module m_FlowCase
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_MeshAndMetries
  use decomp_2d
  use m_Variables,only:mb1,hi1
#ifdef CFDLPT_TwoWay
  use m_Variables,only:FpForce_x,FpForce_y,FpForce_z
#endif
  use m_Tools,only:CalcUxAver
  use iso_c_binding
  implicit none
  private
  include "fftw3.f03"

  ! statistics variabls
  integer,save:: nfstime
  real(RK),save:: PrGradsum
  real(RK),save,allocatable,dimension(:,:):: SumStat 

  ! SpectraOptions
  logical,save:: clc_Spectra 
  integer,save:: ivSpec,jForLCS

  ! Spectra variables
  type(C_PTR),save::fft_plan_x,fft_plan_z
  integer,save:: nxh,nxhp,nzh,nzhp,nSpectime  
  real(RK),save,allocatable,dimension(:,:,:)::EnergySpecX(:,:,:),EnergySpecZ(:,:,:) 
  
  public:: InitVelocity, Update_uy_ym, InitStatVar, clcStat
contains

#if defined(CFDLPT_TwoWay)
#define NCHASTAT 42
#define NEnergySpec 8

#elif defined(ScalarFlow)
#define NCHASTAT 40
#define NEnergySpec 8

#else 
#define NCHASTAT 35
#define NEnergySpec 8
#endif
  !******************************************************************
  ! InitVelocity
  !******************************************************************
  subroutine InitVelocity(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)),intent(inout)::Deviation
  
    ! locals
    integer :: ii,code,i,j,k,m1,m2
    real(RK):: retau_guass,utau_guass,height,rem,wx,wz,xlxPlus,zlzPlus
    real(RK):: xplus,yplus,zplus,yct,ybar,xp,zp,ratiot,uzmean(nyc),uzmeanR(nyc)

    height=yly
    if(FlowType==FT_CH)height=half*yly

    ux=zero; uy=zero; uz=zero
    rem = ubulk * height / xnu
    retau_guass = 0.1538_RK*rem**0.887741_RK
    utau_guass   = retau_guass*xnu/height
    if(nrank==0)print*,'************** retau_gauss=',retau_guass
    if(nrank==0)print*,'************** utau_gauss= ',utau_guass

    call system_clock(count=code); !code=0
    call random_seed(size = ii)
    call random_seed(put = code+63946*(/ (i - 1, i = 1, ii) /))
    call random_number(Deviation)
    Deviation= 0.2_RK* Deviation + 0.9_RK ! [0.8, 1.2]

    !modulation of the random noise + initial velocity profile
    uzmean=zero
    wx=twopi/500.0_RK; wz=twopi/200.0_RK
    xlxPlus=xlx*utau_guass/xnu;   zlzPlus=zlz*utau_guass/xnu;
    m1=floor(xlxPlus*wx/twopi)+1; wx=real(m1,RK)*twopi/xlxPlus
    m2=floor(zlzPlus*wz/twopi)+1; wz=real(m2,RK)*twopi/zlzPlus
    do j=ystart(2),yend(2)
      yct = height-abs(height-yc(j)) 
      ybar= yct/height; yplus=utau_guass*yct/xnu
      do k=ystart(3),yend(3)
        zp   =real(k-1,kind=RK)*dz+dz*half
        zplus=utau_guass*zp/xnu
        do i=ystart(1),yend(1)
          xp   =real(i-1,kind=RK)*dx+dx*half
          xplus=utau_guass*xp/xnu
          !ux(i,j,k) = 0.0052_RK*ubulk*yplus*exp(-yplus*yplus/1800.0_RK)*cos(wz*zplus)*Deviation(i,j,k) ! original expression
          !uz(i,j,k) = 0.0050_RK*ubulk*yplus*exp(-yplus*yplus/1800.0_RK)*sin(wx*xplus)*Deviation(i,j,k) ! original expression
          ux(i,j,k) = ubulk*ybar*exp(-4.5_RK*ybar*ybar)*cos(wz*zplus)*Deviation(i,j,k)
          uz(i,j,k) = ubulk*ybar*exp(-4.5_RK*ybar*ybar)*sin(wx*xplus)*Deviation(i,j,k)
          ux(i,j,k) = ux(i,j,k)+ three*ubulk*(ybar-half*ybar*ybar)
          uzmean(j) = uzmean(j)+ uz(i,j,k)
        enddo
      enddo
      uzmean(j)=uzmean(j)/real(nxc*nzc,RK)
    enddo
    ratiot=ubulk/CalcUxAver(ux)
    call MPI_Bcast(ratiot,1,real_type,0,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uzmean,uzmeanR,nyc,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    ux= ux*ratiot -uCRF
    do j=ystart(2),yend(2)
      uz(:,j,:)=uz(:,j,:)-uzmeanR(j)
    enddo   
  end subroutine InitVelocity

  !******************************************************************
  ! Update_uy_ym
  !******************************************************************   
  subroutine Update_uy_ym(uy_ym, duy_ym)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(out):: uy_ym
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(inout):: duy_ym
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = zero
    duy_ym = uy_ym - duy_ym
  end subroutine Update_uy_ym

  !******************************************************************
  ! InitStatVar
  !******************************************************************
  subroutine InitStatVar(ChannelPrm)
    implicit none 
    character(*),intent(in)::ChannelPrm

    ! locals
    integer::ierror,nUnit
    character(len=50)::filename
    real(RK),dimension(:,:,:),allocatable::Arr1,Arr2
    type(fftw_iodim),dimension(1)::iodim,iodim_howmany
    NAMELIST/SpectraOptions/clc_Spectra,ivSpec,jForLCS

    nUnit = GetFileUnit() 
    open(unit=nUnit, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar", "Cannot open file: "//trim(ChannelPrm))
    read(nUnit, nml=SpectraOptions)
    close(nUnit,IOSTAT=ierror)
    
    if(nrank==0) then
      write(MainLog%nUnit, nml=SpectraOptions)
      if(mod(saveStat,ivstats)/=0 )    call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","ivstats wrong !!!")
      if(clc_Spectra .and. (mod(saveStat,ivSpec)/=0 .or. mod(ivSpec,ivstats)/=0 )) then
        call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","ivSpec wrong !!!")
      endif
      if(IsUxConst)then
        nUnit=GetFileUnit()
        write(filename,'(A,I10.10)')trim(Res_dir)//"PrGrad",ilast
        open(nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
        if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","Cannot open file: "//trim(filename))
        close(nUnit,IOSTAT=ierror)
      endif
    endif
    allocate(SumStat(NCHASTAT,nyp),Stat=ierror)
    if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH: ","Allocation failed")  
    nfstime=0;  nSpectime=0; SumStat=zero; PrGradsum=zero
    if(.not.clc_Spectra) return

    nxh=nxc/2; nxhp=nxh+1
    nzh=nzc/2; nzhp=nzh+1
    allocate(EnergySpecX(nxhp,xsize(2),NEnergySpec));EnergySpecX=zero
    allocate(Arr1(xsize(1),xsize(2),xsize(3)),Arr2(xsize(1),xsize(2),xsize(3)))
    iodim(1)%n  = xsize(1)
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = xsize(2)*xsize(3)
    iodim_howmany(1)%is = xsize(1)
    iodim_howmany(1)%os = xsize(1)
    fft_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,[FFTW_R2HC],FFTW_ESTIMATE)
    deallocate(Arr1,Arr2)
    
    allocate(EnergySpecZ(nzhp,zsize(2),NEnergySpec));EnergySpecZ=zero
    allocate(Arr1(zsize(1),zsize(2),zsize(3)),Arr2(zsize(1),zsize(2),zsize(3)))
    iodim(1)%n  = zsize(3)
    iodim(1)%is = zsize(1)*zsize(2)
    iodim(1)%os = zsize(1)*zsize(2)
    iodim_howmany(1)%n  = zsize(1)*zsize(2)
    iodim_howmany(1)%is = 1
    iodim_howmany(1)%os = 1
    fft_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,[FFTW_R2HC],FFTW_ESTIMATE)
    deallocate(Arr1,Arr2)
  end subroutine InitStatVar

  !******************************************************************
  ! clcStat
  !******************************************************************
#ifdef ScalarFlow
  subroutine clcStat(ux,uy,uz,pressure,scalar,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,pressure,scalar
#else
  subroutine clcStat(ux,uy,uz,pressure,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,pressure
#endif
    real(RK),dimension(ysize(1),ysize(2),ysize(3)),intent(inout)::ArrTemp1,ArrTemp2
   
    ! locals
    character(len=50)::filename
    integer(kind=MPI_OFFSET_KIND)::disp,disp_inc
    real(RK),dimension(:,:),allocatable::arrYplane
    real(RK),dimension(:,:,:),allocatable::arrx1,arrx2,arrz1,arrz2
    real(RK),allocatable,dimension(:,:,:)::EnergySpecXR,EnergySpecZR
    real(RK)::dudxx,dudyy,dudzz,dvdxx,dvdyy,dvdzz,dwdxx,dwdyy,dwdzz,SumStatR(NCHASTAT,nyp),SumVec(NCHASTAT),rdxt
    real(RK)::uxloc,uyloc,uzloc,prloc,uxCell,uyCell,uzCell,prloc2,inxz,infstime,cac,caj,cacU,dudyU,dvdyM
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vor_x,vor_y,vor_z,InterpY1,InterpY2,dudzC,dvdzC,dudxC,dvdxU
    integer::ic,jc,kc,im,jm,km,ip,jp,kp,myistat,ierror,is,ks,iu,ku,coord1,coord2,color,key,SpecWORLDX,SpecWORLDZ,nrankX,nrankZ,nUnit
#ifdef CFDLPT_TwoWay
    real(RK)::uyCells,uyCellm,uyCellp
#endif
#ifdef ScalarFlow
    real(RK)::scloc
#endif
    call myupdate_halo(ux, mb1, hi1)
    call myupdate_halo(uy, mb1, hi1)
    call myupdate_halo(uz, mb1, hi1)
    rdxt=rdx/12.0_RK
    inxz=one/(real(nxc,RK)*real(nzc,RK))
    DO jc=1,nyc
      jp=jc+1; jm=jc-1;
      InterpY1= half*YinterpCoe(jc); InterpY2=half-InterpY1
      cac=rdyc(jc);cacU=rdyc(jp); caj=rdyp(jc); SumVec=zero
      do kc=ystart(3),yend(3)
        ks=kc-2;km=kc-1;kp=kc+1;ku=kc+2
        do ic=ystart(1),yend(1)
          is=ic-2;im=ic-1;ip=ic+1;iu=ic+2

          uxloc= ux(ic,jc,kc)+uCRF
          uyloc= uy(ic,jc,kc)
          uzloc= uz(ic,jc,kc)
          prloc= pressure(ic,jc,kc)
          uxCell= InterpCoe1*ux(im,jc,kc) +InterpCoe2*ux(ic,jc,kc) +InterpCoe3*ux(ip,jc,kc) +InterpCoe4*ux(iu,jc,kc)+uCRF
          uyCell= (uyloc+ uy(ic,jp,kc))*half
          uzCell= InterpCoe1*uz(ic,jc,km) +InterpCoe2*uz(ic,jc,kc) +InterpCoe3*uz(ic,jc,kp) +InterpCoe4*uz(ic,jc,ku) 
          prloc2= InterpY1*(pressure(im,jm,kc)+pressure(ic,jm,kc))+ InterpY2*(pressure(im,jc,kc)+prloc)
 
          dudx= dxCoe1*ux(im,jc,kc) +dxCoe2*ux(ic,jc,kc) +dxCoe3*ux(ip,jc,kc) +dxCoe4*ux(iu,jc,kc)
          dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
          dudz= dzCoe1*ux(ic,jc,ks) +dzCoe2*ux(ic,jc,km) +dzCoe3*ux(ic,jc,kc) +dzCoe4*ux(ic,jc,kp)
          dvdx= dxCoe1*uy(is,jc,kc) +dxCoe2*uy(im,jc,kc) +dxCoe3*uy(ic,jc,kc) +dxCoe4*uy(ip,jc,kc)
          dvdy= (uy(ic,jp,kc)-uyloc)*caj
          dvdz= dzCoe1*uy(ic,jc,ks) +dzCoe2*uy(ic,jc,km) +dzCoe3*uy(ic,jc,kc) +dzCoe4*uy(ic,jc,kp)
          dwdx= dxCoe1*uz(is,jc,kc) +dxCoe2*uz(im,jc,kc) +dxCoe3*uz(ic,jc,kc) +dxCoe4*uz(ip,jc,kc)
          dwdy= (uzloc-uz(ic,jm,kc))*cac
          dwdz= dzCoe1*uz(ic,jc,km) +dzCoe2*uz(ic,jc,kc) +dzCoe3*uz(ic,jc,kp) +dzCoe4*uz(ic,jc,ku)

          dudxx= dxxCoe1*ux(is,jc,kc) +dxxCoe2*ux(im,jc,kc) +dxxCoe3*ux(ic,jc,kc) +dxxCoe4*ux(ip,jc,kc) +dxxCoe5*ux(iu,jc,kc)
          dudyy= ap2c(jc)*ux(ic,jp,kc)+ac2c(jc)*ux(ic,jc,kc)+am2c(jc)*ux(ic,jm,kc)
          dudzz= dzzCoe1*ux(ic,jc,ks) +dzzCoe2*ux(ic,jc,km) +dzzCoe3*ux(ic,jc,kc) +dzzCoe4*ux(ic,jc,kp) +dzzCoe5*ux(ic,jc,ku)
          dvdxx= dxxCoe1*uy(is,jc,kc) +dxxCoe2*uy(im,jc,kc) +dxxCoe3*uy(ic,jc,kc) +dxxCoe4*uy(ip,jc,kc) +dxxCoe5*uy(iu,jc,kc)
          dvdyy= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          dvdzz= dzzCoe1*uy(ic,jc,ks) +dzzCoe2*uy(ic,jc,km) +dzzCoe3*uy(ic,jc,kc) +dzzCoe4*uy(ic,jc,kp) +dzzCoe5*uy(ic,jc,ku)
          dwdxx= dxxCoe1*uz(is,jc,kc) +dxxCoe2*uz(im,jc,kc) +dxxCoe3*uz(ic,jc,kc) +dxxCoe4*uz(ip,jc,kc) +dxxCoe5*uz(iu,jc,kc)
          dwdyy= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)
          dwdzz= dzzCoe1*uz(ic,jc,ks) +dzzCoe2*uz(ic,jc,km) +dzzCoe3*uz(ic,jc,kc) +dzzCoe4*uz(ic,jc,kp) +dzzCoe5*uz(ic,jc,ku)
          vor_x= dwdy-dvdz
          vor_y= dudz-dwdx
          vor_z= dvdx-dudy

          dudxC= (-ux(iu,jc,kc)+8.0_RK*ux(ip,jc,kc)-8.0_RK*ux(im,jc,kc)+ux(is,jc,kc))*rdxt
          dvdxU= dxCoe1*uy(is,jp,kc) +dxCoe2*uy(im,jp,kc) +dxCoe3*uy(ic,jp,kc) +dxCoe4*uy(ip,jp,kc)
          dudyU= (ux(ic,jp,kc)-ux(ic,jc,kc))*cacU
          dvdyM= (uy(im,jp,kc)-uy(im,jc,kc))*caj
          dudzC= (InterpY1*(ux(ic,jm,kp)-ux(ic,jm,km)) +InterpY2*(ux(ic,jc,kp)-ux(ic,jc,km)))*rdz
          dvdzC= (uy(im,jc,kp)-uy(im,jc,km)+ uy(ic,jc,kp)-uy(ic,jc,km))*rdz*quarter

          SumVec( 1)=SumVec( 1)+ uxloc
          SumVec( 2)=SumVec( 2)+ uyloc                     ! yp
          SumVec( 3)=SumVec( 3)+ uzloc
          SumVec( 4)=SumVec( 4)+ prloc
          SumVec( 5)=SumVec( 5)+ uxloc*uxloc
          SumVec( 6)=SumVec( 6)+ uyloc*uyloc               ! yp  
          SumVec( 7)=SumVec( 7)+ uzloc*uzloc
          SumVec( 8)=SumVec( 8)+ prloc*prloc
          SumVec( 9)=SumVec( 9)+ uxCell*uyCell
          SumVec(10)=SumVec(10)+ uyCell*uzCell
          SumVec(11)=SumVec(11)+ uxCell*uzCell
          SumVec(12)=SumVec(12)+ uxCell*prloc
          SumVec(13)=SumVec(13)+ uyCell*prloc
          SumVec(14)=SumVec(14)+ uzCell*prloc
          SumVec(15)=SumVec(15)+ uxCell*uxCell*uyCell 
          SumVec(16)=SumVec(16)+ uyloc *uyloc *uyloc       ! yp
          SumVec(17)=SumVec(17)+ uzCell*uzCell*uyCell
          SumVec(18)=SumVec(18)+ uxCell*uyCell*uyCell
          SumVec(19)=SumVec(19)+ uxloc*uxloc*uxloc
          SumVec(20)=SumVec(20)+ uzloc*uzloc*uzloc
          SumVec(21)=SumVec(21)+ uxloc*uxloc*uxloc*uxloc
          SumVec(22)=SumVec(22)+ uyloc*uyloc*uyloc*uyloc   ! yp
          SumVec(23)=SumVec(23)+ uzloc*uzloc*uzloc*uzloc
          SumVec(24)=SumVec(24)+ prloc*dudx
          SumVec(25)=SumVec(25)+ prloc*dvdy
          SumVec(26)=SumVec(26)+ prloc*dwdz
          SumVec(27)=SumVec(27)+ prloc2*(dudy+dvdx)        ! yp !
          SumVec(28)=SumVec(28)+ uxloc*(dudxx+dudyy+dudzz)
          SumVec(29)=SumVec(29)+ uyloc*(dvdxx+dvdyy+dvdzz) ! yp !
          SumVec(30)=SumVec(30)+ uzloc*(dwdxx+dwdyy+dwdzz)
          SumVec(31)=SumVec(31)+ dudxC*(dvdxU+dvdx)*half+ (dudyU+dudy)*(dvdyM+dvdy)*quarter
          SumVec(32)=SumVec(32)+ dudzC*dvdzC          ! yp !
          SumVec(33)=SumVec(33)+ vor_x*vor_x          ! yp !
          SumVec(34)=SumVec(34)+ vor_y*vor_y
          SumVec(35)=SumVec(35)+ vor_z*vor_z          ! yp !
#ifdef CFDLPT_TwoWay
#ifdef ForceInCell
          SumVec(36)=SumVec(36)+ FpForce_x(ic,jc,kc)
          SumVec(37)=SumVec(37)+ FpForce_y(ic,jc,kc)
          SumVec(38)=SumVec(38)+ FpForce_z(ic,jc,kc)
          SumVec(39)=SumVec(39)+ FpForce_x(ic,jc,kc)*uxCell
          SumVec(40)=SumVec(40)+ FpForce_y(ic,jc,kc)*uyCell
          SumVec(41)=SumVec(41)+ FpForce_z(ic,jc,kc)*uzCell
          SumVec(42)=SumVec(42)+ FpForce_y(ic,jc,kc)*uxCell +FpForce_x(ic,jc,kc)*uyCell
#else
          SumVec(36)=SumVec(36)+ FpForce_x(ic,jc,kc)
          SumVec(37)=SumVec(37)+ FpForce_y(ic,jc,kc)  ! yp !
          SumVec(38)=SumVec(38)+ FpForce_z(ic,jc,kc)
          SumVec(39)=SumVec(39)+ FpForce_x(ic,jc,kc)*uxloc
          SumVec(40)=SumVec(40)+ FpForce_y(ic,jc,kc)*uyloc ! yp !
          SumVec(41)=SumVec(41)+ FpForce_z(ic,jc,kc)*uzloc
          SumVec(42)=SumVec(42)+ half*(FpForce_y(ic,jc,kc)+FpForce_y(ic,jp,kc))*uxCell
          
          uyCells=(uy(is,jc,kc)+ uy(is,jp,kc))*half
          uyCellm=(uy(im,jc,kc)+ uy(im,jp,kc))*half
          uyCellp=(uy(ip,jc,kc)+ uy(ip,jp,kc))*half
          SumVec(42)=SumVec(42)+ FpForce_x(ic,jc,kc)*(InterpCoe1*uyCells +InterpCoe2*uyCellm +InterpCoe3*uyCell +InterpCoe4*uyCellp) 
#endif
#endif
#ifdef ScalarFlow
          scloc=scalar(ic,jc,kc)
          SumVec(36)=SumVec(36)+ scloc
          SumVec(37)=SumVec(37)+ scloc*scloc
          SumVec(38)=SumVec(38)+ scloc*uxCell
          SumVec(39)=SumVec(39)+ scloc*uyCell
          SumVec(40)=SumVec(40)+ scloc*uzCell
#endif
        enddo
      enddo
      do kc=1,NCHASTAT
        SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
      enddo
    ENDDO

    ! nyp only
    jc=nyp; jm=jc-1; cac=rdyc(jc); SumVec=zero
    InterpY1= half*YinterpCoe(jc); InterpY2=half-InterpY1
    do kc=ystart(3),yend(3)
      km=kc-1;kp=kc+1
      do ic=ystart(1),yend(1)
        im=ic-1;ip=ic+1
        prloc2= InterpY1*(pressure(im,jm,kc)+pressure(ic,jm,kc))+ InterpY2*(pressure(im,jc,kc)+pressure(ic,jc,kc))
        dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
        dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
        vor_x=  dwdy
        vor_z= -dudy
        SumVec(27)=SumVec(27)+ prloc2*(dudy+zero)   ! yp !
        SumVec(33)=SumVec(33)+ vor_x*vor_x          ! yp !
        SumVec(35)=SumVec(35)+ vor_z*vor_z          ! yp !
#if defined(CFDLPT_TwoWay) && !defined(ForceInCell)
        SumVec(37)=SumVec(37)+ FpForce_y(ic,jc,kc)              ! yp !
        SumVec(40)=SumVec(40)+ FpForce_y(ic,jc,kc)*uy(ic,jc,kc) ! yp !
#endif
      enddo
    enddo
    do kc=1,NCHASTAT
      SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
    enddo

    ! shear stress and pressure gradient
    PrGradsum   = PrGradsum+ PrGradAver
    if(nrank==0 .and. IsUxConst) then
      nUnit=GetFileUnit()
      write(filename,'(A,I10.10)')trim(Res_dir)//"PrGrad",ilast
      open(nUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=myistat)
      if(myistat/=0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      else
        write(nUnit,'(I7,2ES24.15)')itime,SimTime,PrGradAver
      endif
      close(nUnit,IOSTAT=myistat)
    endif
    nfstime= nfstime + 1

    ! Calculate Energy Spectra
    if(clc_Spectra .and. mod(itime,ivSpec)==0) then

      ! Determine ux in cell center
      Block 
      integer::it,jt,kt
      DO kc=ystart(3),yend(3)
        kt=kc-ystart(3)+1
        do jc=ystart(2),yend(2)
          jt=jc
          do ic=ystart(1),yend(1)
            it=ic-ystart(1)+1
            im=ic-1;ip=ic+1;iu=ic+2
            ArrTemp1(it,jt,kt)=InterpCoe1*ux(im,jc,kc) +InterpCoe2*ux(ic,jc,kc) +InterpCoe3*ux(ip,jc,kc) +InterpCoe4*ux(iu,jc,kc)+uCRF
          enddo
        enddo
      ENDDO
      End block
      
      !=============== Spectra in x-dir ===============
      allocate(arrYplane(ysize(1),ysize(3)))
      allocate(arrx1(xsize(1),xsize(2),xsize(3)))
      allocate(arrx2(xsize(1),xsize(2),xsize(3)))
      call transpose_y_to_x(ux(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrx1); arrx1=arrx1+uCRF
      call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
      call clcEnergySpectraX(arrx1,1)
      
      ! Calculate Linear coherence spectrum (LCS)
      ! Real part
      call transpose_x_to_y(arrx1,ArrTemp2)
      do kc=1,ysize(3)
        do ic=1,ysize(1)
          arrYplane(ic,kc)=ArrTemp2(ic,jForLCS,kc)
        enddo
      enddo
      do kc=1,ysize(3)
        do jc=1,ysize(2)
          do ic=1,ysize(1)
            ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc)
          enddo
        enddo
      enddo
      call transpose_y_to_x(ArrTemp2,arrx2)
      call clcLTS_x_real(arrx2,7)
      
      ! Imaginary part
      arrx2(1,:,:)=arrx1(1,:,:)
      do kc=1,xsize(3)
        do jc=1,xsize(2)
          do ic=2,xsize(1)
            arrx2(ic,jc,kc)=arrx1(nxc+2-ic,jc,kc)
          enddo
        enddo
      enddo
      call transpose_x_to_y(arrx2,ArrTemp2)
      do kc=1,ysize(3)
        do jc=1,ysize(2)
          do ic=1,ysize(1)
            ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc)
          enddo
        enddo
      enddo
      call transpose_y_to_x(ArrTemp2,arrx2)
      call clcLTS_x_imag(arrx2,8)
      
      call transpose_y_to_x(uy(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrx2)
      call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
      call clcEnergySpectraX(arrx2,2)

      call transpose_y_to_x(uz(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrx2)
      call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
      call clcEnergySpectraX(arrx2,3)

      call transpose_y_to_x(pressure(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrx2)
      call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
      call clcEnergySpectraX(arrx2,4)

      ! Determine ux and uy in cell center
      Block 
      integer::it,jt,kt
      DO kc=ystart(3),yend(3)
        kt=kc-ystart(3)+1
        do jc=ystart(2),yend(2)
          jt=jc
          jp=jc+1
          do ic=ystart(1),yend(1)
            it=ic-ystart(1)+1
            ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
          enddo
        enddo
      ENDDO
      End block
      call transpose_y_to_x(ArrTemp2,arrx2)
      call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
      call clcCospectraX(arrx1,arrx2,5)
      call transpose_y_to_x(ArrTemp1,arrx1)
      call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
      call clcCospectraX(arrx1,arrx2,6)
      deallocate(arrx1,arrx2)
       
      !=============== Spectra in z-dir ===============
      allocate(arrz1(zsize(1),zsize(2),zsize(3)))
      allocate(arrz2(zsize(1),zsize(2),zsize(3)))
      call transpose_y_to_z(ux(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrz1); arrz1=arrz1+uCRF
      call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
      call clcEnergySpectraZ(arrz1,1)

      ! Calculate Linear coherence spectrum (LCS)
      ! Real part
      call transpose_z_to_y(arrz1,ArrTemp2)
      do kc=1,ysize(3)
        do ic=1,ysize(1)
          arrYplane(ic,kc)=ArrTemp2(ic,jForLCS,kc)
        enddo
      enddo
      do kc=1,ysize(3)
        do jc=1,ysize(2)
          do ic=1,ysize(1)
            ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc)
          enddo
        enddo
      enddo
      call transpose_y_to_z(ArrTemp2,arrz2)
      call clcLTS_z_real(arrz2,7)
      
      ! Imaginary part
      arrz2(:,:,1)=arrz1(:,:,1)
      do kc=2,zsize(3)
        do jc=1,zsize(2)
          do ic=1,zsize(1)
            arrz2(ic,jc,kc)=arrz1(ic,jc,nzc+2-kc)
          enddo
        enddo
      enddo
      call transpose_z_to_y(arrz2,ArrTemp2)
      do kc=1,ysize(3)
        do jc=1,ysize(2)
          do ic=1,ysize(1)
            ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc)
          enddo
        enddo
      enddo
      call transpose_y_to_z(ArrTemp2,arrz2)
      call clcLTS_z_imag(arrz2,8)
      
      call transpose_y_to_z(uy(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrz2)
      call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
      call clcEnergySpectraZ(arrz2,2)

      call transpose_y_to_z(uz(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrz2)
      call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
      call clcEnergySpectraZ(arrz2,3)

      call transpose_y_to_z(pressure(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrz2)
      call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
      call clcEnergySpectraZ(arrz2,4)

      ! Determine ux and uy in cell center
      Block 
      integer::it,jt,kt
      DO kc=ystart(3),yend(3)
        kt=kc-ystart(3)+1
        do jc=ystart(2),yend(2)
          jt=jc
          jp=jc+1
          do ic=ystart(1),yend(1)
            it=ic-ystart(1)+1
            ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
          enddo
        enddo
      ENDDO
      End block       
      call transpose_y_to_z(ArrTemp2,arrz2)
      call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
      call clcCospectraZ(arrz1,arrz2,5)
      call transpose_y_to_z(ArrTemp1,arrz1)
      call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
      call clcCospectraZ(arrz1,arrz2,6)     
      deallocate(arrz1,arrz2)
      deallocate(arrYplane)

      nSpectime= nSpectime+1
    endif
    if(mod(itime,SaveStat)/=0) return

    ! Write Energy Spectra
    IF(clc_Spectra) THEN
      infstime= one/real(nSpectime,RK)
      coord1=int(nrank/p_col); coord2=mod(nrank,p_col)

      ! Spectra in x-dir
      write(filename,"(A,A,I10.10)") trim(Res_Dir),'SpecX',itime
      call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,nUnit,ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      disp_inc=int(mytype_bytes,8)*int(nyc,8)*int(nxhp,8)
       
      allocate(EnergySpecXR(nxhp,xsize(2),NEnergySpec));EnergySpecXR=zero
      color=coord1; key=nrank
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,SpecWORLDX,ierror)
      call MPI_REDUCE(EnergySpecX,EnergySpecXR,NEnergySpec*xsize(2)*nxhp,real_type,MPI_SUM,0,SpecWORLDX,ierror)
      call MPI_COMM_RANK(SpecWORLDX,nrankX,ierror)
      if(nrankX==0) then
        EnergySpecXR=EnergySpecXR*infstime
        disp=int(mytype_bytes,8)*int(xstart(2)-1,8)*int(nxhp,8)
        do kc=1,NEnergySpec
          call MPI_FILE_WRITE_AT(nUnit,disp,EnergySpecXR(:,:,kc),xsize(2)*nxhp,real_type,MPI_STATUS_IGNORE,ierror);disp=disp+disp_inc
        enddo
      endif
      deallocate(EnergySpecXR)
      call MPI_COMM_FREE(SpecWORLDX,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)

      ! Spectra in z-dir
      write(filename,"(A,A,I10.10)") trim(Res_Dir),'SpecZ',itime
      call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,nUnit,ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      disp_inc=int(mytype_bytes,8)*int(nyc,8)*int(nzhp,8)
            
      allocate(EnergySpecZR(nzhp,zsize(2),NEnergySpec));EnergySpecZR=zero
      color=coord2; key=nrank
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,SpecWORLDZ,ierror)
      call MPI_REDUCE(EnergySpecZ,EnergySpecZR,NEnergySpec*zsize(2)*nzhp,real_type,MPI_SUM,0,SpecWORLDZ,ierror)
      call MPI_COMM_RANK(SpecWORLDZ,nrankZ,ierror)
      if(nrankZ==0) then
        EnergySpecZR=EnergySpecZR*infstime
        disp=int(mytype_bytes,8)*int(zstart(2)-1,8)*int(nzhp,8)
        do kc=1,NEnergySpec
          call MPI_FILE_WRITE_AT(nUnit,disp,EnergySpecZR(:,:,kc),zsize(2)*nzhp,real_type,MPI_STATUS_IGNORE,ierror);disp=disp+disp_inc
        enddo
      endif
      deallocate(EnergySpecZR)
      call MPI_COMM_FREE(SpecWORLDZ,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    ENDIF

    ! Write statistics
    call MPI_REDUCE(SumStat,SumStatR,NCHASTAT*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = one/real(nfstime,RK)
      nUnit=GetFileUnit()
      write(filename,"(A,I10.10)") trim(Res_Dir)//'stats',itime
      open(nUnit,file=filename,status='replace',form='formatted',IOSTAT=myistat)
      IF(myistat/=0) THEN
        call MainLog%CheckForError(ErrT_Abort,"clcStat_CH","Cannot open file: "//trim(filename))
      ELSE
        write(nUnit,'(a,I7,a,I7,a,I7)')'  The time step range for this fluid statistics is ', &
                                    itime-(nfstime-1)*ivstats, ':', ivstats, ':', itime
        write(nUnit,'(A)')'  '
        if(IsUxConst) then
          write(nUnit,'(A)')'  Constant velocity in x-dir by adding a pressure gradient.'
          write(nUnit,'(A, ES24.15)')'    time averaged pressure gradient is: ',abs(PrGradsum)*infstime
        else
          write(nUnit,'(A)')'  Variable velocity in x-dir while adding a constant body force.'
        endif
        write(nUnit,'(A)')'  '
        
        Block 
        character(len=50)::FormatStr
        write(FormatStr,'(A,I3,A)')'(',NCHASTAT,'ES24.15)'
        do jc=1,nyp
          write(nUnit,FormatStr)SumStatR(1:NCHASTAT,jc)*infstime
        enddo
        End block
      ENDIF
      close(nUnit,IOSTAT=myistat)
    endif

    nfstime=0; SumStat=zero; PrGradsum=zero; nSpectime=0; 
    if(clc_Spectra) then
      EnergySpecX=zero; EnergySpecZ=zero
    endif
  end subroutine clcStat

  !******************************************************************
  ! clcEnergySpectraX
  !******************************************************************
  subroutine clcEnergySpectraX(arrx,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(xsize(1),xsize(2),xsize(3)),intent(in)::arrx

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,xsize(2)
      EgyX=zero
      do kc=1,xsize(3)
        EgyX(1)   = EgyX(1)   + arrx(1,jc,kc)   *arrx(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx(nxhp,jc,kc)*arrx(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx(ic,jc,kc)*arrx(ic,jc,kc)+arrx(icCnter,jc,kc)*arrx(icCnter,jc,kc))*two
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraX

  !******************************************************************
  ! clcLTS_x_real
  !******************************************************************
  subroutine clcLTS_x_real(arrx,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(xsize(1),xsize(2),xsize(3)),intent(in)::arrx

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,xsize(2)
      EgyX=zero
      do kc=1,xsize(3)
        EgyX(1)   = EgyX(1)   + arrx(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx(ic,jc,kc)+arrx(icCnter,jc,kc)) ! Note, there is NO "*two" here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcLTS_x_real

  !******************************************************************
  ! clcLTS_x_imag
  !******************************************************************
  subroutine clcLTS_x_imag(arrx,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(xsize(1),xsize(2),xsize(3)),intent(in)::arrx

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,xsize(2)
      EgyX=zero
      do kc=1,xsize(3)
        EgyX(1)   = zero
        EgyX(nxhp)= zero
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx(ic,jc,kc)-arrx(icCnter,jc,kc)) ! Note here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcLTS_x_imag
   
  !******************************************************************
  ! clcCospectraX
  !******************************************************************
  subroutine clcCospectraX(arrx1,arrx2,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(xsize(1),xsize(2),xsize(3)),intent(in)::arrx1,arrx2

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,xsize(2)
      EgyX=zero
      do kc=1,xsize(3)
        EgyX(1)   = EgyX(1)   + arrx1(1,jc,kc)   *arrx2(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx1(nxhp,jc,kc)*arrx2(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx1(ic,jc,kc)*arrx2(ic,jc,kc)+arrx1(icCnter,jc,kc)*arrx2(icCnter,jc,kc))*two
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcCospectraX

  !******************************************************************
  ! clcEnergySpectraZ
  !******************************************************************
  subroutine clcEnergySpectraZ(arrz,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(zsize(1),zsize(2),zsize(3)),intent(in)::arrz

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,zsize(2)
      EgyZ=zero
      do ic=1,zsize(1)
        EgyZ(1)   = EgyZ(1)   + arrz(ic,jc,1)   *arrz(ic,jc,1)
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz(ic,jc,nzhp)*arrz(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz(ic,jc,kc)*arrz(ic,jc,kc)+arrz(ic,jc,kcCnter)*arrz(ic,jc,kcCnter))*two
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraZ

  !******************************************************************
  ! clcLTS_z_real
  !******************************************************************
  subroutine clcLTS_z_real(arrz,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(zsize(1),zsize(2),zsize(3)),intent(in)::arrz

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,zsize(2)
      EgyZ=zero
      do ic=1,zsize(1)
        EgyZ(1)   = EgyZ(1)   + arrz(ic,jc,1)   
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz(ic,jc,kc)+arrz(ic,jc,kcCnter))
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcLTS_z_real

  !******************************************************************
  ! clcLTS_z_imag
  !******************************************************************
  subroutine clcLTS_z_imag(arrz,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(zsize(1),zsize(2),zsize(3)),intent(in)::arrz

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,zsize(2)
      EgyZ=zero
      do ic=1,zsize(1)
        EgyZ(1)   = EgyZ(1)   + zero   
        EgyZ(nzhp)= EgyZ(nzhp)+ zero
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz(ic,jc,kc)-arrz(ic,jc,kcCnter))
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcLTS_z_imag
    
  !******************************************************************
  ! clcCospectraZ
  !******************************************************************
  subroutine clcCospectraZ(arrz1,arrz2,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(zsize(1),zsize(2),zsize(3)),intent(in)::arrz1,arrz2

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,zsize(2)
      EgyZ=zero
      do ic=1,zsize(1)
        EgyZ(1)   = EgyZ(1)   + arrz1(ic,jc,1)   *arrz2(ic,jc,1)
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz1(ic,jc,nzhp)*arrz2(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz1(ic,jc,kc)*arrz2(ic,jc,kc)+arrz1(ic,jc,kcCnter)*arrz2(ic,jc,kcCnter))*two
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcCospectraZ
#undef NCHASTAT    
#undef NEnergySpec
end module m_FlowCase
