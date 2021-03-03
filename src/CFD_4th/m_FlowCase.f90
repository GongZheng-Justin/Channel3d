module m_FlowCase
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_MeshAndMetries
  use decomp_2d
  use m_Variables,only: mb1,hi1
  use m_Tools,only:CalcUxAver
  use mydecomp_2d_extra
  use iso_c_binding
  implicit none
  private
  include "fftw3.f"

  ! statistics variabls
  integer,save:: nfstime
  real(RK),save:: PrGradsum
  real(RK),save,allocatable,dimension(:,:):: SumStat 

  ! Spectra variables
  integer,save:: nSpectime
  integer,save:: njSpecX,njSpecZ
  integer,save:: nxh,nxhp,nzh,nzhp
  integer,save,allocatable,dimension(:)::jSpecX,jSpecZ
  real(RK),save,allocatable,dimension(:,:,:)::EnergySpecX(:,:,:),EnergySpecZ(:,:,:) 

  type,bind(C)::fftw_iodim
    integer(C_INT) n,is,os
  end type fftw_iodim
  interface
    type(C_PTR) function fftw_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,matin,matout,kindt,flags)  &
      bind(C, name='fftw_plan_guru_r2r')
      import
      integer(C_INT), value :: rank
      type(fftw_iodim), dimension(*), intent(in) :: dims
      integer(C_INT), value :: howmany_rank
      type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
      real(C_DOUBLE), dimension(*), intent(out) :: matin,matout
      integer(C_INT) :: kindt
      integer(C_INT), value :: flags
    end function fftw_plan_guru_r2r
  end interface
  type(C_PTR),save:: fft_plan_x
  type(C_PTR),save:: fft_plan_z 

  public:: InitVelocity, Update_uy_ym,InitStatVar,  clcStat
contains

  !******************************************************************
  ! InitVelocity
  !******************************************************************
  subroutine InitVelocity(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)),intent(inout)::Deviation
  
    ! locals
    integer :: ii,code,i,j,k,m,m1,m2
	  real(RK):: retau_guass,utau_guass,height,rem,wx,wz,xlxPlus,zlzPlus
	  real(RK):: xplus,yplus,zplus,yct,ybar,xp,zp,ratiot,uzmean(nyc),uzmeanR(nyc)
      
    select case(FlowType)
    case(FT_CH)
      height=half*yly
    case(FT_HC)
      height=yly
    endselect
    ux=zero; uy=zero; uz=zero
	  rem = ubulk * height / xnu
    retau_guass = 0.1538_RK*rem**0.887741_RK
	  utau_guass 	= retau_guass*xnu/height
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
          ux(i,j,k) = 0.0052_RK*ubulk*yplus*exp(-yplus*yplus/1800.0_RK)*cos(wz*zplus)*Deviation(i,j,k) ! original expression
          uz(i,j,k) = 0.0050_RK*ubulk*yplus*exp(-yplus*yplus/1800.0_RK)*sin(wx*xplus)*Deviation(i,j,k) ! original expression
          !ux(i,j,k) = ubulk*ybar*exp(-4.5_RK*ybar*ybar)*cos(wz*zplus)*Deviation(i,j,k)
          !uz(i,j,k) = ubulk*ybar*exp(-4.5_RK*ybar*ybar)*sin(wx*xplus)*Deviation(i,j,k)
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
  subroutine Update_uy_ym(uy_ym, duy_ym, TimeNew)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(out):: uy_ym
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(inout):: duy_ym    
    real(RK),intent(in):: TimeNew
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = zero
    duy_ym = uy_ym - duy_ym
  end subroutine Update_uy_ym

  !******************************************************************
  ! InitStatVar
  !******************************************************************
  subroutine InitStatVar()
    implicit none

    ! locals
    integer::iErr,jc
    real(RK),dimension(:,:,:),allocatable::arrx,arrz
    type(fftw_iodim),dimension(1)::iodim,iodim_howmany

    if(nrank==0) then
      if(jSpecSet<1 .or. jSpecEnd>nyc) call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","jSpecSet or jSpecEnd wrong !!!")
      if(mod(saveStat,ivstats)/=0 )    call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","ivstats wrong !!!")
      if(clc_Spectra .and. (mod(saveStat,ivSpec)/=0 .or. mod(ivSpec,ivstats)/=0 )) then
        call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","ivSpec wrong !!!")
      endif
      if(IsUxConst)then
        open(79,file='./CFD/Results/PrGrad.txt',status='replace',form='formatted',IOSTAT=iErr)
        if(iErr /= 0) call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","Cannot open file")
        close(79,IOSTAT=iErr)
      endif
    endif
    allocate(SumStat(35,nyp),Stat=iErr)
    if(iErr /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH: ","Allocation failed")  
    nfstime=0;  nSpectime=0; SumStat=zero; PrGradsum=zero
    if(.not.clc_Spectra) return

    nxh=nxc/2; nxhp=nxh+1
    nzh=nzc/2; nzhp=nzh+1
    njSpecX=0; njSpecZ=0
    do jc=jSpecSet,jSpecEnd,jSpecInc
      if(jc>=xstart(2) .and. jc<=xend(2)) njSpecX=njSpecX+1
      if(jc>=zstart(2) .and. jc<=zend(2)) njSpecZ=njSpecZ+1
    enddo
    if(njSpecX>0) then
      allocate(jSpecX(njSpecX), EnergySpecX(4,njSpecX,nxhp));EnergySpecX=zero
      allocate(arrx(xsize(1),njSpecX,xsize(3)))
      iodim(1)%n  = xsize(1)
      iodim(1)%is = 1
      iodim(1)%os = 1
      iodim_howmany(1)%n  = njSpecX*xsize(3)
      iodim_howmany(1)%is = xsize(1)
      iodim_howmany(1)%os = xsize(1)
      fft_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrx,arrx,FFTW_R2HC,FFTW_ESTIMATE)
    endif
    if(njSpecZ>0) then
      allocate(jSpecZ(njSpecZ), EnergySpecZ(4,njSpecZ,nzhp));EnergySpecZ=zero
      allocate(arrz(zsize(1),njSpecZ,zsize(3)))
      iodim(1)%n  = zsize(3)
      iodim(1)%is = zsize(1)*njSpecZ
      iodim(1)%os = zsize(1)*njSpecZ
      iodim_howmany(1)%n  = zsize(1)*njSpecZ
      iodim_howmany(1)%is = 1
      iodim_howmany(1)%os = 1
      fft_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrz,arrz,FFTW_R2HC,FFTW_ESTIMATE)
    endif

    njSpecX=0; njSpecZ=0
    do jc=jSpecSet,jSpecEnd,jSpecInc
      if(jc>=xstart(2) .and. jc<=xend(2)) then
        njSpecX=njSpecX+1; jSpecX(njSpecX)=jc
      endif
      if(jc>=zstart(2) .and. jc<=zend(2)) then
        njSpecZ=njSpecZ+1; jSpecZ(njSpecZ)=jc
      endif
    enddo
  end subroutine InitStatVar

  !******************************************************************
  ! clcStat
  !******************************************************************
  subroutine clcStat(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,pressure
   
    ! locals
    character(len=50)::filename,myformat
    real(RK),dimension(xsize(1),xstart(2):xend(2),xsize(3))::arrx
    real(RK),dimension(zsize(1),zstart(2):zend(2),zsize(3))::arrz
    real(RK),allocatable,dimension(:,:,:)::EnergySpecXR,EnergySpecZR
    real(RK)::dudxx,dudyy,dudzz,dvdxx,dvdyy,dvdzz,dwdxx,dwdyy,dwdzz,SumStatR(35,nyp),SumVec(35),rdxt
    real(RK)::uxloc,uyloc,uzloc,prloc,uxCell,uyCell,uzCell,prloc2,inxz,infstime,cac,caj,cacU,dudyU,dvdyM
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vor_x,vor_y,vor_z,InterpY1,InterpY2,dudzC,dvdzC,dudxC,dvdxU
    integer::ic,jc,kc,im,jm,km,ip,jp,kp,myistat,ierror,is,ks,iu,ku,coord1,coord2,color,key,SpecWORLDX,SpecWORLDZ,nrankX,nrankZ

    call myupdate_halo(ux, mb1, hi1)
    call myupdate_halo(uy, mb1, hi1)
    call myupdate_halo(uz, mb1, hi1)
    rdxt=rdx/12.0_RK
    inxz=one/(real(nxc,RK)*real(nzc,RK))
    DO jc=1,nyc
      jp=jc+1; jm=jc-1;
      InterpY1= half*YinterpCoe(jm); InterpY2=half-InterpY1
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
          SumVec(29)=SumVec(29)+ uyloc*(dvdxx+dvdyy+dvdzz) ! yp
          SumVec(30)=SumVec(30)+ uzloc*(dwdxx+dwdyy+dwdzz)
          SumVec(31)=SumVec(31)+ dudxC*(dvdxU+dvdx)*half+ (dudyU+dudy)*(dvdyM+dvdy)*quarter
          SumVec(32)=SumVec(32)+ dudzC*dvdzC          ! yp
          SumVec(33)=SumVec(33)+ vor_x*vor_x          ! yp !
          SumVec(34)=SumVec(34)+ vor_y*vor_y
          SumVec(35)=SumVec(35)+ vor_z*vor_z          ! yp !
        enddo
      enddo
      do kc=1,35
        SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
      enddo
    ENDDO

    ! nyp only
    jc=nyp; jm=jc-1; cac=rdyc(jc); SumVec=zero
    InterpY1= half*YinterpCoe(jm); InterpY2=half-InterpY1
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
      enddo
    enddo
    do kc=1,35
      SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
    enddo

    ! shear stress and pressure gradient
    PrGradsum   = PrGradsum+ PrGradAver
    if(nrank==0 .and. IsUxConst) then
      open(79,file='./CFD/Results/PrGrad.txt',status='old',position='append',form='formatted',IOSTAT=myistat)
      if(myistat /= 0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file 1")
      else
        write(79,'(I7,2ES24.15)')itime,SimTime,PrGradAver
      endif
      close(79,IOSTAT=myistat)
    endif
    nfstime= nfstime + 1

    ! Calculate Energy Spectra
    if(clc_Spectra .and. mod(itime,ivSpec)==0) then
      ! spectra in x-dir
      call transpose_y_to_x(ux(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrx); arrx=arrx+uCRF
      if(njSpecX>0) call clcEnergySpectraX(arrx,1)

      call transpose_y_to_x(uy(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrx)
      if(njSpecX>0) call clcEnergySpectraX(arrx,2)

      call transpose_y_to_x(uz(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrx)
      if(njSpecX>0) call clcEnergySpectraX(arrx,3)

      call transpose_y_to_x(pressure(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrx)
      if(njSpecX>0) call clcEnergySpectraX(arrx,4)

      ! spectra in z-dir
      call transpose_y_to_z(ux(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrz); arrz=arrz+uCRF
      if(njSpecZ>0) call clcEnergySpectraZ(arrz,1)

      call transpose_y_to_z(uy(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrz)
      if(njSpecZ>0) call clcEnergySpectraZ(arrz,2)

      call transpose_y_to_z(uz(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrz)
      if(njSpecZ>0) call clcEnergySpectraZ(arrz,3)

      call transpose_y_to_z(pressure(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),arrz)
      if(njSpecZ>0) call clcEnergySpectraZ(arrz,4)

      nSpectime= nSpectime+1
    endif
    if(mod(itime,SaveStat)/=0) return

    ! Write Energy Spectra
    if(clc_Spectra) then
      infstime= one/real(nSpectime,RK)
      coord1=int(nrank/p_col); coord2=mod(nrank,p_col)

      ! Spectra in x-dir 
      color=coord1; key=nrank
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,SpecWORLDX,ierror) 
      if(njSpecX>0) then
        allocate(EnergySpecXR(4,njSpecX,nxhp));EnergySpecXR=zero
        call MPI_REDUCE(EnergySpecX,EnergySpecXR,4*njSpecX*nxhp,real_type,MPI_SUM,0,SpecWORLDX,ierror)
        call MPI_COMM_RANK(SpecWORLDX,nrankX,ierror)
        if(nrankX==0) then
          EnergySpecXR=EnergySpecXR*infstime
          write(filename,"(A,I10.10,A,I4.4)") './CFD/Results/SpecX',itime,'_',coord1
          open(70,file=trim(filename),status='replace',form='formatted',IOSTAT=myistat)
          IF(myistat/=0) THEN
            call MainLog%CheckForError(ErrT_Abort,"clcStat_CH","open file wrong 1")
          ELSE
            write(myformat,'(A,I5,A)')'(',4*njSpecX,'ES24.15)'
            do ic=1,nxhp
              write(70,myformat)EnergySpecXR(:,:,ic)
            enddo
          ENDIF
          close(70,IOSTAT=myistat)
        endif
        deallocate(EnergySpecXR)
      endif
      call MPI_COMM_FREE(SpecWORLDX,ierror)

      ! Spectra in z-dir
      color=coord2; key=nrank
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,SpecWORLDZ,ierror)
      if(njSpecZ>0) then
        allocate(EnergySpecZR(4,njSpecZ,nzhp));EnergySpecZR=zero
        call MPI_REDUCE(EnergySpecZ,EnergySpecZR,4*njSpecZ*nzhp,real_type,MPI_SUM,0,SpecWORLDZ,ierror)
        call MPI_COMM_RANK(SpecWORLDZ,nrankZ,ierror)
        if(nrankZ==0) then
          EnergySpecZR=EnergySpecZR*infstime
          write(filename,"(A,I10.10,A,I4.4)") './CFD/Results/SpecZ',itime,'_',coord2
          open(71,file=trim(filename),status='replace',form='formatted',IOSTAT=myistat)
          IF(myistat/=0) THEN
            call MainLog%CheckForError(ErrT_Abort,"clcStat_CH","open file wrong 2")
          ELSE
            write(myformat,'(A,I5,A)')'(',4*njSpecZ,'ES24.15)'
            do kc=1,nzhp
              write(71,myformat)EnergySpecZR(:,:,kc)
            enddo
          ENDIF
          close(71,IOSTAT=myistat)
        endif
        deallocate(EnergySpecZR)
      endif
      call MPI_COMM_FREE(SpecWORLDZ,ierror)
    endif

    ! Write statistics
    call MPI_REDUCE(SumStat,SumStatR,35*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = one/real(nfstime,RK)
      write(filename,"('./CFD/Results/stats',I10.10)") itime
      open(69,file=trim(filename),status='replace',form='formatted',IOSTAT=myistat)
      IF(myistat/=0) THEN
        print*,"open file wrong: ",trim(filename)," please check it! "
      ELSE
        write(69,'(a,I7,a,I7,a,I7)')'  The time step range for this fluid statistics is ', &
                                    itime-(nfstime-1)*ivstats, ':', ivstats, ':', itime
        write(69,'(A)')'  '
        if(IsUxConst) then
          write(69,'(A)')'  Constant velocity in x-dir by adding a pressure gradient.'
          write(69,'(A, ES24.15)')'    time averaged pressure gradient is: ',abs(PrGradsum)*infstime
        else
          write(69,'(A)')'  Variable velocity in x-dir while adding a constant body force.'
        endif
        write(69,'(A)')'  '
        do jc=1,nyp
          write(69,'(40ES24.15)')SumStatR(1:35,jc)*infstime
        enddo
      ENDIF
      close(69,IOSTAT=myistat)
    endif

    nfstime=0; SumStat=zero; PrGradsum=zero; nSpectime=0; 
    if(njSpecX>0) EnergySpecX=zero; if(njSpecZ>0) EnergySpecZ=zero
  end subroutine clcStat

  !******************************************************************
  ! clcEnergySpectraX
  !******************************************************************
  subroutine clcEnergySpectraX(arrx,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(xsize(1),xstart(2):xend(2),xsize(3)),intent(in)::arrx

    ! locals
    integer::ic,jc,kc,jt,icCnter
    real(RK)::normEgy,EgyX(nxhp)
    real(RK),dimension(xsize(1),njSpecX,xsize(3))::arrxSpec

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,njSpecX
      jt=jSpecX(jc); arrxSpec(:,jc,:)=arrx(:,jt,:)
    enddo
    call dfftw_execute_r2r(fft_plan_x,arrxSpec,arrxSpec)
    do jc=1,njSpecX
      EgyX=zero
      do kc=1,xsize(3)
        EgyX(1)   = EgyX(1)   +arrxSpec(1,jc,kc)   *arrxSpec(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+arrxSpec(nxhp,jc,kc)*arrxSpec(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) +(arrxSpec(ic,jc,kc)*arrxSpec(ic,jc,kc)+arrxSpec(icCnter,jc,kc)*arrxSpec(icCnter,jc,kc))*two
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(m,jc,ic)= EnergySpecX(m,jc,ic) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraX

  !******************************************************************
  ! clcEnergySpectraZ
  !******************************************************************
  subroutine clcEnergySpectraZ(arrz,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(zsize(1),zstart(2):zend(2),zsize(3)),intent(in)::arrz

    ! locals
    integer::ic,jc,kc,jt,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)
    real(RK),dimension(zsize(1),njSpecZ,zsize(3))::arrzSpec

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,njSpecZ
      jt=jSpecZ(jc); arrzSpec(:,jc,:)=arrz(:,jt,:)
    enddo
    call dfftw_execute_r2r(fft_plan_z,arrzSpec,arrzSpec)
    do jc=1,njSpecZ
      EgyZ=zero
      do ic=1,zsize(1)
        EgyZ(1)   = EgyZ(1)   +arrzSpec(ic,jc,1)   *arrzSpec(ic,jc,1)
        EgyZ(nzhp)= EgyZ(nzhp)+arrzSpec(ic,jc,nzhp)*arrzSpec(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) +(arrzSpec(ic,jc,kc)*arrzSpec(ic,jc,kc)+arrzSpec(ic,jc,kcCnter)*arrzSpec(ic,jc,kcCnter))*two
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(m,jc,kc)= EnergySpecZ(m,jc,kc) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraZ
    
end module m_FlowCase
