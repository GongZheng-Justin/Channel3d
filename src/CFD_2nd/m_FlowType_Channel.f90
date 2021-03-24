module m_FlowType_Channel
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_MeshAndMetries
  use decomp_2d
  use m_Variables,only: mb1
  use m_Tools,only:CalcUxAver
  implicit none
  private

  ! statistics variabls
  integer,save:: nfstime
  real(RK),save:: PrGradsum
  real(RK),save,allocatable,dimension(:,:):: SumStat 

  public:: InitVelocity_CH,Update_uy_ym_CH,InitStatVar_CH,clcStat_CH
contains

  !******************************************************************
  ! InitVelocity_CH
  !******************************************************************
  subroutine InitVelocity_CH(ux,uy,uz,Deviation)
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
	    rem = ubulk*(half*yly)/xnu
    case(FT_HC)
      height=yly
	    rem = ubulk*yly/xnu
    endselect

    ux=zero; uy=zero; uz=zero
    retau_guass = 0.1538_RK*rem**0.887741_RK
	  utau_guass 	= retau_guass*xnu/height
    if(nrank==0)print*,'************** retau_gauss=',retau_guass
    if(nrank==0)print*,'************** utau_gauss= ',utau_guass

    call system_clock(count=code); !code=0
    call random_seed(size = ii)
    call random_seed(put = code+63946*(/(i-1,i=1,ii)/))
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
    ux= ux*ratiot
    do j=ystart(2),yend(2)
      uz(:,j,:)=uz(:,j,:)-uzmeanR(j)
    enddo
    
  end subroutine InitVelocity_CH

  !******************************************************************
  ! Update_uy_ym_CH
  !******************************************************************   
  subroutine Update_uy_ym_CH(uy_ym, duy_ym, TimeNew)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(out):: uy_ym
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(inout):: duy_ym    
    real(RK),intent(in):: TimeNew
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = zero
    duy_ym = uy_ym - duy_ym
  end subroutine Update_uy_ym_CH

  !******************************************************************
  ! InitStatVar_CH
  !******************************************************************
  subroutine InitStatVar_CH()
    implicit none

    ! locals
    integer::iErr
    character(len=50)::filename

    if(nrank==0) then
      if(mod(saveStat,ivstats)/=0 )    call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","ivstats wrong !!!")
      if(IsUxConst)then
        write(filename,'(A,I10.10)')trim(Res_dir)//"PrGrad",ilast
        open(79,file=filename,status='replace',form='formatted',IOSTAT=iErr)
        if(iErr /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","Cannot open file: "//trim(filename))
        close(79,IOSTAT=iErr)
      endif
    endif
    allocate(SumStat(35,nyp),Stat=iErr)
    if(iErr /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","Allocation failed")  
    nfstime=0;  SumStat=zero; PrGradsum=zero
  end subroutine InitStatVar_CH

  !******************************************************************
  ! clcStat_CH
  !******************************************************************
  subroutine clcStat_CH(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
   
    ! locals
    character(len=50)::filename,myformat
    real(RK)::dudxx,dudyy,dudzz,dvdxx,dvdyy,dvdzz,dwdxx,dwdyy,dwdzz,SumStatR(35,nyp),SumVec(35),rdxh
    real(RK)::uxloc,uyloc,uzloc,prloc,uxCell,uyCell,uzCell,prloc2,inxz,infstime,cac,caj,cacU,dudyU,dvdyM
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vor_x,vor_y,vor_z,InterpY1,InterpY2,dudzC,dvdzC,dudxC,dvdxU
    integer::ic,jc,kc,im,jm,km,ip,jp,kp,myistat,ierror,coord1,coord2,color,key,SpecWORLDX,SpecWORLDZ,nrankX,nrankZ

    rdxh=rdx*half
    inxz = one/(real(nxc,RK)*real(nzc,RK))
    DO jc=1,nyc
      jp=jc+1; jm=jc-1;
      InterpY1= half*YinterpCoe(jm); InterpY2=half-InterpY1
      cac=rdyc(jc);cacU=rdyc(jp); caj=rdyp(jc); SumVec=zero
      do kc=ystart(3),yend(3)
        km=kc-1;kp=kc+1
        do ic=ystart(1),yend(1)
          im=ic-1;ip=ic+1

          uxloc = ux(ic,jc,kc)
          uyloc = uy(ic,jc,kc)
          uzloc = uz(ic,jc,kc)
          prloc = pressure(ic,jc,kc)
          uxCell= half*(ux(ic,jc,kc)+ux(ip,jc,kc))
          uyCell= half*(uy(ic,jc,kc)+uy(ic,jp,kc))
          uzCell= half*(uz(ic,jc,kc)+uz(ic,jc,kp))
          prloc2= InterpY1*(pressure(im,jm,kc)+pressure(ic,jm,kc))+ InterpY2*(pressure(im,jc,kc)+pressure(ic,jc,kc))
 
          dudx= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
          dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
          dudz= (ux(ic,jc,kc)-ux(ic,jc,km))*rdz
          dvdx= (uy(ic,jc,kc)-uy(im,jc,kc))*rdx
          dvdy= (uy(ic,jp,kc)-uy(ic,jc,kc))*caj
          dvdz= (uy(ic,jc,kc)-uy(ic,jc,km))*rdz
          dwdx= (uz(ic,jc,kc)-uz(im,jc,kc))*rdx
          dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
          dwdz= (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz

          dudxx= (ux(ip,jc,kc)-two*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          dudyy= ap2c(jc)*ux(ic,jp,kc)+ac2c(jc)*ux(ic,jc,kc)+am2c(jc)*ux(ic,jm,kc)
          dudzz= (ux(ic,jc,kp)-two*ux(ic,jc,kc)+ux(ic,jc,km))*rdz2
          dvdxx= (uy(ip,jc,kc)-two*uy(ic,jc,kc)+uy(im,jc,kc))*rdx2
          dvdyy= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          dvdzz= (uy(ic,jc,kp)-two*uy(ic,jc,kc)+uy(ic,jc,km))*rdz2
          dwdxx= (uz(ip,jc,kc)-two*uz(ic,jc,kc)+uz(im,jc,kc))*rdx2
          dwdyy= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)
          dwdzz= (uz(ic,jc,kp)-two*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          vor_x= dwdy-dvdz
          vor_y= dudz-dwdx
          vor_z= dvdx-dudy

          dudxC= (ux(ip,jc,kc)-ux(im,jc,kc))*rdxh
          dvdxU= (uy(ic,jp,kc)-uy(im,jp,kc))*rdx
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
      write(filename,'(A,I10.10)')trim(Res_dir)//"PrGrad",ilast
      open(79,file=filename,status='old',position='append',form='formatted',IOSTAT=myistat)
      if(myistat/=0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      else
        write(79,'(I7,2ES24.15)')itime,SimTime,PrGradAver
      endif
      close(79,IOSTAT=myistat)
    endif
    nfstime= nfstime + 1
    if(mod(itime,SaveStat)/=0) return

    ! Write statistics
    call MPI_REDUCE(SumStat,SumStatR,35*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = one/real(nfstime,RK)
      write(filename,"(A,I10.10)") trim(Res_dir)//'stats',itime
      open(69,file=filename,status='replace',form='formatted',IOSTAT=myistat)
      IF(myistat/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
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

    nfstime=0; SumStat=zero; PrGradsum=zero;
  end subroutine clcStat_CH

end module m_FlowType_Channel
