module m_MeshAndMetries
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use decomp_2d
  implicit none
  private 
  
  real(RK),public,save::dx,  dy,  dz   ! average mesh intervals in three directions
  real(RK),public,save::dx2, dy2, dz2  ! square of average mesh intervals in three directions  
  real(RK),public,save::rdx, rdy, rdz  ! inverse average mesh intervals in three directions
  real(RK),public,save::rdx2,rdy2,rdz2 ! square of the inverse average mesh intervals in three directions

  real(RK),public,save,allocatable,dimension(:):: yp   ! point coordinate in y-dir. Suffix 'v' means 'vector'   
  real(RK),public,save,allocatable,dimension(:):: yc   ! center coordinate in y-dir 
  real(RK),public,save,allocatable,dimension(:):: VolCell   ! cell volume
  real(RK),public,save,allocatable,dimension(:):: DeltaCell ! (dx*dy*dz)^(1/3)

  real(RK),public,save,allocatable,dimension(:):: dyp  ! point coordinate interval in y-dir
  real(RK),public,save,allocatable,dimension(:):: dyc  ! center coordinate interval in y-dir
  real(RK),public,save,allocatable,dimension(:):: rdyp ! inverse point coordinate interval in y-dir
  real(RK),public,save,allocatable,dimension(:):: rdyc ! inverse center coordinate interval in y-dir

  ! Pressure Laplacian metries in y-dir
  real(RK),public,save,allocatable,dimension(:):: ap2ph
  real(RK),public,save,allocatable,dimension(:):: ac2ph  
  real(RK),public,save,allocatable,dimension(:):: am2ph 
  
  ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)
  real(RK),public,save,allocatable,dimension(:):: am2c
  real(RK),public,save,allocatable,dimension(:):: ac2c
  real(RK),public,save,allocatable,dimension(:):: ap2c
  ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
  real(RK),public,save,allocatable,dimension(:):: am2cForCN
  real(RK),public,save,allocatable,dimension(:):: ac2cForCN
  real(RK),public,save,allocatable,dimension(:):: ap2cForCN  

  ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
  real(RK),public,save,allocatable,dimension(:):: am2p
  real(RK),public,save,allocatable,dimension(:):: ac2p  
  real(RK),public,save,allocatable,dimension(:):: ap2p

  ! for linear interpolation in y-dir
  real(RK),public,save,allocatable,dimension(:):: YinterpCoe

  ! for 4-order FD
  real(RK),public,save::dxCoe1,dxCoe2,dxCoe3,dxCoe4
  real(RK),public,save::dzCoe1,dzCoe2,dzCoe3,dzCoe4
  real(RK),public,save::dxxCoe1,dxxCoe2,dxxCoe3,dxxCoe4,dxxCoe5
  real(RK),public,save::dzzCoe1,dzzCoe2,dzzCoe3,dzzCoe4,dzzCoe5
  real(RK),public,save::InterpCoe1,InterpCoe2,InterpCoe3,InterpCoe4
  
  public:: InitMeshAndMetries
contains

  !******************************************************************
  ! InitMeshAndMetries
  !****************************************************************** 
  subroutine InitMeshAndMetries()
    implicit none 
    
    ! locals
    real(RK):: tstr,xi,y0
    character(64) :: chFile
    integer:: j,nypt,iErrSum,nUnitFile,myistat
    integer:: iErr01, iErr02, iErr03, iErr04, iErr05, iErr06, iErr07, iErr08, iErr09, iErr10 
    integer:: iErr11, iErr12, iErr13, iErr14, iErr15, iErr16, iErr17, iErr18, iErr19, iErr20, iErr21
    
    dx = xlx/real(nxc,kind=RK)
    dy = yly/real(nyc,kind=RK)
    dz = zlz/real(nzc,kind=RK)
    rdx= real(nxc,kind=RK)/xlx
    rdy= real(nyc,kind=RK)/yly
    rdz= real(nzc,kind=RK)/zlz 
    dx2= dx*dx;  rdx2= rdx*rdx
    dy2= dy*dy;  rdy2= rdy*rdy
    dz2= dz*dz;  rdz2= rdz*rdz

    ! for 4-order FD 
    dxCoe1=  rdx/24.00_RK; dxCoe2= -rdx*1.125_RK; dxCoe3=  rdx*1.125_RK; dxCoe4= -rdx/24.00_RK;
    dzCoe1=  rdz/24.00_RK; dzCoe2= -rdz*1.125_RK; dzCoe3=  rdz*1.125_RK; dzCoe4= -rdz/24.00_RK
    dxxCoe1=-rdx2/12.0_RK; dxxCoe2=4.0_RK/3.0_RK*rdx2; dxxCoe3=-2.5_RK*rdx2; dxxCoe4=4.0_RK/3.0_RK*rdx2; dxxCoe5=-rdx2/12.0_RK;
    dzzCoe1=-rdz2/12.0_RK; dzzCoe2=4.0_RK/3.0_RK*rdz2; dzzCoe3=-2.5_RK*rdz2; dzzCoe4=4.0_RK/3.0_RK*rdz2; dzzCoe5=-rdz2/12.0_RK;
    InterpCoe1= -1.0_RK/16.00_RK; InterpCoe2=  9.0_RK/16.00_RK; InterpCoe3=  9.0_RK/16.00_RK; InterpCoe4= -1.0_RK/16.00_RK
    
    allocate(yp(0:nyp),   Stat = iErr01)
    allocate(yc(0:nyp),   Stat = iErr02)
    allocate(dyp(0:nyp),  Stat = iErr03)
    allocate(rdyp(0:nyp), Stat = iErr04)
    allocate(dyc(0:nyp),  Stat = iErr05)
    allocate(rdyc(0:nyp), Stat = iErr06)    
    
    allocate(ap2ph(1:nyc), Stat=iErr07)
    allocate(ac2ph(1:nyc), Stat=iErr08)    
    allocate(am2ph(1:nyc), Stat=iErr09)
    
    allocate(ap2c(1:nyc), Stat=iErr10)
    allocate(ac2c(1:nyc), Stat=iErr11)
    allocate(am2c(1:nyc), Stat=iErr12)
    allocate(ap2cForCN(1:nyc), Stat=iErr13)
    allocate(ac2cForCN(1:nyc), Stat=iErr14)
    allocate(am2cForCN(1:nyc), Stat=iErr15)
    
    allocate(ap2p(1:nyc), Stat=iErr16)
    allocate(ac2p(1:nyc), Stat=iErr17)    
    allocate(am2p(1:nyc), Stat=iErr18)

    allocate(VolCell(0:nyp),   Stat = iErr19)
    allocate(DeltaCell(0:nyp), Stat = iErr20)
    allocate(YinterpCoe(0:nyc), Stat = iErr21)

    iErrSum=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)+abs(iErr05)+abs(iErr06)+abs(iErr07)+  &
            abs(iErr08)+abs(iErr09)+abs(iErr10)+abs(iErr11)+abs(iErr12)+abs(iErr13)+abs(iErr14)+  &
            abs(iErr15)+abs(iErr16)+abs(iErr17)+abs(iErr18)+abs(iErr19)+abs(iErr20)+abs(iErr21)
    if(iErrSum/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries: ","Allocation failed")

    SELECT CASE(istret)
    CASE (0) ! No clustering
      do j=1,nyp
        yp(j)=dy*real(j-1,kind=RK)
      enddo
    CASE(1) ! Two walls clustering
      tstr=sin(half*cStret*PI)
      do j=1,nyp/2
        xi= two*real(j-1,kind=RK)/real(nyc,kind=RK)-one  ! For j: [1,nyp/2],  xi: [-1,0]
        y0= sin(cStret*xi*PI*half)/tstr+ one             ! For j: [1,nyp/2],  y0: [0, 1]
        yp(j) = y0*yly*half
        yp(nyp+1-j)=yly-yp(j)
      enddo
      yp(1)=zero; yp(nyp)=yly
      if (mod(nyp,2)==1)yp(nyp/2+1)=half*yly
    CASE(2) ! Bottom wall clustering
      tstr=sin(half*cStret*PI)
      do j=1,nyp
        xi= real(j-1,kind=RK)/real(nyc,kind=RK)-one      ! For j: [1,nyp],  xi: [-1,0]
        y0= sin(cStret*xi*PI*half)/tstr + one            ! For j: [1,nyp],  y0: [0, 1]
        yp(j) = y0*yly
      enddo
      yp(1)=zero; yp(nyp)=yly
    END SELECT

    chFile = trim(Res_Dir)//"yMeshFor"//trim(RunName)//".txt"
    nUnitFile = d_writeyMesh_unit
    if(nrank==0) then
      open(unit=nUnitFile, file=chfile,status='replace',form='formatted',IOSTAT=myistat)
      if(myistat /= 0) then
        call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file:  "//trim(adjustl(chFile)))
      endif
      do j=1,nyp
        write(nUnitFile,*)j,yp(j)
      enddo
      close(nUnitFile)
    endif

    ! yc, center coordinate interval in y-dir
    do j=1,nyc
      yc(j) = half*(yp(j)+yp(j+1))
    enddo
    yc(0)=two*yp(1)-yc(1)
    yc(nyp)=two*yp(nyp)-yc(nyc)

    ! dyp, point coordinate interval in y-dir 
    do j=1,nyc
      dyp(j) = yp(j+1) - yp(j)    
    enddo
    dyp(0)  =dyp(1)
    dyp(nyp)=dyp(nyc)

    !VolCell,volume of the cell
    do j=0,nyp
      VolCell(j)=dyp(j)*dx*dz
    enddo
  
    !DeltaCell, (VolCell)^(1/3)
    do j=0,nyp
      DeltaCell(j)=(dyp(j)*dx*dz)**(one/three)
    enddo

    ! dyc, center coordinate interval in y-dir
    do j=2,nyc
      dyc(j) = yc(j)-yc(j-1)
    enddo
    dyc(1)=dyp(1)
    dyc(nyp)=dyp(nyc)
    
    ! rdyp and rdyc, the reverse of the dyp and dyc, respectively.
    do j=0, nyp
      rdyp(j) = one/dyp(j)    
    enddo    
    do j=1,nyp
      rdyc(j) = one/dyc(j)    
    enddo

    ! for linear interpolation in y-dir
    do j=0,nyc
      YinterpCoe(j)= dyp(j)/(dyp(j)+dyp(j+1))
    enddo
    
    ! Pressure Laplacian metries in y-dir
    do j=1,nyc
      am2ph(j)= rdyp(j)*rdyc(j)
      ap2ph(j)= rdyp(j)*rdyc(j+1)
    enddo
    am2ph(1)= zero
    ap2ph(1)= rdyp(1)*rdyc(2) 
    am2ph(nyc)= rdyp(nyc)*rdyc(nyc)
    ap2ph(nyc)= zero 
    ac2ph = -(am2ph+ap2ph)
    
    ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)=====================
    do j=1,nyc
      am2c(j)= rdyp(j)*rdyc(j)
      ap2c(j)= rdyp(j)*rdyc(j+1)
    enddo
    am2c(1)= four*rdyc(1)/( dyc(1)+two*dyc(2) )
    ap2c(1)= four*rdyc(2)/( dyc(1)+two*dyc(2) )
    if(FlowType==FT_CH ) then
      am2c(nyc)= four*rdyc(nyc)/( dyc(nyp)+two*dyc(nyc) )
      ap2c(nyc)= four*rdyc(nyp)/( dyc(nyp)+two*dyc(nyc) )
    elseif(FlowType==FT_HC) then
      am2c(nyc)= rdyp(nyc)*rdyc(nyc)
      ap2c(nyc)= zero
    endif
    ac2c= -(am2c+ap2c)
    
    ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
    ap2cForCN = ap2c;  ac2cForCN = ac2c;  am2cForCN = am2c;
    ap2cForCN(1)= ap2c(1)    
    ac2cForCN(1)= ac2c(1)-am2c(1)
    am2cForCN(1)= two*am2c(1)
    if(FlowType==FT_CH) then
      am2cForCN(nyc)= am2c(nyc)
      ac2cForCN(nyc)= ac2c(nyc)-ap2c(nyc)
      ap2cForCN(nyc)= two*ap2c(nyc)
    endif

    ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
    do j=1,nyc
      am2p(j)= rdyc(j)*rdyp(j-1)        
      ap2p(j)= rdyc(j)*rdyp(j)        
    enddo
    ac2p= -(am2p+ap2p)
      
  end subroutine InitMeshAndMetries

end module m_MeshAndMetries
