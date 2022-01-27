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
#if defined CFDDEM || defined CFDACM
  real(RK),public,save:: dyUniform,rdyUniform
#endif

  real(RK),public,save,allocatable,dimension(:):: yp   ! point coordinate in y-dir. Suffix 'v' means 'vector'
  real(RK),public,save,allocatable,dimension(:):: xc   ! center coordinate in x-dir    
  real(RK),public,save,allocatable,dimension(:):: yc   ! center coordinate in y-dir 
  real(RK),public,save,allocatable,dimension(:):: zc   ! center coordinate in z-dir
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

  ! uy/uz Laplacian metries in x-dir (STAGGERED VARIABLE)
  real(RK),public,save,allocatable,dimension(:):: am1c
  real(RK),public,save,allocatable,dimension(:):: ac1c
  real(RK),public,save,allocatable,dimension(:):: ap1c 
  ! uy/uz Laplacian metries in x-dir (for Crank-Nicolson scheme purpose)
  real(RK),public,save,allocatable,dimension(:):: am1cForCN
  real(RK),public,save,allocatable,dimension(:):: ap1cForCN 
  
  ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)
  real(RK),public,save,allocatable,dimension(:):: am2c
  real(RK),public,save,allocatable,dimension(:):: ac2c
  real(RK),public,save,allocatable,dimension(:):: ap2c
  ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
  real(RK),public,save,allocatable,dimension(:):: am2cForCN
  real(RK),public,save,allocatable,dimension(:):: ap2cForCN  

  ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
  real(RK),public,save,allocatable,dimension(:):: am2p
  real(RK),public,save,allocatable,dimension(:):: ac2p  
  real(RK),public,save,allocatable,dimension(:):: ap2p
  
  ! ux/uy Laplacian metries in z-dir (STAGGERED VARIABLE)
  real(RK),public,save,allocatable,dimension(:):: am3c
  real(RK),public,save,allocatable,dimension(:):: ac3c
  real(RK),public,save,allocatable,dimension(:):: ap3c
  ! ux/uy Laplacian metries in z-dir (for Crank-Nicolson scheme purpose)
  real(RK),public,save,allocatable,dimension(:):: am3cForCN
  real(RK),public,save,allocatable,dimension(:):: ap3cForCN 

  ! for linear interpolation in y-dir
  real(RK),public,save,allocatable,dimension(:):: YinterpCoe
  
  public:: InitMeshAndMetries
contains

  !******************************************************************
  ! InitMeshAndMetries
  !****************************************************************** 
  subroutine InitMeshAndMetries()
    implicit none 
    
    ! locals
    character(64) :: chFile
    real(RK):: tstr,xi,y0,ylyt,ystart  
    integer:: j,nypt,iErrSum,nUnitFile,myistat
    integer:: iErr01,iErr02,iErr03,iErr04,iErr05,iErr06,iErr07,iErr08,iErr09,iErr10,iErr11
    integer:: iErr12,iErr13,iErr14,iErr15,iErr16,iErr17,iErr18,iErr19,iErr20,iErr21,iErr22
    integer:: iErr23,iErr24,iErr25,iErr26,iErr27,iErr28,iErr29,iErr30,iErr31,iErr32

    dx = xlx/real(nxc,kind=RK)
    dy = yly/real(nyc,kind=RK)
    dz = zlz/real(nzc,kind=RK)
    rdx= real(nxc,kind=RK)/xlx
    rdy= real(nyc,kind=RK)/yly
    rdz= real(nzc,kind=RK)/zlz 
    dx2= dx*dx;  rdx2= rdx*rdx
    dy2= dy*dy;  rdy2= rdy*rdy
    dz2= dz*dz;  rdz2= rdz*rdz
    
    allocate(yp(0:nyp),   Stat = iErr01)
    allocate(yc(0:nyp),   Stat = iErr02)
    allocate(dyp(0:nyp),  Stat = iErr03)
    allocate(rdyp(0:nyp), Stat = iErr04)
    allocate(dyc(0:nyp),  Stat = iErr05)
    allocate(rdyc(0:nyp), Stat = iErr06)    
    
    allocate(ap2ph(1:nyc), Stat=iErr07)
    allocate(ac2ph(1:nyc), Stat=iErr08)    
    allocate(am2ph(1:nyc), Stat=iErr09)
    
    allocate(ap1c(1:nxc), Stat=iErr10)
    allocate(ac1c(1:nxc), Stat=iErr11)
    allocate(am1c(1:nxc), Stat=iErr12)
    allocate(ap1cForCN(1:nxc), Stat=iErr13)
    allocate(am1cForCN(1:nxc), Stat=iErr14)  
    
    allocate(ap2c(1:nyc), Stat=iErr15)
    allocate(ac2c(1:nyc), Stat=iErr16)
    allocate(am2c(1:nyc), Stat=iErr17)
    allocate(ap2cForCN(1:nyc), Stat=iErr18)
    allocate(am2cForCN(1:nyc), Stat=iErr19)
    
    allocate(ap2p(1:nyc), Stat=iErr20)
    allocate(ac2p(1:nyc), Stat=iErr21)    
    allocate(am2p(1:nyc), Stat=iErr22)
    
    allocate(ap3c(1:nzc), Stat=iErr23)
    allocate(ac3c(1:nzc), Stat=iErr24)
    allocate(am3c(1:nzc), Stat=iErr25)
    allocate(ap3cForCN(1:nzc), Stat=iErr26)
    allocate(am3cForCN(1:nzc), Stat=iErr27)

    allocate(xc(0:nxp),        Stat = iErr28)
    allocate(zc(0:nzp),        Stat = iErr29)
    allocate(VolCell(0:nyp),   Stat = iErr30)
    allocate(DeltaCell(0:nyp), Stat = iErr31)

    allocate(YinterpCoe(0:nyc), Stat = iErr32)
    iErrSum=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)+abs(iErr05)+abs(iErr06)+abs(iErr07)+abs(iErr08)+  &
            abs(iErr09)+abs(iErr10)+abs(iErr11)+abs(iErr12)+abs(iErr13)+abs(iErr14)+abs(iErr15)+abs(iErr16)+  &
            abs(iErr17)+abs(iErr18)+abs(iErr19)+abs(iErr20)+abs(iErr21)+abs(iErr22)+abs(iErr23)+abs(iErr24)+  &
            abs(iErr25)+abs(iErr26)+abs(iErr27)+abs(iErr28)+abs(iErr29)+abs(iErr30)+abs(iErr31)+abs(iErr32)
    if(iErrSum/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Allocation failed")
    
    ! yp, point coordinate in y-dir
    if(BcOption(ym_dir)==BC_PERIOD .and. istret /= 0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Bc=Peroidic in y-dir, while mesh is also stretching, WRONG!")
    endif

    ! Here:
    !   istret = 0, No clustering;
    !   istret = 1, Clustered in the center;
    !   istret = 2, Two walls clustering;
    !   istret = 3, Bottom wall clustering.
    ! Besides:
    !   If nyUniform>0 and istret /= 0, the first nyUniform grids in [0,yUniform] will be uniform.
    if(nyUniform<0)nyUniform=0
    if(istret<0 .or. istret>3)       call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries"," istret WRONG!")
    if(istret/=0 .and. nrank==0 .and. nyUniform>nyc) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","nyUniform>nyc WRONG!")
    if(istret/=0 .and. nrank==0 .and. yUniform>yly ) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries"," yUniform>yly WRONG!")
#if defined CFDDEM || defined CFDACM
    dyUniform=dy; rdyUniform=rdy
    if(istret/=0 .and. nyUniform>0) then
      dyUniform = yUniform/real(nyUniform,kind=RK)
      rdyUniform= real(nyUniform,kind=RK)/yUniform
    endif
#endif

    if(istret==0) then
      do j=1,nyp
        yp(j)=dy*real(j-1,kind=RK)
      enddo
    else
      if(nyUniform>0) then
        do j=1,nyUniform
          yp(j)=(yUniform/real(nyUniform,RK))*real(j-1,RK)
        enddo
        ystart=yUniform
        ylyt=yly-yUniform
        nypt=nyp-nyUniform
      else
        ystart=zero
        ylyt=yly
        nypt=nyp
      endif

      SELECT CASE(istret)
      CASE(1)
        tstr= tanh(cStret)
        do j= 1,nypt/2
          xi= real(j-1,kind=RK)/real(nypt-1,RK)       ! For j: [1,nypt/2],  xi: [0,1/2]
          y0= half*tanh(two*cStret*xi)/tstr           ! For j: [1,nypt/2],  y0: [0,1/2]
          yp(nyUniform+j)       = ystart +ylyt*y0
          yp(nyUniform+nypt+1-j)= ystart +ylyt*(one-y0)
        enddo
        if(mod(nypt,2)==1)yp(nyUniform+nypt/2+1)=half*ylyt+ystart
      CASE(2)
        tstr= tanh(cStret*half)
        do j= 1,nypt/2
          xi= real(j-1,kind=RK)/real(nypt-1,RK)       ! For j: [1,nypt/2],  xi: [0,1/2]
          y0= half*tanh(cStret*(xi-half))/tstr+ half  ! For j: [1,nypt/2],  y0: [0,1/2]
          yp(nyUniform+j)       = ystart +ylyt*y0
          yp(nyUniform+nypt+1-j)= ystart +ylyt*(one-y0)
        enddo
        if (mod(nypt,2)==1)yp(nyUniform+nypt/2+1)=half*ylyt+ystart
      CASE(3)
        tstr= tanh(cStret)
        do j= 1,nypt
          xi= real(j-1,kind=RK)/real(nypt-1)          ! For j: [1,nypt],    xi: [0,1]
          y0= tanh(cStret*(xi-one))/tstr + one        ! For j: [1,nypt],    y0: [0,1]
          yp(nyUniform+j)       = ystart +ylyt*y0
        enddo
      ENDSELECT
    endif
    yp(1)=zero; yp(nyp)=yly

    write(chFile,"(A)") trim(Res_Dir)//"yMeshFor"//trim(RunName)//".txt"
    nUnitFile=GetFileUnit()
    if(nrank==0) then
      open(unit=nUnitFile, file=chfile,status='replace',form='formatted',IOSTAT=myistat)
      if(myistat/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file: "//trim(chFile))
      do j=1,nyp
        write(nUnitFile,*)j,yp(j)
      enddo
      close(nUnitFile,IOSTAT=myistat)
    endif

    ! yc, center coordinate interval in y-dir
    do j=1,nyc
      yc(j) = half*(yp(j)+yp(j+1))
    enddo
    yc(0)=two*yp(1)-yc(1)
    yc(nyp)=two*yp(nyp)-yc(nyc)

    ! xc,zc, center coordinate interval in x-dir and z-dir
    do j=0,nxp
      xc(j)=dx*(real(j,RK)-half)
    enddo
    do j=0,nzp
      zc(j)=dz*(real(j,RK)-half)
    enddo

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
      YinterpCoe(j)= dyp(j)/(dyp(j)+dyp(j-1))
    enddo
    
    ! Pressure Laplacian metries in y-dir
    do j=1,nyc
      am2ph(j)= rdyp(j)*rdyc(j)
      ap2ph(j)= rdyp(j)*rdyc(j+1)
    enddo
    if(BcOption(ym_dir)==BC_NSLIP .or. BcOption(ym_dir)==BC_FSLIP) then
      am2ph(1)= zero
      ap2ph(1)= rdyp(1)*rdyc(2) 
    endif
    if(BcOption(yp_dir)==BC_NSLIP .or. BcOption(yp_dir)==BC_FSLIP) then
      am2ph(nyc)= rdyp(nyc)*rdyc(nyc)
      ap2ph(nyc)= zero 
    endif    
    ac2ph = -(am2ph+ap2ph)
    
    ! uy/uz Laplacian metries in x-dir (STAGGERED VARIABLE)=====================
    am1c = rdx2;  ap1c=rdx2
    if(BcOption(xm_dir)==BC_NSLIP ) then
      am1c(1)= four/three*rdx2
      ap1c(1)= four/three*rdx2 
    elseif(BcOption(xm_dir)==BC_FSLIP ) then
      am1c(1)= zero
      ap1c(1)= rdx2
    endif    
    if(BcOption(xp_dir)==BC_NSLIP ) then
      am1c(nxc)= four/three*rdx2
      ap1c(nxc)= four/three*rdx2 
    elseif(BcOption(xp_dir)==BC_FSLIP ) then
      am1c(nxc)= rdx2
      ap1c(nxc)= zero
    endif 
    ac1c= -(am1c+ap1c)

    ! uy/uz Laplacian metries in x-dir (for Crank-Nicolson scheme purpose)
    ap1cForCN = ap1c;   am1cForCN = am1c;
    if(BcOption(xm_dir)==BC_NSLIP) then
      ap1cForCN(1)= ap1c(1)
      am1cForCN(1)= two*am1c(1)    
    endif
    if(BcOption(xp_dir)==BC_NSLIP) then
      am1cForCN(nxc)= am1c(nxc)
      ap1cForCN(nxc)= two*ap1c(nxc)
    endif
    
    ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)=====================
    do j=1,nyc
      am2c(j)= rdyp(j)*rdyc(j)
      ap2c(j)= rdyp(j)*rdyc(j+1)
    enddo
    if(BcOption(ym_dir)==BC_NSLIP ) then
      am2c(1)= four*rdyc(1)/( dyc(1)+two*dyc(2) )
      ap2c(1)= four*rdyc(2)/( dyc(1)+two*dyc(2) )
    elseif(BcOption(ym_dir)==BC_FSLIP) then
      am2c(1)= zero
      ap2c(1)= rdyp(1)*rdyc(2)
    endif
    if(BcOption(yp_dir)==BC_NSLIP ) then
      am2c(nyc)= four*rdyc(nyc)/( dyc(nyp)+two*dyc(nyc) )
      ap2c(nyc)= four*rdyc(nyp)/( dyc(nyp)+two*dyc(nyc) )
    elseif(BcOption(yp_dir)==BC_FSLIP) then
      am2c(nyc)= rdyp(nyc)*rdyc(nyc)
      ap2c(nyc)= zero
    endif
    ac2c= -(am2c+ap2c)
    
    ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
    ap2cForCN = ap2c;   am2cForCN = am2c;
    if(BcOption(ym_dir)==BC_NSLIP) then
      ap2cForCN(1)= ap2c(1)
      am2cForCN(1)= two*am2c(1)    
    endif
    if(BcOption(yp_dir)==BC_NSLIP) then
      am2cForCN(nyc)= am2c(nyc)
      ap2cForCN(nyc)= two*ap2c(nyc)
    endif
    
    ! ux/uy Laplacian metries in z-dir (STAGGERED VARIABLE)=====================
    am3c = rdz2;  ap3c=rdz2;
    if(BcOption(zm_dir)==BC_NSLIP ) then
      am3c(1)= four/three*rdz2
      ap3c(1)= four/three*rdz2 
    elseif(BcOption(zm_dir)==BC_FSLIP ) then
      am3c(1)= zero
      ap3c(1)= rdz2
    endif
    if(BcOption(zp_dir)==BC_NSLIP ) then
      am3c(nzc)= four/three*rdz2
      ap3c(nzc)= four/three*rdz2 
    elseif(BcOption(zp_dir)==BC_FSLIP) then
      am3c(nzc)= rdz2
      ap3c(nzc)= zero
    endif
    ac3c= -(am3c+ap3c)

    ! ux/uy Laplacian metries in z-dir (for Crank-Nicolson scheme purpose)
    ap3cForCN = ap3c;   am3cForCN = am3c;
    if(BcOption(zm_dir)==BC_NSLIP) then
      ap3cForCN(1)= ap3c(1)
      am3cForCN(1)= two*am3c(1)    
    endif
    if(BcOption(zp_dir)==BC_NSLIP) then
      am3cForCN(nzc)= am3c(nzc)
      ap3cForCN(nzc)= two*ap3c(nzc)
    endif

    ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
    do j=1,nyc
      am2p(j)= rdyc(j)*rdyp(j-1)        
      ap2p(j)= rdyc(j)*rdyp(j)        
    enddo
    ac2p= -(am2p+ap2p)
      
  end subroutine InitMeshAndMetries

end module m_MeshAndMetries
