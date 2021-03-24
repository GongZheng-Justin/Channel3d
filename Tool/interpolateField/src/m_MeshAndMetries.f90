module m_MeshAndMetries
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  implicit none
  private 
  
  real(RK),public,save,allocatable,dimension(:):: xpOld,xpNew,xcOld,xcNew
  real(RK),public,save,allocatable,dimension(:):: ypOld,ypNew,ycOld,ycNew
  real(RK),public,save,allocatable,dimension(:):: zpOld,zpNew,zcOld,zcNew
  public:: InitMeshAndMetries
contains

  !******************************************************************
  ! InitMeshAndMetries
  !****************************************************************** 
  subroutine InitMeshAndMetries()
    implicit none 
    
    ! locals
    character(64) :: chFile
    integer:: j,iErr01,iErr02,iErr03,iErr04,ierror,iTemp,nUnitFile

    ! x mesh =============
    allocate(xpOld(1:nxpOld), Stat=iErr01)
    allocate(xpNew(1:nxpNew), Stat=iErr02)
    allocate(xcOld(0:nxpOld), Stat=iErr03)
    allocate(xcNew(0:nxpNew), Stat=iErr04)
    ierror=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries: ","Allocation failed: x-dir")
    do j=1,nxpOld
      xpOld(j)=real(j-1,RK)/real(nxcOld,RK)
    enddo
    do j=1,nxpNew
      xpNew(j)=real(j-1,RK)/real(nxcNew,RK)
    enddo
    xpOld(1)=zero; xpOld(nxpOld)=one
    xpNew(1)=zero; xpNew(nxpNew)=one

    ! xc, center coordinate interval in x-dir
    do j=1,nxcOld
      xcOld(j) = half*(xpOld(j)+xpOld(j+1))
    enddo
    do j=1,nxcNew
      xcNew(j) = half*(xpNew(j)+xpNew(j+1))
    enddo
    xcOld(0)     =two*xpOld(1)     -xcOld(1)
    xcNew(0)     =two*xpNew(1)     -xcNew(1)
    xcOld(nxpOld)=two*xpOld(nxpOld)-xcOld(nxcOld)
    xcNew(nxpNew)=two*xpNew(nxpNew)-xcNew(nxcNew)

    ! z mesh =============
    allocate(zpOld(1:nzpOld), Stat=iErr01)
    allocate(zpNew(1:nzpNew), Stat=iErr02)
    allocate(zcOld(0:nzpOld), Stat=iErr03)
    allocate(zcNew(0:nzpNew), Stat=iErr04)
    ierror=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries: ","Allocation failed: z-dir")
    do j=1,nzpOld
      zpOld(j)=real(j-1,RK)/real(nzcOld,RK)
    enddo
    do j=1,nzpNew
      zpNew(j)=real(j-1,RK)/real(nzcNew,RK)
    enddo
    zpOld(1)=zero; zpOld(nzpOld)=one
    zpNew(1)=zero; zpNew(nzpNew)=one

    ! zc, center coordinate interval in z-dir
    do j=1,nzcOld
      zcOld(j) = half*(zpOld(j)+zpOld(j+1))
    enddo
    do j=1,nzcNew
      zcNew(j) = half*(zpNew(j)+zpNew(j+1))
    enddo
    zcOld(0)     =two*zpOld(1)     -zcOld(1)
    zcNew(0)     =two*zpNew(1)     -zcNew(1)
    zcOld(nzpOld)=two*zpOld(nzpOld)-zcOld(nzcOld)
    zcNew(nzpNew)=two*zpNew(nzpNew)-zcNew(nzcNew)

    ! y mesh =============
    allocate(ypOld(1:nypOld), Stat=iErr01)
    allocate(ypNew(1:nypNew), Stat=iErr02)
    allocate(ycOld(0:nypOld), Stat=iErr03)
    allocate(ycNew(0:nypNew), Stat=iErr04)
    ierror=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries: ","Allocation failed: y-dir")

    ! Mesh information
    call clcYp(ypOld,"old")
    call clcYp(ypNew,"new") 
    if(nrank==0) then
      nUnitFile=d_readyMesh_unit
      chFile = trim(Res_Dir)//"yMeshOldFor"//trim(RunName)//".txt"
      open(unit=nUnitFile,file=chFile,status="replace",form='formatted',IOSTAT=ierror)
      if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file:  "//trim(chFile))
      do j=1,nypOld
        write(nUnitFile,*)j,ypOld(j)
      enddo
      close(nUnitFile)
      chFile = trim(Res_Dir)//"yMeshNewFor"//trim(RunName)//".txt"
      open(unit=nUnitFile,file=chFile,status="replace",form='formatted',IOSTAT=ierror)
      if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file:  "//trim(chFile))
      do j=1,nypNew
        write(nUnitFile,*)j,ypNew(j)
      enddo
      close(nUnitFile)
    endif

    ! yc, center coordinate interval in y-dir
    do j=1,nycOld
      ycOld(j) = half*(ypOld(j)+ypOld(j+1))
    enddo
    do j=1,nycNew
      ycNew(j) = half*(ypNew(j)+ypNew(j+1))
    enddo
    ycOld(0)     =two*ypOld(1)     -ycOld(1)
    ycNew(0)     =two*ypNew(1)     -ycNew(1)
    ycOld(nypOld)=two*ypOld(nypOld)-ycOld(nycOld)
    ycNew(nypNew)=two*ypNew(nypNew)-ycNew(nycNew)
  end subroutine InitMeshAndMetries

  !******************************************************************
  ! clcYp
  !******************************************************************
  subroutine clcYp(yp,old_or_new)
    implicit none
    real(RK),dimension(:),intent(out)::yp
    character(len=3),intent(in)::old_or_new
    
    ! locals
    integer::j,nypt,nyp,nyc,StretType,istret,nyUniform
    real(RK)::dy,tstr,xi,y0,ylyt,ystart,cStret,yUniform

    if(old_or_new == "old") then     ! old
      cStret=cStretOld; istret=istretOld; yUniform=yUniformOld
      nyp=nypOld; nyUniform=nyUniformOld; StretType=StretTypeOld
    elseif(old_or_new == "new") then ! new
      cStret=cStretNew; istret=istretNew; yUniform=yUniformNew
      nyp=nypNew; nyUniform=nyUniformNew; StretType=StretTypeNew
    endif
    nyc=nyp-1;      dy=yly/real(nyc,RK)

    if(StretType==0) then
      if(nyUniform<0)nyUniform=0
      if(istret<0 .or. istret>3)       call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries"," istret WRONG!")
      if(istret/=0 .and. nrank==0 .and. nyUniform>nyc) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","nyUniform>nyc WRONG!")
      if(istret/=0 .and. nrank==0 .and. yUniform>yly ) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries"," yUniform>yly WRONG!")

      if(istret==0) then ! 0:tangent hyperbolic function
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

    elseif(StretType==1) then ! 1:sine/cosine function

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
    endif
  end subroutine clcYp

end module m_MeshAndMetries
