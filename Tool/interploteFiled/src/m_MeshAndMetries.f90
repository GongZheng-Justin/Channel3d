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
    character(64)::chFile
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

    ! ypOld
    chFile = trim(OldMeshName); nUnitFile=d_readyMesh_unit
    open(unit=nUnitFile,file=chfile,status='old',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file:  "//trim(adjustl(chFile)))
    do j=1,nypOld
      read(nUnitFile,*)iTemp,ypOld(j)
    enddo
    close(nUnitFile)
    if(abs(ypOld(1))>1.0E-10_RK .or. abs(ypOld(nypOld)-yly)>1.0E-10_RK) then
      call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Old mesh wrong!!")
    else
      ypOld(1)=zero; ypOld(nypOld)=yly
    endif

    ! ypNew
    chFile = trim(NewMeshName); nUnitFile=d_readyMesh_unit
    open(unit=nUnitFile,file=chfile,status='old',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file:  "//trim(adjustl(chFile)))
    do j=1,nypNew
      read(nUnitFile,*)iTemp,ypNew(j)
    enddo
    close(nUnitFile)
    if(abs(ypNew(1))>1.0E-10_RK .or. abs(ypNew(nypNew)-yly)>1.0E-10_RK) then
      call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Old mesh wrong!!")
    else
      ypNew(1)=zero; ypNew(nypNew)=yly
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

end module m_MeshAndMetries
