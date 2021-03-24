module m_Parameters
  use m_TypeDef
  use m_LogInfo
  implicit none  
  private
  
  ! file units for different routines
  integer,parameter,public:: d_read_Param_unit= 1001
  integer,parameter,public:: d_LogInfo_unit   = 1002
  integer,parameter,public:: d_readyMesh_unit = 1003

  integer,parameter,public:: xm_dir=1  ! x- direction
  integer,parameter,public:: xp_dir=2  ! x+ direction
  integer,parameter,public:: ym_dir=3  ! y- direction
  integer,parameter,public:: yp_dir=4  ! y+ direction
  integer,parameter,public:: zm_dir=5  ! z- direction
  integer,parameter,public:: zp_dir=6  ! z+ direction

  logical, public,save:: IsUxConst
  real(RK),public,save:: ubulk

  ! mesh options
  real(RK),public,save::yly
  integer,public,save::nxpOld,nypOld,nzpOld     ! Old grid point number in three directions
  integer,public,save::nxcOld,nycOld,nzcOld     ! Old grid center number in three directions. nxc = nxp-1
  integer,public,save::nxpNew,nypNew,nzpNew     ! New grid point number in three directions
  integer,public,save::nxcNew,nycNew,nzcNew     ! New grid center number in three directions. nxc = nxp-1
  integer,public,save::StretTypeOld,istretOld,nyUniformOld
  integer,public,save::StretTypeNew,istretNew,nyUniformNew
  real(RK),public,save::cStretOld,yUniformOld,cStretNew,yUniformNew
  
  ! Boundary conditions
  !  0 : periodic conditions for all variables
  ! -1: no slip for three velocity components and ZERO GRADIENT for pressure field
  ! -2: free slip for three velocity components and ZERO GRADIENT for pressure field
  integer,parameter,public::BC_PERIOD =  0  
  integer,parameter,public::BC_NSLIP  = -1
  integer,parameter,public::BC_FSLIP  = -2
  integer, public,dimension(6):: BcOption
  real(RK),public,dimension(6):: uxBcValue
  real(RK),public,dimension(6):: uyBcValue
  real(RK),public,dimension(6):: uzBcValue

  character(64),public,save:: RunName      ! Run name
  character(64),public,save:: Res_Dir      ! Result directory
  character(64),public,save:: OldRestartName
  character(64),public,save:: NewRestartName
  integer,public,save:: Cmd_LFile_Freq= 1  ! report frequency in the terminal 
  integer,public,save:: LF_file_lvl   = 5  ! logfile report level      
  integer,public,save:: LF_cmdw_lvl   = 3  ! terminal report level

  ! Decomp2d options
  integer,public,save:: p_row,p_col

  public:: ReadAndInitParameters,DumpReadedParam
contains
    
  !******************************************************************
  ! InitParameters
  !****************************************************************** 
  subroutine ReadAndInitParameters(chFile)
    implicit none 
    character(*),intent(in)::chFile
    
    ! locals
    integer:: nUnitFile,myistat
    NAMELIST/BasicParam/nxpOld,nypOld,nzpOld,nxpNew,nypNew,nzpNew,BcOption,uxBcValue,uyBcValue,uzBcValue,  &
                        RunName,Res_Dir,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,p_row,p_col,yly,IsUxConst,  &
                        ubulk,OldRestartName,NewRestartName,StretTypeOld,istretOld,cStretOld,nyUniformOld, &
                        yUniformOld,StretTypeNew,istretNew,cStretNew,nyUniformNew,yUniformNew

    nUnitFile = d_read_Param_unit 
    open(unit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=myistat )
    if(myistat/=0) call MainLog%CheckForError(ErrT_Abort,"ReadAndInitParameters","Open file wrong: "//trim(chFile))
    read(nUnitFile,nml=BasicParam)
    close(nUnitFile)  

    nxcOld =nxpOld-1; nycOld=nypOld-1; nzcOld=nzpOld-1
    nxcNew =nxpNew-1; nycNew=nypNew-1; nzcNew=nzpNew-1
    if(p_row*p_col /= nproc) call MainLog%CheckForError(ErrT_Abort,"ReadAndInitParameters","p_row*p_col /= nproc")
  end subroutine ReadAndInitParameters

  !******************************************************************
  ! DumpReadedParam
  !****************************************************************** 
  subroutine DumpReadedParam()
    implicit none

    ! locals
    NAMELIST/BasicParam/nxpOld,nypOld,nzpOld,nxpNew,nypNew,nzpNew,BcOption,uxBcValue,uyBcValue,uzBcValue,  &
                        RunName,Res_Dir,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,p_row,p_col,yly,IsUxConst,  &
                        ubulk,OldRestartName,NewRestartName,StretTypeOld,istretOld,cStretOld,nyUniformOld, &
                        yUniformOld,StretTypeNew,istretNew,cStretNew,nyUniformNew,yUniformNew

    write(d_LogInfo_unit, nml=BasicParam)
  end subroutine DumpReadedParam
  
end module m_Parameters
