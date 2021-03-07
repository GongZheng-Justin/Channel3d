module m_Parameters
  use m_TypeDef
  implicit none  
  private
  
  ! file units for different routines
  integer,parameter,public:: d_read_Param_unit= 1001
  integer,parameter,public:: d_LogInfo_unit   = 1002
  integer,parameter,public:: d_writeyMesh_unit= 1003
  integer,parameter,public:: d_write_xdmf_unit= 1004 
  
  ! Flow type option
  integer,parameter,public:: FT_CH=1  ! Channel
  integer,parameter,public:: FT_HC=2  ! Half Channel
  integer,parameter,public:: FT_TG=3  ! Taylor-Green vortex
  integer,parameter,public:: FT_HI=4  ! Homogenerous isotropic turbulence
  integer,parameter,public:: FT_AN=5  ! Added new
  integer, public,save:: FlowType
  real(RK),public,save:: ubulk
  logical,public,save :: IsUxConst
  real(RK),public,save:: PrGradAver
  
  ! mesh options
  real(RK),public,save::xlx, yly, zlz  ! domain length in three directions
  integer,public,save::nxp,nyp,nzp     ! grid point number in three directions
  integer,public,save::nxc,nyc,nzc     ! grid center number in three directions. nxc = nxp-1
  integer, public,save:: istret        ! Stretch y-mesh or not (0:no, 1:center, 2:both sides, 3:bottom)
  real(RK),public,save:: cStret        ! Stretching parameter, if istret=0. this parameter doesn't work. 
  integer,public,save::  nyUniform     ! If nyUniform>0 and istret /= 0, the first nyUniform grids near bottom will set to be uniform.
  real(RK),public,save:: yUniform      ! If nyUniform>0 and istret /= 0, the first nyUniform grids in [0,yUniform] will be uniform.
  
  integer,parameter,public:: x_pencil=1
  integer,parameter,public:: y_pencil=2
  integer,parameter,public:: z_pencil=3
  integer,parameter,public:: xm_dir=1  ! x- direction
  integer,parameter,public:: xp_dir=2  ! x+ direction
  integer,parameter,public:: ym_dir=3  ! y- direction
  integer,parameter,public:: yp_dir=4  ! y+ direction
  integer,parameter,public:: zm_dir=5  ! z- direction
  integer,parameter,public:: zp_dir=6  ! z+ direction  
  
  ! Physical properties
  real(RK),save,public:: xnu                    ! Kinematic viscosity
  real(RK),save,public:: FluidDensity           ! Fluid density 
  real(RK),dimension(3),save,public:: gravity   ! Gravity or  other constant body forces (if any)
  
  ! Time stepping scheme and Projection method options
  integer,parameter,public:: FI_AB2 =1          !  AB2 and treat the viscous term implicitly in y-dir by C-N
  integer,parameter,public:: FI_RK3 =2          !  RK3 and treat the viscous term implicitly in y-dir by C-N
  real(RK),save,public:: dt         ! current time step
  real(RK),save,public:: dtMax      ! Maxium time step
  real(RK),save,public:: SimTime    ! Real simulation time
  integer,save,public :: iCFL       ! Use CFL condition to change time step dynamically( 1: yes, 2:no ).
  real(RK),save,public:: CFLc       ! CFL parameter
  integer,save,public::  itime      ! current time step
  integer,save,public::  ifirst     ! First iteration
  integer,save,public::  ilast      ! Last iteration 
  integer,save,public::  iadvance
  integer,save,public::  ischeme
  integer,save,public::  IsImplicit  !(0=full explicit, 1=partial implicit, 2=full implicit )
  real(RK),dimension(3),save,public :: pmGammaConst
  real(RK),dimension(3),save,public :: pmThetaConst
  real(RK),dimension(3),save,public :: pmAlphaConst
  real(RK),save,public::pmGamma
  real(RK),save,public::pmTheta
  real(RK),save,public::pmAlpha
  real(RK),save,public::pmAlphaC
  real(RK),save,public::pmBeta
  real(RK),save,public::pmBetaT

  ! FFTW_option
  integer,save,public:: FFTW_plan_type
  
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

  ! I/O, Statistics
  integer,public,save:: ivstats            ! time step interval for statistics calculation 
  integer,public,save:: SaveVisu           ! Output visulizing file frequency
  integer,public,save:: BackupFreq         ! Output Restarting file frequency
  integer,public,save:: SaveStat           ! Output Statistics file frequency

  logical,public,save::       RestartFlag  ! restart or not
  character(64),public,save:: RunName      ! Run name
  character(64),public,save:: Res_Dir      ! Result directory
  character(64),public,save:: RestartDir   ! Restart directory
  integer,public,save:: Cmd_LFile_Freq= 1  ! report frequency in the terminal 
  integer,public,save:: LF_file_lvl   = 5  ! logfile report level      
  integer,public,save:: LF_cmdw_lvl   = 3  ! terminal report level

  ! Decomp2d options
  integer,public,save:: p_row,p_col

  ! limited velocity and div
  real(RK),public,save:: vel_limit
  real(RK),public,save:: div_limit

  ! LES options
  ! 0:none; 1:Smagorinsky Model; 2:constant Smagorinsky Model; 3:Dynamic Smagorinsky Model; 4:MTS Model;
  integer,public,save:: LES_type
  integer,public,save:: FilterType   ! 0:trapezoidal type, 1:Simpson type

  public:: ReadAndInitParameters,DumpReadedParam,PMcoeUpdate
contains
    
  !******************************************************************
  ! InitParameters
  !****************************************************************** 
  subroutine ReadAndInitParameters(chFile)
    implicit none 
    character(*),intent(in)::chFile
    
    ! locals
    integer:: nUnitFile,myistat
    NAMELIST/BasicParam/FlowType,ubulk,xlx,yly,zlz,nxp,nyp,nzp,istret,cStret,nyUniform,yUniform,xnu,      &
                        dtMax,iCFL,CFLc,ifirst,ilast,ischeme,IsImplicit,FFTW_plan_type, BcOption,gravity, &
                        uxBcValue, uyBcValue, uzBcValue, ivstats,BackupFreq,SaveStat, SaveVisu,IsUxConst, &
                        RestartFlag,RunName,Res_Dir,RestartDir,Cmd_LFile_Freq, LF_file_lvl,LF_cmdw_lvl,   &
                        p_row,p_col, vel_limit, div_limit, FluidDensity
    NAMELIST/LesOptions/LES_type, FilterType

    nUnitFile = d_read_Param_unit 
    open(unit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=myistat )
    if(myistat /= 0) then
      print*,"Cannot open file: "//trim(adjustl(chFile)); STOP
    endif
    read(nUnitFile, nml=BasicParam)
    read(nUnitFile, nml=LesOptions)
    close(nUnitFile)  

    nxc= nxp-1
    nyc= nyp-1
    nzc= nzp-1
    if(ischeme==FI_AB2) then
      iadvance = 1
      pmGammaConst = (/ three/two, zero, zero/)
      pmThetaConst = (/     -half, zero, zero/)
    elseif(ischeme==FI_RK3) then
      iadvance = 3
      pmGammaConst = (/ eight/fifteen,       five/twelve,    three/four/)
      pmThetaConst = (/          zero,  -seventeen/sixty,  -five/twelve/)
    endif
    pmAlphaConst = pmGammaConst + pmThetaConst   
    
  end subroutine ReadAndInitParameters

  !******************************************************************
  ! DumpReadedParam
  !****************************************************************** 
  subroutine DumpReadedParam()
    implicit none

    ! locals
    NAMELIST/BasicParam/FlowType,ubulk,xlx,yly,zlz,nxp,nyp,nzp,istret,cStret,nyUniform,yUniform,xnu,      &
                        dtMax,iCFL,CFLc,ifirst,ilast,ischeme,IsImplicit,FFTW_plan_type, BcOption,gravity, &
                        uxBcValue, uyBcValue, uzBcValue, ivstats,BackupFreq,SaveStat, SaveVisu,IsUxConst, &
                        RestartFlag,RunName,Res_Dir,RestartDir,Cmd_LFile_Freq, LF_file_lvl,LF_cmdw_lvl,   &
                        p_row,p_col, vel_limit, div_limit, FluidDensity
    NAMELIST/LesOptions/LES_type, FilterType

    write(d_LogInfo_unit, nml=BasicParam)
    write(d_LogInfo_unit, nml=LesOptions)
  end subroutine DumpReadedParam

  !******************************************************************
  ! PMcoeUpdate
  !******************************************************************
  subroutine PMcoeUpdate(ns)
    implicit none
    integer,intent(in)::ns
      
    pmGamma = pmGammaConst(ns) * dt
    pmTheta = pmThetaConst(ns) * dt
    pmAlpha = pmAlphaConst(ns) * dt
    pmBeta  = half * pmAlphaConst(ns) *xnu * dt
    SimTime=SimTime+dt*real(pmAlphaConst(ns),RK)

    pmAlphaC= pmAlphaConst(ns)              ! for pressure gradient purpose only
    pmBetaT = half * pmAlphaConst(ns) * dt  ! for LES purpose only
    if(ns==1) PrGradAver=zero
  end subroutine PMcoeUpdate
  
end module m_Parameters
