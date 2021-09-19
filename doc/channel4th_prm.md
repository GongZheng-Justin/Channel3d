# input parameter for `channel4th` 

&emsp;Input file for `channel4th` can be divided into 3 part. Every part contains one **namelist**, i.e. `BasicParam`, `SpectraOptions` and `IO_Options`

&emsp;**NOTE:**

* Now, **channel4th** can only be used for wall-bounded channel flow (channel or open-channel), and only the partial semi-implicit scheme is available now. Fourth-order scheme for other Bcs and other viscous treatment will be considered in the future.
* The parameter sequence in a **namelist** is not important.
* A **namelist** is started with `&NamelistName` (e.g. `&BasicParam`), and ended with `/End`.
* You can also add blank lines or comments (after sign `!`) in the input file, just like what you do in a free format Fortran code.

&emsp;Consider the following input file:

```
!===================
&BasicParam
!===================

  RestartFlag = F         ! restart or not

  ! Flow type ( 1=Channel, 2=Half channel )
  FlowType  = 1
  ubulk= 0.666666666666666666666666666667
  IsUxConst = T
  IsUseCRF  = T          ! use Converting Reference Frame or not

  ! Mesh options
  xlx = 12.56637061      ! domain length in x-dir
  yly = 2.0              ! domain length in y-dir 
  zlz = 4.188790205      ! domain length in z-dir
  nxp =  385             ! grid point number in x-dir
  nyp =  193             ! grid point number in y-dir
  nzp =  181             ! grid point number in z-dir
  istret = 1             ! Stretch y-mesh or not (0:no, 1:both sides, 2:bottom)
  cStret = 1.0           ! Stretching parameter, if istret=0. this parameter doesn't work.

  ! Physical properties
  xnu = 2.3310E-4        ! kinematic viscosity
  gravity =0.0 0.0 0.0   ! Gravity or other constant body forces (if any)

  ! Time stepping
  dtMax= 0.02            ! Maxium time step
  iCFL = 2               ! Use CFL condition to change time step dynamically( 1: yes, 2:no ).
  CFLc = 1.5             ! CFL parameter
  ifirst= 1              ! First iteration
  ilast = 27000          ! Last iteration

  ! Numerical scheme options
  ischeme = 3                 ! (1=AB2, 2=RK2, 3=RK3)
  FFTW_plan_type = 1          ! (1=FFTW_MEASURE,  2=FFTW_ESTIMATE)

  ! Boundary conditions
  !                    x-,  x+,  y-,  y+,  z-,  z+
  ! And the velocity Bc values will ONLY be used corresponding to the no-slip Bc
  ! While at the same time, you are allow to assign a  transpiration Bc for uy in y- bottom(uyBcValue(3) dosen't work).
  uxBcValue = 0.0, 0.0,  0.0,  0.0,  0.0,  0.0
  uyBcValue = 0.0, 0.0,  0.0,  0.0,  0.0,  0.0
  uzBcValue = 0.0, 0.0,  0.0,  0.0,  0.0,  0.0 

  ! I/O, Statistics
  ivstats  = 10                  ! time step interval for statistics calculation
  saveStat = 3000                ! Output Statistics file frequency
  SaveVisu = 3000                ! Output visualizing file frequency
  BackupFreq = 9000              ! Output Restarting file frequency
  RunName  = "Cha180_"           ! Run name
  Res_Dir  = "./CFD/Results/"    ! Result directory
  RestartDir = "./CFD/Restart/"  ! Restart directory
  Cmd_LFile_Freq  = 5            ! Report frequency in the terminal
  LF_file_lvl     = 5            ! Logfile report level 
  LF_cmdw_lvl     = 3            ! Terminal report level

  ! Decomp2d options
  p_row = 4 !7
  p_col = 2 !7

  ! limited velocity and div
  vel_limit = 2.0
  div_limit = 0.2

/End of NAMELIST "&BasicParam"

!=================
&SpectraOptions
!=================

  clc_Spectra =  T
  ivSpec   = 30
  jSpecSet = 5
  jSpecEnd = 95
  jSpecInc = 5

/End of NAMELIST "SpectraOptions"

!=================
&IO_Options
!=================
  
  iskip   = 1
  jskip   = 1
  kskip   = 1
  XDMF_SET_TYPE=1  ! 0: cell, 1:Node

  save_ux    = T
  save_uy    = F
  save_uz    = F
  save_wx    = F
  save_wy    = F
  save_wz    = F
  save_wMag  = F
  save_pr    = F
  save_Q_vor = F
  save_lamda2= F

  WriteHistOld = T
  ReadHistOld  = T

/End of NAMELIST "&IO_Options"
```
&emsp;As you can see, most parameters are followed by some comments to illustrate its function.

## BasicParam
&emsp;**BasicParam** specifies several basic parameters, e.g. mesh options, physical properties.

* `RestartFlag`: logical type. Restart or not. If RestartFlag=.true., the simulation will start from a prestored restarting file.
* `FlowType`: integer type. Specify the flow type, 1=Channel, 2=Half channel.
* `ubulk`: real type. Mean bulk velocity in x-direction. This parameter is only used for wall-bounded turbulent channel flows (channel or half channel) to set the initial mean streamwise velocity.
* `IsUxConst`:  logical type. Whether ux is constant or not. This parameter is only used for wall-bounded turbulent channel flows (channel or half channel). If IsUxConst=.true., the turbulent channel flow is driven by a uniform pressure gradient, which varies in time to maintain a constant mass flow rate in streamwise (x) direction. If IsUxConst=.false., the uniform pressure gradient will be also constant in time, which is set by parameter `gravity`.
* `IsUseCRF`: logical type. Whether use Converting Reference Frame or not. If so, the computations are performed in a moving reference frame for which the bulk velocity (net streamwise mass flux) is zero.
* `xlx`: real type. Domain length in x-direction.
* `yly`: real type. Domain length in y-direction.
* `zlz`: real type. Domain length in z-direction.
* `nxp`: integer type. Grid point number in x-direction. **Note:** If nxp=129, there will be 129 grid points in x-direction, and the x-domain will be divided into 128 parts.
* `nyp`: integer type. Grid point number in y-direction.
* `nzp`: integer type. Grid point number in z-direction.
* `istret`: integer type. Whether y-mesh is stretched or not. 0:no,  1:both sides, 2:bottom. Sine function is used here. See file `../src/CFD_4th/m_MeshAndMetries.f90` for more information.
* `cStret`: real type. Corresponding stretching parameter, if istret=0. this parameter doesn't work.
* `xnu`: real type. Kinematic viscosity.
* `gravity`: real vector containing 3 components. Gravity or other constant body forces (if any), which can be used to drive the flow.
* `dtMax`: real type. Maxium allowable time step. If iCFL=1, `dt` will adjust according to the allowable CFL number. `dt=min{dtMax, CFLc/max(|u|/dx+|v|/dy+|w|/dz)}`. If iCFL=2, `dt=dtMax`. 
* `iCFL`: integer type. Whether use CFL condition to change time step dynamically or not. 1:yes, 2:no.
* `CFLc`: real type. Allowable CFL parameter. If iCFL=2, this parameter will not work.
* `ifirst`: integer type. First iteration.
* `ilast`:  integer type. Last iteration.
* `ischeme`: integer type. Specify the time integral scheme. 1=AB2, 2=RK2, 3=RK3.
* `FFTW_plan_type`: integer type. 1=FFTW_MEASURE, 2=FFTW_ESTIMATE. **Note:** In my practice, using *FFTW_MEASURE* is faster than *FFTW_ESTIMATE*, while *FFTW_MEASURE* might lead to different *FFTW_plan*, even for the same simulation case. And further, different *FFTW_plan* will result in slight different DFT values (The last several digits of the decimal point might be different). If you want to debug the code, please use *FFTW_ESTIMATE*, which will lead to the same DFT values for a specific case.
* `uxBcValue`: real vector containing 6 components to specify the u-velocity Bc values, and these values will ONLY be used corresponding to the no-slip Bc.
* `uyBcValue`: real vector containing 6 components to specify the v-velocity Bc values.
* `uzBcValue`: real vector containing 6 components to specify the w-velocity Bc values.
* `ivstats`: integer type. Time step interval for statistics calculation. The statistic value will be calculated every `ivstats` time step.
* `SaveStat`: integer type. Output Statistics file frequency. The statistic value will be written to the folder `Res_Dir` every `SaveStat` time step.
* `SaveVisu`: integer type. Output visualizing file frequency. The visualizing data will be written to the folder `Res_Dir` every `SaveVisu` time step.
* `BackupFreq`: integer type. Output restarting file frequency. The restarting file will be be written to the folder `RestartDir` every `BackupFreq` time step, and the previous restarting file will be deleted in order to save memory room.
* `RunName`: character type. Run name string.
* `Res_Dir`: character type. The folder to store result data.
* `RestartDir`: character type. The folder to store restart data.
* `Cmd_LFile_Freq`: integer type. Reporting frequency. Reporting information  will be printed to the terminal and written to the logfile every `Cmd_LFile_Freq` time step.
* `LF_file_lvl`: integer type. Logfile report level. From 1 to 5. `5` means every message will be reported into logfile, smaller `LF_file_lvl` corresponds to less reported message.
* `LF_cmdw_lvl`: integer type. Terminal report level. From 1 to 5.
* `p_row`: integer type. Number of processors in row pencil.
* `p_col`: integer type. Number of processors in column pencil.
* `vel_limit`: real type. Maximum allowable velocity. If `min{|u|,|v|,|w|}>vel_limit`, the program will abort.
* `div_limit`: real type. Maximum allowable divergence. If `max{abs(du/dx +dv/dy + dw/dz)}>div_limit`, the program will abort.

## SpectraOptions
&emsp;**SpectraOptions** designates options for energy spectra calculation.

* `clc_Spectra`: logical type. Calculate energy spectra or not.
* `ivSpec`: integer type. Time step interval for energy spectra calculation. The energy spectr will be calculated every `ivSpec` time step.
* `jSpecSet`: integer type. The first y-index for energy spectra calculation.
* `jSpecEnd`: integer type. The last y-index for energy spectra calculation.
* `jSpecInc`: integer type. The y-index interval for energy spectra calculation. Here I don't calculate energy spectra for every y-index (from 1 to nyp-1), and only energy spectra for [jSpecSet:jSpecEnd:jSpecInc] (from `jSpecSet` to `jSpecEnd`, every `jSpecInc` y-index) are calculated.

## IO_Options
&emsp;**IO_Options** sets input/output options.

* `iskip`: integer type. Flow data will be written into output visualizing file every `iskip` grid mesh in x-dir.
* `jskip`: integer type. Flow data will be written into output visualizing file every `jskip` grid mesh in y-dir.
* `kskip`: integer type. Flow data will be written into output visualizing file every `kskip` grid mesh in z-dir.
* `XDMF_SET_TYPE`: integer type. The output XDMF_SET_TYPE, 0: cell, 1:Node. 
* `save_ux`: logical type. Save u-velocity or not.
* `save_uy`: logical type. Save v-velocity or not.
* `save_uz`: logical type. Save w-velocity or not.
* `save_wx`: logical type. Save x-vorticity or not.
* `save_wy`: logical type. Save y-vorticity or not.
* `save_wz`: logical type. Save z-vorticity or not.
* `save_wMag`: logical type. Save vorticity magnitude or not.
* `save_pr`: logical type, Save pressure or not.
* `save_Q_vor`: logical type. Save Q vortex criterion value or not.
* `save_lamda2`: logical type. Save lamda2 vortex criterion value or not.
* `WriteHistOld`: logical type. Write the old convective data (more exactly, the terms which are treated explicitly) to restart value or not.
* `ReadHistOld`: logical type. Read the old convective data (more exactly, the terms which are treated explicitly) to restart value or not. This parameter, coupled with the previous one `WriteHistOld`, can help us to make the program restarting more flexible. Sometimes, we want to restart a `AB2` case from `RK3` case. `AB2` integral approach need the old convective data, while `RK3` method not. So, if we want to restart a `AB2` case from `RK3` case, we can set `ReadHistOld=F` for `AB2` case. If we want to restart a `AB2` case from a previous `AB2` case, we can set `ReadHistOld=T` for the new `AB2` case, and `WriteHistOld=T` for the previous case. If we want to restart a `RK3` case from `AB2` case, we can set `WriteHistOld=F` for the previous `AB2` case.

