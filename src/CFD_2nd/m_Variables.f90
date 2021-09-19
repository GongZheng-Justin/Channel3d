module m_Variables
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use decomp_2d
  implicit none
  private
  
  ! define all major arrays here 
  real(RK), save,public,allocatable, dimension(:,:,:) :: ux
  real(RK), save,public,allocatable, dimension(:,:,:) :: uy
  real(RK), save,public,allocatable, dimension(:,:,:) :: uz
  real(RK), save,public,allocatable, dimension(:,:,:) :: HistxOld
  real(RK), save,public,allocatable, dimension(:,:,:) :: HistyOld
  real(RK), save,public,allocatable, dimension(:,:,:) :: HistzOld
  real(RK), save,public,allocatable, dimension(:,:,:) :: pressure
  
  real(RK), save,public,allocatable, dimension(:,:,:) :: RealArr1
  real(RK), save,public,allocatable, dimension(:,:,:) :: RealArr2
  real(RK), save,public,allocatable, dimension(:,:,:) :: Realhalo

  real(RK), save,public,allocatable, dimension(:,:,:) :: nut

  real(RK), save,public,allocatable, dimension(:,:) :: uy_ym, duy_ym

  type(MatBound),save,public:: mb1                                   ! matrix bound type 1
  type(HaloInfo),save,public:: hi1,hi_uxPrSrc,hi_uyPrSrc,hi_uzPrSrc  ! halo info type 1

#ifdef CFDDEM
  type(MatBound),save,public:: mb_dist  ! matrix bound for distribution
  type(HaloInfo),save,public:: hi_dist  ! halo info for distribution
  real(RK), save,public,allocatable, dimension(:,:,:)::FpForce_x,FpForce_y,FpForce_z
#endif

#ifdef CFDACM
  type(MatBound),save,public:: mb_dist  ! matrix bound for distribution
  type(HaloInfo),save,public:: hi_dist  ! halo info for distribution
  real(RK), save,public,allocatable, dimension(:,:,:):: FpForce_x,FpForce_y,FpForce_z
  real(RK), save,public,allocatable, dimension(:,:,:):: IBMArr1,IBMArr2,IBMArr3
  character,save,public,allocatable, dimension(:,:,:):: FluidIndicator
#endif
    
  public:: AllocateVariables
contains  

  !******************************************************************
  ! InverseTridiagonal
  !****************************************************************** 
  subroutine AllocateVariables()
    implicit none
      
    ! locals
    integer::iErr01, iErr02, iErr03, iErr04, iErr05, iErr06, iErr07, iErrSum

    mb1%pencil = y_pencil  
    mb1%xme=1;  mb1%xpe=2
    mb1%yme=1;  mb1%ype=2
    mb1%zme=1;  mb1%zpe=2

    hi1%pencil = y_pencil
    hi1%xmh=1;  hi1%xph=1
    hi1%zmh=1;  hi1%zph=1
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi1%ymh=1;  hi1%yph=1
    else
      hi1%ymh=0;  hi1%yph=0
    endif
        
    hi_uxPrSrc%pencil = y_pencil
    hi_uxPrSrc%xmh=0;  hi_uxPrSrc%xph=1
    hi_uxPrSrc%ymh=0;  hi_uxPrSrc%yph=0
    hi_uxPrSrc%zmh=0;  hi_uxPrSrc%zph=0

    hi_uyPrSrc%pencil = y_pencil
    hi_uyPrSrc%xmh=0;  hi_uyPrSrc%xph=0;
    hi_uyPrSrc%zmh=0;  hi_uyPrSrc%zph=0;
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_uyPrSrc%ymh=0;  hi_uyPrSrc%yph=1
    else
      hi_uyPrSrc%ymh=0;  hi_uyPrSrc%yph=0
    endif

    hi_uzPrSrc%pencil = y_pencil
    hi_uzPrSrc%xmh=0;  hi_uzPrSrc%xph=0
    hi_uzPrSrc%ymh=0;  hi_uzPrSrc%yph=0
    hi_uzPrSrc%zmh=0;  hi_uzPrSrc%zph=1

    !-------------------------------------------------
    ! Arrays with ghost cells
    !-------------------------------------------------
    call myallocate(ux, mb1, opt_global=.true.)
    call myallocate(uy, mb1, opt_global=.true.)
    call myallocate(uz, mb1, opt_global=.true.)
    call myallocate(pressure, mb1, opt_global=.true.)
    call myallocate(RealHalo, mb1, opt_global=.true.)
    if(LES_type>0) then 
      call myallocate(nut, mb1, opt_global=.true.)
      nut=zero
    endif
#ifdef CFDACM
    call myallocate(IBMArr1, mb1, opt_global=.true.)
    call myallocate(IBMArr2, mb1, opt_global=.true.)
    call myallocate(IBMArr3, mb1, opt_global=.true.)
#endif
   
    !-------------------------------------------------
    ! Arrays without ghost cells
    !-------------------------------------------------
    allocate(HistxOld(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), Stat = iErr01)
    allocate(HistyOld(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), Stat = iErr02)
    allocate(HistzOld(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), Stat = iErr03)
    allocate(RealArr1(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), Stat = iErr04)
    allocate(RealArr2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), Stat = iErr05)
#ifdef CFDACM
    allocate(FluidIndicator(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
#endif

    allocate(uy_ym(ystart(1):yend(1),  ystart(3):yend(3)), Stat = iErr06 )
    allocate(duy_ym(ystart(1):yend(1), ystart(3):yend(3)), Stat = iErr07 )
      
    iErrSum=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)+abs(iErr05)+abs(iErr06)+abs(iErr07)
    if(iErrSum /= 0) call MainLog%CheckForError(ErrT_Abort,"AllocateVariables","Allocation failed")     
    
    ux=zero;         uy=zero;        uz=zero
    HistxOld=zero;   HistyOld=zero;  HistzOld=zero
    pressure=zero;   RealArr1=zero;  RealArr2=zero;  RealHalo=zero;
    uy_ym=zero;      duy_ym=zero
#ifdef CFDACM
    IBMArr1=zero; IBMArr2=zero; IBMArr3=zero
#endif
      
  end subroutine AllocateVariables
    
end module m_Variables
