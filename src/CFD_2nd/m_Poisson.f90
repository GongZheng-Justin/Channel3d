module m_Poisson
  use MPI
  use decomp_2d
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_MeshAndMetries
  use iso_c_binding
  use m_Variables,only: mb1
  use m_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal
  implicit none
  private
  include "fftw3.f03"
   
  ABSTRACT INTERFACE
    subroutine clcPPE_x(prsrc, prphiHalo)
      use m_TypeDef,  only: RK
      use m_Variables,only: mb1 
      use decomp_2d,  only: ystart,yend
      implicit none
      real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::prsrc
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
    end subroutine clcPPE_x
  END INTERFACE
  procedure(clcPPE_x),pointer::clcPPE
  
  ! engine-specific global variables
  real(RK),save::normfft
  type(C_PTR),save::fwd_plan_x,bwd_plan_x,fwd_plan_z,bwd_plan_z
  real(RK),save,allocatable,dimension(:)::ak1,ak3 ! modified waves in x-dir and z-dir
  
  public:: InitPoissonSolver, clcPPE
contains
    
  !******************************************************************
  ! InitPoissonSolver
  !******************************************************************     
  subroutine InitPoissonSolver()
    implicit none
    
    ! locals
    real(RK)::wa1,wa3,norm
    integer::i,k,iErr01, iErr02,iErrSum,plan_type
    type(fftw_iodim),dimension(1)::iodim,iodim_howmany
    integer(C_FFTW_R2R_KIND),dimension(1)::kind_fwd,kind_bwd
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1,arrx2
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1,arrz2

    if(BcOption(ym_dir)==BC_PERIOD ) then
      clcPPE => clcPPE_0
    else
      clcPPE => clcPPE_1
    endif
    
    if(FFTW_plan_type == 1) then
      plan_type=FFTW_MEASURE
    else
      plan_type=FFTW_ESTIMATE
    endif
    normfft = one

    ! fft in x
    iodim(1)%n  = xsize(1)
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = xsize(2)*xsize(3)
    iodim_howmany(1)%is = xsize(1)
    iodim_howmany(1)%os = xsize(1)
    IF(BcOption(xm_dir)==BC_PERIOD ) THEN
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
      norm = one
    else
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01 
      norm = two
    endif
    fwd_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrx1,arrx2,kind_fwd,plan_type)
    bwd_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrx1,arrx2,kind_bwd,plan_type)
    normfft = normfft*norm*real(nxc,kind=RK)
   
    ! fft in z
    iodim(1)%n  = zsize(3)
    iodim(1)%is = zsize(1)*zsize(2)
    iodim(1)%os = zsize(1)*zsize(2)
    iodim_howmany(1)%n  = zsize(1)*zsize(2)
    iodim_howmany(1)%is = 1
    iodim_howmany(1)%os = 1
    IF(BcOption(zm_dir)==BC_PERIOD) THEN
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
      norm = one
    else
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01 
      norm = two
    endif
    fwd_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrz1,arrz2,kind_fwd,plan_type)
    bwd_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrz1,arrz2,kind_bwd,plan_type)
    normfft = normfft*norm*real(nzc,kind=RK)
    normfft = one/normfft
       
    !   modified wave number in x1 and x3 direct.
    allocate(ak1(nxc), Stat =iErr01 )
    allocate(ak3(nzc), Stat =iErr02 )
    iErrSum=abs(iErr01)+abs(iErr02)
    if(iErrSum /= 0) call MainLog%CheckForError(ErrT_Abort,"InitPoissonSolver","Allocation failed")
    
    ! Here only two types of BCs are considered in x-dir and z-dir respectively.
    !  i.e.   two BC_PERIOD   or   two Neumann BCs
    IF(BcOption(xm_dir)==BC_PERIOD ) THEN
      do i=1,nxc     
        wa1= twopi*real(i-1,RK)/real(nxc,RK)
        ak1(i)=two*rdx2*(cos(wa1)-one)
      enddo
    ELSE
      do i=1,nxc        
        wa1= pi*real(i-1,RK)/real(nxc,RK)
        ak1(i)=two*rdx2*(cos(wa1)-one)
      enddo        
    ENDIF

    IF(BcOption(zm_dir)==BC_PERIOD ) THEN
      do k=1,nzc
        wa3= twopi*real(k-1,RK)/real(nzc,RK)
        ak3(k)=two*rdz2*(cos(wa3)-one)
      enddo
    ELSE
      do k=1,nzc
        wa3= pi*real(k-1,RK)/real(nzc,RK)
        ak3(k)=two*rdz2*(cos(wa3)-one)
      enddo
   ENDIF
    
  end subroutine InitPoissonSolver

  !******************************************************************
  ! clcPPE_0
  !******************************************************************
  subroutine clcPPE_0(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1,arrx2
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1,arrz2
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj

    call transpose_y_to_x(prsrc,arrx1)  
    call dfftw_execute_r2r(fwd_plan_x,arrx1,arrx2)
    call transpose_x_to_z(arrx2,arrz1)
    call dfftw_execute_r2r(fwd_plan_z,arrz1,arrz2)
    call transpose_z_to_y(arrz2,prsrc)

    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridmj(ic,jc)=am2ph(jc)
        tridpj(ic,jc)=ap2ph(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,yend(3)
        do jc=1,nyc
          do ic=1,yend(1)
            tridcj(ic,jc)=ac2ph(jc)+ak1(ic)+ak3(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prsrc(1:yend(1),1:nyc,kc),ysize(1),nyc)
      enddo

      ! For nrank=0, ystart(1)=ystart(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=zero
      tridcj(ic,jc)= one
      tridmj(ic,jc)=zero
      tridfj(ic,jc)=zero
      do ic=2,yend(1)
        tridcj(ic,jc)=ac2ph(jc)+ak1(ic)+ak3(kc)
        tridfj(ic,jc)=prsrc(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,yend(1)
          tridcj(ic,jc)=ac2ph(jc)+ak1(ic)+ak3(kc)
          tridfj(ic,jc)=prsrc(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=1,nyc
        do ic=1,yend(1)
          prsrc(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo

    ELSE
      do kc=ystart(3),yend(3)
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
            tridcj(ic,jc)=ac2ph(jc)+ak1(ic)+ak3(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prsrc(ystart(1):yend(1),1:nyc,kc),ysize(1),nyc)
      enddo       
    ENDIF
    call transpose_y_to_z(prsrc,arrz2)
    call dfftw_execute_r2r(bwd_plan_z,arrz2, arrz1)
    call transpose_z_to_x(arrz1,arrx2)
    call dfftw_execute_r2r(bwd_plan_x,arrx2, arrx1)
    call transpose_x_to_y(arrx1,prsrc)
    do kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_0
  
  !******************************************************************
  ! clcPPE_1
  !******************************************************************
  subroutine clcPPE_1(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1,arrx2
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1,arrz2
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj

    call transpose_y_to_x(prsrc,arrx1)  
    call dfftw_execute_r2r(fwd_plan_x,arrx1,arrx2)
    call transpose_x_to_z(arrx2,arrz1)
    call dfftw_execute_r2r(fwd_plan_z,arrz1,arrz2)
    call transpose_z_to_y(arrz2,prsrc)

    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridmj(ic,jc)=am2ph(jc)
        tridpj(ic,jc)=ap2ph(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,yend(3)
        do jc=1,nyc
          do ic=1,yend(1)
            tridcj(ic,jc)=ac2ph(jc)+ak1(ic)+ak3(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prsrc(1:yend(1),1:nyc,kc),ysize(1),nyc)
      enddo

      ! For nrank=0, ystart(1)=ystart(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=zero
      tridcj(ic,jc)=one
      tridmj(ic,jc)=zero
      tridfj(ic,jc)=zero
      do ic=2,yend(1)
        tridcj(ic,jc)=ac2ph(jc)+ak1(ic)+ak3(kc)
        tridfj(ic,jc)=prsrc(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,yend(1)
          tridcj(ic,jc)=ac2ph(jc)+ak1(ic)+ak3(kc)
          tridfj(ic,jc)=prsrc(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=1,nyc
        do ic=1,yend(1)
          prsrc(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo

    ELSE
      do kc=ystart(3),yend(3)
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
            tridcj(ic,jc)=ac2ph(jc)+ak1(ic)+ak3(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prsrc(ystart(1):yend(1),1:nyc,kc),ysize(1),nyc)
      enddo       
    ENDIF
    call transpose_y_to_z(prsrc,arrz2)
    call dfftw_execute_r2r(bwd_plan_z,arrz2, arrz1)
    call transpose_z_to_x(arrz1,arrx2)
    call dfftw_execute_r2r(bwd_plan_x,arrx2, arrx1)
    call transpose_x_to_y(arrx1,prsrc)
    do kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc) * normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_1

end module m_Poisson
