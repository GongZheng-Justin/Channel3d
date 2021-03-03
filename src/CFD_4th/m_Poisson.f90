module m_Poisson
  use MPI
  use decomp_2d 
  use mydecomp_2d_extra
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_MeshAndMetries
  use iso_c_binding
  use m_Variables,only: mb1
  use m_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal
  implicit none
  private
  include "fftw3.f"
  
  type, bind(C) :: fftw_iodim
    integer(C_INT) n, is, os
  end type fftw_iodim
  interface
    type(C_PTR) function fftw_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,matin,matout,kindt,flags)  &
      bind(C, name='fftw_plan_guru_r2r')
      import
      integer(C_INT), value :: rank
      type(fftw_iodim), dimension(*), intent(in) :: dims
      integer(C_INT), value :: howmany_rank
      type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
      real(C_DOUBLE), dimension(*), intent(out) :: matin,matout
      integer(C_INT) :: kindt
      integer(C_INT), value :: flags
    end function fftw_plan_guru_r2r
  end interface
  
  ! engine-specific global variables
  type(C_PTR),save :: fwd_plan_x, bwd_plan_x
  type(C_PTR),save :: fwd_plan_z, bwd_plan_z
  real(RK),save::normfft
  
  real(RK),save,allocatable,dimension(:)::ak1 ! modified waves in x-dir
  real(RK),save,allocatable,dimension(:)::ak3 ! modified waves in z-dir
  
  public:: InitPoissonSolver, clcPPE
contains
    
  !******************************************************************
  ! InitPoissonSolver
  !******************************************************************     
  subroutine InitPoissonSolver()
    implicit none
    
    ! locals
    real(RK)::wa1,wa3,norm
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1,arrx2
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1,arrz2
    integer::i,k,iErr01,iErr02,iErrSum,kind_fwd,kind_bwd,plan_type
    type(fftw_iodim), dimension(1) :: iodim,iodim_howmany
    
    if(FFTW_plan_type == 1) then
      plan_type=FFTW_MEASURE
    else
      plan_type=FFTW_ESTIMATE
    endif

    ! fft in x
    iodim(1)%n  = xsize(1)
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = xsize(2)*xsize(3)
    iodim_howmany(1)%is = xsize(1)
    iodim_howmany(1)%os = xsize(1)
    kind_fwd = FFTW_R2HC
    kind_bwd = FFTW_HC2R
    fwd_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrx1,arrx2,kind_fwd,plan_type)
    bwd_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrx1,arrx2,kind_bwd,plan_type)
   
    ! fft in z
    iodim(1)%n  = zsize(3)
    iodim(1)%is = zsize(1)*zsize(2)
    iodim(1)%os = zsize(1)*zsize(2)
    iodim_howmany(1)%n  = zsize(1)*zsize(2)
    iodim_howmany(1)%is = 1
    iodim_howmany(1)%os = 1
    kind_fwd = FFTW_R2HC
    kind_bwd = FFTW_HC2R
    fwd_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrz1,arrz2,kind_fwd,plan_type)
    bwd_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,arrz1,arrz2,kind_bwd,plan_type)
    normfft = one/(real(nxc,RK)*real(nzc,RK))
       
    ! modified wave number in x1 and x3 direct.
    allocate(ak1(nxc), Stat =iErr01 )
    allocate(ak3(nzc), Stat =iErr02 )
    iErrSum=abs(iErr01)+abs(iErr02)
    if(iErrSum /= 0) call MainLog%CheckForError(ErrT_Abort,"InitPoissonSolver: ","Allocation failed")
    
    do i=1,nxc
      wa1= twopi*real(i-1,RK)/real(nxc,RK)
      ak1(i)=rdx2/288.0_RK*(cos(three*wa1)-54.0_RK*cos(two*wa1)+783.0_RK*cos(wa1)-730.0_RK)
    enddo
    do k=1,nzc
      wa3= twopi*real(k-1,RK)/real(nzc,RK)
      ak3(k)=rdz2/288.0_RK*(cos(three*wa3)-54.0_RK*cos(two*wa3)+783.0_RK*cos(wa3)-730.0_RK)
    enddo
  end subroutine InitPoissonSolver

  !******************************************************************
  ! clcPPE
  !******************************************************************
  subroutine clcPPE(prsrc, prphiHalo)
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

      ! for nrank=0, ystart(1)=ystart(3)=1
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
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE

end module m_Poisson
