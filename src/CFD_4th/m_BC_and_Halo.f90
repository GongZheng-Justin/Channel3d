module m_BC_and_Halo
  use m_TypeDef
  use m_Parameters
  use decomp_2d
  use m_Variables,only: mb1,hi1,hi_pr,hi_uxPrSrc,hi_uzPrSrc
  implicit none
  private
  
  public:: SetBC_and_UpdateHalo,SetBC_and_UpdateHaloForPrSrc,SetBC_and_UpdateHalo_pr
contains
    
  !******************************************************************
  ! SetBC_and_UpdateHalo
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo(ux,uy,uz,uy_ym)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in)::uy_ym
    
    ! locals
    integer::ic,kc

    ! yp-dir
    SELECT CASE(FlowType)
    CASE(FT_CH)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic,nyp,  kc) = uxBcValue(yp_dir)*two -ux(ic, nyc, kc)
          uy(ic,nyp,  kc) = uyBcValue(yp_dir)
          uz(ic,nyp,  kc) = uzBcValue(yp_dir)*two -uz(ic, nyc, kc)
        enddo
      enddo 
    CASE(FT_HC)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic, nyp,  kc) = ux(ic, nyc, kc)
          uy(ic, nyp,  kc) = zero
          uz(ic, nyp,  kc) = uz(ic, nyc, kc)
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    do kc=ystart(3),yend(3)
      do ic=ystart(1),yend(1)
        ux(ic, 0, kc) = uxBcValue(ym_dir)*two-ux(ic, 1, kc)
        uy(ic, 1, kc) = uy_ym(ic,kc)
        uz(ic, 0, kc) = uzBcValue(ym_dir)*two-uz(ic, 1, kc)
      enddo
    enddo

    ! update halo
    call myupdate_halo(ux, mb1, hi1)
    call myupdate_halo(uy, mb1, hi1)
    call myupdate_halo(uz, mb1, hi1)
  
  end subroutine SetBC_and_UpdateHalo

  !******************************************************************
  ! SetBC_and_UpdateHaloForPrSrc
  !******************************************************************   
  subroutine SetBC_and_UpdateHaloForPrSrc(ux,uy,uz,uy_ym)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in)::uy_ym
    
    ! locals
    integer::ic,kc

    ! yp-dir
    SELECT CASE(FlowType)
    CASE(FT_CH)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic,nyp,  kc) = uyBcValue(yp_dir)
        enddo
      enddo 
    CASE(FT_HC)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic, nyp,  kc) = zero
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    do kc=ystart(3),yend(3)
      do ic=ystart(1),yend(1)
        uy(ic, 1, kc) = uy_ym(ic,kc)
      enddo
    enddo

    ! update halo
    call myupdate_halo(ux, mb1, hi_uxPrSrc)
    call myupdate_halo(uz, mb1, hi_uzPrSrc)
  
  end subroutine SetBC_and_UpdateHaloForPrSrc
  
  !******************************************************************
  ! SetBC_and_UpdateHalo_pr
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo_pr(pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure

    ! locals
    integer::ic,kc    
        
    do kc=ystart(3),yend(3)
      do ic=ystart(1),yend(1)
        pressure(ic, 0,   kc) = pressure(ic, 1,   kc)
        pressure(ic, nyp, kc) = pressure(ic, nyc, kc)
      enddo
    enddo

    call myupdate_halo(pressure, mb1, hi_pr)
  end subroutine SetBC_and_UpdateHalo_pr
  
end module m_BC_and_Halo
