module m_BC_and_Halo
  use m_TypeDef
  use m_Parameters
  use decomp_2d
  use m_Variables,only: mb1,hi1,hi_uxPrSrc,hi_uyPrSrc,hi_uzPrSrc
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
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic,nyp,  kc) = uxBcValue(yp_dir)*two -ux(ic, nyc, kc)
          uy(ic,nyp,  kc) = uyBcValue(yp_dir)
          uy(ic,nyp+1,kc) = uyBcValue(yp_dir)*two -uy(ic, nyc, kc)
          uz(ic,nyp,  kc) = uzBcValue(yp_dir)*two -uz(ic, nyc, kc)
        enddo
      enddo 
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic, nyp,  kc) = ux(ic, nyc, kc)
          uy(ic, nyp,  kc) = zero
          uy(ic, nyp+1,kc) =-uy(ic, nyc, kc)
          uz(ic, nyp,  kc) = uz(ic, nyc, kc)
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic, 0, kc) = uxBcValue(ym_dir)*two-ux(ic, 1, kc)
          uy(ic, 1, kc) = uy_ym(ic,kc)
          uy(ic, 0, kc) = uy_ym(ic,kc)*two -uy(ic,2,kc)
          uz(ic, 0, kc) = uzBcValue(ym_dir)*two-uz(ic, 1, kc)
        enddo
      enddo
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic, 0, kc) = ux(ic, 1, kc)
          uy(ic, 1, kc) = zero
          uy(ic, 0, kc) =-uy(ic, 2, kc)
          uz(ic, 0, kc) = uz(ic, 1, kc)
        enddo
      enddo    
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil, 3))
    CASE(BC_NSLIP)
      ux(nxp,  0:nyp,ystart(3):yend(3))=  uxBcValue(xp_dir)
      ux(nxp+1,0:nyp,ystart(3):yend(3))=  uxBcValue(xp_dir)*two -ux(nxc, 0:nyp,ystart(3):yend(3))
      uy(nxp,  0:nyp,ystart(3):yend(3))=  uyBcValue(xp_dir)*two -uy(nxc, 0:nyp,ystart(3):yend(3))
      uz(nxp,  0:nyp,ystart(3):yend(3))=  uzBcValue(xp_dir)*two -uz(nxc, 0:nyp,ystart(3):yend(3))
    CASE(BC_FSLIP)
      ux(nxp,  0:nyp,ystart(3):yend(3))= zero
      ux(nxp+1,0:nyp,ystart(3):yend(3))=-ux(nxc, 0:nyp,ystart(3):yend(3))
      uy(nxp,  0:nyp,ystart(3):yend(3))= uy(nxc, 0:nyp,ystart(3):yend(3))
      uz(nxp,  0:nyp,ystart(3):yend(3))= uz(nxc, 0:nyp,ystart(3):yend(3))        
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil, 4))
    CASE(BC_NSLIP)
      ux(1,0:nyp,ystart(3):yend(3)) = uxBcValue(xm_dir)
      ux(0,0:nyp,ystart(3):yend(3)) = uxBcValue(xm_dir)*two -ux(2,0:nyp,ystart(3):yend(3)) 
      uy(0,0:nyp,ystart(3):yend(3)) = uyBcValue(xm_dir)*two -uy(1,0:nyp,ystart(3):yend(3))
      uz(0,0:nyp,ystart(3):yend(3)) = uzBcValue(xm_dir)*two -uz(1,0:nyp,ystart(3):yend(3))       
    CASE(BC_FSLIP)
      ux(1,0:nyp,ystart(3):yend(3)) = zero
      ux(0,0:nyp,ystart(3):yend(3)) =-ux(2,0:nyp,ystart(3):yend(3)) 
      uy(0,0:nyp,ystart(3):yend(3)) = uy(1,0:nyp,ystart(3):yend(3))
      uz(0,0:nyp,ystart(3):yend(3)) = uz(1,0:nyp,ystart(3):yend(3))         
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(y_pencil, 1))
    CASE(BC_NSLIP) 
      ux(ystart(1):yend(1), 0:nyp, nzp  ) = uxBcValue(zp_dir)*two -ux(ystart(1):yend(1), 0:nyp, nzc) 
      uy(ystart(1):yend(1), 0:nyp, nzp  ) = uyBcValue(zp_dir)*two -uy(ystart(1):yend(1), 0:nyp, nzc)
      uz(ystart(1):yend(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
      uz(ystart(1):yend(1), 0:nyp, nzp+1) = uzBcValue(zp_dir)*two -uz(ystart(1):yend(1), 0:nyp, nzc)
    CASE(BC_FSLIP)
      ux(ystart(1):yend(1), 0:nyp, nzp  ) = ux(ystart(1):yend(1), 0:nyp, nzc) 
      uy(ystart(1):yend(1), 0:nyp, nzp  ) = uy(ystart(1):yend(1), 0:nyp, nzc)
      uz(ystart(1):yend(1), 0:nyp, nzp  ) = zero
      uz(ystart(1):yend(1), 0:nyp, nzp+1) =-uz(ystart(1):yend(1), 0:nyp, nzc)
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(y_pencil, 2))
    CASE(BC_NSLIP)
      ux(ystart(1):yend(1), 0:nyp, 0) = uxBcValue(zm_dir)*two -ux(ystart(1):yend(1), 0:nyp, 1)
      uy(ystart(1):yend(1), 0:nyp, 0) = uyBcValue(zm_dir)*two -uy(ystart(1):yend(1), 0:nyp, 1)
      uz(ystart(1):yend(1), 0:nyp, 1) = uzBcValue(zm_dir)
      uz(ystart(1):yend(1), 0:nyp, 0) = uzBcValue(zm_dir)*two -uz(ystart(1):yend(1), 0:nyp, 2)
    CASE(BC_FSLIP)
      ux(ystart(1):yend(1), 0:nyp, 0) = ux(ystart(1):yend(1), 0:nyp, 1)
      uy(ystart(1):yend(1), 0:nyp, 0) = uy(ystart(1):yend(1), 0:nyp, 1)
      uz(ystart(1):yend(1), 0:nyp, 1) = zero
      uz(ystart(1):yend(1), 0:nyp, 0) =-uz(ystart(1):yend(1), 0:nyp, 2)      
    END SELECT     

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
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic,nyp,  kc) = uyBcValue(yp_dir)
        enddo
      enddo 
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic, nyp,  kc) = zero
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic, 1, kc) = uy_ym(ic,kc)
        enddo
      enddo
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic, 1, kc) = zero
        enddo
      enddo    
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil, 3))
    CASE(BC_NSLIP)
      ux(nxp,  0:nyp,ystart(3):yend(3))=  uxBcValue(xp_dir)
    CASE(BC_FSLIP)
      ux(nxp,  0:nyp,ystart(3):yend(3))= zero      
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil, 4))
    CASE(BC_NSLIP)
      ux(1,0:nyp,ystart(3):yend(3)) = uxBcValue(xm_dir)      
    CASE(BC_FSLIP)
      ux(1,0:nyp,ystart(3):yend(3)) = zero         
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(y_pencil, 1))
    CASE(BC_NSLIP)
      uz(ystart(1):yend(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
    CASE(BC_FSLIP)
      uz(ystart(1):yend(1), 0:nyp, nzp  ) = zero
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(y_pencil, 2))
    CASE(BC_NSLIP)
      uz(ystart(1):yend(1), 0:nyp, 1) = uzBcValue(zm_dir)
    CASE(BC_FSLIP)
      uz(ystart(1):yend(1), 0:nyp, 1) = zero    
    END SELECT     

    ! update halo
    call myupdate_halo(ux, mb1, hi_uxPrSrc)
    call myupdate_halo(uy, mb1, hi_uyPrSrc)
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

    ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NSLIP, BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          pressure(ic, nyp, kc) = pressure(ic, nyc, kc)
        enddo
      enddo
    END SELECT         
        
    ! ym-dir
    SELECT CASE(BcOption(ym_dir))        
    CASE(BC_NSLIP, BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          pressure(ic, 0, kc) = pressure(ic, 1, kc)
        enddo
      enddo     
    END SELECT

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil, 3))
    CASE(BC_NSLIP, BC_FSLIP)
      pressure(nxp, 0:nyp,ystart(3):yend(3)) = pressure(nxc, 0:nyp,ystart(3):yend(3))        
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil, 4))
    CASE(BC_NSLIP, BC_FSLIP)
      pressure(0,0:nyp,ystart(3):yend(3)) = pressure(1,0:nyp,ystart(3):yend(3))       
    END SELECT      
        
    ! zp-dir
    SELECT CASE(myProcNghBC(y_pencil, 1))
    CASE(BC_NSLIP, BC_FSLIP)
      pressure(ystart(1):yend(1), 0:nyp, nzp) = pressure(ystart(1):yend(1), 0:nyp, nzc)         
    END SELECT         
    
    ! zm-dir
    SELECT CASE(myProcNghBC(y_pencil, 2))
    CASE(BC_NSLIP, BC_FSLIP)
      pressure(ystart(1):yend(1), 0:nyp, 0) = pressure(ystart(1):yend(1), 0:nyp, 1)         
    END SELECT     
  
    call myupdate_halo(pressure, mb1, hi1)
  end subroutine SetBC_and_UpdateHalo_pr
  
end module m_BC_and_Halo
