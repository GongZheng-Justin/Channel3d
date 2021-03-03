module m_FlowCase
  use m_TypeDef
  use m_Parameters
  use m_FlowType_Channel
  use m_FlowType_TGVortex
  use m_FlowType_AddedNew
  use decomp_2d
  use m_Variables,only: mb1
  implicit none
  private
  
  public:: InitVelocity, Update_uy_ym
  public:: InitStatVar,  clcStat
contains

  !******************************************************************
  ! SetBC_and_UpdateHalo
  !******************************************************************
  subroutine InitVelocity(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)),intent(inout)::Deviation
      
    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call InitVelocity_CH(ux,uy,uz,Deviation)
    case(FT_TG) ! Taylor-Green vortex
      call InitVelocity_TG(ux,uy,uz,Deviation)          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
      call InitVelocity_AN(ux,uy,uz,Deviation)          
    end select
    
  end subroutine InitVelocity

  !******************************************************************
  ! Update_uy_ym
  !******************************************************************   
  subroutine Update_uy_ym(uy_ym, duy_ym, TimeNew)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(inout):: uy_ym,duy_ym    
    real(RK),intent(in):: TimeNew
  
    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call Update_uy_ym_CH(uy_ym, duy_ym, TimeNew)
    case(FT_TG) ! Taylor-Green vortex
      call Update_uy_ym_TG(uy_ym, duy_ym, TimeNew)          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
      call Update_uy_ym_AN(uy_ym, duy_ym, TimeNew)          
    end select
    
  end subroutine Update_uy_ym

  !******************************************************************
  ! InitStatVar
  !******************************************************************
  subroutine InitStatVar()
    implicit none

    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call InitStatVar_CH()
    case(FT_TG) ! Taylor-Green vortex
      call InitStatVar_TG()          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
          
    end select

  end subroutine InitStatVar

  !******************************************************************
  ! clcStat
  !****************************************************************** 
  subroutine clcStat(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure

    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call clcStat_CH(ux,uy,uz,pressure)
    case(FT_TG) ! Taylor-Green vortex
      call clcStat_TG(ux,uy,uz,pressure)          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
          
    end select

  end subroutine clcStat 
    
end module m_FlowCase
