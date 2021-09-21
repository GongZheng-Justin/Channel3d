module m_FlowType_AddedNew
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_MeshAndMetries
  use decomp_2d
  use m_Variables,only: mb1
  use m_Tools,only:CalcUxAver
  implicit none
  private    

  public:: InitVelocity_AN, Update_uy_ym_AN
  public:: InitStatVar_AN,  clcStat_AN
contains

  !******************************************************************
  ! InitVelocity_AN
  !******************************************************************
  subroutine InitVelocity_AN(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)),intent(inout)::Deviation
  

    ! locals
    integer :: ic,jc,kc
    real(RK):: VelRef,yct
      
    VelRef=one
    do kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        yct=yc(jc)
        do ic=ystart(1),yend(1)
          ux(ic,jc,kc)=  zero !VelRef*(two*yct/yly-one)
        enddo
      enddo
    enddo
    uy=zero
    uz=zero
  end subroutine InitVelocity_AN

  !******************************************************************
  ! Update_uy_ym_AN
  !******************************************************************   
  subroutine Update_uy_ym_AN(uy_ym, duy_ym, TimeNew)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(out):: uy_ym
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(inout):: duy_ym    
    real(RK),intent(in):: TimeNew
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = zero
     
    duy_ym = uy_ym - duy_ym
    
  end subroutine Update_uy_ym_AN

  !******************************************************************
  ! InitStatVar_AN
  !******************************************************************
  subroutine InitStatVar_AN()
    implicit none


  end subroutine InitStatVar_AN

  !******************************************************************
  ! clcStat_CH
  !******************************************************************
  subroutine clcStat_AN(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
   

  end subroutine clcStat_AN
    
end module m_FlowType_AddedNew
