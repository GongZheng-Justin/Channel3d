module m_TScheme
  use MPI
  use decomp_2d
  use m_LogInfo
  use m_TypeDef  
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1, nut
#ifdef CFDDEM
  use m_Variables,only: FpForce_x,FpForce_y,FpForce_z
#endif
  use m_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal,InversePTriFixedCoe
  implicit none
  private 

  ABSTRACT INTERFACE
    subroutine clcRhsX_xxx(ux,uy,uz,RhsX,HistXold,pressure)
      use m_TypeDef,  only: RK
      use m_Variables,only: mb1 
      use decomp_2d,  only: ystart,yend    
      implicit none
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
      real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::RhsX
      real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistXold
    end subroutine clcRhsX_xxx

    subroutine  clcRhsY_xxx(ux,uy,uz,RhsY,HistYold,pressure)
      use m_TypeDef,  only: RK
      use m_Variables,only: mb1 
      use decomp_2d,  only: ystart,yend
      implicit none
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
      real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::RhsY
      real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistYold
    end subroutine  clcRhsY_xxx

    subroutine  clcRhsZ_xxx(ux,uy,uz,RhsZ,HistZold,pressure)
      use m_TypeDef,  only: RK
      use m_Variables,only: mb1 
      use decomp_2d,  only: ystart,yend
      implicit none
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
      real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistZold
    end subroutine  clcRhsZ_xxx

    subroutine clcU1Hat_xxx(ux,RhsX)
      use m_TypeDef,  only: RK
      use m_Variables,only: mb1 
      use decomp_2d,  only: ystart,yend
      implicit none
      real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsX
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux
    end subroutine clcU1Hat_xxx  

    subroutine clcU2Hat_xxx(uy,RhsY,duy_ym)
      use m_TypeDef,  only: RK
      use m_Variables,only: mb1 
      use decomp_2d,  only: ystart,yend
      implicit none
      real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsY
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
      real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in):: duy_ym
    end subroutine clcU2Hat_xxx  

    subroutine clcU3Hat_xxx(uz,RhsZ)
      use m_TypeDef,  only: RK
      use m_Variables,only: mb1 
      use decomp_2d,  only: ystart,yend
      implicit none
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
    end subroutine clcU3Hat_xxx

    subroutine clcPrSrc_xxx(ux,uy,uz,prsrc,pressure,divmax)
      use m_TypeDef,  only: RK
      use m_Variables,only: mb1 
      use decomp_2d,  only: ystart,yend
      implicit none
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
      real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::prsrc
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
      real(RK),intent(out)::divmax
    end subroutine clcPrSrc_xxx

    subroutine PressureUpdate_xxx(pressure, prphiHalo)
      use m_TypeDef,  only: RK
      use m_Variables,only: mb1 
      implicit none
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
      real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out):: pressure
    end subroutine PressureUpdate_xxx
  END INTERFACE
  procedure(clcRhsX_xxx), pointer,public::clcRhsX
  procedure(clcRhsY_xxx), pointer,public::clcRhsY
  procedure(clcRhsZ_xxx), pointer,public::clcRhsZ
  procedure(clcU1Hat_xxx),pointer,public::clcU1Hat
  procedure(clcU2Hat_xxx),pointer,public::clcU2Hat
  procedure(clcU3Hat_xxx),pointer,public::clcU3Hat
  procedure(clcPrSrc_xxx),pointer,public::clcPrSrc
  procedure(PressureUpdate_xxx),pointer,public:: PressureUpdate
  
  public:: InitTimeScheme,FluidVelUpdate
contains    

#include "m_TSchemeFEXP.inc"
#include "m_TSchemePIMP.inc"
#include "m_TSchemeFIMP.inc"

  !******************************************************************
  ! InitTimeScheme
  !******************************************************************
  subroutine InitTimeScheme()
    implicit none

    ! FEXP 0: full explicit
    ! PIMP 1: partial implicit, only use C-N in y-dir 
    ! FIMP 2: full implicit, use C-N in all 3 dirs.
    if( (BcOption(xm_dir)==BC_PERIOD .and. BcOption(xp_dir)/=BC_PERIOD) .or.  (BcOption(ym_dir)==BC_PERIOD .and. BcOption(yp_dir)/=BC_PERIOD) .or. &
        (BcOption(zm_dir)==BC_PERIOD .and. BcOption(zp_dir)/=BC_PERIOD)   ) then
      call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Periodic Bc Wrong ")
    endif

    ! Note here, if use full implicit time scheme, the potential periodic bc should be in x-dir first,a nd then z-dir, then y-dir
    if(IsImplicit==2) then
      if(BcOption(xm_dir) /=BC_PERIOD .and. (BcOption(ym_dir)==BC_PERIOD .or. BcOption(zm_dir)==BC_PERIOD) ) then
        call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Bc type wrong 1 ")
      endif
      if(BcOption(zm_dir) /=BC_PERIOD .and. BcOption(ym_dir) ==BC_PERIOD ) then
        call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Bc type wrong 2 ")
      endif
    endif

    if(LES_type==0) then
      SELECT CASE(IsImplicit)
      CASE(0)
        clcRhsX => clcRhsX_FEXP
        clcRhsY => clcRhsY_FEXP
        clcRhsZ => clcRhsZ_FEXP
        clcU1Hat=> clcU1Hat_FEXP
        clcU2Hat=> clcU2Hat_FEXP
        clcU3Hat=> clcU3Hat_FEXP
        clcPrSrc=> clcPrSrcOther
        PressureUpdate => PressureUpdate_FEXP
      CASE(1)
        clcRhsX => clcRhsX_PIMP
        clcRhsY => clcRhsY_PIMP
        clcRhsZ => clcRhsZ_PIMP
        if(BcOption(ym_dir)==BC_PERIOD) then
          clcU1Hat=> clcU1Hat_PIMP_0
          clcU2Hat=> clcU2Hat_PIMP_0
          clcU3Hat=> clcU3Hat_PIMP_0
        else
          clcU1Hat=> clcU1Hat_PIMP
          clcU2Hat=> clcU2Hat_PIMP
          clcU3Hat=> clcU3Hat_PIMP
        endif
        clcPrSrc=> clcPrSrcOther
        PressureUpdate => PressureUpdate_PIMP
      CASE(2)
        clcRhsX => clcRhsX_FIMP
        clcRhsY => clcRhsY_FIMP
        clcRhsZ => clcRhsZ_FIMP
        if(BcOption(xm_dir) ==BC_PERIOD) then ! There can be several periodic Bcs
          if(BcOption(zm_dir) ==BC_PERIOD) then
            if(BcOption(ym_dir) ==BC_PERIOD) then !
              clcU1Hat=> clcU1Hat_FIMP_000
              clcU2Hat=> clcU2Hat_FIMP_000 
              clcU3Hat=> clcU3Hat_FIMP_000    
            else
              clcU1Hat=> clcU1Hat_FIMP_010
              clcU2Hat=> clcU2Hat_FIMP_010
              clcU3Hat=> clcU3Hat_FIMP_010
            endif
          else
            clcU1Hat=> clcU1Hat_FIMP_011
            clcU2Hat=> clcU2Hat_FIMP_011
            clcU3Hat=> clcU3Hat_FIMP_011
          endif
        else                                  ! NO periodic Bc exist
          clcU1Hat=> clcU1Hat_FIMP_111
          clcU2Hat=> clcU2Hat_FIMP_111
          clcU3Hat=> clcU3Hat_FIMP_111
        endif
        clcPrSrc=> clcPrSrc_FIMP
        PressureUpdate => PressureUpdate_FIMP
      END SELECT
    else
      SELECT CASE(IsImplicit)
      CASE(0)
        clcRhsX => clcRhsX_FEXP_LES
        clcRhsY => clcRhsY_FEXP_LES
        clcRhsZ => clcRhsZ_FEXP_LES
        clcU1Hat=> clcU1Hat_FEXP_LES
        clcU2Hat=> clcU2Hat_FEXP_LES
        clcU3Hat=> clcU3Hat_FEXP_LES
        clcPrSrc=> clcPrSrcOther
        PressureUpdate => PressureUpdate_FEXP_LES
      CASE(1)
        clcRhsX => clcRhsX_PIMP_LES
        clcRhsY => clcRhsY_PIMP_LES
        clcRhsZ => clcRhsZ_PIMP_LES
        if(BcOption(ym_dir)==BC_PERIOD) then
          clcU1Hat=> clcU1Hat_PIMP_LES_0
          clcU2Hat=> clcU2Hat_PIMP_LES_0
          clcU3Hat=> clcU3Hat_PIMP_LES_0
        else
          clcU1Hat=> clcU1Hat_PIMP_LES
          clcU2Hat=> clcU2Hat_PIMP_LES
          clcU3Hat=> clcU3Hat_PIMP_LES
        endif
        clcPrSrc=> clcPrSrcOther
        PressureUpdate => PressureUpdate_PIMP_LES
      CASE(2)
        clcRhsX => clcRhsX_FIMP_LES
        clcRhsY => clcRhsY_FIMP_LES
        clcRhsZ => clcRhsZ_FIMP_LES
        if(BcOption(xm_dir) ==BC_PERIOD) then ! There can be several periodic Bcs
          if(BcOption(zm_dir) ==BC_PERIOD) then
            if(BcOption(ym_dir) ==BC_PERIOD) then !
              clcU1Hat=> clcU1Hat_FIMP_LES_000
              clcU2Hat=> clcU2Hat_FIMP_LES_000 
              clcU3Hat=> clcU3Hat_FIMP_LES_000
              ! a nut metrix in x-pencil and in z-pencil is needed, and update_halo in x-pencil and z-pencil is also needed!!!
              ! after fineshed those, the LES for full semi-implicit C-N ca be used!
              call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Have not finished yet - 1 !!!")    
            else
              clcU1Hat=> clcU1Hat_FIMP_LES_010
              clcU2Hat=> clcU2Hat_FIMP_LES_010
              clcU3Hat=> clcU3Hat_FIMP_LES_010
              call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Have not finished yet - 2 !!!")
            endif
          else
            !clcU1Hat=> clcU1Hat_FIMP_LES_011
            !clcU2Hat=> clcU2Hat_FIMP_LES_011
            !clcU3Hat=> clcU3Hat_FIMP_LES_011
            call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Have not finished yet - 3 !!!")
          endif
        else                                  ! NO periodic Bc exist
          !clcU1Hat=> clcU1Hat_FIMP_LES_111
          !clcU2Hat=> clcU2Hat_FIMP_LES_111
          !clcU3Hat=> clcU3Hat_FIMP_LES_111
          call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Have not finished yet - 4 !!!")
        endif
        clcPrSrc=> clcPrSrcOther
        PressureUpdate => PressureUpdate_FIMP_LES
      END SELECT 
    endif
  end subroutine InitTimeScheme

  !******************************************************************
  ! clcPrSrcOther
  !******************************************************************  
  subroutine clcPrSrcOther(ux,uy,uz,prsrc,pressure,divmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
    real(RK),intent(out)::divmax
    
    ! locals
    integer::ic,jc,kc,ip,jp,kp,ierror
    real(RK)::sudtal,sucaj,rdiv,divmax1

    divmax1=zero
    sudtal=one/pmAlpha
    DO kc=ystart(3),yend(3)
       kp=kc+1
       do jc=ystart(2),yend(2)
         jp=jc+1
         sucaj=rdyp(jc)
         do ic=ystart(1),yend(1)
           ip=ic+1
           rdiv= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx + (uy(ic,jp,kc)-uy(ic,jc,kc))*sucaj + &
                 (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
           divmax1=max(abs(rdiv),divmax1)
           prsrc(ic,jc,kc)= sudtal * rdiv
         enddo
       enddo
     ENDDO
     call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)
  end subroutine clcPrSrcOther    

  !******************************************************************
  ! clcPrSrc_FIMP
  !******************************************************************  
  subroutine clcPrSrc_FIMP(ux,uy,uz,prsrc,pressure,divmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
    real(RK),intent(out)::divmax
    
    ! locals
    integer::ic,jc,kc,ip,jp,kp,ierror
    real(RK)::sudtal,sucaj,rdiv,divmax1,xnuhm

    divmax1=zero
    xnuhm= -half*xnu
    sudtal=one/pmAlpha
    DO kc=ystart(3),yend(3)
       kp=kc+1
       do jc=ystart(2),yend(2)
         jp=jc+1
         sucaj=rdyp(jc)
         do ic=ystart(1),yend(1)
           ip=ic+1
           rdiv= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx + (uy(ic,jp,kc)-uy(ic,jc,kc))*sucaj + &
                 (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
           divmax1=max(abs(rdiv),divmax1)
           prsrc(ic,jc,kc)= sudtal * rdiv
           pressure(ic,jc,kc)= pressure(ic,jc,kc)+ xnuhm*rdiv         ! new added for full implicit scheme
         enddo
       enddo
     ENDDO
     call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)
  end subroutine clcPrSrc_FIMP
  
  !******************************************************************
  ! FluidVelUpdate
  !******************************************************************
  subroutine FluidVelUpdate(prphiHalo,ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out):: ux,uy,uz
    
    ! locals
    integer::ic,jc,kc,im,jm,km
    real(RK)::rdxEta2,sucacEta2,rdzEta2,locphi

    rdxEta2=rdx*pmAlpha
    rdzEta2=rdz*pmAlpha
    DO kc=ystart(3),yend(3)
      km=kc-1
      do jc=ystart(2),yend(2)
        jm=jc-1
        sucacEta2=rdyc(jc)*pmAlpha
        do ic=ystart(1),yend(1)
          im=ic-1
          locphi=prphiHalo(ic,jc,kc)
          ux(ic,jc,kc)=ux(ic,jc,kc)-  (locphi-prphiHalo(im,jc,kc))*rdxEta2
          uy(ic,jc,kc)=uy(ic,jc,kc)-  (locphi-prphiHalo(ic,jm,kc))*sucacEta2
          uz(ic,jc,kc)=uz(ic,jc,kc)-  (locphi-prphiHalo(ic,jc,km))*rdzEta2                
        enddo
      enddo
    ENDDO
  end subroutine FluidVelUpdate

end module m_TScheme
