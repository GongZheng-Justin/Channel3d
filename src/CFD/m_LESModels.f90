module m_LESModels
  use decomp_2d
  use mydecomp_2d_extra
  use m_TypeDef
  use m_LogInfo
  use m_MeshAndMetries
  use m_Parameters
  use m_Variables,only: mb1,hi1
  implicit none
  private
     
  real(RK),dimension(:),allocatable,save::Fcoe

  public:: InitFilterCoe,Clc_SGS_Vis
contains

  !******************************************************************
  ! InitFilterCoe
  !******************************************************************
  subroutine InitFilterCoe()
    implicit none

    if(nrank==0) then
      !call MainLog%CheckForError(ErrT_Pass,"InitFilterCoe: ","NOTE: now the LES model in the code is ONLY used for channel flows(pbc in x-dir & z-dir, npbc in y-dir)!")
    endif
    allocate(Fcoe(9))
    if(FilterType==0) then     ! 0:trapezoidal type
      Fcoe = (/one, two,  one, two, four,    two,  one, two,  one/)/sixteen
    elseif(FilterType==1) then ! 1:Simpson type
      Fcoe = (/one, four, one, four,sixteen, four, one, four, one/)/36.00_RK
    else
      call MainLog%CheckForError(ErrT_Abort,"InitFilterCoe: ","unrecognized FilterType")
    endif
  end subroutine InitFilterCoe

  !******************************************************************
  ! Clc_SGS_Vis
  !******************************************************************
  subroutine Clc_SGS_Vis(ux,uy,uz,nut)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in) ::ux,uy,uz
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::nut

    SELECT CASE(LES_type)
    CASE(1)  ! 1: Smagorinsky Model

    CASE(2)  ! 2: constant Smagorinsky Model

    CASE(3)  ! 3: Dynamic Smagorinsky Model

    CASE(4)  ! 4: MTS Model
      call Clc_SGS_MTS(ux,uy,uz,nut)
    END SELECT 

    call SetBC_and_UpdateHalo_SGS(nut)
  end subroutine Clc_SGS_Vis

  !******************************************************************
  ! Calculate the SGS kinematic viscosity using Mixed-Time-Scale model 
  ! 
  ! Reference:
  ! Inagaki M, Kondoh T, Nagano Y, et al. 
  !   A Mixed-Time-Scale SGS Model With Fixed Model-Parameters for Practical LES. 
  !     Journal of Fluids Engineering-transactions of The Asme, 2005, 127(1): 1-13.
  !******************************************************************
  subroutine Clc_SGS_MTS(ux,uy,uz,nut)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in) ::ux,uy,uz
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::nut

    ! locals
    type(HaloInfo):: hi_K
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::uxf,uyf,uzf,Ksgs,Ts,rdelta
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))::strain_cha,Kes
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm)::Ktemp

    call ClcVelStrainRate(ux,uy,uz,strain_cha)

    ! filter ux, only in homogeneous directions
    DO kc=ystart(3),yend(3)
      km=kc-1
      kp=kc+1
      do jc=ystart(2),yend(2)    
        do ic=ystart(1),yend(1)
          im=ic-1
          ip=ic+1
          uxf= Fcoe(1)*ux(im,jc,km) + Fcoe(2)     *ux(ic,jc,km) +Fcoe(3)*ux(ip,jc,km) + &
               Fcoe(4)*ux(im,jc,kc) +(Fcoe(5)-one)*ux(ic,jc,kc) +Fcoe(6)*ux(ip,jc,kc) + &
               Fcoe(7)*ux(im,jc,kp) + Fcoe(8)     *ux(ic,jc,kp) +Fcoe(9)*ux(ip,jc,kp)
          Ktemp(ic,jc,kc)= uxf*uxf
        enddo
      enddo
    ENDDO
    hi_K%pencil = y_pencil
    hi_K%xmh=0;  hi_K%xph=1
    hi_K%ymh=0;  hi_K%yph=0
    hi_K%zmh=0;  hi_K%zph=0
    call myupdate_halo(Ktemp, mb1, hi_K)
    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          ip=ic+1
          Kes(ic,jc,kc)= half*(Ktemp(ic,jc,kc) +Ktemp(ip,jc,kc))
        enddo
      enddo
    ENDDO
    
    ! filter uy, only in homogeneous directions
    DO kc=ystart(3),yend(3)
      km=kc-1
      kp=kc+1
      do jc=1,nyp
        do ic=ystart(1),yend(1)
          im=ic-1
          ip=ic+1
          uyf= Fcoe(1)*uy(im,jc,km) + Fcoe(2)     *uy(ic,jc,km) +Fcoe(3)*uy(ip,jc,km) + &
               Fcoe(4)*uy(im,jc,kc) +(Fcoe(5)-one)*uy(ic,jc,kc) +Fcoe(6)*uy(ip,jc,kc) + &
               Fcoe(7)*uy(im,jc,kp) + Fcoe(8)     *uy(ic,jc,kp) +Fcoe(9)*uy(ip,jc,kp)
          Ktemp(ic,jc,kc)= uyf*uyf
        enddo
      enddo
    ENDDO
    DO kc=ystart(3),yend(3)
      do jc=1,nyc
        jp=jc+1
        do ic=ystart(1),yend(1)
          Kes(ic,jc,kc)= Kes(ic,jc,kc)+ half*(Ktemp(ic,jc,kc) +Ktemp(ic,jp,kc))
        enddo
      enddo
    ENDDO

    ! filter uz, only in homogeneous directions
    DO kc=ystart(3),yend(3)
      km=kc-1
      kp=kc+1
      do jc=ystart(2),yend(2)    
        do ic=ystart(1),yend(1)
          im=ic-1
          ip=ic+1
          uzf= Fcoe(1)*uz(im,jc,km) + Fcoe(2)     *uz(ic,jc,km) +Fcoe(3)*uz(ip,jc,km) + &
               Fcoe(4)*uz(im,jc,kc) +(Fcoe(5)-one)*uz(ic,jc,kc) +Fcoe(6)*uz(ip,jc,kc) + &
               Fcoe(7)*uz(im,jc,kp) + Fcoe(8)     *uz(ic,jc,kp) +Fcoe(9)*uz(ip,jc,kp)
          Ktemp(ic,jc,kc)= uzf*uzf
        enddo
      enddo
    ENDDO
    hi_K%pencil = y_pencil
    hi_K%xmh=0;  hi_K%xph=0
    hi_K%ymh=0;  hi_K%yph=0
    hi_K%zmh=0;  hi_K%zph=1
    call myupdate_halo(Ktemp, mb1, hi_K)
    DO kc=ystart(3),yend(3)
      kp=kc+1
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          Kes(ic,jc,kc)= Kes(ic,jc,kc)+ half*(Ktemp(ic,jc,kc) +Ktemp(ic,jc,kp))
        enddo
      enddo
    ENDDO

    ! Finally calculate the nut
    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        rdelta = one/DeltaCell(jc)
        do ic=ystart(1),yend(1)
          Ksgs= Kes(ic,jc,kc)
          Ts = sqrt(Ksgs)*rdelta + strain_cha(ic,jc,kc)/ten
          if(Ts < 1.00E-10_RK) then
            Ts= zero
          else
            Ts= one/Ts
          endif
          nut(ic,jc,kc)=0.05_RK*Ksgs*Ts + xnu
        enddo
      enddo
    ENDDO
     
  end subroutine Clc_SGS_MTS

  !******************************************************************
  ! ClcVelStrainTensor
  !******************************************************************
  subroutine ClcVelStrainTensor(ux,uy,uz,st1,st2,st3,st4,st5,st6)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::st1,st2,st3
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::st4,st5,st6

    ! locals
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::caj,cac1,cac2,cac12

    ! st1 = duxdx                st2 = duydy                st3 = duzdz
    ! st4 =(duxdy + duydx)*0.5   st5 =(duxdz + duzdx)*0.5   st6 =(duydz + duzdy)*0.5
    DO kc=ystart(3),yend(3)
      kp=kc+1
      km=kc-1
      do jc=ystart(2),yend(2)
        jp=jc+1
        jm=jc-1
        caj  = rdyp(jc)
        cac1 = rdyc(jc)
        cac2 = rdyc(jp)    
        cac12= cac1 - cac2
        do ic=ystart(1),yend(1)
          ip=ic+1
          im=ic-1
          st1(ic,jc,kc)=(ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
          st2(ic,jc,kc)=(uy(ic,jp,kc)-uy(ic,jc,kc))*caj
          st3(ic,jc,kc)=(uz(ic,jc,kp)-uz(ic,jc,kc))*rdz

          !rt1=(ux(ic,jc,kc)-ux(ic,jm,kc))*cac1 +(uy(ic,jc,kc)-uy(im,jc,kc))*rdx ! 2*st4 at (i,    j,    k+1/2)
          !rt2=(ux(ic,jp,kc)-ux(ic,jc,kc))*cac2 +(uy(ic,jp,kc)-uy(im,jp,kc))*rdx ! 2*st4 at (i,    j+1,  k+1/2)
          !rt3=(ux(ip,jc,kc)-ux(ip,jm,kc))*cac1 +(uy(ip,jc,kc)-uy(ic,jc,kc))*rdx ! 2*st4 at (i+1,  j,    k+1/2)
          !rt4=(ux(ip,jp,kc)-ux(ip,jc,kc))*cac2 +(uy(ip,jp,kc)-uy(ic,jp,kc))*rdx ! 2*st4 at (i+1,  j+1,  k+1/2)
          !st4(ic,jc,kc)=oneEighth*(rt1+rt2+rt3+rt4)                             !   st4 at (i+1/2,j+1/2,k+1/2)
          st4(ic,jc,kc)=oneEighth*( (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx  &
                                   +(ux(ip,jp,kc) +ux(ic,jp,kc))*cac2                             &
                                   +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12                            &
                                   -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1                             )

          !rt1=(ux(ic,jc,kc)-ux(ic,jc,km))*rdz  +(uz(ic,jc,kc)-uz(im,jc,kc))*rdx ! 2*st5 at (i,    j+1/2,k    )
          !rt2=(ux(ic,jc,kp)-ux(ic,jc,kc))*rdz  +(uz(ic,jc,kp)-uz(im,jc,kp))*rdx ! 2*st5 at (i,    j+1/2,k+1  )
          !rt3=(ux(ip,jc,kc)-ux(ip,jc,km))*rdz  +(uz(ip,jc,kc)-uz(ic,jc,kc))*rdx ! 2*st5 at (i+1,  j+1/2,k    )
          !rt4=(ux(ip,jc,kp)-ux(ip,jc,kc))*rdz  +(uz(ip,jc,kp)-uz(ic,jc,kp))*rdx ! 2*st5 at (i+1,  j+1/2,k+1  )
          !st5(ic,jc,kc)=oneEighth*(rt1+rt2+rt3+rt4)                             !   st5 at (i+1/2,j+1/2,k+1/2)
          st5(ic,jc,kc)=oneEighth*( (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz  &  
                                   +(uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx  )

          !rt1=(uz(ic,jc,kc)-uz(ic,jm,kc))*cac1 +(uy(ic,jc,kc)-uy(ic,jc,km))*rdz ! 2*st6 at (i+1/2,j,    k    )
          !rt2=(uz(ic,jp,kc)-uz(ic,jc,kc))*cac2 +(uy(ic,jp,kc)-uy(ic,jp,km))*rdz ! 2*st6 at (i+1/2,j+1,  k    )
          !rt3=(uz(ic,jc,kp)-uz(ic,jm,kp))*cac1 +(uy(ic,jc,kp)-uy(ic,jc,kc))*rdz ! 2*st6 at (i+1/2,j,    k+1  )
          !rt4=(uz(ic,jp,kp)-uz(ic,jc,kp))*cac2 +(uy(ic,jp,kp)-uy(ic,jp,kc))*rdz ! 2*st6 at (i+1/2,j+1,  k+1  )
          !st6(ic,jc,kc)=oneEighth*(rt1+rt2+rt3+rt4)                             !   st6 at (i+1/2,j+1/2,k+1/2)
          st6(ic,jc,kc)=oneEighth*( (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz  &
                                   +(uz(ic,jp,kp) +uz(ic,jp,kc))*cac2                             &
                                   +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12                            &
                                   -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1                             )
        enddo
      enddo
    ENDDO

  end subroutine ClcVelStrainTensor

  !******************************************************************
  ! ClcVelStrainRate
  !******************************************************************
  subroutine ClcVelStrainRate(ux,uy,uz,strain_cha)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::strain_cha

    ! locals
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::st1,st2,st3,st4,st5,st6,caj,cac1,cac2,cac12

    !  characteristic filtered rate of strain: strain_cha
    !  |S|= sqrt(2*S_ij*S_ij), where S_ij = 0.5*(dui/dxj + duj/dxi)
    DO kc=ystart(3),yend(3)
      kp=kc+1
      km=kc-1
      do jc=ystart(2),yend(2)
        jp=jc+1
        jm=jc-1
        caj  = rdyp(jc)
        cac1 = rdyc(jc)
        cac2 = rdyc(jp)    
        cac12= cac1 - cac2
        do ic=ystart(1),yend(1)
          ip=ic+1
          im=ic-1
          st1=(ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
          st2=(uy(ic,jp,kc)-uy(ic,jc,kc))*caj
          st3=(uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
          st4=oneEighth*( (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx  &
                         +(ux(ip,jp,kc) +ux(ic,jp,kc))*cac2                             &
                         +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12                            &
                         -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1                             )
          st5=oneEighth*( (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz  &  
                         +(uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx  )
          st6=oneEighth*( (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz  &
                         +(uz(ic,jp,kp) +uz(ic,jp,kc))*cac2                             &
                         +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12                            &
                         -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1                             )
          strain_cha(ic,jc,kc)= sqrt(two*(st1*st1+ st2*st2+ st3*st3)+ four*(st4*st4 + st5*st5 + st6*st6))
        enddo
      enddo
    ENDDO

  end subroutine ClcVelStrainRate

  !******************************************************************
  ! ClcVelStrainRate
  !******************************************************************
  subroutine SetBC_and_UpdateHalo_SGS(nut)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::nut

   ! locals
    integer::ic,kc

    ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NSLIP)           
      nut(ystart(1):yend(1), nyp, ystart(3):yend(3)) =  two*xnu-nut(ystart(1):yend(1), nyc, ystart(3):yend(3)) 
    CASE(BC_FSLIP) 
      nut(ystart(1):yend(1), nyp, ystart(3):yend(3)) =          nut(ystart(1):yend(1), nyc, ystart(3):yend(3))  
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NSLIP)        
      nut(ystart(1):yend(1), 0, ystart(3):yend(3)) =  two*xnu-nut(ystart(1):yend(1), 1, ystart(3):yend(3)) 
    CASE(BC_FSLIP)
      nut(ystart(1):yend(1), 0, ystart(3):yend(3)) =          nut(ystart(1):yend(1), 1, ystart(3):yend(3))   
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil, 3))
    CASE(BC_NSLIP)           
      nut(nxp, 0:nyp,ystart(3):yend(3))= two*xnu-nut(nxc, 0:nyp,ystart(3):yend(3))
    CASE(BC_FSLIP) 
      nut(nxp, 0:nyp,ystart(3):yend(3))=         nut(nxc, 0:nyp,ystart(3):yend(3))        
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil, 4))
    CASE(BC_NSLIP)        
      nut(0,0:nyp,ystart(3):yend(3)) = two*xnu-nut(1,0:nyp,ystart(3):yend(3))       
    CASE(BC_FSLIP)
      nut(0,0:nyp,ystart(3):yend(3)) =         nut(1,0:nyp,ystart(3):yend(3))         
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(y_pencil, 1))
    CASE(BC_NSLIP)         
      nut(ystart(1):yend(1), 0:nyp, nzp) = two*xnu-nut(ystart(1):yend(1), 0:nyp, nzc)
    CASE(BC_FSLIP)
      nut(ystart(1):yend(1), 0:nyp, nzp) =         nut(ystart(1):yend(1), 0:nyp, nzc) 
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(y_pencil, 2))
    CASE(BC_NSLIP)         
      nut(ystart(1):yend(1), 0:nyp, 0) = two*xnu-nut(ystart(1):yend(1), 0:nyp, 1)
    CASE(BC_FSLIP) 
      nut(ystart(1):yend(1), 0:nyp, 0) =         nut(ystart(1):yend(1), 0:nyp, 1)
    END SELECT     

    call myupdate_halo(nut, mb1, hi1)
  end subroutine SetBC_and_UpdateHalo_SGS

end module m_LESModels
