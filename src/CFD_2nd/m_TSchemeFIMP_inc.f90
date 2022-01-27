
  ! This file is included in the module m_TScheme

  !******************************************************************
  ! clcRhsX_FIMP
  !******************************************************************    
  subroutine  clcRhsX_FIMP(ux,uy,uz,RhsX,HistXold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::RhsX
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistXold
    
    ! locals    
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,ierror
    real(RK)::qdx1,qdx3,h11,h12,h13,sucaj,s3tot,s3tot1,dp1ns,dpmdxns
    real(RK)::d11q1,d22q1,d33q1,dcq123,convEd1,gradp1,InterpY1,InterpY2,InterpY3,InterpY4
    
    s3tot=zero
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=ystart(3),yend(3)
      km=kc-1; kp=kc+1
      do jc=ystart(2),yend(2)
        jm=jc-1; jp=jc+1
        sucaj=half*rdyp(jc)
        InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=one-InterpY3            
        do ic=ystart(1),yend(1)
          im=ic-1; ip=ic+1

          h11=( (ux(ip,jc,kc)+ux(ic,jc,kc))* (ux(ip,jc,kc)+ux(ic,jc,kc))  &
               -(ux(ic,jc,kc)+ux(im,jc,kc))* (ux(ic,jc,kc)+ux(im,jc,kc)) )*qdx1 
          h12=( (uy(ic,jp,kc)+uy(im,jp,kc))* (InterpY3*ux(ic,jc,kc) +InterpY4*ux(ic,jp,kc))  &
               -(uy(ic,jc,kc)+uy(im,jc,kc))* (InterpY1*ux(ic,jm,kc) +InterpY2*ux(ic,jc,kc)) )*sucaj            
          h13=( (uz(ic,jc,kp)+uz(im,jc,kp))* (ux(ic,jc,kp)+ux(ic,jc,kc))  &
               -(uz(ic,jc,kc)+uz(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jc,km)) )*qdx3
                
          d11q1= (ux(ip,jc,kc)-two*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          d22q1= ap2c(jc)*ux(ic,jp,kc) + ac2c(jc)*ux(ic,jc,kc)+ am2c(jc)*ux(ic,jm,kc)                
          d33q1= ap3c(kc)*ux(ic,jc,kp) + ac3c(kc)*ux(ic,jc,kc)+ am3c(kc)*ux(ic,jc,km)
          dcq123= d11q1+d22q1+d33q1
             
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
#if defined CFDDEM || defined CFDACM || defined CFDLPT_TwoWay
          s3tot=s3tot+ux(ic,jc,kc)*dyp(jc)
#else
          s3tot=s3tot+dcq123*dyp(jc)
#endif
#ifdef CFDDEM
          convEd1= -h11-h12-h13 + gravity(1)+half*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
#elif CFDLPT_TwoWay
#ifdef ForceInCell
          convEd1= -h11-h12-h13 + gravity(1)+half*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
#else
          convEd1= -h11-h12-h13 + gravity(1)+FpForce_x(ic,jc,kc)
#endif
#else
          convEd1= -h11-h12-h13 + gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1+ two*pmBeta*dcq123
          HistXold(ic,jc,kc)=convEd1
        enddo
      enddo
    ENDDO
  
    ! in dp1ns there is the mean pressure gradient to keep constant mass
#if defined CFDDEM || defined CFDACM || defined CFDLPT_TwoWay
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(s3tot,s3tot1,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      s3tot1=s3tot1/(real(nxc*nzc,kind=RK))/yly
      dpmdxns= ubulk - s3tot1

      dp1ns  = pmAlphaC* dpmdxns
      do kc=ystart(3),yend(3)
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+ dp1ns
          enddo
        enddo
      enddo
      PrGradAver = PrGradAver+ dpmdxns * pmAlphaC/dt
    ENDIF
#else
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(s3tot,s3tot1,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      dpmdxns= xnu*s3tot1/(real(nxc*nzc,kind=RK))/yly
      dp1ns  = pmAlpha* dpmdxns
      do kc=ystart(3),yend(3)
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc)- dp1ns
          enddo
        enddo
      enddo
      PrGradAver = PrGradAver+ dpmdxns * pmAlphaC
    ENDIF
#endif

  end subroutine clcRhsX_FIMP

  !******************************************************************
  ! clcRhsX_FIMP_LES
  !******************************************************************    
  subroutine  clcRhsX_FIMP_LES(ux,uy,uz,RhsX,HistXold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::RhsX
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistXold
    
    ! locals    
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,ierror
    real(RK)::qdx1,qdx3,h11,h12,h13,sucaj,s3tot,s3tot1,dp1ns,dpmdxns
    real(RK)::d11q1,d22q1,d33q1,dcq123,convEd1,gradp1
    real(RK)::visa,visb,visc,visd,sgs1,sgs2,sgs3
    
    s3tot=zero
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=ystart(3),yend(3)
      km=kc-1
      kp=kc+1
      do jc=ystart(2),yend(2)
        jm=jc-1
        jp=jc+1
        sucaj=quarter*rdyp(jc)            
        do ic=ystart(1),yend(1)
          im=ic-1
          ip=ic+1
          h11=( (ux(ip,jc,kc)+ux(ic,jc,kc))* (ux(ip,jc,kc)+ux(ic,jc,kc))  &
               -(ux(ic,jc,kc)+ux(im,jc,kc))* (ux(ic,jc,kc)+ux(im,jc,kc)) )*qdx1

          h12=( (uy(ic,jp,kc)+uy(im,jp,kc))* (ux(ic,jp,kc)+ux(ic,jc,kc))  &
               -(uy(ic,jc,kc)+uy(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jm,kc)) )*sucaj            
          h13=( (uz(ic,jc,kp)+uz(im,jc,kp))* (ux(ic,jc,kp)+ux(ic,jc,kc))  &
               -(uz(ic,jc,kc)+uz(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jc,km)) )*qdx3

          visa=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc))
          visb=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jp,kc)+nut(im,jp,kc))
          visc=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,km)+nut(im,jc,km))
          visd=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,kp)+nut(im,jc,kp))
          sgs1=rdx2*(nut(ic,jc,kc)*(ux(ip,jc,kc)-ux(ic,jc,kc))-nut(im,jc,kc)*(ux(ic,jc,kc)-ux(im,jc,kc)))
          sgs2=rdx*rdyp(jc)*(visb*(uy(ic,jp,kc)-uy(im,jp,kc))-visa*(uy(ic,jc,kc)-uy(im,jc,kc)))
          sgs3=rdx*rdz     *(visd*(uz(ic,jc,kp)-uz(im,jc,kp))-visc*(uz(ic,jc,kc)-uz(im,jc,kc)))
                
          d11q1= sgs1
          d22q1= ap2c(jc)*visb*(ux(ic,jp,kc)-ux(ic,jc,kc))-am2c(jc)*visa*(ux(ic,jc,kc)-ux(ic,jm,kc))
          d33q1= ap3c(kc)*visd*(ux(ic,jc,kp)-ux(ic,jc,kc))-am3c(kc)*visc*(ux(ic,jc,kc)-ux(ic,jc,km))
          dcq123= d11q1+d22q1+d33q1
                
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
#if defined CFDDEM || defined CFDACM || defined CFDLPT_TwoWay
          s3tot=s3tot+ux(ic,jc,kc)*dyp(jc)
#else
          s3tot=s3tot+dcq123*dyp(jc)
#endif              
#ifdef CFDDEM
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+gravity(1)+half*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
#elif CFDLPT_TwoWay
#ifdef ForceInCell
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+gravity(1)+half*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
#else
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+gravity(1)+FpForce_x(ic,jc,kc)
#endif
#else
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1+ two*pmBetaT*dcq123
          HistXold(ic,jc,kc)=convEd1
        enddo
      enddo
    ENDDO
  
    ! in dp1ns there is the mean pressure gradient to keep constant mass
#if defined CFDDEM || defined CFDACM || defined CFDLPT_TwoWay
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(s3tot,s3tot1,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      s3tot1=s3tot1/(real(nxc*nzc,kind=RK))/yly
      dpmdxns= ubulk - s3tot1

      dp1ns  = pmAlphaC* dpmdxns
      do kc=ystart(3),yend(3)
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+ dp1ns
          enddo
        enddo
      enddo
      PrGradAver = PrGradAver+ dpmdxns * pmAlphaC/dt
    ENDIF
#else
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(s3tot,s3tot1,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      dpmdxns= s3tot1/(real(nxc*nzc,kind=RK))/yly
      dp1ns  = pmAlpha* dpmdxns
      do kc=ystart(3),yend(3)
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc)- dp1ns
          enddo
        enddo
      enddo
      PrGradAver = PrGradAver+ dpmdxns * pmAlphaC
    ENDIF
#endif
  end subroutine clcRhsX_FIMP_LES

  !******************************************************************
  ! clcRhsY_FIMP
  !******************************************************************    
  subroutine  clcRhsY_FIMP(ux,uy,uz,RhsY,HistYold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure 

    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::RhsY
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistYold
   
    ! locals 
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::hdx1,hdx3,h21,h22,h23,sucac,qsucac
    real(RK)::d11q2,d22q2,d33q2,dcq123,convEd2,gradp2,InterpY1,InterpY2    
    
    hdx1=half*rdx
    hdx3=half*rdz
    DO kc=ystart(3),yend(3)
      km=kc-1; kp=kc+1
      do jc=ystart(2),yend(2)
        jm=jc-1; jp=jc+1
        sucac = rdyc(jc)
        qsucac= quarter*sucac
        InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1
        do ic=ystart(1),yend(1)
          im=ic-1; ip=ic+1

          h21=( (InterpY1*ux(ip,jm,kc)+InterpY2*ux(ip,jc,kc))* (uy(ip,jc,kc)+uy(ic,jc,kc)) &
               -(InterpY1*ux(ic,jm,kc)+InterpY2*ux(ic,jc,kc))* (uy(ic,jc,kc)+uy(im,jc,kc)) )*hdx1
          h22=( (uy(ic,jp,kc)+uy(ic,jc,kc))* (uy(ic,jp,kc)+uy(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jm,kc)) )*qsucac
          h23=( (InterpY1*uz(ic,jm,kp)+InterpY2*uz(ic,jc,kp))* (uy(ic,jc,kp)+uy(ic,jc,kc)) &
               -(InterpY1*uz(ic,jm,kc)+InterpY2*uz(ic,jc,kc))* (uy(ic,jc,kc)+uy(ic,jc,km)) )*hdx3
                
          d11q2= ap1c(ic)*uy(ip,jc,kc)+ac1c(ic)*uy(ic,jc,kc)+am1c(ic)*uy(im,jc,kc)
          d22q2= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          d33q2= ap3c(kc)*uy(ic,jc,kp)+ac3c(kc)*uy(ic,jc,kc)+am3c(kc)*uy(ic,jc,km)
          dcq123= d11q2+d22q2+d33q2
       
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+ gravity(2) +InterpY1*FpForce_y(ic,jm,kc)+InterpY2*FpForce_y(ic,jc,kc)
#elif CFDLPT_TwoWay
#ifdef ForceInCell
          convEd2= -h21-h22-h23+ gravity(2) +InterpY1*FpForce_y(ic,jm,kc)+InterpY2*FpForce_y(ic,jc,kc)
#else
          convEd2= -h21-h22-h23+ gravity(2) +FpForce_y(ic,jc,kc)
#endif
#else
          convEd2= -h21-h22-h23+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2+ two*pmBeta*dcq123
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_FIMP

  !******************************************************************
  ! clcRhsY_FIMP_LES
  !******************************************************************    
  subroutine  clcRhsY_FIMP_LES(ux,uy,uz,RhsY,HistYold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::RhsY
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistYold
   
    ! locals 
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::qdx1,qdx3,h21,h22,h23,sucac,qsucac
    real(RK)::d11q2,d22q2,d33q2,dcq123,convEd2,gradp2  
    real(RK)::visa,visb,visc,visd,sgs1,sgs2,sgs3  
    
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=ystart(3),yend(3)
      km=kc-1
      kp=kc+1
      do jc=ystart(2),yend(2)
        jm=jc-1
        jp=jc+1
        sucac= rdyc(jc)
        qsucac=quarter*sucac
        do ic=ystart(1),yend(1)
          im=ic-1                
          ip=ic+1
          h21=( (ux(ip,jc,kc)+ux(ip,jm,kc))* (uy(ip,jc,kc)+uy(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jm,kc))* (uy(ic,jc,kc)+uy(im,jc,kc)) )*qdx1
          h22=( (uy(ic,jp,kc)+uy(ic,jc,kc))* (uy(ic,jp,kc)+uy(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jm,kc)) )*qsucac
          h23=( (uz(ic,jc,kp)+uz(ic,jm,kp))* (uy(ic,jc,kp)+uy(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jc,km)) )*qdx3

          visa=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc))
          visb=quarter*(nut(ic,jc,kc)+nut(ip,jc,kc)+nut(ic,jm,kc)+nut(ip,jm,kc))
          visc=quarter*(nut(ic,jc,kc)+nut(ic,jm,kc)+nut(ic,jc,km)+nut(ic,jm,km))
          visd=quarter*(nut(ic,jc,kc)+nut(ic,jm,kc)+nut(ic,jc,kp)+nut(ic,jm,kp))
          sgs1=rdx*sucac*(visb*(ux(ip,jc,kc)-ux(ip,jm,kc))-visa*(ux(ic,jc,kc)-ux(ic,jm,kc)))
          sgs2=ap2p(jc)*nut(ic,jc,kc)*(uy(ic,jp,kc)-uy(ic,jc,kc))-am2p(jc)*nut(ic,jm,kc)*(uy(ic,jc,kc)-uy(ic,jm,kc))
          sgs3=rdz*sucac*(visd*(uz(ic,jc,kp)-uz(ic,jm,kp))-visc*(uz(ic,jc,kc)-uz(ic,jm,kc)))

          d11q2= ap1c(ic)*visb*(uy(ip,jc,kc)-uy(ic,jc,kc) ) -am1c(ic)*visa*(uy(ic,jc,kc)- uy(im,jc,kc)) 
          d22q2= sgs2
          d33q2= ap3c(kc)*visd*(uy(ic,jc,kp)-uy(ic,jc,kc) ) -am3c(kc)*visc*(uy(ic,jc,kc)-uy(ic,jc,km))
          dcq123= d11q2+d22q2+d33q2
                
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+ sgs1+sgs2+sgs3+ gravity(2)+half*(FpForce_y(ic,jc,kc)+FpForce_y(ic,jm,kc))
#elif CFDLPT_TwoWay
#ifdef ForceInCell
          convEd2= -h21-h22-h23+ sgs1+sgs2+sgs3+ gravity(2)+half*(FpForce_y(ic,jc,kc)+FpForce_y(ic,jm,kc))
#else
          convEd2= -h21-h22-h23+ sgs1+sgs2+sgs3+ gravity(2)+FpForce_y(ic,jc,kc)
#endif
#else
          convEd2= -h21-h22-h23+ sgs1+sgs2+sgs3+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2+ two*pmBetaT*dcq123
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_FIMP_LES

  !******************************************************************
  ! clcRhsZ_FIMP
  !****************************************************************** 
  subroutine  clcRhsZ_FIMP(ux,uy,uz,RhsZ,HistZold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistZold
   
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::d11q3,d22q3,d33q3,dcq123,convEd3,gradp3
    real(RK)::qdx1,qdx3,h31,h32,h33,sucaj,InterpY1,InterpY2,InterpY3,InterpY4
    
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=ystart(3),yend(3)
      km=kc-1;kp=kc+1
      do jc=ystart(2),yend(2)
        jm=jc-1;jp=jc+1
        sucaj =half*rdyp(jc)
        InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=one-InterpY3
        do ic=ystart(1),yend(1)
          im=ic-1;ip=ic+1

          h31=( (ux(ip,jc,kc)+ux(ip,jc,km))* (uz(ip,jc,kc)+uz(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jc,km))* (uz(ic,jc,kc)+uz(im,jc,kc)) )*qdx1
          h32=( (uy(ic,jp,kc)+uy(ic,jp,km))* (InterpY3*uz(ic,jc,kc) +InterpY4*uz(ic,jp,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jc,km))* (InterpY1*uz(ic,jm,kc) +InterpY2*uz(ic,jc,kc)) )*sucaj                
          h33=( (uz(ic,jc,kp)+uz(ic,jc,kc))* (uz(ic,jc,kp)+uz(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jc,km)) )*qdx3

          d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
          d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
          d33q3= (uz(ic,jc,kp)-two*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          dcq123= d11q3+d22q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz                
#ifdef CFDDEM
          convEd3= -h31-h32-h33+ gravity(3)+half*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
#ifdef ForceInCell
          convEd3= -h31-h32-h33+ gravity(3)+half*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#else
          convEd3= -h31-h32-h33+ gravity(3)+FpForce_z(ic,jc,kc)
#endif
#else
          convEd3= -h31-h32-h33+ gravity(3)
#endif
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3+ two*pmBeta*dcq123
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_FIMP

  !******************************************************************
  ! clcRhsZ_FIMP_LES
  !****************************************************************** 
  subroutine  clcRhsZ_FIMP_LES(ux,uy,uz,RhsZ,HistZold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistZold
   
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::qdx1,qdx3,h31,h32,h33,sucaj
    real(RK)::d11q3,d22q3,d33q3,dcq123,convEd3,gradp3
    real(RK)::visa,visb,visc,visd,sgs1,sgs2,sgs3
    
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=ystart(3),yend(3)
      km=kc-1
      kp=kc+1
      do jc=ystart(2),yend(2)
        jm=jc-1
        jp=jc+1
        sucaj=quarter*rdyp(jc)
        do ic=ystart(1),yend(1)
          im=ic-1
          ip=ic+1
          h31=( (ux(ip,jc,kc)+ux(ip,jc,km))* (uz(ip,jc,kc)+uz(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jc,km))* (uz(ic,jc,kc)+uz(im,jc,kc)) )*qdx1
          h32=( (uy(ic,jp,kc)+uy(ic,jp,km))* (uz(ic,jp,kc)+uz(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jm,kc)) )*sucaj                
          h33=( (uz(ic,jc,kp)+uz(ic,jc,kc))* (uz(ic,jc,kp)+uz(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jc,km)) )*qdx3

          visa=quarter*(nut(ic,jc,km)+nut(im,jc,km)+nut(ic,jc,kc)+nut(im,jc,kc))
          visb=quarter*(nut(ic,jc,km)+nut(ip,jc,km)+nut(ic,jc,kc)+nut(ip,jc,kc))
          visc=quarter*(nut(ic,jc,km)+nut(ic,jm,km)+nut(ic,jm,kc)+nut(ic,jc,kc))
          visd=quarter*(nut(ic,jc,km)+nut(ic,jp,km)+nut(ic,jp,kc)+nut(ic,jc,kc))
          sgs1=rdx*rdz     *(visb*(ux(ip,jc,kc)-ux(ip,jc,km))-visa*(ux(ic,jc,kc)-ux(ic,jc,km)))
          sgs2=rdz*rdyp(jc)*(visd*(uy(ic,jp,kc)-uy(ic,jp,km))-visc*(uy(ic,jc,kc)-uy(ic,jc,km)))
          sgs3=rdz2*(nut(ic,jc,kc)*(uz(ic,jc,kp)-uz(ic,jc,kc))-nut(ic,jc,km)*(uz(ic,jc,kc)-uz(ic,jc,km)))

          d11q3= ap1c(ic)*visb*(uz(ip,jc,kc)-uz(ic,jc,kc)) -am1c(ic)*visa*(uz(ic,jc,kc)-uz(im,jc,kc)) 
          d22q3= ap2c(jc)*visd*(uz(ic,jp,kc)-uz(ic,jc,kc)) -am2c(jc)*visc*(uz(ic,jc,kc)-uz(ic,jm,kc))        
          d33q3= sgs3
          dcq123= d11q3+d22q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz 
#ifdef CFDDEM
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+gravity(3)+half*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
#ifdef ForceInCell
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+gravity(3)+half*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#else
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+gravity(3)+FpForce_z(ic,jc,kc)
#endif
#else
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+gravity(3)
#endif   
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3+ two*pmBetaT*dcq123
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_FIMP_LES

  !******************************************************************
  ! clcU1Hat_FIMP_000
  !******************************************************************    
  subroutine clcU1Hat_FIMP_000(ux,RhsX)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mjc,cjc,pjc,mkc,ckc,pkc
    real(RK),dimension(xsize(2),xsize(1))::tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y_to_x(RhsX,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      dO kc=1,zsize(3)  
        do ic=1,zsize(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO 

    ! compute dq2* sweeping in the y dirextion
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*two +one
    call transpose_z_to_y(arrz1,RhsX)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) =  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_000

  !******************************************************************
  ! clcU1Hat_FIMP_LES_000
  !******************************************************************    
  subroutine clcU1Hat_FIMP_LES_000(ux,RhsX)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::im,jm,km,ic,jc,kc,jp,kp
    real(RK):: rt1,rt2,betax,betay,betaz,visa,visb,visc,visd
    real(RK),dimension(xstart(2):xend(2),xstart(1):xend(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(zstart(1):zend(1),zstart(3):zend(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))::arrx1
    real(RK),dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    betax=pmBetaT*rdx2
    call transpose_y_to_x(RhsX,arrx1)
    DO kc=xstart(3),xend(3) 
      do jc=xstart(2),xend(2)
        do ic=xstart(1),xend(1)
          im=ic-1
          rt1=betax*nut(im,jc,kc)
          rt2=betax*nut(ic,jc,kc)
          tridmi(jc,ic)= -rt1
          tridci(jc,ic)=  rt1+rt2+one      
          tridpi(jc,ic)= -rt2
          tridfi(jc,ic)=  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmi,tridci,tridpi,tridfi,xsize(2),nxc)
      do jc=xstart(2),xend(2)
        do ic=xstart(1),xend(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    betaz=quarter*pmBetaT*rdz2
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=zstart(2),zend(2)
      do kc=zstart(3),zend(3)  
        kp=kc+1
        km=kc-1
        do ic=zstart(1),zend(1)
          im=ic-1
          visa=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,km)+nut(im,jc,km)
          visb=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,kp)+nut(im,jc,kp)
          rt1=betaz*visa
          rt2=betaz*visb
          tridmk(ic,kc)= -rt1        
          tridck(ic,kc)=  rt1+rt2+one      
          tridpk(ic,kc)= -rt2 
          tridfk(ic,kc)=  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmk,tridck,tridpk,tridfk,zsize(1),nzc)
      do kc=zstart(3),zend(3)
        do ic=zstart(1),zend(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO 

    ! compute dq2* sweeping in the y dirextion
    betay=quarter*pmBetaT*rdy2
    call transpose_z_to_y(arrz1,Rhsx)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        jm=jc-1
        jp=jc+1
        do ic=ystart(1),yend(1)
          im=ic-1
          visc=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc)
          visd=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jp,kc)+nut(im,jp,kc)
          rt1=betay*visc
          rt2=betay*visd
          tridmj(ic,jc)= -rt1        
          tridcj(ic,jc)=  rt1+rt2+one      
          tridpj(ic,jc)= -rt2 
          tridfj(ic,jc)=  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_LES_000

  !******************************************************************
  ! clcU1Hat_FIMP_010
  !******************************************************************    
  subroutine clcU1Hat_FIMP_010(ux,RhsX)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc
    real(RK),dimension(xsize(2),xsize(1))::tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y_to_x(RhsX,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      dO kc=1,zsize(3)  
        do ic=1,zsize(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+one
      enddo
    enddo 
    call transpose_z_to_y(arrz1,Rhsx)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) = RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_010

  !******************************************************************
  ! clcU1Hat_FIMP_LES_010
  !******************************************************************    
  subroutine clcU1Hat_FIMP_LES_010(ux,RhsX)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::im,jm,km,ic,jc,kc,jp,kp
    real(RK):: rt1,rt2,betax,betaz,visa,visb,visc,visd
    real(RK),dimension(xstart(2):xend(2),xstart(1):xend(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(zstart(1):zend(1),zstart(3):zend(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))::arrx1
    real(RK),dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    betax=pmBetaT*rdx2
    call transpose_y_to_x(RhsX,arrx1)
    DO kc=xstart(3),xend(3) 
      do jc=xstart(2),xend(2)
        do ic=xstart(1),xend(1)
          im=ic-1
          rt1=betax*nut(im,jc,kc)
          rt2=betax*nut(ic,jc,kc)
          tridmi(jc,ic)= -rt1        
          tridci(jc,ic)=  rt1+rt2+one      
          tridpi(jc,ic)= -rt2
          tridfi(jc,ic)=  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmi,tridci,tridpi,tridfi,xsize(2),nxc)
      do jc=xstart(2),xend(2)
        do ic=xstart(1),xend(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    betaz=quarter*pmBetaT*rdz2
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=zstart(2),zend(2)
      do kc=zstart(3),zend(3)  
        kp=kc+1
        km=kc-1
        do ic=zstart(1),zend(1)
          im=ic-1
          visa=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,km)+nut(im,jc,km)
          visb=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,kp)+nut(im,jc,kp)
          rt1=betaz*visa
          rt2=betaz*visb
          tridmk(ic,kc)= -rt1        
          tridck(ic,kc)=  rt1+rt2+one      
          tridpk(ic,kc)= -rt2 
          tridfk(ic,kc)=  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmk,tridck,tridpk,tridfk,zsize(1),nzc)
      do kc=zstart(3),zend(3)
        do ic=zstart(1),zend(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO 

    ! compute dq2* sweeping in the y dirextion
    call transpose_z_to_y(arrz1,RhsX)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        jm=jc-1
        jp=jc+1
        do ic=ystart(1),yend(1)
          im=ic-1
          visc= quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc))
          visd= quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jp,kc)+nut(im,jp,kc))
          rt1 = pmBetaT*visc*am2cForCN(jc)
          rt2 = pmBetaT*visd*ap2cForCN(jc)
          tridmj(ic,jc) = -rt1 
          tridcj(ic,jc) =  rt1+rt2+one
          tridpj(ic,jc) = -rt2
          tridfj(ic,jc) =  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_LES_010

  !******************************************************************
  ! clcU1Hat_FIMP_011
  !******************************************************************    
  subroutine clcU1Hat_FIMP_011(ux,RhsX)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic
    real(RK),dimension(xsize(2),xsize(1))::tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y_to_x(RhsX,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    do kc=1,zsize(3)  
      do ic=1,zsize(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+one
      enddo
    enddo
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      do kc=1,zsize(3)  
        do ic=1,zsize(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+one
      enddo
    enddo
    call transpose_z_to_y(arrz1,Rhsx)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) = RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_011

  !******************************************************************
  ! clcU1Hat_FIMP_111
  !******************************************************************    
  subroutine clcU1Hat_FIMP_111(ux,RhsX)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic
    real(RK),dimension(xsize(2),xsize(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    do jc=1,xsize(2)
      tridmi(jc,1) = zero
      tridci(jc,1) = one
      tridpi(jc,1) = zero
      do ic=2,xsize(1)
        tridmi(jc,ic) =  mic
        tridci(jc,ic) =  cic
        tridpi(jc,ic) =  pic
      enddo
    enddo
    call transpose_y_to_x(RhsX,arrx1)
    DO kc=1,xsize(3)
      do jc=1,xsize(2)
        tridfi(jc,1) = zero
        do ic=2,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmi,tridci,tridpi,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    do kc=1,zsize(3)  
      do ic=1,zsize(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+one
      enddo
    enddo
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      do kc=1,zsize(3)  
        do ic=1,zsize(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+one
      enddo
    enddo 
    call transpose_z_to_y(arrz1,Rhsx)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) = RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_111

  !******************************************************************
  ! clcU2Hat_FIMP_000
  !******************************************************************    
  subroutine clcU2Hat_FIMP_000(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in):: duy_ym    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mjc,cjc,pjc,mkc,ckc,pkc
    real(RK),dimension(xsize(2),xsize(1))::tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y_to_x(RhsY,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      dO kc=1,zsize(3)  
        do ic=1,zsize(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO 

    ! compute dq2* sweeping in the y dirextion
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*two +one
    call transpose_z_to_y(arrz1,RhsY)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) =  RhsY(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU2Hat_FIMP_000

  !******************************************************************
  ! clcU2Hat_FIMP_LES_000
  !******************************************************************    
  subroutine clcU2Hat_FIMP_LES_000(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in):: duy_ym    
    
    ! locals
    integer:: im,jm,km,ic,jc,kc,ip,kp
    real(RK)::rt1,rt2,betax,betay,betaz,visa,visb,visc,visd
    real(RK),dimension(xstart(2):xend(2),xstart(1):xend(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(zstart(1):zend(1),zstart(3):zend(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))::arrx1
    real(RK),dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    betax=quarter*pmBetaT*rdx2
    call transpose_y_to_x(RhsY,arrx1)
    DO kc=xstart(3),xend(3) 
      do jc=xstart(2),xend(2)
        jm=jc-1
        do ic=xstart(1),xend(1)
          ip=ic+1
          im=ic-1
          visa=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc)
          visb=nut(ic,jc,kc)+nut(ip,jc,kc)+nut(ic,jm,kc)+nut(ip,jm,kc)
          rt1=betax*visa
          rt2=betax*visb
          tridmi(jc,ic)= -rt1        
          tridci(jc,ic)=  rt1+rt2+one      
          tridpi(jc,ic)= -rt2
          tridfi(jc,ic)=  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmi,tridci,tridpi,tridfi,xsize(2),nxc)
      do jc=xstart(2),xend(2)
        do ic=xstart(1),xend(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO

    ! compute dq3* sweeping in the z dirextion
    betaz=quarter*pmBetaT*rdz2
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=zstart(2),zend(2)
      do kc=zstart(3),zend(3)  
        kp=kc+1
        km=kc-1
        do ic=zstart(1),zend(1)
          im=ic-1
          visc=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,km)+nut(im,jc,km)
          visd=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,kp)+nut(im,jc,kp)
          rt1=betaz*visc
          rt2=betaz*visd
          tridmk(ic,kc)= -rt1        
          tridck(ic,kc)=  rt1+rt2+one      
          tridpk(ic,kc)= -rt2 
          tridfk(ic,kc)=  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmk,tridck,tridpk,tridfk,zsize(1),nzc)
      do kc=zstart(3),zend(3)
        do ic=zstart(1),zend(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    betay=pmBetaT*rdy2
    call transpose_z_to_y(arrz1,RhsY)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        jm=jc-1
        do ic=ystart(1),yend(1)
          rt1=betay*nut(ic,jm,kc)
          rt2=betay*nut(ic,jc,kc)
          tridmj(ic,jc)= -rt1        
          tridcj(ic,jc)=  rt1+rt2+one      
          tridpj(ic,jc)= -rt2 
          tridfj(ic,jc)=  RhsY(ic,jc,kc)
        enddo
      enddo  
      call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU2Hat_FIMP_LES_000

  !******************************************************************
  ! clcU2Hat_FIMP_010
  !******************************************************************    
  subroutine clcU2Hat_FIMP_010(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in):: duy_ym    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc
    real(RK),dimension(xsize(2),xsize(1))::tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y_to_x(RhsY,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      dO kc=1,zsize(3)  
        do ic=1,zsize(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    call transpose_z_to_y(arrz1,RhsY)
    do ic=ystart(1),yend(1) 
      tridpj(ic,1)=zero
      tridcj(ic,1)=one
      tridmj(ic,1)=zero
    enddo
    do jc=2,nyc
      do ic=ystart(1),yend(1) 
        tridpj(ic,jc) = -pmBeta*ap2p(jc)
        tridcj(ic,jc) = -pmBeta*ac2p(jc)+one
        tridmj(ic,jc) = -pmBeta*am2p(jc)
      enddo
    enddo
    DO kc=ystart(3),yend(3) 
      do ic=ystart(1),yend(1)
        tridfj(ic,1)=duy_ym(ic,kc)
      enddo
      do jc=2,nyc
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=1,nyc
         do ic=ystart(1),yend(1) 
           uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_FIMP_010

  !******************************************************************
  ! clcU2Hat_FIMP_LES_010
  !******************************************************************    
  subroutine clcU2Hat_FIMP_LES_010(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in):: duy_ym    
    
    ! locals
    integer:: im,jm,km,ic,jc,kc,ip,kp
    real(RK)::rt1,rt2,betax,betaz,visa,visb,visc,visd
    real(RK),dimension(xstart(2):xend(2),xstart(1):xend(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(zstart(1):zend(1),zstart(3):zend(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))::arrx1
    real(RK),dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    betax=quarter*pmBetaT*rdx2
    call transpose_y_to_x(RhsY,arrx1)
    DO kc=xstart(3),xend(3) 
      do jc=xstart(2),xend(2)
        jm=jc-1
        do ic=xstart(1),xend(1)
          ip=ic+1
          im=ic-1
          visa=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc)
          visb=nut(ic,jc,kc)+nut(ip,jc,kc)+nut(ic,jm,kc)+nut(ip,jm,kc)
          rt1=betax*visa
          rt2=betax*visb
          tridmi(jc,ic)= -rt1        
          tridci(jc,ic)=  rt1+rt2+one      
          tridpi(jc,ic)= -rt2
          tridfi(jc,ic)=  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmi,tridci,tridpi,tridfi,xsize(2),nxc)
      do jc=xstart(2),xend(2)
        do ic=xstart(1),xend(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO

    ! compute dq3* sweeping in the z dirextion
    betaz=quarter*pmBetaT*rdz2
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=zstart(2),zend(2)
      do kc=zstart(3),zend(3)  
        kp=kc+1
        km=kc-1
        do ic=zstart(1),zend(1)
          im=ic-1
          visc=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,km)+nut(im,jc,km)
          visd=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,kp)+nut(im,jc,kp)
          rt1=betaz*visc
          rt2=betaz*visd
          tridmk(ic,kc)= -rt1        
          tridck(ic,kc)=  rt1+rt2+one      
          tridpk(ic,kc)= -rt2 
          tridfk(ic,kc)=  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmk,tridck,tridpk,tridfk,zsize(1),nzc)
      do kc=zstart(3),zend(3)
        do ic=zstart(1),zend(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    call transpose_z_to_y(arrz1,RhsY)
    DO kc=ystart(3),yend(3) 
      do ic=ystart(1),yend(1) 
        tridpj(ic,1)=zero
        tridcj(ic,1)=one
        tridmj(ic,1)=zero
        tridfj(ic,1)=duy_ym(ic,kc)
      enddo
      do jc=2,nyc
        jm=jc-1
        do ic=ystart(1),yend(1) 
          rt1= pmBetaT*nut(ic,jm,kc)*am2p(jc)
          rt2= pmBetaT*nut(ic,jc,kc)*ap2p(jc)
          tridmj(ic,jc) = -rt1
          tridcj(ic,jc) =  rt1+rt2+one
          tridpj(ic,jc) = -rt2
          tridfj(ic,jc) =  RhsY(ic,jc,kc)
        enddo
      enddo
     
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=1,nyc
        do ic=ystart(1),yend(1) 
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_FIMP_LES_010

  !******************************************************************
  ! clcU2Hat_FIMP_011
  !******************************************************************    
  subroutine clcU2Hat_FIMP_011(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in):: duy_ym    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic
    real(RK),dimension(xsize(2),xsize(1))::tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y_to_x(RhsY,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO

    ! compute dq3* sweeping in the z dirextion
    do kc=1,zsize(3)  
      do ic=1,zsize(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+one
      enddo
    enddo
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      do kc=1,zsize(3)  
        do ic=1,zsize(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    do ic=ystart(1),yend(1) 
      tridpj(ic,1)=zero
      tridcj(ic,1)=one
      tridmj(ic,1)=zero
    enddo
    do jc=2,nyc
      do ic=ystart(1),yend(1) 
        tridpj(ic,jc) = -pmBeta*ap2p(jc)
        tridcj(ic,jc) = -pmBeta*ac2p(jc)+one
        tridmj(ic,jc) = -pmBeta*am2p(jc)
      enddo
    enddo
    call transpose_z_to_y(arrz1,RhsY)
    DO kc=ystart(3),yend(3) 
      do ic=ystart(1),yend(1)
        tridfj(ic,1)=duy_ym(ic,kc)
      enddo
      do jc=2,nyc
        do ic=ystart(1),yend(1) 
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=1,nyc
        do ic=ystart(1),yend(1) 
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_FIMP_011

  !******************************************************************
  ! clcU2Hat_FIMP_111
  !******************************************************************    
  subroutine clcU2Hat_FIMP_111(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in):: duy_ym    
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(xsize(2),xsize(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    
    ! compute dq1* weeping in the x dirextion
    do jc=1,xsize(2)
      do ic=1,xsize(1)
        tridpi(jc,ic)= -pmBeta*ap1cForCN(ic)
        tridmi(jc,ic)= -pmBeta*am1cForCN(ic)
        tridci(jc,ic)= -tridpi(jc,ic)-tridmi(jc,ic)+one
      enddo
    enddo
    call transpose_y_to_x(RhsY,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic)=  arrx1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmi,tridci,tridpi,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO

    ! compute dq3* sweeping in the z dirextion
    do kc=1,zsize(3)  
      do ic=1,zsize(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+one
      enddo
    enddo
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      do kc=1,zsize(3)  
        do ic=1,zsize(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    do ic=ystart(1),yend(1) 
      tridpj(ic,1)=zero
      tridcj(ic,1)=one
      tridmj(ic,1)=zero
    enddo
    do jc=2,nyc
      do ic=ystart(1),yend(1) 
        tridpj(ic,jc) = -pmBeta*ap2p(jc)
        tridcj(ic,jc) = -pmBeta*ac2p(jc)+one
        tridmj(ic,jc) = -pmBeta*am2p(jc)
      enddo
    enddo
    call transpose_z_to_y(arrz1,RhsY)
    DO kc=ystart(3),yend(3) 
      do ic=ystart(1),yend(1)
        tridfj(ic,1)=duy_ym(ic,kc)
      enddo
      do jc=2,nyc
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=1,nyc
        do ic=ystart(1),yend(1) 
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_FIMP_111

  !******************************************************************
  ! clcU3Hat_FIMP_000
  !******************************************************************    
  subroutine clcU3Hat_FIMP_000(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mjc,cjc,pjc,mkc,ckc,pkc
    real(RK),dimension(xsize(2),xsize(1))::tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))::RhsZ_temp

    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y_to_x(RhsZ_temp,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      dO kc=1,zsize(3)  
        do ic=1,zsize(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO 

    ! compute dq2* sweeping in the y dirextion
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*two +one
    call transpose_z_to_y(arrz1,RhsZ_temp)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) =  RhsZ_temp(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_000

  !******************************************************************
  ! clcU3Hat_FIMP_LES_000
  !******************************************************************    
  subroutine clcU3Hat_FIMP_LES_000(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::im,jm,km,ic,jc,kc,ip,jp
    real(RK):: rt1,rt2,betax,betay,betaz,visa,visb,visc,visd
    real(RK),dimension(xstart(2):xend(2),xstart(1):xend(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(zstart(1):zend(1),zstart(3):zend(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))::arrx1
    real(RK),dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))::arrz1
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))::RhsZ_temp

    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* weeping in the x dirextion
    betax=quarter*pmBetaT*rdx2
    call transpose_y_to_x(RhsZ_temp,arrx1)
    DO kc=xstart(3),xend(3) 
      do jc=xstart(2),xend(2)
        jm=jc-1
        do ic=xstart(1),xend(1)
          ip=ic+1
          im=ic-1
          visa=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc)
          visb=nut(ic,jc,kc)+nut(ip,jc,kc)+nut(ic,jm,kc)+nut(ip,jm,kc)
          rt1=betax*visa
          rt2=betax*visb
          tridmi(jc,ic)= -rt1        
          tridci(jc,ic)=  rt1+rt2+one      
          tridpi(jc,ic)= -rt2
          tridfi(jc,ic)=  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmi,tridci,tridpi,tridfi,xsize(2),nxc)
      do jc=xstart(2),xend(2)
        do ic=xstart(1),xend(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO

    ! compute dq3* sweeping in the z dirextion
    betaz=pmBetaT*rdz2
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=zstart(2),zend(2)
      do kc=zstart(3),zend(3)
        km=kc-1
        do ic=zstart(1),zend(1)
          rt1=betaz*nut(ic,jc,km)
          rt2=betaz*nut(ic,jc,kc)
          tridmk(ic,kc)= -rt1        
          tridck(ic,kc)=  rt1+rt2+one      
          tridpk(ic,kc)= -rt2 
          tridfk(ic,kc)=  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmk,tridck,tridpk,tridfk,zsize(1),nzc)
      do kc=zstart(3),zend(3)
        do ic=zstart(1),zend(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    betay=quarter*pmBetaT*rdy2
    call transpose_z_to_y(arrz1,RhsZ_temp)
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        jm=jc-1
        jp=jc+1
        do ic=ystart(1),yend(1)
          im=ic-1
          visc=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc)
          visd=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jp,kc)+nut(im,jp,kc)
          rt1=betay*visc
          rt2=betay*visd
          tridmj(ic,jc)= -rt1        
          tridcj(ic,jc)=  rt1+rt2+one      
          tridpj(ic,jc)= -rt2 
          tridfj(ic,jc)=  RhsZ_temp(ic,jc,kc)
        enddo
      enddo  
      call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_LES_000

  !******************************************************************
  ! clcU3Hat_FIMP_010
  !******************************************************************    
  subroutine clcU3Hat_FIMP_010(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc
    real(RK),dimension(xsize(2),xsize(1))::tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))::RhsZ_temp

    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y_to_x(RhsZ_temp,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      dO kc=1,zsize(3)  
        do ic=1,zsize(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO 

    ! compute dq2* sweeping in the y dirextion
    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+one
      enddo
    enddo
    call transpose_z_to_y(arrz1,RhsZ_temp)
    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) = RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_010

  !******************************************************************
  ! clcU3Hat_FIMP_LES_010
  !******************************************************************    
  subroutine clcU3Hat_FIMP_LES_010(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::im,jm,km,ic,jc,kc,ip,jp
    real(RK):: rt1,rt2,betax,betaz,visa,visb,visc,visd
    real(RK),dimension(xstart(2):xend(2),xstart(1):xend(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(zstart(1):zend(1),zstart(3):zend(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))::arrx1
    real(RK),dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))::arrz1
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))::RhsZ_temp

    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* weeping in the x dirextion
    betax=quarter*pmBetaT*rdx2
    call transpose_y_to_x(RhsZ_temp,arrx1)
    DO kc=xstart(3),xend(3) 
      do jc=xstart(2),xend(2)
        jm=jc-1
        do ic=xstart(1),xend(1)
          ip=ic+1
          im=ic-1
          visa=nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc)
          visb=nut(ic,jc,kc)+nut(ip,jc,kc)+nut(ic,jm,kc)+nut(ip,jm,kc)
          rt1=betax*visa
          rt2=betax*visb
          tridmi(jc,ic)= -rt1        
          tridci(jc,ic)=  rt1+rt2+one      
          tridpi(jc,ic)= -rt2
          tridfi(jc,ic)=  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmi,tridci,tridpi,tridfi,xsize(2),nxc)
      do jc=xstart(2),xend(2)
        do ic=xstart(1),xend(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO

    ! compute dq3* sweeping in the z dirextion
    betaz=pmBetaT*rdz2
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=zstart(2),zend(2)
      do kc=zstart(3),zend(3)
        km=kc-1
        do ic=zstart(1),zend(1)
          rt1=betaz*nut(ic,jc,km)
          rt2=betaz*nut(ic,jc,kc)
          tridmk(ic,kc)= -rt1        
          tridck(ic,kc)=  rt1+rt2+one      
          tridpk(ic,kc)= -rt2 
          tridfk(ic,kc)=  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmk,tridck,tridpk,tridfk,zsize(1),nzc)
      do kc=zstart(3),zend(3)
        do ic=zstart(1),zend(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    call transpose_z_to_y(arrz1,RhsZ_temp)
    DO kc=ystart(3),yend(3)
      km=kc-1
      do jc=ystart(2),yend(2)
        jm=jc-1
        jp=jc+1
        do ic=ystart(1),yend(1)
          visc=quarter*(nut(ic,jc,km)+nut(ic,jm,km)+nut(ic,jc,kc)+nut(ic,jm,kc))
          visd=quarter*(nut(ic,jc,km)+nut(ic,jp,km)+nut(ic,jc,kc)+nut(ic,jp,kc))
          rt1= pmBetaT*visc*am2cForCN(jc)
          rt2= pmBetaT*visd*ap2cForCN(jc)
          tridmj(ic,jc) = -rt1
          tridcj(ic,jc) =  rt1+rt2+one
          tridpj(ic,jc) = -rt2
          tridfj(ic,jc) =  RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_LES_010

  !******************************************************************
  ! clcU3Hat_FIMP_011
  !******************************************************************    
  subroutine clcU3Hat_FIMP_011(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc
    real(RK),dimension(xsize(2),xsize(1))::tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))::RhsZ_temp

    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* weeping in the x dirextion
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y_to_x(RhsZ_temp,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    do ic=1,zsize(1)
      tridpk(ic,1) =  zero
      tridck(ic,1) =  one
      tridmk(ic,1) =  zero
    enddo
    dO kc=2,zsize(3)  
      do ic=1,zsize(1)
        tridpk(ic,kc) =  pkc
        tridck(ic,kc) =  ckc
        tridmk(ic,kc) =  mkc
      enddo
    enddo
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      do ic=1,zsize(1)
        tridfk(ic,1) =  zero
      enddo
      dO kc=2,zsize(3)  
        do ic=1,zsize(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk,tridck,tridpk,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+one
      enddo
    enddo
    call transpose_z_to_y(arrz1,RhsZ_temp)
    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) = RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_011

  !******************************************************************
  ! clcU3Hat_FIMP_111 (unfinished !!!)
  !******************************************************************    
  subroutine clcU3Hat_FIMP_111(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mkc,ckc,pkc
    real(RK),dimension(xsize(2),xsize(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(zsize(1),zsize(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(xsize(1),xsize(2),xsize(3))::arrx1
    real(RK),dimension(zsize(1),zsize(2),zsize(3))::arrz1
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))::RhsZ_temp

    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* weeping in the x dirextion
    do jc=1,xsize(2)
      do ic=1,xsize(1)
        tridpi(jc,ic)= -pmBeta*ap1cForCN(ic)
        tridmi(jc,ic)= -pmBeta*am1cForCN(ic)
        tridci(jc,ic)= -tridpi(jc,ic)-tridmi(jc,ic)+one
      enddo
    enddo
    call transpose_y_to_x(RhsZ_temp,arrx1)
    DO kc=1,xsize(3)  
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmi,tridci,tridpi,tridfi,xsize(2),nxc)
      do jc=1,xsize(2)
        do ic=1,xsize(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 

    ! compute dq3* sweeping in the z dirextion
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    do ic=1,zsize(1)
      tridpk(ic,1) =  zero
      tridck(ic,1) =  one
      tridmk(ic,1) =  zero
    enddo
    dO kc=2,zsize(3)  
      do ic=1,zsize(1)
        tridpk(ic,kc) =  pkc
        tridck(ic,kc) =  ckc
        tridmk(ic,kc) =  mkc
      enddo
    enddo
    call transpose_x_to_z(arrx1,arrz1)
    Do jc=1,zsize(2)
      do ic=1,zsize(1)
        tridfk(ic,1) =  zero
      enddo
      dO kc=2,zsize(3)  
        do ic=1,zsize(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk,tridck,tridpk,tridfk,zsize(1),nzc)
      do kc=1,zsize(3)
        do ic=1,zsize(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO

    ! compute dq2* sweeping in the y dirextion
    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+one
      enddo
    enddo
    call transpose_z_to_y(arrz1,RhsZ_temp)
    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) = RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_111

  !******************************************************************
  ! PressureUpdate_FIMP
  !******************************************************************
  subroutine PressureUpdate_FIMP(pressure, prphiHalo)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out):: pressure
    
    ! locals
    integer::ic,jc,kc

    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          pressure(ic,jc,kc)= pressure(ic,jc,kc)+ prphiHalo(ic,jc,kc)
        enddo
      enddo
    ENDDO
  end subroutine PressureUpdate_FIMP

  !******************************************************************
  ! PressureUpdate_FIMP_LES
  !******************************************************************
  subroutine PressureUpdate_FIMP_LES(pressure, prphiHalo)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out):: pressure
    
    !locals
    integer::im,jm,km,ic,jc,kc,ip,jp,kp
    real(RK)::visa,visb,visc,visd,vise,visf,metr1,metr2

    DO kc=ystart(3),yend(3)
      kp=kc+1
      km=kc-1
      do jc=ystart(2),yend(2)
        jp=jc+1
        jm=jc-1
        metr1=rdyp(jc)*rdyc(jc)
        metr2=rdyp(jc)*rdyc(jp)
        do ic=ystart(1),yend(1)
          ip=ic+1
          im=ic-1
          visa= half*(nut(ic,jc,kc)+nut(im,jc,kc))*rdx2
          visb= half*(nut(ic,jc,kc)+nut(ip,jc,kc))*rdx2
          visc= half*(nut(ic,jc,kc)+nut(ic,jm,kc))*metr1
          visd= half*(nut(ic,jp,kc)+nut(ic,jc,kc))*metr2
          vise= half*(nut(ic,jc,kc)+nut(ic,jc,km))*rdz2
          visf= half*(nut(ic,jc,kc)+nut(ic,jc,kp))*rdz2
          pressure(ic,jc,kc)=  pressure(ic,jc,kc)+ prphiHalo(ic,jc,kc)- pmBetaT* &
                              (visb*prphiHalo(ip,jc,kc)-(visb+visa)*prphiHalo(ic,jc,kc)+visa*prphiHalo(im,jc,kc)+ &
                               visd*prphiHalo(ic,jp,kc)-(visc+visd)*prphiHalo(ic,jc,kc)+visc*prphiHalo(ic,jm,kc)+ &
                               visf*prphiHalo(ic,jc,kp)-(visf+vise)*prphiHalo(ic,jc,kc)+vise*prphiHalo(ic,jc,km)  )

        enddo
      enddo
    ENDDO
  end subroutine PressureUpdate_FIMP_LES
