module m_TScheme
  use MPI
  use decomp_2d
  use m_LogInfo
  use m_TypeDef  
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1,hi1
  use m_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal,InversePTriFixedCoe
  implicit none
  private

  type(HaloInfo),save::hi_uxhx,hi_uyhx,hi_uzhx,hi_uxhz,hi_uyhz,hi_uzhz
  public:: InitTimeScheme,FluidVelUpdate,clcRhsXYZ
  public:: clcPrSrc,clcU1Hat,clcU2Hat,clcU3Hat,PressureUpdate
contains    

  !******************************************************************
  ! InitTimeScheme
  !******************************************************************
  subroutine InitTimeScheme()
    implicit none
  
    ! Halo info
    hi_uxhx%pencil= y_pencil; hi_uyhx%pencil= y_pencil; hi_uzhx%pencil= y_pencil
    hi_uxhz%pencil= y_pencil; hi_uyhz%pencil= y_pencil; hi_uzhz%pencil= y_pencil
    hi_uxhx%xmh=2;hi_uxhx%xph=1;hi_uxhx%ymh=0;hi_uxhx%yph=0;hi_uxhx%zmh=0;hi_uxhx%zph=0
    hi_uyhx%xmh=1;hi_uyhx%xph=2;hi_uyhx%ymh=0;hi_uyhx%yph=0;hi_uyhx%zmh=0;hi_uyhx%zph=0
    hi_uzhx%xmh=1;hi_uzhx%xph=2;hi_uzhx%ymh=0;hi_uzhx%yph=0;hi_uzhx%zmh=1;hi_uzhx%zph=2
    hi_uxhz%xmh=1;hi_uxhz%xph=2;hi_uxhz%ymh=0;hi_uxhz%yph=0;hi_uxhz%zmh=1;hi_uxhz%zph=2
    hi_uyhz%xmh=0;hi_uyhz%xph=0;hi_uyhz%ymh=0;hi_uyhz%yph=0;hi_uyhz%zmh=1;hi_uyhz%zph=2
    hi_uzhz%xmh=0;hi_uzhz%xph=0;hi_uzhz%ymh=0;hi_uzhz%yph=0;hi_uzhz%zmh=2;hi_uzhz%zph=1
  end subroutine InitTimeScheme

  !******************************************************************
  ! clcPrSrc
  !******************************************************************  
  subroutine clcPrSrc(ux,uy,uz,prsrc,pressure,divmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
    real(RK),intent(out)::divmax
    
    ! locals
    integer::ic,jc,kc,ip,jp,kp,ierror,km,im,ku,iu
    real(RK)::sudtal,sucaj,rdiv,divmax1,dudx,dvdy,dwdz

    divmax1=zero
    sudtal=one/pmAlpha
    DO kc=ystart(3),yend(3)
       km=kc-1;kp=kc+1;ku=kc+2
       do jc=ystart(2),yend(2)
         jp=jc+1;sucaj=rdyp(jc)
         do ic=ystart(1),yend(1)
           im=ic-1;ip=ic+1;iu=ic+2
           dudx=  dxCoe1*ux(im,jc,kc) +dxCoe2*ux(ic,jc,kc) +dxCoe3*ux(ip,jc,kc) +dxCoe4*ux(iu,jc,kc)
           dvdy=  (uy(ic,jp,kc)-uy(ic,jc,kc))*sucaj
           dwdz=  dzCoe1*uz(ic,jc,km) +dzCoe2*uz(ic,jc,kc) +dzCoe3*uz(ic,jc,kp) +dzCoe4*uz(ic,jc,ku)

           rdiv= dudx +dvdy + dwdz
           divmax1=max(abs(rdiv),divmax1)
           prsrc(ic,jc,kc)= sudtal * rdiv
         enddo
       enddo
     ENDDO
     call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)
  end subroutine clcPrSrc

  !******************************************************************
  ! FluidVelUpdate
  !******************************************************************
  subroutine FluidVelUpdate(prphiHalo,ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out):: ux,uy,uz
    
    ! locals
    real(RK)::sucac,dpdx,dpdy,dpdz
    integer::ic,jc,kc,im,jm,km,is,ks,ip,kp

    DO kc=ystart(3),yend(3)
      ks=kc-2;km=kc-1;kp=kc+1
      do jc=ystart(2),yend(2)
        jm=jc-1
        sucac=rdyc(jc)
        do ic=ystart(1),yend(1)
          is=ic-2;im=ic-1;ip=ic+1
          dpdx= dxCoe1*prphiHalo(is,jc,kc) +dxCoe2*prphiHalo(im,jc,kc) +dxCoe3*prphiHalo(ic,jc,kc) +dxCoe4*prphiHalo(ip,jc,kc)
          dpdy=(prphiHalo(ic,jc,kc)-prphiHalo(ic,jm,kc))*sucac
          dpdz= dzCoe1*prphiHalo(ic,jc,ks) +dzCoe2*prphiHalo(ic,jc,km) +dzCoe3*prphiHalo(ic,jc,kc) +dzCoe4*prphiHalo(ic,jc,kp)
          ux(ic,jc,kc)=ux(ic,jc,kc)-  pmAlpha*dpdx
          uy(ic,jc,kc)=uy(ic,jc,kc)-  pmAlpha*dpdy
          uz(ic,jc,kc)=uz(ic,jc,kc)-  pmAlpha*dpdz               
        enddo
      enddo
    ENDDO
  end subroutine FluidVelUpdate

  !******************************************************************
  ! clcRhsXYZ
  !******************************************************************    
  subroutine clcRhsXYZ(ux,uy,uz,RhsX,RhsY,RhsZ,HistXOld,HistYOld,HistZOld,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(out)::RhsX,RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::HistXOld,HistYOld,HistZOld
    
    ! locals    
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,is,ks,iu,ku,ierror
    real(RK)::d11q1,d22q1,d33q1,d11q2,d22q2,d33q2,d11q3,d22q3,d33q3,dcq13
    real(RK)::uxhyp,uxhyc,uxhym,uxhyu,uzhym,uzhyc,uzhyp,uzhyu,gradp1,gradp2,gradp3
    real(RK)::sucaj,s3tot,s3tot1,dp1ns,dpmdxns,convEd1,convEd2,convEd3,sucac,qsucac
    real(RK)::h11,h12,h13,h21,h22,h23,h31,h32,h33,InterpY1,InterpY2,InterpY3,InterpY4
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm)::uxhx,uxhz,uyhx,uyhz,uzhx,uzhz

    ! Update_VelInterp_uvw
    DO kc=ystart(3),yend(3)
      ks=kc-2;km=kc-1;kp=kc+1;ku=kc+2
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
           is=ic-2;im=ic-1;ip=ic+1;iu=ic+2
           uxhx(ic,jc,kc)= InterpCoe1*ux(im,jc,kc) +InterpCoe2*ux(ic,jc,kc) +InterpCoe3*ux(ip,jc,kc) +InterpCoe4*ux(iu,jc,kc)
           uyhx(ic,jc,kc)= InterpCoe1*uy(is,jc,kc) +InterpCoe2*uy(im,jc,kc) +InterpCoe3*uy(ic,jc,kc) +InterpCoe4*uy(ip,jc,kc)
           uzhx(ic,jc,kc)= InterpCoe1*uz(is,jc,kc) +InterpCoe2*uz(im,jc,kc) +InterpCoe3*uz(ic,jc,kc) +InterpCoe4*uz(ip,jc,kc)

           uxhz(ic,jc,kc)= InterpCoe1*ux(ic,jc,ks) +InterpCoe2*ux(ic,jc,km) +InterpCoe3*ux(ic,jc,kc) +InterpCoe4*ux(ic,jc,kp)
           uyhz(ic,jc,kc)= InterpCoe1*uy(ic,jc,ks) +InterpCoe2*uy(ic,jc,km) +InterpCoe3*uy(ic,jc,kc) +InterpCoe4*uy(ic,jc,kp)
           uzhz(ic,jc,kc)= InterpCoe1*uz(ic,jc,km) +InterpCoe2*uz(ic,jc,kc) +InterpCoe3*uz(ic,jc,kp) +InterpCoe4*uz(ic,jc,ku)
        enddo
      enddo
      jc=nyp
      do ic=ystart(1),yend(1)
        uyhx(ic,jc,kc)=zero; uyhz(ic,jc,kc)=zero;
      enddo
    enddo
    call myupdate_halo(uxhx, mb1, hi_uxhx)
    call myupdate_halo(uyhx, mb1, hi_uyhx)
    call myupdate_halo(uzhx, mb1, hi_uzhx)
    call myupdate_halo(uxhz, mb1, hi_uxhz)
    call myupdate_halo(uyhz, mb1, hi_uyhz)
    call myupdate_halo(uzhz, mb1, hi_uzhz)

    ! RhsX Part ===========
    s3tot=zero
    DO kc=ystart(3),yend(3)
      ks=kc-2;km=kc-1;kp=kc+1;ku=kc+2
      do jc=ystart(2),yend(2)
        jm=jc-1;jp=jc+1
        sucaj=rdyp(jc)
        InterpY1= YinterpCoe(jm); InterpY2=one-InterpY1  
        InterpY3= YinterpCoe(jc); InterpY4=one-InterpY3         
        do ic=ystart(1),yend(1)
          is=ic-2;im=ic-1;ip=ic+1;iu=ic+2
 
          uxhyc= InterpY1*ux(ic,jm,kc) +InterpY2*ux(ic,jc,kc)
          uxhyp= InterpY3*ux(ic,jc,kc) +InterpY4*ux(ic,jp,kc)
          h11=  dxCoe1*uxhx(is,jc,kc)*(uxhx(is,jc,kc)+uCRF) +dxCoe2*uxhx(im,jc,kc)*(uxhx(im,jc,kc)+uCRF) &
               +dxCoe3*uxhx(ic,jc,kc)*(uxhx(ic,jc,kc)+uCRF) +dxCoe4*uxhx(ip,jc,kc)*(uxhx(ip,jc,kc)+uCRF)
          h12=( uyhx(ic,jp,kc)*uxhyp -uyhx(ic,jc,kc)*uxhyc )*sucaj
          h13=  dzCoe1*uxhz(ic,jc,km)*uzhx(ic,jc,km) +dzCoe2*uxhz(ic,jc,kc)*uzhx(ic,jc,kc) &
               +dzCoe3*uxhz(ic,jc,kp)*uzhx(ic,jc,kp) +dzCoe4*uxhz(ic,jc,ku)*uzhx(ic,jc,ku)
          
          d11q1= dxxCoe1*ux(is,jc,kc) +dxxCoe2*ux(im,jc,kc) +dxxCoe3*ux(ic,jc,kc) +dxxCoe4*ux(ip,jc,kc) +dxxCoe5*ux(iu,jc,kc)
          d22q1= ap2c(jc)*ux(ic,jp,kc)+ac2c(jc)*ux(ic,jc,kc)+am2c(jc)*ux(ic,jm,kc)                
          d33q1= dzzCoe1*ux(ic,jc,ks) +dzzCoe2*ux(ic,jc,km) +dzzCoe3*ux(ic,jc,kc) +dzzCoe4*ux(ic,jc,kp) +dzzCoe5*ux(ic,jc,ku)
          dcq13= d11q1+d33q1
          s3tot=s3tot+(dcq13+d22q1)*dyp(jc)

          convEd1= -h11-h12-h13 +xnu*dcq13 + gravity(1)
          gradp1= dxCoe1*pressure(is,jc,kc)+dxCoe2*pressure(im,jc,kc) +dxCoe3*pressure(ic,jc,kc) +dxCoe4*pressure(ip,jc,kc)
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1+ two*pmBeta*d22q1
          HistXOld(ic,jc,kc)=convEd1
        enddo
      enddo
    ENDDO

    ! RhsY Part ===========
    DO kc=ystart(3),yend(3)
      ks=kc-2;km=kc-1;kp=kc+1;ku=kc+2
      do jc=ystart(2),yend(2)
        jm=jc-1;jp=jc+1
        sucac= rdyc(jc)
        qsucac=quarter*sucac
        InterpY1= YinterpCoe(jm); InterpY2=one-InterpY1
        do ic=ystart(1),yend(1)
          is=ic-2;im=ic-1;ip=ic+1;iu=ic+2

          uxhym= InterpY1*ux(im,jm,kc)+InterpY2*ux(im,jc,kc)+uCRF;  uzhym=InterpY1*uz(ic,jm,km)+InterpY2*uz(ic,jc,km)
          uxhyc= InterpY1*ux(ic,jm,kc)+InterpY2*ux(ic,jc,kc)+uCRF;  uzhyc=InterpY1*uz(ic,jm,kc)+InterpY2*uz(ic,jc,kc)
          uxhyp= InterpY1*ux(ip,jm,kc)+InterpY2*ux(ip,jc,kc)+uCRF;  uzhyp=InterpY1*uz(ic,jm,kp)+InterpY2*uz(ic,jc,kp)
          uxhyu= InterpY1*ux(iu,jm,kc)+InterpY2*ux(iu,jc,kc)+uCRF;  uzhyu=InterpY1*uz(ic,jm,ku)+InterpY2*uz(ic,jc,ku)
          h21=  dxCoe1*uyhx(im,jc,kc)*uxhym +dxCoe2*uyhx(ic,jc,kc)*uxhyc +dxCoe3*uyhx(ip,jc,kc)*uxhyp +dxCoe4*uyhx(iu,jc,kc)*uxhyu
          h22=( (uy(ic,jp,kc)+uy(ic,jc,kc))* (uy(ic,jp,kc)+uy(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jm,kc)) )*qsucac
          h23=  dzCoe1*uyhz(ic,jc,km)*uzhym +dzCoe2*uyhz(ic,jc,kc)*uzhyc +dzCoe3*uyhz(ic,jc,kp)*uzhyp +dzCoe4*uyhz(ic,jc,ku)*uzhyu
          
          d11q2= dxxCoe1*uy(is,jc,kc) +dxxCoe2*uy(im,jc,kc) +dxxCoe3*uy(ic,jc,kc) +dxxCoe4*uy(ip,jc,kc) +dxxCoe5*uy(iu,jc,kc)
          d22q2= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          d33q2= dzzCoe1*uy(ic,jc,ks) +dzzCoe2*uy(ic,jc,km) +dzzCoe3*uy(ic,jc,kc) +dzzCoe4*uy(ic,jc,kp) +dzzCoe5*uy(ic,jc,ku)
          dcq13= d11q2+d33q2
  
          convEd2= -h21-h22-h23 +xnu*dcq13+ gravity(2)     
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2+ two*pmBeta*d22q2
          HistYOld(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO

    ! RhsZ Part ===========
    DO kc=ystart(3),yend(3)
      ks=kc-2;km=kc-1;kp=kc+1;ku=kc+2
      do jc=ystart(2),yend(2)
        jm=jc-1;jp=jc+1
        sucaj=rdyp(jc)
        InterpY1= YinterpCoe(jm); InterpY2=one-InterpY1  
        InterpY3= YinterpCoe(jc); InterpY4=one-InterpY3
        do ic=ystart(1),yend(1)
          is=ic-2;im=ic-1;ip=ic+1;iu=ic+2

          uzhyc= InterpY1*uz(ic,jm,kc) +InterpY2*uz(ic,jc,kc)
          uzhyp= InterpY3*uz(ic,jc,kc) +InterpY4*uz(ic,jp,kc)
          h31=  dxCoe1*uzhx(im,jc,kc)*(uxhz(im,jc,kc)+uCRF) +dxCoe2*uzhx(ic,jc,kc)*(uxhz(ic,jc,kc)+uCRF) &
               +dxCoe3*uzhx(ip,jc,kc)*(uxhz(ip,jc,kc)+uCRF) +dxCoe4*uzhx(iu,jc,kc)*(uxhz(iu,jc,kc)+uCRF)
          h32=( uyhz(ic,jp,kc)*uzhyp-uyhz(ic,jc,kc)*uzhyc)*sucaj
          h33=  dzCoe1*uzhz(ic,jc,ks)*uzhz(ic,jc,ks) +dzCoe2*uzhz(ic,jc,km)*uzhz(ic,jc,km) &
               +dzCoe3*uzhz(ic,jc,kc)*uzhz(ic,jc,kc) +dzCoe4*uzhz(ic,jc,kp)*uzhz(ic,jc,kp)

          d11q3= dxxCoe1*uz(is,jc,kc) +dxxCoe2*uz(im,jc,kc) +dxxCoe3*uz(ic,jc,kc) +dxxCoe4*uz(ip,jc,kc) +dxxCoe5*uz(iu,jc,kc)
          d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
          d33q3= dzzCoe1*uz(ic,jc,ks) +dzzCoe2*uz(ic,jc,km) +dzzCoe3*uz(ic,jc,kc) +dzzCoe4*uz(ic,jc,kp) +dzzCoe5*uz(ic,jc,ku)
          dcq13= d11q3+d33q3
                   
          convEd3= -h31-h32-h33 +xnu*dcq13+ gravity(3)
          gradp3= dzCoe1*pressure(ic,jc,ks) +dzCoe2*pressure(ic,jc,km) +dzCoe3*pressure(ic,jc,kc) +dzCoe4*pressure(ic,jc,kp)
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3+ two*pmBeta*d22q3
          HistZOld(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO

    ! in dp1ns there is the mean pressure gradient to keep constant mass
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
  end subroutine clcRhsXYZ

  !******************************************************************
  ! clcU1Hat
  !******************************************************************    
  subroutine clcU1Hat(ux,RhsX)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(ystart(1):yend(1), ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj

    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc) 
        tridcj(ic,jc) = -pmBeta*ac2cForCN(jc)+one
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
      enddo
    enddo      
    DO kc=ystart(3),yend(3)  
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1) 
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
  end subroutine clcU1Hat
  
  !******************************************************************
  ! clcU2Hat
  !******************************************************************  
  subroutine clcU2Hat(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(in):: duy_ym
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(ystart(1):yend(1), ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj

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
  end subroutine clcU2Hat
  
  !******************************************************************
  ! clcU3Hat
  !******************************************************************  
  subroutine clcU3Hat(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
  
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2))::tridmj,tridcj,tridpj,tridfj

    do jc=ystart(2),yend(2)
      do ic=ystart(1),yend(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc) 
        tridcj(ic,jc) = -pmBeta*ac2cForCN(jc)+one
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
      enddo
    enddo    
    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          tridfj(ic,jc) = RhsZ(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,ysize(1),nyc)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat

  !******************************************************************
  ! PressureUpdate
  !******************************************************************
  subroutine PressureUpdate(pressure, prphiHalo)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out):: pressure
    
    integer::ic,jc,kc,jp,jm
    real(RK)::pmBetap,pmBetac,pmBetam

    DO kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        jp=jc+1;  jm=jc-1
        pmBetap= -pmBeta*ap2ph(jc)
        pmBetac= -pmBeta*ac2ph(jc) + one
        pmBetam= -pmBeta*am2ph(jc)
        do ic=ystart(1),yend(1)
          pressure(ic,jc,kc)= pressure(ic,jc,kc)+ pmBetap*prphiHalo(ic,jp,kc)+ pmBetac*prphiHalo(ic,jc,kc)+ &
                                                  pmBetam*prphiHalo(ic,jm,kc)
        enddo
      enddo
    ENDDO
  end subroutine PressureUpdate

end module m_TScheme
