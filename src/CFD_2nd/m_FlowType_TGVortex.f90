module m_FlowType_TGVortex
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_MeshAndMetries
  use decomp_2d
  use m_Variables,only: mb1
  use m_Tools,only:CalcDissipationRate
  implicit none
  private

  public:: InitVelocity_TG, Update_uy_ym_TG
  public:: InitStatVar_TG,  clcStat_TG
contains

  !******************************************************************
  ! InitVelocity_TG
  !******************************************************************
  subroutine InitVelocity_TG(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)),intent(inout)::Deviation
  
    ! locals
    integer :: ic,jc,kc
	  real(RK):: VelRef,LenRef,xpt,ypt,zpt,xct,yct,zct
      
    VelRef=one
    LenRef=one
    do kc=ystart(3),yend(3)
      zpt=real(kc-1,kind=RK)*dz
      zct=zpt+half*dz
      do jc=ystart(2),yend(2)
        ypt=yp(jc)
        yct=yc(jc)
        do ic=ystart(1),yend(1)
          xpt=real(ic-1,kind=RK)*dx
          xct=xpt+half*dx
          ux(ic,jc,kc) =  VelRef*sin(xpt/LenRef)*cos(yct/LenRef)*cos(zct/LenRef)              
          uy(ic,jc,kc) = -VelRef*cos(xct/LenRef)*sin(ypt/LenRef)*cos(zct/LenRef)
          uz(ic,jc,kc) =  zero !VelRef*cos(xct/LenRef)*cos(yct/LenRef)*sin(zpt/LenRef)
        enddo
      enddo
    enddo
    
  end subroutine InitVelocity_TG

  !******************************************************************
  ! Update_uy_ym_TG
  !******************************************************************   
  subroutine Update_uy_ym_TG(uy_ym, duy_ym, TimeNew)
    implicit none
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(out):: uy_ym
    real(RK),dimension(ystart(1):yend(1),ystart(3):yend(3)),intent(inout):: duy_ym    
    real(RK),intent(in):: TimeNew
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = zero
     
    duy_ym = uy_ym - duy_ym
    
  end subroutine Update_uy_ym_TG

  !******************************************************************
  ! InitStatVar_TG
  !******************************************************************
  subroutine InitStatVar_TG()
    implicit none

    ! locals
    integer:: myistat
    character(len=50)::filename

    if(nrank/=0) return
    filename=trim(Res_Dir)//'TG_dissp.txt'
    open(69, file=filename,status='replace',form='formatted',IOSTAT=myistat)
    if(myistat /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_TG: ","Cannot open file: "//trim(filename))
    close(69,IOSTAT=myistat)

  end subroutine InitStatVar_TG

  !******************************************************************
  ! clcStat_TG
  !******************************************************************
  subroutine clcStat_TG(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
   
    ! locals
    integer::ic,jc,kc,ierror
    real(Rk):: sum_dissp,sumr
    character(len=50)::filename
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))::dissp
 
    sumr= zero
    call CalcDissipationRate(ux,uy,uz,dissp)
    do kc=ystart(3),yend(3)
      do jc=ystart(2),yend(2)
        do ic=ystart(1),yend(1)
          sumr=sumr+ dissp(ic,jc,kc)
        enddo
      enddo
    enddo
    call MPI_REDUCE(sumr,sum_dissp,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)

    if(nrank==0) then
      filename=trim(Res_Dir)//'TG_dissp.txt'
      open(unit=69, file=filename, status='old',position='append',form='formatted',IOSTAT=ierror )
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_TG: ","Cannot open file")
      ELSE
        write(69,'(2ES24.15)')SimTime,xnu*sum_dissp/real(nxc*nyc*nzc,RK)
      ENDIF
      close(69,IOSTAT=ierror)
    endif
  end subroutine clcStat_TG
    
end module m_FlowType_TGVortex
