!********************************************************************!
!*    file name  : mydecomp_2d_extra.f90                            *!
!*    module name: mydecomp_2d_extra                                *!  
!*                                                                  *!
!*    purpose:                                                      *! 
!*      1)  update halo                                             *!
!*      2)                                                          *!
!*                                                                  *!
!*  Author: Zheng Gong           Date: 23:Feb:2020                  *!
!*                                                                  *!
!********************************************************************!
module mydecomp_2d_extra
  use MPI
  use decomp_2d
  implicit none
  private

  ! Define neighboring blocks (including the Bc info)
  !   first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
  !   second dimension 4 neighbour processors
  integer,save,public,dimension(3,4):: myProcNghBC

  type,public:: MatBound
    integer:: pencil
    integer:: xmr, xmm, xme  ! real index, memory index, expand size in xm_dir respectively, xmm=xmr-xme
    integer:: xpr, xpm, xpe  ! real index, memory index, expand size in xp_dir respectively, xpm=xpr+xpe 
    integer:: ymr, ymm, yme  ! real index, memory index, expand size in ym_dir respectively, ymm=ymr-yme
    integer:: ypr, ypm, ype  ! real index, memory index, expand size in yp_dir respectively, ypm=ypr+ype
    integer:: zmr, zmm, zme  ! real index, memory index, expand size in zm_dir respectively, zmm=zmr-zme
    integer:: zpr, zpm, zpe  ! real index, memory index, expand size in zp_dir respectively, zpm=zpr+zpe
  end type MatBound

  type,public:: HaloInfo
    integer :: pencil
    integer :: xmh, xph   ! xm_dir halo, xp_dir halo
    integer :: ymh, yph   ! ym_dir halo, yp_dir halo
    integer :: zmh, zph   ! zm_dir halo, zp_dir halo
  end type HaloInfo

  public:: myInit_neighbour, myallocate, myupdate_halo, my_best_2d_grid

contains

  !**********************************************************************
  ! myInit_neighbour
  !**********************************************************************
  subroutine myInit_neighbour(row,col,BcOpt)
    implicit none
    integer,intent(in)::row,col,BcOpt(6)

    !locals
    integer :: i,j,coord1,coord2,idTop,idBottom,idLeft,idRight,ierror

    ! ------------------------- neigbor information begins------------------------
    ! 
    ! 2D domain decomposition method in Decomp2D(from left to right:  x-pencil, y-pencil, z-pencil)
    ! 
    !   If we have 6 processors = 3 row * 2 col
    ! 
    !     the arrangement of the smubomains(nrank) is as follow:
    !       y               x               x
    !       |        4 5    |        4 5    |        4 5   
    !       |        2 3    |        2 3    |        2 3
    !       |_ _ _z  0 1    |_ _ _z  0 1    |_ _ _y  0 1
    !
    !     the arrangement of the coord1 is as follow:
    !       y               x               x
    !       |        2 2    |        2 2    |        2 2
    !       |        1 1    |        1 1    |        1 1
    !       |_ _ _z  0 0    |_ _ _z  0 0    |_ _ _y  0 0
    ! 
    !     the arrangement of the coord2 is as follow:
    !       y               x               x
    !       |        0 1    |        0 1    |        0 1
    !       |        0 1    |        0 1    |        0 1
    !       |_ _ _ z 0 1    |_ _ _z  0 1    |_ _ _y  0 1
    ! 
    !     neighbor index:
    ! 
    !       y               x               x
    !       |        6 3 5  |        6 3 5  |        6 3 5
    !       |        2 0 1  |        2 0 1  |        2 0 1
    !       |_ _ _ z 7 4 8  |_ _ _z  7 4 8  |_ _ _y  7 4 8
    ! 
    !        Here 0 means the center smubomain, and 1-8 stands for the reduative location of the eight neighbors

    coord1 = int ( nrank / col)
    coord2 = mod ( nrank,  col)

    ! Firstly, all the boundaries are assumed to be periodic
    idTop     = mod(coord1+1,    row)  ! top
    idBottom  = mod(coord1+row-1,row)  ! bottom 
    idLeft    = mod(coord2+col-1,col)  ! left
    idRight   = mod(coord2+1,    col)  ! right
    do i=1,3
      myProcNghBC(i,1) = coord1   * col + idRight; ! myProcNghBC(i,5) = idTop    * col + idRight
      myProcNghBC(i,2) = coord1   * col + idLeft ; ! myProcNghBC(i,6) = idTop    * col + idLeft
      myProcNghBC(i,3) = idTop    * col + coord2 ; ! myProcNghBC(i,7) = idBottom * col + idLeft
      myProcNghBC(i,4) = idBottom * col + coord2 ; ! myProcNghBC(i,8) = idBottom * col + idRight
    enddo

    ! Secondly, modify the edge neighbour ids
    IF(coord1==0) THEN
      if(BcOpt(3)<0) myProcNghBC(1,4)=BcOpt(3)
      if(BcOpt(1)<0) myProcNghBC(2,4)=BcOpt(1)
      if(BcOpt(1)<0) myProcNghBC(3,4)=BcOpt(1)
    ENDIF
    IF(coord1==row-1) THEN
      if(BcOpt(4)<0) myProcNghBC(1,3)=BcOpt(4)
      if(BcOpt(2)<0) myProcNghBC(2,3)=BcOpt(2)
      if(BcOpt(2)<0) myProcNghBC(3,3)=BcOpt(2)
    ENDIF
    IF(coord2==0) THEN
      if(BcOpt(5)<0) myProcNghBC(1,2)=BcOpt(5)
      if(BcOpt(5)<0) myProcNghBC(2,2)=BcOpt(5)
      if(BcOpt(3)<0) myProcNghBC(3,2)=BcOpt(3)
    ENDIF
    IF(coord2==col-1) THEN
      if(BcOpt(6)<0) myProcNghBC(1,1)=BcOpt(6)
      if(BcOpt(6)<0) myProcNghBC(2,1)=BcOpt(6)
      if(BcOpt(4)<0) myProcNghBC(3,1)=BcOpt(4)
    ENDIF

  end subroutine myInit_neighbour

  !**********************************************************************
  ! myallocate
  !**********************************************************************
  subroutine myallocate(var, mb, opt_decomp, opt_global)
    implicit none
    real(mytype), allocatable, dimension(:,:,:),intent(OUT) :: var
    type(MatBound),intent(INOUT)::mb
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
 
    ! locals
    TYPE(DECOMP_INFO) :: decomp
    logical :: global
    integer :: alloc_stat, errorcode
  
    if(mb%pencil<1 .or. mb%pencil>3) then
      errorcode = 8
      call decomp_2d_abort(errorcode,'wrong input pencil when creating new arrays')
    endif

    if (present(opt_decomp)) then
      decomp = opt_decomp
    else
      decomp = decomp_main
    endif
    if (present(opt_global)) then
      global = opt_global
    else
      global = .false.
    endif

    ! first update the MatBound
    select case(mb%pencil)
    case(1) ! x-pencil
      if (global) then
        mb%xmr = decomp%xst(1); mb%xpr = decomp%xen(1)
        mb%ymr = decomp%xst(2); mb%ypr = decomp%xen(2)
        mb%zmr = decomp%xst(3); mb%zpr = decomp%xen(3)
      else
        mb%xmr = 1; mb%xpr = decomp%xsz(1)
        mb%ymr = 1; mb%ypr = decomp%xsz(2)
        mb%zmr = 1; mb%zpr = decomp%xsz(3)
      end if
    case(2)  ! y-pencil
      if (global) then
        mb%xmr = decomp%yst(1); mb%xpr = decomp%yen(1)
        mb%ymr = decomp%yst(2); mb%ypr = decomp%yen(2)
        mb%zmr = decomp%yst(3); mb%zpr = decomp%yen(3)
      else
        mb%xmr = 1; mb%xpr = decomp%ysz(1)
        mb%ymr = 1; mb%ypr = decomp%ysz(2)
        mb%zmr = 1; mb%zpr = decomp%ysz(3)
      end if
    case(3)  ! z-pencil
      if (global) then
        mb%xmr = decomp%zst(1); mb%xpr = decomp%zen(1)
        mb%ymr = decomp%zst(2); mb%ypr = decomp%zen(2)
        mb%zmr = decomp%zst(3); mb%zpr = decomp%zen(3)
      else
        mb%xmr = 1; mb%xpr = decomp%zsz(1)
        mb%ymr = 1; mb%ypr = decomp%zsz(2)
        mb%zmr = 1; mb%zpr = decomp%zsz(3)
      end if
    end select
    mb%xmm= mb%xmr- mb%xme;  mb%xpm= mb%xpr+ mb%xpe
    mb%ymm= mb%ymr- mb%yme;  mb%ypm= mb%ypr+ mb%ype
    mb%zmm= mb%zmr- mb%zme;  mb%zpm= mb%zpr+ mb%zpe 


    allocate(var(mb%xmm:mb%xpm,mb%ymm:mb%ypm,mb%zmm:mb%zpm),stat=alloc_stat)
    if (alloc_stat /= 0) then
      errorcode = 8
      call decomp_2d_abort(errorcode,'Memory allocation failed when creating new arrays')
    end if
  end subroutine myallocate

  !**********************************************************************
  ! myupdate_halo
  !**********************************************************************
  subroutine myupdate_halo(mat,mb,hi)
    implicit none
    type(MatBound),intent(in)::mb
    type(HaloInfo),intent(in)::hi
    real(mytype), dimension(mb%xmm:mb%xpm,mb%ymm:mb%ypm,mb%zmm:mb%zpm),intent(INOUT):: mat

    ! locals
    integer:: j,data_type,halo_type,ierror,ProcNgh(4) 
    integer:: s1,s2,s3,icount,ilength,ijump,requests(2)
    integer:: xs1,xs2,ys1,ys2,zs1,zs2   ! index for halo sending
    integer:: xr1,xr2,yr1,yr2,zr1,zr2   ! index for halo receiving
    integer,dimension(MPI_STATUS_SIZE,2) :: SRstatus
  
    data_type=real_type
    if(mb%pencil /= hi%pencil ) then
      ierror = 14
      call decomp_2d_abort(ierror,'myupdate_halo wrong: pencil error, Gong Zheng')
    endif
    if(mb%xme<hi%xmh .or. mb%yme<hi%ymh .or. mb%zme<hi%zmh .or.  &
       mb%xpe<hi%xph .or. mb%ype<hi%yph .or. mb%zpe<hi%zph .or.  &
       hi%xmh<0      .or. hi%ymh<0      .or. hi%zmh<0      .or.  &
       hi%xph<0      .or. hi%yph<0      .or. hi%zph<0       ) then
      ierror = 14
      call decomp_2d_abort(ierror,'myupdate_halo wrong: HaloSize error, Gong Zheng')
    endif

    s1= mb%xpm-mb%xmm+1
    s2= mb%ypm-mb%ymm+1
    s3= mb%zpm-mb%zmm+1
    do j=1,4
      if(myProcNghBC(mb%pencil,j)<0) then
        ProcNgh(j)=MPI_PROC_NULL
      else
	    ProcNgh(j)=myProcNghBC(mb%pencil,j)
      endif
    enddo

    SELECT CASE(mb%pencil)
    CASE(1)  ! x-pencil

      ! step1: receive from xm_dir, and send to xp_dir
      IF(hi%xmh >0 ) THEN
        xr1=mb%xmr-hi%xmh;      xr2=mb%xmr-1
        yr1=mb%ymr;             yr2=mb%ypr  
        zr1=mb%zmr;             zr2=mb%zpr

        xs1=mb%xpr-hi%xmh+1;    xs2=mb%xpr
        ys1=yr1;                ys2=yr2   
        zs1=zr1;                zs2=zr2
        mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
      ENDIF

      ! step2: receive from xp_dir, and send to xm_dir
      IF(hi%xph >0 ) THEN
        xr1=mb%xpr+1;           xr2=mb%xpr+hi%xph
        yr1=mb%ymr;             yr2=mb%ypr  
        zr1=mb%zmr;             zr2=mb%zpr

        xs1=mb%xmr;             xs2=mb%xmr+hi%xph-1
        ys1=yr1;                ys2=yr2   
        zs1=zr1;                zs2=zr2
        mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
      ENDIF



    CASE(2)  ! y-pencil

      ! step1: receive from ym_dir, and send to yp_dir
      IF(hi%ymh >0 ) THEN
        xr1=mb%xmr;             xr2=mb%xpr 
        yr1=mb%ymr-hi%ymh;      yr2=mb%ymr-1  
        zr1=mb%zmr;             zr2=mb%zpr

        xs1=xr1;                xs2=xr2
        ys1=mb%ypr-hi%ymh+1;    ys2=mb%ypr    
        zs1=zr1;                zs2=zr2
        mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
      ENDIF

      ! step2: receive from yp_dir, and send to ym_dir
      IF(hi%yph >0 ) THEN
        xr1=mb%xmr;             xr2=mb%xpr 
        yr1=mb%ypr+1;           yr2=mb%ypr+hi%yph  
        zr1=mb%zmr;             zr2=mb%zpr

        xs1=xr1;                xs2=xr2
        ys1=mb%ymr;             ys2=mb%ymr+hi%yph-1   
        zs1=zr1;                zs2=zr2
        mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
      ENDIF

      ! step3: receive from xm_dir, and send to xp_dir
      IF(hi%xmh > 0) THEN
        xr1=mb%xmr-hi%xmh;      xr2=mb%xmr-1
        yr1=mb%ymm;             yr2=mb%ypm
        zr1=mb%zmm;             zr2=mb%zpm
        
        xs1=mb%xpr-hi%xmh+1;    xs2=mb%xpr
        ys1=yr1;                ys2=yr2    
        zs1=zr1;                zs2=zr2
        if(ProcNgh(3)==nrank) then  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
          mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
        else
          icount = s2*s3
          ilength= hi%xmh
          ijump  = s1
          call MPI_TYPE_VECTOR(icount,ilength,ijump, data_type, halo_type, ierror)
          call MPI_TYPE_COMMIT(halo_type, ierror)
          call MPI_IRECV( mat(xr1,yr1,zr1),1,halo_type,ProcNgh(4),1,MPI_COMM_WORLD,requests(1),ierror)
          call MPI_ISSEND(mat(xs1,ys1,zs1),1,halo_type,ProcNgh(3),1,MPI_COMM_WORLD,requests(2),ierror)
          call MPI_WAITALL(2,requests,SRstatus,ierror)
          call MPI_TYPE_FREE(halo_type, ierror)
        endif
      ENDIF

      ! step4: receive from xp_dir, and send to xm_dir
      IF(hi%xph > 0) THEN
        xr1=mb%xpr+1;           xr2=mb%xpr+hi%xph
        yr1=mb%ymm;             yr2=mb%ypm
        zr1=mb%zmm;             zr2=mb%zpm
        
        xs1=mb%xmr;             xs2=mb%xmr+hi%xph-1
        ys1=yr1;                ys2=yr2    
        zs1=zr1;                zs2=zr2
        if(ProcNgh(4)==nrank) then  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
          mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
        else
          icount = s2*s3
          ilength= hi%xph
          ijump  = s1
          call MPI_TYPE_VECTOR(icount, ilength, ijump, data_type, halo_type, ierror)
          call MPI_TYPE_COMMIT(halo_type, ierror)
          call MPI_IRECV( mat(xr1,yr1,zr1),1,halo_type,ProcNgh(3),2,MPI_COMM_WORLD,requests(1),ierror)
          call MPI_ISSEND(mat(xs1,ys1,zs1),1,halo_type,ProcNgh(4),2,MPI_COMM_WORLD,requests(2),ierror)
          call MPI_WAITALL(2,requests,SRstatus,ierror)
          call MPI_TYPE_FREE(halo_type, ierror)
        endif
      ENDIF

      ! step5: receive from zm_dir, and send to zp_dir
      IF(hi%zmh > 0) THEN
        xr1=mb%xmm;             xr2=mb%xpm
        yr1=mb%ymm;             yr2=mb%ypm
        zr1=mb%zmr-hi%zmh;      zr2=mb%zmr-1
        
        xs1=xr1;                xs2=xr2
        ys1=yr1;                ys2=yr2    
        zs1=mb%zpr-hi%zmh+1;    zs2=mb%zpr
        if(ProcNgh(1)==nrank) then  ! neighbour is nrank itself, and ProcNgh(2)==nrank also !!!
          mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
        else
          icount=s1*s2*hi%zmh
          call MPI_IRECV( mat(xr1,yr1,zr1),icount,data_type,ProcNgh(2),3,MPI_COMM_WORLD,requests(1),ierror)
          call MPI_ISSEND(mat(xs1,ys1,zs1),icount,data_type,ProcNgh(1),3,MPI_COMM_WORLD,requests(2),ierror)
          call MPI_WAITALL(2,requests,SRstatus,ierror)
        endif
      ENDIF

      ! step6: receive from zp_dir, and send to zm_dir
      IF(hi%zph > 0) THEN
        xr1=mb%xmm;             xr2=mb%xpm
        yr1=mb%ymm;             yr2=mb%ypm
        zr1=mb%zpr+1;           zr2=mb%zpr+hi%zph
        
        xs1=xr1;                xs2=xr2
        ys1=yr1;                ys2=yr2    
        zs1=mb%zmr;             zs2=mb%zmr+hi%zph-1
        if(ProcNgh(2)==nrank) then  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
          mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
        else
          icount=s1*s2*hi%zph
          call MPI_IRECV( mat(xr1,yr1,zr1),icount,data_type,ProcNgh(1),4,MPI_COMM_WORLD,requests(1),ierror)
          call MPI_ISSEND(mat(xs1,ys1,zs1),icount,data_type,ProcNgh(2),4,MPI_COMM_WORLD,requests(2),ierror)
          call MPI_WAITALL(2,requests,SRstatus,ierror)
        endif
      ENDIF

    CASE(3)  ! z-pencil





    END SELECT

  end subroutine myupdate_halo

#include "factor.inc"
  !**********************************************************************
  ! my_best_2d_grid
  !**********************************************************************
  subroutine my_best_2d_grid(nx,ny,nz,pencil,lvlhalo_test,time_halo,BcOpt,iproc,best_p_row,best_p_col)
    implicit none
    integer,intent(in)::  nx,ny,nz,pencil,lvlhalo_test,time_halo,iproc,BcOpt(6)
    integer,intent(out):: best_p_row, best_p_col

    ! locals
    TYPE(MatBound):: mb_test
    TYPE(HaloInfo):: hi_test
    real(mytype):: t1, t2, best_time
    integer,allocatable, dimension(:) :: factors
    integer:: nfact,i,j,k,row,col,ierror,errorcode
    real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3

    if (nrank==0) write(*,*) 'My auto-tuning mode......'

    best_time = huge(t1)
    best_p_row = -1
    best_p_col = -1

    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if (nrank==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

    mb_test%pencil = pencil  
    mb_test%xme=lvlhalo_test;  mb_test%xpe=lvlhalo_test
    mb_test%yme=lvlhalo_test;  mb_test%ype=lvlhalo_test
    mb_test%zme=lvlhalo_test;  mb_test%zpe=lvlhalo_test

    hi_test%pencil = pencil
    hi_test%xmh=lvlhalo_test;  hi_test%xph=lvlhalo_test
    hi_test%ymh=lvlhalo_test;  hi_test%yph=lvlhalo_test
    hi_test%zmh=lvlhalo_test;  hi_test%zph=lvlhalo_test

    do i=1, nfact

       row = factors(i)
       col = iproc / row

       ! enforce the limitation of 2D decomposition
       if (min(nx,ny)>=row .and. min(ny,nz)>=col) then

           call decomp_2d_init(nx,ny,nz,row,col)
           call myInit_neighbour(row,col,BcOpt)

          ! arrays for X,Y and Z-pencils
          allocate(u1(xsize(1),xsize(2),xsize(3)))
          allocate(u2(ysize(1),ysize(2),ysize(3)))
          allocate(u3(zsize(1),zsize(2),zsize(3)))

          ! timing the transposition routines
          t1 = MPI_WTIME()
          do k=1,5
            call transpose_y_to_x(u2,u1)
            call transpose_x_to_z(u1,u3)
            call transpose_z_to_y(u3,u2)
            call transpose_y_to_z(u2,u3)        
            call transpose_z_to_x(u3,u1)
            call transpose_x_to_y(u1,u2)
          enddo
          t2 = MPI_WTIME() - t1

          deallocate(u1,u2,u3)
          call myallocate(u1, mb_test)

          t1 = MPI_WTIME()
          do k=1,5
            do j=1,time_halo
              call myupdate_halo(u1,mb_test,hi_test)
            enddo
          enddo
          t2 = MPI_WTIME() - t1 + t2

          deallocate(u1)
          call decomp_2d_finalize

          call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
          t1 = t1 / real(iproc,kind=mytype)

          if (nrank==0) then
             write(*,*) 'processor grid', row, ' by ', col, ' time=', t1
          endif

          if (best_time > t1) then
             best_time = t1
             best_p_row = row
             best_p_col = col
          endif

       endif

    enddo ! loop through processer grid
    deallocate(factors)

    if (best_p_row/=-1) then
       if (nrank==0) then
          write(*,*) 'the best processor grid is probably ', &
               best_p_row, ' by ', best_p_col
       endif
    else
       errorcode = 9
       call decomp_2d_abort(errorcode, &
            'The processor-grid auto-tuning code failed. ' // &
            'The number of processes requested is probably too large.')
    endif
  end subroutine my_best_2d_grid

end module mydecomp_2d_extra
