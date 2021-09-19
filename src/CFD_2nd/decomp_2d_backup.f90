module decomp_2d
  use MPI
  implicit none
  private

#ifdef DOUBLE_PREC
  integer,parameter,public:: mytype = KIND(0.0D0)
  integer,parameter,public:: real_type = MPI_DOUBLE_PRECISION
#ifdef SAVE_SINGLE
  integer,parameter,public:: mytype_save = KIND(0.0)
  integer,parameter,public:: real_type_save = MPI_REAL
#else
  integer,parameter,public:: mytype_save = KIND(0.0D0)
  integer,parameter,public:: real_type_save = MPI_DOUBLE_PRECISION
#endif
#else
  integer,parameter,public:: mytype = KIND(0.0)
  integer,parameter,public:: real_type = MPI_REAL
  integer,parameter,public:: mytype_single = KIND(0.0)
  integer,parameter,public:: real_type_single = MPI_REAL
#endif
  integer,save,public:: mytype_bytes

  ! some key global variables
  integer,save,public:: nrank  ! local MPI rank 
  integer,save,public:: nproc  ! total number of processors

  ! parameters for 2D Cartesian topology
  integer,save::np_row,np_col,DECOMP_2D_COMM_ROW,DECOMP_2D_COMM_COL

  ! staring/ending index and size of data held by current processor
  ! duplicate 'decomp_main', needed by apps to define data structure 
  integer,save,dimension(3),public::xstart,xend,xsize  ! x-pencil
  integer,save,dimension(3),public::ystart,yend,ysize  ! y-pencil
  integer,save,dimension(3),public::zstart,zend,zsize  ! z-pencil

  ! These are the buffers used by MPI_ALLTOALL(V) calls
  integer,save :: decomp_buf_size = 0
  real(mytype),allocatable,dimension(:)::work1_r,work2_r

  ! derived type to store decomposition info for a given global data size
  TYPE,public :: DECOMP_INFO
    ! staring/ending index and size of data held by current processor
    integer,dimension(3)::xst,xen,xsz  ! x-pencil
    integer,dimension(3)::yst,yen,ysz  ! y-pencil
    integer,dimension(3)::zst,zen,zsz  ! z-pencil

    ! how each dimension is distributed along pencils
    integer,allocatable, dimension(:)::x1size,y1size,y2size,z2size  

#ifdef EVEN
    logical :: even
    integer:: xycount,yzcount
#else
    ! send/receive buffer counts and displacements for MPI_ALLTOALLV
    integer, allocatable, dimension(:)::x1cnts, y1cnts, y2cnts, z2cnts
    integer, allocatable, dimension(:)::x1disp, y1disp, y2disp, z2disp
#endif   
                                                          
    ! This is only for the real datatype                                                    
    integer,dimension(:),allocatable::zcnts_xz,ztypes_xz,xcnts_xz,xtypes_xz                           
#ifdef MPI3
    integer :: xtozNeighborComm,ztoxNeighborComm
    integer,dimension(:),allocatable::xranks,zranks
    integer(MPI_ADDRESS_KIND),dimension(:),allocatable::zdispls_xz,xdispls_xz
#else
    integer,dimension(:),allocatable::zdispls_xz,xdispls_xz
#endif     
  END TYPE DECOMP_INFO
  type(decomp_info),save::decomp_main ! main (default) decomposition information

  ! public user routines
  public :: decomp_2d_init, decomp_2d_finalize,decomp_info_init, &
       decomp_info_finalize,decomp_2d_abort, &
       transpose_x_to_y, transpose_y_to_z, &
       transpose_z_to_y, transpose_y_to_x, &
       transpose_x_to_z, transpose_z_to_x, &
       decomp_2d_write_var,decomp_2d_read_var,decomp_2d_write_plane,decomp_2d_write_every

  ! Define neighboring blocks (including the Bc info)
  !   first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
  !   second dimension 4 neighbour processors
  integer,save,public,dimension(3,4):: myProcNghBC

  type,public:: MatBound
    integer:: pencil
    integer:: xmr,xmm, xme  ! real index, memory index, expand size in xm_dir respectively, xmm=xmr-xme
    integer:: xpr,xpm, xpe  ! real index, memory index, expand size in xp_dir respectively, xpm=xpr+xpe 
    integer:: ymr,ymm, yme  ! real index, memory index, expand size in ym_dir respectively, ymm=ymr-yme
    integer:: ypr,ypm, ype  ! real index, memory index, expand size in yp_dir respectively, ypm=ypr+ype
    integer:: zmr,zmm, zme  ! real index, memory index, expand size in zm_dir respectively, zmm=zmr-zme
    integer:: zpr,zpm, zpe  ! real index, memory index, expand size in zp_dir respectively, zpm=zpr+zpe
  end type MatBound

  type,public:: HaloInfo
    integer :: pencil
    integer :: xmh, xph   ! xm_dir halo, xp_dir halo
    integer :: ymh, yph   ! ym_dir halo, yp_dir halo
    integer :: zmh, zph   ! zm_dir halo, zp_dir halo
  end type HaloInfo
  public:: myInit_neighbour, myallocate, myupdate_halo, my_best_2d_grid
contains

  !******************************************************************
  ! Routine to be called by applications to initialise this library
  !******************************************************************
  subroutine decomp_2d_init(nx,ny,nz,p_row,p_col)
    implicit none
    integer,intent(in)::nx,ny,nz,p_row,p_col

    ! locals
    integer::errorcode,ierror,row,col

    if(p_row==0 .and. p_col==0)then
      call best_2d_grid(nx,ny,nz,nproc,row,col)
    else
      if(nproc/=p_row*p_col) then
        errorcode = 1
        call decomp_2d_abort(errorcode,'Invalid 2D processor grid - nproc /= p_row*p_col')
      else
        row = p_row
        col = p_col
      endif
    endif
    call decomp_info_init(nx,ny,nz,row,col,decomp_main)

    ! make a copy of the decomposition information 
    np_row = row
    np_col = col
    xstart = decomp_main%xst
    ystart = decomp_main%yst
    zstart = decomp_main%zst
    xend   = decomp_main%xen
    yend   = decomp_main%yen
    zend   = decomp_main%zen
    xsize  = decomp_main%xsz
    ysize  = decomp_main%ysz
    zsize  = decomp_main%zsz
    call MPI_TYPE_SIZE(real_type,mytype_bytes,ierror)
#ifdef EVEN
    if (nrank==0) print*, 'Padded ALLTOALL optimisation on'
#endif
  end subroutine decomp_2d_init

  !******************************************************************
  ! Routine to be called by applications to clean things up
  !******************************************************************
  subroutine decomp_2d_finalize()
    implicit none
    call decomp_info_finalize(decomp_main)
    decomp_buf_size=0
    deallocate(work1_r,work2_r)
  end subroutine decomp_2d_finalize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Advanced Interface allowing applications to define globle domain of
  ! any size, distribute it, and then transpose data among pencils.
  !  - generate 2D decomposition details as defined in DECOMP_INFO
  !  - the default global data size is nx*ny*nz
  !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
  !  - multiple global sizes can co-exist in one application, each
  !    using its own DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_init(nx,ny,nz,p_row,p_col,decomp)
    implicit none
    integer,intent(in)::nx,ny,nz,p_row,p_col
    type(decomp_info),intent(inout)::decomp

    ! locals
    integer::buf_size,ierror1,ierror2, errorcode,coord(2)
    integer,allocatable,dimension(:)::x1st,y1st,y2st,z2st
    integer,allocatable,dimension(:)::x1en,y1en,y2en,z2en
    integer::i,j,k,prank,subsize_y,offset_y,ierror,DECOMP_2D_COMM_CART
#ifdef MPI3
    integer::index_src, index_dest
    integer,dimension(nproc)::xranks,zranks,xweights,zweights
#endif

    ! verify the global size can actually be distributed as pencils
    if(nx<p_row .or. ny<p_row .or. ny<p_col .or. nz<p_col)then
      errorcode=6
      call decomp_2d_abort(errorcode,'Invalid 2D processor grid. ' // &
        'Make sure that min(nx,ny) >= p_row and min(ny,nz) >= p_col')
    endif

    call MPI_CART_CREATE(MPI_COMM_WORLD,2,(/p_row,p_col/),(/.false.,.false./),.false.,DECOMP_2D_COMM_CART, ierror)
    call MPI_CART_COORDS(DECOMP_2D_COMM_CART,nrank,2,coord,ierror)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART,(/.true.,.false./),DECOMP_2D_COMM_COL,ierror)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART,(/.false.,.true./),DECOMP_2D_COMM_ROW,ierror)

#ifdef EVEN
    if(mod(nx,p_row)==0 .and. mod(ny,p_row)==0 .and. mod(ny,p_col)==0 .and. mod(nz,p_col)==0) then
      decomp%even = .true.
    else
      decomp%even = .false.
    endif
#endif

    ! distribute mesh points
    allocate(decomp%x1size(0:p_row-1),decomp%y1size(0:p_row-1))
    allocate(decomp%y2size(0:p_col-1),decomp%z2size(0:p_col-1))
    
    !============  Add following lines by Gong Zheng, transpose_x_z purpose ==========    
    allocate(x1st(0:p_row-1),  x1en(0:p_row-1))
    allocate(y1st(0:p_row-1),  y1en(0:p_row-1))
    allocate(y2st(0:p_col-1),  y2en(0:p_col-1))
    allocate(z2st(0:p_col-1),  z2en(0:p_col-1)) 
    ! ================== adding ends ==================
    
    call distribute(nx,p_row,x1st,x1en); decomp%x1size=x1en-x1st+1
    call distribute(ny,p_row,y1st,y1en); decomp%y1size=y1en-y1st+1
    call distribute(ny,p_col,y2st,y2en); decomp%y2size=y2en-y2st+1
    call distribute(nz,p_col,z2st,z2en); decomp%z2size=z2en-z2st+1

    ! generate partition information - starting/ending index etc.
    call partition(nx,ny,nz,coord,(/1,2,3/),p_row,p_col,decomp%xst,decomp%xen,decomp%xsz)
    call partition(nx,ny,nz,coord,(/2,1,3/),p_row,p_col,decomp%yst,decomp%yen,decomp%ysz)
    call partition(nx,ny,nz,coord,(/2,3,1/),p_row,p_col,decomp%zst,decomp%zen,decomp%zsz)

    !============  Add following lines by Gong Zheng, transpose_x_z purpose ==========
    allocate(decomp%xcnts_xz(nproc),decomp%zcnts_xz(nproc))
    allocate(decomp%xtypes_xz(nproc),decomp%ztypes_xz(nproc))
    allocate(decomp%xdispls_xz(nproc))
    allocate(decomp%zdispls_xz(nproc))
    ! ================== adding ends ==================

#ifdef EVEN
    ! MPI_ALLTOALL buffer information
    ! For unevenly distributed data, pad smaller messages. Note the 
    ! last blocks along pencils always get assigned more mesh points
    decomp%xycount = decomp%x1size(p_row-1)*decomp%y1size(p_row-1)*decomp%xsz(3)
    decomp%yzcount = decomp%y2size(p_col-1)*decomp%z2size(p_col-1)*decomp%zsz(1)
#else
    ! Prepare send/receive buffer displacement and count for ALLTOALLV
    allocate(decomp%x1cnts(0:p_row-1),decomp%y1cnts(0:p_row-1))
    allocate(decomp%y2cnts(0:p_col-1),decomp%z2cnts(0:p_col-1))
    allocate(decomp%x1disp(0:p_row-1),decomp%y1disp(0:p_row-1))
    allocate(decomp%y2disp(0:p_col-1),decomp%z2disp(0:p_col-1))
    
    ! MPI_ALLTOALLV buffer information
    j=0; k=0
    do i=0, p_row-1
      decomp%x1cnts(i)= decomp%x1size(i)*decomp%xsz(2)*decomp%xsz(3)
      decomp%y1cnts(i)= decomp%ysz(1)*decomp%y1size(i)*decomp%ysz(3)
      decomp%x1disp(i)= j; j=j+decomp%x1cnts(i)
      decomp%y1disp(i)= k; k=k+decomp%y1cnts(i)
    enddo

    j=0; k=0
    do i=0,p_col-1
      decomp%y2cnts(i)= decomp%ysz(1)*decomp%y2size(i)*decomp%ysz(3)
      decomp%z2cnts(i)= decomp%zsz(1)*decomp%zsz(2)*decomp%z2size(i)
      decomp%y2disp(i)= j; j=j+decomp%y2cnts(i)
      decomp%z2disp(i)= k; k=k+decomp%z2cnts(i)
    enddo
#endif
    ! Information for MPI_Alltoallw for complex X <=> Z transposes
    decomp%xdispls_xz=0
    decomp%zdispls_xz=0
    decomp%xcnts_xz=0
    decomp%zcnts_xz=0
    decomp%xtypes_xz(:)=MPI_INTEGER
    decomp%ztypes_xz(:)=MPI_INTEGER
#ifdef MPI3
    index_src=0
    index_dest=0
#endif
    do k=0,p_row-1
      do i=0,p_col-1
        call MPI_Cart_rank(DECOMP_2D_COMM_CART,(/k,i/),prank,ierror)
        if(decomp%zst(2)<=y1en(k) .and. decomp%zen(2)>=y1st(k)) then
          decomp%zcnts_xz(prank+1)=1
          subsize_y=min(decomp%zen(2),y1en(k))-max(decomp%zst(2),y1st(k))+1
          offset_y =max(decomp%zst(2),y1st(k))-decomp%zst(2)
#ifdef MPI3
          index_src=index_src+1
          zranks(index_src)=prank
          zweights(index_src)=decomp%zsz(1)*subsize_y*decomp%z2size(i)
#endif
          call MPI_Type_create_subarray(3,decomp%zsz,(/decomp%zsz(1),subsize_y,decomp%z2size(i)/), &
            (/0,offset_y,z2st(i)-decomp%zst(3)/),MPI_ORDER_FORTRAN,real_type,decomp%ztypes_xz(prank+1),ierror)
          call MPI_Type_commit(decomp%ztypes_xz(prank+1),ierror)
!JD send to process with x-pencil defined by (k,i)
!JD x-bounds are taken from the z-pencils
!             send: decomp%zst(1):decomp%zen(1)
!JD y-bounds are the overlapping region of both pencils.
!                   max(decomp%zst(2),decomp%y1st(k)):min(decomp%zen(2),decomp%y1en(k))
!JD z-bounds are taken from the x-pencils.
!                   decomp%z2st(i):decomp%z2en(i)
        endif
        if(decomp%xst(2)<=y2en(i) .and.decomp%xen(2)>=y2st(i)) then
          decomp%xcnts_xz(prank+1)=1
          subsize_y=min(decomp%xen(2),y2en(i))-max(decomp%xst(2),y2st(i))+1
          offset_y =max(decomp%xst(2),y2st(i))-decomp%xst(2)
#ifdef MPI3
          index_dest=index_dest+1
          xranks(index_dest)=prank
          xweights(index_dest)=decomp%x1size(k)*subsize_y*decomp%xsz(3)
#endif
          call MPI_Type_create_subarray(3,decomp%xsz,(/decomp%x1size(k),subsize_y,decomp%xsz(3)/), &
            (/x1st(k)-decomp%xst(1),offset_y,0/),MPI_ORDER_FORTRAN,real_type,decomp%xtypes_xz(prank+1),ierror)
          call MPI_Type_commit(decomp%xtypes_xz(prank+1),ierror)
!JD recv from process with z-pencil defined by (k,i)
!JD x-bounds are taken from the z-pencils
!             send: x1st(k):x1en(k)
!JD y-bounds are the overlapping region of both pencils.
!                   max(decomp%xst(2),decomp%y2st(i)):min(xen(2),y2en(i))
!JD z-bounds are taken from the x-pencils.
!                   decomp%xst(3):decomp%xen(3)
        endif
      enddo
    enddo
#ifdef MPI3
    allocate(decomp%xranks(index_dest))
    allocate(decomp%zranks(index_src))
    decomp%xranks=xranks(1:index_dest)+1
    decomp%zranks=zranks(1:index_src)+1
    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART,index_src,zranks(1:index_src),zweights(1:index_src), &
      index_dest,xranks(1:index_dest),xweights(1:index_dest),MPI_inFO_NULL,.true.,decomp%xtozNeighborComm,ierror)
    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART,index_dest,xranks(1:index_dest),xweights(1:index_dest), &
      index_src,zranks(1:index_src),zweights(1:index_src),MPI_inFO_NULL,.true.,decomp%ztoxNeighborComm,ierror)
#endif

    ! allocate memory for the MPI_ALLTOALL(V) buffers
    ! define the buffers globally for performance reason
    buf_size= max(decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3), &
      max(decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3),decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)))
#ifdef EVEN
    ! padded alltoall optimisation may need larger buffer space
    buf_size=max(buf_size,max(decomp%xycount*p_row,decomp%yzcount*p_col) ) 
#endif

    ! check if additional memory is required
    ! *** TODO: consider how to share the real/complex buffers 
    if(buf_size>decomp_buf_size) then
       decomp_buf_size = buf_size
       if(allocated(work1_r)) deallocate(work1_r)
       if(allocated(work2_r)) deallocate(work2_r)
       allocate(work1_r(buf_size), STAT=ierror1)
       allocate(work2_r(buf_size), STAT=ierror2)
       if(abs(ierror1)+abs(ierror2) /= 0) then
         errorcode = 2
         call decomp_2d_abort(errorcode,'Out of memory when allocating 2DECOMP workspace')
       endif
    endif
  end subroutine decomp_info_init

  !******************************************************************
  ! Release memory associated with a DECOMP_INFO object
  !******************************************************************
  subroutine decomp_info_finalize(decomp)
    implicit none
    type(decomp_info),intent(inout)::decomp

    ! locals
    integer::i,ierror

    do i=1,nproc
      if (decomp%ztypes_xz(i).ne.MPI_INTEGER) then
        call MPI_TYPE_FREE(decomp%ztypes_xz(i),ierror)
      endif
      if (decomp%xtypes_xz(i).ne.MPI_INTEGER) then
        call MPI_TYPE_FREE(decomp%xtypes_xz(i),ierror)
      endif
    enddo
    deallocate(decomp%xcnts_xz,decomp%zcnts_xz,decomp%xtypes_xz)
    deallocate(decomp%ztypes_xz,decomp%xdispls_xz,decomp%zdispls_xz)
    deallocate(decomp%x1size,decomp%y1size,decomp%y2size,decomp%z2size)
#ifndef EVEN
    deallocate(decomp%x1cnts,decomp%y1cnts,decomp%y2cnts,decomp%z2cnts)
    deallocate(decomp%x1disp,decomp%y1disp,decomp%y2disp,decomp%z2disp)
#endif
#ifdef MPI3
    deallocate(decomp%xranks,decomp%zranks)
#endif
  end subroutine decomp_info_finalize

  !******************************************************************
  ! Find sub-domain information held by current processor
  !******************************************************************
  subroutine partition(nx,ny,nz,coord,pdim,p_row,p_col,lstart,lend,lsize)
    implicit none
    integer,intent(in)::nx,ny,nz,coord(2),pdim(3),p_row,p_col
    integer,dimension(3),intent(out)::lstart,lend,lsize

    ! locals
    integer::i,gsize
    integer,allocatable,dimension(:)::st,en

    do i= 1, 3
      if(i==1) then
        gsize= nx
      elseif(i==2) then
        gsize= ny
      elseif (i==3) then
        gsize= nz
      endif

      if(pdim(i) == 1) then        ! all local
        lstart(i)= 1
        lend(i)  = gsize
      elseif(pdim(i) == 2) then    ! distribute across p_row
        allocate(st(0:p_row-1))
        allocate(en(0:p_row-1))
        call distribute(gsize,p_row,st,en)
        lstart(i)= st(coord(1))
        lend(i)  = en(coord(1))
        deallocate(st,en)
      elseif(pdim(i) == 3) then    ! distribute across p_col
        allocate(st(0:p_col-1))
        allocate(en(0:p_col-1))
        call distribute(gsize,p_col,st,en)
        lstart(i)= st(coord(2))
        lend(i)  = en(coord(2))
        deallocate(st,en)
      endif
      lsize(i)=lend(i)-lstart(i)+1
    enddo
  end subroutine partition
 
  !******************************************************************
  ! Distibutes grid points in one dimension
  !******************************************************************
  subroutine distribute(ndata,proc,st,en)
    implicit none
    integer::ndata,proc,st(0:proc-1),en(0:proc-1)

    ! locals
    integer::i,isize,nl

    isize=ndata/proc
    nl=proc-mod(ndata,proc)
    do i=0,proc-1
      st(i)=(i*isize+1)+max(i-nl,0)
      en(i)=(i+1)*isize+max(i+1-nl,0)
    enddo
  end subroutine distribute

  !******************************************************************
  ! Auto-tuning algorithm to select the best 2D processor grid
  !******************************************************************
  subroutine best_2d_grid(nx,ny,nz,iproc,best_p_row,best_p_col)
    implicit none
    integer,intent(in) :: nx,ny,nz,iproc
    integer,intent(out):: best_p_row, best_p_col

    ! locals
    type(decomp_info)::decomp
    real(mytype):: t1,t2,best_time
    integer:: nfact,i,row,col,ierror,errorcode
    integer,allocatable,dimension(:) :: factors
    real(mytype),allocatable,dimension(:,:,:) :: u1, u2, u3

    if(nrank==0) write(*,*) 'In auto-tuning mode......'
    best_time = huge(t1)
    best_p_row= -1
    best_p_col= -1

    i=int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if(nrank==0)write(*,*) 'factors: ', (factors(i), i=1,nfact)

    do i=1, nfact
      row = factors(i)
      col = iproc / row
      if(min(nx,ny)>=row .and.min(ny,nz)>=col) then
        np_row= row
        np_col= col
        call decomp_info_init(nx,ny,nz,row,col,decomp)

        allocate(u1(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
        allocate(u2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
        allocate(u3(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))
        t1 = MPI_WTIME()
        call transpose_x_to_y(u1,u2,decomp)
        call transpose_y_to_z(u2,u3,decomp)        
        call transpose_z_to_y(u3,u2,decomp)
        call transpose_y_to_x(u2,u1,decomp)
        t2 = MPI_WTIME() - t1
        deallocate(u1,u2,u3)
        call decomp_info_finalize(decomp)

        call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
        t1=t1/real(nproc,mytype)
        if(nrank==0)write(*,*) 'processor grid', row, ' by ', col, ' time=', t1

        if(best_time > t1) then
          best_time = t1
          best_p_row= row
          best_p_col= col
        endif
      endif
    enddo
    deallocate(factors)

    if(best_p_row/=-1) then
     if(nrank==0)write(*,*)'the best processor grid is probably ',best_p_row,' by ',best_p_col
    else
     errorcode = 9
     call decomp_2d_abort(errorcode,'The processor-grid auto-tuning code failed. ' // &
       'The number of processes requested is probably too large.')
    endif
  end subroutine best_2d_grid

  !******************************************************************
  ! Error handling
  !******************************************************************
  subroutine decomp_2d_abort(errorcode,msg)
    implicit none
    integer,intent(in)::errorcode
    character(len=*),intent(in)::msg

    ! locals
    integer::ierror

    if(nrank==0) then
      write(*,*) '2DECOMP&FFT ERROR - errorcode: ', errorcode
      write(*,*) 'ERROR MESSAGE: ' // msg
    endif
    call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
  end subroutine decomp_2d_abort

  !******************************************************************
  ! findfactor
  !******************************************************************
  subroutine findfactor(num, factors, nfact)
    implicit none
    integer,intent(in)::num
    integer,intent(out),dimension(*)::factors
    integer,intent(out)::nfact

    ! locals
    integer :: i, m

    ! find the factors <= sqrt(num)
    m = int(sqrt(real(num)))
    nfact = 1
    do i=1,m
      if(num/i*i /= num)cycle
      factors(nfact)= i
      nfact = nfact + 1
    enddo
    nfact = nfact - 1

    ! derive those > sqrt(num)
    if(factors(nfact)**2/=num) then
      do i=nfact+1, 2*nfact
        factors(i) = num/factors(2*nfact-i+1)
      enddo
      nfact = nfact * 2
    else
      do i=nfact+1, 2*nfact-1
        factors(i)= num/factors(2*nfact-i)
      enddo
      nfact = nfact * 2 - 1
    endif
  end subroutine findfactor

  !******************************************************************
  ! transpose_x_to_y
  !******************************************************************  
  subroutine transpose_x_to_y(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) :: src
    real(mytype),dimension(:,:,:),intent(out):: dst
    type(decomp_info),intent(in), optional :: opt_decomp

    ! locals
    type(decomp_info) :: decomp
    integer:: nsize,s2,s3,d1,d3,ierror,i,j,k,m,i1,i2,pos

    if(present(opt_decomp)) then
      decomp = opt_decomp
    else
      decomp = decomp_main
    endif
    s2= SIZE(src,2)
    s3= SIZE(src,3)
    d1= SIZE(dst,1)
    d3= SIZE(dst,3)

    ! Rearrange source array as send buffer
    i1=1;i2=0
    do m=0,np_row-1
      nsize=decomp%x1size(m)
      i2=i2+nsize
#ifdef EVEN
      pos=m*decomp%xycount
#else
      pos=decomp%x1disp(m)
#endif
      do k=1,s3
        do j=1,s2
          do i=i1,i2
            pos = pos + 1
            work1_r(pos) = src(i,j,k)
          enddo
        enddo
      enddo
      i1=i1+nsize
    enddo

    ! Transpose using MPI_ALLTOALL(V)
#ifdef EVEN
    call MPI_ALLTOALL(work1_r, decomp%xycount, &
      real_type, work2_r, decomp%xycount,real_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_r, decomp%x1cnts, decomp%x1disp, &
      real_type, work2_r, decomp%y1cnts, decomp%y1disp,real_type, DECOMP_2D_COMM_COL, ierror)
#endif

    ! Rearrange receive buffer
    i1=1;i2=0
    do m=0,np_row-1
      nsize=decomp%y1size(m)
      i2=i2+nsize
#ifdef EVEN
      pos =m*decomp%xycount
#else
      pos =decomp%y1disp(m)
#endif
      do k=1,d3
        do j=i1,i2
          do i=1,d1
            pos = pos + 1
            dst(i,j,k)= work2_r(pos)
          enddo
        enddo
      enddo
      i1=i1+nsize
    enddo
  end subroutine transpose_x_to_y

  !******************************************************************
  ! transpose_y_to_x
  !******************************************************************
  subroutine transpose_y_to_z(src, dst, opt_decomp)
    implicit none
    real(mytype), dimension(:,:,:), intent(in) :: src
    real(mytype), dimension(:,:,:), intent(out) :: dst
    type(decomp_info), intent(in), optional :: opt_decomp

    ! locals
    type(decomp_info) :: decomp
    integer:: nsize,s1,s3,d1,d2,ierror,i,j,k, m,i1,i2,pos

    if(present(opt_decomp)) then
      decomp= opt_decomp
    else
      decomp= decomp_main
    endif
    s1= SIZE(src,1)
    s3= SIZE(src,3)
    d1= SIZE(dst,1)
    d2= SIZE(dst,2)

    ! Rearrange source array as send buffer
    i1=1;i2=0
    do m=0,np_col-1
      nsize=decomp%y2size(m)
      i2=i2+nsize
#ifdef EVEN
       pos = m*decomp%yzcount
#else
       pos = decomp%y2disp(m)
#endif
       do k=1,s3
         do j=i1,i2
           do i=1,s1
           pos = pos + 1
           work1_r(pos) = src(i,j,k)
          enddo
        enddo
      enddo
      i1=i1+nsize
    enddo

#ifdef EVEN
    if(decomp%even) then
      call MPI_ALLTOALL(work1_r, decomp%yzcount, &
        real_type, dst, decomp%yzcount,real_type, DECOMP_2D_COMM_ROW, ierror)
    else
      call MPI_ALLTOALL(work1_r, decomp%yzcount, &
        real_type, work2_r, decomp%yzcount,real_type, DECOMP_2D_COMM_ROW, ierror)
    endif
#else
    call MPI_ALLTOALLV(work1_r, decomp%y2cnts, decomp%y2disp, &
      real_type, dst, decomp%z2cnts, decomp%z2disp,real_type, DECOMP_2D_COMM_ROW, ierror)
#endif

    ! Rearrange receive buffer
#ifdef EVEN
    if(.not. decomp%even) then
      i1=1;i2=0
      do m=0,np_col-1
        nsize=decomp%z2size(m)
        i2=i2+nsize
        pos = m * decomp%yzcount
        do k=i1,i2
          do j=1,d2
            do i=1,d1
              pos = pos + 1
              dst(i,j,k) = work2_r(pos)
            enddo
          enddo
        enddo
        i1=i1+nsize
      enddo
    endif
#else
    ! Note the receive buffer is already in natural (i,j,k) order
    ! So no merge operation needed
#endif
  end subroutine transpose_y_to_z

  !******************************************************************
  ! transpose_z_to_x
  !******************************************************************
  subroutine transpose_z_to_x(src, dst, opt_decomp)
    implicit none
    real(mytype), dimension(:,:,:), intent(in) :: src
    real(mytype), dimension(:,:,:), intent(out) :: dst
    type(decomp_info), intent(in), optional :: opt_decomp

    ! locals
    integer::ierror
    type(decomp_info):: decomp

    if(present(opt_decomp))then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

#ifdef MPI3
    call MPI_Neighbor_alltoallw( &
      src,decomp%zcnts_xz(decomp%zranks),decomp%zdispls_xz(decomp%zranks),decomp%ztypes_xz(decomp%zranks), &
      dst,decomp%xcnts_xz(decomp%xranks),decomp%xdispls_xz(decomp%xranks),decomp%xtypes_xz(decomp%xranks), &
      decomp%ztoxNeighborComm,ierror)
#else
    call MPI_Alltoallw(src,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xz, &
      dst,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xz,MPI_COMM_WORLD,ierror)
#endif
  end subroutine transpose_z_to_x

  !******************************************************************
  ! transpose_x_to_z
  !******************************************************************
  subroutine transpose_x_to_z(src, dst, opt_decomp)
    implicit none
    real(mytype), dimension(:,:,:), intent(in) :: src
    real(mytype), dimension(:,:,:), intent(out) :: dst
    type(decomp_info), intent(in), optional :: opt_decomp

    ! locals
    integer :: ierror
    type(decomp_info) :: decomp

    if(present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    endif

#ifdef MPI3
    call MPI_Neighbor_alltoallw( &
      src,decomp%xcnts_xz(decomp%xranks),decomp%xdispls_xz(decomp%xranks),decomp%xtypes_xz(decomp%xranks), &
      dst,decomp%zcnts_xz(decomp%zranks),decomp%zdispls_xz(decomp%zranks),decomp%ztypes_xz(decomp%zranks), &
      decomp%xtozNeighborComm,ierror)
#else
    call MPI_Alltoallw(src,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xz, &
      dst,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xz,MPI_COMM_WORLD,ierror)
#endif
  end subroutine transpose_x_to_z

  !******************************************************************
  ! transpose_z_to_y
  !******************************************************************
  subroutine transpose_z_to_y(src, dst, opt_decomp)
    implicit none
    real(mytype), dimension(:,:,:), intent(in) :: src
    real(mytype), dimension(:,:,:), intent(out) :: dst
    type(decomp_info), intent(in), optional :: opt_decomp

    ! locals
    type(decomp_info) :: decomp
    integer::nsize,s1,s2,d1,d3,ierror,i,j,k, m,i1,i2,pos

    if(present(opt_decomp)) then
      decomp = opt_decomp
    else
      decomp = decomp_main
    endif
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    d1 = SIZE(dst,1)
    d3 = SIZE(dst,3)

    ! Rearrange source array as send buffer
#ifdef EVEN
    if(.not. decomp%even)then
      i1=1;i2=0
      do m=0,np_col-1
        nsize=decomp%z2size(m)
        i2=i2+nsize
        pos= m * decomp%yzcount
        do k=i1,i2
          do j=1,s2
            do i=1,s1
              pos = pos + 1
              work1_r(pos) = src(i,j,k)
            enddo
          enddo
        enddo
        i1=i1+nsize
      enddo
    endif
#else
    ! Note the src array is suitable to be a send buffer
    ! So no split operation needed
#endif

#ifdef EVEN
    if(decomp%even) then
      call MPI_ALLTOALL(src, decomp%yzcount, &
        real_type, work2_r, decomp%yzcount,real_type, DECOMP_2D_COMM_ROW, ierror)
    else
      call MPI_ALLTOALL(work1_r, decomp%yzcount, &
        real_type, work2_r, decomp%yzcount,real_type, DECOMP_2D_COMM_ROW, ierror)
    endif
#else
    call MPI_ALLTOALLV(src, decomp%z2cnts, decomp%z2disp, &
      real_type, work2_r, decomp%y2cnts, decomp%y2disp,real_type, DECOMP_2D_COMM_ROW, ierror)
#endif

    ! Rearrange receive buffer
    i1=1;i2=0
    do m=0,np_col-1
      nsize=decomp%y2size(m)
      i2=i2+nsize
#ifdef EVEN
      pos = m * decomp%yzcount
#else
      pos = decomp%y2disp(m)
#endif
      do k=1,d3
        do j=i1,i2
          do i=1,d1
            pos = pos + 1
            dst(i,j,k) = work2_r(pos)
          enddo
        enddo
      enddo
      i1=i1+nsize
    enddo
  end subroutine transpose_z_to_y

  !******************************************************************
  ! transpose_y_to_x
  !******************************************************************
  subroutine transpose_y_to_x(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) :: src
    real(mytype),dimension(:,:,:),intent(out) :: dst
    type(decomp_info),intent(in), optional :: opt_decomp

    ! locals
    type(decomp_info)::decomp
    integer::nsize,s1,s3,d2,d3,ierror,i,j,k,m,i1,i2,pos

    if(present(opt_decomp)) then
      decomp = opt_decomp
    else
      decomp = decomp_main
    endif
    s1= SIZE(src,1)
    s3= SIZE(src,3)
    d2= SIZE(dst,2)
    d3= SIZE(dst,3)

    ! Rearrange source array as send buffer
    i1=1;i2=0
    do m=0,np_row-1
      nsize=decomp%y1size(m)
      i2=i2+nsize
#ifdef EVEN
      pos= m * decomp%xycount + 1
#else
      pos= decomp%y1disp(m) + 1
#endif
      do k=1,s3
        do j=i1,i2
          do i=1,s1
            work1_r(pos) = src(i,j,k)
            pos= pos + 1
          enddo
        enddo
      enddo
      i1=i1+nsize
    enddo

    ! Transpose using MPI_ALLTOALL(V)
#ifdef EVEN
    call MPI_ALLTOALL(work1_r, decomp%xycount, &
      real_type, work2_r, decomp%xycount,real_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
      real_type, work2_r, decomp%x1cnts, decomp%x1disp,real_type, DECOMP_2D_COMM_COL, ierror)
#endif

    ! Rearrange receive buffer
    i1=1;i2=0
    do m=0,np_row-1
      nsize=decomp%x1size(m)
      i2=i2+nsize
#ifdef EVEN
      pos= m * decomp%xycount
#else
      pos= decomp%x1disp(m)
#endif
     do k=1,d3
       do j=1,d2
         do i=i1,i2
           pos = pos + 1
           dst(i,j,k) = work2_r(pos)
          enddo
        enddo
      enddo
      i1=i1+nsize
    enddo
  end subroutine transpose_y_to_x

  !******************************************************************
  ! Write a 3D array as part of a big MPI-IO file.
  !  'disp' will be updated after the writing operation.
  !******************************************************************
  subroutine decomp_2d_write_var(fh,disp,ipencil,var,opt_decomp)
    implicit none
    integer,intent(in)::fh
    integer,intent(in)::ipencil
    integer(MPI_OFFSET_KIND),intent(inout)::disp
    real(mytype),dimension(:,:,:),intent(in)::var
    type(decomp_info),intent(in),optional::opt_decomp

    ! locals
    type(decomp_info)::decomp
    integer::ierror,data_type,newtype
    integer,dimension(3)::sizes,subsizes,starts

    data_type = real_type
    if(present(opt_decomp))then
      decomp = opt_decomp
    else
      decomp = decomp_main
    endif
    sizes(1)= decomp%xsz(1)
    sizes(2)= decomp%ysz(2)
    sizes(3)= decomp%zsz(3)
    if(ipencil==1)then
      subsizes= decomp%xsz
      starts= decomp%xst-1
    elseif(ipencil==2)then
      subsizes= decomp%ysz
      starts= decomp%yst-1
    elseif(ipencil==3)then
      subsizes= decomp%zsz
      starts= decomp%zst-1
    endif
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,data_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,var,subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8) ! Update displacement
  end subroutine decomp_2d_write_var

  !******************************************************************
  ! Read a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
  !******************************************************************
  subroutine decomp_2d_read_var(fh,disp,ipencil,var,opt_decomp)
    implicit none
    integer,intent(in) :: fh
    integer,intent(in)::ipencil
    integer(MPI_OFFSET_KIND),intent(inout)::disp
    real(mytype),dimension(:,:,:),intent(inout)::var
    type(decomp_info),intent(in),optional::opt_decomp

    ! locals
    type(decomp_info)::decomp
    integer::ierror,data_type,newtype
    integer,dimension(3)::sizes,subsizes,starts

    data_type = real_type
    if(present(opt_decomp))then
      decomp = opt_decomp
    else
      decomp = decomp_main
    endif
    sizes(1)= decomp%xsz(1)
    sizes(2)= decomp%ysz(2)
    sizes(3)= decomp%zsz(3)
    if(ipencil==1)then
      subsizes= decomp%xsz
      starts= decomp%xst-1
    elseif(ipencil==2)then
      subsizes= decomp%ysz
      starts= decomp%yst-1
    elseif(ipencil==3)then
      subsizes= decomp%zsz
      starts= decomp%zst-1
    endif
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,data_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh,var,subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp=disp + int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)  ! Update displacement
  end subroutine decomp_2d_read_var

  !******************************************************************
  ! Write a 2D slice of the 3D data to a file
  ! It is much easier to implement if all mpi ranks participate I/O.
  ! Transpose the 3D data if necessary.
  !******************************************************************
  subroutine decomp_2d_write_plane(ipencil,var,iplane,n,filename,opt_decomp)
    implicit none
    integer,intent(in)::n      !which plane to write (global coordinate)
    integer,intent(in)::iplane !x-plane= 1; y-plane= 2; z-plane= 3
    integer,intent(in)::ipencil!x-pencil=1; y-pencil=2; z-pencil=3
    real(mytype),dimension(:,:,:),intent(in) :: var
    character(len=*),intent(in)::filename
    type(decomp_info),intent(in),optional::opt_decomp

    ! lcoals
    type(decomp_info)::decomp
    integer(MPI_OFFSET_KIND)::disp
    real(mytype),allocatable,dimension(:,:,:)::wk
    real(mytype_save),allocatable,dimension(:,:,:)::wk2d
    integer::i,j,k,ierror,newtype,fh,data_type
    integer,dimension(3)::sizes,subsizes,starts

    data_type= real_type_save
    if(present(opt_decomp))then
      decomp= opt_decomp
    else
      decomp= decomp_main
    endif
    if(iplane==1)then
      allocate(wk2d(1,decomp%xsz(2),decomp%xsz(3)))
      IF(ipencil==1)THEN
        do k=1,decomp%xsz(3)
          do j=1,decomp%xsz(2)
            wk2d(1,j,k)=real(var(n,j,k),mytype_save)
          enddo
        enddo
      ELSE
        allocate(wk(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
        if(ipencil==2)then
          call transpose_y_to_x(var,wk,decomp)
        elseif(ipencil==3)then
          call transpose_z_to_x(var,wk,decomp)
        endif
        do k=1,decomp%xsz(3)
          do j=1,decomp%xsz(2)
            wk2d(1,j,k)=real(wk(n,j,k),mytype_save)
          enddo
        enddo
        deallocate(wk)
      ENDIF
      sizes   =(/1,decomp%ysz(2),decomp%zsz(3)/)
      subsizes=(/1,decomp%xsz(2),decomp%xsz(3)/)
      starts  =(/1,decomp%xst(2),decomp%xst(3)/)-1
    elseif(iplane==2) then
      allocate(wk2d(decomp%ysz(1),1,decomp%ysz(3)))
      IF(ipencil==2)THEN
        do k=1,decomp%ysz(3)
          do i=1,decomp%ysz(1)
            wk2d(i,1,k)=real(var(i,n,k),mytype_save)
          enddo
        enddo
      ELSE
        allocate(wk(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
        if(ipencil==1)then
          call transpose_x_to_y(var,wk,decomp)
        elseif(ipencil==3)then
          call transpose_z_to_y(var,wk,decomp)
        endif
        do k=1,decomp%ysz(3)
          do i=1,decomp%ysz(1)
            wk2d(i,1,k)=real(wk(i,n,k),mytype_save)
          enddo
        enddo
        deallocate(wk)
      ENDIF
      sizes   =(/decomp%xsz(1),1,decomp%zsz(3)/)
      subsizes=(/decomp%ysz(1),1,decomp%ysz(3)/)
      starts  =(/decomp%yst(1),1,decomp%yst(3)/)-1
    elseif(iplane==3)then
      allocate(wk2d(decomp%zsz(1),decomp%zsz(2),1))
      IF(ipencil==3)THEN
        do j=1,decomp%zsz(2)
          do i=1,decomp%zsz(1) 
            wk2d(i,j,1)=real(var(i,j,n),mytype_save)
          enddo
        enddo
      ELSE
        allocate(wk(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))
        if(ipencil==1)then
          call transpose_x_to_z(var,wk,decomp)
        elseif(ipencil==2)then
          call transpose_y_to_z(var,wk,decomp)
        endif
        do j=1,decomp%zsz(2)
          do i=1,decomp%zsz(1) 
            wk2d(i,j,1)=real(wk(i,j,n),mytype_save)
          enddo
        enddo
        deallocate(wk)
      ENDIF
      sizes   =(/decomp%xsz(1),decomp%ysz(2),1/)
      subsizes=(/decomp%zsz(1),decomp%zsz(2),1/)
      starts  =(/decomp%zst(1),decomp%zst(2),1/)-1
    endif
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh,ierror)
    call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,wk2d,subsizes(1)*subsizes(2)*subsizes(3),data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    deallocate(wk2d)
  end subroutine decomp_2d_write_plane

  !******************************************************************
  ! Write 3D array data for every specified mesh point
  !******************************************************************
  subroutine decomp_2d_write_every(ipencil,var,iskip,jskip,kskip,filename,from1)
    implicit none
    integer,intent(in)::ipencil
    integer,intent(in)::iskip,jskip,kskip 
    character(len=*),intent(in)::filename
    real(mytype),dimension(:,:,:),intent(in)::var
    logical,intent(in)::from1  ! T:save 1,n+1,2n+1...; F:save n,2n,3n...

    ! locals
    integer(MPI_OFFSET_KIND)::disp
    real(mytype_save),allocatable,dimension(:,:,:)::wk
    integer::i,j,k,id,jd,kd,ierror,newtype,fh,key,color,newcomm,data_type
    integer,dimension(3)::sizes,subsizes,starts,xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    data_type = real_type_save
    skip=(/iskip,jskip,kskip/)
    do i=1,3
      if(from1)then
        xst(i)= (xstart(i)+skip(i)-1)/skip(i)
        yst(i)= (ystart(i)+skip(i)-1)/skip(i)
        zst(i)= (zstart(i)+skip(i)-1)/skip(i)
        if(mod(xstart(i)+skip(i)-1,skip(i))/=0) xst(i)=xst(i)+1
        if(mod(ystart(i)+skip(i)-1,skip(i))/=0) yst(i)=yst(i)+1
        if(mod(zstart(i)+skip(i)-1,skip(i))/=0) zst(i)=zst(i)+1
        xen(i)= (xend(i)+skip(i)-1)/skip(i)
        yen(i)= (yend(i)+skip(i)-1)/skip(i)
        zen(i)= (zend(i)+skip(i)-1)/skip(i)
      else
        xst(i)= xstart(i)/skip(i)
        yst(i)= ystart(i)/skip(i)
        zst(i)= zstart(i)/skip(i)
        if(mod(xstart(i),skip(i))/=0) xst(i)=xst(i)+1
        if(mod(ystart(i),skip(i))/=0) yst(i)=yst(i)+1
        if(mod(zstart(i),skip(i))/=0) zst(i)=zst(i)+1
        xen(i)= xend(i)/skip(i)
        yen(i)= yend(i)/skip(i)
        zen(i)= zend(i)/skip(i)
      endif
      xsz(i)= xen(i)-xst(i)+1
      ysz(i)= yen(i)-yst(i)+1
      zsz(i)= zen(i)-zst(i)+1
    enddo

    ! if 'skip' value is large it is possible that some ranks do not 
    ! contain any points to be written. Subarray constructor requires 
    ! nonzero size so it is not possible to use MPI_COMM_WORLD for IO.
    ! Create a sub communicator for this...
    color=1; key=0  ! rank order doesn't matter
    if(ipencil==1)then
      if(xsz(1)==0.or.xsz(2)==0.or.xsz(3)==0)color=2
    elseif(ipencil==2)then
      if(ysz(1)==0.or.ysz(2)==0.or.ysz(3)==0)color=2
    elseif(ipencil==3)then
      if(zsz(1)==0.or.zsz(2)==0.or.zsz(3)==0)color=2
    endif
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,newcomm,ierror)

    if(color/=1) return ! only ranks in this group do IO collectively
    sizes(1)= xsz(1)
    sizes(2)= ysz(2)
    sizes(3)= zsz(3)
    if(ipencil==1) then
      subsizes= xsz
      starts= xst-1
    elseif(ipencil==2) then
      subsizes= ysz
      starts= yst-1
    elseif(ipencil==3) then
      subsizes= zsz
      starts= zst-1
    endif

    ! Copy data from original array
    IF(ipencil==1)THEN
      allocate(wk(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
      if(from1)then
        do k=xst(3),xen(3)
          kd=(k-1)*kskip-xstart(3)+2
          do j=xst(2),xen(2)
            jd=(j-1)*jskip-xstart(2)+2
            do i=xst(1),xen(1)
              id=(i-1)*iskip-xstart(1)+2
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      else
        do k=xst(3),xen(3)
          kd=k*kskip-xstart(3)+1
          do j=xst(2),xen(2)
            jd=j*jskip-xstart(2)+1
            do i=xst(1),xen(1)
              id=i*iskip-xstart(1)+1
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      endif   
    ELSEIF(ipencil==2)THEN
      allocate(wk(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))
      if(from1)then
        do k=yst(3),yen(3)
          kd=(k-1)*kskip-ystart(3)+2
          do j=yst(2),yen(2)
            jd=(j-1)*jskip-ystart(2)+2
            do i=yst(1),yen(1)
              id=(i-1)*iskip-ystart(1)+2
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      else
        do k=yst(3),yen(3)
          kd=k*kskip-ystart(3)+1
          do j=yst(2),yen(2)
            jd=j*jskip-ystart(2)+1
            do i=yst(1),yen(1)
              id=i*iskip-ystart(1)+1
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      endif
    ELSEIF(ipencil==3)THEN
      allocate(wk(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))
      if(from1)then
        do k=zst(3),zen(3)
          kd=(k-1)*kskip-zstart(3)+2
          do j=zst(2),zen(2)
            jd=(j-1)*jskip-zstart(2)+2
            do i=zst(1),zen(1)
              id=(i-1)*iskip-zstart(1)+2
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
         enddo
      else
        do k=zst(3),zen(3)
          kd=k*kskip-zstart(3)+1
          do j=zst(2),zen(2)
            jd=j*jskip-zstart(2)+1
            do i=zst(1),zen(1)
              id=i*iskip-zstart(1)+1
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      endif
    ENDIF
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(newcomm,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierror)
    call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,wk,subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    deallocate(wk)
  end subroutine decomp_2d_write_every

  !**********************************************************************
  ! myInit_neighbour
  !**********************************************************************
  subroutine myInit_neighbour(row,col,BcOpt)
    implicit none
    integer,intent(in)::row,col,BcOpt(6)

    !locals
    integer :: i,coord1,coord2,idTop,idBottom,idLeft,idRight

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
    !       |        2 * 1  |        2 * 1  |        2 * 1
    !       |_ _ _ z 7 4 8  |_ _ _z  7 4 8  |_ _ _y  7 4 8
    ! 
    !        Here * means the center subdomain, and 1-8 stands for the reduative location of the eight neighbors

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
    real(mytype), allocatable, dimension(:,:,:),intent(out) :: var
    type(MatBound),intent(INout)::mb
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
    real(mytype), dimension(mb%xmm:mb%xpm,mb%ymm:mb%ypm,mb%zmm:mb%zpm),intent(INout):: mat

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
    best_p_row= -1
    best_p_col= -1

    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if(nrank==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

    mb_test%pencil = pencil  
    mb_test%xme=lvlhalo_test;  mb_test%xpe=lvlhalo_test
    mb_test%yme=lvlhalo_test;  mb_test%ype=lvlhalo_test
    mb_test%zme=lvlhalo_test;  mb_test%zpe=lvlhalo_test

    hi_test%pencil = pencil
    hi_test%xmh=lvlhalo_test;  hi_test%xph=lvlhalo_test
    hi_test%ymh=lvlhalo_test;  hi_test%yph=lvlhalo_test
    hi_test%zmh=lvlhalo_test;  hi_test%zph=lvlhalo_test

    do i=1,nfact
      row= factors(i)
      col= iproc / row
      if(min(nx,ny)>=row .and. min(ny,nz)>=col) then
        call decomp_2d_init(nx,ny,nz,row,col)
        call myInit_neighbour(row,col,BcOpt)

        allocate(u1(xsize(1),xsize(2),xsize(3)))
        allocate(u2(ysize(1),ysize(2),ysize(3)))
        allocate(u3(zsize(1),zsize(2),zsize(3)))
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

        call myallocate(u1,mb_test)
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
        t1=t1/real(iproc,kind=mytype)

        if(nrank==0)write(*,*) 'processor grid', row, ' by ', col, ' time=', t1
        if(best_time > t1) then
          best_time = t1
          best_p_row= row
          best_p_col= col
        endif
      endif
    enddo
    deallocate(factors)

    if(best_p_row/=-1) then
      if(nrank==0)write(*,*)'the best processor grid is probably ',best_p_row,' by ',best_p_col
    else
      errorcode = 9
      call decomp_2d_abort(errorcode,'The processor-grid auto-tuning code failed. ' // &
        'The number of processes requested is probably too large.')
    endif
  end subroutine my_best_2d_grid

end module decomp_2d
