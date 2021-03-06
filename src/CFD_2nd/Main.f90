program main_channel3d
  use MPI
  use decomp_2d
  use mydecomp_2d_extra
  use m_ChannelSystem
  use m_Parameters
  use m_Variables
  use m_IOAndVisu
  use m_LogInfo
  implicit none
  character(len=64)::ChannelPrm
  integer:: ierror

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! read Channel options
  if(command_argument_count()/=1 .and. nrank==0) write(*,*)'command argument wrong!'
  call get_command_argument(1,ChannelPrm)
  call ReadAndInitParameters(ChannelPrm)

  !call my_best_2d_grid(nx,ny,nz,pencil,lvlhalo_test,time_halo,BcOpt,iproc,best_p_row,best_p_col)
  if(p_row==0 .and. p_col==0) call my_best_2d_grid(nxc,nyc,nzc,y_pencil,2, 8, BcOption, nproc, p_row, p_col)
  call decomp_2d_init(nxc,nyc,nzc,p_row,p_col)
  call myInit_neighbour(p_row,p_col,BcOption)

  call ChannelInitialize(ChannelPrm)    ! Topest level initialing for Channel body

  asso_Q: associate(Q_vor =>RealArr1)
  call  dump_visu(ifirst-1,ux,uy,uz,pressure,Q_vor)
  end associate asso_Q

  do itime=ifirst, ilast
    call ChannelIterate()
  enddo
  if(nrank==0)call MainLog%OutInfo("Good job! Channel3d finished successfully!",1)

  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
end program main_channel3d
