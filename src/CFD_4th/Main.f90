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
  integer:: ierror,BcOption(6)

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! read Channel options
  ChannelPrm = "./CFD/data_in/channel3d_4th.prm"
  call ReadAndInitParameters(ChannelPrm)

  !call my_best_2d_grid(nx,ny,nz,pencil,lvlhalo_test,time_halo,BcOpt,iproc,best_p_row,best_p_col)
  if(FlowType==FT_CH) then
    BcOption=(/0,0,-1,-1,0,0/)
  elseif(FlowType==FT_HC) then
    BcOption=(/0,0,-1,-2,0,0/)
  endif
  if(p_row==0 .and. p_col==0) call my_best_2d_grid(nxc,nyc,nzc,y_pencil,2,12,BcOption,nproc,p_row,p_col)
  call decomp_2d_init(nxc,nyc,nzc,p_row,p_col)
  call myInit_neighbour(p_row,p_col,BcOption)

  call ChannelInitialize(ChannelPrm)    ! Topest level initialing for Channel body

  asso_Q: associate(Q_vor =>RealArr1)
  call  dump_visu(ifirst-1,ux,uy,uz,pressure,Q_vor)
  end associate asso_Q

  do itime=ifirst, ilast
    call ChannelIterate()
  enddo
  if(nrank==0)call MainLog%OutInfo("Good job! Channel3d_4th finished successfully!",1)

  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
end program main_channel3d
