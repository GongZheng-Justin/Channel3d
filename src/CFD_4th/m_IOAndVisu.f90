module m_IOAndVisu
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_Variables,only:mb1
  use MPI
  use decomp_2d
  use m_MeshAndMetries
  use m_Tools,only: Clc_Q_vor,Clc_lamda2
  implicit none
  private

  ! VisuOption
  integer,save:: iskip,jskip,kskip
  integer,save:: XDMF_SET_TYPE =0  ! 0: Cell, 1: Node
  character(16),save:: XDMF_SET_TYPE_STR
  integer,save:: Prev_BackUp_itime  = 123456789
  logical,save:: save_ux,save_uy,save_uz,save_wx,save_wy,save_wz,save_wMag
  logical,save:: save_pr,save_Q_vor,save_lamda2,WriteHistOld,ReadHistOld

  public:: InitVisu, dump_visu, read_restart, write_restart, Delete_Prev_Restart

contains

  !******************************************************************
  ! InitVisu
  !******************************************************************
  subroutine InitVisu(ChannelPrm)
    implicit none
    character(*),intent(in)::ChannelPrm

    ! locals
    character(256)::XdmfFile
    integer::nUnitFile,ierror,nflds,ifld,iprec,i,j,k,nxOut,nyOut,nzOut
    NAMELIST /IO_Options/ save_ux,save_uy,save_uz,save_pr,save_wx,save_wy,save_wz,save_wMag,save_Q_vor, &
                          save_lamda2,WriteHistOld,ReadHistOld,XDMF_SET_TYPE,iskip,jskip,kskip
 
    nUnitFile = GetFileUnit() 
    open(unit=nUnitFile, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitVisu", "Cannot open file: "//trim(ChannelPrm))
    read(nUnitFile, nml=IO_Options)
    close(nUnitFile,IOSTAT=ierror)
    if(nrank==0)write(MainLog%nUnit, nml=IO_Options)

    ! write XDMF file
    if(nrank/=0) return
    if(XDMF_SET_TYPE==0) then ! 0: Cell
      XDMF_SET_TYPE_STR='" Center="Cell">'
      nxOut=nxp; nyOut=nyp; nzOut=nzp
    else                      ! 1: Node
      XDMF_SET_TYPE_STR='" Center="Node">'
      nxOut=nxc; nyOut=nyc; nzOut=nzc
    endif

    nUnitFile = GetFileUnit()
    XdmfFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//".xmf"
    open(unit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitVisu","Cannot open file:  "//trim(XdmfFile))
    ! XDMF/XMF Title
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '<Domain>'

    ! grid
    iprec=mytype_save
    write(nUnitFile,'(A,3I5,A)')'    <Topology name="TOPO" TopologyType="3DRectMesh" Dimensions="',nzOut,nyOut,nxOut,'"/>'
    write(nUnitFile,'(A)')'    <Geometry name="GEO" GeometryType="VXVYVZ">'
    ! x-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nxOut,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do i=1,nxOut
      if(XDMF_SET_TYPE==0) then
        write(nUnitFile,'(E14.7)',advance='no') (i-1)*dx
      else
        write(nUnitFile,'(E14.7)',advance='no') (i-1)*dx+dx*half
      endif
    enddo
    write(nUnitFile,'(A)')''; write(nUnitFile,'(A)')'        </DataItem>'
    ! y-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nyOut,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do j=1,nyOut
      if(XDMF_SET_TYPE==0) then
        write(nUnitFile,'(E14.7)',advance='no') yp(j)
      else
        write(nUnitFile,'(E14.7)',advance='no') yc(j)
      endif
    enddo
    write(nUnitFile,'(A)')''; write(nUnitFile,'(A)')'        </DataItem>'
    ! z-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nzOut,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do k=1,nzOut
      if(XDMF_SET_TYPE==0) then
        write(nUnitFile,'(E14.7)',advance='no') (k-1)*dz
      else
        write(nUnitFile,'(E14.7)',advance='no') (k-1)*dz+dz*half
      endif
    enddo
    write(nUnitFile,'(A)')''; write(nUnitFile,'(A)')'        </DataItem>'
    write(nUnitFile,'(A)')'    </Geometry>'

    ! Time series
    nflds = (ilast - ifirst +1)/SaveVisu  + 1
    write(nUnitFile,'(A)')'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    write(nUnitFile,'(A)')'        <Time TimeType="List">'
    write(nUnitFile,'(A,I6,A)')'        <DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
    write(nUnitFile,'(A)',advance='no')'        '
    do ifld = ifirst-1,ilast,SaveVisu
      write(nUnitFile,'(I10)',advance='no') ifld
    enddo
    write(nUnitFile,'(A)')'        </DataItem>'
    write(nUnitFile,'(A)')''; write(nUnitFile,'(A)')'       </Time>'

    ! attribute
    do  ifld=ifirst-1,ilast,SaveVisu
      write(nUnitFile,'(A,I10.10,A)')'        <Grid Name="T',ifld,'" GridType="Uniform">'
      write(nUnitFile,'(A)')'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(nUnitFile,'(A)')'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
      if(save_ux)    call Write_XDMF_One(nUnitFile,ifld,'ux')
      if(save_uy)    call Write_XDMF_One(nUnitFile,ifld,'uy')
      if(save_uz)    call Write_XDMF_One(nUnitFile,ifld,'uz')
      if(save_pr)    call Write_XDMF_One(nUnitFile,ifld,'pr')
      if(save_wx)    call Write_XDMF_One(nUnitFile,ifld,'wx')
      if(save_wy)    call Write_XDMF_One(nUnitFile,ifld,'wy')
      if(save_wz)    call Write_XDMF_One(nUnitFile,ifld,'wz')
      if(save_wMag)  call Write_XDMF_One(nUnitFile,ifld,'wMag')
      if(save_Q_vor) call Write_XDMF_One(nUnitFile,ifld,'Q' )
      if(save_lamda2)call Write_XDMF_One(nUnitFile,ifld,'lamda2')
      write(nUnitFile,'(A)')'        </Grid>'
    enddo

    write(nUnitFile,'(A)')'    </Grid>'
    write(nUnitFile,'(A)')'</Domain>'
    write(nUnitFile,'(A)')'</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)
  end subroutine InitVisu

  !******************************************************************
  ! Write_XDMF_One
  !******************************************************************
  subroutine Write_XDMF_One(nUnitFile, ifld,chAttribute)
    implicit none
    integer,intent(in)::nUnitFile,ifld
    character(*),intent(in)::chAttribute

    ! locals
    character(64)::chFile
    integer::iprec=mytype_save
    
    write(chFile,'(A,A,I10.10)')"VisuFor"//trim(RunName),"_"//trim(adjustl(chAttribute))//"_",ifld
    write(nUnitFile,'(A)')'            <Attribute Name="'//trim(chAttribute)//XDMF_SET_TYPE_STR
    write(nUnitFile,'(A,I1,A,3I5,A)')'                <DataItem Format="Binary" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nzc,nyc,nxc,'">'
    write(nUnitFile,'(A)')'                    '//trim(chFile)
    write(nUnitFile,'(A)')'                </DataItem>'
    write(nUnitFile,'(A)')'            </Attribute>'

  end subroutine Write_XDMF_One

  !******************************************************************
  ! dump_visu
  !******************************************************************
  subroutine dump_visu(ntime,ux,uy,uz,pressure,ArrTemp)
    implicit none
    integer,intent(in)::ntime
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),intent(inout)::ArrTemp

    ! locals
    character(24)::ch
    character(256)::chFile
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::dudy,dudz,dvdx,dvdz,dwdx,dwdy
    real(RK)::caj,cac1,cac2,cac12,vor_x,vor_y,vor_z

    ! begin to dump
    write(ch,'(I10.10)')ntime
 
    ! ux
    if(save_ux) then
      do kc=ystart(3),yend(3)
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
             ArrTemp(ic,jc,kc)=half*(ux(ic+1,jc,kc)+ux(ic,jc,kc))+uCRF
          enddo
        enddo
      enddo
      chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_ux_"//trim(adjustl(ch))
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! uy
    if(save_uy) then
      do kc=ystart(3),yend(3)
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
             ArrTemp(ic,jc,kc)=half*(uy(ic,jc+1,kc)+uy(ic,jc,kc))
          enddo
        enddo
      enddo
      chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_uy_"//trim(adjustl(ch))
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! uz
    if(save_uz) then
      do kc=ystart(3),yend(3)
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
             ArrTemp(ic,jc,kc)=half*(uz(ic,jc,kc+1)+uz(ic,jc,kc))
          enddo
        enddo
      enddo
      chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_uz_"//trim(adjustl(ch))
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! pressure
    if(save_pr) then
      chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_pr_"//trim(adjustl(ch))
      call decomp_2d_write_every(y_pencil,pressure(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)),iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wx
    if(save_wx) then
      do kc=ystart(3),yend(3)
        kp=kc+1
        km=kc-1
        do jc=ystart(2),yend(2)
          jp=jc+1
          jm=jc-1
          cac1 = rdyc(jc)
          cac2 = rdyc(jp)    
          cac12= cac1 - cac2
          do ic=ystart(1),yend(1)
            dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *quarter
            dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                  +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                  -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *quarter
            ArrTemp(ic,jc,kc)= dwdy -dvdz
          enddo
        enddo
      enddo
      chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_wx_"//trim(adjustl(ch))
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wy
    if(save_wy) then
      do kc=ystart(3),yend(3)
        kp=kc+1
        km=kc-1
        do jc=ystart(2),yend(2)
          do ic=ystart(1),yend(1)
            ip=ic+1
            im=ic-1
            dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *quarter
            dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *quarter
            ArrTemp(ic,jc,kc)= dudz -dwdx
          enddo
        enddo
      enddo
      chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_wy_"//trim(adjustl(ch))
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wz
    if(save_wz) then
      do kc=ystart(3),yend(3)
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
            dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2  &
                  +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12 &
                  -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1  )  *quarter
            dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *quarter
            ArrTemp(ic,jc,kc)= dvdx -dudy
          enddo
        enddo
      enddo
      chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_wz_"//trim(adjustl(ch))
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wMag
    if(save_wMag) then
      do kc=ystart(3),yend(3)
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

            dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2    &
                  +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12   &
                  -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1    )    *quarter
            dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *quarter
      
            dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *quarter
            dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *quarter

            dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *quarter
            dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                  +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                  -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *quarter

            vor_x= dwdy-dvdz
            vor_y= dudz-dwdx
            vor_z= dvdx-dudy
            ArrTemp(ic,jc,kc)= sqrt(vor_x*vor_x +vor_y*vor_y +vor_z*vor_z)
          enddo
        enddo
      enddo
      chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_wMag_"//trim(adjustl(ch)) 
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! Q
    if(save_Q_vor) then
      call Clc_Q_vor(ux,uy,uz,ArrTemp)
      chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_Q_"//trim(adjustl(ch)) 
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif  
  
   ! lamda2
   if(save_lamda2) then
     call Clc_lamda2(ux,uy,uz,ArrTemp)
     chFile = trim(Res_Dir)//"VisuFor"//trim(RunName)//"_lamda2_"//trim(adjustl(ch)) 
     call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
   endif

  end subroutine dump_visu

  !**********************************************************************
  ! Delete_Prev_Restart
  !**********************************************************************
  subroutine Delete_Prev_Restart(ntime)
    implicit none
    integer,intent(in):: ntime

    ! lcoals
    character(24)::ch
    character(256)::chFile

    if(nrank/=0) return
    write(ch,'(I10.10)')Prev_BackUp_itime   
    chFile = trim(RestartDir)//"RestartFor"//trim(RunName)//trim(adjustl(ch))
    call system("rm "//trim(adjustl(chFile))//" 2> /dev/null")
    
    Prev_BackUp_itime = ntime
  end subroutine Delete_Prev_Restart

  !******************************************************************
  ! write_restart
  !******************************************************************
  subroutine write_restart(ntime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
    implicit none
    integer,intent(in)::ntime
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(ysize(1),ysize(2),ysize(3)),intent(in):: HistXOld,HistYOld,HistZOld

    ! locals
    character(24)::ch
    character(256)::chFile
    integer::fh,code
    integer(kind=MPI_OFFSET_KIND)::disp,filesize

    ! begin to write restart file
    write(ch,'(I10.10)')ntime
    chFile = trim(RestartDir)//"RestartFor"//trim(RunName)//trim(adjustl(ch))
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, code)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,code)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND

    call decomp_2d_write_var(fh,disp,y_pencil,      ux(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))+uCRF)
    call decomp_2d_write_var(fh,disp,y_pencil,      uy(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,      uz(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,pressure(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    if(ischeme==FI_AB2 .and. WriteHistOld) then
      call decomp_2d_write_var(fh,disp,y_pencil,HistXOld)
      call decomp_2d_write_var(fh,disp,y_pencil,HistYOld)
      call decomp_2d_write_var(fh,disp,y_pencil,HistZOld)
    endif
    call MPI_FILE_CLOSE(fh,code)

  end subroutine write_restart

  !******************************************************************
  ! read_restart
  !******************************************************************
  subroutine read_restart(ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz,pressure
    real(RK),dimension(ysize(1),ysize(2),ysize(3)),intent(out):: HistXOld,HistYOld,HistZOld

    ! locals
    character(24)::ch
    character(256)::chFile
    integer::fh,code,ntime
    integer(kind=MPI_OFFSET_KIND)::disp,byte_total,filebyte

    ! begin to write restart file
    ntime= ifirst - 1
    write(ch,'(I10.10)')ntime
    chFile = trim(RestartDir)//"RestartFor"//trim(RunName)//trim(adjustl(ch))
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, code)
    if(code/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Read_Restart","Cannot open file: "//trim(chFile))

    if(ischeme==FI_AB2 .and. ReadHistOld) then
      byte_total=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*7_MPI_OFFSET_KIND
    else
      byte_total=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*4_MPI_OFFSET_KIND
    endif
    call MPI_FILE_GET_SIZE(fh,filebyte,code)
    if(filebyte /= byte_total .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"Read_Restart","file byte wrong")
    endif

    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_read_var(fh,disp,y_pencil,      ux(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)));ux=ux-uCRF
    call decomp_2d_read_var(fh,disp,y_pencil,      uy(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,      uz(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,pressure(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    if(ischeme==FI_AB2 .and. ReadHistOld) then
      call decomp_2d_read_var(fh,disp,y_pencil,HistXOld)
      call decomp_2d_read_var(fh,disp,y_pencil,HistYOld)
      call decomp_2d_read_var(fh,disp,y_pencil,HistZOld)
    endif
    call MPI_FILE_CLOSE(fh,code)

  end subroutine read_restart

end module m_IOAndVisu
