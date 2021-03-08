module m_LogInfo
  use MPI
  use decomp_2d,only: nrank  
  implicit none
  private

  integer,parameter,public:: ErrT_NoError = 0
  integer,parameter,public:: ErrT_Abort   = 1
  integer,parameter,public:: ErrT_Pass    = 2
    
  character(512):: g_Ch_LastReportedError
  integer :: g_ET_LastReportedError
  character(2),dimension(10),parameter::bullet= (/">>"," >","++"," +","--"," -","**"," *","==" ," ="/)
  character(11),parameter,dimension(0:2):: chErrorT = (/"No error   " , "Fatal error" , "Warning    "/)
    
  type LogInfo
    character(256)::  chLogFileName = "file.log"
    integer:: rprt_lvl_file = 2
    integer:: rprt_lvl_cmdw = 3
    integer:: unit_file     = 100 
  contains
    procedure:: InitLog   => LI_InitLog
    procedure:: OpenFile  => LI_OpenFile
    procedure:: CloseFile => LI_CloseFile
    procedure:: LI_OutInfo
    procedure:: LI_OutInfo2
    generic:: OutInfo         => LI_OutInfo, LI_OutInfo2
    procedure:: CheckForError => LI_CheckForError
  end type LogInfo
  type(LogInfo),public :: MainLog
contains

  !********************************************************************************
  !   LI_InitLog
  !********************************************************************************
  subroutine LI_InitLog(this,nunit, Dir_Res, RunName , file_lvl, cmdw_lvl )
    implicit none
    class(LogInfo):: this
    integer, intent(in):: nunit
    character(*), intent(in):: Dir_Res, RunName
    integer,intent(in):: file_lvl, cmdw_lvl
    character(256):: chFile
    
    chFile = trim(Dir_Res)//trim(RunName)//".log"
    call this%OpenFile(nunit, chFile, RunName )
    this%rprt_lvl_file = file_lvl
    this%rprt_lvl_cmdw = cmdw_lvl
  end subroutine LI_InitLog

  !********************************************************************************
  !   LI_OpenFile
  !********************************************************************************
  subroutine LI_OpenFile(this,nUnit,chFile,RunName)
    implicit none
    class(LogInfo)::this
    integer,intent(in)::nUnit
    character(*),intent(in)::chFile,RunName

    ! locals
    integer:: myistat
    
    this%unit_file = nUnit
    this%chLogFileName = chFile
    if(nrank==0) then
      open(file = chFile, unit = nUnit, status='replace', IOSTAT=myistat)
      write(this%unit_file , *, IOSTAT=myistat ) "**************************************************************************"
      write(this%unit_file , *, IOSTAT=myistat ) "Log file for run :"// trim(RunName)
      write(this%unit_file , *, IOSTAT=myistat ) "**************************************************************************"
      write(this%unit_file , *, IOSTAT=myistat )
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,myistat)
    if(nrank/=0) open(file = chFile, unit = nUnit, status='old', IOSTAT=myistat)
  end subroutine 

  !********************************************************************************
  !   LI_CloseFile
  !********************************************************************************
  subroutine LI_CloseFile(this)
    implicit none
    class(LogInfo) this

    ! locals
    integer:: myistat
    close(this%unit_file, IOSTAT=myistat)
  end subroutine LI_CloseFile

  !********************************************************************************
  !   LI_OutInfo
  !********************************************************************************
  subroutine LI_OutInfo( this, chInfo, lvl , no_bull )
    implicit none
    class(LogInfo) this
    character(*),intent(in):: chInfo
    integer,intent(in):: lvl
    logical,optional,intent(in):: no_bull

    ! locals
    logical:: l_no
    integer:: myistat
    character(64) ch100, ch101
    
    l_no = .false.
    if( present(no_bull) ) l_no = no_bull
      if( lvl > 1 ) then    
        write(ch100,"(A,I3,A)")"(",2*(lvl-1), "x , A2, x ,A)" 
        write(ch101,"(A,I3,A)")"(",2*(lvl-1), "x ,A)"
      else
        ch100 = "(x , A2, x ,A)"
	      ch101 = "(x ,A)"
      endif
    
      if( lvl <= this%rprt_lvl_file )then
        if(lvl == 1 )write(this%unit_file, *, IOSTAT=myistat )
        if(l_no)then
          write(this%unit_file, ch101, IOSTAT=myistat ) trim(chInfo)    
        else
          write(this%unit_file, ch100, IOSTAT=myistat ) bullet(lvl), trim(chInfo)    
        end if
      end if
    
      if( lvl <= this%rprt_lvl_cmdw )then
        if(lvl == 1 )write(*,*)
        if(l_no)then
          write(*, ch101 ) trim(chInfo)
        else
          write(*, ch100 ) bullet(lvl), trim(chInfo)
        end if
      end if
  end subroutine  LI_OutInfo

  !********************************************************************************
  !   LI_OutInfo2
  !********************************************************************************
  subroutine LI_OutInfo2( this, chInfo , chInfo2, lvl , no_bull )
    implicit none
    class(LogInfo) this
    character(*), intent(in):: chInfo, chInfo2
    integer,intent(in):: lvl
    logical,optional,intent(in):: no_bull
    logical:: l_no
    
    l_no = .false.
    if( present(no_bull) ) l_no = no_bull
    call this%OutInfo(chInfo, lvl, l_no );
    call this%OutInfo(chInfo2, lvl, l_no );
  end subroutine LI_OutInfo2
            
  !********************************************************************************
  !   LI_CheckForError
  !     checking the input error message (Err_type) and creating a message in command  
  !     window and logfile, then stopping the execution of the program if necessary
  !********************************************************************************
  subroutine LI_CheckForError( this, Err_type, chMethod, chMessage )
    implicit none
    class(LogInfo):: this
    integer, intent(in):: Err_type
    character(*),intent(in):: chMethod, chMessage
    
    select case(Err_type)
    case(ErrT_NoError)
      ! no error occurred, the program will continue its normal execution
      g_ET_LastReportedError = Err_type
      g_ch_LastReportederror = "No error occurred in :"//trim(chMethod)//". Message: "//trim(chMessage)
      return
        
    case(ErrT_Abort)
      ! a severe error occurred and program should be aborted 
      call this%OutInfo( "A severe error occurred in program", 1 )
      call this%OutInfo( "Error occurred in :"//trim(chMethod), "Error message is:"//trim(chMessage) , 2 )
      g_ET_LastReportedError = Err_type
      g_ch_LastReportederror = "A severe error occurred in :"//trim(chMethod)//". Message: "//trim(chMessage)
      stop
        
    case(ErrT_Pass)
      ! a warning occurred in the program, a message will appear on the screen and 
      ! a message will be sent to log file but program continues running
        
      call this%OutInfo( "A warning is reported in program", 1 )
      call this%OutInfo( "A warning occurred in :"//trim(chMethod), "Warning message is:"//trim(chMessage) , 2 )
      g_ET_LastReportedError = Err_type
      g_ch_LastReportederror = "A warning occurred in :"//trim(chMethod)//". Message: "//trim(chMessage)
      return
    end select
  end subroutine LI_CheckForError
end module m_LogInfo
