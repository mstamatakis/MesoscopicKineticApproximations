module err_warn_module

use constants_module, only: iwrite, crestartfname, maxrepwarnings
use parser_module

#ifdef WITH_MPI_
  use mpi
#endif

implicit none

! Error-warning module - modus operandi:
! For large software applications it becomes tedious to have one error module with all the possible 
! error codes and messages. It is much better from a modularity standpoint to have the error codes 
! and descriptions in the module to which they are pertinent. A small problem to be resolved is however
! that the error codes may not be unique; it would be difficult to keep track of error codes in all other
! modules and make sure that they are not duplicated. A way of this is to impose a convention:
!
!       Decide how may errors could be there at maximum in a module; here we allow for 1000. 
!       Then assign an integer in the set 0, 1000, 2000, 3000, ... to each programming unit that uses 
!       the err_warn_module.
!
! This works as follows: each module has a private error subroutine which
!       - defines an integer parameter, the callerID, and a string, the callerName
!       - has a list of all the possible errors that can be issued by this module
!       - calls the issue_error subroutine of the err_warn_module, giving it the callerID and callerName,
!         error number (errorID), description, and more information.
! The issue_error subroutine issues the error with error-code equal to callerID+errorID; 
! the caller name is reported as well to avoid confusion (in fact we still have more than one module 
! use the same callerID to be compatible with the error numbers of older versions of Zacros). 
! Subsequently, the issue_error subroutine either stops execution or not, depending on the running mode
! which is set by the stop_on_error field.

! Warnings work in a similar way; however, we allow for suppressing the output of repeat-warnings, if they
! are issued more than maxrepwarnings (this constant is defined in contants_module). The way this works is 
! as follows:
!       - every time a warning is issued, the callerID, callerName and warnID are logged into fields 
!         warngCallerIDs, warngCallerNames and warngIDs, respectively
!       - in field warngCounters, we also keep track of the number of times the warning is issued
!       - if warngCounters is equal to maxrepwarnings, we append a note in the last warning to be printed
!         (saying that further output will be supressed)
!       - if warngCounters exceeds maxrepwarnings, no further output is generated, but still the warnings 
!         are counted.
! We have introduced memory amortization to the aforementioned fields (warngCallerIDs, warngCallerNames, 
! warngIDs, warngCounters) to avoid problems with filling them up, or overestimating how much memory is needed
! for them in the first place. For this reason, it doesn't really matter what the value of n0 is upon initialization
! (see the arguments of err_warn_initialize).

integer, parameter :: callerNameLen = 64 ! string length of caller names

type err_warn
  
!  logical :: stop_on_error = .true.

  integer :: nWarngs
  integer, allocatable :: warngCallerIDs(:)
  character(len=callerNameLen), allocatable :: warngCallerNames(:) ! This might seem superfluous, but it's there for safety, 
                                                                   ! since the same callerID might be used for more than one modules.
  integer, allocatable :: warngIDs(:)
  integer, allocatable :: warngCounters(:)
  
  contains
    final :: cleanup_err_warn
    procedure :: initialize => err_warn_initialize
    procedure :: error => issue_error
    procedure :: warning => issue_warning
    procedure :: summarize_warnings
    procedure :: index_of_warngID
    procedure :: addNewWarngToLog
    procedure :: updateWarngInLog
    procedure :: print => print_err_warn_structure
  
end type err_warn

! Define an allocatable error-warning reporter object
type(err_warn), allocatable :: errwarn_globj

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine cleanup_err_warn(this)

    implicit none

    type(err_warn), intent(inout) :: this
	
    ! Cleanup err_warn_module allocatable fields
    if(allocated(this%warngCallerIDs)) deallocate(this%warngCallerIDs)
    if(allocated(this%warngCallerNames)) deallocate(this%warngCallerNames)
    if(allocated(this%warngIDs)) deallocate(this%warngIDs)
    if(allocated(this%warngCounters)) deallocate(this%warngCounters)

    this%nWarngs = 0

	return

  end subroutine cleanup_err_warn
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine err_warn_initialize(this,n0)

    implicit none

    class(err_warn), intent(inout) :: this
	
	integer n0

    this%nWarngs = 0

    allocate(this%warngCallerIDs(n0))
	this%warngCallerIDs = -1
    allocate(this%warngCallerNames(n0))
	this%warngCallerNames = ''
    allocate(this%warngIDs(n0))
	this%warngIDs = 0
    allocate(this%warngCounters(n0))
	this%warngCounters = 0

	return

  end subroutine err_warn_initialize

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine issue_error(this, callerID, callerName, errorID, &
                         writeunit, errormessage, moreinfo)

    implicit none

    class(err_warn), intent(inout) :: this
    logical :: stop_on_error = .true.
    integer, intent(in) :: callerID, errorID
    integer, intent(in), optional :: writeunit
    integer writeunit_local
    integer nwords, iw(6), ierror, errorcode
    character(len=callerNameLen) :: callerName
	character(len=*), intent(in) :: errormessage
	character(len=*), intent(in), optional :: moreinfo
   
    if (present(writeunit)) then
       writeunit_local = writeunit
    else
       writeunit_local = iwrite
    endif
	
	write(writeunit_local,'(/,a)') '***************'

	if (errorID == -1) then
		write(writeunit_local,'(/,a)') 'Low level internal error from ' // &
			trim(callerName) // ': attempted to issue an error message ' // &
			'with unknown error code ' // trim(errormessage) // & 
            ! errormessage holds the unknown error code, since errorID was replaced by -1
			'. Please notify the developers about this...'

	elseif (errorID == -10) then
		write(writeunit_local,'(/,a)') 'General error code ' // trim(int2str(-errorID)) // ' from ' // &
			trim(callerName) // ': at least one non-optional input file is missing.'
		write(writeunit_local,'(/,a)') 'More information: '
		write(writeunit_local,'(a)') 'Input file ' // trim(moreinfo) // ' could not be found.'

	elseif (errorID == -11) then
		write(writeunit_local,'(/,a)') 'General error code ' // trim(int2str(-errorID)) // ' from ' // &
			trim(callerName) // ': error reading from restart file ' // trim(crestartfname) // '.'
		write(writeunit_local,'(a)') 'Unable to restart the simulation.'
		write(writeunit_local,'(/,a)') 'More information: '
		write(writeunit_local,'(a)') 'File was last positioned at byte ' // trim(moreinfo) // '.'

	elseif (errorID == -12) then
		write(writeunit_local,'(/,a)') 'General error code ' // trim(int2str(-errorID)) // ' from ' // &
			trim(callerName) // ': error writing to restart file ' // trim(crestartfname) // '.'
		write(writeunit_local,'(a)') 'Simulation will probably be impossible to restart.'
		write(writeunit_local,'(/,a)') 'More information: '
		write(writeunit_local,'(a)') 'File was last positioned at byte ' // trim(moreinfo) // '.'

	elseif (errorID == -100) then
		write(writeunit_local,'(/,a)') 'General error code ' // trim(int2str(-errorID)) // ' from ' // &
			trim(callerName) // ': invalid duplicate specification.'
		write(writeunit_local,'(/,a)') 'More information: '
        call break_words(moreinfo,';','#',nwords,iw)
		write(writeunit_local,'(a)') 'Keyword ' // moreinfo(iw(1):iw(2)) // ' has already been parsed; ' // &
                               'yet, it ' // 'appears again in line ' // moreinfo(iw(3):iw(4)) // &
                               ' of ' // moreinfo(iw(5):iw(6)) // '.'

	else

		if (callerID >= 800000) then
			write(writeunit_local,'(/,a)',advance='no') 'Internal error code '
		else
			write(writeunit_local,'(/,a)',advance='no') 'Error code '
		endif
		
		write(writeunit_local,'(a)') trim(int2str(callerID+errorID)) // ' from ' // &
			trim(callerName) // ': ' // trim(errormessage)
		
                if (present(moreinfo)) then
			write(writeunit_local,'(/,a)') 'More information: '
			write(writeunit_local,'(a)') trim(moreinfo)
		endif
	
	endif
	
    write(writeunit_local,'(/,a)') '***************'
    
    if (stop_on_error) then
		write(writeunit_local,'(/,a)') '> ABNORMAL TERMINATION DUE TO FATAL ERROR <'
		flush(writeunit_local)
#ifdef WITH_MPI_
        call MPI_Abort(MPI_COMM_WORLD, errorcode, ierror)
#endif
		stop
	else
		write(writeunit_local,'(/,a)') '> TESTING MODE - STOPPING DUE TO FATAL ERROR DEFERRED <'
	endif

	return

  end subroutine issue_error

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine issue_warning(this, callerID, callerName, warngID, &
                         writeunit, warngmessage, moreinfo)

    implicit none

    class(err_warn), intent(inout) :: this

    integer, intent(in) :: callerID, warngID
    integer, intent(in), optional :: writeunit
    integer writeunit_local, i
    character(len=callerNameLen) :: callerName
	character(len=*), intent(in) :: warngmessage
	character(len=*), intent(in), optional :: moreinfo

    i = this%index_of_warngID(callerID,callerName,warngID)
    if (i < 0) then
        i = this%nWarngs + 1
        call this%addNewWarngToLog(callerID,callerName,warngID)
    else
        call this%updateWarngInLog(i)
    end if

    if (this%warngCounters(i) > maxrepwarnings) return

    if (present(writeunit)) then
        writeunit_local = writeunit
    else
        writeunit_local = iwrite
    endif

	write(writeunit_local,'(/,a)') '***************'

	if (warngID == -1) then
		write(writeunit_local,'(/,a)') 'Low level internal warning from ' // &
			trim(callerName) // ': attempted to issue a warning message ' // &
			'with unknown warning code ' // trim(warngmessage) // &
            ! warngmessage holds the unknown warning code, since warngID was replaced by -1
			'. Please notify the developers about this...'
    else
        write(writeunit_local,'(/,a)') 'Warning code ' // trim(int2str(callerID+warngID)) // ' from ' // &
            trim(callerName) // ': ' // trim(warngmessage)
    endif

	if (present(moreinfo)) then
		write(writeunit_local,'(/,a)') 'More information: '
		write(writeunit_local,'(a)') trim(moreinfo)
	endif

    if (this%warngCounters(i) == maxrepwarnings) then
        write(writeunit_local,'(/,a)') 'NOTE: Repeat-warnings of this kind will be suppressed from producing further output...'
    end if

    write(writeunit_local,'(/,a)') '***************'

    flush(writeunit_local)
    
	return

  end subroutine issue_warning

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  integer function index_of_warngID(this,callerID,callerName,warngID)

    use parser_module, only: striccompare

    class(err_warn), intent(inout) :: this
    integer, intent(in) :: callerID, warngID
    character(len=callerNameLen) :: callerName

    integer i

    index_of_warngID = -1

    do i = 1,this%nWarngs
        if (this%warngCallerIDs(i) == callerID) then
            if (striccompare(trim(this%warngCallerNames(i)),trim(callerName)) .and. &
                this%warngIDs(i) == warngID) then
                index_of_warngID = i
                return
            endif
        endif
	enddo

	return

  end function index_of_warngID

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine addNewWarngToLog(this,callerID,callerName,warngID)

    class(err_warn), intent(inout) :: this
    integer, intent(in) :: callerID, warngID
    character(len=callerNameLen) :: callerName

    integer, allocatable :: tmpInt(:)
    character(len=callerNameLen), allocatable :: tmpChar(:)

    integer i

    ! Memory amortization: check size and increase if necessary
    if (this%nWarngs == size(this%warngCallerIDs)) then

        ! Take care of integer fields first
        allocate(tmpInt(this%nWarngs))

        ! copy over warngCallerIDs
        do i = 1,this%nWarngs
            tmpInt(i) = this%warngCallerIDs(i)
        end do
        deallocate(this%warngCallerIDs)
        allocate(this%warngCallerIDs(2*this%nWarngs))
        this%warngCallerIDs = -1
        do i = 1,this%nWarngs
            this%warngCallerIDs(i) = tmpInt(i)
        end do

        ! copy over warngIDs
        do i = 1,this%nWarngs
            tmpInt(i) = this%warngIDs(i)
        end do
        deallocate(this%warngIDs)
        allocate(this%warngIDs(2*this%nWarngs))
        this%warngIDs = 0
        do i = 1,this%nWarngs
            this%warngIDs(i) = tmpInt(i)
        end do

        ! copy over warngCounters
        do i = 1,this%nWarngs
            tmpInt(i) = this%warngCounters(i)
        end do
        deallocate(this%warngCounters)
        allocate(this%warngCounters(2*this%nWarngs))
        this%warngCounters = 0
        do i = 1,this%nWarngs
            this%warngCounters(i) = tmpInt(i)
        end do

        deallocate(tmpInt)

        ! Then take care of the character field
        allocate(tmpChar(this%nWarngs))

        ! copy over warngCallerNames
        do i = 1,this%nWarngs
            tmpChar(i) = this%warngCallerNames(i)
        end do
        deallocate(this%warngCallerNames)
        allocate(this%warngCallerNames(2*this%nWarngs))
        this%warngCallerNames = ''
        do i = 1,this%nWarngs
            this%warngCallerNames(i) = tmpChar(i)
        end do

        deallocate(tmpChar)

    end if

    this%nWarngs = this%nWarngs + 1
    this%warngCallerIDs(this%nWarngs) = callerID
    this%warngCallerNames(this%nWarngs) = callerName
    this%warngIDs(this%nWarngs) = warngID
    this%warngCounters(this%nWarngs) = 1

	return

  end subroutine addNewWarngToLog

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine updateWarngInLog(this,i)

    class(err_warn), intent(inout) :: this
    integer, intent(in) :: i

    this%warngCounters(i) = this%warngCounters(i) + 1

    return

  end subroutine updateWarngInLog

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine print_err_warn_structure(this,writeunit)

    class(err_warn), intent(inout) :: this
    integer, optional, intent(in) :: writeunit

    integer i, writeunit_local

    if (present(writeunit)) then
        writeunit_local = writeunit
    else
        writeunit_local = iwrite
    endif

    write(writeunit_local,*) '--------------------------------------'
    write(writeunit_local,*) 'err_warn_structure:'

    do i = 1,this%nWarngs

        write(writeunit_local,'(i4,1x,i4,1x,a15,1x,i6,1x,i4)') i, &
                    this%warngCallerIDs(i), &
                    this%warngCallerNames(i), &
                    this%warngIDs(i), this%warngCounters(i)

    end do

    write(writeunit_local,*) '--------------------------------------'

    return

end subroutine print_err_warn_structure

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine summarize_warnings(this,writeunit)

    class(err_warn), intent(inout) :: this
    integer, optional, intent(in) :: writeunit

    integer i, iw, writeunit_local

    character(100) tmpc

    if (present(writeunit)) then
        writeunit_local = writeunit
    else
        writeunit_local = iwrite
    endif

    write(writeunit_local,'(/,a)') 'Summary of warnings:'
    write(writeunit_local,'(a,/)') '~~~~~~~~~~~~~~~~~~~~'

    do i = 1,this%nWarngs
        if (this%warngCounters(i) > maxrepwarnings) then
            iw = this%warngCounters(i) - maxrepwarnings
            tmpc = 's, out of which ' // trim(int2str(iw)) // ' times in silent mode.'
        elseif (this%warngCounters(i) == 1) then
            tmpc = '.'
        else ! this%warngCounters(i) > 1 .and. this%warngCounters(i) <= maxrepwarnings 
            tmpc = 's.'
        end if
        iw = this%warngCallerIDs(i)+this%warngIDs(i)
        write(writeunit_local,'(a)') 'Warning ' // trim(int2str(iw)) // &
            ' from '// trim(this%warngCallerNames(i)) // ' was triggered ' // &
            trim(int2str(this%warngCounters(i))) // ' time' // tmpc
    end do

end subroutine summarize_warnings

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module err_warn_module
