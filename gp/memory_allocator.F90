!
! $Id$
!
! memory allocation subsystem
!
! initially conceived to be used in the dirac and dalton packages
! as a replacement to the memget/memrel routines and related
!
! design/discussions    : ulf ekstrom, andre gomes, jetze sikkema
! inital implementation : andre gomes
!                         vu amsterdam, winter/spring 2008
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! memory allocator module
!
! allocates memory and optionally (work in progress) keeps track of memory usage etc  
!
! the idea is that it starts replacement of the memget/memrel setup
! for intrinsic types but eventually could be used to allocate/deallocate
! all sorts of derived types, if necessary.
!
! the information in the tracking might be useful to provide extended debugging/statistics
! of memory usage per subroutine and also at a higher level, like SCF or DFT.
! the idea is to integrate its use with the information from qenter/qexit so
! that one can always associate a given allocation to a subroutine, but also
! with the variable used,
!
      module memory_allocator

         use memory_tracker
                  
         implicit none
         private
         public alloc, dealloc, initData, allocator_init

         interface alloc
            module procedure alloc_real
            module procedure alloc_int
         end interface alloc

         interface dealloc
            module procedure dealloc_real
            module procedure dealloc_int
         end interface dealloc

         interface initData
            module procedure initDataReal
            module procedure initDataInt
         end interface initData


! error codes, variable kinds

         integer, parameter :: kchar = 1, kint = 4, kreal = 8, kcomplex = 8

         integer, parameter :: alreadyAllocated = 1, alreadyDeallocated = 2, &
                               allocateFailure  = 3, deallocateFailure  = 4

! actions for getting backtrace information
! have to find a better way to keep this in sync with the C part...
                               
         integer, parameter :: printAllToStdout = 1, printAllocCaller = 2, &
                               printPrettyTrace = 3, prinRawTrace     = 4
         
         contains

!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! error handling in allocation and other goodies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!            subroutine allocErrorHandler(status,size,data_kind,ignoreError)
            subroutine allocErrorHandler(status,size,data_kind)
               integer, parameter             :: charArrayLength = 132
               integer, optional              :: size
               integer, optional              :: data_kind
!               logical, optional              :: ignoreError
              logical          :: ignoreError

               integer                        :: status
               integer                        :: kbytesRequested
               logical                        :: isFatal = .True.
               character(len=charArrayLength) :: ErrorMessage, InfoMessage
               character(len=charArrayLength) :: trace_result
               character(len=charArrayLength) :: kbytesAvailable
               character(len=1)               :: dummy

!               if (.not.present(ignoreError)) then
!                  ignoreError = .False.
!               else
!                  write (*,*) ignoreError
!               endif
!               if (.not.present(size))        size = 0
!               if (.not.present(data_kind))   data_kind = 4
 
               ignoreError = .False.
                             
               InfoMessage  = 'Unknown error!'
               ErrorMessage = InfoMessage
               trace_result = 'unknown (stack trace info unavailable)'

               call get_traceback_info (printAllocCaller,trace_result,charArrayLength)
               
               select case (status)
                  case (alreadyAllocated)      
                     InfoMessage  = 'Trying to allocate variable already allocated. Please submit a bug report.'
                     ErrorMessage = 'Variable is already allocated!'
                     isFatal = .True.

                  case (alreadyDeallocated)     
                     InfoMessage  = 'Trying to deallocate variable already deallocated. Please submit a bug report.'
                     ErrorMessage = 'Trying to deallocate variable already deallocated'
                     isFatal = .False.
                  
                  case (allocateFailure)
                     call get_memory_info(kbytesAvailable,charArrayLength)
                     kbytesRequested = (size*data_kind)/1024

                     InfoMessage  = 'Error! Allocation failed!'
                     write (*,*) kbytesRequested
                     write (*,*) kbytesAvailable
                     write (InfoMessage,'(a,i12,a)') 'Attempted to allocate ',kbytesRequested, &
                                                     ' kB but only '//kbytesAvailable(1:24)//'available'

                     ErrorMessage = 'Memory allocation failure at '//trace_result
                     isFatal = .True.


                  case (deallocateFailure)
                     InfoMessage  = 'Memory deallocation failed at '//trace_result//' !'
                     ErrorMessage = 'deallocation failed!'
                     isFatal = .False.
   
               end select

               if (isFatal) then
                  if (ignoreError) then
                     write (*,*) ' '
                     write (*,*) 'Non-fatal error at '//trace_result//'. The program will continue.'
                     write (*,*) ' '
                  else
                     write (*,*) ' '
                     write (*,*) InfoMessage
                     write (*,*) ' '
! temporarily this traceback call has been enablen in quit_ so we don't need it now...
!                    call get_traceback_info (printPrettyTrace,dummy,1)
                     call quit(ErrorMessage)
                  endif
               else
                  write (*,*) ' '
                  write (*,*) 'Warning! non-fatal error at '//trace_result
                  write (*,*) ' '
                  write (*,*) InfoMessage
                  write (*,*) ' '
                  write (*,*) 'The program will continue...'
                  write (*,*) ' '
               endif
               
            end subroutine allocErrorHandler

            subroutine allocInitAllocated(data,size)
               integer  :: i, size
               real(kind=kreal), allocatable, target :: data(:)

               call initData(data,size)
!               call allocRegisterAllocation(kind(data),size)

            end subroutine allocInitAllocated

            subroutine quit2(message)
               character(len=*) :: message
 
               write (*,*) message
               stop
            end subroutine quit2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialization for data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            subroutine initDataReal(data,size)
               integer :: i,size
               real(kind=kreal) :: data(:)

               forall (i=1:size) data(i) = 0.0D0
            end subroutine initDataReal
  
            subroutine initDataInt(data,size)
               integer :: i, size
               integer :: data(:)

               forall (i=1:size) data(i) = 0
            end subroutine initDataInt
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! memory allocation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
            subroutine alloc_real(data,size)
               integer            :: size, status
               real(kind=kreal), allocatable, target :: data(:)
              
               if (.not.allocated(data)) then
                  allocate(data(size),stat=status) 

                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,size,kind(data))

                  else
                     call initData(data,size)
!                     call registerAllocation(data,size)

                  endif
               else
                  call allocErrorHandler(alreadyAllocated)

               endif

            end subroutine alloc_real

!            subroutine registerAllocationReal(data,size)
!               real(kind=kreal), allocatable, target :: data(:)
!               integer(kind=kint), pointer           :: p_mem
!               character(len=300)                    :: trace_result
               
!               call get_memory_address(data, p_mem)
!               call print_trace_allocfromfnct(trace_result)
!!               call appendAllocation(chunk,p_mem,trace_result)
               
!            end subroutine registerAllocationReal

            subroutine alloc_int(data,size)
               integer  :: status
               integer :: size
               integer, allocatable :: data(:)
               
               if (.not.allocated(data)) then
                  allocate(data(size),stat=status)

                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,size,kind(data))

                  else
                     call initData(data,size)

                  endif
               else
                  call allocErrorHandler(alreadyAllocated)

               endif

            end subroutine alloc_int


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! memory deallocation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

            subroutine dealloc_real(data)
               integer :: status
               real(kind=kreal), allocatable :: data(:)

               if (allocated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure)
                     
                  endif
                  
               else
                  call allocErrorHandler(alreadyDeallocated)
                  
               endif
               
            end subroutine dealloc_real


            subroutine dealloc_int(data)
               integer, allocatable :: data(:)

               deallocate(data)
            end subroutine dealloc_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interface to the tracker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine allocator_init
         call trackInitMemoryStack()
!         call trackCreateGroup('main')          ! always create a main tracking group
!         call trackInitGroup('main')
      end subroutine allocator_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! end of memory allocator module definition
!
      end module memory_allocator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
