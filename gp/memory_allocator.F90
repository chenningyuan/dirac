!*
!*
!* Copyright (c) 2008 Andre Severo Pereira Gomes <gomes@few.vu.nl>
!* All rights reserved.
!*
!* Redistribution and use in source and binary forms, with or without
!* modification, are permitted provided that the following conditions
!* are met:
!* 1. Redistributions of source code must retain the above copyright
!*    notice, this list of conditions and the following disclaimer.
!* 2. Redistributions in binary form must reproduce the above copyright
!*    notice, this list of conditions and the following disclaimer in the
!*    documentation and/or other materials provided with the distribution.
!*
!* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
!* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
!* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
!* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
!* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
!* SUCH DAMAGE.
!*
!*
!
! $Id: memory_allocator.F90 7009 2008-06-07 23:10:40Z gomes $
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
         public alloc, dealloc, allocp, deallocp, initData, alloc_maxsize
         public allocator_init, allocator_setmax

         interface alloc
            module procedure alloc_real_1d
            module procedure alloc_real_2d
            module procedure alloc_real_3d
            module procedure alloc_real_4d
            module procedure alloc_real_5d
            module procedure alloc_real_6d
            module procedure alloc_int
         end interface alloc

         interface allocp
            module procedure alloc_preal
            module procedure alloc_pint
         end interface allocp

         interface dealloc
            module procedure dealloc_real_1d
            module procedure dealloc_real_2d
            module procedure dealloc_real_3d
            module procedure dealloc_real_4d
            module procedure dealloc_real_5d
            module procedure dealloc_real_6d
            module procedure dealloc_int
         end interface dealloc

         interface deallocp
            module procedure dealloc_preal
            module procedure dealloc_pint_0d
            module procedure dealloc_pint_1d
         end interface deallocp

         interface initData
            module procedure initData_real_1d
            module procedure initData_real_2d
            module procedure initData_real_3d
            module procedure initData_real_4d
            module procedure initData_real_5d
            module procedure initData_real_6d
            module procedure initDataInt
         end interface initData

         interface alloc_maxsize
            module procedure alloc_maxsize_k8
         end interface alloc_maxsize 


! error codes, variable kinds

         integer, parameter :: charArrayLength = 132
         integer, parameter :: kreal = 8, kcomplex = 8

         integer, parameter :: alreadyAllocated = 1, alreadyDeallocated = 2, &
                               allocateFailure  = 3, deallocateFailure  = 4

! actions for getting backtrace information
! have to find a better way to keep this in sync with the C part...
                               
         integer, parameter :: printAllToStdout = 1, printAllocCaller = 2, &
                               printPrettyTrace = 3, prinRawTrace     = 4
         
         integer :: max_mem = -1, curr_mem = 0

         contains

!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! error handling in allocation and other goodies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            subroutine allocErrorHandler(status,size,data_kind,var_id,ignoreError)
               integer                        :: status
! optional arguments
               logical, optional              :: ignoreError
               integer, optional              :: size
               integer, optional              :: data_kind
               character(len=*), optional     :: var_id 
! local variables
               logical                        :: ignoreErrorLocal, isFatal
               integer                        :: kbytesRequested, sizeLocal, kindLocal
               character(len=charArrayLength) :: ErrorMessage, InfoMessage
               character(len=charArrayLength) :: trace_result
               character(len=charArrayLength) :: kbytesAvailable
               character(len=1)               :: dummy

               isFatal = .True.

               if (present(ignoreError)) then 
                  ignoreErrorLocal = ignoreError
               else
                  ignoreErrorLocal = .False.
               endif

               if (present(size)) then 
                  sizeLocal = size
               else
                  sizeLocal = 0
               endif

               if (present(data_kind)) then
                  kindLocal = data_kind
               else
                  kindLocal = 4
               endif

               if (present(var_id)) then
                  if (var_id.ne.'None') then
                     write (*,*) ' '
                     write (*,*) ' *** Warning: at allocation for identifier '//var_id(1:42)
!                    write (*,*) ' '
                  endif
               endif
 
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
                     kbytesRequested = (sizeLocal*kindLocal)/1024

                     InfoMessage  = 'Error! Allocation failed!'
!                    write (*,*) kbytesRequested
!                    write (*,*) kbytesAvailable
                     write (InfoMessage,'(a,i9,a)')                       &
                        'Attempted to allocate',kbytesRequested,' Kb '//  &
                        'but memory available is '//kbytesAvailable(1:24)

!                    call catfile()
                     ErrorMessage = 'Memory allocation failure at '//trace_result(1:42)
                     isFatal = .True.


                  case (deallocateFailure)
                     InfoMessage  = 'Memory deallocation failed at '//trace_result(1:42)//' !'
                     ErrorMessage = 'deallocation failed!'
                     isFatal = .False.
   
               end select

               if (isFatal) then
                  if (ignoreErrorLocal) then
!                    write (*,*) ' '
                     write (*,*) 'Non-fatal error at '//trace_result(1:42)//'. The program will continue.'
!                    write (*,*) ' '
                  else
!                    write (*,*) ' '
                     write (*,*) InfoMessage
!                    write (*,*) ' '
                     call quit(ErrorMessage)
                  endif
               else
!                 write (*,*) ' '
                  write (*,*) 'Warning! non-fatal error at '//trace_result(1:42)
                  write (*,*) ' '
                  write (*,*) InfoMessage
                  write (*,*) ' '
                  write (*,*) 'The program will continue...'
!                 write (*,*) ' '
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

            subroutine initData_real_1di8(data,size)
               integer(kind=8) :: i,size
               real(kind=kreal) :: data(:)

               forall (i=1:size) data(i) = 0.0D0
            end subroutine initData_real_1di8

            subroutine initData_real_1d(data,size)
               integer :: i,size
               real(kind=kreal) :: data(:)

               forall (i=1:size) data(i) = 0.0D0
            end subroutine initData_real_1d

            subroutine initData_real_2d(data,size1,size2)
               integer :: i,j,size1,size2
               real(kind=kreal) :: data(:,:)

               forall (i=1:size1, j=1:size2) data(i,j) = 0.0D0
            end subroutine initData_real_2d

            subroutine initData_real_3d(data,size1,size2,size3)
               integer :: i,j,k,size1,size2,size3
               real(kind=kreal) :: data(:,:,:)

               forall (i=1:size1, j=1:size2, k=1:size3) data(i,j,k) = 0.0D0
            end subroutine initData_real_3d

            subroutine initData_real_4d(data,size1,size2,size3,size4)
               integer :: i,j,k,l,size1,size2,size3,size4
               real(kind=kreal) :: data(:,:,:,:)

               forall (i=1:size1, j=1:size2, k=1:size3, l=1:size4) data(i,j,k,l) = 0.0D0
            end subroutine initData_real_4d

            subroutine initData_real_5d(data,size1,size2,size3,size4,size5)
               integer :: i,j,k,l,m,size1,size2,size3,size4,size5
               real(kind=kreal) :: data(:,:,:,:,:)

               forall (i=1:size1, j=1:size2, k=1:size3, l=1:size4, m=1:size5) data(i,j,k,l,m) = 0.0D0
            end subroutine initData_real_5d

            subroutine initData_real_6d(data,size1,size2,size3,size4,size5,size6)
               integer :: i,j,k,l,m,n,size1,size2,size3,size4,size5,size6
               real(kind=kreal) :: data(:,:,:,:,:,:)

               forall (i=1:size1, j=1:size2, k=1:size3, l=1:size4, m=1:size5, n=1:size6) data(i,j,k,l,m,n) = 0.0D0
            end subroutine initData_real_6d
  
            subroutine initDataInt(data,size)
               integer :: i, size
               integer :: data(:)

               forall (i=1:size) data(i) = 0
            end subroutine initDataInt
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! memory allocation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!            subroutine registerAllocationReal(data,size)
!               real(kind=kreal), allocatable, target :: data(:)
!               integer, pointer           :: p_mem
!               character(len=300)                    :: trace_result

!               call get_memory_address(data, p_mem)
!               call print_trace_allocfromfnct(trace_result)
!               call appendAllocation(chunk,p_mem,trace_result)

!            end subroutine registerAllocationReal


!
! using allocatables
!
            subroutine alloc_real_1d(data,size,id)
               integer                               :: size, status
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), allocatable, target :: data(:)
              
               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None' 
               endif

               if (.not.allocated(data)) then
                  allocate(data(size),stat=status) 

                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,size,kind(data),var_id=local_id)
                  else
                     call initData(data,size)
!                    call registerAllocation(data,size)
                  endif
               else
                  call allocErrorHandler(alreadyAllocated,var_id=local_id)
               endif
            end subroutine alloc_real_1d

            subroutine alloc_real_2d(data,size1,size2,id)
               integer                               :: size1,size2, status, totalsize
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), allocatable, target :: data(:,:)
              
               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None' 
               endif

               if (.not.allocated(data)) then
                  allocate(data(size1,size2),stat=status) 

                  totalsize = size1*size2
                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,totalsize,kind(data),var_id=local_id)
                  else
!                    call initData(data,totalsize)
                     call initData(data,size1,size2)
!                    call registerAllocation(data,size)
                  endif
               else
                  call allocErrorHandler(alreadyAllocated,var_id=local_id)
               endif
            end subroutine alloc_real_2d

            subroutine alloc_real_3d(data,size1,size2,size3,id)
               integer                               :: size1,size2, size3, status, totalsize
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), allocatable, target :: data(:,:,:)
              
               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None' 
               endif

               if (.not.allocated(data)) then
                  allocate(data(size1,size2,size3),stat=status) 

                  totalsize = size1*size2*size3
                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,totalsize,kind(data),var_id=local_id)
                  else
!                    call initData(data,totalsize)
                     call initData(data,size1,size2,size3)
                  endif
               else
                  call allocErrorHandler(alreadyAllocated,var_id=local_id)
               endif
            end subroutine alloc_real_3d

            subroutine alloc_real_4d(data,size1,size2,size3,size4,id)
               integer                               :: size1,size2, size3, size4, status, totalsize
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), allocatable, target :: data(:,:,:,:)
              
               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None' 
               endif

               if (.not.allocated(data)) then
                  allocate(data(size1,size2,size3,size4),stat=status) 

                  totalsize = size1*size2*size3*size4
                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,totalsize,kind(data),var_id=local_id)
                  else
!                    call initData(data,totalsize)
                     call initData(data,size1,size2,size3,size4)
                  endif
               else
                  call allocErrorHandler(alreadyAllocated,var_id=local_id)
               endif
            end subroutine alloc_real_4d

            subroutine alloc_real_5d(data,size1,size2,size3,size4,size5,id)
               integer                               :: size1,size2, size3, size4,size5, status, totalsize
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), allocatable, target :: data(:,:,:,:,:)
              
               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None' 
               endif

               if (.not.allocated(data)) then
                  allocate(data(size1,size2,size3,size4,size5),stat=status) 

                  totalsize = size1*size2*size3*size4*size5
                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,totalsize,kind(data),var_id=local_id)
                  else
!                    call initData(data,totalsize)
                     call initData(data,size1,size2,size3,size4,size5)
                  endif
               else
                  call allocErrorHandler(alreadyAllocated,var_id=local_id)
               endif
            end subroutine alloc_real_5d

            subroutine alloc_real_6d(data,size1,size2,size3,size4,size5,size6,id)
               integer                               :: size1,size2,size3,size4,size5,size6
               integer                               :: status, totalsize
               character(len=*), optional            :: id
               character(len=charArrayLength)        :: local_id
               real(kind=kreal), allocatable, target :: data(:,:,:,:,:,:)
              
               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None' 
               endif

               if (.not.allocated(data)) then
                  allocate(data(size1,size2,size3,size4,size5,size6),stat=status) 

                  totalsize = size1*size2*size3*size4*size5*size6
                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,totalsize,kind(data),var_id=local_id)
                  else
!                    call initData(data,totalsize)
                     call initData(data,size1,size2,size3,size4,size5,size6)
                  endif
               else
                  call allocErrorHandler(alreadyAllocated,var_id=local_id)
               endif
            end subroutine alloc_real_6d


            subroutine alloc_int(data,size,id)
               integer                        :: size, status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               integer, allocatable, target   :: data(:)
               
               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (.not.allocated(data)) then
                  allocate(data(size),stat=status)

                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,size,kind(data),var_id=local_id)
                  else
                     call initData(data,size)
                  endif
               else
                  call allocErrorHandler(alreadyAllocated,var_id=local_id)
               endif
            end subroutine alloc_int

!
! using pointers
!

            subroutine alloc_preal(data,size,id)
               integer                        :: size, status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               real(kind=kreal), pointer      :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (.not.associated(data)) then
                  allocate(data(size),stat=status)

                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,size,kind(data),var_id=local_id)
                  else
                     call initData(data,size)
                  endif
               else
                  call allocErrorHandler(alreadyAllocated,var_id=local_id)
               endif
            end subroutine alloc_preal

            subroutine alloc_pint(data,size,id)
               integer                        :: size, status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               integer, pointer               :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (.not.associated(data)) then
                  allocate(data(size),stat=status)

                  if ( status.ne.0 ) then
                     call allocErrorHandler(allocateFailure,size,kind(data),var_id=local_id)
                  else
                     call initData(data,size)
                  endif
               else
                  call allocErrorHandler(alreadyAllocated,var_id=local_id)
               endif
            end subroutine alloc_pint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! memory deallocation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

            subroutine dealloc_real_1d(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               real(kind=kreal), allocatable  :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (allocated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_real_1d

            subroutine dealloc_real_2d(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               real(kind=kreal), allocatable  :: data(:,:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (allocated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_real_2d

            subroutine dealloc_real_3d(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               real(kind=kreal), allocatable  :: data(:,:,:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (allocated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_real_3d

            subroutine dealloc_real_4d(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               real(kind=kreal), allocatable  :: data(:,:,:,:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (allocated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_real_4d

            subroutine dealloc_real_5d(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               real(kind=kreal), allocatable  :: data(:,:,:,:,:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (allocated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_real_5d

            subroutine dealloc_real_6d(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               real(kind=kreal), allocatable  :: data(:,:,:,:,:,:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (allocated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_real_6d


            subroutine dealloc_preal(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               real(kind=kreal), pointer      :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (associated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_preal

            subroutine dealloc_int(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               integer, allocatable :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (allocated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_int

            subroutine dealloc_pint_1d(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               integer, pointer               :: data(:)

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (associated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_pint_1d

            subroutine dealloc_pint_0d(data, id)
               integer                        :: status
               character(len=*), optional     :: id
               character(len=charArrayLength) :: local_id
               integer, pointer               :: data

               if (present(id)) then
                  write (local_id,'(a)') id
               else
                  write (local_id,'(a)') 'None'
               endif

               if (associated(data)) then
                  deallocate(data,stat=status)
                  
                  if ( status.ne.0 ) then
                     call allocErrorHandler(deallocateFailure,var_id=local_id)
                  endif
               else
                  call allocErrorHandler(alreadyDeallocated,var_id=local_id)
               endif
            end subroutine dealloc_pint_0d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interface to the tracker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine allocator_init
         call trackInitMemoryStack()
!         call trackCreateGroup('main')          ! always create a main tracking group
!         call trackInitGroup('main')
      end subroutine allocator_init

      subroutine allocator_setmax(maxmemsize)
         integer :: maxmemsize

         if (maxmemsize.lt.0) then ! there are no limits so we use a large integer
            max_mem = huge(maxmemsize)/32
            maxmemsize = max_mem
         else
            max_mem = maxmemsize
         endif
         
      end subroutine allocator_setmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine alloc_maxsize_k8(maxsize,minsize,verbose)
         integer :: step, i, nstep, maximum, status, minsize
         integer :: maxsize, verbose, minimum, factor
         real(kind=8), allocatable :: data(:)

         factor  = 4096 

         call get_maximum_mwords(maximum,verbose)

         minimum = minsize 

         if      ( (max_mem.gt.0).and.(maximum.gt.0) ) then  ! if  we have set max_mem, we use it as a guess for
            maximum = min(max_mem,maximum*1024*1024)         ! the maximum or have as the ceiling for allocation
         else if ( (max_mem.le.0).and.(maximum.gt.0) ) then     
            maximum = maximum*1024*1024                      ! if max_mem is not set but maximum is, we use the latter
         else if ( (max_mem.gt.0).and.(maximum.le.0) ) then     
            maximum = max_mem                                ! if max_mem is set and maximum not, we use the former
         else 
            maximum = huge(maximum)/32                       ! if both max_mem and maximum are not valid, we use a large integer
         endif


!        write (*,*) 'debug, maximum(a):',maximum
         if (mod(maximum,factor).ne.0) then
            maximum = maximum+mod(maximum,factor) 
         endif
!        write (*,*) 'debug, maximum(b):',maximum

!        write (*,*) 'debug, minimum(a):',minimum
         if (mod(minimum,factor).ne.0) then
            minimum = minimum+mod(minimum,factor) 
         endif
!        write (*,*) 'debug, minimum(b):',minimum

         if (minimum.gt.0) nstep = maximum / minimum

         if (maximum.ge.minimum) then

            if (verbose.gt.0) then
               write (*,*) 'minimal 8-byte words required   :',minimum
               write (*,*) 'number of steps attepmted       :',nstep
               write (*,*) 'maximum 8-byte words start      :',maximum
            endif

            do i=1,nstep
               allocate(data(maximum),stat=status)
               if (status.eq.0) then
                  call initData(data,maximum)
                  if (verbose.gt.0) write (*,*) 'at try',i,': allocation ok for ',maximum,'8-byte words'
                  maxsize = maximum 
                  deallocate(data)
                  return
               else
                  maximum = maximum - maximum/(nstep-i+1) 
               endif
            enddo
         else
            call quit('Minimum space required unavaiable. Exiting.')
         endif

       end subroutine alloc_maxsize_k8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! end of memory allocator module definition
!
      end module memory_allocator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
