!*
!*
!* Copyright (c) 2008-2011 Andre Severo Pereira Gomes <andre.gomes@univ-lille1.fr>
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
! $Id: memory_allocator.F90 9830 2009-10-30 04:37:36Z gomes $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! memory allocator module
!
! allocates memory and (optionally) keeps track of memory usage 
!
      module allocator_eh
      
         use allocator_parameters
         use allocator_cfg
         use allocator_track_if

         implicit none
         
         private 

         public allocator_errorHandler

         public varIsAliveP
         public varIsAliveA

         interface varIsAliveP
            module procedure rPointerIsAlive_1d
            module procedure iPointerIsAlive_1d
            module procedure rPointerIsAlive_2d
            module procedure iPointerIsAlive_2d
            module procedure rPointerIsAlive_3d
            module procedure iPointerIsAlive_3d
            module procedure rPointerIsAlive_4d
            module procedure iPointerIsAlive_4d
         end interface varIsAliveP

         interface varIsAliveA
            module procedure rAllocatableIsAlive_1d
            module procedure iAllocatableIsAlive_1d
            module procedure rAllocatableIsAlive_2d
            module procedure iAllocatableIsAlive_2d
            module procedure rAllocatableIsAlive_3d
            module procedure iAllocatableIsAlive_3d
            module procedure rAllocatableIsAlive_4d
            module procedure iAllocatableIsAlive_4d
         end interface varIsAliveA

         contains


            subroutine allocator_errorHandler(status,size,data_kind,var_id,ignoreError)
               integer                        :: status
! optional arguments
               logical, optional              :: ignoreError
               integer, optional              :: size
               integer, optional              :: data_kind
               character(len=*), optional     :: var_id 
! local variables
               logical                        :: ignoreErrorLocal, isFatal, isVerbose
               integer                        :: kbytesRequested, sizeLocal, kindLocal
               character(len=charArrayLength) :: ErrorMessage, InfoMessage
               character(len=charArrayLength) :: trace_result
               character(len=charArrayLength) :: kbytesAvailable
               character(len=1)               :: dummy
               real(kind=kreal)               :: mbytesInUse, mbytesRequested, mbytesMax

               isFatal = .True.
               call allocator_get_verbose(isVerbose)

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

               if (present(var_id).and.isVerbose) then
                  if (var_id.ne.'None') then
                     write (*,*) 'Warning: at allocation for var. id '//var_id(1:42)
                  endif
               endif
 
               InfoMessage  = 'Unknown error!'
               ErrorMessage = InfoMessage
               trace_result = 'unknown (stack trace info unavailable)'

               call get_traceback_info (printAllocCaller,trace_result,charArrayLength)
               
               select case (status)
                  case (alreadyAllocated)      
                     InfoMessage  = 'Trying to allocate variable already allocated.'
                     ErrorMessage = 'Memory allocation failure at '//trace_result(1:42)
                     isFatal = .True.

                  case (alreadyDeallocated)     
                     InfoMessage  = 'Trying to deallocate variable already deallocated.'
                     ErrorMessage = 'Trying to deallocate variable already deallocated'
                     isFatal = .False.
                  
                  case (allocateFailure)
                     call get_memory_info(kbytesAvailable,charArrayLength)
                     kbytesRequested = (sizeLocal*kindLocal)/b2kbfac

                     InfoMessage  = 'Error! Allocation failed!'
                     write (InfoMessage,'(a,i9,a)')                       &
                        'Asked for ',kbytesRequested,' Kb, '//  &
                        'memory available is '//kbytesAvailable(1:24)
!                        radovan: this error message can be very confusing
!                        real example: "Asked for   1274446 Kb, memory available is 45590396 Kb"
!                        the memory available is the memory on the node but it may not be available
!                        if the batch system limits it

                     write (*,*) 'Please verify:'
                     write (*,*) ' - if the amount of memory installed in the system is sufficient'
                     write (*,*) '   taking into account the total memory requested (master+slaves)'
                     write (*,*) ' - if there are any limitations enforced by the shell, queueing '
                     write (*,*) '   system or OS'

                     ErrorMessage = 'Memory allocation failure at '//trace_result(1:42)
                     isFatal = .True.

                  case (deallocateFailure)
                     InfoMessage  = 'Memory deallocation failed at '//trace_result(1:42)
                     ErrorMessage = 'deallocation failed!'
                     isFatal = .False.

                  case (malformedAllocCall)
                     InfoMessage  = 'Malformed call to alloc() at '//trace_result(1:42)
                     ErrorMessage = 'Please submit a bug report.'
                     isFatal = .True.

                  case (reqMemoryOverLimit)
                     InfoMessage  = 'Error! Allocation failed!'

                     mbytesRequested = real(sizeLocal)*real(kindLocal)/b2mbfac
!                    write (*,*) 'debug, requested :',mbytesRequested

                     call allocator_get_meminuse(mbytesInUse)  
                     mbytesInUse = mbytesInUse/b2mbfac
!                    write (*,*) 'debug, in use :',mbytesInUse

                     call allocator_get_max(mbytesMax)
                     mbytesMax   = mbytesMax/b2mbfac
!                    write (*,*) 'debug, maximum :',mbytesMax

                     write (InfoMessage,'(a,f12.0,a,f12.0,a,f12.0,a)')  &
     &                  'Request (',mbytesRequested,' Mb) +'//          &
     &                  ' memory in use (',mbytesInUse,' Mb)'//         &
     &                  '  >  user-defined limit (',mbytesMax,' Mb)'
#ifdef PRG_DIRAC
                     write (*,*) 'Make -aw',((mbytesRequested+mbytesInUse+1.0D0)/8.0D0)
#endif

                     ErrorMessage = 'Memory allocation failure at '//trace_result(1:42)
                     isFatal = .True.

                  case (negativeSizeCall)
                     InfoMessage  = 'Negative array dimension (e.g. a '//   &
     &                              'value larger than what the variable can hold). '// &
     &                              'Allocation is not possible, program will exit now.'
                     ErrorMessage = 'Memory allocation failure at '//trace_result(1:42)
                     isFatal = .True.

               end select

               if (isFatal) then
                  if (ignoreErrorLocal) then
                     if (isVerbose) &
     &               write (*,*) 'Non-fatal error at '//trace_result(1:42)//'. The program will continue.'
                  else
                     write (*,*) InfoMessage
                     call quit(ErrorMessage)
                  endif
               else
                  if (isVerbose) then
                     write (*,*) 'Non-fatal error at '//trace_result(1:42)
                     write (*,*) InfoMessage
                     write (*,*) 'The program will continue.'
                  endif
               endif
               
            end subroutine allocator_errorHandler


            subroutine quit(errormessage)
               character(len=*) :: errormessage
 
               write (*,*) errormessage
               call exit()

            end subroutine quit

!
! other error-handling related routines
!

!
! 1d-objects
!
            logical function rPointerIsAlive_1d(data)
               real(kind=kreal), pointer :: data(:)

               rPointerIsAlive_1d = .false.
               if ( associated(data) ) rPointerIsAlive_1d = .true.
            end function rPointerIsAlive_1d

            logical function iPointerIsAlive_1d(data)
               integer(kind=kint), pointer :: data(:)

               iPointerIsAlive_1d = .false.
               if ( associated(data) ) iPointerIsAlive_1d = .true.
            end function iPointerIsAlive_1d

            logical function rAllocatableIsAlive_1d(data)
               real(kind=kreal), allocatable :: data(:)

               rAllocatableIsAlive_1d = .false.
               if ( allocated(data) ) rAllocatableIsAlive_1d = .true.
            end function rAllocatableIsAlive_1d

            logical function iAllocatableIsAlive_1d(data)
               integer(kind=kint), allocatable :: data(:)

               iAllocatableIsAlive_1d = .false.
               if ( allocated(data) ) iAllocatableIsAlive_1d = .true.
            end function iAllocatableIsAlive_1d
!
! 2d-objects
!
            logical function rPointerIsAlive_2d(data)
               real(kind=kreal), pointer :: data(:,:)

               rPointerIsAlive_2d = .false.
               if ( associated(data) ) rPointerIsAlive_2d = .true.
            end function rPointerIsAlive_2d

            logical function iPointerIsAlive_2d(data)
               integer(kind=kint), pointer :: data(:,:)

               iPointerIsAlive_2d = .false.
               if ( associated(data) ) iPointerIsAlive_2d = .true.
            end function iPointerIsAlive_2d

            logical function rAllocatableIsAlive_2d(data)
               real(kind=kreal), allocatable :: data(:,:)

               rAllocatableIsAlive_2d = .false.
               if ( allocated(data) ) rAllocatableIsAlive_2d = .true.
            end function rAllocatableIsAlive_2d

            logical function iAllocatableIsAlive_2d(data)
               integer(kind=kint), allocatable :: data(:,:)

               iAllocatableIsAlive_2d = .false.
               if ( allocated(data) ) iAllocatableIsAlive_2d = .true.
            end function iAllocatableIsAlive_2d

!
! 3d-objects
!
            logical function rPointerIsAlive_3d(data)
               real(kind=kreal), pointer :: data(:,:,:)

               rPointerIsAlive_3d = .false.
               if ( associated(data) ) rPointerIsAlive_3d = .true.
            end function rPointerIsAlive_3d

            logical function iPointerIsAlive_3d(data)
               integer(kind=kint), pointer :: data(:,:,:)

               iPointerIsAlive_3d = .false.
               if ( associated(data) ) iPointerIsAlive_3d = .true.
            end function iPointerIsAlive_3d

            logical function rAllocatableIsAlive_3d(data)
               real(kind=kreal), allocatable :: data(:,:,:)

               rAllocatableIsAlive_3d = .false.
               if ( allocated(data) ) rAllocatableIsAlive_3d = .true.
            end function rAllocatableIsAlive_3d

            logical function iAllocatableIsAlive_3d(data)
               integer(kind=kint), allocatable :: data(:,:,:)

               iAllocatableIsAlive_3d = .false.
               if ( allocated(data) ) iAllocatableIsAlive_3d = .true.
            end function iAllocatableIsAlive_3d

!
! 4d-objects
!
            logical function rPointerIsAlive_4d(data)
               real(kind=kreal), pointer :: data(:,:,:,:)

               rPointerIsAlive_4d = .false.
               if ( associated(data) ) rPointerIsAlive_4d = .true.
            end function rPointerIsAlive_4d

            logical function iPointerIsAlive_4d(data)
               integer(kind=kint), pointer :: data(:,:,:,:)

               iPointerIsAlive_4d = .false.
               if ( associated(data) ) iPointerIsAlive_4d = .true.
            end function iPointerIsAlive_4d

            logical function rAllocatableIsAlive_4d(data)
               real(kind=kreal), allocatable :: data(:,:,:,:)

               rAllocatableIsAlive_4d = .false.
               if ( allocated(data) ) rAllocatableIsAlive_4d = .true.
            end function rAllocatableIsAlive_4d

            logical function iAllocatableIsAlive_4d(data)
               integer(kind=kint), allocatable :: data(:,:,:,:)

               iAllocatableIsAlive_4d = .false.
               if ( allocated(data) ) iAllocatableIsAlive_4d = .true.
            end function iAllocatableIsAlive_4d

      end module allocator_eh
