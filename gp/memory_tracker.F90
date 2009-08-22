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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! memory tracker module
!
! this provides functionality for the memory_allocator module to keep
! track of memory comsumption
!
!  MI/5.4.'08:
!  to compile it with GNU, gfortran-4.2 is needed (not the default gfortran-4.1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module memory_tracker

         implicit none

         private

         integer, parameter :: kchar = 1, kint = 4, kreal = 8, kcomplex = 8

         type, public :: fString
            integer            :: size = 0 
            character, pointer :: s(:) 
         end type fString

!
! the next is a kind of hack to keep track of memory addresses of the
! possibla variable types we have, since there is no equivalent to
! a void pointer in fortran... maybe redefine as an interface to pointers is better?
!
         type, public :: memoryContent
            integer           , pointer     :: iContent 
            real            , pointer       :: rContent 
        end type memoryContent
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! type definitions for the memory tracking system
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         type, public :: memChunk
            integer                  :: n_elements = 0                ! size of the memory, in units 
            integer                  :: type       = 0                ! type of variable  
            integer                  :: address = 0             ! address of the memory
            type (memoryContent)     :: content
            type (memChunk), pointer :: next, prev              ! for the linked list
            type (fString), pointer  :: routine                 ! optionally holds the name of the routine where this chunk was allocated 
            type (fString), pointer  :: variable                ! optionally holds the name of the variable this chunk was allocated for 
         end type memChunk


         type, public :: memGroup
            integer                  :: nchunks
            type (fString)           :: label
            type (memChunk), pointer :: firstActive, lastActive          ! list of allocated chunks
            type (memChunk), pointer :: firstFreed,  lastFreed           ! list of deallocated chunks
            type (memChunk), pointer :: largestActive, smallestActive    ! largest and smallest active chunks for this group
            integer                  :: largestEver, smallestEver        ! size of largest and smallest ever for this group, in bytes
         end type memGroup


         type, public :: memTracker
            integer                  :: ngroups
            type (memGroup), pointer :: first, last                      ! first and last active groups
            integer                  :: largestEver, smallestEver        ! size of globally largest and smallest allocations in bytes 
            type (fString)           :: largestLabel, smallestLabel
         end type memTracker
         
         type(memTracker), save :: memoryTracking

         public trackInitMemoryStack
         contains 

! manager

            subroutine trackInitMemoryStack
               memoryTracking%ngroups = 0
               
               memoryTracking%first   => NULL()
               memoryTracking%last    => NULL()
               
               memoryTracking%largestEver  = 0
               memoryTracking%smallestEver = 0
               memoryTracking%first   => NULL()
               
            end subroutine trackInitMemoryStack


! memGroup manipulation routines

            subroutine trackCreateGroup(label,group_ptr)
               character(len=*) :: label
               type (memGroup), pointer :: group_ptr

               allocate(group_ptr)

            end subroutine trackCreateGroup

            subroutine trackDestroyGroup(group_ptr,label)
               character(len=*) :: label
               type (memGroup), pointer :: group_ptr

               group_ptr%largestActive  => NULL()
               group_ptr%smallestActive => NULL()

               ! loop over all active chunks, and destroy all of them
               deallocate(group_ptr%firstActive)

               !gth
               deallocate(group_ptr%largestActive)
               deallocate(group_ptr%largestActive)

            end subroutine trackDestroyGroup




! registrar

            subroutine trackRegisterAllocation(size,type,variableName,trackLevel)
               integer  :: size, type
               character(len=*) :: variableName, trackLevel
               type (memChunk), pointer :: newChunk

               allocate(newChunk) 
            end subroutine trackRegisterAllocation
         
      end module memory_tracker


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! simplified memory tracker module
!
! before the full-fleged beast above is finished, this module will provide
! the most rudimentary functionality
!
      module simple_memory_tracker

         use iso_c_binding

         implicit none

         private

         integer, parameter :: kchar = 1, kint = 4, kreal = 8, kcomplex = 8

         type, public :: memChunk
            integer                  :: n_elements = 0                ! size of the memory, in units
            integer                  :: type       = 0                ! type of variable
            type (memChunk), pointer :: next, prev              ! for the linked list
         end type memChunk

         type, public :: memTracker
            integer                  :: n_created, n_destroyed
            type (memChunk), pointer :: memCreatedFirst, memCreatedLast
            type (memChunk), pointer :: memDestroyedFirst, memDestroyedLast
            integer                  :: largestSize,  largestKind
            integer                  :: smallestSize, smallestKind 
         end type memTracker

         type(memTracker), save :: memTracking

         public trackCreateTracker, trackRegisterAllocation, trackPrintStats

         contains

            subroutine trackCreateTracker

                memTracking%n_created     = 0
                memTracking%n_destroyed   = 0

                memTracking%memCreatedFirst    => NULL()
                memTracking%memCreatedLast     => NULL()
                memTracking%memDestroyedFirst  => NULL()
                memTracking%memDestroyedLast   => NULL()

                memTracking%largestSize   = 0
                memTracking%largestKind   = 0

            end subroutine trackCreateTracker

            subroutine trackPrintStats(mytid)
               integer :: mytid
               real (kind=8) :: totalMemory

               call trackTotalMemoryAllocated(totalMemory)

               if (mytid.eq.0) then
                  write (*,*) ' '
                  write (*,*) ' >>> Dynamical Memory Usage Summary' 
                  write (*,*) ' '
               endif
               write (*,*) ' >>> node',mytid,' : ',memTracking%n_created,' allocations, ',totalMemory/(1024*1024),' Mb'

!              if ( memTracking%n_created .ne. memTracking%n_destroyed ) then
!                 write (*,*) '          warning! memory was left allocated!'
!              endif

            end subroutine trackPrintStats

            subroutine trackTotalMemoryAllocated(totalMemory)
               integer :: pass
               real (kind=8) :: totalMemory
               type (memChunk), pointer :: tmp

               totalMemory = 0
               pass = 0

!              write (*,*) memTracking%memCreatedFirst%n_elements
!              write (*,*) memTracking%memCreatedLast%n_elements
!              write (*,*) memTracking%memCreatedLast%prev%n_elements
!              write (*,*) memTracking%memCreatedLast%prev%prev%n_elements

!              write (*,*) memTracking%memCreatedFirst%next%n_elements

               tmp => memTracking%memCreatedLast

               do
                  if (.not.associated(tmp)) then
                     exit
                  endif

                  pass = pass + 1
!                 write (*,*) 'debug: pass # ',pass
!                 write (*,*) tmp%n_elements

                  totalMemory = totalMemory + (tmp%n_elements * tmp%type) 
                  tmp => tmp%prev
               enddo

            end subroutine trackTotalMemoryAllocated

            subroutine trackRegisterAllocation(size,type,variableID,trackLevel)
               integer  :: size, type, trackLevel
               character(len=*) :: variableID
               type (memChunk), pointer :: newChunk

               allocate(newChunk)

               newChunk%n_elements = size
               newChunk%type       = type

               newChunk%next => NULL()
               newChunk%prev => memTracking%memCreatedLast

               memTracking%memCreatedLast => newChunk

               if ( .not.associated(memTracking%memCreatedFirst) ) then
                  memTracking%memCreatedFirst => newChunk
               else
                  if ( .not.associated(memTracking%memCreatedFirst%next)) then
                     memTracking%memCreatedFirst%next => newChunk 
                  endif 
               endif

               memTracking%n_created = memTracking%n_created + 1

!              write (*,*) '----------'
!              if ( associated(memTracking%memCreatedFirst) ) then
!              write (*,*) 'debug: first', memTracking%memCreatedFirst%n_elements
!              else 
!              write (*,*) 'debug: FIRST is null'
!              endif

!              if ( associated(memTracking%memCreatedLast) ) then
!              write (*,*) 'debug: last',  memTracking%memCreatedLast%n_elements
!              else 
!              write (*,*) 'debug: LAST is null'
!              endif

!              if ( associated(newChunk) ) then
!              write (*,*) 'debug: current', newChunk%n_elements
!              else 
!              write (*,*) 'debug: current is null'
!              endif

!              if ( associated(newChunk%prev) ) then
!              write (*,*) 'debug: current->prev', newChunk%prev%n_elements
!              else 
!              write (*,*) 'debug: current->prev is null'
!              endif

!              if ( associated(newChunk%next) ) then
!              write (*,*) 'debug: current->next', newChunk%next%n_elements
!              else 
!              write (*,*) 'debug: current->next is null'
!              endif
!              write (*,*) '----------'

            end subroutine trackRegisterAllocation 

            subroutine trackRegisterDeallocation(variableID,trackLevel)
               integer          :: trackLevel
               character(len=*) :: variableID
               type (memChunk), pointer :: tmp
 
               tmp = memTracking%memCreatedLast
!              for ( ; tmp%prev .eq. NULL; tmp=tmp%prev) 
!                 if (variableID .eq. tmp%id ) then
!                    if ( .not.associated(memTracking%memDestroyedFirst) ) then
!                    endif
!                 endif
!              enddo

            end subroutine trackRegisterDeallocation               


     end module simple_memory_tracker
