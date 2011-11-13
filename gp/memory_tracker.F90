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
! $Id: memory_tracker.F90 7009 2008-06-07 23:10:40Z gomes $
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
            integer(kind=kint)  :: size = 0 
            character, pointer :: s(:) 
         end type fString

!
! the next is a kind of hack to keep track of memory addresses of the
! possibla variable types we have, since there is no equivalent to
! a void pointer in fortran... maybe redefine as an interface to pointers is better?
!
         type, public :: memoryContent
            integer(kind=kint), pointer     :: iContent 
            real(kind=kreal), pointer       :: rContent 
        end type memoryContent
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! type definitions for the memory tracking system
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         type, public :: memChunk
            integer(kind=kint)       :: size = 0                ! size of the memory, in bytes
            integer(kind=kint)       :: address = 0             ! address of the memory
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
            integer(kind=kint)       :: largestEver, smallestEver        ! size of largest and smallest ever for this group, in bytes
         end type memGroup


         type, public :: memTracker
            integer                  :: ngroups
            type (memGroup), pointer :: first, last                      ! first and last active groups
            integer(kind=kint)       :: largestEver, smallestEver        ! size of globally largest and smallest allocations in bytes 
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

         
end module memory_tracker
