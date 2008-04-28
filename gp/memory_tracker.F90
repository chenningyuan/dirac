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
