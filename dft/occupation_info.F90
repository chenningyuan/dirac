
!> \brief   write/read occupation information to/from file
!> \author  radovan bast
!> \date    2010-05-03
!> \warning no mpi bcasting done on these parameters

module occupation_info

   implicit none

   public occupation_info_write
   public occupation_info_read
   public occupation_info_delete

   private

!  -----------------------------------------------------------------------------
!  public

   integer, public :: nr_open_shells
   integer, public :: nr_electrons_inactive
   integer, public :: nr_electrons_active
   integer, public :: nr_electrons_total

!  -----------------------------------------------------------------------------
!  private

   character(*), parameter :: file_name = 'occupation_info'
   integer,      parameter :: file_unit = 137
   logical                 :: file_exists

contains

!> \brief reset parameters and deallocate arrays
   subroutine occupation_info_delete()

      nr_open_shells        = 0
      nr_electrons_inactive = 0
      nr_electrons_active   = 0
      nr_electrons_total    = 0

   end subroutine

!> \brief write occupation info to file
   subroutine occupation_info_write()

#include "dcborb.h"

      call occupation_info_delete()

      inquire(file  = file_name, exist = file_exists)

      if (file_exists) then
         close(file_unit, status = 'keep')
         return
      else
         open(file_unit,            &
              file   = file_name,   &
              status = 'new',       &
              form   = 'formatted', &
              access = 'sequential')
         rewind(file_unit)
      end if

      nr_open_shells        = nasht
      nr_electrons_inactive = nelect
      nr_electrons_active   = naelec
      nr_electrons_total    = nr_electrons_inactive + nr_electrons_active

      write(file_unit, *) nr_open_shells
      write(file_unit, *) nr_electrons_inactive
      write(file_unit, *) nr_electrons_active
      write(file_unit, *) nr_electrons_total

      close(file_unit, status = 'keep')

      call occupation_info_delete()

   end subroutine

!> \brief read occupation info from file
   subroutine occupation_info_read()

      call occupation_info_delete()

      inquire(file  = file_name, exist = file_exists)

      if (file_exists) then
         open(file_unit,            &
              file   = file_name,   &
              status = 'unknown',   &
              form   = 'formatted', &
              access = 'sequential')
         rewind(file_unit)
      else
         call quit('cannot read from file occupation_info')
      end if

      read(file_unit, *) nr_open_shells
      read(file_unit, *) nr_electrons_inactive
      read(file_unit, *) nr_electrons_active
      read(file_unit, *) nr_electrons_total

      close(file_unit, status = 'keep')

   end subroutine

end module
