module visual_utils

  implicit none

  public newunit
  public visual_quit
  public visual_walltime

contains

  function visual_walltime()
!  ----------------------------------------------------------------------------
    use interface_to_mpi
    real(8) :: visual_walltime
!  ----------------------------------------------------------------------------
#ifdef VAR_MPI
    visual_walltime = interface_MPI_WTIME()
#endif
    call cpu_time(visual_walltime)
   end function visual_walltime


   integer function newunit(unit) result(n)
!  ----------------------------------------------------------------------------
      integer, intent(out), optional:: unit
      integer, parameter:: nmin = 10
      integer, parameter:: nmax = 1000
      character(len = 100):: errmsg
      logical:: is_used

!     returns lowest i/o unit number not in use
!     adapted from http://fortranwiki.org/fortran/show/newunit
!     gosia, May 2022
!  ----------------------------------------------------------------------------
      do n = nmin, nmax
          inquire(unit = n, opened = is_used)
          if (.not. is_used) then
              if (present(unit)) unit = n
              return
          end if
      end do
      write (errmsg, *) 'newunit ERROR: no available file unit in a range ',&
                       nmin, nmax, ' Stopping now!'
      call visual_quit(errmsg)
   end function


   subroutine visual_quit(text)
!  ----------------------------------------------------------------------------
      character(*):: text
      character(len = len_trim(text)):: local_text
!  ----------------------------------------------------------------------------

      write(local_text, '(a)') trim(text)
      write(*,'(a)') local_text
      stop

   end subroutine

end module
