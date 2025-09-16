module exacorr_utils

!  Collection of simple utility functions used in exacorr

  implicit none
  private

  public relative_error
  public get_free_fileunit
  public quicksort
  public shift_orbital_energy
 
  contains

    real(8) function relative_error(quantity_1,quantity_2)
    ! Gives percentage error for values > 1, otherwise absolute error

      real(8),intent(in):: quantity_1,quantity_2
      real(8)           :: error

      error = quantity_1 - quantity_2

      if (abs(quantity_1).gt.1.0) then
        relative_error = abs (error / quantity_1)
      else
        relative_error = abs (error)
      end if

    end function relative_error
       
    subroutine get_free_fileunit(unit)
    ! idea from http://fortranwiki.org, will become obsolete once fortran2008 is supported everywhere
    ! gets free file unit
      integer, intent(out) :: unit
      integer, parameter   :: LUN_MIN=100, LUN_MAX=1000
      logical :: opened
      integer :: lun
      unit=-1
      do lun=LUN_MIN,LUN_MAX
        inquire(unit=lun,opened=opened)
        if (.not. opened) then
          unit=lun
          exit
        end if
      end do
    end subroutine get_free_fileunit

    recursive subroutine quicksort(mo_list, i_beg, i_end)

        implicit none

        integer,intent(inout) :: mo_list(:)
        integer, intent(in) :: i_beg, i_end
        integer :: i,j
        integer :: ref, tem

        i=floor((float(i_end)+float(i_beg))/2)
        ref=mo_list(i)
        i=i_beg
        j=i_end

        do
          do while (mo_list(i)<ref)
            i=i+1
          end do
          do while (mo_list(j)>ref)
            j=j-1
          end do
          if (i<j) then
            tem=mo_list(i)
            mo_list(i)=mo_list(j)
            mo_list(j)=tem
            i=i+1
            j=j-1
          else if(i.eq.j) then
            i=i+1
            EXIT
          else
            EXIT
          end if
        end do

        if (i_beg<j) call quicksort(mo_list, i_beg, j)
        if (i<i_end) call quicksort(mo_list, i, i_end)

    end subroutine quicksort

    subroutine shift_orbital_energy(eps_vir,eps_occ,level_shift)

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
          use interface_to_mpi

          implicit none 

          integer :: my_MPI_rank
          integer :: my_MPI_master=0
#endif

          real(8), intent(inout)   :: eps_vir(:),eps_occ(:)
          real(8), intent(in)      :: level_shift
          
          integer            :: print_level=1
          integer            :: lumo,i
          real(8)            :: e_homo, e_lumo
          real(8), parameter :: THRESHOLD=1.D-8
          
          e_homo = maxval(eps_occ)
          lumo   = minloc(eps_vir,1)

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
          call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
#endif

          if (level_shift.gt.1.0D-14) then
            e_lumo = minval(eps_vir)

            do i=1,size(eps_vir)
              eps_vir(i)=eps_vir(i)+level_shift

              if (print_level.gt.-1) then
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
                if (my_MPI_rank == my_MPI_master) write(*,*) "LUMO shifted from  ", e_lumo," to ",eps_vir(lumo)
#else           
                write(*,*) "LUMO shifted from  ", e_lumo," to ",eps_vir(lumo)
#endif
              end if

            end do
          end if
          
          if (eps_vir(lumo)-e_homo<THRESHOLD) then
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
            if (my_MPI_rank == my_MPI_master) write(*,*) 'WARNING: Negative HOMO-LUMO gap (.LSHIFT possible)'
#else           
            write(*,*) 'WARNING: Negative HOMO-LUMO gap (.LSHIFT possible)'
#endif
          end if

    end subroutine shift_orbital_energy

end module exacorr_utils

