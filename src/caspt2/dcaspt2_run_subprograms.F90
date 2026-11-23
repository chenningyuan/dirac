subroutine dcaspt2_run_subprograms(ivo_mode)
    use dcaspt2_restart_file, only: read_and_validate_restart_file
    use module_file_manager, only: open_formatted_file
    use module_global_variables
#ifdef HAVE_MPI
    use module_mpi, only: bcast_wrapper
#endif
    use module_realonly, only: check_realonly
    use module_validation, only: validate_ndet
    use read_input_module, only: read_input
    implicit none
#include "dcbpsi.h"
    logical, intent(in) :: ivo_mode
    integer :: unit_input, i
    integer :: ciroots_size_1, ciroots_size_2, max_selectroot_list_size
    character(:), allocatable :: filename

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    if (rank == 0) then
        print *, ' '
        print *, 'Reading MRCONEE (1-e integrals)'
    end if
    allocate (filename, source='MRCONEE')
    call check_dirac_integer_size(filename)
    call read_mrconee(filename)
    call check_realonly

    if (docountndet) then
        call search_cas_configuration ! If docountndet is true, only search the number of CASCI configurations and print them
        return ! skip subprograms if docountndet is true
    end if

    if (ivo_mode) then
        if (.not. doivo) return
        call r4divo_co
        call dcaspt2_deallocate
        ! skip the following subprogram if ivo_mode is true
        ! if docasci or docaspt2 is true, this subroutine will be called again with ivo_mode = .false.
        return
    end if

#ifdef HAVE_MPI
    if (rank == 0) then
        ciroots_size_1 = size(caspt2_ciroots, 1)
        ciroots_size_2 = size(caspt2_ciroots, 2)
        max_selectroot_list_size = size(max_selectroot_list)
    end if
    call bcast_wrapper(docasci, 0)
    call bcast_wrapper(docaspt2, 0)
    call bcast_wrapper(ciroots_size_1, 0)
    call bcast_wrapper(ciroots_size_2, 0)
    call bcast_wrapper(max_selectroot_list_size, 0)
    if (.not. allocated(caspt2_ciroots)) then
        allocate(caspt2_ciroots(ciroots_size_1, ciroots_size_2))
    end if
    if (.not. allocated(max_selectroot_list)) then
        allocate(max_selectroot_list(max_selectroot_list_size))
    end if
    call bcast_wrapper(caspt2_ciroots, 0)
    call bcast_wrapper(max_selectroot_list, 0)
#endif
    do i = 1, size(caspt2_ciroots, 1)
        totsym = caspt2_ciroots(i, 1)
        selectroot = caspt2_ciroots(i, 2)
        nroot = max_selectroot_list(totsym)
        e2_subspace = 0
        sumc2 = 0
        sumc2_subspace = 0
        if (enable_restart) call read_and_validate_restart_file

        if (rank == 0) then
            print '(a)', 'start casci/caspt2, rasci/raspt2 calculation'
            print '(a,x,i0)', "totsym:", totsym
            print '(a,x,i0)', "selectroot:", selectroot
            print '(a,x,i0)', "nroot:", nroot
        end if
        if (casci_done(totsym)) call validate_ndet

        if (docasci .and. .not. casci_done(totsym)) then
            call r4dcasci
            casci_done(totsym) = .true.
        end if

        if (docaspt2) then
            call r4dcaspt2_tra
        end if
    end do
    call dcaspt2_deallocate

end subroutine dcaspt2_run_subprograms
