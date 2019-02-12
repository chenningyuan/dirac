!This is a stand-alone version of exacorr that reads all information
!from an interface file.

   program exacorr_main

#if defined (VAR_MPI)
      use interface_to_mpi
#endif
      use exacorr_ao_to_mo,only     : exacorr_compute_mo_integrals
      use talsh_ao_to_mo,only       : talsh_compute_mo_integrals
      use talsh_cc,only             : talsh_cc_driver
      use exacorr_cc,only           : exacorr_cc_driver
      use exacc_cfg
      use exacorr_global

      implicit none

      logical mo_transform, cc_calculation
      integer              :: numproc, mytid
      integer              :: imo, nmo(4)
      logical              :: input_found=.false.
      integer              :: file_unit
      integer, allocatable :: mo_list(:)
      integer, allocatable :: mokr_occ(:), mokr_vir(:)
      integer              :: nesh, kvec
      integer              :: i

!     Initialize MPI variables
      mytid   = 0
      numproc = 0

#if defined (VAR_MPI)
      call interface_mpi_init()
      call interface_mpi_comm_rank(global_communicator,mytid)
      call interface_mpi_comm_size(global_communicator,numproc)
#endif

      mo_transform   = .false.
      cc_calculation = .true.

      if (mytid == 0) then
        call get_free_fileunit(file_unit)
        call read_menu_input('exacorr.inp',file_unit,'**EXACC',input_found) 
        if (.not.input_found) call quit('Input was not found!')
      end if

#if defined (VAR_MPI)
      call interface_mpi_bcast_l0(exa_input%talsh,1,0,global_communicator)
      if (mytid /= 0 .and. exa_input%talsh) call interface_mpi_finalize()
#else
!     serial version of DIRAC
      if (.not.exa_input%talsh) then
!        EXATENSOR version does not work with serial version of DIRAC
         print*, "Attempting to use exatensor with a non MPI version"
         call quit ("error in setup: exatensor requires MPI")
      end if
#endif

!     Initialize global data (read from DFCOEF or similar by master)
      call set_talsh_only(exa_input%talsh)
      call initialize_global_data
      if (mytid == 0) print*," Initialized global data"

      if (mo_transform) then
!        We have no input driver for this test code, so we hardwire this here.
         nmo = 8
         allocate(mo_list(nmo(1)+nmo(2)+nmo(3)+nmo(4)))
         do imo = 1, nmo(1)
             mo_list(imo)                      = imo
             mo_list(imo+nmo(1))               = imo
             mo_list(imo+nmo(1)+nmo(2))        = imo
             mo_list(imo+nmo(1)+nmo(2)+nmo(3)) = imo
         end do

!        call talsh_compute_mo_integrals (nmo, mo_list)
         call exacorr_compute_mo_integrals (nmo, mo_list)
      end if

      if (mytid == 0) then
        nesh = get_nsolutions() / 2

        exa_input%nkr_occ = 0
        call numlst (string_occupied,kvec,1,1,nesh,1,exa_input%nkr_occ)
        allocate (mokr_occ(exa_input%nkr_occ))
        call numlst (string_occupied,mokr_occ,exa_input%nkr_occ,1,nesh,1,exa_input%nkr_occ)
        exa_input%nkr_occ = 2 * exa_input%nkr_occ
        allocate (exa_input%mokr_occ(exa_input%nkr_occ))
        do i = 1, exa_input%nkr_occ
          exa_input%mokr_occ(i) = (mokr_occ(1)-1)*2+i
        end do
        deallocate (mokr_occ)

        exa_input%nkr_vir = 0
        call numlst (string_virtual,kvec,1,1,nesh,1,exa_input%nkr_vir)
        allocate (mokr_vir(exa_input%nkr_vir))
        call numlst (string_virtual,mokr_vir,exa_input%nkr_vir,1,nesh,1,exa_input%nkr_vir)
        exa_input%nkr_vir = 2 * exa_input%nkr_vir
        allocate (exa_input%mokr_vir(exa_input%nkr_vir))
        do i = 1, exa_input%nkr_vir
          exa_input%mokr_vir(i) = exa_input%mokr_occ(exa_input%nkr_occ)+i
        end do
        deallocate (mokr_vir)
      end if

      if (cc_calculation) then
         if (mytid == 0) print*, "a CC calculation with",exa_input%nkr_occ," occupied"," and ",exa_input%nkr_vir," virtuals"
         if (exa_input%talsh) then
            if (mytid == 0) call talsh_cc_driver (exa_input)
         else
#if defined (VAR_MPI)
!           slaves need the same info to get started
            call interface_mpi_bcast(exa_input%nkr_occ,1,0,global_communicator)
            if (mytid /= 0) allocate(exa_input%mokr_occ(exa_input%nkr_occ))
            call interface_mpi_bcast(exa_input%nkr_vir,1,0,global_communicator)
            if (mytid /= 0) allocate(exa_input%mokr_vir(exa_input%nkr_vir))
            call exacc_sync_cw(0)
#endif
            call exacorr_cc_driver (exa_input)
         end if
      end if

      deallocate(exa_input%mokr_occ)
      deallocate(exa_input%mokr_vir)

#if defined (VAR_MPI) 
      call interface_mpi_finalize()
#endif

   end program exacorr_main
