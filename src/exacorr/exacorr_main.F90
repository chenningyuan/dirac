
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
      use exacorr_datatypes,only    : exacc_input

      implicit none

      logical mo_transform, cc_calculation
      integer              :: imo, nmo(4)
      integer, allocatable :: mo_list(:)
      type(exacc_input)    :: exa_input

      mo_transform   = .false.
      cc_calculation = .true.
      exa_input%talsh = .false.
      exa_input%ccd = .false.
      exa_input%lambda = .false.
      exa_input%noints = .false.
      exa_input%exa_blocksize = 75
      exa_input%moint_scheme = 1

#if defined (VAR_MPI)
      call interface_mpi_init()
#endif

      nmo = 8
      allocate(mo_list(nmo(1)+nmo(2)+nmo(3)+nmo(4)))
      do imo = 1, nmo(1)
          mo_list(imo)                      = imo
          mo_list(imo+nmo(1))               = imo
          mo_list(imo+nmo(1)+nmo(2))        = imo
          mo_list(imo+nmo(1)+nmo(2)+nmo(3)) = imo
      end do

      if (mo_transform) then
         call talsh_compute_mo_integrals (nmo, mo_list)
         call exacorr_compute_mo_integrals (nmo, mo_list)
      end if

      !Hard-coded water input to correspond with x2c test
      exa_input%nkr_occ = 3
      exa_input%nkr_vir = 12
      allocate(exa_input%mokr_occ(exa_input%nkr_occ))
      allocate(exa_input%mokr_vir(exa_input%nkr_vir))
      exa_input%mokr_occ = (/3, 4, 5/)
      exa_input%mokr_vir = (/6, 7, 8, 9, 10 , 11, 12, 13, 14, 15, 16, 17/)

      if (cc_calculation) then
         print*, "a CC calculation with",exa_input%nkr_occ," occupied"," and ",exa_input%nkr_vir," virtuals"
         if (exa_input%talsh) then
            call talsh_cc_driver (exa_input)
         else
            call exacorr_cc_driver (exa_input)
         end if
      end if

      deallocate(exa_input%mokr_occ)
      deallocate(exa_input%mokr_vir)

#if defined (VAR_MPI) 
      call interface_mpi_finalize()
#endif

   end program exacorr_main
