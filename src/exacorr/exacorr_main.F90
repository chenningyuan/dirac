
!This is a stand-alone version of exacorr that reads all information
!from an interface file.

   program exacorr_main

#if defined (VAR_MPI) 
      use interface_to_mpi
#endif

      use exacorr_global,only   : get_talsh_only, set_talsh_only
      use exacorr_ao_to_mo,only : exacorr_compute_mo_integrals
      use talsh_ao_to_mo,only   : talsh_compute_mo_integrals
      use talsh_cc,only         : cc_driver
      use exacorr_cc,only       : exacorr_cc_driver

      implicit none

      logical mo_transform, cc_calculation 
      integer              :: imo, nmo(4)
      integer, allocatable :: mo_list(:)
      integer              :: nocc, nvir
      integer, allocatable :: mo_occ(:)
      integer, allocatable :: mo_vir(:)

      mo_transform   = .true.
      cc_calculation = .true.
      call set_talsh_only(.false.)
      

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

      nocc = 2
      nvir  = 6
      allocate(mo_occ(nocc))
      do imo = 1, nocc
          mo_occ(imo) = imo
      end do
      allocate(mo_vir(nvir))
      do imo = 1, nvir
          mo_vir(imo) = nocc+imo
      end do

      if (cc_calculation) then
         print*, "a CC calculation with",nocc," occupied"," and ",nvir," virtuals"
         if (get_talsh_only()) then
            call cc_driver (nocc,nvir,mo_occ,mo_vir)
         else
            call exacorr_cc_driver (nocc,nvir,mo_occ,mo_vir)
            print*, "exacorr cc not yet available"
         end if
      end if

#if defined (VAR_MPI) 
      call interface_mpi_finalize()
#endif

   end program exacorr_main
