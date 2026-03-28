!
!    Interface with ADF. We have three situations to deal with:
!    0. no access to ADF routines
!       - compile a few placeholder routines (needed to avoid undesirable code duplication in ROSE)
!       - Slater type orbitals can not be used
!    1. only access to KF modules
!       - read only geometry and MO coefficients using direct KF read commands
!       - read AO overlap from precomputed TAPE15
!    2. full access to AMS modules
!       - use AMS modules to initialize and read fragments
!       - compute AO overlap with the appropriate ADF subroutines
!
!    The placeholders below are replaced by actual routines if situation 2 applies

    subroutine read_in_ADF (ao_basis,spinors)

       use read_xyz, only : molecule_info
       use rose_ao
       use rose_mo

       implicit none
       type(basis_set_info_t), target, intent(inout) :: ao_basis
       type(mo_basis), intent(inout)       :: spinors
 
    end subroutine read_in_ADF

    subroutine get_ao_overlap_ADF (ao_basis_bra,ao_basis_ket,sao)

       use rose_ao

       implicit none
       type(basis_set_info_t), intent(in), target :: ao_basis_bra, ao_basis_ket
       real(8)                                    :: sao(:,:)

    end subroutine get_ao_overlap_ADF

    subroutine write_ROSE_results (spinors,ribo_occ_per_frag, ribo_vir_per_frag, RIBO_Fock,ncore)

      use rose_mo

      implicit none
      type(mo_basis), intent(inout) :: spinors
      integer, intent(in)           :: ribo_occ_per_frag(*)
      integer, intent(in)           :: ribo_vir_per_frag(*)
      real(8), intent(in)           :: RIBO_Fock(spinors%nmo-ncore,spinors%nmo-ncore,1)
      integer, intent(in)           :: ncore

    end subroutine write_ROSE_results
