!
!    Interface with ADF. Two situations remain:
!    0. no access to ADF routines, which is the standalone build
!       - compile the placeholder routines below, which abort if ever reached
!       - a Slater basis is read from an h5 file produced by rose_kf_bridge.py,
!         and its overlaps are computed by slater_overlap.F90
!    1. full access to AMS modules, which is the build inside AMS
!       - use AMS modules to initialize and read fragments
!       - compute AO overlap with the appropriate ADF subroutines
!
!    The placeholders below are replaced by the real routines of rose_inAMS.f90
!    when situation 1 applies.
!
!    There used to be a third situation, linking only to CoreKF so that ROSE could
!    read an rkf and a precomputed TAPE15 directly. That route, ADF_interface.F90
!    and the SCM library discovery in CMakeLists.txt were removed once the KF
!    bridge and slater_overlap made them unnecessary; see the git history.
!
!    The four adf_* procedures at the end serve the Slater (basis_type 2) pathway
!    of import_fragments.F90 and have no meaning without AMS: there is no way to
!    answer them from a Gaussian basis. They therefore abort rather than return a
!    value. They are bare external procedures on both sides, declared at the call
!    site with "integer, external" or called directly, so nothing checks that these
!    declarations still agree with the ones in rose_inAMS.f90 -- keep them in step
!    by hand whenever that file changes.

    subroutine read_structure_in_ADF (molecule)

       use read_xyzfile, only : molecule_info

       implicit none
       type(molecule_info), intent(inout)  :: molecule

    end subroutine read_structure_in_ADF

    subroutine read_basis_in_ADF (molecule_name,ao_basis,spinors)

       use rose_ao
       use rose_mo

       implicit none
       character(*), intent(in)              :: molecule_name
       type(basis_set_info_t), intent(inout) :: ao_basis
       type(mo_basis), intent(inout)         :: spinors

    end subroutine read_basis_in_ADF

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

    integer function adf_n_atomic_centers (ao_basis)

!      Number of AO centers of an ADF (STO) basis. See rose_inAMS.f90.

       use rose_ao

       implicit none
       type(basis_set_info_t), intent(in) :: ao_basis

       adf_n_atomic_centers = 0
       error stop 'adf_n_atomic_centers: Slater basis sets are only available in the AMS build'

    end function adf_n_atomic_centers

    subroutine adf_atomic_centers (ao_basis,n_centers,center_coord)

!      Coordinates of those centers. See rose_inAMS.f90.

       use rose_ao

       implicit none
       type(basis_set_info_t), intent(in) :: ao_basis
       integer, intent(in)                :: n_centers
       real(8), intent(out)               :: center_coord(3,n_centers)

       center_coord = 0.D0
       error stop 'adf_atomic_centers: Slater basis sets are only available in the AMS build'

    end subroutine adf_atomic_centers

    subroutine adf_place_fragment (ao_basis,x)

!      Move a single-atom ADF fragment basis to the position of its atom.
!      See rose_inAMS.f90.

       use rose_ao

       implicit none
       type(basis_set_info_t), intent(inout) :: ao_basis
       real(8), intent(in)                   :: x(3)

       error stop 'adf_place_fragment: Slater basis sets are only available in the AMS build'

    end subroutine adf_place_fragment

    integer function adf_frozen_core_capacity (atom_number)

!      Number of electrons that ADF has frozen away for this element.
!      See rose_inAMS.f90.

       implicit none
       integer, intent(in) :: atom_number

       adf_frozen_core_capacity = 0
       error stop 'adf_frozen_core_capacity: the AMS frozen core is only available in the AMS build'

    end function adf_frozen_core_capacity
