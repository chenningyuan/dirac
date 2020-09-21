!This module defines the data types used in the ExaCorr module of DIRAC
       module exacorr_datatypes

        use exatensor
        use talsh
        use, intrinsic:: ISO_C_BINDING

        implicit none
        private
!PARAMETERS:
!Generic:
        integer, parameter, private:: REALD=8 ! default real for this module
!DERIVED TYPES:
!       Basis function info:
        type, public:: basis_func_info_t
         integer:: orb_momentum=-1       ! orbital momentum: (0,1,2,3,...)
         integer:: atom_number=-1        ! atom sequential number [1..MAX] (if centered on a nucleus of some atom)
         integer:: atom_element=0        ! atomic element number (if centered on a nucleus of some atom)
         integer:: n_primitives=0        ! number of primitives in the contraction
         real(REALD):: coord(1:3)        ! coordinate of the basis function center (typically a nucleus)
         real(REALD), pointer :: exponent(:)
         real(REALD), pointer :: coefficient(:)
        end type basis_func_info_t
!       Basis set info:
        type, public:: basis_set_info_t
         integer:: nshells=0             ! number of shells in the basis set
         integer:: nao=0                 ! number of basis functions
         integer:: basis_angular         ! 1=cartesian, 2=spherical
         type(basis_func_info_t), pointer:: gtos(:)
         integer, pointer :: shell_indices(:)
        end type basis_set_info_t

        type, public :: spinor_info_t
           integer                 :: nocc, nvir             ! number of spinors
           integer, allocatable    :: occ(:), vir(:)         ! list of active spinors
           real(8), allocatable    :: eps_occ(:), eps_vir(:) ! orbital energies
        end type spinor_info_t

!       Input parameter
        type, public :: exacc_input
            logical                 :: talsh = .true.           ! if false, use ExaTensor
            logical                 :: ccd = .false.            ! if false=>CCSD, if true=>CCD
            logical                 :: lambda = .false.
            logical                 :: noints = .false.
            logical                 :: eom = .false.
            logical                 :: eom_properties = .false.
            character(len=2)        :: eom_flavour = "  "
            logical                 :: beta_occ = .false.
            logical                 :: beta_vir = .false.

            integer                 :: exa_blocksize = 75
            integer                 :: talsh_buff = 50          ! (GB) total memory available for talsh
            integer                 :: moint_scheme = 3
            integer                 :: ncycles = 30
            integer                 :: print_level = 0
            integer                 :: nkr_occ, nkr_vir         ! number of Kramer pairs
            integer                 :: nocc, nvir               ! number of spinors

            integer, allocatable    :: mokr_occ(:), mokr_vir(:) ! list of active Kramer pairs
            integer, allocatable    :: mo_occ(:), mo_vir(:)     ! list of active spinors
            
            real(REALD)             :: t_econv = 1.0D-9         ! convergence criterium
            real(REALD)             :: level_shift = 0.D0       ! level shift
            real(REALD)             :: t_cholesky = 1.0D-9      ! threshold for cholesky
        end type exacc_input

!       All two electron integral tensors
        type, public :: exatns_intg_tens_t
         type(tens_rcrsv_t) :: oooo
         type(tens_rcrsv_t) :: ooov
         type(tens_rcrsv_t) :: oovv, vovo
         type(tens_rcrsv_t) :: vovv
         type(tens_rcrsv_t) :: vvvv
        end type exatns_intg_tens_t

!       All cholesky tensors
        type, public :: exatns_chol_tens_t
         type(tens_rcrsv_t) :: oou
         type(tens_rcrsv_t) :: ovu
         type(tens_rcrsv_t) :: vou
         type(tens_rcrsv_t) :: vvu
        end type exatns_chol_tens_t

!       Common blocks between t and lambda
        type, public :: exatns_comm_tens_t
         type(tens_rcrsv_t) :: foo, fov, fvv
         type(tens_rcrsv_t) :: goo, gvv
         type(tens_rcrsv_t) :: hoo, hov, hvv
         type(tens_rcrsv_t) :: a_int, h_int
         type(tens_rcrsv_t) :: t1, t2, tau
        end type exatns_comm_tens_t

        type, public :: exatns_lambda_tens_t
         type(tens_rcrsv_t) :: fbar_oo, fbar_ov, fbar_vv
         type(tens_rcrsv_t) :: w_oooo
         type(tens_rcrsv_t) :: w_ovoo, w_ooov
         type(tens_rcrsv_t) :: w_vovo, wbar_vovo
         type(tens_rcrsv_t) :: w_vvvo, w_vovv
         type(tens_rcrsv_t) :: w_vvvv
         type(tens_rcrsv_t) :: goo, gvv
        end type exatns_lambda_tens_t

        type, public :: exatns_dm_tens_t
         type(tens_rcrsv_t) :: prime_oo, symm_oo
         type(tens_rcrsv_t) :: prime_ov, prime_vo, symm_ov, symm_vo
         type(tens_rcrsv_t) :: prime_vv, symm_vv
        end type exatns_dm_tens_t

        type, public :: space_dims
         integer       :: occ_space_id, vir_space_id
         integer(INTL) :: occ_space_root,vir_space_root
        end type space_dims

!       All two electron integral tensors
        type, public :: talsh_intg_tens_t
         type(talsh_tens_t) :: oooo
         type(talsh_tens_t) :: ooov
         type(talsh_tens_t) :: oovv, vovo
         type(talsh_tens_t) :: vovv
         type(talsh_tens_t) :: vvvv
        end type talsh_intg_tens_t

!       Common blocks between t and lambda
        type, public :: talsh_comm_tens_t
         type(talsh_tens_t) :: foo, fov, fvv
         type(talsh_tens_t) :: goo, gvv
         type(talsh_tens_t) :: hoo, hov, hvv
         type(talsh_tens_t) :: a_int, h_int
         type(talsh_tens_t) :: t1, t2, tau
        end type talsh_comm_tens_t

        type, public :: talsh_lambda_tens_t
         type(talsh_tens_t) :: fbar_oo, fbar_ov, fbar_vv
         type(talsh_tens_t) :: w_oooo
         type(talsh_tens_t) :: w_ovoo, w_ooov
         type(talsh_tens_t) :: w_vovo, wbar_vovo
         type(talsh_tens_t) :: w_vvvo, w_vovv
         type(talsh_tens_t) :: w_vvvv
         type(talsh_tens_t) :: goo, gvv
        end type talsh_lambda_tens_t

        type, public :: talsh_dm_tens_t
         type(talsh_tens_t) :: prime_oo, symm_oo
         type(talsh_tens_t) :: prime_ov, prime_vo, symm_ov, symm_vo
         type(talsh_tens_t) :: prime_vv, symm_vv
        end type talsh_dm_tens_t

! intermediates for eom calculations. here we should only have things that do not change and are independent of eom flavor
!
        type, extends(talsh_lambda_tens_t), public :: talsh_intermediates_tens_t
           type(talsh_tens_t) :: t1, t2  ! t amplitudes used to construct the Hbar matrix elements, right now to be copied from the contents of comm_t, 
                                         ! but could be pointers to that storage as well. but the point of keeping them here is to keep everything 
                                         ! in one placewe make a copy of what may be pointers to what is  
           procedure(sigma_vector), pointer, nopass :: sigma_rhs_p => NULL()
           procedure(sigma_vector), pointer, nopass :: sigma_lhs_p => NULL()

        end type talsh_intermediates_tens_t

        abstract interface       
           subroutine sigma_vector(exa_input, int_t, interm_t, r1, r2, s1, s2)
           import :: exacc_input
           import :: talsh_intg_tens_t
           import :: talsh_intermediates_tens_t
           import :: talsh_tens_t

              type(exacc_input),     intent(in   ) :: exa_input
              type(talsh_intg_tens_t), intent(inout) :: int_t
              type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
              type(talsh_tens_t),    intent(inout) :: r1
              type(talsh_tens_t),    intent(inout) :: r2
              type(talsh_tens_t),    intent(inout) :: s1
              type(talsh_tens_t),    intent(inout) :: s2
            end subroutine  
        end interface


! this is a generic type for projectors, could be used e.g in eom for cvs, or in the ground state to project out orbitals in solving the cc equations
        type, public :: talsh_projectors_tens_t
         type(talsh_tens_t) :: p1, p2
        end type talsh_projectors_tens_t

!       One electron integrals
        type, public  :: one_el_t
         integer     :: n_spinor=-1          ! number of spinors, used for consistency checks
         real(8)     :: e_core               ! energy of the core electrons plus nuclear repulsion energy
         complex(8), pointer :: h_core(:,:)  ! one-body part of Hamiltonian (kinetic energy + nuclear attraction + core shielding)
        end type one_el_t

       end module exacorr_datatypes
