!This module defines a data type used elsewhere in DIRAC in case ExaCorr is not used.
       module exacorr_datatypes

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

!       Input parameter
        type, public :: exacc_input
            logical                 :: talsh = .true.           ! if false, use ExaTensor
            logical                 :: ccd = .false.            ! if false=>CCSD, if true=>CCD
            logical                 :: cc2 = .false.            !L if true CC2
            logical                 :: mp2no = .false.          !L if true mp2 natural orbital
            logical                 :: lambda = .false.
            logical                 :: beta_occ = .false.
            logical                 :: beta_vir = .false.
            logical                 :: do_triples = .true.
            logical                 :: do_no = .false.          ! if true, generate natural orbital

            integer                 :: exa_blocksize = 75
            integer                 :: talsh_buff = 50          ! (GB) total memory available for talsh
            integer                 :: moint_scheme = 4
            integer                 :: ncycles = 30
            integer                 :: print_level = 0
            integer                 :: nkr_occ, nkr_vir         ! number of Kramer pairs
            integer                 :: nocc, nvir               ! number of spinors
            integer                 :: tripl_block(1:2) = -1    ! for running triples
            integer                 :: nff(1:2) = 0             ! num. prop, num. field strengths
            
            integer, allocatable    :: mokr_occ(:), mokr_vir(:) ! list of active Kramer pairs
            integer, allocatable    :: mo_occ(:), mo_vir(:)     ! list of active spinors

            real(REALD)             :: t_econv = 1.0D-9         ! convergence criterium
            real(REALD)             :: level_shift = 0.D0       ! level shift
            real(REALD)             :: t_cholesky = 1.0D-9      ! threshold for cholesky
            real(REALD)             :: t_mp2no = 1.0D-5         ! threshold for mp2 natural orbital occupation number

            complex(8), allocatable   :: ff(:,:)                  ! field strengths
            character(8), allocatable :: ff_names(:)            ! property labels
        end type exacc_input

!       One electron integrals
        type, public  :: one_el_t
         integer     :: n_spinor=-1           ! number of spinors, used for consistency checks
         integer     :: n_prop=0              ! number of property matrices that were read
         real(8)     :: e_core                ! energy of the core electrons plus nuclear repulsion energy
         complex(8), pointer :: h_core(:,:)   ! one-body part of Hamiltonian (kinetic energy + nuclear attraction + core shielding)
         complex(8), pointer :: h_prop(:,:,:) ! property integrals to enable adding a finite field perturbation
         character*8,pointer :: property_labels(:)  ! strings that should correspond to the label as written to MDPROP (or other interface files)
        end type one_el_t

       end module exacorr_datatypes
