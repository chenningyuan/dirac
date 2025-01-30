!This module defines the data types used in the ExaCorr module of DIRAC
       module exacorr_datatypes

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        use exatensor
#endif
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

        type, public :: eom_config_t
           integer             :: nroots
           character(len=8)    :: flavor
           logical             :: do_cvs_prj
           logical             :: do_expval
           logical             :: do_linrsp
           logical             :: do_quarsp
           logical             :: do_cubrsp
           integer             :: rhs_expval_states
           integer             :: lhs_expval_states
           character(len=8)    :: properties
           complex(8), pointer :: rsp_pert_freqs(:)
        end type eom_config_t

!       Input parameter
        type, public :: exacc_input
            logical                 :: talsh = .true.           ! if false, use ExaTensor
            logical                 :: ccd = .false.            ! if false=>CCSD, if true=>CCD
            logical                 :: cc2 = .false.            !L if true CC2
            logical                 :: do_cis = .false.         ! enables CIS, and branch off away from the CC driver 
            logical                 :: lambda = .false.
            logical                 :: beta_occ = .false.
            logical                 :: beta_vir = .false.
            logical                 :: do_iterative   = .true.  ! this provides a way to switch off iterative methods, e.g. if we want to do only MP2 
            logical                 :: do_perturbative = .true. ! this generalizes the do_triples
            logical                 :: do_triples = .true.
            logical                 :: tripl_lap = .false.
            logical                 :: do_no = .false.          ! if true, generate natural orbital
            logical                 :: eom_lr = .false.
            
            logical                 :: do_1rdm_gs=.false.
            logical                 :: do_2rdm_gs=.false.

            logical                 :: auto_kernels = .false.   ! if true, use automatically generated kernels for t/l equations
            integer                 :: auto_scheme  = 2         ! 1: autogen;  2: tenpi

            integer                 :: exa_blocksize = 75
            integer                 :: talsh_buff = 50          ! (GB) total memory available for talsh
            integer                 :: moint_scheme = 4
            integer                 :: mp2no = 0                ! > 0 mp2 natural orbital, < 0 mp2lap
            integer                 :: ncycles = 30
            integer                 :: print_level = 0
            integer                 :: nkr_occ, nkr_vir         ! number of Kramer pairs
            integer                 :: nocc, nvir               ! number of spinors
            integer                 :: iterative_level = 2      ! at which excitation level iterative amplitudes/lagrange mult. will be calculated: 1, 2, 3, 4...
            integer                 :: iterative_scheme = 21    ! scheme/variant for the iterative determination of amplitudes 
                                                                ! 10: CCS (not currently implemented, but logically should follow as this)
                                                                ! 20: CCD      ; 21: CCSD       ; 22: CC2        ; 
                                                                ! 30: CCSDT    ; 31: CCSDT-1a   ; 32: CCSDT-1b   ; 33: CCSDT-2   ; 34: CCSDT-3   ; 35:  CC3
                                                                ! 40: CCSDTQ   ; 41: CCSDTQ-1a  ; 42: CCSDTQ-1b  ; 43: CCSDTQ-2  ; 44: CCSDTQ-3  ; 45:  CC4
                                                                ! 50: CCSDTQP  ; 51: CCSDTQP-1a ; 52: CCSDTQP-1b ; 53: CCSDTQP-2 ; 54: CCSDTQP-3 ; 55:  CC5
            integer                 :: perturbative_level = 3   ! at which excitation level pertubative corrections will be calculated: 2, 3, 4...
            integer                 :: perturbative_scheme = 30 ! scheme/variant for the iterative determination of energies/amplitudes/... 
                                                                ! 20: MP2 (MP2 is always performed)
                                                                ! 30: CCSD(T)

!  aspg, TODO: harmonize this tripl_block with the more general setup of excitation rank
            integer                 :: tripl_block(1:2) = -1    ! for running triples
            integer                 :: nff(1:2) = 0             ! num. prop, num. field strengths
            integer                 :: n_laplace = 10           ! number of laplace points to use
            integer                 :: N_eom_trial = 1          ! number of trial vectors in CC linear response
            integer                 :: wf_type = 2              ! wave-function type
            integer                 :: N_Prop_A, N_Prop_B    ! number of property
            integer                 :: N_Prop_C, N_Prop_D    ! number of property
            integer                 :: N_omega_lr = 1               ! number of frequency in cclr
            integer                 :: N_omega_qr1 = 1               ! number of frequency1 in ccqr
            integer                 :: N_omega_qr2 = 1               ! number of frequency2 in ccqr
            integer                 :: N_Time_anti = 0          ! Number of time-antisymmetric operators
           
            
            integer, allocatable    :: mokr_occ(:), mokr_vir(:) ! list of active Kramer pairs
            integer, allocatable    :: mo_occ(:), mo_vir(:)     ! list of active spinors
            integer, allocatable    :: Prop_A(:), Prop_B(:), Prop_C(:), Prop_D(:) ! index of Property
            INTEGER                 ::  TPA_state             ! index of final excited state in eom-tpa
            real(REALD), allocatable :: w_lr_re(:)            ! real frequency in CC linear response
            real(REALD), allocatable :: w_lr_im(:)            ! imaginary frequency in CC response

            real(REALD), allocatable :: w_qr1_re(:)            ! 1st real frequency in CC quadratic response
            real(REALD), allocatable :: w_qr1_im(:)            ! 1st imaginary frequency in CC quadratic response
            real(REALD), allocatable :: w_qr2_re(:)            ! 2nd real frequency in CC quadratic response
            real(REALD), allocatable :: w_qr2_im(:)            ! 2nd imaginary frequency in CC quadratic response
           

            real(REALD)             :: t_econv = 1.0D-9         ! convergence criterium
            real(REALD)             :: level_shift = 0.D0       ! level shift
            real(REALD)             :: t_cholesky = 1.0D-9      ! threshold for cholesky / mp2lap
            real(REALD)             :: t_mp2no = 1.0D-5         ! threshold for mp2 natural orbital occupation number
            
           

            real(REALD), allocatable:: rew_energies(:,:) ! restricted excitation window, in terms of energy thresholds (a.u.)

! this just signals we should handle *EOMCC inputs under **EXACC
            logical                 :: eom = .false.
            integer                 :: current_eom_input = -1
        !     logical                 :: eom_lr = .false.
            logical                 :: eom_qr = .false.
            logical                 :: eom_tpa = .false.
! handle one or more *EOMCC blocks 
            type(eom_config_t), allocatable :: eom_input(:)

            complex(8), allocatable   :: ff(:,:)                  ! field strengths
            character(8), allocatable :: ff_names(:)            ! property labels
            
        end type exacc_input

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
!       All two electron integral tensors
        type, public :: exatns_intg_tens_t
         type(tens_rcrsv_t) :: oooo
         type(tens_rcrsv_t) :: ooov
         type(tens_rcrsv_t) :: oovv, vovo
         type(tens_rcrsv_t) :: vovv
         type(tens_rcrsv_t) :: vvvv
         type(tens_rcrsv_t) :: foo, fov, fvv
        end type exatns_intg_tens_t

!       Common blocks between t and lambda
        type, public :: exatns_comm_tens_t
         type(tens_rcrsv_t) :: foo, fov, fvv
         type(tens_rcrsv_t) :: goo, gvv
         type(tens_rcrsv_t) :: hoo, hov, hvv
         type(tens_rcrsv_t) :: a_int, h_int
         type(tens_rcrsv_t) :: t1, t2, t3, t4, t5, tau
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
! one-body density matrices
         type(tens_rcrsv_t) :: prime_oo, symm_oo
         type(tens_rcrsv_t) :: prime_ov, prime_vo, symm_ov, symm_vo
         type(tens_rcrsv_t) :: prime_vv, symm_vv
! two-body density matrices
         type(tens_rcrsv_t) :: prime_vvvv
         type(tens_rcrsv_t) :: prime_vvvo, prime_vovv, prime_vvov, prime_ovvv
         type(tens_rcrsv_t) :: prime_vvoo, prime_oovv
         type(tens_rcrsv_t) :: prime_ovoo, prime_ooov, prime_vooo, prime_oovo
         type(tens_rcrsv_t) :: prime_ovov, prime_voov, prime_vovo, prime_ovvo
         type(tens_rcrsv_t) :: prime_oooo
         type(tens_rcrsv_t) :: symm_vvvv
         type(tens_rcrsv_t) :: symm_vvvo,  symm_vovv,  symm_vvov,  symm_ovvv
         type(tens_rcrsv_t) :: symm_vvoo,  symm_oovv
         type(tens_rcrsv_t) :: symm_vooo,  symm_ooov,  symm_ovoo,  symm_oovo
         type(tens_rcrsv_t) :: symm_vovo,  symm_voov,  symm_ovvo,  symm_ovov
         type(tens_rcrsv_t) :: symm_oooo
        end type exatns_dm_tens_t


        type, public :: space_dims
         integer       :: occ_space_id, vir_space_id
         integer(INTL) :: occ_space_root,vir_space_root
        end type space_dims
#endif

!       All fixed one and two electron integral tensors (defines 2nd quantized Hamiltonian)
        type, public :: talsh_intg_tens_t
         type(talsh_tens_t) :: oooo
         type(talsh_tens_t) :: ooov
         type(talsh_tens_t) :: oovv, vovo
         type(talsh_tens_t) :: vovv
         type(talsh_tens_t) :: vvvv
         type(talsh_tens_t) :: foo, fov, fvv
        end type talsh_intg_tens_t

!       Common blocks between t and lambda
        type, public :: talsh_comm_tens_t
         type(talsh_tens_t) :: foo, fov, fvv ! include finite field perturbations, otherwise identical to f in talsh_intg_tens
         type(talsh_tens_t) :: goo, gvv
         type(talsh_tens_t) :: hoo, hov, hvv
         type(talsh_tens_t) :: a_int, h_int
         type(talsh_tens_t) :: t1, t2, t3, t4, t5, tau
        end type talsh_comm_tens_t

        type, public :: talsh_lambda_tens_t
         type(talsh_tens_t) :: fbar_oo, fbar_ov, fbar_vv
         type(talsh_tens_t) :: w_oooo
         type(talsh_tens_t) :: w_ovoo, w_ooov
         type(talsh_tens_t) :: w_vovo, wbar_vovo
         type(talsh_tens_t) :: w_vvvo, w_vovv
         type(talsh_tens_t) :: w_vvvv
         type(talsh_tens_t) :: goo, gvv
         type(talsh_tens_t) :: eta_ov, eta_oovv
        end type talsh_lambda_tens_t


        type, public :: talsh_orb_prop_integ_tens_t
         type(talsh_tens_t) :: prop_oo, prop_ov, prop_vv, prop_vo
        end type talsh_orb_prop_integ_tens_t

        
! intermediates for eom calculations. here we should only have things that do not change and are independent of eom flavor
!
        type, extends(talsh_lambda_tens_t), public :: talsh_intermediates_tens_t
           logical :: ccd = .false.! aspg: having this here is a hack, in reality the sigma routines should be created by composing with more primitive collections of contractions (so ccsd collects "ccs+ccd+whatever") and then specialized routines pick and choose which approximated contractions to use.
           integer :: nocc, nvir
           type(talsh_tens_t) :: t1, t2  ! t amplitudes used to construct the Hbar matrix elements, right now to be copied from the contents of comm_t, 
                                         ! but could be pointers to that storage as well. but the point of keeping them here is to keep everything 
                                         ! in one placewe make a copy of what may be pointers to what is  
           procedure(sigma_vector), pointer, nopass :: sigma_rhs_p => NULL()
           procedure(sigma_vector), pointer, nopass :: sigma_lhs_p => NULL()

        end type talsh_intermediates_tens_t

        abstract interface       
           subroutine sigma_vector(int_t, interm_t, r1, r2, s1, s2)
           import :: exacc_input
           import :: talsh_intg_tens_t
           import :: talsh_intermediates_tens_t
           import :: talsh_tens_t

              type(talsh_intg_tens_t), intent(inout) :: int_t
              type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
              type(talsh_tens_t),    intent(inout) :: r1
              type(talsh_tens_t),    intent(inout) :: r2
              type(talsh_tens_t),    intent(inout) :: s1
              type(talsh_tens_t),    intent(inout) :: s2
            end subroutine  
        end interface

        ! intermediates used in eom prop and cc response 
        type, public :: talsh_rsp_interm_tens_t
              type(talsh_tens_t) :: prop_bar_vv, prop_bar_oo, prop_bar_ov
              type(talsh_tens_t) :: Q_vvvo, Q_ovoo
              type(talsh_tens_t) :: eta_ov, eta_oovv
              type(talsh_tens_t) :: epsilon_prop_vo, epsilon_prop_vvoo
              type(talsh_tens_t) :: eta_prop_ov, eta_prop_oovv, eta_prop_ov_eom
              type(talsh_tens_t) :: O_vo
       end type talsh_rsp_interm_tens_t

! this is a generic type for projectors, could be used e.g in eom for cvs, or in the ground state to project out orbitals in solving the cc equations
        type, public :: talsh_projectors_tens_t
          integer, pointer :: p1_dims(:), p2_dims(:)        
          integer :: p1_rank, p2_rank                       
         type(talsh_tens_t) :: p1, p2
        end type talsh_projectors_tens_t


        type, public :: talsh_dm_tens_t
! one-body density matrices
         type(talsh_tens_t) :: prime_oo, symm_oo
         type(talsh_tens_t) :: prime_ov, prime_vo, symm_ov, symm_vo
         type(talsh_tens_t) :: prime_vv, symm_vv
! two-body density matrices
         type(talsh_tens_t) :: prime_vvvv
         type(talsh_tens_t) :: prime_vvvo, prime_vovv, prime_vvov, prime_ovvv
         type(talsh_tens_t) :: prime_vvoo, prime_oovv
         type(talsh_tens_t) :: prime_ovoo, prime_ooov, prime_vooo, prime_oovo
         type(talsh_tens_t) :: prime_ovov, prime_voov, prime_vovo, prime_ovvo
         type(talsh_tens_t) :: prime_oooo
         type(talsh_tens_t) :: symm_vvvv
         type(talsh_tens_t) :: symm_vvvo,  symm_vovv,  symm_vvov,  symm_ovvv
         type(talsh_tens_t) :: symm_vvoo,  symm_oovv
         type(talsh_tens_t) :: symm_vooo,  symm_ooov,  symm_ovoo,  symm_oovo
         type(talsh_tens_t) :: symm_vovo,  symm_voov,  symm_ovvo,  symm_ovov
         type(talsh_tens_t) :: symm_oooo
        end type talsh_dm_tens_t

        ! type, public, extends(talsh_comm_tens) :: talsh_eom_fixed_tens
        !     type(talsh_intg_tens)   :: int_t
        !     type(talsh_lambda_tens) :: lambda_t
        !     type(talsh_tens_t)      :: l2_tensor ! lambda amplitude
        ! end type talsh_eom_fixed_tens

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
