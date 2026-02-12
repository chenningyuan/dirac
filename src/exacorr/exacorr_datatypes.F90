!This module defines the tensorial data types used in the ExaCorr module of DIRAC
       module exacorr_datatypes

        use exatensor
        use exacorr_tensor_interface, only: tensor_t
        use exacorr_cfg, only: exacc_input
        use, intrinsic:: ISO_C_BINDING

        implicit none
        private

        integer, parameter, private:: REALD=8 ! default real for this module

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

!       All fixed one and two electron integral tensors (defines 2nd quantized Hamiltonian)
        type, public :: eti_intg_tens_t
         type(tensor_t) :: oooo
         type(tensor_t) :: ooov,ovoo
         type(tensor_t) :: oovv, vovo, ovov,ovvo
         type(tensor_t) :: vovv,ovvv,vvvo
         type(tensor_t) :: vvvv
         type(tensor_t) :: foo, fov, fvv
        end type eti_intg_tens_t

!       Common blocks between t and lambda
        type, public :: eti_comm_tens_t
         type(tensor_t) :: foo, fov, fvv ! include finite field perturbations, otherwise identical to f in eti_intg_tens
         type(tensor_t) :: goo, gvv
         type(tensor_t) :: hoo, hov, hvv
         type(tensor_t) :: a_int, h_int
         type(tensor_t) :: t1, t2, t3, t4, t5, tau
         type(tensor_t) :: l1, l2, l3, l4, l5
        end type eti_comm_tens_t

        type, public :: eti_lambda_tens_t
         type(tensor_t) :: fbar_oo, fbar_ov, fbar_vv
         type(tensor_t) :: w_oooo
         type(tensor_t) :: w_ovoo, w_ooov
         type(tensor_t) :: w_vovo, wbar_vovo
         type(tensor_t) :: w_vvvo, w_vovv
         type(tensor_t) :: w_vvvv
         type(tensor_t) :: goo, gvv
         type(tensor_t) :: eta_ov, eta_oovv
        end type eti_lambda_tens_t


        type, public :: eti_orb_prop_integ_tens_t
         type(tensor_t) :: prop_oo, prop_ov, prop_vv, prop_vo
         character*8 :: property_label
        end type eti_orb_prop_integ_tens_t


! intermediates for eom calculations. here we should only have things that do not change and are independent of eom flavor
!
        type, extends(eti_lambda_tens_t), public :: eti_intermediates_tens_t
           logical :: ccd = .false.! aspg: having this here is a hack, in reality the sigma routines should be created by composing with more primitive collections of contractions (so ccsd collects "ccs+ccd+whatever") and then specialized routines pick and choose which approximated contractions to use.
           integer :: nocc, nvir
           type(tensor_t) :: t1, t2  ! t amplitudes used to construct the Hbar matrix elements, right now to be copied from the contents of comm_t, 
                                         ! but could be pointers to that storage as well. but the point of keeping them here is to keep everything 
                                         ! in one placewe make a copy of what may be pointers to what is  
           procedure(sigma_vector), pointer, nopass :: sigma_rhs_p => NULL()
           procedure(sigma_vector), pointer, nopass :: sigma_lhs_p => NULL()

        end type eti_intermediates_tens_t

        abstract interface
           subroutine sigma_vector(int_t, interm_t, r1, r2, s1, s2)
           import :: exacc_input
           import :: eti_intg_tens_t
           import :: eti_intermediates_tens_t
           import :: tensor_t

              type(eti_intg_tens_t), intent(inout) :: int_t
              type(eti_intermediates_tens_t),  intent(inout) :: interm_t
              type(tensor_t),    intent(inout) :: r1
              type(tensor_t),    intent(inout) :: r2
              type(tensor_t),    intent(inout) :: s1
              type(tensor_t),    intent(inout) :: s2
            end subroutine  
        end interface

        ! intermediates used in eom prop and cc response 
        type, public :: eti_rsp_interm_tens_t
              type(tensor_t) :: prop_bar_vv, prop_bar_oo, prop_bar_ov
              type(tensor_t) :: prop_bar_oovv
              type(tensor_t) :: Q_vvvo, Q_ovoo
              type(tensor_t) :: eta_ov, eta_oovv
              type(tensor_t) :: epsilon_prop_vo, epsilon_prop_vvoo
              type(tensor_t) :: eta_prop_ov, eta_prop_oovv, eta_prop_ov_eom
              type(tensor_t) :: O_vo
              character*8    :: property_label
       end type eti_rsp_interm_tens_t

! this is a generic type for projectors, could be used e.g in eom for cvs, or in the ground state to project out orbitals in solving the cc equations
        type, public :: eti_projectors_tens_t
          integer, pointer :: p1_dims(:), p2_dims(:)        
          integer :: p1_rank, p2_rank                       
         type(tensor_t) :: p1, p2
        end type eti_projectors_tens_t


        type, public :: eti_dm_tens_t
! one-body density matrices
         type(tensor_t) :: prime_oo, symm_oo
         type(tensor_t) :: prime_ov, prime_vo, symm_ov, symm_vo
         type(tensor_t) :: prime_vv, symm_vv
         type(tensor_t) :: asym_oo, asym_vv, asym_ov, asym_vo
! two-body density matrices
         type(tensor_t) :: prime_vvvv
         type(tensor_t) :: prime_vvvo, prime_vovv, prime_vvov, prime_ovvv
         type(tensor_t) :: prime_vvoo, prime_oovv
         type(tensor_t) :: prime_ovoo, prime_ooov, prime_vooo, prime_oovo
         type(tensor_t) :: prime_ovov, prime_voov, prime_vovo, prime_ovvo
         type(tensor_t) :: prime_oooo
         type(tensor_t) :: symm_vvvv
         type(tensor_t) :: symm_vvvo,  symm_vovv,  symm_vvov,  symm_ovvv
         type(tensor_t) :: symm_vvoo,  symm_oovv
         type(tensor_t) :: symm_vooo,  symm_ooov,  symm_ovoo,  symm_oovo
         type(tensor_t) :: symm_vovo,  symm_voov,  symm_ovvo,  symm_ovov
         type(tensor_t) :: symm_oooo
        end type eti_dm_tens_t

        ! type, public, extends(eti_comm_tens) :: eti_eom_fixed_tens
        !     type(eti_intg_tens)   :: int_t
        !     type(eti_lambda_tens) :: lambda_t
        !     type(tensor_t)      :: l2_tensor ! lambda amplitude
        ! end type eti_eom_fixed_tens

       end module exacorr_datatypes
