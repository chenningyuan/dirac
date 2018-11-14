!This module defines the data types used in the ExaCorr module of DIRAC
       module exacorr_datatypes

        use talsh, only     : talsh_tens_t
        use exatensor, only : tens_rcrsv_t

        implicit none
        private
!PARAMETERS:
!Generic:
        integer, parameter, private:: REALD=8 !default real for this module
!DERIVED TYPES:
!Basis function info:
        type, public:: basis_func_info_t
         integer:: orb_momentum=-1       !orbital momentum: (0,1,2,3,...)
         integer:: atom_number=-1        !atom sequential number [1..MAX] (if centered on a nucleus of some atom)
         integer:: atom_element=0        !atomic element number (if centered on a nucleus of some atom)
         integer:: n_primitives=0        !number of primitives in the contraction
         real(REALD):: coord(1:3)        !coordinate of the basis function center (typically a nucleus)
         real(REALD), pointer :: exponent(:)
         real(REALD), pointer :: coefficient(:)
        end type basis_func_info_t
!Basis set info:
        type, public:: basis_set_info_t
         integer:: nshells=0             !number of shells in the basis set
         integer:: nao=0                 !number of basis functions
         type(basis_func_info_t), pointer:: gtos(:)
         integer, pointer :: shell_indices(:)
        end type basis_set_info_t

!Input parameter
        type, public :: exacc_input
            logical                 :: talsh = .true.           !If false, use ExaTensor
            logical                 :: ccd = .false.            !If false=>CCSD, if true=>CCD
            integer                 :: nkr_occ, nkr_vir         !Number of Kramer pairs
            integer, allocatable    :: mokr_occ(:), mokr_vir(:) !List of active Kramer pairs
            real(REALD)             :: t_econv = 1.0D-13        !Convergence criterium
        end type exacc_input

!All two electron integral tensors
        type, public :: talsh_intg_tens
         type(talsh_tens_t) :: oooo
         type(talsh_tens_t) :: oovv, vovo
         type(talsh_tens_t) :: ooov, vovv
         type(talsh_tens_t) :: vvvv
        end type talsh_intg_tens

        type, public :: exatns_intg_tens
         type(tens_rcrsv_t) :: oooo
         type(tens_rcrsv_t) :: oovv, vovo
         type(tens_rcrsv_t) :: ooov, vovv
         type(tens_rcrsv_t) :: vvvv
        end type exatns_intg_tens

       end module exacorr_datatypes
