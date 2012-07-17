!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module adc_cfg

! This module contains all input set by the user. Used for RELADC

  implicit none
! integer, dimension(1:16), public :: relcc_nelec   = 0  
! integer, dimension(1:16), public :: relcc_nelec_f1= 0 
! integer, dimension(1:16), public :: relcc_nelec_f2= 0 
! integer, dimension(1:16), public :: relcc_nfroz   = 0
! save

! General control flags
! logical, public :: nelec_input                   = .false.
! logical, public :: relcc_debug                   = .false.
! integer, public :: relcc_max_dimension_diis      = 8
! integer, public :: relcc_max_iterations          = 30
! integer, public :: relcc_desired_convergence     = 12
!
! Information relevant for the integral sorting and setup stage
! logical, public :: relcc_breit                    = .false.
! logical, public :: relcc_add_finite_field         = .false.
! character(10), public :: relcc_fs_interface
! character( 8), public :: relcc_ff_property_names

! Information relevant for the energy calculation
! logical, public :: relcc_do_mp2                  = .true.
! logical, public :: relcc_do_singles              = .false.
! logical, public :: relcc_do_doubles              = .false.

end module
