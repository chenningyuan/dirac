! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Module:       Kramers unrestricted SCF module
!
! Description:  Public derived data types for the module 
!
! Contains:     
!
! Licensing:    dirac_copyright_start
!                     Copyright (c) 2010 by the authors of DIRAC.
!                     All Rights Reserved.
!                     
!                     This source code is part of the DIRAC program package.
!                     It is provided under a written license and may be used,
!                     copied, transmitted, or stored only in accordance to the
!                     conditions of that written license.
!                     
!                     In particular, no part of the source code or compiled modules may
!                     be distributed outside the research group of the license holder.
!                     This means also that persons (e.g. post-docs) leaving the research
!                     group of the license holder may not take any part of Dirac,
!                     including modified files, with him/her, unless that person has
!                     obtained his/her own license.
!                     
!                     For information on how to get a license, as well as the
!                     author list and the complete list of contributors to the
!                     DIRAC program, see: http://dirac.chem.vu.nl
!               dirac_copyright_end
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Logs:         
!
! Comments:     Note: do not change the order of components in derived data types
!                     unless you appropriately modify the intrinsic structure constructor 
!                     for the derived data type initialisation
!                            
! ------------------------------------------------------------------------------------

MODULE SCF_DATA_TYPES

  implicit none

  public


  !> definition of "scf" data type
  type type_scf 
       integer :: max_scf_iterations
       integer :: max_diis_dimension
       real(8) :: diis_error_treshold
       real(8) :: energy_change_treshold
  end type type_scf
  type (type_scf) :: scf


  !> definition of "element" data type
  type type_element 
       character(len=7) :: name 
       real(8)          :: mass  
       real(8)          :: number
       character(len=4) :: symbol 
  end type type_element 
  type (type_element)   :: element(1:109)


  !> definition of "atom" data type 
  type type_atom
       real(8) :: charge 
       real(8) :: coordinate_x
       real(8) :: coordinate_y
       real(8) :: coordinate_z
       real(8) :: gnu_exponent 
  end  type
  type(type_atom), allocatable :: atom(:)


  !> definition of "gto" data type
  type type_gto
       integer :: index
       integer :: offset
       integer :: lvalue
       integer :: origin
       integer :: sdegen 
       integer :: cdegen 
       real(8) :: norm(1)
       real(8) :: exponent(1)
       real(8) :: coefficient(1)
  end  type
  type(type_gto), allocatable, target :: gto(:)


  !> definition of "real_vector" data type
  type real_vector 
       real(8), pointer :: elms(:) 
  end  type


  !> definition of "real_matrix" data type 
  type real_matrix
       integer          :: nrow 
       integer          :: ncol  
       real(8), pointer :: elms(:,:) 
  end  type


  !> definition of "complex_pointer" data type 
  type complex_pointer
       complex(8), pointer :: p => NULL()
  end  type


  !> definition of "complex_matrix" data type 
  type complex_matrix
       integer                 :: nrow 
       integer                 :: ncol  
       type(complex_pointer)   :: info(4,4)  
       complex(8), allocatable :: elms(:,:,:) 
  end  type

END MODULE
