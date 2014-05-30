!dirac_copyright_start
!      Copyright (c) 2013 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

module pcmmod_cfg
   
  use, intrinsic :: iso_c_binding

  implicit none

  public host_input

! logical block
  logical, public, save :: pcmmod_host_provides_input = .false.
  logical, public, save :: pcmmod_skipss   =  .false.
  logical, public, save :: pcmmod_separate =  .false.
  logical, public, save :: pcmmod_dospf    =  .false.
  logical, public, save :: pcmmod_skipoit  =  .false.
! printlevel
  integer, public, save :: pcmmod_print    = 0
! *PCMSOL section
! cavity specification *PCMSOL section
  character(kind=c_char, len=8), public, save :: pcmmod_cavity_type = 'gepol  '//c_null_char
  integer(c_int), public, save  :: pcmmod_patch_level = 2
  real(c_double), public, save  :: pcmmod_coarsity = 0.5
  real(c_double), public, save  :: pcmmod_cavity_area = 0.3
  real(c_double), public, save  :: pcmmod_min_distance = 0.1
  integer(c_int), public, save  :: pcmmod_der_order = 4
  logical(c_bool), public, save  :: pcmmod_scaling = .true. 
  character(kind=c_char, len=8), public, save :: pcmmod_radii_set = 'bondi  '//c_null_char
  character(kind=c_char, len=20), public, save :: pcmmod_restart_name = 'cavity.npz         '//c_null_char
! solver specification *PCMSOL section
  character(kind=c_char, len=7), public, save :: pcmmod_solver_type = 'iefpcm'//c_null_char
  character(kind=c_char, len=16), public, save :: pcmmod_solvent = 'water          '//c_null_char
! green specification *PCMSOL section
  integer(c_int), public, save :: pcmmod_green_inside_type = 0
  integer(c_int), public, save :: pcmmod_green_outside_type = 0

  type, bind(c) :: cavityInput
          character(kind=c_char) :: cavity_type(8)
          integer(c_int) :: patch_level
          real(c_double) :: coarsity
          real(c_double) :: area
          real(c_double) :: min_distance
          integer(c_int) :: der_order
          logical(c_bool) :: scaling
          character(kind=c_char) :: radii_set(8)
          character(kind=c_char) :: restart_name(20)
  end type
  
  type, bind(c) :: solverInput
          character(kind=c_char) :: solver_type(7)
          character(kind=c_char) :: solvent(16)
  end type
  
  type, bind(c) :: greenInput
          integer(c_int) :: inside_type
          integer(c_int) :: outside_type
  end type
  
  contains

     subroutine host_input(cavity, solver, green) bind(c, name='host_input')
! Performs syntactic checks on PCMSolver input and fills the data structures
! holding input data
                                                           
     type(cavityInput) :: cavity
     type(solverInput) :: solver
     type(greenInput)  :: green

     integer :: i

     do i = 1, size(cavity%cavity_type)
        cavity%cavity_type(i)  = pcmmod_cavity_type(i:i)
     end do
     cavity%patch_level  = pcmmod_patch_level
     cavity%coarsity     = pcmmod_coarsity
     cavity%area         = pcmmod_cavity_area
     cavity%min_distance = pcmmod_min_distance
     cavity%der_order    = pcmmod_der_order
     cavity%scaling      = pcmmod_scaling
     do i = 1, size(cavity%radii_set)
        cavity%radii_set(i)    = pcmmod_radii_set(i:i)
     end do
     do i = 1, size(cavity%restart_name)
        cavity%restart_name(i) = pcmmod_restart_name(i:i)
     end do

     do i = 1, size(solver%solver_type)
        solver%solver_type(i) = pcmmod_solver_type(i:i)
     end do
     do i = 1, size(solver%solvent)
        solver%solvent(i)     = pcmmod_solvent(i:i)
     end do

     green%inside_type  = pcmmod_green_inside_type 
     green%outside_type = pcmmod_green_outside_type
                                                           
     end subroutine host_input

end module
