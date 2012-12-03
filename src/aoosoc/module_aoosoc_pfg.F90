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

module module_aoosoc_pfg

  implicit none

  save

! parameters

  integer, parameter :: aoo_funit                   =  25    ! file unit for X2CMAT

  logical            :: aoo_mdirac                             =  .false. ! .true. if MDIRAC in DIRAC is active

! real(8) block
  real(8), public    :: aoo_cspeed                             = -1.0d0   ! speed of light

  integer            :: dim_pshell(2)                          =  0
  integer            :: dim_eshell(2)                          =  0
  integer            :: dim_oshell(2)                          =  0
  integer            :: nzt_aoo                                =  0
  integer            :: nr_quat                                =  0
  integer            :: nr_fsym                                =  0
  integer            :: mfsym                                  =  0
  integer            :: nrows                                  =  0
  integer            :: ncols                                  =  0
  integer            :: nr_ao_total_aoo                        =  0
  integer            :: nr_ao_large_aoo                        =  0
  integer            :: nr_ao_small_aoo                        =  0
  integer            :: n2bastq_dim_aoo                        =  0
  integer            :: nr_cmo_q                               =  0
  integer            :: nr_ao_all(2)                           =  0
  integer            :: nr_ao_l(2)                             =  0
  integer            :: nr_ao_s(2)                             =  0
  integer            :: ioff_aomat_x(2,2)                      =  0
  integer            :: aoo_cb_pq_to_uq(4, 0:7)                =  0
  integer            :: aoo_cb_uq_to_pq(4, 0:7)                =  0
  integer            :: aoo_bs_to_fs(0:7, 2)                   =  0
  integer            :: aoo_num_nuclei_i                       =  0
  integer            :: aoo_num_nuclei_d                       =  0
  integer            :: aoo_is_nuclei                          =  0
  integer            :: aoo_pointer_quat(0:7, 1:2)             =  0
  integer            :: aoo_bs_irrep_mat(1:4, 0:7)             =  0
  integer            :: aoo_iqmult_trip_q(1:4, 1:4, 1:4)       =  0
  integer            :: aoo_pointer_quat_op(0:7)               =  0
  integer            :: nr_2e_fock_matrices                    =  0
  integer            :: aoo_intflg                             =  0
  integer            :: aoo_mxatom                             =  0
  integer            :: aoo_mxcent                             =  0
  integer            :: aoo_nontyp                             =  0
  integer, allocatable :: aoo_nont(:)
  integer, allocatable :: aoo_nucdeg(:)

  !> definition of "atom" data type 
  type type_atom
        real(8) :: charge 
        real(8) :: coordinate_x
        real(8) :: coordinate_y
        real(8) :: coordinate_z
  end  type
  type(type_atom), allocatable :: atom(:)
  !> end


end module
