!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

module xamfi_internal_parameters

  implicit none

  save

! parameters

  integer, parameter :: aoo_funit                   =  25    ! file unit for X2CMAT
  integer, parameter :: aoo_mxopen                  =  10    ! same as mxopen

  logical            :: aoo_mdirac                             =  .false. ! .true. if MDIRAC in DIRAC is active

! real(8) block
  real(8), public    :: aoo_cspeed                             = -1.0d0   ! speed of light
  real(8)            :: aoo_dfopen(0:aoo_mxopen)               =  0.0d0

  integer            :: dim_pshell(1:2)                        =  0
  integer            :: dim_eshell(1:2)                        =  0
  integer            :: dim_oshell(1:2)                        =  0
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
  integer            :: n2bbasxq_dim_aoo                       =  0
  integer            :: nr_cmo_q                               =  0
  integer            :: aoo_n2tmt                              =  0
  integer            :: nr_ao_all(1:2)                         =  0
  integer            :: nr_ao_l(1:2)                           =  0
  integer            :: nr_ao_s(1:2)                           =  0
  integer            :: ioff_aomat_x(1:2,1:2)                  =  0
  integer            :: aoo_cb_pq_to_uq(1:4, 0:7)              =  0
  integer            :: aoo_cb_uq_to_pq(1:4, 0:7)              =  0
  integer            :: aoo_bs_to_fs(0:7, 1:2)                 =  0
  integer            :: aoo_num_nuclei_i                       =  0
  integer            :: aoo_num_nuclei_d                       =  0
  integer            :: aoo_is_nuclei                          =  0
  integer            :: aoo_pointer_quat(0:7, 1:2)             =  0
  integer            :: aoo_bs_irrep_mat(1:4, 0:7)             =  0
  integer            :: aoo_iqmult_trip_q(1:4, 1:4, 1:4)       =  0
  integer            :: aoo_pointer_quat_op(0:7)               =  0
  integer            :: aoo_ihqmat(1:4,-1:1)                   =  0
  integer            :: nr_2e_fock_matrices                    =  0
  integer            :: aoo_intflg                             =  0
  integer            :: aoo_mxatom                             =  0
  integer            :: aoo_mxcent                             =  0
  integer            :: aoo_nontyp                             =  0
  integer            :: aoo_nopen                              =  0
  integer            :: lutmat_aoo                             =  0
  integer            :: nr_tmo(2)                              =  0
  integer            :: ioff_tmt(2)                            =  0
  integer            :: ioff_tmot(2)                           =  0
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

  !> cb exchange
  integer            :: nr_quat_save            = 0
  integer            :: nzt_aoo_save            = 0
  integer            :: nr_ao_total_aoo_save    = 0
  integer            :: nrows_save              = 0
  integer            :: ncols_save              = 0
  integer            :: nr_ao_large_aoo_save    = 0
  integer            :: nr_ao_small_aoo_save    = 0
  integer            :: n2bastq_dim_aoo_save    = 0
  integer            :: n2bbasxq_dim_aoo_save   = 0
  integer            :: nr_cmo_q_save           = 0
  integer            :: aoo_n2tmt_save          = 0

  integer            :: dim_pshell_save(1:2)    = 0
  integer            :: dim_eshell_save(1:2)    = 0
  integer            :: dim_oshell_save(1:2)    = 0
  integer            :: nr_ao_all_save(1:2)     = 0
  integer            :: nr_ao_l_save(1:2)       = 0
  integer            :: nr_ao_s_save(1:2)       = 0

  integer            :: nr_tmo_save(1:2)        = 0
  integer            :: ioff_tmt_save(1:2)      = 0
  integer            :: ioff_tmot_save(1:2)     = 0

  integer            :: ioff_aomat_x_save(1:2,1:2) = 0


end module
