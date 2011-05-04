!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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

module xc_driver

!  written by T. Helgaker
!  DIRAC adaption by T. Saue

#ifdef PRG_DIRAC
!  these dependencies should be removed
   use fde_mod
#endif
   use num_grid_cfg
   use basis_info
   use occupation_info
   use dft_cfg
   use xc_derv
   use num_grid_io
   use ao_eval
   use density_eval
   use xc_blas3
   use xc_max_block_length
#ifdef MOD_UNRELEASED
   use overlap_diagnostic
   use xc_geometric
   use xc_london_c1
#endif
#ifdef MOD_OPENRSP
   use xc_response
#endif
   use xc_mpi
#ifdef MOD_XCFUN
   use xcfun
#endif

   implicit none

   public xc_integrator
   public get_xc_energy
   public get_xc_mat_energy

   private

   type omega_prefactor
      real(8), allocatable :: n(:)
      real(8), allocatable :: s(:, :)
      real(8), allocatable :: gn(:, :)
      real(8), allocatable :: gs(:, :, :)
      real(8), allocatable :: tn(:)
      real(8), allocatable :: ts(:, :)
!gosia, temporarily for sdft+lao
      real(8), allocatable :: s_lao(:, :)
      real(8), allocatable :: gs_lao(:, :, :)
   end type

   integer, parameter :: bllen = 1 !later more, now slowdft
   integer            :: block_length

   save

   integer, parameter :: file_unit      = 6
   integer, parameter :: file_unit_grid = 44

   real(8), allocatable :: dmat_saop_outer(:)
   real(8), allocatable :: dmat_saop_resp(:)
   real(8), allocatable :: rx(:)
   real(8), allocatable :: ry(:)
   real(8), allocatable :: rz(:)
   real(8), allocatable :: rw(:)
   real(8), allocatable :: fde_prop1(:)
   real(8), allocatable :: fde_prop2(:)
   real(8), allocatable :: ao(:)
   real(8), allocatable :: buffer(:)
   real(8), allocatable :: dmat_sorted(:, :)

   real(8) :: tau_integral
   real(8) :: nr_electrons
   real(8) :: xc_energy
   real(8) :: xc_mat_energy

   logical :: use_gga_qr = .false. !quaternion real
   logical :: use_gga_qi = .false. !quaternion imaginary
   logical :: use_gga    = .false.

   logical :: parallel_xc

!             number of points in current batch
   integer :: nr_points_batch

!             number of points currently looped over on this processor
   integer :: nr_points_on_this_proc

!             total nr of points (sum of nr_points_batch)
   integer :: nr_points_total

!             total nr of points that were actually used
!             we skip points with very small density
!             typically smaller than nr_points_total
   integer :: nr_points_used

   integer :: max_ao_g_order
   integer :: max_ao_m_order
   integer :: max_fun_derv

   real(8), external :: second
   real(8)           :: time_integration, time_integration_start
   real(8)           :: time_ao,          time_ao_start
   real(8)           :: time_matrix_dist, time_matrix_dist_start
   real(8)           :: time_density,     time_density_start
   real(8)           :: time_derv,        time_derv_start

   integer                      :: mat_dim
   integer                      :: nz

   real(8), pointer             :: dmat_0(:, :)
   real(8), allocatable, target :: dmat_0_container(:, :)

   real(8), allocatable         :: tmat_1(:, :)
   real(8), allocatable         :: tmat_2(:, :)

   integer                      :: nr_dmat
   real(8), pointer             :: dmat(:, :, :, :)
   real(8), allocatable, target :: dmat_container(:, :, :, :)

   integer                      :: nr_fmat
   real(8), pointer             :: fmat(:, :, :, :)
   real(8), allocatable, target :: fmat_container(:, :, :, :)

   integer, pointer             :: dmat_pg_sym(:)
   integer, allocatable, target :: dmat_pg_sym_container(:)
   integer, pointer             :: dmat_ih_sym(:)
   integer, allocatable, target :: dmat_ih_sym_container(:)
   integer, pointer             :: fmat_pg_sym(:)
   integer, allocatable, target :: fmat_pg_sym_container(:)

   integer                      :: nr_atoms
   real(8), pointer             :: property_gradient(:)
   real(8), allocatable, target :: property_gradient_container(:)

#ifdef MOD_UNRELEASED
   real(8), allocatable :: fde_froz_n(:)
   real(8), allocatable :: fde_froz_ng(:,:)
   real(8), allocatable :: fde_elpot(:)
#endif

   logical :: do_potential
   integer :: response_order_mo
   integer :: response_order_ao
   logical :: do_london_rhs_direct
   logical :: do_london_rhs_ro
   logical :: do_geo
   logical :: do_geo_0
   logical :: do_geo_1
   logical :: do_geo_2
   logical :: do_london_ks
   logical :: do_london_lr
   logical :: do_overlap_diagnostic
   logical :: do_fde_update_vemb
   logical :: do_fde_response

contains

   subroutine xc_integrator(xc_mat_dim,               &
                            xc_nz,                    &
                            xc_dmat_0,                &
                            xc_nr_dmat,               &
                            xc_nr_fmat,               &
                            xc_dmat,                  &
                            xc_fmat,                  &
                            xc_dmat_ih_sym,           &
                            xc_dmat_pg_sym,           &
                            xc_fmat_pg_sym,           &
                            xc_nr_atoms,              &
                            xc_property_gradient,     &
                            xc_do_potential,          &
                            xc_do_london_rhs_direct,  &
                            xc_do_london_rhs_ro,      &
                            xc_do_geo_0,              &
                            xc_do_geo_1,              &
                            xc_do_geo_2,              &
                            xc_do_london_ks,          &
                            xc_do_london_lr,          &
                            xc_do_overlap_diagnostic, &
                            xc_do_fde_update_vemb,    &
                            xc_do_fde_response,       &
                            xc_response_order_mo,     &
                            xc_response_order_ao)

! response_order_mo/_ao:  0 - do nothing (no response)
!                         1 - ks lr
!                         2 - ks qr
!                         *_mo - old code, *_ao - openrsp code
!
! do_potential            = .true.  - do energy and xc potential
!
! contributions to differentiated kohn-sham matrix:
! do_london_rhs_direct       = .true.  - xc contributions that require mag dervs of ao-s, "direct" (if london)
! do_london_rhs_ro           = .true.  - xc contributions to reotrhonormalization terms (if london)

!     --------------------------------------------------------------------------
      integer,                   intent(in)    :: xc_mat_dim
      integer,                   intent(in)    :: xc_nz
      real(8),           target, intent(in)    :: xc_dmat_0(xc_mat_dim, xc_mat_dim)
      integer,                   intent(in)    :: xc_nr_dmat
      integer,                   intent(in)    :: xc_nr_fmat
      real(8), optional, target, intent(in)    :: xc_dmat(xc_mat_dim, xc_mat_dim, xc_nz, xc_nr_dmat)
!     real(8), optional, target, intent(inout) :: xc_fmat(xc_mat_dim, xc_mat_dim, xc_nz, xc_nr_fmat)
      real(8), optional, target                :: xc_fmat(xc_mat_dim, xc_mat_dim, xc_nz, xc_nr_fmat)
      integer, optional, target, intent(in)    :: xc_dmat_ih_sym(xc_nr_dmat)
      integer, optional, target, intent(in)    :: xc_dmat_pg_sym(xc_nr_dmat)
      integer, optional, target, intent(in)    :: xc_fmat_pg_sym(xc_nr_fmat)

      integer, optional,         intent(in)    :: xc_nr_atoms
      real(8), optional, target, intent(inout) :: xc_property_gradient(*)

      logical, optional,         intent(in)    :: xc_do_potential
      logical, optional,         intent(in)    :: xc_do_london_rhs_direct
      logical, optional,         intent(in)    :: xc_do_london_rhs_ro
      logical, optional,         intent(in)    :: xc_do_geo_0
      logical, optional,         intent(in)    :: xc_do_geo_1
      logical, optional,         intent(in)    :: xc_do_geo_2
      logical, optional,         intent(in)    :: xc_do_london_ks
      logical, optional,         intent(in)    :: xc_do_london_lr
      logical, optional,         intent(in)    :: xc_do_overlap_diagnostic
      logical, optional,         intent(in)    :: xc_do_fde_update_vemb
      logical, optional,         intent(in)    :: xc_do_fde_response

      integer, optional,         intent(in)    :: xc_response_order_mo
      integer, optional,         intent(in)    :: xc_response_order_ao
!     --------------------------------------------------------------------------
    real(8), external :: ddot
   integer, parameter :: max_mat_number = 10
   integer, target    :: default_pg_sym(max_mat_number)
   integer, target    :: default_ih_sym(max_mat_number)
    real(8)             :: save_ddot1, save_ddot2



!   ----------------------------------------------------------------------------

!     start timer
      time_integration_start = second()
!     reset timer
      time_ao          = 0.0d0
      time_matrix_dist = 0.0d0
      time_density     = 0.0d0
      time_derv        = 0.0d0

      nr_electrons = 0.0d0
      xc_energy    = 0.0d0
      tau_integral = 0.0d0

      parallel_xc = (xc_mpi_get_nr_proc() > 1)

#ifdef VAR_MPI
      if (parallel_xc) then
         call sync_dft_cfg()
         call xc_mpi_bcast(fun_is_automatic)
      end if
#endif

      if (xc_mpi_is_master()) then
         mat_dim = xc_mat_dim
         nz      = xc_nz
         nr_dmat = xc_nr_dmat
         nr_fmat = xc_nr_fmat
      end if

!     set defaults
      do_potential          = .false.
      response_order_mo     = 0
      response_order_ao     = 0
      do_london_rhs_direct  = .false.
      do_london_rhs_ro      = .false.
      do_geo_0              = .false.
      do_geo_1              = .false.
      do_geo_2              = .false.
      do_london_ks          = .false.
      do_london_lr          = .false.
      do_overlap_diagnostic = .false.
      do_fde_update_vemb    = .false.
      do_fde_response       = .false.
      nr_atoms              = 0

      if (xc_mpi_is_master()) then
!        change defaults
         if (present(xc_do_potential))          do_potential          = xc_do_potential
         if (present(xc_response_order_mo))     response_order_mo     = xc_response_order_mo
         if (present(xc_response_order_ao))     response_order_ao     = xc_response_order_ao
         if (present(xc_do_london_rhs_direct))  do_london_rhs_direct  = xc_do_london_rhs_direct
         if (present(xc_do_london_rhs_ro))      do_london_rhs_ro      = xc_do_london_rhs_ro
         if (present(xc_do_geo_0))              do_geo_0              = xc_do_geo_0
         if (present(xc_do_geo_1))              do_geo_1              = xc_do_geo_1
         if (present(xc_do_geo_2))              do_geo_2              = xc_do_geo_2
         if (present(xc_do_london_ks))          do_london_ks          = xc_do_london_ks
         if (present(xc_do_london_lr))          do_london_lr          = xc_do_london_lr
         if (present(xc_nr_atoms))              nr_atoms              = xc_nr_atoms
         if (present(xc_do_overlap_diagnostic)) do_overlap_diagnostic = xc_do_overlap_diagnostic
         if (present(xc_do_fde_update_vemb))    do_fde_update_vemb    = xc_do_fde_update_vemb
         if (present(xc_do_fde_response))       do_fde_response       = xc_do_fde_response
      end if

#ifdef VAR_MPI
      if (parallel_xc) then
!        broadcast non-allocatables
         call xc_mpi_bcast(mat_dim)
         call xc_mpi_bcast(nz)
         call xc_mpi_bcast(nr_dmat)
         call xc_mpi_bcast(nr_fmat)
         call xc_mpi_bcast(nr_atoms)
         call xc_mpi_bcast(do_potential)
         call xc_mpi_bcast(do_london_rhs_direct)
         call xc_mpi_bcast(do_london_rhs_ro)
         call xc_mpi_bcast(response_order_mo)
         call xc_mpi_bcast(response_order_ao)
         call xc_mpi_bcast(do_geo_0)
         call xc_mpi_bcast(do_geo_1)
         call xc_mpi_bcast(do_geo_2)
         call xc_mpi_bcast(do_london_ks)
         call xc_mpi_bcast(do_london_lr)
         call xc_mpi_bcast(do_overlap_diagnostic)
         call xc_mpi_bcast(do_fde_update_vemb)
         call xc_mpi_bcast(do_fde_response)
      end if
#endif

!     nullify(tmat_1)
!     nullify(tmat_2)
      if (do_geo_1) then
         allocate(tmat_1(mat_dim, mat_dim))
         tmat_1(:, :) = xc_dmat(:, :, 1, 1)
      end if

      if (do_geo_2) then
         allocate(tmat_1(mat_dim, mat_dim))
         allocate(tmat_2(mat_dim, mat_dim))
         tmat_1(:, :) = xc_dmat(:, :, 1, 1)
         tmat_2(:, :) = xc_dmat(:, :, 1, 2)
      end if

      nullify(dmat_0)
      if (.not. xc_mpi_is_master()) then
         allocate(dmat_0_container(mat_dim, mat_dim))
         dmat_0 => dmat_0_container
      else
         dmat_0 => xc_dmat_0
      end if

      nullify(dmat)
      if (nr_dmat > 0) then
         if (.not. xc_mpi_is_master()) then
            allocate(dmat_container(mat_dim, mat_dim, nz, nr_dmat))
            dmat => dmat_container
         else
            if (present(xc_dmat)) then
               dmat => xc_dmat
            end if
         end if
      end if

      nullify(fmat)
      if (nr_fmat > 0) then
         if (.not. xc_mpi_is_master()) then
            allocate(fmat_container(mat_dim, mat_dim, nz, nr_fmat))
            fmat_container = 0.0d0
            fmat => fmat_container
         else
            if (present(xc_fmat)) then
               fmat => xc_fmat
            end if
         end if
      end if

      default_pg_sym = 1
      default_ih_sym = 0

      nullify(fmat_pg_sym)
      nullify(dmat_pg_sym)
      nullify(dmat_ih_sym)

      if (.not. xc_mpi_is_master()) then
         if (nr_fmat > 0) then
            allocate(fmat_pg_sym_container(nr_fmat))
            fmat_pg_sym => fmat_pg_sym_container
         end if
         if (nr_dmat > 0) then
            allocate(dmat_pg_sym_container(nr_dmat))
            allocate(dmat_ih_sym_container(nr_dmat))
            dmat_pg_sym => dmat_pg_sym_container
            dmat_ih_sym => dmat_ih_sym_container
         end if
      else
         if (present(xc_fmat_pg_sym)) then
            fmat_pg_sym => xc_fmat_pg_sym
         else
            fmat_pg_sym => default_pg_sym(1:nr_fmat)
         end if
         if (present(xc_dmat_pg_sym)) then
            dmat_pg_sym => xc_dmat_pg_sym
         else
            dmat_pg_sym => default_pg_sym(1:nr_dmat)
         end if
         if (present(xc_dmat_ih_sym)) then
            dmat_ih_sym => xc_dmat_ih_sym
         else
            dmat_ih_sym => default_ih_sym(1:nr_dmat)
         end if
      end if

      do_geo = .false.
      if (do_geo_0 .or. do_geo_1 .or. do_geo_2) then
         do_geo = .true.
      end if

      nullify(property_gradient)
      if (.not. xc_mpi_is_master()) then
         if (do_geo) then
            allocate(property_gradient_container(nr_atoms*3))
            property_gradient => property_gradient_container
            property_gradient = 0.0d0
         end if
      else
         if (present(xc_property_gradient)) then
            property_gradient => xc_property_gradient(1:nr_atoms*3)
            property_gradient = 0.0d0
         end if
      end if

!     set and sync functional
#ifdef MOD_XCFUN
      if (fun_is_automatic) then
         call xcfun_set_functional(xc_fun,       response_order_mo, parallel_xc)
         call xcfun_set_functional(xc_fun_alda,  response_order_mo, parallel_xc)
         call xcfun_set_functional(xc_fun_xalda, response_order_mo, parallel_xc)
      else
#endif
#ifdef VAR_MPI
      if (parallel_xc) then
         call dftsyncfunc(xc_mpi_is_master())
      end if
#endif
#ifdef MOD_XCFUN
      end if
#endif

!     check whether gga code structure will be used
      use_gga_qr = (fun_is_gga() .or. fun_is_tau_mgga())
      use_gga_qi = (fun_is_gga() .or. fun_is_tau_mgga())
      if (response_order_mo > 0) then
         if (dft_cfg_alda_hs) use_gga_qr = .false.
         if (dft_cfg_alda_ha) use_gga_qi = .false.
      end if
      if (use_gga_qr .or. use_gga_qi) then
         use_gga = .true.
      end if
!     in qr it is not so easy to separate h+ and h-
!     contributions since they may mix
!     therefore i quit here if the user wants qr
!     and not the full alda or xalda
!     it is really unfair to quit as late as here
!     this can be checked earlier
!     but this is anyway an exotic option
      if (response_order_mo > 1) then
         if (dft_cfg_alda_hs .and. (.not. dft_cfg_alda_ha)) then
            call quit('qr and partial alda is not implemented')
         end if
         if (dft_cfg_alda_ha .and. (.not. dft_cfg_alda_hs)) then
            call quit('qr and partial alda is not implemented')
         end if
      end if

      max_ao_g_order = 0
      if (use_gga) then
        max_ao_g_order = 1
      end if
      if (do_geo) then
        max_ao_g_order = max_ao_g_order + 1
      end if
     
      max_ao_m_order = 0
      if (do_london_rhs_direct) then
        max_ao_m_order = 1
      end if
     
      call ao_eval_init(max_ao_g_order,    &
                        max_ao_m_order,    &
                        dft_cfg_screening, &
                        .true.)

      if (xc_mpi_is_master()) then
!        initial expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot1 = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
         end if

!        get average density from open shells
         if (nr_open_shells > 0) call gtdoav(mat_dim, nz, 1, dmat_0, 1)

         call scale_density_matrices(2.0d0)

         call insert_half_phases()
      end if

#ifdef VAR_MPI
      if (parallel_xc) then
!        broadcast allocatables
         call xc_mpi_bcast(dmat_0)
         if (nr_fmat > 0) then
!           never bcast fmat
!           only master's fmat contains "HF" contribution
!           otherwise it is multiplied
            call xc_mpi_bcast(fmat_pg_sym)
         end if
         if (nr_dmat > 0) then
            call xc_mpi_bcast(dmat)
            call xc_mpi_bcast(dmat_pg_sym)
            call xc_mpi_bcast(dmat_ih_sym)
         end if
      end if
#endif

!     allocate saop density matrices
      if (dft_cfg_saop) then
         dft_cfg_saop_is_active = .true.
         allocate(dmat_saop_outer(mat_dim*mat_dim))
         if (dft_cfg_saop_with_response_part) then
            allocate(dmat_saop_resp(mat_dim*mat_dim))
         end if
!        only master calculates saop density matrices
         if (xc_mpi_is_master()) then
            call dftac_saop_get_dmats(dmat_saop_outer, &
                                      dmat_saop_resp)
         end if
#ifdef VAR_MPI
         if (parallel_xc) then
            call xc_mpi_bcast(dmat_saop_outer)
            if (dft_cfg_saop_with_response_part) then
               call xc_mpi_bcast(dmat_saop_resp)
            end if
         end if
#endif
      end if

#ifdef MOD_UNRELEASED
!     aspg, toying with having the frozen density and adfgrid available
!     for the respose calcs or otherwise...
      if (have_frozen_density) then
         call FDE_testFrozenDensity()
      endif
#endif

!#define PRINT_FMAT
#ifdef PRINT_FMAT
      write(*, *) 'debug: XC potential matrix zeroed out prior to integration'
      call dzero(fmat, mat_dim*mat_dim*nz)
#endif

#ifdef MOD_UNRELEASED
      if (do_overlap_diagnostic) then
         call overlap_diagnostic_init(parallel_xc, xc_mpi_is_master())
      end if
#endif

      call loop_over_batches()

#ifdef MOD_UNRELEASED
      if (do_overlap_diagnostic) then
#ifdef VAR_MPI
         if (parallel_xc) then
            call overlap_diagnostic_collect(xc_mpi_is_master())
         end if
#endif
         if (xc_mpi_is_master()) then
            call o_matrix_to_disc()
         end if
      end if
#endif

#ifdef VAR_MPI
      if (parallel_xc) then
!        collect results from processors
         call xc_mpi_reduce(nr_electrons)
         call xc_mpi_reduce(nr_points_used)
         if (do_potential) then
            call xc_mpi_reduce(xc_energy)
         end if
         if (nr_fmat > 0) then
            call xc_mpi_reduce(fmat)
         end if
         if (do_geo) then
            call xc_mpi_reduce(property_gradient)
         end if
      end if
#endif

!     deallocate saop density matrices
      if (dft_cfg_saop) then
        deallocate(dmat_saop_outer)
        if (dft_cfg_saop_with_response_part) then
          deallocate(dmat_saop_resp)
        end if
      end if

      if (xc_mpi_is_master()) then
         call scale_density_matrices(0.5d0)

!        get average density from open shells
         if (nr_open_shells > 0) call gtdoav(mat_dim, nz, 1, dmat_0, -1)

!        symmetrize matrices
         if (dft_cfg_blocked .or. (use_gga                          &
                                   .and. .not. do_london_rhs_direct &
                                   !.and. .not. do_london_rhs_ro &
                                   .and. .not. do_london_ks         &
!radovan: with do_london_lr matrices have to be symmetrized!
!                                  .and. .not. do_london_lr         &
                                   .and. .not. do_geo)) then
            if (nr_dmat > 0) then
               call gga_sym(mat_dim, nz, nr_fmat, fmat, fmat_pg_sym, dmat_ih_sym)
            else
               call gga_sym(mat_dim, nz, 1, fmat, (/1/), (/1/))
            end if
         end if

! gosia: later i will join it with what is done above (now it is not working if i just
! add another if)
         if (do_london_rhs_direct) then
           call mat_lao_sym(mat_dim, nz, nr_fmat, fmat, fmat_pg_sym)
         end if

         call insert_half_phases()

!        final expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot2    = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
            xc_mat_energy = save_ddot2 - save_ddot1
         end if

         call report()
      end if

      call nullify_and_release()

#ifdef PRINT_FMAT
      write(*, *) 'debug: real part of the XC potential matrix'
      call prqmat(xc_fmat,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  1,              &
                  (/1, 2, 3, 4/), &
                  file_unit)
#endif

   end subroutine

   subroutine nullify_and_release()

      if (associated(dmat)) then
         nullify(dmat)
      end if
      if (associated(fmat)) then
         nullify(fmat)
      end if
      if (associated(dmat_0)) then
         nullify(dmat_0)
      end if
      if (associated(dmat_pg_sym)) then
         nullify(dmat_pg_sym)
      end if
      if (associated(dmat_ih_sym)) then
         nullify(dmat_ih_sym)
      end if
      if (associated(fmat_pg_sym)) then
         nullify(fmat_pg_sym)
      end if
      if (associated(property_gradient)) then
         nullify(property_gradient)
      end if
      if (do_geo_1) then
         deallocate(tmat_1)
      end if
      if (do_geo_2) then
         deallocate(tmat_1)
         deallocate(tmat_2)
      end if
!     if (associated(tmat_1)) then
!        nullify(tmat_1)
!     end if
!     if (associated(tmat_2)) then
!        nullify(tmat_2)
!     end if

      if (.not. xc_mpi_is_master()) then
         if (allocated(dmat_container)) then
            deallocate(dmat_container)
         end if
         if (allocated(fmat_container)) then
            deallocate(fmat_container)
         end if
         if (allocated(dmat_0_container)) then
            deallocate(dmat_0_container)
         end if
         if (allocated(dmat_pg_sym_container)) then
            deallocate(dmat_pg_sym_container)
         end if
         if (allocated(dmat_ih_sym_container)) then
            deallocate(dmat_ih_sym_container)
         end if
         if (allocated(fmat_pg_sym_container)) then
            deallocate(fmat_pg_sym_container)
         end if
         if (allocated(property_gradient_container)) then
            deallocate(property_gradient_container)
         end if
      end if

   end subroutine

   subroutine scale_density_matrices(f)

!     --------------------------------------------------------------------------
      real(8) :: f
!     --------------------------------------------------------------------------

      call dscal(mat_dim*mat_dim, f, dmat_0, 1)
      if (nr_dmat > 0) then
        call dscal(mat_dim*mat_dim*nz*nr_dmat, f, dmat, 1)
      end if

   end subroutine

   subroutine insert_half_phases()

!     --------------------------------------------------------------------------
      integer :: iz, im, iq
!     --------------------------------------------------------------------------

      if (nz < 4) then
         if (nr_dmat > 0) then
            do iz = 1, nz
               do im = 1, nr_dmat
                  iq = pq_to_uq(iz, dmat_pg_sym(im) - 1)
                  call q2bphase('D', iq, 1, dmat(1, 1, iz, im))
               end do
            end do
         end if
         if (nr_fmat > 0) then
            do iz = 1, nz
               do im = 1, nr_fmat
                  iq = pq_to_uq(iz, fmat_pg_sym(im) - 1)
                  call q2bphase('F', iq, 1, fmat(1, 1, iz, im))
               end do
            end do
         end if
      end if

   end subroutine





   subroutine integrate()

!     --------------------------------------------------------------------------
!     unperturbed density
      real(8) :: n_0(max_block_length)
!     unperturbed density gradient
      real(8) :: gn_0_pointwise(3)

      real(8) :: tau
      real(8) :: gnn_0(max_block_length)
      real(8) :: n_0_t(max_block_length)
      real(8) :: gnn_0_t(max_block_length)
      real(8) :: w(max_block_length)
!     --------------------------------------------------------------------------
!     geometric business
      real(8) :: r_t1(max_block_length)
      real(8) :: r_t2(max_block_length)
!     --------------------------------------------------------------------------
!     fde and fde response
      logical :: fde_resp_verbose =.false.
      real(8) :: z0tt(3)
!     --------------------------------------------------------------------------
      integer :: id, m
      integer :: ii, ji
      real(8) :: temp
      real(8) :: block_threshold
      real(8) :: rho13, rho_outer, rho_resp, f
    integer :: ipoint
    integer :: k
      real(8) :: VT(3)
    integer              :: derv_length, im, i
    real(8), allocatable :: derv(:, :)
    real(8), allocatable :: derv_nadd(:, :)

    logical :: print_data_at_point

    real(8) :: n_b(3)
    real(8) :: gnn_b(3)
    real(8) :: gn_b(3, 3)

      logical :: un_above_threshold
      logical :: us_above_threshold
      logical :: ugn_above_threshold
      logical :: ugs_above_threshold
      logical :: utn_above_threshold
      logical :: uts_above_threshold
      logical :: us_lao_above_threshold
      logical :: ugs_lao_above_threshold

      real(8), allocatable :: ao_block(:, :)
      real(8), allocatable :: buffer_block(:, :)

    type(omega_prefactor) :: u

    allocate(u%n (      nr_fmat))
    allocate(u%s (   3, nr_fmat))
    allocate(u%gn(   3, nr_fmat))
    allocate(u%gs(3, 3, nr_fmat))
    allocate(u%tn(      nr_fmat))
    allocate(u%ts(   3, nr_fmat))
    allocate(u%s_lao(3, nr_fmat))
    allocate(u%gs_lao(3, 3, nr_fmat))





    allocate(ao    (total_ao_vector_length))
    ao = 0.0d0
    allocate(buffer(total_ao_vector_length))
    buffer = 0.0d0

    if (dft_cfg_blocked) then
       block_threshold = get_block_threshold(dft_cfg_dfthr0, mat_dim, dmat_0)
    end if


!   check what functional derivatives do we need
!   ============================================

    max_fun_derv = 0

    if (do_potential) then
       max_fun_derv = 1
    end if
    if (response_order_mo > 0) then
       max_fun_derv = response_order_mo + 1
    end if
    if (do_london_rhs_direct) then
       if (dft_cfg_no_sdft) then
          max_fun_derv = 1
       else
          max_fun_derv = 2
       end if
    end if
    if (do_geo_0) then
       max_fun_derv = 1
    end if
    if (do_geo_1) then
       max_fun_derv = 2
    end if
    if (do_geo_2) then
       max_fun_derv = 3
    end if
    if (do_london_ks) then
       max_fun_derv = 1
    end if
    if (do_london_lr) then
       max_fun_derv = 2
    end if

#ifdef MOD_OPENRSP
    if (do_london_ks .or. do_london_lr) then
       call xc_london_c1_init(mat_dim)
    end if
#endif

#ifdef MOD_XCFUN
#ifdef MOD_OPENRSP
    if (response_order_ao > 0) then
      derv_length = xc_output_length(xc_fun%id, response_order_ao)
      allocate(derv(derv_length, max_block_length))
      allocate(derv_nadd(derv_length, max_block_length))
    else
#endif
#endif
      derv_length = nr_nonzero_derv
      allocate(derv(max_block_length, 0:derv_length))
      allocate(derv_nadd(max_block_length, 0:derv_length))
#ifdef MOD_XCFUN
#ifdef MOD_OPENRSP
    end if
#endif
#endif

!       better safe than sorry
        derv      = 0.0d0
        derv_nadd = 0.0d0

!       set it to zero, otherwise it is undefined with lda
        gnn_0 = 0.0d0


#ifdef MOD_UNRELEASED
   if (dft_cfg_blocked) then
      loop_over_points_block: do ipoint = 1, nr_points_on_this_proc, max_block_length
         block_length = min(max_block_length, nr_points_on_this_proc - ipoint + 1)

         allocate(ao_block(block_length,     total_ao_vector_length))
         ao_block = 0.0d0
         allocate(buffer_block(block_length, total_ao_vector_length))
         buffer_block = 0.0d0

         call get_ao(                                  &
                     block_length,                     &
                     rx(ipoint:ipoint+block_length-1), &
                     ry(ipoint:ipoint+block_length-1), &
                     rz(ipoint:ipoint+block_length-1), &
                     ao_block,                         &
                     buffer_block                      &
                    )

         call xc_contribution_blas3(                          &
                                    use_gga,                  &
                                    (response_order_mo == 1), &
                                    block_length,             &
                                    block_threshold,          &
                                    mat_dim,                  &
                                    nz,                       &
                                    buffer,                   &
                                    ao_block,                 &
                                    nr_fmat,                  &
                                    fmat,                     &
                                    dmat_0,                   &
                                    nr_dmat,                  &
                                    dmat,                     &
                                    dmat_pg_sym,              &
                                    dmat_ih_sym,              &
                                    derv_length,              &
                                    derv,                     &
                                    max_fun_derv,             &
                                    rw(ipoint),               &
                                    nr_electrons,             &
                                    xc_energy                 &
                                   )

         deallocate(ao_block)
         deallocate(buffer_block)

      end do loop_over_points_block
   else
#endif /* ifdef MOD_UNRELEASED */

   print_data_at_point = .true.

!  here starts the expensive loop over points
!  do not evaluate anything inside this loop
!  that does not change from point to point
   do ipoint = 1, nr_points_on_this_proc

      time_ao_start = second()
      call get_ao(rx(ipoint), &
                  ry(ipoint), &
                  rz(ipoint), &
                  ao,         &
                  buffer)



      time_ao = time_ao + second() - time_ao_start

      time_density_start = second()
      if (need_ao_order(1, 0)) then
!        density and density gradient
         call get_gn(n_0(1), gn_0_pointwise, 0, mat_dim, dmat_0, buffer, ao, ao(ao_start(1, 0):ao_end(1, 0)))
         gnn_0(1) = gn_0_pointwise(1)*gn_0_pointwise(1) + gn_0_pointwise(2)*gn_0_pointwise(2) + gn_0_pointwise(3)*gn_0_pointwise(3)
      else
!        density
         call get_n(n_0(1), 0, mat_dim, dmat_0, buffer, ao)
!        for mcscf/srdft also evaluate spin density here
!        call get_s(spin_density, 0, mat_dim, dmat_0, buffer, ao)
      end if
      time_density = time_density + second() - time_density_start

#ifdef MOD_UNRELEASED
      if (do_overlap_diagnostic) then
         call integrate_o_matrix(rw(ipoint), ao)
      end if
#endif

      if (n_0(1) > dft_cfg_dfthr0) then
         nr_points_used = nr_points_used + 1

         u%n  = 0.0d0
         u%s  = 0.0d0
         u%gn = 0.0d0
         u%gs = 0.0d0
         u%tn = 0.0d0
         u%ts = 0.0d0
         u%s_lao  = 0.0d0
         u%gs_lao  = 0.0d0

         w(1) = rw(ipoint)

#ifdef MOD_UNRELEASED
         if (fun_is_tau_mgga()) then
            call get_tau(tau,     &
                         0,       &
                         mat_dim, &
                         dmat_0,  &
                         buffer,  &
                         ao(ao_start(1, 0):ao_end(1, 0)))
            tau_integral = tau_integral + w(1)*tau
         end if

        if (do_fde_update_vemb .or. do_fde_response) then
           if ( fde_froz_n(ipoint) > (dft_cfg_dfthr0/10000) ) then
              n_0_t(1) = n_0(1) + fde_froz_n(ipoint)

              z0tt(1)  = gn_0_pointwise(1) + fde_froz_ng(ipoint,1)
              z0tt(2)  = gn_0_pointwise(2) + fde_froz_ng(ipoint,2)
              z0tt(3)  = gn_0_pointwise(3) + fde_froz_ng(ipoint,3)

              gnn_0_t(1) = z0tt(1)*z0tt(1) + z0tt(2)*z0tt(2) + z0tt(3)*z0tt(3)

           else
              n_0_t(1) = n_0(1)
              gnn_0_t(1) = gnn_0(1)
           endif
      endif
#endif /* ifdef MOD_UNRELEASED */

!       number of electrons
        nr_electrons = nr_electrons + w(1)*n_0(1)

        time_derv_start = second()
#ifdef MOD_XCFUN
#ifdef MOD_OPENRSP
!       if we run openrsp with xcfun we evaluate derivatives
!       somewhere else
        if (response_order_ao == 0) then
#endif
#endif
          call get_functional_derv(xc_fun,                    &
                                   xc_fun_alda,               &
                                   xc_fun_xalda,              &
                                   max_fun_derv,              &
                                   bllen,                     &
                                   w,                         &
                                   n_0,                       &
                                   gnn_0,                     &
                                   derv_length,               &
                                   derv,                      &
                                   alda_real=dft_cfg_alda_hs, &
                                   alda_imag=dft_cfg_alda_ha, &
                                   xalda=dft_cfg_xalda,       &
                                   tau=tau)

#ifdef MOD_XCFUN
           if (do_fde_update_vemb .or. do_fde_response) then
!
! aspg: this way we are forcing the same xc functional for the active and non-additive
!       contributions. this is not desired, but we have to overcome xcfun's limitation
!       on using only up to 5 functionals...
!
              call get_fde_nadd_fun_derv(max_fun_derv, &
                                         bllen,        &
                                         w,            &
                                         n_0,          &
                                         n_0_t,        &
                                         gnn_0,        &
                                         gnn_0_t,      &
                                         derv_length,  &
                                         derv_nadd,    &
                                         parallel_xc)

!      if (print_data_at_point) then
!        write (*,*) 'debug: at grid point ',ipoint
!        write (*,*) '       for embedding potential, response construction'
!        write (*,*) '       n_0:',n_0(1),'  n_0_t:',n_0_t(1),"   w:",w(1)
!        write (*,*) '       ep :',fde_elpot(ipoint)
!        ii = 1
!         do ji = 0, derv_length
!           write (*,*) 'ss(',ii,',',ji,'):',derv(ii,ji),'    nadd(',ii,',',ji,'):',derv_nadd(ii,ji)
!         enddo
!        print_data_at_point = .false.
!      endif

              do ii = 1, max_block_length
                do ji = 0, derv_length
                  if (ji >  1 ) &
      &             derv(ii,ji) = derv(ii,ji) + derv_nadd(ii,ji)
! we add the electrostatic potential from the frozen subsystem(s) to the derivative
                  if (ji == 1 .and. update_vemb) &
      &             derv(ii,ji) = derv(ii,ji) + derv_nadd(ii,ji) + w(1)*fde_elpot(ipoint)
                enddo
              enddo
           endif
#endif /*  MOD_XCFUN */


#ifdef MOD_XCFUN
#ifdef MOD_OPENRSP
        end if
#endif
#endif
        time_derv = time_derv + second() - time_derv_start

        time_matrix_dist_start = second()

!       exchange-correlation potential
!       ==============================

        if (do_potential) then

!         asymptotic correction
          if (dft_cfg_saop_is_active) then
            call get_n(rho_outer,       &
                       0,               &
                       mat_dim,         &
                       dmat_saop_outer, &
                       buffer,          &
                       ao)
            if (dft_cfg_saop_with_response_part) then
              call get_n(rho_resp,       &
                         0,              &
                         mat_dim,        &
                         dmat_saop_resp, &
                         buffer,         &
                         ao)
            end if
          end if
          if (dft_cfg_grac_is_active .or. dft_cfg_saop_is_active) then
            rho13 = n_0(1)**(1.0d0/3.0d0)
            call dftac(derv, w(1), n_0(1), rho13, &
     &                 rho_resp, rho_outer, gnn_0(1))
          end if

!         energy density
          xc_energy = xc_energy + derv(1, d0000000)

          u%n(1) = derv(1, d1000000)
          if (use_gga) then
!           this factor has nothing to do
!           with symmetrization later
!           it comes from the fact that derivative of z
!           yields two identical terms (one left, one right)
            u%gn(:, 1) = 2.0d0*derv(1, d0010000)*gn_0_pointwise
          end if
#ifdef MOD_UNRELEASED
         if (fun_is_tau_mgga()) then
            u%tn(1) = derv(1, d0000010)
         end if
#endif

        end if !exchange-correlation potential



!       direct-lao contribution to differentiated kohn-sham matrix
!       ==========================================================

!       gosia: this part is correct only for closed-shell systems
!       some xc contributions to differentiated Fock matrix have not been
!       programmed as for closed-shell systems they are = 0 by definition

        if (do_london_rhs_direct) then
            if (use_gga) then
              call gga_london_distribute_r(derv(1, d1000000),                 &
                                           2.0d0*derv(1, d0010000),           &
                                           gn_0_pointwise,                    &
                                           fmat_pg_sym,                       &
                                           mat_dim,                           &
                                           nz,                                &
                                           fmat,                              &
                                           ao(ao_start(0, 0):ao_end(0, 0)),   &
                                           ao(ao_start(0, 1):ao_end(0, 1)),   &
                                           ao(ao_start(1, 0):ao_end(1, 0)),   &
                                           ao(ao_start(1, 1):ao_end(1, 1)),   &
                                           rx(ipoint), ry(ipoint), rz(ipoint))
            else
              call lda_london_distribute_r(derv(1, d1000000),                 &
                                           fmat_pg_sym,                       &
                                           mat_dim,                           &
                                           nz,                                &
                                           fmat,                              &
                                           ao,                                &
                                           ao(ao_start(0, 1):ao_end(0, 1)),   &
                                           rx(ipoint), ry(ipoint), rz(ipoint))
            end if

          end if


!       linear response
!       ===============

        if (response_order_mo == 1) then
           if (do_london_rhs_direct .and. .not. dft_cfg_no_sdft) then

              !call lr_sdft_lao_direct(bllen,          &
              !                        mat_dim,        &
              !                        nz,             &
              !                        ao,             &
              !                        n_0,            &
              !                        gn_0_pointwise, &
              !                        gnn_0,          &
              !                        w,              &
              !                        use_gga_qi,     &
              !                        buffer,         &
              !                        derv_length,    &
              !                        derv, u,        &
              !                        rx(ipoint), ry(ipoint), rz(ipoint))

           else if (do_london_rhs_ro .and. .not. dft_cfg_no_sdft) then
              call lr_sdft_lao_reorth(bllen,          &
                                      mat_dim,        &
                                      nz,             &
                                      ao,             &
                                      n_0,            &
                                      gn_0_pointwise,           &
                                      gnn_0,            &
                                      w,              &
                                      use_gga_qi,     &
                                      dmat, &
                                      dmat_pg_sym, &
                                      buffer,         &
                                      derv_length,    &
                                      derv, u)

           else

              call lr(bllen,          &
                      mat_dim,        &
                      nz,             &
                      nr_dmat,        &
                      dmat,           &
                      dmat_pg_sym,    &
                      dmat_ih_sym,    &
                      .false.,        &
                      ao,             &
                      n_0,            &
                      gn_0_pointwise,           &
                      gnn_0,            &
                      w,              &
                      use_gga_qr,     &
                      use_gga_qi,     &
                      buffer,         &
                      derv_length,    &
                      derv, u)

          endif
        endif


!       quadratic response
!       ==================
        if (response_order_mo == 2) then

          call sdft_qr(bllen,           &
                       use_gga,                &
                       dft_cfg_no_sdft,        &
                       mat_dim, dmat,        &
                       dmat_pg_sym,dmat_ih_sym,ao,ao(ao_start(1, 0)), &
                       n_0(1),gnn_0(1),gn_0_pointwise,w(1), &
                       buffer,             &
                  derv_length,    &
                  derv, u)

!         the idea here is that   r^[B]               \omega_{ia}^[C]
!                                 r^[C]               \omega_{ia}^[B]
!                                 r^[BC]              \omega_{ia}
!         and the corresponding   q^[B]  \cdot \nabla \omega_{ia}^[C]
!                                 q^[C]  \cdot \nabla \omega_{ia}^[B]
!                                 q^[BC] \cdot \nabla \omega_{ia}

!         can be evaluated using lr code (below)
!         only the remaining part is done in sdft_qr (above)
!
          call lr(bllen,    &
                  mat_dim,         &
                  nz,              &
                  nr_dmat,       &
                  dmat,          &
                  dmat_pg_sym,     &
                  dmat_ih_sym,          &
                  .true.,          &
                  ao,            &
                  n_0,             &
                  gn_0_pointwise,             &
                  gnn_0,             &
                  w,               &
                  use_gga_qr,      &
                  use_gga_qi,      &
                  buffer,          &
                  derv_length,     &
                  derv, u)

        endif


#ifdef MOD_XCFUN
#ifdef MOD_OPENRSP
          if (response_order_ao > 0) then
            call xc_response_contribution(response_order_ao,                 &
                                          mat_dim,   &
                                          nz,        &
                                          fmat,      &
                                          nr_dmat, &
                                          dmat,    &
                                          dmat_pg_sym,      &
                                          dmat_ih_sym,      &
                                          fmat_pg_sym,      &
                                          use_gga,                           &
                                          ao(ao_start(0, 0):ao_end(0, 0)), &
                                          ao(ao_start(1, 0):ao_end(1, 0)), &
                                          n_0(1),                            &
                                          gn_0_pointwise,                               &
                                          rx(ipoint),                      &
                                          ry(ipoint),                      &
                                          rz(ipoint),                      &
                                          w(1),                              &
                                          buffer, &
                                          derv_length,  &
                                          derv)
          end if
#endif
#endif

#ifdef MOD_UNRELEASED
          if (do_geo_0) then
             if (use_gga) then
                call xc_geo_0_gga_c1(property_gradient,               &
                                     derv(1, d1000000),                    &
                                     derv(1, d0010000),                    &
                                     gn_0_pointwise,                             &
                                     mat_dim,                         &
                                     dmat_0,                          &
                                     ao,                              &
                                     ao(ao_start(1, 0):ao_end(1, 0)), &
                                     ao(ao_start(2, 0):ao_end(2, 0)))
             else
                call xc_geo_0_lda_c1(property_gradient, &
                                     derv(1, d1000000),      &
                                     mat_dim,           &
                                     dmat_0,            &
                                     ao,                &
                                     ao(ao_start(1, 0):ao_end(1, 0)))
             end if
          end if

          if (do_geo_1) then
             call get_n(r_t1(1), 0, mat_dim, tmat_1, buffer, ao)
             call xc_geo_1_lda_c1(property_gradient, &
                                  derv(1, d1000000),      &
                                  derv(1, d2000000),     &
                                  r_t1(1),           &
                                  mat_dim,           &
                                  dmat_0,            &
                                  tmat_1,            &
                                  ao,              &
                                  ao(ao_start(1, 0):ao_end(1, 0)))
          end if

          if (do_geo_2) then
             call get_n(r_t1(1), 0, mat_dim, tmat_1, buffer, ao)
             call get_n(r_t2(1), 0, mat_dim, tmat_2, buffer, ao)
             call xc_geo_2_lda_c1(property_gradient, &
                                  derv(1, d2000000),     &
                                  derv(1, d3000000),    &
                                  r_t1(1),           &
                                  r_t2(1),           &
                                  mat_dim,           &
                                  dmat_0,            &
                                  tmat_1,            &
                                  tmat_2,            &
                                  ao,              &
                                  ao(ao_start(1, 0):ao_end(1, 0)))
          end if

          if (do_london_ks) then
             if (use_gga) then
                call xc_london_c1_ks_gga(derv(1, d1000000),               &
                                         derv(1, d0010000),               &
                                         gn_0_pointwise,                             &
                                         mat_dim,                         &
                                         nz,                              &
                                         fmat,                            &
                                         ao,                              &
                                         ao(ao_start(1, 0):ao_end(1, 0)), &
                                         (/rx(ipoint), ry(ipoint), rz(ipoint)/))
             else
                call xc_london_c1_ks_lda(derv(1, d1000000), &
                                         mat_dim,           &
                                         nz,                &
                                         fmat,              &
                                         ao,                &
                                         (/rx(ipoint), ry(ipoint), rz(ipoint)/))
             end if
          end if

          if (do_london_lr) then
             if (use_gga) then
                call xc_london_c1_lr_gga_densities(mat_dim,                         &
                                                   nz,                              &
                                                   dmat,                            &
                                                   ao,                              &
                                                   ao(ao_start(1, 0):ao_end(1, 0)), &
                                                 (/rx(ipoint),                      &
                                                   ry(ipoint),                      &
                                                   rz(ipoint)/),                    &
                                                   gn_0_pointwise,                            &
                                                   n_b,                             &
                                                   gn_b,                            &
                                                   gnn_b)
             else
                call xc_london_c1_lr_lda_densities(mat_dim,      &
                                                   nz,           &
                                                   dmat,         &
                                                   ao,           &
                                                 (/rx(ipoint),   &
                                                   ry(ipoint),   &
                                                   rz(ipoint)/), &
                                                   n_b)
             end if
!radovan:    a short explanation of what is happening here
!            (because this is done differently than in Dalton):
!            the above routines evaluate magnetic field derivatives
!            (there are three: im = 1, 3)
!            of perturbed density variables using the transition density matrix dmat
!            next we form u%n and u%gn and distribute over matrices
!            as if it was a normal linear response run
!            this makes the otherwise very lengthy operation rather compact
             do im = 1, 3
                u%n(im) = derv(1, d2000000)*n_b(im)
                if (use_gga) then
                   u%n(im) = u%n(im) + derv(1, d1010000)*gnn_b(im)
                   u%gn(1:3, im) = 2.0d0*gn_0_pointwise(1:3)*n_b(im)*derv(1, d1010000)   &
                                 + 2.0d0*gn_0_pointwise(1:3)*gnn_b(im)*derv(1, d0020000) &
                                 + 2.0d0*gn_b(1:3, im)*derv(1, d0010000)
                end if
             end do
          end if
#endif /* ifdef MOD_UNRELEASED */



!         distribute over fock matrices
!         =============================

          do im = 1, nr_fmat

             un_above_threshold = .false.
             if (dabs(u%n(im)) > tiny(0.0d0)) then
                un_above_threshold = .true.
             end if

             ugn_above_threshold = .false.
             if (maxval((/dabs(u%gn(1, im)), &
                          dabs(u%gn(2, im)), &
                          dabs(u%gn(3, im))/)) > tiny(0.0d0)) then
                ugn_above_threshold = .true.
             end if

             if (ugn_above_threshold) then
!               F_pq +=   u%n    \chi_p* \chi_q
!                     + 2 u%gn_i \chi_p* \nabla_i \chi_q
                call nabla_omega_real(fmat,                &
                                      im,                  &
                                      fmat_pg_sym(im) - 1, &
                                      mat_dim,             &
                                      nz,                  &
                                      u%n(im),             &
                                      u%gn(1, im),         &
                                      ao,                  &
                                      ao(ao_start(1, 0)))
             else
                if (un_above_threshold) then
!                  F_pq += u%n \chi_p* \chi_q
                   call omega_real(fmat,                &
                                   im,                  &
                                   fmat_pg_sym(im) - 1, &
                                   mat_dim,             &
                                   nz,                  &
                                   u%n(im),             &
                                   ao)
                end if
             end if

             us_above_threshold = .false.
             if (maxval((/dabs(u%s(1, im)), &
                          dabs(u%s(2, im)), &
                          dabs(u%s(3, im))/)) > tiny(0.0d0)) then
                us_above_threshold = .true.
             end if

             ugs_above_threshold = .false.
             if (maxval((/dabs(u%gs(1, 1, im)), &
                          dabs(u%gs(1, 2, im)), &
                          dabs(u%gs(1, 3, im)), &
                          dabs(u%gs(2, 1, im)), &
                          dabs(u%gs(2, 2, im)), &
                          dabs(u%gs(2, 3, im)), &
                          dabs(u%gs(3, 1, im)), &
                          dabs(u%gs(3, 2, im)), &
                          dabs(u%gs(3, 3, im))/)) > tiny(0.0d0)) then
                ugs_above_threshold = .true.
             end if

             if (ugs_above_threshold) then
!               F_pq +=   u%s_j   \chi_p* \Sigma_j \chi_q
!                     + 2 u%gs_ij \chi_p* \Sigma_j \nabla_i \chi_q
                call nabla_omega_imag(fmat,                &
                                      im,                  &
                                      fmat_pg_sym(im) - 1, &
                                      mat_dim,             &
                                      nz,                  &
                                      u%s(1, im),          &
                                      u%gs(1, 1, im),      &
                                      ao,                  &
                                      ao(ao_start(1, 0)))
             else
                if (us_above_threshold) then
!                  F_pq += u%s_j \chi_p* \Sigma_j \chi_q
                   call omega_imag(fmat,                &
                                   im,                  &
                                   fmat_pg_sym(im) - 1, &
                                   mat_dim,             &
                                   nz,                  &
                                   u%s(1, im),          &
                                   ao)
                end if
             end if

#ifdef MOD_UNRELEASED
             utn_above_threshold = .false.
             if (dabs(u%tn(im)) > tiny(0.0d0)) then
                utn_above_threshold = .true.
             end if

             if (utn_above_threshold) then
!               F_pq += 0.5 u%tn (\nabla_i \chi_p*) (\nabla_i \chi_q)
                call omega_real(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                0.5d0*u%tn(im),      &
                                ao(ao_start(1, 0)))
                call omega_real(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                0.5d0*u%tn(im),      &
                                ao(ao_start(1, 0) + mat_dim))
                call omega_real(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                0.5d0*u%tn(im),      &
                                ao(ao_start(1, 0) + mat_dim*2))
             end if

             uts_above_threshold = .false.
             if (maxval((/dabs(u%ts(1, im)), &
                          dabs(u%ts(2, im)), &
                          dabs(u%ts(3, im))/)) > tiny(0.0d0)) then
                uts_above_threshold = .true.
             end if

             if (uts_above_threshold) then
!               F_pq += 0.5 u%ts_j (\nabla_i \chi_p*) \Sigma_j (\nabla_i \chi_q)
                u%ts(:, im) = 0.5d0*u%ts(:, im)
                call omega_imag(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                u%ts(1, im),         &
                                ao(ao_start(1, 0)))
                call omega_imag(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                u%ts(1, im),         &
                                ao(ao_start(1, 0) + mat_dim))
                call omega_imag(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                u%ts(1, im),         &
                                ao(ao_start(1, 0) + mat_dim*2))
             end if

!gosia:
!sdft+lao now separately, later will be moved to u%s and u%gs,
!now i left it separately for debugging

!radovan: this code does not seem to be covered by any test, is this true?

             us_lao_above_threshold = .false.
             if (maxval((/dabs(u%s_lao(1, im)), &
                          dabs(u%s_lao(2, im)), &
                          dabs(u%s_lao(3, im))/)) > tiny(0.0d0)) then
                us_lao_above_threshold = .true.
             end if

             ugs_lao_above_threshold = .false.
             if (maxval((/dabs(u%gs_lao(1, 1, im)), &
                          dabs(u%gs_lao(1, 2, im)), &
                          dabs(u%gs_lao(1, 3, im)), &
                          dabs(u%gs_lao(2, 1, im)), &
                          dabs(u%gs_lao(2, 2, im)), &
                          dabs(u%gs_lao(2, 3, im)), &
                          dabs(u%gs_lao(3, 1, im)), &
                          dabs(u%gs_lao(3, 2, im)), &
                          dabs(u%gs_lao(3, 3, im))/)) > tiny(0.0d0)) then
                ugs_lao_above_threshold = .true.
             end if

             if (ugs_lao_above_threshold) then
                call nabla_omega_imag(fmat,                &
                                      im,                  &
                                      fmat_pg_sym(im) - 1, &
                                      mat_dim,             &
                                      nz,                  &
                                      u%s_lao(1, im),      &
                                      u%gs_lao(1, 1, im),  &
                                      ao,                  &
                                      ao(ao_start(1, 0)))
             else
                if (us_lao_above_threshold) then
!                  F_pq += u%s_j \chi_p* \Sigma_j \chi_q
                   call omega_imag(fmat,                &
                                   im,                  &
                                   fmat_pg_sym(im) - 1, &
                                   mat_dim,             &
                                   nz,                  &
                                   u%s_lao(1:3, im),    &
                                   ao)
                end if
             end if
!end gosia

#endif /* ifdef MOD_UNRELEASED */

          end do

          time_matrix_dist = time_matrix_dist + second() - time_matrix_dist_start

      end if !if (n_0(1) > dft_cfg_dfthr0) then
    end do !do ipoint = 1, nr_points_on_this_proc
#ifdef MOD_UNRELEASED
   end if !if (dft_cfg_blocked) then
#endif /* ifdef MOD_UNRELEASED */

    deallocate(derv)
    deallocate(derv_nadd)

    deallocate(ao)
    deallocate(buffer)

    deallocate(u%n)
    deallocate(u%s)
    deallocate(u%gn)
    deallocate(u%gs)
    deallocate(u%tn)
    deallocate(u%ts)
    deallocate(u%s_lao)
    deallocate(u%gs_lao)

   end subroutine

   subroutine report()
    real(8) :: error_nr_electrons

      call occupation_info_read()

      error_nr_electrons  = nr_electrons &
                          - nr_electrons_total

      write(file_unit, *)

      write(file_unit, '(3x, a, i9)')                            &
            'number of grid points                          = ', &
            nr_points_total
      write(file_unit, '(3x, a, i9)')                            &
            'number of grid points below density threshold  = ', &
            nr_points_total - nr_points_used

      if (do_potential) then
         write(file_unit, '(3x, a, f26.16)')                        &
               'DFT exchange-correlation energy:               = ', &
               xc_energy
      end if

      write(file_unit, '(3x, a, f26.16)')                        &
            'number of electrons from numerical integration = ', &
            nr_electrons

      if (abs(error_nr_electrons) > dft_cfg_accepted_error) then
         write(file_unit, '(3x, a, i9)')                            &
               'number of electrons from orbital occupations   = ', &
               nr_electrons_total
         write(file_unit, '(3x, a, f26.16)')                        &
               'WARNING: error in the number of electrons      = ', &
               error_nr_electrons
         write(file_unit, '(12x, a, f15.6)')                        &
               'is larger than .ERRELS (set input)    =  ',         &
                dft_cfg_accepted_error
         write(file_unit, '(12x, a)')                              &
               'this can happen when starting from coefficients'
         write(file_unit, '(12x, a)')                              &
               'from a different geometry'
         write(file_unit, '(12x, a)')                               &
               'or it can mean that the quadrature '//              &
               'grid is inappropriate'
      end if
      write(file_unit, *)

!     timing
      time_integration = second() - time_integration_start
      call timtxt('  time spent in DFT integration                  =', &
                  time_integration, file_unit)
     !call timtxt('                AO evaluation                    =', &
     !            time_ao, file_unit)
     !call timtxt('                matrix distribution              =', &
     !            time_matrix_dist, file_unit)
     !call timtxt('                density evaluation               =', &
     !            time_density, file_unit)
     !call timtxt('                functional derivatives           =', &
     !            time_derv, file_unit)
      write(file_unit, '(3x, a, i9)')                           &
            'number of processors                           = ', &
            xc_mpi_get_nr_proc()

      write(file_unit, *)

   end subroutine





   function get_xc_energy()

!     --------------------------------------------------------------------------
      real(8) :: get_xc_energy
!     --------------------------------------------------------------------------

      get_xc_energy = xc_energy

   end function

   function get_xc_mat_energy()

!     --------------------------------------------------------------------------
      real(8) :: get_xc_mat_energy
!     --------------------------------------------------------------------------

      get_xc_mat_energy = xc_mat_energy

   end function

   subroutine loop_over_batches()
    real(8) :: dummy
    integer :: idummy
    integer             :: ibatch
    integer             :: i
    real(8), allocatable :: mpi_buffer(:)

      if (xc_mpi_is_master()) then
         if (do_fde_update_vemb .or. do_fde_response) then
!           don't open and rewind anything
!           data is already stored at this point
         else
            rewind file_unit_grid
         end if
      end if

      ibatch                = 0
      nr_points_total       = 0
      nr_points_used        = 0

!     loop over batches
      do

!        get number of points in this batch
         if (xc_mpi_is_master()) then
            if (num_grid_cfg_import_grid) then
               read(file_unit_grid, *) nr_points_batch
#ifdef MOD_UNRELEASED
            else if (do_fde_update_vemb .or. do_fde_response) then
               nr_points_batch = fde_grid_np
#endif
            else
               read(file_unit_grid) nr_points_batch
            end if
         end if

#ifdef VAR_MPI
         if (parallel_xc) then
            call xc_mpi_bcast(nr_points_batch)
         end if
#endif
         if (nr_points_batch < 0) exit

         allocate(rx(nr_points_batch))
         allocate(ry(nr_points_batch))
         allocate(rz(nr_points_batch))
         allocate(rw(nr_points_batch))

!radovan: grid stuff will be moved outside this module
         if (xc_mpi_is_master()) then
            if (num_grid_cfg_import_grid) then
               call num_grid_read(rx, ry, rz, rw, file_unit_grid, nr_points_batch)
#ifdef MOD_UNRELEASED
            else if (do_fde_update_vemb .or. do_fde_response) then
               if (parallel_xc) then
                  call quit('fix fde_response mpi')
               end if

               allocate(fde_froz_n(nr_points_batch))
               fde_froz_n = 0.0d0
               allocate(fde_elpot(nr_points_batch))
               fde_elpot = 0.0d0
               allocate(fde_froz_ng(nr_points_batch, 3))
               fde_froz_ng = 0.0d0

!              copy data over to rx, ry, rz, rw, fde_froz_n, fde_froz_ng
               do i = 1, nr_points_batch
                  rx(i) = fde_grid(1,i)
                  ry(i) = fde_grid(2,i)
                  rz(i) = fde_grid(3,i)
                  rw(i) = fde_grid(4,i)

                  fde_froz_n(i) = fde_frozen_prp(3,i)
                  fde_elpot(i)  = fde_frozen_prp(1,i) + fde_frozen_prp(2,i)

                  if (fde_nprop.eq.6) then
                     fde_froz_ng(i,1) = fde_frozen_prp(4,i)
                     fde_froz_ng(i,2) = fde_frozen_prp(5,i)
                     fde_froz_ng(i,3) = fde_frozen_prp(6,i)
                  else
                     fde_froz_ng(i,1) = 0.0d0
                     fde_froz_ng(i,2) = 0.0d0
                     fde_froz_ng(i,3) = 0.0d0
                  endif
               end do

#endif
            else
               call num_grid_read_unformatted(rx, ry, rz, rw, file_unit_grid, nr_points_batch)
            end if
         end if

         ibatch          = ibatch + 1
         nr_points_total = nr_points_total + nr_points_batch

         nr_points_on_this_proc = nr_points_batch
#ifdef VAR_MPI
         if (parallel_xc) then
!           redefines nr_points_on_this_proc
            call xc_mpi_distribute_points(rx, ry, rz, rw, nr_points_batch, nr_points_on_this_proc)
         end if
#endif

         call integrate()

         deallocate(rx)
         deallocate(ry)
         deallocate(rz)
         deallocate(rw)

#ifdef MOD_UNRELEASED
        if (do_fde_update_vemb .or. do_fde_response) then
           if (allocated(fde_froz_n))  deallocate(fde_froz_n)
           if (allocated(fde_elpot))   deallocate(fde_elpot)
           if (allocated(fde_froz_ng)) deallocate(fde_froz_ng)
!          don't try to read another batch
           goto 999
        end if
#endif

!     end loop over batches
      end do
 999  continue

   end subroutine




  subroutine lr(block_length,       &
                mat_dim,            &
                nz,                 &
                nr_mat,             &
                dmaterturbed,     &
                isym,               &
                ih,                 &
                allow_nonhermitian, &
                ao_p,               &
                n_0,                &
                gn_0_pointwise,               &
                gnn_0,                &
                w,                  &
                is_gga_qr,          &
                is_gga_qi,          &
                buffer,             &
                n,                  &
                derv, u)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(in)    :: dmaterturbed(*)
    integer, intent(in)    :: isym(nr_mat)
    integer, intent(in)    :: ih(nr_mat)
    logical, intent(in)    :: allow_nonhermitian
    real(8), intent(in)    :: ao_p(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gnn_0(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga_qr
    logical, intent(in)    :: is_gga_qi
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv(max_block_length, 0:n)
!   ----------------------------------------------------------------------------
    logical                :: calculate_qr, calculate_qi


    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)
    real(8)                ::  t_b, ts_b(3)

    real(8)                :: z_b
    real(8)                :: y_b(3)

    integer                :: k_start, k, im, im_off, j
    type(omega_prefactor) :: u
!   ----------------------------------------------------------------------------

    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

    do im = 1, nr_mat

      calculate_qr = .false.
      calculate_qi = .false.

      if (ih(im) == +1) calculate_qr = .true.
      if (ih(im) == -1) calculate_qi = .true.

      if (ih(im) == 0 .and. allow_nonhermitian) then
        calculate_qr = .true.
        calculate_qi = .true.
      end if

      if (dft_cfg_no_sdft) then
        calculate_qi = .false.
      end if

      im_off = mat_dim*mat_dim*nz*(im - 1)

!     hermitian contribution (qr = quaternion real)
      if (calculate_qr) then

#ifdef MOD_UNRELEASED
         if (fun_is_tau_mgga()) then
            call get_tau(t_b,     &
                         isym(im) - 1,             &
                         mat_dim, &
                         dmaterturbed(1 + im_off), &
                         buffer,  &
                         ao_p(ao_start(1, 0)))
         end if
#endif /* ifdef MOD_UNRELEASED */

         if (is_gga_qr) then
            if (allow_nonhermitian) then
               call get_gn_nonhermitian(r_b,                      &
                                        gr_b,                     &
                                        isym(im) - 1,             &
                                        mat_dim,                  &
                                        dmaterturbed(1 + im_off), &
                                        buffer,                   &
                                        ao_p,                     &
                                        ao_p(ao_start(1, 0)))
            else
               call get_gn(r_b,                      &
                           gr_b,                     &
                           isym(im) - 1,             &
                           mat_dim,                  &
                           dmaterturbed(1 + im_off), &
                           buffer,                   &
                           ao_p,                     &
                           ao_p(ao_start(1, 0)))
            end if
         else
            call get_n(r_b,                        &
                       isym(im) - 1,               &
                       mat_dim,                    &
                       dmaterturbed(1 + im_off), &
                       buffer,                     &
                       ao_p)
         end if

        u%n(im) = u%n(im) + derv(1, d2000000)*r_b

        if (is_gga_qr) then
           z_b = 2.0d0*(gn_0_pointwise(1)*gr_b(1) &
                      + gn_0_pointwise(2)*gr_b(2) &
                      + gn_0_pointwise(3)*gr_b(3))

           u%n(im) = u%n(im) + derv(1, d1010000)*z_b

           u%gn(1:3, im) = u%gn(1:3, im) + 2.0d0*gn_0_pointwise(1:3)*r_b*derv(1, d1010000) &
                                         + 2.0d0*gn_0_pointwise(1:3)*z_b*derv(1, d0020000) &
                                         + 2.0d0*gr_b(1:3)*derv(1, d0010000)
        end if

#ifdef MOD_UNRELEASED
        if (fun_is_tau_mgga()) then
           u%tn(im) = u%tn(im) + derv(1, d0000020)*t_b
           u%tn(im) = u%tn(im) + derv(1, d1000010)*r_b
           u%tn(im) = u%tn(im) + derv(1, d0010010)*z_b
           u%n(im)  = u%n(im)  + derv(1, d1000010)*t_b
           u%gn(1:3, im) = u%gn(1:3, im) + 2.0d0*gn_0_pointwise(1:3)*t_b*derv(1, d0010010)
        end if
#endif /* ifdef MOD_UNRELEASED */

      end if !end of qr contribution

!     antihermitian contribution (qi = quaternion imaginary)
      if (calculate_qi) then


        if (is_gga_qi) then
          if (allow_nonhermitian) then
            call get_gs_nonhermitian(s_b,                       &
                                     gs_b,                       &
                                     isym(im) - 1,               &
                                     mat_dim,                    &
                                     dmaterturbed(1 + im_off), &
                                     buffer,                     &
                                     ao_p,                       &
                                     ao_p(ao_start(1, 0)))
          else
            call get_gs(s_b,                       &
                        gs_b,                       &
                        isym(im) - 1,               &
                        mat_dim,                    &
                        dmaterturbed(1 + im_off), &
                        buffer,                     &
                        ao_p,                       &
                        ao_p(ao_start(1, 0)))
          end if
        else
        call get_s(s_b,                        &
                   isym(im) - 1,               &
                   mat_dim,                    &
                   dmaterturbed(1 + im_off), &
                   buffer,                     &
                   ao_p)
        end if
#ifdef MOD_UNRELEASED
        if (fun_is_tau_mgga()) then
           call get_tau_imag(ts_b,     &
                             isym(im) - 1,             &
                             mat_dim, &
                             dmaterturbed(1 + im_off), &
                             buffer,  &
                             ao_p(ao_start(1, 0)))
        end if
#endif /* ifdef MOD_UNRELEASED */

        do k = k_start, 3
           u%s(k, im) = u%s(k, im) + derv(1, d0200000)*s_b(k)

           if (is_gga_qi) then
              y_b(k)  = gn_0_pointwise(1)*gs_b(1, k) &
                      + gn_0_pointwise(2)*gs_b(2, k) &
                      + gn_0_pointwise(3)*gs_b(3, k)

              u%s(k, im) = u%s(k, im) + derv(1, d0101000)*y_b(k)

              u%gs(1:3, k, im) = u%gs(1:3, k, im) +       gn_0_pointwise(1:3)*s_b(k)*derv(1, d0101000) &
                                                  +       gn_0_pointwise(1:3)*y_b(k)*derv(1, d0002000) &
                                                  + 2.0d0*gs_b(1:3, k)*derv(1, d0000100)
           end if
#ifdef MOD_UNRELEASED
           if (fun_is_tau_mgga()) then
              u%ts(k, im) = u%ts(k, im) + derv(1, d0000002)*ts_b(k)
              u%ts(k, im) = u%ts(k, im) + derv(1, d0100001)*s_b(k)
              u%ts(k, im) = u%ts(k, im) + derv(1, d0001001)*y_b(k)
              u%s(k, im) = u%s(k, im) + derv(1, d0100001)*ts_b(k)
              u%gs(1:3, k, im) = u%gs(1:3, k, im) + gn_0_pointwise(1:3)*ts_b(k)*derv(1, d0001001)
           end if
#endif /* ifdef MOD_UNRELEASED */
        end do

      end if !end of qi contribution

    end do

  end subroutine


  subroutine lr_sdft_lao_direct(block_length,       &
                                mat_dim,            &
                                nz,                 &
                                ao_p,               &
                                n_0,                &
                                gn_0_pointwise,     &
                                gnn_0,              &
                                w,                  &
                                is_gga_qi,          &
                                buffer,             &
                                n,                  &
                                derv, u,            &
                                px, py, pz)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(in)    :: ao_p(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gnn_0(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga_qi
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv(max_block_length, 0:n)
    real(8), intent(in)    :: px, py, pz
!   ----------------------------------------------------------------------------
    real(8)                ::  s_b_lao(3, 3), gs_b_lao(3, 3, 3)

    real(8)                :: z_b
    real(8)                :: y_b(3, 3)

    integer                :: k_start, k, im, im_off, j
    type(omega_prefactor) :: u
    integer :: i
!   ----------------------------------------------------------------------------

    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation


    if (is_gga_qi) then
      call get_gs_lao(s_b_lao,                           &
                      gs_b_lao,                          &
                      fmat_pg_sym,                       &
                      mat_dim,                           &
                      dmat_0,                            &
                      buffer,                            &
                      ao_p(ao_start(0, 0):ao_end(0, 0)), &
                      ao_p(ao_start(1, 0):ao_end(1, 0)), &
                      ao_p(ao_start(0, 1):ao_end(0, 1)), &
                      ao_p(ao_start(1, 1):ao_end(1, 1)), &
                      px, py, pz)
    else
      call get_s_lao(s_b_lao,                            &
                     fmat_pg_sym,                        &
                     mat_dim,                            &
                     dmat_0,                             &
                     ao_p(ao_start(0, 0):ao_end(0, 0)),  &
                     ao_p(ao_start(0, 1):ao_end(0, 1)),  &
                     px, py, pz,                         &
                     buffer)

    end if

    do im = 1, nr_fmat
        do k = k_start, 3

           u%s_lao(k, im) = u%s_lao(k, im) + derv(1, d0200000)*s_b_lao(k, im)

           if (is_gga_qi) then
              y_b(k, im)  = gn_0_pointwise(1)*gs_b_lao(1, k, im) &
                          + gn_0_pointwise(2)*gs_b_lao(2, k, im) &
                          + gn_0_pointwise(3)*gs_b_lao(3, k, im)

              u%s_lao(k, im) = u%s_lao(k, im) + derv(1, d0101000)*y_b(k, im)

              u%gs_lao(1:3, k, im) = u%gs_lao(1:3, k, im)                               &
                                   +(gn_0_pointwise(1)*s_b_lao(1, k)                    &
                                   + gn_0_pointwise(2)*s_b_lao(2, k)                    &
                                   + gn_0_pointwise(3)*s_b_lao(3, k))*derv(1, d0101000) &
                                   + gn_0_pointwise(1:3)*y_b(k, im)*derv(1, d0002000)   &
                                   + 2.0d0*gs_b_lao(1:3, k, im)*derv(1, d0000100)
           end if
        end do

    end do

  end subroutine


  subroutine lr_sdft_lao_reorth(block_length,       &
                                mat_dim,            &
                                nz,                 &
                                ao_p,               &
                                n_0,                &
                                gn_0_pointwise,     &
                                gnn_0,              &
                                w,                  &
                                is_gga,             &
                                tb_dmat,            &
                                isym,               &
                                buffer,             &
                                n,                  &
                                derv, u)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(in)    :: ao_p(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gnn_0(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga
    real(8), intent(in)    :: tb_dmat(*)
    integer, intent(in)    :: isym(*)
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv(max_block_length, 0:n)
!   ----------------------------------------------------------------------------
    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)

    real(8)                :: z_b
    real(8)                :: y_b(3)

    integer                :: k_start, k, im, im_off, j
    type(omega_prefactor) :: u
    integer :: i
!   ----------------------------------------------------------------------------

    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

    do im = 1, nr_fmat

       im_off = mat_dim*mat_dim*nz*(im - 1)

!      only charge density part
!       if (is_gga) then
!          call get_gn(r_b,                      &
!                      gr_b,                     &
!                      isym(im) - 1,             &
!                      mat_dim,                  &
!                      tb_dmat(1 + im_off),      &
!                      buffer,                   &
!                      ao_p,                     &
!                      ao_p(ao_start(1, 0)))
!       else
!          call get_n(r_b,                       &
!                     isym(im) - 1,              &
!                     mat_dim,                   &
!                     tb_dmat(1 + im_off),       &
!                     buffer,                    &
!                     ao_p)
!       end if
!
!       u%n(im) = u%n(im) + derv(1, d2000000)*r_b
!
!       if (is_gga) then
!          z_b = 2.0d0*(gn_0_pointwise(1)*gr_b(1) &
!                     + gn_0_pointwise(2)*gr_b(2) &
!                     + gn_0_pointwise(3)*gr_b(3))
!
!          u%n(im) = u%n(im) + derv(1, d1010000)*z_b
!
!          u%gn(1:3, im) = u%gn(1:3, im) + 2.0d0*gn_0_pointwise(1:3)*r_b*derv(1, d1010000) &
!                                        + 2.0d0*gn_0_pointwise(1:3)*z_b*derv(1, d0020000) &
!                                        + 2.0d0*gr_b(1:3)*derv(1, d0010000)
!       end if

!      SDFT
       if (.not. dft_cfg_no_sdft) then
          if (is_gga) then
             call get_gs(s_b,                       &
                         gs_b,                      &
                         isym(im) - 1,              &
                         mat_dim,                   &
                         tb_dmat(1 + im_off),       &
                         buffer,                    &
                         ao_p,                      &
                         ao_p(ao_start(1, 0)))
          else
             call get_s(s_b,                        &
                        isym(im) - 1,               &
                        mat_dim,                    &
                        tb_dmat(1 + im_off),        &
                        buffer,                     &
                        ao_p)
          end if

          do k = k_start, 3

             u%s(k, im) = u%s(k, im) + derv(1, d0200000)*s_b(k)

             if (is_gga) then
                y_b(k)  = gn_0_pointwise(1)*gs_b(1, k) &
                        + gn_0_pointwise(2)*gs_b(2, k) &
                        + gn_0_pointwise(3)*gs_b(3, k)

                u%s(k, im) = u%s(k, im) + derv(1, d0101000)*y_b(k)

                u%gs(1:3, k, im) = u%gs(1:3, k, im) +       gn_0_pointwise(1:3)*s_b(k)*derv(1, d0101000) &
                                                    +       gn_0_pointwise(1:3)*y_b(k)*derv(1, d0002000) &
                                                    + 2.0d0*gs_b(1:3, k)*derv(1, d0000100)
             end if
          end do

      end if
    end do

  end subroutine


      subroutine sdft_qr(block_length, &
                         add_grad, &
                         no_sdft, &
                         mat_dim,  &
                         dtrmat, &
                         isymop, &
                         ihrmop, &
                         gao, &
                         gao1, &
                         n_0, &
                         gnn_0, &
                         gn_0_pointwise, &
                         w, &
                         buf, &
                         derv_len, &
                         derv, u)

    type(omega_prefactor) :: u
      integer, intent(in) :: block_length
      integer, intent(in)    :: derv_len
      integer, intent(in)    :: mat_dim
      real(8), intent(in)    :: derv(max_block_length, 0:derv_len)
      logical   add_grad, &
     &          no_sdft
     real(8) ::          dtrmat(mat_dim,mat_dim,nz,3)
     integer ::          isymop(3)
     integer ::          ihrmop(3)
     real(8) ::          gao(*)
     real(8) ::          gao1(*)
     real(8) ::          buf(*)
     real(8) ::          n_0(*)
     real(8) ::          gnn_0(*)
     real(8) ::          w(*)
!
     real(8) ::          gn_0_pointwise (  3)
     real(8) ::          gr_b (  3)
     real(8) ::          gr_c (  3)
!
     real(8) ::          s_b  (  3)
     real(8) ::          s_c  (  3)
!
     real(8) ::          gs_b (3,3)
     real(8) ::          gs_c (3,3)
!
     real(8) ::          y_b  (  3)
     real(8) ::          y_c  (  3)
     real(8) ::          y_bc (  3)
!
     real(8) ::          x_bc (  3)
!
     real(8) ::          v_bc_h_y_0(3)
     real(8) ::          v_bc_h_y_b(3)
     real(8) ::          v_bc_h_y_c(3)
     real(8) ::          v_bc_h_x_b(3)
     real(8) ::          v_bc_h_x_c(3)


      real(8) ::  v_bc_h_gnn_0
      real(8) ::  s_bc_h
      real(8) ::  v_bc_h_z_b
      real(8) ::  v_bc_h_z_c
      real(8) ::  z_b, z_c, z_bc






!
!               here we construct (see thesis for notation):
!
!               k_{pq,0}^(bc) = n_0^(bc) \omega_{pq,0} + q_{0j}^(bc) \nabla_j \omega_{pq,0}
!               k_{pq,z}^(bc) = r_z^(bc) \omega_{pq,z} + q_{zj}^(bc) \nabla_j \omega_{pq,z}
!               k_{pq,y}^(bc) = r_y^(bc) \omega_{pq,y} + q_{yj}^(bc) \nabla_j \omega_{pq,y}
!               k_{pq,x}^(bc) = r_x^(bc) \omega_{pq,x} + q_{xj}^(bc) \nabla_j \omega_{pq,x}
!
!               all the rest is done elsewhere
!
!               scalar h+ part: n_0^[bc]
!               s_bc_h
!               vector h+ part: q_{0j}^[bc] with j = x,y,z
     real(8) ::          v_bc_h (  3)
!               scalar h- part: r_k^[bc]    with k = x,y,z
     real(8) ::          s_bc_ah(  3)
!               vector h- part: q_{kj}^[bc] with k = x,y,z and j = x,y,z
     real(8) ::          v_bc_ah(3,3)
      real(8) :: f, r_c, r_b
      integer :: i
    integer                :: k_start, k, im, im_off, j


      k_start = 1
      if(dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

      if (add_grad) then
         call get_gn_nonhermitian(r_b, &
     &                            gr_b, &
     &                            isymop(1) - 1, &
     &                            mat_dim,            &
     &                            dtrmat(1, 1, 1, 1), &
     &                            buf, &
     &                            gao, &
     &                            gao1)
         call get_gn_nonhermitian(r_c, &
     &                            gr_c, &
     &                            isymop(2) - 1, &
     &                            mat_dim,            &
     &                            dtrmat(1, 1, 1, 2), &
     &                            buf, &
     &                            gao, &
     &                            gao1)
      else
         call get_n(r_b, &
     &              isymop(1) - 1, &
     &              mat_dim,            &
     &              dtrmat(1, 1, 1, 1), &
     &              buf, &
     &              gao)
         call get_n(r_c, &
     &              isymop(2) - 1, &
     &              mat_dim,         &
     &              dtrmat(1, 1, 1, 2), &
     &              buf, &
     &              gao)
      end if

      if (.not. no_sdft) then
         if (add_grad) then
            call get_gs_nonhermitian(s_b, &
     &                               gs_b, &
     &                               isymop(1) - 1, &
     &                               mat_dim,            &
     &                               dtrmat(1, 1, 1, 1), &
     &                               buf, &
     &                               gao, &
     &                               gao1)
            call get_gs_nonhermitian(s_c, &
     &                               gs_c, &
     &                               isymop(2) - 1, &
     &                               mat_dim,            &
     &                               dtrmat(1, 1, 1, 2), &
     &                               buf, &
     &                               gao, &
     &                               gao1)
         else
           call get_s(s_b, &
     &                isymop(1) - 1, &
     &                mat_dim,            &
     &                dtrmat(1, 1, 1, 1), &
     &                buf, &
     &                gao)
           call get_s(s_c, &
     &                isymop(2) - 1, &
     &                mat_dim,            &
     &                dtrmat(1, 1, 1, 2), &
     &                buf, &
     &                gao)
         end if
      endif

      if (add_grad) then
!       z = \nabla r \cdot \nabla r
        z_b  = gn_0_pointwise(1)*gr_b(1) &
             + gn_0_pointwise(2)*gr_b(2) &
             + gn_0_pointwise(3)*gr_b(3)
        z_c  = gn_0_pointwise(1)*gr_c(1) &
             + gn_0_pointwise(2)*gr_c(2) &
             + gn_0_pointwise(3)*gr_c(3)
        z_bc = gr_b(1)*gr_c(1) &
             + gr_b(2)*gr_c(2) &
             + gr_b(3)*gr_c(3)
        z_b  = 2.0d0*z_b
        z_c  = 2.0d0*z_c
        z_bc = 2.0d0*z_bc
      end if

      if(add_grad .and. .not. no_sdft) then
!       y = \nabla r \cdot \nabla s
        do k = k_start,3
        y_b(k)  = gn_0_pointwise(1)*gs_b(1, k) &
                + gn_0_pointwise(2)*gs_b(2, k) &
                + gn_0_pointwise(3)*gs_b(3, k)
        y_c(k)  = gn_0_pointwise(1)*gs_c(1, k) &
                + gn_0_pointwise(2)*gs_c(2, k) &
                + gn_0_pointwise(3)*gs_c(3, k)
        y_bc(k) = gr_b(1)*gs_c(1, k) &
                + gr_b(2)*gs_c(2, k) &
                + gr_b(3)*gs_c(3, k) &
                + gr_c(1)*gs_b(1, k) &
                + gr_c(2)*gs_b(2, k) &
                + gr_c(3)*gs_b(3, k)
        enddo

!       x = \nabla s \cdot \nabla s
        do k = k_start,3
          x_bc(k) = gs_b(1, k)*gs_c(1, k) &
                  + gs_b(2, k)*gs_c(2, k) &
                  + gs_b(3, k)*gs_c(3, k)
          x_bc(k) = 2.0d0*x_bc(k)
        enddo
      endif

!-------------------------------------------------------------------------------
!     collect s_bc_h

      s_bc_h = 0.0d0

      s_bc_h = s_bc_h + derv(1, d3000000)*(r_b*r_c)

      if(add_grad) then
        s_bc_h = s_bc_h + derv(1, d2010000)*(r_b*z_c + z_b*r_c) &
     &                  + derv(1, d1020000)*(z_b*z_c) &
     &                  + derv(1, d1010000)*z_bc
      endif

      if(.not. no_sdft) then
        do k = k_start,3
          s_bc_h = s_bc_h + derv(1, d1200000)*(s_b(k)*s_c(k))

          if(add_grad) then
       s_bc_h = s_bc_h + derv(1, d1101000)*(s_b(k)*y_c(k) + y_b(k)*s_c(k)) &
     &                 + derv(1, d1002000)*(y_b(k)*y_c(k)) &
     &                 + derv(1, d1000100)*x_bc(k)
          endif
        enddo
      endif

!-------------------------------------------------------------------------------
!     collect s_bc_ah

      s_bc_ah = 0.0d0

      if(.not. no_sdft) then
        do k = k_start,3
       s_bc_ah(k) = s_bc_ah(k) &
     &            + derv(1, d1200000)*(r_b*s_c(k) + s_b(k)*r_c)

          if(add_grad) then
       s_bc_ah(k) = s_bc_ah(k)  &
     &            + derv(1, d1101000)*(r_b*y_c(k) + y_b(k)*r_c) &
     &            + derv(1, d0210000)*(z_b*s_c(k) + s_b(k)*z_c) &
     &            + derv(1, d0111000)*(z_b*y_c(k) + y_b(k)*z_c) &
     &            + derv(1, d0101000)*y_bc(k)
          endif
        enddo
      endif

!-------------------------------------------------------------------------------
!     collect vector components

      v_bc_h  = 0.0d0
      v_bc_ah = 0.0d0

      if(add_grad) then

        v_bc_h_gnn_0 = 0.0d0
        v_bc_h_z_b = 0.0d0
        v_bc_h_z_c = 0.0d0

        v_bc_h_gnn_0 = v_bc_h_gnn_0 &
     &             + derv(1, d0030000)*(z_b*z_c) &
     &             + derv(1, d1020000)*(r_b*z_c + z_b*r_c) &
     &             + derv(1, d2010000)*(r_b*r_c)

        v_bc_h_gnn_0 = v_bc_h_gnn_0 &
     &             + derv(1, d0020000)*z_bc
        v_bc_h_z_b = v_bc_h_z_b &
     &             + derv(1, d1010000)*r_c &
     &             + derv(1, d0020000)*z_c
        v_bc_h_z_c = v_bc_h_z_c &
     &             + derv(1, d1010000)*r_b &
     &             + derv(1, d0020000)*z_b

        if(.not. no_sdft) then

          v_bc_h_y_0 = 0.0d0
          v_bc_h_y_b = 0.0d0
          v_bc_h_y_c = 0.0d0

          v_bc_h_x_b = 0.0d0
          v_bc_h_x_c = 0.0d0

          do k = k_start,3
            v_bc_h_gnn_0  = v_bc_h_gnn_0 &
     &                  + derv(1, d0012000)*(y_b(k)*y_c(k)) &
     &                  + derv(1, d0111000)*(s_b(k)*y_c(k) + y_b(k)*s_c(k)) &
     &                  + derv(1, d0210000)*(s_b(k)*s_c(k))

            v_bc_h_y_0(k) = v_bc_h_y_0(k) &
     &                    + derv(1, d0002000)*y_bc(k)
            v_bc_h_y_b(k) = v_bc_h_y_b(k) &
     &                    + derv(1, d0101000)*s_c(k) &
     &                    + derv(1, d0002000)*y_c(k)
            v_bc_h_y_c(k) = v_bc_h_y_c(k) &
     &                    + derv(1, d0101000)*s_b(k) &
     &                    + derv(1, d0002000)*y_b(k)

            v_bc_h_x_b(k) = v_bc_h_x_b(k) &
     &                    + derv(1, d1000100)*r_c &
     &                    + derv(1, d0010100)*z_c
            v_bc_h_x_c(k) = v_bc_h_x_c(k) &
     &                    + derv(1, d1000100)*r_b &
     &                    + derv(1, d0010100)*z_b

            v_bc_h_y_0(k) = v_bc_h_y_0(k) &
     &                    + derv(1, d1101000)*(r_b*s_c(k) + s_b(k)*r_c) &
     &                    + derv(1, d1002000)*(r_b*y_c(k) + y_b(k)*r_c) &
     &                    + derv(1, d0111000)*(z_b*s_c(k) + s_b(k)*z_c) &
     &                    + derv(1, d0012000)*(z_b*y_c(k) + y_b(k)*z_c)
          enddo
        endif

        do j = 1,3
          v_bc_h(j) = v_bc_h(j) + v_bc_h_gnn_0*gn_0_pointwise(j)*2.0d0 &
     &                          + v_bc_h_z_b*gr_b(j)*2.0d0 &
     &                          + v_bc_h_z_c*gr_c(j)*2.0d0
          if(.not. no_sdft) then
            do k = k_start,3
              v_bc_ah(j,k) = v_bc_ah(j,k) + v_bc_h_y_0(k)*gn_0_pointwise(j) &
     &                                    + v_bc_h_y_b(k)*gr_b(j) &
     &                                    + v_bc_h_y_c(k)*gr_c(j)

              v_bc_h(j) = v_bc_h(j) + v_bc_h_y_b(k)*gs_b(j,k) &
     &                              + v_bc_h_y_c(k)*gs_c(j,k)

              v_bc_ah(j,k) = v_bc_ah(j,k) + v_bc_h_x_b(k)*gs_b(j,k)*2.0d0 &
     &                                    + v_bc_h_x_c(k)*gs_c(j,k)*2.0d0
            enddo
          endif
        enddo
      endif

!-------------------------------------------------------------------------------

      s_bc_h =     2.0d0*s_bc_h
      call dscal(3,2.0d0,v_bc_h,1)

      if (add_grad) then
        u%gn(:, 3) = u%gn(:, 3) + v_bc_h
      end if
      u%n(3) = u%n(3) + s_bc_h

      if (.not. no_sdft) then

        s_bc_ah = 2.0d0*s_bc_ah
        v_bc_ah = 2.0d0*v_bc_ah

        if (add_grad) then
          u%gs(:, :, 3) = u%gs(:, :, 3) + v_bc_ah
        end if
        u%s(:, 3) = u%s(:, 3) + s_bc_ah
      end if

      end subroutine

end module
