module xcint_main

   use interface_ao
   use interface_mo
   use interface_file_io
   use dft_cfg
   use interface_functional_read
   use xc_derv
   use xc_ac
   use extra
   use file_units
   use interface_grid
   use ao_eval
   use density_eval
   use xc_blas3
   use xc_max_block_length
   use xc_geometric
#ifdef MOD_UNRELEASED
   use overlap_diagnostic
   use xc_london_c1
#endif
#ifdef MOD_OPENRSP
   use xc_response
#endif
   use xc_mpi
   use xcfun

   implicit none

   public integrate_xc
   public get_xc_energy
   public get_xc_mat_energy
   public get_nr_electrons_integrated

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

   real(8), allocatable :: dmat_saop_outer(:)
   real(8), allocatable :: dmat_saop_resp(:)
   real(8), allocatable :: rx(:)
   real(8), allocatable :: ry(:)
   real(8), allocatable :: rz(:)
   real(8), allocatable :: rw(:)
   real(8), allocatable :: dmat_sorted(:, :)

   real(8) :: tau_integral
   real(8) :: nr_electrons_integrated
   real(8) :: xc_energy
   real(8) :: xc_mat_energy

   logical :: use_gga_qr = .false. !quaternion real
   logical :: use_gga_qi = .false. !quaternion imaginary
   logical :: use_gga    = .false.

   logical :: use_blas3_scheme = .false.

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

contains

   subroutine integrate_xc(xc_mat_dim,               &
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

      integer, optional,         intent(in)    :: xc_response_order_mo
      integer, optional,         intent(in)    :: xc_response_order_ao
!     --------------------------------------------------------------------------
      real(8), external :: ddot
      integer, parameter :: max_mat_number = 10
      integer, target    :: default_pg_sym(max_mat_number)
      integer, target    :: default_ih_sym(max_mat_number)
      real(8)             :: save_ddot1, save_ddot2, h
!     --------------------------------------------------------------------------

      integer, parameter :: max_response_order_mo = 5
      integer            :: i, j, iz, imat

!#define DEBUG_XC
#ifdef DEBUG_XC
      call matrix_to_file('xc_dmat0', xc_mat_dim*xc_mat_dim, xc_dmat_0)
#endif

!     start timer
      time_integration_start = second()
!     reset timer
      time_ao          = 0.0d0
      time_matrix_dist = 0.0d0
      time_density     = 0.0d0
      time_derv        = 0.0d0

      nr_electrons_integrated = 0.0d0
      xc_energy    = 0.0d0
      tau_integral = 0.0d0

      parallel_xc = (xc_mpi_get_nr_proc() > 1)

#ifdef VAR_MPI
      if (parallel_xc) then
         call xc_mpi_bcast(pure_hf_run)
      end if
#endif
      if (pure_hf_run) then
         write (*,*) 'this is a pure hf run?!?'
         return
      end if

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
      end if

#ifndef MOD_UNRELEASED
      if (do_london_rhs_direct .or. &
          do_london_rhs_ro     .or. &
          do_london_ks         .or. &
          do_london_lr) then
         print *, 'error: LAO+DFT not available in this version'
         stop 1
      end if
      if (do_overlap_diagnostic) then
         print *, 'error: overlap diagnostic not available in this version'
         stop 1
      end if
#endif

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

      if (response_order_mo > max_response_order_mo) then
          stop 'response_order_mo > max_response_order_mo'
      end if

!     set and sync functional
      if (fun_is_automatic) then
!        call xcfun_set_functional(xc_fun,       response_order_mo, parallel_xc)
!        call xcfun_set_functional(xc_fun_alda,  response_order_mo, parallel_xc)
!        call xcfun_set_functional(xc_fun_xalda, response_order_mo, parallel_xc)
!radovan: hack - this is needed otherwise with rhs_direct we go out of the array
         call xcfun_set_functional(xc_fun,       max_response_order_mo, parallel_xc)
         call xcfun_set_functional(xc_fun_alda,  max_response_order_mo, parallel_xc)
         call xcfun_set_functional(xc_fun_xalda, max_response_order_mo, parallel_xc)
      else
#ifdef VAR_MPI
      if (parallel_xc) then
         call dftsyncfunc(xc_mpi_is_master())
      end if
#endif
      end if

#ifndef MOD_UNRELEASED
      if (fun_is_tau_mgga(xc_fun)) then
         print *, 'error: meta-GGA not available in this version'
         stop 1
      end if
#endif

!     check whether gga code structure will be used
      use_gga_qr = (fun_is_gga(xc_fun) .or. fun_is_tau_mgga(xc_fun))
      use_gga_qi = (fun_is_gga(xc_fun) .or. fun_is_tau_mgga(xc_fun))
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
            call my_quit('qr and partial alda is not implemented')
         end if
         if (dft_cfg_alda_ha .and. (.not. dft_cfg_alda_hs)) then
            call my_quit('qr and partial alda is not implemented')
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
     
      call ao_eval_init(max_ao_g_order, &
                        max_ao_m_order, &
                        .true.)

      if (xc_mpi_is_master()) then
!        initial expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot1 = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
         end if

!        get average density from open shells
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) call gtdoav(mat_dim, nz, 1, dmat_0, 1)

#ifdef PRG_DIRAC
         call scale_density_matrices(2.0d0)
#endif

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

      call interface_mo_read()

      if (dft_cfg_saop .or. dft_cfg_grac) then
         call dftac_eigenvalues()
#ifdef VAR_MPI
         !fixme ugly
         if (xc_mpi_is_master()) then
            h = hf_exchange_factor
         end if
         if (parallel_xc) then
            call xc_mpi_bcast(h)
            call set_hf_exchange_factor(h)
         end if
#endif
      end if

!     allocate saop density matrices
      if (dft_cfg_saop) then
         allocate(dmat_saop_outer(mat_dim*mat_dim))
         if (dft_cfg_saop_with_response_part) then
            allocate(dmat_saop_resp(mat_dim*mat_dim))
         end if
         if (mo_eigenvalues_available .and. mo_coef_available) then
            dft_cfg_saop_is_active = .true.
!           only master calculates saop density matrices
            if (xc_mpi_is_master()) then
               call dftac_saop_get_dmats(dmat_saop_outer, &
                                         dmat_saop_resp)
            end if
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

#ifdef DEBUG_XC
      write(*, *) 'debug: XC potential matrix zeroed out prior to integration'
      fmat = 0.0d0
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
         call xc_mpi_reduce(nr_electrons_integrated)
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
#ifdef PRG_DIRAC
         call scale_density_matrices(0.5d0)
#endif

!        get average density from open shells
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) call gtdoav(mat_dim, nz, 1, dmat_0, 1)

!        symmetrize matrices
         if (use_blas3_scheme .or. (use_gga                         &
                                   .and. .not. do_london_rhs_direct &
                                   .and. .not. do_london_ks         &
                                   .and. .not. do_geo)) then
            if (nr_dmat > 0) then
               call gga_sym(mat_dim, nz, nr_fmat, fmat, fmat_pg_sym, dmat_ih_sym)
            else
               call gga_sym(mat_dim, nz, 1, fmat, (/1/), (/1/))
            end if
         end if

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

#ifdef DEBUG_XC
      write(*, *) 'debug: real part of the XC potential matrix'
      call prqmat(xc_fmat,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  1,              &
                  (/1, 2, 3, 4/), &
                  file_unit)
      call matrix_to_file('xc_fmat', mat_dim*mat_dim, xc_fmat)
#endif

!     check for nan in xc matrix
      do imat = 1, xc_nr_fmat
         do iz = 1, nz
            do i = 1, mat_dim
               do j = 1, mat_dim
                  if (xc_fmat(j, i, iz, imat) /= xc_fmat(j, i, iz, imat)) then
                     print *, 'error: found NaN in XC matrix, quitting'
                     stop
                  end if
               end do
            end do
         end do
      end do

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
                  call my_q2bphase('D', iq, 1, dmat(1, 1, iz, im))
               end do
            end do
         end if
         if (nr_fmat > 0) then
            do iz = 1, nz
               do im = 1, nr_fmat
                  iq = pq_to_uq(iz, fmat_pg_sym(im) - 1)
                  call my_q2bphase('F', iq, 1, fmat(1, 1, iz, im))
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

      real(8), allocatable :: ao(:, :)
      real(8), allocatable :: buffer(:, :)

    type(omega_prefactor) :: u

    allocate(u%n (      nr_fmat))
    allocate(u%s (   3, nr_fmat))
    allocate(u%gn(   3, nr_fmat))
    allocate(u%gs(3, 3, nr_fmat))
    allocate(u%tn(      nr_fmat))
    allocate(u%ts(   3, nr_fmat))
    allocate(u%s_lao(3, nr_fmat))
    allocate(u%gs_lao(3, 3, nr_fmat))






    if (use_blas3_scheme) then
       block_threshold = get_block_threshold(dft_cfg_tinydens, mat_dim, dmat_0)
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

    if (do_london_rhs_ro) then
       if (max_fun_derv < 2) then
          max_fun_derv = 2
       end if
    end if

#ifdef MOD_OPENRSP
    if (response_order_ao > 0) then
      derv_length = xc_output_length(xc_fun%id, response_order_ao)
      allocate(derv(derv_length, max_block_length))
      allocate(derv_nadd(derv_length, max_block_length))
    else
#endif
      derv_length = nr_nonzero_derv
      allocate(derv(max_block_length, 0:derv_length))
      allocate(derv_nadd(max_block_length, 0:derv_length))
#ifdef MOD_OPENRSP
    end if
#endif

!       better safe than sorry
        derv      = 0.0d0
        derv_nadd = 0.0d0

!       set it to zero, otherwise it is undefined with lda
        gnn_0 = 0.0d0


   if (use_blas3_scheme) then
      loop_over_points_block: do ipoint = 1, nr_points_on_this_proc, max_block_length
         block_length = min(max_block_length, nr_points_on_this_proc - ipoint + 1)

         allocate(ao(block_length,     nr_ao_slices*nr_ao_cartesian))
         ao = 0.0d0
         allocate(buffer(block_length, nr_ao_slices*nr_ao_cartesian))
         buffer = 0.0d0

         call get_ao(                                  &
                     block_length,                     &
                     rx(ipoint:ipoint+block_length-1), &
                     ry(ipoint:ipoint+block_length-1), &
                     rz(ipoint:ipoint+block_length-1), &
                     ao,                               &
                     buffer                            &
                    )

         call xc_contribution_blas3(                          &
                                    use_gga,                  &
                                    use_gga_qr,               &
                                    use_gga_qi,               &
                                    (response_order_mo == 1), &
                                    block_length,             &
                                    block_threshold,          &
                                    mat_dim,                  &
                                    nz,                       &
                                    buffer,                   &
                                    ao,                       &
                                    nr_fmat,                  &
                                    fmat,                     &
                                    dmat_0,                   &
                                    dmat_saop_outer,          &
                                    dmat_saop_resp,           &
                                    nr_dmat,                  &
                                    dmat,                     &
                                    dmat_pg_sym,              &
                                    dmat_ih_sym,              &
                                    derv_length,              &
                                    derv,                     &
                                    max_fun_derv,             &
                                    rw(ipoint),               &
                                    nr_electrons_integrated,  &
                                    xc_energy                 &
                                   )

#ifdef MOD_UNRELEASED
         if (do_overlap_diagnostic) then
            call integrate_o_matrix(block_length, rw(ipoint), ao)
         end if
#endif

         deallocate(ao)
         deallocate(buffer)

      end do loop_over_points_block
   end if !if (use_blas3_scheme) then

   if (.not. use_blas3_scheme) then

      allocate(ao(1,     nr_ao_slices*nr_ao_cartesian))
      ao = 0.0d0
      allocate(buffer(1, nr_ao_slices*nr_ao_cartesian))
      buffer = 0.0d0

   print_data_at_point = .true.

!  here starts the expensive loop over points
!  do not evaluate anything inside this loop
!  that does not change from point to point
   do ipoint = 1, nr_points_on_this_proc
   
      time_ao_start = second()
      call get_ao(1,          &
                  rx(ipoint), &
                  ry(ipoint), &
                  rz(ipoint), &
                  ao,         &
                  buffer)
      time_ao = time_ao + second() - time_ao_start

      time_density_start = second()
      if (need_ao_order(1, 0)) then
!        density and density gradient
         call get_gn(n_0(1), gn_0_pointwise, 0, mat_dim, dmat_0, buffer, ao)
         gnn_0(1) = gn_0_pointwise(1)*gn_0_pointwise(1) + gn_0_pointwise(2)*gn_0_pointwise(2) + gn_0_pointwise(3)*gn_0_pointwise(3)
      else
!        density
         call get_n(n_0(1), 0, mat_dim, dmat_0, buffer, ao)
!        for mcscf/srdft also evaluate spin density here
!        call get_s(spin_density, 0, mat_dim, dmat_0, buffer, ao)
      end if
      time_density = time_density + second() - time_density_start

      if (n_0(1) > dft_cfg_tinydens) then
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
         if (fun_is_tau_mgga(xc_fun)) then
            call get_tau(tau,     &
                         0,       &
                         mat_dim, &
                         dmat_0,  &
                         buffer,  &
                         ao)
            tau_integral = tau_integral + w(1)*tau
         end if
#endif /* ifdef MOD_UNRELEASED */

!       number of electrons
        nr_electrons_integrated = nr_electrons_integrated + w(1)*n_0(1)
        
        time_derv_start = second()
#ifdef MOD_OPENRSP
!       if we run openrsp with xcfun we evaluate derivatives
!       somewhere else
        if (response_order_ao == 0) then
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
                                   tau=(/tau/))

#ifdef MOD_OPENRSP
        end if
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
             call dftac(1, derv, w, n_0, (/rho_resp/), (/rho_outer/), gnn_0)
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
         if (fun_is_tau_mgga(xc_fun)) then
            u%tn(1) = derv(1, d0000010)
         end if
#endif

        end if !exchange-correlation potential



#ifdef MOD_UNRELEASED
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
                                           ao, &
                                           rx(ipoint), ry(ipoint), rz(ipoint))
            else
              call lda_london_distribute_r(derv(1, d1000000),                 &
                                           fmat_pg_sym,                       &
                                           mat_dim,                           &
                                           nz,                                &
                                           fmat,                              &
                                           ao,                                &
                                           rx(ipoint), ry(ipoint), rz(ipoint))
            end if

          end if


!       london reorthogonalization spin density contribution
!       ====================================================

        if (do_london_rhs_ro .and. .not. dft_cfg_no_sdft) then
           call lr_sdft_lao_reorth(bllen,          &
                                   mat_dim,        &
                                   nz,             &
                                   ao,             &
                                   n_0,            &
                                   gn_0_pointwise, &
                                   gnn_0,          &
                                   w,              &
                                   use_gga_qi,     &
                                   dmat,           &
                                   dmat_pg_sym,    &
                                   buffer,         &
                                   derv_length,    &
                                   derv, u)
        end if


!       london direct spin density contribution
!       =======================================

        if (do_london_rhs_direct .and. .not. dft_cfg_no_sdft) then
!radovan: this code is commented and tests pass, strange
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
        end if
#endif /* #ifdef MOD_UNRELEASED */


!       linear response
!       ===============

        if (response_order_mo == 1) then
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
                   gn_0_pointwise, &
                   gnn_0,          &
                   w,              &
                   use_gga_qr,     &
                   use_gga_qi,     &
                   buffer,         &
                   derv_length,    &
                   derv, u)
        endif


!       quadratic response
!       ==================
        if (response_order_mo == 2) then

          call sdft_qr(bllen,           &
                       use_gga,                &
                       dft_cfg_no_sdft,        &
                       mat_dim, dmat,        &
                       dmat_pg_sym,dmat_ih_sym,ao,&
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
                                          ao, &
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

          if (do_geo_0) then
             if (use_gga) then
                call xc_geo_0_gga_c1(property_gradient,               &
                                     derv(1, d1000000),                    &
                                     derv(1, d0010000),                    &
                                     gn_0_pointwise,                             &
                                     mat_dim,                         &
                                     dmat_0,                          &
                                     ao)
             else
                call xc_geo_0_lda_c1(property_gradient, &
                                     derv(1, d1000000),      &
                                     mat_dim,           &
                                     dmat_0,            &
                                     ao)
             end if
          end if

#ifdef MOD_UNRELEASED
          if (do_geo_1) then
             call get_n(r_t1(1), 0, mat_dim, tmat_1, buffer, ao)
             call xc_geo_1_lda_c1(property_gradient, &
                                  derv(1, d1000000),      &
                                  derv(1, d2000000),     &
                                  r_t1(1),           &
                                  mat_dim,           &
                                  dmat_0,            &
                                  tmat_1,            &
                                  ao)
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
                                  ao)
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
                                      ao)
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
                                      ao)
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
                                ao(1, ao_off_g1_m0(1, 0) + 1))
                call omega_real(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                0.5d0*u%tn(im),      &
                                ao(1, ao_off_g1_m0(2, 0) + 1))
                call omega_real(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                0.5d0*u%tn(im),      &
                                ao(1, ao_off_g1_m0(3, 0) + 1))
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
                                ao(1, ao_off_g1_m0(1, 0) + 1))
                call omega_imag(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                u%ts(1, im),         &
                                ao(1, ao_off_g1_m0(2, 0) + 1))
                call omega_imag(fmat,                &
                                im,                  &
                                fmat_pg_sym(im) - 1, &
                                mat_dim,             &
                                nz,                  &
                                u%ts(1, im),         &
                                ao(1, ao_off_g1_m0(3, 0) + 1))
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
                                      ao)
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

      end if !if (n_0(1) > dft_cfg_tinydens) then
    end do !do ipoint = 1, nr_points_on_this_proc

    deallocate(ao)
    deallocate(buffer)

   end if !if (.not. use_blas3_scheme) then

    deallocate(derv)
    deallocate(derv_nadd)

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

      error_nr_electrons = nr_electrons_integrated - nr_electrons

      write(file_unit, *)
      write(file_unit, *) 'DFT integration'
      write(file_unit, *)

      write(file_unit, '(3x, a, i9)')                            &
            'number of grid points                          = ', &
            nr_points_total

      if (do_potential) then
         write(file_unit, '(3x, a, f26.16)')                        &
               'DFT exchange-correlation energy:               = ', &
               xc_energy
      end if

      write(file_unit, '(3x, a, f26.16)')                        &
            'number of electrons from numerical integration = ', &
            nr_electrons_integrated

      if (abs(error_nr_electrons) > 1.0d-3) then
         write(file_unit, '(3x, a, i9)')                            &
               'number of electrons from orbital occupations   = ', &
               nr_electrons
         write(file_unit, '(3x, a, f26.16)')                        &
               'WARNING: error in the number of electrons      = ', &
               error_nr_electrons
         write(file_unit, '(12x, a, f15.6)')                        &
               'is larger than 1.0d-3'
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




   function get_nr_electrons_integrated()
      real(8) :: get_nr_electrons_integrated
      get_nr_electrons_integrated = nr_electrons_integrated
   end function

   function get_xc_energy()
      real(8) :: get_xc_energy
      get_xc_energy = xc_energy
   end function

   function get_xc_mat_energy()
      real(8) :: get_xc_mat_energy
      get_xc_mat_energy = xc_mat_energy
   end function

   subroutine loop_over_batches()
    real(8) :: dummy
    integer :: idummy
    integer             :: ibatch
    integer             :: i
    real(8), allocatable :: mpi_buffer(:)

      if (xc_mpi_is_master()) then
         interface_file_open = .false. 
         open(interface_file_unit,       &
! aspg debug
             file   = 'numerical_grid', &
!            file   = 'ADFGRID',        &
              status = 'unknown',        &
              form   = 'formatted',      &
              access = 'sequential')
         rewind(interface_file_unit)
         interface_file_open = .true. 
      end if

      ibatch                = 0
      nr_points_total       = 0
      nr_points_used        = 0

!     loop over batches
      do

!        get number of points in this batch
         if (xc_mpi_is_master()) then
            read(interface_file_unit, *) nr_points_batch
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
            call num_grid_read(rx, ry, rz, rw, interface_file_unit, nr_points_batch)
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

         if (do_potential .or. (response_order_mo == 1)) then
            use_blas3_scheme = .true.
         else
            use_blas3_scheme = .false.
         end if
         if (fun_is_tau_mgga(xc_fun)) then
            use_blas3_scheme = .false.
         end if
         if (dft_cfg_blocked) then
            use_blas3_scheme = .true.
         end if
         if (dft_cfg_pointwise) then
            use_blas3_scheme = .false.
         end if

         call integrate()

         deallocate(rx)
         deallocate(ry)
         deallocate(rz)
         deallocate(rw)

!     end loop over batches
      end do
 999  continue

      if (xc_mpi_is_master()) then
         if (interface_file_open) then
            close(interface_file_unit, status='keep')
         end if
      end if

   end subroutine




  subroutine lr(block_length,       &
                mat_dim,            &
                nz,                 &
                nr_mat,             &
                dmaterturbed,     &
                isym,               &
                ih,                 &
                allow_nonhermitian, &
                ao,               &
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
    real(8), intent(in)    :: ao(*)
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
         if (fun_is_tau_mgga(xc_fun)) then
            call get_tau(t_b,     &
                         isym(im) - 1,             &
                         mat_dim, &
                         dmaterturbed(1 + im_off), &
                         buffer,  &
                         ao)
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
                                        ao)
            else
               call get_gn(r_b,                      &
                           gr_b,                     &
                           isym(im) - 1,             &
                           mat_dim,                  &
                           dmaterturbed(1 + im_off), &
                           buffer,                   &
                           ao)
            end if
         else
            call get_n(r_b,                        &
                       isym(im) - 1,               &
                       mat_dim,                    &
                       dmaterturbed(1 + im_off), &
                       buffer,                     &
                       ao)
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
        if (fun_is_tau_mgga(xc_fun)) then
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
                                     ao)
          else
            call get_gs(s_b,                       &
                        gs_b,                       &
                        isym(im) - 1,               &
                        mat_dim,                    &
                        dmaterturbed(1 + im_off), &
                        buffer,                     &
                        ao)
          end if
        else
        call get_s(s_b,                        &
                   isym(im) - 1,               &
                   mat_dim,                    &
                   dmaterturbed(1 + im_off), &
                   buffer,                     &
                   ao)
        end if
#ifdef MOD_UNRELEASED
        if (fun_is_tau_mgga(xc_fun)) then
           call get_tau_imag(ts_b,     &
                             isym(im) - 1,             &
                             mat_dim, &
                             dmaterturbed(1 + im_off), &
                             buffer,  &
                             ao)
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
           if (fun_is_tau_mgga(xc_fun)) then
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
                                ao,               &
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
    real(8), intent(in)    :: ao(*)
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
      print *, 'ERROR: get_gs_lao needs to be programmed'
      stop
    else
      call get_s_lao(s_b_lao,                            &
                     fmat_pg_sym,                        &
                     mat_dim,                            &
                     dmat_0,                             &
                     ao, &
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
                                ao,               &
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
    real(8), intent(in)    :: ao(*)
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
!                      ao)
!       else
!          call get_n(r_b,                       &
!                     isym(im) - 1,              &
!                     mat_dim,                   &
!                     tb_dmat(1 + im_off),       &
!                     buffer,                    &
!                     ao)
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
                         ao)
          else
             call get_s(s_b,                        &
                        isym(im) - 1,               &
                        mat_dim,                    &
                        tb_dmat(1 + im_off),        &
                        buffer,                     &
                        ao)
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
                         ao, &
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
     real(8) ::          ao(*)
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
     &                            ao)
         call get_gn_nonhermitian(r_c, &
     &                            gr_c, &
     &                            isymop(2) - 1, &
     &                            mat_dim,            &
     &                            dtrmat(1, 1, 1, 2), &
     &                            buf, &
     &                            ao)
      else
         call get_n(r_b, &
     &              isymop(1) - 1, &
     &              mat_dim,            &
     &              dtrmat(1, 1, 1, 1), &
     &              buf, &
     &              ao)
         call get_n(r_c, &
     &              isymop(2) - 1, &
     &              mat_dim,         &
     &              dtrmat(1, 1, 1, 2), &
     &              buf, &
     &              ao)
      end if

      if (.not. no_sdft) then
         if (add_grad) then
            call get_gs_nonhermitian(s_b, &
     &                               gs_b, &
     &                               isymop(1) - 1, &
     &                               mat_dim,            &
     &                               dtrmat(1, 1, 1, 1), &
     &                               buf, &
     &                               ao)
            call get_gs_nonhermitian(s_c, &
     &                               gs_c, &
     &                               isymop(2) - 1, &
     &                               mat_dim,            &
     &                               dtrmat(1, 1, 1, 2), &
     &                               buf, &
     &                               ao)
         else
           call get_s(s_b, &
     &                isymop(1) - 1, &
     &                mat_dim,            &
     &                dtrmat(1, 1, 1, 1), &
     &                buf, &
     &                ao)
           call get_s(s_c, &
     &                isymop(2) - 1, &
     &                mat_dim,            &
     &                dtrmat(1, 1, 1, 2), &
     &                buf, &
     &                ao)
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
