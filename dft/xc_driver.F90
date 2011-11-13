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

module xc_driver

!   written by T. Helgaker
!   DIRAC adaption by T. Saue

#ifdef MOD_EMB
  use fde_mod
#endif
  use dirac_cfg
  use basis_info
  use occupation_info
  use memory_allocator
  use dft_cfg
  use dft_functional_derv
  use dft_ao_in_mem
  use ao_eval
  use density_eval
  use dirac_mpi


  implicit none

   public xc_integrator
   public xc_mpi_wake_up_coworkers
   public xc_mpi_release_coworkers
   public get_xc_energy
   public get_xc_mat_energy
   public reset_grid

   private



   type omega_prefactor
      real(8), allocatable :: n(:)
      real(8), allocatable :: s(:, :)
      real(8), allocatable :: gn(:, :)
      real(8), allocatable :: gs(:, :, :)
      real(8), allocatable :: tn(:)
      real(8), allocatable :: ts(:, :)
   end type




#include "functionals.h"

#include "infpar.h"
#include "dcbgen.h"
#include "dcborb.h"
#ifdef VAR_MPI
#include "mpif.h"
#ifdef INT_STAR8
  dimension     istat(mpi_status_size)
#else
  integer    :: istat(mpi_status_size)
#endif /* ifdef INT_STAR8 */
#endif /* ifdef VAR_MPI */

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
   real(8), allocatable :: ao_block(:, :)
   integer, allocatable :: ao_sort_index(:)
   integer, allocatable :: ao_sort_index_reverse(:)
   real(8), allocatable :: buffer_block(:, :)
   real(8), allocatable :: dmat_sorted(:)
   real(8), allocatable :: xmat_sorted(:)
   real(8), allocatable :: ao_sorted(:)
   real(8), allocatable :: bmat(:)

   real(8) :: tau_integral
   real(8) :: nr_electrons
   real(8) :: xc_energy
   real(8) :: xc_mat_energy
   logical :: grid_is_done   = .false.

   logical :: use_gga_qr = .false. !quaternion real
   logical :: use_gga_qi = .false. !quaternion imaginary
   logical :: use_gga    = .false.

   logical :: is_tau_mgga = .false.

   logical :: i_am_master
   logical :: i_am_coworker
   logical :: parallel_xc
   logical :: coworkers_are_awake = .false.

   integer :: nr_points_process
   integer :: nr_points_part
   integer :: nr_processes
   integer :: nr_points_looped_over
   integer :: nr_points_used
   integer :: nr_points_total
   integer :: nr_points
   integer :: nr_points_rest

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


   logical                      :: do_potential
   integer                      :: response_order_mo
   integer                      :: response_order_ao
   logical                      :: do_london_rhs_direct
   logical                      :: do_london_rhs_ro
   logical                      :: do_geo
   logical                      :: do_geo_0
   logical                      :: do_geo_1
   logical                      :: do_geo_2
   logical                      :: do_overlap_diagnostic

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
      real(8), optional, target, intent(inout) :: xc_fmat(xc_mat_dim, xc_mat_dim, xc_nz, xc_nr_fmat)
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
      logical, optional,         intent(in)    :: xc_do_overlap_diagnostic

      integer, optional,         intent(in)    :: xc_response_order_mo
      integer, optional,         intent(in)    :: xc_response_order_ao
!     --------------------------------------------------------------------------
    real(8), external :: ddot
   integer, parameter :: max_mat_number = 10
   integer, target    :: default_sym(max_mat_number)
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

      i_am_master   = .true.
      i_am_coworker = .false.
      parallel_xc   = .false.

#ifdef VAR_MPI
      if (parcal) then
         call dirac_mpi_bcast(coworkers_are_awake)
         if (coworkers_are_awake) then
!           parallel_xc makes it possible to run
!           instances of xc_integrator in sequential
!           while rest of the program runs in parallel
!           if coworkers are not ready then do not
!           try to integrate in parallel
            parallel_xc = .true.
         end if
         if (mytid == mparid) then
            i_am_master   = .true.
            i_am_coworker = .false.
         else
            i_am_master   = .false.
            i_am_coworker = .true.
         end if
      end if
#endif

      call basis_info_read(i_am_master, parallel_xc)

#ifdef VAR_MPI
      if (parallel_xc) then
         call sync_dft_cfg()
      end if
#endif

      if (i_am_master) then
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
      do_overlap_diagnostic = .false.
      nr_atoms              = 0

      if (i_am_master) then
!        change defaults
         if (present(xc_do_potential))          do_potential          = xc_do_potential
         if (present(xc_response_order_mo))     response_order_mo     = xc_response_order_mo
         if (present(xc_response_order_ao))     response_order_ao     = xc_response_order_ao
         if (present(xc_do_london_rhs_direct))  do_london_rhs_direct  = xc_do_london_rhs_direct
         if (present(xc_do_london_rhs_ro))      do_london_rhs_ro      = xc_do_london_rhs_ro
         if (present(xc_do_geo_0))              do_geo_0              = xc_do_geo_0
         if (present(xc_do_geo_1))              do_geo_1              = xc_do_geo_1
         if (present(xc_do_geo_2))              do_geo_2              = xc_do_geo_2
         if (present(xc_nr_atoms))              nr_atoms              = xc_nr_atoms
         if (present(xc_do_overlap_diagnostic)) do_overlap_diagnostic = xc_do_overlap_diagnostic
      end if

#ifdef VAR_MPI
      if (parallel_xc) then
!        broadcast non-allocatables
         call dirac_mpi_bcast(mat_dim)
         call dirac_mpi_bcast(nz)
         call dirac_mpi_bcast(nr_dmat)
         call dirac_mpi_bcast(nr_fmat)
         call dirac_mpi_bcast(nr_atoms)
         call dirac_mpi_bcast(do_potential)
         call dirac_mpi_bcast(do_london_rhs_direct)
         call dirac_mpi_bcast(do_london_rhs_ro)
         call dirac_mpi_bcast(response_order_mo)
         call dirac_mpi_bcast(response_order_ao)
         call dirac_mpi_bcast(do_geo_0)
         call dirac_mpi_bcast(do_geo_1)
         call dirac_mpi_bcast(do_geo_2)
         call dirac_mpi_bcast(do_overlap_diagnostic)
      end if
#endif

!     nullify(tmat_1)
!     nullify(tmat_2)
      if (do_geo_1) then
         call alloc(tmat_1, mat_dim, mat_dim)
         tmat_1(:, :) = xc_dmat(:, :, 1, 1)
      end if

      if (do_geo_2) then
         call alloc(tmat_1, mat_dim, mat_dim)
         call alloc(tmat_2, mat_dim, mat_dim)
         tmat_1(:, :) = xc_dmat(:, :, 1, 1)
         tmat_2(:, :) = xc_dmat(:, :, 1, 2)
      end if

      nullify(dmat_0)
      if (i_am_coworker) then
         call alloc(dmat_0_container, mat_dim, mat_dim)
         dmat_0 => dmat_0_container
      else
         dmat_0 => xc_dmat_0
      end if

      nullify(dmat)
      if (nr_dmat > 0) then
         if (i_am_coworker) then
            call alloc(dmat_container, mat_dim, mat_dim, nz, nr_dmat)
            dmat => dmat_container
         else
            if (present(xc_dmat)) then
               dmat => xc_dmat
            end if
         end if
      end if

      nullify(fmat)
      if (nr_fmat > 0) then
         if (i_am_coworker) then
            call alloc(fmat_container, mat_dim, mat_dim, nz, nr_fmat)
            fmat => fmat_container
         else
            if (present(xc_fmat)) then
               fmat => xc_fmat
            end if
         end if
      end if

      default_sym = 1

      nullify(fmat_pg_sym)
      nullify(dmat_pg_sym)
      nullify(dmat_ih_sym)

      if (i_am_coworker) then
         if (nr_fmat > 0) then
            call alloc(fmat_pg_sym_container, nr_fmat)
            fmat_pg_sym => fmat_pg_sym_container
         end if
         if (nr_dmat > 0) then
            call alloc(dmat_pg_sym_container, nr_dmat)
            call alloc(dmat_ih_sym_container, nr_dmat)
            dmat_pg_sym => dmat_pg_sym_container
            dmat_ih_sym => dmat_ih_sym_container
         end if
      else
         if (present(xc_fmat_pg_sym)) then
            fmat_pg_sym => xc_fmat_pg_sym
         else
            fmat_pg_sym => default_sym(1:nr_fmat)
         end if
         if (present(xc_dmat_pg_sym)) then
            dmat_pg_sym => xc_dmat_pg_sym
         else
            dmat_pg_sym => default_sym(1:nr_dmat)
         end if
         if (present(xc_dmat_ih_sym)) then
            dmat_ih_sym => xc_dmat_ih_sym
         else
            dmat_ih_sym => default_sym(1:nr_dmat)
         end if
      end if

      do_geo = .false.
      if (do_geo_0 .or. do_geo_1 .or. do_geo_2) then
         do_geo = .true.
      end if
      nullify(property_gradient)
      if (i_am_coworker) then
         if (do_geo) then
            call alloc(property_gradient_container, nr_atoms*3)
            property_gradient => property_gradient_container
            property_gradient = 0.0d0
         end if
      else
         if (present(xc_property_gradient)) then
            property_gradient => xc_property_gradient(1:nr_atoms*3)
            property_gradient = 0.0d0
         end if
      end if

      if (i_am_master) then
!        initial expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot1 = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
         end if

!        get average density from open shells
         if (nr_open_shells > 0) call gtdoav(mat_dim, nz, 1, dmat_0, 1)

         call scale_density_matrices(2.0d0)
         call insert_half_phases()

!        calculate grid
         if (.not. dirac_cfg_fde_response) then
            if (.not. grid_is_done) then
              call dftgrd(dmat_0, 1)
              grid_is_done = .true.
            end if
         else
            write (*,*) ' '
            write (*,*) '   FDE info: requested FDE response. Grid generation will'
            write (*,*) '             be skipped, and impoted grid will be used.'
            write (*,*) ' '
         end if
      end if

#ifdef VAR_MPI
      if (parallel_xc) then
!        broadcast allocatables
         call dirac_mpi_bcast(dmat_0)
         if (nr_fmat > 0) then
!           never bcast fmat
!           only master's fmat contains "HF" contribution
!           otherwise it is multiplied
            call dirac_mpi_bcast(fmat_pg_sym)
         end if
         if (nr_dmat > 0) then
            call dirac_mpi_bcast(dmat)
            call dirac_mpi_bcast(dmat_pg_sym)
            call dirac_mpi_bcast(dmat_ih_sym)
         end if
      end if
#endif

!     allocate saop density matrices
      if (dft_cfg_saop) then
         dft_cfg_saop_is_active = .true.
         call alloc(dmat_saop_outer, mat_dim*mat_dim)
         if (dft_cfg_saop_with_response_part) then
            call alloc(dmat_saop_resp, mat_dim*mat_dim)
         end if
!        only master calculates saop density matrices
         if (i_am_master) then
            call dftac_saop_get_dmats(dmat_saop_outer, &
                                      dmat_saop_resp)
         end if
#ifdef VAR_MPI
         if (parallel_xc) then
            call dirac_mpi_bcast(dmat_saop_outer)
            if (dft_cfg_saop_with_response_part) then
               call dirac_mpi_bcast(dmat_saop_resp)
            end if
         end if
#endif
      end if

!     set and sync functional
#ifdef VAR_MPI
      if (parallel_xc) then
         call dftsyncfunc(i_am_master)
      end if
#endif

!     check whether gga code structure will be used
         use_gga_qr = dftisgga()
         use_gga_qi = dftisgga()
      if (response_order_mo > 0) then
         if (dft_cfg_alda_hs) use_gga_qr = .false.
         if (dft_cfg_alda_ha) use_gga_qi = .false.
      end if
      if (use_gga_qr .or. use_gga_qi) then
         use_gga = .true.
      end if
!     in QR it is not so easy to separate H+ and H-
!     contributions since they may mix
!     therefore i quit here if the user wants QR
!     and not the full ALDA or XALDA
!     it is really unfair to quit as late as here
!     this can be checked earlier
!     but this is anyway an exotic option
      if (response_order_mo > 1) then
         if (dft_cfg_alda_hs .and. (.not. dft_cfg_alda_ha)) then
            call quit('QR and partial ALDA is not implemented')
         end if
         if (dft_cfg_alda_ha .and. (.not. dft_cfg_alda_hs)) then
            call quit('QR and partial ALDA is not implemented')
         end if
      end if


!#define PRINT_FMAT
#ifdef PRINT_FMAT
      write(*, *) 'debug: XC potential matrix zeroed out prior to integration'
      call dzero(fmat, mat_dim*mat_dim*nz)
#endif


      call loop_over_batches()


#ifdef VAR_MPI
      if (parallel_xc) then
         call dft_collect()
      end if
#endif

!     deallocate saop density matrices
      if (dft_cfg_saop) then
        call dealloc(dmat_saop_outer)
        if (dft_cfg_saop_with_response_part) then
          call dealloc(dmat_saop_resp)
        end if
      end if

      if (i_am_master) then
         call scale_density_matrices(0.5d0)

!        get average density from open shells
         if (nr_open_shells > 0) call gtdoav(mat_dim, nz, 1, dmat_0, -1)

!        symmetrize matrices
         if (dft_cfg_blocked .or. (use_gga .and. .not. do_london_rhs_direct .and. .not. do_geo)) then
            if (nr_dmat > 0) then
               call gga_sym(mat_dim, nz, nr_fmat, fmat, fmat_pg_sym, dmat_ih_sym)
            else
               call gga_sym(mat_dim, nz, 1, fmat, (/1/), (/1/))
            end if
         end if


!#define DO_DAXPY
#ifdef DO_DAXPY
!!gosia fixme     if london-direct then antisymmetrize matrices
         if (do_london_rhs_direct) then
            call antisym_london(mat_dim, nz, nr_fmat, fmat, (/1/), (/-1/))
         end if
#endif

         call insert_half_phases()

!        final expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot2    = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
            xc_mat_energy = save_ddot2 - save_ddot1
         end if

         call report()
      end if

      call basis_info_delete()
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
         call dealloc(tmat_1)
      end if
      if (do_geo_2) then
         call dealloc(tmat_1)
         call dealloc(tmat_2)
      end if
!     if (associated(tmat_1)) then
!        nullify(tmat_1)
!     end if
!     if (associated(tmat_2)) then
!        nullify(tmat_2)
!     end if

      if (i_am_coworker) then
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

   subroutine xc_mpi_wake_up_coworkers()

#ifdef VAR_MPI
      if (parcal) then
         coworkers_are_awake = .true.
         call dirac_parctl(4)
      end if
#endif /* ifdef VAR_MPI */

   end subroutine

   subroutine xc_mpi_release_coworkers()

      coworkers_are_awake = .false.

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
      real(8) :: r_0(mxbllen)
      real(8) :: tau
      real(8) :: z_0(mxbllen)
      real(8) :: r_0_t(mxbllen)
      real(8) :: z_0_t(mxbllen)
      real(8) :: w(mxbllen)
!     --------------------------------------------------------------------------
!     geometric business
      real(8) :: r_t1(mxbllen)
      real(8) :: r_t2(mxbllen)
!     --------------------------------------------------------------------------
!     fde and fde response
      logical :: fde_resp_verbose =.false.
      real(8) :: z0tt(3)
!     --------------------------------------------------------------------------
      integer :: id, m
      real(8) :: temp
      real(8) :: block_threshold
      real(8) :: rho13, rho_outer, rho_resp, f
    integer :: ipoint
    integer :: k
      real(8) :: RHG(3),VT(3)
    integer              :: derv_length, im
    real(8), allocatable :: derv(:, :)

!   gosia: temporarily
    real(8) :: u_qi(3,3), s_b(3,3)
    integer :: l

      logical :: un_above_threshold
      logical :: us_above_threshold
      logical :: ugn_above_threshold
      logical :: ugs_above_threshold
      logical :: utn_above_threshold
      logical :: uts_above_threshold

    type(omega_prefactor) :: u

    allocate(u%n (      nr_fmat))
    allocate(u%s (   3, nr_fmat))
    allocate(u%gn(   3, nr_fmat))
    allocate(u%gs(3, 3, nr_fmat))
    allocate(u%tn(      nr_fmat))
    allocate(u%ts(   3, nr_fmat))




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
                      dft_cfg_screening)

    call alloc(ao,     total_ao_vector_length)
    call alloc(buffer, total_ao_vector_length)


    if (dft_cfg_blocked) then
       call alloc(ao_block,              mxbllen, total_ao_vector_length)
       call alloc(buffer_block,          mxbllen, total_ao_vector_length)
       call alloc(ao_sort_index,         mat_dim)
       call alloc(ao_sort_index_reverse, mat_dim)
       call alloc(dmat_sorted,           mat_dim*mat_dim)
       call alloc(xmat_sorted,           mxbllen*mat_dim)
       call alloc(ao_sorted,             mxbllen*total_ao_vector_length)
       call alloc(bmat,                  mxbllen*total_ao_vector_length)
       block_threshold = get_block_threshold(mat_dim, dmat_0)
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


      derv_length = nr_nonzero_derv
      allocate(derv(mxbllen, 0:derv_length))
!       set it to zero, otherwise it is undefined with lda
        z_0 = 0.0d0




!  here starts the expensive loop over points
!  do not evaluate anything inside this loop
!  that does not change from point to point
   do ipoint = 1, nr_points_process
      nr_points_looped_over = nr_points_looped_over + 1

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
         call get_gn(r_0(1), rhg, 0, mat_dim, dmat_0, buffer, ao, ao(ao_start(1, 0):ao_end(1, 0)))
         z_0(1) = rhg(1)*rhg(1) + rhg(2)*rhg(2) + rhg(3)*rhg(3)
      else
!        density
         call get_n(r_0(1), 0, mat_dim, dmat_0, buffer, ao)
!        for mcscf/srdft also evaluate spin density here
!        call get_s(spin_density, 0, mat_dim, dmat_0, buffer, ao)
      end if
      time_density = time_density + second() - time_density_start


      if (r_0(1) > dft_cfg_dfthr0) then
         nr_points_used = nr_points_used + 1

         u%n  = 0.0d0
         u%s  = 0.0d0
         u%gn = 0.0d0
         u%gs = 0.0d0
         u%tn = 0.0d0
         u%ts = 0.0d0

         w(1) = rw(ipoint)


!       number of electrons
        nr_electrons = nr_electrons + w(1)*r_0(1)

        time_derv_start = second()
              call get_functional_derv(xc_fun,       &
                                       max_fun_derv, &
                                       bllen,        &
                                       w,            &
                                       r_0,          &
                                       z_0,          &
                                       derv_length,  &
                                       derv,         &
                                       alda_real=dft_cfg_alda_hs, &
                                       alda_imag=dft_cfg_alda_ha, &
                                       tau=tau)
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
            rho13 = r_0(1)**(1.0d0/3.0d0)
            call dftac(derv, w(1), r_0(1), rho13, &
     &                 rho_resp, rho_outer, z_0(1))
          end if

!         energy density
          xc_energy = xc_energy + derv(1, d0000000)

          u%n(1) = derv(1, d1000000)
          if (use_gga) then
!           this factor has nothing to do
!           with symmetrization later
!           it comes from the fact that derivative of z
!           yields two identical terms (one left, one right)
            u%gn(:, 1) = 2.0d0*derv(1, d0010000)*rhg
          end if

        end if !exchange-correlation potential



!       direct-lao contribution to differentiated kohn-sham matrix
!       ==========================================================

        if (do_london_rhs_direct .and. response_order_mo == 0) then
          if (dft_cfg_no_sdft) then
!           spin-density contribution for closed-shell systems is 0
            if (use_gga) then
              call gga_london_distribute_r(derv(1, d1000000),                      &
                                           2.0d0*derv(1, d0010000),                      &
                                           rhg,                               &
                                           0,                                 &
                                           mat_dim,                           &
                                           nz,                                &
                                           fmat,                              &
                                           ao(ao_start(0, 0):ao_end(0, 0)), &
                                           ao(ao_start(0, 1):ao_end(0, 1)), &
                                           ao(ao_start(1, 0):ao_end(1, 0)), &
                                           ao(ao_start(1, 1):ao_end(1, 1)), &
                                           rx(ipoint), ry(ipoint), rz(ipoint))
            else
              call lda_london_distribute_r(derv(1, d1000000),                      &
                                           fmat_pg_sym,                       &
                                           mat_dim,                           &
                                           nz,                           &
                                           fmat,                              &
                                           ao,                              &
                                           ao(ao_start(0, 1):ao_end(0, 1)), &
                                           rx(ipoint), ry(ipoint), rz(ipoint))
            end if

          end if

!       gosia: maybe move this to lr() xc_driver_helpers.F90
        else if (do_london_rhs_direct .and. response_order_mo == 1) then
!           spin-density contribution only
          if (.not. dft_cfg_no_sdft)  &
             call quit("error: no_sdft but do_london_rhs_direct and response_order_mo = 1")

!           get "ds_x/dBx ds_y/dBx ds_z/dBx ds_x/dBy ..." to s_b
!           where ds_i/dBj = (R_MN \times r)_j (\chi_mu \Sigma_i \chi_nu) (c_{\mu j}^* c_{\nu j})
          call get_s_lao(s_b,                               &
                         dmat_pg_sym,                       & ! dmat_0_pg_sym
                         mat_dim,                           &
                         dmat_0,                            &
                         ao(ao_start(0, 0):ao_end(0, 0)),   &
                         ao(ao_start(0, 1):ao_end(0, 1)),   &
                         rx(ipoint), ry(ipoint), rz(ipoint))

!         gosia: check whether you multiply correct two numbers (check ordering in derv(,d0200000)
          do l = 1, 3
            do k = 1, 3
              u_qi(k, l) = derv(1, d0200000)*s_b(k, l)
            end do
          end do

          if (use_gga) then
          !  call nabla_omega_imag(...
          else
            do k = 1, 3
               call omega_imag(fmat,                &              
                               1,                   &
                               fmat_pg_sym(k) - 1,  & ! i hope pg sym of "Bx", "By" or "Bz" contribution, gosia: to check
                               mat_dim,             &
                               nz,                  &
                               u_qi(1, k),          &
                               ao)
            end do
          end if
        end if



!       linear response
!       ===============
        if (response_order_mo == 1) then

!         gosia: todo use this call to add lao+sdft
          call lr(bllen,   &
                  mat_dim,        &
                  nz,             &
                  nr_dmat,      &
                  dmat,         &
                  dmat_pg_sym,    &
                  dmat_ih_sym,         &
                  .false.,        &
                  ao,           &
                  r_0,            &
                  rhg,            &
                  z_0,            &
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
                       dmat_pg_sym,dmat_ih_sym,ao,ao(ao_start(1, 0)), &
                       r_0(1),z_0(1),rhg,w(1), &
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
                  r_0,             &
                  rhg,             &
                  z_0,             &
                  w,               &
                  use_gga_qr,      &
                  use_gga_qi,      &
                  buffer,          &
                  derv_length,     &
                  derv, u)

        endif






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
       
       
          end do

          time_matrix_dist = time_matrix_dist + second() - time_matrix_dist_start

      end if !if (r_0(1) > dft_cfg_dfthr0) then
    end do !do ipoint = 1, nr_points_process

     deallocate(derv)

    call dealloc(ao)
    call dealloc(buffer)
    if (dft_cfg_blocked) then
       call dealloc(ao_block)
       call dealloc(buffer_block)
       call dealloc(ao_sort_index)
       call dealloc(ao_sort_index_reverse)
       call dealloc(dmat_sorted)
       call dealloc(xmat_sorted)
       call dealloc(ao_sorted)
       call dealloc(bmat)
    end if

    deallocate(u%n)
    deallocate(u%s)
    deallocate(u%gn)
    deallocate(u%gs)
    deallocate(u%tn)
    deallocate(u%ts)

   end subroutine

   subroutine report()
    real(8) :: error_nr_electrons

      call occupation_info_read()

      error_nr_electrons  = nr_electrons &
                          - nr_electrons_total

      write(file_unit, *)

!     verify that the number of points we have actually looped
!     over indeed matches the expected total number of points
!     if not something went wrong in the distribution
!     if ((nr_points_total - nr_points_looped_over) /= 0) then
!        call quit('parallel dft: error in distribution of points')
!     end if

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
            nr_processes

      write(file_unit, *)

   end subroutine




#ifdef VAR_MPI
   subroutine dft_collect()
      integer :: ierr
    real(8) :: sum_r
    integer :: sum_i

      call mpi_reduce(nr_electrons,         &
                      sum_r,                &
                      1,                    &
                      mpi_double_precision, &
                      mpi_sum,              &
                      0,                    &
                      mpi_comm_world,       &
                      ierr)
      if (i_am_master) then
         nr_electrons = sum_r
      end if

      call mpi_reduce(nr_points_looped_over, &
                      sum_i,                 &
                      1,                     &
                      mpi_integer,           &
                      mpi_sum,               &
                      0,                     &
                      mpi_comm_world,        &
                      ierr)
      if (i_am_master) then
         nr_points_looped_over = sum_i
      end if

      call mpi_reduce(nr_points_used,        &
                      sum_i,                 &
                      1,                     &
                      mpi_integer,           &
                      mpi_sum,               &
                      0,                     &
                      mpi_comm_world,        &
                      ierr)
      if (i_am_master) then
         nr_points_used = sum_i
      end if

      if (do_potential) then
         call mpi_reduce(xc_energy,            &
                         sum_r,                &
                         1,                    &
                         mpi_double_precision, &
                         mpi_sum,              &
                         0,                    &
                         mpi_comm_world,       &
                         ierr)
         if (i_am_master) then
            xc_energy = sum_r
         end if
      end if

      if (nr_fmat > 0) then
         if (i_am_master) then
            call alloc(fmat_container, mat_dim, mat_dim, nz, nr_fmat)
         end if

         call mpi_reduce(fmat,                       &
                         fmat_container,             &
                         mat_dim*mat_dim*nz*nr_fmat, &
                         mpi_double_precision,       &
                         mpi_sum,                    &
                         0,                          &
                         mpi_comm_world,             &
                         ierr)

         if (i_am_master) then
            call dcopy(mat_dim*mat_dim*nz*nr_fmat, fmat_container, 1, fmat, 1)
            call dealloc(fmat_container)
         end if
      end if

      if (do_geo) then
         if (i_am_master) then
            call alloc(property_gradient_container, nr_atoms*3)
         end if

         call mpi_reduce(property_gradient,           &
                         property_gradient_container, &
                         nr_atoms*3,                   &
                         mpi_double_precision,         &
                         mpi_sum,                      &
                         0,                            &
                         mpi_comm_world,               &
                         ierr)

         if (i_am_master) then
            call dcopy(nr_atoms*3, property_gradient_container, 1, property_gradient, 1)
            call dealloc(property_gradient_container)
         end if
      end if

   end subroutine
#endif /* ifdef VAR_MPI */

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

   subroutine reset_grid()

      grid_is_done = .false.

   end subroutine

   subroutine loop_over_batches()
      integer :: ierr
    real(8) :: dummy
    integer :: idummy
    integer             :: ibatch
    integer             :: i
    real(8), allocatable :: mpi_buffer(:)

      if (i_am_master) then
         if (dirac_cfg_fde_response) then
!           don't open and rewind anything
!           data is already stored at this point
         else
            rewind file_unit_grid
            if (dft_cfg_multigrid) then
               rewind dft_cfg_multigrid_file_unit(response_order_ao)
            end if
         end if
      end if

      ibatch                = 0
      nr_points_total       = 0
      nr_points_looped_over = 0
      nr_points_used        = 0

!     loop over batches
      do

!        get number of points in this batch
         if (i_am_master) then
            if (dft_cfg_import_grid) then
               read(file_unit_grid, *) nr_points
            else if (dft_cfg_multigrid) then
               read(dft_cfg_multigrid_file_unit(response_order_ao), *) nr_points
            else
               read(file_unit_grid) nr_points
            end if
         end if

#ifdef VAR_MPI
         if (parallel_xc) then
            call dirac_mpi_bcast(nr_points)
         end if
#endif
         if (nr_points < 0) exit

         nr_points_process = nr_points
         nr_processes      = 1
#ifdef VAR_MPI
         if (parallel_xc) then
             nr_processes      = nodes + 1
             nr_points_part    = nr_points/nr_processes
             nr_points_rest    = mod(nr_points, nr_processes)
             nr_points_process = nr_points_part
!            the last process gets the rest - these points are usually
!            the most diffuse ones and therefore the most screened ones
             if (mytid == nodes) nr_points_process = nr_points_part + nr_points_rest
         end if
#endif

         if (i_am_master) then
            call alloc(rx, nr_points)
            call alloc(ry, nr_points)
            call alloc(rz, nr_points)
            call alloc(rw, nr_points)
         else
            call alloc(rx, nr_points_process)
            call alloc(ry, nr_points_process)
            call alloc(rz, nr_points_process)
            call alloc(rw, nr_points_process)
         end if

         if (i_am_master) then
            if (dft_cfg_import_grid) then
               call import_dft_grid(rx, ry, rz, rw, file_unit_grid, nr_points)
            else if (dft_cfg_multigrid) then
               call import_dft_grid(rx, ry, rz, rw, dft_cfg_multigrid_file_unit(response_order_ao), nr_points)
            else
               call reaqua(rx, ry, rz, rw, file_unit_grid, nr_points)
            end if
         end if

         ibatch          = ibatch + 1
         nr_points_total = nr_points_total + nr_points

#ifdef VAR_MPI
         if (parallel_xc) then
            if (i_am_master) then
               call alloc(mpi_buffer, nr_points_process)
               call mpi_scatter(rx, nr_points_process,                                       &
                                mpi_real8, mpi_buffer, nr_points_process, mpi_real8, mparid, &
                                mpi_comm_world, ierr)
               call mpi_scatter(ry, nr_points_process,                                       &
                                mpi_real8, mpi_buffer, nr_points_process, mpi_real8, mparid, &
                                mpi_comm_world, ierr)
               call mpi_scatter(rz, nr_points_process,                                       &
                                mpi_real8, mpi_buffer, nr_points_process, mpi_real8, mparid, &
                                mpi_comm_world, ierr)
               call mpi_scatter(rw, nr_points_process,                                       &
                                mpi_real8, mpi_buffer, nr_points_process, mpi_real8, mparid, &
                                mpi_comm_world, ierr)
               call dealloc(mpi_buffer)
               if (nr_points_rest /= 0) then
                  i = nr_points - nr_points_rest + 1
                  call mpi_send(rx(i), nr_points_rest, &
                                mpi_real8, nodes, 90, mpi_comm_world, ierr)
                  call mpi_send(ry(i), nr_points_rest, &
                                mpi_real8, nodes, 91, mpi_comm_world, ierr)
                  call mpi_send(rz(i), nr_points_rest, &
                                mpi_real8, nodes, 92, mpi_comm_world, ierr)
                  call mpi_send(rw(i), nr_points_rest, &
                                mpi_real8, nodes, 93, mpi_comm_world, ierr)
               end if
            end if

            if (i_am_coworker) then
               call mpi_scatter(dummy, idummy,                                       &
                                mpi_real8, rx, nr_points_process, mpi_real8, mparid, &
                                mpi_comm_world, ierr)
               call mpi_scatter(dummy, idummy,                                       &
                                mpi_real8, ry, nr_points_process, mpi_real8, mparid, &
                                mpi_comm_world, ierr)
               call mpi_scatter(dummy, idummy,                                       &
                                mpi_real8, rz, nr_points_process, mpi_real8, mparid, &
                                mpi_comm_world, ierr)
               call mpi_scatter(dummy, idummy,                                       &
                                mpi_real8, rw, nr_points_process, mpi_real8, mparid, &
                                mpi_comm_world, ierr)
               if ((nr_points_rest /= 0) .and. (mytid == nodes)) then
                  i = nr_points_part + 1
                  call mpi_recv(rx(i), nr_points_rest, &
                                mpi_real8, mparid, 90, mpi_comm_world, istat, ierr)
                  call mpi_recv(ry(i), nr_points_rest, &
                                mpi_real8, mparid, 91, mpi_comm_world, istat, ierr)
                  call mpi_recv(rz(i), nr_points_rest, &
                                mpi_real8, mparid, 92, mpi_comm_world, istat, ierr)
                  call mpi_recv(rw(i), nr_points_rest, &
                                mpi_real8, mparid, 93, mpi_comm_world, istat, ierr)
               end if
            end if
         end if
#endif

         call integrate()

         call dealloc(rx)
         call dealloc(ry)
         call dealloc(rz)
         call dealloc(rw)


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
                r_0,                &
                gr_0,               &
                z_0,                &
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
    real(8), intent(in)    :: r_0(*)
    real(8), intent(in)    :: gr_0(*)
    real(8), intent(in)    :: z_0(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga_qr
    logical, intent(in)    :: is_gga_qi
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv(mxbllen, 0:n)
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
           z_b = 2.0d0*(gr_0(1)*gr_b(1) &
                      + gr_0(2)*gr_b(2) &
                      + gr_0(3)*gr_b(3))
          
           u%n(im) = u%n(im) + derv(1, d1010000)*z_b
          
           u%gn(1:3, im) = u%gn(1:3, im) + 2.0d0*gr_0(1:3)*r_b*derv(1, d1010000) &
                                         + 2.0d0*gr_0(1:3)*z_b*derv(1, d0020000) &
                                         + 2.0d0*gr_b(1:3)*derv(1, d0010000)
        end if


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

        do k = k_start, 3
           u%s(k, im) = u%s(k, im) + derv(1, d0200000)*s_b(k)
           
           if (is_gga_qi) then
              y_b(k)  = gr_0(1)*gs_b(1, k) &
                      + gr_0(2)*gs_b(2, k) &
                      + gr_0(3)*gs_b(3, k)
           
              u%s(k, im) = u%s(k, im) + derv(1, d0101000)*y_b(k)
           
              u%gs(1:3, k, im) = u%gs(1:3, k, im) +       gr_0(1:3)*s_b(k)*derv(1, d0101000) &
                                                  +       gr_0(1:3)*y_b(k)*derv(1, d0002000) &
                                                  + 2.0d0*gs_b(1:3, k)*derv(1, d0000100)
           end if
        end do

      end if !end of qi contribution

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
                         r_0, &
                         z_0, &
                         gr_0, &
                         w, &
                         buf, &
                         derv_len, &
                         derv, u)

    type(omega_prefactor) :: u
      integer, intent(in) :: block_length
      integer, intent(in)    :: derv_len
      integer, intent(in)    :: mat_dim
      real(8), intent(in)    :: derv(mxbllen, 0:derv_len)
      logical   add_grad, &
     &          no_sdft
     real(8) ::          dtrmat(mat_dim,mat_dim,nz,3)
     integer ::          isymop(3)
     integer ::          ihrmop(3)
     real(8) ::          gao(*)
     real(8) ::          gao1(*)
     real(8) ::          buf(*)
     real(8) ::          r_0(*)
     real(8) ::          z_0(*)
     real(8) ::          w(*)
!
     real(8) ::          gr_0 (  3)
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


      real(8) ::  v_bc_h_z_0
      real(8) ::  s_bc_h
      real(8) ::  v_bc_h_z_b
      real(8) ::  v_bc_h_z_c
      real(8) ::  z_b, z_c, z_bc






!
!               here we construct (see thesis for notation):
!
!               k_{pq,0}^(bc) = r_0^(bc) \omega_{pq,0} + q_{0j}^(bc) \nabla_j \omega_{pq,0}
!               k_{pq,z}^(bc) = r_z^(bc) \omega_{pq,z} + q_{zj}^(bc) \nabla_j \omega_{pq,z}
!               k_{pq,y}^(bc) = r_y^(bc) \omega_{pq,y} + q_{yj}^(bc) \nabla_j \omega_{pq,y}
!               k_{pq,x}^(bc) = r_x^(bc) \omega_{pq,x} + q_{xj}^(bc) \nabla_j \omega_{pq,x}
!
!               all the rest is done elsewhere
!
!               scalar h+ part: r_0^[bc]
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
        z_b  = gr_0(1)*gr_b(1) &
             + gr_0(2)*gr_b(2) &
             + gr_0(3)*gr_b(3)
        z_c  = gr_0(1)*gr_c(1) &
             + gr_0(2)*gr_c(2) &
             + gr_0(3)*gr_c(3)
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
        y_b(k)  = gr_0(1)*gs_b(1, k) &
                + gr_0(2)*gs_b(2, k) &
                + gr_0(3)*gs_b(3, k)
        y_c(k)  = gr_0(1)*gs_c(1, k) &
                + gr_0(2)*gs_c(2, k) &
                + gr_0(3)*gs_c(3, k)
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

        v_bc_h_z_0 = 0.0d0
        v_bc_h_z_b = 0.0d0
        v_bc_h_z_c = 0.0d0

        v_bc_h_z_0 = v_bc_h_z_0 &
     &             + derv(1, d0030000)*(z_b*z_c) &
     &             + derv(1, d1020000)*(r_b*z_c + z_b*r_c) &
     &             + derv(1, d2010000)*(r_b*r_c)

        v_bc_h_z_0 = v_bc_h_z_0 &
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
            v_bc_h_z_0  = v_bc_h_z_0 &
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
          v_bc_h(j) = v_bc_h(j) + v_bc_h_z_0*gr_0(j)*2.0d0 &
     &                          + v_bc_h_z_b*gr_b(j)*2.0d0 &
     &                          + v_bc_h_z_c*gr_c(j)*2.0d0
          if(.not. no_sdft) then
            do k = k_start,3
              v_bc_ah(j,k) = v_bc_ah(j,k) + v_bc_h_y_0(k)*gr_0(j) &
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
