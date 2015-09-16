module fde_dirac_matrices_integration

!embed modules
   use fde_types
   use fde_cfg
   use fde_io
   use fde_data
   use fde_nadd_derv
   use fde_xcfun_interface
   use fde_max_block_length

! dirac-specific modules
   use interface_ao
   use interface_mo
   use interface_file_io
   use interface_grid
   use dirac_ao_eval
   use density_eval
   use fde_mpi

   implicit none

   public fde_dirac_emb_matrices_via_integration
   public fde_dirac_emb_pot_energy

   private

   type fde_omega_prefactor
      real(8), allocatable :: n(:)
      real(8), allocatable :: gn(:, :)
      real(8), allocatable :: s(:, :)
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

   real(8), allocatable :: rx(:)
   real(8), allocatable :: ry(:)
   real(8), allocatable :: rz(:)
   real(8), allocatable :: rw(:)

   real(8), allocatable :: frozen_n(:)
   real(8), allocatable :: frozen_gnx(:)
   real(8), allocatable :: frozen_gny(:)
   real(8), allocatable :: frozen_gnz(:)
   real(8), allocatable :: frozen_elpot(:)
   
   real(8), allocatable :: dmat_sorted(:, :)

   real(8) :: nr_electrons_integrated
   real(8), save :: fde_energy
   real(8), save :: fde_energy_density
   real(8), save :: fde_mat_energy

   logical :: use_gga_qr = .false. !quaternion real
   logical :: use_gga_qi = .false. !quaternion imaginary
   logical :: use_gga    = .false.

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
   logical :: do_fde_update_vemb
   logical :: do_fde_response

   integer, allocatable, target :: default_fmat_pg_sym(:)
   integer, allocatable, target :: default_dmat_pg_sym(:)
   integer, allocatable, target :: default_dmat_ih_sym(:)

contains

    function fde_dirac_emb_pot_energy()
       real(8) :: fde_dirac_emb_pot_energy 
       fde_dirac_emb_pot_energy = fde_mat_energy
    end function


!   ----------------------------------------------------------------------------
    subroutine fde_dirac_emb_matrices_via_integration(  &
                             fde_mat_dim,               &
                             fde_nz,                    &
                             fde_dmat_0,                &
                             fde_nr_dmat,               &
                             fde_nr_fmat,               &
                             fde_dmat,                  &
                             fde_fmat,                  &
                             fde_dmat_ih_sym,           &
                             fde_dmat_pg_sym,           &
                             fde_fmat_pg_sym,           &
                             fde_nr_atoms,              &
                             fde_property_gradient,     &
                             fde_do_potential,          &
                             fde_do_london_rhs_direct,  &
                             fde_do_london_rhs_ro,      &
                             fde_do_geo_0,              &
                             fde_do_geo_1,              &
                             fde_do_geo_2,              &
                             fde_do_london_ks,          &
                             fde_do_london_lr,          &
                             fde_response_order_mo,     &
                             fde_response_order_ao)

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
      integer,                   intent(in)    :: fde_mat_dim
      integer,                   intent(in)    :: fde_nz
      real(8),           target, intent(in)    :: fde_dmat_0(fde_mat_dim, fde_mat_dim)
      integer,                   intent(in)    :: fde_nr_dmat
      integer,                   intent(in)    :: fde_nr_fmat
      real(8), optional, target, intent(in)    :: fde_dmat(fde_mat_dim, fde_mat_dim, fde_nz, fde_nr_dmat)
!     real(8), optional, target, intent(inout) :: fde_fmat(fde_mat_dim, fde_mat_dim, fde_nz, fde_nr_fmat)
      real(8), optional, target                :: fde_fmat(fde_mat_dim, fde_mat_dim, fde_nz, fde_nr_fmat)
      integer, optional, target, intent(in)    :: fde_dmat_ih_sym(fde_nr_dmat)
      integer, optional, target, intent(in)    :: fde_dmat_pg_sym(fde_nr_dmat)
      integer, optional, target, intent(in)    :: fde_fmat_pg_sym(fde_nr_fmat)

      integer, optional,         intent(in)    :: fde_nr_atoms
      real(8), optional, target, intent(inout) :: fde_property_gradient(*)

      logical, optional,         intent(in)    :: fde_do_potential
      logical, optional,         intent(in)    :: fde_do_london_rhs_direct
      logical, optional,         intent(in)    :: fde_do_london_rhs_ro
      logical, optional,         intent(in)    :: fde_do_geo_0
      logical, optional,         intent(in)    :: fde_do_geo_1
      logical, optional,         intent(in)    :: fde_do_geo_2
      logical, optional,         intent(in)    :: fde_do_london_ks
      logical, optional,         intent(in)    :: fde_do_london_lr

      integer, optional,         intent(in)    :: fde_response_order_mo
      integer, optional,         intent(in)    :: fde_response_order_ao
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
      call matrix_to_file('fde_dmat0', fde_mat_dim*fde_mat_dim, fde_dmat_0)
#endif

!     start timer
      time_integration_start = second()
!     reset timer
      time_ao          = 0.0d0
      time_matrix_dist = 0.0d0
      time_density     = 0.0d0
      time_derv        = 0.0d0

      nr_electrons_integrated = 0.0d0
      fde_energy    = 0.0d0
      fde_energy_density = 0.0d0

      parallel_fde = (fde_mpi_get_nr_proc() > 1)

#ifdef VAR_MPI
      if (parallel_fde) then
         call sync_fde_cfg()
      end if
#endif

      if (fde_mpi_is_master()) then
         mat_dim = fde_mat_dim
         nz      = fde_nz
         nr_dmat = fde_nr_dmat
         nr_fmat = fde_nr_fmat
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

      nr_atoms              = 0

     if (fde_mpi_is_master()) then      
!        change defaults
         if (present(fde_do_potential))          do_potential          = fde_do_potential
         if (present(fde_response_order_mo))     response_order_mo     = fde_response_order_mo
         if (present(fde_response_order_ao))     response_order_ao     = fde_response_order_ao
         if (present(fde_do_london_rhs_direct)) do_london_rhs_direct  = fde_do_london_rhs_direct
         if (present(fde_do_london_rhs_ro))      do_london_rhs_ro      = fde_do_london_rhs_ro
         if (present(fde_do_geo_0))              do_geo_0              = fde_do_geo_0
         if (present(fde_do_geo_1))              do_geo_1              = fde_do_geo_1
         if (present(fde_do_geo_2))              do_geo_2              = fde_do_geo_2
         if (present(fde_do_london_ks))          do_london_ks          = fde_do_london_ks
         if (present(fde_do_london_lr))          do_london_lr          = fde_do_london_lr
         if (present(fde_nr_atoms))              nr_atoms              = fde_nr_atoms

      end if

#ifdef VAR_MPI
      if (parallel_fde) then
!        broadcast non-allocatables
         call fde_mpi_bcast(mat_dim)
         call fde_mpi_bcast(nz)
         call fde_mpi_bcast(nr_dmat)
         call fde_mpi_bcast(nr_fmat)
         call fde_mpi_bcast(nr_atoms)
         call fde_mpi_bcast(do_potential)
         call fde_mpi_bcast(do_london_rhs_direct)
         call fde_mpi_bcast(do_london_rhs_ro)
         call fde_mpi_bcast(response_order_mo)
         call fde_mpi_bcast(response_order_ao)
         call fde_mpi_bcast(do_geo_0)
         call fde_mpi_bcast(do_geo_1)
         call fde_mpi_bcast(do_geo_2)
         call fde_mpi_bcast(do_london_ks)
         call fde_mpi_bcast(do_london_lr)

      end if
#endif

!     nullify(tmat_1)
!     nullify(tmat_2)
      if (do_geo_1) then
         allocate(tmat_1(mat_dim, mat_dim))
         tmat_1(:, :) = fde_dmat(:, :, 1, 1)
      end if

      if (do_geo_2) then
         allocate(tmat_1(mat_dim, mat_dim))
         allocate(tmat_2(mat_dim, mat_dim))
         tmat_1(:, :) = fde_dmat(:, :, 1, 1)
         tmat_2(:, :) = fde_dmat(:, :, 1, 2)
      end if

     nullify(dmat_0)
      if (.not. fde_mpi_is_master()) then
         allocate(dmat_0_container(mat_dim, mat_dim))
         dmat_0 => dmat_0_container
      else
         dmat_0 => fde_dmat_0
      end if

      nullify(dmat)
      if (nr_dmat > 0) then
         if (.not. fde_mpi_is_master()) then
            allocate(dmat_container(mat_dim, mat_dim, nz, nr_dmat))
            dmat => dmat_container
         else
            if (present(fde_dmat)) then
               dmat => fde_dmat
            end if
         end if
      end if

      nullify(fmat)
      if (nr_fmat > 0) then
         if (.not. fde_mpi_is_master()) then
            allocate(fmat_container(mat_dim, mat_dim, nz, nr_fmat))
            fmat_container = 0.0d0
            fmat => fmat_container
         else
            if (present(fde_fmat)) then
               fmat => fde_fmat
            end if
         end if
      end if

      if (nr_fmat > 0) then
         allocate(default_fmat_pg_sym(nr_fmat))
         default_fmat_pg_sym = 1
      end if
      if (nr_dmat > 0) then
         allocate(default_dmat_pg_sym(nr_dmat))
         allocate(default_dmat_ih_sym(nr_dmat))
         default_dmat_pg_sym = 1
         default_dmat_ih_sym = 0
      end if

      nullify(fmat_pg_sym)
      nullify(dmat_pg_sym)
      nullify(dmat_ih_sym)

      if (.not. fde_mpi_is_master()) then
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
         if (present(fde_fmat_pg_sym)) then
            fmat_pg_sym => fde_fmat_pg_sym
         else
            fmat_pg_sym => default_fmat_pg_sym(1:nr_fmat)
         end if
         if (present(fde_dmat_pg_sym)) then
            dmat_pg_sym => fde_dmat_pg_sym
         else
            dmat_pg_sym => default_dmat_pg_sym(1:nr_dmat)
         end if
         if (present(fde_dmat_ih_sym)) then
            dmat_ih_sym => fde_dmat_ih_sym
         else
            dmat_ih_sym => default_dmat_ih_sym(1:nr_dmat)
         end if
      end if

      do_geo = .false.
      if (do_geo_0 .or. do_geo_1 .or. do_geo_2) then
         do_geo = .true.
      end if

      nullify(property_gradient)
      if (.not. fde_mpi_is_master()) then
         if (do_geo) then
            allocate(property_gradient_container(nr_atoms*3))
            property_gradient => property_gradient_container
            property_gradient = 0.0d0
         end if
      else
         if (present(fde_property_gradient)) then
            property_gradient => fde_property_gradient(1:nr_atoms*3)
            property_gradient = 0.0d0
         end if
      end if

      if (response_order_mo > max_response_order_mo) then
          stop 'response_order_mo > max_response_order_mo'
      end if

! functionals are set
      call xcfun_set_functional(fde_xcf,       max_response_order_mo, parallel_fde)
      call xcfun_set_functional(fde_xcf_alda,  max_response_order_mo, parallel_fde)
      call xcfun_set_functional(fde_xcf_xalda, max_response_order_mo, parallel_fde)
      call xcfun_set_functional(fde_kef,       max_response_order_mo, parallel_fde)
      call xcfun_set_functional(fde_kef_alda,  max_response_order_mo, parallel_fde)
      call xcfun_set_functional(fde_kef_xalda, max_response_order_mo, parallel_fde)

      use_gga_qr = (fun_is_gga(fde_xcf)) .or. (fun_is_gga(fde_kef))
      use_gga_qi = (fun_is_gga(fde_xcf)) .or. (fun_is_gga(fde_kef))
      if (response_order_mo > 0) then
         if (fde_cfg_alda_hs) use_gga_qr = .false.
         if (fde_cfg_alda_ha) use_gga_qi = .false.
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
         if (fde_cfg_alda_hs .and. (.not. fde_cfg_alda_ha)) then
            call fde_quit('qr and partial alda is not implemented')
         end if
         if (fde_cfg_alda_ha .and. (.not. fde_cfg_alda_hs)) then
            call fde_quit('qr and partial alda is not implemented')
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

      call dirac_ao_eval_init(max_ao_g_order, &
                        max_ao_m_order, &
                        .true.)

      if (fde_mpi_is_master()) then
!        initial expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot1 = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
         end if

!        get average density from open shells
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) &
            call gtdoav(mat_dim, nz, 1, dmat_0, 1)

#ifdef PRG_DIRAC
         call scale_density_matrices(2.0d0)
#endif

         call insert_half_phases()

      end if

#ifdef VAR_MPI
      if (parallel_fde) then
!        broadcast allocatables
         print *, "before dmat_0 bcast for process ",fde_get_my_rank()
         call fde_mpi_bcast(dmat_0)
         print *, "after dmat_0 bcast for process ",fde_get_my_rank()
         if (nr_fmat > 0) then
!           never bcast fmat
!           only master's fmat contains "HF" contribution
!           otherwise it is multiplied
         print *, "before fmat_pg_sym bcast for process ",fde_get_my_rank()
            call fde_mpi_bcast(fmat_pg_sym)
         print *, "after fmat_pg_sym bcast for process ",fde_get_my_rank()
         end if
         if (nr_dmat > 0) then
         print *, "before dmat bcast for process ",fde_get_my_rank()
            call fde_mpi_bcast(dmat)
         print *, "before dmat_pg_sym bcast for process ",fde_get_my_rank()
            call fde_mpi_bcast(dmat_pg_sym)
         print *, "before dmat_ih_sym bcast for process ",fde_get_my_rank()
            call fde_mpi_bcast(dmat_ih_sym)
         end if
      end if
      print *, "after allocatables bcast for process ",fde_get_my_rank()
#endif














































#ifdef DEBUG_XC
      write(*, *) 'debug: XC potential matrix zeroed out prior to integration'
      fmat = 0.0d0
#endif




#ifdef VAR_MPI
!     print *, "before loop_over_batches for process ",fde_get_my_rank()
#endif
      call loop_over_batches()
#ifdef VAR_MPI
!     print *, "before loop_over_batches for process ",fde_get_my_rank()
#endif















#ifdef VAR_MPI
      if (parallel_fde) then
!        collect results from processors
         call fde_mpi_reduce(nr_electrons_integrated)
         call fde_mpi_reduce(nr_points_used)
         if (do_potential) then
            call fde_mpi_reduce(fde_energy)
            call fde_mpi_reduce(fde_energy_density)
         end if
         if (nr_fmat > 0) then
            call fde_mpi_reduce(fmat)
         end if
         if (do_geo) then
            call fde_mpi_reduce(property_gradient)
         end if
      end if
#endif




      if (fde_mpi_is_master()) then
#ifdef PRG_DIRAC
         call scale_density_matrices(0.5d0)
#endif

!        get average density from open shells
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) call gtdoav(mat_dim, nz, 1, dmat_0, 1)

!        symmetrize matrices
         if (fde_use_blas3_scheme .or. (use_gga                         &
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
            fde_mat_energy = save_ddot2 - save_ddot1
         end if

         call report()
      end if

      call nullify_and_release()

#ifdef DEBUG_XC
      write(*, *) 'debug: real part of the XC potential matrix'
      call prqmat(fde_fmat,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  1,              &
                  (/1, 2, 3, 4/), &
                  file_unit)
      call matrix_to_file('fde_fmat', mat_dim*mat_dim, fde_fmat)
#endif

!     check for nan in xc matrix
      do imat = 1, fde_nr_fmat
         do iz = 1, nz
            do i = 1, mat_dim
               do j = 1, mat_dim
                  if (fde_fmat(j, i, iz, imat) /= fde_fmat(j, i, iz, imat)) then
                     print *, 'error: found NaN in XC matrix, quitting'
                     stop
                  end if
               end do
            end do
         end do
      end do

   end subroutine


!   ----------------------------------------------------------------------------
   subroutine nullify_and_release()
!   ----------------------------------------------------------------------------

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

      if (.not. fde_mpi_is_master()) then
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

      if (allocated(default_fmat_pg_sym)) then
         deallocate(default_fmat_pg_sym)
      end if
      if (allocated(default_dmat_pg_sym)) then
         deallocate(default_dmat_pg_sym)
      end if
      if (allocated(default_dmat_ih_sym)) then
         deallocate(default_dmat_ih_sym)
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


! ------------------------------------------------------------------------------
   subroutine integrate()
! ------------------------------------------------------------------------------
!        unperturbed density
   real(8) :: n_0(max_block_length)
!        unperturbed density gradient
   real(8) :: gn_0_pointwise(3)

   real(8) :: gnn_0(max_block_length)
   real(8) :: n_0_t(max_block_length)
   real(8) :: gnn_0_t(max_block_length)
   real(8) :: w(max_block_length)
!        --------------------------------------------------------------------------
!     geometric business
   real(8) :: r_t1(max_block_length)
   real(8) :: r_t2(max_block_length)
!     --------------------------------------------------------------------------
!     fde and fde response
!      logical :: fde_resp_verbose =.false.
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
   real(8), allocatable :: derv_active(:, :)
   real(8), allocatable :: derv_total(:, :)

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

   type(fde_omega_prefactor) :: u


   allocate(u%n (      nr_fmat))
   allocate(u%gn(   3, nr_fmat))
   allocate(u%s (   3, nr_fmat))
   allocate(u%gs(3, 3, nr_fmat))
   allocate(u%tn(      nr_fmat))
   allocate(u%ts(   3, nr_fmat))
   allocate(u%s_lao(3, nr_fmat))
   allocate(u%gs_lao(3, 3, nr_fmat))

!  if (fde_use_blas3_scheme) then
!      block_threshold = get_block_threshold(fde_cfg_tinydens, mat_dim, dmat_0)
!   end if

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
    if (fde_cfg_no_sdft) then
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

   derv_length = nr_nonzero_derv
   allocate(derv_active(max_block_length, 0:derv_length))
   allocate(derv_total(max_block_length, 0:derv_length))

!  better safe than sorry
   derv_active = 0.0d0
   derv_total  = 0.0d0

!  set it to zero, otherwise it is undefined with lda
   gnn_0 = 0.0d0

  if (fde_use_blas3_scheme) then
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

!        call fde_contribution_blas3(                          &
!                                   use_gga,                  &
!                                   use_gga_qr,               &
!                                   use_gga_qi,               &
!                                   (response_order_mo == 1), &
!                                   block_length,             &
!                                   block_threshold,          &
!                                   mat_dim,                  &
!                                   nz,                       &
!                                   buffer,                   &
!                                   ao,                       &
!                                   nr_fmat,                  &
!                                   fmat,                     &
!                                   dmat_0,                   &
!                                   dmat_saop_outer,          &
!                                   dmat_saop_resp,           &
!                                   nr_dmat,                  &
!                                   dmat,                     &
!                                   dmat_pg_sym,              &
!                                   dmat_ih_sym,              &
!                                   derv_length,              &
!                                   derv,                     &
!                                   max_fun_derv,             &
!                                   rw(ipoint),               &
!                                   nr_electrons_integrated,  &
!                                   fde_energy_density        &
!                                  )

         deallocate(ao)
         deallocate(buffer)

      end do loop_over_points_block
   end if 

   if (.not. fde_use_blas3_scheme) then

         allocate(ao(1,     nr_ao_slices*nr_ao_cartesian))
         ao = 0.0d0
         allocate(buffer(1, nr_ao_slices*nr_ao_cartesian))
         buffer = 0.0d0

         print_data_at_point = .true.

!  here starts the expensive loop over points
!  do not evaluate anything inside this loop
!  that does not change from point to point
!  print *, "fde nr points",fde_get_my_rank(),nr_points_on_this_proc
   do ipoint = 1, nr_points_on_this_proc
      
!     print *, "fde check points",ipoint,fde_get_my_rank(),ipoint,rx(ipoint),ry(ipoint),rz(ipoint)
      time_ao_start = second()
      call get_ao(1,    &
            rx(ipoint), &
            ry(ipoint), &
            rz(ipoint), &
            ao,         &
            buffer)
      time_ao = time_ao + second() - time_ao_start

!     print *, "fde check ao",ipoint,fde_get_my_rank(),ao
!     print *, "fde check buffer",ipoint,fde_get_my_rank(),buffer

      time_density_start = second()
      if (need_ao_order(1, 0)) then
!        density and density gradient
         call get_gn(n_0(1), gn_0_pointwise, 0, mat_dim, dmat_0, buffer, ao)
         gnn_0(1) = gn_0_pointwise(1)*gn_0_pointwise(1) &
                  + gn_0_pointwise(2)*gn_0_pointwise(2) &
                  + gn_0_pointwise(3)*gn_0_pointwise(3)
      else
!        density
         call get_n(n_0(1), 0, mat_dim, dmat_0, buffer, ao)
      endif

      time_density = time_density + second() - time_density_start

!     print *, "fde check density ",ipoint,fde_get_my_rank(),fde_cfg_tinydens,n_0(1),frozen_n(ipoint)

      if (n_0(1) > fde_cfg_tinydens) then
         nr_points_used = nr_points_used + 1

!        print *,"inside ",fde_get_my_rank(),fde_cfg_tinydens,nr_points_used

         u%n  = 0.0d0
         u%gn = 0.0d0
         u%s  = 0.0d0
         u%gs = 0.0d0
         u%tn = 0.0d0
         u%ts = 0.0d0
         u%s_lao  = 0.0d0
         u%gs_lao  = 0.0d0

         w(1) = rw(ipoint)


         if ( frozen_n(ipoint) > 0.0d0 ) then
            n_0_t(1) = n_0(1)            + frozen_n(ipoint)
            z0tt(1)  = gn_0_pointwise(1) + frozen_gnx(ipoint)
            z0tt(2)  = gn_0_pointwise(2) + frozen_gny(ipoint)
            z0tt(3)  = gn_0_pointwise(3) + frozen_gnz(ipoint)

            gnn_0_t(1) = z0tt(1)*z0tt(1) + z0tt(2)*z0tt(2) + z0tt(3)*z0tt(3)

         else
            n_0_t(1)   = n_0(1)
            gnn_0_t(1) = gnn_0(1)
         endif
!       number of electrons
         nr_electrons_integrated = nr_electrons_integrated + w(1)*n_0(1)

         time_derv_start = second()

! aspg: this way we are forcing the same xc functional for the active and non-additive
!       contributions. this is not desired, but we have to overcome xcfun's limitation
!       on using only up to 5 functionals...
!
         call get_xcke_fun_derv(max_fun_derv, & 
                                bllen,        &
                                w,            &
                                n_0_t,        &
                                gnn_0_t,      &
                                derv_total)

         call get_xcke_fun_derv(max_fun_derv, & 
                                bllen,        &
                                w,            &
                                n_0,          &
                                gnn_0,        &
                                derv_active)

         time_derv = time_derv + second() - time_derv_start


         time_matrix_dist_start = second()

!       exchange-correlation potential
!       ==============================

         if (do_potential) then

!         energy will not be calcuated here but instead via
!         call calc_fde_interaction_energy(energy,gf_active,gf_frozen)
!         outside of the integration loop  
            fde_energy_density = fde_energy_density + (derv_total(1, d1000000) - derv_active(1, d1000000)) + w(1)*frozen_elpot(ipoint)

            u%n(1) = (derv_total(1, d1000000) - derv_active(1, d1000000)) + w(1)*frozen_elpot(ipoint) 

         if (use_gga) then
!           this factor has nothing to do
!           with symmetrization later
!           it comes from the fact that derivative of z
!           yields two identical terms (one left, one right)
           u%gn(:, 1) = 2.0d0*( derv_total(1, d0010000)*z0tt - derv_active(1, d0010000)*gn_0_pointwise)
         end if

         end if !exchange-correlation potential


!       direct-lao contribution to differentiated kohn-sham matrix
!       ==========================================================
! N/A

!       london reorthogonalization spin density contribution
!       ====================================================
! N/A

!       london direct spin density contribution
!       =======================================
! N/A

!       linear response
!       ===============

        if (response_order_mo == 1) then
           call lr_fde_uncoupled(  &
                   bllen,          &
                   mat_dim,        &
                   nz,             &
                   nr_dmat,        &
                   dmat,           &
                   dmat_pg_sym,    &
                   dmat_ih_sym,    &
                   .false.,        &
                   ao,             &
                   n_0,            &
                   n_0_t,          &
                   gn_0_pointwise, &
                   z0tt,           &
                   gnn_0,          &
                   gnn_0_t,        &
                   w,              &
                   use_gga_qr,     &
                   use_gga_qi,     &
                   buffer,         &
                   derv_length,    &
                   derv_active,    &
                   derv_total,     &
                   u)
      endif

#ifdef FDE_PLACEHOLDER
!       quadratic response
!       ==================
! N/A
! what goes below is strictly what the standard dft integrator has, so here it's only a placeholder

      if (response_order_mo == 2) then

         call sdft_qr(bllen,           &
                    use_gga,                &
                    fde_cfg_no_sdft,        &
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
         call lr_fde_uncoupled( &
                bllen,          &
                mat_dim,        &
                nz,             &
                nr_dmat,        &
                dmat,           &
                dmat_pg_sym,    &
                dmat_ih_sym,    &
                .false.,        &
                ao,             &
                n_0,            &
                n_0_t,          &
                gn_0_pointwise, &
                z0tt,           &
                gnn_0,          &
                gnn_0_t,        &
                w,              &
                use_gga_qr,     &
                use_gga_qi,     &
                buffer,         &
                derv_length,    &
                derv_active,    &
                derv_total,     &
                u)
 
      endif

!
! all the geo, london stuff below are placeholders, roughly taken from the dft integrator
!
      if (do_geo_0) then
         if (use_gga) then
            call xc_geo_0_gga_c1(property_gradient, &
                                 derv_active(1, d1000000), &
                                 derv_active(1, d0010000), &
                                 derv_total(1, d1000000), &
                                 derv_total(1, d0010000), &
                                 gn_0_pointwise,    &
                                 mat_dim,           &
                                 dmat_0,            &
                                 ao)
         else
            call xc_geo_0_lda_c1(property_gradient, &
                                 derv_active(1, d1000000), &
                                 derv_total(1, d1000000), &
                                 mat_dim,           &
                                 dmat_0,            &
                                 ao)
         end if
      end if

      if (do_geo_1) then
         call get_n(r_t1(1), 0, mat_dim, tmat_1, buffer, ao)
         call xc_geo_1_lda_c1(property_gradient, &
                              derv_active(1, d1000000), &
                              derv_active(1, d2000000), &
                              derv_total(1, d1000000), &
                              derv_total(1, d2000000), &
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
                            derv_active(1, d2000000), &
                            derv_active(1, d3000000), &
                            derv_total(1, d2000000),  &
                            derv_total(1, d3000000),  &
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
          call xc_london_c1_ks_gga( &
                                   derv_active(1, d1000000),        &
                                   derv_active(1, d0010000),        &
                                   derv_total(1, d1000000),         &
                                   derv_total(1, d0010000),         &
                                   gn_0_pointwise,                  &
                                   mat_dim,                         &
                                   nz,                              &
                                   fmat,                            &
                                   ao,                              &
                                   (/rx(ipoint), ry(ipoint), rz(ipoint)/))
       else
          call xc_london_c1_ks_lda( &
                                    derv_active(1, d1000000), &
                                    derv_total(1, d1000000), &
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

#endif


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

! spin-density stuff 

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
          end do

! end spin-density stuff
          time_matrix_dist = time_matrix_dist + second() - time_matrix_dist_start

      end if !if (n_0(1) > fde_cfg_small_density_threshold) then
    end do !do ipoint = 1, nr_points_on_this_proc
    end if

    deallocate(ao)
    deallocate(buffer)
    deallocate(derv_total)
    deallocate(derv_active)

    deallocate(u%n)
    deallocate(u%gn)
    deallocate(u%gs)
    deallocate(u%tn)
    deallocate(u%ts)
    deallocate(u%s_lao)
    deallocate(u%gs_lao)

   end subroutine integrate

   
!   ----------------------------------------------------------------------------
      subroutine report()
!   ----------------------------------------------------------------------------
         real(8) :: error_nr_electrons, fde_nr_electrons = 0

         error_nr_electrons = nr_electrons_integrated - fde_nr_electrons

         write(file_unit, *)
         write(file_unit, *) 'FDE integration'
         write(file_unit, *)
         write(file_unit, '(3x, a, i9)')                            &
            'number of grid points                          = ', &
            nr_points_total
         write(file_unit, '(3x, a, i9)')                            &
            'number of grid points below density threshold  = ', &
            nr_points_total - nr_points_used

      write(file_unit, '(3x, a, f26.16)')                        &
            'number of electrons from numerical integration = ', &
            nr_electrons_integrated

!     if (do_potential) then
!        write(file_unit, '(3x, a, f26.16)')                        &
!              'excess subsystem energy:               = ', &
!              fde_mat_energy
!     end if

      write(file_unit, *)

!     timing
      time_integration = second() - time_integration_start
      call timtxt('  time spent in FDE integration                  =', &
                  time_integration, file_unit)
      write(file_unit, *)

   end subroutine report


!   ----------------------------------------------------------------------------
   real(8) function fde_get_nr_electrons_integrated()
!   ----------------------------------------------------------------------------
      fde_get_nr_electrons_integrated = nr_electrons_integrated
   end function fde_get_nr_electrons_integrated


!   ----------------------------------------------------------------------------
      subroutine loop_over_batches()
!   ----------------------------------------------------------------------------
         real(8) :: dummy
         integer :: idummy
         integer :: ibatch
         integer :: i
         real(8), allocatable :: mpi_buffer(:)

!           print *, "debug, entering loop_over_batches()"
#ifdef VAR_MPI
!           print *, "this is process ",fde_get_my_rank()
#endif

         if (fde_mpi_is_master()) then
            nr_points_batch = fde_grid_im%npoints
         end if   
#ifdef VAR_MPI
         if (parallel_fde) then
!           print *, "nr_ponts_batch #1 is ",nr_points_batch,"for is process ",fde_get_my_rank()
            call fde_mpi_bcast(nr_points_batch)
!           print *, "nr_ponts_batch #2 is ",nr_points_batch,"for is process ",fde_get_my_rank()
         end if
#endif
         ibatch          = 0
         nr_points_total = 0
         nr_points_used  = 0


!        loop over batches
         do
            print *, "looping over batches, ibatch =",ibatch

            if (nr_points_batch < 0) exit

            allocate(rx(nr_points_batch))
            allocate(ry(nr_points_batch))
            allocate(rz(nr_points_batch))
            allocate(rw(nr_points_batch))

            allocate(frozen_n(nr_points_batch))
            allocate(frozen_gnx(nr_points_batch))
            allocate(frozen_gny(nr_points_batch))
            allocate(frozen_gnz(nr_points_batch))
            allocate(frozen_elpot(nr_points_batch))

!radovan: grid stuff will be moved outside this module
            if (fde_mpi_is_master()) then
                rx = fde_grid_im%r(1,:)
                ry = fde_grid_im%r(2,:)
                rz = fde_grid_im%r(3,:)
                rw = fde_grid_im%w
 
                frozen_n   = gf_frozen%n
                frozen_gnx = gf_frozen%gn(1,:)
                frozen_gny = gf_frozen%gn(2,:)
                frozen_gnz = gf_frozen%gn(3,:)
                frozen_elpot= gf_frozen%elpot
!               print *, "values at point 1 for master"
!               print *, "(x,y,z, w) = (",rx(1),",",ry(1),",",rz(1),",",rw(1),")"
!               print *, "(n, gn, v) = (",frozen_n(1),",",frozen_gnx(1),",",frozen_gny(1),",",frozen_gnz(1),",",frozen_elpot(1),")"
            end if

            ibatch = ibatch + 1
            nr_points_total = nr_points_total + nr_points_batch
            nr_points_on_this_proc = nr_points_batch

#ifdef VAR_MPI
!        if (.not.fde_mpi_is_master()) then 
!               print *, "values at point 1 for slave ",fde_get_my_rank(),"before distributing"
!               print *, "(x,y,z, w) = (",rx(1),",",ry(1),",",rz(1),",",rw(1),")"
!               print *, "(n, gn, v) = (",frozen_n(1),",",frozen_gnx(1),",",frozen_gny(1),",",frozen_gnz(1),",",frozen_elpot(1),")"
!        endif
         if (parallel_fde) then
!           redefines nr_points_on_this_proc
!           print *, "distributing points"
            call fde_mpi_distribute_points(rx, ry, rz, rw, nr_points_batch, nr_points_on_this_proc) 
!           print *, "distributing frozen density, potential"
            call fde_mpi_distribute_frozen(frozen_n, frozen_gnx, frozen_gny, frozen_gnz, frozen_elpot,nr_points_batch, nr_points_on_this_proc) 
         end if

!        if (.not.fde_mpi_is_master()) then 
!               print *, "values at point 1 for process ",fde_get_my_rank()
!               print *, "#a",fde_get_my_rank(),rx(1),ry(1),rz(1),rw(1)
!               print *, "#b",fde_get_my_rank(),frozen_n(1),frozen_gnx(1),frozen_gny(1),frozen_gnz(1),frozen_elpot(1)
!        endif
#endif
!        if (do_potential .or. (response_order_mo == 1)) then
!            fde_use_blas3_scheme = .true.
!         else
!            fde_use_blas3_scheme = .false.
!         end if
!         if (fun_is_tau_mgga(xc_fun)) then
!            fde_use_blas3_scheme = .false.
!         end if
!         if (fde_cfg_blocked) then
!            fde_use_blas3_scheme = .true.
!         end if
!         if (fde_cfg_pointwise) then
!            fde_use_blas3_scheme = .false.
!         end if

!           print *, "debug, before integrate()"
#ifdef VAR_MPI
!           print *, "this is process ",fde_get_my_rank()
#endif
            call integrate()
!           print *, "debug, after integrate()"
#ifdef VAR_MPI
!           print *, "this is process ",fde_get_my_rank()
#endif

            deallocate(rx)
            deallocate(ry)
            deallocate(rz)
            deallocate(rw)
            deallocate(frozen_n)
            deallocate(frozen_gnx)
            deallocate(frozen_gny)
            deallocate(frozen_gnz)
            deallocate(frozen_elpot)

!          don't try to read another batch
           goto 999

!     end loop over batches
         end do
 999     continue
         
!           print *, "debug, leaving loop_over_batches()"
#ifdef VAR_MPI
!           print *, "this is process ",fde_get_my_rank()
#endif
      end subroutine loop_over_batches


!   ----------------------------------------------------------------------------
      subroutine lr_fde_uncoupled(block_length,       &
                    mat_dim,            &
                    nz,                 &
                    nr_mat,             &
                    dmaterturbed,       &
                    isym,               &
                    ih,                 &
                    allow_nonhermitian, &
                    ao,                 &
                    n_0,                &
                    n_0_t,              &
                    gn_0_pointwise,     &
                    gn_0_t_pointwise,   &
                    gnn_0,              &
                    gnn_0_t,            &
                    w,                  &
                    is_gga_qr,          &
                    is_gga_qi,          &
                    buffer,             &
                    n,                  &
                    derv_active,        &
                    derv_total,         &
                    u)
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
    real(8), intent(in)    :: n_0_t(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gn_0_t_pointwise(*)
    real(8), intent(in)    :: gnn_0(*)
    real(8), intent(in)    :: gnn_0_t(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga_qr
    logical, intent(in)    :: is_gga_qi
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv_active(max_block_length, 0:n)
    real(8), intent(in)    :: derv_total(max_block_length, 0:n)
!   ----------------------------------------------------------------------------
    logical                :: calculate_qr, calculate_qi
    logical, save          :: have_triplet_contrib(100) = .false.
    logical, save          :: have_triplet = .false.
    logical, save          :: print_triplet_contrib = .true.


    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)
    real(8)                ::  t_b, ts_b(3)

    real(8)                :: z_b
    real(8)                :: z_b_t
    real(8)                :: y_b(3)
    real(8)                :: y_b_t(3)

    integer                :: k_start, k, im, im_off, j
    type(fde_omega_prefactor) :: u
!   ----------------------------------------------------------------------------

    k_start = 1

    do im = 1, nr_mat

      calculate_qr = .false.
      calculate_qi = .false.

      if (ih(im) == +1) calculate_qr = .true.
      if (ih(im) == -1) calculate_qi = .true.

      if (ih(im) == 0 .and. allow_nonhermitian) then
        calculate_qr = .true.
        calculate_qi = .true.
      end if

      if (fde_cfg_no_sdft) then
        calculate_qi = .false.
      end if

      if (calculate_qi) have_triplet_contrib(im) = .true.

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

        u%n(im) = u%n(im) + (derv_total(1, d2000000) - derv_active(1, d2000000))*r_b

        if (is_gga_qr) then
           z_b = 2.0d0*(gn_0_pointwise(1)*gr_b(1) &
                      + gn_0_pointwise(2)*gr_b(2) &
                      + gn_0_pointwise(3)*gr_b(3))

           z_b_t = 2.0d0*(gn_0_t_pointwise(1)*gr_b(1) &
                        + gn_0_t_pointwise(2)*gr_b(2) &
                        + gn_0_t_pointwise(3)*gr_b(3))

           u%n(im) = u%n(im) + (derv_total(1, d1010000)*z_b_t - derv_active(1, d1010000)*z_b)

           u%gn(1:3, im) = u%gn(1:3, im) + 2.0d0*gn_0_t_pointwise(1:3)*r_b*derv_total(1, d1010000)   &
                                         + 2.0d0*gn_0_t_pointwise(1:3)*z_b_t*derv_total(1, d0020000) &
                                         + 2.0d0*gr_b(1:3)*derv_total(1, d0010000)                   &
                                         - 2.0d0*gn_0_pointwise(1:3)*r_b*derv_active(1, d1010000)    &
                                         - 2.0d0*gn_0_pointwise(1:3)*z_b*derv_active(1, d0020000)    &
                                         - 2.0d0*gr_b(1:3)*derv_active(1, d0010000)
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

        do k = k_start, 3
           u%s(k, im) = u%s(k, im) + (derv_total(1, d0200000) - derv_active(1, d0200000))*s_b(k)

           if (is_gga_qi) then
              y_b(k)  = gn_0_pointwise(1)*gs_b(1, k) &
                      + gn_0_pointwise(2)*gs_b(2, k) &
                      + gn_0_pointwise(3)*gs_b(3, k)

              y_b_t(k) = gn_0_t_pointwise(1)*gs_b(1, k) &
                       + gn_0_t_pointwise(2)*gs_b(2, k) &
                       + gn_0_t_pointwise(3)*gs_b(3, k)

              u%s(k, im) = u%s(k, im) + (derv_total(1, d0101000)*y_b_t(k) - derv_active(1, d0101000)*y_b(k))

              u%gs(1:3, k, im) = u%gs(1:3, k, im) +       gn_0_t_pointwise(1:3)*s_b(k)*derv_total(1, d0101000)   &
                                                  +       gn_0_t_pointwise(1:3)*y_b_t(k)*derv_total(1, d0002000) &
                                                  + 2.0d0*gs_b(1:3, k)*derv_total(1, d0000100)                   &
                                                  -       gn_0_pointwise(1:3)*s_b(k)*derv_active(1, d0101000)    &
                                                  -       gn_0_pointwise(1:3)*y_b(k)*derv_active(1, d0002000)    &
                                                  - 2.0d0*gs_b(1:3, k)*derv_active(1, d0000100)
           end if
        end do

      end if !end of qi contribution

    end do

! now for some output
    if (print_triplet_contrib) then
       do im = 1, nr_mat
          if (have_triplet_contrib(im)) have_triplet = .true.
       end do

       if (have_triplet) then
          write(*,*) 'Triplet contributions to FDE kernel included'
       else
          write(*,*) 'Triplet contributions to FDE kernel neglected'
       endif
       print_triplet_contrib = .false.
    endif
!

  end subroutine

end module fde_dirac_matrices_integration
