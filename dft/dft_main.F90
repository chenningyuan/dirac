module dft_main

!   written by T. Helgaker
!   DIRAC adaption by T. Saue

  use memory_allocator
  use dft_cfg
  use dft_functional_derv
  use dft_ao_in_mem
  use dft_ao
  use dft_density
! if you really have to use common blocks, then
! please put them in dft_cb_interface
  use dft_cb_interface
#ifdef MOD_AUTODIFF
  use xc_fun
#endif
#ifdef MOD_AORESPONSE
  use xc_response
#endif
#ifdef VAR_MPI
  use dirac_mpi
#endif

  implicit none

  public dftexc
#ifdef VAR_MPI
  public dftnod
#endif

  private

#include "infpar.h"
#ifdef VAR_MPI
#include "mpif.h"
#ifdef INT_STAR8
  dimension     istat(mpi_status_size)
#else
  integer    :: istat(mpi_status_size)
#endif /* ifdef INT_STAR8 */
#endif /* ifdef VAR_MPI */

  real(8), parameter :: tiny = 1.0d-20

  integer, parameter :: block_length = 1 !later more, now slowdft

  save

  logical, public :: london_g_term   = .false.
  logical, public :: dft_temp_grid_is_done   = .false.
  logical, public :: dft_temp_reload_batches = .false.

  real(8), public :: dft_temp_xc_energy     = 0.0d0
  real(8), public :: dft_temp_xc_mat_energy = 0.0d0
  real(8), public :: dft_temp_nr_electrons  = 0.0d0

  logical :: use_gga_qr = .false. !quaternion real
  logical :: use_gga_qi = .false. !quaternion imaginary
  logical :: use_gga    = .false.

contains

  subroutine dftexc(mat_dim,           &
                    nz,                &
                    nr_fmat,           &
                    fmat,              &
                    dmat_0,            &
                    nr_dmat_p,         &
                    dmat_p,            &
                    isymop_dmat,       &
                    isymop_fmat,       &
                    ihrmop,            &
                    response_order,    &
                    response_order_ao, &
                    g_order,           &
                    m_order)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: mat_dim
    integer, intent(in) :: nz
    integer, intent(in) :: nr_fmat
    real(8)             :: fmat(*)
    real(8), intent(in) :: dmat_0(mat_dim, mat_dim, nz)
    integer, intent(in) :: nr_dmat_p
    real(8), intent(in) :: dmat_p(*)
    integer, intent(in) :: isymop_dmat(*)
    integer, intent(in) :: isymop_fmat(*)
    integer, intent(in) :: ihrmop(*)
    integer, intent(in) :: response_order
    integer, intent(in) :: response_order_ao
    integer, intent(in) :: g_order
    integer, intent(in) :: m_order
!   ----------------------------------------------------------------------------
    integer             :: iz, im, iq
    real(8)             :: time_start, time_in_xc, second
    real(8)             :: ddot

    real(8), allocatable :: dmat_saop_outer(:)
    real(8), allocatable :: dmat_saop_resp(:)

    integer              :: max_ao_g_order, max_ao_m_order, ioff
!   ----------------------------------------------------------------------------

      call init_dft_cb()

!     start timer
      time_start = second()

!     initial expectation value of the two-electron matrix
      if ((response_order == 0) .and. (response_order_ao == 0)) then
        dft_temp_xc_mat_energy = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
      end if

!     get average density from open shells
      if (nr_open_shells > 0) call gtdoav(dmat_0, 1)

!     scale
      call dscal(mat_dim*mat_dim, 2.0d0, dmat_0, 1)


!     scale and insert half-phases
      if (nr_dmat_p > 0) then
        call dscal(mat_dim*mat_dim*nz*nr_dmat_p, 2.0d0, dmat_p, 1)
        if (nz < 4) then
          do iz = 1, nz
            do im = 1, nr_dmat_p
              iq = dft_cb_pq_to_uq(iz, isymop_dmat(im) - 1)
              ioff = mat_dim*mat_dim*nz*(im - 1) + mat_dim*mat_dim*(iz - 1) + 1
              call q2bphase('D', iq, 1, dmat_p(ioff))
            end do
            do im = 1, nr_fmat
              iq = dft_cb_pq_to_uq(iz, isymop_fmat(im) - 1)
              ioff = mat_dim*mat_dim*nz*(im - 1) + mat_dim*mat_dim*(iz - 1) + 1
              call q2bphase('F', iq, 1, fmat(ioff))
            end do
          end do
        end if
      end if

!     calculate grid
      if (.not. dft_temp_grid_is_done) then
        call dftgrd(dmat_0, 1)
        dft_temp_grid_is_done = .true.
        if (dft_cfg_keep_ao .and. ao_mem_allocated) dft_temp_reload_batches = .true.
      end if

!     allocate and get saop density matrices
      if (dft_cfg_saop) then
        dft_cfg_saop_is_active = .true.
        call alloc(dmat_saop_outer, mat_dim*mat_dim)
        if (dft_cfg_saop_with_response_part) then
          call alloc(dmat_saop_resp, mat_dim*mat_dim)
        end if
        call dftac_saop_get_dmats(dmat_saop_outer, &
                                  dmat_saop_resp)
      end if


!     check whether gga code structure will be used

      use_gga_qr = is_gga_functional
      use_gga_qi = is_gga_functional

#ifdef MOD_AUTODIFF
      if (dft_cfg_autodiff) then

!       xc_fun_set_functional is called several times
!       this is because gga query is not ideally structured
        call xc_fun_set_functional(response_order, response_order_ao)

        use_gga_qr = dft_cfg_is_gga_autodiff
        use_gga_qi = dft_cfg_is_gga_autodiff

      end if
#endif

      if (response_order > 0) then
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
      if (response_order > 1) then
        if (dft_cfg_alda_hs .and. (.not. dft_cfg_alda_ha)) then
          call quit('QR and partial ALDA is not implemented')
        end if
        if (dft_cfg_alda_ha .and. (.not. dft_cfg_alda_hs)) then
          call quit('QR and partial ALDA is not implemented')
        end if
      end if

#ifdef VAR_MPI
      if (is_parallel_calculation) then

        call dirac_parctl(4)

        if (nodes > maxnod) then
          call quit('number of nodes is too large, nodes > maxnod')
        end if

        write(print_unit, '(/5x, a, i5/)') &
                          'number of nodes in parallel DFT:', nodes

        timing = .false.
          call sync_dft_cb()
          call sync_dft_cfg()
    
      call dirac_mpi_bcast(mat_dim)
      call dirac_mpi_bcast(nz)
      call dirac_mpi_bcast(nr_fmat)
      call dirac_mpi_bcast(nr_dmat_p)

      call sync_dmat_0(mat_dim, nz, dmat_0)

      if (dft_cfg_saop_is_active) then
        call sync_dmat_0(mat_dim, nz, dmat_saop_outer)
        if (dft_cfg_saop_with_response_part) then
          call sync_dmat_0(mat_dim, nz, dmat_saop_resp)
        end if
      end if

      if (nr_dmat_p > 0) then
        call sync_dmat_p(mat_dim, nz, nr_dmat_p, dmat_p, isymop_dmat, nr_fmat, isymop_fmat, ihrmop)
      end if

      call dirac_mpi_bcast(response_order)
      call dirac_mpi_bcast(response_order_ao)
      call dirac_mpi_bcast(g_order)
      call dirac_mpi_bcast(m_order)

    end if
#endif /* ifdef VAR_MPI */

!   at this stage all the variables of the next call
!   should be synced between nodes      
!   and all arrays should be allocated at each node
!   and synced between nodes
!   but do not broadcast any other variables that can be synced
!   more compactly deeper down in dftdrv
    call dftdrv(mat_dim,           &
                nz,                &
                nr_fmat,           &
                fmat,              &
                dmat_0,            &
                dmat_saop_outer,   &
                dmat_saop_resp,    &
                nr_dmat_p,         &
                dmat_p,            &
                isymop_dmat,            &
                isymop_fmat,            &
                ihrmop,            &
                response_order,    &
                response_order_ao, &
                g_order,           &
                m_order)

!     deallocate saop density matrices
      if (dft_cfg_saop) then
        call dealloc(dmat_saop_outer)
        if (dft_cfg_saop_with_response_part) then
          call dealloc(dmat_saop_resp)
        end if
      end if

!     rescale
      call dscal(mat_dim*mat_dim, 0.5d0, dmat_0, 1)

!     get average density from open shells
      if (nr_open_shells > 0) call gtdoav(dmat_0, -1)

!     in case of gga symmetrize matrices
      if (use_gga) then
        if (nr_dmat_p > 0) then
          call gga_sym(mat_dim, nz, nr_fmat, fmat, isymop_fmat, ihrmop)
        else
          call gga_sym(mat_dim, nz, 1, fmat, (/1/), (/1/))
        end if
      end if

!     scale and insert half-phases
      if (nr_dmat_p > 0) then
        call dscal(mat_dim*mat_dim*nz*nr_dmat_p, 0.5d0, dmat_p, 1)
        if (nz < 4) then
          do iz = 1, nz
            do im = 1, nr_dmat_p
              iq = dft_cb_pq_to_uq(iz, isymop_dmat(im) - 1)
              ioff = mat_dim*mat_dim*nz*(im - 1) + mat_dim*mat_dim*(iz - 1) + 1
              call q2bphase('D', iq, 1, dmat_p(ioff))
            end do
            do im = 1, nr_fmat
              iq = dft_cb_pq_to_uq(iz, isymop_fmat(im) - 1)
              ioff = mat_dim*mat_dim*nz*(im - 1) + mat_dim*mat_dim*(iz - 1) + 1
              call q2bphase('F', iq, 1, fmat(ioff))
            end do
          end do
        end if
      end if

!     calculate final expectation value of two-electron matrix
      if ((response_order == 0) .and. (response_order_ao == 0)) then
!       the xc contribution can then be subtracted 
        dft_temp_xc_mat_energy = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1) &
                               - dft_temp_xc_mat_energy
      end if 

!     timing
      time_in_xc = second() - time_start
      write(print_unit, '(1x)')
      call timtxt('>>>> time used in dftexc:', time_in_xc, print_unit)
      write(print_unit, '(1x)')

      end subroutine


  subroutine dftdrv(mat_dim,           &
                    nz,                &
                    nr_fmat,           &
                    fmat,              &
                    dmat_0,            &
                    dmat_saop_outer,   &
                    dmat_saop_resp,    &
                    nr_dmat_p,         &
                    dmat_p,            &
                    isymop_dmat,       &
                    isymop_fmat,       &
                    ihrmop,            &
                    response_order,    &
                    response_order_ao, &
                    g_order,           &
                    m_order)

    use fde_mod
!   ----------------------------------------------------------------------------
    integer, intent(in) :: mat_dim
    integer, intent(in) :: nz
    integer, intent(in) :: nr_fmat
    real(8)             :: fmat(*)
    real(8), intent(in) :: dmat_0(*)
    real(8), intent(in) :: dmat_saop_outer(*)
    real(8), intent(in) :: dmat_saop_resp(*)
    integer, intent(in) :: nr_dmat_p
    real(8), intent(in) :: dmat_p(*)
    integer, intent(in) :: isymop_dmat(*)
    integer, intent(in) :: isymop_fmat(*)
    integer, intent(in) :: ihrmop(*)
    integer, intent(in) :: response_order
    integer, intent(in) :: response_order_ao
    integer, intent(in) :: g_order
    integer, intent(in) :: m_order
!   ----------------------------------------------------------------------------
    integer             :: ibatch







      real(8), allocatable :: kx(:)
      real(8), allocatable :: ky(:)
      real(8), allocatable :: kz(:)
      real(8), allocatable :: kw(:)

      integer, allocatable :: idisp(:)
      integer, allocatable :: isend(:)
      real(8), allocatable :: kxdum(:)
      real(8) :: dummy
      real(8) :: error_nr_electrons, nr_electrons_orbital_total
      integer :: i, idummy, ierr, myNPOINT, nr_points_part, NODESALL
      integer :: npnts, npoint, nrest, IDISP_REST
      logical :: i_am_the_master


    dft_temp_nr_electrons = 0.0d0
    dft_temp_xc_energy    = 0.0d0
    i_am_the_master = .false.
    if (mytid == mparid) i_am_the_master = .true.

    if (i_am_the_master) then
      rewind dft_temp_file_unit
      if (dft_cfg_multigrid) then
        rewind dft_cfg_multigrid_file_unit(response_order_ao)
      end if
      npnts  = 0
      ibatch = 0
    end if

    if (dft_cfg_keep_ao .and. ao_mem_allocated .and. dft_temp_reload_batches) then
      call kill_ao_in_mem()
      ao_mem_allocated        = .false.
      dft_temp_reload_batches = .false.
    end if


! aspg, toying with having the frozen density and adfgrid available
! for the respose calcs or otherwise...
    if (have_frozen_density) then
       call FDE_testFrozenDensity()
    endif

!*    Loop over batches of grid points
  100 CONTINUE

      if (i_am_the_master) then
        if (dft_cfg_import_grid) then
          read(dft_temp_file_unit, *) npoint
        elseif (dft_cfg_multigrid) then
          read(dft_cfg_multigrid_file_unit(response_order_ao), *) npoint
        else
          read(dft_temp_file_unit) npoint
        end if
      end if

#ifdef VAR_MPI
      if (is_parallel_calculation) then
        call dirac_mpi_bcast(npoint)
      end if
#endif

      IF( NPOINT .lt. 0 ) THEN
         GOTO 200
      ELSE IF(  NPOINT == 0 ) THEN
          GOTO 100      
      ELSE
        IF(MYTID.EQ.MPARID) THEN
          npnts  = npnts  + npoint
          ibatch = ibatch + 1

          call alloc(kx, npoint)
          call alloc(ky, npoint)
          call alloc(kz, npoint)
          call alloc(kw, npoint)

          if (dft_cfg_import_grid) then
            call import_dft_grid(kx, ky, kz, kw, dft_temp_file_unit, npoint)
          elseif (dft_cfg_multigrid) then
            call import_dft_grid(kx, ky, kz, kw, dft_cfg_multigrid_file_unit(response_order_ao), npoint)
          else
            call reaqua(kx, ky, kz, kw, dft_temp_file_unit, npoint)
          end if
        end if

#ifdef VAR_MPI
        if (is_parallel_calculation) then
          call dirac_mpi_bcast(ibatch)

          nodesall = nodes + 1

!         if aos are in mem one could skip the broadcast of
!         x, y, and z (still need the w though)
!         but for the moment let's do it anyway

          call alloc(idisp, nodesall)
          call alloc(isend, nodesall)

          nr_points_part = NPOINT/NODESALL
          NREST = MOD(NPOINT,NODESALL)
          myNPOINT = nr_points_part
!.........send nothing to master
          IF(MYTID.EQ.MPARID) THEN ! master
            ISEND(1) = 0
            IDISP(1) = 0
            DO I = 1,NODES-1
              ISEND(I+1) = nr_points_part
              IDISP(I+1) = nr_points_part*I
            ENDDO
            ISEND(NODESALL) = nr_points_part + NREST
            IDISP(NODESALL) = nr_points_part*NODES
!           send data to slaves - all grid points in one send chunk
!           FIXME: scatterv makes problems with INT_STAR8...
!                  SK - Apr 2009
!
#if !defined (INT_STAR8)
!............scatterv data
             CALL MPI_SCATTERV(KX,ISEND,IDISP,MPI_REAL8, &
     &                         DUMMY,ISEND,MPI_REAL8,MPARID, &
     &                         MPI_COMM_WORLD,IERR)
             CALL MPI_SCATTERV(KY,ISEND,IDISP,MPI_REAL8, &
     &                         DUMMY,ISEND,MPI_REAL8,MPARID, &
     &                         MPI_COMM_WORLD,IERR)
             CALL MPI_SCATTERV(KZ,ISEND,IDISP,MPI_REAL8, &
     &                         DUMMY,ISEND,MPI_REAL8,MPARID, &
     &                         MPI_COMM_WORLD,IERR)
             CALL MPI_SCATTERV(KW,ISEND,IDISP,MPI_REAL8, &
     &                         DUMMY,ISEND,MPI_REAL8,MPARID, &
     &                         MPI_COMM_WORLD,IERR)
#else
!............scatter data and send remainings to last node
             call alloc(kxdum, mynpoint)
             CALL MPI_SCATTER(KX,myNPOINT,MPI_REAL8, &
     &                        KXDUM,myNPOINT,MPI_REAL8,MPARID, &
     &                        MPI_COMM_WORLD,IERR)
             CALL MPI_SCATTER(KY,myNPOINT,MPI_REAL8, &
     &                        KXDUM,myNPOINT,MPI_REAL8,MPARID, &
     &                        MPI_COMM_WORLD,IERR)
             CALL MPI_SCATTER(KZ,myNPOINT,MPI_REAL8, &
     &                        KXDUM,myNPOINT,MPI_REAL8,MPARID, &
     &                        MPI_COMM_WORLD,IERR)
             CALL MPI_SCATTER(KW,myNPOINT,MPI_REAL8, &
     &                        KXDUM,myNPOINT,MPI_REAL8,MPARID, &
     &                        MPI_COMM_WORLD,IERR)
             call dealloc(kxdum)
             IF( NREST .ne. 0 )THEN
               IDISP_REST = IDISP(NODESALL) + nr_points_part
               CALL MPI_SEND(KX(1+IDISP_REST),NREST,MPI_REAL8, &
     &                       NODES,90,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(KY(1+IDISP_REST),NREST,MPI_REAL8, &
     &                       NODES,91,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(KZ(1+IDISP_REST),NREST,MPI_REAL8, &
     &                       NODES,92,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(KW(1+IDISP_REST),NREST,MPI_REAL8, &
     &                       NODES,93,MPI_COMM_WORLD,IERR)
             END IF
#endif
          ELSE ! slaves
!
!           the last process gets the rest; these points are usually
!           the most diffuse ones and therefore the most screened ones
            IF(MYTID.eq.NODES) myNPOINT = nr_points_part+NREST
            call alloc(kx, mynpoint)
            call alloc(ky, mynpoint)
            call alloc(kz, mynpoint)
            call alloc(kw, mynpoint)
!           receive data from master - all grid points in one receive chunk
!
!           FIXME: scatterv does makes problems with INT_STAR8...
!                  SK - Apr 2009
#if !defined (INT_STAR8)
!...........receive data via scatterv
            CALL MPI_SCATTERV(DUMMY,IDUMMY,IDUMMY,MPI_REAL8, &
     &                        KX,myNPOINT,MPI_REAL8,MPARID, &
     &                        MPI_COMM_WORLD,IERR)
            CALL MPI_SCATTERV(DUMMY,IDUMMY,IDUMMY,MPI_REAL8, &
     &                        KY,myNPOINT,MPI_REAL8,MPARID, &
     &                        MPI_COMM_WORLD,IERR)
            CALL MPI_SCATTERV(DUMMY,IDUMMY,IDUMMY,MPI_REAL8, &
     &                        KZ,myNPOINT,MPI_REAL8,MPARID, &
     &                        MPI_COMM_WORLD,IERR)
            CALL MPI_SCATTERV(DUMMY,IDUMMY,IDUMMY,MPI_REAL8, &
     &                        KW,myNPOINT,MPI_REAL8,MPARID, &
     &                        MPI_COMM_WORLD,IERR)
#else
!...........reset for scatter length
            IF(MYTID.eq.NODES) myNPOINT = nr_points_part
!...........receive data via scatter
            CALL MPI_SCATTER(DUMMY,IDUMMY,MPI_REAL8, &
     &                       KX,myNPOINT,MPI_REAL8,MPARID, &
     &                       MPI_COMM_WORLD,IERR)
            CALL MPI_SCATTER(DUMMY,IDUMMY,MPI_REAL8, &
     &                       KY,myNPOINT,MPI_REAL8,MPARID, &
     &                       MPI_COMM_WORLD,IERR)
            CALL MPI_SCATTER(DUMMY,IDUMMY,MPI_REAL8, &
     &                       KZ,myNPOINT,MPI_REAL8,MPARID, &
     &                       MPI_COMM_WORLD,IERR)
            CALL MPI_SCATTER(DUMMY,IDUMMY,MPI_REAL8, &
     &                       KW,myNPOINT,MPI_REAL8,MPARID, &
     &                       MPI_COMM_WORLD,IERR)
!...........remainings for the last node
            IF(MYTID .eq. NODES .and. NREST .ne. 0 ) THEN
              IDISP_REST = nr_points_part
              CALL MPI_RECV(KX(1+IDISP_REST),NREST,MPI_REAL8, &
     &                      MPARID,90,MPI_COMM_WORLD,ISTAT,IERR)
              CALL MPI_RECV(KY(1+IDISP_REST),NREST,MPI_REAL8, &
     &                      MPARID,91,MPI_COMM_WORLD,ISTAT,IERR)
              CALL MPI_RECV(KZ(1+IDISP_REST),NREST,MPI_REAL8, &
     &                      MPARID,92,MPI_COMM_WORLD,ISTAT,IERR)
              CALL MPI_RECV(KW(1+IDISP_REST),NREST,MPI_REAL8, &
     &                      MPARID,93,MPI_COMM_WORLD,ISTAT,IERR)
            END IF
#endif
          end if !IF(MYTID.EQ.MPARID) THEN ! master
!.........add rest to total number of myNPOINT on last node
          IF(MYTID.eq.NODES) myNPOINT = nr_points_part+NREST
          call dealloc(idisp)
          call dealloc(isend)
        ELSE ! sequential calculation
#endif
            myNPOINT = NPOINT
#if defined (VAR_MPI) 
        end if !IF(is_parallel_calculation) THEN
#endif
        
        call dftdr1(mat_dim,           &
                    nz,                &
                    nr_fmat,           &
                    fmat,              &
                    dmat_0,            &
                    dmat_saop_outer,   &
                    dmat_saop_resp,    &
                    nr_dmat_p,         &
                    dmat_p,            &
                    isymop_dmat,       &
                    isymop_fmat,       &
                    ihrmop,            &
                    response_order,    &
                    response_order_ao, &
                    g_order,           &
                    m_order,           &
                    kx,                &
                    ky,                &
                    kz,                &
                    kw,                &
                    mynpoint,          &
                    ibatch)

        call dealloc(kx)
        call dealloc(ky)
        call dealloc(kz)
        call dealloc(kw)

        goto 100
        end if !IF( NPOINT .lt. 0 ) THEN
 200    CONTINUE



#if defined (VAR_MPI)
    if (is_parallel_calculation) then
      call dft_collect(mat_dim,        &
                       nz,             &
                       nr_fmat,        &
                       fmat,           &
                       response_order, &
                       response_order_ao)
    end if
#endif


      IF(MYTID.NE.MPARID) RETURN

      nr_electrons_orbital_total = nr_electrons_orbital &
                                 + nr_electrons_orbital_active
      error_nr_electrons  = dft_temp_nr_electrons &
                          - nr_electrons_orbital_total

      write(print_unit, *)
      write(print_unit, '(3x, a, i24)')                          &
              'number of quadrature points:  ',                  &
              npnts
      if ((response_order == 0) .and. (response_order_ao == 0)) then
        write(print_unit, '(3x, a, f24.16)')                     &
              'DFT exchange-correlation energy:               ', &
              dft_temp_xc_energy
      end if
      write(print_unit, '(3x, a, f24.16)')                       &
              'number of electrons from numerical integration:', &
              dft_temp_nr_electrons
      if (abs(error_nr_electrons) > dft_cfg_accepted_error) then
        write(print_unit, '(3x, a, f14.6)')                      &
              'number of electrons from orbital occupations:  ', &
               nr_electrons_orbital_total
        write(print_unit, '(3x, a, f14.6)')                      &
              'WARNING: error in the number of electrons:     ', &
               error_nr_electrons
        write(print_unit, '(12x, a, f14.6)')                     &
              'is larger than .ERRELS (set input):   ',          &
               dft_cfg_accepted_error
        write(print_unit, '(12x, a)') &
              'this can happen when starting from coefficients'
        write(print_unit, '(12x, a)') &
              'from a different geometry'
        write(print_unit, '(12x, a)')                 &
              'or it can mean that the quadrature '// &
              'grid is inappropriate'
      end if
      write(print_unit, *)

      end subroutine

    subroutine dftdr1(mat_dim,           &
                      nz,                &
                      nr_fmat,           &
                      fmat,              &
                      dmat_0,            &
                      dmat_saop_outer,   &
                      dmat_saop_resp,    &
                      nr_dmat_p,         &
                      dmat_p,            &
                      isymop_dmat,       &
                      isymop_fmat,       &
                      ihrmop,            &
                      response_order,    &
                      response_order_ao, &
                      g_order,           &
                      m_order,           &
                      corx,              &
                      cory,              &
                      corz,              &
                      weight,            &
                      npoint,            &
                      ibatch)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: mat_dim
    integer, intent(in) :: nz
    integer, intent(in) :: nr_fmat
    real(8)             :: fmat(*)
    real(8), intent(in) :: dmat_0(*)
    real(8), intent(in) :: dmat_saop_outer(*)
    real(8), intent(in) :: dmat_saop_resp(*)
    integer, intent(in) :: nr_dmat_p
    real(8), intent(in) :: dmat_p(*)
    integer, intent(in) :: isymop_dmat(*)
    integer, intent(in) :: isymop_fmat(*)
    integer, intent(in) :: ihrmop(*)
    integer, intent(in) :: response_order
    integer, intent(in) :: response_order_ao
    integer, intent(in) :: g_order
    integer, intent(in) :: m_order
    real(8), intent(in) :: corx(*)
    real(8), intent(in) :: cory(*)
    real(8), intent(in) :: corz(*)
    real(8), intent(in) :: weight(*)
    integer, intent(in) :: npoint
    integer, intent(in) :: ibatch
!   ----------------------------------------------------------------------------
    integer             :: ioff, i


      real(8) :: RHG(3),VT(3)
      real(8) :: derv(mxbllen, 0:99)
      real(8) :: r_0(mxbllen)
      real(8) :: z_0(mxbllen)
      real(8) :: w(mxbllen)
      real(8) :: rho13, rho_outer, rho_resp

      integer :: ipoint, need_limit
      real(8) :: ddot, f

      integer :: itest, jtest
    integer             :: max_ao_g_order
    integer             :: max_ao_m_order

!   this pointer is useful for .KEEPAO, then it can point to memory
!   otherwise it points to ao
    real(8), pointer             :: ao_p(:)

    real(8), allocatable, target :: ao(:)
    real(8), allocatable         :: buffer(:)


#if defined (VAR_MPI)
    if (is_parallel_calculation) then
      call dirac_mpi_bcast(use_gga_qr)
      call dirac_mpi_bcast(use_gga_qi)
      call dirac_mpi_bcast(use_gga)
      call dftsyncfunc(mytid == mparid)
    end if
#endif

    max_ao_g_order = 0
    if (use_gga) then
      max_ao_g_order = 1
    end if

    max_ao_m_order = m_order

    call dft_ao_init(is_parallel_calculation, &
                     max_ao_g_order,          &
                     max_ao_m_order,          &
                     dft_cfg_screening,       &
                     dft_cfg_dfthri,          &
                     .false.)

    call dft_density_init(is_parallel_calculation)

#ifdef MOD_AUTODIFF
    if (dft_cfg_autodiff) then
      call xc_fun_set_functional(response_order, response_order_ao)
    end if
#endif

    call alloc(ao,     total_ao_vector_length)
    call alloc(buffer, total_ao_vector_length)

    if (dft_cfg_keep_ao) then
      if (dft_cfg_ao_mem_limit < 1) then
!       write(*, *) 'debug will not keep ao in mem, mem too small'
        dft_cfg_keep_ao = .false.
      end if
!     write(*, *) 'debug will try to keep ao in mem', dft_cfg_ao_mem_limit
!   else
!     write(*, *) 'debug will not try to keep ao in mem', dft_cfg_ao_mem_limit
    end if

    if (dft_cfg_keep_ao) then
      if (.not. ao_mem_allocated) then
!       allocate
        write(*, *) 'debug allocating ao in mem'
        call init_ao_in_mem(dft_cfg_ao_mem_limit)
        ao_mem_allocated = .true.
!       reset
        batch_pointer    = 1
        batch_is_in_mem  = .false.
      end if

      ioff = batch_pointer(ibatch)
      if (.not. batch_is_in_mem(ibatch)) then
        do i = ibatch + 1, max_nr_of_batches
          batch_pointer(i) = batch_pointer(i) + npoint*total_ao_vector_length
        end do
      end if

      need_limit = ioff + npoint*total_ao_vector_length
      if (need_limit > dft_cfg_ao_mem_limit) then
        batch_fits_in_mem(ibatch) = .false.
        write(*, *) 'debug cannot keep batch', ibatch
        write(*, *) 'debug need', need_limit*8/1000000
        write(*, *) 'debug have', dft_cfg_ao_mem_limit*8/1000000
      end if
    end if

    if (.not. (dft_cfg_keep_ao .and. batch_is_in_mem(ibatch))) then
!     so we either don't use .KEEPAO or the AOs are not yet stored
!     so let it point to ao
      ao_p => ao
    end if    

    do ipoint = 1, npoint

      if (ipoint == 1) then
        if (dft_cfg_keep_ao .and. batch_is_in_mem(ibatch)) then
          write(*, *) 'debug reading ao from mem', ioff
!       else
!         write(*, *) 'debug recalculating aos'
        end if
      end if
      
      if (.not. (dft_cfg_keep_ao .and. batch_is_in_mem(ibatch))) then
        call get_so(corx(ipoint), cory(ipoint), corz(ipoint), &
                    ao_p, buffer)
      end if
      if (dft_cfg_keep_ao) then

!       copy ao to mem
        if (batch_fits_in_mem(ibatch) .and. .not. batch_is_in_mem(ibatch)) then
          if (ipoint == 1) then
            write(*, *) 'debug storing ao to mem', ioff
          end if
          call dcopy(total_ao_vector_length, ao_p, 1, ao_in_mem(ioff), 1)
        end if

!       if already in mem, set pointer
        if (batch_is_in_mem(ibatch)) then
          ao_p => ao_in_mem(ioff:ioff + total_ao_vector_length)
        end if
        ioff = ioff + total_ao_vector_length
        
        if (ipoint == npoint) then
          if (batch_fits_in_mem(ibatch)) batch_is_in_mem(ibatch) = .true.
        end if
      end if

!     density
      call get_n(r_0(1), 0, mat_dim, dmat_0, buffer, ao_p)

      if (r_0(1) > dft_cfg_dfthr0) then

!       set it to zero, otherwise it is undefined with lda
        z_0 = 0.0d0

!       density gradient
        if (need_ao_order(1, 0)) then
          call get_gn(rhg, 0, mat_dim, dmat_0, buffer, ao_p, ao_p(ao_start(1, 0):ao_end(1, 0)))
          z_0(1) = ddot(3, rhg, 1, rhg, 1)
        end if

        w(1) = weight(ipoint)

!       number of electrons
        dft_temp_nr_electrons = dft_temp_nr_electrons + w(1)*r_0(1)


!       exchange-correlation potential 
!       ==============================

        if ((response_order == 0) .and. (response_order_ao == 0)) then

          call get_functional_derv(1,                  &
                                   block_length,       &
                                   w,                  &
                                   r_0,                &
                                   z_0,                &
                                   nr_nonzero_derv(1), &
                                   derv)

!         energy density
          dft_temp_xc_energy = dft_temp_xc_energy + derv(1, d_0)

          if (use_gga) then
!           this factor has nothing to do
!           with symmetrization later
!           it comes from the fact that derivative of z
!           yields two identical terms (one left, one right)
            derv(1, d_z) = 2.0d0*derv(1, d_z)
          end if

!         asymptotic correction
          if (dft_cfg_saop_is_active) then
            call get_n(rho_outer,       &
                       0,               &
                       mat_dim,         &
                       dmat_saop_outer, &
                       buffer,          &
                       ao_p)
            if (dft_cfg_saop_with_response_part) then
              call get_n(rho_resp,       &
                         0,              &
                         mat_dim,        &
                         dmat_saop_resp, &
                         buffer,         &
                         ao_p)
            end if
          end if
          if (dft_cfg_grac_is_active .or. dft_cfg_saop_is_active) then
            rho13 = r_0(1)**(1.0d0/3.0d0)
            call dftac(derv, w(1), r_0(1), rho13, &
     &                 rho_resp, rho_outer, z_0(1))
          end if

        end if !exchange-correlation potential




!       exchange-correlation contribution to kohn-sham matrix
        if ((response_order == 0) .and. (response_order_ao == 0)) then
!       here we will add LAO contribution to KS matrix

! if london: 
          if (max_ao_m_order == 1) then
            call lda_london_distribute_r(derv(1, d_r),                      &
                                         0,                                 &
                                         mat_dim,                           &
                                         fmat,                              &
                                         ao_p,                              &
                                         ao_p(ao_start(0, 1):ao_end(0, 1)), &
                                         corx(ipoint), cory(ipoint), corz(ipoint))

          else if (max_ao_m_order == 0) then
! if not london:

          if (use_gga) then
            vt = derv(1, d_z)*rhg

#ifdef TEST_GAO1_SYMMETRY
            if (abs(corx(ipoint)) > 1.0d-10) then
              vt(2) = vt(1)*cory(ipoint)/corx(ipoint)
              vt(3) = vt(1)*corz(ipoint)/corx(ipoint)
            end if
#endif

            call gga_distribute_r(derv(1, d_r), &
                                  vt,           &
                                  0,            &
                                  mat_dim,      &
                                  fmat,         &
                                  ao_p,         &
                                  ao_p(ao_start(1, 0):ao_end(1, 0)))
          else
            call lda_distribute_r(derv(1, d_r), &
                                  0,            &
                                  mat_dim,      &
                                  fmat,         &
                                  ao_p)
          end if
        end if ! max_ao_m_order==1
        end if


!       linear response
        if (response_order == 1) then

          if (max_ao_m_order == 1) then
! gosia: temporary fix for london
            call g_lao_contrib(block_length,    &
                               mat_dim,          &
                               nz,          &
                               nr_dmat_p,          &
                               dmat_p,          &
                               fmat,          &
                               isymop_dmat,          &
                               ihrmop,          &
                               .false.,         &
                               ao_p,            &
                               r_0,             &
                               rhg,             &
                               z_0,             &
                               w,               &
                               use_gga_qr, &
                               use_gga_qi, &
                               buffer)
!no london:
          else if (max_ao_m_order == 0) then

          call lr(block_length,    &
                  mat_dim,          &
                  nz,          &
                  nr_dmat_p,          &
                  dmat_p,          &
                  fmat,          &
                  isymop_dmat,     &
                  ihrmop,          &
                  .false.,         &
                  ao_p,            &
                  r_0,             &
                  rhg,             &
                  z_0,             &
                  w,               &
                  use_gga_qr, &
                  use_gga_qi, &
                  buffer)

        endif
        endif


!       quadratic response
        if (response_order == 2) then

          call sdft_qr(block_length,    &
                       use_gga,          &
                       dft_cfg_no_sdft, &
                       mat_dim, dmat_p,fmat,   &
                       isymop_dmat,ihrmop,ao_p,ao_p(ao_start(1, 0)), &
                       r_0(1),z_0(1),rhg,w(1),buffer)

!         the idea here is that   r^[B]               \omega_{ia}^[C]
!                                 r^[C]               \omega_{ia}^[B]
!                                 r^[BC]              \omega_{ia}
!         and the corresponding   q^[B]  \cdot \nabla \omega_{ia}^[C]
!                                 q^[C]  \cdot \nabla \omega_{ia}^[B]
!                                 q^[BC] \cdot \nabla \omega_{ia}
          
!         can be evaluated using lr code (below)
!         only the remaining part is done in sdft_qr (above)
!         
          call lr(block_length,    &
                  mat_dim,          &
                  nz,          &
                  nr_dmat_p,          &
                  dmat_p,          &
                  fmat,          &
                  isymop_dmat,          &
                  ihrmop,          &
                  .true.,          &
                  ao_p,            &
                  r_0,             &
                  rhg,             &
                  z_0,             &
                  w,               &
                  use_gga_qr, &
                  use_gga_qi, &
                  buffer)

        endif


#if defined MOD_AUTODIFF && MOD_AORESPONSE
          if (response_order_ao > 0) then
            call xc_response_contribution(response_order_ao,                 &
                                          mat_dim,   &
                                          nz,        &
                                          fmat,      &
                                          nr_dmat_p, &
                                          dmat_p,    &
                                          isymop_dmat,      &
                                          isymop_fmat,      &
                                          use_gga,                           &
                                          ao_p(ao_start(0, 0):ao_end(0, 0)), &
                                          ao_p(ao_start(1, 0):ao_end(1, 0)), &
                                          r_0(1),                            &
                                          rhg,                               &
                                          corx(ipoint),                      &
                                          cory(ipoint),                      &
                                          corz(ipoint),                      &
                                          w(1),                              &
                                          buffer)
          end if
#endif


      end if !if (r_0(1) > dft_cfg_dfthr0) then
    end do !do ipoint = 1, npoint

    nullify(ao_p)

    call dealloc(ao)
    call dealloc(buffer)
      

  end subroutine

#ifdef VAR_MPI
  subroutine dft_collect(mat_dim,        &
                         nz,             &
                         nr_fmat,        &
                         fmat,           &
                         response_order, &
                         response_order_ao)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: mat_dim
    integer, intent(in)  :: nz
    integer, intent(in)  :: nr_fmat
    real(8)              :: fmat(*)
    integer, intent(in)  :: response_order
    integer, intent(in)  :: response_order_ao
!   ----------------------------------------------------------------------------
    integer              :: ierr
    integer              :: length_fmat
    real(8)              :: sum_nr_electrons
    real(8)              :: sum_xc_energy
    real(8), allocatable :: sum_fmat(:)
!   ----------------------------------------------------------------------------

    call mpi_reduce(dft_temp_nr_electrons, &
                    sum_nr_electrons,      &
                    1,                     &
                    mpi_double_precision,  &
                    mpi_sum,               &
                    0,                     &
                    mpi_comm_world,        &
                    ierr) 

    if ((response_order == 0) .and. (response_order_ao == 0)) then
      call mpi_reduce(dft_temp_xc_energy,   &
                      sum_xc_energy,        &
                      1,                    &
                      mpi_double_precision, &
                      mpi_sum,              &
                      0,                    &
                      mpi_comm_world,       &
                      ierr) 
    end if

    length_fmat = mat_dim*mat_dim*nz*nr_fmat

    if (mytid == mparid) then
      call alloc(sum_fmat, length_fmat)
      sum_fmat = 0.0d0
    end if

    call mpi_reduce(fmat,                 &
                    sum_fmat,             &
                    length_fmat,          &
                    mpi_double_precision, &
                    mpi_sum,              &
                    0,                    &
                    mpi_comm_world,       &
                    ierr) 

    if (mytid == mparid) then
      dft_temp_nr_electrons = sum_nr_electrons
      if ((response_order == 0) .and. (response_order_ao == 0)) then
        dft_temp_xc_energy = sum_xc_energy
      end if
      call dcopy(length_fmat, sum_fmat, 1, fmat, 1)     
      call dealloc(sum_fmat)
    end if

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DFTNOD()
!***********************************************************************
!     Based on HERNOD 
!     written by T. Saue March 4
!
!    This is the node program for DFT                           *
!***********************************************************************



      integer, allocatable :: isymop_dmat(:)
      integer, allocatable :: isymop_fmat(:)
      integer, allocatable :: ihrmop(:)
      real(8), allocatable :: fmat(:)
      real(8), allocatable :: dmat_0(:)
      real(8), allocatable :: dmat_p(:)
      real(8), allocatable :: dmat_saop_outer(:)
      real(8), allocatable :: dmat_saop_resp(:)
      integer              :: nr_dmat_p, nr_fmat
      real(8)              :: wall1, cpu1
      integer              :: iprint
      integer              :: max_ao_g_order, max_ao_m_order
      integer              :: mat_dim, nz
      integer              :: response_order
      integer              :: response_order_ao
      integer              :: g_order
      integer              :: m_order
!
!
!
!     Start timing.
!
      CALL GETTIM(CPU1,WALL1)
!
!     This is a slave in a parallel run
!
      parher = .true.
      slave  = .true.
      debug  = .false.
!
!     Receive initialization data from master 1
!     =========================================
      call sync_dft_cb()
      call sync_dft_cfg()

    call dirac_mpi_bcast(mat_dim)
    call dirac_mpi_bcast(nz)
    call dirac_mpi_bcast(nr_fmat)
    call dirac_mpi_bcast(nr_dmat_p)
       
    call alloc(dmat_0, mat_dim*mat_dim*nz)
    call sync_dmat_0(mat_dim, nz, dmat_0)

    if (dft_cfg_saop_is_active) then
      call alloc(dmat_saop_outer, mat_dim*mat_dim*nz)
      call sync_dmat_0(mat_dim, nz, dmat_saop_outer)
      if (dft_cfg_saop_with_response_part) then
        call alloc(dmat_saop_resp, mat_dim*mat_dim*nz)
        call sync_dmat_0(mat_dim, nz, dmat_saop_resp)
      end if
    end if

    if (nr_dmat_p > 0) then
      call alloc(dmat_p, mat_dim*mat_dim*nz*nr_dmat_p)
      call alloc(isymop_dmat, nr_dmat_p)
      call alloc(isymop_fmat, nr_fmat)
      call alloc(ihrmop, nr_dmat_p)
      call sync_dmat_p(mat_dim, nz, nr_dmat_p, dmat_p, isymop_dmat, nr_fmat, isymop_fmat, ihrmop)
    end if

    call dirac_mpi_bcast(response_order)
    call dirac_mpi_bcast(response_order_ao)
    call dirac_mpi_bcast(g_order)
    call dirac_mpi_bcast(m_order)

    call alloc(fmat, mat_dim*mat_dim*nz*nr_fmat)
    fmat = 0.0d0

    call dftdrv(mat_dim,           &
                nz,                &
                nr_fmat,           &
                fmat,              &
                dmat_0,            &
                dmat_saop_outer,   &
                dmat_saop_resp,    &
                nr_dmat_p,         &
                dmat_p,            &
                isymop_dmat,       &
                isymop_fmat,       &
                ihrmop,            &
                response_order,    &
                response_order_ao, &
                g_order,           &
                m_order)

    call dealloc(fmat)
    call dealloc(dmat_0)

    if (dft_cfg_saop_is_active) then
      call dealloc(dmat_saop_outer)
      if (dft_cfg_saop_with_response_part) then
        call dealloc(dmat_saop_resp)
      end if
    end if

    if (nr_dmat_p > 0) then
      call dealloc(dmat_p)
      call dealloc(isymop_dmat)
      call dealloc(isymop_fmat)
      call dealloc(ihrmop)
    end if

  end subroutine

  subroutine sync_dmat_0(mat_dim, &
                         nz,      &
                         dmat_0)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: mat_dim
    integer, intent(in) :: nz
    real(8)             :: dmat_0(mat_dim, mat_dim, nz)
!   ----------------------------------------------------------------------------

    call dirac_mpi_bcast(dmat_0)

  end subroutine

  subroutine sync_dmat_p(mat_dim,     &
                         nz,          &
                         nr_dmat_p,   &
                         dmat_p,      &
                         isymop_dmat, &
                         nr_fmat,     &
                         isymop_fmat, &
                         ihrmop)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: mat_dim
    integer, intent(in) :: nz
    integer, intent(in) :: nr_dmat_p
    integer, intent(in) :: nr_fmat
    real(8)             :: dmat_p(mat_dim, mat_dim, nz, nr_dmat_p)
    integer             :: isymop_dmat(nr_dmat_p)
    integer             :: isymop_fmat(nr_fmat)
    integer             :: ihrmop(nr_dmat_p)
!   ----------------------------------------------------------------------------

    call dirac_mpi_bcast(dmat_p)
    call dirac_mpi_bcast(isymop_dmat)
    call dirac_mpi_bcast(isymop_fmat)
    call dirac_mpi_bcast(ihrmop)

  end subroutine
#endif /* ifdef VAR_MPI */


  subroutine lr(block_length,       &
                mat_dim,            &
                nz,             &
                nr_mat,             &
                dmat_perturbed,     &
                fmat,               &
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
                buffer)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(in)    :: dmat_perturbed(*)
    real(8), intent(inout) :: fmat(*)
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
!   ----------------------------------------------------------------------------
    logical                :: calculate_qr, calculate_qi

    real(8)                :: u_qr,    u_qi(3)
    real(8)                :: v_qr(3), v_qi(3, 3)

    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)

    real(8)                :: z_b
    real(8)                :: y_b(3)

    real(8)                :: derv(mxbllen, 0:99)
    real(8)                :: ddot
    integer                :: k_start, k, im, im_off, j
!   ----------------------------------------------------------------------------

    call get_functional_derv(2,                  &
                             block_length,       &           
                             w,                  &
                             r_0,                &
                             z_0,                &  
                             nr_nonzero_derv(2), &
                             derv)

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

        call get_n(r_b,                        &
                   isym(im) - 1,               &
                   mat_dim,                    &
                   dmat_perturbed(1 + im_off), &
                   buffer,                     &
                   ao_p)
        
        if (is_gga_qr) then
          if (allow_nonhermitian) then
            call get_gn_nonhermitian(gr_b,                       &
                                     isym(im) - 1,               &
                                     mat_dim,                    &
                                     dmat_perturbed(1 + im_off), &
                                     buffer,                     &
                                     ao_p,                       &
                                     ao_p(ao_start(1, 0)))
          else
            call get_gn(gr_b,                       &
                        isym(im) - 1,               &
                        mat_dim,                    &
                        dmat_perturbed(1 + im_off), &
                        buffer,                     &
                        ao_p,                       &
                        ao_p(ao_start(1, 0)))
          end if
        end if

        u_qr = derv(1, d_rr)*r_b
        
        if (is_gga_qr) then
          z_b  = 2.0d0*ddot(3, gr_0, 1, gr_b, 1)
          u_qr = u_qr + derv(1, d_rz)*z_b

          do j = 1, 3
            v_qr(j) = gr_0(j)*(                   & 
                                r_b*derv(1, d_rz) &
                              + z_b*derv(1, d_zz) &   
                              )                   &
                    + gr_b(j)*derv(1, d_z)
          end do
          v_qr = 2.0d0*v_qr
        end if

        if (is_gga_qr) then
          call gga_distribute_r(u_qr,             &
                                v_qr,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p,             &
                                ao_p(ao_start(1, 0)))
        else
          call lda_distribute_r(u_qr,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p)
        end if

      end if !end of qr contribution

!     antihermitian contribution (qi = quaternion imaginary)
      if (calculate_qi) then

        call get_s(s_b,                        &
                   isym(im) - 1,               &
                   mat_dim,                    &
                   dmat_perturbed(1 + im_off), &
                   buffer,                     &
                   ao_p)
        
        if (is_gga_qi) then
          if (allow_nonhermitian) then
            call get_gs_nonhermitian(gs_b,                       &
                                     isym(im) - 1,               &
                                     mat_dim,                    &
                                     dmat_perturbed(1 + im_off), &
                                     buffer,                     &
                                     ao_p,                       &
                                     ao_p(ao_start(1, 0)))
          else
            call get_gs(gs_b,                       &
                        isym(im) - 1,               &
                        mat_dim,                    &
                        dmat_perturbed(1 + im_off), &
                        buffer,                     &
                        ao_p,                       &
                        ao_p(ao_start(1, 0)))
          end if
        end if

!       initialize otherwise there can be numerical dirt
!       and then collinear numbers are wrong
        u_qi = 0.0d0
        v_qi = 0.0d0

        do k = k_start, 3

          u_qi(k) = derv(1, d_ss)*s_b(k)
         
          if (is_gga_qi) then
            y_b(k)  = ddot(3, gr_0, 1, gs_b(1, k), 1)
            u_qi(k) = u_qi(k) + derv(1, d_sy)*y_b(k)

            do j = 1, 3
              v_qi(j, k) = gr_0(j)*(                   & 
                                  s_b(k)*derv(1, d_sy) &
                                + y_b(k)*derv(1, d_yy) &   
                                )                      &
                      + 2.0d0*gs_b(j, k)*derv(1, d_x)
            end do
          end if
        end do

        if (is_gga_qi) then
          call gga_distribute_i(u_qi,             &
                                v_qi,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p,             &
                                ao_p(ao_start(1, 0)))
        else
          call lda_distribute_i(u_qi,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p)
        end if

      end if !end of qi contribution

    end do

  end subroutine


! gosia: temporary fix:
  subroutine g_lao_contrib(block_length,       &
                mat_dim,            &
                nz,             &
                nr_mat,             &
                dmat_perturbed,     &
                fmat,               &
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
                buffer)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(in)    :: dmat_perturbed(*)
    real(8), intent(inout) :: fmat(*)
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
!   ----------------------------------------------------------------------------
    logical                :: calculate_qr, calculate_qi

    real(8)                :: u_qr,    u_qi(3)
    real(8)                :: v_qr(3), v_qi(3, 3)

    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)

    real(8)                :: z_b
    real(8)                :: y_b(3)

    real(8)                :: derv(mxbllen, 0:99)
    real(8)                :: ddot
    integer                :: k_start, k, im, im_off, j, off_mag
!   ----------------------------------------------------------------------------

    call get_functional_derv(2,                  &
                             block_length,       &           
                             w,                  &
                             r_0,                &
                             z_0,                &  
                             nr_nonzero_derv(2), &
                             derv)

    do im = 1, nr_mat

      im_off = mat_dim*mat_dim*nz*(im - 1)

        call get_n(r_b,                        &
                   isym(im) - 1,               &
                   mat_dim,                    &
                   dmat_perturbed(1 + im_off), &
                   buffer,                     &
                   ao_p)
        
        u_qr = derv(1, d_rr)*r_b
        
          call lda_distribute_r(u_qr,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p)

    end do

  end subroutine



end module
