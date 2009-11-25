module dft_main

  use memory_allocator
  use dft_cfg
  use dft_temp
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

contains

      subroutine dftexc(dmat, &
                        dmat_perturbed, &
                        fmat, &
                        nr_dmat_perturbed, &
                        nr_fmat, &
                        isymop, &
                        ihrmop, &
                        dft_options, &
                        iprint)

!     written by T. Helgaker
!     DIRAC adaption by T. Saue

!     --------------------------------------------------------------------------
      real(8), intent(in) :: dmat(nr_ao_total, nr_ao_total, nr_qu, *)
      real(8)             :: dmat_perturbed(nr_ao_total, nr_ao_total, nr_qu, *)
      real(8)             :: fmat(nr_ao_total, nr_ao_total, nr_qu, *)
      integer, intent(in) :: nr_dmat_perturbed
      integer, intent(in) :: nr_fmat
      integer, intent(in) :: isymop(nr_dmat_perturbed)
      integer, intent(in) :: ihrmop(nr_dmat_perturbed)
      type(dft_options_t) :: dft_options
      integer, intent(in) :: iprint
!     --------------------------------------------------------------------------
      integer             :: iz, im, iq, g_order, m_order
      real(8)             :: time_start, time_in_xc, second
      real(8)             :: ddot

      real(8), allocatable :: dmat_saop_outer(:)
      real(8), allocatable :: dmat_saop_resp(:)

      integer              :: max_ao_g_order, max_ao_m_order
!     --------------------------------------------------------------------------

      call init_dft_cb()

!     start timer
      time_start = second()

!     initial expectation value of the two-electron matrix
      if (dft_options%energy) then
        dft_temp_xc_mat_energy = ddot(nr_ao_total*nr_ao_total, dmat, 1, fmat, 1)
      end if

!     get average density from open shells
      if (nr_open_shells > 0) call gtdoav(dmat, 1)

!     scale
      call dscal(nr_ao_total*nr_ao_total, 2.0d0, dmat, 1)


!     scale and insert half-phases

      if (dft_options%lr  &
     &    .or. dft_options%qr &
     &    .or. (dft_options%aor_xc_order > 0)) then

        call dscal(nr_ao_total*nr_ao_total*nr_qu*nr_dmat_perturbed, &
     &             2.0d0, dmat_perturbed, 1)

        if (nr_qu < 4) then
          do iz = 1, nr_qu
            do im = 1, nr_dmat_perturbed
              iq = dft_cb_pq_to_uq(iz, isymop(im) - 1)
              call q2bphase('D', iq, 1, dmat_perturbed(1, 1, iz, im))
            end do
            do im = 1, nr_fmat
              iq = dft_cb_pq_to_uq(iz, isymop(im) - 1)
              call q2bphase('F', iq, 1, fmat(1, 1, iz, im))
            end do
          end do
        end if

      end if

!     calculate grid
      if (.not. dft_temp_grid_is_done) then
        call dftgrd(dmat, 1)
        dft_temp_grid_is_done = .true.
        if (dft_cfg_keep_ao .and. ao_mem_allocated) dft_temp_reload_batches = .true.
      end if

!     allocate and get saop density matrices
      if (dft_cfg_saop) then
        dft_cfg_saop_is_active = .true.
        call alloc(dmat_saop_outer, nr_ao_total*nr_ao_total)
        if (dft_cfg_saop_with_response_part) then
          call alloc(dmat_saop_resp, nr_ao_total*nr_ao_total)
        end if
        call dftac_saop_get_dmats(dmat_saop_outer, &
                                  dmat_saop_resp)
      end if


!     check whether gga code structure will be used

      dft_temp_gga_hs = is_gga_functional
      dft_temp_gga_ha = is_gga_functional

#ifdef MOD_AUTODIFF
      if (dft_cfg_autodiff) then
!       xc_fun_set_functional is called several times
!       this is because gga query is not ideally structured
        call xc_fun_set_functional(dft_options%aor_xc_order)
        dft_temp_gga_hs = dft_cfg_is_gga_autodiff
        dft_temp_gga_ha = dft_cfg_is_gga_autodiff
      end if
#endif

      if (dft_options%lr .or. dft_options%qr) then
        if (dft_cfg_alda_hs) dft_temp_gga_hs = .false.
        if (dft_cfg_alda_ha) dft_temp_gga_ha = .false.
      end if

!     in QR it is not so easy to separate H+ and H-
!     contributions since they may mix
!     therefore i quit here if the user wants QR
!     and not the full ALDA or XALDA
!     it is really unfair to quit as late as here
!     this can be checked earlier
!     but this is anyway an exotic option
      if (dft_options%qr) then
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
        call dftpar_init1(dft_options, nr_dmat_perturbed, nr_fmat)
        call dftpar_init2(dft_options, nr_dmat_perturbed, isymop, ihrmop, &
     &                    dmat, dmat_perturbed, dmat_saop_resp, dmat_saop_outer)
      end if
#endif /* ifdef VAR_MPI */

    max_ao_g_order = 0
!radovan: need to test vrg first
!   if (dft_temp_gga_hs .or. dft_temp_gga_ha .or. dft_cfg_cdft_vrg) then
    if (dft_temp_gga_hs .or. dft_temp_gga_ha) then
      max_ao_g_order = 1
    end if

    max_ao_m_order = 0
    if (dft_options%london) then
      max_ao_m_order = 1
    end if

#ifdef VAR_MPI
    if (is_parallel_calculation) then
      call dirac_mpi_bcast(max_ao_g_order)
      call dirac_mpi_bcast(max_ao_m_order)
    end if
#endif

    call dft_ao_init(is_parallel_calculation, &
                     max_ao_g_order,          &
                     max_ao_m_order,          &
                     dft_cfg_screening,       &
                     dft_cfg_dfthri)
      
!     integrate
      call dftdrv(dmat,              &
                  dmat_saop_resp,    &
                  dmat_saop_outer,   &
                  dmat_perturbed,    &
                  fmat,              &
                  nr_dmat_perturbed, &
                  nr_fmat,           &
                  isymop,            &
                  ihrmop,            &
                  dft_options,       &
                  dft_cfg_print_level)

!     deallocate saop density matrices
      if (dft_cfg_saop) then
        call dealloc(dmat_saop_outer)
        if (dft_cfg_saop_with_response_part) then
          call dealloc(dmat_saop_resp)
        end if
      end if

!     rescale
      call dscal(nr_ao_total*nr_ao_total, 0.5d0, dmat, 1)

!     get average density from open shells
      if (nr_open_shells > 0) call gtdoav(dmat, -1)

!     in case of gga symmetrize matrices
      if (dft_temp_gga_hs .or. dft_temp_gga_ha) then
        call gga_sym(nr_fmat, fmat, isymop, ihrmop)
      end if

!     rescale and insert half-phases

      if (dft_options%lr  &
     &    .or. dft_options%qr &
     &    .or. (dft_options%aor_xc_order > 0)) then

        if (nr_qu < 4) then
          do iz = 1, nr_qu
            do im = 1, nr_dmat_perturbed
              iq = dft_cb_pq_to_uq(iz, isymop(im) - 1)
              call q2bphase('D', iq, 1, dmat_perturbed(1, 1, iz, im))
            end do
            do im = 1, nr_fmat
              iq = dft_cb_pq_to_uq(iz, isymop(im) - 1)
              call q2bphase('F', iq, 1, fmat(1, 1, iz, im))
            end do
          end do
        end if

        call dscal(nr_ao_total*nr_ao_total*nr_qu*nr_dmat_perturbed, &
     &             0.5d0, dmat_perturbed, 1)

      end if

!     calculate final expectation value of two-electron matrix
      if (dft_options%energy) then
!       the xc contribution can then be subtracted 
        dft_temp_xc_mat_energy = ddot(nr_ao_total*nr_ao_total, dmat, 1, fmat, 1) &
                               - dft_temp_xc_mat_energy
      end if 

!     timing
      time_in_xc = second() - time_start
      write(print_unit, '(1x)')
      call timtxt('>>>> time used in dftexc:', time_in_xc, print_unit)
      write(print_unit, '(1x)')

      end subroutine





!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DFTDRV(DMAT,DMAT_RESP_AO,DMAT_OUTER_AO, &
                        DTRMAT,EXCMAT,nr_dmat_perturbed, nr_fmat, &
                        ISYMOP,IHRMOP,dft_options, &
                        IPRINT)
!
!     T. Helgaker sep 99
!

      real(8) :: EXCMAT(*)
      real(8) :: DMAT(*),DMAT_RESP_AO(*),DMAT_OUTER_AO(*), &
     &          DTRMAT(*)
      integer    ::ISYMOP(nr_dmat_perturbed),IHRMOP(nr_dmat_perturbed)
      integer    ::iprint, nr_dmat_perturbed, nr_fmat, ibatch
      type(dft_options_t) :: dft_options

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

#ifdef MOD_AUTODIFF
    if (dft_cfg_autodiff) then
      call xc_fun_set_functional(dft_options%aor_xc_order)
    end if
#endif

    dft_temp_nr_electrons = 0.0d0
    dft_temp_xc_energy    = 0.0d0
    i_am_the_master = .false.
    if (mytid == mparid) i_am_the_master = .true.

    if (i_am_the_master) then
      rewind dft_temp_file_unit
      npnts  = 0
      ibatch = 0
    end if

    if (dft_cfg_keep_ao .and. ao_mem_allocated .and. dft_temp_reload_batches) then
      call kill_ao_in_mem()
      ao_mem_allocated        = .false.
      dft_temp_reload_batches = .false.
    end if

!*    Loop over batches of grid points
  100 CONTINUE

      if (i_am_the_master) then
        if (dft_cfg_import_grid) then
          read(dft_temp_file_unit, *) npoint
        else
          read(dft_temp_file_unit)    npoint
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
        
        call dftdr1(kx, ky, kz, kw, &
                    mynpoint, ibatch, &
                    DMAT,DMAT_RESP_AO,DMAT_OUTER_AO, &
                    DTRMAT,EXCMAT,nr_dmat_perturbed,nr_fmat,&
                    ISYMOP,IHRMOP, &
                    dft_options, &
                    IPRINT)

        call dealloc(kx)
        call dealloc(ky)
        call dealloc(kz)
        call dealloc(kw)

        GO TO 100
        end if !IF( NPOINT .lt. 0 ) THEN
 200    CONTINUE




#if defined (VAR_MPI)
      if (is_parallel_calculation) then
        call dftres(excmat, nr_fmat, dft_options, iprint)
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
      if (dft_options%energy) then
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
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DFTDR1(CORX,CORY,CORZ, &
     &                  WEIGHT, &
                        NPOINT, ibatch, &
     &                  DMAT,DMAT_RESP_AO,DMAT_OUTER_AO, &
     &                  DTRMAT,EXCMAT,nr_dmat_perturbed, nr_fmat, &
                        ISYMOP,IHRMOP, &
     &                  dft_options, &
     &                  IPRINT)
!***********************************************************************
!
!     Written by Trond Saue and T.Helgaker
!
!***********************************************************************


      integer :: npoint, iprint, nr_dmat_perturbed, nr_fmat, ibatch, ioff, i
      logical :: gga_qr
      real(8) :: corx(*)
      real(8) :: cory(*)
      real(8) :: corz(*)
      real(8) :: weight(*)
      real(8) :: dmat(*)
      real(8) :: dmat_resp_ao(*)
      real(8) :: dmat_outer_ao(*)
      real(8) :: dtrmat(*)
      real(8) :: excmat(*)
      integer :: isymop(nr_dmat_perturbed)
      integer :: ihrmop(nr_dmat_perturbed)
      real(8) :: energy_density

      real(8) :: RHG(3),VT(3)
      real(8) :: derv(mxbllen, 0:99)
      real(8) :: r_0(mxbllen)
      real(8) :: z_0(mxbllen)
      real(8) :: w(mxbllen)
      real(8) :: rho13, rho_outer, rho_resp
      type(dft_options_t) :: dft_options

      integer :: ipoint, need_limit
      real(8) :: ddot, f

      integer :: itest, jtest

!   this pointer is useful for .KEEPAO, then it can point to memory
!   otherwise it points to ao
    real(8), pointer             :: ao_p(:)

    real(8), allocatable, target :: ao(:)
    real(8), allocatable         :: buffer(:)

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
!       here get magnetic dervs of aos
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
      call get_n(r_0(1), 0, dmat, buffer, ao_p)

      if (r_0(1) > dft_cfg_dfthr0) then

!       set it to zero, otherwise it is undefined with lda
        z_0 = 0.0d0

!       density gradient
        if (need_ao_order(1, 0)) then
          call get_gn(rhg, 0, dmat, buffer, ao_p, ao_p(ao_start(1, 0):ao_end(1, 0)))
          z_0(1) = ddot(3, rhg, 1, rhg, 1)
        end if

        w(1) = weight(ipoint)

!       number of electrons
        dft_temp_nr_electrons = dft_temp_nr_electrons + w(1)*r_0(1)

!       energy density
        if (dft_options%energy) then

          call get_functional_derv(0,            &
                                   block_length, &
                                   w,            &
                                   r_0,          &
                                   z_0,          &
                                   0,            &
                                   derv)

          dft_temp_xc_energy = dft_temp_xc_energy + derv(1, d_0)

        end if


!       exchange-correlation potential 
!       ==============================

        if (dft_options%energy &
     &      .or. dft_options%gradient &
     &      .or. dft_options%london) then

          call get_functional_derv(1,                  &
                                   block_length,       &
                                   w,                  &
                                   r_0,                &
                                   z_0,                &
                                   nr_nonzero_derv(1), &
                                   derv)

          if (dft_temp_gga_hs) then
!           this factor has nothing to do
!           with symmetrization later
!           it comes from the fact that derivative of z
!           yields two identical terms (one left, one right)
            derv(1, d_z) = 2.0d0*derv(1, d_z)
          end if

!         asymptotic correction
          if (dft_cfg_saop_is_active) then
            call get_n(rho_outer,     &
                       0,             &
                       dmat_outer_ao, &
                       buffer,           &
                       ao_p)
            if (dft_cfg_saop_with_response_part) then
              call get_n(rho_resp,     &
                         0,            &
                         dmat_resp_ao, &
                         buffer,          &
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
        if (dft_options%energy) then
!       here we will add LAO contribution to KS matrix
          if (dft_temp_gga_hs) then
            vt = derv(1, d_z)*rhg

#ifdef TEST_GAO1_SYMMETRY
            if (abs(corx(ipoint)) > 1.0d-10) then
              vt(2) = vt(1)*cory(ipoint)/corx(ipoint)
              vt(3) = vt(1)*corz(ipoint)/corx(ipoint)
            end if
#endif

            call gga_distribute(derv(1, d_r), vt, 0, excmat, ao_p, ao_p(ao_start(1, 0):ao_end(1, 0)))
          else
            call lda_distribute(derv(1, d_r),     0, excmat, ao_p)
          end if
        end if


!       linear response
        if(dft_options%lr) then

!#define DEBUG_XC_RESPONSE
#ifndef DEBUG_XC_RESPONSE
          call lr(block_length,    &
                  nr_dmat_perturbed,          &
                  dtrmat,          &
                  excmat,          &
                  isymop,          &
                  ihrmop,          &
                  .false.,         &
                  ao_p,            &
                  r_0,             &
                  rhg,             &
                  z_0,             &
                  w,               &
                  dft_temp_gga_hs, &
                  dft_temp_gga_ha, &
                  buffer)
#else
            call xc_response_contribution(2,  &
     &                      nr_dmat_perturbed,  &
     &                      dft_temp_gga_hs, &
     &                      dtrmat,  &
     &                      excmat, &
     &                      isymop, &
     &                      ao_p(ao_start(0, 0):ao_end(0, 0)), &
     &                      ao_p(ao_start(1, 0):ao_end(1, 0)), &
     &                      r_0(1), &
     &                      rhg, &
     &                      w(1), &
     &                      buffer)
#endif
        endif


!       quadratic response
        if(dft_options%qr) then

          gga_qr = dft_temp_gga_hs .or. dft_temp_gga_ha
          call sdft_qr(block_length,    &
                       gga_qr,          &
                       dft_cfg_no_sdft, &
                       dtrmat,excmat,   &
                       isymop,ihrmop,ao_p,ao_p(ao_start(1, 0)), &
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
                  nr_dmat_perturbed,          &
                  dtrmat,          &
                  excmat,          &
                  isymop,          &
                  ihrmop,          &
                  .true.,          &
                  ao_p,            &
                  r_0,             &
                  rhg,             &
                  z_0,             &
                  w,               &
                  dft_temp_gga_hs, &
                  dft_temp_gga_ha, &
                  buffer)

        endif


#ifdef MOD_AUTODIFF
#ifdef MOD_AORESPONSE
          if (dft_options%aor_xc_order > 0) then
            call xc_response_contribution(dft_options%aor_xc_order,  &
     &                      nr_dmat_perturbed,  &
     &                      dft_temp_gga_hs, &
     &                      dtrmat,  &
     &                      excmat, &
     &                      isymop, &
     &                      ao_p(ao_start(0, 0):ao_end(0, 0)), &
     &                      ao_p(ao_start(1, 0):ao_end(1, 0)), &
     &                      r_0(1), &
     &                      rhg, &
     &                      w(1), &
     &                      buffer)
          end if
#endif /* ifdef MOD_AORESPONSE */
#endif /* ifdef MOD_AUTODIFF */


      end if !if (r_0(1) > dft_cfg_dfthr0) then
    end do !do ipoint = 1, npoint

    nullify(ao_p)

    call dealloc(ao)
    call dealloc(buffer)

    end subroutine
      


#ifdef VAR_MPI
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DFTRES(EXCMAT,nr_fmat,dft_options, &
     &                  IPRINT)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!     Send receive DFTresults
!     Written by T.Saue and T.helgaker
!
!***********************************************************************

      integer :: nr_fmat, iprint
      real(8) :: EXCMAT(*)
      type(dft_options_t) :: dft_options
      real(8), allocatable :: kfbuffer(:)
      integer :: ierr, LFMAT
      real(8) :: sum_nr_electrons
      real(8) :: sum_xc_energy
!
!

      call mpi_reduce(dft_temp_nr_electrons, &
                      sum_nr_electrons,      &
                      1,                     &
                      mpi_double_precision,  &
                      mpi_sum,               &
                      0,                     &
                      mpi_comm_world,        &
                      ierr) 
      if (dft_options%energy) then
        call mpi_reduce(dft_temp_xc_energy,   &
                        sum_xc_energy,        &
                        1,                    &
                        mpi_double_precision, &
                        mpi_sum,              &
                        0,                    &
                        mpi_comm_world,       &
                        ierr) 
      end if
      IF(MYTID.EQ.MPARID) THEN
        dft_temp_nr_electrons = sum_nr_electrons
        if (dft_options%energy) then
          dft_temp_xc_energy = sum_xc_energy
        end if
      ENDIF

      LFMAT = nr_ao_total*nr_ao_total*nr_fmat*nr_qu
      IF(MYTID.EQ.MPARID) THEN
        call alloc(kfbuffer, lfmat)
        kfbuffer = 0.0d0
      end if
      call mpi_reduce(excmat,                &
                      kfbuffer,                 &
                      lfmat,                 &
                      mpi_double_precision,  &
                      mpi_sum,               &
                      0,                     &
                      mpi_comm_world,        &
                      ierr) 
      IF(MYTID.EQ.MPARID) THEN
        CALL DCOPY(LFMAT,KFbuffer,1,EXCMAT,1)     
        call dealloc(kfbuffer)
      ENDIF

      end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DFTNOD()
!***********************************************************************
!     Based on HERNOD 
!     written by T. Saue March 4
!
!    This is the node program for DFT                           *
!***********************************************************************


!
      type(dft_options_t) :: dft_options
      integer, allocatable :: ksym(:)
      integer, allocatable :: khrm(:)
      real(8), allocatable :: kexc(:)
      real(8), allocatable :: kden(:)
      real(8), allocatable :: kdtr(:)
      real(8), allocatable :: KDMAT_OUTER_AO(:)
      real(8), allocatable :: KDMAT_RESP_AO(:)
      integer              :: nr_dmat_perturbed, nr_fmat
      real(8)              :: wall1, cpu1
      integer              :: iprint
      integer              :: max_ao_g_order, max_ao_m_order
!
!
      PARHER = .TRUE.
!
!     Start timing.
!
      CALL GETTIM(CPU1,WALL1)
!
!     This is a slave in a parallel run
!
      SLAVE = .TRUE.
      DEBUG  = .FALSE.
!
!     Receive initialization data from master 1
!     =========================================
      call sync_dft_cb()
      call dftpar_init1(dft_options, nr_dmat_perturbed, nr_fmat)
      call alloc(kexc, nr_ao_total*nr_ao_total*nr_qu*nr_fmat)
      kexc = 0.0d0
      call alloc(ksym, nr_dmat_perturbed)
      call alloc(khrm, nr_dmat_perturbed)
      call alloc(kden,nr_ao_total*nr_ao_total)
      if (dft_options%lr .or. dft_options%qr .or. (dft_options%aor_xc_order > 0)) then
        call alloc(kdtr, nr_ao_total*nr_ao_total*nr_qu*nr_dmat_perturbed)
      ENDIF
!.....Allocate density matrices for SAOP
      if (dft_cfg_saop_is_active) then
        call alloc(KDMAT_OUTER_AO,nr_ao_total*nr_ao_total)
        if (dft_cfg_saop_with_response_part) then
          call alloc(KDMAT_RESP_AO,nr_ao_total*nr_ao_total)
        ENDIF
      ENDIF
 
!     Receive initialization data from master 2
!     =========================================
      CALL DFTPAR_INIT2(dft_options,nr_dmat_perturbed, ksym, khrm, &
     &                  KDEN,KDTR, &
     &                  KDMAT_RESP_AO,KDMAT_OUTER_AO)

    call dirac_mpi_bcast(max_ao_g_order)
    call dirac_mpi_bcast(max_ao_m_order)
    call dft_ao_init(.true.,            &
                     max_ao_g_order,    &
                     max_ao_m_order,    &
                     dft_cfg_screening, &
                     dft_cfg_dfthri)

!.....Initialization
      IPRINT = 0
!.....Integrate
      CALL DFTDRV(KDEN,KDMAT_RESP_AO,KDMAT_OUTER_AO, &
     &            KDTR,KEXC,nr_dmat_perturbed, nr_fmat, ksym, khrm, &
     &            dft_options,IPRINT)

      call dealloc(KEXC)
      call dealloc(ksym)
      call dealloc(khrm)
      call dealloc(kden)
      if (dft_options%lr .or. dft_options%qr .or. (dft_options%aor_xc_order > 0)) then
        call dealloc(kdtr)
      ENDIF
      if (dft_cfg_saop_is_active) then
        call dealloc(KDMAT_OUTER_AO)
        if (dft_cfg_saop_with_response_part) then
          call dealloc(KDMAT_RESP_AO)
        ENDIF
      ENDIF

      end subroutine


      subroutine dftpar_init1(dft_options, nr_dmat_perturbed, nr_fmat)

!     written by T. Saue april 7 2009


!     --------------------------------------------------------------------------
      type(dft_options_t)                :: dft_options
      integer                            :: nr_dmat_perturbed, nr_fmat
!     --------------------------------------------------------------------------
      logical                            :: i_am_the_master
      integer                            :: ierr
!     --------------------------------------------------------------------------

      if (mytid == mparid) then
        i_am_the_master = .true.
      else
        i_am_the_master = .false.
      end if

!     sync nr_dmat_perturbed, nr_fmat
      call dirac_mpi_bcast(nr_dmat_perturbed)
      call dirac_mpi_bcast(nr_fmat)

!     sync various parameters
      call sync_dft_cfg()
      call sync_dft_temp()
      call sync_dft_options(dft_options)

!     sync functionals
      call dftsyncfunc(i_am_the_master)
 
      end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DFTPAR_INIT2(dft_options,nr_dmat_perturbed,ISYMOP,IHRMOP, &
     &                       DMAT,DTRMAT,DMAT_RESP_AO,DMAT_OUTER_AO)
!***********************************************************************
!
!     Initialize parallel DFT run
!
!     Written by T. Saue April 7 2009
! 
!***********************************************************************
      real(8) :: DMAT(nr_ao_total, nr_ao_total),DTRMAT(*), &
     &          DMAT_RESP_AO(nr_ao_total, nr_ao_total), &
     &          DMAT_outer_AO(nr_ao_total, nr_ao_total)
      integer ::          ISYMOP(nr_dmat_perturbed),IHRMOP(nr_dmat_perturbed)
      type(dft_options_t) :: dft_options
      integer :: nr_dmat_perturbed, lmat, ierr
      call dirac_mpi_bcast(isymop)
      call dirac_mpi_bcast(ihrmop)
      call dirac_mpi_bcast(dmat)
      if (dft_options%lr .or. dft_options%qr .or. (dft_options%aor_xc_order > 0)) then
        lmat = nr_ao_total*nr_ao_total*nr_qu*nr_dmat_perturbed
        CALL MPI_BCAST(DTRMAT,LMAT,MPI_DOUBLE_PRECISION,MPARID, &
     &                   MPI_COMM_WORLD,IERR)
      ENDIF

      if (dft_cfg_saop_is_active) then
        call dirac_mpi_bcast(dmat_outer_ao)
        if (dft_cfg_saop_with_response_part) then
          call dirac_mpi_bcast(dmat_resp_ao)
        ENDIF
      ENDIF
      end subroutine
#endif /* ifdef VAR_MPI */


  subroutine lr(block_length,       &
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

      im_off = nr_ao_total*nr_ao_total*nr_qu*(im - 1)

!     hermitian contribution (qr = quaternion real)
      if (calculate_qr) then

        call get_n(r_b,                        &
                   isym(im) - 1,               &
                   dmat_perturbed(1 + im_off), &
                   buffer,                     &
                   ao_p)
        
        if (is_gga_qr) then
          if (allow_nonhermitian) then
            call get_gn_nonhermitian(gr_b,                       &
                                     isym(im) - 1,               &
                                     dmat_perturbed(1 + im_off), &
                                     buffer,                     &
                                     ao_p,                       &
                                     ao_p(ao_start(1, 0)))
          else
            call get_gn(gr_b,                       &
                        isym(im) - 1,               &
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
          call gga_distribute(u_qr,             &
                              v_qr,             &
                              isym(im) - 1,     &
                              fmat(1 + im_off), &
                              ao_p,             &
                              ao_p(ao_start(1, 0)))
        else
          call lda_distribute(u_qr,             &
                              isym(im) - 1,     &
                              fmat(1 + im_off), &
                              ao_p)
        end if

      end if !end of qr contribution

!     antihermitian contribution (qi = quaternion imaginary)
      if (calculate_qi) then

        call get_s(s_b,                        &
                   isym(im) - 1,               &
                   dmat_perturbed(1 + im_off), &
                   buffer,                     &
                   ao_p)
        
        if (is_gga_qi) then
          if (allow_nonhermitian) then
            call get_gs_nonhermitian(gs_b,                       &
                                     isym(im) - 1,               &
                                     dmat_perturbed(1 + im_off), &
                                     buffer,                     &
                                     ao_p,                       &
                                     ao_p(ao_start(1, 0)))
          else
            call get_gs(gs_b,                       &
                        isym(im) - 1,               &
                        dmat_perturbed(1 + im_off), &
                        buffer,                     &
                        ao_p,                       &
                        ao_p(ao_start(1, 0)))
          end if
        end if

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
          call gga_distribute_qi(u_qi,             &
                                 v_qi,             &
                                 isym(im) - 1,     &
                                 fmat(1 + im_off), &
                                 ao_p,             &
                                 ao_p(ao_start(1, 0)))
        else
          call lda_distribute_qi(u_qi,             &
                                 isym(im) - 1,     &
                                 fmat(1 + im_off), &
                                 ao_p)
        end if

      end if !end of qi contribution

    end do

  end subroutine

end module
