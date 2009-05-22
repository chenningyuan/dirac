module dft_main

  use memory_allocator
  use dft_cfg
  use dft_temp
  use dft_ao_in_mem
#ifdef MOD_AORESPONSE
  use aor_xc
#endif

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


      implicit none

#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "dcbbas.h"

#ifdef VAR_MPI
! parcal in dcbgen.h
#include "dcbgen.h"
#include "infpar.h"
#include "mpif.h"
#endif

#include "nuclei.h"
#include "dcborb.h"
#include "dgroup.h"
#include "functionals.h"
#include "densfit.h"

!     --------------------------------------------------------------------------
      real(8), intent(in) :: dmat(n2bbasx, nz, *)
      real(8)             :: dmat_perturbed(n2bbasx, nz, *)
      real(8)             :: fmat(n2bbasx, nz, *)
      integer, intent(in) :: nr_dmat_perturbed
      integer, intent(in) :: nr_fmat
      integer, intent(in) :: isymop(nr_dmat_perturbed)
      integer, intent(in) :: ihrmop(nr_dmat_perturbed)
      type(dft_options_t) :: dft_options
      integer, intent(in) :: iprint
!     --------------------------------------------------------------------------
      integer             :: iz, im, iq
      real(8)             :: time_start, time_in_xc, second
      real(8)             :: ddot

      real(8), allocatable :: xc_vector_b(:)
      real(8), allocatable :: dmat_saop_outer(:)
      real(8), allocatable :: dmat_saop_resp(:)

#ifdef VAR_MPI
#ifndef INT_STAR8
      integer             :: istat
#endif /* ifndef INT_STAR8 */
      dimension              istat(mpi_status_size)
#endif /* ifdef VAR_MPI */
!     --------------------------------------------------------------------------


!     start timer
      time_start = second()

!     prepare for accumulating integrals with density fitting
#ifdef MOD_DNF
      if (df_fitvxc) then
        call alloc(xc_vector_b, nbasisaux)
        xc_vector_b = 0.0d0
        if (dft_options%energy) then 
          call save_matrix(xc_vector_b, 1, nbasisaux, 'FITVXC')
        else if (dft_options%lr) then 
          call save_matrix(xc_vector_b, 1, nbasisaux, 'LRDFVX')
        end if
      end if
#endif

!     initial expectation value of the two-electron matrix
#ifdef MOD_DNF
      if (.not. df_fitvxc) then
#endif
      if (dft_options%energy) then
        dft_temp_xc_mat_energy = ddot(n2bbasx, dmat, 1, fmat, 1)
      end if
#ifdef MOD_DNF
      end if
#endif

!     get average density from open shells
      if (nasht > 0) call gtdoav(dmat, 1)

!     scale
      call dscal(n2bbasx, 2.0d0, dmat, 1)


!     scale and insert half-phases

#ifdef MOD_DNF
      if (.not. df_fitvxc) then
#endif
      if (dft_options%lr  &
     &    .or. dft_options%qr &
     &    .or. (dft_options%aor_xc_order > 0)) then

        call dscal(n2bbasxq*nr_dmat_perturbed, &
     &             2.0d0, dmat_perturbed, 1)

        if (nz < 4) then
          do iz = 1, nz
            do im = 1, nr_dmat_perturbed
              iq = ipqtoq(iz, isymop(im) - 1)
              call q2bphase('D', iq, 1, dmat_perturbed(1, iz, im))
            end do
            do im = 1, nr_fmat
              iq = ipqtoq(iz, isymop(im) - 1)
              call q2bphase('F', iq, 1, fmat(1, iz, im))
            end do
          end do
        end if

      end if
#ifdef MOD_DNF
      end if
#endif

!     calculate grid
      if (.not. dft_temp_grid_is_done) then
        call dftgrd(dmat, 1)
        dft_temp_grid_is_done = .true.
      end if

!     allocate and get saop density matrices
      if (dft_cfg_saop) then
        dft_cfg_saop_is_active = .true.
        call alloc(dmat_saop_outer, n2bbasx)
        if (dft_cfg_saop_with_response_part) then
          call alloc(dmat_saop_resp, n2bbasx)
        end if
        call dftac_saop_get_dmats(dmat_saop_outer, &
                                  dmat_saop_resp)
      end if


!     check whether gga code structure will be used

      dft_temp_gga_hs = dftisgga()
      dft_temp_gga_ha = dft_temp_gga_hs

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


!     check whether derivatives of basis functions are needed

!     it should already be set to .false. but just to be sure
      dft_temp_need_gao1 = .false.
      dft_temp_need_gao2 = .false.

      if (dft_temp_gga_hs .or. dft_temp_gga_ha) then
        dft_temp_need_gao1 = .true.
      end if

      if (dft_cfg_cdft_vrg) then
        dft_temp_need_gao1 = .true.
      end if

#ifdef VAR_MPI
      if (parcal) then

        call dirac_parctl(4)

        if (nodes > maxnod) then
          call quit('number of nodes is too large, nodes > maxnod')
        end if

        write(lupri, '(/5x, a, i5/)') &
                     'number of nodes in parallel DFT:', nodes

        timing = .false.

        call dftpar_init1(dft_options, nr_dmat_perturbed)
        call dftpar_init2(dft_options, nr_dmat_perturbed, isymop, ihrmop, &
     &                    dmat, dmat_perturbed, dmat_saop_resp, dmat_saop_outer)

      end if
#endif /* ifdef VAR_MPI */

!     integrate
      call dftdrv(dmat,              &
                  dmat_saop_resp,    &
                  dmat_saop_outer,   &
                  dmat_perturbed,    &
                  fmat,              &
                  nr_dmat_perturbed, &
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
      call dscal(n2bbasx, 0.5d0, dmat, 1)

!     get average density from open shells
      if (nasht > 0) call gtdoav(dmat, -1)

#ifdef MOD_DNF
      if (.not. df_fitvxc) then
#endif

!     in case of gga symmetrize matrices
      if (dft_temp_gga_hs) then
        call sym_hs(fmat, nr_fmat, isymop, ihrmop)
      end if
      if (dft_temp_gga_ha .and. dft_cfg_sdft) then
        call sym_ha(fmat, nr_fmat, isymop, ihrmop)
      end if


!     rescale and insert half-phases

      if (dft_options%lr  &
     &    .or. dft_options%qr &
     &    .or. (dft_options%aor_xc_order > 0)) then

        if (nz < 4) then
          do iz = 1, nz
            do im = 1, nr_dmat_perturbed
              iq = ipqtoq(iz, isymop(im) - 1)
              call q2bphase('D', iq, 1, dmat_perturbed(1, iz, im))
            end do
            do im = 1, nr_fmat
              iq = ipqtoq(iz, isymop(im) - 1)
              call q2bphase('F', iq, 1, fmat(1, iz, im))
            end do
          end do
        end if

        call dscal(n2bbasxq*nr_dmat_perturbed, &
     &             0.5d0, dmat_perturbed, 1)

      end if

!     calculate final expectation value of two-electron matrix
      if (dft_options%energy) then
!       the xc contribution can then be subtracted 
        dft_temp_xc_mat_energy = ddot(n2bbasx, dmat, 1, fmat, 1) &
                               - dft_temp_xc_mat_energy
      end if 

#ifdef MOD_DNF
      else
!       df_fitvxc is .true.
        call dealloc(xc_vector_b)
      end if
#endif

!     timing
      time_in_xc = second() - time_start
      write(lupri, '(1x)')
      call timtxt('>>>> time used in dftexc:', time_in_xc, lupri)
      write(lupri, '(1x)')

      end subroutine





!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! /* Deck dftdrv */
      SUBROUTINE DFTDRV(DMAT,DMAT_RESP_AO,DMAT_OUTER_AO, &
     &                  DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,dft_options, &
     &                  IPRINT)
!
!     T. Helgaker sep 99
!


#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, &
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)
!
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#include "infpar.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#endif
#include "nuclei.h"
#include "orgcom.h"
#include "functionals.h"
#include "dgroup.h"
#include "densfit.h"
!
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
#endif
      DIMENSION EXCMAT(*)
      DIMENSION DMAT(*),DMAT_RESP_AO(*),DMAT_OUTER_AO(*), &
     &          DTRMAT(*),ISYMOP(NXCMAT),IHRMOP(NXCMAT)
      type(dft_options_t) :: dft_options

      real(8), allocatable :: gao(:)
      real(8), allocatable :: gao1(:)
      real(8), allocatable :: gao2(:)
      real(8), allocatable :: gab1(:)
      real(8), allocatable :: buf(:)
      integer, allocatable :: ncnt(:)

      real(8), allocatable :: kx(:)
      real(8), allocatable :: ky(:)
      real(8), allocatable :: kz(:)
      real(8), allocatable :: kw(:)

      integer, allocatable :: idisp(:)
      integer, allocatable :: isend(:)
!
!
!



!
!     Allocate arrays for basis functions and their derivatives at grid points
!
      NDIM = NTBAS(0)
#ifdef MOD_DNF
      IF (DF_FITVXC) NDIM = MAX(NTBAS(0),NBASISAUX)
#endif




      call alloc(gao, ndim)
      if (dft_temp_need_gao1) then
        call alloc(gao1, ndim*3)
      end if
      if (dft_temp_need_gao2) then
        call alloc(gao2, ndim*6)
      end if
      if (dft_options%london) then
        call alloc(gab1, ndim*3)
      end if
      call alloc(ncnt, ndim)
      call alloc(buf,  ndim*4)


!
!     integrate
!
      dft_temp_nr_electrons = 0.0d0
      dft_temp_xc_energy    = 0.0d0
!
#ifndef VAR_MPI
      mytid  = 0
      mparid = 0
#endif

      IF(MYTID.EQ.MPARID) THEN
        REWIND dft_temp_file_unit
        NPNTS = 0
      ENDIF
!*    Loop over batches of grid points
  100 CONTINUE
      IF(MYTID.EQ.MPARID) THEN
        IF(dft_cfg_import_grid) THEN
          READ(dft_temp_file_unit,*) NPOINT
        ELSE
          READ(dft_temp_file_unit) NPOINT
        ENDIF
      ENDIF
#if defined (VAR_MPI) 
      IF(PARCAL) CALL MPI_BCAST(NPOINT,1,my_MPI_INTEGER,MPARID, &
     &                          MPI_COMM_WORLD,IERR)
#endif
      IF( NPOINT .lt. 0 ) THEN
         GOTO 200
      ELSE IF(  NPOINT .eq. 0 ) THEN
          GOTO 100      
      ELSE
        IF(MYTID.EQ.MPARID) THEN
          NPNTS = NPNTS + NPOINT
          call alloc(kx, npoint)
          call alloc(ky, npoint)
          call alloc(kz, npoint)
          call alloc(kw, npoint)

          IF(dft_cfg_import_grid) THEN
            CALL import_dft_grid(kx, ky, kz, kw, &
     &                      dft_temp_file_unit,NPOINT)
          ELSE
            CALL REAQUA(kx, ky, kz, kw, &
     &                 dft_temp_file_unit,NPOINT)
          ENDIF
        ENDIF
#if defined (VAR_MPI) 
        IF(PARCAL) THEN
          NODESALL=NODES+1

          call alloc(idisp, nodesall)
          call alloc(isend, nodesall)

          NBATCH = NPOINT/NODESALL
          NREST = MOD(NPOINT,NODESALL)
          myNPOINT = NBATCH
!.........send nothing to master
          IF(MYTID.EQ.MPARID) THEN ! master
            ISEND(1) = NBATCH
            IDISP(1) = 0
            DO I = 1,NODES-1
              ISEND(I+1) = NBATCH
              IDISP(I+1) = NBATCH*I
            ENDDO
            ISEND(NODESALL) = NBATCH + NREST
            IDISP(NODESALL) = NBATCH*NODES
!           send data to slaves - all grid points in one send chunk
!           FIXME: scatterv makes problems with INT_STAR8...
!                  SK - Apr 2009
!
#if !defined (INT_STAR8)
!............scatterv data
             CALL MPI_SCATTERV(KX,ISEND,IDISP,MPI_REAL8, &
     &                         DUMMY,0,MPI_REAL8,MPARID, &
     &                         MPI_COMM_WORLD,IERR)
             CALL MPI_SCATTERV(KY,ISEND,IDISP,MPI_REAL8, &
     &                         DUMMY,0,MPI_REAL8,MPARID, &
     &                         MPI_COMM_WORLD,IERR)
             CALL MPI_SCATTERV(KZ,ISEND,IDISP,MPI_REAL8, &
     &                         DUMMY,0,MPI_REAL8,MPARID, &
     &                         MPI_COMM_WORLD,IERR)
             CALL MPI_SCATTERV(KW,ISEND,IDISP,MPI_REAL8, &
     &                         DUMMY,0,MPI_REAL8,MPARID, &
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
               IDISP_REST = IDISP(NODESALL) + NBATCH
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
            IF(MYTID.eq.NODES) myNPOINT = NBATCH+NREST
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
            IF(MYTID.eq.NODES) myNPOINT = NBATCH
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
              IDISP_REST = NBATCH
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
          ENDIF
!.........add rest to total number of myNPOINT on last node
          IF(MYTID.eq.NODES) myNPOINT = NBATCH+NREST
          call dealloc(idisp)
          call dealloc(isend)
        ELSE
#endif
            myNPOINT = NPOINT
#if defined (VAR_MPI) 
        ENDIF
#endif
        CALL DFTDR1(KX,Ky,KZ,KW,myNPOINT, &
     &              DMAT,DMAT_RESP_AO,DMAT_OUTER_AO, &
     &              DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP, &
     &              dft_options, &
     &              gao, gao1, gao2, gab1, ncnt, buf, &
     &              IPRINT,0)

          call dealloc(kx)
          call dealloc(ky)
          call dealloc(kz)
          call dealloc(kw)



        GO TO 100
        ENDIF
 200    CONTINUE
#if defined (VAR_MPI)
      IF(PARCAL) CALL DFTRES(EXCMAT,NXCMAT,dft_options, &
     &                       IPRINT)
#endif


      call dealloc(gao)
      if (dft_temp_need_gao1) then
        call dealloc(gao1)
      end if
      if (dft_temp_need_gao2) then
        call dealloc(gao2)
      end if
      if (dft_options%london) then
        call dealloc(gab1)
      end if
      call dealloc(ncnt)
      call dealloc(buf)


      IF(MYTID.NE.MPARID) RETURN


      orbital_nr_electrons = nelect + naelec
      error_nr_electrons   = dft_temp_nr_electrons - orbital_nr_electrons

      write(lupri, *)
      write(lupri, '(3x, a, i24)')                               &
              'number of quadrature points:  ',                  &
              npnts
      if (dft_options%energy) then
        write(lupri, '(3x, a, f24.16)')                          &
              'DFT exchange-correlation energy:               ', &
              dft_temp_xc_energy
      end if
      write(lupri, '(3x, a, f24.16)')                            &
              'number of electrons from numerical integration:', &
              dft_temp_nr_electrons
      if (abs(error_nr_electrons) > dft_cfg_accepted_error) then
        write(lupri, '(3x, a, f14.6)')                           &
              'number of electrons from orbital occupations:  ', &
               orbital_nr_electrons
        write(lupri, '(3x, a, f14.6)')                           &
              'WARNING: error in the number of electrons:     ', &
               error_nr_electrons
        write(lupri, '(12x, a, f14.6)')                          &
              'is larger than .ERRELS (set input):   ',          &
               dft_cfg_accepted_error
        write(lupri, '(12x, a)') &
              'this can happen when starting from coefficients'
        write(lupri, '(12x, a)') &
              'from a different geometry'
        write(lupri, '(12x, a)')                      &
              'or it can mean that the quadrature '// &
              'grid is inappropriate'
      end if
      write(lupri, *)

      end subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck dftdr1 */
      SUBROUTINE DFTDR1(CORX,CORY,CORZ,WEIGHT,NPOINT, &
     &                  DMAT,DMAT_RESP_AO,DMAT_OUTER_AO, &
     &                  DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP, &
     &                  dft_options, &
     &                  GAO,GAO1,GAO2,GAB1, &
     &                  NCNT,BUF,IPRINT,NYTID)
!***********************************************************************
!
!     Written by Trond Saue and T.Helgaker
!
!***********************************************************************

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, &
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "shells.h"
#include "inforb.h"
#include "orgcom.h"
#include "functionals.h"
#include "dgroup.h"
#include "densfit.h"

#ifdef MOD_AUTODIFF
#include "../autodiff/xc_fortran_interface.h"
#endif

      logical gga_qr

      DIMENSION CORX(NPOINT),CORY(NPOINT),CORZ(NPOINT), &
     &            WEIGHT(NPOINT), &
     &            DMAT(NBAST,NBAST),DMAT_RESP_AO(*),DMAT_OUTER_AO(*), &
     &            DTRMAT(*),                 GAO2(*), &
     &            GAB1(*),EXCMAT(*),  &
     &            NCNT(NBAST), BUF(NBAST,4), &
     &            VX(DERVS1),ISYMOP(NXCMAT),IHRMOP(NXCMAT)
      DIMENSION RHG(3),VT(3)
      DIMENSION derv(99)
      type(dft_options_t) :: dft_options

      real(8), target  :: gao(*)
      real(8), target  :: gao1(*)
      real(8), pointer :: gao_p(:)
      real(8), pointer :: gao1_p(:)

      real(8), allocatable :: kcfit(:)
      real(8), allocatable :: kxcvec(:)
      real(8), allocatable :: klrcfit(:)
      real(8), allocatable :: klrxcvec(:)

!
!




!
      NDIM = NBAST
#ifdef MOD_DNF
      ELCTRB = D0
      nskipped = 0
! aspg+rbast 2006-09-21
! memory allocation for variables to hold information related
! to density fitting. used in connection to the scheme of
! Birkenheuer et at (IJQC 102, 743, 2005) to compute the 
! exchange-correlation contibutions together with density fitting 
      IF (DF_FITVXC) THEN ! 
! actions disabled so far...
         dft_options%gradient = .false.
         dft_options%qr  = .false.
         dft_options%london = .false.
         call alloc(KCFIT,NBASISAUX)
         call alloc(KXCVEC,NBASISAUX)
         CALL READ_MATRIX(KXCVEC,1,NBASISAUX,'FITVXC')
         if (using_diffdens) then
            CALL READ_MATRIX(KCFIT,1,NBASISAUX,'FITDDC')
         else
            CALL READ_MATRIX(KCFIT,1,NBASISAUX,'FITCOE')
         endif
         NDIM = MAX(NBAST,NBASISAUX)
! aspg+rb, 2007-01-19
! aspg: will we need this really for linear response?
         if (dft_options%lr) then
!           DF_FITVXC_LR = .true.
         call alloc(KlrCFIT,NBASISAUX)
         call alloc(KlrXCVEC,NBASISAUX)
            CALL READ_MATRIX(KLRXCVEC,1,NBASISAUX,'LRDFVX')
            if (using_diffdens) then
               CALL READ_MATRIX(KLRCFIT,1,NBASISAUX,'LRDFCD')
            else
               CALL READ_MATRIX(KLRCFIT,1,NBASISAUX,'LRDFCE')
            endif
         write (*,*) ' '
         write (*,*) ' DEBUG: perturbed fit coefficients'
         endif
! aspg+rb, 2007-01-19
      ENDIF
!aspg+rbast, 2006-09-21
#endif /* MOD_DNF */

        nder = 0
        if (dft_temp_need_gao1) nder = 1

        gao_p  => gao(1:nbast)
        gao1_p => gao1(1:nbast*3)



!radovan: please don't delete, code under construction
!     if (dft_cfg_keep_ao) then

!       itimes = nbast
!       if (dft_temp_need_gao1) itimes = itimes*4

!       radovan: hardcoded for testing, will fix later ...
!       if (.not. associated(ao_in_mem_p)) then
!         call init_ao_in_mem(17357*itimes)
!       end if

!       if (npoint .eq. 8535) then
!         ao_in_mem_p => ao_in_mem(1:8535*itimes)
!       end if
!       if (npoint .eq. 8822) then
!         ao_in_mem_p => ao_in_mem(8535*itimes+1:17357*itimes)
!       end if

!     end if




      DO 300 IPNT = 1, NPOINT
       
        IF (IPRINT .GE. 100) THEN
          WRITE (LUPRI,'(2X,I6,4F20.12)')  &
     &       IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
        END IF
!
        WGHT = WEIGHT(IPNT)

!radovan: please don't delete, code under construction
!      if (.not. dft_cfg_ao_are_in_mem) then

         call getsos(gao, gao1, gao2, gab1, ncnt, &
     &               corx(ipnt), cory(ipnt), corz(ipnt), &
     &               buf, ndim, &
     &               nder, dft_options%london, iprint)

!radovan: please don't delete, code under construction
!        copy gao to mem
!        if (dft_cfg_keep_ao) then
!          if (dft_temp_need_gao1) then
!            call dcopy(nbast, &
!    &                  gao, 1, &
!    &                  ao_in_mem_p(4*nbast*(ipnt-1) + 1), 1)
!            call dcopy(3*nbast, &
!    &                  gao1, 1, &
!    &                  ao_in_mem_p(4*nbast*(ipnt-1) + nbast + 1), 1)
!          else
!            call dcopy(nbast, &
!    &                  gao, 1, &
!    &                  ao_in_mem_p(nbast*(ipnt-1) + 1), 1)
!          end if
!          if (ipnt .eq. 8822) then
!          dft_cfg_ao_are_in_mem = .true.
!          end if
!        end if
!      end if
!      if (dft_cfg_ao_are_in_mem) then
!        if (dft_temp_need_gao1) then
!          iao  = 4*nbast*(ipnt - 1)         + 1
!          iao1 = 4*nbast*(ipnt - 1) + nbast + 1
!          gao1_p => ao_in_mem_p(iao1:iao1 + nbast*3)
!        else
!          iao = nbast*(ipnt - 1) + 1
!        end if
!        gao_p => ao_in_mem_p(iao:iao + nbast)
!      end if












!
!       Density
!       =======
!
#ifdef MOD_DNF
        IF (DF_FITVXC) THEN  ! 
          RHO = DDOT(NBASISAUX,GAO,1,KCFIT,1)
        ELSE
#endif
          CALL GETRHO(RHO,0,GAO,DMAT,BUF)
#ifdef MOD_DNF
        ENDIF
#endif
!aspg+rbast, 2006-09-21
        IF (RHO.LE.dft_cfg_dfthr0) GO TO 300
!
!       Gradient of density
!       ===================
!

        IF (dft_temp_need_gao1) THEN
#ifdef MOD_DNF
          IF (DF_FITVXC) THEN  ! 
             DO IXYZ = 1,3
              IOFF = NBASISAUX * (IXYZ - 1) + 1
              RHG(IXYZ) =  &
     &          D1*DDOT(NBASISAUX,GAO1(IOFF),1,KCFIT,1)
! enable this later...
!             if (dft_options%lr)  &
!    &           b3(ixyz) =  &
!    &              D1*DDOT(NBASISAUX,GAO1(IOFF),1,KLRCFIT,1)
             ENDDO
          ELSE
#endif
            CALL GETGRHO(RHG,0,GAO,GAO1,DMAT,BUF)      
#ifdef MOD_DNF
          ENDIF
#endif
          RHOGRD = DSQRT(RHG(1)**2 + RHG(2)**2 + RHG(3)**2)
        END IF


#ifdef ENABLE_HESSIAN_OF_DENSITY
!
!       Hessian of density
!       ===================
!
        IF (DOLRHO) THEN
          CALL GETLRHO(RXX,RXY,RXZ,RYY,RYZ,RZZ,GAO,GAO1,GAO2, &
     &                 DMAT,BUF)       
          RHOLAP = RXX + RYY + RZZ
          RHOGHG = RHG(1)*RXX*RHG(1)+RHG(2)*RYY*RHG(2)+RHG(3)*RZZ*RHG(3) &
     &           + D2*(RHG(1)*RXY*RHG(2) + RHG(1)*RXZ*RHG(3)  &
     &               + RHG(2)*RYZ*RHG(3))
        END IF
#endif /* ifdef ENABLE_HESSIAN_OF_DENSITY */



!
!       Number of electrons
!       ===================
!
        dft_temp_nr_electrons = dft_temp_nr_electrons + WGHT*RHO




!       energy density
!       ==============

        if (dft_temp_gga_hs) Z = rhogrd*rhogrd

        if (dft_options%energy) then

          energy_density = 0.0d0

#ifdef MOD_AUTODIFF
          if (dft_cfg_autodiff) then

            if (dft_temp_gga_hs) then
              call lda_eval(LDA_SLATER, 0, rho, derv)
              energy_density = energy_density + derv(1)
              call gga_eval(GGA_BECKE, 0, rho, Z, derv)
              energy_density = energy_density + derv(1)
              call gga_eval(GGA_LYP, 0, rho, Z, derv)
              energy_density = energy_density + derv(1)
            else
              call lda_eval(LDA_SLATER, 0, rho, derv)
              energy_density = energy_density + derv(1)
              call lda_eval(LDA_VWN, 0, rho, derv)
              energy_density = energy_density + derv(1)
            end if

          else
#endif /* ifdef MOD_AUTODIFF */

            energy_density = dftenergy(rho, rhogrd)

#ifdef MOD_AUTODIFF
          end if
#endif

          dft_temp_xc_energy = dft_temp_xc_energy &
                             + wght*energy_density

        end if


!       exchange-correlation potential 
!       ==============================

        if (dft_options%energy &
     &      .or. dft_options%gradient &
     &      .or. dft_options%london) then

#ifdef MOD_AUTODIFF
          if (dft_cfg_autodiff) then

            vx = 0.0d0

            if (dft_temp_gga_hs) then
              Z = rhogrd*rhogrd

!             remember that BLYP = Slater + Becke correction + LYP
              call lda_eval(LDA_SLATER, 1, rho, derv)
              vx(fr0) = vx(fr0) &
     &                + taylor_deriv(derv, 2, (/1, 0/))
              call gga_eval(GGA_BECKE, 1, rho, Z, derv)
              vx(fr0) = vx(fr0) &
     &                + taylor_deriv(derv, 5, (/1, 0, 0, 0, 0/))
              vx(fz0) = vx(fz0) &
     &                + taylor_deriv(derv, 5, (/0, 0, 1, 0, 0/))
              call gga_eval(GGA_LYP, 1, rho, Z, derv)
              vx(fr0) = vx(fr0) &
     &                + taylor_deriv(derv, 5, (/1, 0, 0, 0, 0/))
              vx(fz0) = vx(fz0) &
     &                + taylor_deriv(derv, 5, (/0, 0, 1, 0, 0/))
            else
              call lda_eval(LDA_SLATER, 1, rho, derv)
              vx(fr0) = vx(fr0) &
     &                + taylor_deriv(derv, 2, (/1, 0/))
              call lda_eval(LDA_VWN, 1, rho, derv)
              vx(fr0) = vx(fr0) &
     &                + taylor_deriv(derv, 2, (/1, 0/))
            end if

            vx = vx*wght

          else
#endif /* ifdef MOD_AUTODIFF */

            call dftpot0(vx, wght, rho, rhogrd)

#ifdef MOD_AUTODIFF
          end if
#endif

          if (dft_temp_gga_hs) then
!           this factor has nothing to do
!           with symmetrization later
!           it comes from the fact that derivative of Z
!           yields two identical terms (one left, one right)
            vx(fz0) = 2.0d0*vx(fz0)
          end if

!         asymptotic correction
          if (dft_cfg_saop_is_active) then
            call getrho(rho_outer, 0, gao, dmat_outer_ao, buf)
            if (dft_cfg_saop_with_response_part) then
              call getrho(rho_resp, 0, gao, dmat_resp_ao, buf)
            end if
          end if
          if (dft_cfg_grac_is_active .or. dft_cfg_saop_is_active) then
            rho13 = rho**dp3
            call dftac(vx, wght, rho, rho13, &
     &                 rho_resp, rho_outer, rhogrd)
          end if

        end if !exchange-correlation potential





!
!
!       Exchange-correlation contribution to Kohn-Sham matrix
!       =====================================================
!
        IF (dft_options%energy) THEN
           IREPOP = 0
#ifdef MOD_DNF
! aspg+rb, 2007-01-19, df-vxc modularized code
           if (DF_FITVXC) then
               IF (dft_temp_gga_hs) THEN
                  DO I = 1,3
                     VT(I) = VX(FZ0)*RHG(I)
                  ENDDO
               ENDIF
               CALL fit_DOEXC(dft_temp_gga_hs, &
     &                     VX(FR0),VT,GAO,GAO1, &
     &                     KXCVEC,NBASISAUX)
           else  ! df_fitvxc
#endif
             if (dft_temp_gga_hs) then
                DO I = 1,3
                   VT(I) = VX(FZ0)*RHG(I)
                ENDDO
                CALL DOEXC2(VX(FR0),VT,IREPOP,EXCMAT,GAO,GAO1)      
             ELSE
                CALL DOEXC1(VX(FR0),IREPOP,EXCMAT,GAO)
             END IF
#ifdef MOD_DNF
           endif ! df_fitvxc
#endif
        END IF ! dft_options%energy


!       linear response
        if(dft_options%lr) then

          call dft_lr(dtrmat, &
     &                excmat, &
     &                nxcmat, &
     &                isymop, &
     &                ihrmop, &
     &                .false., &
     &                gao_p, &
     &                gao1_p, &
     &                rho, &
     &                rhogrd, &
     &                rhg, &
     &                wght, &
     &                dft_temp_gga_hs, &
     &                buf)

#ifdef MOD_CDFT
!         spin magnetization contribution
          if (dft_cfg_sdft) then
            call sdft_lr(dtrmat,excmat,nxcmat, &
     &                   isymop,ihrmop,.false., &
     &                   gao,gao1, &
     &                   rho,rhogrd,rhg,wght,dft_temp_gga_ha,buf)
          endif
!         VRG contribution
          if (dft_cfg_cdft_vrg) then
            call vrg_lr(dtrmat,excmat,nxcmat, &
     &                  isymop,ihrmop,gao,gao1, &
     &                  rho,rhogrd,rhg,wght,dft_temp_gga_hs,buf)
          endif
#endif /* ifdef MOD_CDFT */
        endif






!       quadratic response
        if(dft_options%qr) then

          gga_qr = dft_temp_gga_hs .or. dft_temp_gga_ha
          call sdft_qr(gga_qr, &
     &                 sdft, &
     &                 dtrmat,excmat, &
     &                 isymop,ihrmop,gao,gao1, &
     &                 rho,rhogrd,rhg,wght,buf)

!         the idea here is that   r^[B]               \omega_{ia}^[C]
!                                 r^[C]               \omega_{ia}^[B]
!                                 r^[BC]              \omega_{ia}
!         and the corresponding   q^[B]  \cdot \nabla \omega_{ia}^[C]
!                                 q^[C]  \cdot \nabla \omega_{ia}^[B]
!                                 q^[BC] \cdot \nabla \omega_{ia}
          
!         can be evaluated using (s)dft_lr (below)
!         only the remaining part is done in sdft_qr (above)
          
          call dft_lr(dtrmat,excmat,3, &
     &                isymop,ihrmop,.true., &
     &                gao,gao1, &
     &                rho,rhogrd,rhg,wght,dft_temp_gga_hs,buf)
          if (dft_cfg_sdft) then
            call sdft_lr(dtrmat,excmat,3, &
     &                   isymop,ihrmop,.true., &
     &                   gao,gao1, &
     &                   rho,rhogrd,rhg,wght,dft_temp_gga_ha,buf)
          endif

        endif


#ifdef MOD_AUTODIFF
#ifdef MOD_AORESPONSE
          if (dft_options%aor_xc_order > 0) then
            call aor_xc_new(dft_options%aor_xc_order,  &
     &                      nxcmat,  &
     &                      dft_temp_gga_hs, &
     &                      dtrmat,  &
     &                      excmat, &
     &                      isymop, &
     &                      gao_p, &
     &                      gao1_p, &
     &                      rho, &
     &                      rhg, &
     &                      wght, &
     &                      buf)
          end if
#endif /* ifdef MOD_AORESPONSE */
#endif /* ifdef MOD_AUTODIFF */


  300 CONTINUE !end of loop over grid points



#ifdef MOD_DNF
!aspg+rbast, 2006-09-21
! commented code below is still in use! debugging response with dfvxc
!
!     write(LUPRI,*) ' a total of ',nskipped,' out of ',npoint, &
!    &               ' integration were skipped while getting b0'
      IF (DF_FITVXC) THEN ! releasing temp space used for the df code
! aspg+rb, 2007-01-19
! aspg: do we really need this for linear response?
         if (dft_options%lr) then
!radovan: KXLRCVEC not set
          call quit('lr and fitvxc not implemented')
!         CALL SAVE_MATRIX(KXLRCVEC,1,NBASISAUX,'LRDFVX')
         call dealloc(klrcfit)
         call dealloc(klrxcvec)
         endif
! aspg+rb, 2007-01-19
!
         CALL SAVE_MATRIX(KXCVEC,1,NBASISAUX,'FITVXC')
         call dealloc(kcfit)
         call dealloc(kxcvec)
      ENDIF
!aspg+rbast, 2006-09-21
#endif /* MOD_DNF */

      end subroutine
!  /* Deck reaqua */
      SUBROUTINE REAQUA(X,Y,Z,WT,LUQUAD,NPOINT)
#include "implicit.h" 
#include "priunit.h"
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      READ(LUQUAD) X,Y,Z,WT
      end subroutine
      
#ifdef ENABLE_DFT_NUMERICAL_TEST
!  /* Deck dftnum */
      SUBROUTINE DFTNUM(DODRC,DOVWN,DOLYP,DOBCK,RHO,RHO13,RHOGRD)
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,  &
     &           DP3 = D1/D3, ADD = 1.0D-06)
      LOGICAL DODRC,DOVWN,DOLYP,DOBCK
!
      ADD2 = D2*ADD
!
      R0 = RHO
      RP = RHO + ADD
      RM = RHO - ADD
      Z0 = RHOGRD
      ZP = RHOGRD + ADD 
      ZM = RHOGRD - ADD
      R013 = R0**DP3
      RP13 = RP**DP3
      RM13 = RM**DP3
!
      IF (DODRC) THEN
         CALL VDRC(G10,R013)
         CALL V1DRC(F10,F20,R0,R013)
!
         CALL EDRC(EP,RP,RP13)
         CALL VDRC(V1P,RP13)
         CALL EDRC(EM,RM,RM13)
         CALL VDRC(V1M,RM13)
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
!
         WRITE(LUPRI,'(A,3F15.10)')' DR010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' DRC10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' DRC20 test',F20,H20,F20-H20
      END IF
      IF (DOVWN) THEN
         CALL VVWN(G10,R0,R013)
         CALL V1VWN(F10,F20,R0,R013)
!
         CALL EVWN(EP,RP,RP13)
         CALL VVWN(V1P,RP,RP13)
         CALL EVWN(EM,RM,RM13)
         CALL VVWN(V1M,RM,RM13)
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
!
         WRITE(LUPRI,'(A,3F15.10)')' VW010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' VWN10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' VWN20 test',F20,H20,F20-H20
      END IF
!
      IF (DOBCK) THEN
         CALL GBCK(G10,G01,R0,R013,Z0)
         G01 = DP5*G01*Z0 
         CALL V1BCK(F10,F01,F20,F11,F02,R0,Z0)
!
         CALL EBCK(EP,RP,RP13,Z0)
         CALL GBCK(V1P,V2P,RP,RP13,Z0)
         V2P = DP5*Z0*V2P
         CALL EBCK(EM,RM,RM13,Z0)
         CALL GBCK(V1M,V2M,RM,RM13,Z0)
         V2M = DP5*Z0*V2M
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
         H11 = (V2P - V2M)/ADD2
!
         CALL EBCK(EP,R0,R013,ZP)
         CALL GBCK(V1P,V2P,R0,R013,ZP)
         V2P = DP5*ZP*V2P
         CALL EBCK(EM,R0,R013,ZM)
         CALL GBCK(V1M,V2M,R0,R013,ZM)
         V2M = DP5*ZM*V2M
         H01 = (EP - EM)/ADD2
         H02 = (V2P - V2M)/ADD2
!
         WRITE(LUPRI,'(A,3F15.10)')' BC010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' BC001 test',G01,H01,G01-H01
         WRITE(LUPRI,'(A,3F15.10)')' BCK10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' BCK01 test',F01,H01,F01-H01
         WRITE(LUPRI,'(A,3F15.10)')' BCK20 test',F20,H20,F20-H20
         WRITE(LUPRI,'(A,3F15.10)')' BCK11 test',F11,H11,F11-H11
         WRITE(LUPRI,'(A,3F15.10)')' BCK02 test',F02,H02,F02-H02
      END IF
      IF (DOLYP) THEN
         CALL GLYP(G10,G01,R0,R013,Z0)
         G01 = DP5*G01*Z0 
         CALL V1LYP(F10,F01,F20,F11,F02,R0,Z0)
!
         CALL ELYP(EP,RP,RP13,Z0)
         CALL GLYP(V1P,V2P,RP,RP13,Z0)
         V2P = DP5*Z0*V2P
         CALL ELYP(EM,RM,RM13,Z0)
         CALL GLYP(V1M,V2M,RM,RM13,Z0)
         V2M = DP5*Z0*V2M
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
         H11 = (V2P - V2M)/ADD2
!
         CALL ELYP(EP,R0,R013,ZP)
         CALL GLYP(V1P,V2P,R0,R013,ZP)
         V2P = DP5*ZP*V2P
         CALL ELYP(EM,R0,R013,ZM)
         CALL GLYP(V1M,V2M,R0,R013,ZM)
         V2M = DP5*ZM*V2M
         H01 = (EP - EM)/ADD2
         H02 = (V2P - V2M)/ADD2
!
         WRITE(LUPRI,'(A,3F15.10)')' LY010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' LY001 test',G01,H01,G01-H01
         WRITE(LUPRI,'(A,3F15.10)')' LYP10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' LYP01 test',F01,H01,F01-H01
         WRITE(LUPRI,'(A,3F15.10)')' LYP20 test',F20,H20,F20-H20
         WRITE(LUPRI,'(A,3F15.10)')' LYP11 test',F11,H11,F11-H11
         WRITE(LUPRI,'(A,3F15.10)')' LYP02 test',F02,H02,F02-H02
      END IF
      end subroutine
#endif /* ENABLE_DFT_NUMERICAL_TEST */
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck Getavd */
      SUBROUTINE GTDOAV(DMAT,IFAC)
!***********************************************************************
!
!     Get average density from open shells; only real part !
!     Written by T.Saue March 2004
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
!
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION DMAT(N2BBASX,NZ,NFMAT)
      IF(AOC) THEN
        DO ISHELL = 1, NOPEN
          FAC = IFAC*DF(ISHELL)
          CALL DAXPY(N2BBASX,FAC,DMAT(1,1,ISHELL+1),1,DMAT,1)
        ENDDO
      ENDIF
!
      end subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE sym_ha(AMAT, &
     &                  NMAT, &
     &                  ISYMOP, &
     &                  IHRMOP)
!=======================================================================
!     symmetrize only specific symmetry blocks of a matrix
!     written by Trond Saue march 2004
!     adapted to antihermitian matrices
!-----------------------------------------------------------------------
!                                                           radovan bast
!                                           last revision: february 2007
!=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (DP5 = 0.5D0,D2 = 2.0D0)

#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AMAT(NBAST,NBAST,NZ,NMAT)
#else
      DIMENSION AMAT(NTBAS(0),NTBAS(0),NZ,NMAT)
#endif
      DIMENSION ISYMOP(NMAT), &
     &          IHRMOP(NMAT)

#include "ibtfun.h"

      DO IM = 1,NMAT
        IF(.NOT. IHRMOP(IM) .EQ. 1) THEN
          DO IXYZ = 1,3
            iq = jm4pos(4-ixyz)
            iq = iqtopq(iq,isymax(ixyz,2))
            IREPOP = IBTXOR(ISYMAX(IXYZ,2),ISYMOP(IM)-1)
            DO 110 IREP1 = 0,NBSYM-1
              IREP2 = IBTXOR(IREP1,IREPOP)
              IF(IREP2 .LT. IREP1) THEN
                DO 120 IC = 1,MC
                  NB1 = NBBAS(IREP1,IC)
                  IF(NB1 .EQ. 0) GOTO 120
                  NB2 = NBBAS(IREP2,IC)
                  IF(NB2 .EQ. 0) GOTO 120
                  JB1 = IBBAS(IREP1,IC)
                  JB2 = IBBAS(IREP2,IC)
                  DO I = 1,NB1
                    II = JB1 + I
                    DO J = 1,NB2
                      JJ = JB2 + J
                      AVERAGE = DP5*(AMAT(II,JJ,IQ,IM)  &
     &                             + AMAT(JJ,II,IQ,IM))
                      AMAT(II,JJ,IQ,IM) = AVERAGE 
                      AMAT(JJ,II,IQ,IM) = AVERAGE 
                    ENDDO
                  ENDDO
  120           CONTINUE
              ELSEIF(IREP2 .EQ. IREP1) THEN
                DO 130 IC = 1,MC
                  NB = NBBAS(IREP1,IC)
                  IF(NB .EQ. 0) GOTO 130
                  JB = IBBAS(IREP1,IC)
                  DO I = 1,NB
                    II = JB + I
                    DO J = 1,I-1
                      JJ = JB + J
                      AVERAGE = DP5*(AMAT(II,JJ,IQ,IM)  &
     &                             + AMAT(JJ,II,IQ,IM))
                      AMAT(II,JJ,IQ,IM) = AVERAGE 
                      AMAT(JJ,II,IQ,IM) = AVERAGE 
                    ENDDO
                  ENDDO
  130           CONTINUE
              ENDIF
  110       CONTINUE
          ENDDO
        ENDIF
      ENDDO

      end subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE sym_hs(AMAT, &
     &                  NMAT, &
     &                  ISYMOP, &
     &                  IHRMOP)
!=======================================================================
!     symmetrize only specific symmetry blocks of a matrix
!     written by Trond Saue march 2004
!=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (DP5 = 0.5D0,D2 = 2.0D0)

#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AMAT(NBAST,NBAST,NZ,NMAT)
#else
      DIMENSION AMAT(NTBAS(0),NTBAS(0),NZ,NMAT)
#endif
      DIMENSION ISYMOP(NMAT), &
     &          IHRMOP(NMAT)

#include "ibtfun.h"

      DO IM = 1,NMAT
        IF(.NOT. IHRMOP(IM) .EQ. -1) THEN
          IREPOP = ISYMOP(IM) - 1
          DO 110 IREP1 = 0,NBSYM-1
            IREP2 = IBTXOR(IREP1,IREPOP)
            IF(IREP2 .LT. IREP1) THEN
              DO 120 IC = 1,MC
                NB1 = NBBAS(IREP1,IC)
                IF(NB1 .EQ. 0) GOTO 120
                NB2 = NBBAS(IREP2,IC)
                IF(NB2 .EQ. 0) GOTO 120
                JB1 = IBBAS(IREP1,IC)
                JB2 = IBBAS(IREP2,IC)
                DO I = 1,NB1
                  II = JB1 + I
                  DO J = 1,NB2
                    JJ = JB2 + J
                    AVERAGE = DP5*(AMAT(II,JJ,1,IM)  &
     &                           + AMAT(JJ,II,1,IM))
                    AMAT(II,JJ,1,IM) = AVERAGE 
                    AMAT(JJ,II,1,IM) = AVERAGE 
                  ENDDO
                ENDDO
  120         CONTINUE
            ELSEIF(IREP2 .EQ. IREP1) THEN
              DO 130 IC = 1,MC
                NB = NBBAS(IREP1,IC)
                IF(NB .EQ. 0) GOTO 130
                JB = IBBAS(IREP1,IC)
                DO I = 1,NB
                  II = JB + I
                  DO J = 1,I-1
                    JJ = JB + J
                    AVERAGE = DP5*(AMAT(II,JJ,1,IM)  &
     &                           + AMAT(JJ,II,1,IM))
                    AMAT(II,JJ,1,IM) = AVERAGE 
                    AMAT(JJ,II,1,IM) = AVERAGE 
                  ENDDO
                ENDDO
  130         CONTINUE
            ENDIF
  110     CONTINUE
        ENDIF
      ENDDO

      end subroutine
#if defined (VAR_MPI)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck dftres */
      SUBROUTINE DFTRES(EXCMAT,NXCMAT,dft_options, &
     &                  IPRINT)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!     Send receive DFTresults
!     Written by T.Saue and T.helgaker
!
!***********************************************************************

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
!
#include "maxorb.h"
#include "mxcent.h"
#include "infpar.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION EXCMAT(*)
      type(dft_options_t) :: dft_options
      real(8), allocatable :: kfbuf(:)
!
      LFMAT = N2BBASX*NXCMAT*NZ
!
      if (mytid == 0) then
        call alloc(kfbuf, lfmat)
      end if

      call mpi_reduce(dft_temp_nr_electrons, &
                      sum_nr_electrons,      &
                      1,                     &
                      mpi_double_precision,  &
                      mpi_sum,               &
                      0,                     &
                      mpi_comm_world,        &
                      ierr) 

      CALL MPI_REDUCE(EXCMAT,KFBUF,LFMAT,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,0, MPI_COMM_WORLD, IERR) 

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

      IF(MYTID.EQ.0) THEN
        dft_temp_nr_electrons = sum_nr_electrons
        if (dft_options%energy) then
          dft_temp_xc_energy = sum_xc_energy
        end if
        CALL DCOPY(LFMAT,KFBUF,1,EXCMAT,1)     
        call dealloc(kfbuf)
      ENDIF
      end subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck dftprd */
!***********************************************************************
      SUBROUTINE DFTPRD
!
!     Communicate DIRAC common blocks for parallell DFT-calculation
!     Communicate arguments received from main program
!
!     Common blocks initialized : dcborb, dcbbas, dgroup, 
!                                 symmet
!
!     Based on traprd, Written by T. Saue 2000
!
!***********************************************************************
#include "implicit.h"
      PARAMETER(D0=0.0D0)
!
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "symmet.h"
#include "dgroup.h"
#include "infpar.h"
#include "parint.h"
#include "blocks.h"
!TROND
#include "shells.h"
#include "nuclei.h"
!
#if defined (VAR_MPI)

!TROND: Temporary fix !!!
      NSAUE=3*NUCIND   
      CALL MPI_BCAST(CORD,NSAUE,MPI_DOUBLE_PRECISION,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      NSAUE=3*MXSHEL
      CALL MPI_BCAST(CENT,NSAUE,MPI_DOUBLE_PRECISION,MPARID, &
     &               MPI_COMM_WORLD,IERR)
!
!     ********************************
!     ****** Send COMMON blocks ******
!     ********************************
!
!.....dcbham.h
      IF(MYTID.EQ.MPARID) THEN
        NUMELM = ICOMMSIZE(I1_HAM,I2_HAM)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_HAM,NUMELM,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
!.....dcborb.h
      IF(MYTID.EQ.MPARID) THEN
        NUMELM = ICOMMSIZE(I1_DCBORB,I2_DCBORB)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_DCBORB,NUMELM,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
!.....dcbbas.h
      IF(MYTID.EQ.MPARID) THEN
        NUMELM = ICOMMSIZE(I1_DCBBAS,I2_DCBBAS)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_DCBBAS,NUMELM,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
!.....dgroup.h
      IF(MYTID.EQ.MPARID) THEN
        NUMELM = ICOMMSIZE(I1_DGROUP,I2_DGROUP)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_DGROUP,NUMELM,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
!.....symmet.h
      IF(MYTID.EQ.MPARID) THEN
        NUMELM = ICOMMSIZE(I1_SYMMTI,I2_SYMMTI)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_SYMMTI,NUMELM,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NSYMMTR
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FMULT,NUMELM,MPI_DOUBLE_PRECISION,MPARID, &
     &               MPI_COMM_WORLD,IERR)
#endif
!
      end subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck dftnod */
      SUBROUTINE DFTNOD()
!***********************************************************************
!     Based on HERNOD 
!     written by T. Saue March 4
!
!    This is the node program for DFT                           *
!***********************************************************************


#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
!
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
!
#include "dcbbas.h"
#include "infpar.h"
!
      type(dft_options_t) :: dft_options
      integer, allocatable :: ksym(:)
      integer, allocatable :: khrm(:)
      real(8), allocatable :: kexc(:)
      real(8), allocatable :: kden(:)
      real(8), allocatable :: kdtr(:)
      real(8), allocatable :: KDMAT_OUTER_AO(:)
      real(8), allocatable :: KDMAT_RESP_AO(:)
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
      CALL DFTPAR_INIT1(dft_options,NXCMAT)
!.....Memory allocation
      LMAT = N2BBASXQ*NXCMAT
      call alloc(KEXC,LMAT)
      call alloc(ksym, nxcmat)
      call alloc(khrm, nxcmat)
      call alloc(kden,N2BBASX)
      IF(dft_options%lr.OR.dft_options%qr) THEN
        call alloc(kdtr,LMAT)
      ENDIF
!.....Allocate density matrices for SAOP
      if (dft_cfg_saop_is_active) then
        call alloc(KDMAT_OUTER_AO,N2BBASX)
        if (dft_cfg_saop_with_response_part) then
          call alloc(KDMAT_RESP_AO,N2BBASX)
        ENDIF
      ENDIF
!
!     Receive initialization data from master 2
!     =========================================
      CALL DFTPAR_INIT2(dft_options,NXCMAT, ksym, khrm, &
     &                  KDEN,KDTR, &
     &                  KDMAT_RESP_AO,KDMAT_OUTER_AO)
!.....Initialization
      IPRINT = 0
      CALL DZERO(KEXC,LMAT)
!.....Integrate
      CALL DFTDRV(KDEN,KDMAT_RESP_AO,KDMAT_OUTER_AO, &
     &            KDTR,KEXC,NXCMAT, ksym, khrm, &
     &            dft_options,IPRINT)

      call dealloc(KEXC)
      call dealloc(ksym)
      call dealloc(khrm)
      call dealloc(kden)
      IF(dft_options%lr.OR.dft_options%qr) THEN
        call dealloc(kdtr)
      ENDIF
      if (dft_cfg_saop_is_active) then
        call dealloc(KDMAT_OUTER_AO)
        if (dft_cfg_saop_with_response_part) then
          call dealloc(KDMAT_RESP_AO)
        ENDIF
      ENDIF

      end subroutine


#ifdef VAR_MPI
      subroutine dftpar_init1(dft_options, nxcmat)

!     written by T. Saue april 7 2009

      implicit none

#include "mpif.h"
#include "infpar.h"

!     --------------------------------------------------------------------------
      type(dft_options_t)                :: dft_options
      integer                            :: nxcmat
!     --------------------------------------------------------------------------
      logical                            :: i_am_the_master
      integer                            :: ierr
!     --------------------------------------------------------------------------

      if (mytid == mparid) then
        i_am_the_master = .true.
      else
        i_am_the_master = .false.
      end if

!     sync nxcmat
      call mpi_bcast(nxcmat, &
     &               1, my_MPI_INTEGER, mparid, mpi_comm_world, ierr)

!     sync various common blocks
      call dftprd()

!     sync various parameters
      call sync_dft_cfg()
      call sync_dft_temp()
      call sync_dft_options(dft_options)

!     sync functionals
      call dftsyncfunc(i_am_the_master)
 
      end subroutine
#endif /* ifdef VAR_MPI */





!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck dftinit2 */
      SUBROUTINE DFTPAR_INIT2(dft_options,NXCMAT,ISYMOP,IHRMOP, &
     &                       DMAT,DTRMAT,DMAT_RESP_AO,DMAT_OUTER_AO)
!***********************************************************************
!
!     Initialize parallel DFT run
!
!     Written by T. Saue April 7 2009
! 
!***********************************************************************


#include "implicit.h"
#include "priunit.h"
!
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "infpar.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION DMAT(N2BBASX),DTRMAT(*), &
     &          DMAT_RESP_AO(*),DMAT_OUTER_AO(*), &
     &          ISYMOP(NXCMAT),IHRMOP(NXCMAT)
      type(dft_options_t) :: dft_options
      CALL MPI_BCAST(ISYMOP,NXCMAT,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IHRMOP,NXCMAT,my_MPI_INTEGER,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DMAT,N2BBASX,MPI_DOUBLE_PRECISION,MPARID, &
     &               MPI_COMM_WORLD,IERR)
      IF(dft_options%lr.OR.dft_options%qr) THEN
        LMAT = N2BBASXQ*NXCMAT
        CALL MPI_BCAST(DTRMAT,LMAT,MPI_DOUBLE_PRECISION,MPARID, &
     &                   MPI_COMM_WORLD,IERR)
      ENDIF

      if (dft_cfg_saop_is_active) then
        CALL MPI_BCAST(DMAT_OUTER_AO,N2BBASX, &
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        if (dft_cfg_saop_with_response_part) then
!             SAOP with response part
        CALL MPI_BCAST(DMAT_RESP_AO,N2BBASX, &
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        ENDIF
      ENDIF
      end subroutine
#endif /* VAR_MPI */


end module

!radovan: the following routines cannot be part of the module
!(circular dep) 

#ifdef MOD_DNF
! /* deck fit_DOEXC */
      SUBROUTINE fit_DOEXC(ISGGA,VXC_scalar,VXC_vector,GAO,GAO1, &
     &                     VXC_integrals,NDIM)
!***********************************************************************
!
!     Construct the set of integrals of Vxc over fit fucntions,
!     <V_{xc}|\Chi_i>
!
!     Written by Andre Gomes and Radovan Bast
!     Strasbourg, January 2007
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
!
      DIMENSION GAO(*), GAO1(*), VXC_vector(3), VXC_integrals(*)
      LOGICAL ISGGA
!
! ndim should be equal to nbasisaux...
!
      IF (ISGGA) THEN
         DO I = 1, NDIM
            IOFF1 = I
            IOFF2 = IOFF1 + NDIM
            IOFF3 = IOFF2 + NDIM
! be careful with the factor below, compared to dalton and to the 
! evaluation over orbital basis 
            VXC_integrals(I) = VXC_integrals(I) &
     &                       + D1 * VXC_scalar * GAO(I) &
     &                       + D1 * VXC_vector(1) * GAO1(IOFF1) &
     &                       + D1 * VXC_vector(2) * GAO1(IOFF2) &
     &                       + D1 * VXC_vector(3) * GAO1(IOFF3)
         ENDDO
      ELSE ! only lda
         DO I = 1, NDIM 
            VXC_integrals(I) = VXC_integrals(I) &
     &                       + D1 * VXC_scalar * GAO(I)
         ENDDO
      END IF
      end subroutine
#endif /* MOD_DNF */



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck doexc1 */
      SUBROUTINE DOEXC1(VXC,IREPOP,EXCMAT,GAO)      
!***********************************************************************
!
!     Construct Kohn-Sham matrix
!
!     Written by T.Saue and T.Helgaker
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
!
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),EXCMAT(NTBAS(0),NTBAS(0))
#endif
!      
#include "ibtfun.h"
!
      DO 10 IRP1 = 0,NBSYM-1
        IRP2 = IBTXOR(IRP1,IREPOP)
        DO 20 IC = 1,MC
          NB1 = NBBAS(IRP1,IC)
          IF(NB1.EQ.0) GOTO 20
          NB2 = NBBAS(IRP2,IC)
          IF(NB2.EQ.0) GOTO 20
          JB1 = IBBAS(IRP1,IC)
          JB2 = IBBAS(IRP2,IC)
          DO J = 1, NBBAS(IRP2,IC)
            GVXC = VXC*GAO(JB2+J)
            CALL DAXPY(NBBAS(IRP1,IC),GVXC, &
     &                 GAO(JB1+1),1,EXCMAT(JB1+1,JB2+J),1)            
          ENDDO
   20   CONTINUE
   10 CONTINUE
!
      end subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck doexc2 */
      SUBROUTINE DOEXC2(SXC,VXC,IREPOP,EXCMAT,GAO,GAO1)      
!***********************************************************************
!
!     Construct the matrix
!
!     F(i,j) = GAO(I)*SXC*GAO(J)+2*GAO(I)*VXC(K)*nabla[GAO1(J)]
!
!     Upon symmetrization 
!
!     F(i,j) = (1/2)[F(i,j)+F(j,i)]
!
!     later on one obtains
!
!     E(i,j) = GAO(I)*SXC*GAO(J)+VXC(K)*nabla[GAO(I)*GAO(J)]
!
!     Written by T.Saue and T.Helgaker
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
!
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),EXCMAT(NBAST,NBAST), &
     &          GAO1(NBAST,3)
#else
      DIMENSION GAO(NTBAS(0)),EXCMAT(NTBAS(0),NTBAS(0)), &
     &          GAO1(NTBAS(0),3)
#endif
      DIMENSION VXC(3)
#include "ibtfun.h"
!.....Vector part is scaled to get symmetrization right
      DO I = 1,3
        VXC(I)=D2*VXC(I)
      ENDDO
      DO 10 IRP1 = 0,NBSYM-1
        IRP2 = IBTXOR(IRP1,IREPOP)
        DO 20 IC = 1,MC
          NB1 = NBBAS(IRP1,IC)
          IF(NB1.EQ.0) GOTO 20
          NB2 = NBBAS(IRP2,IC)
          IF(NB2.EQ.0) GOTO 20
          JB1 = IBBAS(IRP1,IC)
          JB2 = IBBAS(IRP2,IC)
          DO J = 1, NBBAS(IRP2,IC)
            FC = SXC*GAO(JB2+J)      +VXC(1)*GAO1(JB2+J,1) &
     &         + VXC(2)*GAO1(JB2+J,2)+VXC(3)*GAO1(JB2+J,3)
            CALL DAXPY(NBBAS(IRP1,IC),FC, &
     &                 GAO(JB1+1),1,EXCMAT(JB1+1,JB2+J),1)            
          ENDDO
   20   CONTINUE
   10 CONTINUE
!
      end subroutine
