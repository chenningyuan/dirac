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

module visual

  use memory_allocator
  use visual_cfg
  use codata
  use ao_eval
  use basis_info
  use matrix_defop
  use matrix_genop, only: copy_matrix_sym
  use dirac_interface
  use visual_integrate
  use visual_in_point

  implicit none

  public run_visual

  private

#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "dcbxlr.h"



    integer, parameter :: file_unit = 137

contains

    subroutine run_visual()

!   ----------------------------------------------------------------------------
    type(matrix)         :: D_mo, D(max_nr_xvectors)
    type(matrix)         :: D_con_ao, D_con_mo, D_con(max_nr_xvectors)
    type(matrix)         :: C, Cig, Ciu, Csg, Csu

    real(8), allocatable :: eigenvalues(:)
    real(8), allocatable :: occupation(:)
    real(8), allocatable :: ao(:)
    real(8), allocatable :: buffer(:)

    real(8)              :: energy, time_start, second
    integer              :: max_ao_g_order, max_ao_m_order
    integer              :: ioff, imat, iq, iz, m, nr_dmat, ivector
    integer              :: nr_dmat_con

    logical              :: debug_me

    integer :: i, j, off_ij, off_iz, off_ij_s, off_ij_e
    type(matrix) :: D_mo_dagger, D_1a, D_1b, D_2
!   ----------------------------------------------------------------------------

    max_ao_g_order = 0
    max_ao_m_order = 0

    write(LUPRI,'(/,2x,a)') '*** Visualisation module ***'

    do ivector = 1, visual_cfg_nr_xvectors

      if (visual_cfg_iq_xvector(ivector) == iq_divm     ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_rotm     ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_divj     ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_rotj     ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_jpara    ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_divjpara ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_rotjpara ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_elf      ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_k        ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_h        ) max_ao_g_order = 1
      if (visual_cfg_iq_xvector(ivector) == iq_gradient ) max_ao_g_order = 1

      if (visual_cfg_iq_xvector(ivector) == iq_laplacian) then
!         debug this:
!         call getsos(gao,gao1,gao2,(/0.0d0/),ncnt,px,py,pz, &
!    &                buf,ntbas(0),2,.false.,0)
          call quit('debug laplacian in visual/visual.f90')
      end if

    end do


!   find occupation
!   ===============

    call alloc(eigenvalues, norbt)
    eigenvalues = 0.0d0

    call opnfil(file_unit, 'DFCOEF', 'OLD', 'visual.F90')
    call reacmo_no_work(file_unit, 'DFCOEF', (/0.0d0/), eigenvalues, (/0/), energy, 4)

    close(file_unit, status = 'keep')
    write(LUPRI,'(2x,a)') 'Visualisation : DFCOEF file read'

    call fndocc(eigenvalues)
    call dealloc(eigenvalues)


!   get orbital string
!   ==================

    call alloc(occupation, norbt)
    call setup_orbital_string(occupation)


!   get coefficients
!   ================

    call get_c(C,   i=1.0d0, s=1.0d0, g=1.0d0, u=1.0d0)
    call get_c(Cig, i=1.0d0, s=0.0d0, g=1.0d0, u=0.0d0)
    call get_c(Csg, i=0.0d0, s=1.0d0, g=1.0d0, u=0.0d0)

    if (nfsym == 2) then
      call get_c(Ciu, i=1.0d0, s=0.0d0, g=0.0d0, u=1.0d0)
      call get_c(Csu, i=0.0d0, s=1.0d0, g=0.0d0, u=1.0d0)
    end if

!   gosia debug
!    debug_me = .true.
    debug_me = .false.
    if (debug_me) then
      write(lupri, *) 'nz, nfsym, nbsym = ', nz, nfsym, nbsym
      call print_mat(C, label='C ')
      call print_mat(Cig, label='Cig ')
      call print_mat(Csg, label='Csg ')
      if (nfsym == 2) then
        call print_mat(Ciu, label='Ciu ')
        call print_mat(Csu, label='Csu ')
      end if
    end if

!   get unperturbed density matrix
!   ==============================

    call init_mat(D_mo, norbt, norbt)
    D_mo%elms = 0.0d0

    do m = 1, norbt
      D_mo%elms((m - 1)*norbt + m) = occupation(m)
    end do

    D(1) = C*(D_mo*dag(C))


!   get perturbed matrices
!   ======================

    nr_dmat = 1



    do ivector = 1, visual_cfg_nr_xvectors


      if (.not. visual_cfg_xvector_file(ivector) .eq. 'ZORDER') then
        nr_dmat = nr_dmat + 1
        vector_to_matrix(ivector) = nr_dmat

        call setup_fo_dmat_ao(ivector, occupation, D_mo)

        if (nfsym == 2) then
          if (jbtof(D_mo%irep, 1) == 2) then
!           ungerade perturbation
            D(nr_dmat) = (Cig*(D_mo*dag(Csu))) &
                       + (Ciu*(D_mo*dag(Csg))) &
                       - (Csg*(D_mo*dag(Ciu))) &
                       - (Csu*(D_mo*dag(Cig)))
          else
!           gerade perturbation
            D(nr_dmat) = (Cig*(D_mo*dag(Csg))) &
                       + (Ciu*(D_mo*dag(Csu))) &
                       - (Csg*(D_mo*dag(Cig))) &
                       - (Csu*(D_mo*dag(Ciu)))
          end if
        else
          D(nr_dmat) = (Cig*(D_mo*dag(Csg))) &
                     - (Csg*(D_mo*dag(Cig)))
        end if

        call copy_matrix_sym(D(nr_dmat), D_mo)

      else
        vector_to_matrix(ivector) = 1
      end if



    end do   ! do ivector = 1, visual_cfg_nr_xvectors



!   insert half phases
!   ==================

! gosia: in london direct i need it for D(1) also
! 
    !do imat = 2, nr_dmat
    do imat = 1, nr_dmat
      if (nz < 4) then
        do iz = 1, nz
          iq   = ipqtoq(iz, D(imat)%irep)
          ioff = ntbas(0)*ntbas(0)*(iz - 1)
          call q2bphase('D', iq, 1, D(imat)%elms(1 + ioff))
        end do
      end if
    end do



!   scale density matrices
!   ======================

    do imat = 1, nr_dmat
!gosia: here valgrind gives a warning about uninitialised value
!now D is a "temporary matrix" (function new_tmp in matrix_defop.f90)
!actually a pointer to (initialized) matrix D_mo
!to be checked...  
      D(imat) = 2.0d0*D(imat)
    end do


!   deallocate various things
!   =========================

    call dealloc(occupation)

    D_mo = 0

    C    = 0
    Cig  = 0
    Csg  = 0

    if (nfsym == 2) then
      Ciu = 0
      Csu = 0
    end if



!   initialize ao module
!   ====================

    call basis_info_read(.true., .false.)
    call ao_eval_init(max_ao_g_order, &
                      max_ao_m_order, &
                      .false.)

    call alloc(ao,     total_ao_vector_length)
    call alloc(buffer, total_ao_vector_length)


!   plot and integrate
!   ==================

    time_start = second()

      call visual_driver(nr_dmat, D, ao, buffer, 0)

!   end if
    call timtxt('>>>> time used in plot:', (second() - time_start), lupri)

    time_start = second()
    if (visual_cfg_integrate) then

        call visual_integration(nr_dmat, D, ao, buffer, 0)
    end if
    call timtxt('>>>> time used in integration:', (second() - time_start), lupri)

    call dealloc(ao)
    call dealloc(buffer)


!   deallocate density matrices
!   ===========================

    do imat = 1, nr_dmat
      D(imat) = 0
    end do


    call basis_info_delete()

    end subroutine




      subroutine visual_driver(nr_dmat, D, ao, buffer, nr_dmat_con, D_con)

      integer, intent(in) :: nr_dmat
      type(matrix)        :: D(nr_dmat)
      real(8), intent(in) :: ao(*)
      real(8), intent(in) :: buffer(*)

      integer, intent(in)      :: nr_dmat_con
      type(matrix), optional :: D_con(nr_dmat_con)

      real(8) ::          QUANTITY(0:3,visual_cfg_NCUBE(3))
      real(8) ::          STEP_3D(3,0:3)
      real(8) ::          ATOM_XYZ(4,NUCDEP)


      character(2) :: grid_orientation(1:3) = (/'xy', 'xz', 'yz'/)
      real(8) :: r(3), v(3), px, py, pz
      integer :: ICUBE, jCUBE, kCUBE, IORIENTATION, IPOINT_IN_LIST
      integer :: istep, i, j, k, l, m, ni, nj, nk, nr_steps, IXVECTOR
      integer :: lunit_2d, lunit_3d, lunit_line


!     open files

      lunit_2d   = 13
      lunit_3d   = 14
      lunit_line = 15

      if (visual_cfg_2d) then
        open(lunit_2d, &
     &       file   = 'PLOT_2D', &
     &       status = 'new', &
     &       access = 'sequential', &
     &       form   = 'formatted')
      end if

      if (visual_cfg_3d) then
        open(lunit_3d, &
     &       file   = 'PLOT_3D', &
     &       status = 'new', &
     &       access = 'sequential', &
     &       form   = 'formatted')
      end if

      if (visual_cfg_line) then
        open(lunit_line, &
     &       file   = 'PLOT_LINE', &
     &       status = 'new', &
     &       access = 'sequential', &
     &       form   = 'formatted')
      end if

      if (visual_cfg_list) then

        CALL AROUND('Requested density at point(s)')

        WRITE(LUPRI,'(A)')
        WRITE(LUPRI,'(A8,2A14,A18,A21,2A18)') 'r_x', &
     &                                        'r_y', &
     &                                        'r_z', &
     &                                        'scalar', &
     &                                        'x-component', &
     &                                        'y-component', &
     &                                        'z-component'
        WRITE(LUPRI,'(A)')

        DO IPOINT_IN_LIST = 1, visual_cfg_nr_points_in_list


          CALL GET_QUANTITY_IN_POINT(nr_dmat, D, &
     &                               AO, &
     &            visual_cfg_XYZ_LIST(IPOINT_IN_LIST,1), &
     &            visual_cfg_XYZ_LIST(IPOINT_IN_LIST,2), &
     &            visual_cfg_XYZ_LIST(IPOINT_IN_LIST,3), &
     &                               buffer, &
     &                               QUANTITY)

          WRITE(LUPRI,'(3E14.6,4E18.10)')  &
     &             visual_cfg_XYZ_LIST(IPOINT_IN_LIST,1), &
     &             visual_cfg_XYZ_LIST(IPOINT_IN_LIST,2), &
     &             visual_cfg_XYZ_LIST(IPOINT_IN_LIST,3), &
     &                                    QUANTITY(0,1), &
     &                                    QUANTITY(1,1), &
     &                                    QUANTITY(2,1), &
     &                                    QUANTITY(3,1)

        ENDDO

        WRITE(LUPRI,'(A)')
        CALL PRSYMB(LUPRI,'-',100,0)

        WRITE(LUPRI,'(A)')
        WRITE(LUPRI,'(A)')

      ENDIF

      if (visual_cfg_line) then

        v = visual_cfg_line_to - visual_cfg_line_from
        v = v/visual_cfg_line_nr_steps

        do istep = 0, visual_cfg_line_nr_steps

          r = visual_cfg_line_from + istep*v

          call get_quantity_in_point(nr_dmat, D, &
     &                               ao, &
     &                               r(1), r(2), r(3), &
     &                               buffer, &
     &                               quantity)

          write(lunit_line, '(7e18.10)') r(1),  &
     &                                   r(2), &
     &                                   r(3), &
     &                                   quantity(0, 1), &
     &                                   quantity(1, 1), &
     &                                   quantity(2, 1), &
     &                                   quantity(3, 1)

        end do
      end if

      if (visual_cfg_2d) then

        NR_STEPS = 1 + visual_cfg_ADD_2D*2.0d0/visual_cfg_STEP_2D

        DO IORIENTATION = 1,3

          NI = NR_STEPS
          NJ = NR_STEPS
          NK = NR_STEPS

          IF(IORIENTATION .EQ. 1) NK = 1   ! 'xy'
          IF(IORIENTATION .EQ. 2) NJ = 1   ! 'xz'
          IF(IORIENTATION .EQ. 3) NI = 1   ! 'yz'

          WRITE(LUNIT_2D,'(A)')
          WRITE(LUNIT_2D,'(A)')  &
     &             'plane: ' // grid_orientation(IORIENTATION)

          DO I = 1,NI
            DO J = 1,NJ
              DO K = 1,NK

        PX = visual_cfg_ORIGIN_2D(1) - visual_cfg_ADD_2D  &
     &          + (I-1)*visual_cfg_STEP_2D
        PY = visual_cfg_ORIGIN_2D(2) - visual_cfg_ADD_2D &
     &          + (j-1)*visual_cfg_STEP_2D
        PZ = visual_cfg_ORIGIN_2D(3) - visual_cfg_ADD_2D &
     &          + (k-1)*visual_cfg_STEP_2D

                IF(IORIENTATION .EQ. 1) &                        ! 'xy'
     &                     PZ = visual_cfg_ORIGIN_2D(3)
                IF(IORIENTATION .EQ. 2) &                        ! 'xz'
     &                     PY = visual_cfg_ORIGIN_2D(2)
                IF(IORIENTATION .EQ. 3) &                        ! 'yz'
     &                     PX = visual_cfg_ORIGIN_2D(1)


                  CALL GET_QUANTITY_IN_POINT(nr_dmat, D, &
     &                                       ao, &
     &                                       PX,PY,PZ, &
     &                                       buffer, &
     &                                       QUANTITY)


!               write to file

                WRITE(LUNIT_2D,'(7E18.10)') PX, PY, PZ, &
     &                                      QUANTITY(0,1), &
     &                                      QUANTITY(1,1), &
     &                                      QUANTITY(2,1), &
     &                                      QUANTITY(3,1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF !2D

      if (visual_cfg_3d) then

!       calculate cube parameters

!       NCUBE(0)      number of atoms
!       NCUBE(1)      number of points in the "x"-direction
!       NCUBE(2)      number of points in the "y"-direction
!       NCUBE(3)      number of points in the "z"-direction

!       STEP_3D(*,0)  coordinates of initial point
!       STEP_3D(*,1)  step vector for the "x"-direction
!       STEP_3D(*,2)  step vector for the "y"-direction
!       STEP_3D(*,3)  step vector for the "z"-direction

!       ATOM_XYZ(1,I) x-coordinate of atom I
!       ATOM_XYZ(2,I) y-coordinate of atom I
!       ATOM_XYZ(3,I) z-coordinate of atom I
!       ATOM_XYZ(4,I) charge of atom I

        visual_cfg_NCUBE(0) = NUCDEP
        CALL DZERO(STEP_3D,12)
        DO I = 1,3
          STEP_3D(I,0) = CORD(I,1)
          STEP_3D(I,I) = CORD(I,1)
        ENDDO

        I = 0
        DO J = 1,NUCIND !loop over all centers
          K = ISTBNU(J)
          DO L = 0,MAXOPR
            IF(iand(L,K) .EQ. 0) THEN
              I = I + 1
              DO M = 1,3
                ATOM_XYZ(M,I) = PT(iand(ISYMAX(M,1),L))*CORD(M,J)

                STEP_3D(M,0) = MIN(STEP_3D(M,0),ATOM_XYZ(M,I))
                STEP_3D(M,M) = MAX(STEP_3D(M,M),ATOM_XYZ(M,I))
              ENDDO
              ATOM_XYZ(4,I) = CHARGE(J)
            ENDIF
          ENDDO
        ENDDO

!       add some space, calculate step sizes and convert step vectors to angstrom

        DO I = 1,3
          STEP_3D(I,0) = STEP_3D(I,0) - visual_cfg_ADD_3D
          STEP_3D(I,I) = (STEP_3D(I,I) + (visual_cfg_ADD_3D*2.0d0))*XTANG
      STEP_3D(I,I) = (STEP_3D(I,I) - STEP_3D(I,0)) &
     &        /(visual_cfg_NCUBE(I) - 1)
        ENDDO

!       write header

        WRITE(LUNIT_3D,'(A)') 'some'
        WRITE(LUNIT_3D,'(A)') 'text'
        DO I = 0,3
          WRITE(LUNIT_3D,'(I5,4F12.6)') visual_cfg_NCUBE(I),  &
     &               (STEP_3D(J,I),J=1,3)
        ENDDO
        DO I = 1,NUCDEP
          WRITE(LUNIT_3D,'(I5,4F12.6)') &
     &         NINT(ATOM_XYZ(4,I)), ATOM_XYZ(4,I), (ATOM_XYZ(J,I),J=1,3)
        ENDDO

!       grid walk

        DO ICUBE = 1,visual_cfg_NCUBE(1)
          DO JCUBE = 1,visual_cfg_NCUBE(2)
            DO KCUBE = 1,visual_cfg_NCUBE(3)

              PX = STEP_3D(1,0) + (ICUBE-1)*STEP_3D(1,1) &
     &                          + (JCUBE-1)*STEP_3D(1,2) &
     &                          + (KCUBE-1)*STEP_3D(1,3)
              PY = STEP_3D(2,0) + (ICUBE-1)*STEP_3D(2,1) &
     &                          + (JCUBE-1)*STEP_3D(2,2) &
     &                          + (KCUBE-1)*STEP_3D(2,3)
              PZ = STEP_3D(3,0) + (ICUBE-1)*STEP_3D(3,1) &
     &                          + (JCUBE-1)*STEP_3D(3,2) &
     &                          + (KCUBE-1)*STEP_3D(3,3)

!             write even grid so that it can be fed into XC integrator
!             for subsequent experiments
!             zero weight is so that it does not mess up results
!             write(lunit_3d, '(4d20.12)') px, py, pz, 0.0d0


                CALL GET_QUANTITY_IN_POINT(nr_dmat, D, &
     &                                     ao, &
     &                                     PX,PY,PZ, &
     &                                     buffer, &
     &                                     QUANTITY(0,KCUBE))


            ENDDO

!           write to file

            !write(*, *) 'xyzq', px, py, pz, 
            WRITE(LUNIT_3D,'(6E13.5)')  &
     &         (QUANTITY(0,M),M=1, visual_cfg_NCUBE(3))
          ENDDO
        ENDDO
      ENDIF !3D


!     close files

      if (visual_cfg_2d)   close(lunit_2d,   status = 'keep')
      if (visual_cfg_3d)   close(lunit_3d,   status = 'keep')
      if (visual_cfg_line) close(lunit_line, status = 'keep')

      end subroutine


      subroutine setup_orbital_string(occupation)

      real(8)              :: occupation(norbt)

      integer, allocatable :: kvec(:)
      integer              :: i, ifsym, nvec, ioff, ishift

      occupation = 0.0d0

      ioff = 0

      do ifsym = 1, nfsym

        ioff = ioff + npsh(ifsym)

        if (visual_cfg_use_orbital_string) then

          call alloc(kvec, norb(ifsym))
          kvec = 0
          nvec = 1

          call numlst(visual_cfg_orbital_string(ifsym), &
     &                kvec, &
     &                norb(ifsym), &
     &               -npsh(ifsym), &
     &                nesh(ifsym), &
     &                ifsym, &
     &                nvec)

          do i = 1, nvec
            ishift = kvec(i)
            if (ishift < 0) ishift = ishift + 1
            occupation(ioff + ishift) = 1.0d0
          end do

          call dealloc(kvec)

        else

          do i = 1, nish(ifsym)
            occupation(ioff + i) = 1.0d0
          end do

        end if

        ioff = ioff + nesh(ifsym)

      end do

      end subroutine

      SUBROUTINE SETUP_FO_DMAT_AO(IXVECTOR, &
     &                            STRING, D_mo)

      CHARACTER LAB1*16,&
     &          LAB2*16,&
     &          TYP1*2,&
     &          TYP2*2,&
     &          LBUF*16,&
     &          PBUF*1
      LOGICAL   FILE_EXISTS
      type(matrix) :: D_mo
      real(8) :: STRING(*)

      real(8), allocatable :: bvec(:)

      integer, allocatable :: orbital_rotation_indices(:, :)
      real(8), allocatable :: response_vector(:, :)
      integer              :: length
      character(2)         :: vector_type
      real(8)              :: ih_factor

      integer              :: nr_vectors, LUNIT_XVECTOR
      integer              :: ixvector, NBAS, NORBS, len
      integer              :: JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2, IVECTOR
      real(8)              :: FREQ1, FREQ2, ERGRSP, unused

!gosia, i'll remove it later
      integer :: i, j
      character(10) :: pattern
      character(16) :: text   

      INQUIRE(FILE = visual_cfg_xvector_file(IXVECTOR),EXIST = FILE_EXISTS)
      IF(.NOT. FILE_EXISTS) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A,A6,A)') 'solution vector file ', &
     &    visual_cfg_XVECTOR_FILE(IXVECTOR),' not found'
        CALL QUIT('ERROR in VISUAL')
      ENDIF
      LUNIT_XVECTOR = 17
      OPEN(LUNIT_XVECTOR, &
     &     FILE = visual_cfg_xvector_file(IXVECTOR), &
     &     STATUS = 'OLD', &
     &     ACCESS = 'SEQUENTIAL', &
     &     FORM = 'UNFORMATTED')

!     get number of vectors on file

      NR_VECTORS = 0
      REWIND LUNIT_XVECTOR
 1    CONTINUE
      READ(LUNIT_XVECTOR) LAB1
      IF(LAB1 .EQ. 'END_OF_THIS_FILE') THEN
        GOTO 2
      ELSE
        NR_VECTORS = NR_VECTORS + 1
!       skip solution vector
        READ(LUNIT_XVECTOR)
!       skip orbital rotations
        READ(LUNIT_XVECTOR)
        GOTO 1
      ENDIF
 2    CONTINUE

!     errors

      IF(visual_cfg_NR_xVECTORS .EQ. 0) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'NR_VECTORS .EQ. 0'
        CALL QUIT('ERROR in VISUAL')
      ENDIF

      IF(visual_cfg_IWHICH_XVECTOR(IXVECTOR) .GT. NR_VECTORS) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'IWHICH_XVECTOR(IXVECTOR) .GT. NR_VECTORS'
        CALL QUIT('ERROR in VISUAL')
      ENDIF

      IF(visual_cfg_IWHICH_XVECTOR(IXVECTOR) .LT. 1) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'IWHICH_XVECTOR(IXVECTOR) .LT. 1'
        CALL QUIT('ERROR in VISUAL')
      ENDIF

!     move to the right record

      IVECTOR = 0
      REWIND LUNIT_XVECTOR
 3    CONTINUE
      READ(LUNIT_XVECTOR)
      IVECTOR = IVECTOR + 1
      IF(IVECTOR .EQ. visual_cfg_IWHICH_XVECTOR(IXVECTOR)) THEN
        BACKSPACE(LUNIT_XVECTOR)
        GOTO 4
      ELSE
        READ(LUNIT_XVECTOR)
        READ(LUNIT_XVECTOR)
        GOTO 3
      ENDIF
 4    CONTINUE


!     get type and dimensions of record IWHICH_XVECTOR(IXVECTOR)

      READ(LUNIT_XVECTOR) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2, &
     &                    JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2, &
     &                    UNUSED, LEN, INTFLG, ERGRSP, NBAS, NORBS



      visual_cfg_irep_xvector(IXVECTOR) = JSYMOP1 - 1
      length  = LEN/NZ
      D_mo%irep   = JSYMOP1 - 1


    select case (typ1(1:1))
      case ('E')
        vector_type = 'pp'
      case ('P')
        vector_type = 'pn'
      case default
        call quit('unexpected typ1(1:1)')
    end select

    select case (typ1(2:2))
      case ('+')
        D_mo%ih_sym =  1
        ih_factor   =  1.0d0
      case ('-')
        D_mo%ih_sym = -1
        ih_factor   = -1.0d0
      case default
        call quit('unexpected typ1(2:2)')
    end select


!   allocate and read response vector
    call alloc(response_vector, length, nz)
    call readt(lunit_xvector, length*nz, response_vector)
!   allocate and read orbital rotation array
    call alloc(orbital_rotation_indices, 2, length)
    call readi(lunit_xvector, length*2, orbital_rotation_indices)
!   close file of solution vectors
    close(lunit_xvector, status = 'keep')

!fixme: this is a good place to zero elements according to occupation

    D_mo%elms = 0.0d0
    call scatter_vector(length,                   &
                        orbital_rotation_indices, &
                        ih_factor,                &
                        response_vector,          &
                        D_mo%elms,                &
                        D_mo%irep)


    call dealloc(orbital_rotation_indices)
    call dealloc(response_vector)

  end subroutine



  subroutine scatter_vector(length,                   &
                            orbital_rotation_indices, &
                            h,                        &
                            response_vector,          &
                            matrix,                   &
                            irep)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: length
    integer, intent(in)  :: orbital_rotation_indices(2, length)
    real(8), intent(in)  :: h
    real(8), intent(in)  :: response_vector(length, nz)
    real(8), intent(out) :: matrix(norbt, norbt, nz)
    integer, intent(in)  :: irep
!   ----------------------------------------------------------------------------
    integer              :: i, s, is, iz
    real(8)              :: f
!   ----------------------------------------------------------------------------

    do is = 1, length

      i = orbital_rotation_indices(1, is)
      s = orbital_rotation_indices(2, is)

      do iz = 1, nz

        if (ipqtoq(iz, irep) > 1) then
          f = -1.0
        else
          f =  1.0
        end if

!              row
!              ¦
!              ¦  column
!              ¦  ¦
        matrix(s, i, iz) = matrix(s, i, iz) &
                         +     response_vector(is, iz)
        matrix(i, s, iz) = matrix(i, s, iz) &
                         - f*h*response_vector(is, iz)
      end do
    end do

  end subroutine

end module
