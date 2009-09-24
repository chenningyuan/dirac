module visual

  use memory_allocator
  use visual_cfg
  use codata
  use dft_ao
  use matrix_defop
  use dirac_interface

  implicit none

  public run_visual
  public get_quantity_in_point

  private

#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"

      integer, parameter :: IQ_N           =  1 !                  \psi^\dag        \psi = n
      integer, parameter :: IQ_M           =  2 !                  \psi^\dag \Sigma \psi
      integer, parameter :: IQ_DIVM        =  3 !    \nabla \cdot  \psi^\dag \Sigma \psi
      integer, parameter :: IQ_ROTM        =  4 !    \nabla \times \psi^\dag \Sigma \psi
      integer, parameter :: IQ_J           =  5 !-ec               \psi^\dag \alpha \psi = j
      integer, parameter :: IQ_DIVJ        =  6 !-ec \nabla \cdot  \psi^\dag \alpha \psi
      integer, parameter :: IQ_ROTJ        =  7 !-ec \nabla \times \psi^\dag \alpha \psi
      integer, parameter :: IQ_JPARA       =  8 !                  j_para
      integer, parameter :: IQ_DIVJPARA    =  9 !    \nabla \cdot  j_para
      integer, parameter :: IQ_ROTJPARA    = 10 !    \nabla \times j_para
      integer, parameter :: IQ_JDIA        = 11 !                  j_dia (hardcoded for B along z)
      integer, parameter :: IQ_ELF         = 12 !                  ELF
      integer, parameter :: IQ_K           = 13 !                c \psi^\dag (\alpha \cdot p) \psi
      integer, parameter :: IQ_H           = 14 !                  \psi^\dag (\Sigma \cdot p) \psi
      integer, parameter :: IQ_GAMMA5      = 15 !                  \psi^\dag  \gamma^5        \psi

!     WARNING: the factor 1/c^2 is missing on purpose for compatibility with the
!              operators defined in prp/pamprp.F
      integer, parameter :: IQ_NUCDIPX     = 16 !      -1/(r_K^3) (r_K \times j)_x
      integer, parameter :: IQ_NUCDIPY     = 17 !      -1/(r_K^3) (r_K \times j)_y
      integer, parameter :: IQ_NUCDIPZ     = 18 !      -1/(r_K^3) (r_K \times j)_z

      integer, parameter :: IQ_EDIPX       = 19 !                  r_x n
      integer, parameter :: IQ_EDIPY       = 20 !                  r_y n
      integer, parameter :: IQ_EDIPZ       = 21 !                  r_z n
      integer, parameter :: IQ_MDIPX       = 22 !            -1/2 (r_gauge \times j)_x
      integer, parameter :: IQ_MDIPY       = 23 !            -1/2 (r_gauge \times j)_y
      integer, parameter :: IQ_MDIPZ       = 24 !            -1/2 (r_gauge \times j)_z

      integer, parameter :: IQ_GRADIENT    = 25 !
      integer, parameter :: IQ_LAPLACIAN   = 26 !



      integer, save      :: vector_to_matrix(max_nr_xvectors) = 0

contains

      SUBROUTINE run_VISUAL()

      real(8), allocatable :: occupation(:)

      real(8), allocatable :: cmo(:)
      real(8), allocatable :: eig(:)

      real(8), allocatable :: ao(:)
      real(8), allocatable :: buffer(:)

      integer :: IXVECTOR, lucoef
      real(8) :: toterg, time_start, second
      integer :: max_ao_g_order, max_ao_m_order, i, nr_dmat, io

      type(matrix)         :: D(max_nr_xvectors), D_mo, C

      max_ao_g_order = 0
      max_ao_m_order = 0

      do ixvector = 1, visual_cfg_nr_xvectors
        if (visual_cfg_iq_xvector(ixvector) == iq_divm     ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_rotm     ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_divj     ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_rotj     ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_jpara    ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_divjpara ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_rotjpara ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_elf      ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_k        ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_h        ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_gradient ) max_ao_g_order = 1
        if (visual_cfg_iq_xvector(ixvector) == iq_laplacian) then
!         debug this:
!         call getsos(gao,gao1,gao2,(/0.0d0/),ncnt,px,py,pz, &
!    &                buf,ntbas(0),2,.false.,0)
          call quit('debug laplacian in visual/visual.f90')
        end if
      end do


      call alloc(cmo      ,n2bbasxq)
      call alloc(eig      ,norbt)

!     coefficients

      CALL DZERO(CMO,N2BBASXQ)
      CALL DZERO(EIG,NORBT)

      lucoef = 137
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','VISUAL')
      CALL REACMO_no_work(LUCOEF,'DFCOEF',CMO,EIG,idummy, &
     &            toterg, 6)

      CLOSE(LUCOEF,STATUS = 'KEEP')

      call fndocc(eig)
      call dealloc(eig)

!     orbital string

      call alloc(occupation, norbt)
      call setup_orbital_string(occupation)

!     get coefficients

      call get_C(C, i=1.0d0, s=1.0d0, g=1.0d0, u=1.0d0)

!     get unperturbed MO density matrix

      call init_mat(D_mo, norbt, norbt)
      D_mo%elms = 0.0d0
      do io = 1, norbt
        D_mo%elms((io-1)*norbt + io) = occupation(io)
      end do

!     transform to AO

      D(1) = C*(D_mo*dag(C))
      D(1) = 2.0d0*D(1)
      D_mo = 0
      C    = 0

!     get perturbed matrices
      nr_dmat = 1
      do ixvector = 1, visual_cfg_nr_xvectors
        if (.not. visual_cfg_xvector_file(ixvector) .eq. 'ZORDER') then
          nr_dmat = nr_dmat + 1
          D(nr_dmat) = 1.0d-100*D(1)
          call setup_fo_dmat_ao(cmo, ixvector, D(nr_dmat)%elms, occupation)
          vector_to_matrix(ixvector) = nr_dmat
        else
          vector_to_matrix(ixvector) = 1
        end if
      end do

      call dealloc(occupation)
      call dealloc(cmo)

!   initialize AO module
    call dft_ao_init(.false.,        &
                     max_ao_g_order, &
                     max_ao_m_order, &
                     .false.,        &
                     1.0d-15)

    call alloc(ao,     total_ao_vector_length)
    call alloc(buffer, total_ao_vector_length)

    time_start = second()
!   if (visual_cfg_integrate) then
      call visual_driver(nr_dmat, D, ao, buffer)
!   end if
    call timtxt('>>>> time used in plot:', (second() - time_start), lupri)

    time_start = second()
    if (visual_cfg_integrate) then
      call visual_integration(nr_dmat, D, ao, buffer)
    end if
    call timtxt('>>>> time used in integration:', (second() - time_start), lupri)

    call dealloc(ao)
    call dealloc(buffer)

!   deallocate density matrices
    do i = 1, nr_dmat
      D(i) = 0
    end do

    end subroutine


      subroutine visual_driver(nr_dmat, D, ao, buffer)

      integer, intent(in) :: nr_dmat
      type(matrix)        :: D(nr_dmat)
      real(8), intent(in) :: ao(*)
      real(8), intent(in) :: buffer(*)

      real(8) ::          QUANTITY(0:3,visual_cfg_NCUBE(3))
      real(8) ::          STEP_3D(3,0:3)
      real(8) ::          ATOM_XYZ(4,NUCDEP)

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

          IF(visual_cfg_ORIENTATION(IORIENTATION) .EQ. 'xy') NK = 1
          IF(visual_cfg_ORIENTATION(IORIENTATION) .EQ. 'xz') NJ = 1
          IF(visual_cfg_ORIENTATION(IORIENTATION) .EQ. 'yz') NI = 1

          WRITE(LUNIT_2D,'(A)')
          WRITE(LUNIT_2D,'(A)')  &
     &             'plane: ' // visual_cfg_ORIENTATION(IORIENTATION)

          DO I = 1,NI
            DO J = 1,NJ
              DO K = 1,NK

        PX = visual_cfg_ORIGIN_2D(1) - visual_cfg_ADD_2D  &
     &          + (I-1)*visual_cfg_STEP_2D
        PY = visual_cfg_ORIGIN_2D(2) - visual_cfg_ADD_2D &
     &          + (j-1)*visual_cfg_STEP_2D
        PZ = visual_cfg_ORIGIN_2D(3) - visual_cfg_ADD_2D &
     &          + (k-1)*visual_cfg_STEP_2D

                IF(visual_cfg_ORIENTATION(IORIENTATION) .EQ. 'xy') &
     &                     PZ = visual_cfg_ORIGIN_2D(3)
                IF(visual_cfg_ORIENTATION(IORIENTATION) .EQ. 'xz') &
     &                     PY = visual_cfg_ORIGIN_2D(2)
                IF(visual_cfg_ORIENTATION(IORIENTATION) .EQ. 'yz') &
     &                     PX = visual_cfg_ORIGIN_2D(1)

                CALL GET_QUANTITY_IN_POINT(nr_dmat, D, &
     &                                     ao, &
     &                                     PX,PY,PZ, &
     &                                     buffer, &
     &                                     QUANTITY)

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

              CALL GET_QUANTITY_IN_POINT(nr_dmat, D, &
     &                                   ao, &
     &                                   PX,PY,PZ, &
     &                                   buffer, &
     &                                   QUANTITY(0,KCUBE))
            ENDDO

!           write to file

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
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* GET_QUANTITY_IN_POINT */
      SUBROUTINE GET_QUANTITY_IN_POINT(nr_dmat, D, &
     &                                 ao, &
     &                                 PX,PY,PZ, &
     &                                 buffer, &
     &                                 QUANTITY_ADDMEUP)

      integer, intent(in) :: nr_dmat
      type(matrix)        :: D(nr_dmat)
      real(8) ::          AO(*)
      real(8) ::          buffer(*)
      real(8) ::          QUANTITY(0:3)
      real(8) ::          QUANTITY_ADDMEUP(0:3)
      real(8) ::          VECTOR_A(3)
      real(8) ::          VECTOR_B(3)
      real(8) :: factor, fermi_kin, elf_chi, GRADIENT_SQUARED, density
      real(8) :: RSQ
      real(8) :: rxx, rxy, rxz, ryy, ryz, rzz
      real(8) :: px, py, pz, ddot, scalar
      real(8) :: rx, ry, rz, TF_KIN, BOSON_KIN
      integer :: ic, iq, IXVECTOR, ixyz, imat

      call get_so(px, py, pz, ao, buffer)

      CALL DZERO(QUANTITY_ADDMEUP,4)

      DO IXVECTOR = 1, visual_cfg_NR_XVECTORS
        imat = vector_to_matrix(ixvector)

!       IF(visual_cfg_XVECTOR_FILE(IXVECTOR) .EQ. 'ZORDER') &

        CALL DZERO(QUANTITY,4)

        IQ = IQ_N
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          DENSITY = 0.0d0

          DO IC = 1,MC
            IF(IC .EQ. 1) THEN
              IF(.NOT. visual_cfg_IOPTION_XVECTOR(IXVECTOR) .EQ. 2) THEN
                CALL GETRH1(DENSITY,1,visual_cfg_irep_xvector(IXVECTOR), &
     &                      AO,D(imat)%elms,buffer)
              ENDIF
            ENDIF
            IF(IC .EQ. 2) THEN
              IF(.NOT. visual_cfg_IOPTION_XVECTOR(IXVECTOR) .EQ. 1) THEN
                CALL GETRH1(DENSITY,2,visual_cfg_irep_xvector(IXVECTOR), &
     &                      AO,D(imat)%elms,buffer)
              ENDIF
            ENDIF
          ENDDO

          CALL DAXPY(1,1.0d0,DENSITY,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_GRADIENT
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GETGRHO(VECTOR_A, &
     &                 visual_cfg_irep_xvector(IXVECTOR), &
     &                 AO, &
     &                 ao(ao_start(1, 0):ao_end(1, 0)), &
     &                 D(imat)%elms, &
     &                 buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_LAPLACIAN
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
!         CALL GETLRHO(RXX, &
!    &                 RXY, &
!    &                 RXZ, &
!    &                 RYY, &
!    &                 RYZ, &
!    &                 RZZ, &
!    &                 GAO, &
!    &                 GAO1, &
!    &                 GAO2, &
!    &                 D(imat)%elms, &
!    &                 buffer)

!         CALL DAXPY(1,1.0d0,RXX + RYY + RZZ,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_ELF
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          FERMI_KIN = 0.0d0
          DO IXYZ = 1,3
            DO IC = 1,MC
              SCALAR = 0.0d0

              CALL GETRH1(SCALAR, &
     &                    IC, &
     &                    visual_cfg_irep_xvector(IXVECTOR), &
     &                    AO(ao_start(1, 0)+(IXYZ-1)*NTBAS(0)), &
     &                    D(imat)%elms, &
     &                    buffer)

              FERMI_KIN = FERMI_KIN + SCALAR
            ENDDO
          ENDDO

          FERMI_KIN = FERMI_KIN/2.0d0

          CALL GETRHO(DENSITY, &
     &                visual_cfg_irep_xvector(IXVECTOR), &
     &                AO, &
     &                D(imat)%elms, &
     &                buffer)

          CALL GETGRHO(VECTOR_A, &
     &                 visual_cfg_irep_xvector(IXVECTOR), &
     &                 AO, &
     &                 ao(ao_start(1, 0):ao_end(1, 0)), &
     &                 D(imat)%elms, &
     &                 buffer)

          GRADIENT_SQUARED = DDOT(3,VECTOR_A,1,VECTOR_A,1)

          BOSON_KIN = GRADIENT_SQUARED/(8.0D0*DENSITY)

          TF_KIN = 0.3D0*((3.0d0*PI*PI)**(2.0d0/3.0d0))*(DENSITY**(5.0d0/3.0d0))

          ELF_CHI = (FERMI_KIN - BOSON_KIN)/TF_KIN

          SCALAR = 1.0d0/(1.0d0 + ELF_CHI*ELF_CHI)

          CALL DAXPY(1,1.0d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_JDIA
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GETRHO(DENSITY, &
     &                visual_cfg_irep_xvector(IXVECTOR), &
     &                AO, &
     &                D(imat)%elms, &
     &                buffer)

          VECTOR_A(1) = -0.5d0*PY*DENSITY
          VECTOR_A(2) =  0.5d0*PX*DENSITY
          VECTOR_A(3) =  0.0d0

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_M
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_M(.TRUE., &
     &               VECTOR_A, &
     &               visual_cfg_irep_xvector(IXVECTOR), &
     &               AO, &
     &               D(imat)%elms, &
     &               buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_J
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE., &
     &               VECTOR_A, &
     &               visual_cfg_irep_xvector(IXVECTOR), &
     &               AO, &
     &               D(imat)%elms, &
     &               buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_JPARA
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_JPARA(.TRUE., &
     &                   VECTOR_A, &
     &                   visual_cfg_irep_xvector(IXVECTOR), &
     &                   AO, &
     &                   ao(ao_start(1, 0):ao_end(1, 0)), &
     &                   D(imat)%elms, &
     &                   buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_DIVM
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVM(.TRUE., &
     &                  SCALAR, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(1,1.0d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_DIVJ
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVJ(.TRUE., &
     &                  SCALAR, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(1,1.0d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_H
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVM(.TRUE., &
     &                  SCALAR, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(1,0.25d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_K
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVJ(.TRUE., &
     &                  SCALAR, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(1,0.5d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_ROTM
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_ROTM(.TRUE., &
     &                  VECTOR_A, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_ROTJ
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_ROTJ(.TRUE., &
     &                  VECTOR_A, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_ROTJPARA
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_ROTJPARA(.TRUE., &
     &                      VECTOR_A, &
     &                      visual_cfg_irep_xvector(IXVECTOR), &
     &                      AO, &
     &                      ao(ao_start(1, 0):ao_end(1, 0)), &
     &                      D(imat)%elms, &
     &                      buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_GAMMA5
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_GAMMA5(.TRUE., &
     &                    SCALAR, &
     &                    visual_cfg_irep_xvector(IXVECTOR), &
     &                    AO, &
     &                    D(imat)%elms, &
     &                    buffer)

          CALL DAXPY(1,1.0d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IF((visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPX) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPY) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPZ)) THEN

          CALL GET_J(.TRUE., &
     &               VECTOR_A, &
     &               visual_cfg_irep_xvector(IXVECTOR), &
     &               AO, &
     &               D(imat)%elms, &
     &               buffer)

          RX = PX - CORD(1,visual_cfg_IOPTION_XVECTOR(IXVECTOR))
          RY = PY - CORD(2,visual_cfg_IOPTION_XVECTOR(IXVECTOR))
          RZ = PZ - CORD(3,visual_cfg_IOPTION_XVECTOR(IXVECTOR))

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          RSQ = RX*RX + RY*RY + RZ*RZ
          FACTOR = -1.0d0/(RSQ**1.5d0)

          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPX) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(1),1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPY) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(2),1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPZ) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(3),1,QUANTITY(0),1)
        ENDIF

        IF((visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPX) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPY) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPZ)) THEN

          CALL GETRHO(DENSITY, &
     &                visual_cfg_irep_xvector(IXVECTOR), &
     &                AO, &
     &                D(imat)%elms, &
     &                buffer)

          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPX) &
     &      CALL DAXPY(1,1.0d0,DENSITY*PX,1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPY) &
     &      CALL DAXPY(1,1.0d0,DENSITY*PY,1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPZ) &
     &      CALL DAXPY(1,1.0d0,DENSITY*PZ,1,QUANTITY(0),1)
        ENDIF

        IF((visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPX) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPY) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPZ)) THEN

          CALL GET_J(.TRUE., &
     &               VECTOR_A, &
     &               visual_cfg_irep_xvector(IXVECTOR), &
     &               AO, &
     &               D(imat)%elms, &
     &               buffer)

          RX = PX - visual_cfg_GAUGE_ORIGIN(1)
          RY = PY - visual_cfg_GAUGE_ORIGIN(2)
          RZ = PZ - visual_cfg_GAUGE_ORIGIN(3)

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          FACTOR = -0.5d0

          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPX) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(1),1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPY) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(2),1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPZ) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(3),1,QUANTITY(0),1)
        ENDIF

        CALL DAXPY(4,visual_cfg_COMMON_FACTOR, &
     &       QUANTITY(0),1,QUANTITY_ADDMEUP(0),1)
      ENDDO

      RETURN
      END subroutine

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

end module
