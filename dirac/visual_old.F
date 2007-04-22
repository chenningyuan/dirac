C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* VISUAL_DRIVER */
      SUBROUTINE VISUAL_DRIVER(DMAT_AO,
     &                         DMAT_AO_TEMP,
     &                         GAO,
     &                         GAO1,
     &                         NCNT,
     &                         BUF,
     &                         WORK,LWORK)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "codata.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "dcbgen.h"
#include "visual.h"
#include "dcborb.h"








      CHARACTER*50 FILE_NAME
      LOGICAL   LBIT
      DIMENSION DMAT_AO(*),
     &          DMAT_AO_TEMP(*),
     &          GAO(*),
     &          GAO1(*),
     &          NCNT(*),
     &          BUF(*),
     &          WORK(LWORK),
     &          LUNIT(3,NR_QUANTITIES),
     &          QUANTITY(0:3,NR_QUANTITIES,NCUBE(3)),
     &          QUANTITY_INTEGRAL(0:3,NR_QUANTITIES),
     &          QUANTITY_INTEGRAL_LAST(0:3,NR_QUANTITIES),
     &          STEP_3D(3,0:3),
     &          ATOM_XYZ(4,NUCDEP)

#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "dgroup.h"
#include "symmet.h"
#include "sdpre.h"
#include "ibtfun.h"

#include "memint.h"







!     open files

      LUNIT_BASE = 13
      DO I = 1,3 !1D, 2D, 3D
        IF(LBIT(LBIT_VISUAL,I)) THEN
          DO IQ = 1,NR_QUANTITIES
            IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
              WRITE(FILE_NAME,'(A5,I1,A2,A)') 'PLOT.', I, 'D.',
     &                                        QUANTITY_NAME(IQ)
              OPEN(LUNIT_BASE,
     &             FILE   = FILE_NAME,
     &             STATUS = 'NEW',
     &             ACCESS = 'SEQUENTIAL',
     &             FORM   = 'FORMATTED')

              LUNIT(I,IQ) = LUNIT_BASE
              LUNIT_BASE = LUNIT_BASE + 1
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      IF(LBIT(LBIT_VISUAL,2)) THEN !2D
        DO IORIENTATION = 1,3

          NI = NCUBE(1)
          NJ = NCUBE(2)
          NK = NCUBE(3)

          IF(ORIENTATION(IORIENTATION) .EQ. 'xy') THEN
            NK = 1
          ELSEIF(ORIENTATION(IORIENTATION) .EQ. 'xz') THEN
            NJ = 1
          ELSEIF(ORIENTATION(IORIENTATION) .EQ. 'yz') THEN
            NI = 1
          ENDIF

          DO IQ = 1,NR_QUANTITIES
            IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
              WRITE(LUNIT(2,IQ),'(A)')
              WRITE(LUNIT(2,IQ),'(A)')
              WRITE(LUNIT(2,IQ),'(A2)') ORIENTATION(IORIENTATION)
              WRITE(LUNIT(2,IQ),'(A)')
            ENDIF
          ENDDO

          DO I = 1,NI
            DO J = 1,NJ
              DO K = 1,NK

                PX = -6.0D0 + (I-1)*12.0D0/(NI-1)
                PY = -6.0D0 + (J-1)*12.0D0/(NJ-1)
                PZ = -6.0D0 + (K-1)*12.0D0/(NK-1)

                IF(ORIENTATION(IORIENTATION) .EQ. 'xy') THEN
                  PZ = D0
                ELSEIF(ORIENTATION(IORIENTATION) .EQ. 'xz') THEN
                  PY = D0
                ELSEIF(ORIENTATION(IORIENTATION) .EQ. 'yz') THEN
                  PX = D0
                ENDIF

#if defined (MOD_VISUAL)
                CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                     GAO,GAO1,
     &                                     PX,PY,PZ,
     &                                     NCNT,BUF,
     &                                     QUANTITY,
     &                                     WORK,KFREE,LFREE)
#endif

!               write to file(s)

                DO IQ = 1,NR_QUANTITIES
                  IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
                    IF(QUANTITY_IS_A_VECTOR(IQ)) THEN
                      WRITE(LUNIT(2,IQ),'(7F20.10)')
     &                      PX, PY, PZ,
     &                      QUANTITY(0,IQ,1),
     &                      QUANTITY(1,IQ,1),
     &                      QUANTITY(2,IQ,1),
     &                      QUANTITY(3,IQ,1)
                    ELSE
                      WRITE(LUNIT(2,IQ),'(4F20.10)')
     &                      PX, PY, PZ,
     &                      QUANTITY(0,IQ,1)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF !2D

      IF(LBIT(LBIT_VISUAL,3)) THEN !3D

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

        NCUBE(0) = NUCDEP
        CALL DZERO(STEP_3D,12)
        DO I = 1,3
          TEMP = CORD(I,1)
          STEP_3D(I,0) = TEMP
          STEP_3D(I,I) = TEMP
        ENDDO

        I = 0
        DO J = 1,NUCIND !loop over all centers
          K = ISTBNU(J)
          DO L = 0,MAXOPR
            IF(IBTAND(L,K) .EQ. 0) THEN
              I = I + 1
              DO M = 1,3
                ATOM_XYZ(M,I) = PT(IBTAND(ISYMAX(M,1),L))*CORD(M,J)
                STEP_3D(M,0) = MIN(STEP_3D(M,0),ATOM_XYZ(M,I))
                STEP_3D(M,M) = MAX(STEP_3D(M,M),ATOM_XYZ(M,I))
              ENDDO
              ATOM_XYZ(4,I) = CHARGE(J)
            ENDIF
          ENDDO
        ENDDO

!       add some space, calculate step sizes and convert step vectors to angstrom

        DO I = 1,3
          STEP_3D(I,0) = STEP_3D(I,0) - ADD_3D
          STEP_3D(I,I) = (STEP_3D(I,I) + (ADD_3D*D2))*XTANG
          STEP_3D(I,I) = (STEP_3D(I,I) - STEP_3D(I,0))/(NCUBE(I) - 1)
        ENDDO

!       write header

        DO IQ = 1,NR_QUANTITIES
          IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
            WRITE(LUNIT(3,IQ),'(A)') 'some'
            WRITE(LUNIT(3,IQ),'(A)') 'text'
            DO I = 0,3
              WRITE(LUNIT(3,IQ),'(I5,4F12.6)')
     &                                    NCUBE(I), (STEP_3D(J,I),J=1,3)
            ENDDO
            DO I = 1,NUCDEP
              WRITE(LUNIT(3,IQ),'(I5,4F12.6)')
     &         NINT(ATOM_XYZ(4,I)), ATOM_XYZ(4,I), (ATOM_XYZ(J,I),J=1,3)
            ENDDO
          ENDIF
        ENDDO

!       grid walk

        DO ICUBE = 1,NCUBE(1)
          DO JCUBE = 1,NCUBE(2)
            DO KCUBE = 1,NCUBE(3)

              PX = STEP_3D(1,0) + (ICUBE-1)*STEP_3D(1,1)
     &                          + (JCUBE-1)*STEP_3D(1,2)
     &                          + (KCUBE-1)*STEP_3D(1,3)
              PY = STEP_3D(2,0) + (ICUBE-1)*STEP_3D(2,1)
     &                          + (JCUBE-1)*STEP_3D(2,2)
     &                          + (KCUBE-1)*STEP_3D(2,3)
              PZ = STEP_3D(3,0) + (ICUBE-1)*STEP_3D(3,1)
     &                          + (JCUBE-1)*STEP_3D(3,2)
     &                          + (KCUBE-1)*STEP_3D(3,3)

#if defined (MOD_VISUAL)
              CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                   GAO,GAO1,
     &                                   PX,PY,PZ,
     &                                   NCNT,BUF,
     &                                   QUANTITY(0,1,KCUBE),
     &                                   WORK,KFREE,LFREE)
#endif
            ENDDO

!           write fo file(s)

            DO IQ = 1,NR_QUANTITIES
              IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
                WRITE(LUNIT(3,IQ),'(6E13.5)')
     &               (QUANTITY(0,IQ,M),M=1,NCUBE(3))
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF !3D

!     close files

      DO I = 1,3 !1D, 2D, 3D
        IF(LBIT(LBIT_VISUAL,I)) THEN
          DO IQ = 1,NR_QUANTITIES
            IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
              CLOSE(LUNIT(I,IQ), STATUS = 'KEEP')
            ENDIF
          ENDDO
        ENDIF
      ENDDO


      RETURN
      END























#if defined (oldcode)
      IF(LBIT(LBIT_VISUAL,8)) THEN !integration

        IF(LBIT(LBIT_VISUAL,3)) THEN !3D

          CALL AROUND('3D numerical integration')

          LUNIT_3D_GRID = 44

          CALL DFTGRD(DMAT_AO,WORK(KFREE),LFREE,1)

          REWIND LUNIT_3D_GRID
          READ(LUNIT_3D_GRID) NR_POINTS

          CALL MEMGET('REAL',KX,NR_POINTS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KY,NR_POINTS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KZ,NR_POINTS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KW,NR_POINTS,WORK,KFREE,LFREE)

          CALL REAQUA(WORK(KX),WORK(KY),WORK(KZ),
     &                WORK(KW),
     &                LUNIT_3D_GRID,NR_POINTS)

#if defined (MOD_VISUAL)
          CALL INTEGRATE_QUANTITY(QUANTITY_INTEGRAL,
     &                            DMAT_AO,
     &                            GAO,GAO1,
     &                            WORK(KX),WORK(KY),WORK(KZ),
     &                            WORK(KW),
     &                            NR_POINTS,
     &                            NCNT,BUF,
     &                            WORK,KFREE,LFREE)
#endif

          CALL MEMREL('PLOT_SUB',WORK,KX,KX,KFREE,LFREE)

          CLOSE(LUNIT_3D_GRID,STATUS = 'DELETE')

          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A8,A29,3A20)') 'quantity',
     &                                 'scalar/norm',
     &                                 'x-component',
     &                                 'y-component',
     &                                 'z-component'
          WRITE(LUPRI,'(A)')

          DO IQ = 1,NR_QUANTITIES
            IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
              IF(QUANTITY_IS_A_VECTOR(IQ)) THEN
                WRITE(LUPRI,'(A20,4E20.10)') QUANTITY_NAME(IQ),
     &                                   QUANTITY_INTEGRAL(IQ,0),
     &                                   QUANTITY_INTEGRAL(IQ,1),
     &                                   QUANTITY_INTEGRAL(IQ,2),
     &                                   QUANTITY_INTEGRAL(IQ,3)
              ELSE
                WRITE(LUPRI,'(A20,E20.10)') QUANTITY_NAME(IQ),
     &                                  QUANTITY_INTEGRAL(IQ,0)
              ENDIF
            ENDIF
          ENDDO

          WRITE(LUPRI,'(A)')
          CALL PRSYMB(LUPRI,'-',100,0)
          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A)')

        ENDIF !3D

        IF(LBIT(LBIT_VISUAL,2)) THEN !2D

          CALL AROUND('2D Gauss-Lobatto numerical integration')

          DO IORIENTATION = 2,2
            WRITE(LUPRI,'(A)')
            WRITE(LUPRI,'(A12,A3,A17,A19,A20)') 'orientation:',
     &                                        ORIENTATION(IORIENTATION),
     &                                          'pieces',
     &                                          'order',
     &                                          'error'
            WRITE(LUPRI,'(A)')

            NR_PIECES = 5
            N = 10
            DIMI = 10.0D0
            DIMJ = 10.0D0
            LUNIT_2D_GRID = 45

            CALL DZERO(QUANTITY_INTEGRAL_LAST,NR_QUANTITIES*4)

            DO WHILE(.TRUE.)
#if defined (MOD_VISUAL)
              CALL LOBATTO_2D_GRID(ORIENTATION(IORIENTATION),
     &                             NR_PIECES,N,DIMI,DIMJ,LUNIT_2D_GRID)
#endif

              REWIND LUNIT_2D_GRID
              READ(LUNIT_2D_GRID) NR_POINTS

              CALL MEMGET('REAL',KX,NR_POINTS,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KY,NR_POINTS,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KZ,NR_POINTS,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KW,NR_POINTS,WORK,KFREE,LFREE)

              CALL REAQUA(WORK(KX),WORK(KY),WORK(KZ),
     &                    WORK(KW),
     &                    LUNIT_2D_GRID,NR_POINTS)

#if defined (MOD_VISUAL)
              CALL INTEGRATE_QUANTITY(QUANTITY_INTEGRAL,
     &                                DMAT_AO,
     &                                GAO,GAO1,
     &                                WORK(KX),WORK(KY),WORK(KZ),
     &                                WORK(KW),
     &                                NR_POINTS,
     &                                NCNT,BUF,
     &                                WORK,KFREE,LFREE)
#endif

              CALL MEMREL('PLOT_SUB',WORK,KX,KX,KFREE,LFREE)

              CLOSE(LUNIT_2D_GRID,STATUS = 'DELETE')

              ERROR = D0
              DO IQ = 1,NR_QUANTITIES
                DO IXYZ = 0,3
                  SOMETHING = QUANTITY_INTEGRAL(IQ,IXYZ)
     &                      - QUANTITY_INTEGRAL_LAST(IQ,IXYZ)
                  ERROR = ERROR + DABS(SOMETHING)
                ENDDO
              ENDDO
              WRITE(LUPRI,'(I30,I20,E30.10)') NR_PIECES, N, ERROR
              IF(ERROR .LT. 1.0D-6) GOTO 1

              N = N + 10
              CALL DZERO(QUANTITY_INTEGRAL_LAST,NR_QUANTITIES*4)
              CALL DAXPY(NR_QUANTITIES*4,D1,QUANTITY_INTEGRAL,1,
     &                                      QUANTITY_INTEGRAL_LAST,1)
            ENDDO
 1          CONTINUE

            WRITE(LUPRI,'(A)')
            WRITE(LUPRI,'(A8,A29,3A20)') 'quantity',
     &                                   'scalar/norm',
     &                                   'x-component',
     &                                   'y-component',
     &                                   'z-component'
            WRITE(LUPRI,'(A)')

            DO IQ = 1,NR_QUANTITIES
              IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
                IF(QUANTITY_IS_A_VECTOR(IQ)) THEN
                  WRITE(LUPRI,'(A20,4E20.10)') QUANTITY_NAME(IQ),
     &                                   QUANTITY_INTEGRAL(IQ,0),
     &                                   QUANTITY_INTEGRAL(IQ,1),
     &                                   QUANTITY_INTEGRAL(IQ,2),
     &                                   QUANTITY_INTEGRAL(IQ,3)
                ELSE
                  WRITE(LUPRI,'(A20,E20.10)') QUANTITY_NAME(IQ),
     &                                   QUANTITY_INTEGRAL(IQ,0)
                ENDIF
              ENDIF
            ENDDO
          ENDDO !DO IORIENTATION = 1,3
        ENDIF !2D

        WRITE(LUPRI,'(A)')
        CALL PRSYMB(LUPRI,'-',100,0)

      ENDIF !integration
#endif
#if defined (oldcode)
      IF(LBIT(LBIT_VISUAL,7)) THEN !use orbital string
        IF(ORBSTRING1(4:5) .EQ. '..') THEN
          READ(ORBSTRING1(1:3),'(I3)') ISTRING1_BEGIN
          READ(ORBSTRING1(6:8),'(I3)') ISTRING1_END
          GOTO 1
        ENDIF
        IF(ORBSTRING1(3:4) .EQ. '..') THEN
          READ(ORBSTRING1(1:2),'(I3)') ISTRING1_BEGIN
          READ(ORBSTRING1(5:8),'(I3)') ISTRING1_END
          GOTO 1
        ENDIF
        IF(ORBSTRING1(2:3) .EQ. '..') THEN
          READ(ORBSTRING1(1:1),'(I3)') ISTRING1_BEGIN
          READ(ORBSTRING1(4:8),'(I3)') ISTRING1_END
          GOTO 1
        ENDIF
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'wrong orbital string format'
        CALL QUIT('ERROR in PLOT')
 1      CONTINUE
        IF(
     &          (.NOT. (ISTRING1_BEGIN .LE. ISTRING1_END))
     &     .OR. (.NOT. (ISTRING1_BEGIN .GE. 0))
     &     .OR. (.NOT. (ISTRING1_END .LE. NESH(1)))) THEN
          WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
          WRITE(LUPRI,'(A)') 'wrong orbital string dimension'
          CALL QUIT('ERROR in PLOT')
        ENDIF
        CALL DFILL(WORK(KSTRING1),D0,NESH(1))
        DO I = ISTRING1_BEGIN,ISTRING1_END
          IF(ISTRING1_BEGIN .NE. 0) WORK(KSTRING1-1+I) = D1
        ENDDO

        IF(NFSYM .GT. 1) THEN
          IF(ORBSTRING2(4:5) .EQ. '..') THEN
            READ(ORBSTRING2(1:3),'(I3)') ISTRING2_BEGIN
            READ(ORBSTRING2(6:8),'(I3)') ISTRING2_END
            GOTO 2
          ENDIF
          IF(ORBSTRING2(3:4) .EQ. '..') THEN
            READ(ORBSTRING2(1:2),'(I3)') ISTRING2_BEGIN
            READ(ORBSTRING2(5:8),'(I3)') ISTRING2_END
            GOTO 2
          ENDIF
          IF(ORBSTRING2(2:3) .EQ. '..') THEN
            READ(ORBSTRING2(1:1),'(I3)') ISTRING2_BEGIN
            READ(ORBSTRING2(4:8),'(I3)') ISTRING2_END
            GOTO 2
          ENDIF
          WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
          WRITE(LUPRI,'(A)') 'wrong orbital string format'
          CALL QUIT('ERROR in PLOT')
 2        CONTINUE
          IF(
     &            (.NOT. (ISTRING2_BEGIN .LE. ISTRING2_END))
     &       .OR. (.NOT. (ISTRING2_BEGIN .GE. 0))
     &       .OR. (.NOT. (ISTRING2_END .LE. NESH(2)))) THEN
            WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
            WRITE(LUPRI,'(A)') 'wrong orbital string dimension'
            CALL QUIT('ERROR in PLOT')
          ENDIF
          CALL DFILL(WORK(KSTRING2),D0,NESH(2))
          DO I = ISTRING2_BEGIN,ISTRING2_END
            IF(ISTRING2_BEGIN .NE. 0) WORK(KSTRING2-1+I) = D1
          ENDDO
        ENDIF !IF(NFSYM .GT. 1) THEN
      ENDIF !IF(LBIT(LBIT_VISUAL,7)) THEN !use orbital string

#endif
