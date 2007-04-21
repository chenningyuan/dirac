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
     &          QUANTITY(NR_QUANTITIES,0:3,NCUBE(3)),
     &          QUANTITY_INTEGRAL(NR_QUANTITIES,0:3),
     &          QUANTITY_INTEGRAL_LAST(NR_QUANTITIES,0:3),
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
          DO IQUANTITY = 1,NR_QUANTITIES
            IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
              WRITE(FILE_NAME,'(A5,I1,A2,A)') 'PLOT.', I, 'D.',
     &                                        QUANTITY_NAME(IQUANTITY)
              OPEN(LUNIT_BASE,
     &             FILE   = FILE_NAME,
     &             STATUS = 'NEW',
     &             ACCESS = 'SEQUENTIAL',
     &             FORM   = 'FORMATTED')

              LUNIT(I,IQUANTITY) = LUNIT_BASE
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

          DO IQUANTITY = 1,NR_QUANTITIES
            IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
              WRITE(LUNIT(2,IQUANTITY),'(A)')
              WRITE(LUNIT(2,IQUANTITY),'(A)')
              WRITE(LUNIT(2,IQUANTITY),'(A2)') ORIENTATION(IORIENTATION)
              WRITE(LUNIT(2,IQUANTITY),'(A)')
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

                CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                     GAO,GAO1,
     &                                     PX,PY,PZ,
     &                                     NCNT,BUF,
     &                                     QUANTITY,
     &                                     WORK,KFREE,LFREE)

!               write to file(s)

                DO IQUANTITY = 1,NR_QUANTITIES
                  IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
                    IF(QUANTITY_IS_A_VECTOR(IQUANTITY)) THEN
                      WRITE(LUNIT(2,IQUANTITY),'(7F20.10)')
     &                      PX, PY, PZ,
     &                      QUANTITY(IQUANTITY,0,1),
     &                      QUANTITY(IQUANTITY,1,1),
     &                      QUANTITY(IQUANTITY,2,1),
     &                      QUANTITY(IQUANTITY,3,1)
                    ELSE
                      WRITE(LUNIT(2,IQUANTITY),'(4F20.10)')
     &                      PX, PY, PZ,
     &                      QUANTITY(IQUANTITY,0,1)
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

        DO IQUANTITY = 1,NR_QUANTITIES
          IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
            WRITE(LUNIT(3,IQUANTITY),'(A)') 'some'
            WRITE(LUNIT(3,IQUANTITY),'(A)') 'text'
            DO I = 0,3
              WRITE(LUNIT(3,IQUANTITY),'(I5,4F12.6)')
     &                                    NCUBE(I), (STEP_3D(J,I),J=1,3)
            ENDDO
            DO I = 1,NUCDEP
              WRITE(LUNIT(3,IQUANTITY),'(I5,4F12.6)')
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

              CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                   GAO,GAO1,
     &                                   PX,PY,PZ,
     &                                   NCNT,BUF,
     &                                   QUANTITY(1,0,KCUBE),
     &                                   WORK,KFREE,LFREE)
            ENDDO

!           write fo file(s)

            DO IQUANTITY = 1,NR_QUANTITIES
              IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
                WRITE(LUNIT(3,IQUANTITY),'(6E13.5)')
     &               (QUANTITY(IQUANTITY,0,M),M=1,NCUBE(3))
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF !3D

!     close files

      DO I = 1,3 !1D, 2D, 3D
        IF(LBIT(LBIT_VISUAL,I)) THEN
          DO IQUANTITY = 1,NR_QUANTITIES
            IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
              CLOSE(LUNIT(I,IQUANTITY), STATUS = 'KEEP')
            ENDIF
          ENDDO
        ENDIF
      ENDDO


      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_QUANTITY_IN_POINT */
      SUBROUTINE GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                 GAO,GAO1,
     &                                 PX,PY,PZ,
     &                                 NCNT,BUF,
     &                                 QUANTITY,
     &                                 WORK,KFREE,LFREE)
C=======================================================================
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (DP12 = 1.0D0/2.0D0)
      PARAMETER (DP23 = 2.0D0/3.0D0)
      PARAMETER (DP53 = 5.0D0/3.0D0)

#include "visual.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbxrs.h"

      LOGICAL   LBIT
      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          NCNT(*),
     &          BUF(*),
     &          QUANTITY(NR_QUANTITIES,0:3),
     &          VECTOR(3),
     &          WORK(*)


      IPOINT_N        =  1 !    number density
      IPOINT_M        =  2 !    magnetization
      IPOINT_DIVM     =  3 !div(magnetization)
      IPOINT_ROTM     =  4 !rot(magnetization)
      IPOINT_RXJ      =  5 !    r \times j
      IPOINT_J        =  6 !    current density
      IPOINT_DIVJ     =  7 !div(current density)
      IPOINT_ROTJ     =  8 !rot(current density)
      IPOINT_JPARA    =  9 !    paramagnetic current density
      IPOINT_DIVJPARA = 10 !div(paramagnetic current density)
      IPOINT_ROTJPARA = 11 !rot(paramagnetic current density)
      IPOINT_JDIA     = 12 !    diamagnetic current density
      IPOINT_V        = 13 !    velocity density
      IPOINT_DIVV     = 14 !div(velocity density)
      IPOINT_ROTV     = 15 !rot(velocity density) a.k.a. vorticity
      IPOINT_ELF      = 16 !    ELF
      IPOINT_K        = 17 !    kinetic energy density
      IPOINT_H        = 18 !    helicity density
      IPOINT_GAMMA5   = 19 !    \gamma_5 density


      CALL DZERO(QUANTITY,NR_QUANTITIES*4)

      CALL GETSOS(GAO,GAO1,UNUSED,UNUSED,NCNT,PX,PY,PZ,
     &            BUF,WORK(KFREE),LFREE,NTBAS(0),1,.FALSE.,0)

!     density

      IQUANTITY = IPOINT_N        !    number density
      CALL GETRHO(DENSITY,IREP_VISUAL,GAO,DMAT_AO,BUF)
      QUANTITY(IQUANTITY,0) = DENSITY

!     ELF

      IQUANTITY = IPOINT_ELF      !    ELF
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN

        CALL GETGRHO(VECTOR,IREP_VISUAL,GAO,GAO1,DMAT_AO,BUF)
        GRADIENT_SQUARED = VECTOR(1)**D2 + VECTOR(2)**D2 + VECTOR(3)**D2

        ELFKIN = D0
        DO I = 1,3
          DO IC = 1,MC
            TEMP = D0
            CALL GETRH1(TEMP,IC,IREP_VISUAL,GAO1((I-1)*NTBAS(0)+1),
     &                  DMAT_AO,BUF)
            ELFKIN = ELFKIN + TEMP
          ENDDO
        ENDDO
        ELFKIN = ELFKIN/D2

        PI =  3.1415926536D0
        FACTOR = 0.3D0*((D3*PI*PI)**DP23)

        H = FACTOR*(DENSITY**DP53)
        G = ELFKIN - GRADIENT_SQUARED/(8.0D0*DENSITY)

        ELFCHI = G/H

        QUANTITY(IQUANTITY,0) = D1/(D1 + ELFCHI*ELFCHI)
      ENDIF

!     diamagnetic current density

      IQUANTITY = IPOINT_JDIA     !    diamagnetic current density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
        CALL GETRHO(DENSITY,0,GAO,DMAT_AO,BUF)

        QUANTITY(IQUANTITY,1) = -(D1/D2)*PY*DENSITY
        QUANTITY(IQUANTITY,2) =  (D1/D2)*PX*DENSITY
        QUANTITY(IQUANTITY,3) =   D0
      ENDIF

!     paramagnetic current density

      IQUANTITY = IPOINT_JPARA    !    paramagnetic current density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_JPARA(.TRUE.,
     &                 VECTOR,
     &                 IREP_VISUAL,
     &                 GAO,
     &                 GAO1,
     &                 DMAT_AO,
     &                 BUF)

        QUANTITY(IQUANTITY,1) = VECTOR(1)
        QUANTITY(IQUANTITY,2) = VECTOR(2)
        QUANTITY(IQUANTITY,3) = VECTOR(3)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     magnetization

      IQUANTITY = IPOINT_M        !    magnetization
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_M(.TRUE.,
     &             VECTOR,
     &             IREP_VISUAL,
     &             GAO,
     &             UNUSED,
     &             DMAT_AO,
     &             BUF)

        QUANTITY(IQUANTITY,1) = VECTOR(1)
        QUANTITY(IQUANTITY,2) = VECTOR(2)
        QUANTITY(IQUANTITY,3) = VECTOR(3)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     current density

      IQUANTITY = IPOINT_J        !    current density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_J(.TRUE.,
     &             VECTOR,
     &             IREP_VISUAL,
     &             GAO,
     &             UNUSED,
     &             DMAT_AO,
     &             BUF)

        QUANTITY(IQUANTITY,1) = VECTOR(1)
        QUANTITY(IQUANTITY,2) = VECTOR(2)
        QUANTITY(IQUANTITY,3) = VECTOR(3)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     r \times j

      IQUANTITY = IPOINT_RXJ      !    r \times j
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_J(.TRUE.,
     &             VECTOR,
     &             IREP_VISUAL,
     &             GAO,
     &             UNUSED,
     &             DMAT_AO,
     &             BUF)

        QUANTITY(IQUANTITY,1) = PY*VECTOR(3) - PZ*VECTOR(2)
        QUANTITY(IQUANTITY,2) = PZ*VECTOR(1) - PX*VECTOR(3)
        QUANTITY(IQUANTITY,3) = PX*VECTOR(2) - PY*VECTOR(1)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     divergence of magnetization

      IQUANTITY = IPOINT_DIVM     !div(magnetization)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_DIVM(.TRUE.,
     &                QUANTITY(IQUANTITY,0),
     &                IREP_VISUAL,
     &                GAO,
     &                GAO1,
     &                DMAT_AO,
     &                BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     divergence of current density

      IQUANTITY = IPOINT_DIVJ     !div(current density)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_DIVJ(.TRUE.,
     &                QUANTITY(IQUANTITY,0),
     &                IREP_VISUAL,
     &                GAO,
     &                GAO1,
     &                DMAT_AO,
     &                BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     divergence of paramagnetic current density

      IQUANTITY = IPOINT_DIVJPARA !div(paramagnetic current density)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
        CALL QUIT('not implemented')
      ENDIF

!     curl of magnetization

      IQUANTITY = IPOINT_ROTM     !rot(magnetization)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_ROTM(.TRUE.,
     &                VECTOR,
     &                IREP_VISUAL,
     &                GAO,
     &                GAO1,
     &                DMAT_AO,
     &                BUF)

        QUANTITY(IQUANTITY,1) = VECTOR(1)
        QUANTITY(IQUANTITY,2) = VECTOR(2)
        QUANTITY(IQUANTITY,3) = VECTOR(3)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     curl of current density

      IQUANTITY = IPOINT_ROTJ     !rot(current density)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_ROTJ(.TRUE.,
     &                VECTOR,
     &                IREP_VISUAL,
     &                GAO,
     &                GAO1,
     &                DMAT_AO,
     &                BUF)

        QUANTITY(IQUANTITY,1) = VECTOR(1)
        QUANTITY(IQUANTITY,2) = VECTOR(2)
        QUANTITY(IQUANTITY,3) = VECTOR(3)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     curl of paramagnetic current density

      IQUANTITY = IPOINT_ROTJPARA !rot(paramagnetic current density)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_ROTJPARA(.TRUE.,
     &                    VECTOR,
     &                    IREP_VISUAL,
     &                    GAO,
     &                    GAO1,
     &                    DMAT_AO,
     &                    BUF)

        QUANTITY(IQUANTITY,1) = VECTOR(1)
        QUANTITY(IQUANTITY,2) = VECTOR(2)
        QUANTITY(IQUANTITY,3) = VECTOR(3)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     \gamma_5 density

      IQUANTITY = IPOINT_GAMMA5   !    \gamma_5 density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)

        CALL GET_GAMMA5(.TRUE.,
     &                  QUANTITY(IQUANTITY,0),
     &                  IREP_VISUAL,
     &                  GAO,
     &                  UNUSED,
     &                  DMAT_AO,
     &                  BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     kinetic energy density

      IQUANTITY = IPOINT_K        !    kinetic energy density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
        CALL QUIT('not implemented')
      ENDIF

!     helicity

      IQUANTITY = IPOINT_H        !    helicity density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
        CALL QUIT('not implemented')
      ENDIF

!     IQUANTITY = IPOINT_V        !    velocity density
!     IQUANTITY = IPOINT_DIVV     !div(velocity density)
!     IQUANTITY = IPOINT_ROTV     !rot(velocity density) a.k.a. vorticity

!     get norm for quantities that are vectors

      DO IQUANTITY = 1,NR_QUANTITIES
        IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
          IF(QUANTITY_IS_A_VECTOR(IQUANTITY)) THEN
            QUANTITY(IQUANTITY,0) = QUANTITY(IQUANTITY,1)**D2
     &                            + QUANTITY(IQUANTITY,2)**D2
     &                            + QUANTITY(IQUANTITY,3)**D2

            QUANTITY(IQUANTITY,0) = DSQRT(QUANTITY(IQUANTITY,0))
          ENDIF
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
#endif
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

          DO IQUANTITY = 1,NR_QUANTITIES
            IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
              IF(QUANTITY_IS_A_VECTOR(IQUANTITY)) THEN
                WRITE(LUPRI,'(A20,4E20.10)') QUANTITY_NAME(IQUANTITY),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,0),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,1),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,2),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,3)
              ELSE
                WRITE(LUPRI,'(A20,E20.10)') QUANTITY_NAME(IQUANTITY),
     &                                  QUANTITY_INTEGRAL(IQUANTITY,0)
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
              DO IQUANTITY = 1,NR_QUANTITIES
                DO IXYZ = 0,3
                  SOMETHING = QUANTITY_INTEGRAL(IQUANTITY,IXYZ)
     &                      - QUANTITY_INTEGRAL_LAST(IQUANTITY,IXYZ)
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

            DO IQUANTITY = 1,NR_QUANTITIES
              IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
                IF(QUANTITY_IS_A_VECTOR(IQUANTITY)) THEN
                  WRITE(LUPRI,'(A20,4E20.10)') QUANTITY_NAME(IQUANTITY),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,0),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,1),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,2),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,3)
                ELSE
                  WRITE(LUPRI,'(A20,E20.10)') QUANTITY_NAME(IQUANTITY),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,0)
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
