C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* VISUAL_INPUT */
      SUBROUTINE VISUAL_INPUT(WORD,RESET,WORK,LWORK)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"

      PARAMETER (NTABLE = 60)

#include "dcbgen.h"
#include "visual.h"
#include "dgroup.h"

      LOGICAL   SET,
     &          RESET,
     &          LBIT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD_SAVE*7
      DIMENSION WORK(LWORK)

      SAVE SET

!                   1         2         3         4
      DATA TABLE /'.XXXXXX','.2D    ','.3D    ','.QUANTI',
!                   5         6         7         8
     &            '.XXXXXX','.XXXXXX','.STRING','.INTEGR',
!                   9        10        11        12
     &            '.3DADD ','.CUBE  ','.COMMON','.FREQ+H',
!                  13        14        15        16
     &            '.FREQ-H','.ORIGIN','.XXXXXX','.XXXXXX',
!                  17        18        19        20
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  21        22        23        24
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  25        26        27        28
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  29        30        31        32
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  33        34        35        36
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  37        38        39        40
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  41        42        43        44
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  45        46        47        48
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  49        50        51        52
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  53        54        55        56
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  57        58        59        60
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/

      DATA SET /.FALSE./

      IF(SET) THEN
        IF(RESET) SET = .FALSE.
        RETURN
      ENDIF
      IF(RESET) THEN
        SET = .FALSE.
      ELSE
        SET = .TRUE.
      ENDIF

!     initialize visual.h

      LBIT_VISUAL = 0

      NCUBE(1) = 80
      NCUBE(2) = 80
      NCUBE(3) = 80

      VISUAL_ORIGIN(1) = 0.0D0
      VISUAL_ORIGIN(2) = 0.0D0
      VISUAL_ORIGIN(3) = 0.0D0

      ADD_3D        = 4.0D0
      COMMON_FACTOR = 1.0D0

      DO IXVECTOR = 1,MAX_NR_XVECTORS
        IREP_XVECTOR(IXVECTOR)    = 0
        IWHICH_XVECTOR(IXVECTOR)  = 1
        IOPTION_XVECTOR(IXVECTOR) = 0
      ENDDO

      NR_XVECTORS = 0

      ORIENTATION(1) = 'xy'
      ORIENTATION(2) = 'xz'
      ORIENTATION(3) = 'yz'

!     process input

      IF(WORD .EQ. '*VISUAL') THEN
        WORD_SAVE = WORD
  100   CONTINUE
        READ(LUCMD,'(A7)') WORD
        PROMPT = WORD(1:1)
        IF(PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
          GOTO 100
        ELSEIF(PROMPT .EQ. '.') THEN
          DO 200 I = 1,NTABLE
            IF(TABLE(I) .EQ. WORD) THEN
            GOTO ( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     &            11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     &            21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
     &            31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
     &            41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
     &            51, 52, 53, 54, 55, 56, 57, 58, 59, 60), I
            ENDIF
  200     CONTINUE
          IF(WORD .EQ. '.OPTION') THEN
            CALL PRTAB(NTABLE,TABLE,WORD_SAVE//' input keywords',LUPRI)
            GOTO 100
          ENDIF
          WRITE(LUPRI,'(/,3A,/)') ' Keyword "', WORD,
     &                            '" not recognized in VISUAL_INPUT.'
          CALL PRTAB(NTABLE,TABLE,WORD_SAVE//' input keywords',LUPRI)
          CALL QUIT('Illegal keyword in VISUAL_INPUT.')

    1     CONTINUE
          GOTO 100

    2     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,2)
          GOTO 100

    3     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,3)
          GOTO 100

    4     CONTINUE
          READ(LUCMD,*) NR_XVECTORS
          DO IXVECTOR = 1,NR_XVECTORS
            READ(LUCMD,*) XVECTOR_FILE(IXVECTOR),
     &                    IWHICH_XVECTOR(IXVECTOR),
     &                    IQ_XVECTOR(IXVECTOR),
     &                    IOPTION_XVECTOR(IXVECTOR)
          ENDDO
          GOTO 100

    5     CONTINUE
    6     CONTINUE
          GOTO 100

    7     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,7)
          READ(LUCMD,'(A)') ORBITAL_STRING(1)
          IF(NFSYM .GT. 1) READ(LUCMD,'(A)') ORBITAL_STRING(2)
          GOTO 100

    8     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,8)
          GOTO 100

    9     CONTINUE
          READ(LUCMD,*) ADD_3D
          GOTO 100

   10     CONTINUE
          READ(LUCMD,*) NCUBE(1), NCUBE(2), NCUBE(3)
          GOTO 100

   11     CONTINUE
          READ(LUCMD,*) COMMON_FACTOR
          GOTO 100

   12     CONTINUE
   13     CONTINUE
          GOTO 100

   14     CONTINUE
          READ(LUCMD,*) VISUAL_ORIGIN(1),
     &                  VISUAL_ORIGIN(2),
     &                  VISUAL_ORIGIN(3)
          GOTO 100

   15     CONTINUE
   16     CONTINUE
   17     CONTINUE
   18     CONTINUE
   19     CONTINUE
   20     CONTINUE
   21     CONTINUE
   22     CONTINUE
   23     CONTINUE
   24     CONTINUE
   25     CONTINUE
   26     CONTINUE
   27     CONTINUE
   28     CONTINUE
   29     CONTINUE
   30     CONTINUE
   31     CONTINUE
   32     CONTINUE
   33     CONTINUE
   34     CONTINUE
   35     CONTINUE
   36     CONTINUE
   37     CONTINUE
   38     CONTINUE
   39     CONTINUE
   40     CONTINUE
   41     CONTINUE
   42     CONTINUE
   43     CONTINUE
   44     CONTINUE
   45     CONTINUE
   46     CONTINUE
   47     CONTINUE
   48     CONTINUE
   49     CONTINUE
   50     CONTINUE
   51     CONTINUE
   52     CONTINUE
   53     CONTINUE
   54     CONTINUE
   55     CONTINUE
   56     CONTINUE
   57     CONTINUE
   58     CONTINUE
   59     CONTINUE
   60     CONTINUE
          GOTO 100

        ELSEIF(PROMPT .EQ. '*') THEN
          GOTO 300
        ELSE
          WRITE(LUPRI,'(/,3A,/)') ' Keyword "', WORD,
     &                            '" not recognized in VISUAL_INPUT.'
          CALL PRTAB(NTABLE,TABLE,WORD_SAVE//' input keywords',LUPRI)
          CALL QUIT('Illegal keyword in VISUAL_INPUT.')
        ENDIF
      ENDIF
  300 CONTINUE

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* VISUAL */
      SUBROUTINE VISUAL(WORK,LWORK)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "visual.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dgroup.h"

      LOGICAL   LBIT
      DIMENSION WORK(LWORK)

#include "memint.h"

!     CALL VISUAL_CHECK_CONSISTENCY()
!     CALL VISUAL_REPORT()

!     allocate memory

      CALL MEMGET('REAL',KCMO      ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG      ,NORBT,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_AO  ,N2BBASXQ*(MAX_NR_XVECTORS+1),
     &            WORK,KFREE,LFREE)
!     CALL MEMGET('REAL',KSMAT     ,N2BBASX,
!    &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSTRING   ,NESH(1) + NESH(2),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO      ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO1     ,NTBAS(0)*3,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCNT     ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF      ,NTBAS(0)*4,
     &            WORK,KFREE,LFREE)

!     coefficients

      CALL DZERO(WORK(KCMO),N2BBASXQ)
      CALL DZERO(WORK(KEIG),NORBT)

      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','VISUAL')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),UNUSED,
     &            UNUSED,6,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS = 'KEEP')

      CALL FNDOCC(WORK(KEIG),WORK,KFREE,LFREE)

!     orbital string

      CALL SETUP_ORBITAL_STRING(WORK(KSTRING))

!     density matrices

      IOFF = 0
      DO IXVECTOR = 1,NR_XVECTORS
        IF(.NOT. XVECTOR_FILE(IXVECTOR) .EQ. 'ZORDER') THEN
          CALL SETUP_FO_DMAT_AO(WORK(KCMO),
     &                          IXVECTOR,
     &                          WORK(KDMAT_AO+IOFF),
     &                          WORK(KSTRING),
     &                          WORK(KFREE),LFREE)
        ENDIF
        IOFF = IOFF + N2BBASXQ
      ENDDO
      CALL SETUP_ZO_DMAT_AO(WORK(KCMO),
     &                      WORK(KDMAT_AO+IOFF),
     &                      WORK(KSTRING),
     &                      WORK(KFREE),LFREE)

!     CALL GTOVLX(WORK(KSMAT),SSMTRC,WORK(KFREE),LFREE)
!     addme = DDOT(N2BBASX,WORK(KDMAT_AO),1,WORK(KSMAT),1)
!     write(lupri,*) 'debug', addme
!     addme = DDOT(N2BBASX,WORK(KDMAT_AO+IOFF),1,WORK(KSMAT),1)
!     write(lupri,*) 'debug', addme

      CALL VISUAL_DRIVER(WORK(KDMAT_AO),
     &                   WORK(KGAO),
     &                   WORK(KGAO1),
     &                   WORK(KNCNT),
     &                   WORK(KBUF),
     &                   WORK(KFREE),LFREE)

      IF(LBIT(LBIT_VISUAL,8)) THEN !integration
        CALL VISUAL_INTEGRATION(WORK(KDMAT_AO),
     &                          WORK(KGAO),
     &                          WORK(KGAO1),
     &                          WORK(KNCNT),
     &                          WORK(KBUF),
     &                          WORK(KFREE),LFREE)
      ENDIF

!     release memory

      CALL MEMREL('VISUAL',WORK,KCMO,KCMO,KFREE,LFREE)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* VISUAL_DRIVER */
      SUBROUTINE VISUAL_DRIVER(DMAT_AO,
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
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "dgroup.h"
#include "symmet.h"
#include "sdpre.h"
#include "ibtfun.h"

      LOGICAL   LBIT
      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          NCNT(*),
     &          BUF(*),
     &          WORK(LWORK),
     &          QUANTITY(0:3,NCUBE(3)),
     &          STEP_3D(3,0:3),
     &          ATOM_XYZ(4,NUCDEP)

#include "memint.h"

!     open files

      LUNIT_2D = 13
      LUNIT_3D = 14

      IF(LBIT(LBIT_VISUAL,2)) THEN !2D
        OPEN(LUNIT_2D,
     &       FILE   = 'PLOT_2D',
     &       STATUS = 'NEW',
     &       ACCESS = 'SEQUENTIAL',
     &       FORM   = 'FORMATTED')
      ENDIF
      IF(LBIT(LBIT_VISUAL,3)) THEN !3D
        OPEN(LUNIT_3D,
     &       FILE   = 'PLOT_3D',
     &       STATUS = 'NEW',
     &       ACCESS = 'SEQUENTIAL',
     &       FORM   = 'FORMATTED')
      ENDIF

      IF(LBIT(LBIT_VISUAL,2)) THEN !2D
        DO IORIENTATION = 1,3

          NI = NCUBE(1)
          NJ = NCUBE(2)
          NK = NCUBE(3)

          IF(ORIENTATION(IORIENTATION) .EQ. 'xy') NK = 1
          IF(ORIENTATION(IORIENTATION) .EQ. 'xz') NJ = 1
          IF(ORIENTATION(IORIENTATION) .EQ. 'yz') NI = 1

          WRITE(LUNIT_2D,'(A)')
          WRITE(LUNIT_2D,'(A)')
          WRITE(LUNIT_2D,'(A2)') ORIENTATION(IORIENTATION)
          WRITE(LUNIT_2D,'(A)')

          RANGE_MIN = -6.0D0
          RANGE_MAX = -RANGE_MIN

          DO I = 1,NI
            DO J = 1,NJ
              DO K = 1,NK

                PX = VISUAL_ORIGIN(1)
     &             + RANGE_MIN + (I-1)*(RANGE_MAX - RANGE_MIN)/(NI-1)
                PY = VISUAL_ORIGIN(2)
     &             + RANGE_MIN + (J-1)*(RANGE_MAX - RANGE_MIN)/(NJ-1)
                PZ = VISUAL_ORIGIN(3)
     &             + RANGE_MIN + (K-1)*(RANGE_MAX - RANGE_MIN)/(NK-1)

                IF(ORIENTATION(IORIENTATION) .EQ. 'xy')
     &                                        PZ = VISUAL_ORIGIN(3) + D0
                IF(ORIENTATION(IORIENTATION) .EQ. 'xz')
     &                                        PY = VISUAL_ORIGIN(2) + D0
                IF(ORIENTATION(IORIENTATION) .EQ. 'yz')
     &                                        PX = VISUAL_ORIGIN(1) + D0

                CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                     GAO,
     &                                     GAO1,
     &                                     PX,PY,PZ,
     &                                     NCNT,
     &                                     BUF,
     &                                     QUANTITY,
     &                                     WORK,KFREE,LFREE)

!               write to file

                WRITE(LUNIT_2D,'(7F20.10)') PX, PY, PZ,
     &                                      QUANTITY(0,1),
     &                                      QUANTITY(1,1),
     &                                      QUANTITY(2,1),
     &                                      QUANTITY(3,1)
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
          STEP_3D(I,0) = CORD(I,1)
          STEP_3D(I,I) = CORD(I,1)
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

        WRITE(LUNIT_3D,'(A)') 'some'
        WRITE(LUNIT_3D,'(A)') 'text'
        DO I = 0,3
          WRITE(LUNIT_3D,'(I5,4F12.6)') NCUBE(I), (STEP_3D(J,I),J=1,3)
        ENDDO
        DO I = 1,NUCDEP
          WRITE(LUNIT_3D,'(I5,4F12.6)')
     &         NINT(ATOM_XYZ(4,I)), ATOM_XYZ(4,I), (ATOM_XYZ(J,I),J=1,3)
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
     &                                   GAO,
     &                                   GAO1,
     &                                   PX,PY,PZ,
     &                                   NCNT,
     &                                   BUF,
     &                                   QUANTITY(0,KCUBE),
     &                                   WORK,KFREE,LFREE)
            ENDDO

!           write to file

            WRITE(LUNIT_3D,'(6E13.5)') (QUANTITY(0,M),M=1,NCUBE(3))
          ENDDO
        ENDDO
      ENDIF !3D

!     close files

      IF(LBIT(LBIT_VISUAL,2)) CLOSE(LUNIT_2D, STATUS = 'KEEP')
      IF(LBIT(LBIT_VISUAL,3)) CLOSE(LUNIT_3D, STATUS = 'KEEP')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_QUANTITY_IN_POINT */
      SUBROUTINE GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                 GAO,
     &                                 GAO1,
     &                                 PX,PY,PZ,
     &                                 NCNT,
     &                                 BUF,
     &                                 QUANTITY_ADDMEUP,
     &                                 WORK,KFREE,LFREE)
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

      PARAMETER (D0 = 0.0D0)
      PARAMETER (D1 = 1.0D0)
      PARAMETER (D2 = 2.0D0)
      PARAMETER (D3 = 3.0D0)
      PARAMETER (D4 = 4.0D0)
      PARAMETER (D5 = 5.0D0)

      PARAMETER (DP12 = D1/D2)
      PARAMETER (DP14 = D1/D4)
      PARAMETER (DP23 = D2/D3)
      PARAMETER (DP32 = D3/D2)
      PARAMETER (DP53 = D5/D3)

#include "dcbgen.h"
#include "visual.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "dgroup.h"
#include "symmet.h"
#include "sdpre.h"
#include "ibtfun.h"

      LOGICAL   LBIT
      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          NCNT(*),
     &          BUF(*),
     &          QUANTITY(0:3),
     &          QUANTITY_ADDMEUP(0:3),
     &          WORK(*),
     &          VECTOR_A(3),
     &          VECTOR_B(3)

      IQ_N           =  1 !    number density
      IQ_M           =  2 !    magnetization
      IQ_DIVM        =  3 !div(magnetization)
      IQ_ROTM        =  4 !rot(magnetization)
!                       5
      IQ_J           =  6 !    current density
      IQ_DIVJ        =  7 !div(current density)
      IQ_ROTJ        =  8 !rot(current density)
      IQ_JPARA       =  9 !    paramagnetic current density
      IQ_DIVJPARA    = 10 !div(paramagnetic current density)
      IQ_ROTJPARA    = 11 !rot(paramagnetic current density)
      IQ_JDIA        = 12 !     diamagnetic current density
      IQ_VPARA       = 13 !    paramagnetic velocity density
      IQ_DIVVPARA    = 14 !div(paramagnetic velocity density)
      IQ_ROTVPARA    = 15 !rot(paramagnetic velocity density) a.k.a. vorticity
      IQ_ELF         = 16 !    ELF
      IQ_K           = 17 !    kinetic energy density
      IQ_H           = 18 !    helicity density
      IQ_GAMMA5      = 19 !    \gamma_5 density
      IQ_NUC1SPINX   = 20 !    give me a name ...
      IQ_NUC1SPINY   = 21 !    give me a name ...
      IQ_NUC1SPINZ   = 22 !    give me a name ...
      IQ_NUC2SPINX   = 23 !    give me a name ...
      IQ_NUC2SPINY   = 24 !    give me a name ...
      IQ_NUC2SPINZ   = 25 !    give me a name ...
      IQ_EDIPX       = 26 !    electric x-dipole density
      IQ_EDIPY       = 27 !    electric y-dipole density
      IQ_EDIPZ       = 28 !    electric z-dipole density
      IQ_MDIPX       = 29 !    magnetic x-dipole density
      IQ_MDIPY       = 30 !    magnetic y-dipole density
      IQ_MDIPZ       = 31 !    magnetic z-dipole density

      CALL DZERO(QUANTITY_ADDMEUP,4)

      CALL GETSOS(GAO,GAO1,UNUSED,UNUSED,NCNT,PX,PY,PZ,
     &            BUF,WORK(KFREE),LFREE,NTBAS(0),1,.FALSE.,0)

      DO IXVECTOR = 1,NR_XVECTORS
        IOFF = N2BBASXQ*(IXVECTOR-1)

        IF(XVECTOR_FILE(IXVECTOR) .EQ. 'ZORDER')
     &                                       IOFF = N2BBASXQ*NR_XVECTORS

        CALL DZERO(QUANTITY,4)

        IQ = IQ_N !number density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GETRHO(DENSITY,
     &                IREP_XVECTOR(IXVECTOR),
     &                GAO,
     &                DMAT_AO(1+IOFF),
     &                BUF)

          CALL DAXPY(1,D1,DENSITY,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_ELF !ELF
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          ELF_KIN = D0
          DO IXYZ = 1,3
            DO IC = 1,MC
              SCALAR = D0

              CALL GETRH1(SCALAR,
     &                    IC,
     &                    IREP_XVECTOR(IXVECTOR),
     &                    GAO1(1+(IXYZ-1)*NTBAS(0)),
     &                    DMAT_AO(1+IOFF),
     &                    BUF)

              ELF_KIN = ELF_KIN + SCALAR/D2
            ENDDO
          ENDDO

          FACTOR = 0.3D0*((D3*PI*PI)**DP23)

          CALL GETRHO(DENSITY,
     &                IREP_XVECTOR(IXVECTOR),
     &                GAO,
     &                DMAT_AO(1+IOFF),
     &                BUF)

          CALL GETGRHO(VECTOR_A,
     &                 IREP_XVECTOR(IXVECTOR),
     &                 GAO,
     &                 GAO1,
     &                 DMAT_AO(1+IOFF),
     &                 BUF)

          GRADIENT_SQUARED = VECTOR_A(1)**D2
     &                     + VECTOR_A(2)**D2
     &                     + VECTOR_A(3)**D2

          H = FACTOR*(DENSITY**DP53)
          G = ELF_KIN - GRADIENT_SQUARED/(8.0D0*DENSITY)

          ELF_CHI = G/H

          SCALAR = D1/(D1 + ELF_CHI*ELF_CHI)

          CALL DAXPY(1,D1,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_JDIA !diamagnetic current density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GETRHO(DENSITY,
     &                IREP_XVECTOR(IXVECTOR),
     &                GAO,
     &                DMAT_AO(1+IOFF),
     &                BUF)

          VECTOR_A(1) = -DP12*PY*DENSITY
          VECTOR_A(2) =  DP12*PX*DENSITY
          VECTOR_A(3) =  D0

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_M !magnetization
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_M(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_J !current density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_JPARA !paramagnetic current density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_JPARA(.TRUE.,
     &                   VECTOR_A,
     &                   IREP_XVECTOR(IXVECTOR),
     &                   GAO,
     &                   GAO1,
     &                   DMAT_AO(1+IOFF),
     &                   BUF)

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_DIVM !div(magnetization)
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVM(.TRUE.,
     &                  SCALAR,
     &                  IREP_XVECTOR(IXVECTOR),
     &                  GAO,
     &                  GAO1,
     &                  DMAT_AO(1+IOFF),
     &                  BUF)

          CALL DAXPY(1,D1,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_DIVJ !div(current density)
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVJ(.TRUE.,
     &                  SCALAR,
     &                  IREP_XVECTOR(IXVECTOR),
     &                  GAO,
     &                  GAO1,
     &                  DMAT_AO(1+IOFF),
     &                  BUF)

          CALL DAXPY(1,D1,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_H !helicity density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVM(.TRUE.,
     &                  SCALAR,
     &                  IREP_XVECTOR(IXVECTOR),
     &                  GAO,
     &                  GAO1,
     &                  DMAT_AO(1+IOFF),
     &                  BUF)

          CALL DAXPY(1,DP14,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_K !kinetic energy density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVJ(.TRUE.,
     &                  SCALAR,
     &                  IREP_XVECTOR(IXVECTOR),
     &                  GAO,
     &                  GAO1,
     &                  DMAT_AO(1+IOFF),
     &                  BUF)

          CALL DAXPY(1,DP12,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_ROTM !rot(magnetization)
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_ROTM(.TRUE.,
     &                  VECTOR_A,
     &                  IREP_XVECTOR(IXVECTOR),
     &                  GAO,
     &                  GAO1,
     &                  DMAT_AO(1+IOFF),
     &                  BUF)

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_ROTJ !rot(current density)
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_ROTJ(.TRUE.,
     &                  VECTOR_A,
     &                  IREP_XVECTOR(IXVECTOR),
     &                  GAO,
     &                  GAO1,
     &                  DMAT_AO(1+IOFF),
     &                  BUF)

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_ROTJPARA !rot(paramagnetic current density)
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_ROTJPARA(.TRUE.,
     &                      VECTOR_A,
     &                      IREP_XVECTOR(IXVECTOR),
     &                      GAO,
     &                      GAO1,
     &                      DMAT_AO(1+IOFF),
     &                      BUF)

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_GAMMA5 !\gamma_5 density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_GAMMA5(.TRUE.,
     &                    SCALAR,
     &                    IREP_XVECTOR(IXVECTOR),
     &                    GAO,
     &                    DMAT_AO(1+IOFF),
     &                    BUF)

          CALL DAXPY(1,D1,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_NUC1SPINX !give me a name ...
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - CORD(1,IOPTION_XVECTOR(IXVECTOR))
          RY = PY - CORD(2,IOPTION_XVECTOR(IXVECTOR))
          RZ = PZ - CORD(3,IOPTION_XVECTOR(IXVECTOR))

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          RSQ = RX*RX + RY*RY + RZ*RZ
          FACTOR = -D1/(RSQ**DP32)

          CALL DAXPY(1,FACTOR,VECTOR_B(1),1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_NUC1SPINY !give me a name ...
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - CORD(1,IOPTION_XVECTOR(IXVECTOR))
          RY = PY - CORD(2,IOPTION_XVECTOR(IXVECTOR))
          RZ = PZ - CORD(3,IOPTION_XVECTOR(IXVECTOR))

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          RSQ = RX*RX + RY*RY + RZ*RZ
          FACTOR = -D1/(RSQ**DP32)

          CALL DAXPY(1,FACTOR,VECTOR_B(2),1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_NUC1SPINZ !give me a name ...
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - CORD(1,IOPTION_XVECTOR(IXVECTOR))
          RY = PY - CORD(2,IOPTION_XVECTOR(IXVECTOR))
          RZ = PZ - CORD(3,IOPTION_XVECTOR(IXVECTOR))

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          RSQ = RX*RX + RY*RY + RZ*RZ
          FACTOR = -D1/(RSQ**DP32)

          CALL DAXPY(1,FACTOR,VECTOR_B(3),1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_NUC2SPINX !give me a name ...
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - CORD(1,IOPTION_XVECTOR(IXVECTOR))
          RY = PY - CORD(2,IOPTION_XVECTOR(IXVECTOR))
          RZ = PZ - CORD(3,IOPTION_XVECTOR(IXVECTOR))

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          RSQ = RX*RX + RY*RY + RZ*RZ
          FACTOR = -D1/(RSQ**DP32)

          CALL DAXPY(1,FACTOR,VECTOR_B(1),1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_NUC2SPINY !give me a name ...
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - CORD(1,IOPTION_XVECTOR(IXVECTOR))
          RY = PY - CORD(2,IOPTION_XVECTOR(IXVECTOR))
          RZ = PZ - CORD(3,IOPTION_XVECTOR(IXVECTOR))

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          RSQ = RX*RX + RY*RY + RZ*RZ
          FACTOR = -D1/(RSQ**DP32)

          CALL DAXPY(1,FACTOR,VECTOR_B(2),1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_NUC2SPINZ !give me a name ...
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - CORD(1,IOPTION_XVECTOR(IXVECTOR))
          RY = PY - CORD(2,IOPTION_XVECTOR(IXVECTOR))
          RZ = PZ - CORD(3,IOPTION_XVECTOR(IXVECTOR))

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          RSQ = RX*RX + RY*RY + RZ*RZ
          FACTOR = -D1/(RSQ**DP32)

          CALL DAXPY(1,FACTOR,VECTOR_B(3),1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_EDIPX !electric x-dipole density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GETRHO(DENSITY,
     &                IREP_XVECTOR(IXVECTOR),
     &                GAO,
     &                DMAT_AO(1+IOFF),
     &                BUF)

          CALL DAXPY(1,D1,DENSITY*PX,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_EDIPY !electric y-dipole density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GETRHO(DENSITY,
     &                IREP_XVECTOR(IXVECTOR),
     &                GAO,
     &                DMAT_AO(1+IOFF),
     &                BUF)

          CALL DAXPY(1,D1,DENSITY*PY,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_EDIPZ !electric z-dipole density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GETRHO(DENSITY,
     &                IREP_XVECTOR(IXVECTOR),
     &                GAO,
     &                DMAT_AO(1+IOFF),
     &                BUF)

          CALL DAXPY(1,D1,DENSITY*PZ,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_MDIPX !magnetic x-dipole density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - D0
          RY = PY - D0
          RZ = PZ - D0

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          FACTOR = DP12

          CALL DAXPY(1,FACTOR,VECTOR_B(1),1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_MDIPY !magnetic y-dipole density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - D0
          RY = PY - D0
          RZ = PZ - D0

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          FACTOR = DP12

          CALL DAXPY(1,FACTOR,VECTOR_B(2),1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_MDIPZ !magnetic z-dipole density
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - D0
          RY = PY - D0
          RZ = PZ - D0

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          FACTOR = DP12

          CALL DAXPY(1,FACTOR,VECTOR_B(3),1,QUANTITY(0),1)
        ENDIF

        CALL DAXPY(4,COMMON_FACTOR,QUANTITY(0),1,QUANTITY_ADDMEUP(0),1)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SETUP_ORBITAL_STRING */
      SUBROUTINE SETUP_ORBITAL_STRING(STRING)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dgroup.h"
#include "visual.h"

      LOGICAL   LBIT
      DIMENSION STRING(*)

      DO IFSYM = 1,NFSYM
        CALL DFILL(STRING(1+NESH(1)*(IFSYM-1)),D0,NESH(IFSYM))
        IF(LBIT(LBIT_VISUAL,7)) THEN !use orbital string

          IF(ORBITAL_STRING(IFSYM)(4:5) .EQ. '..') THEN
            READ(ORBITAL_STRING(IFSYM)(1:3),'(I3)') ISTART
            READ(ORBITAL_STRING(IFSYM)(6:8),'(I3)') IEND
            GOTO 1
          ENDIF
          IF(ORBITAL_STRING(IFSYM)(3:4) .EQ. '..') THEN
            READ(ORBITAL_STRING(IFSYM)(1:2),'(I3)') ISTART
            READ(ORBITAL_STRING(IFSYM)(5:8),'(I3)') IEND
            GOTO 1
          ENDIF
          IF(ORBITAL_STRING(IFSYM)(2:3) .EQ. '..') THEN
            READ(ORBITAL_STRING(IFSYM)(1:1),'(I3)') ISTART
            READ(ORBITAL_STRING(IFSYM)(4:8),'(I3)') IEND
            GOTO 1
          ENDIF
          WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
          WRITE(LUPRI,'(A)') 'wrong orbital string format'
          CALL QUIT('ERROR in VISUAL')
 1        CONTINUE

          IF(     (.NOT. (ISTART .LE. IEND))
     &       .OR. (.NOT. (ISTART .GE. 0))
     &       .OR. (.NOT. (IEND .LE. NESH(IFSYM)))) THEN
            WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
            WRITE(LUPRI,'(A)') 'wrong orbital string dimension'
            CALL QUIT('ERROR in VISUAL')
          ENDIF

          DO I = ISTART,IEND
            IF(ISTART .NE. 0) STRING(I+NESH(1)*(IFSYM-1)) = D1
          ENDDO

        ELSE
          CALL DFILL(STRING(1+NESH(1)*(IFSYM-1)),D1,NISH(IFSYM))
        ENDIF
      ENDDO

      RETURN
      END
