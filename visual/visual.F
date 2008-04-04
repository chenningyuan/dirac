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
      DATA TABLE /'.LIST  ','.2D    ','.3D    ','.QUANTI',
!                   5         6         7         8
     &            '.XXXXXX','.XXXXXX','.ORBSTR','.INTEGR',
!                   9        10        11        12
     &            '.3DADD ','.CUBE  ','.COMMON','.XXXXXX',
!                  13        14        15        16
     &            '.XXXXXX','.2DORIG','.GAUGE ','.2DADD ',
!                  17        18        19        20
     &            '.2DSTEP','.XXXXXX','.XXXXXX','.XXXXXX',
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

      ORIGIN_2D(1)    = 0.0D0
      ORIGIN_2D(2)    = 0.0D0
      ORIGIN_2D(3)    = 0.0D0

      GAUGE_ORIGIN(1) = 0.0D0
      GAUGE_ORIGIN(2) = 0.0D0
      GAUGE_ORIGIN(3) = 0.0D0

      ADD_3D        = 4.0D0
      COMMON_FACTOR = 1.0D0
      ADD_2D        = 4.0D0
      STEP_2D       = 0.1D0

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
          CALL SET_LBIT(LBIT_VISUAL,1)
          READ(LUCMD,*) NR_POINTS_IN_LIST
          IF(NR_POINTS_IN_LIST .GT. MAX_NR_POINTS_IN_LIST)
     &       CALL QUIT('increase MAX_NR_POINTS_IN_LIST in visual.h')
          DO IPOINT_IN_LIST = 1,NR_POINTS_IN_LIST
            READ(LUCMD,*) XYZ_LIST(IPOINT_IN_LIST,1),
     &                    XYZ_LIST(IPOINT_IN_LIST,2),
     &                    XYZ_LIST(IPOINT_IN_LIST,3)
          ENDDO
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
          READ(LUCMD,*) ORIGIN_2D(1),
     &                  ORIGIN_2D(2),
     &                  ORIGIN_2D(3)
          GOTO 100

   15     CONTINUE
          READ(LUCMD,*) GAUGE_ORIGIN(1),
     &                  GAUGE_ORIGIN(2),
     &                  GAUGE_ORIGIN(3)
          GOTO 100

   16     CONTINUE
          READ(LUCMD,*) ADD_2D
          GOTO 100

   17     CONTINUE
          READ(LUCMD,*) STEP_2D
          GOTO 100

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
      CALL MEMGET('REAL',KGAO2     ,NTBAS(0)*6,
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
     &                   WORK(KGAO2),
     &                   WORK(KNCNT),
     &                   WORK(KBUF),
     &                   WORK(KFREE),LFREE)

      IF(LBIT(LBIT_VISUAL,8)) THEN !integration
        CALL VISUAL_INTEGRATION(WORK(KDMAT_AO),
     &                          WORK(KGAO),
     &                          WORK(KGAO1),
     &                          WORK(KGAO2),
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
     &                         GAO2,
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
     &          GAO2(*),
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

      IF(LBIT(LBIT_VISUAL,1)) THEN !list of points

        CALL AROUND('Requested density at point(s)')

        WRITE(LUPRI,'(A)')
        WRITE(LUPRI,'(A8,2A14,A18,A21,2A18)') 'r_x',
     &                                        'r_y',
     &                                        'r_z',
     &                                        'scalar',
     &                                        'x-component',
     &                                        'y-component',
     &                                        'z-component'
        WRITE(LUPRI,'(A)')

        DO IPOINT_IN_LIST = 1,NR_POINTS_IN_LIST

          CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                               GAO,
     &                               GAO1,
     &                               GAO2,
     &                               XYZ_LIST(IPOINT_IN_LIST,1),
     &                               XYZ_LIST(IPOINT_IN_LIST,2),
     &                               XYZ_LIST(IPOINT_IN_LIST,3),
     &                               NCNT,
     &                               BUF,
     &                               QUANTITY,
     &                               WORK,KFREE,LFREE)

          WRITE(LUPRI,'(3E14.6,4E18.10)') XYZ_LIST(IPOINT_IN_LIST,1),
     &                                    XYZ_LIST(IPOINT_IN_LIST,2),
     &                                    XYZ_LIST(IPOINT_IN_LIST,3),
     &                                    QUANTITY(0,1),
     &                                    QUANTITY(1,1),
     &                                    QUANTITY(2,1),
     &                                    QUANTITY(3,1)

        ENDDO

        WRITE(LUPRI,'(A)')
        CALL PRSYMB(LUPRI,'-',100,0)

        WRITE(LUPRI,'(A)')
        WRITE(LUPRI,'(A)')

      ENDIF

      IF(LBIT(LBIT_VISUAL,2)) THEN !2D

        NR_STEPS = 1 + ADD_2D*D2/STEP_2D

        DO IORIENTATION = 1,3

          NI = NR_STEPS
          NJ = NR_STEPS
          NK = NR_STEPS

          IF(ORIENTATION(IORIENTATION) .EQ. 'xy') NK = 1
          IF(ORIENTATION(IORIENTATION) .EQ. 'xz') NJ = 1
          IF(ORIENTATION(IORIENTATION) .EQ. 'yz') NI = 1

          WRITE(LUNIT_2D,'(A)')
          WRITE(LUNIT_2D,'(A)') 'plane: ' // ORIENTATION(IORIENTATION)

          DO I = 1,NI
            DO J = 1,NJ
              DO K = 1,NK

                PX = ORIGIN_2D(1) - ADD_2D + (I-1)*STEP_2D
                PY = ORIGIN_2D(2) - ADD_2D + (J-1)*STEP_2D
                PZ = ORIGIN_2D(3) - ADD_2D + (K-1)*STEP_2D

                IF(ORIENTATION(IORIENTATION) .EQ. 'xy')
     &                                        PZ = ORIGIN_2D(3)
                IF(ORIENTATION(IORIENTATION) .EQ. 'xz')
     &                                        PY = ORIGIN_2D(2)
                IF(ORIENTATION(IORIENTATION) .EQ. 'yz')
     &                                        PX = ORIGIN_2D(1)

                CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                     GAO,
     &                                     GAO1,
     &                                     GAO2,
     &                                     PX,PY,PZ,
     &                                     NCNT,
     &                                     BUF,
     &                                     QUANTITY,
     &                                     WORK,KFREE,LFREE)

!               write to file

                WRITE(LUNIT_2D,'(7E18.10)') PX, PY, PZ,
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
     &                                   GAO2,
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
     &                                 GAO2,
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

      LOGICAL   LBIT,
     &          NEED_GAO1,
     &          NEED_GAO2
      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          GAO2(*),
     &          NCNT(*),
     &          BUF(*),
     &          QUANTITY(0:3),
     &          QUANTITY_ADDMEUP(0:3),
     &          WORK(*),
     &          VECTOR_A(3),
     &          VECTOR_B(3)

      IQ_N           =  1 !                  \psi^\dag        \psi = n
      IQ_M           =  2 !                  \psi^\dag \Sigma \psi
      IQ_DIVM        =  3 !    \nabla \cdot  \psi^\dag \Sigma \psi
      IQ_ROTM        =  4 !    \nabla \times \psi^\dag \Sigma \psi
      IQ_J           =  5 !-ec               \psi^\dag \alpha \psi = j
      IQ_DIVJ        =  6 !-ec \nabla \cdot  \psi^\dag \alpha \psi
      IQ_ROTJ        =  7 !-ec \nabla \times \psi^\dag \alpha \psi
      IQ_JPARA       =  8 !                  j_para
      IQ_DIVJPARA    =  9 !    \nabla \cdot  j_para
      IQ_ROTJPARA    = 10 !    \nabla \times j_para
      IQ_JDIA        = 11 !                  j_dia (hardcoded for B along z)
      IQ_ELF         = 12 !                  ELF
      IQ_K           = 13 !                c \psi^\dag (\alpha \cdot p) \psi
      IQ_H           = 14 !                  \psi^\dag (\Sigma \cdot p) \psi
      IQ_GAMMA5      = 15 !                  \psi^\dag  \gamma^5        \psi

!     WARNING: the factor 1/c^2 is missing on purpose for compatibility with the
!              operators defined in prp/pamprp.F
      IQ_NUCDIPX     = 16 !      -1/(r_K^3) (r_K \times j)_x
      IQ_NUCDIPY     = 17 !      -1/(r_K^3) (r_K \times j)_y
      IQ_NUCDIPZ     = 18 !      -1/(r_K^3) (r_K \times j)_z

      IQ_EDIPX       = 19 !                  r_x n
      IQ_EDIPY       = 20 !                  r_y n
      IQ_EDIPZ       = 21 !                  r_z n
      IQ_MDIPX       = 22 !            -1/2 (r_gauge \times j)_x
      IQ_MDIPY       = 23 !            -1/2 (r_gauge \times j)_y
      IQ_MDIPZ       = 24 !            -1/2 (r_gauge \times j)_z

      IQ_GRADIENT    = 25 !
      IQ_LAPLACIAN   = 26 !

      NEED_GAO1 = .FALSE.
      NEED_GAO2 = .FALSE.

      DO IXVECTOR = 1,NR_XVECTORS
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_DIVM     ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_ROTM     ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_DIVJ     ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_ROTJ     ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_JPARA    ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_DIVJPARA ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_ROTJPARA ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_ELF      ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_K        ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_H        ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_GRADIENT ) NEED_GAO1 = .TRUE.
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_LAPLACIAN) NEED_GAO1 = .TRUE.

        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_LAPLACIAN) NEED_GAO2 = .TRUE.
      ENDDO

      IF(NEED_GAO1) THEN
        IF(NEED_GAO2) THEN
          CALL GETSOS(GAO,GAO1,GAO2,UNUSED,NCNT,PX,PY,PZ,
     &                BUF,WORK(KFREE),LFREE,NTBAS(0),2,.FALSE.,0)
        ELSE
          CALL GETSOS(GAO,GAO1,UNUSED,UNUSED,NCNT,PX,PY,PZ,
     &                BUF,WORK(KFREE),LFREE,NTBAS(0),1,.FALSE.,0)
        ENDIF
      ELSE
        CALL GETSOS(GAO,UNUSED,UNUSED,UNUSED,NCNT,PX,PY,PZ,
     &              BUF,WORK(KFREE),LFREE,NTBAS(0),0,.FALSE.,0)
      ENDIF

      CALL DZERO(QUANTITY_ADDMEUP,4)

      DO IXVECTOR = 1,NR_XVECTORS
        IOFF = N2BBASXQ*(IXVECTOR-1)

        IF(XVECTOR_FILE(IXVECTOR) .EQ. 'ZORDER')
     &                                       IOFF = N2BBASXQ*NR_XVECTORS

        CALL DZERO(QUANTITY,4)

        IQ = IQ_N
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          DENSITY = D0

          DO IC = 1,MC
            IF(IC .EQ. 1) THEN
              IF(.NOT. IOPTION_XVECTOR(IXVECTOR) .EQ. 2) THEN
                CALL GETRH1(DENSITY,1,IREP_XVECTOR(IXVECTOR),
     &                      GAO,DMAT_AO(1+IOFF),BUF)
              ENDIF
            ENDIF
            IF(IC .EQ. 2) THEN
              IF(.NOT. IOPTION_XVECTOR(IXVECTOR) .EQ. 1) THEN
                CALL GETRH1(DENSITY,2,IREP_XVECTOR(IXVECTOR),
     &                      GAO,DMAT_AO(1+IOFF),BUF)
              ENDIF
            ENDIF
          ENDDO

          CALL DAXPY(1,D1,DENSITY,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_GRADIENT
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GETGRHO(VECTOR_A,
     &                 IREP_XVECTOR(IXVECTOR),
     &                 GAO,
     &                 GAO1,
     &                 DMAT_AO(1+IOFF),
     &                 BUF)

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_LAPLACIAN
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GETLRHO(RXX,
     &                 RXY,
     &                 RXZ,
     &                 RYY,
     &                 RYZ,
     &                 RZZ,
     &                 GAO,
     &                 GAO1,
     &                 GAO2,
     &                 DMAT_AO(1+IOFF),
     &                 BUF)

          CALL DAXPY(1,D1,RXX + RYY + RZZ,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_ELF
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          FERMI_KIN = D0
          DO IXYZ = 1,3
            DO IC = 1,MC
              SCALAR = D0

              CALL GETRH1(SCALAR,
     &                    IC,
     &                    IREP_XVECTOR(IXVECTOR),
     &                    GAO1(1+(IXYZ-1)*NTBAS(0)),
     &                    DMAT_AO(1+IOFF),
     &                    BUF)

              FERMI_KIN = FERMI_KIN + SCALAR
            ENDDO
          ENDDO

          FERMI_KIN = FERMI_KIN/D2

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

          GRADIENT_SQUARED = DDOT(3,VECTOR_A,1,VECTOR_A,1)

          BOSON_KIN = GRADIENT_SQUARED/(8.0D0*DENSITY)

          TF_KIN = 0.3D0*((D3*PI*PI)**DP23)*(DENSITY**DP53)

          ELF_CHI = (FERMI_KIN - BOSON_KIN)/TF_KIN

          SCALAR = D1/(D1 + ELF_CHI*ELF_CHI)

          CALL DAXPY(1,D1,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IQ = IQ_JDIA
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

        IQ = IQ_M
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_M(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_J
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          CALL DAXPY(3,D1,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

        IQ = IQ_JPARA
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

        IQ = IQ_DIVM
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

        IQ = IQ_DIVJ
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

        IQ = IQ_H
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

        IQ = IQ_K
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

        IQ = IQ_ROTM
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

        IQ = IQ_ROTJ
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

        IQ = IQ_ROTJPARA
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

        IQ = IQ_GAMMA5
        IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ) THEN
          CALL GET_GAMMA5(.TRUE.,
     &                    SCALAR,
     &                    IREP_XVECTOR(IXVECTOR),
     &                    GAO,
     &                    DMAT_AO(1+IOFF),
     &                    BUF)

          CALL DAXPY(1,D1,SCALAR,1,QUANTITY(0),1)
        ENDIF

        IF((IQ_XVECTOR(IXVECTOR) .EQ. IQ_NUCDIPX) .OR.
     &     (IQ_XVECTOR(IXVECTOR) .EQ. IQ_NUCDIPY) .OR.
     &     (IQ_XVECTOR(IXVECTOR) .EQ. IQ_NUCDIPZ)) THEN

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

          IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_NUCDIPX)
     &      CALL DAXPY(1,FACTOR,VECTOR_B(1),1,QUANTITY(0),1)
          IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_NUCDIPY)
     &      CALL DAXPY(1,FACTOR,VECTOR_B(2),1,QUANTITY(0),1)
          IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_NUCDIPZ)
     &      CALL DAXPY(1,FACTOR,VECTOR_B(3),1,QUANTITY(0),1)
        ENDIF

        IF((IQ_XVECTOR(IXVECTOR) .EQ. IQ_EDIPX) .OR.
     &     (IQ_XVECTOR(IXVECTOR) .EQ. IQ_EDIPY) .OR.
     &     (IQ_XVECTOR(IXVECTOR) .EQ. IQ_EDIPZ)) THEN

          CALL GETRHO(DENSITY,
     &                IREP_XVECTOR(IXVECTOR),
     &                GAO,
     &                DMAT_AO(1+IOFF),
     &                BUF)

          IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_EDIPX)
     &      CALL DAXPY(1,D1,DENSITY*PX,1,QUANTITY(0),1)
          IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_EDIPY)
     &      CALL DAXPY(1,D1,DENSITY*PY,1,QUANTITY(0),1)
          IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_EDIPZ)
     &      CALL DAXPY(1,D1,DENSITY*PZ,1,QUANTITY(0),1)
        ENDIF

        IF((IQ_XVECTOR(IXVECTOR) .EQ. IQ_MDIPX) .OR.
     &     (IQ_XVECTOR(IXVECTOR) .EQ. IQ_MDIPY) .OR.
     &     (IQ_XVECTOR(IXVECTOR) .EQ. IQ_MDIPZ)) THEN

          CALL GET_J(.TRUE.,
     &               VECTOR_A,
     &               IREP_XVECTOR(IXVECTOR),
     &               GAO,
     &               DMAT_AO(1+IOFF),
     &               BUF)

          RX = PX - GAUGE_ORIGIN(1)
          RY = PY - GAUGE_ORIGIN(2)
          RZ = PZ - GAUGE_ORIGIN(3)

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          FACTOR = -DP12

          IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_MDIPX)
     &      CALL DAXPY(1,FACTOR,VECTOR_B(1),1,QUANTITY(0),1)
          IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_MDIPY)
     &      CALL DAXPY(1,FACTOR,VECTOR_B(2),1,QUANTITY(0),1)
          IF(IQ_XVECTOR(IXVECTOR) .EQ. IQ_MDIPZ)
     &      CALL DAXPY(1,FACTOR,VECTOR_B(3),1,QUANTITY(0),1)
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
