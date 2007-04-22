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

      CALL MEMGET('REAL',KCMO         ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG         ,NORBT,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_MO     ,N2ORBXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_AO     ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_AO_TEMP,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSTRING1     ,NESH(1),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSTRING2     ,NESH(2),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO         ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO1        ,NTBAS(0)*3,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCNT        ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF         ,NTBAS(0)*4,
     &            WORK,KFREE,LFREE)

      CALL SETUP_ORBITAL_STRINGS(WORK(KSTRING1),WORK(KSTRING2))

!     coefficients

      CALL DZERO(WORK(KCMO),N2BBASXQ)
      CALL DZERO(WORK(KEIG),NORBT)

      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','VISUAL')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),UNUSED,
     &            UNUSED,6,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS = 'KEEP')

      CALL FNDOCC(WORK(KEIG),WORK,KFREE,LFREE)

!     density matrices

      IF(LBIT(LBIT_VISUAL,4)) THEN !using the solution vector
        CALL SETUP_FO_DMAT_AO(WORK(KCMO),
     &                        WORK(KDMAT_AO), !this comes back
     &                        WORK(KDMAT_AO_TEMP),
     &                        WORK(KSTRING1),WORK(KSTRING2),
     &                        WORK(KFREE),LFREE)

!       in some cases the ZO DMAT_AO is needed later
!       so it is set up here as well (and called DMAT_AO_TEMP)
!       could be checked and skipped ...
        CALL SETUP_ZO_DMAT_AO(WORK(KCMO),
     &                        WORK(KDMAT_MO),
     &                        WORK(KDMAT_AO_TEMP),
     &                        WORK(KSTRING1),WORK(KSTRING2),
     &                        WORK(KFREE),LFREE)
      ELSE !"normal"
        CALL SETUP_ZO_DMAT_AO(WORK(KCMO),
     &                        WORK(KDMAT_MO),
     &                        WORK(KDMAT_AO),
     &                        WORK(KSTRING1),WORK(KSTRING2),
     &                        WORK(KFREE),LFREE)
      ENDIF

      CALL VISUAL_DRIVER(WORK(KDMAT_AO),
     &                   WORK(KDMAT_AO_TEMP),
     &                   WORK(KGAO),
     &                   WORK(KGAO1),
     &                   WORK(KNCNT),
     &                   WORK(KBUF),
     &                   WORK(KFREE),LFREE)

!     release memory

      CALL MEMREL('VISUAL',WORK,KCMO,KCMO,KFREE,LFREE)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SETUP_ORBITAL_STRINGS */
      SUBROUTINE SETUP_ORBITAL_STRINGS(STRING1,STRING2)
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

      DIMENSION STRING1(*),
     &          STRING2(*)

      CALL DFILL(STRING1,D0,NESH(1))
      CALL DFILL(STRING1,D1,NISH(1))
      IF(NFSYM .GT. 1) THEN
        CALL DFILL(STRING2,D0,NESH(2))
        CALL DFILL(STRING2,D1,NISH(2))
      ENDIF

      RETURN
      END
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
      DATA TABLE /'.XXXXXX','.2D    ','.3D    ','.READXV',
!                   5         6         7         8
     &            '.SKIPEP','.SKIPEE','.STRING','.INTEGR',
!                   9        10        11        12
     &            '.3DADD ','.CUBE  ','.NXVEC ','.XXXXXX',
!                  13        14        15        16
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
!                  17        18        19        20
     &            '.N     ','.M     ','.DIVM  ','.ROTM  ',
!                  21        22        23        24
     &            '.RXJ   ','.J     ','.DIVJ  ','.ROTJ  ',
!                  25        26        27        28
     &            '.JPARA ','.DIVJPA','.ROTJPA','.JDIA  ',
!                  29        30        31        32
     &            '.VPARA ','.DIVVPA','.ROTVPA','.ELF   ',
!                  33        34        35        36
     &            '.K     ','.H     ','.GAMMA5','.XXXXXX',
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

      NCUBE(1)             = 80
      NCUBE(2)             = 80
      NCUBE(3)             = 80

      NXVEC                = 1

      IREP_VISUAL          = 0

      LBIT_VISUAL          = 0
      LBIT_VISUAL_QUANTITY = 0

      ADD_3D               = 4.0D0

      CALL LZERO(QUANTITY_IS_A_VECTOR,NR_QUANTITIES)

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
          GOTO 100

    2     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,2)
          GOTO 100

    3     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,3)
          GOTO 100

    4     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,4)
          READ(LUCMD,'(A6)') XVECFILE
          GOTO 100

    5     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,5)
          GOTO 100

    6     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,6)
          GOTO 100

    7     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL,7)
          READ(LUCMD,'(A)') ORBSTRING1
          IF(NFSYM .GT. 1) READ(LUCMD,'(A)') ORBSTRING2
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
          READ(LUCMD,*) NXVEC
          GOTO 100

   12     CONTINUE
   13     CONTINUE
   14     CONTINUE
   15     CONTINUE
   16     CONTINUE
          GOTO 100

   17     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY, 1)
          QUANTITY_NAME(1)         = 'n'
          GOTO 100

   18     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY, 2)
          QUANTITY_NAME(2)         = 'm'
          QUANTITY_IS_A_VECTOR(2)  = .TRUE.
          GOTO 100

   19     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY, 3)
          QUANTITY_NAME(3)         = 'divm'
          GOTO 100

   20     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY, 4)
          QUANTITY_NAME(4)         = 'rotm'
          QUANTITY_IS_A_VECTOR(4)  = .TRUE.
          GOTO 100

   21     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY, 5)
          QUANTITY_NAME(5)         = 'rxj'
          QUANTITY_IS_A_VECTOR(5)  = .TRUE.
          GOTO 100

   22     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY, 6)
          QUANTITY_NAME(6)         = 'j'
          QUANTITY_IS_A_VECTOR(6)  = .TRUE.
          GOTO 100

   23     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY, 7)
          QUANTITY_NAME(7)         = 'divj'
          GOTO 100

   24     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY, 8)
          QUANTITY_NAME(8)         = 'rotj'
          QUANTITY_IS_A_VECTOR(8)  = .TRUE.
          GOTO 100

   25     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY, 9)
          QUANTITY_NAME(9)         = 'jpara'
          QUANTITY_IS_A_VECTOR(9)  = .TRUE.
          GOTO 100

   26     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,10)
          QUANTITY_NAME(10)        = 'divjpara'
          GOTO 100

   27     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,11)
          QUANTITY_NAME(11)        = 'rotjpara'
          QUANTITY_IS_A_VECTOR(11) = .TRUE.
          GOTO 100

   28     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,12)
          QUANTITY_NAME(12)        = 'jdia'
          QUANTITY_IS_A_VECTOR(12) = .TRUE.
          GOTO 100

   29     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,13)
          QUANTITY_NAME(13)        = 'vpara'
          QUANTITY_IS_A_VECTOR(13) = .TRUE.
          GOTO 100

   30     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,14)
          QUANTITY_NAME(14)        = 'divvpara'
          GOTO 100

   31     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,15)
          QUANTITY_NAME(15)        = 'rotvpara'
          QUANTITY_IS_A_VECTOR(15) = .TRUE.
          GOTO 100

   32     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,16)
          QUANTITY_NAME(16)        = 'elf'
          GOTO 100

   33     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,17)
          QUANTITY_NAME(17)        = 'k'
          GOTO 100

   34     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,18)
          QUANTITY_NAME(18)        = 'h'
          GOTO 100

   35     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,19)
          QUANTITY_NAME(19)        = 'gamma5'
          GOTO 100

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
C     /* GET_QUANTITY_IN_POINT */
      SUBROUTINE GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                 GAO,
     &                                 GAO1,
     &                                 PX,PY,PZ,
     &                                 NCNT,
     &                                 BUF,
     &                                 QUANTITY,
     &                                 WORK,KFREE,LFREE)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0)
      PARAMETER (D1 = 1.0D0)
      PARAMETER (D2 = 2.0D0)
      PARAMETER (D3 = 3.0D0)
      PARAMETER (D4 = 4.0D0)
      PARAMETER (D5 = 5.0D0)

      PARAMETER (DP12 = D1/D2)
      PARAMETER (DP23 = D2/D3)
      PARAMETER (DP53 = D5/D3)

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
     &          QUANTITY(0:3,NR_QUANTITIES),
     &          WORK(*),
     &          VECTOR(3)

      IQ_N        =  1 !    number density
      IQ_M        =  2 !    magnetization
      IQ_DIVM     =  3 !div(magnetization)
      IQ_ROTM     =  4 !rot(magnetization)
      IQ_RXJ      =  5 !    r \times j
      IQ_J        =  6 !    current density
      IQ_DIVJ     =  7 !div(current density)
      IQ_ROTJ     =  8 !rot(current density)
      IQ_JPARA    =  9 !    paramagnetic current density
      IQ_DIVJPARA = 10 !div(paramagnetic current density)
      IQ_ROTJPARA = 11 !rot(paramagnetic current density)
      IQ_JDIA     = 12 !     diamagnetic current density
      IQ_VPARA    = 13 !    paramagnetic velocity density
      IQ_DIVVPARA = 14 !div(paramagnetic velocity density)
      IQ_ROTVPARA = 15 !rot(paramagnetic velocity density) a.k.a. vorticity
      IQ_ELF      = 16 !    ELF
      IQ_K        = 17 !    kinetic energy density
      IQ_H        = 18 !    helicity density
      IQ_GAMMA5   = 19 !    \gamma_5 density

      CALL DZERO(QUANTITY,NR_QUANTITIES*4)

      CALL GETSOS(GAO,GAO1,UNUSED,UNUSED,NCNT,PX,PY,PZ,
     &            BUF,WORK(KFREE),LFREE,NTBAS(0),1,.FALSE.,0)

      IQ = IQ_N !number density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
        CALL GETRHO(QUANTITY(0,IQ),IREP_VISUAL,GAO,DMAT_AO,BUF)
      ENDIF

      IQ = IQ_ELF !ELF
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
        ELF_KIN = D0
        DO IXYZ = 1,3
          DO IC = 1,MC
            SCALAR = D0
            CALL GETRH1(SCALAR,IC,IREP_VISUAL,GAO1(1+(IXYZ-1)*NTBAS(0)),
     &                  DMAT_AO,BUF)
            ELF_KIN = ELF_KIN + SCALAR/D2
          ENDDO
        ENDDO

        PI =  3.1415926536D0
        FACTOR = 0.3D0*((D3*PI*PI)**DP23)

        CALL GETRHO(DENSITY,IREP_VISUAL,GAO,DMAT_AO,BUF)

        CALL GETGRHO(VECTOR,IREP_VISUAL,GAO,GAO1,DMAT_AO,BUF)
        GRADIENT_SQUARED = VECTOR(1)**D2
     &                   + VECTOR(2)**D2
     &                   + VECTOR(3)**D2

        H = FACTOR*(DENSITY**DP53)
        G = ELF_KIN - GRADIENT_SQUARED/(8.0D0*DENSITY)

        ELF_CHI = G/H

        QUANTITY(0,IQ) = D1/(D1 + ELF_CHI*ELF_CHI)
      ENDIF

      IQ = IQ_JDIA !diamagnetic current density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
        CALL GETRHO(DENSITY,0,GAO,DMAT_AO,BUF)

        QUANTITY(1,IQ) = -(D1/D2)*PY*DENSITY
        QUANTITY(2,IQ) =  (D1/D2)*PX*DENSITY
        QUANTITY(3,IQ) =   D0
      ENDIF

      IQ = IQ_JPARA !paramagnetic current density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_JPARA(.TRUE.,
     &                 QUANTITY(1,IQ),
     &                 IREP_VISUAL,
     &                 GAO,
     &                 GAO1,
     &                 DMAT_AO,
     &                 BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

      IQ = IQ_M !magnetization
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_M(.TRUE.,
     &             QUANTITY(1,IQ),
     &             IREP_VISUAL,
     &             GAO,
     &             GAO1,
     &             DMAT_AO,
     &             BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

      IQ = IQ_J !current density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_J(.TRUE.,
     &             QUANTITY(1,IQ),
     &             IREP_VISUAL,
     &             GAO,
     &             GAO1,
     &             DMAT_AO,
     &             BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

      IQ = IQ_RXJ !r \times j
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_J(.TRUE.,
     &             VECTOR,
     &             IREP_VISUAL,
     &             GAO,
     &             GAO1,
     &             DMAT_AO,
     &             BUF)

        QUANTITY(1,IQ) = PY*VECTOR(3) - PZ*VECTOR(2)
        QUANTITY(2,IQ) = PZ*VECTOR(1) - PX*VECTOR(3)
        QUANTITY(3,IQ) = PX*VECTOR(2) - PY*VECTOR(1)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

      IQ = IQ_DIVM !div(magnetization)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_DIVM(.TRUE.,
     &                QUANTITY(0,IQ),
     &                IREP_VISUAL,
     &                GAO,
     &                GAO1,
     &                DMAT_AO,
     &                BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

      IQ = IQ_DIVJ !div(current density)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_DIVJ(.TRUE.,
     &                QUANTITY(0,IQ),
     &                IREP_VISUAL,
     &                GAO,
     &                GAO1,
     &                DMAT_AO,
     &                BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

      IQ = IQ_ROTM !rot(magnetization)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_ROTM(.TRUE.,
     &                QUANTITY(1,IQ),
     &                IREP_VISUAL,
     &                GAO,
     &                GAO1,
     &                DMAT_AO,
     &                BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

      IQ = IQ_ROTJ !rot(current density)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_ROTJ(.TRUE.,
     &                QUANTITY(1,IQ),
     &                IREP_VISUAL,
     &                GAO,
     &                GAO1,
     &                DMAT_AO,
     &                BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

      IQ = IQ_ROTJPARA !rot(paramagnetic current density)
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_ROTJPARA(.TRUE.,
     &                    QUANTITY(1,IQ),
     &                    IREP_VISUAL,
     &                    GAO,
     &                    GAO1,
     &                    DMAT_AO,
     &                    BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

      IQ = IQ_GAMMA5 !\gamma_5 density
      IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
#if defined (MOD_CDFT)
        CALL GET_GAMMA5(.TRUE.,
     &                  QUANTITY(0,IQ),
     &                  IREP_VISUAL,
     &                  GAO,
     &                  GAO1,
     &                  DMAT_AO,
     &                  BUF)
#else
        CALL QUIT('not implemented')
#endif
      ENDIF

!     calculate norm for quantities that are vectors

      DO IQ = 1,NR_QUANTITIES
        IF(LBIT(LBIT_VISUAL_QUANTITY,IQ)) THEN
          IF(QUANTITY_IS_A_VECTOR(IQ)) THEN
            QUANTITY(0,IQ) = QUANTITY(1,IQ)**D2
     &                     + QUANTITY(2,IQ)**D2
     &                     + QUANTITY(3,IQ)**D2

            QUANTITY(0,IQ) = DSQRT(QUANTITY(0,IQ))
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
