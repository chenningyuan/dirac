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
     &            '.V     ','.DIVV  ','.ROTV  ','.ELF   ',
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

      DO IQUANTITY = 1,NR_QUANTITIES
        QUANTITY_IS_A_VECTOR(IQUANTITY) = .FALSE.
      ENDDO

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
          QUANTITY_NAME(13)        = 'v'
          QUANTITY_IS_A_VECTOR(13) = .TRUE.
          GOTO 100

   30     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,14)
          QUANTITY_NAME(14)        = 'divv'
          GOTO 100

   31     CONTINUE
          CALL SET_LBIT(LBIT_VISUAL_QUANTITY,15)
          QUANTITY_NAME(15)        = 'rotv'
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
