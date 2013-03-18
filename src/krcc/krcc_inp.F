C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRCC_INP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C Input reader for Kramers restricted coupled cluster with intermediates
C
C Written by Lasse Soerensen 2010
C
C***********************************************************************
#include "implicit.inc"
#include "ipoist8.inc"
#include "priunit.h"
#include "mxpdim.inc"
#include "cgas.inc"
#include "crun.inc"
#include "lucinp.inc"
* other blocks
#include "dgroup.h"
      PARAMETER (NTABLE = 25)
      PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0,D2=2.0D00)
C
      LOGICAL SET, NEWDEF, RESET, LBIT
      CHARACTER*4 REPNA(64)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, TMPLAB*8
      CHARACTER CTEMP*72, LINE*80, TEXT*20, REPA(8)*4
      DIMENSION WORK(LWORK), MULTB_TMP(64,64)
      INTEGER IKWSET(NTABLE)
C
C
      SAVE SET
      DATA TABLE /'.MAX CC','.MK2REF','.MK2DEL','.GASSH ','.GASSPC',
     &            '.MXCCVE','.PRINT ','.NACTEL','.NOCOVA','.E_CONV',
     &            '.GASSHK','.INACTI','.IBATCH','.NO4IDX','.IREADB',
     &            '.CC_EXC','.CCROOT','.MXCIV ','.MXITLR','.CONVEX',
     &            '.ONLYTR','.******','.******','.******','.******'/
      DATA SET/.FALSE./
C
#include "memint.h"
C
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
!     print*,'LWORK,LFREE',LWORK,LFREE
C
      CALL IZERO(IKWSET,NTABLE)
C
C     Initialize /CGAS/
C     ===================
C

      NGAS      = 0
      INOCOVA   = 0
C
C     Initialize /CRUN/
C     ===================
C
      MK2REF    = 0
      MK2DEL    = 0
      NO4IDX    = 0 ! Skip 4-index transformation
      I_DO_CC_EXC_E = 0
      MXCIV     = 0
      IONLYTRANS = 0
C
C     Initialize /DCOKRCI/
C     ===================
C
      DO I = 1,2
        DO J = 1, MXPNGAS
          IGSOCCX(J,I,1) = 0
        END DO
      END DO
C
      DO I = 1, MXPNGAS
        DO J = 1, 2
          DO K =1, NFSYM
            NGSSH(K,I) = 0
            IF(J.EQ.1) THEN
              NGSSHK(K,I,1) = 0
            ELSE
              NGSSHK(K,I,2) = 0
            END IF
          END DO
        END DO
      END DO
C
      DO I =1,2
        KRCC_NISH(I) = 0
      END DO
C
C
C     Process input for *KRCCCALC
C     ===========================
C
      NEWDEF = (WORD .EQ. '*KRCCCA')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25)
     &                     , I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in KRCCCA.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in KRCICA.')
    1          CONTINUE
C&&& .MAX CC ITER...: maximum number of initial CC iterations.
                  READ(LUCMD,*) MAXIT
                  IKWSET(1) = 1
               GO TO 100
    2          CONTINUE
C&&& MK2REFERENCE: 2 * M_K reference
                  READ(LUCMD,*) MK2REF
                  IKWSET(2) = 1
               GO TO 100
    3          CONTINUE
C&&& MK2DELTA: 2 * DELTA M_K
                  READ(LUCMD,*) MK2DEL
                  IKWSET(3) = 1
               GO TO 100
    4          CONTINUE
C&&& GASSH: GAS setup
                  READ(LUCMD,*) NGAS
                  IF (NGAS .LT. 1 .OR. NGAS .GT. MXPNGAS) THEN
                     CALL QUIT('*** ERROR in KRCC *** ' //
     &                         'illegal value of NGAS')
                  END IF
                  DO I = 1, NGAS
                     READ(LUCMD,*) (NGSSH(J,I),J=1,NFSYM)
                  END DO
                  IKWSET(4) = 1
               GO TO 100
    5          CONTINUE
C&&& GASSPC: GAS space constraints
                  IF (NGAS .EQ. 0) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in KRCC ***'
     &                  //' .GASSH must be specified before .GASSPC!'
                     CALL QUIT('*** ERROR in KRCC ***')
                  END IF
C Make loop over ISPC for multiple calculations
                  ISPC = 1
                  DO I = 1, NGAS
                     READ(LUCMD,*) (IGSOCCX(I,J,ISPC),J=1,2)
                  END DO
                  IKWSET(5) = 1
               GO TO 100
    6          CONTINUE
C&&& .MXCCVE: Number of CC vectors allowed for subspace
                  READ(LUCMD,*) MXCCV
                  IKWSET(6) = 1
               GO TO 100
    7          CONTINUE
C&&& .PRINT : General print level in KRCI-module
                  READ(LUCMD,*) IPRNT
                  IKWSET(7) = 1
               GO TO 100
    8          CONTINUE
C&&& .NACTEL: number of electrons
                  READ(LUCMD,*) NACTEL
                  IKWSET(8) = 1
               GO TO 100
    9          CONTINUE
C&&& .NOCOVA: Max excitation level with any external indices
                  READ(LUCMD,*) INOCOVA
                  IKWSET(9) = 1
               GO TO 100
   10          CONTINUE
C&&& .E_CONV: Convergence of cc vector function
                  READ(LUCMD,*) CCCONV 
                  IKWSET(10) = 1
               GO TO 100
   11          CONTINUE
C&&& GASSHK: GAS setup for spinors
                  READ(LUCMD,*) NGAS
                  IF (NGAS .LT. 1 .OR. NGAS .GT. MXPNGAS) THEN
                     CALL QUIT('*** ERROR in KRCC *** ' //
     &                         'illegal value of NGAS')
                  END IF
                  DO I = 1, NGAS
                    DO J = 1, 2
                      DO K =1, NFSYM
                        IF(J.EQ.1) THEN
                          READ(LUCMD,*) NGSSHK(K,I,1)
                        ELSE
                          READ(LUCMD,*) NGSSHK(K,I,2)
                        END IF
                      END DO
                    END DO
                  END DO
                  IKWSET(11) = 1
               GO TO 100
   12          CONTINUE
               READ(LUCMD,*) (KRCC_NISH(I),I=1,NFSYM)
               GO TO 100
   13          CONTINUE
C&&& IBATCH: Change default batch length according to memory and speed
                 READ(LUCMD,*) LCCB
                 IKWSET(13) = 1
               GO TO 100
   14          CONTINUE
C&&& Skip the 4-index transformation - used for restarrted calculations
                 IKWSET(14) = 1
                 NO4IDX = 1
               GO TO 100
   15          CONTINUE
C&&& The length of the blocks written and read from disc
                 READ(LUCMD,*) I_READ_BLOCK
                 IKWSET(15) = 1
               GO TO 100
   16          CONTINUE
C&&& Do CC excitation energies
                 I_DO_CC_EXC_E = 1
                 IKWSET(16) = 1
               GO TO 100
   17          CONTINUE
C&&& Number of roots per symmetry
                 READ(LUCMD,*) (NEXC_PER_SYM(I),I=1,NBSYM)
                 IKWSET(17) = 1
               GO TO 100
   18          CONTINUE
C&&& Number of vectors per root
                 READ(LUCMD,*) MXCIV
                 IKWSET(18) = 1
               GO TO 100
   19          CONTINUE
C&&& Max number of iterations for linear response
                 READ(LUCMD,*) MXITLR
                 IKWSET(19) = 1
               GO TO 100
   20          CONTINUE
C&&& Convergence treshold for excited states
                 READ(LUCMD,*) CCCONV_EX
                 IKWSET(20) = 1
               GO TO 100
   21          CONTINUE
C&&& Only transformation no CC. Will just exit after transformation
                 IONLYTRANS = 1
                 IKWSET(21) = 1
               GO TO 100
   22          CONTINUE
C&&& 
                 IKWSET(22) = 1
               GO TO 100
   23          CONTINUE
C&&& 
                 IKWSET(23) = 1
               GO TO 100
   24          CONTINUE
C&&& 
                 IKWSET(24) = 1
               GO TO 100
   25          CONTINUE
C&&& 
                 IKWSET(25) = 1
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in KRCCCALC.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in KRCCCALC.')
            END IF
      END IF
  300 CONTINUE
C
      IF (.NOT. NEWDEF) GOTO 400
C
C SET DEFAULTS AND MAKE SMALL CHECKS
C
      IF(IKWSET(2).EQ.0) THEN
        MAXDEFIT = 50
        MAXIT = MAXDEFIT
        WRITE(LUPRI,'(A)') 
     &     ' Setting max number of iterations to ', MAXDEFIT
      END IF
             
C
C Calculate number of orbitals per gasspace and total number of orbitals
C
      print*,'NFSYM',NFSYM
      IF(IKWSET(4).EQ.1) THEN
        IF (NGAS .GT. 0) THEN
          NTOTORB = 0
          DO J = 1, NGAS
            NGSOBT(J) = 0
            DO I = 1, NFSYM
              NGSOBT(J) = NGSOBT(J) + NGSSH(I,J)
            END DO
            NTOTORB = NTOTORB + NGSOBT(J)
            IF(NGSOBT(J).GT.MXTSOB) THEN
              WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &        ' *** ERROR in KRCC *** Too many orbitals: ',
     &        NGSOBT(J),' in GAS ',J, ' max is ',MXPORB
              CALL QUIT('*** ERROR in KRCC *** Too many orbitals!')
            END IF
          END DO
        END IF
        IF(NTOTORB.GT.MXPORB) THEN
          WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &       ' *** ERROR in KRCC *** Too many orbitals: ',
     &       NTOTORB, ' -- max is ',MXPORB
          CALL QUIT('*** ERROR in KRCC *** Too many orbitals!')
        END IF
      END IF
C
C Change of convergence threshold
C
      IF(IKWSET(10).EQ.1) THEN
        WRITE(6,*) ' Convergence threshold changed to: ',CCCONV
      ELSE
        CCCONV = 1.0D-10
        WRITE(6,*) ' Default convergence threshold of ',CCCONV,' used'
      END IF
C
C Calculate number of orbitals per gasspace and total number of orbitals
C Using new input type
C
      IF(IKWSET(11).EQ.1) THEN
        IF (NGAS .GT. 0) THEN
          NTOTORB = 0
          DO J = 1, NGAS
            NGSOBTK(J,1) = 0
            NGSOBTK(J,2) = 0
            DO I = 1, NFSYM
              DO K = 1, 2
                NGSOBTK(J,K) = NGSOBTK(J,K) + NGSSHK(I,J,K)
              END DO
            END DO
            NTOTORB = NTOTORB + NGSOBTK(J,1) + NGSOBTK(J,2)
            IF(NGSOBTK(J,1) + NGSOBTK(J,2).GT.MXTSOB) THEN
              WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &        ' *** ERROR in KRCC *** Too many orbitals: ',
     &        NGSOBTK(J,1)+NGSOBTK(J,2),' in GAS ',J, ' max is ',MXPORB
              CALL QUIT('*** ERROR in KRCC *** Too many orbitals!')
            END IF
          END DO
        END IF
        IF(NTOTORB.GT.MXPORB) THEN
          WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &       ' *** ERROR in KRCC *** Too many orbitals: ',
     &       NTOTORB, ' -- max is ',MXPORB
          CALL QUIT('*** ERROR in KRCC *** Too many orbitals!')
        END IF
      END IF
C
C Default batch length
C
      IF(IKWSET(13).EQ.0) THEN
        WRITE(6,*) ' Batch lengt set to default 1000'
        WRITE(6,*) ' Decrease to reduce memory or increase for speed'
        LCCB = 1000
      ELSE
        WRITE(6,*) ' You have changed the default batch length to',LCCB
      END IF
C
C No 4-index transformation
C
      IF(IKWSET(14).EQ.0) THEN
        WRITE(6,*) ' Will do the 4-index transformation'
      ELSE
        WRITE(6,*) ' This is a restarted calculation from integrals'
        WRITE(6,*) ' No 4-index transformation will be performed'
      END IF
C
C The block size in which files are written and read => smaller blocks =
C less memory but more reads
C
      IF(IKWSET(15).EQ.0) THEN
        I_READ_BLOCK = 100000
        WRITE(6,*) ' The default block sizes are : ',I_READ_BLOCK
      ELSE
        WRITE(6,*) ' You have changed the default block size to : ',
     &               I_READ_BLOCK
      END IF
C
C Consistency check for CC excitation energies
C
      IF(I_DO_CC_EXC_E.EQ.1) THEN
        IF(IKWSET(16).EQ.0) THEN
          WRITE(6,*) ' You need to specify the number of roots '
          CALL QUIT(' You need to specify the number of roots ')
        END IF
        IF(IKWSET(18).EQ.0) THEN
C Default number of vectors per root is 2
          MXCIV = 2
        END IF
      ELSE 
        IF(IKWSET(16).EQ.1) THEN
          CALL QUIT(' I want you to tell me to do excitation energies ')
        END IF
      END IF
      IF(IKWSET(19).EQ.0) THEN
        MXITLR = MAXIT
      END IF
      IF(IKWSET(20).EQ.0) THEN
        CCCONV_EX = CCCONV
      END IF
C
C Setting some defaults that i am too lazy to change in the code
C
        NCMBSPC = ISPC
        NCISPC = ISPC
        DO ICISPC = 1, NCISPC
          LCMBSPC(ICISPC) = 1
          ICMBSPC(1,ICISPC) = ICISPC
        END DO
        LCSBLK = 6400000
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &   ' KRCC: General set-up for KR-CC calculation:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
      WRITE(LUPRI,'(A,I4)')
     &   ' Number of electrons ', NACTEL
      print*,NACTEL
      IF(IKWSET(11).EQ.1) THEN
      WRITE(LUPRI,'(A)')
     &   ' Total number of spinors ', NTOTORB
        IF (NGAS .GT. 0) THEN
          WRITE(LUPRI,'(A,I3,A)')
     &    ' * GAS space setup for ',NGAS,' GAS space(s) : '
          DO I = 1, NGAS
            WRITE(LUPRI,'(A,I3,A,2I4)')
     &           '   - GAS space ',I,'       : ',
     &           ' Unbarred spinors ', (NGSSHK(J,I,1),J=1,NFSYM),
     &           ' Barred spinors ', (NGSSHK(J,I,2),J=1,NFSYM)
            WRITE(LUPRI,'(A,I3,A1,I3,A)')
     &      '    (constraints: min/max active electrons after space : ',
     &           IGSOCCX(I,1,1),'/',IGSOCCX(I,2,1),')'
          END DO
        END IF
      END IF
      IF(IKWSET(9).EQ.1) THEN
        WRITE(6,*) ' Using NOCOVA approximation '
        WRITE(6,*) ' External indices restricted to ',INOCOVA
      END IF
      IF(IKWSET(4).EQ.1) THEN
        WRITE(LUPRI,'(A,I4)')
     &  ' Total number of Kramers pairs ', NTOTORB
        IF (NGAS .GT. 0) THEN
          WRITE(LUPRI,'(A,I3,A)')
     &    ' * GAS space setup for ',NGAS,' GAS space(s) : '
          DO I = 1, NGAS
            WRITE(LUPRI,'(A,I3,A,2I4)')
     &           '   - GAS space ',I,'       : ',
     &           (NGSSH(J,I),J=1,NFSYM)
            WRITE(LUPRI,'(A,I3,A1,I3,A)')
     &      '    (constraints: min/max active electrons after space : ',
     &           IGSOCCX(I,1,1),'/',IGSOCCX(I,2,1),')'
          END DO
        END IF
      END IF
      WRITE(LUPRI,'(A,I4)')
     &   ' Kramers projection of the reference ', MK2REF
      WRITE(LUPRI,'(A,I4)')
     &   ' Maximum flipping of Kramers projection ', MK2DEL
      WRITE(LUPRI,'(A,I4)')
     &   ' Number of inactive shells ',(KRCC_NISH(I),I=1,NFSYM)
      IF(I_DO_CC_EXC_E.EQ.1) THEN
        WRITE(6,*) '                                       '
        WRITE(6,*) ' Will calculate excited state energies '
        WRITE(6,*) '                                       '
        WRITE(6,*) ' ------------------------------------- '
        WRITE(6,*) '                                       '
        WRITE(6,*) ' Number of excited states in each symmetry '
        WRITE(6,*) '                                       '
        DO I = 1,NBSYM
          WRITE(6,*) 'Number of roots : ',NEXC_PER_SYM(I),
     &               ' in symmetry : ',I
        END DO
      END IF
C
C
C
  400 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

