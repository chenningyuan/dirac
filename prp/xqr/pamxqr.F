      SUBROUTINE XQRINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for quadratic response  <<A;B,C>>
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <dcbxqr.h>
#include <mxcent.h>
#include <dcbham.h>
#include <dcbprp.h>
#include <mxorb.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <pgroup.h>
#include <dcbxpr.h>
C
      PARAMETER ( D1=1.0D0, D0=0.0D0 )
      PARAMETER ( NTABLE=28, DMIL=1.0D3 )
C
      LOGICAL SET,NEWDEF,RESET
      CHARACTER PROMPT*1,WORD*7,TABLE(NTABLE)*7,WORD1*7,PNAME*16
      DIMENSION WORK(LWORK),PNAME(3)
C
      SAVE SET
      DATA TABLE /
     &     '.PRINT ','.A OPER','.B OPER','.C OPER',
     &     '.B FREQ','.C FREQ','.THRESH','.MAXITR',
     &     '.SKIPEP','.MAXRED','.QED   ','.XLRNRM',
     &     '.XXXXXX','.XXXXXX','.INTFLG','.NOPREC',
     &     '.DIPLEN','.XXXXXX','.XXXXXX','.XXXXXX',
     &     '.ITRINT','.CNVINT','.XXXXXX','.XXXXXX',
     &     '.RESFAC','.ORBANA','.XXXXXX','.SKIPEE'/
      DATA SET/.FALSE./
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
C
C     Local initialization
C
      ILLINT = 1
      ISLINT = 1
      ISSINT = 1
      INPERR = 0
C
C     Initialize /XCBQR/
C
      NAQROP = 0
      NBQROP = 0
      NCQROP = 0
      NBQRFR = 1
      NCQRFR = 1
      BQRFR(1) = D0
      CQRFR(1) = D0
      IPRXQR = 0
      THCQR = 1.0D-5
      RESXQR = DMIL
      CNVXQR(1) = DUMMY
      CNVXQR(2) = DUMMY
      MAXQRM = 100
      ITRXQR = 30
      INTXQR = ILLINT + 2*ISLINT + 4*ISSINT
      ITRIQR(1) = 1
      ITRIQR(2) = 1
      XQRNRM = .FALSE.
      XQRDIH = .TRUE.
      QREP = .TRUE.
      QREE = .TRUE.
C
      CALL IZERO(LQROP,3*MAXLQR*2)
      IF (LEVYLE) THEN
         INTXQR = ILLINT
         QREP = .FALSE.
      END IF
C
C     Process input from CBIXQR
C
      NEWDEF = (WORD .EQ. '*QUADRA')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
         READ (LUCMD, '(A7)') WORD
 110     CONTINUE
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '.') THEN
            ICHANG = ICHANG + 1
            DO 200 I = 1, NTABLE
               IF (TABLE(I) .EQ. WORD) THEN
                  GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                 11,12,13,14,15,16,17,18,19,20,
     &                 21,22,23,24,25,26,27,28), I
               END IF
 200        CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &           '" not recognized in XQRINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in XQRINP.')
 1          CONTINUE
C     Print level
            READ(LUCMD,*) IPRXQR
            GO TO 100
 2          CONTINUE
C     A operator
            CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &           WORK,LWORK,IPRXQR)
            CALL OP1IND('NAQROP',IND1OP,LAQROP,NAQROP,INDXPR,MAXLQR)
            GO TO 100
 3          CONTINUE
C     B operator
            CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &           WORK,LWORK,IPRXQR)
            CALL OP1IND('NBQROP',IND1OP,LBQROP,NBQROP,INDXPR,MAXLQR)
            GO TO 100
 4          CONTINUE
C     C operator
            CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &           WORK,LWORK,IPRXQR)
            CALL OP1IND('NCQROP',IND1OP,LCQROP,NCQROP,INDXPR,MAXLQR)
            GO TO 100
 5          CONTINUE
C     B Frequencies
            READ (LUCMD, *) NBQRFR
            IF (NBQRFR.GT.0) THEN
               IF (NBQRFR.LE.MAXFQR) THEN
                  READ (LUCMD,*) (BQRFR(J),J=1,NBQRFR)
               ELSE
                  INPERR = INPERR + 1
                  WRITE (LUPRI,'(3(/,A,I5),/)')
     &                 '@Number of B frequencies specified :',
     &                 NBQRFR,
     &                 '@is greater than the allowed number :',
     &                 MAXFQR
                  READ (LUCMD,*) (BQRFR(J),J=1,MAXFQR),
     &                 (FFFF,J=MAXFQR+1,NBQRFR)
                  NBQRFR = MAXFQR
               END IF
            END IF
            GO TO 100
 6          CONTINUE
C     C Frequencies
            READ (LUCMD, *) NCQRFR
            IF (NCQRFR.GT.0) THEN
               IF (NCQRFR.LE.MAXFQR) THEN
                  READ (LUCMD,*) (CQRFR(J),J=1,NCQRFR)
               ELSE
                  INPERR = INPERR + 1
                  WRITE (LUPRI,'(3(/,A,I5),/)')
     &                 '@Number of C frequencies specified :',
     &                 NBQRFR,
     &                 '@is greater than the allowed number :',
     &                 MAXFQR
                  READ (LUCMD,*) (CQRFR(J),J=1,MAXFQR),
     &                 (FFFF,J=MAXFQR+1,NCQRFR)
                  NBQRFR = MAXFQR
               END IF
            END IF
            GO TO 100
 7          CONTINUE
C     Threshold for convergence:
            READ (LUCMD,*) THCQR
            GO TO 100
 8          CONTINUE
C     Maximum number of iterations in LR solver
            READ (LUCMD, '(I5)') ITRXQR
            GO TO 100
 9          CONTINUE
C     Exclude all electron-positron rotations
            QREP = .FALSE.
            GOTO 100
 10         CONTINUE
C     Max dimension of matrix in reduced system
            READ(LUCMD,*) MAXQRM
            GOTO 100
 11         CONTINUE
C
            GOTO 100
 12         CONTINUE
C     Normalize trial vectors
            XQRNRM = .TRUE.
            GOTO 100
C
 13         CONTINUE
            GOTO 100
C
 14         CONTINUE
            GOTO 100
C     Specify what two-electron integrals to include
 15         CONTINUE
            READ(LUCMD,*) ILLINT,ISLINT,ISSINT
            INTXQR = ILLINT + 2*ISLINT + 4*ISSINT
            GOTO 100
 16         CONTINUE
C     No preconditioning of initial trial vectors
            XQRDIH = .FALSE.
            GOTO 100
 17         CONTINUE
C     Specification of a complete beta tensor of dipole operators
            PNAME(1) = 'XDIPLEN'
            PNAME(2) = 'YDIPLEN'
            PNAME(3) = 'ZDIPLEN'
            DO I=1,3
               CALL XPRIND(PNAME(I),1,1,D1,PNAME(I),INDXPR,ISYXPR,
     &              ITRXPR,WORK,LWORK,IPRXQR)
               CALL OP1IND('NAQROP',IND1OP,LAQROP,NAQROP,INDXPR,MAXLQR)
               CALL OP1IND('NBQROP',IND1OP,LBQROP,NBQROP,INDXPR,MAXLQR)
               CALL OP1IND('NCQROP',IND1OP,LCQROP,NCQROP,INDXPR,MAXLQR)
            END DO
            GO TO 100
 18         CONTINUE
C
            GO TO 100
 19         CONTINUE
C
            GO TO 100
 20         CONTINUE
            GO TO 100
 21         CONTINUE
C     Number of iterations before adding SL- and SS-integrals
            READ(LUCMD,*) ITRIQR(1),ITRIQR(2)
            GO TO 100
 22         CONTINUE
C     Convergence thresholds for adding SL- and SS-integrals
            READ(LUCMD,*) CNVXQR(1),CNVXQR(2)
            GO TO 100
 23         CONTINUE
C
            GO TO 100
 24         CONTINUE
C
            GO TO 100
 25         CONTINUE
C     New trial vectors only for parameters whose residual
C     has a norm that is a fraction RESFAC of the max norm
            READ(LUCMD,*) RESXQR
            GO TO 100
 26         CONTINUE
C
            GO TO 100
 27         CONTINUE
C
            GO TO 100
 28         CONTINUE
C     Skip all e-e rotations
            QREE = .FALSE.
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 300
         ELSE
            WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &           '" not recognized in XQRINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal prompt in XQRINP.')
         END IF
      END IF
 300  CONTINUE
C
C     Process section
C
      IF (CNVXQR(1).LT.DUMMY) ITRIQR(1) = 1
      IF (CNVXQR(2).LT.DUMMY) ITRIQR(2) = 1
      MXREDZ = MAXQRM/2
      MAXQRM  = 2*MXREDZ
C
C     Print section
C      
      DOXQR = DOXQR.OR.(NAQROP.GT.0.AND.NBQROP.GT.0.AND.NCQROP.GT.0)
      IF (.NOT.DOXQR) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 
     &    'XQRINP: Electronic quadratic response properties'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I5)') 
     &   '* Print level                 :',IPRXQR
      WRITE(LUPRI,'(1X,A,I5)') 
     &   '* Maximum number of iterations:',ITRXQR
      WRITE(LUPRI,'(1X,A,I5)') 
     &   '* Max. size of reduced matrix :',MAXQRM
      WRITE(LUPRI,'(1X,A,D9.2)') 
     &   '* Threshold for convergence   :',THCQR
C
      WRITE(LUPRI,'(1X,A)') '* Include contributions from '//
     &    'the following two-electron integrals:'
      IF (ILLINT.EQ.1) WRITE(LUPRI,'(3X,A)') '- LL-integrals'
      IF (ISLINT.EQ.1) WRITE(LUPRI,'(3X,A)') '- LS-integrals'
      IF (ISSINT.EQ.1) WRITE(LUPRI,'(3X,A)') '- SS-integrals'
      IF (XQRNRM) WRITE(LUPRI,'(1X,A)')
     &     '* Trial vectors will be normalized.'
      IF (XQRDIH) THEN
        WRITE(LUPRI,'(1X,A)')
     &  '* Initial trial vectors will be preconditioned.'
      ELSE
        WRITE(LUPRI,'(1X,A)')
     &  '* No preconditioning of initial trial vectors.'
      ENDIF
      IF (.NOT.QREE) WRITE(LUPRI,'(1X,A)')
     &     '* No e-e rotation considered.'
      IF (.NOT.QREP) WRITE(LUPRI,'(1X,A)')
     &     '* No e-p rotation considered.'
C
      WRITE(LUPRI,'(1X,A)') 
     &     '* Frequencies in quadratic response in a.u.'
      WRITE(LUPRI,'(3X,A16,6(4F12.8,/,19X))') '- B-frequencies:',
     &     (BQRFR(I),I=1,NBQRFR)
      WRITE(LUPRI,'(3X,A16,6(4F12.8,/,19X))') '- C-frequencies:',
     &     (CQRFR(I),I=1,NCQRFR)
C
      WRITE(LUPRI,'(1X,A)') '  +-------------+'
      WRITE(LUPRI,'(1X,A)') '* | A-operators |'
      WRITE(LUPRI,'(1X,A)') '  +-------------+'
      DO I = 1,NAQROP
        INDXPR = LAQROP(I)
        CALL WRIXPR(I,INDXPR)
      ENDDO
      WRITE(LUPRI,'(1X,A)') '  +-------------+'
      WRITE(LUPRI,'(1X,A)') '* | B-operators |'
      WRITE(LUPRI,'(1X,A)') '  +-------------+'
      DO I = 1,NBQROP
        INDXPR = LBQROP(I)
        CALL WRIXPR(I,INDXPR)
      ENDDO
      WRITE(LUPRI,'(1X,A)') '  +-------------+'
      WRITE(LUPRI,'(1X,A)') '* | C-operators |'
      WRITE(LUPRI,'(1X,A)') '  +-------------+'
      DO I = 1,NCQROP
        INDXPR = LCQROP(I)
        CALL WRIXPR(I,INDXPR)
      ENDDO
C
  999 CONTINUE
C
      RETURN
      END
      SUBROUTINE PRPXQR(WORK,LWRK)
C*****************************************************************************
C
C     Driver routine for computing the quadratic response properties
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxqr.h>
#include <dgroup.h>
#include <dcborb.h>
C
      PARAMETER ( LUCMO = 22 )
      DIMENSION RES_LR(3,MAXQR),RES_QR(3,MAXQR)
      DIMENSION WORK(*)
C
      CALL QENTER('PRPXQR')
C
      KFREE = 1
      LFREE = LWRK
C
C     Get the MO coefficients
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL OPNFIL(LUCMO,'DFCOEF','OLD','LRCALC')
      CALL REACMO(LUCMO,'DFCOEF',WORK(KCMO),DUM,IDUM,TOTERG,4)
      CLOSE(LUCMO,STATUS='KEEP')
C
C     Set up the linear response equations to be solved.
C
      CALL QR_SETUP(WORK(KFREE),LFREE)
C
C     Solve the linear response equations to obtain the necessary 
C     linear response vectors. The vectors are stored in a direct
C     access file on disk with unit numbers LUQREE and LUQREP for
C     the electronic and positronic solutions, respectively. 
C
      CALL QRVEC(WORK(KFREE),LFREE)
C
C     Perform the vector dot product to obtain the linear
C     response function values.
C
      CALL MEMGET('REAL',KVECAEE,MZYEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECAEP,MZYEP*NZ,WORK,KFREE,LFREE)
      CALL LRCALC(RES_LR,WORK(KVECAEE),WORK(KVECAEP),WORK(KCMO),
     &     WORK(KFREE),LFREE)
C
C     Perform the matrix contraction to obtain the quadratic
C     response function values.
C
      CALL MEMGET('REAL',KVECBEE,MZYEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECBEP,MZYEP*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECCEE,MZYEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECCEP,MZYEP*NZ,WORK,KFREE,LFREE)
      CALL QRHYP(RES_QR,WORK(KVECAEE),WORK(KVECBEE),WORK(KVECCEE),
     &     WORK(KVECAEP),WORK(KVECBEP),WORK(KVECCEP),
     &     WORK(KCMO),WORK(KFREE),LFREE)
C
C     Formatted print out of result
C
C     CALL QRPRINT(LRVAL,QRVAL)
C
      CALL QEXIT('PRPXQR')
      RETURN
      END
      SUBROUTINE QR_SETUP(WORK,LWRK)
C*****************************************************************************
C
C     Determine the linear response equations to be solved.
C     For a given quadratic response function <<A;B,C>>_{wb,wc}
C     we need to solve for Na(-(wb+wc)), Nb(wb), and Nc(wc).
C     Note, however, that we make use of that there is a relation
C     between response vectors with opposite frequencies, so only
C     positive frequencies are solved for.
C
C     +-------------------------+
C     | Hermitian perturbations |
C     +-------------------------+
C     N(w) = / Z \  ==>  N(-w) = / Y \
C            \ Y*/               \ Z*/
C
C     +------------------------------+
C     | Anti-hermitian perturbations |
C     +------------------------------+
C     N(w) = / Z \  ==>  N(-w) = / Y \
C            \-Y*/               \-Z*/
C
C     NQROP  = number of unique operator to solve linear response for
C     NQRHYP = number of quadratic resonse functions
C     MZYEE  = max length of electronic response vectors
C     MZYEP  = max length of positronic response vectors
C     LQRHYP(NQRHYP,3)  = the A/B/C operator pointers
C     QRFRHYP(NQRHYP,2) = the B/C frequencies
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxpr.h>
#include <dcbxlr.h>
#include <dcbxqr.h>
#include <dgroup.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dcbxrs.h>
C
      DIMENSION WORK(*)
      LOGICAL DOHYP
C
      CALL QENTER('QR_SETUP')
C
      KFREE = 1
      LFREE = LWRK
      CALL OSTRING(WORK,KFREE,LFREE)
C
      NQROP  = 0
      NQRHYP = 0
C
      DO ICOP=1,NCQROP
      DO ICFR=1,NCQRFR
         DO IBOP=1,NBQROP
         DO IBFR=1,NBQRFR
            ISYMC = IPRPSYM(LCQROP(ICOP))
            ISYMB = IPRPSYM(LBQROP(IBOP))
            DO IAOP=1,NAQROP
               AFREQ=-(BQRFR(IBFR)+CQRFR(ICFR))
               ISYMA = IPRPSYM(LAQROP(IAOP))
C
C     Check if the QR function is to be evaluated. Considerations are made
C     to operator symmetries, permutation symmetries for the QR function,
C     and user specified optical processes. 
C
               CALL QRHYPCHK(LAQROP(IAOP),LBQROP(IBOP),LCQROP(ICOP),
     &              ISYMA,ISYMB,ISYMC,AFREQ,BQRFR(IBFR),CQRFR(ICFR),
     &              DOHYP)
C
               IF (DOHYP) THEN
                  NQRHYP = NQRHYP+1
                  IF (NQRHYP.LE.MAXQR) THEN
                     CALL INDQR(LAQROP(IAOP),ABS(AFREQ),INDEX)
                     CALL INDQR(LBQROP(IBOP),ABS(BQRFR(IBFR)),INDEX)
                     CALL INDQR(LCQROP(ICOP),ABS(CQRFR(ICFR)),INDEX)
                     LQRHYP(NQRHYP,1) = LAQROP(IAOP)
                     LQRHYP(NQRHYP,2) = LBQROP(IBOP)
                     LQRHYP(NQRHYP,3) = LCQROP(ICOP)
                     QRFRHYP(NQRHYP,1) = AFREQ
                     QRFRHYP(NQRHYP,2) = BQRFR(IBFR)
                     QRFRHYP(NQRHYP,3) = CQRFR(ICFR)
                  ELSE
                     NQRHYP = MAXQR
                     WRITE(LUPRI,'(//A,/A,I5)')
     & '@ WARNING!!! Quadratic response function not evaluated since',
     & '             the total number exceeds MAXQR=',MAXQR
                  END IF
               END IF
            END DO
         END DO
         END DO 
      END DO
      END DO
C
C     Get the maximum length of the response vectors
C
      MZYEE  = 0
      MZYEP  = 0
      SKIPEP = .NOT.QREP
      SKIPEE = .NOT.QREE
C
      DO IOP=1,NQROP
         IPT=LQROP(IOP,1)
         JSYMOP = IPRPSYM(IPT)
         JTIMOP = IPRPTIM(IPT)
         JOPSY  = JBTOF(JSYMOP-1,1)
         CALL XLRPAR(WORK,KFREE,LFREE,IPRXQR)
         MZYEE = MAX(MZYEE,2*NZXOPE)
         MZYEP = MAX(MZYEP,2*NZXOPP)
      END DO
C
C     Print out the defined quadratic response functions that remain after
C     reduction.
C
      IF (NQRHYP.GT.0) THEN
         WRITE(LUPRI,'(//A)') 
     &   ' The following quadratic response equations will be computed:'
         WRITE(LUPRI,'(/2A)') 
     &        '----------------------------------------',
     &        '-------------------------'
         WRITE(LUPRI,'(2A)')
     &        ' A operator   B operator   C operator   ',
     &        'B frequency   C frequency'
         WRITE(LUPRI,'(2A)') 
     &        '----------------------------------------',
     &        '-------------------------'
         DO I=1,NQRHYP
            WRITE(LUPRI,'(1X,3A13,F11.6,F14.6)') 
     &           PRPNAM(LQRHYP(I,1)),PRPNAM(LQRHYP(I,2)), 
     &           PRPNAM(LQRHYP(I,3)),QRFRHYP(I,2),QRFRHYP(I,3)
         END DO
         WRITE(LUPRI,'(2A,//)') 
     &        '----------------------------------------',
     &        '-------------------------'
      ELSE
         WRITE(LUPRI,'(4(/A),/)') ' *** WARNING ***',
     &        ' No nonzero Quadratic Response Functions exist for',
     &        ' the specified input. Check the input file DIRAC.INP',
     &        ' and then you are welcome back.'
         CALL QUIT('Quadratic Response Functions are zero by symmetry')
      END IF
C
      CALL QEXIT('QR_SETUP')
      RETURN
      END
      SUBROUTINE INDQR(IPOINT,FREQ,INDEX)
C*****************************************************************************
C
C     Return index to linear equation to solve. The index is also used
C     as record number in the direct access file storing the linear
C     response vectors. A unique (IPOINT,FREQ) pair gives rise to an 
C     linear response equation. If unique store in list. 
C
C     Input:
C     IPOINT points to an operator
C     FREQ is the frequency associated with the operator
C
C     Output:
C     INDEX position of linear response equation in list, return value
C           zero if new in list
C
C     Used from common:
C     LQROP list of unique operator pointers and the respective number 
C           of frequencies
C     NQROP  number of unique operator pointers
C     QRFREQ list of frequencies
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxqr.h>
C
      PARAMETER ( THRES = 1.0D-6 )
      CALL QENTER('INDQR')
C
      INDEX=0
C
      K=0
      DO IOP=1,NQROP
         IF (IPOINT.EQ.LQROP(IOP,1)) THEN
            DO IFR=1,LQROP(IOP,2)
               K=K+1
               IF ( ABS( FREQ-QRFREQ(IFR,IOP) ).LT.THRES ) THEN
                  INDEX = K
                  GOTO 100
               END IF
            END DO         
            LQROP(IOP,2)=LQROP(IOP,2)+1
            QRFREQ(LQROP(IOP,2),IOP)=FREQ
            GOTO 100 
         ELSE
            K=K+LQROP(IOP,2)
         END IF
      END DO
C
      NQROP=NQROP+1
      LQROP(NQROP,1)=IPOINT
      LQROP(NQROP,2)=1
      QRFREQ(1,NQROP)=FREQ
C
 100  CONTINUE
C
      CALL QEXIT('INDQR')
      RETURN
      END
      SUBROUTINE QRHYPCHK(IAOP,IBOP,ICOP,ISYMA,ISYMB,ISYMC,
     &     AFR,BFR,CFR,DOHYP)
C*****************************************************************************
C
C     Check if the QR function is to be evaluated. Considerations are made
C     to operator symmetries, permutation symmetries for the QR function,
C     and user specified optical processes. 
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <dcbxqr.h>
C
      PARAMETER (THD = 1.0D-8)
      LOGICAL DOHYP
      DIMENSION IOP(MAXQR,3),FRQ(MAXQR,3),ISYM(MAXQR,3)
C
      SAVE IOP,FRQ,N,ISYM
C
      DATA N/0/
C
#include <muld2h.h>
C
      CALL QENTER('QRHYPCHK')
C
      IF (MULD2H(ISYMA,MULD2H(ISYMB,ISYMC)).NE.1) THEN
         DOHYP=.FALSE.
         GOTO 100
      END IF
C
      DOHYP = .TRUE.
      DO I = 1,N
      DO J = 1,3
      DO K = 1,3
      IF (K.NE.J) THEN
         DO L = 1,3
         IF (L.NE.K .AND. L.NE.J) THEN
C
            IF (IAOP.EQ.IOP(I,J) .AND.
     &          IBOP.EQ.IOP(I,K) .AND.
     &          ICOP.EQ.IOP(I,L) .AND.
     &          ISYMA.EQ.ISYM(I,J) .AND.
     &          ISYMB.EQ.ISYM(I,K) .AND.
     &          ISYMC.EQ.ISYM(I,L) .AND.
     &          ABS(ABS(AFR)-ABS(FRQ(I,J))).LT.THD .AND.
     &          ABS(ABS(BFR)-ABS(FRQ(I,K))).LT.THD .AND.
     &          ABS(ABS(CFR)-ABS(FRQ(I,L))).LT.THD ) THEN
               DOHYP = .FALSE.
            END IF
C
         END IF
         END DO
      END IF
      END DO
      END DO
      END DO
C
      DOHYP = .TRUE.
C
C     If not done before put in list with unique calculations
C
      IF (DOHYP) THEN
         N = N+1
         IOP(N,1) = IAOP
         IOP(N,2) = IBOP
         IOP(N,3) = ICOP
         ISYM(N,1) = ISYMA
         ISYM(N,2) = ISYMB
         ISYM(N,3) = ISYMC
         FRQ(N,1) = AFR
         FRQ(N,2) = BFR
         FRQ(N,3) = CFR
      END IF
C
 100  CONTINUE
C     
      CALL QEXIT('QRHYPCHK')
C     
      RETURN
      END
      SUBROUTINE QRVEC(WORK,LWRK)
C*****************************************************************************
C
C     Solve the linear response equations associated with quadratic
C     response. 
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxpr.h>
#include <dcbxrs.h>
#include <dcbxqr.h>
#include <pgroup.h>
#include <dgroup.h>
#include <cbiqr.h>
C
      PARAMETER ( D0 = 0.0D0 )
C
      LOGICAL EXST,ZERO
      CHARACTER TTYP(-1:1)*2
      SAVE TTYP
      DATA TTYP /'T-','T0','T+'/
      DIMENSION WORK(*)
C
      CALL QENTER('QRVEC')
C
      KFREE = 1
      LFREE = LWRK
C
C     Initialize symmetry independent variables that are used in the
C     linear response solver routines.
C
      CALL SETXQR
C
      CALL MEMGET('INTE',KIBTYP,2*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBCVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBPVC,NREDM  ,WORK,KFREE,LFREE)
C
C---------------------------------------------------
C     Loop over number of defined properties.
C---------------------------------------------------
C
      DO IOP=1,NQROP
C
C     Initialize configurational and orbital parameters that
C     depend on the symmetry at hand.
C
         KSAVE  = KFREE
         INDPRP = LQROP(IOP,1)
         JSYMOP = IPRPSYM(INDPRP)
         JTIMOP = IPRPTIM(INDPRP)
         JOPSY  = JBTOF(JSYMOP-1,1)
         NFREQ  = LQROP(IOP,2)
         NEXSIM = NFREQ
         NEXSTV = NFREQ
         NEXCNV = NFREQ
         NEVECR = NREDM*NEXSIM
         NCRED  = 0
         NERED  = 0
         NPRED  = 0
         NZRED  = NCRED+NERED+NPRED
         STATIC = (NFREQ.EQ.1 .AND. QRFREQ(1,IOP).EQ.D0)
C
         CALL XLRPAR(WORK,KFREE,LFREE,IPRXQR)
C
C     Check if response vectors exist on restart file
C
         CALL READQR(EXST,PRPNAM(INDPRP),NFREQ,QRFREQ(1,IOP))
         IF (EXST) THEN
            WRITE(LUPRI,'(A,2(/A,A),/A,/)') ' *** WARNING ***',
     &           ' Response vector found in file QRINFO, and the',
     &           ' program will continue',
     &           ' assuming you have supplied the corresponding',
     &           ' direct access files',' QRVEC.EE and QRVEC.PP.'
            GO TO 30
         END IF
C
C     Print header
C
         WRITE(LUPRI,'(//A,//2A,A,A3,2X,A2/)')
     &        ' <<<  SOLVING SETS OF LINEAR EQUATIONS '//
     &        'FOR QUADRATIC RESPONSE PROPERTY >>>',
     &        ' Label: ',PRPNAM(INDPRP),
     &        'Symmetry: ',REP(JSYMOP-1),TTYP(JTIMOP)
         WRITE(LUPRI,'(1X,A,I2,A,3(/5F12.8))') 'at the following ',
     &        NFREQ,' frequencies:',(QRFREQ(I,IOP), I=1,NFREQ)
         CALL FLSHFO(LUPRI)
C
C     Get property gradient
C
         CALL MEMGET('REAL',KGPCI,NZCONFQ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KGPOE,NZXOPEQ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KGPOP,NZXOPPQ,WORK,KFREE,LFREE)
C
         print *, "indprp,jopsy,nzxope,nzxopp,kfree,lfree,iprxqr=",
     &        indprp,jopsy,nzxope,nzxopp,kfree,lfree,iprxqr
         CALL PAMGRD(INDPRP,WORK(KGPCI),WORK(KGPOE),WORK(KGPOP),
     &        JOPSY,WORK(KZCONF),WORK(KZXOPE),WORK(KZXOPP),
     &        NZCONF,NZXOPE,NZXOPP,WORK,KFREE,LFREE,IPRXQR)
C
C     Check norm of property gradient, skip if below threshold and
C     write response vector equal to zero vector to file.
C
         CALL GRDCHK(ZERO,INDPRP,NFREQ,QRFREQ(1,IOP),
     &        NZCONFQ,WORK(KGPCI),NZXOPEQ,WORK(KGPOE),
     &        NZXOPPQ,WORK(KGPOP),WORK(KFREE),LFREE)
C
         IF (ZERO) GOTO 30
C
C     Solve response equation in the reduced space spanned by the 
C     trial vectors.
C     
         CALL MEMGET('REAL',KEVALR,NEXSIM,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEVECR,NEVECR,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KLRF,NFREQ,WORK,KFREE,LFREE)
C
         CALL DCOPY(NFREQ,QRFREQ(1,IOP),1,WORK(KEVALR),1)
         CALL XRSCTL(WORK(KGPCI),WORK(KGPOE),WORK(KGPOP),
     &        WORK(KIBTYP),WORK(KIBCVC),WORK(KIBEVC),
     &        WORK(KIBPVC),WORK(KLRF),WORK(KEVALR),WORK(KEVECR),
     &        WORK,KFREE,LFREE)
C
C     Construct the response vector from the solution in the reduced
C     space and write them to direct access files LUQREE and LUQREP
C     for the electronic and positronic parts, respectively. The 
C     string 'QRVEC' indicates that the response vectors are to be used
C     to evaluate QR functions.
C
         CALL LRVEC(INDPRP,WORK(KIBTYP),WORK(KIBCVC),
     &        WORK(KIBEVC),WORK(KIBPVC),WORK(KEVALR),
     &        WORK(KEVECR),WORK(KFREE),LFREE)
C
C     Release memory used for this response property.
C
 30      CONTINUE
         CALL MEMREL('QRVEC',WORK,1,KSAVE,KFREE,LFREE)
C
C     End loop over linear response equations
C
      END DO 
C
      CALL QEXIT('QRVEC')
      RETURN
      END
      SUBROUTINE SETXQR
C*****************************************************************************
C
C     Initialize variables that are used by the linear response solver
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
C
#include <dcbxlr.h>
#include <dcbxqr.h>
#include <dcbxrs.h>
C
      CALL QENTER('SETXQR')
C
C     Transfer information from XCBXQR to other common blocks 
C     to be used in solving the linear response equations.
C
      LINEQ = .TRUE.
      LSVCFG(1) = LSFG(1)
      LSVCFG(2) = LSFG(2)
      TKNORM = XQRNRM
      DIAGHE = XQRDIH
      IPRXRS = IPRXQR
      THCXRS = THCQR
      RESFAC = RESXQR
      MAXITR = ITRXQR
      NREDM  = MAXQRM
      N2REDM = MAXQRM*MAXQRM
      LOFFTY = 0
      CNVINT(1) = CNVXQR(1)
      CNVINT(2) = CNVXQR(2)
      ITRINT(1) = ITRIQR(1)
      ITRINT(2) = ITRIQR(2)
      INTDEF = INTXQR
      STERNH = XSTERN
C
      CALL QEXIT('SETXQR')
C
      RETURN
      END
      SUBROUTINE LRVEC(INDPRP,IBTYP,IBCVC,IBEVC,IBPVC,
     &     EVALR,EVECR,WORK,LWRK)
C*****************************************************************************
C
C     Form the response vectors and write them to file given the 
C     solutions from the reduced system of trial vectors. 
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <dcbxrs.h>
#include <dcbxlr.h>
#include <dcbxqr.h>
#include <dcbibn.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),
     &          IBPVC(*),EVALR(*),EVECR(NREDM,*)
      DIMENSION WORK(*)
C
      CALL QENTER('LRVEC')
C
      KFREE = 1
      LFREE = LWRK
      NESIM = 0
      NPSIM = 0
      NCSIM = 0
      IFAC  = 1
      IF (.NOT.STATIC) IFAC = 2
C
C     Solution vector - orbital (e-e) part
C
      IF(NZXOPE.GT.0) THEN
        NESIM = NFREQ*IFAC
        CALL MEMGET('REAL',KVECZ,NZXOPE*NZ*NFREQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVECY,NZXOPE*NZ*NFREQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRSPVC,MZYEE*NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBUF,NZXOPE*NZ,WORK,KFREE,LFREE)
        CALL LRVEC1(INDPRP,JBENDX,IBTYP,IBEVC,EVALR,
     &       EVECR,WORK(KVECZ),WORK(KVECY),WORK(KRSPVC),WORK(KBUF),
     &       NZXOPE,MZYEE,WORK(KFREE),LFREE)
        CALL MEMREL('LRVEC',WORK,1,KVECZ,KFREE,LFREE)
      END IF
C
C     Solution vector - orbital (e-p) part
C
      IF(NZXOPP.GT.0) THEN
        NPSIM = NFREQ*IFAC
        CALL MEMGET('REAL',KVECZ,NZXOPP*NZ*NFREQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVECY,NZXOPP*NZ*NFREQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRSPVC,MZYEP*NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBUF,NZXOPP*NZ,WORK,KFREE,LFREE)
        CALL LRVEC1(INDPRP,JBPNDX,IBTYP,IBPVC,EVALR,
     &       EVECR,WORK(KVECZ),WORK(KVECY),WORK(KRSPVC),WORK(KBUF),
     &       NZXOPP,MZYEP,WORK(KFREE),LFREE)
        CALL MEMREL('LRVEC',WORK,1,KVECZ,KFREE,LFREE)
      END IF
C
      CALL QEXIT('LRVEC')
      RETURN
      END
      SUBROUTINE LRVEC1(INDPRP,NBTYP,IBTYP,IBVEC,
     &     EVALR,EVECR,XPOZ,XPOY,RSPVEC,BUF,NZPAR,MZYPAR,WORK,LWRK)
C*****************************************************************************
C
C     Form the response vectors and write them to file given the 
C     solutions from the reduced system of trial vectors. 
C
C     NZPAR  = half-length of current response vector
C     MZYPAR = max full-length of response vectors
C     NBTYP  = 1 for e-e and 2 e-p part, respectively (input)
C     XPOZ   = space for (Z+Y)
C     XPOY   = space for (Z-Y)
C     RSPVEC = space for a full response vector
C     EVALR  = frequencies in the reduced system (E[2] - w*S[2])*N = GP
C     EVECR  = space for eigenvectors of reduced system (read from file)
C     NFREQ  = number of frequencies (input)
C     IBTYP  = type of B trial vectors, dim=2*NREDM (input)
C     IBVEC  = pointers to electron/positron like trial vectors (input)
C     BUF    = extra space of response vector half-length
C     LUBVEC = unit number of direct access file for the trial vectors 
C     INDPRP = pointer to the property for which response is evaluated
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <dcbxrs.h>
#include <dcbxlr.h>
#include <dcbxqr.h>
#include <dcbibn.h>
#include <dcborb.h>
#include <dgroup.h>
C
#include <dcbxpr.h>
C
      PARAMETER ( D1=1.0D0, DM1=-1.0D0, DH=0.5D0 )
      CHARACTER FILNAM*6,FILRSP*8,TYP*1
      DIMENSION IBTYP(2,*),IBVEC(*),EVALR(*),EVECR(NREDM,*),
     &     XPOZ(NZPAR,NZ,NFREQ),XPOY(NZPAR,NZ,NFREQ),
     &     RSPVEC(MZYPAR,NZ),BUF(MZYPAR*NZ),WORK(*)
C
      CALL QENTER('LRVEC1')
C
      KFREE = 1
      LFREE = LWRK
C
      IF (NBTYP.EQ.JBENDX) THEN
        NBRED  = NERED
        LUBVEC = LUBOE
        FILNAM = 'PAMBOE'
        TYP    = 'E'
        LUFILE = LUQREE
        FILRSP = 'QRVEC.EE'
      ELSEIF (NBTYP.EQ.JBPNDX) THEN
        NBRED  = NPRED
        LUBVEC = LUBOP
        FILNAM = 'PAMBOP'
        TYP    = 'P'
        LUFILE = LUQREP
        FILRSP = 'QRVEC.EP'
      ELSE
        CALL QUIT('LRVEC: Unknown NBTYP!')
      ENDIF
C
C     Construct response solution vectors by summing the trial
C     vectors multiplied by the eigenvectors in the reduced
C     system. Note that the trial vectors are of the form 
C     (Z+Y) and (Z-Y).
C
      OPEN(LUBVEC,FILE=FILNAM,FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*NZPAR*NZ,STATUS='OLD')
      CALL DZERO(XPOZ,NZPAR*NZ*NFREQ)
      CALL DZERO(XPOY,NZPAR*NZ*NFREQ)
CTROND: quickfix...allocating IVECS here...
      CALL MEMGET('INTE',KIVECS,NFREQ,WORK,KFREE,LFREE)
      CALL XRSXV1(1,NBTYP,XPOZ,EVECR,NFREQ,IBTYP,IBVEC,
     &            WORK(KIVECS),BUF)
      CALL XRSXV1(-1,NBTYP,XPOY,EVECR,NFREQ,IBTYP,IBVEC,
     &            WORK(KIVECS),BUF)
      CALL MEMREL('LRVEC1',WORK,KIVECS,KIVECS,KFREE,LFREE)
      CLOSE(LUBVEC,STATUS='DELETE')
C
C     Write response vector to file: / Z \
C                                    \ Y*/
C
      OPEN(LUFILE,FILE=FILRSP,FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*MZYPAR*NZ,STATUS='UNKNOWN')
      DO IFR=1,NFREQ
C
         DO IZ=1,NZ
            CALL DCOPY(NZPAR,XPOZ(1,IZ,IFR),1,RSPVEC(1,IZ),1)
            CALL DCOPY(NZPAR,XPOZ(1,IZ,IFR),1,RSPVEC(1+NZPAR,IZ),1)
            CALL DAXPY(NZPAR,D1,XPOY(1,IZ,IFR),1,RSPVEC(1,IZ),1)
            CALL DAXPY(NZPAR,DM1,XPOY(1,IZ,IFR),1,
     &           RSPVEC(1+NZPAR,IZ),1)
            IF (IPQTOQ(IZ,JSYMOP-1).GT.1) 
     &           CALL DSCAL(NZPAR,DM1,RSPVEC(1+NZPAR,IZ),1)
         END DO
         CALL DSCAL(MZYPAR*NZ,DH,RSPVEC,1)
         CALL INDQR(INDPRP,EVALR(IFR),IREC)
         CALL WRTDAC(LUFILE,MZYPAR*NZ,RSPVEC,IREC)
C
C     Write info to QRINFO file to be used on restart
C
         DNORM=DNORM2(MZYPAR*NZ,RSPVEC,1)
         CALL WRTQR(PRPNAM(INDPRP),TYP,EVALR(IFR),IREC,DNORM)
C
      END DO
      CLOSE(LUFILE,STATUS='KEEP')
C
C     Print section
C
      IF(IPRXQR.GE.10) THEN
         WRITE(LUPRI,'(/2A,/A,5(/5F12.8))') 
     &        ' Writing response vectors with label: ',
     &        PRPNAM(INDPRP),' at frequencies: ',
     &        (EVALR(I), I=1,NFREQ)
         WRITE(LUPRI,'(2A)') ' to file: ',FILRSP
C
         CALL HEADER('X(+) vectors in LRVEC1',-1)
         CALL PRBVEC(LUPRI,XPOZ,NFREQ,NZPAR)
         CALL HEADER('X(-) vectors in LRVEC1',-1)
         CALL PRBVEC(LUPRI,XPOY,NFREQ,NZPAR)
      ENDIF
C
      CALL QEXIT('LRVEC1')
      RETURN
      END
      SUBROUTINE LRCALC(RESULT,VECBEE,VECBEP,CMO,WORK,LWRK)
C*****************************************************************************
C
C     ------------------------------------------------------------
C     Routine for computing the linear response function that are
C     solved for in the quadratic response calculation.
C     ------------------------------------------------------------
C
C     Linear response functions are defined as:
C
C     <<A;B>>_{wb} = - A[1]_j * Nb_j(w1)
C
C     Each operator A and B is assumed to span a given boson irrep
C     and be either symmetric(+) or antisymmetric(-) with respect
C     to time reversal.
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbibn.h>
#include <dcbxrs.h>
#include <dcbxpr.h>
#include <dgroup.h>
#include <dcborb.h>
#include <dcbxqr.h>
#include <dcbxlr.h>
#include <dcbgen.h>
C
      PARAMETER ( D4 = 4.0D0 )
C
      CHARACTER ALAB*8,BLAB*8
      DIMENSION VECBEE(MZYEE,NZ),VECBEP(MZYEP,NZ)
      DIMENSION RESULT(3,MAXQR),CMO(*),WORK(*)
C
#include <muld2h.h>
C
      CALL QENTER('LRCALC')
C
      KFREE = 1
      LFREE = LWRK
C
C     Open the property file
C
#if defined (SYS_IBM)
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
C     Zero the result vector
C
      CALL DZERO(RESULT,3*MAXQR)
C
      WRITE(LUPRI,'(A)') ' '
      CALL PRSYMB(LUPRI,'=',70,1)
      WRITE(LUPRI,'(2A)') 
     &     ' >>>>>>>>    L I N E A R   R E S P O N S E   ',
     &     'F U N C T I O N   <<<<<<<<'
      CALL PRSYMB(LUPRI,'=',70,1)
C
C     Open the files with the response vectors
C
      IF (QREE) OPEN(LUQREE,FILE='QRVEC.EE',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*MZYEE*NZ,STATUS='OLD')
      IF (QREP) OPEN(LUQREP,FILE='QRVEC.EP',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*MZYEP*NZ,STATUS='OLD')
C
C     Start loop over A operators
C
      ILRF = 0
      DO IAOP=1,NQROP
C
         KSAVE = KFREE
C
C     Initialize variables for A operator
C
         IAPT  = LQROP(IAOP,1)
         ALAB  = PRPNAM(IAPT)(1:8)
         ISYMA = IPRPSYM(IAPT)
         IFERA = JBTOF(ISYMA-1,1)
         ITIMA = IPRPTIM(IAPT)
C
         JSYMOP = ISYMA
         JOPSY  = IFERA
         JTIMOP = ITIMA
         CALL XLRPAR(WORK,KFREE,LFREE,IPRXQR)
C
         NZYAEE = 2*NZXOPE
         NZYAEP = 2*NZXOPP
C
C     Get the A operator matrix in MO basis
C
         CALL MEMGET('REAL',KA,N2ORBXQ,WORK,KFREE,LFREE)
         CALL PRPMAT(IAPT,IFERA,WORK(KA),CMO,ICMOQ,NORB,
     &               WORK,KFREE,LFREE,IPRXQR)
C
C     Get the A operator gradient:  / g*\
C                                   \ g /
C
         CALL MEMGET('REAL',KGEE,NZYAEE*NZ,WORK,KFREE,LFREE)         
         CALL MEMGET('REAL',KGEP,NZYAEP*NZ,WORK,KFREE,LFREE)         
         CALL RSP1GR(WORK(KGEE),WORK(KGEP),NZYAEE,NZYAEP,ISYMA,
     &        WORK(KA),IPRXQR,WORK(KFREE),LFREE)
C
C     Start loop over B operators
C
         DO IBOP=IAOP,NQROP
            IBPT = LQROP(IBOP,1)
            ISYMB = IPRPSYM(IBPT)
            BLAB = PRPNAM(IBPT)(1:8)
C
C     Boson symmetries of operators A and B must be equal
C
            IF (MULD2H(ISYMA,ISYMB).NE.1) GOTO 99
C
C     Get the B response vector:  / Z \
C                                 \ Y*/
C
            DO IFR=1,LQROP(IBOP,2)
               ILRF = ILRF + 1
               FREQ = QRFREQ(IFR,IBOP)
               CALL INDQR(IBPT,ABS(FREQ),IREC)
               IF (QREE) CALL READAC(LUQREE,MZYEE*NZ,VECBEE,IREC)
               IF (QREP) CALL READAC(LUQREP,MZYEP*NZ,VECBEP,IREC)
C     
C     Compute the LR function value. The response vector Nb is
C     swapped in QRSPDOT.
C     
               CALL QRSPDOT('A[1] Nb',ITIMA*D4,VECBEE,VECBEP,
     &              WORK(KGEE),WORK(KGEP),NZYAEE,NZYAEP,ISYMA,
     &              RESULT(1,ILRF),IPRXQR)
C     
C     Print out information about present LR function
C     
               WRITE(LUPRI,'(/A,2(/A,A10,I4,F10.6),/,3(/A,F20.12))')
     &              ' Linear response function in a.u.',
     &      ' A operator, boson symmetry, frequency: ',ALAB,ISYMA,-FREQ,
     &      ' B operator, boson symmetry, frequency: ',BLAB,ISYMB, FREQ,
     &              ' Value of electronic part: ',RESULT(2,ILRF),
     &              ' Value of positronic part: ',RESULT(3,ILRF),
     &              ' Value of total response : ',RESULT(1,ILRF)
C     
C     Print the Nb response vector
C     
               IF (IPRXQR.GE.10) THEN
                  CALL HEADER('VECBEE '//BLAB//
     &                 ' in LRCALC (response vector)',-1)
                  CALL PRBVEC(LUPRI,VECBEE,1,MZYEE)
                  CALL HEADER('VECBEP '//BLAB//
     &                 ' in LRCALC (response vector)',-1)
                  CALL PRBVEC(LUPRI,VECBEP,1,MZYEP)
               END IF
C
C     End loop over B frequencies
C
            END DO
C
C     End symmetry check
C
 99         CONTINUE
C
C     End loop over B operators
C     
         END DO
C
C     Print the A[1] gradient vector
C
         IF (IPRXQR.GE.10) THEN
            CALL HEADER('VECAEE '//ALAB//
     &           ' in LRCALC (property gradient)',-1)
            CALL PRBVEC(LUPRI,WORK(KGEE),1,NZYAEE)
            CALL HEADER('VECAEP '//ALAB//
     &           ' in LRCALC (property gradient)',-1)
            CALL PRBVEC(LUPRI,WORK(KGEP),1,NZYAEP)
         END IF
C
         CALL MEMREL('LRCALC',WORK,1,KSAVE,KFREE,LFREE)
C
C     End loop over A operators
C
      END DO
C
C     Close files with response vectors
C
      IF (QREE) CLOSE(LUQREE)
      IF (QREP) CLOSE(LUQREP)
C
C     Close the property file
C
      CLOSE(LU1INT,STATUS='KEEP')
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('LRCALC')
C
      RETURN
      END
      SUBROUTINE QRHYP(RESULT,VECAEE,VECBEE,VECCEE,VECAEP,VECBEP,VECCEP,
     &     CMO,WORK,LWRK)
C*****************************************************************************
C
C     ------------------------------------------------------------
C     Driver routine for computing the quadratic response function
C     ------------------------------------------------------------
C
C     Quadratic response functions are defined as:
C
C     <<A;B,C>>_{wb,wc} = 
C          Na_j(-(w1+w2))* B[2]_jk *Nc_k(w2)
C        + Na_j(-(w1+w2))* C[2]_jk *Nb_k(w1)
C        + Nb_j(w1)*( A[2]_jk + A[2]_kj )*Nc_k(w2) 
C        - Na_j(w1+w2)*( E[3]_jkl + E[3]_jlk - 
C                  w1*S[2]_jkl - w2*S[2]_jlk )*Nb_k(w1)*Nc_l(w2)
C
C     Each operator A,B, and C is assumed to span a given boson irrep
C     and be either symmetric(+) or antisymmetric(-) with respect
C     to time reversal.
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxpr.h>
#include <dcbxqr.h>
#include <dgroup.h>
#include <dcbgen.h>
C
      PARAMETER ( D0=0.0D0 )
      CHARACTER ALAB*8,BLAB*8,CLAB*8
      DIMENSION VECAEE(MZYEE,NZ),VECBEE(MZYEE,NZ),VECCEE(MZYEE,NZ),
     &          VECAEP(MZYEP,NZ),VECBEP(MZYEP,NZ),VECCEP(MZYEP,NZ)
      DIMENSION CMO(*),WORK(*)
      DIMENSION RESULT(3,MAXQR)
C
      CALL QENTER('QRHYP')
C
      KFREE = 1
      LFREE = LWRK
C
C     Open the property file
C
#if defined (SYS_IBM)
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
C     Zero the result vector
C     
      CALL DZERO(RESULT,3*MAXQR)
C
C     Print header
C
      WRITE(LUPRI,'(A)') ' '
      CALL PRSYMB(LUPRI,'=',70,1)
      WRITE(LUPRI,'(2A)') 
     &     ' >>>>>    Q U A D R A T I C   R E S P O N S E   ',
     &     'F U N C T I O N   <<<<<'
      CALL PRSYMB(LUPRI,'=',70,1)
C
C     Start loop over nontrivial QR functions
C
      DO IQRF=1,NQRHYP
C
C     Initialize variables
C
         CALL HYPINIT(IQRF,1,AFREQ,IAPT,ISYMA,IFERA,ITIMA,ALAB,
     &        NZYAEE,NZYAEP,VECAEE,VECAEP,WORK(KFREE),LFREE)
         CALL HYPINIT(IQRF,2,BFREQ,IBPT,ISYMB,IFERB,ITIMB,BLAB,
     &        NZYBEE,NZYBEP,VECBEE,VECBEP,WORK(KFREE),LFREE)
         CALL HYPINIT(IQRF,3,CFREQ,ICPT,ISYMC,IFERC,ITIMC,CLAB,
     &        NZYCEE,NZYCEP,VECCEE,VECCEP,WORK(KFREE),LFREE)
C
C     Print out information about present QR function
C
         WRITE(LUPRI,'(/A,2(I2,A),3(/A,A10,I4,F10.6))')
     &        ' Quadratic response function no. ',
     &        IQRF,' out of ',NQRHYP,' in a.u.',
     &      ' A operator, boson symmetry, frequency: ',ALAB,ISYMA,AFREQ,
     &      ' B operator, boson symmetry, frequency: ',BLAB,ISYMB,BFREQ,
     &      ' C operator, boson symmetry, frequency: ',CLAB,ISYMC,CFREQ
         CALL FLSHFO(LUPRI)
C
         CALL A2DRV(IAPT,VECBEE,VECBEP,VECCEE,VECCEP,
     &        ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &        NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,RESULT(1,IQRF),
     &        CMO,WORK(KFREE),LFREE)
C
         CALL A2DRV(IAPT,VECCEE,VECCEP,VECBEE,VECBEP,
     &        ISYMA,ISYMC,ISYMB,ITIMA,ITIMC,ITIMB,IFERA,IFERC,IFERB,
     &        NZYAEE,NZYAEP,NZYCEE,NZYCEP,NZYBEE,NZYBEP,RESULT(1,IQRF),
     &        CMO,WORK(KFREE),LFREE)
C
         CALL X2DRV(VECAEE,VECAEP,IBPT,VECCEE,VECCEP,
     &        ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &        NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,RESULT(1,IQRF),
     &        CMO,WORK(KFREE),LFREE)
C
         CALL X2DRV(VECAEE,VECAEP,ICPT,VECBEE,VECBEP,
     &        ISYMA,ISYMC,ISYMB,ITIMA,ITIMC,ITIMB,IFERA,IFERC,IFERB,
     &        NZYAEE,NZYAEP,NZYCEE,NZYCEP,NZYBEE,NZYBEP,RESULT(1,IQRF),
     &        CMO,WORK(KFREE),LFREE)
C
         CALL T3DRV(VECAEE,VECAEP,VECBEE,VECBEP,VECCEE,VECCEP,
     &        ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &        NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,BFREQ,CFREQ,
     &        RESULT(1,IQRF),CMO,WORK(KFREE),LFREE)
C
         WRITE(LUPRI,'(3(/A,F20.12))')
     &        ' Value of electronic part:',RESULT(2,IQRF),
     &        ' Value of positronic part:',RESULT(3,IQRF),
     &        ' Value of total response :',RESULT(1,IQRF)
C
C     End loop over nontrivial QR functions
C
      END DO
C
C     Close the property file
C
      CLOSE(LU1INT,STATUS='KEEP')
C
      CALL QEXIT('QRHYP')
C
      RETURN
      END
      SUBROUTINE HYPINIT(IQRF,IOP,FREQ,IPT,ISYM,IFER,ITIM,LAB,
     &     NZYEE,NZYEP,VECEE,VECEP,WORK,LWRK)
C*****************************************************************************
C
C     Initialize variables and read response vectors from file
C     for given QR function
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxpr.h>
#include <dcbxqr.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      PARAMETER ( D0=0.0D0, D1=1.0D0, DM1=-1.0D0 )
      CHARACTER LAB*8
      DIMENSION VECEE(MZYEE,NZ),VECEP(MZYEP,NZ)
      DIMENSION WORK(*)
C
      CALL QENTER('HYPINIT')
C
      KFREE = 1
      LFREE = LWRK
C
C     Open file with response vectors
C
      IF (QREE) OPEN(LUQREE,FILE='QRVEC.EE',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*MZYEE*NZ,STATUS='OLD')
      IF (QREP) OPEN(LUQREP,FILE='QRVEC.EP',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*MZYEP*NZ,STATUS='OLD')
C
C     Read VECEE and VECEP
C
      FREQ = QRFRHYP(IQRF,IOP)
      IPT    = LQRHYP(IQRF,IOP)
      ISYM   = IPRPSYM(IPT)
      IFER  =  JBTOF(ISYM-1,1)
      ITIM   = IPRPTIM(IPT)
      LAB    = PRPNAM(IPT)(1:8)
      JSYMOP = ISYM
      JOPSY  = IFER
      JTIMOP = ITIM
      CALL XLRPAR(WORK,KFREE,LFREE,IPRXQR)
      NZYEE = 2*NZXOPE
      NZYEP = 2*NZXOPP
      CALL INDQR(IPT,ABS(FREQ),IREC)
      IF (QREE) CALL READAC(LUQREE,MZYEE*NZ,VECEE,IREC)
      IF (QREP) CALL READAC(LUQREP,MZYEP*NZ,VECEP,IREC)
C
C     The linear response equations have been solved for positive
C     frequencies and we retain the correct response vector
C     with the relations:
C
C     +-------------------------+
C     | Hermitian perturbations |
C     +-------------------------+
C     N(w) = / Z \  ==>  N(-w) = / Y \
C            \ Y*/               \ Z*/
C
C     +------------------------------+
C     | Anti-hermitian perturbations |
C     +------------------------------+
C     N(w) = / Z \  ==>  N(-w) = / Y \
C            \-Y*/               \-Z*/
C
      IF (FREQ.LT.D0) THEN
         IF (ITIM.EQ.1) THEN
            FAC = D1
         ELSE IF (ITIM.EQ.-1) THEN
            FAC = DM1
         END IF
         CALL QSWAP(NZYEE,VECEE,ISYM,FAC,.TRUE.)
         CALL QSWAP(NZYEP,VECEP,ISYM,FAC,.TRUE.)
      END IF
C
C     Print
C
      IF (IPRXQR.GE.10) THEN
         IF (FREQ.LT.D0) WRITE(LUPRI,'(/A,A12,I5,2(/A)))') 
     &     ' HYPINIT: Note that response vector (label/time-symmetry):',
     &     LAB,ITIM,
     &     ' has negative frequency, so solution vector is swapped and',
     &     ' conjugated (ITIM=1) and sign-reversed (ITIM=-1).'
C
         IF (QREE) CALL HEADER('VECEE '//LAB//' in HYPINIT',-1)
         IF (QREE) CALL PRBVEC(LUPRI,VECEE,1,NZYEE)
         IF (QREP) CALL HEADER('VECEP '//LAB//' in HYPINIT',-1)
         IF (QREP) CALL PRBVEC(LUPRI,VECEP,1,NZYEP)
      END IF
C
C     Close files
C
      IF (QREE) CLOSE(LUQREE,STATUS='KEEP')
      IF (QREP) CLOSE(LUQREP,STATUS='KEEP')
C
      CALL QEXIT('HYPINIT')
C
      RETURN
      END
      SUBROUTINE A2DRV(IAPT,VECBEE,VECBEP,VECCEE,VECCEP,
     &     ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &     NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,
     &     RESULT,CMO,WORK,LWRK)
C*****************************************************************************
C
C     Driver routine for computing the contribution:
C          Nb_j(w1)* A[2]_jk *Nc_k(w2)
C
C     and adding the result to the variable:
C          RESULT(1) = total contribution
C          RESULT(2) = electronic part
C          RESULT(3) = positronic part
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxqr.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      PARAMETER ( D1 = 1.0D0, D4 = 4.0D0 )
C
      DIMENSION VECBEE(MZYEE,NZ),VECCEE(MZYEE,NZ),
     &          VECBEP(MZYEP,NZ),VECCEP(MZYEP,NZ)
      DIMENSION RESULT(3),CMO(*),WORK(*)
C
      CALL QENTER('A2DRV')
C
      KFREE = 1
      LFREE = LWRK
C
C     Allocate space for:
C     1. A operator matrix
C     2. Unpacked C vector
C     3. One-index transformed A matrix
C     4. Gradient of transformed A matrix
C
      CALL MEMGET('REAL',KA,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEE,NZYAEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEP,NZYAEP*NZ,WORK,KFREE,LFREE)
C
C     Get the A operator matrix in MO basis
C
      CALL PRPMAT(IAPT,IFERA,WORK(KA),CMO,ICMOQ,NORB,
     &            WORK,KFREE,LFREE,IPRXQR)
C
C     Unpack the C vector into a matrix 
C
      CALL GTZYMT(VECCEE,NZYCEE,VECCEP,NZYCEP,ISYMC,WORK(KC),
     &     IPRXQR,WORK(KFREE),LFREE)
C
C     One-index transform A matrix with unpacked C response vector
C
      CALL DZERO(WORK(KAC),N2ORBXQ)
      CALL OITMAT(WORK(KA),ISYMA,ITIMA,IFERA,WORK(KC),ISYMC,ITIMC,IFERC,
     &     WORK(KAC),ISYMB,ITIMA*ITIMC,IFERB,IPRXQR)
C
C     Make the gradient
C
      CALL RSP1GR(WORK(KGEE),WORK(KGEP),NZYBEE,NZYBEP,ISYMB,
     &     WORK(KAC),IPRXQR,WORK(KFREE),LFREE)
C
C     The gradient for A[2] is normally swapped at this point to conform
C     with the definition (q+ and q are interchanged in the definition
C     of A[2] and X[2] respectively). However, this swap is indirectly 
C     transfered to QRSPDOT by instead swapping the response vector Nx in 
C     that routine.
C
      CALL QRSPDOT('Nx A[2] Ny',D4,VECBEE,VECBEP,
     &     WORK(KGEE),WORK(KGEP),NZYBEE,NZYBEP,ISYMB,RESULT,IPRXQR)
C
      CALL QEXIT('A2DRV')
      RETURN
      END
      SUBROUTINE X2DRV(VECAEE,VECAEP,IBPT,VECCEE,VECCEP,
     &     ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &     NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,
     &     RESULT,CMO,WORK,LWRK)
C*****************************************************************************
C
C     Driver routine for computing the contribution:
C          Na_j(-(w1+w2))* B[2]_jk *Nc_k(w2)
C
C     and adding the result to the variable:
C          RESULT(1) = total contribution
C          RESULT(2) = electronic part
C          RESULT(3) = positronic part
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxqr.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      PARAMETER ( D1 = 1.0D0, D8 = 8.0D0 )
C
      DIMENSION VECAEE(MZYEE,NZ),VECCEE(MZYEE,NZ),
     &          VECAEP(MZYEP,NZ),VECCEP(MZYEP,NZ)
      DIMENSION RESULT(3),CMO(*),WORK(*)
C
      CALL QENTER('X2DRV')
C
      KFREE = 1
      LFREE = LWRK
C
C     Allocate space for:
C     1. B operator matrix
C     2. Unpacked C vector
C     3. One-index transformed B matrix
C     4. Gradient of transformed B matrix
C
      CALL MEMGET('REAL',KB,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEE,NZYAEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEP,NZYAEP*NZ,WORK,KFREE,LFREE)
C
C     Get the B operator matrix in MO basis
C
      CALL PRPMAT(IBPT,IFERB,WORK(KB),CMO,ICMOQ,NORB,
     &            WORK,KFREE,LFREE,IPRXQR)
C
      IF(IPRXQR.GE.10) THEN
         WRITE(LUPRI,'(A)') 
     &        'Property matrix in X2DRV'
         CALL PRQMAT(WORK(KB),NORBT,NORBT,NORBT,NORBT,NZ,
     &        IPQTOQ(1,ISYMB-1),LUPRI)
      ENDIF
C
C     Unpack the C vector into a matrix 
C
      CALL GTZYMT(VECCEE,NZYCEE,VECCEP,NZYCEP,ISYMC,WORK(KC),
     &     IPRXQR,WORK(KFREE),LFREE)
C
C     One-index transform B matrix with unpacked C response vector
C
      CALL DZERO(WORK(KBC),N2ORBXQ)
      CALL OITMAT(WORK(KB),ISYMB,ITIMB,IFERB,WORK(KC),ISYMC,ITIMC,IFERC,
     &     WORK(KBC),ISYMA,ITIMB*ITIMC,IFERA,IPRXQR)
C
C     Make the gradient
C
      CALL RSP1GR(WORK(KGEE),WORK(KGEP),NZYAEE,NZYAEP,ISYMA,
     &     WORK(KBC),IPRXQR,WORK(KFREE),LFREE)
C
C     Perform the dot product with the Na response vector. The response
C     vector Na is swapped in QRSPDOT in accordance with the definition.
C
      CALL QRSPDOT('Na X[2] Ny',D8,VECAEE,VECAEP,
     &     WORK(KGEE),WORK(KGEP),NZYAEE,NZYAEP,ISYMA,RESULT,IPRXQR)
C
      CALL QEXIT('X2DRV')
      RETURN
      END
      SUBROUTINE GTZYMT(VECEE,NZYEE,VECEP,NZYEP,ISYM,ZYMAT,IPRINT,
     &                 WORK,LWRK)
C*****************************************************************************
C
C     Unpack response vector N_{rs} into a matrix W(r,s) to be used in 
C     one-index transformations of operators according to::
C
C                             /          |  Z  |          \
C            /   \            |          |(e-p)|          | 
C            | Z |            |---------------------------| 
C        N = |   | -----> W = | -Y+ (e-p)|     | -Y+ (e-e)| 
C            | Y*|            |---------------------------| 
C            \   /            |          |  Z  |          | 
C                             \          |(e-e)|          /
C                      
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcborb.h>
#include <dcbxrs.h>
#include <dgroup.h>
#include <dcbxqr.h>
C
      DIMENSION VECEE(MZYEE,NZ),VECEP(MZYEP,NZ),ZYMAT(NORBT,NORBT,NZ)
      DIMENSION WORK(*)
C
      CALL QENTER('GTZYMT')
C
      KFREE = 1
      LFREE = LWRK
      NZEE = NZYEE/2
      NZEP = NZYEP/2
C
C     Get the index of the excitations for RPA.
C
      JSYMOP = ISYM
      JOPSY  = JBTOF(JSYMOP-1,1)
      CALL XLRPAR(WORK,KFREE,LFREE,IPRXQR)
      IF (NZXOPE.NE.NZEE .OR. NZXOPP.NE.NZEP) 
     &     CALL QUIT('GTZYMT: Inconsistent orbital rotations')
C
C     Scatter the response vector into matrix using index for RPA. 
C
      CALL DZERO(ZYMAT,N2ORBXQ)
      CALL VECTOMAT(ZYMAT,VECEE,MZYEE,NZEE,WORK(KZXOPE))
      CALL VECTOMAT(ZYMAT,VECEP,MZYEP,NZEP,WORK(KZXOPP))
C
C     Print section
C
      IF(IPRINT.GE.10) THEN
         WRITE(LUPRI,'(A)') 
     &        'Unpacked response vector in GTZYMT'
         CALL PRQMAT(ZYMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &        IPQTOQ(1,ISYM-1),LUPRI)
      ENDIF
C
      CALL QEXIT('GTZYMT')
C
      RETURN
      END
      SUBROUTINE VECTOMAT(ZYMAT,VEC,LVEC,N,INDEX)
C*****************************************************************************
C
C     Scatter the response vector into a  matrix.
C     Called from GTZYMT.
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
#include <dcborb.h>
C                          
      DIMENSION VEC(LVEC,NZ),ZYMAT(NORBT,NORBT,NZ),INDEX(2,N)
C
C     Scatter the response vector into matrix, where I1 is the 
C     inactive index and I2 is the secondary index in RPA. 
C
      IF (N.GT.0) THEN
         DO I=1,N
            I1 = INDEX(1,I)
            I2 = INDEX(2,I)
            DO IZ = 1,NZ
               ZYMAT(I2,I1,IZ) = VEC(I,IZ)
               ZYMAT(I1,I2,IZ) = -VEC(N+I,IZ)
            END DO 
         END DO
      ENDIF
C
      RETURN
      END
      SUBROUTINE OITMAT(A,ISYMA,ITIMA,IFERA,W,ISYMW,ITIMW,IFERW,
     &     C,ISYMC,ITIMC,IFERC,IPRINT)
C*****************************************************************************
C
C     One-index quaternion matrix QA with unpacked quaternion response vector 
C     in matrix QW and add the result to matrix QC. Value in {} parenthesis 
C     denote time reversal symmetries, orbitals within [] parenthesis 
C     denote Kramer pair orbitals, and indices within () denote matrix 
C     rows and columns. 
C     
C     QC(r,s){+1} = QC(r,s){+1} + 
C                   ( QW(r,p){tA}*QA(p,s){+1} - QA(r,p){tW}*QW(p,s){+1} ) , 
C
C     where
C
C     QA(r,s){tA} = A(r,s) + tA*A(r,[s]) j 
C     QW(r,s){tW} = W(r,s) + tW*W(r,[s]) j .
C
C     The structure of QA is general and the unpacked response vector
C     is of the form:
C
C           /          |  Z  |          \
C           |          |(e-p)|          | 
C           |---------------------------| 
C      QW = | -Y+ (e-p)|     | -Y+ (e-e)| 
C           |---------------------------| 
C           |          |  Z  |          | 
C           \          |(e-e)|          /
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcborb.h>
#include <dgroup.h>
C
      PARAMETER ( D1 = 1.0D0, DM1 = -1.0D0 )
      CHARACTER T*1
      DIMENSION A(NORBT,NORBT,NZ),W(NORBT,NORBT,NZ),C(NORBT,NORBT,NZ)
      DIMENSION T(-1:1)
C
#include <muld2h.h>
C
      CALL QENTER('OITMAT')
C
      IF (ISYMC.NE.MULD2H(ISYMA,ISYMW)) THEN
         CALL QUIT('Error with operator symmetry in OITMAT')
      END IF
C
      T(1)  = 'N'
C      T(-1) = 'I'
      T(-1) = 'N'
C
      DO I1 = 1,NFSYM
C
C     Start loop over fermion symmetries
C
        I2 = MOD(I1+IFERA,2) + 1
        I3 = MOD(I2+IFERW,2) + 1
C
        NORB1 = NISH(I1)
        IORB1 = IORB(I1) + NPSH(I1) + 1
        NORB2 = NORB(I2)
        IORB2 = IORB(I2) + 1
        NORB3 = NORB(I3)
        IORB3 = IORB(I3) + 1
        IF (NORB1.GT.0 .AND. NORB2.GT.0 .AND. NORB3.GT.0) THEN
           CALL QGEMM(NORB1,NORB3,NORB2,D1,
     &          'N',T(ITIMA),IPQTOQ(1,ISYMW-1),
     &          W(IORB1,IORB2,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,ISYMA-1),
     &          A(IORB2,IORB3,1),NORBT,NORBT,NZ,
     &          D1,IPQTOQ(1,ISYMC-1),C(IORB1,IORB3,1),NORBT,NORBT,NZ)
        END IF
        NORB1 = NORB(I1)
        IORB1 = IORB(I1) + 1
        NORB2 = NISH(I2)
        IORB2 = IORB(I2) + NPSH(I2) + 1
        NORB3 = NORB(I3)
        IORB3 = IORB(I3) + 1
        IF (NORB1.GT.0 .AND. NORB2.GT.0 .AND. NORB3.GT.0) THEN
           CALL QGEMM(NORB1,NORB3,NORB2,D1,
     &          'N',T(ITIMA),IPQTOQ(1,ISYMW-1),
     &          W(IORB1,IORB2,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,ISYMA-1),
     &          A(IORB2,IORB3,1),NORBT,NORBT,NZ,
     &          D1,IPQTOQ(1,ISYMC-1),C(IORB1,IORB3,1),NORBT,NORBT,NZ)
           CALL QGEMM(NORB1,NORB3,NORB2,DM1,
     &          'N',T(ITIMW),IPQTOQ(1,ISYMA-1),
     &          A(IORB1,IORB2,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,ISYMW-1),
     &          W(IORB2,IORB3,1),NORBT,NORBT,NZ,
     &          D1,IPQTOQ(1,ISYMC-1),C(IORB1,IORB3,1),NORBT,NORBT,NZ)
        END IF
C
        NORB1 = NORB(I1)
        IORB1 = IORB(I1) + 1
        NORB2 = NORB(I2)
        IORB2 = IORB(I2) + 1
        NORB3 = NISH(I3)
        IORB3 = IORB(I3) + NPSH(I3) + 1
        IF (NORB1.GT.0 .AND. NORB2.GT.0 .AND. NORB3.GT.0) THEN
           CALL QGEMM(NORB1,NORB3,NORB2,DM1,
     &          'N',T(ITIMW),IPQTOQ(1,ISYMA-1),
     &          A(IORB1,IORB2,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,ISYMW-1),
     &          W(IORB2,IORB3,1),NORBT,NORBT,NZ,
     &          D1,IPQTOQ(1,ISYMC-1),C(IORB1,IORB3,1),NORBT,NORBT,NZ)
        END IF
C
C     End loop over fermion symmetries
C
      END DO
C
C     Print section
C
      IF(IPRINT.GE.10) THEN
         WRITE(LUPRI,'(A)') 
     &        'One-index transformed matrix in OITMAT'
         CALL PRQMAT(C,NORBT,NORBT,NORBT,NORBT,NZ,
     &        IPQTOQ(1,ISYMC-1),LUPRI)
      ENDIF
C
      CALL QEXIT('OITMAT')
C
      RETURN
      END
      SUBROUTINE RSP1GR(GEE,GEP,NZYEE,NZYEP,ISYMGR,OP,IPRINT,WORK,LWRK)
C*****************************************************************************
C
C     Compute the gradent of a one-electron operator in MO basis
C     according to the formula.
C
C       / g*_{xi} \   /  <0|[q,B]|0>  \   / B(x,i) \
C       |         | = |               | = |        |
C       \ g_{xi}  /   \ -<0|[q+,B]|0> /   \ B(i,x) /
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      DIMENSION GEE(NZYEE,NZ),GEP(NZYEP,NZ),OP(NORBT,NORBT,NZ),WORK(*)
C
      CALL QENTER('RSP1GR')
C
      KFREE = 1 
      LFREE = LWRK
      NZEE = NZYEE/2
      NZEP = NZYEP/2
C
C     Get the index to the excitations.
C
      JSYMOP = ISYMGR
      JOPSY  = JBTOF(JSYMOP-1,1)
      CALL XLRPAR(WORK,KFREE,LFREE,IPRXQR)
      IF (NZXOPE.NE.NZEE .OR. NZXOPP.NE.NZEP) 
     &     CALL QUIT('RSP1GR: Inconsistent orbital rotations')
C
C     Scatter the operator into the gradient using index.
C
      CALL OPTOGR(OP,GEE,NZEE,WORK(KZXOPE))
      CALL OPTOGR(OP,GEP,NZEP,WORK(KZXOPP))
C
C     Print section
C
      IF (IPRINT.GE.10) THEN
         WRITE(LUPRI,'(/A)') 'Gradient (e-e) in RSP1GR'
         CALL PRQMAT(GEE,NZYEE,1,NZYEE,1,NZ,IPQTOQ(1,ISYMGR-1),LUPRI)
         WRITE(LUPRI,'(/A)') 'Gradient (e-p) in RSP1GR'
         CALL PRQMAT(GEP,NZYEP,1,NZYEP,1,NZ,IPQTOQ(1,ISYMGR-1),LUPRI)
      ENDIF
C
      CALL QEXIT('RSP1GR')
C
      RETURN
      END
      SUBROUTINE OPTOGR(OP,G,N,INDEX)
C*****************************************************************************
C
C     Scatter the operator matrix into the gradient for RPA.
C     Called from RSP1GR.
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <dgroup.h>
#include <dcborb.h>
C                          
      DIMENSION OP(NORBT,NORBT,NZ),G(2*N,NZ),INDEX(2,N)
C
C     Compute the gradient, where I1 is the inactive index 
C     and I2 is the secondary index in RPA. 
C
      IF (N.GT.0) THEN
         DO I=1,N
            I1 = INDEX(1,I)
            I2 = INDEX(2,I)
            DO IZ = 1,NZ
               G(I,IZ) = OP(I2,I1,IZ)
               G(N+I,IZ) = OP(I1,I2,IZ)
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
      SUBROUTINE T3DRV(VECAEE,VECAEP,VECBEE,VECBEP,VECCEE,VECCEP,
     &     ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &     NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,BFREQ,CFREQ,
     &     RESULT,CMO,WORK,LWRK)
C*****************************************************************************
C
C     Driver routine for computing the contribution:
C          Na_j(-(w1+w2))* (E[3]_jkl + E[3]_jlk
C                          - wb*S[3]_jkl - wc*S[3]_jlk ) *Nb_k(w1)*Nc_l(w2)
C
C     and adding the result to the variable:
C          RESULT(1) = total contribution
C          RESULT(2) = electronic part
C          RESULT(3) = positronic part
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dcbxqr.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      PARAMETER ( D0 = 0.0D0 )
C
      DIMENSION VECAEE(MZYEE,NZ),VECBEE(MZYEE,NZ),VECCEE(MZYEE,NZ),
     &          VECAEP(MZYEP,NZ),VECBEP(MZYEP,NZ),VECCEP(MZYEP,NZ)
      DIMENSION RESULT(3),CMO(*),WORK(*)
C
      CALL QENTER('T3DRV')
C
      KFREE = 1 
      LFREE = LWRK
C
      CALL  E3DRV(VECAEE,VECAEP,VECBEE,VECBEP,VECCEE,VECCEP,
     &     ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &     NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,
     &     RESULT,CMO,WORK(KFREE),LFREE)
C
      IF (BFREQ.NE.D0) THEN
         CALL  S3DRV(VECAEE,VECAEP,VECBEE,VECBEP,VECCEE,VECCEP,
     &        ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &        NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,BFREQ,
     &        RESULT,CMO,WORK(KFREE),LFREE)
      END IF
C
      IF (CFREQ.NE.D0) THEN
         CALL  S3DRV(VECAEE,VECAEP,VECCEE,VECCEP,VECBEE,VECBEP,
     &        ISYMA,ISYMC,ISYMB,ITIMA,ITIMC,ITIMB,IFERA,IFERC,IFERB,
     &        NZYAEE,NZYAEP,NZYCEE,NZYCEP,NZYBEE,NZYBEP,CFREQ,
     &        RESULT,CMO,WORK(KFREE),LFREE)
      END IF
C
      CALL QEXIT('T3DRV')
C     
      RETURN
      END
      SUBROUTINE E3DRV(VECAEE,VECAEP,VECBEE,VECBEP,VECCEE,VECCEP,
     &     ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &     NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,
     &     RESULT,CMO,WORK,LWRK)
C*****************************************************************************
C
C     Driver routine for computing the contribution:
C          Na_j(-(w1+w2))* (E[3]_jkl + E[3]_jlk ) *Nb_k(w1)*Nc_l(w2)
C
C     and adding the result to the variable:
C          RESULT(1) = total contribution
C          RESULT(2) = electronic part
C          RESULT(3) = positronic part
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dcbxqr.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0 )
C
      DIMENSION VECAEE(MZYEE,NZ),VECBEE(MZYEE,NZ),VECCEE(MZYEE,NZ),
     &          VECAEP(MZYEP,NZ),VECBEP(MZYEP,NZ),VECCEP(MZYEP,NZ)
      DIMENSION ISYMOP(3),IHRMOP(3),IFCKOP(3)
      DIMENSION RESULT(3),CMO(*),WORK(*)
C
      CALL QENTER('E3DRV')
C
      KFREE = 1 
      LFREE = LWRK
C
      CALL MEMGET('REAL',KB,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KC,N2ORBXQ,WORK,KFREE,LFREE)
C
      CALL MEMGET('REAL',KF1AO,3*N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KD1AO,3*N2BBASXQ,WORK,KFREE,LFREE)
C
      KF2AO = KF1AO + N2BBASXQ
      KF3AO = KF2AO + N2BBASXQ
      KD2AO = KD1AO + N2BBASXQ
      KD3AO = KD2AO + N2BBASXQ
C
C     Unpack the response vectors into matrices 
C
      IF (IPRXQR.GE.10) THEN
         WRITE(LUPRI,'(/A)') 
     &        'CMO in e3drv'
         CALL PRQMAT(CMO,NTBAS(0),NORBT,NTBAS(0),NORBT,NZ,
     &        IPQTOQ(1,0),LUPRI)
      END IF
C     
      CALL GTZYMT(VECBEE,NZYBEE,VECBEP,NZYBEP,ISYMB,WORK(KB),
     &     IPRXQR,WORK(KFREE),LFREE)
      CALL GTZYMT(VECCEE,NZYCEE,VECCEP,NZYCEP,ISYMC,WORK(KC),
     &     0,WORK(KFREE),LFREE)
C
C     Construct the modified density matrices
C
      CALL DZERO(WORK(KD1AO),3*N2BBASXQ)
      CALL QDENS1(ISYMB,ITIMB,IFERB,WORK(KB),WORK(KD1AO),
     &     IFCKOP(1),CMO,WORK(KFREE),LFREE)
      CALL QDENS1(ISYMC,ITIMC,IFERC,WORK(KC),WORK(KD2AO),
     &     IFCKOP(2),CMO,WORK(KFREE),LFREE)
      CALL QDENS2(ISYMB,ITIMB,IFERB,WORK(KB),ISYMC,ITIMC,IFERC,WORK(KC),
     &     WORK(KD3AO),ISYMA,ITIMB*ITIMC,IFERA,IFCKOP(3),
     &     CMO,WORK(KFREE),LFREE)
      CALL QDENS2(ISYMC,ITIMC,IFERC,WORK(KC),ISYMB,ITIMB,IFERB,WORK(KB),
     &     WORK(KD3AO),ISYMA,ITIMB*ITIMC,IFERA,IFCKOP(3),
     &     CMO,WORK(KFREE),LFREE)
C
C     Construct the Fock matrices in AO basis. 
C
      CALL E3FCK(ISYMOP,IHRMOP,IFCKOP,ISYMB,ISYMC,ISYMA,
     &     ITIMB,ITIMC,ITIMB*ITIMC,IFERB,IFERC,IFERA)
      CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,WORK(KF1AO),WORK(KD1AO),3,
     &     IDUM,INTFLG,IPRXQR,WORK(KFREE),LFREE)
C
      IF (IPRXQR.GE.10) THEN
         WRITE(LUPRI,'(A)') 
     &        'F1AO in e3drv'
         CALL PRQMAT(work(kf1ao),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &        IPQTOQ(1,isymb-1),LUPRI)
         WRITE(LUPRI,'(A)') 
     &        'F12AO in e3drv'
         CALL PRQMAT(work(kf3ao),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &        IPQTOQ(1,isyma-1),LUPRI)
      ENDIF
C
C     Reuse the memory of the density matrices for the 
C     Fock matrices in MO basis.
C
      CALL MEMREL('E3DRV',WORK,1,KD1AO,KFREE,LFREE)
C
      CALL MEMGET('REAL',KF3MO,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KF2MO,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KF1MO,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFMO,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Read in the Fock matrix in MO basis from file, and transform the 
C     modified Fock matrices constructed above from AO to MO basis
C
      CALL GETFCK(WORK(KFMO),IPRXQR,WORK,KFREE,LFREE)
      CALL FAOMO(WORK(KF1MO),WORK(KF1AO),ISYMB,ITIMB,IFERB,IPRXQR,
     &     CMO,WORK(KFREE),LFREE)
      CALL FAOMO(WORK(KF2MO),WORK(KF2AO),ISYMC,ITIMC,IFERC,IPRXQR,
     &     CMO,WORK(KFREE),LFREE)
      CALL FAOMO(WORK(KF3MO),WORK(KF3AO),ISYMA,ITIMB*ITIMC,IFERA,IPRXQR,
     &     CMO,WORK(KFREE),LFREE)

C
C     Compute the final doubly one-index transformed Fock matrix.
C     
      CALL DSCAL(N2ORBXQ,D2,WORK(KB),1)
      CALL DSCAL(N2ORBXQ,D2,WORK(KC),1)
C
C     [Wb,F] ---> F1
C     [Wc,F] ---> F2
C
      CALL OITMAT(WORK(KFMO),1,1,1,WORK(KB),ISYMB,ITIMB,IFERB,
     &     WORK(KF1MO),ISYMB,ITIMB,IFERB,IPRXQR)
      CALL OITMAT(WORK(KFMO),1,1,1,WORK(KC),ISYMC,ITIMC,IFERC,
     &     WORK(KF2MO),ISYMC,ITIMC,IFERC,IPRXQR)
C
C     [Wc, F1+[Wb,F]] ---> F3
C     [Wb, F2+[Wc,F]] ---> F3
C
      CALL OITMAT(WORK(KF1MO),ISYMB,ITIMB,IFERB,
     &     WORK(KC),ISYMC,ITIMC,IFERC,
     &     WORK(KF3MO),ISYMA,ITIMB*ITIMC,IFERA,IPRXQR)
      CALL OITMAT(WORK(KF2MO),ISYMC,ITIMC,IFERC,
     &     WORK(KB),ISYMB,ITIMB,IFERB,
     &     WORK(KF3MO),ISYMA,ITIMB*ITIMC,IFERA,IPRXQR)     
C
C     Reuse the memory of the intermediate Fock
C     matrices in MO basis for the gradient.  
C
      CALL MEMREL('E3DRV',WORK,1,KF2MO,KFREE,LFREE)
C
      CALL MEMGET('REAL',KGEE,NZYAEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEP,NZYAEP*NZ,WORK,KFREE,LFREE)
C
C     Make the gradient.
C
      CALL RSP1GR(WORK(KGEE),WORK(KGEP),NZYAEE,NZYAEP,ISYMA,
     &     WORK(KF3MO),IPRXQR,WORK(KFREE),LFREE)
C
C     Perform the dot product with the Na response vector. The response
C     vector Na is swapped in QRSPDOT in accordance with the definition.
C
      CALL QRSPDOT('Na E[3] Nb Nc',D2,VECAEE,VECAEP,
     &     WORK(KGEE),WORK(KGEP),NZYAEE,NZYAEP,ISYMA,RESULT,IPRXQR)
C
      CALL QEXIT('E3DRV')
C     
      RETURN
      END
      SUBROUTINE S3DRV(VECAEE,VECAEP,VECBEE,VECBEP,VECCEE,VECCEP,
     &     ISYMA,ISYMB,ISYMC,ITIMA,ITIMB,ITIMC,IFERA,IFERB,IFERC,
     &     NZYAEE,NZYAEP,NZYBEE,NZYBEP,NZYCEE,NZYCEP,BFREQ,
     &     RESULT,CMO,WORK,LWRK)
C*****************************************************************************
C
C     Driver routine for computing the contribution:
C          Na_j(-(w1+w2))* (w1*S[3]_jkl + w2*S[3]_jlk ) *Nb_k(w1)*Nc_l(w2)
C
C     and adding the result to the variable:
C          RESULT(1) = total contribution
C          RESULT(2) = electronic part
C          RESULT(3) = positronic part
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dcbxqr.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      PARAMETER ( D2 = 2.0D0 )
C
      DIMENSION VECAEE(MZYEE,NZ),VECBEE(MZYEE,NZ),VECCEE(MZYEE,NZ),
     &          VECAEP(MZYEP,NZ),VECBEP(MZYEP,NZ),VECCEP(MZYEP,NZ)
      DIMENSION RESULT(3),CMO(*),WORK(*)
C
      CALL QENTER('S3DRV')
C
      KFREE = 1 
      LFREE = LWRK
C
      CALL MEMGET('REAL',KB,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBC,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Unpack the response vectors into matrices 
C
      CALL GTZYMT(VECBEE,NZYBEE,VECBEP,NZYBEP,ISYMB,WORK(KB),
     &     IPRXQR,WORK(KFREE),LFREE)
      CALL GTZYMT(VECCEE,NZYCEE,VECCEP,NZYCEP,ISYMC,WORK(KC),
     &     IPRXQR,WORK(KFREE),LFREE)
C
C     One-index transform B with C
C
      CALL DZERO(WORK(KBC),N2ORBXQ)
      CALL OITMAT(WORK(KB),ISYMB,ITIMB,IFERB,WORK(KC),ISYMC,ITIMC,IFERC,
     &     WORK(KBC),ISYMA,ITIMB*ITIMC,IFERA,IPRXQR)
C
C     Make the gradient.
C
      CALL MEMGET('REAL',KGEE,NZYAEE*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEP,NZYAEP*NZ,WORK,KFREE,LFREE)
C
      CALL RSP1GR(WORK(KGEE),WORK(KGEP),NZYAEE,NZYAEP,ISYMA,
     &     WORK(KBC),IPRXQR,WORK(KFREE),LFREE)
C
C     Perform the dot product with the Na response vector. The response
C     vector Na is swapped in QRSPDOT in accordance with the definition.
C
      CALL QRSPDOT('Na S[3] Nb Nc',BFREQ,VECAEE,VECAEP,
     &     WORK(KGEE),WORK(KGEP),NZYAEE,NZYAEP,ISYMA,RESULT,IPRXQR)
C
      CALL QEXIT('S3DRV')
C     
      RETURN
      END
      SUBROUTINE QDENS1(ISYM1,ITIM1,IFER1,W1,DEN,IFCK,
     &     CMO,WORK,LWRK)
C*****************************************************************************
C
C     Construct the modified density matrix
C
C     IFCK -- used in TWOFCK call, determines where to add the Coulomb 
C             integrals. For time-symmetric operators the Coloumb
C             contribution is real, and for time-antisymmetric operators
C             the contribution is i-imaginary.
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbqal.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
#include <dcbxqr.h>
C
      PARAMETER ( D1 = 1.0D0, DM1 = -1.0D0 )
C
      DIMENSION DEN(N2BBASX,NZ),CMO(*),W1(NORBT,NORBT,NZ)
      DIMENSION WORK(*)
C
      CALL QENTER('QDENS1')
C
      KFREE = 1 
      LFREE = LWRK
C
      CALL MEMGET('REAL',KTMPAO,N2BBASXQ,WORK,KFREE,LFREE)
C
C     D_{ab} = C_{ai} W1_{iu} C+_{bu} - C_{au} W1_{ui} C+_{bi}
C
      CALL QDEN1(1,D1,ISYM1,ITIM1,IFER1,W1,DEN,WORK(KTMPAO),
     &     CMO,WORK(KFREE),LWRK)
      CALL QDEN1(2,DM1,ISYM1,ITIM1,IFER1,W1,DEN,WORK(KTMPAO),
     &     CMO,WORK(KFREE),LWRK)
C
C     Print
C
      IF (IPRXQR.GE.10) THEN
        CALL HEADER('QDENS1: Density matrix',-1)
        CALL PRQMAT(DEN,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,ISYM1-1),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
C
      CALL QEXIT('QDENS1')
C     
      RETURN
      END
      SUBROUTINE QDEN1(JINDX,FAC,ISYM1,ITIM1,IFER1,W1,DEN,DTMP,
     &     CMO,WORK,LWRK)
C*****************************************************************************
C
C     Construct the modified density matrix
C
C     If JINDX=1 then D_{ab} = D_{ab} + FAC * C_{ai} W1_{iu} C+_{bu}
C     If JINDX=2 then D_{ab} = D_{ab} + FAC * C_{au} W1_{ui} C+_{bi}
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbqal.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
#include <dcbxqr.h>
C
      PARAMETER ( D0 = 0.0D0, DM1 = -1.0D0 )
C
      CHARACTER T*1
      DIMENSION DEN(N2BBASX,NZ),DTMP(N2BBASX,NZ),W1(NORBT,NORBT,NZ)
      DIMENSION CMO(*),WORK(*)
C
      CALL QENTER('QDEN1')
C
      KFREE = 1 
      LFREE = LWRK
C
      IF (ITIM1.EQ.1) THEN
         T='S'
      ELSEIF (ITIM1.EQ.-1) THEN
C         T='A'
         T='S'
      ELSE
         CALL QUIT('Unknown time-reversal symmetry of matrix in QDEN1')
      END IF
C
C     Zero temporary storage DTMP
C
      CALL DZERO(DTMP,N2BBASXQ)
C
C     Generate modified density matrix as transformation MO to AO
C
      DO I1 = 1,NFSYM
C
C     Start loop over fermion symmetries
C
        I2 = MOD(I1+IFER1,2) + 1
        IF (JINDX.EQ.1) THEN
           NORB1 = NISH(I1)
           NORB2 = NORB(I2)
           IORB1 = IORB(I1) + NPSH(I1) + 1
           IORB2 = IORB(I2) + 1
           ICMO1 = ICMOQ(I1) + NPSH(I1)*NFBAS(I1,0) + 1
           ICMO2 = ICMOQ(I2) + 1
        ELSEIF (JINDX.EQ.2) THEN
           NORB1 = NORB(I1)
           NORB2 = NISH(I2)
           IORB1 = IORB(I1) + 1
           IORB2 = IORB(I2) + NPSH(I2) + 1
           ICMO1 = ICMOQ(I1) + 1
           ICMO2 = ICMOQ(I2) + NPSH(I2)*NFBAS(I2,0) + 1
        ELSE
           CALL QUIT('Unknown JINDX in QDENS1')
        END IF
C
        IF (NORB1.GT.0 .AND. NORB2.GT.0) THEN
           CALL QTRANS('MOAO',T,D1, 
     &          NFBAS(I1,0),NFBAS(I2,0),NORB1,NORB2,
     &          DTMP(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),
     &          NZ,IPQTOQ(1,ISYM1-1),
     &          W1(IORB1,IORB2,1),NORBT,NORBT,NZ,IPQTOQ(1,ISYM1-1),
     &          CMO(ICMO1),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &          CMO(ICMO2),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &          WORK(KFREE),LFREE,IPRXQR)
        END IF
C
C     End loop over fermion symmetries
C
      END DO
C
C     Add result to DEN
C
      CALL QAXPY(N2BBASX,'N','N',FAC,1,1,
     &     DTMP,N2BBASX,1,IPQTOQ(1,ISYM1-1),NZ,
     &     DEN,N2BBASX,1,IPQTOQ(1,ISYM1-1),NZ)
C
      CALL QEXIT('QDEN1')
C     
      RETURN
      END
      SUBROUTINE QDENS2(ISYM1,ITIM1,IFER1,W1,ISYM2,ITIM2,IFER2,W2,
     &     DEN,ISYM3,ITIM3,IFER3,IFCK,CMO,WORK,LWRK)
C*****************************************************************************
C
C     Construct the modified density matrix
C
C     IFCK -- used in TWOFCK call, determines where to add the Coulomb 
C             integrals. For time-symmetric operators the Coloumb
C             contribution is real, and for time-antisymmetric operators
C             the contribution is i-imaginary.
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbqal.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
#include <dcbxqr.h>
C
      PARAMETER ( D1 = 1.0D0, DM1 = -1.0D0 )
C
      DIMENSION DEN(N2BBASX,NZ),W1(NORBT,NORBT,NZ),W2(NORBT,NORBT,NZ)
      DIMENSION CMO(*),WORK(*)
C
#include <muld2h.h>
C
      CALL QENTER('QDENS2')
C
      KFREE = 1 
      LFREE = LWRK
C
      IF (ITIM3.EQ.1) THEN
         IFCK = 1
      ELSEIF (ITIM3.EQ.-1) THEN
         IFCK = 2
      ELSE
         CALL QUIT('Unknown time-reversal symmetry of matrix in QDENS1')
      END IF
C
      CALL MEMGET('REAL',KTMPAO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMPMO,N2ORBXQ,WORK,KFREE,LFREE)
C
C     D_{ab} = C_{ai} W1_{iu} W2_{uv} C+_{bv} 
C            + C_{av} W2_{vu} W1_{ui} C+_{bi} 
C            - C_{av} W2_{vi} W1_{iu} C+_{bu} 
C            - C_{au} W1_{ui} W2_{iv} C+_{bv} 
C
      CALL QDEN2(1,D1,ISYM1,ITIM1,IFER1,W1,DEN,
     &     ISYM2,ITIM2,IFER2,W2,WORK(KTMPAO),
     &     ISYM3,ITIM3,IFER3,WORK(KTMPMO),
     &     CMO,WORK(KFREE),LWRK)
      CALL QDEN2(3,D1,ISYM2,ITIM2,IFER2,W2,DEN,
     &     ISYM1,ITIM1,IFER1,W1,WORK(KTMPAO),
     &     ISYM3,ITIM3,IFER3,WORK(KTMPMO),
     &     CMO,WORK(KFREE),LWRK)
      CALL QDEN2(2,DM1,ISYM2,ITIM2,IFER2,W2,DEN,
     &     ISYM1,ITIM1,IFER1,W1,WORK(KTMPAO),
     &     ISYM3,ITIM3,IFER3,WORK(KTMPMO),
     &     CMO,WORK(KFREE),LWRK)
      CALL QDEN2(2,DM1,ISYM1,ITIM1,IFER1,W1,DEN,
     &     ISYM2,ITIM2,IFER2,W2,WORK(KTMPAO),
     &     ISYM3,ITIM3,IFER3,WORK(KTMPMO),
     &     CMO,WORK(KFREE),LWRK)
C
C     Print
C
      IF (IPRXQR.GE.10) THEN
        CALL HEADER('QDENS2: Density matrix',-1)
        CALL PRQMAT(DEN,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,ISYM3-1),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
C
      CALL QEXIT('QDENS2')
C     
      RETURN
      END
      SUBROUTINE QDEN2(JINDX,FAC,ISYM1,ITIM1,IFER1,W1,DEN,
     &     ISYM2,ITIM2,IFER2,W2,DTMP,
     &     ISYM3,ITIM3,IFER3,WTMP,
     &     CMO,WORK,LWRK)
C*****************************************************************************
C
C     Construct the modified density matrix
C
C     If JINDX=1 then D_{ab} = D_{ab} + FAC * C_{ai} W1_{iu} W2_{uv} C+_{bv}
C     If JINDX=2 then D_{ab} = D_{ab} + FAC * C_{au} W1_{ui} W2_{iv} C+_{bv}
C     If JINDX=3 then D_{ab} = D_{ab} + FAC * C_{au} W1_{uv} W2_{vi} C+_{bi}
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbqal.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
#include <dcbxqr.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0 )
C
      CHARACTER T*1
      DIMENSION DEN(N2BBASX,NZ),DTMP(N2BBASX,NZ),WTMP(NORBT,NORBT,NZ)
      DIMENSION W1(NORBT,NORBT,NZ),W2(NORBT,NORBT,NZ)
      DIMENSION T(-1:1),IQBF(4),CMO(*),WORK(*)
C
      CALL QENTER('QDEN2')
C
      KFREE = 1 
      LFREE = LWRK
C
C     Zero temporary storage DTMP
C
      CALL DZERO(DTMP,N2BBASXQ)
C
C     Generate modified density matrix as transformation MO to AO
C
      DO I1 = 1,NFSYM
C
C     Start loop over fermion symmetries
C
        I2 = MOD(I1+IFER1,2) + 1
        I3 = MOD(I2+IFER2,2) + 1
        NORB1 = NORB(I1)
        NORB2 = NORB(I2)
        NORB3 = NORB(I3)
        IORB1 = IORB(I1) + 1
        IORB2 = IORB(I2) + 1
        IORB3 = IORB(I3) + 1
        ICMO1 = ICMOQ(I1) + 1
        ICMO2 = ICMOQ(I2) + 1
        ICMO3 = ICMOQ(I3) + 1
        IF (JINDX.EQ.1) THEN
           NORB1 = NISH(I1)
           IORB1 = IORB(I1)  + NPSH(I1) + 1
           ICMO1 = ICMOQ(I1) + NPSH(I1)*NFBAS(I1,0) + 1
        ELSEIF (JINDX.EQ.2) THEN
           NORB2 = NISH(I2)
           IORB2 = IORB(I2)  + NPSH(I2) + 1
           ICMO2 = ICMOQ(I2) + NPSH(I2)*NFBAS(I2,0) + 1
        ELSEIF (JINDX.EQ.3) THEN
           NORB3 = NISH(I3)
           IORB3 = IORB(I3)  + NPSH(I3) + 1
           ICMO3 = ICMOQ(I3) + NPSH(I3)*NFBAS(I3,0) + 1
        ELSE
           CALL QUIT('Unknown JINDX in QDENS1')
        END IF
C
        IF (NORB1.GT.0 .AND. NORB2.GT.0 .AND. NORB3.GT.0) THEN
C
           T(1)  = 'N'
C           T(-1) = 'I'
           T(-1) = 'N'
C
           CALL DZERO(WTMP,N2ORBXQ)
           CALL QGEMM(NORB1,NORB3,NORB2,D1,
     &          'N',T(ITIM2),IPQTOQ(1,ISYM1-1),
     &          W1(IORB1,IORB2,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,ISYM2-1),
     &          W2(IORB2,IORB3,1),NORBT,NORBT,NZ,
     &          D0,IPQTOQ(1,ISYM3-1),
     &          WTMP(IORB1,IORB3,1),NORBT,NORBT,NZ)
C           
           T(1)  = 'S'
C           T(-1) = 'A'
           T(-1) = 'S'
C
           CALL QTRANS('MOAO',T(ITIM3),D1, 
     &          NFBAS(I1,0),NFBAS(I3,0),NORB1,NORB3,
     &          DTMP(I2BASX(I1,I3)+1,1),NTBAS(0),NTBAS(0),
     &          NZ,IPQTOQ(1,ISYM3-1),
     &          WTMP(IORB1,IORB3,1),NORBT,NORBT,NZ,IPQTOQ(1,ISYM3-1),
     &          CMO(ICMO1),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &          CMO(ICMO3),NFBAS(I3,0),NORB(I3),NZ,IPQTOQ(1,0),
     &          WORK(KFREE),LFREE,IPRXQR)
C
        END IF
C
C     End loop over fermion symmetries
C
      END DO
C
C     Add result to DEN
C
      CALL QAXPY(N2BBASX,'N','N',FAC,1,1,
     &     DTMP,N2BBASX,1,IPQTOQ(1,ISYM3-1),NZ,
     &     DEN,N2BBASX,1,IPQTOQ(1,ISYM3-1),NZ)
C
      CALL QEXIT('QDEN2')
C     
      RETURN
      END
      SUBROUTINE E3FCK(ISYMOP,IHRMOP,IFCKOP,ISYM1,ISYM2,ISYM3,
     &     ITIM1,ITIM2,ITIM3,IFER1,IFER2,IFER3)
C*****************************************************************************
C
C     Set up the integer arrays that are used by TWOFCK in the 
C     construction of Fock matrices.
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxrs.h>
#include <dcbxqr.h>
C
      DIMENSION ISYMOP(3), IHRMOP(3), IFCKOP(3)
C
      CALL QENTER('E3FCK')
C
      ISYMOP(1) = ISYM1
      ISYMOP(2) = ISYM2
      ISYMOP(3) = ISYM3
      IHRMOP(1) = 0
      IHRMOP(2) = 0
      IHRMOP(3) = 0
      IF (ITIM1.EQ.1) THEN
         IFCKOP(1) = 1
      ELSE
         IFCKOP(1) = 2
      END IF
      IF (ITIM2.EQ.1) THEN
         IFCKOP(2) = 1
      ELSE
         IFCKOP(2) = 2
      END IF
      IF (ITIM3.EQ.1) THEN
         IFCKOP(3) = 1
      ELSE
         IFCKOP(3) = 2
      END IF
C
C     Detrmine which integrals to include. Default is to include
C     all types, but other alternatives can be specified in the 
C     input file. 
C
      INTFLG = INTXQR
C
      IF (INTFLG.LE.0) THEN
         WRITE(LUPRI,'(2(/A),/)') ' *** WARNING ****',
     &        ' No integrals included in Fock matrix construction'
      END IF
C
      CALL QEXIT('E3FCK')
C     
      RETURN
      END
      SUBROUTINE FAOMO(FMO,FAO,ISYM,ITIM,IFER,IPRINT,CMO,WORK,LWRK)
C*****************************************************************************
C
C     Transform Fock matrix from AO to MO basis using QTRANS
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbqal.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
#include <dcbxqr.h>
C
      PARAMETER ( D0 = 0.0D0 )
C
      CHARACTER T*1
C
      DIMENSION FAO(N2BBASXQ),FMO(N2ORBXQ)
      DIMENSION IQBF(4),T(-1:1),CMO(*),WORK(*)
C
      CALL QENTER('FAOMO')
C
      KFREE = 1 
      LFREE = LWRK
C
      T(1)='S'
      T(-1)='A'
C
      CALL DZERO(FMO,N2ORBXQ)
C
      DO I1=1,NFSYM
         I2 = MOD(I1+IFER,2) + 1
         IF(NORB(I1).GT.0 .AND. NORB(I2).GT.0) THEN
            CALL QTRANS('AOMO',T(ITIM),
     &           D0,NFBAS(I1,0),NFBAS(I2,0),NORB(I1),NORB(I2),
     &           FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &           IPQTOQ(1,ISYM-1),
     &           FMO(I2ORBX(I1,I2)+1),NORBT,NORBT,NZ,
     &           IPQTOQ(1,ISYM-1),
     &           CMO(ICMOQ(I1)+1),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &           CMO(ICMOQ(I2)+1),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
         END IF
      END DO
C
C     Print
C
      IF (IPRINT.GE.10) THEN
         CALL HEADER('FAOMO: Fock matrix in MO-basis',-1)
         CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &        IPQTOQ(1,ISYM-1),LUPRI)
      ENDIF
C
      CALL QEXIT('FAOMO')
C     
      RETURN
      END
      SUBROUTINE QRSPDOT(TEXT,FAC,VECAEE,VECAEP,GEE,GEP,
     &     NZYAEE,NZYAEP,ISYMGR,RESULT,IPRINT)
C*****************************************************************************
C
C     Perform the dot product with the response vector and the gradient
C     and add the result times a factor to the variable RESULT. The dot
C     product is performed according to the formula:
C
C             RESULT = -(Y Z*)/ g*\ = -(Z Y*)/ g \
C                             \ g /          \ g*/
C
C     The former is the explicit form of the X[n] and E[n] type of 
C     contributions, whereas the latter is the explicit form of the 
C     A[n] type of contributions. 
C
C     Input to QRSPDOT:
C     -----------------
C     FAC is a factor scaling the result
C
C     GEE and GEP are gradients for (e-e) and (e-p) rotations of 
C     the form: / g*\
C               \ g /
C
C     VECAEE and VECAEP are response vectors  for (e-e) and (e-p) rotations
C     of the form: / Z \
C                  \ Y*/
C
C     Output returned by QRSPDOT:
C     ---------------------------
C     RESULT(1) = total contribution
C     RESULT(2) = electronic part
C     RESULT(3) = positronic part
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxqr.h>
#include <dgroup.h>
C
      PARAMETER ( D1 = 1.0D0, DM1 = -1.0D0 )
C
      CHARACTER*(*) TEXT
      DIMENSION VECAEE(MZYEE,NZ),VECAEP(MZYEP,NZ),
     &          GEE(NZYAEE,NZ),GEP(NZYAEP,NZ)
      DIMENSION RESULT(3),TEMP(3,3)
C
      CALL QENTER('QRSPDOT')
C
      NZAEE = NZYAEE/2
      NZAEP = NZYAEP/2
C
C     TEMP stores the contributions individually
C     ------------------------------------------
C     E+P:Z+Y      E+P:Z     E+P:Y
C       E:Z+Y        E:Z       E:Y
C       P:Z+Y        P:Z       P:Y
C
      CALL DZERO(TEMP,3*3)
      DO IZ=1,NZ
         IF (IPQTOQ(IZ,ISYMGR-1).EQ.1) THEN
            SIGN = DM1
         ELSE 
            SIGN = D1
         END IF
         TEMP(2,2) = TEMP(2,2) +
     &       SIGN*FAC*DDOT(NZAEE,VECAEE(1,IZ),1,GEE(1+NZAEE,IZ),1)
         TEMP(2,3) = TEMP(2,3) + 
     &       SIGN*FAC*DDOT(NZAEE,VECAEE(1+NZAEE,IZ),1,GEE(1,IZ),1)
         TEMP(3,2) = TEMP(3,2) +
     &       SIGN*FAC*DDOT(NZAEP,VECAEP(1,IZ),1,GEP(1+NZAEP,IZ),1)
         TEMP(3,3) = TEMP(3,3) +
     &       SIGN*FAC*DDOT(NZAEP,VECAEP(1+NZAEP,IZ),1,GEP(1,IZ),1)
      END DO
C
      TEMP(1,2) = TEMP(2,2)+TEMP(3,2)
      TEMP(1,3) = TEMP(2,3)+TEMP(3,3)
      TEMP(2,1) = TEMP(2,2)+TEMP(2,3)
      TEMP(3,1) = TEMP(3,2)+TEMP(3,3)
      TEMP(1,1) = TEMP(1,2)+TEMP(1,3)
C
      CALL DAXPY(3,D1,TEMP,1,RESULT,1)
C
C     Print
C
      IF (IPRINT.GE.0) THEN
         WRITE(LUPRI,'(/A16,A22,3F16.8,2(/A38,3F16.8))')
     &        TEXT,'contribution (T,Z,Y):',(TEMP(1,I),I=1,3),
     &        'electronic (T,Z,Y):',(TEMP(2,I),I=1,3),
     &        'positronic (T,Z,Y):',(TEMP(3,I),I=1,3)
      END IF
C
      CALL QEXIT('QRSPDOT')
C     
      RETURN
      END
      SUBROUTINE GRDCHK(ZERO,INDPRP,NFR,FREQ,
     &     NZCQ,GPC,NZEQ,GPE,NZPQ,GPP,WORK,LWRK)
C*****************************************************************************
C
C     Check norm of gradient vector. If norm less than threshold then
C     assume solution to linear response equation to be the zero vector.
C     Write zero vector to file and return ZERO=.TRUE.
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxrs.h>
#include <dcbxqr.h>
#include <dgroup.h>
#include <dcbxpr.h>
C
      PARAMETER ( D0 = 0.0D0, THD = 1.0D-9 )
C
      LOGICAL ZERO
      DIMENSION GPC(NZCQ),GPE(NZEQ),GPP(NZPQ)
      DIMENSION FREQ(NFR),WORK(*)
C
      KFREE = 1
      LFREE = LWRK
C
      ZERO = .FALSE.
C     
      GPNRM = DNRM2(NZCQ,GPC,1)+DNRM2(NZEQ,GPE,1)+DNRM2(NZPQ,GPP,1)
      IF (GPNRM.LT.THD) THEN
C
         ZERO = .TRUE.
C
         WRITE(LUPRI,'(2(/A),A16,2(/A,D7.1),//2A)') 
     &        ' *** WARNING ****',
     &        ' Skipping property     : ',PRPNAM(INDPRP),
     &        ' since norm of gradient: ',GPNRM,
     &        ' is below threshold    : ',THD,
     &        ' The program will continue assuming response vector',
     &        ' equal to zero vector.'
C
         IF (NZXOPE.GT.0) THEN
            OPEN(LUQREE,FILE='QRVEC.EE',FORM='UNFORMATTED',
     &           ACCESS='DIRECT',RECL=8*MZYEE*NZ,STATUS='UNKNOWN')
            CALL MEMGET('REAL',KDUM,MZYEE*NZ,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KDUM),MZYEE*NZ)
            DO IFR=1,NFR
               CALL INDQR(INDPRP,FREQ(IFR),IREC)
               CALL WRTDAC(LUQREE,MZYEE*NZ,WORK(KDUM),IREC)
               CALL WRTQR(PRPNAM(INDPRP),'E',FREQ(IFR),IREC,D0)
            END DO
            CALL MEMREL('GRDCHK',WORK,1,KDUM,KFREE,LFREE)
            CLOSE(LUQREE,STATUS='KEEP')
         END IF
C
         IF (NZXOPP.GT.0) THEN
            OPEN(LUQREP,FILE='QRVEC.EP',FORM='UNFORMATTED',
     &           ACCESS='DIRECT',RECL=8*MZYEP*NZ,STATUS='UNKNOWN')
            CALL MEMGET('REAL',KDUM,MZYEP*NZ,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KDUM),MZYEP*NZ)
            DO IFR=1,NFR
               CALL INDQR(INDPRP,FREQ(IFR),IREC)
               CALL WRTDAC(LUQREP,MZYEP*NZ,WORK(KDUM),IREC)
               CALL WRTQR(PRPNAM(INDPRP),'E',FREQ(IFR),IREC,D0)
            END DO
            CALL MEMREL('GRDCHK',WORK,1,KDUM,KFREE,LFREE)
            CLOSE(LUQREP,STATUS='KEEP')
         END IF
C
      END IF
C     
      RETURN
      END
      SUBROUTINE READQR(EXST,LABEL,NFREQ,FREQ)
C*****************************************************************************
C
C     Check if the response vectors already exist on files
C     QRVEC.EE and QRVEC.EP. QRINFO is an ASCII file containing 
C     information about the response vectors that have been solved
C     for. See also WRTQR which writes information to QRINFO.
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxqr.h>
C
      LOGICAL EXST
      CHARACTER*1 TYP
      CHARACTER*16 LABEL,LAB
      DIMENSION FREQ(NFREQ)
C
      CALL QENTER('READQR')
C
C     If the QRINFO file is old then open else write EOF symbol
C
      INQUIRE(FILE='QRINFO',EXIST=EXST) 
      IF (EXST) THEN 
         OPEN (UNIT=LUQRINFO,FILE='QRINFO',STATUS='OLD')
      ELSE
         OPEN (UNIT=LUQRINFO,FILE='QRINFO',STATUS='NEW')
         WRITE(LUQRINFO,'(A16)') ' **End of File**'
         GO TO 300
      END IF
C
C     Check if the response vectors exist on file already
C
      EXST = .TRUE.
      DO IFREQ=1,NFREQ
C
         REWIND (LUQRINFO)
 100     READ(LUQRINFO,'(A16,A1,F20.12,I3,F20.12)',END=200,ERR=200)
     &        LAB,TYP,FRQ,IREC,DNORM
         IF (LAB .EQ. ' **End of File**') THEN
            EXST = .FALSE.
            GO TO 300
         ELSEIF (LAB.NE.LABEL .OR. FRQ.NE.FREQ(IFREQ)) THEN
            GO TO 100
         END IF
C
      END DO
C
C     All response vectors found on file, EXST=TRUE on return
C
      GO TO 300
C
C     Error in restart
C
 200  CONTINUE
      WRITE(LUPRI,'(/A,/2A)') 
     &     ' READQR: Error in restart module of QR',
     &     ' Label read from QRINFO: LAB=',LAB 
      CALL QUIT('READQR: Error in restart module of QR')
C
 300  CONTINUE
      CLOSE(LUQRINFO,STATUS='KEEP')
C
      CALL QEXIT('READQR')
      RETURN
      END
      SUBROUTINE WRTQR(LABEL,TYP,FREQ,IREC,DNORM)
C*****************************************************************************
C
C     Write info about response vectors to file QRINFO used for restart.
C     See also READQR which reads through the QRINFO file.
C
C     Written by panor 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbxqr.h>
C
      LOGICAL EXST
      CHARACTER*1 TYP,TYPDUM
      CHARACTER*16 LABEL,LABDUM
C
      CALL QENTER('WRTQR')
C
C     Open file QRINFO
C
      INQUIRE(FILE='QRINFO',EXIST=EXST) 
      IF (EXST) THEN 
         OPEN (UNIT=LUQRINFO,FILE='QRINFO',STATUS='OLD')
      ELSE
         LABDUM='**QRINFO**'
         GO TO 200
      END IF
C
      REWIND (LUQRINFO)
 100  READ(LUQRINFO,'(A16,A1,F20.12,I3,F20.12)',END=200,ERR=200)
     &        LABDUM,TYPDUM,FRQDUM,IDUM,DUM
      IF (LABDUM.EQ.' **End of File**') THEN
         BACKSPACE(LUQRINFO)
         WRITE(LUQRINFO,'(A16,A1,F20.12,I3,F20.12)')
     &        LABEL,TYP,FREQ,IREC,DNORM
         WRITE(LUQRINFO,'(A16)') ' **End of File**'
         GO TO 300
      ELSE
         GO TO 100
      END IF
C
C     Error in restart
C
 200  CONTINUE
      WRITE(LUPRI,'(/A,/2A)') 
     &     ' WRTQR: Error in restart module of QR',
     &     ' Label read from QRINFO: LABDUM=',LABDUM 
      CALL QUIT('WRTQR: Error in restart module of QR')
C
 300  CONTINUE
      CLOSE(LUQRINFO,STATUS='KEEP')
C
      CALL QEXIT('WRTQR')
      RETURN
      END
