      SUBROUTINE KRMCVEC(IBTYP,IBCVC,IBEVC,IBPVC,EVALR,EVECR,
     &     DKAPPA,CREF,WORK,LWORK)
C*****************************************************************************
C     
C     Form the response vectors and write them to file given the 
C     solutions from the reduced system of trial vectors. 
C     
C     Originally written by panor 1998
C     Modified by jth Nov 24 1998: write to LUKRMC instead.
C     
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C     
      PARAMETER (DP5 = 0.500, D2 = 2.0D00, D4 = 4.0D00)
C     
#include <dummy.h>
#include <dcbxrs.h>
#include <dcbibn.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbgen.h>
C     
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),
     &     IBPVC(*),EVALR(*),EVECR(NREDM,*),DKAPPA(*)
      DIMENSION CREF(*)
      DIMENSION WORK(*)
C     
#include <memint.h>
C     
      CALL QENTER('KRMCVEC')
C     
C     Solution vector - orbital (e-e) part
C     
      IF(NZXOPE.GT.0) THEN
         write(6,*) 'e-e'
         CALL MEMGET('REAL',KBBUF,NZXOPEQ,WORK,KFREE,LFREE)
         OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
         CALL XRSXV1(JTIMOP,JBENDX,DKAPPA,EVECR,NFREQ,
     &        IBTYP,IBEVC,IVEC,WORK(KBBUF))
         CLOSE(LUBOE,STATUS='KEEP')
         CALL MEMREL('KRMCVEC.e-e',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF
C     
C     Solution vector - orbital (e-p) part
C     
      IF(NZXOPP.GT.0) THEN
         write(6,*) 'e-p'
         CALL MEMGET('REAL',KBBUF,NZXOPPQ,WORK,KFREE,LFREE)
         OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
         CALL XRSXV1(JTIMOP,JBPNDX,DKAPPA(1+NZXOPEQ),EVECR,NFREQ,
     &        IBTYP,IBPVC,IVEC,WORK(KBBUF))
         CLOSE(LUBOP,STATUS='KEEP')
         CALL MEMREL('KRMCVEC.e-p',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF
C     
C     Solution vector - conf. part
C     
      IF ( NZCONF .GT. 0 ) THEN
         write(6,*) 'conf'
         CALL MEMGET('REAL',KBBUF,NZCONFQ,WORK,KFREE,LFREE)
         OPEN(LUBCI,FILE='PAMBCI',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*NZCONFQ,STATUS='UNKNOWN')
         CALL XRSXV1(JTIMOP,JBCNDX,CREF,EVECR,NFREQ,
     &        IBTYP,IBCVC,IVEC,WORK(KBBUF))
         CLOSE(LUBCI,STATUS='KEEP')
         CALL MEMREL('KRMCVEC.conf',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF
C     
C     Write to KAPPA to LUKRM1:
C     
      REWIND(LUKRM1)
C     
C     Scale with .5 (Y+Z=2Y is saved from xrs)
C     ...and scale with 4 from new RSP code.
C     
      write(6,*) 'krmcvec: writing kappa:',nzxoptq
c     CALL DSCAL(NZXOPTQ,D4,DKAPPA,1)
      CALL WRITT(LUKRM1,NZXOPTQ,DKAPPA)
C     
C     
      CALL QEXIT('KRMCVEC')
      RETURN
      END
      SUBROUTINE krmcvec1(NBTYP,IBTYP,IBVEC,
     &     EVALR,EVECR,XPO,BBUF,NZXPAR)
C*****************************************************************************
C
C     Form the response vectors and write them to file given the 
C     solutions from the reduced system of trial vectors. 
C
C     LRF    = linear response funtion values
C     NZXPAR = length of response vector
C     NBTYP  = 1 for e-e and 2 e-p part, respectively (input)
C     XPO    = space for response vectors (output)
C     EVALR  = frequencies in the reduced system (E[2] - w*S[2])*N = GP
C     EVECR  = space for eigenvectors of reduced system (read from file)
C     NFREQ  = number of frequencies (input)
C     IBTYP  = type of B trial vectors, dim=2*NREDM (input)
C     IBVEC  = pointers to electron/positron like trial vectors (input)
C     BBUF   = extra space of response vector length
C     LUBVEC = unit number of direct access file for the trial vectors 
C     INDPRP = pointer to the property for which response is evaluated
C     RSPTYP = string indicating what order response function the 
C              linear response vectors are used for, either equal to
C              'QRVEC' or 'CRVEC' for quadratic or cubic response (to be 
C              implemnted).
C
C     Written by panor 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <dcbxrs.h>
#include <dcbibn.h>
#include <dcborb.h>
#include <dgroup.h>
C
#include <dcbxpr.h>
C
      CHARACTER FILNAM*6,FILRSP*8,TYP*1,RSPTYP*5
      DIMENSION IBTYP(2,*),IBVEC(*),EVALR(*),EVECR(NREDM,*),
     &          XPO(NZXPAR,NZ,NFREQ),BBUF(*)
C
      CALL QENTER('KRMCVEC1')
C
      KFREE = 1
      LFREE = LWRK
C
      IF (NBTYP.EQ.JBENDX) THEN
        NBRED  = NERED
        MZXPAR = NZXOPE
        LUBVEC = LUBOE
        FILNAM = 'PAMBOE'
        TYP    = 'E'
      ELSEIF (NBTYP.EQ.JBPNDX) THEN
        NBRED  = NPRED
        MZXPAR = NZXOPP
        LUBVEC = LUBOP
        FILNAM = 'PAMBOP'
        TYP    = 'P'
      ELSEIF (NBTYP.EQ.JBCNDX) THEN
        NBRED  = NCRED
        MZXPAR = NZCONF
        LUBVEC = LUSCI
        FILNAM = 'PAMBCI'
        TYP    = 'C'
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'KRMCVEC1 ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('KRMCVEC1: Unknown NBTYP!')
      ENDIF
C
C     Construct response solution vectors by summing the trial
C     vectors multiplied by the eigenvectors in the reduced
C     system.
C
      OPEN(LUBVEC,FILE=FILNAM,FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*NZXPAR*NZ,STATUS='UNKNOWN')
      CALL DZERO(XPO,NZXPAR*NZ*NFREQ)
      CALL XRSXV1(JTIMOP,NBTYP,XPO,EVECR,NFREQ,IBTYP,IBVEC,IVEC,BBUF)
      CLOSE(LUBVEC,STATUS='KEEP')
C
C
C     Print section
C
      IF(IPRXLR.GE.0) THEN
         CALL HEADER('KRMCVEC1: Solution vectors in KRMCVEC1',-1)
         CALL PRBVEC(LUPRI,XPO,NFREQ,NZXPAR)
         CALL FLSHFO(LUPRI)
      ENDIF
C
      CALL QEXIT('KRMCVEC1')
      RETURN
      END
      SUBROUTINE PRINT_REAL(N,VEC)
C*****************************************************************************
C
C     Formatted print of vector of real numbers to LUPRI.
C
C     Written by panor 1998
C
C*****************************************************************************
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( NCOL=5 )
      DIMENSION VEC(*)
C
      NROW=INT((N-1)/NCOL)
      DO I=1,NROW
         WRITE(LUPRI,'(5(F12.6))') (VEC(1+I*(NROW-1)+J), J=1,NCOL)
      END DO
C
      WRITE(LUPRI,'(5(F12.6))') (VEC(NROW*NCOL+I),I=1,N-NROW*NCOL)
C
      RETURN
      END

