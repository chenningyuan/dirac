!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

C VERSION : $Revision: 11757 $
C DATE    : $Date: 2010-10-10 09:06:37 +0200 (Sun, 10 Oct 2010) $
C FILE    : vcraypack.F
C
C
#ifdef UNDEF
/**COMDECK log
Define labels used: ESSL for CALL DGEMUL in MXMA
941223-hjaaj: merged with aauc32:~qcprg/gp/vcraypac.f
941208-hjaaj:
  From Lund ; Alliant arhfx8 file 'vcraypack.ibm' was dated Jan 05 1988
  The lund tarfile is now deleted, this was the only of possible interest.
*DECK VCRAYPACK
*/
#endif
      SUBROUTINE VSADD(A,LA,S,C,LC,N)
#include "implicit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION A(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=A(1+LA*I)+S
10    CONTINUE
      RETURN
      END
      SUBROUTINE VSMA(A,LA,S,B,LB,C,LC,N)
#include "implicit.h"
      DIMENSION A(*),B(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=A(1+LA*I)*S+B(1+LB*I)
10    CONTINUE
      RETURN
      END
      SUBROUTINE VSMUL(A,LA,S,C,LC,N)
#include "implicit.h"
      DIMENSION A(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=A(1+LA*I)*S
10    CONTINUE
      RETURN
      END
      SUBROUTINE VSQRT(A,LA,C,LC,N)
#include "implicit.h"
      DIMENSION A(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=SQRT(A(1+LA*I))
10    CONTINUE
      RETURN
      END
      SUBROUTINE VSUB(A,LA,B,LB,C,LC,N)
#include "implicit.h"
      DIMENSION A(*),B(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=B(1+LB*I)-A(1+LA*I)
10    CONTINUE
      RETURN
      END
      SUBROUTINE VMOV(A,LA,B,LB,N)
#include "implicit.h"
      DIMENSION A(*),B(*)
      DO 10 I=0,N-1
         B(1+I*LB)=A(1+I*LA)
10    CONTINUE
      RETURN
      END
      SUBROUTINE VMSA(A,LA,B,LB,S,C,LC,N)
#include "implicit.h"
      DIMENSION A(*),B(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=A(1+LA*I)*B(1+LB*I)+S
10    CONTINUE
      RETURN
      END
      SUBROUTINE VMUL(A,LA,B,LB,C,LC,N)
#include "implicit.h"
      DIMENSION A(*),B(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=A(1+LA*I)*B(1+LB*I)
10    CONTINUE
      RETURN
      END
      SUBROUTINE VDIV(A,LA,B,LB,C,LC,N)
#include "implicit.h"
      DIMENSION A(*),B(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=B(1+LB*I)/A(1+LA*I)
10    CONTINUE
      RETURN
      END
      SUBROUTINE VADD(A,LA,B,LB,C,LC,N)
#include "implicit.h"
      DIMENSION A(*),B(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=A(1+LA*I)+B(1+LB*I)
10    CONTINUE
      RETURN
      END
#if defined (SYS_CRAY)
      SUBROUTINE VCLR(A,LA,N)
#include "implicit.h"
      DIMENSION A(*)
      PARAMETER ( D0 = 0.0D0 )
      DO 10 I=0,N-1
         A(1+I*LA)=D0
10    CONTINUE
      RETURN
      END
#endif
      SUBROUTINE MXMA(A,ICA,IRA,B,ICB,IRB,C,ICC,IRC,NROW,NSUM,NCOL)
#include "implicit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION A(*),B(*),C(*)
#if !defined (VAR_ESSL) 
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,NCOL-1
         SUM=D0
         DO 110 ISUM=0,NSUM-1
            SUM=SUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         INDEX=1+ICC*IROW+IRC*ICOL
         C(INDEX)=SUM
100   CONTINUE
#else
      IF(NSUM.EQ.0) GOTO 600
      IF((NROW.EQ.0).OR.(NCOL.EQ.0)) RETURN
C  TRY TO FIT TO A STANDARD CASE AND SET SELECTOR:
      ISEL=1
      IF ((ICA.EQ.1).AND.(IRA.GE.NROW)) GOTO 101
      IF ((IRA.EQ.1).AND.(ICA.GE.NSUM)) GOTO 102
      GOTO 400
 101  ISEL=ISEL+4
 102  IF ((ICB.EQ.1).AND.(IRB.GE.NSUM)) GOTO 201
      IF ((IRB.EQ.1).AND.(ICB.GE.NCOL)) GOTO 202
      GOTO 400
 201  ISEL=ISEL+2
 202  IF ((ICC.EQ.1).AND.(IRC.GE.NROW)) GOTO 301
      IF ((IRC.EQ.1).AND.(ICC.GE.NCOL)) GOTO 302
      GOTO 400
 301  ISEL=ISEL+1
C  NOTE -- CASES 2 OR 7 IN COMPUTED GOTO STATEMENT BELOW ARE OF
C  THE T-T TYPE, WHICH SOMETIMES GOES WRONG IN DGEMUL.
C  FOR THE TIME BEING, THIS MUST BE TREATED AS A GENERAL CASE.
C  CHANGE BACK AS SOON AS ESSL HAS BEEN PATCHED!! P-AA M 87-03-23.
C302  GOTO (1,2,3,4,5,6,7,8) ISEL
 302  GOTO (1,400,3,4,5,6,400,8) ISEL
      WRITE(*,*)' SOMETHING IS VERY WRONG -- THIS SPOT SHOULD NOT BE'
      WRITE(*,*)' REACHABLE!! MXMA STOPS.'

   1  CALL DGEMUL(B,ICB,'N',A,ICA,'N',C,ICC,NCOL,NSUM,NROW)
      RETURN
 
C  2  CALL DGEMUL(A,ICA,'T',B,ICB,'T',C,IRC,NROW,NSUM,NCOL)
C     RETURN
 
   3  CALL DGEMUL(B,IRB,'T',A,ICA,'N',C,ICC,NCOL,NSUM,NROW)
      RETURN
 
   4  CALL DGEMUL(A,ICA,'T',B,IRB,'N',C,IRC,NROW,NSUM,NCOL)
      RETURN
 
   5  CALL DGEMUL(B,ICB,'N',A,IRA,'T',C,ICC,NCOL,NSUM,NROW)
      RETURN
 
   6  CALL DGEMUL(A,IRA,'N',B,ICB,'T',C,IRC,NROW,NSUM,NCOL)
      RETURN
 
C  7  CALL DGEMUL(B,IRB,'T',A,IRA,'T',C,ICC,NCOL,NSUM,NROW)
C     RETURN
 
   8  CALL DGEMUL(A,IRA,'N',B,IRB,'N',C,IRC,NROW,NSUM,NCOL)
      RETURN
 
 400  DO 520 IROW=0,NROW-1
        DO 520 ICOL=0,NCOL-1
          T=D0
          DO 510 ISUM=0,NSUM-1
            T=T+A(1+IROW*ICA+ISUM*IRA)*B(1+ISUM*ICB+ICOL*IRB)
 510      CONTINUE
          C(1+IROW*ICC+ICOL*IRC)=T
 520  CONTINUE
      RETURN
 
 600  DO 620 IROW=0,NROW-1
        DO 620 ICOL=0,NCOL-1
          C(1+IROW*ICC+ICOL*IRC)=0.0D00
 620  CONTINUE
#endif
      RETURN
      END
      SUBROUTINE MXMAA(A,ICA,IRA, B,ICB,IRB, C,ICC,IRC, NROW,NSUM,NCOL)
#include "implicit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION A(*),B(*),C(*)
 
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,NCOL-1
         SUM=D0
         DO 110 ISUM=0,NSUM-1
            SUM=SUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         INDEX=1+ICC*IROW+IRC*ICOL
         C(INDEX)=C(INDEX)+SUM
100   CONTINUE
      RETURN
      END
      SUBROUTINE MXMAM(A,ICA,IRA, B,ICB,IRB, C,ICC,IRC, NROW,NSUM,NCOL)
#include "implicit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION A(*),B(*),C(*)
 
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,NCOL-1
         SUM=D0
         DO 110 ISUM=0,NSUM-1
            SUM=SUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         INDEX=1+ICC*IROW+IRC*ICOL
         C(INDEX)=C(INDEX)-SUM
100   CONTINUE
      RETURN
      END
      SUBROUTINE MXMT(A,ICA,IRA, B,ICB,IRB, C, NROW,NSUM)
#include "implicit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION A(*),B(*),C(*)
      IND=0
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,IROW
         SUM=D0
         DO 110 ISUM=0,NSUM-1
            SUM=SUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         IND=IND+1
         C(IND)=SUM
100   CONTINUE
      RETURN
      END
