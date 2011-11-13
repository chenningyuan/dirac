C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
C
C
C VERSION : $Revision: 1.19.4.1 $
C DATE    : $Date: 2004/09/01 08:17:56 $
C FILE    : ccxblas.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XAXPY (N,A,X,INCX,Y,INCY)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SAXPY, DAXPY, CAXPY, ZAXPY
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX,INCY
      REAL*8 A,X(*),Y(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CALL ZAXPY (N,A,X,INCX,Y,INCY)
C64B     CALL CAXPY (N,A,X,INCX,Y,INCY)
      ELSE
         CALL DAXPY (N,A,X,INCX,Y,INCY)
C64B     CALL SAXPY (N,A,X,INCX,Y,INCY)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOPY (N,X,INCX,Y,INCY)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SCOPY, DCOPY, CCOPY, ZCOPY
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX,INCY
      REAL*8 X(*),Y(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CALL ZCOPY (N,X,INCX,Y,INCY)
C64B     CALL CCOPY (N,X,INCX,Y,INCY)
      ELSE
         CALL DCOPY (N,X,INCX,Y,INCY)
C64B     CALL SCOPY (N,X,INCX,Y,INCY)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     &                  BETA,C,LDC)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SGEMM, DGEMM, CGEMM, ZGEMM
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*1 TRANSA,TRANSB
      INTEGER M,N,K,LDA,LDB,LDC
      REAL*8 ALPHA,BETA,A(*),B(*),C(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*1 TRANSAL,TRANSBL
C
C---------------Executable code--------------------------------------
C
      IF (M.EQ.0.OR.N.EQ.0.OR.K.EQ.0) RETURN
C
      IF (CARITH) THEN
         CALL ZGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     &               BETA,C,LDC)
C64B     CALL CGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
      ELSE
         TRANSAL = TRANSA
         TRANSBL = TRANSB
         IF (TRANSA.EQ.'C') TRANSAL = 'T'
         IF (TRANSB.EQ.'C') TRANSBL = 'T'
         IF (TRANSA.EQ.'c') TRANSAL = 't'
         IF (TRANSB.EQ.'c') TRANSBL = 't'
         CALL DGEMM (TRANSAL,TRANSBL,M,N,K,ALPHA,A,LDA,B,LDB,
     &               BETA,C,LDC)
C64B     CALL SGEMM (TRANSAL,TRANSBL,M,N,K,ALPHA,A,LDA,B,LDB,
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XGEMV (TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SGEMV, DGEMV, CGEMV, ZGEMV
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*1 TRANS
      INTEGER M,N,LDA,INCX,INCY
      REAL*8 ALPHA,BETA,A(*),X(*),Y(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*1 TRANSL
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CALL ZGEMV (TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
C64B     CALL CGEMV (TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      ELSE
         TRANSL = TRANS
         IF (TRANS.EQ.'C') TRANSL = 'T'
         IF (TRANS.EQ.'c') TRANSL = 't'
         CALL DGEMV (TRANSL,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
C64B     CALL SGEMV (TRANSL,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XSCAL (N,A,X,INCX)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SSCAL, DSCAL, CSCAL, ZSCAL
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX
      REAL*8 A,X(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CALL ZSCAL (N,A,X,INCX)
C64B     CALL CSCAL (N,A,X,INCX)
      ELSE
         CALL DSCAL (N,A,X,INCX)
C64B     CALL SSCAL (N,A,X,INCX)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XSUM (RESULT,N,X,INCX)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines SSUM, DSUM, CSUM, ZSUM
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX
      REAL*8 X(*)
      REAL*8 RESULT(*),DSUM
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         RESULT(1) = DSUM (N,X(1),2*INCX) 
         RESULT(2) = DSUM (N,X(2),2*INCX) 
      ELSE
         RESULT(1) = DSUM (N,X,INCX) 
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined (SYS_AIX)
      REAL*8 FUNCTION DSUM (N,X,INCX)
      REAL*8 X(*),SUM
      II = 1
      SUM = 0.D0
      DO I = 1, N
         SUM = SUM + X(II)
         II = II + INCX
      ENDDO
      DSUM = SUM
      RETURN
      END
#endif
C
#if defined (SYS_AIX)  || defined (SYS_IRIX) || defined (SYS_LINUX) \
 || defined (SYS_HPUX) || defined (SYS_T3D)  || defined (SYS_SUN) || defined (SYS_DARWIN)
C
      COMPLEX*16 FUNCTION ZSUM (N,X,INCX)
      COMPLEX*16 X(*),SUM
      SUM = DCMPLX(0.D0,0.D0)
      II = 1
      DO I = 1, N
         SUM = SUM + X(II)
         II = II + INCX
      ENDDO
      ZSUM = SUM
      RETURN
      END
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XTRANSP (M,N,AR,BR,AC,BC)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Generic transpose routine.
C     Copies matrix A(M,N) into matrix B(N,M)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER M,N
      REAL*8 AR(M,N),BR(N,M)
      COMPLEX*16 AC(M,N),BC(N,M)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER I,J
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         DO I = 1, M
            DO J = 1, N
               BC(J,I) = AC(I,J)
            ENDDO
         ENDDO
      ELSE
         DO I = 1, M
            DO J = 1, N
               BR(J,I) = AR(I,J)
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XDOTC (RESULT,N,X,INCX,Y,INCY)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Generic version of BLAS routines DOTC and DOTU
C     Depending of the value of CARITH and the machine specification
C     the appropriate routine is called
C
C     Note : The routine was written due to problems with the BLAS 
C            ZDOTC and ZDOTU functions, be careful if reintroducing them
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C     Date : 25-9-2003
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX,INCY
      REAL*8 X(*),Y(*)
      REAL*8 DDOT, RESULT(*)
C
C---------------Common Blocks--------------------------------------
C
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         RESULT(1) = DDOT (N,X(1),2*INCX,Y(1),2*INCY) 
     &             + DDOT (N,X(2),2*INCX,Y(2),2*INCY)
         RESULT(2) = DDOT (N,X(1),2*INCX,Y(2),2*INCY) 
     &             - DDOT (N,X(2),2*INCX,Y(1),2*INCY)
      ELSE
         RESULT(1) = DDOT (N,X,INCX,Y,INCY)
      ENDIF
C
      RETURN
C
      ENTRY XDOTU (RESULT,N,X,INCX,Y,INCY)
C
      IF (CARITH) THEN
         RESULT(1) = DDOT (N,X(1),2*INCX,Y(1),2*INCY) 
     &             - DDOT (N,X(2),2*INCX,Y(2),2*INCY)
         RESULT(2) = DDOT (N,X(1),2*INCX,Y(2),2*INCY) 
     &             + DDOT (N,X(2),2*INCX,Y(1),2*INCY)
      ELSE
         RESULT(1) = DDOT (N,X,INCX,Y,INCY)
      ENDIF
C
      RETURN
C
      END
