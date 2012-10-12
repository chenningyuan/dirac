!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

#ifdef MOD_DNF
#ifndef MOD_LINEARS
* Routines for inversion
* From Jeppe - late november 03
      SUBROUTINE INVERT_BY_DIAG2(A,B,SCR,VEC,NDIM)
*
* Invert a symmetric  - hopefully nonsingular - matrix A 
* by diagonalization. 
*
* Jeppe Olsen, Oct 97 to check INVMAT
*              March 00 : Scale initial matrix to obtain unit diagonal
*
*. The input matrix is delivered in A in complete form, and 
*. the output matrix is returned in A, also in complete form
      IMPLICIT REAL*8(A-H,O-Z)
*. Input and output matrix
      DIMENSION A(NDIM*NDIM)       
*. Scratch matrices and vector
      DIMENSION B(NDIM*NDIM),SCR(NDIM*NDIM),VEC(NDIM)
*
      NTEST = 00
*. Reform a to symmetric packed form
      CALL TRIPAK(A,SCR,1,NDIM,NDIM)
*. Extract diagonal
      CALL COPDIA(SCR,VEC,NDIM,1)
*
*.scale 
*
      DO I = 1, NDIM
*. Scaling vector
        IF(VEC(I).EQ.0.0D0) THEN
          VEC(I) = 1.0D0
        ELSE
          VEC(I) = 1.0D0/SQRT(ABS(VEC(I)))
        END IF
      END DO
*. Scale matrix
      IJ = 0
      DO I = 1, NDIM
        DO J = 1, I
          IJ = IJ + 1
          SCR(IJ) = SCR(IJ)*VEC(I)*VEC(J)
        END DO
      END DO
C     DO I = 1, NDIM
C       VEC(I) = 1.0D0/VEC(I)
C     END DO
*. Diagonalize
      CALL EIGEN(SCR,B,NDIM,0,1)
*. Scale eigenvectors
      DO IVEC = 1, NDIM 
        IOFF = 1 + (IVEC-1)*NDIM
        CALL VVTOV(B(IOFF),VEC,B(IOFF),NDIM)
      END DO
*. 
      CALL COPDIA(SCR,VEC,NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of scaled matrix : '
        CALL WRTMAT(VEC,NDIM,1,NDIM,1)
      END IF
*. Invert diagonal elements 
      DO I = 1, NDIM
       IF(ABS(VEC(I)).GT.1.0D-15) THEN
         VEC(I) = 1.0D0/VEC(I)
       ELSE
         VEC(I) = 0.0D0
         WRITE(6,*) ' Singular mode activated '
       END IF
      END DO
*. and obtain inverse matrix by transformation
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(A,B,VEC,NDIM,SCR)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverse matrix from INVERSE_BY_DIAG'
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END 
#endif /* ! MOD_LINEARS */
#endif /* MOD_DNF */


