***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
*   Routine for (complex) transposing rectangular matrix A
*
      subroutine adjoint(B,A,NROWS,NCOLUMNS,IMRECO)
      implicit real*8 (A-H,O-Z)
      dimension A(NROWS,NCOLUMNS),B(NROWS,NCOLUMNS)
*
*   Write zeroes to matrix B
*
      do K=1,NCOLUMNS,1
         do L=1,NROWS,1
            B(L,K) = 0.
         end do
      end do
*
*  transpose
*
      do I=1,NCOLUMNS,1
         do J=1,NROWS,1
            B(I,J) = IMRECO * A(J,I)
         end do
      end do
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine clsone
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
*   Routine(1) for error message
*
      subroutine hammer(LAST,IMAX)
      implicit none
      integer LAST,IMAX
      write(*,*) 'The dimensions have exceeded the maximum work space!'
      write(*,*) 'Dimensions are  ',LAST
      write(*,*) 'Maximum space for this part is ',IMAX
      write(*,*) 'Picasso will put down the brush.'
      write(*,*) 'Bye-bye.'
      Call Abend()
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION I_GET_ADJSYM(ISYM)
*
* Obtain conjugate symmetry of ISYM
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "symm.inc"
#include "lucinp.inc"
*
      IF(PNTGRP.EQ.1) THEN
*. Usual pointgroup, symmetry and complex symmetry identical
        I_GET_ADJSYM = ISYM
      ELSE
*. Currently only working alternative is double group symmetry
        I_GET_ADJSYM = IADJSYM(ISYM)
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' I_GET_ADJSYM , ISYM, IADJSYM',
     &  ISYM, IADJSYM(ISYM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION IFRMR_REL(WORK,IROFF,IELMNT)
*
* An integer array is stored in real array WORK,
* starting from WORK(IROFF). Obtain element
* IELMNT of this array
*
      INTEGER WORK(*)
*
#include "irat.inc"
C     PARAMETER(IRAT = 2 )
*. offset when work is integer array
      IIOFF = 1 + IRAT * (IROFF-1)
      IFRMR_REL = WORK(IIOFF-1+IELMNT)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION IMNMX_REL(IVEC,NDIM,MINMAX)
*
*     Find smallest (MINMAX=1) or largest (MINMAX=2)
*     absolute value of elements in integer vector IVEC
*
      DIMENSION IVEC(NDIM)
*
      IX = 0
      IF(MINMAX.EQ.1) THEN
        IF(NDIM.GT.0) IX=ABS(IVEC(1))
        DO 100 I=2,NDIM
          IX=MIN(IX,ABS(IVEC(I)))
  100   CONTINUE
      END IF
*
      IF(MINMAX.EQ.2) THEN
        IF(NDIM.GT.0) IX=ABS(IVEC(1))
        DO 200 I=2,NDIM
          IX=MAX(IX,ABS(IVEC(I)))
  200   CONTINUE
      END IF
*
      IMNMX_REL = IX
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Use Kramers symmetry (if flag is activated) on a particle
* of an integral.
* Generalized to all kinds of index exchange, including complex
* conjugation symmetry.
*
      subroutine index_exch(I,J,IFAC,INDFAC,IFLOC)
*
*  I,J are indices to be interchanged (by some symmetry)
*  IFAC: Total factor multiplied to integral
*  INDFAC: Multiply indices by this factor (usually +-1)
*  IFLOC: Input symmetry factor for this exchange
*
      include 'implicit.inc'
*
      NTESTL = 0
      if (NTESTL.ge.1) then
        write(6,*)
        write(6,*) 'INDEX_EXCH : input indices ',I,J
        write(6,*) '             input factor  ',IFAC
      end if
*
      IH = I
      JH = J
      I = JH * INDFAC
      J = IH * INDFAC
      IFAC = IFAC * IFLOC
*
      if (NTESTL.ge.1) then
        write(6,*) 'INDEX_EXCH : output indices ',I,J
        write(6,*) '             output factor  ',IFAC
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE INPRDDC(VEC1,VEC2,LU1,LU2,
     &                   I1,I2,IREW,LBLK,XR,XI)
*
* Two complex vectors are given : vec1 = vector number I1 on LU1
*                                 Vec2 = vector number I2 on LU2
*
* Obtain inner product of these two vectors  <Vec(LU1)!Vec(Lu2)>
*
* output given as XR : real part of inner product
*                 XI : Imaginary part of inner product
* If LU1 = LU2, the imaginary part is set to zero
*
* If LU1.NE.LU2, both the real and imaginary parts are
* calculated. This requires rewinding of LU1 and LU2
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL*8 INPRDD
*
      NTEST = 0
*
      IF(IREW.NE.0) THEN
*. Position files at start of complex vectors I1, I2
        CALL SKPVCD(LU1,2*(I1-1),VEC1,1,LBLK)
C              SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
        IF(LU2.NE.LU1) THEN
          CALL SKPVCD(LU2,2*(I2-1),VEC1,1,LBLK)
        END IF
      END IF
*
      IF(LU1.EQ.LU2) THEN
*
* Norm of vector
*
        XRR = INPRDD(VEC1,VEC2,LU1,LU2,0,LBLK)
        IF(NTEST.NE.0)
     &  WRITE(6,*) ' Real part * Real part ', XRR
        XII =  INPRDD(VEC1,VEC2,LU1,LU2,0,LBLK)
        IF(NTEST.NE.0)
     &  WRITE(6,*) ' Imaginary part * Imaginary part ', XII
        XIR = 0.0D0
        XRI = 0.0D0
      ELSE
*
*. General inner product of two vectors,
*
*. Real * Real
        XRR = INPRDD(VEC1,VEC2,LU1,LU2,0,LBLK)
        IF(NTEST.NE.0)
     &  WRITE(6,*) ' Real part * Real part ', XRR
*. Real * Imag
        CALL SKPVCD(LU1,2*(I1-1),VEC1,1,LBLK)
        XRI = INPRDD(VEC1,VEC2,LU1,LU2,0,LBLK)
        IF(NTEST.NE.0)
     &  WRITE(6,*) ' Real part * Imag part ', XRI
*. Imag*real
        CALL SKPVCD(LU2,2*(I2-1),VEC1,1,LBLK)
        XIR = INPRDD(VEC1,VEC2,LU1,LU2,0,LBLK)
        IF(NTEST.NE.0)
     &  WRITE(6,*) ' Imag part * Real part ', XIR
*. Imag*imag
        CALL SKPVCD(LU1,2*(I1-1)+1,VEC1,1,LBLK)
        CALL SKPVCD(LU2,2*(I2-1)+1,VEC1,1,LBLK)
        XII = INPRDD(VEC1,VEC2,LU1,LU2,0,LBLK)
        IF(NTEST.NE.0)
     &  WRITE(6,*) ' Imag part * Imag part ', XII
      END IF
*
      XR = XRR + XII
      XI = XRI - XIR
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Unfold lower triangle square matrix to full size;
*  I suppose that the input matrix is a one-dim. field with
*  N(N/2 + 1) elements
*  1) write vector to matrix
*  2) unfold considering symmetry (e.g. hermitian or not)
*
      subroutine unfold(ARRAY,VECTOR,IDIM,ISYM,NTRIAEL,IFAC)
      implicit real*8 (A-H,O-Z)
      dimension ARRAY(IDIM,IDIM),VECTOR(NTRIAEL)
*
*  unfold
*
      ICOUNT = 0
      do N = 1,IDIM,1
         do M = 1,N,1
            ARRAY(N,M) = IFAC * VECTOR(ICOUNT+M)
         end do
         ICOUNT = ICOUNT + N
      end do
*
*  multiply with sign
*
      do I=1,IDIM,1
         do J=I+1,IDIM,1
            ARRAY(I,J) = ISYM * ARRAY(J,I)
         end do
      end do
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Copy vector to defined partition of matrix
* (Simply the way the vector is arranged in memory)
*
      subroutine vectomat(VEC,WMAT,ILEN,NROW,NCOL,
     &                    NPARTROW,NPARTCOL,IROWOFF,ICOLOFF)
*
      implicit real*8 (A-H,O-Z)
*
      dimension WMAT(NROW,NCOL)
      dimension VEC(ILEN)
*
*   VEC: Input matrix in vector format
*   WMAT: Output matrix (can be a partition of a bigger matrix!)
*   ILEN: Length of input vector
*   NROW: Rows of WMAT
*   NCOL: Columns of WMAT
*   NPARTROW: Number of rows in output partition
*   NPARTCOL: Number of columns in output partition
*   IROWOFF: Output row offset
*   ICOLOFF: Output column offset
*
      K = 1
      do I=ICOLOFF,(ICOLOFF+NPARTCOL-1),1
         do J=IROWOFF,(IROWOFF+NPARTROW-1),1
            WMAT(J,I) = VEC(K)
            K = K + 1
         end do
      end do
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Read indices from file INDICES
*
      subroutine read_indices(IN,NINTS)
#include "clunit.inc"
      dimension IN(4,NINTS)
*
      open(unit=LUSC95,file='INDICES',
     &       status='OLD',form='FORMATTED')
      rewind(LUSC95)
      do I=1,NINTS,1
        read(LUSC95,'(4I4)') (IN(J,I),J=1,4)
      end do
      close(LUSC95)
*
      return
      end
