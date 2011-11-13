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

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
*
*   Routine for (complex) transposing rectangular matrix A
*
      subroutine adjoint(B,A,NROWS,NCOLUMNS,IMRECO)
*
#include "implicit.inc"
*
      dimension A(NROWS,NCOLUMNS),B(NROWS,NCOLUMNS)
*
*   Write zeroes to matrix B
*
      ZERO = 0.D0
      do K=1,NCOLUMNS,1
         do L=1,NROWS,1
            B(L,K) = ZERO
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
*           parallelization by Stefan Knecht                          *
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
*           parallelization by Stefan Knecht                          *
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
*           parallelization by Stefan Knecht                          *
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
*           parallelization by Stefan Knecht                          *
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
*     offset when work is integer array
      IIOFF = 1 + IRAT * (IROFF-1)
      IFRMR_REL = WORK(IIOFF-1+IELMNT)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER FUNCTION ISACTIVE_CPU(NPARBLOCK,NDIM,ISPROC)
*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
*     determine whether CPU is active or not
*
***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION NPARBLOCK(NDIM)
*
      ISACTIVE_DUM = 0
C
      DO I = 1, NDIM
        IF( NPARBLOCK(I) .eq. ISPROC ) ISACTIVE_DUM = 1
      END DO
C
      ISACTIVE_CPU = ISACTIVE_DUM
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
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
#include "implicit.inc"
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
C     IFAC = IFAC * IFLOC
*
      if (NTESTL.ge.1) then
        write(6,*) 'INDEX_EXCH : output indices ',I,J
        write(6,*) '             output factor  ',IFAC
      end if
*
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INPRDD_REAL_CPLX(VEC1,VEC2,LU1,LU2,I1,I2,
     &                            IREW,LBLK,XR,XI,MZ)
C***********************************************************************     
C
C     two in general complex vectors are given : 
C                    vec1 = vector number I1 on LU1
C                    vec2 = vector number I2 on LU2
C
C     obtain inner product of these two vectors  
C     <Vec(LU1)!Vec(Lu2)>
C
C     output given as XR : real part of inner product
C                     XI : imaginary part of inner product
C
C     if LU1 = LU2, the imaginary part is set to zero
C
C     based on INPRDDC.
C
C     written by S. Knecht - Oct 2008
C
C***********************************************************************     
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL*8 INPRDD
      integer LU1, LU2, I1, I2, IREW, LBLK, MZ
      integer internal_rew
C
      NTEST = 0
      XR  = 0.0D0
      XI  = 0.0D0
      XRR = 0.0D0
      XRI = 0.0D0
      XIR = 0.0D0
      XII = 0.0D0
      internal_rew = 0
C
C     position files at start of complex vectors I1, I2
      IF(IREW.NE.0) THEN
        CALL SKPVCD(LU1,MZ*(I1-1),VEC1,IREW,LBLK)
        IF(LU2.NE.LU1) THEN
          CALL SKPVCD(LU2,MZ*(I2-1),VEC1,IREW,LBLK)
        END IF
      END IF
C
      IF(LU1.EQ.LU2) THEN
C
C       Norm of vector
        XRR = INPRDD(VEC1,VEC2,LU1,LU2,internal_rew,LBLK)
        IF( MZ .eq. 2 ) XII =  
     &  INPRDD(VEC1,VEC2,LU1,LU2,internal_rew,LBLK)
        XIR = 0.0D0
        XRI = 0.0D0
      ELSE
C
C       general inner product of two vectors
C       ------------------------------------
C
C       real * real
        XRR = INPRDD(VEC1,VEC2,LU1,LU2,internal_rew,LBLK)
        IF(NTEST.NE.0)
     &  WRITE(6,*) ' Real part * Real part ', XRR
        IF( MZ .eq. 2 )THEN
C
C         real * imag
          CALL SKPVCD(LU1,2*(I1-1),VEC1,1,LBLK)
          XRI = INPRDD(VEC1,VEC2,LU1,LU2,0,LBLK)
          IF(NTEST.NE.0)
     &    WRITE(6,*) ' Real part * Imag part ', XRI
C
C         imag*real
          CALL SKPVCD(LU2,2*(I2-1),VEC1,1,LBLK)
          XIR = INPRDD(VEC1,VEC2,LU1,LU2,0,LBLK)
          IF(NTEST.NE.0)
     &    WRITE(6,*) ' Imag part * Real part ', XIR
C
C         imag*imag
          CALL SKPVCD(LU1,2*(I1-1)+1,VEC1,1,LBLK)
          CALL SKPVCD(LU2,2*(I2-1)+1,VEC1,1,LBLK)
          XII = INPRDD(VEC1,VEC2,LU1,LU2,0,LBLK)
          IF(NTEST.NE.0)
     &    WRITE(6,*) ' Imag part * Imag part ', XII
        END IF
      END IF
C
      XR = XRR + XII
      XI = XRI - XIR
!     print *, 'set XR and XI',XR,XI
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
C**********************************************************************
      INTEGER FUNCTION IPICK_T_ACTIVE(IT_CTALL,IT_CTONE,LPROC_NM,
     &                                NDIM,IMULT)
C
C     find all active T-blocks for all (LPROC_NM) CPU's 
C
C     Written by S. Knecht - Feb 01 2008
C
C     Last revision :
C
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#if defined (VAR_MPI2)
#include "mpif.h"
#include "infpar.h"
#endif
#include "parluci.h"
      DIMENSION IT_CTALL(*), IT_CTONE(*)
      INTEGER ICOUNT_TBLK, IOFF_ONE
C
      IPICK_T_ACTIVE = 0
      ICOUNT_TBLK    = 0
      NTEST          = 0
C
      IF( NTEST .gt. 0 )THEN
        IF( IMULT .eq. 1 ) WRITE(LUWRT,*) ' real part of T'
        IF( IMULT .eq. 2 ) WRITE(LUWRT,*) ' imag part of T'
        WRITE(LUWRT,*) ' T dimension: ',NDIM
        WRITE(LUWRT,*) ' T group dimension: ', LPROC_NM
        IF( IMULT .eq. 1 ) THEN 
          WRITE(LUWRT,*) ' TOTAL T BLOCK ARRAY    '
          CALL IWRTMAMN(IT_CTALL,LPROC_NM,NDIM*IRC_SAVE,LPROC_NM,
     &                  NDIM*IRC_SAVE,LUWRT)
        END IF
      END IF
C
C     start search...
C
      DO 200 ICT = 1, NDIM
C
        INUM_P = 0
 100    CONTINUE
        INUM_P = INUM_P + 1
        IOFF   = NDIM * IRC_SAVE*(INUM_P-1)+(NDIM*(IMULT-1)) + ICT
        IF( IT_CTALL(IOFF) .gt. 0 )THEN
          ICOUNT_TBLK = ICOUNT_TBLK + 1
          IOFF_ONE = ICOUNT_TBLK + NDIM * ( IMULT - 1 )
csk          IT_CTONE( IOFF_ONE ) = ICT + NDIM * ( IMULT - 1 )
          IT_CTONE( IOFF_ONE ) = ICT
C         this block is active, go for the next one
          GOTO 200
        END IF
        IF( INUM_P .lt. LPROC_NM ) GOTO 100
C
 200  CONTINUE
C
csk      WRITE(LUWRT,'(A,I7)') ' total active T blocks: ',ICOUNT_TBLK
C
csk      WRITE(LUWRT,*) ' ACTIVE T BLOCKS'
csk      CALL IWRTMAMN(IT_CTONE,1,NDIM*IRC_SAVE,1,NDIM*IRC_SAVE,LUWRT)
C
      IPICK_T_ACTIVE = ICOUNT_TBLK
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
C**********************************************************************
      INTEGER*8 FUNCTION IGET_MY_T_LEN(IT_TTPL,IDIM,N1INT,IMULT)
C
C     find all active T-blocks for all CPU's 
C
C     Written by S. Knecht - Feb 01 2008
C
C     Last revision :
C
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#if defined (VAR_MPI2)
#include "mpif.h"
#include "infpar.h"
#endif
#include "parluci.h"
#include "ipoist8.inc"
#include "ctcc.inc"
      DIMENSION IT_TTPL(*)
      INTEGER*8 J_COUNT_LEN, IILENGTH
C
      IGET_MY_T_LEN = 0
      J_COUNT_LEN   = 0
      IILENGTH      = 0
C
C     set total number of integrals
C
      NINT_TOTAL = N1INT + N2ELINT
      NINT_TOTAL = NINT_TOTAL
csk   WRITE(LUWRT,*) ' NINT_TOTAL, N1INT, N2ELINT:', 
csk  &                 NINT_TOTAL, N1INT, N2ELINT
csk   WRITE(LUWRT,*) ' IDIM, IMULT:', 
csk  &                 IDIM, IMULT
csk   CALL IWRTMAMN(IT_TTPL,1,IDIM*IMULT,1,
csk  &              IDIM*IMULT,LUWRT)
C
      IF( IT_SHL .ge. 0 )THEN
        DO JMULT = 1, IMULT
          DO IBLK = 1, IDIM
            IF( IT_TTPL( IBLK + ( IDIM * (JMULT - 1 ) ) ) 
     &          .eq. MYNEW_ID_SM )THEN
              IILENGTH    = NELM_CC(IBLK,IDIM,NINT_TOTAL)
csk           WRITE(LUWRT,*) ' IILENGTH:', IILENGTH
              J_COUNT_LEN = J_COUNT_LEN + IILENGTH
            END IF
          END DO
        END DO
      ELSE IF( SPLIT_IJKL )THEN
        DO JMULT = 1, IMULT
          DO IBLK = 1, IDIM
            IF( IT_TTPL( IBLK + ( IDIM * (JMULT - 1 ) ) ) 
     &          .gt. 0 )THEN
              IILENGTH    = NELM_CC(IBLK,IDIM,NINT_TOTAL)
csk           WRITE(LUWRT,*) ' IILENGTH:', IILENGTH
              J_COUNT_LEN = J_COUNT_LEN + IILENGTH
            END IF
          END DO
        END DO
      END IF
C
      IGET_MY_T_LEN = J_COUNT_LEN
csk   WRITE(LUWRT,*) ' IGET_MY_T_LEN:', IGET_MY_T_LEN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine unfold(ARRAY,VECTOR,IDIM,ISYM,NTRIAEL,IFAC)
***********************************************************************
*  Unfold lower triangle square matrix to full size;
*  I suppose that the input matrix is a one-dim. field with
*  N(N/2 + 1) elements
*  1) write vector to matrix
*  2) unfold considering symmetry (e.g. hermitian or not)
***********************************************************************
*
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
*           parallelization by Stefan Knecht                          *
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
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Read indices from file INDICES
*
      subroutine read_indices(IN,NINTS_I8)
#include "clunit.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "files_r.inc"
      integer*8 NINTS_I8
      dimension IN(4,NINTS_I8)
*
      open(unit=LUSC95,file=NLUSC95_MPI(1:LUSC95_MPIL),
     &       status='OLD',form='FORMATTED')
      rewind(LUSC95)
      do I=1,NINTS_I8,1
        read(LUSC95,'(4I4)') (IN(J,I),J=1,4)
      end do
      close(LUSC95)
*
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WRITE_T_DIAG(T_BUFF,ISPINFREE,ENVIRO,SIGDEN_ROUTE,
     &                        INTSPACE,IPRNT,IOFF_T_IJKL)
C**********************************************************************
C
C     write sorted integrals back on file - driver routine                         
C
C     NOTE: (ij|kl) sorting is GAS-specific
C
C     Written by S. Knecht - Jan 10 2008
C     Last revision :
C
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "symm.inc"
#include "oper.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET, IOFF_T_IJKL
#if defined (INT_STAR8)
      INTEGER(KIND=MPI_OFFSET_KIND)  IWRITE_BYTE
#else
      INTEGER IWRITE_BYTE
#endif
      INTEGER(KIND=MPI_OFFSET_KIND) IS_BYTE8 
      INTEGER(KIND=MPI_OFFSET_KIND) IS_LENGHT8
#endif
#include "parluci.h"
#include "typesz_mpi2.h"
      CHARACTER*6 ENVIRO, SIGDEN_ROUTE
      DIMENSION T_BUFF(*)
C

      IS_BYTE8   = 0
      IS_LENGHT8 = 0
C     local memory
      IDUM = 0
      call memmar(KDUM,  IDUM,    'MARK  ',IDUM,'T_DIAG')
C
CTF   Total number of orbitals and active orbitals (spinors, actually)
CTF   from BOTH sets, barred and unbarred !
      NTOOBTOT = NTOOB + NTOOB2
      NACOBTOT = NACOB + NACOB2
C
      call memmar(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      call memmar(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      call memmar(KLK2  ,NTOOB**2,'ADDL  ',2,'KLK   ')
      call memmar(KLH1DA,NTOOB, 'ADDL  ',2,'KLH1DA')
      call memmar(KLH1DB,NTOOB, 'ADDL  ',2,'KLH1DB')
C
      CALL DZERO(WORK(KLH1DA),NTOOB)
      CALL DZERO(WORK(KLH1DB),NTOOB)
      CALL DZERO(WORK(KLJ)   ,NTOOB**2)
      CALL DZERO(WORK(KLK)   ,NTOOB**2)
      CALL DZERO(WORK(KLK2)  ,NTOOB**2)
C
      if (IHAM12.eq.1.or.IHAM12.eq.2) then
        call geth1dia4(WORK(KLH1DA),WORK(KLH1DB),DUMMY,
     &                   INTSPACE,T_BUFF,ENVIRO,SIGDEN_ROUTE)
      end if
C     obtain 2e integrals needed for diagonal
      if (IHAM12.EQ.2.or.IHAM12.eq.3) then
        call get_dia_2e(WORK(KLJ),WORK(KLK),WORK(KLK2),T_BUFF,0,
     &                  ISPINFREE,NGAS,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
     &                  ISPSF_MO,ENVIRO,SIGDEN_ROUTE,IPRNT)
      end if
C
C     write out to disk
C
#if defined (VAR_MPI2)
C
      IOFFSET     = 0
      IWRITE_BYTE = 0 
      IOFFSET     = IOFF_T_IJKL
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB
#if defined (INT_STAR8)
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
#else
      IWRITE_BYTE = ISIZE_dp * NTOOB
#endif
C     WRITE(LUWRT,*)' writing WORK(KLH1DA): offset + bytes',
C    &                IOFFSET, IWRITE_BYTE
C     CALL WRTMAT(WORK(KLH1DA),1,NTOOB,1,NTOOB)
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,WORK(KLH1DA),
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB
#if defined (INT_STAR8)
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
#else
      IWRITE_BYTE = ISIZE_dp * NTOOB
#endif
C     WRITE(LUWRT,*)' writing WORK(KLH1DB): offset + bytes',
C    &                IOFFSET, IWRITE_BYTE
C     CALL WRTMAT(WORK(KLH1DB),1,NTOOB,1,NTOOB)
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,WORK(KLH1DB),
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB**2
#if defined (INT_STAR8)
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
#else
      IWRITE_BYTE = ISIZE_dp * NTOOB**2
#endif
C     WRITE(LUWRT,*)' writing WORK(KLJ): offset + bytes',
C    &                IOFFSET, IWRITE_BYTE
C     CALL WRTMAT(WORK(KLJ),1,NTOOB**2,1,NTOOB**2)
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,WORK(KLJ),
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB**2
#if defined (INT_STAR8)
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
#else
      IWRITE_BYTE = ISIZE_dp * NTOOB**2
#endif
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,WORK(KLK),
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C     
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB**2
#if defined (INT_STAR8)
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
#else
      IWRITE_BYTE = ISIZE_dp * NTOOB**2
#endif
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,WORK(KLK2),
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C
#endif
C     flush local memory
      call memmar(KDUM,  IDUM,    'FLUSM ',IDUM,'T_DIAG')
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WRITE_IJKL_GAS(T,T_BUFF,INTIMP,ISPINFREE,ISYM_T,
     &                          ISOBEX,ISTR_CA,ISTR_CB,ISTR_AA,
     &                          ISTR_AB,NASHT,F1,MDIASV,IBOSYM, 
     &                          LABEL,CIRUN,IPRNT,NDIAMINT,ENVIRO,
     &                          ibeig,SIGDEN_ROUTE,INTSPACE)
C**********************************************************************
C
C     write sorted integrals back on file - driver routine                         
C
C     NOTE: (ij|kl) sorting is GAS-specific
C
C     Written by S. Knecht - Jan 10 2008
C     Last revision :
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) IDISP, IOFFSET
      INTEGER(KIND=MPI_OFFSET_KIND) IS_BYTE8 
      INTEGER(KIND=MPI_OFFSET_KIND) IWRITE_MAX,IMOD_WRITE
      INTEGER(KIND=MPI_OFFSET_KIND) IT_TOFFW,IWRITE_ELM
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_ADD
#if defined (INT_STAR8)
      INTEGER(KIND=MPI_OFFSET_KIND) IWRITE_BYTE
#else
      INTEGER IWRITE_BYTE
#endif
#endif
#include "parluci.h"
#include "mxpdim.inc"
#include "symm.inc"
#include "glbbas.inc"
#include "input.inc"
#include "clunit.inc"
#include "oper.inc"
#include "ctcc.inc"
#include "integrals_off.inc"
*
#include "dgroup.h"
      CHARACTER*6 CIRUN
      DIMENSION   T_BUFF(*), T(*), F1(*), ibeig(*)
#if defined (VAR_MPI2)
#include "krmc_shmem.h"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
      real*8,  allocatable :: t_buffx(:)
#endif
!     scratch memory
      real(8), allocatable :: t_unsort(:)
      real(8), allocatable :: t_buff_lowsort(:)
      real(8), allocatable :: hmat(:)
      real(8), allocatable :: h2xy(:)
      real(8), allocatable :: h2yx(:)
      integer, allocatable :: indxb12(:,:)
      integer, allocatable :: indxb34(:,:)
      integer, allocatable :: idxkr12(:,:)
C
      allocate(t_unsort(len_all_int))
C
      if(iiomod_rel.eq.1 .or. myproc.eq.master)then

!       get dimensions for scratch matrices
        call rgeth2_dim(ndim_ij,ndim_kl,ndim_h12_tmp)
 
!       allocate scratch matrices
        allocate(hmat(ndim_h12_tmp))
        allocate(h2xy(ndim_h12_tmp))
        allocate(h2yx(ndim_h12_tmp))

        allocate(indxb12(2,ndim_ij))
        allocate(indxb34(2,ndim_kl))
        allocate(idxkr12(2,ndim_ij))
        
!       get integrals on t_unsort
        call rgeth21(.false.,.true.,dummy,dummy,t_unsort,
     &               hmat,h2xy,h2yx,idxkr12,indxb12,indxb34,
     &               ibeig,iprnt,luwrt)

!       release scratch space
        deallocate(idxkr12)
        deallocate(indxb34)
        deallocate(indxb12)

        deallocate(hmat)
        deallocate(h2xy)
        deallocate(h2yx)
      end if
C     ^ iiomod_rel == 1 or myproc == master

      IF( IIOMOD_REL .eq. 0 )THEN
C
C       LEN_ALL_INT is INTEGER*8: this might not work with 
C       every MPI implementation (IBM?); anyway bcast should 
C       NOT be our first choice...
C
#if defined (VAR_MPI2)
        CALL MPI_BCAST(t_unsort,LEN_ALL_INT,MPI_REAL8,MASTER,
     &                 IXCOMM,IERR)
#endif
      END IF
!
!     write back to disk pre-sorted integrals
      if( lowsrt_ijkl )then
#if defined (VAR_MPI2)
        CALL MPI_FILE_OPEN(MPI_COMM_SELF,"IJKL_PRESORT",
     &       MPI_MODE_CREATE + MPI_MODE_WRONLY,MPI_INFO_NULL,
     &                     IJKL_PRE,IERR)
        IDISP = 0
        CALL MPI_FILE_SET_VIEW(IJKL_PRE,IDISP,MPI_BYTE,MPI_BYTE,
     &                         "native",MPI_INFO_NULL,IERR)
!
!       check for total number of elements - do not write more than 
!       1.0*10^8 elements at each cycle
        IWRITE_MAX  = 100 000 000
        ILOOP_WRITE = 1
        IMOD_WRITE  = 0
        IT_TOFFW    = 1
        IS_BYTE8    = ISIZE_dp
        IWRITE_ELM  = LEN_ALL_INT
        IOFFSET     = 0
!       more integrals than max per cycle
        IF( IWRITE_ELM .gt. IWRITE_MAX )THEN
          ILOOP_WRITE = IWRITE_ELM / IWRITE_MAX
          IMOD_WRITE  = MOD( IWRITE_ELM, IWRITE_MAX )
          IWRITE_ELM  = IWRITE_MAX
        END IF
#if defined (INT_STAR8)
        IWRITE_BYTE = IWRITE_ELM * IS_BYTE8
#else
        XLENALL = REAL(IWRITE_ELM)
        IWRITE_ELM4 = IGIVE_I_B(XLENALL)
        IWRITE_BYTE = IWRITE_ELM4 * ISIZE_dp
#endif
!       bytes to add in each cycle
        IOFF_ADD      = IWRITE_ELM * IS_BYTE8
!
!
        DO IJ = 1, ILOOP_WRITE
           CALL MPI_FILE_WRITE_AT(IJKL_PRE,IOFFSET,t_unsort(IT_TOFFW),
     &                            IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
!          new offsets
           IT_TOFFW    = IT_TOFFW + IWRITE_ELM
           IOFFSET     = IOFFSET  + IOFF_ADD
        ENDDO
!
!       ... remaining elements
        IF( IMOD_WRITE .gt. 0) THEN
!          bytes to add in remaining step
           IOFF_ADD      = IMOD_WRITE * IS_BYTE8
#if defined (INT_STAR8)
           IWRITE_BYTE   = IMOD_WRITE * IS_BYTE8
#else
           XLENALL = REAL(IMOD_WRITE)
           IWRITE_ELM4 = IGIVE_I_B(XLENALL)
           IWRITE_BYTE = IWRITE_ELM4 * ISIZE_dp
#endif
           CALL MPI_FILE_WRITE_AT(IJKL_PRE,IOFFSET,t_unsort(IT_TOFFW),
     &                            IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
           IOFFSET     = IOFFSET  + IOFF_ADD
        END IF
#endif
!       debug print
!       write(6,*) ' unsorted pre-fetched array of integrals'
!       call wrtmatmn8(t_unsort,1,LEN_ALL_INT,1,LEN_ALL_INT,6)
!
!       free memory and close scratch file
        deallocate(t_unsort)
#if defined (VAR_MPI2)
        call mpi_file_close(ijkl_pre,ierr)
#endif
      end if ! low-memory sorting branch
C
      if( lowsrt_ijkl )then
! 
!       allocate new array for all sorted integrals
        allocate(t_buff_lowsort(nallint))
        t_buff_lowsort  = 0
        call dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,
     &                        ISOBEX,t_buff_lowsort,T_BUFF,
     &                        ISTR_CA,ISTR_CB,ISTR_AA,
     &                        ISTR_AB,NASHT,F1,dummy,
     &                        MDIASV,IBOSYM,LABEL,CIRUN,IPRNT)
      else
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
        allocate(t_buffx(nallint))
        t_buffx = 0
#endif
!
        call dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,ISOBEX,
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
     &                        t_buffx,
#else
     &                        T,
#endif
     &                        T_BUFF,ISTR_CA,ISTR_CB,ISTR_AA,
     &                        ISTR_AB,NASHT,F1,t_unsort,
     &                        MDIASV,IBOSYM,LABEL,CIRUN,IPRNT)
      end if
!
      WRITE(LUWRT,'(/A)')
     & '  =========================================='
      WRITE(LUWRT,'(A)')
     & '  Storing sorted integrals on file IJKL_REOD'
      WRITE(LUWRT,'(A/)')
     & '  =========================================='

!
!     WRITE(LUWRT,*) 'N1DIAINT and N2DIAINT are (2)',
!    &                N1DIAINT,N2DIAINT
      if( lowsrt_ijkl )then
        CALL IJKL_OUTPUT( t_buff_lowsort,T_BUFF,MDIASV,
     &                    NDIAMINT,ENVIRO,SIGDEN_ROUTE,ISPINFREE, 
     &                    INTSPACE, IPRNT,N1DIAINT )
      else
        CALL IJKL_OUTPUT( 
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
     &                    t_buffx, 
#else
     &                    T,
#endif
     &                    T_BUFF, MDIASV, NDIAMINT, ENVIRO, 
     &                    SIGDEN_ROUTE, ISPINFREE, INTSPACE, IPRNT,
     &                    N1DIAINT )
      end if
!
!     flush local memory
      if( lowsrt_ijkl )then
        deallocate(t_buff_lowsort)
      else 
        deallocate(t_unsort)
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
        deallocate(t_buffx)
#endif
      end if
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE IJKL_OUTPUT( T, T_BUFF, MDIASV, NDIAMINT, ENVIRO, 
     &                        SIGDEN_ROUTE, ISPINFREE,INTSPACE,IPRNT,
     &                        NTEST_DIA)
C**********************************************************************
C
C     write sorted integrals on file
C
C     NOTE: (ij|kl) sorting is GAS-specific
C
C     Written by S. Knecht - Jan 11 2008
C     Last revision :
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) IDISP, IOFFSET
      INTEGER(KIND=MPI_OFFSET_KIND) IS_BYTE8 
      INTEGER(KIND=MPI_OFFSET_KIND) IWRITE_MAX,IMOD_WRITE
      INTEGER(KIND=MPI_OFFSET_KIND) IT_TOFFW,IWRITE_ELM
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_T_IJKL, IOFF_ADD
#if defined (INT_STAR8)
      INTEGER(KIND=MPI_OFFSET_KIND) IWRITE_BYTE
#else
      INTEGER IWRITE_BYTE
#endif
#endif
#include "mxpdim.inc"
#include "parluci.h"
#include "integrals_off.inc"
#include "typesz_mpi2.h"
#include "ipoist8.inc"
#include "ctcc.inc"
      DIMENSION T_BUFF(*), T(*), MDIASV(*)
      INTEGER*8 NINT_TP_I8, I_START, I_END, I_OUT
C
      NINT_TOT = N1ELINT + N2ELINT
C     WRITE(LUWRT,*) 'NINT_TOT, N1ELINT, N2ELINT',
C    &                NINT_TOT, N1ELINT, N2ELINT
      IWRITE_ELM4 = 0
      IOFF_ADD    = 0
      IS_BYTE8    = 0
      IWRITE_BYTE = 0
      I_START     = 0
      I_END       = 0
      I_OUT       = 0
#if defined (VAR_MPI2)
      IDISP    = 0
      IOFFSET  = 0
C
#if defined (VAR_PFS) && defined (SYS_AIX)
C
C     provide useful I/O performance information on IBMs GPFS
C
      CALL MPI_INFO_CREATE(FILE_INFO_OBJ, IERR)
      CALL MPI_INFO_SET(FILE_INFO_OBJ, "IBM_largeblock_io",
     &                  "true",IERR)
      CALL MPI_FILE_OPEN(MPI_COMM_SELF,NIIJKL_ROD(1:IIJKL_ROD_LL),
     &     MPI_MODE_CREATE + MPI_MODE_WRONLY,FILE_INFO_OBJ,
     &                       IIJKL_ROD,IERR)
C
      CALL MPI_INFO_FREE(FILE_INFO_OBJ, IERR)
#else
      CALL MPI_FILE_OPEN(MPI_COMM_SELF,NIIJKL_ROD(1:IIJKL_ROD_LL),
     &     MPI_MODE_CREATE + MPI_MODE_WRONLY,MPI_INFO_NULL,
     &                       IIJKL_ROD,IERR)
#endif
C
      CALL MPI_FILE_SET_VIEW(IIJKL_ROD,IDISP,MPI_BYTE,MPI_BYTE,
     &                       "native",MPI_INFO_NULL,IERR)
C
C     WRITE(LUWRT,*) ' NTEST_DIA',NTEST_DIA
      IOFFSET     = 0 
      IWRITE_BYTE = ISIZE_int8
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,NALLINT,
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IWRITE_BYTE = ISIZE_my_int
C     WRITE(LUWRT,*) ' NDIAINT',NDIAINT
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,NDIAINT,
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IWRITE_BYTE = ISIZE_my_int
C     WRITE(LUWRT,*) ' N1DIAINT',N1DIAINT
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,N1DIAINT,
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IWRITE_BYTE = ISIZE_my_int
C     WRITE(LUWRT,*) ' N2DIAINT',N2DIAINT
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,N2DIAINT,
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C
C     T in SIGDEN ordered format
C
C     WRITE(LUWRT,*) ' writing T array with # elements = ',NALLINT
C     WRITE(LUWRT,*) ' # elements in bytes/records = ',
C    &                 IWRITE_BYTE, NALLINT
csk   WRITE(LUWRT,*) ' writing T 100 elements starting at 1 000 000'
csk   IOFF_T_WRITE = 100 000 0
csk   CALL WRTMATMN8(T(IOFF_T_WRITE),1,100,1,100,LUWRT)
C     CALL WRTMATMN8(T,1,NALLINT,1,NALLINT,LUWRT)
C
C     new offset - starting point for T
C
      IOFFSET     = IOFFSET + IWRITE_BYTE
C
C     check for total number of elements - do not write more than 
C     1.0*10⁸ elements at each cycle
C
      IWRITE_MAX  = 100 000 000
      ILOOP_WRITE = 1
      IMOD_WRITE  = 0
      IT_TOFFW    = 1
      IS_BYTE8    = ISIZE_dp
      IWRITE_ELM  = NALLINT
C     more integrals than max per cycle
      IF( NALLINT .gt. IWRITE_MAX )THEN
        ILOOP_WRITE = NALLINT / IWRITE_MAX
        IMOD_WRITE  = MOD( NALLINT, IWRITE_MAX )
        IWRITE_ELM  = IWRITE_MAX
      END IF
#if defined (INT_STAR8)
      IWRITE_BYTE = IWRITE_ELM * IS_BYTE8
#else
      XLENALL = REAL(IWRITE_ELM)
      IWRITE_ELM4 = IGIVE_I_B(XLENALL)
      IWRITE_BYTE = IWRITE_ELM4 * ISIZE_dp
#endif
C     bytes to add in each cycle
      IOFF_ADD      = IWRITE_ELM * IS_BYTE8
C
      DO IJ = 1, ILOOP_WRITE
C        WRITE(LUWRT,*) ' loop index IJ',IJ
C        CALL WRTMATMN8(T(IT_TOFFW),1,IWRITE_ELM,1,
C    &   IWRITE_ELM,LUWRT)
         CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,T(IT_TOFFW),
     &                          IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C        new offsets
         IT_TOFFW    = IT_TOFFW + IWRITE_ELM
         IOFFSET     = IOFFSET  + IOFF_ADD
      ENDDO
C
C     ... remaining elements
C
      IF( IMOD_WRITE .gt. 0) THEN
C        bytes to add in remaining step
         IOFF_ADD      = IMOD_WRITE * IS_BYTE8
#if defined (INT_STAR8)
          IWRITE_BYTE   = IMOD_WRITE * IS_BYTE8
#else
          XLENALL = REAL(IMOD_WRITE)
          IWRITE_ELM4 = IGIVE_I_B(XLENALL)
          IWRITE_BYTE = IWRITE_ELM4 * ISIZE_dp
#endif
          CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,T(IT_TOFFW),
     &                           IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
          IOFFSET     = IOFFSET  + IOFF_ADD
      END IF
C
C     T diagonal in GASDIAT_REL ordered format
C
      IOFF_T_IJKL = 0
      IOFF_T_IJKL = IOFFSET
C     WRITE(LUWRT,*) ' IOFFSET, IWRITE_BYTE, IOFF_T_IJKL',
C    &                 IOFFSET, IWRITE_BYTE, IOFF_T_IJKL
C
      CALL WRITE_T_DIAG(T,ISPINFREE,ENVIRO,SIGDEN_ROUTE,
     &                  INTSPACE,IPRNT,IOFF_T_IJKL)
C
      CALL MPI_FILE_CLOSE(IIJKL_ROD,IERR)
C
#endif
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE READ_IJKL_GAS(T,T_BUFF,T_BUFF_D,IT_TTPL,
     &                         MDIASV,NDIAMINT,IDIM_TYP,N1_INT)
C**********************************************************************
C
C     read sorted integrals from file
C
C     NOTE: (ij|kl) sorting is GAS-specific
C
C     Written by S. Knecht - Jan 11 2008
C     Last revision :
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) IDISP,IOFFSET,IREAD_MAX,IMOD_READ
      INTEGER(KIND=MPI_OFFSET_KIND) IS_BYTE8,IOFFSET_T, IT_TOFFW
      INTEGER(KIND=MPI_OFFSET_KIND) IREAD_ELM,IOFF_ADD
#if defined (INT_STAR8)
      INTEGER(KIND=MPI_OFFSET_KIND)  IREAD_BYTE
#else
      INTEGER IREAD_BYTE
#endif
#endif
#include "parluci.h"
#include "mxpdim.inc"
#include "ipoist8.inc"
#include "ctcc.inc"
#include "integrals_off.inc"
#include "orbinp.inc"
#include "typesz_mpi2.h"
      DIMENSION T_BUFF(*),T(NALLINT),MDIASV(*),IT_TTPL(*)
      DIMENSION T_BUFF_D(*)
      INTEGER*8 NINT_TP_I8, I_START, I_END, I_IN
C
      NINT_TOT = N1_INT + N2ELINT
      IWRITE_ELM4 = 0
      IOFF_ADD    = 0
      I_START  = 0
      I_END    = 0
      I_OUT    = 0
#if defined (VAR_MPI2)
C
      IDISP    = 0
      IOFFSET  = 0
C
#if defined (VAR_PFS) && defined (SYS_AIX)
C
C     provide useful information on IBMs GPFS to enhance 
C     I/O performance
C
      CALL MPI_INFO_CREATE(FILE_INFO_OBJ, IERR)
      CALL MPI_INFO_SET(FILE_INFO_OBJ, "IBM_largeblock_io",
     &                  "true",IERR)
      CALL MPI_FILE_OPEN(MPI_COMM_SELF,NIIJKL_ROD(1:IIJKL_ROD_LL),
     &                   MPI_MODE_RDONLY,FILE_INFO_OBJ,
     &                   IIJKL_ROD,IERR)
C
      CALL MPI_INFO_FREE(FILE_INFO_OBJ, IERR)
#else

      CALL MPI_FILE_OPEN(MPI_COMM_SELF,NIIJKL_ROD(1:IIJKL_ROD_LL),
     &                   MPI_MODE_RDONLY,MPI_INFO_NULL,
     &                   IIJKL_ROD,IERR)
C
#endif
      CALL MPI_FILE_SET_VIEW(IIJKL_ROD,IDISP,MPI_BYTE,MPI_BYTE,
     &                       "native",MPI_INFO_NULL,IERR)
C
      IOFFSET     = 0
      IREAD_BYTE  = ISIZE_int8
      CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,NALLINT,
     &                     IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
C     WRITE(LUWRT,*) ' read-in at ',IOFFSET
C     WRITE(LUWRT,*) ' read-in result: NALLINT ',NALLINT
      IOFFSET     = IOFFSET + IREAD_BYTE
      IREAD_BYTE  = ISIZE_my_int
      CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,NDIAINT,
     &                       IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
C     WRITE(LUWRT,*) ' read-in at ',IOFFSET
C     WRITE(LUWRT,*) ' read-in result: NDIAINT ',NDIAINT
      IOFFSET     = IOFFSET + IREAD_BYTE
      IREAD_BYTE  = ISIZE_my_int
      CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,N1DIAINT,
     &                       IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
C     WRITE(LUWRT,*) ' read-in at ',IOFFSET
C     WRITE(LUWRT,*) ' read-in result: N1DIAINT',N1DIAINT
      IOFFSET     = IOFFSET + IREAD_BYTE
      IREAD_BYTE  = ISIZE_my_int
      CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,N2DIAINT,
     &                       IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
C     WRITE(LUWRT,*) ' read-in at ',IOFFSET
C     WRITE(LUWRT,*) ' read-in result: N2DIAINT',N2DIAINT
C
C     T in SIGDEN ordered format
C
      IOFFSET     = IOFFSET + IREAD_BYTE
      IOFFSET_T   = IOFFSET
C     WRITE(LUWRT,*) ' (T_BUFF) read-in at ',IOFFSET
C
C     read-in integrals to T_BUFF according to type list IT_TTPL
C
C     WRITE(LUWRT,*) ' calling READ_INT_TTYPES: NSP...,N1...',
C    &                 IDIM_TYP, N1_INT
C     CALL IWRTMAMN(IT_TTPL,1,IDIM_TYP,1,IDIM_TYP,LUWRT)
C     ... shared memory modus for integrals / reduced integral list
C         modus
      IF( IT_SHL .ge. 0 .or. SPLIT_IJKL )THEN
         CALL READ_INT_TTYPES(T_BUFF,IT_TTPL,IDIM_TYP,N1_INT,
     &                        IOFFSET_T)
      ELSE
C     ... old-school way of integral handling
C
C       check for total number of elements - do not read more than 
C       1.0*10⁸ elements at each cycle
C
        IREAD_MAX  = 100 000 000
        ILOOP_READ = 1
        IMOD_READ  = 0
        IT_TOFFW   = 1
        IS_BYTE8   = ISIZE_dp
        IREAD_ELM = NALLINT
C       more integrals than max per cycle
        IF( NALLINT .gt. IREAD_MAX )THEN
          ILOOP_READ = NALLINT / IREAD_MAX
          IMOD_READ  = MOD( NALLINT, IREAD_MAX )
          IREAD_ELM  = IREAD_MAX
        END IF
#if defined (INT_STAR8)
        IREAD_BYTE = IREAD_ELM * IS_BYTE8
#else
        XLENALL    = REAL(IREAD_ELM)
        IREAD_ELM4 = IGIVE_I_B(XLENALL)
        IREAD_BYTE = IREAD_ELM4 * ISIZE_dp
#endif
C       bytes to add in each cycle
        IOFF_ADD      = IREAD_ELM * IS_BYTE8
C       WRITE(LUWRT,*) ' IREAD_ELM,IREAD_BYTE,IT_TOFFW',
C    &                   IREAD_ELM,IREAD_BYTE,IT_TOFFW,
C    &                   IMOD_READ, IOFFSET
C
C       WRITE(LUWRT,*) ' ILOOP_READ',ILOOP_READ
        DO IJ = 1, ILOOP_READ
           CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,T(IT_TOFFW),
     &                           IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
C          WRITE(LUWRT,*) ' read cycle IJ',IJ
C          new offsets
           IT_TOFFW    = IT_TOFFW + IREAD_ELM
           IOFFSET     = IOFFSET  + IOFF_ADD
        ENDDO
C
C       ... remaining elements
C
        IF( IMOD_READ .gt. 0) THEN
C         bytes to add in remaining step
          IOFF_ADD   = IMOD_READ * IS_BYTE8
#if defined (INT_STAR8)
          IREAD_BYTE = IMOD_READ * IS_BYTE8
#else
          XLENALL    = REAL(IMOD_READ)
          IREAD_ELM4 = IGIVE_I_B(XLENALL)
          IREAD_BYTE = IREAD_ELM4 * ISIZE_dp
#endif
          CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,T(IT_TOFFW),
     &                           IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
          IOFFSET     = IOFFSET  + IOFF_ADD
        END IF
      
      END IF
C     WRITE(LUWRT,*) ' printing T_BUFF with # elements =',
C    &                                          NALLINT
C     WRITE(LUWRT,*) ' printing first 100 elements'
C     CALL WRTMATMN(T,1,100,1,100,LUWRT)
C
C     T diagonal in GASDIAT_REL ordered format
C
      IS_BYTE8    = ISIZE_dp
      IOFFSET     = IOFFSET_T  + ( NALLINT * IS_BYTE8 )
#if defined (INT_STAR8)
      IREAD_BYTE  = MY_T_LEN_D * IS_BYTE8
#else
      XLENTD        = REAL(MY_T_LEN_D)
      MY_T_LEN_D_I4 = IGIVE_I_B(XLENTD)
      IREAD_BYTE    = MY_T_LEN_D_I4 * ISIZE_dp
#endif
C     WRITE(LUWRT,*) ' read-in T_BUFF_D at',IOFFSET
C     WRITE(LUWRT,*) ' read-in T_BUFF_D # bytes',IREAD_BYTE
      CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,T_BUFF_D,
     &                      IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
C
C     WRITE(LUWRT,*) ' printing T_BUFF_D with # elements =',
C    &                                          MY_T_LEN_D
C     WRITE(LUWRT,*) ' printing first 100 elements'
C     CALL WRTMATMN8(T_BUFF_D,1,100,1,100,LUWRT)
C
      CALL MPI_FILE_CLOSE(IIJKL_ROD,IERR)
C
#endif
C
      END
#if defined (VAR_MPI2)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CHECK_CONSIST( ICHECK_VEC, NDIM )
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     check if logical is consistent on each CPU
C
C     Written by S. Knecht - Feb 22 2008
C
C     Last revision :
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mpif.h"
#include "infpar.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION ICHECK_VEC(*)
C
      ICHECK_SCR   = ICHECK_VEC(1)
      ICHECK_SCR_2 = 0
      DO II = 1, NDIM
         ICHECK_SCR_2 = ICHECK_VEC(II)
         IF( ICHECK_SCR_2 .ne. ICHECK_SCR )THEN
           WRITE(LUWRT,*) ' *** ERROR in CHECK_CONSIST: unconsistency'//
     &                    ' detected ***:', ICHECK_SCR, ICHECK_SCR_2
           CALL QUIT(' *** ERROR in CHECK_CONSIST: unconsistency 
     & detected *** ' )
          END IF
      END DO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DIST_XPROP(XPROP,LUPROP_MAT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     distribute property integrals from MASTER to node-MASTER.
C
C     Written by S. Knecht - Oct 2008
C
C     Last revision :
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mpif.h"
#include "infpar.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
C information for property run
#include "krciprop.h"
#include "dcbxpr.h"
      CHARACTER FILE_INFO_GROUPSZ*4,XPROP_INT*13,FILELAB_REOD*8,RODLAB*2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_XPRP
      DIMENSION XPROP(*)
C
      RODLAB    = 'RO'
      IOFF_XPRP = 0
C
C     distribute from MASTER to node-MASTER(s)
C     ---------------------------------------
C
C     file info object - provide useful hints for the MPI implementation
C
      CALL MPI_INFO_CREATE(FILE_INFO_OBJ, IERR)
C
C     ... number of CPUs sharing the following MPI-I/O file
      WRITE (FILE_INFO_GROUPSZ,'(I4)') NEWCOMM_PROC
      CALL MPI_INFO_SET(FILE_INFO_OBJ,"nb_proc",FILE_INFO_GROUPSZ,IERR)
C
#if defined (VAR_PFS)
C
C     special information on IBMs GPFS to enhance I/O performance
C
      CALL MPI_INFO_SET(FILE_INFO_OBJ, "IBM_largeblock_io", "true",IERR)
#endif
C
C     ... set correct name
      LXPROP_INT = 13
      IF( MY_GROUPN .lt. 10 ) THEN    ! MPI ID has one digit
        WRITE (XPROP_INT,'(A10,I1)') 'XPROP_INT.',MY_GROUPN
        LXPROP_INT = 11
      ELSE IF( MY_GROUPN .lt. 100 ) THEN ! MPI ID has two digits
        WRITE (XPROP_INT,'(A10,I2)') 'XPROP_INT.',MY_GROUPN
        LXPROP_INT = 12
      ELSE ! MPI ID has three digits
        WRITE (XPROP_INT,'(A10,I3)') 'XPROP_INT.',MY_GROUPN
      END IF
C
C     ... open file
      CALL MPI_FILE_OPEN(MYNEW_COMM,XPROP_INT(1:LXPROP_INT),
     &     MPI_MODE_DELETE_ON_CLOSE + MPI_MODE_CREATE + MPI_MODE_RDWR,
     &     FILE_INFO_OBJ,ILPRP_X,IERR)
C
      CALL MPI_INFO_FREE(FILE_INFO_OBJ,IERR)
C
C     ... set fileview
      CALL MPI_FILE_SET_VIEW(ILPRP_X,IOFF_XPRP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
C
C     XPROP distribution
      DO IPROP_X = 1, NPROP_KRCI
         LEN_XPRP_INT = IN1ELPRP_KRCI(IPROP_X)
         IF( LEN_XPRP_INT .gt. 0 ) THEN
           IF( MYPROC .eq. MASTER )THEN
             INDXPR  = LPROP_KRCI(IPROP_X)
             WRITE(FILELAB_REOD,'(A6,A2)') PRPNAM(INDXPR)(1:6),RODLAB
             CALL KRCI_PRPFILE(LUPROP_MAT,FILELAB_REOD,XPROP,
     &                         LEN_XPRP_INT,1)
           END IF
           IF( MYNEW_ID .eq. N_MASTER ) THEN
             CALL MPI_BCAST(XPROP,LEN_XPRP_INT,MPI_REAL8,MASTER,
     &                      ICOMM,IERR)
             CALL MPI_FILE_WRITE_AT(ILPRP_X,IOFF_XPRP,XPROP,
     &                              LEN_XPRP_INT,MPI_REAL8,ISTAT,IERR)
C            new offset
             IOFF_XPRP = IOFF_XPRP + LEN_XPRP_INT
           END IF
         END IF
      END DO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_BLOCKFILE_DATA(BLOCKFLE,DIST_COMPLETE,NMPROC_FILE,
     &                              SPLIT_IJKL_FILE,SHARED_M_FILE,
     &                              ICCTOS,RCCTOS,NPARBLOCK,ICWEIGHTF,
     &                              IT_TTPL,LEN_TTPL)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     read distribution data from formatted file KRCI_BLOCKDIST.x.
C
C     Written by S. Knecht - Dec 2008
C
C     Last revision :
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mpif.h"
#include "infpar.h"
#include "parluci.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER   RCCTOS(NUM_BLOCKS2)
      DIMENSION ICCTOS(I_NZERO_LEN_C,I_NZERO_LEN_S),IT_TTPL(LEN_TTPL)
      DIMENSION NPARBLOCK(NUM_BLOCKS),ICWEIGHTF(I_NZERO_LEN_S)
      CHARACTER BLOCKFLE*16
      LOGICAL SPLIT_IJKL_FILE, SHARED_M_FILE, DIST_COMPLETE
C
C     initialize
      NMPROC_FILE  = 0
      ISPLIT_FILE  = 0
      ISHARED_FILE = -1
      SPLIT_IJKL_FILE = .FALSE.
      SHARED_M_FILE   = .FALSE.
      LUKRCIBLK       = 81
      OPEN(LUKRCIBLK,FILE=BLOCKFLE,STATUS='OLD',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL')
      REWIND LUKRCIBLK
C
C     read data from file
      READ(LUKRCIBLK,'(3I5)') NMPROC_FILE,ISPLIT_FILE,ISHARED_FILE
      READ(LUKRCIBLK,'(66I6)') 
     &    ((ICCTOS(I,J),I = 1,I_NZERO_LEN_C),J=1,I_NZERO_LEN_S)
      READ(LUKRCIBLK,'(66I6)') (RCCTOS(I), I = 1,NUM_BLOCKS2)
      READ(LUKRCIBLK,'(66I6)') (NPARBLOCK(I), I = 1,NUM_BLOCKS)
      READ(LUKRCIBLK,'(66I6)') (ICWEIGHTF(I), I = 1,I_NZERO_LEN_S)
      IF( ISPLIT_FILE .gt.0 ) SPLIT_IJKL_FILE = .TRUE.
C
C     read IT_TTPL only if allocated...
      IF( SPLIT_IJKL_FILE .and. .NOT. SPLIT_IJKL ) 
     &    SPLIT_IJKL_FILE = .FALSE.
      IF( ISHARED_FILE .eq. IT_SHL .and. IT_SHL .ge. 0 )
     &    SHARED_M_FILE  = .TRUE.
C
C     read IT_TTPL
      IF( SPLIT_IJKL_FILE .or. SHARED_M_FILE )THEN
        READ(LUKRCIBLK,'(66I6)') (IT_TTPL(I),I=1,LEN_TTPL)
      END IF
      CLOSE(LUKRCIBLK,STATUS='KEEP')
C
C     check whether all information was available and correct
      IF( NMPROC .eq. NMPROC_FILE )THEN
        IF( (SPLIT_IJKL_FILE .eqv. SPLIT_IJKL) .and.
     &      ( ISHARED_FILE .eq. IT_SHL ) ) DIST_COMPLETE = .TRUE. 
      END IF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PUT_BLOCKFILE_DATA(BLOCKFLE,ICCTOS,RCCTOS,NPARBLOCK,
     &                              ICWEIGHTF,IT_TTPL,LEN_TTPL)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     put distribution data to formatted file KRCI_BLOCKDIST.x.
C
C     Written by S. Knecht - Dec 2008
C
C     Last revision :
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mpif.h"
#include "infpar.h"
#include "parluci.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER   RCCTOS(NUM_BLOCKS2)
      DIMENSION ICCTOS(I_NZERO_LEN_C,I_NZERO_LEN_S),IT_TTPL(LEN_TTPL)
      DIMENSION NPARBLOCK(NUM_BLOCKS),ICWEIGHTF(I_NZERO_LEN_S)
      CHARACTER BLOCKFLE*16
      LOGICAL EX
C
C     initialize
      ISPLIT_FILE  = 0
      ISHARED_FILE = -1
      IF( SPLIT_IJKL ) ISPLIT_FILE = 1
      IF( SHARED_M .and. IT_SHL .ge. 0) ISHARED_FILE = IT_SHL
      LUKRCIBLK       = 81
      INQUIRE(FILE=BLOCKFLE,EXIST=EX)
      IF(EX)THEN 
        OPEN(LUKRCIBLK,FILE=BLOCKFLE,STATUS='OLD',FORM='FORMATTED',
     &       ACCESS='SEQUENTIAL')
        CLOSE(LUKRCIBLK,STATUS='DELETE')
      END IF
      OPEN(LUKRCIBLK,FILE=BLOCKFLE,STATUS='UNKNOWN',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL')
C
C     put data to file
      WRITE(LUKRCIBLK,'(3I5)') NMPROC,ISPLIT_FILE,ISHARED_FILE
      WRITE(LUKRCIBLK,'(66I6)') 
     &    ((ICCTOS(I,J),I = 1,I_NZERO_LEN_C),J=1,I_NZERO_LEN_S)
      WRITE(LUKRCIBLK,'(66I6)') (RCCTOS(I), I = 1,NUM_BLOCKS2)
      WRITE(LUKRCIBLK,'(66I6)') (NPARBLOCK(I), I = 1,NUM_BLOCKS)
      WRITE(LUKRCIBLK,'(66I6)') (ICWEIGHTF(I), I = 1,I_NZERO_LEN_S)
      IF( SPLIT_IJKL .or. ISHARED_FILE .ne. -1 )THEN
        WRITE(LUKRCIBLK,'(66I6)') (IT_TTPL(I),I=1,LEN_TTPL)
      END IF
      CLOSE(LUKRCIBLK,STATUS='KEEP')
      END
#endif
