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
*. offset when work is integer array
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
      return
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
#include "lucipar_rel.inc"
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
      INTEGER*8 FUNCTION IGET_MY_T_LEN(IT_TTPL,IMULT)
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
#include "lucipar_rel.inc"
#include "ctcc.inc"
      DIMENSION IT_TTPL(*)
      INTEGER*8 J_COUNT_LEN, ITLTND, IILENGTH
C
      IGET_MY_T_LEN = 0
      J_COUNT_LEN   = 0
      IILENGTH      = 0
C
      DO JMULT = 1, IMULT
        DO IBLK = 1, NSPOBEX_TP
          IF( IT_TTPL( IBLK + ( NSPOBEX_TP * (JMULT - 1 ) ) ) 
     &        .eq. MYNEW_ID_SM )THEN
            IILENGTH    = NELM_CC(IBLK,NSPOBEX_TP,NINT_TOTAL)
            J_COUNT_LEN = J_COUNT_LEN + IILENGTH
          END IF
        END DO
      END DO
C
      IGET_MY_T_LEN = J_COUNT_LEN
csk      WRITE(LUWRT,*) ' IGET_MY_T_LEN:', IGET_MY_T_LEN
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
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
#include "lucipar_rel.inc"
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
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET, IWRITE_BYTE, IOFF_T_IJKL
      INTEGER(KIND=MPI_OFFSET_KIND) IS_BYTE8, IS_LENGHT8
#endif
#include "lucipar_rel.inc"
#include "typesz_mpi2.h"
      CHARACTER*6 ENVIRO, SIGDEN_ROUTE
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
#if defined(VAR_MPI2)
C
      IOFFSET     = 0
      IWRITE_BYTE = 0 
      IOFFSET     = IOFF_T_IJKL
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,WORK(KLH1DA),
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,WORK(KLH1DB),
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB**2
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,WORK(KLJ),
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB**2
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
      CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,WORK(KLK),
     &                       IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C     
      IOFFSET     = IOFFSET + IWRITE_BYTE
      IS_BYTE8    = ISIZE_dp
      IS_LENGHT8  = NTOOB**2
      IWRITE_BYTE = IS_BYTE8 * IS_LENGHT8
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
     &                          SIGDEN_ROUTE,INTSPACE)
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
#include "ipoist8.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "lucipar_rel.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
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
      DIMENSION   T_BUFF(*), T(NALLINT), F1(*)
C
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'IJKLRO')
C
      CALL MEMMAR_I8(KINTI12,LEN_ALL_INT,'ADDL  ',2,'XTIT12')
C
      IF( IIOMOD_REL .eq. 1 .or. MYPROC .eq. MASTER )THEN
C
        CALL RGETH21_LUCIAREL(.FALSE.,.TRUE.,DUMMY,DUMMY,
     &                         WORK(KINTI12),IPRNT)
      END IF
C     ^ IIOMOD_REL == 1 or MYPROC ==. MASTER
      IF( IIOMOD_REL .eq. 0 )THEN
C
C       LEN_ALL_INT is INTEGER*8: this might not work with 
C       every MPI implementation (IBM?); anyway bcast should 
C       NOT be our first choice...
C
#if defined (VAR_MPI2)
        CALL MPI_BCAST(WORK(KINTI12),LEN_ALL_INT,MPI_REAL8,
     &                 MASTER,IXCOMM,IERR)
#endif
C
      END IF
C
      call dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,
     &                      ISOBEX,T,T_BUFF,
     &                      ISTR_CA,ISTR_CB,ISTR_AA,
     &                      ISTR_AB,NASHT,F1,WORK(KINTI12),
     &                      MDIASV,IBOSYM,LABEL,CIRUN,IPRNT)
C
      WRITE(LUWRT,'(A)') ' '
      WRITE(LUWRT,'(A)')
     & '  =========================================='
      WRITE(LUWRT,'(A)')
     & '  Storing sorted integrals on file IJKL_REOD'
      IF( IT_SHL .eq. - 2 ) THEN
        WRITE(LUWRT,'(A)')
     & '     file type  : (ij|kl)-sharing format    '
      ELSE
        WRITE(LUWRT,'(A)')
     & '     file type  : standard (ij|kl) format   '
      ENDIF
      WRITE(LUWRT,'(A)')
     & '  =========================================='
      WRITE(LUWRT,'(A)') ' ' 

C
      CALL IJKL_OUTPUT( T, T_BUFF, MDIASV, NDIAMINT, ENVIRO, 
     &                  SIGDEN_ROUTE, ISPINFREE, INTSPACE, IPRNT )
C
C     flush local memory
C
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'IJKLRO')
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE IJKL_OUTPUT( T, T_BUFF, MDIASV, NDIAMINT, ENVIRO, 
     &                        SIGDEN_ROUTE, ISPINFREE, INTSPACE, IPRNT)
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
      INTEGER(KIND=MPI_OFFSET_KIND) IS_BYTE8, IWRITE_BYTE
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_T_IJKL
#endif
#include "mxpdim.inc"
#include "lucipar_rel.inc"
#include "ctcc.inc"
#include "integrals_off.inc"
#include "typesz_mpi2.h"
      DIMENSION T_BUFF(*), T(NALLINT), MDIASV(*)
      INTEGER*8 NINT_TP_I8, I_START, I_END, I_OUT
C
      NINT_TOT = N1ELINT + N2ELINT
      IS_BYTE8    = 0
      IWRITE_BYTE = 0
      I_START     = 0
      I_END       = 0
      I_OUT       = 0
      IF( IT_SHL .eq. - 2 )THEN
#if defined (VAR_MPI2)
        IDISP    = 0
        IOFFSET  = 0
C
        CALL MPI_FILE_OPEN(MPI_COMM_SELF,NIIJKL_ROD(1:IIJKL_ROD_LL),
     &       MPI_MODE_CREATE + MPI_MODE_WRONLY,MPI_INFO_NULL,
     &                         IIJKL_ROD,IERR)
C
        CALL MPI_FILE_SET_VIEW(IIJKL_ROD,IDISP,MPI_BYTE,MPI_BYTE,
     &                         "native",MPI_INFO_NULL,IERR)
C
        IOFFSET     = 0 
        IWRITE_BYTE = ISIZE_int8
        CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,NALLINT,
     &                         IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
        IOFFSET     = IOFFSET + IWRITE_BYTE
        IWRITE_BYTE = ISIZE_my_int
        CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,NDIAINT,
     &                         IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
        IOFFSET     = IOFFSET + IWRITE_BYTE
        IWRITE_BYTE = ISIZE_my_int
        CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,N1DIAINT,
     &                         IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
        IOFFSET     = IOFFSET + IWRITE_BYTE
        IWRITE_BYTE = ISIZE_my_int
        CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,N2DIAINT,
     &                         IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
C
C       T in SIGDEN ordered format
C
        IOFFSET     = IOFFSET + IWRITE_BYTE
        IS_BYTE8    = ISIZE_dp
        IWRITE_BYTE = NALLINT * IS_BYTE8
csk        WRITE(LUWRT,*) ' writing T array with # elements = ',NALLINT
csk        WRITE(LUWRT,*) ' # elements in bytes = ',IWRITE_BYTE
csk        CALL WRTMATMN(T,1,NALLINT,1,NALLINT,LUWRT)
        CALL MPI_FILE_WRITE_AT(IIJKL_ROD,IOFFSET,T,
     &                         IWRITE_BYTE,MPI_BYTE,ISTAT,IERR)
csk        WRITE(LUWRT,*) ' write-out at ',IOFFSET
C
C       T diagonal in GASDIAT_REL ordered format
C
        IOFF_T_IJKL = 0
        IOFF_T_IJKL = IOFFSET + IWRITE_BYTE
csk        WRITE(LUWRT,*) ' IOFFSET, IWRITE_BYTE, IOFF_T_IJKL',
csk     &                   IOFFSET, IWRITE_BYTE, IOFF_T_IJKL
C
        CALL WRITE_T_DIAG(T_BUFF,ISPINFREE,ENVIRO,SIGDEN_ROUTE,
     &                    INTSPACE,IPRNT,IOFF_T_IJKL)
C
        CALL MPI_FILE_CLOSE(IIJKL_ROD,IERR)
C
#endif
      ELSE
C
C       open integral file
C
        OPEN(Unit=IIJKL_ROD,File=NIIJKL_ROD(1:IIJKL_ROD_LL),
     &       Status='UNKNOWN',Form='UNFORMATTED')
C
C       write all information that is needed
C
        WRITE(IIJKL_ROD) NALLINT,NDIAINT,N1DIAINT,N2DIAINT
csk      CALL IWRTMAMN8(IPOSDIA,1,NDIAMINT,1,NDIAMINT,LUWRT)
csk      CALL IWRTMAMN(MDIASV,1,4*NDIAMINT,1,4*NDIAMINT,LUWRT)
        WRITE(IIJKL_ROD) (IPOSDIA(I),I=1,  NDIAMINT)
        WRITE(IIJKL_ROD) (MDIASV(I) ,I=1,4*NDIAMINT)
C       ... integrals
        WRITE(IIJKL_ROD) (T(I_OUT), I_OUT = 1, NALLINT)
C
C       close
C
        CLOSE(UNIT=IIJKL_ROD,STATUS='KEEP')
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE READ_IJKL_GAS(T,T_BUFF,T_BUFF_D,IT_TTPL,
     &                         MDIASV,NDIAMINT)
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
      INTEGER(KIND=MPI_OFFSET_KIND) IDISP, IOFFSET
      INTEGER(KIND=MPI_OFFSET_KIND) IS_BYTE8, IREAD_BYTE, IOFFSET_T
#endif
#include "mxpdim.inc"
#include "lucipar_rel.inc"
#include "ctcc.inc"
#include "integrals_off.inc"
#include "orbinp.inc"
#include "typesz_mpi2.h"
      DIMENSION T_BUFF(*),T(NALLINT),MDIASV(*),IT_TTPL(*)
      DIMENSION T_BUFF_D(*)
      INTEGER*8 NINT_TP_I8, I_START, I_END, I_IN
C
      NINT_TOT = N1ELINT + N2ELINT
      I_START  = 0
      I_END    = 0
      I_OUT    = 0
      IF( IT_SHL .ge. 0 )THEN
#if defined (VAR_MPI2)
C
        IDISP    = 0
        IOFFSET  = 0
C
        CALL MPI_FILE_OPEN(MPI_COMM_SELF,NIIJKL_ROD(1:IIJKL_ROD_LL),
     &                     MPI_MODE_RDONLY,MPI_INFO_NULL,
     &                     IIJKL_ROD,IERR)
C
        CALL MPI_FILE_SET_VIEW(IIJKL_ROD,IDISP,MPI_BYTE,MPI_BYTE,
     &                         "native",MPI_INFO_NULL,IERR)
C
        IOFFSET     = 0
        IREAD_BYTE  = ISIZE_int8
        CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,NALLINT,
     &                       IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
csk        WRITE(LUWRT,*) ' read-in at ',IOFFSET
csk        WRITE(LUWRT,*) ' read-in result: NALLINT ',NALLINT
        IOFFSET     = IOFFSET + IREAD_BYTE
        IREAD_BYTE  = ISIZE_my_int
        CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,NDIAINT,
     &                         IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
csk        WRITE(LUWRT,*) ' read-in at ',IOFFSET
csk        WRITE(LUWRT,*) ' read-in result: NDIAINT ',NDIAINT
        IOFFSET     = IOFFSET + IREAD_BYTE
        IREAD_BYTE  = ISIZE_my_int
        CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,N1DIAINT,
     &                         IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
csk        WRITE(LUWRT,*) ' read-in at ',IOFFSET
csk        WRITE(LUWRT,*) ' read-in result: N1DIAINT',N1DIAINT
        IOFFSET     = IOFFSET + IREAD_BYTE
        IREAD_BYTE  = ISIZE_my_int
        CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,N2DIAINT,
     &                         IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
csk        WRITE(LUWRT,*) ' read-in at ',IOFFSET
csk        WRITE(LUWRT,*) ' read-in result: N2DIAINT',N2DIAINT
C
C       T in SIGDEN ordered format
C
        IOFFSET     = IOFFSET + IREAD_BYTE
        IOFFSET_T   = IOFFSET
csk        WRITE(LUWRT,*) ' (T_BUFF) read-in at ',IOFFSET
C
C       raed-in integrals to T_BUFF according to type list IT_TTPL
C
        CALL READ_INT_TTYPES(T_BUFF,IT_TTPL,NSPOBEX_TP,IOFFSET_T)
csk        WRITE(LUWRT,*) ' printing T_BUFF with # elements =',
csk                                                  MY_T_LEN
csk        CALL WRTMATMN(T_BUFF,1,MY_T_LEN,1,MY_T_LEN,LUWRT)
C
C       T diagonal in GASDIAT_REL ordered format
C
        IS_BYTE8    = ISIZE_dp
        IREAD_BYTE  = NALLINT    * IS_BYTE8
        IOFFSET     = IOFFSET_T  + IREAD_BYTE
        IREAD_BYTE  = MY_T_LEN_D * IS_BYTE8
csk        WRITE(LUWRT,*) ' read-in T_BUFF_D at',IOFFSET
csk        WRITE(LUWRT,*) ' read-in T_BUFF_D # bytes',IREAD_BYTE
        CALL MPI_FILE_READ_AT(IIJKL_ROD,IOFFSET,T_BUFF_D,
     &                        IREAD_BYTE,MPI_BYTE,ISTAT,IERR)
C
csk        WRITE(LUWRT,*) ' printing T_BUFF_D with # elements =',
csk                                                  MY_T_LEN_D
csk        CALL WRTMATMN8(T_BUFF_D,1,MY_T_LEN_D,1,MY_T_LEN_D,LUWRT)
C
        CALL MPI_FILE_CLOSE(IIJKL_ROD,IERR)
C
#endif
      ELSE
C
C       open integral file
C
        OPEN(Unit=IIJKL_ROD,File=NIIJKL_ROD(1:IIJKL_ROD_LL),
     &       Status='OLD',Form='UNFORMATTED')
C
        REWIND(IIJKL_ROD)
C
C       read all integrals in
C
        READ(IIJKL_ROD) NALLINT,NDIAINT,N1DIAINT,N2DIAINT
csk      WRITE(LUWRT,*) ' NALLINT ',NALLINT
csk      WRITE(LUWRT,*) ' NDIAMINT',NDIAMINT
        READ(IIJKL_ROD) (IPOSDIA(I),I=1,  NDIAMINT)
        READ(IIJKL_ROD) (MDIASV(I) ,I=1,4*NDIAMINT)
csk      CALL IWRTMAMN8(IPOSDIA,1,NDIAMINT,1,NDIAMINT,LUWRT)
csk      CALL IWRTMAMN(MDIASV,1,4*NDIAMINT,1,4*NDIAMINT,LUWRT)
C
        READ(IIJKL_ROD) (T(I_IN), I_IN = 1, NALLINT)
C
C       close
C
        CLOSE(UNIT=IIJKL_ROD,STATUS='KEEP')
C
      END IF
C
      END
#if defined (VAR_MPI2)
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_IJKL_UPDATE( T, T_BUFF, MDIASV, NDIAMINT,
     &                            UPD_GROUP1, UPD_GROUP2, UPD_GROUP3, 
     &                            UPD_GROUP4 )
C**********************************************************************
C
C     perform (ij|kl) update according to logical flags UPD_GROUPX
C
C     NOTE: (ij|kl) sorting is GAS-specific
C
C     Written by S. Knecht - Jan 14 2008
C     Last revision :
C
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mpif.h"
#include "infpar.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "mxpdim.inc"
#include "lucipar_rel.inc"
#include "ctcc.inc"
#include "integrals_off.inc"
      DIMENSION T_BUFF(*), T(NALLINT), MDIASV(*)
      LOGICAL UPD_GROUP1, UPD_GROUP2, UPD_GROUP3, UPD_GROUP4
C
      IF( UPD_GROUP1 ) THEN
C
C       NODE-MASTER(s) --> slaves
C
        CALL MPI_BCAST(NALLINT, 1,  MPI_INTEGER8, N_MASTER_SM,
     &                 MYNEW_COMM_SM,IERR)
        CALL MPI_BCAST(NDIAINT, 1,my_MPI_INTEGER, N_MASTER_SM,
     &                 MYNEW_COMM_SM,IERR)
        CALL MPI_BCAST(N1DIAINT,1,my_MPI_INTEGER, N_MASTER_SM,
     &                 MYNEW_COMM_SM,IERR)
        CALL MPI_BCAST(N2DIAINT,1,my_MPI_INTEGER, N_MASTER_SM,
     &                 MYNEW_COMM_SM,IERR)
        CALL MPI_BCAST(MDIASV,4*NDIAMINT,my_MPI_INTEGER,N_MASTER_SM,
     &                 MYNEW_COMM_SM,IERR)
        CALL MPI_BCAST(IPOSDIA,NDIAMINT,MPI_INTEGER8,N_MASTER_SM,
     &                 MYNEW_COMM_SM,IERR)
C
      ELSE IF( UPD_GROUP2 ) THEN
C
C       MASTER --> NODE-MASTER(s) + slaves 
C
        CALL MPI_BCAST(NALLINT, 1,  MPI_INTEGER8, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(NDIAINT, 1,my_MPI_INTEGER, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(N1DIAINT,1,my_MPI_INTEGER, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(N2DIAINT,1,my_MPI_INTEGER, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(MDIASV,4*NDIAMINT,my_MPI_INTEGER, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(IPOSDIA,NDIAMINT,MPI_INTEGER8, MASTER,
     &                 MPI_COMM_WORLD,IERR)
C
C       MASTER --> NODE-MASTER(s) broadcast of all (ij|kl)
C
C
C       divide ...
C
        IXINT_COLOR = 1
C
        IF( MYNEW_ID_SM .eq. N_MASTER_SM ) IXINT_COLOR = 4
C
        IXINT_KEY   = MYPROC
C
C       split MPI_COMM_WORLD into N_MASTER_SM and slaves
C       MASTER is always included...
C
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,IXINT_COLOR,
     &                      IXINT_KEY,IXINT_COMM,IERR)
C
        IF( IXINT_COLOR .eq. 1 ) GOTO 200
C
C       communicate integrals
C
        CALL MPI_BCAST(T_BUFF,NALLINT,MPI_REAL8,
     &                 MASTER,IXINT_COMM,IERR)
 200              CONTINUE
C
C       reunite
C
        CALL MPI_COMM_FREE(IXINT_COMM,IERR)
C
      ELSE IF( UPD_GROUP3 ) THEN
C
C       MASTER --> all slaves
C
        CALL MPI_BCAST(NALLINT, 1,  MPI_INTEGER8, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(NDIAINT, 1,my_MPI_INTEGER, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(N1DIAINT,1,my_MPI_INTEGER, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(N2DIAINT,1,my_MPI_INTEGER, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(MDIASV,4*NDIAMINT,my_MPI_INTEGER, MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(IPOSDIA,NDIAMINT,MPI_INTEGER8, MASTER,
     &                 MPI_COMM_WORLD,IERR)
C
C       communicate integrals
C
        CALL MPI_BCAST(T,NALLINT,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
         
C
      ELSE IF( UPD_GROUP4 ) THEN
C
        GOTO 500
      ELSE
        WRITE(LUWRT,*) ' GET_IJKL_UPDATE: ERROR! WHO SHOULD GET AN
     & UPDATE? '
        CALL Abend2( 'error in GET_IJKL_UPDATE! Check logical flags! ')
      ENDIF
C
C     finally ...
C
 500  CONTINUE
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE CHECK_CONSIST( CHECK_VAL, CHECK_VEC, NDIM )
C**********************************************************************
C
C     check if logical is consistent on each CPU
C
C     Written by S. Knecht - Feb 22 2008
C
C     Last revision :
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mpif.h"
#include "infpar.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar_rel.inc"
      LOGICAL CHECK_VAL, CHECK_VEC(*), CHECK_SCR
C
      CHECK_SCR = .FALSE.
C
csk      WRITE(6,*) ' subroutine CHECK_CONSIST: CHECK_VAL is ',   
csk     &             CHECK_VAL
C
      IF( CHECK_VAL ) THEN
        DO II = 1, NDIM
          CHECK_SCR = CHECK_VEC(II)
          IF( .NOT. CHECK_SCR )THEN
             WRITE(LUWRT,*) ' subroutine CHECK_CONSIST: unconsistency 
     &   detected ', CHECK_SCR
            CALL Abend2(' subroutine CHECK_CONSIST: unconsistency 
     &   detected' )
          END IF
        END DO
      ELSE
        DO II = 1, NDIM
          CHECK_SCR = CHECK_VEC(II)
          IF( CHECK_SCR )THEN
             WRITE(LUWRT,*) ' subroutine CHECK_CONSIST: unconsistency 
     &   detected ', CHECK_SCR
            CALL Abend2(' subroutine CHECK_CONSIST: unconsistency 
     &   detected' )
          END IF
        END DO
      END IF
C
      END
#endif
