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
      SUBROUTINE WRITE_IJKL_GAS(T,T_BUFF,INTIMP,ISPINFREE,ISYM_T,
     &                          ISOBEX,ISTR_CA,ISTR_CB,ISTR_AA,
     &                          ISTR_AB,NASHT,F1,MDIASV,IBOSYM, 
     &                          LABEL,CIRUN,IPRNT,NDIAMINT)
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
#include "krmc_shmem.h"
#include "mxdim_mpi2.h"
      POINTER (MY_XINT_SCR_PTR, XINT_SCR(IDUMMY_BUFF_SZ))
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
#if defined (VAR_MPI2)
C
C     dynamic allocation of scratch XINT block (only N_MASTER_SM)
C      --> use in "shared memory" MPI runs, MPI info available
C     if IS_DUMMY > 0 and info im[lemented
C
      IF( SHARED_M ) THEN
C
        IF( MYNEW_ID_SM .eq. N_MASTER_SM )THEN
          IS_DUMMY = 0
          CALL GET_MEM_XINT_SCR( IS_DUMMY )
        END IF
C
      END IF
C
#endif
      IF( .NOT. SHARED_M ) THEN
         CALL MEMMAR_I8(KINTI12,LEN_ALL_INT,'ADDL  ',2,'XTIT12')
      END IF
C
C
      IF( SHARED_M ) THEN
         CALL RGETH21_LUCIAREL(.FALSE.,.TRUE.,DUMMY,DUMMY,
     &                         XINT_SCR,IPRNT)
      ELSE
         CALL RGETH21_LUCIAREL(.FALSE.,.TRUE.,DUMMY,DUMMY,
     &                         WORK(KINTI12),IPRNT)
      END IF
C
C
      IF( SHARED_M ) THEN
        call dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,
     &                        ISOBEX,T,T_BUFF,
     &                        ISTR_CA,ISTR_CB,ISTR_AA,
     &                        ISTR_AB,NASHT,F1,XINT_SCR, 
     &                        MDIASV,IBOSYM,LABEL,CIRUN,IPRNT)
C
      ELSE
C
        call dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,
     &                        ISOBEX,T,T_BUFF,
     &                        ISTR_CA,ISTR_CB,ISTR_AA,
     &                        ISTR_AB,NASHT,F1,WORK(KINTI12), 
     &                        MDIASV,IBOSYM,LABEL,CIRUN,IPRNT)
      END IF
C
      WRITE(LUWRT,'(A)') ' ' 
      WRITE(LUWRT,'(A)')
     & '  =========================================='
      WRITE(LUWRT,'(A)')
     & '  Storing sorted integrals on file IJKL_REOD'
      WRITE(LUWRT,'(A)')
     & '  =========================================='
      WRITE(LUWRT,'(A)') ' ' 

C
      CALL IJKL_OUTPUT( T, T_BUFF, MDIASV, NDIAMINT )
C
C     flush local memory
C
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'IJKLRO')
#if defined (VAR_MPI2)
C
C     release scratch memory
C
      IF( SHARED_M ) THEN
        CALL REL_MEM_XINT_SCR( XINT_SCR )
      END IF
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
      SUBROUTINE IJKL_OUTPUT( T, T_BUFF, MDIASV, NDIAMINT )
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
#include "mxpdim.inc"
#include "lucipar_rel.inc"
#include "ctcc.inc"
#include "integrals_off.inc"
      DIMENSION T_BUFF(*), T(NALLINT), MDIASV(*)
      INTEGER*8 NINT_TP_I8, I_START, I_END, I_OUT
C
      NINT_TOT = N1ELINT + N2ELINT
      I_START  = 0
      I_END    = 0
      I_OUT    = 0
C     
C     open integral file
C
csk      WRITE(LUWRT,*) ' NALLINT ',NALLINT
      OPEN(Unit=IIJKL_ROD,File=NIIJKL_ROD(1:IIJKL_ROD_LL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
C
C     write all information that is needed
C
      WRITE(IIJKL_ROD) NALLINT
C
      IF( .NOT. SHARED_M ) THEN
C
        WRITE(IIJKL_ROD) (T(I_OUT), I_OUT = 1, NALLINT)
C
      ELSE
        WRITE(IIJKL_ROD) (T_BUFF(I_OUT), I_OUT = 1, NALLINT)
csk        CALL WRTMATMN8(T_BUFF,1,NALLINT,1,NALLINT,LUWRT)
      ENDIF
C
C     some more information ...
C
csk      WRITE(LUWRT,*) ' NDIAMINT',NDIAMINT
      WRITE(IIJKL_ROD) NDIAINT,N1DIAINT,N2DIAINT
csk      CALL IWRTMAMN8(IPOSDIA,1,NDIAMINT,1,NDIAMINT,LUWRT)
csk      CALL IWRTMAMN(MDIASV,1,4*NDIAMINT,1,4*NDIAMINT,LUWRT)
      WRITE(IIJKL_ROD) (IPOSDIA(I),I=1,  NDIAMINT)
      WRITE(IIJKL_ROD) (MDIASV(I) ,I=1,4*NDIAMINT)
C
C     close
C
      CLOSE(UNIT=IIJKL_ROD,STATUS='KEEP')
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE READ_IJKL_GAS( T, T_BUFF, MDIASV, NDIAMINT )
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
#include "mxpdim.inc"
#include "lucipar_rel.inc"
#include "ctcc.inc"
#include "integrals_off.inc"
      DIMENSION T_BUFF(*), T(NALLINT), MDIASV(*)
      INTEGER*8 NINT_TP_I8, I_START, I_END, I_IN
C
      NINT_TOT = N1ELINT + N2ELINT
      I_START  = 0
      I_END    = 0
      I_OUT    = 0
C
C     open integral file
C
      OPEN(Unit=IIJKL_ROD,File=NIIJKL_ROD(1:IIJKL_ROD_LL),
     &     Status='OLD',Form='UNFORMATTED')
      REWIND(IIJKL_ROD)
C
C     read all integrals in
C
      READ(IIJKL_ROD) NALLINT
csk      WRITE(LUWRT,*) ' NALLINT ',NALLINT
C
      IF( .NOT. SHARED_M ) THEN
        READ(IIJKL_ROD) (T(I_IN), I_IN = 1, NALLINT)
      ELSE
        READ(IIJKL_ROD) (T_BUFF(I_IN), I_IN = 1, NALLINT)
      ENDIF
C
csk        CALL WRTMATMN8 (T_BUFF,1,NALLINT,1,NALLINT,LUWRT)
C
C     some more information ...
C
      READ(IIJKL_ROD) NDIAINT,N1DIAINT,N2DIAINT
csk      WRITE(LUWRT,*) ' NDIAMINT',NDIAMINT
      READ(IIJKL_ROD) (IPOSDIA(I),I=1,  NDIAMINT)
      READ(IIJKL_ROD) (MDIASV(I) ,I=1,4*NDIAMINT)
csk      CALL IWRTMAMN8(IPOSDIA,1,NDIAMINT,1,NDIAMINT,LUWRT)
csk      CALL IWRTMAMN(MDIASV,1,4*NDIAMINT,1,4*NDIAMINT,LUWRT)
C
C     close
C
      CLOSE(UNIT=IIJKL_ROD,STATUS='KEEP')
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
#endif
