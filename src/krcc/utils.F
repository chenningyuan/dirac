      SUBROUTINE REWINO( LU )
C
C REWIND SEQ FILE LU WITH FASTIO ROUTINES
C
C?    WRITE(6,*) ' TO REWIND FILE ',LU
      IDUM = 1
C     CALL SQFILE(LU,5,IDUM,IDUM)
      REWIND (LU)
C?    WRITE(6,*) ' FILE REWOUND '
C
      RETURN
      END
*
      SUBROUTINE VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
*
* Read vector VEC of length LENGTH from discfile LU in
* standard LUCIA form, as for example written by VEC_TO_DISC
*
* Jeppe Olsen, March 2000
*
#include "implicit.inc"
*. Input
      DIMENSION VEC(LENGTH)
*
      IF(IREW.EQ.1) THEN
        CALL REWINO(LU)
      END IF
      CALL IFRMDS(LENGTH2,1,LBLK,LU)
      CALL FRMDSC(VEC,LENGTH,LBLK,LU,IAMZERO,IAMPACKED)
*. Skip the -1 at end
      CALL IFRMDS(IONEM,1,LBLK,LU)
*
      RETURN
      END
*
      SUBROUTINE VEC_TO_DISC(VEC,LENGTH,IREW,LBLK,LU)
*
* Write vector VEC of length LENGTH to discfile LU in
* standard LUCIA form
*
* 1 : Length as written by ITODSC
* 2 : The vector as written by TODSC
* 3 : End of vector, -1 written by ITODSC
*
* Jeppe Olsen, March 2000
*
#include "implicit.inc"
*. Input
      DIMENSION VEC(LENGTH)
*
      IF(IREW.EQ.1) THEN
        CALL REWINO(LU)
      END IF
      CALL ITODS(LENGTH,1,LBLK,LU)
      CALL TODSC(VEC,LENGTH,LBLK,LU)
      IONEM = -1
      CALL ITODS(IONEM,1,LBLK,LU)
*
      RETURN
      END
*
      SUBROUTINE WRT_SPOX_TP_CC_KRCC(IEX_TP,NEX_TP)
*
* Print types of spin-orbital excitations
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
      INTEGER IEX_TP(4*NGAS,NEX_TP)
*
      WRITE(6,*)
      WRITE(6,*) ' ***************************************** '
      WRITE(6,*) ' Information about spinorbital excitations '
      WRITE(6,*) ' ***************************************** '
      WRITE(6,*)
*
      DO JEX_TP = 1, NEX_TP
        WRITE(6,*)
        WRITE(6,*) ' Included spinorbitalexcitation ', JEX_TP
        WRITE(6,'(A,16I4)')
     &  ' Creation of alpha     :',
     &  (IEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)')
     &  ' Creation of beta      :',
     &  (IEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)')
     &  ' Annihilation of alpha :',
     &  (IEX_TP(I+2*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)')
     &  ' Annihilation of beta  :',
     &  (IEX_TP(I+3*NGAS,JEX_TP),I=1,NGAS)
      END DO
*
      RETURN
      END
*
      SUBROUTINE ISTVC3(IVEC,IOFF,IVAL,NDIM)
*
* IVEC(IOFF-1+I) = IVAL, I = 1, NDIM
*
* Jeppe Olsen, Oct 2000
*
#include "implicit.inc"
*. Output
      INTEGER IVEC(*)
      DO IELMNT = IOFF, IOFF-1+NDIM
        IVEC(IELMNT) = IVAL
      END DO
*
      RETURN
      END
*
      SUBROUTINE SUM_FDIA_FOR_STR(FSUM,FDIA,NSTR,NEL,IOCC)
*
* F(ISTR) = SUM(IEL) F(IOCC(IEL,ISTR)
*
* Jeppe Olsen, Summer of 99
*
#include "implicit.inc"
*. Input
      DIMENSION FDIA(*)
      INTEGER IOCC(NEL,NSTR)
*. Output
      DIMENSION FSUM(NSTR)
*
      DO ISTR = 1, NSTR
        X = 0.0D0
        DO IEL = 1, NEL
          X = X + FDIA(IOCC(IEL,ISTR))
        END DO
        FSUM(ISTR) = X
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' First two elements of Fdia in Orb basis '
        CALL WRTMAT(FDIA,1,2,1,2)
        WRITE(6,*) ' Occ of Strings '
        CALL IWRTMA(IOCC,NEL,NSTR,NEL,NSTR)
        WRITE(6,*) ' FSUM in Str basis '
        CALL WRTMAT(FSUM,1,NSTR,1,NSTR)
      END IF

*
      RETURN
      END
*
      SUBROUTINE COP_SYMMAT(AIN,AOUT,I1OUT,NOUT)
*
* Extract symmetric packed submatrix AOUT from symmetric
* packed submatrix AIN.
*
*. First row/column in output matrix is I1OUT, and
*. number of rows/columns in AOUT is NOUT
*
* symmetric matrix is packed in conventional form  with
* row wise packing
*
* Jeppe Olsen, March 4 2000
*
#include "implicit.inc"
*. Input
      DIMENSION AIN(*)
*. Output
      DIMENSION AOUT(*)
*
      IADD = I1OUT-1
      DO I = 1, NOUT
        DO J = 1, I
          IJ_OUT = I*(I-1)/2 + J
          IJ_IN  = (I+IADD)*(I+IADD-1)/2 + J + IADD
          AOUT(IJ_OUT) = AIN(IJ_IN)
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output matrix from COP_SYMMAT'
        CALL PRSYM(AOUT,NOUT)
      END IF
*
      RETURN
      END
*
      SUBROUTINE COP_SCA_VEC(VECOUT,VECIN,FACTOR,NDIM)
*
* Copy and scale a vector
*
* VECOUT = FACTOR*VECIN
*
#include "implicit.inc"
*
      DIMENSION VECOUT(NDIM),VECIN(NDIM)
*
      DO I=1,NDIM
        VECOUT(I) = FACTOR*VECIN(I)
      END DO
*
      RETURN
      END
