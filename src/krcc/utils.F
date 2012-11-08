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
      SUBROUTINE VEC_FROM_DISC_KRCC(VEC,LENGTH,IREW,LBLK,LU)
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
      CALL IFRMDS_KRCC(LENGTH2,LBLK2,LU)
      IF(LENGTH2.NE.LENGTH.OR.LBLK.NE.LBLK2) THEN
        print*,'LENGTH,LBLK,LENGTH2,LBLK2',LENGTH,LBLK,LENGTH2,LBLK2
        STOP 'CHECK VEC_FROM_DISC_KRCC'
      END IF
C     CALL FRMDSC(VEC,LENGTH,LBLK,LU,IAMZERO,IAMPACKED)
      READ(LU) (VEC(I),I=1,LENGTH)
*. Skip the -1 at end
      IONEM = -1
      CALL IFRMDS_KRCC(IONEM,LBLK,LU)
*
      RETURN
      END
*
      SUBROUTINE ADD_VEC_FROM_DISC_TO_VEC(VEC,LENGTH,IREW,LBLK,LU)
*
* Read and add to vector VEC of length LENGTH from discfile LU in
* standard KRCC form, as for example written by VEC_TO_DISC_KRCC
*
* Lasse 2012
*
#include "implicit.inc"
*. Input
      DIMENSION VEC(LENGTH)
* Scratch
      DIMENSION SCR(LBLK)
*
      IF(IREW.EQ.1) THEN
        CALL REWINO(LU)
      END IF
*
      CALL IFRMDS_KRCC(LENGTH2,LBLK2,LU)
      print*,'LENGTH2,LBLK2',LENGTH2,LBLK2
      IF(LENGTH2.NE.LENGTH.OR.LBLK.NE.LBLK2) THEN
        print*,'LENGTH,LBLK,LENGTH2,LBLK2',LENGTH,LBLK,LENGTH2,LBLK2
        STOP 'CHECK ADD_VEC_FROM_DISC_TO_VEC'
      END IF
*
C     print*,'LENGTH,LBLK',LENGTH,LBLK
C     IF(LENGTH.GT.LBLK) THEN
C       IDIVIDE = LENGTH/LBLK
C       ICHECK = IDIVIDE*LBLK
C       IDIFF = LENGTH - ICHECK
C       ISTART = 1
C       DO I=1,IDIVIDE
C         IEND = ISTART + LBLK
C         READ(LU) (SCR(J),J=1,LBLK)
C         DO J =1,LBLK
C           print*,'VEC(ISTART +J),SCR(J),J',VEC(ISTART +J),SCR(J),J
C           VEC(ISTART +J) = VEC(ISTART +J) + SCR(J)
C         END DO
C         ISTART = IEND + 1
C       END DO
C       IEND = ISTART + IDIFF
C       READ(LU) (SCR(I),I=1,IDIFF)
C       DO J =1,IDIFF
C         VEC(ISTART +J) = VEC(ISTART +J) + SCR(J)
C       END DO
C     ELSE
C       READ(LU) (SCR(I),I=1,LENGTH)
C       DO I =1,LENGTH
C         print*,'VEC(I),SCR(I),I',VEC(I),SCR(I),I
C         VEC(I) = SCR(I)
C         VEC(I) = VEC(I) + SCR(I)
C       END DO
C     END IF
      CALL ADD_FRMDSC(VEC,LENGTH,LBLK,LU)
*. Skip the -1 at end
      IONEM = -1
      CALL IFRMDS_KRCC(IONEM,LBLK,LU)
*
      RETURN
      END
*
      SUBROUTINE IFRMDS_KRCC(NDIM,LBLK,LU)
* Will just read the dimension (NDIM) and the blocksize (LBLK) to check
      IMPLICIT REAL*8(A-H,O-Z)
*
      IF(NDIM.NE.-1) THEN
        READ(LU) NDIM
        READ(LU) LBLK
        print*,'NDIM,LBLK',NDIM,LBLK
      ELSE
        READ(LU) NDIM
      END IF
*
      RETURN
      END
*
      SUBROUTINE ADD_FRMDSC(VEC,NDIM,MBLOCK,IFIL)
C
C     ADD VEC FROM DISK FILE IFIL TO VEC
C
*. Version allowing zero or more blocks 
*
* Lasse 2012
* 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(NDIM)
* Scratch
      DIMENSION SCR(MBLOCK)
*
C     print*,'NDIM,MBLOCK',NDIM,MBLOCK
      IF(NDIM.GT.MBLOCK) THEN
        IDIVIDE = NDIM/MBLOCK
        ICHECK = IDIVIDE*MBLOCK
        IDIFF = NDIM - ICHECK
        ISTART = 1
C       print*,'IDIVIDE,ICHECK,IDIFF',IDIVIDE,ICHECK,IDIFF
        DO I=1,IDIVIDE
          IEND = ISTART + MBLOCK - 1
C         print*,'ISTART,IEND',ISTART,IEND
          READ(IFIL) (SCR(J),J=1,MBLOCK)
          DO J =1,MBLOCK
C           print*,'VEC(ISTART +J),SCR(J),J',VEC(ISTART +J),SCR(J),J
            VEC(ISTART +J -1) = VEC(ISTART +J -1) + SCR(J)
C        print*,'VEC(ISTART +J),SCR(J),J after',VEC(ISTART +J),SCR(J),J
          END DO
          ISTART = IEND + 1
        END DO
        IEND = ISTART + IDIFF - 1
C       print*,'ISTART,IEND',ISTART,IEND
        READ(IFIL) (SCR(I),I=1,IDIFF)
        DO J =1,IDIFF
          VEC(ISTART +J -1) = VEC(ISTART +J -1) + SCR(J)
C         print*,'ISTART +J -1',ISTART +J -1
        END DO
      ELSE
        READ(IFIL) (SCR(I),I=1,NDIM)
        DO I =1,NDIM
C         print*,'VEC(I),SCR(I),I',VEC(I),SCR(I),I
C         VEC(I) = SCR(I)
          VEC(I) = VEC(I) + SCR(I)
C         print*,'VEC(I),SCR(I),I after',VEC(I),SCR(I),I
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE VEC_TO_DISC_KRCC(VEC,LENGTH,IREW,LBLK,LU)
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
* New routine by Lasse 2012 with slightly different format and where the
* block structure is written in the file. This should give more
* possibilities for reading in files.
*
*
#include "implicit.inc"
*. Input
      DIMENSION VEC(LENGTH)
*
      IF(IREW.EQ.1) THEN
        CALL REWINO(LU)
      END IF
      print*,'LENGTH,LBLK',LENGTH,LBLK
      CALL ITODS_KRCC(LENGTH,LBLK,LU)
      CALL TODSC_KRCC(VEC,LENGTH,LBLK,LU)
      IONEM = -1
      CALL ITODS_KRCC(IONEM,LBLK,LU)
*
      RETURN
      END
*
      SUBROUTINE ITODS_KRCC(NDIM,MBLOCK,IFIL)
C Write the length (NDIM) of the file (IFIL) and the blocksize (MBLOCK)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NDIM,MBLOCK
*
      IF(NDIM.GE.0) THEN
        WRITE(IFIL) NDIM
        WRITE(IFIL) MBLOCK
      ELSE
        WRITE(IFIL) NDIM
      END IF
C
      RETURN
      END
*
      SUBROUTINE TODSC_KRCC(VEC,NDIM,MBLOCK,IFIL)
C Write the vector (VEC) of length (NDIM) of the file (IFIL) and the blocksize (MBLOCK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(NDIM)
      INTEGER NDIM,MBLOCK
*
C     print*,'NDIM,MBLOCK',NDIM,MBLOCK
      IF(NDIM.GT.MBLOCK) THEN
        IDIVIDE = NDIM/MBLOCK
        ICHECK = IDIVIDE*MBLOCK
        IDIFF = NDIM - ICHECK
C       print*,'IDIVIDE,ICHECK,IDIFF',IDIVIDE,ICHECK,IDIFF
        ISTART = 1
        DO I=1,IDIVIDE
          IEND = ISTART + MBLOCK - 1
C         print*,'ISTART,IEND',ISTART,IEND
          WRITE(IFIL) (VEC(J),J=ISTART,IEND)
          ISTART = IEND + 1
        END DO
        IEND = ISTART + IDIFF - 1
C       print*,'ISTART,IEND',ISTART,IEND
        WRITE(IFIL) (VEC(I),I=ISTART,IEND)
      ELSE
        WRITE(IFIL) (VEC(I),I=1,NDIM)
      END IF
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
