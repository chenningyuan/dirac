***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine abend(CHARINP)
      implicit none
      character*(*) CHARINP
      write(6,*) CHARINP
      stop 'Abend called'
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine abend1(CHARINP)
      implicit none
      character*(*) CHARINP
      write(6,*) CHARINP
      stop 'Abend1 called'
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine abend2(CHARINP)
      implicit none
      character*(*) CHARINP
      write(6,*) CHARINP
      stop 'Abend2 called'
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine abend3(CHARINP)
      implicit none
      character*(*) CHARINP
      write(6,*) CHARINP
      stop 'Abend3 called'
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine
      subroutine add_info(ECI,FINEIG,NROOT)
      implicit none
      character*4 ECI
      real*8 FINEIG
      integer NROOT
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck bndinv $Revision$
        SUBROUTINE BNDINV(A,EL,N,DETERM,EPSIL,ITEST,NSIZE)
C
C       DOUBLE PRECISION MATRIX INVERSION SUBROUTINE
C       FROM "DLYTAP".
C
C*      DOUBLE PRECISION E,F
C*      DOUBLE PRECISION A,EL,D,DSQRT,C,S,DETERP
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(NSIZE,1),EL(NSIZE,1)
        IF(N.LT.2)GO TO 140
        ISL2=0
        K000FX=2
        IF(ISL2.EQ.0)INDSNL=2
        IF(ISL2.EQ.1)INDSNL=1
C       CALL SLITET(2,INDSNL)
C       CALL OVERFL(K000FX)
C       CALL DVCHK(K000FX)
C
C       SET EL = IDENTITY MATRIX
        DO 30 I=1,N
        DO 10 J=1,N
 10     EL(I,J)=0.0D0
 30     EL(I,I)=1.0D0
C
C       TRIANGULARIZE A, FORM EL
C
        N1=N-1
        M=2
        DO 50 J=1,N1
        DO 45 I=M,N
        IF(A(I,J).EQ.0.0D0)GO TO 45
        D=DSQRT(A(J,J)*A(J,J)+A(I,J)*A(I,J))
        C=A(J,J)/D
        S=A(I,J)/D
 38     DO 39 K=J,N
        D=C*A(J,K)+S*A(I,K)
        A(I,K)=C*A(I,K)-S*A(J,K)
        A(J,K)=D
 39     CONTINUE
        DO 40 K=1,N
        D=C*EL(J,K)+S*EL(I,K)
        EL(I,K)=C*EL(I,K)-S*EL(J,K)
        EL(J,K)=D
 40     CONTINUE
 45     CONTINUE
 50     M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,51),K000FX
C
C       CALCULATE THE DETERMINANT
 51     DETERP=A(1,1)
        DO 52 I=2,N
 52     DETERP=DETERP*A(I,I)
        DETERM=DETERP
C       CALL OVERFL(K000FX)
C       GO TO (140,520,520),K000FX
C
C       IS MATRIX SINGULAR
 520    F=A(1,1)
        E=A(1,1)
        DO 58 I=2,N
        IF(DABS(F).LT.DABS(A(I,I)))F=A(I,I)
        IF(DABS(E).GT.DABS(A(I,I)))E=A(I,I)
 58     CONTINUE
        EPSILP=EPSIL
        IF(EPSILP.LE.0)EPSILP=1.0E-8
        RAT=E/F
        IF(ABS(RAT).LT.EPSILP)GO TO 130
C
C       INVERT TRIANGULAR MATRIX
        J=N
        DO 100 J1=1,N
C       CALL SLITE(2)
        I=J
        ISL2=1
        DO 90 I1=1,J
C       CALL SLITET(2,K000FX)
        IF(ISL2.EQ.0)K000FX=2
        IF(ISL2.EQ.1)K000FX=1
        IF(ISL2.EQ.1)ISL2=0
        GO TO (70,75),K000FX
 70     A(I,J)=1.0D0/A(I,I)
        GO TO 90
 75     KS=I+1
        D=0.0D0
        DO 80 K=KS,J
 80     D=D+A(I,K)*A(K,J)
        A(I,J)=-D/A(I,I)
 90     I=I-1
 100    J=J-1
C       CALL OVERFL(K000FX)
C       GO TO (140,103,103),K000FX

C103    CALL DVCHK(K000FX)
C       GO TO (140,105),K000FX
C
C       PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
 105    M=1
        DO 120 I=1,N
        DO 118 J=1,N
        D=0.0D0
        DO 107 K=M,N
 107    D=D+A(I,K)*EL(K,J)
        EL(I,J)=D
 118    CONTINUE
 120    M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,123,123),K000FX
C
C       RECOPY EL TO A
 123    DO 124 I=1,N
        DO 124 J=1,N
 124    A(I,J)=EL(I,J)
        ITEST=0
C126    IF(INDSNL.EQ.1)CALL SLITE(2)
 126    IF(INDSNL.EQ.1)ISL2=1
        RETURN
C
 130    ITEST=1
        GO TO 126
 140    ITEST=-1
        GO TO 126
        END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine center
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck readin $Revision$
      SUBROUTINE CHAR_TO_INTEGER(CHAR_X,INT_X)

*     Convert an integer number given as a string into
*     an integer variable

      CHARACTER*72 CHAR_X
      INTEGER      INT_X
      CHARACTER*12 STRING

      CALL STDFMT(CHAR_X,STRING)
      CALL RIGHTAD(STRING)
      READ (STRING,'(I12)',IOSTAT=IOS,ERR=100,END=200) INT_X
      RETURN

100   WRITE (*,*)
      WRITE (*,*) ' *** ERROR IN INPUT ***'
      WRITE (*,*) ' THE PROGRAM EXPECTED AN INTEGER NUMBER'
      WRITE (*,*) ' BUT WAS SUPPLIED WITH A REAL NUMBER OR'
      WRITE (*,*) ' A CHARACTER STRING'
      WRITE (*,*)
      WRITE (*,*) ' FORTRAN I/O STATUS = ',IOS
      WRITE (*,*)
      WRITE (*,*) ' THE LAST LINE PROCESSED WAS:'
      WRITE (*,*) CHAR_X
      WRITE (*,*)
      CALL QUIT(16)

200   WRITE (*,*)
      WRITE (*,*) ' *** ERROR IN INPUT ***'
      WRITE (*,*) ' THE PROGRAM EXPECTED AN INTEGER NUMBER'
      WRITE (*,*) ' BUT WAS SUPPLIED WITH AN EMPTY CHARACTER STRING'
      WRITE (*,*)
      WRITE (*,*) ' FORTRAN I/O STATUS = ',IOS
      WRITE (*,*)
      CALL QUIT(16)

      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck char_to_integer_moluc $Revision$
      SUBROUTINE CHAR_TO_INTEGER_MOLUC(CHAR_X,INT_X)
*
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*72 CHAR_X
*
#include "units.inc"
*
      NTESTL = 0
      if (NTESTL.ge.1) then
        write(6,*) ' INFO FROM CHAR_TO_INTEGER'
        write(6,'(A)') ' CHARACTER INPUT ', CHAR_X
      end if
*
      LU_INTERNAL = LUSC_INT
*
      Open(Unit=LU_INTERNAL,
     &     File='LUSCDL',
     &     Status='NEW',
     &     Form='FORMATTED')
*
      REWIND(LU_INTERNAL)
      WRITE(LU_INTERNAL,'(A)') CHAR_X
      REWIND(LU_INTERNAL)
      READ(LU_INTERNAL,*) INT_X
*
      Close(Unit=LU_INTERNAL,
     &      Status='DELETE')
*
      if (NTESTL.ge.1) then
        write(6,*) ' INTEGER OUTPUT ', INT_X
      end if
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck readin $Revision$
      SUBROUTINE CHAR_TO_REAL(CHAR_X,REAL_X)

*     Convert a real number given as a string into
*     a real variable

      CHARACTER*72 CHAR_X
      REAL*8       REAL_X
      CHARACTER*12 STRING

      CALL STDFMT(CHAR_X,STRING)
      CALL RIGHTAD(STRING)
      READ (STRING,'(F18.10)',IOSTAT=IOS,ERR=100,END=200) REAL_X
      RETURN

100   WRITE (*,*)
      WRITE (*,*) ' *** ERROR IN INPUT ***'
      WRITE (*,*) ' THE PROGRAM EXPECTED A REAL NUMBER'
      WRITE (*,*) ' BUT WAS SUPPLIED WITH AN INTEGER NUMBER OR'
      WRITE (*,*) ' A CHARACTER STRING'
      WRITE (*,*)
      WRITE (*,*) ' FORTRAN I/O STATUS = ',IOS
      WRITE (*,*)
      WRITE (*,*) ' THE LAST LINE PROCESSED WAS:'
      WRITE (*,*) CHAR_X
      WRITE (*,*)
      CALL QUIT(16)

200   WRITE (*,*)
      WRITE (*,*) ' *** ERROR IN INPUT ***'
      WRITE (*,*) ' THE PROGRAM EXPECTED A REAL NUMBER'
      WRITE (*,*) ' BUT WAS SUPPLIED WITH AN EMPTY CHARACTER STRING'
      WRITE (*,*)
      WRITE (*,*) ' FORTRAN I/O STATUS = ',IOS
      WRITE (*,*)
      CALL QUIT(16)

      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine clscom
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine daclos
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine dafile
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine daname
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine datimx
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck decode_line $Revision$
      SUBROUTINE DECODE_LINE(LINE,NCHAR,NENTRY,IENTRY,MXENTRY)
*
* A CHAR line is given.
* Find number of separate items, with each item
* being separated by a ,
*
* Jeppe Olsen, Someday in 97 where I really should be doing more
* important things
*
*. Entry
      CHARACTER*(*) LINE
*. Output
      CHARACTER*72 IENTRY(MXENTRY)
*. Local scratch
      CHARACTER*72 CSCR
*
*
      JITEM=0
      JEFF = 0
      DO ICHAR = 0, NCHAR
        IF(ICHAR.EQ.0.OR.LINE(ICHAR:ICHAR).EQ.',') THEN
*Start of new item, make sure there is space and clean up
          JITEM = JITEM + 1
          IF(JITEM .GT.MXENTRY) THEN
            WRITE(6,*) 'DECODE_LINE:MXENTRY too small'
            WRITE(6,*) ' Number of entries larger than MXENTRY'
            WRITE(6,*) ' JITEM, MXENTRY', JITEM, MXENTRY
            Call Abend2( 'DECODE_LINE:MXENTRY too small' )
          END IF
*. Copy previous entry to permanent
          IF(JITEM.NE.1) THEN
            IENTRY(JITEM-1) = CSCR
          END IF
*. and clean
          DO JCHAR = 1, NCHAR
            CSCR(JCHAR:JCHAR) = ' '
          END DO
          JEFF = 0
        ELSE
*. Continuation of previous item
          JEFF = JEFF + 1
          CSCR(JEFF:JEFF) = LINE(ICHAR:ICHAR)
        END IF
*
      END DO
*. Transfer last item to permanant residence
      IF(JEFF.NE.0) IENTRY(JITEM) = CSCR
*
      NENTRY = JITEM
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from Decode line '
        WRITE(6,*) ' ========================'
        WRITE(6,*)
        WRITE(6,*) ' Number of separate entries', NENTRY
        WRITE(6,*)
        DO JENTRY = 1, NENTRY
          WRITE(6,'(A,I3,72A)') 'Entry ',JENTRY,IENTRY(JENTRY)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine errmsg
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck evenodd $Revision$
*  This (stupid) little subroutine merely checks whether an
*  integer is odd or even and passes this info to the calling
*  routine.
*
      subroutine evenodd(IEVOD,INTTEST)
      implicit real*8 (A-H,O-Z)
* add 1
      INTTESTP = INTTEST + 1
* put to real value
      A=INTTEST
      B=INTTESTP
* divide both values by 2
      FIRST = A/2
      SECOND = B/2
* cut off places right of point
      ICUT1 = FIRST
      ICUT2 = SECOND
* compare results
      if (ICUT1.eq.ICUT2) then
* INTTEST is even
         IEVOD = 2
      else
* INTTEST is odd
         IEVOD = 1
      end if
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine findhigh(IARRAY,LENGTH,MVALUE)
*
*  Find highest integer on IARRAY
*  Positive integers assumed
*
      implicit real*8(A-H,O-Z)
*
      dimension IARRAY(LENGTH)
*
      MVALUE = 0
      do I=1,LENGTH,1
         if (IARRAY(I).gt.MVALUE) MVALUE = IARRAY(I)
      end do
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine finish(IDUMMY)
      implicit none
      integer IDUMMY
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck fndmnd $Revision$
      SUBROUTINE FNDMND(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,
     &SCR,ISCAT,SUBVAL,NTESTG)
*
* FIND NSUB LOWEST ELEMENTS OF VECTOR RESIDING ON FILE
* LU. ENSURE THAT NO DEGENERENCIES ARE SPLIT
*
*
* INPUT
*======
* LU : FILE WHERE VECTOR OF INTEREST IS RESIDING, REWINDED
* LBLK : DEFINES FILE STRUCTURE ON FILE LU
* NSUBMX: LARGEST ALLOWED NUMBER OF SORTED ELEMENTS
*
*OUTPUT
*======
* NSUB : ACTUAL NUMBER OF ELEMENTS OBTAINED. CAN BE SMALLER
*       THAN NSUBMX IF THE LATST ELEMENT BELONGS TO A DEGENERATED
*       SET
*ISCAT: SCATTERING ARRAY , ISCAT(I) GIVED FULL ADRESS OF SORTED
*       ELEMENT I
*       SUBVAL : VALUE OF SORTED ELEMENTS

      IMPLICIT REAL*8           ( A-H,O-Z)
      DIMENSION SEGMNT(*), ISCAT(*),SUBVAL(*),SCR(*),ISCR(*)
*
      CALL REWINE(LU,-1)
*
*.LOOP OVER BLOCKS
*
C     write(6,*) ' FNDMND NSUBMX = ', NSUBMX
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
      IBLK = 0
      IBASE = 1
      LSUB = 0
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          CALL IFRMDS(LBL,1,-1,LU)
        END IF
        IBLK = IBLK + 1
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Info about block ',IBLK
          WRITE(6,*) ' Number of elements ', LBL
        END IF
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(SEGMNT,LBL ,-1,LU,IMZERO,IAMPACK)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Elements read in '
             CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
           END IF
           IF(LBL .GE. 0 ) THEN
*. LOWEST ELEMENTS IN SEGMNT  ( ADD TO PREVIOUS LIST )
             MSUBMX = MIN(NSUBMX,LBL)
C     SUBROUTINE SORLOW(WRK,STVAL,ISTART,KZVAR,KEXSTV,JEXSTV,IPRT)
             IF(LBL.GE.1) THEN
               CALL SORLOW(SEGMNT,SCR(1+LSUB),ISCR(1+LSUB),LBL,
     &                     MSUBMX,MSUB,NTEST)
C              WRITE(6,*)
C    &       ' After SORLOW 1 , Scatter array for combined list '
C            call iwrtma(ISCR(1),1,LSUB+MSUB,1,LSUP+MSUP)
             ELSE
              MSUB = 0
             END IF
             DO 10 I = 1, MSUB
   10        ISCR(LSUB+I) = ISCR(LSUB+I) + IBASE - 1
C            Write(6,*)
C    &       ' After 10 , Scatter array for combined list '
C            call iwrtma(ISCR,1,LSUB+MSUB,1,LSUP+MSUP)
* SORT COMBINED LIST
             MSUBMX = MIN(NSUBMX,LSUB+MSUB)
             IF(MSUBMX.GT.0) THEN
             CALL SORLOW(SCR,SUBVAL,ISCAT,LSUB+MSUB,MSUBMX,LLSUB,
     &         NTEST)
C              WRITE(6,*)
C    &       ' After SORLOW 2 , Scatter array for combined list '
C            call iwrtma(ISCR(1),1,LSUB+MSUB,1,LSUP+MSUP)
             ELSE
               LLSUB = 0
             END IF
             LSUB = LLSUB
             DO 20 I = 1, LSUB
               ISCR(I+2*NSUBMX) = ISCR(ISCAT(I))
   20        CONTINUE
*
             CALL ICOPVE(ISCR(1+2*NSUBMX),ISCR(1),LSUB)
             CALL COPVEC(SUBVAL,SCR,LSUB)
             IBASE = IBASE + LBL
             IF(NTEST .GE. 20 ) THEN
               WRITE(6,*) ' Lowest elements and their original place '
               WRITE(6,*) ' Number of elements obtained ', LSUB
               CALL WRTMAT(SUBVAL,1,LSUB,1,LSUB)
               CALL IWRTMA(ISCR,1,LSUB,1,LSUB)
             END IF
          END IF
*
        END IF
*
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      NSUB = LSUB
      CALL ICOPVE(ISCR,ISCAT,NSUB)
      IF(NTEST .GE. 5) THEN
        WRITE(6,*) ' Lowest elements and their original place '
        WRITE(6,*) ' Number of elements obtained ', NSUB
       CALL WRTMAT(SUBVAL,1,NSUB,1,NSUB)
       CALL IWRTMA(ISCAT,1,NSUB,1,NSUB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Fold square matrix to lower triangle form
*  The output matrix is a one-dim. field with
*  N(N+1)/2 elements
*
      subroutine fold(VECTOR,ARRAY,IDIM,ISYM,NTRIAEL)
      implicit real*8 (A-H,O-Z)
      dimension ARRAY(IDIM,IDIM),VECTOR(NTRIAEL)
*
      NTEST = 0
      do M=1,IDIM,1
         J = 0
         do L=1,M,1
            J = J + L
         end do
         ICOUNT = J
         do N=M,IDIM,1
            VECTOR(ICOUNT) = ARRAY(N,M)
            if (N.eq.IDIM.and.M.eq.IDIM) goto 100
            ICOUNT = ICOUNT + N
         end do
      end do
100   continue
      if (NTEST.ne.0) then
        if (ICOUNT.eq.NTRIAEL) then
          write(6,*) 'I have folded successfully.'
        else
          write(6,*) 'Fold unsuccessful.'
          call abend('Quitting.')
        end if
        write(6,*) 'Vector is:'
        do I=1,NTRIAEL,1
           write(6,*) VECTOR(I)
        end do
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck gflush $Revision$
      SUBROUTINE GFLUSH(LU)
c     CALL xFLUSH(LU)
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Save an integer to a real array
*
      subroutine i_save_r(WORK,IVALUE,IELMNT)
*
      implicit none
      integer WORK(*),IVALUE,IELMNT
*
      WORK(IELMNT) = IVALUE
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck icopmt $Revision$
      SUBROUTINE ICOPMT(MATI,NRI,NCI,MATO,NRO,NCO)
*
* Copy integer matrix MATI to MATO
*
*. Input
      INTEGER  MATI(NRI,NCI)
*. Output
      INTEGER MATO(NRO,NCO)
*
      NREFF = MIN(NRI,NRO)
      NCEFF = MIN(NCI,NCO)
*
      DO IC = 1, NCEFF
        DO IR = 1, NREFF
          MATO(IR,IC) = MATI(IR,IC)
        END DO
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck icopve2 $Revision$
      SUBROUTINE ICOPVE2(IIN,IOFF,NDIM,IOUT)
*
* IOUT(I) = IIN(IOFF-1+I),I = 1, NDIM
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION IIN(*)
*. Output
      DIMENSION IOUT(*)
*
      DO I = 1, NDIM
        IOUT(I) = IIN(IOFF-1+I)
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck icpmt2 $Revision$
      SUBROUTINE ICPMT2(AIN,AOUT,NINR,NINC,NOUTR,NOUTC,IZERO)
*
* Copy INTEGER matrix AIN to AOUT . Dimensions can differ
*
* If IZERO .ne. 0 , AOUT is zeroed  first
      IMPLICIT REAL*8           (A-H,O-Z)
*. Input
      INTEGER AIN(NINR,NINC)
*. Output
      INTEGER AOUT(NOUTR,NOUTC)
*
      IF(IZERO.NE.0) CALL ISETVC(AOUT,0,NOUTR*NOUTC)
      DO 100 J = 1, NINC
       CALL ICOPVE(AIN(1,J),AOUT(1,J),NINR)
  100 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck ielsum $Revision$
      FUNCTION IELSUM(IVEC,NELMNT)
*
* Sum elements of integer vector IVEC
*
      DIMENSION IVEC(*)
*
      ISUM = 0
      DO 100 IELMNT = 1, NELMNT
        ISUM = ISUM + IVEC(IELMNT)
  100 CONTINUE
*
      IELSUM = ISUM
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck ifreq $Revision$
      FUNCTION IFREQ(IVEC,IVAL,NDIM)
*
* Number of times IVAL occurs in IVEC
*
      DIMENSION IVEC(*)
*
      NTIME = 0
      DO 100 I = 1, NDIM
        IF(IVEC(I).EQ.IVAL) NTIME = NTIME + 1
  100 CONTINUE
*
      IFREQ = NTIME
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck ifrmr $Revision$
      FUNCTION IFRMR(WORK,IROFF,IELMNT)
*
* An integer array is stored in real array WORK,
* starting from WORK(IROFF). Obtain element
* IELMNT of this array
*
      INTEGER WORK(*)
*
#include "irat.inc"
*. offset when work is integer array
      IIOFF = 1 + IRAT * (IROFF-1)
      IFRMR = WORK(IIOFF-1+IELMNT)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck imnmx $Revision$
      FUNCTION IMNMX(IVEC,NDIM,MINMAX)
*
*     Find smallest (MINMAX=1) or largest (MINMAX=2)
*     absolute value of elements in integer vector IVEC
*
      DIMENSION IVEC(1)
*
      IF(NDIM.GT.0) THEN
        IX = -1
        IF(MINMAX.EQ.1) THEN
          IX=ABS(IVEC(1))
          DO I=2,NDIM
            IX=MIN(IX,ABS(IVEC(I)))
          END DO
        END IF
*
        IF(MINMAX.EQ.2) THEN
          IX=ABS(IVEC(1))
          DO I=2,NDIM
            IX=MAX(IX,ABS(IVEC(I)))
          END DO
        END IF
*
        IMNMX = IX
*
      ELSE IF(NDIM.EQ.0) THEN
*. No components : set to zero and write a warning
        IMNMX = 0
        WRITE(6,*) ' Min/Max taken zero length vector set to zero'
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck inprodb $Revision$
      REAL*8 FUNCTION INPRODB(VEC1,VEC2,NBLK,LBLK,I0BLK)
*
* Inner products between blocked vectors with check of
* zero blocks. Zero blocks are flagged by a unit entry
* in I0BLK or a negative block length
*
*Jeppe Olsen, May 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC1(*),VEC2(*), LBLK(NBLK),I0BLK(NBLK)
*
      REAL*8 INPROD
*
      X = 0.0D0
      IOFF = 1
      DO IBLK = 1, NBLK
        LLBLK = LBLK(IBLK)
C?      WRITE(6,*) ' INPRODB IBLK LLBLK IOFF ',IBLK,LLBLK,IOFF
        IF(I0BLK(IBLK).EQ.0.AND.LLBLK.GT.0) THEN
          X = X + INPROD(VEC1(IOFF),VEC2(IOFF),LLBLK)
C?        WRITE(6,*) ' Vec1 and Vec2 blocks '
C?        CALL WRTMAT(VEC1(IOFF),1,LLBLK,1,LLBLK)
C?        CALL WRTMAT(VEC2(IOFF),1,LLBLK,1,LLBLK)
C?        write(6,*) ' Updated x', X
        END IF
        IF(LLBLK.GT.0) IOFF = IOFF + LLBLK
      END DO
*
      INPRODB = X
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck invmat $Revision$
      SUBROUTINE INVMAT(A,B,MATDIM,NDIM)
C FIND INVERSE OF MATRIX A
C INPUT :
C        A : MATRIX TO BE INVERTED
C        B : SCRATCH ARRAY
C        MATDIM : PHYSICAL DIMENSION OF MATRICES
C        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
C
C OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
C WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MATDIM,MATDIM),B(MATDIM,MATDIM)
C
      ITEST=0
      IF(NDIM.EQ.1)THEN
        IF(A(1,1) .NE. 0.0D0 ) THEN
           A(1,1) = 1.0D0/A(1,1)
        ELSE
           ITEST = 1
        END IF
      ELSE
        DETERM=0.0D0
        EPSIL=0.0D0
        CALL BNDINV(A,B,NDIM,DETERM,EPSIL,ITEST,MATDIM)
      END IF
C
      IF( ITEST .NE. 0 ) THEN
        WRITE (6,'(A,I3)') ' INVERSION PROBLEM NUMBER..',ITEST
      END IF
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' INVERTED MATRIX '
        CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck ivcsum $Revision$
      SUBROUTINE IVCSUM(IA,IB,IC,IFACB,IFACC,NDIM)
*
* Add two (scaled) integer vectors
*
*        IA(*) = IFACB*IB(*) + IFACC*IC(*)
*
      DIMENSION IA(*),IB(*),IC(*)
*
      DO 100 I = 1, NDIM
        IA(I) = IFACB * IB(I) + IFACC * IC(I)
  100 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck iwrtma10 $Revision$
      SUBROUTINE IWRTMA10(IMAT,NROW,NCOL,MAXROW,MAXCOL)
* I10 format
      DIMENSION IMAT(MAXROW,MAXCOL)
C
      DO 100 I = 1, NROW
        WRITE(6,1110) (IMAT(I,J),J= 1,NCOL)
  100 CONTINUE
 1110 FORMAT(/,1X,8I10,/,(1X,8I10))
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine leftad
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck lftpos $Revision$

      SUBROUTINE LFTPOS(CARD,LENGTH)
*
* Eliminate blanks and left position character string CARD
*
      CHARACTER*72 CARD
*
      IEFF = 0
C     WRITE(6,*) ' INPUT string to LFTPOS '
C     WRITE(6,'(1H ,A)') CARD

      DO 100 IPOS = 1, LENGTH
       IF(CARD(IPOS:IPOS).NE.' ') THEN
         IEFF = IEFF + 1
         IF(IEFF.NE.IPOS) CARD(IEFF:IEFF) = CARD(IPOS:IPOS)
       END IF
  100 CONTINUE
*.Fill end with trailing blanks
      DO 200 IPOS = IEFF+1,LENGTH
        CARD(IPOS:IPOS) = ' '
  200 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Left adjusted character string '
        WRITE(6,'(1H ,A)') CARD
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck matml7 $Revision$


      SUBROUTINE MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C =  FACTORC*C + FACTORAB* A * B             FOR ITRNSP = 0
C
C     C =  FACTORC*C + FACTORAB* A(TRANSPOSED) * B FOR ITRNSP = 1
C
C     C =  FACTORC*C + FACTORAB* A * B(TRANSPOSED) FOR ITRNSP = 2
C
C... JEPPE OLSEN,
C
*. Notice : If the summation index has dimension zero nothing
*           is performed
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
      COMMON/MATMLST/XNFLOP,XNCALL,XLCROW,XLCCOL,XLCROWCOL,TMULT
*
c     CALL QENTER('MATML')
*
      NTEST = 00
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) ' A, B and C MATRIX FROM MATML7 '
      WRITE(6,*)
      CALL WRTMAT(A,NAROW,NACOL,NAROW,NACOL)
      CALL WRTMAT(B,NBROW,NBCOL,NBROW,NBCOL)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      WRITE(6,*)      ' NCROW NCCOL NAROW NACOL NBROW NBCOL '
      WRITE(6,'(6I6)')  NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL
      WRITE(6,*) ' ITRNSP : ', ITRNSP
      END IF
*. Statistics
      XNCALL = XNCALL + 1
      XLCROW = XLCROW + NCROW
      XLCCOL = XLCCOL + NCCOL
      XLCROWCOL = XLCROWCOL + NCROW*NCCOL
      T_INI = Second()
*
      IF(ITRNSP.EQ.1) THEN
        XNFLOP = XNFLOP + 2*NCROW*NCCOL*NAROW
      ELSE
        XNFLOP = XNFLOP + 2*NCROW*NCCOL*NACOL
      END IF
*
      IF(NAROW*NACOL*NBROW*NBCOL*NCROW*NCCOL .EQ. 0 ) THEN
        IZERO = 1
      ELSE
        IZERO = 0
      END IF

      IESSL = 1
      ICONVEX = 0
      IF ((ICONVEX.EQ.1 .OR.IESSL.EQ.1).AND. IZERO.EQ.0 ) THEN
*. DGEMM from CONVEX/ESSL  lib
        LDA = MAX(1,NAROW)
        LDB = MAX(1,NBROW)
*
        LDC = MAX(1,NCROW)
        IF(ITRNSP.EQ.0) THEN
        CALL DGEMM('N','N',NAROW,NBCOL,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF (ITRNSP.EQ.1) THEN
        CALL DGEMM('T','N',NACOL,NBCOL,NAROW,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF(ITRNSP.EQ.2) THEN
        CALL DGEMM('N','T',NAROW,NBROW,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        END IF
      ELSE
* Use Jeppes version ( it should be working )
        IF( ITRNSP .EQ. 0 ) THEN
* ======
* C=A*B
* ======
*
C         CALL SCALVE(C,FACTORC,NCROW*NCCOL)
          DO J =1, NCCOL
*. Initialize with FACTORC*C(I,J) + FACTORAB*A(I,1)*B(1,J)
            IF(NBROW.GE.1) THEN
              B1J = FACTORAB*B(1,J)
              DO I = 1, NCROW
                C(I,J) = FACTORC*C(I,J) + B1J*A(I,1)
              END DO
            END IF
*. and the major part
            DO K =2, NBROW
              BKJ = FACTORAB*B(K,J)
              DO I = 1, NCROW
                C(I,J) = C(I,J) + BKJ*A(I,K)
              END DO
            END DO
          END DO
*
        END IF
        IF ( ITRNSP .EQ. 1 ) THEN
*
* =========
* C=A(T)*B
* =========
*
          DO J = 1, NCCOL
            DO I = 1, NCROW
              T = 0.0D0
              DO K = 1, NBROW
                T = T  + A(K,I)*B(K,J)
              END DO
              C(I,J) = FACTORC*C(I,J) + FACTORAB*T
            END DO
          END DO
        END IF
C
        IF ( ITRNSP .EQ. 2 ) THEN
C ===========
C. C = A*B(T)
C ===========
          DO J = 1,NCCOL
*. Initialization
            IF(NBCOL.GE.1) THEN
              BJ1 = FACTORAB*B(J,1)
              DO I = 1, NCROW
                C(I,J) = FACTORC*C(I,J) + BJ1*A(I,1)
              END DO
            END IF
*. And the rest
            DO K = 2,NBCOL
              BJK = FACTORAB*B(J,K)
              DO I = 1, NCROW
                C(I,J) = C(I,J) + BJK*A(I,K)
              END DO
            END DO
          END DO
        END IF
      END IF
C
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) ' C MATRIX FROM MATML7 '
      WRITE(6,*)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
C
      T_END = Second()
      TMULT = TMULT + T_END - T_INI
C
C     CALL QEXIT('MATML')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck matvcb $Revision$
      SUBROUTINE MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION   MATRIX(MATDIM,MATDIM),VECIN(2),VECOUT(2)
C
C     VECOUT=MATRIX*VECIN FOR ITRNSP=0
C     VECOUT=MATRIX(TRANSPOSED)*VECIN FOR ITRNSP .NE. 0
C
      DO 10 I=1,NDIM
   10 VECOUT(I)=0.0D0
      IF(ITRNSP.EQ.0) THEN
C
       DO 100 J=1,NDIM
        VECINJ=VECIN(J)
        DO 90 I=1,NDIM
         VECOUT(I)=VECOUT(I)+MATRIX(I,J)*VECINJ
   90   CONTINUE
  100  CONTINUE
      END IF
C
      IF(ITRNSP.NE.0) THEN
       DO 200 I=1,NDIM
        X=0.0D0
        DO 190 J=1,NDIM
         X=X+MATRIX(J,I)*VECIN(J)
  190   CONTINUE
        VECOUT(I)=X
  200  CONTINUE
      END IF
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck matvcc $Revision$
      SUBROUTINE MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
*
* ITRNS = 0 : VOUT(I) = A(I,J)*VIN(J)
* ITRNS = 1 : VOUT(I) = A(J,I)*VIN(J)
*
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION A(NROW,NCOL)
       DIMENSION VIN(*),VOUT(*)
*
      IF(ITRNS.EQ.0) THEN
        DO 10 I = 1, NROW
   10   VOUT(I) = 0.0D0
        DO 100 J = 1, NCOL
         VINJ = VIN(J)
         DO 60 I = 1, NROW
           VOUT(I) = VOUT(I) + A(I,J)*VINJ
   60    CONTINUE
  100   CONTINUE
      ELSE IF( ITRNS.EQ.1) THEN
        DO 200 I = 1, NCOL
          X = 0.0D0
          DO 160 J = 1, NROW
            X = X + A(J,I)*VIN(J)
  160     CONTINUE
          VOUT(I) = X
  200   CONTINUE
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck nxtij $Revision$
* Output
      SUBROUTINE NXTIJ(I,J,NI,NJ,IJSM,NONEW)
*
* An ordered pair (I,J) is given ,I.LE.NI,J.LE.NJ
*
* Find next pair, if IJSM .ne. 0 ,I .ge. J
*
      NONEW = 0
  100 CONTINUE
      IF(I.LT.NI) THEN
        I = I + 1
      ELSE
        IF(J.LT.NJ) THEN
          I = 1
          J = J+1
        ELSE
          NONEW = 1
          GOTO 101
        END IF
      END IF
      IF(IJSM.NE.0.AND.I.LT.J) GOTO 100
  101 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' next (i,j) pair ', I,J
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck nxtnum2 $Revision$
      SUBROUTINE NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An set of numbers INUM(I),I=1,NELMNT is
* given. Find next compund number.
* Digit I must be in the range MINVAL,MAXVAL(I).
*
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen Oct 1994
*
*. Input
      DIMENSION MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTNUM '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IF(NELMNT.EQ.0) THEN
       NONEW = 1
       GOTO 1001
      END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
          DO JPLACE = 1, IPLACE
            INUM(JPLACE) = 1
          END DO
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck nxtnum2_rev $Revision$
      SUBROUTINE NXTNUM2_REV(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An set of numbers INUM(I),I=1,NELMNT is
* given. Find next compund number.
* Digit I must be in the range MINVAL,MAXVAL(I).
*
* In this version rightmost digits are increased first
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen March 1998
*
*. Input
      DIMENSION MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTNUM '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IF(NELMNT.EQ.0) THEN
       NONEW = 1
       GOTO 1001
      END IF
*
      IPLACE = NELMNT + 1
 1000 CONTINUE
        IPLACE = IPLACE - 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.GT.1 ) THEN
          INUM(IPLACE) = 1
        ELSE IF ( IPLACE. EQ. 1 ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck nxtnum3 $Revision$
      SUBROUTINE NXTNUM3(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An set of numbers INUM(I),I=1,NELMNT is
* given. Find next compund number.
* Digit I must be in the range MINVAL(I),MAXVAL(I).
*
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen Oct 1994
*
*. Input
      DIMENSION MINVAL(*),MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTNUM '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IF(NELMNT.EQ.0) THEN
        NONEW = 1
        GOTO 1001
      END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
          DO JPLACE = 1, IPLACE
            INUM(JPLACE) = MINVAL(JPLACE)
          END DO
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck nxtord $Revision$
      SUBROUTINE NXTORD(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An ordered set of numbers INUM(I),I=1,NELMNT is
* given in strictly ascending order. Values of INUM(*) is
* restricted to the interval MINVAL,MAXVAL .
*
* Find next higher number.
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen May 1989
*
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTORD '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF( IPLACE .LT. NELMNT .AND.
     &      INUM(IPLACE)+1 .LT. INUM(IPLACE+1)
     &  .OR.IPLACE.EQ. NELMNT .AND.
     &      INUM(IPLACE)+1.LE.MAXVAL) THEN
              INUM(IPLACE) = INUM(IPLACE) + 1
              NONEW = 0
              GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
              IF(IPLACE .EQ. 1 ) THEN
                INUM(IPLACE) = MINVAL
              ELSE
                INUM(IPLACE) = INUM(IPLACE-1) + 1
              END IF
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
              NONEW = 1
              GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine opncom
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck ordint $Revision$
      SUBROUTINE ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
*
* ORDER A STRING OF INTEGERS TO ASCENDING ORDER
*
* IINST : INPUT STRING
* IOUTST : OUTPUT STRING
* NELMNT : NUMBER OF INTEGERS
* INO : Mapping array from new to old order
*
* THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
* ( HE IS HEREBY AKNOWLEDGED , AND I AM EXCUSED )
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IINST(NELMNT),IOUTST(NELMNT),INO(NELMNT)
*
      IF(NELMNT.EQ.0) GOTO 1001
      CALL ICOPVE(IINST,IOUTST,NELMNT)
      DO  5 I = 1, NELMNT
        INO(I) = I
    5 CONTINUE
C
C       BEGIN TO ORDER
C
        JOE = 1
  10    I = JOE
  20    CONTINUE
        IF(I.EQ.NELMNT) GO TO 50
        IF(IOUTST(I).LE.IOUTST(I+1)) GO TO 40
        JOE = I + 1
  30    SWAP = IOUTST(I)
        IOUTST(I) = IOUTST(I+1)
        IOUTST(I+1) = SWAP
        ISWAP = INO(I)
        INO(I) = INO(I+1)
        INO(I+1) = ISWAP
        IF(I.EQ.1) GO TO 10
        I = I - 1
        IF(IOUTST(I).GT.IOUTST(I+1)) GO TO 30
        GO TO 10
 40     I = I + 1
      GO TO 20
C
C     END ORDER
C
 50   CONTINUE
*
 1001 CONTINUE
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GE.200) THEN
        WRITE(6,*) ' Result from ORDINT '
        WRITE(6,*)
        WRITE(6,*)  ' Input string '
        CALL IWRTMA(IINST,1,NELMNT,1,NELMNT)
        WRITE(6,*)  ' Ordered string '
        CALL IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
        WRITE(6,*) ' New to old order '
        CALL IWRTMA(INO,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE QSTAT
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine rdnlst
      read(5,*)
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine rightad
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck scalar_product_of_b $Revision$
      FUNCTION SCALAR_PRODUCT_OF_B
     &(VEC1,VEC2,NBLOCK1,IBLOCK1,I2EQ1,NBLOCK2,I2TO1,IOFF)
*
* Two blocked vectors VEC1 and VEC2 are given.
* VEC1 contains NBLOCK1 blocks defined by IBLOCK1.
* If I2EQ1.ne.0  VEC2 has the same block structure as VEC1
* IF I2EQ1 .eq.0 VEC2 contains NBLOCK2 blocks, and I2TO1 gives
* the mapping from blocks of VEC2 to blocks of vec1
*
* Find the scalar product between these two vectors
*
* Jeppe Olsen, October 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION VEC1(*),VEC2(*),IBLOCK1(8,*)
      REAL*8 INPROD
*. In use of I2EQ1 = 0
      DIMENSION I2TO1(*)
*
      X = 0.0D0
      IOFF2 = 1
      DO JBLOCK = 1, NBLOCK2
        IF(I2EQ1.NE.0) THEN
          JJBLOCK = JBLOCK-1+IOFF
          IOFF1   = IBLOCK1(6,JJBLOCK)
          IOFF2   = IOFF1
        ELSE
          JJBLOCK = I2TO1(JBLOCK)
          IOFF1 = IBLOCK1(6,JJBLOCK)
        END IF
        NELMNT = IBLOCK1(8,JJBLOCK)
        X = X + INPROD(VEC1(IOFF1),VEC2(IOFF2),NELMNT)
        IF(I2EQ1.EQ.0) IOFF2 = IOFF2 + NELMNT
      END DO
      SCALAR_PRODUCT_OF_B = X
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Output from INNER_PRODUCT_OF_BLOCKED_VECTORS'
         WRITE(6,*)
         WRITE(6,*) ' Number of blocks included ', NBLOCK2
         WRITE(6,*) ' Inner product ', X
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck scldia $Revision$
      SUBROUTINE SCLDIA(A,FACTOR,NDIM,IPACK)
*
* scale diagonal of square matrix A
*
* IPACK = 0 : full matrix
* IPACK .NE. 0 : Lower triangular packed matrix
*                assumed packed columnwise !!!!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*)
*
      IF( IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          II = (I-1)*NDIM + I
          A(II) = A(II) * FACTOR
  100   CONTINUE
      ELSE
        II = 1
        DO 200 I = 1, NDIM
          A(II) = A(II) * FACTOR
          II = II + NDIM - I + 1
  200   CONTINUE
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine sqfile
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck sqrtmt $Revision$
      SUBROUTINE SQRTMT(A,NDIM,ITASK,ASQRT,AMSQRT,SCR)
*
* Calculate square root of positive definite symmetric matrix A
* if(ITASK .EQ. 2 ) Inverted square root matrix is also calculated
* In case of singularities in A A -1/2 is defined to have the same
* singularity
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
*
      DIMENSION A(NDIM,NDIM)
      DIMENSION ASQRT(NDIM,NDIM),AMSQRT(NDIM,NDIM)
      DIMENSION SCR(*)
* Length of SCR should at least be 2 * NDIM ** 2 + NDIM*(NDIM+1)/2
      KLFREE = 1
*
      KLASYM = KLFREE
      KLAVAL = KLASYM
      KLFREE = KLASYM + NDIM*(NDIM+1)/2
*
      KLAVEC = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
      NTEST = 0
*
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(A,SCR(KLASYM),1,NDIM,NDIM)
      CALL EIGEN(SCR(KLASYM),SCR(KLAVEC),NDIM,0,1)
      CALL COPDIA(SCR(KLASYM),SCR(KLAVAL),NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of matrix : '
        CALL WRTMAT(SCR(KLAVAL),NDIM,1,NDIM,1)
      END IF
*. Check for negative eigenvalues
      DO I = 1, NDIM
       IF(SCR(KLAVAL-1+I).LT.0.0D0) THEN
         WRITE(6,*) ' SQRTMT : Negative eigenvalue ', SCR(KLAVAL-1+I)
         WRITE(6,*) ' SQRTMT : I will STOP '
         Call Abend2(       ' SQRTMT : Negative eigenvalue ' )
       END IF
      END DO
*
      DO 100 I = 1,NDIM
        SCR(KLAVAL-1+I) = SQRT(SCR(KLAVAL-1+I))
  100 CONTINUE
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(ASQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
*
      IF(ITASK .EQ. 2 ) THEN
        DO 200 I = 1,NDIM
          IF(SCR(KLAVAL-1+I) .GT. 1.0D-13 ) then
            SCR(KLAVAL-1+I) = 1.0D0/SCR(KLAVAL-1+I)
          ELSE
            SCR(KLAVAL-1+I) = SCR(KLAVAL-1+I)
          END IF
  200   CONTINUE
        CALL XDIAXT(AMSQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
      END IF
*
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Info from SQRTMT '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' Input matrix to SQRTMT '
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' Square root of matrix '
        CALL WRTMAT(ASQRT,NDIM,NDIM,NDIM,NDIM)
        IF(ITASK .EQ. 2 ) THEN
          WRITE(6,*) ' Inverse square root of matrix '
          CALL WRTMAT(AMSQRT,NDIM,NDIM,NDIM,NDIM)
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine start(CDUMMY)
      implicit none
      character*6 CDUMMY
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine stdfmt
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine stka
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine strnln
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck swapve $Revision$

      SUBROUTINE SWAPVE(VEC1,VEC2,NDIM)
C
C      SWAP ELEMENTS OF VECTORS VEC1 AND VEC2
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC1(1   ) ,VEC2(1   )
      DO 100 I=1,NDIM
       BUF=VEC1(I)
       VEC1(I)=VEC2(I)
       VEC2(I)=BUF
  100 CONTINUE
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck tripak $Revision$
      SUBROUTINE TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
C
C ( NOT A SIMPLIFIED VERSION OF TETRAPAK )
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C   IWAY = 2 : PACKED TO FULL FORM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
C
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO 100 I = 1,NDIM
          DO 50  J = 1, I
           APAK(IJ+J) = AUTPAK(J,I)
   50     CONTINUE
          IJ = IJ + I
  100   CONTINUE
      END IF
C
      IF( IWAY .EQ. 2 ) THEN
        IJ = 0
        DO 200 I = 1,NDIM
          DO 150  J = 1, I
           AUTPAK(I,J) = APAK(IJ+J)
           AUTPAK(J,I) = APAK(IJ+J)
  150     CONTINUE
          IJ = IJ + I
  200   CONTINUE
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPAK '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSYM(APAK,NDIM)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck tripak_blkm $Revision$
*
      SUBROUTINE TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
*
* Switch between packed and unpacked form of a blocked matrix
*
* IWAY = 1 => Unpacked to packed
* IWAY = 2 => Packed to unpacked
*
* Jeppe Olsen, February 1, 1998 (Moensted Kalkgrubber )
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input/output
      DIMENSION AUTPAK(*),APAK(*)
*. Input
      DIMENSION LBLOCK(NBLOCK)
*
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IOFFU = 1
          IOFFP = 1
        ELSE
          IOFFU = IOFFU + LBLOCK(IBLOCK-1)**2
          IOFFP = IOFFP + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
        END IF
        L = LBLOCK(IBLOCK)
        CALL TRIPAK(AUTPAK(IOFFU),APAK(IOFFP),IWAY,L,L)
C            TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from TRIPAK_BLKM'
        WRITE(6,*)
        WRITE(6,*) ' Unpacked matrix'
        WRITE(6,*) ' ==============='
        CALL APRBLM2(AUTPAK,LBLOCK,LBLOCK,NBLOCK,0)
C            APRBLM2(A,LROW,LCOL,NBLK,ISYM)
        WRITE(6,*) ' Packed matrix '
        WRITE(6,*) ' =============='
        CALL APRBLM2(  APAK,LBLOCK,LBLOCK,NBLOCK,1)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck tripk2 $Revision$
      SUBROUTINE TRIPK2(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
C
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C              LOWER HALF OF AUTPAK IS STORED IN APAK
C   IWAY = 2 : PACKED TO FULL FORM
C              APAK STORED IN LOWER HALF
C               SIGN * APAK TRANSPOSED IS STORED IN UPPPER PART
C.. NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
C
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO 100 J = 1,NDIM
          DO 50  I = J , NDIM
           APAK(IJ+I) = AUTPAK(I,J)
   50     CONTINUE
          IJ = IJ +NDIM-J
  100   CONTINUE
      END IF
C
      IF( IWAY .EQ. 2 ) THEN
        IJ = 0
        DO 200 J = 1,NDIM
          DO 150  I = J,NDIM
           AUTPAK(J,I) = SIGN*APAK(IJ+I)
           AUTPAK(I,J) = APAK(IJ+I)
  150     CONTINUE
          IJ = IJ + NDIM-J
  200   CONTINUE
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPK2 '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSM2(APAK,NDIM)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck tripk3 $Revision$
*
      SUBROUTINE TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
C
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C              LOWER HALF OF AUTPAK IS STORED IN APAK
C   IWAY = 2 : PACKED TO FULL FORM
C              APAK STORED IN LOWER HALF
C               SIGN * APAK TRANSPOSED IS STORED IN UPPPER PART
C.. NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS
*
* Some considerations on cache minimization used for IMET = 2 Loop
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
*
*. Packing : No problem with cache misses
*
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO J = 1,NDIM
          DO I = J , NDIM
            APAK(IJ+I) = AUTPAK(I,J)
          END DO
          IJ = IJ +NDIM-J
        END DO
      END IF
*
* Unpacking : cache misses can occur so two routes
*
      IF( IWAY .EQ. 2 ) THEN
*. Use blocked algorithm
      IMET = 2
C     IMET =2
      IF(IMET.EQ.1) THEN
*. No blocking
        IJ = 0
        DO J = 1,NDIM
          DO I = J,NDIM
           AUTPAK(J,I) = SIGN*APAK(IJ+I)
           AUTPAK(I,J) = APAK(IJ+I)
          END DO
          IJ = IJ + NDIM-J
        END DO
      ELSE IF (IMET .EQ. 2 ) THEN
*. Blocking
C       LBLK = 40
        LBLK = 40
        NBLK = MATDIM/LBLK
        IF(LBLK*NBLK.LT.MATDIM) NBLK = NBLK + 1
        DO JBLK = 1, NBLK
          IF(JBLK.EQ.1) THEN
            JOFF = 1
          ELSE
            JOFF = JOFF + LBLK
          END IF
          JEND = MIN(JOFF+LBLK-1,MATDIM)
          DO IBLK = JBLK, NBLK
            IF(IBLK.EQ.JBLK) THEN
              IOFF = JOFF
            ELSE
              IOFF = IOFF + LBLK
            END IF
            IEND = MIN(IOFF+LBLK-1,MATDIM)
              DO J = JOFF,JEND
                IF(IBLK.EQ.JBLK) THEN
                  IOFF2 = J
                ELSE
                  IOFF2 = IOFF
                END IF
                IJOFF = (J-1)*MATDIM-J*(J-1)/2
                DO I = IOFF2,IEND
                  AUTPAK(J,I) = SIGN*APAK(IJOFF+I)
                  AUTPAK(I,J) = APAK(IJOFF+I)
                END DO
              END DO
*. End of loop over I and J
            END DO
          END DO
*. End of loop over blocks of I and J
        END IF
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPK3 '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSM2(APAK,NDIM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck trnmat $Revision$
      SUBROUTINE TRNMAT(A,X,SCRA,NDIM,MATDIM)
C
C TRANSFORM MATRIX A : X(TRANS)*A*X
C A IS OVERWRITTREN BY TRANSFORMED A
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MATDIM,1),X(MATDIM,1),
     +          SCRA(MATDIM,1)
C
C
      NTEST=0
C
      IF(NTEST.GE.3) THEN
       WRITE(6,1020)
 1020  FORMAT(1H0,'*** OUTPUT FROM TRANMAT')
       WRITE(6,1030)
 1030  FORMAT(1H0,'A- AND X-MATRIX')
       CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
       CALL WRTMAT(X,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
C A*X
      DO 1000 I=1,NDIM
       DO 900 J=1,NDIM
       AX=0.0D0
        DO 800 K=1,NDIM
         AX=AX+A(I,K)*X(K,J)
  800   CONTINUE
        SCRA(I,J)=AX
  900  CONTINUE
 1000 CONTINUE
C
      IF(NTEST.GE.2) THEN
       WRITE(6,1040)
 1040  FORMAT(1H0,' AX MATRIX')
       CALL WRTMAT(SCRA,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
C X(TRANS)*(A*X)
      DO 600 I=1,NDIM
       DO 500 J=1,NDIM
        XAX=0.0D0
        DO 400 K=1,NDIM
         XAX=XAX+X(K,I)*SCRA(K,J)
  400   CONTINUE
        A(I,J)=XAX
  500  CONTINUE
  600 CONTINUE
C
      IF(NTEST.GE.2) THEN
       WRITE(6,1010)
 1010  FORMAT(1H0,' TRANSFORMED MATRIX')
       CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck trpad3 $Revision$
      SUBROUTINE TRPAD3(MAT,FACTOR,NDIM)
*
*  MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
*
*. With some considerations of effective cache use for large
*  matrices
*
      IMPLICIT REAL*8           (A-H,O-Z)
      REAL*8            MAT(NDIM,NDIM)
              FAC2 = 1.0D0 - FACTOR**2
*
C     IWAY = 1
      IWAY = 2
      IF(IWAY.EQ.1) THEN
*
*. No blocking
*
*. Lower half
        DO  J = 1, NDIM
          DO  I = J, NDIM
            MAT(I,J) =MAT(I,J) + FACTOR * MAT(J,I)
          END DO
        END DO
*. Upper half
        IF( ABS(FACTOR) .NE. 1.0D0 ) THEN
          FAC2 = 1.0D0 - FACTOR**2
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) = FACTOR*MAT(I,J ) + FAC2 * MAT(J,I)
            END DO
          END DO
        ELSE IF(FACTOR .EQ. 1.0D0) THEN
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) = MAT(I,J )
            END DO
          END DO
        ELSE IF(FACTOR .EQ. -1.0D0) THEN
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) =-MAT(I,J )
            END DO
          END DO
        END IF
      ELSE IF(IWAY .EQ. 2 ) THEN
*. Simple blocking of matrix
        LBLK = 40
        NBLK = NDIM/LBLK
        IF(NBLK*LBLK.LT.NDIM) NBLK = NBLK + 1
        IOFF = 1-LBLK
C?      write(6,*) 'NBLK ',nblk
        DO IBLK = 1, NBLK
          IF(IBLK.EQ.-1) write(6,*) 'IBLK = ',IBLK
          IOFF = IOFF + LBLK
          IEND = MIN(IOFF+LBLK-1,NDIM)
          JOFF = 1 - LBLK
          DO JBLK = 1, IBLK
            JOFF = JOFF + LBLK
            JEND = MIN(JOFF+LBLK-1,NDIM)
*. Lower half
            DO  I = IOFF,IEND
              IF(IBLK.EQ.JBLK) JEND = I
              DO J = JOFF,JEND
                MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
              END DO
            END DO
*. Upper half
            IF( ABS(FACTOR) .NE. 1.0D0 ) THEN
              FAC2 = 1.0D0 - FACTOR**2
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I
                DO J = JOFF, JEND
                  MAT(J,I) = FACTOR*MAT(I,J ) + FAC2 * MAT(J,I)
                 END DO
               END DO
            ELSE IF(FACTOR .EQ. 1.0D0) THEN
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I -1
                DO J = JOFF, JEND
                  MAT(J,I) = MAT(I,J )
                END DO
              END DO
            ELSE IF(FACTOR .EQ. -1.0D0) THEN
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I
                DO J = JOFF, JEND
                  MAT(J,I) = - MAT(I,J )
                END DO
              END DO
            END IF
*. ENd of loop over blocks
          END DO
        END DO
*. End of IWAY branching
      END IF
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck trpmt3 $Revision$
      SUBROUTINE TRPMT3(XIN,NROW,NCOL,XOUT)
*
* XOUT(I,J) = XIN(J,I)
*
*. With a few considerations for large scale cases with cache minimization
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XIN(NROW,NCOL),XOUT(NCOL,NROW)
*
      IWAY = 2
      IF(IWAY.EQ.1) THEN
*. Straightforward, no blocking
        DO IROW =1, NROW
          DO ICOL = 1, NCOL
            XOUT(ICOL,IROW) = XIN(IROW,ICOL)
          END DO
        END DO
      ELSE IF(IWAY.EQ.2) THEN
*. Simple blocking of matrix
        LRBLK = 40
        LCBLK = 40
        NRBLK = NROW/LRBLK
        NCBLK = NCOL/LCBLK
        IF(LRBLK*NRBLK.NE.NROW) NRBLK = NRBLK + 1
        IF(LCBLK*NCBLK.NE.NCOL) NCBLK = NCBLK + 1
*
        DO IRBLK = 1,NRBLK
          IF(IRBLK.EQ.1) THEN
            IROFF = 1
          ELSE
            IROFF = IROFF + LRBLK
          END IF
          IREND = MIN(NROW,IROFF+LRBLK-1)
          DO ICBLK = 1, NCBLK
            IF(ICBLK.EQ.1) THEN
              ICOFF = 1
            ELSE
              ICOFF = ICOFF + LCBLK
            END IF
            ICEND = MIN(NCOL,ICOFF+LCBLK-1)
*
            DO IROW = IROFF,IREND
              DO ICOL = ICOFF,ICEND
                XOUT(ICOL,IROW) = XIN(IROW,ICOL)
              END DO
            END DO
*
          END DO
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck unpcpc $Revision$
       SUBROUTINE UNPCPC(IRC,IR,NR,IC)
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
* A compound index IRC is given for a given element in
* a column packed lower half of a matrix with NR rows.
*
* Find the corrsponding row and column numbers, IR, IC
*
* The relation to be fulfilled is
* IRC = (IC-1)*NR+IR -IC*(IC-1)/2
*
*. Start by solving equation assuming IR = IC
* the corresponding second equation  is
* -XC**2/2 + XC(NR+1.5) -(NR +IRC)
*
* Well according to Mister Bechmann
      A = -0.5D0
      B = FLOAT(NR) + 1.5D0
      C = -FLOAT(NR + IRC)
      XC = -B/(2.0D0*A)
     &   + 1.0D0/(2.0D0*A)* SQRT( B ** 2 - 4.0D0*A*C)
*. Round down to get column number
      IC = XC
*. Row number
      IR = IRC - (IC-1)*NR + IC*(IC-1)/2
*. Check
      IF(IC.GT.IR  .OR. IR.GT.NR .OR. IR .LE. 0 ) THEN
        WRITE(6,*) ' Dear Sir '
        WRITE(6,*) ' I am a subroutine called UNPCPC '
        WRITE(6,*) ' trying to do my best, but my programmer '
        WRITE(6,*) ' ( Jeppe, if you know him ) '
        WRITE(6,*) ' has left me with some problems '
        WRITE(6,*) ' So I quit in a microsecond or two '
        WRITE(6,*)
        WRITE(6,*) ' IRC IR IC NR ', IRC,IR,IC,NR
        Call Abend2( 'UNPCPC ' )
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine upcase
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck uppcas $Revision$
      SUBROUTINE UPPCAS(LINE,LENGTH)
*
* Convert letters in character string LINE to upper case
*
* very stupid and not vectorized !
*
      CHARACTER*(*) LINE
      PARAMETER (NCHAR = 41)
      CHARACTER*1 LOWER(NCHAR)
      CHARACTER*1 UPPER(NCHAR)
*
      DATA LOWER/'a','b','c','d','e',
     &           'f','g','h','i','j',
     &           'k','l','m','n','o',
     &           'p','q','r','s','t',
     &           'u','v','w','x','y',
     &           'z','+','-','<','>',
     &           '=','0','1','2','3',
     &           '4','5','6','7','8',
     &           '9'/
      DATA UPPER/'A','B','C','D','E',
     &           'F','G','H','I','J',
     &           'K','L','M','N','O',
     &           'P','Q','R','S','T',
     &           'U','V','W','X','Y',
     &           'Z','+','-','<','>',
     &           '=','0','1','2','3',
     &           '4','5','6','7','8',
     &           '9'/
*
      DO 100 ICHA = 1, LENGTH
        DO 50 I = 1,NCHAR
          IF(LINE(ICHA:ICHA).EQ.LOWER(I))
     &    LINE(ICHA:ICHA) = UPPER(I)
   50   CONTINUE
  100 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
* Copy vector to defined partition of matrix
* (Simply the way the vector is arranged in memory)
*
      subroutine lvectomat(VEC,WMAT,ILEN,NROW,NCOL,
     &                     NPARTROW,NPARTCOL,IROWOFF,ICOLOFF)
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
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck writve $Revision$
      SUBROUTINE WRITVE(VEC,NDIM)
      DOUBLE PRECISION  VEC
      DIMENSION VEC(1   )
C
      WRITE(6,1010) (VEC(I),I=1,NDIM)
 1010 FORMAT(1H0,2X,4(2X,E15.8),/,(1H ,2X,4(2X,E15.8)))
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck wrtmat $Revision$
      SUBROUTINE WRTMAT(A,NROW,NCOL,NMROW,NMCOL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
C
      DO 100 I=1,NROW
      WRITE(6,1010) I,(A(I,J),J=1,NCOL)
 1010 FORMAT(1H0,I3,2X,4(1X,E16.8),/,(1H ,5X,4(1X,E16.8)))
  100 CONTINUE
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck xquit $Revision$
      Subroutine xquit(rc)
      Integer rc

      Open(66,File='returncode',Form='Formatted',Status='unknown')
      Write(66,*) rc
      Close(66)
      Stop
      End
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck xuflow $Revision$
      SUBROUTINE XUFLOW
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck zero_offdiag $Revision$
      SUBROUTINE ZERO_OFFDIAG(A,NDIM,IPACK)
*
* Zero off-diagonal elements in matrix A
*
* Jeppe Olsen, Oct 1997
*
* IF IPACK.NE.0, the matrix is assumed packed (inner loop over columns)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*)
*
      IF(IPACK.EQ.0) THEN
        DO I = 1, NDIM
          DO J = 1, NDIM
            IF(I.NE.J) THEN
              IJ = (J-1)*NDIM+I
              A(IJ) = 0.0D0
            END IF
          END DO
        END DO
      ELSE
        DO I = 1, NDIM
          DO J = 1, I-1
            IJ = I*(I-1)/2 + J
            A(IJ) = 0.0D0
          END DO
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck zero_offdiag_blm $Revision$
* Version of Oct 21
      SUBROUTINE ZERO_OFFDIAG_BLM(A,NBLOCK,LBLOCK,IPACK)
*
* Zero off-diagonal elements in a block-diagonal matrix
* if IPACK .ne. 0, the matrix is packed (inner loop over columns)
*
* Jeppe Olsen, October 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(*),LBLOCK(NBLOCK)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFF = 1
       ELSE
         LPREV = LBLOCK(IBLOCK-1)
         IF(IPACK.EQ.0) THEN
           IOFF = IOFF + LPREV*(LPREV+1)/2
         ELSE
           IOFF = IOFF + LPREV*LPREV
         END IF
       END IF
       L = LBLOCK(IBLOCK)
       CALL ZERO_OFFDIAG(A(IOFF),L,IPACK)
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output matrix from ZERO_OFFDIAG_BLM'
        CALL APRBLM2(A,LBLOCK,LBLOCK,NBLOCK,IPACK)
C       CALL APRBLM2(WORK(KFOCK),NTOOBS,NTOOBS,NSMOB,ISM)
      END IF
*
      RETURN
      END
