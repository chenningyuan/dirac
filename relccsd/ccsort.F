C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                 OFF,OFF1,OFF2,RBUF1,RBUF2,CBUF1,CBUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I>J,K>L:KLREP) to array BUF2(IK,JL:JLREP=1)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 CBUF1(*),CBUF2(*)
      REAL*8 RBUF1(*),RBUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(2),A0,0,RBUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKJL0 = OFF(JLREP) 
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           ILJK0 = OFF(JKREP) 
     &          + (OFF2(JREP,KREP)+(K-1)*JFIE(JREP)) * NPAIR(ILREP)
           ILJK0 = ILJK0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           JKIL0 = OFF(ILREP) 
     &          + (OFF2(IREP,LREP)+(L-1)*IFIE(IREP)) * NPAIR(JKREP)
           JKIL0 = JKIL0 + OFF1(JREP,KREP) + (K-1)*JFIE(JREP) + 1
           JLIK0 = OFF(IKREP) 
     &          + (OFF2(IREP,KREP)+(K-1)*IFIE(IREP)) * NPAIR(JLREP)
           JLIK0 = JLIK0 + OFF1(JREP,LREP) + (L-1)*JFIE(JREP) + 1
            DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             DO I = IMIN, IFIE(IREP)
              IJKL = IJKL + 1
              IKJL = IKJL0 + I
              ILJK = ILJK0 + I
              JKIL = JKIL0 + (I-1)*NPAIR(JKREP)
              JLIK = JLIK0 + (I-1)*NPAIR(JLREP)
              IF (CARITH) THEN
                 IF (JLREP.EQ.1) CBUF2(IKJL) = CBUF1(IJKL)
                 IF (JKREP.EQ.1) CBUF2(ILJK) = - CBUF1(IJKL)
                 IF (ILREP.EQ.1) CBUF2(JKIL) = - CBUF1(IJKL)
                 IF (IKREP.EQ.1) CBUF2(JLIK) = CBUF1(IJKL)
              ELSE
                 IF (JLREP.EQ.1) RBUF2(IKJL) = RBUF1(IJKL)
                 IF (JKREP.EQ.1) RBUF2(ILJK) = - RBUF1(IJKL)
                 IF (ILREP.EQ.1) RBUF2(JKIL) = - RBUF1(IJKL)
                 IF (IKREP.EQ.1) RBUF2(JLIK) = RBUF1(IJKL)
              ENDIF
            ENDDO
            IKJL0 = IKJL0 + NPAIR(IKREP)
            ILJK0 = ILJK0 + NPAIR(ILREP)
            JKIL0 = JKIL0 + 1
            JLIK0 = JLIK0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT2 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                 OFF,OFF1,OFF2,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IK,JL:JLREP=1)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(2),A0,0,BUF2,1)
      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKJL0 = OFF(JLREP) 
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           IKJL0 = IKJL0 * RCW + 1
           ILJK0 = OFF(JKREP) 
     &          + (OFF2(JREP,KREP)+(K-1)*JFIE(JREP)) * NPAIR(ILREP)
           ILJK0 = ILJK0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           ILJK0 = ILJK0 * RCW + 1
           DO J = 1, JFIE(JREP)
            IF (JLREP.EQ.1) 
     &         CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKJL0),1)
            IF (JKREP.EQ.1) THEN 
               CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(ILJK0),1)
               CALL XSCAL (IFIE(IREP),-A1,BUF2(ILJK0),1)
            ENDIF
            IJKL = IJKL + IFIE(IREP) * RCW
            IKJL0 = IKJL0 + NPAIR(IKREP) * RCW
            ILJK0 = ILJK0 + NPAIR(ILREP) * RCW
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT3 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                 OFF,OFF1,OFF2,RBUF1,RBUF2,CBUF1,CBUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IK,JL:JLREP=1)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 CBUF1(*),CBUF2(*)
      REAL*8 RBUF1(*),RBUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(2),A0,0,RBUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           IKJL0 = OFF(JLREP) 
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
            DO J = 1, JFIE(JREP)
             DO I = 1, IFIE(IREP)
              IJKL = IJKL + 1
              IKJL = IKJL0 + I
              IF (CARITH) THEN
                 IF (JLREP.EQ.1) CBUF2(IKJL) = CBUF1(IJKL)
              ELSE
                 IF (JLREP.EQ.1) RBUF2(IKJL) = RBUF1(IJKL)
              ENDIF
            ENDDO
            IKJL0 = IKJL0 + NPAIR(IKREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT4 (NREP,MULTB,NPAIR,KFIE,LFIE,NTRIPL,OFF,OFF1,
     &                 BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IJ K,L:LREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       N = NPAIR(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        LOFF = (OFF(LREP) + OFF1(IJREP,KREP)) * RCW
        KOFF = (OFF(KREP) + OFF1(IJREP,LREP)) * RCW
        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
            IJK  = (K-1) * N * RCW + 1
            IJKL1 = LOFF + (L-1) * NTRIPL(LREP) * RCW + IJK 
            IJLK1 = KOFF + (K-1) * NTRIPL(KREP) * RCW + IJL 
            CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJKL1),1)
            CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJLK1),1)
            CALL XSCAL (N,-A1,BUF2(IJLK1),1)
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT5 (NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,OFF,OFF2,
     &                 BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I>J,KL:KLREP) to array BUF2(I,J KL:JKLREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF (IREP.LT.JREP) GOTO 10
        JKLREP = IREP
        IKLREP = JREP
        JKLOFF = (OFF(JKLREP)+OFF2(JREP,KLREP)*IFIE(IREP)) * RCW
        IKLOFF = (OFF(IKLREP)+OFF2(IREP,KLREP)*JFIE(JREP)) * RCW
        M2 = IFIE(IREP) * JFIE(JREP)
        DO J = 1, JFIE(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, IFIE(IREP)
            IJKL1 = JKLOFF + ((J-1) * IFIE(IREP) + I - 1) * RCW + 1
            JIKL1 = IKLOFF + ((I-1) * JFIE(JREP) + J - 1) * RCW + 1
            CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(IJKL1),M2)
            CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(JIKL1),M2)
            CALL XSCAL(N,-A1,BUF2(JIKL1),M2)
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT6 (NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,NTRIPL,
     &                 OFF,OFF1,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I,J;KL:KLREP) to array BUF2(I KL, J:JREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        JOFF = (OFF(JREP) + OFF1(IREP,KLREP)) * RCW
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            IKLJ1 = JOFF + ((J-1)*NTRIPL(JREP)+I-1) * RCW + 1
            CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(IKLJ1),IFIE(IREP))
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT7 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                 OFF,OFF1,OFF2,RBUF1,RBUF2,CBUF1,CBUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IL,KJ:KJREP=1)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 CBUF1(*),CBUF2(*)
      REAL*8 RBUF1(*),RBUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(2),A0,0,RBUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           KJREP = MULTB(KREP,JREP,2)
           ILREP = MULTB(IREP,LREP,2)
           ILKJ0 = OFF(KJREP) 
     &          + (OFF2(KREP,JREP)+K-1)*NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
            DO J = 1, JFIE(JREP)
             DO I = 1, IFIE(IREP)
              IJKL = IJKL + 1
              ILKJ = ILKJ0 + I
              IF (CARITH) THEN
                 IF (KJREP.EQ.1) CBUF2(ILKJ) = CBUF1(IJKL)
              ELSE
                 IF (KJREP.EQ.1) RBUF2(ILKJ) = RBUF1(IJKL)
              ENDIF
            ENDDO
            ILKJ0 = ILKJ0 + NPAIR(ILREP) * KFIE(KREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT8 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                 OFF,OFF1,OFF2,RBUF1,RBUF2,CBUF1,CBUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IL,KJ:KJREP=1)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 CBUF1(*),CBUF2(*)
      REAL*8 RBUF1(*),RBUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(2),A0,0,RBUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           KJREP = MULTB(KREP,JREP,2)
           LJREP = MULTB(LREP,JREP,2)
           ILREP = MULTB(IREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILKJ0 = OFF(KJREP) 
     &          + (OFF2(KREP,JREP)+K-1)*NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           IKLJ0 = OFF(LJREP) 
     &          + (OFF2(LREP,JREP)+L-1)*NPAIR(IKREP)
           IKLJ0 = IKLJ0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
            DO J = 1, JFIE(JREP)
             DO I = 1, IFIE(IREP)
              IJKL = IJKL + 1
              ILKJ = ILKJ0 + I
              IKLJ = IKLJ0 + I
              IF (CARITH) THEN
                 IF (KJREP.EQ.1) CBUF2(ILKJ) = CBUF1(IJKL)
                 IF (LJREP.EQ.1) CBUF2(IKLJ) = - CBUF1(IJKL)
              ELSE
                 IF (KJREP.EQ.1) RBUF2(ILKJ) = RBUF1(IJKL)
                 IF (LJREP.EQ.1) RBUF2(IKLJ) = - RBUF1(IJKL)
              ENDIF
            ENDDO
            ILKJ0 = ILKJ0 + NPAIR(ILREP) * KFIE(KREP)
            IKLJ0 = IKLJ0 + NPAIR(IKREP) * LFIE(LREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT9 (NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,NTRIPL,
     &                 OFF,OFF1,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I,J;KL:KLREP) to array -BUF2*(KL J,I:IREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IOFF = (OFF(IREP)+OFF1(KLREP,JREP)) * RCW + 1
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            KLJI = IOFF+((I-1)*NTRIPL(IREP)+(J-1)*N)*RCW
            CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(KLJI),1)
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO
      CALL XSCAL (OFF(NREP+1),-A1,BUF2,1)
      IF (CARITH) CALL CONJUGA (OFF(NREP+1),BUF2,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT10 (NREP,MULTB,NPAIR,KFIE,LFIE,NTRIPL,OFF,OFF1,
     &                 BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF2(IJ K,L:LREP) to array BUF1(IJ,K>L:KLREP) and
C     antisymmetrize K and L. Inverse of SRT4
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       N = NPAIR(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        LOFF = (OFF(LREP) + OFF1(IJREP,KREP)) * RCW
        KOFF = (OFF(KREP) + OFF1(IJREP,LREP)) * RCW
        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
            IJK  = (K-1) * N * RCW + 1
            IJKL1 = LOFF + (L-1) * NTRIPL(LREP) * RCW + IJK 
            IJLK1 = KOFF + (K-1) * NTRIPL(KREP) * RCW + IJL 
            CALL XCOPY (N,BUF1(IJKL1),1,BUF2(IJKL),1)
            CALL XAXPY (N,-A1,BUF1(IJLK1),1,BUF2(IJKL),1)
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT11 (NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,OFF,OFF1,
     &                  BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,KL:KLREP) to array -BUF2(JI,KL:KLREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      JIKLOFF = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJKL = OFF(KLREP) * RCW + 1
       JIKLOFF = OFF(KLREP)
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        NIS = IFIE(IREP)
        NJS = JFIE(JREP)
        DO J = 1, NJS
         DO I = 1, NIS
           JIKL = (JIKLOFF+OFF1(JREP,IREP)+(I-1)*NJS+J-1)*RCW + 1
           CALL XCOPY (N,BUF1(IJKL),M,BUF2(JIKL),M)
           IJKL = IJKL + RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      CALL XSCAL (OFF(NREP+1),-A1,BUF2,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT12 (NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,OFF,OFF2,
     &                  BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I,J KL:JKLREP) to array BUF2(I>J,KL:JKLREP)
C     Inverse of SRT5
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF (IREP.LT.JREP) GOTO 10
        JKLREP = IREP
        IKLREP = JREP
        JKLOFF = (OFF(JKLREP)+OFF2(JREP,KLREP)*IFIE(IREP)) * RCW
        IKLOFF = (OFF(IKLREP)+OFF2(IREP,KLREP)*JFIE(JREP)) * RCW
        M2 = IFIE(IREP) * JFIE(JREP)
        DO J = 1, JFIE(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, IFIE(IREP)
            IJKL1 = JKLOFF + ((J-1) * IFIE(IREP) + I - 1) * RCW + 1
            JIKL1 = IKLOFF + ((I-1) * JFIE(JREP) + J - 1) * RCW + 1
            CALL XCOPY(N,BUF1(IJKL1),M2,BUF2(IJKL+IJ),M)
            CALL XAXPY(N,-A1,BUF1(JIKL1),M2,BUF2(IJKL+IJ),M)
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT13(NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,OFF,OFF2,
     &                 BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I J,KL:KLREP) to array BUF2(I,J KL:JKLREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        JKLREP = IREP
        JKLOFF = (OFF(JKLREP)+OFF2(JREP,KLREP)*IFIE(IREP)) * RCW
        M2 = IFIE(IREP) * JFIE(JREP)
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            IJKL1 = JKLOFF + ((J-1) * IFIE(IREP) + I - 1) * RCW + 1
            CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(IJKL1),M2)
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT14 (NREP,NPAIR1,NPAIR2,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,KL:KLREP) to array BUF2*(KL,IJ:KLREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP,NPAIR1(NREP),NPAIR2(NREP)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IJKL = 1
      KLIJ = 1
      DO KLREP = 1, NREP
       M = NPAIR1(KLREP)
       N = NPAIR2(KLREP)
       DO IJ = 1, M
         CALL XCOPY (N,BUF1(IJKL),M,BUF2(KLIJ),1)
         IJKL = IJKL + RCW
         KLIJ = KLIJ + N * RCW
       ENDDO
       IJKL = KLIJ
      ENDDO
C
      N = (KLIJ - 1) / RCW
      IF (CARITH) CALL CONJUGA (N,BUF2,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT15(NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,OFF,OFF2,
     &                 BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I,JKL:JKLREP) to array BUF2(I J,KL:KLREP)
C     Inverse of SRT 13
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IJKL = 0
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        JKLREP = IREP
        JKLOFF = (OFF(JKLREP)+OFF2(JREP,KLREP)*IFIE(IREP)) * RCW
        M2 = IFIE(IREP) * JFIE(JREP)
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            IJKL1 = JKLOFF + ((J-1) * IFIE(IREP) + I - 1) * RCW + 1
            CALL XCOPY(N,BUF1(IJKL1),M2,BUF2(IJKL+IJ),M)
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT16 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                  OFF,OFF1,OFF2,RBUF1,RBUF2,CBUF1,CBUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IL,KJ:KJREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 CBUF1(*),CBUF2(*)
      REAL*8 RBUF1(*),RBUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,RBUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        DO L = 1, LFIE(LREP)
         DO K = 1, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           KJREP = MULTB(KREP,JREP,2)
           ILREP = MULTB(IREP,LREP,2)
           ILKJ0 = OFF(KJREP) 
     &          + (OFF2(KREP,JREP)+K-1)*NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
            DO J = 1, JFIE(JREP)
             DO I = 1, IFIE(IREP)
              IJKL = IJKL + 1
              ILKJ = ILKJ0 + I
              IF (CARITH) THEN
                 CBUF2(ILKJ) = CBUF1(IJKL)
              ELSE
                 RBUF2(ILKJ) = RBUF1(IJKL)
              ENDIF
            ENDDO
            ILKJ0 = ILKJ0 + NPAIR(ILREP) * KFIE(KREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT17 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                  OFF,OFF1,OFF2,RBUF1,RBUF2,CBUF1,CBUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I>J,K>L:KLREP) to array BUF2(IK,JL:JLREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 CBUF1(*),CBUF2(*)
      REAL*8 RBUF1(*),RBUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,RBUF2,1)
      IJKL = 0
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKJL0 = OFF(IKREP) 
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           ILJK0 = OFF(ILREP) 
     &          + (OFF2(JREP,KREP)+(K-1)*JFIE(JREP)) * NPAIR(ILREP)
           ILJK0 = ILJK0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           JKIL0 = OFF(JKREP) 
     &          + (OFF2(IREP,LREP)+(L-1)*IFIE(IREP)) * NPAIR(JKREP)
           JKIL0 = JKIL0 + OFF1(JREP,KREP) + (K-1)*JFIE(JREP) + 1
           JLIK0 = OFF(JLREP) 
     &          + (OFF2(IREP,KREP)+(K-1)*IFIE(IREP)) * NPAIR(JLREP)
           JLIK0 = JLIK0 + OFF1(JREP,LREP) + (L-1)*JFIE(JREP) + 1
            DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             DO I = IMIN, IFIE(IREP)
              IJKL = IJKL + 1
              IKJL = IKJL0 + I
              ILJK = ILJK0 + I
              JKIL = JKIL0 + (I-1)*NPAIR(JKREP)
              JLIK = JLIK0 + (I-1)*NPAIR(JLREP)
              IF (CARITH) THEN
                 CBUF2(IKJL) = CBUF1(IJKL)
                 CBUF2(JLIK) = CBUF1(IJKL)
                 CBUF2(ILJK) = - CBUF1(IJKL)
                 CBUF2(JKIL) = - CBUF1(IJKL)
              ELSE
                 RBUF2(IKJL) = RBUF1(IJKL)
                 RBUF2(JLIK) = RBUF1(IJKL)
                 RBUF2(ILJK) = - RBUF1(IJKL)
                 RBUF2(JKIL) = - RBUF1(IJKL)
              ENDIF
            ENDDO
            IKJL0 = IKJL0 + NPAIR(IKREP)
            ILJK0 = ILJK0 + NPAIR(ILREP)
            JKIL0 = JKIL0 + 1
            JLIK0 = JLIK0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT18 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                  OFF,OFF1,OFF2,RBUF1,RBUF2,CBUF1,CBUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IK,JL:JLREP) to array BUF2(I>J,K>L:KLREP)
C     Inverse of SRT17
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 CBUF1(*),CBUF2(*)
      REAL*8 RBUF1(*),RBUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IJKL = 0
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           IF (IREP.LT.JREP) GOTO 20
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKJL0 = OFF(IKREP) 
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           ILJK0 = OFF(ILREP) 
     &          + (OFF2(JREP,KREP)+(K-1)*JFIE(JREP)) * NPAIR(ILREP)
           ILJK0 = ILJK0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           JKIL0 = OFF(JKREP) 
     &          + (OFF2(IREP,LREP)+(L-1)*IFIE(IREP)) * NPAIR(JKREP)
           JKIL0 = JKIL0 + OFF1(JREP,KREP) + (K-1)*JFIE(JREP) + 1
           JLIK0 = OFF(JLREP) 
     &          + (OFF2(IREP,KREP)+(K-1)*IFIE(IREP)) * NPAIR(JLREP)
           JLIK0 = JLIK0 + OFF1(JREP,LREP) + (L-1)*JFIE(JREP) + 1
            DO J = 1, JFIE(JREP)
             IMIN = 1
             IF (IREP.EQ.JREP) IMIN = J + 1
             DO I = IMIN, IFIE(IREP)
              IJKL = IJKL + 1
              IKJL = IKJL0 + I
              ILJK = ILJK0 + I
              JKIL = JKIL0 + (I-1)*NPAIR(JKREP)
              JLIK = JLIK0 + (I-1)*NPAIR(JLREP)
              IF (CARITH) THEN
                 CBUF2(IJKL) = CBUF1(IKJL) + CBUF1(JLIK) 
     &                       - CBUF1(JKIL) - CBUF1(ILJK)
              ELSE
                 RBUF2(IJKL) = RBUF1(IKJL) + RBUF1(JLIK) 
     &                       - RBUF1(JKIL) - RBUF1(ILJK)
              ENDIF
            ENDDO
            IKJL0 = IKJL0 + NPAIR(IKREP)
            ILJK0 = ILJK0 + NPAIR(ILREP)
            JKIL0 = JKIL0 + 1
            JLIK0 = JLIK0 + 1
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT19 (NREP,MULTB,NPAIR,KFIE,LFIE,NTRIPL,OFF,OFF1,
     &                  BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,KL:KLREP) to array BUF2(IJ L,K:KREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       N = NPAIR(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        KOFF = (OFF(KREP) + OFF1(IJREP,LREP)) * RCW
        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         DO K = 1, KFIE(KREP)
            IJLK1 = KOFF + (K-1) * NTRIPL(KREP) * RCW + IJL 
            CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJLK1),1)
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT20 (NREP,MULTB,IJPAIR,IFIE,JFIE,KFIE,LFIE,
     &                  NPAIR,OFF,OFF1,OFF2,BUF1,NBUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IK,LJ:LJREP=1)
C     Integrals are fetched from disk
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF1
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP,MULTB(16,16,2),IJPAIR(NREP)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (OFF(2),A0,0,BUF2,1)
      DO 10 KLREP = 1, NREP
       IF (IJPAIR(KLREP).EQ.0) GOTO 10
       IJKL = 1
       KL = 0
       KLMAX = NBUF1 / IJPAIR(KLREP)
       KLSTAR = 0
       KLINCR = KLMAX
       CALL GETVOVV (KLREP,KLSTAR,KLINCR,DONE,BUF1)
       DO 20 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 20
        DO 30 L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO 40 K = KMIN, KFIE(KREP)
          KL = KL + 1
          IF (KL.EQ.KLINCR+1) THEN
C -----------------------
C We need the next buffer
C -----------------------
             KLSTAR = KLSTAR + KLINCR
             KLINCR = KLMAX
             CALL GETVOVV (KLREP,KLSTAR,KLINCR,DONE,BUF1)
             KL = 1
             IJKL = 1
          ENDIF
          DO 50 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           LJREP = MULTB(LREP,JREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           KJREP = MULTB(KREP,JREP,2)
           IKLJ0 = OFF(LJREP) + (OFF2(LREP,JREP)+L-1) * NPAIR(IKREP)
           IKLJ0 = IKLJ0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           IKLJ0 = IKLJ0 * RCW + 1
           ILKJ0 = OFF(KJREP) + (OFF2(KREP,JREP)+K-1) * NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           ILKJ0 = ILKJ0 * RCW + 1
           DO 60 J = 1, JFIE(JREP)
            IF (LJREP.EQ.1) THEN
               CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKLJ0),1)
            ENDIF
            IF (KJREP.EQ.1) THEN 
               CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(ILKJ0),1)
               CALL XSCAL (IFIE(IREP),-A1,BUF2(ILKJ0),1)
            ENDIF
            IJKL = IJKL + IFIE(IREP) * RCW
            IKLJ0 = IKLJ0 + NPAIR(IKREP) * LFIE(LREP) * RCW
            ILKJ0 = ILKJ0 + NPAIR(ILREP) * KFIE(KREP) * RCW
 60        CONTINUE
 50       CONTINUE
 40      CONTINUE
 30     CONTINUE
 20    CONTINUE
 10   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT21(NREP,MULTB,NPAIR1A,NPAIR1B,IFIE,JFIE,NPAIR2,
     &                 OFF,OFF1,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(I>J,KL:KLREP) to array BUF2(I J,KL:KLREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR1A(NREP),NPAIR1B(NREP),NPAIR2(NREP)
      INTEGER IFIE(NREP),JFIE(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL1 = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       M1 = NPAIR1A(IJREP)
       M2 = NPAIR1B(IJREP)
       N = NPAIR2(KLREP)
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF (IREP.LT.JREP) GOTO 10
        DO J = 1, JFIE(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, IFIE(IREP)
            IJKL2 = (OFF(KLREP)+OFF1(IREP,JREP)
     &               +(J-1)*IFIE(IREP)+I-1)*RCW + 1
            JIKL2 = (OFF(KLREP)+OFF1(JREP,IREP)
     &               +(I-1)*JFIE(JREP)+J-1)*RCW + 1
            CALL XCOPY(N,BUF1(IJKL1),M1,BUF2(IJKL2),M2)
            CALL XCOPY(N,BUF1(IJKL1),M1,BUF2(JIKL2),M2)
            CALL XSCAL(N,-A1,BUF2(JIKL2),M2)
            IJKL1 = IJKL1 + RCW
         ENDDO
        ENDDO
 10    CONTINUE
      IJKL1 = IJKL1 + M1 * (N-1) * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT22 (NREP,MULTB,NPAIR,KFIE,LFIE,NTRIPL,OFF,OFF1,
     &                  BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(L IJ,K:KREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       N = NPAIR(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        KOFF = OFF(KREP) + OFF1(LREP,IJREP)
        LOFF = OFF(LREP) + OFF1(KREP,IJREP)
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
            LIJK = (KOFF + (K-1) * NTRIPL(KREP) + L - 1) * RCW + 1 
            KIJL = (LOFF + (L-1) * NTRIPL(LREP) + K - 1) * RCW + 1 
            CALL XCOPY (N,BUF1(IJKL),1,BUF2(LIJK),LFIE(LREP))
            CALL XCOPY (N,BUF1(IJKL),1,BUF2(KIJL),KFIE(KREP))
            CALL XSCAL (N,-A1,BUF2(KIJL),KFIE(KREP))
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT23 (NREP,MULTB,NPAIR1,IFIE,JFIE,NPAIR2,NTRIPL,
     &                  OFF,OFF1,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array VOVV(I,J;KL:KLREP) to array BUF2(I KL, J:JREP)
C     Read first array from file, representationwise
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
      INCLUDE "files.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      INTEGER ISTART,NINT
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      DO KLREP = 1, NREP
       IJREP = KLREP
       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       ISTART = 0
       NINT = N
       CALL GETVOVV (IJREP,ISTART,NINT,DONE,BUF1)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        JOFF = (OFF(JREP) + OFF1(IREP,KLREP)) * RCW
        DO J = 1, JFIE(JREP)
         DO I = 1, IFIE(IREP)
            IKLJ1 = JOFF + ((J-1)*NTRIPL(JREP)+I-1) * RCW + 1
            CALL XCOPY(N,BUF1(IJ),M,BUF2(IKLJ1),IFIE(IREP))
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT24 (NREP,MULTB,NPAIR,KFIE,LFIE,NTRIPL,OFF,OFF1,
     &                 BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF2(IJ K,L:LREP) to array BUF1(IJ,K L:KLREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER NPAIR(NREP),KFIE(NREP),LFIE(NREP),NTRIPL(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      IJKL = 1
      DO KLREP = 1, NREP
       IJREP = KLREP
       N = NPAIR(KLREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        LOFF = (OFF(LREP) + OFF1(IJREP,KREP)) * RCW
        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         DO K = 1, KFIE(KREP)
            IJK  = (K-1) * N * RCW + 1
            IJKL1 = LOFF + (L-1) * NTRIPL(LREP) * RCW + IJK 
            CALL XCOPY (N,BUF1(IJKL1),1,BUF2(IJKL),1)
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT25 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                  OFF,OFF1,OFF2,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF1(IJ,K>L:KLREP) to array BUF2(IK,JL:JLREP)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY(OFF(2),A0,0,BUF2,1)
      IJKL = 1
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           JLREP = MULTB(JREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILREP = MULTB(IREP,LREP,2)
           JKREP = MULTB(JREP,KREP,2)
           IKJL0 = OFF(JLREP) 
     &          + (OFF2(JREP,LREP)+(L-1)*JFIE(JREP)) * NPAIR(IKREP)
           IKJL0 = IKJL0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
           IKJL0 = IKJL0 * RCW + 1
           ILJK0 = OFF(JKREP) 
     &          + (OFF2(JREP,KREP)+(K-1)*JFIE(JREP)) * NPAIR(ILREP)
           ILJK0 = ILJK0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           ILJK0 = ILJK0 * RCW + 1
           DO J = 1, JFIE(JREP)
            CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(IKJL0),1)
            CALL XCOPY (IFIE(IREP),BUF1(IJKL),1,BUF2(ILJK0),1)
            CALL XSCAL (IFIE(IREP),-A1,BUF2(ILJK0),1)
            IJKL = IJKL + IFIE(IREP) * RCW
            IKJL0 = IKJL0 + NPAIR(IKREP) * RCW
            ILJK0 = ILJK0 + NPAIR(ILREP) * RCW
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT26 (NREP,MULTB,IFIE,JFIE,KFIE,LFIE,NPAIR,
     &                  OFF,OFF1,OFF2,RBUF1,RBUF2,CBUF1,CBUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Sort array BUF2(IL,KJ:KJREP=1) to array BUF1(IJ,K>L:KLREP) and
C     anti-symmetrize. Inverse of SRT8
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 CBUF1(*),CBUF2(*)
      REAL*8 RBUF1(*),RBUF2(*)
      INTEGER NREP, MULTB(16,16,2)
      INTEGER IFIE(NREP),JFIE(NREP),KFIE(NREP),LFIE(NREP),NPAIR(NREP)
      INTEGER OFF(NREP+1),OFF1(8,8),OFF2(8,8)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IJKL = 0
      DO KLREP = 1, NREP
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10
        DO L = 1, LFIE(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
          DO 20 JREP = 1, NREP
           IREP = MULTB(JREP,KLREP+NREP,2)
           KJREP = MULTB(KREP,JREP,2)
           LJREP = MULTB(LREP,JREP,2)
           ILREP = MULTB(IREP,LREP,2)
           IKREP = MULTB(IREP,KREP,2)
           ILKJ0 = OFF(KJREP) 
     &          + (OFF2(KREP,JREP)+K-1)*NPAIR(ILREP)
           ILKJ0 = ILKJ0 + OFF1(IREP,LREP) + (L-1)*IFIE(IREP)
           IKLJ0 = OFF(LJREP) 
     &          + (OFF2(LREP,JREP)+L-1)*NPAIR(IKREP)
           IKLJ0 = IKLJ0 + OFF1(IREP,KREP) + (K-1)*IFIE(IREP)
            DO J = 1, JFIE(JREP)
             DO I = 1, IFIE(IREP)
              IJKL = IJKL + 1
              ILKJ = ILKJ0 + I
              IKLJ = IKLJ0 + I
              IF (CARITH) THEN
                 CBUF2(IJKL) = CBUF1(ILKJ) - CBUF1(IKLJ)
              ELSE
                 RBUF2(IJKL) = RBUF1(ILKJ) - RBUF1(IKLJ)
              ENDIF
            ENDDO
            ILKJ0 = ILKJ0 + NPAIR(ILREP) * KFIE(KREP)
            IKLJ0 = IKLJ0 + NPAIR(IKREP) * LFIE(LREP)
           ENDDO
 20       CONTINUE
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
C
      RETURN
      END
