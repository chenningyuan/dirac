      SUBROUTINE ADS1_GAS(NK,I1,XI1S,LI1,IORB,LORB,
     &                KSTR,NKEL,NKSTR,IREO,IZ,
     &                NOCOB,KMAX,KMIN,IEND,KACT,SCLFAC)
*
* Obtain I1(KSTR) = +/- A+ IORB !KSTR>
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* IORB : First orbital to be added
* LORB : Number of orbitals to be added : IORB to IORB-1+LORB
*        are used. They must all be in the same TS group
* KMAX : Largest allowed relative number for K strings
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR,JORB) : ne. 0 => a + JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR,JORB) : above +/-
*          : eq. 0    a + JORB !KSTR> = 0
* Offset is KMIN
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER KSTR(NKEL,NKSTR)
      INTEGER IREO(*), IZ(NOCOB,*)
*.Output
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*
      NTEST = 000
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ================= '
       WRITE(6,*) ' ADS1_GAS speaking '
       WRITE(6,*) ' ================= '
       WRITE(6,*) ' IORB,LORB ', IORB,LORB
      END IF
C?    WRITE(6,*) ' IZ matrix in ADS1 '
C?    CALL IWRTMA(IZ,NOCOB,NKEL+1,NOCOB,NKEL+1)
*
      KACT = 0
      IORBMIN = IORB
      IORBMAX = IORB + LORB - 1
*
      KEND = MIN(NKSTR,KMAX)
      NK = MAX(0,KEND-KMIN+1)
      IF(KEND.GE.NKSTR) THEN
        IEND = 1
      ELSE
        IEND = 0
      END IF
*. Start by zeroing
      ZERO = 0.0D0
      IZERO = 0
      DO IIORB = 1, LORB
C       CALL SETVEC(XI1S(1,IIORB),ZERO,NK)
        CALL ISETVC(I1(1,IIORB),IZERO,NK)
      END DO
*
      DO KKSTR = KMIN,KEND
C?    write(6,'(a,I12)') 'KKSTR    ',KKSTR
C?    WRITE(6,*) ' Occ of K ', (KSTR(JEL,KKSTR),JEL=1,NKEL)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Occupation of string ', KKSTR
         CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
       END IF
*
       SGN = -1.0D0*SCLFAC
       DO KEL = 1, NKEL+1
         SGN = - SGN
         IF(KEL.LE.NKEL) THEN
           KORB = KSTR(KEL,KKSTR)
         ELSE
           KORB = IORB+LORB
         END IF
*
         IF(KORB.GT.IORBMIN) THEN
*. Orbitals can be added before electron KEL
           IORB2 = MIN(KORB-1,IORBMAX)
           IF(KEL.EQ.1) THEN
             IORB1 = IORBMIN
           ELSE
             IORB1 = MAX(IORBMIN,KSTR(KEL-1,KKSTR)+1)
           END IF
           IF(IORB1.GT.IORBMAX) GOTO 2912
           IF(IORB2.GE.IORB1) THEN
*
*. common piece of lexical ordering
C            IF(KEL.EQ.1) THEN
C              IREVLEX = 1
C              DO JEL = 1, NKEL
C                 IREVLEX = IREVLEX + IZ(KSTR(JEL,KKSTR),JEL+1)
C              END DO
C            ELSE
C                 IREVLEX = IREVLEX + IZ(KSTR(KEL-1,KKSTR),KEL-1)
C    &                              - IZ(KSTR(KEL-1,KKSTR),KEL)
C            END IF
C          IF(IORB2.GE.IORB1) THEN
*
             IREVLEX = 1
             DO JEL = 1, KEL-1
               IREVLEX = IREVLEX + IZ(KSTR(JEL,KKSTR),JEL)
             END DO
             DO JEL = KEL, NKEL
               IREVLEX = IREVLEX + IZ(KSTR(JEL,KKSTR),JEL+1)
             END DO
             IF(NTEST.GE.1000) WRITE(6,*)
     &       'KEL IORB1,IORB2,IREVLEX ',KEL,IORB1,IORB2,IREVLEX
             DO IIORB = IORB1, IORB2
               XI1S(KKSTR-KMIN+1,IIORB-IORBMIN+1) = SGN
               I1(KKSTR-KMIN+1,IIORB-IORBMIN+1) =
     &         IREO(IREVLEX  + IZ(IIORB,KEL))
*
               KACT = 1
*
               IF( I1(KKSTR-KMIN+1,IIORB-IORBMIN+1).LE.0) THEN
                 WRITE(6,*)
     &           ' Problems with map, stringnumber .le.0'
                 write(6,*) 'KKSTR         ',KKSTR
                 WRITE(6,*) ' Occupation of K string '
                 CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
                 WRITE(6,*) ' Occupation of previous strings '
                 CALL IWRTMA(KSTR,NKEL,KKSTR-1,NKEL,KKSTR-1)
                 write(6,*) 'KMIN          ',KMIN
                 write(6,*) 'IIORB         ',IIORB
                 write(6,*) 'IORBMIN       ',IORBMIN
                 write(6,*) 'KEL           ',KEL
                 write(6,*) 'IREVLEX       ',IREVLEX
                 write(6,*) 'IZ(IIORB,KEL) ',IZ(IIORB,KEL)
                 write(6,*) 'IREO(IREVLEX  + IZ(IIORB,KEL)) =',
     &                       IREO(IREVLEX  + IZ(IIORB,KEL))
                 Call Abend
     &           ( ' Problems with map, stringnumber .le.0')
               END IF
             END DO
           END IF
         END IF
       END DO
 2912  CONTINUE
      END DO
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADS1_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO 200 IIORB = IORB,IORB+LORB-1
            IIORBR = IIORB-IORB+1
            WRITE(6,*) ' Info for orbital ', IIORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(I1(1,IIORBR),1,NK,1,NK)
            CALL WRTMAT(XI1S(1,IIORBR),1,NK,1,NK)
  200     CONTINUE
        END IF
      END IF
*
C     WRITE(6,*) 'JEPPE forced me to stop in ADS1   '
C     STOP'JEPPE forced me to stop in ADS1   '
      RETURN
      END
