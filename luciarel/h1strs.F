      SUBROUTINE H1STRS(H,ISTROF,NIST,KSTROF,
     &                   NI,IOFF,NJ,JOFF,
     &                   IAMAPO,IAMAPS,LAMAP,IAMPFL,
     &                   ICMAPO,ICMAPS,ICMPO,ICMPL,LCMAP,ICMPFL,
     &                   MXSXST,IEXSTR,FACSTR,NEX,
     &                   IH2TRM,IPIJKK,XIJKK,IOCSTR,
     &                   NEL,NTOOB,NTEST)
*
* Slave routine mastered by H1STR
* ( See my master for further information about my role in life )
*
* ==================
*. Additional input ( compared to SXSTR)
* ==================
*
* ISTROF : Absolute number of first string to be excited from
* NIST   : Number of strings to be excited from
* KSTROF : Offset of strings in resulting type-symmetry block
* N*,*OFF,*=I,J : Number and offset for each orbital set
*
* IAMAPO : Annihilation mapping, orbital part
* IAMAPS : Annihilation mapping, string part
* LAMAP  : Row dimension of Annihilation map
* IAMPFL : Annihilation map complete ?
*
* ICMAPO : Creation     mapping, orbital part
* ICMAPS : Creation     mapping, string part
* LCMAP  : Row dimension of Creation     map
* ICMPFL : Creation     map complete ?
*
* IH2TRM : ne. 0 : include zero and one-electron excitations from
*                  twobody operator
* IPIJKK : Pointer to symmetry adapted integral h(ij)
* XIJKK  : List of integrals  (ij!kk) - (ik!kj)
* NTOOB  : Number of orbitals
* IOCSTR : Occupation of input strings
* NEL    : Number of electrons in input string
*
* Jeppe Olsen, March 1994
      IMPLICIT REAL*8(A-H,O-Z)
*
*. Input
*
      INTEGER IAMAPO(LAMAP,*), IAMAPS(LAMAP,*)
C     INTEGER ICMAPO(LCMAP,*), ICMAPS(LCMAP,*)
      INTEGER ICMAPO(*),ICMAPS(*)
      INTEGER ICMPL(*),ICMPO(*)
      DIMENSION H(NI,NJ)
      DIMENSION IPIJKK(NTOOB,NTOOB),XIJKK(NTOOB,*)
      DIMENSION IOCSTR(NEL,*)
*. Output
      INTEGER IEXSTR(MXSXST,*)
      INTEGER NEX(*)
      DIMENSION FACSTR(MXSXST,*)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' H1STRS : NTEST = ', NTEST
        WRITE(6,*) ' LCMAP LAMAP ', LCMAP,LAMAP
        WRITE(6,*) ' IAMPFL, ICMPFL ',IAMPFL, ICMPFL
        WRITE(6,*) ' MXSXST ', MXSXST
      END IF
*
      DO 1100 ISTR = ISTROF,ISTROF+NIST-1
C?      write(6,*) ' ISTR = ', ISTR
        LEX = 0
        DO 1002 JORB = JOFF,JOFF+NJ-1
C?        write(6,*) ' JORB = ', JORB
*
* =====================================
* 1 :  Remove orbital JORB from ISTR
* =====================================
*
          JOCC = 0
          IF(IAMPFL.EQ.1) THEN
*. Read from full map
            IF(IAMAPO(JORB,ISTR).EQ.-JORB) THEN
              JOCC = 1
              IF(IAMAPS(JORB,ISTR).GT.0) THEN
                JISTR = IAMAPS(JORB,ISTR)
                SJ = 1.0D0
               ELSE
                JISTR = -IAMAPS(JORB,ISTR)
                SJ = -1.0D0
              END IF
            END IF
          ELSE
*. Read from compact map
            DO JELEC = 1, LAMAP
              IF(IAMAPO(JELEC,ISTR).EQ.-JORB) THEN
                JOCC = 1
                IF(IAMAPS(JELEC,ISTR).GT.0) THEN
                   JISTR = IAMAPS(JELEC,ISTR)
                   SJ = 1.0D0
                ELSE
                   JISTR = -IAMAPS(JELEC,ISTR)
                   SJ = -1.0D0
                END IF
              END IF
            END DO
          END IF
C?        WRITE(6,*) ' JOCC = ', JOCC
          IF(JOCC.EQ.0) GOTO 1002
*
* ==================================
*. 2 : Add orbital I to string JISTR
* ==================================
*
          IF(ICMPFL.EQ.1) THEN
            JISTRO = (JISTR-1)*LCMAP
            NIEFF = NI
          ELSE
            JISTRO = ICMPO(JISTR)
            NIEFF =  ICMPL(JISTR)
          END IF
C         DO 1001 IORB = IOFF,IOFF+NI-1
          DO 1001 IIORB = 1, NIEFF
            IF(ICMPFL.EQ.1) THEN
              IORB = IIORB-1+IOFF
              IIIORB = IORB
            ELSE
              IORB = ABS(ICMAPO(JISTRO-1+IIORB))
              IIIORB = IIORB
              IF(IORB.LT.IOFF.OR.IORB.GT.IOFF+NI-1)
     &        GOTO 1001
            END IF
C?          WRITE(6,*) ' IORB = ', IORB
            IJACT = 0
            IF(ICMAPO(JISTRO-1+IIIORB).EQ.+IORB) THEN
              IJACT  = 1
              IF(ICMAPS(JISTRO-1+IIIORB).GT.0) THEN
                IJISTR = ICMAPS(JISTRO-1+IIIORB)
                SIJ = SJ
              ELSE
                IJISTR = -ICMAPS(JISTRO-1+IIIORB)
                SIJ = -SJ
              END IF
            END IF
C?          WRITE(6,*) ' IJACT = ', IJACT
            IF(IJACT.EQ.0) GOTO 1001
*
*. A new excitation has been born, enlist it !
*
            LEX = LEX + 1
            IEXSTR(LEX,ISTR-ISTROF+1) = IJISTR-KSTROF+1
            FACSTR(LEX,ISTR-ISTROF+1) =
     &      SIJ*H(IORB-IOFF+1,JORB-JOFF+1)
*. If IH2TRM .ne. 0 add
*
* sum (k)  a+i aj a+k a k /1+delta(i,j) (ij!kk)-(ik!kj)
           IF(IORB.EQ.JORB) THEN
             FACIJ = 0.5D0*SIJ
           ELSE
             FACIJ = 1.0D0*SIJ
           END IF
           IJEFF = IPIJKK(IORB,JORB)
*
C          WRITE(6,*) ' TESTING in H1STRS '
C          WRITE(6,*) ' IORB JORB FACIJ IJEFF ',
C    &                  IORB,JORB,FACIJ,IJEFF
           DO KEL = 1, NEL
             KORB = IOCSTR(KEL,ISTR)
             FACSTR(LEX,ISTR-ISTROF+1) =
     &       FACSTR(LEX,ISTR-ISTROF+1) +
     &       FACIJ*XIJKK(KORB,IJEFF)
*
C            WRITE(6,*) ' TESTING in H1STRS '
C            WRITE(6,*) ' ISTR KEL KORB',ISTR,KEL,KORB
C            WRITE(6,*) ' XIJKK(K,IJEFF) and explicit '
C            XIJKK2 = GTIJKL(IORB,JORB,KORB,KORB)-
C    &                GTIJKL(IORB,KORB,KORB,JORB)
C            WRITE(6,*) XIJKK2, XIJKK(KORB,IJEFF)
*
           END DO


 1001     CONTINUE
 1002   CONTINUE
        NEX(ISTR-ISTROF+1) = LEX
 1100 CONTINUE
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' **************** '
         WRITE(6,*) ' H1STRS reporting '
         WRITE(6,*) ' **************** '
         WRITE(6,*)
         DO ISTR = ISTROF,ISTROF+NIST-1
           WRITE(6,*) ' excitations from string ',ISTR
           WRITE(6,*)
           WRITE(6,*) '    exc.string      factor '
           WRITE(6,*) ' =============================='
           DO LEX = 1, NEX(ISTR-ISTROF+1)
             WRITE(6,'(3X,I8,F13.8)')
     &       IEXSTR(LEX,ISTR-ISTROF+1),FACSTR(LEX,ISTR-ISTROF+1)
           END DO
           WRITE(6,*)
         END DO
      END IF
*
      RETURN
      END
