      SUBROUTINE SXSTRS(ISTROF,NIST,KSTROF,
     &                   NI,IOFF,NJ,JOFF,
     &                   IAMAPO,IAMAPS,LAMAP,IAMPFL,
     &                   ICMAPO,ICMAPS,LCMAP,ICMPFL,
     &                   ICMPO,ICMPL,
     &                   MXSXST,ISXSTR,JSXSTR,IEXSTR,FACSTR,NEX,
     &                   IOFFDG,NTEST)
*
* Obtain single excitations from string ISTROF-ISTROF+NIST-1,
* Slave routine mastered by SXSTR
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
* ICMPFL : Creation map complete ?
* LCMAP  : Row dimension of Creation     map
*
* Jeppe Olsen, March 1994
      IMPLICIT REAL*8(A-H,O-Z)
*
*. Input
*
      INTEGER IAMAPO(LAMAP,*), IAMAPS(LAMAP,*)
C     INTEGER ICMAPO(LCMAP,*), ICMAPS(LCMAP,*)
      INTEGER ICMAPO(*),ICMAPS(*)
      INTEGER ICMPO(*),ICMPL(*)
*. Output
      INTEGER ISXSTR(MXSXST,*),JSXSTR(MXSXST,*)
      INTEGER IEXSTR(MXSXST,*)
      DIMENSION FACSTR(MXSXST,*)
      INTEGER NEX(*)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' SXFSTS : NTEST = ', NTEST
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
            JISTRO = (JISTR-1)*LCMAP+1
            NIEFF = LCMAP
          ELSE
            JISTRO = ICMPO(JISTR)
            NIEFF = ICMPL(JISTR)
          END IF
C?        WRITE(6,*) ' JISTRO NIEFF ', JISTRO,NIEFF
C         DO 1001 IORB = IOFF,IOFF+NI-1
          DO 1001 IIORB = 1, NIEFF
*
            IF(ICMPFL.EQ.1) THEN
              IORB = IIORB-1+IOFF
              IIIORB = IORB
            ELSE
              IORB = ABS(ICMAPO(JISTRO-1+IIORB))
              IIIORB = IIORB
              IF(IORB.LT.IOFF .OR. IORB .GT. IOFF + NI -1 )
     &        GOTO 1001
            END IF
*
C?          WRITE(6,*) ' IORB = ', IORB
            IF(IOFFDG.NE.0 .AND. IORB.EQ.JORB) GOTO 1001
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
C?          WRITE(6,*) ' IIIORB IJISTR ', IIIORB, IJISTR
C?          WRITE(6,*) ' IJACT = ', IJACT
            IF(IJACT.EQ.0) GOTO 1001
*
*. A new excitation has been born, enlist it !
*
            LEX = LEX + 1
            ISXSTR(LEX,ISTR-ISTROF+1) = IORB
            JSXSTR(LEX,ISTR-ISTROF+1) = JORB
            IEXSTR(LEX,ISTR-ISTROF+1) = IJISTR-KSTROF+1
            FACSTR(LEX,ISTR-ISTROF+1) = SIJ
C?          WRITE(6,*) 'FACSTR = ',  FACSTR(LEX,ISTR-ISTROF+1)
 1001     CONTINUE
 1002   CONTINUE
        NEX(ISTR-ISTROF+1) = LEX
 1100 CONTINUE
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' **************** '
         WRITE(6,*) ' SXSTRS reporting '
         WRITE(6,*) ' **************** '
         WRITE(6,*)
         DO ISTR = ISTROF,ISTROF+NIST-1
           WRITE(6,*) ' excitations from string ',ISTR
           WRITE(6,*)
           WRITE(6,*) ' iorb jorb exc.string phase '
           WRITE(6,*) ' =========================='
           DO LEX = 1, NEX(ISTR-ISTROF+1)
             WRITE(6,'(2I4,I8,F8.3)')
     &       ISXSTR(LEX,ISTR-ISTROF+1),JSXSTR(LEX,ISTR-ISTROF+1),
     &       IEXSTR(LEX,ISTR-ISTROF+1),FACSTR(LEX,ISTR-ISTROF+1)
           END DO
           WRITE(6,*)
         END DO
      END IF
*
      RETURN
      END
