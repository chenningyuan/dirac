      SUBROUTINE SKICKJ_MK1(SKII,CKJJ,NKA,NIB,NJB,NKB,NMB,XIJKL,
     &                  NI,NJ,NK,NL,MAXI,MAXK,MAXM,
     &                  KBIB,XKBIB,MBKB,XMBKB,JBMB,XJBMB,IKORD,
     &                  FACS,MXTSOB,XIJILS)
*
*
* Calculate S(Ka,Ib) = FACS*S(Ka,Ib)
*          +SUM(i,j,k,l,Kb) <Ib!a  jb!Kb><Kb!a kb !Mb>*<Mb! a lb!Jb>
*           (ij!kl)*C(Ka,i,Jb)
*
* as
*
*          +SUM(i,j,k,l,Kb) <Kb!a  jb!Ib><Mb!a kb !Kb>*<Ib! a lb!Mb>
*           (ij!kl)*C(Ka,i,Jb)
*
* Jeppe Olsen, Winter of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CKJJ(NKA,*)
      DIMENSION XIJKL(*)
*
      DIMENSION KBIB(MAXI,*),XKBIB(MAXI,*)
      DIMENSION MBKB(MAXK,*),XMBKB(MAXK,*)
      DIMENSION JBMB(MAXM,*),XJBMB(MAXM,*)
*
*. Input and output
      DIMENSION SKII(NKA,*)
*. Scratch
      DIMENSION XIJILS(MXTSOB)
      IF(NI.GT.MXTSOB.OR.NJ.GT.MXTSOB.OR.NK.GT.MXTSOB
     &   .OR.NL.GT.MXTSOB) THEN
         WRITE(6,*) ' SKICKJ : Too many orbs : NI > MXTSOB '
         WRITE(6,*) ' NI, MXTSOB ',MAX(NI,NJ,NK,NL),MXTSOB
         Call Abend2( ' Redim MXTSOB in SKICKJ' )
      END IF
*
      DO IB = 1, NIB
        DO J = 1, NJ
        IF(KBIB(IB,J).NE.0) THEN
          SIGNJ = XKBIB(IB,J)
          KB = KBIB(IB,J)
          DO K = 1, NK
          IF(MBKB(KB,K).NE.0) THEN
            MB = MBKB(KB,K)
            SIGNK = XMBKB(KB,K)
            DO L = 1, NL
            IF(JBMB(MB,L).NE.0) THEN
              JB = JBMB(MB,L)
              SIGNL = XJBMB(JB,L)
*. All connections have been made and
*. We have now a IB and Jb string, let's do it
C             ISOFF = (IB-1)*NKA + 1
C             ICOFF = (JB-1)*NI*NKA + 1
              ISOFF = IB
*. C(1,1,JB)
              ICOFF = (JB-1)*NI+ 1
*. Integrals are assumed stored as (ij!kl), offset to (1j!kl)
              INTOFF = ((L-1)*NK + K - 1 )*NI*NJ
     &               +  (J-1)*NI + 1
              CALL MATML7(SKII(1,ISOFF),CKJJ(1,ICOFF),XIJKL(INTOFF),
     &                    NKA,1,NKA,NI,NI,1,0,1.0D0,1.0D0)
            END IF
            END DO
          END IF
          END DO
        END IF
        END DO
      END DO
*
      RETURN
      END
