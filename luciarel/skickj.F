      SUBROUTINE SKICKJ(SKII,CKJJ,NKA,NIB,NJB,NKB,XIJKL,
     &     NI,NJ,NK,NL,MAXK,
     &     KBIB,XKBIB,KBJB,XKBJB,IKORD,
     &     FACS,IROUTE,MXTSOB)
*
*
* Calculate S(Ka,Ib,i) = FACS*S(Ka,Ib,i)
*          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
*
* Jeppe Olsen, Spring of 94
*
* : Note : Route 1 has retired, March 97
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CKJJ(*)
      DIMENSION XIJKL(*)
      DIMENSION KBIB(MAXK,*),XKBIB(MAXK,*)
      DIMENSION KBJB(MAXK,*),XKBJB(MAXK,*)
*. Input and output
      DIMENSION SKII(*)
*. Scratch
C     COMMON/SOMESCR/SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
C     DIMENSION IBOFF(MXTSOB*MXTSOB),JBOFF(MXTSOB*MXTSOB)
C     DIMENSION SCR(NI*NL*NJ),IBOFF(NI),JBOFF(NL*NJ)
      DIMENSION XIJILS(MXTSOB)
      IF(NI.GT.MXTSOB.OR.NJ.GT.MXTSOB.OR.NK.GT.MXTSOB
     &     .OR.NL.GT.MXTSOB) THEN
         WRITE(6,*) ' SKICKJ : Too many orbs : NI > MXTSOB '
         WRITE(6,*) ' NI, MXTSOB ',MAX(NI,NJ,NK,NL),MXTSOB
         Call Abend2( ' Redim MXTSOB in SKICKJ' )
      END IF
*
C     CALL QENTER('SKICK')
      IF(IROUTE.EQ.3) THEN
* S(Ka,i,Ib) = S(Ka,i,Ib) + sum(j) (ji!kl) C(Ka,j,Jb)
         DO KB = 1, NKB
*. Number of nonvanishing connections from KB
            LL = 0
            KK = 0
            DO L = 1, NL
               IF(KBJB(KB,L).NE.0) LL = LL + 1
            END DO
            DO K = 1, NK
               IF(KBIB(KB,K).NE.0) KK = KK + 1
            END DO
*
            IF(KK.NE.0.AND.LL.NE.0) THEN
               DO K = 1, NK
                  IB = KBIB(KB,K)
                  IF(IB.NE.0) THEN
                     SGNK = XKBIB(KB,K)
                     DO L = 1, NL
                        JB = KBJB(KB,L)
                        IF(JB.NE.0) THEN
                           SGNL = XKBJB(KB,L)
                           FACTOR = SGNK*SGNL
*. We have now an IB and Jb string, let us do it
                           ISOFF = (IB-1)*NI*NKA + 1
                           ICOFF = (JB-1)*NJ*NKA + 1
C?                         write(6,*) 'ICOFF is now =  ',ICOFF
                           INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
                           IMAX = NI
*
                           IF(IKORD.NE.0) THEN
*. Restrict so (ij) .le. (kl)
                              IMAX  = K
                              JKINTOF = INTOF + (K-1)*NJ
C     CALL COPVEC(XIJKL(JKINTOF),XIJILS,NJ)
                              DO J = L,NL
                                 XIJILS(J) = XIJKL(JKINTOF-1+J)
                              END DO
                              XIJKL(JKINTOF-1+L)
     &                          = 0.5D0*XIJKL(JKINTOF-1+L)
                              DO J = L+1, NL
                                 XIJKL(JKINTOF-1+J) = 0.0D0
                              END DO
                           END IF
C     ONE = 1.0D0
                           CALL MATML7(SKII(ISOFF),CKJJ(ICOFF),
     &                          XIJKL(INTOF),NKA,IMAX,NKA,NJ,
     &                          NJ,IMAX,FACS,FACTOR ,0)
                           IF(IKORD.NE.0) THEN
                              DO J = L,NL
                                 XIJKL(JKINTOF-1+J) =  XIJILS(J)
                              END DO
C     CALL COPVEC(XIJILS,XIJKL(JKINTOF),NJ)
                           END IF
*
                        END IF
*     ^ JB.eq.0 or not
                     END DO
*                    ^ End of L to NL loop
*
                  END IF
*                 ^ IB.eq.0 or not
               END DO
*              ^ End of K to NK loop
            END IF
*           ^ KK, LL.eq.0 or not
         END DO
*        ^ End of loop over KB strings
*
      ELSE IF (IROUTE.EQ.1.OR.IROUTE.EQ.2) THEN

         WRITE(6,*) ' Sorry route 1and 2 have retired'
         Call Abend2( 'SKICKJ:Invalid route=1,2' )
      END IF
*. End of IROUTE branchning
*
C     CALL QEXIT('SKICK')
      RETURN
      END
