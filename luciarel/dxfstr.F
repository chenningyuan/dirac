      SUBROUTINE DXFSTR(IKIN1,IKIN2,IKST,IKSG,NIK,
     &                  KJIN1,KJIN2,KJST,KJSG,NKJ,
     &                  NI,NJ,NK,NL,
     &                  IOFF,JOFF,KOFF,LOFF,
     &                  ISTR,MXSXST,XIJKL,
     *                  IEXSTR,FACSTR,NEX,NTESTG)
*. Two list of single excitations occurs
*
*  IK : String I => String K (IKIN1,IKIN2,IKST,IKSG)
*  KJ : String K => String J (KJIN1,KJIN2,KJST,KJSG)
*
* Obtain effect of two-electron excitation part of two-electron operator
*
*
* g = 1/2sum(i,j,k,l) a+ia+ka l a j (ij!kl)
*   = sum(i>k,j>l)    a+ia+ka l a j ((ij!kl) -(il!jk))
*   = sum(i>k,j>l,j.ne.k)    a+ia j a+k a l ((ij!kl) -(il!jk))
*   + sum(i>k,k>l)    a+ia k a+k a l ((ik!kl) -(il!kk))
*   = sum(i>k,j>l,j.ne.k)    a+ia j a+k a l ((ij!kl) -(il!jk))
*   - sum(i>k,k>l)    a+ia+k a k a l ((ik!kl) -(il!kk))
*   = sum(i>k,j>l,j.ne.k)    a+ia j a+k a l ((ij!kl) -(il!jk))
*   - sum(i>k,k>j)    a+ia j a+k a k ((ik!kj) -(ij!kk))
*
* COLD :
*   = sum(i>k,j.ne.k, j.gt.l)
*     a+ia j a+k a l ((ij!kl) -(il!jk))
*   + sum(i>k, j>l)
*     a+ia j a+j a l ((ij!jl) -(il!jj))
*
* as
*   G !ISTR> = sum(i>k,j.ne.k,if(k.ne.l ) j.gt.l)
*     a+ia ji!KSTR><KSTR! a+k a l!ISTR> ((ij!kl) -(il!jk))
*
* Version only employing single annihilation tables
*
* =====
* Input
* =====
*
* IKIN1 : Index one of <K STR! SX ! ISTR > i.e k
* IKIN2 : Index two of <K STR! SX ! ISTR > i.e l
* IKST  : Corresponding K string
* IKSG :  Phase shift of excitation
* NIK  : Number of excitations per string I=> K
*
* KJIN1 : Index one of <J STR! SX ! KSTR > i.e i
* KJIN2 : Index two of <J STR! SX ! KSTR > i.e j
* KJST  : Corresponding J string
* KJSG :  Phase shift of excitation
* NKJ  : Number of excitations per string K => J
*
* NI,NJ,NK,NL : Number of orbitals in each group
* IOFF,JOFF,KOFF,LOFF : Offset for each group of orbitals.
*
*
* XIJKL : Integrals stored as (j,l,i,k)
* ======
* Output
* ======
*
*     NEX    : Number of excitations
*
*     IEXSTR : IEXSTR(K) is string number of excited string K
*              relative to start of TS block
*     FACSTR : FACSTR(K) is the factor that string INDSTR(K) is multiplied
*              with
*
*     Jeppe Olsen, March 1994, LUCIA version of DXFSTR
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IKIN1(MXSXST,*),IKIN2(MXSXST,*)
      INTEGER IKST(MXSXST,*),NIK(*)
      REAL*8  IKSG(MXSXST,*)
*
      INTEGER KJIN1(MXSXST,*),KJIN2(MXSXST,*)
      INTEGER KJST(MXSXST,*),NKJ(*)
      REAL*8  KJSG(MXSXST,*)
*
      DIMENSION XIJKL(NJ,NL,NI,NK)
*
*. Output
      INTEGER IEXSTR(*)
      DIMENSION FACSTR(*)
*
      CALL QENTER('DXFST')
*
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ****************************************'
        WRITE(6,*) ' DXFSTR at your command : ISTR = ',ISTR
        WRITE(6,*) ' ****************************************'
      END IF
*
      NEX = 0
C?    WRITE(6,*) ' Number of excitations from string ', NIK(ISTR)
      DO 1000 KLEX = 1, NIK(ISTR)
        K = IKIN1(KLEX,ISTR)
        L = IKIN2(KLEX,ISTR)
        KSTR = IKST(KLEX,ISTR)
        SIGN1 = IKSG(KLEX,ISTR)
C?      WRITE(6,*) ' KLEX K L ', KLEX,K,L
C?      WRITE(6,*) ' KSTR SIGN1 ',KSTR,SIGN1
        DO 900 IJEX = 1, NKJ(KSTR)
          I = KJIN1(IJEX,KSTR)
          J = KJIN2(IJEX,KSTR)
C?        WRITE(6,*) ' I J ', I,J
          ICASE = 0
          IF(I.GT.K .AND. J.GT.L .AND.
     &       I.NE.L.AND.I.NE.J.AND.K.NE.L.AND.K.NE.J) THEN
            ICASE = 1
          END IF
          IF(ICASE.EQ.0) GOTO 900
          JSTR = KJST(IJEX,KSTR)
C?        write(6,*) ' KSTR = ', KSTR
C?        write(6,'(A,5I3)')
C?   &   ' Allowed excitation i,j,k,l  JSTR ',I,J,K,L,JSTR
          IF(KJSG(IJEX,KSTR).EQ.1.0D0) THEN
            SIGN =   SIGN1
          ELSE
            SIGN = - SIGN1
          END IF
C?        write(6,*) ' jstr,kstr sign ', jstr,kstr,sign
*
          NEX = NEX + 1
          IEXSTR(NEX) = JSTR
          IF(ICASE.EQ.1) THEN
C    &      SIGN*XIJKL(J-JOFF+1,L-LOFF+1,I-IOFF+1,K-KOFF+1)
            FACSTR(NEX) =
     &      SIGN*(GTIJKL(I,J,K,L)-GTIJKL(I,L,K,J))
C    &      SIGN*XIJKL(J-JOFF+1,L-LOFF+1,I-IOFF+1,K-KOFF+1)
          END IF
C         IF(FORM1.NE.FACSTR(NEX) ) THEN
C           write(6,*) ' problemo in DXFSTR'
C           write(6,*) ' i j k l Form1 FACSTR '
C           write(6,*) i,j,k,l,FORM1,FACSTR(NEX)
C           write(6,*) ' ijkl and ilkj '
C           write(6,*) GTIJKL(I,J,K,L),GTIJKL(I,L,K,J)
C         END IF
  900   CONTINUE
 1000 CONTINUE
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 'Number of excitations obtained ', NEX
        WRITE(6,*)
        WRITE(6,*) ' iex      exc.string         factor '
        WRITE(6,*) ' =================================='
        DO IEX = 1, NEX
          WRITE(6,'(I4,4X,I8,6X,F13.5)')
     &    IEX,IEXSTR(IEX),FACSTR(IEX)
        END DO
      END IF
*
      CALL QEXIT ('DXFST')
      RETURN
      END
