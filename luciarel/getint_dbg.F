      SUBROUTINE GETINT_DBG(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &                      KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &                      IRI,IXCHNG,IKSM,JLSM,ICOUL)
*
*
* Obtain block of double group integrals
*     Jeppe Olsen
*
*  New list structure incorporated
*     Timo Fleig
*
      IMPLICIT REAL*8(A-H,O-Z)
*. output
      DIMENSION XINT(*)
*. Scratch (should be able to hold one integral block)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
*
      INTEGER IREO(4)
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Godmorgen velkommen til GETINT_DBG '
        WRITE(6,*) ' Godmorgen velkommen til GETINT_DBG '
        WRITE(6,*) ' Godmorgen velkommen til GETINT_DBG '
        WRITE(6,*) ' Godmorgen velkommen til GETINT_DBG '
*
        WRITE(6,*) ' ITYP,ISYM,IBU',  ITYP,ISYM,IBU
        WRITE(6,*) ' JTYP,JSYM,JBU',  JTYP,JSYM,JBU
        WRITE(6,*) ' KTYP,KSYM,KBU',  KTYP,KSYM,KBU
        WRITE(6,*) ' LTYP,LSYM,LBU',  LTYP,LSYM,LBU
        WRITE(6,*)
        WRITE(6,*) ' JLSM = ', JLSM
      END IF
*
*
* The lists of integrals are
*
* 1 : (u u | u u)
*--------------------
* 2 : (u b | u u)
* 3 : (b u | u b)
* 4 : (b u | b u)
*--------------------
* 5 : (u b | u b)
* 6 : (b u | u u)
*--------------------
*
* What list of integrals should be accessed
*
      IF(IBU.NE.JBU.AND.KBU.NE.LBU) THEN
*. List 3 or list 4 or list 5
         IF(IBU.EQ.KBU) THEN
*. List 4 or 5, two possibilities : (u b ! u b) (c.c)
*                                   (b u ! b u) (c.c)
            IF(IBU.EQ.2) THEN
* ( U B ! U B) : complex conjugation
*. Output => input index
               IREO(1) = 2
               IREO(2) = 1
               IREO(3) = 3
               IREO(4) = 4
*
               SCALE = 1.0D0
C?             WRITE(6,*) ' GETINT_DBG_S 1 '
               CALL GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &              KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,3)),WORK(KPINT2D(3)),
     &              IREO,SCALE)
*
C                   GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
C    &                           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
C    &                           IXCHNG,IKSM,JLSM,ICOUL,
C    &                           XINTLIST,PINT2D,IREOINDX,SCALE)
            else if (IBU.eq.1) then
*. (B U ! B U) : complex conjugation
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 4
               IREO(4) = 3
               IF(IRI.EQ.1) THEN
                  SCALE = 1.0D0
               ELSE
                  SCALE = -1.0D0
               END IF
C?             WRITE(6,*) ' GETINT_DBG_S 2 '
               CALL GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &              KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,3)),WORK(KPINT2D(3)),
     &              IREO,SCALE)
            END IF
*           ^ End if IBU = 2 or 1
         ELSE
*. List 4 : Two possibilities : (b u ! u b) (directly obtained)
*                               (u b ! b u) (c.c)
            IF(IBU.EQ.1) THEN
*     . Output => input index
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 3
               IREO(4) = 4
*
               SCALE = 1.0D0
C?             WRITE(6,*) ' GETINT_DBG_S 3 '
               CALL GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &                       KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &                       IXCHNG,IKSM,JLSM,ICOUL,
     &                       WORK(KINT2D(IRI,3)),WORK(KPINT2D(3)),
     &                       IREO,SCALE)
            ELSE
*. (U B ! B U) : particle 1-2 interchange
               IREO(1) = 3
               IREO(2) = 4
               IREO(3) = 1
               IREO(4) = 2
*
               SCALE = 1.0D0
C?             WRITE(6,*) ' GETINT_DBG_S 4 '
               CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &                        KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &                        IXCHNG,IKSM,JLSM,ICOUL,
     &                        WORK(KINT2D(IRI,3)),WORK(KPINT2D(3)),
     &                        IREO,SCALE)
            END IF
*         ^ End of the two possibilities
         END IF
*       ^ End of switch list3/list4,5
      END IF
*
*
      IF(IBU.EQ.JBU.AND.KBU.EQ.LBU) THEN
*
* Derived from list 1 : ( U U ! U U )
         IF(JBU.EQ.KBU) THEN
*
* List 1 : ( U U ! U U) or ( B B ! B B)
*
            IF(IBU.EQ.2) THEN
*. ( U U ! U U ) case, simple
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 3
               IREO(4) = 4
*
               SCALE = 1.0D0
C?             WRITE(6,*) ' GETINT_DBG_S 5 '
               CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &              KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,1)),WORK(KPINT2D(1)),
     &              IREO,SCALE)
            ELSE
* ( B B ! B B) case, use kramers symmetry and c.c
               CALL KRA_SYM(ISYM,IBU,ISYMK,IBUK)
               CALL KRA_SYM(JSYM,JBU,JSYMK,JBUK)
               CALL KRA_SYM(KSYM,KBU,KSYMK,KBUK)
               CALL KRA_SYM(LSYM,LBU,LSYMK,LBUK)
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 3
               IREO(4) = 4
               IF(IRI.EQ.1) THEN
                  SCALE = 1.0D0
               ELSE
                  SCALE = -1.0D0
               END IF
C?             WRITE(6,*) ' GETINT_DBG_S 6 '
               CALL  GETINT_DBG_S(XINT,ITYP,ISYMK,IBUK,JTYP,JSYMK,JBUK,
     &              KTYP,KSYMK,KBUK,LTYP,LSYMK,LBUK,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,1)),WORK(KPINT2D(1)),
     &              IREO,SCALE)
            END IF
         ELSE
*. ( U U ! B B ) or ( B B ! U U )

            IF(IBU.EQ.2) THEN
* ( U U ! B1 B2 ) = ( U U ! U2 U1)
               CALL KRA_SYM(KSYM,KBU,KSYMK,KBUK)
               CALL KRA_SYM(LSYM,LBU,LSYMK,LBUK)
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 4
               IREO(4) = 3
C?             WRITE(6,*) ' GETINT_DBG_S 7 '
               SCALE = 1.0D0
               CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &              KTYP,KSYMK,KBUK,LTYP,LSYMK,LBUK,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,1)),WORK(KPINT2D(1)),
     &              IREO,SCALE)
            ELSE IF (IBU.EQ.1) THEN
* ( B1 B2 ! U U) = (U2 U1 ! U U )
               CALL KRA_SYM(ISYM,IBU,ISYMK,IBUK)
               CALL KRA_SYM(JSYM,JBU,JSYMK,JBUK)
               IREO(1) = 2
               IREO(2) = 1
               IREO(3) = 3
               IREO(4) = 4
C?             WRITE(6,*) ' GETINT_DBG_S 8 '
               WRITE(6,*) ' IRI = ', IRI
               SCALE = 1.0D0
               WRITE(6,*) ' KINT2D(IRI,1)', KINT2D(IRI,1)
               CALL  GETINT_DBG_S(XINT,ITYP,ISYMK,IBUK,JTYP,JSYMK,JBUK,
     &              KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,1)),WORK(KPINT2D(1)),
     &              IREO,SCALE)
            END IF
*        ^ End of ( U U ! B B ) or ( B B ! U U ) switch
         END IF
*       ^ End of the cases for list 1,
      END IF
*     ^ End of the use of list 1
*
* List 2 :  (u b ! u u)
*
      IF((IBU.NE.JBU.AND.KBU.EQ.LBU).OR.
     &     (IBU.EQ.JBU.AND.KBU.NE.LBU)) THEN
*
* 8 cases
*         ( U B ! U U ) Simple
*         ( B U ! U U ) use c.c and Kramers symmetry
*         ( U B ! B B ) Use Kramers symmetry
*         ( B U ! B B ) use c.c
*
* The remaining four cases have different symmetries at particle 2,
* use particle 1-2 particle two symmetry to get these in place
*
         IF(IBU.EQ.2.AND.JBU.EQ.1.AND.KBU.EQ.2.AND.LBU.EQ.2) THEN
*. the simple one
            IREO(1) = 1
            IREO(2) = 2
            IREO(3) = 3
            IREO(4) = 4
            SCALE = 1.0D0
C?          WRITE(6,*) ' GETINT_DBG_S 9 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF(IBU.EQ.1.AND.JBU.EQ.2.AND.KBU.EQ.2.AND.LBU.EQ.2)THEN
* (B U ! U U ) => ( U B ! U U) by c.c
            IREO(1) = 2
            IREO(2) = 1
            IREO(3) = 4
            IREO(4) = 3
            IF(IRI.EQ.1) THEN
               SCALE = 1.0D0
            ELSE
               SCALE = -1.0D0
            END IF
C?          WRITE(6,*) ' GETINT_DBG_S 10 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF(IBU.EQ.2.AND.JBU.EQ.1.AND.KBU.EQ.1.AND.LBU.EQ.1) THEN
*. ( U  B1 ! B2 B3) => (U B1 ! U3 U2) by kramers
            IREO(1) = 1
            IREO(2) = 2
            IREO(3) = 4
            IREO(4) = 3
            CALL KRA_SYM(KSYM,KBU,KSYMK,KBUK)
            CALL KRA_SYM(LSYM,LBU,LSYMK,LBUK)
            SCALE = 1.0D0
C?          WRITE(6,*) ' GETINT_DBG_S 11 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYMK,KBUK,LTYP,LSYMK,LBUK,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF (IBU.EQ.1.AND.JBU.EQ.2.AND.KBU.EQ.1.AND.LBU.EQ.1) THEN
*. (B1 U ! B2 B3 ) => ( B1 U ! U3 U2 ) (kr) => (U B1 ! B3 B2) * (c.c)
            IREO(1) = 2
            IREO(2) = 1
            IREO(3) = 3
            IREO(4) = 4
            CALL KRA_SYM(KSYM,KBU,KSYMK,KBUK)
            CALL KRA_SYM(LSYM,LBU,LSYMK,LBUK)
            IF(IRI.EQ.1) THEN
               SCALE = 1.0D0
            ELSE
               SCALE = -1.0D0
            END IF
C?          WRITE(6,*) ' GETINT_DBG_S 12 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYMK,KBUK,LTYP,LSYMK,LBUK,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
*. and then the transpose cases
         ELSE IF(IBU.EQ.2.AND.JBU.EQ.2.AND.KBU.EQ.2.AND.LBU.EQ.1) THEN
* ( U U ! U B ) => ( U B ! U U ) by particle one -  two interchange
            IREO(1) = 3
            IREO(2) = 4
            IREO(3) = 1
            IREO(4) = 2
            SCALE = 1.0D0
C?          WRITE(6,*) ' GETINT_DBG_S 13 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF(IBU.EQ.2.AND.JBU.EQ.2.AND.KBU.EQ.1.AND.LBU.EQ.2) THEN
* (U U ! B U ) => (1-2) (B U ! U U ) => ( U B ! U U) by c.c
*
            IREO(1) = 4
            IREO(2) = 3
            IREO(3) = 2
            IREO(4) = 1
            IF(IRI.EQ.1) THEN
               SCALE = 1.0D0
            ELSE
               SCALE = -1.0D0
            END IF
C?          WRITE(6,*) ' GETINT_DBG_S 14 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF( IBU.EQ.1.AND.JBU.EQ.1.AND.KBU.EQ.2.AND.LBU.EQ.1) THEN
*. (B2 B3 ! U B1 ) =>(1-2) ( U  B1 ! B2 B3) => (U B1 ! U3 U2) by kramers
            CALL KRA_SYM(ISYM,IBU,ISYMK,IBUK)
            CALL KRA_SYM(JSYM,JBU,JSYMK,JBUK)
            SCALE = 1.0D0
            IREO(1) = 4
            IREO(2) = 3
            IREO(3) = 1
            IREO(4) = 2
C?          WRITE(6,*) ' GETINT_DBG_S 15 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYMK,IBUK,JTYP,JSYMK,JBUK,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF (IBU.EQ.1.AND.JBU.EQ.1.AND.KBU.EQ.1.AND.LBU.EQ.2) THEN
*. (B2 B3 ! B1 U) => (1-2) (B1 U ! B2 B3 ) => (Kr)
*  ( B1 U ! U3 U2 )  => (c.c) (U B1 ! U3 U2) *
            CALL KRA_SYM(ISYM,IBU,ISYMK,IBUK)
            CALL KRA_SYM(JSYM,JBU,JSYMK,JBUK)
            IREO(1) = 4
            IREO(2) = 3
            IREO(3) = 2
            IREO(4) = 1
            IF(IRI.EQ.1) THEN
               SCALE = 1.0D0
            ELSE
               SCALE = -1.0D0
            END IF
C?          WRITE(6,*) ' GETINT_DBG_S 16 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYMK,IBUK,JTYP,JSYMK,JBUK,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         END IF
*      ^ End of the eight choices for List 2
      END IF
*     ^ End of list 2
*
      RETURN
      END
