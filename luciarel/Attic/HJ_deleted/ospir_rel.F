      SUBROUTINE OSPIR_REL(NOSPIR,IOSPIR,PNTGRP,NIRREP,IPRNT,
     &                     NORB_SET)
*
* Number and symmetries of orbitals corresponding to a given shell
*
* =====
* Input
* =====
*
*   PNTGRP  : type of pointgroup
*         = 1 => D2h or a subgroup of D2H
*         = 2 => C inf v
*         = 3 => D inf h
*         = 4 => O 3
*         = 5 to 9 => different double groups :
*             9 :   C2H
*             8 :   C2
*             7 :   CS
*             6 :   CI
*             5 :   C1
*   NIRREP : Number of irreducible representations per point group
*   MXPIRR : Largest allowed number of shell irreps
*   MXPOBS : Largest allowed number of orbital symmetries
*
* ======
* Output
* ======
*
*   NOSPIR : Number of orbital symmetries per irrep
*   IOSPIR : Orbital symmetries corresponding to a given irrep
*
* Jeppe Olsen , Winter of 1991
*              Updated, July 97 : Double group symmetry added
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
      INTEGER PNTGRP
*. Double group info
#include "symm.inc"
*. Output
      DIMENSION NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR)
*
*. In general -except for double groups we have just a single
*  orbitals set
      NORB_SET = 1
      IF(PNTGRP.GE.5) NORB_SET = 2
      if (PNTGRP.eq.1) then
*=====
*.D2h
*=====
        NSMOB = 0
        DO 10 IRREP = 1, 8
          NOSPIR(IRREP) = 1
          IOSPIR(1,IRREP) = IRREP
   10   CONTINUE
      ELSE IF(PNTGRP.EQ.2) THEN
* =========
*. C inf V
* =========
* orbital symmetry is numbered as IML - MNMLOB + 1
        MNMLOB = -(NIRREP-1)
        DO 20 IRREP = 1, NIRREP
*.Irrep I contains orbitals with ML = -(IRREP-1),+(IRREP-1)
          IF(IRREP.EQ.1) THEN
            NOSPIR(IRREP) = 1
            IOSPIR(1,IRREP) = IRREP - 1 - MNMLOB + 1
          ELSE
            NOSPIR(IRREP) = 2
            IOSPIR(1,IRREP) = -(IRREP - 1) - MNMLOB + 1
            IOSPIR(2,IRREP) =  (IRREP - 1) - MNMLOB + 1
          END IF
   20   CONTINUE
      ELSE IF(PNTGRP.EQ.3) THEN
* ========
*. D inf H
* ========
* orbital symmetry is numbered as (PARITY-1) * NMLOB + IML - MNMLOB + 1
        MXMLOB =  NIRREP/2-1
        MNMLOB = -MXMLOB
        NMLOB =   NIRREP - 1
        IRREP = 0
        DO 35 IPARI = 1, 2
          IADD = (IPARI-1)*NMLOB
          DO 30 ML = 0,MXMLOB
            IRREP = IRREP + 1
            IF(ML.EQ.0) THEN
              NOSPIR(IRREP) = 1
              IOSPIR(1,IRREP) = IADD + ML - MNMLOB + 1
            ELSE
              NOSPIR(IRREP) = 2
              IOSPIR(1,IRREP) = IADD - ML - MNMLOB + 1
              IOSPIR(2,IRREP) = IADD + ML - MNMLOB + 1
            END IF
   30     CONTINUE
   35   CONTINUE

      ELSE IF(PNTGRP.EQ.4) THEN
* =====
*. O 3
* =====
* orbital symmetry is numbered as (PARITY-1) * NMLOB + IML - MNMLOB + 1
        MXMLOB =  NIRREP/2-1
        MNMLOB = -MXMLOB
        NMLOB =   NIRREP - 1
        DO 45 L = 0, NIRREP - 1
          IF(MOD(L,2).EQ.0) THEN
            IPARI = 1
          ELSE
            IPARI = 2
          END IF
          IADD = (IPARI-1)*NMLOB
          IRREP = L + 1
          NOSPIR(IRREP) = 2 * L + 1
          ICOMP = 0
          DO 40 ML = MNMLOB,MXMLOB
            ICOMP = ICOMP + 1
            IOSPIR(ICOMP,IRREP) = IADD + ML - MNMLOB + 1
   40     CONTINUE
   45   CONTINUE
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
* ========================
*. Double group symmetry
* ========================
*
* We are handling the transformation from orbitals to spinors
* by having two sets of orbitals. This turns out to be
* significant more convenient ( we hope ) than having
* a single set with double dimension
        NORB_SET = 2
        DO  IRREP = 1, NIRR_DG
          NOSPIR(IRREP) = 1
        END DO
*. IOSPIR become the symmetry of the spinors
        DO  IRREP = 1, NIRR_DG
*. The first set corresponds to alpha spin functions
          IOSPIR(1,IRREP) = IMOSF_SP(IRREP,1)
*. The second set corresponds to beta spin functions
          IOSPIR(2,IRREP) = IMOSF_SP(IRREP,2)
        END DO
*. The first NSMOB/2 symmetries are now the fermion irreps, whereas
*  the boson irreps are missing. Add these for completeness
        DO ISMOB = 1, NIRR_PN
         IOSPIR(1,NIRR_PN+ISMOB) = ISMOB
         IOSPIR(2,NIRR_PN+ISMOB) = ISMOB
        END DO
      ELSE
        WRITE(6,*) ' Sorry  PNTGRP out of range , PNTGRP = ', PNTGRP
        WRITE(6,*) ' OSPIR fatally wounded '
        Call Abend1( 5 )
      END IF
*
      NTEST = 00
      NTEST = MAX(IPRNT,NTEST)
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' OSPIR speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' Number of orbital sets ',NORB_SET
        WRITE(6,*) ' Number of orbitals per irrep '
*
        IF(.NOT.(PNTGRP.GE.5.AND.PNTGRP.LE.9)) THEN
*. Info for orbitals
          CALL IWRTMA(NOSPIR,1,NIRREP,1,NIRREP)
          WRITE(6,*) ' Orbital symmetries per irrep '
          DO 100 IRREP = 1, NIRREP
          CALL IWRTMA(IOSPIR(1,IRREP),1,NOSPIR(IRREP),1,NOSPIR(IRREP))
  100     CONTINUE
        ELSE
*. Info for spinors
          WRITE(6,*) ' IOSPIR for alpha spinors '
          WRITE(6,'(10I4)') (IOSPIR(1,I),I=1,NIRR_DG)
          WRITE(6,*) ' IOSPIR for beta spinors '
          WRITE(6,'(10I4)') (IOSPIR(2,I),I=1,NIRR_DG)
        END IF
      END IF
*
      RETURN
      END
