      SUBROUTINE H1STR(H,ISTRSM,ISTRTP,ISTRGP,
     &                  IOBSM,IOBTP,JOBSM,JOBTP,MXSXST,
     *                  IEXSTR,FACSTR,NEX,IH2TRM,NTESTG)
*
* Obtain one-electron string matrix elements sum(ij) h(ij) <J!e(ij)!I>
* with orbital i and j restricted to given orbital TS subsets
*
* from strings of symmetry ISTRSM and class ISTRTP
*
*
* =====
* Input
* =====
*

*   ISTRSM,ISTTP,ISTRGP :  symmetry, type and group  of input string
*   IOBSM,IOBTP : symmetry and type of orbital i
*   JOBSM,JOBTP : symmetry and type of orbital j
*   MXSXST      : Max number of single excitations for given string
*   H           : one-electron integrals for subsets of orbitals
*   IH2TRM      : ne 0 => include diagonal and one-electron
*                 operators from two body operator
*
* ======
* Output
* ======
*
*     NEX(ISTR)        : Number of excitations
*     IEXSTR(IEX,ISTR) : Number of excited string ( relative to offset )
*     FACSTR(IEX,ISTR)     : Phase factor of excitation.
*
*     Jeppe Olsen, March 1994, LUCIA version
*
      IMPLICIT REAL*8(A-H,O-Z)
*. a few include blocks
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1,
C    &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR,
C    &              KSBIA,KSBIB,KVEC3,KPNIJ,KIJKK

      DIMENSION H(*)
*. A bit of local scratch
      DIMENSION IACAR(2),ITPAR(2)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ***************** '
        WRITE(6,*) ' Output from H1STR '
        WRITE(6,*) ' ***************** '
      END IF
*
*. Number of orbitals and offsets
      NI   = NOBPTS(IOBTP,IOBSM)
      IOFF = IBTSOB(IOBTP,IOBSM)
*
      NJ   = NOBPTS(JOBTP,JOBSM)
      JOFF = IBTSOB(JOBTP,JOBSM)
      IF(NTEST.GE.1000) THEN
        write(6,*) ' iobtp,iobsm,jobtp,jobsm',
     &               iobtp,iobsm,jobtp,jobsm
        write(6,*) ' NI NJ IOFF JOFF ', NI,NJ,IOFF,JOFF
      END IF
*
*. Offset for input string
      ISTROF = IFRMR(WORK,KISTSO(ISTRGP),
     &         (ISTRSM-1)*NOCTYP(ISTRGP)+ISTRTP)
*. Number of input strings
      NIST = IFRMR(WORK,KNSTSO(ISTRGP),
     &         (ISTRSM-1)*NOCTYP(ISTRGP)+ISTRTP)
*. Group of output strings
      IACAR(1) = 2
      IACAR(2) = 1
*
      ITPAR(1) = IOBTP
      ITPAR(2) = JOBTP
      CALL NEWTYP(ISTRGP,ISTRTP,IACAR,ITPAR,2,KSTRGP,KSTRTP)
*. Symmetry of output strings
      CALL SYMCOM(3,0,IOBSM,JOBSM,IJSXSM)
      CALL SYMCOM(3,0,ISTRSM,IJSXSM,KSTRSM)
C     WRITE(6,*) ' ISTRSM,IJSXSM KSTRSM',ISTRSM,IJSXSM,KSTRSM
C     SYMCOM(ITASK,IOBJ,I1,I2,I12)
*. Offset for output strings
      KSTROF = IFRMR(WORK,KISTSO(ISTRGP),
     &         (KSTRSM-1)*NOCTYP(ISTRGP)+KSTRTP)
C     WRITE(6,*) ' off set for output strings ', KSTROF
*
*. Type of mappings
*
*. Nel => Nel - 1 electrons
      IF(ISTAC(ISTRGP,1).NE.0.AND.ISTAC(ISTRGP,2).NE.0) THEN
*. full list
        IM1FUL = 1
        LM1    = NACOB
      ELSE
*. Truncated list
        IM1FUL = 0
        LM1 = NELEC(ISTRGP)
      END IF
*. Nel -1 => Nel mapping, not so simple anymore
      NIEL = NELEC(ISTRGP)
      IF(ISTAC(ISTRGP+1,1).NE.0.AND.ISTAC(ISTRGP+1,2).NE.0) THEN
        IP1FUL = 1
        LP1 = NACOB
      ELSE
        IP1FUL = 0
        LP1 = -1
      END IF

      CALL  H1STRS(H,ISTROF,NIST,KSTROF,
     &              NI,IOFF,NJ,JOFF,
     &              WORK(KSTSTM(ISTRGP,1)),WORK(KSTSTM(ISTRGP,2)),
     &              LM1,IM1FUL,
     &              WORK(KSTSTM(ISTRGP+1,1)),WORK(KSTSTM(ISTRGP+1,2)),
     &              WORK(KSTSTMI(ISTRGP+1)), WORK(KSTSTMN(ISTRGP+1)),
     &              LP1,IP1FUL,
     &              MXSXST,IEXSTR,FACSTR,NEX,IH2TRM,WORK(KPNIJ),
     &              WORK(KIJKK),WORK(KOCSTR(ISTRGP)),NIEL,
     &              NTOOB,NTEST)

*
      RETURN
      END
