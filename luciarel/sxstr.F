      SUBROUTINE SXSTR(ISTRSM,ISTRTP,ISTRGP,
     &                  IOBSM,IOBTP,JOBSM,JOBTP,MXSXST,
     *                  ISXSTR,JSXSTR,IEXSTR,FACSTR,NEX,
     &                  IOFFDG,NTESTG)
*
* Obtain all excitations from a set of strings of given sym, type,
* and group that can be obtained by applying
* single excitations where each operator has a given sym and type
*
* IF IOFFDG .ne. 0 , only excitation a+i a j with i.ne.j are generated
*
*
* =====
* Input
* =====
*

*    ISTRSM,ISTTP,ISTRGP :  symmetry, type and group  of input string
*
*   IOBSM,IOBTP : symmetry and type of orbital i
*   JOBSM,JOBTP : symmetry and type of orbital j
*   MXSXST      : Max number of single excitations for given string
*
* ======
* Output
* ======
*     ISXSTR(IEX,ISTR) : I orbital indeces of SX
*     JSXSTR(IEX,ISTR) : J orbital indeces of SX
*     IEXSTR(IEX,ISTR) : Number of excited string ( relative to offset )
*     NEX(ISTR)        : Number of excitations
*     FACSTR(ISTR)     : Phase factor of excitation.
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
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPOBS,MXPR4T),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),
C    &              NOBPTS(6+MXPR4T,MXPOBS),IOBPTS(6+MXPR4T,MXPOBS),
C    &              ITOOBS(MXPOBS),ITPFTO(MXPORB),ISMFSO(MXPORB)
C     COMMON/STRINP/NSTTYP,MNRS1(MXPSTT),MXRS1(MXPSTT),
C    &              MNRS3(MXPSTT),MXRS3(MXPSTT),NELEC(MXPSTT),
C    &              IZORR(MXPSTT),IAZTP,IBZTP,IARTP(3,10),IBRTP(3,10),
C    &              NZSTTP,NRSTTP,ISTTP(MXPSTT)
C     COMMON/STRBAS/KSTINF,KOCSTR(MXPSTT),KNSTSO(MXPSTT),KISTSO(MXPSTT),
C    &              KSTSTM(MXPSTT,2),KZ(MXPSTT),
C    &              KSTREO(MXPSTT),KSTSM(MXPSTT),KSTCL(MXPSTT),
C    &              KEL1(MXPSTT),KEL3(MXPSTT),KACTP(MXPSTT),
C    &              KCOBSM,KNIFSJ,KIFSJ,KIFSJO,KSTSTX
C    &             ,KNDMAP(MXPSTT),KNUMAP(MXPSTT)
*. A bit of local scratch
      DIMENSION IACAR(2),ITPAR(2)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ***************** '
        WRITE(6,*) ' Output from SXSTR '
        WRITE(6,*) ' ***************** '
      END IF
*
*. Number of orbitals and offsets
      NI   = NOBPTS(IOBTP,IOBSM)
      IOFF = IBTSOB(IOBTP,IOBSM)
*
      NJ   = NOBPTS(JOBTP,JOBSM)
      JOFF = IBTSOB(JOBTP,JOBSM)
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
*. Offset for output strings
      KSTROF = IFRMR(WORK,KISTSO(ISTRGP),
     &         (KSTRSM-1)*NOCTYP(ISTRGP)+KSTRTP)
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
*. Nel -1 => Nel mapping, things are not so simple anymore so
      IF(ISTAC(ISTRGP+1,1).NE.0.AND.ISTAC(ISTRGP+1,2).NE.0 ) THEN
        IP1FUL = 1
        LP1 = NACOB
      ELSE
*. Contains only creations, compact with indeces
        IP1FUL = 0
        LP1 = -1
      END IF
C?    write(6,*) ' SXSTR : ISTRGP,IP1FUL LP1' ,
C?   &                     ISTRGP,IP1FUL,LP1
      CALL  SXSTRS(ISTROF,NIST,KSTROF,
     &              NI,IOFF,NJ,JOFF,
     &              WORK(KSTSTM(ISTRGP,1)),WORK(KSTSTM(ISTRGP,2)),
     &              LM1,IM1FUL,
     &              WORK(KSTSTM(ISTRGP+1,1)),WORK(KSTSTM(ISTRGP+1,2)),
     &              LP1,IP1FUL,
     &              WORK(KSTSTMI(ISTRGP+1)), WORK(KSTSTMN(ISTRGP+1)),
     &              MXSXST,ISXSTR,JSXSTR,IEXSTR,FACSTR,NEX,IOFFDG,
     &              NTEST)
*
      RETURN
      END
