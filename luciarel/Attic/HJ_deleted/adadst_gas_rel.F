      SUBROUTINE ADADST_GAS_REL(IOB,IOBSM,IOBTP,NIOB,
     &                      JOB,JOBSM,JOBTP,NJOB,
     &                      ISPGP,ISM,ITP,KMIN,KMAX,
     &                      I1,XI1S,LI1,NK,IEND,IFRST,KFRST,I12,K12,
     &                      SCLFAC)
*
*
* Obtain mappings
* a+IORB a+ JORB !KSTR> = +/-!ISTR>
* In the form
* I1(KSTR) =  ISTR if a+IORB a+ JORB !KSTR> = +/-!ISTR> , ISTR is in
* ISPGP,ISM,IGRP.
* (numbering relative to TS start)
*. Only excitations IOB. GE. JOB are included
* The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
* and IOB (JOB) is the first orbital to be used, and the number of orbitals
* to be checked is NIOB ( NJOB).
*
* Only orbital pairs IOB .gt. JOB are included
*
* The output is given in I1(KSTR,I,J) = I1 ((KSTR,(J-1)*NIOB + I)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* Only Kstrings with relative numbers from KMIN to KMAX are included
*
* If IEND .ne. 0 last string has been checked
*
* Jeppe Olsen , August of 95
*
* Note : Works for spinor symmetries also : Switching is based upon
*        Value of ITP
*
* ======
*. Input
* ======
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*./BIGGY
#include "wrkspc.inc"
*./ORBINP/
#include "orbinp.inc"
#include "strinp.inc"
#include "strbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
      COMMON/HIDSCRR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      COMMON/SSAVER/NELIS(4), NSTRKS(4)
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' ADADST_GAS in service '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        WRITE(6,*) ' IOB,IOBSM,IOBTP ', IOB,IOBSM,IOBTP
        WRITE(6,*) ' JOB,JOBSM,JOBTP ', JOB,JOBSM,JOBTP
      END IF
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) 'Problemo, ADADST '
C?      WRITE(6,*) ' SCLFAC = ',SCLFAC
C?    END IF
CTF
*   determine IUBSTR, at best somehow internally ...
*       IUBSTR = 1  :   adadst_gas for unbarred strings
*       IUBSTR = 2  :   adadst_gas for barred strings
* Barred or unbarred strings : Defined by ITP
        IF(MOD(ITP,2).EQ.0) THEN
*. Barred strings
          IUB = 2
        ELSE
*. Unbarred strings
          IUB = 1
        END IF
CTF
*
*. Internal affairs
*
      IF(I12.LE.4.AND.K12.LE.2) THEN
        KLLOC = KLOCSTR(K12)
        KLLZ = KLZ(I12)
        KLLREO = KLREO(I12)
      ELSE
        WRITE(6,*) ' ADST_GAS : Illegal value of I12 = ', I12
        Call Abend2( ' ADST_GAS : Illegal value of I12  ' )
      END IF

*
*. Supergroup and symmetry of K strings
*
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      CALL NEWTYP_REL(ISPGPABS,1,IOBTP,1,K1SPGPABS)
      CALL NEWTYP_REL(K1SPGPABS,1,JOBTP,1,KSPGPABS)
      CALL SYMCOM_REL(2,0,IOBSM,K1SM,ISM)
      CALL SYMCOM_REL(2,0,JOBSM,KSM,K1SM)
      IF(NTEST.GE.100) WRITE(6,*)
     & ' K1SM,K1SPGPABS,KSM,KSPGPABS : ',
     &   K1SM,K1SPGPABS,KSM,KSPGPABS
* In ADADS1_GAS we need : Occupation of KSTRINGS
*                         lexical => Actual order for I strings
* Generate if required
*
      IF(IFRST.NE.0) THEN
*.. Generate information about I strings
*. Arc weights for ISPGP
        NTEST2 = NTEST
        CALL WEIGHT_SPGP(WORK(KLLZ),NGAS,
     &                  NELFSPGP(1,ISPGPABS),
     &                  NOBPT,WORK(KLZSCR),NTEST2)
        NELI = NELFTP(ITP)
        NELIS(I12) = NELI
*. Reorder array for I strings
        CALL GETSTR_TOTSM_SPGP_REL(ITP,ISPGP,ISM,NELI,NSTRI,
     &                         WORK(KLLOC),NOCOB,1,
     &                         WORK(KLLZ),WORK(KLLREO))
      END IF
      NELK = NELIS(I12) - 2
      IF(KFRST.NE.0) THEN
*. Generate occupation of K STRINGS
*. Type of K : two electrons moved from I so
      KTP = ITP + 4
      KSPGPREL = KSPGPABS- IBSPGPFTP(KTP)+1
       CALL GETSTR_TOTSM_SPGP_REL(KTP,KSPGPREL,KSM,NELK,NSTRK,
     &                        WORK(KLLOC),NOCOB,0,
     &                        IDUM,IDUM)
       NSTRKS(K12) = NSTRK
      END IF
*
      NSTRK = NSTRKS(K12)
*
      IF(IUB.EQ.1) THEN
        IIOB = IOBPTS(IOBTP,IOBSM) + IOB - 1
        JJOB = IOBPTS(JOBTP,JOBSM) + JOB - 1
      ELSE
        IIOB = IOBPTS2(IOBTP,IOBSM) + IOB - 1
        JJOB = IOBPTS2(JOBTP,JOBSM) + JOB - 1
      END IF
*
C?    WRITE(6,*) '  IOBPTS(JOBTP,JOBSM) ',  IOBPTS(JOBTP,JOBSM)
      CALL ADADS1_GAS(NK,I1,XI1S,LI1,IIOB,NIOB,JJOB,NJOB,
     &          WORK(KLLOC),NELK,NSTRK,WORK(KLLREO),WORK(KLLZ),
     &          NOCOB,KMAX,KMIN,IEND,SCLFAC)
*
      RETURN
      END
