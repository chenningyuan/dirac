      SUBROUTINE ADAST_GAS_REL(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
     &                     I1,XI1S,NKSTR,IEND,
     &                     IFRST,KFRST,KACT,SCLFAC,IAC,IAB)
*
*
* Obtain creation or annihilation mapping
*
* IAC = 2 : Creation map
* a+IORB !KSTR> = +/-!ISTR>
*
* IAC = 1 : Annihilation map
* a IORB !KSTR> = +/-!ISTR>
*
* for orbitals of symmetry IOBSM and type IOBTP
* and Istrings defined by the NIGRP groups IGRP and symmetry ISPGPSM
*
* The results are given in the form
* I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR>
* (numbering relative to TS start)
* Above +/- is stored in XI1S
*
* if some nonvanishing excitations were found, KACT is set to 1,
* else it is zero
*
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version
*               October 96   : Improved version
*               September 97 : annihilation mappings added
*                              I groups defined by IGRP
*
* Timo Fleig, extensions for relativistic calculations,
*         - IAB added (1=alpha/unbarred, 2=beta/barred)
*     August 13, 2001
*
*
* ======
*. Input
* ======
*
*./BIGGY
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "csm.inc"
#include "lucinp.inc"
*. Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER ISMFGS(MXPNGAS)
      INTEGER MXVLI(MXPNGAS),MNVLI(MXPNGAS)
      INTEGER MXVLK(MXPNGAS),MNVLK(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      INTEGER KGRP(MXPNGAS)
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
      PARAMETER(LOFFI=8*8*8*8*8)
CTF   DIMENSION IOFFI(LOFFI)
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Will be stored as a matrix of dimension
* (NKSTR,*), Where NKSTR is the number of K-strings of
*  correct symmetry . Nk is provided by this routine.
*
      CALL QENTER('ADAST ')
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADAST_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) '  IOBTP IOBSM : ', IOBTP,IOBSM
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        WRITE(6,*) '  Symmetry of supergroup : ', ISPGPSM
        WRITE(6,*) ' SCLFAC = ', SCLFAC
*
        IF(IAC.EQ.1) THEN
          WRITE(6,*) ' Annihilation mapping '
        ELSE IF(IAC.EQ.2) THEN
          WRITE(6,*) ' Creation mapping '
        ELSE
          WRITE(6,*) ' Unknown IAC parameter in ADAST ',IAC
          Call Abend2( ' Unknown IAC parameter in ADAST ' )
        END IF
*
      END IF
*
      call memmar(IDUM,IDUM,'MARK  ',IDUM,'ADASTG')
*. A few preparations
      if (IAB.eq.1) then
        NORBTS= NOBPTS(IOBTP,IOBSM)
      else if (IAB.eq.2) then
        NORBTS= NOBPTS2(IOBTP,IOBSM)
      end if
      NORBT= NOBPT(IOBTP)
      IACGAS = IOBTP
*. First orbital of given GASpace
*. First orbital of given GASpace and Symmetry
      if (IAB.eq.1) then
        IBORBSPS = IOBPTS(IOBTP,IOBSM)
      else
        IBORBSPS = IOBPTS2(IOBTP,IOBSM)
      end if
      IBORBSP = IELSUM(NOBPT,IOBTP-1)+1
*
*====================================================
*. K strings : Supergroup, symmetry and distributions
*====================================================
*
      IF(IAC.EQ.1) THEN
       IDELTA = +1
      ELSE
       IDELTA = -1
      END IF
*. Is required mapping contained within current set of maps?
*. a:) Is active GASpace included in IGRP - must be
      IACGRP = 0
      DO JGRP = 1, NIGRP
        IF(IGSFGP(IGRP(JGRP)).EQ.IACGAS) IACGRP = JGRP
      END DO
*. Note : IACGRP is not the actual active group, 
*         it is the address of the active group in IGRP
      IF(IACGRP.EQ.0) THEN
        WRITE(6,*) ' ADAST in problems '
        WRITE(6,*) ' Active GASpace not included in IGRP '
        WRITE(6,*) ' Active GASpace : ', IACGAS
        Call Abend2
     &  ( ' ADAST : Active GASpace not included in IGRP ' )
      END IF
*. b:) active group in K strings
      NIEL = NELFGP(IGRP(IACGRP))
      NKEL = NIEL + IDELTA
      IF(NTEST.GE.1000) WRITE(6,*) ' NIEL and NKEL ',NIEL,NKEL
      NOBIACP1 = NOBPT(IACGAS) + 1
      IF(NKEL.EQ.-1.OR.NKEL.EQ.NOBIACP1) THEN
*. No strings with this number of elecs - be happy : No work
        NKSTR = 0
        KACT = 0
        KACGRP = 0
        GOTO 9999
      ELSE
*. Find group with NKEL electrons in IACGAS
        KACGRP = 0
        DO JGRP = IBGPSTR(IACGAS),IBGPSTR(IACGAS)+NGPSTR(IACGAS)-1
          IF(NELFGP(JGRP).EQ.NKEL) KACGRP = JGRP
        END DO
        IF(NTEST.GE.1000) WRITE(6,*) ' KACGRP = ',KACGRP
*. KACGRP is the Active group itself
        IF(KACGRP.EQ.0) THEN
          WRITE(6,*)' ADAST : cul de sac, active K group not found'
          WRITE(6,*)' GAS space and number of electrons ',
     &               IACGAS,NKEL
          Call Abend2
     &    ( ' ADAST : cul de sac, active K group not found' )
        END IF
      END IF
*. Okay active K group was found and is nontrivial
      CALL SYMCOM_REL(2,0,IOBSM,KSM,ISPGPSM)
*. The K supergroup
      CALL ICOPVE(IGRP,KGRP,NIGRP)
      KGRP(IACGRP) = KACGRP
*. Number of strings and symmetry distributions of K strings
      call nst_spgrp_dbg(NIGRP,KGRP,KSM,NSMST,NKSTR,NKDIST,IAB)
CTF   CALL NST_SPGRP(NIGRP,KGRP,KSM,WORK(KNSTSGP(1)),
CTF  &               NSMST,NKSTR,NKDIST)
      IF(NTEST.GE.1000) WRITE(6,*)
     & ' KSM, NKSTR : ', KSM, NKSTR
      IF(NKSTR.EQ.0) GOTO 9999
*. Last active space in K strings and number of strings 
*  per group and sym
      NGASL = 1
      DO JGRP = 1, NIGRP
        IF(NELFGP(KGRP(JGRP)).GT.0) NGASL = JGRP
        if (IAB.eq.1) then
          call icopve2(WORK(KNSTSGP(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,JGRP))
          call icopve2(WORK(KISTSGP(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &                 IISTSGP(1,JGRP))
        else if (IAB.eq.2) then
          call icopve2(WORK(KNSTSGP2(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,JGRP))
          call icopve2(WORK(KISTSGP2(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &                 IISTSGP(1,JGRP))
        end if
      END DO
*. MIN/MAX for Kstrings
      call minmax_for_sym_dist_dbg(NIGRP,KGRP,MNVLK,MXVLK,
     &                             NKDIST_TOT,IAB)
      IF(NTEST.GE.1000) THEN
        write(6,*) 'MNVLK and MXVLK '
        CALL IWRTMA(MNVLK,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MXVLK,1,NIGRP,1,NIGRP)
      END IF
*. (NKDIST_TOT is number of distributions, all symmetries )
* ==============
*. I Strings
* ==============
*. Generate symmetry distributions of I strings with given symmetry
      call memmar(KIOFFI,LOFFI,'ADDL  ',1,'IOFFI ')
      CALL TS_SYM_PNT2_REL(IGRP,NIGRP,MXVLI,MNVLI,ISPGPSM,
     &                     WORK(KIOFFI),LOFFI,IAB)
*. Offset and dimension for active group in I strings
      if (IAB.eq.1) then
        CALL ICOPVE2(WORK(KISTSGP(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               IACIST)
        CALL ICOPVE2(WORK(KNSTSGP(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               NACIST)
      else if (IAB.eq.2) then
        CALL ICOPVE2(WORK(KISTSGP2(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               IACIST)
        CALL ICOPVE2(WORK(KNSTSGP2(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               NACIST)
      end if
*. Last entry in IGRP with a nonvanisking number of strings
      NIGASL = 1
      DO JGRP = 1, NIGRP
        IF(NELFGP(IGRP(JGRP)).GT.0) NIGASL = JGRP
      END DO
C?    WRITE(6,*) ' NIGASL = ', NIGASL

*. Number of electrons before active space
      NELB = 0
      DO JGRP = 1, IACGRP-1
        NELB = NELB + NELFGP(IGRP(JGRP))
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) ' NELB = ', NELB
*
      ZERO =0.0D0
      IZERO = 0    
COLD  CALL SETVEC(XI1S,ZERO,NORBTS*NKSTR)
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
      DO IGAS = 1, NIGRP
        ISMFGS(IGAS) = 1
      END DO
 1000 CONTINUE
*. Next distribution
        CALL NEXT_SYM_DISTR_REL(NGASL,MNVLK,MXVLK,ISMFGS,KSM,KFIRST,
     &                          NONEW)
        IF(NTEST.GE.1000) THEN
          write(6,*) ' Symmetry distribution '
          call iwrtma(ISMFGS,1,NIGRP,1,NIGRP)
        END IF
        IF(NONEW.EQ.1) GOTO 9999
* NONEW is not equal to 1, treat this distribution:
        KFIRST = 0
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
          NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IACGRP)
C       CALL SYMCOM_REL(3,1,IOBSM,ISMFGS(IOBTP),IACSM)
        CALL SYMCOM_REL(3,1,IOBSM,ISMFGS(IACGRP),IACSM)
        ISMFGS(IACGRP) = IACSM
        IBSTRINI = IOFF_SYM_DIST(ISMFGS,NIGASL,WORK(KIOFFI),MXVLI,MNVLI)
        ISMFGS(IACGRP) = ISAVE
C?      WRITE(6,*) ' IOBSM, KACSM IACSM ',
C?   &               IOBSM,ISAVE, IACSM  
C?      WRITE(6,*) 'IBSTRINI ', IBSTRINI
*. Number of strings before active GAS space
        NSTB = 1
C       DO IGAS = 1, IOBTP-1
        DO IGAS = 1, IACGRP-1
          NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number of strings After active GAS space
        NSTA = 1
C       DO IGAS =  IOBTP +1, NIGRP
        DO IGAS =  IACGRP+1, NIGRP
          NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number and offset for active group
        NIAC = NACIST(IACSM)
        IIAC = IACIST(IACSM)
*
        NKAC = NNSTSGP(ISMFGS(IACGRP),IACGRP)
        IKAC = IISTSGP(ISMFGS(IACGRP),IACGRP)
*. I and K strings of given symmetry distribution
        NISD = NSTB*NIAC*NSTA
        NKSD = NSTB*NKAC*NSTA
        IF(NTEST.GE.1000) THEN
          write(6,*) ' nstb nsta niac nkac ',
     &                 nstb,nsta,niac,nkac
        END IF
*. Obtain annihilation/creation mapping for all strings of this type
*. Are group mappings in expanded or compact form
        IF(IAC.EQ.1.AND.ISTAC(KACGRP,2).EQ.0) THEN
          IEC = 2
          LROW_IN = NKEL
        ELSE
          IEC = 1
C         LROW_IN = NORBTS
          LROW_IN = NORBT
        END IF
        NKACT = NSTFGP(KACGRP)
*
CTF     MXAADST = MAX(MXAADST,NKSTR*NORBTS)
C     COMMON/COMJEP/MXACJ,MXACIJ,MXAADST
        if (NSTA*NSTB*NIAC*NKAC.NE.0) then
          if (IAB.eq.1) then
            CALL ADAST_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,
     &                       WORK(KSTSTM(KACGRP,1)),
     &                       WORK(KSTSTM(KACGRP,2)),
     &                       IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &                       NKSTR,KBSTRIN,NELB,
     &                       I1,XI1S,SCLFAC,IAC,
     &                       LROW_IN,IEC)
          else if (IAB.eq.2) then
            CALL ADAST_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,
     &                       WORK(KSTSTM2(KACGRP,1)),
     &                       WORK(KSTSTM2(KACGRP,2)),
     &                       IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &                       NKSTR,KBSTRIN,NELB,
     &                       I1,XI1S,SCLFAC,IAC,
     &                       LROW_IN,IEC)
          end if
        end if
        KSTRBS = KSTRBS + NKSD
C       IF(NGASL-1.GT.0) GOTO 1000
        GOTO 1000
 1001 CONTINUE
*
 9999 call memmar(IDUM,IDUM,'FLUSM ',IDUM,'ADASTG')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADAST_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Total number of K strings ', NKSTR
        IF(NKSTR.NE.0) THEN
          DO IORB = IBORBSPS,IBORBSPS + NORBTS  - 1
            IORBR = IORB-IBORBSPS +1
            WRITE(6,*) ' Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(  I1((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
            CALL WRTMAT(XI1S((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
          END DO
        END IF
      END IF
*
      CALL QEXIT('ADAST ')
      RETURN
      END
