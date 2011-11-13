!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE AATYP(NAATP,ITYP,JTYP,LOCC,ROCC,NGAS)
*
* Obtain types of I,J so
* <L!aI aJ !R> is nonvanishing
* only combinations with type(I) .ge. type(J) * are included
*
*
#include "ipoist8.inc"
      INTEGER LOCC(NGAS),ROCC(NGAS)
*. First ensure that the types differs by annihilation of two elecs
      ICREA = 0
      IANNI = 0
      DO IGAS = 1, NGAS
        IF(LOCC(IGAS).LT.ROCC(IGAS)) THEN
         IANNI = IANNI + ROCC(IGAS)-LOCC(IGAS)
        ELSE IF (LOCC(IGAS).GT.ROCC(IGAS) ) THEN
         ICREA = ICREA + LOCC(IGAS) - ROCC(IGAS)
        END IF
      END DO
*
      IF(ICREA.EQ.0 .AND. IANNI.EQ.2) THEN
       NAATP = 1
      ELSE
       NAATP = 0
      END IF
*
      IANNI = 0
      IF(NAATP.EQ.1) THEN
        DO IGAS = 1, NGAS
          IF(LOCC(IGAS)-ROCC(IGAS) .EQ. -2 ) THEN
            ITYP = IGAS
            JTYP = IGAS
          ELSE IF (LOCC(IGAS)-ROCC(IGAS) .EQ. -1 ) THEN
            IF(IANNI.EQ.1) THEN
             ITYP = IGAS
             IANNI = 2
            ELSE IF( IANNI.EQ. 0 ) THEN
             JTYP = IGAS
             IANNI = 1
            END IF
          END IF
        END DO
      END IF
*
#if defined LUCI_DEBUG
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Double Annihilation connecting  '
        WRITE(6,*) ' Left string '
        CALL IWRTMA(LOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Right string '
        CALL IWRTMA(ROCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        IF(NAATP.EQ.0) THEN
          WRITE(6,*) ' No connecting double annihilation'
        ELSE IF(NAATP.EQ.1) THEN
          WRITE(6,*) ' Types of connecting annihilation ops', ITYP,JTYP
        END IF
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
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
#include "ipoist8.inc"
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
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' ADADST_GAS in service '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        WRITE(6,*) ' IOB,IOBSM,IOBTP ', IOB,IOBSM,IOBTP
        WRITE(6,*) ' JOB,JOBSM,JOBTP ', JOB,JOBSM,JOBTP
      END IF
#endif
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
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) WRITE(6,*)
     & ' K1SM,K1SPGPABS,KSM,KSPGPABS : ',
     &   K1SM,K1SPGPABS,KSM,KSPGPABS
#endif
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
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ADAST_GAS_REL(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
     &                     I1,XI1S,NKSTR,IEND,
     &                     IFRST,KFRST,KACT,SCLFAC,IAC,ISMDST,IAB)
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
*               Aug  03      : ISMDST from the outside added 
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
#include "ipoist8.inc"
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
*. Scratch through the argument list
      INTEGER ISMDST(*)
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
C     CALL QENTER('ADAST ')
*
      NONEW = -1
      NTEST = 0000
#if defined LUCI_DEBUG
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
#endif
*
C!    call memmar(KDUM,IDUM,'MARK  ',IDUM,'ADASTG')
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
#if defined LUCI_DEBUG
      IF(NTEST.GE.1000) WRITE(6,*) ' NIEL and NKEL ',NIEL,NKEL
#endif
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
#if defined LUCI_DEBUG
        IF(NTEST.GE.1000) WRITE(6,*) ' KACGRP = ',KACGRP
#endif
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
#if defined LUCI_DEBUG
      IF(NTEST.GE.1000) WRITE(6,*)
     & ' KSM, NKSTR : ', KSM, NKSTR
#endif
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
#if defined LUCI_DEBUG
      IF(NTEST.GE.1000) THEN
        write(6,*) 'MNVLK and MXVLK '
        CALL IWRTMA(MNVLK,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MXVLK,1,NIGRP,1,NIGRP)
      END IF
#endif
*. (NKDIST_TOT is number of distributions, all symmetries )
* ==============
*. I Strings
* ==============
*. Generate symmetry distributions of I strings with given symmetry
C!!   call memmar(KIOFFI,LOFFI,'ADDL  ',1,'IOFFI ')
C     CALL TS_SYM_PNT2_REL(IGRP,NIGRP,MXVLI,MNVLI,ISPGPSM,
C    &                     WORK(KIOFFI),LOFFI,IAB)
      CALL TS_SYM_PNT2_REL(IGRP,NIGRP,MXVLI,MNVLI,ISPGPSM,
     &                     ISMDST,LOFFI,IAB)
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
#if defined LUCI_DEBUG
      IF(NTEST.GE.1000) WRITE(6,*) ' NELB = ', NELB
#endif
*
      CALL IZERO(I1,NORBTS*NKSTR)
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
#if defined LUCI_DEBUG
        IF(NTEST.GE.1000) THEN
          write(6,*) ' Symmetry distribution '
          call iwrtma(ISMFGS,1,NIGRP,1,NIGRP)
        END IF
#endif
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
        IBSTRINI = IOFF_SYM_DIST(ISMFGS,NIGASL,ISMDST,MXVLI,MNVLI)
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
#if defined LUCI_DEBUG
        IF(NTEST.GE.1000) THEN
          write(6,*) ' nstb nsta niac nkac ',
     &                 nstb,nsta,niac,nkac
        END IF
#endif
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
 9999 CONTINUE
C!    call memmar(KDUM,IDUM,'FLUSM ',IDUM,'ADASTG')
*
#if defined LUCI_DEBUG
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
#endif
*
C     CALL QEXIT('ADAST ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ADD_STR_GROUP_REL(NSTADD,IOFADD,ISTADD,NSTB,NSTA,
     &                         ISTRING,IELOF,NELADD,NELTOT)
*
* Part of assembling strings in individual types to
* super group of strings
*
*. Copying strings belonging to a given type to supergroup of strings
*
* Jeppe Olsen, for once improving performance of LUCIA
*
*.Input
* =====
* NSTADD : Number of strings to be added
* IOFADD : First string to be added
* ISTADD : Strings to be added
* NSTB   : Number of strings belonging to lower gasspaces
* NSTA   : Number of strings belonging to higher gasspaces
* ISTRING: Supergroup of strings under construction
* IELOF  : Place of first electron to be added
* NELADD : Number of electrons to be added
* NELTOT : Total number of electrons
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION ISTADD(NELADD,*)
*. Input and output
      DIMENSION ISTRING(NELTOT,*)
*
      IF(NSTA.GT.1) THEN
        DO IISTR = 1,NSTADD
*. Address of A(1,IISTR,1)
*. A(I(after),Igas,I(before))
          IOFF1 = (IISTR-1)*NSTA + 1
          DO ISTB = 1, NSTB
*. Address of A(1,IISTR,ISTB)
            IOFF2 = IOFF1 + (ISTB-1)*NSTADD*NSTA
            DO ISTA = 1, NSTA
              DO JEL = 1, NELADD
                ISTRING(IELOF-1+JEL,IOFF2-1+ISTA)
     &        = ISTADD(JEL,IOFADD-1+IISTR)
              END DO
            END DO
          END DO
        END DO
      ELSE IF (NSTA .EQ. 1 ) THEN
*. Address of A(1,IISTR,1)
*. A(I(after),Igas,I(before))
        DO ISTB = 1, NSTB
          IOFF0 = (ISTB-1)*NSTADD
          DO IISTR = 1,NSTADD
*. Address of A(1,IISTR,ISTB)
            IOFF2 = IISTR  + IOFF0
            DO JEL = 1, NELADD
              ISTRING(IELOF-1+JEL,IOFF2)
     &      = ISTADD(JEL,IOFADD-1+IISTR)
            END DO
          END DO
        END DO
      END IF
*
#if defined LUCI_DEBUG
      NTEST = 0000
      if (NTEST.ge.1000) then
        write(6,*) 'ADD_STR_GROUP speaking: Input string array'
        call iwrtma(ISTADD,NELADD,IOFADD-1+NSTADD,
     &                     NELADD,IOFADD-1+NSTADD)
      end if
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ADS1_GAS(NK,I1,XI1S,LI1,IORB,LORB,
     &                KSTR,NKEL,NKSTR,IREO,IZ,
     &                NOCOB,KMAX,KMIN,IEND,KACT,SCLFAC)
*
* Obtain I1(KSTR) = +/- A+ IORB !KSTR>
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* IORB : First orbital to be added
* LORB : Number of orbitals to be added : IORB to IORB-1+LORB
*        are used. They must all be in the same TS group
* KMAX : Largest allowed relative number for K strings
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR,JORB) : ne. 0 => a + JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR,JORB) : above +/-
*          : eq. 0    a + JORB !KSTR> = 0
* Offset is KMIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*.Input
      INTEGER KSTR(NKEL,NKSTR)
      INTEGER IREO(*), IZ(NOCOB,*)
*.Output
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ================= '
       WRITE(6,*) ' ADS1_GAS speaking '
       WRITE(6,*) ' ================= '
       WRITE(6,*) ' IORB,LORB ', IORB,LORB
      END IF
#endif
C?    WRITE(6,*) ' IZ matrix in ADS1 '
C?    CALL IWRTMA(IZ,NOCOB,NKEL+1,NOCOB,NKEL+1)
*
      KACT = 0
      IORBMIN = IORB
      IORBMAX = IORB + LORB - 1
*
      KEND = MIN(NKSTR,KMAX)
      NK = MAX(0,KEND-KMIN+1)
      IF(KEND.GE.NKSTR) THEN
        IEND = 1
      ELSE
        IEND = 0
      END IF
*. Start by zeroing
      DO IIORB = 1, LORB
C       CALL SETVEC(XI1S(1,IIORB),ZERO,NK)
        CALL IZERO(I1(1,IIORB),NK)
      END DO
*
      DO KKSTR = KMIN,KEND
C?    write(6,'(a,I12)') 'KKSTR    ',KKSTR
C?    WRITE(6,*) ' Occ of K ', (KSTR(JEL,KKSTR),JEL=1,NKEL)
#if defined LUCI_DEBUG
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Occupation of string ', KKSTR
         CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
       END IF
#endif
*
       SGN = -1.0D0*SCLFAC
       DO KEL = 1, NKEL+1
         SGN = - SGN
         IF(KEL.LE.NKEL) THEN
           KORB = KSTR(KEL,KKSTR)
         ELSE
           KORB = IORB+LORB
         END IF
*
         IF(KORB.GT.IORBMIN) THEN
*. Orbitals can be added before electron KEL
           IORB2 = MIN(KORB-1,IORBMAX)
           IF(KEL.EQ.1) THEN
             IORB1 = IORBMIN
           ELSE
             IORB1 = MAX(IORBMIN,KSTR(KEL-1,KKSTR)+1)
           END IF
           IF(IORB1.GT.IORBMAX) GOTO 2912
           IF(IORB2.GE.IORB1) THEN
*
*. common piece of lexical ordering
C            IF(KEL.EQ.1) THEN
C              IREVLEX = 1
C              DO JEL = 1, NKEL
C                 IREVLEX = IREVLEX + IZ(KSTR(JEL,KKSTR),JEL+1)
C              END DO
C            ELSE
C                 IREVLEX = IREVLEX + IZ(KSTR(KEL-1,KKSTR),KEL-1)
C    &                              - IZ(KSTR(KEL-1,KKSTR),KEL)
C            END IF
C          IF(IORB2.GE.IORB1) THEN
*
             IREVLEX = 1
             DO JEL = 1, KEL-1
               IREVLEX = IREVLEX + IZ(KSTR(JEL,KKSTR),JEL)
             END DO
             DO JEL = KEL, NKEL
               IREVLEX = IREVLEX + IZ(KSTR(JEL,KKSTR),JEL+1)
             END DO
#if defined LUCI_DEBUG
             IF(NTEST.GE.1000) WRITE(6,*)
     &       'KEL IORB1,IORB2,IREVLEX ',KEL,IORB1,IORB2,IREVLEX
#endif
             DO IIORB = IORB1, IORB2
               XI1S(KKSTR-KMIN+1,IIORB-IORBMIN+1) = SGN
               I1(KKSTR-KMIN+1,IIORB-IORBMIN+1) =
     &         IREO(IREVLEX  + IZ(IIORB,KEL))
*
               KACT = 1
*
               IF( I1(KKSTR-KMIN+1,IIORB-IORBMIN+1).LE.0) THEN
                 WRITE(6,*)
     &           ' Problems with map, stringnumber .le.0'
                 write(6,*) 'KKSTR         ',KKSTR
                 WRITE(6,*) ' Occupation of K string '
                 CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
                 WRITE(6,*) ' Occupation of previous strings '
                 CALL IWRTMA(KSTR,NKEL,KKSTR-1,NKEL,KKSTR-1)
                 write(6,*) 'KMIN          ',KMIN
                 write(6,*) 'IIORB         ',IIORB
                 write(6,*) 'IORBMIN       ',IORBMIN
                 write(6,*) 'KEL           ',KEL
                 write(6,*) 'IREVLEX       ',IREVLEX
                 write(6,*) 'IZ(IIORB,KEL) ',IZ(IIORB,KEL)
                 write(6,*) 'IREO(IREVLEX  + IZ(IIORB,KEL)) =',
     &                       IREO(IREVLEX  + IZ(IIORB,KEL))
                 Call Abend
     &           ( ' Problems with map, stringnumber .le.0')
               END IF
             END DO
           END IF
         END IF
       END DO
 2912  CONTINUE
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADS1_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO 200 IIORB = IORB,IORB+LORB-1
            IIORBR = IIORB-IORB+1
            WRITE(6,*) ' Info for orbital ', IIORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(I1(1,IIORBR),1,NK,1,NK)
            CALL WRTMAT(XI1S(1,IIORBR),1,NK,1,NK)
  200     CONTINUE
        END IF
      END IF
#endif
*
C     WRITE(6,*) 'JEPPE forced me to stop in ADS1   '
C     STOP'JEPPE forced me to stop in ADS1   '
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ADST_GAS(IORB,NORB,ISPGP,ISM,ITP,KMIN,KMAX,
     &                    I1,XI1S,LI1,NK,IEND,IFRST,KFRST,I12,K12,
     &                    KACT,SCLFAC)
*
*
* Obtain mappings
* a+JORB !KSTR> = +/-!ISTR> for orbitals IORB - IORB+NORB-1
* and KSTR in the range KMIN-KMAX
*. All orbitals are assumed to belong to the same TS class
*. and are type-symmetry ordered
*
* The results are given in the form
* I1(KSTR,JORB) =  ISTR if A+JORB !KSTR> = +/-!ISTR> , ISTR is in
* ISPGP,ISM,ITP .
*
* if some nonvanishing excitations were found, KACT is set to 1,
* else it is zero

* (numbering relative to TS start)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* If all K strings have been searched IEND is set to 1
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version
*
* ======
*. Input
* ======
*
*./BIGGY
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
*. Local scratch
      COMMON/HIDSCRR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*. save internally
      COMMON/SSAVER/NELIS(4), NSTRKS(4)
*
* =======
*. Output
* =======
*
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*
C     CALL QENTER('ADST  ')
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' ADST_GAS in service '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' ISM ISPGP ITP',ISM,ISPGP,ITP
      END IF
#endif
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) ' Problemo, ADST_GAS'
C?      WRITE(6,*) 'SCLFAC = ',SCLFAC
C?    END IF
*. Internal affairs
      IF(I12.LE.4.AND.K12.LE.2 ) THEN
        KLLOC = KLOCSTR(K12)
        KLLZ = KLZ(I12)
        KLLREO = KLREO(I12)
      ELSE
        WRITE(6,*)
     &   ' ADST_GAS : Illegal value of I12,K12 = ', I12,K12
        Call Abend2( ' ADST_GAS : Illegal value of I12,K12  ' )
      END IF
CTF
*   determine IUBSTR, at best somehow internally ...
*       IUBSTR = 1  :   adst_gas for unbarred strings
*       IUBSTR = 2  :   adst_gas for barred strings
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
*. Supergroup and symmetry of K strings
*
      if (IUB.eq.1) then
        IORBTP = ITPFTO(IORB)
        IORBSM = ISMFTO(IORB)
      else if (IUB.eq.2) then
        IORBTP = ITPFTO2(IORB)
        IORBSM = ISMFTO2(IORB)
      else
        write(6,*) 'IUBSTR out of range in adst_gas ...  ',IUBSTR
        Call Abend('Quitting')
      end if
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) WRITE(6,*)
     &' ADST : IORB IORBTP IORBSM : ', IORB,IORBTP,IORBSM
#endif
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      CALL NEWTYP_REL(ISPGPABS,1,IORBTP,1,KSPGPABS)
      KTP = ITP + 2
      KSPGPREL = KSPGPABS - IBSPGPFTP(KTP) + 1
*.
      CALL SYMCOM_REL(2,0,IORBSM,KSM,ISM)
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) WRITE(6,*)
     & ' KSM and KSPGPABS : ', KSM,KSPGPABS
#endif
*. In order to construct mapping we need to generate
*.      Reorder array : lexical order => actual order for I strings
*.      Occupation of K strings
*
       IF(IFRST.NE.0) THEN
*. Arc weights for ISPGP
          NTEST2 = NTEST
          CALL WEIGHT_SPGP(WORK(KLLZ),NGAS,
     &                     NELFSPGP(1,ISPGPABS),
     &                     NOBPT,WORK(KLZSCR),NTEST2)
*. Reorder array for I strings
         NELI = NELFTP(ITP)
         NELIS(I12) = NELI
C?       CALL QENTER('GENSTR')
         CALL GETSTR_TOTSM_SPGP_REL(ITP,ISPGP,ISM,NELI,NSTRI,
     &                          WORK(KLLOC),NOCOB,
     &                          1,WORK(KLLZ),WORK(KLLREO))
C?       CALL QEXIT('GENSTR')
       END IF
       NELK = NELIS(I12) - 1
C?     WRITE(6,*) ' ADST : NELK, NELI ', NELK,NELI
       IF(KFRST.NE.0) THEN
*. Occupation of K strings
CT       CALL QENTER('GENSTR')
         CALL GETSTR_TOTSM_SPGP_REL(KTP,KSPGPREL,KSM,NELK,NSTRK,
     &                          WORK(KLLOC),NOCOB,
     &                          0,ISUM,IDUM)
CT       CALL QEXIT('GENSTR')
         NSTRKS(K12) = NSTRK
      END IF
*. a+i !kstr> = +/- !istr>
      NSTRK = NSTRKS(K12)
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) WRITE(6,*) ' NSTRK = ', NSTRK
#endif
      CALL ADS1_GAS(NK,I1,XI1S,LI1,IORB,NORB,
     &          WORK(KLLOC),NELK,NSTRK,WORK(KLLREO),WORK(KLLZ),
     &          NOCOB,KMAX,KMIN,IEND,KACT,SCLFAC)
*
C     SAVE NELK , NSTRK
C     CALL QEXIT('ADST  ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ADSTN_GAS_REL(IOBSM,IOBTP,ISPGP,ISPGPSM,ISPGPTP,
     &                    I1,XI1S,NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC)
*
*
* Obtain mappings
* a+IORB !KSTR> = +/-!ISTR> for orbitals of symmetry IOBSM and type IOBTP
* and I strings belonging to supergroup ISPGP wih symmetry ISPGPSM
* and type ISPGPTP(=1=>alpha,=2=>beta)
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
*
* ======
*. Input
* ======
*
*./BIGGY
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
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
*. Local scratch
      INTEGER NELFGS(MXPNGAS),ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
      PARAMETER(LOFFI=8*8*8*8*8*8*8)
      PARAMETER(MXLNGAS=7)
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Will be stored as an matrix of dimension
* (NKSTR,*), Where NKSTR is the number of K-strings of
*  correct symmetry . Nk is provided by this routine.
*
      DIMENSION IOFFI(LOFFI)
      call memmar(KLOFFI,LOFFI,'ADDL',1,'LLOFFI')
C     CALL QENTER('ADSTN ')
*
      IF(NGAS.GT.MXLNGAS) THEN
        WRITE(6,*) ' Ad hoc programming in ADSTN (IOFFI)'
        WRITE(6,*) ' Must be changed - or redimensioned '
        Call Abend2( 'ADST : IOFFI problem ' )
      END IF
*
      NTEST = 00
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADSTN_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) '  IOBTP IOBSM : ', IOBTP,IOBSM
        WRITE(6,*) '  ISPGP ISPGPSM ISPGPTP :  ',
     &                ISPGP,ISPGPSM,ISPGPTP
      END IF
#endif
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) ' Problemo : ADSTN_GAS'
C?      WRITE(6,*) ' SCLFAC .ne. 1 '
C?    END IF
*
CTF
*   determine IUBSTR, at best somehow internally ...
*       IUBSTR = 1  :   adst_gas for unbarred strings
*       IUBSTR = 2  :   adst_gas for barred strings
* Barred or unbarred strings : Defined by ITP
        IF(MOD(ISPGPTP,2).EQ.0) THEN
*. Barred strings
          IUB = 2
        ELSE
*. Unbarred strings
          IUB = 1
        END IF
CTF
*
*. Supergroup and symmetry of K strings
*
      ISPGRPABS = IBSPGPFTP(ISPGPTP)-1+ISPGP
      CALL NEWTYP_REL(ISPGRPABS,1,IOBTP,1,KSPGRPABS)
      CALL SYMCOM_REL(2,0,IOBSM,KSM,ISPGPSM)
CSTART
      IF(IUB.EQ.1) THEN
         NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
      ELSE
         NKSTR = NSTFSMSPGP2(KSM,KSPGRPABS)
      END IF
CEND
#if defined LUCI_DEBUG
      IF(NTEST.GE.200) WRITE(6,*)
     & ' KSM, KSPGPRABS, NKSTR : ', KSM,KSPGRPABS, NKSTR
#endif
      IF(NKSTR.EQ.0) GOTO 9999
*
      IF(IUB.EQ.1) THEN
        NORBTS= NOBPTS(IOBTP,IOBSM)
      ELSE
        NORBTS= NOBPTS2(IOBTP,IOBSM)
      END IF
      CALL DZERO(XI1S,NORBTS*NKSTR)
      CALL IZERO(I1,NORBTS*NKSTR)
*
*. First orbital of given GASSpace
       IBORBSP = IELSUM(NOBPT,IOBTP-1)+1
*. First orbital of given GASSPace and Symmetry
       IF(IUB.EQ.1) THEN
         IBORBSPS = IOBPTS(IOBTP,IOBSM)
       ELSE
         IBORBSPS = IOBPTS2(IOBTP,IOBSM)
       END IF
*
*. Information about I strings
* =============================
*
*. structure of group of strings defining I strings
      NGASL = 1
      DO IGAS = 1, NGAS
       ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPABS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
*. Number of electrons before active type
      NELB = 0
      DO IGAS = 1, IOBTP -1
        NELB = NELB + NELFGS(IGAS)
      END DO
*. Number of electrons in active space
      NACGSOB = NOBPT(IOBTP)

*. Number of strings per symmetry for each symmetry
      DO IGAS = 1, NGAS
        IF(IUB.EQ.1) THEN
          CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,IGAS))
        ELSE
          CALL ICOPVE2(WORK(KNSTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,IGAS))
        END IF
      END DO
*. Offset and dimension for active group in I strings
      IF(IUB.EQ.1) THEN
       CALL ICOPVE2(WORK(KISTSGP(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               IACIST)
       CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               NACIST)
      ELSE
       CALL ICOPVE2(WORK(KISTSGP2(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               IACIST)
       CALL ICOPVE2(WORK(KNSTSGP2(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               NACIST)
      END IF
C?     WRITE(6,*) ' IACIST and NACIST arrays '
C?     CALL IWRTMA(IACIST,1,NSMST,1,NSMST)
C?     CALL IWRTMA(NACIST,1,NSMST,1,NSMST)
*
*. Generate offsets for I strings with given symmetry in
*  each space
*
      DO IGAS = 1, NGAS
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
      IFIRST = 1
      ISTRBS = 1
      NSTRINT = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 2001
        END IF
        IFIRST = 0
*. Symmetry of NGASL -1 spaces given, symmetry of full space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL SYMCOM_REL(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*.  sym of SPACE NGASL
        CALL SYMCOM_REL(2,1,ISTSMM1,ISMGSN,ISPGPSM)
        ISMFGS(NGASL) = ISMGSN
#if defined LUCI_DEBUG
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
#endif
*. Number of strings with this symmetry combination
        NSTRII = 1
        DO IGAS = 1, NGASL
          NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for this symmetry distribution in IOFFI
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGASL
          IOFF = IOFF + (ISMFGS(IGAS)-1)*MULT
          MULT = MULT * NSMST
        END DO
*
        call i_save_r(WORK(KLOFFI),NSTRINT+1,IOFF)
CTF     IOFFI(IOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
#if defined LUCI_DEBUG
        IF(NTEST.GE.200) THEN
          WRITE(6,'(A22,I8,I8,I8)')
     &    'IOFF,NSTRINT+1 NSTRII ',
     &     IOFF,IFRMR_REL(WORK(KLOFFI),1,IOFF),NSTRII
        END IF
#endif
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE
*
*. Gas structure of K strings
*
      NGASL = 1
      DO IGAS = 1, NGAS
       ITPFGS(IGAS) = ISPGPFTP(IGAS,KSPGRPABS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
*. Active group of K-strings
      KACGRP = ITPFGS(IOBTP)
*. Number of strings per symmetry distribution
      DO IGAS = 1, NGAS
        IF(IUB.EQ.1) THEN
        CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(WORK(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
        ELSE
        CALL ICOPVE2(WORK(KNSTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(WORK(KISTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
        END IF
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
 1000 CONTINUE
        IF(KFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        KFIRST = 0
#if defined LUCI_DEBUG
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
#endif
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL SYMCOM_REL(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE NGASL
        CALL SYMCOM_REL(2,1,ISTSMM1,ISMGSN,KSM)
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
#if defined LUCI_DEBUG
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
#endif
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
          NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IOBTP)
        CALL SYMCOM_REL(3,1,IOBSM,ISMFGS(IOBTP),IACSM)
        ISMFGS(IOBTP) = IACSM
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGAS
          IOFF = IOFF + (ISMFGS(IGAS)-1)*MULT
          MULT = MULT * NSMST
        END DO
        ISMFGS(IOBTP) = ISAVE
        IBSTRINI = IFRMR_REL(WORK(KLOFFI),1,IOFF)
*. Number of strings before active GAS space
        NSTB = 1
        DO IGAS = 1, IOBTP-1
          NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number of strings before active GAS space
        NSTA = 1
        DO IGAS =  IOBTP+1, NGAS
          NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number and offset for active group
C?      write(6,*) ' IACSM = ', IACSM
        NIAC  = NACIST(IACSM)
        IIAC =  IACIST(IACSM)
*
        NKAC = NNSTSGP(ISMFGS(IOBTP),IOBTP)
        IKAC = IISTSGP(ISMFGS(IOBTP),IOBTP)
*. I and K strings of given symmetry distribution
        NISD = NSTB*NIAC*NSTA
        NKSD = NSTB*NKAC*NSTA
C?      write(6,*) ' nstb nsta niac nkac ',
C?   &               nstb,nsta,niac,nkac
*. Obtain annihilation n mapping for all strings of this type
*
        IF(IUB.EQ.1) THEN
          NORBTS= NOBPTS(IOBTP,IOBSM)
        ELSE
          NORBTS= NOBPTS2(IOBTP,IOBSM)
        END IF
*
        NKACT = NSTFGP(KACGRP)
C?      write(6,*) ' KACGRP ', KACGRP
        IF(IUB.EQ.2) THEN
        CALL ADSTN_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,
     &                 WORK(KSTSTM2(KACGRP,1)),WORK(KSTSTM2(KACGRP,2)),
     &                 IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &                 NKSTR,KBSTRIN,NELB,NACGSOB,I1,XI1S,SCLFAC)
        ELSE
        CALL ADSTN_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,
     &                 WORK(KSTSTM(KACGRP,1)),WORK(KSTSTM(KACGRP,2)),
     &                 IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &                 NKSTR,KBSTRIN,NELB,NACGSOB,I1,XI1S,SCLFAC)
        END IF
*
        KSTRBS = KSTRBS + NKSD
        IF(NGASL-1.GT.0) GOTO 1000
 1001 CONTINUE
*
 9999 CONTINUE
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADSTN_GAS '
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
#endif
*
      call memmar(KLOFFI,LOFFI,'FREE',1,'LLOFFI')
C     CALL QEXIT('ADSTN ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ADTYP(NATP,JTYP,NGAS,LOCC,ROCC)
*
* Obtain type of J so
* <L! a+J !R> is nonvanishing
*
      INTEGER LOCC(NGAS),ROCC(NGAS)
#include "ipoist8.inc"
*. First ensure that the types differ by creation of one electron
      ICREA = 0
      IANNI = 0
      DO IGAS = 1,NGAS,1
        IF(LOCC(IGAS).LT.ROCC(IGAS)) THEN
         IANNI = IANNI + ROCC(IGAS)-LOCC(IGAS)
        ELSE IF (LOCC(IGAS).GT.ROCC(IGAS) ) THEN
         ICREA = ICREA + LOCC(IGAS) - ROCC(IGAS)
        END IF
      END DO
*
      IF(ICREA.EQ.1 .AND. IANNI.EQ.0) THEN
       NATP = 1
      ELSE
       NATP = 0
      END IF
*
      IF(NATP.EQ.1) THEN
        DO IGAS = 1,NGAS,1
          IF (LOCC(IGAS)-ROCC(IGAS) .EQ. 1) THEN
            IF(ICREA.EQ.1) THEN
             JTYP = IGAS
             ICREA = 1
            end if
          END IF
        END DO
      END IF
*
      NTEST = 0
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Single Creation connecting  '
        WRITE(6,*) ' Left string '
        CALL IWRTMA(LOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Right string '
        CALL IWRTMA(ROCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        IF(NATP.EQ.0) THEN
          WRITE(6,*) ' No connecting single creation'
        ELSE IF(NATP.EQ.1) THEN
          WRITE(6,*) ' Type of connecting creation op', JTYP
        END IF
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ATYP(NATP,JTYP,NGAS,LOCC,ROCC)
*
* Obtain type of J so
* <L! aJ !R> is nonvanishing
*
      INTEGER LOCC(NGAS),ROCC(NGAS)
#include "ipoist8.inc"
*. First ensure that the types differ by annihilation of one electron
      ICREA = 0
      IANNI = 0
      DO IGAS = 1,NGAS,1
        IF(LOCC(IGAS).LT.ROCC(IGAS)) THEN
         IANNI = IANNI + ROCC(IGAS)-LOCC(IGAS)
        ELSE IF (LOCC(IGAS).GT.ROCC(IGAS) ) THEN
         ICREA = ICREA + LOCC(IGAS) - ROCC(IGAS)
        END IF
      END DO
*
      IF(ICREA.EQ.0 .AND. IANNI.EQ.1) THEN
       NATP = 1
      ELSE
       NATP = 0
      END IF
*
      IF(NATP.EQ.1) THEN
        DO IGAS = 1,NGAS,1
          IF (LOCC(IGAS)-ROCC(IGAS) .EQ. -1 ) THEN
            IF(IANNI.EQ.1) THEN
             JTYP = IGAS
             IANNI = 1
            end if
          END IF
        END DO
      END IF
*
      NTEST = 0
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Single Annihilation connecting  '
        WRITE(6,*) ' Left string '
        CALL IWRTMA(LOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Right string '
        CALL IWRTMA(ROCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        IF(NATP.EQ.0) THEN
          WRITE(6,*) ' No connecting single annihilation'
        ELSE IF(NATP.EQ.1) THEN
          WRITE(6,*) ' Type of connecting annihilation op', JTYP
        END IF
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE CRESTR_GAS_REL(STRING,NSTINI,NSTINO,NEL,NORB,IORBOF,
     &                  Z,NEWORD,LSGSTR,ISGSTI,ISGSTO,TI,TTO,NACOB,
     &                  IPRNT)
*
* A type of strings containing NEL electrons are given
* set up all possible ways of adding an electron to this type of strings
*
*------------------------------------------
*  Jeppe Olsen
*
*   Fixed machine dependent problem with IF blocks.
*   They are now split up so that array dimensions are
*   not violated.
*  Timo Fleig, January 2001
*------------------------------------------
*
*========
* Input :
*========
* STRING : Input strings containing NEL electrons
* NSTINI : Number of input  strings
* NSTINO : Number of output strings
* NEL    : Number of electrons in input strings
* NORB   : Number of orbitals
* IORBOF : Number of first orbital
* Z      : Lexical ordering matrix for output strings containing
*          NEL + 1 electrons
* NEWORD : Reordering array for N+1 strings
* LSGSTR : .NE.0 => Include sign arrays ISGSTI,ISGSTO of strings
* ISGSTI : Sign array for NEL   strings
* ISGSTO : Sign array for NEL+1 strings
*
*=========
* Output :
*=========
*
*TI      : TI(I,ISTRIN) .gt. 0 indicates that orbital I can be added
*          to string ISTRIN .
*TTO     : Resulting NEL + 1 strings
*          if the string have a negative sign
*          then the phase equals - 1
      IMPLICIT REAL*8           (A-H,O-Z)
#include "ipoist8.inc"
C     INTEGER STRING,TI,TTO,OCC,UNOCC,STRIN2,Z
      INTEGER STRING,TI,TTO,STRIN2,Z
*.Input
      DIMENSION STRING(NEL,NSTINI),NEWORD(NSTINO),Z(NORB,NEL+1)
      DIMENSION ISGSTI(NSTINI),ISGSTO(NSTINO)
*.Output
      DIMENSION TI(NORB,NSTINI),TTO(NORB,NSTINI)
*.Scratch
      DIMENSION STRIN2(500)
*
      NTEST0 =  000
      NTEST = MAX(IPRNT,NTEST0)
#if defined LUCI_DEBUG
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' CRESTR speaking '
        WRITE(6,*)  ' =============== '
        WRITE(6,*)
         WRITE(6,*) ' Number of input electrons ', NEL
      END IF
#endif
C     WRITE(6,*) ' Reorder array NEWORD '
C     CALL IWRTMA(NEWORD,1,NSTINO,1,NSTINO)
*
      DO 1000 ISTRIN = 1,NSTINI
C?    write(6,*) ' Input string ',istrin,(string(i,istrin),i=1,nel)
        DO 100 IORB = IORBOF, IORBOF-1+NORB
C?      write(6,*) ' orbital ',iorb
           IPLACE = 0
           IF(NEL.EQ.0) THEN
             IPLACE = 1
             GOTO 11
           ELSE IF ( NEL .NE. 0 ) THEN
            DO 10 IEL = 1, NEL
              IF (IEL.eq.1.AND.IORB.lt.STRING(1,ISTRIN)) THEN
                IPLACE = 1
                GOTO 11
              ELSE IF (IEL.EQ.NEL.AND.IORB.GT.STRING(IEL,ISTRIN)) then
                IPLACE = IEL + 1
                goto 11
              else if (IORB.eq.STRING(IEL,ISTRIN)) THEN
                IPLACE = 0
                GOTO 11
              end if
              if (IEL.LT.NEL) then
                if (IORB.GT.STRING(IEL,ISTRIN).AND.
     &              IORB.LT.STRING(IEL+1,ISTRIN)) THEN
                  IPLACE = IEL+1
                  GOTO 11
                end if
              end if
   10       CONTINUE
           END IF
   11     CONTINUE
*
C?        write(6,*) ' iplace = ', iplace
          IF(IPLACE.NE.0) THEN
*. Generate next string
            DO 30 I = 1, IPLACE-1
   30       STRIN2(I) = STRING(I,ISTRIN)
            STRIN2(IPLACE) = IORB
            DO 40 I = IPLACE,NEL
   40       STRIN2(I+1) = STRING(I,ISTRIN)
C?          write(6,*) ' updated string (STRIN2) '
C?          call iwrtma(STRIN2,1,NEL+1,1,NEL+1)
            JSTRIN = ISTRNM(STRIN2,NACOB,NEL+1,Z,NEWORD,1)
C?          write(6,*) ' corresponding number ', JSTRIN
*
            TTO(IORB-IORBOF+1,ISTRIN) = JSTRIN
            IIISGN = (-1)**(IPLACE-1)
            IF(LSGSTR.NE.0)
     &      IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
            IF(IIISGN .EQ. -1 )
     &      TTO(IORB-IORBOF+1,ISTRIN) = - TTO(IORB-IORBOF+1,ISTRIN)
            TI(IORB-IORBOF+1,ISTRIN ) = IORB
          END IF
  100   CONTINUE
*
 1000 CONTINUE
*
#if defined LUCI_DEBUG
      IF ( NTEST .GE. 20) THEN
        MAXPR = 60
        NPR = MIN(NSTINI,MAXPR)
        WRITE(6,*) ' Output from CRESTR : '
        WRITE(6,*) '==================='
*
        WRITE(6,*)
        WRITE(6,*) ' Strings with an electron added  '
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' New strings.. ',
     &     (TTO(I,ISTRIN),I = 1,NORB)
        END DO
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' orbitals added or removed ' ,
     &     (TI(I,ISTRIN),I = 1,NORB)
        END DO
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GENSTR_GAS_REL(NEL,NELMN1,NELMX1,NELMN3,NELMX3,
     &                  ISTASO,IGRP,NOCTYP,NSMST,Z,LSTASO,
     &                  IREORD,STRING,IOC,IOTYP,IPRNT)
*
* Generate strings consisting of NEL electrons fulfilling
*   1 : Between NELMN1 AND NELMX1 electrons in the first NORB1 orbitals
*   2 : Between NELMN3 AND NELMX3 electrons in the last  NORB3 orbitals
*
* In the present version the strings are directly ordered into
* symmetry and occupation type .
*
* Jeppe Olsen Winter of 1990
*
* Special GAS version, Winter of 94 All strings of group IGRP
*
* ========
* Output :
* ========
* STRING(IEL,ISTRIN) : Occupation of strings.
* IREORD             : Reordering array going from lexical
*                      order to symmetry and occupation type order.
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION ISTASO(NSMST,*)
      INTEGER Z(NACOB,NEL)
*.Orbinp
#include "mxpdim.inc"
#include "orbinp.inc"
*
*.Output
      INTEGER STRING(NEL,*),IREORD(*)
*.Scratch arrays
      DIMENSION IOC(*),LSTASO(NOCTYP,NSMST)
*
      CALL ISETVC(LSTASO,0,NOCTYP*NSMST)
*
      NTEST0 = 000
      NTEST = MAX(NTEST0,IPRNT)
*
#if defined LUCI_DEBUG
      IF( NTEST .GE. 10) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' GENSTR speaking '
        WRITE(6,*)  ' =============== '
      END IF
#endif
*
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
* Loop over possible partitionings between RAS1,RAS2,RAS3
      DO 1001 IEL1 = NELMX1,NELMN1,-1
      DO 1003 IEL3 = NELMN3,NELMX3, 1
!      print *, 'IEL1, IEL3, NORB1, NORB3',IEL1, IEL3,NORB1, NORB3
       IF(IEL1.GT. NORB1 ) GOTO 1001
       IF(IEL3.GT. NORB3 ) GOTO 1003
       IEL2 = NEL - IEL1-IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
* Loop over RAS 1 occupancies
  901  CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
!           print *, ' calling ISTVC2',IOC(1)
            CALL ISTVC2(IOC(1),0,1,IEL1)
            IFRST1 = 0
           ELSE
             CALL NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
#if defined LUCI_DEBUG
         IF( NTEST .GE. 500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
#endif
         IFRST2 = 1
         IFRST3 = 1
* Loop over RAS 2 occupancies
  902    CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              CALL ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
              IFRST2 = 0
             ELSE
               CALL NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
#if defined LUCI_DEBUG
           IF( NTEST .GE. 500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
#endif
           IFRST3 = 1
* Loop over RAS 3 occupancies
  903      CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                CALL ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
                IFRST3 = 0
               ELSE
                 CALL NXTORD(IOC(IEL1+IEL2+1),
     &           IEL3,IORB3F,IORB3L,NONEW3)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
#if defined LUCI_DEBUG
             IF( NTEST .GE. 500 ) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
#endif
* Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
*. Symmetry
*                   ISYMST_REL(STRING,NEL)
!            WRITE(6,*) ' get ISYM'
             ISYM = ISYMST_REL(IOC,NEL)
!            WRITE(6,*) ' got ISYM',ISYM
C            IF(NTEST.GE.500)  WRITE(6,*) ' ISYM= ',ISYM
*. Occupation type
C            ITYP = IOCTP2(IOC,NEL,IOTYP)
             ITYP = 1
*
             IF(ITYP.NE.0) THEN
               LSTASO(ITYP,ISYM) = LSTASO(ITYP,ISYM)+ 1
C                      ISTRNM(IOCC,NACTOB,NEL,Z,NEWORD,IREORD)
               LEXCI = ISTRNM(IOC,NACOB,NEL,Z,IREORD,0)
               LACTU = ISTASO(ISYM,IGRP)-1+LSTASO(ITYP,ISYM)
               IREORD(LEXCI) = LACTU
!#if defined LUCI_DEBUG
               IF(NTEST.GT.50) WRITE(6,*) ' LEXCI,LACTU',
     &         LEXCI,LACTU
!#endif
               CALL ICOPVE(IOC,STRING(1,LACTU),NEL)
             END IF
*
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)
     &     GOTO 901
 1003 CONTINUE
 1001 CONTINUE
*
!#if defined LUCI_DEBUG
      IF(NTEST.GE.1 ) THEN
        WRITE(6,*) ' Number of strings generated   ', NSTRIN
      END IF
      IF(NTEST.GE.10)  THEN
        IF(NTEST.GE.100) THEN
          NPR = NSTRIN
        ELSE
          NPR = MIN(NSTRIN,50)
        END IF
        WRITE(6,*) ' Strings generated '
        WRITE(6,*) ' =================='
        ISTRIN = 0
        DO 100 ISYM = 1, NSMST
        DO 100 ITYP = 1,NOCTYP
          LSTRIN = MIN(LSTASO(ITYP,ISYM),NPR-ISTRIN)
          IF(LSTRIN.GT.0) THEN
            WRITE(6,*) ' Strings of type and symmetry ',ITYP,ISYM
            DO 90 KSTRIN = 1,LSTRIN
              ISTRIN = ISTRIN + 1
              WRITE(6,'(2X,I4,8X,(10I5))')
     &        ISTRIN,(STRING(IEL,ISTRIN),IEL = 1,NEL)
   90       CONTINUE
          END IF
  100   CONTINUE
*
        WRITE(6,*) ' Array giving actual place from lexical place'
        WRITE(6,*) ' ============================================'
        CALL IWRTMA(IREORD,1,NPR,1,NPR)
      END IF
!#endif

      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETSTR2_TOTSM_SPGP_REL(IUB,IGRP,NIGRP,ISPGRPSM,NEL,
     &                                  NSTR,ISTR,
     &                                  NORBT,IDOREO,IZ,IREO)
*
* Obtain all super-strings of given total symmetry and given
* occupation in each GAS space 
*
*.If  IDOREO .NE. 0 THEN reordering array : 
*                   lexical => actual order is obtained
*
* Nomenclature of the day : superstring : string in complete 
*                           orbital space, product of strings in
*                           each GAS space 
*
* Compared to GETSTR2_TOTSM_SPGP : Based upon IGRP(NIGRP)
*                                  (Just a few changes in the beginning)
*
* =====
* Input 
* =====
*
* IGRP :  supergroup, here as an array of GAS space 
* NIGRP : Number of active groups 
* ISPGRPSM : Total symmetry of superstrings 
* NEL : Number of electrons 
* IZ  : Reverse lexical ordering array for this supergroup 
*       (IF IDOREO.NE.0)
* 
*
* ======
* Output 
* ======
*
* NSTR : Number of superstrings generated
* ISTR : Occupation of superstring
* IREO : Reorder array ( if IDOREO.NE.0) 
*
*
* Jeppe Olsen, Written  July 1995
*              Version of Dec 1997
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
*. Input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "strbas.inc"
#include "csm.inc"
      INTEGER IZ(NORBT,NEL)
      INTEGER IGRP(NIGRP)
*. output
      INTEGER ISTR(NEL,*), IREO(*)
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
C     CALL QENTER('GETST')
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Welcome to GETSTR_TOTSM_SPGP '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
        WRITE(6,'(A)')  ' Strings to be obtained : ' 
        WRITE(6,'(A)')  ' ======================== '
        WRITE(6,'(A)')
        WRITE(6,'(A,I2)') '   Symmetry : ', ISPGRPSM
        WRITE(6,'(A,16I3)') ' Groups : ', (IGRP(I),I=1,NIGRP)
        WRITE(6,*) ' NEL = ', NEL
        IF(IDOREO.NE.0) THEN 
          WRITE(6,*)
          WRITE(6,*) ' ============= '
          WRITE(6,*) ' The Z array : '
          WRITE(6,*) ' ============= '
          WRITE(6,*)
          WRITE(6,*) ' NORBT,NEL = ',NORBT,NEL
          CALL IWRTMA(IZ,NORBT,NEL,NORBT,NEL)
        END IF
      END IF
#endif
*. Absolut number of this supergroup
*. Occupation per gasspace
*. Largest occupied space 
      NGASL = 0
*. Largest and lowest symmetries active in each GAS space
      DO IGAS = 1, NGAS
        ITPFGS(IGAS) = IGRP(IGAS)
        NELFGS(IGAS) = NELFGP(IGRP(IGAS))          
        IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
      IF(NGASL.EQ.0) NGASL = 1
*. Number of strings per GAS space and offsets for strings of given sym
      DO IGAS = 1, NGAS
        if (IUB.eq.1) then
          CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,IGAS))
          CALL ICOPVE2(WORK(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 IISTSGP(1,IGAS))
        else
          CALL ICOPVE2(WORK(KNSTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,IGAS))
          CALL ICOPVE2(WORK(KISTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 IISTSGP(1,IGAS))
        end if
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST =1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
* Largest and lowest active symmetries for each GAS space
#if defined LUCI_DEBUG
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Type of each GAS space '
         CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of elecs per GAS space '
         CALL IWRTMA(NELFGS,1,NGAS,1,NGAS)
      END IF 
#endif
*
*. Loop over symmetries of each GAS
*
      MAXLEX = 0
      IFIRST = 1
      ISTRBS = 1
 1000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1 
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        IFIRST = 0
#if defined LUCI_DEBUG
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
#endif
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL SYMCOM_REL(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE NGASL 
        CALL SYMCOM_REL(2,1,ISTSMM1,ISMGSN,ISPGRPSM)
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
#if defined LUCI_DEBUG
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
#endif
*. Obtain all strings of this symmetry 
        CALL GETSTRN_GASSM_SPGP_REL(ISMFGS,ITPFGS,ISTR(1,ISTRBS),
     &                              NSTR,NEL,
     &                              NNSTSGP,IISTSGP,IUB)
*. Reorder Info : Lexical => actual number 
        IF(IDOREO.NE.0) THEN
*. Lexical number of NEL electrons
*. Can be made smart by using common factor for first NGAS-1 spaces 
          DO JSTR = ISTRBS, ISTRBS+NSTR-1
            LEX = 1
            DO IEL = 1, NEL 
              LEX = LEX + IZ(ISTR(IEL,JSTR),IEL)
            END DO
C?          WRITE(6,*) ' string '  
C?          CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
C?          WRITE(6,*) ' JSTR and LEX ', JSTR,LEX
*
            MAXLEX = MAX(MAXLEX,LEX)
            IREO(LEX) = JSTR
          END DO
        END IF
*
        ISTRBS = ISTRBS + NSTR 
*. ready for next symmetry distribution 
        IF(NGAS-1.NE.0) GOTO 1000
 1001 CONTINUE
*. End of loop over symmetry distributions
      NSTR = ISTRBS - 1
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NEL(b) = ', NEL
        WRITE(6,*) ' Number of strings generated ', NSTR
        WRITE(6,*)
        WRITE(6,*) ' Strings : '
        WRITE(6,*)
        CALL PRTSTR(ISTR,NEL,NSTR)
*
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) 'Largest Lexical number obtained ', MAXLEX
          WRITE(6,*) ' Reorder array '
          CALL IWRTMA(IREO,1,NSTR,1,NSTR)
        END IF
      END IF
#endif
*
C     CALL QEXIT('GETST')
      RETURN
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETSTR_TOTSM_SPGP_REL(ISTRTP,ISPGRP,ISPGRPSM,
     &                                 NEL,NSTR,ISTR,
     &                                 NORBT,IDOREO,IZ,IREO)
*
* Obtain all super-strings of given total symmetry and given
* occupation in each GAS space
*
*.If  IDOREO .NE. 0 THEN reordering array : lexical => actual order is obtained
*
* Nomenclature of the day : superstring : string in complete
*                           orbital space, product of strings in
*                           each GAS space
* =====
* Input
* =====
*
* ISTRTP  : Type of of superstrings ( alpha => 1, beta => 2 )
* ISPGRP :  supergroup number, (relative to start of this type )
* ISPGRPSM : Total symmetry of superstrings
* NEL : Number of electrons
* IZ  : Reverse lexical ordering array for this supergroup
*
*
* ======
* Output
* ======
*
* NSTR : Number of superstrings generated
* ISTR : Occupation of superstring
* IREO : Reorder array ( if IDOREO.NE.0)
*
*
* Jeppe Olsen, July 1995
*
* Double group modifications (april 98 ) : Separate treatment of alpha and beta
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
*. Input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "strbas.inc"
#include "csm.inc"
*
#include "krmcluci_inf.h"
*
      INTEGER IZ(NORBT,NEL)
*. output
      INTEGER ISTR(NEL,*), IREO(*)
*. Local scratch
      INTEGER NELFGS(MXPNGAS),ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
C     CALL QENTER('GETST')
*
      IF(MOD(ISTRTP,2).EQ.0) THEN
*. Barred spinor
       IUB = 2
      ELSE
       IUB = 1
      END IF
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Welcome to GETSTR_TOTSM_SPGP '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
        WRITE(6,'(A,3I3)')
     & ' Strings to be obtained : Type, supergroup, symmetry ',
     &   ISTRTP,ISPGRP,ISPGRPSM
        WRITE(6,*)
        WRITE(6,*) ' IUB = ', IUB
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) ' IZ matrix '
          CALL IWRTMA(IZ,NORBT,NEL,NORBT,NEL)
        END IF
      END IF
#endif
*. Absolut number of this supergroup
      ISPGRPA = IBSPGPFTP(ISTRTP) - 1 + ISPGRP
*. Occupation per gasspace
*. Largest occupied space
      NGASL = 0
*. Largest and lowest symmetries active in each GAS space
      DO IGAS = 1, NGAS
        ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPA)
        NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
        IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
      IF(NGASL.EQ.0) NGASL = 1
*. Number of strings per GAS space and offsets for strings of given sym
      DO IGAS = 1, NGAS
        IF(IUB.EQ.1) THEN
         CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                NNSTSGP(1,IGAS))
         CALL ICOPVE2(WORK(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                IISTSGP(1,IGAS))
        ELSE
         CALL ICOPVE2(WORK(KNSTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                NNSTSGP(1,IGAS))
         CALL ICOPVE2(WORK(KISTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                IISTSGP(1,IGAS))
        END IF
*
      END DO
*
* Initialize and set MINVAL and MAXVAL
      IZERO = 0
      DO IGAS = 1, NGAS
        MINVAL(IGAS) = IZERO
        MAXVAL(IGAS) = IZERO
        DO ISMST =1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) then
            MAXVAL(IGAS) = ISMST
          end if
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) then
            MINVAL(IGAS) = ISMST
          end if
        END DO
      END DO
#if defined LUCI_DEBUG
      if (NTEST.ge.300) then
        write(6,*) 'NSMST :  ',NSMST
        do IGAS = 1, NGAS
          do ISMST =1, NSMST
            write(6,*) 'NNSTSGP(',ISMST,IGAS,') : ',NNSTSGP(ISMST,IGAS)
          end do
        end do
      end if
#endif
*
CTF
* If no strings are found in a given GAS, we set the symmetry
* value in MINVAL to 1, just to prevent crashes (only a test,
* waiting for more elegant solution ...)
      IONE = 1
      do IGAS = 1,NGAS-1,1
        if (MINVAL(IGAS).eq.0) then
          MINVAL(IGAS) = IONE
        end if
      end do
CTF
*
* Largest and lowest active symmetries for each GAS space
#if defined LUCI_DEBUG
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Type of each GAS space '
         CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of elecs per GAS space '
         CALL IWRTMA(NELFGS,1,NGAS,1,NGAS)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' MINVAL and MAXVAL : '
        WRITE(6,*) ' =================== '
        DO I = 1, NGAS
          WRITE(6,'(2I5)') MINVAL(I),MAXVAL(I)
        END DO
      END IF
#endif
*
*. Loop over symmetries of each GAS
*
      MAXLEX = 0
      IFIRST = 1
      ISTRBS = 1
      NGASCT = NGAS
 1000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
C        NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         do IGS=1,NGAS,1
            if (ISMFGS(IGS).gt.NIRR_DG_LOC) then
#if defined LUCI_DEBUG
              if (NTEST.ge.1) then
                write(6,*) 'Problem in getstr_totsm_spgp'
                write(6,*) 'NIRR_DG_LOC = ',NIRR_DG_LOC
                write(6,*) 'ISMFGS out of range.'
                write(6,*) 'ISMFGS(',IGS,') : ',ISMFGS(IGS)
                write(6,*) 'setting symmetry to 1'
              end if
#endif
              ISMFGS(IGS) = 1
            end if
         end do
         IF(NONEW.NE.0) GOTO 1001
        END IF
        IFIRST = 0
#if defined LUCI_DEBUG
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
#endif
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
#if defined LUCI_DEBUG
          if (NTEST.ge.300) then
            write(6,*) 'ISMFGS(',IGAS,') in getstr_totsm_spgp',
     &                  ISMFGS(IGAS)
          end if
#endif
C?        write(6,*) ' ISTSMM1 : ', ISTSMM1
          CALL SYMCOM_REL(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
C              SYMCOM(ITASK,IOBJ,I1,I2,I12)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE NGASL
#if defined LUCI_DEBUG
        if (NTEST.ge.300) then
          write(6,*) 'ISMGSN in getstr_totsm_spgp',ISMGSN
        end if
#endif
        CALL SYMCOM_REL(2,1,ISTSMM1,ISMGSN,ISPGRPSM)
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
#if defined LUCI_DEBUG
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
#endif
*. Obtain all strings of this symmetry
CT      CALL QENTER('GASSM')
        CALL GETSTRN_GASSM_SPGP_REL(ISMFGS,ITPFGS,ISTR(1,ISTRBS),
     &                              NSTR,NEL,
     &                              NNSTSGP,IISTSGP,IUB)
CT      CALL QEXIT('GASSM')
*. Reorder Info : Lexical => actual number
        IF(IDOREO.NE.0) THEN
*. Lexical number of NEL electrons
*. Can be made smart by using common factor for first NGAS-1 spaces
          DO JSTR = ISTRBS, ISTRBS+NSTR-1
            LEX = 1
            DO IEL = 1, NEL
              LEX = LEX + IZ(ISTR(IEL,JSTR),IEL)
            END DO
C           WRITE(6,*) ' string '
C           CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
C           WRITE(6,*) ' JSTR and LEX ', JSTR,LEX
*
            MAXLEX = MAX(MAXLEX,LEX)
            IREO(LEX) = JSTR
          END DO
        END IF
*
        ISTRBS = ISTRBS + NSTR
*. ready for next symmetry distribution
#if defined LUCI_DEBUG
        if (NTEST.ge.1) then
          write(6,*) 'Loop status in getstr_totsm_spgp'
          write(6,*) 'NGAS-1 is  ',NGAS-1
        end if
#endif
        NGASCT = NGASCT - 1
        if (NGAS-1.ne.0) then
c         if (NGASCT.eq.0) goto 1001
          GOTO 1000
        else
          goto 1001
        end if
 1001 CONTINUE
*. End of loop over symmetry distributions
      NSTR = ISTRBS - 1
*. Check for zeroes
       I_DO_ZERO_CHECK = 0
       IF(I_DO_ZERO_CHECK.EQ.1) THEN
         DO JSTR = 1, NSTR
           DO JEL = 1, NEL
             IF(ISTR(JEL,JSTR).EQ.0) THEN
                WRITE(6,*) ' Occupation string with zeroes '
                CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
                Call Abend2(       ' Occupation string with zeroes ' )
              END IF
            END DO
         END DO
       END IF
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings generated ', NSTR
        WRITE(6,*)
        WRITE(6,*) ' Strings : '
        WRITE(6,*)
        if (NSTR.ne.0.and.NEL.ne.0) CALL PRTSTR(ISTR,NEL,NSTR)
*
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) 'Largest Lexical number obtained ', MAXLEX
          WRITE(6,*) ' Reorder array '
          CALL IWRTMA(IREO,1,NSTR,1,NSTR)
        END IF
      END IF
#endif
*
 1002 continue
*
C     CALL QEXIT('GETST')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETSTRN_GASSM_SPGP_REL(ISMFGS,ITPFGS,ISTROC,NSTR,NEL,
     &                                  NNSTSGP,IISTSGP,IUB)
*
* Obtain all superstrings containing  strings of given sym and type
*
* ( Superstring :contains electrons belonging to all gasspaces
*        string :contains electrons belonging to a given GAS space
* A super string is thus a product of NGAS strings )
*
* Jeppe Olsen, Summer of 95
*              Optimized version, october 1995
*              Option for separate treatment of up- and down- spinors added
*              April 98
*
*. In this subroutine the ordering of strings belonging to a given type
*  is defined !!
* Currently we are using the order
* Loop over GAS 1 strings
*  Loop over GAS 2 strings
*   Loop over GAS 3 strings --
*
*     Loop over gas N strings
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "csm.inc"
*. Specific input
      INTEGER ITPFGS(*), ISMFGS(*)
      INTEGER NNSTSGP(MXPNSMST,*), IISTSGP(MXPNSMST,*)
*. Local scratch
C     INTEGER NSTFGS(MXPNGAS), IBSTFGS(MXPNGAS), ISTRNM(MXPNGAS)
      INTEGER NSTFGS(MXPNGAS), IBSTFGS(MXPNGAS)
*. Output
      INTEGER ISTROC(NEL,*)
*. Number of strings per GAS space
      DO IGAS = 1, NGAS
        NSTFGS(IGAS)  = NNSTSGP(ISMFGS(IGAS),IGAS)
        IBSTFGS(IGAS) = IISTSGP(ISMFGS(IGAS),IGAS)
      END DO
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GETSTR_GASSM_SPGP speaking '
        WRITE(6,*) '  =========================== '
        WRITE(6,*) ' ISMFGS,ITPFGS (input) '
        CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' NSTFGS, IBSTFGS ( intermediate results ) '
        CALL IWRTMA(NSTFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(IBSTFGS,1,NGAS,1,NGAS)
      END IF
#endif
*. Last gasspace with a nonvanishing number of electrons
      IGASL = 0
      DO IGAS = 1, NGAS
        IF( NELFGP(ITPFGS(IGAS)) .NE. 0 ) IGASL = IGAS
      END DO
*
CTI  In certain cases, types are generated that result in all
CTI  zero valued elements of NELFGP. 
CTI  Set following stuff to zero in such cases.
      if (IGASL.ne.0) then
         NELL = NELFGP(ITPFGS(IGASL))
         NELML = NEL - NELL
         NSTRGASL = NSTFGS(IGASL)
         IBGASL = IBSTFGS(IGASL)
      else
         NELL = 0
         NELML = NEL - NELL
         NSTRGASL = 0
         IBGASL = 0
      end if
*
      NSTRTOT = 1
      DO IGAS = 1, NGAS
        NSTRTOT = NSTRTOT*NSTFGS(IGAS)
      END DO
*
      IF(NSTRTOT.EQ.0) GOTO 1001
*. Loop over GAS spaces
      DO IGAS = 1, IGASL
*. Number of electrons in GAS = 1, IGAS - 1
        IF(IGAS.EQ.1) THEN
          NELB = 0
        ELSE
          NELB = NELB +  NELFGP(ITPFGS(IGAS-1))
        END IF
*. Number of electron in IGAS
        NELI = NELFGP(ITPFGS(IGAS))
        IF(NELI.GT.0) THEN

*. The order of strings corresponds to a matrix A(I(after),Igas,I(before))
*. where I(after) loops over strings in IGAS+1 - IGASL and
*  I(before) loop over strings in 1 - IGAS -1
          NSTA = 1
          DO JGAS = IGAS+1, IGASL
            NSTA = NSTA * NSTFGS(JGAS)
          END DO
*
          NSTB =  1
          DO JGAS = 1, IGAS-1
            NSTB = NSTB * NSTFGS(JGAS)
          END DO
*
          IF(IUB.EQ.1) THEN
            KOCSTR_ADD = KOCSTR(ITPFGS(IGAS))
          ELSE
            KOCSTR_ADD = KOCSTR2(ITPFGS(IGAS))
          END IF
*
C?        call memchk
*
          NSTI = NSTFGS(IGAS)
          CALL ADD_STR_GROUP_REL(NSTI,IBSTFGS(IGAS),WORK(KOCSTR_ADD),
     &                       NSTB,NSTA,ISTROC,NELB+1,NELI,NEL)
*. Loop over strings in IGAS
        END IF
      END DO
 1001 CONTINUE
      NSTR = NSTRTOT
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from  GETSTR_GASSM_SPGP '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry and type strings : '
        WRITE(6,*)
        WRITE(6,*) '   AS    Sym  Type '
        WRITE(6,*) ' =================='
        DO IGAS = 1, NGAS
          WRITE(6,'(3I6)') IGAS,ISMFGS(IGAS),ITPFGS(IGAS)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Number of strings generated : ', NSTR
        WRITE(6,*) ' Strings generated '
        CALL PRTSTR(ISTROC,NEL,NSTR)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GXFSTR(ILA,ILB,IRA,IRB,ICA,IAA,ICB,IAB,NGAS,ICON)
*
* Left  type defined by ILA, ILB
* Right type defined by IRA, IRB
* 
* Excitation operator defined by ICA, IAA, ICB,IAB
*
* Does this operator connect the two types ?
*
* ICON = 1 => Connection
* ICON = 0 => No connection
*
*
* Jeppe Olsen, Summer of 99
*
#include "implicit.inc"
*. Input
      INTEGER ILA(NGAS),ILB(NGAS),IRA(NGAS),IRB(NGAS)
      INTEGER ICA(NGAS),ICB(NGAS),IAA(NGAS),IAB(NGAS)
*
      ICON = 1
      DO IGAS = 1, NGAS
        IF(.NOT.( ICA(IGAS) - IAA(IGAS) + IRA(IGAS).EQ.ILA(IGAS)))
     &  ICON = 0
        IF(.NOT.( ICB(IGAS) - IAB(IGAS) + IRB(IGAS).EQ.ILB(IGAS)))
     &  ICON = 0
      END DO
*
#if defined LUCI_DEBUG
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'GXFSTR : IRA, IRB, ILA, ILB '
        CALL IWRTMA(IRA,1,NGAS,1,NGAS)
        CALL IWRTMA(IRB,1,NGAS,1,NGAS)
        CALL IWRTMA(ILA,1,NGAS,1,NGAS)
        CALL IWRTMA(ILB,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' ICA, IAA, ICB, IAB '
        CALL IWRTMA(ICA,1,NGAS,1,NGAS)
        CALL IWRTMA(IAA,1,NGAS,1,NGAS)
        CALL IWRTMA(ICB,1,NGAS,1,NGAS)
        CALL IWRTMA(IAB,1,NGAS,1,NGAS)
        IF(ICON.EQ.1) THEN
          WRITE(6,*) ' We have contact '
        ELSE
          WRITE(6,*) ' No contact '
        END IF
      END IF
#endif
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION ISYMS1_REL(STRING,NEL)
*
* Symmmetry of string, D2H version
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
C     PARAMETER(MXPORB = 500 )
C     PARAMETER(MXPR4T = 10)
C     PARAMETER(MXPOBS = 20)
C     PARAMETER(MXPIRR = 20 )
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB)
#include "mxpdim.inc"
#include "orbinp.inc"
*
      INTEGER SYMPRO(8,8)
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
*. Specific input
      INTEGER STRING(*)
*
      ISYM = 1
      DO 100 IEL = 1, NEL
        ISYM = SYMPRO(ISYM,ISMFTO(STRING(IEL)))
  100 CONTINUE
*
      ISYMS1_REL = ISYM
*
#if defined LUCI_DEBUG
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ISYMS1, String and symmetry '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,*) ISYM
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION ISYMS2_REL(STRING,NEL)
*
* Symmetry of string STRING, D inf h, C inf v, O3 version
*
*
      IMPLICIT REAL*8(A-H,O-Z)
C     PARAMETER(MXPORB = 500 )
c     PARAMETER(MXPR4T = 10)
C     PARAMETER(MXPOBS = 20)
C     PARAMETER(MXPIRR = 20)
c     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB)
#include "mxpdim.inc"
#include "orbinp.inc"
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONABR/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*
      INTEGER STRING(NEL)
*. ML and parity of string
      MLSTR = 0
      IPARI = 1
      DO 10 IEL = 1, NEL
        IF(ISMFTO(STRING(IEL)).LE.NMLOB) THEN
          MLSTR = MLSTR + ISMFTO(STRING(IEL))-1+MNMLOB
        ELSE
          MLSTR = MLSTR + ISMFTO(STRING(IEL))-1+MNMLOB-NMLOB
          IPARI = - IPARI
        END IF
   10 CONTINUE
*
      IF(IPARI.EQ.-1) IPARI = 2
      ISYM  = (IPARI-1) * NMLST+ MLSTR - MNMLST + 1
      ISYMS2_REL = ISYM
*
      NTEST = 0
#if defined LUCI_DEBUG
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' STRING '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,'(A,3I3)') ' MLSTR, IPARI ISYMS2 ', MLSTR,IPARI,ISYM
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION ISYMSDG(STRING,NEL)
*
* Symmmetry of string, Double group  version
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "symm.inc"
*

*. Specific input
      INTEGER STRING(*)
*. Note : How to switch between up and down spinors
      ISYM = 1
      DO IEL = 1, NEL
        ISYM = IDBGMULT(ISYM,ISMFTO(STRING(IEL)))
      END DO
*
      ISYMSDG = ISYM
*
      NTEST = 0
#if defined LUCI_DEBUG
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ISYMS1, String and symmetry '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,*) ISYM
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION ISYMST_REL(STRING,NEL)
*
* Master routine for symmetry of string
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input ( PNTGRP is used )
C     INTEGER PNTGRP,EXTSPC
C     PARAMETER(MXPIRR = 20)
C     PARAMETER ( MXPOBS = 20 )
C     PARAMETER (MXPR4T = 10 )
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C
#include "mxpdim.inc"
#include "lucinp.inc"
*. Specific input
      INTEGER STRING(*)
      IF(PNTGRP.EQ.1) THEN
*.D2h
        ISYMST_REL = ISYMS1_REL(STRING,NEL)
      ELSE IF(PNTGRP.GE.2.AND.PNTGRP.LE.4) THEN
*.Cinfv Dinfh O3
        ISYMST_REL = ISYMS2_REL(STRING,NEL)
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
        ISYMST_REL = ISYMSDG(STRING,NEL)
      ELSE
        WRITE(6,*) ' Sorry PNTGRP option not programmed ', PNTGRP
        WRITE(6,*) ' Enforced stop in ISYMST '
        Call Abend1( 5 )
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MEMSTR_GAS_REL
*
*
* Construct pointers for saving information about strings and
* their mappings
*
* GAS version
*
*========
* Input :
*========
* Number and groups of strings defined by /GASSTR/
* Symmetry information stored in         /CSM/
* String information stored in           /STINF/
*=========
* Output
*=========
* Pointers stored in common block /STRBAS/
*
* Jeppe Olsen , Winter of 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "strbas.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "stinf.inc"

*
*. Start of string information
      call memmar(KSTINF,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      NTEST = 0
#if defined LUCI_DEBUG
      IF(NTEST.NE.0)
     &WRITE(6,*) ' First word with string information',KSTINF
#endif
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        NSTRIN = NSTFGP(IGRP)
        LSTRIN = NSTRIN*NELFGP(IGRP)
        call memmar(KOCSTR(IGRP),LSTRIN,'ADDS  ',1,'OCSTR ')
        call memmar(KSTREO(IGRP),NSTRIN,'ADDS  ',1,'STREO ')
        IF(PNTGRP.GE.5) THEN
        call memmar(KOCSTR2(IGRP),LSTRIN,'ADDS  ',1,'OCSTR2')
        call memmar(KSTREO2(IGRP),NSTRIN,'ADDS  ',1,'STREO2')
        END IF
*. For spin down strings
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
      call memmar(KNSTSGP(1),NSMST*NGRP,'ADDS  ',1,'NSTSGP')
      call memmar(KISTSGP(1),NSMST*NGRP,'ADDS  ',1,'ISTSGP')
      IF(PNTGRP.GE.5) THEN
        call memmar(KNSTSGP2(1),NSMST*NGRP,'ADDS  ',1,'NSTSGP')
        call memmar(KISTSGP2(1),NSMST*NGRP,'ADDS  ',1,'ISTSGP')
      END IF


*
*. Number of strings per symmetry and offset for strings of given sym
*. for types
*
      DO ITP = 1, NSTTP
        call memmar(KNSTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'NSTSO ')
        call memmar(KISTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'ISTSO ')
      END DO
      IF(PNTGRP.GE.5) THEN
       DO ITP = 1, NSTTP
         LEN_T = NSPGPFTP(ITP)*NSMST
         call memmar(KNSTSO2(ITP),LEN_T,'ADDS  ',1,'NSTSO2')
         call memmar(KISTSO2(ITP),LEN_T,'ADDS  ',1,'ISTSO2')
       END DO
      END IF
*
**. Lexical adressing of arrays : use array indeces for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP
        call memmar(KZ(IGRP),NACOB*NELFGP(IGRP),'ADDS  ',1,'Zmat  ')
      END DO
*
*. Mappings between different groups
*
      DO  IGRP = 1, NGRP
        IEL = NELFGP(IGRP)
        IGAS = IGSFGP(IGRP)
        IORB = NOBPT(IGAS)
        ISTRIN = NSTFGP(IGRP)
*. Use full orbital notation
*. Note Only creation maps are currently installed,
*. Generalize when/if annihilation maps are included
        IF(ISTAC(IGRP,2).NE.0) THEN
          LENGTH = IORB*ISTRIN
          call memmar(KSTSTM(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
          call memmar(KSTSTM(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
          IF(PNTGRP.GE.5) THEN
            call memmar(KSTSTM2(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
            call memmar(KSTSTM2(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
          END IF
        ELSE
          KSTSTM(IGRP,1) = -1
          KSTSTM(IGRP,2) = -1
          IF(PNTGRP.GE.5) THEN
            KSTSTM2(IGRP,1) = -1
            KSTSTM2(IGRP,2) = -1
          END IF
        END IF
      END DO
*
*. Symmetry of conjugated orbitals and orbital excitations
*
*     KCOBSM,KNIFSJ,KIFSJ,KIFSJO
      call memmar(KCOBSM,NACOB,'ADDS  ',1,'Cobsm ')
      call memmar(KNIFSJ,NACOB*NSMSX,'ADDS  ',1,'Nifsj ')
      call memmar(KIFSJ,NACOB**2 ,'ADDS  ',1,'Ifsj  ')
      call memmar(KIFSJO,NACOB*NSMSX,'ADDS  ',1,'Ifsjo ')
*
*. Symmetry of excitation connecting  strings of given symmetry
*
      call memmar(KSTSTX,NSMST*NSMST,'ADDS  ',1,'Ststx ')
*
*. Occupation classes
*
      call memmar(KIOCLS,NMXOCCLS*NGAS,'ADDS  ',1,'IOCLS ')
*. Annihilation/Creation map of supergroup types
      call memmar(KSPGPAN,NTSPGP*NGAS,'ADDS  ',1,'SPGPAN')
      call memmar(KSPGPCR,NTSPGP*NGAS,'ADDS  ',1,'SPGPCR')
*
*. Last word of string information
*
      call memmar(KSTINE,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
#if defined LUCI_DEBUG
      IF(NTEST.NE.0)
     &WRITE(6,*) ' Last word with string information',KSTINE-1
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' R*8 words used for storing Strinformation ',
     &               KSTINE - KSTINF
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE NST_SPGRP_DBG(NIGRP,IGRP,ISM_TOT,NSMST,NSTRIN,NDIST,
     &                         IAB)
      use memory_allocator
*
* Number of strings for given combination of groups and 
* symmetry.
*
*. Input
*        
*
*   NIGRP : Number of active groups 
*   IGRP : The active groups
*   ISM_TOT : Total symmetry of supergroup
*   NSMST   : Number of string symmetries
*   IAB : = 1 => unbarred operators
*         = 2 => barred operators 
*
*. Output
*
*  NSTRIN : Number of strings with symmetry ISM_TOT
*  NDIST  : Number of symmetry distributions
*
* Jeppe Olsen, September 1997
*
*. Double group version
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
*. Specific Input
      DIMENSION IGRP(NIGRP)
*. General input
#include "gasstr.inc"
C
C     some scratch
CSK   INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
C     SK - Jan 2010: introduced use of dynamic memory allocation using 
C     the alloc interface.
C     This removes a spurious crash for CI calculations with a 
C     larger number of GA spaces.
C     As a rule, static allocations have to be avoided if a routine 
C     is called billion of times...
      integer, allocatable :: ISM(:)
      integer, allocatable :: MNSM(:)
      integer, allocatable :: MXSM(:)
C
C     allocate scratch memory
      call alloc(ISM, MXPNGAS)
      call alloc(MNSM,MXPNGAS)
      call alloc(MXSM,MXPNGAS)
#if defined LUCI_DEBUG
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' NST_SPGRP_DBG is speaking '
        WRITE(6,*) ' =========================='
*
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
      END IF
#endif
*
      IF(NIGRP.EQ.0) THEN
        IF(ISM_TOT.EQ.1) THEN 
          LENGTH = 1
        ELSE
          LENGTH = 0
        END IF
        GOTO 1001
      END IF
*
*. Set up min and max values for symmetries
      CALL MINMAX_FOR_SYM_DIST_DBG(NIGRP,IGRP,MNSM,MXSM,NDISTX,IAB)
*. Loop over symmetry distributions
      IFIRST = 1
      LENGTH = 0 
      NDIST = 0
*. Last group with symmetry differing from total symmetric
      NIGRPL = 1
      DO JGRP = 1, NIGRP
        IF(MXSM(JGRP).GT.1) NIGRPL = JGRP
      END DO
*. Number of strings in groups after NIGRPL
      NSTRL = 1
      DO JGRP = NIGRPL+1,NIGRP
        IF(IAB.EQ.1) THEN 
          NSTRL = NSTRL*NSTFSMGP(1,IGRP(JGRP))
        ELSE  
          NSTRL = NSTRL*NSTFSMGP2(1,IGRP(JGRP))
        END IF
      END DO
 1000 CONTINUE
*. Next symmetry distribution
* \/  NEXT_SYM_DISTR_REL should work for double group symmetry
        CALL NEXT_SYM_DISTR_REL(NIGRPL,MNSM,MXSM,ISM,ISM_TOT,IFIRST,
     &                          NONEW)
        IF(NONEW.EQ.0) THEN
          LDIST = NSTRL
          DO JGRP = 1, NIGRPL
            IF(IAB.EQ.1) THEN
              LDIST = LDIST*NSTFSMGP(ISM(JGRP),IGRP(JGRP))
            ELSE 
              LDIST = LDIST*NSTFSMGP2(ISM(JGRP),IGRP(JGRP))
            END IF
          END DO
          LENGTH = LENGTH + LDIST
          NDIST = NDIST + 1
      GOTO 1000
        END IF
*
 1001 CONTINUE
      NSTRIN = LENGTH
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings obtained ', LENGTH
        WRITE(6,*) ' Number of symmetry-distributions',NDIST
      END IF
#endif
C
C     release scratch memory
      call dealloc(ISM)
      call dealloc(MNSM)
      call dealloc(MXSM)
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE NSTPTP_GAS_REL(NGAS,ISPGRP,NSTSGP,NSMST,
     &                          NSTSSPGP,IGRP,MXNSTR)
*
* From number of strings per group and sym to number of strings
* per supergroup and sym for given super group.
*
* Jeppe Olsen , Fall of 94
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION ISPGRP(NGAS),NSTSGP(NSMST,*)
*. Output
      DIMENSION NSTSSPGP(NSMST,IGRP)
*. Scratch
#include "mxpdim.inc"
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NSTPTP_GAS is speaking '
        WRITE(6,*) ' ======================='
*
        WRITE(6,*) ' Supergroup in action '
        CALL IWRTMA(ISPGRP,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' NSMST = ' , NSMST
        write(6,*) 'First 3 columns of NSTSGP'
        CALL IWRTMA(NSTSGP,NSMST,3,NSMST,3)
      END IF
#endif
*
      CALL ISETVC(NSTSSPGP(1,IGRP),0,NSMST)
*. Max and Min for allowed symmetries
      DO IGAS = 1, NGAS
        MXSM(IGAS) = 1
        DO ISYM = 1, NSMST
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
        END DO
        MNSM(IGAS) = NSMST
        DO ISYM = NSMST,1, -1
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MNSM(IGAS) = ISYM
        END DO
      END DO
*. First symmetry combination
      DO IGAS = 1, NGAS
         ISM(IGAS) = MNSM(IGAS)
      END DO
*. Loop over symmetries in each gas space
      IFIRST = 1
 1000 CONTINUE
      IF(IFIRST.EQ.1) THEN
        CALL ISETVC(ISM,1,NGAS)
        IFIRST = 0
        NONEW = 0
      ELSE
C       NXTNUM3(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
        CALL NXTNUM3(ISM,NGAS,MNSM,MXSM,NONEW)
      END IF
      IF(NONEW.EQ.0) THEN
*. Symmetry of current combination and number of strings in this supergroup
        ISMSPGP = ISM(1)
        NST = NSTSGP(ISM(1),ISPGRP(1))
        DO JGRP = 2, NGAS
          CALL SYMCOM_REL(3,7,ISMSPGP,ISM(JGRP),ISMSPGPO)
          ISMSPGP = ISMSPGPO
          NST = NST * NSTSGP(ISM(JGRP),ISPGRP(JGRP))
        END DO
C?      WRITE(6,*) ' Symmetry of groups '
C?      CALL IWRTMA(ISM,1,NGAS,1,NGAS)
C?      WRITE(6,*) ' Symmetry of super group and number of strings '
C?      WRITE(6,*) ISMSPGP , NST
        NSTSSPGP(ISMSPGP,IGRP) = NSTSSPGP(ISMSPGP,IGRP) + NST
        GOTO 1000
      END IF
*
      MXNSTR = -1
      DO ISTRSM = 1, NSMST
        MXNSTR = MAX(MXNSTR,NSTSSPGP(ISTRSM,IGRP))
      END DO
*
      NTEST = 00
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
     &  ' Number of strings per symmetry for supergroup',IGRP
        CALL IWRTMA10(NSTSSPGP(1,IGRP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Largest number of strings of given sym ',MXNSTR
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE NSTRSO_GAS_REL(NEL,NORB1,NORB2,NORB3,
     &                          NELMN1,NELMX1,NELMN3,NELMX3,
     &                          IOC,NORB,NSTASO,ISTASO,
     &                          NOCTYP,NSMST,IOTYP,IPRNT)
*
* Number of strings per symmetry for group IOTYP
*
* Gas version, no check of type : set to 1
*
* Jeppe Olsen Winter of 1994
*
      IMPLICIT REAL*8           ( A-H,O-Z)
#include "ipoist8.inc"
C     DIMENSION IOC(*),NSTASO(NOCTYP,NSMST)
      DIMENSION IOC(*),NSTASO(NSMST,*),ISTASO(NSMST,*)
*
      CALL ISETVC(NSTASO(1,IOTYP),0,NSMST)
      NTESTL = 0
      NTEST = MAX(IPRNT,NTESTL)
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
* Loop over possible partitionings between RAS1,RAS2,RAS3
      DO 1001 IEL1 = NELMX1,NELMN1,-1
      DO 1003 IEL3 = NELMN3,NELMX3, 1
       IF(IEL1.GT. NORB1 ) GOTO 1001
       IF(IEL3.GT. NORB3 ) GOTO 1003
       IEL2 = NEL - IEL1-IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
* Loop over RAS 1 occupancies
  901  CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
            CALL ISTVC2(IOC(1),0,1,IEL1)
            IFRST1 = 0
           ELSE
             CALL NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
#if defined LUCI_DEBUG
         IF( NTEST .GE.500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
#endif
         IFRST2 = 1
         IFRST3 = 1
* Loop over RAS 2 occupancies
  902    CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              CALL ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
              IFRST2 = 0
             ELSE
               CALL NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
#if defined LUCI_DEBUG
           IF( NTEST .GE.500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
#endif
           IFRST3 = 1
* Loop over RAS 3 occupancies
  903      CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                CALL ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
                IFRST3 = 0
               ELSE
                 CALL NXTORD(IOC(IEL1+IEL2+1),
     &           IEL3,IORB3F,IORB3L,NONEW3)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
#if defined LUCI_DEBUG
             IF( NTEST .GE. 500) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
#endif
* Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
*. Symmetry of string
             ISYM = ISYMST_REL(IOC,NEL)
C                   ISYMST_REL(STRING,NEL)
*. occupation type of string
COLD         ITYP = IOCTP2(IOC,NEL,IOTYP)
C                   IOCTP2(STRING,NEL)
*
             NSTASO(ISYM,IOTYP) = NSTASO(ISYM,IOTYP)+ 1
*
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)
     &     GOTO 901
 1003 CONTINUE
 1001 CONTINUE
*
*. The corresponding offset
*
      DO ISM = 1, NSMST
        IF(ISM.EQ.1) THEN
          ISTASO(ISM,IOTYP) = 1
        ELSE
          ISTASO(ISM,IOTYP) = ISTASO(ISM-1,IOTYP)+NSTASO(ISM-1,IOTYP)
        END IF
      END DO
 
#if defined LUCI_DEBUG
      IF(NTEST.GE.5)
     &WRITE(6,*) ' Number of strings generated   ', NSTRIN
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Number of strings per sym for group = ', IOTYP 
        WRITE(6,*) '================================================'
        CALL IWRTMA(NSTASO(1,IOTYP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Offset for given symmetry for group = ', IOTYP 
        WRITE(6,*) '================================================'
        CALL IWRTMA(ISTASO(1,IOTYP),1,NSMST,1,NSMST)
      END IF
#endif
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE STRINF_GAS_REL(STIN,IPRNT)
*
* Obtain string information for GAS expansion
*
*--------------------------------------------
*   Jeppe Olsen
*
*   Only calling genstr_gas and crestr_gas
*   if number of electrons .gt. 0
*
*   Creation maps only if ISTTYPC.GE.1
*
*    Timo Fleig, Jan. 2001
*
* call to MAX_NSYMDIST_DBG added, Aug. 2003
*
*--------------------------------------------
*
* =====
*.Input
* =====
*
* /LUCINP/,/ORBINP/,/CSM/, /CGAS/, /GASSTR/
*
* =====
*.Output
* =====
*
* /STRINP/,/STINF/,/STRBAS/ and string information in STIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
*     (and /LUCINP/ not occuring here )
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "strbas.inc"
#include "csm.inc"
#include "cstate.inc"
#include "lucinp.inc"
#include "stinf.inc"
#include "strinp.inc"
#include "irat.inc"
*
#include "wrkspc.inc"
*
      dimension STIN(*)
*. A bit of scratch
C     DIMENSION IOCTYP(MXPNGAS)
*
C     CALL QENTER('STRIN')
*
      NTESTL = 000
      NTEST = MAX(NTESTL,IPRNT)
*
**.2 : Number of classes per string type and mappings between
**.    string types (/STINF/)
*
      CALL ZSTINF_GAS_REL(IPRNT)
*
**.3 : Static memory for string information
*
      CALL MEMSTR_GAS_REL
*
** 4 : Info about group of strings
*
*.First free address
      call memmar(KFREEL,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
#if defined LUCI_DEBUG
      if (NTEST.ge.1) then
        write(6,*) 'Getting first free address in strinf_gas'
        write(6,*) '----------------------------------------'
        write(6,*) 'first free address is   ',KFREEL
      end if
*
*. For spin up spinors :
      if (NTEST.ge.5) then
        write(6,*) '======================================='
        write(6,*)
        write(6,*) ' Info about UNBARRED spinor strings '
        write(6,*)
        write(6,*) '======================================='
      end if
#endif
      DO IGRP = 1, NGRP
*. A gas group can be considered as a RAS group with 0 electrons in
*  RAS1, RAS3 !
        IGAS = IGSFGP(IGRP)
        IF(IGAS.EQ.1) THEN
          NORB1 = 0
        ELSE
          NORB1 = IELSUM(NOBPT,IGAS-1)
        END IF
        NORB2 = NOBPT(IGAS)
        NORB3 = NACOB-NORB1-NORB2
        MNRS1X = 0
        MXRS1X = 0
        MNRS3X = 0
        MXRS3X = 0
        IEL = NELFGP(IGRP)
        IOCTYPX = 1
*. Reverse lexical adresing schemes for each group of string
        CALL WEIGHT_REL(STIN(KZ(IGRP)),IEL,NORB1,NORB2,NORB3,
     &                  MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                  STIN(KFREEL),IPRNT)
!       print *,'blub NSTRSO_GAS_REL'
*. Number of strings per symmetry in a given group
        CALL NSTRSO_GAS_REL(IEL,NORB1,NORB2,NORB3,
     &                      MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                      STIN(KFREEL),NACOB,
     &                      STIN(KNSTSGP(1)),
     &                      STIN(KISTSGP(1)),
     &                      IOCTYPX,NSMST,IGRP,IPRNT)
*
!       print *,'blub blub NSTRSO_GAS_REL'
       CALL ICOPVE2(STIN(KNSTSGP(1)),1+(IGRP-1)*NSMST,NSMST,
     &              NSTFSMGP(1,IGRP))
       CALL ICOPVE2(STIN(KISTSGP(1)),1+(IGRP-1)*NSMST,NSMST,
     &              ISTFSMGP(1,IGRP))
*. Construct the strings ordered by symmetry
        if (IEL.gt.0) then
!       print *,'blub blub GENSTR_GAS_REL',IEL
          CALL GENSTR_GAS_REL(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                        STIN(KISTSGP(1)),IGRP,
     &                        IOCTYPX,NSMST,STIN(KZ(IGRP)),STIN(KFREEL),
     &                        STIN(KSTREO(IGRP)),STIN(KOCSTR(IGRP)),
     &                        STIN(KFREEL+IOCTYPX*NSMST),IGRP,IPRNT)
        end if
!       print *,'blub blub blub GENSTR_GAS_REL'
      END DO
*
*. Min and max of sym for each group
*
#if defined LUCI_DEBUG
      if (NTEST.ge.50) then
        write(6,*) 'NSTFSMGP array:'
        call iwrtma(NSTFSMGP,NSMST,NGRP,MXPNSMST,MXPSTT)
      end if
#endif
*
      DO IGP = 1, NGRP
       MX = 1
       DO ISM = 1, NSMST
         IF(NSTFSMGP(ISM,IGP).GT.0) MX = ISM
       END DO
*
       MN = NSMST
       DO ISM = NSMST,1,-1
         IF(NSTFSMGP(ISM,IGP).GT.0) MN = ISM
       END DO
*
       MINMAX_SM_GP(1,IGP) = MN
       MINMAX_SM_GP(2,IGP) = MX
      END DO
*
#if defined LUCI_DEBUG
      if (NTEST.ge.30) THEN
        WRITE(6,*) ' MINMAX array for sym of groups '
        WRITE(6,*) ' =============================='
        CALL IWRTMA(MINMAX_SM_GP,2,NGRP,2,MXPSTT)
      end if
#endif
*. In double group CI, the spin down spinors may have different symmetry,
*  obtain the corresponding symmetries
      IF(PNTGRP.GE.5) THEN
*. Change symmetry array ISMFTO into ISMFTO2 (Dirty )
      DO IORB = 1, NTOOB
       ISWAP = ISMFTO(IORB)
       ISMFTO(IORB) = ISMFTO2(IORB)
       ISMFTO2(IORB) = ISWAP
      END DO
*
#if defined LUCI_DEBUG
      if (NTEST.ge.50) then
        WRITE(6,*) ' Updated ISMFTO array '
        CALL IWRTMA(ISMFTO,1,NTOOB,1,NTOOB)
      end if
*
      if (NTEST.ge.5) then
        write(6,*) '======================================='
        write(6,*)
        write(6,*) ' Info about BARRED spinor strings '
        write(6,*)
        write(6,*) '======================================='
      end if
#endif
      DO IGRP = 1, NGRP
*. A gas group can be considered as a RAS group with 0 electrons in
*  RAS1, RAS3 !
        IGAS = IGSFGP(IGRP)
        IF(IGAS.EQ.1) THEN
          NORB1 = 0
        ELSE
          NORB1 = IELSUM(NOBPT,IGAS-1)
        END IF
        NORB2 = NOBPT(IGAS)
        NORB3 = NACOB-NORB1-NORB2
        MNRS1X = 0
        MXRS1X = 0
        MNRS3X = 0
        MXRS3X = 0
        IEL = NELFGP(IGRP)
        IOCTYPX = 1
*. Reverse lexical adresing schemes for each group of string
        CALL WEIGHT_REL(STIN(KZ(IGRP)),IEL,NORB1,NORB2,NORB3,
     &                  MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                  STIN(KFREEL),IPRNT)
*. Number of strings per symmetry in a given group
        CALL NSTRSO_GAS_REL(IEL,NORB1,NORB2,NORB3,
     &                      MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                      STIN(KFREEL),NACOB,
     &                      STIN(KNSTSGP2(1)),
     &                      STIN(KISTSGP2(1)),
     &                      IOCTYPX,NSMST,IGRP,IPRNT)
*. Construct the strings ordered by symmetry
        CALL ICOPVE2(STIN(KNSTSGP2(1)),1+(IGRP-1)*NSMST,NSMST,
     &               NSTFSMGP2(1,IGRP))
        CALL ICOPVE2(STIN(KISTSGP2(1)),1+(IGRP-1)*NSMST,NSMST,
     &               ISTFSMGP2(1,IGRP))
        if (IEL.gt.0) then
          CALL GENSTR_GAS_REL(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                        STIN(KISTSGP2(1)),IGRP,
     &                        IOCTYPX,NSMST,STIN(KZ(IGRP)),STIN(KFREEL),
     &                        STIN(KSTREO2(IGRP)),STIN(KOCSTR2(IGRP)),
     &                        STIN(KFREEL+IOCTYPX*NSMST),IGRP,IPRNT)
        end if
      END DO
*
*. Min and max of sym for each group
*
#if defined LUCI_DEBUG
      if (NTEST.ge.50) then
        write(6,*) 'NSTFSMGP2 array:'
        call iwrtma(NSTFSMGP2,NSMST,NGRP,MXPNSMST,MXPSTT)
      end if
#endif
*
      DO IGP = 1, NGRP
       MX = 1
       DO ISM = 1, NSMST
         IF(NSTFSMGP2(ISM,IGP).GT.0) MX = ISM
       END DO
*
       MN = NSMST
       DO ISM = NSMST,1,-1
         IF(NSTFSMGP2(ISM,IGP).GT.0) MN = ISM
       END DO
*
       MINMAX_SM_GP2(1,IGP) = MN
       MINMAX_SM_GP2(2,IGP) = MX
      END DO
*
#if defined LUCI_DEBUG
      if (NTEST.ge.3) THEN
        WRITE(6,*) ' MINMAX array for sym of groups '
        WRITE(6,*) ' =============================='
        CALL IWRTMA(MINMAX_SM_GP2,2,NGRP,2,MXPSTT)
      end if
#endif

*. Restore ORDER !
      DO IORB = 1, NTOOB
       ISWAP = ISMFTO(IORB)
       ISMFTO(IORB) = ISMFTO2(IORB)
       ISMFTO2(IORB) = ISWAP
      END DO
      ELSE
*. Normal pointgroup, alpha and beta spin both totally symmetric
        DO IGRP = 1, NGRP
          KNSTSGP2(IGRP) = KNSTSGP(IGRP)
          KISTSGP2(IGRP) = KISTSGP(IGRP)
          KOCSTR2(IGRP) = KOCSTR(IGRP)
          KSTREO2(IGRP) = KSTREO(IGRP)
        END DO
      END IF
*     ^ End of switch double grop / point group
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Number of strings per group and symmetry '
        CALL IWRTMA10(STIN(KNSTSGP(1)),NSMST,NGRP,NSMST,NGRP)
*
        IF(PNTGRP.GE.5) THEN
        WRITE(6,*)
     &  ' Number of strings per group and symmetry (spindown)'
        CALL IWRTMA10(STIN(KNSTSGP2(1)),NSMST,NGRP,NSMST,NGRP)
        END IF
*
      END IF
#endif
*
*
* 4.5 : Creation/Annihilation mappings between different
*       types of strings (Annihilation mapping not installed)
*
      DO IGRP = 1, NGRP
        IF(ISTAC(IGRP,2).NE.0) THEN
          IGAS = IGSFGP(IGRP)
          NGSOBP = NOBPT(IGAS)
*. First orbital in GAS spacce
          IGSOB = IELSUM(NOBPT,IGAS-1)+1
          IEL = NELFGP(IGRP)
          NSTINI = NSTFGP(IGRP)
          CALL IZERO(STIN(KSTSTM(IGRP,1)),NGSOBP*NSTINI)
          CALL IZERO(STIN(KSTSTM(IGRP,2)),NGSOBP*NSTINI)
          JGRP = ISTAC(IGRP,2)
          CALL CRESTR_GAS_REL(STIN(KOCSTR(IGRP)),NSTFGP(IGRP),
     &         NSTFGP(JGRP),IEL,NGSOBP,IGSOB,STIN(KZ(JGRP)),
     &         STIN(KSTREO(JGRP)),0,IDUM,IDUM,
     &         STIN(KSTSTM(IGRP,1)),STIN(KSTSTM(IGRP,2)),NACOB,IPRNT)
        END IF
      END DO
*
      IF(PNTGRP.GE.5) THEN
*. Mappings for spin-down groups
      DO IGRP = 1, NGRP
        IF(ISTAC(IGRP,2).NE.0) THEN
          IGAS = IGSFGP(IGRP)
          NGSOBP = NOBPT(IGAS)
*. First orbital in GAS space
          IGSOB = IELSUM(NOBPT,IGAS-1)+1
          IEL = NELFGP(IGRP)
          NSTINI = NSTFGP(IGRP)
          CALL IZERO(STIN(KSTSTM2(IGRP,1)),NGSOBP*NSTINI)
          CALL IZERO(STIN(KSTSTM2(IGRP,2)),NGSOBP*NSTINI)
          JGRP = ISTAC(IGRP,2)
          CALL CRESTR_GAS_REL(STIN(KOCSTR2(IGRP)),NSTFGP(IGRP),
     &         NSTFGP(JGRP),IEL,NGSOBP,IGSOB,STIN(KZ(JGRP)),
     &         STIN(KSTREO2(JGRP)),0,IDUM,IDUM,
     &         STIN(KSTSTM2(IGRP,1)),STIN(KSTSTM2(IGRP,2)),NACOB,IPRNT)
        END IF
      END DO
      ELSE
*. Normal pointgroup calculation
         DO IGRP = 1, NGRP
          KSTSTM2(IGRP,1) = KSTSTM(IGRP,1)
          KSTSTM2(IGRP,2) = KSTSTM(IGRP,2)
         END DO
      END IF
*     ^ End of switch Double group/POintgroup

*
*. Now to string types, i.e. strings of with given number of elecs
*
      CALL IZERO(NSTFSMSPGP,MXPNSMST*NTSPGP)
      MXNSTR = -1
      DO ITP = 1, NSTTYP
*. Loop over supergroups of given type . i.e. strings
*  with given occupation in each GAS space
        DO IGRP = 1, NSPGPFTP(ITP)
          IGRPABS = IGRP-1 + IBSPGPFTP(ITP)
          CALL NSTPTP_GAS_REL(NGAS,ISPGPFTP(1,IGRPABS),
     &                        WORK(KNSTSGP(1)),NSMST,
     &                        WORK(KNSTSO(ITP)),IGRP,MXNSTRFSG)
          MXNSTR = MAX(MXNSTR,MXNSTRFSG)
        END DO
*
        CALL ICOPMT(WORK(KNSTSO(ITP)),NSMST,NSPGPFTP(ITP),
     &              NSTFSMSPGP(1,IBSPGPFTP(ITP)),MXPNSMST,NSPGPFTP(ITP))
        CALL ZSPGPIB(WORK(KNSTSO(ITP)),WORK(KISTSO(ITP)),
     &                NSPGPFTP(ITP),NSMST)
*
        IF(PNTGRP.GE.5) THEN
*. Double group CI, symmetry for spin down strings
        DO IGRP = 1, NSPGPFTP(ITP)
          IGRPABS = IGRP-1 + IBSPGPFTP(ITP)
          CALL NSTPTP_GAS_REL(NGAS,ISPGPFTP(1,IGRPABS),
     &                        WORK(KNSTSGP2(1)),NSMST,
     &                        WORK(KNSTSO2(ITP)),IGRP,MXNSTRFSG)
          MXNSTR = MAX(MXNSTR,MXNSTRFSG)
        END DO
*
        CALL ICOPMT(WORK(KNSTSO2(ITP)),NSMST,NSPGPFTP(ITP),
     &       NSTFSMSPGP2(1,IBSPGPFTP(ITP)),MXPNSMST,NSPGPFTP(ITP))
        CALL ZSPGPIB(WORK(KNSTSO2(ITP)),WORK(KISTSO2(ITP)),
     &                NSPGPFTP(ITP),NSMST)
        ELSE
          KISTSO2(ITP) = KISTSO(ITP)
          CALL ICOPMT(WORK(KNSTSO(ITP)),NSMST,NSPGPFTP(ITP),
     &         NSTFSMSPGP2(1,IBSPGPFTP(ITP)),MXPNSMST,NSPGPFTP(ITP))
        END IF
*       ^ End of double group/point group switch
#if defined LUCI_DEBUG
        IF(NTEST.GE.5) THEN
          WRITE(6,*)
     &    ' Number of strings per sym (row) and supergroup(column)',
     &    ' for type = ', ITP
          CALL IWRTMA(WORK(KNSTSO(ITP)),NSMST,NSPGPFTP(ITP),
     &                NSMST,NSPGPFTP(ITP))
          IF(PNTGRP.GE.5) THEN
          WRITE(6,*)
     &    ' Number of strings(spin down) per sym (row) ',
     &    ' and supergroup(column) for type = ', ITP
          CALL IWRTMA(WORK(KNSTSO2(ITP)),NSMST,NSPGPFTP(ITP),
     &                NSMST,NSPGPFTP(ITP))
          WRITE(6,*)
          END IF

        END IF
#endif
*
      END DO
*. Number of electrons in each AS for each supergroup
      CALL ZNELFSPGP_REL(IPRNT)
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' =============================='
        WRITE(6,*) ' Occupation of each supergroup '
        WRITE(6,*) ' =============================='
        WRITE(6,*)
        WRITE(6,*) ' Supergroup   Occupation '
        WRITE(6,*) ' =========================================='
        DO ISPGP = 1, NTSPGP
          WRITE(6,'(1H ,I6,6X,12(I4))')
     &    ISPGP, (NELFSPGP(IGS, ISPGP), IGS = 1, NGAS)
        END DO
      END IF
#endif
*. Largest number of strings belonging to given supergroup
*. Largest Occupation block for given supergroup and sym
      MAX_STR_OC_BLK = -1
      MAX_STR_SPGP = 0
      DO ISPGP = 1, NTSPGP
        NSTR = IELSUM(NSTFSMSPGP(1,ISPGP),NSMST)
        MAX_STR_SPGP = MAX(MAX_STR_SPGP,NSTR)
        NEL = IELSUM(NELFSPGP(1,ISPGP),NGAS)
        DO ISTSM = 1, NSMST
          MAX_STR_OC_BLK
     &  = MAX(MAX_STR_OC_BLK,NEL*NSTFSMSPGP(ISTSM,ISPGP))
          IF(PNTGRP.GE.5) THEN
            MAX_STR_OC_BLK
     &    = MAX(MAX_STR_OC_BLK,NEL*NSTFSMSPGP2(ISTSM,ISPGP))
          END IF
        END DO
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.2) THEN
      WRITE(6,*)
     & ' Largest number of strings of given supergroup        ',
     & MAX_STR_SPGP
      WRITE(6,*) ' Largest block of string occupations ',
     &              MAX_STR_OC_BLK
      WRITE(6,*)
     & ' Largest number of strings of given supergroup and sym', MXNSTR
      END IF
*
      if (NTEST.ge.50) then
        WRITE(6,*) ' Matrix NSTFSMSPGP '
        CALL IWRTMA(NSTFSMSPGP,NSMST,NTSPGP,MXPNSMST,NTSPGP)
      end if
#endif
*
* Largest number of symmetrydistributions for any given supergroup
C           FUNCTION MAX_NSYMDIST_DBG(NSPGRP,ISPGRP,NGAS,NSMST)
      MAX_NSYMDIST = MAX_NSYMDIST_DBG(NTSPGP,ISPGPFTP,NGAS,NSMST)
*
*
* Possible occupation classes
*
      CALL OCCLS_REL(2,NMXOCCLS,WORK(KIOCLS),NACTEL,NGAS,
     &               IGSOCC(1,1),IGSOCC(1,2))
*
* Maps creation/annihilation of given gas orb from given supergroup
* gives new supergroup.
      CALL IZERO(WORK(KSPGPCR),NGAS*NTSPGP)
      CALL IZERO(WORK(KSPGPAN),NGAS*NTSPGP)
*
      DO ISTTYP = 1,NSTTYP
*. Creation map from this type
*. Adding one electron corresponds to reduce the
* string type number with 2
*(very primitive, should be generalized)
        IF(NSPGPFTP(ISTTYP).GT.0) THEN
          ISTTYPC = ISTTYP - 2
          IF (ISTTYPC.GE.1) then
            if (NSPGPFTP(ISTTYPC).GT.0.AND.
     &          NELFTP_REF(ISTTYP).EQ.NELFTP_REF(ISTTYPC)) THEN
C             write(6,*)
C    &          ' creation map from ',ISTTYP,' to ',ISTTYPC
C             write(6,*) ' Bases IBSPGPFTP '
C             write(6,*)   IBSPGPFTP(ISTTYP),IBSPGPFTP(ISTTYPC)
              CALL SPGP_AC(NELFSPGP(1,1), NSPGPFTP(ISTTYP),
     &                     NELFSPGP(1,1),NSPGPFTP(ISTTYPC),
     &                     NGAS,MXPNGAS,2,WORK(KSPGPCR),
     &                     IBSPGPFTP(ISTTYP),IBSPGPFTP(ISTTYPC))
            end if
          END IF
*. Annihilation maps
          ISTTYPA = ISTTYP+2
          IF(ISTTYPA.LE.NSTTYP.AND.NSPGPFTP(ISTTYPA).GT.0.AND.
     &       NELFTP_REF(ISTTYP).EQ.NELFTP_REF(ISTTYPA)) THEN
C          write(6,*)
C    &       ' annihilation map from ',ISTTYP,' to ',ISTTYPA
             CALL SPGP_AC(NELFSPGP(1,1), NSPGPFTP(ISTTYP),
     &                    NELFSPGP(1,1),NSPGPFTP(ISTTYPA),
     &                    NGAS,MXPNGAS,1,WORK(KSPGPAN),
     &                    IBSPGPFTP(ISTTYP),IBSPGPFTP(ISTTYPA))
          END IF
        END IF
      END DO
*
C?    WRITE(6,*) ' Memory Check at end of STRINF_GAS  '
*
C     CALL QEXIT('STRIN')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE STRTYP_GAS_REL(IPRNT)
*
* Find groups of strings in each GA space
*
* Output : /GASSTR/
*
* Jeppe Olsen, Oct 1994
*              July 1997, modified for relativistic calculations
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cstate.inc"
#include "gasstr.inc"
#include "strinp.inc"
#include "stinf.inc"
*. Local scratch
      DIMENSION IOCTYP(MXPNGAS)
*
C     CALL QENTER('STRTYP_REL')
      NTESTL = 000
      NTEST = MAX(IPRNT,NTESTL)
*. As input NCISPC GAS spaces IGSOCCX are given.
* Obtain space that cantains all these as special cases
*
      DO IGAS = 1, NGAS
       MINI = IGSOCCX(IGAS,1,1)
       MAXI = IGSOCCX(IGAS,2,1)
       DO ICISPC = 2, NCISPC
        MINI = MIN(MINI,IGSOCCX(IGAS,1,ICISPC))
        MAXI = MAX(MAXI,IGSOCCX(IGAS,2,ICISPC))
       END DO
       IGSOCC(IGAS,1) = MINI
       IGSOCC(IGAS,2) = MAXI
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Compound GAS space : '
        WRITE(6,*) ' ====================='
        WRITE(6,'(A)')
        WRITE(6,'(A)') '         Min. occ    Max. occ '
        WRITE(6,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(6,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCC(IGAS,1),IGSOCC(IGAS,2)
        END DO
      END IF
#endif
*
*. Find min and max number of elecs in each subspace
*
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          MNGSOC(IGAS) = IGSOCC(IGAS,1)
          MXGSOC(IGAS) = IGSOCC(IGAS,2)
        ELSE
          MXGSOC(IGAS) = IGSOCC(IGAS,2)-IGSOCC(IGAS-1,1)
          MNGSOC(IGAS) = MAX(0,IGSOCC(IGAS,1)-IGSOCC(IGAS-1,2))
        END IF
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,'(A)') ' Min and Max occupation in each GAS space: '
        WRITE(6,'(A)') ' ========================================= '
        WRITE(6,*)
        DO IGAS = 1,  NGAS
          WRITE(6,'(A,I2,4X,2I3)')
     &    '  GAS',IGAS,MNGSOC(IGAS),MXGSOC(IGAS)
        END DO
       END IF
#endif
*.
*. Occupation classes corresponding to largest CI space
*
C     CALL QENTER('OCCLS')
      CALL OCCLS_REL(1,NOCCLS,IOCCLS,NACTEL,NGAS,
     &               IGSOCC(1,1),IGSOCC(1,2))
      NMXOCCLS = NOCCLS
C     CALL QEXIT('OCCLS')
*
*. Loop over MS2 spaces
*
*. Obtain the required range of electrons in each gasspace
      DO IMS2 = 1, NMS2VAL
*. Number of alpha. and beta electrons
        MS2 = MS2VAL(IMS2)
        NAEL = (MS2 + NACTEL ) / 2
        NBEL = (NACTEL - MS2 ) / 2
*
        IF(NAEL + NBEL .NE. NACTEL ) THEN
            WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
          WRITE(6,'(5I4)')   MS2,NACTEL,NAEL,NBEL
          WRITE(6,*)
     &    ' STOP : NUMBER OF ELECTRONS AND MULTIPLICITY INCONSISTENT '
            Call Abend2
     &      ( ' NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY ' )
        END IF
*
#if defined LUCI_DEBUG
        IF(NTEST.GE.5) THEN
          WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
          WRITE(6,'(5I6)')   MS2,NACTEL,NAEL,NBEL
        END IF
#endif
*. Number of electrons to be subtracted or added
        MAXSUB = 2
        MAXADD = 2
*. electrons are only added for systems that atleast have halffilled
*. shells
        IGRP = 0
        MXAL = NAEL
        MNAL = NAEL
        MXBL = NBEL
        MNBL = NBEL
        NORBL = NTOOB
        DO IGAS = 1, NGAS
*. occupation constraints 1
         MXA1 = MIN(MXGSOC(IGAS),2*NOBPT(IGAS),MXAL)
         MXB1 = MIN(MXGSOC(IGAS),2*NOBPT(IGAS),MXBL)
         MNA1 = MAX(0,MNGSOC(IGAS)-MXA1)
         MNB1 = MAX(0,MNGSOC(IGAS)-MXB1)
*. Additional checks can be made here
         MXA = MXA1
         MXB = MXB1
         MNA = MNA1
         MNB = MNB1
*
         MXAL = MXAL - MNA
         MNAL = MAX(0,MNAL-MXA)
         MXBL = MXBL - MNB
         MNBL = MAX(0,MNBL-MXB)
*
#if defined LUCI_DEBUG
         IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Occupation numbers for IGAS = ', IGAS
          WRITE(6,*) ' MXAL MNAL MXBL MNBL ',MXAL,MNAL,MXBL,MNBL
          WRITE(6,*) ' MXA MNA MXB MNB ',MXA,MNA,MXB,MNB
         END IF
#endif
*
         MNAB = MIN(MNA,MNB)
         MXAB = MAX(MXA,MXB)
         MNAB = MAX(0,MNAB-MAXSUB)
         IF(MNAB*2.GE.2*NOBPT(IGAS)) MXAB = MIN(MXAB + 2,2*NOBPT(IGAS))
*
#if defined LUCI_DEBUG
         IF(NTEST.GE.100) WRITE(6,*) ' MNAB,MXAB',MNAB,MXAB
#endif
         IF(IMS2.EQ.1) THEN
           MNELFGP(IGAS) = MNAB
           MXELFGP(IGAS) = MXAB
           MXELFGP(IGAS) = MIN(2*NOBPT(IGAS),MXELFGP(IGAS))
         ELSE
           MNELFGP(IGAS) = MIN(MNELFGP(IGAS),MNAB)
           MXELFGP(IGAS) = MAX(MXELFGP(IGAS),MXAB)
           MXELFGP(IGAS) = MIN(2*NOBPT(IGAS),MXELFGP(IGAS))
         END IF
        END DO
C       ^ End of loop over gasspaces
      END DO
C     ^ End of loop over IMS2 values
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,'(A)')
     &  '   Range of orbital occupations in each gasspace : '
        WRITE(6,*)
        WRITE(6,'(A)')
     &  '   Gas space     Min occ.     Max occ.  '
        WRITE(6,'(A)')
     &  '   ============================================='
        DO IGAS = 1, NGAS
          WRITE(6,'(3X,I3,8X,I3,10X,I3)')
     &    IGAS, MNELFGP(IGAS),MXELFGP(IGAS)
        END DO
      END IF
#endif
*. Total number of groups
      NGRP = 0
      DO IGAS = 1, NGAS
        NGRP = NGRP + MXELFGP(IGAS)-MNELFGP(IGAS)+1
      END DO
#if defined LUCI_DEBUG
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Total number of groups ', NGRP
      END IF
#endif
      IF(NGRP.GT.MXPSTT) THEN
        WRITE(6,*) ' Too many string groups '
        WRITE(6,*) ' Current limit ', MXPSTT
        WRITE(6,*) ' STOP : STRTYP_GAS, Too many string groups'
                     Call Abend2
     &               ( ' STRTYP_GAS, Too many string groups' )
      END IF
*. And then the groups
      LGRP = 0
      DO IGAS = 1, NGAS
        IBGPSTR(IGAS) = LGRP + 1
        NGPSTR(IGAS) = MXELFGP(IGAS)-MNELFGP(IGAS)+1
        IEL = MNELFGP(IGAS)
        DO IGRP = IBGPSTR(IGAS),IBGPSTR(IGAS)+NGPSTR(IGAS)
          IGSFGP(IGRP) = IGAS
          NELFGP(IGRP) = IEL
          IF(IEL.GT.NOBPT(IGAS)) THEN
            NSTFGP(IGRP) = 0
          ELSE
            NSTFGP(IGRP) = IBION(NOBPT(IGAS),IEL)
          END IF
          IEL = IEL + 1
        END DO
        LGRP = LGRP + NGPSTR(IGAS)
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,'(A)') ' Information about Groups of strings '
        WRITE(6,'(A)') ' =================================== '
        WRITE(6,*)
        WRITE(6,*) '     GAS  MNEL  MXEL IBGRP  NGRP'
        WRITE(6,*) '    ============================'
        DO IGAS = 1, NGAS
          WRITE(6,'(5(2X,I4))') IGAS,MNELFGP(IGAS),
     &          MXELFGP(IGAS),IBGPSTR(IGAS),NGPSTR(IGAS)
        END DO
        WRITE(6,'(A,I3)')
     &  ' Total number of groups generated ', NGRP
*
        WRITE(6,'(A)') ' Information about each string group '
        WRITE(6,'(A)') ' ===================================='
        WRITE(6,*)
        IITYPE = 0
        WRITE(6,'(A)') ' GROUP  GAS   NEL      NSTR '
        WRITE(6,'(A)') ' ==========================='
        DO IGRP = 1, NGRP
          IITYPE = IITYPE + 1
          WRITE(6,'(3(2X,I4),2X,I8)')
     &    IITYPE,IGSFGP(IGRP),NELFGP(IGRP),NSTFGP(IGRP)
        END DO
      END IF
#endif
*
*. Creation-annihilation connections between groups
*
      DO IGRP = 1, NGRP
        ISTAC(IGRP,1) = 0
        ISTAC(IGRP,2) = 0
        DO JGRP = 1, NGRP
          IF(IGSFGP(IGRP).EQ.IGSFGP(JGRP).AND.
     &       NELFGP(IGRP).EQ.NELFGP(JGRP)-1) ISTAC(IGRP,2) = JGRP
          IF(IGSFGP(IGRP).EQ.IGSFGP(JGRP).AND.
     &       NELFGP(IGRP).EQ.NELFGP(JGRP)+1) ISTAC(IGRP,1) = JGRP
        END DO
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================================'
        WRITE(6,*) ' Annihilation / Creation connections'
        WRITE(6,*) ' ======================================'
        WRITE(6,*)
        CALL IWRTMA(ISTAC,NGRP,2,MXPSTT,2)
      END IF
#endif
*.
*
* Range of nael and nbel for wave function types.
*. Types with two less electrons must be constructed for
*. Intermediate states
*. range of NAEL
      NAELMX  = (MS2VAL(1) + NACTEL)/2
      NAELMN  = (MS2VAL(1) + NACTEL)/2
      DO IMS2 = 2, NMS2VAL
        NAEL    = (MS2VAL(IMS2) + NACTEL)/2
        NAELMX = MAX(NAELMX,NAEL)
        NAELMN = MIN(NAELMN,NAEL)
      END DO
      NBELMX = NACTEL - NAELMN
      NBELMN = NACTEL - NAELMX
*. Range we should include
      NAELMN = MAX(0,NAELMN-2)
      NBELMN = MAX(0,NBELMN-2)
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Range of alpha electrons',NAELMN,' to ', NAELMX
        WRITE(6,*) ' Range of beta  electrons',NBELMN,' to ', NBELMX
      END IF
#endif
*
* Type 1 : NAEL electrons
*      2 : NBEL ELECTRONS
*      3 : NAEL -1 ELECTRONS
*      4 : NBEL -1 ELECTRONS
*      5 : NAEL -2 ELECTRONS
*      6 : NBEL -2 ELECTRONS
*
*. For each determinant type (given nalpha, nbeta) will we
*  introduce these strings, single and double annihilated strings.
*
* We will order the types as
*    Loop over det types
*      Loop over annihilations 0,1,2 ..
*         alpha type
*         beta type
*      End of loop over annihilations
*    End of loop over types
*
* This maximized symmetry with previous code
      MXANNI = 2
      ISTTYP = 0
      DO IDTTP = 1, NMS2VAL
        NAELRF  = (MS2VAL(IDTTP) + NACTEL)/2
        NBELRF  = NACTEL - NAELRF
        DO IANNI = 0, MXANNI
*. Alpha terms
          ISTTYP = ISTTYP + 1
          NELEC(ISTTYP) = NAELRF-IANNI
          NELFTP(ISTTYP) = NAELRF-IANNI
          NELFTP_REF(ISTTYP) = NAELRF
          IF(IANNI.EQ.0) IST_FOR_DT(1,IDTTP) = ISTTYP
          IF(IANNI.EQ.1) ITPM1(1,IDTTP) = ISTTYP
          IF(IANNI.EQ.2) ITPM2(1,IDTTP) = ISTTYP
*. Beta terms
          ISTTYP = ISTTYP + 1
          NELEC(ISTTYP) = NBELRF-IANNI
          NELFTP(ISTTYP) = NBELRF-IANNI
          NELFTP_REF(ISTTYP) = NBELRF
          IF(IANNI.EQ.0) IST_FOR_DT(2,IDTTP) = ISTTYP
          IF(IANNI.EQ.1) ITPM1(2,IDTTP) = ISTTYP
          IF(IANNI.EQ.2) ITPM2(2,IDTTP) = ISTTYP
        END DO
      END DO
      NSTTYP = ISTTYP
      NSTTP = ISTTYP
*
      DO ITP = 1, NSTTYP
        NOCTYP(ITP) = 0
        NSPGPFTP(ITP) = 0
      END DO
*
* Loop over types, i.e.  given number of electrons
*
      IOFF = 1
      NABEL = NACTEL
      NSPGP_TOT = 0
      DO 2000 ITYP = 1, NSTTYP
*. Number of electrons in reference space ( alpha or beta )
        NELEC_REF = NELFTP_REF(ITYP)
        IDEL = NELEC(ITYP) - NELEC_REF
*. If we are studying beta type, and number of alpha and beta
* electrons are identical, just refer to alpha
C       IF(NAEL.EQ.NBEL.AND.MOD(ITYP,2).EQ.0) THEN
C         IBSPGPFTP(ITYP) =  IBSPGPFTP(ITYP-1)
C         NOCTYP(ITYP) =   NOCTYP(ITYP-1)
C         NSPGPFTP(ITYP) =  NSPGPFTP(ITYP-1)
*. Commented out, for double groups we can have identical
*. number of elecs but different symmetries so we play it safe.
C       ELSE
*. Number of electrons removed compared to reference
C?      WRITE(6,*) '  GASSPC : ITYP IDEL ', ITYP,IDEL
*. Generate all supergroups that belongs to given type of string.
*. Initial type of strings, relative to offset for given group
        DO IGAS = 1, NGAS
          IOCTYP(IGAS) = 1
        END DO
        NSPGP = 0
        IBSPGPFTP(ITYP) = IOFF
        IF(NELEC(ITYP).LT.0) THEN
          NOCTYP(ITYP) = 0
          NSPGPFTP(ITYP) =  0
          GOTO 2000
        END IF
*. Number of electrons in present type
*. Loop over  SUPER GROUPS with current nomenclature!
*. Temp max for loop
        MXLOOP = 10000
        IONE = 1
        NLOOP = 0
 1000   CONTINUE
*. Number of electrons in present supergroup
          NEL = 0
          DO IGAS = 1, NGAS
            NEL = NEL + NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
          END DO
*
          IF(NEL.GT.NELEC(ITYP)) THEN
*. If the number of electrons is to large find next number that
* can be correct.
* The following uses that within a given GAS space
* the number of elecs increases as the type number increases
*
*. First integer  that can be reduced
            IRED = 0
            DO IGAS = 1, NGAS
              IF(IOCTYP(IGAS).NE.1) THEN
                IRED = IGAS
                GOTO 888
              END IF
            END DO
  888       CONTINUE
            IF(IRED.EQ.NGAS) THEN
              NONEW = 1
            ELSE IF(IRED.LT.NGAS) THEN
              IOCTYP(IRED) = 1
*. Increase remanining part
              CALL NXTNUM2(
     &        IOCTYP(IRED+1),NGAS-IRED,IONE,NGPSTR(IRED+1),NONEW)
            END IF
            GOTO 2803
          END IF

          IF(NEL.EQ.NELEC(ITYP)) THEN
*. test 1 has been passed, check additional occupation constraints
*
*. Check from above
           I_AM_OKAY = 1
           DO IGAS = NGAS, 1, -1
*. Number of electrons when all electrons of AS IGAS have been added
             IF(IGAS.EQ.NGAS ) THEN
               IEL = MAX(NABEL,NABEL+IDEL)
             ELSE
               IEL = IEL-NELFGP(IOCTYP(IGAS+1)+IBGPSTR(IGAS+1)-1)
               IF(IEL.LT.MAX(IGSOCC(IGAS,1),IGSOCC(IGAS,1)+IDEL))
     &         I_AM_OKAY = 0
             END IF
           END DO
*
           IF(I_AM_OKAY.EQ.1) THEN
*. passed !!!
             NSPGP = NSPGP + 1
*. Copy supergroup to ISPGPFTP with absolute group numbers
             DO IGAS = 1, NGAS
               ISPGPFTP(IGAS,IOFF-1+NSPGP)
     &       = IOCTYP(IGAS)+IBGPSTR(IGAS)-1
             END DO
           END IF
*
          END IF
*. Next type of strings
          IONE = 1
          CALL NXTNUM2(IOCTYP,NGAS,IONE,NGPSTR,NONEW)
 2803   CONTINUE
        IF(NONEW.EQ.0) GOTO 1000
*. End of loop over possible supergroups, save information about current type
        IOFF = IOFF + NSPGP
        NOCTYP(ITYP) = NSPGP
        NSPGPFTP(ITYP) =  NSPGP
        NSPGP_TOT =  NSPGP_TOT +  NSPGP
C     END IF
 2000 CONTINUE
*     ^ End of loop over types (#s of electrons)
*
       NTSPGP = NSPGP_TOT
       IF(NSPGP_TOT .GT. MXPSTT ) THEN
         WRITE(6,*) ' Too many super groups = ', NSPGP_TOT
         WRITE(6,*) ' Increase MXPSTT to this value'
         WRITE(6,*) ' See you later '
         WRITE(6,*)
         WRITE(6,*) ' STOP Increase MXPSTT '
         Call Abend2( ' Increase MXPSTT' )
       END IF
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.2) THEN
       WRITE(6,*) ' Total number of super groups ', NTSPGP
       WRITE(6,*) ' Number of alpha supergroups  ', NSPGPFTP(1)
       WRITE(6,*) ' Number of beta  supergroups  ', NSPGPFTP(2)
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Information about types of strings'
        WRITE(6,*) ' =================================='
        WRITE(6,*)
        DO ITYP = 1, NSTTYP
          WRITE(6,*)
          WRITE(6,*) '      Type : ', ITYP
          WRITE(6,*) '      ==============='
          WRITE(6,*) '      Number of electrons  ',NELFTP(ITYP)
          WRITE(6,*) '      Number of super groups ', NSPGPFTP(ITYP)
          WRITE(6,*) '      Supergroups '
          DO ISPGP = 1, NSPGPFTP(ITYP)
            IOFF = IBSPGPFTP(ITYP)
            CALL IWRTMA(ISPGPFTP(1,IOFF-1+ISPGP),1,NGAS,1,NGAS)
          END DO
        END DO
        write(6,*) 'NOCTYP vector:'
        CALL IWRTMA(NOCTYP,1,NSTTYP,1,NSTTYP)
      END IF
#endif
*
C     CALL QEXIT('STRTYP_REL')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE TS_SYM_PNT2_REL(IGRP,NIGRP,MAXVAL,MINVAL,ISYM,
     &                           IPNT,LPNT,IAB)
*
* Construct pointers to start of symmetry distributions
* for supergroup of strings with given symmetry
*
* The start of symmetry block ISYM1 ISYM2 ISYM3 .... ISYMN
* is given as
*     1
*     + (ISM1-MINVAL(1))
*     + (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
*     + (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     + (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)
*
* Where L is the last group of strings with nonvanishing occupation
*
* Jeppe Olsen, September 1997
* Same kind of bug fix inserted as in getstr_totsm_spgp. (July 2001)
* IAB inserted for relativistic calculations. (August 2001)
*   Timo Fleig
*
* Version 2 : Uses IGRP and NIGRP to define supergroup
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "csm.inc"
*. Specific Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER NELFGS(MXPNGAS),ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
*. Output
      INTEGER MINVAL(*),MAXVAL(*),IPNT(*)
*
      NTEST = 000
*. Info on groups of strings in supergroup
      NGASL = 1
      DO IGAS = 1, NIGRP
       ITPFGS(IGAS) = IGRP(IGAS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
*. Number of strings per symmetry in each gasspace
        if (IAB.eq.1) then
          CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,IGAS))
        else if (IAB.eq.2) then
          CALL ICOPVE2(WORK(KNSTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,IGAS))
        end if
      END DO
*
      IZERO = 0
      DO IGAS = 1, NIGRP
        MINVAL(IGAS) = IZERO
        MAXVAL(IGAS) = IZERO
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
CTF
* If no strings are found in a given GAS, we set the symmetry
* value in MINVAL to 1, just to prevent crashes (only a test,
* waiting for more elegant solution ...)
      IONE = 1
      do IGAS = 1,NGAS-1,1
        if (MINVAL(IGAS).eq.0) then
          MINVAL(IGAS) = IONE
        end if
      end do
CTF
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' MINVAL and MAXVAL '
        CALL IWRTMA(MINVAL,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MAXVAL,1,NIGRP,1,NIGRP)
        WRITE(6,*) ' NIGRP = ', NIGRP
      END IF
#endif
*. Total number of strings that will be generated
      NBLKS = 1
      DO IGAS = 1, NGASL-1
       NBLKS = NBLKS*(MAXVAL(IGAS)-MINVAL(IGAS)+1)
      END DO
      IF(NBLKS.GT.LPNT) THEN
        WRITE(6,*) ' Problem in TS_SYM_PNT'
        WRITE(6,*) ' Dimension of IPNT too small'
        WRITE(6,*) ' Actual and required length',NBLKS,LPNT
        WRITE(6,*)
        WRITE(6,*) ' I will Stop and wait for instructions'
        Call Abend2( ' TS_SYM_PNT too small' )
      END IF
*. Loop over symmetry blocks in standard order
      IFIRST = 1
      ISTRBS = 1
      NSTRINT = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 2001
        END IF
        IFIRST = 0
*. Symmetry of NGASL -1 spaces given, symmetry of full space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL SYMCOM_REL(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*.  sym of SPACE NGASL
        CALL SYMCOM_REL(2,1,ISTSMM1,ISMGSN,ISYM)
        ISMFGS(NGASL) = ISMGSN
#if defined LUCI_DEBUG
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
#endif
*. Number of strings with this symmetry combination
        NSTRII = 1
        DO IGAS = 1, NGASL
          NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for this symmetry distribution in IOFFI
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGASL-1
          IOFF = IOFF + (ISMFGS(IGAS)-MINVAL(IGAS))*MULT
          MULT = MULT * (MAXVAL(IGAS)-MINVAL(IGAS)+1)
        END DO
*
        IPNT(IOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
#if defined LUCI_DEBUG
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IOFF, IPNT(IOFF) NSTRII ',
     &                 IOFF, IPNT(IOFF),NSTRII
        END IF
#endif
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from TS_SYM_PNT'
        WRITE(6,*) ' Required total symmetry',ISYM
        WRITE(6,*) ' Number of symmetry blocks ', NBLKS
        WRITE(6,*)
        WRITE(6,*) ' Offset array  for symmetry blocks'
        CALL IWRTMA(IPNT,1,NBLKS,1,NBLKS)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE WEIGHT_REL(Z,NEL,NORB1,NORB2,NORB3,
     &                      MNRS1,MXRS1,MNRS3,MXRS3,ISCR,NTEST)
*
* construct vertex weights
*
* Reverse lexical ordering is used for restricted space
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
      INTEGER Z(*),ISCR(*)
*
      NORB = NORB1 + NORB2 + NORB3
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' >>>> WEIGHT <<<<< '
        WRITE(6,*) ' NORB1 NORB2 NORB3 ',NORB1,NORB2,NORB3
        WRITE(6,*) ' NEL MNRS1 MXRS1 MNRS3 MXRS3 '
        WRITE(6,*)   NEL,MNRS1,MXRS1,MNRS3,MXRS3
      END IF
#endif
*
      KLFREE = 1
      KLMAX = KLFREE
      KLFREE = KLFREE + NORB
*
      KLMIN = KLFREE
      KLFREE = KLFREE + NORB
*
      KW = KLFREE
      KLFREE = KW + (NEL+1)*(NORB+1)
*.Max and min arrays for strings
      CALL RSMXMN_REL(ISCR(KLMAX),ISCR(KLMIN),NORB1,NORB2,NORB3,
     &                NEL,MNRS1,MXRS1,MNRS3,MXRS3,NTEST)
*. Arc weights
      CALL GRAPW(ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NEL,NTEST)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ZSTINF_GAS_REL(IPRNT)
*
* Set up common block /STINF/ from information in /STINP/
*
*=========
* Input
*=========
* Information in /CGAS/ and /GASSTR/
*
*======================
* Output ( in /STINF/ )
*======================
* ISTAC (MXPSTT,2) : string type obtained by creating (ISTAC(ITYP,2))
*                    or annihilating (ISTAC(ITYP,1)) an electron
*                    from a string of type  ITYP . A zero indicates
*                    that this mapping is not included
*                    Only strings belonging to the same
*                    Orbital group are mapped
*                    mapped
*. Input
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
C     COMMON/CGAS/IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),
C    &            NGSOB(MXPOBS,MXPNGAS),
C    &            NGSOBT(MXPNGAS),IGSOCC(2,MXPNGAS),IGSINA,IGSDEL
C     COMMON/GASSTR/MNGSOC(MXPNGAS),MXGSOC(MXPNGAS),NGPSTR(MXPNGAS),
C    &              IBGPSTR(MXPNGAS),NELFGP(MXPSTT),IGSFGP(MXPSTT),
C    &              NSTFGP(MXPSTT),MNELFGP(MXPNGAS),MXELFGP(MXPNGAS)
*. Output
#include "stinf.inc"
*./STINF/
C     COMMON/STINF/ISTAC(MXPSTT,2),NOCTYP(MXPSTT),NSTFTP(MXPSTT),
C    &             INUMAP(MXPSTT),INDMAP(MXPSTT)
*. Only the first element, i.e. ISTAC  is defined

*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
* ******************************************************************
* Mappings between strings with the same type ISTTP index , +/- 1 el
* ******************************************************************
      CALL IZERO(ISTAC,2*MXPSTT)
      DO  IGAS = 1, NGAS
*. groups for a given gas spaces goes with increasing number of orbitals,
*  so the first space does not have any creation mapping
*  and the last space does not have any annihilation mapping
*
        MGRP = NGPSTR(IGAS)
        DO IGRP = 1, MGRP
          IIGRP = IGRP + IBGPSTR(IGAS) -1
*. Annihilation map is present : IIGRP => IIGRP - 1
          IF(IGRP.NE.1) ISTAC(IIGRP,1) = IIGRP -1
*. Creation map is present : IIGRP => IIGRP + 1
          IF(IGRP.NE.MGRP) ISTAC(IIGRP,2) = IIGRP + 1
        END DO
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Type - type mapping array ISTAC '
        WRITE(6,*) ' =============================== '
        CALL IWRTMA(ISTAC,NGRP  ,2,MXPSTT,2)
      END IF
#endif
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION MAX_NSYMDIST_DBG(NSPGRP,ISPGRP,NGAS,NSMST)
*
* Obtain largest number of symmetrydistributions of 
* any given supergroup.
*
* Timo and Jeppe, Odense, August 2003 
* ( to avoid that memmar takes 78 % of all time in a test ....)
*
*. Input
* NSPGRP : Number of supergroups
* ISPGRP : The groups of the various supergroups
* NGAS   : Number of GAS spaces
* NSTFSMGP : Number of strings of given sym for given group
* NSMST    : Number of string-symmetries 
*
*
#include "implicit.inc"
*. Input
      INTEGER ISPGRP(NGAS,NSPGRP)
*
      MAX_NSYMDIST_DBG = 0
      MAXX = -1
      NJDIST = -1
      DO IAB = 1, 2
       DO JSPGRP = 1, NSPGRP
        CALL NST_SPGRP_DBG(NGAS,ISPGRP(1,JSPGRP),1,NSMST,
     &                     NJSTR,NJDIST,IAB)
C       call nst_spgrp_dbg(NIGRP,KGRP,KSM,NSMST,NKSTR,NKDIST,IAB)
        MAXX = MAX(MAXX,NJDIST)
       END DO
      END DO
*
      MAX_NSYMDIST_DBG = MAXX
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Largest number of symmetrydistributions ',
     &  MAXX
      END IF
#endif
      END
