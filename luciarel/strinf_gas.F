      SUBROUTINE STRINF_GAS(STIN,IPRNT)
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
      CALL QENTER('STRIN')
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
*
**.2 : Number of classes per string type and mappings between
**.    string types (/STINF/)
*
      CALL ZSTINF_GAS(IPRNT)
*
**.3 : Static memory for string information
*
      CALL MEMSTR_GAS
*
** 4 : Info about group of strings
*
*.First free address
      write(6,*) 'Getting first free address in strinf_gas'
      write(6,*) '----------------------------------------'
      CALL MEMMAN(KFREEL,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      write(6,*) 'first free address is   ',KFREEL
*. For spin up spinors :
      if (NTEST.ge.5) then
        write(6,*) '======================================='
        write(6,*)
        write(6,*) ' Info about UNBARRED spinor strings '
        write(6,*)
        write(6,*) '======================================='
      end if
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
        CALL WEIGHT(STIN(KZ(IGRP)),IEL,NORB1,NORB2,NORB3,
     &              MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &              STIN(KFREEL),IPRNT )
*. Number of strings per symmetry in a given group
        CALL NSTRSO_GAS(IEL,NORB1,NORB2,NORB3,
     &              MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &              STIN(KFREEL),NACOB,
     &              STIN(KNSTSGP(1)),
     &              STIN(KISTSGP(1)),
     &              IOCTYPX,NSMST,IGRP,IPRNT)
*
       CALL ICOPVE2(STIN(KNSTSGP(1)),1+(IGRP-1)*NSMST,NSMST,
     &              NSTFSMGP(1,IGRP))
       CALL ICOPVE2(STIN(KISTSGP(1)),1+(IGRP-1)*NSMST,NSMST,
     &              ISTFSMGP(1,IGRP))
*. Construct the strings ordered by symmetry
        if (IEL.gt.0) then
          CALL GENSTR_GAS(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                STIN(KISTSGP(1)),IGRP,
     &                IOCTYPX,NSMST,STIN(KZ(IGRP)),STIN(KFREEL),
     &                STIN(KSTREO(IGRP)),STIN(KOCSTR(IGRP)),
     &                STIN(KFREEL+IOCTYPX*NSMST),IGRP,IPRNT)
        end if
      END DO
*
*. Min and max of sym for each group
*
C?    write(6,*) 'NSTFSMGP array:'
C?    call iwrtma(NSTFSMGP,NSMST,NGRP,MXPNSMST,MXPSTT)
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
      if (NTEST.ge.30) THEN
        WRITE(6,*) ' MINMAX array for sym of groups '
        WRITE(6,*) ' =============================='
        CALL IWRTMA(MINMAX_SM_GP,2,NGRP,2,MXPSTT)
      end if

*. In double group CI, the spin down spinors may have different symmetry,
*  obtain the corresponding symmetries
      IF(PNTGRP.GE.5) THEN
*. Change symmetry array ISMFTO into ISMFTO2 (Dirty )
      DO IORB = 1, NTOOB
       ISWAP = ISMFTO(IORB)
       ISMFTO(IORB) = ISMFTO2(IORB)
       ISMFTO2(IORB) = ISWAP
      END DO
C     WRITE(6,*) ' Updated ISMFTO array '
C     CALL IWRTMA(ISMFTO,1,NTOOB,1,NTOOB)
*
      if (NTEST.ge.5) then
        write(6,*) '======================================='
        write(6,*)
        write(6,*) ' Info about BARRED spinor strings '
        write(6,*)
        write(6,*) '======================================='
      end if
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
        CALL WEIGHT(STIN(KZ(IGRP)),IEL,NORB1,NORB2,NORB3,
     &              MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &              STIN(KFREEL),IPRNT )
*. Number of strings per symmetry in a given group
        CALL NSTRSO_GAS(IEL,NORB1,NORB2,NORB3,
     &              MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &              STIN(KFREEL),NACOB,
     &              STIN(KNSTSGP2(1)),
     &              STIN(KISTSGP2(1)),
     &              IOCTYPX,NSMST,IGRP,IPRNT)
*. Construct the strings ordered by symmetry
        CALL ICOPVE2(STIN(KNSTSGP2(1)),1+(IGRP-1)*NSMST,NSMST,
     &               NSTFSMGP2(1,IGRP))
        CALL ICOPVE2(STIN(KISTSGP2(1)),1+(IGRP-1)*NSMST,NSMST,
     &               ISTFSMGP2(1,IGRP))
        if (IEL.gt.0) then
          CALL GENSTR_GAS(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                STIN(KISTSGP2(1)),IGRP,
     &                IOCTYPX,NSMST,STIN(KZ(IGRP)),STIN(KFREEL),
     &                STIN(KSTREO2(IGRP)),STIN(KOCSTR2(IGRP)),
     &                STIN(KFREEL+IOCTYPX*NSMST),IGRP,IPRNT)
        end if
      END DO
*
*. Min and max of sym for each group
*
      if (NTEST.ge.5) then
        write(6,*) 'NSTFSMGP2 array:'
        call iwrtma(NSTFSMGP2,NSMST,NGRP,MXPNSMST,MXPSTT)
      end if
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
      if (NTEST.ge.3) THEN
        WRITE(6,*) ' MINMAX array for sym of groups '
        WRITE(6,*) ' =============================='
        CALL IWRTMA(MINMAX_SM_GP2,2,NGRP,2,MXPSTT)
      end if

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
          IZERO = 0
          CALL ISETVC(STIN(KSTSTM(IGRP,1)),IZERO,NGSOBP*NSTINI)
          CALL ISETVC(STIN(KSTSTM(IGRP,2)),IZERO,NGSOBP*NSTINI)
          JGRP = ISTAC(IGRP,2)
          CALL CRESTR_GAS(STIN(KOCSTR(IGRP)),NSTFGP(IGRP),
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
          IZERO = 0
          CALL ISETVC(STIN(KSTSTM2(IGRP,1)),IZERO,NGSOBP*NSTINI)
          CALL ISETVC(STIN(KSTSTM2(IGRP,2)),IZERO,NGSOBP*NSTINI)
          JGRP = ISTAC(IGRP,2)
          CALL CRESTR_GAS(STIN(KOCSTR2(IGRP)),NSTFGP(IGRP),
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
      CALL ISETVC(NSTFSMSPGP,0,MXPNSMST*NTSPGP)
      MXNSTR = -1
      DO ITP = 1, NSTTYP
*. Loop over supergroups of given type . i.e. strings
*  with given occupation in each GAS space
        DO IGRP = 1, NSPGPFTP(ITP)
          IGRPABS = IGRP-1 + IBSPGPFTP(ITP)
          CALL NSTPTP_GAS(NGAS,ISPGPFTP(1,IGRPABS),
     &                    WORK(KNSTSGP(1)),NSMST,
     &                    WORK(KNSTSO(ITP)),IGRP,MXNSTRFSG)
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
          CALL NSTPTP_GAS(NGAS,ISPGPFTP(1,IGRPABS),
     &                    WORK(KNSTSGP2(1)),NSMST,
     &                    WORK(KNSTSO2(ITP)),IGRP,MXNSTRFSG)
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
*
      END DO
*. Number of electrons in each AS for each supergroup
      CALL ZNELFSPGP(IPRNT)
*
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
C?    WRITE(6,*) ' Matrix NSTFSMSPGP '
C?    CALL IWRTMA(NSTFSMSPGP,NSMST,NTSPGP,MXPNSMST,NTSPGP)
*
*
* Possible occupation classes
*
      CALL OCCLS(2,NMXOCCLS,WORK(KIOCLS),NACTEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2))
*
* Maps creation/annihilation of given gas orb from given supergroup
* gives new supergroup.
      IZERO = 0
      CALL ISETVC(WORK(KSPGPCR),IZERO,NGAS*NTSPGP)
      CALL ISETVC(WORK(KSPGPAN),IZERO,NGAS*NTSPGP)
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
      WRITE(6,*) ' Memory Check at end of STRINF_GAS  '
      CALL LMEMCHK
*
      CALL QEXIT('STRIN')
      RETURN
      END
