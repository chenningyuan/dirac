      SUBROUTINE STRTYP_GAS_KRCC
*
* Find groups of strings in each GA space
*
* Output : /GASSTR/
*
* Jeppe Olsen, Oct 1994
*              July 1997, modified for relativistic calculations
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "crun.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cstate.inc"
#include "gasstr.inc"
#include "strinp.inc"
#include "stinf.inc"
*. Local scratch
      DIMENSION IOCTYP(MXPNGAS)
*
      CALL QENTER('STRTY')
      NTEST = 00
*. As input NCISPC GAS spaces IGSOCCX are given.
* Obtain space that contains all these as special cases
*
* Next line added by Lasse for several calc input
      NCISPC = NCMBSPC
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
*
* Lasse addition start
* Notice IUSE_PH harwired and not on common block
*. Particle and hole spaces  :
*  Hole spaces are always more than half occupied
*
      IPHGASL = 0
      NPHGAS = 0
      IUSE_PH = 1
      DO IGAS = 1, NGAS
        IF(IUSE_PH.EQ.1) THEN
*. P/H separation for compound space
C         IF(  MNGSOC(IGAS).GT.NOBPT(IGAS)) THEN
          IF(  MNGSOC(IGAS).GE.NOBPT(IGAS)) THEN
            IPHGAS(IGAS) = 2
            IPHGASL = IGAS
            NPHGAS = NPHGAS + 1
          ELSE
             IPHGAS(IGAS) = 1
          END IF
*. P/H separation for initial  space
          IF(IGAS.EQ.1) THEN
            MIN_OC1 = IGSOCCX(IGAS,1,1)
          ELSE
            MIN_OC1 = IGSOCCX(IGAS,1,1)-IGSOCCX(IGAS-1,2,1)
          END IF
          IF(MIN_OC1.GT.NOBPT(IGAS)) THEN
            IPHGAS1(IGAS) = 2
          ELSE
             IPHGAS1(IGAS) = 1
          END IF
        ELSE IF(IUSE_PH.EQ.0) THEN
          IPHGAS(IGAS) = 1
          IPHGAS1(IGAS) = 1
        END IF
      END DO
*. Large number of particle and hole orbitals of given type
      MXTSOB_P = 0
      MXTSOB_H = 0
      DO IGAS = 1, NGAS
        IF(IPHGAS1(IGAS).EQ.1) THEN
          MXTSOB_P = MAX(MXTSOB_P,NOBPT(IGAS))
        ELSE
          MXTSOB_H = MAX(MXTSOB_H,NOBPT(IGAS))
        END IF
      END DO
*
* Lasse addition end
*
*.
*. Occupation classes corresponding to largest CI space
*
      CALL QENTER('OCCLS')
      CALL OCCLS_KRCC(1,NOCCLS,IOCCLS,NACTEL,NGAS,
     &               IGSOCC(1,1),IGSOCC(1,2))
      NMXOCCLS = NOCCLS
      CALL QEXIT('OCCLS')
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
        IF(NTEST.GE.5) THEN
          WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
          WRITE(6,'(5I6)')   MS2,NACTEL,NAEL,NBEL
        END IF
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
         IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Occupation numbers for IGAS = ', IGAS
          WRITE(6,*) ' MXAL MNAL MXBL MNBL ',MXAL,MNAL,MXBL,MNBL
          WRITE(6,*) ' MXA MNA MXB MNB ',MXA,MNA,MXB,MNB
         END IF
*
         MNAB = MIN(MNA,MNB)
         MXAB = MAX(MXA,MXB)
         MNAB = MAX(0,MNAB-MAXSUB)
         IF(MNAB*2.GE.2*NOBPT(IGAS)) MXAB = MIN(MXAB + 2,2*NOBPT(IGAS))
*
         IF(NTEST.GE.100) WRITE(6,*) ' MNAB,MXAB',MNAB,MXAB
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
        print*,'now max'
        DO IGAS = 1, NGAS
          WRITE(6,'(3X,I3,8X,I3,10X,I3)')
     &    IGAS, MXELFGP(IGAS),MXELFGP(IGAS)
        END DO
      END IF
*. Total number of groups
      NGRP = 0
      DO IGAS = 1, NGAS
        NGRP = NGRP + MXELFGP(IGAS)-MNELFGP(IGAS)+1
      END DO
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Total number of groups ', NGRP
      END IF
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
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================================'
        WRITE(6,*) ' Annihilation / Creation connections'
        WRITE(6,*) ' ======================================'
        WRITE(6,*)
        CALL IWRTMA(ISTAC,NGRP,2,MXPSTT,2)
      END IF
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
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Range of alpha electrons',NAELMN,' to ', NAELMX
        WRITE(6,*) ' Range of beta  electrons',NBELMN,' to ', NBELMX
      END IF
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
      IF(NTEST.GE.2) THEN
       WRITE(6,*) ' Total number of super groups ', NTSPGP
       WRITE(6,*) ' Number of alpha supergroups  ', NSPGPFTP(1)
       WRITE(6,*) ' Number of beta  supergroups  ', NSPGPFTP(2)
      END IF
*
      WRITE(6,*)
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
*
      CALL QEXIT('STRTY')
*
      RETURN
      END
*
      SUBROUTINE STRINF_GAS_KRCC(WORK,KFREE,LFREE)
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
*
      dimension WORK(*)
C
C#include "memint.h"
C
*. A bit of scratch
C     DIMENSION IOCTYP(MXPNGAS)
*
C     CALL MEMGET('INTE',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
C     CALL MEMGET('REAL',KDAX,NORBT*NORBT*NZ,WORK,KFREE,LFREE)
C     SUBROUTINE MEMGET(DATA_TYPE,KBASE,LENGTH,WORK,KFREE,LFREE)
      CALL QENTER('STRIN')
*
      NTEST = 00
      IPRNT = NTEST
*
**.2 : Number of classes per string type and mappings between
**.    string types (/STINF/)
*
      CALL ZSTINF_GAS_KRCC(IPRNT)
*
**.3 : Static memory for string information
*
      CALL MEMSTR_GAS_KRCC(WORK,KFREE,LFREE)
*
** 4 : Info about group of strings
*
*.First free address
      IDUMMY  = 60000
      KFREEL = 0
      CALL MEMGET('INTE',KFREEL,IDUMMY,WORK,KFREE,LFREE)
      if (NTEST.ge.1) then
        write(6,*) 'Getting first free address in strinf_gas'
        write(6,*) '----------------------------------------'
        write(6,*) 'first free address is   ',KFREEL
      end if
* will try to make dummy arrays
C     IDUMMY  = 2000
C     CALL MEMGET('INTE',KFREEL2,IDUMMY,WORK,KFREE,LFREE)
C     CALL MEMGET('INTE',KFREEL3,IDUMMY,WORK,KFREE,LFREE)
C     KFREEL2=KFREEL!+2
*
*. For spin up spinors :
      if (NTEST.ge.5) then
        write(6,*) '======================================='
        write(6,*)
        write(6,*) ' Info about UNBARRED spinor strings '
        write(6,*)
        write(6,*) '======================================='
      end if
*
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
        CALL WEIGHT_KRCC(WORK(KZ(IGRP)),IEL,NORB1,NORB2,NORB3,
     &                  MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                  WORK(KFREEL),IPRNT)
*. Number of strings per symmetry in a given group
        CALL NSTRSO_GAS_KRCC(IEL,NORB1,NORB2,NORB3,
     &                      MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                      WORK(KFREEL),NACOB,
     &                      WORK(KNSTSGP(1)),
     &                      WORK(KISTSGP(1)),
     &                      IOCTYPX,NSMST,IGRP,IPRNT)
*
*. Construct the strings ordered by symmetry
        if (IEL.gt.0) then
          CALL GENSTR_GAS_KRCC(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                        WORK(KISTSGP(1)),IGRP,
     &                        IOCTYPX,NSMST,WORK(KZ(IGRP)),
     &                        WORK(KFREEL),
     &                        WORK(KSTREO(IGRP)),WORK(KOCSTR(IGRP)),
     &                        WORK(KFREEL+IOCTYPX*NSMST),IGRP,IPRNT)
        end if
       CALL ICOPVE2(WORK(KNSTSGP(1)),1+(IGRP-1)*NSMST,NSMST,
     &              NSTFSMGP(1,IGRP))
       CALL ICOPVE2(WORK(KISTSGP(1)),1+(IGRP-1)*NSMST,NSMST,
     &              ISTFSMGP(1,IGRP))
      END DO
*
*. Min and max of sym for each group
*
      if (NTEST.ge.50) then
        write(6,*) 'NSTFSMGP array:'
        call iwrtma(NSTFSMGP,NSMST,NGRP,MXPNSMST,MXPSTT)
        call iwrtma(ISTFSMGP,NSMST,NGRP,MXPNSMST,MXPSTT)
      end if
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
*
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
C     CALL MEMCHK_KRCC(WORK)
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
        CALL WEIGHT_KRCC(WORK(KZ(IGRP)),IEL,NORB1,NORB2,NORB3,
     &                  MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                  WORK(KFREEL),IPRNT)
C       CALL MEMCHK_KRCC(WORK)
*. Number of strings per symmetry in a given group
        CALL NSTRSO_GAS_KRCC(IEL,NORB1,NORB2,NORB3,
     &                      MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                      WORK(KFREEL),NACOB,
     &                      WORK(KNSTSGP2(1)),
     &                      WORK(KISTSGP2(1)),
     &                      IOCTYPX,NSMST,IGRP,IPRNT)
*. Construct the strings ordered by symmetry
        if (IEL.gt.0) then
          CALL GENSTR_GAS_KRCC(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                        WORK(KISTSGP2(1)),IGRP,
     &                        IOCTYPX,NSMST,WORK(KZ(IGRP)),
     &                        WORK(KFREEL),
     &                        WORK(KSTREO2(IGRP)),WORK(KOCSTR2(IGRP)),
     &                        WORK(KFREEL+IOCTYPX*NSMST),IGRP,IPRNT)
        end if
        CALL ICOPVE2(WORK(KNSTSGP2(1)),1+(IGRP-1)*NSMST,NSMST,
     &               NSTFSMGP2(1,IGRP))
        CALL ICOPVE2(WORK(KISTSGP2(1)),1+(IGRP-1)*NSMST,NSMST,
     &               ISTFSMGP2(1,IGRP))
      END DO
*
*. Min and max of sym for each group
*
      if (NTEST.ge.50) then
        write(6,*) 'NSTFSMGP2 array:'
        call iwrtma(NSTFSMGP2,NSMST,NGRP,MXPNSMST,MXPSTT)
        call iwrtma(ISTFSMGP2,NSMST,NGRP,MXPNSMST,MXPSTT)
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
        CALL IWRTMA10(WORK(KNSTSGP(1)),NSMST,NGRP,NSMST,NGRP)
*
        IF(PNTGRP.GE.5) THEN
        WRITE(6,*)
     &  ' Number of strings per group and symmetry (spindown)'
        CALL IWRTMA10(WORK(KNSTSGP2(1)),NSMST,NGRP,NSMST,NGRP)
        END IF
*
      END IF
*
*
* 4.5 : Creation/Annihilation mappings between different
*       types of strings (Annihilation mapping not installed)
*
      DO IGRP = 1, NGRP
        print*,' IGRP ',IGRP,' ISTAC(IGRP,2) ',ISTAC(IGRP,2)
        IF(ISTAC(IGRP,2).NE.0) THEN
          IGAS = IGSFGP(IGRP)
          NGSOBP = NOBPT(IGAS)
*. First orbital in GAS spacce
          IGSOB = IELSUM(NOBPT,IGAS-1)+1
          IEL = NELFGP(IGRP)
          NSTINI = NSTFGP(IGRP)
          IZERO = 0
          CALL ISETVC(WORK(KSTSTM(IGRP,1)),IZERO,NGSOBP*NSTINI)
          CALL ISETVC(WORK(KSTSTM(IGRP,2)),IZERO,NGSOBP*NSTINI)
          JGRP = ISTAC(IGRP,2)
          CALL CRESTR_GAS_KRCC(WORK(KOCSTR(IGRP)),NSTFGP(IGRP),
     &         NSTFGP(JGRP),IEL,NGSOBP,IGSOB,WORK(KZ(JGRP)),
     &         WORK(KSTREO(JGRP)),0,IDUM,IDUM,
     &         WORK(KSTSTM(IGRP,1)),WORK(KSTSTM(IGRP,2)),NACOB,IPRNT)
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
          CALL ISETVC(WORK(KSTSTM2(IGRP,1)),IZERO,NGSOBP*NSTINI)
          CALL ISETVC(WORK(KSTSTM2(IGRP,2)),IZERO,NGSOBP*NSTINI)
          JGRP = ISTAC(IGRP,2)
          CALL CRESTR_GAS_KRCC(WORK(KOCSTR2(IGRP)),NSTFGP(IGRP),
     &         NSTFGP(JGRP),IEL,NGSOBP,IGSOB,WORK(KZ(JGRP)),
     &         WORK(KSTREO2(JGRP)),0,IDUM,IDUM,
     &         WORK(KSTSTM2(IGRP,1)),WORK(KSTSTM2(IGRP,2)),NACOB,IPRNT)
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
          CALL NSTPTP_GAS_KRCC(NGAS,ISPGPFTP(1,IGRPABS),
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
          CALL NSTPTP_GAS_KRCC(NGAS,ISPGPFTP(1,IGRPABS),
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
      CALL ZNELFSPGP_KRCC(IPRNT)
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
      if (NTEST.ge.50) then
        WRITE(6,*) ' Matrix NSTFSMSPGP '
        CALL IWRTMA(NSTFSMSPGP,NSMST,NTSPGP,MXPNSMST,NTSPGP)
      end if
*
*
* Possible occupation classes
*
      CALL OCCLS_KRCC(2,NMXOCCLS,WORK(KIOCLS),NACTEL,NGAS,
     &               IGSOCC(1,1),IGSOCC(1,2))
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
C?    WRITE(6,*) ' Memory Check at end of STRINF_GAS  '
*
      CALL QEXIT('STRIN')
      RETURN
      END
*
      SUBROUTINE ZSTINF_GAS_KRCC(IPRNT)
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
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Output
#include "stinf.inc"
*. Only the first element, i.e. ISTAC  is defined

*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
* ******************************************************************
      CALL ISETVC(ISTAC,0,2*MXPSTT)
      DO  IGAS = 1, NGAS
*. groups for a given gas spaces goes with increasing number of orbitals,
*  so the first space does not have any creation mapping
*  and the last space does not have any annihilation mapping
*
        MGRP = NGPSTR(IGAS)
        DO IGRP = 1, MGRP
          IIGRP = IGRP + IBGPSTR(IGAS) -1
          IF(IGRP.NE.1) THEN
*. Annihilation map is present : IIGRP => IIGRP - 1
            ISTAC(IIGRP,1) = IIGRP -1
          END IF
          IF(IGRP.NE.MGRP) THEN
*. Creation map is present : IIGRP => IIGRP + 1
             ISTAC(IIGRP,2) = IIGRP + 1
          END IF
        END DO
      END DO
*
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Type - type mapping array ISTAC '
        WRITE(6,*) ' =============================== '
        CALL IWRTMA(ISTAC,NGRP  ,2,MXPSTT,2)
      END IF
*
      RETURN
      END
*
      SUBROUTINE MEMSTR_GAS_KRCC(WORK,KFREE,LFREE)
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
#include "cgas.inc"
#include "gasstr.inc"
#include "stinf.inc"
*
      dimension WORK(*)
C
C#include "memint.h"
C
*
*. Start of string information
      IDUMMY = 0
      CALL MEMGET('INTE',KSTINF,IDUMMY,WORK,KFREE,LFREE)

      NTEST = 0
      IF(NTEST.NE.0)
     &WRITE(6,*) ' First word with string information',KSTINF
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        NSTRIN = NSTFGP(IGRP)
        LSTRIN = NSTRIN*NELFGP(IGRP)
        CALL MEMGET('INTE',KOCSTR(IGRP),LSTRIN,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KSTREO(IGRP),NSTRIN,WORK,KFREE,LFREE)
        IF(PNTGRP.GE.5) THEN
        CALL MEMGET('INTE',KOCSTR2(IGRP),LSTRIN,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KSTREO2(IGRP),NSTRIN,WORK,KFREE,LFREE)
        END IF
        CALL ISETVC(WORK(KOCSTR(IGRP)),0,LSTRIN)
        CALL ISETVC(WORK(KSTREO(IGRP)),0,NSTRIN)
        IF(PNTGRP.GE.5) THEN
        CALL ISETVC(WORK(KOCSTR2(IGRP)),0,LSTRIN)
        CALL ISETVC(WORK(KSTREO2(IGRP)),0,NSTRIN)
        END IF
*. For spin down strings
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
        CALL MEMGET('INTE',KNSTSGP(1),NSMST*NGRP,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KISTSGP(1),NSMST*NGRP,WORK,KFREE,LFREE)
      IF(PNTGRP.GE.5) THEN
        CALL MEMGET('INTE',KNSTSGP2(1),NSMST*NGRP,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KISTSGP2(1),NSMST*NGRP,WORK,KFREE,LFREE)
      END IF
        CALL ISETVC(WORK(KNSTSGP(1)),0,NSMST*NGRP)
        CALL ISETVC(WORK(KISTSGP(1)),0,NSMST*NGRP)
      IF(PNTGRP.GE.5) THEN
        CALL ISETVC(WORK(KNSTSGP2(1)),0,NSMST*NGRP)
        CALL ISETVC(WORK(KISTSGP2(1)),0,NSMST*NGRP)
      END IF

*
*. Number of strings per symmetry and offset for strings of given sym
*. for types
*
      DO ITP = 1, NSTTP
        CALL MEMGET('INTE',KNSTSO(ITP),NSPGPFTP(ITP)*NSMST,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KISTSO(ITP),NSPGPFTP(ITP)*NSMST,
     &              WORK,KFREE,LFREE)
      END DO
      IF(PNTGRP.GE.5) THEN
       DO ITP = 1, NSTTP
         LEN_T = NSPGPFTP(ITP)*NSMST
        CALL MEMGET('INTE',KNSTSO2(ITP),NSPGPFTP(ITP)*NSMST,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KISTSO2(ITP),NSPGPFTP(ITP)*NSMST,
     &              WORK,KFREE,LFREE)
       END DO
      END IF
*
**. Lexical adressing of arrays : use array indeces for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP
        CALL MEMGET('INTE',KZ(IGRP),NACOB*NELFGP(IGRP),WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KZ(IGRP)),0,NACOB*NELFGP(IGRP))
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
          CALL MEMGET('INTE',KSTSTM(IGRP,1),LENGTH,WORK,KFREE,LFREE)
          CALL MEMGET('INTE',KSTSTM(IGRP,2),LENGTH,WORK,KFREE,LFREE)
          IF(PNTGRP.GE.5) THEN
            CALL MEMGET('INTE',KSTSTM2(IGRP,1),LENGTH,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KSTSTM2(IGRP,2),LENGTH,WORK,KFREE,LFREE)
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
*
*. Symmetry of excitation connecting  strings of given symmetry
*
*
*. Occupation classes
*
      CALL MEMGET('INTE',KIOCLS,NMXOCCLS*NGAS,WORK,KFREE,LFREE)
*. Annihilation/Creation map of supergroup types
      CALL MEMGET('INTE',KSPGPAN,NTSPGP*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KSPGPCR,NTSPGP*NGAS,WORK,KFREE,LFREE)
*
*. Last word of string information
*
      KSTINE = 0
      CALL MEMGET('INTE',KSTINE,IDUMMY,WORK,KFREE,LFREE)
      IF(NTEST.NE.0)
     &WRITE(6,*) ' Last word with string information',KSTINE-1
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' R*8 words used for storing Strinformation ',
     &               KSTINE - KSTINF - 1
      END IF
*
      RETURN
      END
*
      SUBROUTINE WEIGHT_KRCC(Z,NEL,NORB1,NORB2,NORB3,
     &                      MNRS1,MXRS1,MNRS3,MXRS3,ISCR,NTEST)
*
* construct vertex weights
*
* Reverse lexical ordering is used for restricted space
*
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER Z(*),ISCR(*)
*
      NORB = NORB1 + NORB2 + NORB3
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' >>>> WEIGHT <<<<< '
        WRITE(6,*) ' NORB1 NORB2 NORB3 ',NORB1,NORB2,NORB3
        WRITE(6,*) ' NEL MNRS1 MXRS1 MNRS3 MXRS3 '
        WRITE(6,*)   NEL,MNRS1,MXRS1,MNRS3,MXRS3
      END IF
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
      CALL RSMXMN_KRCC(ISCR(KLMAX),ISCR(KLMIN),NORB1,NORB2,NORB3,
     &                NEL,MNRS1,MXRS1,MNRS3,MXRS3,NTEST)
*. Arc weights
      CALL GRAPW(ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NEL,NTEST)
*
      RETURN
      END
*
      SUBROUTINE RSMXMN_KRCC(MAXEL,MINEL,NORB1,NORB2,NORB3,NEL,
     &                      MIN1_T,MAX1,MIN3,MAX3_T,NTEST)
*
* Construct accumulated MAX and MIN arrays for a RAS set of strings
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION MINEL(*),MAXEL(*)
*
      NTEST = 00
*
      NORB = NORB1 + NORB2 + NORB3
*. accumulated max and min in each of the three spaces
*. ( required max and min at final orbital in each space )
COLD  MIN1A = MIN1_T
      MIN1A = MAX(MIN1_T,NEL-MAX3_T-NORB2)
      MAX1A = MAX1
*
      MIN2A = NEL - MAX3_T
      MAX2A = NEL - MIN3
*
      MIN3A = NEL
      MAX3A = NEL
*
      DO 100 IORB = 1, NORB
        IF(IORB .LE. NORB1 ) THEN
          MINEL(IORB) = MAX(MIN1A+IORB-NORB1,0)
          MAXEL(IORB) = MIN(IORB,MAX1A)
        ELSE IF ( NORB1.LT.IORB .AND. IORB.LE.(NORB1+NORB2)) THEN
          MINEL(IORB) = MAX(MIN2A+IORB-NORB1-NORB2,0)
          IF (NORB1.GT.0) MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1))
          MAXEL(IORB) = MIN(IORB,MAX2A)
        ELSE IF ( IORB .GT. NORB1 + NORB2 ) THEN
          MINEL(IORB) = MAX(MIN3A+IORB-NORB,0)
          IF (NORB1+NORB2.GT.0)
     &        MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1+NORB2))
          MAXEL(IORB) = MIN(IORB,MAX3A)
        END IF
  100 CONTINUE
*
      IF( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from RSMXMN '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
*
      SUBROUTINE NSTRSO_GAS_KRCC(NEL,NORB1,NORB2,NORB3,
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
         IF( NTEST .GE.500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
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
           IF( NTEST .GE.500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
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
             IF( NTEST .GE. 500) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
* Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
*. Symmetry of string
             ISYM = ISYMST_KRCC(IOC,NEL)
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
C
      RETURN
      END
*
      SUBROUTINE GENSTR_GAS_KRCC(NEL,NELMN1,NELMX1,NELMN3,NELMX3,
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
*
      CALL ISETVC(LSTASO,0,NOCTYP*NSMST)
*
      NTEST0 = 000
      NTEST = MAX(NTEST0,IPRNT)
*
      IF( NTEST .GE. 10) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' GENSTR speaking '
        WRITE(6,*)  ' =============== '
      END IF
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
         IF( NTEST .GE. 500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
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
           IF( NTEST .GE. 500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
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
             IF( NTEST .GE. 500 ) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
* Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
*. Symmetry
*                   ISYMST_REL(STRING,NEL)
             ISYM = ISYMST_KRCC(IOC,NEL)
             IF(NTEST.GE.500)  WRITE(6,*) ' ISYM= ',ISYM
             IF(NTEST.GE.500)  WRITE(6,*) ' NEL= ',NEL
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
               IF(NTEST.GT.500) WRITE(6,*) ' LEXCI,LACTU',
     &         LEXCI,LACTU
               IF(NTEST.GT.500) WRITE(6,*) ' NEL,IOC ',NEL
               IF(NTEST.GT.500) CALL IWRTMA(IOC,1,NEL,1,NEL)
               CALL ICOPVE(IOC,STRING(1,LACTU),NEL)
               IF(NTEST.GT.500) CALL IWRTMA(STRING,1,NEL,1,LACTU+1)
             END IF
*
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)
     &     GOTO 901
 1003 CONTINUE
 1001 CONTINUE
*
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
*
      RETURN
      END
*
      FUNCTION ISYMST_KRCC(STRING,NEL)
*
* Master routine for symmetry of string
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input ( PNTGRP is used )
C
#include "mxpdim.inc"
#include "lucinp.inc"
*. Specific input
      INTEGER STRING(*)
      IF(PNTGRP.EQ.1) THEN
*.D2h
        ISYMST_KRCC = ISYMS1_KRCC(STRING,NEL)
      ELSE IF(PNTGRP.GE.2.AND.PNTGRP.LE.4) THEN
*.Cinfv Dinfh O3
        ISYMST_KRCC = ISYMS2_KRCC(STRING,NEL)
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
        ISYMST_KRCC = ISYMSDG_KRCC(STRING,NEL)
      ELSE
        WRITE(6,*) ' Sorry PNTGRP option not programmed ', PNTGRP
        WRITE(6,*) ' Enforced stop in ISYMST '
        Call Abend1( 5 )
      END IF
*
      RETURN
      END
*
      FUNCTION ISYMS1_KRCC(STRING,NEL)
*
* Symmmetry of string, D2H version
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
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
      ISYMS1_KRCC = ISYM
*
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ISYMS1, String and symmetry '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,*) ISYM
      END IF
*
      RETURN
      END
*
      FUNCTION ISYMS2_KRCC(STRING,NEL)
*
* Symmetry of string STRING, D inf h, C inf v, O3 version
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
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
      ISYMS2_KRCC = ISYM
*
      NTEST = 0
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' STRING '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,'(A,3I3)') ' MLSTR, IPARI ISYMS2 ', MLSTR,IPARI,ISYM
      END IF
*
      RETURN
      END
*
      FUNCTION ISYMSDG_KRCC(STRING,NEL)
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
      ISYMSDG_KRCC = ISYM
*
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ISYMS1, String and symmetry '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,*) ISYM
      END IF
*
      RETURN
      END
*
      SUBROUTINE CRESTR_GAS_KRCC(STRING,NSTINI,NSTINO,NEL,NORB,IORBOF,
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
      NTEST0 =  00
      NTEST = MAX(IPRNT,NTEST0)
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*)  ' ==================== '
        WRITE(6,*)  ' CRESTR_KRCC speaking '
        WRITE(6,*)  ' ==================== '
        WRITE(6,*)
         WRITE(6,*) ' Number of input electrons ', NEL
      END IF
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
*
      RETURN
      END
*
      SUBROUTINE NSTPTP_GAS_KRCC(NGAS,ISPGRP,NSTSGP,NSMST,
     &                          NSTSSPGP,IGRP,MXNSTR)
*
* From number of strings per group and sym to number of strings
* per supergroup and sym for given super group.
*
* Jeppe Olsen , Fall of 94
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ISPGRP(NGAS),NSTSGP(NSMST,*)
*. Output
      DIMENSION NSTSSPGP(NSMST,IGRP)
*. Scratch
#include "mxpdim.inc"
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' NSTPTP_GAS_KRCC is speaking '
        WRITE(6,*) ' ============================'
*
        WRITE(6,*) ' Supergroup in action '
        CALL IWRTMA(ISPGRP,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' NSMST = ' , NSMST
        write(6,*) 'First 3 columns of NSTSGP'
        CALL IWRTMA(NSTSGP,NSMST,3,NSMST,3)
      END IF
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
      NCOUNT = 0
 1000 CONTINUE
      IF(IFIRST.EQ.1) THEN
        CALL ISETVC(ISM,1,NGAS)
        IFIRST = 0
        NONEW = 0
        NCOUNT=NCOUNT + 1
        IF(NCOUNT.GE.2) STOP'NCOUNT ERROR'
      ELSE
C       NXTNUM3(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
        CALL NXTNUM3(ISM,NGAS,MNSM,MXSM,NONEW)
      END IF
      IF(NONEW.EQ.0) THEN
*. Symmetry of current combination and number of strings in this supergroup
        ISMSPGP = ISM(1)
        NST = NSTSGP(ISM(1),ISPGRP(1))
        DO JGRP = 2, NGAS
          CALL SYMCOM_KRCC(3,7,ISMSPGP,ISM(JGRP),ISMSPGPO)
          ISMSPGP = ISMSPGPO
          NST = NST * NSTSGP(ISM(JGRP),ISPGRP(JGRP))
        END DO
C       WRITE(6,*) ' Symmetry of groups '
C       CALL IWRTMA(ISM,1,NGAS,1,NGAS)
C       WRITE(6,*) ' Symmetry of super group and number of strings '
C       WRITE(6,*) ISMSPGP , NST
        NSTSSPGP(ISMSPGP,IGRP) = NSTSSPGP(ISMSPGP,IGRP) + NST
        GOTO 1000
      END IF
*
      MXNSTR = -1
      DO ISTRSM = 1, NSMST
        MXNSTR = MAX(MXNSTR,NSTSSPGP(ISTRSM,IGRP))
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
     &  ' Number of strings per symmetry for supergroup',IGRP
        CALL IWRTMA10(NSTSSPGP(1,IGRP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Largest number of strings of given sym ',MXNSTR
      END IF
*
      RETURN
      END
*
      SUBROUTINE ZNELFSPGP_KRCC(NTESTG)
*
* Generate for each supergroup the number of electrons in each active
* orbital space and store in NELFSPGP
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. input
#include "mxpdim.inc"
#include "cgas.inc"
*. Input and Output ( NELFSPGP(MXPNGAS,MXPSTT) )
#include "gasstr.inc"
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
*
      DO ITP = 1, NSTTP
        NSPGP = NSPGPFTP(ITP)
        IBSPGP = IBSPGPFTP(ITP)
        DO ISPGP = IBSPGP,IBSPGP + NSPGP - 1
          DO IGAS = 1, NGAS
            NELFSPGP(IGAS,ISPGP) = NELFGP(ISPGPFTP(IGAS,ISPGP))
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Distribution of electrons in Active spaces '
        DO ITP = 1, NSTTP
          WRITE(6,*) ' String type ', ITP
          WRITE(6,*) ' Row : active space, Column: supergroup '
          NSPGP = NSPGPFTP(ITP)
          IBSPGP = IBSPGPFTP(ITP)
          CALL IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE SET_HOP_DBG_KRCC(IHTYPE,JCMBSPC,WORK,LWORK)
*
* Construct orbital excitation types for different choices
* of Hamiltonian operator, double group version
*
* This routine was taken over by Lasse since we have different
* number of H-type operators for the CI-driven CC and the
* commutator driven CC. This is only for the commutator driven version.
*
* IHTYPE = 1 : Include only Delta MK = 0 terms
* IHTYPE = 2 : Include Delta Mk = 0 for two-electron part
*              and complete one-electron part
* IHTYPE = 3 : Include all terms for general relativistic Hamiltonian
*
#include "implicit.inc"
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
*. Output
#include "ctcc.inc"
#include "ctccp.inc"
C
      DIMENSION WORK(LWORK)
C
C
#include "memint.h"
*
*. Number of spin-orbital exctation operators in H
      CALL MEMGET('INTE',KLSOBEX_TEMP,4*NGAS,WORK,KFREE,LFREE)
      IFLAG = 1
      CALL GET_HX_RELA_KRCC(IFLAG,IHTYPE,NOBPT,NHX,
     &                 IFTONE,IFTTWO,ILTTWO,WORK(KLSOBEX_TEMP),
     &                 DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                 JCMBSPC,IPRNT)
      NSPOBEX_TP = NHX
*
*. and the actual spin-orbital excitations
C
C Needs to consider names of arrayes since these have all been allocated before!!!!
C so just comment out and overwrite to particle indexing (from integral idexing)
C Will be allocated here again for KRCC since there is a new structure
C All arrays will be keept since they will be needed
      CALL MEMGET('INTE',KLSOBEX,4*NGAS*NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIHIND,4*NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNNHX,NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLABEXTP,4*NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLOP_REO,4*NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIOPREO,NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINC_HX,NHX,WORK,KFREE,LFREE)
      IZERO = 0
      CALL ISETVC(WORK(KINC_HX),IZERO,NHX)
      IFLAG = 0
      CALL GET_HX_RELA_KRCC(IFLAG,IHTYPE,NOBPT,NHX,
     &                 IFTONE,IFTTWO,ILTTWO,WORK(KLSOBEX_TEMP),
     &                 WORK(KLSOBEX),WORK(KLIHIND),WORK(KSIGNNHX),
     &                 WORK(KLABEXTP),WORK(KLOP_REO),
     &                 WORK(KSIOPREO),JCMBSPC,IPRNT)
C     CALL STUPID_ROUTINE(WORK(KLSOBEX),NSPOBEX_TP,NGAS)
*
*. Space for info in the different parts of Hamiltonian
* KLLSOBEX : Length of this block
* KLIBSOBEX : Offset for this block
*
      CALL MEMGET('INTE',KLLSOBEX,NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIBSOBEX,NSPOBEX_TP,WORK,KFREE,LFREE)
      ITOP_SM = 1
      CALL IDIM_TCC_KRCC(WORK(KLSOBEX),NSPOBEX_TP,ITOP_SM,
     &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &     MX_SBSTR,
     &     IFTONE,IFTTWO,-1,N1ELINT,N2ELINT)
*
      MX_ST_TSOSO_BLK_MX = MX_ST_TSOSO_BLK
      LEN_T_VEC_MX = LEN_T_VEC
*
      RETURN
      END
*
      SUBROUTINE SET_HOP_DBG_FOR_NEWCCV(IHTYPE,JCMBSPC,WORK,KFREE,LFREE)
*
* Construct orbital excitation types for different choices
* of Hamiltonian operator, double group version
*
*
* IHTYPE = 1 : Include only Delta MK = 0 terms
* IHTYPE = 2 : Include Delta Mk = 0 for two-electron part
*              and complete one-electron part
* IHTYPE = 3 : Include all terms for general relativistic Hamiltonian
*
#include "implicit.inc"
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
*. Output
#include "ctcc.inc"
#include "ctccp.inc"
*
      DIMENSION WORK(*)
C
C
C#include "memint.h"
*
*. Number of spin-orbital exctation operators in H
      CALL MEMGET('INTE',KLSOBEX_TEMP,4*NGAS,WORK,KFREE,LFREE)
      IFLAG = 1
      CALL GET_HX_RELA_FOR_NEWCCV(IFLAG,IHTYPE,NOBPT,NHX,
     &                 IFTONE,IFTTWO,ILTTWO,WORK(KLSOBEX_TEMP),
     &                 DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                 DUMMY,JCMBSPC,IPRNT)
      NSPOBEX_TP = NHX
*
*. and the actual spin-orbital excitations
      CALL MEMGET('INTE',KLSOBEX,4*NGAS*NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIHIND,4*NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIHINDOP,4*NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNNHX,NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLABEXTP,4*NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLOP_REO,4*NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIOPREO,NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINC_HX,NHX,WORK,KFREE,LFREE)
      IZERO = 0
      CALL ISETVC(WORK(KINC_HX),IZERO,NHX)
      IFLAG = 0
* Notice KLIHINDOP becomes IHINDX and KLIHIND IHINDXOP
* It is KLIHIND that is used later in the code!!!
      CALL GET_HX_RELA_FOR_NEWCCV(IFLAG,IHTYPE,NOBPT,NHX,
     &                 IFTONE,IFTTWO,ILTTWO,WORK(KLSOBEX_TEMP),
     &                 WORK(KLSOBEX),WORK(KLIHINDOP),WORK(KLIHIND),
     &                 WORK(KSIGNNHX),
     &                 WORK(KLABEXTP),WORK(KLOP_REO),
     &                 WORK(KSIOPREO),WORK(KINC_HX),JCMBSPC,IPRNT)
*
*. Space for info in the different parts of Hamiltonian
* KLLSOBEX : Length of this block
* KLIBSOBEX : Offset for this block
*
      CALL MEMGET('INTE',KLLSOBEX,NSPOBEX_TP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIBSOBEX,NSPOBEX_TP,WORK,KFREE,LFREE)
*. Dimensions for total symmetric operators, sym = 1
      ITOTSYM = 1
      MSCOMB_T = 0
      CALL IDIM_TCC_KRCC(WORK(KLSOBEX),NSPOBEX_TP,ITOTSYM,
     &           MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &           WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &           MX_TBLK_AS,IFTONE,IFTTWO,ILTTWO,
     &           N1ELINT,N2ELINT)
*
      MX_ST_TSOSO_BLK_MX = MX_ST_TSOSO_BLK
      LEN_T_VEC_MX = LEN_T_VEC
*
      RETURN
      END
*
      SUBROUTINE GET_HX_RELA_KRCC(IFLAG,IHTYPE,NOBPT,NHX,
     &                       IFTONE,IFTTWO,ILTTWO,IHX_LOC,IHX,IHINDX,
     &                       SIGN_NHX,LABEL,IOP_REO,SIGN_OPREO,
     &                       INCLUDED_HX,JCMBSPC,IPRNT)
*
*
*
* Obtain number (IFLAG = 1)
* or number and types (IFLAG.NE.1) of excitations in
* relativistic Hamiltonian
*
*--------------------------------
* Jeppe Olsen and Timo Fleig
*   March/April 2001
* Orbital and electron number limitations
*   Timo Fleig, June 2001
* Updated to exclude unwanted components, Oct 2001
*   JO and TF
* Included LABEL specification of excitation type for DIRAC
*   (and other interfaces)
* Changed to accomedate for several calc
* Lasse Nov 04
*
* Finally changed to operator form in 08 by Lasse with operator index
*
* mapping from integral to operator done as following
* sign     1     -1
*        1->1   1->1
*        3->2   3->2
*        4->3   4->4
*        2->4   2->3
*
*--------------------------------
*
*
* IHTYPE = 1 : Include only Delta MK = 0 terms
* IHTYPE = 2 : Include Delta Mk = 0 for two-electron part
*              and complete one-electron part
* IHTYPE = 3 : Include all terms for relativistic Hamiltonian
*              in a basis of non-relativistic orbitals
* IHTYPE = 4 : Include all terms for relativistic Hamiltonian
*              in a basis of relativistic spinors
*
* LABEL : 1  for unbarred (alpha)
*         -1 for barred   (beta)
*         and refers to integral (not operator!) indexing !!
*
* IOP_REO : Reordering array for output density elements.
*           Relate definition of operator in this routine with
*           general output CA AA CB AB (or corresponding scheme
*           in DMK .ne. 0 cases). Not used in 1-el. case, but
*           due to NHX type counting need to be defined.
*
* SIGN_OPREO : Sign for this operator reordering
*          1.0 or -1.0 : Non-redundant density and sign
*          0.0         : Redundant density type
*
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. input
      dimension NOBPT(NGAS)
*. output
      INTEGER IHX(NGAS,4,*),IHINDX(4,*),LABEL(4,*),IOP_REO(4,*)
      dimension SIGN_NHX(*),SIGN_OPREO(*)
*. Local scratch
      INTEGER IHX_LOC(4*NGAS)
*. An excitation is stored as :
*  Alpha creation,beta creation,alpha annihilation,beta annihilation.
*
      NTESTL = 100
      NTEST = max(NTESTL,IPRNT)
      NTEST = 00
*
      IZERO = 0
      ZERO = 0.D0
      ONE = 1.0D0
*
      NHX = 0
*

* Determine maximum number of electrons per GAS first:
* FIXME  : Hardwired only one type of GASSPC calculation for the
*          moment. This would require other changes as well.
*          Making those at the moment - Lasse
      MXNELGS(1) = IGSOCCX(1,2,JCMBSPC)
      do IGAS=2,NGAS,1
        MXNELGS(IGAS)= IGSOCCX(IGAS,2,JCMBSPC)-IGSOCCX(IGAS-1,1,JCMBSPC)
        MXNELGS(IGAS)= min(MXNELGS(IGAS),2*NOBPT(IGAS))
      end do
*
*. Delta MK = 0 terms
*
* operator:
* a+i alpha a j alpha
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          if (NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(1-1)*NGAS) = 1
              IHX_LOC(JGAS+(3-1)*NGAS) = 1
*. Is this excitation type included ?
              CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(IHINDX(1,NHX),NHX) = 1
                  LABEL(IHINDX(2,NHX),NHX) = 1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            end if
          end if
        END DO
      END DO
*
* operator:
* a+i beta a j beta
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          if (NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(2-1)*NGAS) = 1
              IHX_LOC(JGAS+(4-1)*NGAS) = 1
              CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(IHINDX(1,NHX),NHX) = -1
                  LABEL(IHINDX(2,NHX),NHX) = -1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            end if
          end if
        END DO
      END DO
*
* For general purposes: First type with
*                       1-electron and 2-electron integrals:
C     print*,'NHX one electron',NHX
      IFTONE = 1
      IFTTWO = NHX + 1
*
* operator:
*. a+i alpha  a+j alpha  a k alpha  a l alpha, ordered so i.ge.j, k.ge.l
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1))) then
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS) + 1
                IHX_LOC(KGAS+(3-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS) + 1
C               IHX(IGAS,1,NHX) = 1
C               IHX(JGAS,1,NHX) = IHX(JGAS,1,NHX)+1
* means IGAS=JGAS, so the number of creators of this type is 2!
C               IHX(KGAS,3,NHX) = 1
C               IHX(LGAS,3,NHX) = IHX(LGAS,3,NHX)+1
                CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN
                  WRITE(6,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  a | +a  a)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 2
                    IHINDX(3,NHX) = 3
                    IHINDX(4,NHX) = 4
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = 1
                    LABEL(IHINDX(3,NHX),NHX) = 1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 2
                    SIGN_OPREO(NHX) = -ONE
                  end if
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
* operator:
*. a+i beta a+j beta a k beta a l beta, ordered so i.ge.j, k.ge.l
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1))) then
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(2-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS) + 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS) + 1
                CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN
                  WRITE(6,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  b | +b  b)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 2
                    IHINDX(3,NHX) = 3
                    IHINDX(4,NHX) = 4
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = -1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = -1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 2
                    SIGN_OPREO(NHX) = -ONE
                  end if
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
* operator:
*. a+i alpha  a+j beta  a k beta  a l alpha
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, NGAS
              if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &            MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &            NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &            NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1) then
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = 1
                CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN
                  WRITE(6,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  a | +b  b)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 2
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 3
                    SIGN_NHX(NHX) = -1.D0 * ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 2
                    IOP_REO(4,NHX) = 3
                    SIGN_OPREO(NHX) = -1.D0 * ONE
                  end if
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
C     print*,'NHX before ub|bu ',NHX
      if (IHTYPE.ge.40) then
*
* operator:
*. a+i alpha  a+j beta  a k alpha  a l beta
*   Generates type (ub|bu) of integrals.
*
* Changes made by Lasse !!!
        do IGAS = 1,NGAS,1
          do JGAS = 1,NGAS,1
            do KGAS = 1,NGAS,1
              do LGAS = 1,NGAS,1
                if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &              MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1) then
                  call isetvc(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = 1
                  call caab_included(IHX_LOC,INCLUDED)
                  if (INCLUDED.eq.0) then
                    write(6,*) 'Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  else
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      call icopve(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  b | +b  a)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 2
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 4
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 3
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  end if
                end if
              end do
            end do
          end do
        end do
*
      end if
*
* Last type of non-spin-flip integrals
      ILTTWO = NHX
C     print*,'NHX after non-spin-flip',NHX
*
* one-electron spin orbit
*
      IF (IHTYPE.GT.20) THEN
*
* operator:
* a+i alpha  j beta
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &          NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(1-1)*NGAS) = 1
              IHX_LOC(JGAS+(4-1)*NGAS) = 1
              CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(1,NHX) = 1
                  LABEL(2,NHX) = -1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            END IF
          END DO
        END DO
*
* operator:
* a+i beta a j alpha
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &          NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(2-1)*NGAS) = 1
              IHX_LOC(JGAS+(3-1)*NGAS) = 1
              CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  SIGN_NHX(NHX) = ONE
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  LABEL(1,NHX) = -1
                  LABEL(2,NHX) = 1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            END IF
          END DO
*    quaternionic matrix groups only!
        END DO
      END IF
C     print*,'NHX after one-electron spinflip',NHX
*
* Two-electron terms with Delta mk .ne. 0
*
      if (IHTYPE.eq.3.or.IHTYPE.eq.5) THEN
*
**********************************************************
*  Spin-dependent spinor basis. E.g. DIRAC environment   *
*    (ub|uu), (bu|uu), (bu|bb), (ub|bb) integral classes *
**********************************************************
*
* operator:
*. a+i alpha a+j alpha a k alpha a l beta, i .ge. j
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS
                if ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &               MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &               NOBPT(IGAS).ge.2.and.NOBPT(KGAS).ge.1.and.
     &               NOBPT(LGAS).ge.1).or.
     &               (IGAS.ne.JGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS)+1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = 1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  a | b  a+)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 2
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 4
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = -1.D0 * ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
* operator:
*. a+i beta a+j alpha a k alpha a l alpha, k .ge. l
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if ((KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &               MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &               NOBPT(KGAS).ge.2.and.NOBPT(IGAS).ge.1.and.
     &               NOBPT(JGAS).ge.1).or.
     &               (KGAS.ne.LGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS)+1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  a | a  a+)
* so reorder accordingly:
                      IHINDX(1,NHX) = 2
                      IHINDX(2,NHX) = 1
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 4
Cori                  SIGN_NHX(NHX) = ONE
                      SIGN_NHX(NHX) = -1.D0 * ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = 1
                      IOP_REO(1,NHX) = 3
                      IOP_REO(2,NHX) = 1
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
* operator:
*. a+i beta  a+j beta  a k beta  a l alpha , i .ge. j
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS
                if ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &               MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &               NOBPT(IGAS).ge.2.and.NOBPT(KGAS).ge.1.and.
     &               NOBPT(LGAS).ge.1).or.
     &               (IGAS.ne.JGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(2-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS)+1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = 1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  a | b  b+)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 2
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 3
                      SIGN_NHX(NHX) = -1.0D0 * ONE
Cori                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = 1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
* operator:
*. a+i alpha  a+j beta  a k beta  a l beta , k .ge. l
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if ((KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &               MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &               NOBPT(KGAS).ge.2.and.NOBPT(IGAS).ge.1.and.
     &               NOBPT(JGAS).ge.1).or.
     &               (KGAS.ne.LGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS)+1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  b | +b  b)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 2
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 4
C No the sign is not so confusing
                      SIGN_NHX(NHX) = ONE  !!??
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
      end if
*
**********************************************************
*  Spin-dependent spinor basis. E.g. DIRAC environment   *
*      (bu|bu), (ub|ub) integral classes                 *
**********************************************************
*
C     print*,'NHX after IHTYPE 3 or 5',NHX
      if (IHTYPE.eq.4) then
*
* operator:
*. a+i alpha  a+j alpha  a k beta  a l beta, i.ge.j, k.ge.l
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1))) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS)+1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS)+1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  b | +a  b )
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 2
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 4
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 3
                      IOP_REO(4,NHX) = 2
                      SIGN_OPREO(NHX) = -ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
* operator:
*. a+i beta  a+j beta  a k alpha  a l alpha, i.ge.j, k.ge.l
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1))) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(2-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS)+1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS)+1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  a | +b  a)
* so reorder accordingly:
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 2
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 4
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = 1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 3
                      IOP_REO(4,NHX) = 2
                      SIGN_OPREO(NHX) = -ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
      end if
C     print*,'NHX at the end',NHX
*
      if (IFLAG.ne.1.and.NTEST.ge.500) then
        do I=1,NHX
          write(6,*) 'I,SIGN_NHX ',I,SIGN_NHX(I)
        end do
      end if
*
      NTEST = 00
      if (IFLAG.ne.1) then
        if (NTEST.ge.100) then
          WRITE(6,*) ' Information from  GET_HX_RELA_CC '
          WRITE(6,*)
          write(6,*) ' Max. no of electrons per GAS:'
          call iwrtma(MXNELGS,1,NGAS,1,NGAS)
          WRITE(6,*) ' The excitation types '
          CALL WRT_SPOX_TP_CC_KRCC(IHX,NHX)
        end if
        write(6,*)
        WRITE(6,*) ' Number of excitation types ', NHX
      end if
        WRITE(6,*) ' Number of excitation types ', NHX
*
      RETURN
      END
*
      SUBROUTINE GET_HX_RELA_FOR_NEWCCV(IFLAG,IHTYPE,NOBPT,NHX,
     &                       IFTONE,IFTTWO,ILTTWO,IHX_LOC,IHX,
     &                       IHINDX,IHINDXOP,
     &                       SIGN_NHX,LABEL,IOP_REO,SIGN_OPREO,
     &                       INCLUDED_HX,JCMBSPC,IPRNT)
*
*
*
* Obtain number (IFLAG = 1)
* or number and types (IFLAG.NE.1) of excitations in
* relativistic Hamiltonian
*
*--------------------------------
* Jeppe Olsen and Timo Fleig
*   March/April 2001
* Orbital and electron number limitations
*   Timo Fleig, June 2001
* Updated to exclude unwanted components, Oct 2001
*   JO and TF
* Included LABEL specification of excitation type for DIRAC
*   (and other interfaces)
* Changed to accomedate for several calc
* Lasse Nov 04
*--------------------------------
*
* IHTYPE = 1 : Include only Delta MK = 0 terms
* IHTYPE = 2 : Include Delta Mk = 0 for two-electron part
*              and complete one-electron part
* IHTYPE = 3 : Include all terms for relativistic Hamiltonian
*              in a basis of non-relativistic orbitals
* IHTYPE = 4 : Include all terms for relativistic Hamiltonian
*              in a basis of relativistic spinors
*
* LABEL : 1  for unbarred (alpha)
*         -1 for barred   (beta)
*         and refers to integral (not operator!) indexing !!
*
* IOP_REO : Reordering array for output density elements.
*           Relate definition of operator in this routine with
*           general output CA AA CB AB (or corresponding scheme
*           in DMK .ne. 0 cases). Not used in 1-el. case, but
*           due to NHX type counting need to be defined.
*
* SIGN_OPREO : Sign for this operator reordering
*          1.0 or -1.0 : Non-redundant density and sign
*          0.0         : Redundant density type
*
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. input
      dimension NOBPT(NGAS)
*. output
      INTEGER IHX(NGAS,4,*),IHINDX(4,*),IHINDXOP(4,*)
      INTEGER LABEL(4,*),IOP_REO(4,*)
      INTEGER INCLUDED_HX(*)
      dimension SIGN_NHX(*),SIGN_OPREO(*)
*. Local scratch
      INTEGER IHX_LOC(4*NGAS)
*. An excitation is stored as :
*  Alpha creation,beta creation,alpha annihilation,beta annihilation.
*
      NTESTL = 100
      NTEST = max(NTESTL,IPRNT)
*
      IZERO = 0
      ZERO = 0.D0
      ONE = 1.0D0
*
      NHX = 0
*

* Determine maximum number of electrons per GAS first:
* FIXME  : Hardwired only one type of GASSPC calculation for the
*          moment. This would require other changes as well.
*          Making those at the moment - Lasse
      MXNELGS(1) = IGSOCCX(1,2,JCMBSPC)
      do IGAS=2,NGAS,1
        MXNELGS(IGAS)= IGSOCCX(IGAS,2,JCMBSPC)-IGSOCCX(IGAS-1,1,JCMBSPC)
        MXNELGS(IGAS)= min(MXNELGS(IGAS),2*NOBPT(IGAS))
      end do
*
*. Delta MK = 0 terms
*
* operator:
* a+i alpha a j alpha
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          if (NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(1-1)*NGAS) = 1
              IHX_LOC(JGAS+(3-1)*NGAS) = 1
*. Is this excitation type included ?
              CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  IHINDXOP(1,NHX) = 1
                  IHINDXOP(2,NHX) = 2
                  IHINDXOP(3,NHX) = 0
                  IHINDXOP(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(IHINDX(1,NHX),NHX) = 1
                  LABEL(IHINDX(2,NHX),NHX) = 1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            end if
          end if
        END DO
      END DO
*
* operator:
* a+i beta a j beta
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          if (NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(2-1)*NGAS) = 1
              IHX_LOC(JGAS+(4-1)*NGAS) = 1
              CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  IHINDXOP(1,NHX) = 1
                  IHINDXOP(2,NHX) = 2
                  IHINDXOP(3,NHX) = 0
                  IHINDXOP(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(IHINDX(1,NHX),NHX) = -1
                  LABEL(IHINDX(2,NHX),NHX) = -1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            end if
          end if
        END DO
      END DO
*
* For general purposes: First type with
*                       1-electron and 2-electron integrals:
      IFTONE = 1
      IFTTWO = NHX + 1
*
* operator:
*. a+i alpha  a+j alpha  a k alpha  a l alpha, ordered so i.ge.j, k.ge.l
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1))) then
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS) + 1
                IHX_LOC(KGAS+(3-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS) + 1
C               IHX(IGAS,1,NHX) = 1
C               IHX(JGAS,1,NHX) = IHX(JGAS,1,NHX)+1
* means IGAS=JGAS, so the number of creators of this type is 2!
C               IHX(KGAS,3,NHX) = 1
C               IHX(LGAS,3,NHX) = IHX(LGAS,3,NHX)+1
                CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN
                  WRITE(6,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  a | +a  a)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    IHINDXOP(1,NHX) = 1
                    IHINDXOP(2,NHX) = 2
                    IHINDXOP(3,NHX) = 3
                    IHINDXOP(4,NHX) = 4
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = 1
                    LABEL(IHINDX(3,NHX),NHX) = 1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 2
                    SIGN_OPREO(NHX) = -ONE
                  end if
                END IF
              ELSE
C These Hamiltonian classes will never contribute to anything but we will include
C them because it makes it easier to fetch the integrals. These are marked in
C INCLUDED_HX and then sorted out later without being calculated. So no extra effort is needed.
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS) + 1
                IHX_LOC(KGAS+(3-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS) + 1
C               IHX(IGAS,1,NHX) = 1
C               IHX(JGAS,1,NHX) = IHX(JGAS,1,NHX)+1
* means IGAS=JGAS, so the number of creators of this type is 2!
C               IHX(KGAS,3,NHX) = 1
C               IHX(LGAS,3,NHX) = IHX(LGAS,3,NHX)+1
                CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                    INCLUDED_HX(NHX) = 1
* This excitation corresponds to the integral
*  (+a  a | +a  a)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    IHINDXOP(1,NHX) = 1
                    IHINDXOP(2,NHX) = 2
                    IHINDXOP(3,NHX) = 3
                    IHINDXOP(4,NHX) = 4
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = 1
                    LABEL(IHINDX(3,NHX),NHX) = 1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 2
                    SIGN_OPREO(NHX) = -ONE
                  end if
              END IF
            END DO
          END DO
        END DO
      END DO
*
* operator:
*. a+i beta a+j beta a k beta a l beta, ordered so i.ge.j, k.ge.l
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1))) then
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(2-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS) + 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS) + 1
                CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN
                  WRITE(6,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  b | +b  b)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    IHINDXOP(1,NHX) = 1
                    IHINDXOP(2,NHX) = 2
                    IHINDXOP(3,NHX) = 3
                    IHINDXOP(4,NHX) = 4
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = -1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = -1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 2
                    SIGN_OPREO(NHX) = -ONE
                  end if
                END IF
              ELSE
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(2-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS) + 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS) + 1
                CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                    INCLUDED_HX(NHX) = 1
* This excitation corresponds to the integral
*  (+b  b | +b  b)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    IHINDXOP(1,NHX) = 1
                    IHINDXOP(2,NHX) = 2
                    IHINDXOP(3,NHX) = 3
                    IHINDXOP(4,NHX) = 4
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = -1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = -1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 2
                    SIGN_OPREO(NHX) = -ONE
                  end if
              END IF
            END DO
          END DO
        END DO
      END DO
*
* operator:
*. a+i alpha  a+j beta  a k beta  a l alpha
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, NGAS
              if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &            MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &            NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &            NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1) then
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = 1
                CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN
                  WRITE(6,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  a | +b  b)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    IHINDXOP(1,NHX) = 1
                    IHINDXOP(2,NHX) = 2
                    IHINDXOP(3,NHX) = 4
                    IHINDXOP(4,NHX) = 3
                    SIGN_NHX(NHX) = -1.D0 * ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 2
                    IOP_REO(4,NHX) = 3
                    SIGN_OPREO(NHX) = -1.D0 * ONE
                  end if
                END IF
              ELSE
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = 1
                CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                    INCLUDED_HX(NHX) = 1
* This excitation corresponds to the integral
*  (+a  a | +b  b)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    IHINDXOP(1,NHX) = 1
                    IHINDXOP(2,NHX) = 2
                    IHINDXOP(3,NHX) = 4
                    IHINDXOP(4,NHX) = 3
                    SIGN_NHX(NHX) = -1.D0 * ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 2
                    IOP_REO(4,NHX) = 3
                    SIGN_OPREO(NHX) = -1.D0 * ONE
                  end if
              END IF
            END DO
          END DO
        END DO
      END DO
*
      if (IHTYPE.ge.40) then
*
* operator:
*. a+i alpha  a+j beta  a k alpha  a l beta
*   Generates type (ub|bu) of integrals.
*
* Changes made by Lasse !!!
        do IGAS = 1,NGAS,1
          do JGAS = 1,NGAS,1
            do KGAS = 1,NGAS,1
              do LGAS = 1,NGAS,1
                if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &              MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1) then
                  call isetvc(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = 1
                  call caab_included(IHX_LOC,INCLUDED)
                  if (INCLUDED.eq.0) then
                    write(6,*) 'Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  else
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      call icopve(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  b | +b  a)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 3
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  end if
                ELSE
                  call isetvc(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = 1
                  call caab_included(IHX_LOC,INCLUDED)
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      call icopve(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                      INCLUDED_HX(NHX) = 1
* This excitation corresponds to the integral
*  (+a  b | +b  a)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 3
                      SIGN_OPREO(NHX) = ZERO
                    end if
                end if
              end do
            end do
          end do
        end do
*
      end if
*
* Last type of non-spin-flip integrals
      ILTTWO = NHX
*
* one-electron spin orbit
*
      IF (IHTYPE.GT.20) THEN
*
* operator:
* a+i alpha  j beta
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &          NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(1-1)*NGAS) = 1
              IHX_LOC(JGAS+(4-1)*NGAS) = 1
              CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(1,NHX) = 1
                  LABEL(2,NHX) = -1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            END IF
          END DO
        END DO
*
* operator:
* a+i beta a j alpha
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &          NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(2-1)*NGAS) = 1
              IHX_LOC(JGAS+(3-1)*NGAS) = 1
              CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  SIGN_NHX(NHX) = ONE
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  LABEL(1,NHX) = -1
                  LABEL(2,NHX) = 1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            END IF
          END DO
*    quaternionic matrix groups only!
        END DO
      END IF
*
* Two-electron terms with Delta mk .ne. 0
*
**********************************************************
*  Spin-dependent spinor basis. E.g. DIRAC environment   *
*      (bu|bu), (ub|ub) integral classes                 *
**********************************************************
*
      if (IHTYPE.eq.4) then
*
* operator:
*. a+i alpha  a+j alpha  a k beta  a l beta, i.ge.j, k.ge.l
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1))) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS)+1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS)+1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  b | +a  b )
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      IHINDXOP(1,NHX) = 1
                      IHINDXOP(2,NHX) = 2
                      IHINDXOP(3,NHX) = 3
                      IHINDXOP(4,NHX) = 4
                      SIGN_NHX(NHX) = 1 * ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 3
                      IOP_REO(4,NHX) = 2
                      SIGN_OPREO(NHX) = -ONE
                    end if
                  END IF
                ELSE
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS)+1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS)+1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                      INCLUDED_HX(NHX) = 1
* This excitation corresponds to the integral
*  (+a  b | +a  b )
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      IHINDXOP(1,NHX) = 1
                      IHINDXOP(2,NHX) = 2
                      IHINDXOP(3,NHX) = 3
                      IHINDXOP(4,NHX) = 4
                      SIGN_NHX(NHX) = 1 * ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 3
                      IOP_REO(4,NHX) = 2
                      SIGN_OPREO(NHX) = -ONE
                    end if
                END IF
              END DO
            END DO
          END DO
        END DO
*
* operator:
*. a+i beta  a+j beta  a k alpha  a l alpha, i.ge.j, k.ge.l
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1))) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(2-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS)+1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS)+1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP_CC_KRCC(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  a | +b  a)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 2
                      IHINDX(4,NHX) = 4
                      IHINDXOP(1,NHX) = 1
                      IHINDXOP(2,NHX) = 2
                      IHINDXOP(3,NHX) = 3
                      IHINDXOP(4,NHX) = 4
                      SIGN_NHX(NHX) = 1 * ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = 1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 3
                      IOP_REO(4,NHX) = 2
                      SIGN_OPREO(NHX) = -ONE
                    end if
                  END IF
                ELSE
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(2-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS)+1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS)+1
                  CALL CAAB_INCLUDED_KRCC(IHX_LOC,INCLUDED)
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                      INCLUDED_HX(NHX) = 1
* This excitation corresponds to the integral
*  (+b  a | +b  a)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 2
                      IHINDX(4,NHX) = 4
                      IHINDXOP(1,NHX) = 1
                      IHINDXOP(2,NHX) = 2
                      IHINDXOP(3,NHX) = 3
                      IHINDXOP(4,NHX) = 4
                      SIGN_NHX(NHX) = 1 * ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = 1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 3
                      IOP_REO(4,NHX) = 2
                      SIGN_OPREO(NHX) = -ONE
                    end if
                END IF
              END DO
            END DO
          END DO
        END DO
      end if
*
      if (IFLAG.ne.1.and.NTEST.ge.500) then
        do I=1,NHX
          write(6,*) 'I,SIGN_NHX ',I,SIGN_NHX(I)
        end do
      end if
*
      if (IFLAG.ne.1) then
        if (NTEST.ge.100) then
          WRITE(6,*) ' Information from  GET_HX_RELA '
          WRITE(6,*)
          write(6,*) ' Max. no of electrons per GAS:'
          call iwrtma(MXNELGS,1,NGAS,1,NGAS)
          WRITE(6,*) ' The excitation types '
          CALL WRT_SPOX_TP_CC_KRCC(IHX,NHX)
        write(6,*)
        WRITE(6,*) ' Number of excitation types ', NHX
*
        WRITE(6,*) ' Hamiltonian classes included '
        DO I = 1,NHX
          WRITE(6,*) INCLUDED_HX(I), ' Hamiltonian number = ',I 
        END DO
        end if
      end if
        WRITE(6,*) ' Number of excitation types ', NHX
*
      RETURN
      END
*
      SUBROUTINE GETSTRN_GASSM_SPGP_KRCC(ISMFGS,ITPFGS,ISTROC,NSTR,NEL,
     &                                  NNSTSGP,IISTSGP,IUB,
     &                                  WORK,KFREE,LFREE)
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
*
      DIMENSION WORK(*)
*
*. Number of strings per GAS space
C     CALL MEMCHK_KRCC(WORK)
      DO IGAS = 1, NGAS
        NSTFGS(IGAS)  = NNSTSGP(ISMFGS(IGAS),IGAS)
        IBSTFGS(IGAS) = IISTSGP(ISMFGS(IGAS),IGAS)
      END DO
*
      NTEST = 000
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
C         call memchk_KRCC(WORK)
*
          NSTI = NSTFGS(IGAS)
          CALL ADD_STR_GROUP_KRCC(NSTI,IBSTFGS(IGAS),WORK(KOCSTR_ADD),
     &                       NSTB,NSTA,ISTROC,NELB+1,NELI,NEL)
*. Loop over strings in IGAS
        END IF
      END DO
 1001 CONTINUE
      NSTR = NSTRTOT
*
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
C     CALL MEMCHK_KRCC(WORK)
*
      RETURN
      END
*
      SUBROUTINE SOME_DUMB_PRINT_FOR_WORK(AARRAY,NDIM)
* Dumb routine to print work
*
#include "implicit.inc"
*
      DIMENSION AARRAY(*)
*
      WRITE(6,*) ' Will write something from work '
      CALL IWRTMA(AARRAY,1,NDIM,1,NDIM)
*
      RETURN
      END
*
      SUBROUTINE ADD_STR_GROUP_KRCC(NSTADD,IOFADD,ISTADD,NSTB,NSTA,
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
      NTEST = 000
      if (NTEST.ge.1000) then
        write(6,*) 'ADD_STR_GROUP speaking: Input string array'
        call iwrtma(ISTADD,NELADD,IOFADD-1+NSTADD,
     &                     NELADD,IOFADD-1+NSTADD)
      end if
*
      RETURN
      END
*
      SUBROUTINE ABFLIP_SPOXTP(ICAAB_IN,ICAAB_OUT,NGAS)
*
* Obtain spin-orbital excitation type ICAAB_OUT by
* spinflipping ICAAB_IN
*
* Jeppe Olsen, July 11, 2001
*
#include "implicit.inc"
*. Input
      INTEGER ICAAB_IN(NGAS,4)
*. Output
      INTEGER ICAAB_OUT(NGAS,4)
*
      CALL ICOPVE(ICAAB_IN(1,1),ICAAB_OUT(1,2),NGAS)
      CALL ICOPVE(ICAAB_IN(1,2),ICAAB_OUT(1,1),NGAS)
      CALL ICOPVE(ICAAB_IN(1,3),ICAAB_OUT(1,4),NGAS)
      CALL ICOPVE(ICAAB_IN(1,4),ICAAB_OUT(1,3),NGAS)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Spin-orbital excitation type '
        CALL WRT_SPOX_TP_CC_KRCC(ICAAB_IN,1)
        WRITE(6,*) ' Spinflipped orbital excitation type '
        CALL WRT_SPOX_TP_CC_KRCC(ICAAB_OUT,1)
      END IF
*
      RETURN
      END
*
      SUBROUTINE GET_REF_ALBE_OCC_KRCC(IREFSPC,IREF_AL,IREF_BE,IMS2_REF)
*
* Obtain alpha and beta occupations for reference space
*
* Reference space is assumed to be a single pair of occupations of
* alpha and beta strings ( this includes closed shell HF,
* Highspin open shell and CAS reference)
*
* Only a single valence orbital space is assumed
*
* Jeppe Olsen, March 2000
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "strinp.inc"
*. Output : Alpha and beta occupations for each GAS space
      INTEGER IREF_AL(NGAS),IREF_BE(NGAS)
*
*. Total number of Hole orbitals
      NHOLE = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.1) THEN
          NHOLE = NHOLE + NGSOBT(IGAS)
        END IF
      END DO
*. Number of orbitals in valence space
      NELEC_AL = NELEC(1)
      NELEC_BE = NELEC(2)
      NVAL_AL = NELEC_AL - NHOLE
      NVAL_BE = NELEC_BE - NHOLE
*
      DO IGAS = 1, NGAS
        NORB = NGSOBT(IGAS)
        IF(IHPVGAS(IGAS).EQ.1) THEN
          IREF_AL(IGAS) = NORB
          IREF_BE(IGAS) = NORB
        ELSE IF( IHPVGAS(IGAS).EQ.2) THEN
          IREF_AL(IGAS) = 0
          IREF_BE(IGAS) = 0
        ELSE IF( IHPVGAS(IGAS).EQ.3) THEN
          IREF_AL(IGAS) = NVAL_AL
          IREF_BE(IGAS) = NVAL_BE
        END IF
      END DO
*
* Added by Lasse
*
      NAEL_TOT = 0
      NBEL_TOT = 0
      DO IGAS = 1,NGAS
         NAEL_TOT = IREF_AL(IGAS) + NAEL_TOT
         NBEL_TOT = IREF_BE(IGAS) + NBEL_TOT
      END DO
*
      IMS2_REF = NAEL_TOT - NBEL_TOT
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Alpha and Beta Occupation for reference space '
        WRITE(6,*)
        CALL IWRTMA(IREF_AL,1,NGAS,1,NGAS)
        CALL IWRTMA(IREF_BE,1,NGAS,1,NGAS)
        WRITE(6,*) ' Total spinprojection for reference',IMS2_REF
      END IF
*
      RETURN
      END
*
      SUBROUTINE GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP,NIGRP,ISPGRPSM,NEL,
     &                                  NSTR,ISTR,
     &                                  NORBT,IDOREO,IZ,IREO,
     &                                  WORK,KFREE,LFREE)
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
      DIMENSION WORK(*)
*
      CALL QENTER('GETST')
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================================= '
        WRITE(6,*) ' Welcome to GETSTR_TOTSM_SPGP_KRCC '
        WRITE(6,*) ' ================================= '
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
C     CALL MEMCHK_KRCC(WORK)
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
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' Write out NNSTSGP '
        CALL IWRTMA(NNSTSGP,1,NSMST,1,NGAS)
        DO I=1,NSMST
         DO J=1,NGAS
           print*,'NNSTSGP,I,J',NNSTSGP(I,J),I,J
         END DO
        END DO
        call memchk_KRCC(WORK)
      END IF
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
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Type of each GAS space '
         CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of elecs per GAS space '
         CALL IWRTMA(NELFGS,1,NGAS,1,NGAS)
         call memchk_KRCC(WORK)
      END IF
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
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL SYMCOM_KRCC(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE NGASL
        CALL SYMCOM_KRCC(2,1,ISTSMM1,ISMGSN,ISPGRPSM)
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
          call memchk_KRCC(WORK)
        END IF
*. Obtain all strings of this symmetry
        CALL GETSTRN_GASSM_SPGP_KRCC(ISMFGS,ITPFGS,ISTR(1,ISTRBS),
     &                              NSTR,NEL,
     &                              NNSTSGP,IISTSGP,IUB,
     &                              WORK,KFREE,LFREE)
C       print*,'right arfter'
C     CALL MEMCHK_KRCC(WORK)
C     print*,'jijiji'
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
*
C     CALL MEMCHK_KRCC(WORK)
      CALL QEXIT('GETST')
      RETURN
      END
*
      SUBROUTINE CAAB_INCLUDED_KRCC(IEXOP,INCLUDED_CAAB)
*
* Is excitation operator CAAB included in the
* list of excitations
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
*. Specific input
      INTEGER IEXOP(NGAS,4)
*
      INCLUDED_CAAB = 1
* Creation alpha
      CALL OCC_INCLUDED_KRCC(IEXOP(1,1),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
* Creation beta
      CALL OCC_INCLUDED_KRCC(IEXOP(1,2),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
* Annihilation alpha
      CALL OCC_INCLUDED_KRCC(IEXOP(1,3),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
* Annihilation beta
      CALL OCC_INCLUDED_KRCC(IEXOP(1,4),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' CAAB excitation operator '
       CALL WRT_SPOX_TP_CC_KRCC(IEXOP,1)
       WRITE(6,*) ' INCLUDED_CAAB = ', INCLUDED_CAAB
      END IF
*
      RETURN
      END
*
      SUBROUTINE OCC_INCLUDED_KRCC(IOCC,INCLUDED)
*
* A supergroup is given as occupation for each gas space
*
* Check if all these groups are included, and return
* answer in INCLUDED
*
* Jeppe and Timo, Oct 2001
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
*. Specific input
      INTEGER IOCC(NGAS)
*
      NTEST = 00
*
      INCLUDED = 1
      DO IOBTP = 1, NGAS
        JJGRP = 0
        DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
          IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
        END DO
        IF(JJGRP.EQ.0 ) INCLUDED = 0
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation  '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Included = ',INCLUDED
      END IF
*
      RETURN
      END
*
      SUBROUTINE SORT_OP_MASTER(WORK,KFREE,LFREE)
*
* Routine for sorting all operators according to N,Mk,Mub
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "lucinp.inc"
#include "ctcc.inc"
#include "ctccp.inc"
#include "cgas.inc"
*
       DIMENSION WORK(*)
*
       IZERO = 0
*
* Choosing same dimension for all cluster and intermediate operators
*
      IMAXKRFLIP = NACTEL !MIN(3*MK2DEL,MX_EXC_LEVEL)
      IMUBMAX = NACTEL !MX_EXC_LEVEL
      IDIMINDX = (MX_EXC_LEVEL+1)*(2*IMAXKRFLIP+1)*(2*IMUBMAX+1)
      print*,'IDIMINDX',IDIMINDX
*
* First cluster operator (could perhaps be done before (see if this
* needs to mve later on)
      CALL MEMGET('INTE',KT_IDXS,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KT_IDXS),IZERO,IDIMINDX)
      CALL MEMGET('INTE',KT_IDXF,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KT_IDXF),IZERO,IDIMINDX)
*
* Allocate help arrays
*
      CALL MEMGET('INTE',KHELP_OP,4*NGAS*NSPOBEX_TPE,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHELP_IDX,NSPOBEX_TPE,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KT_IDXH,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KT_IDXH),IZERO,IDIMINDX)
*
      CALL SORT_OP_3_INDEX(WORK(KLSOBEX_CC),NSPOBEX_TPE,NGAS,
     &                     MX_EXC_LEVEL,IMAXKRFLIP,IMUBMAX,
     &                     WORK(KT_IDXS),WORK(KT_IDXF),
     &                     WORK(KT_IDXH),WORK(KHELP_OP),WORK(KHELP_IDX))
*
* release memory
*
      CALL MEMREL('T_IDX',WORK,KHELP_OP,KHELP_OP,KFREE,LFREE)
*
* Sort intermediate 01
*
      CALL MEMGET('INTE',KINT01_IDXS,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT01_IDXS),IZERO,IDIMINDX)
      CALL MEMGET('INTE',KINT01_IDXF,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT01_IDXF),IZERO,IDIMINDX)
*
* Allocate help arrays
*
      CALL MEMGET('INTE',KHELP_OP,4*NGAS*NINTER01,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHELP_IDX,NINTER01,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KT_IDXH,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KT_IDXH),IZERO,IDIMINDX)
*
      CALL SORT_OP_3_INDEX(WORK(KINTM01),NINTER01,NGAS,
     &                     MX_EXC_LEVEL,IMAXKRFLIP,IMUBMAX,
     &                     WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                     WORK(KT_IDXH),WORK(KHELP_OP),WORK(KHELP_IDX))
*
* release memory
*
      CALL MEMREL('INT01_IDX',WORK,KHELP_OP,KHELP_OP,KFREE,LFREE)
*
* Sort intermediate 10
*
      CALL MEMGET('INTE',KINT10_IDXS,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT10_IDXS),IZERO,IDIMINDX)
      CALL MEMGET('INTE',KINT10_IDXF,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT10_IDXF),IZERO,IDIMINDX)
*
* Allocate help arrays
*
      CALL MEMGET('INTE',KHELP_OP,4*NGAS*NINTER10,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHELP_IDX,NINTER10,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KT_IDXH,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KT_IDXH),IZERO,IDIMINDX)
*
      CALL SORT_OP_3_INDEX(WORK(KINTM10),NINTER10,NGAS,
     &                     MX_EXC_LEVEL,IMAXKRFLIP,IMUBMAX,
     &                     WORK(KINT10_IDXS),WORK(KINT10_IDXF),
     &                     WORK(KT_IDXH),WORK(KHELP_OP),WORK(KHELP_IDX))
*
* release memory
*
      CALL MEMREL('INT10_IDX',WORK,KHELP_OP,KHELP_OP,KFREE,LFREE)
*
* Sort intermediate 11
*
      CALL MEMGET('INTE',KINT11_IDXS,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT11_IDXS),IZERO,IDIMINDX)
      CALL MEMGET('INTE',KINT11_IDXF,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT11_IDXF),IZERO,IDIMINDX)
*
* Allocate help arrays
*
      CALL MEMGET('INTE',KHELP_OP,4*NGAS*NINTER11,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHELP_IDX,NINTER11,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KT_IDXH,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KT_IDXH),IZERO,IDIMINDX)
*
      print*,'sort m11'
      CALL SORT_OP_3_INDEX(WORK(KINTM11),NINTER11,NGAS,
     &                     MX_EXC_LEVEL,IMAXKRFLIP,IMUBMAX,
     &                     WORK(KINT11_IDXS),WORK(KINT11_IDXF),
     &                     WORK(KT_IDXH),WORK(KHELP_OP),WORK(KHELP_IDX))
*
* release memory
*
      CALL MEMREL('INT11_IDX',WORK,KHELP_OP,KHELP_OP,KFREE,LFREE)
*
* Sort intermediate 02
*
      CALL MEMGET('INTE',KINT02_IDXS,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT02_IDXS),IZERO,IDIMINDX)
      CALL MEMGET('INTE',KINT02_IDXF,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT02_IDXF),IZERO,IDIMINDX)
*
* Allocate help arrays
*
      CALL MEMGET('INTE',KHELP_OP,4*NGAS*NINTER02,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHELP_IDX,NINTER02,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KT_IDXH,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KT_IDXH),IZERO,IDIMINDX)
*
      print*,'sort m02'
      CALL SORT_OP_3_INDEX(WORK(KINTM02),NINTER02,NGAS,
     &                     MX_EXC_LEVEL,IMAXKRFLIP,IMUBMAX,
     &                     WORK(KINT02_IDXS),WORK(KINT02_IDXF),
     &                     WORK(KT_IDXH),WORK(KHELP_OP),WORK(KHELP_IDX))
*
* release memory
*
      CALL MEMREL('INT02_IDX',WORK,KHELP_OP,KHELP_OP,KFREE,LFREE)
*
* Sort intermediate 20
*
      CALL MEMGET('INTE',KINT20_IDXS,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT20_IDXS),IZERO,IDIMINDX)
      CALL MEMGET('INTE',KINT20_IDXF,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT20_IDXF),IZERO,IDIMINDX)
*
* Allocate help arrays
*
      CALL MEMGET('INTE',KHELP_OP,4*NGAS*NINTER20,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHELP_IDX,NINTER20,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KT_IDXH,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KT_IDXH),IZERO,IDIMINDX)
*
      print*,'sort m20'
      CALL SORT_OP_3_INDEX(WORK(KINTM20),NINTER20,NGAS,
     &                     MX_EXC_LEVEL,IMAXKRFLIP,IMUBMAX,
     &                     WORK(KINT20_IDXS),WORK(KINT20_IDXF),
     &                     WORK(KT_IDXH),WORK(KHELP_OP),WORK(KHELP_IDX))
*
* release memory
*
      CALL MEMREL('INT20_IDX',WORK,KHELP_OP,KHELP_OP,KFREE,LFREE)
*
* Sort intermediate 12
*
      CALL MEMGET('INTE',KINT12_IDXS,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT12_IDXS),IZERO,IDIMINDX)
      CALL MEMGET('INTE',KINT12_IDXF,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KINT12_IDXF),IZERO,IDIMINDX)
*
* Allocate help arrays
*
      CALL MEMGET('INTE',KHELP_OP,4*NGAS*NINTER12,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHELP_IDX,NINTER12,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KT_IDXH,IDIMINDX,WORK,KFREE,LFREE)
      CALL ISETVC(WORK(KT_IDXH),IZERO,IDIMINDX)
      print*,'sort m12'
*
      CALL SORT_OP_3_INDEX(WORK(KINTM12),NINTER12,NGAS,
     &                     MX_EXC_LEVEL,IMAXKRFLIP,IMUBMAX,
     &                     WORK(KINT12_IDXS),WORK(KINT12_IDXF),
     &                     WORK(KT_IDXH),WORK(KHELP_OP),WORK(KHELP_IDX))
*
* release memory
*
      CALL MEMREL('INT12_IDX',WORK,KHELP_OP,KHELP_OP,KFREE,LFREE)
*
      RETURN
      END
*
      SUBROUTINE HTYPE_PH_CONTRACTIONS_MASTER(WORK,KFREE,LFREE)
*
* Master routine to find the number of holes and particles for given
* HTYPE
*
#include "implicit.inc"
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "ctcc.inc"
#include "ctccp.inc"
      DIMENSION WORK(*)
*
* Initialize
*
      NH22 = 0
      NH21 = 0
      NH12 = 0
      NH20 = 0
      NH11 = 0
      NH02 = 0
      NH10 = 0
      NH01 = 0
      NH00 = 0
*
      INI = 0
      CALL HTYPE_PH_CONTRACTIONS(INI,WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &                                 NH22,NH21,NH12,NH20,
     &                                 NH11,NH02,NH10,NH01,NH00,
     &                                 IDUMMY,IDUMMY,
     &                                 IDUMMY,IDUMMY,
     &                                 IDUMMY,IDUMMY,
     &                                 IDUMMY,IDUMMY,
     &                                 IDUMMY)
*
      CALL MEMGET('INTE',KNH22,NH22,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNH21,NH21,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNH12,NH12,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNH20,NH20,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNH11,NH11,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNH02,NH02,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNH10,NH10,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNH01,NH01,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNH00,NH00,WORK,KFREE,LFREE)
*
      INI = 1
      CALL HTYPE_PH_CONTRACTIONS(INI,WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &                                 NH22,NH21,NH12,NH20,
     &                                 NH11,NH02,NH10,NH01,NH00,
     &                                 WORK(KNH22),WORK(KNH21),
     &                                 WORK(KNH12),WORK(KNH20),
     &                                 WORK(KNH11),WORK(KNH02),
     &                                 WORK(KNH10),WORK(KNH01),
     &                                 WORK(KNH00))
*
      CALL MEMCHK_KRCC(WORK)
*
      RETURN
      END
*
      SUBROUTINE HTYPE_PH_CONTRACTIONS(INI,IHTP,NHTP,NGAS,
     &                                 NH22,NH21,NH12,NH20,
     &                                 NH11,NH02,NH10,NH01,NH00,
     &                                 INH22,INH21,INH12,INH20,
     &                                 INH11,INH02,INH10,INH01,INH00)
*
#include "implicit.inc"
*
* Find the number of holes and particles a given HTYPE has.
* Extended to also include an active space
*
      INTEGER IHTP(4*NGAS,NHTP)
      INTEGER INH22(NH22), INH21(NH21), INH12(NH12), INH20(NH20)
      INTEGER INH11(NH11), INH02(NH02), INH10(NH10), INH01(NH01)
      INTEGER INH00(NH00)
*
* Initialize
*
      NH22 = 0
      NH21 = 0
      NH12 = 0
      NH20 = 0
      NH11 = 0
      NH02 = 0
      NH10 = 0
      NH01 = 0
      NH00 = 0
*
      DO ITP = 1,NHTP
        CALL RANK_FOR_CAAB_KRCC(IHTP(1,ITP),NEX,NDEEX,
     &                          NHOLE,NPART,NHOLEV,NPARTV)
* This HTYPE will now have a minimum of NPART,NHOLE and maximum of
* NPART+NPARTV,NHOLE+NHOLEV indices to be contracted
        DO IPART = NPART,MIN(NPART+NPARTV,2)
          DO IHOLE = NHOLE,MIN(NHOLE+NHOLEV,2)
            IF(IPART.EQ.2.AND.IHOLE.EQ.2) THEN
              NH22 = NH22 + 1
              IF(INI.EQ.1) THEN
                INH22(NH22) = ITP
              END IF
            ELSE IF(IPART.EQ.2.AND.IHOLE.EQ.1) THEN
              NH21 = NH21 + 1
              IF(INI.EQ.1) THEN
                INH21(NH21) = ITP
              END IF
            ELSE IF(IPART.EQ.1.AND.IHOLE.EQ.2) THEN
              NH12 = NH12 + 1
              IF(INI.EQ.1) THEN
                INH12(NH12) = ITP
              END IF
            ELSE IF(IPART.EQ.2.AND.IHOLE.EQ.0) THEN
              NH20 = NH20 + 1
              IF(INI.EQ.1) THEN
                INH20(NH20) = ITP
              END IF
            ELSE IF(IPART.EQ.1.AND.IHOLE.EQ.1) THEN
              NH11 = NH11 + 1
              IF(INI.EQ.1) THEN
                INH11(NH11) = ITP
              END IF
            ELSE IF(IPART.EQ.0.AND.IHOLE.EQ.2) THEN
              NH02 = NH02 + 1
              IF(INI.EQ.1) THEN
                INH02(NH02) = ITP
              END IF
            ELSE IF(IPART.EQ.1.AND.IHOLE.EQ.0) THEN
              NH10 = NH10 + 1
              IF(INI.EQ.1) THEN
                INH10(NH10) = ITP
              END IF
            ELSE IF(IPART.EQ.0.AND.IHOLE.EQ.1) THEN
              NH01 = NH01 + 1
              IF(INI.EQ.1) THEN
                INH01(NH01) = ITP
              END IF
            ELSE IF(IPART.EQ.0.AND.IHOLE.EQ.0) THEN
              NH00 = NH00 + 1
              IF(INI.EQ.1) THEN
                INH00(NH00) = ITP
              END IF
            END IF
          END DO
        END DO
      END DO
*
      NHTOT=NH22+NH21+NH12+NH20+NH11+NH02+NH10+NH01+NH00
      IF(NHTOT.NE.NHTP) STOP ' Error in HTYPE_PH_CONTRACTIONS '
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of Hamiltonian Types ',NHTOT
        WRITE(6,*) ' 2 Particles 2 Holes ',NH22
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The Specific types are '
          DO I = 1,NH22
            WRITE(6,*) INH22(I)
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,INH22(I)),1)
          END DO
        END IF
        WRITE(6,*) ' 2 Particles 1 Holes ',NH21
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The Specific types are '
          DO I = 1,NH21
            WRITE(6,*) INH21(I)
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,INH21(I)),1)
          END DO
        END IF
        WRITE(6,*) ' 1 Particles 2 Holes ',NH12
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The Specific types are '
          DO I = 1,NH12
            WRITE(6,*) INH12(I)
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,INH12(I)),1)
          END DO
        END IF
        WRITE(6,*) ' 2 Particles 0 Holes ',NH20
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The Specific types are '
          DO I = 1,NH20
            WRITE(6,*) INH20(I)
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,INH20(I)),1)
          END DO
        END IF
        WRITE(6,*) ' 1 Particles 1 Holes ',NH11
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The Specific types are '
          DO I = 1,NH11
            WRITE(6,*) INH11(I)
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,INH11(I)),1)
          END DO
        END IF
        WRITE(6,*) ' 0 Particles 2 Holes ',NH02
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The Specific types are '
          DO I = 1,NH02
            WRITE(6,*) INH02(I)
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,INH02(I)),1)
          END DO
        END IF
        WRITE(6,*) ' 1 Particles 0 Holes ',NH10
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The Specific types are '
          DO I = 1,NH10
            WRITE(6,*) INH10(I)
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,INH10(I)),1)
          END DO
        END IF
        WRITE(6,*) ' 0 Particles 1 Holes ',NH01
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The Specific types are '
          DO I = 1,NH01
            WRITE(6,*) INH01(I)
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,INH01(I)),1)
          END DO
        END IF
        WRITE(6,*) ' 0 Particles 0 Holes ',NH00
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The Specific types are '
          DO I = 1,NH00
            WRITE(6,*) INH00(I)
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,INH00(I)),1)
          END DO
        END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE RANK_FOR_CAAB_KRCC(IHTP,NEX,NDEEX,
     &                              NHOLE,NPART,NHOLEV,NPARTV)
*
* An operator IHTP is given in CAAB form. 
* Obtain number of excitation and deexcitation operators
*
* IHPVGAS_AB is used to determine the hp nature of spinorbitals
*
* Jeppe Olsen, March 2003
*
* Extended to include NHOLE and NPART and Valence space 2010 Lasse
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
*. Specific input
      INTEGER IHTP(NGAS,4)
*
      NEX = 0
      NHOLE = 0
      NPART = 0
      NHOLEV = 0
      NPARTV = 0
*
      DO ICA = 1, 2
        DO IAB = 1, 2
          ICAAB = (ICA-1)*2+ IAB
          DO IGAS = 1, NGAS
            IF((ICA.EQ.1.AND.IHPVGAS_AB(IGAS,IAB).EQ.2).OR.
     &         (ICA.EQ.2.AND.IHPVGAS_AB(IGAS,IAB).EQ.1)    ) THEN
                NEX = NEX  + IHTP(IGAS,ICAAB)
            ELSE IF(ICA.EQ.2.AND.IHPVGAS_AB(IGAS,IAB).EQ.2) THEN
                NPART = NPART + IHTP(IGAS,ICAAB)
            ELSE IF(ICA.EQ.1.AND.IHPVGAS_AB(IGAS,IAB).EQ.1) THEN
                NHOLE = NHOLE + IHTP(IGAS,ICAAB)
            ELSE IF(ICA.EQ.2.AND.IHPVGAS_AB(IGAS,IAB).EQ.3) THEN
                NPARTV = NPARTV + IHTP(IGAS,ICAAB)
            ELSE IF(ICA.EQ.1.AND.IHPVGAS_AB(IGAS,IAB).EQ.3) THEN
                NHOLEV = NHOLEV + IHTP(IGAS,ICAAB)
            END IF
          END DO
        END DO
      END DO
*
      NDEEX = NPART + NHOLE + NPARTV + NHOLEV
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CAAB operator : '
        CALL WRT_SPOX_TP_CC_KRCC(IHTP,1)
        WRITE(6,*) ' Number of excitation and max deexcitation ops : ',
     &  NEX,NDEEX
        WRITE(6,*) ' The certain deexcitation ops  '
        WRITE(6,*) ' The Holes                     ',NHOLE
        WRITE(6,*) ' The Particles                 ',NPART
        WRITE(6,*) ' The possible deexcitation ops '
        WRITE(6,*) ' The Valence Holes             ',NHOLEV
        WRITE(6,*) ' The Valence Particles         ',NPARTV
      END IF
*
      RETURN
      END
*
      SUBROUTINE FIND_ALL_ABOBEX(INI,NGAS,IOBEX_TP,
     &                           NOBEX_TP,IOBEX_TP_KRCC,NOBEX_TP_KRCC)
*
*
#include "implicit.inc"
*
*. Input
      INTEGER IOBEX_TP(2*NGAS,NOBEX_TP)
*. Output
      INTEGER IOBEX_TP_KRCC(2*NGAS,*)
*. Scratch
      INTEGER ITEMP(2*NGAS,NOBEX_TP)
*
      NTEST = 00
*
* First find number of different creation and annihilation types
*
      CALL FIND_CA_ABOBEX(INI,NGAS,IOBEX_TP,NOBEX_TP,ITEMP,
     &                    NOBEX_TP_KRCC,NC,NA)
*
* Now construct the OBEX for the KRCC
*
      IF(INI.EQ.1) THEN
*
* First copy alpha creators and annihilators
*
        ICOUNT = 0
*
C       DO I=1,4
C         WRITE(6,*) (ITEMP(J,I),J=1,NGAS)
C         WRITE(6,*) (ITEMP(J+NGAS,I),J=1,NGAS)
C       END DO
        
        DO I = 1,NC
          DO J = 1,NA
            ICOUNT = ICOUNT + 1
* Copy creators
C           WRITE(6,*) (ITEMP(K,I),K=1,NGAS)
            CALL ICOPVE(ITEMP(1,I),IOBEX_TP_KRCC(1,ICOUNT),NGAS)
C           WRITE(6,*) (IOBEX_TP_KRCC(K,ICOUNT),K=1,NGAS)
* Copy annihilators
C           WRITE(6,*) (ITEMP(K+NGAS,J),K=1,NGAS)
            CALL ICOPVE(ITEMP(1+NGAS,J),IOBEX_TP_KRCC(1+NGAS,ICOUNT),
     &                  NGAS)
C           WRITE(6,*) (IOBEX_TP_KRCC(K+NGAS,ICOUNT),K=1,NGAS)
          END DO
        END DO
      END IF
*
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Number of combinations ',NOBEX_TP_KRCC
        WRITE(6,*) ' From which ',NC,' is creators and ',NA,' anni'
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The combinations '
          DO I =1,NOBEX_TP_KRCC
            WRITE(6,*) ' Combination ',I
            WRITE(6,*) (IOBEX_TP_KRCC(J,I),J=1,NGAS)
            WRITE(6,*) (IOBEX_TP_KRCC(J+NGAS,I),J=1,NGAS)
          END DO
        END IF
      END IF
*
*      
      RETURN
      END
*
      SUBROUTINE FIND_CA_ABOBEX(INI,NGAS,IABOBEX,NABOBEX,IOUT,
     &                          NDIM,NC,NA)
*
* Finds the number of unique types of creators and annihilators from the
* ABOBEX types
*
#include "implicit.inc"
*
*. Input
      INTEGER IABOBEX(2*NGAS,NABOBEX)
*. Output
      INTEGER IOUT(2*NGAS,NABOBEX)

      INTEGER ISCR(2,NABOBEX)
*
      NTEST = 00
*
      CALL ISETVC(ISCR,1,2*NABOBEX)
*
      ICOUNT = 1
*
* First check creators (in case there is only two we assume they are
* different since one will be the identity)
*      
      DO
        ICOUNT = ICOUNT + 1
        IF(ICOUNT.GT.NABOBEX) EXIT
        IF(ISCR(1,ICOUNT-1).EQ.0) CYCLE
        DO I =ICOUNT,NABOBEX
          IF(ISCR(1,I).EQ.0) CYCLE
* Compare
          CALL COMPARE_OP(IABOBEX(1,I),IABOBEX(1,ICOUNT-1),
     &                           NGAS,IDIFF)
          IF(IDIFF.EQ.0) THEN
            ISCR(1,I) = 0
          END IF        
        END DO
      END DO
*
* Now annihilators
*
      ICOUNT = 1
*      
      DO
        ICOUNT = ICOUNT + 1
        IF(ICOUNT.GT.NABOBEX) EXIT
        IF(ISCR(2,ICOUNT-1).EQ.0) CYCLE
        DO I =ICOUNT,NABOBEX
          IF(ISCR(2,I).EQ.0) CYCLE
* Compare
          CALL COMPARE_OP(IABOBEX(1+NGAS,I),IABOBEX(1+NGAS,ICOUNT-1),
     &                           NGAS,IDIFF)
          IF(IDIFF.EQ.0) THEN
            ISCR(2,I) = 0
          END IF
        END DO
      END DO
*
      NC = 0
      NA = 0
*
      DO I =1,NABOBEX
        IF(ISCR(1,I).EQ.1) THEN
          NC = NC + 1
          IF(INI.EQ.1) THEN
            CALL ICOPVE(IABOBEX(1,I),IOUT(1,NC),NGAS)
          END IF        
        END IF
        IF(ISCR(2,I).EQ.1) THEN
          NA = NA + 1
          IF(INI.EQ.1) THEN
            CALL ICOPVE(IABOBEX(1+NGAS,I),IOUT(1+NGAS,NA),NGAS)
          END IF
        END IF
      END DO
*
      NDIM = NC*NA
*
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Total dimension of combination ',NDIM
        WRITE(6,*) '================================'
        WRITE(6,*) ' Number of unique creators ',NC
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The creator types '
          DO I =1,NC
            WRITE(6,*) (IABOBEX(J,I),J=1,NGAS)
          END DO
        END IF
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The creator types '
          DO I =1,NC
            WRITE(6,*) (IOUT(J,I),J=1,NGAS)
          END DO
        END IF
        WRITE(6,*) ' Number of unique annihilators ',NA
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The annihilator types '
          DO I =1,NA
            WRITE(6,*) (IABOBEX(J+NGAS,I),J=1,NGAS)
          END DO
        END IF  
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' The annihilator types '
          DO I =1,NA
            WRITE(6,*) (IOUT(J+NGAS,I),J=1,NGAS)
          END DO
        END IF  
      END IF
*      
      RETURN
      END
*
      SUBROUTINE COMPARE_OP(IOP,JOP,NGAS,IDIFF)
*
* Compare two operators
*
      INTEGER IOP(NGAS),JOP(NGAS)
*
      IDIFF = 0
*      
      DO IGAS =1,NGAS
        IF(IOP(IGAS).NE.JOP(IGAS)) THEN
          IDIFF = 1
          EXIT
        END IF        
      END DO
*
      RETURN
      END
*
      SUBROUTINE STRING_IN_CCEXCIT_REL(ISTRING,ICAAB,INSPC,
     &                             NSPOBEX_TP,ISPOBEX_TP)
*
* Is String ISTRING ( given as occ in each gas space)
* include as CAAB string ICAAB in any of the NSPOBEX_TP TCC excitations 
* blocks given in ISPOBEX_TP 
*
* Jeppe Olsen, May 2000
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER ISTRING(NGAS)
*
      INSPC = 0
      IB = 1 + (ICAAB-1)*NGAS
      DO KTP = 1, NSPOBEX_TP
        IDEL = 0
        DO IGAS = 1, NGAS
          IDEL = IDEL + ABS(ISPOBEX_TP(IB-1+IGAS,KTP)-ISTRING(IGAS))
        END DO
        IF(IDEL.EQ.0) THEN
          INSPC = 1
          GOTO 1001 
        END IF
      END DO
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ICAAB of string = ', ICAAB
        WRITE(6,*) ' Occupation : '
        CALL IWRTMA(ISTRING,1,NGAS,1,NGAS)
        IF(INSPC.EQ.1) THEN
           WRITE(6,*) ' String is in space '
        ELSE
           WRITE(6,*) ' String is out of space '
        END IF
      END IF
*
      RETURN 
      END
