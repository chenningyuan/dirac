      SUBROUTINE STRTYP_GAS(IPRNT)
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
*.
*. Occupation classes corresponding to largest CI space
*
      CALL QENTER('OCCLS')
      CALL OCCLS(1,NOCCLS,IOCCLS,NACTEL,NGAS,IGSOCC(1,1),IGSOCC(1,2))
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
*. Beta terms
          ISTTYP = ISTTYP + 1
          NELEC(ISTTYP) = NBELRF-IANNI
          NELFTP(ISTTYP) = NBELRF-IANNI
          NELFTP_REF(ISTTYP) = NBELRF
          IF(IANNI.EQ.0) IST_FOR_DT(2,IDTTP) = ISTTYP
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
