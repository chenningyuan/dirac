      SUBROUTINE ADSTN_GAS(IOBSM,IOBTP,ISPGP,ISPGPSM,ISPGPTP,
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
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
      PARAMETER(LOFFI=8*8*8*8*8*8*8)
      DIMENSION IOFFI(LOFFI)
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
      CALL QENTER('ADSTN ')
*
      IF(NGAS.GT.MXLNGAS) THEN
        WRITE(6,*) ' Ad hoc programming in ADSTN (IOFFI)'
        WRITE(6,*) ' Must be changed - or redimensioned '
        Call Abend2( 'ADST : IOFFI problem ' )
      END IF
*
      NTEST = 00
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
      CALL NEWTYP(ISPGRPABS,1,IOBTP,1,KSPGRPABS)
      CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
CSTART
      IF(IUB.EQ.1) THEN
         NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
      ELSE
         NKSTR = NSTFSMSPGP2(KSM,KSPGRPABS)
      END IF
CEND
      IF(NTEST.GE.200) WRITE(6,*)
     & ' KSM, KSPGPRABS, NKSTR : ', KSM,KSPGRPABS, NKSTR
      IF(NKSTR.EQ.0) GOTO 9999
*
      IF(IUB.EQ.1) THEN
        NORBTS= NOBPTS(IOBTP,IOBSM)
      ELSE
        NORBTS= NOBPTS2(IOBTP,IOBSM)
      END IF
      ZERO =0.0D0
      CALL SETVEC(XI1S,ZERO,NORBTS*NKSTR)
      IZERO = 0
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
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
          CALL SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*.  sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISPGPSM)
        ISMFGS(NGASL) = ISMGSN
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
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
        IOFFI(IOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' IOFF, IOFFI(IOFF) NSTRII ',
     &                 IOFF, IOFFI(IOFF),NSTRII
        END IF
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
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,KSM)
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
          NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IOBTP)
        CALL SYMCOM(3,1,IOBSM,ISMFGS(IOBTP),IACSM)
        ISMFGS(IOBTP) = IACSM
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGAS
          IOFF = IOFF + (ISMFGS(IGAS)-1)*MULT
          MULT = MULT * NSMST
        END DO
        ISMFGS(IOBTP) = ISAVE
        IBSTRINI = IOFFI(IOFF)
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
*
      CALL QEXIT('ADSTN ')
      RETURN
      END
