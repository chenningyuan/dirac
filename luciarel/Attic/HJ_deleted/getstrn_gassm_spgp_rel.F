      SUBROUTINE GETSTRN_GASSM_SPGP_REL(ISMFGS,ITPFGS,ISTROC,NSTR,NEL,
     &                          NNSTSGP,IISTSGP,IUB)
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
*
      RETURN
      END
