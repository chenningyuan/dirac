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
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
*. Local scratch
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*. save internally
      COMMON/SSAVE/NELIS(4), NSTRKS(4)
*
* =======
*. Output
* =======
*
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*
      CALL QENTER('ADST  ')
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' ADST_GAS in service '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' ISM ISPGP ITP',ISM,ISPGP,ITP
      END IF
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
        Call Abend()
      end if
*
      IF(NTEST.GE.100) WRITE(6,*)
     &' ADST : IORB IORBTP IORBSM : ', IORB,IORBTP,IORBSM
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      CALL NEWTYP(ISPGPABS,1,IORBTP,1,KSPGPABS)
      KTP = ITP + 2
      KSPGPREL = KSPGPABS - IBSPGPFTP(KTP) + 1
*.
      CALL SYMCOM(2,0,IORBSM,KSM,ISM)
      IF(NTEST.GE.100) WRITE(6,*)
     & ' KSM and KSPGPABS : ', KSM,KSPGPABS
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
         CALL GETSTR_TOTSM_SPGP(ITP,ISPGP,ISM,NELI,NSTRI,
     &                          WORK(KLLOC),NOCOB,
     &                          1,WORK(KLLZ),WORK(KLLREO))
C?       CALL QEXIT('GENSTR')
       END IF
       NELK = NELIS(I12) - 1
C?     WRITE(6,*) ' ADST : NELK, NELI ', NELK,NELI
       IF(KFRST.NE.0) THEN
*. Occupation of K strings
CT       CALL QENTER('GENSTR')
         CALL GETSTR_TOTSM_SPGP(KTP,KSPGPREL,KSM,NELK,NSTRK,
     &                          WORK(KLLOC),NOCOB,
     &                          0,ISUM,IDUM)
CT       CALL QEXIT('GENSTR')
         NSTRKS(K12) = NSTRK
      END IF
*. a+i !kstr> = +/- !istr>
      NSTRK = NSTRKS(K12)
      IF(NTEST.GE.100) WRITE(6,*) ' NSTRK = ', NSTRK
      CALL ADS1_GAS(NK,I1,XI1S,LI1,IORB,NORB,
     &          WORK(KLLOC),NELK,NSTRK,WORK(KLLREO),WORK(KLLZ),
     &          NOCOB,KMAX,KMIN,IEND,KACT,SCLFAC)
*
C     SAVE NELK , NSTRK
      CALL QEXIT('ADST  ')
      RETURN
      END
