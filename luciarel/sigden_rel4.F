      SUBROUTINE SIGDEN_REL4(NSBLOCK,ISBLOCK,CB,SB,C2,
     &           NSSOA_S,NSSOB_S,NSSOA_C,NSSOB_C,
     &           IASTP,IBSTP,IACTP,IBCTP,MAXK,
     &           XINT,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           IPRS,IPRD,LUC,CJRES,SIRES,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,
     &           ISSM,ICSM,ISIGDEN,NCBATCH,
     &           T,
     &           IDIAG,IRIS,IRIC,
     &           NSPOBEX_TP,ITSPOBEX_TP,IBSPOBEX_TP,
     &           ISPOBEX_AC,
     &           TOCC1,TOCC2,TOCC3,TOCC4,
     &           I1G,XI1G,I2G,XI2G,I3G,XI3G,I4G,XI4G,
     &           I1GE,XI1GE,I2GE,XI2GE,I3GE,XI3GE,I4GE,XI4GE,
     &           SIGN_NHX,ISPOBEX_FRZ)
*
*
* Jeppe Olsen, Summer of 1999
*              July 2000, Combinations of T-amplitudes and 
*                         determinants introduced
*              September 2000, ISPOBEX_FRZ added
*
* Timo Fleig, 2001/2002
*                  Modifications for relativistic CI
*                                 
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included 
*
* NSSOA : Number of strings per type and symmetry for alpha strings
* NSSOB : Number of strings per type and symmetry for beta strings
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
*
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "gasstr.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "cstate.inc"
#include "strbas.inc"
*. Specific input
      INTEGER ISBLOCK(8,*)
      INTEGER LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*.General input
      INTEGER NSSOA_S(NSMST,*),NSSOB_S(NSMST,*)
      INTEGER NSSOA_C(NSMST,*),NSSOB_C(NSMST,*)
      dimension SIGN_NHX(NSPOBEX_TP)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*)
      DIMENSION I1(*),I2(*),I3(*),I4(*)
      DIMENSION XI1S(*),XI2S(*),XI3S(*)
*
      DIMENSION CJRES(*),SIRES(*)
*. T-coefficients
      DIMENSION T(*)
      dimension ITSPOBEX_TP(4*NGAS,*),IBSPOBEX_TP(*)
      INTEGER ISPOBEX_AC(*),ISPOBEX_FRZ(*)
*
      CALL QENTER('SIGD4')
*
* Set undefined array to zero to avoid possible initialization 
* crashes. (JEPPE: What is this array supposed for?) FIXME.
* I assume no densities from frozen orbitals!
      IZERO = 0
      call isetvc(ISPOBEX_FRZ,IZERO,NSPOBEX_TP)
*
      NTESTL = 000
      NTEST = MAX(NTESTL,max(IPRS,IPRD))
*
      IF(NTEST.GE.800) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' SIGDEN4 speaking :'
        WRITE(6,*) ' ================'
        WRITE(6,*)  
        WRITE(6,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(6,*) ' TTSS and IOFF for each ACTIVE sigma block'
          DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0) 
     &      WRITE(6,'(10X,I6,4I3,I8)') 
     &      IBLOCK,(ISBLOCK(II,IBLOCK),II=1,5)
          END DO
        WRITE(6,*) ' List of active T blocks : '
        CALL IWRTMA(ISPOBEX_AC,1,NSPOBEX_TP,1,NSPOBEX_TP)
      END IF
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' Initial C vector '
        CALL WRTVCD(CB,LUC,1,-1)
      END IF
*
      IOCTPA = IBSPGPFTP(IASTP)
      IOCTPB = IBSPGPFTP(IBSTP)
      JOCTPA = IBSPGPFTP(IACTP)
      JOCTPB = IBSPGPFTP(IBCTP)
*
      NSTT_BLK = NSBLOCK/NSMST
* Loop over batches over C blocks      
      REWIND LUC
*
      DO 20000 JCBATCH = 1, NCBATCH             
*
*. Read C blocks into core
*
        NJBLOCK = LCBLOCK(JCBATCH)
        NCTT_BLK = NJBLOCK/NSMST
        I1C = I1CBLOCK(JCBATCH)
*
*. Loop over TT blocks of sigma and C in batches and 
*  obtain  contributions 
        DO 9000 ICTT_BLK = 1, NCTT_BLK
          IREAD = 0
*. first block of next TT block of C    
          ICBLK = I1C + (ICTT_BLK-1)*NSMST
          JATP = ICBLOCK(1,ICBLK)
          JBTP = ICBLOCK(2,ICBLK)
          ICOFF = ICBLOCK(5,ICBLK)
*
          IF(NTEST.GE.500) THEN
            WRITE(6,*) ' Next block of C, ICBLK,JATP,JBTP,ICOFF ',
     &      ICBLK,JATP,JBTP, ICOFF
          END IF
*
          IRATP = JATP + JOCTPA - 1
          IRBTP = JBTP + JOCTPB - 1
*
          DO 10000 ISTT_BLK = 1, NSTT_BLK
*. first block of next TT block of Sigma
            ISBLK = (ISTT_BLK-1)*NSMST + 1
*
            IATP = ISBLOCK(1,ISBLK)
            IBTP = ISBLOCK(2,ISBLK)
            ISOFF = ISBLOCK(5,ISBLK)
*
            ILATP = IATP + IOCTPA - 1
            ILBTP = IBTP + IOCTPB - 1
*. Connections ?
            DO 8000 ITTP = 1, NSPOBEX_TP
*
              ICA_OFF = 1
              ICB_OFF = 1 +  NGAS
              IAA_OFF = 1 +2*NGAS
              IAB_OFF = 1 +3*NGAS
*. Connections ?
               CALL GXFSTR(NELFSPGP(1,ILATP),NELFSPGP(1,ILBTP),
     &                     NELFSPGP(1,IRATP),NELFSPGP(1,IRBTP),
     &                     ITSPOBEX_TP(ICA_OFF,ITTP), 
     &                     ITSPOBEX_TP(IAA_OFF,ITTP),
     &                     ITSPOBEX_TP(ICB_OFF,ITTP), 
     &                     ITSPOBEX_TP(IAB_OFF,ITTP),
     &                     NGAS,ICON)
*
               IF(ICON.EQ.1.AND.ISPOBEX_AC(ITTP).EQ.1.AND.
     &            .NOT.(ISPOBEX_FRZ(ITTP).EQ.1.AND.ISIGDEN.EQ.2)) THEN
*. Read block in   
                 IF(IREAD.EQ.0) THEN
*. Number of blocks
                   NCBLK = NSMST
*. And fetch the c-blocks and transform into unpacked blocks in
*. determinant normalization
                   DO JBLK = ICBLK, ICBLK + NCBLK - 1           
                     JASM = ICBLOCK(3,JBLK)
                     JBSM = ICBLOCK(4,JBLK)
                     JOFF = ICBLOCK(5,JBLK)
                     ISCALE = 1
                     CALL GSTTBLD(CB(JOFF),JATP,JASM,JBTP,JBSM,
     &                    WORK(KNSTSO(IACTP)),WORK(KNSTSO2(IBCTP)), 
     &                    PSSIGN,IDC,PLSIGN,
     &                    LUC,C2,NSMST,ISCALE,XSCALE)
*
                     IF(NTEST.GE.500) THEN
                       WRITE(6,*) ' C block read in, JBLK, JOFF : ',
     &                 JBLK,JOFF
                       LEN = ICBLOCK(7,JBLK)
                       CALL WRTMAT(CB(JOFF),1,LEN,1,LEN)
                     END IF
*
                   END DO
                   IREAD = 1
                 END IF
*                ^ End if block should be read in 
                 ITB = IBSPOBEX_TP(ITTP)
*
* Setting undefined variables to zero:
                 ITABTRNSP = 0
                 ICRESTRICT = 0
                 ISRESTRICT = 0
                 ITDIAG = 0
                 SIGN_EXTERNAL = SIGN_NHX(ITTP)
                 CALL GNSIDE_REL(ISIGDEN,IRIS,IRIC,
     &                  ITSPOBEX_TP(ICA_OFF,ITTP),
     &                  ITSPOBEX_TP(ICB_OFF,ITTP),
     &                  ITSPOBEX_TP(IAA_OFF,ITTP),
     &                  ITSPOBEX_TP(IAB_OFF,ITTP),
     &                  NELFSPGP(1,ILATP),NELFSPGP(1,ILBTP),
     &                  NELFSPGP(1,IRATP),NELFSPGP(1,IRBTP),
     &                  NSSOA_S(1,IATP),NSSOB_S(1,IBTP),
     &                  NSSOA_C(1,JATP),NSSOB_C(1,JBTP),
     &                  T(ITB),ITTP,
     &                  SB(ISOFF),CB(ICOFF),ISSM,ICSM,
     &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  XINT,CJRES,SIRES,MAXK,
     &                  TOCC1,TOCC2,TOCC3,TOCC4,
     &                  I1G,XI1G,I2G,XI2G,I3G,XI3G,I4G,XI4G,
     &                  I1GE,XI1GE,I2GE,XI2GE,I3GE,XI3GE,I4GE,XI4GE,
     &                  ITDIAG,ITABTRNSP,ICRESTRICT,ISRESTRICT,
     &                  SIGN_EXTERNAL,IPRS,IPRD)
               END IF
*              ^ End if connection 
 8000       CONTINUE
*.          ^  End of loop over TT sigma blocks 
            IF(NTEST.GE.300) THEN
              WRITE(6,*) ' Updated sigma blocks ' 
              CALL WRTVH1(SB(ISOFF),ISSM,NSSOA_S(1,IATP),
     &                    NSSOB_S(1,IBTP),NSMST,0) 
            END IF
10000     CONTINUE
*         ^ End of loop over blocks of sigma 
          IF(IREAD.EQ.0) THEN
*. These C blocks were not needed, skip them
            DO JBLK = 1, NSMST
              CALL IFRMDS(LBL,-1,1,LUC)
              CALL SKPRCD2(LBL,-1,LUC)
            END DO
          END IF
*
 9000   CONTINUE
*.      ^ End of loop over TT C blocks in Batch
20000 CONTINUE
*.    ^End of loop over batches of C blocks
*
      CALL QEXIT('SIGD4')
      RETURN
      END
