      SUBROUTINE GSBBD1DG(RHO1,NACOB,ISCSM,ISCTP,
     &     ICCSM,ICCTP,IGRP,NROW,
     &     NGAS,ISEL,ICEL,
     &     SB,CB,
     &     NOBPTS,IOBPTS,MAXI,MAXK,
     &     SSCR,CSCR,I1,XI1S,I2,XI2S,
     &     NSMOB,RHO1S,
     &     SCALEF,SCLFAC,IPRNT)
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen
*
*  modified for relativistic double group calculation
*  Real/imaginary logic taken care of outside
*  cleaned and shaved ...
*  summer of '97 and april '98, Timo Fleig
*
* =====
* Input
* =====
* RHO1  : One body density matrix to be calculated/updated
* NACOB : Number of active orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active spaces
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* CB   : Input C block
* SB   : Input S block
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of spinors per type and symmetry (this set)
* IOBPTS : base for spinors of given type and symmetry (this set)
* NSMOB : Number of symmetries of spinors
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO1 : Calculated/updated density block
*        This is real or imaginary part of either the
*        barred-barred or the unbarred-unbarred partition
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
* I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
*              type and symmetry
* RHO1S : Scratch space for one electron density
*
* Jeppe Olsen, Winter of 1991
* Updated for GAS , August '95
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
*. General input
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
*
* Double group symmetry info
#include "symm.inc"
*
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO1(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),RHO1S(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*)
      DIMENSION I2(MAXK,*),XI2S(MAXK,*)
*
*.Local arrays
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS)
*
      NTESTL = 0000
      NTEST = max(NTESTL,IPRNT)
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' ==================='
         WRITE(6,*) ' GSBBD1DG in action '
         WRITE(6,*) ' ==================='
         WRITE(6,*)
         WRITE(6,*) ' Occupation of active left strings '
         CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
         WRITE(6,*) ' Occupation of active Right strings '
         CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
*
         WRITE(6,*) ' GSBBD1DG, sclfac ',SCLFAC
         WRITE(6,*) ' GSBBD1DG, scalef for r/i: ',SCALEF
      END IF
*
      IFRST = 1
      JFRST = 1
*
* Type of single excitations that connects the two column strings
      CALL SXTYP_GAS(NSXTP,ITP,JTP,NGAS,MXPNGAS,ISEL,ICEL)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = IDBGMULT(ISCSM,IADJSYM(ICCSM))
      IF(NTEST.GE.1000)
     &     WRITE(6,*) ' ISCSM,ICCSM,IJSM ', ISCSM,ICCSM,IJSM
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP = 1,NSXTP,1
         ITYP = ITP(IJTP)
         JTYP = JTP(IJTP)
         IF(NTEST.GE.1000) write(6,*) ' ITYP JTYP ', ITYP,JTYP
         DO 800 ISM = 1,NSMOB,1
*. new i and j so new intermediate strings
            KFRST = 1
*
            JSM = IDBGMULT(ISM,IADJSYM(IJSM))
            IF(JSM.EQ.0) GOTO 800
            IF(NTEST.GE.1000) write(6,*) ' ISM JSM ', ISM,JSM
            NIORB = NOBPTS(ITYP,ISM)
            NJORB = NOBPTS(JTYP,JSM)
            IBIORB = IOBPTS(ITYP,ISM)
            IBJORB = IOBPTS(JTYP,JSM)
            if (NTEST.ge.1000) then
               write(6,*) ' NIORB NJORB ', NIORB,NJORB
               write(6,*) ' offset left orb ', IBIORB
               write(6,*) ' offset right orb ', IBJORB
            end if
            IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*
COLD. Loop over partitionings of the row strings
            NIPART = NROW/MAXI
            IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
*. Loop over partitionings of N-1 strings (K strings)
            KBOT = 1-MAXK
            KTOP = 0
 700        CONTINUE
            KBOT = KBOT + MAXK
            KTOP = KTOP + MAXK
*. Single excitation information independent of I strings
*
* Determine mapping <Jb|j_+|Kb>
            J12 = 1
            K12 = 1
* What are these good for? (relevant for info on common block HIDSCR)
            IF(JFRST.EQ.1) KFRST=1
            CALL ADST_GAS(IBJORB,NJORB,ICCTP,ICCSM,IGRP,KBOT,KTOP,
     &           I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,J12,K12,
     &           KKJACT,SCLFAC)
            JFRST = 0
            KFRST = 0
*. Determine mapping <Ib|i_+|Kb>
            I12 = 2
            IF(IFRST.EQ.1) KFRST = 1
            CALL ADST_GAS(IBIORB,NIORB,ISCTP,ISCSM,IGRP,KBOT,KTOP,
     &           I2,XI2S,MAXK,NKBTC,KEND,IFRST,KFRST,I12,K12,
     &           KKIACT,SCLFAC)
            IFRST = 0
            KFRST = 0
*. Appropriate place to start partitioning over I strings
*. Loop over partitionings of the row strings
            DO 701 IPART = 1, NIPART
               IBOT = (IPART-1)*MAXI+1
               ITOP = MIN(IBOT+MAXI-1,NROW)
               NIBTC = ITOP - IBOT + 1
*
* Obtain CSCR(I,K,JORB) = SUM(J) <K!j_+!J> . C(I,J)
*.Gather C Block
               DO JJORB = 1,NJORB
                  ICGOFF = 1 + (JJORB-1)*NKBTC*NIBTC
                  CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                 NKBTC,I1(1,JJORB),XI1S(1,JJORB))
               END DO
*
* Obtain SSCR(I,K,IORB) = SUM(I) <I|i_+|K> . S(I,J)
*.Gather S Block
               DO IIORB = 1,NIORB
                  ISGOFF = 1 + (IIORB-1)*NKBTC*NIBTC
                  CALL MATCG(SB,SSCR(ISGOFF),NROW,NIBTC,IBOT,
     &                 NKBTC,I2(1,IIORB),XI2S(1,IIORB))
               END DO
               NKI = NKBTC*NIBTC
*
               FACTORC = 0.0D0
               FACTORAB = SCALEF
               CALL MATML7(RHO1S,SSCR,CSCR,NIORB,NJORB,NKI,NIORB,
     &              NKI,NJORB,FACTORC,FACTORAB,1)
* "1" means transpose matrix SSCR, so summation is over NKI,NKI as required!
*. Scatter out to complete matrix
               DO 610 JJORB = 1, NJORB
                  JORB = IBJORB-1+JJORB
                  DO 605 IIORB = 1, NIORB
                    IORB = IBIORB-1+IIORB
                    if (NTEST.ge.1000) then
                    write(6,'(A13,1X,F10.8,2X,A16,I5,A3,F10.8)')
     &                    'Adding RHO1S =',
     &                    RHO1S((JJORB-1)*NIORB+IIORB),
     &                    'to element RHO1(',
     &                    (JORB-1)*NACOB+IORB,')= ',
     &                    RHO1((JORB-1)*NACOB+IORB)
                    end if
                    RHO1((JORB-1)*NACOB+IORB) =
     &                   RHO1((JORB-1)*NACOB+IORB) +
     &                   RHO1S((JJORB-1)*NIORB+IIORB)
                    if (NTEST.ge.1000) then
                      write(6,'(3X,A3,F16.14)') ' = ',
     &                   RHO1((JORB-1)*NACOB+IORB)
                    end if
 605              CONTINUE
 610           CONTINUE

 701        CONTINUE
*. /\ end of this I partitioning
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
 800     CONTINUE
*.(end of loop over symmetries)
 900  CONTINUE
 1001 CONTINUE
*
C!    stop ' enforced stop in RSBBD1 '
      RETURN
      END
