C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : traone.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRAONE */
      SUBROUTINE TRAONE(WORK,KFREE,LFREE,NSTR,NSTRT,NSPC,KQ,KQC,KE,KIBE)
C***********************************************************************
C
C     Driver routine for 2-index transformations
C
C     Written by L. Visscher Mar 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxorb.h>
#include <dcbtra.h>
#include <dcbtr3.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbgen.h>
      LOGICAL TOBE,TRASAM
      DIMENSION WORK(*),NSTR(2,0:2,2),NSTRT(2),KQ(2,2),
     &          NSPC(2,0:2),KQC(2),KE(2,2),KIBE(2,4)
C
      CALL QENTER('TRAONE')
      KFRSAV = KFREE
C
C     Get hold of the slaves if we run in parallel.
C        ( ITASK = 1 for Fock matrices )
C
      IF (PARCAL) CALL DIRAC_PARCTL( 1 )
C
      CALL TRAHI(IPRTRA,2)
C
C     Transform core-fock matrix to active MO basis
C     ---------------------------------------------
C
      CALL MEMGET('REAL',KFMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL GTCFCK(WORK(KFMO),IPRTRA,WORK(KQ(1,1)),WORK(KQ(1,2)),
     &            WORK(KQC(1)),NDMOQR,ICMOQR,NSTR,NDMOQC,ICMOQC,NSPC,
     &            ITRA_INTFL2,WORK,KFREE,LFREE)
C
C     Write MOLFDIR file MRCONEE
C     --------------------------
C
      CALL MEMGET('REAL',KORBE,NSTRT(1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBORB,NSTRT(1),WORK,KFREE,LFREE)
C
      CALL MRCONE(WORK(KFMO),WORK(KORBE),WORK(KE(1,1)),WORK(KE(2,1)),
     &  WORK(KIBORB),WORK(KIBE(1,1)),WORK(KIBE(2,1)),
     &  NSPC,NSTR,NSTRT,IPRTRA,WORK,KFREE,LFREE)
C
C     Release the slaves if we run in parallel.
C        ( ITASK = -1 )
C
      IF (PARCAL) CALL DIRAC_PARCTL( -1 )
C
      CALL MEMREL('TRAONE',WORK,KFRSAV,KFMO,KFREE,LFREE)
      CALL QEXIT('TRAONE')
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gtcfck */
      SUBROUTINE GTCFCK(FMO,IPRINT,Q1,Q2,QCO,NDMOQR,ICMOQR,NSTR,
     &                  NDMOQC,ICMOQC,NSPC,INTFLG2,WORK,KFREE,LFREE)
C*****************************************************************************
C     
C     Get core Fock-matrix in MO-basis; 
C     use active coefficients for MO-transformation
C
C     Based on Trond Saue's GETFCK, Luuk Visscher, Mar 17 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbdhf.h>
#include <dgroup.h>
C
      LOGICAL TOBE
      DIMENSION FMO(*),WORK(*)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2),NSTR(2,0:2,2)
      DIMENSION NDMOQC(2,2),ICMOQC(2),NSPC(2,0:2)
      DIMENSION Q1(*),Q2(*),QCO(*)
C
      CALL QENTER('GTCFCK')
      KFRSAV = KFREE
C
C     Set INTFLG in dcbdhf
C
      INTFLG = INTFLG2
C
      CALL MEMGET('REAL',KFAO ,N2BBASXQ,WORK,KFREE,LFREE)      
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)      
C
      CALL GTCFC1(FMO,WORK(KFAO),WORK(KDMAT),Q1,Q2,QCO,NDMOQR,ICMOQR,
     &            NSTR,NDMOQC,ICMOQC,NSPC,IPRINT,WORK,KFREE,LFREE)
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
          CALL HEADER('GTCFCK:Total Fock matrix in MO-basis',-1)
          CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,
     &                NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      WRITE (LUPRI,1000) DHFERG,ELERGY,E1PART,E2PART
C
      CALL MEMREL('GTCFCK',WORK,KFRSAV,KFAO,KFREE,LFREE)
      CALL QEXIT('GTCFCK')
C
 1000 FORMAT (//' Core energy (includes nuclear repulsion) :',T50,F20.10
     & /' - Electronic part :',T50,F20.10
     & /'   - One-electron terms :',T50,F20.10
     & /'   - Two-electron terms :',T50,F20.10)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GTCFC1 */
      SUBROUTINE GTCFC1(FMO,FAO,DMAT,Q1,Q2,QCO,NDMOQR,ICMOQR,NSTR,
     &                  NDMOQC,ICMOQC,NSPC,IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C     
C     Get core Fock-matrix in MO-basis; 
C     use active coefficients for MO-transformation
C
C     Based on Trond Saue's GETFC1, Luuk Visscher, Mar 17 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <aovec.h>
#include <mxorb.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include <dcbgen.h>
#include <dcbfir.h>
#include <dcbdhf.h>
#include <blocks.h>
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbham.h>
#include <cbihr2.h>
      LOGICAL TOBE,BF1,BF2,BF3
      DIMENSION FMO(*),FAO(*),DMAT(*),WORK(*)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2),NSTR(2,0:2,2)
      DIMENSION NDMOQC(2,2),ICMOQC(2),NSPC(2,0:2)
      DIMENSION Q1(*),Q2(*),QCO(*)
C
      CALL QENTER('GTCFC1')
      KFRSAV = KFREE
      LUBUF = 22
C
C     Initialize
C
      NFMAT = 1
      NOPEN = 0
C
C     Check if we need two-electron part
C
      NTOT = 0
      DO IFRP = 1,NFSYM
        NTOT = NTOT + NSPC(IFRP,0)
      ENDDO
      IF(NTOT.GT.0) THEN
C
C       1. Get two-electron Fock matrix in AO-basis
C       ===========================================
C
        NPOS = 0
        IF(PARCAL) THEN
           BF1  = FIRST1
           BF2  = FIRST2
           BF3  = FIRST3
           FIRST1 = .TRUE.
           FIRST2 = .TRUE.
           FIRST3 = .TRUE.
           NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
           NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
           NPOS = NPOS + NSMLBL*(NSMLBL+1)/2
        ENDIF
        CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
C
C       Closed shell Hartree - Fock
C       ===========================
C       Totally symmetric operator
        ISYMOP(1) = 1
C       Fock matrix type
        IFCKOP(1) = 1
C       Hermitian operator
        IHRMOP(1) = 1
C       Get core density matrix
        CALL CDENS(DMAT,QCO,NDMOQC,ICMOQC,NSPC,IPRINT)        
        CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FMO,DMAT,NFMAT,
     &       WORK(KPOS),INTFLG,IPRINT,WORK(KFREE),LFREE)
        IF(PARCAL) THEN
          FIRST1 = BF1
          FIRST2 = BF2
          FIRST3 = BF3
        ENDIF
        CALL MEMREL('GTCFC1.TWOFCK',WORK,1,KFRSAV,KFREE,LFREE)
      ELSE
        CALL DZERO(FMO, N2BBASXQ)        
        CALL DZERO(DMAT,N2BBASXQ)        
      ENDIF
C
C     3. Get one-electron Fock matrix
C     ===============================
C
      INQUIRE(FILE='DFFCK1',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK1','OLD','GTCFC1')
        CALL REAFCK(LUBUF,FAO,.TRUE.,1)
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        WRITE(LUPRI,'(A)')
     &    '*** WARNING *** GTCFC1: No 1-Fock found. Regenerating.'
        CALL ONEFCK(FAO,IPRINT,WORK(KFREE),LFREE)
      ENDIF
C
C     4. Calculate core-core interaction.
C     ===================================
C
      CALL ERGCAL(FAO,FMO,DMAT)
C
C     5. Add 1- and 2-electron matrix and transform to MO-basis
C     =========================================================
C
      IF(NTOT.GT.0) CALL DAXPY(N2BBASXQ,D1,FMO,1,FAO,1) 
C
C     Make the matrix spinfree if desired.
C
C     Note that we did not make the matrices spinfree prior to the energy
C     calculation since this would require one extra large (fock-matrix size)
C     array. For closed shell cores we thus use that the first order spin-orbit
C     effect is zero. This assumption should be checked when implementing
C     non-closed shell cores or unrestricted DHF algorithms !!
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GTCFC1:Total Fock matrix in AO-basis',-1)
        CALL PRQMAT(FAO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF
C
      IF (SPINFR) CALL SPFAO (FMO,FAO,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.5.AND.SPINFR) THEN
        CALL HEADER('GTCFC1:Total spinfree Fock matrix in AO-basis',-1)
        CALL PRQMAT(FAO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF
C
      CALL DZERO(FMO,N2ORBXQ)
      DO 10 I = 1,NFSYM
        IF(NSTR(I,0,1).EQ.0.OR.NSTR(I,0,2).EQ.0) GOTO 10
        NRQ1 = NDMOQR(1,I,1)
        NCQ1 = NDMOQR(2,I,1)
        NRQ2 = NDMOQR(1,I,2)
        NCQ2 = NDMOQR(2,I,2)
        IF(IPRINT.GE.5) THEN
          CALL HEADER('GTCFC1:First transformation matrix',-1)
          CALL PRQMAT(Q1(ICMOQR(I,1)),NRQ1,NCQ1,NRQ1,NCQ1,NZ,
     &                IPQTOQ(1,0),LUPRI)
          CALL HEADER('GTCFC1:Second transformation matrix',-1)
          CALL PRQMAT(Q2(ICMOQR(I,2)),NRQ2,NCQ2,NRQ2,NCQ2,NZ,
     &                IPQTOQ(1,0),LUPRI)
        ENDIF
        CALL QTRANS('AOMO','S',D0,
     &            NFBAS(I,0),NFBAS(I,0),NSTR(I,0,1),NSTR(I,0,2),
     &            FAO(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            FMO(I2ORBX(I,I)+1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &            Q1(ICMOQR(I,1)),NRQ1,NCQ1,NZ,IPQTOQ(1,0),
     &            Q2(ICMOQR(I,2)),NRQ2,NCQ2,NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)
 10   CONTINUE
C
      CALL QEXIT('GTCFC1')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck CDENS */
      SUBROUTINE CDENS(DMAT,QCO,NDMOQC,ICMOQC,NSPC,IPRINT)
C*****************************************************************************
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION DMAT(*),QCO(*)
      DIMENSION NSPC(2,0:2),NDMOQC(2,2),ICMOQC(2)
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
C
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NSPC(I,0).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),
     &              QCO(ICMOQC(I)),NDMOQC(1,I),NDMOQC(2,I),
     &              1,NSPC(I,0),NFBAS(I,0))
   10 CONTINUE
C
C     Print section
C
      IF(IPRINT.GE.4) THEN
        CALL TITLER('Output from CDENS','*',103)
        DO 20 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 20
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(DMAT(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   20   CONTINUE
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mrcone */
      SUBROUTINE MRCONE(FMO,ORBE,ORBE1,ORBE2,IBORB,IBORB1,IBORB2,
     &                  NSPC,NSTR,NSTRT,IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C     
C     Write core Fock-matrix and other information in MRCONEE format 
C     A better interface needs to be defined, this will serve for the
C     moment and has the advantage that RELCCSD and DIRRCI need not be
C     modified.
C
C     Luuk Visscher, Mar 18 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0=0.D0)
C
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbdhf.h>
#include <dgroup.h>
C
      DIMENSION FMO(NORBT,NORBT,NZ)
      DIMENSION ORBE(*),IBORB(*),WORK(*)
      DIMENSION ORBE1(*),ORBE2(*),IBORB1(*),IBORB2(*)
      DIMENSION NSPC(2,0:2),NSTR(2,0:2,2),NSTRT(2)
C
      CALL QENTER('MRCONE')
      KFRSAV = KFREE
C
      IF (NSTRT(1).NE.NSTRT(2)) THEN
         WRITE (LUPRI,*) ' Ranges for bra and ket are different'
         WRITE (LUPRI,*) ' MRCONEE file is not written !!!'
         RETURN
      ENDIF
C
      CALL MEMGET('REAL',KFMOM,NSTRT(1)*NSTRT(1)*8,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIRPMO,NSTRT(1)*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIRPAMO,NSTRT(1)*2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KORBMO,NSTRT(1)*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBSPI,NSTRT(1)*2,WORK,KFREE,LFREE)
C
      II = 0
      DO I = 1, NFSYM
       DO J = 1, NSTR(I,0,1)
         II = II + 1
         IF (I.EQ.1) THEN
            ORBE(II) = ORBE1(J)
            IBORB(II) = IBORB1(J)
         ELSEIF (I.EQ.2) THEN
            ORBE(II) = ORBE2(J)
            IBORB(II) = IBORB2(J)
         ELSE
            CALL QUIT ('More than 2 fermion irreps not implemented') 
         ENDIF
       ENDDO
      ENDDO
C
      CALL MRCON1(FMO,ORBE,IBORB,WORK(KFMOM),NSTR,NSTRT,WORK(KIRPMO),
     &            WORK(KIRPAMO),WORK(KORBMO),WORK(KIBSPI),NSPC,
     &            IPRINT,WORK,KFREE,LFREE)
C
      CALL MEMREL('TRAONE',WORK,KFRSAV,KFMOM,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mrcon1 */
      SUBROUTINE MRCON1(FMO,ORBE,IBORB,FMOM,NSTR,NSTRT,IRPMO,IRPAMO,
     &                  ORBMO,IBSPI,NSPC,IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C     
C     Luuk Visscher, Mar 18 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbdhf.h>
#include <dcbham.h>
#include <dcbtra.h>
#include <dgroup.h>
C
      LOGICAL BREIT
      CHARACTER*14 REPN(8)
      CHARACTER*4 REPA(64)
      DIMENSION FMO(NORBT,NORBT,NZ),ORBE(NSTRT),IBORB(NSTRT)
      DIMENSION FMOM(NSTRT*2,NSTRT*2,2),WORK(*)
      DIMENSION IRPMO(NSTRT*2),IRPAMO(NSTRT*2)
      DIMENSION ORBMO(NSTRT*2),IBSPI(NSTRT*2)
      DIMENSION NSTR(2),IPHS(4),NAEL(8),NSPC(2,0:2)
      DIMENSION MULTB(64,64)
      DATA IPHS/ 1,-1,-1, 1/
C
C
C     Fixed information. At the moment we cannot treat the Breit 
C     interaction and the fermion representations are limited to 
C     gerade and ungerade.
C     We treat the quaternion case as pseudocomplex to avoid ordering
C     problems in the CCSD case.
C
      BREIT = .FALSE.
      IF (NFSYM.EQ.1.AND.NZ.EQ.4) THEN
         NSYMRP = 2
         REPN(1) = '  1A          '
         REPN(2) = '  2A          '
      ELSEIF (NFSYM.EQ.1.AND.NZ.EQ.2) THEN
         NSYMRP = 2
         REPN(1) = '  1E          '
         REPN(2) = '  2E          '
      ELSEIF (NFSYM.EQ.1.AND.NZ.EQ.1) THEN
         NSYMRP = 2
         REPN(1) = '   E  1E      '
         REPN(2) = '   E  2E      '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.4) THEN
         NSYMRP = 4
         REPN(1) = ' 1AG          '
         REPN(2) = ' 2AG          '
         REPN(3) = ' 1AU          '
         REPN(4) = ' 2AU          '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.2) THEN
         NSYMRP = 4
         REPN(1) = ' 1Eg          '
         REPN(2) = ' 2Eg          '
         REPN(3) = ' 1Eu          '
         REPN(4) = ' 2Eu          '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.1) THEN
         NSYMRP = 4
         REPN(1) = '  Eg 1Eg      '
         REPN(2) = '  Eg 2Eg      '
         REPN(3) = '  Eu 1Eu      '
         REPN(4) = '  Eu 2Eu      '
      ELSE
         PRINT*, 'Warning representation name is not defined'
      ENDIF
C
C     Get multiplication table for the abelian subgroup
C
      IF (SPINFR.OR.LEVYLE) THEN
         CALL GMULTSF (NSYMRPA,REPA,MULTB)
      ELSE
         CALL GMULTA (NZ,NFSYM,NSYMRPA,REPA,MULTB)
      ENDIF
C
C     Set the number of active electrons
C
      DO IREP = 1, NSYMRP
         IFRP = (IREP+1) / 2
         NAEL(IREP) = NOCC(IFRP) - NSPC(IFRP,1)
      ENDDO
C
C     The correlation codes do not use TR symmetry, write Kramers
C     pairs as spinors
C
      II    = 0
      IK    = 0
      IREP  = 0
      IREPA = 0
      DO I = 1, NFSYM
         IREP  = IREP  + 1
         IREPA = IREPA + 1
         DO ISTR = 1, NSTR(I)
            II = II + 1
            IK = IK + 1
            IRPMO(II)  = IREP
            IF (SPINFR.OR.LEVYLE) THEN
               IRPAMO(II) = IBORB(IK) + 1
            ELSE
               IRPAMO(II) = IREPA
            ENDIF
            ORBMO(II)  = ORBE(IK)
            IBSPI(II)  = IBORB(IK)
         ENDDO
         IREP = IREP + 1
         IF (NZ.NE.4) IREPA = IREPA + 1
         IK = IK - NSTR(I)
         DO ISTR = 1, NSTR(I)
            II = II + 1
            IK = IK + 1
            IRPMO(II)  = IREP
             IF (SPINFR.OR.LEVYLE) THEN
               IRPAMO(II) = IBORB(IK) + NBSYM + 1
            ELSE
               IRPAMO(II) = IREPA
            ENDIF
            ORBMO(II)  = ORBE(IK)
            IBSPI(II)  = IBORB(IK)
         ENDDO
      ENDDO
C
C     Fill in complete core-Fock matrix in complex form
C
      CALL DZERO(FMOM,NSTRT*NSTRT*8)
      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      DO IFSYM = 1, NFSYM
        IOFF3 = IOFF3 + NSTR(IFSYM)
        DO J = 1, NSTR(IFSYM)
          CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,1),1,
     &                           FMOM(1+IOFF2,J+IOFF2,1),1)
          CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,1),1,
     &                           FMOM(1+IOFF3,J+IOFF3,1),1)
        ENDDO
        IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN
          DO J = 1, NSTR(IFSYM)
            CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,2),1,
     &                             FMOM(1+IOFF2,J+IOFF2,2),1)
            CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,2),1,
     &                             FMOM(1+IOFF3,J+IOFF3,2),1)
            CALL DSCAL(NSTR(IFSYM),DM1,FMOM(1+IOFF3,J+IOFF3,2),1)
          ENDDO
          IF (NZ.GE.4.AND.(.NOT.SPINFR)) THEN
            DO J = 1, NSTR(IFSYM)
              CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,3),1,
     &                               FMOM(1+IOFF2,J+IOFF3,1),1)
              CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,3),1,
     &                               FMOM(1+IOFF3,J+IOFF2,1),1)
              CALL DSCAL(NSTR(IFSYM),DM1,FMOM(1+IOFF3,J+IOFF2,1),1)
              CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,4),1,
     &                               FMOM(1+IOFF2,J+IOFF3,2),1)
              CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,4),1,
     &                               FMOM(1+IOFF3,J+IOFF2,2),1)
            ENDDO
          ENDIF
        ENDIF
        IOFF1 = IOFF1 + NORB(IFSYM)
        IOFF2 = IOFF2 + 2*NSTR(IFSYM)
        IOFF3 = IOFF3 + NSTR(IFSYM)
      ENDDO
C
      IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Core Fock matrix ***'
          CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZ,2),
     &                IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Write MRCONEE file
C
      OPEN (LUMLF1,FILE='MRCONEE',FORM='UNFORMATTED')
      WRITE (LUMLF1) 2*NSTRT,BREIT,DHFERG
      WRITE (LUMLF1) NSYMRP,(REPN(IRP),IRP=1,NSYMRP),
     &              (NAEL(IRP),IRP=1,NSYMRP)
      WRITE (LUMLF1) NSYMRPA,(REPA(IRPA),IRPA=1,NSYMRPA*2)
      WRITE (LUMLF1) ((MULTB(I,J),I=1,2*NSYMRPA),J=1,2*NSYMRPA)
      WRITE (LUMLF1) (IRPMO(I),IRPAMO(I),ORBMO(I),I=1,NSTRT*2),
     &               (IBSPI(I),I=1,2*NSTRT),(NORB(I),I=1,NFSYM)
      WRITE (LUMLF1) ((FMOM(I,J,1),FMOM(I,J,2),
     &                  J=1,NSTRT*2),I=1,NSTRT*2)
      CLOSE (LUMLF1)
C
      WRITE (LUPRI,1000) 
C
      CALL QEXIT('MRCONE')
C
 1000 FORMAT (//' MOLFDIR file MRCONEE is written')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmulta */
      SUBROUTINE GMULTA (NZ,NFSYM,NSYMRPA,REPA,MULTB)
C*****************************************************************************
C     
C     Hard-wired multiplication table for the abelian subgroups of D2h
C     Defines representation names as well.
C
C     Luuk Visscher, Mar 18 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      CHARACTER*4 REPA(64)
      DIMENSION MULTB(64,64)
C
      IF (NFSYM.EQ.1.AND.NZ.EQ.4) THEN
         NSYMRPA = 1
         REPA(1) = '   A'
         REPA(2) = '   a'
         MULTB(1,1) = 2
         MULTB(2,1) = 1
         MULTB(2,2) = 2
      ELSEIF (NFSYM.EQ.1.AND.NZ.LE.2) THEN
         NSYMRPA = 2
         REPA(1) = '  1E'
         REPA(2) = '  2E'
         REPA(3) = '   a'
         REPA(4) = '   b'
         MULTB(1,1) = 4
         MULTB(2,1) = 3
         MULTB(3,1) = 1
         MULTB(4,1) = 2
         MULTB(2,2) = 4
         MULTB(3,2) = 2
         MULTB(4,2) = 1
         MULTB(3,3) = 3
         MULTB(4,3) = 4
         MULTB(4,4) = 3
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.4) THEN
         NSYMRPA = 2
         REPA(1) = '  AG'
         REPA(2) = '  AU'
         REPA(3) = '  ag'
         REPA(4) = '  au'
         MULTB(1,1) = 3
         MULTB(2,1) = 4
         MULTB(3,1) = 1
         MULTB(4,1) = 2
         MULTB(2,2) = 3
         MULTB(3,2) = 2
         MULTB(4,2) = 1
         MULTB(3,3) = 3
         MULTB(4,3) = 4
         MULTB(4,4) = 3
      ELSEIF (NFSYM.EQ.2.AND.NZ.LE.2) THEN
         NSYMRPA = 4
         REPA(1) = ' 1Eg'
         REPA(2) = ' 2Eg'
         REPA(3) = ' 1Eu'
         REPA(4) = ' 2Eu'
         REPA(5) = '  ag'
         REPA(6) = '  bg'
         REPA(7) = '  au'
         REPA(8) = '  bu'
C        The gg-block
         MULTB(1,1) = 6
         MULTB(2,1) = 5
         MULTB(5,1) = 1
         MULTB(6,1) = 2
         MULTB(2,2) = 6
         MULTB(5,2) = 2
         MULTB(6,2) = 1
         MULTB(5,5) = 5
         MULTB(6,5) = 6
         MULTB(6,6) = 5
C        The uu-block
         MULTB(3,3) = 6
         MULTB(4,3) = 5
         MULTB(7,3) = 1
         MULTB(8,3) = 2
         MULTB(4,4) = 6
         MULTB(7,4) = 2
         MULTB(8,4) = 1
         MULTB(7,7) = 5
         MULTB(8,7) = 6
         MULTB(8,8) = 5
C        The ug-block
         MULTB(3,1) = 8
         MULTB(4,1) = 7
         MULTB(7,1) = 3
         MULTB(8,1) = 4
         MULTB(3,2) = 7
         MULTB(4,2) = 8
         MULTB(7,2) = 4
         MULTB(8,2) = 3
         MULTB(7,5) = 7
         MULTB(8,5) = 8
         MULTB(7,6) = 8
         MULTB(8,6) = 7
C        The gu-block
         MULTB(5,3) = 3
         MULTB(6,3) = 4
         MULTB(5,4) = 4
         MULTB(6,4) = 3
      ELSE
         CALL QUIT('Multiplication table is not defined')
      ENDIF
C
      DO I = 1, 2*NSYMRPA
         DO J = 1, I-1
            MULTB(J,I) = MULTB(I,J)
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmultsf */
      SUBROUTINE GMULTSF (NSYMRPA,REPA,MULTB)
C*****************************************************************************
C
C     We make the multiplication table as the direct product of the
C     spatial and spin multiplication tables. For the latter we have
C     two irreps a(lpha) and b(eta) that can couple to give mS=0 or
C     mS = +/- 1. The last class is treated as giving one-dimensional
C     irreps, i.e. we do not employ the full spin symmetry group SU(2)
C     but take the subgroup C4 for this to get an Abelian group that
C     can distinguish the different spin-projection up to |mS|=3/2.
C     An easy extension is to take the subgroup C8 so that one can go
C     up to |mS|=7/2, but this increases the amount of irreps by 2 and
C     has to be taken into account in the multiplication tables.
C
C     Luuk Visscher, Mar 18 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
#include <pgroup.h>
C
      CHARACTER*4 REPA(64)
      CHARACTER*1 REPAS(8)
      DIMENSION MULTB(64,64),MULTBS(8,8)
C
#include <dcbibt.h>
C
C     Start by making the spin-multiplication table
C
      REPAS(1) = 'a' ! mS = + 1/2
      REPAS(2) = 'b' ! mS = - 1/2
      REPAS(3) = '3' ! mS = + 3/2
      REPAS(4) = '3' ! mS = - 3/2
      REPAS(5) = '0' ! mS = 0
      REPAS(6) = '4' ! mS = +/- 4
      REPAS(7) = '2' ! mS = + 2
      REPAS(8) = '2' ! mS = - 2
      MULTBS(1,1) = 7
      MULTBS(2,1) = 5
      MULTBS(3,1) = 8
      MULTBS(4,1) = 6
      MULTBS(5,1) = 1
      MULTBS(6,1) = 3
      MULTBS(7,1) = 4
      MULTBS(8,1) = 2
      MULTBS(2,2) = 8
      MULTBS(3,2) = 6
      MULTBS(4,2) = 7
      MULTBS(5,2) = 2
      MULTBS(6,2) = 4
      MULTBS(7,2) = 1
      MULTBS(8,2) = 3
      MULTBS(3,3) = 7
      MULTBS(4,3) = 5
      MULTBS(5,3) = 3
      MULTBS(6,3) = 1
      MULTBS(7,3) = 2
      MULTBS(8,3) = 4
      MULTBS(4,4) = 8
      MULTBS(5,4) = 4
      MULTBS(6,4) = 2
      MULTBS(7,4) = 3
      MULTBS(8,4) = 1
      MULTBS(5,5) = 5
      MULTBS(6,5) = 6
      MULTBS(7,5) = 7
      MULTBS(8,5) = 8
      MULTBS(6,6) = 5
      MULTBS(7,6) = 8
      MULTBS(8,6) = 7
      MULTBS(7,7) = 6
      MULTBS(8,7) = 5
      MULTBS(8,8) = 6
C
      DO I = 1, 8
         DO J = 1, I-1
            MULTBS(J,I) = MULTBS(I,J)
         ENDDO
      ENDDO
C
C     The multiplication table of the boson irreps can be
C     obtained by bit operations.
C     Make the direct product table
C 
      NSYMRPA = NBSYM * 4
C
      JSR = 0
      DO JS = 1, 8
         DO JREP = 0, NBSYM-1
            JSR = JSR + 1
            ISR = 0
            DO IS = 1, 8
               DO IREP = 0, NBSYM-1
                  ISR = ISR + 1
                  IJS = MULTBS(IS,JS)
                  IJREP = IBTXOR(IREP,JREP)
                  IJSR = (IJS-1)*NBSYM + IJREP + 1
                  MULTB(ISR,JSR) = IJSR
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     Give the representations a name.
C
      ISR = 0
      DO IS = 1, 8
         DO IREP = 0, NBSYM-1
            ISR = ISR + 1
            REPA(ISR) = REP(IREP)//REPAS(IS)
         ENDDO
      ENDDO
C
      RETURN
      END
