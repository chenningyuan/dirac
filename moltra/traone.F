C 
!
!      Copyright (c) 2008 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
C 
C
C
C VERSION : $Revision: 7732 $
C DATE    : $Date: 2008-09-16 22:30:19 +0200 (Tue, 16 Sep 2008) $
C FILE    : traone.F
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
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gtcfck */
      SUBROUTINE GTCFCK(FMO,IPRINT,Q1,Q2,QCO,NDMOQR,ICMOQR,NSTR,
     &                  NDMOQC,ICMOQC,NSPC,INTFLG2,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Get core Fock-matrix in MO-basis; 
C     use active coefficients for MO-transformation
C
C     Called from traone.F/TRAONE
C                 trapam.F/RESOLV
C
C     Based on Trond Saue's GETFCK, Luuk Visscher, Mar 17 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
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
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "aovec.h"
#include "maxorb.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dcbfir.h"
#include "dcbdhf.h"
#include "blocks.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "cbihr2.h"
      LOGICAL TOBE, SAVEFLAGS(4)
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
      NOPENBUF = NOPEN
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
        call SaveTaskDistribFlags(saveflags)
       call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. , .TRUE. /))
        call SetIntTaskArrayDimension(NPOS,PARCAL)
        if (NPOS.GT.0) THEN
           CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
        else
           KPOS = KFREE
        endif
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

        IF(PARCAL) call SetTaskDistribFlags(saveflags)

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
      CALL ERGCAL(FAO,FMO,DMAT,WORK(KFREE),LFREE)
C
C     IF LVCON is used we should take the repulsion between all charge while
C     for the core-core interaction we leave out the valence part. Add on the
C     valence-valence SS-repulsion by recalculating CORRLV
C
C     LV+OF 20-8-2002 : To be tested !
C 
C     IF (LVCON) THEN
C        CORRLV_CORE = CORRLV
C        CALL GENDEN(DMAT,CMO,1,IPRSCF)
C        CALL LVCORR(DMAT,WORK,LWORK,IPRSCF)
C        DHFERG = DHFERG + CORRLV - CORRLV_CORE
C     ENDIF
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
CMI   ... in the BSS mode the spin-free AO matrix is ready on ONEFCK
      IF (SPINFR.AND.(.NOT.BSS)) CALL SPFAO (FMO,FAO,WORK(KFREE),LFREE)
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
C     Restore NOPEN
C
      NOPEN = NOPENBUF
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
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION DMAT(*),QCO(*)
      DIMENSION NSPC(2,0:2),NDMOQC(2,2),ICMOQC(2)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NSPC(I,0).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),D1,D0,
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mrcone */
      SUBROUTINE MRCONE(FMO,ORBE,ORBE1,ORBE2,IBORB,IBORB1,IBORB2,
     &                  NSPC,NSTR,NSTRT,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Write core Fock-matrix and other information in MRCONEE format 
C     A better interface needs to be defined, this will serve for the
C     moment and has the advantage that RELCCSD and DIRRCI need not be
C     modified.
C
C     Luuk Visscher, Mar 18 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
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
      CALL QEXIT('MRCONE')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mrcon1 */
      SUBROUTINE MRCON1(FMO,ORBE,IBORB,FMOM,NSTR,NSTRT,IRPMO,IRPAMO,
     &                  ORBMO,IBSPI,NSPC,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Luuk Visscher, Mar 18 1997
C
C     Added number of boson symmetry reps (of parent group).
C      Timo Fleig, Mar 28, 2002
C
C     Miro Ilias, Febr.2007   
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbtra.h"
#include "dgroup.h"
C
C     NORDER gives the order of the Abelian subgroup used linear molecules. 
C    NORDER=16 (default) uses C_16_h for D_inf_h and C_32 for C_inf_v.
C
      PARAMETER (NORDER=16)
      LOGICAL BREIT,FORWARD
      CHARACTER*14 REPN(8)
      CHARACTER*4 REPA(4*NORDER)
      DIMENSION FMO(NORBT,NORBT,NZ),ORBE(NSTRT),IBORB(NSTRT)
      DIMENSION FMOM(NSTRT*2,NSTRT*2,2),WORK(*)
      DIMENSION IRPMO(NSTRT*2),IRPAMO(NSTRT*2)
      DIMENSION ORBMO(NSTRT*2),IBSPI(NSTRT*2)
      DIMENSION NSTR(2),IPHS(4),NAEL(8),NSPC(2,0:2)
      DIMENSION MULTB(64,64), MJ_TO_REP(-2*NORDER:2*NORDER,2)
      DATA IPHS/ 1,-1,-1, 1/

      CALL QENTER('MRCON1')

      KFRSAV = KFREE
C
C
C     Fixed information. At the moment we cannot treat the Breit 
C     interaction and the fermion representations are limited to 
C     gerade and ungerade.
C     We treat the quaternion case as pseudocomplex to avoid ordering
C     problems in the CCSD case.
C
      BREIT = .FALSE.
      FORWARD = .TRUE.
      IF (NFSYM.EQ.1.AND.NZ.EQ.4) THEN
         NSYMRP = 2
         NBSYMRP = 1
         REPN(1) = '  1A          '
         REPN(2) = '  2A          '
      ELSEIF (NFSYM.EQ.1.AND.NZ.EQ.2) THEN
         NSYMRP = 2
         NBSYMRP = 2
         REPN(1) = '  1E          '
         REPN(2) = '  2E          '
      ELSEIF (NFSYM.EQ.1.AND.NZ.EQ.1) THEN
         NSYMRP = 2
         NBSYMRP = 4
         REPN(1) = '   E  1E      '
         REPN(2) = '   E  2E      '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.4) THEN
         NSYMRP = 4
         NBSYMRP = 2
         REPN(1) = ' 1AG          '
         REPN(2) = ' 2AG          '
         REPN(3) = ' 1AU          '
         REPN(4) = ' 2AU          '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.2) THEN
         NSYMRP = 4
         NBSYMRP = 4
         REPN(1) = ' 1Eg          '
         REPN(2) = ' 2Eg          '
         REPN(3) = ' 1Eu          '
         REPN(4) = ' 2Eu          '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.1) THEN
         NSYMRP = 4
         NBSYMRP = 8
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
      ELSEIF (LINEAR) THEN
         CALL GMULTLN (NSYMRPA,REPA,MJ_TO_REP,MULTB,NORDER)
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
C     The correlation codes (RELCCSD, DIRRCI, LUCITA) do not use 
C     TR symmetry, write Kramers pairs as spinors
C
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
            ELSEIF (LINEAR) THEN
               MJ = IBORB(IK)
               IF (MJ.GT. 2*NORDER/NFSYM) MJ = MJ - 4*NORDER/NFSYM
               IF (MJ.LT.-2*NORDER/NFSYM) MJ = MJ + 4*NORDER/NFSYM
               IRPAMO(II) = MJ_TO_REP(MJ,I)
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
            ELSEIF (LINEAR) THEN
               MJ = -IBORB(IK)
               IF (MJ.GT. 2*NORDER/NFSYM) MJ = MJ - 4*NORDER/NFSYM
               IF (MJ.LT.-2*NORDER/NFSYM) MJ = MJ + 4*NORDER/NFSYM
               IRPAMO(II) = MJ_TO_REP(MJ,I)
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
      CALL MOLTODIR(FMO,FMOM,NSTR,NSTRT,FORWARD,IPRTRA)

C
C     Write MRCONEE file
C
C     Pass on all the information that the post-dhf
C     electron correlation codes might need.
C
      OPEN (LUMLF1,FILE='MRCONEE',FORM='UNFORMATTED')
C     Total number of spinors, breit active in DHF ?,
C     DHF energy, inversion symmetry (yes : 2; no : 1),
C     group type (1 real, 2 complex, 4 quat.),spinfree formalism,
C     total number of orbitals (so including frozen or deleted orbitals)
      WRITE (LUMLF1) 2*NSTRT,BREIT,DHFERG,NFSYM,NZ,SPINFR.OR.LEVYLE,
     &               NORBT
C     NSYMRP        - Number of fermion irreps in parent group
C     REPN          - names of these irreps (gerade, ungerade)
C     NAEL          - number of spinors active in the transformation
C     NSPC          - number of occupied frozen (core) spinors,
C                     0: total
C                     1: positive energy
C                     2: negative energy
C     NESH-NST-NSPC -number of deleted spinors.
      WRITE (LUMLF1) NSYMRP,(REPN(IRP),IRP=1,NSYMRP),
     &                      (NAEL(IRP),IRP=1,NSYMRP),
     &                      (NSTR(IRP),IRP=1,NFSYM),
     &                      (NSPC(IRP,0),IRP=1,NFSYM),
     &                      (NSPC(IRP,1),IRP=1,NFSYM),
     &                      (NSPC(IRP,2),IRP=1,NFSYM),
     &  (NESH(IRP)-NSTR(IRP)-NSPC(IRP,0),IRP=1,NFSYM)
C     Number of fermion irreps in the Abelian subgroup,
C     names of these irreps.
      WRITE (LUMLF1) NSYMRPA,(REPA(IRPA),IRPA=1,NSYMRPA*2)
C     Multiplication table for direct products in the
C     Abelian subgroup.
      WRITE (LUMLF1) ((MULTB(I,J),I=1,2*NSYMRPA),J=1,2*NSYMRPA)
C     Information for each spinor :
C     - Irrep in parent group (1:gerade, 2:ungerade)
C     - Irrep in Abelian subgroup
C     - Orbital energy taken from DHF
C     - Approximate boson irrep identification (needed in LUCIAREL
C       and LUCITA)
C     - Number of boson symmetry reps (for LUCITA)
      WRITE (LUMLF1) (IRPMO(I),IRPAMO(I),ORBMO(I),I=1,NSTRT*2),
     &               (IBSPI(I),I=1,2*NSTRT),(NORB(I),I=1,NFSYM),NBSYMRP
C     Matrix elements over the core Fock operator.
      WRITE (LUMLF1) ((FMOM(I,J,1),FMOM(I,J,2),
     &                 J=1,NSTRT*2),I=1,NSTRT*2)
      CLOSE (LUMLF1)
C
      WRITE (LUPRI,1000) 
C
      CALL QEXIT('MRCON1')
C
 1000 FORMAT (//' MOLFDIR file MRCONEE is written')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmulta */
      SUBROUTINE GMULTA (NZ,NFSYM,NSYMRPA,REPA,MULTB)
C***********************************************************************
C     
C     Hard-wired multiplication table for the abelian subgroups of D2h
C     Defines representation names as well.
C
C     Luuk Visscher, Mar 18 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
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
         CALL QUIT('GMULTA: Multiplication table is not defined')
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmultln */
      SUBROUTINE GMULTLN (NSYMRPA,REPA,MJ_TO_REP,MULTB,NORDER)
C***********************************************************************
C
C     We make the multiplication table for linear molecules where we 
C     have used the mj quantum number to identify blocks. Here the 
C     full point group is the non-Abelian Dinfh or Cinfv and we use 
C     the Abelian Cinfh or Cinf, respectively. Since the Cinfh and 
C     Cinv have an infinite number of irreps and we use a finite 
C     multiplication table to store the coupling information we need 
C     to truncate at some point. We do this at mj=15/2 in Cinfh and 
C     mj=31/2 in Cinf so that we can keep MULTB dimensioned at 64 like 
C     before. In fact this means that we'll be using the subgroups C8h 
C     and C16 respectively. This can be easily extended or decreased 
C     using the input parameter NORDER that indicates how large MULTB 
C     may become.
C
C     Luuk Visscher, Nov 03 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
C
      CHARACTER*4 REPA(4*NORDER)
      CHARACTER PARITY
      DIMENSION MULTB(64,64), MJ_TO_REP(-2*NORDER:2*NORDER,2)
C
#include "dcbibt.h"
C
      IF (4*NORDER.GT.64) CALL QUIT ("increase dimension of multb")
      NSYMRPA = 2*NORDER
C
      IF (NFSYM.EQ.1) THEN
C
C        Start by specifying the names of the irreps
C        and making the array MJ_TO_REP.
C
         IREP = 0
C
C        Fermion (half-integer) irreps
C
         DO MJ = 1, 2*NORDER-1, 2
            IREP = IREP+1
            WRITE (REPA(IREP),1000) MJ
            MJ_TO_REP(MJ,1) = IREP
            IREP = IREP+1
            WRITE (REPA(IREP),1000) -MJ
            MJ_TO_REP(-MJ,1) = IREP
         ENDDO
C
C        Boson (integer) irreps
C
         MJ = 0
         IREP = IREP + 1
         WRITE (REPA(IREP),1000) MJ
         MJ_TO_REP(MJ,1) = IREP
         DO MJ = 2, 2*NORDER-2, 2
            IREP = IREP+1
            WRITE (REPA(IREP),1000) MJ
            MJ_TO_REP(MJ,1) = IREP
            IREP = IREP+1
            WRITE (REPA(IREP),1000) -MJ
            MJ_TO_REP(-MJ,1) = IREP
         ENDDO
         MJ = 2*NORDER
         IREP = IREP + 1
         WRITE (REPA(IREP),1000) MJ
         MJ_TO_REP( MJ,1) = IREP
         MJ_TO_REP(-MJ,1) = IREP
C
C        The multiplication table is straigthforward since 
C        the product irrep is given by the sum of Mj values.
C        The only complication is the finite range of irreps so
C        we may need to add or subtract something, e.g. :
C        MJ = 33/2 -> MJ = -31/2 in the C16 subgroup.
C        (formulated in another way MJ = +/- 16 transforms
C         according to  the totally symmetric irrep in C16)
C        
          DO MJ = -2*NORDER, 2*NORDER, 1
             DO MI = -2*NORDER, 2*NORDER, 1
                MIJ = MI + MJ 
C               Map Cinf to C16
                IF (MIJ.LT.-2*NORDER) MIJ = MIJ + 4*NORDER
                IF (MIJ.GT. 2*NORDER) MIJ = MIJ - 4*NORDER
                IF (MIJ.EQ.-2*NORDER) MIJ = 2*NORDER
                IREP = MJ_TO_REP(MI,1)
                JREP = MJ_TO_REP(MJ,1)
                IJREP = MJ_TO_REP(MIJ,1)
                MULTB(IREP,JREP) = IJREP
             ENDDO
          ENDDO
C
      ELSE
C
C        Same procedure but we now also have inversion symmetry
C
C        Fermion (half-integer) irreps
C
         IREP = 0
         DO IFSYM = 1, NFSYM
            IF (IFSYM.EQ.1) THEN
               PARITY = 'g'
            ELSE
               PARITY = 'u'
            ENDIF
            DO MJ = 1, NORDER-1, 2
               IREP = IREP+1
               WRITE (REPA(IREP),1001) MJ,PARITY
               MJ_TO_REP(MJ,IFSYM) = IREP
               IREP = IREP+1
               WRITE (REPA(IREP),1001) -MJ,PARITY
               MJ_TO_REP(-MJ,IFSYM) = IREP
            ENDDO
         ENDDO
C
C        Boson (integer) irreps
C
         DO IFSYM = 1, NFSYM
            IF (IFSYM.EQ.1) THEN
               PARITY = 'g'
            ELSE
               PARITY = 'u'
            ENDIF
            MJ = 0
            IREP = IREP + 1
            WRITE (REPA(IREP),1001) MJ,PARITY
            MJ_TO_REP(MJ,IFSYM) = IREP
            DO MJ = 2, NORDER-2, 2
               IREP = IREP+1
               WRITE (REPA(IREP),1001) MJ,PARITY
               MJ_TO_REP(MJ,IFSYM) = IREP
               IREP = IREP+1
               WRITE (REPA(IREP),1001) -MJ,PARITY
               IF (MJ.NE.0) MJ_TO_REP(-MJ,IFSYM) = IREP
            ENDDO
C           For MJ = 16 we will not distinguish between positive and 
C           negative MJ
            MJ = NORDER
            IREP = IREP + 1
            WRITE (REPA(IREP),1001) MJ,PARITY
            MJ_TO_REP( MJ,IFSYM) = IREP
            MJ_TO_REP(-MJ,IFSYM) = IREP
         ENDDO
C
C        The multiplication table is rather straigthforward since 
C        the product irrep is the one with the sum of Mj values.
C        The only complication is the finite range of irreps so
C        we may need to add or subtract something, e.g. :
C        MJ = 17/2 -> MJ = -15/2 in the C8h subgroup.
C        (formulated in another way MJ = +/- 8 transforms
C         according to Ag or Au in C8h)
C        
          DO MJ = -NORDER, NORDER, 1
             DO MI = -NORDER, NORDER, 1
                MIJ = MI + MJ 
C               Map Cinfh to C8h
                IF (MIJ.LT.-NORDER) MIJ = MIJ + 2*NORDER
                IF (MIJ.GT. NORDER) MIJ = MIJ - 2*NORDER
                IF (MIJ.EQ.-NORDER) MIJ = NORDER
C               The gg-part
                IREP = MJ_TO_REP(MI,1)
                JREP = MJ_TO_REP(MJ,1)
                IJREP = MJ_TO_REP(MIJ,1)
                MULTB(IREP,JREP) = IJREP
C               The ug-part
                IREP = MJ_TO_REP(MI,2)
                JREP = MJ_TO_REP(MJ,1)
                IJREP = MJ_TO_REP(MIJ,2)
                MULTB(IREP,JREP) = IJREP
C               The gu-part
                IREP = MJ_TO_REP(MI,1)
                JREP = MJ_TO_REP(MJ,2)
                IJREP = MJ_TO_REP(MIJ,2)
                MULTB(IREP,JREP) = IJREP
C               The uu-part
                IREP = MJ_TO_REP(MI,2)
                JREP = MJ_TO_REP(MJ,2)
                IJREP = MJ_TO_REP(MIJ,1)
                MULTB(IREP,JREP) = IJREP
             ENDDO
          ENDDO
      ENDIF
C
 1000 FORMAT (I4)
 1001 FORMAT (I3,A)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmultsf */
      SUBROUTINE GMULTSF (NSYMRPA,REPA,MULTB)
C***********************************************************************
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "pgroup.h"
C
      CHARACTER*4 REPA(64)
      CHARACTER*1 REPAS(8)
      DIMENSION MULTB(64,64),MULTBS(8,8)
C
#include "dcbibt.h"
C
C     Start by making the spin-multiplication table
C
      REPAS(1) = 'a' ! mS = + 1/2
      REPAS(2) = 'b' ! mS = - 1/2
      REPAS(3) = '3' ! mS = + 3/2
      REPAS(4) = '3' ! mS = - 3/2
      REPAS(5) = '0' ! mS = 0
      REPAS(6) = '4' ! mS = +/- 4/2 = +/- 2
      REPAS(7) = '2' ! mS = + 2/2 = + 1
      REPAS(8) = '2' ! mS = - 2/2 = - 1
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck moltodir */
      SUBROUTINE MOLTODIR(FMO,FMOM,NSTR,NSTRT,FORWARD,IPRINT)
C***********************************************************************
C     
C Purpose: Transform DIRAC Fock matrix, FMO, into MOLFDIR format, FMOM
C          (if FORWARD=.true, otherwise transform back)     
C     
C     Luuk Visscher, Mar 26 2002
C         Miro Ilias, June 2007 - added IPRINT
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbtr3.h"
C
      LOGICAL FORWARD
      DIMENSION FMO(NORBT,NORBT,NZ)
      DIMENSION FMOM(NSTRT*2,NSTRT*2,2)
      DIMENSION NSTR(2)

      CALL QENTER('MOLTODIR')
C
C
C     Transform from matrix in quaternion DIRAC form to complex MOLFDIR form 
C     ... or backwards if we want to go from RELCCSD back to DIRAC
C
      IF (FORWARD) THEN
         CALL DZERO(FMOM,NSTRT*NSTRT*8)
      ELSE
         CALL DZERO(FMO ,NORBT*NORBT*NZ)
      ENDIF
C
      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      DO IFSYM = 1, NFSYM
        IOFF3 = IOFF3 + NSTR(IFSYM)
        DO J = 1, NSTR(IFSYM)
          IF (FORWARD) THEN
             CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,1),1,
     &                              FMOM(1+IOFF2,J+IOFF2,1),1)
             CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,1),1,
     &                              FMOM(1+IOFF3,J+IOFF3,1),1)
          ELSE
             CALL DCOPY(NSTR(IFSYM),FMOM(1+IOFF2,J+IOFF2,1),1,
     &                              FMO (1+IOFF1,J+IOFF1,1),1)
          ENDIF
        ENDDO
        IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN
          DO J = 1, NSTR(IFSYM)
            IF (FORWARD) THEN
               CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,2),1,
     &                                FMOM(1+IOFF2,J+IOFF2,2),1)
               CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,2),1,
     &                                FMOM(1+IOFF3,J+IOFF3,2),1)
               CALL DSCAL(NSTR(IFSYM),DM1,FMOM(1+IOFF3,J+IOFF3,2),1)
             ELSE
               CALL DCOPY(NSTR(IFSYM),FMOM(1+IOFF2,J+IOFF2,2),1,
     &                                FMO (1+IOFF1,J+IOFF1,2),1)
             ENDIF
          ENDDO
          IF (NZ.GE.4.AND.(.NOT.SPINFR)) THEN
            DO J = 1, NSTR(IFSYM)
              IF (FORWARD) THEN
                CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,3),1,
     &                                 FMOM(1+IOFF2,J+IOFF3,1),1)
                CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,3),1,
     &                                 FMOM(1+IOFF3,J+IOFF2,1),1)
                CALL DSCAL(NSTR(IFSYM),DM1,FMOM(1+IOFF3,J+IOFF2,1),1)
                CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,4),1,
     &                                 FMOM(1+IOFF2,J+IOFF3,2),1)
                CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,4),1,
     &                                 FMOM(1+IOFF3,J+IOFF2,2),1)
               ELSE
                CALL DCOPY(NSTR(IFSYM),FMOM(1+IOFF2,J+IOFF3,1),1,
     &                                 FMO (1+IOFF1,J+IOFF1,3),1)
                CALL DCOPY(NSTR(IFSYM),FMOM(1+IOFF2,J+IOFF3,2),1,
     &                                 FMO (1+IOFF1,J+IOFF1,4),1)
               ENDIF
            ENDDO
          ENDIF
        ENDIF
        IOFF1 = IOFF1 + NORB(IFSYM)
        IOFF2 = IOFF2 + 2*NSTR(IFSYM)
        IOFF3 = IOFF3 + NSTR(IFSYM)
      ENDDO
C
      IF (IPRINT.GE.7) THEN
        IF (FORWARD) THEN
         CALL HEADER(
     &   'MOLTODIR: transformation from DIRAC into MOLFDIR format',-1)
        ELSE
         CALL HEADER(
     &   'MOLTODIR: transformation from MOLFDIR into DIRAC format',-1)
        ENDIF
C       ....  
        WRITE(LUPRI,'(/3X,A)')
     &  'MOLTODIR:  *** matrix in DIRAC format, FMO ***'
        CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &               IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/3X,A)')
     &  'MOLTODIR:  *** matrix in MOLFDIR format, FMOM ***'
C       CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MINO(NZ,2),
        CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,2,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF
    
      CALL QEXIT('MOLTODIR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck moltodir */
      SUBROUTINE RELTODIR(XMO,XMOM,NFVO,IRCW)
C***********************************************************************
C     
C     Joost van Stralen, april 2002
C
C     XMO - The Lagrangian in Dirac format, this matrix is not sym. packed
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbtr3.h"
#include "dcbmp2.h"
C
      DIMENSION XMO(NORBT,NORBT,NZ)
      DIMENSION XMOM(NFVO,IRCW)
C
C
C
      CALL DZERO(XMO,NORBT*NORBT*NZ)
C
C
      IF (nz.NE.4) THEN
      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      DO IFSYM = 1, NFSYM
         IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
         IOFF2 = IOFF2 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
         IAI = 1
         DO I = 1, NAOCC(IFSYM)
            CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,1),1,
     &                              XMO(1+IOFF1,I+IOFF2,1),1)
         IAI = IAI + NAVIR(IFSYM)
         ENDDO
         IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN
          IAI = 1
          DO I = 1, NAOCC(IFSYM)
              CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,2),1,
     &                                XMO(1+IOFF1,I+IOFF2,2),1)
          IAI = IAI + NAVIR(IFSYM)
          ENDDO
         ENDIF
        IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) + NORB(IFSYM)
        IOFF2 = IOFF2 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) + NORB(IFSYM)
        IOFF3 = IOFF3 + 2*NAVIR(IFSYM)*NAOCC(IFSYM)
      ENDDO
C
      ELSE
C
      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      DO IFSYM = 1, NFSYM
        IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
        IOFF2 = IOFF2 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
        IAI = 1
        DO I = 1, NAOCC(IFSYM)
C    Re
           CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,1),1,
     &                             XMO(1+IOFF1,I+IOFF2,1),1)
C    i
           CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,2),1,
     &                             XMO(1+IOFF1,I+IOFF2,2),1)
           IAI = IAI + 2*NAVIR(IFSYM)
        ENDDO
        IAI = 2*NAVIR(IFSYM)*NAOCC(IFSYM) + 1
        DO I = 1, NAOCC(IFSYM)
C    j
             CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,1),1,
     &                               XMO(1+IOFF1,I+IOFF2,3),1)
C    k
             CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,2),1,
     &                               XMO(1+IOFF1,I+IOFF2,4),1)
             IAI = IAI + 2*NAVIR(IFSYM)
        ENDDO
        IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) + NORB(IFSYM)
        IOFF2 = IOFF2 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) + NORB(IFSYM)
        IOFF3 = IOFF3 + 4*NAVIR(IFSYM)*NAOCC(IFSYM)
      ENDDO
C
      ENDIF
C
      ioff1 = 0
      ioff2 = 0
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck dootodir */
      SUBROUTINE DOOTODIR(DMO,DMOM,NO,NFOO,NREP,IRCW)
C***********************************************************************
C     
C     Joost van Stralen, may 2002
C
C     Transform occ-occ part of 2nd order Density matrix from RELCCSD format 
C     to Dirac format
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbtr3.h"
C
      DIMENSION DMO(NORBT,NORBT,NZ)
      DIMENSION DMOM(NFOO,IRCW)
      DIMENSION NO(NREP)
C
C
      CALL DZERO(DMO,NORBT*NORBT*NZ)
C
      IF(NZ.NE.4) THEN
         IOFF1 = 0
         IOFF2 = 0
         IOFF4 = 0
         DO IFSYM = 1, NFSYM
            IOFF1 = IOFF1 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
            IJ = 1
            DO I = 1, NO(NFSYM+IOFF4)
               CALL DCOPY(NO(IFSYM+IOFF4),       DMOM(IJ+ioff2,1),1,
     &                                     DMO(1+IOFF1,I+IOFF1,1),1)
               IJ = IJ + NO(IFSYM+IOFF4)
            ENDDO
            IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN
               IJ = 1
               DO I = 1, NO(IFSYM+IOFF4)
                  CALL DCOPY(NO(IFSYM+IOFF4),       DMOM(IJ+IOFF2,2),1,
     &                                        DMO(1+IOFF1,I+IOFF1,2),1)
                  IJ = IJ + NO(IFSYM+IOFF4)
               ENDDO
            ENDIF
            IOFF1 = IOFF1 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) +NORB(IFSYM)
            IOFF2 = IOFF2 + 2*NO(IFSYM)*NO(IFSYM)
            IOFF4 = IOFF4 + 1
         ENDDO
      ELSE
         IOFF1 = 0
         IOFF2 = 0
         DO IFSYM = 1, NFSYM
            IOFF1 = IOFF1 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
C     Re       
            IJ = 1
            DO I = 1, (NO(IFSYM)/2)
               CALL DCOPY((NO(IFSYM)/2),      DMOM(IJ+IOFF2,1),1,
     &                                  DMO(1+IOFF1,I+IOFF1,1),1)
C     i
               CALL DCOPY((NO(IFSYM)/2),      DMOM(IJ+IOFF2,2),1,
     &                                  DMO(1+IOFF1,I+IOFF1,2),1)
            IJ = IJ + NO(IFSYM)
            ENDDO
C     j
            IJ = (NO(IFSYM)/2)*NO(IFSYM) + 1
            DO I = 1, (NO(IFSYM)/2)
               CALL DCOPY((NO(IFSYM)/2),      DMOM(IJ+IOFF2,1),1,
     &                                  DMO(1+IOFF1,I+IOFF1,3),1)
C     k
               CALL DCOPY((NO(IFSYM)/2),      DMOM(IJ+IOFF2,2),1,
     &                                  DMO(1+IOFF1,I+IOFF1,4),1)
            IJ = IJ + NO(IFSYM)
            ENDDO
            IOFF1 = IOFF1 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) +NORB(IFSYM)
            IOFF2 = IOFF2 + NO(IFSYM)*NO(IFSYM)
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck dvvtodir */
      SUBROUTINE DVVTODIR(DMO,DMOM,NV,NFVV,NREP,IRCW)
C*****************************************************************************
C     
C     Joost van Stralen, may 2002
C
C     Transform vir-vir part of 2nd order Density matrix from RELCCSD format 
C     to Dirac format
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbtr3.h"
C
      DIMENSION DMO(NORBT,NORBT,NZ)
      DIMENSION DMOM(NFVV,IRCW)
      DIMENSION NV(NREP)
C
      IF(NZ.NE.4) THEN
         IOFF1 = 0
         IOFF2 = 0
         IOFF4 = 0
         DO IFSYM = 1, NFSYM
            IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
            IJ = 1
            DO I = 1, NV(NFSYM+IOFF4)
               CALL DCOPY(NV(IFSYM+IOFF4),       DMOM(IJ+ioff2,1),1,
     &                                     DMO(1+IOFF1,I+IOFF1,1),1)
               IJ = IJ + NV(IFSYM+IOFF4)
            ENDDO
            IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN
               IJ = 1
               DO I = 1, NV(IFSYM+IOFF4)
                  CALL DCOPY(NV(IFSYM+IOFF4),       DMOM(IJ+IOFF2,2),1,
     &                                        DMO(1+IOFF1,I+IOFF1,2),1)
                  IJ = IJ + NV(IFSYM+IOFF4)
               ENDDO
            ENDIF
            IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) +NORB(IFSYM)
            IOFF2 = IOFF2 + 2*NV(IFSYM)*NV(IFSYM)
            IOFF4 = IOFF4 + 1
         ENDDO
      ELSE
         IOFF1 = 0
         IOFF2 = 0
         DO IFSYM = 1, NFSYM
            IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
C     Re       
            IJ = 1
            DO I = 1, (NV(IFSYM)/2)
               CALL DCOPY((NV(IFSYM)/2),      DMOM(IJ+IOFF2,1),1,
     &                                  DMO(1+IOFF1,I+IOFF1,1),1)
C     i
               CALL DCOPY((NV(IFSYM)/2),      DMOM(IJ+IOFF2,2),1,
     &                                  DMO(1+IOFF1,I+IOFF1,2),1)
            IJ = IJ + NV(IFSYM)
            ENDDO
C     j
            IJ = (NV(IFSYM)/2)*NV(IFSYM) + 1
            DO I = 1, (NV(IFSYM)/2)
               CALL DCOPY((NV(IFSYM)/2),      DMOM(IJ+IOFF2,1),1,
     &                                  DMO(1+IOFF1,I+IOFF1,3),1)
C     k
               CALL DCOPY((NV(IFSYM)/2),      DMOM(IJ+IOFF2,2),1,
     &                                  DMO(1+IOFF1,I+IOFF1,4),1)
            IJ = IJ + NV(IFSYM)
            ENDDO
            IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) +NORB(IFSYM)
            IOFF2 = IOFF2 + NV(IFSYM)*NV(IFSYM)
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dmp2ao */
      SUBROUTINE DMP2AO(DMO,DMAT,CMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C     
C     PURPOSE: Transform MO MP2 densitry matrix to quaternionic
C              symmetry-packed AO matrix 
C              
C     INPUT:   DMO  - MO density matrix
C              CMO  - spinor (MO) coefficients
C              
C     OUTPUT:  DMAT - AO MP2 denstity matrix
C     
C     REMARKS:
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D0 = 0.0D0)
C         
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION DMAT(N2BBASX,NZ),CMO(*),DMO(NORBT,NORBT,NZ),WORK(*)
C           
      EXTERNAL BDQALGEBRA
C     
      CALL DZERO(DMAT,N2BBASXQ)
C     
      DO IFSYM = 1,NFSYM
C     
C        Parameters for first index
C
          NS = NORB(IFSYM)
          IS = IORB(IFSYM) + 1
          JS = ICMOQ(IFSYM) + 1
C
C        Parameters for first index
C
          NI  = NORB(IFSYM) 
          II  = IORB(IFSYM) + 1
          JI  = ICMOQ(IFSYM) + 1
C         
          CALL QTRANS('MOAO','S',D0,
     &                NFBAS(IFSYM,0),NFBAS(IFSYM,0),NS,NI,
     &                DMAT(I2BASX(IFSYM,IFSYM)+1,1),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),
     &                DMO(IS,II,1),NORBT,NORBT,
     &                NZ,IPQTOQ(1,0),
     &                CMO(JS),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,1)
C     
      ENDDO
c       CALL HEADER('DMP2AO: AO density matrix',-1)
c       CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
c    &              NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
c       CALL FLSHFO(LUPRI)
C     
      RETURN
      END
