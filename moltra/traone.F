C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRAONE */
      SUBROUTINE TRAONE(WORK,LWORK)
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
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#else
#include <fpvm3h.h>
#endif
#include <mxorb.h>
#include <infpar.h>
#include <dcbtra.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbgen.h>
      LOGICAL TOBE,TRASAM
      DIMENSION WORK(*),NSTR(2,0:2,2),NSTRT(2),KVEC(2,2),KQ(2,2),
     &          IR(2,2),NR(2),NQ(2,2),NQT(2),
     &          NSPC(2,0:2),KQC(2),KE(2,2)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2),ICMOQC(2)
C
      CALL QENTER('TRAONE')
#include <memint.h>
      KFRSAV = KFREE
C
C     Get hold of the slaves if we run in parallel.
C
      IF (PARCAL) THEN
         ITASK = 1
#if defined (VAR_MPI)
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
         CALL PVMFPACK(INTEGER4,ITASK,1,1,INFO)
#endif
      ENDIF
C
      CALL TRAHI(IPRTRA,2)
C
C     Set up index arrays for core orbitals
C     =====================================      
C
      IC   = 1
      DO I = 1,NFSYM
        NR(I)   = NESH(I)
        IR(1,I) = 1
        IR(2,I) = NESH(I)
      ENDDO
C
      DO IFRP = 1,NFSYM
        CALL MEMGET('INTE',KVEC(IFRP,1),NR(IFRP),WORK,KFREE,LFREE)
        NSPC(IFRP,0) = 1
        CALL  NUMLST(TRA_CORSTR(IFRP),WORK(KVEC(IFRP,1)),
     &               NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &               NSPC(IFRP,0),WORK(KFREE),LFREE)
        CALL ORBCNT(WORK(KVEC(IFRP,1)),NSPC(IFRP,0),IFRP,
     &              NSPC(IFRP,2),NSPC(IFRP,1))
        NSPC(IFRP,0) = NSPC(IFRP,1) + NSPC(IFRP,2)
        CALL MEMREL('TRAONE',WORK,KFRSAV,KVEC(IFRP,1),KFREE,LFREE)
        CALL MEMGET('INTE',KVEC(IFRP,1),NSPC(IFRP,0),
     &              WORK,KFREE,LFREE)
      ENDDO
C
C     Print section
C     =============
C
      DO I = 1, NFSYM
         CALL TRAPRI(1,I,WORK(KVEC(I,1)),DUM,DUM,DUM,NSPC)
      ENDDO
C
C     Calculate dimensions of the coefficient array
C
C     ICMOQC : Index to start of representation in the coefficient array
C
      NQT(1) = 0
      ICMOQC(1) = 1
      DO IFRP = 1, NFSYM
         NQ(IFRP,1) = NFBAS(IFRP,0)*NSPC(IFRP,0)*NZ
         NQT(1) = NQT(1) + NQ(IFRP,1)
         IF (IFRP.LT.NFSYM) ICMOQC(IFRP+1) = ICMOQC(IFRP) + NQ(IFRP,1)
      ENDDO
C
C     Allocate the memory for the coefficients
C
      CALL MEMGET('REAL',KQC(1),NQT(1),WORK,KFREE,LFREE)
      KQC(2) = KQC(1) + NQ(1,1)
C
C     Get all coefficients
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(A)') 'Coefficient file not present'
         CALL QUIT('TRAPAM: No coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','TRAONE')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,TOTERG,3)
      ENDIF
C
C     Select the set that we need
C
      DO IFRP = 1, NFSYM
         CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),IFRP,WORK(KQC(IFRP)),
     &                WORK(KVEC(IFRP,1)),
     &                NSPC(IFRP,2),NSPC(IFRP,1),
     &                NFBAS(IFRP,0),NORB(IFRP))
      ENDDO
C
C     Throw away the complete set, keep only the selected sets
C
      CALL MEMREL('TRAONE',WORK,KFRSAV,KCMO,KFREE,LFREE)
C
C     =======================================      
C     Set up index arrays for active orbitals
C     =======================================      
C
      IF(NOPAIR) THEN
        IC   = 1
        DO I = 1,NFSYM
          NR(I)   = NESH(I)
          IR(1,I) = 1
          IR(2,I) = NESH(I)
        ENDDO
      ELSE
        IC   = 0
        DO I = 1,NFSYM
          NR(I)   =  NORB(I)
          IR(1,I) = -NPSH(I)
          IR(2,I) =  NESH(I)
        ENDDO
      ENDIF
C
      DO I = 1,2
        NSTRT(I) = 0
        DO IFRP = 1,NFSYM
          CALL MEMGET('INTE',KVEC(IFRP,I),NR(IFRP),WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = 1
          CALL  NUMLST(TRA2_INDSTR(I,IFRP),WORK(KVEC(IFRP,I)),
     &                 NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                 NSTR(IFRP,0,I),WORK(KFREE),LFREE)
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),IFRP,
     &              NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NSTRT(I) = NSTRT(I) +  NSTR(IFRP,0,I)
          CALL MEMREL('TRAONE',WORK,KFRSAV,KVEC(IFRP,I),KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &                WORK,KFREE,LFREE)
        ENDDO
      ENDDO      
C
C     Print section
C     =============
C
      DO I = 1, NFSYM
         CALL TRAPRI(2,I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),DUM,DUM,NSTR)
      ENDDO
C
C     Calculate dimensions of coefficient arrays
C     The coefficients are stored as 2 arrays of two matrices each
C     They may share the same memory if the transformation ranges
C     are identical
C
C     NDMOQR : 1st dimension row/column dimension of coefficient array
C              2nd dimension fermion irrep
C              3rd dimension transformation index
C
C     ICMOQR : Index to start of representation in the coefficient array
C
      DO I = 1, 2
         NQT(I) = 0
         ICMOQR(1,I) = 1
         DO IFRP = 1, NFSYM
            NDMOQR(1,IFRP,I) = NFBAS(IFRP,0)
            NDMOQR(2,IFRP,I) = NSTR(IFRP,0,I)
            NQ(IFRP,I) = NFBAS(IFRP,0)*NSTR(IFRP,0,I)*NZ
            NQT(I) = NQT(I) + NQ(IFRP,I)
            IF (IFRP.LT.NFSYM) 
     &         ICMOQR(IFRP+1,I) = ICMOQR(IFRP,I) + NQ(IFRP,I)
         ENDDO
      ENDDO
C
C     Allocate the memory for the coefficients, check their relations
C
      DO I = 1, 2
C
         DO J = I, 1, -1
            IF (TRASAM(WORK(KVEC(1,I)),WORK(KVEC(2,I)),
     &                 WORK(KVEC(1,J)),WORK(KVEC(2,J)),
     &          NSTR(1,0,I),NSTR(1,0,J))) ISAME(I) = J
         ENDDO
C
         IF (ISAME(I).NE.I) THEN
            KQ(1,I) = KQ(1,ISAME(I))
            KE(1,I) = KE(1,ISAME(I))
         ELSE
            CALL MEMGET('REAL',KQ(1,I),NQT(I),WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KE(1,I),NSTRT(I),WORK,KFREE,LFREE)
         ENDIF
         KQ(2,I) = KQ(1,I) + NQ(1,I)
         KE(2,I) = KE(1,I) + NSTR(1,0,I)
C
      ENDDO
C
C     Get all coefficients and eigenvalues
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(A)') 'Coefficient file not present'
         CALL QUIT('TRAPAM: No coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','TRAONE')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),TOTERG,2)
      ENDIF
C
C     Select the set that we need
C
      DO I = 1, 2
         DO IFRP = 1, NFSYM
            CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),IFRP,WORK(KQ(IFRP,I)),
     &                   WORK(KVEC(IFRP,I)),
     &                   NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &                   NFBAS(IFRP,0),NORB(IFRP))
            CALL SELEIG (WORK(KEIG+IORB(IFRP)),IFRP,WORK(KE(IFRP,I)),
     &                   WORK(KVEC(IFRP,I)),
     &                   NSTR(IFRP,2,I),NSTR(IFRP,1,I))
         ENDDO
      ENDDO
C
C     Throw away the complete set, keep only the selected sets
C
      CALL MEMREL('TRAONE',WORK,KFRSAV,KCMO,KFREE,LFREE)
C
C     Transform core-fock matrix to active MO basis
C     ---------------------------------------------
C
      CALL MEMGET('REAL',KFMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL GTCFCK(WORK(KFMO),IPRTRA,WORK(KQ(1,1)),WORK(KQ(1,2)),
     &            WORK(KQC(1)),NDMOQR,ICMOQR,NSTR,ICMOQC,NSPC,
     &            ITRA_INTFL2,WORK,KFREE,LFREE)
C
C     Write MOLFDIR file MRCONEE
C     --------------------------
C
      CALL MRCONE(WORK(KFMO),WORK(KE(1,1)),NSTR,NSTRT,IPRTRA,
     &            WORK,KFREE,LFREE)
C
C     Release the slaves if we run in parallel.
C
      IF (PARCAL) THEN
         ITASK = -1
#if defined (VAR_MPI)
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
         CALL PVMFPACK(INTEGER4,ITASK,1,1,INFO)
#endif
      ENDIF
C
      CALL MEMREL('TRAONE',WORK,KFRSAV,KFMO,KFREE,LFREE)
      CALL QEXIT('TRAONE')
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gtcfck */
      SUBROUTINE GTCFCK(FMO,IPRINT,Q1,Q2,QCO,NDMOQR,ICMOQR,NSTR,ICMOQC,
     &                  NSPC,INTFLG2,WORK,KFREE,LFREE)
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
      DIMENSION ICMOQC(2),NSPC(2,0:2)
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
     &            NSTR,ICMOQC,NSPC,IPRINT,WORK,KFREE,LFREE)
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
     &                  ICMOQC,NSPC,IPRINT,WORK,KFREE,LFREE)
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
#include <cbihr2.h>
      LOGICAL TOBE,BF1,BF2,BF3
      DIMENSION FMO(*),FAO(*),DMAT(*),WORK(*)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2),NSTR(2,0:2,2)
      DIMENSION ICMOQC(2),NSPC(2,0:2),IQBF(4)
      DIMENSION Q1(*),Q2(*),QCO(*)
C
      CALL QENTER('GTCFC1')
      KFRSAV = KFREE
      LUBUF = 22
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
        NFMAT = 1
C       Totally symmetric operator
        ISYMOP(1) = 1
C       Fock matrix type
        IFCKOP(1) = 1
C       Hermitian operator
        IHRMOP(1) = 1
C       Get core density matrix
        CALL CDENS(DMAT,QCO,ICMOQC,NSPC,IPRINT)        
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
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GTCFC1:Total Fock matrix in AO-basis',-1)
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
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &            NSTR(I,0,1),NSTR(I,0,2),
     &            FAO(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            FMO(I2ORBX(I,I)+1),NORBT,NORBT,NBF,IQBF,
     &            Q1(ICMOQR(I,1)),NRQ1,NCQ1,NZ,IPQTOQ(1,0),
     &            Q2(ICMOQR(I,2)),NRQ2,NCQ2,NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)
        IF(NBF.NE.NZ) WRITE(LUPRI,'(A)')
     &     '*** WARNING *** GTCFC1: NBF.NE.NZ !'
 10   CONTINUE
C
      CALL QEXIT('GTCFC1')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck CDENS */
      SUBROUTINE CDENS(DMAT,QCO,ICMOQC,NSPC,IPRINT)
C*****************************************************************************
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION DMAT(*),QCO(*)
      DIMENSION NSPC(2,0:2),ICMOQC(2)
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
C
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NSPC(I,0).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),
     &              QCO(ICMOQC(I)),NFBAS(I,0),NSPC(I,0),
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
      SUBROUTINE MRCONE(FMO,ORBE,NSTR,NSTRT,IPRINT,WORK,KFREE,LFREE)
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
      DIMENSION FMO(NORBT,NORBT,NZ),ORBE(*),WORK(*)
      DIMENSION NSTR(2,0:2,2),NSTRT(2)
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
C
      CALL MRCON1(FMO,ORBE,WORK(KFMOM),NSTR,NSTRT,WORK(KIRPMO),
     &            WORK(KIRPAMO),WORK(KORBMO),IPRINT,WORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mrcon1 */
      SUBROUTINE MRCON1(FMO,ORBE,FMOM,NSTR,NSTRT,IRPMO,IRPAMO,ORBMO,
     &                  IPRINT,WORK,KFREE,LFREE)
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
#include <dcbtra.h>
#include <dgroup.h>
C
      LOGICAL BREIT
      CHARACTER*14 REPN(8)
      CHARACTER*4 REPA(8)
      DIMENSION FMO(NORBT,NORBT,NZ),ORBE(NSTRT)
      DIMENSION FMOM(NSTRT*2,NSTRT*2,2),WORK(*)
      DIMENSION IRPMO(NSTRT*2),IRPAMO(NSTRT*2),ORBMO(NSTRT*2)
      DIMENSION NSTR(2),IPHS(4)
      DIMENSION MULTB(8,8)
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
         REPN(1) = '   E  1E          '
         REPN(2) = '   E  2E          '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.4) THEN
         NSYMRP = 4
         REPN(1) = ' 1AG          '
         REPN(2) = ' 1AG          '
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
         REPN(1) = '  Eg 1Eg          '
         REPN(2) = '  Eg 2Eg          '
         REPN(3) = '  Eu 1Eu          '
         REPN(4) = '  Eu 2Eu          '
      ELSE
         PRINT*, 'Warning representation name is not defined'
      ENDIF
C
C     Get multiplication table for the abelian subgroup
C
      CALL GMULTA (NZ,NFSYM,NSYMRPA,REPA,MULTB)
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
            IRPAMO(II) = IREPA
            ORBMO(II)  = ORBE(IK)
         ENDDO
         IREP = IREP + 1
         IF (NZ.NE.4) IREPA = IREPA + 1
         IK = IK - NSTR(I)
         DO ISTR = 1, NSTR(I)
            II = II + 1
            IK = IK + 1
            IRPMO(II)  = IREP
            IRPAMO(II) = IREPA
            ORBMO(II)  = ORBE(IK)
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
        IF (NZ.GE.2) THEN
          DO J = 1, NSTR(IFSYM)
            CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,2),1,
     &                             FMOM(1+IOFF2,J+IOFF2,2),1)
            CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,2),1,
     &                             FMOM(1+IOFF3,J+IOFF3,2),1)
            CALL DSCAL(NSTR(IFSYM),DM1,FMOM(1+IOFF3,J+IOFF3,2),1)
          ENDDO
          IF (NZ.GE.4) THEN
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
      IF (IPRINT.GE.3) THEN
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Core Fock matrix ***'
          CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZ,2),
     &                IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Write MRCONEE file
C
      OPEN (LUMLF1,FILE='MRCONEE',FORM='UNFORMATTED')
      WRITE (LUMLF1) 2*NSTRT,BREIT,DHFERG
      WRITE (LUMLF1) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
      WRITE (LUMLF1) NSYMRPA,(REPA(IRPA),IRPA=1,NSYMRPA*2)
      WRITE (LUMLF1) ((MULTB(I,J),I=1,2*NSYMRPA),J=1,2*NSYMRPA)
      WRITE (LUMLF1) (IRPMO(I),IRPAMO(I),ORBMO(I),I=1,NSTRT*2)
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
      CHARACTER*4 REPA(8)
      DIMENSION MULTB(8,8)
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
         STOP 'Multiplication table is not defined'
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
