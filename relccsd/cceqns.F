!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MP2EQN(EPS,FVO,T1,S1,T2,S2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Solves MP2 amplitude or lambda equations.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C     Intermediaries etc.
C
      REAL*8 EPS(*),FVO(*),T1(*),S1(*),T2(*),S2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "results.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 DIAG
      LOGICAL TEQ
C
C---------------Executable code--------------------------------------
C
      TEQ = .TRUE.
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
      ELSEIF (EQNS.EQ.'FOCKSP') THEN
         TEQ = .TRUE.
      ELSE
         CALL QUIT('Wrong keyword in call to MP2EQNS')
      ENDIF
C
C  INITIALIZE T2 AND T1 AMPLITUDES
C
      CALL XCOPY (NDIMT1,FVO,1,T1,1)
      CALL GETVVOO (T2)
      CALL DENOM (EPS,T1,T2,T1,T2)
      CALL ZCORE (T1,T2)
C
      IF (TEQ) THEN
C
C  EVALUATE CORRELATION ENERGY 
C
         CALL MP2ENG (T1,T2,FVO,S2,EMP2)
         CALL T1DIAG (DIAG,T1,S1)
C
         WRITE(IW,6071)
         WRITE(IW,6072) ESCF
         WRITE(IW,6073) EMP2
         ETOT = ESCF + EMP2
         WRITE(IW,6075) ESCF + EMP2
         WRITE(IW,6077) DIAG
C
      ELSE
C
C  INITIALIZE LAMDA AMPLITUDES
C
         CALL XCOPY (NDIMT1,T1,1,S1,1)
         CALL XCOPY (NDIMT2,T2,1,S2,1)
C
      ENDIF
C
      RETURN
 6071 FORMAT(//,'  MP2 results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 7000 FORMAT(//' Timing of routine MP2EQN :'
     &/' Before MP2EQN :',T30,F12.3,' seconds'
     &/' Energy, T1-Diag :',T30,F12.3,' seconds'
     &/' Untimed parts :',T30,F12.3,' seconds'
     &/' Total time in MP2EQN :',T30,F12.3,' seconds')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCEQN(MAXIT,EPS,FOO,FVO,FVV,T1,S1,T2,S2,AZ,CZ,H,
     &                 HO,HV,HOV,GOO,GVV,
     &                 BUF1,BUF2,BUF3,CC1,BB1,BB2,CONV,
     &                 NBUF2,NBUF3,HEFF,EIGENVALUES,EIGENVECTORS)
C
      use interface_to_mpi
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Solves CCSD amplitude or lambda equations.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C     Intermediaries etc.
C
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*),T1(*),S1(*),T2(*),S2(*)
      REAL*8 AZ(*),H(*),HO(*),HV(*),HOV(*),BUF1(*),BUF2(*),BUF3(*)
      REAL*8 CZ(*) ! is only allocated for Lambda equations
      REAL*8 GOO(*),GVV(*)
C
C     for diis ...
C
      REAL*8 CC1(*)
      REAL*8 BB1(*),BB2(*)
C
      REAL*8 RMS,CONV
      REAL*8 DIAG
C
C     for Fock space ...
C
      REAL*8 HEFF(*),EIGENVECTORS(*),EIGENVALUES(0:*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "ihm.inc"
#include "results.inc"
#include "complex.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
      DIMENSION ISTAT(df_MPI_STATUS_SIZE)
#endif
C
C---------------Local variables--------------------------------------
C
      REAL*8 CPU,CPU4V,CPUH1,CPUH2
      LOGICAL TEQ,FOCKSP,PRENE
      CHARACTER*2 DELS1
      CHARACTER*4 DELS2
      REAL*8 DDOT
C
C---------------Executable code--------------------------------------
C
C
C     Get initial values of timers for the most important contractions
C     This is important in Fock space calculations where we enter this
C     routine more than once. We do it in this manner rather than resetting
C     the timer because we want the overall timing at the end.
C
      CALL XTIMRD (CPU4V,'-- BINTM                       ')
      CALL XTIMRD (CPUH1,'--- HINTM: VVOO contribution   ')
      CALL XTIMRD (CPUH2,'--- T2EQNS: HINTM*T2           ')
C
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
         FOCKSP = .FALSE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
         FOCKSP = .FALSE.
      ELSEIF (EQNS.EQ.'FOCKSP') THEN
         TEQ = .TRUE.
         FOCKSP = .TRUE.
         DELS1 = 'DD'
         IF (ICURFSS.EQ. 2) DELS1 = 'DK'
         IF (ICURFSS.EQ. 3) DELS1 = 'KD'
         IF (ICURFSS.EQ. 4) DELS1 = 'KK'
         IF (ICURFSS.EQ. 5) DELS1 = 'DK'
         IF (ICURFSS.EQ. 6) DELS1 = 'KD'
         DELS2 = 'DDDD'
         IF (ICURFSS.EQ. 2) DELS2 = 'DDDK'
         IF (ICURFSS.EQ. 3) DELS2 = 'KDDD'
         IF (ICURFSS.EQ. 4) DELS2 = 'KDDK'
         IF (ICURFSS.EQ. 5) DELS2 = 'DDKK'
         IF (ICURFSS.EQ. 6) DELS2 = 'KKDD'
C
C        Initialize S1 and S2 in case we want to do zero
C        iterations (tests). In normal runs S1 and S2 will
C        be initialized in T1EQNS and T2EQNS.
C
         CALL XCOPY(NDIMT1,A0,0,S1,1)
         CALL XCOPY(NDIMT2,A0,0,S2,1)
         CALL S1S2_ADDFO (FVO,BUF1,BUF2,S1,S2)
C
      ELSE
         CALL QUIT('Wrong keyword in call to CCEQNS')
      ENDIF
C
C     Zero amplitudes if the NOCCX options are specified
C
      IF (TEQ.AND.NOCCS) CALL XCOPY (NDIMT1,A0,0,T1,1)
      IF (TEQ.AND.NOCCD) CALL XCOPY (NDIMT2,A0,0,T2,1)
C
C     Initialize the work arrays to zero.
C
      CALL XCOPY (NDIMT2,A0,0,BUF1,1)
      CALL XCOPY (NBUF2,A0,0,BUF2,1)
      CALL XCOPY (NBUF3,A0,0,BUF3,1)
C
C     Initialize DIIS.
C
      IT = -1
      CALL DIISX(NDIMT1,NDIMT2,BUF1,T1,BUF2,T2,
     &           IT,CC1,BB1,ITC,BB2,RMS)
C
      ECCSD = 0.D0
      NIT = 0
C
      CALL T1DIAG(DIAG,T1,BUF1)
      CALL CCENG (T1,T2,FVO,BUF1,BUF2,ECCSD,ECCSDIM)
      IF (IPRNT.GE.3) CALL PAIREN (T1,T2,FVO,BUF1,BUF2)
C
      IF (IPRNT.GE.1) THEN
         WRITE(IW,6003)
         WRITE(IW,6001) NIT,ECCSD,RMS,DIAG
      ENDIF
C
C     Calculate fixed intermediaries for Lambda equations
C
      IF (.NOT.TEQ) THEN
         CALL GETAMPT (S1,S2)
         CALL HOINTM (FOO,EPS,S1,S2,BUF1,BUF2,BUF3,HO)
         CALL HVINTM (FVV,EPS,S1,S2,BUF1,BUF2,BUF3,HV)
         CALL HOVINTM (FVO,S1,BUF1,BUF2,HOV)
         CALL GOINTM (FVO,S1,BUF1,BUF2,HO)
         CALL GVINTM (FVO,S1,BUF3,NBUF3,BUF2,NBUF2,HV)
         CALL AINTM (S1,S2,BUF1,BUF2,AZ)
         CALL CINTM (T1,T2,FVO,AZ,BUF1,BUF2,BUF3,NBUF3,CZ)
         CALL DINTM (S1,S2,HOV,BUF1,BUF2,NBUF2,BUF3,NBUF3,H)
         CALL HINTM (S1,S2,BUF1,BUF2,BUF3,NBUF3,H)
      ENDIF
C
C  START ITERATIVE PROCEDURE
C
 1000 CONTINUE
      NIT = NIT + 1
      IF (NIT.GT.MAXIT) THEN
         NIT = NIT - 1
         WRITE(IW,6004) RMS
         GO TO 2000
      ENDIF
C
      IF (.NOT.TEQ) THEN
C        
C        Make Gpq intermediates for lambda equations
C
         CALL GETAMPT (S1,S2)
         CALL MKGPP (S2,T2,BUF1,BUF2,GOO,GVV)
      ENDIF

C  EVALUATE T1 EQUATIONS
C  CHECK FOR CASES WITH NO NON_ZERO T1'S TO AVOID CRASHES
C
      CALL XTIME(5,1,'- T1 equations                  ')
      IF (TEQ) THEN
           CALL HOINTM (FOO,EPS,T1,T2,BUF1,BUF2,BUF3,HO)
      ENDIF
      IF (TEQ) THEN
           CALL HVINTM (FVV,EPS,T1,T2,BUF1,BUF2,BUF3,HV)
      ENDIF
      IF (NDIMT1.GT.0) THEN
         IF (TEQ) THEN
          CALL HOVINTM (FVO,T1,BUF1,BUF2,HOV)
          CALL T1EQNS (FVO,EPS,T1,T2,HO,HV,HOV,H,CZ,GOO,GVV,
     &                 BUF1,BUF2,NBUF2,BUF3,NBUF3,S1)
         ELSE
          CALL T1EQNS (HOV,EPS,T1,T2,HO,HV,HOV,H,CZ,GOO,GVV,
     &                 BUF1,BUF2,NBUF2,BUF3,NBUF3,S1)
         ENDIF
      ENDIF
      CALL XTIME(5,2,'- T1 equations                  ')

C
C     In case of NOCCD we can skip the calculation of T2 amplitudes 
C
      IF (TEQ.AND.NOCCD) GOTO 100
      CALL XTIME(5,1,'- T2 equations                  ')

      
CMI  ... debugging output ...
      IF (DEBUG) THEN
        WRITE (*,111) "ho ",MYPROC,DDOT(RCW*NFOO,HO,1,HO,1)
        WRITE (*,111) "hv ",MYPROC,DDOT(RCW*NFVV,HV,1,HV,1)
        WRITE (*,111) "hov",MYPROC,DDOT(RCW*NDIMT1,HOV,1,HOV,1)
        WRITE (*,111) "s1 ",MYPROC,DDOT(RCW*NDIMT1,S1,1,S1,1)
      ENDIF
C
C  EVALUATE T2 INTERMEDIATES
C
      IF (TEQ) THEN
         CALL XTIME(0,1,'-- GOINTM                        ')
         CALL GOINTM (FVO,T1,BUF1,BUF2,HO)
         CALL XTIME(0,2,'-- GOINTM                        ')
         CALL XTIME(0,1,'-- GVINTM                        ')
         CALL GVINTM (FVO,T1,BUF3,NBUF3,BUF2,NBUF2,HV)
         CALL XTIME(0,2,'-- GVINTM                        ')
         CALL XTIME(0,1,'-- AINTM                         ')
         CALL AINTM (T1,T2,BUF1,BUF2,AZ)
         CALL XTIME(0,2,'-- AINTM                         ')
         CALL XTIME(0,1,'-- HINTM                         ')
         CALL HINTM (T1,T2,BUF1,BUF2,BUF3,NBUF3,H)
         CALL XTIME(0,2,'-- HINTM                         ')

#if defined (VAR_MPI)
C *******************
C ** parallel case: the H-intermediate was partially calculated
C ** on the separate slaves. In order to perform the H*T2
C ** contraction we need the complete H array ==> we complete it
C ** here and synchronize all the nodes. They then have the complete
C ** H-array and this lengthy contraction can also be parallelized.
C
      IF (NMPROC .GT. 1) THEN

        IF(NBUF3.LT.IVOVO(NREP+1)) 
     &    CALL QUIT('BUF3 can not keep the VOVO integrals !')

C ... we collect the H array

        CALL XTIME(4,1,'-- combining HINTM via MPI_ALLREDUCE')
        CALL XCOPY (NBUF3,A0,0,BUF3,1)
CSK        WRITE(IW,*) 'THIS IS IERR before',IERR
        IERR = 0
        call interface_mpi_allreduce_r1_work_f77(H,BUF3,
     &       RCW*IVOVO(NREP+1),
     &       op_MPI_SUM,global_communicator)
        IF(IERR.GT.0) THEN
           CALL QUIT('MPI_REDUCE ERROR in collecting HINTM !')
        ENDIF
        CALL XCOPY(IVOVO(NREP+1),BUF3,1,H,1)
        CALL XTIME(4,2,'-- combining HINTM via MPI_ALLREDUCE')

      ENDIF

#endif

      ENDIF

C
      IF (DEBUG) THEN
         WRITE (*,111)"go ",MYPROC,DDOT(RCW*NFOO,HO,1,HO,1)
         WRITE (*,111)"gv ",MYPROC,DDOT(RCW*NFVV,HV,1,HV,1)
         LEN = IOOOOTT(NREP+1)
         WRITE (*,111)"az ",MYPROC,4.D0*DDOT(RCW*LEN,AZ,1,AZ,1)
         LEN = IVOVO(NREP+1)
         WRITE (*,111)"h  ",MYPROC,DDOT(RCW*LEN,H,1,H,1)
      ENDIF
C
C  INTERMEDIATES FORMED --- START T2 EQUATIONS
C
      CALL XTIME(5,1,'-- T2 EQNS                     ')
      CALL T2EQNS (EPS,T1,T2,HO,HV,HOV,H,AZ,CZ,GOO,GVV,
     &             BUF1,BUF2,BUF3,NBUF3,S2)
      CALL XTIME(5,2,'-- T2 EQNS                     ')
C
C
C  MUST NOW ADD IN B INTERMEDIATE CONTRIBUTION
C
      IF (DEBUG) THEN
         WRITE (*,111)"s2 ",MYPROC,DDOT(RCW*NDIMT2,S2,1,S2,1)
      ENDIF
C      
      CALL XTIME(5,1,'-- BINTM                       ')
      CALL BINTM (T1,T2,BUF1,BUF2,NBUF2,BUF3,NBUF3,S2,IW)
      CALL XTIME(5,2,'-- BINTM                       ')

      IF (DEBUG) THEN
         WRITE (*,111)"s2 after bintm",MYPROC,DDOT(RCW*NDIMT2,S2,1,S2,1)
      ENDIF
C
      CALL XTIME(5,2,'- T2 equations                  ')
C
 100  CONTINUE
C
C     Delete the Pi->Q excitations in case of Int. Ham. CCFS scheme 2
      IF (FOCKSP.AND.DOIH.AND.IHSCHEME.EQ.2) THEN
         CALL DENOMF (EPS,S1,S2,S1,S2,3)
      ENDIF
C
C     Add the first order wave function contribution
C
      CALL S1S2_ADDFO (FVO,BUF1,BUF2,S1,S2)
C
C  ADD PARTIAL S1 AND S2 (PARALLEL RUNS)
C
      CALL XTIME(4,1,'- adding partial T1/T2 amplitudes ')
      CALL PARTS (1,RMS,S1,S2,BUF1,BUF2)
      CALL XTIME(4,2,'- adding partial T1/T2 amplitudes ')
C
      IF (DEBUG) THEN
        WRITE (*,111)"s2 after parts",MYPROC,DDOT(RCW*NDIMT2,S2,1,S2,1)
        WRITE (*,111)"s1 after parts",MYPROC,DDOT(RCW*NDIMT1,S1,1,S1,1)
      ENDIF
C
      IF (FOCKSP) THEN
C       Adding folded diagrams, use HOV and H as scratch arrays
        CALL FOLDED (ICURFSS,IHEQ,T1,T2,BUF1,BUF2,BUF3,HOV,H,
     &               S1,S2)     
        IF (DEBUG) THEN
        WRITE (*,111)"s1 after folded",MYPROC,DDOT(RCW*NDIMT1,S1,1,S1,1)
        WRITE (*,111)"s2 after folded",MYPROC,DDOT(RCW*NDIMT2,S2,1,S2,1)
        ENDIF
      ENDIF
C
C  Make new T1s and T2 by dividing S1s and S2 by eigenvalues
C
      CALL XCOPY (NDIMT1,S1,1,T1,1)
      CALL XCOPY (NDIMT2,S2,1,T2,1)
C
      IF (FOCKSP) THEN
        CALL DENOMF (EPS,T1,T2,T1,T2,1)
C       Remove the Hamiltonian elements, T1 and T2 should only contain
C       the wave function (amplitudes) appropriate for this sector.
        CALL DENOMF (EPS,T1,T2,T1,T2,2)
      ELSE
        CALL DENOM (EPS,T1,T2,T1,T2)
      ENDIF
C
C     Zero amplitudes if the NOCCX options are specified
C
      IF (TEQ.AND.NOCCS) CALL XCOPY (NDIMT1,A0,0,T1,1)
      IF (TEQ.AND.NOCCD) CALL XCOPY (NDIMT2,A0,0,T2,1)
C
      CALL ZCORE (T1,T2)
C
C.....when debugging print out energies in every iteration
      IF (FOCKSP.AND.DEBUG) 
     &   CALL HEFF1(ICURFSS,1,S1,S2,FVO,EPS,BUF1,BUF2,BUF3,
     &              HEFF,EIGENVALUES,EIGENVECTORS)
C
C  EXTRAPOLATE THE AMPLITUDES
C
      CALL XTIME(5,1,'- DIIS extrapolation            ')
      CALL DELFCK ('VO',DELS1,T1)
      CALL DELINT ('VVOO',DELS2,T2,0,0,0)
      CALL DIISX(NDIMT1,NDIMT2,BUF1,T1,BUF2,T2,
     &           IT,CC1,BB1,ITC,BB2,RMS)
      CALL XTIME(5,2,'- DIIS extrapolation            ')
C
C  SYNCHRONIZE T1 AND T2 (PARALLEL RUNS)
C
      CALL XTIME(4,1,'- synchronizing T1 & T2 amplitudes')
      CALL PARTS (2,RMS,T1,T2,BUF1,BUF2)
      CALL XTIME(4,2,'- synchronizing T1 & T2 amplitudes')
C
C  EVALUATE CORRELATION ENERGY AND WRITE OUT RESULTS
C
      IF (TEQ) THEN
         IF (IPRNT.GE.3) CALL PAIREN (T1,T2,FVO,BUF1,BUF2)
         CALL CCENG  (T1,T2,FVO,BUF1,BUF2,ECCSD,ECCSDIM)
      ENDIF
C
      CALL T1DIAG (DIAG,T1,BUF1)
      IF (IPRNT.GE.1) WRITE(IW,6001) NIT,ECCSD,RMS,DIAG
      CALL FLSHFO(IW)
C
C
      IF(RMS.LT.CONV) GO TO 2000
      GO TO 1000
C
 2000 CONTINUE
      CALL T1DIAG (DIAG,T1,BUF1)
      WRITE(IW,6071)
C
C     Only print energy when doing amplitudes for the reference
C     wave function.
C
      IF (.NOT.FOCKSP) THEN
         PRENE =  TEQ
      ELSE
         PRENE = ICURFSS .EQ. 1
      ENDIF
C
      IF (PRENE) THEN
         ETOT = ESCF + ECCSD
         WRITE(IW,6072) ESCF
         IF (.NOT.(NOCCD.OR.NOCCS)) THEN
            WRITE(IW,6074) ECCSD
         ELSE
            IF (NOCCS) WRITE(IW,6084) ECCSD
            IF (NOCCD) WRITE(IW,6085) ECCSD
         ENDIF
         WRITE(IW,6076) ETOT
         IF (IPRNT.GE.2) CALL PAIREN (T1,T2,FVO,BUF1,BUF2)
      ENDIF
C
      WRITE(IW,6077) DIAG
      WRITE(IW,6078) RMS
      WRITE(IW,6079) NIT
C
      IF (FOCKSP) THEN
C        Calculate the energies for this sector
         CALL HEFF1(ICURFSS,IPRNT,S1,S2,FVO,EPS,BUF1,BUF2,BUF3,
     &              HEFF,EIGENVALUES,EIGENVECTORS)
      ENDIF
C
C     Read out timers for the most important contractions
C
      CALL XTIMRD (CPU,'-- BINTM                       ')
      CPU4V = CPU - CPU4V
      CALL XTIMRD (CPU,'--- HINTM: VVOO contribution   ')
      CPUH1 = CPU - CPUH1
      CALL XTIMRD (CPU,'--- T2EQNS: HINTM*T2           ')
      CPUH2 = CPU - CPUH2
C
C     Compute & print the performance in Mega/Giga/TeraFlops
C
      IF (TIMING) CALL BFLOPS (CPU4V,CPUH1,CPUH2,NIT)
C
C     Close and delete DIIS scratch files
C
      IT = -2
      CALL DIISX(NDIMT1,NDIMT2,BUF1,T1,BUF2,T2,
     &           IT,CC1,BB1,ITC,BB2,RMS)
C
 111  FORMAT(A20,I4,2f25.15)
 6001 FORMAT(2X,I3,5X,2F20.15,F10.5)
 6002 FORMAT(3X,'MP2   =',F20.15/3X,'EMP2  =',F20.15/)
 6003 FORMAT(//3X,'NIT',10X,'ENERGY',20X,'RMS',7X,'T1-DIAGN')
 6004 FORMAT(//' WARNING: Requested convergence not reached in CCSD :',
     &         'RMS = ',F20.15)
 6071 FORMAT(//,'  CCSD results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6074 FORMAT(' CCSD correlation energy :',T40,F25.15)
 6084 FORMAT(' CCD correlation energy :',T40,F25.15)
 6085 FORMAT(' CCS correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6076 FORMAT(' Total CCSD energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 6078 FORMAT(' Convergence :',T40,F25.15)
 6079 FORMAT(' Number or iterations used :',T60,I5)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T1EQNS (FVO,EPS,T1,T2,HO,HV,HOV,H,CZ,GOO,GVV,
     &                   BUF1,BUF2,NBUF2,BUF3,NBUF3,S1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate T1 or Lambda1 equations
C     Here we have to take care of the fact that certain integral classes
C     are now split over individual nodes. The corresponding contributions
C     which are linear in the integrals are added at the end of the T1
C     routines. Therefore we have to take care of calculating the 
C     non-distributed contributions to the T1 amplitudes only on the
C     master node in order to avoid double-counting in the PARTS routine !
C     all direct VOVV/VVVV contributions and VOVV/VVVV contributions stemming 
C     from intermediates have to be calculated on EACH node.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF2,NBUF3
      REAL*8 EPS(*),FVO(*)
      REAL*8 T1(*),T2(*)
      REAL*8 S1(*)
      REAL*8 HO(*),HV(*),HOV(*),H(*)
      REAL*8 GOO(*),GVV(*),CZ(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ihm.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ,USEDZ,RIGHT
      REAL*8 SDOT,DDOT
C
C---------------Executable code--------------------------------------
C
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
      ELSEIF (EQNS.EQ.'FOCKSP') THEN
         TEQ = .TRUE.
      ELSE
         CALL QUIT('Wrong keyword in call to T1EQNS')
      ENDIF
C-------------------
C Initialize to zero
C-------------------
      CALL XCOPY (NDIMT1,A0,0,S1,1)
C--------------------------------
C Terms specific for T1 equations
C--------------------------------
      IF (TEQ) THEN
C ********************
C ** TASK SPLITTING **
C ********************
        IF (MYPROC.EQ.MASTER) THEN
C------------------------------------------------------------
C S(A,I) = S(A,I) + T(A,K) * [HOV(K,C) - 2 * F(K,C)] * T(C,I)
C------------------------------------------------------------
         CALL XTIME(0,1,'--- T1EQNS T*[HOV - F]*T      ')
         CALL XCOPY (NFVO,HOV,1,BUF1,1)
         IF (CARITH) CALL CONJUGA (NFVO,BUF1,1)
         CALL XAXPY (NFVO,-A2,FVO,1,BUF1,1)
         CALL DELFCK ('VO','DD',BUF1)
         OFF1 = 1
         DO IRP = 1, NREP
            M = NV(IRP)
            N = NO(IRP)
            K1 = NO(IRP)
            K2 = NV(IRP)
            CALL XGEMM ('C','N',K1,N,K2,A1,BUF1(OFF1),K2,T1(OFF1),K2,
     &                   A0,BUF2,K1)
            CALL XGEMM ('N','N',M,N,K1,A1,T1(OFF1),M,BUF2,K1,
     &                   A1,S1(OFF1),M)
            OFF1 = OFF1 + M * N * RCW
         ENDDO
         CALL XTIME(0,2,'--- T1EQNS T*[HOV - F]*T      ')
        ENDIF
C ****************************************************
C ** END OF TASK SPLITTING. Parallel execution      **
C ****************************************************
C-----------------------------------------------------------------
C S(A,I) = S(A,I) + HOV(CK) * T2(A,C,I,K)
C-----------------------------------------------------------------
         CALL XTIME(0,1,'--- T1EQNS HOV*T2(A,C,I,K       ')
         CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,T2,BUF1)

         K    = IDIST(1,3,1)
         OFF1 = IDIST(2,3,1)
         OFF2 = IDIST(3,3,1)
         M    = MVO(1)
         CALL XGEMV ('N',M,K,A1,BUF1(OFF1),M,HOV(OFF2),1,A1,S1,1)
         CALL XTIME(0,2,'--- T1EQNS HOV*T2(A,C,I,K       ')
      ENDIF

C ********************
C ** TASK SPLITTING **
C ********************

      IF (MYPROC.EQ.MASTER) THEN

      CALL XTIME(0,1,'--- T1EQNS   HV*T / T*HO              ')
C-----------------------------------------------
C S(A,I) = S(A,I) + HV(A,C) * T(C,I)
C-----------------------------------------------
      CALL CNTRCT ('N','N',NV,NO,NV,A1,HV,T1,A1,S1,NREP)
C-----------------------------------------------
C S(A,I) = S(A,I) - T(A,K) * HO(K,I)
C-----------------------------------------------
      CALL CNTRCT ('N','N',NV,NO,NO,-A1,T1,HO,A1,S1,NREP)
      CALL XTIME(0,2,'--- T1EQNS   HV*T / T*HO            ')

C--------------------------------------------
C S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------
      CALL XTIME(0,1,'--- T1EQNS VOOO*TAU                ')
      IF (TEQ) THEN
         CALL GETTAU (T1,T2,BUF3)
      ELSE
         CALL XCOPY (NDIMT2,T2,1,BUF3,1)
      ENDIF
      CALL T1EQN1 (BUF1,BUF2,BUF3,NBUF3,S1)
      CALL XTIME(0,2,'--- T1EQNS VOOO*TAU               ')

C ****************************************************
C ** END OF TASK SPLITTING. THE COMMON STUFF STARTS **
C ****************************************************

      ENDIF

C--------------------------------------------
C S(A,I) = S(A,I) + V(A,K,C,D) * TAU(C,D,I,K)   ! VOVV contributions !
C--------------------------------------------
      IF (TEQ) THEN
         CALL GETTAU (T1,T2,BUF1)
      ELSE
C -------------------------------------------
C TAU'(CD,IK) = L(CD,IK) +  G(C,B) * T(BD,IK)   <<<<<<------- LAMBDA EQNS
C -------------------------------------------
         CALL GETAMPT (BUF1,BUF2)
         CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &                BUF2,BUF1)
         CALL CNTRCT ('N','N',NV,NVOOT,NV,A1,GVV,BUF1,A0,BUF2,NREP)
         CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &                BUF2,BUF1)
         CALL XAXPY (NDIMT2,A1,T2,1,BUF1,1)
      ENDIF
      CALL XTIME(0,1,'--- T1EQNS VOVV contribution  ')
      CALL T1EQN2 (BUF1,BUF2,BUF3,NBUF3,S1)   ! VOVV contractions, term 6 !
      CALL XTIME(0,2,'--- T1EQNS VOVV contribution  ')

C--------------------------------------
C S(A,I) = S(A,I) - W(A,K,C,I) * T(C,K)    ! VOVO contributions
C--------------------------------------
      CALL XTIME(0,1,'--- T1EQNS VOVO * T(C,K)         ')
      IF (TEQ) THEN
         CALL GETVOVO (BUF1)
         CALL DELINT ('VOVO','KDDK',BUF1,0,0,0)
         CALL SRT16 (NREP,MULTB,LTR,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &              BUF1,BUF2)
      ELSE
         CALL XCOPY (JVOVO(2),H,1,BUF2,1)
      ENDIF

      M    = MVO(1)
      K    = IDIST(1,3,1)
      OFF1 = IDIST(2,3,1)
      OFF2 = IDIST(3,3,1)
      CALL XGEMV ('N',M,K,-A1,BUF2(OFF1),M,T1(OFF2),1,A1,S1,1)
      CALL XTIME(0,2,'--- T1EQNS VOVO * T(C,K)        ')
C------------------------------------
C Terms specific for Lambda1 equations
C------------------------------------
      IF (.NOT.TEQ) THEN
C----------------------------------------------------
C S(A,I) = S(A,I) - T(A,L) * V(K,L,C,D) * L(C,D,K,I)
C        = S(A,I) - T(A,L) * V*(C,D,K,L) * L(C,D,K,I)
C----------------------------------------------------
      CALL GETVVOO(BUF2)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             BUF2,BUF1)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             T2,BUF2)
      CALL CNTRCT ('C','N',NO,NO,NVVOT,A1,BUF1,BUF2,A0,BUF3,NREP)
      CALL GETAMPT (BUF1,BUF2)
      CALL CNTRCT ('N','N',NV,NO,NO,-A1,BUF1,BUF3,A1,S1,NREP)
C----------------------------------------------------
C S(A,I) = S(A,I) - L(A,C,K,L) * V(K,L,D,C) * T(D,I)
C        = S(A,I) - L(A,C,K,L) * V*(D,C,K,L) * T(D,I)
C----------------------------------------------------
      CALL GETVVOO(BUF2)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF2,BUF1)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             T2,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,A1,BUF2,BUF1,A0,BUF3,NREP)
      CALL GETAMPT (BUF1,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP =1 , NREP
         M = NV(IRP)
         N = NO(IRP)
         K = NV(IRP)
         CALL XGEMM ('N','N',M,N,K,-A1,BUF3(OFF1),M,BUF1(OFF2),K,
     &               A1,S1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C--------------------------------------
C S(A,I) = S(A,I) - W(C,A,D,I) * G(C,D)
C--------------------------------------
      M = NFVV
      N = NFVO
      IF (NBUF3.LT.JVVVO(2)) THEN
         CALL QUIT(" OUT OF CORE STEP IN GV NOT IMPLEMENTED")
      ENDIF
      USEDZ = .TRUE.
      RIGHT = .FALSE.
clv : changed to SRT20D, untested as the rest of the lambda code
clv : addition should be supplied with a minus sign
      CALL SRT20D (NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &             JJVO,BUF2,NBUF2,S1,GVV,USEDZ,RIGHT)
C--------------------------------------
C S(A,I) = S(A,I) + W(A,K,L,I) * G(K,L)
C--------------------------------------
      CALL SRT1ST4 (NREP,MULTB,LTR,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &              CZ,BUF2)
      M = NFVO
      K = NFOO
      CALL XGEMV ('N',M,K,A1,BUF2,M,GOO,1,A1,S1,1)
C---------------------------------------------------
C S(A,I) = S(A,I) + G(C,B) * W(B,K,C,D) * T(A,D,K,I)
C---------------------------------------------------
      USEDZ = .FALSE.
      RIGHT = .FALSE.
      CALL SRT20D (NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &            JJVO,BUF2,NBUF2,BUF1,GVV,USEDZ,RIGHT)
      CALL GETAMPT (BUF3,BUF2)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,BUF2,BUF3)
      M = MVO(1)
      K = MVO(1)
      CALL XGEMV ('N',M,K,A1,BUF3,M,BUF1,1,A1,S1,1)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T2EQNS (EPS,T1,T2,GO,GV,HOV,H,AZ,CZ,GOO,GVV,
     &                   BUF1,BUF2,BUF3,NBUF3,S2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate T2 equations
C
C---------------Routines called----------------------------------------
C
C     BLAS routines :                 XGEMM, XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF3
      REAL*8 EPS(*)
      REAL*8 T1(*),T2(*)
      REAL*8 S2(*)
      REAL*8 GO(*),GV(*),HOV(*)
      REAL*8 GOO(*),GVV(*)
      REAL*8 AZ(*),H(*),CZ(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "eqns.inc"
#include "ihm.inc"
#include "symm.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ
      REAL*8 T1VAL(2)
      DATA T1VAL /2*0.0D0/
      complex*16 zdotc
      real*8 ddot
C
C---------------Executable code--------------------------------------
C
C
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
      ELSEIF (EQNS.EQ.'FOCKSP') THEN
         TEQ = .TRUE.
      ELSE
         CALL QUIT('Wrong keyword in call to T2EQNS')
      ENDIF

C---- Initialize to zero...
      CALL XCOPY (NDIMT2,A0,0,S2,1)
C---------------------------------------------
C S(AB,IJ) = S(AB,IJ) + TAU(AB,KL) * AZ(KL,IJ)
C---------------------------------------------
      CALL XTIME(0,1,'--- T2EQNS: TAU*AINTM contraction')
      IF (TEQ) THEN
         CALL GETTAU (T1,T2,BUF1)
      ELSE
         CALL XCOPY (NDIMT2,T2,1,BUF1,1)
      ENDIF

      OFF3 = 1
      DO IREP = 1, NREP
        M    = NVVT(IREP)
        N    = NOOT(IREP)
        KLD  = NOOT(IREP)
        K    = IDIST(1,5,IREP)
        OFF1 = IDIST(2,5,IREP)
        OFF2 = IDIST(3,5,IREP)
        CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,AZ(OFF2),KLD,
     &              A1,S2(OFF3),M)
        OFF3 = OFF3 + M * N * RCW
      ENDDO

C ***
C *** original CNTRCT call in serial version
C ***
C *** CALL CNTRCT ('N','N',NVVT,NOOT,NOOT,A1,BUF1,AZ,A1,S2,NREP)
C ***

      CALL XTIME(0,2,'--- T2EQNS: TAU*AINTM contraction')

C ** we will not parallelize a summation over just one occupied index
C ** if the system is small this can lead to zero batches. Only master
C ** does the T2(AB,IK) * GO(K,J) contraction

C ********************
C ** TASK SPLITTING **
C ********************

      IF(MYPROC.EQ.MASTER) THEN
C----------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) - PIJ [ T2(AB,IK) * GO(K,J) ]
C----------------------------------------------------------------
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             T2,BUF1)
      CALL CNTRCT ('N','N',NVVOT,NO,NO,A1,BUF1,GO,A0,BUF2,NREP)
      CALL SRT1T3 (NREP,MULTB,LTR,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             BUF2,BUF1)
      CALL XAXPY (NDIMT2,-A1,BUF1,1,S2,1)

C ***************************
C ** END OF TASK SPLITTING **
C ***************************

      ENDIF
C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PIJ [ W(AB,CJ) * T1(C,I) ]
C-----------------------------------------------------------------------
      CALL XTIME(0,1,'--- T2EQNS: VOVV*T1              ')
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         MINT = NBUF3/NVO(IRP)  ! max. # of NVO batches that fit in mem.
  10     IF (TEQ) THEN
            CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
            CALL DELINT ('VOVV','DKKK',BUF3,IRP,ISTART,NINT)
         ELSE
            CALL GETDZ   (IRP,ISTART,NINT,DONE,BUF3)
         ENDIF
         IF (CARITH) CALL CONJUGA (NVO(IRP)*NINT,BUF3,1)
         IJ = 0
         DO 20 JREP = 1, NREP
           CREP = MULTB(JREP,IRP+NREP,2)
           IREP = CREP
           IF (IREP.LT.JREP) GOTO 20
           DO J = 1, NO(JREP)
              IMIN = 1
              IF (JREP.EQ.IREP) IMIN = J + 1
              DO I = IMIN, NO(IREP)
                 IJ = IJ + 1
                 ABIJ = (IVVOOTT(IRP)+(IJ-1)*NVVT(IRP)+ISTART)*RCW+1
                 DO C = 1, NV(CREP)
                    CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = T1(T1CI)
                    IF (CARITH) T1VAL(2) = T1(T1CI+1)
                    CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABIJ),1)
                 ENDDO
              ENDDO
           ENDDO
 20      CONTINUE
         JI = 0
         DO 21 IREP = 1, NREP
           JREP = MULTB(IREP,IRP+NREP,2)
           CREP = IREP
           IF (JREP.LT.IREP) GOTO 21
           DO I = 1, NO(IREP)
              JMIN = 1
              IF (JREP.EQ.IREP) JMIN = I + 1
              DO J = JMIN, NO(JREP)
                 JI = JI + 1
                 ABJI = (IVVOOTT(IRP)+(JI-1)*NVVT(IRP)+ISTART)*RCW+1
                 DO C = 1, NV(CREP)
                    CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = - T1(T1CI)
                    IF (CARITH) T1VAL(2) = - T1(T1CI+1)
                    CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABJI),1)
                 ENDDO
              ENDDO
           ENDDO
 21      CONTINUE
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
 30    CONTINUE
      CALL XTIME(0,2,'--- T2EQNS: VOVV*T1              ')

C ********************
C ** TASK SPLITTING **
C ********************
C ** here a VOOV dependent contribution (Term 7 in the paper) is to be
C ** calculated only by the master

      IF(MYPROC.EQ.MASTER) THEN

C--------------------------------
C Terms specific for T2 equations
C--------------------------------
      IF (TEQ) THEN
C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ [T1(A,K)*W(BK,IC)*T1(C,J)]
C-----------------------------------------------------------------------
         CALL GETVOVO (BUF1)
         CALL DELINT ('VOVO','KDDK',BUF1,0,0,0)
         CALL SRT19 (NREP,MULTB,NVO,NV,NO,NVOO,KVOOV,KKVOO,
     &               BUF1,BUF2)
         CALL CNTRCT ('N','N',NVOO,NO,NV,A1,BUF2,T1,A0,BUF1,NREP)
         CALL SRT1T3 (NREP,MULTB,LTR,NVO,NO,NO,NVOO,KVOOO,KKVOO,
     &                BUF1,BUF3)
C-------------------------------------------------------------------
C        Integrals are now in VO,OOT format : combine with next term
C-------------------------------------------------------------------
      ENDIF

      ENDIF

C ****************************
C ** END OF TASK SPLITTING. **
C ****************************
C ** this contribution coming from the intermediate G(A,C) has to be
C ** calculated on each node. The fact that G(A,C) contains LOCAL +
C ** DISTRIBUTED parts as well is accounted for by the original calculation
C ** of the G(A,C) array before entering T2EQNS

C-------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB [ GV(A,C) * T2(CB,IJ)]
C-------------------------------------------------
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF2)
      CALL CNTRCT ('N','N',NV,NVOOT,NV,A1,GV,BUF2,A0,BUF1,NREP)

C ** the contraction of G(A,C) with the T2 amplitudes now is in BUF1 !
C ** these are the partial contributions available on EACH node


C ***
C ****** task splitting
C ***
      IF(MYPROC.EQ.MASTER) THEN

C------------------------------------------------
C S(AB,IJ) = S(AB,IJ) - PAB [ T1(A,K) *W(KB,IJ) ]
C------------------------------------------------
      IF (TEQ) THEN
         CALL GETVOOO (BUF2)
         CALL DELINT ('VOOO','KDKK',BUF2,0,0,0)
         CALL XAXPY (NV2,-A1,BUF3,1,BUF2,1)
      ELSE
         CALL XCOPY (NV2,CZ,1,BUF2,1)
      ENDIF
      CALL SRT1L1 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,IVOOOT,IIOV,BUF2,BUF3)
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NO,NV,NOOT,LOVOOT,LLVOOT,
     &             BUF3,BUF2)
      CALL CNTRCT ('N','N',NV,NVOOT,NO,-A1,T1,BUF2,A1,BUF1,NREP)

C ***
C *** end of task splitting
C ***
      ENDIF


      CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF1,BUF2)
      CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)


C ** now the final H*T2 contraction. H is complete since we communicated
C ** it already.

C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ [H(AI,CK)*T2(CK,BJ)]
C-----------------------------------------------------------------------

      CALL XTIME(0,1,'--- T2EQNS: HINTM*T2           ')
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              T2,BUF1)
      CALL XCOPY (J2VOVO(NREP+1),A0,0,BUF2,1)

      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M    = MVO(IREP)
         N    = MVO(JREP)
         KLD  = MVO(IREP)
         K    = IDIST(1,10,JREP)
         OFF1 = IDIST(2,10,JREP)
         OFF2 = IDIST(3,10,JREP)
         CALL XGEMM ('N','N',M,N,K,A1,H(OFF1),M,BUF1(OFF2),KLD,
     &                A0,BUF2(OFF3),M)
         OFF3 = OFF3 + M * N * RCW
      ENDDO

      CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              BUF2,BUF1)
      CALL XAXPY (NDIMT2,A1,BUF1,1,S2,1)
      CALL XTIME(0,2,'--- T2EQNS: HINTM*T2           ')



C==========================================================
C============= L A M B D A   E Q U A T I O N S ============
C==========================================================


      IF (.NOT.TEQ) THEN
C-------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + TAU(AB,KL) * V(KL,CD) * L(CD,IJ)
C-------------------------------------------------------
      CALL GETAMPT (BUF1,BUF2)
      CALL GETTAU (BUF1,BUF2,BUF3)
      CALL GETVVOO (BUF1)
      CALL CNTRCT ('C','N',NOOT,NOOT,NVVT,A1,BUF1,T2,A0,BUF2,NREP)
      CALL CNTRCT ('N','N',NOOT,NOOT,NVVT,A1,BUF3,BUF2,A1,S2,NREP)
C-----------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ. [T(A,I)*H(B,J)]
C-----------------------------------------------
      CALL LAMBD21 (T1,HOV,S2)
C---------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB. [G(A,C)*W(CB,IJ)]
C---------------------------------------------
      CALL GETVVOO (BUF1)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF1,BUF2)
      CALL CNTRCT ('N','N',NV,NVOOT,NV,A1,GVV,BUF2,A0,BUF1,NREP)
      CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF1,BUF2)
      CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)
C---------------------------------------------
C S(AB,IJ) = S(AB,IJ) - PIJ. [W(AB,IK)*G(K,J)]
C---------------------------------------------
      CALL GETVVOO(BUF2)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             BUF2,BUF1)
      CALL CNTRCT ('N','N',NVVOT,NO,NO,A1,BUF1,GOO,A0,BUF2,NREP)
      CALL SRT1T3 (NREP,MULTB,LTR,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             BUF2,BUF1)
      CALL XAXPY (NDIMT2,-A1,BUF1,1,S2,1)
C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ [L(D,I)*V(AK,CD)*T2(BC,KJ)]
C-----------------------------------------------------------------------
      CALL GETAMPT (BUF1,BUF3)
      CALL DINTM2 (T1,BUF1,BUF3,NBUF3,BUF2)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              BUF3,BUF1)
      CALL XCOPY (J2VOVO(NREP+1),A0,0,BUF3,1)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP =1 , NREP
         JRP = MULTB(IRP+NREP,1+NREP,2)
         M = MVO(IRP)
         N = MVO(JRP)
         K = MVO(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,BUF2(OFF1),M,BUF1(OFF2),K,
     &               A0,BUF3(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              BUF3,BUF1)
      CALL XAXPY (NDIMT2,A1,BUF1,1,S2,1)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T3CORR (FVO,EPS,T1,T2,T2S,T2T,BUF1,BUF2,BUF3,TAU3,
     &                   VVVO,VVOO,VOOO)
C
      use interface_to_mpi
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate non-iterative T3 corrections DE(+T), DE((T)), DE(-T)
C     Use formulation of Deegan and Knowles, CPL 227 (1994), 321.
C     Note that their W is called in the description because we use
C     W for the anti-symmetrized integrals.
C     This version keeps everything in core.
C
C---------------Routines called----------------------------------------
C
C     BLAS routines :                 XGEMM, XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVO(*),EPS(*)
      REAL*8 T1(*),T2(*),T2S(*),T2T(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*),TAU3(*)
      REAL*8 VVVO(*)
      REAL*8 VOOO(*),VVOO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "files.inc"
#include "symm.inc"
#include "inpt.inc"
#include "results.inc"
#include "complex.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
      DIMENSION ISTAT(df_MPI_STATUS_SIZE)
#endif

C
C---------------Local variables--------------------------------------
C
      REAL*8 EPSIJK,EIJKABC,VMIN(2)
      INTEGER*8 I8,IT3OFF1,JT3OFF1,KT3OFF1,EBCI,EBCJ,EBCK
C---------------Executable code--------------------------------------
C
C     Initialize the timing.
C
      CALL XTIME(0,1,'-- T3CORR: Integral resorting      ')
C
      ET1 = D0
      ET2 = D0
      ET3 = D0

      CALL GETVOOO(VOOO)
      CALL GETVVOO(VVOO)
C
C     Bring arrays in correct order
C
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,T2,T2S)
      CALL SRT22 (NREP,MULTB,NVVT,NO,NO,NVVOT,LOVVOT,LLOVVT,T2,T2T)

C call SRT23 with appropriate arguments. We use the TAU3 array as scratch space for the integral
C sort. This array is dimensioned like T2.

      CALL SRT23 (NREP,MULTB,NVO,NV,NO,NT3VVT,NT3VVVT,KT3VVVOT,
     &            LLT3VVVT,TAU3,VVVO,NDIMT2)
      IF (CARITH) THEN
C
C      Take complex conjugate. We can not use CONJUGA as the length of this
C      array may exceed the I4 adressing space.
C
         DO I8 = 2, 2*KT3VVVOT(NREP+1), 2
            VVVO(I8) = - VVVO(I8)
         ENDDO
      ENDIF
C
C     Compute TAU3 that is needed in evaluating the -T correction
C
      CALL GETTAU3 (T1,T2,TAU3)
C
      CALL XTIME(0,2,'-- T3CORR: Integral resorting      ')
C
      DO 300 IJKREP = 1, NREP
      DO 200 KREP = 1, NREP
      IJREP = MULTB(KREP,IJKREP,2)
      DO 100 JREP = KREP, NREP
      IREP = MULTB(JREP,IJREP+NREP,2)
      IF (IREP.LT.JREP) GOTO 100
      IKREP = MULTB(IREP,KREP,1)
      JKREP = MULTB(JREP,KREP,1)
      DO K = 1, NO(KREP)
      KK = IO(KREP) + K
      JMIN = 1
      IF (KREP.EQ.JREP) JMIN = K + 1
      DO J = JMIN, NO(JREP)
      JJ = IO(JREP) + J
      IMIN = 1
      IF (JREP.EQ.IREP) IMIN = J + 1
      DO I = IMIN, NO(IREP)
      II = IO(IREP) + I
      EPSIJK = - EPS(II) - EPS(JJ) - EPS(KK)
C
      IJ = IIOOT(IREP,JREP) + (J-1) * NO(IREP) + I
      IK = IIOOT(IREP,KREP) + (K-1) * NO(IREP) + I
      JK = IIOOT(JREP,KREP) + (K-1) * NO(JREP) + J
      IF (IREP.EQ.JREP) IJ = IJ - J*(J+1)/2
      IF (IREP.EQ.KREP) IK = IK - K*(K+1)/2
      IF (JREP.EQ.KREP) JK = JK - K*(K+1)/2
      IJOFF1 = IVVOOT(IJREP) + (IJ-1) * NVV(IJREP)
      IKOFF1 = IVVOOT(IKREP) + (IK-1) * NVV(IKREP)
      JKOFF1 = IVVOOT(JKREP) + (JK-1) * NVV(JKREP)
      IJOFF2 = IVOOOT(IJREP) + (IJ-1) * NVO(IJREP)
      IKOFF2 = IVOOOT(IKREP) + (IK-1) * NVO(IKREP)
      JKOFF2 = IVOOOT(JKREP) + (JK-1) * NVO(JKREP)
      IJOFF3 = IVVOOTT(IJREP) + (IJ-1) * NVVT(IJREP)
      IKOFF3 = IVVOOTT(IKREP) + (IK-1) * NVVT(IKREP)
      JKOFF3 = IVVOOTT(JKREP) + (JK-1) * NVVT(JKREP)

      IT3OFF1  = KT3VVVOT(IREP) + (I-1) * NT3VVVT(IREP)
      JT3OFF1  = KT3VVVOT(JREP) + (J-1) * NT3VVVT(JREP)
      KT3OFF1  = KT3VVVOT(KREP) + (K-1) * NT3VVVT(KREP)

      IOFF2  = LOVVOT(IREP) + (I-1) * NOVVT(IREP)
      JOFF2  = LOVVOT(JREP) + (J-1) * NOVVT(JREP)
      KOFF2  = LOVVOT(KREP) + (K-1) * NOVVT(KREP)
C
      CALL XTIME(0,1,'-- T3CORR: VOVV contraction        ')
C
      CALL XCOPY (IVVVT(NREP+1,IJKREP),A0,0,BUF1,1)
      OFF1 = 1
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      M = NV(AREP)
      N = NVVT(BCREP)
      NN = NT3VVT(BCREP)
      OFF2 = M * IDIST(1,2,BCREP) * RCW
C-------------------------------------------
C U'(ABC) = T(A,E:IJ) * W(E,BC:K)
C-------------------------------------------
      EREP = MULTB(KREP,BCREP+NREP,2)
      O = NV(EREP)
      AEIJ = (IJOFF1 + IIVV(AREP,EREP)) * RCW + 1
      EBCK = (KT3OFF1 + LLT3VVVT(EREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,NN,O,A1,T2S(AEIJ),M,
     &     VVVO(EBCK),O,A1,BUF1(OFF1+OFF2),M)
C-------------------------------------------
C U'(ABC) = U'(ABC) + T(A,E:JK) * W(E,BC:I)
C-------------------------------------------
      EREP = MULTB(IREP,BCREP+NREP,2)
      O = NV(EREP)
      AEJK = (JKOFF1 + IIVV(AREP,EREP)) * RCW + 1
      EBCI = (IT3OFF1 + LLT3VVVT(EREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,NN,O,A1,T2S(AEJK),M,
     &     VVVO(EBCI),O,A1,BUF1(OFF1+OFF2),M)
C-------------------------------------------
C U'(ABC) = U'(ABC) - T(A,E:IK) * W(E,BC:J)
C-------------------------------------------
      EREP = MULTB(JREP,BCREP+NREP,2)
      O = NV(EREP)
      AEIK = (IKOFF1 + IIVV(AREP,EREP)) * RCW + 1
      EBCJ = (JT3OFF1 + LLT3VVVT(EREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,NN,O,-A1,T2S(AEIK),M,
     &     VVVO(EBCJ),O,A1,BUF1(OFF1+OFF2),M)

C ** here the distributed VOVV integrals are used. we 
C ** communicate the partial U(ABC) and proceed on one 
C ** node only. BUF1, BUF2 and BUF3 have the same size so we
C ** can uns them as send/receive buffers here

#if defined (VAR_MPI)
      IF(NMPROC.GT.1) THEN
C       Switch off the overall timer for this part and focus on
C       the communication.
        CALL XTIME(0,2,'-- T3CORR: VOVV contraction        ')
        CALL XTIME(0,1,'-- T3CORR: addition of partial U   ')
        call interface_mpi_reduce_r1_work_f77(BUF1(OFF1),BUF2(1),
     &       RCW*M*N,
     &       op_MPI_SUM,MASTER,
     &       global_communicator)
        CALL XCOPY(M*N,BUF2(1),1,BUF1(OFF1),1)
C       Communication completed, switch back to operations timer
        CALL XTIME(0,2,'-- T3CORR: addition of partial U   ')
        CALL XTIME(0,1,'-- T3CORR: VOVV contraction        ')
      ENDIF
#endif

C *******************
C *** TASK SPLITTING
C *******************

      IF(MYPROC.EQ.MASTER) THEN

C-------------------------------------------
C U'(ABC) = U'(ABC) + W(A,L:IJ) * T(L,BC:K)
C-------------------------------------------
      LREP = MULTB(KREP,BCREP+NREP,2)
      O = NO(LREP)
      ALIJ = (IJOFF2 + IIVO(AREP,LREP)) * RCW + 1
      LBCK = (KOFF2 + LLOVVT(LREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,N,O,A1,VOOO(ALIJ),M,
     &            T2T(LBCK),O,A1,BUF1(OFF1),M)
C-------------------------------------------
C U'(ABC) = U'(ABC) + W(A,L:JK) * T(L,BC:I)
C-------------------------------------------
      LREP = MULTB(IREP,BCREP+NREP,2)
      O = NO(LREP)
      ALJK = (JKOFF2 + IIVO(AREP,LREP)) * RCW + 1
      LBCI = (IOFF2 + LLOVVT(LREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,N,O,A1,VOOO(ALJK),M,
     &            T2T(LBCI),O,A1,BUF1(OFF1),M)
C-------------------------------------------
C U'(ABC) = U'(ABC) - W(A,L:IK) * T(L,BC:J)
C-------------------------------------------
      LREP = MULTB(JREP,BCREP+NREP,2)
      O = NO(LREP)
      ALIK = (IKOFF2 + IIVO(AREP,LREP)) * RCW + 1
      LBCJ = (JOFF2 + LLOVVT(LREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,N,O,-A1,VOOO(ALIK),M,
     &            T2T(LBCJ),O,A1,BUF1(OFF1),M)
C
C ************************************
C **** end of task splitting 
C **** offsets have to be corrected
C ************************************
      ENDIF

      OFF1 = OFF1 + M * N * RCW
      ENDDO
C
      CALL XTIME(0,2,'-- T3CORR: VOVV contraction        ')
      CALL XTIME(0,1,'-- T3CORR: energy calculation      ')
C
C ** the U array is complete on the master node. 
C ** we can permute the indices and compute on the
C ** master only...


      IF(MYPROC.EQ.MASTER) THEN

C-------------------------------------
C U(ABC) = U'(ABC) - U'(BAC) + U'(CAB)
C-------------------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 10 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 10
      ABREP = MULTB(AREP,BREP,1)
      ACREP = MULTB(AREP,CREP,1)
      DO C = 1, NV(CREP)
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      AB = IIVVT(AREP,BREP) + (B-1) * NV(AREP) + AMIN
      AC = IIVVT(AREP,CREP) + (C-1) * NV(AREP) + AMIN
      IF (AREP.EQ.BREP) AB = AB - B*(B+1)/2
      IF (AREP.EQ.CREP) AC = AC - C*(C+1)/2
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC = (IVVVT(BCREP,IJKREP) + (BC-1)*NV(AREP) + AMIN-1) * RCW + 1
      BAC = (IVVVT(ACREP,IJKREP) + (AC-1)*NV(BREP) + B - 1) * RCW + 1
      CAB = (IVVVT(ABREP,IJKREP) + (AB-1)*NV(CREP) + C - 1) * RCW + 1
      DO A = AMIN, NV(AREP)
       BUF3(ABC) = BUF1(ABC) - BUF1(BAC) + BUF1(CAB)
       IF (CARITH) BUF3(ABC+1) = BUF1(ABC+1) - BUF1(BAC+1) + BUF1(CAB+1)
       ABC = ABC + RCW
       BAC = BAC + NV(BREP) * RCW
       CAB = CAB + NV(CREP) * RCW
      ENDDO
      ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C-------------------------------------------------
C Y2'(ABC) = V(AB:IJ) * T(C:K) + T(AB:IJ) * F(C:K)
C          + V(AB:JK) * T(C:I) + T(AB:JK) * F(C:I)
C          - V(AB:IK) * T(C:J) - T(AB:IK) * F(C:J)
C-------------------------------------------------
      CALL XCOPY (IVVVT(NREP+1,IJKREP),A0,0,BUF2,1)
      ABIJ = IJOFF3 * RCW + 1
      ABJK = JKOFF3 * RCW + 1
      ABIK = IKOFF3 * RCW + 1
      DO C = 1, NV(KREP)
         CK = (IVO(KREP) + (K-1) * NV(KREP) + C - 1) * RCW + 1
         ABC = (JVVVT(KREP,IJKREP)+(C-1)*NVVT(IJREP)) * RCW + 1
         CALL XAXPY (NVVT(IJREP),T1 (CK),VVOO(ABIJ),1,BUF2(ABC),1)
         CALL XAXPY (NVVT(IJREP),FVO(CK),T2  (ABIJ),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(IREP)
         CI = (IVO(IREP) + (I-1) * NV(IREP) + C - 1) * RCW + 1
         ABC = (JVVVT(IREP,IJKREP)+(C-1)*NVVT(JKREP)) * RCW + 1
         CALL XAXPY (NVVT(JKREP),T1 (CI),VVOO(ABJK),1,BUF2(ABC),1)
         CALL XAXPY (NVVT(JKREP),FVO(CI),T2  (ABJK),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(JREP)
         CJ = (IVO(JREP) + (J-1) * NV(JREP) + C - 1) * RCW + 1
         ABC = (JVVVT(JREP,IJKREP)+(C-1)*NVVT(IKREP)) * RCW + 1
         VMIN(1) = - T1(CJ)
         IF (CARITH) VMIN(2) = - T1(CJ+1)
         CALL XAXPY (NVVT(IKREP),VMIN,VVOO(ABIK),1,BUF2(ABC),1)
         VMIN(1) = - FVO(CJ)
         IF (CARITH) VMIN(2) = - FVO(CJ+1)
         CALL XAXPY (NVVT(IKREP),VMIN,T2  (ABIK),1,BUF2(ABC),1)
      ENDDO
C-----------------------------------------
C Y2(ABC) = Y2'(ABC) + Y2'(BCA) - Y2'(ACB)
C-----------------------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 20 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 20
      ABREP = MULTB(AREP,BREP,1)
      ACREP = MULTB(AREP,CREP,1)
      DO C = 1, NV(CREP)
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      AB = IIVVT(AREP,BREP) + (B-1) * NV(AREP) + AMIN
      AC = IIVVT(AREP,CREP) + (C-1) * NV(AREP) + AMIN
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      IF (AREP.EQ.BREP) AB = AB - B*(B+1)/2
      IF (AREP.EQ.CREP) AC = AC - C*(C+1)/2
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC1 = (IVVVT(BCREP,IJKREP) + (BC-1)*NV(AREP)+AMIN-1)  * RCW + 1
      ABC2 = (JVVVT(CREP,IJKREP) + (C-1)*NVVT(ABREP)+AB-1)   * RCW + 1
      BCA  = (JVVVT(AREP,IJKREP) + (AMIN-1)*NVVT(BCREP)+BC-1)* RCW + 1
      ACB  = (JVVVT(BREP,IJKREP) + (B-1)*NVVT(ACREP)+AC-1)   * RCW + 1
      DO A = AMIN, NV(AREP)
       BUF1(ABC1) = BUF2(ABC2) + BUF2(BCA) - BUF2(ACB)
       IF (CARITH) BUF1(ABC1+1) = BUF2(ABC2+1)+BUF2(BCA+1)-BUF2(ACB+1)
       ABC1 = ABC1 + RCW
       ABC2 = ABC2 + RCW
       ACB = ACB + RCW
       BCA = BCA + NVVT(BCREP) * RCW
      ENDDO
      ENDDO
      ENDDO
 20   CONTINUE
      ENDDO
C-----------------------------------------
C DE(+T)   = - U(ABC) * U(ABC) / EPS(ABC)
C DE'((T)) = - U(ABC) * Y2(ABC) / EPS(ABC)
C-----------------------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 30 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 30
      DO C = 1, NV(CREP)
      CC = IO(NREP+1) + IV(CREP) + C
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      BB = IO(NREP+1) + IV(BREP) + B
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC = (IVVVT(BCREP,IJKREP)+(BC-1)*NV(AREP)+AMIN-1) * RCW + 1
      DO A = AMIN, NV(AREP)
       AA = IO(NREP+1) + IV(AREP) + A
       EIJKABC = EPSIJK + EPS(AA) + EPS(BB) + EPS(CC)
       ET1 = ET1 - BUF3(ABC) * BUF3(ABC) / EIJKABC
       ET2 = ET2 - BUF3(ABC) * BUF1(ABC) / EIJKABC
       IF (CARITH) THEN
         ET1 = ET1 - BUF3(ABC+1) * BUF3(ABC+1) / EIJKABC
         ET2 = ET2 - BUF3(ABC+1) * BUF1(ABC+1) / EIJKABC
       ENDIF
       ABC = ABC + RCW
      ENDDO
      ENDDO
      ENDDO
 30   CONTINUE
      ENDDO
C--------------------------------
C Y3'(ABC) = TAU3(AB:IJ) * T(C:K)
C          + TAU3(AB:JK) * T(C:I)
C          - TAU3(AB:IK) * T(C:J)
C--------------------------------
      CALL XCOPY (IVVVT(NREP+1,IJKREP),A0,0,BUF2,1)
      ABIJ = IJOFF3 * RCW + 1
      ABJK = JKOFF3 * RCW + 1
      ABIK = IKOFF3 * RCW + 1
      DO C = 1, NV(KREP)
         CK = (IVO(KREP) + (K-1) * NV(KREP) + C - 1) * RCW + 1
         ABC = (JVVVT(KREP,IJKREP)+(C-1)*NVVT(IJREP)) * RCW + 1
         CALL XAXPY (NVVT(IJREP),T1 (CK),TAU3(ABIJ),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(IREP)
         CI = (IVO(IREP) + (I-1) * NV(IREP) + C - 1) * RCW + 1
         ABC = (JVVVT(IREP,IJKREP)+(C-1)*NVVT(JKREP)) * RCW + 1
         CALL XAXPY (NVVT(JKREP),T1 (CI),TAU3(ABJK),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(JREP)
         CJ = (IVO(JREP) + (J-1) * NV(JREP) + C - 1) * RCW + 1
         ABC = (JVVVT(JREP,IJKREP)+(C-1)*NVVT(IKREP)) * RCW + 1
         VMIN(1) = - T1(CJ)
         IF (CARITH) VMIN(2) = - T1(CJ+1)
         CALL XAXPY (NVVT(IKREP),VMIN,TAU3(ABIK),1,BUF2(ABC),1)
      ENDDO
C-----------------------------------------
C Y3(ABC) = Y3'(ABC) + Y3'(BCA) - Y3'(ACB)
C-----------------------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 40 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 40
      ABREP = MULTB(AREP,BREP,1)
      ACREP = MULTB(AREP,CREP,1)
      DO C = 1, NV(CREP)
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      AB = IIVVT(AREP,BREP) + (B-1) * NV(AREP) + AMIN
      AC = IIVVT(AREP,CREP) + (C-1) * NV(AREP) + AMIN
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      IF (AREP.EQ.BREP) AB = AB - B*(B+1)/2
      IF (AREP.EQ.CREP) AC = AC - C*(C+1)/2
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC1 = (IVVVT(BCREP,IJKREP) + (BC-1)*NV(AREP)+AMIN-1)  * RCW + 1
      ABC2 = (JVVVT(CREP,IJKREP) + (C-1)*NVVT(ABREP)+AB-1)   * RCW + 1
      BCA  = (JVVVT(AREP,IJKREP) + (AMIN-1)*NVVT(BCREP)+BC-1)* RCW + 1
      ACB  = (JVVVT(BREP,IJKREP) + (B-1)*NVVT(ACREP)+AC-1)   * RCW + 1
      DO A = AMIN, NV(AREP)
       BUF1(ABC1) = BUF2(ABC2) + BUF2(BCA) - BUF2(ACB)
       IF (CARITH) BUF1(ABC1+1) = BUF2(ABC2+1)+BUF2(BCA+1)-BUF2(ACB+1)
       ABC1 = ABC1 + RCW
       ABC2 = ABC2 + RCW
       ACB = ACB + RCW
       BCA = BCA + NVVT(BCREP) * RCW
      ENDDO
      ENDDO
      ENDDO
 40   CONTINUE
      ENDDO
C-----------------------------
C DE'(-T) = + U(ABC) * Y3(ABC)
C-----------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 50 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 50
      DO C = 1, NV(CREP)
      CC = IO(NREP+1) + IV(CREP) + C
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      BB = IO(NREP+1) + IV(BREP) + B
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC = (IVVVT(BCREP,IJKREP)+(BC-1)*NV(AREP)+AMIN-1) * RCW + 1
      DO A = AMIN, NV(AREP)
       AA = IO(NREP+1) + IV(AREP) + A
       ET3 = ET3 + BUF3(ABC) * BUF1(ABC) 
       IF (CARITH) THEN
         ET3 = ET3 + BUF3(ABC+1) * BUF1(ABC+1)
       ENDIF
       ABC = ABC + RCW
      ENDDO
      ENDDO
      ENDDO
 50   CONTINUE
      ENDDO
C

C ******** end of task splitting *******
      ENDIF
C ******** end of task splitting *******

      CALL XTIME(0,2,'-- T3CORR: energy calculation      ')
      ENDDO
      ENDDO
      ENDDO
  100 CONTINUE
  200 CONTINUE
  300 CONTINUE
C
      WRITE (IW,1071)
      WRITE (IW,1072) ESCF
      WRITE (IW,1073) ECCSD
      WRITE (IW,1074) ET1
      WRITE (IW,1075) ET2
      WRITE (IW,1076) ET3
      WRITE (IW,1077) ESCF+ECCSD+ET1
      WRITE (IW,1078) ESCF+ECCSD+ET1+ET2
      WRITE (IW,1079) ESCF+ECCSD+ET1+ET3
      ETOT = ESCF+ECCSD+ET1+ET2
C
 1071 FORMAT(//,'  Perturbative treatment of triple excitations',/)
 1072 FORMAT(' SCF energy :',T40,F25.15)
 1073 FORMAT(' CCSD correlation energy :',T40,F25.15)
 1074 FORMAT(' 4th order triples correction :',T40,F25.15)
 1075 FORMAT(' 5th order triples (T) correction :',T40,F25.15)
 1076 FORMAT(' 5th order triples -T  correction :',T40,F25.15)
 1077 FORMAT(' Total CCSD+T  energy :',T40,F25.15)
 1078 FORMAT(' Total CCSD(T) energy :',T40,F25.15)
 1079 FORMAT(' Total CCSD-T  energy :',T40,F25.15)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZCORE(T1,T2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Zeroes out core contributions (frozen spinors)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "freeze.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      OFF = 1
      DO IRP = 1, NREP
         DO I = 1, NO(IRP)
            IF (IFROZ(I,IRP).NE.0)
     &      CALL XCOPY (NV(IRP),A0,0,T1(OFF),1)
            OFF = OFF + NV(IRP) * RCW
         ENDDO
      ENDDO
C
      OFF = 1
      DO IJRP = 1, NREP
      N = NVVT(IJRP)
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            IF (IFROZ(I,IRP).NE.0.OR.IFROZ(J,JRP).NE.0)
     &         CALL XCOPY (N,A0,0,T2(OFF),1)
            OFF = OFF + N * RCW
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PARTS(IMODE,RMS,T1,T2,BUF1,BUF2)
C
      use interface_to_mpi
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     IMODE : 
C     1) Sums partial contributions to T1 and T2
C     2) Synchronizes T1 and T2
C
C     BUF1 and BUF2 have the length of T1 and T2 respectively and are
C     not used with imode 2. RMS is also updated because it controls
C     the iterative process
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*),BUF1(*),BUF2(*)
      INTEGER IMODE
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
      DIMENSION ISTAT(df_MPI_STATUS_SIZE)
#endif
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
#if defined (VAR_MPI)
      IF (NMPROC.EQ.1) RETURN
      IF (IMODE.EQ.1) THEN
          call interface_mpi_reduce_r1_work_f77(T1,BUF1,NDIMT1*RCW,
     &                    op_MPI_SUM,MASTER,global_communicator)
          call interface_mpi_reduce_r1_work_f77(T2,BUF2,NDIMT2*RCW,
     &                    op_MPI_SUM,MASTER,global_communicator)
          IF (MYPROC.EQ.MASTER) CALL XCOPY (NDIMT1,BUF1,1,T1,1)
          IF (MYPROC.EQ.MASTER) CALL XCOPY (NDIMT2,BUF2,1,T2,1)
      ENDIF
C
      call interface_mpi_bcast_r1_work_f77(T1,NDIMT1*RCW,MASTER,
     &               global_communicator)
      call interface_mpi_bcast_r1_work_f77(T2,NDIMT2*RCW,MASTER,
     &               global_communicator)
      call interface_mpi_BCAST(RMS,1,MASTER,
     &               global_communicator)
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DIISX (NDIMT1,NDIMT2,T1O,T1N,T2O,T2N,IT,
     &                  CC,BB,ITC,BB2,RMS)
C    
      IMPLICIT INTEGER (A-Z)
C     
C---------------Description--------------------------------------------
C
C     DIIS extrapolation of the t1 and t2 amplitudes
C     
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher (based on a routine by T.J. Lee)
C     
C---------------Calling variables--------------------------------------
C
      REAL*8 T2O(*),T2N(*),T1O(*),T1N(*)
      REAL*8 CC(MAXDIM),BB(MAXDIM+1,MAXDIM+2),BB2(MAXDIM+1,MAXDIM+2)
      REAL*8 RMS
C     
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "diis.inc"
#include "files.inc"
#include "waio.h"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C 
      REAL*8 XM,XMAX,XEN,XFAC,XADD,DET
      REAL*8 SUM1,SUM2,SSUM,DSUM,SDOT,DDOT
      COMPLEX*16 CFAC
C     
C---------------Executable code--------------------------------------
C
      IF (MYPROC.NE.MASTER) RETURN
C     
C     Handle special cases of initialization or finalization first
C
      IF (IT.EQ.-1) THEN
C
C       Initialization of the DIIS procedure;
C       compute record lengths, open files, write first set of amplitudes
C       This information is stored in the common block diis.inc
C
        NST1 = INTOWP(RCW*NDIMT1)/ NWORDS_ON_RECORD
        IF(NST1*NWORDS_ON_RECORD.NE.INTOWP(RCW*NDIMT1)) NST1 = NST1 + 1
        NST2 = INTOWP(RCW*NDIMT2)/ NWORDS_ON_RECORD
        IF(NST2*NWORDS_ON_RECORD.NE.INTOWP(RCW*NDIMT2)) NST2 = NST2 + 1
        CALL WAIO_OPEN(ITAMPL)
        CALL WAIO_OPEN(ITAMPD)
        SECT1 = MAXDIM*(NST1+NST2) + 1
        SECT2 = SECT1 + NST1
C       Place first input amplitudes at record maxdim+1 of ITAMPL
        CALL RWRIT(ITAMPL,T1N,INTOWP(NDIMT1*RCW),SECT1)
        CALL RWRIT(ITAMPL,T2N,INTOWP(NDIMT2*RCW),SECT2)
        IT = 0
        ITC = 0
        RMS = 1.0
C       Return to caller
        RETURN
      ENDIF
C
      IF (IT.EQ.-2) THEN
C
C        Stop DIIS and clean up files
C
        JCODE = 4
        CALL WAIO_CLOSE(ITAMPL,JCODE)
        CALL WAIO_CLOSE(ITAMPD,JCODE)
        RETURN
      ENDIF
C
C     Normal procedure : start by reading the old amplitudes that were used to produce T1N and T2N
C
      SECT1 = MAXDIM*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RREAD(ITAMPL,T1O,INTOWP(NDIMT1*RCW),SECT1)
      CALL RREAD(ITAMPL,T2O,INTOWP(NDIMT2*RCW),SECT2)
C
      IT = IT + 1
      IF(IT.GT.MAXDIM) IT = MAXDIM
      ITC = ITC + 1
      IF(ITC.GT.MAXDIM) ITC = 1
C
C     Calculate the difference between the input and output amplitudes
C     This is our error vector. Below we will consider the real and 
C     imaginary differences as independent so that we can work with
C     a real error vector that is twice as long in complex as in real
C     cases.
C
      CALL XAXPY (NDIMT1,-A1,T1N,1,T1O,1)
      CALL XAXPY (NDIMT2,-A1,T2N,1,T2O,1)
      CALL XSCAL (NDIMT1,-A1,T1O,1)
      CALL XSCAL (NDIMT2,-A1,T2O,1)
C
      IF (CARITH) THEN
         SUM1 = DSUM (RCW*NDIMT1,T1O,1)
         SUM2 = DDOT (RCW*NDIMT1,T1O,1,T1O,1)
         SUM1 = SUM1 + DSUM (RCW*NDIMT2,T2O,1)
         SUM2 = SUM2 + DDOT (RCW*NDIMT2,T2O,1,T2O,1)
         RMS = SUM1 * SUM1
      ELSE
         SUM1 = DSUM (NDIMT1,T1O,1)
         SUM2 = DDOT (NDIMT1,T1O,1,T1O,1)
         SUM1 = SUM1 + DSUM (NDIMT2,T2O,1)
         SUM2 = SUM2 + DDOT (NDIMT2,T2O,1,T2O,1)
         RMS = SUM1 * SUM1
      ENDIF
C
      RMS = SQRT(ABS(SUM2-RMS))
C
      SECT1 = (ITC-1)*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RWRIT(ITAMPD,T1O,INTOWP(NDIMT1*RCW),SECT1)
      CALL RWRIT(ITAMPL,T1N,INTOWP(NDIMT1*RCW),SECT1)
      CALL RWRIT(ITAMPD,T2O,INTOWP(NDIMT2*RCW),SECT2)
      CALL RWRIT(ITAMPL,T2N,INTOWP(NDIMT2*RCW),SECT2)
C
      DO 30 N = 1,IT
C
      SECT1 = (N-1)*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RREAD(ITAMPD,T1N,INTOWP(NDIMT1*RCW),SECT1)
      CALL RREAD(ITAMPD,T2N,INTOWP(NDIMT2*RCW),SECT2)
C
      IF (CARITH) THEN
         SUM2 = DDOT (RCW*NDIMT1,T1N,1,T1O,1)
         SUM2 = SUM2 + DDOT (RCW*NDIMT2,T2N,1,T2O,1)
      ELSE
         SUM2 = DDOT (NDIMT1,T1N,1,T1O,1)
         SUM2 = SUM2 + DDOT (NDIMT2,T2N,1,T2O,1)
      ENDIF
C
      BB2(N,ITC) = SUM2
      BB2(ITC,N) = BB2(N,ITC)
C
   30 CONTINUE
C
C  TRANSFER INTO THE BB ARRAY
C
      DO 60 N = 1,IT
      DO 70 M = 1,N-1
      BB(N,M) = BB2(N,M)
      BB(M,N) = BB2(M,N)
  70  CONTINUE
      BB(N,N) = BB2(N,N)
  60  CONTINUE
C
C  FIND THE MAXIMUM AND SCALE
C
      XM =  ABS(BB(1,1))
      DO 80 N = 1,IT
      DO 90 M = 1,N
      XFAC =  ABS(BB(N,M))
      XMAX =  DMAX1(XM,XFAC)
      XM = XMAX
   90 CONTINUE
   80 CONTINUE
C
      XM = AR1/XM
      DO 100 N = 1,IT
      DO 110 M = 1,IT
      BB(N,M) = BB(N,M)*XM
  110 CONTINUE
  100 CONTINUE
C
      IT1 = IT+1
      IT2 = IT+2
      DO 120 N = 1,IT
      BB(N,IT1) = -AR1
      BB(IT1,N) = -AR1
  120 CONTINUE
C
      BB(IT1,IT1) = AR0
      DO 130 N = 1,IT
      BB(N,IT2) = AR0
  130 CONTINUE
      BB(IT1,IT2) = -AR1
C
      CALL FLIN(BB,MAXDIM+1,IT1,1,DET)
      XADD = AR0
      DO 140 N = 1,IT
      CC(N) = BB(N,IT2)
      XADD = XADD + CC(N)
  140 CONTINUE
C
      CALL XCOPY(NDIMT1,A0,0,T1N,1)
      CALL XCOPY(NDIMT2,A0,0,T2N,1)
C
      DO 150 N = 1,IT
C
      SECT1 = (N-1)*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RREAD(ITAMPL,T1O,INTOWP(NDIMT1*RCW),SECT1)
      CALL RREAD(ITAMPL,T2O,INTOWP(NDIMT2*RCW),SECT2)
C
      CFAC = DCMPLX(CC(N),AR0)
      CALL XAXPY (NDIMT1,CFAC,T1O,1,T1N,1)
      CALL XAXPY (NDIMT2,CFAC,T2O,1,T2N,1)
C
  150 CONTINUE
C
C     Write extrapolated amplitudes (which become input in the next iteration)
C     at the end of the amplitude file.
C
      SECT1 = MAXDIM*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RWRIT(ITAMPL,T1N,INTOWP(NDIMT1*RCW),SECT1)
      CALL RWRIT(ITAMPL,T2N,INTOWP(NDIMT2*RCW),SECT2)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCENG(T1,T2,FVO,BUF1,BUF2,ECC,ECCIM)        
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates CC energy. Checks result (should be real).
C
C---------------Routines called----------------------------------------
C
C     BLAS routines
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C              Miro Ilias,July 2007 - added imaginary part of CCSD energy
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
      REAL*8 FVO(*)
      REAL*8 BUF1(*),BUF2(*)
      REAL*8 ECC,ECCIM
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "files.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CECC1(2),CECC2(2)
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NDIMT1,FVO,1,BUF1,1)
      CALL DELFCK ('VO','DD',BUF1)
C
      CALL XDOTC(CECC1,NFVO,BUF1,1,T1,1)
      ECC = CECC1(1)
      IF (CARITH) THEN
         IF (ABS(CECC1(2)).GT.ACCUR) 
     &      WRITE (IW,1000) 'T1',CECC1(2)
      ENDIF
C
c     print*,'T1 contribution to CC energy',ecc
C
      CALL GETVVOO (BUF1)
      CALL DELINT ('VVOO','DDDD',BUF1,0,0,0)
      CALL GETTAU  (T1,T2,BUF2)
C
      CALL XDOTC(CECC2,NDIMT2,BUF1,1,BUF2,1)
      ECC = ECC + CECC2(1)
      IF (CARITH) THEN
         IF (ABS(CECC2(2)).GT.ACCUR) 
     &      WRITE (IW,1000) 'T2',CECC2(2)
      ENDIF

CMI  ... sum imaginary energies
      IF (CARITH) THEN
        ECCIM = CECC1(2) + CECC2(2)
      ENDIF
C
      RETURN
1000  FORMAT (/' WARNING : Imaginary contribution to ',A2,
     &' part of CC energy :',F20.15)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MP2ENG(T1,T2,FVO,BUF1,EMP2)        
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates MP2 energy. Checks result (should be real).
C
C---------------Routines called----------------------------------------
C
C     BLAS routines
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
      REAL*8 FVO(*)
      REAL*8 BUF1(*)
      REAL*8 EMP2
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "files.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CEMP21(2),CEMP22(2)
C
C---------------Executable code--------------------------------------
C
      CALL XDOTC(CEMP21,NFVO,FVO,1,T1,1)
      EMP2 = CEMP21(1)
      IF (CARITH) THEN
         IF (ABS(CEMP21(2)).GT.ACCUR)    
     &      WRITE (IW,1000) 'T1',CEMP21(2)
      ENDIF
C
c     print*,'T1 contribution to MP2 energy',emp2
C
      CALL GETVVOO(BUF1)
      CALL XDOTC(CEMP22,NDIMT2,BUF1,1,T2,1)
      EMP2 = EMP2 + CEMP22(1)
      IF (CARITH) THEN
         IF (ABS(CEMP22(2)).GT.ACCUR)    
     &      WRITE (IW,1000) 'T2',CEMP22(2)
      ENDIF
C
      RETURN
1000  FORMAT (/' WARNING : Imaginary contribution to ',A2,
     &' part of MP2 energy :',F20.15)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOM (EPS,T1,T2,CT1,CT2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Divide T1 and T2 by denominators
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 T1(*),T2(*)
      COMPLEX*16 CT1(*),CT2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC,FAC1,FAC2,FAC3,ATVAL,ABS,DHOLU
C
C---------------Executable code--------------------------------------
C
      DHOLU = 1.D-6
      II = 0
      AI = 0
      DO IRP = 1, NREP
      DO I = 1, NO(IRP)
         II = II + 1
         FAC1 = EPS(II)
         AOFF = IO(NREP+1) + IV(IRP)
         DO A = 1, NV(IRP)
             AA = AOFF + A
             FAC = FAC1 - EPS(AA)
             AI = AI + 1
C
C            Check for small denominators. All denominators should be 
C            larger than the minimum HOMO-LUMO gap DHOLU to avoid 
C            convergence problems. An exception is the case in which
C            the numerator is zero (e.g. by symmetry), then we allow
C            execution to be continued. 
C
             IF (-FAC.GT.DHOLU) THEN
                IF (CARITH) THEN
                   CT1(AI) = CT1(AI)/FAC
                ELSE
                   T1(AI) = T1(AI)/FAC
                ENDIF
             ELSE
                IF (CARITH) THEN
                   ATVAL = ABS(CT1(AI))
                ELSE
                   ATVAL = ABS(T1(AI))
                ENDIF
                IF (ATVAL.GT.DHOLU) THEN
                   CALL QUIT ('Zero or negative HOMO-LUMO gap')
                ELSE
                  IF (CARITH) THEN
                     CT1(AI) = A0
                  ELSE
                     T1(AI) = AR0
                  ENDIF
                ENDIF
             ENDIF
         ENDDO
      ENDDO
      ENDDO
C
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP) 
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 - EPS(AA)
                  ABIJ = ABIJ + 1
C                 Check for small denominators, see comment above.
                  IF (-FAC.GT.DHOLU) THEN
                     IF (CARITH) THEN
                        CT2(ABIJ) = CT2(ABIJ)/FAC
                     ELSE
                        T2(ABIJ) = T2(ABIJ)/FAC
                     ENDIF
                  ELSE
                     IF (CARITH) THEN
                        ATVAL = ABS(CT2(ABIJ))
                     ELSE
                        ATVAL = ABS(T2(ABIJ))
                     ENDIF
                     IF (ATVAL.GT.DHOLU) THEN
                        CALL QUIT ('Zero or negative HOMO-LUMO gap')
                     ELSE
                       IF (CARITH) THEN
                          CT2(ABIJ) = A0
                       ELSE
                          T2(ABIJ) = AR0
                       ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T1DIAG(DIAG,T1,BUF1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculate T1-diagnostic
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),BUF1(*),DIAG
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER SUM
      REAL*8 RSUM,DDOT,SDOT
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NDIMT1,T1,1,BUF1,1)
      CALL DELFCK ('VO','DD',BUF1)
C
      SUM = 0
      DO I = 1, NREP
         SUM = SUM + NE(I)
      ENDDO
      RSUM = DBLE(SUM)
C
      DIAG = DDOT (RCW*NDIMT1,BUF1,1,BUF1,1)
      DIAG =  SQRT(DIAG/RSUM)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T1EQN1 (BUF1,BUF2,BUF3,NBUF3,S1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate part of T1 or Lambda1 equations
C     The same contraction is used in the construction of the
C     X-intermediate.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF3
      REAL*8 S1(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
C--------------------------------------------
C S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------

C ** TAU now reorganised in BUF1: NV*NVOOT

      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF3,BUF1)
      CALL GETVOOO (BUF3)
      CALL DELINT ('VOOO','DKDD',BUF3,0,0,0)

C     --------------------------------
C     Sort from 12,34 to 134,2 (VOO,O)
C     --------------------------------
      CALL SRT6 (NREP,MULTB,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           BUF3,BUF2)

C ** VOOO integrals now reorganised in BUF2: NVOOT * NO

C     ------------------------
C     V(K,L,C,I) = V*(C,I,K,L)
C     ------------------------
      IF (CARITH) CALL CONJUGA (LOVOOT(NREP+1),BUF2,1)

C     OFF3 = 1
C     DO IREP = 1, NREP
C       M    = NV(IREP)
C       N    = NO(IREP)
C       KLD  = NVOOT(IREP)
C       K    = IDIST(1,6,IREP)
C       OFF1 = IDIST(2,6,IREP)
C       OFF2 = IDIST(3,6,IREP)
C       CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,BUF2(OFF2),KLD,
C    &              A1,S1(OFF3),M)
C       OFF3 = OFF3 + M * N * RCW
C     ENDDO

C***
C*** original CNTRCT call in serial version
C*** since parallel code (above) is not debugged
C***
C
      CALL CNTRCT ('N','N',NV,NO,NVOOT,A1,BUF1,BUF2,A1,S1,NREP)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T1EQN2 (BUF1,BUF2,BUF3,NBUF3,S1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate part of T1 or Lambda1 equations
C     The same contraction is used in the construction of the
C     X-intermediate.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF3
      REAL*8 S1(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      REAL*8 SDOT,DDOT
C
C---------------Executable code--------------------------------------
C
C--------------------------------------------
C S(A,I) = S(A,I) + V(A,K,C,D) * TAU(C,D,I,K)
C--------------------------------------------
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         NINT = NBUF3/NVO(IRP)  ! maximum # of NVO batches

C ** NINT is updated on exit of GETVOVV !

  10     CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,NINT)
         CALL DELINT ('VOVV','KDDD',BUF3,IRP,ISTART,NINT)
         IK = 0
         DO 20 KREP = 1, NREP
           IREP = MULTB(KREP,IRP+NREP,2)
           IF (KREP.GT.IREP) GOTO 20
           DO K = 1, NO(KREP)
              IMIN = 1
              IF (KREP.EQ.IREP) IMIN = K + 1
              DO I = IMIN, NO(IREP)
                 IK = IK + 1
                 CDIK = (IVVOOTT(IRP)+(IK-1)*NVVT(IRP)+ISTART)*RCW+1
                 DO A = 1, NV(IREP)
                    AK = (IIVO(IREP,KREP)+(K-1)*NV(IREP)+A-1) * RCW + 1
                    T1AI = (IVO(IREP)+(I-1)*NV(IREP)+A-1) * RCW + 1
                    IF (CARITH) THEN
                       S1(T1AI) = S1(T1AI)
     &            + DDOT(NINT,BUF1(CDIK),RCW,BUF3(AK),RCW*NVO(IRP))
     &            - DDOT(NINT,BUF1(CDIK+1),RCW,BUF3(AK+1),RCW*NVO(IRP))
                       S1(T1AI+1) = S1(T1AI+1)
     &            + DDOT(NINT,BUF1(CDIK+1),RCW,BUF3(AK),RCW*NVO(IRP))
     &            + DDOT(NINT,BUF1(CDIK),RCW,BUF3(AK+1),RCW*NVO(IRP))
                    ELSE
                       S1(T1AI) = S1(T1AI)
     &            + DDOT(NINT,BUF1(CDIK),1,BUF3(AK),NVO(IRP))
                    ENDIF
                 ENDDO
                 DO A = 1, NV(KREP)
                    AI = (IIVO(KREP,IREP)+(I-1)*NV(KREP)+A-1) * RCW + 1
                    T1AK = (IVO(KREP)+(K-1)*NV(KREP)+A-1) * RCW + 1
                    IF (CARITH) THEN
                       S1(T1AK) = S1(T1AK)
     &            - DDOT(NINT,BUF1(CDIK),RCW,BUF3(AI),RCW*NVO(IRP))
     &            + DDOT(NINT,BUF1(CDIK+1),RCW,BUF3(AI+1),RCW*NVO(IRP))
                       S1(T1AK+1) = S1(T1AK+1)
     &            - DDOT(NINT,BUF1(CDIK+1),RCW,BUF3(AI),RCW*NVO(IRP))
     &            - DDOT(NINT,BUF1(CDIK),RCW,BUF3(AI+1),RCW*NVO(IRP))
                    ELSE
                       S1(T1AK) = S1(T1AK)
     &            - DDOT(NINT,BUF1(CDIK),1,BUF3(AI),NVO(IRP))
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
 20      CONTINUE
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
 30   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BFLOPS(CPU4V,CPUH1,CPUH2,NIT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculate performance of largest contraction
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 CPU4V,CPUH1,CPUH2
      INTEGER NIT
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 MFLOPS,DRCW2,DRCW,DNVVT,DNOOT,DNVO,DNVOOT,DNV,DNO
      REAL*8 MFLOPB,MFLOPH,DMVO
      CHARACTER*30 CONTR
C
C---------------Executable code--------------------------------------
C
C     Prevent division by zero in case somebody could not figure out
C     how the cpu timer works on the system.
C
      IF (CPU4V.EQ.0.0.OR.CPUH1.EQ.0.0.OR.CPUH2.EQ.0.0) RETURN
      MFLOPB = 0.0
      MFLOPH = 0.0
C
      DRCW2 = DBLE(2*RCW**2)
      DRCW  = DBLE(RCW)
      DO IRP = 1, NREP

         DNVVT = DBLE(NVVT(IRP))
         DNOOT = DBLE(NOOT(IRP))
         DNVO  = DBLE(NVO(IRP))
         DNVOOT= DBLE(NVOOT(IRP))
         DNV   = DBLE(NV(IRP))
         DNO   = DBLE(NO(IRP))
         DMVO  = DBLE(MVO(IRP))

C        Compute number of operations for BINTM
         MFLOPB = MFLOPB + DNVVT  * DNVVT * DNOOT * DRCW2
         MFLOPB = MFLOPB + DNVVT  * DNVO  * DNOOT * DRCW2
         MFLOPB = MFLOPB + DNVOOT * DNV   * DNO   * DRCW2
         MFLOPB = MFLOPB + DNVVT  * DNOOT * DRCW

C        Compute number of operations for HINTM
         MFLOPH = MFLOPH + DMVO  * DMVO * DMVO * DRCW2

      ENDDO 
C
C     Convert to MFlops
C
      MFLOPB = MFLOPB * 0.000001
      MFLOPH = MFLOPH * 0.000001
C
      WRITE (IW,1000)
C
      DO I = 1, 3
         IF (I.EQ.1) THEN
            MFLOPS = MFLOPB * DBLE(NIT) / CPU4V
            CONTR = 'VVVV+VOVV (in B: includes I/O)'
         ELSEIF (I.EQ.2) THEN
            MFLOPS = MFLOPH * DBLE(NIT) / CPUH1
            CONTR = 'VOVO (in H: only XGEMM)       '
         ELSE
            MFLOPS = MFLOPH * DBLE(NIT) / CPUH2
            CONTR = 'VOVO (in T2EQN: includes sort)'
         ENDIF
C
         IF (MFLOPS.GE.1.0E6) THEN
            WRITE (IW,1001) CONTR,MFLOPS/1.0E6,'T'
         ELSEIF (MFLOPS.GE.1.0E5) THEN
            WRITE (IW,1003) CONTR,MFLOPS/1.0E3,'G'
         ELSEIF (MFLOPS.GE.1.0E4) THEN
            WRITE (IW,1002) CONTR,MFLOPS/1.0E3,'G'
         ELSEIF (MFLOPS.GE.1.0E3) THEN
            WRITE (IW,1001) CONTR,MFLOPS/1.0E3,'G'
         ELSEIF (MFLOPS.GE.1.0E2) THEN
            WRITE (IW,1003) CONTR,MFLOPS,'M'
         ELSEIF (MFLOPS.GE.1.0E1) THEN
            WRITE (IW,1002) CONTR,MFLOPS,'M'
         ELSE
            WRITE (IW,1001) CONTR,MFLOPS,'M'
         ENDIF
C
      ENDDO
C
 1000 FORMAT (//' Performance of BLAS GEMM in the largest contractions',
     &//' Contraction type',T45,'Performance')
 1001 FORMAT (1X,A,T40,F10.3,3X,A1,'flop/s') 
 1002 FORMAT (1X,A,T40,F10.2,3X,A1,'flop/s') 
 1003 FORMAT (1X,A,T40,F10.1,3X,A1,'flop/s') 
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DELFCK (CLASS,DODEL,FXX)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Delete integrals that contribute to unphysical excitations
C     in the Fock space calculations.
C
C     The spinors space is arranged as
C
C
C      #       O         V
C     -----------------------
C      1      IO        AV
C      2      AO        IV
C      3      AV        AO
C  
C      This routine deletes the AV integrals in the occupied space if
C      the integrals is to be used in a contraction over truly O-type
C      spinor, while in the virtual space it will delete the AO 
C      class, if desired. The parameter CLASS specifies the type
C      of the incoming integrals, the parameter DODEL specifies the
C      indices that are to be cleaned.
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FXX(*)
      CHARACTER*2 CLASS,DODEL
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
      DIMENSION NINDX(0:2,MXREP,2)
C
C---------------Executable code--------------------------------------
C
C     Return when we do not do Fock space calculations
C
      IF (EQNS.NE.'FOCKSP') RETURN
C
C     NINDX(1,IREP,INDEX) : Spinors to be kept
C     NINDX(2,IREP,INDEX) : Additional spinors to be zeroed.
C
      DO INDEX = 1, 2
         IF (CLASS(INDEX:INDEX).EQ.'O') THEN
            DO IREP = 1, NREP
               IF (DODEL(INDEX:INDEX).EQ.'D') THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ELSE
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ENDIF
            ENDDO
         ELSEIF (CLASS(INDEX:INDEX).EQ.'V') THEN
            DO IREP = 1, NREP
               IF (DODEL(INDEX:INDEX).EQ.'D') THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ELSE
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ENDIF
            ENDDO
         ELSE
            CALL QUIT ('DELFCK integral class ')
         ENDIF
      ENDDO
C
      IJ = 1
C
C     Loop over J : irrep, spinor class, spinors
C
      DO 10 JREP = 1, NREP
         IREP = JREP
         DO JC =  1, 2
         DO J =  NINDX(JC-1,JREP,2)+1, NINDX(JC,JREP,2)
C
C    Loop over I : (irrep fixed via JREP), spinor class
C    The loop over spinors is done within XCOPY
C
            DO IC = 1, 2
               IMIN = NINDX(IC-1,IREP,1)+1
               LNGTH=NINDX(IC,IREP,1)-IMIN+1
               IDEL = IC + JC
C
C              Only when IDEL = 2 all indices are allowed
C
               IF (LNGTH.GT.0.AND.IDEL.GT.2) 
     &            CALL XCOPY (LNGTH,A0,0,FXX(IJ),1)
               IF (LNGTH.GT.0) IJ = IJ + RCW * LNGTH
            ENDDO    ! IC
C
         ENDDO    ! J
         ENDDO    ! JC
 10   CONTINUE ! JREP
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DELINT (CLASS,DODEL,VXXXX,IREPVV,ISTART,NINT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Delete integrals that contribute to unphysical excitations
C     in the Fock space calculations.
C
C     The spinors space is arranged as
C
C
C      #       O         V
C     -----------------------
C      1      IO        AV
C      2      AO        IV
C      3      AV        AO
C  
C      This routine deletes the AV integrals in the occupied space if
C      the integrals is to be used in a contraction over truly O-type
C      spinor, while in the virtual space it will delete the AO 
C      class, if desired. The parameter CLASS specifies the type
C      of the incoming integrals, the parameter DODEL specifies the
C      indices that are to be cleaned.
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VXXXX(*)
      CHARACTER*4 CLASS,DODEL
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
      DIMENSION NINDX(0:2,MXREP,4)
      LOGICAL TRIANG(2),ALLOWED,ALLOW1,ALLOW2
C
C---------------Executable code--------------------------------------
C
C     Return when we do not do Fock space calculations
C
      IF (EQNS.NE.'FOCKSP') RETURN
C
      TRIANG(1) = CLASS(1:1).EQ.CLASS(2:2)
      TRIANG(2) = CLASS(3:3).EQ.CLASS(4:4)
C
C     NINDX(1,IREP,INDEX) : Spinors to be kept
C     NINDX(2,IREP,INDEX) : Additional spinors to be zeroed.
C
      INDEXS = 1
      INDEXE = 2
      DO INDEXP = 1, 2
      DO INDEX = INDEXS, INDEXE
         IF (CLASS(INDEX:INDEX).EQ.'O') THEN
            DO IREP = 1, NREP
               IF (DODEL(INDEX:INDEX).EQ.'D') THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ELSE
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ENDIF
C
C              Triangular arrays need special treatment, we can have the
C              case that the to be deleted index is related by anti-
C              symmetry to the index to be kept. 
C              We treat this as the condition 'DD'to have the proper
C              loop structure but modify the delete condition somewhat
C    
               IF (TRIANG(INDEXP).AND.
     &             (DODEL(INDEXS:INDEXE).EQ.'DK'.OR.
     &              DODEL(INDEXS:INDEXE).EQ.'KD')) THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSO(IREP)
                  NINDX(2,IREP,INDEX) = NO(IREP)
               ENDIF
            ENDDO
         ELSEIF (CLASS(INDEX:INDEX).EQ.'V') THEN
            DO IREP = 1, NREP
               IF (DODEL(INDEX:INDEX).EQ.'D') THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ELSE
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ENDIF
C
C              The same treatment for the triangular condition as above
C    
               IF (TRIANG(INDEXP).AND.
     &             (DODEL(INDEXS:INDEXE).EQ.'DK'.OR.
     &              DODEL(INDEXS:INDEXE).EQ.'KD')) THEN
                  NINDX(0,IREP,INDEX) = 0
                  NINDX(1,IREP,INDEX) = NSV(IREP)
                  NINDX(2,IREP,INDEX) = NV(IREP)
               ENDIF
            ENDDO
         ELSE
            CALL QUIT ('DELINT integral class ')
         ENDIF
      ENDDO
      INDEXS = INDEXS + 2
      INDEXE = INDEXE + 2
      ENDDO
C
      IJKL = 1
C
C     Loop over compound irreps
C
      DO 40 KLREP = 1, NREP
         KL = 0
C
C        Loop over L : irrep, spinor classes, spinors
C
         DO 30 LREP = 1, NREP
         KREP = MULTB(LREP,KLREP+NREP,2)
         IF (TRIANG(2).AND.KREP.LT.LREP) GOTO 30
            
         DO LC = 1, 2
         DO L =  NINDX(LC-1,LREP,4)+1, NINDX(LC,LREP,4)
C
C           Loop over K: irrep(fixed by KLREP,JREP),spinor class,spinors
C
            DO KC =  1, 2
            KMIN = NINDX(KC-1,KREP,3)+1
            IF (KREP.EQ.LREP.AND.TRIANG(2)) KMIN = MAX0(KMIN,L+1)
C
C           Condition for keep/delete of this index pair
C
            IF (KC+LC.EQ.2) THEN
               ALLOW2 = .TRUE.
            ELSEIF (KC+LC.EQ.3.AND.TRIANG(2).AND.
     &             (DODEL(3:4).EQ.'DK'.OR.
     &              DODEL(3:4).EQ.'KD')) THEN
               ALLOW2 = .TRUE.
            ELSE
               ALLOW2 = .FALSE.
            ENDIF
C
            DO 20 K = KMIN, NINDX(KC,KREP,3)
C
C           For the VOVV and VVVV integrals we have only a batch in
C           memory, jump out if we are outside this area
C
            KL = KL + 1
            IF (CLASS.EQ.'VOVV'.OR.CLASS.EQ.'VVVV') THEN
               IF (KLREP.NE.IREPVV) GOTO 40
               IF (KL.LE.ISTART.OR.KL.GT.ISTART+NINT) GOTO 20
            ENDIF
C
C    Loop over J : irrep, spinor class, spinors
C
               DO 10 JREP = 1, NREP
               IREP = MULTB(JREP,KLREP+NREP,2)
               IF (IREP.LT.JREP.AND.TRIANG(1)) GOTO 10
C
               DO JC =  1, 2
               DO J =  NINDX(JC-1,JREP,2)+1, NINDX(JC,JREP,2)
C
C    Loop over I : (irrep fixed via KLREP and JREP), spinor class
C    The loop over spinors is done within XCOPY
C
                  DO IC = 1, 2
                     IMIN = NINDX(IC-1,IREP,1)+1
                     IF (IREP.EQ.JREP.AND.TRIANG(1)) IMIN=MAX0(IMIN,J+1)
                     LNGTH=NINDX(IC,IREP,1)-IMIN+1
C
C                    Condition for keep/delete of this index pair
C
                     IF (IC+JC.EQ.2) THEN
                        ALLOW1 = .TRUE.
                     ELSEIF (IC+JC.EQ.3.AND.TRIANG(1).AND.
     &                      (DODEL(1:2).EQ.'DK'.OR.
     &                       DODEL(1:2).EQ.'KD')) THEN
                        ALLOW1 = .TRUE.
                     ELSE
                        ALLOW1 = .FALSE.
                     ENDIF
C
                     ALLOWED = ALLOW1.AND.ALLOW2
                     IF (LNGTH.GT.0.AND.(.NOT.ALLOWED))
     &                  CALL XCOPY (LNGTH,A0,0,VXXXX(IJKL),1)
                     IF (LNGTH.GT.0) IJKL = IJKL + RCW * LNGTH
                  ENDDO    ! IC
C
               ENDDO    ! J
               ENDDO    ! JC
 10            CONTINUE ! JREP
C                    
 20         CONTINUE ! K
            ENDDO ! KC
C
         ENDDO    ! L
         ENDDO    ! LC
 30      CONTINUE ! LREP
C
 40   CONTINUE ! KLREP
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PAIREN (T1,T2,FVO,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Attribute correlation energy to electron pairs
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2(*)
      REAL*8 T1(*),T2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "param.inc"
#include "files.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 EI,EIJ
C
C---------------Executable code--------------------------------------
C
      WRITE (IW,1000)
C
      CALL XCOPY (NDIMT1,FVO,1,BUF1,1)
      CALL DELFCK ('VO','DD',BUF1)
C
      WRITE (IW,2000)
      AI = 1
      DO IRP = 1, NREP
      DO I = 1, NO(IRP)
         EI = AR0
         DO A = 1, NV(IRP)
             EI = EI + BUF1(AI) * T1(AI)
             IF (CARITH) THEN
                EI = EI + BUF1(AI+1) * T1(AI+1)
             ENDIF
             AI = AI + RCW
         ENDDO
         IF (ABS(EI).GT.1.D-6) WRITE (IW,2001) IRP,I,EI
      ENDDO
      ENDDO
C
      CALL GETVVOO (BUF1)
      CALL DELINT ('VVOO','DDDD',BUF1,0,0,0)
      CALL GETTAU  (T1,T2,BUF2)
C
      WRITE (IW,3000)
      ABIJ = 1
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            EIJ = AR0
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            DO B = 1, NV(BRP)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  EIJ = EIJ + BUF1(ABIJ) * BUF2(ABIJ)
                  IF (CARITH) THEN
                     EIJ = EIJ + BUF1(ABIJ+1) * BUF2(ABIJ+1)
                  ENDIF
                  ABIJ = ABIJ + RCW
               ENDDO
            ENDDO
 20         CONTINUE
            IF (ABS(EIJ).GT.1.D-6) WRITE (IW,3001) IRP,I,JRP,J,EIJ
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      WRITE (IW,4000)
C
 1000 FORMAT (//' Individual terms in the correlation energy')
 2000 FORMAT (/'   Irrep Spinor',5X,'Relaxation energy (T1-amplitudes)')
 2001 FORMAT (I7,I8,F30.15)
 3000 FORMAT (/'   Irrep1 Spinor1 Irrep2 Spinor2',5X,
     &'Correlation energy (T2-amplitudes)')
 3001 FORMAT (2(I7,I8),F30.15)
 4000 FORMAT (//)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE S1S2_ADDFO (FVO,BUF1,BUF2,S1,S2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Add the first order wavefunction contribution to S1 and S2.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 BUF1(*),BUF2,FVO(*)
      REAL*8 S1(*),S2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C---------------------------------------
C Only the master gets this contribution
C---------------------------------------
      IF (MYPROC.NE.MASTER) RETURN
C-------------------------------
C S(A,I) = S(A,I) + F(A,I)
C-------------------------------
      IF (EQNS.EQ.'FOCKSP') THEN
C       Delete first order Hamiltonian parts
C       first for Fockspace
        CALL XCOPY (NDIMT1,FVO,1,BUF1,1)
        CALL DENOMF (BUF1,BUF1,BUF2,BUF1,BUF2,2)
        CALL XAXPY (NDIMT1,A1,BUF1,1,S1,1)
      ELSE
        CALL XAXPY (NDIMT1,A1,FVO,1,S1,1)
      ENDIF
C-------------------------------
C S(AB,IJ) = S(AB,IJ) + W(AB,IJ)
C-------------------------------
      CALL GETVVOO (BUF2)
      CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)
C
      RETURN
      END
