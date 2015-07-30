!dirac_copyright_start
!      Copyright (c) 2014 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FOLDED (ICURFSS,IHEQ,T1,T2,BUF1,BUF2,BUF3,
     &                   DS1,DS2,S1,S2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate T*H  terms ("folded diagrams") for FS CC
C     (0,1) and (1,0) complete, others to be done..
C
C     Notation for the folded diagrams :
C     - All Occupied      : I,J,K,L
C     - All Virtual       : A,B,C,D
C     - Active Occupied   : M,N,O,P
C     - Active Virtual    : E,F,G,H
C     - Inactive Occupied : Q,R,S,T
C     - Inactive Virtual  : V,W,X,Y
C
C---------------Routines called----------------------------------------
C
C     BLAS routines :                 XGEMM, XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Ephraim Eliav
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*),S1(*),S2(*),DS1(*),DS2(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "ihm.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL TEQ
C
C---------------Executable code--------------------------------------
C
C  All terms can be derived the regular CCSD contraction by making appropriate substitutions I->E, A->M, etc.
C
C  The Hamiltonian and wave functions elements are all contained in the S1, S2 and T1, T2 arrays, respectively.
C  For each term we also indicate which part of the array needs to be extracted (e.g. a=2, for the
C  second class of virtual indices which are the inactive virtuals).
C
C--------------------------------------------------------------------
      CALL XCOPY (NDIMT1,A0,0,DS1,1)
      CALL XCOPY (NDIMT2,A0,0,DS2,1)
C--------------------------------------------------------------------
C
C  1. Contractions of the 1-electron Hamiltonian
C
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C     Term 1 : S1(V,E) <---  - T1(V,F) * H1(F,E) 
C              contributes in sector (0,1) only
C
C     S1=S1(V,E) (a=2; i=3)
C     H1=S1(E,F) (a=1; i=3)
C
C     Extract the Hamiltonian
      CALL EXTRT1 (S1,BUF1,1,1,3,3,LFA)
C     Extract the wave function
      CALL EXTRT1 (T1,BUF2,2,2,3,3,LFA)
C     Calculate the diagram
      CALL CNTRCT ('N','N',NIV,NAV,NAV,-A1,BUF2,BUF1,A0,BUF3,NREP)
C     Insert the updated wave function
      CALL EXTRT1 (DS1,BUF3,2,2,3,3,LTR)
C--------------------------------------------------------------------
C     Term 2 : S2(AB,IE) <---  - P(34) T2(AB,IF) * H1(F,E) 
C              contributes in sectors (0,1), (1,1) and (0,2)
C
C     S2=T2(AB,IE) (a=1,3; b=1,3; i=1,2,3; j=3)
C     H1=S1(E,F)   (a=1; i=3)
C
C     Extract the wave function, the Hamiltonian is still in BUF1
      CALL EXTRT2R(T2,BUF3,1,3,1,3,3,3,LFA)
C     Sort and contract, after contraction put updated wf back.
      CALL SRT1S3 (NREP,MULTB,LFA,NVVT,NO,NAV,NVVOT,
     &             KVVOAVT,KKVVOT,BUF3,BUF2)
      CALL CNTRCT('N','N',NVVOT,NAV,NAV,-A1,BUF2,BUF1,A0,BUF3,NREP)
      CALL SRT1S3(NREP,MULTB,LTR,NVVT,NO,NAV,NVVOT,
     &            KVVOAVT,KKVVOT,BUF3,BUF2)
      CALL EXTRT2R(DS2,BUF2,1,3,1,3,3,3,LTR)
C--------------------------------------------------------------------
C     Term 3 : S1(M,Q) <---  H1(M,N) * T1(N,Q) 
C              contributes in sector (1,0) only
C
C     S1=S1(M,Q) (a=3; i=1)
C     H1=S1(M,N) (a=3; i=2)
C
C     Extract the Hamiltonian
      CALL EXTRT1 (S1,BUF1,3,3,2,2,LFA)
C     Extract the wave function
      CALL EXTRT1 (T1,BUF2,3,3,1,1,LFA)
      CALL CNTRCT ('N','N',NAO,NIO,NAO,A1,BUF1,BUF2,A0,BUF3,NREP)
C     Insert the updated wave function
      CALL EXTRT1 (DS1,BUF3,3,3,1,1,LTR)
C--------------------------------------------------------------------
C     Term 4 : S2(MB,IJ) <---  P(12) H1(M,N) * T2(NB,IJ)
C              contributes in sectors (1,0), (1,1) and (2,0)
C
C     S2=T2(MB,IJ) (a=3; b=1,2,3; i=1,3; j=1,3)
C     H1=S1(M,N)   (a=3; i=2)
C
C     Extract the wave function, the Hamiltonian is still in BUF1
      CALL EXTRT2L(T2,BUF3,3,3,1,3,1,3,LFA)
C     Sort and contract, after contraction put updated wf back.
      CALL SRT1S2 (NREP,MULTB,LFA,NAOV,NAO,NV,NOOT,LAOVOOT,
     &             LLVOOT,BUF3,BUF2)
      CALL CNTRCT ('N','N',NAO,NVOOT,NAO,A1,BUF1,BUF2,A0,BUF3,NREP)
      CALL SRT1S2 (NREP,MULTB,LTR,NAOV,NAO,NV,NOOT,LAOVOOT,
     &             LLVOOT,BUF3,BUF2)
      CALL EXTRT2L(DS2,BUF2,3,3,1,3,1,3,LTR)
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C
C  2. Contractions of the 2-electron Hamiltonian
C
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C     Term 1 : S2(AB,EF) <---  - T2TAU(AB,GH) * H2(GH,EF) 
C              contributes in sector (0,2)
C
C     S2=S2(AB,EF) (a=1,2; b=1,2 {-11}; i=3; j=3)
C     H2=S2(GH,EF) (a=1, b=1, i=3, j=3)
C     FSTAU(AB,GH)=T2(AB,GH)+T1(A,G)*T1(B,H)-T1(A,H)*T1(B,G)
C                 +DELTA(A,G)*T1(B,H)-DELTA(A,H)*T1(B,G)
C                 +T1(A,G)*DELTA(B,H)-T1(A,H)*DELTA(B,G)
C
C     Make FSTAU
      CALL GET_FSTAU (T1,T2,BUF1,BUF2,BUF3)
C     Extract the wave function

      CALL EXTRT2T(BUF2,BUF3,1,2,3,3,LFA)
C     Extract the 2-particle Hamiltonian
      CALL EXTRT2T(S2,BUF1,1,1,3,3,LFA)
      CALL CNTRCT ('N','N',NSVSVT,NAVAVT,NAVAVT,-A1,BUF3,BUF1,A0,
     &             BUF2,NREP)
C     Insert the updated wave function
      CALL EXTRT2T(DS2,BUF2,1,2,3,3,LTR)
C---------------------------------------------
C     Term 2 : S2(MN,IJ) <---  - H2(MN,OP) * T2TAU(OP,IJ) 
C              contributes in sector (2,0)
C
C     S2=S2(MN,IJ) (a=3; b=3; i=1,2; j=1,2 {-22})
C     H2=T2(MN,OP) (a=3; b=3 ; i=2; j=2)
C
C     FSTAU(OP,IJ)=T2(OP,IJ)+T1(O,I)*T1(P,J)-T1(O,J)*T1(P,I)
C                 +DELTA(O,I)*T1(P,J)-DELTA(O,J)*T1(P,I)
C                 +T1(O,I)*DELTA(P,J)-T1(O,J)*DELTA(P,I)
C
C     Make FSTAU
      CALL GET_FSTAU (T1,T2,BUF1,BUF2,BUF3)
C     Extract the wave function
      CALL EXTRT2T(BUF2,BUF3,3,3,1,2,LFA)
C     Extract the 2-particle Hamiltonian
      CALL EXTRT2T(S2,BUF1,3,3,2,2,LFA)
      CALL CNTRCT ('N','N',NAOAOT,NSOSOT,NAOAOT,-A1,BUF1,BUF3,A0,
     &             BUF2,NREP)
C     Insert the updated wave function
      CALL EXTRT2T(DS2,BUF2,3,3,1,2,LTR)
C-------------------------------------------
C     Term 3 : S2(MB,IE) <---  H2(MH,OE) * T2TAU(OB,IH) 
C              contributes in sector (1,1)
C
C     This term is evaluated after a sort of the second and fourth indices :
C              S2'(ME,IB) <--- H2'(ME,OH) * T2TAU'(OH,IB)
C
C     S2=S2(MB,IE) (a=3; b=1,2; i=1,2; j=3 {exclude 3123, achieved in denomf})
C     H2=T2(MH,OE) (a=3; b=1 ; i=2; j=3)
C
C     FSTAU(OB,IH)=T2(OB,IH)+T1(O,I)*T1(B,H)-T1(O,H)*T1(B,I)
C                 +DELTA(O,I)*T1(B,H)
C                 +T1(O,I)*DELTA(B,H)
C
C     Note that T2TAU has fewer non-zero contributions than in terms 1 and 2 (due to the delta-funs)
C     The term -T1(O,H)*T1(B,I) is deleted because it will not contribute
C
C For the record : Ephraim's original formula's
C (1,1) Sector (the new definition)
c
c  2-electron diagram
c
C  S2=T2(AF,IN) (a=1,2; b=1 ; i=1,2; j=2  {exclude-1122})
C  H2=T2(EF,MN) (a=1; b=1; i=2; j=2)
C
C-------------------------------------------
C (1) S2(AF,IN)=S2(AF,IN)+S2TAU(AE,IM)*H2(EF,MN)  ! this one is implemented below
C (2) S2(AF,IN)=S2(AF,IN)+P(34).S2(AF,IM)*H1(M,N) ! this diagram is contained in term 4 of the 1-e eff. Hamiltonian
C (3) S2(AF,IN)=S2(AF,IN)-P(12).S2(AE,IN)*H1(E,F) ! this diagram is contained in term 2 of the 1-e eff. Hamiltonian
C (4) S2(AF,MN)=S2(AF,MN)-P(12).S1(A,E)*H2(EF,MN) ! these terms are contained in
C (5) S2(EF,IN)=S2(EF,IN)+P(34).S1(I,M)*H2(EF,MN) ! the t2tau wf by using gettau2
C
C S2TAU(AF,IM)=S2(AF,IM)+S1(A,F)*S1(I,M)
C-------------------------------------------
C     Make FSTAU
      CALL GET_FSTAU (T1,T2,BUF1,BUF2,BUF3)
C     Extract the wave function
      CALL EXTRT2S(BUF2,BUF1,3,3,1,2,1,2,3,3,LFA)
C     Sort the wave function to ordering 4.
      CALL SRT16 (NREP,MULTB,LFA,LFA,NAO,NSV,NSO,NAV,MAOAV,
     &            JAOAVOV,JJAOAV,JJOV,BUF1,BUF2)
C     Extract the Hamiltonian
      CALL EXTRT2S(S2,BUF3,3,3,1,1,2,2,3,3,LFA)
C     Sort the Hamiltonian to ordering 4.
      CALL SRT16 (NREP,MULTB,LFA,LFA,NAO,NAV,NAO,NAV,MAOAV,
     &            JAOAVAOAV,JJAOAV,JJAOAV,BUF3,BUF1)
C     Do the contraction.
      CALL CNTRCT ('N','N',MAOAV,MOV,MAOAV,A1,BUF1,BUF2,A0,
     &              BUF3,NREP)
C     Sort the update to S2 to the standard order
      CALL SRT16 (NREP,MULTB,LFA,LTR,NAO,NSV,NSO,NAV,MAOAV,
     &            JAOAVOV,JJAOAV,JJOV,BUF1,BUF3)
C     Insert the updated wave function
      CALL EXTRT2S(DS2,BUF1,3,3,1,2,1,2,3,3,LTR)
C
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C
C  3. Deleting contributions to the effective Hamiltonian that were 
C     computed but should not be included. After this we can add the
C     contribution from the folded diagrams to the total
C
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C
      CALL DENOMF (BUF1,DS1,DS2,DS1,DS2,2)
      CALL XAXPY (NDIMT1,A1,DS1,1,S1,1)
      CALL XAXPY (NDIMT2,A1,DS2,1,S2,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HEFF1(ICURFSS,IPR,S1,S2,FVO,EPS,BUF1,BUF2,BUF3,
     &                 HEFF,EIGENVALUES,EIGENVECTORS)

C
      implicit real*8 (a-h,o-z)
C---------------Description--------------------------------------------
C
C     Extract the Hamiltonian dressing terms from the S1 and S2 arrays
C     and add the zeroth order (CI matrix) terms. Then find the eigen
C     values and eigenvectors by solving a generalized eigenvalue 
C     problem. Finish by analyzing the eigenvectors and computing the
C     transition moments.
C---------------Calling variables--------------------------------------
C
      REAL*8 HEFF(*),EIGENVALUES(0:*),EIGENVECTORS(*)
      REAL*8 S1(*),FVO(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
      CHARACTER*40 TIT
      INTEGER ICURFSS
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "results.inc"
#include "ccpar.inc"
#include "files.inc"
C
C---------------Local variables -----------------------------------
      PARAMETER (NDIMH=600,NDIMH2=NDIMH*NDIMH)
      DIMENSION INDX_ALL(NDIMH*MXREP),IRPS_ALL(NDIMH*MXREP)
      DIMENSION IND(MXREP),NDT_ALL(MXREP)
      REAL*8 E_OFF
C---------------Executable code--------------------------------------
      IF (ICURFSS.EQ.1) THEN
C-------------------------------------------
C (0,0)  Sector
C
C  Nothing to be done (only one eigenvalue)
C
         RETURN
C-------------------------------------------
C (0,1)  Sector
C
C  1-electron diagram
C
C  H01 = FVO(E,F) + S1(E,F) (a=1; i=3)
C-------------------------------------------
      ELSEIF (ICURFSS.EQ.2) THEN
        TIT = '(0h,1p) space, electron affinities    '
        CALL XCOPY (IAVAV(NREP+1),A0,0,HEFF,1)
C       Extract the Hamiltonian from S1
        CALL EXTRT1 (S1,HEFF,1,1,3,3,LFA)
C       Extract and add the Fock matrix elements
        CALL EXTRT1 (FVO,BUF1,1,1,3,3,LFA)
        CALL XAXPY (IAVAV(NREP+1),A1,BUF1,1,HEFF,1)
C-------------------------------------------
C (1,0) Sector
C
C  1-electron diagram
C
C  H10 = - FVO(M,N) - S1(M,N) (a=3; i=2)
C-------------------------------------------
      ELSEIF(ICURFSS.EQ.3) then
        TIT = '(1h,0p) space, ionization energies    '
        CALL XCOPY (IAOAO(NREP+1),A0,0,HEFF,1)
C       Extract the Hamiltonian
        CALL EXTRT1 (S1,HEFF,3,3,2,2,LFA)
C       Extract and add the Fock matrix elements
        CALL EXTRT1 (FVO,BUF1,3,3,2,2,LFA)
        CALL XAXPY (IAOAO(NREP+1),A1,BUF1,1,HEFF,1)
        CALL XSCAL (IAOAO(NREP+1),-A1,HEFF,1)
C-------------------------------------------
C (1,1) Sector
C
C  Slater rules for 2-e determinants !
C
C  General formula
C  H11(MF,OH) = H01(F,H)*DELTA(M,O) - H10(M,O)*DELTA(F,H)
C             - S2(MF,OH)
C
C  We need to extract the (a=3, b=1; i=2, j=3) part
C
C-------------------------------------------
      ELSEIF(ICURFSS.EQ.4) THEN
        TIT = '(1h,1p) space, excitation energies    '
C       Make the one electron Hamiltonian parts
        CALL XCOPY (NDIMT1,S1,1,BUF1,1)
        CALL XAXPY (NDIMT1,A1,FVO,1,BUF1,1)
        CALL SCALE_S1 (BUF1)
C       Make the effective Hamiltonian in the full array.
        CALL XCOPY (NDIMT2,S2,1,BUF2,1)
        CALL XSCAL (NDIMT2,-A1,BUF2,1)
        CALL GETTAU2 (BUF1,BUF3,BUF2)
C       Extract the Hamiltonian
        CALL EXTRT2S(BUF2,BUF3,3,3,1,1,2,2,3,3,LFA)
C       Sort the matrix to ordering 4 to diagonalize it.
        CALL SRT16 (NREP,MULTB,LFA,LFA,NAO,NAV,NAO,NAV,MAOAV,
     &              JAOAVAOAV,JJAOAV,JJAOAV,BUF3,HEFF)
C-------------------------------------------
C (0,2) Sector
C
C  Slater rules for 2-e determinants !
C
C  General formula
C  H02(EF,GH) = S1(E,G)*DELTA(F,H) + S1(F,H)*DELTA(E,G)
C             - S1(E,H)*DELTA(F,G) - S1(F,G)*DELTA(E,H)
C             + S2(EF,GH)
C
C  We need to extract the (a=1, b=1; i=3, j=3) part
C
C-------------------------------------------
      ELSEIF(ICURFSS.EQ.5) THEN
        TIT = '(0h,2p) space                         '
C       Add the one electron contributions
        CALL XCOPY (NDIMT1,S1,1,BUF1,1)
        CALL XAXPY (NDIMT1,A1,FVO,1,BUF1,1)
C       Make the effective Hamiltonian in the full array.
        CALL XCOPY (NDIMT2,S2,1,BUF2,1)
        CALL GETTAU2 (BUF1,BUF3,BUF2)
C       Extract the Hamiltonian
        CALL EXTRT2T(BUF2,HEFF,1,1,3,3,LFA)
C-------------------------------------------
C (2,0) Sector
C
C  Slater rules for 2-e determinants !
C
C  General formula
C  H20(MN,OP) = - S1(M,O)*DELTA(N,P) - S1(N,P)*DELTA(M,O)
C               + S1(M,P)*DELTA(N,O) + S1(N,O)*DELTA(M,P)
C               + S2(MN,OP)
C
C  We need to extract the (a=3, b=3; i=2, j=2) part
C
C-------------------------------------------
      ELSEIF(ICURFSS.EQ.6) THEN
        TIT = '(2h,0p) space'
C       Add the one electron contributions
        CALL XCOPY (NDIMT1,S1,1,BUF1,1)
        CALL XAXPY (NDIMT1,A1,FVO,1,BUF1,1)
        CALL XSCAL (NDIMT1,-A1,BUF1,1)
C       Make the effective Hamiltonian in the full array.
        CALL XCOPY (NDIMT2,S2,1,BUF2,1)
        CALL GETTAU2 (BUF1,BUF3,BUF2)
C       Extract the Hamiltonian
        CALL EXTRT2T(BUF2,HEFF,3,3,2,2,LFA)
      ENDIF
C
C.... Diagonalize the effective Hamiltonian (blocked by symmetry).
C
      DO 100 IRP=1,NREP
         IF (NFS_HE(IRP,ICURFSS).EQ.0) GOTO 100
         CALL SOLVE_CG (NFS_HE(IRP,ICURFSS),
     &                  HEFF(IFS_HE2(IRP,ICURFSS)*RCW+1),
     &                  EIGENVECTORS(IFS_HE2(IRP,ICURFSS)*RCW+1),
     &                  EIGENVALUES(IFS_HE(IRP,ICURFSS)*RCW+RCW))
100   CONTINUE

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SOLVE_CG (NDT,H,EVEC,EV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Find eigenvalues and eigenvectors of a complex or real general matrix
C     NDT : Dimension of the matrix
C     H   : Matrix (complex or real, square)
C     EV  : Eigenvalues (real part, if carith=T also imaginary)
C     EVEC: Eigenvectors (complex)
C
C---------------Routines called----------------------------------------
C
C      CG or RG (from EISPACK) if complex or real case, respectively
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C            Miro Ilias - complex eingenenergies carried in EV, more output
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "inpt.inc"
C
C---------------Calling variables--------------------------------------

      REAL*8 H(RCW,NDT,NDT),EVEC(RCW,NDT,NDT),EV(RCW,NDT)
C---------------Local variables -----------------------------------
      REAL*8, ALLOCATABLE :: EVR(:),EVI(:)
      REAL*8, ALLOCATABLE :: HR(:,:),HI(:,:)
      REAL*8, ALLOCATABLE :: VECR(:,:),VECI(:,:)
      REAL*8, ALLOCATABLE :: SCR(:,:)
      INTEGER , ALLOCATABLE :: INDX(:)
      INTEGER :: IERR
      REAL*8 FNOR(RCW),E_OFF
C---------------Executable code--------------------------------------
C
C....... Allocate memory for the diagonalization routine
C
      ALLOCATE (EVR(NDT),STAT=IERR)
      IF (IERR.NE.0) THEN
         CALL QUIT('Error in EVR allocation !')
      ENDIF
      ALLOCATE (EVI(NDT),STAT=IERR)
      IF (IERR.NE.0) THEN
        CALL QUIT('Error in EVI allocation !')
      ENDIF
      ALLOCATE (VECR(NDT,NDT),STAT=IERR)
      IF (IERR.NE.0) THEN
        CALL QUIT('Error in VECR allocation !')
      ENDIF
      ALLOCATE (SCR(NDT,3),STAT=IERR)
      IF (IERR.NE.0) THEN
        CALL QUIT('Error in SCR allocation !')
      ENDIF
      ALLOCATE (INDX(NDT),STAT=IERR)
      IF (IERR.NE.0) THEN
        CALL QUIT('Error in INDX allocation !')
      ENDIF

C
C....... Diagonalize effective Hamiltonian.
C
      IF  (CARITH) THEN
C        The diagonalization routine expects separate real and imaginary blocks,
C        make a local copy
         NDT2=NDT*NDT
         ALLOCATE (HR(NDT,NDT),STAT=IERR)
         IF (IERR.NE.0) THEN
           CALL QUIT('Error in HR allocation !')
         ENDIF
         ALLOCATE (HI(NDT,NDT),STAT=IERR)
         IF (IERR.NE.0) THEN
           CALL QUIT('Error in HI allocation !')
         ENDIF
         ALLOCATE (VECI(NDT,NDT),STAT=IERR)
         IF (IERR.NE.0) THEN
           CALL QUIT('Error in VECI allocation !')
         ENDIF
         CALL DCOPY (NDT2,H,2,HR,1)
         CALL DCOPY (NDT2,H(2,1,1),2,HI,1)

         IF (IPRNT.GE.10) THEN
CMI    ... control print out
           WRITE(IW,'(2X,A)')
     &'Output from SOLVE_CG - complex matrix (RCW=2) to be diagonalized'
           CALL OUTPUT(HR,1,NDT,1,NDT,NDT,NDT,1,IW)
           WRITE(IW,'(/2X,A)') 'Imaginary part'
           CALL OUTPUT(HI,1,NDT,1,NDT,NDT,NDT,1,IW)
         ENDIF

         CALL CG(NDT,NDT,HR,HI,EVR,EVI,1,VECR,VECI,
     &           SCR,SCR(1,2),SCR(1,3),IERR)
         DEALLOCATE (HR,STAT=IERR)
         IF (IERR.NE.0) THEN
            CALL QUIT('Error in HR deallocation !')
         ENDIF
         DEALLOCATE (HI,STAT=IERR)
         IF (IERR.NE.0) THEN
           CALL QUIT('Error in HI deallocation !')
         ENDIF
      ELSE
            IF (IPRNT.GE.10) THEN
CMI    ... control print out
           WRITE(IW,'(2X,A)')
     &'Output from SOLVE_CG - real matrix to be diagonalized'
             CALL OUTPUT(H,1,NDT,1,NDT,NDT,NDT,1,IW)
            ENDIF
            CALL RG(NDT,NDT,H,EVR,EVI,1,VECR,SCR,SCR(1,2),IERR)
      ENDIF
C
C.....Check for errors
C
      IF (IERR.NE.0) THEN
         WRITE (IW,'('' CG/RG IERR='',I4)') IERR
         CALL QUIT('Error in diagonalization of matrix')
      ENDIF
C
C.....Sort eigenvalues
C
      CALL INDEXX (NDT,EVR,INDX)
C
C.....Store the eigenvalues 
      DO I = 1, NDT
       EV(1,I) = EVR(INDX(I))
       IF (CARITH) EV(2,I) = EVI(INDX(I))
      ENDDO
C
C.....Store the eigenvectors in EVEC
C     
      DO I = 1, NDT
         DO J = 1, NDT
            EVEC(1,J,I) = VECR(J,INDX(I))
            IF (CARITH) EVEC(2,J,I) = VECI(J,INDX(I))
         ENDDO
      ENDDO
C
C.....Normalize the eigenvectors
C
      DO I = 1, NDT
         CALL XDOTC(FNOR,NDT,EVEC(1,1,I),1,EVEC(1,1,I),1)
         FNOR(1) = SQRT(1.D0 / FNOR(1))
         IF (CARITH) FNOR(2) = 0.D0
         CALL XSCAL(NDT,FNOR,EVEC(1,1,I),1)
      ENDDO
C
C.....Release the memory
C
      DEALLOCATE (EVR,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in EVR deallocation !')
      ENDIF
      DEALLOCATE (EVI,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in EVI deallocation !')
      ENDIF
      DEALLOCATE (VECR,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in VECR deallocation !')
      ENDIF
      IF (CARITH) THEN
        DEALLOCATE (VECI,STAT=IERR)
        IF (IERR.NE.0) THEN
          CALL QUIT('Error in VECI deallocation !')
        ENDIF
      ENDIF
      DEALLOCATE (SCR,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in SCR deallocation !')
      ENDIF
      DEALLOCATE (INDX,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in INDX deallocation !')
      ENDIF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FS_ANALYSIS1(ICURFSS,EPS,EIGENVALUES,EIGENVECTORS,
     &  NSP,IREPSPI)      
C
      implicit real*8 (a-h,o-z)
C---------------Description--------------------------------------------
C
C     Analyse the eigenvectors of the effective Hamiltonian and
C     compute oscilator strengths.
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "inpt.inc"
#include "complex.inc"
#include "symm.inc"
#include "results.inc"
#include "ccpar.inc"
#include "files.inc"
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*),EIGENVALUES(RCW,0:*),EIGENVECTORS(*) 
      INTEGER ICURFSS,NSP
C
C---------------Local variables -----------------------------------
      INTEGER , ALLOCATABLE :: INDX_ALL(:), IRPS_ALL(:)
      INTEGER :: IERR, ISIZE
      DIMENSION IND(MXREP)
      DIMENSION IREPSPI(NSP,MXREP,2)
      REAL*8 E_OFF
C---------------Executable code--------------------------------------

      ISIZE = IFS_HE(NREP+1,ICURFSS)+1
      ALLOCATE (INDX_ALL(ISIZE),STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in INDX_ALL allocation !')
      ENDIF
      ALLOCATE (IRPS_ALL(ISIZE),STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in IRPS_ALL allocation !')
      ENDIF

C     Initialize the reference energy

      E_OFF = ESCF + ECCSD
C
C     Print an ordered list of eigenvalues
C
      IF (ICURFSS.EQ.2.OR.ICURFSS.EQ.3) THEN
C        Odd electron case : fermion irreps
         CALL PRT_EV (IW,NREP,NFS_HE(1,ICURFSS),REPNA(1),
     &                E_OFF,ECCSDIM,RCW,
     &                EIGENVALUES(1:RCW,1:ISIZE),ISIZE,
     &                INDX_ALL,IRPS_ALL,IND,EFS,STATE_E_FSCC)
      ELSE IF (ICURFSS.EQ.4) THEN
C        Add energy of (0,0) sector to list of printed energies for the (1,1) sector
C        This makes the routine print excitation energies
         NFS_HE(1,ICURFSS) = NFS_HE(1,ICURFSS) + 1
         EIGENVALUES(1:RCW,0) = 0.D0
C        Even electron case : boson irreps
         CALL PRT_EV (IW,NREP,NFS_HE(1,ICURFSS),REPNA(NREP+1),
     &                E_OFF,ECCSDIM,RCW,
     &                EIGENVALUES(1:RCW,0:ISIZE),ISIZE,
     &                INDX_ALL,IRPS_ALL,IND,EFS,STATE_E_FSCC)
C        One should not forget this line : the real dimension of the first irrep should be retained !
         NFS_HE(1,ICURFSS) = NFS_HE(1,ICURFSS) - 1
      ELSE
C        Even electron case : boson irreps
         CALL PRT_EV (IW,NREP,NFS_HE(1,ICURFSS),REPNA(NREP+1),
     &                E_OFF,ECCSDIM,RCW,
     &                EIGENVALUES(1:RCW,1:ISIZE),ISIZE,
     &                INDX_ALL,IRPS_ALL,IND,EFS,STATE_E_FSCC)
      ENDIF
C     Detailed print-out including oscillator strengths
      CALL FS_ANALYSIS (ICURFSS,EPS,IPRNT,IFS_HE(NREP+1,ICURFSS),
     &                E_OFF,EIGENVALUES(1,1),EIGENVECTORS,NSP,IREPSPI)
      DEALLOCATE (INDX_ALL,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in INDX_ALL deallocation !')
      ENDIF
      DEALLOCATE (IRPS_ALL,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in IRPS_ALL deallocation !')
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMF (EPS,T1,T2,CT1,CT2,OPTION)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Option 1) :
C     Divide parts of T1 and T2 that correspond to WF by denominators
C     Parts that correspond to the effective Hamiltonian are untouched.
C     In the intermediate Hamiltonian IH-2 formalism we zero out Pi->Q
C     excitations (but keep Pm->Q) and corresponding Hamiltonian parts.
C
C     Option 2) :
C     Delete parts of T1 and T2 that correspond to the Hamiltonian.
C
C     Option 3) :
C     In the intermediate Hamiltonian (IH-2) formalism we zero out Pi->Q
C     excitations (but keep Pm->Q) and corresponding Hamiltonian parts.
C     NO division of the T1 and T2 by dets!
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher  & Ephraim Eliav
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
#include "files.inc"
#include "complex.inc"
#include "param.inc"
#include "ihm.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC,FAC1,FAC2,FAC3,VALUE,FACIH
      INTEGER IOI(0:3), IVA(0:3)
      INTEGER IOJ(0:3), IVB(0:3)
      CHARACTER*2 CLASS(3,2)
      LOGICAL DENO,ZERO,Z1,Z2,Z3,Z4,IHZERO
      LOGICAL DOIHDN
C
C---------------Executable code--------------------------------------
C
      CLASS(1,1) = 'Oi'
      CLASS(2,1) = 'Oa'
      CLASS(3,1) = 'Va'
      CLASS(1,2) = 'Va'
      CLASS(2,2) = 'Vi'
      CLASS(3,2) = 'Oa'
c
c
c The Fock Space blocks constraction:
c
c
c      #       O         V
c     -----------------------
c      1      IO        AV
c      2      AO        IV
c      3      AV        AO
C
C T1 part:
C
c
      II = 0
      AI = 0
C
      DO IRP = 1, NREP
c
      IOI(0)=0
      IOI(1)=NIO(IRP)+IOI(0)
      IOI(2)=NAO(IRP)+IOI(1)
      IOI(3)=NAV(IRP)+IOI(2)
c
      do ico=1,3
      DO I = IOI(ico-1)+1, IOI(ico)
         Z1 = .FALSE.
         II = II + 1
         FAC1 = EPS(II)
         ISECTP=0
         IF (DOIH.AND.ICO.EQ.3) THEN
             ISECTP=1
             IF (IPIORB(II).EQ.1) Z1=.TRUE.
         ENDIF
         AOFF = IO(NREP+1) + IV(IRP)
c
         IVA(0)=0
         IVA(1)=NAV(IRP)+IVA(0)
         IVA(2)=NIV(IRP)+IVA(1)
         IVA(3)=NAO(IRP)+IVA(2)
c
         do acv=1,3
C        Find out whether we need to divide by the denominators
C        leave Hamiltonian matrix elements untouched.
         IF (OPTION.EQ.1) THEN
            IF ((ICO.EQ.3.AND.ACV.EQ.3).OR.
     &          (ICO.EQ.3.AND.ACV.EQ.1).OR.
     &          (ICO.EQ.2.AND.ACV.EQ.3)) THEN
                DENO = .FALSE.
            ELSE
                DENO = .TRUE.
            ENDIF
         ELSE
            DENO = .FALSE. 
         ENDIF
C        Find out whether we need to zero non-wf parts
C        This will zero Hamiltonian elements
         IF (OPTION.EQ.2) THEN
            IF ((ICO.EQ.3.AND.ACV.EQ.3).OR.
     &          (ICO.EQ.3.AND.ACV.EQ.1).OR.
     &          (ICO.EQ.2.AND.ACV.EQ.3)) THEN
                ZERO = .TRUE.
            ELSE
                ZERO = .FALSE.
            ENDIF
         ELSE
            ZERO = .FALSE. 
         ENDIF
C
         DO A = IVA(acv-1)+1, IVA(acv)
             Z2 = .FALSE.
             AA = AOFF + A
             FAC = FAC1 - EPS(AA)
             FACIH = FAC
             ISECTH=0
             IF (DOIH.AND.ACV.EQ.3) THEN
                ISECTH=1
                IF (IPIORB(AA).EQ.1) Z2=.TRUE.
             ENDIF
             AI = AI + 1
             DOIHDN=.FALSE.
             IHZERO=.FALSE.
             DOIHDN = (Z1.OR.Z2).AND.IHSCHEME.EQ.1
             IHZERO = (Z1.OR.Z2).AND.IHSCHEME.EQ.2
             IF (DENO) THEN
               IF (DOIHDN) THEN
                IF (ISECTP.EQ.1.AND.ISECTH.EQ.0) FACIH=FAC+SHIFT_IH(1,1)
                IF (ISECTP.EQ.0.AND.ISECTH.EQ.1) FACIH=FAC+SHIFT_IH(1,2)
               ENDIF
                IF (CARITH) THEN
                    IF (IHZERO) CT1(AI) = A0
                    IF (DOIHDN.AND.DABS(FACIH-FAC).GT.1.D-7) THEN
                CT1(AI)=CT1(AI)/FACIH
                      IF (DABS(AIH).GT.1.D-9) THEN
                      IF (NIH.GT.0.AND.NIH.LT.100) THEN
                CT1(AI)=CT1(AI)*
     x                  (1.d0 - (AIH*(FACIH-FAC)/FACIH)**NIH)/
     x                  (1.d0 - AIH*(FACIH-FAC)/FACIH)
                       ELSEIF (NIH.GE.100) THEN
                CT1(AI)=CT1(AI)/
     x           (1.d0 - AIH*(FACIH-FAC)/FACIH)
                      ENDIF
                      ENDIF
                    ELSE
                    CT1(AI) = CT1(AI)/FAC
                    ENDIF
                 ELSE
                    IF (IHZERO) T1(AI) = A0
                    IF (DOIHDN.AND.DABS(FACIH-FAC).GT.1.D-7) THEN
                T1(AI)=T1(AI)/FACIH
                      IF (DABS(AIH).GT.1.D-9) THEN
                      IF (NIH.GT.0.AND.NIH.LT.100) THEN
                T1(AI)=T1(AI)*
     x                  (1.d0 - (AIH*(FACIH-FAC)/FACIH)**NIH)/
     x                  (1.d0 - AIH*(FACIH-FAC)/FACIH)
                       ELSEIF (NIH.GE.100) THEN
                T1(AI)=T1(AI)/
     x           (1.d0 - AIH*(FACIH-FAC)/FACIH)
                      ENDIF
                      ENDIF
                     ELSE
                    T1(AI) = T1(AI)/FAC
                    ENDIF
                ENDIF
             ENDIF
C
C
             IF (ZERO.OR.(OPTION.EQ.3.AND.IHZERO)) THEN
                IF (CARITH) THEN
                   CT1(AI) = A0
                ELSE
                   T1(AI) = A0
                ENDIF
             ENDIF
          ISECTP=0
          ISECTH=0
         ENDDO
         enddo
      ENDDO
      enddo
      ENDDO
C

c T2 part:
      ABIJ = 0
c
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
c
      IOJ(0)=0
      IOJ(1)=NIO(JRP)+IOJ(0)
      IOJ(2)=NAO(JRP)+IOJ(1)
      IOJ(3)=NAV(JRP)+IOJ(2)
c
      IOFF = IO(IRP) 
c
      IOI(0)=0
      IOI(1)=NIO(IRP)+IOI(0)
      IOI(2)=NAO(IRP)+IOI(1)
      IOI(3)=NAV(IRP)+IOI(2)
c

          do jco=1,3
          DO J = IOJ(jco-1)+1, IOJ(jco)
            Z1 = .FALSE.
            JJ = JJ + 1
            FAC1 = EPS(JJ)
            ISECTP1=0
            IF (DOIH.AND.JCO.EQ.3) THEN
                 ISECTP1=1
                 IF (IPIORB(JJ).EQ.1) Z1 = .TRUE.
            ENDIF
            icomin=1
            IF (IRP.EQ.JRP) icoMIN = jco
            do ico=icoMIN,3
            IMIN = IOI(ico-1)+1
            IF (IRP.EQ.JRP.and.ico.eq.jco) IMIN = J + 1
            DO I = IMIN, IOI(ico)
               Z2 = .FALSE.
               II = IOFF + I
               FAC2 = EPS(II) + FAC1
               ISECTP2=0
               IF (DOIH.AND.ICO.EQ.3) THEN
                  ISECTP2=1
                  IF (IPIORB(II).EQ.1) Z2 = .TRUE.
               ENDIF
c
             DO 20 BRP = 1, NREP
             BB = IV(BRP) + IO(NREP+1)
             ARP = MULTB(BRP,IJRP+NREP,2)
             IF (ARP.LT.BRP) GOTO 20
c
             IVB(0)=0
             IVB(1)=NAV(BRP)+IVB(0)
             IVB(2)=NIV(BRP)+IVB(1)
             IVB(3)=NAO(BRP)+IVB(2)
c
             AOFF = IV(ARP) + IO(NREP+1)
c
             IVA(0)=0
             IVA(1)=NAV(ARP)+IVA(0)
             IVA(2)=NIV(ARP)+IVA(1)
             IVA(3)=NAO(ARP)+IVA(2)
c
                  do bcv=1,3
                  DO B = IVB(bcv-1)+1, IVB(bcv)
                     Z3 = .FALSE.
                     BB = BB + 1
                     FAC3 = FAC2 - EPS(BB)
                     ISECTH1=0
                     IF (DOIH.AND.BCV.EQ.3) THEN
                       ISECTH1=1
                       IF (IPIORB(BB).EQ.1) Z3 = .TRUE.
                     ENDIF
                     acvMIN=1
                     IF (ARP.EQ.BRP) acvMIN = Bcv
                     do acv=acvMIN,3
C                    Find out whether we need to divide by the denoms
C                    leave Hamiltonian matrix elements untouched.
                     IF (OPTION.EQ.1) THEN
                 IF ( (ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.1.AND.BCV.EQ.1)
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.2.AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.2.AND.ACV.EQ.3.AND.BCV.EQ.1)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.2.AND.ACV.EQ.1.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.1)
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.3.AND.ACV.EQ.1.AND.BCV.EQ.3)
ceph &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.3))
     &            .OR.(             JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.             AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.             AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.3             ))
     &                  THEN
                            DENO = .FALSE.
                        ELSE
                            DENO = .TRUE.
                        ENDIF
                     ELSE
                        DENO = .FALSE. 
                     ENDIF
C                    Find out whether we need to zero non-wf parts
C                    This will zero out Hamiltonian elements
                     IF (OPTION.EQ.2) THEN
                 IF ( (ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.1.AND.BCV.EQ.1)      ! H2 (EF,GH) (0,2) sector
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.2.AND.ACV.EQ.3.AND.BCV.EQ.3)      ! H2 (MN,OP) (2,0) sector
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.2.AND.ACV.EQ.3.AND.BCV.EQ.1)      ! H2 (MF,GP) (1,1) sector
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.2.AND.ACV.EQ.1.AND.BCV.EQ.3)      ! H2 (EN,GP) (1,1) sector
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.1)      ! H2 (MF,OH) (1,1) sector
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.3.AND.ACV.EQ.1.AND.BCV.EQ.3)      ! H2 (EN,OH) (1,1) sector
ceph &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.3))     ! H2 (MN,GH)  ?
     &            .OR.(             JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.             AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.             AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.3             ))
     &                  THEN
                            ZERO = .TRUE.
                        ELSE
                            ZERO = .FALSE.
                        ENDIF
                     ELSE
                        ZERO = .FALSE. 
                     ENDIF
                     AMIN = IVA(acv-1)+1
                     IF (ARP.EQ.BRP.and.acv.eq.bcv) AMIN = B + 1
                        DO A = AMIN, IVA(acv)
                           Z4 = .FALSE.
                           AA = AOFF + A
                           FAC = FAC3 - EPS(AA)
                           FACIH = FAC
                           ISECTH2=0
                           IF (DOIH.AND.ACV.EQ.3) THEN
                            ISECTH2=1
                            IF (IPIORB(AA).EQ.1) Z4 = .TRUE.
                           ENDIF
                           ABIJ = ABIJ + 1
                           DOIHDN = .FALSE.
                           IHZERO = .FALSE.
                           DOIHDN = (Z1.OR.Z2.OR.Z3.OR.Z4)
     &                              .AND.IHSCHEME.EQ.1
                           IHZERO = (Z1.OR.Z2.OR.Z3.OR.Z4)
     &                              .AND.IHSCHEME.EQ.2
              ISECTP=ISECTP1+ISECTP2
              ISECTH=ISECTH1+ISECTH2
              IF (DENO) THEN
               IF (DOIHDN) THEN
                IF (ISECTP.EQ.1.AND.ISECTH.EQ.0) THEN
                    FACIH=FAC+SHIFT_IH(2,1)
                ELSEIF (ISECTP.EQ.0.AND.ISECTH.EQ.1) THEN
                    FACIH=FAC+SHIFT_IH(2,2)
                ELSEIF (ISECTP.EQ.2.AND.ISECTH.EQ.0) THEN
                   FACIH=FAC+SHIFT_IH(2,5)
                   IF (Z1.AND.Z2) FACIH=FACIH+SHIFT_IH(2,5)
                ELSEIF (ISECTP.EQ.0.AND.ISECTH.EQ.2) THEN
                   FACIH=FAC+SHIFT_IH(2,6)
                   IF (Z3.AND.Z4) FACIH=FACIH+SHIFT_IH(2,6)
                ELSEIF (ISECTP.EQ.1.AND.ISECTH.EQ.1) THEN
                   FACIH=FAC+SHIFT_IH(2,4)
                   IF ((Z1.AND.Z4).OR.(Z1.AND.Z3).OR.
     x                 (Z2.AND.Z4).OR.(Z2.AND.Z3)) 
     x             FACIH=FACIH+SHIFT_IH(2,4)
                ENDIF
               ENDIF
C                            Check for zero denominators
                             IF (DABS(FACIH).LT.1.D-9) THEN
                               WRITE (IW,1000) IRP,JRP,ARP,BRP,
     &                         I,J,A,B,
     &                         CLASS(ICO,1),CLASS(JCO,1),
     &                         CLASS(ACV,2),CLASS(BCV,2),
     &                         EPS(II),EPS(JJ),EPS(AA),EPS(BB),FAC,FACIH
                               CALL QUIT('Zero denominator in DENOMF')
                             ENDIF
               IF (CARITH) THEN
                IF (DOIHDN.AND.DABS(FACIH-FAC).GT.1.D-7) THEN
                CT2(ABIJ)=CT2(ABIJ)/FACIH
                      IF (DABS(AIH).GT.1.D-9) THEN
                      IF (NIH.GT.0.AND.NIH.LT.100) THEN
                CT2(ABIJ)=CT2(ABIJ)*
     x           (1.d0 - (AIH*(FACIH-FAC)/FACIH)**NIH)/
     x           (1.d0 - AIH*(FACIH-FAC)/FACIH)
                       ELSEIF (NIH.GE.100) THEN
                CT2(ABIJ)=CT2(ABIJ)/
     x           (1.d0 - AIH*(FACIH-FAC)/FACIH)
                      ENDIF
                      ENDIF
                 ELSE
                  CT2(ABIJ) = CT2(ABIJ)/FAC
                 ENDIF
                ELSE
                 IF (IHZERO) T2(ABIJ) = A0
                 IF (DOIHDN.AND.DABS(FACIH-FAC).GT.1.D-7) THEN
                 T2(ABIJ)=T2(ABIJ)/FACIH
                      IF (DABS(AIH).GT.1.D-9) THEN
                      IF (NIH.GT.0.AND.NIH.LT.100) THEN
                 T2(ABIJ)=T2(ABIJ)*
     x           (1.d0 - (AIH*(FACIH-FAC)/FACIH)**NIH)/
     x           (1.d0 - AIH*(FACIH-FAC)/FACIH)
                       ELSEIF (NIH.GE.100) THEN
                 T2(ABIJ)=T2(ABIJ)/
     x           (1.d0 - AIH*(FACIH-FAC)/FACIH)
                      ENDIF
                      ENDIF
                  ELSE
                   T2(ABIJ) = T2(ABIJ)/FAC
                 ENDIF
                ENDIF
C-----------------------------------------------------------------------
CThis may cause large output : activate when debugging..
C                            Check for large amplitudes
C                            IF (CARITH) THEN
C                               VALUE = CDABS(CT2(ABIJ))
C                            ELSE
C                               VALUE = ABS(T2(ABIJ))
C                            ENDIF
C                            IF (VALUE.GT.0.1D0) THEN
C                               WRITE (IW,1001) IRP,JRP,ARP,BRP,
C    &                          I,J,A,B,
C    &                          CLASS(ICO,1),CLASS(JCO,1),
C    &                          CLASS(ACV,2),CLASS(BCV,2),
C    &                          EPS(II),EPS(JJ),EPS(AA),EPS(BB),VALUE
C                            ENDIF
C-----------------------------------------------------------------------
              ENDIF
C
C
                           IF (ZERO.OR.(OPTION.EQ.3.AND.IHZERO)) THEN
                             IF (CARITH) THEN
                                CT2(ABIJ) = A0
                             ELSE
                                T2(ABIJ) = A0
                             ENDIF
                           ENDIF
                        ENDDO
                        enddo
                  ENDDO
                  enddo
 20          CONTINUE
               ENDDO
               enddo
          ENDDO
         enddo
 10   CONTINUE
      ENDDO
C
      IF (ABIJ.NE.NDIMT2) CALL QUIT ('ERROR in DENOMF')
C
      RETURN
 1000 FORMAT (//' Zero denominator found in DENOMF',
     &/T30,'I',14X,'J',14x,'A',14X,'B',
     &/' Irreps      : ',4I15,
     &/' Indices     : ',4I15,
     &/' Classes     : ',4(13X,A2),
     &/' Energies    : ',4F15.5,
     &/' Denominators: ',2E15.5,
     &//' Please choose a better model space ! We will stop now.')
 1001 FORMAT (//' Large amplitude in DENOMF',
     &/T30,'I',14X,'J',14x,'A',14X,'B',
     &/' Irreps      : ',4I15,
     &/' Indices     : ',4I15,
     &/' Classes     : ',4(13X,A2),
     &/' Energies    : ',4F15.5,
     &/' Absolute value ',E15.5)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXTRT1 (T1,FT1,ACVMIN,ACVMAX,ICOMIN,ICOMAX,DOINV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Extract the different subsets defined in Fock space CC from
C     the total T1 array.
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Ephraim Eliav
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),FT1(*)
      LOGICAL DOINV
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      dimension IOI(0:3), IVA(0:3)
C
C---------------Executable code--------------------------------------
c
c
c The Fock Space blocks construction:
c
c
c      #       O         V
c     -----------------------
c      1      IO        AV
c      2      AO        IV
c      3      AV        AO

C
      ioffs1=1
c
      DO IRP = 1, NREP
c
c
      IOI(0)=0
      IOI(1)=NIO(IRP)+IOI(0)
      IOI(2)=NAO(IRP)+IOI(1)
      IOI(3)=NAV(IRP)+IOI(2)
c
c
      IVA(0)=0
      IVA(1)=NAV(IRP)+IVA(0)
      IVA(2)=NIV(IRP)+IVA(1)
      IVA(3)=NAO(IRP)+IVA(2)
c
      AI= IVO(IRP) + IOI(icomin-1)*NV(IRP) +IVA(acvmin-1) + 1
      AI = (AI - 1) * RCW + 1
      LNGTH=IVA(ACVMAX)-IVA(ACVMIN-1)
c
      do ico=ICOMIN,ICOMAX
      DO I = IOI(ico-1)+1, IOI(ico)
         IF (DOINV) THEN
            CALL XAXPY (LNGTH,A1,FT1(ioffs1),1,T1(AI),1)
         ELSE  
            CALL XCOPY (LNGTH,T1(AI),1,FT1(ioffs1),1)
         ENDIF
         ioffs1=ioffs1 + RCW * LNGTH
         AI = AI+ RCW * NV(IRP)
      ENDDO
      enddo 
      ENDDO
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXTRT2R(T2,FT2,ABCVMIN,ABCVMAX,ICOMIN,ICOMAX,JCOMIN,
     &                   JCOMAX,DOINV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Extract or insert a subset defined in Fock space CC from (to)
C     the total T2 array.
C
C     This routine takes a subset of the indices 3 and 4 and 
C     puts it in a square array.
C
C     Coding of the different classes of spinors
C     This is used below when we indicate  ACVMIN, ACVMAX, etc.
c
c      #       O         V
c     -----------------------
c      1      IO        AV
c      2      AO        IV
c      3      AV        AO
C
C
C---------------Routines called----------------------------------------
C
C     XCOPY,XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T2(*),FT2(*)
      INTEGER ICOMIN,ICOMAX,JCOMIN,JCOMAX,ABCVMIN,ABCVMAX
      LOGICAL DOINV
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      dimension IOI(0:3), IVA(0:3)
      dimension IOJ(0:3), IVB(0:3)
      COMPLEX*16 FACTOR 
C
C---------------Executable code--------------------------------------
C
      IOFF2 = 1
C
C     Loop over compound irreps
C
      DO IJRP = 1, NREP
C
C    Loop over J : irrep, spinor space, spinors
C
         IJRPOFF = IVVOOTT(IJRP)
C
         DO 10 JRP = 1, NREP
         IRP = MULTB(JRP,IJRP+NREP,2)
C
         IOJ(0)=0
         IOJ(1)=NIO(JRP)+IOJ(0)
         IOJ(2)=NAO(JRP)+IOJ(1)
         IOJ(3)=NAV(JRP)+IOJ(2)
C
         IOI(0)=0
         IOI(1)=NIO(IRP)+IOI(0)
         IOI(2)=NAO(IRP)+IOI(1)
         IOI(3)=NAV(IRP)+IOI(2)
C
         DO JCO = JCOMIN, JCOMAX
         DO J = IOJ(JCO-1)+1, IOJ(JCO)
C
C    Loop over I : (irrep fixed via IJRP and JRP), spinor space, spinors
C
            DO ICO = ICOMIN, ICOMAX
            DO I = IOI(ICO-1)+1, IOI(ICO)
               IF (IRP.GT.JRP) THEN
C                 The amplitude ij is canonical
                  IJ = IIOOT(IRP,JRP) + (J-1)*NO(IRP) + I
                  FACTOR = A1
               ELSEIF (IRP.LT.JRP) THEN
C                 We have to take the amplitude ji
                  IJ = IIOOT(JRP,IRP) + (I-1)*NO(JRP) + J
                  FACTOR = - A1
               ELSEIF (I.GT.J) THEN
C                 The amplitude ij is canonical
                  IJ = IIOOT(IRP,JRP) + (J-1)*NO(IRP) + I - J*(J+1)/2
                  FACTOR = A1
               ELSEIF (I.LT.J) THEN
C                 We have to take the amplitude ji
                  IJ = IIOOT(JRP,IRP) + (I-1)*NO(JRP) + J - I*(I+1)/2
                  FACTOR = - A1
               ELSE
C                 The amplitude is zero, start at arbitrary position
                  IJ = 1
                  FACTOR = A0
               ENDIF
               IJOFF = IJRPOFF + (IJ-1) * NVVT(IJRP)
C
C    Loop over B : irrep, spinor space, spinors
C
               DO 20 BRP = 1, NREP
               ARP = MULTB(BRP,IJRP+NREP,2)
               IF (ARP.LT.BRP) GOTO 20
C
               IVB(0)=0
               IVB(1)=NAV(BRP)+IVB(0)
               IVB(2)=NIV(BRP)+IVB(1)
               IVB(3)=NAO(BRP)+IVB(2)
C
               IVA(0)=0
               IVA(1)=NAV(ARP)+IVA(0)
               IVA(2)=NIV(ARP)+IVA(1)
               IVA(3)=NAO(ARP)+IVA(2)
C
               DO BCV = ABCVMIN, ABCVMAX
               DO B = IVB(BCV-1)+1, IVB(BCV)
C
C    Loop over A : (irrep fixed via IJRP and BRP), spinor space
C    The loop over spinors is done within XCOPY/XAXPY
C
                  ACVSTRT = ABCVMIN
                  IF (ARP.EQ.BRP) ACVSTRT = BCV
                  DO ACV = ACVSTRT, ABCVMAX
                     AMIN = IVA(ACV-1)+1
                     IF (ARP.EQ.BRP.AND.ACV.EQ.BCV) AMIN = B + 1 
                     LNGTH=IVA(ACV)-AMIN+1
                     IF (LNGTH.LT.0) stop 'extrt2r'
                     AB = IIVVT(ARP,BRP) + (B-1)*NV(ARP) + AMIN
                     IF (ARP.EQ.BRP) AB = AB - B*(B+1)/2
                     ABIJ = (IJOFF + AB - 1) * RCW + 1
                     IF (DOINV) THEN
                        CALL XAXPY (LNGTH,FACTOR,FT2(IOFF2),1,
     &                                       T2(ABIJ),1)
                     ELSE
                        CALL XCOPY (LNGTH,T2(ABIJ),1,
     &                                    FT2(IOFF2),1)
                        CALL XSCAL (LNGTH,FACTOR,FT2(IOFF2),1)
                     ENDIF
                     IOFF2 = IOFF2 + RCW * LNGTH
                  ENDDO
C
               ENDDO
               ENDDO
 20            CONTINUE
C                    
            ENDDO
            ENDDO
C
         ENDDO
         ENDDO
 10      CONTINUE
C
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXTRT2L(T2,FT2,ACVMIN,ACVMAX,BCVMIN,BCVMAX,IJCOMIN,
     &                   IJCOMAX,DOINV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Extract or insert a subset defined in Fock space CC from (to)
C     the total T2 array.
C
C     This routine takes a subset of the indices 3 and 4 and 
C     puts it in a square array.
C
C     Coding of the different classes of spinors
C     This is used below when we indicate  ACVMIN, ACVMAX, etc.
c
c      #       O         V
c     -----------------------
c      1      IO        AV
c      2      AO        IV
c      3      AV        AO
C
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
      REAL*8 T2(*),FT2(*)
      INTEGER IJCOMIN,IJCOMAX,ACVMIN,ACVMAX,BCVMIN,BCVMAX
      LOGICAL DOINV
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      dimension IOI(0:3), IVA(0:3)
      dimension IOJ(0:3), IVB(0:3)
      COMPLEX*16 FACTOR
C
C---------------Executable code--------------------------------------
C
      IOFF2 = 1
C
C     Loop over compound irreps
C
      DO IJRP = 1, NREP
C
C    Loop over J : irrep, spinor space, spinors
C
         IJRPOFF = IVVOOTT(IJRP)
C
         DO 10 JRP = 1, NREP
         IRP = MULTB(JRP,IJRP+NREP,2)
         IF (IRP.LT.JRP) GOTO 10
C
         IOJ(0)=0
         IOJ(1)=NIO(JRP)+IOJ(0)
         IOJ(2)=NAO(JRP)+IOJ(1)
         IOJ(3)=NAV(JRP)+IOJ(2)
C
         IOFF = IO(IRP) 
C
         IOI(0)=0
         IOI(1)=NIO(IRP)+IOI(0)
         IOI(2)=NAO(IRP)+IOI(1)
         IOI(3)=NAV(IRP)+IOI(2)
C
         DO JCO = IJCOMIN, IJCOMAX
         DO J = IOJ(JCO-1)+1, IOJ(JCO)
C
C    Loop over I : (irrep fixed via IJRP and JRP), spinor space, spinors
C
            ICOSTRT = IJCOMIN
            IF (IRP.EQ.JRP) ICOSTRT = JCO
            DO ICO = ICOSTRT, IJCOMAX
            IMIN = IOI(ICO-1)+1
            IF (IRP.EQ.JRP.AND.ICO.EQ.JCO) IMIN = J + 1
            DO I = IMIN, IOI(ICO)
               IJ = IIOOT(IRP,JRP) + (J-1)*NO(IRP) + I
               IF (IRP.EQ.JRP) IJ = IJ - J*(J+1)/2
               IJOFF = IJRPOFF + (IJ-1) * NVVT(IJRP)
C
C    Loop over B : irrep, spinor space, spinors
C
               DO 20 BRP = 1, NREP
               ARP = MULTB(BRP,IJRP+NREP,2)
C
               IVB(0)=0
               IVB(1)=NAV(BRP)+IVB(0)
               IVB(2)=NIV(BRP)+IVB(1)
               IVB(3)=NAO(BRP)+IVB(2)
C
               IVA(0)=0
               IVA(1)=NAV(ARP)+IVA(0)
               IVA(2)=NIV(ARP)+IVA(1)
               IVA(3)=NAO(ARP)+IVA(2)
C
               DO BCV = BCVMIN, BCVMAX
               DO B = IVB(BCV-1)+1, IVB(BCV)
C
C    Loop over A : (irrep fixed via IJRP and BRP), spinor space
C    The loop over spinors is done within XCOPY/XAXPY
C
                  ACVSTRT = ACVMIN
                  DO ACV = ACVSTRT, ACVMAX
                  DO A = IVA(ACV-1)+1, IVA(ACV)
                     IF (ARP.GT.BRP) THEN
C                       The amplitude ab is canonical
                        AB = IIVVT(ARP,BRP) + (B-1)*NV(ARP) + A
                        FACTOR = A1
                     ELSEIF (ARP.LT.BRP) THEN
C                       We have to take the amplitude ba
                        AB = IIVVT(BRP,ARP) + (A-1)*NV(BRP) + B
                        FACTOR = - A1
                     ELSEIF (A.GT.B) THEN
C                       The amplitude ab is canonical
                        AB = IIVVT(ARP,BRP) + (B-1)*NV(ARP)+A-B*(B+1)/2
                        FACTOR = A1
                     ELSEIF (A.LT.B) THEN
C                       We have to take the amplitude ba
                        AB = IIVVT(BRP,ARP) + (A-1)*NV(BRP)+B-A*(A+1)/2
                        FACTOR = - A1
                     ELSE
C                       The amplitude is zero, start at arbitr. position
                        AB = 1
                        FACTOR = A0
                     ENDIF
                     ABIJ = (IJOFF + AB - 1) * RCW + 1
                     IF (DOINV) THEN
                        T2(ABIJ) = T2(ABIJ) + FACTOR * FT2(IOFF2)
                        IF (CARITH) 
     &                  T2(ABIJ+1) = T2(ABIJ+1) + FACTOR * FT2(IOFF2+1)
                     ELSE
                        FT2(IOFF2) = FACTOR * T2(ABIJ)
                        IF (CARITH) 
     &                  FT2(IOFF2+1) = FACTOR * T2(ABIJ+1)
                     ENDIF
                     IOFF2 = IOFF2 + RCW
                  ENDDO
                  ENDDO
C
               ENDDO
               ENDDO
 20            CONTINUE
C                    
            ENDDO
            ENDDO
C
         ENDDO
         ENDDO
 10      CONTINUE
C
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRAMPF (T1,T2,CT1,CT2,TRES)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Print out T1 and T2 amplitudes or hamiltonian matrix elements
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher  & Ephraim Eliav
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
      COMPLEX*16 CT1(*),CT2(*)
      REAL*8 TRES,ABS
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "files.inc"
#include "complex.inc"
#include "param.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER IOI(0:3), IVA(0:3)
      INTEGER IOJ(0:3), IVB(0:3)
      CHARACTER*2 CLASS(5)
C
C---------------Executable code--------------------------------------
C
      IF (EQNS.NE.'FOCKSP') RETURN
      CLASS(1) = 'Oi'
      CLASS(2) = 'Oa'
      CLASS(3) = 'Va'
      CLASS(4) = 'Vi'
      CLASS(5) = 'Oa'
c
c
c The Fock Space blocks constraction:
c
c
c      #       O         V
c     -----------------------
c      1      IO        AV
c      2      AO        IV
c      3      AV        AO
C
C
C T1 part:
C
c
      WRITE (IW,1000) 'T1'
      II = 0
      AI = 0
C
      DO IRP = 1, NREP
c
      IOI(0)=0
      IOI(1)=NIO(IRP)+IOI(0)
      IOI(2)=NAO(IRP)+IOI(1)
      IOI(3)=NAV(IRP)+IOI(2)
c
      do ico=1,3
      DO I = IOI(ico-1)+1, IOI(ico)
         II = II + 1
         AOFF = IO(NREP+1) + IV(IRP)
c
         IVA(0)=0
         IVA(1)=NAV(IRP)+IVA(0)
         IVA(2)=NIV(IRP)+IVA(1)
         IVA(3)=NAO(IRP)+IVA(2)
c
         do acv=1,3
         DO A = IVA(acv-1)+1, IVA(acv)
             AA = AOFF + A
             AI = AI + 1
             IF (CARITH.AND.ABS(CT1(AI)).GE.TRES) THEN
               WRITE (IW,1001) IRP,
     &                         A,I,
     &                         CLASS(ACV+2),
     &                         CLASS(ICO),CT1(AI)
             ELSEIF (.NOT.CARITH.AND.ABS(T1(AI)).GE.TRES) THEN
               WRITE (IW,1001) IRP,
     &                         A,I,
     &                         CLASS(ACV+2),
     &                         CLASS(ICO),T1(AI)
             ENDIF
         ENDDO
         enddo
      ENDDO
      enddo
      ENDDO
C
      WRITE (IW,1000) 'T2'
c T2 part:
      ABIJ = 0
c
      DO IJRP = 1, NREP
          DO 10 JRP = 1, NREP
          JJ = IO(JRP)
          IRP = MULTB(JRP,IJRP+NREP,2)
          IF (IRP.LT.JRP) GOTO 10
c
          IOJ(0)=0
          IOJ(1)=NIO(JRP)+IOJ(0)
          IOJ(2)=NAO(JRP)+IOJ(1)
          IOJ(3)=NAV(JRP)+IOJ(2)
c
          IOFF = IO(IRP) 
c
          IOI(0)=0
          IOI(1)=NIO(IRP)+IOI(0)
          IOI(2)=NAO(IRP)+IOI(1)
          IOI(3)=NAV(IRP)+IOI(2)
c

          DO JCO=1,3
          DO J = IOJ(JCO-1)+1, IOJ(JCO)
            JJ = JJ + 1
            ICOMIN=1
            IF (IRP.EQ.JRP) ICOMIN = JCO
            DO ICO=ICOMIN,3
            IMIN = IOI(ICO-1)+1
            IF (IRP.EQ.JRP.AND.ICO.EQ.JCO) IMIN = J + 1
               DO I = IMIN, IOI(ICO)
                  II = IOFF + I
c
                  DO 20 BRP = 1, NREP
                  B = IV(BRP) + IO(NREP+1)
                  ARP = MULTB(BRP,IJRP+NREP,2)
                  IF (ARP.LT.BRP) GOTO 20
c
                  IVB(0)=0
                  IVB(1)=NAV(BRP)+IVB(0)
                  IVB(2)=NIV(BRP)+IVB(1)
                  IVB(3)=NAO(BRP)+IVB(2)
c
                  AOFF = IV(ARP) + IO(NREP+1)
c
                  IVA(0)=0
                  IVA(1)=NAV(ARP)+IVA(0)
                  IVA(2)=NIV(ARP)+IVA(1)
                  IVA(3)=NAO(ARP)+IVA(2)
c
                  DO BCV=1,3
                  DO B = IVB(BCV-1)+1, IVB(BCV)
                     ACVMIN=1
                     IF (ARP.EQ.BRP) ACVMIN = BCV
                     DO ACV=ACVMIN,3
                     AMIN = IVA(ACV-1)+1
                     IF (ARP.EQ.BRP.AND.ACV.EQ.BCV) AMIN = B + 1
                        DO A = AMIN, IVA(ACV)
                           AA = AOFF + A
                           ABIJ = ABIJ + 1
                           IF (CARITH.AND.ABS(CT2(ABIJ)).GT.TRES) THEN
                             WRITE (IW,1002) ARP,BRP,IRP,JRP,
     &                          A,B,I,J,
     &                          CLASS(ACV+2),
     &                          CLASS(BCV+2),
     &                          CLASS(ICO),CLASS(JCO),
     &                          CT2(ABIJ)
                           ELSEIF (.NOT.CARITH.AND.
     &                             ABS(T2(ABIJ)).GT.TRES) THEN
                             WRITE (IW,1002) ARP,BRP,IRP,JRP,
     &                          A,B,I,J,
     &                          CLASS(ACV+2),
     &                          CLASS(BCV+2),
     &                          CLASS(ICO),CLASS(JCO),
     &                          T2(ABIJ) 
                           ENDIF
                        ENDDO
                        ENDDO
                  ENDDO
                  ENDDO
 20               CONTINUE
               ENDDO
               ENDDO
          ENDDO
          ENDDO
C
 10   CONTINUE
      ENDDO

C
      IF (ABIJ.NE.NDIMT2) CALL QUIT ('ERROR in PRAMPF')
C
      RETURN
 1000 FORMAT (//' Debug output : The amplitudes stored in the ',A2,
     &' array'/' Notation : irreps,indices,classes,value')
 1001 FORMAT (3I5,2(1X,A2),2G20.8)
 1002 FORMAT (8I5,4(1X,A2),2G20.8)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXTRT2T(T2,FT2,ABCVMIN,ABCVMAX,IJCOMIN,IJCOMAX,DOINV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Extract or insert a subset defined in Fock space CC from (to)
C     the total T2 array.
C
C     This routine takes a subset of the indexpairs 12 and 34 and
C     puts it in a triangular array.
C
C     Coding of the different classes of spinors
C     This is used below when we indicate  ACVMIN, ACVMAX, etc.
c
c      #       O         V
c     -----------------------
c      1      IO        AV
c      2      AO        IV
c      3      AV        AO
C
C
C---------------Routines called----------------------------------------
C
C     XCOPY,XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T2(*),FT2(*)
      INTEGER IJCOMIN,IJCOMAX,ABCVMIN,ABCVMAX
      LOGICAL DOINV
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      dimension IOI(0:3), IVA(0:3)
      dimension IOJ(0:3), IVB(0:3)
C
C---------------Executable code--------------------------------------
C
      IOFF2 = 1
C
C     Loop over compound irreps
C
      DO IJRP = 1, NREP
C
C    Loop over J : irrep, spinor space, spinors
C
         IJRPOFF = IVVOOTT(IJRP)
C
         DO 10 JRP = 1, NREP
         IRP = MULTB(JRP,IJRP+NREP,2)
         IF (IRP.LT.JRP) GOTO 10
C
         IOJ(0)=0
         IOJ(1)=NIO(JRP)+IOJ(0)
         IOJ(2)=NAO(JRP)+IOJ(1)
         IOJ(3)=NAV(JRP)+IOJ(2)
C
         IOI(0)=0
         IOI(1)=NIO(IRP)+IOI(0)
         IOI(2)=NAO(IRP)+IOI(1)
         IOI(3)=NAV(IRP)+IOI(2)
C
         DO JCO = IJCOMIN, IJCOMAX
         DO J = IOJ(JCO-1)+1, IOJ(JCO)
C
C    Loop over I : (irrep fixed via IJRP and JRP), spinor space, spinors
C
            ICOSTRT = IJCOMIN
            IF (IRP.EQ.JRP) ICOSTRT = JCO
            DO ICO = ICOSTRT, IJCOMAX
            IMIN = IOI(ICO-1)+1
            IF (IRP.EQ.JRP.AND.ICO.EQ.JCO) IMIN = J + 1
            DO I = IMIN, IOI(ICO)
               IJ = IIOOT(IRP,JRP) + (J-1)*NO(IRP) + I
               IF (IRP.EQ.JRP) IJ = IJ - J*(J+1)/2
               IJOFF = IJRPOFF + (IJ-1) * NVVT(IJRP)
C
C    Loop over B : irrep, spinor space, spinors
C
               DO 20 BRP = 1, NREP
               ARP = MULTB(BRP,IJRP+NREP,2)
               IF (ARP.LT.BRP) GOTO 20
C
               IVB(0)=0
               IVB(1)=NAV(BRP)+IVB(0)
               IVB(2)=NIV(BRP)+IVB(1)
               IVB(3)=NAO(BRP)+IVB(2)
C
               IVA(0)=0
               IVA(1)=NAV(ARP)+IVA(0)
               IVA(2)=NIV(ARP)+IVA(1)
               IVA(3)=NAO(ARP)+IVA(2)
C
               DO BCV = ABCVMIN, ABCVMAX
               DO B = IVB(BCV-1)+1, IVB(BCV)
C
C    Loop over A : (irrep fixed via IJRP and BRP), spinor space
C    The loop over spinors is done within XCOPY/XAXPY
C
                  ACVSTRT = ABCVMIN
                  IF (ARP.EQ.BRP) ACVSTRT = BCV
                  DO ACV = ACVSTRT, ABCVMAX
                     AMIN = IVA(ACV-1)+1
                     IF (ARP.EQ.BRP.AND.ACV.EQ.BCV) AMIN = B + 1 
                     LNGTH=IVA(ACV)-AMIN+1
                     IF (LNGTH.LT.0) stop 'extrt2t'
                     AB = IIVVT(ARP,BRP) + (B-1)*NV(ARP) + AMIN
                     IF (ARP.EQ.BRP) AB = AB - B*(B+1)/2
                     ABIJ = (IJOFF + AB - 1) * RCW + 1
                     IF (DOINV) THEN
                        CALL XAXPY (LNGTH,A1,FT2(IOFF2),1,T2(ABIJ),1)
                     ELSE
                        CALL XCOPY (LNGTH,T2(ABIJ),1,FT2(IOFF2),1)
                     ENDIF
                     IOFF2 = IOFF2 + RCW * LNGTH
                  ENDDO
C
               ENDDO
               ENDDO
 20            CONTINUE
C                    
            ENDDO
            ENDDO
C
         ENDDO
         ENDDO
 10      CONTINUE
C
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXTRT2S(T2,FT2,ACVMIN,ACVMAX,BCVMIN,BCVMAX,
     &                   ICOMIN,ICOMAX,JCOMIN,JCOMAX,DOINV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Extract or insert a subset defined in Fock space CC from (to)
C     the total T2 array.
C
C     This routine takes a subset of all indices and 
C     puts it in a square array.
C
C     Coding of the different classes of spinors
C     This is used below when we indicate  ACVMIN, ACVMAX, etc.
c
c      #       O         V
c     -----------------------
c      1      IO        AV
c      2      AO        IV
c      3      AV        AO
C
C
C---------------Routines called----------------------------------------
C
C     XCOPY,XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T2(*),FT2(*)
      INTEGER ICOMIN,ICOMAX,JCOMIN,JCOMAX,ACVMIN,ACVMAX,BCVMIN,BCVMAX
      LOGICAL DOINV
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
      dimension IOI(0:3), IVA(0:3)
      dimension IOJ(0:3), IVB(0:3)
      COMPLEX*16 FACTOR,FACTOR1
C
C---------------Executable code--------------------------------------
C
      IOFF2 = 1
C
C     Loop over compound irreps
C
      DO IJRP = 1, NREP
C
C    Loop over J : irrep, spinor space, spinors
C
         IJRPOFF = IVVOOTT(IJRP)
C
         DO 10 JRP = 1, NREP
         IRP = MULTB(JRP,IJRP+NREP,2)
C
         IOJ(0)=0
         IOJ(1)=NIO(JRP)+IOJ(0)
         IOJ(2)=NAO(JRP)+IOJ(1)
         IOJ(3)=NAV(JRP)+IOJ(2)
C
         IOI(0)=0
         IOI(1)=NIO(IRP)+IOI(0)
         IOI(2)=NAO(IRP)+IOI(1)
         IOI(3)=NAV(IRP)+IOI(2)
C
         DO JCO = JCOMIN, JCOMAX
         DO J = IOJ(JCO-1)+1, IOJ(JCO)
C
C    Loop over I : (irrep fixed via IJRP and JRP), spinor space, spinors
C
            DO ICO = ICOMIN, ICOMAX
            DO I = IOI(ICO-1)+1, IOI(ICO)
               IF (IRP.GT.JRP) THEN
C                 The amplitude ij is canonical
                  IJ = IIOOT(IRP,JRP) + (J-1)*NO(IRP) + I
                  FACTOR1 = A1
               ELSEIF (IRP.LT.JRP) THEN
C                 We have to take the amplitude ji
                  IJ = IIOOT(JRP,IRP) + (I-1)*NO(JRP) + J
                  FACTOR1 = - A1
               ELSEIF (I.GT.J) THEN
C                 The amplitude ij is canonical
                  IJ = IIOOT(IRP,JRP) + (J-1)*NO(IRP) + I - J*(J+1)/2
                  FACTOR1 = A1
               ELSEIF (I.LT.J) THEN
C                 We have to take the amplitude ji
                  IJ = IIOOT(JRP,IRP) + (I-1)*NO(JRP) + J - I*(I+1)/2
                  FACTOR1 = - A1
               ELSE
C                 The amplitude is zero, start at arbitrary position
                  IJ = 1
                  FACTOR1 = A0
               ENDIF
               IJOFF = IJRPOFF + (IJ-1) * NVVT(IJRP)
C
C    Loop over B : irrep, spinor space, spinors
C
               DO 20 BRP = 1, NREP
               ARP = MULTB(BRP,IJRP+NREP,2)
C
               IVB(0)=0
               IVB(1)=NAV(BRP)+IVB(0)
               IVB(2)=NIV(BRP)+IVB(1)
               IVB(3)=NAO(BRP)+IVB(2)
C
               IVA(0)=0
               IVA(1)=NAV(ARP)+IVA(0)
               IVA(2)=NIV(ARP)+IVA(1)
               IVA(3)=NAO(ARP)+IVA(2)
C
               DO BCV = BCVMIN, BCVMAX
               DO B = IVB(BCV-1)+1, IVB(BCV)
C
C    Loop over A : (irrep fixed via IJRP and BRP), spinor space
C    The loop over spinors is done within XCOPY/XAXPY
C
                  ACVSTRT = ACVMIN
                  DO ACV = ACVSTRT, ACVMAX
                  DO A = IVA(ACV-1)+1, IVA(ACV)
                     IF (ARP.GT.BRP) THEN
C                       The amplitude ab is canonical
                        AB = IIVVT(ARP,BRP) + (B-1)*NV(ARP) + A
                        FACTOR = A1 * FACTOR1
                     ELSEIF (ARP.LT.BRP) THEN
C                       We have to take the amplitude ba
                        AB = IIVVT(BRP,ARP) + (A-1)*NV(BRP) + B
                        FACTOR = - A1 * FACTOR1
                     ELSEIF (A.GT.B) THEN
C                       The amplitude ab is canonical
                        AB = IIVVT(ARP,BRP) + (B-1)*NV(ARP)+A-B*(B+1)/2
                        FACTOR = A1 * FACTOR1
                     ELSEIF (A.LT.B) THEN
C                       We have to take the amplitude ba
                        AB = IIVVT(BRP,ARP) + (A-1)*NV(BRP)+B-A*(A+1)/2
                        FACTOR = - A1 * FACTOR1
                     ELSE
C                       The amplitude is zero, start at arbitr. position
                        AB = 1
                        FACTOR = A0
                     ENDIF
                     ABIJ = (IJOFF + AB - 1) * RCW + 1
                     IF (DOINV) THEN
                        T2(ABIJ) = T2(ABIJ) + FACTOR * FT2(IOFF2)
                        IF (CARITH) 
     &                  T2(ABIJ+1) = T2(ABIJ+1) + FACTOR * FT2(IOFF2+1)
                     ELSE
                        FT2(IOFF2) = FACTOR * T2(ABIJ)
                        IF (CARITH) 
     &                  FT2(IOFF2+1) = FACTOR * T2(ABIJ+1)
                     ENDIF
                     IOFF2 = IOFF2 + RCW
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
 20            CONTINUE
C                    
            ENDDO
            ENDDO
C
         ENDDO
         ENDDO
 10      CONTINUE
C
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SCALE_S1 (S1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Scale the different subsets defined in Fock space CC with the 
C     appropriate sign : +1 for EA, -1 for IE.
C     Of course only the negative scaling is required.
C
C---------------Routines called----------------------------------------
C
C     XSCAL
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 S1(*)
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
c
c
c The Fock Space blocks construction:
c
c
c      #       O         V
c     -----------------------
c      1      IO        AV
c      2      AO        IV
c      3      AV        AO

C
      DO IRP = 1, NREP
C
C     We need the second class of occupied spinors : these are the active occupied spinors
C
      ISTART = NIO(IRP)
      IEND = ISTART + NAO(IRP)
C
C     Compute the address of the first amplitude in this irrep that needs scaling
C     The address consists of the pointer to the start of the irrep, the length of the blocks of
C     unscaled inactive occupieds and the skipped true virtuals (we need the third block of V)
C
      AI= IVO(IRP) + ISTART * NV(IRP) + NAV(IRP) + NIV(IRP) + 1
      AI = (AI - 1) * RCW + 1
C
      DO I = ISTART+1, IEND
         CALL XSCAL (NAO(IRP),-A1,S1(AI),1)
         AI = AI+ RCW * NV(IRP)
      ENDDO
C
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_FSTAU (T1,T2,BUF1,BUF2,BUF3)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Deletes Hamiltonian parts and forms T2TAU needed in folded diagrams
C     T2TAU will be placed in BUF2, BUF1 and BUF3 are scratch arrays
C
C---------------Routines called----------------------------------------
C
C      XCOPY, DEOMF, GETTAU and GETTAU2
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*), T2(*)
      REAL*8 BUF1(*), BUF2(*), BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
C     Copy the wave function plus Hamiltonian in T1 to BUF1
      CALL XCOPY (NDIMT1,T1,1,BUF1,1)
C     Store T2 first in BUF3, GETTAU will copy it to BUF2
      CALL XCOPY (NDIMT2,T2,1,BUF3,1)
C     Delete the Hamiltonian parts from these arrays
      CALL DENOMF (BUF2,BUF1,BUF3,BUF1,BUF3,2)
C     Make T2TAU
      CALL GETTAU (BUF1,BUF3,BUF2)
      CALL GETTAU2(BUF1,BUF3,BUF2)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FS_ANALYSIS (ICURFSS,EPS,IPR,NDT_TOT,E_OFF,
     &                        EV_ALL,H,NSP,IREPSPI)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Analyze the excited states and calculate transition moments
C
C     On input ICURFSS : currently active sector
C                  EPS : orbital energies
C                  IPR : print level
C              NDT_TOT : dimension of the CI matrix
C                E_OFF : CCSD energy to be added to eigenvalue
C               EV_ALL : eigenvalues for all irreps
C                    H : CI coefficients for all irreps (packed)
C   output   INTENSITY : intensity of the transition
C         TRANS_DIPOLE : transition dipoles (x,y,z)
C
C---------------Routines called----------------------------------------
C
C     XCOPY, FS_COMPUTE_INTENSITY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
CMI/July 2007 works only with real ECCSD
C
C---------------Calling variables--------------------------------------
C
#include "complex.inc"
#include "symm.inc"
      REAL*8 TRES,EPS(*),E_OFF
      REAL*8 EV_ALL(RCW,NDT_TOT),H(RCW,*) 
      REAL*8 INTENSITY(NDT_TOT),
     &       TRANS_DIPOLE(3,RCW,NDT_TOT)
      INTEGER ICURFSS,IPR,NSP
      DIMENSION IREPSPI(NSP,MXREP,2)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "param.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      PARAMETER (NDIMH=5000)
      INTEGER I_EIGENVALUE, J_EIGENVALUE, I_COEFF
      CHARACTER*50 LABEL(NDIMH)
      REAL*8 VALUE
C
C---------------Executable code--------------------------------------
 
#ifdef MOD_UNRELEASED
!     At the moment only sector 4 has intensities
      IF (ICURFSS.EQ.4) THEN
      CALL FS_COMPUTE_INTENSITY (ICURFSS,H,INTENSITY,TRANS_DIPOLE,
     &                           NDT_TOT,IREPSPI,NSP,E_OFF,EV_ALL)
      WRITE (IW,1100)
      IF (.NOT.CARITH) THEN
         WRITE (IW,1110)
      ELSE
         WRITE (IW,1111)
      ENDIF
      I_EIGENVALUE = 0
      DO IRP=1,NREP
!        Not all irreps may be spanned
         NDT = NFS_HE(IRP,ICURFSS)
         IF (NDT.EQ.0) CYCLE
         DO J_EIGENVALUE = 1, NDT
              I_EIGENVALUE = I_EIGENVALUE + 1
              IF (INTENSITY(I_EIGENVALUE).GT.1E-05) 
     &        WRITE (IW,1112) REPNA(NREP+IRP),J_EIGENVALUE,
     &            EV_ALL(1,I_EIGENVALUE),
     &            INTENSITY(I_EIGENVALUE),
     &            ((TRANS_DIPOLE(IX,IY,I_EIGENVALUE),IX=1,3),IY=1,RCW)
         END DO
      END DO
      END IF ! Intensities for excitation energies
#endif
C
C     Now print an analysis of the vectors
C     Initialize and get labels for the excited state coefficiens
      IF (IPR.LE.0) THEN
         TRES = 0.1
      ELSEIF (IPR.LE.2) THEN
         TRES = 0.001
      ELSE
         TRES = 0.0
      END IF

      I_EIGENVALUE = 0
      I_COEFF = 0
      IJ = 1
      CALL GET_FSLABEL (ICURFSS,EPS,LABEL,NDIMH)

      WRITE (IW,1000) TRES
      DO IRP=1,NREP
 
C     Get the dimension of the space
 
         NDT = NFS_HE(IRP,ICURFSS)
         IF (NDT.EQ.0) CYCLE

C     Add an offset for irrep names for bosonic irreps

         IF (ICURFSS.LE.3) THEN
          IRPOFF = 0
         ELSE
          IRPOFF = NREP
         ENDIF

C        Loop over the eigenvalues

         DO J_EIGENVALUE = 1, NDT
            I_EIGENVALUE = I_EIGENVALUE + 1
C           I have seen crazy ouputs with large Fock space calcs, put a limit
            IF (J_EIGENVALUE.LT.100) THEN
            WRITE (IW,1010) REPNA(IRPOFF+IRP),J_EIGENVALUE,
     &             E_OFF+EV_ALL(1,I_EIGENVALUE),EV_ALL(1,I_EIGENVALUE)
C           Loop over the coefficients and print the interesting ones
            DO J_COEFF = 1, NDT
               K_COEFF = I_COEFF + J_COEFF
               VALUE = ABS(H(1,IJ))
               IF (CARITH) VALUE = VALUE + ABS(H(2,IJ))
               IF (CARITH .AND. VALUE.GT.TRES) THEN
                   WRITE (IW,'(2F15.5,3X,A50)') 
     &             H(1,IJ),H(2,IJ),LABEL(K_COEFF)
               ELSE IF (VALUE.GT.TRES) THEN              
               WRITE (IW,'(F15.5,3X,A50)') H(1,IJ),LABEL(K_COEFF)
               END IF
               IJ = IJ + 1
            END DO
            ELSE
               IJ = IJ + NDT
            END IF
         END DO
         I_COEFF = I_COEFF + NDT
      END DO
 
 1000 FORMAT (//' Analysis of eigenvectors of the effective Hamiltonian'
     &/' First line  :  Energy, Eigenvalue'
     &/' other lines : Coefficients of contributing determinants',
     & ' (above a treshold of ',E6.1,')')
 1010 FORMAT (/' Irrep ',A4,' State ',I5,2F16.8)
#ifdef MOD_UNRELEASED
 1100 FORMAT (//' Excitation energies and intensities'/)
 1110 FORMAT (/' Irrep   State   Ex. Energy     f        ',
     & '  mu(x)     mu(y)     mu(z) '  )
 1111 FORMAT (/' Irrep   State   Ex. Energy   f     ',
     & 'mu_r(x)   mu_r(y)   mu_r(z)    mu_i(x)    mu_i(y)    mu_i(z)')
 1112 FORMAT (2X,A4,I8,2f12.6,6f10.5)
#endif
      END
#ifdef MOD_UNRELEASED
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FS_COMPUTE_INTENSITY (ICURFSS,EIGENVECTORS,INTENSITY,
     &                                 TRANS_DIPOLE,NDT_TOT,IREPSPI,NSP,
     &                                 E_OFF,EV_ALL) 
c
      IMPLICIT INTEGER (A-Z)   
C
C---------------Description--------------------------------------------
C
C     Read dipole moment in MO basis and contract with CI vectors
C
C     Author : Luuk Visscher
C
C---------------Routines called--------------------------------------
C---------------Last modified-----------------------------------------
C---------------Common Blocks------------------------------------
C
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "dcbxpr.h"
c
c---------------Calling Variables-----------------------------------
c 
      INTEGER ICURFSS,NSP,NDT_TOT,IREPSPI(NSP,MXREP,2)
      REAL*8 EIGENVECTORS(RCW,IFS_HE2(NREP+1,ICURFSS))
      REAL*8 E_OFF,EV_ALL(RCW,NDT_TOT)
!output
      REAL*8 INTENSITY(NDT_TOT),
     &       TRANS_DIPOLE(3,RCW,NDT_TOT)
c
c---------------Local Variables-----------------------------------------------
c
c    NB: note that this scheme is treating a H-eff as were it a Hermitian CI matrix, a probably better approach
c        is to use the left and righthand side eigenvectors, but this needs to be worked out more carefully.
c
      REAL*8, ALLOCATABLE:: PROP(:,:,:)
c
c---------------Executable code--------------------------------------------------
c

      ALLOCATE(PROP(RCW,NSP,NSP),STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in PROP allocation!')
      ENDIF

      CALL DZERO (INTENSITY,NDT_TOT)
      CALL DZERO (TRANS_DIPOLE,3*RCW*NDT_TOT)

       DO EOPER=1, 3 !loop over components of dipole vector
      CALL DZERO (PROP,NSP*NSP*RCW)
c
c Read in dipole matrices in MO-basis
c
      CALL READ_PROP(NAMEE(EOPER),NSP,PROP)
c
c  Transition dipole moments for (0,1) or (1,0) sector
c
      IF ((ICURFSS.EQ.2).OR.(ICURFSS.EQ.3)) THEN 
! 
!      CALL TRANSITION_ZEROONE_ONEZERO (ICURFSS,EOPER,PROP,
!     &               TRANS_DIPOLE,EIGENVECTORS,NDT_TOT,NSP)
c
c  Transition dipole moments for (1,1) sector
c 
      ELSEIF (ICURFSS.EQ.4) THEN         

      CALL TRANSITION_ONEONE (EOPER,PROP,TRANS_DIPOLE,
     &                        EIGENVECTORS,NDT_TOT,IREPSPI,NSP)
      ENDIF
c
c Calculate intensity, gauge dependence? I = 2/3 * |mu|^2 * excitation energy
c
       IF (CARITH) THEN
        DO I=1,NDT_TOT
         INTENSITY(I)= INTENSITY(I) + 2.0D0/3.0D0 *
     &        (TRANS_DIPOLE(EOPER,1,I)**2
     &        + TRANS_DIPOLE(EOPER,2,I)**2)
     &        * EV_ALL(1,I)
        ENDDO
       ELSE
        DO I=1,NDT_TOT
          INTENSITY(I) = INTENSITY(I) + 2.0D0/3.0D0 *
     &     (TRANS_DIPOLE(EOPER,1,I)**2)*EV_ALL(1,I)
       ENDDO
      ENDIF 

       ENDDO !loop over components of dipole vector

      DEALLOCATE(PROP,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in PROP deallocation!')
      ENDIF

!!    Write excitation energies and transition dipoles

!      OPEN(20,FILE='FSCC_Spectra.out',STATUS='NEW')
!      WRITE(20,*) EV_ALL,TRANS_DIPOLE
!      CLOSE(20)

!!    PRINT SPECTRA

!      DO I=1,NDT_TOT
!       IF (ABS(INTENSITY(EOPER,I)).GT.0.0D0) THEN
!         CALL OUTPUT_SPECTRA(ICURFSS,IPR,NDT_TOT,EOPER,
!     &                      NFS_HE(NREP+1,ICURFSS))
!       ELSE 
!          WRITE(*,*) 'No spectra output as intensities are zero'
!       ENDIF
!      ENDDO

      RETURN
      END 
!C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TRANSITION_ONEONE (EOPER,PROP,TRANS_DIPOLE,
     &                       EIGENVECTORS,NDT_TOT,IREPSPI,NSP)
      IMPLICIT INTEGER (A-Z)
c
c     Contraction of dipole integrals with coefficients for sector (1,1)
c
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "dcbxpr.h"
c
c-------------------Calling Variables-------------------------------------------------
c
      INTEGER NSP,NDT_TOT,EOPER,IREPSPI(NSP,MXREP,2)
      REAL*8 EIGENVECTORS(RCW,IFS_HE2(NREP+1,4)),PROP(RCW,NSP,NSP)
!output
      !REAL*8 TRANS_DIPOLE(NEOPER,RCW,NDT_TOT)
      REAL*8 TRANS_DIPOLE(3,RCW,NDT_TOT) ! miro: out-of-bounds fix
c
c-------------------Local Variables--------------------------------------------------
c
      INTEGER, ALLOCATABLE:: IND(:,:)
      REAL*8 DDOT
c
c------------------------------------------------------------------------------------
c
c
c Contraction
c
      K = 0
      JI = 0
      DO IJREP = 1, NREP
      NDT = NFS_HE(IJREP,4) 
      IF (NDT.EQ.0) CYCLE
      ALLOCATE (IND(2,NDT))
      CALL PROPINDEX(IND,IJREP,NDT,IREPSPI,NSP)

      DO I = 1,NDT  !loop over eigenvalues
         K = K + 1
         DO J = 1,NDT !loop over states in the model space
           JI = JI + 1
           IF (CARITH) THEN
       TRANS_DIPOLE(EOPER,1,K) = TRANS_DIPOLE(EOPER,1,K)
     &        + PROP(1,IND(1,J),IND(2,J))*EIGENVECTORS(1,JI)
     &        + PROP(2,IND(1,J),IND(2,J))*EIGENVECTORS(2,JI)

       TRANS_DIPOLE(EOPER,2,K) = TRANS_DIPOLE(EOPER,2,K)
     &        + PROP(1,IND(1,J),IND(2,J))*EIGENVECTORS(2,JI)
     &        - PROP(2,IND(1,J),IND(2,J))*EIGENVECTORS(1,JI)
           ELSE

       TRANS_DIPOLE(EOPER,1,K) = TRANS_DIPOLE(EOPER,1,K)
     &       + PROP(1,IND(1,J),IND(2,J))*EIGENVECTORS(1,JI)
           ENDIF
             
          ENDDO !model space 
        ENDDO !eigenvalues
        DEALLOCATE(IND)
       enddo !nrep

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
       SUBROUTINE READ_PROP(DIPLEN,NSP,PROP)
C
       IMPLICIT INTEGER (A-Z)
C
C This routine is analogous to RDPROP in relccsd/ccints.F (without the sorting).
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "inpt.inc"
C
C---------------Calling variables------------------------------------
C
       INTEGER NSP  
       REAL*8 PROP(RCW,NSP,NSP)
       CHARACTER*8 DIPLEN
C
C---------------Local variables--------------------------------------
C
      LOGICAL TOBE
      CHARACTER*32 ACHAR
      REAL*8 DUMPROP(2,NSP,NSP)
      REAL*8 DUMMY,AMAXR,AMAXI,ANORM,APHASE(2)
C
C---------------Executable code--------------------------------------
C
C Property integrals are read in in MOLFDIR (RELCCSD) format
C
      INQUIRE(FILE='MDPROP',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
!        We need a warning here once this gets further automatized and the file should be present
         PROP = 0.0
         RETURN
      END IF

      OPEN (MDPROP,FILE='MDPROP',FORM='UNFORMATTED')
        READ (MDPROP,IOSTAT=IERR) ACHAR
        IF (IERR.NE.0) THEN
          WRITE (6,*) ' Error reading property ',DIPLEN,'on file MDPROP'
          CALL QUIT(' Error reading property integrals')
        ENDIF
        DO
          IF (ACHAR(1:8).NE.'********'.OR.ACHAR(25:32).NE.DIPLEN) THEN
            READ (MDPROP,IOSTAT=IERR) ACHAR
          ELSE
            EXIT
          ENDIF
        ENDDO
        IF (IERR.NE.0) THEN
          WRITE (6,*) ' Property ',DIPLEN, ' not found on file MDPROP'
          CALL QUIT(' Property integrals missing')
        ENDIF
        IF (IPRNT.GE.1) THEN
          WRITE (IW,*) 'Read integral type ', DIPLEN,
     &    'created ', ACHAR(9:16), 'storage info : ', ACHAR(17:24)
        ENDIF
        READ (MDPROP) DUMPROP
      CLOSE (MDPROP,STATUS='KEEP')
C
C Multiplication with phase factors
C For real symmetries the property integrals have to be multiplied with a phase factor
C
      IF (CARITH) THEN
        PROP = DUMPROP
      ELSE
C  
      APHASE(1) = D1
      APHASE(2) = D0
      AMAXI = D0
      IMAX = 0
      JMAX = 0
      DO I = 1,NSP
       DO J = 1,NSP
        IF (ABS(DUMPROP(2,I,J)).GT.ABS(AMAXI)) THEN
         AMAXR = DUMPROP(1,I,J) 
         AMAXI = DUMPROP(2,I,J)
         IMAX = I
         JMAX = J
        ENDIF
       ENDDO
      ENDDO
       
      IF (ABS(AMAXI).GT.ACCUR) THEN
       APHASE(1) = AMAXR
       APHASE(2) = AMAXI
       ANORM = SQRT(AMAXR*AMAXR + AMAXI*AMAXI)
       APHASE(1) = APHASE(1) / ANORM
       APHASE(2) = APHASE(2) / ANORM     
      ENDIF

      DO I = 1,NSP
        DO J = 1,NSP
          PROP(1,I,J) = DUMPROP(1,I,J)*APHASE(1) +
     &                  DUMPROP(2,I,J)*APHASE(2)
          AMAXI = D0
          DUMMY = DUMPROP(1,I,J)*APHASE(2) +
     &            DUMPROP(2,I,J)*APHASE(1)
           IF (ABS(DUMMY).GT.ABS(AMAXI)) THEN
            AMAXI = DUMMY
C consistency check
           IF (ABS(AMAXI).GT.ACCUR) THEN
              WRITE(*,*) I,J,AMAXI
              CALL QUIT('USE COMPLEX ARITHMETICS')
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      ENDIF !CARITH

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PROPINDEX(IND,IJREP,NDT,IREPSPI,NSP) 

      IMPLICIT INTEGER (A-Z)
c
c define offsets of property matrix needed in subroutine Transition_OneOne
c                                                      
c----------------Common Blocks-----------------------------------------
c 
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "dcbxpr.h"
c
c----------------Calling variables---------------------------------------------
c
      INTEGER NDT, IJREP
      INTEGER IND(2,NDT)
      DIMENSION IREPSPI(NSP,MXREP,2)

c
c  Determine indices for property matrix (in MOLFDIR format!)
c 
      II = 0
      DO JREP = 1, NREP
         IREP = MULTB(JREP,IJREP+NREP,1)
          DO J = 1, NAV(JREP)
             DO I = 1, NAO(IREP)
                II = II + 1
                ind(1,ii) = IREPSPI(nio(irep)+i,irep,1) ! irepspi points to the orginal list of spinors
                ind(2,ii) = IREPSPI(j,jrep,2) ! we need no offset here because active virtual come first
             END DO
          END DO
      END DO

      RETURN
      END
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_FSLABEL (ICURFSS,EPS,LABEL,NDIMH)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Define the labels of the model space vector
C
C     Author : Luuk Visscher
C
C
C---------------Common Blocks--------------------------------------
C
#include "symm.inc"
#include "files.inc"
C
C---------------Variables-----------------------------------------------
C
      CHARACTER*50 LABEL(NDIMH),NOT_DEFINED
      INTEGER ICURFSS,I
      REAL*8 EPS(*)
      DIMENSION IAO(MXREP),IAV(MXREP)

C     Determine offsets for start of orbital energies

      IAO(1) = NIO(1)
      IAV(1) = NIO(1) + NAO(1)
      DO IREP = 1, NREP-1
         IAO(IREP+1) = IAO(IREP) + NAO(IREP) + NAV(IREP)   + NIO(IREP+1)
         IAV(IREP+1) = IAV(IREP) + NAV(IREP) + NIO(IREP+1) + NAO(IREP+1)
      END DO

      II = 0
C     II counts the active orbitals
      IF (ICURFSS .EQ. 2) THEN
         DO IREP = 1, NREP
            DO I = 1, NAV(IREP)
               II = II + 1
               WRITE(LABEL(II),1001) REPNA(IREP),I,EPS(IAV(IREP)+I)
            END DO
         END DO
      
      ELSEIF (ICURFSS .EQ. 3) THEN
         DO IREP = 1, NREP
            DO I = 1, NAO(IREP)
               II = II + 1
               WRITE(LABEL(II),1001) REPNA(IREP),I,EPS(IAO(IREP)+I)
            END DO
         END DO

C     II counts single excitations
      ELSEIF (ICURFSS .EQ. 4) THEN
         DO IJREP =1, NREP
            DO JREP = 1, NREP
               IREP = MULTB(JREP,IJREP+NREP,1)
               DO J = 1, NAV(JREP)
                  DO I = 1, NAO(IREP)
                     II = II + 1
                     WRITE(LABEL(II),1011) 
     &                  REPNA(IREP),I,EPS(IAO(IREP)+I),
     &                  REPNA(JREP),J,EPS(IAV(JREP)+J)
                  END DO
               END DO
            END DO
         END DO

C     II counts the double creations
      ELSEIF (ICURFSS .EQ. 5) THEN
         DO IJREP =1, NREP
            DO JREP = 1, NREP
               IREP = MULTB(JREP,IJREP+NREP,2)
               IF (IREP.EQ.JREP) THEN
                  DO J = 1, NAV(JREP)
                     DO I = J+1, NAV(IREP)
                       II = II + 1
                       WRITE(LABEL(II),1002) 
     &                       REPNA(IREP),I,EPS(IAV(IREP)+I),
     &                       REPNA(JREP),J,EPS(IAV(JREP)+J)
                     END DO
                  END DO

               ELSEIF (IREP.GT.JREP) THEN
                  DO J = 1, NAV(JREP)
                     DO I = 1, NAV(IREP)
                       II = II + 1
                       WRITE(LABEL(II),1002) 
     &                       REPNA(IREP),I,EPS(IAV(IREP)+I),
     &                       REPNA(JREP),J,EPS(IAV(JREP)+J)
                     END DO
                  END DO
               END IF
            END DO
         END DO
C     II counts the double annihilations
      ELSEIF (ICURFSS .EQ. 6) THEN
         DO IJREP =1, NREP
            DO JREP = 1, NREP
               IREP = MULTB(JREP,IJREP+NREP,2)
               IF (IREP.EQ.JREP) THEN
                  DO J = 1, NAO(JREP)
                     DO I = J+1, NAO(IREP)
                       II = II + 1
                       WRITE(LABEL(II),1002) 
     &                       REPNA(IREP),I,EPS(IAO(IREP)+I),
     &                       REPNA(JREP),J,EPS(IAO(JREP)+J)
                     END DO
                  END DO
               ELSEIF (IREP.GT.JREP) THEN
                  DO J = 1, NAO(JREP)
                     DO I = 1, NAO(IREP)
                       II = II + 1
                       WRITE(LABEL(II),1002) 
     &                       REPNA(IREP),I,EPS(IAO(IREP)+I),
     &                       REPNA(JREP),J,EPS(IAO(JREP)+J)
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF

 1001 FORMAT ('| ',A4,' #',I4,' (',F8.3,') |',25X)
 1002 FORMAT ('| ',A4,' #',I4,' (',F8.3,'), ',
     &             A4,' #',I4,' (',F8.3,') |')
 1011 FORMAT (1X,A4,' #',I4,' (',F8.3,') -> ',
     &             A4,' #',I4,' (',F8.3,')  ')
      END
