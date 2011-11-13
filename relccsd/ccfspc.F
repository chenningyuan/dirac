C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
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
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
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
      real*8 ddot
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
      SUBROUTINE HEFF1(ICURFSS,IPR,S1,S2,FVO,H,BUF1,BUF2,BUF3)
C
      implicit real*8 (a-h,o-z)
C---------------Calling variables--------------------------------------
C
      REAL*8 S1(*),FVO(*),H(*)
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
      PARAMETER (NDIMH=300,NDIMH2=NDIMH*NDIMH)
      REAL*8 EV(NDIMH),EVI(NDIMH),HR(NDIMH2),HI(NDIMH2)
      REAL*8 SCR(3*NDIMH),VEC(NDIMH2),VECI(NDIMH2)
      REAL*8 EV_ALL(NDIMH*MXREP)
      DIMENSION INDX_ALL(NDIMH*MXREP),IRPS_ALL(NDIMH*MXREP)
      DIMENSION NDT_ALL(MXREP),IND(MXREP)
      REAL*8 FNOR,E_OFF
      INTEGER INDX(NDIMH)
C---------------Executable code--------------------------------------
C-------------------------------------------
C (0,1)  Sector
C
C  1-electron diagram
C
C  H01 = FVO(E,F) + S1(E,F) (a=1; i=3)
C-------------------------------------------
      IF (ICURFSS.EQ.1) then
        TIT = '(0h,1p) space, electron affinities    '
        CALL XCOPY (IAVAV(NREP+1),A0,0,H,1)
C       Extract the Hamiltonian
        CALL EXTRT1 (S1,H,1,1,3,3,LFA)
C       Extract and add the Fock matrix elements
        CALL EXTRT1 (FVO,BUF1,1,1,3,3,LFA)
        CALL XAXPY (IAVAV(NREP+1),A1,BUF1,1,H,1)
C-------------------------------------------
C (1,0) Sector
C
C  1-electron diagram
C
C  H10 = - FVO(M,N) - S1(M,N) (a=3; i=2)
C-------------------------------------------
      ELSEIF(ICURFSS.EQ.10) then
        TIT = '(1h,0p) space, ionization energies    '
        CALL XCOPY (IAOAO(NREP+1),A0,0,H,1)
C       Extract the Hamiltonian
        CALL EXTRT1 (S1,H,3,3,2,2,LFA)
C       Extract and add the Fock matrix elements
        CALL EXTRT1 (FVO,BUF1,3,3,2,2,LFA)
        CALL XAXPY (IAOAO(NREP+1),A1,BUF1,1,H,1)
        CALL XSCAL (IAOAO(NREP+1),-A1,H,1)
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
      ELSEIF(ICURFSS.EQ.2) then
        TIT = '(0h,2p) space                         '
C       Add the one electron contributions
        CALL XCOPY (NDIMT1,S1,1,BUF1,1)
        CALL XAXPY (NDIMT1,A1,FVO,1,BUF1,1)
C       Make the effective Hamiltonian in the full array.
        CALL XCOPY (NDIMT2,S2,1,BUF2,1)
        CALL GETTAU2 (BUF1,BUF3,BUF2)
C       Extract the Hamiltonian
        CALL EXTRT2T(BUF2,H,1,1,3,3,LFA)
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
      ELSEIF(ICURFSS.EQ.20) THEN
        TIT = '(2h,0p) space'
C       Add the one electron contributions
        CALL XCOPY (NDIMT1,S1,1,BUF1,1)
        CALL XAXPY (NDIMT1,A1,FVO,1,BUF1,1)
        CALL XSCAL (NDIMT1,-A1,BUF1,1)
C       Make the effective Hamiltonian in the full array.
        CALL XCOPY (NDIMT2,S2,1,BUF2,1)
        CALL GETTAU2 (BUF1,BUF3,BUF2)
C       Extract the Hamiltonian
        CALL EXTRT2T(BUF2,H,3,3,2,2,LFA)
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
      ELSEIF(ICURFSS.EQ.11) THEN
        TIT = '(1h,1p) space WARNING: MAYBE INCORRECT'
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
     &              JAOAVAOAV,JJAOAV,JJAOAV,BUF3,H)
      ENDIF
C   
      NDT_TOT = 0
      IOFF = 1
      DO 100 IRP=1,NREP
         IF (ICURFSS.EQ.1) THEN
          NDT = NAV(IRP)
         ELSEIF (ICURFSS.EQ.10) THEN
          NDT = NAO(IRP)
         ELSEIF (ICURFSS.EQ.2) THEN
          NDT = NAVAVT(IRP)
         ELSEIF (ICURFSS.EQ.20) THEN
          NDT = NAOAOT(IRP)
         ELSEIF (ICURFSS.EQ.11) THEN
          NDT = MAOAV(IRP)
         ELSE
          NDT = 0
         ENDIF
         NDT_ALL(IRP) = NDT
         IF (NDT.EQ.0) GOTO 100
C
C....... Diagonalize effective Hamiltonian.
C
         IF (NDT.GT.NDIMH) CALL QUIT ('Increase NDIMH in HEFF1')
         NDT2=NDT*NDT
         IF  (CARITH) THEN
            CALL DCOPY (NDT2,H(IOFF),2,HR,1)
            CALL DCOPY (NDT2,H(IOFF+1),2,HI,1)
            CALL CG(NDT,NDT,HR,HI,EV,EVI,1,VEC,VECI,
     &               SCR,SCR(NDT+1),SCR(2*NDT+1),IERR)
Cluuk       uncomment when debugging, gives huge output
Cluuk       call output (hr,1,ndt,1,ndt,ndt,ndt,1,iw)
Cluuk       call output (hi,1,ndt,1,ndt,ndt,ndt,1,iw)
         ELSE
Cluuk       uncomment when debugging, gives huge output
Cluuk       call output (h(ioff),1,ndt,1,ndt,ndt,ndt,1,iw)
            CALL RG(NDT,NDT,H(IOFF),EV,EVI,1,VEC,SCR,SCR(NDT+1),IERR)
         ENDIF
C
C........Check for errors
C
         IF (IERR.NE.0) THEN                
            WRITE (IW,'('' RG IERR='',I4)') IERR
            CALL QUIT('Error in diagonalization of eff. Hamiltonian')
         ENDIF
C
C........Sort eigenvalues 
C
         CALL INDEXX (NDT,EV,INDX)
C
         IF (IPR.GT.0) WRITE (IW,8300) TIT,IRP,(EV(INDX(I)),I=1,NDT)
C
C........Write  eigenvalues & vectors if desired.
C
         DO I = 1, NDT
C          Calculate pointer to start of vector
           IVEC = (INDX(I)-1)*NDT + 1
C          Normalize
           FNOR = DDOT (NDT,VEC(IVEC),1,VEC(IVEC),1)
           IF (CARITH) FNOR = FNOR + DDOT(NDT,VECI(IVEC),1,VECI(IVEC),1)
           FNOR = A1 / DSQRT(FNOR)
           CALL DSCAL (NDT,FNOR,VEC(IVEC),1)
           IF (CARITH) CALL DSCAL (NDT,FNOR,VECI(IVEC),1)
           IF (IPR.GT.1) 
     &         WRITE (IW,8400) I,EV(INDX(I)),(VEC(IVEC+J-1),J=1,NDT)
           IF (CARITH.AND.IPR.GT.1) 
     &         WRITE (IW,8401) I,EVI(INDX(I)),(VECI(IVEC+J-1),J=1,NDT)
         ENDDO
         IOFF=IOFF+NDT2*RCW
C
C........Store the eigenvalues for later processing by PRT_EV
C
         CALL DCOPY (NDT,EV,1,EV_ALL(NDT_TOT+1),1)
         NDT_TOT = NDT_TOT + NDT
 100  CONTINUE
C
C     Print an ordered list of eigenvalues to which the ccsd energy is
C     added.
C
      E_OFF = ESCF + ECCSD
      e_off = 0.d0
      IF (ICURFSS.EQ.1.OR.ICURFSS.EQ.10) THEN
         CALL PRT_EV (IW,NREP,NDT_ALL,REPNA,E_OFF,EV_ALL,
     &                INDX_ALL,IRPS_ALL,IND)
      ELSE
         CALL PRT_EV (IW,NREP,NDT_ALL,REPNA(NREP+1),E_OFF,EV_ALL,
     &                INDX_ALL,IRPS_ALL,IND)
      ENDIF
C
8300  FORMAT (//1x,'Diagonalization of ',A40,' irrep =',I2/
     &' Energies in atomic units:'
     &/(5F16.12))
8400  FORMAT (/1x,' Energy',i3,1h:,F20.12,' Real Vector:',/(1x,6f12.7))
8401  FORMAT (' Imaginary energy',i3,1h:,F20.12,' Imag Vector:',
     & (1x,6F12.7))
C
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
C
C     Option 2) :
C     Delete parts of T1 and T2 that correspond to the Hamiltonian.
C
C     Option 3) :
C     In the intermediate Hamiltonian formalism we zero out Pi->Q
C     excitations (but keep Pm->Q) and corresponding Hamiltonian parts.
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
      REAL*8 FAC,FAC1,FAC2,FAC3,VALUE
      INTEGER IOI(0:3), IVA(0:3)
      INTEGER IOJ(0:3), IVB(0:3)
      CHARACTER*2 CLASS(3,2)
      LOGICAL DENO,ZERO,Z1,Z2,Z3,Z4,IHZERO
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
         IF (DOIH.AND.IPIORB(II).EQ.1) THEN
             IF (ICO.EQ.3) Z1=.TRUE.
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
             IF (DOIH.AND.IPIORB(AA).EQ.1) THEN
                IF (ACV.EQ.3) Z2=.TRUE.
             ENDIF
             AI = AI + 1
             IF (DENO) THEN
                IF (CARITH) THEN
                       CT1(AI) = CT1(AI)/FAC
                 ELSE
                       T1(AI) = T1(AI)/FAC
                ENDIF
             ENDIF
C
             IHZERO = (Z1.OR.Z2).AND.OPTION.EQ.3
             IF (DOIH.AND.IHZERO) THEN
                IF (CARITH) THEN
                       CT1(AI) = A0
                 ELSE
                       T1(AI) = A0
                ENDIF
             ENDIF
C
             IF (ZERO) THEN
                IF (CARITH) THEN
                   CT1(AI) = A0
                ELSE
                   T1(AI) = A0
                ENDIF
             ENDIF
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
            IF (DOIH.AND.IPIORB(JJ).EQ.1) THEN
               IF (JCO.EQ.3) Z1 = .TRUE.
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
               IF (DOIH.AND.IPIORB(II).EQ.1) THEN
                  IF (ICO.EQ.3) Z2 = .TRUE.
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
                     IF (DOIH.AND.IPIORB(BB).EQ.1) THEN
                       IF (BCV.EQ.3) Z3 = .TRUE.
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
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.3))
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
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.3))     ! H2 (MN,GH)  ?
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
                           IF (DOIH.AND.IPIORB(AA).EQ.1) THEN
                              IF (ACV.EQ.3) Z4 = .TRUE.
                           ENDIF
                           ABIJ = ABIJ + 1
                           IF (DENO) THEN
C                            Check for zero denominators
                             IF (ABS(FAC).LT.1.D-12) THEN
                                WRITE (IW,1000) IRP,JRP,ARP,BRP,
     &                          I,J,A,B,
     &                          CLASS(ICO,1),CLASS(JCO,1),
     &                          CLASS(ACV,2),CLASS(BCV,2),
     &                          EPS(II),EPS(JJ),EPS(AA),EPS(BB),FAC
                                CALL QUIT('Zero denominator in DENOMF')
                             ENDIF
                             IF (CARITH) THEN
                                   CT2(ABIJ) = CT2(ABIJ)/FAC
                              ELSE
                                   T2(ABIJ) = T2(ABIJ)/FAC
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
                           IHZERO = (Z1.OR.Z2.OR.Z3.OR.Z4)
     &                              .AND.OPTION.EQ.3
                           IF (DOIH.AND.IHZERO) THEN
                             IF (CARITH) THEN
                                CT2(ABIJ) = A0
                             ELSE
                                T2(ABIJ) = A0
                             ENDIF
                           ENDIF
C
                           IF (ZERO) THEN
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
     &/' Denominator : ',E15.5,
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
c The Fock Space blocks constraction:
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
      REAL*8 TRES
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
             IF (CARITH.AND.ABS(CT1(AI)).GT.TRES) THEN
               WRITE (IW,1001) IRP,
     &                         A,I,
     &                         CLASS(ACV+2),
     &                         CLASS(ICO),CT1(AI)
             ELSEIF (.NOT.CARITH.AND.ABS(T1(AI)).GT.TRES) THEN
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
 1001 FORMAT (3I5,2(1X,A2),2F20.16)
 1002 FORMAT (8I5,4(1X,A2),2F20.16)
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
