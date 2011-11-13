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
C
C
C VERSION : $Revision: 1.32.2.2 $
C DATE    : $Date: 2004/09/01 08:17:51 $
C FILE    : ccgrad.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCDENSC(IOPT,T1,S1,T2,S2,
     &                   BUF1,BUF2,BUF3,DOO,DVV,DVO,DOV,DMO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Menu driven routine. Depending on value of IOPT
C
C     Calculate unperturbed density matrix at the
C       - SCF level (0)
C       - MP2 level (1)
C       - CCSD level (2)
C       - CCSD(T) level (3)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER IOPT,nps
      REAL*8 T1(*),S1(*),T2(*),S2(*)
      REAL*8 DOO(*),DVV(*),DVO(*),DOV(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
      REAL*8 DMO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NFOO,A0,0,DOO,1)
      CALL XCOPY (NFVV,A0,0,DVV,1)
      IF (IOPT.LT.1) RETURN
C----------------------------------------------------------------------
C D(I,J) = D(I,J) - T(CD,IL) * L(CD,JL) = D(I,J) - T(CDL,I) * L(CDL,J)
C----------------------------------------------------------------------
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,T2,BUF1)
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,S2,BUF2)
      CALL CNTRCT ('C','N',NO,NO,NVVOT,-A1,BUF1,BUF2,A1,DOO,NREP)
C----------------------------------------------------------------------
C D(I,J) = D(I,J) - T(C,I) * L(C,J) 
C----------------------------------------------------------------------
      CALL CNTRCT ('C','N',NO,NO,NV,-A1,T1,S1,A1,DOO,NREP)
C----------------------------------------------------------------------
C D(A,B) = D(A,B) + L(AD,KL) * T(BD,KL) = H(A,C) + L(A,DKL) * T(B,DKL)
C----------------------------------------------------------------------
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF1)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,S2,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,A1,BUF1,BUF2,A1,DVV,NREP)
C----------------------------------------------------------------------
C D(A,B) = D(A,B) + L(A,K) * T(B,K)
C----------------------------------------------------------------------
      CALL CNTRCT ('N','C',NV,NV,NO,A1,T1,S1,A1,DVV,NREP)
C
      IF (CARITH) THEN
         CALL CONJUGA (NFOO,DOO,1)
         CALL CONJUGA (NFVV,DVV,1)
      ENDIF
C
C###########################################
C convert the density matrix to DIRAC format
C###########################################
      CALL CONVDM(BUF1,BUF2,DOO,DVV,DMO)
C
      IF (IOPT.LT.2) RETURN
C----------------------------------------------------------------------
C DVO(A,I) = L(A,I)
C----------------------------------------------------------------------
      CALL XCOPY (NDIMT1,S1,1,DVO,1)         
C----------------------------------------------------------------------
C DOV(A,I) = T(A,I)
C----------------------------------------------------------------------
      CALL XCOPY (NDIMT1,T1,1,DOV,1)         
C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) + T(AC,IK) * L(C,K)
C----------------------------------------------------------------------
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              T2,BUF1)
      CALL XGEMV ('N',NDIMT1,NDIMT1,A1,BUF1,M,S1,1,A1,DOV,1)
C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) - T(A,K) * L(C,K) * T(C,I)
C----------------------------------------------------------------------
      CALL CNTRCT ('C','N',NO,NO,NV,-A1,S1,T1,A0,BUF1,NREP)
      CALL CNTRCT ('N','N',NV,NO,NO,A1,T1,BUF1,A1,DOV,NREP)
C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) - T(A,K) * L(CD,LK) * T(CD,LI)
C----------------------------------------------------------------------
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,S2,BUF1)
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,T2,BUF2)
      CALL CNTRCT ('C','N',NO,NO,NVVOT,-A1,BUF1,BUF2,A0,BUF3,NREP)
      CALL CNTRCT ('N','N',NV,NO,NO,A1,T1,BUF3,A1,DOV,NREP)
C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) - L(CD,KL) * T(AD,KL) * T(D,I)
C----------------------------------------------------------------------
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,S2,BUF1)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,-A1,BUF1,BUF2,A0,BUF3,NREP)
      CALL CNTRCT ('N','N',NV,NO,NV,A1,BUF3,T1,A1,DOV,NREP)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CONVDM(BUF1,BUF2,DOO,DVV,dmo)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Convert the density matrix to DIRAC format
C
C
C---------------Routines called----------------------------------------
C
C     DCOPY
C     DOOTODIR
C     DVVTODIR
C---------------Last modified------------------------------------------
C
C     Author : Joost van Stralen - March 2002
C
C---------------Calling variables--------------------------------------
C
      REAL*8 DOO(*),DVV(*)
      REAL*8 BUF1(*),BUF2(*)
      REAL*8 DMO(*)
C
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
C
C     Reorder density matrix
C
C######################################################################
C  STEP 1:
C
C   OO       VO         OV      VV            Re          Im
C +-----+   +-----+   +------+   +-----+    +----------+----------+
C |     |   |     |   |      |   |     |    |          |          |
C |     | , |     | , |      | , |     | -> |          |          |
C +-----+   +-----+   +------+   +-----+    |          |          |
C                                           |          |          |
C                                           +----------+----------+
C                                           <---O+V---->
C                                            unbar + bar
C######################################################################
C
C##############
C     DOO
C##############
      IJ = 0
      DO IREP = 1, NREP
         DO J = 1, NO(IREP)
            DO I = 1, NO(IREP)
               IJ = IJ + 1
               IF (CARITH) THEN
                  BUF1(IJ) = DOO(RCW*IJ-1)
                  BUF2(IJ) = -DOO(RCW*IJ)
               ELSE
                  BUF1(IJ)= DOO(IJ)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C Put imaginary part after real part
C
      IF (CARITH) THEN
         CALL DCOPY(NFOO,BUF2,1,BUF1(NFOO+1),1)
      ELSE
      ENDIF
c     do i = 1, nfoo
c        print*,'real DOO ',i,buf1(i)
c        print*,'imag DOO ',i,buf1(i+nfoo)
c     enddo
C
c     print*,'laat oo deel dmo op 0'
      CALL DOOTODIR(dmo,BUF1,NO,nfoo,NREP,RCW)
C##############
C     DVV
C##############
      AB = 0
      DO IREP = 1, NREP
         DO B = 1, NV(IREP)
            DO A = 1, NV(IREP)
               AB = AB + 1
               IF (CARITH) THEN
                  BUF1(AB) = DVV(RCW*AB-1)
                  BUF2(AB) = -DVV(RCW*AB)
               ELSE
                  BUF1(AB)= DVV(AB)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C Put imaginary part after real part
C
      IF (CARITH) THEN
         CALL DCOPY(NFVV,BUF2,1,BUF1(NFVV+1),1)
      ELSE
      ENDIF
c     do i = 1, nfvv
c        print*,'real DVV ',i,buf1(i)
c        print*,'imag DVV ',i,buf1(i+nfvv)
c     enddo
C
c     print*,'laat vv deel dmo op 0'
      CALL DVVTODIR(dmo,BUF1,NV,nfvv,NREP,RCW)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCDENSZ(IOPT,NSP,IREPSPI,EPS,T1,S1,T2,S2,
     &                   BUF1,BUF2,BUF3,NBUF2,NBUF3,CBUF,
     &                   XVO,XOV,DOO,DVV,DVO,DOV,
     &                   EVO,EOV,D2VO,D2OV,DG,CC1,BB1,BB2,CONV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Menu driven routine. Depending on value of IOPT
C
C     Solve Z-vector equations at the
C       - SCF level (0)
C       - MP2 level (1)
C       - CCSD level (2)
C       - CCSD(T) level (3)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER IOPT
      REAL*8 T1(*),S1(*),T2(*),S2(*)
      REAL*8 DOO(*),DVO(*),DOV(*),DVV(*)
      REAL*8 D2VO(*),D2OV(*),DG(*)
      REAL*8 XVO(*),XOV(*)
      REAL*8 EVO(*),EOV(*)
      REAL*8 EPS(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*),CBUF(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
C
C---------------Local variables--------------------------------------
C
      CHARACTER*8 NAME
      COMPLEX*16 EPHASE(MAXOP)
      LOGICAL BDONE,USEDZ,RIGHT
C
C---------------Executable code--------------------------------------
C
      IF (IOPT.LT.1) RETURN
C
C     No relaxation if there are no virtuals within the right irrep
C
      IF (NFVO.EQ.0) RETURN
C
C     Calculate X-intermediate
C
      CALL XCOPY (NDIMX,A0,0,XVO,1)
C--------------------------------------------
C X(A,I) = X(A,I) + V(A,K,C,D) * TAU(C,D,I,K)
C--------------------------------------------
      CALL XCOPY (NDIMT2,S2,1,BUF1,1)
      CALL T1EQN2 (BUF1,BUF2,BUF3,NBUF3,XVO)
C
C------------------------------------
C X(A,I) = X(A,I) + W(CI,DA) * D(C,D)
C------------------------------------
      M = NFVV
      N = NFVO
      USEDZ = .FALSE.
      RIGHT = .FALSE.
      CALL SRT20D(NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &            JJVO,BUF2,NBUF2,XVO,DVV,USEDZ,RIGHT)
C
C In parallel runs we add all distributed contribution together here and
C continue on the master.
C
#if defined (VAR_MPI)
      IF(NMPROC.GT.1) THEN
C
C     LV 16-10-2002
C     The next statements might be combined to an MPI_ALLREDUCE
C     but I found the SGI implementation to give a runtime error.
C
        CALL MPI_REDUCE(XVO,BUF1,RCW*NDIMX,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,
     &    MPI_COMM_WORLD,IERR)
        CALL XCOPY(NDIMX,BUF1,1,XVO,1)
        CALL MPI_BCAST (XVO,RCW*NDIMX,
     &    MPI_DOUBLE_PRECISION,MASTER,
     &    MPI_COMM_WORLD,IERR)
      ENDIF
#endif
C--------------------------------------------
C X(A,I) = X(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------
      CALL XCOPY (NDIMT2,S2,1,BUF3,1)
      CALL T1EQN1 (BUF1,BUF2,BUF3,NBUF3,XVO)
C---------------------------------------------------------------------
C X(A,I) = X(A,I) + V(IK,AL) * D(K,L)
C---------------------------------------------------------------------
      CALL GETVOOO (BUF1)
      CALL SRT1ST4 (NREP,MULTB,LTR,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &              BUF1,BUF2)
      M = NFVO
      K = NFOO
      CALL XGEMV ('N',M,K,A1,BUF2,M,DOO,1,A1,XVO,1)
C
      CALL XSCAL(NDIMX,-A1,XVO,1)
      CALL XCOPY (NDIMX,XVO,1,XOV,1)
      IF (CARITH) CALL CONJUGA (NDIMX,XOV,1)
C
C     Solve Z-vector equations by calling the RPA routine
C
      NAME = 'MP2-GRAD'
      DO EOPER = 1, NEOPER
         EOFF = (EOPER-1)*NDIMX*RCW+1
         CALL RDPROP(IPRNT,NAMEE(EOPER),NSP,LFA,CBUF,IREPSPI,
     &               EPHASE(EOPER),EVO(EOFF),EOV(EOFF),DUM,DUM)
      ENDDO
C
C     Silence RPA routine when the print level is smaller than one
C
      IF (IPRNT.GE.1) THEN
         JPRNT = IPRNT
      ELSE
         JPRNT = -1
      ENDIF
C
      BDONE = .FALSE.
      MAXIT = MXITGR
      CALL RPAA (MAXIT,NEOPER,NAMEE,NAME,EPHASE,A1,BDONE,EPS,A0,
     &           EVO,EOV,XVO,XOV,DVO,D2VO,DOV,D2OV,DG,BUF1,BUF2,
     &           CC1,BB1,BB2,CONV,JPRNT,DEBUG,TIMING)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MP2ZVEC(T2,BUF1,BUF2,BUF3,NBUF3,XMO,DMO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Solve Z-vector equations at the MP2 level, using the AO based 
C     implementation.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T2(*)
      REAL*8 DMO(*)
      REAL*8 XMO(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "inpt.inc"
#include "symm.inc"
#include "param.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
C     +------------------------------------------------+
C     | Generate the T's in the right order for the    |
C     | contraction with the 3/4 transformed integrals |
C     +------------------------------------------------+
C
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              T2,BUF1)
      CALL SRTDIR (NREP,MULTB,MVO,J2VOVO,JJVO,JJVO,BUF1,BUF2)
      CALL SRTT1T2(BUF2,BUF1)
C
C     +---------------------------------------------------+
C     | Generate the tau dependent part of the Lagrangian |
C     +---------------------------------------------------+
C
      LFREE=NBUF3*RCW
      CALL TR5INI(XMO,BUF2,BUF1,BUF3,LFREE)
C
C     +-------------------------------------------------------+
C     | Calculate the ji' and b'a parts of the density matrix |
C     +-------------------------------------------------------+
C
      KFREE=1
      LFREE=NBUF3*RCW
      CALL CALCDJI(DMO,XMO,BUF3,KFREE,LFREE)
C
C     +-------------------------------------------------+
C     | Generate the D dependent part of the Lagrangian |
C     +-------------------------------------------------+
C
      KFREE=1 
      LFREE=NBUF3*RCW
      CALL MPGLVC(DMO,XMO,BUF3,KFREE,LFREE,1)
C
C     +----------------------------------------+
C     | Solve Z-vector equations by using XRPA |
C     +----------------------------------------+
C
      KFREE=1
      LFREE=NBUF3*RCW
      CALL PRPXLR_MPG(DMO,XMO,BUF3,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CONVXM(BUF1,BUF2,XVO,XMO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Convert x-matrix(Lagrangian) to DIRAC format
C
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Joost van Stralen - April 2002
C
C---------------Calling variables--------------------------------------
C
      REAL*8 XVO(*)
      REAL*8 BUF1(*),BUF2(*)
      REAL*8 XMO(*)
C
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
C
C###########################################
C STEP 1: Put imaginary part after real part
C###########################################
      AI = 0
      DO IREP = 1,NREP
         DO I = 1, NO(IREP)
            DO A = 1, NV(IREP)
               AI = AI + 1
               IF (CARITH) THEN
                  BUF1(AI) = XVO(AI*RCW-1)
                  BUF2(AI) = XVO(AI*RCW)
c                 print*,'real XVO',ai,BUF1(AI)
c                 print*,'imag XVO',ai,BUF2(AI)
               ELSE
                  BUF1(AI) = XVO(AI)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
      IF (CARITH) THEN
         CALL DCOPY(NFVO,BUF2,1,BUF1(NFVO+1),1)
      ELSE
      ENDIF
C######################################################################
C  STEP 2:
C
C           Re      Im
C     -  +-------+-------+
C   V |  |*******|*******|
C  u+b|  |*******|*******|
C     -  +-------+-------+
C        <---O--->
C         unbar+bar
C                     -------------->
C
C          Re        I         J       K
C       +--------+--------+--------+--------+
C     P |        |        |        |        |
C     O |        |        |        |        |
C     V |  ***   |  ***   |  ***   |  ***   |
C       +--------+--------+--------+--------+
C        P  O  V
C######################################################################
C
      CALL RELTODIR(XMO,BUF1,NFVO,RCW)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXPVAL (DOO,DVV,DVO,DOV,EOO,EVV,EVO,EOV,EPHASE,CEXPA)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates expectation value of an operator.
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
      REAL*8 DOO(*),DVV(*),DVO(*),DOV(*)
      REAL*8 EOO(*),EVV(*),EVO(*),EOV(*)
      COMPLEX*16 CEXPA,EPHASE
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 EXPA,EXPA_SUM
      REAL*8 DDOT,ZDOT
C
C---------------Executable code--------------------------------------
C
      EXPA = (0.0, 0.0)
      CALL XDOTU (EXPA,NFVO,EOV,1,DVO,1)
      EXPA_SUM = EXPA
      CALL XDOTU (EXPA,NFVO,EVO,1,DOV,1)
      EXPA_SUM = EXPA_SUM + EXPA
c     print*,' Relaxation contribution ',expa_sum
      CALL XDOTU (EXPA,NFOO,EOO,1,DOO,1)
c     print*,' OO contribution ',EXPA
      EXPA_SUM = EXPA_SUM + EXPA
      CALL XDOTU (EXPA,NFVV,EVV,1,DVV,1)
c     print*,' VV contribution ',EXPA
      EXPA_SUM = EXPA_SUM + EXPA
C
C     This is the correlation contribution, add the reference value
C
      II = 1
      DO IREP = 1, NREP
         CALL XDOTU (EXPA,NO(IREP),EOO(II),NO(IREP)+1,A1,0)
         EXPA_SUM = EXPA_SUM + EXPA
         II = II + NO(IREP)*NO(IREP)*RCW
      ENDDO
C
      CEXPA = EPHASE * EXPA_SUM
C
      RETURN
      END
