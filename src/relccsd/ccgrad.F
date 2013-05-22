!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
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
      use interface_to_mpi
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
        CALL XCOPY(NDIMX,XVO,1,BUF1,1)
        CALL DZERO(XVO,RCW*NDIMX)
        call interface_MPI_allREDUCE_r1_work_f77(BUF1,XVO,RCW*NDIMX,
     &                                           op_mpi_sum,
     &                                           global_communicator)
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
C     Called from: CCFOPR/ccdriv.F
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
CMI 
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
      CALL FLSHFO(IW)

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
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE  MP2_NATORB(DMO,NOPRNT,NOZG)
C***********************************************************************
C
C     Calculate MP2 nat. orb. and nat. orb. occ. num.
C
C     Called from: CCFOPR/ccdriv.F
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Stefan Knecht and Hans Joergen Aa. Jensen, Feb. 2008
C
C     Last modifications: 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbmp2no.h"
#include "infpar.h"
C
      DIMENSION DMO(NORBT,NORBT,NZ)
      LOGICAL NOZG
      real(8), allocatable :: CMO_NAT(:), OCC_NAT(:)
C
C---------------Local variables--------------------------------------
C
      CHARACTER*8 LABEL
C
      CALL QENTER('MP2_NATORB')
C
C     memory allocation
C     -----------------
      allocate(CMO_NAT(NORBT*NESHT*NZ))
      allocate(OCC_NAT(NESHT))
C
C     initialize
C
      CALL DZERO(CMO_NAT,NORBT*NESHT*NZ)
      CALL DZERO(OCC_NAT,NESHT)
C
C     build diagonal according to Fock matrix
C     ---------------------------------------
      CALL ADD_TWO_FOCK(DMO,NOZG,NOPRNT)
C
C     calculate NO occ. num. 
C     ----------------------
      CALL MP2_NATORB_GEN(DMO,CMO_NAT,OCC_NAT,NOPRNT)
C
C     save nat. orb. occ. num. and transform. matrix to file 
C     ------------------------------------------------------
      IF( MYTID .eq. MPARID )THEN
C
C       ... KRMCOLD file
C
C       sk + hjaaj: we use KRMCOLD and not KRMCSCF for the MP2 nat.orb.
C       coefficients, because the KRMC module starts by deleting the
C       KRMCSCF file !!! (And it uses CMO from KRMCOLD as first priority
C       for starting orbitals.)
C
        CALL OPNFIL(LUKRMC,'KRMCOLD','UNKNOWN','MP2_NA')
        CALL NEWLAB('*KRMCMP2',LUKRMC,LUPRI)
C
        CALL WRTKRMC(LUKRMC,'MP2NATOB',CMO_NAT,NORBT*NESHT*NZ)
        CALL WRTKRMC(LUKRMC,'MP2NOOCC',OCC_NAT,NESHT)
        CLOSE(LUKRMC,STATUS='KEEP')
C
      END IF
C
C     memory deallocation
C     -------------------
      deallocate(OCC_NAT)
      deallocate(CMO_NAT)
C
      CALL QEXIT('MP2_NATORB')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MP2_NATORB_GEN(DMO,CMO_NAT,OCC_NAT,NOPRNT)
C***********************************************************************
C
C     diagonalize second order density matrix DMO
C     and calculate nat. orb. occ. num.
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Stefan Knecht and Hans Joergen Aa. Jensen, Feb. 2008
C
C     Last modifications: 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
C
      DIMENSION OCC_NAT(*), CMO_NAT(NORBT,NESHT,NZ)
      DIMENSION DMO(NORBT,NORBT,*)
      real(8), allocatable :: scratch(:)
C
      IF ( NOPRNT .ge. 5 )THEN
         CALL HEADER(
     &               'MP2_NATORB_GEN: DMO :',-1)
         CALL PRQMAT(DMO,
     &               NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF

      IOFF2 = 1
C
      DO IFSYM = 1, NFSYM  
C
        IOFF1 = IORB(IFSYM) + NPSH(IFSYM) + 1
C
        IF(DOJACO) THEN
           lfree = ((NESH(IFSYM)*(NESH(IFSYM)+1))/2)+2*NESH(IFSYM)+16
           allocate(scratch(lfree))
           CALL RSJACO(NORBT,NESH(IFSYM),NESH(IFSYM),DMO(IOFF1,IOFF1,1),
     &                 OCC_NAT(IOFF2),1,-1,0,CMO_NAT(IOFF1,IOFF2,1),
     &                 scratch,LFREE)
           deallocate(scratch)
        ELSE
           CALL QDIAG90(NZ,NESH(IFSYM),DMO(IOFF1,IOFF1,1),NORBT,NORBT,
     &                OCC_NAT(IOFF2),1,
     &                CMO_NAT(IOFF1,IOFF2,1),NORBT,NESH(IFSYM),IERR)
CSK        WRITE(LUPRI,*) ' Fermionic symmetry ', IFSYM
CSK        CALL HEADER(
CSK  &     'MP2_NATORB_GEN: Eigenvalues unsorted:',-1)
CSK        WRITE(LUPRI,'(I5,F25.10)')
CSK  &     (J,OCC_NAT(IOFF2-1+J),J=1,NESH(IFSYM)) 
C          reorder eigenvalues and eigenvectors: descending order ==>
C          IORDER = - 1
           IORDER = - 1
           CALL QORDER3(NESH(IFSYM),OCC_NAT(IOFF2),
     &                  CMO_NAT(IOFF1,IOFF2,1),
     &                  NORBT,NESH(IFSYM),NESH(IFSYM),NZ,IORDER)
        ENDIF
        IOFF2 = IOFF2 + NESH(IFSYM)
      END DO
C
C     debug print section
C
      IOFF2 = 1
      IF ( NOPRNT .ge. 5 ) THEN
         DO IFSYM = 1, NFSYM
           IOFF1 = IORB(IFSYM) + NPSH(IFSYM) + 1
           CALL HEADER(
     &     'MP2_NATORB_GEN: Eigenvectors :',-1)
           WRITE(LUPRI,*) ' fermionic irrep', IFSYM
           CALL PRQMAT(CMO_NAT(IOFF1,IOFF2,1),
     &                 NORBT,NESH(IFSYM),NORBT,NESH(IFSYM),
     &                 NZ,IPQTOQ(1,0),LUPRI)
           IOFF2 = IOFF2 + NESH(IFSYM)
         END DO
      ENDIF
C
C     MP2 - nat. orb. occ. num.
C
      IOFF2 = 1
C
      IF ( NOPRNT .ge. 3 ) THEN
        DO I = 1, NFSYM
           NOCCI = NESH(I)
           IF (NOCCI .EQ. 0) THEN
              WRITE(LUPRI,9001) I
           ELSE
              OCCSUM = DSUM(NOCCI,OCC_NAT(IOFF2),1)
              WRITE(LUPRI,9002) I, (OCC_NAT(IOFF2-1+J),J=1,NESH(I))
              WRITE(LUPRI,9003) OCCSUM
           END IF
           IOFF2 = IOFF2 + NESH(I)
        END DO
      END IF
C
 9001 FORMAT(/'   Symmetry',I3,//,'   No occupied orbitals')
 9002 FORMAT(/'   Natural orbital occupation numbers, symmetry',
     &            I3,//,(5F14.9))
 9003 FORMAT(/'   Sum =',T15,F14.9)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADD_TWO_FOCK(DMO,NOZG,NOPRNT)
C***********************************************************************
C
C     build proper diagonal part in second order density matrix DMO
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Stefan Knecht, Feb. 2008
C
C     Last modifications: 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "../moltra/dcbtr3.h"
C
      DIMENSION DMO(NORBT,NORBT,*)
      LOGICAL NOZG
C
      D2 = 2.0D0
C
C     add "alpha" and "beta" by scaling with 2.0
C
      CALL DSCAL(NZ*NORBT*NORBT,D2,DMO,1)
C
      IF ( NOPRNT .ge. 10 )THEN
         CALL HEADER(
     &               'ADD_TWO_FOCK input: DMO :',-1)
         CALL PRQMAT(DMO,
     &               NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      IOFF2 = 1
C
      DO IFSYM = 1, NFSYM  
        IOFF1 = IORB(IFSYM) + NPSH(IFSYM) + 1
C
!       WRITE(LUPRI,*) ' NDMOQC(2,IFSYM,1)',NDMOQC(2,IFSYM,1)
C       zero 'wrong' offdiagonal terms for frozen orbs
        IF( .NOT. NOZG )THEN
          DO IZ = 1, NZ
            CALL DZERO(DMO(IOFF1,IOFF1,IZ),NDMOQC(2,IFSYM,1)*NORBT)
          END DO
        END IF
C
        IOCC_A = 0
        DO I = 1, NOCC(IFSYM)
           DMO(IOFF1+IOCC_A,IOFF1+IOCC_A,1) = 
     &     DMO(IOFF1+IOCC_A,IOFF1+IOCC_A,1) + D2
           IOCC_A = IOCC_A + 1
        END DO
      END DO
C
      IF ( NOPRNT .ge. 10 )THEN
         CALL HEADER(
     &               'ADD_TWO_FOCK output: DMO :',-1)
         CALL PRQMAT(DMO,
     &               NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_CC_DENSITY (WF,BASIS,DMAT,IERR)
C
C---------------Description--------------------------------------------
C
C     Retrieve relaxed density matrix for processing by other modules
C     of DIRAC (currenty only first-order properties)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Common Blocks--------------------------------------
C
#include "implicit.h"
#include "files.inc"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
C
C---------------Calling variables--------------------------------------
C
      PARAMETER (HALF=0.5D0)
      CHARACTER*4 WF                     ! Wave function type
      CHARACTER*2 BASIS                  ! Representation
      REAL*8 DMAT(*)                     ! Density matrix
      INTEGER IERR                       ! Error code
C
C---------------Local variables--------------------------------------
C
      REAL*8, ALLOCATABLE:: DMO(:,:,:)   ! Density matrix in MO basis
      LOGICAL FND,FNDLAB,ADD_2_FOCK
      CHARACTER*8 LABEL
      CHARACTER*6 FILE_NAME
C
C---------------Executable code--------------------------------------
 
      ADD_2_FOCK = .FALSE.
      SELECT CASE (WF)
         CASE ("MP2 ")
            WRITE (LUPRI,*) " Retrieving MP2 density matrix from file"
            LABEL        = 'MP2 DENS'
            FILE_NAME    = 'CCDENS'
            ADD_2_FOCK   = .TRUE.
            IERR         = 0
         CASE ("MC  ")
            WRITE (LUPRI,*) " Retrieving MCSCF density matrix from file"
            LABEL        = 'MC  DENS'
            FILE_NAME    = 'MCDENS'
            ADD_2_FOCK   = .FALSE.
            IERR         = 0
            CALL QUIT('MCSCF density read not fully implemented yet!')
         CASE ("CI  ")
            WRITE (LUPRI,*) " Retrieving CI density matrix from file"
            LABEL        = 'CI  DENS'
            FILE_NAME    = 'MCDENS'
            ADD_2_FOCK   = .FALSE.
            IERR         = 0
            CALL QUIT('CI density read not fully implemented yet!')
         CASE ("HF  ")
            WRITE (LUPRI,*) " Generating Hartree-Fock density matrix"
            FILE_NAME    = 'CCDENS'
            ADD_2_FOCK   = .TRUE.
            IERR         = 0
         CASE DEFAULT
            WRITE (LUPRI,*) " Unknown or unimplemented density: ",WF
            WRITE (LUPRI,*) " Will use Hartree-Fock density matrix"
            IERR = 1
      END SELECT

C     Read the appropriate correlated density from file

      ALLOCATE (DMO(NORBT,NORBT,NZ),STAT=IERR)
      IF (IERR.NE.0) THEN
        CALL QUIT('Error in DMO allocation !')
      ENDIF
      CALL OPNFIL(MCCDENS,FILE_NAME,'OLD','GET_CC')
      FND =  FNDLAB(LABEL,MCCDENS)
      CALL READT(MCCDENS,NORBT*NORBT*NZ,DMO)
      CLOSE(Unit=MCCDENS,Status='KEEP')

      
C     Add the closed shell Hartree-Fock contribution
      IF(ADD_2_FOCK)THEN
C       The add_two routine scales by 2, scale back to the orginal afterwards
        CALL ADD_TWO_FOCK(DMO,.TRUE.,0)
        CALL DSCAL(NZ*NORBT*NORBT,HALF,DMO,1)
      END IF

      SELECT CASE (BASIS)
         CASE ("MO")
            CALL DCOPY (N2ORBTQ,DMO,1,DMAT,1)
         CASE ("AO")
            CALL DENSMO_TO_DENSAO (DMO,DMAT)
         CASE DEFAULT
            WRITE (LUPRI,*) " Unknown or unimplemented basis: ",BASIS
            WRITE (LUPRI,*) " Use density matrix in AO basis "
            IERR = IERR + 2
            CALL DENSMO_TO_DENSAO (DMO,DMAT)
      END SELECT

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE STORE_CC_DENSITY (WF,DMO)
C
C---------------Description--------------------------------------------
C
C     Retrieve relaxed density matrix for processing by other modules
C     of DIRAC (currenty only first-order properties)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Common Blocks--------------------------------------
C
#include "implicit.h"
#include "files.inc"
#include "priunit.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*4 WF                     ! Wave function type
      REAL*8 DMO(*)                      ! Density matrix (MO-basis)
C
C---------------Local variables--------------------------------------
C
      CHARACTER*8 LABEL
C---------------Executable code--------------------------------------
 
      WRITE (LUPRI,*) "Storing the ",WF," density matrix on file"
      LABEL = WF//'DENS'

C     Write the correlated density to file

      CALL OPNFIL(MCCDENS,'CCDENS','UNKNOWN','STORE_')
      CALL NEWLAB(LABEL,MCCDENS,LUPRI)
      CALL WRITT(MCCDENS,NORBT*NORBT*NZ,DMO)
      CLOSE(Unit=MCCDENS,Status='KEEP')

      END
