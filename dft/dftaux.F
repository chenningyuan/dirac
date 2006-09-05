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
C /* Deck dft_linres */
      SUBROUTINE DFT_LINRES(DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,GAO,GAO1,
     &                      RHO,RHOGRD,RHG,WGHT,DOGGA,BUF)
C
C     T. Helgaker sep 99/oct 00
C     adapted to DIRAC by T.Saue May 3 2001
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,DP5 = 0.5D0,D2=2.0D0)
C
#include "inforb.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dgroup.h"
#else
#include "wrkrsp.h"
#include "dftinf.h"
#endif
#include "dftcom.h"
#include "cdft.h"
#include "functionals.h"
C
      INTEGER A, B
      LOGICAL DOGGA, TRPLET
      DIMENSION DTRMAT(NBAST,NBAST,NZ,NXCMAT),
     &          EXCMAT(NBAST,NBAST,NZ,NXCMAT),
     &          ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &          GAO(NBAST), GAO1(NBAST,3),
     &          RHG(3),VT(3), VXC(DERVS2), BUF(*)
      DIMENSION B3(3)
C
      DO IM = 1,NXCMAT
      IREPOP = ISYMOP(IM)-1
      IF(IHRMOP(IM).EQ.1) THEN
C
C       Get perturbed density
C       =====================
C
        CALL GETRHO(B0,IREPOP,GAO,DTRMAT(1,1,1,IM),BUF)
C
C       ***************
C       ***** LDA *****
C       ***************
C
        IF(.NOT.DOGGA) THEN
          IF(DABS(B0).GT.DFTHRL) THEN
C
C           Calculated VT
C
            CALL DFTPOT1(VXC, WGHT,RHO, 0D0)
            SXC = VXC(FRR)*B0
C
C           Linear transformation
C
            CALL DOEXC1(SXC,IREPOP,EXCMAT(1,1,1,IM),GAO)
          END IF
        ELSE
C
C         ***************
C         ***** GGA *****
C         ***************
C
          CALL GETGRHO(B3,IREPOP,GAO,GAO1,DTRMAT(1,1,1,IM),BUF)
          BMAX = MAX(ABS(B0),ABS(B3(1)),ABS(B3(2)),ABS(B3(3)))
C
          IF (BMAX.GT.DFTHRL) THEN
            CALL DFTPOT1(VXC,WGHT,RHO,RHOGRD)
            BR = (B3(1)*RHG(1)+ B3(2)*RHG(2)+ B3(3)*RHG(3))/RHOGRD
            FACR = VXC(FRZ)*B0 + (VXC(FZZ)-VXC(FZ0)/RHOGRD)*BR
            SXC = VXC(FRR)*B0 + VXC(FRZ)*BR
            DO I = 1,3
              VT(I) = (VXC(FZ0)*B3(I)+FACR*RHG(I))/RHOGRD
            ENDDO

C           Linear transformation
C
            CALL DOEXC2(SXC,VT,IREPOP,EXCMAT(1,1,1,IM),GAO,GAO1)
CTROND            CALL DOLRG(VT,IREPOP,EXCMAT(1,1,1,IM),GAO,GAO1)
          END IF
        END IF
#if defined (MOD_CDFT)
      ELSEIF(IHRMOP(IM).EQ.-1) THEN
        IF(SDFT) THEN
          CALL SDFT_LINEAR_TRANS(DTRMAT(1,1,1,IM),
     &                           EXCMAT(1,1,1,IM),
     &                           IREPOP,
     &                           GAO,RHO,WGHT,DOGGA,BUF)
        ENDIF
      ELSE
        CALL QUIT('DTRMAT has to be hermitian or antihermitian')
#endif
      ENDIF
      ENDDO
      RETURN
      END
      SUBROUTINE DFT_QUADRES(DTRMAT,EXCMAT,ISYMOP,IHRMOP,GAO,GAO1,
     &                       RHO,RHOGRD,RHG,WGHT,DOGGA,BUF)
C***********************************************************************
C
C     Based on DFT_LINRES
C
C     Driver for DFT contributions to quadratic response.
C
C     On input:
C     DTRMAT(*,*,*,1) --- B density matrix
C     DTRMAT(*,*,*,2) --- C density matrix
C     DTRMAT(*,*,*,3) --- BC density matrix
C     GAO             --- 
C
C     On output:
C     EXCMAT(*,*,*,1) --- G_B matrix
C     EXCMAT(*,*,*,2) --- G_C matrix
C     EXCMAT(*,*,*,3) --- H_BC matrix
C
C     tsaue & johhe Aug. 2006
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,DP5 = 0.5D0,D2=2.0D0)
C
#include "inforb.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dgroup.h"
#else
#include "wrkrsp.h"
#include "dftinf.h"
#endif
#include "dftcom.h"
#include "functionals.h"
C
      INTEGER A, B
      LOGICAL DOGGA, TRPLET
      DIMENSION DTRMAT(NBAST,NBAST,NZ,3),
     &          EXCMAT(NBAST,NBAST,NZ,3),
     &          ISYMOP(3),IHRMOP(3),
     &          GAO(NBAST), GAO1(NBAST,3), 
     &          RHG(3),VT(3), VXC(DERVS3), BUF(*)
!johhe Check DERVS3
      DIMENSION B3(3)
C
C     Get perturbed densities
C     =======================
C
      IREPBOP = ISYMOP(1)-1
      IF(IHRMOP(1).EQ.1) THEN
         CALL GETRHO(B0,IREPBOP,GAO,DTRMAT(1,1,1,1),BUF)
      ELSE
         CALL QUIT('DFTQR program general case!')
      END IF
C
      IREPCOP = ISYMOP(2)-1
      IF(IHRMOP(2).EQ.1) THEN
         CALL GETRHO(C0,IREPCOP,GAO,DTRMAT(1,1,1,2),BUF)
      ELSE
         CALL QUIT('DFTQR program general case!')
      END IF
C
      IREPBCOP = ISYMOP(3)-1
      IF(IHRMOP(3).EQ.1) THEN
         CALL GETRHO(BC0,IREPBCOP,GAO,DTRMAT(1,1,1,3),BUF)
      ELSE
         CALL QUIT('DFTQR program general case!')
      END IF
C
C       ***************
C       ***** LDA *****
C       ***************
C
      IF(.NOT.DOGGA) THEN
C
C     Possibly check sizes of B0, C0, and BC0, if all small, skip this.
         CALL DFTPOT2(VXC,WGHT,RHO,0D0)
C
C     Constructing G_B matrix
         IF(DABS(B0).GT.DFTHRL) THEN 
C
C           Calculated VT
C
            SXC = VXC(FRR)*B0
C
C           Linear transformation
C  
            CALL DOEXC1(SXC,IREPBOP,EXCMAT(1,1,1,1),GAO)
         END IF
C
C     Constructing G_C matrix
         IF(DABS(C0).GT.DFTHRL) THEN 
C
C           Calculated VT
C
            SXC = VXC(FRR)*C0
C
C           Linear transformation
C  
            CALL DOEXC1(SXC,IREPCOP,EXCMAT(1,1,1,2),GAO)
         END IF
C
C     Constructing H_BC matrix
         IF(DABS(BC0).GT.DFTHRL) THEN 
C
C           Calculated VT
C
            SXC = VXC(FRRR)*B0*C0 + VXC(FRR)*BC0
C
C           Linear transformation
C  
            CALL DOEXC1(SXC,IREPBCOP,EXCMAT(1,1,1,3),GAO)
         END IF
      ELSE  
C
C         ***************
C         ***** GGA *****
C         ***************
C
C     Write me!
C
      END IF
!      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dolrg */
      SUBROUTINE DOLRG(VT,IREPOP,EXCMAT,GAO,GAO1)
C***********************************************************************
C
C     Construct Kohn-Sham matrix
C
C     Written by T.Saue and T.Helgaker and P.Salek
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),GAO1(NBAST,3),EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),GAO1(NTBAS(0),3),EXCMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION VT(0:3)
C
#include "ibtfun.h"
C
      DO 10 IRP1 = 0,NBSYM-1
        DO 20 IRP2 = 0,NBSYM-1
          IF(IBTXOR(IRP1,IRP2).NE.IREPOP) GOTO 20
          DO 30 IC = 1,MC
            NB1 = NBBAS(IRP1,IC)
            NB2 = NBBAS(IRP2,IC)
            JB1 = IBBAS(IRP1,IC)
            JB2 = IBBAS(IRP2,IC)
            DO J = 1, NB2
               G0 = GAO(JB2+J)
               GX = GAO1(JB2+J,1)
               GY = GAO1(JB2+J,2)
               GZ = GAO1(JB2+J,3)
               DO I = 1, NB1
                  A0 = G0*GAO(JB1+I)
                  AX = GX*GAO(JB1+I) + G0*GAO1(JB1+I,1)
                  AY = GY*GAO(JB1+I) + G0*GAO1(JB1+I,2)
                  AZ = GZ*GAO(JB1+I) + G0*GAO1(JB1+I,3)
                  EXCMAT(JB1+I,JB2+J) = EXCMAT(JB1+I,JB2+J)
     &                 + VT(0)*A0 + VT(1)*AX
     &                 + VT(2)*AY + VT(3)*AZ
              END DO
            ENDDO
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      END
c
      SUBROUTINE FORTWRT(str, slength)
#include "priunit.h"
      INTEGER unit, slength
      CHARACTER str(slength)
      WRITE(LUPRI, *) str
      END

      SUBROUTINE DFTSETCAM(W,BETA)
#include "implicit.h"
#include "dcbham.h"
      HFXMU = W
      HFXATT = BETA 
      END
