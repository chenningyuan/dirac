C /* Deck dftrpa */
      SUBROUTINE DFTRPA(DMAT,DTRMAT,EXCMAT,NXCMAT,IREPOP,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Prepare for DFT response....
C
C     Written byu T.Saue May 11 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2=2.0D0)
C
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
        DIMENSION DMAT(N2BBASX,NZ),DTRMAT(N2BBASX,NZ,NXCMAT),
     &            EXCMAT(N2BBASX,NZ,NXCMAT),WORK(*)
CTROND        CALL DCOPY(N2BBASX,DMAT,1,CMO,1)           
CTROND        IOFF = 1
CTROND        DO ISHELL = 1, NOPEN
CTROND          IOFF = IOFF + N2BBASXQ
CTROND          CALL DAXPY(N2BBASX,DF(ISHELL),DMAT(IOFF),1,CMO,1)
CTROND        ENDDO
CTROND
C       Insert half-phases
        IF(NZ.LT.4) THEN
           DO IM = 1,NXCMAT
            DO IZ = 1,NZ
              IQ = IPQTOQ(IZ,IREPOP)
              CALL Q2BPHASE('D',IQ,1,DTRMAT(1,IZ,IM))
              CALL Q2BPHASE('F',IQ,1,EXCMAT(1,IZ,IM))
            ENDDO
          ENDDO
        ENDIF
        CALL DSCAL(N2BBASX,D2,DMAT,1)
        CALL DSCAL(N2BBASXQ*NXCMAT,D2,DTRMAT,1)
#if defined (VAR_MPI)
        IF(PARCAL) THEN
CTROND: This must be fixed !
          CALL QUIT('Not yet implemented !')
          IRNTYP = 7
          CALL HER_PARDRV(WORK(KFREE),LFREE,FOCK,CMO,NFMAT,
     &          IDUM,IDUM,
     &          IRNTYP,IDUM,IDUM,.TRUE.,.TRUE.,.FALSE.,
     &          .FALSE.,.FALSE.,.FALSE.,IDUM,IDUM,IDUM,
     &          IDUM,DUM,DUM,DUM,DUM,DUM)
        ELSE
#endif
          CALL DFTEXC(DMAT,DTRMAT,EXCMAT,NXCMAT,IREPOP,
     &               .FALSE.,.FALSE.,.FALSE.,.TRUE.,
     &               WORK(KFREE),LFREE,IPRINT)
#if defined (VAR_MPI)
        ENDIF
#endif
C       Insert half-phases
        IF(NZ.LT.4) THEN
           DO IM = 1,NXCMAT
            DO IZ = 1,NZ
              IQ = IPQTOQ(IZ,IREPOP)
              CALL Q2BPHASE('D',IQ,1,DTRMAT(1,IZ,IM))
              CALL Q2BPHASE('F',IQ,1,EXCMAT(1,IZ,IM))
            ENDDO
          ENDDO
        ENDIF
      RETURN
      END       
C /* Deck dftltr */
      SUBROUTINE DFTLTR(DTRMAT,EXCMAT,NXCMAT,IREPOP,GAO,GAO1,RHO,
     &                  RHOGRD,RHG,WGHT,DOGGA,BUF)
C
C     T. Helgaker sep 99/oct 00
C     adapted to DIRAC by T.Saue May 3 2001
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,DP5 = 0.5D0)
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
      DIMENSION DTRMAT(NBAST,NBAST,NZ,NXCMAT),
     &          EXCMAT(NBAST,NBAST,NZ,NXCMAT),
     &          GAO(NBAST), GAO1(NBAST,3), 
     &          RHG(3),VT(0:3), VXC(DERVS2), BUF(*)
      DIMENSION B3(3)
C
      DO IM = 1,NXCMAT      
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
        IF (.NOT.DOGGA) THEN
CTROND           IF (DABS(B0).GT.DFTHRL) THEN 
           IF (DABS(B0).GT.1.0D-10) THEN 
C
C             Calculated VT
C
              CALL DFTPOT1(VXC, WGHT,RHO, 0D0)
              VT(0) = VXC(FRR)*B0
C
C             Linear transformation
C  
              CALL DOEXC1(VT(0),IREPOP,EXCMAT(1,1,1,IM),GAO)
           END IF
        ELSE  
C
C          ***************
C          ***** GGA *****
C          ***************
C
           CALL GETGRHO(B3,IREPOP,GAO,GAO1,DTRMAT(1,1,1,IM),BUF)
           BMAX = MAX(ABS(B0),ABS(B3(1)),ABS(B3(2)),ABS(B3(3)))
C          
           
           IF (.TRUE..OR.BMAX.GT.DFTHRL) THEN
              CALL DFTPOT1(VXC, WGHT,RHO, RHOGRD)
              BR = (B3(1)*RHG(1)+ B3(2)*RHG(2)+ B3(3)*RHG(3))/RHOGRD
              FACR = VXC(FRZ)*B0 + (VXC(FZZ)-VXC(FZ0)/RHOGRD)*BR
              VT(0) = VXC(FRR)*B0 + VXC(FRZ)*BR
              VT(1) = (RHG(1)*FACR + VXC(FZ0)*B3(1))/RHOGRD
              VT(2) = (RHG(2)*FACR + VXC(FZ0)*B3(2))/RHOGRD
              VT(3) = (RHG(3)*FACR + VXC(FZ0)*B3(3))/RHOGRD

C             Linear transformation
C
              CALL DOLRG(VT,IREPOP,EXCMAT(1,1,1,IM),GAO,GAO1)
         END IF
      END IF
      ENDDO
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
      MC = 2
      IF(LEVYLE) MC = 1      
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
