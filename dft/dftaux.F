C /* Deck dftrpa */
      SUBROUTINE DFTRPA(DMAT,DTRMAT,EXCMAT,NXCMAT,IREPOP,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Prepare for DFT respons....
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
      SUBROUTINE DFTLTR(DTRMAT,EXCMAT,NXCMAT,IREPOP,GAO,GAO1,RHO,RHO13,
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
     &          RHG(3),VXC(DERVS2), BUF(*)
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
              VT = VXC(FRR)*B0
C
C             Linear transformation
C  
              CALL DOEXC1(VT,IREPOP,EXCMAT(1,1,1,IM),GAO)
           END IF
        ELSE  
C
C          ***************
C          ***** GGA *****
C          ***************
CTRONDC
CTRONDC          B0, BX=B3(1), BY=B3(2), BZ=B3(3), BMAX
CTRONDC
CTRONDC             B3 = GAO1'*DTGAO
CTROND              CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D0,B3,1)
CTRONDC             DTGAO= DMAT'*GAO
CTROND            CALL DGEMV('T',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
CTRONDC             B3 = B3 + GAO1'*DTGAO
CTROND              CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D1,B3,1)
CTROND              BMAX = MAX(DABS(B0),DABS(B3(1)),DABS(B3(2)),DABS(B3(3)))
CTRONDC
CTROND           IF (BMAX.GT.DFTHRL) THEN
CTRONDC
CTRONDC             ZNV, FZ0, FRR, FRZ, FZZ
CTRONDC
CTROND              IF (DODRC) CALL V1DRC(VDRC0,FRRDRC,RHO,RHO13)
CTROND              IF (DOBCK) CALL V1BCK(FR0BCK,FZ0BCK,FRRBCK,FRZBCK,
CTROND     &                              FZZBCK,RHO,RHOGRD)
CTROND              IF (DOVWN) THEN
CTROND                 IF (.NOT.TRPLET) THEN
CTROND                    CALL V1VWN(VVWN0,FRRVWN,RHO,RHO13)
CTROND                 ELSE
CTROND                    CALL VTVWN(FRRVWN,RHO,RHO13)
CTROND                    FRRVWN = DP5*FRRVWN 
CTROND                 END IF
CTROND              END IF 
CTROND              IF (DOLYP) THEN
CTROND                 RHOA = DP5*RHO
CTROND                 RHGA = (DP5*RHOGRD)**2
CTRONDC                CALL V1LYP(FR0LYP,FZ0LYP,FRRLYP,FRZLYP,
CTRONDC    &                      FZZLYP,RHO,RHOGRD)
CTROND                 CALL GLYPCO(DF1000,DF0100,DF0010,DF0001,
CTROND     &                       DF00001,RHO,RHO13,RHOGRD,.TRUE.)
CTROND                 CALL VTLYP (DF2000,DF0200,DF1100,DF1010,
CTROND     &                       DF0101,DF1001,DF0110,DF10001,
CTROND     &                       DF01001,RHOA,RHOA,RHGA,RHGA,RHGA)
CTROND                 IF (.NOT.TRPLET) THEN
CTROND                    FZ0LYP = DP5*(DF0010 + DF00001)*RHOGRD
CTROND                    FRRLYP = DP5*(DF2000 + DF1100)
CTROND                    FRZLYP = DP5*(DF1010 + DF1001+DF10001)*RHOGRD
CTROND                    FZZLYP = FZ0LYP/RHOGRD 
CTROND                 ELSE
CTROND                    FZ0LYP = DP5*(DF0010 - DF00001)*RHOGRD
CTROND                    FRRLYP = DP5*(DF2000 - DF1100)
CTROND                    FRZLYP = DP5*(DF1010 - DF1001)*RHOGRD
CTROND                    FZZLYP = FZ0LYP/RHOGRD 
CTROND                 END IF
CTROND              ELSE 
CTROND                 FZ0LYP = 0D0
CTROND                 FRRLYP = 0D0
CTROND                 FRZLYP = 0D0
CTROND                 FZZLYP = 0D0
CTROND              END IF  
CTRONDC
CTROND              ZNV = D1/RHOGRD
CTROND              FZ0 = ZNV*(WBCK*FZ0BCK + WLYP*FZ0LYP)
CTROND              FRR = WDRC*FRRDRC + WVWN*FRRVWN 
CTROND     &            + WBCK*FRRBCK + WLYP*FRRLYP 
CTROND              FRZ = WBCK*FRZBCK + WLYP*FRZLYP 
CTROND              FZZ = WBCK*FZZBCK + WLYP*FZZLYP 
CTRONDC 
CTROND              RX = ZNV*RHG(1) 
CTROND              RY = ZNV*RHG(2)
CTROND              RZ = ZNV*RHG(3)
CTRONDC
CTRONDC             Linear transformation
CTRONDC
CTROND                 BR = B3(1)*RX + B3(2)*RY + B3(3)*RZ
CTROND                 FAC0 = FRR*B0 + FRZ*BR
CTROND                 FACR = FRZ*B0 + FZZ*BR
CTROND                 IF (NOSYM) THEN
CTROND                    DO I = 1, NBAST
CTROND                       G0 = GAO(I)
CTROND                       GX = GAO1(I,1)
CTROND                       GY = GAO1(I,2)
CTROND                       GZ = GAO1(I,3)
CTROND                       DO J = 1, I 
CTROND                          A0 = G0*GAO(J)
CTROND                          AX = GX*GAO(J) + G0*GAO1(J,1)
CTROND                          AY = GY*GAO(J) + G0*GAO1(J,2)
CTROND                          AZ = GZ*GAO(J) + G0*GAO1(J,3)
CTROND                          AR = AX*RX + AY*RY + AZ*RZ
CTROND                          AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) - AR*BR
CTROND                          EXCMAT(J,I) = EXCMAT(J,I)
CTROND     &                                +FAC0*A0+FACR*AR+FZ0*AB
CTROND                         END DO
CTROND                      END DO
CTROND                 ELSE
CTROND                    DO ISYM = 1, NSYM
CTROND                       ISTR = IBAS(ISYM) + 1
CTROND                       IEND = IBAS(ISYM) + NBAS(ISYM)
CTROND                       JSYM = MULD2H(ISYM,KSYMOP)
CTROND                       IF (ISYM.GE.JSYM) THEN
CTROND                          JSTR = IBAS(JSYM) + 1
CTROND                          JEND = IBAS(JSYM) + NBAS(JSYM)
CTROND                          DO I = ISTR, IEND
CTROND                           G0 = GAO(I)
CTROND                           GX = GAO1(I,1)
CTROND                           GY = GAO1(I,2)
CTROND                           GZ = GAO1(I,3)
CTROND                           DO J = JSTR, MIN(I,JEND) 
CTROND                              A0 = G0*GAO(J)
CTROND                              AX = GX*GAO(J) + G0*GAO1(J,1)
CTROND                              AY = GY*GAO(J) + G0*GAO1(J,2)
CTROND                              AZ = GZ*GAO(J) + G0*GAO1(J,3)
CTROND                              AR = AX*RX + AY*RY + AZ*RZ
CTROND                              AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) -AR*BR
CTROND                              EXCMAT(J,I) = EXCMAT(J,I) + FAC0*A0 
CTROND     &                                    + FACR*AR + FZ0*AB
CTROND                           END DO
CTROND                        END DO
CTROND                     END IF
CTROND                  END DO
CTROND               END IF
CTROND         END IF
      END IF
      ENDDO
      RETURN
      END

      SUBROUTINE FORTWRT(str, slength)
#include "priunit.h"
      INTEGER unit, slength
      CHARACTER str(slength)
      WRITE(LUPRI, *) str
      END
