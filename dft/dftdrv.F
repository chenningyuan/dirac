C /* Deck dftexc */
      SUBROUTINE DFTEXC(DMAT,DTRMAT,EXCMAT,NXCMAT,
     &                  IREPOP,DOERG,DOGRD,DOLND,DOATR,
     &                  WORK,LWORK,JPRINT)
C
C     T. Helgaker sep 99
C     DIRAC adaption by T.Saue
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C 
      PARAMETER (DP5 = 0.5D0)
      LOGICAL DOERG, DOGRD, DOLND, DOATR
      DIMENSION DMAT(*),DTRMAT(*),EXCMAT(*),WORK(LWORK)
C
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbham.h"
#endif
#if defined (VAR_MPI)
#include "infmpi.h"
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "infpar.h"
#endif
#include "inforb.h"
#include "nuclei.h"
#include "dftcom.h"
#include "memint.h"
C
CTROND      LUQUAD = 29
      LUQUAD = 44
C
C     Calculate grid
C
      IF (.NOT.GRDONE) THEN
         CALL DFTGRD(WORK,LWORK,1,.FALSE.)
         GRDONE = .TRUE.
      END IF
#if defined (VAR_MPI)
      IF(PARCAL) THEN
        CALL DFTDRV(DMAT,DTRMAT,EXCMAT,NXCMAT,
     &              IREPOP,DOERG,DOGRD,DOLND,DOATR,
     &              DUM,DUM,DUM,DUM,DUM,DUM,DUM,
     &              WORK(KFREE),LFREE,LUQUAD,IPRDFT)
      ELSE
#endif
C
C     Integrate
C
      CALL MEMGET('REAL',KGAO,NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGA1,3*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGA2,6*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGB1,3*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCNT,NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDST,NATOMS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF,NBAST*4,WORK,KFREE,LFREE)
      CALL DFTDRV(DMAT,DTRMAT,EXCMAT,NXCMAT,
     &            IREPOP,DOERG,DOGRD,DOLND,DOATR,
     &            WORK(KGAO),WORK(KGA1),WORK(KGA2),WORK(KGB1),
     &            WORK(KCNT),WORK(KDST),WORK(KBUF),
     &            WORK(KFREE),LFREE,LUQUAD,IPRDFT)
      CALL MEMREL('DFTDRV',WORK,KWORK,KWORK,KFREE,LFREE)
#if defined (VAR_MPI)
      ENDIF
#endif
      RETURN
      END
C /* Deck dftdrv */
      SUBROUTINE DFTDRV(DMAT,DTRMAT,EXCMAT,NXCMAT,
     &                  IREPOP,DOERG,DOGRD,DOLND,DOATR,
     &                  GAO,GAO1,GAO2,GAB1,
     &                  NCNT,DST,BUF,WORK,LWORK,LUQUAD,IPRINT)
C
C     T. Helgaker sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)
CTROND      PARAMETER (THRSH = 1.0D-8, THRSH2 = 1.0D-10) 
      PARAMETER (THRSH = 1.0D-10, THRSH2 = 1.0D-10) 
C
#if defined (VAR_MPI)
#include "infmpi.h"
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "infpar.h"
#endif
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbham.h"
#endif
#include "inforb.h"
#include "nuclei.h"
#include "dfterg.h"
#include "dftcom.h"
#include "orgcom.h"
#include "dftbuf.h"
#include "functionals.h"
C
      LOGICAL DOERG, DOGRD, DOLND, DOATR
      DIMENSION DMAT(*),DTRMAT(*),
     &          GAO(*), GAO1(*), GAO2(*),GAB1(*),
     &          WORK(LWORK),EXCMAT(NBAST,NBAST,NXCMAT), 
     &          NCNT(*), DST(*),BUF(*)
C
CTROND
      TIMSTR = SECOND()
#include "memint.h"
C
      IF (IPRINT .GE. 100) THEN
         CALL HEADER('DMAT in DFTDRV ',-1)
         CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      IF (DOERG) THEN
         NDER = 0
#if defined (PRG_DIRAC)
         EXCTRO = DP5*DDOT(NBAST*NBAST,DMAT,1,EXCMAT,1)
#else
         EXCTRO = DP5*DDOT(NBAST*NBAST,DMAT,1,EXCMAT,1)
#endif
      END IF 
      IF (DOGRD) THEN
         NDER = 1
         CALL DZERO(GRADFT,MXCOOR)
      END IF
      IF (DOLND) THEN
         NDER = 0
      END IF
      IF (DOATR) THEN
         NDER = 0
      END IF
C
      DOGGA = DFT_IS_GGA()
CTROND
      DOGRHO = .FALSE.
      DOLRHO = .FALSE.
C
      IF (DOGGA) THEN
         FROMVX = DFTPOT
         DOGRHO = .TRUE.
         NDER = NDER + 1 
         IF (FROMVX) THEN
            DOLRHO = .TRUE.
            NDER = NDER + 1 
         END IF
      ELSE
         FROMVX = .FALSE.
      END IF
C
C     integrate
C
      ELCTRN = D0
      EDFTY  = D0
C
      REWIND LUQUAD
      R2TOT = D0
      NPNTS = 0
  100 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         NPNTS = NPNTS + NPOINT
         CALL MEMGET('REAL',KX,NPOINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KY,NPOINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KZ,NPOINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KW,NPOINT,WORK,KFREE,LFREE)
         CALL REAQUA(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),LUQUAD,NPOINT)
#if defined (VAR_MPI) 
         IF(PARCAL) THEN
           NBATCH = NPOINT/NODES
           NREST = NPOINT
           KXB = KX
           KYB = KY
           KZB = KZ
           KWB = KW
           NTEST = 1
           DO NODE = 1,NODES
             CALL MPI_RECV(ITASK,1,MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                  MPI_COMM_WORLD,ISTAT,IERR)
             NWHO = ISTAT(MPI_SOURCE)
             IF(NODE.EQ.NODES) NBATCH = NREST
             CALL MPI_SEND(NBATCH,1,MPI_INTEGER,
     &                     NWHO,30,MPI_COMM_WORLD,IERR)
             CALL MPI_SEND(WORK(KXB),NBATCH,MPI_DOUBLE_PRECISION,
     &                 NWHO,32,MPI_COMM_WORLD,IERR)
             KXB = KXB + NBATCH
             CALL MPI_SEND(WORK(KYB),NBATCH,MPI_DOUBLE_PRECISION,
     &                 NWHO,33,MPI_COMM_WORLD,IERR)
             KYB = KYB + NBATCH
             CALL MPI_SEND(WORK(KZB),NBATCH,MPI_DOUBLE_PRECISION,
     &                 NWHO,34,MPI_COMM_WORLD,IERR)
             KZB = KZB + NBATCH
             CALL MPI_SEND(WORK(KWB),NBATCH,MPI_DOUBLE_PRECISION,
     &                 NWHO,35,MPI_COMM_WORLD,IERR)
             KWB = KWB + NBATCH
             NREST = NREST - NBATCH
           ENDDO
         ELSE
#endif
           CALL DFTDR1(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NPOINT,
     &               DMAT,DTRMAT,EXCMAT,NXCMAT,IREPOP,
     &               DOERG,DOGRD,DOLND,DOATR,DFTASC,
     &               GAO,GAO1,GAO2,GAB1,
     &               NCNT,DST,BUF,WORK(KFREE),LFREE,IPRINT)
#if defined (VAR_MPI) 
         ENDIF
#endif
         CALL MEMREL('DFTDR1',WORK,KWORK,KWORK,KFREE,LFREE)
         GO TO 100
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 100
      ELSE
         GO TO 200
      END IF
  200 CONTINUE
#if defined (VAR_MPI)
      IF(PARCAL) CALL DFTRES(EXCMAT,NXCMAT,WORK(KFREE),LFREE,IPRINT)
#endif
      IF (DOERG) THEN
        DENERGY = EDFTY
        EDFTY = EDFTY +EXCTRO-DP5*DDOT(NBAST*NBAST,DMAT,1,EXCMAT,1)
        DO I = 1, NBAST
          DO J = 1, I - 1 
            AVERAG = DP5*(EXCMAT(I,J,1) + EXCMAT(J,I,1))
            EXCMAT(I,J,1) = AVERAG 
            EXCMAT(J,I,1) = AVERAG 
          ENDDO
        ENDDO
      END IF 
CTROND      IF (DOATR) THEN
CTROND        DO I = 1, NBAST
CTROND          DO J = 1, I - 1 
CTROND            AVERAG = EXCMAT(I,J,1) + EXCMAT(J,I,1)
CTROND            EXCMAT(I,J,1) = AVERAG 
CTROND            EXCMAT(J,I,1) = AVERAG 
CTROND          ENDDO
CTROND        ENDDO
CTROND      END IF 
      IF (DOLND) THEN
        DO I = 1, NBAST
          DO J = 1, I - 1 
            AVERAG = DP5*(EXCMAT(I,J,1) - EXCMAT(J,I,1))
            EXCMAT(I,J,1) =   AVERAG 
            EXCMAT(J,I,1) = - AVERAG 
          ENDDO
        ENDDO
      END IF 
C
C     Print section
C
      IF (IPRINT .GE. 2) THEN 
         WRITE (LUPRI,'(/2X,A,I24)')
     &      ' Number of quadrature points:       ', NPNTS
         WRITE (LUPRI,'( 2X,A,F24.16)')
     &      ' Number of electrons:               ', ELCTRN
         IF (DOERG) THEN
CPAWSA            CALL DFT_ENE_REPORT()
            WRITE (LUPRI,'(2X,A,13X,F14.6)') 
     &          ' DFT exchange-correlation energy:', 
     &          DENERGY
         END IF
      END IF
      IF (IPRINT .GE. 100) THEN
         IF (DOERG) THEN
            CALL HEADER('EXCMAT in DFTDRV ',-1)
            CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (DOLND) THEN
            CALL HEADER('EXCMAT in DFTDRV (x direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,1),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL HEADER('EXCMAT in DFTDRV (y direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,2),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL HEADER('EXCMAT in DFTDRV (z direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,3),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (DOGRD) THEN
            CALL HEADER('GRADFT in DFTDRV ',-1)
            CALL OUTPUT(GRADFT,1,3,1,NATOMS,3,NATOMS,1,LUPRI)
         END IF
         IF (DOATR) THEN
            CALL HEADER('EXCMAT in DFTDRV ',-1)
            CALL OUTPUT(EXCMAT(1,1,1),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
      TIMDRV = SECOND() - TIMSTR
      WRITE (LUPRI,'(1X)')
      CALL TIMTXT('>>>> Time used in DFTDRV:',TIMDRV,LUPRI)
      WRITE (LUPRI,'(1X)')
      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck dftdr1 */
      SUBROUTINE DFTDR1(CORX,CORY,CORZ,WEIGHT,NPOINT,
     &               DMAT,DTRMAT,EXCMAT,NXCMAT,IREPOP,
     &               DOERG,DOGRD,DOLND,DOATR,DFTASC,
     &               GAO,GAO1,GAO2,GAB1,
     &               NCNT,DST,BUF,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Written by Trond Saue and T.Helgaker
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)
      PARAMETER (THRSH = 1.0D-8, THRSH2 = 1.0D-10) 
CTROND      PARAMETER (THRSH = 1.0D-12, THRSH2 = 1.0D-10) 
#include "mxcent.h"
C
#include "nuclei.h"
#include "inforb.h"
#include "orgcom.h"
#include "dftbuf.h"
#include "dfterg.h"
#include "functionals.h"
      LOGICAL DOERG,DOGRD,DOLND,DOATR, DFTASC
      DIMENSION CORX(NPOINT),CORY(NPOINT),CORZ(NPOINT),
     &            WEIGHT(NPOINT),
     &            DMAT(NBAST,NBAST),DTRMAT(*),
     &            GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6),
     &            GAB1(NBAST,3),EXCMAT(NBAST,NBAST,NXCMAT), 
     &            NCNT(NBAST), DST(NATOMS),BUF(NBAST,4),
     &            VX(DERVS1),
     &            WORK(LWORK)
      DIMENSION RHG(3)
      DO 300 IPNT = 1, NPOINT
       
        IF (IPRINT .GE. 100) THEN
          WRITE (LUPRI,'(2X,I6,4F20.12)') 
     &       IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
        END IF
C
        WGHT = WEIGHT(IPNT)
C
C       AOs
C       ===
C
        CALL GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,
     &              CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &              BUF,WORK,LWORK,NBAST,NDER,DOLND,IPRINT)
C
C       Density
C       =======
C
        CALL GETRHO(RHO,0,GAO,DMAT,BUF)
CTROND        IF (RHO.LE.THRSH) GO TO 300
        RHO13  = RHO**DP3
C
C       Gradient of density
C       ===================
C
        IF (DOGRHO) THEN
          CALL GETGRHO(RHG,GAO,GAO1,DMAT,BUF)      
          RHOGRD = DSQRT(RHG(1)**2 + RHG(2)**2 + RHG(3)**2)
        END IF
C
C       Numerical test
C       ==============
C
CPAWSA  IF (.FALSE.) THEN
CPAWSA     CALL DFTNUM(DODRC,DOVWN,DOLYP,DOBCK,RHO,RHO13,RHOGRD)
CPAWSA  END IF
C
C       Hessian of density
C       ===================
C
        IF (DOLRHO) THEN
          CALL GETLRHO(RXX,RXY,RXZ,RYY,RYZ,RZZ,GAO,GAO1,GAO2,
     &                 DMAT,BUF)       
          RHOLAP = RXX + RYY + RZZ
          RHOGHG = RHG(1)*RXX*RHG(1)+RHG(2)*RYY*RHG(2)+RHG(3)*RZZ*RHG(3)
     &           + D2*(RHG(1)*RXY*RHG(2) + RHG(1)*RXZ*RHG(3) 
     &               + RHG(2)*RYZ*RHG(3))
        END IF
C
C       Number of electrons
C       ===================
C
        ELCTRN = ELCTRN + WGHT*RHO
C
C       Energy
C       ======
C
        IF(DOERG) THEN
           EDFTY = EDFTY + DFTENERGY(RHO,RHOGRD)*WGHT
        ENDIF
C
C       Exchange-correlation potential 
C       ==============================
C
        IF (DOERG.OR.DOGRD.OR.DOLND) THEN
           CALL DFTPOT0(VX, WGHT, RHO, RHOGRD)
           IF(DOGGA) VX(FZ0) = 2D0*VX(FZ0)/RHOGRD
        ENDIF
C
C       Exchange-correlation contribution to Kohn-Sham matrix
C       =====================================================
C
        IF (DOERG) THEN
          IF (DOGGA .AND. .NOT.FROMVX) THEN
            CALL DOEXC2(VX(FR0),VX(FZ0),RHG,EXCMAT,GAO,GAO1)      
          ELSE
            CALL DOEXC1(VX(FR0),0,EXCMAT,GAO)
          END IF
        END IF
C
C       Exchange-correlation contribution to molecular gradient
C       =======================================================
C
        IF(DOGRD) THEN
CTROND           CALL DFTFRC
        ENDIF
C
C       London contributions to Kohn-Sham matrix 
C       ========================================
C
        IF (DOLND) THEN
CTROND          CALL DFTMAG         
        ENDIF
C
C       Hessian transformation
C       ======================
C
        IF (DOATR) THEN
          CALL DFTLTR(DTRMAT,EXCMAT,NXCMAT,IREPOP,GAO,GAO1,RHO,RHO13,
     &                RHOGRD,RHG,WGHT,DOGGA,BUF)
CTROND     &                  DOBCK,DOLYP,DOGGA,TRPLET,DTGAO)
        END IF
  300 CONTINUE
C
      RETURN
      END
#ifndef PRG_DIRAC
C  /* Deck dftxcm */
      SUBROUTINE DFTXCM(FMAT,CMO,WORK,LWORK)
C
C     T. Helgaker Sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0)
C
      DIMENSION CMO(*), FMAT(NORBT,NORBT), WORK(LWORK)
C
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
C
      KDAO = 1
      KFAO = KDAO + N2BASX
      KLST = KFAO + N2BASX
      LWRK = LWORK - KLST + 1
C
C     CALL FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,IPRINT,WORK,LWORK)
C
      CALL FCKDEN(.TRUE.,.FALSE.,WORK(KDAO),DUMMY,CMO,DUMMY,
     &            IPRINT,WORK(KLST),LWRK)
      CALL DZERO(WORK(KFAO),N2BASX)
      CALL DFTEXC(WORK(KDAO),1,WORK(KFAO),1,
     &            .TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,
     &            WORK(KLST),LWORK,IPRINT)
C
C     Transform Fock matrices to MO basis
C
      CALL DZERO(FMAT,NORBT*NORBT)
      CALL MPATB(CMO,NBAST,NORBT,NBAST,NORBT,
     &           WORK(KFAO),NBAST,NBAST,NBAST,NBAST,
     &           WORK(KDAO),NORBT,NBAST)
      CALL AMPAB(WORK(KDAO),NORBT,NBAST,NORBT,NBAST,
     &           CMO,NBAST,NORBT,NBAST,NORBT,
     &           FMAT,NORBT,NORBT)
C
      RETURN
      END
#endif
C  /* Deck reaqua */
      SUBROUTINE REAQUA(X,Y,Z,WT,LUQUAD,NPOINT)
#include "implicit.h" 
#include "priunit.h"
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      READ(LUQUAD) X,Y,Z,WT
      RETURN
      END
C /* Deck dftest */
      SUBROUTINE DFTEST(WORK,LWORK,NBUF,LUQUAD,IPRINT)
C
C     T. Helgaker sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C 
      DIMENSION WORK(LWORK)
C
#include "inforb.h"
#include "nuclei.h"
C
      KX    = 1
      KY    = KX   + NBUF
      KZ    = KY   + NBUF
      KW    = KZ   + NBUF
      KGAO  = KW   + NBUF
      KGA1  = KW   + NBUF
      KDEN  = KGA1 + 3*NBAST
      KCNT  = KDEN + NBAST*NBAST
      KLST  = KCNT + NBAST
      LWRK  = LWORK - KLST + 1
      CALL DFTES1(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NBUF,
     &            WORK(KDEN),WORK(KGAO),WORK(KGA1),WORK(KCNT),
     &            WORK(KLST),LWRK,IPRINT,LUQUAD)
      RETURN
      END
C /* Deck dftest */
      SUBROUTINE DFTES1(CORX,CORY,CORZ,WEIGHT,NBUF,DMAT,GAO,GAO1,NCNT,
     &                  WORK,LWORK,IPRINT,LUQUAD)
C
C     T. Helgaker sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
C
#include "inforb.h"
C
      DIMENSION CORX(NBUF), CORY(NBUF), CORZ(NBUF), WEIGHT(NBUF),
     &          DMAT(NBAST,NBAST), GAO(NBAST), GAO1(NBAST,3),
     &          NCNT(NBAST), WORK(LWORK)
C
      LOGICAL NODV, NODC, DOLND
C
      DOLND = .FALSE.
C
C     Calculate density matrix
C     
      NODC = .FALSE.
      NODV = .FALSE. 
      CALL GETDEN(DMAT,WORK,LWORK,NODC,NODV,IPRINT)
      ELCTRN = D0
C
      NDER = 0
C
C     Integrate
C
      REWIND LUQUAD
  100 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         CALL REAQUA(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
         DO 300 IPNT = 1, NPOINT

C HJAaJ TODO June 2001: wrong argument list to GETAOS
            call quit('Correct argument list to GETAOS in DFTES1')
chj         CALL GETAOS(GAO,GAO1,NCNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),
chj  &                  WORK,LWORK,NBAST,NDER,DOLND,IPRINT)
chj from diraos.F:
chj   SUBROUTINE GETAOS(GAO,GAO1,GAO2,GAB1,NCNT,CORPX,CORPY,CORPZ,
chj  &                  WORK,LWORK,NBAST,NDER,DOLND,IPRINT)
C HJAaJ TODO end

            DO 400 I = 1, NBAST
            DO 400 J = 1, NBAST
               ELCTRN = ELCTRN + WEIGHT(IPNT)*DMAT(I,J)*GAO(I)*GAO(J)
  400       CONTINUE 
  300    CONTINUE
         GO TO 100
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 100
      ELSE
         GO TO 200
      END IF
  200 CONTINUE
C
      WRITE (LUPRI,'(2X,A,F20.12)') ' Number of electrons ',ELCTRN
C
      RETURN 
      END
#ifndef PRG_DIRAC
C  /* Deck dftsol */
      SUBROUTINE DFTSOL(WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      DIMENSION WORK(*)
#include "infvar.h"
#include "inflin.h"
C
      KVEC = 1
      KFCK = KVEC + NWOPPT
      KLST = KFCK + NNORBT - 1
      IF (KLST .GT. LWORK) CALL STOPIT('DFTSOL',' ',KLST,LWORK)
      CALL DFTSO1(WORK(KVEC),WORK(KFCK))
      RETURN
      END
C  /* Deck dftso1 */
      SUBROUTINE DFTSO1(VEC,FCK)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
      PARAMETER (FAC = -1.0D0/4.0D0)
      DIMENSION VEC(*), FCK(*)
#include "nuctap.h"
#include "infvar.h"
#include "inflin.h"
#include "inforb.h"
      IA(I) = JWOP(1,I)
      IB(I) = JWOP(2,I)
      IAA(I) = IA(I)*(IA(I)+1)/2
      IBB(I) = IB(I)*(IB(I)+1)/2
C
      REWIND (LUMC)
      CALL MOLLAB('SIR IPH ',LUMC,LUPRI)
      READ (LUMC) 
      READ (LUMC)
      READ (LUMC) 
      READ (LUMC)
      READ (LUMC)
      READ (LUMC)
      CALL READSQ(LUMC,IRAT*NNORBT,FCK)
      DO 100 I = 1,NWOPPT
         VEC(I) = FAC*VEC(I)/(FCK(IAA(I)) - FCK(IBB(I)))
  100 CONTINUE
      IF (.FALSE.) THEN
         WRITE(LUPRI,'(//,A)')' Orbital part of solution'
         PRFAC = 0.1D0
         CALL PRKAP(NWOPPT,VEC,PRFAC,LUPRI)
      END IF
      RETURN
      END
#endif /* ifndef PRG_DIRAC */
#ifdef ENABLE_DFT_NUMERICAL_TEST
C  /* Deck dftnum */
      SUBROUTINE DFTNUM(DODRC,DOVWN,DOLYP,DOBCK,RHO,RHO13,RHOGRD)
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, 
     &           DP3 = D1/D3, ADD = 1.0D-06)
      LOGICAL DODRC,DOVWN,DOLYP,DOBCK
C
      ADD2 = D2*ADD
C
      R0 = RHO
      RP = RHO + ADD
      RM = RHO - ADD
      Z0 = RHOGRD
      ZP = RHOGRD + ADD 
      ZM = RHOGRD - ADD
      R013 = R0**DP3
      RP13 = RP**DP3
      RM13 = RM**DP3
C
      IF (DODRC) THEN
         CALL VDRC(G10,R013)
         CALL V1DRC(F10,F20,R0,R013)
C
         CALL EDRC(EP,RP,RP13)
         CALL VDRC(V1P,RP13)
         CALL EDRC(EM,RM,RM13)
         CALL VDRC(V1M,RM13)
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' DR010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' DRC10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' DRC20 test',F20,H20,F20-H20
      END IF
      IF (DOVWN) THEN
         CALL VVWN(G10,R0,R013)
         CALL V1VWN(F10,F20,R0,R013)
C
         CALL EVWN(EP,RP,RP13)
         CALL VVWN(V1P,RP,RP13)
         CALL EVWN(EM,RM,RM13)
         CALL VVWN(V1M,RM,RM13)
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' VW010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' VWN10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' VWN20 test',F20,H20,F20-H20
      END IF
C
      IF (DOBCK) THEN
         CALL GBCK(G10,G01,R0,R013,Z0)
         G01 = DP5*G01*Z0 
         CALL V1BCK(F10,F01,F20,F11,F02,R0,Z0)
C
         CALL EBCK(EP,RP,RP13,Z0)
         CALL GBCK(V1P,V2P,RP,RP13,Z0)
         V2P = DP5*Z0*V2P
         CALL EBCK(EM,RM,RM13,Z0)
         CALL GBCK(V1M,V2M,RM,RM13,Z0)
         V2M = DP5*Z0*V2M
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
         H11 = (V2P - V2M)/ADD2
C
         CALL EBCK(EP,R0,R013,ZP)
         CALL GBCK(V1P,V2P,R0,R013,ZP)
         V2P = DP5*ZP*V2P
         CALL EBCK(EM,R0,R013,ZM)
         CALL GBCK(V1M,V2M,R0,R013,ZM)
         V2M = DP5*ZM*V2M
         H01 = (EP - EM)/ADD2
         H02 = (V2P - V2M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' BC010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' BC001 test',G01,H01,G01-H01
         WRITE(LUPRI,'(A,3F15.10)')' BCK10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' BCK01 test',F01,H01,F01-H01
         WRITE(LUPRI,'(A,3F15.10)')' BCK20 test',F20,H20,F20-H20
         WRITE(LUPRI,'(A,3F15.10)')' BCK11 test',F11,H11,F11-H11
         WRITE(LUPRI,'(A,3F15.10)')' BCK02 test',F02,H02,F02-H02
      END IF
      IF (DOLYP) THEN
         CALL GLYP(G10,G01,R0,R013,Z0)
         G01 = DP5*G01*Z0 
         CALL V1LYP(F10,F01,F20,F11,F02,R0,Z0)
C
         CALL ELYP(EP,RP,RP13,Z0)
         CALL GLYP(V1P,V2P,RP,RP13,Z0)
         V2P = DP5*Z0*V2P
         CALL ELYP(EM,RM,RM13,Z0)
         CALL GLYP(V1M,V2M,RM,RM13,Z0)
         V2M = DP5*Z0*V2M
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
         H11 = (V2P - V2M)/ADD2
C
         CALL ELYP(EP,R0,R013,ZP)
         CALL GLYP(V1P,V2P,R0,R013,ZP)
         V2P = DP5*ZP*V2P
         CALL ELYP(EM,R0,R013,ZM)
         CALL GLYP(V1M,V2M,R0,R013,ZM)
         V2M = DP5*ZM*V2M
         H01 = (EP - EM)/ADD2
         H02 = (V2P - V2M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' LY010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' LY001 test',G01,H01,G01-H01
         WRITE(LUPRI,'(A,3F15.10)')' LYP10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' LYP01 test',F01,H01,F01-H01
         WRITE(LUPRI,'(A,3F15.10)')' LYP20 test',F20,H20,F20-H20
         WRITE(LUPRI,'(A,3F15.10)')' LYP11 test',F11,H11,F11-H11
         WRITE(LUPRI,'(A,3F15.10)')' LYP02 test',F02,H02,F02-H02
      END IF
      RETURN
      END
#endif /* ENABLE_DFT_NUMERICAL_TEST */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doexc1 */
      SUBROUTINE DOEXC1(VXC,IREPOP,EXCMAT,GAO)      
C***********************************************************************
C
C     Construct Kohn-Sham matrix
C
C     Written by T.Saue and T.Helgaker
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
      DIMENSION GAO(NBAST),EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),EXCMAT(NTBAS(0),NTBAS(0))
#endif
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
            IF(NB1.EQ.0) GOTO 30
             NB2 = NBBAS(IRP2,IC)
            IF(NB2.EQ.0) GOTO 30
            JB1 = IBBAS(IRP1,IC)
            JB2 = IBBAS(IRP2,IC)
            DO J = 1, NBBAS(IRP2,IC)
              GVXC = VXC*GAO(JB2+J)
              CALL DAXPY(NBBAS(IRP1,IC),GVXC,
     &                   GAO(JB1+1),1,EXCMAT(JB1+1,JB2+J),1)            
            ENDDO
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doexc2 */
      SUBROUTINE DOEXC2(VXC,VXB,RHG,EXCMAT,GAO,GAO1)      
C***********************************************************************
C
C     Get density in point (PX,PY,PZ).
C     Gets values of AOs in this point from GAO
C
C     Written by T.Saue and T.Helgaker
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
      DIMENSION GAO(NBAST),EXCMAT(NBAST,NBAST),
     &          GAO1(NBAST,3)
#else
      DIMENSION GAO(NTBAS(0)),EXCMAT(NTBAS(0),NTBAS(0)),
     &          GAO1(NTBAS(0),3)
#endif
      DIMENSION RHG(3)
C      
      MC = 2
      IF(LEVYLE) MC = 1      
      NBRP = 4/NZ
      DO IFRP = 1,NFSYM
        DO IC = 1,MC
          IP = MOD(IFRP+IC,NFSYM) + 1
          DO JSYM = 1,NBRP
            IREP = JFSYM(JSYM,IP) - 1
            IF(NBBAS(IREP,IC).GT.0) THEN
              JJ = IBBAS(IREP,IC)
              DO J = 1, NBBAS(IREP,IC)
                FC = VXC*GAO(JJ+J) + VXB*
     &               (RHG(1)*GAO1(JJ+J,1)+RHG(2)*GAO1(JJ+J,2)
     &               +RHG(3)*GAO1(JJ+J,3))
                CALL DAXPY(NBBAS(IREP,IC),FC,
     &                     GAO(JJ+1),1,EXCMAT(JJ+1,JJ+J),1)            
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck dftres */
      SUBROUTINE DFTRES(EXCMAT,NFMAT,WORK,LWORK,IPRINT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Send receive DFTresults
C     Written by T.Saue and T.helgaker
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
C
#include "mxorb.h"
#include "mxcent.h"
#if defined (VAR_MPI)
#include "infmpi.h"
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#else
#include "fpvm3h.h"
#endif
#include "infpar.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dftbuf.h"
#include "dftcom.h"
#include "dfterg.h"
#include "dgroup.h"
      DIMENSION EXCMAT(*),WORK(*)
C
#include "memint.h"
      LFMAT = N2BBASX*NFMAT*NZ
      IF (MYTID.EQ.0) THEN
        CALL MEMGET('REAL',KFBUF,LFMAT,WORK,KFREE,LFREE)
C
C     **************************
C     ***** Master receive *****
C     **************************
C
        NTEST = -1
        DO I = 1,NODES
          CALL MPI_RECV(ITASK,1,MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                  MPI_COMM_WORLD,ISTAT,IERR)
          NWHO = ISTAT(MPI_SOURCE)
          CALL MPI_SEND(NTEST,1,MPI_INTEGER,NWHO,30,
     &                    MPI_COMM_WORLD,IERR)
          CALL MPI_RECV(TEMP,1,MPI_DOUBLE_PRECISION,NWHO,40,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
          ELCTRN = ELCTRN + TEMP
          CALL MPI_RECV(TEMP,1,MPI_DOUBLE_PRECISION,NWHO,41,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
          ERGDRC = ERGDRC + TEMP
          CALL MPI_RECV(TEMP,1,MPI_DOUBLE_PRECISION,NWHO,42,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
          ERGVWN = ERGVWN + TEMP
          CALL MPI_RECV(TEMP,1,MPI_DOUBLE_PRECISION,NWHO,43,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
          ERGBCK = ERGBCK + TEMP
          CALL MPI_RECV(TEMP,1,MPI_DOUBLE_PRECISION,NWHO,44,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
          ERGLYP = ERGLYP + TEMP
          CALL MPI_RECV(WORK(KFBUF),LFMAT,MPI_DOUBLE_PRECISION,NWHO,45,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
          CALL DAXPY(LFMAT,D1,WORK(KFBUF),1,EXCMAT,1)
        ENDDO
        CALL MEMREL('DFTRES',WORK,KWORK,KWORK,KFREE,LFREE)
      ELSE 
C
C     ***********************
C     ***** Slave send  *****
C     ***********************
C
         CALL MPI_SEND(ELCTRN,1,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(ERGDRC,1,MPI_DOUBLE_PRECISION,MPARID,41,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(ERGVWN,1,MPI_DOUBLE_PRECISION,MPARID,42,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(ERGBCK,1,MPI_DOUBLE_PRECISION,MPARID,43,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(ERGLYP,1,MPI_DOUBLE_PRECISION,MPARID,44,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(EXCMAT,LFMAT,MPI_DOUBLE_PRECISION,MPARID,45,
     &                 MPI_COMM_WORLD,IERR)
      ENDIF
      RETURN
      END
C /* Deck dftpar */
      SUBROUTINE DFTPAR(DMAT,DTRMAT,EXCMAT,NXCMAT,
     &                  DOERG,DOGRD,DOLND,DOATR,
     &                  WORK,LWORK,IPRINT)
C
C     T.Saue and T. Helgaker 
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)
      PARAMETER (THRSH = 1.0D-8, THRSH2 = 1.0D-10) 
C
#if defined (VAR_MPI)
#include "infmpi.h"
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "infpar.h"
#endif
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbham.h"
#endif
#include "inforb.h"
#include "nuclei.h"
#include "dgroup.h"
#include "dfterg.h"
#include "dftcom.h"
#include "orgcom.h"
#include "dftbuf.h"
#include "functionals.h"
C
      LOGICAL DOERG, DOGRD, DOLND, DOATR
      DIMENSION DMAT(*), DTRMAT(*),
     &          WORK(LWORK),EXCMAT(*)

C
CTROND
      DIMENSION RHOB(0:2)
      TIMSTR = SECOND()
#include "memint.h"
C
        LFMAT = NBAST*NBAST*NZ*NXCMAT
        CALL DZERO(EXCMAT,LFMAT)
C
CTROND: This is temporary fix to have DFTenergy work in parallell...
      IREPOP = 0
CTROND
      CALL MEMGET('REAL',KGAO,NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGA1,3*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGA2,6*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGB1,3*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCNT,NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDST,NATOMS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF,NBAST*6,WORK,KFREE,LFREE)
      IF (DOERG) THEN
         NDER = 0
      END IF 
      IF (DOGRD) THEN
         NDER = 1
         CALL DZERO(GRADFT,MXCOOR)
      END IF
      IF (DOLND) THEN
         NDER = 0
      END IF
      IF (DOATR) THEN
         NDER = 0
      END IF
C
      DOGGA = DFT_IS_GGA()
CTROND
      DOGRHO = .FALSE.
      DOLRHO = .FALSE.
C
      IF (DOGGA) THEN
         FROMVX = DFTPOT
         DOGRHO = .TRUE.
         NDER = NDER + 1 
         IF (FROMVX) THEN
            DOLRHO = .TRUE.
            NDER = NDER + 1 
         END IF
      ELSE
         FROMVX = .FALSE.
      END IF
C
C     integrate
C
      ELCTRN = D0
C
      KFRSAV = KFREE
  100 CONTINUE
      CALL MPI_SEND(MYTID,1,MPI_INTEGER,MPARID,20,MPI_COMM_WORLD,IERR)
C
      CALL MPI_RECV(NPOINT,1,MPI_INTEGER,MPARID,30,
     &               MPI_COMM_WORLD,ISTAT,IERR)
      IF (NPOINT .EQ. -1) GOTO 300
      CALL MEMGET('REAL',KX,NPOINT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KY,NPOINT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KZ,NPOINT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KW,NPOINT,WORK,KFREE,LFREE)
      CALL MPI_RECV(WORK(KX),NPOINT,MPI_DOUBLE_PRECISION,
     &               MPARID,32,MPI_COMM_WORLD,ISTAT,IERR)
      CALL MPI_RECV(WORK(KY),NPOINT,MPI_DOUBLE_PRECISION,
     &                 MPARID,33,MPI_COMM_WORLD,ISTAT,IERR)
      CALL MPI_RECV(WORK(KZ),NPOINT,MPI_DOUBLE_PRECISION,
     &                 MPARID,34,MPI_COMM_WORLD,ISTAT,IERR)
      CALL MPI_RECV(WORK(KW),NPOINT,MPI_DOUBLE_PRECISION,
     &                 MPARID,35,MPI_COMM_WORLD,ISTAT,IERR)
      CALL DFTDR1(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NPOINT,
     &            DMAT,DTRMAT,EXCMAT,NXCMAT,IREPOP,
     &            DOERG,DOGRD,DOLND,DOATR,DFTASC,
     &            WORK(KGAO),WORK(KGA1),WORK(KGA2),WORK(KGB1),
     &            WORK(KCNT),WORK(KDST),WORK(KBUF),
     &            WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DFPAR1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      GOTO 100
C
 300  CONTINUE
      CALL MEMREL('DFPAR1',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dftprd */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DFTPRD
C
C     Communicate DIRAC common blocks for parallell DFT-calculation
C     Communicate arguments received from main program
C
C     Common blocks initialized : dcborb, dcbbas, dgroup, 
C                                 symmet
C
C     Based on traprd, Written by T. Saue 2000
C
#include "implicit.h"
      PARAMETER(D0=0.0D0)
C
#if defined (VAR_MPI)
#include "infmpi.h"
      INCLUDE 'mpif.h'
#else
#include "fpvm3h.h"
#endif
#include "priunit.h"
#include "mxorb.h"
#include "aovec.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "symmet.h"
#include "dgroup.h"
#include "infpar.h"
#include "parint.h"
#include "blocks.h"
C
#if defined (VAR_MPI)
C
C     ********************************
C     ****** Send COMMON blocks ******
C     ********************************
C
      NUMELM = NDCBORB
      CALL MPI_BCAST(NORB,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NDCBBAS
      CALL MPI_BCAST(NBBAS,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NDGROUP
      CALL MPI_BCAST(NZ,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NSYMMTI
      CALL MPI_BCAST(MAXREP,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NSYMMTR
      CALL MPI_BCAST(FMULT,NUMELM,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
#endif
C
      RETURN
      END
#endif
