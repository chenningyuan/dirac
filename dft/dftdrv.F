C 
!
!      Copyright (c) 2008 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
C 
C /* Deck dftexc */
      SUBROUTINE DFTEXC(DMAT,DTRMAT,EXCMAT,NXCMAT,
     &                  ISYMOP,IHRMOP,IOPTS,
     &                  WORK,LWORK,IPRINT)
C
C     T. Helgaker sep 99
C     DIRAC adaption by T.Saue
C
C     00001: DOERG, calculate DFT energy
C     00010: DOGRD, calculate molecular gradient
C     00100: DOLR,  linear response contribution
C     01000: DOQR,  quadratic response contribution
C     10000: DOLND, include London orbitals
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
C 
      PARAMETER (DP5 = 0.5D0,D2=2.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "infpar.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "nuclei.h"
#include "dftcom.h"
#include "dftac.h"
#include "cdft.h"
#include "dfterg.h"
#include "dgroup.h"
#include "functionals.h"
#include "densfit.h"
      LOGICAL DOERG, 
     &        DOGRD, 
     &        DOLND, 
     &        DOLR,
     &        LBIT,
     &        DOQR
      DIMENSION DMAT(N2BBASX,NZ,*),EXCMAT(N2BBASX,NZ,*),
     &          ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &          DTRMAT(N2BBASX,NZ,*),WORK(LWORK)
C
#include "memint.h"
C
      DOERG = LBIT(IOPTS,1)
      DOGRD = LBIT(IOPTS,2)
      DOLR  = LBIT(IOPTS,3)
      DOQR  = LBIT(IOPTS,4)
      DOLND = LBIT(IOPTS,5)
      LUQUAD = 44
C
C
C     Initial expectation value of two-electron matrix
      IF (DOERG) THEN
         ERRKS = DDOT(N2BBASX,DMAT,1,EXCMAT,1)
      END IF 


!     check hermiticity of density matrices

      NR_HERMITIAN     = 0
      NR_ANTIHERMITIAN = 0
      NR_GENERAL       = 0

      DO I = 1,NXCMAT
        IF(IHRMOP(I) .EQ. -1) NR_ANTIHERMITIAN = NR_ANTIHERMITIAN + 1
        IF(IHRMOP(I) .EQ.  0) NR_GENERAL       = NR_GENERAL       + 1
        IF(IHRMOP(I) .EQ.  1) NR_HERMITIAN     = NR_HERMITIAN     + 1
      ENDDO

      IF(DOLR) THEN
        IF(NR_GENERAL .GT. 0)
     &           CALL QUIT('DFT LR: general density matrix or matrices')
        IF(NR_HERMITIAN .EQ. 0) GOTO 10
      ENDIF

      IF(NASHT.GT.0) CALL GTDOAV(DMAT,1)
C
C     Scale density matrices
C
      CALL DSCAL(N2BBASX,D2,DMAT,1)
      IF (DOLR.OR.DOQR) THEN
        CALL DSCAL(N2BBASXQ*NXCMAT,D2,DTRMAT,1)
C       Insert half-phases:
C       The symmetry scheme tries to pack the real symemtric and three 
C       imaginary anti-symmetric matrices of the original quaternion 
C       matrix. The DFT contribution goes into the real part. 
C       For a totally symmetric density (IREPOP=0) the DFT contribution 
C       will go into the first packed matrix (IZ=1) having the correct
C       symmetry about the diagonal, but this is not generally true.
C       The insertion of half-phases assures the correcy symmetry about
C       the diagonal for the appropriate symmetry blocks...
        IF(NZ.LT.4) THEN
          DO IM = 1,NXCMAT
            IREPOP = ISYMOP(IM)-1
            DO IZ = 1,NZ
              IQ = IPQTOQ(IZ,IREPOP)
              CALL Q2BPHASE('D',IQ,1,DTRMAT(1,IZ,IM))
              CALL Q2BPHASE('F',IQ,1,EXCMAT(1,IZ,IM))
            ENDDO
          ENDDO
        ENDIF
      END IF
C
C
C     get density matrices for SAOP
C     =============================
C
      KDMAT_OUTER_AO = KFREE
      KDMAT_RESP_AO  = KFREE
      IF(LBIT(LBIT_DFTAC,3)) THEN !SAOP
        CALL SET_LBIT(LBIT_DFTAC,4) !SAOP is active

        CALL MEMGET('REAL',KDMAT_OUTER_AO,N2BBASX,
     &              WORK,KFREE,LFREE)
        IF(LBIT(LBIT_DFTAC,5)) THEN !SAOP with response part
          CALL MEMGET('REAL',KDMAT_RESP_AO,N2BBASX,
     &                WORK,KFREE,LFREE)
        ENDIF
C
        CALL DFTAC_SAOP_GET_DMATS(WORK(KDMAT_OUTER_AO),
     &                            WORK(KDMAT_RESP_AO),
     &                            WORK(KFREE),LFREE)
      ENDIF
C
C
C     Calculate grid
C
      IF (.NOT.GRDONE) THEN
         CALL DFTGRD(DMAT,WORK(KFREE),LFREE,1)
         GRDONE = .TRUE.
      END IF
#if defined (VAR_MPI)
      IF(PARCAL) THEN
C       Get hold of slaves
        CALL DIRAC_PARCTL( 4 )
C
        CALL DFT_PARDRV(DMAT,WORK(KDMAT_RESP_AO),
     &                  WORK(KDMAT_OUTER_AO),
     &                  DTRMAT,EXCMAT,NXCMAT,
     &                  ISYMOP,IHRMOP,IOPTS,
     &                  WORK(KFREE),LFREE,LUQUAD,IPRDFT)
C
      ELSE
#endif
C
C     Integrate
        NDIMENSION = NTBAS(0)
C
        CALL MEMGET('REAL',KGAO,NDIMENSION,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KGA1,3*NDIMENSION,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KGA2,6*NDIMENSION,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KGB1,3*NDIMENSION,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KCNT,NDIMENSION,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KDST,NATOMS,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBUF,NDIMENSION*4,WORK,KFREE,LFREE)
C        
        CALL DFTDRV(DMAT,WORK(KDMAT_RESP_AO),WORK(KDMAT_OUTER_AO),
     &              DTRMAT,EXCMAT,NXCMAT,
     &              ISYMOP,IHRMOP,IOPTS,
     &              WORK(KGAO),WORK(KGA1),WORK(KGA2),WORK(KGB1),
     &              WORK(KCNT),WORK(KDST),WORK(KBUF),
     &              WORK(KFREE),LFREE,LUQUAD,IPRDFT)
C             SAOP is active
      IF(LBIT(LBIT_DFTAC,4)) THEN
        CALL MEMREL('DFTDRV',WORK,KDMAT_OUTER_AO,KDMAT_OUTER_AO,
     &              KFREE,LFREE)
      ELSE
        CALL MEMREL('DFTDRV-buf',WORK,KBUF,KBUF,KFREE,LFREE)
        CALL MEMREL('DFTDRV-dst',WORK,KDST,KDST,KFREE,LFREE)
        CALL MEMREL('DFTDRV-cnt',WORK,KCNT,KCNT,KFREE,LFREE)
        CALL MEMREL('DFTDRV-gb1',WORK,KGB1,KGB1,KFREE,LFREE)
        CALL MEMREL('DFTDRV-ga2',WORK,KGA2,KGA2,KFREE,LFREE)
        CALL MEMREL('DFTDRV-ga1',WORK,KGA1,KGA1,KFREE,LFREE)
        CALL MEMREL('DFTDRV-gao',WORK,KGAO,KGAO,KFREE,LFREE)
C       CALL MEMREL('DFTDRV-all',WORK,KGAO,KGAO,
C    &              KFREE,LFREE)
      ENDIF
#if defined (VAR_MPI)
      ENDIF
#endif


!     in case of GGA symmetrize matrices
!     see DFT RPA notes for more info

      IF(DFTISGGA() .AND. .NOT. DF_FITVXC) THEN
        CALL SYM_HERMITIAN(EXCMAT,NXCMAT,ISYMOP,IHRMOP)
      ENDIF

      IF (DOLR.OR.DOQR) THEN
C       Insert half-phases
        IF(NZ.LT.4) THEN
           DO IM = 1,NXCMAT
            IREPOP = ISYMOP(IM)-1
            DO IZ = 1,NZ
              IQ = IPQTOQ(IZ,IREPOP)
              CALL Q2BPHASE('D',IQ,1,DTRMAT(1,IZ,IM))
              CALL Q2BPHASE('F',IQ,1,EXCMAT(1,IZ,IM))
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C     Rescale density matrices
C
      CALL DSCAL(N2BBASX,DP5,DMAT,1)
      IF(NASHT.GT.0) CALL GTDOAV(DMAT,-1)
      IF (DOLR.OR.DOQR) THEN
        CALL DSCAL(N2BBASXQ*NXCMAT,DP5,DTRMAT,1)
      END IF
C
C     Calculate final expectation value of two-electron matrix;
C     the Kohn-Sham contribution can then be subtracted 
      IF (DOERG) THEN
        ERRKS = DDOT(N2BBASX,DMAT,1,EXCMAT,1) - ERRKS
      END IF 
C
 10   RETURN
      END
C /* Deck dftdrv */
      SUBROUTINE DFTDRV(DMAT,DMAT_RESP_AO,DMAT_OUTER_AO,
     &                  DTRMAT,EXCMAT,NXCMAT,
     &                  ISYMOP,IHRMOP,IOPTS,
     &                  GAO,GAO1,GAO2,GAB1,NCNT,DST,BUF,
     &                  WORK,LWORK,LUQUAD,IPRINT)
C
C     T. Helgaker sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#include "infpar.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#endif
C#include "inforb.h"
#include "nuclei.h"
#include "dfterg.h"
#include "dftcom.h"
#include "orgcom.h"
#include "dftbuf.h"
#include "functionals.h"
#include "dgroup.h"
C
      LOGICAL DOERG, DOGRD, DOLND, DOLR, DOQR, LBIT
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION EXCMAT(NBAST,NBAST,NZ,NXCMAT)
#else
      DIMENSION EXCMAT(NTBAS(0),NTBAS(0),NZ,NXCMAT)
#endif
      DIMENSION DMAT(*),DMAT_RESP_AO(*),DMAT_OUTER_AO(*),
     &          DTRMAT(*),NCNT(*),DST(*),BUF(*),
     &          GAO(*), GAO1(*), GAO2(*),GAB1(*),
     &          ISYMOP(NXCMAT),IHRMOP(NXCMAT),WORK(LWORK)
C
CTROND
      TIMSTR = SECOND()
#include "memint.h"
C
      DOERG = LBIT(IOPTS,1)
      DOGRD = LBIT(IOPTS,2)
      DOLR  = LBIT(IOPTS,3)
      DOQR  = LBIT(IOPTS,4)
      DOLND = LBIT(IOPTS,5)
C
      IF (IPRINT .GE. 100) THEN
         CALL HEADER('DMAT in DFTDRV ',-1)
         CALL OUTPUT(DMAT,1,NTBAS(0),1,NTBAS(0),
     &               NTBAS(0),NTBAS(0),1,LUPRI)
      END IF
C
      DOGGA = DFTISGGA()
      IF (DOGRD) THEN
         NDER = 1
         CALL DZERO(GRADFT,MXCOOR)
      END IF
      IF (DOLND) THEN
         NDER = 0
      END IF
      IF (DOLR.OR.DOQR) THEN
         NDER = 0
         IF(ALDAON) DOGGA = .FALSE.
      END IF
C
      DOGRHO = .FALSE.
      DOLRHO = .FALSE.
#if defined (VAR_MPI) 
      DO NODE = 1,NODES
        CALL MPI_SEND(DOGGA,1,MPI_LOGICAL,NODE,11,
     &                MPI_COMM_WORLD,IERR)
      END DO             
#endif
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

      IF(import_grid) THEN
        READ(LUQUAD,*) NPOINT
      ELSE
        READ(LUQUAD) NPOINT
      ENDIF

      IF (NPOINT.GT.0) THEN
         NPNTS = NPNTS + NPOINT
         CALL MEMGET('REAL',KX,NPOINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KY,NPOINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KZ,NPOINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KW,NPOINT,WORK,KFREE,LFREE)

         IF(import_grid) THEN
           CALL import_dft_grid(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),
     &                   LUQUAD,NPOINT)
         ELSE
           CALL REAQUA(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),
     &                 LUQUAD,NPOINT)
         ENDIF

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
             CALL MPI_RECV(ITASK,1,my_MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                  MPI_COMM_WORLD,ISTAT,IERR)
             NWHO = ISTAT(MPI_SOURCE)
             IF(NODE.EQ.NODES) NBATCH = NREST
             CALL MPI_SEND(NBATCH,1,my_MPI_INTEGER,
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
     &                 DMAT,DMAT_RESP_AO,DMAT_OUTER_AO,
     &                 DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,IOPTS,
     &                 GAO,GAO1,GAO2,GAB1,
     &                 NCNT,DST,BUF,WORK(KFREE),LFREE,IPRINT)
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
      IF(PARCAL) CALL DFTRES(EXCMAT,NXCMAT,IOPTS,WORK,KFREE,LFREE,
     &                       IPRINT)
#endif
      IF (DOLND) THEN
        DO I = 1, NTBAS(0)
          DO J = 1, I - 1 
            AVERAG = DP5*(EXCMAT(I,J,1,1) - EXCMAT(J,I,1,1))
            EXCMAT(I,J,1,1) =   AVERAG 
            EXCMAT(J,I,1,1) = - AVERAG 
          ENDDO
        ENDDO
      END IF 
C
C     Test on the number of electrons
C
      ELCTRX = NELECT+NAELEC
      ERROR  = ELCTRN - ELCTRX
      IF (ABS(ERROR) .GT. DFTELS) THEN
         WRITE (LUPRI,'(4(/2X,A,F14.6),/2X,A)')
     &   ' Number of electrons from numerical integration:',ELCTRN,
     &   ' Number of electrons from orbital occupations:  ',ELCTRX,
     &   ' Error in the number of electrons:              ',ERROR,
     &   ' Error larger than DFTELS (set input):          ',DFTELS,
     &   ' Calculation aborted.'
CTROND         CALL QUIT
CTROND     &    ('Wrong number of electrons in DFTDRV. Calculation aborted.')
        WRITE(LUPRI,'(A)') 'WARNING !!!'
      END IF
C
C
C     Print section
C
      IF(IPRINT .GE. 2) THEN 
        WRITE(LUPRI,'(2X,A,I24)')
     &              ' Number of quadrature points:       ', NPNTS
        WRITE(LUPRI,'(2X,A,F24.16)')
     &              ' Number of electrons:               ', ELCTRN
        IF(DOERG) THEN
          WRITE(LUPRI,'(2X,A,F24.16)')
     &                ' DFT exchange-correlation energy:   ', EDFTY
        ENDIF
      ENDIF


      IF (IPRINT .GE. 100) THEN
         IF (DOERG) THEN
            CALL HEADER('EXCMAT in DFTDRV ',-1)
            CALL OUTPUT(EXCMAT,1,NTBAS(0),1,NTBAS(0),NTBAS(0),
     &                  NTBAS(0),1,LUPRI)
         END IF
         IF (DOLND) THEN
            CALL HEADER('EXCMAT in DFTDRV (x direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,1,1),
     &                  1,NTBAS(0),1,NTBAS(0),NTBAS(0),NTBAS(0),1,LUPRI)
            CALL HEADER('EXCMAT in DFTDRV (y direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,1,2),
     &                  1,NTBAS(0),1,NTBAS(0),NTBAS(0),NTBAS(0),1,LUPRI)
            CALL HEADER('EXCMAT in DFTDRV (z direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,1,3),
     &                  1,NTBAS(0),1,NTBAS(0),NTBAS(0),NTBAS(0),1,LUPRI)
         END IF
         IF (DOGRD) THEN
            CALL HEADER('GRADFT in DFTDRV ',-1)
            CALL OUTPUT(GRADFT,1,3,1,NATOMS,3,NATOMS,1,LUPRI)
         END IF
         IF (DOLR.OR.DOQR) THEN
            CALL HEADER('EXCMAT in DFTDRV ',-1)
            CALL OUTPUT(EXCMAT(1,1,1,1),
     &                  1,NTBAS(0),1,NTBAS(0),NTBAS(0),NTBAS(0),1,LUPRI)
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
     &                  DMAT,DMAT_RESP_AO,DMAT_OUTER_AO,
     &                  DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,IOPTS,
     &                  GAO,GAO1,GAO2,GAB1,
     &                  NCNT,DST,BUF,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Written by Trond Saue and T.Helgaker
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)
#include "mxcent.h"
C
#include "nuclei.h"
#include "inforb.h"
#include "orgcom.h"
#include "dftbuf.h"
#include "dftcom.h"
#include "dftac.h"
#include "cdft.h"
#include "dfterg.h"
#include "functionals.h"
#include "dgroup.h"
#include "densfit.h"
      LOGICAL DOERG,DOGRD,DOLND,DOLR,DOQR,LBIT
      DIMENSION CORX(NPOINT),CORY(NPOINT),CORZ(NPOINT),
     &            WEIGHT(NPOINT),
     &            DMAT(NBAST,NBAST),DMAT_RESP_AO(*),DMAT_OUTER_AO(*),
     &            DTRMAT(*),GAO(*), GAO1(*), GAO2(NBAST,6),
     &            GAB1(NBAST,3),EXCMAT(NBAST,NBAST,NZ,NXCMAT), 
     &            NCNT(NBAST), DST(NATOMS),BUF(NBAST,4),
     &            VX(DERVS1),ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &            WORK(LWORK)
      DIMENSION RHG(3),VT(3)
C
#include "memint.h"
C
      DOERG = LBIT(IOPTS,1)
      DOGRD = LBIT(IOPTS,2)
      DOLR  = LBIT(IOPTS,3)
      DOQR  = LBIT(IOPTS,4)
      DOLND = LBIT(IOPTS,5)
C
      
      NDIMENSION = NBAST



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
C
        CALL GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,
     &              CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &              BUF,WORK,LWORK,NDIMENSION,NDER,DOLND,IPRINT)
C
C       Density
C       =======
C
          CALL GETRHO(RHO,0,GAO,DMAT,BUF)
Caspg+rbast, 2006-09-21
        IF (RHO.LE.DFTHR0) GO TO 300
        RHO13  = RHO**DP3
C
C       Gradient of density
C       ===================
C
        IF (DOGRHO) THEN
            CALL GETGRHO(RHG,0,GAO,GAO1,DMAT,BUF)      
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
           
!          WRITE(lupri,*) 'DOGGA2', DOGGA
           IF (DOGGA) THEN
             VX(FZ0) = VX(FZ0)/RHOGRD
           ENDIF
C
C
C         asymptotic correction
C         =====================
C
C                 SAOP is active
          IF(LBIT(LBIT_DFTAC,4)) THEN
            CALL GETRHO(RHO_OUTER,0,GAO,DMAT_OUTER_AO,BUF)
C                   SAOP with response part
            IF(LBIT(LBIT_DFTAC,5)) THEN
              CALL GETRHO(RHO_RESP,0,GAO,DMAT_RESP_AO,BUF)
            ENDIF
          ENDIF
C                 GRAC is active          SAOP is active
          IF(LBIT(LBIT_DFTAC,2) .OR. LBIT(LBIT_DFTAC,4)) THEN
            CALL DFTAC(VX,WGHT,
     &                 RHO,RHO13,RHO_RESP,RHO_OUTER,RHOGRD)
          ENDIF
C
C
        ENDIF
C
C
C       Exchange-correlation contribution to Kohn-Sham matrix
C       =====================================================
C
        IF (DOERG) THEN
           IREPOP = 0
C conventional code...
!            WRITE(lupri,*) 'DOGGA4', DOGGA
             IF (DOGGA .AND. .NOT.FROMVX) THEN
                DO I = 1,3
                   VT(I) = VX(FZ0)*RHG(I)
                ENDDO
                CALL DOEXC2(VX(FR0),VT,IREPOP,EXCMAT,GAO,GAO1)      
             ELSE
                CALL DOEXC1(VX(FR0),IREPOP,EXCMAT,GAO)
             END IF
        END IF ! DOERG

!
!       Exchange-correlation contribution to molecular gradient
!       =======================================================
!
!       IF(DOGRD) THEN
!TROND           CALL DFTFRC
!       ENDIF
!
!       London contributions to Kohn-Sham matrix 
!       ========================================
!
!       IF (DOLND) THEN
!TROND          CALL DFTMAG         
!       ENDIF




!       linear response
        if(dolr) then
          call dft_lr(dtrmat,excmat,nxcmat,
     &                isymop,ihrmop,.false.,
     &                gao,gao1,
     &                rho,rhogrd,rhg,wght,dogga,buf)

        endif

!       quadratic response
        if(doqr) then
          call dft_qr(dtrmat,excmat,
     &                isymop,ihrmop,gao,gao1,
     &                rho,rhogrd,rhg,wght,dogga,buf)
        endif

  300 CONTINUE !end of loop over grid points




      RETURN
      END
C  /* Deck reaqua */
      SUBROUTINE REAQUA(X,Y,Z,WT,LUQUAD,NPOINT)
#include "implicit.h" 
#include "priunit.h"
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      READ(LUQUAD) X,Y,Z,WT
      RETURN
      END
#ifndef PRG_DIRAC
C  /* Deck dftsol */
      SUBROUTINE DFTSOL(WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
      DIMENSION WORK(*)
#include "mx2orb.h"
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
      PARAMETER (FAC = -1.0D0/4.0D0)
      DIMENSION VEC(*), FCK(*)
#include "nuctap.h"
#include "mx2orb.h"
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
      DO 10 IRP1 = 0,NBSYM-1
        IRP2 = IBTXOR(IRP1,IREPOP)
        DO 20 IC = 1,MC
          NB1 = NBBAS(IRP1,IC)
          IF(NB1.EQ.0) GOTO 20
          NB2 = NBBAS(IRP2,IC)
          IF(NB2.EQ.0) GOTO 20
          JB1 = IBBAS(IRP1,IC)
          JB2 = IBBAS(IRP2,IC)
          DO J = 1, NBBAS(IRP2,IC)
            GVXC = VXC*GAO(JB2+J)
            CALL DAXPY(NBBAS(IRP1,IC),GVXC,
     &                 GAO(JB1+1),1,EXCMAT(JB1+1,JB2+J),1)            
          ENDDO
   20   CONTINUE
   10 CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doexc2 */
      SUBROUTINE DOEXC2(SXC,VXC,IREPOP,EXCMAT,GAO,GAO1)      
C***********************************************************************
C
C     Construct the matrix
C
C     F(i,j) = GAO(I)*SXC*GAO(J)+2*GAO(I)*VXC(K)*nabla[GAO1(J)]
C
C     Upon symmetrization 
C
C     F(i,j) = (1/2)[F(i,j)+F(j,i)]
C
C     later on one obtains
C
C     E(i,j) = GAO(I)*SXC*GAO(J)+VXC(K)*nabla[GAO(I)*GAO(J)]
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
      DIMENSION VXC(3)
#include "ibtfun.h"
C.....Vector part is scaled to get symmetrization right
      DO I = 1,3
        VXC(I)=D2*VXC(I)
      ENDDO
      DO 10 IRP1 = 0,NBSYM-1
        IRP2 = IBTXOR(IRP1,IREPOP)
        DO 20 IC = 1,MC
          NB1 = NBBAS(IRP1,IC)
          IF(NB1.EQ.0) GOTO 20
          NB2 = NBBAS(IRP2,IC)
          IF(NB2.EQ.0) GOTO 20
          JB1 = IBBAS(IRP1,IC)
          JB2 = IBBAS(IRP2,IC)
          DO J = 1, NBBAS(IRP2,IC)
            FC = SXC*GAO(JB2+J)      +VXC(1)*GAO1(JB2+J,1)
     &         + VXC(2)*GAO1(JB2+J,2)+VXC(3)*GAO1(JB2+J,3)
            CALL DAXPY(NBBAS(IRP1,IC),FC,
     &                 GAO(JB1+1),1,EXCMAT(JB1+1,JB2+J),1)            
          ENDDO
   20   CONTINUE
   10 CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Getavd */
      SUBROUTINE GTDOAV(DMAT,IFAC)
C***********************************************************************
C
C     Get average density from open shells; only real part !
C     Written by T.Saue March 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION DMAT(N2BBASX,NZ,NFMAT)
      IF(AOC) THEN
        DO ISHELL = 1, NOPEN
          FAC = IFAC*DF(ISHELL)
          CALL DAXPY(N2BBASX,FAC,DMAT(1,1,ISHELL+1),1,DMAT,1)
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SYM_ANTIHERMITIAN */
      SUBROUTINE SYM_ANTIHERMITIAN(AMAT,
     &                             NMAT,
     &                             ISYMOP,
     &                             IHRMOP)
C=======================================================================
C     symmetrize only specific symmetry blocks of a matrix
C     written by Trond Saue march 2004
C     adapted to antihermitian matrices
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                           last revision: february 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (DP5 = 0.5D0,D2 = 2.0D0)

#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AMAT(NBAST,NBAST,NZ,NMAT)
#else
      DIMENSION AMAT(NTBAS(0),NTBAS(0),NZ,NMAT)
#endif
      DIMENSION ISYMOP(NMAT),
     &          IHRMOP(NMAT)

#include "ibtfun.h"

      DO IM = 1,NMAT
        IF(.NOT. IHRMOP(IM) .EQ. 1) THEN
          DO IXYZ = 1,3
            iq = jm4pos(4-ixyz)
            iq = iqtopq(iq,isymax(ixyz,2))
            IREPOP = IBTXOR(ISYMAX(IXYZ,2),ISYMOP(IM)-1)
            DO 110 IREP1 = 0,NBSYM-1
              IREP2 = IBTXOR(IREP1,IREPOP)
              IF(IREP2 .LT. IREP1) THEN
                DO 120 IC = 1,MC
                  NB1 = NBBAS(IREP1,IC)
                  IF(NB1 .EQ. 0) GOTO 120
                  NB2 = NBBAS(IREP2,IC)
                  IF(NB2 .EQ. 0) GOTO 120
                  JB1 = IBBAS(IREP1,IC)
                  JB2 = IBBAS(IREP2,IC)
                  DO I = 1,NB1
                    II = JB1 + I
                    DO J = 1,NB2
                      JJ = JB2 + J
                      AVERAGE = DP5*(AMAT(II,JJ,IQ,IM) 
     &                             + AMAT(JJ,II,IQ,IM))
                      AMAT(II,JJ,IQ,IM) = AVERAGE 
                      AMAT(JJ,II,IQ,IM) = AVERAGE 
                    ENDDO
                  ENDDO
  120           CONTINUE
              ELSEIF(IREP2 .EQ. IREP1) THEN
                DO 130 IC = 1,MC
                  NB = NBBAS(IREP1,IC)
                  IF(NB .EQ. 0) GOTO 130
                  JB = IBBAS(IREP1,IC)
                  DO I = 1,NB
                    II = JB + I
                    DO J = 1,I-1
                      JJ = JB + J
                      AVERAGE = DP5*(AMAT(II,JJ,IQ,IM) 
     &                             + AMAT(JJ,II,IQ,IM))
                      AMAT(II,JJ,IQ,IM) = AVERAGE 
                      AMAT(JJ,II,IQ,IM) = AVERAGE 
                    ENDDO
                  ENDDO
  130           CONTINUE
              ENDIF
  110       CONTINUE
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SYM_HERMITIAN */
      SUBROUTINE SYM_HERMITIAN(AMAT,
     &                         NMAT,
     &                         ISYMOP,
     &                         IHRMOP)
C=======================================================================
C     symmetrize only specific symmetry blocks of a matrix
C     written by Trond Saue march 2004
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (DP5 = 0.5D0,D2 = 2.0D0)

#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AMAT(NBAST,NBAST,NZ,NMAT)
#else
      DIMENSION AMAT(NTBAS(0),NTBAS(0),NZ,NMAT)
#endif
      DIMENSION ISYMOP(NMAT),
     &          IHRMOP(NMAT)

#include "ibtfun.h"

      DO IM = 1,NMAT
        IF(.NOT. IHRMOP(IM) .EQ. -1) THEN
          IREPOP = ISYMOP(IM) - 1
          DO 110 IREP1 = 0,NBSYM-1
            IREP2 = IBTXOR(IREP1,IREPOP)
            IF(IREP2 .LT. IREP1) THEN
              DO 120 IC = 1,MC
                NB1 = NBBAS(IREP1,IC)
                IF(NB1 .EQ. 0) GOTO 120
                NB2 = NBBAS(IREP2,IC)
                IF(NB2 .EQ. 0) GOTO 120
                JB1 = IBBAS(IREP1,IC)
                JB2 = IBBAS(IREP2,IC)
                DO I = 1,NB1
                  II = JB1 + I
                  DO J = 1,NB2
                    JJ = JB2 + J
                    AVERAGE = DP5*(AMAT(II,JJ,1,IM) 
     &                           + AMAT(JJ,II,1,IM))
                    AMAT(II,JJ,1,IM) = AVERAGE 
                    AMAT(JJ,II,1,IM) = AVERAGE 
                  ENDDO
                ENDDO
  120         CONTINUE
            ELSEIF(IREP2 .EQ. IREP1) THEN
              DO 130 IC = 1,MC
                NB = NBBAS(IREP1,IC)
                IF(NB .EQ. 0) GOTO 130
                JB = IBBAS(IREP1,IC)
                DO I = 1,NB
                  II = JB + I
                  DO J = 1,I-1
                    JJ = JB + J
                    AVERAGE = DP5*(AMAT(II,JJ,1,IM) 
     &                           + AMAT(JJ,II,1,IM))
                    AMAT(II,JJ,1,IM) = AVERAGE 
                    AMAT(JJ,II,1,IM) = AVERAGE 
                  ENDDO
                ENDDO
  130         CONTINUE
            ENDIF
  110     CONTINUE
        ENDIF
      ENDDO

      RETURN
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck dftres */
      SUBROUTINE DFTRES(EXCMAT,NXCMAT,IOPTS,WORK,KFREE,LFREE,IPRINT)
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
#include "maxorb.h"
#include "mxcent.h"
#include "infpar.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dftbuf.h"
#include "dftcom.h"
#include "dfterg.h"
#include "dgroup.h"
      LOGICAL DOERG, DOGRD, DOLND, DOLR, DOQR, LBIT
      DIMENSION EXCMAT(*),WORK(*)
C
      KWORK = KFREE
      LFMAT = N2BBASX*NXCMAT*NZ
C
      DOERG = LBIT(IOPTS,1)
      DOGRD = LBIT(IOPTS,2)
      DOLR  = LBIT(IOPTS,3)
      DOQR  = LBIT(IOPTS,4)
      DOLND = LBIT(IOPTS,5)
C
      IF (MYTID.EQ.0) THEN
        CALL MEMGET('REAL',KFBUF,LFMAT,WORK,KFREE,LFREE)
C
C     **************************
C     ***** Master receive *****
C     **************************
C
        NTEST = -1
        DO I = 1,NODES
          CALL MPI_RECV(ITASK,1,my_MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                  MPI_COMM_WORLD,ISTAT,IERR)
          NWHO = ISTAT(MPI_SOURCE)
          CALL MPI_SEND(NTEST,1,my_MPI_INTEGER,NWHO,30,
     &                    MPI_COMM_WORLD,IERR)
          CALL MPI_RECV(TEMP,1,MPI_DOUBLE_PRECISION,NWHO,40,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
          ELCTRN = ELCTRN + TEMP
          IF(DOERG) THEN
            CALL MPI_RECV(TEMP,1,MPI_DOUBLE_PRECISION,NWHO,40,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
            EDFTY = EDFTY + TEMP
          ENDIF
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
         IF(DOERG) THEN
           CALL MPI_SEND(EDFTY,1,MPI_DOUBLE_PRECISION,MPARID,40,
     &                   MPI_COMM_WORLD,IERR)
         ENDIF
         CALL MPI_SEND(EXCMAT,LFMAT,MPI_DOUBLE_PRECISION,MPARID,45,
     &                 MPI_COMM_WORLD,IERR)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck dftpar */
      SUBROUTINE DFTPAR(DMAT,DMAT_RESP_AO,DMAT_OUTER_AO,
     &                  DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,IOPTS,
     &                  WORK,KFREE,LFREE,IPRINT)
C
C     T.Saue and T. Helgaker 
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)
C
#if defined (VAR_MPI)
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
      LOGICAL DOERG, DOGRD, DOLND, DOLR, DOQR, LBIT
      DIMENSION DMAT(*),DMAT_RESP_AO(*),DMAT_OUTER_AO(*),
     &          DTRMAT(*),EXCMAT(*),ISYMOP(NXCMAT),
     &          IHRMOP(NXCMAT),WORK(*)
C
      KWORK = KFREE
      TIMSTR = SECOND()
C
      DOERG = LBIT(IOPTS,1)
      DOGRD = LBIT(IOPTS,2)
      DOLR  = LBIT(IOPTS,3)
      DOQR  = LBIT(IOPTS,4)
      DOLND = LBIT(IOPTS,5)
C
      LEXCMAT = NBAST*NBAST*NZ*NXCMAT
      CALL DZERO(EXCMAT,LEXCMAT)
C
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
      IF (DOLR.OR.DOQR) THEN
         NDER = 0
      END IF
C
      DOGGA = DFTISGGA()
      CALL MPI_RECV(DOGGA,1,MPI_LOGICAL,MPARID,11,
     &              MPI_COMM_WORLD,ISTAT,IERR)
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
      KFRSAV = KFREE
  100 CONTINUE
      CALL MPI_SEND(MYTID,1,my_MPI_INTEGER,MPARID,20,
     &              MPI_COMM_WORLD,IERR)
C
      CALL MPI_RECV(NPOINT,1,my_MPI_INTEGER,MPARID,30,
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
     &            DMAT,DMAT_RESP_AO,DMAT_OUTER_AO,
     &            DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,IOPTS,
     &            WORK(KGAO),WORK(KGA1),WORK(KGA2),WORK(KGB1),
     &            WORK(KCNT),WORK(KDST),WORK(KBUF),
     &            WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DFTPAR.loop',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      GOTO 100
C
 300  CONTINUE
      CALL MEMREL('DFTPAR.end',WORK,KWORK,KWORK,KFREE,LFREE)
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
      INCLUDE 'mpif.h'
#endif
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
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
C.....dcbham.h
      IF(MYTID.EQ.0) THEN
        NUMELM = ICOMMSIZE(I1_HAM,I2_HAM)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_HAM,NUMELM,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C.....dcborb.h
      IF(MYTID.EQ.0) THEN
        NUMELM = ICOMMSIZE(I1_DCBORB,I2_DCBORB)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_DCBORB,NUMELM,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C.....dcbbas.h
      IF(MYTID.EQ.0) THEN
        NUMELM = ICOMMSIZE(I1_DCBBAS,I2_DCBBAS)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_DCBBAS,NUMELM,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C.....dgroup.h
      IF(MYTID.EQ.0) THEN
        NUMELM = ICOMMSIZE(I1_DGROUP,I2_DGROUP)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_DGROUP,NUMELM,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C.....symmet.h
      IF(MYTID.EQ.0) THEN
        NUMELM = ICOMMSIZE(I1_SYMMTI,I2_SYMMTI)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_SYMMTI,NUMELM,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NSYMMTR
      CALL MPI_BCAST(NUMELM,1,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FMULT,NUMELM,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pardrv */
      SUBROUTINE DFT_PARDRV(DMAT,DMAT_RESP_AO,DMAT_OUTER_AO,
     &                      DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,IOPTS,
     &                      WORK,LWORK,LUQUAD,IPRINT)
C
C     A cut and paste job by T. Saue Sep 25 2003
C     Based on HER_PARDRV which was
C     written by Henrik Koch and Paal Dahle Nov-1993 -> Okt-1994
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
C
#include "infpar.h"
#include "dcbbas.h"
#include "dftcom.h"
#include "dftac.h"
#include "cdft.h"
#include "dgroup.h"
C
      LOGICAL DOERG, DOGRD, DOLND, DOLR,DOQR,LBIT, LSEND(6)
C
      DIMENSION EXCMAT(*),DMAT(N2BBASX),
     &          DMAT_RESP_AO(*),DMAT_OUTER_AO(*),
     &          DTRMAT(N2BBASX,NZ,NXCMAT),
     &          ISYMOP(NXCMAT),IHRMOP(NXCMAT),WORK(LWORK)
C
#include "memint.h"
      CALL QENTER('DFT_PARDRV')
C
      DOERG = LBIT(IOPTS,1)
      DOGRD = LBIT(IOPTS,2)
      DOLR  = LBIT(IOPTS,3)
      DOQR  = LBIT(IOPTS,4)
      DOLND = LBIT(IOPTS,5)
C
C     Ensure that NODES .LE. MAXNOD
C
      IF (NODES .GT. MAXNOD) THEN
         WRITE(LUPRI,'(/5X,A,/12X,A,/5X,A16,I8,/5X,A16,I8)')
     &    'ERROR: Number of nodes allocated has exceeded hardcoded',
     &    'limit. Reduce number of nodes or recompile program. ',
     &    'Allocated (NODES) :',NODES,'Limit (MAXNOD) :',MAXNOD
         CALL QUIT('Number of nodes requested/allocated is too large.')
      ELSE IF (IPRINT .GT. 0) THEN
         WRITE(LUPRI,'(/5X,A,I5/)')
     &    'DFT_PARDRV: Number of nodes allocated:', NODES
      ENDIF
      TIMING = .FALSE.
C
C-----------------------
C     Send data
C-----------------------
C
      CALL MPI_BCAST(NXCMAT,1,my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
C     Transfer various COMMON blocks
      CALL DFTPRD
C     Transfer dftcom
      LSEND(1)  = DFTADD
      LSEND(2)  = GRDONE
      LSEND(3)  = DFTRUN
      LSEND(4)  = DFTPOT
      LSEND(5)  = DFTORD
      LSEND(6)  = DFTASC
      CALL MPI_BCAST(LSEND,6,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
C


!     broadcast thresholds
      CALL MPI_BCAST(DFTHR0,1,
     &               MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTHRL,1,
     &               MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)



C     transfer parameters for asymptotic correction
C     =============================================
C
      CALL MPI_BCAST(LBIT_DFTAC,1,
     &               my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
C             GRAC is active
      IF(LBIT(LBIT_DFTAC,2)) THEN
        CALL MPI_BCAST(GRAC_ALPHA,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(GRAC_BETA,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(DFTAC_IP,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(HOMO_EIGENVALUE,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(DFTAC_THRESHOLD,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C             SAOP is active
      IF(LBIT(LBIT_DFTAC,4)) THEN
        CALL MPI_BCAST(DMAT_OUTER_AO,N2BBASX,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C             SAOP with response part
      IF(LBIT(LBIT_DFTAC,5)) THEN
        CALL MPI_BCAST(DMAT_RESP_AO,N2BBASX,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C
C
C     Transfer functionals
      CALL DFTSYNCFUNC(1)
C     Transfer control parameters
C     write(lupri,*) 'Before IOPTS broadcast'
      CALL MPI_BCAST(IOPTS,1,my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
C     write(lupri,*) 'After IOPTS broadcast'
C
      CALL MPI_BCAST(ISYMOP,NXCMAT,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IHRMOP,NXCMAT,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(DMAT,N2BBASX,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      IF(DOLR.OR.DOQR) THEN
        LMAT = N2BBASXQ*NXCMAT
        CALL MPI_BCAST(DTRMAT,LMAT,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      ENDIF

C     write(lupri,*) 'Broadcasted'
C
C-----------------------
C     Start calculation.
C-----------------------
C
      CALL DFTDRV(DMAT,DMAT_RESP_AO,DMAT_OUTER_AO,
     &            DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,IOPTS,
     &            DUM,DUM,DUM,DUM,IDUM,DUM,DUM,
     &            WORK(KFREE),LFREE,LUQUAD,IPRINT)
C     write(lupri,*) 'DFTDRV done'
C
      CALL MEMREL('DFT_PARDRV',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('DFT_PARDRV')
C     write(lupri,*) 'Exit DFT_PARDRV'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dftnod */
      SUBROUTINE DFTNOD(WORK,LWORK)
C***********************************************************************
C     Based on HERNOD 
C     written by T. Saue March 4
C
C    This is the node program for DFT                           *
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
C
#include "dcbbas.h"
#include "dftcom.h"
#include "dftac.h"
#include "cdft.h"
#include "infpar.h"
C
      DIMENSION WORK(LWORK)
      LOGICAL DOERG,DOGRD,DOLND,DOLR,LSEND(6),LBIT
C
      CALL QENTER('DFTNOD')
#include "memint.h"
C
      PARHER = .TRUE.
C
C     Start timing.
C
      CALL GETTIM(CPU1,WALL1)
C
C     This is a slave in a parallel run
C
      SLAVE = .TRUE.
      DEBUG  = .FALSE.
C
C     ***************************************************
C     ***** Receive initialization data from master *****
C     ***************************************************
C
      CALL MPI_BCAST(NXCMAT,1,my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
C     Various common blocks
      CALL DFTPRD
C     Receive some dftcom variables
      CALL MPI_BCAST(LSEND,6,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
      DFTADD = LSEND(1)
      GRDONE = LSEND(2)
      DFTRUN = LSEND(3)
      DFTPOT = LSEND(4)
      DFTORD = LSEND(5)
      DFTASC = LSEND(6)
C


!     broadcast thresholds
      CALL MPI_BCAST(DFTHR0,1,
     &               MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTHRL,1,
     &               MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)

C     receive parameters for asymptotic correction
C     ============================================
C
      CALL MPI_BCAST(LBIT_DFTAC,1,
     &               my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
C             GRAC is active
      IF(LBIT(LBIT_DFTAC,2)) THEN
        CALL MPI_BCAST(GRAC_ALPHA,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(GRAC_BETA,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(DFTAC_IP,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(HOMO_EIGENVALUE,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(DFTAC_THRESHOLD,1,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C             SAOP is active
      IF(LBIT(LBIT_DFTAC,4)) THEN
        CALL MEMGET('REAL',KDMAT_OUTER_AO,N2BBASX,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KDMAT_OUTER_AO),N2BBASX,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ELSE
        KSAOPDMATRESP = KFREE
      ENDIF
C             SAOP with response part
      IF(LBIT(LBIT_DFTAC,5)) THEN
        CALL MEMGET('REAL',KDMAT_RESP_AO,N2BBASX,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KDMAT_RESP_AO),N2BBASX,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ELSE
        KSAOPDMATOUTER = KFREE
      ENDIF
C
C
C     Receive functionals
      CALL DFTSYNCFUNC(0)
C     Transfer control parameters
      CALL MPI_BCAST(IOPTS,1,my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      CALL MEMGET('INTE',KSYM,NXCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHRM,NXCMAT,WORK,KFREE,LFREE)
      CALL DFTNO1(NXCMAT,WORK(KSYM),WORK(KHRM),IOPTS,
     &            WORK(KDMAT_RESP_AO),WORK(KDMAT_OUTER_AO),
     &            WORK,KFREE,LFREE)
C
      CALL MEMREL('DFTNOD',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('DFTNOD')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck dftno1 */
      SUBROUTINE DFTNO1(NXCMAT,ISYMOP,IHRMOP,IOPTS,
     &                  DMAT_RESP_AO,DMAT_OUTER_AO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Next part of DFT node program.
C     T. Saue March 4 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "infpar.h"
#include "dcbbas.h"
      LOGICAL DOERG, DOGRD, DOLND, DOLR, DOQR, LBIT
      DIMENSION ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &          DMAT_RESP_AO(*),DMAT_OUTER_AO(*),WORK(*)
C
      DOERG = LBIT(IOPTS,1)
      DOGRD = LBIT(IOPTS,2)
      DOLR  = LBIT(IOPTS,3)
      DOQR  = LBIT(IOPTS,4)
      DOLND = LBIT(IOPTS,5)
C
      CALL MPI_BCAST(ISYMOP,NXCMAT,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IHRMOP,NXCMAT,my_MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C
      LMAT = N2BBASXQ*NXCMAT
      CALL MEMGET('REAL',KDEN,N2BBASX,WORK,KFREE,LFREE)
      CALL MPI_BCAST(WORK(KDEN),N2BBASX,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      IF(DOLR.OR.DOQR) THEN
        CALL MEMGET('REAL',KDTR,LMAT,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KDTR),LMAT,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      ELSE
        KDTR = KFREE
      ENDIF
      CALL MEMGET('REAL',KEXC,LMAT,WORK,KFREE,LFREE)
C
C
      IPRINT = 0
      CALL DFTPAR(WORK(KDEN),
     &            DMAT_RESP_AO,DMAT_OUTER_AO,
     &            WORK(KDTR),WORK(KEXC),NXCMAT,ISYMOP,IHRMOP,IOPTS,
     &            WORK,KFREE,LFREE,IPRINT)
      CALL DFTRES(WORK(KEXC),NXCMAT,IOPTS,WORK,KFREE,LFREE,IPRINT)
      RETURN
      END
#endif

