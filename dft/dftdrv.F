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
C /* Deck dftexc */
      SUBROUTINE DFTEXC(DMAT,DTRMAT,EXCMAT,NXCMAT,
     &                  ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
     &                  WORK,LWORK,IPRINT)
C
C     T. Helgaker sep 99
C     DIRAC adaption by T.Saue
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
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "infpar.h"
#endif
#include "nuclei.h"
#include "dftcom.h"
#include "dftac.h"
#include "dfterg.h"
#include "dgroup.h"
#include "functionals.h"
      LOGICAL DOERG, DOGRD, DOLND, DOATR
      LOGICAL FILEISTHERE, FILEISOPENED
      DIMENSION DMAT(N2BBASX,NZ,*),EXCMAT(N2BBASX,NZ,*),
     &          ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &          DTRMAT(N2BBASX,NZ,*),WORK(LWORK)
C
#include "memint.h"
C
      LUQUAD = 44
C     Initial expectation value of two-electron matrix
      IF (DOERG) THEN
         ERRKS = DDOT(N2BBASX,DMAT,1,EXCMAT,1)
      END IF 
C     Check for Hermitian density matrices;
C     The present implementation uses only total density and no
C     spin polarization nor current density so DFT contributes
C     only for Hermitian densities
      NHM = 0
      DO IM = 1,NXCMAT
        IF(IHRMOP(IM).EQ.1) NHM = NHM + 1
      ENDDO
      IF(NHM.EQ.0) GOTO 10
C
      IF(NASHT.GT.0) CALL GTDOAV(DMAT,1)
C
C     Scale density matrices
C
      CALL DSCAL(N2BBASX,D2,DMAT,1)
      IF (DOATR) THEN
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
C     get orbital energies and density matrices for SAOP
C     ==================================================
C
      IF(DOSAOP.AND.SAOPON) THEN
C
C
C       allocate memory
C       ===============
C
        CALL MEMGET('REAL',KSAOPCMO          ,N2BBASXQ,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KSAOPDMATOUTER    ,N2BBASX,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KSAOPDMATORBOUTER ,N2ORBX,
     &              WORK,KFREE,LFREE)
        IF(DORESP) THEN
          CALL MEMGET('REAL',KSAOPDMATRESP   ,N2BBASX,
     &                WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KSAOPDMATORBRESP,N2ORBX,
     &                WORK,KFREE,LFREE)
        ENDIF
        CALL MEMGET('REAL',KEIG              ,NORBT,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KSAOPENERGIES     ,NORBT,
     &              WORK,KFREE,LFREE)
C
C
C       read coefficients and give info
C       ===============================
C
        INQUIRE(FILE = 'DFCOEF',EXIST = FILEISTHERE)
        IF(.NOT.FILEISTHERE) THEN
          WRITE(LUPRI,'(/A)') 'SAOP: Coefficient file DFCOEF not found!'
          CALL QUIT('SAOP: Coefficient file DFCOEF not found!')
        ENDIF
C
        INQUIRE(FILE = 'DFCOEF',OPENED = FILEISOPENED)
        IF(FILEISOPENED) THEN
          INQUIRE(FILE = 'DFCOEF',NUMBER = LUSAOP)
          REWIND LUSAOP
        ELSE
          LUSAOP = 42
          CALL OPNFIL(LUSAOP,'DFCOEF','UNKNOWN','SAOP')
        ENDIF
C
        CALL DZERO(WORK(KSAOPCMO),N2BBASXQ)
        CALL DZERO(WORK(KEIG),NORBT)
        CALL REACMO(LUSAOP,'DFCOEF',WORK(KSAOPCMO),WORK(KEIG),TEMP,
     &              TOTERG,6,WORK,KFREE,LFREE)
        IF(.NOT.FILEISOPENED) THEN
          CLOSE(LUSAOP,STATUS = 'KEEP')
        ENDIF
C
C
C       calculate energy differences
C       ============================
C
        K = 0
        CALL DZERO(WORK(KSAOPENERGIES),NORBT)
        DO I = 1,NFSYM
          IF(NISH(I).GT.0) THEN
            DO J = 1,NISH(I)
              K = K+1
              WORK(KSAOPENERGIES+K-1) =
     &                       EIGMAX-WORK(KEIG+J-1+(I-1)*NORB(1)+NPSH(I))
            ENDDO
          ENDIF
        ENDDO
C
C
C       get orbital density matrices
C       ============================
C
        CALL DZERO(WORK(KSAOPDMATOUTER),N2BBASX)
        CALL DZERO(WORK(KSAOPDMATORBOUTER),N2ORBX)
        IF(DORESP) THEN
          CALL DZERO(WORK(KSAOPDMATRESP),N2BBASX)
          CALL DZERO(WORK(KSAOPDMATORBRESP),N2ORBX)
        ENDIF
        K = 0
        DO I = 1,NFSYM
          IF(NISH(I).GT.0) THEN
            DO J = 1,NISH(I)
              K = K+1
              JOFF = (NORB(1)*(I-1)+NPSH(I))*NORBT
              JOFF = JOFF+NORB(1)*(I-1)+NPSH(I)
              JOFF = JOFF+(J-1)*NORBT+J-1
              IF(DORESP) THEN
                WORK(KSAOPDMATORBRESP+JOFF)
     &                                   = SQRT(WORK(KSAOPENERGIES+K-1))
              ENDIF
              WORK(KSAOPDMATORBOUTER+JOFF)
     &                          = EXP(-D2*(WORK(KSAOPENERGIES+K-1)**D2))
            ENDDO
          ENDIF
        ENDDO
C
C
C       to avoid numerical problems for degenerate orbitals
C       ===================================================
C
        DO I = 1,NFSYM
          IF(NISH(I).GT.0) THEN
            DO JHELP = 1,(NISH(I)-1)
              J = NISH(I)-JHELP
              KCHECK = NISH(1)*(I-1)+J
C
              JOFF = (NORB(1)*(I-1)+NPSH(I))*NORBT
              JOFF = JOFF+NORB(1)*(I-1)+NPSH(I)
              JOFF = JOFF+(J-1)*NORBT+J-1
C
              JOFFNEXT = JOFF+NORBT+1
C
              SOMEDIFF = ABS(WORK(KSAOPENERGIES+KCHECK)
     &                                    -WORK(KSAOPENERGIES+KCHECK-1))
              IF(SOMEDIFF.LT.1.0D-5) THEN
                IF(DORESP) THEN
                  WORK(KSAOPDMATORBRESP+JOFF)
     &                                 = WORK(KSAOPDMATORBRESP+JOFFNEXT)
                ENDIF
                WORK(KSAOPDMATORBOUTER+JOFF)
     &                                = WORK(KSAOPDMATORBOUTER+JOFFNEXT)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
C
C
C       transform density matrices
C       ==========================
C
        DO I = 1,NFSYM
          IF(NISH(I).GT.0) THEN
            IF(DORESP) THEN
              CALL QTRANS('MOAO','S',D0,
     &                    NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &                    WORK(KSAOPDMATRESP+I2BASX(I,I)),
     &                    NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &                    WORK(KSAOPDMATORBRESP+I2ORBX(I,I)),
     &                    NORBT,NORBT,1,IPQTOQ(1,0),
     &                    WORK(KSAOPCMO+ICMOQ(I)),
     &                    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                    WORK(KSAOPCMO+ICMOQ(I)),
     &                    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                    WORK(KFREE),LFREE,IPRINT)
            ENDIF
            CALL QTRANS('MOAO','S',D0,
     &                  NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &                  WORK(KSAOPDMATOUTER+I2BASX(I,I)),
     &                  NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &                  WORK(KSAOPDMATORBOUTER+I2ORBX(I,I)),
     &                  NORBT,NORBT,1,IPQTOQ(1,0),
     &                  WORK(KSAOPCMO+ICMOQ(I)),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                  WORK(KSAOPCMO+ICMOQ(I)),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                  WORK(KFREE),LFREE,IPRINT)
          ENDIF
        ENDDO
        IF(DORESP) THEN
          CALL DSCAL(N2BBASX,D2,WORK(KSAOPDMATRESP),1)
        ENDIF
        CALL DSCAL(N2BBASX,D2,WORK(KSAOPDMATOUTER),1)
C
C
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
        IF(DOSAOP.AND.SAOPON) THEN
          IF(DORESP) THEN
            CALL DFT_PARDRV(DMAT,WORK(KSAOPDMATRESP),
     &                      WORK(KSAOPDMATOUTER),
     &                      DTRMAT,EXCMAT,NXCMAT,
     &                      ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
     &                      WORK(KFREE),LFREE,LUQUAD,IPRDFT)
          ELSE
            CALL DFT_PARDRV(DMAT,DUMMY,
     &                      WORK(KSAOPDMATOUTER),
     &                      DTRMAT,EXCMAT,NXCMAT,
     &                      ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
     &                      WORK(KFREE),LFREE,LUQUAD,IPRDFT)
          ENDIF
C
          CALL MEMREL('DFTDRV',WORK,KSAOPCMO,KSAOPCMO,KFREE,LFREE)
        ELSE
          CALL DFT_PARDRV(DMAT,DUMMY,DUMMY,
     &                    DTRMAT,EXCMAT,NXCMAT,
     &                    ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
     &                    WORK(KFREE),LFREE,LUQUAD,IPRDFT)
        ENDIF
C
      ELSE
#endif
C
C     Integrate
C
      CALL MEMGET('REAL',KGAO,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGA1,3*NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGA2,6*NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGB1,3*NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCNT,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDST,NATOMS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF,NTBAS(0)*4,WORK,KFREE,LFREE)
C
      IF(DOSAOP.AND.SAOPON) THEN
        IF(DORESP) THEN
          CALL DFTDRV(DMAT,WORK(KSAOPDMATRESP),WORK(KSAOPDMATOUTER),
     &                DTRMAT,EXCMAT,NXCMAT,
     &                ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
     &                WORK(KGAO),WORK(KGA1),WORK(KGA2),WORK(KGB1),
     &                WORK(KCNT),WORK(KDST),WORK(KBUF),
     &                WORK(KFREE),LFREE,LUQUAD,IPRDFT)
        ELSE
          CALL DFTDRV(DMAT,DUMMY,WORK(KSAOPDMATOUTER),
     &                DTRMAT,EXCMAT,NXCMAT,
     &                ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
     &                WORK(KGAO),WORK(KGA1),WORK(KGA2),WORK(KGB1),
     &                WORK(KCNT),WORK(KDST),WORK(KBUF),
     &                WORK(KFREE),LFREE,LUQUAD,IPRDFT)
        ENDIF
        CALL MEMREL('DFTDRV',WORK,KSAOPCMO,KSAOPCMO,KFREE,LFREE)
      ELSE
        CALL DFTDRV(DMAT,DUMMY,DUMMY,
     &              DTRMAT,EXCMAT,NXCMAT,
     &              ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
     &              WORK(KGAO),WORK(KGA1),WORK(KGA2),WORK(KGB1),
     &              WORK(KCNT),WORK(KDST),WORK(KBUF),
     &              WORK(KFREE),LFREE,LUQUAD,IPRDFT)
        CALL MEMREL('DFTDRV',WORK,KGAO,KGAO,KFREE,LFREE)
      ENDIF
#if defined (VAR_MPI)
      ENDIF
#endif
      IF(DFTISGGA()) CALL SYMSYM(EXCMAT,NXCMAT,ISYMOP,IHRMOP)
      IF (DOATR) THEN
C       Insert half-phases
        IF(NZ.LT.4) THEN
           DO IM = 1,NXCMAT
           IF(IHRMOP(IM).EQ.1) THEN
            IREPOP = ISYMOP(IM)-1
            DO IZ = 1,NZ
              IQ = IPQTOQ(IZ,IREPOP)
              CALL Q2BPHASE('D',IQ,1,DTRMAT(1,IZ,IM))
              CALL Q2BPHASE('F',IQ,1,EXCMAT(1,IZ,IM))
            ENDDO
          ENDIF
          ENDDO
        ENDIF
      ENDIF
C
C     Symmetrize matrices
C     ..include DOGGA check !
C
C     Rescale density matrices
C
      CALL DSCAL(N2BBASX,DP5,DMAT,1)
      IF(NASHT.GT.0) CALL GTDOAV(DMAT,-1)
      IF (DOATR) THEN
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
      SUBROUTINE DFTDRV(DMAT,SAOPDMATRESP,SAOPDMATOUTER,
     &                  DTRMAT,EXCMAT,NXCMAT,
     &                  ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
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
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "infpar.h"
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
      LOGICAL DOERG, DOGRD, DOLND, DOATR
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION EXCMAT(NBAST,NBAST,NZ,NXCMAT)
#else
      DIMENSION EXCMAT(NTBAS(0),NTBAS(0),NZ,NXCMAT)
#endif
      DIMENSION DMAT(*),SAOPDMATRESP(*),SAOPDMATOUTER(*),
     &          DTRMAT(*),NCNT(*),DST(*),BUF(*),
     &          GAO(*), GAO1(*), GAO2(*),GAB1(*),
     &          ISYMOP(NXCMAT),IHRMOP(NXCMAT),WORK(LWORK)
C
CTROND
      TIMSTR = SECOND()
#include "memint.h"
C
      IF (IPRINT .GE. 100) THEN
         CALL HEADER('DMAT in DFTDRV ',-1)
         CALL OUTPUT(DMAT,1,NTBAS(0),1,NTBAS(0),
     &               NTBAS(0),NTBAS(0),1,LUPRI)
      END IF
C
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
      DOGGA = DFTISGGA()
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
     &                 DMAT,SAOPDMATRESP,SAOPDMATOUTER,
     &                 DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,
     &                 DOERG,DOGRD,DOLND,DOATR,
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
      IF(PARCAL) CALL DFTRES(EXCMAT,NXCMAT,DOERG,DOGRD,DOLND,DOATR,
     &                       WORK,KFREE,LFREE,IPRINT)
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
      IF (IPRINT .GE. 2) THEN 
         WRITE (LUPRI,'(/2X,A,I24)')
     &      ' Number of quadrature points:       ', NPNTS
         WRITE (LUPRI,'( 2X,A,F24.16)')
     &      ' Number of electrons:               ', ELCTRN
         IF (DOERG) THEN
CTROND            CALL DFT_ENE_REPORT()
            WRITE (LUPRI,'(2X,A,13X,F24.16)') 
     &          ' DFT exchange-correlation energy:', 
     &          EDFTY
         END IF
      END IF
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
         IF (DOATR) THEN
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
     &                  DMAT,SAOPDMATRESP,SAOPDMATOUTER,
     &                  DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,
     &                  DOERG,DOGRD,DOLND,DOATR,
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
#include "dfterg.h"
#include "functionals.h"
#include "dgroup.h"
      LOGICAL DOERG,DOGRD,DOLND,DOATR
      DIMENSION CORX(NPOINT),CORY(NPOINT),CORZ(NPOINT),
     &            WEIGHT(NPOINT),
     &            DMAT(NBAST,NBAST),SAOPDMATRESP(*),SAOPDMATOUTER(*),
     &            DTRMAT(*),GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6),
     &            GAB1(NBAST,3),EXCMAT(NBAST,NBAST,NZ,NXCMAT), 
     &            NCNT(NBAST), DST(NATOMS),BUF(NBAST,4),
     &            VX(DERVS1),ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &            WORK(LWORK)
      DIMENSION RHG(3),VT(3),RHOOUTERGRD(3)
C
C
C     radovan:
C   { debug asymptotic correction
C     IF(DODUMP) THEN
C       NPOINT = 200
C       DO I = 1,NPOINT
C         WEIGHT(I) = D1
C         CORX(I) = D0
C         CORY(I) = D0
C         CORZ(I) = D0+(I-1)*(10.0D0/(NPOINT-1))
C       ENDDO
C     ENDIF
C     debug asymptotic correction }
C
C
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
           IF(DOGGA) THEN
             VX(FZ0) = VX(FZ0)/RHOGRD
           ENDIF
C
C
C         asymptotic correction
C         =====================
C
C         IF(DODUMP) THEN
C           WRITE(0,*) CORZ(IPNT)
C         ENDIF
C
          IF(DOGRAC.AND.GRACON) THEN
C           R = SQRT(CORX(IPNT)**2+CORY(IPNT)**2+CORZ(IPNT)**2)
            CALL DFTAC(VX,WGHT,RHO,RHO13,DUMMY,DUMMY,RHOGRD,WORK,LWORK)
          ENDIF
C
          IF(DOSAOP.AND.SAOPON) THEN
            IF(DORESP) THEN
              CALL GETRHO(RHORESP,0,GAO,SAOPDMATRESP,BUF)
            ENDIF
C
            CALL GETRHO(RHOOUTER,0,GAO,SAOPDMATOUTER,BUF)
C
            VXZSAVE = VX(FZ0)
C
            IF(DORESP) THEN
              CALL DFTAC(VX,WGHT,RHO,DUMMY,RHORESP,RHOOUTER,RHOGRD,
     &                   WORK,LWORK)
            ELSE
              CALL DFTAC(VX,WGHT,RHO,DUMMY,DUMMY,RHOOUTER,RHOGRD,
     &                   WORK,LWORK)
            ENDIF
C
            CALL GETGRHO(RHOOUTERGRD,0,GAO,GAO1,SAOPDMATOUTER,BUF)
C
            GGAOUTER = RHG(1)*(RHOOUTER*RHG(1)-RHO*RHOOUTERGRD(1))
     &               + RHG(2)*(RHOOUTER*RHG(2)-RHO*RHOOUTERGRD(2))
     &               + RHG(3)*(RHOOUTER*RHG(3)-RHO*RHOOUTERGRD(3))
            GGAOUTER = GGAOUTER/(RHO*RHO)
C
C           VX(FR0) = VX(FR0)+VXZSAVE*GGAOUTER
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
          IF (DOGGA .AND. .NOT.FROMVX) THEN
            DO I = 1,3
              VT(I) = VX(FZ0)*RHG(I)
            ENDDO
            CALL DOEXC2(VX(FR0),VT,IREPOP,EXCMAT,GAO,GAO1)      
          ELSE
            CALL DOEXC1(VX(FR0),IREPOP,EXCMAT,GAO)
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
          CALL DFTLTR(DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,GAO,GAO1,RHO,
     &                RHOGRD,RHG,WGHT,DOGGA,BUF)
        END IF
  300 CONTINUE
C
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
C/* Deck Symsym */
      SUBROUTINE SYMSYM(AMAT,NMAT,ISYMOP,IHRMOP)
C***********************************************************************
C
C     Symmetrize only specific symmetry blocks of a matrix
C
C     Written by T. Saue Mar 10 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP5 = 0.5D0,D2=2.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AMAT(NBAST,NBAST,NZ,NMAT)
#else
      DIMENSION AMAT(NTBAS(0),NTBAS(0),NZ,NMAT)
#endif

      DIMENSION ISYMOP(NMAT),IHRMOP(NMAT)
#include "ibtfun.h"
C
      DO IM = 1,NMAT
      IF(IHRMOP(IM).EQ.1) THEN
        IREPOP = ISYMOP(IM)-1
        DO 10 IRP1 = 0,NBSYM-1
          IRP2 = IBTXOR(IRP1,IREPOP)
          IF(IRP2.LT.IRP1) THEN
          DO 20 IC = 1,MC
            NB1 = NBBAS(IRP1,IC)
            IF(NB1.EQ.0) GOTO 20
            NB2 = NBBAS(IRP2,IC)
            IF(NB2.EQ.0) GOTO 20
            JB1 = IBBAS(IRP1,IC)
            JB2 = IBBAS(IRP2,IC)
            DO I = 1,NB1
              II = JB1+I
              DO J = 1,NB2
                JJ = JB2+J
                AVERAG = DP5*(AMAT(II,JJ,1,IM) + AMAT(JJ,II,1,IM))
                AMAT(II,JJ,1,IM) = AVERAG 
                AMAT(JJ,II,1,IM) = AVERAG 
              ENDDO
            ENDDO
   20     CONTINUE
          ELSEIF(IRP2.EQ.IRP1) THEN
          DO 30 IC = 1,MC
            NB = NBBAS(IRP1,IC)
            IF(NB.EQ.0) GOTO 30
            JB = IBBAS(IRP1,IC)
            DO I = 1,NB
              II = JB+I
              DO J = 1,I-1
                JJ = JB+J
                AVERAG = DP5*(AMAT(II,JJ,1,IM) + AMAT(JJ,II,1,IM))
                AMAT(II,JJ,1,IM) = AVERAG 
                AMAT(JJ,II,1,IM) = AVERAG 
              ENDDO
            ENDDO
   30     CONTINUE
          ENDIF
   10   CONTINUE
      ENDIF
      ENDDO
C
      RETURN
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck dftres */
      SUBROUTINE DFTRES(EXCMAT,NXCMAT,DOERG,DOGRD,DOLND,DOATR,
     &                  WORK,KFREE,LFREE,IPRINT)
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
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "infpar.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dftbuf.h"
#include "dftcom.h"
#include "dfterg.h"
#include "dgroup.h"
      LOGICAL DOERG, DOGRD, DOLND, DOATR
      DIMENSION EXCMAT(*),WORK(*)
C
      KWORK = KFREE
      LFMAT = N2BBASX*NXCMAT*NZ
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
      SUBROUTINE DFTPAR(DMAT,SAOPDMATRESP,SAOPDMATOUTER,
     &                  DTRMAT,EXCMAT,NXCMAT,
     &                  ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
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
      LOGICAL DOERG, DOGRD, DOLND, DOATR
      DIMENSION DMAT(*),SAOPDMATRESP(*),SAOPDMATOUTER(*),
     &          DTRMAT(*),EXCMAT(*),ISYMOP(NXCMAT),
     &          IHRMOP(NXCMAT),WORK(*)
C
      KWORK = KFREE
      TIMSTR = SECOND()
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
      IF (DOATR) THEN
         NDER = 0
      END IF
C
      DOGGA = DFTISGGA()
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
     &            DMAT,SAOPDMATRESP,SAOPDMATOUTER,
     &            DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,
     &            DOERG,DOGRD,DOLND,DOATR,
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
      CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_HAM,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C.....dcborb.h
      IF(MYTID.EQ.0) THEN
        NUMELM = ICOMMSIZE(I1_DCBORB,I2_DCBORB)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_DCBORB,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C.....dcbbas.h
      IF(MYTID.EQ.0) THEN
        NUMELM = ICOMMSIZE(I1_DCBBAS,I2_DCBBAS)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_DCBBAS,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C.....dgroup.h
      IF(MYTID.EQ.0) THEN
        NUMELM = ICOMMSIZE(I1_DGROUP,I2_DGROUP)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_DGROUP,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C.....symmet.h
      IF(MYTID.EQ.0) THEN
        NUMELM = ICOMMSIZE(I1_SYMMTI,I2_SYMMTI)
      ENDIF
      CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_SYMMTI,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NSYMMTR
      CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FMULT,NUMELM,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pardrv */
      SUBROUTINE DFT_PARDRV(DMAT,SAOPDMATRESP,SAOPDMATOUTER,
     &                      DTRMAT,EXCMAT,NXCMAT,
     &                      ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
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
#include "dgroup.h"
C
      LOGICAL DOERG, DOGRD, DOLND, DOATR, LSEND(6), LACBCAST(8)
      DIMENSION EXCMAT(*),DMAT(N2BBASX),
     &          SAOPDMATRESP(*),SAOPDMATOUTER(*),
     &          DTRMAT(N2BBASX,NZ,NXCMAT),
     &          ISYMOP(NXCMAT),IHRMOP(NXCMAT),WORK(LWORK)
C
#include "memint.h"
      CALL QENTER('DFT_PARDRV')
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
      CALL MPI_BCAST(NXCMAT,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
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
C
C     transfer booleans for asymptotic correction
      LACBCAST(1) = DOGRAC
      LACBCAST(2) = GRACON
      LACBCAST(3) = DOSAOP
      LACBCAST(4) = SAOPON
      LACBCAST(5) = DORESP
      LACBCAST(6) = DOLB94
      LACBCAST(7) = DOLBALPHA
      LACBCAST(8) = DODUMP
      CALL MPI_BCAST(LACBCAST,8,
     &               MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
C
      IF(GRACON) THEN
        CALL MPI_BCAST(GRACALPHA,1,MPI_DOUBLE_PRECISION,
     &                             MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(GRACBETA,1,MPI_DOUBLE_PRECISION,
     &                            MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(ACIP,1,MPI_DOUBLE_PRECISION,
     &                        MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(EIGMAX,1,MPI_DOUBLE_PRECISION,
     &                          MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(ACTHRESHOLD,1,MPI_DOUBLE_PRECISION,
     &                               MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C
      IF(GRACON.OR.SAOPON) THEN
        CALL MPI_BCAST(ICOUNTAC,1,
     &                 MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C
      IF(DORESP) THEN
        CALL MPI_BCAST(SAOPDMATRESP,N2BBASX,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C
      IF(DOSAOP) THEN
        CALL MPI_BCAST(SAOPDMATOUTER,N2BBASX,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C
C
C     Transfer functionals
      CALL DFTSYNCFUNC(1)
C     Transfer control parameters
      CALL MPI_BCAST(DOERG,1,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOGRD,1,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOLND,1,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOATR,1,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(ISYMOP,NXCMAT,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IHRMOP,NXCMAT,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(DMAT,N2BBASX,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      IF(DOATR) THEN
        LMAT = N2BBASXQ*NXCMAT
        CALL MPI_BCAST(DTRMAT,LMAT,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      ENDIF
C
C-----------------------
C     Start calculation.
C-----------------------
C
      CALL DFTDRV(DMAT,SAOPDMATRESP,SAOPDMATOUTER,
     &            DTRMAT,EXCMAT,NXCMAT,
     &            ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
     &            DUM,DUM,DUM,DUM,IDUM,DUM,DUM,
     &            WORK(KFREE),LFREE,LUQUAD,IPRINT)
C
      CALL MEMREL('DFT_PARDRV',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('DFT_PARDRV')
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
#include "infpar.h"
C
      DIMENSION WORK(LWORK)
      LOGICAL DOERG,DOGRD,DOLND,DOATR,LSEND(6),LACBCAST(8)
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
      CALL MPI_BCAST(NXCMAT,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
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
C
C     receive parameters for asymptotic correction
      CALL MPI_BCAST(LACBCAST,8,
     &               MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
      DOGRAC     = LACBCAST(1)
      GRACON     = LACBCAST(2)
      DOSAOP     = LACBCAST(3)
      SAOPON     = LACBCAST(4)
      DORESP     = LACBCAST(5)
      DOLB94     = LACBCAST(6)
      DOLBALPHA  = LACBCAST(7)
      DODUMP     = LACBCAST(8)
C
      IF(GRACON) THEN
        CALL MPI_BCAST(GRACALPHA,1,MPI_DOUBLE_PRECISION,
     &                             MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(GRACBETA,1,MPI_DOUBLE_PRECISION,
     &                            MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(ACIP,1,MPI_DOUBLE_PRECISION,
     &                        MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(EIGMAX,1,MPI_DOUBLE_PRECISION,
     &                          MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(ACTHRESHOLD,1,MPI_DOUBLE_PRECISION,
     &                               MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C
      IF(GRACON.OR.SAOPON) THEN
        CALL MPI_BCAST(ICOUNTAC,1,
     &                 MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C
      IF(DORESP) THEN
        CALL MEMGET('REAL',KSAOPDMATRESP,N2BBASX,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KSAOPDMATRESP),N2BBASX,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C
      IF(DOSAOP) THEN
        CALL MEMGET('REAL',KSAOPDMATOUTER,N2BBASX,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KSAOPDMATOUTER),N2BBASX,
     &                 MPI_DOUBLE_PRECISION,MPARID,MPI_COMM_WORLD,IERR)
      ENDIF
C
C
C     Receive functionals
      CALL DFTSYNCFUNC(0)
C     Transfer control parameters
      CALL MPI_BCAST(DOERG,1,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOGRD,1,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOLND,1,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOATR,1,MPI_LOGICAL,MPARID,MPI_COMM_WORLD,IERR)
      CALL MEMGET('INTE',KSYM,NXCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHRM,NXCMAT,WORK,KFREE,LFREE)
      IF(DOSAOP) THEN
        IF(DORESP) THEN
          CALL DFTNO1(NXCMAT,WORK(KSYM),WORK(KHRM),
     &                DOERG,DOGRD,DOLND,DOATR,
     &                WORK(KSAOPDMATRESP),WORK(KSAOPDMATOUTER),
     &                WORK,KFREE,LFREE)
        ELSE
          CALL DFTNO1(NXCMAT,WORK(KSYM),WORK(KHRM),
     &                DOERG,DOGRD,DOLND,DOATR,
     &                DUMMY,WORK(KSAOPDMATOUTER),
     &                WORK,KFREE,LFREE)
        ENDIF
      ELSE
        CALL DFTNO1(NXCMAT,WORK(KSYM),WORK(KHRM),
     &              DOERG,DOGRD,DOLND,DOATR,
     &              DUMMY,DUMMY,
     &              WORK,KFREE,LFREE)
      ENDIF
C
      CALL MEMREL('DFTNOD',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('DFTNOD')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck dftno1 */
      SUBROUTINE DFTNO1(NXCMAT,ISYMOP,IHRMOP,
     &                  DOERG,DOGRD,DOLND,DOATR,
     &                  SAOPDMATRESP,SAOPDMATOUTER,WORK,KFREE,LFREE)
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
      LOGICAL DOERG, DOGRD, DOLND, DOATR
      DIMENSION ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &          SAOPDMATRESP(*),SAOPDMATOUTER(*),WORK(*)
C
      CALL MPI_BCAST(ISYMOP,NXCMAT,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IHRMOP,NXCMAT,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C
      LMAT = N2BBASXQ*NXCMAT
      CALL MEMGET('REAL',KDEN,N2BBASX,WORK,KFREE,LFREE)
      CALL MPI_BCAST(WORK(KDEN),N2BBASX,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      IF(DOATR) THEN
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
     &            SAOPDMATRESP,SAOPDMATOUTER,
     &            WORK(KDTR),WORK(KEXC),NXCMAT,
     &            ISYMOP,IHRMOP,DOERG,DOGRD,DOLND,DOATR,
     &            WORK,KFREE,LFREE,IPRINT)
      CALL DFTRES(WORK(KEXC),NXCMAT,DOERG,DOGRD,DOLND,DOATR,
     &            WORK,KFREE,LFREE,IPRINT)
      RETURN
      END
#endif
