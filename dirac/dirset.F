C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
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
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dirset.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamset */
      SUBROUTINE PAMSET(WORK,LWORK)
C*****************************************************************************
C
C     Subroutine for setting dimensions and generate overlap and Lowdin
C     matrices
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
C
#include <dcbham.h>
#include <dcblab.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbgen.h>
#include <cbihr1.h>
#include <infpar.h>
C
      DIMENSION WORK(LWORK)
C
C     Set DCB common blocks - first time
C     ==================================
C
C      CALL SETDC1(IPREAD)
C
C
C     Make primitive labels
C     =====================
C
      CALL GETLAB(WORK,LWORK,IPRGEN)
C
C     Generate all one-electron integrals
C     ===================================
C
      CALL ONEGEN(WORK,LWORK)
      IF(DOHRM) RETURN
C
C     Generate MO-transformation matrix and write to file
C     ===================================================
C
      CALL GMOTRA(WORK,LWORK)
C
C     Set DCB common blocks - second time
C     ===================================
C
      CALL SETDC2(IPRGEN)
C
C     Prepare for two-electron integrals
C     ==================================
C
      IF(.NOT.ONESYS) CALL TWOPRP(WORK,LWORK)
C
C     Prepare for parallel calculation (set inforb)
C     ==============================================
C
      CALL RELINF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setdc1 */
      SUBROUTINE SETDC1(IPRINT)
C*****************************************************************************
C
C     Set DIRAC common blocks (DCB)
C
C     The general variable syntax is:
C
C     prefix1+prefix2+body+suffix1+suffix2 (irp,ic)
C
C     prefix1:
C       N   - length of vector
C       N2  - length of full matrix
C       NN  - length of row-packed matrix(lower triangle)
C       I   - offset of vector
C       I2  - offset of full matrix
C       IN  - offset of row-packed matrix(lower triangle)
C       J   - pointer of vector
C       J2  - pointer of full matrix
C       JJ  - pointer of row-packed matrix(lower triangle)
C
C     prefix2:
C       B - boson irreps
C       F - fermion irreps
C
C     body:
C       BAS - AO-basis
C       ORB - MO-basis
C       ISH - inactive spinors
C
C     suffix1:
C       T - symmetry-packing
C       X - no symmetry-packing
C       Q - quaternionic symmetry-packed matrix (NZ matrices adjoined)
C
C     indices:
C       irp - irrep (boson or fermion)
C       ic  - component:
C              0 - total
C              1 - large
C              2 - small
C
C     Written by T.Saue - October 1994
C     Last revision: Nov 23 1994 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
C
#include <dcbgen.h>
#include <dcbbas.h>
#include <symmet.h>
#include <pgroup.h>
#include <dgroup.h>
#include <dcbfir.h>
      CHARACTER COMP(2)*1
C
      IF(IPRINT.GE.2) CALL TITLER('Output from SETDC1','*',103)
C
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     *********************************************
C     ***** COMMON/DCBBAS/ - Basis dimensions *****
C     *********************************************
C
C     NTBAS(ic)      - total number of basis functions for component ic
C     NBBAS(ibrp,ic) - number of SO basis functions
C                      for boson irrep ibrp and component ic
C     NFBAS(ifrp,ic) - number of SO basis functions
C                      for fermion irrep ifrp and component ic
C     MXBBAS         - max.number of SO basis functions in a boson irrep
C
C     Matrix dimensions:
C        N2BBAST,NNBBAST,N2BBASX,N2BBASXQ
C        N2BAS(ifrp),N2BASQ(ifrp),N2BAST,
C     Offsets:
C        IBBAS(ibrp,ic),
C        I2BBASX(ibrp,jbrp,ic,jc)
C        IBAS(ifrp),I2BASX(ifrp,ifrp),I2BAST(ifrp)
C
      NNBAST = 0
      DO IBREP = 0,MAXREP
        NBBAS(IBREP,0) = NAOS(IBREP+1)
        NNBBAST        = NNBBAST + NBBAS(IBREP,0)*(NBBAS(IBREP,0)+1)/2
      ENDDO
      CALL IZERO(NTBAS,3)
      CALL IZERO(NFBAS,6)
      MXBBAS   = 0
      N2BBAST  = 0
      IBASI    = 0
      NBRP     = 4/NZ
C
      IF(NFSYM.EQ.1)     THEN
        IBAS(1)    = IBASI
        DO IC = 1,2
          DO ISYM = 1,NBRP
            IBSYM           = JFSYM(ISYM,1)
            IBREP           = IBSYM - 1
            IBBAS(IBREP,IC) = IBASI
            NBBAS(IBREP,IC) = NCOS(IBSYM,IC)
            IBASI           = NCOS(IBSYM,IC) + IBASI
            NFBAS(1,IC)     = NCOS(IBSYM,IC) + NFBAS(1,IC)
            MXBBAS          = MAX(MXBBAS,NBBAS(IBREP,IC))
          ENDDO
          NTBAS(IC) = NTBAS(IC) + NFBAS(1,IC)
        ENDDO
        NFBAS(1,0) = NFBAS(1,1)  + NFBAS(1,2)
        N2BAS(1)   = NFBAS(1,0)*NFBAS(1,0)
        N2BAST     = N2BAS(1)
        N2BASQ(1)  = N2BAS(1)*NZ
        I2BAST(1)  = 0
        N2BASTQ    = N2BASQ(1)
      ELSEIF(NFSYM.EQ.2) THEN
        DO IFRP = 1,NFSYM
          IBAS(IFRP)    = IBASI
          DO IC = 1,2
C           IP is gerade(1)/ungerade(2) !
            IP    = MOD(IFRP+IC,2) + 1
            DO ISYM = 1,NBRP
              IBSYM           = JFSYM(ISYM,IP)
              IBREP           = IBSYM - 1
              IBBAS(IBREP,IC) = IBASI
              NBBAS(IBREP,IC) = NCOS(IBSYM,IC)
              IBASI           = NCOS(IBSYM,IC) + IBASI
              NFBAS(IFRP,IC)  = NCOS(IBSYM,IC) + NFBAS(IFRP,IC)
              MXBBAS          = MAX(MXBBAS,NBBAS(IBREP,IC))
            ENDDO
            NTBAS(IC)       = NTBAS(IC) + NFBAS(IFRP,IC)
          ENDDO
          NFBAS(IFRP,0) = NFBAS(IFRP,1)  + NFBAS(IFRP,2)
          N2BAS(IFRP)   = NFBAS(IFRP,0)*NFBAS(IFRP,0)
          N2BAST        = N2BAST  + N2BAS(IFRP)
          N2BASQ(IFRP)  = N2BAS(IFRP)*NZ
          I2BAST(IFRP)  = N2BASTQ
          N2BASTQ       = N2BASTQ + N2BASQ(IFRP)
        ENDDO
      ENDIF
      NTBAS(0)   = NTBAS(1) + NTBAS(2)
      NNBBASX    = NTBAS(0)*(NTBAS(0)+1)/2
      N2BBASX    = NTBAS(0)*NTBAS(0)
      N2BBASXQ   = N2BBASX*NZ
      MXFBAS = MAX(NFBAS(1,0),NFBAS(2,0))
C
      DO JC = 1,2
        DO IC = 1,2
          DO JBRP = 0,NBSYM-1
            DO IBREP = 0,NBSYM-1
              I2BBASX(IBREP,JBRP,IC,JC) =
     &               NTBAS(0)*IBBAS(JBRP,JC)+IBBAS(IBREP,IC)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C
      IOFF = 0
      DO J = 1,NFSYM
        DO I = 1,NFSYM
          I2BASX(I,J) = IOFF + IBAS(I)
        ENDDO
        IOFF = IOFF + NTBAS(0)*NFBAS(J,0)
      ENDDO
C
C     Print section
C     =============
C
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,'(2X,6(1X,A6))')
     &        'Irrep ','NTBAS ','NBASL ','NBASS ','IBASL ','IBASS '
        WRITE(LUPRI,'(3X,41A1)') ('-',I=1,41)
        WRITE(LUPRI,'(3X,I1,2X,A3,5I7)')
     &    ((I+1),REP(I),(NBBAS(I,J),J=0,2),(IBBAS(I,J),J=1,2),
     &     I=0,NBSYM-1)
        WRITE(LUPRI,'(3X,41A1)') ('-',I=1,41)
        WRITE(LUPRI,'(3X,A4,I2,3I7)')
     &       ('Frep',I,(NFBAS(I,J),J=0,2),I = 1,NFSYM)
        WRITE(LUPRI,'(3X,41A1)') ('-',I=1,41)
        WRITE(LUPRI,'(3X,A6,3I7,4X,A3,I7)') 'Total ',(NTBAS(J),J=0,2),
     &        'Max',MXBBAS
        WRITE(LUPRI,'(/3X,2(A,I6,2X))') 'NNBBAST = ',NNBBAST,
     &        'NNBBASX ',NNBBASX
C        WRITE(LUPRI,'(/A)') '* I2BBASX:'
C        WRITE(LUPRI,'(4(3X,2(A3,2X),2(A1,2X),2X,I8))')
C     &     ((((REP(I),REP(J),COMP(IC),COMP(JC),
C     &        I2BBASX(I,J,IC,JC),IC = 1,2),JC = 1,2)
C     &        ,I=0,NBSYM-1),J=0,NBSYM-1)
      ENDIF
C
C     *****************************************
C     ***** COMMON/PAMIOU/ - Unit numbers *****
C     *****************************************
C
      LUCOEF = 1
      LUOVLP = 2
      LUTMAT = 3
      LU1INT = 19
C
C     *****************************************
C     ***** COMMON/FIRSTN/ -              *****
C     *****************************************
C
C     
      FIRST1 = .TRUE.
      FIRST2 = .TRUE.
      FIRST3 = .TRUE.
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setdc2 */
      SUBROUTINE SETDC2(IPRINT)
C*****************************************************************************
C
C     Set DIRAC common blocks (DCB) - second run
C
C     The general variable syntax is:
C
C     prefix1+prefix2+body+suffix1+suffix2 (irp,ic)
C
C     prefix1:
C       N   - length of vector
C       N2  - length of full matrix
C       NN  - length of row-packed matrix(lower triangle)
C       I   - offset of vector
C       I2  - offset of full matrix
C       IN  - offset of row-packed matrix(lower triangle)
C       J   - pointer of vector
C       J2  - pointer of full matrix
C       JJ  - pointer of row-packed matrix(lower triangle)
C
C     prefix2:
C       B - boson irreps
C       F - fermion irreps
C
C     body:
C       BAS - AO-basis
C       ORB - MO-basis
C       ISH - inactive spinors
C
C     suffix1:
C       T - symmetry-packing
C       X - no symmetry-packing
C       Q - quaternionic symmetry-packed matrix (NZ matrices adjoined)
C
C     indices:
C       irp - irrep (boson or fermion)
C       ic  - component:
C              0 - total
C              1 - large
C              2 - small
C
C     Written by T.Saue - November 1994
C     Last revision: Nov 23 1994 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcbham.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <pgroup.h>
#include <dgroup.h>
      CHARACTER COMP(2)*1
C
C
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     *********************************************
C     ***** COMMON/DCBBAS/ - Basis dimensions *****
C     *********************************************
C
C     The following dimensions are set in LOWDIN:
C       NFORB(ifrp,ic) - number of MO basis functions
C                        for fermion irrep ifrp and component ic
C     Here we set:
C       NTORB(ic)      - number of MO basis functions for a component ic
C
      IF(.NOT.URKBAL) NTORB(2) = NTORB(1)
C
C     **********************************************
C     ***** COMMON/DCBORB - Orbital dimensions *****
C     **********************************************
C
C     A shell is occupied by a Kramers pair
C
C     NORB(ifrp) - total number of shells
C     NISH(ifrp) - inactive (electron) shells
C     NASH(ifrp) - active (electron) shells
C     NOCC(ifrp) - total number of partially/fully occupied shells
C     NSSH(ifrp) - secondary (electron) shells
C     NESH(ifrp) - electron shells
C     NPSH(ifrp) - positron shells
C
C     NCMO(ifrp) - dimension of MO - coefficient matrix
C     NCAO(ifrp) - dimension of occupied part of MO-coef. matrix
C
      NORBT  = 0
      N2ORBT = 0
      NFROT  = 0
      NISHT  = 0
      N2ISHT = 0
      NASHT  = 0
      N2ASHT = 0
      NOCCT  = 0
      N2OCCT = 0
      NSSHT  = 0
      N2SSHT = 0
      NPSHT  = 0
      N2PSHT = 0
      NCMOT  = 0
      NCAOT  = 0
      N2TMT  = 0
      DO 100 IFRP = 1,NFSYM
        IORB(IFRP)  = NORBT
        IASH(IFRP)  = NASHT
        NESH(IFRP)  = NFORB(IFRP,1)
        IF(URKBAL) THEN
          NZT         = 1
          NPSH(IFRP)  = NFORB(IFRP,2)
        ELSEIF(VEXTPJ.OR.FREEPJ.OR.LEVYLE) THEN
          NZT         = NZ
          NPSH(IFRP)  = 0
        ELSE
          NZT         = NZ
          NPSH(IFRP)  = MIN(NFORB(IFRP,1),NFORB(IFRP,2))
        ENDIF
        NORB(IFRP)  = NESH(IFRP) + NPSH(IFRP)
        N2ORB(IFRP) = NORB(IFRP)*NORB(IFRP)
        NORBT       = NORBT  + NORB(IFRP)
        I2ORBT(IFRP)= N2ORBT*NZ
        N2ORBT      = N2ORBT + N2ORB(IFRP)
        N2ISH(IFRP) = NISH(IFRP)*NISH(IFRP)
        NFROT       = NFROT  + NFRO(IFRP)
        NISHT       = NISHT  + NISH(IFRP)
        N2ISHT      = N2ISHT + N2ISH(IFRP)
        N2ASH(IFRP) = NASH(IFRP)*NASH(IFRP)
        NASHT       = NASHT  + NASH(IFRP)
        I2ASHT(IFRP)= N2ASHT*NZ
        N2ASHT      = N2ASHT + N2ASH(IFRP)
        N2OCC(IFRP) = NOCC(IFRP)*NOCC(IFRP)
        NOCCT       = NOCCT  + NOCC(IFRP)
        N2OCCT      = N2OCCT + N2OCC(IFRP)
        NSSH(IFRP)  = NESH(IFRP) - NOCC(IFRP)
        N2SSH(IFRP) = NSSH(IFRP)*NSSH(IFRP)
        NSSHT       = NSSHT  + NSSH(IFRP)
        N2SSHT      = N2SSHT + N2SSH(IFRP)
        N2ESH(IFRP) = NESH(IFRP)*N2ESH(IFRP)
        NESHT       = NESHT  + NESH(IFRP)
        N2ESHT      = N2ESHT + N2ESH(IFRP)
        N2PSH(IFRP) = NPSH(IFRP)*N2PSH(IFRP)
        NPSHT       = NPSHT  + NPSH(IFRP)
        N2PSHT      = N2PSHT + N2PSH(IFRP)
        ICMO(IFRP)  = NCMOT
        ICMOQ(IFRP) = NCMOT*NZ
        NCMO(IFRP)  = NFBAS(IFRP,0)*NORB(IFRP)
        NCMOQ(IFRP) = NCMO(IFRP)*NZ
        NCMOT       = NCMOT + NCMO(IFRP)
        I2TMT(IFRP) = N2TMT
        N2TMT       = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZT
  100 CONTINUE
      N2ORBTQ = N2ORBT*NZ
      NCMOTQ  = NCMOT*NZ
      IOFF = 0
      DO 200 J = 1,NFSYM
        DO 300 I = 1,NFSYM
          I2ORBX(I,J) = IOFF + IORB(I)
 300    CONTINUE
        IOFF = IOFF + NORB(J)*NORBT
 200  CONTINUE
      IOFF = 0
      DO J = 1,NFSYM
         DO I = 1,NFSYM
            I2ASHX(I,J) = IOFF + IASH(I)
         END DO
         IOFF = IOFF + NASH(J)*NASHT
      END DO
      N2ORBX  = NORBT*NORBT
      N2ORBXQ = NZ*N2ORBX
      N2ISHX  = NISHT*NISHT
      N2ASHX  = NASHT*NASHT
      N2OCCX  = NOCCT*NOCCT
      N2SSHX  = NSSHT*NSSHT
      N2ESHX  = NESHT*NESHT
      N2PSHX  = NPSHT*NPSHT
C
      IF(IPRINT.GE.1) THEN
        CALL HEADER('Output from SETDC2',-1)
        WRITE(LUPRI,'(A7,9X,A5,14X,A8)') 
     &       'Fermion','Basis','Orbitals'
        WRITE(LUPRI,'(A6,2(5X,3A5))')
     &       'ircop  ','Total',' Lbas',' Sbas',
     &                 'Total',' Eorb',' Porb'
        WRITE(LUPRI,'(2X,A3,1X,5X,3I5,5X,3I5)')
     &    (FREP(I),(NFBAS(I,J),J = 0,2),(NFORB(I,J),J = 0,2),
     &     I = 1,NFSYM)
      ENDIF
C
C     Sanety check if NOCC > NORB
C
      DO I = 1,NFSYM
         IF (NOCC(I) .GT. NORB(I)) THEN
            WRITE(LUPRI,*) 
     &         'Number of occ. orbitals (',NOCC(I),
     &         ') greater than orbitals (',NORB(I),') for symmetry ',I
            CALL QUIT('NOCC>NORB')
         END IF
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setdhf */
      SUBROUTINE SETDHF(IPRINT)
C*****************************************************************************
C
C     Set dimensions for COMMON blocks in SCF-module
C
C     Written by T.Saue November 1994
C     Last revision : Nov 23 1994 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbdhf.h>
C
C     *****************************************
C     ***** COMMON/DHFIOU/ - Unit numbers *****
C     *****************************************
C
      LUCYCL    =  9
      LUFCK1    = 10
      LUFCK2    = 11
      LUFOCK    = 12
      LUDIIS    = 13
      LU2HER    = 20
      LUDENS    = 21
      LUEVEC    = 22
      LUCMOS    = 23
      LUSMOS    = 24
      N2EVECTQ  = 0
      NSMOTQ    = 0
      DO IFRP = 1,NFSYM
         ISMOQ(IFRP) = NSMOTQ
         N2EVECTQ = N2EVECTQ + N2ORB(IFRP)*NZ
         NSMOTQ   = NSMOTQ   + NOCC(IFRP)*NORB(IFRP)*NZ
      ENDDO
C
C     Closed shell Hartree - Fock
C     ===========================
C
      IF (NOPEN.GT.0) THEN
         NFMAT = NOPEN+1
      ELSE
         IF (NASHT.GT.0) THEN
            NFMAT = 2
         ELSE
            NFMAT = 1
         END IF
      END IF
      DO I = 1,NFMAT
C        Totally symmetric operator
         ISYMOP(I) = 1
C        Fock matrix type
         IFCKOP(I) = 1
C        Hermitian operator
         IHRMOP(I) = 1
      END DO
C
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SCDENS */
      SUBROUTINE SCDENS(Q,SCQ,N,IPRINT)
C***********************************************************************
C
C     INPUT:  Charge of atom
C     OUTPUT: Small component charge
C
C     Small component charges have been obtained from calculaltions
C     using the relativistic atomic code GRASP.
C     see T.SAUE, Ph.D. thesis, University of Oslo 1997
C
C     Quadratic fit:
C       A[0] = 0.0299562717
C       A[1] = -0.00349803581
C       A[2] = 0.000122746761
C       R square = 0.9985796
C     Cubic fit:
C       A[0] = -0.000873131732
C       A[1] = -2.4427764e-05
C       A[2] = 3.96486383e-05
C       A[3] = 5.32680272e-07
C       R square = 0.9999936
C     4th degree fit:
C       A[0] = 0.00137797077
C       A[1] = -0.000439315949
C       A[2] = 5.73988753e-05
C       A[3] = 2.68047049e-07
C       A[4] = 1.27227511e-09
C       R square = 0.999999
C     5th degree fit:
C       A[0] = 0.000348425837
C       A[1] = -0.000161165557
C       A[2] = 3.90375579e-05
C       A[3] = 7.35391124e-07
C       A[4] = -3.77065743e-09
C       A[5] = 1.93958944e-11
C       R square = 0.9999999
C     Written by T.Saue June 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(DTOL = 1.0D-3,D0=0.0D0)
      LOGICAL XPOL,IPOL
      DIMENSION SCCHRG(103),BFIT(6),Q(N),SCQ(N)
#include <sdpre.h>
      DATA (SCCHRG(I),I = 1,103)
     1/ 0.0000D0, 0.0001D0, 0.0002D0, 0.0004D0, 0.0007D0, 
     2  0.0010D0, 0.0014D0, 0.0020D0, 0.0026D0, 0.0034D0, 
     3  0.0043D0, 0.0053D0, 0.0065D0, 0.0077D0, 0.0091D0, 
     4  0.0106D0, 0.0123D0, 0.0141D0, 0.0160D0, 0.0181D0, 
     5  0.0203D0, 0.0227D0, 0.0252D0, 0.0279D0, 0.0308D0, 
     6  0.0339D0, 0.0371D0, 0.0405D0, 0.0440D0, 0.0478D0, 
     7  0.0517D0, 0.0558D0, 0.0602D0, 0.0647D0, 0.0694D0, 
     8  0.0743D0, 0.0793D0, 0.0846D0, 0.0901D0, 0.0958D0, 
     9  0.1017D0, 0.1078D0, 0.1141D0, 0.1206D0, 0.1273D0, 
     O  0.1343D0, 0.1415D0, 0.1489D0, 0.1566D0, 0.1644D0, 
     1  0.1725D0, 0.1809D0, 0.1895D0, 0.1983D0, 0.2073D0, 
     2  0.2166D0, 0.2261D0, 0.2359D0, 0.2459D0, 0.2563D0, 
     3  0.2668D0, 0.2777D0, 0.2888D0, 0.3003D0, 0.3118D0, 
     4  0.3238D0, 0.3360D0, 0.3485D0, 0.3614D0, 0.3746D0, 
     5  0.3879D0, 0.4017D0, 0.4157D0, 0.4301D0, 0.4447D0, 
     6  0.4597D0, 0.4750D0, 0.4907D0, 0.5067D0, 0.5230D0, 
     7  0.5396D0, 0.5566D0, 0.5739D0, 0.5916D0, 0.6097D0, 
     8  0.6281D0, 0.6468D0, 0.6659D0, 0.6854D0, 0.7053D0, 
     9  0.7256D0, 0.7464D0, 0.7675D0, 0.7888D0, 0.8107D0, 
     O  0.8332D0, 0.8560D0, 0.8790D0, 0.9027D0, 0.9268D0, 
     1  0.9514D0, 0.9765D0, 1.0020D0 /
      DATA BFIT/0.000348425837D0,
     &         -0.000161165557D0,
     &          3.90375579D-05,
     &          7.35391124D-07,
     &         -3.77065743D-09,
     &          1.93958944D-11/
C
C     Check that charge is on list
C            
      DO I = 1,N
        IQ   = NINT(Q(I))
        XPOL = IQ.GT.103
        IPOL = (Q(I)-SDPRE(IQ)).GT.DTOL.AND.(.NOT.XPOL)
        IF(XPOL.OR.IPOL) THEN
          SCQ(I) = - POLVAL(5,BFIT,Q(I))
        ELSEIF(IQ.GT.0) THEN
          SCQ(I) = - SCCHRG(IQ)
        ELSE
          SCQ(I) = D0
        ENDIF
        IF(IPRINT.GE.1) THEN
          IF    (XPOL) THEN
            WRITE(LUPRI,'(A,F9.4,A,F6.4,A)')
     &        'Charge: ',Q(I),' --> Small component charge:',SCQ(I),
     &        '(extrapolation)'
          ELSEIF(IPOL) THEN
            WRITE(LUPRI,'(A,F9.4,A,F6.4,A)')
     &        'Charge: ',Q(I),' --> Small component charge:',SCQ(I),
     &        '(interpolation)'
          ELSE
            WRITE(LUPRI,'(A,F9.4,A,F6.4)')
     &        'Charge: ',Q(I),' --> Small component charge:',SCQ(I)
          ENDIF
        ENDIF            
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Lvcorr */
      SUBROUTINE LVCORR(WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Calculate classic repulsive potential of 
C     small component charges
C
C     Written by T.Saue June 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
C
#include <nuclei.h>
      DIMENSION WORK(LWORK)
C
#include <memint.h>
      CALL MEMGET('REAL',KSCQ,NUCIND,WORK,KFREE,LFREE)
      CALL LVCOR1(WORK(KSCQ),IPRINT)
      CALL MEMREL('LVCORR',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Lvcor1 */
      SUBROUTINE LVCOR1(SCQ,IPRINT)
C***********************************************************************
C
C     Calculate classic repulsive potential of 
C     small component charges
C
C     Written by T.Saue June 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
C
#include <nuclei.h>
#include <frame.h>
      DIMENSION SCQ(*)
C
C     Get small component charges
C
      CALL SCDENS(CHARGE,SCQ,NUCIND,IPRINT)
C
C     Calculate nuclear repulsion energy
C
      CORRLV = GETPOT(SCQ)
C
C     Print section
C      
      CALL HEADER('Output from LVCORR',-1)
      WRITE(LUPRI,'(1X,A)')
     &  '* Small component atomic charges:'
      WRITE(LUPRI,'(4X,A4,F10.4)') (NAMN(I),SCQ(I),I = 1,NUCIND)
      WRITE(LUPRI,'(1X,A,F18.12)')
     &  '* Classical repulsion of small component atomic charges:',
     &  CORRLV
      RETURN
C
      END

