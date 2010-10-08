!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
!      DIRAC program, see: http://wiki.chem.vu.nl/dirac
!dirac_copyright_end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck esr_init */
      SUBROUTINE ESR_INIT(WF,IPRINT,NCIESR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Initialize ESR calculation: find out what to do
C
C     Input : WF - wave function type
C
C     Output: NCIESR - dimension of CI space
C
C     Author: Hans Joergen Aa. Jensen
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcborb.h>
#include <dcbdhf.h>
#include <dcbopt.h>
#include <dgroup.h>
C
      DIMENSION WORK(*)
      CHARACTER WF*4
C
      CALL QENTER('ESR_INIT')
C
      IPROPT = IPRINT
      IF (WF(1:3) .NE. 'DHF' .AND. WF .NE. 'DFT') THEN
         CALL QUIT('ESR only implemented for DHF and DFT so far')
      END IF
C
      IF (WF(1:4) .EQ. 'DHF1') THEN
C        ... special input for using old code for a single open shell
         IF (NOPEN .NE. 1) THEN
            CALL QUIT('"DHF1" is only for NOPEN .eq. 1')
         END IF
         NCIESR = 2
         GO TO 9999
      END IF
C
      MZ = MIN(2,NZ)
C
C     Define dcbidx.h
C
      CALL SETDCBIDX()
C
C     Initialize IOPT_*(2) for generation of FC and FV matrices
C     (used in GASCIP_FULLCI for getting FC)
C
      DO I = 1,2
         IOPT_ISYMOP(I) = 1
         IOPT_IFCKOP(I) = 1
         IOPT_IHRMOP(I) = 1
      END DO
C
C     We reuse LUKRMC=50 and RSET_KRMC from krmcscf module
C     for esr_resolve information (but with file name
C     ESR_RESOLVE instead of KRMCSCF).
C
      LUKRMC = 50
      CALL OPNFIL(LUKRMC,'ESR_RESOLVE','UNKNOWN','ESR_INIT')
      CALL NEWLAB('*ESR_RES',LUKRMC,LUPRI)
C
C     **** CI problem is defined by transfer of information in
C          CALL TDHFORB in beginning of PRP module
C
C     **** Find necessary CI information:
C          - find number of determinants,
C          - generate list of determinants.
C          Note: only GASCIP program implemented for ESR_RESOLVE !
C
      IOPT_SYMMETRY = -1
C     ... code for including all symmetries
      CALL RSETCI('GASCIP',WORK,KFREE,LFREE)
C
C     This is CI, so no orbital rotations:
C
      NZXOPE = 0
      NZHOPE = 0
      LZXOPE = NZHOPE*2
      NZXOPP = 0
      LZXOPP = NZXOPP*2
      CALL MEMGET('INTE',KZXOPE,LZXOPE,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KZXOPP,LZXOPP,WORK,KFREE,LFREE)
C
C     Write information to LUKRMC
C
      CALL RSET_KRMC(WORK,KFREE,LFREE)
C
C     Transfer NZCONF to NCIESR - dimension of ESR matrix
C
      NCIESR = NZCONF
      IF (NCIESR .EQ. 0) CALL QUIT('NCIESR = 0 !')
C
 9999 CALL QEXIT('ESR_INIT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck esr_resolve */
      SUBROUTINE ESR_RESOLVE(CMO,IPRINT,WORK,KEIGVL,KEIGVC,KFREE,LFREE)
C***********************************************************************
C
C     Resolve CI for ESR with open shell DHF code
C
C     Input : CMO - MO coefficents
C
C     Output: KEIGVL - WORK(KEIGVL) has eigenvalues
C             KEIGVC - WORK(KEIGVC) has eigenvectors
C
C     Author: Hans Joergen Aa. Jensen
C     (Based on GASCIP_STARTVEC by Joern Thyssen June 2000)
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D0 = 0.0D0)
C
#include <dcbgen.h>
#include <dcborb.h>
#include <dcbopt.h>
#include <dgroup.h>
C
      DIMENSION CMO(*), WORK(*)
      CHARACTER SZTRA(2,4)*200
C
      CALL QENTER('ESR_RESOLVE')
C
      MZ = MIN(2,NZ)
C
C     ***************************************************
C     *** Which integrals to calculate ******************
C     ***************************************************
C
      CALL INTCON(IOPT_INTFLG,IOPT_INTBUF,IOPT_INTDEF,D0,
     &        OPT_CNVINT,9999,IOPT_ITRINT,OPT_INTTYP)
C
      WRITE(LUPRI,9020) OPT_INTTYP
C
 9020 FORMAT(/' (ESR_RESOLVE) Integrals for CI: ',A)
C
C
C     ***************************************************
C     *** Transform two-electron integrals **************
C     ***************************************************
C
      CALL GETTIM(CPUITR1,WALLITR1)
C
      JTRLVL = 0
C     ... only (aa|aa) integrals
      CALL INISZT(JTRLVL,SZTRA,.FALSE.,IDUMMY)
      OPT_NOPFQ = .TRUE.
      CALL RTRACTL(CMO,SZTRA,WORK(KFREE),LFREE)
C
C     **** Resolve CI problem -- find eigenvalues and eigenvectors
C
      N2ZCONF = NZCONF*NZCONF
      CALL MEMGET('REAL',KEIGVL,NZCONF    ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVC,N2ZCONF*MZ,WORK,KFREE,LFREE)
Chj   IF (MAXCIT .GE. 0) THEN
C
C        MAXCIT Davidson CI iterations.
C        -- not implemented yet for ESR_RESOLVE
C
Chj      CALL GASCIP_RCISTD(WORK(KZCONF),CMO,CREF,ECI,
Chj  &        WORK(KFREE),LFREE)
Chj      WRITE(LUPRI,'(/A,F20.8)')
Chj  &        ' (RCIST) CI energy:',ECI
C           
Chj   ELSE
C
C        Full CI.
C           
C debug TODO remove again: 
         IPROPT =  MAX(IPRINT,5)
         CALL GASCIP_FULLCI(WORK(KZCONF),CMO,WORK(KEIGVL),WORK(KEIGVC),
     &                     ECORE,WORK,KFREE,LFREE)
C debug TODO remove again: 
         IPROPT = IPRINT
C        CALL GASCIP_FULLCI(IDET,CMO,EIGVAL,EIGVEC,ECORE,WORK,KFREE,LFREE)
C           
Chj   END IF
C
Chj TODO: calculate and save an index vector with symmetry of each CI eigenvector
C
      CALL WRTKRMC(LUKRMC,'CI energ',WORK(KEIGVL),NZCONF)
      CALL WRTKRMC(LUKRMC,'CI vecs ',WORK(KEIGVC),N2ZCONF*MZ)
      CLOSE(LUKRMC,STATUS='KEEP')
C
      CALL QEXIT('ESR_RESOLVE')
C
      RETURN
      END
