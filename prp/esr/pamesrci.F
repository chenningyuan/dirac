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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck esr_init */
      SUBROUTINE ESR_INIT(WF,IPRINT,N_CIESR,SPINPOLESR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Initialize ESR calculation: find out what to do
C
C     Input : WF - wave function type
C
C     Output: N_CIESR - dimension of CI space
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
      LOGICAL   SPINPOLESR
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
         N_CIESR = 2
         GO TO 9999
      END IF
C
C
C     Include spin polarization in CI ?? If yes redefine GAS spaces
C     NASH(:), NASHT, NACTEL, NAELEC etc. are all reset in THDFORB
C
      IF (SPINPOLESR) THEN
         CALL TDHFORB('ESRCI')
         CALL SETDC2(0)
      END IF

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
C     Transfer NZCONF to N_CIESR - dimension of ESR matrix
C
      N_CIESR = NZCONF
      IF (N_CIESR .EQ. 0) CALL QUIT('N_CIESR = 0 !')
C
 9999 CALL QEXIT('ESR_INIT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck esr_resolve */
      SUBROUTINE ESR_RESOLVE(CMO,IPRINT,N_CIROOTS,MAX_CI_IT,THR_CI,
     &           WORK,K_EIGVAL,K_EIGVEC,KFREE,LFREE)

C***********************************************************************
C
C     Resolve CI for ESR with open shell DHF code
C
C     Input : CMO - MO coefficents
C
C     Output: K_EIGVAL - WORK(K_EIGVAL) has eigenvalues
C             K_EIGVEC - WORK(K_EIGVEC) has eigenvectors
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
!#include <krmcluci_inf.h>
C
      DIMENSION CMO(*), WORK(*)
      CHARACTER SZTRA(2,4)*200
C
      CALL QENTER('ESR_RESOLVE')
      KFRSAV = KFREE
C
      IPROPT =  MAX(IPRINT,3)
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
      CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
      CALL IZERO(WORK(KIBEIG),NORBT)
      CALL RTRACTL(CMO,WORK(KIBEIG),SZTRA,WORK(KFREE),LFREE)
      CALL MEMREL('ESR_RESOLVE AFTER RTRACTL',
     &             WORK,1,KFRSAV,KFREE,LFREE) 
C
C     **** Resolve CI problem -- find eigenvalues and eigenvectors
C
      IF (MAX_CI_IT .GT. 0) THEN
C
C        MAX_CI_IT Davidson CI iterations.
C
         L_EIGVAL = N_CIROOTS
         L_EIGVEC = N_CIROOTS*NZCONF*NZ_in_CI ! final CI eigen vectors
         L_EIGVXX = N_CIROOTS*NZCONF*NZ_in_CI ! CI eigen vectors from GASCIP_RCISTD
         CALL MEMGET('REAL',K_EIGVAL,L_EIGVAL,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',K_EIGVXX,L_EIGVEC,WORK,KFREE,LFREE)
!        assign pointer to CI vectors which is in the argument list of this subroutine and needs to be
!        passed upwards to the calling routine.
         K_EIGVEC = K_EIGVXX
      call memchk('before RCISTD',work,1)
         CALL GASCIP_RCISTD(N_CIROOTS,THR_CI,MAX_CI_IT,
     &                     WORK(K_EIGVAL),WORK(K_EIGVXX),NZCONF,
     &                     WORK(KZCONF),CMO,
     &                     IOPT_ISTATE,OPT_CIPROGAM,
     &                     IPROPT, WORK(KFREE),LFREE)
      call memchk('after RCISTD',work,1)
C      SUBROUTINE GASCIP_RCISTD(N_ROOTS,RTHRCI,MAX_CI_IT,
C     &                        ECI,CVECS,NDET,IDET,CMO,
C     &                        ISTATE,CIPROG,IPRINT,
C     &                        WORK,LWORK)

       
C
C        Sort CI eigenvectors from (real C1, imag C1, real C2, etc)
C        to (real C1, real C2, ..., imag C1, imag c2, ...) for later use in QTRANS
C
         IF (NZ_in_CI .EQ. 2) THEN
            CALL MEMGET('REAL',K_EIGVX2,L_EIGVEC,WORK,KFREE,LFREE)
            J_EIGVEC_R = K_EIGVX2 - 1
            J_EIGVEC_I = J_EIGVEC_R + N_CIROOTS*NZCONF
            J_EIGVXX_R = K_EIGVXX - 1
            J_EIGVXX_I = J_EIGVXX_R + NZCONF
            DO J = 1,N_CIROOTS
               DO I = 1,NZCONF
                  WORK(J_EIGVEC_R+I+(J-1)*NZCONF)=
     &            WORK(J_EIGVXX_R+I+(J-1)*NZCONFQ)
                  WORK(J_EIGVEC_I+I+(J-1)*NZCONF)=
     &            WORK(J_EIGVXX_I+I+(J-1)*NZCONFQ)
               END DO
            END DO
            call dcopy(L_EIGVEC,work(K_EIGVX2),1,work(K_EIGVXX),1)

            CALL MEMREL('ESR_RESOLVE AFTER CI STORAGE SORT',
     &               WORK,1,K_EIGVX2,KFREE,LFREE) 

         END IF
C
      ELSE
C
C        Full CI.
C           
         L_EIGVAL = NZCONF
         L_EIGVEC = NZCONF*NZCONF * NZ_in_CI
         CALL MEMGET('REAL',K_EIGVAL,L_EIGVAL,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',K_EIGVEC,L_EIGVEC,WORK,KFREE,LFREE)
C debug TODO remove again: 
         CALL GASCIP_FULLCI(NZCONF,WORK(KZCONF),CMO,
     &                     WORK(K_EIGVAL),WORK(K_EIGVEC),ECORE,
     &                     IPROPT,WORK(KFREE),LFREE)
C         CALL GASCIP_FULLCI(NDET,IDET,CMO,EIGVAL,EIGVEC,ECORE,
C     &                         IPRINT,WORK,LWORK)

      END IF
C           
C
Chj TODO: calculate and save an index vector with symmetry of each CI eigenvector
C
      CALL WRTKRMC(LUKRMC,'CI energ',WORK(K_EIGVAL),L_EIGVAL)
      CALL WRTKRMC(LUKRMC,'CI vecs ',WORK(K_EIGVEC),L_EIGVEC)
      CLOSE(LUKRMC,STATUS='KEEP')

#ifdef ESR_DEBUG
        WRITE(LUPRI,*)
     &   '(ESR_RESOLVE): Real part of the N_ROOTS CI vectors'
        CALL OUTPUT(WORK(K_EIGVEC),
     &              1,NZCONF,1,N_CIROOTS,NZCONF,N_CIROOTS,-1,LUPRI)
      IF (NZ_in_CI .eq. 2) THEN
        WRITE(LUPRI,*)
     &   '(ESR_RESOLVE): Imaginary part of the N_ROOTS CI vectors'
        CALL OUTPUT(WORK(K_EIGVEC+NZCONF*N_CIROOTS),
     &              1,NZCONF,1,N_CIROOTS,NZCONF,N_CIROOTS,-1,LUPRI)
      END IF
      DO I = 1,N_CIROOTS
         WRITE(lupri,*) 'I,WORK(K_EIGVAL-1+1)', I, WORK(K_EIGVAL-1+I) 
      END DO
#endif
C
      CALL QEXIT('ESR_RESOLVE')
C
      RETURN
      END
