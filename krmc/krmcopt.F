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
C...   group of the licence holder may not take any part of Dalton,
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
C FILE    : diropt.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiopt */
      SUBROUTINE PSIOPT(WORK,LWORK)
C*****************************************************************************
C
C     Driver routine to PSIOP1
C     Get second order optimized wave function.
C     Either DHF or KR-MCSCF...
C
C     Written by J. Thyssen and H. J. Aa. Jensen - Oct 22 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION WORK(LWORK)
C
#include <dcborb.h>
#include <dgroup.h>
#include <dcbopt.h>
#include <dcbgen.h>
C
#include <memint.h>
C
      CALL QENTER('PSIOPT')
C
      IF (NASHT .GT. 1) THEN
         CALL TITLER('KR-MCSCF wave function optimization','*',110)
         WRITE(LUPRI,'(/3X,A/)')
     &      'Kramers Restricted Multiconfigurational '//
     &      'Self-Consistent Field'
      ELSE
         CALL TITLER('QC-DHF wave function optimization','*',110)
         WRITE(LUPRI,'(/3X,A/)')
     &      'Quadratic convergent Dirac-Hartree-fock '//
     &      'wave function optimization'
      END IF
      write(6,*) 'ipropt in PSIOPT',ipropt
C
C     ***************************************************
C     *** Open files needed *****************************
C     ***************************************************
C
      LUKRMC = 20
      CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','PSIOPT')
      CALL NEWLAB('SODLABEL',LUKRMC,LUERR)
      REWIND LUKRMC
C
      LUKRM1 = 101
      CALL OPNFIL(LUKRM1,'KRMC KAPPA','UNKNOWN','PSIOPT')
C
      LUKRM2 = 102
      CALL OPNFIL(LUKRM2,'KRMC CMO','UNKNOWN','PSIOPT')
      CALL NEWLAB('SODLABEL',LUKRM2,LUERR)
      REWIND LUKRM2
C
      LUKRM3 = 103
      CALL OPNFIL(LUKRM3,'KRMC FOCK','UNKNOWN','PSIOPT')
      CALL NEWLAB('SODLABEL',LUKRM3,LUERR)
      REWIND LUKRM3
C
      LUITFO = 104
      CALL OPNFIL(LUITFO,'KRMC ITINFO','UNKNOWN','PSIOPT')
C
C
C     ***************************************************
C     *** Initialize configurational and orbital      ***
C     *** parameters.                                 ***
C     ***************************************************
C
C
      CALL RSETWOP(WORK,KFREE,LFREE)
C
C
C     ***************************************************
C     *** Memory allocation *****************************
C     ***************************************************
C
C
      CALL MEMGET('REAL',KCMO  ,NCMOTQ        ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDV   ,N2ASHXQ       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC   ,N2ORBXQ       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFV   ,N2ORBXQ       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFQ   ,NZ*NASHT*NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAC ,N2ASHXQ       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCREF ,NZCONFQ       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2AC ,N2ASHXQ*N2ASHT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGCI  ,NZCONFQ       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGOEE ,NZXOPEQ       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGOEP ,NZXOPPQ       ,WORK,KFREE,LFREE)
      IF ( OPT_CIPROGRAM .EQ. 'GOSCI' ) THEN
         LENPV = ( 2 * NASHT ) ** 4 * MIN0(NZ,2)
      ELSE
         IF ( NASHT .GE. 1 ) THEN
            WRITE(LUPRI,'(/,A,/,A,I3,2A)')
     $           '*** ERROR in PISOPT ***',
     $           'NASHT = ',NASHT,' not programmed for CI program ',
     $           OPT_CIPROGRAM
            CALL QUIT('*** ERROR in PSIOPT ***')
         END IF
         LENPV = 0
      END IF

C
C     ***************************************************
C     *** Call 2nd order optimization *******************
C     ***************************************************
C
      CALL PSIOP1(WORK(KCMO),WORK(KDV),WORK(KFC),
     &            WORK(KFV),WORK(KFQ),WORK(KFCAC),WORK(KCREF),
     &            WORK(KH2AC),WORK(KGCI),WORK(KGOEE),WORK(KGOEP),
     &            WORK,KFREE,LFREE)
C
C
C     ***************************************************
C     *** Memory deallocation ***************************
C     ***************************************************
C
      CALL MEMREL('PSIOPT',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
C     ***************************************************
C     *** Call output routines **************************
C     ***************************************************
C
C
      CALL KRMCOUT(WORK(KFREE),LFREE)
C
C
C     ***************************************************
C     *** Close files ***********************************
C     ***************************************************
C
      CLOSE(LUKRMC,STATUS='KEEP')
      CLOSE(LUKRM1,STATUS='KEEP')
      CLOSE(LUKRM2,STATUS='KEEP')
      CLOSE(LUKRM3,STATUS='KEEP')
      CLOSE(LUITFO,STATUS='KEEP')
c     CLOSE(LUKRM1,STATUS='DELETE')
c     CLOSE(LUKRM2,STATUS='DELETE')
c     CLOSE(LUKRM3,STATUS='DELETE')
c     CLOSE(LUITFO,STATUS='DELETE')
C
      CALL QEXIT('PSIOPT')
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiopt */
      SUBROUTINE PSIOP1(CMO,DV,FC,FV,FQ,FCAC,CREF,H2AC,GCI,GOEE,GOEP,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Get second order optimized wave function.
C     Either DHF or KR-MCSCF...
C
C     Written by J. Thyssen and H. J. Aa. Jensen - Oct 22 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( D0 = 0.0D00, D1 = 1.0D00 )
C
#include <dcbpsi.h>
#include <dcbkrmc_itinfo.h>
#include <dcbgen.h>
#include <dcbopt.h>
#include <dcborb.h>
#include <frame.h>
C
      DIMENSION WORK(*)
      DIMENSION CMO(*),DV(*),FC(*),FV(*),FQ(*),FCAC(*),CREF(*),
     &          H2AC(*),GCI(*),GOEE(*),GOEP(*)
      LOGICAL   CHKSTP,DO4INDEX
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
C
C
      CALL QENTER('PSIOP1')
      CALL GETTIM(CPUTOT1,WALLTOT1)
      KFRSAV = KFREE
CDEBUG OPTIONS
      OPT_THRGRD = 1.0D-7
CEND DEBUG
      write(6,*) 'iopt_intdef psiop1',iopt_intdef,ipropt
C
C
C     ***************************************************
C     *** Initialize some parameters ********************
C     ***************************************************
C
C
      DO I = 1,2
         IOPT_ISYMOP(I) = 1
         IOPT_IFCKOP(I) = 1
         IOPT_IHRMOP(I) = 1
      END DO
      ITMAC  = 0
      KTRLVL = -1
      JTRLVL = 0
      DO4INDEX = ( NASHT .GT. 1 )
      EMCOLD = D0
      DEPRED = D0
      STPLEN = D0
      CALL DZERO(GNORM,5)
C
C     Integrals in macro iteration zero.
C
      CALL INTCON(IOPT_INTFLG,IOPT_INTBUF,IOPT_INTDEF,GNORM(5),
     &     OPT_CNVINT,ITMAC,IOPT_ITRINT,OPT_INTTYP)
C
C
C     ***************************************************
C     *** Get start vectors *****************************
C     ***************************************************
C
C
#ifdef UNDEF
      write(6,*) '(1) iopt_intflg=',iopt_intflg
      write(6,*) '(1) iopt_intbuf=',iopt_intflg
      write(6,*) '(1) iopt_intdef=',iopt_intflg
      write(6,*) '(1) gnorm(5)=',gnorm(5)
      write(6,*) '(1) opt_cnvint=',opt_cnvint
      write(6,*) '(1) itmac=',itmac
      write(6,*) '(1) iopt_itrint=',iopt_itrint
      write(6,*) '(1) opt_inttyp=',opt_inttyp
#endif
C
C     Get start orbitals and CI vectors:
C
      CALL ROPTST(WORK(KFREE),LFREE)
C
C     Initialize some more parameters:
C
C     If DEPRED < 0, we have restarted a calculation.
C
      IF (DEPRED .NE. D0) THEN
         CHKSTP = .TRUE.
         ISTEP  = -1
      ELSE
         CHKSTP = .FALSE.
         ISTEP  = 0
      END IF
C
C     Check for legal values of Beta:
C
      BETA = ABS(BETA)
      IF (BETA .LE. EPSILON)  BETA = D1
      IF (BETMIN .LE. D0 .OR.
     &    BETMIN .GT. BETMAX) BETMIN = D1
      IF (BETMAX .LE. D1)     BETMAX = D1
C
C
C     ***************************************************
C     *** Start of macro iteration loop *****************
C     ***************************************************
C
C
  100 CONTINUE
         CALL DZERO(DINFO,LDINFO)
         CALL IZERO(IINFO,LIINFO)
         ITMAC = ITMAC+1
         CALL GETTIM(CPUMAC1,WALLMAC1)
         WRITE(LUPRI,9000) ITMAC
C
 9000    FORMAT(//' <<< MACRO ITERATION',I3,' >>>'
     &           /' ---------------------------')
         CALL FLSHFO(LUPRI)
C
C        ***************************************************
C        *** Which integrals to calculate ******************
C        ***************************************************
C
         CALL INTCON(IOPT_INTFLG,IOPT_INTBUF,IOPT_INTDEF,GNORM(5),
     &        OPT_CNVINT,ITMAC,IOPT_ITRINT,OPT_INTTYP)
C
         WRITE(LUPRI,9001) OPT_INTTYP
C
 9001    FORMAT(/1X,'* Integrals in this macro iteration: ',A)
C
C        ***************************************************
C        *** Get MO-coefficients ***************************
C        ***************************************************
C
C
         JRDCMO = 0
         CALL RREADMO(CMO,JRDCMO,1)
C
C
C        ***************************************************
C        *** Transform two-electron integrals **************
C        ***************************************************
C
         CALL GETTIM(CPUITR1,WALLITR1)
         IF ( DO4INDEX ) THEN
C
C           Do 4-index transformation.
C           KTRLVL is the order of the transformation.
C
            JTRLVL = KTRLVL
            CALL RTRACTL(JTRLVL,WORK(KFREE),LFREE)
C
         END IF
         CALL GETTIM(CPUITR2,WALLITR2)
C
C
C
C        ***************************************************
C        *** Get Gradient **********************************
C        ***************************************************
C
C
         CALL RGRAD(CMO,WORK(KZXOPE),WORK(KZXOPP),
     &              GCI,GOEE,GOEP,EMY,EACTIV,DV,FC,FV,FQ,
     &              FCAC,H2AC,WORK(KFREE),LFREE,IPROPT)
         CALL RGRDINF(GCI,GOEE,GOEP)
C
C
         EMCSCF = POTNUC + EMY + EACTIV
         WRITE(LUPRI,9005) POTNUC,EMY,EACTIV,EMCSCF
C
 9005    FORMAT(/' Nuclear repulsion            : ',F25.15,/,
     &           ' Inactive energy              : ',F25.15,/,
     &           ' Active energy                : ',F25.15,/,
     &           ' Total KR-MCSCF energy        : ',F25.15)
C
         WRITE(LUPRI,9010) (ITMAC,GNORM(I),I=1,5)
         CALL FLSHFO(LUPRI)
C
 9010    FORMAT(/' CI gradient        (It. ',I3,') : ',F25.15,/,
     &           ' Orbital e-e grad.  (It. ',I3,') : ',F25.15,/,
     &           ' Orbital e-p grad.  (It. ',I3,') : ',F25.15,/,
     &           ' Total orbital grad.(It. ',I3,') : ',F25.15,/,
     &           ' Total gradient     (It. ',I3,') : ',F25.15)
C
C
C        ***************************************************
C        *** Converged ??? *********************************
C        ***************************************************
C
         write(6,*) 'jth: THRGRD',OPT_THRGRD
         IF (GNORM(5) .LE. OPT_THRGRD) THEN
            ICONV = 1
            ISTEP = -2
         END IF
C
C
C        ***************************************************
C        *** Check step ************************************
C        ***************************************************
C
         IF (CHKSTP) THEN
            CALL DIRSTP(ISTEP,CMO,WORK(KFREE),LFREE)
         ELSE
            EMCOLD = EMCSCF
            ISTEP  = 0
         END IF
C
         DINFO(9) = RTRUST
         IF (ICONV .EQ. 1) GOTO 1000
C
         CHKSTP = .TRUE.
C
C
C
C        ***************************************************
C        *** Second order step *****************************
C        ***************************************************
C
C
         CALL DIRNR(CREF,GCI,GOEE,GOEP,CMO,DV,FC,FV,FCAC,
     &              H2AC,WORK,KFREE,LFREE)
C
C
C
C        ***************************************************
C        *** Step accepted: save info on LUKRMC ************
C        ***************************************************
C
C
         CALL KRMCSAV(0,WORK(KFREE),LFREE)
C
C
 1000    CONTINUE
         CALL GETTIM(CPUMAC2,WALLMAC2)
         CPUTID = SECTID(CPUMAC2-CPUMAC1)
         WALLTID = SECTID(WALLMAC2-WALLMAC1)
         WRITE(LUPRI,9100) CPUTID,WALLTID
C
C        Save information for final summary output
C
         IINFO(1) = ITMAC
         DINFO(1) = EMY
         DINFO(2) = EACTIV
         DINFO(3) = EMCSCF
         DINFO(16) = CPUMAC2-CPUMAC1
         DINFO(18) = CPUITR2-CPUITR1
         WRITE(LUITFO) DINFO,IINFO
C
 9100    FORMAT(/'>>> CPU (WALL) TIME FOR THIS MACRO IT. :',A,'(',A,')')
         write(6,*) 'iconv=',iconv
         write(6,*) 'itmac,mopt_mxmacro=',itmac,mopt_mxmacro
         CALL FLSHFO(LUPRI)
         IF (ICONV .GT. 0) GOTO 2000 
         IF (ITMAC .LT. MOPT_MXMACRO) THEN
            IOPT_INTBUF = IOPT_INTFLG
            GOTO 100
         END IF
C
C
C        ***************************************************
C        *** End of macro iteration loop *******************
C        ***************************************************
C
C
C
C
C
C     ***************************************************
C     *** Not converged *********************************
C     ***************************************************
C
C
      WRITE(LUPRI,9200) 
 9200 FORMAT(/'WARNING: Wave function did not converge...')
      GOTO 3000
 2000 CONTINUE
C
C
C     ***************************************************
C     *** We have converged *****************************
C     ***************************************************
C
C
 3000 CONTINUE
      CALL GETTIM(CPUTOT2,WALLTOT2)
      CPUTID = SECTID(CPUTOT2-CPUTOT1)
      WALLTID = SECTID(WALLTOT2-WALLTOT1)
      WRITE(LUPRI,9300) CPUTID,WALLTID
C
 9300 FORMAT(/'>>> TOTAL CPU (WALL) TIME IN KR-MCSCF: ',A,'(',A,')')
C
      CALL MEMREL('PSIOP1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PSIOP1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgrad */
      SUBROUTINE RGRAD(CMO,JXOPE,JXOPP,
     &                 GCI,GOEE,GOEP,EMCMY,EMCACT,DV,FC,FV,FQ,FCAC,
     &                 H2AC,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Calculate the gradient.
C
C     Input:
C        CMO     - the orbitals for the CEP
C        CREF    - CI vector
C        JXOPE   - orbital e-e rotations
C        JXOPP   - orbital e-p rotations
C
C
C     Output:
C        EMCACT  - active energy
C        GCI     - the CI gradient of the CEP
C        GOEE    - the e-e orbital gradient of the CEP
C        GOEP    - the e-p gradient of the CEP
C        EMCMY   - inactive energy 
C        DV      - one-electron active density matrix
C        H2AC    - two-electron integrals with active indices
C        FC      - the inactive Fock matrix
C        FCAC    - the inactive Fock matrix with active indices
C        FV      - the active Fock matrix
C        FQ      - the Q Fock matrix
C     
C
C     Written by J. Thyssen and H. J. Aa. Jensen - Oct 22 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D0 = 0.0D00, D1 = 1.0D00, DM1 = -D1)
C
      DIMENSION GCI(*),DV(*),FC(*),FV(*),FQ(*)
      DIMENSION JXOPE(*),JXOPP(*),GOEE(*),GOEP(*)
      DIMENSION WORK(LWORK)
C
C     Local variable
C
      DIMENSION DFAC(2)
C
      LOGICAL   FND,FNDLAB
C
#include <dcbopt.h>
#include <dcbgen.h>
#include <dcborb.h>
#include <dgroup.h>
C
#include <memint.h>
C
      CALL QENTER('RGRAD')
C
      DFAC(1) = D1
      DFAC(2) = DM1
C
C     ***************************************************
C     *** Clean file KRMC FOCK **************************
C     ***************************************************
C
      REWIND LUKRM3 
      FND =  FNDLAB('SODLABEL',LUKRM3)
      BACKSPACE LUKRM3
      CALL NEWLAB('SODLABEL',LUKRM3,LUERR)
      REWIND(LUKRM3)
C
C     ***************************************************
C     *** Calculate active density matrices *************
C     ***************************************************
C
      IF (NASHT .NE. 0) THEN
         CALL RMAKDM(DV,WORK(KFREE),LFREE)
      END IF
C
C     ***************************************************
C     *** Calculate inactive and active Fock matrices ***
C     ***************************************************
C
      CALL RFCKMAT(DV,CMO,EMCMY,EMCACT,FC,FV,WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** Calculate FQ matrix and two-electron inte-  ***
C     *** grals with active orbital indices.          ***
C     ***************************************************
C
      IF (NASHT .GE. 1) THEN
C
C        Get Fock matrix with both indeces active.
C        FCAC is written to KRMC FOCK.
C
         CALL RGETAC(FC,FCAC,IPROPT)
C
      END IF
C
      IF (NASHT .GT. 1) THEN
C
C        Calculate FQ matrix and two-electron integrals
C        with active orbital indeces (H2AC).     
C        FQ and H2AC are written to KRMC FOCK.
C
         CALL RGETH2(DV,FQ,WORK(KFREE),LFREE)
C
C        Active energy.
C
         EMCACT = RENRACT(DV,FCAC,FQ,IPROPT)
C
      ELSE IF (NASHT. EQ. 1) THEN
         CALL DZERO(FQ,NZ*NASHT*NORBT)
         EMCACT = RENRACT(DV,FCAC,FQ,IPROPT)
      ELSE
         EMCACT = D0
      END IF
C
C
C     ***************************************************
C     *** Calculate Fock matrix contributions to the  ***
C     *** diagonal of the orbital Hessian             ***
C     ***************************************************
C
C
      IF (NZXOPT .GT. 1) THEN
c        CALL RDIAG
      END IF
C
C
C     ***************************************************
C     *** Calculate the orbital gradient ****************
C     ***************************************************
C
C
      IF (NZXOPE .GT. 0) THEN
         CALL DZERO(GOEE,NZXOPEQ)
         LFPV = NZ*NASHT*NORBT
         write(6,*) 'jth: lfree before rorbgrd',LFREE
         CALL MEMGET('REAL',KFPV,LFPV,WORK,KFREE,LFREE)
         CALL RORBGRD(DV,FC,FV,FQ,GOEE,
     &                JXOPE,NZXOPE,WORK(KFPV),DFAC,
     &                IPQTOQ(1,0),IPRINT)
         CALL MEMREL('RORBGRD.KFPV',WORK,KFPV,KFPV,KFREE,LFREE)
         write(6,*) 'jth: lfree after rorbgrd',LFREE
      END IF
      IF (NZXOPP .GT. 0) THEN
         CALL DZERO(GOEP,NZXOPPQ)
         LFPV = NZ*NASHT*NORBT
         write(6,*) 'jth: lfree before rorbgrd',LFREE
         CALL MEMGET('REAL',KFPV,LFPV,WORK,KFREE,LFREE)
         CALL RORBGRD(DV,FC,FV,FQ,GOEP,
     &                JXOPP,NZXOPP,WORK(KFPV),DFAC,
     &                IPQTOQ(1,0),IPRINT)
         CALL MEMREL('RORBGRD.KFPV',WORK,KFPV,KFPV,KFREE,LFREE)
         write(6,*) 'jth: lfree after rorbgrd',LFREE
      END IF
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('rGRAD: gradient:',-1)
         CALL RPRKAP(GOEE,GOEP,JXOPE,NZXOPE,
     &               JXOPP,NZXOPP,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C
C     ***************************************************
C     *** Remove redundant orbital rotations ************
C     ***************************************************
C
C
C     se REDPAR i ./prp/pamset.F
C
C
C     ***************************************************
C     *** Calculate the CI gradient *********************
C     ***************************************************
C
      IF (NZCONF .GT. 0) THEN
c        CALL RCIDIAG
c        CALL RCIGRAD
      END IF
C
C
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('RGRAD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rmakdm */
      SUBROUTINE RMAKDM(DV,WORK,LWORK)
C*****************************************************************************
C
C     Construct one- and two-electron active density matrix
C
C     Input:
C        CMO     - the orbitals for the CEP
C        CREF    - CI vector
C
C     Output:
C        DV      - one electron active density matrix
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D00, D0 = 0.0D00, DEPS = 1.0D-10 )
C
#include <dcborb.h>
#include <dcbgen.h>
#include <dcbopt.h>
#include <dgroup.h>
#include <dcbqal.h>
C
      DIMENSION DV(*)
      DIMENSION WORK(*)
      DIMENSION DTRCE(4)
C
C
#include <memint.h>
C
      CALL QENTER('RMAKDM')
C
      IF ( NASHT .GT. 1 )
     $     CALL MEMGET('REAL',KPV,LENPV,WORK,KFREE,LFREE)
      write(6,*) 'lenpv ',lenpv
C
      IF ( OPT_CIPROGRAM .EQ. 'GOSCI' ) THEN
C
         IF ( NASHT .EQ. 1 ) THEN
            DV(1) = D1
            DO IZ = 2, NZ
               DV(IZ) = D0
            END DO
         ELSE
C
C           Calculate density matrices
C
            CALL MEMGET('REAL',KDVMOLF,(2*NASHT)**2*MIN0(2,NZ),
     $           WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KCREF,NZCONFQ,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KIKR,(2*NASHT)*2,WORK,KFREE,LFREE)
            CALL REAKRMC(LUKRMC,'CISTART0',WORK(KCREF),NZCONFQ)
            CALL GOSCI_RMAKDM(WORK(KCREF),WORK(KDVMOLF),WORK(KPV),
     $           WORK(KIKR))
            CALL MFC2QFC(WORK(KDVMOLF),DV,IPROPT)
            CALL MEMREL('RMAKDM.GOSCI',WORK,KDVMOLF,KDVMOLF,KFREE,LFREE)
C
C           Check trace of DV
C
            DO IZ = 1, NZ
               DTRCE(IZ) = D0
            END DO
            DO I = 1,NASHT
               DO IZ = 1, NZ
                  DTRCE(IZ) = DTRCE(IZ) +
     $                 DV( I + ( I - 1 ) * NASHT +
     $                 ( IZ - 1 ) * NASHT*NASHT )
               END DO
            END DO
            IF ( IPROPT .GE. 10 ) THEN
               CALL HEADER('RMAKDM: trace of DV',-1)
               WRITE(LUPRI,'(A,4X,3A)')
     $              'Matrix','    Actual    ','   Expected   ',
     $              '     Diff     '
               DO IZ = 1, NZ
                  IF ( IZ .EQ. 1 ) THEN
                     DTREXP = NAELEC
                  ELSE
                     DTREXP = D0
                  END IF
                  WRITE(LUPRI,'(1P,A1,A,6X,D10.3,4X,D10.3,4X,D10.3)')
     $                 QUNIT(IPQTOQ(IZ,0)),'-part',
     $                 DTRCE(IZ),DTREXP,ABS(DTRCE(IZ)-DTREXP)
                  IF ( ABS ( DTRCE(IZ) - DTREXP ) .GE. DEPS )
     $                 CALL QUIT('*** ERROR in RMAKDM ***')
               END DO
            END IF
         END IF
      ELSE
         IF (NASHT .EQ. 1) THEN
            DV(1) = D1
            DO IZ = 2,NZ
               DV(IZ) = D0
            END DO
         END IF
      END IF
C
      IF ( IPROPT .GE. 5 ) THEN
         CALL HEADER('RMAKDM: DV density matrix',-1)
         CALL PRQMAT(DV,NASHT,NASHT,NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Save active density matrices on file KRMC FOCK
C
      CALL WRTKRMC(LUKRM3,'DVMO    ',DV,N2ASHXQ)
      IF ( NASHT .GT. 1 ) THEN
         CALL WRTKRMC(LUKRM3,'PVMO    ',WORK(KPV),LENPV)
      END IF
      CALL MEMREL('RMAKDM.PV',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RMAKDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rfckmat */
      SUBROUTINE RFCKMAT(DV,CMO,EMCMY,EMCACT,FC,FV,WORK,LWORK)
C*****************************************************************************
C
C     Construct active and inactive Fock matrices
C
C     Input:
C        CMO     - the orbitals for the CEP
C        DV      - the active one electron density matrix
C        EMCACT  - the energy of CEP
C
C     Output:
C
C        EMCMY   - inactive energy 
C        FC      - the inactive Fock matrix
C        FV      - the active Fock matrix
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <aovec.h>
C
      PARAMETER (DP5 = 0.50D00, D2 = 2.00D00, D1 = 1.00D0)
C
#include <blocks.h>
#include <dcborb.h>
#include <dcbopt.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dgroup.h>
      DIMENSION DV(*),CMO(*),FC(*),FV(*)
      DIMENSION WORK(LWORK)
      LOGICAL NOFV
      DIMENSION IQBF(4)
C
C
#include <memint.h>
C
      CALL QENTER('RFCKMAT')
C
      CALL MEMGET('REAL',KH1AO,N2BBASXQ,WORK,KFREE,LFREE)
      IF (NASHT .GT. 0) THEN
         IF (NISHT .GT. 0) THEN
C
            CALL MEMGET('REAL',KFCAO,2*N2BBASXQ,WORK,KFREE,LFREE)
            KFAO  = KFCAO
            KFVAO = KFCAO + N2BBASXQ
C
            CALL MEMGET('REAL',KDCAO,2*N2BBASXQ,WORK,KFREE,LFREE)
            KDAO  = KDCAO
            KDVAO = KDCAO + N2BBASXQ
C
            NFMAT = 2
C
         ELSE
C
            CALL MEMGET('REAL',KFCAO,2*N2BBASXQ,WORK,KFREE,LFREE)
            KFVAO = KFCAO + N2BBASXQ
            KFAO  = KFVAO
C
            CALL MEMGET('REAL',KDVAO,N2BBASXQ,WORK,KFREE,LFREE)
            KFCAO = KDVAO
            KDAO  = KDVAO
C
            NFMAT = 1
         END IF
      ELSE
         IF (NISHT .GT. 0) THEN
C
            CALL MEMGET('REAL',KFCAO,N2BBASXQ,WORK,KFREE,LFREE)
            KFAO = KFCAO
C
            CALL MEMGET('REAL',KDCAO,N2BBASXQ,WORK,KFREE,LFREE)
            KDAO = KDCAO
C
            NFMAT = 1
C
         ELSE
            CALL QUIT('no inactive nor active orbitals...')
         END IF
      END IF
      NPOS = 0
      IF (PARCAL) THEN
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NSMLBL*(NSMLBL+1)/2
      END IF
      CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
C
C     ***************************************************
C     *** Get one-electron density matrices *************
C     ***************************************************
C
C     Construct DCAO
C
      IF (NISHT .GT. 0) THEN
         CALL DENSTY(WORK(KDCAO),CMO,IPROPT)
      END IF
C
C     Construct DVAO
C
      IF (NASHT .GT. 0) THEN
         CALL DENSDV(WORK(KDVAO),CMO,DV,IPROPT,WORK(KFREE),LFREE)
C
C        Scale DVao with 0.5 due to scaling with 2.0 in TWOFCK.
C        (TWOFCK assumes inactive density matrices normalized to 1)
C
         CALL DSCAL(N2BBASXQ,DP5,WORK(KDVAO),1)
      END IF
C
C     ***************************************************
C     *** Get 2-electron Fock matrices ******************
C     ***************************************************
C
      CALL TWOFCK(IOPT_ISYMOP,IOPT_IHRMOP,IOPT_IFCKOP,
     &            WORK(KFAO),WORK(KDAO),
     &            NFMAT,WORK(KPOS),IOPT_INTFLG,IPROPT,WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** Get 1-electron Fock matrix ********************
C     ***************************************************
C
      CALL DIRH1(WORK(KH1AO),IPROPT,WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** Calculate inactive energy *********************
C     ***************************************************
C
      IF (NISHT .GT. 0) THEN
C
C        1-electron inactive energy
C        (factor D2: DCao normed to 1 (not 2))
C
         EMCMY = D2*DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KH1AO),1)
         IF (IPROPT .GE. 4) 
     &      WRITE(LUPRI,9000) 1,EMCMY
C
C        Add 2-electron inactive energy
C        (factor D2: DCao normed to 1 (not 2))
C
         EMCMY2E = DP5 * D2 * DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KFCAO),1)
         IF (IPROPT .GE. 4) 
     &      WRITE(LUPRI,9000) 2,EMCMY2E
         EMCMY = EMCMY + EMCMY2E
         IF (IPROPT .GE. 4) 
     &      WRITE(LUPRI,9001) EMCMY
      ELSE
         EMCMY = D0
      END IF
 9000 FORMAT (/2X'rFCKMAT: ',I1,'-electron inactive energy: ',F25.15)
 9001 FORMAT (/2X'rFCKMAT: Total inactive energy     : ',F25.15)
C
C     ***************************************************
C     *** Add one-electron contribution to FCAO *********
C     ***************************************************
C
      CALL DAXPY(N2BBASXQ,D1,WORK(KH1AO),1,WORK(KFCAO),1)
C
      IF ( IPROPT .GE. 5 ) THEN
         WRITE(LUPRI,9005)
         IF ( NISHT .GT. 0 ) THEN
            WRITE(LUPRI,9006) 'FCao'
            CALL PRQMAT(WORK(KFCAO),NASHT,NASHT,NASHT,NASHT,
     $           NZ,IPQTOQ(1,0),LUPRI)
         END IF
         IF ( NASHT .GT. 0 ) THEN
            WRITE(LUPRI,9006) 'FVao'
            CALL PRQMAT(WORK(KFCAO),NASHT,NASHT,NASHT,NASHT,
     $           NZ,IPQTOQ(1,0),LUPRI)
         END IF
      END IF
 9005 FORMAT(/'rFCKMAT: AO Fock matrices')
 9006 FORMAT(/'rFCKMAT: ',A,' fock matrix')
C
C     ***************************************************
C     *** Transform to MO basis *************************
C     ***************************************************
C
      CALL DZERO(FC,N2ORBXQ)
      IF (NASHT .GT. 0) CALL DZERO(FV,N2ORBXQ)
      DO 100 I = 1,NFSYM
         IF (NORB(I) .EQ. 0) GOTO 100
         IF (IPROPT .GE. 5) WRITE(LUPRI,9010) I
         CALL QTRANS('AOMO','S',D0,
     &               NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &               WORK(KFCAO+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),
     &               FC(I2ORBX(I,I)+1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               WORK(KFREE),LFREE,IPROPT)
         IF (IPROPT .GE. 5) THEN
            WRITE(LUPRI,9011) 'FC',I
            CALL PRQMAT(FC(I2ORBX(I,I)+1),NORB(I),NORB(I),
     &                  NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
         IF (NASHT .EQ. 0) GOTO 100
         CALL QTRANS('AOMO','S',D0,
     &               NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &               WORK(KFVAO+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),
     &               FV(I2ORBX(I,I)+1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               WORK(KFREE),LFREE,IPROPT)
         IF (IPROPT .GE. 5) THEN
            WRITE(LUPRI,9011) 'FV',I
            CALL PRQMAT(FV(I2ORBX(I,I)+1),NORB(I),NORB(I),
     &                  NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
  100 CONTINUE
C
C     Write to file KRMC FOCK
C
c      write(6,*) 'writing files to KRMC FOCK'
      CALL WRTKRMC(LUKRM3,'FCMO    ',FC,N2ORBXQ)
      CALL WRTKRMC(LUKRM3,'FVMO    ',FV,N2ORBXQ)
c      call getmat('FCMO',FC,10,.TRUE.,WORK(KFREE),LFREE)
c      call prqmat(fc,norbt,norbt,norbt,norbt,nz,ipqtoq(1,0),lupri)
c      call getmat('FVMO',FV,10,.TRUE.,WORK(KFREE),LFREE)
c      call prqmat(fv,norbt,norbt,norbt,norbt,nz,ipqtoq(1,0),lupri)
C
C
C
      CALL MEMREL('RFCKMAT',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('RFCKMAT')
      RETURN
 9010 FORMAT(/'rFCKMAT: MO Fock matrices, irrep ',I1)
 9011 FORMAT(/'rFCKMAT: ',A,' fock matrix, irrep ',I1)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgetfc */
      SUBROUTINE RGETFC(CMO,FC,ECORE,WORK,LWORK)
C*****************************************************************************
C
C     Construct inactive Fock matrices
C
C     Input:
C        CMO     - the orbitals for the CEP
C
C     Output:
C
C        FC      - the inactive Fock matrix
C
C     Written by J. Thyssen - June 30 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <aovec.h>
C
      PARAMETER (DP5 = 0.50D00, D2 = 2.00D00, D1 = 1.00D0)
C
#include <blocks.h>
#include <dcborb.h>
#include <dcbopt.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <frame.h>
      DIMENSION CMO(*),FC(*)
      DIMENSION WORK(LWORK)
      DIMENSION IQBF(4)
C
C
#include <memint.h>
C
      CALL QENTER('RGETFC')
C
C
C     ***************************************************
C     *** Get 1-electron Fock matrix ********************
C     ***************************************************
C
      CALL MEMGET('REAL',KH1AO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL DIRH1(WORK(KH1AO),IPROPT,WORK(KFREE),LFREE)
C
      IF (NISHT .GT. 0) THEN
C
         CALL MEMGET('REAL',KFCAO,N2BBASXQ,WORK,KFREE,LFREE)
C
         CALL MEMGET('REAL',KDCAO,N2BBASXQ,WORK,KFREE,LFREE)
C
         NFMAT = 1
C
         NPOS = 0
         IF (PARCAL) THEN
            NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
            NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
            NPOS = NPOS + NSMLBL*(NSMLBL+1)/2
         END IF
         CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
C     
C        ***************************************************
C        *** Get one-electron density matrices *************
C        ***************************************************
C
C        Construct DCAO
C
         CALL DENSTY(WORK(KDCAO),CMO,IPROPT)
C     
C        ***************************************************
C        *** Get 2-electron Fock matrices ******************
C        ***************************************************
C
         CALL TWOFCK(IOPT_ISYMOP,IOPT_IHRMOP,IOPT_IFCKOP,
     &        WORK(KFCAO),WORK(KDCAO),
     &        NFMAT,WORK(KPOS),IOPT_INTFLG,IPROPT,WORK(KFREE),LFREE)
C
C
C        1-electron inactive energy
C        (factor D2: DCao normed to 1 (not 2))
C
         ECORE = D2*DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KH1AO),1)
         IF (IPROPT .GE. 4) 
     &        WRITE(LUPRI,9000) 1,ECORE
C
C        Calculate 2-electron inactive energy
C        (factor D2: DCao normed to 1 (not 2))
C     
         ECORE2E = DP5 * D2 * DDOT(N2BBASXQ,WORK(KDCAO),1,WORK(KFCAO),1)
         IF (IPROPT .GE. 4) 
     &        WRITE(LUPRI,9000) 2,ECORE2E
         ECORE = ECORE + EMCMY2E
C
C        ***************************************************
C        *** Add two-electron contribution to H1AO *********
C        ***************************************************
C     
         CALL DAXPY(N2BBASXQ,D1,WORK(KFCAO),1,WORK(KH1AO),1)
C
      ELSE
         EMCMY = D0
      END IF
C
      IF (IPROPT .GE. 4) 
     &     WRITE(LUPRI,9001) POTNUC
      ECORE = ECORE + POTNUC
      IF (IPROPT .GE. 4) 
     &     WRITE(LUPRI,9002) ECORE
C
      IF ( IPROPT .GE. 5 ) THEN
         WRITE(LUPRI,9006) 'FCao'
         CALL PRQMAT(WORK(KH1AO),NASHT,NASHT,NASHT,NASHT,
     $        NZ,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
      END IF
 9000 FORMAT (/2X'rRGETFC: ',I1,'-electron core energy: ',F25.15)
 9001 FORMAT (/2X'rRGETFC: Nuclear repulsion     : ',F25.15)
 9002 FORMAT (/2X'rRGETFC: Total core energy     : ',F25.15)
 9006 FORMAT(/'RGETFC: ',A,' fock matrix')
C
C     ***************************************************
C     *** Transform to MO basis *************************
C     ***************************************************
C
      CALL DZERO(FC,N2ORBXQ)
      DO 100 I = 1,NFSYM
         write(6,*) 'rGETFC: norb = ',i,norb(i),nfbas(i,0),i2basx(i,i),
     &        ntbas(0),i2orbx(i,i),norbt,icmoq(i)
         IF (NORB(I) .EQ. 0) GOTO 100
         CALL QTRANS('AOMO','S',D0,
     &               NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &               WORK(KH1AO+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),
     &               FC(I2ORBX(I,I)+1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),
     &               WORK(KFREE),LFREE,IPROPT)
         IF (IPROPT .GE. 5) THEN
            WRITE(LUPRI,9011) 'FC',I
            CALL PRQMAT(FC(I2ORBX(I,I)+1),NORB(I),NORB(I),
     &                  NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
            call flshfo(lupri)
         END IF
  100 CONTINUE
C
C     Write to file KRMC FOCK
C
      write(6,*) 'writing files to KRMC FOCK'
      CALL WRTKRMC(LUKRM3,'FCMO    ',FC,N2ORBXQ)
            call flshfo(lupri)
c      call getmat('FCMO',FC,10,.TRUE.,WORK(KFREE),LFREE)
c      call prqmat(fc,norbt,norbt,norbt,norbt,nz,ipqtoq(1,0),lupri)
c      call getmat('FVMO',FV,10,.TRUE.,WORK(KFREE),LFREE)
c      call prqmat(fv,norbt,norbt,norbt,norbt,nz,ipqtoq(1,0),lupri)
C
C
C
      CALL MEMREL('RGETFC',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('RGETFC')
      RETURN
 9010 FORMAT(/'RGETFC: MO Fock matrices, irrep ',I1)
 9011 FORMAT(/'RGETFC: ',A,' fock matrix, irrep ',I1)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirh1 */
      SUBROUTINE DIRH1(H1AO,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Get one-electron Fock matrix in AO-basis
C
C     Input:
C
C     Output:
C        H1AO    - one-electron Fock matrix
C
C     Written by J. Thyssen - Oct 26 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
C
#include <dcbdhf.h>
      DIMENSION H1AO(*), WORK(*)
      LOGICAL   EX
C
      CALL QENTER('DIRH1')
C
C     (1) If file DFFCK1 is present, then read H1AO from DFFCK1.
C     ----------------------------------------------------------
C
      INQUIRE (FILE='DFFCK1', EXIST=EX)
      IF (EX) THEN
         CALL OPNFIL(LUFCK1,'DFFCK1','OLD','DIRH1')
         CALL REAFCK(LUFCK1,H1AO,.FALSE.,1)
         CLOSE(LUFCK1)
      ELSE
C
C     (2) Generate H1AO and save on DFFCK1.
C     -------------------------------------
C
         CALL ONEFCK(H1AO,IPRONE,WORK,LWORK)
         CALL OPNFIL(LUFCK1,'DFFCK1','UNKNOWN','ONEFC1')
         CALL WRIFCK(LUFCK1,H1AO,1)
         CLOSE (LUFCK1,STATUS='KEEP')
C
      END IF
C
      CALL QEXIT('DIRH1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsetwop */
      SUBROUTINE RSETWOP(WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Initialize configurational and orbital parameters for 
C     (KR-MC)SCF calculation.
C
C     Input:
C        OPT_SKIPEE  - skip e-e rotations
C        OPT_SKIPEP  - skip e-e rotations
C
C     Output:
C        KXOPE   - electron-electron rotations
C        NXOPE   - number of -"- 
C        KXOPP   - electron-positron rotations
C        NXOPP   - number of -"- 
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcbopt.h>
#include <dcborb.h>
C
      DIMENSION WORK(*)
C
C
      CALL QENTER('RSETWOP')
C
C     ***************************************************
C     *** Configurational parameters ********************
C     ***************************************************
C
      IF ( NASHT .GT. 1 ) THEN
C
C        A non-trivial number of CI rotations.
C         
         IF ( OPT_CIPROGRAM .EQ. 'GOSCI' ) THEN
C
C           Use GOSCI:
C
            CALL GOSCI_SETCB()
            CALL GOSCI_GENDET(NZCONF)
            NZCONFQ = NZCONF
            IF ( NZ .GE. 2 ) NZCONFQ = 2 * NZCONF
C
            WRITE(LUPRI,'(/A,I5,A,I5,A/)')
     $           ' * RSETWOP: Number of determinants: ',NZCONF,
     $           ' (dimension ',NZCONFQ,')'
C
         ELSE
C
            WRITE(LUPRI,'(A/3A)')
     $           '*** ERROR in RSETWOP ***',
     $           'Unknown CI program: "',OPT_CIPROGRAM,'"'
            CALL QUIT('*** ERROR in RSETWOP ***')
C
         END IF
      ELSE
C
C        Dirac-Hartree-Fock: only 1 CI coefficient
C     
         NZCONF = 1
         IF ( NZ .GE. 2 ) THEN
            NZCONFQ = 2
         ELSE
            NZCONFQ = 1
         END IF
C
      END IF
C
C     Allocate memory of CI rotations
C
      write(6,*) 'kfree before CI rot: ',kfree
      CALL MEMGET('INTE',KZCONF,2*NZCONFQ,WORK,KFREE,LFREE)
      write(6,*) 'kfree after CI rot: ',kfree
C
C     ***************************************************
C     *** Orbital parameters ****************************
C     ***************************************************
C
      WRITE(LUPRI,'(1X,A)') 'Orbitals classes:'
      DO I = 1,NFSYM
         WRITE(LUPRI,'(1X,A,I1)') '* Fermion irrep ',I
         IF (NISH(I).GT.0) THEN
            WRITE(OPT_INDSTR(1,I),'(I4,A2,I4)') 
     &         1,'..',NISH(I)
         ELSE
            OPT_INDSTR(1,I) = ''
         END IF
         IF (NASH(I).GT.0) THEN
            WRITE(OPT_INDSTR(2,I),'(I4,A2,I4)') 
     &         NISH(I)+1,'..',NISH(I)+NASH(I)
         ELSE
            OPT_INDSTR(2,I) = ''
         END IF
         IF (NPSH(I) .GT. 0 .AND. NSSH(I) .GT. 0) THEN
            WRITE(OPT_INDSTR(3,I),'(I4,A2,I4,A1,I4,A2,I4)') 
     &         -NPSH(I),'..',-1,',',NISH(I)+NASH(I)+1,'..',NESH(I)
         ELSE IF (NPSH(I) .GT. 0) THEN
            WRITE(OPT_INDSTR(3,I),'(I4,A2,I4)') -NPSH(I),'..',-1
         ELSE IF (NSSH(I) .GT. 0) THEN
            WRITE(OPT_INDSTR(3,I),'(I4,A2,I4)') 
     &         NISH(I)+NASH(I)+1,'..',NESH(I)
         ELSE
            OPT_INDSTR(3,I) = ''
         END IF
C
         WRITE(LUPRI,'(3X,A,A)')
     &      'Inactive orbitals   : ',OPT_INDSTR(1,I)
         WRITE(LUPRI,'(3X,A,A)')
     &      'Active orbitals     : ',OPT_INDSTR(2,I)
         WRITE(LUPRI,'(3X,A,A)')
     &      'Secondary orbitals  : ',OPT_INDSTR(3,I)
      END DO
      CALL FLSHFO(LUPRI)
C
C     Always totally symmetric rotations (=> JOPSY = 1)
C
      JOPSY = 1
      CALL SETXOP(OPT_INDSTR,JOPSY,NZXOPE,KZXOPE,NZXOPP,KZXOPP,
     &            WORK,KFREE,LFREE,OPT_SKIPEP,OPT_SKIPEE,.FALSE.,IPROPT)
      NZXOPEQ  = NZXOPE *NZ
      NZXOPPQ  = NZXOPP *NZ
      NZXOPT   = NZXOPE + NZXOPP
      NZXOPTQ  = NZXOPT *NZ
C
C     ***************************************************
C     *** Total number of parameters ********************
C     ***************************************************
C
      NZVAR    = NZCONF + NZXOPT
      NZVARQ   = NZVAR*NZ
C     
C
C
      CALL QEXIT('RSETWOP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgetac */
      SUBROUTINE RGETAC(FC,FCAC,IPRINT)
C*****************************************************************************
C
C     Extract block with active-active orbital indices out of FC
C
C     Input:
C        FC      - Inactive Fock matrix
C
C     Output:
C        FCAC    - Inactive Fock matrix with active-active indices
C
C     Written by J. Thyssen - Oct 29 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dcbgen.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION FC(*),FCAC(*)
C
      CALL QENTER('RGETAC')
      IF (IPRINT .GE. 5) CALL HEADER('Output from rGETAC',-1)
C
      IF (NASHT .EQ. 0) GOTO 9999
C
C     Print FC 
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('rGETAC: FC matrix',-1)
         CALL PRQMAT(FC,NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL DZERO(FCAC,N2ASHXQ)
      IOFF = 1
      DO 100 IFSYM = 1,NFSYM
         IF (NASH(IFSYM) .EQ. 0) GOTO 100
C
         ISTART = IORB(IFSYM) + NPSH(IFSYM)+NISH(IFSYM)+1
         IEND   = ISTART-1+NASH(IFSYM)
         CALL RGETA1(FC,NORBT,
     *               FCAC,NASHT,
     &               ISTART,IEND,IASH(IFSYM),NZ)
C
  100 CONTINUE
C
C     Print FCAC
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('rGETAC: FCAC matrix',-1)
         CALL PRQMAT(FCAC,NASHT,NASHT,NASHT,NASHT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Save on file LUKRM3
C
      CALL WRTKRMC(LUKRM3,'FCAC    ',FCAC,N2ASHXQ)
C
 9999 CALL QEXIT('RGETAC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgeta1 */
      SUBROUTINE RGETA1(FC,NORBT,FCAC,NASHT,ISTART,IEND,IASH,NZ)
C*****************************************************************************
C
C     Extract block with active-active orbital indices out of FC
C
C     Input:
C        FC      - Inactive Fock matrix for a given symmetry
C
C     Output:
C        FCAC    - Inactive Fock matrix with active-active indices
C                  for a given symmetry
C
C     Written by J. Thyssen - Oct 29 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
C
      DIMENSION FCAC(NASHT,NASHT,NZ),FC(NORBT,NORBT,NZ)
C
      CALL QENTER('RGETA1')
C
      DO IZ = 1,NZ
         K = IASH
         DO I = ISTART,IEND
            K = K + 1
            L = IASH
            DO J = ISTART,IEND
               L = L + 1
               FCAC(K,L,IZ) = FC(I,J,IZ)
            END DO
         END DO
      END DO
C
      CALL QEXIT('RGETA1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck renract */
      FUNCTION RENRACT(DV,FCAC,FQ,IPRINT)
C*****************************************************************************
C
C     Calculate active electronic energy of the KR-MCSCF state
C
C     Input:
C        DV      - one-electron active density matrix
C        FQ      - the Q Fock matrix
C        FCAC    - the inactive Fock matrix with active indices
C
C     Returns:
C        RENRACT - active electronic energy
C
C     Written by J. Thyssen - Oct 29 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (DP5 = 0.50D 00)
C
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION DV(*),FCAC(*),FQ(*)
C
      CALL QENTER('RENRACT')
C
C     Eactive = \sum_{uv} DV_{uv} FC_{uv} + 
C               \sum_{uvxy} P_{uv,xy} (uv|xy)
C             = \sum_{uv} DV_{uv} FC_{uv} + 
C               1/2 \sum_{u} FQ_{uu}
C
C
C
C     Print FCAC
C
      IF (IPRINT .GE. 5) THEN
         write(6,*) 'n2asht',N2ASHT
         CALL HEADER('rEnract: FCAC matrix',-1)
         CALL PRQMAT(FCAC,NASHT,NASHT,NASHT,NASHT,
     &               NZ,IPQTOQ(1,0),LUPRI)
         CALL HEADER('rEnract: DV matrix',-1)
         CALL PRQMAT(DV,NASHT,NASHT,NASHT,NASHT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF


      EACT = DDOT(N2ASHXQ,DV,1,FCAC,1)
      IF (IPRINT .GE. 4) 
     &   WRITE(LUPRI,9000) EACT
C
C
      IOFF = 1
      DO IFSYM = 1,NFSYM
         EACT = EACT + DP5 * RENRAC1(FQ(IOFF),NORBT,NASH(IFSYM),
     &                               NPSH(IFSYM)+NISH(IFSYM))
         IOFF = IOFF + NZ*NORBT*NASH(IFSYM)
      END DO
      IF (IPRINT .GE. 4) 
     &   WRITE(LUPRI,9001) EACT
      RENRACT = EACT
      CALL QEXIT('RENRACT')
      RETURN
 9000 FORMAT(/2X'rENRACT: sum_{uv} DV_{uv} FC_{uv} : ',F25.15)
 9001 FORMAT(/2X'rENRACT: Total active energy      : ',F25.15)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck renrac1 */
      FUNCTION RENRAC1(FQ,NORBT,INASH,IOFF)
C*****************************************************************************
C
C     Calculate \sum_{u} FQ_{u,u}
C
C     Input:
C        FQ      - the Q Fock matrix for a given symmetry
C
C     Returns:
C        RENRAC1 - \sum_{u} FQ_{u,u}
C
C     Written by J. Thyssen - Oct 29 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION FQ(NORBT,INASH)
C
      CALL QENTER('RENRAC1')
C
      DO I = 1,INASH
         E = E + FQ(I+IOFF,I)
      END DO
      RENRAC1 = E
C
      CALL QEXIT('RENRAC1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rreadmo */
      SUBROUTINE RREADMO(CMO,JRDMO,IOPT)
C*****************************************************************************
C
C     Retrieve MO coefficients from file
C
C     Input:
C
C     Output:
C        CMO     - MO coefficients
C
C     Written by J. Thyssen - Nov 3 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcborb.h>
#include <dcbopt.h>
#include <dgroup.h>
#include <dcbbas.h>
C
      DIMENSION CMO(*)
      LOGICAL NOSTOP
      CHARACTER LABEL*8
C
      CALL QENTER('RREADMO')
      IF (JRDMO .LT. 0) THEN
         JRDMO = 0
         NOSTOP = .TRUE.
      ELSE
         NOSTOP = .FALSE.
      END IF
      IF (IOPT .LT. 1 .OR. IOPT .GT. 2) THEN
          WRITE(LUPRI,'(/1X,A,I10)') 
     &       'rREADMO: ERROR: Unknown option ',IOPT
          CALL QUIT('rREADMO: ERROR: Unknown option.')
      END IF
      GOTO (100,100), IOPT
  100 CONTINUE
C
C     Read orbitals from label 'NEWORB'
C
      REWIND(LUKRMC)
      IERR = -1
      IF ( IOPT .EQ. 1 ) THEN
         LABEL = 'NEWORB  '
      ELSE
         LABEL = 'OLDORB  '
      END IF
C
      CALL MOLLAB(LABEL,LUKRMC,IERR)
      IF (IERR .GE. 0) THEN
         CALL READT(LUKRMC,NCMOTQ,CMO)
      ELSE IF (IERR .EQ. -1) THEN
         IF (NOSTOP) THEN
            JRDMO = JRDMO + 1
         ELSE
            CALL QUIT('rREADMO: ERROR: Label ' // LABEL // 
     &                'not found on file KRMCSCF')
         END IF
      ELSE
         IF (NOSTOP) THEN
            JRDMO = JRDMO + 1
         ELSE
            CALL QUIT('rREADMO: ERROR: error reading file "KRMCSCF"')
         END IF
      END IF
      IF ( JRDMO .GT. 0 ) THEN
         GO TO 1010
      ELSE
         GO TO 1000
      END IF
C
C     End:
C
 1000 CONTINUE
      IF (IPROPT .GE. 5) THEN
         CALL HEADER('Output from rREADMO:',-1)
         DO I = 1,NFSYM
            CALL HEADER('Coefficients for irrep '//CHAR(I+48),-1)
            CALL PRQMAT(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
 1010 CONTINUE
      CALL QEXIT('RREADMO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck roptst */
      SUBROUTINE ROPTST(WORK,LWORK)
C*****************************************************************************
C
C     Generate start vectors
C
C     Input: None
C
C     Output: None
C
C
C     Written by J. Thyssen - Nov 17 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcborb.h>
#include <dcbopt.h>
C
      LOGICAL EX, FNDLAB
      DIMENSION WORK(*)
C
#include <memint.h>
C
      CALL QENTER('ROPTST')
C
C     (1) MO Orbitals
C     ----------------
C
C     Get MO orbitals:
C
C     (1.1) From a previous MCSCF calculation.
C           Read orbitals from KRMCSCF file
C
      CALL MEMGET('REAL',KCMO  ,NCMOTQ        ,WORK,KFREE,LFREE)
      INQUIRE ( FILE = 'KRMCSCF', EXIST = EX )
      IF ( EX ) THEN
C
C     (1.1.1) Read from label NEWORB
C
         JRDMO = -1
         CALL RREADMO(WORK(KCMO),JRDMO,1)
         IF (JRDMO .EQ. 0) THEN
            WRITE(LUPRI,'(A)')
     $           '*** Start orbitals read from' //
     $           ' label NEWORB on file KRMCSCF'
            GO TO 1000
         END IF
C
C     (1.1.2) Read from label OLDORD
C     
         JRDMO = -1
         CALL RREADMO(WORK(KCMO),JRDMO,2)
         IF (JRDMO .EQ. 0) THEN
            WRITE(LUPRI,'(A)')
     $           '*** Start orbitals read from' //
     $           ' label OLDORB on file KRMCSCF'
            GO TO 1000
         END IF
      END IF
C
C     (1.2) From previous DHF calculation.
C           Read orbitals from DFCOEF
C
      INQUIRE ( FILE = 'DFCOEF', EXIST = EX )
      IF ( EX ) THEN
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','RREADMO')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,DUMMY,DUMMY,4)
         CLOSE(LUCOEF,STATUS='KEEP')
         WRITE(LUPRI,'(A)')
     $        '*** Start orbitals read from' //
     $        ' file DFCOEF'
         GO TO 1000
      END IF
C
C     (1.3) From a usual DHF start guess, i.e. bare nucleus guess.
C           FIXME: to be implemented
C 
C
C
 1000 CONTINUE
C
C     Write coefficients back to file as label OLDORB
C
      CALL WRTKRMC(LUKRMC,'OLDORB  ',WORK(KCMO),NCMOTQ)
C
C
C     (2) CI Coefficients
C     -------------------
C
C
      IF ( NZCONF .GT. 1 ) THEN
C
C        Do a 0th order 4-index transformation.
C
         JTRLVL = 0
         CALL RTRACTL(JTRLVL,WORK(KFREE),LFREE)
C
         CALL MEMGET('REAL',KCREF ,NZCONFQ     ,WORK,KFREE,LFREE)
         CALL RCIST(WORK(KCMO),WORK(KCREF),WORK(KFREE),LFREE)
C     
C        Write CI vector to file as 'CISTART0'
C     
         CALL WRTKRMC(LUKRMC,'CISTART0',WORK(KCREF),NZCONFQ)
C
      END IF
C     
C     Write coefficients back to file as label NEWORB
C
      CALL WRTKRMC(LUKRMC,'NEWORB  ',WORK(KCMO),NCMOTQ)
C      
      CALL MEMREL('ROPTST',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('ROPTST')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rcist */
      SUBROUTINE RCIST(CMO,CREF,WORK,LWORK)
C*****************************************************************************
C
C     Generate CI start vectors.
C
C     Input : CMO - MO coefficents
C
C     Output: CREF - CI vector
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcborb.h>
C
      CALL QENTER('RCIST')
C
C     Do some (e.g., 3) Davidson CI iterations to get a start vector:
C
C     (2.1) Contruct diagonal of CI matrix
C     (2.2) Pick start vector with zeros everywhere except for the lowest 
C           value of the diagonal of the CI matrix. If 2 or more diagonal
C           value are close, pick more start vectors.
C     (2.3) Do a couple of Davidson CI iterations.
C
C     ...or (the quick hack implemented here):
C
C     Diagonalize the full CI matrix now. This is of course only possible
C     for sufficiently small CI's.
C
      CALL GOSCI_STARTVEC(CMO,CREF,WORK,LWORK)
C
      CALL QEXIT('RCIST')
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirnr */
      SUBROUTINE DIRNR(CREF,GCI,GOEE,GOEP,CMO,DV,FC,FV,FCAC,
     &                 H2AC,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Perform NR step
C
C     Input:
C        CMO     - the orbitals for the CEP
C        CREF    - CI vector
C        G       - the gradient
C
C     Output:
C
C     Written by J. Thyssen and H. J. Aa. Jensen - Nov 17 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
      PARAMETER (D0 = 0.00D00, D1 = 1.00D00)
C
#include <dcbopt.h>
#include <dgroup.h>
#include <dcborb.h>
#include <dcbkrmc_itinfo.h>
C
      LOGICAL OPTRST,OPTNRM,OPTDIH,OPTSTERN
      DIMENSION CREF(*),GCI(*),DV(*),FC(*),FV(*),
     &          H2AC(*),WORK(*),GOEE(*),GOEP(*)
      DIMENSION ITRINT_MIC(2),CNVINT_MIC(2)
C
      CALL QENTER('DIRNR')
      KFRSAV = KFREE
C
C
C     ***************************************************
C     *** Transfer info to Response common blocks *******
C     ***************************************************
C
C
      OPTRST = .FALSE.
      OPTNRM = .FALSE.
CDEBUG
      OPTNRM = .TRUE.
CEND DEBUG
      OPTDIH = .TRUE.
      OPTSTERN = .FALSE.
C
C     Calculate with the same integrals used for gradient and energy.
C     In a refined version one might skip some of the integrals in
C     the micro iterations.
C
      INTDEF_MIC = IOPT_INTDEF
      ITRINT_MIC(1) = 1
      ITRINT_MIC(2) = 1
      CNVINT_MIC(1) = DUMMY
      CNVINT_MIC(2) = DUMMY
C
C     Heuristics ripped from Dalton (CNVTHR, sirneo.F)
C
      THR = MIN (0.2 * GNORM(5), 0.8 * GNORM(5)**2 )
      IF ( THR .LE. OPT_THRGRD ) THEN
         THR = 0.98D00 * OPT_THRGRD
      ELSE IF ( THR .LE. 4.0D00 * OPT_THRGRD ) THEN
         THR = 0.7D00 * OPT_THRGRD
      END IF
      write(6,*) 'thr in microit',thr
      CALL TDCBXRS(OPTRST,OPTNRM,OPTDIH,OPTSTERN,IPROPT,MOPT_MXMICRO,
     &             INTDEF_MIC,ITRINT_MIC,CNVINT_MIC,
     &             NZVAR, NZVARQ ,NZCONF,NZCONFQ,
     &             NZXOPT,NZXOPTQ,NZXOPE,NZXOPEQ,NZXOPP,NZXOPPQ,
     &             KZCONF,KZXOPE,KZXOPP,THR)
                   
C
C
C
      NREDM = 50
      NEXSIM = 1
      CALL MEMGET('REAL',KIBTYP,2*NREDM     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KIBCVC,NREDM       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KIBEVC,NREDM       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KIBPVC,NREDM       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVALR,NEXSIM      ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVECR,NEXSIM*NREDM,WORK,KFREE,LFREE)
C
C
C     ***************************************************
C     *** Perform micro iterations **********************
C     ***************************************************
C
C
      WORK(KEVALR) = D0
      CALL XRSCTL(GCI,GOEE,GOEP,
     &            WORK(KIBTYP),WORK(KIBCVC),WORK(KIBEVC),WORK(KIBPVC),
     &            RCNV,WORK(KEVALR),WORK(KEVECR),WORK,KFREE,LFREE)
C
C
      CALL GTNORB(CMO,WORK(KIBTYP),
     &            WORK(KIBCVC),WORK(KIBEVC),WORK(KIBPVC),
     &            WORK(KEVALR),WORK(KEVECR),NREDM,WORK,KFREE,LFREE)
C
C
C     Save information for final summary output
C
      DINFO(7) = D0
      DINFO(8) = STPLEN
      DINFO(10) = D1
C
      CALL MEMREL('DIRNR.after XRSCTL',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DIRNR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tdcbxlr */
      SUBROUTINE TDCBXRS(OPTRST,OPTNRM,OPTDIH,OPTSTERN,IPROPT,MXMICRO,
     &                   INTDEF_MIC,ITRINT_MIC,CNVINT_MIC,
     &                   NZVAR_OPT,NZVARQ_OPT,NZCONF_OPT,NZCONFQ_OPT,
     &                   NZXOPT_OPT,NZXOPTQ_OPT,NZXOPE_OPT,NZXOPEQ_OPT,
     &                   NZXOPP_OPT,NZXOPPQ_OPT,
     &                   KZCONF_OPT,KZXOPE_OPT,KZXOPP_OPT,THR)
C*****************************************************************************
C
C     Transfer information to response common blocks (dcbxnr.h)
C
C     Input:
C        Information from optimization module
C
C     Output:
C
C     Written by J. Thyssen - Nov 19 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
C
      LOGICAL OPTRST,OPTNRM,OPTDIH,OPTSTERN
      DIMENSION ITRINT_MIC(2),CNVINT_MIC(2)
C
C     Logicals...
C
      RSREST = OPTRST
      STATIC = .TRUE.
      LINEQ  = .TRUE.
      LSVCFG(1) = .TRUE.
      LSVCFG(2) = .TRUE.
      TKNORM = OPTNRM
      DIAGHE = OPTDIH
      STERNH = OPTSTERN
      FKRMC  = .TRUE.
C
C     integers...
C
      IPRXRS = IPROPT
      JSYMOP = 1
      JTIMOP = 1
      JOPSY  = 1
      MAXITR = MXMICRO
      NFREQ  = 1
      NEXSIM = NFREQ
      NEXSTV = NFREQ
      NEXCNV = NFREQ
      LOFFTY = 0
      NREDM  = 50
      NCRED  = 0
      NERED  = 0
      NPRED  = 0
      NZRED  = NCRED + NERED + NPRED
      N2REDM  = NREDM**2
      INTDEF = INTDEF_MIC
      ITRINT(1) = ITRINT_MIC(1)
      ITRINT(2) = ITRINT_MIC(2)
      NZVAR     = NZVAR_OPT
      NZVARQ    = NZVARQ_OPT
      NZCONF    = NZCONF_OPT
      NZCONFQ   = NZCONFQ_OPT
      NZXOPT    = NZXOPT_OPT
      NZXOPTQ   = NZXOPTQ_OPT
      NZXOPE    = NZXOPE_OPT
      NZXOPEQ   = NZXOPEQ_OPT
      NZXOPP    = NZXOPP_OPT
      NZXOPPQ   = NZXOPPQ_OPT
C
C     reals...
C
      THCXRS = THR
      CNVINT(1) = CNVINT_MIC(1)
      CNVINT(1) = CNVINT_MIC(2)
      RESFAC    = 1.0D3
C
C     pointers...
C
      KZCONF = KZCONF_OPT
      KZXOPE = KZXOPE_OPT
      KZXOPP = KZXOPP_OPT
C
C     units...
C
      LUXVC = 29
      LURST = 29
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtnorb */
      SUBROUTINE GTNORB(CMO,IBTYP,IBCVC,IBEVC,IBPVC,
     &                  EVALR,EVECR,NREDM,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Calculate new orbitals ( CMO_{new} = exp ( -\kappa) CMO_{old} )
C
C     Input:
C        CMO - old coefficients
C
C     Output:
C        CMO - new coefficients
C
C     Written by J. Thyssen - Nov 24 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <dcbopt.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbgen.h>
C
      DIMENSION CMO(*), WORK(*), IBTYP(2,*), IBCVC(*),
     &          IBEVC(*),IBPVC(*), EVECR(NREDM,*)
C
      CALL QENTER('GTNORB')
      KFRSAV = KFREE
C
C
C     ***************************************************
C     *** Get solution vectors = kappa ******************
C     ***************************************************
C
C
      CALL MEMGET('REAL',KKAPPA,NZXOPTQ,WORK,KFREE,LFREE)
      CALL KRMCVEC(IBTYP,IBCVC,IBEVC,IBPVC,
     &           EVALR,EVECR,WORK(KKAPPA),WORK(KFREE),LFREE)
C
      IF (IPROPT .GE. 5) THEN
         CALL HEADER('GTNORB: Kappa:',-1)
         CALL RPRKAP(WORK(KKAPPA),WORK(KKAPPA+NZXOPEQ),
     &               WORK(KZXOPE),NZXOPE,
     &               WORK(KZXOPP),NZXOPP,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C
C     ***************************************************
C     *** Rotate old coefficients ***********************
C     ***************************************************
C
C     Unpack KAPPA
C
      CALL MEMGET('REAL',KUKAPPA,N2ORBXQ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KUKAPPA),N2ORBXQ)
      CALL RUPKWOP(WORK(KKAPPA),WORK(KKAPPA+NZXOPEQ),
     &             WORK(KUKAPPA),NORBT,
     &             WORK(KZXOPE),NZXOPE,
     &             WORK(KZXOPP),NZXOPP,NZ)
      IF (IPROPT .GE. 5) THEN
         CALL HEADER('GTNORB: Unpacked Kappa:',-1)
         CALL PRQMAT(WORK(KUKAPPA),NORBT,NORBT,NORBT,NORBT,NZ,
     &               IPQTOQ(1,0),LUPRI)
      END IF
C
C     Save old coefficients
C
      CALL WRTKRMC(LUKRM2,'OLDORB  ',CMO,NCMOTQ)
C
C     Get new coefficients
C
      IOFF = 0
      DO I = 1,NFSYM
         IF (NOCC(I) .GT. 0) THEN
CJTH     IF (NORB(I) .GT. 0) THEN
            CALL HEADER('Fermion symmetry '//CHAR(I+48),-1)
            CALL ROTCMO(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,
     &                  WORK(KUKAPPA+IOFF),
     &                  NORBT,IPQTOQ(1,0),IPROPT,WORK(KFREE),LFREE)
         END IF
         IOFF = IOFF + (NORBT+1)*NORB(I)
      END DO
C
C     Write new coefficients
C
      CALL WRTKRMC(LUKRM2,'NEWORB  ',CMO,NCMOTQ)
C
      CALL MEMREL('GTNORB',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('GTNORB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rprkap */
      SUBROUTINE RPRKAP(XKAPEE,XKAPEP,JXOPE,NZXOPE,JXOPP,NZXOPP,
     &                  NZ,IP,LUPRI)
C*****************************************************************************
C
C     Print Kappa vector
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Nov 24 1998
C
C*****************************************************************************
#include <implicit.h>
C
      DIMENSION XKAPEE(NZXOPE,NZ),XKAPEP(NZXOPP,NZ),
     &          JXOPE(2,NZXOPE),JXOPP(2,NZXOPP),IP(4)
      CHARACTER M(4)*1
      DATA M /'1','i','j','k'/
C
      IF (NZXOPE .EQ. 0) GOTO 100
      WRITE(LUPRI,9000) 'e-e'
      WRITE(LUPRI,9005) (M(IP(I)),I=1,NZ)
      WRITE(LUPRI,9006) ('---------------  ',I=1,NZ)
      DO IG = 1,NZXOPE
         WRITE(LUPRI,9010) IG,JXOPE(1,IG),JXOPE(2,IG),
     &                     (XKAPEE(IG,I),I=1,NZ)
      END DO
  100 CONTINUE
      IF (NZXOPP .EQ. 0) GOTO 200
      WRITE(LUPRI,9000) 'e-p'
      WRITE(LUPRI,9005) (M(I),I=1,NZ)
      WRITE(LUPRI,9006) ('---------------  ',I=1,NZ)
      DO IG = 1,NZXOPP
         WRITE(LUPRI,9010) IG,JXOPP(1,IG),JXOPP(2,IG),
     &                     (XKAPEP(IG,I),I=1,NZ)
      END DO
 9000 FORMAT(/A,' part of matrix:')
 9005 FORMAT(1P,/'Index(r,s)   r    s   ',4(A1,'-part           '))
 9006 FORMAT(    '----------  ---  ---  ',4A17)
 9010 FORMAT(1P,4X,I4,4X,I3,2X,I3,2X,4(D15.8,'  '))
C
  200 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrtkrmc */
      SUBROUTINE WRTKRMC(LUNIT,LABEL,VECTOR,ISIZE)
C*****************************************************************************
C
C     Write VECTOR to LUKRMC under label LABEL.
C
C     Input:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Output:
C
C     Written by J. Thyssen - Nov 24 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <iratdef.h>
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8,C(4)*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C     Label exists.
C     Seek to label SODLABEL and write from there.
C
         REWIND(LUNIT)
         IF (.NOT. FNDLAB('SODLABEL',LUNIT)) THEN
            CALL QUIT('WRTKRMC: wrong structure in file')
         END IF
      ELSE
C
C     Label do not exist
C     Seek to label EOFLABEL and write from there
C
         REWIND(LUNIT)
         IF (FNDLAB('EOFLABEL',LUNIT))
     &      BACKSPACE LUNIT
C
      END IF
C
C     Write label
C
      CALL NEWLAB(LABEL,LUNIT,LUERR)
C
C     Write vector
C
      LEN = MAX(4,ISIZE)
      CALL WRITSQ(LUNIT,IRAT*LEN,VECTOR)
      CALL NEWLAB('EOFLABEL',LUNIT,LUERR)
C
      REWIND LUNIT
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reakrmc */
      SUBROUTINE REAKRMC(LUNIT,LABEL,VECTOR,ISIZE)
C*****************************************************************************
C
C     Read VECTOR from LUKRMC under label LABEL.
C
C     Input:
C
C     Output:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Written by J. Thyssen - Jul 4 2000
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <iratdef.h>
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8,C(4)*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C        Label exists.
C
         LEN = MAX(4,ISIZE)
         CALL READT(LUNIT,ISIZE,VECTOR)
C
      ELSE
         WRITE(LUPRI,'(//A/,3A,I3)')
     $        '*** ERROR in REAKRMC ***',
     $        'Label <',LABEL,'> does not exists on unit ',LUNIT
         CALL QUIT('*** ERROR in REAKRMC ***')
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rupkwop */
      SUBROUTINE RUPKWOP(WOPEE,WOPEP,UNWOP,NORBT,
     &                   JXOPE,NZXOPE,JXOPP,NZXOPP,NZ)
C*****************************************************************************
C
C     Unpack WOP into UNWOP
C     WOP is a rotation type vector and UNWOP is a square matrix.
C
C     Input:
C        WOPEE   - e-e orbital rotation type vector
C        WOPEP   - e-p orbital rotation type vector
C        JXOPE   - e-e rotations
C        JXOPP   - e-p rotations
C
C     Output:
C        UNWOP   - unpacked matrix
C
C     Written by J. Thyssen - Nov 24 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION WOPEE(NZXOPE,NZ),WOPEP(NZXOPP,NZ),
     &          UNWOP(NORBT,NORBT,NZ),JXOPE(2,*),JXOPP(2,*)
C
C     e-e rotations
C
      DO IG = 1,NZXOPE
         K = JXOPE(1,IG)
         L = JXOPE(2,IG)
         DO IZ = 1,NZ
            UNWOP(K,L,IZ) =  WOPEE(IG,IZ)
            UNWOP(L,K,IZ) = -WOPEE(IG,IZ)
         END DO
      END DO
C
C     e-p rotations
C
      DO IG = 1,NZXOPP
         K = JXOPP(1,IG)
         L = JXOPP(2,IG)
         DO IZ = 1,NZ
            UNWOP(K,L,IZ) =  WOPEP(IG,IZ)
            UNWOP(L,K,IZ) = -WOPEP(IG,IZ)
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rotcmo */
      SUBROUTINE ROTCMO(CMO,NBASI,NORBI,NZ,XKAPPA,NORBT,IP,IPRINT,
     &                  WORK,LWORK)
C*****************************************************************************
C
C     Rotate coefficients
C
C     CMO_{new} = CMO_{old} exp( -\kappa ) 
C     We use a trick a la Sostrup notes, chapter 3.1.5
C
C     Input:
C        CMO     - old coefficients
C        XKAPPA  - Kappa matrix (Kappa is overwritten with junk)
C        IPRINT  - print level 
C
C     Output:
C        CMO     - new coefficients
C
C     Written by J. Thyssen - Nov 24 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D1 = 1.00 D00, D0 = 0.00 D00, DM1 = -D1)
      PARAMETER (DTHRS = 1.0D-8, DSIXTH = D1/6.00D00, DP5 = 0.50D00)
      PARAMETER (D1P5 = 1.5D00, D2 = 2.0D00, THREQL = 1.0D-12)
C
      DIMENSION CMO(NBASI,NORBI,NZ), XKAPPA(NORBT,NORBT,NZ), 
     &          IP(4), WORK(*)
      DIMENSION IPK(4)
      CHARACTER SECTID*12,CPUTID*12,WALLTID*12
C
#include <memint.h>
C
      CALL QENTER('ROTCMO')
      CALL GETTIM(CPU1,WALL1)
      IF (IPRINT .GE. 5) CALL HEADER('Output from ROTCMO',-1)
C
#ifdef OLDCODE
C
      NORBXQ = NORBI*NORBI*NZ
      CALL MEMGET('REAL',KEIG,NORBI,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KKSQ,NORBI*NBASI*NZ,WORK,KFREE,LFREE)
      KCOS = KKSQ
      CALL MEMGET('REAL',KSIN,NORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KV  ,NORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP,NORBXQ,WORK,KFREE,LFREE)
C
C     Calculate -\kappa^2 = \kappa^{\dagger} \kappa
C     ---------------------------------------------
C
      CALL QGEMM(NORBI,NORBI,NORBI,D1,
     &           'H','N',IP,XKAPPA,NORBT,NORBT,NZ,
     &           'N','N',IP,XKAPPA,NORBT,NORBT,NZ,
     &           D0,IP,WORK(KKSQ),NORBI,NORBI,NZ)
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: -kappa^2',-1)
         CALL PRQMAT(WORK(KKSQ),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Diagonalize -\kappa^2
C     ---------------------------------------------
C
      CALL QDIAG(NZ,NORBI,WORK(KKSQ),NORBI,NORBI,
     &           WORK(KEIG),1,WORK(KV),NORBI,NORBI,
     &           WORK(KFREE),LFREE,IERR)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: V matrix (eigenvector)',-1)
         CALL PRQMAT(WORK(KV),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO: eigenvalues',-1)
         CALL PRQMAT(WORK(keig),1,NORBI,1,NORBI,1,IP,LUPRI)
      END IF
C
C     Make cos(d) and sin(d)/d ; d = sqrt(eigenvalues)
C     ------------------------------------------------
C
      CALL GETTIM(CPU2,WALL2)
      CALL DZERO(WORK(KCOS),NORBXQ)
      CALL DZERO(WORK(KSIN),NORBXQ)
      DO I = 1,NORBI
         IADDR = (I-1) + (I-1)*NORBI
         EIG = WORK(KEIG+I-1)
         IF (EIG .LT. D0) THEN
            IF (IPRINT .GE. 3)
     &         WRITE(LUPRI,9010) EIG
            EIG = D0
            WORK(KEIG+I-1)   = D0
         END IF
         IF (EIG .LT. DTHRS) THEN
C
C           Use Taylor expansion to avoid division by zero and/or
C           to speed up evaluation (avoid call to sin/cos/sqrt).
C           For EIG < 1.0D-8 the 3. order Taylor expansion is correct
C             to machine precision.
C
            WORK(KCOS+IADDR) = D1-DP5*EIG
            WORK(KSIN+IADDR) = D1-DSIXTH*EIG
         ELSE
            SQEIG = SQRT(EIG)
            WORK(KCOS+IADDR) = COS(SQEIG)
            WORK(KSIN+IADDR) = SIN(SQEIG)/SQEIG
         END IF
      END DO
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: sin(d)/d',-1)
         CALL PRQMAT(WORK(KSIN),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO: cos(d)',-1)
         CALL PRQMAT(WORK(KCOS),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU3,WALL3)
 9010 FORMAT('WARNING from ROTCMO: Negative eigenvalue ',1P,D15.6)
C
C     Make V cos(d) V^{\dagger}
C     (use KTMP as temporatory storage)
C     ------------------------------------
C
      CALL QTRANS('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI,
     &            WORK(KTMP),NORBI,NORBI,NZ,IP,
     &            WORK(KCOS),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KFREE),LFREE,IPROPT)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: V cos(d) V^{dagger}',-1)
         CALL PRQMAT(WORK(KTMP),NORBI,NORBI,NORBI,NORBI,NZK,IPK,LUPRI)
      END IF
C
C     Make V sin(d)/d V^{\dagger}
C     (use KCOS as temporatory storage)
C     ------------------------------------
C
      CALL QTRANS('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI,
     &            WORK(KCOS),NORBI,NORBI,NZ,IP,
     &            WORK(KSIN),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KFREE),LFREE,IPROPT)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: V sin(d)/d V^{dagger}',-1)
         CALL PRQMAT(WORK(KCOS),NORBI,NORBI,NORBI,NORBI,NZK,IPK,LUPRI)
      END IF
C
C     Make V cos(d) V^{\dagger} - V sin(d)/d V^{\dagger} \kappa
C     ---------------------------------------------------------
C
      CALL QGEMM(NORBI,NORBI,NORBI,DM1,
     &           'N','N',IP,WORK(KCOS),NORBI,NORBI,NZ,
     &           'N','N',IP,XKAPPA,NORBT,NORBT,NZ,
     &           D1,IP,WORK(KTMP),NORBI,NORBI,NZ)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: exp(-kappa)',-1)
         CALL PRQMAT(WORK(KTMP),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Check if unitary
C
C#ifdef KRMC_DEBUG
      CALL QGEMM(NORBI,NORBI,NORBI,D1,
     &           'N','N',IP,WORK(KTMP),NORBI,NORBI,NZ,
     &           'H','N',IP,WORK(KTMP),NORBI,NORBI,NZ,
     &           D0,IP,WORK(KCOS),NORBI,NORBI,NZ)
      ioff = 0
      i_max = 0
      j_max = 0
      d_max = D0
      do i = 1,norbi
      do j = 1,norbi
         if (i.eq.j) then
            if (work(kcos+ioff)-d1.gt.d_max) then
               i_max = i
               j_max = j
               d_max = work(kcos+ioff)-d1
            end if
         else
            if (work(kcos+ioff).gt.d_max) then
               i_max = i
               j_max = j
               d_max = work(kcos+ioff)
            end if
         end if
         ioff = ioff + 1
      end do
      end do
      if (d_max .gt. d0) then
         write(6,*) 'ROTCMO: maximum deviation from unitary matrix'
         write(6,*) 'Element ',i_max,j_max,' value :',d_max
      end if
C#endif
C
C
C     Calculate new orbitals: CMO_{old} exp(-kappa)
C     ---------------------------------------------
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: old orbitals',-1)
         CALL PRQMAT(CMO,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL QGEMM(NBASI,NORBI,NORBI,D1,
     &           'N','N',IP,CMO,NBASI,NORBI,NZ,
     &           'N','N',IP,WORK(KTMP),NORBI,NORBI,NZ,
     &           D0,IP,WORK(KCOS),NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,WORK(KCOS),1,CMO,1)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('ROTCMO: new orbitals',-1)
         CALL PRQMAT(cmo,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU4,WALL4)
      IF (IPRINT .GE. 2) THEN
         CPUTID = SECTID(CPU4-CPU1)
         WALLTID = SECTID(WALL4-WALL1)
         WRITE(LUPRI,9020) CPUTID,WALLTID
         CPUTID = SECTID(CPU3-CPU2)
         WALLTID = SECTID(WALL3-WALL2)
         WRITE(LUPRI,9021) CPUTID,WALLTID 
      END IF
 9020 FORMAT(/'ROTCMO: total CPU (WALL) time       : ',A12,' (',A12,')')
 9021 FORMAT( 'ROTCMO: sin/sqrt/cos CPU (WALL) time: ',A12,' (',A12,')')
C
#else
C
C     Memory alloc
C
      CALL MEMGET('REAL',KXKOLD,NORBI*NORBI*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXKNEW,NORBI*NORBI*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP  ,NBASI*NBASI*NZ,WORK,KFREE,LFREE)
C
C
C     Calculate exp(-K) = 1 - K
C
      K = 0
      DO IZ = 1,NZ
         DO J = 1,NORBI
            DO I = 1,NORBI
               K = K + 1
               WORK(KXKNEW + K - 1) = -XKAPPA(I,J,IZ)
               IF ( I .EQ. J .AND. IZ .EQ. 1)
     $              WORK(KXKNEW + K - 1 ) = WORK(KXKNEW + K - 1 ) +  D1
            END DO
         END DO
      END DO
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('ROTCMO: 1 - Kappa',-1)
         CALL PRQMAT(WORK(KXKNEW),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
      DMX = D0
      ITS = 0
C
C     Symmetric Orthogonalization with N-R iterative algorithm:
C
C     C(p+1) = 1/2( 3C(p) - C(p) S Ct(p) C(p))
C     
C            = 3/2 C(p) + ( -1/2 C(p) Ct(p) ) C(p)
C
C     (here S = unit matrix)
C
 100  CONTINUE
         ITS = ITS + 1
         CALL DCOPY(NORBI*NORBI*NZ,WORK(KXKNEW),1,WORK(KXKOLD),1)
         CALL QGEMM(NORBI,NORBI,NORBI,D1,
     &        'N','N',IP,WORK(KXKOLD),NORBI,NORBI,NZ,
     &        'T','N',IP,WORK(KXKOLD),NORBI,NORBI,NZ,
     &        D0,     IP,WORK(KTMP),  NORBI,NORBI,NZ)
         CALL QGEMM(NORBI,NORBI,NORBI,-DP5,
     &        'N','N',IP,WORK(KTMP),  NORBI,NORBI,NZ,
     &        'N','N',IP,WORK(KXKOLD),NORBI,NORBI,NZ,
     &        D1P5,   IP,WORK(KXKNEW),NORBI,NORBI,NZ)
C
         DMX1 = DMX
         DMX = DABSMAXDIFF(WORK(KXKNEW),WORK(KXKOLD),NORBI*NORBI*NZ)
         IF (DMX .GT. D2 * DMX1 .AND. DMX1 .GT. D0)
     $        CALL QUIT('divergence...')
         IF (DMX .GT. THREQL) GOTO 100
C
C                  
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: symm. ortho. rot. mat.',-1)
         CALL PRQMAT(WORK(KXKNEW),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C
C     Calculate new orbitals: CMO_{old} * (symm.ortho.rot.mat)
C     --------------------------------------------------------
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: old orbitals',-1)
         CALL PRQMAT(CMO,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL QGEMM(NBASI,NORBI,NORBI,D1,
     &           'N','N',IP,CMO,NBASI,NORBI,NZ,
     &           'N','N',IP,WORK(KXKNEW),NORBI,NORBI,NZ,
     &           D0,IP,WORK(KTMP),NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,WORK(KTMP),1,CMO,1)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('ROTCMO: new orbitals',-1)
         CALL PRQMAT(cmo,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
#endif
      CALL MEMREL('ROTCMO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('ROTCMO')
      RETURN
      END
      DOUBLE PRECISION FUNCTION DABSMAXDIFF(A,B,N)
#include <implicit.h>
      PARAMETER ( D0 = 0.0D00 )
      DOUBLE PRECISION A(N),B(N)
      DMX = D0
      DO I = 1,N
         DEL = ABS(A(N) - B(N))
         DMX = MAX(DMX,DEL)
      END DO
      DABSMAXDIFF = DMX
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirstp */
      SUBROUTINE DIRSTP(ISTEP,CMO,WORK,LWORK)
C*****************************************************************************
C
C     Step control
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Nov 27 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (THDE = 1.0D-10, D1 = 1.00D00)
C
#include <dcbopt.h>
#include <dcbkrmc_itinfo.h>
C
      DIMENSION CMO(*), WORK(*)
      LOGICAL   NOREJT
C
      CALL QENTER('DIRSTP')
C
      NOREJT = .FALSE.
      IF (ISTEP .LT. 0) THEN
         NOREJT = .TRUE.
      END IF
C
      ISTEP = 0
      DEACT = EMCSCF - EMCOLD
      IF ( ABS(DEPRED) .GT. THDE ) THEN
         RATIO = DEACT / DEPRED
         WRITE(LUPRI,9000) DEACT,DEPRED,RATIO
      ELSE
         RATIO = D1
         ISTEP = -1
         WRITE(LUPRI,9002) DEACT,DEPRED
      END IF
 9000 FORMAT(/'DIRSTP: Energy difference;',
     &       /' actual, predicted, and ratio: ',1P,2D15.5,0P,F12.6)
 9002 FORMAT(/'DIRSTP: Close to convergence, ratio set to one.'
     &       /' Energy difference; actual and predicted: ',1P,2D15.5)
C
      IF (ISTEP .LE. 0) THEN
         EMCOLD = EMCSCF
         GOTO 1000
      END IF
C
 1000 CONTINUE
C
C     Save information for final summary outout
C
      DINFO(4) = DEPRED
      DINFO(5) = DEACT
      DINFO(6) = RATIO
      CALL QEXIT('DIRSTP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmc_getfck */
      SUBROUTINE KRMC_GETFCK(FMO,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Read Fock matrix from file KRMC FOCK
C
C     Input:
C
C     Output:
C        FMO    - Fock matrix in MO basis
C
C     Written by J. Thyssen - Nov 27 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
#include <dcbgen.h>
C
      DIMENSION WORK(*),FMO(*)
C
#include <memint.h>
C
      REWIND(LUKRM3)
      CALL MOLLAB('FC      ',LUKRM3,LUERR)
      CALL READT(LUKRM3,N2ORBXQ,FMO)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Output from KRMC_GETFCK',-1)
         CALL HEADER('KRMC_GETFCK: FMO',-1)
         CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C   
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgeth2 */
      SUBROUTINE RGETH2(FQ,WORK,LWORK)
C*****************************************************************************
C
C     Calculate FQ and get H2AC (two-electron integrals
C     with active indeces).
C
C     Input:
C
C     Output:
C        FQ   -   aux. FQ fock matrix
C        H2AC -   (written to disk)
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbopt.h>
#include <dcborb.h>
#include <dcbgen.h>
C
      DIMENSION FQ(*)
      DIMENSION WORK(*)
C
C
#include <memint.h>
C
      CALL QENTER('RGETH2')
C
      IF ( OPT_CIPROGRAM .EQ. 'GOSCI' ) THEN
         write(6,*) 'rgeth2: lenpv',lenpv,(2*nasht)**4
         CALL MEMGET('REAL',KH2AC,LENPV,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KPV,LENPV,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KIDXAC,8*NORBT,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KIKR,4*NASHT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KFQ,2*NASHT*2*NORBT,WORK,KFREE,LFREE)
         CALL REAKRMC(LUKRM3,'PVMO    ',WORK(KPV),LENPV)
         CALL GOSCI_RGETH2(WORK(KFQ),2*NORBT,WORK(KH2AC),
     $        WORK(KIDXAC),WORK(KIKR),WORK(KPV))
         CALL MEMCHK('after gosci_rgeth2',WORK,KWORK)
         CALL MFQ2QFQ(WORK(KFQ),FQ,IPROPT)
         CALL MEMCHK('after mfq2qfq',WORK,KWORK)
         CALL WRTKRMC(LUKRM3,'FQMO    ',FQ,N2ASHXQ)
         CALL WRTKRMC(LUKRM3,'H2ACMO  ',WORK(H2AC),LENPV)
         CALL MEMREL('RGETH2.GOSCI',WORK,KWORK,KWORK,KFREE,LFREE)
         call quit('ralle ralle')
      ELSE
         WRITE(LUPRI,'(/,A,/,2A)')
     $        '*** ERROR in RGETH2 ***',
     $        'FQ and H2AC not implemented for CI program ',
     $        OPT_CIPROGRAM
         CALL QUIT('*** ERROR in RGETH2 ***')
      END IF
      CALL QEXIT('RGETH2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsnrprd */
      SUBROUTINE RSNRPRD(REDGP,REDE,EVECR,NREDM)
C*****************************************************************************
C
C     Calculate second order change in energy
C
C     Written by J. Thyssen - Nov 27 1998
C     Last revision : 
C       
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C         
      PARAMETER (DP5 = 0.50D00, D0 = 0.00D00, DP25 = 0.25D00)
C         
#include <dcbopt.h>
C           
      DIMENSION REDGP(NREDM),REDE(NREDM,NREDM),EVECR(NREDM)
C         
      CALL HEADER('Output from RSNRPRD',-1)
      CALL HEADER('solution vector',-1)
      CALL OUTPUT(EVECR,1,NREDM,1,1,NREDM,1,1,LUPRI)
      CALL HEADER('gradient',-1)
      CALL OUTPUT(REDGP,1,NREDM,1,1,NREDM,1,1,LUPRI)
      CALL HEADER('hessian',-1)
      CALL OUTPUT(REDE,1,NREDM,1,NREDM,NREDM,NREDM,1,LUPRI)
      SUM = D0
      DO I = 1,NREDM
         SUM1 = D0
         DO J = 1,NREDM
            SUM1 = SUM1 + REDE(I,J)*EVECR(J)
         END DO
         SUM = SUM + EVECR(I)*(-REDGP(I) + DP5 * SUM1)
      END DO
      write(6,*) 'rsnrprd: predicted energy change: ',sum
      DEPRED = DP5*SUM
C       
      RETURN
      END

