C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgrdinf */
      SUBROUTINE RGRDINF(GCI,GOEE,GOEP)
C*****************************************************************************
C
C     Calculate norm of gradient
C
C     Input:
C        GCI     - the CI gradient of the CEP
C        GOEE    - the e-e orbital gradient of the CEP
C        GOEP    - the e-p gradient of the CEP
C
C     Output:
C        GNORM   - (1) CI-gradient norm
C                  (2) orbital e-e gradient norm
C                  (3) orbital e-p gradient norm
C                  (4) total orbital gradient norm
C                  (5) total norm
C
C     Written by J. Thyssen - Nov 3 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D0 = 0.0D00)
C
#include <dcbopt.h>
#include <dcbkrmc_itinfo.h>
C
      DIMENSION GCI(*),GOEE(*),GOEP(*)
C
      CALL QENTER('RGRDINF')
C
      IF (NZCONF .GT. 0) THEN
         GCINRM = DNORM2(NZCONFQ,GCI,1)
      ELSE
         GCINRM = D0
      END IF
      IF (NZXOPT .GT. 0) THEN
         IF (NZXOPE .GT. 0) THEN
            GOBEENRM = DNORM2(NZXOPEQ,GOEE,1)
         ELSE
            GOBEENRM = D0
         END IF
         IF (NZXOPP .GT. 0) THEN
            GOBEPNRM = DNORM2(NZXOPPQ,GOEP,1)
         ELSE
            GOBEPNRM = D0
         END IF
      ELSE
         GOBEENRM = D0
         GOBEPNRM = D0
      END IF
      GNORM(1) = GCINRM
      GNORM(2) = GOBEENRM
      GNORM(3) = GOBEPNRM
      GNORM(4) = SQRT(GOBEENRM**2+GOBEPNRM**2)
      GNORM(5) = SQRT(GCINRM**2+GNORM(4)**2)
C
C     Save gradient information for final summary output
C
C
      CALL DCOPY(5,GNORM,1,DINFO(11),1)
C
      CALL QEXIT('RGRDINF')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rcigrad */
      SUBROUTINE RCIGRAD(CREF,FCAC,H2AC,GCI,EMCACT,WORK,LWORK)
C*****************************************************************************
C
C     Calculate CI gradient
C
C     Input:
C        CREF    - CI vector
C        FCAC    - FC Fock matrix with active indeces
C        H2AC    - Two-electron integrals with active indeces
C        EMCACT  - active energy
C
C     Output:
C        GCI     - CI gradient
C
C     Written by J. Thyssen - Nov 20 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( THRSML = 1.0D-8 )
      PARAMETER ( THRZER = 1.0D-12 )
C
#include <dcbopt.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION CREF(*)
      DIMENSION FCAC(*)
      DIMENSION H2AC(*)
      DIMENSION GCI(*)
      DIMENSION WORK(LWORK)
C
#include <memint.h>
C
      CALL QENTER('RCIGRAD')
C
C     Calculate CI gradient
C     ---------------------
C
      IF ( OPT_CIPROGRAM .EQ. 'GOSCI' ) THEN
C
         CALL RCIGRAD_GOSCI(CREF,FCAC,H2AC,GCI,EMCACX,IPROPT,
     &        WORK(KFREE),LFREE)
C
      ELSE
         WRITE(LUPRI,'(/,A,/,2A)')
     $        '*** ERROR in RCIGRAD ***',
     $        'gradient not implemented for CI program ',
     $        OPT_CIPROGRAM
         CALL QUIT('*** ERROR in RCIGRAD ***')
      END IF
C
C     Consistency check of active energy
C     ----------------------------------
C         
      DNZCONF = NZCONF
      DIFFOK = SQRT(DNZCONF) * THRZER
      IF ( ABS (EMCACT) .GT. THRZER)
     &     DIFFOK = DIFFOK * ABS(EMCACT)
      IF (NZXOPT .GT. 0 .AND. ABS(EMCACT-EMCACX) .GT. DIFFOK) THEN
         IF (ABS(EMCACT-EMCACX) .GT. THRSML) THEN
            WRITE(LUPRI,9000)
     &           '*** ERROR in RCIGRAD ***',
     &           EMCACT,EMCACX
            CALL QUIT('*** ERROR in RCIGRAD ***')
         ELSE
            WRITE(LUPRI,9000)
     &           '*** WARNING in RCIGRAD ***',
     &           EMCACT,EMCACX
         END IF
      END IF
 9000 FORMAT(1P,//,A,/,3X,'Active energy from DV and PV:',D25.15,
     &     /,3X,'Active energy from RCIGRAD:  ',D25.15,/)
C
      CALL QEXIT('RCIGRAD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rcigrad_gosci */
      SUBROUTINE RCIGRAD_GOSCI(CREF,FCAC,H2AC,GCI,EMCACT,
     &     IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Calculate CI gradient with CI program GOSCI.
C
C     Input:
C        CREF    - CI vector
C        FCAC    - FC Fock matrix with active indeces
C        H2AC    - Two-electron integrals with active indeces
C        EMCACT  - active energy
C
C     Output:
C        GCI     - CI gradient
C
C     Written by J. Thyssen - Nov 20 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( D1 = 1.0D00, DM1 = -D1, D0 = 0.0D00, D2 = 2.0D00 )
C
#include <dcbopt.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbgen.h>
C
      DIMENSION CREF(NZCONF,2)
      DIMENSION FCAC(*)
      DIMENSION H2AC(*)
      DIMENSION GCI(NZCONF,2)
      DIMENSION WORK(LWORK)
      DIMENSION EMC(2)
C
#include <memint.h>
C
      CALL QENTER('RCIGRAD_GOSCI')
C
C     ***************************************************
C     *** Calculate CI matrix times reference vector. ***
C     ***************************************************
C
C     In princinple we would do this direct, i.e. calculating a
C     CI sigma vector. Instead we do it the ``easy'' way:
C
C     - calculate CI matrix
C     - do a matrix x vector multiplication
C
C     Calculate CI matrix
C     ===================
C
C
C     Transform quarternion FCAC (Dirac MO-basis) to Molfdir MO-basis
C     ---------------------------------------------------------------
C
      CALL MEMGET('REAL',KMOLFFC,(2*NASHT)*(2*NASHT)*MIN0(NZ,2),
     $     WORK,KFREE,LFREE)
      CALL QFC2MFC(FCAC,WORK(KMOLFFC),IPROPT)
      CALL GOSCI_SETONE(WORK(KMOLFFC),2*NASHT)
      CALL MEMREL('after RGETFC',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Transform integrals from Dirac (NZ,3) basis to Molfdir basis
C     ------------------------------------------------------------
C
      CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * 2, WORK,KFREE,LFREE)
      CALL DZERO(WORK(KMUUUU),(2*NASHT)**4 * 2)
C
      CALL DNZ32M(H2AC,WORK(KMUUUU),IPROPT)
C
C     Construct Hamiltonian:
C     ----------------------
C
      NZCONF2 = NZCONF ** 2
      CALL MEMGET('REAL',KHR,NZCONF2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHI,NZCONF2,WORK,KFREE,LFREE)
C
      CALL DZERO(WORK(KHR),NZCONF2)
      CALL DZERO(WORK(KHI),NZCONF2)
      CALL GOSCI_MAKEH(WORK(KHR),WORK(KHI), WORK(KMUUUU))
C
C     Calculate CI matrix times reference vector
C     ------------------------------------------
C
C     Real part: real(H) * real(CREF) - imag(H) * imag(CREF)
C
C     real(H) * real(CREF):
C
      CALL DGEMM('N','N',NZCONF,1,NZCONF,D1,
     &     WORK(KHR),NZCONF,
     &     CREF(1,1),NZCONF,
     &     D0,GCI(1,1),NZCONF)
C
      IF (NZ .GE. 2) THEN
C
C        imag(H) * imag(CREF):
C
         CALL DGEMM('N','N',NZCONF,1,NZCONF,DM1,
     &        WORK(KHI),NZCONF,
     &        CREF(1,2),NZCONF,
     &        D1,GCI(1,1),NZCONF)
C
C
C        Imaginary part: real(H) * imag(CREF) + imag(H) * real(CREF)
C
C        real(H) * imag(CREF):
C
         CALL DGEMM('N','N',NZCONF,1,NZCONF,D1,
     &        WORK(KHR),NZCONF,
     &        CREF(1,2),NZCONF,
     &        D0,GCI(1,2),NZCONF)
C
C        imag(H) * real(CREF):
C
         CALL DGEMM('N','N',NZCONF,1,NZCONF,D1,
     &        WORK(KHI),NZCONF,
     &        CREF(1,1),NZCONF,
     &        D1,GCI(1,2),NZCONF)
C
      END IF
C
      IF ( IPRINT .GE. 0 ) THEN
         CALL HEADER('Output from RCIGRAD_GOSCI:',-1)
         WRITE(LUPRI,'(1X,A)') 'H * CREF, real part'
         CALL OUTPUT(GCI(1,1),1,NZCONF,1,1,NZCONF,1,1,LUPRI)
         IF (NZ .GE. 2) THEN
            WRITE(LUPRI,'(1X,A)') 'H * CREF, imag part'
            CALL OUTPUT(GCI(1,2),1,NZCONF,1,1,NZCONF,1,1,LUPRI)
         END IF
      END IF
C
C
C     Calculate active energy
C     -----------------------
C
C     Real part of energy:
C
      EMC(1) = DDOT(NZCONF,CREF(1,1),1,GCI(1,1),1)
C
      IF (NZ .GE. 2) THEN
         EMC(1) = EMC(1) - DDOT(NZCONF,CREF(1,2),1,GCI(1,2),1)
         EMC(2) =
     &        DDOT(NZCONF,CREF(1,1),1,GCI(1,2),1)
     &        + DDOT(NZCONF,CREF(1,2),1,GCI(1,1),1)
      ELSE
         EMC(2) = D0
      END IF
      IF (IPRINT .GE. 0 .OR. EMC(2) .NE. D0) THEN
         WRITE(LUPRI,'(1P,2/,1X,A,2(/,3X,A,D25.15),2/)')
     &        'Active energy from CI gradient module:',
     &        'Real part: ', EMC(1),
     &        'Imag part: ', EMC(2)
         IF ( EMC(2) .NE. D0) THEN
            WRITE(LUPRI,'(2/,A,/,A,2/)')
     &           '*** WARNING ***',
     &           'Active energy is imaginary!'
         END IF
      END IF
C
C
      EMCACT = EMC(1)
C
C     Finish the CI gradient
C     ----------------------
C
C     GCI -= EMCACT * CREF
C
      DO IZ = 1, MIN0(NZ,2)
         DO I = 1, NZCONF
            GCI(I,IZ) = GCI(I,IZ) - EMCACT * CREF(I,IZ)
         END DO
      END DO
C
      IF (IPRINT .GE. 0) THEN
         WRITE(LUPRI,'(1X,A)') 'CI gradient, real part'
         CALL OUTPUT(GCI(1,1),1,NZCONF,1,1,NZCONF,1,1,LUPRI)
         IF (NZ .GE. 2) THEN
            WRITE(LUPRI,'(1X,A)') 'CI gradient, imag part'
            CALL OUTPUT(GCI(1,2),1,NZCONF,1,1,NZCONF,1,1,LUPRI)
         END IF
      END IF
C
C     Save diagonal of Hessian
C     ------------------------
C
C     E^{[2]c^*c}_{\mu,\mu} = <\mu|H|\mu> - E[0]
C        - ( g_{\mu}^* C_{\mu}^* + C_{\mu} g_{\mu} ) =
C
C     E^{[2]c^*c}_{\mu,\mu} = <\mu|H|\mu> - E[0] 
C        - 2 Re ( g_{\mu} C_{\mu} )
C
      CALL MEMGET('REAL',KHCDIAG,NZCONF,WORK,KFREE,LFREE)
C
      DO I = 1, NZCONF
         WORK(KHCDIAG + I - 1) =
     &        WORK(KHR + (I-1) + (I-1) * NZCONF) - EMCACT
     &        - D2 * GCI(I,1) * CREF(I,1)
         IF (NZ .GE. 2) 
     &        WORK(KHCDIAG + I - 1) = 
     &        WORK(KHCDIAG + I - 1) + D2 * GCI(I,2) * CREF(I,2)
      END DO
      IF (IPRINT .GE. 0) THEN
         WRITE(LUPRI,'(1X,A)') 'E2 conf. diagonal'
         CALL OUTPUT(WORK(KHCDIAG),1,NZCONF,1,1,NZCONF,1,1,LUPRI)
      END IF
      CALL WRTKRMC(LUKRM3,'E2CDIAG ',WORK(KHCDIAG),NZCONF)
C
      CALL MEMREL('E2 diag',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RCIGRAD_GOSCI')
      RETURN
      END
