#ifdef UNDEF
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rorbgrd */
      SUBROUTINE RORBGRD(DV,FC,FV,FQ,GRD,JXOP,NXOP,FPV,DFAC,
     &                   IPQ,IPRINT)
C*****************************************************************************
C
C     Calculate the orbital gradient
C
C     Input:
C        DV      - one-electron active density matrix
C                  (dimensioned N2ASHTQ)
C        FC      - the inactive Fock matrix
C                  (dimensioned full matrix: N2ORBXQ)
C        FV      - the active Fock matrix
C                  (dimensioned full matrix: N2ORBXQ)
C        FQ      - the Q Fock matrix
C                  (dimensioned ?)
C
C     Output:
C        GRD     - the e-p gradient of the CEP
C     
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
C
      DIMENSION DV(*),FC(*),FV(*),FQ(*),GRD(NXOP,NZ),
     &          JXOP(2,*),FPV(*),DFAC(2),IPQ(*)
      LOGICAL   NOFV
C
      CALL QENTER('RORBGRD')
C
C
C
C     ***************************************************
C     *** Calculate F_{pv}                          *****
C     ***    = \sum_{u} D_{vu} F_^C_{pu} + F^Q_{pv} *****
C     ***************************************************
C
C
      CALL FGEN_PV(FC,FQ,DV,FPV,IPQ,IPRINT)
C
C
C     ***************************************************
C     *** Calculate gradient elements *******************
C     ***************************************************
C
C
      NOFV = (NASHT .EQ. 0)
      DO IG = 1,NXOP
         K = JXOP(1,IG)
         L = JXOP(2,IG)
         CALL IOBTYP(K,ITYPK,ISYMK)
         CALL IOBTYP(L,ITYPL,ISYML)
C        
C        g_{ip} = F_{pi} - F_{ip}^*
C
C        g_{vp} = F_{pv} - F_{vp}^*
C
C
C
         IOFFKL = 0
         IOFF_FPV = 1
         IF (ISYMK .EQ. 2) THEN
            IOFFKL = NORB(1)
            IOFF_FPV = IOFF_FPV + NZ*NASH(1)*NORB(1)
         END IF
         CALL CALCGR(FC,FV,FPV,NORBT,NASHT,NZ,K,L,NOFV,
     &               GRD(IG,1),NXOP,ITYPK,ITYPL,
     &               IORB(ISYMK)+NPSH(ISYMK)+NISH(ISYMK),
     &               IORB(ISYML)+NPSH(ISYML)+NISH(ISYML),
     &               DFAC,IPQ)
         write(6,*) 'gradient for rotation ',K,L
         write(6,*) (grd(ig,iz),iz=1,nz)
      END DO
C
C
  200 CALL QEXIT('RORBGRD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck iobtyp */
      SUBROUTINE IOBTYP(I,ITYP,ISYM)
C*****************************************************************************
C
C     Return orbital type of orbital I
C
C     Input:
C        I       - orbital number
C
C     Output:
C        ITYP    - orbital type of orbital I
C        ISYM    - symmetry of orbital I
C
C     Written by J. Thyssen - Oct 29 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dgroup.h>
C
      II = I
      ITYP = -1
      DO IFSYM = 1,NFSYM
         IF (II .LE. NPSH(IFSYM)) THEN
            ITYP = 1
         ELSE IF (II .LE. NPSH(IFSYM)+NISH(IFSYM)) THEN
            ITYP = 2
         ELSE IF (II .LE. NPSH(IFSYM)+NISH(IFSYM)+NASH(IFSYM)) THEN
            ITYP = 3
         ELSE IF (II .LE. NPSH(IFSYM)+NISH(IFSYM)+
     &                    NASH(IFSYM)+NSSH(IFSYM)) THEN
            ITYP = 4
         END IF
         IF (ITYP .EQ. -1) THEN
            II = II-NORB(IFSYM)
         ELSE
            ISYM = IFSYM
            RETURN
         END IF
      END DO
      WRITE(LUPRI,9000) I
      WRITE(LUPRI,*) I
 9000 FORMAT(/'Orbital number: ',I3)
      CALL QUIT('IOBTYP: unknown orbital type')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fgen_pv */
      SUBROUTINE FGEN_PV(FC,FQ,DV,FGEN,IPQ,IPRINT)
C*****************************************************************************
C
C     Calculate \sum_u F^C_{pu} D_{uv}^* + F^Q_{pv}
C               = (F^C * D^\dagger + F^Q)_{pv}
C
C     Input:
C        FC      - inactive Fock matrix
C        DV      - active density matrix
C        FQ      - the Q Fock matrix
C
C     Output:
C        FGEN    - generalized Fock matrix with second index active
C
C     Written by J. Thyssen - Nov 16 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D1 = 1.00D00, D0 = 0.00D00 )
C
#include <dgroup.h>
#include <dcborb.h>
      DIMENSION FC(NORBT,NORBT,NZ),FQ(NORBT,NORBT,NZ),DV(*),
     &          FGEN(NORBT,NASHT,NZ),IPQ(*)
C
      CALL QENTER('FGEN_PV')
      IF (IPRINT .GE. 5) CALL HEADER('Output from FGEN_PV',-1)
C
      IF (NASHT .GT. 1) THEN
         CALL QUIT('FGEN_PV: not implemented...')
         FAC = D1
      ELSE
         FAC = D0
      END IF
C
      IOFF_DV   = 1
      IOFF_FGEN = 1
C
CJTHMAERKE: Skal det mon vaere IPQ eller IPQTOQ(1,0)
C
      DO 10 I = 1,NFSYM
         IF (NASH(I) .GT. 0) THEN

         CALL QGEMM(NORBT,NASH(I),NASH(I),D1,
     &      'N','N',IPQTOQ(1,0),FC(1,1+IORB(I)+NPSH(I)+NISH(I),1),
     &         NORBT,NORBT,NZ,
     &      'H','N',IPQTOQ(1,0),DV(IOFF_DV),NASH(I),NASH(I),NZ,
     &      FAC,IPQTOQ(1,0),FGEN(1,1+IASH(I),1),NORBT,NASH(I),NZ)
         IF (IPRINT .GE. 5) THEN
            WRITE(LUPRI,9000) I
            CALL PRQMAT(FGEN(1,1+IASH(I),1),NORBT,NASH(I),NORBT,NASH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
         END IF
         IOFF_DV   = IOFF_DV   + N2ASH(I)*NZ
c        IOFF_FGEN = IOFF_FGEN + NORB(I)*NASH(I)*NZ
  10  CONTINUE
C
      CALL QEXIT('FGEN_PV')
 9000 FORMAT(/'Generalized Fock matrix irrep ',I1,
     &        ' with second index active')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck calcgr */
      SUBROUTINE CALCGR(FC,FV,FPV,NORBT,NASHT,NZ,K,L,NOFV,GORB,NZX,
     &                  ITYPK,ITYPL,IOFFK,IOFFL,DFAC,IPQ)
C*****************************************************************************
C
C     Calculate gradient element
C
C     Input:
C        FC      - the inactive Fock matrix
C        FV      - the active Fock matrix
C        FPV     - generalized Fock matrix with second index active
C
C     Output:
C        GORB    - gradient element
C
C     Written by J. Thyssen - Nov 17 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D0 = 0.0D00, D2 = 2.00D00)
C
      DIMENSION FC(NORBT,NORBT,NZ), FV(NORBT,NORBT,NZ),
     &          FPV(NORBT,NASHT,NZ), GORB(NZX,NZ), DFAC(*),
     &          IPQ(*)
      LOGICAL   NOFV
C
C        
C        g_{ip} = F_{pi} - F_{ip}^*
C
C        g_{vp} = F_{pv} - F_{vp}^*
C
C
C     First index:
C     
      IF (ITYPK .EQ. 2) THEN
C
C        ...inactive 
C
C        Add F_{pi}
C
         IF (NOFV) THEN
            DO IZ = 1,NZ 
               GORB(1,IZ) = DFAC(1) * D2 * FC(L,K,IZ)
            END DO
         ELSE
            DO IZ = 1,NZ 
               GORB(1,IZ) = DFAC(1) * D2 * (FC(L,K,IZ) + FV(L,K,IZ))
            END DO
         END IF
      ELSE
C
C        ...active
C
C        Add F_{pv}
C
         DO IZ = 1,NZ 
            GORB(1,IZ) = DFAC(1) * FPV(L,K-IOFFK,IZ)
         END DO
      END IF
C
C     Second index:
C     
      IF (ITYPL .EQ. 3) THEN
C
C        ...active
C
C        Add F_{pv}^*
C
         DO IZ = 1,NZ
            IP = IPQ(IZ)
            IF (IP .EQ. 1) THEN
               GORB(1,IZ) = GORB(1,IZ) + DFAC(2) * FPV(K,L-IOFFL,IZ)
            ELSE
               GORB(1,IZ) = GORB(1,IZ) - DFAC(2) * FPV(K,L-IOFFL,IZ)
            END IF
         END DO
C     ELSE
C
C        ...secondary/positronic
C
C        Subtract 0
C
C
      END IF
C
      RETURN
      END
#endif
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
