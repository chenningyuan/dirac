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
C /* Deck dftac */ 
C      SUBROUTINE DFTAC(VA,DST,CRX,CRY,CRZ)
C
C     T. Helgaker
C     based on code by Handy and Tozer
C
C#include "implicit.h"
C#include "mxcent.h"
C      PARAMETER (D0 = 0.0D0) 
C      LOGICAL INSIDE 
C      DIMENSION DST(NATOMS)
C#include "dfterg.h"
C#include "dftcom.h"
C#include "dftacb.h"
C#include "nuclei.h"
C#include "codata.h"
C
C      DIMENSION RADIUS(0:18)
C      DATA RADIUS /0.0D0,0.35D0,0.35D0,
C     &      1.45D0,1.05D0,0.85D0,0.70D0,0.65D0,0.60D0,0.50D0,0.45D0,
C     &      1.80D0,1.50D0,1.25D0,1.10D0,1.00D0,1.00D0,1.00D0,1.00D0/
C
C      SHIFT = DFTIPT + DFTHMO
C
C      INSIDE = .FALSE.
C      DO 100 I = 1, NATOMS
C         BRAG   = DFTBR1*RADIUS(NINT(CHARGE(I)))/XTANG
C         DST(I) = DSQRT((CRX - CORD(1,I))**2
C     &                + (CRY - CORD(2,I))**2
C     &                + (CRZ - CORD(3,I))**2)
C         IF (DST(I).LT.BRAG) INSIDE = .TRUE.
C  100 CONTINUE
C
C      IF (.NOT.INSIDE) THEN
C         RHOT = D0
C         RINV = D0 
C         DO 200 I = 1, NATOMS
C            RHOT = RHOT + DFTCHR(I)
C            RINV = RINV + DFTCHR(I)/DST(I)
C  200    CONTINUE
C         ASYMP = -WDFTX*RINV/RHOT + SHIFT
C
C         NFRC = 0
C         FRAC = D0 
C         DO 300 I = 1, NATOMS
C            BRAG  = RADIUS(NINT(CHARGE(I)))/XTANG
C            BRAG1 = DFTBR1*BRAG
C            BRAG2 = DFTBR2*BRAG
C            IF (DST(I) .LT. BRAG2) THEN
C               NFRC = NFRC + 1
C               FRAC = FRAC + (DST(I) - BRAG1)/(BRAG2 - BRAG1)
C            END IF
C  300    CONTINUE
C
C         IF (NFRC.EQ.0) THEN
C            VA = ASYMP
C         ELSE
C            VA = VA + FRAC*(ASYMP - VA)/NFRC
C         END IF
C      END IF
C
C      VA = VA - SHIFT
C      RETURN
C      END
#ifdef oldcode
c---------------------------------------------------------------------
      SUBROUTINE dft_grac_old(vx,wght,rho,rhogrd,r)
c---------------------------------------------------------------------
c gradient regulated asymptotic correction (closed-shell version)
c see: M. Gr\"uning et al. JCP 114 (2001) p. 652
c remark: this version only does the correction for the
c         dF/drho quantity WHILE dF/dzeta is scaled with
c         switch function
#include "implicit.h"
#include "dcbham.h"
#include "dftcom.h"
#include "dftbuf.h"
#include "functionals.h"
#include "priunit.h"
      PARAMETER (b=0.05d0)
      PARAMETER(pi=3.14159265358979312d0)
      PARAMETER(p1=0.5d0,p2=40d0)
      DIMENSION vx(dervs1)

      gp(x)=-b*x**2/(1+3d0*b*x*asinh(x))

      IF(gracalpha.EQ.0d0) gracalpha=p1
      IF(gracbeta.EQ.0d0)  gracbeta=p2

      rho13=rho**(1d0/3d0)
      rho43=rho*rho13

      sigma=2d0*rhogrd*rhogrd
      x=0d0
      IF(rho43.NE.0d0) x=SQRT(sigma)/rho43

      rhoa=0.5d0*rho
      rhoa13=rhoa**(1d0/3d0)
      rhoa43=rhoa*rhoa13
      grada=0.5d0*rhogrd
      sigmaaa=grada*grada
      xa=SQRT(sigmaaa)/rhoa43

c...asymptotic potential
      vac=0d0
      CALL vdrc(v,rho13)
      vac=vac+v
      CALL vvwn(v,rho,rho13)
      vac=vac+v
      vac=vac+rhoa13*gp(xa)

c...shape correction
      d=1d0/(1d0+EXP(-gracalpha*(x-gracbeta)))

      vx(fr0)=vx(fr0)/wght
c      write(lupri,'(5e13.5)') r,d,vac*wght,
c     >      (vx(fr0)-gracdelt)*wght,vx(fr0)*wght
      vx(fr0)=(vx(fr0)-gracdelt)*(1d0-d) + vac*(1d0-hfxfac)*d
      vx(fr0)=vx(fr0)*wght
      vx(fz0)=vx(fz0)*(1d0-d)
c      write(lupri,'(5e13.5)') r,vx(fr0),vx(fz0)

      RETURN
      END
#endif
c---------------------------------------------------------------------
      REAL*8 FUNCTION asinh(x)
c---------------------------------------------------------------------
#include "implicit.h"
      asinh=dlog(x+dsqrt((x*x) + 1d0))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck dftac */
      SUBROUTINE DFTAC(VX,WGHT,RHO,RHO13,RHORESP,RHOOUTER,RHOGRD,
     &                 WORK,LWORK)
C  .....................................................................
C  .                                                                   .
C  .  DFTAC                                                            .
C  .  common routine for GRAC and SAOP                                 .
C  .                                                                   .
C  .  SAOP                                                             .
C  .  statistical averaging of (model) orbital potentials              .
C  .  P. R. T. Schipper, O. V. Gritsenko, S. J. A. van Gisbergen, E. J..
C  .  Baerends, JCP 112, 1344 (2000).                                  .
C  .                                                                   .
C  .  GRAC                                                             .
C  .  gradient regulated asymptotic correction                         .
C  .  M. Gr\"uning et al. JCP 114, 652 (2001).                         .
C  .                                                                   .
C  .  closed-shell version                                             .
C  .                                                                   .
C  .  modelled after Andreas Hesselman's DFT_GRAC                      .
C  .  which was finally eaten by this routine,                         .
C  .  so many thanks to Andreas!                                       .
C  .                                                                   .
C  .....................................................................
C  .  author:                                          last revision:  .
C  .  radovan bast                                     october 2005    .
C  .....................................................................
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
#include "pi.h"
C
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dftcom.h"
#include "dftac.h"
#include "dftbuf.h"
#include "dgroup.h"
#include "functionals.h"
C
      DIMENSION VX(DERVS1),WORK(LWORK)
C
#include "memint.h"
C
C
C     memory allocation
C     =================
C
      CALL MEMGET('REAL',KVXBULK     ,DERVS1,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVXOUTER    ,DERVS1,
     &            WORK,KFREE,LFREE)
C
      IF(DORESP) THEN
        CALL MEMGET('REAL',KVXSLATER ,DERVS1,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVXRESP   ,DERVS1,
     &              WORK,KFREE,LFREE)
      ENDIF
C
C
C     GRAC parameters
C     ===============
C
      IF(DOGRAC) THEN
        RHO43 = RHO*RHO13
        SIGMA = D2*RHOGRD*RHOGRD
        X = D0
C
        IF(RHO43.NE.D0) THEN
           X = SQRT(SIGMA)/RHO43
        ELSE
          RETURN
        ENDIF
      ENDIF
C
C
C     shift bulk?
C     ===========
C
      IF(DOGRAC) THEN
        SHIFTBULK = ACIP+EIGMAX
      ELSE
        SHIFTBULK = D0
      ENDIF
C
C
C     bulk potential
C     ==============
C
      IF(DORESP) THEN
C
        CALL DZERO(WORK(KVXSLATER),DERVS1)
        CALL GLLBHOLEPOT(WORK(KVXSLATER),WGHT,RHO,RHOGRD)  
C
        CALL DZERO(WORK(KVXRESP),DERVS1)
        DK = 0.42D0
        WORK(KVXRESP) = DK*WGHT*RHORESP/RHO
C
        CALL DZERO(WORK(KVXBULK),DERVS1)
        DO I = 1,DERVS1
          WORK(KVXBULK-1+I) = WORK(KVXSLATER-1+I)+WORK(KVXRESP-1+I)
        ENDDO
        WORK(KVXBULK+1) = VX(FZ0)
C
      ELSE
C
        CALL DZERO(WORK(KVXBULK),DERVS1)
        DO I = 1,DERVS1
          WORK(KVXBULK-1+I) = VX(I)
        ENDDO
C
      ENDIF
C
C
C     asymptotic potential
C     ====================
C
      CALL DZERO(WORK(KVXOUTER),DERVS1)
      IF(DOLBALPHA) THEN
        CALL LBALPHAPOT(WORK(KVXOUTER),WGHT,RHO,RHOGRD)
      ENDIF
      IF(DOLB94) THEN
        CALL LB94POT(WORK(KVXOUTER),WGHT,RHO,RHOGRD)
      ENDIF
C
C
C     howto FADE (with D0 .LE. FADE .LE. D1)
C     ======================================
C
      IF(DOSAOP) THEN
        FADE = RHOOUTER/RHO
      ENDIF
      IF(DOGRAC) THEN
        FADE = D1/(D1+EXP(-GRACALPHA*(X-GRACBETA)))
      ENDIF
C     IF(DODUMP) WRITE(0,*) FADE
C
C
C     connection
C     ==========
C
      CALL DZERO(VX,DERVS1)
C
      VX(FR0) = (D1-FADE)*(WORK(KVXBULK)-WGHT*SHIFTBULK)
     &        +     FADE *(D1-HFXFAC)*WORK(KVXOUTER)
C
      VX(FZ0) = (D1-FADE)*WORK(KVXBULK+1)
C
C
C     memory release
C     ==============
C
      CALL MEMREL('DFTAC',WORK,KVXBULK,KVXBULK,KFREE,LFREE)
C
C
      RETURN
      END
