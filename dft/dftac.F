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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck dftac */
      SUBROUTINE DFTAC(VX,WEIGHT,
     &                 RHO,RHO13,RHO_RESP,RHO_OUTER,GRADIENT_NORM)
C=======================================================================
C     common routine for GRAC and SAOP
C
C     SAOP
C     statistical averaging of (model) orbital potentials
C     P. R. T. Schipper, O. V. Gritsenko, S. J. A. van Gisbergen, E. J.
C     Baerends, JCP 112, 1344 (2000).
C
C     GRAC
C     gradient regulated asymptotic correction
C     M. Gr\"uning et al. JCP 114, 652 (2001).
C
C     closed-shell version
C
C     modelled after Andreas Hesselman's DFT_GRAC
C     which was finally eaten by this routine,
C     so many thanks to Andreas!
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "pi.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dftcom.h"
#include "dftac.h"
#include "dftbuf.h"
#include "dgroup.h"
#include "functionals.h"
C
      LOGICAL   LBIT
      DIMENSION VX(DERVS1),
     &          VX_BULK(DERVS1),
     &          VX_OUTER(DERVS1),
     &          VX_SLATER(DERVS1),
     &          VX_RESP(DERVS1)
C
C
C     GRAC parameters
C     ===============
C
      IF(LBIT(LBIT_DFTAC,1)) THEN !GRAC
        RHO43 = RHO*RHO13
        X = D0
        IF(RHO43 .NE. D0) THEN
          X = SQRT(GRADIENT_NORM*GRADIENT_NORM)/RHO43
        ELSE
          RETURN
        ENDIF
      ENDIF
C
C
C     shift bulk?
C     ===========
C
      SHIFTBULK = D0
      IF(LBIT(LBIT_DFTAC,1)) SHIFTBULK = DFTAC_IP + HOMO_EIGENVALUE !GRAC
C
C
C     bulk potential
C     ==============
C
      CALL DZERO(VX_BULK,DERVS1)
C
      IF(LBIT(LBIT_DFTAC,4)) THEN !SAOP with response part
        CALL DZERO(VX_SLATER,DERVS1)
        CALL DZERO(VX_RESP,DERVS1)
C
        CALL GLLBHOLEPOT(VX_SLATER,WEIGHT,RHO,GRADIENT_NORM)
        VX_RESP(FR0) = 0.42D0*WEIGHT*RHO_RESP/RHO
        VX_BULK(FR0) = VX_SLATER(FR0) + VX_RESP(FR0)
      ELSE !GRAC
        VX_BULK(FR0) = VX(FR0)
      ENDIF
C
      VX_BULK(FZ0) = VX(FZ0)
C
C
C     asymptotic potential
C     ====================
C
      CALL DZERO(VX_OUTER,DERVS1)
C
      IF(LBIT(LBIT_DFTAC,6))
     &                CALL LB94POT(VX_OUTER,WEIGHT,RHO,GRADIENT_NORM)
      IF(LBIT(LBIT_DFTAC,7))
     &                CALL LBALPHAPOT(VX_OUTER,WEIGHT,RHO,GRADIENT_NORM)
C
C
C     howto FADE (with D0 .LE. FADE .LE. D1)
C     ======================================
C
!     SAOP
      IF(LBIT(LBIT_DFTAC,3))
     &                     FADE = RHO_OUTER/RHO
!     GRAC
      IF(LBIT(LBIT_DFTAC,1))
     &                     FADE = D1/(D1+EXP(-GRAC_ALPHA*(X-GRAC_BETA)))
C
C
C     connection
C     ==========
C
      CALL DZERO(VX,DERVS1)
C
      VX(FZ0) = (D1-FADE)* VX_BULK(FZ0)
      VX(FR0) = (D1-FADE)*(VX_BULK(FR0)-WEIGHT*SHIFTBULK)
     &        +     FADE * VX_OUTER(FR0)*(D1-HFXFAC)
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck dftac_eigenvalues */
      SUBROUTINE DFTAC_EIGENVALUES(EIGENVALUES,
     &                             OLD_HOMO_EIGENVALUE)
C=======================================================================
C     check orbital energy convergence
C     for asymptotic correction
C     modelled after Andreas Hesselmann's routine HOMO
C     which has been incorporated into this one, so thanks Andreas!
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dftcom.h"
#include "dftac.h"
C
      LOGICAL LBIT
      DIMENSION EIGENVALUES(*)


!     get the homo eigenvalue
!     =======================

      if (.not. lbit(lbit_dftac, 2)) then !grac is not active
        homo_eigenvalue = -1.0d20
        do i = 1, nfsym
          if(nish(i) .gt. 0) then
            do j = 1, nish(i)
              e = eigenvalues((i-1)*norb(1) + npsh(i) + j)
              if(e .gt. homo_eigenvalue) homo_eigenvalue = e
            end do
          end if
        end do
      end if
C
C
C     turn asymptotic correction on ?
C     ===============================
C
      IF(LBIT(LBIT_DFTAC,3)) CALL SET_LBIT(LBIT_DFTAC,4)
C
      DELTA_HOMO = ABS(OLD_HOMO_EIGENVALUE - HOMO_EIGENVALUE)
      OLD_HOMO_EIGENVALUE = HOMO_EIGENVALUE
C
      IF(DELTA_HOMO .LE. DFTAC_THRESHOLD) THEN
        IF(LBIT(LBIT_DFTAC,1)) CALL SET_LBIT(LBIT_DFTAC,2)
      ENDIF
C
!     GRAC is active
      IF(LBIT(LBIT_DFTAC,2))
     &                WRITE(LUPRI,*) 'GRAC asymptotic correction active'
!     SAOP is active
      IF(LBIT(LBIT_DFTAC,4))
     &                WRITE(LUPRI,*) 'SAOP asymptotic correction active'
C
C     ICOUNTER_DFTAC_DEBUG = ICOUNTER_DFTAC_DEBUG + 1 !for debug only
C     IF(ICOUNTER_DFTAC_DEBUG .GT. 15) LBIT_DFTAC = LBIT_DFTAC + 128
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck dftac_saop_get_dmats */
      SUBROUTINE DFTAC_SAOP_GET_DMATS(DMAT_OUTER_AO,
     &                                DMAT_RESP_AO,
     &                                WORK,LWORK)
C=======================================================================
C     get orbital energies and density matrices for SAOP
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "nuclei.h"
#include "dftcom.h"
#include "dftac.h"
#include "dfterg.h"
#include "dgroup.h"
#include "functionals.h"
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LOGICAL   LBIT,
     &          FILE_IS_THERE,
     &          FILE_IS_OPENED
      DIMENSION DMAT_OUTER_AO(*),
     &          DMAT_RESP_AO(*),
     &          WORK(LWORK)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "memint.h"
      CALL QENTER('DFTAC_SAOP_GET_DMATS')
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C
C     allocate memory
C     ===============
C
      CALL MEMGET('REAL',KSAOPCMO          ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_OUTER_MO ,N2ORBX,
     &            WORK,KFREE,LFREE)
C             SAOP with response part
      IF(LBIT(LBIT_DFTAC,5)) THEN
        CALL MEMGET('REAL',KDMAT_RESP_MO,N2ORBX,
     &              WORK,KFREE,LFREE)
      ENDIF
      CALL MEMGET('REAL',KEIG              ,NORBT,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSAOPENERGIES     ,NORBT,
     &            WORK,KFREE,LFREE)
C
C
C     read coefficients and give info
C     ===============================
C
      INQUIRE(FILE = 'DFCOEF',EXIST = FILE_IS_THERE)
      IF(.NOT.FILE_IS_THERE) THEN
        WRITE(LUPRI,'(/A)') 'SAOP: Coefficient file DFCOEF not found!'
        CALL QUIT('SAOP: Coefficient file DFCOEF not found!')
      ENDIF
C
      INQUIRE(FILE = 'DFCOEF',OPENED = FILE_IS_OPENED)
      IF(FILE_IS_OPENED) THEN
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
     &            TOTERG,6,WORK,KFREE,LFREE)
      IF(.NOT.FILE_IS_OPENED) THEN
        CLOSE(LUSAOP,STATUS = 'KEEP')
      ENDIF


!     get the homo eigenvalue
!     =======================

      homo_eigenvalue = -1.0d20
      do i = 1, nfsym
        if(nish(i) .gt. 0) then
          do j = 1, nish(i)
            e = work(keig + j - 1 + (i-1)*norb(1) + npsh(i))
            if(e .gt. homo_eigenvalue) homo_eigenvalue = e
          end do
        end if
      end do
C
C
C     calculate energy differences
C     ============================
C
      K = 0
      CALL DZERO(WORK(KSAOPENERGIES),NORBT)
      DO I = 1,NFSYM
        IF(NISH(I).GT.0) THEN
          DO J = 1,NISH(I)
            K = K+1
            WORK(KSAOPENERGIES+K-1) =
     &            HOMO_EIGENVALUE-WORK(KEIG+J-1+(I-1)*NORB(1)+NPSH(I))
          ENDDO
        ENDIF
      ENDDO
C
C
C     get orbital density matrices
C     ============================
C
      CALL DZERO(DMAT_OUTER_AO,N2BBASX)
      CALL DZERO(WORK(KDMAT_OUTER_MO),N2ORBX)
C             SAOP with response part
      IF(LBIT(LBIT_DFTAC,5)) THEN
        CALL DZERO(DMAT_RESP_AO,N2BBASX)
        CALL DZERO(WORK(KDMAT_RESP_MO),N2ORBX)
      ENDIF
      K = 0
      DO I = 1,NFSYM
        IF(NISH(I).GT.0) THEN
          DO J = 1,NISH(I)
            K = K+1
            JOFF = (NORB(1)*(I-1)+NPSH(I))*NORBT
            JOFF = JOFF+NORB(1)*(I-1)+NPSH(I)
            JOFF = JOFF+(J-1)*NORBT+J-1
C                   SAOP with response part
            IF(LBIT(LBIT_DFTAC,5)) THEN
              WORK(KDMAT_RESP_MO+JOFF)
     &                                 = SQRT(WORK(KSAOPENERGIES+K-1))
            ENDIF
            WORK(KDMAT_OUTER_MO+JOFF)
     &                        = EXP(-D2*(WORK(KSAOPENERGIES+K-1)**D2))
          ENDDO
        ENDIF
      ENDDO
C
C
C     to avoid numerical problems for degenerate orbitals
C     ===================================================
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
     &                                  -WORK(KSAOPENERGIES+KCHECK-1))
            IF(SOMEDIFF.LT.1.0D-5) THEN
C                     SAOP with response part
              IF(LBIT(LBIT_DFTAC,5)) THEN
                WORK(KDMAT_RESP_MO+JOFF)
     &                               = WORK(KDMAT_RESP_MO+JOFFNEXT)
              ENDIF
              WORK(KDMAT_OUTER_MO+JOFF)
     &                              = WORK(KDMAT_OUTER_MO+JOFFNEXT)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C
C     transform density matrices
C     ==========================
C
      DO I = 1,NFSYM
        IF(NISH(I).GT.0) THEN
C                 SAOP with response part
          IF(LBIT(LBIT_DFTAC,5)) THEN
            CALL QTRANS('MOAO','S',D0,
     &                  NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &                  DMAT_RESP_AO(1+I2BASX(I,I)),
     &                  NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &                  WORK(KDMAT_RESP_MO+I2ORBX(I,I)),
     &                  NORBT,NORBT,1,IPQTOQ(1,0),
     &                  WORK(KSAOPCMO+ICMOQ(I)),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                  WORK(KSAOPCMO+ICMOQ(I)),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                  WORK(KFREE),LFREE,0)
          ENDIF
          CALL QTRANS('MOAO','S',D0,
     &                NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &                DMAT_OUTER_AO(1+I2BASX(I,I)),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &                WORK(KDMAT_OUTER_MO+I2ORBX(I,I)),
     &                NORBT,NORBT,1,IPQTOQ(1,0),
     &                WORK(KSAOPCMO+ICMOQ(I)),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                WORK(KSAOPCMO+ICMOQ(I)),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,0)
        ENDIF
      ENDDO
C             SAOP with response part
      IF(LBIT(LBIT_DFTAC,5)) THEN
        CALL DSCAL(N2BBASX,D2,DMAT_RESP_AO,1)
      ENDIF
      CALL DSCAL(N2BBASX,D2,DMAT_OUTER_AO,1)
C
C
C     memory release
C     ==============
C
      CALL MEMREL('DFTAC_SAOP_GET_DMATS',WORK,KSAOPCMO,KSAOPCMO,
     &            KFREE,LFREE)
      CALL QEXIT('DFTAC_SAOP_GET_DMATS')
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
