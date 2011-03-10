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

      SUBROUTINE DFTAC(derv,WEIGHT,
     &                 RHO,RHO13,RHO_RESP,RHO_OUTER,z_0)
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

       use dft_cfg
       use dft_functional_derv

#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "pi.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "../../dft/functionals.h"
C
      real(8) :: derv(mxbllen, 0:99)
      DIMENSION VX_BULK(DERVS1),
     &          VX_OUTER(DERVS1),
     &          VX_SLATER(DERVS1),
     &          VX_RESP(DERVS1)
C
C
C     GRAC parameters
C     ===============
C
      if (dft_cfg_grac) then
        RHO43 = RHO*RHO13
        X = D0
        IF(RHO43 .NE. D0) THEN
          X = dsqrt(z_0)/RHO43
        ELSE
          RETURN
        ENDIF
      end if
C
C
C     shift bulk?
C     ===========
C
      shiftbulk = 0.0d0
      if (dft_cfg_grac) then
        shiftbulk = dft_cfg_ac_ip + dft_cfg_homo_eigenvalue
      end if
C
C
C     bulk potential
C     ==============
C
      vx_bulk = 0.0d0
C
      IF (dft_cfg_saop_with_response_part) THEN
        vx_slater = 0.0d0
        vx_resp = 0.0d0
C
        CALL GLLBHOLEPOT(VX_SLATER,WEIGHT,RHO,dsqrt(z_0))
        VX_RESP(FR0) = 0.42D0*WEIGHT*RHO_RESP/RHO
        VX_BULK(FR0) = VX_SLATER(FR0) + VX_RESP(FR0)
      ELSE
        VX_BULK(FR0) = derv(1, d1000000)
      END IF
C
      VX_BULK(FZ0) = derv(1, d0010000)
C
C
C     asymptotic potential
C     ====================
C
      VX_OUTER = 0.0d0
C
      IF(dft_cfg_asymptote_is_lb94)
     &                CALL LB94POT(VX_OUTER,WEIGHT,RHO,dsqrt(z_0))
      IF(dft_cfg_asymptote_is_lbalpha)
     &                CALL LBALPHAPOT(VX_OUTER,WEIGHT,RHO,dsqrt(z_0))
C
C
C     howto FADE (with D0 .LE. FADE .LE. D1)
C     ======================================
C
!     SAOP
      IF(dft_cfg_saop)
     & FADE = RHO_OUTER/RHO
!     GRAC
      IF(dft_cfg_grac)
     & FADE = D1/(D1+EXP(-dft_cfg_GRAC_ALPHA*(X-dft_cfg_GRAC_BETA)))
C
C
C     connection
C     ==========
C
!     VX = 0.0d0
C
      derv(1, d0010000) = (D1-FADE)* VX_BULK(FZ0)
      derv(1, d1000000) = (D1-FADE)*(VX_BULK(FR0)-WEIGHT*SHIFTBULK)
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
       use dft_cfg
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
C
      DIMENSION EIGENVALUES(*)


!     get the homo eigenvalue
!     =======================

      if (.not. dft_cfg_grac_is_active) then
        dft_cfg_homo_eigenvalue = -1.0d20
        do i = 1, nfsym
          if(nish(i) .gt. 0) then
            do j = 1, nish(i)
              e = eigenvalues((i-1)*norb(1) + npsh(i) + j)
              if(e .gt. dft_cfg_homo_eigenvalue) 
     &            dft_cfg_homo_eigenvalue = e
            end do
          end if
        end do
      end if
C
C
C     turn asymptotic correction on ?
C     ===============================
C
      if (dft_cfg_saop) dft_cfg_saop_is_active = .true.
C
      DELTA_HOMO = ABS(OLD_HOMO_EIGENVALUE - dft_cfg_homo_eigenvalue)
      OLD_HOMO_EIGENVALUE = dft_cfg_homo_eigenvalue
C
      if (delta_homo .le. dft_cfg_ac_threshold) then
        if (dft_cfg_grac) dft_cfg_grac_is_active = .true.
      end if

!     GRAC is active
      if (dft_cfg_grac_is_active)
     &                WRITE(LUPRI,*) 'GRAC asymptotic correction active'
!     SAOP is active
      if (dft_cfg_saop_is_active)
     &                WRITE(LUPRI,*) 'SAOP asymptotic correction active'

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck dftac_saop_get_dmats */
      SUBROUTINE DFTAC_SAOP_GET_DMATS(DMAT_OUTER_AO,
     &                                DMAT_RESP_AO)
C=======================================================================
C     get orbital energies and density matrices for SAOP
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
       use dft_cfg
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
#include "dgroup.h"
#include "../../dft/functionals.h"
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LOGICAL   FILE_IS_THERE,
     &          FILE_IS_OPENED
      DIMENSION DMAT_OUTER_AO(n2bbasx),
     &          DMAT_RESP_AO(n2bbasx)
      real(8), allocatable :: ksaopcmo(:)
      real(8), allocatable :: dmat_outer_mo(:)
      real(8), allocatable :: dmat_resp_mo(:)
      real(8), allocatable :: keig(:)
      real(8), allocatable :: ksaopenergies(:)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C

!     allocate memory
      allocate(ksaopcmo(      n2bbasxq))
      allocate(dmat_outer_mo( n2orbx))
      if (dft_cfg_saop_with_response_part) then
        allocate(dmat_resp_mo( n2orbx))
      end if
      allocate(keig(          norbt))
      allocate(ksaopenergies( norbt))
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
      KSAOPCMO = 0.0d0
      keig = 0.0d0
      CALL REACMO_no_work(LUSAOP,'DFCOEF',KSAOPCMO,KEIG,TEMP,
     &            TOTERG,6)
      IF(.NOT.FILE_IS_OPENED) THEN
        CLOSE(LUSAOP,STATUS = 'KEEP')
      ENDIF


!     get the homo eigenvalue
!     =======================

      dft_cfg_homo_eigenvalue = -1.0d20
      do i = 1, nfsym
        if(nish(i) .gt. 0) then
          do j = 1, nish(i)
            e = keig(1 + j - 1 + (i-1)*norb(1) + npsh(i))
            if(e .gt. dft_cfg_homo_eigenvalue) 
     &          dft_cfg_homo_eigenvalue = e
          end do
        end if
      end do
C
C
C     calculate energy differences
C     ============================
C
      K = 0
      KSAOPENERGIES = 0.0d0
      DO I = 1,NFSYM
        IF(NISH(I).GT.0) THEN
          DO J = 1,NISH(I)
            K = K+1
            KSAOPENERGIES(k) =
     &            dft_cfg_homo_eigenvalue
     &           -KEIG(J+(I-1)*NORB(1)+NPSH(I))
          ENDDO
        ENDIF
      ENDDO
C
C
C     get orbital density matrices
C     ============================
C
      dmat_outer_ao = 0.0d0
      dmat_outer_mo = 0.0d0
      if (dft_cfg_saop_with_response_part) then
        dmat_resp_ao = 0.0d0
        dmat_resp_mo = 0.0d0
      end if
      K = 0
      DO I = 1,NFSYM
        IF(NISH(I).GT.0) THEN
          DO J = 1,NISH(I)
            K = K+1
            JOFF = (NORB(1)*(I-1)+NPSH(I))*NORBT
            JOFF = JOFF+NORB(1)*(I-1)+NPSH(I)
            JOFF = JOFF+(J-1)*NORBT+J-1
            IF(dft_cfg_saop_with_response_part) THEN
              DMAT_RESP_MO(1+JOFF)
     &                                 = SQRT(KSAOPENERGIES(K))
            ENDIF
            DMAT_OUTER_MO(1+JOFF)
     &                        = EXP(-D2*(KSAOPENERGIES(K)**D2))
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
            SOMEDIFF = ABS(KSAOPENERGIES(1+KCHECK)
     &                                  -KSAOPENERGIES(KCHECK))
            IF(SOMEDIFF.LT.1.0D-5) THEN
              IF(dft_cfg_saop_with_response_part) THEN
                dmat_resp_mo(1+joff) = dmat_resp_mo(1+joffnext)
              ENDIF
           dmat_outer_mo(1+joff) = dmat_outer_mo(1+joffnext)
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
          IF(dft_cfg_saop_with_response_part) THEN
            CALL QTRANS90('MOAO','S',D0,
     &                  NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &                  DMAT_RESP_AO(1+I2BASX(I,I)),
     &                  NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &                  DMAT_RESP_MO(1+I2ORBX(I,I)),
     &                  NORBT,NORBT,1,IPQTOQ(1,0),
     &                  KSAOPCMO(1+ICMOQ(I)),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                  KSAOPCMO(1+ICMOQ(I)),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                  0)
          ENDIF
          CALL QTRANS90('MOAO','S',D0,
     &                NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &                DMAT_OUTER_AO(1+I2BASX(I,I)),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &                DMAT_OUTER_MO(1+I2ORBX(I,I)),
     &                NORBT,NORBT,1,IPQTOQ(1,0),
     &                KSAOPCMO(1+ICMOQ(I)),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                KSAOPCMO(1+ICMOQ(I)),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &                0)
        ENDIF
      ENDDO
      IF(dft_cfg_saop_with_response_part) THEN
        CALL DSCAL(N2BBASX,D2,DMAT_RESP_AO,1)
      ENDIF
      CALL DSCAL(N2BBASX,D2,DMAT_OUTER_AO,1)
C
C


!     release memory
      deallocate(ksaopcmo)
      deallocate(dmat_outer_mo)
      if (dft_cfg_saop_with_response_part) then
        deallocate(dmat_resp_mo)
      end if
      deallocate(keig)
      deallocate(ksaopenergies)

      end subroutine
