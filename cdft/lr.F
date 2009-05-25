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
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     /* sdft_lr */
      subroutine sdft_lr(dtrmat,
     &                   excmat,
     &                   nxcmat,
     &                   isymop,
     &                   ihrmop,
!
!                        allow_nonhermitian allows us to use
!                        this routine also for the evaluation
!                        of quadratic response contributions
     &                   allow_nonhermitian,
!
     &                   gao,
     &                   gao1,
     &                   r_0,
     &                   gr_0_norm,
     &                   gr_0,
     &                   w,
     &                   is_gga,
     &                   buf)
c=======================================================================
c     linear transformation, spin density contribution
c-----------------------------------------------------------------------
c                                                           radovan bast
c                                            last revision: january 2008
c=======================================================================

      use dft_cfg
      use dft_density

#include "implicit.h"
#include "priunit.h"

      parameter (d0 = 0.0d0,d1 = 1.0d0,d2 = 2.0d0,d3 = 3.0d0,d4 = 4.0d0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dcbham.h"
#include "dgroup.h"
#include "inforb.h"
#include "../dft/functionals.h"

#ifdef MOD_AUTODIFF
#include "../autodiff/xc_fortran_interface.h"
#endif

      logical   is_gga,
     &          lbit,
     &          allow_nonhermitian
      dimension dtrmat(*),
     &          excmat(*),
     &          isymop(*),
     &          ihrmop(*),
     &          gao(*),
     &          gao1(*),
     &          buf(*),
!
     &          gr_0(  3),
     &          s_b (  3),
     &          gs_b(3,3),
!
     &          vxc(99),
     &          derv(99)

      k_start = 1
      if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

      do im = 1,nxcmat
        irepop = isymop(im) - 1

        if(ihrmop(im) .eq.  1) goto 99
        if((ihrmop(im) .gt.-1) .and. (.not. allow_nonhermitian)) goto 99

        im_off = nbast*nbast*nz*(im-1)

        call get_m(.true.,
     &             s_b,
     &             isymop(im)-1,
     &             gao,
     &             dtrmat(1+im_off),
     &             buf)

        if(is_gga) then
          if(allow_nonhermitian) then
            call get_m_dervs_nh(.true.,
     &                          gs_b,
     &                          isymop(im)-1,
     &                          gao,
     &                          gao1,
     &                          dtrmat(1+im_off),
     &                          buf)
          else
            call get_gs(gs_b,
     &                  isymop(im) - 1,
     &                  dtrmat(1 + im_off),
     &                  buf,
     &                  gao,
     &                  gao1)
          endif
        endif

        worth_it = dsqrt(ddot(3,s_b,1,s_b,1))

        if(is_gga) then
          worth_it = max(
     &                   worth_it,
     &                   dabs(gs_b(1,1)),
     &                   dabs(gs_b(1,2)),
     &                   dabs(gs_b(1,3)),
     &                   dabs(gs_b(2,1)),
     &                   dabs(gs_b(2,2)),
     &                   dabs(gs_b(2,3)),
     &                   dabs(gs_b(3,1)),
     &                   dabs(gs_b(3,2)),
     &                   dabs(gs_b(3,3))
     &                  )
        endif

        if(worth_it .gt. dft_cfg_dfthrl) then

          if (is_gga) then

#ifdef MOD_AUTODIFF
            if (dft_cfg_autodiff) then

              vxc = 0.0d0
              Z   = gr_0_norm*gr_0_norm

!             remember that BLYP = Slater + Becke correction + LYP
              call lda_eval(LDA_SLATER, 2, r_0, derv)
              vxc(fss) = vxc(fss)
     &                 + taylor_deriv(derv, 2, (/0, 2/))

              call gga_eval(GGA_BECKE, 2, r_0, Z, derv)
              vxc(fss) = vxc(fss)
     &                 + taylor_deriv(derv, 5, (/0, 2, 0, 0, 0/))
              vxc(fsy) = vxc(fsy)
     &                 + taylor_deriv(derv, 5, (/0, 1, 0, 0, 1/))
              vxc(fyy) = vxc(fyy)
     &                 + taylor_deriv(derv, 5, (/0, 0, 0, 0, 2/))
              vxc(fx)  = vxc(fx)
     &                 + taylor_deriv(derv, 5, (/0, 0, 0, 1, 0/))

              call gga_eval(GGA_LYP, 2, r_0, Z, derv)
              vxc(fss) = vxc(fss)
     &                 + taylor_deriv(derv, 5, (/0, 2, 0, 0, 0/))
              vxc(fsy) = vxc(fsy)
     &                 + taylor_deriv(derv, 5, (/0, 1, 0, 0, 1/))
              vxc(fyy) = vxc(fyy)
     &                 + taylor_deriv(derv, 5, (/0, 0, 0, 0, 2/))
              vxc(fx)  = vxc(fx)
     &                 + taylor_deriv(derv, 5, (/0, 0, 0, 1, 0/))

              vxc = vxc*w

            else
#endif

              call sdftpot1(vxc, w, r_0, gr_0_norm)

#ifdef MOD_AUTODIFF
            end if
#endif

          else

#ifdef MOD_AUTODIFF
            if (dft_cfg_autodiff) then

              vxc = 0.0d0

              call lda_eval(LDA_SLATER, 2, r_0, derv)
              vxc(fss) = vxc(fss)
     &                 + taylor_deriv(derv, 2, (/0, 2/))

              call lda_eval(LDA_VWN, 2, r_0, derv)
              vxc(fss) = vxc(fss)
     &                 + taylor_deriv(derv, 2, (/0, 2/))

              vxc = vxc*w

            else
#endif

              call sdftpot1(vxc, w, r_0, 0.0d0)

#ifdef MOD_AUTODIFF
            end if
#endif

          end if

          do k = k_start,3
            call sdft_lr_collinear(vxc,
     &                             excmat(1+im_off),
     &                             isymop(im),
     &                             gao,
     &                             gao1,
     &                             k,
     &                             s_b(k),
     &                             gs_b(1,k),
     &                             gr_0,
     &                             is_gga)
          enddo
        endif
 99     continue
      enddo

      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     /* sdft_lr_collinear */
      subroutine sdft_lr_collinear(vxc,
     &                             excmat,
     &                             isymop,
     &                             gao,
     &                             gao1,
     &                             k,
     &                             s_b,
     &                             gs_b,
     &                             gr_0,
     &                             is_gga)
c=======================================================================
c     describe me ...
c-----------------------------------------------------------------------
c                                                           radovan bast
c                                               last revision: june 2007
c=======================================================================

#include "implicit.h"
#include "priunit.h"

      parameter (d0 = 0.0d0,d1 = 1.0d0,d2 = 2.0d0,d3 = 3.0d0,d4 = 4.0d0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dcbham.h"
#include "dgroup.h"
#include "inforb.h"
#include "../dft/functionals.h"

      logical   is_gga
      dimension vxc(*),
     &          excmat(*),
     &          gao(*),
     &          gao1(*),
     &          gs_b(3),
     &          gr_0(3),
     &          y_pq(3),
     &          x_pq_b(3),
     &          vxc_vector(3)


      vxc_scalar = vxc(fss)*s_b

      if(is_gga) then
        y_b = ddot(3,gr_0,1,gs_b,1)

        do j = 1,3
          x_pq_b(j) = vxc(fx)

          y_pq(j) = s_b*vxc(fsy) + y_b*vxc(fyy)
        enddo

        vxc_scalar = vxc_scalar + vxc(fsy)*y_b
      endif

      if(is_gga) then
        do j = 1,3
          vxc_vector(j) = y_pq(j)*gr_0(j)
     &                  + x_pq_b(j)*gs_b(j)*d2
        enddo

        call dscal(3,d2,vxc_vector,1) !symmetrization later in the code
      endif

      call exc_m(is_gga,
     &           vxc_scalar,
     &           vxc_vector,
     &           k,
     &           isymop,
     &           excmat,
     &           gao,
     &           gao1)

      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* VRG_LR */
      SUBROUTINE VRG_LR(DTRMAT,
     &                  EXCMAT,
     &                  NXCMAT,
     &                  ISYMOP,
     &                  IHRMOP,
     &                  GAO,
     &                  GAO1,
     &                  DENSITY,
     &                  DENSITY_GRADIENT_NORM,
     &                  DENSITY_GRADIENT,
     &                  WEIGHT,
     &                  IS_GGA,
     &                  BUF)
C=======================================================================
C                                                           radovan bast
C                                              last revision: march 2007
C=======================================================================

      use dft_cfg

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (DP13 = D1/D3)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dcbham.h"
#include "dgroup.h"
#include "inforb.h"
#include "../dft/functionals.h"

      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          ISYMOP(*),
     &          IHRMOP(*),
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          FO_VORTICITY(3)

      DO IMAT = 1,NXCMAT !loop over all frequencies
        IF(IHRMOP(IMAT) .EQ. -1) THEN

          JOFF = NBAST*NBAST*NZ*(IMAT-1)

          CALL GET_ROTVPARA(.TRUE.,
     &                      FO_VORTICITY,
     &                      DENSITY,
     &                      DENSITY_GRADIENT,
     &                      ISYMOP(IMAT)-1,
     &                      GAO,
     &                      GAO1,
     &                      DTRMAT(1+JOFF),
     &                      BUF)

          RS  = (4.1887902D0*DENSITY)**(-DP13)

          IF(RS .GT. 10.0D0) GOTO 999

          DKF = (29.6088132D0*DENSITY)**DP13

          CHI_RATIO = (D1 + 0.028D0*RS)*DEXP(-0.042D0*RS)

          G = D2*WEIGHT*(DKF/236.870506D0)*(CHI_RATIO - D1)

          CALL DSCAL(3,G,FO_VORTICITY,1)

          FO_VORTICITY_NORM = DSQRT(DDOT(3,FO_VORTICITY,1,
     &                                     FO_VORTICITY,1))

          IF(FO_VORTICITY_NORM .LT. dft_cfg_DFTHRL) GOTO 999

          CALL GET_ROTVPARA(.FALSE.,
     &                      FO_VORTICITY,
     &                      DENSITY,
     &                      DENSITY_GRADIENT,
     &                      ISYMOP(IMAT)-1,
     &                      GAO,
     &                      GAO1,
     &                      EXCMAT(1+JOFF),
     &                      BUF)

 999      CONTINUE

        ENDIF
      ENDDO

      RETURN
      END
