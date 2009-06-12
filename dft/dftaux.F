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
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     /* dft_lr */
      subroutine dft_lr(dtrmat,
     &                  excmat,
     &                  nxcmat,
     &                  isymop,
     &                  ihrmop,
!
!                       allow_nonhermitian allows us to use
!                       this routine also for the evaluation
!                       of quadratic response contributions
     &                  allow_nonhermitian,
!
     &                  gao,
     &                  gao1,
     &                  r_0,
     &                  gr_0_norm,
     &                  gr_0,
     &                  w,
     &                  is_gga,
     &                  buf)
c=======================================================================
c     linear transformation, spin density contribution is in sdft_lr
c-----------------------------------------------------------------------
c     T. Helgaker, T. Saue
c                                  last revision: january 2008 (radovan)
c=======================================================================

      use dft_main
      use dft_cfg
      use dft_density
      use dft_cb_interface

      implicit none

#include "functionals.h"
#ifdef MOD_AUTODIFF
#include "../autodiff/xc_fortran_interface.h"
#endif

      logical   is_gga,
     &          allow_nonhermitian
      integer:: nxcmat, im, j, irepop
      integer:: isymop(nxcmat),
     &          ihrmop(nxcmat)
      real(8):: r_0, w, gr_0_norm, dabs, Z_b, vxc_scalar, ddot, r_b
      real(8):: dtrmat(nr_ao_total,nr_ao_total,nr_qu,nxcmat),
     &          excmat(nr_ao_total,nr_ao_total,nr_qu,nxcmat),
     &          gao(nr_ao_total),
     &          gao1(nr_ao_total,3),
     &          buf(*),
!
     &          gr_0(3),
     &          gr_b(3),
!
     &          vxc_vector(3),
!
     &          vxc(99),
     &          derv(99), taylor_deriv, z

      do im = 1,nxcmat
        irepop = isymop(im) - 1

        if(ihrmop(im) .eq. -1) goto 99
        if((ihrmop(im) .lt. 1) .and. (.not. allow_nonhermitian)) goto 99

        call get_n(r_b, 
     &             irepop, 
     &             dtrmat(1, 1, 1, im),
     &             buf, 
     &             gao)

        if(.not. is_gga) then

!-------------------------------------------------------------------------------
!         LDA

          if(dabs(r_b) .gt. dft_cfg_dfthrl) then

#ifdef MOD_AUTODIFF
            if (dft_cfg_autodiff) then

              vxc = 0.0d0
             
              call lda_eval(LDA_SLATER, 2, r_0, derv)
              vxc(frr) = vxc(frr)
     &                 + taylor_deriv(derv, 2, (/2, 0/))
             
              call lda_eval(LDA_VWN, 2, r_0, derv)
              vxc(frr) = vxc(frr)
     &                 + taylor_deriv(derv, 2, (/2, 0/))
             
              vxc = vxc*w

            else
#endif

              call dftpot1(vxc, w, r_0, 0.0d0)

#ifdef MOD_AUTODIFF
            end if
#endif

            vxc_scalar = vxc(fRR)*r_b

            call lda_distribute(vxc_scalar,
     &                          irepop,
     &                          excmat(1, 1, 1, im),
     &                          gao)
          endif
        else

!-------------------------------------------------------------------------------
!         GGA

          if (allow_nonhermitian) then
            call get_gn_nonhermitian(gr_b,
     &                               irepop,
     &                               dtrmat(1, 1, 1, im),
     &                               buf,
     &                               gao,
     &                               gao1)
          else
            call get_gn(gr_b,
     &                  irepop,
     &                  dtrmat(1, 1, 1, im),
     &                  buf,
     &                  gao,
     &                  gao1)
          end if

!         worth_it = max(
!    &                   dabs(r_b),
!    &                   dabs(gr_b(1)),
!    &                   dabs(gr_b(2)),
!    &                   dabs(gr_b(3))
!    &                  )

!         if(worth_it .gt. dft_cfg_dfthrl) then

#ifdef MOD_AUTODIFF
            if (dft_cfg_autodiff) then

              vxc = 0.0d0
              Z   = gr_0_norm*gr_0_norm
             
!             remember that BLYP = Slater + Becke correction + LYP
              call lda_eval(LDA_SLATER, 2, r_0, derv)
              vxc(frr) = vxc(frr)
     &                 + taylor_deriv(derv, 2, (/2, 0/))
             
              call gga_eval(GGA_BECKE, 2, r_0, Z, derv)
              vxc(frr) = vxc(frr)
     &                 + taylor_deriv(derv, 5, (/2, 0, 0, 0, 0/))
              vxc(frz) = vxc(frz)
     &                 + taylor_deriv(derv, 5, (/1, 0, 1, 0, 0/))
              vxc(fzz) = vxc(fzz)
     &                 + taylor_deriv(derv, 5, (/0, 0, 2, 0, 0/))
              vxc(fz0) = vxc(fz0)
     &                 + taylor_deriv(derv, 5, (/0, 0, 1, 0, 0/))
             
              call gga_eval(GGA_LYP, 2, r_0, Z, derv)
              vxc(frr) = vxc(frr)
     &                 + taylor_deriv(derv, 5, (/2, 0, 0, 0, 0/))
              vxc(frz) = vxc(frz)
     &                 + taylor_deriv(derv, 5, (/1, 0, 1, 0, 0/))
              vxc(fzz) = vxc(fzz)
     &                 + taylor_deriv(derv, 5, (/0, 0, 2, 0, 0/))
              vxc(fz0) = vxc(fz0)
     &                 + taylor_deriv(derv, 5, (/0, 0, 1, 0, 0/))
             
              vxc = vxc*w

            else
#endif

              call dftpot1(vxc, w, r_0, gr_0_norm)

#ifdef MOD_AUTODIFF
            end if
#endif

            Z_b = 2.0d0*ddot(3,gr_0,1,gr_b,1)

            vxc_scalar = vxc(fRR)*r_b
     &                 + vxc(fRZ)*Z_b

            do j = 1,3
              vxc_vector(j) = vxc(fZZ)*Z_b*gr_0(j)
     &                      + vxc(fRZ)*r_b*gr_0(j)
     &                      + vxc(fZ0)*gr_b(j)
            enddo

!           this factor is *not* for the symmetrization later
            call dscal(3,2.0d0,vxc_vector,1)

            call gga_distribute(vxc_scalar,
     &                          vxc_vector,
     &                          irepop,
     &                          excmat(1, 1, 1, im),
     &                          gao,
     &                          gao1)
!         endif
        endif
 99     continue
      enddo

      end
c
      subroutine fortwrt(str, slength)

      implicit none

#include "priunit.h"

#ifdef INT_STAR8
      integer(8) :: slength
#else
      integer    :: slength
#endif /* ifdef INT_STAR8 */
      character  :: str(slength)

      write(lupri, *) str

      end subroutine

      SUBROUTINE DFTSETCAM(W,BETA)
#include "implicit.h"
#include "dcbham.h"
      HFXMU = W
      HFXATT = BETA
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GTDOAV(DMAT,IFAC)
!***********************************************************************
!
!     Get average density from open shells; only real part !
!     Written by T.Saue March 2004
!
!***********************************************************************
      use dft_cb_interface
      implicit none
#include"dcbdhf.h"
      integer :: ifac, ishell
      real(8) :: DMAT(nr_ao_total, nr_ao_total, nr_qu, nfmat)
      real(8) :: fac
      IF(AOC) THEN
        DO ISHELL = 1, NOPEN
          FAC = IFAC*DF(ISHELL)
          CALL DAXPY(nr_ao_total*nr_ao_total,FAC,
     &        DMAT(1,1,1,ISHELL+1),1,DMAT,1)
        ENDDO
      ENDIF
!
      end subroutine
