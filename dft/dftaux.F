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

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"

      parameter (d0 = 0.0d0,d1 = 1.0d0,d2 = 2.0d0,d3 = 3.0d0,d4 = 4.0d0)

#include "inforb.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dgroup.h"
#else
#include "wrkrsp.h"
#include "dftinf.h"
#endif
#include "dftcom.h"
#include "functionals.h"

      logical   is_gga,
     &          lbit,
     &          allow_nonhermitian
      dimension dtrmat(nbast,nbast,nz,nxcmat),
     &          excmat(nbast,nbast,nz,nxcmat),
     &          isymop(nxcmat),
     &          ihrmop(nxcmat),
     &          gao(nbast),
     &          gao1(nbast,3),
     &          buf(*),
!
     &          gr_0(3),
     &          gr_b(3),
!
     &          vxc_vector(3),
!
     &          vxc(99)

      do im = 1,nxcmat
        irepop = isymop(im) - 1

        if(ihrmop(im) .eq. -1) goto 99
        if((ihrmop(im) .lt. 1) .and. (.not. allow_nonhermitian)) goto 99

        call getrho(r_b,irepop,gao,
     &              dtrmat(1,1,1,im),buf)

        if(.not. is_gga) then

!-------------------------------------------------------------------------------
!         LDA

          if(dabs(r_b) .gt. dfthrl) then

            call dftpot1(vxc,w,r_0,d0)
            vxc_scalar = vxc(fRR)*r_b

            call doexc1(vxc_scalar,irepop,excmat(1,1,1,im),gao)
          endif
        else

!-------------------------------------------------------------------------------
!         GGA

            call getgrho(gr_b,irepop,gao,gao1,
     &                   dtrmat(1,1,1,im),buf)

          worth_it = max(
     &                   dabs(r_b),
     &                   dabs(gr_b(1)),
     &                   dabs(gr_b(2)),
     &                   dabs(gr_b(3))
     &                  )

          if(worth_it .gt. dfthrl) then
            call dftpot1(vxc,w,r_0,gr_0_norm)

            Z_b = d2*ddot(3,gr_0,1,gr_b,1)

            vxc_scalar = vxc(fRR)*r_b
     &                 + vxc(fRZ)*Z_b

            do j = 1,3
              vxc_vector(j) = vxc(fZZ)*Z_b*gr_0(j)
     &                      + vxc(fRZ)*r_b*gr_0(j)
     &                      + vxc(fZ0)*gr_b(j)
            enddo

!           this factor is *not* for the symmetrization later
            call dscal(3,d2,vxc_vector,1)

            call doexc2(vxc_scalar,vxc_vector,
     &                  irepop,excmat(1,1,1,im),gao,gao1)
          endif
        endif
 99     continue
      enddo

      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dft_qr(dtrmat,excmat,isymop,ihrmop,gao,gao1,
     &                  rho,nrmgrd,grdrho,wght,dogga,buf)
!***********************************************************************
!
!     Calculate the DFT contributions the quadratic response.
!
!     Input:
!     dtrmat(*,*,*,1) --- B density matrix
!     dtrmat(*,*,*,2) --- C density matrix
!     dtrmat(*,*,*,3) --- BC density matrix
!     isymop          --- operator symmetry
!     ihrmop          --- operator hermiticity
!     gao             ---
!     gao1            ---
!     rho             --- \rho in current grid point
!     nrmgrd          --- |\nabla\rho| in current grid point
!                         (corresponds to rhogrd in LR)
!     grdrho          --- \nabla\rho in current grid point
!                         (corresponds to rhg in LR)
!     wght            ---
!     dogga           --- LDA or GGA
!     buf             --- scratch memory
!
!     Output:
!     excmat(*,*,*,1) --- G^{[B]C} matrix
!     excmat(*,*,*,2) --- G^{[C]B} matrix
!     excmat(*,*,*,3) --- H^{[BC]} matrix
!
!     Written by johhe Jan 2007
!
!***********************************************************************
!
      implicit none
!
#include "mxcent.h"
#include "inforb.h"
#if defined (PRG_DIRAC)
#include "dgroup.h"
#else
#include "wrkrsp.h"
#include "dftinf.h"
#endif
#include "dftcom.h"
#include "functionals.h"
!
!     I/O VARIABLES
!
      double precision dtrmat(nbast,nbast,nz,3)
      double precision excmat(nbast,nbast,nz,3)
      integer isymop(3), ihrmop(3)
      double precision gao(nbast), gao1(nbast,3)
      double precision rho, nrmgrd, grdrho(3), wght, buf(*)
      logical dogga
!
!     INTERNAL VARIABLES
!
!     Symmetry of the different operators
      integer irepbop, irepcop, irepbcop
!
!     One and two times modified densities w.r.t. perturbations B and C
!     (rho_B corresponds to b0 in LR)
      double precision rho_B, rho_C, rho_BC
!
!     Gradients of the modified densities
!     (grdrho_B corresponds to b3 in LR)
      double precision grdrho_B(3), grdrho_C(3), grdrho_BC(3)
!
!     \zeta^X=|\nabla\rho^X|^2 for the different modified densities
!     (zeta_B corresponds to br in LR)
      double precision zeta_B, zeta_C, zeta_BC
!
!     xc-potential (functional derivatives of e_{xc} w.r.t. \rho and \zeta)
      double precision vxc(dervs3)
!
!     s and v in the equation
!       f_{pq} = s \Omega_{pq} + \mathbf{v} \cdot \nabla\Omega_{pq}
      double precision sxc, vt(3)
!
!     Pre-factors associated with different parts of \mathbf{v}
      double precision facr, facb, facc
!
!     To check if a point needs to be included or not.
      double precision bmax, cmax, bcmax
!
!     Loop index
      integer i
!
!     Calculate the ireps for the different perturbing operators
      irepbop  = isymop(1)-1
      irepcop  = isymop(2)-1
      irepbcop = isymop(3)-1
!
!     Get perturbed densities
!     =======================
!
!     Modified density w.r.t. B
      call getrho(rho_B,irepbop,gao,dtrmat(1,1,1,1),buf)
!
!     Modified density w.r.t. C
      call getrho(rho_C,irepcop,gao,dtrmat(1,1,1,2),buf)
!
!     Two times modified density w.r.t. B and C
      call getrho(rho_BC,irepbcop,gao,dtrmat(1,1,1,3),buf)
!
!       ***************
!       ***** LDA *****
!       ***************
!
      if(.not.dogga) then
!
!     Get the functional derivatives
         call dftpot2(vxc,wght,rho,0D0)
!
!     Constructing g^{[C]B}
         if(dabs(rho_B).gt.dfthrl) then
!
!           Calculate sxc
            sxc = vxc(frr)*rho_B
!
!           Linear transformation
            call doexc1(sxc,irepbop,excmat(1,1,1,1),gao)
         end if
!
!     Constructing g^{[B]C}
         if(dabs(rho_C).gt.dfthrl) then
!
!           Calculate sxc
            sxc = vxc(frr)*rho_C
!
!           Linear transformation
            call doexc1(sxc,irepcop,excmat(1,1,1,2),gao)
         END IF
!
!     Constructing h^{[BC]}
         if(dabs(rho_BC).gt.dfthrl) then
!
!           Calculate sxc
            sxc = 2.0D0*vxc(frrr)*rho_B*rho_C + vxc(frr)*rho_BC
!
!           Linear transformation
            call doexc1(sxc,irepbcop,excmat(1,1,1,3),gao)
         end if
!
      else
!
!         ***************
!         ***** GGA *****
!         ***************
!
!     Get gradients of perturbed densities
!     ====================================
!
!     Gradients of the different modified densities
      if(ihrmop(1) .eq. 1) then
        call getgrho(grdrho_b,irepbop,
     &               gao,gao1,
     &               dtrmat(1,1,1,1),buf)
      else
        call get_gradient_nh(grdrho_b,
     &                       irepbop,
     &                       gao,
     &                       gao1,
     &                       dtrmat(1,1,1,1),
     &                       buf)
      endif
      if(ihrmop(2) .eq. 1) then
        call getgrho(grdrho_c,irepcop,
     &               gao,gao1,
     &               dtrmat(1,1,1,2),buf)
      else
        call get_gradient_nh(grdrho_c,
     &                       irepcop,
     &                       gao,
     &                       gao1,
     &                       dtrmat(1,1,1,2),
     &                       buf)
      endif
      if(ihrmop(3) .eq. 1) then
        call getgrho(grdrho_bc,irepbcop,
     &               gao,gao1,
     &               dtrmat(1,1,1,3),buf)
      else
        call get_gradient_nh(grdrho_bc,
     &                       irepbcop,
     &                       gao,
     &                       gao1,
     &                       dtrmat(1,1,1,3),
     &                       buf)
      endif
!
!     \zeta^B = \frac{\partial\zeta}{\partial\epsilon_B}
         zeta_B = 2.0D0*(grdrho_B(1)*grdrho(1) + grdrho_B(2)*grdrho(2) +
     &                   grdrho_B(3)*grdrho(3))
!
!     \zeta^B = \frac{\partial\zeta}{\partial\epsilon_C}
         zeta_C = 2.0D0*(grdrho_C(1)*grdrho(1) + grdrho_C(2)*grdrho(2) +
     &                   grdrho_C(3)*grdrho(3))
!
!     \zeta^{BC} = \frac{\partial^2\zeta}{\partial\epsilon_B\partial\epsilon_C}
!     Note the factor 4 which comes from that \rho^{BC} is already symmetrized,
!     and, hence, a factor 2 is enough on the first line.
         zeta_BC = 2.0D0*(grdrho_BC(1)*grdrho(1) +
     &                    grdrho_BC(2)*grdrho(2) +
     &                    grdrho_BC(3)*grdrho(3)) +
     &             4.0d0*(grdrho_B(1)*grdrho_C(1) +
     &                    grdrho_B(2)*grdrho_C(2) +
     &                    grdrho_B(3)*grdrho_C(3))
!
!     Get the functional derivatives
         call dftpot2(vxc,wght,rho,nrmgrd)
!
!     Constructing g^{[B]C}
!
         bmax = max(abs(rho_B),
     &              abs(grdrho_B(1)),abs(grdrho_B(2)),abs(grdrho_B(3)))
         if (bmax.gt.dfthrl) then
!
!           Calculate sxc
            sxc  = vxc(frr)*rho_B + vxc(frz)*zeta_B
!
!           Calculate vt
            facr = vxc(frz)*rho_B + vxc(fzz)*zeta_B
            do i=1,3
               vt(i) = 2.0D0 * (facr*grdrho(i) + vxc(fz0)*grdrho_B(i))
            end do
!
!           Linear transformation
            call doexc2(sxc,vt,irepbop,excmat(1,1,1,1),gao,gao1)
!
         end if
!
!     Constructing g^{[C]B}
!
         cmax = max(abs(rho_C),
     &              abs(grdrho_C(1)),abs(grdrho_C(2)),abs(grdrho_C(3)))
         if (cmax.gt.dfthrl) then
!
!           Calculate sxc
            sxc  = vxc(frr)*rho_C + vxc(frz)*zeta_C
!
!           Calculate vt
            facr = vxc(frz)*rho_C + vxc(fzz)*zeta_C
            do I=1,3
               vt(i) = 2.0D0 * (facr*grdrho(i) + vxc(fz0)*grdrho_C(i))
            end do
!
!           Linear transformation
            call doexc2(sxc,vt,irepcop,excmat(1,1,1,2),gao,gao1)
!
         end if
!
!     Constructing h^{[BC]}
!
         bcmax = max(bmax,cmax,abs(rho_BC),
     &            abs(grdrho_BC(1)),abs(grdrho_BC(2)),abs(grdrho_BC(3)))
         if (bcmax.gt.dfthrl) then
!
!           Calculate sxc
            sxc =  2.0D0*vxc(frrr)*rho_B*rho_C +
     &             2.0d0*vxc(frrz)*(rho_B*zeta_C+rho_C*zeta_B) +
     &             2.0D0*vxc(frzz)*zeta_B*zeta_C+vxc(frr)*rho_BC +
     &             vxc(frz)*zeta_BC
!
!           Calculate vt
            facr = 2.0D0*vxc(fzzz)*zeta_B*zeta_C +
     &             2.0d0*vxc(frzz)*(rho_B*zeta_C+rho_C*zeta_B) +
     &             2.0D0*vxc(frrz)*rho_B*rho_C + vxc(frz)*rho_BC +
     &             vxc(fzz)*zeta_BC
            facb = 2.0d0*vxc(frz)*rho_C + 2.0d0*vxc(fzz)*zeta_C
            facc = 2.0d0*vxc(frz)*rho_B + 2.0d0*vxc(fzz)*zeta_B
            do i=1,3
               vt(i) = 2.0D0 * ( facr*grdrho(i) +
     &                           facb*grdrho_B(i) + facc*grdrho_C(i) +
     &                           vxc(fz0) * grdrho_BC(i) )
            end do
!
!           Linear transformation
            call doexc2(sxc,vt,irepbcop,excmat(1,1,1,3),gao,gao1)
!
         end if
!
      end if
      end
c
      SUBROUTINE FORTWRT(str, slength)
#include "priunit.h"
      INTEGER unit, slength
      CHARACTER str(slength)
      WRITE(LUPRI, *) str
      END

      SUBROUTINE DFTSETCAM(W,BETA)
#include "implicit.h"
#include "dcbham.h"
      HFXMU = W
      HFXATT = BETA
      END
C /* Deck dft_df_lr */
