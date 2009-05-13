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

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"

      parameter (d0 = 0.0d0,d1 = 1.0d0,d2 = 2.0d0,d3 = 3.0d0,d4 = 4.0d0)

#include "inforb.h"
#include "nuclei.h"
#include "dgroup.h"
#include "functionals.h"

#ifdef MOD_AUTODIFF
#include "../autodiff/xc_fortran_interface.h"
#endif

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
     &          vxc(99),
     &          derv(99)

      do im = 1,nxcmat
        irepop = isymop(im) - 1

        if(ihrmop(im) .eq. -1) goto 99
        if((ihrmop(im) .lt. 1) .and. (.not. allow_nonhermitian)) goto 99

        call getrho(r_b,irepop,gao,
     &              dtrmat(1,1,1,im),buf)

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

            call doexc1(vxc_scalar,irepop,excmat(1,1,1,im),gao)
          endif
        else

!-------------------------------------------------------------------------------
!         GGA

#ifdef MOD_CDFT
          if(allow_nonhermitian) then
            call get_gr_nh(.true.,gr_b,irepop,gao,gao1,
     &                     dtrmat(1,1,1,im),buf)
          else
#endif /* MOD_CDFT */
            call getgrho(gr_b,irepop,gao,gao1,
     &                   dtrmat(1,1,1,im),buf)
#ifdef MOD_CDFT
          endif
#endif /* MOD_CDFT */

          worth_it = max(
     &                   dabs(r_b),
     &                   dabs(gr_b(1)),
     &                   dabs(gr_b(2)),
     &                   dabs(gr_b(3))
     &                  )

          if(worth_it .gt. dft_cfg_dfthrl) then

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
C /* Deck dft_df_lr */
#ifdef MOD_DNF
      SUBROUTINE fit_DFT_LR(pCFIT,pVXCint,NXCMAT,ISYMOP,IHRMOP,
     &                      GAO,GAO1,RHO,RHOGRD,RHG,WGHT,DOGGA,B0,
     &                      nskipped)
C
C adapted from DFT_LR for density fitting
C andre gomes and radovan bast, january 2007
C
      use dft_cfg
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,DP5 = 0.5D0,D2=2.0D0)
C
#include "inforb.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dgroup.h"
#else
#include "wrkrsp.h"
#include "dftinf.h"
#endif
#include "functionals.h"
C
      INTEGER A, B
      LOGICAL DOGGA, TRPLET, LBIT
      DIMENSION pCFIT(*),pVXCint(*),ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &          GAO(*), GAO1(*),
     &          RHG(3),VT(3), VXC(DERVS2)
      DIMENSION B3(3)
C
      if (NXCMAT.GT.1)
     &   CALL QUIT ('DF-VXC not working for multiple densities yet...')

      DO IM = 1,NXCMAT
      IREPOP = ISYMOP(IM)-1
      IF(IHRMOP(IM).EQ.1) THEN
C
C       Get perturbed density
C       =====================
C
        B0 = DDOT(NBASISAUX,GAO,1,pCFIT,1)
C
C       ***************
C       ***** LDA *****
C       ***************
C
        IF(.NOT.DOGGA) THEN
          IF(DABS(B0).GT.dft_cfg_DFTHRL) THEN
C
C            Calculated VT
C
             CALL DFTPOT1(VXC, WGHT,RHO, 0D0)
             SXC = VXC(FRR)*B0
C
C            Linear transformation, saving the integrals
C
             CALL fit_DOEXC(DOGGA,SXC,DUMMY,GAO,DUMMY,
     &                      pVXCint,NBASISAUX)
          else
             nskipped=nskipped+1
          END IF
        ELSE
C
C         ***************
C         ***** GGA *****
C         ***************
C
C  getting gradient from perturbed fitted density
C
           DO IXYZ = 1,3
              IOFF = NBASISAUX * (IXYZ - 1) + 1
              b3(ixyz) =
     &              D1*DDOT(NBASISAUX,GAO1(IOFF),1,pCFIT,1)
           ENDDO
           BMAX = MAX(ABS(B0),ABS(B3(1)),ABS(B3(2)),ABS(B3(3)))
C
           IF (BMAX.GT.dft_cfg_DFTHRL) THEN
              CALL DFTPOT1(VXC,WGHT,RHO,RHOGRD)
              BR = (B3(1)*RHG(1)+ B3(2)*RHG(2)+ B3(3)*RHG(3))/RHOGRD
              FACR = VXC(FRZ)*B0 + (VXC(FZZ)-VXC(FZ0)/RHOGRD)*BR
              SXC = VXC(FRR)*B0 + VXC(FRZ)*BR
              DO I = 1,3
                 VT(I) = (VXC(FZ0)*B3(I)+FACR*RHG(I))/RHOGRD
              ENDDO
c        write (*,*) 'debug: '
c        write (*,'(A3,E20.12)') '       b3(i,j,k)',b3(1),b3(2),b3(3)
c        write (*,'(A3,E20.12)') '      rhg(i,j,k)',rhg(1),rhg(2),rhg(3)
c        write (*,'(A3,E20.12)') '       vt(i,j,k)',vt(1),vt(2),vt(3)
c        write (*,'(A2,E20.12)') ' sxc,b0,br,f1,f2',sxc,b0,br,
c    &                              VXC(FRR)*B0,VXC(FRZ)*BR
C
C           Linear transformation, saving the integrals
C
              CALL fit_DOEXC(DOGGA,SXC,VT,GAO,GAO1,pVXCint,NBASISAUX)
          else
             nskipped=nskipped+1
          END IF
        END IF

      ENDIF
      ENDDO
      RETURN
      END
#endif /* MOD_DNF */
