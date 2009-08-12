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
      subroutine dft_lr(block_length,
     &                  dtrmat,
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
     &                  z_0,
     &                  gr_0,
     &                  w,
     &                  is_gga,
     &                  buf)
c=======================================================================
!fixme sdft_lr and dft_lr should really be only one routine
c     linear transformation, spin density contribution is in sdft_lr
c-----------------------------------------------------------------------
c     T. Helgaker, T. Saue
c                                  last revision: january 2008 (radovan)
c=======================================================================

      use dft_main
      use dft_cfg
      use dft_density
      use dft_cb_interface
      use dft_temp
      use dft_functional_derv
      use dft_basis_info

      implicit none

      integer:: block_length
      logical   is_gga,
     &          allow_nonhermitian
      integer:: nxcmat, im, j, irepop, i
      integer:: isymop(nxcmat),
     &          ihrmop(nxcmat)
      real(8):: r_0(*), w(*), z_0(*), dabs, Z_b, vxc_scalar, ddot, r_b
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
     &          derv(mxbllen, 0:99)
      real(8), parameter :: tiny = 1.0d-20

! gosia: added
      real(8)             :: cx, cy, cz
      type(dft_options_t) :: dft_options
      real(8)             :: gaob1(nr_ao_total,3)

      do im = 1,nxcmat
        irepop = isymop(im) - 1

        if(ihrmop(im) .eq. -1) goto 99
        if((ihrmop(im) .lt. 1) .and. (.not. allow_nonhermitian)) goto 99

        call get_n(r_b, 
     &             irepop, 
     &             dtrmat(1, 1, 1, im),
     &             buf, 
     &             gao)

        if (is_gga) then
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
        end if

        call get_functional_derv(2,
     &                           block_length,
     &                           w,
     &                           r_0,
     &                           z_0,
     &                           nr_nonzero_derv(2),
     &                           derv)

        if (is_gga) then

          z_b = 2.0d0*ddot(3,gr_0,1,gr_b,1)

          vxc_scalar = derv(1, d_rr)*r_b
     &               + derv(1, d_rz)*z_b

          do j = 1,3
            vxc_vector(j) = derv(1, d_zz)*z_b*gr_0(j)
     &                    + derv(1, d_rz)*r_b*gr_0(j)
     &                    + derv(1, d_z)*gr_b(j)
          enddo

!         this factor is *not* for the symmetrization later
          call dscal(3, 2.0d0, vxc_vector, 1)

          call gga_distribute(vxc_scalar,
     &                        vxc_vector,
     &                        irepop,
     &                        excmat(1, 1, 1, im),
     &                        gao,
     &                        gao1)

        else
          vxc_scalar = derv(1, d_rr)*r_b
          call lda_distribute(vxc_scalar,
     &                        irepop,
     &                        excmat(1, 1, 1, im),
     &                        gao)
        end if

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
