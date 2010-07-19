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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!dirac_copyright_end

      subroutine sdft_qr(block_length,
     &                   add_grad,
     &                   no_sdft,
     &                   mat_dim, 
     &                   dtrmat,
     &                   excmat,
     &                   isymop,
     &                   ihrmop,
     &                   gao,
     &                   gao1,
     &                   r_0,
     &                   z_0,
     &                   gr_0,
     &                   w,
     &                   buf,
     &                   derv_len,
     &                   derv)
c=======================================================================
c     describe me ...
c-----------------------------------------------------------------------
c                                                           radovan bast
c                                               last revision: june 2007
c=======================================================================

      use xc_driver
      use dft_cfg
      use density_eval
      use dft_functional_derv

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

      integer, intent(in) :: block_length
      integer, intent(in)    :: derv_len
      real(8), intent(in)    :: derv(mxbllen, 0:derv_len)
      logical   add_grad,
     &          no_sdft
      dimension dtrmat(nbast,nbast,nz,3),
     &          excmat(nbast,nbast,nz,3),
     &          isymop(3),
     &          ihrmop(3),
     &          gao(*),
     &          gao1(*),
     &          buf(*),
     &          r_0(*),
     &          z_0(*),
     &          w(*),
!
     &          gr_0 (  3),
     &          gr_b (  3),
     &          gr_c (  3),
!
     &          s_b  (  3),
     &          s_c  (  3),
!
     &          gs_b (3,3),
     &          gs_c (3,3),
!
     &          Y_b  (  3),
     &          Y_c  (  3),
     &          Y_bc (  3),
!
     &          X_bc (  3),
!
     &          v_bc_h_Y_0(3),
     &          v_bc_h_Y_b(3),
     &          v_bc_h_Y_c(3),
     &          v_bc_h_X_b(3),
     &          v_bc_h_X_c(3),
!
!               here we construct (see thesis for notation):
!
!               K_{pq,0}^(BC) = r_0^(BC) \Omega_{pq,0} + q_{0j}^(BC) \nabla_j \Omega_{pq,0}
!               K_{pq,z}^(BC) = r_z^(BC) \Omega_{pq,z} + q_{zj}^(BC) \nabla_j \Omega_{pq,z}
!               K_{pq,y}^(BC) = r_y^(BC) \Omega_{pq,y} + q_{yj}^(BC) \nabla_j \Omega_{pq,y}
!               K_{pq,x}^(BC) = r_x^(BC) \Omega_{pq,x} + q_{xj}^(BC) \nabla_j \Omega_{pq,x}
!
!               all the rest is done elsewhere
!
!               scalar h+ part: r_0^[BC]
!               s_bc_h
!               vector h+ part: q_{0j}^[BC] with j = x,y,z
     &          v_bc_h (  3),
!               scalar h- part: r_k^[BC]    with k = x,y,z
     &          s_bc_ah(  3),
!               vector h- part: q_{kj}^[BC] with k = x,y,z and j = x,y,z
     &          v_bc_ah(3,3)
      real(8), parameter :: tiny = 1.0d-20
      real(8) :: f
      integer :: i

#include "ibtfun.h"

      k_start = 1
      if(dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

      if (add_grad) then
         call get_gn_nonhermitian(r_b,
     &                            gr_b,
     &                            isymop(1) - 1,
     &                            mat_dim,           
     &                            dtrmat(1, 1, 1, 1),
     &                            buf,
     &                            gao,
     &                            gao1)
         call get_gn_nonhermitian(r_c,
     &                            gr_c,
     &                            isymop(2) - 1,
     &                            mat_dim,           
     &                            dtrmat(1, 1, 1, 2),
     &                            buf,
     &                            gao,
     &                            gao1)
      else
         call get_n(r_b,
     &              isymop(1) - 1,
     &              mat_dim,           
     &              dtrmat(1, 1, 1, 1),
     &              buf,
     &              gao)
         call get_n(r_c,
     &              isymop(2) - 1,
     &              mat_dim,           
     &              dtrmat(1, 1, 1, 2),
     &              buf,
     &              gao)
      end if

      if (.not. no_sdft) then
         if (add_grad) then
            call get_gs_nonhermitian(s_b,
     &                               gs_b,
     &                               isymop(1) - 1,
     &                               mat_dim,           
     &                               dtrmat(1, 1, 1, 1),
     &                               buf,
     &                               gao,
     &                               gao1)
            call get_gs_nonhermitian(s_c,
     &                               gs_c,
     &                               isymop(2) - 1,
     &                               mat_dim,           
     &                               dtrmat(1, 1, 1, 2),
     &                               buf,
     &                               gao,
     &                               gao1)
         else
           call get_s(s_b,
     &                isymop(1) - 1,
     &                mat_dim,           
     &                dtrmat(1, 1, 1, 1),
     &                buf,
     &                gao)
           call get_s(s_c,
     &                isymop(2) - 1,
     &                mat_dim,           
     &                dtrmat(1, 1, 1, 2),
     &                buf,
     &                gao)
         end if
      endif

      if (add_grad) then
!       Z = \nabla R \cdot \nabla R
        Z_b  = d2*ddot(3,gr_0,1,gr_b,1)
        Z_c  = d2*ddot(3,gr_0,1,gr_c,1)
        Z_bc = d2*ddot(3,gr_b,1,gr_c,1)
      end if

      if(add_grad .and. .not. no_sdft) then
!       Y = \nabla R \cdot \nabla S
        do k = k_start,3
          Y_b (k) = ddot(3,gr_0,1,gs_b(1,k),1)
          Y_c (k) = ddot(3,gr_0,1,gs_c(1,k),1)
          Y_bc(k) = ddot(3,gr_b,1,gs_c(1,k),1)
     &            + ddot(3,gr_c,1,gs_b(1,k),1)
        enddo

!       X = \nabla S \cdot \nabla S
        do k = k_start,3
          X_bc(k) = d2*ddot(3,gs_b(1,k),1,gs_c(1,k),1)
        enddo
      endif

!-------------------------------------------------------------------------------
!     collect s_bc_h

      s_bc_h = d0

      s_bc_h = s_bc_h + derv(1, d_RRR)*(r_b*r_c)

      if(add_grad) then
        s_bc_h = s_bc_h + derv(1, d_RRZ)*(r_b*Z_c + Z_b*r_c)
     &                  + derv(1, d_RZZ)*(Z_b*Z_c)
     &                  + derv(1, d_RZ)*Z_bc
      endif

      if(.not. no_sdft) then
        do k = k_start,3
          s_bc_h = s_bc_h + derv(1, d_RSS)*(s_b(k)*s_c(k))

          if(add_grad) then
       s_bc_h = s_bc_h + derv(1, d_RSY)*(s_b(k)*Y_c(k) + Y_b(k)*s_c(k))
     &                 + derv(1, d_RYY)*(Y_b(k)*Y_c(k))
     &                 + derv(1, d_RX)*X_bc(k)
          endif
        enddo
      endif

!-------------------------------------------------------------------------------
!     collect s_bc_ah

      s_bc_ah = 0.0d0

      if(.not. no_sdft) then
        do k = k_start,3
       s_bc_ah(k) = s_bc_ah(k)
     &            + derv(1, d_RSS)*(r_b*s_c(k) + s_b(k)*r_c)

          if(add_grad) then
       s_bc_ah(k) = s_bc_ah(k) 
     &            + derv(1, d_RSY)*(r_b*Y_c(k) + Y_b(k)*r_c)
     &            + derv(1, d_SSZ)*(Z_b*s_c(k) + s_b(k)*Z_c)
     &            + derv(1, d_SZY)*(Z_b*Y_c(k) + Y_b(k)*Z_c)
     &            + derv(1, d_SY)*Y_bc(k)
          endif
        enddo
      endif

!-------------------------------------------------------------------------------
!     collect vector components

      v_bc_h  = 0.0d0
      v_bc_ah = 0.0d0

      if(add_grad) then

        v_bc_h_Z_0 = d0
        v_bc_h_Z_b = d0
        v_bc_h_Z_c = d0

        v_bc_h_Z_0 = v_bc_h_Z_0
     &             + derv(1, d_ZZZ)*(Z_b*Z_c)
     &             + derv(1, d_RZZ)*(r_b*Z_c + Z_b*r_c)
     &             + derv(1, d_RRZ)*(r_b*r_c)

        v_bc_h_Z_0 = v_bc_h_Z_0
     &             + derv(1, d_ZZ)*Z_bc
        v_bc_h_Z_b = v_bc_h_Z_b
     &             + derv(1, d_RZ)*r_c
     &             + derv(1, d_ZZ)*Z_c
        v_bc_h_Z_c = v_bc_h_Z_c
     &             + derv(1, d_RZ)*r_b
     &             + derv(1, d_ZZ)*Z_b

        if(.not. no_sdft) then

          v_bc_h_Y_0 = 0.0d0
          v_bc_h_Y_b = 0.0d0
          v_bc_h_Y_c = 0.0d0

          v_bc_h_X_b = 0.0d0
          v_bc_h_X_c = 0.0d0

          do k = k_start,3
            v_bc_h_Z_0  = v_bc_h_Z_0
     &                  + derv(1, d_ZYY)*(Y_b(k)*Y_c(k))
     &                  + derv(1, d_SZY)*(s_b(k)*Y_c(k) + Y_b(k)*s_c(k))
     &                  + derv(1, d_SSZ)*(s_b(k)*s_c(k))

            v_bc_h_Y_0(k) = v_bc_h_Y_0(k)
     &                    + derv(1, d_YY)*Y_bc(k)
            v_bc_h_Y_b(k) = v_bc_h_Y_b(k)
     &                    + derv(1, d_SY)*s_c(k)
     &                    + derv(1, d_YY)*Y_c(k)
            v_bc_h_Y_c(k) = v_bc_h_Y_c(k)
     &                    + derv(1, d_SY)*s_b(k)
     &                    + derv(1, d_YY)*Y_b(k)

            v_bc_h_X_b(k) = v_bc_h_X_b(k)
     &                    + derv(1, d_RX)*r_c
     &                    + derv(1, d_ZX)*Z_c
            v_bc_h_X_c(k) = v_bc_h_X_c(k)
     &                    + derv(1, d_RX)*r_b
     &                    + derv(1, d_ZX)*Z_b

            v_bc_h_Y_0(k) = v_bc_h_Y_0(k)
     &                    + derv(1, d_RSY)*(r_b*s_c(k) + s_b(k)*r_c)
     &                    + derv(1, d_RYY)*(r_b*Y_c(k) + Y_b(k)*r_c)
     &                    + derv(1, d_SZY)*(Z_b*s_c(k) + s_b(k)*Z_c)
     &                    + derv(1, d_ZYY)*(Z_b*Y_c(k) + Y_b(k)*Z_c)
          enddo
        endif

        do j = 1,3
          v_bc_h(j) = v_bc_h(j) + v_bc_h_Z_0*gr_0(j)*d2
     &                          + v_bc_h_Z_b*gr_b(j)*d2
     &                          + v_bc_h_Z_c*gr_c(j)*d2
          if(.not. no_sdft) then
            do k = k_start,3
              v_bc_ah(j,k) = v_bc_ah(j,k) + v_bc_h_Y_0(k)*gr_0(j)
     &                                    + v_bc_h_Y_b(k)*gr_b(j)
     &                                    + v_bc_h_Y_c(k)*gr_c(j)

              v_bc_h(j) = v_bc_h(j) + v_bc_h_Y_b(k)*gs_b(j,k)
     &                              + v_bc_h_Y_c(k)*gs_c(j,k)

              v_bc_ah(j,k) = v_bc_ah(j,k) + v_bc_h_X_b(k)*gs_b(j,k)*d2
     &                                    + v_bc_h_X_c(k)*gs_c(j,k)*d2
            enddo
          endif
        enddo
      endif

!-------------------------------------------------------------------------------

      s_bc_h =     d2*s_bc_h
      call dscal(3,d2,v_bc_h,1)

      if (add_grad) then
        call gga_distribute_r(s_bc_h,
     &                        v_bc_h,
     &                        isymop(3) - 1,
     &                        mat_dim, 
     &                        excmat(1, 1, 1, 3),
     &                        gao,
     &                        gao1)
      else
        call lda_distribute_r(s_bc_h,
     &                        isymop(3) - 1,
     &                        mat_dim, 
     &                        excmat(1, 1, 1, 3),
     &                        gao)
      end if

      if (.not. no_sdft) then

        s_bc_ah = 2.0d0*s_bc_ah
        v_bc_ah = 2.0d0*v_bc_ah

        if (add_grad) then
          call gga_distribute_i(s_bc_ah,
     &                          v_bc_ah,
     &                          isymop(3) - 1,
     &                          mat_dim, 
     &                          excmat(1, 1, 1, 3),
     &                          gao,
     &                          gao1)
        else
          call lda_distribute_i(s_bc_ah,
     &                          isymop(3) - 1,
     &                          mat_dim, 
     &                          excmat(1, 1, 1, 3),
     &                          gao)
        end if
      end if

      end subroutine
