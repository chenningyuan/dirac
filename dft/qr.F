      subroutine sdft_qr(add_grad,
     &                   add_spin,
     &                   dtrmat,
     &                   excmat,
     &                   isymop,
     &                   ihrmop,
     &                   gao,
     &                   gao1,
     &                   r_0,
     &                   gr_0_norm,
     &                   gr_0,
     &                   w,
     &                   buf)
c=======================================================================
c     describe me ...
c-----------------------------------------------------------------------
c                                                           radovan bast
c                                               last revision: june 2007
c=======================================================================

      use dft_main
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

      logical   add_grad,
     &          add_spin
      dimension dtrmat(nbast,nbast,nz,3),
     &          excmat(nbast,nbast,nz,3),
     &          isymop(3),
     &          ihrmop(3),
     &          gao(*),
     &          gao1(*),
     &          buf(*),
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
     &          v_bc_ah(3,3),
!
     &          vxc(99),
     &          derv(99)
      real(8), parameter :: tiny = 1.0d-20
      real(8) :: f
      integer :: i

#include "ibtfun.h"

      k_start = 1
      if(dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

      call get_n(r_b,
     &           isymop(1) - 1,
     &           dtrmat(1, 1, 1, 1),
     &           buf,
     &           gao)
      call get_n(r_c,
     &           isymop(2) - 1,
     &           dtrmat(1, 1, 1, 2),
     &           buf,
     &           gao)

      if (add_spin) then
        call get_s(s_b,
     &             isymop(1) - 1,
     &             dtrmat(1, 1, 1, 1),
     &             buf,
     &             gao)
        call get_s(s_c,
     &             isymop(2) - 1,
     &             dtrmat(1, 1, 1, 2),
     &             buf,
     &             gao)
      endif

      if (add_grad) then
        call get_gn_nonhermitian(gr_b,
     &                           isymop(1) - 1,
     &                           dtrmat(1, 1, 1, 1),
     &                           buf,
     &                           gao,
     &                           gao1)
        call get_gn_nonhermitian(gr_c,
     &                           isymop(2) - 1,
     &                           dtrmat(1, 1, 1, 2),
     &                           buf,
     &                           gao,
     &                           gao1)

!       Z = \nabla R \cdot \nabla R
        Z_b  = d2*ddot(3,gr_0,1,gr_b,1)
        Z_c  = d2*ddot(3,gr_0,1,gr_c,1)
        Z_bc = d2*ddot(3,gr_b,1,gr_c,1)
      end if

      if(add_grad .and. add_spin) then
        call get_gs_nonhermitian(gs_b,
     &                           isymop(1) - 1,
     &                           dtrmat(1, 1, 1, 1),
     &                           buf,
     &                           gao,
     &                           gao1)
        call get_gs_nonhermitian(gs_c,
     &                           isymop(2) - 1,
     &                           dtrmat(1, 1, 1, 2),
     &                           buf,
     &                           gao,
     &                           gao1)

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

#ifdef MOD_AUTODIFF
      if (dft_cfg_autodiff) then

        vxc = 0.0d0
        Z   = gr_0_norm*gr_0_norm

!       lda
        do i = 1, size(dft_cfg_lda_functional_weight_array)
          f = dft_cfg_lda_functional_weight_array(i)
          if (dabs(f) > tiny) then
            call lda_eval(i, 3, r_0, derv)

      vxc(frrr) = vxc(frrr) + f*taylor_deriv(derv, 2, (/3, 0/))
      vxc(frss) = vxc(frss) + f*taylor_deriv(derv, 2, (/1, 2/))

          end if
        end do

!       gga
        do i = 1, size(dft_cfg_gga_functional_weight_array)
          f = dft_cfg_gga_functional_weight_array(i)
          if (dabs(f) > tiny) then
            call gga_eval(i, 3, r_0, z, derv)

      vxc(frrr) = vxc(frrr) + f*taylor_deriv(derv, 5, (/3, 0, 0, 0, 0/))
      vxc(frss) = vxc(frss) + f*taylor_deriv(derv, 5, (/1, 2, 0, 0, 0/))

      vxc(frrz) = vxc(frrz) + f*taylor_deriv(derv, 5, (/2, 0, 1, 0, 0/))
      vxc(frsy) = vxc(frsy) + f*taylor_deriv(derv, 5, (/1, 1, 0, 0, 1/))
      vxc(frx)  = vxc(frx)  + f*taylor_deriv(derv, 5, (/1, 0, 0, 1, 0/))
      vxc(fryy) = vxc(fryy) + f*taylor_deriv(derv, 5, (/1, 0, 0, 0, 2/))
      vxc(frz)  = vxc(frz)  + f*taylor_deriv(derv, 5, (/1, 0, 1, 0, 0/))
      vxc(frzz) = vxc(frzz) + f*taylor_deriv(derv, 5, (/1, 0, 2, 0, 0/))
      vxc(fssz) = vxc(fssz) + f*taylor_deriv(derv, 5, (/0, 2, 1, 0, 0/))
      vxc(fsy)  = vxc(fsy)  + f*taylor_deriv(derv, 5, (/0, 1, 0, 0, 1/))
      vxc(fszy) = vxc(fszy) + f*taylor_deriv(derv, 5, (/0, 1, 1, 0, 1/))
      vxc(fyy)  = vxc(fyy)  + f*taylor_deriv(derv, 5, (/0, 0, 0, 0, 2/))
      vxc(fzx)  = vxc(fzx)  + f*taylor_deriv(derv, 5, (/0, 0, 1, 1, 0/))
      vxc(fzyy) = vxc(fzyy) + f*taylor_deriv(derv, 5, (/0, 0, 1, 0, 2/))
      vxc(fzz)  = vxc(fzz)  + f*taylor_deriv(derv, 5, (/0, 0, 2, 0, 0/))
      vxc(fzzz) = vxc(fzzz) + f*taylor_deriv(derv, 5, (/0, 0, 3, 0, 0/))

          end if
        end do

        vxc = vxc*w

      else
#endif /* ifdef MOD_AUTODIFF */

      if (add_grad) then
        call dftpot2(vxc, w, r_0, gr_0_norm)
      else
        call dftpot2(vxc, w, r_0, 0.0d0)
      end if

      if (add_spin) then
        if (add_grad) then
          call sdftpot2(vxc, w, r_0, gr_0_norm)
        else
          call sdftpot2(vxc, w, r_0, 0.0d0)
        end if
      end if

#ifdef MOD_AUTODIFF
      end if
#endif

!-------------------------------------------------------------------------------
!     collect s_bc_h

      s_bc_h = d0

      s_bc_h = s_bc_h + vxc(fRRR)*(r_b*r_c)

      if(add_grad) then
        s_bc_h = s_bc_h + vxc(fRRZ)*(r_b*Z_c + Z_b*r_c)
     &                  + vxc(fRZZ)*(Z_b*Z_c)
     &                  + vxc(fRZ)*Z_bc
      endif

      if(add_spin) then
        do k = k_start,3
          s_bc_h = s_bc_h + vxc(fRSS)*(s_b(k)*s_c(k))

          if(add_grad) then
            s_bc_h = s_bc_h + vxc(fRSY)*(s_b(k)*Y_c(k) + Y_b(k)*s_c(k))
     &                      + vxc(fRYY)*(Y_b(k)*Y_c(k))
     &                      + vxc(fRX)*X_bc(k)
          endif
        enddo
      endif

!-------------------------------------------------------------------------------
!     collect s_bc_ah

      s_bc_ah = 0.0d0

      if(add_spin) then
        do k = k_start,3
          s_bc_ah(k) = s_bc_ah(k) + vxc(fRSS)*(r_b*s_c(k) + s_b(k)*r_c)

          if(add_grad) then
           s_bc_ah(k) = s_bc_ah(k) + vxc(fRSY)*(r_b*Y_c(k) + Y_b(k)*r_c)
     &                             + vxc(fSSZ)*(Z_b*s_c(k) + s_b(k)*Z_c)
     &                             + vxc(fSZY)*(Z_b*Y_c(k) + Y_b(k)*Z_c)
     &                             + vxc(fSY)*Y_bc(k)
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
     &             + vxc(fZZZ)*(Z_b*Z_c)
     &             + vxc(fRZZ)*(r_b*Z_c + Z_b*r_c)
     &             + vxc(fRRZ)*(r_b*r_c)

        v_bc_h_Z_0 = v_bc_h_Z_0
     &             + vxc(fZZ)*Z_bc
        v_bc_h_Z_b = v_bc_h_Z_b
     &             + vxc(fRZ)*r_c
     &             + vxc(fZZ)*Z_c
        v_bc_h_Z_c = v_bc_h_Z_c
     &             + vxc(fRZ)*r_b
     &             + vxc(fZZ)*Z_b

        if(add_spin) then

          v_bc_h_Y_0 = 0.0d0
          v_bc_h_Y_b = 0.0d0
          v_bc_h_Y_c = 0.0d0

          v_bc_h_X_b = 0.0d0
          v_bc_h_X_c = 0.0d0

          do k = k_start,3
            v_bc_h_Z_0  = v_bc_h_Z_0
     &                  + vxc(fZYY)*(Y_b(k)*Y_c(k))
     &                  + vxc(fSZY)*(s_b(k)*Y_c(k) + Y_b(k)*s_c(k))
     &                  + vxc(fSSZ)*(s_b(k)*s_c(k))

            v_bc_h_Y_0(k) = v_bc_h_Y_0(k)
     &                    + vxc(fYY)*Y_bc(k)
            v_bc_h_Y_b(k) = v_bc_h_Y_b(k)
     &                    + vxc(fSY)*s_c(k)
     &                    + vxc(fYY)*Y_c(k)
            v_bc_h_Y_c(k) = v_bc_h_Y_c(k)
     &                    + vxc(fSY)*s_b(k)
     &                    + vxc(fYY)*Y_b(k)

            v_bc_h_X_b(k) = v_bc_h_X_b(k)
     &                    + vxc(fRX)*r_c
     &                    + vxc(fZX)*Z_c
            v_bc_h_X_c(k) = v_bc_h_X_c(k)
     &                    + vxc(fRX)*r_b
     &                    + vxc(fZX)*Z_b

            v_bc_h_Y_0(k) = v_bc_h_Y_0(k)
     &                    + vxc(fRSY)*(r_b*s_c(k) + s_b(k)*r_c)
     &                    + vxc(fRYY)*(r_b*Y_c(k) + Y_b(k)*r_c)
     &                    + vxc(fSZY)*(Z_b*s_c(k) + s_b(k)*Z_c)
     &                    + vxc(fZYY)*(Z_b*Y_c(k) + Y_b(k)*Z_c)
          enddo
        endif

        do j = 1,3
          v_bc_h(j) = v_bc_h(j) + v_bc_h_Z_0*gr_0(j)*d2
     &                          + v_bc_h_Z_b*gr_b(j)*d2
     &                          + v_bc_h_Z_c*gr_c(j)*d2
          if(add_spin) then
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

      worth_it = max(
     &               dabs(s_bc_h),
     &               dabs(v_bc_h(1)),
     &               dabs(v_bc_h(2)),
     &               dabs(v_bc_h(3))
     &              )
      if(worth_it .gt. dft_cfg_dfthrl) then
        if (add_grad) then
          call gga_distribute(s_bc_h,
     &                        v_bc_h,
     &                        isymop(3) - 1,
     &                        excmat(1, 1, 1, 3),
     &                        gao,
     &                        gao1)
        else
          call lda_distribute(s_bc_h,
     &                        isymop(3) - 1,
     &                        excmat(1, 1, 1, 3),
     &                        gao)
        end if
      endif

      if (add_spin) then

        s_bc_ah = 2.0d0*s_bc_ah
        v_bc_ah = 2.0d0*v_bc_ah

        if (add_grad) then
          call gga_distribute_qi(s_bc_ah,
     &                           v_bc_ah,
     &                           isymop(3) - 1,
     &                           excmat(1, 1, 1, 3),
     &                           gao,
     &                           gao1)
        else
          call lda_distribute_qi(s_bc_ah,
     &                           isymop(3) - 1,
     &                           excmat(1, 1, 1, 3),
     &                           gao)
        end if
      end if

      end subroutine
