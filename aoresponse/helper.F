      subroutine fill_work(do_pp,
     &                     do_pn,
     &                     length_pp,
     &                     length_pn,
     &                     orbital_rotation_vector_pp,
     &                     orbital_rotation_vector_pn,
     &                     jxop_pp,
     &                     jxop_pn)

      implicit none

      logical   do_pp
      logical   do_pn
      integer   length_pp
      integer   length_pn
      integer   orbital_rotation_vector_pp
      integer   orbital_rotation_vector_pn
      integer   jxop_pp
      integer   jxop_pn

      dimension orbital_rotation_vector_pp(2, *)
      dimension orbital_rotation_vector_pn(2, *)
      dimension jxop_pp(2, *)
      dimension jxop_pn(2, *)

      integer   i

      if (do_pp) then
        do i = 1, length_pp
          jxop_pp(1, i) = orbital_rotation_vector_pp(1, i)
          jxop_pp(2, i) = orbital_rotation_vector_pp(2, i)
        end do
      end if

      if (do_pn) then
        do i = 1, length_pn
          jxop_pn(1, i) = orbital_rotation_vector_pn(1, i)
          jxop_pn(2, i) = orbital_rotation_vector_pn(2, i)
        end do
      end if

      end

      subroutine general_density_block(density,
     &                                 block,
     &                                 get_density,
     &                                 irepop,
     &                                 gao_slice_left,
     &                                 gao_slice_right,
     &                                 buffer,
     &                                 component_left,
     &                                 component_right,
     &                                 scale_by)

      implicit none

      double precision density
      double precision block
      logical          get_density
      integer          irepop
      double precision gao_slice_left
      double precision gao_slice_right
      double precision buffer
      integer          component_left
      integer          component_right
      double precision scale_by

#include "dgroup.h"
#include "dcbham.h"
#include "dcbbas.h"

      dimension        block(ntbas(0), ntbas(0))
      dimension        gao_slice_left(ntbas(0))
      dimension        gao_slice_right(ntbas(0))
      dimension        buffer(ntbas(0))

      integer          i1, i2, n1, n2, j, irep1, irep2
      double precision f, ddot

      if (.not.((component_left + component_right > 2)
     &    .and. (mc == 1))) then

        do irep1 = 0, nbsym - 1
          n1 = nbbas(irep1, component_left)
          if (n1 > 0) then

            irep2 = ieor(irep1, irepop)
            n2 = nbbas(irep2, component_right)
            if (n2 > 0) then

              if (get_density) then
                i1 = ibbas(irep1, component_left)  + 1
                i2 = ibbas(irep2, component_right) + 1

                call dgemv('N',
     &                     nbbas(irep1, component_left),
     &                     nbbas(irep2, component_right),
     &                     1.0d0,
     &                     block(i1, i2),
     &                     ntbas(0),
     &                     gao_slice_right(i2),
     &                     1,
     &                     0.0d0,
     &                     buffer,
     &                     1)

                density = density
     &                  + scale_by*ddot(nbbas(irep1, component_left),
     &                                  gao_slice_left(i1),
     &                                  1,
     &                                  buffer,
     &                                  1)

              else
                i1 = ibbas(irep1, component_left)
                i2 = ibbas(irep2, component_right)

                do j = 1, n2
                  f = scale_by*density*gao_slice_right(i2 + j)
                  call daxpy(n1,
     &                       f,
     &                       gao_slice_left(i1 + 1),
     &                       1,
     &                       block(i1 + 1, i2 + j),
     &                       1)
                end do
              end if
            end if
          end if
        end do
      end if

      end


!     XC stuff
!     ========

      subroutine di_get_geomderiv_fxd_dft()
        call quit('di_get_geomderiv_fxd_dft not there')
      end

      subroutine di_get_magderiv_gxd_dft()
        call quit('di_get_magderiv_gxd_dft not there')
      end

      subroutine di_get_magderiv_fxd_dft()
        call quit('di_get_magderiv_fxd_dft not there')
      end

      subroutine di_get_sigma_xc_cont()
        call quit('di_get_sigma_xc_cont not there')
      end

      subroutine di_get_t_xc_cont()
        call quit('di_get_t_xc_cont not there')
      end

      subroutine di_get_geomderiv_molgrad_dft()
        call quit('di_get_geomderiv_molgrad_dft not there')
      end

      subroutine di_get_geomderiv_gxd_dft()
        call quit('di_get_geomderiv_gxd_dft not there')
      end
