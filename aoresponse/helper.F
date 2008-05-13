      subroutine fill_work(length_pp,
     &                     length_pn,
     &                     iorbital_rotation_vector_pp,
     &                     iorbital_rotation_vector_pn,
     &                     jxop_pp,
     &                     jxop_pn)

#include "implicit.h"

      dimension iorbital_rotation_vector_pp(2, *)
      dimension iorbital_rotation_vector_pn(2, *)
      dimension jxop_pp(2, *)
      dimension jxop_pn(2, *)

      do i = 1, length_pp
        jxop_pp(1, i) = iorbital_rotation_vector_pp(1, i)
        jxop_pp(2, i) = iorbital_rotation_vector_pp(2, i)
      end do

      do i = 1, length_pn
        jxop_pn(1, i) = iorbital_rotation_vector_pn(1, i)
        jxop_pn(2, i) = iorbital_rotation_vector_pn(2, i)
      end do

      end
