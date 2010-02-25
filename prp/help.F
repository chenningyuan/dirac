      subroutine getfck_wrap(tab, print_prp_level, work, lwork)

#include "implicit.h"

      real(8) :: tab(*)
      integer :: print_prp_level
      real(8) :: work(*)
      integer :: lwork

#include "memint.h"

      call getfck(tab, print_prp_level, work, kfree, lfree)

      end
