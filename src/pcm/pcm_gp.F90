module pcm_gp

   use, intrinsic :: iso_c_binding

   implicit none

   public collect_nctot
   public collect_atoms 

   private

   contains

      subroutine collect_nctot(nr_nuclei) bind(c, name='collect_nctot') 

#include "mxcent.h"
#include "nuclei.h"
      
      integer(c_int), intent(out) :: nr_nuclei
      
      nr_nuclei = nucdep

      end subroutine collect_nctot
      
      subroutine collect_atoms(atomic_charges, atomic_centers) bind(c, name='collect_atoms')

#include "mxcent.h"
#include "nuclei.h"
      
      real(c_double), intent(out) :: atomic_charges(*)
      real(c_double), intent(out) :: atomic_centers(3,*)
     
      integer :: i, j, k 

! Get coordinates
      call getacord(atomic_centers)
! Get charges      
      i = 0
      do j = 1, nucind
        do k = 1, nucdeg(j)
          i = i + 1
          atomic_charges(i) = charge(j)
        enddo
      enddo
   
      end subroutine collect_atoms

end module pcm_gp
