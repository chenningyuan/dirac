module pcm_gp

   use, intrinsic :: iso_c_binding

   implicit none

   public collect_nctot
   public collect_atoms 
   public set_point_group
   public timing

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

      subroutine set_point_group(nr_gen, gen1, gen2, gen3) bind(c, name='set_point_group')

#include "mxcent.h"                     
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"

      integer(c_int), intent(inout) :: nr_gen 
      integer(c_int), intent(inout) :: gen1, gen2, gen3

      nr_gen = 0 !pcm_igen(1)   
      gen1   = 0 !pcm_igen(2) 
      gen2   = 0 !pcm_igen(3)
      gen3   = 0 !pcm_igen(4)

      end subroutine set_point_group

      real(8) function timing()

      use interface_to_mpi
      
#if defined (VAR_MPI)
      timing = interface_MPI_WTIME()
#else
      call cpu_time(timing)
#endif                      

      end function timing

end module pcm_gp
