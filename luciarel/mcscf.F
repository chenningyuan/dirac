* Controlling routine for relativistic spinor optimization with
* time-reversal adapted spinors in double group symmetry
*
      subroutine mcscf(SPININTS,HARD,SOFT,MCHARD,MCSOFT,ITRAF)
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "input.inc"
*
      dimension SPININTS(ITRAF)
      dimension HARD(MCHARD),SOFT(MCSOFT)
*
*  Calculate 1- and 2-particle density matrices for complete spinor
*  space (including inactive and virtual)
      call caldens()


*
*  Construct Fock matrices (inactive and active)
      call fockmat(FOCKMAT,SPININTS,ITRAF)


*
*  Construct spinor gradient; there are 9 different types of
*  elements: uu,ub,bu and i-a,a-v,i-v for each
      call gradvec()



*
*  Construct Hessian matrix; in principle, there are 81 different
*  types (9x9), but quite a number might be redundant
      call hessmat()



*
*  Start QNR iterative solution of MCSCF equations
*  (What method are we going to use here?)
      call quanew()






      return
      end
