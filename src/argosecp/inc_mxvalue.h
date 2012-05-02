!     # the following parameters should agree exactly with those in the
!     # seg1mn() routine.  (rls)

!     # msu    = max number of symmetry inequivalent types of atoms.
      integer    msup
      parameter (msup=203)

!     # mstu   = max number of irreps.
      integer    mstup
      parameter (mstup=8)

!     # mrcru  = max number of members in a contraction set.
      integer    mrcrup
      parameter (mrcrup=50) !v3 20->50

!     # mconu  = max number of primitives in a contraction.
      integer    mconup
      parameter (mconup=50) !v3 24->50

!     # mcu    = max number of symmetry equivalent atoms.
      integer    mcup
      parameter (mcup=50)   !v2 12->50

!     # mconsu = max number of contraction sets.
      integer    mconsp
      parameter (mconsp=200) !v2 65->100

!     # msfu   = max number of function sets.
      integer    msfup
      parameter (msfup=200)  ! 52->v1(70),v2(100),1000

      integer    mnsfup
      parameter (mnsfup=(msfup*(msfup+1))/2)

!     # mgu    = max number of operators in the
!     #          nuclear interchange group.
      integer    mgup
      parameter (mgup=49)

!     # mccu   = max number of symmetry unique center
!     #          combinations from 4 nuclei.
      integer    mccup
      parameter (mccup=182)
 
! ### from (subroutine seg1mn) ###

!     # kaords = max number of ao reduction sets.
      integer    kaordp
      parameter (kaordp=240)

!     # mgcsu  = max number of symmetry orbital transformation sets.
      integer    mgcsup
      parameter (mgcsup=240)

!     # mru    = max number of irreps in an ao reduction set.
      integer    mrup
      parameter (mrup=100) !v2 24->50

!     # mcru   = max number of core potential expansion functions.
      integer    mcrup
      parameter (mcrup=140)

!     # msfru  = max number of so sets.
      integer    msfrup
      parameter (msfrup=1000) !312,500

!     # mnru   = max number of charge dists. from a function set pair.
      integer    mnrup
      parameter (mnrup=1000) !576->1000

