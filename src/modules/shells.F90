!-----------------------------------------
! Module with variables previously defined
! at aovec.h include file
!-----------------------------------------
   module shells_module
   use primitives_module, only: MXSHEL
   implicit none
   private

   !---------------------------------------------------------------
   ! Public
   !---------------------------------------------------------------
   public :: CENT, NHKT, KHKT, KCKT, ISTBAO, NUCO, JSTRT, NSTRT
   public :: NCENT, NRCO, NUMCF, NBCH, KSTRT, LCLASS, IPTSHL
   public :: NUMCFT, KMAX, NLRGSH, NSMLSH, NORBS
   public :: MST, MBSID, MBIDBT, SEGM, SPHR

   !*****************************************************************************
   !     Data on the AO-block associated with a given shell(ISHELL):
   !       NUCO    - number of uncontracted functions
   !       NRCO    - number of contracted functions
   !     Data on a given shell in an AO-block:
   !       NCENT   - index of symmetry independent center
   !       NUMCF   - index of shell in AO-block
   !       NBCH    - index of block in AO-vector
   !       LAMN    - name of center
   !       ISTBAO  - stabiliser: basic sym. op. that do not move center
   !       NHKT    - angular quantum number (s=1,p=2,d=3 etc.)
   !       KHKT    - number of spherical (Cartesian) components
   !       KCKT    - number of Cartesian components
   !       LCLASS  - class: large component (1), small or Huckel(2), density fitting (0)
   !       CENT    - coordinates of center
   !       JSTRT   - offset for the shell
   !       MBSID and MBIDBT have been added for multiple basis sets (WK/UniKA/31-10-2002).
   !       SEGM    - segmented contraction
   !*****************************************************************************

   real(8) :: CENT(MXSHEL, 3, 2)
   integer :: NHKT(MXSHEL), KHKT(MXSHEL), KCKT(MXSHEL), ISTBAO(MXSHEL)
   integer :: NUCO(MXSHEL), JSTRT(MXSHEL), NSTRT(MXSHEL)
   integer :: NCENT(MXSHEL), NRCO(MXSHEL), NUMCF(MXSHEL)
   integer :: NBCH(MXSHEL), KSTRT(MXSHEL), LCLASS(MXSHEL)
   integer :: IPTSHL(MXSHEL), NUMCFT(MXSHEL)
   integer :: KMAX, NLRGSH, NSMLSH, NORBS
   integer :: MST(MXSHEL), MBSID(MXSHEL), MBIDBT(MXSHEL)

   logical :: SEGM(MXSHEL), SPHR(MXSHEL)

   end module shells_module
