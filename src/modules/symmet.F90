!-----------------------------------------
! Module with variables previously defined
! at symmet.h include file
!-----------------------------------------
   module symmet_mod
   use maxaqn_mod, only: MXQN, MXAQN
   use mxcent_mod, only: MXCENT
   implicit none
!
   ! ---- Variables from old COMMON /SYMMTI/ ----
   integer :: MAXREP, MAXOPR
   integer :: MULT(0:7)
   integer :: ISYMAX(3, 2)
   integer :: ISYMAO(MXQN, MXAQN)
   integer :: NPARSU(8)
   integer :: NAOS(8)
   integer :: NPARNU(8, 8)
   integer, allocatable :: IPTSYM(:, :)   ! (AO components, 0:7)
   integer :: IPTCNT(3*MXCENT, 0:7, 2)
   integer :: NCRREP(0:7, 2)
   integer :: IPTCOR(3*MXCENT, 2)
   integer :: NAXREP(0:7, 2)
   integer :: IPTAX(3, 2)
   integer :: IPTXYZ(3, 0:7, 2)
   integer :: IPTNUC(MXCENT, 0:7)
   integer :: ISOP(0:7)
   integer :: NROTS, NINVC, NREFL
   integer :: IXVAL(0:7, 0:7)
   integer :: NCOS(8, 2)
   integer, allocatable :: ICLASS(:)      ! (SOs)
   integer :: IPAR(0:7)
   integer :: JSOP(0:7)
   integer :: ICOS(8, 2)
   integer :: I2COSX(8, 8, 2, 2)
!
   ! ---- Variables from old COMMON /SYMMTR/ ----
   real*8 :: FMULT(0:7)
   real*8 :: PT(0:7)
!
   contains

   subroutine alloc_symmet_orbitals(nao, nso)
!     (Re)allocate the AO/SO arrays, zeroed as the old COMMON /SYMMTI/ was.
      integer, intent(in) :: nao, nso
      if (allocated(IPTSYM)) deallocate(IPTSYM, ICLASS)
      allocate(IPTSYM(nao, 0:7), ICLASS(nso), source=0)
   end subroutine alloc_symmet_orbitals

   end module symmet_mod
