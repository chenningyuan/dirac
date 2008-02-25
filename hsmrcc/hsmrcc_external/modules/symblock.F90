!module SymBlock
module SymBlocks
!------------------------------------------------------------------
!
!   Gfortran does not accept identical names of module and of
!  a derived data type inside.
!
!------------------------------------------------------------------
   implicit none

   integer NSYM_BLOCK

   type symblock ! gfortran complains due to identical module/type names 
      integer :: ms2
      integer :: irp
      integer :: nvac
      integer, allocatable :: vacs(:)
      integer, allocatable :: occp(:,:)
      double precision, allocatable :: s2mat(:,:)
   end type symblock

   type op_symblock
      double precision, allocatable :: op(:,:)
      double precision, allocatable :: eval(:)
      double precision, allocatable :: rvec(:,:), lvec(:,:)
      integer, allocatable :: s2p1(:)
   end type op_symblock

   type (symblock), allocatable, target :: symblock_array(:)

end module SymBlocks
!end module SymBlock
