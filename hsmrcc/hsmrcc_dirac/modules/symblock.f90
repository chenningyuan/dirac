module SymBlock

   implicit none

   integer NSYM_BLOCK

   type symblock
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

end module SymBlock
