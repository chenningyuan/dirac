!-----------------------------------------
! Module with the parameters and the statement functions previously
! defined at dcbibt.h include file: four indices packed in one integer
!-----------------------------------------
   module dcbibt_mod
   implicit none
   private
   public :: IBIT2, IBIT2M, IBIT4, IBIT4M, INDGET, LGET, KGET, JGET, IGET, IPACK
#if defined (INT_STAR8)
   integer, parameter :: IBIT2 = 65536, IBIT2M = 65535, IBIT4 = 4294967296, IBIT4M = 4294967295
   integer, parameter :: NSHIFT = 16   ! bits per packed index
#else
   integer, parameter :: IBIT2 = 256, IBIT2M = 255, IBIT4 = 65536, IBIT4M = 65535
   integer, parameter :: NSHIFT = 8
#endif
   contains
   pure elemental integer function INDGET(I, J)
   integer, intent(in) :: I, J
   INDGET = IAND(ISHFT(I, -J), IBIT2M)
   end function INDGET
   pure elemental integer function LGET(I)
   integer, intent(in) :: I
   LGET = INDGET(I, 0)
   end function LGET
   pure elemental integer function KGET(I)
   integer, intent(in) :: I
   KGET = INDGET(I, NSHIFT)
   end function KGET
   pure elemental integer function JGET(I)
   integer, intent(in) :: I
   JGET = INDGET(I, 2*NSHIFT)
   end function JGET
   pure elemental integer function IGET(I)
   integer, intent(in) :: I
   IGET = INDGET(I, 3*NSHIFT)
   end function IGET
   pure elemental integer function IPACK(I, J, K, L)
   integer, intent(in) :: I, J, K, L
   IPACK = IOR(IOR(IOR(ISHFT(L, 0), ISHFT(K, NSHIFT)), ISHFT(J, 2*NSHIFT)), ISHFT(I, 3*NSHIFT))
   end function IPACK
   end module dcbibt_mod
