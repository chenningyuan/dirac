! This subroutine zeros the elements
! of TNEW/TOLD(i->a;vac) and TNEW/TOLD(i,j->a,b;vac)
! which are not valid excitations
! for the given vaccum. Elements which
! are valid excitations are kept untouched.
! this routine is required if indexing
! program has to be vaccum independent and
! keep it very simple. proper use of this
! routine can minimize errors due to
! unnecessary multiplications which have to be
! carried out, but do not(should not) contribute.

subroutine MakeZero(T,TSize,vac)

   use Parameters
   use VacMod
   implicit none
   integer, intent(in) :: vac
   integer, intent(in) :: TSize
   double precision, intent(inout) :: T(0:TSize)
   integer :: i, j, a, b
   integer :: eq, ti, ij, ab, ix

   T(0) = dble(0)
   ! ... Singles ...
   do i = 1, NHMO
      ix = (i-1)*NPMO
      if (.not.IsHole(i,vac)) then
         T(ix+1:ix+NPMO) = dble(0)
         cycle
      endif
      do a = 1, NPMO
         eq = ix + a
         if (.not.IsPtcl(a,vac)) then
            T(eq) = dble(0)
            cycle
         endif
         if (IsActHole(i-NCMO,vac).and.IsActPtcl(a,vac)) T(eq) = dble(0)
      enddo
   enddo
   ! ... Doubles ...
   ij = 0
   do i = 2, NHMO
   do j = 1, (i-1)
      ij = ij + 1
      ti = TIND(ij)
      if ((.not.IsHole(i,vac)).or.(.not.IsHole(j,vac))) then
         T(ti+1:ti+NPMOSQ) = dble(0)
         cycle
      endif
      ab = 0
      do a = 2, NPMO
      do b = 1, (a-1)
         ab = ab + 1
         eq = ti + ab
         if ((.not.IsPtcl(a,vac)).or.(.not.IsPtcl(b,vac))) then
            T(eq) = dble(0)
            cycle
         endif
         if (IsActHole(i-NCMO,vac).and.IsActPtcl(a,vac).and. &
             IsActHole(j-NCMO,vac).and.IsActPtcl(b,vac))     &
             T(eq) = dble(0)
      enddo
      enddo
   enddo
   enddo

   return

!contains
!include "tind.h"

end subroutine MakeZero
