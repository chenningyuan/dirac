!   This subroutine applies compensation terms to residual in TRES

subroutine Compensate(TRES,TOLD,TSize,vac)

   use Parameters
   use OrbMod
   use VacMod
   use Shift
   implicit none
   integer, intent(in) :: vac
   integer, intent(in) :: TSize
   double precision, intent(inout) :: TRES(0:TSize)
   double precision, intent(out) :: TOLD(0:TSize)
   integer :: i, j, a, b
   integer :: eq, ti, ij, ab, ix
   double precision :: comp

   if (IsMainVac(vac)) return

   ! ... Singles part ...

   do i = 1, NHMO
      if (.not.IsHole(i,vac)) cycle
      ix = (i-1)*NPMO
      do a = 1, NPMO
         if (.not.IsPtcl(a,vac)) cycle
         if (IsActHole(i-NCMO,vac).and.IsActPtcl(a,vac)) cycle
         eq = ix + a
         comp = ShiftEn(i) - ShiftEn(a+NCMO)
         TRES(eq) = TRES(eq) + Beeta*comp*TOLD(eq)
      enddo
   enddo

   ! ... Doubles part of energy denominator ...

   do i = 2, NHMO
      if (.not.IsHole(i,vac)) cycle
      do j = 1, (i-1)
         if (.not.IsHole(j,vac)) cycle
         ij = MAXMIN(i,j)
         ti = TIND(ij)
         do a = 2, NPMO
            if (.not.IsPtcl(a,vac)) cycle
            do b = 1, (a-1)
               if (.not.IsPtcl(b,vac)) cycle
               if (IsActHole(i-NCMO,vac).and.IsActHole(j-NCMO,vac)    &
                   .and.IsActPtcl(a,vac).and.IsActPtcl(b,vac)) cycle
               ab = MAXMIN(a,b)
               eq = ti + ab
               comp = ShiftEn(i) + ShiftEn(j) - ShiftEn(a+NCMO) - ShiftEn(b+NCMO)
               TRES(eq) = TRES(eq) + Beeta*comp*TOLD(eq)
            enddo
         enddo
      enddo
   enddo

   return

contains
#include "maxmin.h"
!include "tind.h"

end subroutine Compensate
