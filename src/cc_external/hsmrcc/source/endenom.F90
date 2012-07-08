!   This subroutine adds the energy
!   denominator to TNEW/TOLD for a given
!   vaccum. All valence parts r not
!   affected.

subroutine EnDenom(T,TSize,vac)

   use Parameters
   use OrbMod
   use VacMod
   implicit none
   integer, intent(in) :: vac
   integer, intent(in) :: TSize
   double precision, intent(inout) :: T(0:TSize)
   integer :: i, j, a, b
   integer :: eq, ti, ij, ab, ix
   double precision :: denom

   ! ... residual is already compensated, so here we apply shifts ...

   GenEn = OrbEn
   if (.not.IsMainVac(vac)) GenEn = GenEn - ShiftEn

   ! ... Singles part of energy denominator ...

   do i = 1, NHMO
      if (.not.IsHole(i,vac)) cycle
      ix = (i-1)*NPMO
      do a = 1, NPMO
         if (.not.IsPtcl(a,vac)) cycle
         if (IsActHole(i-NCMO,vac).and.IsActPtcl(a,vac)) cycle
         eq = ix + a
         denom = GenEn(i) - GenEn(a+NCMO)
         T(eq) = T(eq)/denom
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
               denom = GenEn(i) + GenEn(j) - GenEn(a+NCMO) - GenEn(b+NCMO)
               T(eq) = T(eq)/denom
            enddo
         enddo
      enddo
   enddo

   return

contains

!include "maxmin.h"
#include "../modules/maxmin.h"
!include "tind.h"

end subroutine EnDenom
