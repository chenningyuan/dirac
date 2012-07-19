!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine EndHSMR
!---------------------------------------------------------------
!
!   Ending procedure of HSMRCC
!
!   Called from HSMRCC
!
!---------------------------------------------------------------
   use Parameters
   use Options
   use IOModule
   use S2MatEle
   use VacMod
   use HEFFModule
   use VINFMod
   use TAmps
   use ConvMod
   use Shift
   use DirectDiags
   use GMSMod
   implicit none
   integer :: vac
   integer :: i, a, j, b, ij, ab, ix, ti, ia
   character*7 :: flag_valid

   ! ... For Integral Generation only runs

   if (OPTION(NOHSMR)) then
      call CloseIntFiles
      call UpdateVINF
      call EndGMSMod
      call EndVacMod
      call EndIOMod
      call CloseVINF
      return
   endif

   ! ... Update data in FileVINF       ...
   ! ... Save Final cluster amplitudes ...

   call UpdateVINF
   call SaveTAMP

   ! ... Printout Final T-Amplitudes

   do vac = 1, NSTATE

      read(UTAMP,rec=vac) TNEW
      ia = 0
      write(90,*) "Vaccum # ",vac

      do i = 1, NHMO
      do a = 1, NPMO
         ia = ia + 1
         flag_valid = " VALID "
         if (.not.IsHole(i,vac)) flag_valid = "INVALID"
         if (.not.IsPtcl(a,vac)) flag_valid = "INVALID"
         if (IsActHole(i-NCMO,vac).and.IsActPtcl(a,vac)) flag_valid = "INVALID"
         write(90,10) i, a, ia, "SINGLES", flag_valid, TNEW(ia)
         10 format(1x,2(i3,1x),8x,i10,1x,A7,1x,A7,1x,F14.9)
      enddo
      enddo
      ij = 0
      do i = 2, NHMO
      do j = 1, (i-1)
         ij = ij + 1
         ti = TIND(ij)
         ab = 0
         do a = 2, NPMO
         do b = 1, (a-1)
            ab = ab + 1
            flag_valid = " VALID "

            if (.not.IsHole(i,vac)) flag_valid = "INVALID"
            if (.not.IsHole(j,vac)) flag_valid = "INVALID"
            if (.not.IsPtcl(a,vac)) flag_valid = "INVALID"
            if (.not.IsPtcl(b,vac)) flag_valid = "INVALID"

            if (IsActHole(i-NCMO,vac).and.IsActHole(j-NCMO,vac)   &
                .and.IsActPtcl(a,vac).and.IsActPtcl(b,vac)) flag_valid = "INVALID"

            write(90,20) i, j, a, b, ti+ab, "DOUBLES", flag_valid, TNEW(ti+ab)
         20 format(1x,4(i3,1x),i10,1x,A7,1x,A7,1x,F14.9)

         enddo
         enddo
      enddo
      enddo
      write(90,*)
   enddo

   ! ... End Module Routines & Deallocate Data ...

   call DeAllocTAMP
   call CloseOtherFiles
   call EndDirect
   call EndHEFF
   call EndVacMod
   call EndS2MatEle
   call EndConvMod
   call EndShift
   call DeAllocIntBufs
   call CloseIntFiles
   call CloseVINF
   call EndIOMod

   ! ... Write Ending Information before closing Output ...

   write(UOUT,'(/,4x,a,/)') '--End of HSMR Run--'
   close(UOUT) ! mi: is it necessary ?

   return

!contains
!include "tind.h"

end subroutine EndHSMR
