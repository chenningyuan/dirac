module Shift

   use Parameters
   use IOModule
   use Options
   use OrbMod
   use GMSMod
   use ConvMod
   implicit none
   double precision :: CoreShiftVal, VirtShiftVal, ActvShiftVal
   double precision :: Beeta
   logical, private :: Init = .FALSE.

contains

subroutine InitShift

   implicit none
   integer i
   double precision :: AvgEn

   if (Init) then
      write(UOUT,*) 'Shift Module already initialized'
      call flush(UOUT)
      call abort
   endif
   Init = .TRUE.

   ! ... Allocation for Storage of Orbital Energies ...

   if (.not.allocated(OrbEn))   allocate(OrbEn(NTMO))
   if (.not.allocated(ShiftEn)) allocate(ShiftEn(NTMO))
   if (.not.allocated(GenEn))   allocate(GenEn(NTMO))

   ! ... Prepare Orbital enegies & Shifts ...

   OrbEn = OrbEnergy

   ! ... Do not shift core and virtuals ... 
   ! ... Shift Core Type Intermediate Orbitals upwards by +CoreShiftVal
   ! ... Shift Virt Type Intermediate Orbitals upwards by -VirtShiftVal
   ! ... Shift Main Orbitals by ActvShiftVal, and optionally average their energies

   ShiftEn(1:NCMO) = dble(0)
   ShiftEn(NHMO+1:NTMO) = dble(0)

   AvgEn = dble(0)
   do i = 1, NVMO
      if (IsMainOrb(i)) AvgEn = AvgEn + OrbEn(i)
   enddo
   AvgEn = AvgEn/dble(NMAINORB)

   do i = 1, NVMO
      if (IsIntmCore(i)) ShiftEn(NCMO+i) = -CoreShiftVal
      if (IsIntmVirt(i)) ShiftEn(NCMO+i) = +VirtShiftVal
      if (IsMainOrb(i)) then
         if (OPTION(AVGVAL)) then
            ShiftEn(NCMO+i) = OrbEn(NCMO+i) - AvgEn
         else
            ShiftEn(NCMO+i) = dble(0)
         endif
         ShiftEn(NCMO+i) = ShiftEn(NCMO+i) + ActvShiftVal
      endif
   enddo

!  ShiftEn(1:NCMO)      = +CoreShiftVal
!  ShiftEn(NHMO+1:NTMO) = -VirtShiftVal

   return

end subroutine InitShift

subroutine EndShift

   implicit none

   if (.not.Init) then
      write(UOUT,*) 'Shift Module not initialized'
      call flush(UOUT)
      call abort
   endif
   Init = .FALSE.

   if (allocated(OrbEn))   deallocate(OrbEn)
   if (allocated(ShiftEn)) deallocate(ShiftEn)
   if (allocated(GenEn))   deallocate(GenEn)

   return

end subroutine EndShift

end module Shift
