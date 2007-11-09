module OrbMod

   use Parameters
   use IOModule
   implicit none
   double precision, allocatable, dimension(:) :: OrbEn
   double precision, allocatable, dimension(:) :: GenEn
   double precision, allocatable, dimension(:) :: ShiftEn
   integer, allocatable :: SpinOfOrb(:), SpinOfActvOrb(:)
   integer, allocatable :: OccGMSRef(:), OccRef(:)
   integer, allocatable :: IrrepOrb(:), IrpLabOrb(:)
   integer, allocatable :: OrbType(:), OrbSpace(:)
   double precision, allocatable :: OrbEnergy(:), OrbCoeff(:,:)
   double precision, allocatable :: HCoreMatrix(:,:)

contains

function IsMainOrb(i)

   implicit none
   integer i
   logical IsMainOrb

   if ((i.le.0).or.(i.gt.NVMO)) then
      write(UOUT,*) "Invalid value for IsMainOrb"
      call flush(UOUT)
      call abort
   endif

   IsMainOrb = OrbSpace(i+NCMO).eq.PM_SPACE

   return

end function IsMainOrb

function IsIntmOrb(i)

   implicit none
   integer i
   logical IsIntmOrb

   if ((i.le.0).or.(i.gt.NVMO)) then
      write(UOUT,*) "Invalid value for IsIntmOrb"
      call flush(UOUT)
      call abort
   endif

   IsIntmOrb = OrbSpace(i+NCMO).eq.PI_SPACE

   return

end function IsIntmOrb

function IsIntmCore(i)

   implicit none
   integer i
   logical IsIntmCore

   if ((i.le.0).or.(i.gt.NVMO)) then
      write(UOUT,*) "Invalid value for IsIntmCore"
      call flush(UOUT)
      call abort
   endif

   IsIntmCore = IsIntmOrb(i)
   if (IsIntmCore) IsIntmCore = OccRef(i+NCMO).eq.1

   return

end function IsIntmCore

function IsIntmVirt(i)

   implicit none
   integer i
   logical IsIntmVirt

   if ((i.le.0).or.(i.gt.NVMO)) then
      write(UOUT,*) "Invalid value for IsIntmVirt"
      call flush(UOUT)
      call abort
   endif

   IsIntmVirt = IsIntmOrb(i)
   if (IsIntmVirt) IsIntmVirt = OccRef(i+NCMO).eq.0

   return

end function IsIntmVirt

end module OrbMod
