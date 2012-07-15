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

   integer, allocatable :: MOLTRA_indx(:) ! for MOLTRA indexes

   double precision, allocatable :: OrbEnergy(:), OrbCoeff(:,:)
   double precision, allocatable :: HCoreMatrix(:,:)

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function IsMainOrb(i)
   implicit none
   integer i
   logical IsMainOrb

   if ((i.le.0).or.(i.gt.NVMO)) then
      write(UOUT,*) "Invalid value for IsMainOrb"
      write(UOUT,*) "i,NVMO:",i,NVMO
      call flushout
      call quit('Invalid value for IsMainOrb')
   endif

  ! if (.not.allocated(OrbSpace).or.sizeof(OrbSpace) < i+NCMO ) then
   if (.not.allocated(OrbSpace).or.size(OrbSpace) < i+NCMO ) then
     write(UOUT,*) 'OrbSpace error !...size(OrbSpace)=',size(OrbSpace)
     write(UOUT,*) 'IsMainOrb(i): i=',i,' i+NCMO=',i+NCMO
     call quit('IsMainOrb: error with OrbSpace!')
   endif

   IsMainOrb = OrbSpace(i+NCMO).eq.PM_SPACE

   return
end function IsMainOrb

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function IsIntmOrb(i)
   implicit none
   integer i
   logical IsIntmOrb

   if ((i.le.0).or.(i.gt.NVMO)) then
      write(UOUT,*) "Invalid value for IsIntmOrb"
      write(UOUT,*) "i,NVMO:",i,NVMO
      call flushout
      call quit('Invalid value for IsIntmOrbInvalid value for IsIntmOrb')
   endif

   IsIntmOrb = OrbSpace(i+NCMO).eq.PI_SPACE

   return
end function IsIntmOrb

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function IsIntmCore(i)
   implicit none
   integer i
   logical IsIntmCore

   if ((i.le.0).or.(i.gt.NVMO)) then
     write(UOUT,*) "Invalid value for IsIntmCore"
     write(UOUT,*) "i,NVMO:",i,NVMO
     call flushout
     call quit('Invalid value for IsIntmCore')
   endif

   IsIntmCore = IsIntmOrb(i)
   if (IsIntmCore) IsIntmCore = OccRef(i+NCMO).eq.1

   return
end function IsIntmCore

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function IsIntmVirt(i)
   implicit none
   integer i
   logical IsIntmVirt

   if ((i.le.0).or.(i.gt.NVMO)) then
     write(UOUT,*) "Invalid value for IsIntmVirt:"
     write(UOUT,*) "i,NVMO:",i,NVMO
     call flushout
     call quit('Invalid value for IsIntmVirt')
   endif

   IsIntmVirt = IsIntmOrb(i)
   if (IsIntmVirt) IsIntmVirt = OccRef(i+NCMO).eq.0

   return
end function IsIntmVirt

subroutine print_orbmod_functions
 implicit none
 integer i

 write(UOUT,'(/,6x,a,i3,a,i3,a)') '*** Output from OrbMod ***,(1-NVMO=',NVMO,' NCMO=',NCMO,')'
 do i=1, NVMO
   write(UOUT,'(2x,i3,a,l1,a,l1,a,l1,a,l1)')  & 
    i,'. IsMainOrb:',IsMainOrb(i), ' IsIntmOrb:',IsIntmOrb(i),  & 
    ' IsIntmCore:',IsIntmCore(i),' IsIntmVirt:',IsIntmVirt(i)
 enddo
 call flushout

end subroutine print_orbmod_functions

end module OrbMod
