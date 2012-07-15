module VacMod

   use Parameters
   use IOModule
   use GMSMod
   use OrbMod
   !use SymBlock
   use SymBlocks
   use s2matele
   implicit none
   public
   logical, private :: Init = .FALSE.
   integer, allocatable :: VacString(:), VacTable(:,:)
   integer, allocatable :: OccVac(:,:)
   integer, allocatable :: IrpVac(:), Ms2Vac(:)
   integer, allocatable :: OrderVac(:), RevOrderVac(:)
   logical, allocatable :: IsMainVac(:)
   interface ExcVaccum
      module procedure ExcVaccum1
      module procedure ExcVaccum2
   end interface
   interface Exct
      module procedure Exct1
      module procedure Exct2
   end interface
   private :: HoleNum, PtclNum
   private :: GenVac
   private :: MAXMIN

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine InitVacMod
!---------------------------------------------------------------------------
!
!                  Initialize vacuum module
!
!     called from source/inihsbase.F90 - routine InitHSMRBase
!
!---------------------------------------------------------------------------
   implicit none


   if (get_iprint() >= 7) then
    call header('Output from "InitVacMod & GenVac"',-1)
    write(UOUT,'(2x,a,i4)')  &
    '...output from InitVacMod - parameter for allocations, NSTATE=',NSTATE
   endif

   if (Init) then
     write(UOUT,'(/,2x,a)') 'InitVacMod: VacMod already initialized !'
     call flushout
     call Quit('InitVacMod: VacMod already initialized !')
   endif
   Init = .TRUE.

!   ... MI: TODO - add allocation checks   

   if (.not.allocated(VacString))   allocate(VacString(NSTATE))
   if (.not.allocated(VacTable))    allocate(VacTable(NSTATE,NSTATE))
   if (.not.allocated(OccVac))      allocate(OccVac(NVMO,NSTATE))
   if (.not.allocated(IrpVac))      allocate(IrpVac(NSTATE))
   if (.not.allocated(Ms2Vac))      allocate(Ms2Vac(NSTATE))
   if (.not.allocated(OrderVac))    allocate(OrderVac(NSTATE))
   if (.not.allocated(RevOrderVac)) allocate(RevOrderVac(NSTATE))
   if (.not.allocated(IsMainVac))   allocate(IsMainVac(NSTATE))

   call GenVac

   return

end subroutine InitVacMod

subroutine EndVacMod

   implicit none

   if (.not.Init) then
      write(UOUT,*) 'EndVacMod: VacMod not initialized !'
      !call flush(UOUT)
      call flushout
      !call abort
      call Quit('EndVacMod: VacMod not initialized !')
   endif
   Init = .FALSE.

 ! ... TODO: include proper mem check
   if (allocated(VacString))   deallocate(VacString)
   if (allocated(VacTable))    deallocate(VacTable)
   if (allocated(OccVac))      deallocate(OccVac)
   if (allocated(IrpVac))      deallocate(IrpVac)
   if (allocated(Ms2Vac))      deallocate(Ms2Vac)
   if (allocated(OrderVac))    deallocate(OrderVac)
   if (allocated(RevOrderVac)) deallocate(RevOrderVac)
   if (allocated(IsMainVac))   deallocate(IsMainVac)

   return

end subroutine EndVacMod

function IsHole(i,vac)

   implicit none
   integer, intent(in) :: i, vac
   logical :: IsHole

   IsHole = (i<=NCMO)
   if (.not.IsHole) IsHole = BTEST(VacString(vac),(i-NCMO))

   return

end function IsHole

function IsPtcl(a,vac)

   implicit none
   integer, intent(in) :: a, vac
   logical :: IsPtcl

   IsPtcl = (a>NVMO)
   if (.not.IsPtcl) IsPtcl = .NOT.BTEST(VacString(vac),a)

   return

end function IsPtcl

function IsActHole(i,vac)

   implicit none
   integer, intent(in) :: i, vac
   logical :: IsActHole

   IsActHole = (i>0).and.(BTEST(VacString(vac),i))

   return

end function IsActHole

function IsActPtcl(a,vac)

   implicit none
   integer, intent(in) :: a, vac
   logical :: IsActPtcl

   IsActPtcl = (a<=NVMO).and.(.NOT.BTEST(VacString(vac),a))

   return

end function IsActPtcl

function HoleNum(i,vac)

   implicit none
   integer, intent(in) :: i, vac
   integer :: HoleNum
   integer :: j, k

   !if (.not.IsActHole(i,vac)) stop 'HoleNum : index not act hole'
   if (.not.IsActHole(i,vac)) then
     write(UOUT,*) 'HoleNum error ! i=',i,'  vac=',vac, &
     ' IsActHole(i,vac):',IsActHole(i,vac)
     call flushout
     call Quit('HoleNum : index not active hole !')
   endif

   k = 0
   do j = 1, i
      if (IsActHole(j,vac)) k = k + 1
   enddo
   HoleNum = k

   return

end function HoleNum

function PtclNum(a,vac)

   implicit none
   integer, intent(in) :: a, vac
   integer :: PtclNum
   integer :: b, c


  ! if (.not.IsActPtcl(a,vac)) stop 'PtclNum : index not act ptcl'
   if (.not.IsActPtcl(a,vac)) then
     write(UOUT,*) 'PctlNum error ! a=',a,'  vac=',vac, &
     ' and IsActPtcl(a,vac):',IsActPtcl(a,vac)
     call flushout
     call quit('PtclNum : index not act ptcl')
   endif

   c = 0
   do b = 1, a
      if (IsActPtcl(b,vac)) c = c + 1
   enddo
   PtclNum = c

   return

end function PtclNum

function VacNum(tvac)
!
!  returns the index of vacuum string for which 
! the string is equal with the entering value
!
   implicit none
   integer, intent(in) :: tvac
   integer :: VacNum
   integer :: vac

   do vac = 1, NSTATE
     if ( VacString(vac) == tvac ) exit
   enddo
   VacNum = vac

   return

end function VacNum

#include "maxmin.h"

!   These subroutines return the vaccum
!   number of a given excitation from
!   the given vaccum. The number is
!   extracted from VacTable by finding
!   the exc index of the given excitation.

function ExcVaccum1(i,a,vac)

   implicit none
   integer, intent(in) :: i, a, vac
   integer :: ExcVaccum1, exc

   if (.not.IsActHole(i,vac)) stop 'ExcVaccum : i not act hole'
   if (.not.IsActPtcl(a,vac)) stop 'ExcVaccum : a not act ptcl'

   exc = (HoleNum(i,vac)-1)*NPTCL + PtclNum(a,vac) + 1
   ExcVaccum1 = VacTable(exc,vac)

   return

end function ExcVaccum1

function ExcVaccum2(i,j,a,b,vac)

   implicit none
   integer, intent(in) :: i, j, a, b, vac
   integer :: ExcVaccum2, ij, ab, exc

  ! if (.not.IsActHole(i,vac)) stop 'ExcVaccum : i not act hole'
   if (.not.IsActHole(i,vac)) call quit('ExcVaccum2 : i not act hole')
  ! if (.not.IsActHole(j,vac)) stop 'ExcVaccum : j not act hole'
   if (.not.IsActHole(j,vac)) call quit('ExcVaccum2 : j not act hole')
  ! if (.not.IsActPtcl(a,vac)) stop 'ExcVaccum : a not act ptcl'
   if (.not.IsActPtcl(a,vac)) call quit('ExcVaccum2 : a not act ptcl')
  ! if (.not.IsActPtcl(b,vac)) stop 'ExcVaccum : b not act ptcl'
   if (.not.IsActPtcl(b,vac)) call quit('ExcVaccum2 : b not act ptcl')

   ij = MAXMIN(HoleNum(i,vac),HoleNum(j,vac))
   ab = MAXMIN(PtclNum(a,vac),PtclNum(b,vac))
   exc = (ij-1)*((NPTCL*(NPTCL-1))/2) + ab + NHOLE*NPTCL + 1
   ExcVaccum2 = VacTable(exc,vac)

   return

end function ExcVaccum2

function Exct1(i,a,vac)

   implicit none
   integer, intent(in) :: i, a, vac
   integer :: Exct1

  ! if (.not.IsActHole(i,vac)) stop 'ExcVaccum : i not act hole'
   if (.not.IsActHole(i,vac)) call quit('Exct1 : i not act hole')

  ! if (.not.IsActPtcl(a,vac)) stop 'ExcVaccum : a not act ptcl'
   if (.not.IsActPtcl(a,vac)) call quit('Exct1 : a not act ptcl')

   Exct1 = (HoleNum(i,vac)-1)*NPTCL + PtclNum(a,vac) + 1

   return

end function Exct1

function Exct2(i,j,a,b,vac)

   implicit none
   integer, intent(in) :: i, j, a, b, vac
   integer :: Exct2, ij, ab

  ! if (.not.IsActHole(i,vac)) stop 'ExcVaccum : i not act hole'
   if (.not.IsActHole(i,vac)) call quit('Exct2 : i not act hole')
  ! if (.not.IsActHole(j,vac)) stop 'ExcVaccum : j not act hole'
   if (.not.IsActHole(j,vac)) call quit('Exct2 : j not act hole')
  !if (.not.IsActPtcl(a,vac)) stop 'ExcVaccum : a not act ptcl'
   if (.not.IsActPtcl(a,vac)) call quit('Exct2 : a not act ptcl')
  !if (.not.IsActPtcl(b,vac)) stop 'ExcVaccum : b not act ptcl'
   if (.not.IsActPtcl(b,vac)) call quit('Exct2 : b not act ptcl')

   ij = MAXMIN(HoleNum(i,vac),HoleNum(j,vac))
   ab = MAXMIN(PtclNum(a,vac),PtclNum(b,vac))
   Exct2 = (ij-1)*((NPTCL*(NPTCL-1))/2) + ab + NHOLE*NPTCL + 1

   return

end function Exct2

!------------------------------------------------------------------------------
! This program is upto quadruple (double) excitations
! within model space. Higher excitations, if required,
! can be introduced in a similar fashion.
!------------------------------------------------------------------------------

subroutine GenVac
!
!    called from InitVacMod
!
   implicit none
   integer :: i, j, k, l, a, b, c, d
   integer :: vac, tvac, exc
   integer :: CurrVac, RefVac
   integer :: iocc, irp, ims2, ms2, ms2max, ovac
   logical :: flag
   integer :: nvac, ivac, nvac_sum
   integer :: isym_block

!-----------------------------------------------------------------------------------------
!
!  RefVac = 2**(NHOLE+1) - 2
!
!  RefVac is read from OccRef (which at the moment is read from gamess/dirac OccGMSRef)
!
!------------------------------------------------------------------------------------------

   if (get_iprint() >= 10) then
      write(UOUT,'(2x,a,i4,2x,a,i4)') 'GenVac..NCMO=',NCMO,'NVMO=',NVMO
      write(UOUT,'(2x,a,i4,a,/,20(1x,i2))') '-- OccRef(1-',NSO,'):',(OccRef(i),i=1,NSO)
   endif

   RefVac = 0
   do i = 1, NVMO
      if (OccRef(NCMO+i).eq.0) cycle
      if (OccRef(NCMO+i).eq.1) then
         RefVac = IBSET(RefVac,i)
      else
         write(UOUT,*) 'GenVac: incorrect occupation in OccRef array...see array OccRef'
         write(UOUT,*) '* OccRef: ',OccRef
         call flushout
         call Quit('GenVac: Incorrect Occupation in OccRef array !')
      endif
   enddo

   if (get_iprint() >= 10) then
      write(UOUT,'(2x,a,i4)') 'GenVac...set (by using NCMO,NVMO,OccRef) RefVac=',RefVac
   endif

   ! Generating other vaccums
   ! and storing them.

   VacString(1) = RefVac
   vac = 1

   ! Single Excitations

   do i = 1, NHOLE
   do a = 1, NPTCL
      vac = vac + 1
      CurrVac = RefVac
      CurrVac = IBCLR(CurrVac,i)
      CurrVac = IBSET(CurrVac,a+NHOLE)
      VacString(vac) = CurrVac

      if (get_iprint() >= 10) then
        write(UOUT,'(2x,a,2i5,a,i5,a,i8)')   & 
        'Genvac - single ext; i,a:',i,a,' VacString(',vac,')=',VacString(vac)
      endif

   enddo
   enddo

   ! Double Excitations

   do i = 2, NHOLE
   do j = 1, (i-1)
   do a = 2, NPTCL
   do b = 1, (a-1)
      vac = vac + 1
      CurrVac = RefVac
      CurrVac = IBCLR(CurrVac,i)
      CurrVac = IBCLR(CurrVac,j)
      CurrVac = IBSET(CurrVac,a+NHOLE)
      CurrVac = IBSET(CurrVac,b+NHOLE)
      VacString(vac) = CurrVac
      if ( get_iprint() >= 10 ) then
       write(UOUT,'(2x,a,4i5,a,i5,a,i8)')   & 
       'GenVac - double excit.; i,j,a,b:',i,j,a,b,' VacString(',vac,')=',VacString(vac)
      endif
   enddo
   enddo
   enddo
   enddo

   ! Triple Excitations

   do i = 2, NHOLE
   do j = 1, (i-1)
   do k = 1, (j-1)
   do a = 2, NPTCL
   do b = 1, (a-1)
   do c = 1, (b-1)
      vac = vac + 1
      CurrVac = RefVac
      CurrVac = IBCLR(CurrVac,i)
      CurrVac = IBCLR(CurrVac,j)
      CurrVac = IBCLR(CurrVac,k)
      CurrVac = IBSET(CurrVac,a+NHOLE)
      CurrVac = IBSET(CurrVac,b+NHOLE)
      CurrVac = IBSET(CurrVac,c+NHOLE)
      VacString(vac) = CurrVac
      if ( get_iprint() >= 10 ) then
        write(UOUT,'(2x,a,6i5,a,i5,a,i8)')   & 
        'GenVac - triple excit.; i,j,k,a,b,c:',i,j,k,a,b,c,' VacString(',vac,')=',VacString(vac)
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo

   ! Quadruple Excitation

   do i = 2, NHOLE
   do j = 1, (i-1)
   do k = 1, (j-1)
   do l = 1, (k-1)
   do a = 2, NPTCL
   do b = 1, (a-1)
   do c = 1, (b-1)
   do d = 1, (c-1)
      vac = vac + 1
      CurrVac = RefVac
      CurrVac = IBCLR(CurrVac,i)
      CurrVac = IBCLR(CurrVac,j)
      CurrVac = IBCLR(CurrVac,k)
      CurrVac = IBCLR(CurrVac,l)
      CurrVac = IBSET(CurrVac,a+NHOLE)
      CurrVac = IBSET(CurrVac,b+NHOLE)
      CurrVac = IBSET(CurrVac,c+NHOLE)
      CurrVac = IBSET(CurrVac,d+NHOLE)
      VacString(vac) = CurrVac
      if ( get_iprint() >= 10 ) then
       write(UOUT,'(2x,a,8i5,a,i5,a,i8)')   & 
       'GenVac - quadruple excit.; i,j,k,l,a,b,c,d:',i,j,k,l,a,b,c,d,' VacString(',vac,')=',VacString(vac)
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo

   ! -------------------------------------------------------
   !       Higher Excitations not yet implemented
   ! -------------------------------------------------------

   ! Generation of vaccum table to be used
   ! in renormalization diagrams.

   do vac = 1, NSTATE

      exc = 1
      VacTable(exc,vac) = vac

      ! single excitations

      do i = 1, NVMO
         if (.not.IsActHole(i,vac)) cycle
         do a = 1, NVMO
            if (.not.IsActPtcl(a,vac)) cycle
            exc = exc + 1
            tvac = VacString(vac)
            tvac = IBCLR(tvac,i)
            tvac = IBSET(tvac,a)
            ! VacNum function returns the vaccum
            ! number of Vaccum String tvac
            VacTable(exc,vac) = VacNum(tvac)
         enddo
      enddo

      ! double excitations

      do i = 1, NVMO
         if (.not.IsActHole(i,vac)) cycle
         do j = 1, (i-1)
            if (.not.IsActHole(j,vac)) cycle
            do a = 1, NVMO
               if (.not.IsActPtcl(a,vac)) cycle
               do b = 1, (a-1)
                  if (.not.IsActPtcl(b,vac)) cycle
                  exc = exc + 1
                  tvac = VacString(vac)
                  tvac = IBCLR(tvac,i)
                  tvac = IBCLR(tvac,j)
                  tvac = IBSET(tvac,a)
                  tvac = IBSET(tvac,b)
                  ! VacNum function returns the vaccum
                  ! number of Vaccum String tvac
                  VacTable(exc,vac) = VacNum(tvac)
               enddo
            enddo
         enddo
      enddo

      ! triple excitations

      do i = 1, NVMO
         if (.not.IsActHole(i,vac)) cycle
         do j = 1, (i-1)
            if (.not.IsActHole(j,vac)) cycle
            do k = 1, (j-1)
               if (.not.IsActHole(k,vac)) cycle
               do a = 1, NVMO
                  if (.not.IsActPtcl(a,vac)) cycle
                  do b = 1, (a-1)
                     if (.not.IsActPtcl(b,vac)) cycle
                     do c = 1, (b-1)
                        if (.not.IsActPtcl(c,vac)) cycle
                        exc = exc + 1
                        tvac = VacString(vac)
                        tvac = IBCLR(tvac,i)
                        tvac = IBCLR(tvac,j)
                        tvac = IBCLR(tvac,k)
                        tvac = IBSET(tvac,a)
                        tvac = IBSET(tvac,b)
                        tvac = IBSET(tvac,c)
                        ! VacNum function returns the vaccum
                        ! number of Vaccum String tvac
                        VacTable(exc,vac) = VacNum(tvac)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      ! quadruple excitations

      do i = 1, NVMO
         if (.not.IsActHole(i,vac)) cycle
         do j = 1, (i-1)
            if (.not.IsActHole(j,vac)) cycle
            do k = 1, (j-1)
               if (.not.IsActHole(k,vac)) cycle
               do l = 1, (k-1)
                  if (.not.IsActHole(l,vac)) cycle
                  do a = 1, NVMO
                     if (.not.IsActPtcl(a,vac)) cycle
                     do b = 1, (a-1)
                        if (.not.IsActPtcl(b,vac)) cycle
                        do c = 1, (b-1)
                           if (.not.IsActPtcl(c,vac)) cycle
                           do d = 1, (c-1)
                              if (.not.IsActPtcl(d,vac)) cycle
                              exc = exc + 1
                              tvac = VacString(vac)
                              tvac = IBCLR(tvac,i)
                              tvac = IBCLR(tvac,j)
                              tvac = IBCLR(tvac,k)
                              tvac = IBCLR(tvac,l)
                              tvac = IBSET(tvac,a)
                              tvac = IBSET(tvac,b)
                              tvac = IBSET(tvac,c)
                              tvac = IBSET(tvac,d)
                              ! VacNum function returns the vaccum
                              ! number of Vaccum String tvac
                              VacTable(exc,vac) = VacNum(tvac)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      ! -------------------------------------------
      !    higher excitations not yet implemented
      ! -------------------------------------------

   enddo

   if ( get_iprint() >= 10 ) then
     write(UOUT,'(/,2x,a,i5)') 'GenVac: totally_symmetric_irrep=',totally_symmetric_irrep
     write(UOUT,*) '**** irrep_product table ****:'
     write(UOUT,*) irrep_product
     write(UOUT,*) 
   endif

   ! miro: I have to check whether this condition is sufficient for read in symmetry
   if ( totally_symmetric_irrep /= 1 .or. numirrep == 0) then
     write(UOUT,'(/,2x,a,i5)') 'GenVac: totally_symmetric_irrep=',totally_symmetric_irrep
     write(UOUT,*) '**** irrep_product table ****:'
     write(UOUT,*) irrep_product
     write(UOUT,*) 
     call flushout
     call quit('GenVac: symmetry data NOT ready !')
   endif

   ! vaccum table generation is complete
   ! ------------------------------------------------------
   !
   ! For each vacuum, create the following
   !
   !    OccVac => Occupancy
   !    IrpVac => Irrep
   !    Ms2Vac => 2*M_s (excess number of alpha orbitals)
   !
   !    Needs "irrep_product" !
   !

   ms2max = 0

   do vac = 1, NSTATE

      irp = totally_symmetric_irrep
      ms2 = 0
      flag = .true.
      do i = 1, NVMO
         if (IsActHole(i,vac)) then
            iocc = 1
            irp = irrep_product(irp,IrrepOrb(i+NCMO))
            ms2 = ms2 + SpinOfOrb(i+NCMO)
            !if (get_iprint() >= 10) then
            !endif
         else
            iocc = 0
         endif
         if (flag) then
            if ((IsIntmCore(i)).and.(iocc.eq.0)) flag = flag.and..false.
            if ((IsIntmVirt(i)).and.(iocc.eq.1)) flag = flag.and..false.
         endif
         OccVac(i,vac) = iocc
      enddo
      IrpVac(vac) = irp
      Ms2Vac(vac) = ms2
      IsMainVac(vac) = flag
      ms2max = MAX(ms2,ms2max)
   enddo ! of vac = 1, NSTATE

   ! a symmetry block is defined by [ms2, irrep]
   ! find number of symmetry blocks in given model space
   ! also setup OrderVac and RevOrderVac arrays which
   ! help in reordering vacua according to symmetry blocks

   if ( get_iprint() >= 10 ) then
     write(UOUT,'(/,4x,a,i3,1x,a,i3)') 'ms2max=',ms2max,'numirrep=',numirrep
     do i=1,NSTATE
       write(UOUT,'(a,i2,a,i2,1x,a,i2,a,i2,1x,a,i2,a,i2)')    & 
       'IrpVac(',i,')=',IrpVac(i),'Ms2Vac(',i,')=',Ms2Vac(i), &
       'IsMainVac(',i,')=',IsMainVac(vac)
     enddo
   endif

   isym_block = 0
   nvac_sum = 0
   ovac = 0

   do ims2 = 0, ms2max

      ms2 = +ims2
      do irp = 1, numirrep

         nvac = 0
         do vac = 1, NSTATE
            if ((IrpVac(vac).eq.irp).and.(Ms2Vac(vac).eq.ms2)) then
               nvac = nvac + 1
               ovac = ovac + 1
               OrderVac(vac) = ovac
               RevOrderVac(ovac) = vac
            endif
         enddo
         if (nvac.ne.0) isym_block = isym_block + 1
         if (nvac.ne.0) nvac_sum = nvac_sum + nvac

      enddo ! of irp = 1, numirrep

      if (ims2.eq.0) cycle

      ms2 = -ims2
      do irp = 1, numirrep

         nvac = 0
         do vac = 1, NSTATE
            if ((IrpVac(vac).eq.irp).and.(Ms2Vac(vac).eq.ms2)) then
               nvac = nvac + 1
               ovac = ovac + 1
               OrderVac(vac) = ovac
               RevOrderVac(ovac) = vac
            endif
         enddo

         if (nvac.ne.0) isym_block = isym_block + 1
         if (nvac.ne.0) nvac_sum = nvac_sum + nvac

      enddo  ! of irp = 1, numirrep

   enddo ! of ims2 = 0, ms2max

   NSYM_BLOCK = isym_block

   if ( get_iprint() >= 10 ) then
     write(UOUT,*) 
     do i=1,NSTATE
       write(UOUT,'(2x,a,i3,a,i3,2x,a,i3)')   & 
         'OrderVac(',i,')=',OrderVac(i),      & 
         'RevOrderVac[OrderVac(i)]=',RevOrderVac(OrderVac(i))
     enddo
   endif

   ! ... error catching
   if (nvac_sum.ne.NSTATE) then
     write(UOUT,'(//,2x,a)')  &
     "GenVac: nvac_sum <> NSTATE; something wrong in finding symmetry blocks !"
     write(UOUT,'(2x,a,i4,1x,a,i4,/)') 'nvac_sum=',nvac_sum,' NSTATE=',NSTATE
     call flushout
     call quit('GenVac: something wrong in finding symmetry blocks !')
   endif

   ! ... allocate symblock_array and fill its elements
   allocate(symblock_array(NSYM_BLOCK))

   isym_block = 0
   ovac = 0
   do ims2 = 0, ms2max

      ms2 = +ims2
      do irp = 1, numirrep

         nvac = 0
         do vac = 1, NSTATE
            if ((IrpVac(vac).eq.irp).and.(Ms2Vac(vac).eq.ms2)) nvac = nvac + 1
         enddo

         if (nvac.ne.0) then
            isym_block = isym_block + 1
            symblock_array(isym_block)%ms2  = ms2
            symblock_array(isym_block)%irp  = irp
            symblock_array(isym_block)%nvac = nvac
            allocate(symblock_array(isym_block)%vacs(nvac))
            allocate(symblock_array(isym_block)%occp(NVMO,nvac))
            allocate(symblock_array(isym_block)%s2mat(nvac,nvac))
            ivac = 0
            do vac = 1, NSTATE
               if ((IrpVac(vac).eq.irp).and.(Ms2Vac(vac).eq.ms2)) then
                  ivac = ivac + 1
                  symblock_array(isym_block)%vacs(ivac) = vac
                  symblock_array(isym_block)%occp(:,ivac) = OccVac(:,vac)
               endif
            enddo

            ! s2mat evaluation done here
            call hmatrix(nvac,NVMO,NVMOSQ,nvac,            &
                         symblock_array(isym_block)%s2mat, &
                         symblock_array(isym_block)%occp,  &
                         s2_1, s2_2)
         endif

      enddo ! of irp = 1, numirrep

      if (ims2.eq.0) cycle

      ms2 = -ims2
      do irp = 1, numirrep

         nvac = 0
         do vac = 1, NSTATE
            if ((IrpVac(vac).eq.irp).and.(Ms2Vac(vac).eq.ms2)) nvac = nvac + 1
         enddo

         if (nvac.ne.0) then
            isym_block = isym_block + 1
            symblock_array(isym_block)%ms2  = ms2
            symblock_array(isym_block)%irp  = irp
            symblock_array(isym_block)%nvac = nvac
            allocate(symblock_array(isym_block)%vacs(nvac))
            allocate(symblock_array(isym_block)%occp(NVMO,nvac))
            allocate(symblock_array(isym_block)%s2mat(nvac,nvac))
            ivac = 0
            do vac = 1, NSTATE
               if ((IrpVac(vac).eq.irp).and.(Ms2Vac(vac).eq.ms2)) then
                  ivac = ivac + 1
                  symblock_array(isym_block)%vacs(ivac) = vac
                  symblock_array(isym_block)%occp(:,ivac) = OccVac(:,vac)
               endif
            enddo
            ! s2mat evaluation done here
         endif

      enddo ! of irp = 1, numirrep

   enddo ! of ims2 = 0, ms2max

   if (isym_block /= NSYM_BLOCK) then
     write(UOUT,'(/,2x,a)')   &
     'GenVac: isym_block <> NSYM_BLOCK; something wrong in creating sym_blocks !'
     write(UOUT,*) 'isym_block=',isym_block,' NSYM_BLOCK=',NSYM_BLOCK
     call flushout
     call Quit('GenVac')
   endif

   return

end subroutine GenVac

end module VacMod
