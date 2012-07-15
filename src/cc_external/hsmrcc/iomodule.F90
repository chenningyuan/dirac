module IOmodule

   use Parameters
   use Options
   implicit none

   ! ... Initialization Status Flag ...
   logical, private :: Init = .FALSE.

   ! ... IO Unit Numbers ...

   integer, parameter :: NZVINTBUFSIZE = 2048

   integer, parameter :: UMR1EINT = 10
   integer, parameter :: UMR2EINT = 11
   integer, parameter :: USYMINF  = 12

   integer, parameter :: UVINF = 21
   integer, parameter :: UVINT = 22
   integer, parameter :: UVIND = 23
   integer, parameter :: UFHH  = 24
   integer, parameter :: UFHP  = 25
   integer, parameter :: UFPP  = 26
   integer, parameter :: UHEFF = 27
   integer, parameter :: UCOEF = 28
   integer, parameter :: UCNVR = 29

   !integer :: UOUT  = 51 ! ... moved to parameters as this module looks more general 

   integer :: UTOLD = 52
   integer :: UTNEW = 53
   integer :: UTAMP = 54

   ! ... File Names ...

   character(len=80) :: FileHSMR
   character(len=80) :: FileMR1EINT
   character(len=80) :: FileMR2EINT
   character(len=80) :: FileSYMINF
   character(len=80) :: FileVINF
   character(len=80) :: FileVINT
   character(len=80) :: FileVIND
   character(len=80) :: FileFHH
   character(len=80) :: FileFHP
   character(len=80) :: FileFPP
   character(len=80) :: FileTAMP
   character(len=80) :: FileTOLD
   character(len=80) :: FileTNEW
   character(len=80) :: FileHEFF
   character(len=80) :: FileCOEF
   character(len=80) :: FileCNVR

   ! ... File Attributes ...

   integer, dimension(NVFIL) :: NumRec, RecLen
   integer, dimension(NVFIL) :: StartRec, EndngRec
   integer, dimension(NVFIL) :: TotNZInt, NumBlock

   ! ... Array containing info on Non-Zero Integrals in a record ...
   ! ... NZIntRec is a generic pointer ...

   integer, dimension(:), pointer :: NZIntRec
   integer, allocatable, target :: NZIntRec_V2H2P(:)
   integer, allocatable, target :: NZIntRec_VPPPP(:)
   integer, allocatable, target :: NZIntRec_VHHHH(:)
   integer, allocatable, target :: NZIntRec_VHHPP(:)
   integer, allocatable, target :: NZIntRec_V2HPP(:)
   integer, allocatable, target :: NZIntRec_V2PHH(:)
   integer, allocatable, target :: NZIntRec_VPHHP(:)
   integer, allocatable, target :: NZIntRec_VPH2P(:)
   integer, allocatable, target :: NZIntRec_VPHPP(:)
   integer, allocatable, target :: NZIntRec_V2HHP(:)
   integer, allocatable, target :: NZIntRec_VHHHP(:)

   ! ... Integral Data Buffers and their pointers ...

   integer         , dimension(:), pointer :: IndBuf
   double precision, dimension(:), pointer :: IntBuf
   integer         , private, allocatable, target :: BufIndex(:)
   double precision, private, allocatable, target :: BufIntgrl(:)

   ! ... Data for Block Read Operations ...

   integer, private, allocatable, dimension(:) :: StBlkRec_VPPPP, StBlkRec_VPH2P, StBlkRec_VPHPP
   integer, private, allocatable, dimension(:) :: EnBlkRec_VPPPP, EnBlkRec_VPH2P, EnBlkRec_VPHPP
   integer, private, allocatable, dimension(:) :: BlkStInd_VPPPP, BlkStInd_VPH2P, BlkStInd_VPHPP
   integer, private, allocatable, dimension(:) :: BlkEnInd_VPPPP, BlkEnInd_VPH2P, BlkEnInd_VPHPP
   integer, private, allocatable, dimension(:) :: BlkSt_VPPPP   , BlkSt_VPH2P   , BlkSt_VPHPP
   integer, private, allocatable, dimension(:) :: BlkEn_VPPPP   , BlkEn_VPH2P   , BlkEn_VPHPP

   ! ... Block Read Data Setup Routines ...

   private :: SetupBlkRead_VPPPP
   private :: SetupBlkRead_VPH2P
   private :: SetupBlkRead_VPHPP
   private :: NumBlock_VPPPP
   private :: NumBlock_VPH2P
   private :: NumBlock_VPHPP

   ! ... Generic Interface for Integral FileIO ...

   interface ReadFile
      module procedure ReadFileFull
      module procedure ReadFileBlock
   end interface

   ! ... External Routines ... mi: these are replaced by dirac/own counterparts ...

!   external :: flush
!   external :: abort

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_dirac_hscc_file_names
!
! assign file names for 6 Dirac-prepared files 
!
! called from Dirac_HSMRCC interface 
!

! write(UOUT,*) 'set_dirac_hscc_file_names'
! write(UOUT,*) 'get_iprint=',get_iprint()

FileVINT='Dirac_HSMRCC.VINT'
FileVIND='Dirac_HSMRCC.VIND'
FileVINF='Dirac_HSMRCC.VINF'
FileFHH ='Dirac_HSMRCC.FHH'
FileFHP ='Dirac_HSMRCC.FHP'
FileFPP ='Dirac_HSMRCC.FPP'

FileTAMP = 'Dirac_HSMRCC.TAMP'
FileTOLD  ='Dirac_HSMRCC.TOLD'
FileTNEW  ='Dirac_HSMRCC.TNEW'
FileHEFF  ='Dirac_HSMRCC.HEFF'
FileCOEF  ='Dirac_HSMRCC.COEF'
FileCNVR  ='Dirac_HSMRCC.CNVR'


if (get_iprint() >= 5) then
 write(UOUT,'(/,4x,a,/)') '=== Output from "set_dirac_hscc_file_names" ==='
 write(UOUT,'(2x,a,i2,a,a20)') 'lu,UVINF=',UVINF,' dirac assigned  FileVINF=',FileVINF
 write(UOUT,'(2x,a,i2,a,a20)') 'lu,UVINT=',UVINT,' dirac assigned  FileVINT=',FileVINT
 write(UOUT,'(2x,a,i2,a,a20)') 'lu,UVIND=',UVIND,' dirac assigned  FileVIND=',FileVIND
 write(UOUT,'(2x,a,i2,a,a20)') 'lu, UFHH=',UFHH, ' dirac assigned   FileFHH=',FileFHH
 write(UOUT,'(2x,a,i2,a,a20)') 'lu, UFHP=',UFHP, ' dirac assigned   FileFHH=',FileFHP
 write(UOUT,'(2x,a,i2,a,a20)') 'lu, UFPP=',UFPP, ' dirac assigned   FileFHH=',FileFPP
 write(UOUT,*) 
endif

end subroutine set_dirac_hscc_file_names

subroutine reset_Init
! ... set to FALSE initialization variable (needed in DIRAC HSCC preparing files part)
 Init=.FALSE.
end subroutine reset_Init

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine InitIOMod
!---------------------------------------------------------------------------------
!
! calculate attributes NumRec & RecLen and allocate arrays NZIntRec_xxxx
! for reading Gamess/DIRAC-made integrals & for preparing DIRAC-HSMRCC integrals
!
!  called from : 
!         InitHSMRBase ("source/inihsbase.F90")
!         dirac HSCC interface ("dir_hscc_mods.F90")
!   
!---------------------------------------------------------------------------------

 implicit none

   ! ... Initialization Check ...

   if (Init) then
      write(UOUT,*) 'Error Initializing IOMod, Init=',Init
      call quit('InitIOMod:  Error Initializing IOMod ! It was already set to true !')
   endif
   Init = .TRUE. ! it means it's initialized

   ! first check NumRec/RecLen dimensions (ie nonzero NHMO, NPMO parameters)
   call check_FileRec_param

   ! ... Setup File Attributes ...
   NumRec(V2H2P) = NHMO*NPMO
   NumRec(VPPPP) = NPMOSQ
   NumRec(VHHHH) = NHMOSQ
   NumRec(VHHPP) = NHMOSQ
   NumRec(V2HPP) = NHMO*NHMO
   NumRec(V2PHH) = NPMO*NPMO
   NumRec(VPHHP) = NHMO*NPMO
   NumRec(VPH2P) = NPMO*NPMO
   NumRec(VPHPP) = NPMO*NHMO
   NumRec(V2HHP) = NHMO*NHMO
   NumRec(VHHHP) = NHMOSQ

   NumRec(FHH)   = NHMO
   NumRec(FHP)   = NHMO
   NumRec(FPP)   = NPMO

   RecLen(V2H2P) = NHMO*NPMO
   RecLen(VPPPP) = NPMOSQ
   RecLen(VHHHH) = NHMOSQ
   RecLen(VHHPP) = NPMOSQ
   RecLen(V2HPP) = NPMOSQ
   RecLen(V2PHH) = NHMOSQ
   RecLen(VPHHP) = NHMO*NPMO
   RecLen(VPH2P) = NHMO*NPMO
   RecLen(VPHPP) = NPMOSQ
   RecLen(V2HHP) = NHMO*NPMO
   RecLen(VHHHP) = NHMO*NPMO
   RecLen(FHH)   = NVMO*NHMO
   RecLen(FHP)   = NVMO*NPMO
   RecLen(FPP)   = NVMO*NPMO

   if (get_iprint() >= 5) then
     call print_int_file_attributes
   endif

   call InitIOMod_alloc_NZIntRec

   return

end subroutine InitIOMod

subroutine print_int_file_attributes

  write(UOUT,'(/,4x,a,/)') '*** NumRec & RecLen data ***'

  write(UOUT,9234) 'NumRec(V2H2P)=',NumRec(V2H2P),' RecLen(V2H2P)=',RecLen(V2H2P)
  write(UOUT,9234) 'NumRec(VPPPP)=',NumRec(VPPPP),' RecLen(V2H2P)=',RecLen(VPPPP)
  write(UOUT,9234) 'NumRec(VHHHH)=',NumRec(VHHHH),' RecLen(VHHHH)=',RecLen(VHHHH)
  write(UOUT,9234) 'NumRec(VHHPP)=',NumRec(VHHPP),' RecLen(VHHPP)=',RecLen(VHHPP)
  write(UOUT,9234) 'NumRec(V2HPP)=',NumRec(V2HPP),' RecLen(V2HPP)=',RecLen(V2HPP)
  write(UOUT,9234) 'NumRec(V2PHH)=',NumRec(V2PHH),' RecLen(V2PHH)=',RecLen(V2PHH)
  write(UOUT,9234) 'NumRec(VPHHP)=',NumRec(VPHHP),' RecLen(VPHHP)=',RecLen(VPHHP)
  write(UOUT,9234) 'NumRec(VPH2P)=',NumRec(VPH2P),' RecLen(VPH2P)=',RecLen(VPH2P)
  write(UOUT,9234) 'NumRec(VPHPP)=',NumRec(VPHPP),' RecLen(VPHPP)=',RecLen(VPHPP)
  write(UOUT,9234) 'NumRec(V2HHP)=',NumRec(V2HHP),' RecLen(V2HHP)=',RecLen(V2HHP)
  write(UOUT,9234) 'NumRec(VHHHP)=',NumRec(VHHHP),' RecLen(VHHHP)=',RecLen(VHHHP)
  write(UOUT,*) 
  write(UOUT,9235) 'NumRec(FHH)=',NumRec(FHH),' RecLen(FHH)=',RecLen(FHH)
  write(UOUT,9235) 'NumRec(FHP)=',NumRec(FHP),' RecLen(FHP)=',RecLen(FHP)
  write(UOUT,9235) 'NumRec(FPP)=',NumRec(FPP),' RecLen(FPP)=',RecLen(FPP)

   9234 format(2x,a,i5,a,i5)
   9235 format(2x,a,i4,a,i4)

end subroutine print_int_file_attributes 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine InitIOMod_alloc_NZIntRec
!-------------------------------------------------------------------------------------
!
! allocation of NZIntRec in separate routine
!
! called from InitIOMod
!
!-------------------------------------------------------------------------------------

   ! ... Allocate & Setup NZIntRec Arrays ...
   if (.not.allocated(NZIntRec_V2H2P)) allocate(NZIntRec_V2H2P(NumRec(V2H2P)))
   if (.not.allocated(NZIntRec_VPPPP)) allocate(NZIntRec_VPPPP(NumRec(VPPPP)))
   if (.not.allocated(NZIntRec_VHHHH)) allocate(NZIntRec_VHHHH(NumRec(VHHHH)))
   if (.not.allocated(NZIntRec_VHHPP)) allocate(NZIntRec_VHHPP(NumRec(VHHPP)))
   if (.not.allocated(NZIntRec_V2HPP)) allocate(NZIntRec_V2HPP(NumRec(V2HPP)))
   if (.not.allocated(NZIntRec_V2PHH)) allocate(NZIntRec_V2PHH(NumRec(V2PHH)))
   if (.not.allocated(NZIntRec_VPHHP)) allocate(NZIntRec_VPHHP(NumRec(VPHHP)))
   if (.not.allocated(NZIntRec_VPH2P)) allocate(NZIntRec_VPH2P(NumRec(VPH2P)))
   if (.not.allocated(NZIntRec_VPHPP)) allocate(NZIntRec_VPHPP(NumRec(VPHPP)))
   if (.not.allocated(NZIntRec_V2HHP)) allocate(NZIntRec_V2HHP(NumRec(V2HHP)))
   if (.not.allocated(NZIntRec_VHHHP)) allocate(NZIntRec_VHHHP(NumRec(VHHHP)))

   if (get_iprint() > 5) then
     write(UOUT,'(/,4x,a)') '*** Output from InitIOMod_alloc_NZIntRec ***'
     write(UOUT,'(2x,a,/)') '11 NZIntRec_Vxxxx arrays allocated'
   endif

end subroutine InitIOMod_alloc_NZIntRec

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine ReadNZIntRecs
!-------------------------------------------------------------------------------------------------------
!
!
!  Read data on records on FileVINF file: StartRec,EndngRec,TotNZInt, NZIntRec_Vxxxx for 11 int.types
!
!
!-------------------------------------------------------------------------------------------------------
   implicit none

   if (.not.Init) then
      write(UOUT,*) 'ReadNZIntRecs: IOMod not initialized'
      call flushout
      call quit('ReadNZIntRecs: IOMod not initialized!')
   endif

   rewind(UVINF)

   read(UVINF)
   read(UVINF)
   read(UVINF)
   read(UVINF)
   read(UVINF)
   read(UVINF)
   read(UVINF)                
   read(UVINF)                
   read(UVINF)                

   read(UVINF) StartRec(V2H2P), EndngRec(V2H2P), TotNZInt(V2H2P), NZIntRec_V2H2P
   read(UVINF) StartRec(VPPPP), EndngRec(VPPPP), TotNZInt(VPPPP), NZIntRec_VPPPP
   read(UVINF) StartRec(VHHHH), EndngRec(VHHHH), TotNZInt(VHHHH), NZIntRec_VHHHH
   read(UVINF) StartRec(VHHPP), EndngRec(VHHPP), TotNZInt(VHHPP), NZIntRec_VHHPP
   read(UVINF) StartRec(V2HPP), EndngRec(V2HPP), TotNZInt(V2HPP), NZIntRec_V2HPP
   read(UVINF) StartRec(V2PHH), EndngRec(V2PHH), TotNZInt(V2PHH), NZIntRec_V2PHH
   read(UVINF) StartRec(VPHHP), EndngRec(VPHHP), TotNZInt(VPHHP), NZIntRec_VPHHP
   read(UVINF) StartRec(VPH2P), EndngRec(VPH2P), TotNZInt(VPH2P), NZIntRec_VPH2P
   read(UVINF) StartRec(VPHPP), EndngRec(VPHPP), TotNZInt(VPHPP), NZIntRec_VPHPP
   read(UVINF) StartRec(V2HHP), EndngRec(V2HHP), TotNZInt(V2HHP), NZIntRec_V2HHP
   read(UVINF) StartRec(VHHHP), EndngRec(VHHHP), TotNZInt(VHHHP), NZIntRec_VHHHP
   read(UVINF) StartRec(FHH)  , EndngRec(FHH)  , TotNZInt(FHH)
   read(UVINF) StartRec(FHP)  , EndngRec(FHP)  , TotNZInt(FHP)
   read(UVINF) StartRec(FPP)  , EndngRec(FPP)  , TotNZInt(FPP)

   if (get_iprint() > 5) then
     write(UOUT,'(/,4x,a)') '*** Output from ReadNZIntRecs ***'
     write(UOUT,'(2x,a)') 'StartRec/EndngRec/TotNZInt/NZIntRec_Vxxxx record arrays read from FileVINF'
   endif

   return
end subroutine ReadNZIntRecs

subroutine InitBlockIO

   implicit none

   if (.not.Init) then
      write(UOUT,*) 'InitBlockIO: IOMod not initialized'
      call flushout
      call quit('InitBlockIO: IOMod not initialized')
   endif

   ! ... Setup Data for Block IO   ...
   ! ... Block IO is available     ...
   ! ... only for following files  ...
   ! ...   VPPPP, VPH2P, VPHPP     ...

   NumBlock(V2H2P) = 1
   NumBlock(VHHHH) = 1
   NumBlock(VHHPP) = 1
   NumBlock(V2HPP) = 1
   NumBlock(V2PHH) = 1
   NumBlock(VPHHP) = 1
   NumBlock(V2HHP) = 1
   NumBlock(VHHHP) = 1
   NumBlock(FHH)   = 1

   ! ... block reading procedure ...
   NumBlock(VPPPP) = NumBlock_VPPPP()
   NumBlock(VPH2P) = NumBlock_VPH2P()
   NumBlock(VPHPP) = NumBlock_VPHPP()

   if (.not.allocated(StBlkRec_VPPPP)) allocate(StBlkRec_VPPPP(NumBlock(VPPPP)))
   if (.not.allocated(EnBlkRec_VPPPP)) allocate(EnBlkRec_VPPPP(NumBlock(VPPPP)))
   if (.not.allocated(BlkStInd_VPPPP)) allocate(BlkStInd_VPPPP(NumBlock(VPPPP)))
   if (.not.allocated(BlkEnInd_VPPPP)) allocate(BlkEnInd_VPPPP(NumBlock(VPPPP)))
   if (.not.allocated(BlkSt_VPPPP))    allocate(BlkSt_VPPPP(NumBlock(VPPPP)))
   if (.not.allocated(BlkEn_VPPPP))    allocate(BlkEn_VPPPP(NumBlock(VPPPP)))

   if (.not.allocated(StBlkRec_VPH2P)) allocate(StBlkRec_VPH2P(NumBlock(VPH2P)))
   if (.not.allocated(EnBlkRec_VPH2P)) allocate(EnBlkRec_VPH2P(NumBlock(VPH2P)))
   if (.not.allocated(BlkStInd_VPH2P)) allocate(BlkStInd_VPH2P(NumBlock(VPH2P)))
   if (.not.allocated(BlkEnInd_VPH2P)) allocate(BlkEnInd_VPH2P(NumBlock(VPH2P)))
   if (.not.allocated(BlkSt_VPH2P))    allocate(BlkSt_VPH2P(NumBlock(VPH2P)))
   if (.not.allocated(BlkEn_VPH2P))    allocate(BlkEn_VPH2P(NumBlock(VPH2P)))

   if (.not.allocated(StBlkRec_VPHPP)) allocate(StBlkRec_VPHPP(NumBlock(VPHPP)))
   if (.not.allocated(EnBlkRec_VPHPP)) allocate(EnBlkRec_VPHPP(NumBlock(VPHPP)))
   if (.not.allocated(BlkStInd_VPHPP)) allocate(BlkStInd_VPHPP(NumBlock(VPHPP)))
   if (.not.allocated(BlkEnInd_VPHPP)) allocate(BlkEnInd_VPHPP(NumBlock(VPHPP)))
   if (.not.allocated(BlkSt_VPHPP))    allocate(BlkSt_VPHPP(NumBlock(VPHPP)))
   if (.not.allocated(BlkEn_VPHPP))    allocate(BlkEn_VPHPP(NumBlock(VPHPP)))

   call SetupBlkRead_VPPPP
   call SetupBlkRead_VPH2P
   call SetupBlkRead_VPHPP

   return

end subroutine InitBlockIO

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine AllocIntBufs
!-------------------------------------------------------------------
!
!                Allocate integrals/indexes buffers
!
!-------------------------------------------------------------------
   implicit none
   integer :: IntBufSize

   if (.not.Init) then
      write(UOUT,*) 'AllocIntBufs: IOMod not initialized'
      call flushout
      call quit('AllocIntBufs : IOMod not initialized!')
   endif

   ! ... Determine IntBufSize ...

   IntBufSize = MAX(1,(EndngRec(VHHHH)-StartRec(VHHHH)+1))
   IntBufSize = MAX(IntBufSize,(EndngRec(VHHPP)-StartRec(VHHPP)+1))
   IntBufSize = MAX(IntBufSize,(EndngRec(V2H2P)-StartRec(V2H2P)+1))
   IntBufSize = MAX(IntBufSize,(EndngRec(V2HPP)-StartRec(V2HPP)+1))
   IntBufSize = MAX(IntBufSize,(EndngRec(V2PHH)-StartRec(V2PHH)+1))
   IntBufSize = MAX(IntBufSize,(EndngRec(VPHHP)-StartRec(VPHHP)+1))
   IntBufSize = MAX(IntBufSize,(EndngRec(V2HHP)-StartRec(V2HHP)+1))
   IntBufSize = MAX(IntBufSize,(EndngRec(VHHHP)-StartRec(VHHHP)+1))
   IntBufSize = IntBufSize*NZVINTBUFSIZE

   ! ... Setup Integral Data Buffers ...

   Nullify(IndBuf)
   Nullify(IntBuf)
   if (.not.allocated(BufIndex)) then
      allocate(BufIndex(IntBufSize))
      if (get_iprint() >= 12) then
        write(UOUT,'(/,2x,a,i8)')  & 
        'AllocIntBufs: indexes (integer) BufIndex allocated, size=',IntBufSize
      endif
   endif

   if (.not.allocated(BufIntgrl)) then 
      allocate(BufIntgrl(IntBufSize))
      if (get_iprint() >= 12) then
        write(UOUT,'(2x,a,i8)')  & 
        'AllocIntBufs: integrals (double) BufIntgrl allocated, size=',IntBufSize
      endif
   endif

   return

end subroutine AllocIntBufs

subroutine DeAllocIntBufs

   implicit none

   if (.not.Init) then
      write(UOUT,*) 'DeAllocIntBufs: IOMod not initialized'
      call flushout
      call quit('DeAllocIntBufs: IOMod not initialized!')
   endif

   nullify(IndBuf)
   nullify(IntBuf)
   if (allocated(BufIndex))  deallocate(BufIndex)
   if (allocated(BufIntgrl)) deallocate(BufIntgrl)

   return

end subroutine DeAllocIntBufs

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine check_VIND_VINT_readiness
!  ... check whether files VIND,VINT are open,
!       and arrays BufIntgrl,BufIndex allocated
logical :: is_VIND_file, is_VIND_open, is_VINT_file, is_VINT_open

inquire(file=FileVIND,exist=is_VIND_file, opened=is_VIND_open)
inquire(file=FileVINT,exist=is_VINT_file, opened=is_VINT_open)

if (.not.(allocated(BufIntgrl).and.allocated(BufIndex).and.is_VIND_open.and.is_VINT_open)) then
  write(UOUT,*) 'is_VIND_file=',is_VIND_file
  write(UOUT,*) 'is_VINT_file=',is_VINT_file
  write(UOUT,*) 'is_VIND_open=',is_VIND_open
  write(UOUT,*) 'is_VINT_open=',is_VINT_open
  write(UOUT,*) 'allocated(BufIntgrl):',allocated(BufIntgrl)
  write(UOUT,*) ' allocated(BufIndex):',allocated(BufIndex)
  call quit('check_VIND_VINT_readiness: VINT/VIND not ready for reading !')
endif

end subroutine check_VIND_VINT_readiness

function NumBlock_VPPPP()

   implicit none
   integer :: NumBlock_VPPPP
   integer :: NumFileRec, BufSize, MaxBufRec
   integer :: IntRec, NRead, NunRead, NZInt
   integer :: NBlk

   NZIntRec => NZIntRec_VPPPP
   NumFileRec = EndngRec(VPPPP) - StartRec(VPPPP) + 1
   BufSize = MIN(SIZE(BufIntgrl),SIZE(BufIndex))
   MaxBufRec  = MIN(BufSize/NZVINTBUFSIZE,NumFileRec) 
   if ( MaxBufRec < 1 ) then
      write(UOUT,*) 'Integral Buffer too small for VPPPP'
      call flushout
      call quit('NumBlock_VPPPP: Integral Buffer too small for VPPPP !')
   endif
   ! BufSize is reset to effective buffer size
   BufSize = MaxBufRec*NZVINTBUFSIZE

   IntRec = 1
   NRead = 0
   NBlk = 0
   do
      NBlk = NBlk + 1
      NZInt = 0
      do
         if ( IntRec > NumRec(VPPPP) ) exit
         if ( (NZInt + NZIntRec(IntRec)) > (BufSize-NRead) ) exit
         NZInt = NZInt + NZIntRec(IntRec)
         IntRec = IntRec + 1
      enddo
      if (NZInt == 0) then 
      ! ... more debug info...
        if ( IntRec > NumRec(VPPPP) ) then
          write(UOUT,*) 'IntRec=',IntRec,' NumRec(VPPPP)=',NumRec(VPPPP)
          write(UOUT,*) 'exit due to IntRec > NumRec(VPPPP) !'
        endif
        if ( (NZInt + NZIntRec(IntRec)) > (BufSize-NRead) ) then
          write(UOUT,*) 'exit due to  (NZInt + NZIntRec(IntRec)) > (BufSize-NRead) !'
          write(UOUT,*) ' (NZInt + NZIntRec(IntRec)):',NZInt + NZIntRec(IntRec)
          write(UOUT,*) ' (BufSize-NRead)=',(BufSize-NRead) 
        endif
       ! write(UOUT,*) 'NBlk=',NBlk
       ! write(UOUT,*) 'NZIntRec(IntRec)=',NZIntRec(IntRec)
        call Quit('NumBlock_VPPPP: Error NZInt=0 !')
      endif
      NunRead = BufSize - NRead - NZInt
      NRead = NZVINTBUFSIZE - MOD(NunRead,NZVINTBUFSIZE)
      if (NRead == NZVINTBUFSIZE) NRead = 0
      if (IntRec > NumRec(VPPPP)) exit
   enddo
   NumBlock_VPPPP = NBlk

   return

end function NumBlock_VPPPP

function NumBlock_VPH2P()

   implicit none
   integer :: NumBlock_VPH2P
   integer :: NumFileRec, BufSize, MaxBufRec
   integer :: IntRec, NRead, NunRead, NZInt
   integer :: NBlk, st, en, ni

   NZIntRec => NZIntRec_VPH2P
   NumFileRec = EndngRec(VPH2P) - StartRec(VPH2P) + 1
   BufSize = MIN(SIZE(BufIntgrl),SIZE(BufIndex))
   MaxBufRec  = MIN(BufSize/NZVINTBUFSIZE,NumFileRec) 
   if ( MaxBufRec < 1 ) then
      write(UOUT,*) 'NumBlock_VPH2P: Integral Buffer too small for VPH2P'
      call flushout
      call quit('NumBlock_VPH2P: Integral Buffer too small for VPH2P!')
   endif

   ! BufSize is reset to effective buffer size
   BufSize = MaxBufRec*NZVINTBUFSIZE

   IntRec = 1
   NRead = 0
   NBlk = 0
   do
      NBlk = NBlk + 1
      NZInt = 0
      do
         if ( IntRec > NPMO ) exit
         st = (IntRec-1)*NPMO+1
         en = st + NPMO - 1
         ni = SUM(NZIntRec(st:en))
         if ((NZInt + ni) > (BufSize-NRead)) exit
         NZInt = NZInt + ni
         IntRec = IntRec + 1
      enddo
      if (NZInt == 0) then
        call Quit('NumBlock_VPH2P: Error NZInt=0')
      endif
      NunRead = BufSize - NRead - NZInt
      NRead = NZVINTBUFSIZE - MOD(NunRead,NZVINTBUFSIZE)
      if (NRead == NZVINTBUFSIZE) NRead = 0
      if (IntRec > NPMO) exit
   enddo
   NumBlock_VPH2P = NBlk

   return

end function NumBlock_VPH2P

function NumBlock_VPHPP()

   implicit none
   integer :: NumBlock_VPHPP
   integer :: NumFileRec, BufSize, MaxBufRec
   integer :: IntRec, NRead, NunRead, NZInt
   integer :: NBlk, st, en, ni

   NZIntRec => NZIntRec_VPHPP
   NumFileRec = EndngRec(VPHPP) - StartRec(VPHPP) + 1
   BufSize = MIN(SIZE(BufIntgrl),SIZE(BufIndex))
   MaxBufRec  = MIN(BufSize/NZVINTBUFSIZE,NumFileRec) 
   if ( MaxBufRec < 1 ) then
      write(UOUT,*) 'NumBlock_VPHPP: Integral Buffer too small for VPHPP'
      call flushout
      call quit('NumBlock_VPHPP: Integral Buffer too small for VPHPP')
   endif
   ! BufSize is reset to effective buffer size
   BufSize = MaxBufRec*NZVINTBUFSIZE

   IntRec = 1
   NRead = 0
   NBlk = 0
   do
      NBlk = NBlk + 1
      NZInt = 0
      do
         if ( IntRec > NPMO ) exit
         st = (IntRec-1)*NHMO+1
         en = st + NHMO - 1
         ni = SUM(NZIntRec(st:en))
         if ((NZInt + ni) > (BufSize-NRead)) exit
         NZInt = NZInt + ni
         IntRec = IntRec + 1
      enddo
      if (NZInt == 0) then 
        call quit('NumBlock_VPHPP: Error NZInt=0')
      endif
      NunRead = BufSize - NRead - NZInt
      NRead = NZVINTBUFSIZE - MOD(NunRead,NZVINTBUFSIZE)
      if (NRead == NZVINTBUFSIZE) NRead = 0
      if (IntRec > NPMO) exit
   enddo
   NumBlock_VPHPP = NBlk

   return

end function NumBlock_VPHPP

subroutine SetupBlkRead_VPPPP

   implicit none
   integer :: NumBlock_VPPPP
   integer :: NumFileRec, BufSize, MaxBufRec
   integer :: IntRec, NRead, NunRead, NZInt
   integer :: SRec, StR, EnR, StRec, EnRec, NBRec
   integer :: NBlk

   NZIntRec => NZIntRec_VPPPP
   NumFileRec = EndngRec(VPPPP) - StartRec(VPPPP) + 1
   BufSize = MIN(SIZE(BufIntgrl),SIZE(BufIndex))
   MaxBufRec  = MIN(BufSize/NZVINTBUFSIZE,NumFileRec) 
   if ( MaxBufRec < 1 ) then
      write(UOUT,*) 'SetupBlkRead_VPPPP: Integral Buffer too small for VPPPP'
      call flushout
      call quit('SetupBlkRead_VPPPP: Integral Buffer too small for VPPPP')
   endif
   ! BufSize is reset to effective buffer size
   BufSize = MaxBufRec*NZVINTBUFSIZE

   IntRec = 1
   NRead = 0
   SRec = StartRec(VPPPP)
   NBlk = 0
   do
      NBlk = NBlk + 1
      StBlkRec_VPPPP(NBlk) = SRec
      StRec = IntRec
      NZInt = 0
      do
         if ( IntRec > NumRec(VPPPP) ) exit
         if ( (NZInt + NZIntRec(IntRec)) > (BufSize-NRead) ) exit
         NZInt = NZInt + NZIntRec(IntRec)
         IntRec = IntRec + 1
      enddo
      if ( NZInt == 0 ) then
         call quit('SetupBlkRead_VPPPP: Error NZInt=0')
      endif
      EnRec = IntRec - 1
      NBRec = (NRead+NZInt)/NZVINTBUFSIZE
      if (NBRec*NZVINTBUFSIZE/=(NRead+NZInt)) NBRec = NBRec + 1
      EnBlkRec_VPPPP(NBlk) = SRec + NBRec - 1
      BlkStInd_VPPPP(NBlk) = NRead + 1
      BlkEnInd_VPPPP(NBlk) = NRead + NZInt
      BlkSt_VPPPP(NBlk) = StRec
      BlkEn_VPPPP(NBlk) = EnRec
      SRec = SRec + NBRec
      NunRead = BufSize - NRead - NZInt
      NRead = NZVINTBUFSIZE - MOD(NunRead,NZVINTBUFSIZE)
      if ( NRead == NZVINTBUFSIZE ) NRead = 0
      if ( NRead /= 0 ) SRec = SRec - 1
      if ( IntRec > NumRec(VPPPP)) exit
    enddo

end subroutine SetupBlkRead_VPPPP

subroutine SetupBlkRead_VPH2P

   implicit none
   integer :: NumBlock_VPH2P
   integer :: NumFileRec, BufSize, MaxBufRec
   integer :: IntRec, NRead, NunRead, NZInt
   integer :: SRec, StR, EnR, StRec, EnRec, NBRec
   integer :: NBlk, st, en, ni

   NZIntRec => NZIntRec_VPH2P
   NumFileRec = EndngRec(VPH2P) - StartRec(VPH2P) + 1
   BufSize = MIN(SIZE(BufIntgrl),SIZE(BufIndex))
   MaxBufRec  = MIN(BufSize/NZVINTBUFSIZE,NumFileRec) 
   if ( MaxBufRec < 1 ) then
      write(UOUT,*) 'SetupBlkRead_VPH2P: Integral Buffer too small for VPH2P'
      call flushout
      call quit('SetupBlkRead_VPH2P: Integral Buffer too small for VPH2P')
   endif
   ! BufSize is reset to effective buffer size
   BufSize = MaxBufRec*NZVINTBUFSIZE

   IntRec = 1
   NRead = 0
   SRec = StartRec(VPH2P)
   NBlk = 0
   do
      NBlk = NBlk + 1
      StBlkRec_VPH2P(NBlk) = SRec
      StRec = IntRec
      NZInt = 0
      do
         if ( IntRec > NPMO ) exit
         st = (IntRec-1)*NPMO+1
         en = st + NPMO - 1
         ni = SUM(NZIntRec(st:en))
         if ((NZInt + ni) > (BufSize-NRead)) exit
         NZInt = NZInt + ni
         IntRec = IntRec + 1
      enddo
      !if ( NZInt == 0 ) stop 'Error NZInt=0'
      if ( NZInt == 0 ) then 
        call quit('SetupBlkRead_VPH2P: Error NZInt=0')
      endif
      EnRec = IntRec - 1
      NBRec = (NRead+NZInt)/NZVINTBUFSIZE
      if (NBRec*NZVINTBUFSIZE/=(NRead+NZInt)) NBRec = NBRec + 1
      EnBlkRec_VPH2P(NBlk) = SRec + NBRec - 1
      BlkStInd_VPH2P(NBlk) = NRead + 1
      BlkEnInd_VPH2P(NBlk) = NRead + NZInt
      BlkSt_VPH2P(NBlk) = StRec
      BlkEn_VPH2P(NBlk) = EnRec
      SRec = SRec + NBRec
      NunRead = BufSize - NRead - NZInt
      NRead = NZVINTBUFSIZE - MOD(NunRead,NZVINTBUFSIZE)
      if ( NRead == NZVINTBUFSIZE ) NRead = 0
      if ( NRead /= 0 ) SRec = SRec - 1
      if ( IntRec > NPMO) exit
    enddo

end subroutine SetupBlkRead_VPH2P

subroutine SetupBlkRead_VPHPP

   implicit none
   integer :: NumBlock_VPHPP
   integer :: NumFileRec, BufSize, MaxBufRec
   integer :: IntRec, NRead, NunRead, NZInt
   integer :: SRec, StR, EnR, StRec, EnRec, NBRec
   integer :: NBlk, st, en, ni

   NZIntRec => NZIntRec_VPHPP
   NumFileRec = EndngRec(VPHPP) - StartRec(VPHPP) + 1
   BufSize = MIN(SIZE(BufIntgrl),SIZE(BufIndex))
   MaxBufRec  = MIN(BufSize/NZVINTBUFSIZE,NumFileRec) 
   if ( MaxBufRec < 1 ) then
      write(UOUT,*) 'SetupBlkRead_VPHPP:Integral Buffer too small for VPHPP'
      call flushout
      call quit('SetupBlkRead_VPHPP:Integral Buffer too small for VPHPP')
   endif
   ! BufSize is reset to effective buffer size
   BufSize = MaxBufRec*NZVINTBUFSIZE

   IntRec = 1
   NRead = 0
   SRec = StartRec(VPHPP)
   NBlk = 0
   do
      NBlk = NBlk + 1
      StBlkRec_VPHPP(NBlk) = SRec
      StRec = IntRec
      NZInt = 0
      do
         if ( IntRec > NPMO ) exit
         st = (IntRec-1)*NHMO+1
         en = st + NHMO - 1
         ni = SUM(NZIntRec(st:en))
         if ((NZInt + ni) > (BufSize-NRead)) exit
         NZInt = NZInt + ni
         IntRec = IntRec + 1
      enddo
      !if ( NZInt == 0 ) stop 'Error NZInt=0'
      if ( NZInt == 0 ) call quit('SetupBlkRead_VPHPP: Error NZInt=0')
      EnRec = IntRec - 1
      NBRec = (NRead+NZInt)/NZVINTBUFSIZE
      if (NBRec*NZVINTBUFSIZE/=(NRead+NZInt)) NBRec = NBRec + 1
      EnBlkRec_VPHPP(NBlk) = SRec + NBRec - 1
      BlkStInd_VPHPP(NBlk) = NRead + 1
      BlkEnInd_VPHPP(NBlk) = NRead + NZInt
      BlkSt_VPHPP(NBlk) = StRec
      BlkEn_VPHPP(NBlk) = EnRec
      SRec = SRec + NBRec
      NunRead = BufSize - NRead - NZInt
      NRead = NZVINTBUFSIZE - MOD(NunRead,NZVINTBUFSIZE)
      if ( NRead == NZVINTBUFSIZE ) NRead = 0
      if ( NRead /= 0 ) SRec = SRec - 1
      if ( IntRec > NPMO) exit
    enddo

end subroutine SetupBlkRead_VPHPP

subroutine EndIOMod

   implicit none

   ! ... Initilization Check ...

   if (.not.Init) then
      write(UOUT,*) 'EndIOMod: EndIOMod called Without Initilization'
      call flushout
      call quit('EndIOMod: EndIOMod called Without Initilization')
   endif
   Init = .FALSE.

   ! ... deallocate only if they are allocated ...
   ! ... when NOHSMR the following entities ...
   ! ... are not allocated                  ...

   if (allocated(NZIntRec_V2H2P)) deallocate(NZIntRec_V2H2P)
   if (allocated(NZIntRec_VPPPP)) deallocate(NZIntRec_VPPPP)
   if (allocated(NZIntRec_VHHHH)) deallocate(NZIntRec_VHHHH)
   if (allocated(NZIntRec_VHHPP)) deallocate(NZIntRec_VHHPP)
   if (allocated(NZIntRec_V2HPP)) deallocate(NZIntRec_V2HPP)
   if (allocated(NZIntRec_V2PHH)) deallocate(NZIntRec_V2PHH)
   if (allocated(NZIntRec_VPHHP)) deallocate(NZIntRec_VPHHP)
   if (allocated(NZIntRec_VPH2P)) deallocate(NZIntRec_VPH2P)
   if (allocated(NZIntRec_VPHPP)) deallocate(NZIntRec_VPHPP)
   if (allocated(NZIntRec_V2HHP)) deallocate(NZIntRec_V2HHP)
   if (allocated(NZIntRec_VHHHP)) deallocate(NZIntRec_VHHHP)

   if (allocated(StBlkRec_VPPPP)) deallocate(StBlkRec_VPPPP)
   if (allocated(EnBlkRec_VPPPP)) deallocate(EnBlkRec_VPPPP)
   if (allocated(BlkStInd_VPPPP)) deallocate(BlkStInd_VPPPP)
   if (allocated(BlkEnInd_VPPPP)) deallocate(BlkEnInd_VPPPP)
   if (allocated(BlkSt_VPPPP))    deallocate(BlkSt_VPPPP)
   if (allocated(BlkEn_VPPPP))    deallocate(BlkEn_VPPPP)

   if (allocated(StBlkRec_VPH2P)) deallocate(StBlkRec_VPH2P)
   if (allocated(EnBlkRec_VPH2P)) deallocate(EnBlkRec_VPH2P)
   if (allocated(BlkStInd_VPH2P)) deallocate(BlkStInd_VPH2P)
   if (allocated(BlkEnInd_VPH2P)) deallocate(BlkEnInd_VPH2P)
   if (allocated(BlkSt_VPH2P))    deallocate(BlkSt_VPH2P)
   if (allocated(BlkEn_VPH2P))    deallocate(BlkEn_VPH2P)

   if (allocated(StBlkRec_VPHPP)) deallocate(StBlkRec_VPHPP)
   if (allocated(EnBlkRec_VPHPP)) deallocate(EnBlkRec_VPHPP)
   if (allocated(BlkStInd_VPHPP)) deallocate(BlkStInd_VPHPP)
   if (allocated(BlkEnInd_VPHPP)) deallocate(BlkEnInd_VPHPP)
   if (allocated(BlkSt_VPHPP))    deallocate(BlkSt_VPHPP)
   if (allocated(BlkEn_VPHPP))    deallocate(BlkEn_VPHPP)

   return

end subroutine EndIOMod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine ReadFileFull(FileName)
!----------------------------------------------------------------------------------
!
! read big chunk of seleced integrals (accroding to FileName) together with indexes
!
!-------------------------------------------------------------------------------

   implicit none
   integer, intent(in) :: FileName
   integer :: StRec, EnRec
   integer :: StInd, EnInd
   integer :: RecInd

   select case (FileName)
      case (V2H2P)
         NZIntRec => NZIntRec_V2H2P
      case (VHHHH)
         NZIntRec => NZIntRec_VHHHH
      case (VHHPP)
         NZIntRec => NZIntRec_VHHPP
      case (V2HPP)
         NZIntRec => NZIntRec_V2HPP
      case (V2PHH)
         NZIntRec => NZIntRec_V2PHH
      case (VPHHP)
         NZIntRec => NZIntRec_VPHHP
      case (V2HHP)
         NZIntRec => NZIntRec_V2HHP
      case (VHHHP)
         NZIntRec => NZIntRec_VHHHP
      case default
         write(UOUT,*) 'ReadFileFull: Unable to Read Full File, wrong case FileName=',FileName
         call flushout
         call quit('ReadFileFull: Unable to Read Full File, wrong case !')
   end select

   StRec = StartRec(FileName)
   EnRec = EndngRec(FileName)
   do RecInd = StRec, EnRec
      StInd = (RecInd-StRec)*NZVINTBUFSIZE+1
      EnInd = StInd+NZVINTBUFSIZE-1
      read(UVINT,rec=RecInd) BufIntgrl(StInd:EnInd)
      read(UVIND,rec=RecInd) BufIndex(StInd:EnInd)
   enddo

   IntBuf => BufIntgrl(1:TotNZInt(FileName))
   IndBuf => BufIndex(1:TotNZInt(FileName))

   return

end subroutine ReadFileFull

subroutine ReadFileBlock(FileName,Block,St,En)
!---------------------------------------------------------------------------------
!
! read block of integrals & indexes from files UVINT, UVIND  
!
! reads VPPPP, VPH2P, VPHPP types of integrals (large - thus block reading necessary)
!
! returns St, En
!
!----------------------------------------------------------------------------------
   implicit none
   integer, intent(in)  :: FileName, Block
   integer, intent(out) :: St, En
   integer :: StRec, EnRec
   integer :: StInd, EnInd
   integer :: RecInd

   select case (FileName)
      case (VPPPP)
         NZIntRec => NZIntRec_VPPPP
         StRec = StBlkRec_VPPPP(Block)
         EnRec = EnBlkRec_VPPPP(Block)
         do RecInd = StRec, EnRec
            StInd = (RecInd-StRec)*NZVINTBUFSIZE+1
            EnInd = StInd+NZVINTBUFSIZE-1
            read(UVINT,rec=RecInd) BufIntgrl(StInd:EnInd)
            read(UVIND,rec=RecInd) BufIndex(StInd:EnInd)
         enddo
         StInd = BlkStInd_VPPPP(Block)
         EnInd = BlkEnInd_VPPPP(Block)
         IntBuf => BufIntgrl(StInd:EnInd)
         IndBuf => BufIndex(StInd:EnInd)
         St = BlkSt_VPPPP(Block)
         En = BlkEn_VPPPP(Block)
      case (VPH2P)
         NZIntRec => NZIntRec_VPH2P
         StRec = StBlkRec_VPH2P(Block)
         EnRec = EnBlkRec_VPH2P(Block)
         do RecInd = StRec, EnRec
            StInd = (RecInd-StRec)*NZVINTBUFSIZE+1
            EnInd = StInd+NZVINTBUFSIZE-1
            read(UVINT,rec=RecInd) BufIntgrl(StInd:EnInd)
            read(UVIND,rec=RecInd) BufIndex(StInd:EnInd)
         enddo
         StInd = BlkStInd_VPH2P(Block)
         EnInd = BlkEnInd_VPH2P(Block)
         IntBuf => BufIntgrl(StInd:EnInd)
         IndBuf => BufIndex(StInd:EnInd)
         St = BlkSt_VPH2P(Block)
         En = BlkEn_VPH2P(Block)
      case (VPHPP)
         NZIntRec => NZIntRec_VPHPP
         StRec = StBlkRec_VPHPP(Block)
         EnRec = EnBlkRec_VPHPP(Block)
         do RecInd = StRec, EnRec
            StInd = (RecInd-StRec)*NZVINTBUFSIZE+1
            EnInd = StInd+NZVINTBUFSIZE-1
            read(UVINT,rec=RecInd) BufIntgrl(StInd:EnInd)
            read(UVIND,rec=RecInd) BufIndex(StInd:EnInd)
         enddo
         StInd = BlkStInd_VPHPP(Block)
         EnInd = BlkEnInd_VPHPP(Block)
         IntBuf => BufIntgrl(StInd:EnInd)
         IndBuf => BufIndex(StInd:EnInd)
         St = BlkSt_VPHPP(Block)
         En = BlkEn_VPHPP(Block)
      case default
         write(UOUT,*) 'ReadFileBlock: Unable to Perform Block IO, wrong case=',FileName
         call flushout
         call quit('ReadFileBlock: Unable to Perform Block IO, wrong case !')
   end select

   return

end subroutine ReadFileBlock

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine OpenIntFiles
!----------------------------------------------------------------------------------------
!
! Open HSMRCC intgeral files:
!
!     - opening as new when running HSMRCC from Gamess (own HSMRCC integral files generation)
!        or when called from dirac_hsmrcc interface
!     - HSMRCC integral files MUST BE READY when running HSMRCC from DIRAC
!
! called from:  
!              InitHSMRBase (source/inihsbase.F90)
!                              dir_hscc_mods.F90
!
!-----------------------------------------------------------------------------------------
   implicit none
   character(len=10) :: stat
   integer:: recl_vind,recl_vint,recl_fhh,recl_fhp,recl_fpp,ios 
   logical:: exist_file

   if (get_iprint() >= 7) then
     call header('OpenIntFiles: going to open UVIND, UVINT, UFHH, UFHP, UFPP files',-1)
   endif

   call check_FileRec_param

   if (.not.Init) then
      write(UOUT,'(/,2x,a)') 'OpenIntFiles: IOMod not initialized !'
      call flushout
      call quit('OpenIntFiles: IOMod not initialized')
   endif

   if (OPTION(GENINT).or.get_dirac_hscc_interface_mode()) then
     stat = 'NEW'
     ! stat = 'UNKNOWN'
   else
     stat = 'OLD'
   endif

   if (get_iprint() >= 7) then
     write(UOUT,'(2X,A,A)') '...going to open files as :',stat
   endif

! ----------------------------------------------------
   inquire(file=FileVIND,exist=exist_file)
   if (.not.exist_file.and.OPTION(HSCC_DIRAC_MODE)) then
    call quit('OpenIntFiles: FileVIND does not exist!')
   endif
   
   if (exist_file.and.stat=='NEW') then
     call quit('OpenIntFiles: FileVIND exist, can not ne open as NEW !')
   endif

   recl_vind=(4*2*NZVINTBUFSIZE)/INTSFP
   open(UVIND,file=FileVIND,access='DIRECT',form='UNFORMATTED',status=stat,recl=recl_vind,iostat=ios)

   if (ios.ne.0) then
     write(UOUT,*) 'UVIND=',UVIND
     write(UOUT,*) 'FileVIND=',FileVIND
     write(UOUT,*) 'recl_vind=',recl_vind 
     call quit('OpenIntFiles: FileVIND opening failure !')
   endif

! ----------------------------------------------------
   inquire(file=FileVINT,exist=exist_file)
   if (.not.exist_file.and.OPTION(HSCC_DIRAC_MODE)) then
    call quit('OpenIntFiles: FileVINT does not exist!')
   endif

   if (exist_file.and.stat=='NEW') then
     call quit('OpenIntFiles: FileVINT exist, can not ne open as NEW !')
   endif

   recl_vint=4*2*NZVINTBUFSIZE
   open(UVINT,file=FileVINT,access='DIRECT',form='UNFORMATTED',status=stat,recl=recl_vint,iostat=ios)
! ----------------------------------------------------
   inquire(file=FileFHH,exist=exist_file)
   if (.not.exist_file.and.OPTION(HSCC_DIRAC_MODE)) then
    call quit('OpenIntFiles: FileFHH does not exist!')
   endif

   if (exist_file.and.stat=='NEW') then
     call quit('OpenIntFiles: FileFHH exist, can not ne open as NEW !')
   endif

   recl_fhh=4*2*NHMO*NHMO
   open(UFHH ,file=FileFHH ,access='DIRECT',form='UNFORMATTED',status=stat,recl=recl_fhh,iostat=ios)

! ----------------------------------------------------
   inquire(file=FileFHP,exist=exist_file)
   if (.not.exist_file.and.OPTION(HSCC_DIRAC_MODE)) then
    call quit('OpenIntFiles: FileFHP does not exist!')
   endif

   if (exist_file.and.stat=='NEW') then
     call quit('OpenIntFiles: FileFHP exist, can not ne open as NEW !')
   endif

   recl_fhp=4*2*NHMO*NPMO
   open(UFHP ,file=FileFHP ,access='DIRECT',form='UNFORMATTED',status=stat,recl=recl_fhp,iostat=ios)

! ----------------------------------------------------
   inquire(file=FileFPP,exist=exist_file)
   if (.not.exist_file.and.OPTION(HSCC_DIRAC_MODE)) then
    call quit('OpenIntFiles: FileFPP does not exist!')
   endif

   if (exist_file.and.stat=='NEW') then
     call quit('OpenIntFiles: FileFPP exist, can not ne open as NEW !')
   endif

   recl_fpp=4*2*NPMO*NPMO
   open(UFPP ,file=FileFPP ,access='DIRECT',form='UNFORMATTED',status=stat,recl=recl_fpp,iostat=ios)

!  ... print out info ...
   if (get_iprint() >= 7) then
     write(UOUT,'(/,2x,a,a)')   & 
     'UVIND, UVINT, UFHH, UFHP, UFPP files were open with status=',stat
     write(UOUT,'(a,i2,a,a18,a,i8)') 'UVIND=',UVIND,' 2-el.integrals, FileVIND=',FileVIND,' recl=',recl_vind
     write(UOUT,'(a,i2,a,a18,a,i8)') 'UVINT=',UVINT,' 2-el.integ. indexes, FileVINT=',FileVINT,' recl=',recl_vint
     write(UOUT,'(a,i2,a,a18,a,i8)') ' UFHH=',UFHH, ' 1-el. integ., FileFHH=',FileFHH,' recl=',recl_fhh
     write(UOUT,'(a,i2,a,a18,a,i8)') ' UFHP=',UFHP, ' 1-el. integ., FileFHP=',FileFHP,' recl=',recl_fhp
     write(UOUT,'(a,i2,a,a18,a,i8)') ' UFPP=',UFPP, ' 1-el. integ., FileFPP=',FileFPP,' recl=',recl_fpp
     write(UOUT,*) 
     call flushout
   endif

   return
end subroutine OpenIntFiles

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine OpenOtherFiles
!
!       called from InitHSMRIter (source/inihsiter.F90)
!
   implicit none
   integer :: recl_ntsize, ios

   if (.not.Init) then
      write(UOUT,*) 'OpenOtherFiles: IOMod not initialized'
      call flushout
      call quit('OpenOtherFiles: IOMod not initialized')
   endif

   if (get_iprint() >= 7) then
     write(UOUT,*) 'OpenOtherFiles: going to open FileHEFF,FileCOEF with lu: ',UHEFF,UCOEF
     write(UOUT,*) 'file names:',FileHEFF,FileCOEF
   endif

   open(UHEFF,file=FileHEFF,access='SEQUENTIAL',form='UNFORMATTED',status='UNKNOWN',iostat=ios)
   if (ios.ne.0) then
     call quit('error in opening FileHEFF !')
   endif

   open(UCOEF,file=FileCOEF,access='SEQUENTIAL',form='UNFORMATTED',status='UNKNOWN',iostat=ios)
   if (ios.ne.0) then
     call quit('error in opening FileCOEF !')
   endif

! .. direct access files
   recl_ntsize=(4*2*(1+NTSIZE))

   if (get_iprint() >= 7) then
!    write(UOUT,*) 'OpenOtherFiles: going to open FileTNEW,FileTOLD,FileTAMP files with recl_ntsize=',\\
!    recl_ntsize
   endif

   open(UTNEW,file=FileTNEW,access='DIRECT',form='UNFORMATTED',status='UNKNOWN',recl=recl_ntsize,iostat=ios)
   open(UTOLD,file=FileTOLD,access='DIRECT',form='UNFORMATTED',status='UNKNOWN',recl=recl_ntsize,iostat=ios)
   open(UTAMP,file=FileTAMP,access='DIRECT',form='UNFORMATTED',status='UNKNOWN',recl=recl_ntsize,iostat=ios)

   if (get_iprint() >= 7) then
     call header('OpenOtherFiles: UHEFF,UCOEF,UTNEW,UTOLD,UTAMP files open.',-1)
     write(UOUT,'(a,i2,a,a)') 'UHEFF=',UHEFF,' sequent. FileHEFF=',FileHEFF
     write(UOUT,'(a,i2,a,a)') 'UCOEF=',UCOEF,' sequent. FileCOEF=',FileCOEF
     write(UOUT,'(a,i2,a,a,a,i7)') 'UTNEW=',UTNEW,' direct FileTNEW=',FileTNEW,' recl=',recl_ntsize
     write(UOUT,'(a,i2,a,a,a,i7)') 'UTOLD=',UTOLD,' direct FileTOLD=',FileTOLD,' recl=',recl_ntsize
     write(UOUT,'(a,i2,a,a,a,i7)') 'UTAMP=',UTAMP,' direct FileTAMP=',FileTAMP,' recl=',recl_ntsize
     write(UOUT,*)  
     call flushout
   endif

   return
end subroutine OpenOtherFiles

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine CloseIntFiles
!-----------------------------------------------------------------------------
!                close VINT,VIND,FHH,FHP,FPP files
!-----------------------------------------------------------------------------

   implicit none

   if (.not.Init) then
     write(UOUT,*) 'CloseIntFiles: IOMod not initialized'
     call flushout
     call quit('CloseIntFiles: IOMod not initialized')
   endif

   close(UVINT)
   close(UVIND)
   close(UFHH)
   close(UFHP)
   close(UFHP)
   close(UFPP)

   if (get_iprint() >= 7) then
     write(UOUT,'(/,2x,a,/)') 'CloseIntFiles: UVINT, UVIND, UFHH, UFHP, UFHP, UFPP files closed.'
   endif

   return

end subroutine CloseIntFiles

subroutine CloseOtherFiles

   implicit none

   if (.not.Init) then
      write(UOUT,*) 'CloseOtherFiles: IOMod not initialized'
      call flushout
      call quit('CloseOtherFiles: IOMod not initialized')
   endif

   close(UTAMP)
   close(UHEFF)
   close(UCOEF)
   close(UCNVR)
   close(UTNEW,status='DELETE')
   close(UTOLD,status='DELETE')

   if (get_iprint() >= 7) then
     write(UOUT,'(/,2x,a,/)') 'CloseOtherFiles: UTAMP,UHEFF,UCOEF,UCNVR files closed, UTNEW,UTOLD deleted.'
   endif

   return

end subroutine CloseOtherFiles

end module IOmodule

