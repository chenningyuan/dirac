!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module GMSMod
!------------------------------------------------------------------------
!
!     Necessary stuff read from Gamess (integral) files
!
! also for DIRAC mode
!
!------------------------------------------------------------------------
   use Parameters
   use Options
   use IOModule
   use SymInfHSCC
   use OrbMod
   implicit none

  ! integer :: NAO, NSO, NFRZV
  ! integer :: NAO=0, NSO=0, NFRZV=0
  ! integer :: NAO=-1, NSO=-1, NFRZV=-1
   integer :: NAO=-1, NSO=-1, NFRZV=0

  ! integer :: NVMO !  it's in Parameters
   double precision :: TotEn, NucEn, TrfTol

   integer :: NTMO_GMS, NCMO_GMS, NVMO_GMS, NEMO_GMS, NHOLE_GMS, NPTCL_GMS

   integer :: IUNIQ, LABELSIZE, NINTMX

   ! ... arrays for GAMESS integral files ...
   double precision, allocatable :: XXPREAD(:)
   integer, allocatable :: IXPREAD(:)

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine InitGMSMod
!---------------------------------------------------------------------------
!
! Initialize Gamess stuff 
!
! Called from InitHSMRBase (source/inihsbase.F90) - if not in the DIRAC mode !
!
!---------------------------------------------------------------------------

   implicit none

   ! ... Obtain information for Variables of GMSMod       ...
   ! ... There are two ways to obtain this info           ...
   ! ... if the run involves GENINT, this info has to     ...
   ! ... read from FileMR2EINT, otherwise it can be read  ...
   ! ... from FileVINF which has a copy of this info      ...

   if (get_iprint() >= 6) then
    if (.not.OPTION(HSCC_DIRAC_MODE)) then
     Call Header('Output from InitGMSMod routine - GAMESS mode',-1)
    else
     Call Header('Output from InitGMSMod routine - DIRAC mode',-1)
    endif
   endif

   call InitReadGMS

!  ... final print-out of what has been read inside - simple variables & arrays
   if (get_iprint() >=7) then
    call print_gms_data
   endif

   return
end subroutine InitGMSMod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine EndGMSMod
!--------------------------------------------------------------------
!
! Deallocate some arrays.
!
! Called from : EndHSMR - source/endhsmr.F90
!               read_moltra_files (dir_hscc_mods.F90)
!
! TODO: include status checking of deallocations
!
!--------------------------------------------------------------------
   implicit none

   if (get_iprint() >=8) then
      write(UOUT,'(/,2x,a)') 'EndGMSMod: going to deallocate HSCC specific data.'
   endif

   if (allocated(SpinOfOrb))     deallocate(SpinOfOrb)
   if (allocated(SpinOfActvOrb)) deallocate(SpinOfActvOrb)
   if (allocated(OccGMSRef))     deallocate(OccGMSRef)
   if (allocated(OccRef))        deallocate(OccRef)
   if (allocated(IrrepOrb ))     deallocate(IrrepOrb )
   if (allocated(IrpLabOrb ))    deallocate(IrpLabOrb)
   if (allocated(OrbType))       deallocate(OrbType)
   if (allocated(OrbSpace))      deallocate(OrbSpace)
   if (allocated(OrbEnergy))     deallocate(OrbEnergy)
   if (allocated(OrbCoeff))      deallocate(OrbCoeff)
   if (allocated(HCoreMatrix))   deallocate(HCoreMatrix)

 !  if (allocated(MOLTRA_indx))   deallocate(MOLTRA_indx)

   if (get_iprint() >=8) then
     write(UOUT,'(/,2x,a,/)') 'EndGMSMod: HSCC specific arrays (SpinOfOrb,SpinOfActvOrb,..) deallocated'
   endif

   if (allocated(XXPREAD)) then
      deallocate(XXPREAD)
      if (get_iprint() >=8) then
         write(UOUT,'(/,2x,a)') 'EndGMSMod: GMS specific array XXPREAD deallocated.'
      endif
   endif

   if (allocated(IXPREAD)) then
     deallocate(IXPREAD)
     if (get_iprint() >=8) then
        write(UOUT,'(/,2x,a)') 'EndGMSMod: GMS specific array IXPREAD deallocated.'
     endif
   endif

   return
end subroutine EndGMSMod


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine alloc_VINF_arrays
!----------------------------------------------------------------------
!
! allocate 1-electron data
! called from InitReadVINF, InitReadMREINT ....
!
!  ... needs: NSO, NVMO, NAO
!
!----------------------------------------------------------------------
  implicit none

  if (get_iprint() >=7) then
    write(UOUT,'(/,2x,a,/)') '=== Output from alloc_VINF_arrays - before allocations ==='
    write(UOUT,'(2x,a,3i4)') 'used variables - NSO,NVMO,NAO:', NSO,NVMO,NAO
    call flushout
  endif

  if (NSO <=0 .or. NVMO <=0 .or. NAO <=0) then
    write(UOUT,'(2x,a,3i4)')  & 
    'alloc_bunch: used bad variables ! - NSO,NVMO,NAO:',NSO,NVMO,NAO
    call flushout
    call quit('alloc_bunch: wrong sizef for allocation !')
  endif

  if (.not.allocated(SpinOfOrb))     allocate(SpinOfOrb(NSO))
  if (.not.allocated(OccGMSRef))     allocate(OccGMSRef(NSO))
  if (.not.allocated(IrrepOrb))      allocate(IrrepOrb(NSO))
  if (.not.allocated(IrpLabOrb))     allocate(IrpLabOrb(NSO))
  if (.not.allocated(OrbType))       allocate(OrbType(NSO))
  if (.not.allocated(OrbSpace))      allocate(OrbSpace(NSO))
  if (.not.allocated(OrbEnergy))     allocate(OrbEnergy(NSO))
  if (.not.allocated(OrbCoeff))      allocate(OrbCoeff(NAO,NSO))
  if (.not.allocated(HCoreMatrix))   allocate(HCoreMatrix(NSO,NSO))

  if (.not.allocated(SpinOfActvOrb)) allocate(SpinOfActvOrb(NVMO))
  if (.not.allocated(OccRef))        allocate(OccRef(NSO))

  if ( OPTION(DIRAC_MODE) ) then
  !   ... allocate only when in DIRAC preparing mode for HSMRCC ...
     if (.not.allocated(MOLTRA_indx)) allocate(MOLTRA_indx(NSO))
  endif

  if (get_iprint() >=7) then
    write(UOUT,'(2x,a,/)') 'alloc_bunch: SpinOfOrb, SpinOfActvOrb and other data allocated...'
    call flushout
  endif

end subroutine alloc_VINF_arrays

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine InitReadGMS
!--------------------------------------------------------------------------
!
!       read either Gamess integral files; or VINF file when present
!
!--------------------------------------------------------------------------
   implicit none
   integer NXX, IXX

   if (get_iprint() >=7) then
     call header('Output from InitReadGMS',-1)
   endif

   ! ...
  ! if (OPTION(GENINT)) then
   if (OPTION(GENINT).and..not.OPTION(HSCC_DIRAC_MODE)) then
     if (get_iprint() >= 7) then
       write(UOUT,'(2x,a,/)')   & 
       '...going to call InitReadMREINT - read core info from GMS integral files; not DIRAC mode'
     endif
     call InitReadMREINT
   else
     if (get_iprint() >= 7) then
       write(UOUT,'(2x,a,/)') '...going to call InitReadVINF - read core info from VINF file'
     endif
     call InitReadVINF
   endif

   ! ... allocate gamess integral buffers 

! Miro:  we do need it when OPTION(GENINT)=.true. !
   ! ... set up size for gamess integral labels
   NXX = NINTMX
   IF (LABELSIZE.EQ.1) THEN
      IF(INTSFP.EQ.2) IXX = NXX
      IF(INTSFP.EQ.1) IXX = (NXX+1)/2
   ELSE IF (LABELSIZE.EQ.2) THEN
      IF(INTSFP.EQ.2) IXX = NXX*2
      IF(INTSFP.EQ.1) IXX = NXX
   ENDIF

   ! ... double precision array for Gamess integrals ...
   if (.not.allocated(XXPREAD)) allocate(XXPREAD(NXX))

   ! ... integer array for Gamess integral labels ...
   if (.not.allocated(IXPREAD)) allocate(IXPREAD(IXX))

   !             ... control print-out ... 
   write(UOUT,'(/,a,/,a,i2,a,i2,/,2(a,i7))')    &
   '----------------------------------------------------------------------------',    &
   'InitReadGMS: read LABELSIZE=',LABELSIZE,'; set INTSFP (parameters.F90) ',INTSFP,   & 
   '...read NXX=',NXX,' and set IXX=',IXX

#if defined (ARCH32BIT) 
   write(UOUT,'(a)') '...we are on the 32-bit architecture '
#else
   write(UOUT,'(a)') '...we are on the 64-bit architecture '
#endif
   write(UOUT,'(a)')  &
   '-----------------------------------------------------------------------------'

end subroutine  InitReadGMS

subroutine InitReadDIRAC  ! ... NOT ACTIVE !!!
!> read DIRAC file with data

character*12 :: dirac_file='Dirac_HSMRCC'
real*8    :: ENSCF_DIR
integer   :: ludir = 18, iprint_dirac, NELEC_DIR, NSP, NORBT, NCORE_DIR, NACT_DIR,  & 
             numirrep_dirac,irrep_product_dirac(1:8,1:8)
logical   :: is_file, dir_spfr

inquire(file=dirac_file, exist=is_file)
if (.not.is_file) then
  write(UOUT,*) 'InitReadDIRAC: problem with file "dirac_file"=',dirac_file
  call quit('InitReadDIRAC: no dirac input file !...')
endif

! ... trunk/hsmrcc/dir_hscc_mods.F90
open(ludir,file=dirac_file,form='unformatted',status='old', access='sequential')

read(ludir) iprint_dirac

read(ludir) NELEC_DIR

read(ludir) dir_spfr

read(ludir) ENSCF_DIR

read(ludir) NSP
read(ludir) NORBT
read(ludir) NCORE_DIR
read(ludir) NACT_DIR
read(ludir) numirrep_dirac
read(ludir) irrep_product_dirac

write(UOUT,*) 'read iprint_dirac=',iprint_dirac
write(UOUT,*) 'NELEC_DIR=',NELEC_DIR
write(UOUT,*) 'dir_spfrdir_spfr=',dir_spfr
write(UOUT,*) 'ENSCF_DIR=',ENSCF_DIR
write(UOUT,*) 'NSP=',NSP
write(UOUT,*) 'NORBT=',NORBT
write(UOUT,*) 'NCORE_DIR=',NCORE_DIR
write(UOUT,*) 'NACT_DIR=',NACT_DIR
write(UOUT,*) 'numirrep_dirac=',numirrep_dirac
write(UOUT,*) 'irrep_product_dirac',irrep_product_dirac

! ... assign DIRAC HSCC variables ...
NAO = NORBT
NEMO = NAO

! ... active space - double excitations
NPTCL = 2
NHOLE = 2

!NACT = NACT_DIR
NCMO = NCORE_DIR
write(UOUT,*) 'assigned NCMO=',NCMO

NSO = NSP
TotEn = ENSCF_DIR

! ... irrep_product - spin-free !
 irrep_product = irrep_product_dirac


! ...  prepare variable  ...
  NVMO = NHOLE + NPTCL

 call alloc_VINF_arrays

! ... read arrays ...

 read(ludir) SpinOfOrb

 read(ludir) OccGMSRef

 read(ludir) OrbType

 read(ludir) OrbSpace

 read(ludir) IrrepOrb

 read(ludir) IrpLaborb

 read(ludir) OrbEnergy

close(ludir) 

!       ... update ...
      SpinOfActvOrb = SpinOfOrb(NCMO+1:NCMO+NVMO)
      OccRef = OccGMSRef

! ... other  symmetry stuff ... will it be sufficient ?
      totally_symmetric_irrep = 1

      numirrep = numirrep_dirac

end subroutine  InitReadDIRAC


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine InitReadVINF
!---------------------------------------------------------------------
!
! read existing VINF file (when restarting)
!
! called from InitReadGMS
!
!---------------------------------------------------------------------

implicit none

      if (get_iprint() >=7) then
        write(UOUT,'(/,6x,a,/)') 'InitReadVINF: going to read FileVINF file...'
      endif

      !******************************************************
      ! ... Read the information from FileVINF ...
      ! ... This has been already opened       ...
      !******************************************************

      rewind(UVINF)

      write(UOUT,*) 'rewind of of FileVINF done'

      read(UVINF)
      read(UVINF)
      read(UVINF)
      read(UVINF)

      read(UVINF) point_group, naxis, norder, numirrep, irrep_product, &
                  totally_symmetric_irrep, irrep_label
      read(UVINF)

      read(UVINF) NAO, NSO, NFRZV, TotEn, NucEn, TrfTol, IUNIQ, LABELSIZE, NINTMX
 
 !      write(UOUT,*) 'till here FileVINF read'

      ! ... Check for Validity of Data ...

      if ((NAO <= 0).or.(NSO <= 0).or.(NFRZV < 0).or.(NSO > 2*(NAO - NFRZV)))then
         write(UOUT,*) 'Invalid Data in File UVINF'
         call flushout
         call quit(': Invalid Data in File UVINF !')
      endif

      if (((IUNIQ.NE.0).and.(IUNIQ.NE.1)).or.              &
          ((LABELSIZE.NE.1).AND.(LABELSIZE.NE.2)).or.      &
           (NINTMX.LE.0)) then
         write(UOUT,*) ' IUNIQ=',IUNIQ,' LABELSIZE=',LABELSIZE,' NINTMX=',NINTMX
         write(UOUT,*) 'Invalid Data in File UVINF, 2! '
         call flushout
         call quit(': Invalid Data in File UVINF 2 !')
      endif

      ! mi: NTMO is zero 
      if (NTMO /= NSO.and..not.OPTION(HSCC_DIRAC_MODE)) then
         write(UOUT,*) ' NTMO=',NTMO,' NSO=',NSO
         write(UOUT,*) 'Current NTMO is not valid for integral file'
         call flushout
         call quit('InitReadVINF: Current NTMO is not valid for integral file')
      endif

      ! ... Allocate for next read ...
      call alloc_VINF_arrays

      ! ... Read Orbital Energies and MO Coeffeints ... variables must be allocated !

      read(UVINF) SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb, &
                  OrbType, OrbSpace, OrbEnergy, OrbCoeff, HCoreMatrix

      SpinOfActvOrb = SpinOfOrb(NCMO+1:NCMO+NVMO)
      OccRef = OccGMSRef

      if (get_iprint() >=7) then
          call print_gms_data
      endif

end subroutine InitReadVINF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_gms_data
!-------------------------------------------------------------------------------------------------
!
!> print allocated data read from GAMESS/DIRAC files
!
!  Prints  SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb, OrbType, OrbSpace, OrbEnergy, OrbCoeff, 
!
!  and HCoreMatrix if wished
!
!  called from both DIRAC/HSCC
!
!-------------------------------------------------------------------------------------------------
integer :: i
logical :: print_hcore_mtx = .false.

call header('Output from "print_gms_data"',-1)

    if (OPTION(HSCC_DIRAC_MODE)) then
      write(UOUT,'(/,2x,a,/)') '===== InitGMSMod: obtained variables - from DIRAC ====='
    else if (.not.OPTION(HSCC_DIRAC_MODE).and..not.OPTION(DIRAC_MODE)) then
      write(UOUT,'(/,2x,a,/)') '===== InitGMSMod: obtained variables - from GAMESS ====='
    else if (OPTION(DIRAC_MODE)) then
      write(UOUT,'(/,2x,a,/)') '===== InitGMSMod: in DIRAC_mode ====='
    else
     call quit('dead brach ! fix it !') 
    endif

!   -----
    write(UOUT,'(8x,a,i4)')     '   NAO=',NAO
    write(UOUT,'(8x,a,i4)')     '   NSO=',NSO
    write(UOUT,'(8x,a,i4)')     ' NFRZV=',NFRZV
    write(UOUT,'(8x,a,i4)')     '  NCMO=',NCMO
    write(UOUT,'(8x,a,i4)')     '  NEMO=',NEMO
    write(UOUT,'(8x,a,i4)')     ' NHOLE=',NHOLE
    write(UOUT,'(8x,a,i4)')     ' NPTCL=',NPTCL
    write(UOUT,'(8x,a,f18.13)') ' TotEn=',TotEn
    write(UOUT,'(8x,a,f18.13)') ' NucEn=',NucEn
    write(UOUT,'(8x,a,f18.13)') 'TrfTol=',TrfTol

write(UOUT,'(/,2x,a,/,2x,a)') 'Data arrays(i, SpinOfOrb, IrrepOrb, IrpLabOrb, OrbType, OrbSpace, OrbEnergy)', &
'--------------------------------------------------------------------------------------------------------------'

write(UOUT,'(4x,a,i5)') 'range of arrays, NSO=',NSO 

if (NSO==0) then 
 write(UOUT,*) 'Error ! NSO=',NSO
 call quit('InitReadVINF: NSO==0 !')
endif

do i=1, NSO
  call print_gms_data_line(i)
enddo


if (print_hcore_mtx) then
! ... print HCoreMatrix
  write(UOUT,'(/,2x,a)') '*** HCoreMatrix ***'
  call prqmat(HCoreMatrix,NSO,NSO,NSO,NSO,1,1,UOUT)
endif

end subroutine print_gms_data

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine  print_gms_data_line(i)
!
! ... prints one line of indexed data
! 
integer, intent(in) :: i

  if (.not.allocated(MOLTRA_indx)) then
  ! ... print without MOLTRA  spinor index

   write(UOUT,1298)   i,  & 
           SpinOfOrb(i),  & 
           OccGMSRef(i),  & 
            IrrepOrb(i),  & 
           IrpLabOrb(i),  &
             OrbType(i),  &
            OrbSpace(i),  &
           OrbEnergy(i)
  else
  ! ... print with MOLTRA spinor index

   write(UOUT,1299)   i,  & 
           SpinOfOrb(i),  & 
           OccGMSRef(i),  & 
            IrrepOrb(i),  & 
           IrpLabOrb(i),  &
             OrbType(i),  &
            OrbSpace(i),  &
           OrbEnergy(i),  & 
         MOLTRA_indx(i)

  endif

  1298 format (i3,1x,'Spin=',i2,1x,'Occ=',i1,1x,'Irp=',i1,1x,      & 
     'IrpLab=',a3,'Type=',i2,1x,'Space=',i2,1x,'OrbEn=',f12.6)

  1299 format (i3,1x,'Spin=',i2,1x,'Occ=',i1,1x,'Irp=',i1,1x,      & 
     'IrpLab=',a3,'Type=',i2,1x,'Space=',i2,1x,'OrbEn=',f12.6,1x,'MOLTRA indx=',i4)

end subroutine print_gms_data_line


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine InitReadMREINT
!----------------------------------------------------------------------------
!
! read data from GAMESS integral files, FileMR1EINT and FileMR2EINT
!
! called from InitReadGMS when OPTION(GENINT)=T
!
!----------------------------------------------------------------------------
   implicit none

   integer j1, j2, j3, j4, j5, j6, j7, j8, j9
   double precision :: d1, d2, d3
   logical :: errflag, existMR1,existMR2
  ! integer NXX, IXX

      ! ...****************************************************************...
      ! ... Open & Read Information Present in FileMR1EINT and FileMR2EINT ...
      ! ...****************************************************************...
      ! ...                                                                ...
      ! ...                                                                ...
      ! ... FileMR1EINT contains following information                     ...
      ! ... ---------------------------------------------------------------...
      ! ... Record #1 : NAO, NSO, NFRZV, TotEn, NucEn, TrfTol,             ...
      ! ...             NTMO, NCMO, NVMO, NEMO, NHOLE, NPTCL               ...
      ! ... Record #2 : SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb          ...
      ! ...             OrbType, OrbSpace, OrbEnergy, OrbCoeff             ...
      ! ... Record #3 : HCoreMatrix (1-e integrals)                        ...
      ! ...                                                                ...
      ! ...                                                                ...
      ! ... FileMR2EINT contains following information                     ...
      ! ... ---------------------------------------------------------------...
      ! ... Record #1 : NAO, NSO, NFRZV, TotEn, NucEn, TrfTol,             ...
      ! ...             NTMO, NCMO, NVMO, NEMO, NHOLE, NPTCL               ...
      ! ... Record #2 : SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb          ...
      ! ...             OrbType, OrbSpace, OrbEnergy, OrbCoeff             ...
      ! ... Record #3 : IUNIQ, LABELSIZE, NINTMX                           ...
      ! ... Record #4 onwards anti-symmetric 2-e integrals                 ...
      ! ...                                                                ...

      if (get_iprint() >=7) then
       write(UOUT,*) 'InitReadMRINT: going to read GAMESS integral files...'
        write(UOUT,'(2x,a,a)') ' Going to open FileMR1EINT=',FileMR1EINT
        write(UOUT,'(2x,a,a)') ' Going to open FileMR2EINT=',FileMR2EINT
      endif

      ! ... check the presence ...
      inquire(file=FileMR1EINT, exist=existMR1)
      inquire(file=FileMR2EINT, exist=existMR2)
      if (.not.(existMR1.or.existMR2)) then
        call quit('Gamess integral file(s) not exist...')
      endif

      open(UMR1EINT,file=FileMR1EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')
      open(UMR2EINT,file=FileMR2EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')

      if (get_iprint() >=8) then
        write(UOUT,*) 'InitGMSMod: files FileMR1EINT, FileMR2EINT are open'
      endif

      rewind(UMR1EINT)
      rewind(UMR2EINT)

      ! ... Read First Record ...

      read(UMR1EINT) j1, j2, j3, d1, d2, d3, j4, j5, j6, j7, j8, j9

      ! ... Check for Validity of Data ...

      errflag = (j1 <= 0).or.(j2 <= 0).or.(j3 < 0)
      errflag = errflag.or.(j2 > 2*(j1 - j3))
      errflag = errflag.or.(j2.ne.j4)
      errflag = errflag.or.(j4.ne.(j5+j6+j7))
      errflag = errflag.or.(j6.ne.(j8+j9))

      if (errflag) then
         write(UOUT,*) 'Invalid Data in File UMR1EINT'
         write(UOUT,*) "file = ",FileMR1EINT
         write(UOUT,*) j1, j2, j3, j4, j5, j6, j7, j8, j9
         call flushout
         call quit('Invalid Data in File UMR1EINT')
      endif

      read(UMR2EINT) j1, j2, j3, d1, d2, d3, j4, j5, j6, j7, j8, j9

      ! ... Check for Validity of Data ...

      errflag = (j1 <= 0).or.(j2 <= 0).or.(j3 < 0)
      errflag = errflag.or.(j2 > 2*(j1 - j3))
      errflag = errflag.or.(j2.ne.j4)
      errflag = errflag.or.(j4.ne.(j5+j6+j7))
      errflag = errflag.or.(j6.ne.(j8+j9))

      if (errflag) then
         write(UOUT,*) 'Invalid Data in File UMR1EINT'
         write(UOUT,*) "file = ",FileMR1EINT
         write(UOUT,*) j1, j2, j3, j4, j5, j6, j7, j8, j9
         call flushout
         call quit('errflag....Invalid Data in File UMR1EINT')
      endif

      ! ... set the variables in module ...

      NAO = j1
      NSO = j2
      NFRZV = j3
      TotEn = d1
      NucEn = d2
      TrfTol = d3

      NTMO_GMS = j4
      NCMO_GMS = j5
      NVMO_GMS = j6
      NEMO_GMS = j7
      NHOLE_GMS = j8
      NPTCL_GMS = j9

!
!      ...  Set up from gamess
!
      if (NCMO .eq.-1) NCMO = NCMO_GMS
      if (NEMO .eq.-1) NEMO = NEMO_GMS
      if (NHOLE.eq.-1) NHOLE = NHOLE_GMS
      if (NPTCL.eq.-1) NPTCL = NPTCL_GMS

      NVMO = NHOLE + NPTCL

      ! ... Allocate for next read ...

      call alloc_VINF_arrays

      ! ... Read next Record ...

      read(UMR1EINT)
      read(UMR2EINT) SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb,   &
                     OrbType, OrbSpace, OrbEnergy, OrbCoeff

      SpinOfActvOrb = SpinOfOrb(NCMO+1:NCMO+NVMO)
      OccRef = OccGMSRef

      ! ... Read Third Record ...

      read(UMR1EINT) HCoreMatrix

      read(UMR2EINT) IUNIQ, LABELSIZE, NINTMX

      if ( get_iprint() >=8 ) then
        write(UOUT,'(2x,a,i2,a,i2,a,i8)')  &
        ' InitGMSMod: from UMR2EINT read IUNIQ=',IUNIQ,  &
        ' LABELSIZE=',LABELSIZE,' NINTMX=',NINTMX
      endif

      if (((IUNIQ.NE.0).and.(IUNIQ.NE.1)).or.              &
          ((LABELSIZE.NE.1).AND.(LABELSIZE.NE.2)).or.      &
           (NINTMX.LE.0)) then
         write(UOUT,*) 'Invalid Data in File UMR2EINT'
         call flushout
         call quit('Invalid Data in File UMR2EINT')
      endif

      ! ... close files ...
      close(UMR1EINT) 
      close(UMR2EINT)

      if (get_iprint() >=7) then
        call print_gms_data
      endif

      ! read gamess symmetry info ...
      call read_gamess_syminf

end subroutine InitReadMREINT

end module GMSMod
