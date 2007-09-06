!
!
!
!
!
module GMSMod

   use Parameters
   use Options
   use IOModule
   use SymInf
   use OrbMod
   implicit none
   integer :: NAO, NSO, NFRZV
   double precision :: TotEn, NucEn, TrfTol
   integer :: NTMO_GMS, NCMO_GMS, NVMO_GMS, NEMO_GMS, NHOLE_GMS, NPTCL_GMS
   integer :: IUNIQ, LABELSIZE, NINTMX
   double precision, allocatable :: XXPREAD(:)
   integer, allocatable :: IXPREAD(:)


contains

subroutine InitGMSMod

   implicit none
   integer j1, j2, j3, j4, j5, j6, j7, j8, j9
   double precision :: d1, d2, d3
   logical :: errflag
   integer NXX, IXX

   ! ... Obtain information for Variables of GMSMod       ...
   ! ... There are two ways to obtain this info           ...
   ! ... if the run involves GENINT, this info has to     ...
   ! ... read from FileMR2EINT, otherwise it can be read  ...
   ! ... from FileVINF which has a copy of this info      ...


  print *,'in InitGMSMod routine !'

   if (OPTION(GENINT)) then

      ! ... Open & Read Information Present in FileMR1EINT and FileMR2EINT ...
      ! ...                                                                ...
      ! ... FileMR1EINT contains following information                     ...
      ! ... Record #1 : NAO, NSO, NFRZV, TotEn, NucEn, TrfTol,             ...
      ! ...             NTMO, NCMO, NVMO, NEMO, NHOLE, NPTCL               ...
      ! ... Record #2 : SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb          ...
      ! ...             OrbType, OrbSpace, OrbEnergy, OrbCoeff             ...
      ! ... Record #3 : HCoreMatrix (1-e integrals)                        ...
      ! ...                                                                ...
      ! ... FileMR2EINT contains following information                     ...
      ! ... Record #1 : NAO, NSO, NFRZV, TotEn, NucEn, TrfTol,             ...
      ! ...             NTMO, NCMO, NVMO, NEMO, NHOLE, NPTCL               ...
      ! ... Record #2 : SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb          ...
      ! ...             OrbType, OrbSpace, OrbEnergy, OrbCoeff             ...
      ! ... Record #3 : IUNIQ, LABELSIZE, NINTMX                           ...
      ! ... Record #4 onwards anti-symmetric 2-e integrals                 ...

      open(UMR1EINT,file=FileMR1EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')
      open(UMR2EINT,file=FileMR2EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')

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
         call flush(UOUT)
         call abort
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
         call flush(UOUT)
         call abort
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
!           Set up from gamess
!
      if (NCMO .eq.-1) NCMO = NCMO_GMS
      if (NEMO .eq.-1) NEMO = NEMO_GMS
      if (NHOLE.eq.-1) NHOLE = NHOLE_GMS
      if (NPTCL.eq.-1) NPTCL = NPTCL_GMS

      NVMO = NHOLE + NPTCL

      ! ... Allocate for next read ...

      if (.not.allocated(SpinOfOrb))     allocate(SpinOfOrb(NSO))
      if (.not.allocated(SpinOfActvOrb)) allocate(SpinOfActvOrb(NVMO))
      if (.not.allocated(OccGMSRef))     allocate(OccGMSRef(NSO))
      if (.not.allocated(OccRef))        allocate(OccRef(NSO))
      if (.not.allocated(IrrepOrb))      allocate(IrrepOrb(NSO))
      if (.not.allocated(IrpLabOrb))     allocate(IrpLabOrb(NSO))
      if (.not.allocated(OrbType))       allocate(OrbType(NSO))
      if (.not.allocated(OrbSpace))      allocate(OrbSpace(NSO))
      if (.not.allocated(OrbEnergy))     allocate(OrbEnergy(NSO))
      if (.not.allocated(OrbCoeff))      allocate(OrbCoeff(NAO,NSO))
      if (.not.allocated(HCoreMatrix))   allocate(HCoreMatrix(NSO,NSO))

      ! ... Read Second Record ...

      read(UMR1EINT)
      read(UMR2EINT) SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb,   &
                     OrbType, OrbSpace, OrbEnergy, OrbCoeff
      SpinOfActvOrb = SpinOfOrb(NCMO+1:NCMO+NVMO)
      OccRef = OccGMSRef

      ! ... Read Third Record ...

      read(UMR1EINT) HCoreMatrix
      read(UMR2EINT) IUNIQ, LABELSIZE, NINTMX

      if (((IUNIQ.NE.0).and.(IUNIQ.NE.1)).or.              &
          ((LABELSIZE.NE.1).AND.(LABELSIZE.NE.2)).or.      &
           (NINTMX.LE.0)) then
         write(UOUT,*) 'Invalid Data in File UMR2EINT'
         call flush(UOUT)
         call abort
      endif

      close(UMR1EINT)
      close(UMR2EINT)

      ! read gamess symmetry info ...
      call read_gamess_syminf

   else

      ! ... Read the information from FileVINF ...
      ! ... This has been already opened       ...

      rewind(UVINF)

      read(UVINF)
      read(UVINF)
      read(UVINF)
      read(UVINF)

      read(UVINF) point_group, naxis, norder, numirrep, irrep_product, &
                  totally_symmetric_irrep, irrep_label
      read(UVINF)

      read(UVINF) NAO, NSO, NFRZV, TotEn, NucEn, TrfTol, IUNIQ, LABELSIZE, NINTMX

      ! ... Check for Validity of Data ...

      if ((NAO <= 0).or.(NSO <= 0).or.(NFRZV < 0).or.(NSO > 2*(NAO - NFRZV)))then
         write(UOUT,*) 'Invalid Data in File UVINF'
         call flush(UOUT)
         call abort
      endif

      if (((IUNIQ.NE.0).and.(IUNIQ.NE.1)).or.              &
          ((LABELSIZE.NE.1).AND.(LABELSIZE.NE.2)).or.      &
           (NINTMX.LE.0)) then
         write(UOUT,*) 'Invalid Data in File UVINF'
         call flush(UOUT)
         call abort
      endif

      if (NTMO /= NSO) then
         write(UOUT,*) 'Current NTMO is not valid for integral file'
         call flush(UOUT)
         call abort
      endif

      ! ... Allocate for next read ...

      if (.not.allocated(SpinOfOrb))     allocate(SpinOfOrb(NSO))
      if (.not.allocated(SpinOfActvOrb)) allocate(SpinOfActvOrb(NVMO))
      if (.not.allocated(OccGMSRef))     allocate(OccGMSRef(NSO))
      if (.not.allocated(OccRef))        allocate(OccRef(NSO))
      if (.not.allocated(IrrepOrb))      allocate(IrrepOrb(NSO))
      if (.not.allocated(IrpLabOrb))     allocate(IrpLabOrb(NSO))
      if (.not.allocated(OrbType))       allocate(OrbType(NSO))
      if (.not.allocated(OrbSpace))      allocate(OrbSpace(NSO))
      if (.not.allocated(OrbEnergy))     allocate(OrbEnergy(NSO))
      if (.not.allocated(OrbCoeff))      allocate(OrbCoeff(NAO,NSO))
      if (.not.allocated(HCoreMatrix))   allocate(HCoreMatrix(NSO,NSO))

      ! ... Read Orbital Energies and MO Coeffeints ...

      read(UVINF) SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb, &
                  OrbType, OrbSpace, OrbEnergy, OrbCoeff, HCoreMatrix
      SpinOfActvOrb = SpinOfOrb(NCMO+1:NCMO+NVMO)
      OccRef = OccGMSRef

      ! ... Set other variables in GMSMod

   endif

   ! ... allocate gamess integral buffers

   NXX = NINTMX
   IF (LABELSIZE.EQ.1) THEN
      IF(INTSFP.EQ.2) IXX = NXX
      IF(INTSFP.EQ.1) IXX = (NXX+1)/2
   ELSE IF (LABELSIZE.EQ.2) THEN
      IF(INTSFP.EQ.2) IXX = NXX*2
      IF(INTSFP.EQ.1) IXX = NXX
   ENDIF

   if (.not.allocated(XXPREAD)) allocate(XXPREAD(NXX))
   if (.not.allocated(IXPREAD)) allocate(IXPREAD(IXX))

   return

end subroutine InitGMSMod

subroutine EndGMSMod

   implicit none

   if (allocated(SpinOfOrb))     deallocate(SpinOfOrb)
   if (allocated(SpinOfActvOrb)) deallocate(SpinOfOrb)
   if (allocated(OccGMSRef))     deallocate(OccGMSRef)
   if (allocated(OccRef))        deallocate(OccRef)
   if (allocated(IrrepOrb ))     deallocate(IrrepOrb )
   if (allocated(IrpLabOrb ))    deallocate(IrpLabOrb)
   if (allocated(OrbType))       deallocate(OrbType)
   if (allocated(OrbSpace))      deallocate(OrbSpace)
   if (allocated(OrbEnergy))     deallocate(OrbEnergy)
   if (allocated(OrbCoeff))      deallocate(OrbCoeff)
   if (allocated(HCoreMatrix))   deallocate(HCoreMatrix)

   if (allocated(XXPREAD)) deallocate(XXPREAD)
   if (allocated(IXPREAD)) deallocate(IXPREAD)

   return

end subroutine EndGMSMod

end module GMSMod
