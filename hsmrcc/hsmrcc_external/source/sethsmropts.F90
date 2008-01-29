!*****************************************************************
!
!
!   Set of routines called from the C-driver program: hsmrmain.c
!
!
!*****************************************************************

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine SetHSMR_IPRINT(ipr)
!-----------------------------------------------------------------
!
!        ... set up the general control print-out variable 
!
!-----------------------------------------------------------------
  use Parameters
  implicit none

  integer(kind=4), intent(in) :: ipr
  integer :: ipr_temp 
   
  ipr_temp = ipr ! needed proper integer data type 
  call set_iprint (ipr_temp) ! call module routine

end subroutine SetHSMR_IPRINT

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine SetHSMRflags(FlagGENINT , FlagPRNTINF , FlagVERBOSE , &
                        FlagRESTART, FlagTTYRUN  , FlagAVGVAL  , &
                        FlagREPTIME, FlagQUAD    , FlagCONVTAMP, &
                        FlagNOHSMR , FlagMAXIT   , FlagINTTOL  , &
                        FlagCNVTOL , FlagSHFTCORE, FlagSHFTVIRT, &
                        FlagNOAPND , FlagRDGUESS , FlagDIRAC     )
!-----------------------------------------------------------------
!
!  set up of various HSMRCC flags 
!
!-----------------------------------------------------------------

   use Options
   use Parameters

   implicit none
   integer*4, intent(in) :: FlagGENINT , FlagPRNTINF , FlagVERBOSE
   integer*4, intent(in) :: FlagRESTART, FlagTTYRUN  , FlagAVGVAL
   integer*4, intent(in) :: FlagREPTIME, FlagQUAD    , FlagCONVTAMP
   integer*4, intent(in) :: FlagNOHSMR , FlagMAXIT   , FlagINTTOL
   integer*4, intent(in) :: FlagCNVTOL , FlagSHFTCORE, FlagSHFTVIRT
   integer*4, intent(in) :: FlagNOAPND , FlagRDGUESS , FlagDIRAC

   if ( FlagGENINT   == 1 ) OPTION(GENINT)   = .TRUE.
   if ( FlagPRNTINF  == 1 ) OPTION(PRNTINF)  = .TRUE.
   if ( FlagVERBOSE  == 1 ) OPTION(VERBOSE)  = .TRUE.
   if ( FlagRESTART  == 1 ) OPTION(RESTART)  = .TRUE.
   if ( FlagTTYRUN   == 1 ) OPTION(TTYRUN)   = .TRUE.
   if ( FlagAVGVAL   == 1 ) OPTION(AVGVAL)   = .TRUE.
   if ( FlagREPTIME  == 1 ) OPTION(REPTIME)  = .TRUE.
   if ( FlagQUAD     == 1 ) OPTION(QUAD)     = .TRUE.
   if ( FlagCONVTAMP == 1 ) OPTION(CONVTAMP) = .TRUE.
   if ( FlagNOHSMR   == 1 ) OPTION(NOHSMR)   = .TRUE.
   if ( FlagMAXIT    == 1 ) OPTION(MAXIT)    = .TRUE.
   if ( FlagINTTOL   == 1 ) OPTION(INTTOL)   = .TRUE.
   if ( FlagCNVTOL   == 1 ) OPTION(CNVTOL)   = .TRUE.
   if ( FlagSHFTCORE == 1 ) OPTION(SHFTCORE) = .TRUE.
   if ( FlagSHFTVIRT == 1 ) OPTION(SHFTVIRT) = .TRUE.
   if ( FlagNOAPND   == 1 ) OPTION(NOAPND)   = .TRUE.
   if ( FlagRDGUESS  == 1 ) OPTION(RDGUESS)  = .TRUE.
   if ( FlagDIRAC    == 1 ) OPTION(DIRACp)   = .TRUE.  ! ... DIRAC on

   if (get_iprint() > 5) then
      call print_Options 
   endif

   return

end subroutine SetHSMRflags

subroutine SetHSMRpars(ValNCMO,ValNHOLE,ValNPTCL,ValNEMO)
!
!        set important parameters for the code
!
   use Parameters
   implicit none
   integer*4, intent(in) :: ValNCMO,  ValNEMO
   integer*4, intent(in) :: ValNHOLE, ValNPTCL

   NCMO  = ValNCMO
   NHOLE = ValNHOLE
   NPTCL = ValNPTCL
   NEMO  = ValNEMO

   return

end subroutine SetHSMRpars

subroutine SetHSMRshiftpars(SCVal,SVVal,SAVal,BeetaVal)

   use Shift
   implicit none
   double precision, intent(in) :: SCVal
   double precision, intent(in) :: SVVal
   double precision, intent(in) :: SAVal
   double precision, intent(in) :: BeetaVal

   CoreShiftVal = SCVal
   VirtShiftVal = SVVal
   ActvShiftVal = SAVal
   Beeta = BeetaVal

   return

end subroutine SetHSMRshiftpars

subroutine SetHSMRfiles(IFile,OFile)
!----------------------------------------------------------
! 
!       Routine for setting files for HSMRCC run.
!
!----------------------------------------------------------
   use IOModule
   implicit none
   character(len=*), intent(in) :: IFile
   character(len=*), intent(in) :: OFile
   integer :: ilen, olen

   Call QEnter('SetHSMRfiles')

   if (get_iprint() > 5) then
    write(UOUT,*) '>>> SetHSMRfiles called form hsmrmain.c !'
    write(UOUT,*) 'SetHSMRfiles: IFile=',IFile
    write(UOUT,*) 'SetHSMRfiles: OFile=',OFile
   endif

   FileHSMR = REPEAT(' ',LEN(FileHSMR))

   ilen = LEN_TRIM(IFile)
   if (ilen == 0) then
      Write(*,*) 'Input File of Zero Length'
      Write(UOUT,*) 'Input File of Zero Length'
      !call abort
      call Quit('SetHSMRfiles: Input File of Zero Length !')
   endif

   if (TRIM(IFile)==REPEAT(' ',ilen)) then
      !Write(*,*) 'blank Input File'
      Write(UOUT,*) 'blank Input File'
      !call abort
      call Quit('SetHSMRfiles: blank Input File!')
   endif

! =================================================
!       Gamess produced files read by HSMRCC 
! =================================================
   FileMR1EINT  = TRIM(IFile)//'.F295'   ! binary
   FileMR2EINT  = TRIM(IFile)//'.F296'   ! binary
   FileSYMINF   = TRIM(IFile)//'.F297'   ! text file...

   FileVINF  = TRIM(IFile)//'.VINF'

   FileVINT  = TRIM(IFile)//'.VINT'
   FileVIND  = TRIM(IFile)//'.VIND'
   FileFHH   = TRIM(IFile)//'.FHH'
   FileFHP   = TRIM(IFile)//'.FHP'
   FileFPP   = TRIM(IFile)//'.FPP'
   FileTAMP  = TRIM(IFile)//'.TAMP'
   FileTOLD  = TRIM(IFile)//'.TOLD'
   FileTNEW  = TRIM(IFile)//'.TNEW'
   FileHEFF  = TRIM(IFile)//'.HEFF'
   FileCOEF  = TRIM(IFile)//'.COEF'
   FileCNVR  = TRIM(IFile)//'.CNVR'

   olen = LEN_TRIM(OFile)
   if (olen == 0) then
      FileHSMR = TRIM(IFile)//'.HSMR'
   else
      if (TRIM(OFile)==REPEAT(' ',olen)) then
         Write(*,*) 'blank Output File'
         Write(UOUT,*) 'blank Output File'
         !call abort
         call Quit('SetHSMRfiles: blank Output File')
      endif
      FileHSMR = TRIM(OFile)
   endif

   if (get_iprint()>5) then
    write(UOUT,*) 'HSNRCC output file is ',FileHSMR
   endif

   Call QExit('SetHSMRfiles')
   return

end subroutine SetHSMRfiles

subroutine SetHSMRConvTol(Tol)

   use ConvMod
   implicit none
   double precision, intent(in) :: Tol

   CONVTOL = Tol

   return

end subroutine SetHSMRConvTol

subroutine SetHSMRIntTol(Tol)

   use IntGenMod
   implicit none
   double precision, intent(in) :: Tol

   CUTOFF = Tol
   return

end subroutine SetHSMRIntTol

subroutine SetHSMRMaxIt(It)

   use ConvMod
   implicit none
   integer*4, intent(in) :: It

   MAXITER = It

end subroutine SetHSMRMaxIt

