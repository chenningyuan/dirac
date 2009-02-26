!*************************************************************************
!
!
!  Set of routines called from the main C (driver) program: hsmrmain.c
!
!
!*************************************************************************

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine SetHSMR_IPRINT(ipr)
!-----------------------------------------------------------------
!
!     ... set up the general control print-out variable 
!
!  Called from main C driver "hsmrmain.c"
!
!-----------------------------------------------------------------
  use Parameters
  implicit none

  integer(kind=4), intent(in) :: ipr
  integer :: ipr_temp 
   
  ipr_temp = ipr ! needed proper integer data type 
  call set_hscc_iprint (ipr_temp) ! call module routine

end subroutine SetHSMR_IPRINT

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine SetHSMRflags(FlagGENINT , FlagPRNTINF , FlagVERBOSE , &
                        FlagRESTART, FlagTTYRUN  , FlagAVGVAL  , &
                        FlagREPTIME, FlagQUAD    , FlagCONVTAMP, &
                        FlagNOHSMR , FlagMAXIT   , FlagINTTOL  , &
                        FlagCNVTOL , FlagSHFTCORE, FlagSHFTVIRT, &
                        FlagNOAPND , FlagRDGUESS , FlagDIRAC   , &
                        FlagGMSCHK ) 
!------------------------------------------------------------------------
!
!  set up of various HSMRCC flags 
!
!  special note for the DIRAC flag ...
!
!  called as the first F90 routine from main C program in "hsmrmain.c"
!
!------------------------------------------------------------------------

   use Options
   use Parameters

   implicit none
   integer*4, intent(in) :: FlagGENINT , FlagPRNTINF , FlagVERBOSE
   integer*4, intent(in) :: FlagRESTART, FlagTTYRUN  , FlagAVGVAL
   integer*4, intent(in) :: FlagREPTIME, FlagQUAD    , FlagCONVTAMP
   integer*4, intent(in) :: FlagNOHSMR , FlagMAXIT   , FlagINTTOL
   integer*4, intent(in) :: FlagCNVTOL , FlagSHFTCORE, FlagSHFTVIRT
   integer*4, intent(in) :: FlagNOAPND , FlagRDGUESS , FlagDIRAC
   integer*4, intent(in) :: FlagGMSCHK

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

!  ... when DIRAC is on
   if ( FlagDIRAC    == 1 ) then
      OPTION(HSCC_DIRAC_MODE)  = .TRUE. 
      OPTION(GENINT)  = .FALSE.
      write(*,'(/,2x,a,/)') 'SetHSMRflags: DIRAC is ON, therefore OPTION(GENINT) set to false !'
   endif

!  .... when GAMESS integrals testing is on
   if ( FlagGMSCHK   == 1 ) OPTION(GMSCHK)   = .TRUE. 

   if (get_iprint() > 5) then
      call print_Options 
   endif

   if (OPTION(HSCC_DIRAC_MODE)) then
  !  write(UOUT,'(/,2x,a,/)') 'SetHSMRflags: DIRAC is ON!'
     write(*,'(/,2x,a,/)') 'SetHSMRflags 2: DIRAC is ON!'
     call flushout
     !call Quit('SetHSMRflags')
   endif

   return

end subroutine SetHSMRflags

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine SetHSMRpars(ValNCMO,ValNHOLE,ValNPTCL,ValNEMO)
!------------------------------------------------------------------------
!
!    set constants for the code :
!
!                     NCMO
!                     NHOLE
!                     NPCTL
!                     NEMO
!
! TODO: Check their values ...
!
!  Called as the last F90 routine from the main C-program "hsmrmain.c".
!
!------------------------------------------------------------------------
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine SetHSMRfiles(IFile,OFile)
!--------------------------------------------------------------------------------
! 
!  Routine for setting input/output files for HSMRCC run.
!
!  Called as the last F90 routine from the main C-program "hsmrmain.c", also from ..
!
!  NOTE: print level is already set 
!
!--------------------------------------------------------------------------------
   use IOModule
   use Options

   implicit none

   character(len=*), intent(in) :: IFile
   character(len=*), intent(in) :: OFile
   integer :: ilen, olen


   if (OPTION(HSCC_DIRAC_MODE)) then
     Write(*,'(2x,a)') 'SetHSMRfiles: SetHSMRfiles: DIRAC mode is on!'
   else
     Write(*,'(2x,a)') 'SetHSMRfiles: SetHSMRfiles: GAMESS mode is on!'
   endif

   if (get_iprint() >= 5) then
    !write(*,'(/,2x,a)') '>>> SetHSMRfiles called form hsmrmain.c !'
    write(*,'(2x,a,a)') 'SetHSMRfiles: IFile=',IFile
    write(*,'(2x,a,a)') 'SetHSMRfiles: OFile=',OFile
   endif

   if (OPTION(HSCC_DIRAC_MODE)) then
     Write(*,'(2x,a)') 'SetHSMRfiles: DIRAC mod is on!'
     !Call Quit('SetHSMRfiles xx: DIRAC is on!')
   endif

   if (.not. OPTION(HSCC_DIRAC_MODE) ) then

     FileHSMR = REPEAT(' ',LEN(FileHSMR))
     ilen = LEN_TRIM(IFile)

     if (ilen == 0) then
        Write(*,*) 'Input File of Zero Length'
        call Quit('SetHSMRfiles: Input File of Zero Length !')
     endif

     if (TRIM(IFile)==REPEAT(' ',ilen)) then
       Write(*,*) 'SetHSMRfiles:blank Input File'
       !Write(UOUT,*) 'blank Input File'
       call Quit('SetHSMRfiles: blank Input File!')
     endif

!      Gamess produced files read by HSMRCC - read only for Gamess reading
     FileMR1EINT  = TRIM(IFile)//'.F295'   ! binary
     FileMR2EINT  = TRIM(IFile)//'.F296'   ! binary
     FileSYMINF   = TRIM(IFile)//'.F297'   ! text file...

   endif

!---------------------------------------------------------
!  ... the rest are own HSMRCC files
!---------------------------------------------------------

!   ... these file are provided by DIRAC in the case of DIRAC run
   FileVINF  = TRIM(IFile)//'.VINF'
   FileVINT  = TRIM(IFile)//'.VINT'
   FileVIND  = TRIM(IFile)//'.VIND'
   FileFHH   = TRIM(IFile)//'.FHH'
   FileFHP   = TRIM(IFile)//'.FHP'
   FileFPP   = TRIM(IFile)//'.FPP'

!  ... files generated by HSMRCC code, but names must be provided in the case of DIRAC run !
   FileTAMP  = TRIM(IFile)//'.TAMP'
   FileTOLD  = TRIM(IFile)//'.TOLD'
   FileTNEW  = TRIM(IFile)//'.TNEW'
   FileHEFF  = TRIM(IFile)//'.HEFF'
   FileCOEF  = TRIM(IFile)//'.COEF'
   FileCNVR  = TRIM(IFile)//'.CNVR'

if (.not. OPTION(HSCC_DIRAC_MODE)) then

   olen = LEN_TRIM(OFile)
   if (olen == 0) then
      FileHSMR = TRIM(IFile)//'.HSMR'
      if (get_iprint() >= 6) then
          write(*,'(/,2x,a,a)') 'HSMRCC output file is ',FileHSMR
      endif
   else
      if (TRIM(OFile)==REPEAT(' ',olen)) then
         Write(*,*) 'blank Output File'
         Write(*,'(2x,a)') 'SetHSMRfiles: blank Output File'
         call Quit('SetHSMRfiles: blank Output File')
      endif
      FileHSMR = TRIM(OFile)
   endif

else ! ie DIRAC is on
 continue
endif

! ... some output ...
 if (get_iprint() >= 8) then
    !write(*,'(/,2x,a)') '*** SetHSMRfiles: HSMRCC working files '
    call Header ('*** SetHSMRfiles: HSMRCC working files ',-1)

    if (OPTION(HSCC_DIRAC_MODE)) then
     write(*,'(a)') 'first six files are produced by DIRAC:'
    endif

    write(*,'(2x,a,a)') 'FileVINF=',FileVINF
    write(*,'(2x,a,a)') 'FileVINT=',FileVINT
    write(*,'(2x,a,a)') 'FileVIND=',FileVIND
    write(*,'(2x,a,a)') ' FileFHH=',FileFHH
    write(*,'(2x,a,a)') ' FileFHP=',FileFHP
    write(*,'(2x,a,a)') ' FileFPP=',FileFPP

    if (OPTION(HSCC_DIRAC_MODE)) then
     write(*,'(/,a)') 'these files are produced by the HSMRCC module:'
    endif

    write(*,'(2x,a,a)') 'FileTAMP=',FileTAMP
    write(*,'(2x,a,a)') 'FileTOLD=',FileTOLD
    write(*,'(2x,a,a)') 'FileTNEW=',FileTNEW
    write(*,'(2x,a,a)') 'FileTOLD=',FileTOLD
 endif

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
  ! integer*4, intent(in) :: It
   integer, intent(in) :: It

   MAXITER = It

end subroutine SetHSMRMaxIt

subroutine QEnter_C(RoutineName)
!
!  C substitute of DIRAC QENTER
!
character(len=*), intent(in) :: RoutineName
!write(*,*) 'QEnter_C...RoutineName=',RoutineName
Call QEnter(RoutineName)
end subroutine QEnter_C

subroutine QExit_C(RoutineName)
!
!  C substitute of DIRAC QEXIT
!
character(len=*), intent(in) :: RoutineName
!write(*,*) 'QExit_C...RoutineName=',RoutineName
Call QExit(RoutineName)
end subroutine QExit_C

