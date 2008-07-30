!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine InitHSMRBase
!---------------------------------------------------------
!
!  ... Basic Initialization of HSMR Run ...
!
!      Called from "source/hsmrcc.F90"/HSMRCC
!
!---------------------------------------------------------

   use Parameters
   use Options
   use IOModule
   use VacMod
   use IntGenMod
   use ConvMod
   use Shift
   use VINFMod
   use TimingMod
   use GMSMod
   implicit none
   integer :: SvNCMO, SvNHOLE, SvNPTCL, SvNEMO
   integer :: mxit
   logical :: CFlag, VFlag
   double precision :: Core, Virt
   double precision :: cnv, cut
   integer :: fseek, ierrno
   external :: fseek

   if ( get_iprint() > 5 ) then
     call HEADER('Output from InitHSMRBase',-1)
   endif

   ! ... Start of Timing ...

   call StartTiming

   ! ... Setting up the output stream ...

   if (OPTION(TTYRUN)) then
      UOUT = 6 ! ... terminal output ...
   else
      close(6)
      close(5)
      close(0)
      open(UOUT,file=FileHSMR) ! ... output to the file ...
      if ( get_iprint() >= 15 ) then
          write(UOUT,'(2x,a,i3)') 'InitHSMRBase: UOUT variable set to ',UOUT
          write(UOUT,'(2x,a,a)')  &
          'InitHSMRBase: UOUT assigned file is ',FileHSMR
      endif
      if (OPTION(RESTART).and.(.not.OPTION(NOAPND))) then
         ierrno = fseek(UOUT,0,2)
         if (ierrno/=0) then
           write(UOUT,*) 'InitHSMRBase: fseek error on UOUT'
           call flushout
           call Quit('InitHSMRBase: fseek error on UOUT!')
         endif
      endif
   endif

   ! ... Opening FileVINF file ...

   call OpenVINF

   ! ... Read Options from FileVINF unless it is fresh run ...
   ! ... In case of fresh run, read it from gamess files   ...

   if (.not.OPTION(GENINT)) then
      rewind(UVINF)
      read(UVINF) SvNCMO, SvNHOLE, SvNPTCL, SvNEMO
!
!  .... these variables are in the module GMSmod
!
      NCMO  = SvNCMO
      NEMO  = SvNEMO
      NHOLE = SvNHOLE
      NPTCL = SvNPTCL
      NVMO = NHOLE + NPTCL

      if (OPTION(RESTART)) then
         read(UVINF) OPTPRV
         OPTION(AVGVAL) = OPTPRV(AVGVAL)
         OPTION(QUAD)   = OPTPRV(QUAD)
         read(UVINF) mxit, cnv
         read(UVINF) cut
         if (.not.OPTION(MAXIT))  MAXITER = mxit
         if (.not.OPTION(CNVTOL)) CONVTOL = cnv
         if (.not.OPTION(INTTOL)) CUTOFF  = cut
         read(UVINF) 
         read(UVINF) Core, Virt
         CoreShiftVal = Core
         VirtShiftVal = Virt
      endif
   endif

   ! ... Initialize gamess module, GMSMod ...
   call InitGMSMod

   if ((.not.OPTION(RESTART)).and.(.not.OPTION(INTTOL))) CUTOFF  = TrfTol

   ! ... Generate HSMRCC Parameters ...
   call GenHSMRpars

   ! ... Initialize IOModule   ...
   ! ... Setup NZIntRec Arrays ...

   call InitIOMod

   if (.not.OPTION(GENINT)) call ReadNZIntRecs

   ! ... Initialize VacMod ...

   call InitS2MatEle
   call InitVacMod

   ! ... Open Integral Files ...
   call OpenIntFiles

   return

end subroutine InitHSMRBase
