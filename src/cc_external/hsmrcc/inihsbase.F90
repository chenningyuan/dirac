!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine InitHSMRBase
!---------------------------------------------------------
!
!     ... initialization for the HSMRCC run ...
!
!    Called from "source/hsmrcc.F90"/HSMRCC
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
  ! external :: fseek

   if ( get_iprint() > 5 ) then
     call HEADER('Output from InitHSMRBase',-1)
   endif

   ! ... Start of Timing ...

   call StartTiming

   ! ... Open FileVINF file ...either new or existing
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

   ! ... set CUTOFF 
   if ((.not.OPTION(RESTART)).and.(.not.OPTION(INTTOL))) CUTOFF  = TrfTol

   ! ... Generate HSMRCC Parameters ...
   call GenHSMRpars

   ! ... Initialize IOModule   ...
   ! ... Setup NZIntRec Arrays ...

   call InitIOMod

   ! ... read 
   if (.not.OPTION(GENINT)) call ReadNZIntRecs

   ! ... Initialize S2MatrixElem...
   call InitS2MatEle

   ! ... Initialize VacMod ...
   call InitVacMod

   ! ... Open Integral Files ...
   call OpenIntFiles

   return

end subroutine InitHSMRBase


subroutine set_output_stream()
!
! ... sets the output stream ...
!    called from source/hsmrcc.F90 HSMRCC
!
use IOModule

   !integer :: fseek, ierrno
   integer(kind=4) :: ierrno ! kind=4 is gfortran demand
   ! external :: fseek

   print '(a,l1)', 'set_output_stream, OPTION(TTYRUN)=',OPTION(TTYRUN)

   if (OPTION(TTYRUN)) then
      UOUT = 6 ! ... terminal output ...
      if ( get_iprint() >= 6 ) then
        write(UOUT,'(/,2x,a,i3,$)') 'set_output_stream: UOUT=',UOUT
        write(UOUT,'(2x,a,i3)')   'print level is ',get_iprint()
      endif
   else
   ! ... output to the files...
      close(6)
      close(5)
      close(0)
      open(UOUT,file=FileHSMR) ! ... output to the file ...
      if ( get_iprint() >= 6 ) then
         write(UOUT,'(2x,a,i3)') ': UOUT variable set to ',UOUT
         write(UOUT,'(2x,a,a)')  &
         ': UOUT assigned file is ',FileHSMR
      endif
      if (OPTION(RESTART).and.(.not.OPTION(NOAPND))) then

         !ierrno = fseek(UOUT,0,2) ! mi: original 
         !call fseek(UOUT,0,2)

        ! call fseek(UOUT,0,2,ierrno) ! GNU/gfortran substitute

         if (ierrno/=0) then
           write(UOUT,*) ': fseek error on UOUT'
           call flushout
           call Quit(': fseek error on UOUT!')
         endif
         ! call quit('fseek deactivated for the current opteron gfortran, which has not fseek')
      endif
   endif

end subroutine set_output_stream
