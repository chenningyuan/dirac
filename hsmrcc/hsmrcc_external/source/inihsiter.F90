subroutine InitHSMRIter

   use IOModule
   use TAmps
   use HEFFModule
   use ConvMod
   use Shift
   use DirectDiags
   implicit none

   ! ... Setup Integral Data Buffers & Cluster Amplitudes ...
   ! ... Initilize Direct Diagram Driver                  ...
   ! ... Open Files for Storing Cluster Amplitues & HEFF  ...

   call AllocIntBufs
   call AllocTAMP
   call InitDirect
   call OpenOtherFiles

   ! ... Setup BlockIO ...

   call InitBlockIO

   ! ... Initialize Shift, HEFFModule & ConvMod Modules ...

   call InitShift
   call InitHEFF
   call InitConvMod

   ! ... Generate Initial Guess ...

   call InitGuess

   return

end subroutine InitHSMRIter
