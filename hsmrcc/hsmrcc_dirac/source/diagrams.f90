!  this subroutine produces TNEW from TOLD
!  by calculating direct, heff and renormalization
!  diagrams of HSMRCC.

!  Notes
!  -----
!
!  Before calculating renormalization
!  terms, copy the singles part of
!  TOLD for current vaccum into TSING
!  this is used for calculating disconnected
!  diagrams in HSMRCC doubles equations. this
!  needs to be copied because TOLD is
!  overwritten with TOLD(tvac) in
!  renormalization routines.

subroutine Diagrams

   use Parameters
   use IOModule
   use TAmps
   use DirectDiags
   use HEFFModule
   use RenormDiags
   use ConvMod
   use TimingMod
   implicit none
   integer :: vac
   ! ... Interface for External Subroutines ...
   interface
      subroutine EnDenom(T,TSize,vac)
         implicit none
         integer, intent(in) :: vac
         integer, intent(in) :: TSize
         double precision, intent(inout) :: T(0:TSize)
      end subroutine EnDenom
      subroutine MakeZero(T,TSize,vac)
         implicit none
         integer, intent(in) :: vac
         integer, intent(in) :: TSize
         double precision, intent(inout) :: T(0:TSize)
      end subroutine MakeZero
      subroutine Compensate(Tr,Ta,TSize,vac)
         implicit none
         integer, intent(in) :: vac
         integer, intent(in) :: TSize
         double precision, intent(inout) :: Tr(0:TSize)
         double precision, intent(in) :: Ta(0:TSize)
      end subroutine Compensate
   end interface

   ! ... Diagram calculation for each vaccum ...

   do vac = 1, NSTATE
      TNEW = dble(0)
      read(UTOLD,rec=vac) TOLD
      call DrvDirect(vac)
      call DrvHEFF(vac)
      TSING = TOLD(1:NHMO*NPMO)
      call DrvRenorm(vac)
      ! ... TOLD would have been modified in DrvRenorm, restore it ...
      read(UTOLD,rec=vac) TOLD
      ! ... TNEW has uncompensated residual ...
      ! ... Apply compensation to residual  ...
      call Compensate(TNEW,TOLD,NTSIZE,vac)
      ! ... Apply shifted energy denominators to compensated residual ...
      call EnDenom(TNEW,NTSIZE,vac)
      call MakeZero(TNEW,NTSIZE,vac)
      TDIFF(vac) = SUM(DABS(TNEW-TOLD))
      write(UTNEW,rec=vac) TNEW
   enddo

   return

end subroutine Diagrams
