!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module Parameters
!------------------------------------------------------------------------------
!
!              important data parameters for the HSMRCC run
!
!------------------------------------------------------------------------------

   implicit none

!----------------------------------------------------------------------------
!       number of integers in floating point number (machine dependent)
!----------------------------------------------------------------------------
!
!  ... copied from GAMESS-US
!
!C     THIS IS THE NUMBER OF INTEGERS IN A WORKING PRECISION
!C     FLOATING POINT NUMBER, WHICH IS MOST COMMONLY 64 BITS.
!C     THUS SYSTEMS WITH 64 BIT INTEGERS USE A VALUE OF 1, WHILE
!C          SYSTEMS WITH 32 BIT INTEGERS USE A VALUE OF 2.
!
!----------------------------------------------------------------------------
!     select 32/64-bit code by the preprocessor statement
!----------------------------------------------------------------------------
#if defined (ARCH32BIT)
   integer, parameter :: INTSFP = 2 
#else
   integer, parameter :: INTSFP = 1
#endif

   ! number of logical files and subroutines

   integer, parameter :: NVFIL  = 14
   integer, parameter :: NSUBR  = 17

   integer, parameter :: V2H2P  = 01
   integer, parameter :: VPPPP  = 02
   integer, parameter :: VHHHH  = 03
   integer, parameter :: VHHPP  = 04
   integer, parameter :: V2HPP  = 05
   integer, parameter :: V2PHH  = 06
   integer, parameter :: VPHHP  = 07
   integer, parameter :: VPH2P  = 08
   integer, parameter :: VPHPP  = 09
   integer, parameter :: V2HHP  = 10
   integer, parameter :: VHHHP  = 11
   integer, parameter :: FHH    = 12
   integer, parameter :: FHP    = 13
   integer, parameter :: FPP    = 14

   integer, parameter :: HEFFEC = 15
   integer, parameter :: RENORM = 16
   integer, parameter :: INTGEN = 17

   ! define alpha spin and beta spin

   integer, parameter :: ALPHA_SPIN = +1
   integer, parameter :: BEETA_SPIN = -1

   ! define core, active and virtual spinorbitals

   integer, parameter :: CORE_ORB = +1
   integer, parameter :: ACTV_ORB =  0
   integer, parameter :: VIRT_ORB = -1

   ! define Q, P-M and P-I spaces

   integer, parameter :: Q_SPACE  = -1
   integer, parameter :: PM_SPACE = +1
   integer, parameter :: PI_SPACE =  0

   ! some double precision parameters

   double precision, parameter :: zero = 0.0d00
   double precision, parameter :: one  = 1.0d00
   double precision, parameter :: half = 0.5d00

   double precision, parameter :: NaN    = Z'7FF8000000000000'

  ! ... check it for different compilers !!!

  ! double precision, parameter :: NegInf = Z'FFF0000000000000'
 !  gfortran 4-1
  ! double precision, parameter :: NegInf = -1.0/0.0
 !  gfortran 4-2
   double precision, parameter :: NegInf = transfer(Z'FFF0000000000000', 0d0)

  ! double precision, parameter :: PosInf = Z'7FF0000000000000'
  !double precision, parameter :: PosInf = 1.0/0.0
  double precision, parameter :: PosInf = transfer(Z'7FF0000000000000', 0d0)

   ! HSMRCC parameters

   integer :: NCMO
   integer :: NHOLE
   integer :: NPTCL
   integer :: NEMO
   integer :: NELE

   integer :: NVMO
   integer :: NHMO
   integer :: NPMO
   integer :: NTMO

   integer :: NHMOSQ
   integer :: NPMOSQ
   integer :: NVMOSQ

   integer :: NSING
   integer :: NDBLE
   integer :: NTSIZE

   integer :: MAXEXC

   integer :: NMAINORB
   integer :: NINTMORB
   integer :: NINTMORB_CORE
   integer :: NINTMORB_VIRT

   integer :: NSTATE
   integer :: NSTATE_MAIN
   integer :: NSTATE_INTM

   ! ... control print-out parameter (introduced by MI)
   integer, private :: IPRINT = 0

   ! ... print-out channel - moved here from the module "IOmodule"
   integer :: UOUT  = 51

   ! ------------------------------
   ! ... DIRAC input parameters ...
   ! ------------------------------
   type dirac 
   ! ... flag for the DIRAC input 
     logical :: is = .false.
   ! ... flag for boson/fermion symmetry
     logical :: spinfree = .true.  
   end type dirac
 
  ! private :: dirac

   contains

   subroutine set_iprint (ipr)
   ! ... set the IPRINT variable
    integer, intent(in) :: ipr
    integer, parameter :: ipr_max = 100
    call QEnter('set_iprint')
    if (ipr < 0 .or. ipr > ipr_max) then
     write(UOUT,*) 'set_iprint: entering iprint = ',ipr
     write(UOUT,*) 'set_iprint: wrong entering variable ! must be <0,100>'
     call Quit('set_iprint: wrong entering variable ! must be <0,100>')
    else 
      IPRINT = ipr
    endif
    call QExit('set_iprint')
   end subroutine set_iprint

   function get_iprint() result(ipr)
   ! ... returns the IPRINT variable
     integer :: ipr
     ipr = IPRINT
   end function get_iprint

   subroutine flushout
   ! ... short version for flushing out the buffer
    call flush(UOUT)
   end subroutine flushout

end module Parameters
