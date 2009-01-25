!> important data parameters for the HSMRCC run
module Parameters

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
!
! TODO: generalize for ILP32, I32LP64 and ILP64 data types !
!
#if defined (ARCH32BIT)
    !> ...32 bit architecture on
   integer, parameter :: INTSFP = 2 
#else
    !> ... 64-bit architecture on (x86_64) 
   integer, parameter :: INTSFP = 1 
#endif

   ! number of logical files and subroutines

   integer, parameter :: NVFIL  = 14
   integer, parameter :: NSUBR  = 17

    !> V2H2P label index ...
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

     !> zero number
   double precision, parameter :: zero = 0.0d00

   double precision, parameter :: one  = 1.0d00
   double precision, parameter :: half = 0.5d00

!  double precision, parameter :: NaN    = Z'7FF8000000000000'
   double precision, parameter :: NaN    = transfer(Z'7FF8000000000000',0d0)


  ! double precision, parameter :: NegInf = Z'FFF0000000000000'
 !  gfortran 4-1
  ! double precision, parameter :: NegInf = -1.0/0.0
 !  gfortran 4-2
   double precision, parameter :: NegInf = transfer(Z'FFF0000000000000', 0d0)

  ! double precision, parameter :: PosInf = Z'7FF0000000000000'
  !double precision, parameter :: PosInf = 1.0/0.0
  double precision, parameter :: PosInf = transfer(Z'7FF0000000000000', 0d0)

   ! HSMRCC parameters

   integer :: NCMO=0
   !> number of holes
   integer :: NHOLE=0
   !> number of particled
   integer :: NPTCL=0
   integer :: NEMO=0

 !> number of active electrons ...
   integer :: NELE=0

   integer :: NVMO=0
   integer :: NHMO=0
   integer :: NPMO=0
   integer :: NTMO=0

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

   !> ... control print-out parameter (introduced by MI)
   integer, private :: IPRINT = 0

   !> ... print-out channel - moved here from the module "IOmodule"
   integer :: UOUT  = 51

   !> ... DIRAC mode parameters ...
   type dirac 
   !> ... logical flag for the DIRAC mode
     logical :: is = .false.
   !> ... is boson/fermion symmetry
     logical :: spinfree = .true.  
   !> ... is complex algebra (i.e. when complex relativ.group used)
     logical :: is_complex = .false.  
   end type dirac

   !> ... is dirac-like print-out channel
   logical :: dirac_UOUT=.false.

   !> ... is dirac_hscc_interface on
   logical, private :: dirac_hscc_interface_mode=.false.
 
   contains

   !> set print-out channel
   subroutine set_dirac_UOUT
#include "priunit.h"
   UOUT=LUPRI
   write(UOUT,'(/,2x,a,i2)') 'UOUT set according to DIRAC LUPRI to value ',UOUT
   dirac_UOUT=.true.
   end subroutine set_dirac_UOUT

  !> ... tell we are in the Dirac_HSMRCC interface mode, where we are reusing a lot of HSMRCC stuff
   subroutine set_dirac_hscc_interface_mode
     dirac_hscc_interface_mode=.true.
   end subroutine set_dirac_hscc_interface_mode

   subroutine unset_dirac_hscc_interface_mode
     dirac_hscc_interface_mode=.false.
   end subroutine unset_dirac_hscc_interface_mode

   !> ... returns the protected dirac_hscc_interface_mode
   function get_dirac_hscc_interface_mode() result(is_mode)
     logical :: is_mode
     is_mode = dirac_hscc_interface_mode
   end function get_dirac_hscc_interface_mode

   subroutine  check_dirac_hscc_interface_mode
    if (.not.get_dirac_hscc_interface_mode()) then
      call quit('check_dirac_hscc_interface_mode: NOT in the DIRAC_HSMRCC_interface mode !')
    endif
   end  subroutine  check_dirac_hscc_interface_mode

   !> ... set the IPRINT variable
   subroutine set_hscc_iprint (ipr)
   !>  entering general printing variable (between 0-100, veryfied inside)
    integer, intent(in) :: ipr
    integer, parameter :: ipr_max = 100
    if (ipr < 0 .or. ipr > ipr_max) then
     write(UOUT,*) 'set_hscc_iprint: entering iprint = ',ipr
     write(UOUT,*) 'set_hscc_iprint: wrong entering variable ! must be <0,100>'
     call quit('set_hscc_iprint: wrong entering variable ! must be <0,100>')
    else 
      IPRINT = ipr
    endif
    print *,'set_hscc_iprint: print parameter variable set to IPRINT=',IPRINT
    !call quit('set_hscc_iprint') ! test crash ...
   end subroutine set_hscc_iprint

   !> ... returns the protected IPRINT variable
   function get_iprint() result(ipr)
     integer :: ipr
     ipr = IPRINT
   end function get_iprint

   !> ... short version of flushing out the buffer 
   subroutine flushout
      call flush(UOUT)
   end subroutine flushout

!include "tind.h"
function TIND(ij)
!--------------------------------------------------------------------------------------
! ... inline integer function inserted in various places of the code
! ... depends on module variables NPMOSQ, NSING stored in module Parameters
!--------------------------------------------------------------------------------------
   integer, intent(in) :: ij
   integer :: TIND
   TIND = (ij-1)*NPMOSQ + NSING
   return
end function TIND

subroutine check_FileRec_param
!
! check variables for dimensions for NumRec/RecLen arrays
!

! NHMO, NPMO, NVMO, NHMOSQ, NPMOSQ
if (NHMO==0.or.NPMO==0.or.NVMO==0.or.NHMOSQ==0.or.NPMOSQ==0) then
    write(UOUT,'(/,5x,a)') 'check_FileRec_param: some parameters are zero !'
   write(UOUT,*) 'NHMO=',NHMO
   write(UOUT,*) 'NPMO=',NPMO
   write(UOUT,*) 'NVMO=',NVMO
   write(UOUT,*) 'NHMSQ=',NHMOSQ
   write(UOUT,*) 'NPMOSQ=',NPMOSQ
   write(UOUT,*) 
   call quit('check_FileRec_param : some NumRec param. are zero !')
endif
end subroutine check_FileRec_param

subroutine print_basic_params
! ... print basic parameters
  write(UOUT,'(/,8x,a,/)') '**** (Parameters) basis parameters: ****'
  write(UOUT,'(2x,a,i4)') ' NCMO=',NCMO
  write(UOUT,'(2x,a,i4)') ' NVMO=',NVMO
  write(UOUT,'(2x,a,i4)') 'NHOLE=',NHOLE
  write(UOUT,'(2x,a,i4)') ' NEMO=',NEMO
  write(UOUT,'(2x,a,i4)') 'NPTCL=',NPTCL
  call flushout
end subroutine print_basic_params

subroutine print_gms_parameters
! ....
    write(UOUT,'(/,8x,a,/)') '**** (GenGSMRpars) calculated parameters: ****'

    write(UOUT,'(2x,a,i7)') '       NSTATE=',NSTATE
    write(UOUT,'(2x,a,i7)') '  NSTATE_MAIN=',NSTATE_MAIN
    write(UOUT,'(2x,a,i7)') '  NSTATE_INTM=',NSTATE_INTM
    write(UOUT,'(2x,a,i7)') '       NAXEXC=',MAXEXC

    write(UOUT,'(2x,a,i7)') '     NMAINORB=',NMAINORB
    write(UOUT,'(2x,a,i7)') '     NINTMORB=',NINTMORB

    write(UOUT,'(2x,a,i7)') 'NINTMORB_CORE=',NINTMORB_CORE
    write(UOUT,'(2x,a,i7)') 'NINTMORB_VIRT=',NINTMORB_VIRT

    write(UOUT,'(2x,a,i7)') 'NHMO=',NHMO
    write(UOUT,'(2x,a,i7)') 'NPMO=',NPMO
    write(UOUT,'(2x,a,i7)') 'NELE=',NELE
    write(UOUT,'(2x,a,i7)') 'NTMO=',NTMO
    write(UOUT,'(2x,a,i7)') 'NVMO=',NVMO

end subroutine print_gms_parameters

end module Parameters
