!================
!   File_Units
!================
Module File_Units
  Integer :: FUPrint, FUDal, FUError, FUHSA, FUPhaseSpc, FileStatus
End Module File_Units

!--------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module Dirac_HSMRCC
!
! contains all data from DIRAC for standalone HSMRCC code 
!
integer (kind=2), parameter :: ludh = 12 ! lu unit of formatted DIRAC_for_HSMRCC file
integer (kind=4) ::  NTMO=0, NCMO=0, NVMO=0, NEMO=0, NHOLE=0, NPTCL=0, NELE=0

integer(kind=4) :: iprint_hscc

contains

subroutine open_dirac_HSMRCC_file
open(ludh,file='Dirac_for_HSMRCC',status='formatted', access='sequential') 
end subroutine open_dirac_HSMRCC_file

end module Dirac_HSMRCC

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module io_moltra

! ... logical numbers ... 
integer (kind=2) :: MDINT=20, MRCONEE=21

CHARACTER*10 DATEX,TIMEX*8
real*8 :: ECORE

!integer (kind=4) :: NSP
integer  :: NSP, NFSYM, NZ, NORBT, NSYMRP, NELEC(16), NREP 
CHARACTER*14 REPNT(16)

integer, parameter :: MXREP=32,MAXSP=1000
character*4 REPN(16), REPNA(MXREP*2)

integer :: IREPSP(MAXSP*4),MULTB(MXREP*2,MXREP*2,2)
real*8 EPS1(MAXSP)

LOGICAL BREIT,SPFR

end module io_moltra

