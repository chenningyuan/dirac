module hsfs_parameters
! ... to replace param_ccsd.h (the old 'parameters.f')

! ... initialize all to -1
integer :: nbig =  -1
integer :: ndij =  -1
integer :: ndint = -1
integer :: ndt1 = - 1
integer :: ndt2 = -1
integer :: ndp3 = -1
integer :: ndp4 = -1
integer :: ndp5 = -1
integer :: ndp6 = -1
integer :: nrle1 = -1
integer :: nrle2 = -1
integer :: nrle3 = -1
integer :: nrle4 = -1

contains

subroutine assign_hsfs_parametersXXX
#include "priunit.h"
#include "param_ccsd.h"
write(LUPRI,*) 'HSFS parameters from param_ccsd.h...'
 call print_hsfs_parameters
end subroutine assign_hsfs_parametersXXX

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine assign_hsfs_parameters
!---------------------------------------------------------------------------------------------
! very important - assign dimensions of arrays for various allocations
!
! MUST be the very first calling 
!---------------------------------------------------------------------------------------------
#include "priunit.h"

write(LUPRI,'(/,2x,a)') 'assigning HSFS dimensions parameters'

! nbig =  9000000
!  nbig = 20000000
 !nbig =  90000000
 !nbig  = 240000000
 !nbig  = 600000000
 nbig  = 300000000

 ndij = 960000

! ndint =150000
  ndint =15000000

 ndt1 =1000
 ndt2 =500000

! ndp3 = 700000
 !ndp3 = 7000000
  ndp3 = 20000000

! ndp4 =2000000
! ndp4 = 3000000
! ndp4 = 5000000
! ndp4 = 20000000
!  ndp4 = 35000000
  ndp4 = 95000000

 ndp5 = 325000

 ndp6 = ndp3+(3*ndp5)

 nrle1= 8*ndt1
 nrle2 = 8*ndt2 

! nrle3 = 1
 nrle3 = 8*ndp3  ! fix for CCSDT run

 nrle4 = 1000

 call print_hsfs_parameters
  
end subroutine assign_hsfs_parameters

subroutine print_hsfs_parameters
#include "priunit.h"
write(LUPRI,*) 'assigned HSFS parameters...'

write(LUPRI,*) 'nbig=',nbig
write(LUPRI,*) 'ndij=',ndij
write(LUPRI,*) 'ndint=',ndint
write(LUPRI,*) 'ndt1=',ndt1
write(LUPRI,*) 'ndt2=',ndt2
write(LUPRI,*) 'ndp3=',ndp3
write(LUPRI,*) 'ndp4=',ndp4
write(LUPRI,*) 'ndp5=',ndp5
write(LUPRI,*) 'ndp6=',ndp6
write(LUPRI,*) 'nrle1=',nrle1
write(LUPRI,*) 'nrle2=',nrle2
write(LUPRI,*) 'nrle3=',nrle3
write(LUPRI,*) 'nrle4=',nrle4

end subroutine print_hsfs_parameters

end module hsfs_parameters
