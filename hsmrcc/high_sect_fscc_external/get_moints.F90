!
! module for getting MO integrals into higher sector FSCC code
!
! here these are coming in from Sham's modified GAMESS
!


 module     get_gamess_mo_ints

   integer, parameter :: LuHCRMO=18
   integer, parameter :: LuMULMO=19

   character(len=5) :: FileHCRMO="HCRMO"
   character(len=5) :: FileMULMO="MULMO"

   integer,allocatable :: IOCCMO(:)
   real*8,allocatable  :: ORBENMO(:)

 contains

 subroutine open_gamess_files
#include "implicit.h"
#include "priunit.h"

  open (unit=LuHCRMO,file=FileHCRMO,form="unformatted",status="old",iostat=ios_hcrmo)
  open (unit=LuMULMO,file=FileHCRMO,form="unformatted",status="old",iostat=ios_mulmo)

  if (ios_hcrmo.eq.0.and.ios_mulmo.eq.0) then
   write(LUPRI,*) 'GAMESS files open'
  else
   write(LUPRI,*) 'GAMESS files not open !'
   call quit('GAMESS files not open !')

  endif

 end subroutine open_gamess_files

 subroutine get_gamess_ints
#include "implicit.h"
#include "priunit.h"

   call open_gamess_files

   call rd_hcrmo

   call quit('get_gamess_ints')

 end subroutine get_gamess_ints

 subroutine rd_hcrmo
  use memory_allocator
#include "implicit.h"
#include "priunit.h"
!C     THE FIRST RECORD CONTAINS INFORMATION                             
!     WRITE(IHCRMO) NMO                                                 
!     WRITE(IHCRMO) IOCCMO                                  
!     WRITE(IHCRMO) ORBENMO          

   read(LuHCRMO) NMO
   write(LUPRI,*) 'read NMO=',NMO

   call alloc(IOCCMO,NMO)
   call alloc(ORBENMO,NMO)

   read(LuHCRMO) IOCCMO
   read(LuHCRMO) ORBENMO

   write(LUPRI,*) 'IOCCMO:', (IOCCMO(I),I=1,NMO)
   write(LUPRI,*) 'ORBENMO:',(ORBENMO(I),I=1,NMO)

 end subroutine rd_hcrmo 

 end module get_gamess_mo_ints


