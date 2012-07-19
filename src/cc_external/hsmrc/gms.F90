!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module  get_gamess_mo_ints
!
! module for getting MO integrals into higher sector FSCC code
!
! here these are coming in from Sham's modified GAMESS
!
   integer, parameter :: LuHCRMO=18
   integer, parameter :: LuMULMO=19

   character(len=5) :: FileHCRMO="HCRMO"
   character(len=5) :: FileMULMO="MULMO"

   integer,allocatable :: IOCCMO(:)
   real*8,allocatable  :: ORBENMO(:)
   real*8,allocatable  :: HCOREMO(:,:) 

 contains

 subroutine open_gamess_files
#include "implicit.h"
#include "priunit.h"

  open (unit=LuHCRMO,file=FileHCRMO,form="unformatted",status="old",iostat=ios_hcrmo)
  open (unit=LuMULMO,file=FileMULMO,form="unformatted",status="old",iostat=ios_mulmo)

  if (ios_hcrmo.eq.0.and.ios_mulmo.eq.0) then
   write(LUPRI,*) 'GAMESS  produced files open'
  else
   write(LUPRI,*) 'GAMESS files not open !'
   call quit('GAMESS files not open !')
  endif

 end subroutine open_gamess_files

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine get_gamess_ints(nort,eps,nn,ind,v,nint)
!
!    provide data for high-sect CCSD from GAMESS 
!  
! called from readv(ccsd3.F)
!
#include "implicit.h"
#include "priunit.h"
   real*8, intent(inout) :: eps(*), v(*)
   integer, intent(inout) :: ind(*),nint, nort

   call open_gamess_files

   call rd_hcrmo(eps, nort)
   call rd_mulmo(v,ind,nint)

   !call quit('get_gamess_ints')

 end subroutine get_gamess_ints

 subroutine rd_hcrmo(eps,nort)
  use memory_allocator ! use Andre's stuff
#include "implicit.h"
#include "priunit.h"
   real*8, intent(inout) :: eps(*)
   integer, intent(inout) :: nort


!C     THE FIRST RECORD CONTAINS INFORMATION                             
!     WRITE(IHCRMO) NMO                                                 
!     WRITE(IHCRMO) IOCCMO                                  
!     WRITE(IHCRMO) ORBENMO          

   read(LuHCRMO) NMO

   write(LUPRI,*) 'read NMO=',NMO

   call alloc(IOCCMO,NMO,id="IOCCMO")
   call alloc(ORBENMO,NMO,id="ORBENMO")
   call alloc(HCOREMO,NMO,NMO,id="HCOREMO")

   read(LuHCRMO) IOCCMO
   read(LuHCRMO) ORBENMO

  ! if (iopt(1).ge.5) then
     write(LUPRI,*) 'IOCCMO:', (IOCCMO(I),I=1,NMO)
     write(LUPRI,*) 'ORBENMO:',(ORBENMO(I),I=1,NMO)
  ! endif

 ! now read Hcore
 ! WRITE(IHCRMO) HCOREMO

   read(LuHCRMO) HCOREMO

  ! if (iopt(1).ge.5) then
     write(LUPRI,'(/,4x,a,2i3)') 'Hcore matrix',NMO,NMO
     call OUTPUT(HCOREMO,1,NMO,1,NMO,NMO,NMO,1,LUPRI)
  ! endif

   nort = NMO
   call dcopy (NMO,ORBENMO,1,eps,1)
 end subroutine rd_hcrmo 

 subroutine rd_mulmo(v,ind,nint)
  use memory_allocator
  use Parameters ! from HSMRCC
  use IntGenMod ! from SHam's HSMRCC
#include "implicit.h"
#include "priunit.h"
   real*8, intent(inout) :: v(*)
   integer, intent(inout) :: ind(*),nint
!  C     THE FIRST RECORD CONTAINS INFORMATION                             
!      WRITE(IMULMO) NMO, IUNIQ, LABELSIZE, NINTMX
!      WRITE(IMULMO) IOCCMO
!      WRITE(IMULMO) ORBENMO
  rewind(LuMULMO)
  read(LuMULMO,iostat=ios) NMO, IUNIQ, LABELSIZE, NINTMX
  if (ios.ne.0) then
    write(LUPRI,'(/,2x,a)') 'error in MULMO reading of NMO, IUNIQ, LABELSIZE, NINTMX '
    call quit('error in MULMO reading of NMO, IUNIQ, LABELSIZE, NINTMX ')
  endif
  read(LuMULMO) ! skip IOCCMO
  read(LuMULMO) ! skip ORBENMO
  write(LUPRI,*) 'NMO=',NMO
  write(LUPRI,*) 'IUNIQ=',IUNIQ
  write(LUPRI,*) 'LABELSIZE=',LABELSIZE
  write(LUPRI,*) 'NINTMX=',NINTMX
! ... continue with reading of integrals&indexes
 ! iprint=iopt(1); 
  !call set_hscc_iprint(50)
  call set_hscc_iprint(1)
 ! ... routine in HSMRCC directories !
  call get_gms_spacial_integrals(LuMULMO,IUNIQ,LABELSIZE,NINTMX,v,ind,nint)
  !call get_gms_mo_integrals(LuMULMO,IUNIQ,LABELSIZE,NINTMX)
  !call get_gms_mo_integrals
 end subroutine rd_mulmo 
end module get_gamess_mo_ints
