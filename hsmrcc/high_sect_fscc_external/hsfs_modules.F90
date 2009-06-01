!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module print_level

 integer, private :: iprint=0
 logical, private :: print_level_read = .false.

contains

subroutine extract_print_level
! extract print level from common block variabla
#include "dat.h"
 iprint = iopt(1)
 print_level_read = .true.! print variable is extracted
end subroutine extract_print_level

function getiprint() result(ipr)
! get private module value
  integer :: ipr

 if (.not.print_level_read) then
   call quit('getiprint: print variable NOT extracted !')
 else
  ipr = iprint
  return
 endif

end function getiprint

end module print_level


module high_sect_execs

integer :: num_methods = 0

character(len=3), allocatable :: method_name(:)

character(len=3):: hs_method_name(3) = (/'csd','csu','csv'/)

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_line_arguments
!------------------------------------------------------------------
!
!    get arguments for order of methods 
!
!------------------------------------------------------------------
#include "priunit.h"
integer :: nargs,i
character(len=3):: this_arg
nargs = command_argument_count()
if (nargs == 0) then
  call quit('no argument !')
else
 allocate (method_name(nargs))
 do i=1,nargs
    CALL get_command_argument(i, this_arg)
    if  (.not.check_method_name(this_arg)) then
      write(LUPRI,*) i,'.argument - method',this_arg, ' not found in the list !'
      call PRTAB(size(hs_method_name),hs_method_name,this_arg//'input keywords',LUPRI)
      call quit('Illegal keyword in HSFS argument list ')
    endif
    method_name(i) = this_arg
 enddo
endif
num_methods=nargs
end subroutine get_line_arguments

function check_method_name(what_method) result (found_method)
!--------------------------------------------------------------
!             checks method name
!--------------------------------------------------------------
character(len=3), intent(in) :: what_method
integer :: i
logical :: found_method

found_method=.false.
do i=1,size(hs_method_name)
 if (what_method==hs_method_name(i)) then
   found_method=.true.
   Exit
 endif
enddo

end function check_method_name


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine launch_hsfs_methods
!--------------------------------------------------------------------
!
!      launches various HSFS methods according to array list
!
!--------------------------------------------------------------------
#include "priunit.h"
integer :: i

write(LUPRI,*) 'num_methods=',num_methods

do i=1, num_methods

  write(LUPRI,*) 'going to launch:', method_name(i)

  if      ( method_name(i) == 'csd' ) then
      call csd
  else if ( method_name(i) == 'csu' ) then
      call csu
  else if ( method_name(i) == 'csv' ) then
      call csv
  else
      write(LUPRI,*) method_name(i),' is not in the list !'
      call quit('launch_hsfs_methods')
  endif

enddo

end subroutine launch_hsfs_methods

end module high_sect_execs


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


