!
! ===  Modules to replace individual common blocks with dynamical allocations &
! other object oriented functionalities
!

module hsfs_parameters
end module hsfs_parameters

module gen_mod
! to replace gen.h
real*8, allocatable :: gen(:)
integer, allocatable :: igen(:)
       
contains

subroutine init_gen_mod
#include "param_ccsd.h"
integer :: istat
 allocate(gen(nbig),stat=istat)
 allocate(igen(nbig),stat=istat)
end subroutine init_gen_mod

end module gen_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module cint_mod
! to replace cint.h
  integer:: nint,idum
  integer, allocatable :: iint(:)
  real*8, allocatable ::  vint(:)
contains
subroutine init_cint_mod
#include "param_ccsd.h"
integer :: istat
 allocate(iint(ndint),stat=istat)
 call alloc_stat_check("iint(ndint)",istat)
 allocate(vint(ndint),stat=istat)
 call alloc_stat_check("vint(ndint)",istat)
end subroutine init_cint_mod

end module cint_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module fock_mod
! to replace fock.h
 integer :: nfock,idum3
 integer, allocatable :: ifock(:)
 real*8,  allocatable :: fock(:)
contains
subroutine init_fock_mod
 allocate(ifock(ndt1),stat=istat)
 call alloc_stat_check("ifock(ndt1)",istat)
 allocate(fock(ndt1),stat=istat)
 call alloc_stat_check("fock(ndt1)",istat)
end subroutine init_fock_mod

end module fock_mod

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

character(len=4), allocatable :: method_name(:)

character(len=4):: hs_method_name(4) = (/'csd','csu','csv','sixh' /)

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_line_arguments
!------------------------------------------------------------------
!
!    get arguments for order of methods 
!
! NOTE: be aware not all F90 are supporting argument
!      line reading intrinsic function !
!
!------------------------------------------------------------------
#include "priunit.h"
integer :: nargs,i
character(len=4):: this_arg

! get number of arguments
!nargs = command_argument_count()
nargs = iargc()

if (nargs == 0) then
  call quit('no argument !')
else
 allocate (method_name(nargs))
 do i=1,nargs

   ! CALL get_command_argument(i, this_arg)
    CALL getarg(i, this_arg)

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
character(len=4), intent(in) :: what_method
integer :: i
logical :: found_method

found_method=.false.
do i=1,size(hs_method_name)
 if (what_method==hs_method_name(i)) then
   found_method=.true.
   Exit
 endif
enddo

! mi: now skip checking
found_method=.true.
return

end function check_method_name


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine launch_hsfs_methods
!--------------------------------------------------------------------
!
!       launches various HSFS methods according to array list   
!
!--------------------------------------------------------------------
#include "priunit.h"
integer :: i

write(LUPRI,*) 'num_methods=',num_methods

do i=1, num_methods

  call header('going to launch higher sector method:',-1) 
  write(LUPRI,'(4x,a)') method_name(i)

  if      ( method_name(i) == 'csd' ) then
      call csd
  else if ( method_name(i) == 'csv' ) then
      call csv
  else if ( method_name(i) == 'csu' ) then
      call csu
  else if ( method_name(i) == 'csw' ) then
      call csw
  else if ( method_name(i) == 'csx' ) then
      call csx
  else if ( method_name(i) == 'csy' ) then
      call csy
  else if ( method_name(i) == 'csz' ) then
      call csz
  else if ( method_name(i) == 'chv' ) then
      call chv
  else if ( method_name(i) == 'chu' ) then
      call chu
  else if ( method_name(i) == 'chw' ) then
      call chw
  else if ( method_name(i) == 'chx' ) then
      call chx
  else if ( method_name(i) == 'chy' ) then
      call chy
  else if ( method_name(i) == 'chz' ) then
      call chz
  else if ( method_name(i) == 'chp' ) then
      call chp
  else if ( method_name(i) == 'chq' ) then
      call chq
  else if ( method_name(i) == 'cst' ) then
      call cst
  else if ( method_name(i) == 'cht' ) then
      call cht
  else if ( method_name(i) == 'sixh' ) then
      call sixh
  else if ( method_name(i) == 'sixp' ) then
      call sixp
  else if ( method_name(i) == 'sixv' ) then
      call sixv
  else if ( method_name(i) == 'm12b' ) then
      call m12b
  else if ( method_name(i) == 'm21b' ) then
      call m21b
  else
      write(LUPRI,*) method_name(i),' is not in the list !'
      call quit('launch_hsfs_methods')
  endif

enddo

end subroutine launch_hsfs_methods

end module high_sect_execs


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 module     get_gamess_mo_ints

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

