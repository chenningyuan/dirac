!
! ===  Modules to replace individual common blocks with dynamical allocations &
!     other object oriented functionalities
!

module gen_mod
! to replace gen.h
real*8, allocatable :: gen(:)
integer, allocatable :: igen(:)
! note: igen is storing last 3/4/5 indexes, the 1st(2nd,3rd) are stored in ij[] array
contains
subroutine init_gen_mod
 use hsfs_parameters
implicit none
#include "priunit.h"
integer :: istat
 allocate(gen(nbig),stat=istat)
 call alloc_stat_check("gen(nbig)",istat)
 allocate(igen(nbig),stat=istat)
 call alloc_stat_check("igen(nbig)",istat)
 write(LUPRI,*) "gen_mod allocated (gen,igen), size=",nbig
end subroutine init_gen_mod

subroutine dealloc_gen_mod
 deallocate(gen)
 deallocate(igen)
end subroutine dealloc_gen_mod

subroutine  check_gen_mod
 use hsfs_parameters
implicit none
#include "priunit.h"
integer :: isize_gen, isize_igen, itemp1, itemp2
real*8 :: temp1,temp2
if (allocated(gen).and.allocated(igen)) then
isize_gen=size(gen); isize_igen=size(igen)
temp1=gen(1);temp2=gen(nbig)
itemp1=igen(1);itemp2=igen(nbig)
endif
end subroutine  check_gen_mod

end module gen_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module cint_mod
! to replace cint.h
  integer:: nint,idum
  integer, allocatable :: iint(:)
  real*8, allocatable ::  vint(:)
contains
subroutine init_cint_mod
 use hsfs_parameters
#include "priunit.h"
!#include "param_ccsd.h"
integer :: istat
 allocate(iint(ndint),stat=istat)
 call alloc_stat_check("iint(ndint)",istat)
 allocate(vint(ndint),stat=istat)
 call alloc_stat_check("vint(ndint)",istat)
 write(LUPRI,*) 'cint_mod initialized with the parameter ndint=',ndint
end subroutine init_cint_mod

subroutine dealloc_cint_mod
 deallocate(iint)
 deallocate(vint)
end subroutine dealloc_cint_mod
end module cint_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module fock_mod
! to replace fock.h
 integer :: nfock,idum3
 integer, allocatable :: ifock(:)
 real*8,  allocatable :: fock(:)
contains
subroutine init_fock_mod
 use hsfs_parameters
!#include "param_ccsd.h"
#include "priunit.h"
integer :: istat
 allocate(ifock(ndt1),stat=istat)
 call alloc_stat_check("ifock(ndt1)",istat)
 allocate(fock(ndt1),stat=istat)
 call alloc_stat_check("fock(ndt1)",istat)
 write(LUPRI,*) 'fock_mod initialized'
end subroutine init_fock_mod

subroutine dealloc_fock_mod
 deallocate(ifock)
 deallocate(fock)
end subroutine dealloc_fock_mod

end module fock_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module cfiles_mod
! to replace "cfiles.h"
 integer, parameter :: cfiles_size=200 ! to be changed
 integer :: nflcur
 integer, allocatable :: namefl(:),nstrt(:,:),iiia(:),iiva(:)
contains
subroutine init_cfiles_mod
 use hsfs_parameters
#include "priunit.h"
integer :: istat
 allocate(namefl(cfiles_size),stat=istat)
 call alloc_stat_check("namefl(cfiles_size)",istat)
 allocate(nstrt(2,0:cfiles_size),stat=istat)
 call alloc_stat_check("nstrt(2,0:cfiles_size",istat)
 allocate(iiia(cfiles_size),stat=istat)
 call alloc_stat_check("iiia(cfiles_size)",istat)
 allocate(iiva(cfiles_size),stat=istat)
 call alloc_stat_check("iiva(cfiles_size)",istat) 

 write(LUPRI,*) 'cfiles_mod initialized'

end subroutine init_cfiles_mod

subroutine dealloc_cfiles_mod
 deallocate(namefl)
 deallocate(nstrt)
 deallocate(iiia)
 deallocate(iiva)
end subroutine dealloc_cfiles_mod
end module cfiles_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module cij1_mod
! ... to replace cij1.h
integer :: nij1,iii1,iiv1
integer, allocatable ::  ij1(:),ije1(:)
! ij1-stores (packed) first indexes,
! ije1 - stored information about ending of records
contains
subroutine init_cij1_mod
 use hsfs_parameters
#include "priunit.h"
!#include "param_ccsd.h"
 allocate(ij1(1:ndij))
 allocate(ije1(0:ndij))
 write(LUPRI,*) 'cij1_mod initialized'
end subroutine init_cij1_mod

subroutine check_cij1_mod
 use hsfs_parameters
implicit none
#include "priunit.h"
integer ::  itemp1,itemp2,itemp3,itemp4
itemp1=ij1(1);itemp2=ij1(ndij)
itemp3=ije1(0);itemp4=ije1(ndij)
end subroutine check_cij1_mod

subroutine dealloc_cij1_mod
 deallocate(ij1)
 deallocate(ije1)
end subroutine dealloc_cij1_mod
end module cij1_mod

module cij2_mod
! ... to replace cij2.h
integer :: nij2,iii2,iiv2
integer, allocatable ::  ij2(:),ije2(:)
contains
subroutine init_cij2_mod
 use hsfs_parameters
#include "priunit.h"
!#include "param_ccsd.h"
 allocate(ij2(1:ndij))
 allocate(ije2(0:ndij))
 write(LUPRI,*) 'cij2_mod initialized'
end subroutine init_cij2_mod

subroutine dealloc_cij2_mod
 deallocate(ij2)
 deallocate(ije2)
end subroutine dealloc_cij2_mod
end module cij2_mod

module cij3_mod
! ... to replace cij3.h
integer :: nij3,iii3,iiv3
integer, allocatable ::  ij3(:),ije3(:)
contains
subroutine init_cij3_mod
 use hsfs_parameters
#include "priunit.h"
!#include "param_ccsd.h"
 allocate(ij3(1:ndij))
 allocate(ije3(0:ndij))
 write(LUPRI,*) 'cij3_mod initialized'
end subroutine init_cij3_mod

subroutine dealloc_cij3_mod
 deallocate(ij3)
 deallocate(ije3)
end subroutine dealloc_cij3_mod
end module cij3_mod

module cij4_mod
! ... to replace cij4.h
integer :: nij4,iii4,iiv4
integer, allocatable :: ij4(:),ije4(:)
contains
subroutine init_cij4_mod
 use hsfs_parameters
#include "priunit.h"
!#include "param_ccsd.h"
 allocate(ij4(1:ndij))
 allocate(ije4(0:ndij))
 write(LUPRI,*) 'cij4_mod initialized'
end subroutine init_cij4_mod

subroutine dealloc_cij4_mod
 deallocate(ij4)
 deallocate(ije4)
end subroutine dealloc_cij4_mod
end module cij4_mod

module cij5_mod
! ... to replace cij5.h
integer :: nij5,iii5,iiv5
integer, allocatable :: ij5(:),ije5(:)
contains
subroutine init_cij5_mod
 use hsfs_parameters
#include "priunit.h"
!#include "param_ccsd.h"
 allocate(ij5(1:ndij))
 allocate(ije5(0:ndij))
 write(LUPRI,*) 'cij5_mod initialized'
end subroutine init_cij5_mod

subroutine dealloc_cij5_mod
 deallocate(ij5)
 deallocate(ije5)
end subroutine dealloc_cij5_mod
end module cij5_mod

module cij6_mod
! ... to replace cij6.h
integer :: nij6,iii6,iiv6
integer, allocatable :: ij6(:),ije6(:)
contains
subroutine init_cij6_mod
 use hsfs_parameters
#include "priunit.h"
!#include "param_ccsd.h"
 allocate(ij6(1:ndij))
 allocate(ije6(0:ndij))
 write(LUPRI,*) 'cij6_mod initialized'
end subroutine init_cij6_mod

subroutine dealloc_cij6_mod
 deallocate(ij6)
 deallocate(ije6)
end subroutine dealloc_cij6_mod
end module cij6_mod

module rsp1_mod
! .. to replace "rsp1.h"
integer :: nrsp1
integer, allocatable ::  irsp1(:)
contains
subroutine init_rsp1_mod
 use hsfs_parameters
#include "priunit.h"
integer, parameter :: rsp1_size=100
allocate(irsp1(rsp1_size))
 write(LUPRI,*) 'rsp1_mod initialized'
end subroutine init_rsp1_mod

subroutine dealloc_rsp1_mod
 deallocate(irsp1)
end subroutine dealloc_rsp1_mod
end module rsp1_mod

module rsp2_mod
! .. to replace "rsp2.h"
integer :: nrsp2
integer, allocatable ::  irsp2(:)
contains
subroutine init_rsp2_mod
 use hsfs_parameters
#include "priunit.h"
integer, parameter :: rsp2_size=500
allocate(irsp2(rsp2_size))
 write(LUPRI,*) 'rsp2_mod initialized'
end subroutine init_rsp2_mod

subroutine dealloc_rsp2_mod
 deallocate(irsp2)
end subroutine dealloc_rsp2_mod
end module rsp2_mod

module dpd_mod
! .. to replace dpd.h
integer :: nir
integer, allocatable :: igmt(:,:),ici(:)
logical :: dpd_mod_initialized = .false.
contains
subroutine init_dpd_mod
 use hsfs_parameters
#include "priunit.h"
integer,parameter :: dpd_size=200
if (.not.dpd_mod_initialized) then
allocate(ici(1:dpd_size))
allocate(igmt(1:8,1:8))
 write(LUPRI,*) 'dpd_mod initialized'
 dpd_mod_initialized = .true.
else 
 write(LUPRI,*) 'dpd_mod was already initialized, leaving'
endif
end subroutine init_dpd_mod

subroutine dealloc_dpd_mod
 deallocate(ici)
 deallocate(igmt)
end subroutine dealloc_dpd_mod
end module dpd_mod

module ctmp_mod
! .. to replace ctmp.h
real(kind=8), allocatable  ::    tmp(:)
     integer, allocatable  ::   itmp(:)
contains
subroutine init_ctmp_mod
 use hsfs_parameters
#include "priunit.h"
!#include "param_ccsd.h"
allocate(tmp(1:ndp4))
allocate(itmp(1:ndp4))
write(LUPRI,*) 'ctmp_mod initialized, ndp4=',ndp4
!write(LUPRI,*) 'ctmp_mod initialized, itmp:',itmp(1),itmp(2),itmp(ndp4-1),itmp(ndp4)
!write(LUPRI,*) 'ctmp_mod: size(itmp)=',size(itmp)
end subroutine init_ctmp_mod

subroutine check_ctmp_mod
! for run-time check...
 use hsfs_parameters
#include "priunit.h"
if (allocated(tmp).and.allocated(itmp)) then
 temp1=tmp(1);temp2=tmp(ndp4);
 if (size(tmp).ne.ndp4) call quit('tmp_mod: size <> ndp4!')
 itemp1=itmp(1);itemp2=itmp(ndp4);
 if (size(itmp).ne.ndp4) call quit('size(itmp) <> ndp4!')
else
 call FLSHFO(LUPRI)
 call quit('ctmp_mod not allocated!')
endif
end subroutine check_ctmp_mod

subroutine dealloc_ctmp_mod
 deallocate(tmp)
 deallocate(itmp)
end subroutine dealloc_ctmp_mod
end module ctmp_mod

module crle_mod
! .. to replace crle.h
 real(kind=8), allocatable  :: trle1(:),trle2(:),trle3(:),trle4(:) 
 integer,      allocatable  :: isize(:)
contains
subroutine init_crle_mod
 use hsfs_parameters
 integer, parameter :: isize_size = 4
#include "priunit.h"
!#include "param_ccsd.h"
allocate(trle1(1:nrle1))
allocate(trle2(1:nrle2))
allocate(trle3(1:nrle3))
allocate(trle4(1:nrle4))
allocate(isize(1:isize_size))
 write(LUPRI,*) 'crle_mod initialized'
end subroutine init_crle_mod

subroutine dealloc_crle_mod
 deallocate(trle1)
 deallocate(trle2)
 deallocate(trle3)
 deallocate(trle4)
 deallocate(isize)
end subroutine dealloc_crle_mod
end module crle_mod

module det_mod
! .. to replace det.h
integer :: nvlspo, ndets,isplit
integer, allocatable  :: lstvso(:),lstdet(:)
contains

subroutine init_det_mod
 use hsfs_parameters
#include "priunit.h"
integer, parameter :: lstvso_size=100, lstdet_size=1000
allocate(lstvso(1:lstvso_size))
allocate(lstdet(1:lstdet_size))
 write(LUPRI,*) 'det_mod initialized'
end subroutine init_det_mod

subroutine dealloc_det_mod
 deallocate(lstvso)
 deallocate(lstdet)
end subroutine dealloc_det_mod
end module det_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module dat_mod
! .. to replace dat.h
integer :: norb, idum2
integer, allocatable  :: ic(:), icv(:), iopt(:)
real*8, allocatable  :: eps(:)
logical :: dat_mod_initialized = .false.
integer, parameter :: dat_size=200, iopt_size=40
contains
subroutine init_dat_mod
 use hsfs_parameters
#include "priunit.h"

if (.not.dat_mod_initialized) then
!
  allocate(ic(1:dat_size))
  allocate(icv(1:dat_size))
  allocate(eps(1:dat_size))
  allocate(iopt(1:iopt_size))
!
  dat_mod_initialized = .true.
  write(LUPRI,*) 'dat_mod initialized'
!
else
  write(LUPRI,*) 'dat_mod was already initialized'
endif
end subroutine init_dat_mod

subroutine print_dat_iopt
! ... prints some iopt(:) data ...
#include "priunit.h"
write(LUPRI,*) 'iopt(1,3,17,19):',iopt(1),iopt(3),iopt(17),iopt(19)
end subroutine print_dat_iopt

subroutine check_dat_mod
implicit none
#include "priunit.h"
integer :: itemp1,itemp2,itemp3,itemp4
real*8 :: temp1,temp2
if (dat_mod_initialized) then
itemp1=ic(1);itemp2=ic(int(dat_size))
temp1=eps(1); temp2=eps(int(dat_size))
itemp3=iopt(1);itemp4=iopt(int(iopt_size))
endif
end subroutine check_dat_mod

subroutine dealloc_dat_mod
 deallocate(ic)
 deallocate(icv)
 deallocate(eps)
 deallocate(iopt)
end subroutine dealloc_dat_mod

end module dat_mod

module idv_mod
! .. to replace idv.h
integer :: nvalo
integer, allocatable  :: iv(:)
integer, parameter :: idv_size=200
contains
subroutine init_idv_mod
! use hsfs_parameters
#include "priunit.h"
allocate(iv(1:idv_size))
write(LUPRI,*) 'idv_mod initialized'
end subroutine init_idv_mod

subroutine check_idv_mod
implicit none
#include "priunit.h"
integer :: itemp1,itemp2
if (allocated(iv)) then
!itemp1=iv(1);itemp2=iv(int(idv_size))
itemp1=iv(1)
itemp2=iv(idv_size)
if (size(iv).ne.idv_size) then
  write(LUPRI,*) 'idv_mod: size(iv),idv_size',size(iv),idv_size
  call quit('idv_mod: size(iv)<>idv_size!')
endif
else
call quit('iv not allocated!')
endif
end subroutine check_idv_mod

subroutine dealloc_idv_mod
 deallocate(iv)
end subroutine dealloc_idv_mod
end module idv_mod

module prm_mod
! .. to replace prm.h
integer, allocatable  :: iperm(:,:)
contains
subroutine init_prm_mod
! use hsfs_parameters
#include "priunit.h"
allocate(iperm(0:6,1:720))
write(LUPRI,*) 'prm_mod initialized - iperm(0:6,1:720)'
end subroutine init_prm_mod

subroutine dealloc_prm_mod
 deallocate(iperm)
end subroutine dealloc_prm_mod

subroutine check_prm_mod
implicit none
#include "priunit.h"
integer :: itemp1,itemp2,itemp3,itemp4
if (allocated(iperm)) then
itemp1=iperm(0,1); itemp2=iperm(0,720)
itemp3=iperm(6,1); itemp4=iperm(6,720)
else
call quit('iperm NOT allocated ! No reason to check !')
endif
end subroutine check_prm_mod

end module prm_mod

module sym_mod
! .. to replace sym.h
integer ::  nsym,nclass
integer, allocatable  :: iclass(:),isym(:,:),ibias(:)
logical :: sym_mod_initialized = .false.
contains
subroutine init_sym_mod
#include "priunit.h"
if (.not.sym_mod_initialized) then
allocate(iclass(1:16))
allocate(isym(1:16,1:8))
allocate(ibias(1:17))
 sym_mod_initialized = .true.
 write(LUPRI,*) 'sym_mod initialized'
else
  write(LUPRI,*) 'sym_mod was not initialized, leaving'
  !call quit('sym_mod was not initialized, leaving!')
endif
end subroutine init_sym_mod

subroutine check_sym_mod
implicit none
#include "priunit.h"
integer ::  itemp1,itemp2,itemp3,itemp4,itemp5,itemp6,itemp7,itemp8
if (sym_mod_initialized) then
  itemp1=iclass(1); itemp2=iclass(16)
  itemp3=isym(1,1); itemp4=isym(16,1);itemp5=isym(1,8);itemp6=isym(16,8)
  itemp7=ibias(1);  itemp8=ibias(17)
else
 call quit('sym_mod: not initialized !')
endif
end subroutine check_sym_mod

subroutine dealloc_sym_mod
 deallocate(iclass)
 deallocate(isym)
 deallocate(ibias)
end subroutine dealloc_sym_mod
end module sym_mod

module ijtemp_mod
! .. to replace ijtemp.h
integer, allocatable  :: ijtmp1(:)
contains
subroutine init_ijtemp_mod
 use hsfs_parameters
 implicit none
#include "priunit.h"
!#include "param_ccsd.h"
allocate(ijtmp1(1:ndp6))
 write(LUPRI,*) 'ijtemp_mod initialized'
end subroutine init_ijtemp_mod

subroutine check_ijtemp_mod
 use hsfs_parameters
implicit none
#include "priunit.h"
integer :: itemp1,itemp2
itemp1=ijtmp1(1);itemp2=ijtmp1(ndp6)
end subroutine check_ijtemp_mod

subroutine dealloc_ijtemp_mod
 deallocate(ijtmp1)
end subroutine dealloc_ijtemp_mod
end module ijtemp_mod

module tv2_mod
! ... to replace tv2.h
real*8, allocatable :: a(:,:)
contains
subroutine init_tv2_mod
implicit none
#include "priunit.h"
integer, parameter :: tv2_size=500
 allocate (a(tv2_size,tv2_size))
 write(LUPRI,*) 'tv2_mod initialized'
end subroutine init_tv2_mod

subroutine dealloc_tv2_mod
 deallocate(a)
end subroutine dealloc_tv2_mod
end module tv2_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module print_level
! module for handling print level, and thus amount of printed data
 integer, private :: iprint=0
 logical(kind=1), private :: print_level_read=.false.
 !logical :: print_level_read = .false.
contains

subroutine extract_print_level
! extract print level from the common block variable
#if defined HSFS_DYN_ALLOC
use dat_mod
#else
#include "dat.h"
#endif
implicit none
#include "priunit.h"
 if (.not.print_level_read) then
  iprint = iopt(1)
  write(LUPRI,'(2X,A,I3)') 'print level, iopt(1), set to:',iprint
 else
  write(LUPRI,'(2X,A,I3)') 'print level was already set before to:',iprint
 endif
 print_level_read = .true.! print variable is extracted
end subroutine extract_print_level

subroutine set_print_level(ipr)
implicit none
#include "priunit.h"
 integer, intent(in) :: ipr
 iprint=ipr
 print_level_read = .true.! print variable is extracted
 write(LUPRI,'(2x,a,i3)') 'set_print_level: new print level set, iprint=',iprint
 write(LUPRI,*) 'set_print_level: print_level_read=',print_level_read
end subroutine set_print_level

function getiprint() result(ipr)
! get private module value
#include "priunit.h"
 integer :: ipr
! write(LUPRI,*) 'getiprint: print_level_read=',print_level_read,' iprint=',iprint
! call FLSHFO(LUPRI)
 if (.not.print_level_read) then
   call quit('getiprint: print variable NOT extracted !')
 else
   ipr = iprint
  return
 endif
 return
end function getiprint

end module print_level

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module high_sect_execs

integer :: num_methods = 0

!character(len=4), allocatable :: method_name(:)
integer, parameter :: max_methods=30
character(len=4) :: method_name(max_methods)

! ... all method-names
character(len=4) :: hs_method_name(19) = (/'csd ','csv ','csu ','csw ','csx ','csy ','csz ',  & 
                                           'chv ','chu ','chw ','chx ','chy ','chz ','chp ',  & 
                                           'chq ','cst ','cht ','m21b','m12b' /)

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_line_arguments
!------------------------------------------------------------------
!
! get arguments for order of methods for the standalone HighSect FSCC code
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
else if (nargs > max_methods) then
  call quit('number of methods exceeded !')
else
 !allocate (method_name(nargs))
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
! ... save number of methods
num_methods=nargs
end subroutine get_line_arguments

function check_method_name(what_method) result (found_method)
!--------------------------------------------------------------
!   checks method name entering in "what_method"
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module ipack_mod
!--------------------------------------------------------------------
! module for integer indexes packing - ipack3, ipack4, ipack5, ipack6
! for retrieving individual indexes from packed indexes i1upack2,i2unpack2
! and with some checking procedures
!--------------------------------------------------------------------
  implicit none
  logical, private :: check_ipack_indexes=.false.
  integer, parameter, private :: min_indx = 1, max_indx = 255
  !private :: set_check_ipack_indexes, unset_check_ipack_indexes
  contains

  subroutine set_check_ipack_indexes
    check_ipack_indexes=.true. ! sets checking of ipack indexes
  end subroutine 

  subroutine unset_check_ipack_indexes
    check_ipack_indexes=.false. ! unsets checking of ipack indexes
  end subroutine 

  integer function ipack2(i,j)
    integer, intent(in) :: i,j
    if (check_ipack_indexes) call check_ipack2_indexes(i,j)
    ipack2=ior(ishft(i,8),j)
  end function 

  integer function i1unpack2(ii)
  ! retreive 1st index from the packed index of two
    integer, intent(in) :: ii
    i1unpack2=iand(ishft(ii,-8),z'ff')
    if (check_ipack_indexes) call check_ipack1_indexes(i1unpack2)
  end function 

  integer function i2unpack2(ii)
  ! retrieve 2nd index from the packed index of two
    integer, intent(in) :: ii
    i2unpack2=iand(ii,z'ff')
    if (check_ipack_indexes) call check_ipack1_indexes(i2unpack2)
  end function 

  integer function ipack3(i,j,k)
    integer, intent(in) :: i,j,k
    if (check_ipack_indexes) call check_ipack3_indexes(i,j,k)
    ipack3=ior(ior(ishft(i,16),ishft(j,8)),k)
  end function 

  integer function ipack4(i,j,k,l)
    integer, intent(in) :: i,j,k,l
    if (check_ipack_indexes) call check_ipack4_indexes(i,j,k,l)
    ipack4=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))
  end function 

  integer function i1unpack4(ii)
  ! unpack 1st index of packed index of four
    integer, intent(in) :: ii
    i1unpack4=iand(ishft(ii,-24),z'ff')
    if (check_ipack_indexes) call check_ipack1_indexes(i1unpack4)
  end function 

  integer function i2unpack4(ii)
  ! unpack 2. index of packed index of four
    integer, intent(in) :: ii
    i2unpack4=iand(ishft(ii,-16),z'ff')
    if (check_ipack_indexes) call check_ipack1_indexes(i2unpack4)
  end function 

  integer function i3unpack4(ii)
  ! unpack 3. index
    integer, intent(in) :: ii
    i3unpack4=iand(ishft(ii,-8),z'ff')
    if (check_ipack_indexes) call check_ipack1_indexes(i3unpack4)
  end function 

  integer function i4unpack4(ii)
  ! retreive 4.index
    integer, intent(in) :: ii
    i4unpack4=iand(ii,z'ff')
    if (check_ipack_indexes) call check_ipack1_indexes(i4unpack4)
  end function 

  integer function ipack5(i,j,k,l,m)
    integer, intent(in) :: i,j,k,l,m
    if (check_ipack_indexes) call check_ipack5_indexes(i,j,k,l,m)
    ipack5=ior(ishft(i,16),ior(ior(ishft(j,12),ishft(k,8)),ior(ishft(l,4),m)))
  end function 

  integer function ipack6(i,j,k,l,m,n)
    integer, intent(in) :: i,j,k,l,m,n
    if (check_ipack_indexes) call check_ipack6_indexes(i,j,k,l,m,n)
    ipack6=ior(ior(ishft(i,20),ishft(j,16)),     &
           ior(ior(ishft(k,12),ishft(l,8)),      & 
               ior(ishft(m,4),n)))
  end function 

  subroutine check_ipack1_indexes(i)
#include "priunit.h"
    integer, intent(in) :: i
    if (.not.(i >= min_indx .and. i <= max_indx )) then
        write(LUPRI,*) 'ipack1: i=',i
        call quit('ipack1 index out of range !')
    endif
  end subroutine 

  subroutine check_ipack2_indexes(i,j)
#include "priunit.h"
    integer, intent(in) :: i,j
    if (.not.(i >= min_indx .and. i <= max_indx  .and.  & 
              j >= min_indx .and. j <= max_indx)) then
        write(LUPRI,*) 'ipack2: i=',i,' j=',j
        call quit('ipack2 index out of range !')
    endif
  end subroutine 

  subroutine check_ipack3_indexes(i,j,k)
#include "priunit.h"
    integer, intent(in) :: i,j,k
    if (.not.(i >= min_indx .and. i <= max_indx  .and.  & 
              j >= min_indx .and. j <= max_indx  .and.  & 
              k >= min_indx .and. k <= max_indx  )) then
        write(LUPRI,*) 'ipack3: i=',i,' j=',j,' k=',k
        call quit('ipack3 index out of range !')
    endif
  end subroutine 

  subroutine check_ipack4_indexes(i,j,k,l)
#include "priunit.h"
    integer, intent(in) :: i,j,k,l
    if (.not.(i >= min_indx .and. i <= max_indx  .and.  & 
              j >= min_indx .and. j <= max_indx  .and.  & 
              k >= min_indx .and. k <= max_indx  .and.  & 
              l >= min_indx .and. l <= max_indx     )) then
        write(LUPRI,*) 'ipack4: i=',i,' j=',j,' k=',k,' l=',l
        call quit('ipack4 index out of range !')
    endif
  end subroutine 

  subroutine check_ipack5_indexes(i,j,k,l,m)
#include "priunit.h"
    integer, intent(in) :: i,j,k,l,m
    if (.not.(i >= min_indx .and. i <= max_indx  .and.  & 
              j >= min_indx .and. j <= max_indx  .and.  & 
              k >= min_indx .and. k <= max_indx  .and.  & 
              l >= min_indx .and. l <= max_indx  .and.  &
              m >= min_indx .and. m <= max_indx   )) then
        write(LUPRI,*) 'ipack5: i=',i,' j=',j,' k=',k,' l=',l,' m=',m
        call quit('ipack5 index out of range !')
    endif
  end subroutine 

  subroutine check_ipack6_indexes(i,j,k,l,m,n)
#include "priunit.h"
    integer, intent(in) :: i,j,k,l,m,n
    if (.not.(i >= min_indx .and. i <= max_indx  .and.  & 
              j >= min_indx .and. j <= max_indx  .and.  & 
              k >= min_indx .and. k <= max_indx  .and.  & 
              l >= min_indx .and. l <= max_indx  .and.  &
              m >= min_indx .and. m <= max_indx  .and.  & 
              n >= min_indx .and. n <= max_indx  )) then
       write(LUPRI,*) 'ipack6: i=',i,' j=',j,' k=',k,' l=',l,' m=',m,' n=',n
       call quit('ipack6 index out of range !')
    endif
  end subroutine 

  subroutine test_ipack_functions
    ! subroutine for demonstration and testing of ipack<2-6> functions
    ! it packs and subsequently unpack set of indexes
#include "priunit.h"
    integer :: i,j,k,l,m,n
    integer :: i0,j0,k0,l0,m0,n0
    integer :: ip2,ip3,ip4,ip5,ip6

    write(LUPRI,'(/,2x,a)') '=== Testing of ipack<2-6> index packing functions ==='
    call set_check_ipack_indexes

    i0=1;j0=255
    ip2=ipack2(i0,j0) 
    i=iand(ishft(ip2,-8),z'ff')
    j=iand(ip2,z'ff')
    write(LUPRI,'(1x,2(a,i3),a,i8,2(a,i3))') 'ipack2(',i0,',',j0,')=',ip2,'=> i=',i,' j=',j

    i0=2;j0=232;k0=6
    ip3=ipack3(i0,j0,k0)
    i=iand(ishft(ip3,-16),z'ff')
    j=iand(ishft(ip3,-8),z'ff')
    k=iand(ip3,z'ff')
    write(LUPRI,'(1x,3(a,i3),a,i8,3(a,i3))')  & 
     'ipack3(',i0,',',j0,',',k0,')=',ip3,'=> i=',i,' j=',j,' k=',k

    i0=2;j0=8;k0=12;l0=137
    ip4=ipack4(i0,j0,k0,l0)
    i=iand(ishft(ip4,-24),z'ff')
    j=iand(ishft(ip4,-16),z'ff')
    k=iand(ishft(ip4, -8),z'ff')
    l=iand(      ip4,     z'ff')
    write(LUPRI,'(1x,4(a,i3),a,i8,4(a,i3))')  & 
    'ipack4(',i0,',',j0,',',k0,',',l0,')=',ip4,'=> i=',i,' j=',j,' k=',k,' l=',l

    i0=201;j0=255;k0=12;l0=48;m0=112
    ip5=ipack5(i0,j0,k0,l0,m0)
    j=iand(ishft(ip5,-24),z'ff')
    k=iand(ishft(ip5,-16),z'ff')
    l=iand(ishft(ip5,-8), z'ff')
    m=iand(      ip5,     z'ff')
    write(LUPRI,'(1x,5(a,i3),a,i8,5(a,i3))')  & 
    'ipack5(',i0,',',j0,',',k0,',',l0,',',m0,')=',ip5,'=> last 4 indexes, i=',0,' j=',j,' k=',k,' l=',l,' m=',m

    i0=251;j0=8;k0=199;l0=20;m0=11;n0=178
    ip6=ipack6(i0,j0,k0,l0,m0,n0)
    write(LUPRI,*) 'ipack6(',i0,',',j0,',',k0,',',l0,',',m0,',',n0,')=',ip6

    write(LUPRI,*)
    call unset_check_ipack_indexes
  end subroutine test_ipack_functions

end module ipack_mod


subroutine check_modules
!##########################################################################
!
!   General checking routine for all dynamically allocated arrays
!
!##########################################################################
#if defined HSFS_DYN_ALLOC
use gen_mod
use cint_mod
use prm_mod
use ctmp_mod
use dat_mod
use idv_mod
use cij1_mod
use sym_mod
use ijtemp_mod
#include "priunit.h"

call check_gen_mod
call check_ctmp_mod
call check_prm_mod
call check_dat_mod
call check_idv_mod
call check_cij1_mod
call check_ijtemp_mod
#else
#include "priunit.h"
 write(LUPRI,*) 'no checking of modules, this is nasty STATIC version.'
#endif
end subroutine check_modules

subroutine space_usage
!> at the end of some method print out info 
!! about the consumed space of gen & igen arrays
!! TODO: compare to the size of allocations
#if defined HSFS_DYN_ALLOC
use cfiles_mod
implicit none
#else
implicit none
#include "cfiles.h"
!common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,iiia(200),iiva(200)
#endif
#include "priunit.h"
write (LUPRI,'(/,i8,a,i8,a)') nstrt(1,nflcur),' words used in igen,', & 
                              nstrt(2,nflcur),' words used in gen.'
end subroutine
