!==============================================================================!
subroutine impmos(qcpack,filimp)
!------------------------------------------------------------------------------!
!
! import MO coefficients and energies from ASCII files of another QC package
!
!   (1) MOLCAS InpOrb type of files (ScfOrb, RasOrb, ...) for non-relativistic
!       or scalar-relativic calculations
!
!       MOs need to be converted from the spherical Harmonic into the
!       Cartesian AO basis
!
!
! Benjamin Helmich-Paris, fall 2016
!
!------------------------------------------------------------------------------!

 use car2sph
 use file_handler, only : getunit, releaseunit

 implicit none

! parameters
#include "priunit.h"

! common blocks
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

! parameter:
 character(len=*), parameter :: chrdbg = "impmos>"
 logical, parameter :: locdbg = .true.

! input:
 character(len=64), intent(in) :: qcpack, filimp

! local scalars:
 character(len=74) :: text
 integer :: idim(3,2), lucmo, luimp
 real(8) :: toterg

! local allocatables:
 integer, allocatable :: ibeig(:), norb0(:), nbbas0(:)
 real(8), allocatable :: cmosph(:), cmocar(:), eig(:)

! functions:

 call qenter("impmos")

 allocate(cmocar(ncmotq),cmosph(norbt*norbt),eig(norbt),ibeig(norbt),&
          norb0(nbsym),nbbas0(nbsym))

 ! import MOLCAS information
 call rdinporb(cmosph,eig,norb0,nbbas0,filimp,nbsym,norbt,ntbas(0))

 if (locdbg) then
  write(lupri,"(a)") chrdbg//" imported MOs per irrep:"
  write(lupri,"(5I5)") norb0

  write(lupri,"(a)") chrdbg//" imported AOs per irrep:"
  write(lupri,"(5I5)") nbbas0

  write(lupri,"(a)") chrdbg//" imported MO energies:"
  call output(eig,1,norbt,1,1,norbt,1,2,lupri)

  write(lupri,"(a)") chrdbg//" imported MO coeffs.:"
  call output(cmosph,1,norbt,1,norbt,norbt,norbt,2,lupri)
 end if

 ! save imported and converted MOs to DFPCMO file
 lucmo = getunit("DFPCMO")
 call wripcmo(lucmo,cmocar,eig,ibeig,toterg)
 call releaseunit(lucmo)

 ! clean up
 deallocate(cmocar,cmosph,eig,ibeig,norb0,nbbas0)

 call qexit("impmos")

 stop 'panic'

end subroutine impmos
!==============================================================================!
!==============================================================================!
subroutine rdinporb(cmo,eig,norb,nbbas, &! out
                    filimp,&! inp
                    nbsym,ntorb,ntbas) ! dim
!------------------------------------------------------------------------------!
!
! read from MOLCAS' InpOrb type of files
!
!  (1) MOs in spherical AO basis 
!
!  (2) MO energies
!
! inspired by J. Olsen's / T. Fleig's RDVEC routine
!
!------------------------------------------------------------------------------!

 use file_handler, only : open_seq, releaseunit

 implicit none

#include "priunit.h"

! parameter
 character(len=*), parameter :: chrdbg = "rdinporb>"
 logical, parameter :: locdbg = .true.

 integer, parameter :: ndiv  = 4

! dimensions:
 integer, intent(in) :: nbsym, ntorb, ntbas

! input:
 character(len=64), intent(in) :: filimp

! output:
 integer, intent(out) :: norb(nbsym), nbbas(nbsym)
 real(8), intent(out) :: cmo(ntbas*ntorb), eig(ntorb)

! local scalars
 character(len=5)  :: fmtorb
 character(len=19) :: chrref
 character(len=72) :: line
 logical :: found
 integer :: kcmo, idum1, idum2, isym, iorb, jorb, luimp, ierr, ibf, jbf
 integer :: keig, nbsym0

 call qenter("rdinporb")

 norb(:) = -999

 luimp = open_seq(filimp,cform="FORMATTED")

 ! --------------------------------------------- !
 ! read header information and check consistency
 ! --------------------------------------------- !
 rewind(luimp)

 found = .false.
 do while (.not.found)
  read(luimp,'(a72)',iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')
  found = (index(line,"#INFO").gt.0)
 end do

 read(luimp,'(a72)',iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')

 read(luimp,'(a72)',iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')

 read(line,"(3I8)") idum1,nbsym0,idum2
 
 ! check if number of irreps is consistent
 if (nbsym0.ne.nbsym) then
  call quit(chrdbg//" number of irreps from import file inconsistent with Dirac!")
 end if

 write(fmtorb,'("(",I1,"I8)")')  nbsym0

! TODO check that
 read(luimp,'(a72)',iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')
 read(line,fmtorb) (norb(isym),isym=1,nbsym0)

!! check if number of MOs per irrep is consistent
!do isym = 1,nbsym
! if (norb0(isym).ne.norb(isym)) then
!  call quit(chrdbg//" number of MOs per irrep from import file inconsistent with Dirac!")
! end if
!end do

 read(luimp,'(a72)',iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')
 read(line,fmtorb) (nbbas(isym),isym=1,nbsym0)

!! check if number of basis functions per irrep are consistent
!do isym = 1,nbsym
! if (nbbas0(isym).ne.nbbas(isym)) then
!  call quit(chrdbg//" number of MOs per irrep from import file inconsistent with Dirac!")
! end if
!end do
! TODO check that

 if (locdbg) then
  write(lupri,*) chrdbg//" number of irreps:",nbsym0
  write(lupri,*) chrdbg//" number of MOs:",( norb(isym),isym=1,nbsym0)
  write(lupri,*) chrdbg//" number of AOs:",(nbbas(isym),isym=1,nbsym0)
 end if

 ! -------------------- !
 ! read MO coefficients
 ! -------------------- !
 rewind(luimp)

 found = .false.
 do while (.not.found)
  read(luimp,'(a72)',iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')
  found = (index(line,"#ORB").gt.0)
 end do

 kcmo  = 0
 do isym=1,nbsym
  do iorb=1,norb(isym)

   ! reference string to look for
   chrref(1:9) = "* ORBITAL"
   write(chrref(10:14),"(I5)") isym
   write(chrref(15:19),"(I5)") iorb

   ! find proper orbital string
   found = .false.
   do while (.not.found)
    read(luimp,'(a72)',iostat=ierr) line
    if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')
    found = (index(line,chrref).gt.0)
   end do

   do ibf=1,nbbas(isym),ndiv
    read(luimp,'(a72)',iostat=ierr) line
    if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')
    read(line,"(4E18.12)") (cmo(jbf+kcmo),jbf=ibf,min(ibf+3,nbbas(isym)))
   end do
   kcmo=kcmo+nbbas(isym)
  end do
 end do

 ! ---------------- !
 ! read MO energies
 ! ---------------- !
 rewind(luimp)

 found = .false.
 do while (.not.found)
  read(luimp,'(a72)',iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')
  found = (index(line,"#ONE").gt.0)
 end do

 found = .false.
 do while (.not.found)
  read(luimp,'(a72)',iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')
  found = (index(line,"* ONE ELECTRON ENERGIES").gt.0)
 end do

 keig=0
 do isym=1,nbsym
  do iorb=1,norb(isym),ndiv
   read(luimp,'(a72)',iostat=ierr) line
   if (ierr.ne.0) call quit(chrdbg//' error while reading vector source file!')
   read(line,"(4E18.12)") (eig(jorb+keig),jorb=iorb,min(iorb+3,norb(isym)))
  end do
  keig=keig+norb(isym)
 end do

 ! close file
 call releaseunit(luimp)

 call qexit("rdinporb")

end subroutine rdinporb
!==============================================================================!
