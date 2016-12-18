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
! NOTE: works so far only in C1 symmetry since MOLCAS export MOs in SAO basis
!
!------------------------------------------------------------------------------!

 use carsph, only : sph2car
 use file_handler, only : open_seq, releaseunit

 implicit none

! parameters
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "consts.h"

! common blocks
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "shells.h"
#include "nuclei.h"
#include "blocks.h"

! parameter:
 character(len=*), parameter :: chrdbg = "impmos>"
 logical, parameter :: locdbg = .false.

! input:
 character(len=64), intent(in) :: qcpack, filimp

! local scalars:
 character(len=74) :: text
 integer :: idim(3,2), lucmo, luimp, ioff, ish, nhkti, norbi, khkti
 integer :: idx, ideg, jorb, itmp, mxscr, iofcar, iofsph, nscr, imo, ibf
 real(8) :: toterg, orbsum, chksum

! local allocatables:
 integer, allocatable :: ibeig(:), norb0(:), nbbas0(:)
 real(8), allocatable :: cmosph(:), cmocar(:), eig(:), ctmp(:)

 call qenter("impmos")

 write(lupri,"(/,a)") "    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
 write(lupri,"(  a)") "    !!!        W A R N I N G     !!!"
 write(lupri,"(  a)") "    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
 write(lupri,"(/,a)") "    import (N/S)R MO energies and coefficients from "//qcpack

 ! check if we have C1 symmetry
 if (nbsym.gt.1) then
  call quit(chrdbg//" MO import works so far only in C1 symmetry!")
 end if

 if (qcpack.ne.'MOLCAS') then
  call quit(chrdbg//" MO import works so far only for MOLCAS' InpOrb files!")
 end if
 
 ! determine dimension of scratch array
 mxscr = 0
 do ish = 1,maxshl
  itmp = (2*nhktsh(ish)-1)*norbsh(ish)
  mxscr = max(mxscr,itmp)
 end do

 allocate(cmocar(ncmotq),cmosph(norbt*norbt),eig(norbt),ibeig(norbt),&
          norb0(nbsym),nbbas0(nbsym),ctmp(mxscr))

 ! initialze CMOs with 0
 call dzero(cmocar,ncmotq)

 ! ------------------------------------------------------ !
 ! import MOLCAS information
 ! ------------------------------------------------------ !
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
  !call matprt(cmosph,norbt,norbt)
 end if

 ! ------------------------------------------------------ !
 ! reorder spherical Harm. AOs from MOLCAS to DIRAC order
 ! ------------------------------------------------------ !
 ioff = 0
 do imo = 1,norbt
  do ish = 1,maxshl
   norbi = norbsh(ish)
   nhkti = nhktsh(ish)
   khkti = 2*nhkti-1
   idx = 0
   do ideg = 1,khkti
    do jorb = 1,norbi
     idx = idx+1
     ctmp(ideg+(jorb-1)*khkti) = cmosph(ioff+idx)
    end do
   end do
   nscr = khkti*norbi
   call dcopy(nscr,ctmp,1,cmosph(ioff+1),1)
   ioff = ioff+nscr
  end do
 end do

 if (locdbg) then
  write(lupri,"(a)") chrdbg//" reordered MO coeffs in sph. Harm. basis.:"
  call output(cmosph,1,norbt,1,norbt,norbt,norbt,1,lupri)
  !call matprt(cmosph,norbt,norbt)
 end if

 ! ------------------------------------------------------ !
 ! do the spherical Harm. to Cartesian transformation
 ! ------------------------------------------------------ !
 iofsph = 0
 iofcar = 0
 do imo = 1,norbt
  do ish = 1,kmax
   nhkti = nhkt(ish)
   call sph2car(cmocar(iofcar+1),cmosph(iofsph+1),nhkti)
   iofsph = iofsph+(2*nhkti-1)
   iofcar = iofcar+(nhkti*(nhkti+1)/2)
  end do
 end do

 if (locdbg) then
  write(lupri,"(a)") chrdbg//" reordered MO coeffs. in Cart. basis:"
  call output(cmocar,1,ntbas(0),1,norbt,ntbas(0),norbt,1,lupri)
  !call matprt(cmocar,ntbas(0),norbt)
 end if

 ! ------------------------------------------------------ !
 ! compute check sum that is printed / tested
 ! ------------------------------------------------------ !
 chksum = d0
 iofcar = 0
 do imo = 1,norbt
  orbsum = d0
  do ibf = 1,ntbas(0)
   orbsum = orbsum+cmocar(iofcar+ibf)
  end do
  chksum = chksum+abs(orbsum)
  iofcar = iofcar+ntbas(0)
 end do

 write(lupri,"(/,a,E20.13,/)") "checksum of imported and converted MOs:",chksum

 ! ------------------------------------------------------ !
 ! save imported and converted MOs to DFPCMO file
 ! ------------------------------------------------------ !
 toterg = d0
 ibeig(:) = 0

 lucmo = open_seq("DFPCMO",cform="FORMATTED")
 call wripcmo(lucmo,cmocar,eig,ibeig,toterg)
 call releaseunit(lucmo)

 ! clean up
 deallocate(cmocar,cmosph,eig,ibeig,norb0,nbbas0)

 call qexit("impmos")

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
 logical, parameter :: locdbg = .false.

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
 integer :: kcmo, iuhf, idum, isym, iorb, jorb, luimp, ierr, ibf, jbf
 integer :: keig, nbsym0

 call qenter("rdinporb")

 luimp = open_seq(filimp,cform="FORMATTED")

 ! --------------------------------------------- !
 ! read header information and check consistency
 ! --------------------------------------------- !
 rewind(luimp)

 found = .false.
 do while (.not.found)
  read(luimp,'(a72)',iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error (1) while reading vector source file!')
  found = (index(line,"#INFO").gt.0)
 end do

 read(luimp,'(a72)',iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error (2) while reading vector source file!')

 read(luimp,'(a72)',iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error (3) while reading vector source file!')

 read(line,"(3I8)") iuhf,nbsym0,idum

 ! check if we have restricted wave function
 if (iuhf.ne.0) then
  call quit(chrdbg//" You try to import unrestricted orbitals. Dirac only works in restricted mode!")
 end if

 ! check if number of irreps is consistent
 if (nbsym0.ne.nbsym) then
  call quit(chrdbg//" number of irreps from import file inconsistent with Dirac!")
 end if

 write(fmtorb,'("(",I1,"I8)")')  nbsym0

 read(luimp,'(a72)',iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error (4) while reading vector source file!')
 read(line,fmtorb) (nbbas(isym),isym=1,nbsym0)

 read(luimp,'(a72)',iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error (5) while reading vector source file!')
 read(line,fmtorb) ( norb(isym),isym=1,nbsym0)

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
  if (ierr.ne.0) call quit(chrdbg//' error (6) while reading vector source file!')
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
    if (ierr.ne.0) call quit(chrdbg//' error (7) while reading vector source file!')
    found = (index(line,chrref).gt.0)
   end do

   do ibf=1,nbbas(isym),ndiv
    read(luimp,'(a72)',iostat=ierr) line
    if (ierr.ne.0) call quit(chrdbg//' error (8) while reading vector source file!')
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
  if (ierr.ne.0) call quit(chrdbg//' error (9) while reading vector source file!')
  found = (index(line,"#ONE").gt.0)
 end do

 found = .false.
 do while (.not.found)
  read(luimp,'(a72)',iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error (10) while reading vector source file!')
  found = (index(line,"* ONE ELECTRON ENERGIES").gt.0)
 end do

 keig=0
 do isym=1,nbsym
  do iorb=1,norb(isym),ndiv
   read(luimp,'(a72)',iostat=ierr) line
   if (ierr.ne.0) call quit(chrdbg//' error (11) while reading vector source file!')
   read(line,"(4E18.12)") (eig(jorb+keig),jorb=iorb,min(iorb+3,norb(isym)))
  end do
  keig=keig+norb(isym)
 end do

 ! close file
 call releaseunit(luimp)

 call qexit("rdinporb")

end subroutine rdinporb
!==============================================================================!
