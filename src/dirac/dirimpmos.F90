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
!   (2) ReSpect MOs for 1C and 2C Hamilonians. So far, the MO export is not
!       released at ReSpect side
!
!
! Benjamin Helmich-Paris, fall 2016
!
! NOTE: works so far only in C1 symmetry since MOLCAS export MOs in SAO basis
!
!------------------------------------------------------------------------------!

 use carsph, only : sph2car, carsph_init, carsph_clear
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
#include "dcbham.h"
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
 integer :: idim(3,2), lucmo, luimp, ioff, ish, nhkti, nhktmx
 integer :: itmp, mxscr, iofcar, iofsph, imo, ibf, nz0, minteg
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

 if (x2cmmf .or. x2c .or. bss .or. dkh1 .or.dkh2) then
  nz0 = 4
 else if (nonrel .or. levyle .or. spinfr) then
  nz0 = 1
 end if
 
 ! determine dimension of scratch array
 mxscr = 0
 nhktmx = 0
 do ish = 1,maxshl
  itmp = (2*nhktsh(ish)-1)*norbsh(ish)
  nhktmx = max(nhktmx,nhktsh(ish))
  mxscr = max(mxscr,itmp)
 end do

 allocate(cmocar(ncmotq),cmosph(norbt*norbt),eig(norbt),ibeig(norbt),&
          norb0(nbsym),nbbas0(nbsym),ctmp(mxscr))

 if (index(qcpack(1:6),'MOLCAS').gt.0) then

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
  end if
 
  ! ------------------------------------------------------ !
  ! reorder spherical Harm. AOs from MOLCAS to DIRAC order
  ! ------------------------------------------------------ !
  call reorder_molcas(cmosph,ctmp,& ! i/o
                      norbsh,nhktsh,& ! inp
                      maxshl,norbt,mxscr)
 
 else if (index(qcpack(1:7),'RESPECT').gt.0) then

  call rdrespect(cmosph,eig, &! out
                 filimp,kmax,natoms,&! inp
                 norbt,nz0) ! dim

  if (locdbg) then
   write(lupri,"(a)") chrdbg//" imported MO energies:"
   call output(eig,1,norbt,1,1,norbt,1,2,lupri)

   write(lupri,"(a)") chrdbg//" imported MO coeffs in sph. Harm. basis.:"
   call output(cmosph,1,norbt,1,norbt,norbt,norbt,1,lupri)
  end if
 
  call reorder_respect(cmosph,cmocar,& ! i/o
                       ncntsh,nhktsh,norbsh,& ! inp
                       maxshl,norbt,nz0,natoms)

 else 
  call quit(chrdbg//" MO import works so far only for MOLCAS and ReSpect!")
 end if

 ! ------------------------------------------------------ !
 ! do the spherical Harm. to Cartesian transformation
 ! ------------------------------------------------------ !

 ! different normalizations of spherical Harmonics
 if      (index(qcpack(1:8),'MOLCAS  ').gt.0) then
  minteg = 2
 else if (index(qcpack(1:8),'RESPECT ').gt.0) then
  minteg = 1
 end if

 ! initialze carsph module
 call carsph_init(nhktmx,minteg)

 ! initialze CMOs with 0
 call dzero(cmocar,ncmotq)

 iofsph = 0
 iofcar = 0
 do imo = 1,norbt
  do ish = 1,kmax
   nhkti = nhkt(ish)
   call sph2car(cmocar(iofcar+1),cmosph(iofsph+1),1,nhkti)
   iofsph = iofsph+(2*nhkti-1)
   iofcar = iofcar+(nhkti*(nhkti+1)/2)
  end do
 end do

 ! destroy carsph module
 call carsph_clear()

 if (locdbg) then
  write(lupri,"(a)") chrdbg//" reordered MO coeffs. in Cart. basis:"
  call output(cmocar,1,ntbas(0),1,norbt,ntbas(0),norbt,1,lupri)
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
! Benjamin Helmich-Paris, fall 2016
!
! NOTE: Only InpOrb formats 1.1 and 2.0 are recognized so far.
!
!------------------------------------------------------------------------------!

 use file_handler, only : open_seq, releaseunit

 implicit none

#include "priunit.h"

! parameter
 character(len=*), parameter :: chrdbg = "rdinporb>"
 logical, parameter :: locdbg = .false.

! dimensions:
 integer, intent(in) :: nbsym, ntorb, ntbas

! input:
 character(len=64), intent(in) :: filimp

! output:
 integer, intent(out) :: norb(nbsym), nbbas(nbsym)
 real(8), intent(out) :: cmo(ntbas*ntorb), eig(ntorb)

! local scalars
 character(len=3)   :: chrver
 character(len=5)   :: fmtsym
 character(len=6)   :: fmtlin
 character(len=9)   :: fmtone
 character(len=9)   :: fmtorb
 character(len=19)  :: chrref
 character(len=120) :: line
 logical :: found
 integer :: kcmo, iuhf, idum, isym, iorb, jorb, luimp, ierr, ibf, jbf
 integer :: keig, nbsym0, ncorb, ncone

 call qenter("rdinporb")

 luimp = open_seq(filimp,cform="FORMATTED")

 rewind(luimp)

 ! --------------------------------------------- !
 ! figure out which INPORB format is used
 ! --------------------------------------------- !
 found = .false.
 do while (.not.found)
  read(luimp,'(a11)',iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error (1) while reading vector source file!')
  found = (index(line,"#INPORB").gt.0)
 end do

 chrver = line(9:11)

 select case (chrver)
  case ('1.1')
   write(lupri,"(a/)") "read IMPORB type of file "//trim(filimp)//" with format version "//chrver//" ..."
   fmtlin = "(a72) "
   fmtorb = "(4E18.12)"
   fmtone = "(4E18.12)"
   ncorb = 4
   ncone = 4
  case ('2.0')
   write(lupri,"(a/)") "read IMPORB type of file "//trim(filimp)//" with format version "//chrver//" ..."
   fmtlin = "(a110)"
   fmtorb = "(5E22.14)"
   fmtone = "(10E12.4)"
   ncorb = 5
   ncone = 10
  case default
   call quit('Unsupported MOLCAS INPORB file format (so far only 1.1 and 2.0)!')
 end select

 ! --------------------------------------------- !
 ! read header information and check consistency
 ! --------------------------------------------- !
 found = .false.
 do while (.not.found)
  read(luimp,fmtlin,iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error (1) while reading vector source file!')
  found = (index(line,"#INFO").gt.0)
 end do

 read(luimp,fmtlin,iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error (2) while reading vector source file!')

 read(luimp,fmtlin,iostat=ierr) line
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

 write(fmtsym,'("(",I1,"I8)")')  nbsym0

 read(luimp,fmtlin,iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error (4) while reading vector source file!')
 read(line,fmtsym) (nbbas(isym),isym=1,nbsym0)

 read(luimp,fmtlin,iostat=ierr) line
 if (ierr.ne.0) call quit(chrdbg//' error (5) while reading vector source file!')
 read(line,fmtsym) ( norb(isym),isym=1,nbsym0)

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
  read(luimp,fmtlin,iostat=ierr) line
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
    read(luimp,fmtlin,iostat=ierr) line
    if (ierr.ne.0) call quit(chrdbg//' error (7) while reading vector source file!')
    found = (index(line,chrref).gt.0)
   end do

   do ibf=1,nbbas(isym),ncorb
    read(luimp,fmtlin,iostat=ierr) line
    if (ierr.ne.0) call quit(chrdbg//' error (8) while reading vector source file!')
    read(line,fmtorb) (cmo(jbf+kcmo),jbf=ibf,min(ibf+ncorb-1,nbbas(isym)))
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
  read(luimp,fmtlin,iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error (9) while reading vector source file!')
  found = (index(line,"#ONE").gt.0)
 end do

 found = .false.
 do while (.not.found)
  read(luimp,fmtlin,iostat=ierr) line
  if (ierr.ne.0) call quit(chrdbg//' error (10) while reading vector source file!')
  found = (index(line,"* ONE ELECTRON ENERGIES").gt.0)
 end do

 keig=0
 do isym=1,nbsym
  do iorb=1,norb(isym),ncone
   read(luimp,fmtlin,iostat=ierr) line
   if (ierr.ne.0) call quit(chrdbg//' error (11) while reading vector source file!')
   read(line,fmtone) (eig(jorb+keig),jorb=iorb,min(iorb+ncone-1,norb(isym)))
  end do
  keig=keig+norb(isym)
 end do

 ! close file
 call releaseunit(luimp)

 call qexit("rdinporb")

end subroutine rdinporb
!==============================================================================!
!==============================================================================!
subroutine rdrespect(cmo,eig, &! out
                     filimp,kmax,ncent,&! inp
                     ntorb,nz) ! dim
!------------------------------------------------------------------------------!
!
! read from ReSpect temporary file provided by Michal Repisky
!
!  (1) MOs in spherical AO basis 
!
!  (2) MO energies
!
! Benjamin Helmich-Paris, spring 2017
!
!------------------------------------------------------------------------------!

 use file_handler, only : open_seq, releaseunit

 implicit none

#include "priunit.h"

! parameter
 character(len=*), parameter :: chrdbg = "rdinporb>"
 logical, parameter :: locdbg = .false.

! dimensions:
 integer, intent(in) :: ntorb, nz

! input:
 character(len=64), intent(in) :: filimp
 integer, intent(in) :: kmax, ncent

! output:
 real(8), intent(out) :: cmo(ntorb*ntorb,nz), eig(ntorb)

! local scalars:
 integer(4) :: i4dum, ndim4, ncent4, nbbas4
 integer :: idx, luimp, iz
 real(8) :: dummy

 call qenter("rdrespect")

 luimp = open_seq(filimp,cform="FORMATTED")

 rewind(luimp)

 !atomic data 
 !-----------
 read(luimp,*)  ncent4

 if (int(ncent4).ne.ncent) then
  call quit("Number of centers read from ReSpect do not match!")
 end if

 read(luimp,*) (i4dum,idx=1,ncent4)
 read(luimp,*) (dummy,idx=1,ncent4*3)
 
 !basis set data (uncontacted basis)
 !----------------------------------
 read(luimp,*) nbbas4

 if (int(nbbas4).ne.kmax) then
  call quit("Number of shells read from ReSpect do not match!")
 end if

 read(luimp,*) (i4dum,idx = 1,nbbas4)
 read(luimp,*) (i4dum,idx = 1,nbbas4)
 read(luimp,*) (dummy,idx = 1,nbbas4)
 read(luimp,*) (dummy,idx = 1,nbbas4)
 
 !MO energies
 !-----------------
 read(luimp,*) ndim4
 read(luimp,*) (eig(idx),idx = 1,ndim4)

 !MO coefficients
 !-----------------
 read(luimp,*) ndim4
 do iz = 1,nz
  read(luimp,*) (cmo(idx,iz),idx = 1,ndim4)
 end do

 ! close file
 call releaseunit(luimp)

 call qexit("rdrespect")

end subroutine rdrespect
!==============================================================================!
!==============================================================================!
subroutine reorder_molcas(cmosph,ctmp,& ! i/o
                          norbsh,nhktsh,& ! inp
                          maxshl,norbt,mxscr)
!------------------------------------------------------------------------------!
! reorder MOLCAS MOs to Dirac MO order
!------------------------------------------------------------------------------!

 implicit none

#include "priunit.h"

! local parameter:
 character(len=*), parameter :: chrdbg = "reorder_molcas>"
 logical, parameter :: locdbg = .false.

! dimensions:
 integer, intent(in) :: maxshl, norbt, mxscr

! input:
 integer, intent(in) :: norbsh(maxshl), nhktsh(maxshl)

! in/output:
 real(8), intent(inout) :: ctmp(mxscr), cmosph(norbt*norbt)

! local scalars:
 integer :: ioff, imo, ish, norbi, nhkti, khkti, idx, ideg
 integer :: jorb, nscr

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
 end if

end subroutine reorder_molcas
!==============================================================================!
!==============================================================================!
subroutine reorder_respect(cmosph,scratch,& ! i/o
                           ncntsh,nhktsh,norbsh,& ! inp
                           maxshl,norbt,nz,natoms)
!------------------------------------------------------------------------------!
! reorder ReSpect MOs to Dirac MO order
!------------------------------------------------------------------------------!

 implicit none

#include "priunit.h"
#include "heapsort.h"

! local parameter:
 character(len=*), parameter :: chrdbg = "reorder_respect>"
 logical, parameter :: locdbg = .false.

! dimensions:
 integer, intent(in) :: maxshl, norbt, nz, natoms

! input:
 integer, intent(in) :: nhktsh(maxshl), ncntsh(maxshl), norbsh(maxshl)

! in/output:
 real(8), intent(inout) :: scratch(norbt,norbt,nz), cmosph(norbt,norbt,nz)

! local scalars:
 integer :: ishel, jshel, iat, isum, itmp, ihkt, nktmx, lhkt, iz
 integer :: imo, khkti, nhkti, norbi, iorb, ioff

! local alloctables:
 integer, allocatable :: imask(:), ncent0(:), nhkt0(:), iofkt(:)

 allocate(imask(maxshl),ncent0(maxshl),nhkt0(maxshl),iofkt(maxshl+1))

 ! --------------------------------------------------------- !
 ! 1.) determine order of shell blocks with ReSpect ordering
 ! --------------------------------------------------------- !
 call icopy(maxshl,nhktsh,1,nhkt0,1)

 do ishel = 1,maxshl
  imask(ishel) = ishel
 end do

 call heapsorti(nhkt0,maxshl,.false.,imask)

! write(6,*) " nhkt:", nhktsh(1:maxshl)
! write(6,*) "nhkt0:", nhkt0
! write(6,*) "imask:", imask

 itmp = 0
 ihkt = 0
 do ishel = 1,maxshl
  if (nhkt0(ishel).ne.itmp) then
   itmp = nhkt0(ishel)
   ihkt = ihkt+1
   iofkt(ihkt) = ishel-1
  end if
 end do
 nktmx = ihkt
 iofkt(nktmx+1) = maxshl

! write(6,*) "iofkt:", iofkt(1:nktmx+1)

 do ihkt = 1,nktmx
  lhkt = iofkt(ihkt+1)-iofkt(ihkt)
  do jshel = 1,lhkt
   ncent0(jshel) = ncntsh(imask(jshel+iofkt(ihkt)))
  end do
  call heapsorti(ncent0,lhkt,.false.,imask(iofkt(ihkt)+1))
 end do

! write(6,*) "imask:", imask

 ! --------------------------------------------------------- !
 ! 2.) re-order basis functions in Dirac order for each MO
 ! --------------------------------------------------------- !

 ioff = 0
 do ishel = 1,maxshl
  nhkti = nhktsh(ishel)
  khkti = 2*nhkti-1
  norbi = norbsh(ishel)
  iofkt(ishel) = ioff
  ioff = ioff+khkti*norbi
 end do
 iofkt(maxshl+1) = ioff

! write(6,*) "iofkt:", iofkt

 do iz = 1,nz
  do imo = 1,norbt
   ioff = 0
   do jshel = 1,maxshl
    nhkti = nhkt0(jshel)
    khkti = 2*nhkti-1
    norbi = norbsh(imask(jshel))
!    write(6,*) "debug:",jshel,imask(jshel),nhkti,norbi
    call dcopy(khkti*norbi,cmosph(1+ioff,imo,iz),1,&
                           scratch(1+iofkt(imask(jshel)),imo,iz),1)
    ioff = ioff+khkti*norbi
   end do
  end do
 end do

 if (locdbg) then
  write(lupri,"(a)") chrdbg//" reordered MO coeffs in sph. Harm. basis.:"
  call output(scratch,1,norbt,1,norbt,norbt,norbt,1,lupri)
 end if

 call dcopy(nz*norbt*norbt,scratch,1,cmosph,1)

 deallocate(imask,ncent0,nhkt0,iofkt)

end subroutine reorder_respect
!==============================================================================!
