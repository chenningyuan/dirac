!==============================================================================!
subroutine impmos()
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
!   (2) ReSpect MOs for 1C Hamilonians. So far, the MO export is not
!       released at ReSpect side
!
!   (3) Turbomole MOs for 1C Hamiltonians
!
!
! Benjamin Helmich-Paris, fall 2016
!
! NOTE: works so far only in C1 symmetry
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
#include "dcbimp.h"

! parameter:
 character(len=*), parameter :: chrdbg = "impmos>"
 logical, parameter :: locdbg = .false.

! local scalars:
 integer :: lucmo, ish, nhkti, nhktmx, iz
 integer :: itmp, mxscr, iofcar, iofsph, imo, ibf, nz0, minteg, ndeg
 real(8) :: toterg, orbsum, chksum

! local allocatables:
 integer, allocatable :: ibeig(:), norb0(:), nbbas0(:)
 real(8), allocatable :: cmosph(:,:,:), cmocar(:,:,:), eig(:,:), ctmp(:)

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
  nz0  = 4
  ndeg = 2
 else if (nonrel .or. levyle .or. spinfr) then
  nz0  = 1
  ndeg = 1
 end if
 
 ! determine dimension of scratch array
 mxscr = 0
 nhktmx = 0
 do ish = 1,maxshl
  itmp = (2*nhktsh(ish)-1)*norbsh(ish)
  nhktmx = max(nhktmx,nhktsh(ish))
  mxscr = max(mxscr,itmp)
 end do

 allocate(cmocar(ntbas(0),ntbas(0),4),cmosph(norbt,norbt,nz0),eig(norbt,ndeg),ibeig(norbt),&
          norb0(nbsym),nbbas0(nbsym),ctmp(mxscr))

 if (index(qcpack(1:6),'MOLCAS').gt.0) then

  ! ------------------------------------------------------ !
  ! import MOLCAS information
  ! ------------------------------------------------------ !
  call rdinporb(cmosph,eig,norb0,nbbas0,filimp(1),nbsym,norbt,ntbas(0))
 
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
                 filimp(1),kmax,natoms,&! inp
                 ndeg,norbt,nz0) ! dim

  !if (nz0.gt.1) then
  ! do imo = 1,norbt-1
  !  eig(imo+1) = eig(ndeg*imo+1)
  ! end do
  !end if

  if (locdbg) then
   write(lupri,"(a)") chrdbg//" imported MO energies:"
   call output(eig,1,norbt,1,1,ndeg*norbt,1,2,lupri)

   do iz = 1,nz0
    write(lupri,"(a,i1)") chrdbg//" imported MO coeffs in sph. Harm. basis.:",iz
    call output(cmosph(1,1,iz),1,norbt,1,norbt,norbt,norbt,1,lupri)
   end do
  end if
 
  call reorder_respect(cmosph,cmocar,& ! i/o
                       ncntsh,nhktsh,norbsh,& ! inp
                       maxshl,norbt,nz0,natoms)

 else if (index(qcpack(1:9),'TURBOMOLE').gt.0) then

  ! read MO coefficients from mos or spinor.r/i type files
  call rdturbomole(cmosph,eig,filimp,nkrblk,norbt,nz0) ! dim


  if (locdbg) then
   write(lupri,"(a)") chrdbg//" imported MO energies:"
   call output(eig,1,norbt,1,ndeg,norbt,ndeg,1,lupri)

   do iz = 1,nz0
    write(lupri,"(a,i1)") chrdbg//" imported MO coeffs in sph. Harm. basis.:",iz
    call output(cmosph(1,1,iz),1,norbt,1,norbt,norbt,norbt,1,lupri)
   end do
  end if

  ! reorder TODO
  call reorder_turbomole(cmosph,& ! i/o
                         ncntsh,nhktsh,norbsh,nhktmx,& ! inp
                         maxshl,norbt,nz0)

   do iz = 1,nz0
    write(lupri,"(a,i1)") chrdbg//" resorted MO coeffs in sph. Harm. basis.:",iz
    call output(cmosph(1,1,iz),1,norbt,1,norbt,norbt,norbt,1,lupri)
   end do

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
 else if (index(qcpack(1:9),'TURBOMOLE').gt.0) then
  minteg = 3
 else
  minteg = -1
 end if

 ! initialze carsph module
 call carsph_init(nhktmx,minteg)

 ! initialze CMOs with 0
 call dzero(cmocar,ncmotq)

 do iz = 1,nz0
  do imo = 1,norbt
   iofsph = 0
   iofcar = 0
   do ish = 1,kmax
    nhkti = nhkt(ish)
    call sph2car(cmocar(iofcar+1,imo,iz),cmosph(iofsph+1,imo,iz),1,nhkti)
    iofsph = iofsph+(2*nhkti-1)
    iofcar = iofcar+(nhkti*(nhkti+1)/2)
   end do
  end do
 end do

 if (nz0.eq.4) then
  call dscal(ntbas(0)*norbt,-1.d0,cmocar(1,1,3),1)
 end if

 ! destroy carsph module
 call carsph_clear()

 if (locdbg) then
  do iz = 1,nz0
   write(lupri,"(a,i1)") chrdbg//" reordered MO coeffs. in Cart. basis:",iz
   call output(cmocar(1,1,iz),1,ntbas(0),1,norbt,ntbas(0),norbt,1,lupri)
  end do
 end if

 ! ------------------------------------------------------ !
 ! compute check sum that is printed / tested
 ! ------------------------------------------------------ !
 chksum = d0
 do iz = 1,nz0
  do imo = 1,norbt
   orbsum = d0
   do ibf = 1,ntbas(0)
    orbsum = orbsum+cmocar(ibf,imo,iz)
   end do
   chksum = chksum+abs(orbsum)
  end do
 end do

 write(lupri,"(/,a,E20.13,/)") "checksum of imported and converted MOs:",chksum
 call flush(lupri)

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
                    nbsym,norbt,ntbas) ! dim
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
 integer, intent(in) :: nbsym, norbt, ntbas

! input:
 character(len=64), intent(in) :: filimp

! output:
 integer, intent(out) :: norb(nbsym), nbbas(nbsym)
 real(8), intent(out) :: cmo(ntbas*norbt), eig(norbt)

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
                     ndeg,norbt,nz) ! dim
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
 integer, intent(in) :: norbt, nz, ndeg

! input:
 character(len=64), intent(in) :: filimp
 integer, intent(in) :: kmax, ncent

! output:
 real(8), intent(out) :: cmo(norbt*norbt,nz), eig(ndeg,norbt)

! local scalars:
 integer(4) :: i4dum, ndim4, ncent4, nbbas4
 integer :: idx, luimp, iz, ideg
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

 read(luimp,*) (i4dum,idx=1,ncent)
 read(luimp,*) (dummy,idx=1,ncent*3)
 
 !basis set data (uncontacted basis)
 !----------------------------------
 read(luimp,*) nbbas4

 if (int(nbbas4).ne.kmax) then
  call quit("Number of shells read from ReSpect do not match!")
 end if

 read(luimp,*) (i4dum,idx = 1,kmax)
 read(luimp,*) (i4dum,idx = 1,kmax)
 read(luimp,*) (dummy,idx = 1,kmax)
 read(luimp,*) (dummy,idx = 1,kmax)
 
 !MO energies
 !-----------------
 read(luimp,*) ndim4

 if (int(ndim4).ne.ndeg*norbt) then
  call quit("Number of orbitals read from ReSpect do not match!")
 end if

 read(luimp,*) ((eig(ideg,idx),ideg=1,ndeg),idx = 1,norbt)

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
subroutine rdturbomole(cmo,eig,filimp,ndeg,norbt,nz)
!------------------------------------------------------------------------------!
!
! read from TURBOMOLE MO coefficients
!
!  (1) MOs in spherical AO basis 
!
!  (2) MO energies
!
! Benjamin Helmich-Paris, fall 2017
!
!------------------------------------------------------------------------------!

 use file_handler, only : open_seq, releaseunit

 implicit none

#include "priunit.h"

! parameter
 character(len=*), parameter :: chrdbg = "rdturbomole>"
 logical, parameter :: locdbg = .false.
 character(len=7), parameter :: cheader(2) = (/ "$scfmo " , "$spinor" /)
 character(len=4), parameter :: creim(2) = (/ "real", "imag" /)

! dimensions:
 integer, intent(in) :: norbt, nz, ndeg

! input:
 character(len=64), intent(in) :: filimp(ndeg)

! output:
 real(8), intent(out) :: cmo(norbt,norbt*ndeg,nz), eig(norbt,ndeg)

! local scalars:
 character(len=80) :: line
 character(len=12) :: check
 integer :: ideg, luimp, nbscf, imo, idx, kdx, nlen

! local alloctable:
 real(8), allocatable :: suxx(:)

 call qenter("rdturbomole")

 allocate(suxx(norbt*ndeg))

 do ideg = 1,ndeg

  luimp = open_seq(filimp(ideg),cform="FORMATTED")

  rewind(luimp)

  read(luimp,'(a)') line 

  check = adjustl(cheader(ndeg))
  nlen = len_trim(check)
  if (ndeg.gt.1) then
   check(nlen+1:nlen+5) = "_"//creim(ideg)
   nlen = nlen+5
  end if
  
  if (line(1:nlen).ne.check) call quit("Header line wrong!")

  imo = 0
  do 

   read(luimp,'(a)') line 
   if(line(1:4).eq.'$end') exit
   if (line(1:1).eq.'#')   cycle
   kdx = index(line,'eigenvalue=')+11
  
   imo = imo+1

   read(line(kdx:80),*) eig(imo,ideg)
   kdx = index(line,'nsaos=')+6
   read(line(kdx:kdx+6),*) nbscf
   if (norbt*ndeg.ne.nbscf) &
    call quit("Number of basis functions do not match!")

   read(luimp,'(4d20.14)') (suxx(idx),idx=1,norbt*ndeg)

   ! Re(a) / Im(a)
   call dcopy(norbt,suxx         ,1,cmo(1,imo,  ideg),1)

   ! Re(b) / Im(b)
   if (ndeg.gt.1) then
    call dcopy(norbt,suxx(1+norbt),1,cmo(1,imo,2+ideg),1)
   end if

  enddo

  ! close file
  call releaseunit(luimp)

 end do

 deallocate(suxx)

 call qexit("rdturbomole")

end subroutine rdturbomole
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

 call qenter("reorder_molcas")

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

 call qexit("reorder_molcas")

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
 integer :: ishel, jshel, itmp, ihkt, nktmx, lhkt, iz
 integer :: imo, khkti, nhkti, norbi, ioff

! local alloctables:
 integer, allocatable :: imask(:), ncent0(:), nhkt0(:), iofkt(:)

 call qenter("reorder_respect")

 allocate(imask(maxshl),ncent0(maxshl),nhkt0(maxshl),iofkt(maxshl+1))

 ! --------------------------------------------------------- !
 ! 1.) determine order of shell blocks with ReSpect ordering
 ! --------------------------------------------------------- !
 call icopy(maxshl,nhktsh,1,nhkt0,1)

 do ishel = 1,maxshl
  imask(ishel) = ishel
 end do

 call heapsorti(nhkt0,maxshl,.false.,imask)

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

 do ihkt = 1,nktmx
  lhkt = iofkt(ihkt+1)-iofkt(ihkt)
  do jshel = 1,lhkt
   ncent0(jshel) = ncntsh(imask(jshel+iofkt(ihkt)))
  end do
  call heapsorti(ncent0,lhkt,.false.,imask(iofkt(ihkt)+1))
 end do

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

 do iz = 1,nz
  do imo = 1,norbt
   ioff = 0
   do jshel = 1,maxshl
    nhkti = nhkt0(jshel)
    khkti = 2*nhkti-1
    norbi = norbsh(imask(jshel))
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

 call qexit("reorder_respect")

end subroutine reorder_respect
!==============================================================================!
!==============================================================================!
subroutine reorder_turbomole(cmosph,& ! i/o
                             ncntsh,nhktsh,norbsh,nhktmx,& ! inp
                             maxshl,norbt,nz)
!------------------------------------------------------------------------------!
! reorder Turbomole MOs to Dirac MO order
!
!  spherical Harmonics are order in 
!
!   Dirac :  -l, -l+1, ... , 0, 1, ..., +l (linear)
!
!      TM :  0, 1, -1, -2, +2, +3, -3, ... (sinuidal)
!
!------------------------------------------------------------------------------!

 implicit none

#include "priunit.h"
#include "consts.h"
#include "heapsort.h"

! local parameter:
 character(len=*), parameter :: chrdbg = "reorder_respect>"
 logical, parameter :: locdbg = .false.

! dimensions:
 integer, intent(in) :: maxshl, norbt, nz

! input:
 integer, intent(in) :: nhktmx
 integer, intent(in) :: nhktsh(maxshl), ncntsh(maxshl), norbsh(maxshl)

! in/output:
 real(8), intent(inout) :: cmosph(norbt,norbt,nz)

! local scalars:
 integer :: ish, iz, imo, ioff, ioff2, norbi, nhkti
 integer :: khkti, jorb, idx, jdx, mval

! local alloctables:
 real(8), allocatable :: ctmp(:)

 call qenter("reorder_turbomole")

 allocate(ctmp(2*nhktmx-1))

 do iz = 1,nz

  do imo = 1,norbt

   ioff2 = 0

   do ish = 1,maxshl
    norbi = norbsh(ish)
    nhkti = nhktsh(ish)
    khkti = 2*nhkti-1

    ! only for d, f, g, .... functions
    if (nhkti.ge.3) then

     ! --------------------- !
     !
     ! flip phase for ... 
     !
     !  l = 3, m_l = -3     (position 7)
     !
     !  l = 4, m_l = +2, -3 (position 5, 7)
     !
     ! --------------------- !
!     if (nhkti.eq.3) then
!
!      ioff = 0
!      do jorb = 1,norbi
!       cmosph(ioff2+ioff+7,imo,iz) = - cmosph(ioff2+ioff+7,imo,iz)
!       ioff = ioff+khkti
!      end do
!
!     else if (nhkti.eq.4) then
!
!      ioff = 0
!      do jorb = 1,norbi
!       cmosph(ioff2+ioff+5,imo,iz) = - cmosph(ioff2+ioff+5,imo,iz)
!       cmosph(ioff2+ioff+7,imo,iz) = - cmosph(ioff2+ioff+7,imo,iz)
!       ioff = ioff+khkti
!      end do
!
!     end if

     ioff = 0
     do jorb = 1,norbi

      jdx = 0
      ! -l, -l+1, ..., 0
      do mval = -nhkti+1,-1
       jdx = jdx+1
       idx = 2*abs(mval)+mod(abs(mval),2)
       ctmp(jdx) = cmosph(ioff2+ioff+idx,imo,iz)
      end do
      ! 0, 1, ..., l
      do mval = 0,nhkti-1
       jdx = jdx+1
       idx = 2*mval+mod(mval+1,2)
       ctmp(jdx) = cmosph(ioff2+ioff+idx,imo,iz)
      end do

      call dcopy(khkti,ctmp,1,cmosph(ioff2+ioff+1,imo,iz),1)
      ioff = ioff+khkti

     end do

    end if

    ioff2 = ioff2+khkti*norbi

   end do

  end do

 end do

 deallocate(ctmp)

 call qexit("reorder_turbomole")

end subroutine reorder_turbomole
!==============================================================================!
