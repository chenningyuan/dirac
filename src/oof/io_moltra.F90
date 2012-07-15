!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module io_moltra
!----------------------------------------------------------------------------------
!
!     contains variables and 'member' routines for reading 
!  moltra integral files MRCONEE & MDCINT
!
!----------------------------------------------------------------------------------
 use Memory
 use Vectors

! ... file numbers ... 
integer :: MDINT=20, MRCONEE=21
CHARACTER*10 DATEX,TIMEX*8
real*8 :: ECORE
integer :: NSP, NFSYM, NZ, NORBT, NSYMRP, NELEC(16), NREP, NKR
integer :: ilenght=0
!--------------------------------------------------------
integer, parameter :: NSP2 = 300 ! static allocation 
INTEGER KR(-NSP2:NSP2)
!--------------------------------------------------------
integer :: irrep_product_dirac(1:8,1:8)
CHARACTER*14 REPNT(16)
!--------------------------------------------------------
      integer, parameter :: MXREP=32, MAXSP=1000
!--------------------------------------------------------
character*4 REPN(16), REPNA(MXREP*2)
integer :: IREPSP(MAXSP*4),MULTB(MXREP*2,MXREP*2,2)
real*8 EPS1(MAXSP) ! orbital energies ...
LOGICAL BREIT,SPFR,CARITH ! Breit interaction, spin-free, complex-arithmetic
CHARACTER*10 NEXTFIL, INTERFACE
logical :: is_complex, mrconee_read = .false.

! either HSFS mode, or 
logical :: HSFS_mode=.false.

integer, public :: nelec_tot=0
integer :: nnzi=0, nintbl=0, nintmaxbl=0
integer :: RCW=0

! store also cbuf, pos_ij, pos_kl
type(Array2Dim), save :: cbuf
type(Array3Dim), save :: pos_ij, pos_kl

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine alloc_mdcint_data
!------------------------------------------------------------------------------------
!
!  Allocate arrays for MDCINT 2-el.integrals & their indexes
!
!  These are:  cbuf (integrals), and pos_ij, pos_kl  - indexes
!
!------------------------------------------------------------------------------------
implicit none
#include "priunit.h"

if (.not.IsMemorySticked() ) then
 write(LUPRI,*) ' memory array not present !' 
 call Quit(' memory array not present !  ')
else
 write(LUPRI,*) ' alloc_mdcint_data: memory sticked !' 
endif

if (ilenght.eq.0.or.RCW.eq.0.or.nnzi.eq.0) then
  call quit('alloc_mdcint_data: some parameters are ready !')
endif

write(LUPRI,*) 'alloc_mdcint_data: calculated array size, ilenght=',ilenght

! ... note: "ccdriv.F":1923 is using RCW !!!
cbuf = Make2DimArray('REAL',RCW,nnzi,'MDCINT two-electron integrals') ! big array !!!!

! third index : 1. - position, 2. - block no
pos_ij = Make3DimArray('INTE',2*NKR,2*NKR,2,'position of (IJ| int.blocks on integral array') 

! ... later take care of +1 index in better way
pos_kl = Make3DimArray('INTE',2*NKR,2*NKR,nintbl+1,'position of |KL) integrals in given block on integral array') ! huge array!

end subroutine alloc_mdcint_data

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine dalloc_mdcint_data

! dealloc cbuf, pos_ij, pos_kl - reverse order !
call DeAllocObject(pos_kl)
call DeAllocObject(pos_ij)
call DeAllocObject(cbuf)

end subroutine dalloc_mdcint_data

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_nonzero_integral_position(i,j,k,l) result(int_pos)
!----------------------------------------------------------------------------
!
!     returns position of NONZERO integral (ij|kl) on big array,
!   or -1 if not found
!
!----------------------------------------------------------------------------
implicit none
#include "priunit.h"
integer, intent(in) :: i,j,k,l
integer :: int_pos

int_pos=-1

! pos_ij, pos_kl must be allocated !
if (.not.pos_ij%is_alloc.or..not.pos_kl%is_alloc) then
  write(LUPRI,*) "get_nonzero_integral_position: arrays are not allocated !"
  call quit('get_nonzero_integral_position: arrays not allocated !')
endif 

if ( pos_ij%p_iA3D(i,j,1).ne.0 .and. pos_ij%p_iA3D(i,j,2).ne.0     &
     .and. pos_kl%p_iA3D(k,l,pos_ij%p_iA3D(i,j,2)).ne.0 )    then

  int_pos=(pos_ij%p_iA3D(i,j,1) + pos_kl%p_iA3D(k,l,pos_ij%p_iA3D(i,j,2)) - 1)
  
 ! write(LUPRI,*) 'pos_ij%p_iA3D(i,j,1) =',pos_ij%p_iA3D(i,j,1) 
 ! write(LUPRI,*) 'pos_ij%p_iA3D(i,j,2),nintbl =',pos_ij%p_iA3D(i,j,2),nintbl
 ! write(LUPRI,*) 'pos_kl%p_iA3D(k,l,pos_ij%p_iA3D(i,j,2)) =',pos_kl%p_iA3D(k,l,pos_ij%p_iA3D(i,j,2))
 ! write(LUPRI,*) 'int_pos=',int_pos

endif

end function get_nonzero_integral_position


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine count_nelec_tot
!-----------------------------------------------------------------------------------
! .. count total number of electrons through DIRAC irreps
!
! the NELEC(irp) must be read before and be nonzero !
!-----------------------------------------------------------------------------------
implicit none
#include "priunit.h"
integer :: irp

if (.not.mrconee_read) call quit('count_nelec_tot: MRCONEE not read !')
nelec_tot = 0
do irp=1,NSYMRP
 write(LUPRI,'(2x,i2,a,i2,a,i2,a)') irp,'/',NSYMRP,' -> ',NELEC(irp),' electrons'
 nelec_tot = nelec_tot + NELEC(irp)
enddo

write(LUPRI,'(/,2x,a,i2,/)')  & 
'count_nelec_tot: total number of electrons is:',nelec_tot

if (nelec_tot.le.0) call quit('count_nelec_tot: wrong number of electrons!')

end subroutine count_nelec_tot

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_mrconee
!----------------------------------------------------------------------------------
!
!  read data from MRCONEE 
!
!----------------------------------------------------------------------------------
implicit none
#include "priunit.h"
integer :: IRP,I,J,IREP,JREP,KREP,ISP,NBSYM
integer, external :: GetNBSYM
logical :: is_mrc, adjust, print_mt = .false.

 inquire(FILE='MRCONEE',EXIST=is_mrc)
 if (.not.is_mrc) call QUIT('mrconee file not present !')

  OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED',status = 'old')
  READ (MRCONEE,END=10,ERR=10) NSP,BREIT,ECORE,NFSYM,NZ,SPFR,NORBT
  READ (MRCONEE,END=10,ERR=10) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP),(NELEC(IRP),IRP=1,NSYMRP)
  READ (MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)
! ... mult.table is incomplete, later modified ...
  READ (MRCONEE) ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP)
! ... IREPSP is to be used as 2-dimensional array, read it like that
  READ (MRCONEE) (IREPSP(ISP),IREPSP(NSP+ISP),EPS1(ISP),ISP=1,NSP)
  CLOSE (MRCONEE,STATUS='KEEP')

  mrconee_read = .true.

  adjust=.true.
  if (adjust) then 
! ... copied from "ccdriv.F"
!C
!C      Right adjust representation names
!C
       DO IRP = 1, NSYMRP
          REPN(IRP) = REPNT(IRP)(1:4)
          CALL BLANKL (4,REPN(IRP))
       ENDDO
       DO IRP = 1, NREP
          CALL BLANKL (4,REPNA(IRP))
       ENDDO
!C
!C     Make inverse of multiplication table
!C
      DO JREP = 1, 2*NREP
         DO IREP = 1, 2*NREP
            KREP = MULTB(IREP,JREP,1)
            MULTB(IREP,KREP,2) = JREP
         ENDDO
      ENDDO
!C
!C     Get rid of offset for boson reps
!C
      DO JREP = 1, NREP
         DO IREP = 1, NREP
            MULTB(IREP,JREP,1) = MULTB(IREP,JREP,1) - NREP
            MULTB(IREP,JREP,2) = MULTB(IREP,JREP,2) - NREP
         ENDDO
      ENDDO

      DO JREP = NREP + 1, 2*NREP
         DO IREP = NREP + 1, 2*NREP
            MULTB(IREP,JREP,1) = MULTB(IREP,JREP,1) - NREP
            MULTB(IREP,JREP,2) = MULTB(IREP,JREP,2) - NREP
         ENDDO
      ENDDO

endif

!-----------------------------------------------------------


! fill irrep_product_dirac !

  irrep_product_dirac=0

 NBSYM = GetNBSYM()
  write(LUPRI,*) 'NBSYM=',NBSYM
! write(LUPRI,*) 'Multiplication table'
  do i=1,NBSYM
  do j=1,NBSYM
   if (print_mt) then
    write(LUPRI,*) i,j,'...',REPNA(MULTB(i,j,1)),REPNA(MULTB(i,j,1)-2*NBSYM)
   endif
   irrep_product_dirac(i,j)= MULTB(i,j,1)-2*NBSYM
  enddo
  enddo

!---------------------------------------------------------------------
!     if HSFS mode, allocate array for HSFS_indexes and fill it 
!---------------------------------------------------------------------
 if (HSFS_mode) call get_HSFS_indexes() 

  return

  10 WRITE (LUPRI,*) "Missing information on MRCONEE, I assume"//  & 
                     " this is a MOLFDIR file ... QUIT !"
  Call Quit('read_mrconee: some error in reading MRCONEE file!')

end subroutine read_mrconee

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_mrconee_data
!----------------------------------------------------------------------------------
!
!       Print read data from the MRCONEE file if these were read 
!
!----------------------------------------------------------------------------------
implicit none
#include "priunit.h"
integer :: multbs(1:8,1:8)
integer :: IRP,ISP,i,j
logical :: print_mt = .false.


if (.not.mrconee_read) then
  Write(LUPRI,*) 'MRCONEE file was not read...'
  call quit('print_mrconee_data:  MRCONEE file was not read!')
endif

Call Header('MRCONEE print-out',-1)

write(LUPRI,'(2x,a,i3)')     'Total # of KU spinors, NSP=',NSP
write(LUPRI,'(2x,a,l1)')      ' is BREIT correction=',BREIT
write(LUPRI,'(2x,a,f12.6)')  ' ECORE=',ECORE
write(LUPRI,'(2x,a,i3)')     ' NFSYM=',NFSYM
write(LUPRI,'(2x,a,i3)')     ' NZ=',NZ

write(LUPRI,'(2x,a,l1)')      ' is spin-free, SPFR=',SPFR

write(LUPRI,'(2x,a,i3)')     ' NORBT=',NORBT
write(LUPRI,'(2x,a,i3)')     ' NSP=',NSP

write(LUPRI,'(/,2x,a,i2)')  'double group irreps, NSYMRP=',NSYMRP

write(LUPRI,'(2x,a)')  '** REPNT(1-NSYMRP) **'
DO IRP = 1, NSYMRP
  write(LUPRI,'(2x,a,i2,a,a)') 'IRP=',IRP,' REPNT(IRP)=',REPNT(IRP)
ENDDO

write(LUPRI,'(/,2x,a)')  '** NELEC(1-NSYMRP) **'
DO IRP = 1, NSYMRP
  write(LUPRI,'(2x,a,i3,a,i3)') 'IRP=',IRP,' NELEC(IRP)=',NELEC(IRP)
ENDDO

write(LUPRI,'(/,2x,a,i3)')  'reducible representations, up to 2*NREP=',2*NREP
write(LUPRI,'(2x,a)')  '** REPNA: **'
DO IRP = 1, 2*NREP
 ! write(LUPRI,'(2x,a,i2,a,a)') 'IRP=',IRP,' REPNA(IRP)=',REPNA(IRP)
  write(LUPRI,'(2x,a,i2,a,i2,a,a)') ' REPNA(',IRP,'/',2*NREP,')=',REPNA(IRP)
ENDDO

write (LUPRI,'(/,2x,a)') '*** EPS1, IREPSP:'
DO ISP = 1, NSP
   write (LUPRI,'(2x,a,i4,a,2i3,a,f18.13,a,a4)') 'ISP=',ISP,' IREPSP(ISP/NSP+ISP):', &
   IREPSP(ISP),IREPSP(NSP+ISP), &
   ' EPS1(ISP)=',EPS1(ISP),'-->',REPNA(IREPSP(NSP+ISP))
ENDDO

if  (print_mt) then

write (LUPRI,'(/,2x,a)') 'Multiplication table; *** ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP) ***'
write (LUPRI,'(2x,a,2i3)') ' NREP, 2*NREP:',NREP,(2*NREP)

do i=1,2*NREP
do j=1,2*NREP
  ! write(LUPRI,'(2X,A,2I3,A,I3)') 'I,J:',I,J,' MULTB(I,J,1):',MULTB(I,J,1)
  write(LUPRI,'(4x,2i3,a,2i4,1x,a4,1x,a4)')   & 
  i,j,' - multb(i,j,1),multb(i,j,2):',MULTB(i,j,1),MULTB(i,j,2),REPNA(MULTB(i,j,1)),REPNA(MULTB(i,j,2))
enddo
enddo

endif

!do i=1,8
!do j=1,8
! write(LUPRI,'(2x,2i3,a,i4,1x,a4)')    & 
!  i,j,' - ',MULTBS(i,j),REPNA(MULTBS(i,j))
!enddo
!enddo

end subroutine print_mrconee_data

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_mdcint
!-------------------------------------------------------------------------
!  read the MDCINT two-electron transformed-MO-integrals file
!
!  we have ready cbuf, pos_ij, pos_kl arrays
!-------------------------------------------------------------------------
implicit none
#include "priunit.h"
type(Vect) :: indk_vect, indl_vect
logical :: is_mdc, check_mem, print_int, print_batch_info 
double precision :: Iijkl
integer :: MOINT,icbuf,iblock,INZ,NONZR,i,j,k,l
integer :: IKR,JKR,int_pos,icount

Call Header('MDCINT control print-out',-1)

! ... check if memory block is present
if (.not.IsMemorySticked() ) then
 write(LUPRI,*) 'read_mdcint: memory array not present !' 
 call Quit('read_mdcint: memory array not present !  ')
else
 write(LUPRI,*) 'read_mdcint: memory sticked !' 
endif

! ... check if MRCONEE was read 
if (.NOT.mrconee_read) then
 write(LUPRI,*) 'read_mdcint: MRCONEE file was not read before !' 
 call Quit('read_mdcint: MRCONEE file was not read before !  ')
!else
! write(LUPRI,*) 'read_mdcint: memory sticked !' 
endif

inquire(file="MDCINT",EXIST=is_mdc); if (.not.is_mdc) call Quit('No MDCINT file!')

open(MDINT,file="MDCINT",form='UNFORMATTED',status='old')
rewind(MDINT)

write(LUPRI,'(//2x,a,i3)') 'The MDCINT file open and rewinded...lu=',MDINT

READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,(KR(I),KR(-I), I=1,NKR )


!---------------------------------------------------------------------
!     if HSFS mode, allocate array for HSFS_indexes and fill it 
!---------------------------------------------------------------------
! if (HSFS_mode) call get_HSFS_indexes()


if (NKR > NSP2) then
 call QUIT('(NKR > NSP2 ! Error quit when reading MDCINT !')
endif

if (NKR /= NORBT .and. 2*NKR /= NORBT) then
 write(LUPRI,*) 'CHECK: NKR=',NKR,' NORBT=',NORBT
 !write(LUPRI,*) 'MRCONEE and MDCINT are incompatible !'
 !call Quit('MRCONEE and MDCINT are incompatible !')
endif

write(LUPRI,*) 'Info from MDCINT file'

! ... determine whether entering 2-electron integrals are real or complex 
is_complex = .NOT.SPFR .AND. (NZ==2 .or. NZ==4)

 !call scan_mdcint(nnzi, nintbl)
 call scan_mdcint

 write(LUPRI,'(2x,a,i8)') 'number of ALL nonzero 2-el MO-integrals on MDCINT=',nnzi
 write(LUPRI,'(2x,a,i7)') 'number of (I J|-integral BLOCKS on MDCINT=',nintbl
 write(LUPRI,'(2x,A,$)') 'Integrals have been calculated on the'
 write(LUPRI,'(A,A12,1X,A,A)') '      ',DATEX,' at ',TIMEX

! -----------------------------------------------------------------
  write(LUPRI,'(2x,a,i3)') 'NKR is ',NKR
  write(LUPRI,'(/,2x,a)') 'KR array, KR(+I),KR(-I), I=1,NKR :'
  do I=1,NKR
   write(LUPRI,'(8X,2I8)') KR(I),KR(-I)
  enddo
!
! ... read the rest ... integrals and indexes  ...
! ... allocate memmory for integrals/indexes   ...

if (is_complex) then
  RCW = 2
else
  RCW = 1
endif
write(LUPRI,'(/2x,a,i1/)') 'the parameter RCW =',RCW

ilenght = (2*NKR)*(2*NKR) ! 

write(LUPRI,*) 'calculated array size, ilenght=NKR*NKR=',ilenght

call alloc_mdcint_data

! call print_mem_param
! ... note: "ccdriv.F":1923 is using RCW !!!
!cbuf = Make2DimArray('REAL',RCW,nnzi,'MDCINT two-electron integrals') ! big array !!!!
!pos_ij = Make3DimArray('INTE',NKR,NKR,2,'position of <ij| int.blocks on integral array') 
!pos_kl = Make3DimArray('INTE',NKR,NKR,nintbl,'position of |kl> integrals in bloks on integral array') ! huge array 
indk_vect = MakeVector('INTE',ilenght,'MDCINT "k" indexes')
indl_vect = MakeVector('INTE',ilenght,'MDCINT "l" indexes')

write(LUPRI,*) 'stuff allocated: indk_vect, indl_vect, cbuf, pos_ij, pos_kl'
    CALL FLSHFO(LUPRI)

! call print_mem_param

!check_mem = .true.
check_mem = .false.
if (check_mem) then
call set_vect_iprint(15)
Call ZeroObject(indk_vect)
Call ZeroObject(indl_vect)
Call ZeroObject(cbuf)
 write(LUPRI,*) 'going to deallocate all objects ...'
 call DeAllocObject(cbuf)
 call DeAllocObject(indl_vect)
 call DeAllocObject(indk_vect)
 write(LUPRI,*) 'going to quit here due to memory check...'
 call Quit('read_mdcint: quit here due to memory check !')
endif


 !print_int = .false.
 print_int = .true.
! print_batch_info = .false.
 print_batch_info = .true.

  rewind(MDINT)
  read(MDINT)

 ! call  set_vect_iprint(30)
  call ZeroObject(cbuf)
  call ZeroObject(pos_ij); call ZeroObject(pos_kl)

  MOINT = 0
  icbuf = 1
  iblock = 0

  icount=0 ! counter for read integrals 

 1 continue

  if (is_complex) then
   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,                             &
   ( indk_vect%p_iA1D(INZ), indl_vect%p_iA1D(INZ),           INZ=1,NONZR,1),  & 
   ( cbuf%p_rA2D(1,INZ+icbuf-1), cbuf%p_rA2D(2,INZ+icbuf-1), INZ=1,NONZR,1) 
  else
   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,                       &
   ( indk_vect%p_iA1D(INZ), indl_vect%p_iA1D(INZ), INZ=1,NONZR,1),      &
   (                   cbuf%p_rA2D(1,INZ+icbuf-1), INZ=1,NONZR,1) 
  endif

  if (NONZR > ilenght) then
    write(LUPRI,*) 'NONZR :',NONZR,' ilenght :',ilenght
    write(LUPRI,*) 'NONZR > ilenght; out of memory allocation !'
    call Quit('NONZR > ilenght; out of memory allocation !')
  endif

  MOINT = MOINT + NONZR

  if (print_batch_info) then
    write(LUPRI,'(/,2x,a,i5)')  '...read (IJ|-integral block:',iblock
    write(LUPRI,'(2x,a,3i4)') 'read IKR,JKR,NONZR:',IKR,JKR,NONZR
    write(LUPRI,'(2x,a,i6)')  ' actual number of read integrals :',MOINT
    write(LUPRI,'(2x,a,i6)')  ' position of cbuf_vect array :',icbuf
    CALL FLSHFO(LUPRI)
  endif

if (print_int) then

 write(LUPRI,'(2x,a,2i4)') 'read ( I J |... indexes:',KR(IKR),KR(JKR)

 do INZ=1, NONZR
   icount = icount + 1
   if (.NOT.is_complex) then
   !  ... unique, KU indexes
      write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12)') icount,'. (',KR(IKR),' ',KR(JKR),'|',  & 
           KR(indk_vect%p_iA1D(INZ)),' ',KR(indl_vect%p_iA1D(INZ)),')=', &
           cbuf%p_rA2D(1,INZ+icbuf-1)
      CALL FLSHFO(LUPRI)
   else
      write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,2f18.12)') icount,'. (',KR(IKR),' ',KR(JKR),'|',  & 
           KR(indk_vect%p_iA1D(INZ)),' ',KR(indl_vect%p_iA1D(INZ)),')=', &
           cbuf%p_rA2D(1,INZ+icbuf-1), cbuf%p_rA2D(2,INZ+icbuf-1)
      CALL FLSHFO(LUPRI)
   endif
 enddo

endif

  pos_ij%p_iA3D(KR(IKR),KR(JKR),1) = icbuf

  icbuf = icbuf + NONZR

  iblock = iblock + 1    ! increment integral block

  pos_ij%p_iA3D(KR(IKR),KR(JKR),2) = iblock

! save KL-indexes of a given block on a big array !
  do INZ=1,NONZR
   ! ... store position on block array
    pos_kl%p_iA3D(KR(indk_vect%p_iA1D(INZ)),KR(indl_vect%p_iA1D(INZ)),iblock) = INZ
  enddo

    IF (IKR.EQ.0) THEN
         IF (JKR.EQ.0) THEN
            NONZR = -1
         ELSE
            READ (MDINT) NEXTFIL
            CLOSE (MDINT,STATUS='KEEP')
            OPEN(MDINT,FILE=NEXTFIL,FORM='UNFORMATTED')
            write(LUPRI,'(2x,a,a)') 'HSCC: proceeding to the next MDCINTxx file ',NEXTFIL
            GOTO 1
         ENDIF
         !GOTO 1
    ENDIF

 if (NONZR == -1) then
  ! call Quit('NONZR == -1!')
 !  RETURN
   GoTo 2
 endif

 write(LUPRI,*) '... jumping to read the next integral batch in the file ... '
 CALL FLSHFO(LUPRI)
 GoTo 1

 2 continue

 write(LUPRI,'(/,2x,a,/)')  &
 ' all integrals from the MDCINT file read ...closing the file'
 CALL FLSHFO(LUPRI)
 CLOSE(MDINT,STATUS='KEEP')

 call DeAllocObject(indl_vect)
 call DeAllocObject(indk_vect)
 write(LUPRI,*) 'indl_vect, indk_vect released ...'; CALL FLSHFO(LUPRI)

! ... print all NONZERO integrals in proper order
 write(LUPRI,'(/,2x,a)') '*** All nonzero MOLTRA integrals in one big array ***'

 do i=1,2*NKR
 do j=1,2*NKR
 do k=1,2*NKR
 do l=1,2*NKR  
    ! must be nonzero integral 
    if ( pos_ij%p_iA3D(i,j,1).ne.0 .and.  & 
         pos_ij%p_iA3D(i,j,2).ne.0 .and.  & 
         pos_kl%p_iA3D(k,l,pos_ij%p_iA3D(i,j,2)).ne.0 ) then 

      int_pos= ( pos_ij%p_iA3D(i,j,1)+pos_kl%p_iA3D(k,l,pos_ij%p_iA3D(i,j,2)) - 1) 

      if (int_pos.gt.nnzi) then ! handle error indexes
       write(LUPRI,*) 'int_pos=',int_pos,' while nnzi=',nnzi
       write(LUPRI,*) 'i,j,pos_ij%p_iA3D(i,j,1):',i,j, pos_ij%p_iA3D(i,j,1)
       write(LUPRI,*) 'pos_ij%p_iA3D(i,j,2):',pos_ij%p_iA3D(i,j,2)
       write(LUPRI,*) 'pos_kl%p_iA3D(i,j,p_iA3D(i,j,2)):',pos_kl%p_iA3D(k,l,pos_ij%p_iA3D(i,j,2))
      endif

      Iijkl=cbuf%p_rA2D(1,int_pos)

      write(LUPRI,'(2x,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,i7)') '(IJ|KL), KU indexes, (',i,' ',j,'|',  & 
                    k,' ',l,')=',Iijkl,' array pos:',int_pos

    endif

 enddo
 enddo
 enddo
 enddo

! call prepare_HSMRCC_Vfiles(cbuf,pos_ij,pos_kl)

 write(LUPRI,*)  &
 ' ... going to deallocate objects in LIFO order ... '
 CALL FLSHFO(LUPRI)

 write(LUPRI,*) 'leavig "read_mdcint"...bye bye'; CALL FLSHFO(LUPRI)

 return

!     ... reading file error outputs ...
 10000 WRITE(LUPRI,'(2x,a)') 'error reading MDCINT' 
      CALL QUIT('ERROR READING HEADER MDCINT')

 10010 WRITE(6,'(/A/A,L10,2A)') & 
      'END OF FILE branch taken when reading 2-e integrals', &
      'CARITH =',CARITH,', INTERFACE =',INTERFACE
       CALL QUIT('END OF FILE READING 2-E INTEGRALS')

 10020 WRITE(6,'(/A/A,L10,2A)')   &
      'ERROR branch taken when reading 2-e integrals',  &
      'CARITH =',CARITH,', INTERFACE =',INTERFACE

end subroutine read_mdcint

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine scan_mdcint
!--------------------------------------------------------------------------------
!
!   Scan (read once) the integral file MDCINT: 
! 
!  -  count total number of nonzero integrals   :      nnzi
!
!  -  count total number of IJ-blocks           :      nintbl
!
!  -  find out the largest integral block with number of integrals : nintmaxbl
!
!--------------------------------------------------------------------------------
implicit none
#include "priunit.h"
integer :: IKR,JKR,NONZR,MOINT,INZ
logical :: is_open, is_file
integer :: myIDUMMY, myIDUMMY2
double precision :: myDUMMY, myDUMMY2

write(LUPRI,'(///,2x,a)') ' Going to count MDCINT integrals & integral blocks:'
write(LUPRI,*) 'is_complex=',is_complex
CALL FLSHFO(LUPRI)

nintbl=0; nnzi=0; nintmaxbl=0

inquire(unit=MDINT,opened=is_open,exist=is_file)

if (.not.is_file) then
  call quit('scan_mdcint: MDCINT file is not present !')
else
  if (.not.is_open) then
    write(LUPRI,'(2x,a)') 'going to open MDCINT file'
    open(MDINT,file="MDCINT",form='UNFORMATTED',status='old')
  else
    write(LUPRI,'(2x,a)') 'going MDCINT file is already open.'
  endif
endif


rewind(MDINT)
read(MDINT) ! skip record: DATEX,TIMEX,NKR,(KR(I),KR(-I), I=1,NKR ) 

 MOINT=0

 1 continue

  if (is_complex) then

   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,           &
    ( myIDUMMY, myIDUMMY2, INZ=1,NONZR,1),                   &
    ( myDUMMY , myDUMMY2, INZ=1,NONZR,1)
  else

   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,           &
   (  myIDUMMY, myIDUMMY2,  INZ=1,NONZR,1),                  &
   (            myDUMMY,    INZ=1,NONZR,1) 

  endif

  MOINT = MOINT + NONZR

  if (IKR.ne.0.and.JKR.ne.0.and.NONZR.ne.0) then
    nintbl = nintbl + 1 ! increment number of integral blocks
    if (NONZR.gt.nintmaxbl) nintmaxbl=NONZR
  endif

  IF (IKR.EQ.0) THEN
         IF (JKR.EQ.0) THEN
            NONZR = -1
         ELSE
            READ (MDINT) NEXTFIL
            CLOSE (MDINT,STATUS='KEEP')
            OPEN(MDINT,FILE=NEXTFIL,FORM='UNFORMATTED')
            write(LUPRI,'(2x,a,a)') '...proceeding to the next MDCINTxx file ',NEXTFIL
            GOTO 1
         ENDIF
         !GOTO 1
  ENDIF

 if (NONZR == -1) then
   GoTo 2
 endif

! write(LUPRI,*) '... jumping to read the next integral batch in the file ... '
! CALL FLSHFO(LUPRI)

 GoTo 1

 2 continue

! finally store total number of integrals
 nnzi=MOINT

 write(LUPRI,'(2x,a,i7)')  &
 'all integrals from the MDCINT file counted...number=',MOINT

 write(LUPRI,'(2x,a,i7)')  &
 'maximum number of integrals from the MDCINT integral block:',nintmaxbl

! write(LUPRI,*) ' while (2*NKR)**4=',16*NKR*NKR*NKR*NKR,' ratio=',DFLOAT(MOINT)/DFLOAT(16*NKR*NKR*NKR*NKR)

 write(LUPRI,*)  &
 'total number of integral blocks:',nintbl
 
return

!     ... reading file error outputs ...
 10000 WRITE(LUPRI,'(2x,a)') 'error reading MDCINT' 
      CALL QUIT('ERROR READING HEADER MDCINT')

 10010 WRITE(6,'(/A/A,L10,2A)') & 
      'END OF FILE branch taken when reading 2-e integrals', &
      'CARITH =',CARITH,', INTERFACE =',INTERFACE
       CALL QUIT('END OF FILE READING 2-E INTEGRALS')

 10020 WRITE(6,'(/A/A,L10,2A)')   &
      'ERROR branch taken when reading 2-e integrals',  &
      'CARITH =',CARITH,', INTERFACE =',INTERFACE

  call quit('scan_mdcint: error brach - wrong reading of integrals ...')

end subroutine scan_mdcint

end module io_moltra

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function GetNBSYM() result(nbsym_dirac)
!--------------------------------------------
!    get number of DIRAC bosonic irreps
!--------------------------------------------
implicit none
#include "priunit.h"
#include "dgroup.h"
 integer:: nbsym_dirac
 if (NBSYM.gt.0) then
   nbsym_dirac = NBSYM
 else
   write(LUPRI,*) 'GetNBSYM: NBSYM from dgroup.h=',NBSYM
   call quit('GetNBSYM: wrong number of DIRAC boson irreps')
 endif
end function GetNBSYM
