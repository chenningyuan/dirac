!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module io_moltra
!----------------------------------------------------------------------------------
!
!        contains variables and 'member' routines for reading 
!  moltra integral files MRCONEE & MDCINT
!
!----------------------------------------------------------------------------------

! ... file numbers ... 
integer :: MDINT=20, MRCONEE=21

CHARACTER*10 DATEX,TIMEX*8
real*8 :: ECORE

integer :: NSP, NFSYM, NZ, NORBT, NSYMRP, NELEC(16), NREP 

!--------------------------------------------------------
integer, parameter :: NSP2 = 300 ! static allocation 
!--------------------------------------------------------

INTEGER KR(-NSP2:NSP2)

CHARACTER*14 REPNT(16)

!--------------------------------------------------------
     integer, parameter :: MXREP=32, MAXSP=1000
!--------------------------------------------------------
character*4 REPN(16), REPNA(MXREP*2)

integer :: IREPSP(MAXSP*4),MULTB(MXREP*2,MXREP*2,2)
real*8 EPS1(MAXSP)

LOGICAL BREIT,SPFR

CHARACTER*10 NEXTFIL

logical :: is_complex, mrconee_read = .false.

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_mrconee
!----------------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"

logical :: is_mrc

 inquire(FILE='MRCONEE',EXIST=is_mrc)
 if (.not.is_mrc) call QUIT('mrconee file not present !')

  OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED',status = 'old')

  READ (MRCONEE,END=10,ERR=10) NSP,BREIT,ECORE,NFSYM,NZ,SPFR,NORBT

  READ (MRCONEE,END=10,ERR=10) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP),  & 
                               (NELEC(IRP),IRP=1,NSYMRP)

  READ (MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)

  READ (MRCONEE) ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP)

! ... IREPSP is to be used as 2-dimensional array, read it like that
  READ (MRCONEE) (IREPSP(ISP),IREPSP(NSP+ISP),EPS1(ISP),ISP=1,NSP)


  CLOSE (MRCONEE,STATUS='KEEP')
  mrconee_read = .true.

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
#include "implicit.h"
#include "priunit.h"

logical :: print_mt = .false.

if (.not.mrconee_read) then
  Write(LUPRI,*) 'MRCONEE file was not read...'
  Return
endif

Call Header('MRCONEE print-out',-1)

write(LUPRI,'(2x,a,i3)')     'Total # of KU spinors, NSP=',NSP
write(LUPRI,'(2x,a,l)')      ' is BREIT correction=',BREIT
write(LUPRI,'(2x,a,f12.6)')  ' ECORE=',ECORE
write(LUPRI,'(2x,a,i3)')     ' NFSYM=',NFSYM
write(LUPRI,'(2x,a,i3)')     ' NZ=',NZ
write(LUPRI,'(2x,a,l)')      ' is spin-free, SPFR=',SPFR
write(LUPRI,'(2x,a,i3)')     ' NORBT=',NORBT

write(LUPRI,'(/,2x,a,i2)')  ' NSYMRP=',NSYMRP

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
   write (LUPRI,'(2x,a,i2,a,2i3,a,f18.13)') 'ISP=',ISP,' IREPSP(ISP/NSP+ISP):', &
         IREPSP(ISP),IREPSP(NSP+ISP), &
         ' EPS1(ISP)=',EPS1(ISP)
ENDDO

write (LUPRI,'(/,2x,a)') '*** ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP) ***'
write (LUPRI,'(2x,a,2i3)') ' NREP, 2*NREP:',NREP,(2*NREP)

if (print_mt) then
  DO I=1,2*NREP
  DO J=1,2*NREP
    write(LUPRI,'(2X,A,2I3,A,I3)') 'I,J:',I,J,' MULTB(I,J,1):',MULTB(I,J,1)
  ENDDO
  ENDDO
endif

end subroutine print_mrconee_data

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_mdcint
!-------------------------------------------------------------------------
!
!    ... read the MDCINT file - testing purposes
!
!-------------------------------------------------------------------------
  use Memory
  use Vectors
#include "implicit.h"
#include "priunit.h"

type(Vect) :: indk_vect, indl_vect
type(Array2Dim) :: cbuf_vect

logical :: is_mdc, check_mem, print_int, print_batch_info 

integer :: RCW

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
!open(MDINT,file="MDCINT",form='UNFORMATTED')
rewind(MDINT)

write(LUPRI,'(//2x,a,i3)') 'The MDCINT file open and rewinded...lu=',MDINT

READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,(KR(I),KR(-I), I=1,NKR )

if (NKR > NSP2) then
 call QUIT('(NKR > NSP2 ! Error quit when reading MDCINT !')
endif

if (NKR /= NORBT .and. 2*NKR /= NORBT) then
 write(LUPRI,*) 'NKR=',NKR,' NORBT=',NORBT
 write(LUPRI,*) 'MRCONEE and MDCINT are incompatible !'
! call Quit('MRCONEE and MDCINT are incompatible !')
endif

 write(LUPRI,*) 'Info from MDCINT file'

 write(LUPRI,'(A,$)') 'Integrals have been calculated on the'
 write(LUPRI,'(A,A12,1X,A,A)') '      ',DATEX,' at ',TIMEX

!do i = 1, NKR
!  write(LUPRI,'(2x,a,i3,a,i3,a,i3)') 'i=',i,' kr(-i)=',KR(-i),' kr(i)=',KR(i)
!enddo

  write(LUPRI,'(2x,a,i3)') 'NKR is ',NKR
  write(LUPRI,'(/,2x,a)') 'KR array, KR(+I),KR(-I), I=1,NKR :'
  do I=1,NKR
     write(LUPRI,'(2I8)') KR(I),KR(-I)
  enddo
!
! ... read the rest ... integrals and indexes  ...
! ... allocate memmory for integrals/indexes   ...

! ... determine whether entering 2-electron integrals are real or complex 
is_complex = .NOT.SPFR .AND. (NZ==2 .or. NZ==4)

if (is_complex) then
  RCW = 2
else
  RCW = 1
endif

write(LUPRI,'(/2x,a,i1/)') 'the parameter RCW =',RCW

ilenght = NKR*NKR ! 

 write(LUPRI,*) 'calculated array size, ilenght=',ilenght

 call print_mem_param

indk_vect = MakeVector('INTE',ilenght,'MDCINT "k" indexes')

 call print_mem_param

indl_vect = MakeVector('INTE',ilenght,'MDCINT "l" indexes')

 call print_mem_param

! ... note: "ccdriv.F":1923 is using RCW !!!
cbuf_vect = Make2DimArray('REAL', RCW, ilenght, 'MDCINT two-electron integrals')

 call print_mem_param

!check_mem = .true.
check_mem = .false.
if (check_mem) then
call set_vect_iprint(15)
Call ZeroObject(indk_vect)
Call ZeroObject(indl_vect)
Call ZeroObject(cbuf_vect)
 write(LUPRI,*) 'going to deallocate all objects ...'
 call DeAllocObject(cbuf_vect)
 call DeAllocObject(indl_vect)
 call DeAllocObject(indk_vect)
 write(LUPRI,*) 'going to quit here due to memory check...'
 call Quit('read_mdcint: quit here due to memory check !')
endif

!Call ZeroObject(indk_vect)
!Call ZeroObject(indl_vect)
!Call ZeroObject(cbuf_vect)

 print_int = .false.
 print_batch_info = .false.

  MOINT = 0

  rewind(MDINT)
  read(MDINT)

 1 continue

  if (is_complex) then

   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,                       &
   ( indk_vect%p_iA1D(INZ), indl_vect%p_iA1D(INZ),     INZ=1,NONZR,1),  & 
   ( cbuf_vect%p_rA2D(1,INZ), cbuf_vect%p_rA2D(2,INZ), INZ=1,NONZR,1) 

  else

   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,                       &
   ( indk_vect%p_iA1D(INZ), indl_vect%p_iA1D(INZ), INZ=1,NONZR,1),      & 
   ( cbuf_vect%p_rA2D(1,INZ), INZ=1,NONZR,1) 

  endif

  if (NONZR > ilenght) then
    write(LUPRI,*) 'NONZR :',NONZR,' ilenght :',ilenght
    write(LUPRI,*) 'NONZR > ilenght; out of memory allocation !'
    call Quit('NONZR > ilenght; out of memory allocation !')
  endif

   MOINT = MOINT + NONZR

  ! call print_vect_attrib(indk_vect)
  ! call print_vect_attrib(indl_vect)

  if (print_batch_info) then
    write(LUPRI,'(2x,a,3i4)') 'read IKR,JKR,NONZR:',IKR,JKR,NONZR
    write(LUPRI,'(2x,a,i6)') ' actual number of read integrals:',MOINT
    CALL FLSHFO(LUPRI)
  endif

if (print_int) then

 do i=1,NONZR
   if (.NOT.is_complex) then
     write(LUPRI,'(2x,a,2i5,f18.12)') 'k,l, real integral:',  & 
           indk_vect%p_iA1D(i),indl_vect%p_iA1D(i),cbuf_vect%p_rA2D(1,i)
     write(LUPRI,'(2x,a,4i4)') 'I,J,K,L:',  & 
           KR(IKR),KR(JKR),KR(indk_vect%p_iA1D(i)),KR(indl_vect%p_iA1D(i))
    ! CALL FLSHFO(LUPRI)
   else
     write(LUPRI,*) 'k,l,cmplx int:',indk_vect%p_iA1D(i),indl_vect%p_iA1D(i),  &
                     cbuf_vect%p_rA2D(1,i), cbuf_vect%p_rA2D(2,i)
    ! CALL FLSHFO(LUPRI)
   endif

 ! ... unique ,KU indexes
   write(LUPRI,'(2x,a,4i5)') 'kramers unrestricted, unique(>0) indexes, I,J,K,L:',  &
     KR(IKR),KR(JKR),KR(indk_vect%p_iA1D(i)),KR(indl_vect%p_iA1D(i))
   CALL FLSHFO(LUPRI)

 enddo

endif
   
    IF (IKR.EQ.0) THEN
         IF (JKR.EQ.0) THEN
            NONZR = -1
         ELSE
            READ (MDINT) NEXTFIL
            CLOSE (MDINT,STATUS='KEEP')
            OPEN(MDINT,FILE=NEXTFIL,FORM='UNFORMATTED')
            write(LUPRI,*) 'HSCC: proceeding to the next file MDCINT.'
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

 write(LUPRI,*)  &
 ' all integrals from the MDCINT file read ...closing the file'
 CALL FLSHFO(LUPRI)

 CLOSE(MDINT,STATUS='KEEP')

 write(LUPRI,*)  &
 ' ... going to deallocate objects in LIFO order ... '
 CALL FLSHFO(LUPRI)

! ... deallocate objects in LIFO order ...
 call DeAllocObject(cbuf_vect)
 write(LUPRI,*) 'cbuf_vect released ...'; CALL FLSHFO(LUPRI)

 call DeAllocObject(indl_vect)
 write(LUPRI,*) 'indl_vect released ...'; CALL FLSHFO(LUPRI)

 call DeAllocObject(indk_vect)
 write(LUPRI,*) 'indk_vect released ...'; CALL FLSHFO(LUPRI)

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

end module io_moltra


module Dirac_HSMRCC
!----------------------------------------------------------------------------------
!
!    handles all necessary data from DIRAC for the standalone HSMRCC code 
!
!----------------------------------------------------------------------------------
integer, parameter :: ludh = 12 ! lu unit of formatted DIRAC_for_HSMRCC file
integer ::  NTMO=0, NCMO=0, NVMO=0, NEMO=0, NHOLE=0, NPTCL=0, NELE=0

integer :: NFRZV, NACT

real(kind=8) :: TotEn, NucEn, TrfTol

integer, public :: iprint_hscc = 0

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine prepare_dirac_HSMRCC_file
!----------------------------------------------------------------------------------
!
! prepare file with DIRAC data (orbitals, symmetry, number of electrons..etc) for the 
! standalone HSMRCC run
!
!----------------------------------------------------------------------------------
use io_moltra
#include "implicit.h"
#include "priunit.h"
!character(len=22) :: dirac_file='Dirac_for_HSMRCC.txt'
!character(len=*) :: dirac_file='Dirac_HSMRCC'
character :: dirac_file='Dirac_HSMRCC'
 
if (iprint_hscc >= 5) then
  call Header('Output from "prepare_dirac_HSMRCC_file"',-1)
endif

nelec_tot = 0
do irp=1,NSYMRP
  nelec_tot = nelec_tot + NELEC(irp)
enddo
write(LUPRI,*) 'nelec_tot=',nelec_tot

!open(ludh,file='Dirac_for_HSMRCC.txt',form='formatted',status='unknown', access='sequential') 
!open(ludh,file=dirac_file,form='formatted',status='unknown', access='sequential') 
open(ludh,file=dirac_file,form='unformatted',status='unknown', access='sequential') 

write(ludh) iprint_hscc

write(ludh) nelec_tot

write(ludh) SPFR

write(ludh) ECORE

write(ludh) NSP

write(ludh) NORBT

write(ludh) NCMO

write(ludh) NACT

close(ludh,status="keep")

write(LUPRI,*) 'saved NELEC_TOT=',nelec_tot
write(LUPRI,*) 'saved NCORE=',NCORE
write(LUPRI,*) 'saved NACT=',NACT
write(LUPRI,*) 'saved NORBT=',NACT
write(LUPRI,*) 'saved NSP=',NSP

if (iprint_hscc >= 5) then
  write(LUPRI,'(/,2x,a,a)') '... (unformatted) Dirac input file for the HSMRCC run ready under name:',dirac_file
endif

!call Quit('prepare_dirac_HSMRCC_file')

end subroutine prepare_dirac_HSMRCC_file

end module Dirac_HSMRCC

