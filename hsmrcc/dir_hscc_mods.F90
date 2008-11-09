
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

integer :: irrep_product_dirac(1:8,1:8)

INTEGER KR(-NSP2:NSP2)

CHARACTER*14 REPNT(16)

!--------------------------------------------------------
      integer, parameter :: MXREP=32, MAXSP=1000
!--------------------------------------------------------
character*4 REPN(16), REPNA(MXREP*2)

integer :: IREPSP(MAXSP*4),MULTB(MXREP*2,MXREP*2,2)

real*8 EPS1(MAXSP) ! orbital energies ...

LOGICAL BREIT,SPFR

CHARACTER*10 NEXTFIL

logical :: is_complex, mrconee_read = .false.

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function GetNBSYM() result(nbsym_dirac)
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
 nbsym_dirac = NBSYM
 write(LUPRI,*) 'get NBSYM=',NBSYM
end function GetNBSYM

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_mrconee
!----------------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"

logical :: is_mrc
!integer, external :: GetNBSYM

 inquire(FILE='MRCONEE',EXIST=is_mrc)
 if (.not.is_mrc) call QUIT('mrconee file not present !')

  OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED',status = 'old')

  READ (MRCONEE,END=10,ERR=10) NSP,BREIT,ECORE,NFSYM,NZ,SPFR,NORBT

  READ (MRCONEE,END=10,ERR=10) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP),  & 
                               (NELEC(IRP),IRP=1,NSYMRP)

  READ (MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)

! ... mult.table is incomplete, later modified ...
  READ (MRCONEE) ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP)

! ... IREPSP is to be used as 2-dimensional array, read it like that
  READ (MRCONEE) (IREPSP(ISP),IREPSP(NSP+ISP),EPS1(ISP),ISP=1,NSP)

  CLOSE (MRCONEE,STATUS='KEEP')

  mrconee_read = .true.


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

!-----------------------------------------------------------

! fill irrep_product_dirac !

  irrep_product_dirac=0

  write(LUPRI,*) 'GetNBSYM=',NBSYM
  do i=1,NBSYM
  do j=1,NBSYM
   ! write(LUPRI,*) i,j,'...',REPNA(MULTB(i,j,1))
    write(LUPRI,*) i,j,'...',REPNA(MULTB(i,j,1)),REPNA(MULTB(i,j,1)-2*NBSYM)
    irrep_product_dirac(i,j)= MULTB(i,j,1)-2*NBSYM
  enddo
  enddo

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

integer :: multbs(1:8,1:8)

!logical :: print_mt = .false.

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
   write (LUPRI,'(2x,a,i2,a,2i3,a,f18.13,a,a4)') 'ISP=',ISP,' IREPSP(ISP/NSP+ISP):', &
   IREPSP(ISP),IREPSP(NSP+ISP), &
   ' EPS1(ISP)=',EPS1(ISP),'-->',REPNA(IREPSP(NSP+ISP))
ENDDO

write (LUPRI,'(/,2x,a)') 'Multiplication table; *** ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP) ***'
write (LUPRI,'(2x,a,2i3)') ' NREP, 2*NREP:',NREP,(2*NREP)

do i=1,2*NREP
do j=1,2*NREP
  ! write(LUPRI,'(2X,A,2I3,A,I3)') 'I,J:',I,J,' MULTB(I,J,1):',MULTB(I,J,1)
  write(LUPRI,'(4x,2i3,a,2i4,1x,a4,1x,a4)')   & 
  i,j,' - multb(i,j,1),multb(i,j,2):',MULTB(i,j,1),MULTB(i,j,2),REPNA(MULTB(i,j,1)),REPNA(MULTB(i,j,2))
enddo
enddo

write(LUPRI,*)

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

! call print_mem_param

indk_vect = MakeVector('INTE',ilenght,'MDCINT "k" indexes')

! call print_mem_param

indl_vect = MakeVector('INTE',ilenght,'MDCINT "l" indexes')

! call print_mem_param

! ... note: "ccdriv.F":1923 is using RCW !!!
cbuf_vect = Make2DimArray('REAL', RCW, ilenght, 'MDCINT two-electron integrals')

! call print_mem_param

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

 !print_int = .false.
 print_int = .true.
! print_batch_info = .false.
 print_batch_info = .true.

  MOINT = 0

  rewind(MDINT)
  read(MDINT)

 1 continue

  if (is_complex) then

   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,                       &
   ( indk_vect%p_iA1D(INZ), indl_vect%p_iA1D(INZ),     INZ=1,NONZR,1),  & 
   ( cbuf_vect%p_rA2D(1,INZ), cbuf_vect%p_rA2D(2,INZ), INZ=1,NONZR,1) 

  else

!   write(LUPRI,*) 'going to read indk_vect%p_iA1D/indl_vect%p_iA1D/cbuf_vect%p_rA2D ...'

   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,                       &
   ( indk_vect%p_iA1D(INZ), indl_vect%p_iA1D(INZ), INZ=1,NONZR,1),      & 
   (                      cbuf_vect%p_rA2D(1,INZ), INZ=1,NONZR,1) 

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
    write(LUPRI,'(2x,a,i6)') ' actual number of read integrals :',MOINT
    CALL FLSHFO(LUPRI)
  endif

if (print_int) then

 do i=1,NONZR
   if (.NOT.is_complex) then
   !  write(LUPRI,'(2x,a,2i5,f18.12)') 'k,l, real integral:',  & 
   !        indk_vect%p_iA1D(i),indl_vect%p_iA1D(i),cbuf_vect%p_rA2D(1,i)

   !  write(LUPRI,'(2x,a,4i4)') 'I,J,K,L:',  & 
   !        KR(IKR),KR(JKR),KR(indk_vect%p_iA1D(i)),KR(indl_vect%p_iA1D(i))
   else
   !  write(LUPRI,*) 'k,l,cmplx int:',indk_vect%p_iA1D(i),indl_vect%p_iA1D(i),  &
   !                  cbuf_vect%p_rA2D(1,i), cbuf_vect%p_rA2D(2,i)
   endif

   if (.NOT.is_complex) then
 ! ... unique ,KU indexes
      write(LUPRI,'(2x,a,i3,a,i3,a,i3,a,i3,a,f18.12)') '<',KR(IKR),' ',KR(JKR),'|',  & 
           KR(indk_vect%p_iA1D(i)),' ',KR(indl_vect%p_iA1D(i)),'>=', &
           cbuf_vect%p_rA2D(1,i)
      CALL FLSHFO(LUPRI)
   else
      write(LUPRI,'(2x,a,i3,a,i3,a,i3,a,i3,a,2f18.12)') '<',KR(IKR),' ',KR(JKR),'|',  & 
           KR(indk_vect%p_iA1D(i)),' ',KR(indl_vect%p_iA1D(i)),'>=', &
           cbuf_vect%p_rA2D(1,i), cbuf_vect%p_rA2D(2,i)
      CALL FLSHFO(LUPRI)
   endif

 enddo

endif
   
    IF (IKR.EQ.0) THEN
         IF (JKR.EQ.0) THEN
            NONZR = -1
         ELSE
            READ (MDINT) NEXTFIL
            CLOSE (MDINT,STATUS='KEEP')
            OPEN(MDINT,FILE=NEXTFIL,FORM='UNFORMATTED')
            write(LUPRI,'(2x,a,a)') 'HSCC: proceeding to the next file ',FILE
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
  !  ... lu units of dirac prepared files for HSMRCC
integer, parameter :: ludh = 12, lu_vind=13, lu_vint=14, lu_fhh=15, lu_fhp=16, lu_fpp=17 

integer ::  NTMO=0, NCMO=0, NVMO=0, NEMO=0, NHOLE=0, NPTCL=0, NELE=0
integer :: NHMO=0, NPMO=0

integer :: NFRZV, NACT

real(kind=8) :: TotEn, NucEn, TrfTol

integer, public :: iprint_hscc = 0

! ... stuff from HSMRCC code for preparing VINF,VIND,VINT, FHH, FHP, FPP files ...
integer, parameter :: NVFIL  = 14

integer, dimension(NVFIL) :: NumRec, RecLen
integer, dimension(NVFIL) :: StartRec, EndngRec
integer, dimension(NVFIL) :: TotNZInt, NumBlock

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine prepare_dirac_HSMRCC_file
!----------------------------------------------------------------------------------
!
! prepare files with DIRAC data (orbitals, symmetry, number of electrons..etc) for the 
! standalone HSMRCC run
!
!                       content of VINF file  ....
!
   ! ...............................................................................................
   ! ...  Rec     Module    : Variables present in the record                                    ...
   ! ...............................................................................................
   ! ...  01 [ Parameters ] : NCMO, NHOLE, NPTCL, NEMO                                           ...
   ! ...  02 [ Options    ] : OPTION                                                             ...
   ! ...  03 [ ConvMod    ] : MAXITER, CONVTOL                                                   ...
   ! ...  04 [ IntGenMod  ] : CUTOFF                                                             ...
   ! ...  05 [ SymInf     ] : point_group, naxis, norder, numirrep, irrep_product,               ...
   ! ...                      totally_symmetric_irrep, irrep_label                               ...
   ! ...  06 [ Shift      ] : CoreShiftVal, VirtShiftVal                                         ...
   ! ...  07 [ GMSMod     ] : NAO, NSO, NFRZV, TotEn, NucEn,                                     ...
   !                          TrfTol, IUNIQ, LABELSIZE, NINTMX                                   ...
   ! ...  08 [ GMSMod     ] : SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb,                         ...
   !                          OrbType, OrbSpace, OrbEnergy, OrbCoeff, HCoreMatrix                ...
   ! ...  09 [ IOModule   ] : NumRec, RecLen                                                     ...
   ! ...  10 [ IOModule   ] : StartRec(V2H2P), EndngRec(V2H2P), TotNZInt(V2H2P), NZIntRec_V2H2P  ...
   ! ...  11 [ IOModule   ] : StartRec(VPPPP), EndngRec(VPPPP), TotNZInt(VPPPP), NZIntRec_VPPPP  ...
   ! ...  12 [ IOModule   ] : StartRec(VHHHH), EndngRec(VHHHH), TotNZInt(VHHHH), NZIntRec_VHHHH  ...
   ! ...  13 [ IOModule   ] : StartRec(VHHPP), EndngRec(VHHPP), TotNZInt(VHHPP), NZIntRec_VHHPP  ...
   ! ...  14 [ IOModule   ] : StartRec(V2HPP), EndngRec(V2HPP), TotNZInt(V2HPP), NZIntRec_V2HPP  ...
   ! ...  15 [ IOModule   ] : StartRec(V2PHH), EndngRec(V2PHH), TotNZInt(V2PHH), NZIntRec_V2PHH  ...
   ! ...  16 [ IOModule   ] : StartRec(VPHHP), EndngRec(VPHHP), TotNZInt(VPHHP), NZIntRec_VPHHP  ...
   ! ...  17 [ IOModule   ] : StartRec(VPH2P), EndngRec(VPH2P), TotNZInt(VPH2P), NZIntRec_VPH2P  ...
   ! ...  18 [ IOModule   ] : StartRec(VPHPP), EndngRec(VPHPP), TotNZInt(VPHPP), NZIntRec_VPHPP  ...
   ! ...  19 [ IOModule   ] : StartRec(V2HHP), EndngRec(V2HHP), TotNZInt(V2HHP), NZIntRec_V2HHP  ...
   ! ...  20 [ IOModule   ] : StartRec(VHHHP), EndngRec(VHHHP), TotNZInt(VHHHP), NZIntRec_VHHHP  ...
   ! ...  21 [ IOModule   ] : StartRec(FHH)  , EndngRec(FHH)  , TotNZInt(FHH)                    ...
   ! ...  22 [ IOModule   ] : StartRec(FHP)  , EndngRec(FHP)  , TotNZInt(FHP)                    ...
   ! ...  23 [ IOModule   ] : StartRec(FPP)  , EndngRec(FPP)  , TotNZInt(FPP)                    ...
   ! ...............................................................................................
!
!   
!
!modules/iomodule.F90:   integer, dimension(NVFIL) :: NumRec, RecLen
!modules/iomodule.F90:   NumRec(V2H2P) = NHMO*NPMO
!modules/iomodule.F90:   NumRec(VPPPP) = NPMOSQ
!modules/iomodule.F90:   NumRec(VHHHH) = NHMOSQ
!modules/iomodule.F90:   NumRec(VHHPP) = NHMOSQ
!modules/iomodule.F90:   NumRec(V2HPP) = NHMO*NHMO
!modules/iomodule.F90:   NumRec(V2PHH) = NPMO*NPMO
!modules/iomodule.F90:   NumRec(VPHHP) = NHMO*NPMO
!modules/iomodule.F90:   NumRec(VPH2P) = NPMO*NPMO
!modules/iomodule.F90:   NumRec(VPHPP) = NPMO*NHMO
!modules/iomodule.F90:   NumRec(V2HHP) = NHMO*NHMO
!modules/iomodule.F90:   NumRec(VHHHP) = NHMOSQ
!modules/iomodule.F90:   NumRec(FHH)   = NHMO
!modules/iomodule.F90:   NumRec(FHP)   = NHMO
!modules/iomodule.F90:   NumRec(FPP)   = NPMO
!
!modules/iomodule.F90:   integer, dimension(NVFIL) :: NumRec, RecLen
!modules/iomodule.F90:   integer, dimension(NVFIL) :: StartRec, EndngRec
!modules/iomodule.F90:   integer, dimension(NVFIL) :: TotNZInt, NumBlock
!modules/parameters.F90:   integer, parameter :: NVFIL  = 14
!
!----------------------------------------------------------------------------------
use Memory
use Vectors
use io_moltra
#include "implicit.h"
#include "priunit.h"
!character(len=22) :: dirac_file='Dirac_for_HSMRCC.txt'
!character(len=*) :: dirac_file='Dirac_HSMRCC'
!character :: dirac_file='Dirac_HSMRCC'
character(len=12) :: dirac_file='Dirac_HSMRCC'
character(len=4) :: repna_temp_a, repna_temp_b
integer :: multbs(1:8,1:8)

type(Vect) :: SpinOfOrb,OccGMSRef,IrrepOrb,IrpLabOrb,OrbType,OrbSpace,OrbEnergy
 
if (iprint_hscc >= 5) then
  call Header('Output from "prepare_dirac_HSMRCC_file"',-1)
endif

! ... check if memory block is present
if (.not.IsMemorySticked() ) then
 write(LUPRI,*) 'memory array not present !' 
 call Quit(': memory array not present !  ')
else
 write(LUPRI,*) ' memory sticked !' 
endif

SpinOfOrb  = MakeVector('INTE',NSP,'GAMESS-like SpinOfOrb')
OccGMSRef  = MakeVector('INTE',NSP,'GAMESS-like OccGMSRef')
IrrepOrb   = MakeVector('INTE',NSP,'GAMESS-like IrrepOrb')
IrpLabOrb  = MakeVector('INTE',NSP,'GAMESS-like IrpLabOrb')
OrbType    = MakeVector('INTE',NSP,'GAMESS-like OrbType')
OrbSpace   = MakeVector('INTE',NSP,'GAMESS-like OrbSpace')

OrbEnergy  = MakeVector('REAL',NSP,'GAMESS-like OrbEnergy')

write(LUPRI,*)

nelec_tot = 0
do irp=1,NSYMRP
  write(LUPRI,'(2x,i2,a,i2,a,i2,a)') irp,'/',NSYMRP,' -> ',NELEC(irp),' electrons'
  nelec_tot = nelec_tot + NELEC(irp)
enddo

write(LUPRI,*) 'nelec_tot=',nelec_tot
write(LUPRI,*)

numirrep_dirac = GetNBSYM()
write(LUPRI,*) 'numirrep_dirac=',numirrep_dirac

! ... fill irreps now ...
write(LUPRI,'(/,2x,a,i4,a)') 'Irreps (1-',NSP,' )'
do i=1, NSP
  write(LUPRI,'(2x,i3,a,2i2,a,a4,a,f12.7)')   &  
  i,': irrep no ',IREPSP(i),IREPSP(i+NSP),'- irrep label:',REPNA(IREPSP(i+NSP)),' orb.energy:',EPS1(i) 
enddo
write(LUPRI,*)

! once again...
do i=1, NSP/2
 write(LUPRI,'(2x,i3,1x,a4,1x,a4)') i,REPNA(IREPSP(i+NSP)),REPNA(IREPSP(i+NSP/2+NSP))
enddo

write(LUPRI,*)

! fill orb.en.
ii = 1
do i=1, NSP/2

! write(LUPRI,*) 'EPS1(',i,')=',EPS1(i),EPS1(i+NSP/2)
 OrbEnergy%p_rA1D(ii)   = EPS1(i) 
 OrbEnergy%p_rA1D(ii+1) = EPS1(i+NSP/2)
 
 repna_temp_a = REPNA(IREPSP(i+NSP))
 repna_temp_b = REPNA(IREPSP(i+NSP/2+NSP))

 !write(LUPRI,'(2x,i3,1x,a4,1x,a4)') i,REPNA(IREPSP(i+NSP)),REPNA(IREPSP(i+NSP/2+NSP))
 write(LUPRI,'(2x,i3,1x,a4,i2,1x,a4,i2)') i, repna_temp_a,IREPSP(i), repna_temp_b, IREPSP(i+NSP/2)

 IrrepOrb%p_iA1D(ii)    = IREPSP(i+NSP)
 read(repna_temp_a,'(a4)') IrpLabOrb%p_iA1D(ii)

! IrrepOrb%p_iA1D(ii+1)  = IREPSP(i+NSP/2+NSP)
 IrrepOrb%p_iA1D(ii+1)  = IREPSP(i+NSP) ! previously saved irrep...

 read(repna_temp_b,'(a4)') IrpLabOrb%p_iA1D(ii+1)

 SpinOfOrb%p_iA1D(ii)   = +1
 SpinOfOrb%p_iA1D(ii+1) = -1

! write(LUPRI,'(a,2i4,2f12.7)') 'orb_en:',ii,ii+1,OrbEnergy%p_rA1D(ii),OrbEnergy%p_rA1D(ii+1)

   ii = ii + 2

enddo

do i=1,NSP
 if (i<=nelec_tot) then
  OccGMSRef%p_iA1D(i) = 1
 else
  OccGMSRef%p_iA1D(i) = 0
 endif
!
 if (i<=NCMO) then
  OrbType%p_iA1D(i) = -1
  OrbSpace%p_iA1D(i) = -1
 else if (i>NCMO .and. i<=NCMO+NACT) then
  OrbType%p_iA1D(i) = 0
  OrbSpace%p_iA1D(i) = 1
 else if (i>NCMO+NACT) then
  OrbType%p_iA1D(i) = 1
  OrbSpace%p_iA1D(i) = -1
 else
  call quit('dead branch !')
 endif

 if (iprint_hscc >= 5) then
   write(LUPRI,876) i,SpinOfOrb%p_iA1D(i),OccGMSRef%p_iA1D(i),  & 
                    OrbType%p_iA1D(i),OrbSpace%p_iA1D(i),       & 
                    IrrepOrb%p_iA1D(i),IrpLabOrb%p_iA1D(i),OrbEnergy%p_rA1D(i)

   876 format(i3,1x,'spin=',i2,1x,'occ=',i2,1x,'orbtype=',i2,1x,'orbspace=',i2,1x, &
                   'irp=',i2,1x,' irplab=',a4,1x,'orb.en=',f12.7)
 endif

enddo

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

write(ludh) numirrep_dirac

write(ludh) irrep_product_dirac

!write(ludh) ( SpinOfOrb%p_iA1D(i), i=1, NSP )
write(ludh) ( SpinOfOrb%p_iA1D(1:NSP) )

write(ludh) ( OccGMSRef%p_iA1D(1:NSP) )

write(ludh) ( OrbType%p_iA1D(1:NSP) )

write(ludh) ( OrbSpace%p_iA1D(1:NSP) )

write(ludh) ( IrrepOrb%p_iA1D(1:NSP) )

write(ludh) ( IrpLabOrb%p_iA1D(1:NSP) )

!write(ludh) ( OrbEnergy%p_rA1D(i), i=1, NSP )
write(ludh) ( OrbEnergy%p_rA1D(1:NSP) )

close(ludh,status="keep")

write(LUPRI,*) 
write(LUPRI,*) 'saved NELEC_TOT=',nelec_tot
write(LUPRI,*) 'saved NCORE=',NCORE
write(LUPRI,*) 'saved NACT=',NACT
write(LUPRI,*) 'saved NORBT=',NACT
write(LUPRI,*) 'saved NSP=',NSP
write(LUPRI,*) 'saved numirrep_dirac=',numirrep_dirac

if (iprint_hscc >= 5) then
    write(LUPRI,'(/,2x,a,a)')   & 
    '... (unformatted) Dirac input file for the HSMRCC run ready under name:',dirac_file
endif

call DeAllocObject(OrbEnergy)
call DeAllocObject(OrbSpace)
call DeAllocObject(OrbType)
call DeAllocObject(IrpLabOrb)
call DeAllocObject(IrrepOrb)
call DeAllocObject(OccGMSRef)
call DeAllocObject(SpinOfOrb)

!call Quit('prepare_dirac_HSMRCC_file')
end subroutine prepare_dirac_HSMRCC_file


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine prepare_dirac_HSMRCC_int_files
!-------------------------------------------------------------------------
!
!
   ! ... Generate Integrals ...

   ! ... The order of integral
   ! ... generation is the order
   ! ... in which NZIntRec info
   ! ... is written into FileVINF
   ! ... is the same as order in which
   ! ... Integrals Info is written to
   ! ... FileVIND & FileVINT
   ! ... FileFHH, FileFHP & FileFPP
   !
   ! ... The order of generation
   !
   !   Two-body integral files
   !
   !       01. V2H2P
   !       02. VPPPP
   !       03. VHHHH
   !       04. VHHPP
   !       05. V2HPP
   !       06. V2PHH
   !       07. VPHHP
   !       08. VPH2P
   !       09. VPHPP
   !       10. V2HHP
   !       11. VHHHP
   !
   !   One-body itegral files
   !
   !       01. FHH
   !       02. FHP
   !       03. FPP
   !     

!
!
!
!-------------------------------------------------------------------------
use Memory
use Vectors
use io_moltra
#include "implicit.h"
#include "priunit.h"
!character(len=22) :: dirac_file='Dirac_for_HSMRCC.txt'
!character(len=*) :: dirac_file='Dirac_HSMRCC'
!character :: dirac_file='Dirac_HSMRCC'
character(len=17) :: dirac_file_VIND='Dirac_HSMRCC.VIND'
character(len=17) :: dirac_file_VINT='Dirac_HSMRCC.VINT'
character(len=17) :: dirac_file_FHH='Dirac_HSMRCC.FHH'
character(len=17) :: dirac_file_FHP='Dirac_HSMRCC.FHP'
character(len=17) :: dirac_file_FPP='Dirac_HSMRCC.FPP'

!logical doprepare=.false.
logical:: doprepare=.true.


if (doprepare) then

  write(LUPRI,*) 'preparing DIRAC integral files for HSMRCC !!!'

open(lu_vind,file=dirac_file_VIND,form='unformatted',status='unknown', access='sequential') 
open(lu_vint,file=dirac_file_VINT,form='unformatted',status='unknown', access='sequential') 
open(lu_fhh,file=dirac_file_FHH,form='unformatted',status='unknown', access='sequential') 
open(lu_fhp,file=dirac_file_FHP,form='unformatted',status='unknown', access='sequential') 
open(lu_fpp,file=dirac_file_FPP,form='unformatted',status='unknown', access='sequential') 

close(lu_vind)
close(lu_vint)
close(lu_fhh)
close(lu_fhp)
close(lu_fpp)

endif

end subroutine prepare_dirac_HSMRCC_int_files

end module Dirac_HSMRCC

