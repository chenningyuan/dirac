!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end


module DIRAC_HSFS
 ! use Memory
 ! use Vectors
! HSFS parameters read in through the DIRAC input
! convergence
integer :: iconv_hsfs=8
! print level
integer :: iprint_hsfs=0
integer, parameter :: luhsfs=18
! ... depends only on integral file (no transformation etc..)
logical :: i_file_based=.false.
! array for HSFS indexes
!type(Vect), save :: hsfs_indx, boson_irps_indx, hsfs_orben
integer, save, allocatable :: hsfs_indx(:), boson_irps_indx(:)
real*8,  save, allocatable :: hsfs_orben(:)
contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine open_dirac_hsfs_file
!----------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------
implicit none
#include "priunit.h"
open(luhsfs,file="DIRAC_hsfs",form="unformatted",status="unknown")
write(LUPRI,*) 'DIRAC file was open.'
end subroutine open_dirac_hsfs_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine close_dirac_hsfs_file
!----------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------
implicit none
#include "priunit.h"
close(luhsfs,status="keep")
write(LUPRI,*) 'DIRAC file was closed.'
end subroutine close_dirac_hsfs_file
end module DIRAC_HSFS

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine PAMHSFSCC
!----------------------------------------------------------------------------
! Main launching subroutine for HSFSCC methods from DIRAC
!
! called from main/dirac.F
!----------------------------------------------------------------------------
use memory_allocator
use allocator_cfg

implicit none
real(8), allocatable :: WORK(:)
integer :: lwork, WORK_memsize

LWORK = WORK_memsize
call alloc(WORK,LWORK,id='WORK in PAMHSFSCC')

call HSFS_STUFF(WORK,LWORK)
call HIGH_SECTOR_FSCC_METHODS
call dealloc(WORK,id='WORK in PAMHSFSCC')

end subroutine PAMHSFSCC

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSFS_STUFF(WORK,LWORK)
!----------------------------------------------------------------------------------
!
!  Read DIRAC-MOLTRA produced files MRCONEE/MDCINT and prepare
!  DIRAC_hsfs file containing all data for the subsequent HSFS run.
! 
!  If MOLTRA section is skipped, or if .I_FILE is active in HSFS input, do
!  not prepare the DIRAC_hsfs file
! 
!  Written by Miro Ilias, Tel Aviv/Prievidza 2009
!----------------------------------------------------------------------------------
!use Memory
use io_moltra
use DIRAC_HSFS
implicit none
#include "priunit.h"
integer, intent(in) :: LWORK
real*8, intent(in) :: WORK(LWORK)
! ... skip if specified
if (i_file_based) then
 write(LUPRI,'(/,2x,a)') '*** No HSFS data preparation from DIRAC-MOLTRA files. ***'
 write(LUPRI,'(/,2x,a,/)') 'Assuming interface file DIRAC_hsfs is provided. '
 return
endif
call print_HSFS_Dirac_info
!call StickMemory(WORK,LWORK)
HSFS_mode=.true.
call read_mrconee;  call print_mrconee_data
call scan_mdcint
call read_mdcint_hsfs
!call UnstickMemory
end subroutine HSFS_STUFF


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_HSFS_Dirac_info
!-----------------------------------------------------------------------------
!
!  Print-out basic info on the Dirac_HSMRCC interface run
!
!  Called from:
!
! Written by Miro ILIAS, Tel Aviv 2009
!
!-----------------------------------------------------------------------------
implicit none
#include "priunit.h"

write(LUPRI,9001)

9001 format('                                                                 ',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '|                    I N T E R F A C E                          |',/, &
            '|                                                               |',/, &
            '| FOR HIGH-SECTORS FOCK SPACE COUPLED CLUSTER METHODS           |',/, &
            '|    (in spinfree mode, boson symmetry)                         |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '| Written by Miro ILIAS, Tel Aviv 2009,2012                     |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' )

call FLSHFO(LUPRI)

end subroutine print_HSFS_Dirac_info

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSFSINP(WORD)
!-----------------------------------------------------------------------------
!
!  Read parameters for the external HSFS programs, which are passed 
!  through the formatted input file, together with other DIRAC input
! 
!  Called from PSIINP (dirrdn.F), entering WORD is label currently read
!
!  Written by Miro ILIAS, Tel Aviv, 2009, Banska Bystrica, 2010
!
!-----------------------------------------------------------------------------
  use DIRAC_HSFS
  use print_level
  use high_sect_execs
implicit none
#include "priunit.h"
CHARACTER(len=7), INTENT(INOUT) :: WORD
Logical :: PROPER_POS, iprint_read=.false., rdin_read=.false.
CHARACTER:: PROMPT*1=' ',PROMPT_M*1=' ' 
CHARACTER(len=7) :: METHOD, KEYWORD=' '
! ... nice F90 version 
character(len=7) :: TABLE_TEMP, TABLE(10) =  & 
  (/'.PRINT ','.NACT  ','.NCORE ','.NOHSD ','.DRYRUN', & 
    '.CONV  ','.METHOD','.I_FILE','.INPUT ','.XXXXXX' /)
integer :: i
!
!  Initializations...
!

! ... entering label MUST be the beginning of *HSFS input list
PROPER_POS = (WORD .EQ. '*HSFS')

IF (PROPER_POS) THEN
Do 
  Read(LUCMD,'(A7)',ERR=910, END=930) WORD; PROMPT = WORD(1:1);  ! print *,'WORD=',WORD
  KEYWORD = WORD 
! Now decide !
  If (PROMPT == '*' ) Then
    !write(LUPRI,*) 'HSFSINP: leaving reading cycle !'
    Exit
  Else If ((PROMPT == '#') .or. (PROMPT == '!')) Then
   ! print *,'comment catched !', WORD
    Cycle ! ... continue the cycle..
  Else If (PROMPT == '.') Then
    ! ... keywords are beginning ...
    If (Any(TABLE(1:size(TABLE)) == WORD )) Then
     Select Case (WORD)
      Case ('.OPTION')
        Call PRTAB(size(TABLE), TABLE, KEYWORD//' input keywords',LUPRI)
      Case ('.PRINT')
        READ(LUCMD,*) iprint_hsfs
        write(LUPRI,*) ' read IPRINT_hsfs=',iprint_hsfs
        iprint_read = .true.
       Case ('.CONV')
         READ(LUCMD,*) iconv_hsfs
         write(LUPRI,*) ' read iconv_hsfs=',iconv_hsfs
       Case('.METHOD')
         write(LUPRI,*) ' going to read HSFS methdos !'
          call read_hsfs_methods(WORD)
!      Case ('.DRYRUN')
!        nohsccfiles=.true.
      Case('.INPUT')
         write(LUPRI,*) ' going to read HSFS "input" data ...'
         call read_rdin_input(WORD)
         rdin_read = .true.
      Case ('.I_FILE')
          i_file_based=.true.
      Case Default
        Write(LUPRI,'(/,3A,/)') ' Keyword "',KEYWORD, &
        '" is not yet implemented in HSFSINP.'
     End Select
    EndIf
  Else If (rdin_read) Then
    write(LUPRI,*) 'rdin was read...'
  Else
    WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,'" not recognized in HSFSINP.'
    CALL PRTAB(size(TABLE),TABLE,WORD//' input keywords',LUPRI)
    CALL QUIT('Illegal keyword in HSFSINP.')
  EndIf

EndDo
ELSE
 CALL QUIT('HSFSINP: entering label MUST be *HSFS !')
ENDIF

!  ... print out read parameters ... 

 WRITE (LUPRI,'(//)')
 CALL PRSYMB(LUPRI,'=',75,0)
 WRITE(LUPRI,'(1X,A)') 'HSFSINP: DIRAC set-up for the external standalone HSFS codes'
 CALL PRSYMB(LUPRI,'=',75,0)

   WRITE(LUPRI,'(/,1X,A,I3)') '* General print level   : ',iprint_hsfs
 if (iprint_read) then
   call set_print_level(iprint_hsfs)
   WRITE(LUPRI,'(/,2X,A)') 'print level set in this procedure'
 else
   WRITE(LUPRI,'(/,2X,A)') 'print level to be set after reading "input"'
 endif

! if (nohsccfiles) then

 WRITE(LUPRI,'(1X,A,I3)')   '* Convergence           : ',iconv_hsfs

 write(LUPRI,'(1X,A,I3)') '* Number of individual HSFS methods: ',num_methods
 if (num_methods == 0) then
  write(LUPRI,'(2X,A)') '- no method called'
 else
  write(LUPRI,'(2X,A)') '- list of methods:'
   do i=1,num_methods
     write(LUPRI,*) i,method_name(i)
   enddo
 endif

 if (iprint_read) then
   call set_print_level(iprint_hsfs)
 else

 endif

! if (nohsccfiles) then
!  WRITE(LUPRI,'(1X,A)') '* Do not prepare HSCC files (dry run). You might take them from GAMESS.'
! else
!  WRITE(LUPRI,'(1X,A)') '* Prepare HSCC files by Dirac.'
! endif

! if (nohsccdirac) then
!  WRITE(LUPRI,'(1X,A)') '* Only prepare HSCC files by DIRAC, do not continue with "afterDIRAC" HSCC.'
! endif

! call quit('hsfs input reading')

 Return

910  CONTINUE
        CALL QUIT('Error reading LUCMD, no "*HSFS" input found !')
920  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSFS" input some error occured...')
930  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSFS" input end of file reached...')

end subroutine HSFSINP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_hsfs_methods(WORDIN)
!-----------------------------------------------------------------------------
!
!    read the sequence of individual HSFS methods which are to be run
!
!    called from HSFSINP
!
!-----------------------------------------------------------------------------
 use high_sect_execs
implicit none
#include "priunit.h"
CHARACTER(len=7) :: WORDIN, METHOD, KEYWORD=' '
CHARACTER:: PROMPT*1=' ',PROMPT_M*1=' ' 
integer :: im

IF (.NOT.WORDIN == '.METHOD') THEN
   WRITE(LUPRI,*) 'WORDIN=',WORDIN
   CALL QUIT('read_hsfs_methods: wrong entering word ! Must be .METHOD ')
ELSE
  im = 0
  Do
           Read(LUCMD,'(A7)') METHOD;  PROMPT_M = METHOD(1:1)
           !write(LUPRI,*) 'HSFS methods; METHOD=',METHOD
           !If (PROMPT_M == '*' .or. PROMPT_M == '.' ) Then
           If (METHOD == '.END') Then
             !write(LUPRI,*) 'HSFS methods leaving'
             Exit
           Else IF (PROMPT_M == '#' .or. PROMPT_M == '!') Then
             Cycle
           Else If (Any(hs_method_name(1:size(hs_method_name)) == METHOD )) Then
               im = im + 1
               if (im > max_methods) then
                write(LUPRI,*) 'im, max_methods=',im, max_methods
                call quit('number of methods exceeded ! increase max_methods' &
               //' parameter !')
               endif
              ! write(LUPRI,*) 'valid method name=',METHOD
               method_name(im)=METHOD 
               write(LUPRI,*) im,method_name(im)
           Else
             WRITE (LUPRI,'(/,3A,/)') ' Keyword "',METHOD,'" not recognized in HSFSMETHODS.'
             CALL PRTAB(size(hs_method_name),hs_method_name,METHOD//' input keywords',LUPRI)
             CALL QUIT('Illegal keyword in HSFSMETHODS.')
           EndIf
  EndDo
  ! store number of methods
   num_methods = im
ENDIF
end subroutine read_hsfs_methods

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_rdin_input(WORD)
!------------------------------------------------------------------------------
!      ... for reading "input" from within DIRAC input file ...
!------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
      use sym_mod
      use dpd_mod
#else
!     implicit real*8 (a-h,o-z)
#include "dat.h"
!      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
#include "dpd.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#endif

implicit none
#include "priunit.h"
CHARACTER(len=7), INTENT(IN) :: WORD
CHARACTER(len=7) :: WORDIN, METHOD, KEYWORD=' '
CHARACTER:: PROMPT*1=' ',PROMPT_M*1=' '
#if defined HSFS_DYN_ALLOC
   call init_dat_mod
   call init_sym_mod
   call init_dpd_mod
#endif
   call rdin
   call rdsym
end subroutine read_rdin_input


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_HSFS_indexes()
!----------------------------------------------------------------------------------
!
!                 assumes MRCONEE is read !
!
!----------------------------------------------------------------------------------
!use Memory
!use Vectors

use DIRAC_HSFS
use io_moltra
implicit none
#include "priunit.h"
logical :: indx_found
  integer, external :: GetNBSYM
  integer :: isp, icheck,i,ik,j,k,nirp,ix
  integer :: istat

!  ... allocate vectors ...
!hsfs_indx=MakeVector('INTE',NSP,'HSFS indexes')
  allocate (hsfs_indx(NSP),stat=istat);  call alloc_stat_check("hsfs_indx(NSP)",istat)
  write(LUPRI,*)  'get_HSFS_indexes: hsfs_indx vector allocated'
!boson_irps_indx=MakeVector('INTE',GetNBSYM(),'array for boson irreps')
allocate (boson_irps_indx(GetNBSYM()),stat=istat); call alloc_stat_check("boson_irps_indx(GetNBSYM())",istat)
write(LUPRI,*)  'get_HSFS_indexes: boson_irps_indx vector allocated, GetNBSYM()=',GetNBSYM()
call izero(boson_irps_indx,GetNBSYM())

do isp=1, nsp
  if (iprint_hsfs >= 2) then
   write (LUPRI,'(2x,a,i4,a,2i3,a,f18.10,a,a4)') 'ISP=',ISP,' IREPSP(ISP/NSP+ISP):', &
   IREPSP(ISP),IREPSP(NSP+ISP),                                                      &
   ' EPS1(ISP)=',EPS1(ISP),'--> ',REPNA(IREPSP(NSP+ISP))
  endif

  if (IREPSP(ISP).eq.1.or.IREPSP(ISP).eq.3) then
     !boson_irps_indx%p_iA1D( IREPSP(NSP+ISP) ) = boson_irps_indx%p_iA1D( IREPSP(NSP+ISP) ) + 1 
     boson_irps_indx( IREPSP(NSP+ISP) ) = boson_irps_indx( IREPSP(NSP+ISP) ) + 1 

     if (iprint_hsfs >= 2) then
       !write(LUPRI,*) '- irrep:',boson_irps_indx%p_iA1D( IREPSP(NSP+ISP) )
       write(LUPRI,*) '- irrep:',boson_irps_indx( IREPSP(NSP+ISP) )
      endif
  endif
enddo

! do the check
icheck=0
do i=1,GetNBSYM()
 if (iprint_hsfs >= 6) then
    !write(LUPRI,*) i,'irrep:',boson_irps_indx%p_iA1D(i)
    write(LUPRI,*) i,'irrep:',boson_irps_indx(i)
 endif
 !icheck=icheck+boson_irps_indx%p_iA1D(i)
 icheck=icheck+boson_irps_indx(i)
enddo

if (icheck /= nsp/2) then
 write(LUPRI,*) 'sum of all spinors, icheck=',icheck
 write(LUPRI,*) '               while nsp/2=',nsp/2
 call print_mrconee_data
 call quit('wrong number of spinors !')
endif

if (iprint_hsfs >= 2) then
write(LUPRI,*) 'filled irreps'
do i=1,GetNBSYM()
! boson_irps_indx%p_iA1D(i)=0
 !write(LUPRI,*) i,REPNA(i),' irrep:',boson_irps_indx%p_iA1D(i)
 write(LUPRI,*) i,REPNA(i),' irrep:',boson_irps_indx(i)
enddo
endif

! ... fill "hsfs_indx" vector: map order MRCONEE -> HSFS
ik=0
do i=1,GetNBSYM()
do j=1,nsp

 if ( i == IREPSP(NSP+j) ) then
 !  write(LUPRI,'(a,3i3)') 'i,j,IREPSP(NSP+j):',i,j,IREPSP(NSP+j)
   k=ik
   k=k+1
      !hsfs_indx%p_iA1D(j) = k
      hsfs_indx(j) = k
   ik=k
 endif

enddo
enddo

write(LUPRI,'(/,2x,a)') '* HSFS - MOLTRA MO-integrals index mapping done.'

do isp=1,nsp
if (IREPSP(ISP).eq.1.or.IREPSP(ISP).eq.3) then
  if ( iprint_hsfs >= 2) then
    write (LUPRI,'(2x,a,i4,a,2i3,a,f18.13,a,a4,a,i4)') 'ISP=',ISP,' IREPSP(ISP/NSP+ISP):', &
            IREPSP(ISP),IREPSP(NSP+ISP), &
    !' EPS1(ISP)=',EPS1(ISP),'-->',REPNA(IREPSP(NSP+ISP)),' hsfs_indx:',hsfs_indx%p_iA1D(isp)
    ' EPS1(ISP)=',EPS1(ISP),'-->',REPNA(IREPSP(NSP+ISP)),' hsfs_indx:',hsfs_indx(isp)
  endif
endif
enddo

call open_dirac_hsfs_file
rewind (luhsfs)
write(LUPRI,'(/,2x,a)') 'going to write luhsfs data...'
 ! ... save some stuff 
nirp=GetNBSYM()

write(luhsfs) nirp, nsp/2
write(LUPRI,'(2x,a,i2,a,i4)') 'written nirp=',nirp,' nort=',nsp/2
!write(luhsfs) ( boson_irps_indx%p_iA1D(i),i=1,nirp )
write(luhsfs) ( boson_irps_indx(i),i=1,nirp )
!write(LUPRI,'(1x,a,8i4,/)') 'written boson_irps:',( boson_irps_indx%p_iA1D(i),i=1,nirp )
write(LUPRI,'(1x,a,8i4,/)') 'written boson_irps:',( boson_irps_indx(i),i=1,nirp )

ix=0
do isp=1,nsp
  indx_found=.false.
  do j=1,nsp
     !if (isp == hsfs_indx%p_iA1D(j)) then
     if (isp == hsfs_indx(j)) then
       if ( iprint_hsfs >= 10 ) then
         !write(LUPRI,*) 'found:',isp,hsfs_indx%p_iA1D(j),j  
         write(LUPRI,*) 'found:',isp,hsfs_indx(j),j  
       endif
       indx_found=.true.
       Exit
     endif
  enddo
  if (indx_found) then
   ix = ix + 1
   write(luhsfs) EPS1( j )
   if ( iprint_hsfs >= 10 ) then
     write(LUPRI,*) isp,j,REPNA(IREPSP(NSP+j)),EPS1(j)
   endif
 endif
enddo

write(LUPRI,'(/,2x,a)') 'one-electron HSFS stuff written into luhsfs file'
if (ix.ne.nsp/2) then
write(LUPRI,'(/,2x,a,2i4)') 'one-electron HSFS stuff written into luhsfs file, numbers ix,nsp/2:',ix,nsp/2
  call quit('wrongly extracted HSFS indexes !')
endif

end subroutine get_HSFS_indexes

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_mdcint_hsfs
!-------------------------------------------------------------------------
!
!  Read of MDCINT integrals for HSFSCC
!
!-------------------------------------------------------------------------
use DIRAC_HSFS
use io_moltra
implicit none
#include "priunit.h"
logical :: print_batch_info, print_int, is_mdc, duplicity, still_duplicity
!type(Vect) :: indk_vect, indl_vect, cbuf_vect, mdcint_vect 
integer, allocatable :: indk_vect(:), indl_vect(:)
real*8, allocatable  ::  cbuf_vect(:), mdcint_vect(:)
 
!type(Array2Dim) :: ij_store
integer, allocatable :: ij_store(:,:)

integer :: i,j,k,l,nn,nnzi_save,ij,MOINT,iblock,icount,IKR,JKR,NONZR
integer :: nnnn, INZ, i1,i2,i3,i4,i1hsfs,i2hsfs,i3hsfs,i4hsfs
integer :: ipos,i1hsfs_p,i2hsfs_p,i3hsfs_p,i4hsfs_p
integer :: icbuf,ic
!integer, external :: Find_ijkl_pos
 interface
   integer function Find_ijkl_pos(i0,j0,k0,l0,ij_store)
     integer, intent (in) :: i0,j0,k0,l0
     integer, intent(in) :: ij_store(:,:)
   end function Find_ijkl_pos
 end interface
real*8 :: xint
integer :: istat

Call Header('MDCINT control print-out',-1)

! ... check if memory block is present
!if (.not.IsMemorySticked() ) then
! write(LUPRI,*) 'read_mdcint: memory array not present !' 
! call Quit('read_mdcint: memory array not present !  ')
!else
! write(LUPRI,*) 'read_mdcint: memory sticked !' 
!endif

print_batch_info = ( iprint_hsfs >= 8 )
       print_int = ( iprint_hsfs >= 14 )

!print_batch_info = .true.
!print_int = .true.

write(LUPRI,*) ' print MDCINT integrals batches info:',print_batch_info
write(LUPRI,*) '          print all MDCINT integrals:',print_int

! ... check if MRCONEE was read 
if (.NOT.mrconee_read) then
 write(LUPRI,*) 'read_mdcint: MRCONEE file was not read before !' 
 call Quit('read_mdcint: MRCONEE file was not read before !  ')
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
 !call scan_mdcint

 nnzi_save = nnzi
 write(LUPRI,*) ' number of ALL nonzero 2-el MO-integrals on MDCINT=',nnzi
 write(LUPRI,*) ' number of <I J|-integral BLOCKS on MDCINT=',nintbl

! 
!!! write(luhsfs) nnzi ! write total number of integrals

 write(LUPRI,'(A,$)') 'Integrals have been calculated on the'
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

!call alloc_mdcint_data
! call print_mem_param
! ... note: "ccdriv.F":1923 is using RCW !!!
!cbuf = Make2DimArray('REAL',RCW,nnzi,'MDCINT two-electron integrals') ! big array !!!!
!pos_ij = Make3DimArray('INTE',NKR,NKR,2,'position of <ij| int.blocks on integral array') 
!pos_kl = Make3DimArray('INTE',NKR,NKR,nintbl,'position of |kl> integrals in bloks on integral array') ! huge array 

!indk_vect = MakeVector('INTE',ilenght,'MDCINT "k" indexes')
  allocate(indk_vect(ilenght),stat=istat);  call alloc_stat_check("indk_vect(ilenght)",istat)
!indl_vect = MakeVector('INTE',ilenght,'MDCINT "l" indexes')
  allocate(indl_vect(ilenght),stat=istat);  call alloc_stat_check("indl_vect(ilenght)",istat)

nn=NKR*(NKR+1)/2
!ij_store = Make2DimArray('INTE',nn,2,'i,j-indexes of common ij')
  allocate(ij_store(nn,2),stat=istat);  call alloc_stat_check("ij_store(nn,2)",istat)

ij=0
do i=1,NKR
do j=1,i
 ij=ij+1
 !ij_store%p_iA2D(ij,1)=i; ij_store%p_iA2D(ij,2)=j
 ij_store(ij,1)=i; ij_store(ij,2)=j
enddo
enddo

nnnn=Find_ijkl_pos(NKR,NKR,NKR,NKR,ij_store)

write(LUPRI,*) 'ijkl batch, all,nnnn=',nnnn

!cbuf_vect   =   MakeVector('REAL',nintmaxbl,'MDCINT largest integral block')
allocate(cbuf_vect(nintmaxbl),stat=istat); call alloc_stat_check("cbuf_vect(nintmaxbl)",istat)

!mdcint_vect = MakeVector('REAL',nnnn,'all spinfree MDCINT integrals block')
allocate(mdcint_vect(nnnn),stat=istat);  call alloc_stat_check("mdcint_vect(nnnn)",istat)


  rewind(MDINT)
  read(MDINT)

 ! call ZeroObject(cbuf_vect)
 ! call ZeroObject(mdcint_vect)

  write(LUPRI,*) 'Going to read MDCINT transformed integrals'

  MOINT = 0
  iblock = 0
  icount = 0
 1 continue
! ... only real integrals !!!
   !READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,                               &
   !( indk_vect%p_iA1D(INZ), indl_vect%p_iA1D(INZ),          INZ=1,NONZR,1),     &
   !(                        cbuf_vect%p_rA1D(INZ),          INZ=1,NONZR,1) 

   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,     &
   ( indk_vect(INZ), indl_vect(INZ),          INZ=1,NONZR,1),     &
   (                        cbuf_vect(INZ),          INZ=1,NONZR,1) 



  if (NONZR > ilenght) then
    write(LUPRI,*) 'NONZR :',NONZR,' ilenght :',ilenght
    write(LUPRI,*) 'NONZR > ilenght; out of memory allocation !'
    call Quit('NONZR > ilenght; out of memory allocation !')
  endif

  MOINT = MOINT + NONZR

  if (print_batch_info) then
    write(LUPRI,'(2x,a,i5)')  '...read (IJ|-integral block:',iblock
    write(LUPRI,'(2x,a,3i4)') 'read IKR,JKR,NONZR:',IKR,JKR,NONZR
    write(LUPRI,'(2x,a,i6)')  ' actual number of read integrals :',MOINT
    !write(LUPRI,'(2x,a,i6)')  ' position of cbuf_vect array :',icbuf
    CALL FLSHFO(LUPRI)
  endif

if (print_int) then
 write(LUPRI,'(/,2x,a,2i4)') 'read ( I J |... indexes:',KR(IKR),KR(JKR)
 do INZ=1, NONZR
   i1=KR(IKR); i2=KR(JKR)
   !i3=KR(indk_vect%p_iA1D(INZ)); i4=KR(indl_vect%p_iA1D(INZ))
   i3=KR(indk_vect(INZ)); i4=KR(indl_vect(INZ))
   icount = icount + 1
   !  ... unique, KU indexes ... (i1 i2|i3 i4)=<i1 i3|i2 i4>
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,i3,a,i3,a,i3,a,i3,a)')  & 
           icount,'. (',KR(IKR),' ',KR(JKR),'|',                              & 
           !KR(indk_vect%p_iA1D(INZ)),' ',KR(indl_vect%p_iA1D(INZ)),')=',           &
           KR(indk_vect(INZ)),' ',KR(indl_vect(INZ)),')=',           &
           !cbuf_vect%p_rA1D(INZ),' = <',i1,' ',i3,'|',i2,' ',i4,'>'
           cbuf_vect(INZ),' = <',i1,' ',i3,'|',i2,' ',i4,'>'
     CALL FLSHFO(LUPRI)
 enddo
endif

! ... Extract integrals, indexes determined by array position ...

  do INZ=1, NONZR
   icount = icount + 1
   i1=KR(IKR); i2=KR(JKR)
   !i3=KR(indk_vect%p_iA1D(INZ)); i4=KR(indl_vect%p_iA1D(INZ))
   i3=KR(indk_vect(INZ)); i4=KR(indl_vect(INZ))
!  (i1 i2|i3 i4)=<i1 i3|i2 i4>
   !i1hsfs=hsfs_indx%p_iA1D(i1); i2hsfs=hsfs_indx%p_iA1D(i3)
   i1hsfs=hsfs_indx(i1); i2hsfs=hsfs_indx(i3)
   !i3hsfs=hsfs_indx%p_iA1D(i2); i4hsfs=hsfs_indx%p_iA1D(i4)
   i3hsfs=hsfs_indx(i2); i4hsfs=hsfs_indx(i4)
   call permut_indxs(i1hsfs  ,i2hsfs  ,i3hsfs   ,i4hsfs  ,  & 
                     i1hsfs_p,i2hsfs_p,i3hsfs_p ,i4hsfs_p)
    ipos= Find_ijkl_pos(i1hsfs_p,i2hsfs_p,i3hsfs_p,i4hsfs_p,ij_store) 

    if (ipos > nnnn) then
    ! write(LUPRI,*) 'i,j,ij:',i1hsfs_p,i2hsfs_p,ij
    ! write(LUPRI,*) 'k,l,kl:',i3hsfs_p,i4hsfs_p,kl
     write(LUPRI,*) 'ipos,nnnn=',ipos,nnnn
     call quit('wrong , ipos > nnnn !')
    endif

   ! store integral according to its unique indexes
    !mdcint_vect%p_rA1D(ipos)= cbuf_vect%p_rA1D(INZ)
    mdcint_vect(ipos)= cbuf_vect(INZ)

  if (print_int) then
    write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3,a,2i3,a,2i3,a,i4)')  & 
          icount,'. (',i1,' ',i2,'|',  & 
          i3,' ',i4,')=', &
          !cbuf_vect%p_rA1D(INZ),'-',i1hsfs,i2hsfs,i3hsfs,i4hsfs,' :p <',   &  
          cbuf_vect(INZ),'-',i1hsfs,i2hsfs,i3hsfs,i4hsfs,' :p <',   &  
          i1hsfs_p,i2hsfs_p,'|',i3hsfs_p,i4hsfs_p,'>; pos',ipos
  endif

  enddo

  iblock = iblock + 1    ! increment integral block

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
  ! call Quit('NONZR == -1!')
 !  RETURN
   GoTo 2
 endif

 if (print_int) then
  write(LUPRI,*) '... jumping to read the next integral batch in the file ... '
 endif

 CALL FLSHFO(LUPRI)
 GoTo 1

 2 continue

 write(LUPRI,'(/,2x,a,/)')  &
 ' all integrals from the MDCINT file read ...closing the file'
 CLOSE(MDINT,STATUS='KEEP')

 write(LUPRI,'(2x,a,/)') 'Going to to store unique MO 2-el. integrals into the DIRAC_hsfs file...'
 CALL FLSHFO(LUPRI)
 write(luhsfs) NKR
 ic=0
 do i=1,NKR
 do j=1,i
 do k=1,i
 do l=1,j
   ipos= Find_ijkl_pos(i,j,k,l,ij_store)
   !xint=mdcint_vect%p_rA1D(ipos)
   xint=mdcint_vect(ipos)
   if (dabs(xint).ge.1.0d-15) then
     ic = ic + 1
     write(luhsfs) i,j,k,l,xint
     if (print_int) then
       write(LUPRI,'(i4,a,4i4,f18.12)') ic,'.<ij|kl>:',i,j,k,l,xint
     endif
   endif
 enddo
 enddo
 enddo
 enddo

 write(luhsfs) 0,0,0,0,0.0d0
 write(LUPRI,*) 'saved nonzero integrals & their indexes :',ic
 CALL FLSHFO(LUPRI)

! ... deallocate ...
 !call DeAllocObject(mdcint_vect)
 !call DeAllocObject(cbuf_vect)
 !call DeAllocObject(ij_store)
 !call DeAllocObject(indl_vect)
 !call DeAllocObject(indk_vect)
 deallocate(mdcint_vect)
 deallocate(cbuf_vect)
 deallocate(ij_store)
 deallocate(indl_vect)
 deallocate(indk_vect)

 write(LUPRI,*) 'mdcint_vect, cbuf_vect, indl_vect, indk_vect released ...'; CALL FLSHFO(LUPRI)

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

end subroutine read_mdcint_hsfs


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine permut_indxs(i1hsfs  ,i2hsfs  ,i3hsfs  ,i4hsfs,  & 
                        i1hsfs_p,i2hsfs_p,i3hsfs_p,i4hsfs_p)
!----------------------------------------------------------------------------------------------
!       Put two-electron <pq|rs >integrals indexes into this order:
!
!    <p q|r s> indexes to be permuted so that:
!
!       p >= q and (p,q) >= (r,s)
!
! <p q|r s> = <p s|r q> = <r q|p s> = <r s|p q>
!
! <p q|r s> = <q p|s r> = <q r|s p> = <s r|q p> = <s p|q r>
!
! p,q,r,s=i1hsfs,i2hsfs,i3hsfs,i4hsfs (on entry)
!----------------------------------------------------------------------------------------------
implicit none
#include "priunit.h"
integer, intent(in)  :: i1hsfs  ,i2hsfs  ,i3hsfs  ,i4hsfs
integer, intent(out) :: i1hsfs_p,i2hsfs_p,i3hsfs_p,i4hsfs_p
integer :: ih

 i1hsfs_p=i1hsfs; i2hsfs_p=i2hsfs; i3hsfs_p=i3hsfs; i4hsfs_p=i4hsfs

if (i1hsfs_p < i3hsfs_p) then  ! <pq|rs>=<rq|ps>
  ih =i1hsfs_p; i1hsfs_p=i3hsfs_p; i3hsfs_p=ih
endif

if (i2hsfs_p < i4hsfs_p) then ! <pq|rs>=<ps|rq>
  ih =i2hsfs_p; i2hsfs_p=i4hsfs_p; i4hsfs_p=ih
endif

if (  (i1hsfs_p < i2hsfs_p .or. i3hsfs_p < i4hsfs_p)    ) then
 if ( (max(i1hsfs_p,i2hsfs_p) < max(i3hsfs_p,i4hsfs_p) ) & 
    .or.(i1hsfs_p < i2hsfs_p).or.((i1hsfs_p.eq.i2hsfs_p).and.(i3hsfs_p<i4hsfs_p))  ) then ! <pq|rs>=<qp|sr>
   ih = i1hsfs_p; i1hsfs_p=i2hsfs_p; i2hsfs_p=ih
   ih = i3hsfs_p; i3hsfs_p=i4hsfs_p; i4hsfs_p=ih
 endif
endif

! Check permuted indexes of <i1 i2|i3 i4>
if (.not.((i1hsfs_p>=i2hsfs_p).and.(max(i1hsfs_p,i2hsfs_p) >= max(i3hsfs_p,i4hsfs_p)))) then
  write(LUPRI,'(a,4i3,a,4i3)') 'permut:',i1hsfs,i2hsfs,i3hsfs,i4hsfs,'->',i1hsfs_p,i2hsfs_p,i3hsfs_p,i4hsfs_p
  call quit('permut_indxs: inproper permutations !')
endif

end subroutine permut_indxs

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_indxs(pq,p,q)
!--------------------------------------------------------
! returns p,q indexes, where pq = p(p-1)/2+q, p>=q
!
! finds them in cycle
!--------------------------------------------------------
implicit none
integer, intent(in) :: pq
integer, intent(out) :: p,q
   p = 0
   do
     p = p + 1
     q = pq - ( (p*(p-1))/2 )
     if (q<=p) exit
   enddo
return
end subroutine get_indxs

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function Find_ijkl_pos(i0,j0,k0,l0,ij_store) result(ipos)
!----------------------------------------------------------------
!
! based on <ij|kl> indexes find position on array (ipos)
!
! cycles:  i=1,n
!          j=1,i  (j<=i)
!          k=1,i  (k<=i)
!          l=1,j  
!
!      ...  ipos of <ij|kl> 
!
!----------------------------------------------------------------
 ! use Memory
 ! use Vectors
  implicit none
  integer, intent (in) :: i0,j0,k0,l0
  integer :: ipos,kl,i,j,is,ij,ij0
  !type(Array2Dim), intent(in) :: ij_store
  integer, intent(in) :: ij_store(:,:)

  kl=(k0-1)*j0 + l0

  ij0 = (i0*(i0-1)/2) + j0
  is =0
  do ij=1,ij0-1
    ! call get_indxs(ij,i,j)
    !i=ij_store%p_iA2D(ij,1); j=ij_store%p_iA2D(ij,2)
    i=ij_store(ij,1); j=ij_store(ij,2)
    is = is +(i*j)
  enddo
  ipos=is + kl ! resulting position on array depending of i,j,k,l indexes
end function Find_ijkl_pos

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_dirac_ints(nort,eps,nn,ind,v,nint)
!------------------------------------------------------------------------------
!
!  provide data (integrals,orbital energies) for high-sect CCSD from DIRAC
!  
!  nort - number of orbitals  
!  
!  nint - number of 2-el.integrals
!  
!  called from readv(ccsd3.F) -  HighSectFSCC method
!
!------------------------------------------------------------------------------
use print_level
use ipack_mod, only: ipack4
use DIRAC_HSFS
#if defined HSFS_DYN_ALLOC
use hsfs_parameters
#else
#include "param_ccsd.h"
#endif
implicit none
#include "priunit.h"
   real*8, intent(inout) :: eps(*), v(*)
   integer, intent(inout) :: ind(*),nint, nort
   integer, dimension(8) :: irreps
   integer :: ios,nirp,i,nintx,i1,i2,i3,i4,nn,nkr
   logical :: unread_file

   call header('Output from "get_dirac_ints"',-1)

!  call set_print_level(iprint_hsfs)
!  write(LUPRI,'(2x,a,i5)') 'getiprint=',getiprint()
   iprint_hsfs = getiprint()
!  write(LUPRI,'(2x,a,i5)') 'iprint_hsfs=',iprint_hsfs

   call open_dirac_hsfs_file; rewind(luhsfs)
   write(LUPRI,*) luhsfs,' file open'

   read(luhsfs,iostat=ios) nirp,nort
   if (ios.ne.0) then
     write(LUPRI,*) 'broken reading, nirp=',nirp,' nort=',nort
     call quit('error in luhsfs reading 1!')
   endif

   write(LUPRI,'(2x,a,i2,a,i4)') 'read nirp=',nirp,' nort=',nort

   read(luhsfs,iostat=ios) (irreps(i),i=1,nirp)
   if (ios.ne.0) call quit('error in luhsfs reading 2!')
   write(LUPRI,'(2x,a,8i4)') 'read irreps:',(irreps(i),i=1,nirp)

   write(LUPRI,'(2x,a,i5)') 'going to read orbital energies, nort=',nort

   do i=1, nort
    read(luhsfs,iostat=ios) eps(i)
    if (ios.ne.0) call quit('error in luhsfs reading 3 !')

    if (iprint_hsfs >= 3) then
      write(LUPRI,*) i,eps(i)
    endif

   enddo 

   read(luhsfs,iostat=ios) nkr
   if (ios.ne.0) call quit('error in luhsfs nint reading!')
   write(LUPRI,*) 'nkr=',nkr

   write(LUPRI,'(2x,a)') 'Going to read integrals & indexes:'

   nintx=0
   unread_file=.true.
   do while (unread_file) 
    nintx = nintx + 1
    if (nintx > ndint) then
      write(LUPRI,*) 'nintx,ndint:',nintx,ndint
      call quit('get_dirac_ints: number of integrals exceeds allocated arrays !')
    endif

    read(luhsfs,iostat=ios) i1,i2,i3,i4,v(nintx)
    if (ios.ne.0) call quit('error in luhsfs i1,i2,i3,i4,v(i) reading !')
    ind(nintx)=ipack4(i1,i2,i3,i4) 
    
    unread_file=.not.( (i1.eq.0).and.(i2.eq.0)) 
    if (iprint_hsfs >= 8) then
      write(LUPRI,'(i5,a,4i4,1x,a,i11,a,f18.12)') nintx,'.',i1,i2,i3,i4,'pck in:',ind(nintx),' int:',v(nintx)
    endif
   enddo
    
   nint=nintx-1
   write(LUPRI,*) 'read all integrals&indexes from the file, their number,  nint=',nint
   call close_dirac_hsfs_file

end subroutine get_dirac_ints

