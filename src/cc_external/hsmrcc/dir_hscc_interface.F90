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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSCC_STUFF(WORK,LWORK)
!----------------------------------------------------------------------------------
!
!  Called from PAMTRA, if defined (MOD_HSCC)
!
!  Written by Miro Ilias, Tel Aviv/Prievidza 2007,2008,2009
!
!----------------------------------------------------------------------------------
use Memory
use Dirac_HSMRCC
use io_moltra

implicit none
#include "priunit.h"
#include "dcbpsi.h"

integer,intent(in) :: LWORK
real*8, intent(in) :: WORK(LWORK)
logical :: check_point = .false.
integer :: ipos1,ipos2,ipos3,ipos4
logical :: check_MOints = .true.

! ... first run test_oof features
if (iprint_hscc >= 10) then
  call TEST_OOF(WORK,LWORK)
  call TEST_OO_MEMORY(WORK,LWORK)
endif

if (.not.DOHSCC) then
 write(LUPRI,'(2x,a)') 'DOHSCC=.false., no HSCC data preparation'
 return
endif

write(LUPRI,*)
write(LUPRI,*) 'in HSCC_STUFF !! LWORK=',LWORK

call print_HSMRCC_Dirac_info

call StickMemory(WORK,LWORK)

!if (check_MOints.and.iprint_hscc >= 10) then
! call read_mrconee; call print_mrconee_data
! call read_mdcint ! read MDCINT file
!endif

! read MRCONEE & MDCINT, prepare integral files
call read_moltra_files 

!call prepare_dirac_HSMRCC_file
!call prepare_dirac_HSMRCC_int_files

call UnstickMemory

end subroutine HSCC_STUFF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSCCINP(WORD,reset)
!-----------------------------------------------------------------------------
!
!  Read parameters for the external HSMRCC program, which are passed 
! through the formatted input file, together with other DIRAC input
! 
! Called from TRAINP if DOHSCC=.true., i.e. at the beginning
!
! Written by Miro ILIAS, Tel Aviv, Prievidza 2007, 2008
!
!-----------------------------------------------------------------------------

 use Dirac_HSMRCC  ! ... parameters for the external HSMRCC code
 use Parameters ! .. from hsmrcc_external/modules

implicit none
#include "priunit.h"

Logical, save :: set = .false.
Logical, intent(in) :: reset

CHARACTER:: PROMPT*1=' '
CHARACTER(len=7) :: WORD, KEYWORD=' '

! ... nice F90 version 
character(len=7) :: TABLE_TEMP, TABLE(10) =  (/'.PRINT ','.NACT  ','.NCORE ','.NOHSD ','.DRYRUN', & 
                                              '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX' /)

 if (set) then
  if (reset) set = .false.
  write(lupri,*) 'HSCCINP was called .... leaving at the beginning ... '
  Return
 endif
 if (reset) then
   set = .false.
 else
   set = .true.
 endif

 write(LUPRI,*) 'in HSCCINP, size of TABLE=',size(TABLE)

! REWIND (LUCMD,IOSTAT=IOS); print *,'input file rewind...' 

! Do While (KEYWORD /= '*HSCC  ')
!   READ (LUCMD,'(A7)',ERR=910, END=910) KEYWORD
!   print *,'read KEYWORD=',KEYWORD
! EndDo

 call FLSHFO(LUPRI)

!
!  Initializations...
!

Do 

  Read(LUCMD,'(A7)',ERR=910, END=930) WORD; PROMPT = WORD(1:1);  print *,'WORD=',WORD
  KEYWORD = WORD 

! Now decide !

  If (PROMPT == '*' ) Then
    print *,'WORD with "*" catched !', WORD
    Exit
  EndIf

  If ((PROMPT == '#') .or. (PROMPT == '!')) Then
     print *,'comment catched !', WORD
     Cycle ! ... continue the cycle..
  EndIf
  
  If (PROMPT == '.') Then
    ! ... keywords are beginning ...
    If (Any(TABLE(1:size(TABLE)) == WORD )) Then
     Select Case (WORD)
      Case ('.OPTION')
        Call PRTAB(size(TABLE), TABLE, KEYWORD//' input keywords',LUPRI)
      Case ('.NCORE')
        READ(LUCMD,*) NCMO_DIRAC
        write(LUPRI,*) ' read NCMO_DIRAC=',NCMO_DIRAC
      Case ('.NACT')
        READ(LUCMD,*) NACT_DIRAC
        write(LUPRI,*) ' read NACT_DIRAC=',NACT_DIRAC
      Case ('.NOHSDIR')
        nohsccdirac=.true. 
      Case ('.PRINT')
        READ(LUCMD,*) iprint_hscc
        write(LUPRI,*) ' read IPRINT_hscc=',iprint_hscc
      Case ('.DRYRUN')
        nohsccfiles=.true.
      Case Default
        Write(LUPRI,'(/,3A,/)') ' Keyword "',KEYWORD, &
        '" is not yet implemented in HSCCINP.'
      End Select
    EndIf

  Else
     WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,'" not recognized in HSCCINP.'
     CALL PRTAB(size(TABLE),TABLE,WORD//' input keywords',LUPRI)
     CALL QUIT('Illegal keyword in HSCCINP.')

  EndIf

EndDo

!  ... print out read parameters ... 

 WRITE (LUPRI,'(//)')
 CALL PRSYMB(LUPRI,'=',75,0)
 WRITE(LUPRI,'(1X,A)') 'HSCCINP: Set-up for the external standalone HSMRCC code'
 CALL PRSYMB(LUPRI,'=',75,0)
 WRITE(LUPRI,'(/1X,A,I3)') '* General print level   : ',iprint_hscc
! WRITE(LUPRI,'( 1X,A,I3)') '*                NELE   : ',NELE

 WRITE(LUPRI,'(1X,A,I3)')  '*             NCMO_DIRAC=',NCMO_DIRAC
 WRITE(LUPRI,'(1X,A,I3)')  '*             NACT_DIRAC=',NACT_DIRAC
! WRITE(LUPRI,'(1X,A,I3)')  '*            NFRZV_DIRAC=',NFRZV_DIRAC

 if (nohsccfiles) then
  WRITE(LUPRI,'(1X,A)') '* Do not prepare HSCC files (dry run). Take them from GAMESS.'
 else
  WRITE(LUPRI,'(1X,A)') '* Prepare HSCC files by DIRAC.'
 endif

 if (nohsccdirac) then
  WRITE(LUPRI,'(1X,A)') '* Only prepare HSCC files by DIRAC, do not continue with "afterDIRAC" HSCC.'
 endif

 call set_dirac_UOUT

 call set_dirac_hscc_interface_mode

 call check_dirac_hscc_interface_mode

 if (.not.get_dirac_hscc_interface_mode()) call quit('NOT in dirac_hscc_interface_mode !')

 call set_hscc_iprint(iprint_hscc) ! call hsmrcc_external routine 

 WRITE(LUPRI,*) 'HSCC get_iprint()=',get_iprint()
 WRITE(LUPRI,*) 'HSCC variable UOUT=',UOUT
 
 print *,' returning from HSCCINP ...'

 Return

910  CONTINUE
        CALL QUIT('Error reading LUCMD, no "*HSCC" input found !')
920  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSCC" input some error occured...')
930  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSCC" input end of file reached...')

end subroutine HSCCINP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_HSMRCC_Dirac_info
!-----------------------------------------------------------------------------
!
! Print-out basic info on the Dirac_HSMRCC interface run
!
! Called from:
!
! Written by Miro ILIAS, Tel Aviv 2007, Prievidza 2008
!
!-----------------------------------------------------------------------------
#include "priunit.h"

write(LUPRI,9001)

9001 format('                                                                 ',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '|        I N T E R F A C E                                      |',/, &
            '|                                                               |',/, &
            '| FOR THE STANDALONE PROGRAM FOR THE EFFECTIVE AND INTERMEDIATE |',/, &
            '| HAMILTONIAN HILBERT SPACE MULTI-REFERENCE COUPLED CLUSTER     |',/, &
            '| METHODS                                                       |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '| Written by Miro ILIAS, Shamasunar K.R.                        |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' )

end subroutine print_HSMRCC_Dirac_info

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7
function get_moltra_index(is_hole,i) result(moi)
!
! returns MOLTRA spinor index of HSMRCC NHMO/NPMO  index i
!
use io_moltra
use Parameters
use GMSMod

implicit none
#include "priunit.h"

integer, intent(in) :: i
logical, intent(in) :: is_hole
integer :: moi

!  write(LUPRI,*) 'NCMO:',NCMO
!  write(LUPRI,*) 'NVMO:',NVMO
!  write(LUPRI,*) 'NEMO:',NEMO

if (is_hole) then ! ... is hole (1-NHMO)
 ! write(LUPRI,*) 'hole index, i=',i
 !call print_gms_data_line(i)
 moi=MOLTRA_indx(i)
else ! .. is particle (index 1-NPMO)
 ! write(LUPRI,*) 'particle index, i=',i
 ! call print_gms_data_line(i+NCMO)
  moi=MOLTRA_indx(i+NCMO)
endif

!write(LUPRI,*) 'moi=',moi

end function get_moltra_index

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_moltra_integral(i,j,k,l) result(moint)
!-----------------------------------------------------------------------------
!
!  get MOLTRA integral (IJ|KL) from the MOLTRA integral list
!
!  if none found, simply returns zero
!
!------------------------------------------------------------------------------
!C     Real groups       : classes 1-4 are used
!C     Complex groups    : classes 1-8 are used
!C     Quaternion groups : classes 1-16 are used
!C
!C     Class 1  : Real part of integrals      (i    j    |k    l   )
!C     Class 2  : Real part of integrals      (i    j    |kbar lbar)
!C     Class 3  : Real part of integrals      (i    jbar |k    lbar)
!C     Class 4  : Real part of integrals      (i    jbar |kbar l   )
!C     Class 5  : Imaginary part of integrals (i    j    |k    l   )
!C     Class 6  : Imaginary part of integrals (i    j    |kbar lbar)
!C     Class 7  : Imaginary part of integrals (i    jbar |k    lbar)
!C     Class 8  : Imaginary part of integrals (i    jbar |kbar l   )
!C     Class 9  : Real part of integrals      (i    j    |kbar l   )
!C     Class 10 : Real part of integrals      (i    j    |k    lbar)
!C     Class 11 : Real part of integrals      (ibar j    |k    l   )
!C     Class 12 : Real part of integrals      (i    jbar |k    l   )
!C     Class 13 : Imaginary part of integrals (i    j    |kbar l   )
!C     Class 14 : Imaginary part of integrals (i    j    |k    lbar)
!C     Class 15 : Imaginary part of integrals (ibar j    |k    l   )
!C     Class 16 : Imaginary part of integrals (i    jbar |k    l   )
!-----------------------------------------------------------------------------
use io_moltra
use Parameters
use GMSMod

implicit none
#include "priunit.h"

integer, intent(in) :: i,j,k,l ! KU indexes
!integer, intent(in) :: i_kr,j_kr,k_kr,l_kr  ! ...
integer :: i_kr,j_kr,k_kr,l_kr
integer :: ix,jx,kx,lx,ipos_ijkl,ipos_jilk,ipos_klij,ipos_lkji,int_jilk
real*8 :: moint

moint=0.0d0

!write(LUPRI,'(4x,a,4i3)') 'get_moltra_integral: (i,j,k,l):',i,j,k,l

i_kr=+1;j_kr=+1;k_kr=+1;l_kr=+1

if (i>NSP/2) i_kr=-1
if (j>NSP/2) j_kr=-1
if (k>NSP/2) k_kr=-1
if (l>NSP/2) l_kr=-1


! TODO: print-out 
!write(LUPRI,'(4x,a,a,i3,a,i2,a,i3,a,i2,a,i3,a,i2,a,i3,a,i2,a)')  & 
!          'get_moltra_integral: (i,j,k,l):',                     & 
!          '(',i,'[',i_kr,']',j,'[',j_kr,']|',k,'[',k_kr,']',l,'[',l_kr,'])'

! ... take into account Kramers symmetries of integrals ...
if (i_kr.eq.+1.and.j_kr.eq.+1.and.k_kr.eq.+1.and.l_kr.eq.+1) then      ! (aa|aa)
 ix=i; jx=j; kx=k;     lx=l
else if (i_kr.eq.+1.and.j_kr.eq.+1.and.k_kr.eq.-1.and.l_kr.eq.-1) then ! (aa|bb)
 ix=i; jx=j; kx=k-NSP/2; lx=l-NSP/2
else if (i_kr.eq.-1.and.j_kr.eq.-1.and.k_kr.eq.+1.and.l_kr.eq.+1) then ! (bb|aa)
 ix=i-NSP/2; jx=j-NSP/2; kx=k; lx=l
else if (i_kr.eq.-1.and.j_kr.eq.-1.and.k_kr.eq.-1.and.l_kr.eq.-1) then ! (bb|bb)
 ix=i-NSP/2; jx=j-NSP/2; kx=k-NSP/2; lx=l-NSP/2
else
 return ! no suitable combination of 'bars' for nonzero integral
endif

! write(LUPRI,'(4x,a,4i3)') '...new indexes, i,j,k,l:',ix,jx,kx,lx 

!
! ... take into account general symmetry properties of (IJ|KL) 
!
!  (IJ|KL) = (JI|LK)=(KL|IJ)=(LK|JI) for Coulomb integrals
!

! (IJ|KL)

  ipos_ijkl=get_nonzero_integral_position(ix,jx,kx,lx)
  ipos_jilk=get_nonzero_integral_position(jx,ix,lx,kx)
  ipos_klij=get_nonzero_integral_position(kx,lx,ix,jx)
  ipos_lkji=get_nonzero_integral_position(lx,kx,jx,ix)
 
 if (ipos_ijkl.ne.-1) then
    moint=cbuf%p_rA2D(1,ipos_ijkl)
 else if (ipos_jilk.ne.-1) then
    moint=cbuf%p_rA2D(1,int_jilk)
 else if (ipos_klij.ne.-1) then
    moint=cbuf%p_rA2D(1,ipos_klij)
 else if (ipos_lkji.ne.-1) then
    moint=cbuf%p_rA2D(1,ipos_lkji)
 endif

 return

end function get_moltra_integral


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine make_dirac_syminf
!-----------------------------------------------------------------------------
! prepare file for the gamess-like input file
!
! $GMS/"syminf.f90"
!-----------------------------------------------------------------------------
!  end subroutine make_dirac_syminf
