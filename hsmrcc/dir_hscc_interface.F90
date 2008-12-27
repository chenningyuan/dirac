!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSCC_STUFF(WORK,LWORK)
!----------------------------------------------------------------------------------
!
!  Called from PAMTRA, if defined (MOD_HSCC)
!
!  Written by Miro Ilias, Tel Aviv/Prievidza 2007,2008
!
!----------------------------------------------------------------------------------
use Memory
use Dirac_HSMRCC
use io_moltra

#include "implicit.h"
#include "priunit.h"
#include "dcbpsi.h"

!integer,intent(in) :: LWORK
real*8, intent(in) :: WORK(LWORK)
logical :: check_point = .false.
!logical :: check_point = .true.
!integer (kind=4) :: ipos1,ipos2,ipos3,ipos4
integer :: ipos1,ipos2,ipos3,ipos4

!logical :: check_MOints = .false.
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
 call read_mrconee; call print_mrconee_data
 call read_mdcint ! read MDCINT file
!endif

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

#include "implicit.h"
#include "priunit.h"

Logical, save :: set = .false.
Logical, intent(in) :: reset

CHARACTER:: PROMPT*1=''
CHARACTER(len=7) :: WORD, KEYWORD=''

! ... nice F90 version 
!character(len=7) :: TABLE_TEMP, TABLE(5) =  (/'.PRINT ','.NELEC ','.NACT  ','.XXXXXX','.XXXXXX' /)
!character(len=7) :: TABLE_TEMP, TABLE(5) =  (/'.PRINT ','.NELEC ','.XXXXXX','.XXXXXX','.XXXXXX' /)
character(len=7) :: TABLE_TEMP, TABLE(5) =  (/'.PRINT ','.NACT  ','.NCORE ','.DRYRUN','.XXXXXX' /)

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
        READ(LUCMD,*) NCMO
        write(LUPRI,*) ' read NCMO=',NCMO
      Case ('.NACT')
        READ(LUCMD,*) NACT
        write(LUPRI,*) ' read NACT=',NACT
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
 WRITE(LUPRI,'( 1X,A,I3)') '*                NELE   : ',NELE

 if (nohsccfiles) then
  WRITE(LUPRI,'(1X,A)') '* Do not prepare HSCC files (dry run).'
 else
  WRITE(LUPRI,'(1X,A)') '* Prepare HSCC files by Dirac.'
 endif

! Call QUIT('in HSCCINP !')

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

!  write(LUPRI,'(/2x,a/)') 'interface for eff. and intermed.Ham. HSMRCC method'

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine make_dirac_syminf
!-----------------------------------------------------------------------------
! prepare file for the gamess-like input file
!
! $GMS/"syminf.f90"
!-----------------------------------------------------------------------------
!  end subroutine make_dirac_syminf

