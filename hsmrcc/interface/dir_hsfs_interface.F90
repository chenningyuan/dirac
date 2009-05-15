!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSFS_STUFF(WORK,LWORK)
!----------------------------------------------------------------------------------
!
!  Called from PAMTRA, if defined (MOD_HSFS)
!
!  Written by Miro Ilias, Tel Aviv/Prievidza 2009
!
!----------------------------------------------------------------------------------
use Memory
use io_moltra

#include "implicit.h"
#include "priunit.h"
#include "dcbpsi.h"

real*8, intent(in) :: WORK(LWORK)
!logical :: check_point = .false.
!logical :: check_point = .true.
!integer (kind=4) :: ipos1,ipos2,ipos3,ipos4
!integer :: ipos1,ipos2,ipos3,ipos4

!logical :: check_MOints = .false.
!logical :: check_MOints = .true.

if (.not.DOHSFS) then
 write(LUPRI,'(2x,a)') 'DOHSFS=.false., no HSFS data preparation'
 return
endif

call print_HSFS_Dirac_info

call StickMemory(WORK,LWORK)

HSFS_mode=.true.

call read_mrconee;  call print_mrconee_data

call scan_mdcint

call read_mdcint_hsfs

call UnstickMemory

end subroutine HSFS_STUFF


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_HSFS_Dirac_info
!-----------------------------------------------------------------------------
!
! Print-out basic info on the Dirac_HSMRCC interface run
!
! Called from:
!
! Written by Miro ILIAS, Tel Aviv 2009
!
!-----------------------------------------------------------------------------
#include "priunit.h"

write(LUPRI,9001)

9001 format('                                                                 ',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '|        I N T E R F A C E                                      |',/, &
            '|                                                               |',/, &
            '| FOR THE STANDALONE PROGRAMS FOR HIGHER-SECTOR COUPLED         |',/, &
            '| CLUSER CODES                                                  |',/, &
            '|                                                               |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '| Written by Miro ILIAS, Tel Aviv 2009                          |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' )

call FLSHFO(LUPRI)

end subroutine print_HSFS_Dirac_info

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSFSINP(WORD,reset)
!-----------------------------------------------------------------------------
!
!  Read parameters for the external HSFS programs, which are passed 
!  through the formatted input file, together with other DIRAC input
! 
!  Called from TRAINP if DOHFS=.true., i.e. at the beginning
!
!  Written by Miro ILIAS, Tel Aviv, 2009
!
!-----------------------------------------------------------------------------
  use DIRAC_HSFS
#include "implicit.h"
#include "priunit.h"

Logical, save :: set = .false.
Logical, intent(in) :: reset

CHARACTER:: PROMPT*1=' '
CHARACTER(len=7) :: WORD, KEYWORD=' '

! ... nice F90 version 
character(len=7) :: TABLE_TEMP, TABLE(10) =  (/'.PRINT ','.NACT  ','.NCORE ','.NOHSD ','.DRYRUN', & 
                                              '.CONV  ','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX' /)

 if (set) then
  if (reset) set = .false.
  write(lupri,*) 'HSFSINP was called .... leaving at the beginning ... '
  Return
 endif
 if (reset) then
   set = .false.
 else
   set = .true.
 endif


 write(LUPRI,*) 'in HSCCINP, size of TABLE=',size(TABLE)
 call FLSHFO(LUPRI)

!
!  Initializations...
!

Do 

  Read(LUCMD,'(A7)',ERR=910, END=930) WORD; PROMPT = WORD(1:1);  ! print *,'WORD=',WORD
  KEYWORD = WORD 

! Now decide !

  If (PROMPT == '*' ) Then
   ! print *,'WORD with "*" catched !', WORD
    Exit
  EndIf

  If ((PROMPT == '#') .or. (PROMPT == '!')) Then
   !  print *,'comment catched !', WORD
     Cycle ! ... continue the cycle..
  EndIf
  
  If (PROMPT == '.') Then
    ! ... keywords are beginning ...
    If (Any(TABLE(1:size(TABLE)) == WORD )) Then
     Select Case (WORD)
      Case ('.OPTION')
        Call PRTAB(size(TABLE), TABLE, KEYWORD//' input keywords',LUPRI)
!      Case ('.NCORE')
!        READ(LUCMD,*) NCMO_DIRAC
!        write(LUPRI,*) ' read NCMO_DIRAC=',NCMO_DIRAC
!      Case ('.NACT')
!        READ(LUCMD,*) NACT_DIRAC
!        write(LUPRI,*) ' read NACT_DIRAC=',NACT_DIRAC
!      Case ('.NOHSDIR')
!        nohsccdirac=.true. 
      Case ('.PRINT')
        READ(LUCMD,*) iprint_hsfs
        write(LUPRI,*) ' read IPRINT_hsfs=',iprint_hsfs
       Case ('.CONV')
         READ(LUCMD,*) iconv_hsfs
         write(LUPRI,*) ' read iconv_hsfs=',iconv_hsfs
!      Case ('.DRYRUN')
!        nohsccfiles=.true.
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
 WRITE(LUPRI,'(1X,A)') 'HSFSINP: Set-up for the external standalone HSFS codes'
 CALL PRSYMB(LUPRI,'=',75,0)
 WRITE(LUPRI,'(/1X,A,I3)') '* General print level   : ',iprint_hsfs
 WRITE(LUPRI,'(/1X,A,I3)') '* Convergence           : ',iconv_hsfs

! WRITE(LUPRI,'( 1X,A,I3)') '*                NELE   : ',NELE

! WRITE(LUPRI,'(1X,A,I3)')  '*             NCMO_DIRAC=',NCMO_DIRAC
! WRITE(LUPRI,'(1X,A,I3)')  '*             NACT_DIRAC=',NACT_DIRAC
! WRITE(LUPRI,'(1X,A,I3)')  '*            NFRZV_DIRAC=',NFRZV_DIRAC

! if (nohsccfiles) then
!  WRITE(LUPRI,'(1X,A)') '* Do not prepare HSCC files (dry run). You might take them from GAMESS.'
! else
!  WRITE(LUPRI,'(1X,A)') '* Prepare HSCC files by Dirac.'
! endif

! if (nohsccdirac) then
!  WRITE(LUPRI,'(1X,A)') '* Only prepare HSCC files by DIRAC, do not continue with "afterDIRAC" HSCC.'
! endif


 Return

910  CONTINUE
        CALL QUIT('Error reading LUCMD, no "*HSFS" input found !')
920  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSFS" input some error occured...')
930  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSFS" input end of file reached...')

end subroutine HSFSINP

