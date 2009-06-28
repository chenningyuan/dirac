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

!call read_mdcint_hsfs
!call read_mdcint_hsfs2
call read_mdcint_hsfs3

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
subroutine HSFSINP(WORD)
!-----------------------------------------------------------------------------
!
!  Read parameters for the external HSFS programs, which are passed 
!  through the formatted input file, together with other DIRAC input
! 
!  Called from PSIINP (dirrdn.F), entering WORD is label currently read
!
!  Written by Miro ILIAS, Tel Aviv, 2009
!
!-----------------------------------------------------------------------------
  use DIRAC_HSFS
  use high_sect_execs
#include "implicit.h"
#include "priunit.h"

CHARACTER(len=7), INTENT(INOUT) :: WORD

!Logical, save :: set = .false.
!Logical, intent(in) :: reset

Logical :: PROPER_POS

CHARACTER:: PROMPT*1=' ',PROMPT_M*1=' ' 
CHARACTER(len=7) :: METHOD, KEYWORD=' '

! ... nice F90 version 
character(len=7) :: TABLE_TEMP, TABLE(10) =  & 
  (/'.PRINT ','.NACT  ','.NCORE ','.NOHSD ','.DRYRUN', & 
    '.CONV  ','.METHOD','.XXXXXX','.XXXXXX','.XXXXXX' /)
!
!  Initializations...
!

!write(LUPRI,*) 'in HSFSINP, entering WORD=',WORDIN

! ... entering label MUST be the beginning of *HSFS input list
PROPER_POS = (WORD .EQ. '*HSFS')

IF (PROPER_POS) THEN

Do 

  Read(LUCMD,'(A7)',ERR=910, END=930) WORD; PROMPT = WORD(1:1);  ! print *,'WORD=',WORD
  KEYWORD = WORD 

  !write(LUPRI,*) 'HSFSINP: read WORD=',WORD,' PROMPT=',PROMPT

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
       Case ('.CONV')
         READ(LUCMD,*) iconv_hsfs
         write(LUPRI,*) ' read iconv_hsfs=',iconv_hsfs
       Case('.METHOD')
         write(LUPRI,*) ' going to read methdos !'
          call read_hsfs_methods(WORD)
!      Case ('.DRYRUN')
!        nohsccfiles=.true.
      Case Default
        Write(LUPRI,'(/,3A,/)') ' Keyword "',KEYWORD, &
        '" is not yet implemented in HSFSINP.'
     End Select
    EndIf
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
!-----------------------------------------------------------------------------
 use high_sect_execs
#include "implicit.h"
#include "priunit.h"
CHARACTER(len=7) :: WORDIN, METHOD, KEYWORD=' '
CHARACTER:: PROMPT*1=' ',PROMPT_M*1=' ' 

!write(LUPRI,*) 'HSFS methods input reading, WORDIN=',WORDIN

IF (.NOT.WORDIN == '.METHOD') THEN
   WRITE(LUPRI,*) 'WORDIN=',WORDIN
   CALL QUIT('read_hsfs_methods: worng entering word ! Must be .WORDIN ')
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
                call quit('number of methods exceeded !')
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

