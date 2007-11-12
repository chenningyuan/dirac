!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSCCINP(WORD,reset)
!-----------------------------------------------------------------------------
!
! Read parameters for the external HSMRCC program, which are passed through
! the file.
! 
! Called from TRAINP
!
!
!
!-----------------------------------------------------------------------------

 use Dirac_HSMRCC  ! ... contains all parameters for the external HSMRCC code

#include "implicit.h"
#include "priunit.h"

Logical, save :: set = .false.
Logical, intent(in) :: reset

CHARACTER:: PROMPT*1=''
CHARACTER(len=7) :: WORD, KEYWORD=''

! ... nice F90 version 
!character(len=7) :: TABLE_TEMP, TABLE(5) =  (/'.PRINT ','.NELEC ','.NACT  ','.XXXXXX','.XXXXXX' /)
character(len=7) :: TABLE_TEMP, TABLE(5) =  (/'.PRINT ','.NELEC ','.XXXXXX','.XXXXXX','.XXXXXX' /)

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

 Call QENTER('HSCCINP')

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
      Case ('.NELEC')
        READ(LUCMD,*) NELE
        print  *,' read NELE=',NELE
      Case ('.NACT')
        READ(LUCMD,*) NELE
        print  *,' read NELE=',NELE
      Case ('.PRINT')
        READ(LUCMD,*) iprint_hscc
        print  *,' read IPRINT_hscc=',iprint_hscc
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

! Call QUIT('in HSCCINP !')

 print *,' returning from HSCCINP ...'

 Call QEXIT('HSCCINP')
 Return

910  CONTINUE
        CALL QUIT('Error reading LUCMD, no "*HSCC" input found !')
920  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSCC" input some error occured...')
930  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSCC" input end of file reached...')

end subroutine HSCCINP


subroutine HSMRCC_Dirac_input


end subroutine HSMRCC_Dirac_input

!
!
!
  subroutine print_HSMRCC_Dirac_info
#include "priunit.h"

  write(LUPRI,9001)

9001 format('                                                                 ',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|        I N T E R F A C E                                      |',/, &
            '|      PROGRAM FOR EFFECTIVE AND INTERMEDIATE HAMILTONIAN       |',/, &
            '|     HILBERT SPACE MULTI-REFERENCE COUPLED CLUSTER METHODS     |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '|           Shamasundar. K. R. and Sourav Pal                   |',/, &
            '|           Theoretical Chemistry Group                         |',/, &
            '|           Physical Chemistry Division                         |',/, &
            '|           National Chemical Laboratory (NCL)                  |',/, &
            '|           Pashon Road, Pune - 411 008                         |',/, &
            '|           INDIA                                               |',/, &
            '|                                                               |',/, &
            '|                       AND                                     |',/, &
            '|                                                               |',/, &
            '|           Ephraim Eliav and Uzi Kaldor                        |',/, &
            '|           School of Chemistry                                 |',/, &
            '|           Tel Aviv University                                 |',/, &
            '|           Tel Aviv, 69978                                     |',/, &
            '|           ISRAEL                                              |',/, &
            '|                                                               |',/, &
            '|    Note : This version is ported for opteron machines at      |',/, &
            '|           University of Waterloo and Tel Aviv University.     |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                INFORMATION ON CURRENT RUN                     |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' )


  end subroutine print_HSMRCC_Dirac_info


!
! prepare file for the gamess-like input file
!
! $GMS/"syminf.f90"
  subroutine make_gamess_syminf

  end subroutine make_gamess_syminf

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSCC_STUFF(WORK,LWORK)

use Memory
use Dirac_HSMRCC

#include "implicit.h"
#include "priunit.h"
!integer,intent(in) :: LWORK
real*8, intent(in) :: WORK(LWORK)
logical :: check_point = .false.
!logical :: check_point = .true.
integer (kind=4) :: ipos1,ipos2,ipos3,ipos4

write(LUPRI,*)
write(LUPRI,*) 'in HSCC_STUFF !! LWORK=',LWORK

  call print_HSMRCC_Dirac_info

  call StickMemory(WORK,LWORK)

if (check_point) then
 N=4;NNZ=2; ISIZE=N*N*NNZ
 call mymemget('INTE',ipos1,ISIZE)
 call mymemget('REAL',ipos2,ISIZE)
 call mymemget('LOGI',ipos3,ISIZE)
 call mymemget('INTE',ipos4,ISIZE)
! call assign_pointers77( PW(ipos1), PW(ipos2), PW(ipos3), PW(ipos4),N,N,NNZ)
! call check_pointers
endif

!call read_moltra_integrals
call read_mrconee; call print_mrconee_data
call read_mdcint ! read MDCINT file

call UnstickMemory

end subroutine HSCC_STUFF

! ----------------------------------

subroutine read_moltra_integrals
 use io_moltra
#include "priunit.h"
!logical is_mdcint, is_mrconee
open(MDINT,file="MDCINT",form='UNFORMATTED',status='old')
write(LUPRI,*) 'MDCINT open...'
! READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,  & 
!                                       (KR(I),KR(-I),I=1,NKR)
close(MDINT)
end subroutine read_moltra_integrals

subroutine read_mrconee

 use io_moltra

#include "implicit.h"
#include "priunit.h"

 logical :: is_mrc

     inquire(FILE='MRCONEE',EXIST=is_mrc)

     if (.not.is_mrc) call QUIT('mrconee file not present !')

       OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED',status = 'old')
       READ (MRCONEE,END=10,ERR=10) NSP,BREIT,ECORE, &
                                    NFSYM,NZ,SPFR,NORBT
       READ (MRCONEE,END=10,ERR=10) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP)  & 
                      ,(NELEC(IRP),IRP=1,NSYMRP)
       GOTO 11

!C---- Keep OLD (MOLFDIR) format here for compatibility
   10     WRITE (LUPRI,*) "Missing information on MRCONEE, I assume"//  & 
          " this is a MOLFDIR file and try to continue.."
!         INTERFACE='MOLFDIR'
!         REWIND(MRCONEE)
!         READ (MRCONEE) NSP,BREIT,ECORE
!         READ (MRCONEE) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP) 
!         DO IRP = 1, NSYMRP
!            NELEC(IRP) = 0
!         ENDDO

   11  CONTINUE
       READ (MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)
       READ (MRCONEE) ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP)
!C      IREPSP is to be used as 2-dimensional array, read it like that
       READ (MRCONEE) (IREPSP(ISP),IREPSP(NSP+ISP),   & 
                       EPS1(ISP),ISP=1,NSP)
       CLOSE (MRCONEE,STATUS='KEEP')

end subroutine read_mrconee

subroutine print_mrconee_data
 use io_moltra
#include "implicit.h"
#include "priunit.h"

write(LUPRI,*) 'NREP=',NREP,' NSP=',NSP
write(LUPRI,*) 'REPNA: ',(REPNA(IRP),IRP=1,2*NREP)
write (LUPRI,*) 'EPS1: ',(EPS1(ISP),ISP=1,NSP)

end subroutine print_mrconee_data

!
!  ... read MDCINT 
!
subroutine  read_mdcint
   use io_moltra
#include "implicit.h"
#include "priunit.h"

logical :: is_mdc

inquire(file="MDCINT",EXIST=is_mdc); if (.not.is_mdc) call QUIT('No MDCINT file!')

open(MDINT,file="MDCINT",form='UNFORMATTED',status='old')
write(LUPRI,*) 'MDCINT open...'
!READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,( KR(I),KR(-I), I=1,NKR )

CLOSE (MDINT,STATUS='KEEP')

 10000 CALL QUIT('ERROR READING HEADER MDCINT')


end subroutine  read_mdcint

