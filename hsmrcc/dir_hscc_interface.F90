!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSCCINP(WORD,reset)
!-----------------------------------------------------------------------------
!
! Read parameters for the external HSMRCC program, which are passed through
! the file.
! 
! Called from TRAINP if DOHSCC=.true.
!
! Written by Miro ILIAS, Tel Aviv, Prievidza, 2007, 2008
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
      ! print  *,' read NELE=',NELE
      Case ('.NACT')
        READ(LUCMD,*) NACT
      Case ('.PRINT')
        READ(LUCMD,*) iprint_hscc
      !  print  *,' read IPRINT_hscc=',iprint_hscc
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

  subroutine print_HSMRCC_Dirac_info
!------------------------------------------------------------
!
!
!------------------------------------------------------------
#include "priunit.h"

   write(LUPRI,'(/2x,a/)') 'interface for eff. and intermed.Ham. HSMRCC method'

!  write(LUPRI,9001)

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
!----------------------------------------------------------------------------------
!
!  Called from PAMTRA, if defined (MOD_HSCC)
!
!  Written by Miro Ilias, Tel Aviv/Prievidza 2007,2008
!
!----------------------------------------------------------------------------------
use Memory
use Dirac_HSMRCC

#include "implicit.h"
#include "priunit.h"
#include "dcbpsi.h"

!integer,intent(in) :: LWORK
real*8, intent(in) :: WORK(LWORK)
logical :: check_point = .false.
!logical :: check_point = .true.
integer (kind=4) :: ipos1,ipos2,ipos3,ipos4

!logical :: check_MOints = .false.
logical :: check_MOints = .true.

Call QEnter('HSCC_STUFF')

if (.not.DOHSCC) then
 write(LUPRI,'(2x,a)') 'DOHSCC=.false., no HSCC data preparation'
 Call QExit('HSCC_STUFF')
 return
endif

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

if (check_MOints) then
!call read_moltra_integrals
 call read_mrconee; call print_mrconee_data
 call read_mdcint ! read MDCINT file
endif

call UnstickMemory

Call QExit('HSCC_STUFF')
end subroutine HSCC_STUFF

! ----------------------------------


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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_mdcint
!----------------------------------------------------
!
!    ... read MDCINT - testing purposes
!
!----------------------------------------------------
  ! use Memory
   use Vectors
   use io_moltra
#include "implicit.h"
#include "priunit.h"

type(Vect) :: indk_vect, indl_vect, indmo_vect
type(Array2Dim) :: cbuf_vect

logical :: is_mdc

inquire(file="MDCINT",EXIST=is_mdc); if (.not.is_mdc) call QUIT('No MDCINT file!')

!open(MDINT,file="MDCINT",form='UNFORMATTED',status='old')
open(MDINT,file="MDCINT",form='UNFORMATTED')
rewind(MDINT)

write(LUPRI,'(//2x,a,i3)') 'MDCINT open and rewinded...lu=',MDINT

READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,(KR(I),KR(-I), I=1,NKR )

!READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR
!READ (MDINT,ERR=10000,END=10000) (  KR(I),KR(-I), I=1,NKR )

if (NKR > NSP2) then
 call QUIT('(NKR > NSP2 ! Error quit when reading MDCINT !')
endif

 write(LUPRI,*) 'Info from MDCINT file'

 write(LUPRI,'(A,$)') 'Integrals have been calculated on the'
 write(LUPRI,'(A,A12,1X,A,A)') '      ',DATEX,' at ',TIMEX

!do i = 1, NKR
!  write(LUPRI,'(2x,a,i3,a,i3,a,i3)') 'i=',i,' kr(-i)=',KR(-i),' kr(i)=',KR(i)
!enddo

        write(LUPRI,*) 'NKR is ',NKR
        write(LUPRI,*) 'KR array:'
        do I=1,NKR,1
          write(LUPRI,'(2I8)') KR(-I),KR(I)
        end do
!
! ... read the rest ... integrals and indexes  ...
! ... allocate memmory for integrals/indexes   ...
!

is_complex = .false.

ilenght = NKR*NKR ! note: "ccdriv.F":1923 is using RCW !!!

indk_vect = MakeVector('INTE',ilenght,'MDCINT "k" indexes')

indl_vect = MakeVector('INTE',ilenght,'MDCINT "l" indexes')

cbuf_vect = Make2DimArray('REAL',ilenght,1,'MDCINT two-el. integrals')

indmo_vect = MakeVector('INTE',NKR,'for MDCINT - indmo array')

Call ZeroObject(indmo_vect)
Call ZeroObject(indk_vect)
Call ZeroObject(indl_vect)
Call ZeroObject(cbuf_vect)

!do i=1,NKR ! fill...
!  indmo_vect%p_iA1D(i) = i
!enddo

  MOINT = 0

  rewind(MDINT)
  read(MDINT)

 1 continue
   READ(MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,                       &
   ( indk_vect%p_iA1D(INZ), indl_vect%p_iA1D(INZ), INZ=1,NONZR,1),      & 
   ( cbuf_vect%p_rA2D(1,INZ), INZ=1,NONZR,1) 

  ! call print_vect_attrib(indk_vect)
  ! call print_vect_attrib(indl_vect)

    write(LUPRI,'(2x,a,3i4)') 'read IKR,JKR,NONZR:',IKR,JKR,NONZR

    if (NONZR > ilenght) then
     call Quit('HSCC: NONZR > ilenght !!!')
    endif

 do i=1,NONZR
   write(LUPRI,*) indk_vect%p_iA1D(i),indl_vect%p_iA1D(i),cbuf_vect%p_rA2D(1,i)
 enddo
   
    MOINT = MOINT + NONZR

    IF (IKR.EQ.0) THEN
         IF (JKR.EQ.0) THEN
            NONZR = -1
         ELSE
            READ (MDINT) NEXTFIL
            CLOSE (MDINT,STATUS='KEEP')
            OPEN(MDINT,FILE=NEXTFIL,FORM='UNFORMATTED')
            GOTO 1
         ENDIF
        write(LUPRI,*) 'HSCC: End of file MDCINT.'
        GoTo 1
    ENDIF

 if (NONZR == -1) then
  ! call Quit('NONZR == -1!')
 !  RETURN
   GoTo 2
 endif

! do i=1,NONZR
!   write(LUPRI,*) indk_vect%p_iA1D(i),indl_vect%p_iA1D(i),cbuf_vect%p_rA2D(1,i)
! enddo

!  GoTo 1

 !write(LUPRI,*) 'all nonzero integrals from MDCINT together with indexes read ... '

 i = indmo_vect%p_iA1D(KR(IKR))
 j = indmo_vect%p_iA1D(KR(JKR))

 itr = indmo_vect%p_iA1D(KR(-IKR))
 jtr = indmo_vect%p_iA1D(KR(-JKR))

 SignIJ = SIGN(1,IKR) * SIGN(1,JKR)

 do INZ = 1, NONZR
  kkr = indk_vect%p_iA1D(INZ)
  lkr = indl_vect%p_iA1D(INZ)
  
  k = indmo_vect%p_iA1D(KR(kkr))
  l = indmo_vect%p_iA1D(KR(lkr))

  ktr = indmo_vect%p_iA1D(KR(-kkr))
  ltr = indmo_vect%p_iA1D(KR(-lkr))

  SignKL = SIGN(1,kkr) * SIGN(1,lkr)

! --> Original integral
    IF (i.NE.0.AND.j.NE.0.AND.k.NE.0.AND.l.NE.0) THEN
      write(LUPRI,'(4i3,1x,a,f22.14)')  i,j,k,l,':',cbuf_vect%p_rA2D(1,INZ)
    ENDIF

 enddo

 2 continue
 CLOSE (MDINT,STATUS='KEEP')

 return

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

