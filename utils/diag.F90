!************************************************************************ 
!
! Standalone program for comprehensive testing of all diagonalization
! routines used by DIRAC. 
!
! These are Householder, RSJACOBI, QJACOBI and
! DSYEVR.
!
! The program either reads the matrix file, or generates a random matrix.
!
!
! Written by Miro Ilias, 07,08-2014 Banska Bystrica, Slovakia
! with the help of P.S. Bagus and T.Saue
!
!************************************************************************ 

module time
  implicit none
  character*10 start_date, end_date ! format of date: "dd_mmm_yy "
  character*8  start_time, end_time ! format of time: "hh:mm:ss"
end module time

Program Test_DIRAC_Diagonalization_Routines
 use diagmod
#include "priunit.h"
  write(LUPRI,'(/,2x,a,/)') "**** Welcome to the testing program of DIRAC diagonalization routines ! ****"
  call DBLGRP !> initialize include file variable for printing matrices (PRQMAT)
  call read_input_file
  call read_matrix
  if (do_rs)      call test_dirac_rs
  if (do_dsyevr)  call test_lapack_dsyevr
  if (do_dsyevr_ulf)  call test_lapack_dsyevr_ulf
  if (do_qjaco)   call QJACOBITEST
  if (do_householder_psb) call test_PSB_householder
End Program Test_DIRAC_Diagonalization_Routines

!--------------------------------------------------------------
! Reads the **DIAG section of the text input file "DIAGONALIZE_TESTS.INP".
! This small self-contained reader is used instead of the DIRAC menu reader
! (read_menu_input), which dispatches to the input readers of all DIRAC
! modules and would thereby link nearly all of DIRAC into this program.
! Keywords are case-insensitive; .TITLE, .PRINT and .MTXFIL take their
! value from the next line.
!--------------------------------------------------------------
subroutine read_input_file
  use diagmod, only: title_text, print_level, matrix_file_name, do_dsyevr, do_dsyevr_ulf, &
                     do_rs, do_rsjaco, do_qjaco, do_householder_psb, do_eigv_check
  implicit none
  character*21 :: inpfile="DIAGONALIZE_TESTS.INP"
  logical :: input_found
  character(len=7) :: word, section
  integer :: i
#include "priunit.h"
    inquire(file=inpfile, exist=input_found)
    if (.not.input_found) then
      write(LUPRI,"(2x,a,a)") "read_input_file: the input file not found ! input_file=",inpfile
      call quit("read_input_file: the input file not found !")
    endif
    open(LUSTDIN, file=inpfile, status='old', action='read')
    section = ' '
    do
      read(LUSTDIN, '(a7)', end=1) word
      do i = 1, len(word)   ! to upper case
        if (word(i:i) >= 'a' .and. word(i:i) <= 'z') word(i:i) = achar(iachar(word(i:i)) - 32)
      end do
      if (word == ' ' .or. word(1:1) == '!' .or. word(1:1) == '#') cycle
      if (word(1:1) == '*') then
        if (word == '*END OF' .or. word == '**END O') exit
        section = word
        cycle
      end if
      if (section /= '**DIAG') cycle
      select case (word)
        case ('.TITLE ')
          read(LUSTDIN, '(a)', err=2, end=2) title_text
        case ('.PRINT ')
          read(LUSTDIN, *, err=2, end=2) print_level
        case ('.MTXFIL')
          read(LUSTDIN, '(a)', err=2, end=2) matrix_file_name
        case ('.DSYEVR')
          do_dsyevr = .true.
        case ('.ULF   ')
          do_dsyevr_ulf = .true.   ! Ulf's simple routine
        case ('.RS    ')
          do_rs = .true.
        case ('.RSJACO')
          do_rsjaco = .true.
        case ('.QJACO ')
          do_qjaco = .true.
        case ('.PAUL  ')
          do_householder_psb = .true.
        case ('.EIGVTS')
          do_eigv_check = .true.
        case default
          write(LUPRI,"(2x,a,a)") "read_input_file: unknown keyword in **DIAG: ",word
          call quit("read_input_file: unknown keyword in **DIAG")
      end select
    end do
1   close(LUSTDIN)
    return
2   write(LUPRI,"(2x,a,a)") "read_input_file: could not read the value of keyword ",word
    call quit("read_input_file: could not read keyword value")
end subroutine read_input_file
