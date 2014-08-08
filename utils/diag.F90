!************************************************************************ 
! Standalone program for comprehensive testing of all diagonalization
! routines used by DIRAC. 
!
! These are Householder, RSJACOBI, QJACOBI and
! DSYEVR.
!
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
 
   write(LUPRI,'(2x,a)') "Welcome to the testing program for DIRAC's diagonalization routines !"

  call DBLGRP

  call read_input_file

  call read_matrix

  if (do_rs)      call test_dirac_rs
  if (do_dsyevr)  call test_lapack_dsyevr
  if (do_dsyevr_ulf)  call test_lapack_dsyevr_ulf
  if (do_qjaco)   call QJACOBITEST
  

End Program Test_DIRAC_Diagonalization_Routines
