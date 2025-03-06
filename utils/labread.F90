!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

!
      PROGRAM READLABL
!
! Written by Hans Jorgen Aa. Jensen, Nov. 1983
! CRAY version 10. Oct. 1986
! Alliant version 30. Dec. 1987
! Linux/IBM, +updates, Joern Thyssen, Nov. 1999
!
!#include "implicit.h"
#include "../include/implicit.h"
!TROND      DOUBLE PRECISION TIME, SECOND
      CHARACTER*8 B(4), STARS
      PARAMETER (STARS = '********')
      PARAMETER (LUINP = 1, LUPRI = 6)
      INTEGER   NARG,IARG,JARG
      CHARACTER*80 ARGC
!
      NARG = command_argument_count()
      IF (NARG .EQ. 0) THEN
         CALL get_command_argument(0,ARGC)
         JARG = LEN_TRIM(ARGC)
         WRITE(6,'(/3A/)') 'Usage: ',ARGC(1:JARG),' file1 file2 ...'
         STOP ' '
      END IF
      DO IARG = 1,NARG
         CALL get_command_argument(IARG,ARGC)
         JARG = LEN_TRIM(ARGC)
         WRITE (LUPRI, '(/2A/)')                                        &
       ' MOLECULE labels found on file :',ARGC(1:JARG)
!
         OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
         REWIND LUINP
         IREC = 0
         IERR = 0
    1    READ (LUINP,END=3,ERR=2) B
            IREC = IREC + 1
         IF (B(1) .NE. STARS) GO TO 1
            WRITE (LUPRI, '(5X,I5,3X,4(2X,A8))')  IREC, B
         GO TO 1
!
    2    CONTINUE
         IREC = IREC + 1
         IERR = IERR + 1
         WRITE (LUPRI, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
         REWIND LUINP
         DO 102 I = 1,IREC
            READ (LUINP) J
  102    CONTINUE
         IF (IERR .LE. 2) GO TO 1
  202    CONTINUE
            READ (LUINP,END=3) J
            IREC = IREC + 1
         GO TO 202
!
    3    CONTINUE
!      TIME = SECOND() - TIME
         WRITE (LUPRI,'(/I10,A)') IREC,                                 &
       ' records read before EOF on file.'
!        WRITE (LUPRI,'(F10.2,A)') TIME,' CPU seconds used.'
         CLOSE(LUINP)
      END DO
      STOP ' '
!
      END
