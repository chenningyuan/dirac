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

      PROGRAM RSPREAD
!
! Written by Trond Saue July 4 2001
! based on program LABREAD
!
#include "implicit.h"
#include "priunit.h"
      PARAMETER (LUINP = 1)
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2,LBUF*16
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
       ' Solution vectors found on file : ',ARGC(1:JARG)
!
         OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
         REWIND LUINP
         IREC = 0
         IERR = 0
         IVEC = 0
         LBUF = '                '
    1    CONTINUE
         READ (LUINP,END=3,ERR=2) LAB1
         IF(LAB1.EQ.'END_OF_THIS_FILE') GOTO 10
         BACKSPACE LUINP
         READ (LUINP,END=3,ERR=2) LAB1,LAB2,TYP1,TYP2,                  &
           FREQ1,FREQ2,JSYMOP1,JSYMOP2,                                 &
           JTIMOP1,JTIMOP2,RNORM,LEN,INTFLG,                            &
           ERGRSP,NBAS,NORB
         IVEC = IVEC + 1
         IREC = IREC + 1
         IF(LAB1.NE.LBUF) THEN
           WRITE(LUPRI,'(A,A16,2(A,I3))')                               &
       '** Solution vector : ',LAB1,                                    &
            ' Irrep: ',JSYMOP1,' Trev: ',JTIMOP1
           LBUF = LAB1
         ENDIF
         WRITE(LUPRI,'(I3,A,A2,A,F12.6,A,E10.2,A,I10)')                 &
      IVEC, ' Type : ',TYP1,' Freq.: ',FREQ1,                           &
      ' Rnorm: ',RNORM,' Length:', LEN
         READ(LUINP,END=3,ERR=2)
         READ(LUINP,END=3,ERR=2)
         IREC = IREC + 1
         GO TO 1
!
    2    CONTINUE
         IREC = IREC + 1
         IERR = IERR + 1
         WRITE (LUPRI, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
         GOTO 10
!
    3    CONTINUE
         WRITE (LUPRI,'(/I10,A)') IREC,                                 &
       ' records read before EOF on file.'
         GOTO 10
 10      CONTINUE
         CLOSE(LUINP)
      ENDDO
!
      END
