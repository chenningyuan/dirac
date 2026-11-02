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

C
C AUTHOR : Juan J. Aucar, 2023.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE NQCC_INPUT(WORD,WORK,LWORK)
      use priunit_mod, only: LUCMD, LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
!      use dcbprp_mod, only: NQCCCART
#include "implicit.h"
      INTEGER :: I
#include "dcbprp.h"
      PARAMETER (NTABLE = 1)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
      DATA TABLE /'.CARTES'/
C
C      
      WRITE (LUPRI,'(/,A)') ' * NQCC input section:'

C     Process input
C     =========================
C
      NEWDEF = (WORD .EQ. '*NQCC  ')
      ICHANG = 0
      INPERR = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in NQCCINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in NQCCINP.')
C
    1          CONTINUE
               NQCCcart = .TRUE.
               WRITE(LUPRI,'(3X,A)') 'The NQCC will be calculated'//
     &                               ' in the cartesian axes system.'
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in NQCCINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in NQCCINP.')
            END IF
      END IF
  300 CONTINUE
  999 CONTINUE           

      END
