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

      PROGRAM DIAGRID
!***********************************************************************
!
!     Utility program for diatomics:
!       generates even-spaced grid within classical turning points
!
!***********************************************************************
      use codata
#include "implicit.h"
#include "priunit.h"
      INTEGER ILOGG
      CHARACTER PROJECT*40,REPLY*1
      ILOGG = 1
!
!     Heading
!
      CALL TITLER('DIAGRID for diatomics : Written by T. Saue ','*',110)
      WRITE(LUPRI,'(A)') 'Give project name (A40)'
      READ(LUSTDIN,*) PROJECT
      i=LNBLNK(PROJECT)
      OPEN(ILOGG,FILE=PROJECT(1:i)//'.diagrid',STATUS='UNKNOWN',          &
      FORM='FORMATTED',ACCESS='SEQUENTIAL')
      WRITE(ILOGG,*) '* Input data:'
!
!     Equilibrium bond length
!     =======================
!
      WRITE(LUPRI,'(A)') 'Give equilibrium bond length in Angstroms: '
      READ(LUSTDIN,*) REQ
      WRITE(LUPRI,'(A)') 'Give harmonic frequency in cm-1'
      READ(LUSTDIN,*) WEQ
      WRITE(LUPRI,'(A)') 'Select vibrational quantum number'
      READ(LUSTDIN,*) NU
      WRITE(ILOGG,*) ' - Equilibrium bond length (Angstrom)   ',REQ
      WRITE(ILOGG,*) ' - Harmonic frequency (cm-1)            ',WEQ
      WRITE(ILOGG,*) ' - Selected vibrational quantum number  ',NU
!
!     Get masses
!     ==========
!
 10   CONTINUE
      WRITE(LUPRI,'(A)')                                                &
   'Select one of the following:',                                      &
   '  1. Select masses of the most abundant isotopes.',                 &
   '  2. Employ user-defined atomic masses.',                           &
   '  3. Normal modes: Give reduced mass.'
      READ(LUSTDIN,*) I
      IF(I.EQ.1) THEN
        WRITE(LUPRI,'(A)') '* Give charge of atom A:'
        READ(LUSTDIN,*) IZ
        CALL GETMSS(IZ,AM,ABUND,INFO)
        IF(INFO.EQ.0) THEN
          WRITE(LUPRI,'(A,F12.4)')                                      &
     '* Mass     :', AM,                                                &
     '* Abundance:', ABUND
        ELSE
          WRITE(LUPRI,'(A)') 'Routine failed. Give mass in Daltons:'
          READ(LUSTDIN,*) AM
        ENDIF
        WRITE(LUPRI,'(A)') '* Give charge of atom B :'
        READ(LUSTDIN,*) IZ
        CALL GETMSS(IZ,BM,ABUND,INFO)
        IF(INFO.EQ.0) THEN
          WRITE(LUPRI,'(A,F12.4)')                                      &
     '* Mass     :', BM,                                                &
     '* Abundance:', ABUND
        ELSE
          WRITE(LUPRI,'(A)') 'Routine failed. Give mass in Daltons:'
          READ(LUSTDIN,*) BM
        ENDIF
      ELSEIF(I.EQ.2) THEN
        WRITE(LUPRI,'(A)') '* Give mass of atom A in Daltons:'
        READ(LUSTDIN,*) AM
        WRITE(LUPRI,'(A)') '* Give mass of atom B in Daltons:'
        READ(LUSTDIN,*) BM
      ELSEIF(I.EQ.3) THEN
        WRITE(LUPRI,'(A)') '* Give reduced mass in Daltons:'
        READ(LUSTDIN,*) UM
        GOTO 50
      ELSE
        WRITE(LUPRI,'(A)') ' You stupid fart ! Try again !'
        GOTO 10
      ENDIF
      UM = D1/(D1/AM + D1/BM)
      WRITE(ILOGG,'(A)') '* MASSES:'
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom A      :',AM
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom B      :',BM
 50   CONTINUE
      WRITE(ILOGG,*) ' - Reduced mass:',UM
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      UM = XFAMU*UM
      CLOSE(ILOGG,STATUS='KEEP')
      END
