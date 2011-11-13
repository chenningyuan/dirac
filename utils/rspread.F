C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
      PROGRAM RSPREAD
C
C Written by Trond Saue July 4 2001
C based on program LABREAD
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (LUINP = 1)
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2,LBUF*16
      INTEGER   NARG,IARG,JARG
      CHARACTER*80 ARGC
C
      NARG = IARGC()
      IF (NARG .EQ. 0) THEN
         CALL GETARG(0,ARGC)
         JARG = LNBLNK(ARGC)
         WRITE(6,'(/3A/)') 'Usage: ',ARGC(1:JARG),' file1 file2 ...'
         STOP ' '
      END IF
      DO IARG = 1,NARG
         CALL GETARG(IARG,ARGC)
         JARG = LNBLNK(ARGC)
         WRITE (LUPRI, '(/2A/)') 
     &      ' Solution vectors found on file : ',ARGC(1:JARG)
C
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
         READ (LUINP,END=3,ERR=2) LAB1,LAB2,TYP1,TYP2,
     &          FREQ1,FREQ2,JSYMOP1,JSYMOP2,
     &          JTIMOP1,JTIMOP2,RNORM,LEN,INTFLG,
     &          ERGRSP,NBAS,NORB
         IVEC = IVEC + 1
         IREC = IREC + 1
         IF(LAB1.NE.LBUF) THEN
           WRITE(LUPRI,'(A,A16,2(A,I3))') 
     &      '** Solution vector : ',LAB1,
     &           ' Irrep: ',JSYMOP1,' Trev: ',JTIMOP1
           LBUF = LAB1
         ENDIF
         WRITE(LUPRI,'(I3,A,A2,A,F12.6,A,E10.2,A,I10)')
     &     IVEC, ' Type : ',TYP1,' Freq.: ',FREQ1,
     &     ' Rnorm: ',RNORM,' Length:', LEN
         READ(LUINP,END=3,ERR=2)
         IREC = IREC + 1
         GO TO 1
C
    2    CONTINUE
         IREC = IREC + 1
         IERR = IERR + 1
         WRITE (LUPRI, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
         GOTO 10
C
    3    CONTINUE
         WRITE (LUPRI,'(/I10,A)') IREC,
     *      ' records read before EOF on file.'
         GOTO 10
 10      CONTINUE
         CLOSE(LUINP)
      ENDDO
C
      END
