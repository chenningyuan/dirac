C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
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
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : labread.f
C SOURCE  : $Source$
C
      PROGRAM READLABL
C
C Written by Hans Jorgen Aa. Jensen, Nov. 1983
C CRAY version 10. Oct. 1986
C Alliant version 30. Dec. 1987
C Linux/IBM, +updates, Joern Thyssen, Nov. 1999
C
#include <implicit.h>
CTROND      DOUBLE PRECISION TIME, SECOND
      CHARACTER*8 B(4), STARS
      PARAMETER (STARS = '********')
      PARAMETER (LUINP = 1, LUPRI = 6)
      INTEGER   NARG,IARG
      CHARACTER*80 ARGC
C
C      TIME = SECOND()
      NARG = IARGC()
      IF (NARG .EQ. 0) THEN
         CALL GETARG(0,ARGC)
         WRITE(6,'(/A,A,A/)') 'Usage: ',ARGC,' file1 file2 ...'
         STOP ' '
      END IF
      DO IARG = 1,NARG
         CALL GETARG(IARG,ARGC)
         WRITE (LUPRI, '(/A,A/)') 
     &      ' MOLECULE labels found on file :',ARGC
C
         OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
         REWIND LUINP
         IREC = 0
         IERR = 0
    1    READ (LUINP,END=3,ERR=2) B
            IREC = IREC + 1
         IF (B(1) .NE. STARS) GO TO 1
            WRITE (LUPRI, '(5X,I5,3X,4(2X,A8))')  IREC, B
         GO TO 1
C
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
C
    3    CONTINUE
C      TIME = SECOND() - TIME
         WRITE (LUPRI,'(/I10,A)') IREC,
     *      ' records read before EOF on file.'
         WRITE (LUPRI,'(F10.2,A)') TIME,' CPU seconds used.'
         CLOSE(LUINP)
      END DO
      STOP ' '
C
      END
