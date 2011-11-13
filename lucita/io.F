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
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DISKUN
      IMPLICIT REAL*8           (A-H,O-Z)
*
* Assign logical unit numbers for LUCIA:
*==================================================
*  Complete revision motivated by MCSCF development
*==================================================
*
*   1 -  9 :  Retained for I/O
*  ---------
*        5 :  Standard IN
*        6 :  Standard OUT
*
*  60 - 69 :  Important files retained (DIAG,C,HC etc.)
*  ---------
*       60 :    CI diagonal
*       61 :    CI reference vector
*       62 :    CI sigma vector
*
*  70 - 79 :  Scratch files
*  ---------
*
*  80 - 99 :  Files in use in different context
*             (e.g. standalone, INS, MOLCAS etc.)
*  ---------
*
*  Timo Fleig    Mar 28, 2002
*==================================================
*
*  Old setup:
*
* All file with some kind of input information  :  10 - 19
* All files containing final results            :  90 - 99
* Scratch files                                 :  30 - 40
* Internal files (retained through job)         :  20 - 29
*
#include "clunit.inc"
* =========================
* Standard input and output
* =========================
*. Input file
      LUIN = 5
*. Special unit for input file from MOLCAS conversion routine
*    Timo Fleig
      LUCMOL = 91
*. Output file
      LUOUT = 6
* =================
* Input information
* =================
* Input file containing MO-AO transformation matrix
      LUMOIN = 80
*. Input file for CI-vectors
*. restart from file 21 is assumed
*. Input , two electron integrals - MOLCAS
      LU2INT = 81
*. Input , one electron integrals - MOLCAS
      LU1INT = 82
*. Input , property one-electron integral files
      LUPRP  = 83
*. Sirius interface file
      LUSIR1 = 84
*. File containing additional states for transition densities
      LUEXC = 85
* =================
* Internal files
* =================
*. CI diagonal
      LUDIA = 60
      Open(Unit=LUDIA,File='CIDIA',Status='UNKNOWN',Form='UNFORMATTED')
*. CI vector
      LUC = 61
      Open(Unit=LUC,File='LUCVECT',Status='UNKNOWN',Form='UNFORMATTED')
*. Sigma vector file
      LUHC = 62
      Open(Unit=LUHC,File='HCFILE',Status='UNKNOWN',Form='UNFORMATTED')
* =================
* Scratch files
* =================
      LUSC1 = 70
      LUSC2 = 71
      LUSC3 = 72
      Open(Unit=LUSC1,File='LUSC1',Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC2,File='LUSC2',Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC3,File='LUSC3',Status='UNKNOWN',Form='UNFORMATTED')
*. Scratch space for subspace handling
      LUSC34 = 73
      LUSC35 = 74
      LUSC36 = 75
      LUSC37 = 76
      LUSC38 = 77
      LUSC39 = 78
      LUSC40 = 79
      Open(Unit=LUSC34,File='LUSC34',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC35,File='LUSC35',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC36,File='LUSC36',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC37,File='LUSC37',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC38,File='LUSC38',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC39,File='LUSC39',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC40,File='LUSC40',Status='UNKNOWN',
     &     Form='UNFORMATTED')
* =================
* Output files
* =================
*. output file for CI-vectors
*. Not in use
      LUCIVO = 86
*. Natural orbitals in terms of input orbitals
*.
      LUMOUT = 87
*     Open(Unit=LUMOUT,File='CIAORB',Status='UNKNOWN',Form='FORMATTED')
*. Dumping 1- and 2- electron integrals in formatted form
* ( LU90 just defined here, it is not exported )
      LU90  = 88
*. Dumping symmmetry info, MO-AO expansion matrix and property integral
      LU91 = 69
*. CC amplitudes in formatted form
      LU_CCAMP = 90
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      subroutine setunits
*
#include "clunit.inc"
#include "units.inc"
*
* Scratch file
      LUSC_INT = LU91
* Input for LUCITA, converted
      LUINP = LUSC1
* COMFILE from MOLCAS-SCF
      LUCOM = LUSC36
* JOBIPH file from MOLCAS-RASSCF
      LUIPH = LUSC37
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TODSCN(VEC,NREC,LREC,LBLK,LU)
*
* Write VEC as multiple record file accordin to NREC and LREC
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC(*)
      INTEGER LREC(NREC)
*
      IOFF = 1
      DO IREC = 1, NREC
C?      WRITE(6,*) ' TODSCN: IREC, LREC ',IREC,LREC(IREC)
C?      WRITE(6,*) ' Input record '
C?      CALL WRTMAT(VEC(IOFF),1,LREC(IREC),1,LREC(IREC))
        IF(LREC(IREC).GE.0) THEN
          CALL ITODS(LREC(IREC),1,LBLK,LU)
          CALL TODSC(VEC(IOFF),LREC(IREC),LBLK,LU)
          IOFF = IOFF + LREC(IREC)
        ELSE
          CALL ITODS(-LREC(IREC),1,LBLK,LU)
          CALL ZERORC(IDUMMY,LU,0)
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TODSCNP(VEC,NREC,LREC,LBLK,LU)
*
* Write VEC as multiple record file accordin to NREC and LREC
* Only nonzero elements are written to DISC
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC(*)
      INTEGER LREC(NREC)
*
      IOFF = 1
      DO IREC = 1, NREC
C?      WRITE(6,*) ' TODSCN: IREC, LREC ',IREC,LREC(IREC)
C?      WRITE(6,*) ' Input record '
C?      CALL WRTMAT(VEC(IOFF),1,LREC(IREC),1,LREC(IREC))
        IF(LREC(IREC).GE.0) THEN
*. Normal complete record
          CALL ITODS(LREC(IREC),1,LBLK,LU)
          CALL TODSCP(VEC(IOFF),LREC(IREC),LBLK,LU)
          IOFF = IOFF + LREC(IREC)
        ELSE
*. zero record of length -LREC(IREC)
          CALL ITODS(-LREC(IREC),1,LBLK,LU)
          CALL ZERORC(IDUMMY,LU,1)
        END IF

*
      END DO
*
      RETURN
      END
