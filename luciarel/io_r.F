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
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE COPNBLKD(LUIN,LUOUT,SEGMNT,NBLOCK,IREW,LBLK)
*
* Copy Nblocks from FILE LUIN to file LUOUT
*
*
* Type of file LUOUT is inherited from LUIN
*
* If less than NBLOCK blocks are left in the given vector
* on LUIN a warning is printed,
* but normal return
*
*
* NO EOV mark written
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        CALL REWINE( LUIN ,LBLK)
        CALL REWINE( LUOUT ,LBLK)
      END IF

      DO IIBLK = 1, NBLOCK
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          IF(LBL.GE.0) CALL ITODS (LBL,1,-1,LUOUT)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IAMPACK.EQ.0) THEN
            CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
          ELSE
            CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
          END IF
        END IF
        IF(LBL.EQ.-1) THEN
           WRITE(6,*) ' END OF VECTOR mark for BLOCK'
           write(6,*) 'IBLOCK not defined!!'
           call quit('Stop in COPNBLKD')
           GOTO 3333
        END IF
      END DO
*     ^End of loop over blocks
 3333 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_REL(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C
C COPY VECTOR ON FILE LUIN TO FILE LUOUT
C
C
C LBLK DEFINES STRUCTURE OF FILE
*
* Structure of output file is inherited by output file,
* if input file is packed, so is output file
*
*
C Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        CALL REWINE( LUIN ,LBLK)
        CALL REWINE( LUOUT ,LBLK)
      END IF
 
C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUOUT)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &         NO_ZEROING)
          IF(IAMPACK.NE.0) THEN
C?          WRITE(6,*) ' COPVCD, IAMPACK,FILE = ', IAMPACK,LUIN
          END IF
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
            ELSE  
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            CALL ZERORC_REL(LBL,LUOUT,IAMPACK)
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCDC(LUIN,LUOUT,VEC1,IREW,IRI,LBLK)
*
* Copy complex vector from LUIN to LUOUT
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION VEC1(*)
*
      IF(IREW.NE.0) THEN
        CALL REWINE(LUIN,LBLK)
        CALL REWINE(LUOUT,LBLK)
      END IF
*
*. Real part
      CALL COPVCD_REL(LUIN,LUOUT,VEC1,0,LBLK)
*
      if (IRI.eq.2) then
*. Imaginary part
        CALL COPVCD_REL(LUIN,LUOUT,VEC1,0,LBLK)
      end if
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE DISKUN_REL
      IMPLICIT REAL*8           (A-H,O-Z)
#include "ipoist8.inc"
*
* Assign logical unit numbers for LUCIA:
*
*==================================================
*  Complete revision motivated by MCSCF development
*==================================================
*
*   1 -  9 :  Retained for I/O
*  ---------
*        5 :  Standard IN
*        6 :  Standard OUT
*        7 :  Standard OUT (Slaves)
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
*  Timo Fleig    Mar 19, 2002
*  changes for parallel business: Stefan Knecht   Apr 24, 2007
*==================================================
*
*  Old setup:
*
* All files with some kind of input information  :  10 - 19
* Internal files (retained through job)          :  20 - 29
* Scratch files                                  :  30 - 40
* Files for machine information			 :  50 - 59
* Additional scratch files		       	 :  41,61,62
* Files for relativistic business (SO ints etc.) :  70 - 79
* All files containing final results             :  90 - 99
*
*
#include "clunit.inc"
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "lucipar_rel.inc"
#include "files_r.inc"
      INTEGER LUSCR_MPIL,LUDIA_MPIL,LUC_MPIL, LUHC_MPIL, LUSC123_MPIL
      INTEGER LUSC344_MPIL, LUDUMPFLL
      CHARACTER*9 LUDUMPFN
      CHARACTER*13 NLUDUMPF
*
CSK      WRITE(6,*) 'I have been in LUCI_SETUP',MYPROC
*
* =========================
* Standard input and output
* =========================
*. Input file
      LUIN = 5
      LU5 = LUIN
*. Output file
      LUOUT = 6
* =================
* Input information
* =================
* Input file containing MO-AO transformation matrix
      LUMOIN = 80
*. Input file for CI-vectors
*. restart from file 21 is assumed
*. Input , two electron integrals
      LU2INT = 81
      LU40 = LU2INT
*. Input , one electron integrals
      LU1SO = 82
      LU1INT = 83
      LU2 = LU1INT        ! for MOLCAS environment
      LU9 = LU1INT        ! for MOLECULE environment
*. Input , property one-electron integral files
      LUPRP  = 84
*. Sirius interface file
      LUSIR1 = 85
*. Input MOAO coefficients
      LU19 = 86
*. Input Spin-orbit integrals (X)
      LU69 = 87
*. Input Spin-orbit integrals (X)
      LU70 = 88
*. Input Spin-orbit integrals (X)
      LU71 = 89
* =================
* Internal files
* =================
*. CI diagonal
      LUDIA = 60
      LUDIA_MPIN ="CIDIA_REL"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUDIA_MPI,'(A9,A1,I1)') LUDIA_MPIN,'.',MYPROC
         LUDIA_MPIL=11
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUDIA_MPI,'(A9,A1,I2)') LUDIA_MPIN,'.',MYPROC
         LUDIA_MPIL=12
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUDIA_MPI,'(A9,A1,I3)') LUDIA_MPIN,'.',MYPROC
         LUDIA_MPIL=13
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUDIA_MPI=LUDIA_MPIN
      LUDIA_MPIL=9
#endif
*
      Open(Unit=LUDIA,File=NLUDIA_MPI(1:LUDIA_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
*. CI vector
      LUC =   61
      LUC_MPIN ="LUCVECT_REL"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUC_MPI,'(A11,A1,I1)') LUC_MPIN,'.',MYPROC
         LUC_MPIL=13
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUC_MPI,'(A11,A1,I2)') LUC_MPIN,'.',MYPROC
         LUC_MPIL=14
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUC_MPI,'(A11,A1,I3)') LUC_MPIN,'.',MYPROC
         LUC_MPIL=15
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUC_MPI=LUC_MPIN
      LUC_MPIL=11
#endif
      Open(Unit=LUC,File=NLUC_MPI(1:LUC_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
*. Sigma vector file
      LUHC =  62
      LUHC_MPIN ="HCFILE_REL"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUHC_MPI,'(A10,A1,I1)') LUHC_MPIN,'.',MYPROC
         LUHC_MPIL=12
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUHC_MPI,'(A10,A1,I2)') LUHC_MPIN,'.',MYPROC
         LUHC_MPIL=13
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUHC_MPI,'(A10,A1,I3)') LUHC_MPIN,'.',MYPROC
         LUHC_MPIL=14
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUHC_MPI=LUHC_MPIN
      LUHC_MPIL=10
#endif
      Open(Unit=LUHC,File=NLUHC_MPI(1:LUHC_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
* =================
* Scratch files
* =================
      LUSC1 = 70
      LUSC1_MPIN ="LUSC1_REL"
      LUSC2 = 71
      LUSC2_MPIN ="LUSC2_REL"
      LUSC3 = 72
      LUSC3_MPIN ="LUSC3_REL"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC1_MPI,'(A9,A1,I1)') LUSC1_MPIN,'.',MYPROC
         WRITE (NLUSC2_MPI,'(A9,A1,I1)') LUSC2_MPIN,'.',MYPROC
         WRITE (NLUSC3_MPI,'(A9,A1,I1)') LUSC3_MPIN,'.',MYPROC
         LUSC123_MPIL=11
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC1_MPI,'(A9,A1,I2)') LUSC1_MPIN,'.',MYPROC
         WRITE (NLUSC2_MPI,'(A9,A1,I2)') LUSC2_MPIN,'.',MYPROC
         WRITE (NLUSC3_MPI,'(A9,A1,I2)') LUSC3_MPIN,'.',MYPROC
         LUSC123_MPIL=12
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC1_MPI,'(A9,A1,I3)') LUSC1_MPIN,'.',MYPROC
         WRITE (NLUSC2_MPI,'(A9,A1,I3)') LUSC2_MPIN,'.',MYPROC
         WRITE (NLUSC3_MPI,'(A9,A1,I3)') LUSC3_MPIN,'.',MYPROC
         LUSC123_MPIL=13
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC1_MPI=LUSC1_MPIN
      NLUSC2_MPI=LUSC2_MPIN
      NLUSC3_MPI=LUSC3_MPIN
      LUSC123_MPIL=9
#endif
      Open(Unit=LUSC1,File=NLUSC1_MPI(1:LUSC123_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC2,File=NLUSC2_MPI(1:LUSC123_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC3,File=NLUSC3_MPI(1:LUSC123_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
*. Scratch space for subspace handling
      LUSC34 = 73
      LUSC34_MPIN ="LUSC34_REL"
      LUSC35 = 74
      LUSC35_MPIN ="LUSC35_REL"
      LUSC36 = 75
      LUSC36_MPIN ="LUSC36_REL"
      LUSC37 = 76
      LUSC37_MPIN ="LUSC37_REL"
      LUSC38 = 77
      LUSC38_MPIN ="LUSC38_REL"
      LUSC39 = 78
      LUSC39_MPIN ="LUSC39_REL"
      LUSC41 = 79
      LUSC41_MPIN ="LUSC41_REL"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC34_MPI,'(A10,A1,I1)') LUSC34_MPIN,'.',MYPROC
         WRITE (NLUSC35_MPI,'(A10,A1,I1)') LUSC35_MPIN,'.',MYPROC
         WRITE (NLUSC36_MPI,'(A10,A1,I1)') LUSC36_MPIN,'.',MYPROC
         WRITE (NLUSC37_MPI,'(A10,A1,I1)') LUSC37_MPIN,'.',MYPROC
         WRITE (NLUSC38_MPI,'(A10,A1,I1)') LUSC38_MPIN,'.',MYPROC
         WRITE (NLUSC39_MPI,'(A10,A1,I1)') LUSC39_MPIN,'.',MYPROC
         WRITE (NLUSC41_MPI,'(A10,A1,I1)') LUSC41_MPIN,'.',MYPROC
         LUSC344_MPIL=12
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC34_MPI,'(A10,A1,I2)') LUSC34_MPIN,'.',MYPROC
         WRITE (NLUSC35_MPI,'(A10,A1,I2)') LUSC35_MPIN,'.',MYPROC
         WRITE (NLUSC36_MPI,'(A10,A1,I2)') LUSC36_MPIN,'.',MYPROC
         WRITE (NLUSC37_MPI,'(A10,A1,I2)') LUSC37_MPIN,'.',MYPROC
         WRITE (NLUSC38_MPI,'(A10,A1,I2)') LUSC38_MPIN,'.',MYPROC
         WRITE (NLUSC39_MPI,'(A10,A1,I2)') LUSC39_MPIN,'.',MYPROC
         WRITE (NLUSC41_MPI,'(A10,A1,I2)') LUSC41_MPIN,'.',MYPROC
         LUSC344_MPIL=13
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC34_MPI,'(A10,A1,I3)') LUSC34_MPIN,'.',MYPROC
         WRITE (NLUSC35_MPI,'(A10,A1,I3)') LUSC35_MPIN,'.',MYPROC
         WRITE (NLUSC36_MPI,'(A10,A1,I3)') LUSC36_MPIN,'.',MYPROC
         WRITE (NLUSC37_MPI,'(A10,A1,I3)') LUSC37_MPIN,'.',MYPROC
         WRITE (NLUSC38_MPI,'(A10,A1,I3)') LUSC38_MPIN,'.',MYPROC
         WRITE (NLUSC39_MPI,'(A10,A1,I3)') LUSC39_MPIN,'.',MYPROC
         WRITE (NLUSC41_MPI,'(A10,A1,I3)') LUSC41_MPIN,'.',MYPROC
         LUSC344_MPIL=14
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC34_MPI=LUSC34_MPIN
      NLUSC35_MPI=LUSC35_MPIN
      NLUSC36_MPI=LUSC36_MPIN
      NLUSC37_MPI=LUSC37_MPIN
      NLUSC38_MPI=LUSC38_MPIN
      NLUSC39_MPI=LUSC39_MPIN
      NLUSC41_MPI=LUSC41_MPIN
      LUSC344_MPIL=10
#endif
      Open(Unit=LUSC34,File=NLUSC34_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC35,File=NLUSC35_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC36,File=NLUSC36_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC37,File=NLUSC37_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC38,File=NLUSC38_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC39,File=NLUSC39_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC41,File=NLUSC41_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
*
*. Reference vectors for (transition) densities
*
      LUSC61 = 65
C     LUSC62 = 66
C  avoid unit 66. Used in DALTON
CSK      LUSC62 = 95    ! conflict might occur with LU18 in genh1.F
      LUSC62 = 67    ! might this work? SK April 24, 2007
      LUSC61_MPIN ="LUSC61_REL"
      LUSC62_MPIN ="LUSC62_REL"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC61_MPI,'(A10,A1,I1)') LUSC61_MPIN,'.',MYPROC
         WRITE (NLUSC62_MPI,'(A10,A1,I1)') LUSC62_MPIN,'.',MYPROC
         LUSC344_MPIL=12
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC61_MPI,'(A10,A1,I2)') LUSC61_MPIN,'.',MYPROC
         WRITE (NLUSC62_MPI,'(A10,A1,I2)') LUSC62_MPIN,'.',MYPROC
         LUSC344_MPIL=13
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC61_MPI,'(A10,A1,I3)') LUSC61_MPIN,'.',MYPROC
         WRITE (NLUSC62_MPI,'(A10,A1,I3)') LUSC62_MPIN,'.',MYPROC
         LUSC344_MPIL=14
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC61_MPI=LUSC61_MPIN
      NLUSC62_MPI=LUSC62_MPIN
      LUSC344_MPIL=10
#endif
      Open(Unit=LUSC61,File=NLUSC61_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN', Form='UNFORMATTED')
      Open(Unit=LUSC62,File=NLUSC62_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN', Form='UNFORMATTED')
*. Some more general purpose scratch
      LUSC92 = 92
      LUSC93 = 93
      LUSC92_MPIN ="LUSC92_REL"
      LUSC93_MPIN ="LUSC93_REL"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC92_MPI,'(A10,A1,I1)') LUSC92_MPIN,'.',MYPROC
         WRITE (NLUSC93_MPI,'(A10,A1,I1)') LUSC93_MPIN,'.',MYPROC
         LUSC929_MPIL=12
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC92_MPI,'(A10,A1,I2)') LUSC92_MPIN,'.',MYPROC
         WRITE (NLUSC93_MPI,'(A10,A1,I2)') LUSC93_MPIN,'.',MYPROC
         LUSC929_MPIL=13
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC92_MPI,'(A10,A1,I3)') LUSC92_MPIN,'.',MYPROC
         WRITE (NLUSC93_MPI,'(A10,A1,I3)') LUSC93_MPIN,'.',MYPROC
         LUSC929_MPIL=14
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC92_MPI=LUSC92_MPIN
      NLUSC93_MPI=LUSC93_MPIN
      LUSC929_MPIL=10
#endif
*
      LUSC94 = 94
      LUSC95 = 95
      LUSC94_MPIN ="INTEGRALS"
      LUSC95_MPIN ="INDICES"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC94_MPI,'(A9,A1,I1)') LUSC94_MPIN,'.',MYPROC
         WRITE (NLUSC95_MPI,'(A7,A1,I1)') LUSC95_MPIN,'.',MYPROC
         LUSC94_MPIL=11
         LUSC95_MPIL=9
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC94_MPI,'(A9,A1,I2)') LUSC94_MPIN,'.',MYPROC
         WRITE (NLUSC95_MPI,'(A7,A1,I2)') LUSC95_MPIN,'.',MYPROC
         LUSC94_MPIL=12
         LUSC95_MPIL=10
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC94_MPI,'(A9,A1,I3)') LUSC94_MPIN,'.',MYPROC
         WRITE (NLUSC95_MPI,'(A7,A1,I3)') LUSC95_MPIN,'.',MYPROC
         LUSC94_MPIL=13
         LUSC95_MPIL=11
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUSC94_MPI=LUSC94_MPIN
      NLUSC95_MPI=LUSC95_MPIN
      LUSC94_MPIL=9
      LUSC95_MPIL=7
#endif
*. open and close of LUSC92 - LUSC95 at relevant places ...
      LUSC96 = 96     ! File MRCONEE
      LUSC97 = 97     ! File MDCINT
      LUSC98 = 98
      LUDUMPF = 99
      LUDUMPFN = "LDUMPFILE"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUDUMPF,'(A9,A1,I1)') LUDUMPFN,'.',MYPROC
         LUDUMPFLL=11
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUDUMPF,'(A9,A1,I2)') LUDUMPFN,'.',MYPROC
         LUDUMPFLL=12
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUDUMPF,'(A9,A1,I3)') LUDUMPFN,'.',MYPROC
         LUDUMPFLL=13
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUDUMPF = LUDUMPFN
      LUDUMPFLL=9
#endif
      Open(Unit=LUDUMPF,File=NLUDUMPF(1:LUDUMPFLL),
     &     Status='UNKNOWN', Form='FORMATTED')
*           
* =================
* Output files
* =================
*. Natural spinors in terms of double group adapted atomic spin-orbitals
*.
      LUMOUTR = 90
      LUMOUTI = 91
      LUMOUTR_MPIN ="NATSPIR" 
      LUMOUTI_MPIN ="NATSPII" 
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUMOUTR_MPI,'(A7,A1,I1)') LUMOUTR_MPIN,'.',MYPROC
         WRITE (NLUMOUTI_MPI,'(A7,A1,I1)') LUMOUTI_MPIN,'.',MYPROC
         LUMOUTRI_MPIL=9
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUMOUTR_MPI,'(A7,A1,I2)') LUMOUTR_MPIN,'.',MYPROC
         WRITE (NLUMOUTI_MPI,'(A7,A1,I2)') LUMOUTI_MPIN,'.',MYPROC
         LUMOUTRI_MPIL=10
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUMOUTR_MPI,'(A7,A1,I3)') LUMOUTR_MPIN,'.',MYPROC
         WRITE (NLUMOUTI_MPI,'(A7,A1,I3)') LUMOUTI_MPIN,'.',MYPROC
         LUMOUTRI_MPIL=11
      ELSE
        CALL QUIT("NMPROC.GT.1000! DISKUN_REL: That is too much for me")
      ENDIF
#else
      NLUMOUTR_MPI = LUMOUTR_MPIN
      NLUMOUTI_MPI = LUMOUTI_MPIN
      LUMOUTRI_MPIL = 7
#endif
*. open and close of LUMOUTR and LUMOUTI at relevant places ...
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE LUCI_CLOSEDOWN()
*
*     close LUCIAREL files and keep them
*
*     Stefan Knecht   April 24, 2007
*
      IMPLICIT REAL*8           (A-H,O-Z)
#include "ipoist8.inc"
#include "clunit.inc"
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "lucipar_rel.inc"
*
CSK      WRITE(6,*) 'I have been in LUCI_CLOSEDOWN',MYPROC
      CLOSE(UNIT=LUDIA,STATUS='KEEP')
      CLOSE(UNIT=LUC,STATUS='KEEP')
      CLOSE(UNIT=LUHC,STATUS='KEEP')
      CLOSE(UNIT=LUSC1,STATUS='KEEP')
      CLOSE(UNIT=LUSC2,STATUS='KEEP')
      CLOSE(UNIT=LUSC3,STATUS='KEEP')
      CLOSE(UNIT=LUSC34,STATUS='KEEP')
      CLOSE(UNIT=LUSC35,STATUS='KEEP')
      CLOSE(UNIT=LUSC36,STATUS='KEEP')
      CLOSE(UNIT=LUSC37,STATUS='KEEP')
      CLOSE(UNIT=LUSC38,STATUS='KEEP')
      CLOSE(UNIT=LUSC39,STATUS='KEEP')
      CLOSE(UNIT=LUSC41,STATUS='KEEP')
      CLOSE(UNIT=LUSC61,STATUS='DELETE')
      CLOSE(UNIT=LUSC62,STATUS='DELETE')
      CLOSE(UNIT=LUDUMPF,STATUS='KEEP')
CSK      CLOSE(UNIT=LUDUMPF,STATUS='DELETE')
*
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine opnone
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Use vector partitioning information to save a vector
* to LUCIAREL file.
*
*    Timo Fleig, March 16, 2002
*
      subroutine savevcdc(VECT,SCR,LU,IRILP,IPRNT)
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
*
#include "krmcluci_inf.h"
*
      dimension VECT(*),SCR(*)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      ZERO = 0.D0
*
      if (NTEST.ge.5) then
        write(6,*)
        write(6,*) ' +++++++++++++++++++++++++++++++ '
        write(6,*) '  SAVEVCDC called.               '
        write(6,*) '   (Save complex vec. to file.)  '
        write(6,*) ' +++++++++++++++++++++++++++++++ '
        write(6,*) '  File structure information:'
        write(6,*) '  Number of MK values :',NMK2VAL
        JBATABS = 0
        do IMK2 = 1,NMK2VAL,1
          write(6,*) '  MK value no.: ',IMK2
          do JBATS = 1,NBAT_MS2_MC(IMK2)
            JBATABS = JBATABS + 1
            write(6,*) '   Absolute batch no.: ',JBATABS
            write(6,*) '   No. of elements in batch :',
     &                     LEBATS_MC(JBATABS)
          end do
        end do 
      end if
*
      JOFF = 1
      LENGTH = 0
*
      DO IRI = 1,IRILP,1
* Write end of vector mark between real and imaginary part.
        if (IRI.eq.2) call itods(-1,1,-1,LU)
*
        JBATABS = 0
        if (NTEST.ge.10) write(6,*) '  REAL/IMAG = ',IRI
        DO IMK2 = 1, NMK2VAL
*
* Loop over batches of sigma for this MK2 projection value
          do JBATS = 1,NBAT_MS2_MC(IMK2)
            JBATABS = JBATABS + 1
*
            if (NTEST.ge.10) then
               WRITE(6,*)
               WRITE(6,*) ' Start of writing of batch ', JBATABS
               write(6,'(2X,A,I3)') ' for MK2 value number ',IMK2
               WRITE(6,*) 
            end if
*
            LS = LEBATS_MC(JBATABS)
            LENGTH = LENGTH + LS
            if (NTEST.ge.10.and.LS.gt.0) then
              WRITE(6,*) ' **********************************'
              WRITE(6,*) ' * Array containing current batch *'
              WRITE(6,*) ' **********************************'
              CALL WRTMAT(VECT(JOFF),1,LS,1,LS)
            end if
            KBATSEND = I1BATS_MC(JBATABS) + LBATS_MC(JBATABS)
            do ISBLK = I1BATS_MC(JBATABS),KBATSEND-1,1
              IOFF = IBATS_MC(6,ISBLK)
              LEN_V = IBATS_MC(8,ISBLK)
              if (NTEST.ge.100) then
                write(6,*) '    Saving block no. ',ISBLK
                write(6,*) '     Block offset:   ',IOFF
                write(6,*) '     Block length:   ',LEN_V
              end if
              call itods(LEN_V,1,-1,LU)
              call todsc(VECT(JOFF+IOFF-1),LEN_V,-1,LU)
            end do
            JOFF = JOFF + LS
*
          end do
*         ^ End loop over S batches of given MK2 value
        end do
*       ^ End of loop over IMK2
*
      end do
*     ^ End of loop over IRI
* Write end of vector mark.
      call itods(-1,1,-1,LU)
*
      RETURN
      END
