
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck copvcd $Revision$
      SUBROUTINE COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
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
            CALL ZERORC(LBL,LUOUT,IAMPACK)
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck copvcdp $Revision$
      SUBROUTINE COPVCDP(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C
C COPY VECTOR ON FILE LUIN TO FILE LUOUT
*
* Packed version
C
C
C LBLK DEFINES STRUCTURE OF FILE
C Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
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
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck diskun $Revision$
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
*. Dumping symmmetry info, MO-AO expansion matrix and property integrals
      LU91 = 89
*. CC amplitudes in formatted form
      LU_CCAMP = 90
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck frmdsc $Revision$
      SUBROUTINE FRMDSC(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
*. Version allowing zero and packed blocks
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRAY(*)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)

C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        MMBLOCK = MBLOCK
C       IF(MMBLOCK.GE.2) MMBLOCK = 2
C       CALL IFRMDS(ISCR,2,MMBLOCK,IFILE)
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(I_AM_PACKED.NE.0) THEN
C       WRITE(6,*) ' File is packed, file number = ', IFILE
        END IF
        IF(IMZERO.EQ.1) THEN
          ZERO = 0.0D0
C?        write(6,*) ' frmdsc, length of zero block',NDIM
          CALL SETVEC(ARRAY,ZERO,NDIM)
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF(I_AM_PACKED.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(ARRAY,ZERO,NDIM)
*. Loop over packed records of dimension LPBLK
      NBATCH = 0
C1000 CONTINUE
*. The next LPBLK elements
  999   CONTINUE
          NBATCH = NBATCH + 1
          IF(NBATCH.NE.1) THEN
            LBATCHP = LBATCH
          END IF
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE) (IPAK(I),I=1, LBATCH)
            READ(IFILE) (XPAK(I),I=1, LBATCH)
          END IF
          READ(IFILE) ISTOP
          DO IELMNT = 1, LBATCH
            IF(IPAK(IELMNT).LE.0.OR.IPAK(IELMNT).GT.NDIM) THEN
              WRITE(6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT
              WRITE(6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT )
              WRITE(6,*) ' LBATCH IFILE  = ',LBATCH,IFILE
              IF(NBATCH.EQ.1) THEN
               WRITE(6,*) ' NBATCH = 1 '
              ELSE
               WRITE(6,*) ' NBATCH, LBATCHP', NBATCH,LBATCHP
              END IF
              WRITE(6,*) ' NDIM,IMZERO = ', NDIM,IMZERO
              Call Abend2(' problem in FRMDSC ')
            END IF
            ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
          END DO
        IF(ISTOP.EQ.0) GOTO 999
*. End of loop over records of truncated elements
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
         ELSE
          READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
        END IF
C
        IF( MBLOCK.LT.0.AND.NDIM.GT.0.AND.ICRAY.EQ.0 ) THEN
         CALL SQFILE(IFILE,2,ARRAY,2*NDIM)
        END IF
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck frmdsc2 $Revision$
      SUBROUTINE FRMDSC2(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED,
     &                   NO_ZEROING)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
*. Version allowing zero and packed blocks
*
* If NO_ZEROING = 1, the elements of zero blocks
*    are not set to zero, the routine just returns with
*    IMZERO = 1
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRAY(*)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)

      IMZERO = 0
C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        MMBLOCK = MBLOCK
C       IF(MMBLOCK.GE.2) MMBLOCK = 2
C       CALL IFRMDS(ISCR,2,MMBLOCK,IFILE)
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(IMZERO.EQ.1) THEN
          IF(NO_ZEROING.EQ.0) THEN
            ZERO = 0.0D0
            CALL SETVEC(ARRAY,ZERO,NDIM)
          END IF
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF(I_AM_PACKED.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(ARRAY,ZERO,NDIM)
*. Loop over packed records of dimension LPBLK
      NBATCH = 0
C1000 CONTINUE
*. The next LPBLK elements
  999   CONTINUE
          NBATCH = NBATCH + 1
          IF(NBATCH.NE.1) THEN
            LBATCHP = LBATCH
          END IF
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE) (IPAK(I),I=1, LBATCH)
            READ(IFILE) (XPAK(I),I=1, LBATCH)
          END IF
          READ(IFILE) ISTOP
          DO IELMNT = 1, LBATCH
            IF(IPAK(IELMNT).LE.0.OR.IPAK(IELMNT).GT.NDIM) THEN
              WRITE(6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT
              WRITE(6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT )
              WRITE(6,*) ' LBATCH IFILE  = ',LBATCH,IFILE
              IF(NBATCH.EQ.1) THEN
               WRITE(6,*) ' NBATCH = 1 '
              ELSE
               WRITE(6,*) ' NBATCH, LBATCHP', NBATCH,LBATCHP
              END IF
              WRITE(6,*) ' NDIM,IMZERO = ', NDIM,IMZERO
              Call Abend2(' problem in FRMDSC ')
            END IF
            ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
          END DO
        IF(ISTOP.EQ.0) GOTO 999
*. End of loop over records of truncated elements
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
         ELSE
          READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
        END IF
C
        IF( MBLOCK.LT.0.AND.NDIM.GT.0.AND.ICRAY.EQ.0 ) THEN
         CALL SQFILE(IFILE,2,ARRAY,2*NDIM)
        END IF
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck frmdscn3 $Revision$
      SUBROUTINE FRMDSCN3(VEC,NREC,LBLK,LU,NO_ZEROING,I_AM_ZERO,
     &                    LBLOCK)
*
* Read  VEC as multiple record file, NREC records read
* If NO_ZEROING.EQ.0 then zero blocks are not
* set to zero;  a 1 is instead flagged in the relevant block
* of I_AM_ZERO
*
* IF LBLOCK(IBLK) is less than zero, then no space for this
* block is made
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION LBLOCK(*)
*. Output
      DIMENSION VEC(*),I_AM_ZERO(*)
*
      IOFF = 1
C?    write(6,*) ' FRMDSCN3: Number of records to be read', NREC
      DO IREC = 1, NREC
        CALL IFRMDS(LREC,1,LBLK,LU)
        IF(LBLOCK(IREC).GE.0) THEN
          CALL FRMDSC2(VEC(IOFF),LREC,LBLK,LU,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          I_AM_ZERO(IREC) = IMZERO
          IOFF = IOFF + LREC
        ELSE
          I_AM_ZERO(IREC) = 1
          CALL SKPRCD2(LBLK,-1,LU)
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
*fordeck ifrmds $Revision$
      SUBROUTINE IFRMDS(IARRAY,NDIM,MBLOCK,IFILE)
C
C     TRANSFER INTEGER ARRAY FROM DISC FILE IFILE
C
C NBLOCK .LT. 0 INDICATES USE OF FASTIO
C
C If nblock .eq. 0 NBLOCK = NDIM
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IARRAY(1)
C
      ICRAY = 1
      NBLOCK = MBLOCK

      IF( ICRAY.EQ.1.OR.NBLOCK .GE. 0 ) THEN
C       DO NOT USE FASTIO
        IF(NBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
          IF(IREST.GT.NBLOCK) THEN
            READ(IFILE) (IARRAY(IBASE+I),I=1,NBLOCK)
            IBASE=IBASE+NBLOCK
            IREST=IREST-NBLOCK
          ELSE
            READ(IFILE) (IARRAY(IBASE+I),I=1,IREST)
            IREST=0
          END IF
        IF( IREST .GT. 0 ) GOTO 100
      ELSE
C       USE FAST IO
        CALL SQFILE(IFILE,2,IARRAY,NDIM)
      END IF
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck itods $Revision$
      SUBROUTINE ITODS(IA,NDIM,MBLOCK,IFIL)
C TRANSFER ARRAY INTEGER IA(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IA(1)
      INTEGER START,STOP
*
      ICRAY = 1
      NBLOCK = MBLOCK
      IF( NBLOCK .GE.0.OR.ICRAY.EQ.1 ) THEN
C
      IF(NBLOCK .LE. 0 ) NBLOCK = NDIM
      STOP=0
      NBACK=NDIM
C LOOP OVER RECORDS
  100 CONTINUE
       IF(NBACK.LE.NBLOCK) THEN
         NTRANS=NBACK
         NLABEL=-NTRANS
       ELSE
         NTRANS=NBLOCK
         NLABEL=NTRANS
       END IF
       START=STOP+1
       STOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       WRITE(IFIL) (IA(I),I=START,STOP),NLABEL
      IF(NBACK.NE.0) GOTO 100
      END IF
C
      IF(ICRAY.EQ.0.AND. NBLOCK .LT. 0 ) THEN
       CALL SQFILE(IFIL,1,IA,NDIM)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck nexrec $Revision$
      SUBROUTINE NEXREC(LU,LBLK,REC,IEND,LENGTH)
*
* OBTAIN NEXT RECORD ON FILE LU, IF
* AN END OF FILE IS ISSUED THE RECORD IS EMPTY
* AND IEND IS SET TO 1
*
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION REC(*)
*
      IF(LBLK .GT. 0 ) THEN
        IEND = 0
        LENGTH = LBLK
COLD    CALL FRMDSC(REC,LENGTH,LBLK,LU)
        CALL FRMDSC(REC,LENGTH,-1,LU,IMZERO,IAMPACK)
      ELSE
        CALL IFRMDS(LENGTH,1,LBLK,LU)
C?    write(6,*) ' Length in NEXREC ',LENGTH
        IF(LENGTH.GE.0) THEN
          IEND = 0
          CALL FRMDSC(REC,LENGTH,LBLK,LU,IMZERO,IAMPACK)
C?     write(6,*) ' Record read in '
C?     CALL WRTMAT(REC,1,LENGTH,1,LENGTH)
        ELSE
          IEND = 1
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck rewine $Revision$
      SUBROUTINE REWINE( LU ,LBLK )
*
* LBLK .LT. 0 :  REWIND SEQ FILE LU WITH FASTIO ROUTINES
* LBLK .GE. 0 :  rewinf seq file LU with normal REWIND
      ICRAY = 1
      IF ( ICRAY.EQ.0.AND.LBLK .LT. 0 ) THEN
        IDUM = 1
        CALL SQFILE(LU,5,IDUM,IDUM)
      ELSE
        REWIND LU
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck rewino $Revision$
      SUBROUTINE REWINO( LU )
C
C REWIND SEQ FILE LU WITH FASTIO ROUTINES
C
C?    WRITE(6,*) ' TO REWIND FILE ',LU
      IDUM = 1
C     CALL SQFILE(LU,5,IDUM,IDUM)
      REWIND (LU)
C?    WRITE(6,*) ' FILE REWOUND '
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck setunits $Revision$
*
      subroutine setunits
*
#include "clunit.inc"
#include "units.inc"
*
* Scratch file
      LUSC_INT = LU91
      Open(Unit=LUSC_INT,
     &     File='LUSCDL',
     &     Status='UNKNOWN',
     &     Form='FORMATTED')
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
*fordeck skprcd2 $Revision$
      SUBROUTINE SKPRCD2(NDIM,MBLOCK,IFILE)
C
C     Skip record in file IFILE
C
*. Version allowing zero and packed blocks
*
* Dos not work with FASTIO - I expect
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)

C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(IMZERO.EQ.1) THEN
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF(I_AM_PACKED.EQ.1) THEN
*. Loop over packed records of dimension LPBLK
*. The next LPBLK elements
  999   CONTINUE
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE)
            READ(IFILE)
          END IF
          READ(IFILE) ISTOP
        IF(ISTOP.EQ.0) GOTO 999
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
         ELSE
          READ(IFILE)
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
        END IF
C
C       IF( MBLOCK.LT.0.AND.NDIM.GT.0.AND.ICRAY.EQ.0 ) THEN
C        CALL SQFILE(IFILE,2,ARRAY,2*NDIM)
C       END IF
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck skpvcd $Revision$
      SUBROUTINE SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
C
C SKIP OVER NVEC VECTORS ON FILE LUIN
C
C LBLK DEFINES STRUCTURE OF FILE
C (see note on structure of files )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        CALL REWINE(LU ,LBLK)
      END IF
      DO 1001 IVEC = 1, NVEC
C
C LOOP OVER BLOCKS OF GIVEN VECTOR
C
 1000   CONTINUE
C
          IF( LBLK .GT. 0 ) THEN
            LBL = LBLK
          ELSE IF (LBLK .EQ. 0 ) THEN
            READ(LU) LBL
          ELSE IF (LBLK .LT. 0 ) THEN
            CALL IFRMDS(LBL,1,-1,LU)
          END IF
C
          IF( LBL .GE. 0 ) THEN
            IF(LBLK .GE.0 ) THEN
              KBLK = LBLK
            ELSE
              KBLK = -1
            END IF
            CALL FRMDSC(SEGMNT,LBL,KBLK,LU,IMZERO,IAMPACK)
          END IF
        IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
 1001 CONTINUE
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck todsc $Revision$
      SUBROUTINE TODSC(A,NDIM,MBLOCK,IFIL)
C TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1)
      INTEGER START,STOP
      REAL*8 INPROD
      INTEGER ISCR(2)
*
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = INPROD(A,A,NDIM)
        IF(XNORM.EQ.0.0D0) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. No packing
         ISCR(2) = 0
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
      ICRAY = 1
      IF( MBLOCK .GE.0 .OR.ICRAY .EQ. 1 ) THEN
C
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
      STOP=0
      NBACK=NDIM
C LOOP OVER RECORDS
  100 CONTINUE
       IF(NBACK.LE.NBLOCK) THEN
         NTRANS=NBACK
         NLABEL=-NTRANS
       ELSE
         NTRANS=NBLOCK
         NLABEL=NTRANS
       END IF
       START=STOP+1
       STOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       WRITE(IFIL) (A(I),I=START,STOP),NLABEL
      IF(NBACK.NE.0) GOTO 100
      END IF
C
      IF( ICRAY.EQ.0.AND.MBLOCK.LT.0.AND.NDIM.GT.0) THEN
       CALL SQFILE(IFIL,1,A,2*NDIM)
      END IF
*
 1001 CONTINUE
C
C?    write(6,*) ' leaving TODSC '
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck todscn $Revision$
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
*fordeck todscnp $Revision$
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
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*fordeck todscp $Revision$
      SUBROUTINE TODSCP(A,NDIM,MBLOCK,IFIL)
*
C TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
*
* Packed version : Store only nonzero elements
*. Small elements should be xeroed outside
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1)
      INTEGER START,STOP
      REAL*8 INPROD
      INTEGER ISCR(2)
*
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)
*
*
C?    write(6,*) ' entering TODSCP, file = ', IFIL
C?    CALL XFLUSH(6)
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = INPROD(A,A,NDIM)
        IF(XNORM.EQ.0.0D0) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. Packing
        ISCR(2) = 1
C       CALL ITODS(ISCR,2,MMBLOCK,IFIL)
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
      ICRAY = 1
      IF( MBLOCK .GE.0 .OR.ICRAY .EQ. 1 ) THEN
C
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
*. Loop over packed records of dimension LPBLK
      IELMNT = 0
 1000 CONTINUE
*. The next LPBLK elements
      LBATCH = 0
*. Obtain next batch of elemnts
  999 CONTINUE
       IF(NDIM.GE.1) THEN
       IELMNT = IELMNT+1
       IF(A(IELMNT).NE.0.0D0) THEN
         LBATCH=LBATCH+1
         IPAK(LBATCH) = IELMNT
         XPAK(LBATCH) = A(IELMNT)
       END IF
       END IF
       IF(LBATCH.EQ.LPBLK.OR.IELMNT.EQ.NDIM) goto 998
       GOTO 999
*. Send to DISC
 998   CONTINUE
       WRITE(IFIL) LBATCH
       IF(LBATCH.GT.0) THEN
         WRITE(IFIL) (IPAK(I),I=1, LBATCH)
         WRITE(IFIL) (XPAK(I),I=1, LBATCH)
       END IF
       IF(IELMNT.EQ.NDIM) THEN
         WRITE(IFIL) -1
       ELSE
         WRITE(IFIL) 0
         GOTO 1000
       END IF
*. End of loop over records of truncated elements
      END IF
 1001 CONTINUE
*
C?    CALL XFLUSH(6)
      RETURN
      END
