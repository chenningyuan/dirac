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
      SUBROUTINE COPVCDC(LUIN,LUOUT,VEC1,IREW,LBLK)
*
* Copy complex vector from LUIN to LUOUT
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*)
*
      IF(IREW.NE.0) THEN
        CALL REWINE(LUIN,LBLK)
        CALL REWINE(LUOUT,LBLK)
      END IF
*
*. Real part
      CALL COPVCD_REL(LUIN,LUOUT,VEC1,0,LBLK)
*. Imaginary part
      CALL COPVCD_REL(LUIN,LUOUT,VEC1,0,LBLK)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE DISKUN_REL
      IMPLICIT REAL*8           (A-H,O-Z)
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
#include "clunit.inc"
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
*. CI vector
      LUC =   61
*. Sigma vector file
      LUHC =  62
* =================
* Scratch files
* =================
      LUSC1 = 70
      LUSC2 = 71
      LUSC3 = 72
*. Scratch space for subspace handling
      LUSC34 = 73
      LUSC35 = 74
      LUSC36 = 75
      LUSC37 = 76
      LUSC38 = 77
      LUSC39 = 78
      LUSC41 = 79
*. Reference vectors for (transition) densities
      LUSC61 = 65
C     LUSC62 = 66
C  avoid unit 66. Used in DALTON
      LUSC62 = 95    ! conflict might occur with LU18 in genh1.F
*. Some more general purpose scratch
      LUSC92 = 92
      LUSC93 = 93
      LUSC94 = 94
      LUSC95 = 95
      LUSC96 = 96
      LUSC97 = 97
      LUSC98 = 98
      LUSC99 = 99
* =================
* Output files
* =================
*. Natural spinors in terms of double group adapted atomic spin-orbitals
*.
      LUMOUTR = 90
      LUMOUTI = 91
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine opnone
      return
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
      subroutine savevcdc(VECT,LU,IPRNT)
*
#include "implicit.inc"
#include "mxpdim.inc"
*
#include "civecinf.h"
*
      dimension VECT(*)
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
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
      IRILP = 2
      DO IRI = 1,IRILP,1
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
* Write end of vector mark between real and imaginary part.
        call itods(-1,1,-1,LU)
*
      end do
*     ^ End of loop over IRI
* Write end of vector mark.
      call itods(-1,1,-1,LU)
*
      RETURN
      END
