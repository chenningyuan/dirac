      SUBROUTINE FNDMND(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,
     &                  SCR,ISCAT,SUBVAL,NTESTG)
*
* FIND NSUB LOWEST ELEMENTS OF VECTOR RESIDING ON FILE
* LU. ENSURE THAT NO DEGENERENCIES ARE SPLIT
*
*
* INPUT
*======
* LU : FILE WHERE VECTOR OF INTEREST IS RESIDING, REWOUND
* LBLK : DEFINES FILE STRUCTURE ON FILE LU
* NSUBMX: LARGEST ALLOWED NUMBER OF SORTED ELEMENTS
*
*OUTPUT
*======
* NSUB : ACTUAL NUMBER OF ELEMENTS OBTAINED. CAN BE SMALLER
*       THAN NSUBMX IF THE LATST ELEMENT BELONGS TO A DEGENERATE
*       SET
*ISCAT: SCATTERING ARRAY , ISCAT(I) GIVED FULL ADRESS OF SORTED
*       ELEMENT I
*       SUBVAL : VALUE OF SORTED ELEMENTS

      IMPLICIT REAL*8           ( A-H,O-Z)
      DIMENSION SEGMNT(*), ISCAT(*),SUBVAL(*),SCR(*),ISCR(*)
*
      CALL REWINE(LU,-1)
*
*.LOOP OVER BLOCKS
*
C     write(6,*) ' FNDMND NSUBMX = ', NSUBMX
*
      NTESTL = 0000
      NTEST = MAX(NTESTG,NTESTL)
*
      IBLK = 0
      IBASE = 1
      LSUB = 0
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          CALL IFRMDS(LBL,1,-1,LU)
        END IF
        IBLK = IBLK + 1
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Info about block ',IBLK
          WRITE(6,*) ' Number of elements ', LBL
        END IF
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(SEGMNT,LBL ,-1,LU,IMZERO,IAMPACK)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Elements read in '
             CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
           END IF
           IF(LBL .GE. 0 ) THEN
*. LOWEST ELEMENTS IN SEGMNT  ( ADD TO PREVIOUS LIST )
             MSUBMX = MIN(NSUBMX,LBL)
C     SUBROUTINE SORLOW(WRK,STVAL,ISTART,KZVAR,KEXSTV,JEXSTV,IPRT)
             IF(LBL.GE.1) THEN
               CALL SORLOW(SEGMNT,SCR(1+LSUB),ISCR(1+LSUB),LBL,
     &                     MSUBMX,MSUB,NTEST)
C              WRITE(6,*)
C    &       ' After SORLOW 1 , Scatter array for combined list '
C            call iwrtma(ISCR(1),1,LSUB+MSUB,1,LSUP+MSUP)
             ELSE
              MSUB = 0
             END IF
             DO 10 I = 1, MSUB
   10        ISCR(LSUB+I) = ISCR(LSUB+I) + IBASE - 1
C            Write(6,*)
C    &       ' After 10 , Scatter array for combined list '
C            call iwrtma(ISCR,1,LSUB+MSUB,1,LSUP+MSUP)
* SORT COMBINED LIST
             MSUBMX = MIN(NSUBMX,LSUB+MSUB)
             IF(MSUBMX.GT.0) THEN
             CALL SORLOW(SCR,SUBVAL,ISCAT,LSUB+MSUB,MSUBMX,LLSUB,
     &         NTEST)
C              WRITE(6,*)
C    &       ' After SORLOW 2 , Scatter array for combined list '
C            call iwrtma(ISCR(1),1,LSUB+MSUB,1,LSUP+MSUP)
             ELSE
               LLSUB = 0
             END IF
             LSUB = LLSUB
             DO 20 I = 1, LSUB
               ISCR(I+2*NSUBMX) = ISCR(ISCAT(I))
   20        CONTINUE
*
             CALL ICOPVE(ISCR(1+2*NSUBMX),ISCR(1),LSUB)
             CALL COPVEC(SUBVAL,SCR,LSUB)
             IBASE = IBASE + LBL
             IF(NTEST .GE. 20 ) THEN
               WRITE(6,*) ' Lowest elements and their original place '
               WRITE(6,*) ' Number of elements obtained ', LSUB
               CALL WRTMAT(SUBVAL,1,LSUB,1,LSUB)
               CALL IWRTMA(ISCR,1,LSUB,1,LSUB)
             END IF
          END IF
*
        END IF
*
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      NSUB = LSUB
      CALL ICOPVE(ISCR,ISCAT,NSUB)
      IF(NTEST .GE. 10) THEN
        WRITE(6,*) ' Lowest elements and their original place '
        WRITE(6,*) ' Number of elements obtained ', NSUB
       CALL WRTMAT(SUBVAL,1,NSUB,1,NSUB)
       CALL IWRTMA(ISCAT,1,NSUB,1,NSUB)
      END IF
*
      RETURN
      END
