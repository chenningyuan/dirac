C
C$Id$
C
C
C*************************************************************
C* This is the procedures that controls the reading from     *
C* file, and calls all the other subroutines.                *
C* It returns a bunch of variables. It returns:              *
C*     IQM - (Highest anglar quantumnumber) + 1              *
C*     JCO - Redundant variable that gives how many blocks   *
C*           the primitives and contractioncoeffesients are  *
C*           given in (will be set to one) (array)           *
C*     NUC - Number of primitives for a given                *
C*           (quantumnuber + 1)                              *
C*     NRC - Number of coloumns with contraction-coeffesients*
C*           for a given (quantumnumber + 1).                *
C*     SEG - Gives whether the contraction-coeffesients are  *
C*           segmented (for a given (quantumnumber + 1),     *
C*           (array).                                        *
C*     ALPHA - Array with the primitives (for a given        *
C*             (quantumnumber + 1)).                         *
C*     CPRIMU - Matrix (3D) with the NOT normalized          *
C*              contraction coeffesients (for a given        *
C*              (quantumnumber + 1)).                        *
C*     CPRIM - Like CPRIMU, but this time the contraction    *
C*             -coeffesients are normalized.                 *
*     NBLOCK - For some reason is this defined as the same  *
C*              as IQM                                       *
C*     KAOVEC, KPRIM, Q, DSM are max (angular quantumnumber  *
C*             + 1), max number of primitives, nuclear       *
C*             charge, highest number to be accepted as zero *
C*************************************************************
C
      SUBROUTINE BASLIB (IQM, JCO, NUC, NRC, SEG, ALPHA, CPRIM,
     &                   CPRIMU, NBLOCK, KAOVEC, KPRIM, Q, ISET,
     &                   DSM)
#include <implicit.h>
#include <priunit.h>
      LOGICAL SEG
      DIMENSION JCO(IQM), NUC(KAOVEC), NRC(KAOVEC), SEG(KAOVEC),
     +          ALPHA(KPRIM, KAOVEC), CPRIMU(KPRIM, KPRIM, KAOVEC), 
     +          CPRIM(KPRIM,KPRIM, KAOVEC)
#include <cbirea.h>
C*************************************************************
C* Variable declarations:                                    *
C*     BASNAM - A character variable that contains the name  *
C*              of the file.                                 *
C*     NEWEL - Logical variable, that is returned from       *
C*             the subroutine finpos, gives whether the next *
C*             line in basis-file is a new element or not.   *
C*     SEGEJ - This gives whether the contraction            *
C*             -coeffesients are segmented or not (for       *
C*             (angular quantum-number + 1))                 *
C* External procedures:                                      *
C*     SEGORB - A procedure that checks if the contraction   *
C*              -coeffesients (for a given (angular quantum  *
C*              number + 1)) are segmented or not.           *
C*     NRMORB - Normalizing a matrix with contraction        *
C*              -coeffesients                                *
C*************************************************************
C
      CHARACTER*13 BASTMP
      LOGICAL NEWEL, SEGEJ
C
      INTQ = INT(Q)
      NBLOCK = 0
      IAUG = 0
      SEGEJ = .TRUE.
      IF (ISET .EQ. 2) THEN
         BASTMP = 'STO-3G       '
      ELSE
         IF (BASNAM(2:10) .EQ. 'aug-cc-pV') THEN
            BASTMP(1:11) = BASNAM(2:12)
            IF (BASNAM(1:1) .EQ. 'd') THEN
               IAUG = 1
            ELSE IF (BASNAM(1:1) .EQ. 't') THEN
               IAUG = 2
            ELSE IF (BASNAM(1:1) .EQ. 'q') THEN
               IAUG = 3
            ELSE
               WRITE (LUPRI,'(/A1,A)') BASNAM(1:1),' is an unknown'//
     &              'augmentation level'
               CALL QUIT('Illegal augmentation level in HERBAS')
            END IF
         ELSE
            BASTMP = BASNAM
         END IF
      END IF
C      
C Finds the right element in the file. 
C
      CALL FINDEL(BASTMP, INTQ)
C
C Finds the number of primitives and coloumns of contraction coeffesients
C
      CALL FINPOS ( NEWEL, INTEXP, INTORB)
C
 10   CONTINUE
      IF ( .NOT. NEWEL) THEN
         NBLOCK = NBLOCK + 1
C
C Setting ALPHA, CPRIMU AND CPRIM to zero. 
C
         CALL DZERO(ALPHA(1, NBLOCK), KPRIM)
         CALL DZERO(CPRIMU(1,1, NBLOCK), KPRIM*KPRIM)
         CALL DZERO(CPRIM (1, 1, NBLOCK), KPRIM*KPRIM)
C
C Reading the primitives and contraction-coeffesients from file in READNU
C
         CALL READNU( INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU, ALPHA,
     +                KPRIM,IAUG)
         IF (INTORB .EQ. 0) INTORB = INTEXP
C
C Checking if the matrix, whith the contraction-coefesients, are segmented 
C
         CALL SEGORB( SEGEJ, INTEXP, INTORB, CPRIMU(1,1,NBLOCK), KPRIM, 
     +                DSM)
         NUC(NBLOCK) = INTEXP
         NRC(NBLOCK) = INTORB
         SEG(NBLOCK) = SEGEJ
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
         CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +               ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +               CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
C
C Define the type of the next relevant line in file
C
         CALL FINPOS ( NEWEL, INTEXP, INTORB)
         GOTO 10
      END IF
C
      IQM = NBLOCK
C
C Setting JCO to one.
C
      CALL SETONE(JCO, IQM)
C
C Closing file
C
      CLOSE (20, STATUS = 'KEEP', IOSTAT = IOERR, ERR = 1000)
      RETURN
 1000 CONTINUE
      WRITE (LUPRI,'(/A,I4)') 'Error in closing file. Errorcode: ',IOERR
C
      END
C
C*************************************************************
C*            The subroutine that sets JCO as one            *
C*                 (remnant from the past)                   *
C*************************************************************
C
C/* Deck setone */
      SUBROUTINE SETONE (JCO, IQM)
#include <implicit.h>
      DIMENSION JCO(IQM)
C
      DO 100 I = 1, IQM
         JCO(I) = 1
 100  CONTINUE
      END
C
C*************************************************************
C* This is the subroutine that takes care of the reading of  *
C* the primitives og the contraction coeffesients from the   *
C* file (for a given format). The variables that are         *
C* transferred out are:                                      *
C*     CPRIMU - The matrix that the contraction coeffesients *
C*              are put in (not normalized).                 *
C*     ALPHA - The matrix where the primitives are put in    *
C*     (Both of the variables are given for a given (quantum-*
C*     number + 1))                                          *
C*************************************************************
C
C/* Deck readnu */
      SUBROUTINE READNU(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU, ALPHA, 
     +                  KPRIM,IAUG)
#include <implicit.h>
#include <priunit.h>
C
C Reading the exponents
C
      DIMENSION CPRIMU( KPRIM, KPRIM, KAOVEC), ALPHA(KPRIM, KAOVEC)
C*************************************************************
C* Variable declarations:                                    *
C*     SIGN - Character that gives what kind of sentence it  *
C*            is.                                            *
C*     LINE - A character array that helps bulletproofing the*
C*            sucbroutine.                                   *
C*     STRING - A character string that helps bulletproofing *
C*             the subroutine                                *
C*     BLANK - A logical variable that gives whether the line*
C*             is blank or not.                              *
C*************************************************************
      CHARACTER SIGN, LINE
      CHARACTER*88 STRING
      DIMENSION LINE(88)
      LOGICAL BLANK 
C
C CONEOR is the maximum number of coloumns of contraction coeffesients on 
C the 1. line in file. If INTORB is .GT. CONEOR then the contraction coeffesients
C will continue on the next line. CONORB is the max. number of contraction
C coeffesients on the next lines.
C
      PARAMETER(CONEOR = 6, CONORB = 7)
C
C Initializing a counting variable that must be 0 before loop
C
      J = 0
C
 20   CONTINUE
      IF ( J .LT. INTEXP) THEN
C
C Reading the primitive and contracted coeffesients
C
         READ(20, '(A88)', IOSTAT = IOERR, ERR = 2000)
     +        STRING
         READ (STRING, '(A1)') SIGN
         IF (SIGN .EQ. ' ') THEN
C
C Calling a subroutine that checks whether the line is a blank one. If
C STRING is a blank line BLANK will return as .TRUE.
C
            CALL CHBLANK(BLANK, STRING)
            IF (.NOT. BLANK) THEN
C
C We have found a line with a primitive and coeffesients               
C
               J = J + 1
C
C CONLIN returns the number of lines we have the contraction coeffesients written on.
C
               IF (INTORB .EQ. 0) THEN
                  READ (STRING, '(F16.9)') ALPHA(J,NBLOCK)
                  CPRIMU(J,J,NBLOCK) = 1.0D0
               ELSE
                  CALL CONLIN(INTORB, NUMLIN)
C     
C     Getting the format for the read statement right.
C     
                  IF (NUMLIN .EQ. 1) THEN
                     KNTORB = INTORB
                  ELSE
                     KNTORB = CONEOR
                  END IF
C     
C     Reading the first line with exponents and contraction- coeffesients
C
                  READ (STRING, '(F16.9, 6F12.9)')
     +                 ALPHA(J, NBLOCK), (CPRIMU(J,I,NBLOCK),
     +                 I = 1, KNTORB)
C     
C     If there are more lines with contraction-coeffesients, 
C     they will be read here.
C
                  DO 100 I = 2, NUMLIN
                     NUMNUM = CONEOR + (I-1)*CONORB
C     
C Getting the format for the read-stat right.
C
                     KNTORB = MIN0(NUMNUM, INTORB)
 30                  CONTINUE
C
C     Making the usual safety-precautions before we read the 
C     contraction-coeffesients.
C
                     READ(20, '(A88)', IOSTAT = IOERR, ERR = 2000)
     +                    STRING
                     READ (STRING, '(A1)') SIGN
                     IF (SIGN .EQ. ' ') THEN
                        CALL CHBLANK(BLANK, STRING)
                        IF (.NOT. BLANK) THEN
C     
C     We now have a line with contraction-coeffesients.
C     
                           READ (STRING,'(7F12.9)') (CPRIMU(J,M,NBLOCK),
     +                          M = CONEOR + (I-2)*CONORB +1, KNTORB)
                        END IF
C     
C     Blank line read the next one
C
                     ELSE
                        GOTO 30
                     END IF
 100              CONTINUE
                  GOTO 20
               END IF
            ELSE 
C
C blank line, read the next one.
C
               GOTO 20
            END IF         
         ELSE
C
C Found a line with nothing sensible in it       
C
            GOTO 20
         END IF
      ELSE
C     
C Read all of the primitives, and the coeffesients
C If this is some kind of augmented basis set, we augment it here, kr-96
C
         IF (IAUG .GT. 0) THEN
            DO 199 KAUG = 1, IAUG
               INTEXP = INTEXP + 1
               INTORB = INTORB + 1
               ALPHA(INTEXP,NBLOCK) = ALPHA(INTEXP-1,NBLOCK)
     &                               *ALPHA(INTEXP-1,NBLOCK)
     &                               /ALPHA(INTEXP-2,NBLOCK)
               CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0
 199        CONTINUE 
         END IF
         RETURN
      END IF
C
 2000 CONTINUE
      WRITE (LUPRI,'(/A,I4)') 'Error in reading from file. Errorcode: ',
     &     IOERR
      END
C
C*************************************************************
C* This is a subroutine that checks whether a line is a blank*
C* line. It returns a logical variable if it is a blank one. *
C* The in/out-going variables are:                          *
C*     STRING - A character string (the string checked by the*
C*              subroutine).                                 *
C*     BLANK - The logical variable. It will be .TRUE. if the*
C*             line is a blank one.                          *
C*************************************************************
C/* Deck chblank */
      SUBROUTINE CHBLANK(BLANK, STRING)
C
      PARAMETER(LENGTH = 88)
      CHARACTER*88 STRING
      LOGICAL BLANK
C*************************************************************
C* Variable declarations:                                    *
C*     LINE - A character array that contains all the        *
C*            characters in STRING (one by one).             *
C*************************************************************
      CHARACTER LINE
      DIMENSION LINE(LENGTH)
C
      BLANK = .TRUE.
C
      READ (STRING, '(88A1)') (LINE(I), I = 1, LENGTH)
      DO 100 J = 1, LENGTH
         BLANK = BLANK .AND. (LINE(J) .EQ. ' ')
 100  CONTINUE
      END
C
C*************************************************************
C* This is a subroutine that finds out how many lines the    *
C* contraction- coeffesients are written on.  There are 6    *
C* cc (contraction-coeffesients) on the 1. line, and then    *
C* there are seven on the rest of the lines. In out going    *
C* variables are :                                           *
C*      INTORB - Total number of coloumns of contraction     *
C*               -coeffesients.                              *
C*      NUMLIN - Returns the number of lines the contaction  *
C*               -coeffesients are written on.               *
C*************************************************************
C/* Deck conlin */
C
      SUBROUTINE CONLIN(INTORB, NUMLIN)
#include <implicit.h>
C      
C This is the parameters for the number of contraction-coeffesients on 1. and 2. line.
C
      PARAMETER(CONEOR = 6, CONORB = 7)
C
C The intrisic functions DBLE makes a souble precision reak number of an integer.
C
      B = DBLE(CONORB)
      C = DBLE(INTORB) - DBLE(CONEOR) 
C
C This finds out how many lines we have, and puts it into numlin.
C
      IF ((INTORB - CONEOR) .LE. 0) THEN
         NUMLIN = 1
      ELSE IF (DMOD(C,B) .LT. DSM) THEN
         NUMLIN = (INTORB - CONEOR)/CONORB + 1
      ELSE
         NUMLIN = (INTORB - CONEOR)/CONORB + 2
      END IF
      END
C********************************************************
C* This are the subroutines that  reads  the primitives *
C* and the contracted coeffesients from file.           *
C*********************************************************      
C* This is the subroutine that searches through the file *
C* and finds the element in question.                    *
C*********************************************************
C/* Deck findel */
      SUBROUTINE FINDEL(BASNAM, INTQ)
#include <implicit.h>
#include <priunit.h>
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basisfile                *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*     EXST - Logical variable that helps inquire if     *
C*            there exists a file with that name.        *
C*********************************************************
C
      CHARACTER*(*) BASNAM
      CHARACTER*88 STRING
      CHARACTER SIGN
      LOGICAL EXST
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
C
      STRING = BASDIR(1:LENBAS)//BASNAM
      INQUIRE (FILE = STRING, EXIST = EXST)
      IF (EXST) THEN
         OPEN(20, FILE = STRING, IOSTAT = IOERR,
     +        ERR = 1000)
C
C
C Searching the file for the element.
C
 20      CONTINUE
         READ(20,'(A88)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)') SIGN
C
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
            READ (STRING, 110) TEGN, NUCEL
 110        FORMAT ( A1, I4)
            IF (INTQ .EQ. NUCEL) THEN
               RETURN
            ELSE
               GOTO 20
            END IF
         ELSE 
            GOTO 20
         END IF   
      ELSE
         WRITE (LUPRI,'(/,3A)') 'Basis ',STRING,' doesn''t exist'
         CALL QUIT('Non-existing basis set in HERBAS')
      END IF
C
C Errormessages 
C
 1000 CONTINUE
      WRITE (LUPRI,'(/A)') 'Error when opening or closing file'
      CALL QUIT('I/O error in HERBAS')
 2000 CONTINUE
      WRITE (LUPRI,'(/A)') 'Error when reading from file'
      CALL QUIT('I/O error in HERBAS')
 200  CONTINUE
      WRITE (LUPRI,'(/I3,2A)') INTQ,
     &     ' is an unsupported element for basis ',BASNAM
      CALL QUIT('Unsupported elements in HERBAS')
      END 
C               
C
C*********************************************************	
C* This is the subroutine that finds out wheter it is a  *
C* new block with primitives and contraction-            *
C* coeffesients, or if it is a new element. It returns   *
C* three variables:                                      *
C*          NEWEL -Logical variable that gives whether it*
C*                 is a new element or not.              *
C*          INEXP, INTORB - Number of primitives and     *
C*                          coloumns of contraction      *
C*                          coeffesients                 *
C*********************************************************
C
C/* Deck finpos */
      SUBROUTINE FINPOS(NEWEL, INTEXP, INTORB)
#include <implicit.h>
#include <priunit.h>
C*********************************************************
C* Variable declarations:                                *
C*     SIGN - The first sign in a sentence, to find out  *
C*            what kind of sentence it is.               *
C*     STRING - A character string to help bulletproofing*
C*              the subroutine.                          *
C*     LINE - Character-array to help bulletproofing the *
C*            subroutine.                                *
C*     BLANK - A logical variable to help bulletproofing *
C*             the subroutine.
C*********************************************************
      LOGICAL NEWEL
      CHARACTER SIGN
      CHARACTER*88 STRING
      LOGICAL BLANK
      CHARACTER LINE
      DIMENSION LINE(88)
C
C Initialising NEWEL, because one is innocent until proven guilty
C
      NEWEL = .FALSE.
C
 20   CONTINUE
      READ (20, '(A88)', IOSTAT = IOERR, ERR = 1000, END = 200)
     +     STRING
      READ (STRING, '(A1)') SIGN
C
      IF (SIGN .EQ. ' ') THEN
C
C Calling on a subroutine that checks whether the line is a blank one or not.
C If STRING is a blank line, then BLANK will return as .TRUE.
         CALL CHBLANK(BLANK, STRING)
         IF (BLANK) THEN
            GOTO 20
         ELSE
            READ (STRING, '(2I5)') INTEXP, INTORB
            RETURN
         END IF
      ELSE IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN  
         NEWEL = .TRUE.
         RETURN
      ELSE 
         GOTO 20
      END IF
C
C No more orbitals when end of file, so the orbitals for the last element 
C are read
C
 200  CONTINUE
      NEWEL = .TRUE.
C
C NEWEL must be .TRUE. to break the if-loop.
C
      RETURN
C
 1000 CONTINUE
      WRITE (LUPRI,'(/A)') 'Error in reading file, your basis '//
     &     'is not complete'
      CALL QUIT('Incomplete basis set file in HERBAS')
      END
