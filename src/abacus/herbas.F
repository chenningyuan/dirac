!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

C
C
C
C FILE : herbas.F
C
#include "single.h"
C
C*************************************************************
C* This is the procedures that controls the reading from     *
C* file, and calls all the other subroutines.                *
C* It returns a bunch of variables. It returns:              *
C*     IQM - (Highest angular quantum number) + 1            *
C*     JCO - Variable that gives how many blocks             *
C*           the primitives and contraction coefficients are *
C*           given in (array)                                *
C*     NUC - Number of primitives for a given                *
C*           (quantum number + 1)                            *
C*     NRC - Number of columns with contraction coefficients *
C*           for a given (quantum number + 1).               *
C*     SEG - Gives whether the contraction coefficients are  *
C*           segmented (for a given (quantum number + 1),    *
C*           (array).                                        *
C*     ALPHA - Array with the primitives (for a given        *
C*             (quantum number + 1)).                        *
C*     CPRIMU - Matrix (3D) with the NOT normalized          *
C*              contraction coefficients (for a given        *
C*              (quantum number + 1)).                       *
C*     CPRIM - Like CPRIMU, but this time the contraction    *
C*             coefficients are normalized.                  *
C*     NBLOCK - Number of function blocks                    *
C*     KAOVEC, KPRIM, Q, DSM are max (angular quantum number *
C*             + 1), max number of primitives, nuclear       *
C*             charge, highest number to be accepted as zero *
C*     BASREF - documentation lines for basis set            *
C*     QEFF is eff. nuclear charge which is .le. Q,          *
C*             Q-QEFF will be described with ECP             *
C*************************************************************
C
      SUBROUTINE BASLIB (IQM, JCO, NUC, NRC, SEG, ALPHA, CPRIM,
     &                   CPRIMU, NBLOCK, ISGEN, KAOVEC, KPRIM, 
     &                   Q, QEFF, ISET, DSM, BASREF)
#include "implicit.h"
#include "priunit.h"
      LOGICAL   SEG, READNU_CALLED
      DIMENSION JCO(*), NUC(KAOVEC), NRC(KAOVEC), SEG(KAOVEC),
     +          ALPHA(KPRIM, KAOVEC), CPRIMU(KPRIM, KPRIM, KAOVEC), 
     +          CPRIM(KPRIM,KPRIM, KAOVEC), ISGEN(KAOVEC)
#include "cbirea.h"
#include "dcbgen.h"
#include "dcbham.h"
C*************************************************************
C* Variable declarations:                                    *
C*     BASNAM - A character variable that contains the name  *
C*              of the file with the basis set.              *
C*     NEWEL - Logical variable, that is returned from       *
C*             the subroutine finpos, gives whether the next *
C*             line in basis-file is a new element or not.   *
C*     SEGEJ - This gives whether the contraction            *
C*             coefficients are segmented or not (for        *
C*             (angular quantum-number + 1))                 *
C* External procedures:                                      *
C*     SEGORB - A procedure that checks if the contraction   *
C*              coefficients (for a given (angular quantum   *
C*              number + 1)) are segmented or not.           *
C*     NRMORB - Normalizing a matrix with contraction        *
C*              coefficients                                 *
C*************************************************************
C
      CHARACTER BASTMP*80, BASSAV*80
      CHARACTER*80 BASREF(10)
      LOGICAL NEWEL, SEGEJ, ANO, SADLEJ, POLFUN, FOUND, BLANK, NQVD
#ifdef PRG_DIRAC
C     Define a maximum number of primitives that can go in one block
C     for UNCONTRACTED basis sets. The value given is for s-functions
C     that generate 3 small component p-functions. The maximum for
C     other types of functions is scaled down by KHK(IQM+1) that
C     gives the number of functions in the small component for that type.
      PARAMETER (MAX_IN_BLOCK=3*15)
      LOGICAL SPLIT,TEST
      integer :: icount = 0
#endif /* PRG_DIRAC */
C
      CALL QENTER('BASLIB')
C

      BASTMP = ' '
      BASSAV = BASNAM
      INTQ   = NINT(Q)
      LUBAS = 20
      IQM    = 0
      NBLOCK = 0
      IAUG   = 0
      ANO    = .FALSE.
      SADLEJ = .FALSE.
      NQVD   = .FALSE.
      POLFUN = .FALSE.
CDIRAC:..allow contraction only for elements Z.LE.36
C      except for 2c NR calculations
      TEST = (TWOCOMP.AND.(.NOT.TWOCOMPBSS)).OR.LEVYLE
      TEST = UNCONT.OR.((INTQ.GT.36).AND.(.NOT.TEST))
      IF(TEST.AND.(.NOT.UNCONT)) THEN
        WRITE(LUPRI,'(/A/A,I5/2A/A/A)')
     &   '*** WARNING (BASLIB) : Decontracting basis set',
     &   ' - nuclear charge: ',INTQ,
     &   ' - basis file    : ',BASNAM,
     &   'Reason: DIRAC cannot yet create proper small component basis',
     &   '        for contracted large component basis when Z > 36'
      ENDIF

      IF (ISET .EQ. 2) THEN
C jan 2000 hjaaj:
#ifdef STO3G_HUCKEL
         BASTMP = 'Huckel       '
#else
C
C hjaaj jan 2000:
C        New Huckel based on ano-4 instead of STO-3G.
C        - Advantages: atomic orbitals are orthonormal,
C                      better description of atomic shells
C        - The number of contracted functions is based
C        on the available Huckel parameters in huckel.F
C        (it would be good to extend with 2p parameters for Li, Be etc.)
C        BASNAM is needed here but contains true basis set name,
C        this is restored from BASSAV just before RETURN.
C
         ANO = .TRUE.
         BASTMP = 'ano-4'
         ISTART = 6
         IF (INTQ .LE. 2) THEN
C           ... H and He
            BASNAM = 'ano-4   1 0 0 0'
         ELSE IF (INTQ .LE. 4) THEN
C           ... Li, Be
            BASNAM = 'ano-4   2 0 0 0'
         ELSE IF (INTQ .LE. 10) THEN
C           ... B, C, N, O, F, Ne
            BASNAM = 'ano-4   2 1 0 0'
         ELSE IF (INTQ .LE. 12) THEN
C           ... Na, Mg
            BASNAM = 'ano-4   3 1 0 0'
         ELSE IF (INTQ .LE. 18) THEN
C           ... Al, Si, P, S, Cl, Ar
            BASNAM = 'ano-4   3 2 0 0'
         ELSE IF (INTQ .LE. 20) THEN
C           ...  K, Ca
            BASNAM = 'ano-4   4 2 0 0'
         ELSE IF (INTQ .LE. 30) THEN
C           ...  Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
            BASNAM = 'ano-4   4 2 1 0'
         ELSE IF (INTQ .LE. 36) THEN
C           ...  Ga, Ge, As, Se, Br, Kr
            BASNAM = 'ano-4   4 3 1 0'
         ELSE
            CALL QUIT('BASLIB error: Huckel not defined for Z.gt.36')
         END IF
#endif   /* STO3G-HUCKEL */
C djw Feb 2005
C notation changed from 'd' to 'd-' for consistency with EMSL
      ELSE
         IF (BASNAM(3:10) .EQ. 'aug-cc-p') THEN
            BASTMP(1:78) = BASNAM(3:80)
            IF (BASNAM(1:2) .EQ. 'd-') THEN
               IAUG = 1
            ELSE IF (BASNAM(1:2) .EQ. 't-') THEN
               IAUG = 2
            ELSE IF (BASNAM(1:2) .EQ. 'q-') THEN
               IAUG = 3

!radovan: added some nonstandard augmentation prefixes
!         i need them to test high-order response
!         use at own risk
            else if (basnam(1:2) .eq. '5-') then
               iaug = 4
            else if (basnam(1:2) .eq. '6-') then
               iaug = 5
            else if (basnam(1:2) .eq. '7-') then
               iaug = 6

            ELSE
               WRITE (LUPRI,'(/A2,A)') BASNAM(1:2),' is an unknown'//
     &              'augmentation level'
               CALL QUIT('Illegal augmentation level in HERBAS')
            END IF

#define ALLOW_NONSTANDARD_AUGMENTATION
#ifdef ALLOW_NONSTANDARD_AUGMENTATION
!        else if (basnam(3:15) == 'aug-Turbomole') then
!radovan: this does not work because basis set is max 20 characters long
         else if (basnam(3:13) == 'a-Turbomole') then
           bastmp(1:76) = basnam(5:80)
           if      (basnam(1:2) == 's-') then
             iaug = 1
           else if (basnam(1:2) == 'd-') then
             iaug = 2
           else if (basnam(1:2) == 't-') then
             iaug = 3
           else if (basnam(1:2) == 'q-') then
             iaug = 4
           else
             call quit('illegal augmentation in herbas')
           end if
         else if (basnam(3:11) == 'aug-dyall') then
           bastmp(1:74) = basnam(7:80)
           if      (basnam(1:2) == 's-') then
             iaug = 1
           else if (basnam(1:2) == 'd-') then
             iaug = 2
           else if (basnam(1:2) == 't-') then
             iaug = 3
           else if (basnam(1:2) == 'q-') then
             iaug = 4
           else
             call quit('illegal augmentation in herbas')
           end if
#endif /* ifdef ALLOW_NONSTANDARD_AUGMENTATION */

C for backward compatibility the old notation is retained
         ELSE IF (BASNAM(2:9) .EQ. 'aug-cc-p') THEN
            BASTMP(1:79) = BASNAM(2:80)
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
         ELSE IF (BASNAM(1:3) .EQ. 'ano') THEN
            ANO = .TRUE.
            BASTMP(1:5) = BASNAM(1:5)
            ISTART = 6
            !ISTART = 1
         ELSE IF (BASNAM(1:7) .EQ. 'ANO-RCC') THEN
            ANO = .TRUE.
            BASTMP(1:7) = BASNAM(1:7)
            ISTART = 8
         ELSE IF (BASNAM(1:6) .EQ. 'sadlej') THEN
            SADLEJ = .TRUE.
            BASTMP = 'sadlej'
         ELSE IF (BASNAM(1:4) .EQ. 'NQvD') THEN
            NQVD = .TRUE.
            BASTMP = 'NQvD'
         ELSE
            ISTART = INDEX(BASNAM,' ')           
            BASTMP = BASNAM(1:ISTART - 1)
         END IF
      END IF

      iaug_save = iaug
C
C     Determine if there are any user added basis functions/polarization 
C     functions
C
      IF (ISET .EQ. 1) THEN
         IPOLST = INDEX(BASNAM,'Pol')
      ELSE
         POLFUN = .FALSE.
         IPOLST = 0
      END IF
      IF (IPOLST .GT. 0) THEN
         POLFUN = .TRUE.
         IPOLST = IPOLST + 3
      END IF

      IF (IPREAD.GE.10) THEN
        CALL HEADER('Output from BASLIB',-1)
        WRITE(LUPRI,*) 'ISTART=',ISTART
        WRITE(LUPRI,'(2A)') ' BASTMP=',BASTMP
        WRITE(LUPRI,*) 'IPOLST=',IPOLST
        WRITE(LUPRI,*) 'POLFUN=',POLFUN
        WRITE(LUPRI,*) 'ANO   =',ANO
        WRITE(LUPRI,*) 'SADLEJ=',SADLEJ
        WRITE(LUPRI,*) 'NQVD  =',NQVD
      ENDIF
C      
C Finds the right element in the file. 
C
      IF (NQVD) THEN
         CALL FINNQD(BASNAM,INTQ,NBLOCK,ALPHA,CPRIMU,CPRIM,NUC,NRC,SEG,
     &               KPRIM,KAOVEC,DSM,POLFUN,IPOLST,LUBAS)
      ELSE
         CALL FINDEL(BASTMP, INTQ,LUBAS)
C
C Finds the number of primitives and columns of contraction coeffecients
C
         CALL FINPOS ( NEWEL, INTEXP, INTORB, INTISG, LUBAS,IPREAD)
C     
 10      CONTINUE
C
         IF (TEST) THEN
            INTORB = -INTORB
            IF (ANO .OR. SADLEJ) THEN
               WRITE (LUPRI,'(/A)') ' FATAL ERROR: '//
     &            '.UNCONT not implemented for ano and sadlej yet.'
               CALL QUIT(
     &            '.UNCONT not implemented for ano and sadlej yet')
Chj-aug99: MAERKE find out how to implement UNCONT for ano or sadlej !
            END IF
         END IF
C
         IF ( .NOT. NEWEL) THEN

            READNU_CALLED=.FALSE.
            
            IQM = IQM + 1
            JCO(IQM) = 1
  11        CONTINUE
CMI           ... to have INTORB fixed for repeating calls of READNU
            NBLOCK = NBLOCK + 1
            IF (NBLOCK .GT. KAOVEC) THEN
               CALL QUIT('Too many AO blocks'//
     &            ' for this program version; increase WRKMEM')
            END IF
C
C           We should not make blocks of primitives too large in case of uncontracted sets
C           Compute maximum size (should of course always be one, hence the MAX in the 
C           determination step).
C
            KHK_SC = (IQM+1)*(IQM+2)/2
            MAXPRIM_IN_BLOCK = MAX(1,MAX_IN_BLOCK/KHK_SC)
C
            IF (INTORB.LE.0.AND.INTEXP.GT.MAXPRIM_IN_BLOCK) THEN
               IF (ANO .OR. SADLEJ) THEN
                  SPLIT = .FALSE.
               ELSE
                  INTEXP_2 = INTEXP - MAXPRIM_IN_BLOCK
                  INTEXP = MAXPRIM_IN_BLOCK
                  JCO(IQM) = JCO(IQM) + 1
                  SPLIT = .TRUE.
               ENDIF
            ELSE
               SPLIT = .FALSE.
            ENDIF
            
C     
C Setting ALPHA, CPRIMU AND CPRIM to zero. 
C
            CALL DZERO(ALPHA(1, NBLOCK), KPRIM)
            CALL DZERO(CPRIMU(1,1, NBLOCK), KPRIM*KPRIM)
            CALL DZERO(CPRIM (1, 1, NBLOCK), KPRIM*KPRIM)
C
C Reading the primitives and contraction-coeffecients from file in READNU
C
            IF (ANO .OR. SADLEJ) THEN
!              write(lupri,*) 'NBLOCK,ISTART=',NBLOCK,ISTART
               CALL REAANO(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU,ALPHA,
     +                     KPRIM,ANO,SADLEJ,BASNAM,ISTART,POLFUN,IPOLST,
     +                     LUBAS,IPREAD)
               IF (INTORB .EQ. 0) THEN
                  NBLOCK = NBLOCK - 1
                  CALL FINPOS ( NEWEL, INTEXP, INTORB,
     &                           INTISG, LUBAS, IPREAD)
                  GOTO 10
               END IF
            ELSE

               IF (READNU_CALLED) THEN
C                INTEXP = INTEXP_SAVE
                 INTORB = INTORB_SAVE
                 IF (IPREAD.GE.7) THEN
                 write(lupri,*) 'restored INTORB=',INTORB
                 ENDIF
               END IF

               if (split) then
!                 block has been split and this is not the last part
!                 augment only last part
!                 do not augment this part
                  iaug = 0
               else
                  iaug = iaug_save
               end if

               INTEXP_SAVE = INTEXP
               INTORB_SAVE = INTORB
               CALL READNU(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU,ALPHA,
     +                KPRIM,IAUG,DSM,BASNAM,POLFUN,IPOLST,LUBAS,IPREAD)

               READNU_CALLED = .TRUE.

            END IF
C
C Checking if the matrix, with the contraction-coeffecients, is segmented 
C
            SEGEJ  = .TRUE.
            CALL SEGORB(SEGEJ,INTEXP,INTORB,CPRIMU(1,1,NBLOCK),KPRIM,
     +                  DSM)
            NUC(NBLOCK) = INTEXP
            NRC(NBLOCK) = INTORB
            SEG(NBLOCK) = SEGEJ
            ISGEN(NBLOCK) = INTISG
            
C
C           Reorder primitive orbitals
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUC(NBLOCK),
     &                  NRC(NBLOCK),SEG(NBLOCK),KPRIM,DSM)

!#define PRINT_GAUSSIAN_INP
#ifdef PRINT_GAUSSIAN_INP
            CALL print_gaussian_inp(icount, iqm, ALPHA(1,NBLOCK),
     &                  CPRIMU(1,1,NBLOCK),NUC(NBLOCK),
     &                  NRC(NBLOCK),SEG(NBLOCK),KPRIM,DSM)
#endif
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
             CALL NRMORB(IQM, NRC(NBLOCK), NUC(NBLOCK), 
     +                   ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                   CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
C
C Go back to read additional exponents in case of block splitting
C
            IF (SPLIT) THEN
               INTORB = 0
               INTEXP = INTEXP_2
               GOTO 11
            ENDIF
C
C Define the type of the next relevant line in file
C
            CALL FINPOS ( NEWEL, INTEXP, INTORB, INTISG, LUBAS, IPREAD)
            GOTO 10
         END IF
      END IF
C
      IF (POLFUN) THEN
 32      CONTINUE
         INTEXP = 0
         INTORB = 0
         IQM = IQM + 1
         NBLOCK = NBLOCK + 1
         FOUND = .FALSE.
         ISTART = IPOLST
 33      CONTINUE
         CALL CHBLANK(BLANK,BASNAM(ISTART:))
         IF (.NOT. BLANK) THEN
            CALL FREFRM(BASNAM,ISTART,IQUANT,DUMMY,'INT')
            CALL FREFRM(BASNAM,ISTART,IDUMMY,EXPON,'REA')
            IF (IQUANT .EQ. IQM) THEN
               FOUND  = .TRUE.
               INTEXP = INTEXP + 1
               INTORB = INTORB + 1
               ALPHA(INTEXP,NBLOCK) = EXPON
               CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0
            END IF
            GOTO 33
         END IF
         IF (FOUND) THEN
            SEGEJ  = .TRUE.
            CALL SEGORB( SEGEJ, INTEXP, INTORB, CPRIMU(1,1,NBLOCK),
     +                   KPRIM, DSM)
            NUC(NBLOCK) = INTEXP
            NRC(NBLOCK) = INTORB
            SEG(NBLOCK) = SEGEJ
C
C           Reorder primitive orbitals
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUC(NBLOCK),
     &                  NRC(NBLOCK),SEG(NBLOCK),KPRIM,DSM)
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
            CALL NRMORB(IQM, NRC(NBLOCK), NUC(NBLOCK), 
     +                  ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                  CPRIMU(1,1, NBLOCK), KPRIM, NBLOCK)
            GOTO 32
         ELSE
            NBLOCK = NBLOCK - 1
         END IF
      END IF
C
C Find reference information for the basis set. 
C
      CALL FINREF (BASREF,LUBAS,BASSAV)
C
C Closing file
C

      BASNAM = BASSAV
      CLOSE (LUBAS, STATUS = 'KEEP')
      CALL QEXIT('BASLIB')
      RETURN
C
      END
C
C*************************************************************
C* This is the subroutine that takes care of the reading of  *
C* the primitives and the contraction coefficients from the  *
C* file (for a given format). The variables that are         *
C* transferred out are:                                      *
C*     CPRIMU - The matrix that the contraction coefficients *
C*              are put in (not normalized).                 *
C*     ALPHA - The matrix where the primitives are put in    *
C*     (Both of the variables are given for a given (quantum *
C*     number + 1))                                          *
C*************************************************************
C
C/* Deck readnu */
      SUBROUTINE READNU(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU, ALPHA, 
     +            KPRIM,IAUG,DSM,BASNAM,POLFUN,IPOLST,LUBAS,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
C Reading the exponents
C
      DIMENSION CPRIMU( KPRIM, KPRIM, KAOVEC), ALPHA(KPRIM, KAOVEC)
C*************************************************************
C* Variable declarations:                                    *
C*     SIGN - Character that gives what kind of sentence it  *
C*            is.                                            *
C*     STRING - A character string that helps bulletproofing *
C*             the subroutine                                *
C*     BLANK - A logical variable that gives whether the line*
C*             is blank or not.                              *
C*************************************************************
      CHARACTER SIGN
      CHARACTER STRING*88, BASNAM*(*)
      LOGICAL BLANK, POLFUN
C
C CONEOR is the maximum number of columns of contraction coefficients on 
C the 1. line in file. If INTORB is .GT. CONEOR then the contraction
C coefficients will continue on the next line. CONORB is the max. number
C of contraction coefficients on the next lines.
C
      INTEGER CONEOR, CONORB
      PARAMETER(CONEOR = 6, CONORB = 7)
      IF (IPRINT.GE.5) THEN
       CALL HEADER('Output from READNU',-1)
       write(lupri,*) 'entering READNU, IPREAD=',IPRINT
       write(lupri,*) 'NBLOCK,INTORB:',NBLOCK,INTORB
       write(lupri,*) 'DSM=',DSM
      ENDIF

C
C Initializing a counting variable that must be 0 before loop
C
      J = 0
C
 20   CONTINUE
      IF ( J .LT. INTEXP) THEN
C
C Reading the primitive and contracted coefficients
C
         READ(LUBAS, '(A88)', IOSTAT = IOERR, ERR = 2000)
     +        STRING

         IF (IPRINT.GE.7)
     &      WRITE(LUPRI,'(A,A)') '  READNU: read STRING=',STRING

         READ (STRING, '(A1)') SIGN

         IF (IPRINT.GE.7)
     &      WRITE(LUPRI,'(A,A)') '  READNU: read SIGN=',SIGN

         IF (SIGN .EQ. ' ') THEN
C
C Calling a subroutine that checks whether the line is a blank one. If
C STRING is a blank line BLANK will return as .TRUE.
C
            CALL CHBLANK(BLANK, STRING)
            IF (.NOT. BLANK) THEN
C
C We have found a line with a primitive and coefficients             
C
               J = J + 1
C
C CONLIN returns the number of lines we have the contraction
C coefficients written on.
C
               JNTORB = ABS(INTORB)

               CALL CONLIN(JNTORB, NUMLIN, IPRINT)

               IF (INTORB .EQ. 0) THEN
C                 uncontracted basis set: allow unformatted read                                                                                      
                  READ (STRING, *,IOSTAT=IOS) ALPHA(J,NBLOCK)
                  CPRIMU(J,J,NBLOCK) = 1.0D0

                  IF (IPRINT.GE.5) write(lupri,*)
     &              'READNU read 1 - NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)',
     &              NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)
                  GOTO 20
               ELSEIF (INTORB .LT. 0) THEN
C                 forced uncontracted with .UNCONT   
                  READ (STRING, '(F16.9)',IOSTAT=IOS) ALPHA(J,NBLOCK)
                  IF (IOS.NE.0) CALL QUIT('Error in reading alpha')

C                 skip continuation lines for contraction coefficients
                  DO I = 2,NUMLIN
                     READ (LUBAS,*)
                  END DO
                   
                  CPRIMU(J,J,NBLOCK) = 1.0D0

                  IF (IPRINT.GE.5) write(lupri,*)
     &              'READNU read 1 - NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)',
     &              NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)

                  GOTO 20
               ELSE
C
C     Getting the format for the read statement right.
C
                  IF (NUMLIN .EQ. 1) THEN
                     KNTORB = INTORB
                  ELSE
                     KNTORB = CONEOR
                  END IF
C
C     Reading the first line with exponents and contraction coefficients
C
                  READ (STRING, '(F16.9, 6F12.9)',IOSTAT=IOS)
     +                 ALPHA(J, NBLOCK), (CPRIMU(J,I,NBLOCK),
     +                 I = 1, KNTORB)
                  IF (IOS.NE.0) THEN
                    WRITE(LUPRI,*) 'KNTORB=',KNTORB
                    WRITE(LUPRI,'(3A)') ' STRING = "',STRING,'"'
                    WRITE(LUPRI,*) 'J,NBLOCK:',J,NBLOCK
                    WRITE(LUPRI,*) 'ALPHA(J, NBLOCK)=',ALPHA(J, NBLOCK)
                    WRITE(LUPRI,*) 'CPRIMU(J,I,NBLOCK), I = 1, KNTORB:',
     &              ( CPRIMU(J,I,NBLOCK),I = 1, KNTORB )
                    CALL QUIT('Error in reading basis set...')
                  ENDIF

                  IF (IPRINT.GE.6) THEN
                    write(lupri,*)
     &              'READNU read 2 - J,NBLOCK,ALPHA(J,NBLOCK)',
     &              J,NBLOCK,ALPHA(J,NBLOCK)
                    write(lupri,*)
     &              'READNU read 2 - CPRIMU:',
     &              (CPRIMU(J,I,NBLOCK),I=1,KNTORB)
                  ENDIF
C   
C     If there are more lines with contraction coefficients, 
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
C     contraction coefficients.
C
                     READ(LUBAS, '(A88)', IOSTAT = IOERR, ERR = 2000)
     +                    STRING
                     READ (STRING, '(A1)') SIGN
                     IF (SIGN .EQ. ' ') THEN
                        CALL CHBLANK(BLANK, STRING)
                        IF (.NOT. BLANK) THEN
C   
C     We now have a line with contraction coefficients.
C   
                           READ (STRING,'(F16.9,6F12.9)',IOSTAT=IOS) 
     &                          (CPRIMU(J,M,NBLOCK),
     +                          M = CONEOR + (I-2)*CONORB +1, KNTORB)
                           IF (IOS.NE.0) THEN
                             CALL QUIT('Error in reading sjjd393jdj')
                           ENDIF
                           IF (IPRINT.GE.6) THEN
                             WRITE(LUPRI,*) 
     &                       J,NLBLOCK,' CPRIMU(J,M,NBLOCK):'
                             WRITE(LUPRI,*) ( CPRIMU(J,M,NBLOCK),
     &                         M = CONEOR + (I-2)*CONORB +1, KNTORB )
                           ENDIF
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
         ELSE    ! IF (SIGN .EQ. ' ') THEN
C
C Found a line with nothing sensible in it
C
            GOTO 20
         END IF
      ELSE   ! IF ( J .LT. INTEXP) THEN
C
C        For uncontracted basis sets, update INTORB.
C
         IF ( INTORB .LE. 0 ) INTORB = INTEXP
C
C   
C Read all of the primitives, and the coefficients
C If this is some kind of augmented basis set, we augment it here, kr-96
C
         IF (IAUG .GT. 0) THEN
            do kaug = 1, iaug

!             find lowest nonzero alpha
!             note that they do not have to be monotonically decreasing!
              alpha_lowest = maxval(alpha(1:intexp, nblock))
              do iprim = 1, intexp
                if (
     &              (alpha(iprim, nblock) < alpha_lowest)
     &              .and.
     &              (alpha(iprim, nblock) > 0.0d0)
     &             ) then
                  alpha_lowest = alpha(iprim, nblock)
                end if
              end do

!             find second lowest nonzero alpha
              alpha_second_lowest = maxval(alpha(1:intexp, nblock))
              do iprim = 1, intexp
                if (
     &              (alpha(iprim, nblock) < alpha_second_lowest)
     &              .and.
     &              (alpha(iprim, nblock) > alpha_lowest)
     &             ) then
                  alpha_second_lowest = alpha(iprim, nblock)
                end if
              end do

              if (dabs(alpha_second_lowest 
     &                    - alpha_lowest) < 1.0d-20) then
#ifdef ALLOW_NONSTANDARD_AUGMENTATION
!               if there is only one primitive in the original set
!               then alpha_second_lowest equals alpha_lowest
!               in this case alpha_second_lowest is multiplied by 3.5
!               (factor 3.5 is my choice, no really good reason for it)
                alpha_second_lowest = 3.5d0*alpha_lowest
#else
                call quit('error in automatic augmentation')
#endif
              end if

              intexp = intexp + 1
              intorb = intorb + 1

              alpha(intexp, nblock) 
     &            = alpha_lowest*alpha_lowest/alpha_second_lowest
              cprimu(intexp, intorb, nblock) = 1.0d0
            end do

             IF (IPRINT.GE.5) THEN
              write(lupri,*) 'READNU:  INTEXP,NBLOCK,ALPHA',
     &         INTEXP,NBLOCK,ALPHA(INTEXP,NBLOCK)
              write(LUPRI,*) 'READNU: CPRIMU(INTEXP,INTORB,NBLOCK) =',
     &         CPRIMU(INTEXP,INTORB,NBLOCK)
             ENDIF

         END IF

         IF (POLFUN) THEN
            ISTART = IPOLST
 33         CONTINUE
            CALL CHBLANK(BLANK,BASNAM(ISTART:))
            IF (.NOT. BLANK) THEN
               CALL FREFRM(BASNAM,ISTART,IQUANT,DUMMY,'INT')
               CALL FREFRM(BASNAM,ISTART,IDUMMY,EXPON,'REA')
               IF (IQUANT .EQ. NBLOCK) THEN
                  INTEXP = INTEXP + 1
                  INTORB = INTORB + 1
                  ALPHA(INTEXP,NBLOCK) = EXPON
                  CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0

              IF (IPRINT.GE.5) THEN
             write(lupri,*) 'INTEXP,NBLOCK,ALPHA',
     &         INTEXP,NBLOCK,ALPHA(INTEXP,NBLOCK)
              ENDIF

               END IF
               GOTO 33
            END IF
         END IF

         IF (IPRINT.GE.4) THEN
            write(lupri,*) 'READNU...NBLOCK,J=',NBLOCK,J
            write(lupri,*) (ALPHA(IX,NBLOCK),IX=1,J)
         ENDIF

      END IF
C
      RETURN
C
 2000 CONTINUE
      WRITE (LUPRI,'(/A,I5)') 'Error in reading from basis file.'//
     &     ' IOSTAT =',IOERR
      END
C/* Deck reaano */
      SUBROUTINE REAANO(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU, ALPHA, 
     +                  KPRIM,ANO,SADLEJ,BASNAM,ISTART,POLFUN,IPOLST,
     +                  LUBAS,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
C
C Reading the exponents. Format from MOLCAS files
C
      LOGICAL ANO, SADLEJ, POLFUN, BLANK
      DIMENSION CPRIMU( KPRIM, KPRIM, KAOVEC), ALPHA(KPRIM, KAOVEC)
      CHARACTER SIGN
      CHARACTER BASNAM*80
C
C Determine number of contracted functions to be read
C
      IF (IPRINT.GE.8) THEN
        CALL HEADER('Output from REAANO',-1)
        wRITE(LUPRI,'(A,A)') 'BASNAM=',BASNAM
        wRITE(LUPRI,*) 'ISTART=',ISTART
        wRITE(LUPRI,*) 'INTEXP,INTORB=',INTEXP,INTORB
        WRITE(LUPRI,*) 'NBLOCK=',NBLOCK
        WRITE(LUPRI,*) ' KPRIM=',KPRIM
      ENDIF

!
!      ISTART_ANO=10
!     IF (ANO) CALL FREFRM(BASNAM,ISTART,INTORB,DUMMY,'INT')
!     IF (ANO) CALL FREFRM(BASNAM,ISTART_ANO,INTORB,DUMMY,'INT')

      IF (IPRINT.GE.8) THEN
        WRITE(LUPRI,'(/,8x,a)') 'REAANO - After FREFRM:'
        WRITE(LUPRI,'(A,A)') 'BASNAM=',BASNAM
        WRITE(LUPRI,*) 'ISTART=',ISTART
        wRITE(LUPRI,*) 'INTEXP,INTORB=',INTEXP,INTORB
      ENDIF

C
 20   CONTINUE
C
C Reading the primitive and contracted coeffecients
C
      MAXEXP = 8
      NUMEXP = (INTEXP - 1)/MAXEXP
      J = 1
      IF (INTORB .GT. 0) THEN
         DO I = 0, NUMEXP
            IF (I .EQ. NUMEXP) THEN
               IF (INTEXP .EQ. MAXEXP) THEN
                  READ (LUBAS,*,IOSTAT=IOS) (ALPHA(K,NBLOCK), K = J, 
     &                                              J - 1 + MAXEXP)
                  IF (IOS.NE.0) THEN
                    CALL QUIT('error in reading skdsjdj')
                  ENDIF
               ELSE
                  READ (LUBAS,*,IOSTAT=IOS) (ALPHA(K,NBLOCK), K = J, 
     &                                    J - 1 + MOD(INTEXP,MAXEXP))
                  IF (IOS.NE.0) THEN
                    CALL QUIT('error in reading s233jdj')
                  ENDIF
               END IF
            ELSE
               READ (LUBAS,*,IOSTAT=IOS) 
     &         (ALPHA(K,NBLOCK), K = J, J + MAXEXP - 1)
               IF (IOS.NE.0) THEN
                  CALL QUIT('error in reading 2s2adjj')
               ENDIF
            END IF
            J = J + MAXEXP
         END DO
         DO I = 1, INTEXP
            READ (LUBAS,*,IOSTAT=IOS) 
     &       (CPRIMU(I,J,NBLOCK), J = 1, INTORB)
            IF (IOS.NE.0) THEN
              CALL QUIT('error in reading assjj')
            ENDIF
         END DO
      ELSE
         DO I = 0, NUMEXP
            IF (I .EQ. NUMEXP) THEN
               READ (LUBAS,*)
            ELSE
               READ (LUBAS,*)
            END IF
            J = J + MAXEXP
         END DO
         DO I = 1, INTEXP
            READ (LUBAS,*)
         END DO
      END IF
      IF (POLFUN) THEN
         IREAD = IPOLST
 33      CONTINUE
         CALL CHBLANK(BLANK,BASNAM(IREAD:))
         IF (.NOT. BLANK) THEN
            CALL FREFRM(BASNAM,IREAD,IQUANT,DUMMY,'INT')
            CALL FREFRM(BASNAM,IREAD,IDUMMY,EXPON,'REA')
            IF (IQUANT .EQ. NBLOCK) THEN
               INTEXP = INTEXP + 1
               INTORB = INTORB + 1
               ALPHA(INTEXP,NBLOCK) = EXPON
               CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0
            END IF
            GOTO 33
         END IF
      END IF

      IF (IPRINT.GE.10) THEN
        WRITE(LUPRI,*) 'REAANO-read coefficients...'
        WRITE(LUPRI,*) 'NBLOCK,INTEXP:',NBLOCK,INTEXP
        WRITE(LUPRI,*) 'INTORB:',INTORB
        DO I = 1, INTEXP
          WRITE(LUPRI,*)
     &   (CPRIMU(I,J,NBLOCK), J = 1, INTORB)
        ENDDO
      ENDIF

      RETURN
      END
C
C*************************************************************
C* CHBLANK is a subroutine that checks whether a line is a blank*
C* line. It returns a logical variable if it is a blank one. *
C* The in/out-going variables are:                           *
C*     STRING - A character string (the string checked by the*
C*              subroutine).                                 *
C*     BLANK - The logical variable. It will be .TRUE. if the*
C*             line is a blank one.                          *
C*************************************************************
      SUBROUTINE CHBLANK(BLANK,STRING)
C
      CHARACTER*(*) STRING
      LOGICAL BLANK
C
      BLANK = .TRUE.
C
      DO 100 J = 1, LEN(STRING)
         BLANK = BLANK .AND. (STRING(J:J) .EQ. ' ')
 100  CONTINUE
      RETURN
      END
C
C*************************************************************
C* CONLIN is a subroutine that finds out how many lines the  *
C* contraction coefficients are written on.  There are 6     *
C* cc (contraction coefficients) on the 1. line, and then    *
C* there are seven on the rest of the lines. In out going    *
C* variables are :                                           *
C*      INTORB - Total number of columns of contraction      *
C*               coefficients.                               *
C*      NUMLIN - Returns the number of lines the contraction *
C*               coefficients are written on.                *
C*************************************************************
C
      SUBROUTINE CONLIN(INTORB, NUMLIN, IPRINT)
#include "implicit.h"
#include "priunit.h"
C      
C This is the parameters for the number of contraction coefficients
C on 1. and 2. line.
C
      INTEGER    CONORB
      PARAMETER (CONORB = 7)
C
C This finds out how many lines we have, and puts it into numlin.
C
      NUMLIN = INTORB/CONORB + 1
!     truncation of: NUMLIN = ( (INTORB + 1) - 1)/CONORB + 1
!     as we have INTORB + 1 numbers: exponent + INTORB contraction
!     coefficients.

      IF (IPRINT.GE.10) THEN
       write(lupri,'(/A)') '  *** CONLIN output ***'
       write(lupri,*) 'input INTORB =',INTORB
       write(lupri,*) 'output NUMLIN=',NUMLIN
      ENDIF

      RETURN
      END
C********************************************************
C* This are the subroutines that  reads  the primitives *
C* and the contracted coefficients from file.           *
C*********************************************************
C* This is the subroutine that searches through the file *
C* and finds the element in question.                    *
C*********************************************************
C/* Deck findel */
      SUBROUTINE FINDEL(BASNAM, INTQ,LUBAS)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basis file               *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*     EXST - Logical variable that helps inquire if     *
C*            there exists a file with that name.        *
C*********************************************************
      CHARACTER*(*) BASNAM
      CHARACTER*200 STRING  ! long because contains the full directory path
      CHARACTER SIGN
      LOGICAL EXST, EXST2, USE_BASNAM
#include "inftap.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
C     The environment variable is a : separated string, e.g.:
C     /programs/Dirac/basis:/programs/Dalton/basis:.
C
C     We loop over these directories in search for the basis set and
c     element given.
C
C     Depending on the machine, the string has trailing blanks
C     (ASCII char=32) or null elements (ASCII=0).
      ISTART = 1
      IBLANK=INDEX(BASDIR,CHAR(32))
      INULL =INDEX(BASDIR,CHAR(0))
      ISTOP=MAX(IBLANK,INULL)
      IF(ISTOP.EQ.0) THEN
C       If string is empty add current directory 
        BASDIR(1:1) = '.'
        ISTOP = 1
      ELSE
        ISTOP=ISTOP-1
      ENDIF
C     
 10   CONTINUE
C     
C     Find next colon in BASDIR
C     
      ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:ISTOP), ':' )
C     
C     In there are no colons, we are the last directory in BASDIR
C     
      IF ( ICOLON .EQ. ISTART - 1 ) ICOLON = ISTOP + 1
      IEND = ICOLON - 1
      STRING = BASDIR(ISTART:IEND)//' '
C
C     Add trailing slash if not present.
C
      IEND = INDEX ( STRING(1:), ' ') - 1
!     stefan: make check with idb indicates a possible access penalty if
!     BASDIR is not set...
      if(iend.lt.1)then
        write(lupri,'(/a)') '  *** error in FINDEL: environment '//
     &        'variable BASDIR is not set'
        call quit('*** error in FINDEL: BASDIR is not set') 
      end if
      IF ( STRING(IEND:IEND) .NE. '/' ) THEN
         IEND = IEND + 1
         STRING(IEND:IEND) = '/'
      END IF
      STRING = STRING(1:IEND)//BASNAM
      IEND = INDEX(STRING(1:),' ') - 1
Cluuk Activate when debugging basis set search (should not be on by default)
CLuuk WRITE(LUPRI,'(A/,A1,A,A1)')
CLuuk&     ' Searching basis : trying file: ','"',STRING(1:IEND),'"'
C
C     Inquire if "/path/basis-set" exists.
C
      INQUIRE (FILE = STRING, EXIST = EXST)
CMI  ... on certain systems (like PC-linux,OpenMPI-gfortran) string
C     beginning with '/' is NOT the proper file name ! therefore
C     continue in inquiring:
      EXST2=EXST
      USE_BASNAM=.FALSE.
      IF (.NOT.EXST2) THEN
        INQUIRE (FILE = BASNAM, EXIST = EXST)
        IF (EXST) THEN  
          WRITE(LUPRI,"(/2A/2A)")
     &     "Info: used proper basis set name: ",BASNAM,
     &     "Info: instead of                : ",STRING
          USE_BASNAM=.TRUE.
        ENDIF
      ENDIF
        
      IF (EXST) THEN
        IF (.NOT.USE_BASNAM) THEN
         OPEN(LUBAS, FILE=STRING)
        ELSE
         OPEN(LUBAS, FILE=BASNAM)
        ENDIF
C     
C
C        Searching the file for the element.
C
 20      CONTINUE
         READ(LUBAS,'(A88)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)',IOSTAT=IOS) SIGN
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
            READ (STRING, 110,IOSTAT=IOS) SIGN, NUCEL
            IF (IOS.NE.0) THEN
               CALL QUIT('FINDEL: Error in reading NUCEL')
            ENDIF
 110        FORMAT (BN, A1, I4)
            IF (INTQ .EQ. NUCEL) THEN
               RETURN
            ELSE
               GOTO 20
            END IF
         ELSE 
            GOTO 20
         END IF   
      END IF
      ISTART = ICOLON + 1
      IF(ISTART.GT.ISTOP) GOTO 300
      GOTO 10
 300  CONTINUE
      LSTR = 1
      DO J = LEN(BASNAM), 1, -1
         IF ( BASNAM(J:J) .NE. ' ' ) THEN
            LSTR = J
            GOTO 301
         END IF
      END DO
 301  CONTINUE
      WRITE (LUPRI,'(/,3A)')
     +     'Basis "',BASNAM(1:LSTR),'" doesn''t exist'
      WRITE(LUPRI,'(2X,A,A)') 'basis name, STRING=',STRING
      WRITE(LUPRI,'(2X,A,A)') 'basis name, BASNAM=',BASNAM
      IF (BASNAM(1:6) .EQ. 'Huckel') THEN
         WRITE (LUPRI,'(A)')
     +        'NB! Basis Huckel must be available to use Huckel guess!'
      END IF
      CALL QUIT('Non-existing basis set in HERBAS')
C
C Error messages 
C
 2000 CONTINUE
      WRITE (LUPRI,'(/2A/A,I5)')
     &   'FINDEL: Error when reading from basis file ',STRING(1:IEND),
     &   ' IOSTAT =',IOERR
      CALL QUIT('I/O error in FINDEL in HERBAS')
 200  CONTINUE
      WRITE (LUPRI,'(/I3,2A)') INTQ,
     &     ' is an unsupported element for basis ',BASNAM
      CALL QUIT('Unsupported element in HERBAS')
      END 
C/* Deck finnqd */
      SUBROUTINE FINNQD(BASNAM, INTQ,NBLOCK,ALPHA,CPRIMU,CPRIM,NUC,
     &     NRC,SEG,KPRIM,KAOVEC,DSM,POLFUN,IPOLST,
     $     LUBAS)
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basis file               *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*     EXST - Logical variable that helps inquire if     *
C*            there exists a file with that name.        *
C*********************************************************
C
C
      CHARACTER*(*) BASNAM
      CHARACTER*140 STRING
      CHARACTER SIGN
      LOGICAL EXST, POLFUN, SEG, SEGEJ, BLANK
      DIMENSION ALPHA(KPRIM,KAOVEC), CPRIMU(KPRIM,KPRIM,KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC), NUC(KAOVEC), NRC(KAOVEC),
     &          SEG(KAOVEC)
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
C     Determine primitive basis
C
      ISTART = 5
      CALL FREFRM(BASNAM,ISTART,ISPRIM,DUMMY,'INT')
      IF (INTQ .GT. 4) THEN
         CALL FREFRM(BASNAM,ISTART,IPPRIM,DUMMY,'INT')
      END IF
      CALL FREFRM(BASNAM,ISTART,ISCONT,DUMMY,'INT')
      IF (INTQ .GT. 4) THEN
         CALL FREFRM(BASNAM,ISTART,IPCONT,DUMMY,'INT')
      END IF
      STRING = BASDIR(1:LENBAS)//'NQvD'
      INQUIRE (FILE = STRING, EXIST = EXST)
      IF (EXST) THEN
         OPEN(LUBAS, FILE = STRING)
C
C
C Searching the file for the element.
C
 20      CONTINUE
         READ(LUBAS,'(A88)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)', IOSTAT = IOS) SIGN
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
            ISTART = 2
            CALL FREFRM(STRING,ISTART,NUCEL,DUMMY,'INT')
            CALL FREFRM(STRING,ISTART,IS,DUMMY,'INT')
            IF (NUCEL .EQ. INTQ .AND. IS .EQ. ISPRIM) THEN
               IF (NUCEL .GT. 4) THEN
                  CALL FREFRM(STRING,ISTART,IP,DUMMY,'INT')
                  IF (IP .NE. IPPRIM) GO TO 20
               END IF
               NBLOCK = NBLOCK + 1
               IF (NBLOCK .GT. KAOVEC) THEN
                  CALL QUIT('Too many AO blocks'//
     &               ' for this program version; increase WRKMEM')
            END IF
C
C Setting ALPHA, CPRIMU AND CPRIM to zero. 
C
               CALL DZERO(ALPHA(1, NBLOCK), KPRIM)
               CALL DZERO(CPRIMU(1,1, NBLOCK), KPRIM*KPRIM)
               CALL DZERO(CPRIM (1,1, NBLOCK), KPRIM*KPRIM)
               IF (NUCEL .GT. 4) THEN
                  CALL DZERO(ALPHA(1, NBLOCK + 1), KPRIM)
                  CALL DZERO(CPRIMU(1,1, NBLOCK + 1), KPRIM*KPRIM)
                  CALL DZERO(CPRIM (1,1, NBLOCK + 1), KPRIM*KPRIM)
               END IF
               IF (NUCEL .LE. 4) THEN
                  ICONT = 1
                  DO ISLOOP = 1, ISPRIM
                     READ (LUBAS,*,IOSTAT=IOS) 
     &                ALPHA(ISLOOP,NBLOCK), CONT
                     IF (IOS.NE.0) THEN
                        CALL QUIT('FINNQD: Error reading ALPHA, CONTS')
                     ENDIF
                     IF (NUCEL .EQ. 1) ALPHA(ISLOOP,NBLOCK) =
     &                    ALPHA(ISLOOP,NBLOCK)*1.44D0
                     IF (ICONT .LT. ISCONT) THEN
                        CPRIMU(ISLOOP,ICONT,NBLOCK) = 1.0D0
                        ICONT = ICONT + 1
                     ELSE
                        CPRIMU(ISLOOP,ICONT,NBLOCK) = CONT
                     END IF
                  END DO
               ELSE
                  ICONTS = 1
                  ICONTP = 1
                  DO IPLOOP = 1, IPPRIM
                     READ (LUBAS,*,IOSTAT=IOS) 
     &                     ALPHA(IPLOOP,NBLOCK), CONTS,
     &                     ALPHA(IPLOOP,NBLOCK + 1), CONTP
                     IF (IOS.NE.0) THEN
                       CALL QUIT('FINNQD Error reading ALPHA,CONTS,...')
                     ENDIF
                     IF (ICONTP .LT. IPCONT) THEN
                        CPRIMU(IPLOOP,ICONTP,NBLOCK + 1) = 1.0D0
                        ICONTP = ICONTP + 1
                     ELSE
                        CPRIMU(IPLOOP,ICONTP,NBLOCK + 1) = CONTP
                     END IF
                     IF (ICONTS .LT. ISCONT) THEN
                        CPRIMU(IPLOOP,ICONTS,NBLOCK) = 1.0D0
                        ICONTS = ICONTS + 1
                     ELSE
                        CPRIMU(IPLOOP,ICONTS,NBLOCK) = CONTS
                     END IF
                  END DO
                  DO ISLOOP = IPLOOP, ISPRIM
                     READ (LUBAS,*,IOSTAT=IOS)
     &                 ALPHA(ISLOOP,NBLOCK), CONT
                     IF (IOS.NE.0) THEN
                        CALL QUIT('FINNQD Error reading ALPHA,CONTS')
                     ENDIF
                     IF (ICONTS .LT. ISCONT) THEN
                        CPRIMU(ISLOOP,ICONTS,NBLOCK) = 1.0D0
                        ICONTS = ICONTS + 1
                     ELSE
                        CPRIMU(ISLOOP,ICONTS,NBLOCK) = CONT
                     END IF
                  END DO
               END IF
               IF (POLFUN) THEN
                  ISTART = IPOLST
 33               CONTINUE
                  CALL CHBLANK(BLANK,BASNAM(ISTART:))
                  IF (.NOT. BLANK) THEN
                     CALL FREFRM(BASNAM,ISTART,IQUANT,DUMMY,'INT')
                     CALL FREFRM(BASNAM,ISTART,IDUMMY,EXPON,'REA')
                     IF (IQUANT .EQ. NBLOCK) THEN
                        ISPRIM = ISPRIM + 1
                        ISCONT = ISCONT + 1
                        ALPHA(ISPRIM,NBLOCK) = EXPON
                        CPRIMU(ISPRIM,ISCONT,NBLOCK) = 1.0D0
                     ELSE IF (IQUANT .EQ. (NBLOCK + 1) 
     &                       .AND. NUCEL .GT. 4) THEN
                        IPPRIM = IPPRIM + 1
                        IPCONT = IPCONT + 1
                        ALPHA(IPPRIM,NBLOCK + 1) = EXPON
                        CPRIMU(IPPRIM,IPCONT,NBLOCK + 1) = 1.0D0
                     END IF
                     GOTO 33
                  END IF
               END IF
C
C Checking if the matrix, with the contraction coefficients, is segmented 
C
               SEGEJ  = .TRUE.
               CALL SEGORB(SEGEJ, ISPRIM, ISCONT, CPRIMU(1,1,NBLOCK), 
     +                     KPRIM, DSM)
               NUC(NBLOCK) = ISPRIM
               NRC(NBLOCK) = ISCONT
               SEG(NBLOCK) = SEGEJ
C
C              Reorder primitive orbitals
C
               CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),
     &                     NUC(NBLOCK),NRC(NBLOCK),SEG(NBLOCK),
     &                     KPRIM,DSM)
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
               CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +                     ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                     CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
               IF (NUCEL .GT. 4) THEN
C
C Checking if the matrix, with the contraction coefficients, is segmented 
C
                  NBLOCK = NBLOCK + 1
                  SEGEJ = .TRUE.
                  CALL SEGORB(SEGEJ, IPPRIM, IPCONT,
     +                        CPRIMU(1,1,NBLOCK),KPRIM, DSM)
                  NUC(NBLOCK) = IPPRIM
                  NRC(NBLOCK) = IPCONT
                  SEG(NBLOCK) = SEGEJ
C
C                 Reorder primitive orbitals
C
                  CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),
     &                        NUC(NBLOCK),NRC(NBLOCK),SEG(NBLOCK),
     &                        KPRIM,DSM)
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
                  CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +                        ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                        CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
               END IF
               RETURN
            ELSE
               GOTO 20
            END IF
         ELSE 
            GOTO 20
         END IF
      ELSE
         LSTR=1
         DO J = 140,1,-1
            IF (STRING(J:J) .NE. ' ') THEN
               LSTR = J
               GO TO 300
            END IF
         END DO
  300    WRITE (LUPRI,'(/,3A)')
     +      'Basis "',STRING(1:LSTR),'" doesn''t exist'
         CALL QUIT('Non-existing basis set in HERBAS')
      END IF
C
C Error messages 
C
 2000 CONTINUE
      WRITE (LUPRI,'(/2A/A,I5)')
     &   'FINNQD: Error when reading from basis file ',
     &   BASDIR(1:LENBAS)//'NQvD',
     &   ' IOSTAT =',IOERR
      CALL QUIT('I/O error in FINNQD in HERBAS')
 200  CONTINUE
      WRITE (LUPRI,'(/I3,2A)') INTQ,
     &     ' is an unsupported element for basis ',BASNAM
      CALL QUIT('Unsupported element in HERBAS')
      END 
C/* Deck finref */
      SUBROUTINE FINREF(BASREF, LUBAS, BASNAM)
C*********************************************************
C* This is the subroutine that searches through the file *
C* and finds reference info if available.                *
C* Based on FINPOS, LV, 14-5-2003                         *
C*********************************************************
#include "implicit.h"
#include "priunit.h"
#ifdef VAR_MPI
#include "infpar.h"
#endif

      CHARACTER*80 LINE,BASREF(10), BASNAM
C*********************************************************
C* Variable declarations:                                *
C*     SIGN - The first sign in a sentence, to find out  *
C*            what kind of sentence it is.               *
C*     STRING - A character string to help bulletproofing*
C*              the subroutine.                          *
C*********************************************************
      CHARACTER SIGN
      CHARACTER*88 STRING
      LOGICAL BLANK
C
      REWIND (LUBAS)
C
C Initialising IREF (number of documentation lines), set to one
C once we find a line containing the string REFERENCE
C
      IREF = 0
C
C     The loops are written f77 style since f90 may still not be
C     available everywhere, better to use while loops here..
C
 20   CONTINUE
      READ (LUBAS, '(A88)', IOSTAT = IOERR, ERR = 1000, END = 200)
     +     STRING
      IF (IOERR .NE. 0) GO TO 200
      READ (STRING, '(A1)',IOSTAT=IOS) SIGN
C
      IF (SIGN .EQ. '$') THEN
         IF (IREF.EQ.0) THEN
C           Check for start of documentation 
            DO I = 2, 80
               IF (STRING(I:I+8).EQ.'REFERENCE') IREF = 1
            ENDDO
         ELSE
C           We have documentation, omit blank lines and store in BASREF
            CALL CHBLANK(BLANK, STRING(2:81))
            IF (.NOT.BLANK) THEN
               BASREF(IREF) = STRING(2:81)
               IREF = IREF + 1
            ENDIF
         END IF
C
C        Maximum number of documentation lines is 10
C
         IF (IREF.GT.10) RETURN
      ELSE 
C        We return if we found a REFERENCE line but no more $ lines.
         IF (IREF.GT.0) RETURN
      END IF
C
      GOTO 20
C
C No documentation when end of file is found
C
 200  CONTINUE
 1000 CONTINUE
      END
C             
C
C*********************************************************
C* This is the subroutine that finds out whether it is a *
C* new block with primitives and contraction             *
C* coefficients, or if it is a new element. It returns   *
C* three variables:                                      *
C*          NEWEL -Logical variable that gives whether it*
C*                 is a new element or not.              *
C*          INEXP, INTORB - Number of primitives and     *
C*                          columns of contraction       *
C*                          coefficients                 *
C*********************************************************
C
C/* Deck finpos */
      SUBROUTINE FINPOS(NEWEL, INTEXP, INTORB, INTISG, LUBAS, IPRINT)
#include "implicit.h"
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     SIGN - The first sign in a sentence, to find out  *
C*            what kind of sentence it is.               *
C*     STRING - A character string to help bulletproofing*
C*              the subroutine.                          *
C*     BLANK - A logical variable to help bulletproofing *
C*             the subroutine.
C*********************************************************
      LOGICAL NEWEL
      CHARACTER SIGN
      CHARACTER*88 STRING
      LOGICAL BLANK
C
C Initializing NEWEL, because one is innocent until proven guilty
C
      NEWEL = .FALSE.
C
 20   CONTINUE
      READ (LUBAS, '(A88)', IOSTAT = IOERR, ERR = 1000, END = 200)
     +     STRING
      READ (STRING, '(A1)',IOSTAT=IOS) SIGN
C
      IF (SIGN .EQ. ' ') THEN
C
C Calling on a subroutine that checks whether the line is a blank one or not.
C If STRING is a blank line, then BLANK will return as .TRUE.
         CALL CHBLANK(BLANK, STRING)
         IF (BLANK) THEN
            GOTO 20
         ELSE
            READ (STRING, '(BN,3I5)',IOSTAT=IOS) INTEXP, INTORB, INTISG
            IF (IOS.NE.0) THEN
              CALL QUIT('Error in reading ... ssjdjncncbhd')
            ENDIF
            IF (IPRINT.GE.10) THEN
             write(lupri,'(/,2X,A)') 'Output from FINPOS...reading'
             write(lupri,*) 'INTEXP, INTORB, INTISG',
     &          INTEXP, INTORB, INTISG
            ENDIF
            RETURN
         END IF
      ELSE IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
         NEWEL = .TRUE.
         IF (IPRINT.GE.10) THEN
            write(lupri,'(/,2X,A)') 'Output from FINPOS'
            write(lupri,*) 'INTEXP, INTORB, INTISG',
     &         INTEXP, INTORB, INTISG
         ENDIF
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

      IF (IPRINT.GE.10) THEN
        write(lupri,'(/,2X,A)') 'Output from FINPOS'
        write(lupri,*) 'INTEXP, INTORB, INTISG',
     &          INTEXP, INTORB, INTISG
      ENDIF

      RETURN
C
 1000 CONTINUE
      WRITE (LUPRI,'(/A,I5)') 'Error in reading file, your basis '//
     &     'is not complete. IOSTAT code:',IOERR
      CALL QUIT('Incomplete basis set file in HERBAS')
      END

      SUBROUTINE print_gaussian_inp(icount, iqm,ALPHA,CPRIMU,NPRI,
     &                              NRCI,SEG,KPRIM,DSM)
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
      LOGICAL SEG
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
      character(1) :: ch
      if (iqm == 1) then
         icount = icount + 1
         if (icount > 1) then
            write(*, '(a)') '****'
         end if
         write(*, '(a2, i2)') namn(icount), 0
      end if
C
C     *****************************
C     *** Segmented contraction ***
C     *****************************
C
      IDONE = 0
      IF (SEG) THEN
         DO 100 ICONTR = 1, NRCI
            ISTART = IDONE + 1
            NLEFT  = NPRI - IDONE
C
C           Find first primitive
C           ====================
C
            IMXA = ISTART + IDAMAX(NLEFT,ALPHA(ISTART),1) - 1
            CALL DSWAP(1,ALPHA(ISTART),1,ALPHA(IMXA),1)
            CALL DSWAP(NRCI,CPRIMU(ISTART,1),KPRIM,CPRIMU(IMXA,1),
     &                 KPRIM)
C
C           Find corresponding contracted function
C           ======================================
C
            IMXC = IDAMAX(NRCI,CPRIMU(ISTART,1),KPRIM)
            CALL DSWAP(NPRI,CPRIMU(1,ICONTR),1,CPRIMU(1,IMXC),1)
C
C           Collect other primitives contributing to same contracted
C           ========================================================
C
            IPRI = 1
            DO 200 I = ISTART + 1, NPRI
               IF (ABS(CPRIMU(I,ICONTR)) .GT. DSM) THEN
                  CALL DSWAP(1,ALPHA(I),1,ALPHA(ISTART+IPRI),1)
                  CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                            CPRIMU(ISTART+IPRI,1),KPRIM)
                  IPRI = IPRI + 1
               END IF
  200       CONTINUE
C
C           Sort primitives
C           ===============
C
            IF (IPRI .GT. 2) THEN
               DO 300 I = ISTART + 1, ISTART + IPRI - 2
                  DO 400 J = I + 1, ISTART + IPRI - 1
                     IF (ALPHA(J) .GT. ALPHA(I)) THEN
                        CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
                        CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                                  CPRIMU(J,1),KPRIM)
                     END IF
  400             CONTINUE
  300          CONTINUE
            END IF
C
            IDONE = IDONE + IPRI
  100    CONTINUE
C
C     ***************************
C     *** General contraction ***
C     ***************************
C
      ELSE
         DO 500 I = 1, NPRI - 1
            DO 600 J = I + 1, NPRI
            IF(ALPHA(J) .GT. ALPHA(I)) THEN
               CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
               CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,CPRIMU(J,1),KPRIM)
           END IF
  600      CONTINUE
  500   CONTINUE
      END IF

      select case (iqm)
         case (1)
            ch = 'S'
         case (2)
            ch = 'P'
         case (3)
            ch = 'D'
         case (4)
            ch = 'F'
         case (5)
            ch = 'G'
         case (6)
            ch = 'H'
         case default
            call quit('too high iqm')
      end select

      do i = 1, npri
         write(*, '(a2, a)') ch, '   1 1.00       0.000000000000'
         write(*, '(e22.10, e18.10)') alpha(i), 1.0d0
      end do

      END
C  -- end of herbas.F --
