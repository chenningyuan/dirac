C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
C     M. Watson
C
C  /* Deck fmminp */
      SUBROUTINE FMMINP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "cbifmm.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 23)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "clsfmm.h"
#include "cbifck3.h"
      DATA TABLE /'.SKIP  ','.PRINT ','.STOP  ','.LMAX  ',
     &            '.NLEVEL','.FQUAD ','.NNQUAD','.BQUAD ',
     &            '.NLOGN ','.FMM   ','.GRSRCH','.TLMAX ',
     &            '.UMAT  ','.BRFREE','.GRAIN ','.DYNLMA',
     &            '.ALLSQR','.TRSRCH','.NOBOXP','.PRSTAT',
     &            '.TCON  ','.SCREEN','.SKIPNN'/
C
      LOGICAL USEUMAT,MMBRFREE,DYNLMAX,ALLSQR,TRSRCH,
     &        NOBOXP,GRSRCH,INCNN
      INTEGER NLEVEL,TLMAX,ALGORITHM,LEXTRA,TCON,
     &        DO_FQ,DO_NN,DO_BQ,DO_NLOGN,DO_FMM
      DOUBLE PRECISION GRAIN,RPQMIN,SCREEN
C
c      THRCLS = 1D-15
C
      NEWDEF = (WORD .EQ. '*FMM   ')
      WORD1 = WORD
C
C     initialise default parameters from MM code which are used internally;
C     external parameters (e.g. LMAX for multipole generation) should be
C     initialised elsewhere (earlier);
C     this is done so that MM code is reset with the defaults if in fact
C     no keywords are specified here.
C
      CALL mm_get_named_paras(DO_FQ,DO_NN,DO_BQ,DO_NLOGN,DO_FMM)
      CALL mm_init_defs(NLEVEL,GRAIN,TLMAX,ALGORITHM,
     &                  USEUMAT,MMBRFREE,RPQMIN,DYNLMAX,LEXTRA,
     &                  ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,TCON,SCREEN)
C
C     Now check for input changes to both internal and external
C     MM run-type parameters
C
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,
     &                       13,14,15,16,17,18,19,20,21,22,23), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in FMMINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in FMMINP.')
    1          CONTINUE
                  MMSKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. MMPRDF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  MMCUT  = .TRUE.
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) LMAX 
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) NLEVEL 
               GO TO 100
    6          CONTINUE
                  ALGORITHM = DO_FQ
               GO TO 100
    7          CONTINUE
                  ALGORITHM = DO_NN
               GO TO 100
    8          CONTINUE
                  ALGORITHM = DO_BQ
               GO TO 100
    9          CONTINUE
                  ALGORITHM = DO_NLOGN
               GO TO 100
   10          CONTINUE
                  ALGORITHM = DO_FMM
               GO TO 100
   11          CONTINUE
                  GRSRCH = .TRUE.
               GO TO 100
   12          CONTINUE
                  READ (LUCMD,*) TLMAX 
               GO TO 100
   13          CONTINUE
                  USEUMAT = .TRUE.
               GO TO 100
   14          CONTINUE
                  MMBRFREE = .TRUE. 
                  READ (LUCMD,*) RPQMIN 
                  RPQMIN2 = RPQMIN*RPQMIN
               GO TO 100
   15          CONTINUE
                  READ (LUCMD,*) GRAIN 
               GO TO 100
   16          CONTINUE
                  DYNLMAX = .TRUE. 
                  READ (LUCMD,*) LEXTRA 
               GO TO 100
   17          CONTINUE
                  ALLSQR = .TRUE.
               GO TO 100
   18          CONTINUE
                  TRSRCH = .TRUE.
               GO TO 100
   19          CONTINUE
                  NOBOXP = .TRUE.
               GO TO 100
   20          CONTINUE
               CALL mm_enable_stats
               GO TO 100
   21          CONTINUE
                  READ (LUCMD,*) TCON 
               GO TO 100
   22          CONTINUE
                  READ (LUCMD,*) SCREEN 
               GO TO 100
   23          CONTINUE
                  INCNN = .FALSE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in FMMINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in FMMINP.')
            END IF
      END IF
  300 CONTINUE
C
C The MM code has it's own input validation and processing
C but here we use standard Dalton procedure for the input
C elements which have influence outside the MM f90 modules.
C
      CALL HEADER('Changes of defaults for FMMDRV:',0)
      IF (MMSKIP) THEN
         WRITE (LUPRI,'(A)') ' FMMDRV skipped in this run.'
      ELSE
         IF (IPRINT .NE. MMPRDF) THEN
            WRITE (LUPRI,'(A,I5)') ' Print level in FMMDRV:',IPRINT
         END IF
         IF (LMAX .NE. LMAXDF) THEN
            WRITE (LUPRI,'(A,I5)') 
     &           ' Generation of multipole moments to order: ',LMAX
         END IF
         IF (MMCUT) THEN
            WRITE (LUPRI,'(/,A)') ' Program is stopped after FMMDRV.'
         END IF
      END IF
C
C     pass input to MM code, validate, and store in f90 code for
C     use on eventual main entry into MM driver.
C
      CALL mm_init_scheme(NLEVEL,GRAIN,LMAX,TLMAX,ALGORITHM,
     &                   USEUMAT,MMBRFREE,RPQMIN,DYNLMAX,LEXTRA,
     &                   ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,TCON,
     &                   SCREEN,CONTRACTED)
      RETURN
      END
C
C  /* Deck fmmini */
      SUBROUTINE FMMINI
C
C     M. Watson
C
#include "implicit.h"
#include "cbifmm.h"
#include "clsfmm.h"
C
C This initialises both internal and external
C  MM run-type parameters
C
      THRCLS = 1D-15
      IPRINT = MMPRDF
      MMSKIP = MMSKIPDF
      MMCUT  = MMCUTDF
C     order of multipole moment generation
      LMAX   = LMAXDF
C
      CALL mm_set_def_scheme(LMAX)
      RETURN
      END
C
      SUBROUTINE FMMONE(HMAT,NBAST,NNBASX,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "maxmom.h"
#include "iratdef.h"
#include "cbifmm.h" 
#include "cbifck3.h"
C

      DIMENSION HMAT(NNBASX), WORK(LWORK)
C      
      IF (MMSKIP) RETURN                           

      IPRINT = 0

      T1 = SECOND()
 
      IF ( OLDMOM ) THEN
       NCOMP = 2*LMAX + 1
       CALL PUT2INFORB(NBAST,NNBASX)
       KINTS = 1
       KDMAT = KINTS + NCOMP*NNBASX
       KHMSQ = KDMAT + NNBASX
       KLAST = KHMSQ + NBAST*NBAST
       IF (KLAST .GT. LWORK) CALL STOPIT('FMMFCK',' ',KLAST,LWORK)
       LWRK  = LWORK - KLAST + 1
       CALL FMMONE2(HMAT,WORK(KHMSQ),WORK(KINTS),WORK(KDMAT),
     &              NCOMP,WORK(KLAST),LWRK)
      ELSE
       KFREE = 1
       LFREE = LWORK

       CALL SETAOBTCH(.FALSE.,.TRUE.,WORK,LWORK,IPRINT)   

       CALL SETODBTCH(KIODBCLS,KIODBTCH,KRODBTCH,
     &                KFREE,LFREE,WORK,LWORK,IPRINT,.false.)
C      WORK(KFREE) is here just a dummy array for SETMM
       CALL SETMM_OR_JFF(WORK(KFREE),0,dummy,WORK(KIODBCLS),
     &            WORK(KIODBTCH),WORK(KRODBTCH),
     &            WORK(KFREE),LFREE,IPRFCK3,'FMM',dummy)

       KHMSQ = KFREE
       KLAST = KHMSQ + NBAST*NBAST
       IF (KLAST .GT. LWORK) CALL STOPIT('FMMFCK',' ',KLAST,LWORK)
       LWRK  = LWORK - KLAST + 1

       CALL FMMONE1(HMAT,WORK(KHMSQ),NBAST,WORK(KLAST),LWRK)

      END IF

      TTOT = SECOND()-T1
c      CALL TIMTXT('>>> TIME USED in mm_get_J_matrix(1el)',TTOT,LUPRI)

      RETURN
      END

      SUBROUTINE PUT2INFORB(NBAST_VAL,NNBASX_VAL)
#include "implicit.h"
#include "inforb.h"
      NBAST  = NBAST_VAL
      NNBASX = NNBASX_VAL
      RETURN
      END

      SUBROUTINE FMMONE1(HMAT,HMSQ,NBAST,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "clsfmm.h"
#include "cbifmm.h" 
C
      DIMENSION HMAT(*),HMSQ(NBAST,NBAST),WORK(LWORK)
C
      CALL DZERO(HMSQ,NBAST*NBAST)
      CALL mm_get_J_matrix('ONE_EL',HMSQ,NBAST,NBAST,WORK,LWORK)

      IJ = 0
      DO I = 1, NBAST
         DO J = 1, I 
            IJ = IJ + 1
            IF (I.EQ.J) THEN
               HMAT(IJ) = HMAT(IJ) + HMSQ(I,I)
            ELSE
               HMAT(IJ) = HMAT(IJ) + 0.5d0*(HMSQ(I,J) + HMSQ(J,I))
            END IF
         END DO
      END DO

      RETURN
      END

      SUBROUTINE FMMONE2(HMAT,HMSQ,SPHINT,DENMAT,NCOMP,WORK,LWORK)
C
C     M. Watson & T. Helgaker 
C
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "clsfmm.h"
#include "cbifmm.h" 
      DIMENSION HMAT(NNBASX),HMSQ(NBAST,NBAST), 
     &          SPHINT(NNBASX,NCOMP),DENMAT(NNBASX),WORK(LWORK) 
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMDR1','*',103)
      WRITE(6,*) NBAST,NNBASX
C
C     Dummy nonzero density matrix
C     ============================
C
      DO I = 1, NNBASX
       DENMAT(I) = 1.0D0
      END DO
C
C     Multipole moment integrals 
C     ==========================
C
      PRINT *,"generating MM integrals to order", LMAX
      CALL FMMGEN(SPHINT,DENMAT,LMAX,WORK,LWORK)
C
      CALL DZERO(HMSQ,NBAST*NBAST)
C
      CALL mm_get_J_matrix('ONE_EL',HMSQ,NBAST,NBAST,WORK,LWORK)

      IJ = 0
      DO I = 1, NBAST
       DO J = 1, I 
        IJ = IJ + 1
        IF (I.EQ.J) THEN
           HMAT(IJ) = HMAT(IJ) + HMSQ(I,I)
        ELSE
           HMAT(IJ) = HMAT(IJ) + 0.5d0*(HMSQ(I,J) + HMSQ(J,I))
        END IF
       END DO
      END DO

      RETURN
      END

C  /* Deck fmmfck */
      SUBROUTINE FMMFCK(FMAT,DMAT,WORK,LWORK)
C
C     M. Watson & T. Helgaker 
C
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "cbifmm.h" 
#include "cbifck3.h" 
#include "clsfmm.h"
#include "densfit.h"
C
      DIMENSION FMAT(*), DMAT(*), WORK(LWORK)
C
      IF (MMSKIP) RETURN                           
      IF (NSYM .GT. 1) THEN
         WRITE (LUPRI,'(2X,2A)')
     &      ' Symmetry not implemented in FMMFCK.',
     &      ' Program stopped.'
         CALL QUIT('Symmetry not implemented in FMMFCK.')
      END IF

      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMFCK','*',103)
C
      T1 = SECOND()
C
      IF (OLDMOM) THEN
       NCOMP = 2*LMAX + 1
C
C     Allocations
C    
       KINTS = 1
       KDMAT = KINTS + NCOMP*NNBASX
       KLAST = KDMAT + NNBASX
       IF (KLAST .GT. LWORK) CALL STOPIT('FMMFCK',' ',KLAST,LWORK)
       LWRK  = LWORK - KLAST + 1
       CALL FMMFC1(FMAT,DMAT,WORK(KINTS),WORK(KDMAT),
     &             NCOMP,WORK(KLAST),LWRK)
      ELSE
       T1 = SECOND()
       NDIM1 = NBAST
       NDIM2 = NBAST
       IF (DOFIT.AND.LHSFIT) THEN
         NDIM1 = 1
         NDIM2 = NBASEAUX
       ENDIF
       CALL mm_get_J_matrix('TWO_EL',FMAT,NDIM1,NDIM2,WORK,LWORK)
       TTOT = SECOND()-T1
c       CALL TIMTXT('>>> TIME USED in mm_get_J_matrix(2el)', TTOT, LUPRI)
      END IF

      TTOT = SECOND()-T1
C     CALL TIMTXT('>>> TIME USED in _mmfck', TTOT, LUPRI)
      WRITE(LUPRI,*)

      RETURN
      END
C

C  /* Deck fmmdrv */
      SUBROUTINE FMMFC1(FMAT,DMAT,SPHINT,DENMAT,NCOMP,WORK,LWORK)
     &                  
C
C     T. Helgaker 
C
#include "implicit.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxmom.h"
#include "maxorb.h"
#include "aovec.h"
C
      LOGICAL   FIRST
      PARAMETER (D2 = 2.0D0)
      DIMENSION FMAT(*), DMAT(NBAST,NBAST),
     &          SPHINT(NNBASX,NCOMP), DENMAT(NNBASX), WORK(LWORK)
C
#include "inforb.h"
#include "inftap.h"
#include "infinp.h"
#include "nuclei.h"
#include "cbiher.h"
#include "cbieri.h"
#include "clsfmm.h"
#include "cbifmm.h" 
#include "densfit.h" 
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMDR1','*',103)
C
C     AO density matrix
C     =================
C
      IJ = 0
      DO I = 1, NBAST
         DO J = 1, I
            IJ = IJ + 1
            IF (I.EQ.J) THEN
               DENMAT(IJ) = DMAT(J,I)
            ELSE
               DENMAT(IJ) = D2*DMAT(J,I)
            END IF
         END DO
      END DO
C
C     Multipole moment integrals 
C     ==========================
C
      PRINT *,"generating MM integrals to order", LMAX
      CALL FMMGEN(SPHINT,DENMAT,LMAX,WORK,LWORK)
C
      T1 = SECOND()
       NDIM1 = NBAST
       NDIM2 = NBAST
       IF (DOFIT.AND.LHSFIT) THEN
         NDIM1 = 1
         NDIM2 = NBASEAUX
       ENDIF
      CALL mm_get_J_matrix('TWO_EL',FMAT,NDIM1,NDIM2,WORK,LWORK)
      TTOT = SECOND()-T1
      CALL TIMTXT('>>> TIME USED in mm_get_J_matrix(2el)', TTOT, LUPRI)
C
      RETURN
      END

C  /* Deck fmmdrv */
      SUBROUTINE FMMDRV(WORK,LWORK)
C
C     T. Helgaker 
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "maxmom.h"
#include "iratdef.h"
#include "inforb.h"
#include "cbifmm.h" 
#include "cbifck3.h" 
      DIMENSION WORK(LWORK)
C
      IF (MMSKIP) RETURN                           
      IF (NSYM .GT. 1) THEN
         WRITE (LUPRI,'(2X,2A)')
     &      ' Symmetry not implemented in FMMDRV.',
     &      ' Program stopped.'
         CALL QUIT('Symmetry not implemented in FMMDRV.')
      END IF
      IPRINT = 0
C
      NCOMP = 2*LMAX + 1
C
C     Allocations
C    
      KINTS = 1
      KDMAT = KINTS + NCOMP*NNBASX
      KLAST = KDMAT + NNBASX
      IF (KLAST .GT. LWORK) CALL STOPIT('FMMDRV',' ',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
C
      T1 = SECOND()
      CALL FMMDR1(WORK(KINTS),WORK(KDMAT),NCOMP,WORK(KLAST),LWRK)
      TTOT = SECOND()-T1
      CALL TIMTXT('>>> TIME USED in FMMDR1', TTOT, LUPRI)
C
      RETURN
      END
C
C  /* Deck fmmdr1 */
      SUBROUTINE FMMDR1(SPHINT,DENMAT,NCOMP,WORK,LWORK)
C
C     T. Helgaker 
C
#include "implicit.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxmom.h"
#include "maxorb.h"
#include "aovec.h"
C
      DIMENSION SPHINT(NNBASX,NCOMP),DENMAT(NNBASX),WORK(LWORK)
C
#include "inforb.h"
#include "inftap.h"
#include "infinp.h"
#include "nuclei.h"
#include "cbiher.h"
#include "cbieri.h"
#include "clsfmm.h"
#include "cbifmm.h" 
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMDR1','*',103)
C
C     Calculate nonclassical energy
C     =============================
C
      T1 = SECOND()
      CALL FMMNCL(WORK,LWORK)
      TTOT = SECOND()-T1
      WRITE (LUPRI,*)
      CALL TIMTXT('>>> TIME USED in FMMNCL', TTOT, LUPRI)
      T1 = SECOND()
C
C     AO density matrix
C     =================
C
      CALL FMMDNS(DENMAT,WORK,LWORK,IPRINT)
C
C     Multipole moment integrals 
C     ==========================
C
      PRINT *,"Calculating MM integrals to order", LMAX
      CALL FMMGEN(SPHINT,DENMAT,LMAX,WORK,LWORK)
      TTOT = SECOND()-T1
      WRITE (LUPRI,*)
      CALL TIMTXT('>>> TIME USED in FMMINT', TTOT, LUPRI)
C
      T1 = SECOND()
      CALL mm_get_energy(EL2CLS,WORK,LWORK)
      TTOT = SECOND()-T1
      WRITE (LUPRI,*)
      CALL TIMTXT('>>> TIME USED in mm_get_energy', TTOT, LUPRI)
      WRITE(*,'(">>> TIME USED in mm_get_energy  =",F14.2,/)') TTOT
C
      WRITE (LUPRI,'(/," Returned to main DALTON code.",/)')
      RETURN
      END
C  /* Deck fmmncl */
      SUBROUTINE FMMNCL(WORK,LWORK)
C
C     T. Helgaker 
C
#include "implicit.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxmom.h"
#include "maxorb.h"
#include "aovec.h"
C
      PARAMETER (D0 = 0.0D0, D100 = 1.0D2)
      DIMENSION WORK(LWORK)
C
#include "inforb.h"
#include "inftap.h"
#include "infinp.h"
#include "nuclei.h"
#include "cbiher.h"
#include "cbieri.h"
#include "erienr.h"
#include "clsfmm.h"
#include "cbifmm.h"
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMNCL','*',103)
C
      KDMAT  = 1
      KHMAT  = KDMAT  + N2BASX
      KCCFBT = KHMAT  + NNBAST
      KINDXB = KCCFBT + MXCONT*MXPRIM
      KLAST  = KINDXB + 8*MXSHEL*MXCONT
      IF (KLAST .GT. LWORK) CALL STOPIT('FMMDR1','2EL',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
C
C     Density matrix
C     ==============
C
!radovan:
      call quit('subroutine GETDMT missing')
!     CALL GETDMT(WORK(KDMAT),1,WORK(KLAST),LWRK,.FALSE.,.TRUE.,.TRUE.,
!    &            IPRINT)
C
C     One-electron interactions
C     =========================
C
C     total
C
      LUONEL = -1
      CALL GPOPEN(LUONEL,'AOONEINT',' ',' ',' ',IDUMMY,.FALSE.)
      CALL READIN(WORK(KLAST),LWRK,.TRUE.)
      CALL ONEDRV(WORK(KLAST),LWRK,IPRINT,.FALSE.,0,.FALSE.,.TRUE.,
     &            .TRUE.,.FALSE.,.FALSE.,.FALSE.)
!radovan:
      call quit('subroutine GETDMT missing')
!     CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KHMAT),NNBAST)
      CALL GET1EL(ERGONE,WORK(KDMAT),WORK(KHMAT),IPRINT)
!     CALL RDONEL('KINETINT',.TRUE.,WORK(KHMAT),NNBAST)
      CALL GET1EL(ERGKIN,WORK(KDMAT),WORK(KHMAT),IPRINT)
C
C     nonclassical 
C
      LUONEL = -1
      CALL GPOPEN(LUONEL,'AOONEINT',' ',' ',' ',IDUMMY,.FALSE.)
      CALL READIN(WORK(KLAST),LWRK,.TRUE.)
      CALL ONEDRV(WORK(KLAST),LWRK,IPRINT,.FALSE.,0,.FALSE.,.TRUE.,
     &            .TRUE.,.FALSE.,.TRUE.,.TRUE.)
!radovan:
      call quit('subroutine GETDMT missing')
!     CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KHMAT),NNBAST)
      CALL GET1EL(ERGNCL,WORK(KDMAT),WORK(KHMAT),IPRINT)
C
C     contributions
C
      EL1TOT = ERGONE - ERGKIN
      EL1NON = ERGNCL - ERGKIN 
      EL1CLS = EL1TOT - EL1NON 
C
C     Two-electron interactions
C     =========================
C
C     total

      CALL HRFCK2(DUMMY,WORK(KDMAT),1,IDUMMY,0,IPRINT,WORK(KCCFBT),
     &            WORK(KINDXB),DUMMY,WORK(KLAST),LWRK)
      ALLTOT = ERGTOT
      COUTOT = ERGCOU
C
C     nonclassical
C
      NCLERI = .TRUE.
      CALL HRFCK2(DUMMY,WORK(KDMAT),1,IDUMMY,0,IPRINT,WORK(KCCFBT),
     &            WORK(KINDXB),DUMMY,WORK(KLAST),LWRK)
      NCLERI = .FALSE.
      ALLNCL = ERGTOT
      COUNCL = ERGCOU
C
C     contributions
C
      ERGEXC = ALLTOT - COUTOT
      EL2TOT = COUTOT
      EL2NON = COUNCL
      EL2CLS = EL2TOT - EL2NON
C
C     Print
C     =====
C
      ELNNON = EL1NON + EL2NON
      ELNCLS = EL1CLS + EL2CLS + POTNUC
      ELNTOT = ELNNON + ELNCLS 
C
      ETOT = ELNTOT + ERGKIN + ERGEXC
C
      write(*,'(/,A,/)') " ERI analysis of Hartree energy:"
C
      WRITE (*,'(2X,17X,2A,/)') '       non-classical',
     &             '         classical           total'
      WRITE (*,'(2X,A,E18.9,E21.12,F16.8)') 
     &       ' nuclear attraction',EL1NON,EL1CLS,EL1TOT
      WRITE (*,'(2X,A,E18.9,E21.12,F16.8)') 
     &       ' electron repulsion',EL2NON,EL2CLS,EL2TOT
      WRITE (*,'(2X,A,E18.9,E21.12,F16.8)') 
     &       ' nuclear repulsion ',D0,POTNUC,POTNUC
      WRITE (*,'(2X,A,E18.9,E21.12,F16.8,/)') 
     &       ' total Coulomb     ',ELNNON,ELNCLS,ELNTOT
C
      CALL HEADER('ERI decomposition of Hartree energy',-1)
      WRITE (LUPRI,'(2X,17X,2A,/)') '       non-classical',
     &             '         classical           total'
      WRITE (LUPRI,'(2X,A,E18.9,E21.12,F16.8)') 
     &       ' nuclear attraction',EL1NON,EL1CLS,EL1TOT
      WRITE (LUPRI,'(2X,A,E18.9,E21.12,F16.8)') 
     &       ' electron repulsion',EL2NON,EL2CLS,EL2TOT
      WRITE (LUPRI,'(2X,A,E18.9,E21.12,F16.8)') 
     &       ' nuclear repulsion ',D0,POTNUC,POTNUC
      WRITE (LUPRI,'(2X,A,E18.9,E21.12,F16.8)') 
     &       ' total Coulomb     ',ELNNON,ELNCLS,ELNTOT
C
C     WRITE (LUPRI,'(2X,A,F18.10,F21.14,F16.8)') 
C    &       ' nuclear attraction',EL1NON,EL1CLS,EL1TOT
C     WRITE (LUPRI,'(2X,A,F18.10,F21.14,F16.8)') 
C    &       ' electron repulsion',EL2NON,EL2CLS,EL2TOT
C     WRITE (LUPRI,'(2X,A,F18.10,F21.14,F16.8)') 
C    &       ' nuclear repulsion ',D0,POTNUC,POTNUC
C     WRITE (LUPRI,'(2X,A,F18.10,F21.14,F16.8)') 
C    &       ' total Coulomb     ',ELNNON,ELNCLS,ELNTOT
C      CALL HEADER('Decomposition of total electronic energy',-1)
C     WRITE (LUPRI,'(2X,A,F20.13,4X,F10.2,A)') 
C    &       ' kinetic energy    ',ERGKIN, D100*ERGKIN/ETOT,' %'
C     WRITE (LUPRI,'(2X,A,F20.13,4X,F10.2,A)') 
C    &       ' Coulomb energy    ',ELNTOT, D100*ELNTOT/ETOT,' %'
C     WRITE (LUPRI,'(2X,A,F20.13,4X,F10.2,A)') 
C    &       ' exchange energy   ',ERGEXC, D100*ERGEXC/ETOT,' %'
C     WRITE (LUPRI,'(2X,A,F20.13,4X,F10.2,A)') 
C    &       ' total energy      ',ETOT, D100,' %'
C
C
      RETURN
      END
C  /* Deck fmmgen */
      SUBROUTINE FMMGEN(SPHINT,DENMAT,LORDER,WORK,LWORK)
C
C     Pawel Salek & T. Helgaker 
C
#include "implicit.h"
#include "mxcent.h"
c
      DIMENSION SPHINT(NNBASX,*), DENMAT(*), WORK(LWORK)
c
#include "priunit.h"
#include "inforb.h"
#include "cbifmm.h"
#include "nuclei.h"
c
      KCNTS  = 1
      KDISM  = KCNTS  + NNBASX
      KFMAT  = KDISM  + 14*(NNBASX + NUCIND)
      KINDEX = KFMAT  +  2*NBAST*NBAST 
      KLAST  = KINDEX +  2*(NNBASX + NUCIND)
      IF (KLAST .GT. LWORK) CALL STOPIT('SPHAVE',' ',KLAST,LWORK)
      LWRK = LWORK - KLAST +1
      CALL DZERO(WORK(KCNTS),NNBASX)
c
      LUINTM = -1
      CALL GPOPEN(LUINTM,'MM_DATA','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
C
C     Electronic information
C     ======================
C
      NMOMNT = 0
      NDISTR = 0
      DO LCUR=0, LORDER
         CALL FMMINT(SPHINT,LCUR,WORK(KLAST),LWRK,IPRINT)
C
C     Expectation values
C     ==================
C
         CALL SPHAVE(LUINTM,SPHINT,DENMAT,LCUR,NDISTR,NMOMNT,
     &               WORK(KCNTS),WORK(KDISM))
      END DO
c
c     write to file the count of generated moments
c
      WRITE (LUINTM) -1,0, 0,0, 0,0,0, 0, 0D0, 0D0,0D0,0D0, 0D0,0D0
C
C     Nuclear information
C     =================== 
C
      DO I = 1, NUCIND
         WRITE (LUINTM) CHARGE(I), CORD(1,I), CORD(2,I), CORD(3,I)
      END DO
C
      CALL GPCLOSE(LUINTM,'KEEP')
c
c     Write counts to a separate file for fast access.
      LUINTM = -1
      CALL GPOPEN(LUINTM,'MM_CNTS','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      WRITE (LUINTM) NDISTR, NUCIND, NMOMNT
      WRITE (LUPRI,'(/,2X,A,2I8,2X,F12.6)')
     &      ' Number of nonzero distributions in FMMGEN: ', NDISTR
      WRITE (LUPRI,'(2X,A,2I12,2X,F12.6)')
     &      ' Number of moments written in FMMGEN:   ', NMOMNT
C
      CALL GPCLOSE(LUINTM,'KEEP')
c
      END
C  /* Deck fmmdns */
      SUBROUTINE FMMDNS(DENMAT,WORK,LWORK,IPRINT)
C
C     T. Helgaker 
C
#include "implicit.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
      LOGICAL NODC, NODV
      DIMENSION DENMAT(NNBASX), WORK(LWORK)
#include "inforb.h"
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMDNS','*',103)
C
      KFOCK = 1
      KLAST = KFOCK + NNBASX
      IF (KLAST .GT. LWORK) CALL STOPIT('FMMDNS',' ',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
C
      NODC = .FALSE.
      NODV = .TRUE.
      CALL DZERO(DENMAT,NNBASX)
      CALL DZERO(WORK(KFOCK),NNBASX)
      CALL DSOFSO(DENMAT,WORK(KFOCK),WORK(KLAST),LWRK,IPRINT,NODC,NODV)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Density matrix in FMMDNS',-1)
         CALL OUTPAK(DENMAT,NBAST,1,LUPRI)
      END IF
      RETURN
      END
      SUBROUTINE FMMINT(SPHINT,LMAX,WORK,LWORK,IPRINT)
C
C     T. Helgaker 
C
#include "implicit.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxmom.h"
      LOGICAL TOFILE, TRIMAT, DOINT(2,2)
      DIMENSION SPHINT(NNBASX,*), WORK(LWORK)
#include "inforb.h"
#include "nuclei.h"
#include "cbiher.h"
#include "orgcom.h"
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMINT','*',103)
C
      NCMPS = (LMAX + 1)**2
C
      KLABL = 1 
      KINTR = KLABL + NCMPS
      KINTA = KINTR + (NCMPS + IRAT - 1)/IRAT
      KLAST = KINTA + (NCMPS + IRAT - 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('FMMINT',' ',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
C
      NCOMP = 1
      DOINT(1,1) = .TRUE.
      DOINT(2,1) = .FALSE.
      DOINT(1,2) = .FALSE.
      DOINT(2,2) = .FALSE.
      TRIANG = .TRUE.
      TOFILE = .FALSE.
      PROPRI = IPRINT .GT. 5
C
      CALL PR1IN1(WORK(KLAST),LWRK,WORK(KINTR),WORK(KINTA),
     &            WORK(KLABL),'SPHMOML',LMAX,MPQUAD,TRIANG,PROPRI,
     &            IPRINT,SPHINT,NCOMP,TOFILE,'TRIANG',
     &            DOINT)
      IF (IPRINT .GT. 5) CALL FMMSPR(SPHINT,WORK(KLABL),
     &                               WORK(KINTR),NCMPS)

C      CALL FMMSPR(SPHINT,WORK(KLABL),
C     &                               WORK(KINTR),NCMPS)
      RETURN
      END
C
C  /* Deck fmmspr */
      SUBROUTINE FMMSPR(SPHINT,LABINT,INTREP,NCOMPS)
C
C     T. Helgaker 
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
      CHARACTER LABINT(*)*8
      DIMENSION SPHINT(NNBASX,*), INTREP(*) 
#include "inforb.h"
C
      DO ICOMP = 1, NCOMPS
         CALL AROUND('Integrals of operator: '//LABINT(ICOMP))
         WRITE (LUPRI,'(A,I2)') ' Symmetry of operator:',
     &          INTREP(ICOMP) + 1
         CALL OUTPAK(SPHINT(1,ICOMP),NBAST,1,LUPRI)
      END DO
      RETURN
      END
C  /* Deck sphave */
      SUBROUTINE SPHAVE(LUINTM,SPHINT,DENMAT,L,NDISTR,NMOM,ICNTS,DIS)
C     compute moments for given moment L
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxmom.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "inforb.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL NODSPH
      DIMENSION SPHINT(NNBASX,*),DENMAT(NNBASX),ICNTS(NNBASX),DIS(*)
#include "ccom.h"
#include "primit.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "onecom.h"
#include "clsfmm.h"
#include "cbifmm.h"
C
c      CALL HEADER('Output from SPHAVE',-1)
C
      SPHMAX = D0
      DSPHMAX = D0
      IORBA = 0
      IJ  = 0
      DO ISHELA = 1,KMAX
      DO ICOMPA = 1, KHKT(ISHELA)
         IORBA = IORBA + 1
         IORBB = 0
         DO ISHELB = 1, KMAX
         DO ICOMPB = 1, KHKT(ISHELB)
            IORBB = IORBB + 1
            IF (IORBB .LE. IORBA) THEN
               IJ = IJ + 1
               NUCA   = NUCO (ISHELA)
               NUCB   = NUCO (ISHELB)
               JSTA   = JSTRT(ISHELA)
               JSTB   = JSTRT(ISHELB)
               NUMCFA = NUMCF(ISHELA)
               NUMCFB = NUMCF(ISHELB)
               CORAX  = CENT(ISHELA,1,1)
               CORAY  = CENT(ISHELA,2,1)
               CORAZ  = CENT(ISHELA,3,1)
               CORBX  = CENT(ISHELB,1,1)
               CORBY  = CENT(ISHELB,2,1)
               CORBZ  = CENT(ISHELB,3,1)
               LA     = NHKT(ISHELA) - 1
               LB     = NHKT(ISHELB) - 1
!radovan:
               call quit('subroutine GAUDST missing')
!              CALL GAUDST(PX,PY,PZ,PM,IPRINT)
CMAW FIXME: why are some of these centres broken?
               DO M = -L, L
c                 columns are ordered as 0, 1, -1, 2, -2, etc.
                  IF(M.GT.0) THEN
                     ICOMP = +2*M
                  ELSE
                     ICOMP = -2*M+1
                  END IF
                  DSPH = -DENMAT(IJ)*SPHINT(IJ,ICOMP)
		  IF (ABS(SPHINT(IJ,ICOMP)) .GT. SPHMAX) THEN
                     SPHMAX = ABS(SPHINT(IJ,ICOMP)) 
                     DSPHMAX = ABS(DSPH) 
                  END IF
                  IF (ABS(DSPH).GT.THRCLS) THEN
                   IF (ABS(SPHINT(IJ,ICOMP)).GT.THRCLS) THEN
                     NMOM = NMOM + 1
                     IF(ICNTS(IJ).EQ.0) THEN
                        NDISTR = NDISTR + 1
                        ICNTS(IJ) = NDISTR
                     END IF
                     WRITE (LUINTM) L,M, IORBA,IORBB, LA,LB, IJ, 
     &                             ICNTS(IJ), PM, PX,PY,PZ,
     &                            -SPHINT(IJ,ICOMP), DENMAT(IJ)
                     WRITE (6,'(A,8I4,A,4F10.6,E10.4,F10.6)') 
     &                            'L M IA IB LA LB IB ICNT',
     &                             L,M, IORBA,IORBB, LA,LB, IJ,
     &                 ICNTS(IJ),' P(M,X,Y,Z) S D', PM, PX,PY,PZ,
     &                            -SPHINT(IJ,ICOMP), DENMAT(IJ)
C                     WRITE (LUPRI,'(A,7I4,A,6F10.6)') 
C     &                            'L M IA IB LA LB ICNT',
C     &                             L,M,IORBA,IORBB,LA,LB,ICNTS(IJ),
C     &                            ' P(M,X,Y,Z) S D', PM, PX,PY,PZ,
C     &                            -SPHINT(IJ,ICOMP), DENMAT(IJ)
         	  END IF
                 END IF
	       END DO
	    END IF
	 END DO
	 END DO
      END DO
      END DO
C      
C     Screening check for significant contributions to J_matrix
C      
CFIXME: This will highlight possible errors, but doesn't solve
C       the problem should they be there!
C      
      IF ((SPHMAX .GT. 1d0) .AND. (DSPHMAX .LT.THRCLS)) THEN
	 WRITE(LUPRI,*) "SPHMAX ", SPHMAX
	 WRITE(LUPRI,*) "DSPHMAX", DSPHMAX
	 WRITE(LUPRI,*) "THRCLS ", THRCLS
	 WRITE(LUPRI,*)    
     &  "*** WARNING! moments screened that could be sig. for J_matrix!"
      END IF
      RETURN
      END
CC  /* Deck get1el */
      SUBROUTINE GET1EL(ERGONE,DMAT,HMAT,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.D0, D2 = 2.0D0)
      DIMENSION DMAT(NBAST,NBAST), HMAT(NNBAST)
#include "inforb.h"
      ERGONE = D0
      IJ = 0
      DO I = 1, NBAST
         DO J = 1, I - 1
            IJ = IJ + 1
            ERGONE = ERGONE + D2*DMAT(I,J)*HMAT(IJ)
         END DO
         IJ = IJ + 1
         ERGONE = ERGONE + DMAT(I,I)*HMAT(IJ)
      END DO
      RETURN
      END
c ---------------------------------------------------------------------------
#if defined(VAR_MPI)
      SUBROUTINE mm_getj_wake_up_slaves
      use interface_to_mpi
#include "maxorb.h"
#include "inforb.h"
#include "infpar.h"
      INTEGER IPRTYP, IPRINT, IERR
      IPRTYP = 6
      IPRINT = 0
      IF(mynum.EQ.MASTER) THEN
c        IPRTYP and IPRINT are received by eri slave dispatch routines
         CALL interface_MPI_Bcast(IPRTYP,1,0,global_communicator)
         CALL interface_MPI_Bcast(IPRINT,1,0,global_communicator)
c        nbast is received by fmm_slave_fock 
         CALL interface_MPI_Bcast(NBAST, 2,0,global_communicator)
      END IF
      END
!
      SUBROUTINE FMM_SLAVE_FOCK(WORK,LWORK,IPRINT)
!     FMM parallelization: Pawel Salek
      use interface_to_mpi
      DIMENSION WORK(LWORK)
      INTEGER MWORK, IPRINT
      CHARACTER MMTYPE(6)
#include "maxorb.h"
#include "infpar.h"
#include "priunit.h"
#include "densfit.h"
      CALL interface_MPI_Bcast(NBAST, 2,master,global_communicator)
      KHMSQ = 1
      LWRK  = KHMSQ+NBAST*NBAST
!     check mem here.
      CALL DZERO(WORK(KHMSQ), NBAST*NBAST)
       NDIM1 = NBAST
       NDIM2 = NBAST
       IF (DOFIT.AND.LHSFIT) THEN
         NDIM1 = 1
         NDIM2 = NBASEAUX
       ENDIF
      CALL mm_get_J_matrix(MMTYPE,WORK(KHMSQ),NDIM1,NDIM2,WORK,LWORK)
      END
#endif
