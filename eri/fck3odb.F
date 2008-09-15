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
C    Overlap Distribution classes :
C
C    IODBCLS(5,*) :
C
C      1 - ICENT1
C      2 - ICENT2
C      3 - NHKT12 ( = J1 + J2 + 1; J = ang. momentum )
C      4 - KSTBTCH - number of the first OD batch of 
C          this class ( used in IODBTCH )
C      5 - NBTCH   - number of OD batches in this class
C
C    RODBCLS(2,*) :
C
C      1 - Max CS screening parameter (Gab) in class
C      2 - Max density contracted CS Gab*Dab in class
C
C    --------------------------------------------------
C
C    Overlap Distribution batches :
C
C    IODBTCH(2,*) :
C
C      1 - IAO1 - First AO batch
C      2 - IAO2 - First AO batch
C
C    ( Note that AO batch is not the number of atomic 
C      orbital ( that is accessed via INDXBTCH array 
C      from fck3aob.h common block )  ! )
C  
C    RODBTCH(6+NDMT,*) 
C
C      1 - x-coodinate of the center
C      2 - y-coodinate of the center
C      3 - z-coodinate of the center
C      4 - non-classical extent (from the center)
C      5 - Cauchy screening parameter
C      6 - Max CS screening parameter (Gab) in batch
C      7-(6+NDMT) - Max density contracted CS Gab*Dab(IDMT) in batch

#if defined (VAR_SINGLE)
#include "single.h"
#endif
      SUBROUTINE SETODBTCH(KIODBCLS,KRODBCLS,KIODBTCH,KRODBTCH,
     &                KFREE,LFREE,WORK,LWORK,GAB,
     &                GABFIT,DMAT,CFIT,NDMT,IPRINT,tvecflag)
#include "implicit.h"
#include "priunit.h"
      DIMENSION WORK(*)
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#include "clsfmm.h"
#include "pi.h"
C
#include "pbc.h"
C
      DIMENSION GAB(*), GABFIT(*), DMAT(*), CFIT(*)
      logical tvecflag

      IF (IPRINT .GT. 5) CALL TITLER('Output from SETODBTCH','*',103)

      CALL GET_ODPAR(IPRINT)  

      CALL GETMAXODBCLS(MAXODBCLS,MAXODBTCH,IPRINT,tvecflag)

      KIODBCLS = KFREE
      KRODBCLS = KIODBCLS + 5*MAXODBCLS/IRAT + 1
      KIODBTCH = KRODBCLS + 2*MAXODBCLS
      KRODBTCH = KIODBTCH + 2*MAXODBTCH/IRAT + 1 
      KLAST    = KRODBTCH + (6+NDMT)*MAXODBTCH

      IF (IPRINT.GT.35) THEN
        WRITE(LUPRI,*) 'Temporary memory: '
        WRITE(LUPRI,*) 'KIODBCLS, KRODBCLS   ', KIODBCLS, KRODBCLS
        WRITE(LUPRI,*) 'KIODBTCH, KRODBTCH   ', KIODBTCH, KRODBTCH
        WRITE(LUPRI,*) 'Total temp. memory   ', KLAST - KFREE
        WRITE(LUPRI,*) 'Last work elem. used ', KLAST
        WRITE(LUPRI,*) 'Free mem. pre alloc. ', LFREE
        WRITE(LUPRI,*) 'Memory left in work', LFREE - (KLAST - KFREE)
      ENDIF

      LWRK     = LFREE - ( KLAST-KFREE) 

      IF ( LWRK .LT. 0 ) THEN 
        CALL STOPIT('SETODBTCH','FCK3ODB',KLAST-KFREE,LWORK)
      END IF 

      CALL SETODBTCH1(WORK(KIODBCLS),WORK(KRODBCLS),WORK(KIODBTCH),
     &                WORK(KRODBTCH), WORK(KLAST),LWRK,GAB,GABFIT,
     &                DMAT,CFIT,NDMT,IPRINT,tvecflag)

      IF ( ( NODBCLS .GT. MAXODBCLS ) .OR. ( NODBTCH .GT. MAXODBTCH ) ) 
     &   THEN 
        WRITE(LUPRI,*) 
     &   ' Problem in estimation of the maximum OD batches or classes !'
        WRITE(LUPRI,*) ' NODBCLS, MAXODBCLS :', NODBCLS, MAXODBCLS
        WRITE(LUPRI,*) ' NODBTCH, MAXODBTCH :', NODBTCH, MAXODBTCH
        CALL QUIT('Errors in SETODBTCH')
      END IF
 
C
C     Compress the class iformation
C     =============================
C 
C     sizes of class and batch informations
C
      LICLS  = 5*NODBCLS/IRAT + 1
      LRCLS  = 2*NODBCLS
      LIBTCH = 2*NODBTCH/IRAT + 1
      LRBTCH = (6+NDMT)*NODBTCH

      LTMP  = MAX(LIBTCH,LRBTCH)
      LTMP  = MAX(LTMP,LRCLS)

      KTMP  = KLAST 
      KLAST = KTMP  + LTMP

      IF (IPRINT.GT.35) THEN
        WRITE(LUPRI,*) 'Allocations before contraction: '
        WRITE(LUPRI,*) 'KIODBCLS, KRODBCLS  ', KIODBCLS, KRODBCLS
        WRITE(LUPRI,*) 'KIODBTCH, KRODBTCH  ', KIODBTCH, KRODBTCH
        WRITE(LUPRI,*) 'KIODBTCH, KRODBTCH  ', KIODBTCH, KRODBTCH
        WRITE(LUPRI,*) 'KTMP, KLAST,  LTMP  ', KTMP, KLAST, LTMP
        WRITE(LUPRI,*) 'LICLS, LRCLS        ', LICLS, LRCLS
        WRITE(LUPRI,*) 'LIBTCH, LRBTCH      ', LIBTCH, LRBTCH
        WRITE(LUPRI,*) 'Total temp. memory  ', KLAST - KFREE
        WRITE(LUPRI,*) 'Memory left in work ', LFREE - (KLAST - KFREE)
      ENDIF

c     if (.not. (dopbc .and. pbcactive)) then
c        WRITE(LUPRI,*) 'Max memory usage in SETODBTCH', KLAST - KFREE
c      CALL FLSHFO(LUPRI)
c     end if

      IF ( LTMP .GT. LWRK ) THEN 
        CALL STOPIT('SETODBTCH','ODBTCHTMP',KLAST,LWORK)
      END IF 

C
C     copy the batch information to the temporary place 
C     and than to its final place.
C   ( tmp copying just to avoid any possible problems
C     with copying array to overlaping place )    
C

      IF (IPRINT.GT.35) THEN
        WRITE(LUPRI,*) 'Allocations after contraction: '
        WRITE(LUPRI,*) 'KTMP + 2*NODBCLS ', KTMP + 2*NODBCLS
        WRITE(LUPRI,*) 'KRODBCLS         ', KIODBCLS + LICLS
        WRITE(LUPRI,*) 'KTMP + 2*NODBTCH ', KTMP + 2*NODBTCH/IRAT + 1
        WRITE(LUPRI,*) 'KIODBTCH         ', KRODBCLS + LRCLS
        WRITE(LUPRI,*) 'KTMP + (6+NDMT)*NODBTCH ', 
     &                  KTMP + (6+NDMT)*NODBTCH
        WRITE(LUPRI,*) 'KRODBTCH         ', KIODBTCH + LIBTCH
        WRITE(LUPRI,*) 'KFREE            ', KRODBTCH + LRBTCH
        WRITE(LUPRI,*) 'Total memory OD  ', KRODBTCH + LRBTCH - KFREE
        WRITE(LUPRI,*) 'Memory left      ', 
     &   LFREE - ( LICLS + LRCLS + LIBTCH + LRBTCH )
      ENDIF

      CALL DCOPY(2*NODBCLS,WORK(KRODBCLS),1,WORK(KTMP),1)
      KRODBCLS = KIODBCLS + LICLS
      CALL DCOPY(2*NODBCLS,WORK(KTMP),1,WORK(KRODBCLS),1 )

      CALL ICOPY(2*NODBTCH,WORK(KIODBTCH),1,WORK(KTMP),1)
      KIODBTCH = KRODBCLS + LRCLS
      CALL ICOPY(2*NODBTCH,WORK(KTMP),1,WORK(KIODBTCH),1 )

      CALL DCOPY((6+NDMT)*NODBTCH,WORK(KRODBTCH),1,WORK(KTMP),1)
      KRODBTCH = KIODBTCH + LIBTCH
      CALL DCOPY((6+NDMT)*NODBTCH,WORK(KTMP),1,WORK(KRODBTCH),1 )

c     if (.not. (dopbc .and. pbcactive)) then
c        WRITE(LUPRI,*) 'Memory requirement in SETODBTCH', 
c    &        KRODBTCH + LRBTCH - KFREE
c      CALL FLSHFO(LUPRI)
c     end if

      KFREE = KRODBTCH + LRBTCH
      LFREE = LFREE - ( LICLS + LRCLS + LIBTCH + LRBTCH )

C
C     Print section
C     =============

      CALL PRINTODBTCH( WORK(KIODBCLS), WORK(KIODBTCH), 
     &                  IPRINT )

      RETURN 
      END 

C
C ------------------------------------------------------
C
      SUBROUTINE SETODBTCH1(IODBCLS,RODBCLS,IODBTCH,RODBTCH,
     &                      WORK,LWORK,GAB,GABFIT,DMAT,CFIT,
     &                      NDMT,IPRINT,tvecflag)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#include "densfit.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbifck3.h"
C
C aspg, 2006-06-27 making ftnchek happy
C
#include "erithr.h"
C
#include "pbc.h"
      DIMENSION IODBCLS(5,*),RODBCLS(2,*),IODBTCH(2,*),
     &          RODBTCH(6+NDMT,*),
     &          WORK(*), KCCF(2,2),GAB(NBASE,NBASE),GABFIT(NBASEAUX), 
     &          DMAT(NBASE,NBASE), CFIT(NBASEAUX),
     &          tvec1(3), tvec2(3)
      LOGICAL AUX, CSTRESH, tvecflag
C
C     In the construction of the classes 
C     we use the fact that AO batches are  
C     ordered based on the center of the basis 
C     funtion
C     
      DOSCREENCS = ((RHSFIT .OR. LHSFIT) .AND. DOFIT)
     &             .AND. (.NOT. NOSCRE) .AND. .NOT. (RHSFIT.AND.LHSFIT)
cSimen NOSCRE must be put to a commonblock
      IF (IPRINT.GT.35) THEN
        WRITE(LUPRI,*) 'SETODBTCH1 entered'
        WRITE(LUPRI,*) 'DOSCREENCS      ', DOSCREENCS
        WRITE(LUPRI,*) 'LWORK           ', LWORK
        WRITE(LUPRI,*) 'NBASE, NBASEAUX ', NBASE, NBASEAUX
      ENDIF

      IF (DOSCREENCS) 
     &  CALL SET_MAXSCREEN_CS(DMAXDGAB,DMAXDENS,DMAXGAB,
     &                        GAB,GABFIT,DMAT,CFIT,NDMT)

      IF (IPRINT.GT.35) THEN
        WRITE(LUPRI,*) 'Maximum Screening values'
        WRITE(LUPRI,*) 'DMAXGAB  ', DMAXGAB
        WRITE(LUPRI,*) 'DMAXDGAB ', DMAXDGAB
        WRITE(LUPRI,*) 'DMAXDENS ', DMAXDENS
      ENDIF
 
      NPRFMAX = 0 

      ICLS   = 0 
      IBTCH  = 1
      IBTCH0 = 1

      NODCLS11 = 0
      NODCLS03 = 0
      NODBT11 = 0
      NODBT03 = 0

      DO ICLS1 = 1, NAOBCLS
       ICENT1   = IAOCLS(1,ICLS1)
       KSTAOB1  = IAOCLS(2,ICLS1)
       NAOB1    = IAOCLS(3,ICLS1)
       LCLASS1  = IAOCLS(4,ICLS1)
       KENDAOB1 = KSTAOB1 + NAOB1 - 1

       if (dopbc .and. pbcactive) then
          icls2_max = NAOBCLS
       else
          icls2_max = ICLS1
       end if

       DO ICLS2 = 1, icls2_max

        ICENT2   = IAOCLS(1,ICLS2)
        KSTAOB2  = IAOCLS(2,ICLS2)
        NAOB2    = IAOCLS(3,ICLS2)
        LCLASS2  = IAOCLS(4,ICLS2)
        KENDAOB2 = KSTAOB2 + NAOB2 - 1

C       Cauchy-Schwartz screening factor
        CSCLSMAX = 0.d0

        IF ( (LCLASS1.EQ.1.AND.LCLASS2.EQ.1) .OR.
     &       (LCLASS1.EQ.0.AND.LCLASS2.EQ.3) ) THEN

         AUX = LCLASS1.EQ.0.AND.LCLASS2.EQ.3
         if (dopbc .and. pbcactive) then
C           Here we add translations to the two orbitals that
C           form the overlap distribution.
C           Essentially we want to compute J_ab = \sum_cd (ab|cd) D_cd
C           and the translations are different for ab and cd. Thus,
C           we need the flag tvecflag to tell us whether we have ab
C           or cd overlaps.
            if (.not. tvecflag) then
C              We have cd overlap and should add nf_vec to center 1
C              and nf_vec+nf2_vec to center 2.
               tvec1(1) = nf_vec_std(1)
               tvec1(2) = nf_vec_std(2)
               tvec1(3) = nf_vec_std(3)
               tvec2(1) = nf_vec_std(1) + nf2_vec_std(1)
               tvec2(2) = nf_vec_std(2) + nf2_vec_std(2)
               tvec2(3) = nf_vec_std(3) + nf2_vec_std(3)
            else
C              We have ab overlap. Center 1 is always in the central cell
C              and rspc_vec should be added to center 2.
               tvec1(1) = 0.0
               tvec1(2) = 0.0
               tvec1(3) = 0.0
               tvec2(1) = rspc_vec_std(1)
               tvec2(2) = rspc_vec_std(2)
               tvec2(3) = rspc_vec_std(3)
            end if
            CALL GETR12( NSYMOPR, NSYMOPS, ICENT1, ICENT2, 2,
     &                   tvec1, tvec2 )
         else
            tvec1(1) = 0.0
            tvec1(2) = 0.0
            tvec1(3) = 0.0
            tvec2(1) = 0.0
            tvec2(2) = 0.0
            tvec2(3) = 0.0
            CALL GETR12( NSYMOPR, NSYMOPS, ICENT1, ICENT2, 2,
     &                   tvec1, tvec2 )
         end if

         DO NHKT12 = 1, MAXODNHK
          NBTCH = 0 

          SCREENMAXCS  = 0.0d0
          SCREENMAXDCS = 0.0d0
          DO IAOB1 = KSTAOB1, KENDAOB1
           CALL GETAOB(NPRF1,NCTF1,NTYP1,NHKT1,JMAX1,KHKT1,KCKT1,
     &                 KCSQ1,KEXP1,KCCF(1,1),IAOB1)

           if (.not. (dopbc .and. pbcactive)) then
              IF ( ICLS1 .EQ. ICLS2 ) KENDAOB2=IAOB1
C          else
C             ... do nothing ...
C             we require all pairs (IAOB1,IAOB2) for PBC /ErikT
           end if

           DO IAOB2 = KSTAOB2, KENDAOB2
            CALL GETAOB(NPRF2,NCTF2,NTYP2,NHKT2,JMAX2,KHKT2,KCKT2,
     &                  KCSQ2,KEXP2,KCCF(1,2),IAOB2)          
            IF ( NHKT1 + NHKT2 .EQ. NHKT12 + 1 ) THEN  

             CSTRESH = .TRUE.
             IF (DOSCREENCS) THEN
              if (dopbc .and. pbcactive) then
C                Things to check before combining C-S screening and PBC:
C                1. Assumed that matrices are triangular?
C                2. Do we have the PBC density matrix and not just, say,
C                   the D^{00} component?
C                3. We screen out an infinite series and not just a few
C                   terms.
C                4. Do we always send in real GAB and not just a
C                   dummy? / ErikT
                 write(LUPRI,*) 'SETODBTCH1: Cauchy-Schwarz screening',
     &                ' for PBC need to be thought through.'
                 call quit('SETODBTCH1: Cauchy-Schwarz used with PBC.')
              end if
              CALL SETMAXCSBTCH(RODBTCH(1,IBTCH), NPRF1, NCTF1, NTYP1, 
     &              NPRF2, NCTF2, NTYP2, KAOBTCH(1,IAOB1), 
     &              KAOBTCH(1,IAOB2), KHKT1, KHKT2,
     &              GAB, GABFIT, DMAT, CFIT, NDMT, AUX, IPRINT)
              IF (RHSFIT .AND. LHSFIT) THEN
                IF (AUX) THEN
                  THRESHOLD = THRCS/(MAX(DMAXDENS*DMAXGAB,DMAXDGAB))
                ELSE
                  THRESHOLD = 0.0d0
                ENDIF
              ELSE IF (RHSFIT) THEN
                IF (AUX) THEN
                  THRESHOLD = THRCS/DMAXDENS/DMAXGAB
                ELSE
                  THRESHOLD = THRCS/DMAXDGAB
                ENDIF
              ELSE IF (LHSFIT) THEN
                IF (AUX) THEN
                  THRESHOLD = THRCS/DMAXDGAB
                ELSE
                  THRESHOLD = THRCS/DMAXDENS/DMAXGAB
                ENDIF
              ELSE
                WRITE(LUPRI,*) 'Error in setting up OD-batches with ',
     &                         'DOSCREENCS = .TRUE.'
                CALL QUIT('Either RHSFIT or LHSFIT must be true')
              ENDIF

              SCREENCS     = RODBTCH(6,IBTCH)
              CSTRESH      = SCREENCS .GT. THRESHOLD .OR. NCSODB

              SCREENMAXCS  = MAX(SCREENMAXCS, SCREENCS)
              SCREENMAXDCS = MAX(SCREENMAXDCS,RODBTCH(7,IBTCH))
             ENDIF

             IF (CSTRESH) THEN
C      PBC comment: SCREENPRIM12 uses coordinates and a distance set
C        by GETR12, so the necessary translations have already
C        been put in. Hence, SCREENPRIM12 works without modification for
C        PBC. One thing, though: In the molecular case, one screens out
C        just a few terms, but in the PBC case screening means neglecting
C        a rapidly converging infinite series. This means that the error
C        estimates may be a little different. / ErikT
               CALL SCREENPRIM12(NPRF12,RODBTCH(1,IBTCH),IPR12,.FALSE.,
     &              NPRF1,NCTF1,NTYP1,NPRF2,NCTF2,NTYP2,
     &              EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &              CCFBTCH(KCCF(1,1)),CCFBTCH(KCCF(2,1)),
     &              CCFBTCH(KCCF(1,2)),CCFBTCH(KCCF(2,2)),
     &              R122(0,0),COOR1(1,0),COOR2(1,0),
     &              NDMT,IPRINT)

               IF ( NPRF12 .GT. 0 ) THEN                         
                 NBTCH = NBTCH + 1 
C Counts number of OD-batches of normal (11) and 
C auxilliary (03) type
                 IF (LCLASS1.EQ.1)  NODBT11 = NODBT11 + 1
                 IF (LCLASS1.EQ.0)  NODBT03 = NODBT03 + 1
                 IODBTCH(1,IBTCH) = IAOB1
                 IODBTCH(2,IBTCH) = IAOB2
                 IBTCH = IBTCH + 1
                 NPRFMAX = MAX( NPRFMAX, NPRF12 ) 
                 CSCLSMAX = MAX( CSCLSMAX, RODBTCH(5,IBTCH) )
               END IF
             ENDIF
            END IF

           END DO
          END DO

          IF ( NBTCH .GT. 0 ) THEN 
           ICLS = ICLS + 1 
           IODBCLS(1,ICLS) = ICENT1
           IODBCLS(2,ICLS) = ICENT2
           IODBCLS(3,ICLS) = NHKT12
           IODBCLS(4,ICLS) = IBTCH0
           IODBCLS(5,ICLS) = NBTCH

           IF (DOSCREENCS) THEN
             RODBCLS(1,ICLS) = SCREENMAXCS
             RODBCLS(2,ICLS) = SCREENMAXDCS
           ENDIF
           IBTCH0 = IBTCH
C          Counts number of OD-classes of normal (11) and 
C          auxilliary (03) type
           IF (LCLASS1.EQ.1)  NODCLS11 = NODCLS11 + 1
           IF (LCLASS1.EQ.0)  NODCLS03 = NODCLS03 + 1
          END IF

         END DO

        END IF !Test of LCLASS1, LCLASS2

       END DO
      END DO

      NODBCLS = ICLS 
      NODBTCH = IBTCH - 1

      NODCLSNORM = NODCLS11
      NODCLSAUX = NODCLS03

      NODBTNORM = NODBT11
      NODBTAUX = NODBT03
      RETURN
      END

C
C -----------------------------------------------------
C

      SUBROUTINE GETMAXODBCLS(MAXODBCLS,MAXODBTCH,IPRINT,tvecflag)
#include "implicit.h"
#include "priunit.h" 
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#include "pi.h"
#include "densfit.h"
C
#include "pbc.h"
      dimension tvec1(3), tvec2(3)
      logical tvecflag

      IF (IPRINT.GT.5)
     &  CALL HEADER( 'OUTPUT FROM GETMAXODBCLS', -1 )

      EXPMUMIN = EXPAOMIN/2.0d0
      COEFMAX = R2PI52 / EXPAOMIN * CCFAOMAX * CCFAOMAX 
      R2AOCNTMAX  =  ABS( LOG( THRCOD/COEFMAX ) )/EXPMUMIN

      MAXODBTCH = 0
      MAXODBCLS = 0

C
C aspg, 2006-06-27 making ftnchek happy
C
      NSYMOPR = 0
      NSYMOPS = 0
C

      IF (IPRINT.GT.5) THEN
        WRITE(LUPRI,'(A,3F14.6,/,A,F14.6,2I12,/,A,I12,F14.6,/,A,F14.6)') 
     &       'R2PI52, EXPAOMIN, CCFAOMAX   : ',
     &       R2PI52, EXPAOMIN, CCFAOMAX,
     &       'THRCOD, NAOBCLS, NSYMOPR     :  ',
     &       THRCOD, NAOBCLS, NSYMOPR,
     &       'NSYMOPS,  R2AOCNTMAX         :  ',
     &       NSYMOPS, R2AOCNTMAX,
     &       'COEFMAX                      :  ',
     &       COEFMAX
      ENDIF

      DO ICLS1 = 1, NAOBCLS
       ICENT1   = IAOCLS(1,ICLS1)
       NAOB1    = IAOCLS(3,ICLS1)
       LCLASS1  = IAOCLS(4,ICLS1)
       if (dopbc .and. pbcactive) then
          icls2_max = NAOBCLS
       else
          icls2_max = ICLS1
       end if
       DO ICLS2 = 1, icls2_max
        ICENT2   = IAOCLS(1,ICLS2)
        NAOB2    = IAOCLS(3,ICLS2)
        LCLASS2  = IAOCLS(4,ICLS2)

        IF (LCLASS1.EQ.1.AND.LCLASS2.EQ.1) THEN
           if (dopbc .and. pbcactive) then
            if (.not. tvecflag) then
C              We have cd overlap and should add nf_vec to center 1
C              and nf_vec+nf2_vec to center 2.
               tvec1(1) = nf_vec_std(1)
               tvec1(2) = nf_vec_std(2)
               tvec1(3) = nf_vec_std(3)
               tvec2(1) = nf_vec_std(1) + nf2_vec_std(1)
               tvec2(2) = nf_vec_std(2) + nf2_vec_std(2)
               tvec2(3) = nf_vec_std(3) + nf2_vec_std(3)
            else
C              We have ab overlap. Center 1 is always in the central cell
C              and rspc_vec should be added to center 2.
               tvec1(1) = 0.0
               tvec1(2) = 0.0
               tvec1(3) = 0.0
               tvec2(1) = rspc_vec_std(1)
               tvec2(2) = rspc_vec_std(2)
               tvec2(3) = rspc_vec_std(3)
            end if
           else
              tvec1(1) = 0.0
              tvec1(2) = 0.0
              tvec1(3) = 0.0
              tvec2(1) = 0.0
              tvec2(2) = 0.0
              tvec2(3) = 0.0
           end if
           CALL GETR12( NSYMOPR, NSYMOPS, ICENT1, ICENT2, 2,
     &                  tvec1, tvec2 )

          IF ( R122(0,0) .LE. R2AOCNTMAX ) THEN
cSimen      Should differentiate MAXAONHK for normal- and aux-basis
            MAXODBCLS = MAXODBCLS + ( 2*MAXAONHK - 1 )
            MAXODBTCH = MAXODBTCH + NAOB1*NAOB2 
          END IF
        ELSEIF (LCLASS1.EQ.0.AND.LCLASS2.EQ.3) THEN
          MAXODBCLS = MAXODBCLS + ( 2*MAXAONHK - 1 )
          MAXODBTCH = MAXODBTCH + NAOB1*NAOB2
        END IF
       END DO
      END DO

      IF (IPRINT.GT.5) THEN
        WRITE(LUPRI,'(A,3F14.6,/,A,2I12)') 
     &       'RAOCNTMAX, EXPAOMIN, COEFMAX : ',
     &       SQRT(R2AOCNTMAX), EXPAOMIN, COEFMAX,
     &       'MAXODBCLS, MAXODBTCH           ',
     &       MAXODBCLS, MAXODBTCH   
      ENDIF

      RETURN
      END

C
C -----------------------------------------------------
C

      SUBROUTINE PRINTODBTCH( IODBCLS, IODBTCH, IPRINT )
#include "implicit.h"
#include "priunit.h" 
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3odb.h"

      DIMENSION IODBCLS( 5, * ), IODBTCH( 2, * )

      IF ( IPRINT .GT. 2 ) THEN 

       CALL HEADER( 'OUTPUT FROM SETODBTCH', -1 )

       WRITE(LUPRI,'(2(3X,A,I9,/))') 
     &    'Number of ODB classes  : ', NODBCLS,
     &    'Total number of ODBTCHs: ', NODBTCH

       DO ICLS = 1, NODBCLS
        KBTCH = IODBCLS(4,ICLS)
        NBTCH = IODBCLS(5,ICLS)
        WRITE(LUPRI,'(3(3X,A,I5),3X,A,2I3)')
     &        'Class number : ', ICLS,
     &        'NHKT12 : ', IODBCLS(3,ICLS),
     &        '#ODBTCH: ', NBTCH,
     &        'Centers :', 
     &        IODBCLS(1,ICLS), IODBCLS(2,ICLS)

 
       IF ( IPRINT .GT. 4 ) THEN 

        DO IBTCH = KBTCH, (KBTCH+NBTCH-1)
         WRITE( LUPRI, '(I10,2I6)') IBTCH,
     &          IODBTCH(1,IBTCH),IODBTCH(2,IBTCH)
        END DO
       END IF 

       END DO
      END IF 

      RETURN
      END

C
C -----------------------------------------------------
C

      SUBROUTINE SETMAXCSBTCH(RODBTCH, NPRF1, NCTF1, NTYP1, 
     &     NPRF2, NCTF2, NTYP2, IAO1, IAO2, KHKT1, KHKT2,
     &     GAB, GABFIT, DMAT, CFIT, NDMT, AUX, IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "fck3thr.h"
#include "pi.h"
#include "densfit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
      DIMENSION RODBTCH((6+NDMT)) ,GAB(NBASE,NBASE), GABFIT(NBASEAUX), 
     &          DMAT(NBASE,NBASE), CFIT(NBASEAUX), IAO1(2), IAO2(2)
      LOGICAL AUX

C     *******************************************************
C     *          Picks up the largest GAB element         *
C     *          and largest Dab x GAB element            *
C     *******************************************************

      IF (IPRINT.GT.45) THEN
        WRITE(LUPRI,*) '*****Output from SETMAXCSBTCH*****'
        WRITE(LUPRI,*) 'NPRF1, NCTF1, NTYP1 ', NPRF1, NCTF1, NTYP1
        WRITE(LUPRI,*) 'NPRF2, NCTF2, NTYP2 ', NPRF2, NCTF2, NTYP2
        WRITE(LUPRI,*) 'IAO1, IAO2 (', IAO1(1), IAO1(2), '), (', 
     &                  IAO2(1), IAO2(2), ')'
        WRITE(LUPRI,*) 'KHKT1, KHKT2, AUX   ', KHKT1, KHKT2, AUX
        CALL FLSHFO(LUPRI)
      ENDIF
      IF (DOSCREENCS) THEN
        GABMAX  = 0.0d0
        DGABMAX = 0.0d0
C       *** Auilliary batch ***
        IF (AUX) THEN
          IF (LHSFIT .AND.(.NOT.RHSFIT)) THEN
C           Shared sp-exponents
            IF (NTYP2 .EQ. 2) THEN
              IORBSTART = IAO2(1) - NBASE
              DO IORB2 = IORBSTART, IORBSTART + NCTF2 - 1
                GABMAX  = MAX( GABMAX,  GABFIT (IORB2) )
              ENDDO
              IORBSTART = IAO2(2) - NBASE
              DO IORB2 = IORBSTART, IORBSTART + NCTF2*3 - 1
                GABMAX  = MAX( GABMAX,  GABFIT (IORB2) )
              ENDDO
C           Not shared sp-exponents
            ELSE
              IORBSTART = IAO2(1) - NBASE
              DO IORB2 = IORBSTART, IORBSTART + NCTF2*KHKT2 - 1
                GABMAX  = MAX( GABMAX,  GABFIT (IORB2) )
              ENDDO
            ENDIF
          ELSE ! IF RHSFIT
C           Shared sp-exponents
            IF (NTYP2 .EQ. 2) THEN
              IORBSTART = IAO2(1) - NBASE
              DO IORB2 = IORBSTART, IORBSTART + NCTF2 - 1
                GABMAX  = MAX( GABMAX,  GABFIT (IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(CFIT(IORB2)) 
     &                                 * GABFIT (IORB2) )
              ENDDO
              IORBSTART = IAO2(2) - NBASE
              DO IORB2 = IORBSTART, IORBSTART + NCTF2*3 - 1
                GABMAX  = MAX( GABMAX,  GABFIT (IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(CFIT(IORB2)) 
     &                                 * GABFIT (IORB2) )
              ENDDO
C           Not shared sp-exponents
            ELSE
              IORBSTART = IAO2(1) - NBASE
              DO IORB2 = IORBSTART, IORBSTART + NCTF2*KHKT2 - 1
                GABMAX  = MAX( GABMAX,  GABFIT (IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(CFIT(IORB2)) 
     &                                 * GABFIT (IORB2) )
              ENDDO
            ENDIF
          ENDIF
C       *** Normal batch ***
        ELSE
C         Both AO-batches share sp-exponents
          IF (NTYP1 .EQ. 2 .AND. NTYP2 .EQ. 2) THEN
            DO IORB1 = IAO1(1), IAO1(1) + NCTF1 - 1
              DO IORB2 = IAO2(1), IAO2(1) + NCTF2 - 1
                GABMAX  = MAX( GABMAX,  GAB (IORB1,IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(DMAT(IORB1,IORB2) )
     &                                 * GAB (IORB1,IORB2) )
                IF (.NOT.SYMDENS) DGABMAX = 
     &            MAX(DGABMAX,ABS(DMAT(IORB2,IORB1))*GAB(IORB1,IORB2))
              ENDDO
              DO IORB2 = IAO2(2), IAO2(2) + NCTF2*3 - 1
                GABMAX  = MAX( GABMAX,  GAB (IORB1,IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(DMAT(IORB1,IORB2) )
     &                                 * GAB (IORB1,IORB2) )
                IF (.NOT.SYMDENS) DGABMAX = 
     &            MAX(DGABMAX,ABS(DMAT(IORB2,IORB1))*GAB(IORB1,IORB2))
              ENDDO
            ENDDO
            DO IORB1 = IAO1(2), IAO1(2) + NCTF1*3 - 1
              DO IORB2 = IAO2(1), IAO2(1) + NCTF2 - 1
                GABMAX  = MAX( GABMAX,  GAB (IORB1,IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(DMAT(IORB1,IORB2) )
     &                                 * GAB (IORB1,IORB2) )
                IF (.NOT.SYMDENS) DGABMAX = 
     &            MAX(DGABMAX,ABS(DMAT(IORB2,IORB1))*GAB(IORB1,IORB2))
              ENDDO
              DO IORB2 = IAO2(2), IAO2(2) + NCTF2*3 - 1
                GABMAX  = MAX( GABMAX,  GAB (IORB1,IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(DMAT(IORB1,IORB2) )
     &                                 * GAB (IORB1,IORB2) )
                IF (.NOT.SYMDENS) DGABMAX = 
     &            MAX(DGABMAX,ABS(DMAT(IORB2,IORB1))*GAB(IORB1,IORB2))
              ENDDO
            ENDDO
C       First AO-batch share sp-exponents
          ELSE IF (NTYP1 .EQ. 2) THEN
            DO IORB1 = IAO1(1), IAO1(1) + NCTF1 - 1
              DO IORB2 = IAO2(1), IAO2(1) + NCTF2*KHKT2 - 1
                GABMAX  = MAX( GABMAX,  GAB (IORB1,IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(DMAT(IORB1,IORB2) )
     &                                 * GAB (IORB1,IORB2) )
                IF (.NOT.SYMDENS) DGABMAX = 
     &            MAX(DGABMAX,ABS(DMAT(IORB2,IORB1))*GAB(IORB1,IORB2))
              ENDDO
            ENDDO
            DO IORB1 = IAO1(2), IAO1(2) + NCTF1*3 - 1
              DO IORB2 = IAO2(1), IAO2(1) + NCTF2*KHKT2 - 1
                GABMAX  = MAX( GABMAX,  GAB (IORB1,IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(DMAT(IORB1,IORB2) )
     &                                 * GAB (IORB1,IORB2) )
                IF (.NOT.SYMDENS) DGABMAX = 
     &            MAX(DGABMAX,ABS(DMAT(IORB2,IORB1))*GAB(IORB1,IORB2))
              ENDDO
            ENDDO
C     Second AO-batch share sp-exponents
          ELSE IF (NTYP2 .EQ. 2) THEN
            DO IORB1 = IAO1(1), IAO1(1) + NCTF1*KHKT1 - 1
              DO IORB2 = IAO2(1), IAO2(1) + NCTF2 - 1
                GABMAX  = MAX( GABMAX,  GAB (IORB1,IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(DMAT(IORB1,IORB2) )
     &                                 * GAB (IORB1,IORB2) )
                IF (.NOT.SYMDENS) DGABMAX = 
     &            MAX(DGABMAX,ABS(DMAT(IORB2,IORB1))*GAB(IORB1,IORB2))
              ENDDO
              DO IORB2 = IAO2(2), IAO2(2) + NCTF2*3 - 1
                GABMAX  = MAX( GABMAX,  GAB (IORB1,IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(DMAT(IORB1,IORB2) )
     &                                 * GAB (IORB1,IORB2) )
                IF (.NOT.SYMDENS) DGABMAX = 
     &            MAX(DGABMAX,ABS(DMAT(IORB2,IORB1))*GAB(IORB1,IORB2))
              ENDDO
            ENDDO
C     No AO-batches share sp-exponents
          ELSE
            DO IORB1 = IAO1(1), IAO1(1) + NCTF1*KHKT1 - 1
              DO IORB2 = IAO2(1), IAO2(1) + NCTF2*KHKT2 - 1
                GABMAX  = MAX( GABMAX,  GAB (IORB1,IORB2) )
                DGABMAX = MAX( DGABMAX, ABS(DMAT(IORB1,IORB2) )
     &                                 * GAB (IORB1,IORB2) )
                IF (.NOT.SYMDENS) DGABMAX = 
     &            MAX(DGABMAX,ABS(DMAT(IORB2,IORB1))*GAB(IORB1,IORB2))
              ENDDO
            ENDDO
          ENDIF ! NTYP
        ENDIF ! AUX
      ENDIF ! DOSCREENCS

      RODBTCH(6) = GABMAX
      RODBTCH(7) = DGABMAX

      RETURN
      END 


      SUBROUTINE SCREENPRIM12(NPRF12,RODBTCH,IPR12,SETIPR,
     &                        NPRF1,NCTF1,NTYP1,NPRF2,NCTF2,NTYP2,
     &                        EXPBT1,EXPBT2,CCF11,CCF21,CCF12,CCF22,
     &                        DIST122,CORD1,CORD2,NDMT,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "fck3thr.h"
#include "pi.h"
#include "densfit.h"
      DIMENSION IPR12(2,NPRF1*NPRF2),RODBTCH(6+NDMT),
     &          CORD1(3),CORD2(3),EXPBT1(NPRF1),EXPBT2(NPRF2),
     &          CCF11(NPRF1,NCTF1),CCF21(NPRF1,NCTF1),
     &          CCF12(NPRF2,NCTF2),CCF22(NPRF2,NCTF2)


      LOGICAL SETIPR

      SM=0.0D0
      AM=0.0D0
      BM=0.0D0
      EXTPMAX=0.0D0
      FACSCREEN=0.0D0

      IPRF12 = 0 
      DO IPRF1=1, NPRF1
       EXP1   = EXPBT1(IPRF1)
        
       DO IPRF2=1, NPRF2
        EXP2   = EXPBT2(IPRF2)
        EXP12  = EXP1 + EXP2 
         
        EXP12I = 1.0D0 / EXP12
        EXPMU  = EXP1 * EXP2 / EXP12
          
        CCFSUM11=0.0D0
        CCFSUM21=0.0D0
        CCFSUM12=0.0D0
        CCFSUM22=0.0D0
          
        DO ICTF1 = 1, NCTF1
         CCFSUM11 = CCFSUM11 + ABS(CCF11(IPRF1,ICTF1))
        END DO
          
        DO ICTF2 = 1, NCTF2
         CCFSUM12 = CCFSUM12 + ABS(CCF12(IPRF2,ICTF2))
        END DO
          
        IF ( NTYP1 .EQ. 2 ) THEN 
         DO ICTF1 = 1, NCTF1
          CCFSUM21 = CCFSUM21 + ABS(CCF21(IPRF1,ICTF1))
         END DO
        END IF
          
        IF ( NTYP2 .EQ. 2 ) THEN 
         DO ICTF2 = 1, NCTF2
          CCFSUM22 = CCFSUM22 + ABS(CCF22(IPRF2,ICTF2))
         END DO
        END IF
          
        CCFSUM = MAX(CCFSUM11,CCFSUM21) * MAX(CCFSUM12,CCFSUM22)


C     
C     This factor for screening is a very  simple
C     approach based on exponencial prefactor.
C     It is not rigorous, and the best would be
C     to put GAB for real Cauchy-Schwarz screening
C     
            
        FAC = R2PI52*EXP(-EXPMU*DIST122)*EXP12I * CCFSUM
C     
C     What is better way of setting the Cauchy-Schwarz -
C     making sum of FAC, or just maximum ???
C     
C     FACSCREEN = FACSCREEN + FAC
      
c        IF ( FAC .GT. THRPOD .OR. (RHSFIT.OR.LHSFIT) ) THEN

        IF ( FAC .GT. THRPOD ) THEN
            
         FACSCREEN = MAX(FACSCREEN,FAC)
            
         WT = ( CCFSUM11+CCFSUM21 ) * ( CCFSUM12 + CCFSUM22 )
            
         IPRF12 = IPRF12 + 1
            
         IF ( SETIPR ) THEN
          IPR12(1,IPRF12) = IPRF1
          IPR12(2,IPRF12) = IPRF2
         END IF
            
         WT = WT * EXP(-EXPMU*DIST122) 
         SM = SM + WT
         AM = AM + WT * EXP1 * EXP12I
         BM = BM + WT * EXP2 * EXP12I
            
         EXTP = SQRT(EXP12I) * FACCLS
         EXTPMAX = MAX( EXTPMAX,EXTP )
            
        END IF
          
       END DO
      END DO
      
      NPRF12 = IPRF12
      
      IF ( IPRF12 .GT. 0) THEN
       AM = AM/SM
       BM = BM/SM
        
       RODBTCH(1) = AM * CORD1(1) + BM * CORD2(1)
       RODBTCH(2) = AM * CORD1(2) + BM * CORD2(2)
       RODBTCH(3) = AM * CORD1(3) + BM * CORD2(3)
       RODBTCH(4) = EXTPMAX
       RODBTCH(5) = FACSCREEN
          
      END IF

      RETURN
      END 

C
C -----------------------------------------------------
C

      SUBROUTINE GETR12( NSYMOP1, NSYMOP2, 
     &                   ICENT1, ICENT2, NACTCENT,
     &                   TRANSL1, TRANSL2 )
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "nuclei.h"
      dimension TRANSL1(*), TRANSL2(*)

      IF (ICENT1.EQ.0) THEN
        COOR1(1,0) = 0.0d0
        COOR1(2,0) = 0.0d0
        COOR1(3,0) = 0.0d0
      ELSE
        COOR1(1,0) = CORD(1,ICENT1) + TRANSL1(1)
        COOR1(2,0) = CORD(2,ICENT1) + TRANSL1(2)
        COOR1(3,0) = CORD(3,ICENT1) + TRANSL1(3)
      ENDIF

      COOR2(1,0) = CORD(1,ICENT2) + TRANSL2(1)
      COOR2(2,0) = CORD(2,ICENT2) + TRANSL2(2)
      COOR2(3,0) = CORD(3,ICENT2) + TRANSL2(3)

      R12(1,0,0) = COOR1(1,0) - COOR2(1,0)
      R12(2,0,0) = COOR1(2,0) - COOR2(2,0)
      R12(3,0,0) = COOR1(3,0) - COOR2(3,0)

      R122(0,0)  = R12(1,0,0)**2+R12(2,0,0)**2+R12(3,0,0)**2 

      RETURN
      END 

C
C -----------------------------------------------------
C

      SUBROUTINE PUT2ERICOM(MAXODNHK,MAXODKCK,NCSQ1VAL,NCSQ2VAL)
#include "implicit.h"
#include "ericom.h"

      NHKMAX = MAXODNHK
      KCKMAX = MAXODKCK
      NCSQ1  = NCSQ1VAL
      NCSQ2  = NCSQ2VAL

      RETURN
      END

      SUBROUTINE GET_ODPAR(IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3rel.h"
#include "fck3thr.h"
#include "clsfmm.h"
#include "priunit.h"
#include "maxaqn.h"
      MAXODNHK = 0
      MAXODPRF = 0
      MAXODCTF = 0
      DO ICLS1 = 1, NAOBCLS
        KSTAOB1  = IAOCLS(2,ICLS1)
        NAOB1    = IAOCLS(3,ICLS1)
        LCLASS1  = IAOCLS(4,ICLS1)
        KENDAOB1 = KSTAOB1 + NAOB1 - 1
        DO ICLS2 = 1, ICLS1
          KSTAOB2  = IAOCLS(2,ICLS2)
          NAOB2    = IAOCLS(3,ICLS2)
          LCLASS2  = IAOCLS(4,ICLS2)
          KENDAOB2 = KSTAOB2 + NAOB2 - 1
          IF ( (LCLASS1.EQ.1.AND.LCLASS2.EQ.1) .OR.
     &         (LCLASS1.EQ.0.AND.LCLASS2.EQ.3) ) THEN
            DO IAOB1 = KSTAOB1, KENDAOB1
              NPRF1  = NPRFBTCH( IAOB1 )
              NCTF1  = NCTFBTCH( IAOB1 )
              NHKT1  = NHKTBTCH( IAOB1 )
              DO IAOB2 = KSTAOB2, KENDAOB2
                NPRF2  = NPRFBTCH( IAOB2 )
                NCTF2  = NCTFBTCH( IAOB2 )
                NHKT2  = NHKTBTCH( IAOB2 )
                MAXODNHK = MAX(MAXODNHK, NHKT1 + NHKT2 - 1)
                MAXODPRF = MAX(MAXODPRF, NPRF1 * NPRF2)
                MAXODCTF = MAX(MAXODCTF, NCTF1 * NCTF2)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      MAXODNHK = MAXODNHK + 2*NDER
      IF (MAXODNHK.GT.MXQN) THEN
        WRITE(LUPRI,'(1X,A,I2,A,I2,A)') 'MAXODNHK =',MAXODNHK, 
     &    ' is greater than MXQN =', MXQN, ' in GET_ODPAR'
        CALL QUIT('Maximum angular momentum greater than MXQN')
      ENDIF
      MAXODKCK = MAXAOKCK * MAXAOKCK
      MAXODKHK = MAXAOKHK * MAXAOKHK

      NCSQ1 = NSPHER(0,1,0,0) - 1
      NCSQ2 = 1  

      FACCLS = ERFCIV(THRCLS)

      CALL PUT2ERICOM( MAXODNHK, MAXODKCK, NCSQ1, NCSQ2 )

      IF (IPRINT .GT. 25) THEN
        WRITE (LUPRI,'(7X,A,2I5)')'MAXODNHK,MAXODPRF ',MAXODNHK,MAXODPRF
        WRITE (LUPRI,'(7X,A, I5)')'MAXODCTF          ',MAXODCTF
      ENDIF
      RETURN
      END

      SUBROUTINE SET_MAXSCREEN_CS(DMAXDGAB,DMAXDENS,DMAXGAB,
     &                            GAB,GABFIT,DMAT,CFIT,NDMT)
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "densfit.h"
      DIMENSION GAB(NBASE,NBASE),GABFIT(NBASEAUX), 
     &          DMAT(NBASE,NBASE,NDMT), CFIT(NBASEAUX,NDMT)
      DMAXDGAB = 0.0d0
      DMAXDENS = 0.0d0
      DMAXGAB  = 0.0d0
      DO IDMT=1,NDMT
        IF (LHSFIT .AND. RHSFIT) THEN
          DO I=1,NBASEAUX
            DMAXGAB  = MAX(DMAXGAB,  GABFIT(I))
            DMAXDGAB = MAX(DMAXDGAB, GABFIT(I)*ABS(CFIT(I,IDMT)) )
            DMAXDENS = MAX(DMAXDENS, ABS(CFIT(I,IDMT)) )
          ENDDO
        ELSE IF (LHSFIT) THEN
          DO I=1,NBASE
            DO J=1,NBASE
              DMAXDGAB = MAX(DMAXDGAB, GAB(I,J)*ABS(DMAT(I,J,IDMT)))
              DMAXDENS = MAX(DMAXDENS, ABS(DMAT(I,J,IDMT)) )
            ENDDO
          ENDDO
          DO I=1,NBASEAUX
            DMAXGAB  = MAX(DMAXGAB,  GABFIT(I))
          ENDDO
        ELSE IF (RHSFIT) THEN
          DO I=1,NBASE
            DO J=1,NBASE
              DMAXGAB = MAX(DMAXGAB, GAB(I,J)  )
            ENDDO
          ENDDO
          DO I=1,NBASEAUX
            DMAXDGAB = MAX(DMAXDGAB, GABFIT(I)*ABS(CFIT(I,IDMT)) )
            DMAXDENS = MAX(DMAXDENS, ABS(CFIT(I,IDMT)) )
          ENDDO
        ELSE
          DO I=1,NBASEAUX
            DMAXGAB  = MAX(DMAXGAB,  GABFIT(I))
            DMAXDGAB = MAX(DMAXDGAB, GABFIT(I)*ABS(CFIT(I,IDMT)) )
            DMAXDENS = MAX(DMAXDENS, ABS(CFIT(I,IDMT)) )
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
