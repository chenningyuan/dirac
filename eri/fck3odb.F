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
C    Overlap Distribution classes :
C
C    IODBCLS(5,*) :
C
C      1 -  ICENT1
C      2 -  ICENT2
C      3 -  NHKT12 ( = J1 + J2 + 1; J = ang. momentum )
C      4 -  KSTBTCH - number of the first OD batch of 
C           this class ( used in IODBTCH )
C      5 -  NBTCH   - number of OD batches in this class
C
C    --------------------------------------------------
C
C    Overlap Distribution batches :
C
C    IODBTCH(2,*) :
C
C      1 - IAO1 - First AO batch
C      2 - IAO2 - Second AO batch
C
C    ( Note that AO batch is not the number of atomic 
C      orbital ( that is accessed via INDXBTCH array 
C      from fck3aob.h common block )  ! )
C  
C    RODBTCH(5,*) 
C
C      1 - x-coodinate of the center
C      2 - y-coodinate of the center
C      3 - z-coodinate of the center
C      4 - non-classical extent (from the center)
C      5 - Cauchy screening parameter

#if defined (VAR_SINGLE)
#include "single.h"
#endif
      SUBROUTINE SETODBTCH(KIODBCLS,KIODBTCH,KRODBTCH,
     &                KFREE,LFREE,WORK,LWORK,IPRINT)
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
#include "pi.h"

C
C  MAXODNHK is the JMAX12+1 for overlap distribution
C  while MAXAONHK is JMAX1+1 for the atomic orbitals
C

      MAXODNHK = 2 * MAXAONHK - 1
      MAXODKCK =     MAXODNHK * ( MAXODNHK + 1 ) / 2
      MAXODKHK = 2 * MAXODNHK -  1
      MAXODPRF = MAXAOPRF * MAXAOPRF
      MAXODCTF = MAXAOCTF * MAXAOCTF

      NCSQ1 = NSPHER(0,1,0,0) - 1
      NCSQ2 = 1  

      FACCLS = ERFCIV(THRCLS)

C
C    Assign these variables also to ericom common block
C   

      CALL PUT2ERICOM( MAXODNHK, MAXODKCK, NCSQ1, NCSQ2 )

      CALL GETMAXODBCLS(MAXODBCLS,MAXODBTCH,IPRINT)

      KIODBCLS = KFREE
      KIODBTCH = KIODBCLS + 5*MAXODBCLS/IRAT + 1
      KRODBTCH = KIODBTCH + 2*MAXODBTCH/IRAT + 1 
      KLAST    = KRODBTCH + 5*MAXODBTCH

      LWRK     = LFREE - ( KLAST-KFREE) 

      IF ( LWRK .LT. 0 ) THEN 
        CALL STOPIT('SETODBTCH','FCK3ODB',KLAST,LWORK)
      END IF 

      CALL SETODBTCH1(WORK(KIODBCLS),WORK(KIODBTCH),WORK(KRODBTCH),
     &                WORK(KLAST),LWRK,IPRINT)

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
      LIBTCH = 2*NODBTCH/IRAT + 1
      LRBTCH = 5*NODBTCH

      LTMP  = MAX(LIBTCH,LRBTCH)

      KTMP  = KLAST 
      KLAST = KTMP  + LTMP

      IF ( LTMP .GT. LWRK ) THEN 
        CALL STOPIT('SETODBTCH','ODBTCHTMP',KLAST,LWORK)
      END IF 

C
C     copy the batch information to the temporary place 
C     and than to its final place.
C  ( tmp copying just to avoind any possible problems
C     with copying array to  overlaping place )    
C
   
      CALL ICOPY(2*NODBTCH,WORK(KIODBTCH),1,WORK(KTMP),1)
      KIODBTCH = KIODBCLS + LICLS
      CALL ICOPY(2*NODBTCH,WORK(KTMP),1,WORK(KIODBTCH),1 )

      CALL DCOPY(5*NODBTCH,WORK(KRODBTCH),1,WORK(KTMP),1)
      KRODBTCH = KIODBTCH + LIBTCH
      CALL DCOPY(5*NODBTCH,WORK(KTMP),1,WORK(KRODBTCH),1 )

      KFREE = KRODBTCH + LRBTCH
      LFREE = LFREE - ( LICLS + LIBTCH + LRBTCH )

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
      SUBROUTINE SETODBTCH1(IODBCLS,IODBTCH,RODBTCH,
     &                      WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
      DIMENSION IODBCLS(5,*),IODBTCH(2,*),RODBTCH(5,*),WORK(*),
     &          KCCF(2,2)
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
C
C     In the construction of the classes 
C     we use the fact that AO batches are  
C     ordered based on the center of the basis 
C     funtion
C      
      NPRFMAX = 0 

      ICLS   = 0 
      IBTCH  = 1
      IBTCH0 = 1

      DO ICLS1 = 1, NAOBCLS
       ICENT1   = IAOCLS(1,ICLS1)
       KSTAOB1  = IAOCLS(2,ICLS1)
       NAOB1    = IAOCLS(3,ICLS1)
       KENDAOB1 = KSTAOB1 + NAOB1 - 1

       DO ICLS2 = 1, ICLS1

        ICENT2   = IAOCLS(1,ICLS2)
        KSTAOB2  = IAOCLS(2,ICLS2)
        NAOB2    = IAOCLS(3,ICLS2)
        KENDAOB2 = KSTAOB2 + NAOB2 - 1

        CALL GETR12( NSYMOPR, NSYMOPS, ICENT1, ICENT2, 2 )

        DO NHKT12 = 1, MAXODNHK
         NBTCH = 0 

         DO IAOB1 = KSTAOB1, KENDAOB1
          CALL GETAOB(NPRF1,NCTF1,NTYP1,NHKT1,JMAX1,KHKT1,KCKT1,
     &                KCSQ1,KEXP1,KCCF(1,1),IAOB1)

          IF ( ICLS1 .EQ. ICLS2 ) KENDAOB2=IAOB1

          DO IAOB2 = KSTAOB2, KENDAOB2
           CALL GETAOB(NPRF2,NCTF2,NTYP2,NHKT2,JMAX2,KHKT2,KCKT2,
     &                 KCSQ2,KEXP2,KCCF(1,2),IAOB2)          

           IF ( NHKT1 + NHKT2 .EQ. NHKT12 + 1 ) THEN  

            CALL SCREENPRIM12(NPRF12,RODBTCH(1,IBTCH),IPR12,.FALSE.,
     &                       NPRF1,NCTF1,NTYP1,NPRF2,NCTF2,NTYP2,
     &                       EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &                       CCFBTCH(KCCF(1,1)),CCFBTCH(KCCF(2,1)),
     &                       CCFBTCH(KCCF(1,2)),CCFBTCH(KCCF(2,2)),
     &                       R122(0,0),COOR1(1,0),COOR2(1,0),IPRINT)
            IF ( NPRF12 .GT. 0 ) THEN                         
             NBTCH = NBTCH + 1 
             IODBTCH(1,IBTCH) = IAOB1
             IODBTCH(2,IBTCH) = IAOB2
             IBTCH = IBTCH + 1
             NPRFMAX = MAX( NPRFMAX, NPRF12 ) 
            END IF
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
          IBTCH0 = IBTCH
        END IF

        END DO

       END DO
 
      END DO

      NODBCLS = ICLS 
      NODBTCH = IBTCH - 1

      RETURN
      END

C
C -----------------------------------------------------
C

      SUBROUTINE GETMAXODBCLS(MAXODBCLS,MAXODBTCH,IPRINT)
#include "implicit.h"
#include "priunit.h" 
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#include "pi.h"

      EXPMUMIN = EXPAOMIN/2 
      COEFMAX = R2PI52 / EXPAOMIN * CCFAOMAX * CCFAOMAX 
      R2AOCNTMAX  =  ABS( LOG( THRCOD/COEFMAX ) )/EXPMUMIN

      MAXODBTCH = 0 
      MAXODBCLS = 0 

      DO ICLS1 = 1, NAOBCLS
       ICENT1   = IAOCLS(1,ICLS1)
       NAOB1    = IAOCLS(3,ICLS1)
       DO ICLS2 = 1, ICLS1
        ICENT2   = IAOCLS(1,ICLS2)
        NAOB2    = IAOCLS(3,ICLS2)
        CALL GETR12( NSYMOPR, NSYMOPS, ICENT1, ICENT2, 2 )

        IF ( R122(0,0) .LE. R2AOCNTMAX ) THEN
          MAXODBCLS = MAXODBCLS + ( 2*MAXAONHK - 1 )
          MAXODBTCH = MAXODBTCH + NAOB1*NAOB2 
        END IF

       END DO
      END DO
      IF(iprint.GT.0)
     &WRITE(LUPRI,'(A,3F12.6,/,A,2I12)') 
     &        'RAOCNTMAX, EXPAOMIN, COEFMAX : ',
     &          SQRT(R2AOCNTMAX), EXPAOMIN, COEFMAX,
     &        'MAXODBCLS, MAXODBTCH           ',
     &         MAXODBCLS, MAXODBTCH   

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

      SUBROUTINE SCREENPRIM12(NPRF12,RODBTCH,IPR12,SETIPR,
     &                        NPRF1,NCTF1,NTYP1,NPRF2,NCTF2,NTYP2,
     &                        EXPBT1,EXPBT2,CCF11,CCF21,CCF12,CCF22,
     &                        R122,COOR1,COOR2,IPRINT)
#include"implicit.h"
#include"priunit.h"
#include"fck3thr.h"
#include"pi.h"
      DIMENSION IPR12(2,NPRF1*NPRF2),RODBTCH(5),
     &          COOR1(3),COOR2(3),EXPBT1(NPRF1),EXPBT2(NPRF2),
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
C    This factor for screening is a very  simple
C    approach based on exponencial prefactor.
C    It is not rigorous, and the best would be
C    to put GAB for real Cauchy-Schwarz screening
C

        FAC = R2PI52*EXP(-EXPMU*R122)*EXP12I * CCFSUM
C
C    What is better way of setting the Cauchy-Schwarz -
C       making sum of FAC, or just maximum ???
C
C        FACSCREEN = FACSCREEN + FAC

        IF ( FAC .GT. THRPOD ) THEN

         FACSCREEN = MAX(FACSCREEN,FAC)

         WT = ( CCFSUM11+CCFSUM21 ) * ( CCFSUM12 + CCFSUM22 )

         IPRF12 = IPRF12 + 1

         IF ( SETIPR ) THEN
          IPR12(1,IPRF12) = IPRF1
          IPR12(2,IPRF12) = IPRF2
         END IF

         WT = WT * EXP(-EXPMU*R122) 
         SM = SM + WT
         AM = AM + WT * EXP1 * EXP12I
         BM = BM + WT * EXP2 * EXP12I

         EXTP = SQRT(EXP12I) * FACCLS
         EXTPMAX = MAX( EXTPMAX,EXTP )

        END IF

       END DO
      END DO

      NPRF12 = IPRF12

      IF ( IPRF12 .GT. 0 ) THEN
       AM = AM/SM
       BM = BM/SM

       RODBTCH(1) = AM * COOR1(1) + BM * COOR2(1)
       RODBTCH(2) = AM * COOR1(2) + BM * COOR2(2)
       RODBTCH(3) = AM * COOR1(3) + BM * COOR2(3)

       RODBTCH(4) = EXTPMAX
       RODBTCH(5) = FACSCREEN

      END IF

      RETURN
      END 

C
C -----------------------------------------------------
C

      SUBROUTINE GETR12( NSYMOP1, NSYMOP2, 
     &                   ICENT1, ICENT2, NACTCENT )
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "nuclei.h"

      COOR1(1,0) = CORD(1,ICENT1)
      COOR1(2,0) = CORD(2,ICENT1)
      COOR1(3,0) = CORD(3,ICENT1)

      COOR2(1,0) = CORD(1,ICENT2)
      COOR2(2,0) = CORD(2,ICENT2)
      COOR2(3,0) = CORD(3,ICENT2)

      R12(1,0,0) = CORD(1,ICENT1) - CORD(1,ICENT2)
      R12(2,0,0) = CORD(2,ICENT1) - CORD(2,ICENT2)
      R12(3,0,0) = CORD(3,ICENT1) - CORD(3,ICENT2)

      R122(0,0)  = R12(1,0,0)**2+R12(2,0,0)**2+R12(3,0,0)**2 

      RETURN
      END 

C
C -----------------------------------------------------
C

      SUBROUTINE PUT2ERICOM(MAXODNHK,MAXODKCK,NCSQ1VAL,NCSQ2VAL)
#include"implicit.h"
#include"ericom.h"

      NHKMAX = MAXODNHK
      KCKMAX = MAXODKCK
      NCSQ1  = NCSQ1VAL
      NCSQ2  = NCSQ2VAL

      RETURN
      END
