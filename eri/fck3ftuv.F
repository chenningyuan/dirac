
C --------------------------------------------------------------
C
C             Ftuv
C           ========
C
C  Ftuv = Etij * Eukl * Evmn * Dcd
C
C  - summed for all the angular components of the primitive 
C integrals (ikm), (jln) with the corresponding density matrix 
C elements.
C
C  - Since we are constructing it for the second electron, 
C  there is a (-1)^(t+u+v) factor included for it in Etij.
C  ( In EXPCFT ). However for the nuclear attraction check
C  these signs are excluded
C
C  - for each total angular momentum in FTUV there is 
C   maximum NATOMS^2 classes: NFTUVCL <= NATOMS^2*MAXODNHK
C   based on one-center/two-center and center of the C
C   criteria
C 
C  - Class information is in IFTUVCL( 10, NFTUVCL )
C      1 - NTFTUVBTCH = number of FTUV's batches in the class
C      2 - NTFTUVPRIM = number of primitive FTUV's in the clas
C      3 - KSIFTUV   = pointer to the integer information about
C                    the first FTUV of this class in IFTUV
C                     ( IFTUV(KSIFTUV) )
C      4 - KSRFTUV   = pointer to the real information about
C                    the first FTUV of this class in RFTUV
C                     ( RFTUV(KSRFTUV) )
C      5 - KSFTUV   = pointer to the FTUV of this 
C                     class in FTUV ( FTUV(KSFTUV) )
C      6 - JMAXF    = total angular momentum in this class
C      7 - KSBTCH   = Number of the first OD batch in the class
C      8 - NTYPF    =  type of s-p combined overlap distribition
C                      - 1 if both contributing AO's are single AO's
C                      - 2 if at least on of the contributing AO's
C                        is a combined one from s and p functions 
C                        with the same exponents
C      9 - ICENTC   = center of the C basis functions 
C     10 - ICENTD   = center of the D basis functions 
C                     
C  - Real class infromation RFTUVCL( 4, NFTUVCL )
C      1 - distance between centers CD
C      2 - max. extent of q
C      3 - max. Ftuv
C      4 - max Cauchy screen parama
C
C  - Integer information about each FTUV:  IFTUV( 2, NTFTUV )
C
C      1 - the number of the originating OD 
C          ( for screening out the whole contracted OD 
C            together at the nonclassical test )
C      2 - number of primitives from CD in the same OD batch
C
C  - Real information about each FTUV in class is 
C   in RFTUV(7,NTFTUV)
C      1 - exponent q
C      2 - inverse of q ( 1/q )
C      3 - x-coordinate of Q
C      4 - y-coordinate of Q
C      5 - z-coordinate of Q
C      6 - max ftuv
C      7 - max Cauchy screen. param
C
C
C  - For each class the FTUV(KFTUV) is the beginning 
C    of array FTUV(NTFTUV,NDMAT,0:JMAXF,0:JMAXF,0:JMAXF)
C
C  - Real information about contracted OD corresponding 
C        to FTUV:  RFBTCH(5,NFTUVBTCH)
C  
C      1 - x-coordinate of the center
C      2 - y-coordinate of the center
C      3 - z-coordinate of the center
C      4 - non-classical extent (from the center)
C      5 - Cauchy screening parameter
C  - Integer information about contracted OD corresponding 
C        to FTUV:  IFBTCH(2,NFTUVBTCH)
C
C      1 - number of primitives for the OD
C      2 - relative index of the first primitive of the OD
C          in the class
C
C -------------------------------------------------------------

      SUBROUTINE SETFTUV( KIFTUVCL,KRFTUVCL,KIFTUV,KRFTUV,KFTUV,
     &                    KRFBTCH,KIFBTCH,IODBCLS,IODBTCH,DMAT,
     &                    ISYMDM,NDMAT,KFREE,LFREE,IELECTR,
     &                    WORK,LWORK,NPARFTUV_KSTARTODCLS,
     &                    NPARFTUV_NSTEPODCLS,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3mem.h"
C
      DIMENSION DMAT(NBASE,NBASE,NDMAT),ISYMDM(NDMAT),WORK(LWORK),
     &          IODBCLS(5,*),IODBTCH(2,*)

      CALL GETMAXSIZEFTUV(MAXSUMPRF,MAXSIZEFTUV,IODBCLS,IODBTCH,
     &                    NPARFTUV_KSTARTODCLS,NPARFTUV_NSTEPODCLS)

      MAXFTUVCLS  = NODBCLS
      MAXFTUVBTCH = NODBTCH

      NFTUVCLS    = MAXFTUVCLS

      KIFBTCH  = KFREE
      KRFBTCH  = KIFBTCH  +  2*MAXFTUVBTCH/IRAT + 1 
      KIFTUVCL = KRFBTCH  +  5*MAXFTUVBTCH 
      KRFTUVCL = KIFTUVCL + 10*MAXFTUVCLS/IRAT + 1
      KIFTUV   = KRFTUVCL +  4*MAXFTUVCLS 
      KRFTUV   = KIFTUV   +  2*MAXSUMPRF/IRAT + 1
      KFTUV    = KRFTUV   +  7*MAXSUMPRF
      KLAST1   = KFTUV    +    MAXSIZEFTUV

      MEMFTUV1 = KLAST1 - KFREE

C
C     Local allocations 
C     =================
 
      JTOP  = 4*(MAXAONHK - 1)
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6

      KINDHR  = KLAST1
      KINDSQ  = KINDHR + ( JTOP3 - 1 ) / IRAT + 1
      KIODHR  = KINDSQ + ( NRTOP - 1 ) / IRAT + 1
      KLMNPWR = KIODHR + ( 8 * NRTOP - 1) / IRAT + 1
      KLMNSYM = KLMNPWR + ( 3*MAXODNHK*MAXODKCK - 1)/IRAT + 1
      KCSQ    = KLMNSYM + (16*MAXODNHK*MAXODKCK - 1)/IRAT + 1
      KLAST   = KCSQ    + NCSQ1*NCSQ2

      MEMFTUV2 = KLAST - KLAST1

      LNEED = KLAST - KFREE

      IF ( LNEED .GT. LFREE ) THEN
         WRITE (LUPRI,'(/,2(2X,A,I15,/))')
     &      'Memory requirements total     : ',LNEED,
     &      'Memory available              : ',LFREE
         CALL QUIT('Insufficient memory in SETFTUV')
      END IF

      CALL DZERO(WORK(KFREE),LNEED)

      KWRK = KLAST
      LWRK = LFREE - LNEED 

      IF ( IPRINT .GT. 2 ) THEN 
         CALL HEADER('Memory allocations in SETFTUV:',-1)
         WRITE(LUPRI,'(9(3X,A,2I7,/))')
     &        'KIFBTCH  , LIFBTCH  :', KIFBTCH, KRFBTCH  - KIFBTCH,
     &        'KRFBTCH  , LRFBTCH  :', KRFBTCH, KIFTUVCL - KRFBTCH,
     &        'KIFTUVCL , LIFTUVCL :', KIFTUVCL,KRFTUVCL - KIFTUVCL,
     &        'KRFTUVCL , LRFTUVCL :', KRFTUVCL,KIFTUV   - KRFTUVCL,
     &        'KIFTUV   , LIFTUV   :', KIFTUV,  KRFTUV   - KIFTUV,
     &        'KRFTUV   , LRFTUV   :', KRFTUV,  KFTUV    - KRFTUV,
     &        'KFTUV    , LFTUV    :', KFTUV,   KLMNPWR  - KFTUV ,
     &        'KLMNPWR  , LLMNPWR  :', KLMNPWR, KLMNSYM  - KLMNPWR, 
     &        'KLMNSYM  , LLMNSYM  :', KLMNSYM, KCSQ     - KLMNSYM, 
     &        'KCSQ     , LCSQ     :', KCSQ   , KLAST    - KCSQ 
         WRITE(LUPRI,'(/,A,2I4)') 'NCSQ1, NCSQ2:',NCSQ1,NCSQ2
      END IF 


C
C     Spherical-harmonic transformation matrices
C     ==========================================
C

      CALL SPHMAT(WORK(KCSQ),NCSQ1,NCSQ2,0,.FALSE.,
     &            WORK(KWRK),LWRK,IPRINT)
C
C     Prepare index array INDHER and LMNPWR 
C     =====================================
C  ( for FTUV packing and to get carteesian indexes ijk 
C   corresponding to a certain spherical component )
C

      CALL HERPRP(WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR),
     &            IPRINT)
      CALL LMNPRP(WORK(KLMNPWR),WORK(KLMNSYM),IPRINT)

      CALL SETFTUV1(WORK(KIFTUVCL),WORK(KRFTUVCL),WORK(KIFTUV),
     &              WORK(KRFTUV),WORK(KFTUV),WORK(KRFBTCH),
     &              WORK(KIFBTCH),IODBCLS,IODBTCH,DMAT,ISYMDM,NDMAT,
     &              WORK(KINDHR),WORK(KLMNPWR),WORK(KLMNSYM),WORK(KCSQ),
     &              NPARFTUV_KSTARTODCLS,NPARFTUV_NSTEPODCLS,
     &              IELECTR,WORK(KWRK),LWRK,IPRINT)      

      LFREE = LFREE - ( KLAST1-KFREE )      
      KFREE = KLAST1

      RETURN
      END

      SUBROUTINE SETFTUV1(IFTUVCL,RFTUVCL,IFTUV,RFTUV,FTUV,
     &                    RFBTCH,IFBTCH,IODBCLS,IODBTCH,
     &                    DMAT,ISYMDM,NDMAT,INDHER,LMNPWR,LMNSYM,CSQ, 
     &                    NPARFTUV_KSTARTODCLS,NPARFTUV_NSTEPODCLS,
     &                    IELECTR,WORK,LWORK,IPRINT                    )
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#include "fck3mem.h"
#include "cbifck3.h"
      DIMENSION IFTUVCL(10,*),RFTUVCL(4,*),IFTUV(*),RFTUV(*),FTUV(*),
     &          RFBTCH(5,*),IFBTCH(2,*),IODBCLS(5,*),IODBTCH(2,*),
     &          DMAT(*),CSQ(NCSQ1,NCSQ2),WORK(*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          LMNPWR(MAXODKCK,MAXODNHK,3),
     &          LMNSYM(MAXODKCK,MAXODNHK,2,0:7),
     &          KCCF(2,2)

      KETIJ    = 1
      KPREF12  = KETIJ    + MAXODPRF * (MAXODNHK**3)*3*2
      KIPR12   = KPREF12  + MAXODPRF
      KCOOR12  = KIPR12   + MAXODPRF * 2 / IRAT + 1 
      KEXP12   = KCOOR12  + MAXODPRF * 3
      KHPI     = KEXP12   + MAXODPRF * 3
      KP1      = KHPI     + MAXODPRF
      KP2      = KP1      + MAXODPRF * 3
      KDMATSP  = KP2      + MAXODPRF * 3
      KDMATCP  = KDMATSP  + MAXODPRF * (MAXODKHK**2)
      KDTMP    = KDMATCP  + MAXODPRF * (MAXODKCK**2)
      KRNAI    = KDTMP    + MAXODKCK *  MAXODKHK
      KAHGTF   = KRNAI    + MAXODPRF * (MAXODNHK**3)
      KFTUVTMP = KAHGTF   +        2 * (MAXODNHK**3) 
      KRFTUVTMP= KFTUVTMP + MAXFCLSPRF *(MAXODNHK**3)
      KLAST    = KRFTUVTMP+ MAXFCLSPRF * 7

      LNEED    = KLAST    - 1 
      MEMFTUV2 = MEMFTUV2 + LNEED

      IF ( LNEED .GT. LWORK ) THEN         
         WRITE (LUPRI,'(/,4(2X,A,I15,/))')
     &      'MAXODPRF                      : ',MAXODPRF,
     &      'MAXFCLSPRF                    : ',MAXFCLSPRF,
     &      'MAXFTUVTMPSIZE                : ',MAXFCLSPRF*(MAXODNHK**3),
     &      'Memory requirements total     : ',LNEED,
     &      'Memory available              : ',LWORK
         CALL QUIT('Insufficient memory in SETFTUV1')
      END IF

      IFCLS  = 1
      KSRFTUV = 1
      KSIFTUV = 1
      KSFTUV  = 1

      NFTUVBTCH = 0

      DO ICLS = NPARFTUV_KSTARTODCLS, NODBCLS, NPARFTUV_NSTEPODCLS 

       EXTQMAX = 0.0D0
       FTUVMAX = 0.0D0
       SCRCD   = 0.0D0
       NTYP12  = 1

       NTFTUVBTCH = 0
       NTFTUVPRIM = 0

       ICENT1  = IODBCLS( 1, ICLS ) 
       ICENT2  = IODBCLS( 2, ICLS ) 
       JMAX12  = IODBCLS( 3, ICLS ) - 1 
       NHKT12  = IODBCLS( 3, ICLS ) 
       KSTBTCH = IODBCLS( 4, ICLS ) 
       NBTCH   = IODBCLS( 5, ICLS ) 

       IFTUVCL( 3,ICLS) = KSIFTUV
       IFTUVCL( 4,ICLS) = KSRFTUV
       IFTUVCL( 5,ICLS) = KSFTUV
       IFTUVCL( 6,ICLS) = JMAX12
       IFTUVCL( 7,ICLS) = KSTBTCH
       IFTUVCL( 9,ICLS) = ICENT1
       IFTUVCL(10,ICLS) = ICENT2

       KSRFTUV0 = KSRFTUV
       KSIFTUV0 = KSIFTUV
       KSFTUV0  = KSFTUV

       IF ( PACKFTUV ) THEN
         NTUV = NHKT12*(NHKT12+1)*(NHKT12+2)/6
       ELSE 
         NTUV = NHKT12 * NHKT12 * NHKT12
       END IF

       CALL GETR12( NSYMOPT, NSYMOPU, ICENT1, ICENT2, 2 )

       DO IBTCH = KSTBTCH, ( KSTBTCH+NBTCH-1)

        IAOB1 = IODBTCH( 1, IBTCH ) 
        IAOB2 = IODBTCH( 2, IBTCH ) 

        CALL GETAOB(NPRF1,NCTF1,NTYP1,NHKT1,JMAX1,KHKT1,KCKT1,
     &              KCSQ1,KEXP1,KCCF(1,1),IAOB1)

        CALL GETAOB(NPRF2,NCTF2,NTYP2,NHKT2,JMAX2,KHKT2,KCKT2,
     &              KCSQ2,KEXP2,KCCF(1,2),IAOB2)

        CALL SCREENPRIM12(NPRF12,RFBTCH(1,IBTCH),WORK(KIPR12),.TRUE.,
     &                    NPRF1,NCTF1,NTYP1,NPRF2,NCTF2,NTYP2,
     &                    EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &                    CCFBTCH(KCCF(1,1)),CCFBTCH(KCCF(2,1)),
     &                    CCFBTCH(KCCF(1,2)),CCFBTCH(KCCF(2,2)),
     &                    R122(0,0),COOR1(1,0),COOR2(1,0),IPRINT)

        EXTQMAX = MAX( RFBTCH(4,IBTCH), EXTQMAX )
        SCRCD   = MAX( RFBTCH(5,IBTCH), SCRCD   )

        IF ( NTYP1+NTYP2 .GT. 2 ) NTYP12 = 2

        CALL SETETIJ(WORK(KETIJ),WORK(KPREF12),
     &               R12(1,0,0),R122(0,0),
     &               WORK(KIPR12),EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &               WORK(KCOOR12),WORK(KEXP12),WORK(KHPI),WORK(KP1),
     &               WORK(KP2),NPRF1,NHKT1,JMAX1,
     &               NPRF2,NHKT2,JM AX2,NPRF12,NHKT12,IELECTR,IPRINT)

        J1 = JMAX1 
        J2 = JMAX2 
        JMAX12 = NHKT12 - 1 

        CALL DZERO(WORK(KFTUVTMP+KSFTUV-KSFTUV0),NPRF12*NTUV)

        DO ICOMP1 = 1, NTYP1
         IF ( NTYP1 .EQ. 2 ) CALL SETJKKK(J1,KCSQ1,KCKT1,KHKT1,ICOMP1)

         DO ICOMP2 = 1, NTYP2 
          IF ( NTYP2 .EQ. 2 ) CALL SETJKKK(J2,KCSQ2,KCKT2,KHKT2,ICOMP2)

          CALL SETDMATCP(WORK(KDMATCP),DMAT,IAOB1,IAOB2,ICOMP1,ICOMP2,
     &        CSQ(KCSQ1,1),CSQ(KCSQ2,1),WORK(KIPR12),
     &        CCFBTCH(KCCF(ICOMP1,1)),CCFBTCH(KCCF(ICOMP2,2)),
     &        WORK(KDMATSP), WORK(KDTMP),KCKT1,KCKT2,KHKT1,KHKT2,
     &        NPRF1,NPRF2,NCTF1,NCTF2,NPRF12,IPRINT)

         CALL CONTRACTDMATE(WORK(KFTUVTMP+KSFTUV-KSFTUV0),
     &        WORK(KDMATCP),WORK(KETIJ),WORK(KPREF12),INDHER,LMNPWR,
     &        KCKT1,KCKT2,J1,J2,JMAX1,JMAX2,JMAX12,NPRF12,NTUV,IPRINT)

         END DO
        END DO

        CALL SETRFTUV(IFTUV(KSIFTUV),RFTUV(KSRFTUV),
     &         COOR1(1,0),COOR2(1,0),EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &         WORK(KIPR12),NPRF12,IBTCH,IPRINT)

        IFBTCH(1,IBTCH) = NPRF12
        IFBTCH(2,IBTCH) = NTFTUVPRIM+1

        KSIFTUV = KSIFTUV + NPRF12 * 2
        KSRFTUV = KSRFTUV + NPRF12 * 7 
        KSFTUV  = KSFTUV  + NPRF12 * NTUV

        NTFTUVBTCH = NTFTUVBTCH + 1
        NTFTUVPRIM = NTFTUVPRIM + NPRF12

       END DO

C       CALL SORTFTUV(FTUV(KSFTUV0),RFTUV(KSRFTUV0),
C     &               WORK(KFTUVTMP),WORK(KRFTUVTMP),
C     &               RFBTCH(1,KSTBTCH),IFBTCH(1,KSTBTCH),
C     &               NTUV,NTFTUVPRIM,NTFTUVBTCH,IPRINT)

       CALL RTRMAT(NTUV,NTFTUVPRIM,WORK(KFTUVTMP),FTUV(KSFTUV0))

       IF ( IPRINT .GE. 5 ) THEN 
        CALL PRIFTUV2(FTUV(KSFTUV0),ICLS,NTFTUVPRIM,JMAX12)

        DO I=1,NTFTUVPRIM
         WRITE(LUPRI,'(3X,I5,7F12.7,I4)')
     &         I,(RFTUV(KSRFTUV0+7*(I-1)+J-1),J=1,7),
     &        IFTUV(KSIFTUV0+2*(I-1)+1)
        END DO      
       END IF
  
       NFTUVBTCH = NFTUVBTCH + NTFTUVBTCH

       IFTUVCL( 1, ICLS ) = NTFTUVBTCH
       IFTUVCL( 2, ICLS ) = NTFTUVPRIM
       IFTUVCL( 8, ICLS)  = NTYP12

       RFTUVCL( 1, ICLS ) = SQRT( R122(0,0) ) 
       RFTUVCL( 2, ICLS ) = EXTQMAX
       RFTUVCL( 3, ICLS ) = FTUVMAX
       RFTUVCL( 4, ICLS ) = SCRCD

      END DO

      RETURN 
      END 

      SUBROUTINE GETMAXSIZEFTUV(MAXSUMPRF,MAXSIZEFTUV,IODBCLS,IODBTCH,
     &                    NPARFTUV_KSTARTODCLS,NPARFTUV_NSTEPODCLS)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "cbifck3.h"
      DIMENSION IODBCLS(5,*),IODBTCH(2,*) 

      MAXFCLSPRF  = 0
      MAXSUMPRF   = 0
      MAXSIZEFTUV = 0

      DO ICLS = NPARFTUV_KSTARTODCLS, NODBCLS, NPARFTUV_NSTEPODCLS

       NCLSPRF = 0

       KSTBTCH = IODBCLS(4,ICLS)
       NBTCH   = IODBCLS(5,ICLS)
       
       DO IBTCH = KSTBTCH, (KSTBTCH + NBTCH -1 )
        IAOB1 = IODBTCH( 1, IBTCH ) 
        IAOB2 = IODBTCH( 2, IBTCH ) 

        NHKT1 = NHKTBTCH( IAOB1 )
        NHKT2 = NHKTBTCH( IAOB2 )

        NHKT = NHKT1 + NHKT2 - 1 
 
        NPRF1 = NPRFBTCH( IAOB1 ) 
        NPRF2 = NPRFBTCH( IAOB2 ) 

        NPRF = NPRF1 * NPRF2

        IF ( PACKFTUV ) THEN
         NTUV = NHKT * ( NHKT + 1 ) * ( NHKT + 2 ) / 6
        ELSE 
         NTUV = NHKT * NHKT * NHKT
        END IF

        NCLSPRF     = NCLSPRF     + NPRF
        MAXSUMPRF   = MAXSUMPRF   + NPRF
        MAXSIZEFTUV = MAXSIZEFTUV + NPRF * NTUV 
       END DO

       MAXFCLSPRF = MAX( MAXFCLSPRF, NCLSPRF ) 

      END DO 

      RETURN 
      END 

      SUBROUTINE SETETIJ(ETIJ,PREF12,R12,R122,IPR12,EXPBT1,EXPBT2,
     &                   COOR12,EXP12,HPI,P1,P2,
     &                   NPRF1,NHKT1,JMAX1,NPRF2,NHKT2,JMAX2,
     &                   NPRF12,NHKT12,IELECTR,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
      DIMENSION ETIJ(*),PREF12(NPRF12),
     &          COOR12(NPRF12,3),EXP12(NPRF12,3), 
     &          HPI(*), P1(NPRF12,3), P2(NPRF12,3),R12(3),I120(3),
     &          EXPBT1(NPRF1),EXPBT2(NPRF2),IPR12(2,NPRF12)

      INTEGER X

      CALL DZERO(ETIJ,NPRF12*NHKT12*NHKT1*NHKT2*3*2)

      CALL LOCSYM(R12,1,1,I120,IXYZ)

      DO X=1,3
       DO I=1, NPRF12
        COOR12(I,X) = R12(X)
       END DO       
      END DO

      DO I=1,NPRF12
      
       IPRF1 = IPR12(1,I)
       IPRF2 = IPR12(2,I)

       EXP1  = EXPBT1(IPRF1)
       EXP2  = EXPBT2(IPRF2)
 
       EXPPI = 1.0D0 / ( EXP1 + EXP2 )

       EXP12(I,1) = EXP1 * EXPPI
       EXP12(I,2) = EXP2 * EXPPI
       EXP12(I,3) = EXPPI

       PREF12(I) = exp( -EXP1*EXP2*EXPPI*R122 )

       DO X = 1,3
        P1(I,X) = - EXP12(I,2) * COOR12(I,X) 
        P2(I,X) =   EXP12(I,1) * COOR12(I,X)
       END DO 

      END DO

C
C  IRUTIN - here just to print it
C
      IRUTIN = 1

      ISCAL1 = 0 
      ISCAL2 = 0 
C
C  level of derivatives calculated
C
      IDL    = 0 
      NCNT12 = 0

      CALL EXPCFT(ETIJ,NPRF12,JMAX1,JMAX2,COOR12,EXP12,I120,HPI,
     &            P1,P2,IRUTIN,IELECTR,JMAX1,JMAX2,NCNT12,IDL,
     &            ISCAL1,ISCAL2,IPRINT,DUM,DUM,DUM,DUM)

      IF ( IPRINT .GT. 10 ) THEN
       CALL HEADER('Exponents in SETETIJ',-1)
       DO I=1,NPRF12
        EXP1 = EXP12(I,1) / EXP12(I,3)
        EXP2 = EXP12(I,2) / EXP12(I,3)
        WRITE(LUPRI,'(A,I4,A,I4,3X,2F12.6)')
     &     'IPRF1:',IPR12(1,I),',   IPRF2:',IPR12(2,I),
     &     EXP1,EXP2
       END DO
      END IF 

      RETURN 
      END 

      SUBROUTINE SETDMATCP(DMATCP,DMAT,IAOB1,IAOB2,ISPCOMP1,ISPCOMP2,
     &              CSQ1,CSQ2,IPR12,CCF1,CCF2,DMATSP,DTMP,KCKT1,KCKT2,
     &              KHKT1,KHKT2,NPRF1,NPRF2,NCTF1,NCTF2,NPRF12,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"

      DIMENSION DMAT(NBASE,NBASE),
     &          DMATCP(KCKT1,KCKT2,NPRF12),
     &          DMATSP(KHKT1,KHKT2,NPRF12),
     &          DTMP(KCKT1,KHKT2),
     &          CSQ1(KHKT1,KCKT1), CSQ2(KHKT2,KCKT2),
     &          CCF1(NPRF1,NCTF1), CCF2(NPRF2,NCTF2), 
     &          IPR12(2,NPRF12)

      CALL DZERO(DMATCP, KCKT1*KCKT2*NPRF12 )
      CALL DZERO(DMATSP, KHKT1*KHKT2*NPRF12 )
      CALL DZERO(DTMP, KCKT1*KHKT2 )

      IPRF12 = 1

      FAC = 1.0D0
      IF ( IAOB1 .EQ. IAOB2 ) FAC = 0.5D0

      DO IPRF12 = 1, NPRF12
       IPRF1 = IPR12(1,IPRF12)
       IPRF2 = IPR12(2,IPRF12)

        DO ICTF1 = 1, NCTF1
         ISTAO1 = INDXBTCH(KNDXBTCH(ISPCOMP1,IAOB1) + ICTF1 - 1 )

         DO ICTF2 = 1, NCTF2
          ISTAO2 = INDXBTCH(KNDXBTCH(ISPCOMP2,IAOB2) + ICTF2 - 1 )

          DO ISP1=1, KHKT1
           IAO1 = ISTAO1+ISP1-1
           DO ISP2=1, KHKT2
            IAO2 = ISTAO2 + ISP2 - 1
            DMATSP(ISP1,ISP2,IPRF12) =  DMATSP(ISP1,ISP2,IPRF12) 
     &                   + FAC * CCF1(IPRF1,ICTF1)*CCF2(IPRF2,ICTF2)
     &                     *  ( DMAT(IAO1,IAO2) + DMAT( IAO2, IAO1 ) )     
           END DO
          END DO

         END DO
        END DO

        CALL DGEMM('T','N', KCKT1, KHKT2, KHKT1, 1.0D0,
     &              CSQ1, KHKT1, DMATSP(1,1,IPRF12), KHKT1, 
     &              0.0D0, DTMP, KCKT1  )      

        CALL DGEMM('N','N', KCKT1, KCKT2, KHKT2, 1.0D0,
     &              DTMP, KCKT1, CSQ2, KHKT2, 
     &              0.0D0, DMATCP(1,1,IPRF12), KCKT1  )  

      END DO

      IF ( IPRINT .GE. 10 ) THEN   
        CALL HEADER('DMAT (AO)',-1)
        CALL OUTPUT(DMAT,1,NBASE,1,NBASE,NBASE,NBASE,1,LUPRI)

        CALL HEADER('Contraction Matrices',-1)
        WRITE (LUPRI,1000) 'NPRF1, NCTF1',NPRF1,NCTF1,
     &                     'NPRF2, NCTF2',NPRF2,NCTF2

        CALL HEADER('CCF1',5)
        CALL OUTPUT(CCF1,1,NPRF1,1,NCTF1,NPRF1,NCTF1,1,LUPRI)

        CALL HEADER('CCF2',5)
        CALL OUTPUT(CCF2,1,NPRF2,1,NCTF2,NPRF2,NCTF2,1,LUPRI)

        CALL HEADER('DMATSP - Spherical primitives',-1)
        WRITE (LUPRI,1100) 
     &       'KHKT1, KHKT2, KHKT12', KHKT1, KHKT2, KHKT1*KHKT2,
     &       'KPRF1, KPRF2, NPRF12', NPRF1, NPRF2, NPRF1*NPRF2

        CALL OUTPUT(DMATSP,1,KHKT1*KHKT2,1,NPRF1*NPRF2,
     &             KHKT1*KHKT2, NPRF1*NPRF2,1,LUPRI)
  
        CALL HEADER('Cartesian-spherical Matrices',-1)

        WRITE (LUPRI,1000) 'KHKT1, KCKT1', KHKT1, KCKT1,
     &                     'KHKT2, KCKT2', KHKT2, KCKT2

        CALL HEADER('CSQ1',5)
        CALL OUTPUT(CSQ1,1,KHKT1,1,KCKT1,KHKT1,KCKT1,1,LUPRI)

        CALL HEADER('CSQ2',5)
        CALL OUTPUT(CSQ2,1,KHKT2,1,KCKT2,KHKT2,KCKT2,1,LUPRI)

        CALL HEADER('DMATCP - Cartesian primitives',-1)
          
        WRITE (LUPRI,1100) 
     &       'KCKT1, KCKT2, KCKT12',KCKT1, KCKT2, KCKT1*KCKT2,
     &       'KPRF1, KPRF2, NPRF12',NPRF1, NPRF2, NPRF1*NPRF2

        CALL OUTPUT(DMATCP,1,KCKT1*KCKT2,1,NPRF1*NPRF2,
     &             KCKT1*KCKT2, NPRF1*NPRF2,1,LUPRI)

      END IF 

 1000 FORMAT (5X,'Dimensions : ',A,2I10,/,18X,A,2I10,/)
 1100 FORMAT (5X,'Dimensions : ',A,3I10,/,18X,A,3I10,/)

      RETURN
      END

      SUBROUTINE CONTRACTDMATE(FTUVTMP,DMATCP,ETIJ,PREF12,INDHER,LMNPWR,
     &                      KCKT1,KCKT2,J1,J2,JMAX1,JMAX2,JMAX12,NPRF12,
     &                      NTUV,IPRINT)
#include"implicit.h"
#include"priunit.h"
#include"fck3odb.h"
#include"cbifck3.h"
 
      DIMENSION FTUVTMP(NTUV,NPRF12),
     &          ETIJ(NPRF12,0:JMAX12,0:JMAX1,0:JMAX2,3),
     &          DMATCP(KCKT1,KCKT2,*),PREF12(NPRF12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),LMNPWR(MAXODKCK,MAXODNHK,3)

      INTEGER T,U,V,TUV

      NHKT12  = JMAX12 + 1
      NHKT122 = NHKT12 * NHKT12

      DO T = 0, JMAX12
       DO U = 0, JMAX12 - T 
        DO V = 0, JMAX12 - T - U  

         IF ( PACKFTUV ) THEN
          TUV = INDHER(T,U,V)
         ELSE
          TUV = 1 + T + NHKT12 * U + V * NHKT122
         END IF

         DO KC1 = 1, KCKT1

          L1 = LMNPWR(KC1,J1+1,1)
          M1 = LMNPWR(KC1,J1+1,2)
          N1 = LMNPWR(KC1,J1+1,3)

          DO KC2 = 1, KCKT2
 
           L2 = LMNPWR(KC2,J2+1,1)
           M2 = LMNPWR(KC2,J2+1,2)
           N2 = LMNPWR(KC2,J2+1,3)

           DO I=1,NPRF12
            FTUVTMP(TUV,I) = FTUVTMP(TUV,I)+
     &             DMATCP(KC1,KC2,I) * PREF12(I) * ETIJ(I,T,L1,L2,1)
     &           * ETIJ(I,U,M1,M2,2) * ETIJ(I,V,N1,N2,3) 
           END DO

          END DO

         END DO

        END DO
       END DO
      END DO

      IF ( IPRINT .GT. 10 ) THEN 
       CALL HEADER('Output from CONTRACTDMATE',-1)

       DO I=1,NPRF12
        DO T = 0, JMAX12
         DO U = 0, JMAX12 - T 
          DO V = 0, JMAX12 - T - U  
           IF ( PACKFTUV ) THEN
            TUV = INDHER(T,U,V)
           ELSE
            TUV = 1 + T + NHKT12 * U + V * NHKT122            
           END IF
           WRITE(LUPRI,'(3X,4I4,1F16.8)') I,T,U,V,FTUVTMP(TUV,I) 
          END DO
         END DO
        END DO
       END DO

      END IF

      RETURN
      END

      SUBROUTINE SETJKKK( J, KCSQ, KCKT, KHKT, ICOMP )
#include "implicit.h"
   
      J    = ICOMP - 1
      KCSQ = KSQADR(J,0,0,0,0)
      KCKT = ( J + 1 ) * ( J + 2 ) / 2
      KHKT = 2 * J + 1

      RETURN 
      END 

      SUBROUTINE SETRFTUV(IFTUV,RFTUV,COOR1,COOR2,EXPBT1,EXPBT2,
     &                    IPR12,NPRF12,IBTCH,IPRINT)
#include "implicit.h"
#include "priunit.h"

      DIMENSION IFTUV(2,NPRF12),RFTUV(7,NPRF12),COOR1(3),COOR2(3),
     &          EXPBT1(*),EXPBT2(*),IPR12(2,*)

      DO I=1,NPRF12

       IFTUV(1,I) = IBTCH
       IFTUV(2,I) = NPRF12
      
       IPRF1 = IPR12(1,I)
       IPRF2 = IPR12(2,I)

       EXP1  = EXPBT1(IPRF1)
       EXP2  = EXPBT2(IPRF2)

       EXP12  = EXP1 + EXP2         
       EXP12I = 1.0D0 / ( EXP1 + EXP2 )

       RFTUV(1,I) = EXP12
       RFTUV(2,I) = EXP12I
       RFTUV(3,I) = ( COOR1(1)*EXP1 + COOR2(1)*EXP2 ) * EXP12I
       RFTUV(4,I) = ( COOR1(2)*EXP1 + COOR2(2)*EXP2 ) * EXP12I
       RFTUV(5,I) = ( COOR1(3)*EXP1 + COOR2(3)*EXP2 ) * EXP12I

C
C      These two parameters are not used at the moment :
C
       RFTUV(6,I) = 0.0D0
       RFTUV(7,I) = 0.0D0

      END DO
      
      RETURN
      END

      SUBROUTINE SORTFTUV(FTUV0,RFTUV,FTUVTMP,RFTUVTMP,RFBTCH,
     &                    IFBTCH,NTUV,NTFTUVPRIM,NTFTUVBTCH,IPRINT)
#include"implicit.h"
#include"priunit.h"
      DIMENSION FTUV0(NTFTUVPRIM,NTUV),FTUVTMP(NTUV,NTFTUVPRIM),
     &          RFTUV(7,*),RFTUVTMP(7,*),RFBTCH(5,*),IFBTCH(2,*)

C
C  Here we sort the FTUV decreasingly according to 
C Cauchy screening factor for FTUV batches
C
C  Surely here one could think of a faster algorithm for 
C sorting ... ( althogh NTFTUVBTCH is only the number 
C of OD batches in the class, it wouldn't hurt to 
C improve it )
C

      IPRF0 = 1

      DO I=1,NTFTUVBTCH-1

       SCR = RFBTCH(5,I)
       ISCR = I

C
C   Find the nect element with the biggest screening factor
C

       DO J = I+1, NTFTUVBTCH 
        IF ( RFBTCH(5,J) .GT. SCR ) THEN
          SCR  = RFBTCH(5,J)
          ISCR = J
        END IF
       END DO
     
C
C   ... and interchange the OD information to the possition starting 
C     at I, and put the FTUV batch to the FTUV0 starting at IPRF

       IF ( ISCR .NE. I ) THEN 
        DO K=1,5
         RTMP           = RFBTCH(K,I)
         RFBTCH(K,I)    = RFBTCH(K,ISCR)   
         RFBTCH(K,ISCR) = RTMP
        END DO

        DO K=1,2
         ITMP           = IFBTCH(K,I)
         IFBTCH(K,I)    = IFBTCH(K,ISCR)   
         IFBTCH(K,ISCR) = ITMP
        END DO
       END IF

C
C  If ISCR was different from I we have
C  copied IFBTCH(K,ISCR) to IFBTCH(K,I) already
C  --> we read the data from IFBTCH(K,I) here.
C

       NPRF    = IFBTCH(1,I)
       IPRFTMP = IFBTCH(2,I)
       IFBTCH(2,I) = IPRF0

       DO IPRFOFF = 0, NPRF-1
        IPRF1 = IPRF0   + IPRFOFF
        IPRF2 = IPRFTMP + IPRFOFF

        DO K=1,7
         RFTUV(K,IPRF1) = RFTUVTMP(K,IPRF2)   
        END DO
        DO J=1,NTUV
         FTUV0(IPRF1,J) = FTUVTMP(J,IPRF2)
        END DO
       END DO

       IPRF0 = IPRF0 + NPRF
      END DO

      NPRF    = IFBTCH(1,NTFTUVBTCH)
      IPRFTMP = IFBTCH(2,NTFTUVBTCH)

      DO IPRFOFF = 0, NPRF-1
       IPRF1 = IPRF0   + IPRFOFF
       IPRF2 = IPRFTMP + IPRFOFF
       DO K=1,7
        RFTUV(K,IPRF1) = RFTUVTMP(K,IPRF2)   
       END DO
       DO J=1,NTUV
        FTUV0(IPRF1,J) = FTUVTMP(J,IPRF2)
       END DO
      END DO

      RETURN
      END 

C
C ---------------------------------------------------------------------
C

      SUBROUTINE PRIFTUV2(FTUV,ICL,NPRFAB,JMAX)
#include "implicit.h"
#include "priunit.h"       
      DIMENSION FTUV(NPRFAB,0:JMAX,0:JMAX,0:JMAX)
      INTEGER T,U,V

      CALL HEADER('FTUV',-1)

      WRITE(LUPRI,'(/,7X,3(A,I5,A,3X),/)')
     &            'FTUV class number :',ICL,',',
     &            'NTFTUV :', NPRFAB,',',
     &            'JMAXF  :', JMAX,' '

      WRITE(LUPRI,'(7X,A)') 
     &            'I    T    U    V     FTUV(I,T,U,V)'

      DO I = 1, NPRFAB

       DO T = 0, JMAX
        DO U = 0, JMAX
         DO V = 0, JMAX

          IF ( T+U+V .LE. JMAX ) THEN
            WRITE(LUPRI,'(3X,4I5,F16.8)') I,T,U,V,FTUV(I,T,U,V)
          ELSE
            IF ( FTUV(I,T,U,V) .GT. 0D0 ) THEN 
             WRITE(LUPRI,'(3X,4I5,F16.8,A)') I,T,U,V,FTUV(I,T,U,V),
     &                   ' ERROR - T+U+V > JMAX ! '
            END IF 
          END IF 
         END DO
        END DO
       END DO

      END DO

      RETURN
      END

