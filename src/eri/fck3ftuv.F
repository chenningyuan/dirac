!dalton_copyright_start
!      Copyright (c) 2011 by the authors of Dalton (see below).
!      All Rights Reserved.
!   
!      The source code in this file is part of
!      "Dalton, a molecular electronic structure program,
!       Release DALTON2011 (2011), see http://daltonprogram.org"
!   
!      This source code is provided under a written licence and may be
!      used, copied, transmitted, or stored only in accord with that
!      written licence.
!   
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the licence holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the licence holder may not take any part of Dalton,
!      including modified files, with him/her, unless that person has
!      obtained his/her own licence.
!   
!      For further information, including how to get a licence, see:
!         http://daltonprogram.org
!dalton_copyright_end

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
C      2 - NTFTUVPRIM = number of primitive FTUV's in the class
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
C      8 - NTYPF    =  type of s-p combined overlap distribution
C                      - 1 if both contributing AO's are single AO's
C                      - 2 if at least on of the contributing AO's
C                        is a combined one from s and p functions 
C                        with the same exponents
C      9 - ICENTC   = center of the C basis functions 
C     10 - ICENTD   = center of the D basis functions 
C                     
C  - Real class infromation RFTUVCL( 6, NFTUVCL )
C      1 - distance between centers CD
C      2 - max. extent of q
C      3 - max. Ftuv
C      4 - max Cauchy screen parama
C      5 - max. exponent of q
C      6 - max. contracted density and CS element, Dab*Gab
C
C  - Integer information about each FTUV:  IFTUV( 2, NTFTUV )
C
C      1 - the number of the originating OD 
C          ( for screening out the whole contracted OD 
C            together at the nonclassical test )
C      2 - number of primitives from CD in the same OD batch
C
C  - Real information about each FTUV in class is 
C   in RFTUV(8,NTFTUV)
C      1 - exponent q
C      2 - inverse of q ( 1/q )
C      3 - x-coordinate of Q
C      4 - y-coordinate of Q
C      5 - z-coordinate of Q
C      6 - max ftuv
C      7 - max Cauchy screen. param
C      8 - weight for Yang's partitioning of the density
C
C
C  - For each class the FTUV(KFTUV) is the beginning 
C    of array FTUV(NTFTUV,NDMAT,0:JMAXF,0:JMAXF,0:JMAXF)
C
C  - Real information about contracted OD corresponding 
C        to FTUV:  RFBTCH(6,NFTUVBTCH)
C  
C      1 - x-coordinate of the center
C      2 - y-coordinate of the center
C      3 - z-coordinate of the center
C      4 - non-classical extent (from the center)
C      5 - Cauchy screening parameter
C      6 - max. contracted density and CS element, Dab*Gab
C
C  - Integer information about contracted OD corresponding 
C        to FTUV:  IFBTCH(2,NFTUVBTCH)
C
C      1 - number of primitives for the OD
C      2 - relative index of the first primitive of the OD
C          in the class
C
C -------------------------------------------------------------

      SUBROUTINE SETFTUV( KIFTUVCL,KRFTUVCL,KIFTUV,KRFTUV,KFTUV,
     &                    KRFBTCH,KIFBTCH,RODBCLS,IODBCLS,RODBTCH,
     &                    IODBTCH,DMAT,
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
#include "densfit.h"
C
#include "pbc.h"
C
      DIMENSION DMAT(NDIMC,NDIMD),ISYMDM(NDMAT),WORK(LWORK),
     &          IODBCLS(5,*),IODBTCH(2,*),RODBCLS(2,*),RODBTCH(7,*)

      CALL GETMAXSIZEFTUV(MAXSUMPRF,MAXSIZEFTUV,IODBCLS,IODBTCH,
     &                    NPARFTUV_KSTARTODCLS,NPARFTUV_NSTEPODCLS)

      MAXFTUVCLS  = NODBCLS
      MAXFTUVBTCH = NODBTCH

      NFTUVCLS    = MAXFTUVCLS

      KIFBTCH  = KFREE
      KRFBTCH  = KIFBTCH  +  2*MAXFTUVBTCH/IRAT + 1 
      KIFTUVCL = KRFBTCH  +  6*MAXFTUVBTCH 
      KRFTUVCL = KIFTUVCL + 10*MAXFTUVCLS/IRAT + 1
      KIFTUV   = KRFTUVCL +  6*MAXFTUVCLS 
      KRFTUV   = KIFTUV   +  2*MAXSUMPRF/IRAT + 1
      KFTUV    = KRFTUV   +  8*MAXSUMPRF
      KLAST1   = KFTUV    +    MAXSIZEFTUV

      MEMFTUV1 = KLAST1 - KFREE

C
C     Local allocations 
C     =================
 
      JTOP  = 4*(MAXAONHK - 1)
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6

      CALL PUT2HERTOP(JTOP,NRTOP)
C
C     Dimensions for CSQ
C     ==================
C     
      CALL SPHDIM


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

      CALL HERPRP(WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR))
      CALL LMNPRP(WORK(KLMNPWR),WORK(KLMNSYM),IPRINT)

      CALL SETFTUV1(WORK(KIFTUVCL),WORK(KRFTUVCL),WORK(KIFTUV),
     &              WORK(KRFTUV),WORK(KFTUV),WORK(KRFBTCH),
     &              WORK(KIFBTCH),RODBCLS,IODBCLS,RODBTCH,IODBTCH,DMAT,
     &              ISYMDM,
     &              NDMAT,WORK(KINDHR),WORK(KLMNPWR),WORK(KLMNSYM),
     &              WORK(KCSQ),NPARFTUV_KSTARTODCLS,NPARFTUV_NSTEPODCLS,
     &              IELECTR,WORK(KWRK),LWRK,IPRINT)      


      IF ( IELECTR .EQ. 1) THEN 
        if (dopbc .and. pbcactive) then
           write(LUPRI,*) 'SETFTUV: No call to CHECKNUCATR with PBC.'
           call quit('SETFTUV: No call to CHECKNUCATR with PBC.')
        end if
        CALL CHECKNUCATR(WORK(KIFTUVCL),WORK(KRFTUV),WORK(KFTUV),
     &               DMAT(1,1),WORK(KINDHR),JTOP,NDMAT,NFTUVCLS,
     &               WORK(KWRK),LWRK,IPRINT)
      END IF 

      LFREE = LFREE - ( KLAST1-KFREE )      
      KFREE = KLAST1

      RETURN
      END

      SUBROUTINE SETFTUV1(IFTUVCL,RFTUVCL,IFTUV,RFTUV,FTUV,
     &                    RFBTCH,IFBTCH,RODBCLS,IODBCLS,RODBTCH,IODBTCH,
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
#include "densfit.h"
C
#include "pbc.h"
      DIMENSION IFTUVCL(10,*),RFTUVCL(6,*),IFTUV(*),RFTUV(*),FTUV(*),
     &          RFBTCH(6,*),IFBTCH(2,*),IODBCLS(5,*),IODBTCH(2,*),
     &          DMAT(*),CSQ(NCSQ1,NCSQ2),WORK(*), RODBCLS(2,*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          LMNPWR(MAXODKCK,MAXODNHK,3),
     &          LMNSYM(MAXODKCK,MAXODNHK,2,0:7),
     &          KCCF(2,2),RODBTCH(7,*)
      dimension tvec1(3), tvec2(3)
      LOGICAL ACTIVE
      PARAMETER (ONE=1.0d0, HALF=0.5d0)

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
      KLAST    = KRFTUVTMP+ MAXFCLSPRF * 8

      LNEED    = KLAST    - 1 
      MEMFTUV2 = MEMFTUV2 + LNEED
C
C aspg, 2006-06-27 making ftnchek happy
C
      IF ( LNEED .GT. LWORK ) THEN         
         WRITE (LUPRI,'(/,4(2X,A,I15,/))')
     &      'MAXODPRF                      : ',MAXODPRF,
     &      'MAXCLSPRF                     : ',MAXFCLSPRF,
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

      NCLSLAST = NODCLSNORM
      IF (DOFIT.AND.RHSFIT) NCLSLAST = NODBCLS

      DO ICLS = NPARFTUV_KSTARTODCLS, NCLSLAST, NPARFTUV_NSTEPODCLS 

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

       if (dopbc .and. pbcactive) then
C         We ultimately want J_ab = \sum_mj \sum_cd (ab|cd) D_cd
C         Here we add m to the c-center and m+j to the d-center. /ErikT
          if (pbc_sep_Rtensor) then
             tvec1(1) = nf_vec_std(1)
             tvec1(2) = nf_vec_std(2)
             tvec1(3) = nf_vec_std(3)
             tvec2(1) = nf_vec_std(1) + nf2_vec_std(1)
             tvec2(2) = nf_vec_std(2) + nf2_vec_std(2)
             tvec2(3) = nf_vec_std(3) + nf2_vec_std(3)
          else
C            The false value of pbc_set_Rtensor indicates that we
C            should ignore the m-translation (nf_vec_std) here, since
C            this is instead taken into account later in the summation
C            of R tensors.
             tvec1(1) = 0.0
             tvec1(2) = 0.0
             tvec1(3) = 0.0
             tvec2(1) = nf2_vec_std(1)
             tvec2(2) = nf2_vec_std(2)
             tvec2(3) = nf2_vec_std(3)
          end if
          CALL GETR12( NSYMOPT, NSYMOPU, ICENT1, ICENT2, 2,
     &                 tvec1, tvec2 )
       else
          tvec1(1) = 0.0
          tvec1(2) = 0.0
          tvec1(3) = 0.0
          tvec2(1) = 0.0
          tvec2(2) = 0.0
          tvec2(3) = 0.0
          CALL GETR12( NSYMOPT, NSYMOPU, ICENT1, ICENT2, 2,
     &                 tvec1, tvec2 )
       end if

       DO IBTCH = KSTBTCH, ( KSTBTCH+NBTCH-1)

        IAOB1 = IODBTCH( 1, IBTCH ) 
        IAOB2 = IODBTCH( 2, IBTCH ) 

c       IF (.TRUE.) THEN
        ACTIVE = ACTVAOBTCH(IAOB1).OR.ACTVAOBTCH(IAOB2)
        IF ((.NOT.DOBOX).OR.ACTIVE) THEN
C       Weight added for Yang's partitioning of the density according to
C       Gallant, R.T., Chem. Phys. Lett. 256(1996) p. 569.
        DWEIGHT = ONE
        IF (DOBOX.AND.(.NOT.(ACTVAOBTCH(IAOB1).AND.ACTVAOBTCH(IAOB2))))
     &    DWEIGHT = HALF

        CALL GETAOB(NPRF1,NCTF1,NTYP1,NHKT1,JMAX1,KHKT1,KCKT1,
     &              KCSQ1,KEXP1,KCCF(1,1),IAOB1)

        CALL GETAOB(NPRF2,NCTF2,NTYP2,NHKT2,JMAX2,KHKT2,KCKT2,
     &              KCSQ2,KEXP2,KCCF(1,2),IAOB2)

C       PBC comment: SCREENPRIM12 uses coordinates and a distance set
C        by GETR12, so the necessary translations have already
C        been put in. Hence, SCREENPRIM12 works without modification for
C        PBC. One thing, though: In the molecular case, one screens out
C        just a few terms, but in the PBC case screening means neglecting
C        a rapidly converging infinite series. This means that the error
C        estimates may be a little different. / ErikT
        CALL SCREENPRIM12(NPRF12,RFBTCH(1,IBTCH),WORK(KIPR12),.TRUE.,
     &                    NPRF1,NCTF1,NTYP1,NPRF2,NCTF2,NTYP2,
     &                    EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &                    CCFBTCH(KCCF(1,1)),CCFBTCH(KCCF(2,1)),
     &                    CCFBTCH(KCCF(1,2)),CCFBTCH(KCCF(2,2)),
     &                    R122(0,0),COOR1(1,0),COOR2(1,0),1,IPRINT)

        EXTQMAX = MAX( RFBTCH(4,IBTCH), EXTQMAX )
        SCRCD   = MAX( RFBTCH(5,IBTCH), SCRCD   )

        IF ( NTYP1+NTYP2 .GT. 2 ) NTYP12 = 2

        CALL SETETIJ(WORK(KETIJ),WORK(KPREF12),
     &               R12(1,0,0),R122(0,0),
     &               WORK(KIPR12),EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &               WORK(KCOOR12),WORK(KEXP12),WORK(KHPI),WORK(KP1),
     &               WORK(KP2),NPRF1,NHKT1,JMAX1,
     &               NPRF2,NHKT2,JMAX2,NPRF12,NHKT12,IELECTR,IPRINT)

c        find largest exponent of the ab distribution
         F=1D30
         DO I=1,NPRF12
            F = MIN(F,WORK(KEXP12+I-1))
         END DO
         EXPMAX = 1d0/F
c
        J1 = JMAX1 
        J2 = JMAX2 
        JMAX12 = NHKT12 - 1 

        CALL DZERO(WORK(KFTUVTMP+KSFTUV-KSFTUV0),NPRF12*NTUV)

        DO ICOMP1 = 1, NTYP1
         IF ( NTYP1 .EQ. 2 ) CALL SETJKKK(J1,KCSQ1,KCKT1,KHKT1,ICOMP1)

         DO ICOMP2 = 1, NTYP2 
          IF ( NTYP2 .EQ. 2 ) CALL SETJKKK(J2,KCSQ2,KCKT2,KHKT2,ICOMP2)

          CALL SETDMATCP(WORK(KDMATCP),DMAT,IAOB1,IAOB2,ICOMP1,ICOMP2,
     &         CSQ(KCSQ1,1),CSQ(KCSQ2,1),WORK(KIPR12),
     &         CCFBTCH(KCCF(ICOMP1,1)),CCFBTCH(KCCF(ICOMP2,2)),
     &         WORK(KDMATSP), WORK(KDTMP),KCKT1,KCKT2,KHKT1,KHKT2,
     &         NPRF1,NPRF2,NCTF1,NCTF2,NPRF12,IPRINT)
          
          CALL CONTRACTDMATE(WORK(KFTUVTMP+KSFTUV-KSFTUV0),
     &         WORK(KDMATCP),WORK(KETIJ),WORK(KPREF12),INDHER,LMNPWR,
     &         KCKT1,KCKT2,J1,J2,JMAX1,JMAX2,JMAX12,NPRF12,NTUV,IPRINT)

         END DO ! ICOMP2
        END DO ! ICOMP1

        CALL SETRFTUV(IFTUV(KSIFTUV),RFTUV(KSRFTUV),
     &         COOR1(1,0),COOR2(1,0),EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &         WORK(KIPR12),NPRF12,IBTCH,DWEIGHT,IPRINT)

        IFBTCH(1,IBTCH) = NPRF12
        IFBTCH(2,IBTCH) = NTFTUVPRIM+1

        RFBTCH(6,IBTCH) = RODBTCH( 7, IBTCH ) 

        KSIFTUV = KSIFTUV + NPRF12 * 2
        KSRFTUV = KSRFTUV + NPRF12 * 8 
        KSFTUV  = KSFTUV  + NPRF12 * NTUV

        NTFTUVBTCH = NTFTUVBTCH + 1
        NTFTUVPRIM = NTFTUVPRIM + NPRF12

       ELSE
         IFBTCH(1,IBTCH) = 0
       ENDIF ! Active batches
       END DO ! IBTCH

C       CALL SORTFTUV(FTUV(KSFTUV0),RFTUV(KSRFTUV0),
C     &               WORK(KFTUVTMP),WORK(KRFTUVTMP),
C     &               RFBTCH(1,KSTBTCH),IFBTCH(1,KSTBTCH),
C     &               NTUV,NTFTUVPRIM,NTFTUVBTCH,IPRINT)

       CALL RTRMAT(NTUV,NTFTUVPRIM,WORK(KFTUVTMP),FTUV(KSFTUV0))

       IF ( IPRINT .GE. 5 ) THEN 
        IF (PACKFTUV) THEN
          CALL PRIFTUV3(FTUV(KSFTUV0),INDHER,
     &                  ICLS,NTFTUVPRIM,JMAX12,NTUV)
        ELSE
          CALL PRIFTUV2(FTUV(KSFTUV0),ICLS,NTFTUVPRIM,JMAX12)
        END IF

        DO I=1,NTFTUVPRIM
         WRITE(LUPRI,'(3X,I5,8F12.7,I4)')
     &         I,(RFTUV(KSRFTUV0+8*(I-1)+J-1),J=1,8),
     &        IFTUV(KSIFTUV0+2*(I-1)+1)
        END DO      
       END IF
       DO I = 1, NTFTUVPRIM * NTUV
          FTUVMAX = MAX(FTUVMAX,ABS(FTUV(KSFTUV0+I-1)))
       END DO
  
       NFTUVBTCH = NFTUVBTCH + NTFTUVBTCH

       IFTUVCL( 1, ICLS ) = NTFTUVBTCH
       IFTUVCL( 2, ICLS ) = NTFTUVPRIM
       IFTUVCL( 8, ICLS)  = NTYP12

       RFTUVCL( 1, ICLS ) = SQRT( R122(0,0) ) 
       RFTUVCL( 2, ICLS ) = EXTQMAX
       RFTUVCL( 3, ICLS ) = FTUVMAX
       RFTUVCL( 4, ICLS ) = SCRCD
       RFTUVCL( 5, ICLS ) = EXPMAX
       RFTUVCL( 6, ICLS ) = RODBCLS ( 2, ICLS )

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
#include "cbieri.h"
#include "densfit.h"
#include "fck3rel.h"
#include "pbc.h"
      DIMENSION ETIJ(*),PREF12(NPRF12),
     &          COOR12(NPRF12,3),EXP12(NPRF12,3), 
     &          HPI(*), P1(NPRF12,3), P2(NPRF12,3),R12(3),I120(3),
     &          EXPBT1(NPRF1),EXPBT2(NPRF2),IPR12(2,NPRF12)
      REAL*8, ALLOCATABLE :: TEXP1(:), TEXP2(:)
      INTEGER X
      LOGICAL NOLOCS_SAVE

C
C     level of derivatives calculated
      IDL   = 0
      IF (IELECTR.EQ.1) THEN
        IDER  = NDERLHS
      ELSE
        IDER  = NDERRHS
      ENDIF

      JM1 = JMAX1 + IDER
      JM2 = JMAX2 + IDER

      CALL DZERO(ETIJ,NPRF12*(NHKT12+2*IDER)*(NHKT1+IDER)
     &           *(NHKT2+IDER)*3*2)

      IF (DOFIT .AND. IELECTR.EQ.2 .AND. RHSFIT) THEN
        NOLOCS_SAVE = NOLOCS
        NOLOCS = .TRUE.
      else if (dopbc .and. pbcactive .and. .not. pbc_sep_Rtensor) then
        NOLOCS_SAVE = NOLOCS
        NOLOCS = .TRUE.
      ENDIF

      CALL LOCSYM(R12,1,1,I120,IXYZ)

      IF (DOFIT .AND. IELECTR.EQ.2 .AND. RHSFIT) THEN
         NOLOCS = NOLOCS_SAVE
      else if (dopbc .and. pbcactive .and. .not. pbc_sep_Rtensor) then
         NOLOCS = NOLOCS_SAVE
      end if

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

      NCNT12 = 0

      ALLOCATE(TEXP1(NPRF12))
      ALLOCATE(TEXP2(NPRF12))

      CALL EXPCFT(ETIJ,NPRF12,JM1,JM2,COOR12,EXP12,I120,HPI,
     &            P1,P2,IRUTIN,IELECTR,JM1,JM2,NCNT12,IDL,
     &            ISCAL1,ISCAL2,IPRINT,
     &            TEXP1,TEXP2,AOVERP,BOVERP)

      DEALLOCATE(TEXP1)
      DEALLOCATE(TEXP2)
C

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
#include "densfit.h"
C
#include "pbc.h"
C
      DIMENSION DMAT(NDIMC,NDIMD),
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

      IOFFSET = 0
      IF (DOFIT.AND.RHSFIT) IOFFSET = NBASE

      DO IPRF12 = 1, NPRF12
       IPRF1 = IPR12(1,IPRF12)
       IPRF2 = IPR12(2,IPRF12)

        DO ICTF1 = 1, NCTF1
C        The last index of KNDXBTCH is hardwired to 1 (L.C.) at the
C        moment (june 2007), S.C. density contributions come later
         ISTAO1 = INDXBTCH(KNDXBTCH(ISPCOMP1,IAOB1,1) + ICTF1 - 1, 1 )

         DO ICTF2 = 1, NCTF2
          ISTAO2 = INDXBTCH(KNDXBTCH(ISPCOMP2,IAOB2,1) + ICTF2 - 1, 1 ) 
     &             - IOFFSET
          DO ISP1=1, KHKT1
           IAO1 = ISTAO1+ISP1-1
           DO ISP2=1, KHKT2
            IAO2 = ISTAO2 + ISP2 - 1
            if (dopbc .and. pbcactive) then
              FAC = 1.0D0
              DMATSP(ISP1,ISP2,IPRF12) =  DMATSP(ISP1,ISP2,IPRF12) 
     &             + FAC * CCF1(IPRF1,ICTF1)*CCF2(IPRF2,ICTF2)
     &             * DMAT(IAO1,IAO2)    
               
            else IF (DOFIT.AND.RHSFIT) THEN
            DMATSP(ISP1,ISP2,IPRF12) =  DMATSP(ISP1,ISP2,IPRF12) 
     &             + FAC * CCF1(IPRF1,ICTF1)*CCF2(IPRF2,ICTF2)
     &             *  DMAT(1,IAO2)
            ELSE
              DMATSP(ISP1,ISP2,IPRF12) =  DMATSP(ISP1,ISP2,IPRF12) 
     &             + FAC * CCF1(IPRF1,ICTF1)*CCF2(IPRF2,ICTF2)
     &             *  ( DMAT(IAO1,IAO2) + DMAT( IAO2, IAO1 ) )     
            ENDIF
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
        CALL OUTPUT(DMAT,1,NDIMC,1,NDIMD,NDIMC,NDIMD,1,LUPRI)

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
#include "implicit.h"
#include "priunit.h"
#include "fck3odb.h"
#include "cbifck3.h"
 
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
     &                    IPR12,NPRF12,IBTCH,DWEIGHT,IPRINT)
#include "implicit.h"
#include "priunit.h"

      DIMENSION IFTUV(2,NPRF12),RFTUV(8,NPRF12),COOR1(3),COOR2(3),
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
       RFTUV(8,I) = DWEIGHT

      END DO
      
      RETURN
      END

      SUBROUTINE SORTFTUV(FTUV0,RFTUV,FTUVTMP,RFTUVTMP,RFBTCH,
     &                    IFBTCH,NTUV,NTFTUVPRIM,NTFTUVBTCH,IPRINT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION FTUV0(NTFTUVPRIM,NTUV),FTUVTMP(NTUV,NTFTUVPRIM),
     &          RFTUV(8,*),RFTUVTMP(8,*),RFBTCH(6,*),IFBTCH(2,*)

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

        DO K=1,8
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
       DO K=1,8
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

C
C  Check nuclar attraction energy
C
      SUBROUTINE CHECKNUCATR(IFTUVCL,RFTUV,FTUV,DMAT,INDHER,JTOP,
     &                       NDMT,NCLS,WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "aovec.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "ericom.h"
#include "cbieri.h"
#include "aobtch.h"
#include "symmet.h"
#include "odbtch.h"
#include "odclss.h"
#include "eriao.h"

      DIMENSION IFTUVCL(10,*),RFTUV(*),FTUV(*),DMAT(*),WORK(*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP)

      RNA = 0.0D0

      DO ICL=1,NCLS

       NTFTUVPRIM = IFTUVCL(2,ICL)
       KSRFTUV    = IFTUVCL(4,ICL)
       KSFTUV     = IFTUVCL(5,ICL)

       JMAXF     = IFTUVCL(6,ICL)
       NTUV      = ( JMAXF + 1 ) * ( JMAXF + 2) * ( JMAXF + 3 ) / 6
       NTUV3     = ( JMAXF + 1 )**3

       KRNAI  = 1 
       KAHGTF = KRNAI + NTFTUVPRIM*NTUV3

       CALL DZERO(WORK(KRNAI),NTFTUVPRIM*NTUV3*NDMT)

       CALL GETRNAI(WORK(KRNAI),WORK(KAHGTF),NTFTUVPRIM,JMAXF, 
     &              RFTUV(KSRFTUV),8,IPRINT)

       RNA = RNA + ADDRNAI2(WORK(KRNAI),
     &               FTUV(KSFTUV),INDHER,JTOP,NTFTUVPRIM, JMAXF )

      END DO

      CALL HEADER('Nuclear attraction using FTUV',-1)
      WRITE (LUPRI,'(F30.16)') RNA

      CALL REFERENCE_NA(REFRNA,DMAT,NDMT,WORK,LWORK, IPRINT)

      WRITE (LUPRI,'(A,5X,E20.14)') 'difference:',REFRNA-RNA

      RETURN
      END

      FUNCTION ADDRNAI2(RNAI, FTUV, INDHER, JTOP, NPP12, JMAX12 )
#include "implicit.h"
C
      INTEGER T,U,V,TUV
      DIMENSION RNAI(NPP12,0:JMAX12,0:JMAX12,0:JMAX12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),FTUV(NPP12,*)

      ADDRNAI2=0.0D0

      DO T=0,JMAX12
       DO U=0,JMAX12
        DO V=0,JMAX12
         TUV = INDHER(T,U,V)
          DO I=1,NPP12
           ADDRNAI2 = ADDRNAI2 + RNAI(I,T,U,V)*FTUV(I,TUV)
         END DO
        END DO
       END DO
      END DO

      RETURN
      END 

      SUBROUTINE PRIFTUV3(FTUV,INDHER,ICL,NFTUV,JMAX,NTUV)
#include "implicit.h"
#include "priunit.h"       
#include "fck3odb.h"
      DIMENSION FTUV(NFTUV,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP)
      INTEGER T,U,V, TUV

      CALL HEADER('FTUV',-1)

      WRITE(LUPRI,'(/,7X,3(A,I5,A,3X),/)')
     &            'FTUV class number :',ICL,',',
     &            'NTFTUV :', NFTUV,',',
     &            'JMAXF  :', JMAX,' '

      WRITE(LUPRI,'(7X,A)') 
     &            'I    T    U    V     FTUV(I,T,U,V)'

      DO I = 1, NFTUV

       DO T = 0, JMAX
        DO U = 0, JMAX - T 
         DO V = 0, JMAX - T - U  
          TUV = INDHER(T,U,V)
          WRITE(LUPRI,'(3X,4I5,F16.8)') I,T,U,V,FTUV(I,TUV)
         END DO
        END DO
       END DO
      END DO

      RETURN
      END
