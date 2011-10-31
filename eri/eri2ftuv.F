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
C   maximum NATOMS^2 classes: NFTUVCL <= NATOMS^2*NKHMAX
C   based on one-center/two-center and center of the C
C   criteria
C 
C  - Class information is in IFTUVCL( 6, NFTUVCL )
C      1 - NTFTUV   = number of FTUV's in the class
C      2 - KRFTUV   = pointer to the real information about
C                    the first FTUV of this class in RFTUV
C                     ( RFTUV(KRFTUV) )
C      3 - KSFTUV   = pointer to the FTUV of this 
C                     class in FTUV ( FTUV(KSFTUV) )
C      4 - JMAXF    = total angular momentum in this class
C      5 - NHKTF    = JMAXF + 1
C      6 - NCNTF    = 1 for one-center (original) OD
C                   = 2 for two-center (original) OD
C      7 - ICNTFC   = center of the C basis functions 
C      8 - ICNTFD   = center of the D basis functions 
C      9 - KIFTUV   =  pointer to the integer information about
C                    the first FTUV of this class in IFTUV
C                     ( IFTUV(KIFTUV) )
C     10 - NODF   = number of OD contributing (contracted)
C                   to this class
C                     
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
C    of array FTUV(NTFTUV,NDMT,0:JMAXF,0:JMAXF,0:JMAXF)
C
C  - Real information about contracted OD corresponding 
C        to FTUV:  RODQ(4,KODBC2)
C  
C      1 - x-coordinate of the center
C      2 - y-coordinate of the center
C      3 - z-coordinate of the center
C      4 - extent (from the center)
C
C -------------------------------------------------------------

      SUBROUTINE GETFTUV( KIFTUVCL,KRFTUVCL,KIFTUV,KRFTUV,KFTUV,
     &                    KRODQ,DMAT,ISYMDM,NDMT,
     &                    IODCL2,IODBC2,RODBC2,IODPP2,RODPP2,
     &                    KFREE,LFREE,CCFBT,INDXBT,
     &                    IELECTR,WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "cbieri.h"
#include "odbtch.h"
#include "odclss.h"
#include "ericom.h"
#include "inforb.h"
#include "nuclei.h"
C
      DIMENSION DMAT(NBAST,NBAST,NDMT),ISYMDM(NDMT),WORK(LWORK),
     &          CCFBT(*), INDXBT(*),  IODCL2(NODCL2,NITCL),
     &          IODBC2(NODBC2,NITBC), RODBC2(NODBC2,NRTBC),
     &          IODPP2(NODPP2,NITPP), RODPP2(NODPP2,NRTPP)

C
C   Ensure the density matrix is symmetric
C   ======================================
C
      DO IDMT=1,NDMT
       DO I=1,NBAST
        DO J=1,I-1
         DMAT(I,J,IDMT)=0.5D0*(DMAT(I,J,IDMT)+DMAT(J,I,IDMT))
         DMAT(J,I,IDMT)=DMAT(I,J,IDMT)
        END DO
       END DO
      END DO
C
C   Find out the number of different FTUV classes
C   =============================================
C   ( which equals to the (max J) + 1 = NHKMAX )
 

      NHKMAX  = 0 
      NPRFMAX = 0 

      DO 100 KODCL2 = 1, NODCL2

       NHKTC  = IODCL2(KODCL2,3)
       NHKTD  = IODCL2(KODCL2,4)
       NHKTCD = NHKTC + NHKTD - 1
 
       NPRFC  = IODCL2(KODCL2,15)
       NPRFD  = IODCL2(KODCL2,16)
       NPRFCD = NPRFC*NPRFD

       NHKMAX = MAX( NHKMAX, NHKTCD )
       NPRFMAX = MAX(NPRFMAX,NPRFCD)

 100  CONTINUE

C
C
C     Allocate space for FTUV class information 
C     =========================================
C     ( and the maximum space needed for 
C       RFTUV (but would be reduced due to the 
C       screening later );
C       

      NQLEVEL=1

      NFTUVCL = NQLEVEL*(NATOMS**2)*NHKMAX 

      KQLEVEL  = KFREE
      KRODQ    = KQLEVEL  +  NQLEVEL
      KIFTUVCL = KRODQ    +  4*NODBC2
      KRFTUVCL = KIFTUVCL + 10*NFTUVCL/IRAT + 1
      KIFTUV   = KRFTUVCL +  4*NFTUVCL 
      KRFTUV   = KIFTUV   +  2*NODBC2*NPRFMAX/IRAT + 1
      KLAST    = KRFTUV   +  7*NODBC2*NPRFMAX

C
C     Allocate space for LMNPWR,LMNSYM
C     ================================

C      NHKMAX = NHTYP + MAXDER
      KCKMAX = NHKMAX*(NHKMAX + 1)/2

      KLMNPWR = KLAST
      KLMNSYM = KLMNPWR + ( 3*NHKMAX*KCKMAX - 1)/IRAT + 1
      KLAST   = KLMNSYM + (16*NHKMAX*KCKMAX - 1)/IRAT + 1

C      WRITE(LUPRI,*) 'KIFTUVCL...',KIFTUVCL,KRFTUV,
C     &                KLMNPWR,KLMNSYM

C
C     Allocate space for CSQ
C     =======================
C 
      CALL SPHDIM

      KCSQ = KLAST
      KLAST = KCSQ   + NCSQ1*NCSQ2

      LNEED = KLAST - KFREE

      IF ( LNEED .GT. LFREE ) THEN
         WRITE (LUPRI,'(/,2(2X,A,I5,/))')
     &      'Memory requirements total     : ',LNEED,
     &      'Memory available              : ',LWORK
         CALL QUIT('Insufficient memory in GETFTUV, 0')
      END IF

      KWRK=KLAST
      LWRK=LFREE-LNEED

      CALL DZERO(WORK(KFREE),KWRK-KFREE)

C
C     Spherical-harmonic transformation matrices
C     ==========================================
C
      CALL SPHMAT(WORK(KCSQ),NCSQ1,NCSQ2,MAXDER,GDER,
     &            WORK(KWRK),LWRK,IPRINT)

C
C     Prepare index array LMNPWR 
C     ==========================
C     ( to get carteesian indexes ijk corresponding to 
C       a certain spherical component )
C

      CALL LMNPRP(WORK(KLMNPWR),WORK(KLMNSYM),IPRINT)

      CALL GETFTUV1( WORK(KIFTUVCL),WORK(KRFTUVCL),
     &               WORK(KIFTUV),WORK(KRFTUV),WORK(KRODQ),
     &               KFTUV_TMP,DMAT,NDMT,
     &               IODCL2,IODBC2,RODBC2,IODPP2,RODPP2,
     &               WORK(KLMNPWR),WORK(KLMNSYM),WORK(KCSQ),
     &               CCFBT,INDXBT,NPRFMAX,
     &               WORK(KQLEVEL),NQLEVEL,
     &               IELECTR,WORK,KWRK,LWRK,IPRINT)

C      CALL PRIFULLIFTUV(WORK(KIFTUVCL),WORK(KIFTUV),NFTUVCL)

      CALL REORDERFTUV(WORK(KIFTUVCL),KIFTUVCL,
     &                KRFTUVCL,KIFTUV,KRFTUV,KFTUV,
     &                KFTUV_TMP,NDMT,KFREE,LFREE,WORK,IPRINT)

      IF ( IELECTR .EQ. 1 ) CALL CHECK_RNAI(WORK(KIFTUVCL),
     &                         WORK(KRFTUV),WORK(KFTUV),
     &                         DMAT,NDMT,WORK(KFREE),LFREE,IPRINT)

C      CALL PRIFULLIFTUV(WORK(KIFTUVCL),WORK(KIFTUV),NFTUVCL)

      CALL SUMMARY_FTUV(WORK(KIFTUVCL),WORK(KRFTUVCL),
     &                  WORK(KIFTUV),WORK(KRFTUV),WORK(KFTUV),
     &                  NDMT,KIFTUVCL,KRFTUVCL,KIFTUV,KRFTUV,KFTUV,
     &                  IODCL2,KFREE,LFREE,IPRINT)

      RETURN 
      END 

      SUBROUTINE GETFTUV1( IFTUVCL,RFTUVCL,IFTUV,RFTUV,RODQ,
     &                     KFTUV_TMP,DMAT,NDMT,
     &                     IODCL2,IODBC2,RODBC2,IODPP2,RODPP2,
     &                     LMNPWR,LMNSYM,CSQ,
     &                     CCFBT,INDXBT,NPRFMAX,
     &                     QLEVEL,NQLEVEL,
     &                     IELECTR,WORK,KWRK,LWORK,IPRINT)
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
#include "erithr.h"
#include "nuclei.h"

C
      DIMENSION IFTUVCL(10,*),RFTUVCL(4,*),IFTUV(*),RFTUV(*),RODQ(4,*),
     &          DMAT(NBASE,NBASE,NDMT),
     &          IODCL2(NODCL2,NITCL),
     &          IODBC2(NODBC2,NITBC), RODBC2(NODBC2,NRTBC),
     &          IODPP2(NODPP2,NITPP), RODPP2(*),
     &          LMNPWR(KCKMAX,NHKMAX,3),LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          CSQ(NCSQ1,NCSQ2),CCFBT(*),INDXBT(MXSHEL*MXCONT,0:7),
     &          QLEVEL(*),WORK(LWORK)    

      LOGICAL ONECENTNEW

      QLEVEL(1) =  0.0D0
C      QLEVEL(2) =  2.0D0
C      QLEVEL(3) =  4.0D0


      KHKMAX   = 2*NHKMAX + 1

C
C     Make the temporary allocations 
C     ==============================
C

      KWEIGHT  = KWRK
      KECOEF   = KWEIGHT  + NPRFMAX
      KCOOR12  = KECOEF   + NPRFMAX * (NHKMAX**3)*3*2
      KEXP12   = KCOOR12  + NPRFMAX * 3
      KPREF12  = KEXP12   + NPRFMAX * 3
      KHPI     = KPREF12  + NPRFMAX
      KPC      = KHPI     + NPRFMAX
      KPD      = KPC      + NPRFMAX * 3
      KDMATSP  = KPD      + NPRFMAX * 3
      KDMATCP  = KDMATSP  + NPRFMAX * (KHKMAX**2)
      KDTMP    = KDMATCP  + NPRFMAX * (KCKMAX**2)
      KRNAI    = KDTMP    + KCKMAX  * KHKMAX
      KAHGTF   = KRNAI    + NPRFMAX * (NHKMAX**3)

      KFTUV1   = KAHGTF   +       2 * (NHKMAX**3) 
      KFTUV_TMP  = KFTUV1     + NPRFMAX * (NHKMAX**3)
      KFTUV_TMP2 = KFTUV_TMP  + NPRFMAX * (NHKMAX**3)*NDMT
      KLAST      = KFTUV_TMP2 + NPRFMAX * (NHKMAX**3)*NDMT

C
C  (  LNEED2 to store NDMT new FTUV's )
C     temporary locations for FTUV
C
      LWORK       =  LWORK - (KLAST - KWRK)
      LNEED2      =  NPRFMAX * (NHKMAX**3)*NDMT

C
C     Loop over all classes of OD's for electron 2 
C     ============================================
C      
      KCLSF = 0

      DO 800 IQLEVEL=1,NQLEVEL

      QMIN=QLEVEL(IQLEVEL)
      IF (IQLEVEL .LT. NQLEVEL) THEN 
       QMAX = QLEVEL(IQLEVEL+1)
      ELSE
       QMAX = 1.0D10
      END IF

C      WRITE(LUPRI,*) 'IQLEVEL',IQLEVEL,QMIN,QMAX

      DO 700 ICNTFC=1,NATOMS
      DO 600 ICNTFD=ICNTFC,NATOMS

       NHKTF = 0

       DO 500 KODCL2 = 1, NODCL2

        NODTCD = IODCL2(KODCL2,1)
        KODSCD = IODCL2(KODCL2,2)

        NHKTC  = IODCL2(KODCL2,3)
        NHKTD  = IODCL2(KODCL2,4)
        NHKTCD = NHKTC + NHKTD - 1

        KCSQC = KSQADR(NHKTC-1,0,0,0,0)
        KCSQD = KSQADR(NHKTD-1,0,0,0,0)
      
        KHKTC  = IODCL2(KODCL2,5)
        KHKTD  = IODCL2(KODCL2,6)
        KCKTC  = IODCL2(KODCL2,7)
        KCKTD  = IODCL2(KODCL2,8)

        NCNTCD = IODCL2(KODCL2,13)

        ICMATC = IODCL2(KODCL2,19)
        ICMATD = IODCL2(KODCL2,20)

        JMAXC  = NHKTC - 1 
        JMAXD  = NHKTD - 1 
        JMAXCD = JMAXC + JMAXD 

        JMAX1  = NHKTC - 1 
        JMAX2  = NHKTD - 1 
        JMAX12 = JMAX1 + JMAX2 

        NPRFC  = IODCL2(KODCL2,15)
        NPRFD  = IODCL2(KODCL2,16)
        NPRFCD = NPRFC*NPRFD

        NCTFC  = IODCL2(KODCL2,17)
        NCTFD  = IODCL2(KODCL2,18)

        NPP12  = NPRFCD

        NTUV   = NHKTCD**3

C
C aspg, 2006-06-27 making ftnchek happy
C
        LWRK = LWORK
C
C     Loop over all of OD's in the class
C     ==================================

         DO 400 KODBC2=KODSCD,(KODSCD+NODTCD-1)

          IF ( LNEED2 .GT. LWORK ) THEN
           WRITE (LUPRI,'(/,2(2X,A,I5,/))')
     &      'Memory requirements total LNEED2 : ',LNEED2,
     &      'Memory available                 : ',LWRK
           CALL QUIT('Insufficient memory in GETFTUV1, 2')
          END IF

          SCR  = RODBC2(KODBC2,1)

          IC  = IODBC2(KODBC2,2)
          ID  = IODBC2(KODBC2,3)

          FAC = 1.0D0 
          IF ( IC .NE. ID ) FAC = 2.0D0

          IF ( IPRINT .GE. 4 ) THEN 
           WRITE(LUPRI,*) '  KODCL2 :', KODCL2,
     &                   '  NODTCD :', NODTCD,
     &                   '  KODBC2 :', KODBC2,
     &                   '  NPRFCD :', NPRFCD,
     &                   '  NHKTC  :', NHKTC,
     &                   '  NHKTD  :', NHKTD,
     &                   '  IC,ID  :', IC,ID

          END IF
 
          ICNTC = NCNTBT(IC)
          ICNTD = NCNTBT(ID)

C          WRITE(LUPRI,'(3X,A,4I5,6F12.6)')
C     &         '@COND1: ',ICNTFC,ICNTC,ICNTFD,ICNTD,
C     &               CORXBT(IC), CORYBT(IC), CORZBT(IC),
C     &               CORXBT(ID), CORYBT(ID), CORZBT(ID)
          IF (((ICNTFC.EQ.ICNTC).AND.(ICNTFD.EQ.ICNTD)) .OR.
     &        ((ICNTFC.EQ.ICNTD).AND.(ICNTFD.EQ.ICNTC))) THEN

C          WRITE(LUPRI,*) '@KODBC2:',KODBC2

C
C     Get the hermite expansion coefficients Etij
C     ===========================================

          CALL  GETETIJ(KODBC2,RCD,WORK(KECOEF),
     &            WORK(KCOOR12),WORK(KEXP12),WORK(KPREF12),
     &            WORK(KHPI),WORK(KPC),WORK(KPD),
     &            IODBC2,RODBC2,IODPP2,RODPP2,
     &            NPRFC, NPRFD, IELECTR, IPRINT)

C
C    Loop over the density matrices
C    ==============================
C

          DO 200, IDMT = 1, NDMT

C
C    Transform density matrix to cartesian primitive basis
C    =====================================================
C
C    The corresponding submatrix from the density matrix.
C     We assume DMAT is full without symmetry reduction.
C

           CALL GETDMATCP( DMAT(1,1,IDMT),
     &          WORK(KDMATCP), WORK(KDMATSP), WORK(KDTMP),
     &          CSQ(KCSQC,1), CSQ(KCSQD,1),
     &          CCFBT(ICMATC),CCFBT(ICMATD),INDXBT,
     &          IC,ID,IPRINT )            
  
C
C    Contract the density in cartesian basis with Etij 
C    ================================================
C

           CALL GETFTUV2( WORK(KFTUV1),WORK(KDMATCP),LMNPWR,
     &                   WORK(KECOEF),WORK(KPREF12),
     &                   IELECTR,FAC,IPRINT)

C
C  First copy it to a temporary location ...
C

           KOFF = (IDMT-1)*NPRFCD*NTUV

           CALL DCOPY(NPRFCD*NTUV,WORK(KFTUV1),1,
     &                           WORK(KFTUV_TMP2+KOFF),1)

           IF ( IPRINT .GE. 15 ) THEN
            CALL PRIFTUV(WORK(KFTUV1),KODBC2,NPRFCD,JMAX12)
            CALL PRIFTUV(WORK(KFTUV_TMP2+KOFF),KODBC2,NPRFCD,JMAX12)
           END IF

 200      CONTINUE

C
C     Check if this is a new FTUV class
C     =================================

C          WRITE(LUPRI,*) '@COND2: ',NHKTF,NHKTCD,ICNTC,ICNTD

          IF (NHKTF.NE.NHKTCD) THEN

           NHKTF = NHKTCD

           ITFTUV = 0 

            IF ( KCLSF .EQ. 0 ) THEN
             KRFTUV = 1 
             KSFTUV = 1
             KIFTUV = 1 
             KCLSF  = 1 
            ELSE
             NTFTUV_PREV = IFTUVCL( 1, KCLSF )
             KRFTUV_PREV = IFTUVCL( 2, KCLSF )
             KSFTUV_PREV = IFTUVCL( 3, KCLSF )
             KIFTUV_PREV = IFTUVCL( 9, KCLSF )
             NHKTF_PREV  = IFTUVCL( 5, KCLSF )

             NTUV_PREV = (NHKTF_PREV)**3
             KRFTUV = KRFTUV_PREV + NTFTUV_PREV * 7
             KSFTUV = KSFTUV_PREV + NTFTUV_PREV * NDMT * NTUV_PREV
             KIFTUV = KIFTUV_PREV + NTFTUV_PREV * 2 
             IF ( NTFTUV_PREV .GT. 0 ) THEN
                KCLSF = KCLSF+1
             END IF
            END IF

            IFTUVCL( 2, KCLSF ) = KRFTUV
            IFTUVCL( 3, KCLSF ) = KSFTUV
         
            IFTUVCL( 4, KCLSF ) = NHKTF - 1 
            IFTUVCL( 5, KCLSF ) = NHKTF

            NCNTF=1
            IF (ICNTFC.NE.ICNTFD) NCNTF=2

            IFTUVCL( 6, KCLSF ) = NCNTF
            IFTUVCL( 7, KCLSF ) = ICNTFC
            IFTUVCL( 8, KCLSF ) = ICNTFD
            IFTUVCL( 9, KCLSF ) = KIFTUV
            IFTUVCL(10, KCLSF ) = 0

            RFTUVCL( 1, KCLSF ) = RCD

C            WRITE(LUPRI,*) '@FTUV-class:', KCLSF,
C     &                  ' KRFTUV:',KRFTUV,
C     &                  ' KSFTUV:',KSFTUV, 
C     &                  ' NDMT:',NDMT, 
C     &                  ' NHKTF:', NHKTF,
C     &                  ' NCENT:', NCNTF,
C     &                  ' ICENT:', ICNTF
          END IF

C
C    Screen and save only if higher than threshold
C    =============================================
C



          CALL SCREENFTUV(RODQ,RFTUV(KRFTUV),IFTUV(KIFTUV),
     &                  WORK(KFTUV_TMP+KSFTUV-1),WORK(KFTUV_TMP2),
     &                  ITFTUV,WORK(KWEIGHT),KODBC2,NDMT,NSAVED,
     &                  WORK(KEXP12),IC,ID,QMIN,QMAX,SCR,
     &                  CCFBT(ICMATC),CCFBT(ICMATD),
     &                  DMAT,INDXBT,WORK,LWRK)


C          WRITE(LUPRI,*) '@KRFTUV,ITFTUV,ICL',NSAVED, KRFTUV, ITFTUV,NHKTCD

C          CALL PRIIFTUV(IFTUV(KIFTUV))

          IFTUVCL( 1, KCLSF ) = ITFTUV
 
          IF (NSAVED .GT. 0) IFTUVCL(10,KCLSF)=IFTUVCL(10,KCLSF)+1 
           

          LSAVED = NSAVED*NDMT*NTUV
          LWRK   = LWRK - LSAVED

          KFTUV_TMP2 = KFTUV_TMP2 + LSAVED

          IF ( IPRINT .GE. 15 ) THEN 
           WRITE(LUPRI,'(3X,A,I5)') 
     &      'Screening for FTUV class :', NHKTCD
           WRITE(LUPRI,'(3X,A,2I5)') 
     &      'Screening OD class and batch :', KODCL2,KODBC2
           WRITE(LUPRI,'(3X,A,I10)') 
     &      'Total number of primitive FTUV s :',
     &       NPRFCD*NODTCD
           WRITE(LUPRI,'(3X,A,I10,/)')
     &      'Number of saved FTUV s           :',
     &       NSAVED 
           WRITE(LUPRI,*) 'kftuv_tmp : ',KFTUV_TMP,KSFTUV,ITFTUV
           CALL PRIFTUV(WORK(KFTUV_TMP+KSFTUV-1+ITFTUV-NSAVED),
     &                 KODBC2,NSAVED,JMAX12)
           DO I=1,NSAVED
            WRITE(LUPRI,'(3X,A,I5,7F10.6,/)') '@COORQ:',
     &       ITFTUV-NSAVED+I, 
     &       (RFTUV(KRFTUV+7*(ITFTUV-NSAVED+I-1)+J-1),J=1,7)
           END DO
          END IF

         END IF

 400    CONTINUE

 500   CONTINUE

 600  CONTINUE
 700  CONTINUE
 800  CONTINUE

      NFTUVCL = KCLSF

      RETURN
      END 

      SUBROUTINE GETETIJ(KODBC2,RCD,ETIJ,
     &                   COOR12,EXP12,PREF12,HPI,P1,P2,
     &                   IODBC2,RODBC2,IODPP2,RODPP2,
     &                   NPRFA, NPRFB, IELECTR,IPRINT)
C 
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "eriao.h"
#include "cbieri.h"
#include "aobtch.h"
#include "symmet.h"
#include "iratdef.h"
#include "odbtch.h"
#include "odclss.h"
C
      DIMENSION IODCL2(NODCL2,NITCL),
     &          IODBC2(NODBC2,NITBC), RODBC2(NODBC2,NRTBC),
     &          IODPP2(NODPP2,NITPP), RODPP2(NODPP2,NRTPP),
     &          ETIJ(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3,*),
     &          COOR12(NPP12,3), EXP12(NPP12,3), PREF12(NPP12), 
     &          HPI(NPP12), P1(NPP12,3), P2(NPP12,3)

      INTEGER X

      CALL DZERO(ETIJ(1,0,0,0,1,1),NPP12*(JMAX1+JMAX2)
     &            *(JMAX1+1)*(JMAX2+1)*3*2)

      IF ( IPRINT .GE. 5 ) THEN
        CALL HEADER('OUTPUT from GETETIJ',-1)
      END IF
C
C     Collect coordinates and exponents
C     =================================
C
         
      I1 = IODBC2(KODBC2,2)
      I2 = IODBC2(KODBC2,3)

      COOR12(1,1) = CORXBT(I1) -  CORXBT(I2) 
      COOR12(1,2) = CORYBT(I1) -  CORYBT(I2) 
      COOR12(1,3) = CORZBT(I1) -  CORZBT(I2) 

      RCD = SQRT(COOR12(1,1)**2+COOR12(1,2)**2+COOR12(1,3)**2)

      CALL LOCSYM(COOR12,1,NPP12,I120,IXYZ)

      DO 100 I=1, NPP12
        DO X=1,3 
          COOR12(I,X) = COOR12(1,X)
        END DO
 100  CONTINUE         

      IF ( IPRINT .GE. 5 ) THEN
        CALL HEADER('AO batches numbers :',5)
        WRITE(LUPRI,'(2(5X,A,I5,/))') 'AO batch I1: ', I1,
     &                              'AO batch I2: ', I2
        CALL HEADER('Coordinates : ',5)
        WRITE(LUPRI,'(3(5X,A,3F12.6,/))') 
     &          'I1   :', CORXBT(I1), CORYBT(I1), CORZBT(I1),
     &          'I2   :', CORXBT(I2), CORYBT(I2), CORZBT(I2),
     &          'I1-I2:', ( COOR12(1,I), I=1,3)
        WRITE(LUPRI,*) 'NPRF1, NPRF2, NPP12: ', 
     &                  NPRFA, NPRFB, NPP12
      END IF

      I = 1
      DO 300 IPR1=1, NPRFA
       DO 200 IPR2=1, NPRFB
        
         KEXP1 = KEXPBT(I1) - 1
         KEXP2 = KEXPBT(I2) - 1 

         EXP1  = EXPBT(KEXP1 + IPR1)
         EXP2  = EXPBT(KEXP2 + IPR2)

         IF ( IPRINT .GE. 5 ) THEN
           WRITE(LUPRI,*) 'IPR1, IPR2 : ', IPR1, IPR2
           WRITE(LUPRI,*) 'KEXP1, KEXP2 :',KEXP1, KEXP2
           WRITE(LUPRI,*) 'EXP1,  EXP2  :',EXP1, EXP2
         END IF

         EXPPI = 1.0D0 / ( EXP1 + EXP2 )

         EXP12(I,1) = EXP1 * EXPPI
         EXP12(I,2) = EXP2 * EXPPI
         EXP12(I,3) = EXPPI

         R2 = 0.0D0

         DO 150 X = 1,3
           P1(I,X) = - EXP12(I,2) * COOR12(I,X) 
           P2(I,X) =   EXP12(I,1) * COOR12(I,X)
           R2 = R2 + COOR12(I,X)**2
 150     CONTINUE

         PREF12(I) = exp( -EXP1*EXP2*EXPPI*R2 )

         I = I + 1
 200    CONTINUE
 300   CONTINUE

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

      IF ( IPRINT .GE. 20 ) THEN
        CALL HEADER( 'In GETETIJ',-1 )

        WRITE (LUPRI,'(A,I12)')
     &             'CALLING EXPCFT for ODBACH: ', KODBC2
        WRITE (LUPRI,'(7X,A,I6,3X,A,I6)')  
     &               'I1 =',I1,'I2 = ', I2
        WRITE (LUPRI,'(20X,A,2(/),2X,3(10X,A,6X),/)')
     &               'Coordinates','X','Y','Z'
        WRITE (LUPRI,'(2(2X,3F16.8,/))')
     &                CORXBT(I1), CORYBT(I1), CORZBT(I1),
     &                CORXBT(I2), CORYBT(I2), CORZBT(I2)

        CALL HEADER( 'In GETETIJ: Exponents',-1 )

        DO I=1,NPP12
          WRITE( LUPRI, '(2F16.8)') EXP12(I,1)/EXP12(I,3), 
     &                              EXP12(I,2)/EXP12(I,3)
        END DO
      END IF

C
C aspg, 2006-06-27 making ftnchek happy
C
C      CALL EXPCFT(ETIJ,NPP12,JMAX1,JMAX2,COOR12,EXP12,I120,HPI,
C     &                  P1,P2,IRUTIN,IELECTR,JMAX1,JMAX2,NCNT12,IDL,
C     &                  ISCAL1,ISCAL2,IPRINT)
      CALL EXPCFT(ETIJ,NPP12,JMAX1,JMAX2,COOR12,EXP12,I120,HPI,
     &                  P1,P2,IRUTIN,IELECTR,JMAX1,JMAX2,NCNT12,IDL,
     &                  ISCAL1,ISCAL2,IPRINT,
     &                  TEXPA,TEXPB,AOVERP,BOVERP)
C

      IF ( IPRINT .GE. 20 ) THEN
        CALL HEADER( 'After GETETIJ',-1 )

        WRITE (LUPRI,'(7X,A,I6,3X,A,I6)')  
     &               'I1 =',I1,'I2 = ', I2
        WRITE (LUPRI,'(20X,A,2(/),2X,3(10X,A,6X),/)')
     &               'Coordinates I1-I2','X','Y','Z'
        WRITE (LUPRI,'(2X,3F16.8,/)')
     &           ( COOR12(1,I), I=1,3)
      END IF

      RETURN
      END

      SUBROUTINE GETFTUV2(FTUV,DMATCP,LMNPWR,ETIJ,
     &                    PREFCD,IELECTR,FAC,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "ericom.h"
#include "eriao.h"
#include "cbieri.h"
#include "aobtch.h"
#include "symmet.h"
#include "odbtch.h"
#include "odclss.h"
C
      DIMENSION LMNPWR(KCKMAX,NHKMAX,3),
     &  DMATCP((JMAXC+1)*(JMAXC+2)/2,(JMAXD+1)*(JMAXD+2)/2,*),
     &          FTUV(0:JMAXCD,0:JMAXCD,0:JMAXCD,NPRFCD),
     &          ETIJ(NPRFCD,0:JMAXC+JMAXD,0:JMAXC,0:JMAXD,3,*),
     &          PREFCD(NPRFCD)

      INTEGER T,U,V

      CALL DZERO(FTUV(0,0,0,1),NPRFCD*(JMAXCD+1)**3)

      IF ( IPRINT .GE. 20 ) THEN
         CALL HEADER('Exponencial prefactors',-1)
      
         DO I=1,NPRFCD
           WRITE(LUPRI,'(5X,I3,F16.8)') I, PREFCD(I)
         END DO
         WRITE(LUPRI,'(/)')
         WRITE(LUPRI,*) 'JMAXC, JMAXD, JMAXCD : ',JMAXC,JMAXD,JMAXCD
      END IF
      
      KHKTC = 2*JMAXC+1
      KCKTC = (JMAXC+1)*(JMAXC+2)/2
      KHKTD = 2*JMAXD+1
      KCKTD = (JMAXD+1)*(JMAXD+2)/2

      DO T = 0, JMAXCD
       DO U = 0, JMAXCD - T 
        DO V = 0, JMAXCD - T - U  

         DO KCC = 1, KCKTC

          LC = LMNPWR(KCC,JMAXC+1,1)
          MC = LMNPWR(KCC,JMAXC+1,2)
          NC = LMNPWR(KCC,JMAXC+1,3)

          IF ( IPRINT .GE. 20 ) THEN
            WRITE(LUPRI,'(7X,A,5I5)') 'KCC, JMAXC, LC, MC, NC', 
     &                                 KCC, JMAXC, LC, MC, NC
          END IF 

          DO KCD = 1, KCKTD
 
           LD = LMNPWR(KCD,JMAXD+1,1)
           MD = LMNPWR(KCD,JMAXD+1,2)
           ND = LMNPWR(KCD,JMAXD+1,3)

           IF ( IPRINT .GE. 20 ) THEN
             WRITE(LUPRI,'(7X,A,5I5)') 'KCD, JMAXD, LD, MD, ND', 
     &                                  KCD, JMAXD, LD, MD, ND
           END IF
 
C           IF ( T+U+V .LE. JMAXCD ) THEN
            DO I=1,NPRFCD
             F_TMP = FAC* DMATCP(KCC,KCD,I)* PREFCD(I)
     &                                 * ETIJ(I,T,LC,LD,1,1)
     &                                 * ETIJ(I,U,MC,MD,2,1)
     &                                 * ETIJ(I,V,NC,ND,3,1) 

 
             FTUV(T,U,V,I) = FTUV(T,U,V,I) + F_TMP 

C
C  It seems the correct sign (-1)^(T+U+V) is 
C   included already in the generation of 
C   the Etij coefficients, so we skip it here.

             IF ( IPRINT .GE. 20 ) THEN
              WRITE(LUPRI,'(/,5X,A,3I6,/)') ' T U V : ', T, U, V 
              WRITE(LUPRI,'(3(5XA,F16.10,/))') 
     &         'ETIJ(1,T,LC,LD,1,1) = ', ETIJ(1,T,LC,LD,1,1),
     &         'ETIJ(1,U,MC,MD,2,1) = ', ETIJ(1,U,MC,MD,2,1),
     &         'ETIJ(1,V,NC,ND,3,1) = ', ETIJ(1,V,NC,ND,3,1)
              WRITE(LUPRI,'(2(3X,A,F16.10),/)') 'F_TMP :', F_TMP,
     &                'DMATCP(KCC,KCD,I)',DMATCP(KCC,KCD,I)
             END IF 
            END DO
C          END IF 
          END DO
         END DO
        END DO
       END DO
      END DO

      RETURN
      END

C
C  DMAT:    - AO is in Spherical Contracted
C        CP - Cartesian Primitive
C        SP - Spherical Primitive
C
C  So far we do not take care about ISYMDM
C  DMAT is assumed to be full dimension without any 
C  symmetry reduction
C
      SUBROUTINE GETDMATCP( DMAT, DMATCP, DMATSP, DTMP,
     &                      CSQC, CSQD, CCFC, CCFD, INDXBT,
     &                      IC, ID, IPRINT)   
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "aovec.h"
#include "aobtch.h"
#include "symmet.h"
#include "ericom.h"
C
      DIMENSION DMAT(NBASE,NBASE),
     &          DMATCP(KCKTC,KCKTD,NPRFCD),
     &          DMATSP(KHKTC,KHKTD,NPRFCD),
     &          DTMP(KCKTC,KHKTD),
     &          CSQC(KHKTC,KCKTC), CSQD(KHKTD,KCKTD),
     &          CCFC(NPRFC,NCTFC), CCFD(NPRFD,NCTFD),
     &          INDXBT(MXSHEL*MXCONT,0:7)

      IF ( IPRINT .GE. 10 ) THEN
        CALL HEADER('Output from GETDMATCP',-1)
      END IF 

C
C ( For the moment we do not consider the symmetry )
C

      CALL DZERO(DMATCP, KCKTC*KCKTD*NPRFCD )
      CALL DZERO(DMATSP, KHKTC*KHKTD*NPRFCD )
      CALL DZERO(DTMP, KHKTC*KHKTD )

      IREPC = 0
      IREPD = 0

      IPRFCD = 1 

      DO 500 IPRFC = 1, NPRFC
       DO 400 IPRFD = 1, NPRFD

        DO 300 ICTFC = 1, NCTFC
         DO 200 ICTFD = 1, NCTFD

          IAOC = INDXBT(KNDXBT(IC) - 1 + ICTFC,IREPC)
          IAOD = INDXBT(KNDXBT(ID) - 1 + ICTFD,IREPD)

          DO 150 ISPC=1, KHKTC
C            CALL DAXPY( KHKTD,CCFC(IPRFC,ICTFC)*CCFD(IPRFD,ICTFD),
C     &                  DMAT(IAOC+ISPC-1,IAOD),1, 
C     &                  DMATSP(ISPC,1,IPRFCD), 1 )
           DO 120 ISPD=1, KHKTD
             DMATSP(ISPC,ISPD,IPRFCD) =  DMATSP(ISPC,ISPD,IPRFCD)
     &            + CCFC(IPRFC,ICTFC)*CCFD(IPRFD,ICTFD)
     &            * DMAT(IAOC+ISPC-1,IAOD+ISPD-1)
 120       CONTINUE
 150      CONTINUE

 200     CONTINUE
 300    CONTINUE

        CALL DGEMM('T','N', KCKTC, KHKTD, KHKTC, 1.0D0,
     &              CSQC, KHKTC, DMATSP(1,1,IPRFCD), KHKTC, 
     &              0.0D0, DTMP, KCKTC  )      

        CALL DGEMM('N','N', KCKTC, KCKTD, KHKTD, 1.0D0,
     &              DTMP, KCKTC, CSQD, KHKTD, 
     &              0.0D0, DMATCP(1,1,IPRFCD), KCKTC  )  

        IPRFCD = IPRFCD + 1 

 400   CONTINUE
 500  CONTINUE


      IF ( IPRINT .GE. 20 ) THEN   
        CALL HEADER('DMAT (AO)',-1)
        CALL OUTPUT(DMAT,1,NBASE,1,NBASE,NBASE,NBASE,1,LUPRI)

        CALL PRINTDMAT(1,NBASE,DMAT)

        CALL HEADER('Contraction Matrices',-1)

        WRITE (LUPRI,1000) 'NPRFC, NCTFC',NPRFC,NCTFC,
     &                     'NPRFD, NCTFD',NPRFD,NCTFD

        CALL HEADER('CCFC',5)
        CALL OUTPUT(CCFC,1,NPRFC,1,NCTFC,NPRFC,NCTFC,1,LUPRI)

        CALL HEADER('CCFD',5)
        CALL OUTPUT(CCFD,1,NPRFD,1,NCTFD,NPRFD,NCTFD,1,LUPRI)

        CALL HEADER('DMATSP - Spherical primitives',-1)
          

        WRITE (LUPRI,1100) 
     &       'KHKTC, KHKTD, KHKTCD', KHKTC, KHKTD, KHKTC*KHKTD,
     &       'KPRFC, KPRFD, NPRFCD', NPRFC, NPRFD, NPRFC*NPRFD

        CALL OUTPUT(DMATSP,1,KHKTC*KHKTD,1,NPRFC*NPRFD,
     &             KHKTC*KHKTD, NPRFC*NPRFD,1,LUPRI)
  
        CALL HEADER('Cartesian-spherical Matrices',-1)

        WRITE (LUPRI,1000) 'KHKTC, KCKTC', KHKTC, KCKTC,
     &                     'KHKTD, KCKTD', KHKTD, KCKTD

        CALL HEADER('CSQC',5)
        CALL OUTPUT(CSQC,1,KHKTC,1,KCKTC,KHKTC,KCKTC,1,LUPRI)

        CALL HEADER('CSQD',5)
        CALL OUTPUT(CSQD,1,KHKTD,1,KCKTD,KHKTD,KCKTD,1,LUPRI)

        CALL HEADER('DMATCP - Cartesian primitives',-1)
          
        WRITE (LUPRI,1100) 
     &       'KCKTC, KCKTD, KCKTCD',KCKTC, KCKTD, KCKTC*KCKTD,
     &       'KPRFC, KPRFD, NPRFCD',NPRFC, NPRFD, NPRFC*NPRFD

        CALL OUTPUT(DMATCP,1,KCKTC*KCKTD,1,NPRFC*NPRFD,
     &             KCKTC*KCKTD, NPRFC*NPRFD,1,LUPRI)
    
        CALL HEADER('End of GETDMATCP',-1)

      END IF 

 1000 FORMAT (5X,'Dimensions : ',A,2I10,/,18X,A,2I10,/)
 1100 FORMAT (5X,'Dimensions : ',A,3I10,/,18X,A,3I10,/)

      RETURN
      END 

      SUBROUTINE SCREENFTUV(RODQ,RFTUV,IFTUV,FTUVTMP1,FTUVTMP2,
     &                      ITFTUV,WEIGHT,KODBC2,NDMT,NSAVED,
     &                      EXPCD,IC,ID,QMIN,QMAX,SCR,CCFC,CCFD,
     &                      DMAT,INDXBT,WORK,LWRK) 
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
#include "aobtch.h"
#include "odclss.h"
#include "clsfmm.h"
#include "erithr.h"
#include "pi.h"
#include "ericom.h"

      DIMENSION RODQ(4,*),RFTUV(7,*),IFTUV(2,*),FTUVTMP1(NTUV,NDMT,*),
     &          FTUVTMP2(NTUV,NDMT,NPRFC*NPRFD),
     &          CCFC(NPRFC,NCTFD),CCFD(NPRFC,NCTFD),DMAT(NBASE,NBASE),
     &          EXPCD(NPRFC*NPRFD,3),WEIGHT(NPRFC*NPRFD),WORK(*),
     &          INDXBT(MXSHEL*MXCONT,0:7)
      INTEGER TUV           

      FACCLS = ERFCIV(THRCLS) 

      FACMAX = 0.0D0

      ITFTUV0 = ITFTUV

      IREPC = 0 
      IREPD = 0      

C
C   ( screen the same FTUV's for all the density matrices )
C
C
C       EXPQI=EXPCD(IPRFCD,3)
C       EXPQ =1.0D0/EXPQI
C       
C       WRITE(LUPRI,*) 'EXPQ,QMIN,QMAX)',EXPQ,QMIN,QMAX
C

       EXPQMIN=0.0D0

       IPRFCD=0
       DO IPRFC=1, NPRFC
        DO IPRFD=1, NPRFD

         IPRFCD = IPRFCD + 1

         KEXPC = KEXPBT(IC) - 1 + IPRFC
         KEXPD = KEXPBT(ID) - 1 + IPRFD

         EXPC  = EXPBT(KEXPC)
         EXPD  = EXPBT(KEXPD)

         EXPQ  = EXPC + EXPD

         EXPQMIN = MIN(EXPQMIN,EXPQ)

        END DO
       END DO

C
C   We create different classes for OD with minimum exponent
C   in certain ranges - that should be important for 
C   screening for nonclassical integrals !!!
C

      IF (( EXPQMIN .GE. QMIN ) .AND. ( EXPQMIN .LT. QMAX )) THEN


C
C  We put screening consistent with OD2 distribution to 
C select later the same nonclassical AO's
C

       CX = CORXBT(IC)
       CY = CORYBT(IC)
       CZ = CORZBT(IC)

       DX = CORXBT(ID)
       DY = CORYBT(ID)
       DZ = CORZBT(ID)

       RCD2 = (CX-DX)**2 + (CY-DY)**2 + (CZ-DZ)**2

C       WRITE(LUPRI,'(A,5I3,6F9.4)') 'ScreenFTUV:', KODBC2,IC,ID,
C     &   NPRFC,NPRFD,CX,CY,CZ,DX,DY,DZ
C
C      WRITE(LUPRI,'(A,9I5)') 'ScreenD:', KODBC2,IC,ID,
C     &   NCTFC,NCTFD,IREPC,IREPD,NBASE,NBASE


   
       DSUM = 0.0D0
     
       DO ICTFC = 1, NCTFC
        DO ICTFD = 1, NCTFD

         IAOC = INDXBT(KNDXBT(IC) - 1 + ICTFC,IREPC)
         IAOD = INDXBT(KNDXBT(ID) - 1 + ICTFD,IREPD)
C         WRITE(LUPRI,'(A,4I5)') 'ScreenD1', IAOC,IAOD,KHKTC,KHKTD
C         WRITE(LUPRI,'(A,3I5)') 'ScreenD2', KNDXBT(IC),KNDXBT(ID)
C         WRITE(LUPRI,'(A,4I5)') 'ScreenD2', KNDXBT(IC)-1+ICTFC,
C     &                INDXBT(KNDXBT(IC) - 1 + ICTFC,IREPC),
C     &                KNDXBT(ID)-1+ICTFD,
C     &                INDXBT(KNDXBT(ID) - 1 + ICTFD,IREPD)
        
         DSUMSP = 0.0D0
         DO ISPC=1, KHKTC
          DO ISPD=1, KHKTD
           DSUMSP = DSUMSP + ABS(DMAT(IAOC+ISPC-1,IAOD+ISPD-1))
          END DO
         END DO

         CCSUM  = 0.0D0
         DO IPRFC=1, NPRFC
          DO IPRFD=1, NPRFD
           CCSUM=CCSUM+ ABS(CCFC(IPRFC,ICTFC)*CCFD(IPRFD,ICTFD))
          END DO
         END DO

         DSUM=DSUM+DSUMSP*CCSUM

        END DO
       END DO

C      WRITE(LUPRI,'(A,9I5,F12.8)') 'ScreenD:', KODBC2,IC,ID,
C     &   NCTFC,NCTFD,IREPC,IREPD,KHKTC,KHKTD,DSUM

C       DSUM=1.0D0

       SM=0.0D0
       CM=0.0D0
       DM=0.0D0
       QIM=0.0D0

       IPRFCD=0
       DO 200 IPRFC=1, NPRFC
        DO 100 IPRFD=1, NPRFD

         IPRFCD = IPRFCD + 1

         KEXPC = KEXPBT(IC) - 1 + IPRFC
         KEXPD = KEXPBT(ID) - 1 + IPRFD

         KCCFC = KCCFBT(IC) - 1 + IPRFC
         KCCFD = KCCFBT(ID) - 1 + IPRFD

         EXPC  = EXPBT(KEXPC)
         EXPD  = EXPBT(KEXPD)

         EXPQ  = EXPC + EXPD
         EXPQI = 1.0D0/EXPQ
         EXPMU = EXPC*EXPD/EXPQ

         FMAX = 0.0D0
         DO IDMT=1,NDMT
          DO TUV=1,NTUV
           FMAX=MAX(FMAX,ABS(FTUVTMP2(TUV,IDMT,IPRFCD))*EXPQI)
          END DO
         END DO

         FMAX = FMAX * R2PI52

         FAC=0.0D0

         DO ICTFC = 1, NCTFC
          DO ICTFD = 1, NCTFD
           FAC = FAC+ABS(CCFC(IPRFC,ICTFC)*CCFD(IPRFD,ICTFD))
          END DO
         END DO      

         FAC = R2PI52*EXP(-EXPMU*RCD2)*EXPQI * FAC

         IF ( ABS(FAC) .GT. THRSH ) THEN 

          QIM = MAX(QIM,EXPQI)
          ITFTUV = ITFTUV+1                 

          WT=0.0D0
          DO ICTFC = 1, NCTFC
           DO ICTFD = 1, NCTFD
            WT = WT + ABS(CCFC(IPRFC,ICTFC)*CCFD(IPRFD,ICTFD))
           END DO
          END DO      

          WT = WT * EXP(-EXPMU*RCD2)
          SM = SM + WT
          CM = CM + WT * EXPC * EXPQI
          DM = DM + WT * EXPD * EXPQI
C          WRITE(LUPRI,'(4F12.8)') WT,SM,CM,DM
          CALL DCOPY(NDMT*NTUV,FTUVTMP2(1,1,IPRFCD),1
     &                      ,FTUVTMP1(1,1,ITFTUV),1)

          RFTUV(1,ITFTUV) = EXPQ
          RFTUV(2,ITFTUV) = EXPQI

          RFTUV(3,ITFTUV) = (EXPC*CORXBT(IC)+EXPD*CORXBT(ID))*EXPQI
          RFTUV(4,ITFTUV) = (EXPC*CORYBT(IC)+EXPD*CORYBT(ID))*EXPQI
          RFTUV(5,ITFTUV) = (EXPC*CORZBT(IC)+EXPD*CORZBT(ID))*EXPQI 
          RFTUV(6,ITFTUV) = FMAX
          RFTUV(7,ITFTUV) = SCR


          FACMAX = MAX( FACMAX, FAC ) 

         END IF

 100    CONTINUE
 200   CONTINUE

C       WRITE(6,*) 'Screenfac:', FACMAX, SCR

       CM = CM/SM
       DM = DM/SM

       RODQ(1,KODBC2) = CM * CX + DM * DX
       RODQ(2,KODBC2) = CM * CY + DM * DY
       RODQ(3,KODBC2) = CM * CZ + DM * DZ
       RODQ(4,KODBC2) = FACCLS * SQRT(QIM)
 
C       WRITE(LUPRI,'((3X,A,I5,3F12.7))') 
C     &        'CM,DM,SM',KODBC2,CM,DM,SM
C
C       WRITE(LUPRI,'((3X,A,6F7.3))') 
C     &        'C,D', CX,CY,CZ,DX,DY,DZ
C       WRITE(LUPRI,'((3X,A,4F7.3))') 
C     &        'RODQ', RODQ(1,KODBC2),RODQ(2,KODBC2),
C     &                RODQ(3,KODBC2),RODQ(4,KODBC2)


       END IF

       NSAVED = ITFTUV - ITFTUV0
 
       DO I=ITFTUV0+1,ITFTUV
        IFTUV(1,I) = KODBC2
        IFTUV(2,I) = NSAVED

C        WRITE(LUPRI,'(A,6I6)') 'SCR-IFTUV:',
C     &     I, KODBC2, NSAVED, IFTUV(1,I), IFTUV(2,I), NPRFCD
       END DO

      RETURN
      END

      SUBROUTINE REORDERFTUV(IFTUVCL,KIFTUVCL,
     &                       KRFTUVCL,KIFTUV,KRFTUV,KFTUV,KFTUVTMP,
     &                       NDMT,KFREE,LFREE,WORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aobtch.h"
#include "odclss.h"
#include "iratdef.h"
C
C aspg, 2006-06-27 making ftnchek happy
C
C#include "fck3odb.h"
#include "ericom.h"
C
      DIMENSION IFTUVCL(10,*),WORK(*)
      INTEGER SIZEIFTUVCL,SIZERFTUVCL,SIZEIFTUV,SIZERFTUV,SIZEFTUV

C      CALL PRIIFTUV(WORK(KIFTUV))

      SIZEIFTUVCL= 10 * NFTUVCL/IRAT + 1
      SIZERFTUVCL=  4 * NFTUVCL

      SIZEIFTUV  = 0
      SIZERFTUV  = 0
      SIZEFTUV   = 0
 
      NFTUV = 0 
 
      DO ICL=1,NFTUVCL
       NTFTUV = IFTUVCL(1,ICL)
       NHKTF  = IFTUVCL(5,ICL)

       NFTUV  = NFTUV + NTFTUV
       NTUV   = NHKTF**3 
       
       SIZEIFTUV   = SIZEIFTUV + 2 * NTFTUV / IRAT + 1
       SIZERFTUV   = SIZERFTUV + 7 * NTFTUV 
       SIZEFTUV    = SIZEFTUV  + NTFTUV * NTUV * NDMT
      END DO      

      KTMP = KFTUVTMP + SIZEFTUV

      CALL DCOPY(4*NFTUVCL,WORK(KRFTUVCL),1,WORK(KTMP),1)
      KRFTUVCL=KIFTUVCL+SIZEIFTUVCL
      CALL DCOPY(4*NFTUVCL,WORK(KTMP),1,WORK(KRFTUVCL),1)

      CALL ICOPY(2*NFTUV,WORK(KIFTUV),1,WORK(KTMP),1)
      KIFTUV=KRFTUVCL+SIZERFTUVCL
      CALL ICOPY(2*NFTUV,WORK(KTMP),1,WORK(KIFTUV),1)

      CALL DCOPY(7*NFTUV,WORK(KRFTUV),1,WORK(KTMP),1)
      KRFTUV=KIFTUV+SIZEIFTUV
      CALL DCOPY(7*NFTUV,WORK(KTMP),1,WORK(KRFTUV),1)

      KFTUV    = KRFTUV + SIZERFTUV

C
C  first copy to the safe location 
C
      IF ( LFREE .LT. SIZEFTUV ) THEN
        CALL QUIT('Insufficient memory in GETFTUV,reorder')      
      END IF

      CALL DCOPY( SIZEFTUV, WORK(KFTUVTMP),1,
     &                     WORK(KFTUVTMP+SIZEFTUV),1 )

C
C   than reorder
C
      DO ICL=1,NFTUVCL
       NTFTUV = IFTUVCL(1,ICL)
       KSFTUV = IFTUVCL(3,ICL)
       NHKTF  = IFTUVCL(5,ICL)

       NTUV   = NHKTF**3
       
       CALL REORDER_FTUV2(WORK(KFTUVTMP+SIZEFTUV+KSFTUV-1),
     &                   WORK(KFTUV+KSFTUV-1),NTFTUV,NTUV,NDMT)
       
      END DO  

      LFREE = LFREE - (KFTUV+SIZEFTUV-KFREE)
      KFREE = KFTUV + SIZEFTUV

      RETURN
      END

      SUBROUTINE REORDER_FTUV2(FTUVTMP,FTUV,NTFTUV,NTUV,NDMT)
#include "implicit.h"
#include "priunit.h"

      DIMENSION FTUVTMP(NTUV*NDMT,*),FTUV(NTFTUV,NTUV*NDMT)

      N=NDMT*NTUV

      DO I=1,NTFTUV
       DO J=1,N
        FTUV(I,J)=FTUVTMP(J,I)
       END DO
      END DO

      RETURN
      END

      SUBROUTINE SUMMARY_FTUV(IFTUVCL,RFTUVCL,IFTUV,RFTUV,FTUV,NDMT,
     &                        KIFTUVCL,KRFTUVCL,KIFTUV,KRFTUV,KFTUV,
     &                        IODCL2,KFREE,LFREE,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "ericom.h"
#include "odclss.h"
#include "clsfmm.h"
       DIMENSION IODCL2(NODCL2,NITCL),
     &           IFTUVCL(10,*),RFTUVCL(4,*),
     &           IFTUV(*),RFTUV(*),FTUV(*)

C      WRITE(LUPRI,'(A,6I6)') 'IFTUV-beg.: ', (IFTUV(I),I=1,6)

      FACCLS = ERFCIV(THRCLS) 

      IF ( IPRINT .GE. 2 ) THEN
       CALL HEADER('FTUV summary:',-1)
       WRITE(LUPRI,'(3X,11(1X,A))') 'Class','NTFTUV', 
     &        'KSFTUV','JMAXF','NCNTF','ICNTC',' ICNTD', 
     &        '  RCD ',' max extent','   min F',' Cauchy'
      END IF

       NFTUV=0

       DO ICL=1,NFTUVCL
 
        NTFTUV=IFTUVCL( 1, ICL )
        KRFTUV_REL=IFTUVCL( 2, ICL )
        KSFTUV=IFTUVCL( 3, ICL )
        JMAXF =IFTUVCL( 4, ICL )
        NHKTF =IFTUVCL( 5, ICL )
        NCNTF =IFTUVCL( 6, ICL )
        ICNTC =IFTUVCL( 7, ICL )
        ICNTD =IFTUVCL( 8, ICL )
        KIFTUV_REL=IFTUVCL( 9, ICL )

        NTUV   = NHKTF**3
 
        NFTUV = NFTUV + NTFTUV 

        EXPQIMAX  = 0.0D0
        FTUV_MAX  = 0.0D0
        CAUCH_MAX = 0.0D0
        DO I=1,NTFTUV
         EXPQIMAX  = MAX( EXPQIMAX,  RFTUV( KRFTUV_REL + 7*(I-1) + 1 ) )
         FTUV_MAX  = MAX( FTUV_MAX,  RFTUV( KRFTUV_REL + 7*(I-1) + 5 ) )
         CAUCH_MAX = MAX( CAUCH_MAX, RFTUV( KRFTUV_REL + 7*(I-1) + 6 ) )
        END DO

        EXTQMAX = FACCLS*SQRT(EXPQIMAX)

        RFTUVCL(2,ICL)=EXTQMAX
        RFTUVCL(3,ICL)=FTUV_MAX
        RFTUVCL(4,ICL)=CAUCH_MAX

        IF ( IPRINT .GE. 2 ) THEN
         WRITE(LUPRI,'(3X,2I5,I8,4I6,2X,2F8.4,2F12.6)') ICL,
     &       NTFTUV, KSFTUV, JMAXF, NCNTF, ICNTC, ICNTD,
     &       RFTUVCL(1,ICL),RFTUVCL(2,ICL),RFTUVCL(3,ICL),RFTUVCL(4,ICL)
         IF ( IPRINT .GE. 5 ) THEN
C
C aspg, 2006-06-27 making ftnchek happy
C
C          CALL PRIFTUV2(FTUV(KSFTUV),ICL,NTFTUV,JMAXF,NDMT)
          CALL PRIFTUV2(FTUV(KSFTUV),ICL,NTFTUV,JMAXF)
C
          DO I=1,NTFTUV
           WRITE(LUPRI,'(3X,I5,7F12.7,I4)')
     &        I,(RFTUV(KRFTUV_REL+7*(I-1)+J-1),J=1,7),
     &        IFTUV(KIFTUV_REL+2*(I-1)+1)
          END DO      
         END IF
        END IF
       END DO 

       NODBC2P=0
       DO KODCL2=1,NODCL2
        NODTCD = IODCL2(KODCL2,1)
        NPRFC  = IODCL2(KODCL2,15)
        NPRFD  = IODCL2(KODCL2,16)        
        NPRFCD = NPRFC*NPRFD
        NODBC2P = NODBC2P + NPRFCD * NODTCD
       END DO

       IF ( IPRINT .GE. 2 ) THEN
        CALL HEADER('Screening of FTUVs',-1)

        WRITE(LUPRI,'(3X,A,I6,3X,A,I6,F12.6,A)')
     &         'NFTUV(total):',NFTUV,
     &         '/  NOD2primitivs  :',NODBC2P,
     &         100.0D0*NFTUV/NODBC2P, ' % '

        CALL HEADER('Memory requirements for storing FTUV',-1)

        WRITE(LUPRI,'(6(3X,A,I12,3X,A,I12,/))')
     &      'KIFTUVCL :',KIFTUVCL,'LIFTUVCL',KRFTUVCL - KIFTUVCL,
     &      'KRFTUVCL :',KRFTUVCL,'LRFTUVCL',KIFTUV   - KRFTUVCL,
     &      'KIFTUV   :',KRFTUV,  'LRFTUV  ',KRFTUV   - KIFTUV,
     &      'KRFTUV   :',KRFTUV,  'LRFTUV  ',KFTUV    - KRFTUV,
     &      'KFTUV    :',KFTUV,   'LFTUV   ',KFREE    - KFTUV,
     &      'KFREE    :',KFREE,   'LFREE   ',LFREE 

       END IF

      RETURN
      END

      SUBROUTINE CHECK_RNAI(IFTUVCL,RFTUV,FTUV,DMAT,NDMT,
     &                      WORK,LWORK,IPRINT)
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

      DIMENSION IFTUVCL(10,*),RFTUV(*),FTUV(*),DMAT(*),WORK(*)

      RNA = 0.0D0

      DO ICL=1,NFTUVCL
       NTFTUV = IFTUVCL(1,ICL)
       KRFTUV = IFTUVCL(2,ICL)
       KSFTUV = IFTUVCL(3,ICL)
       JMAXF  = IFTUVCL(4,ICL)
       NHKTF  = IFTUVCL(5,ICL)

       NTUV   = NHKTF**3

       KRNAI  = 1 
       KAHGTF = KRNAI + NTFTUV*NTUV

       CALL DZERO(WORK(KRNAI),NTFTUV*NTUV*NDMT)

       CALL GETRNAI(WORK(KRNAI),WORK(KAHGTF),NTFTUV,JMAXF, 
     &              RFTUV(KRFTUV),7,IPRINT)

       IF ( IPRINT .GE. 20 ) THEN
        CALL PRIRNAI(WORK(KRNAI),ICL,NTFTUV,JMAXF)
        CALL PRIFTUV2(FTUV(KSFTUV),ICL,NTFTUV,JMAXF)
       END IF

       RNA = RNA + ADDRNAI(WORK(KRNAI),
     &               FTUV(KSFTUV),NTFTUV, JMAXF )

      END DO

      CALL HEADER('Nuclear attraction using FTUV',-1)
      WRITE (LUPRI,'(F16.10)') RNA

      CALL REFERENCE_NA(REF_NA,DMAT,NDMT,WORK,LWORK, IPRINT)

    
      RETURN
      END

      SUBROUTINE GETRNAI( RNAI, AHGTF, NPP12, JMAX12,
     &                    RFTUV,NRFTUV,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "ericom.h"
#include "cbieri.h"
#include "aobtch.h"
#include "symmet.h"
#include "odbtch.h"
#include "odclss.h"
#include "nuclei.h"
C
      DIMENSION RNAI(NPP12,0:JMAX12,0:JMAX12,0:JMAX12),
     &          RFTUV(NRFTUV,*),AHGTF((1+1)*(JMAX12+1)**3)

      ISTEPU = JMAX12 + 1
      ISTEPV = ISTEPU*ISTEPU
      NAHGTF = ISTEPU*ISTEPV

      DO 200 IPP12 = 1, NPP12

        IADR = 1

        EXPP = RFTUV(1,IPP12)
        EXPPI= RFTUV(2,IPP12)

        COORPX = RFTUV(3,IPP12)
        COORPY = RFTUV(4,IPP12)
        COORPZ = RFTUV(5,IPP12)

        CALL DZERO(AHGTF,(1+1)*(JMAX12+1)**3) 

        DO 100 IATOMC = 1, NATOMS

         FACTOR = - CHARGE(IATOMC)

         DIFCPX = CORD(1,IATOMC) - COORPX
         DIFCPY = CORD(2,IATOMC) - COORPY
         DIFCPZ = CORD(3,IATOMC) - COORPZ

         CALL HERNAI(AHGTF,JMAX12,EXPP,
     &               DIFCPX,DIFCPY,DIFCPZ,FACTOR,IADR,
     &               ISTEPU,ISTEPV,NAHGTF,IPRINT)      

         IA = 0 
         DO K=0,JMAX12
          DO J=0,JMAX12
           DO I=0,JMAX12
            IA = IA + 1
            RNAI(IPP12,I,J,K) = RNAI(IPP12,I,J,K) + AHGTF(IA)
           END DO
          END DO 
         END DO
 
 100    CONTINUE
 200  CONTINUE

      RETURN
      END 

      FUNCTION ADDRNAI(RNAI, FTUV, NPP12, JMAX12 )
#include "implicit.h"
C
      INTEGER T,U,V
      DIMENSION RNAI(NPP12,0:JMAX12,0:JMAX12,0:JMAX12),
     &          FTUV(NPP12,0:JMAX12,0:JMAX12,0:JMAX12)

      ADDRNAI=0.0D0

      DO I=1,NPP12
       DO T=0,JMAX12
        DO U=0,JMAX12
         DO V=0,JMAX12

           ADDRNAI = ADDRNAI + RNAI(I,T,U,V)*FTUV(I,T,U,V)
     &                 
         END DO
        END DO
       END DO
      END DO

      RETURN
      END 

      SUBROUTINE PRINTDMAT( NSIM, NORBT, D )   
#include "implicit.h"
#include "priunit.h" 

      DIMENSION D(*)

      DO 200 ISIM=1,NSIM 
        WRITE( LUPRI, * ) 'MATRIX NUMBER : ', ISIM
          DO 100 I=1,NORBT
            WRITE( LUPRI, 400 ) (D((I-1)*NORBT+J),J=1,NORBT)
 100      CONTINUE
 200  CONTINUE

 400  FORMAT(6F16.12)

      RETURN
      END
C
C *** END OF PRINTDMAT
C

      SUBROUTINE REFERENCE_NA(REF_NA, DMAT,NDMT,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h" 
#include "inforb.h"

C
C So far only for 1 density matrix ...
C (NDMT=1)
C
      LOGICAL FOUND1, FOUND2
      DIMENSION DMAT(NBAST,NBAST),WORK(*)

#ifdef PRG_DIRAC

      WRITE(LUPRI,*)
     &  'WARNING: REFERENCE_NA not implemented for Dirac, nothing done.'

#else

C aspg, 2006-06-29, making ftnchek happy
C
      FOUND1 = .FALSE.
      FOUND2 = .FALSE.
C
      KH1AO  = 1
      KKINAO = KH1AO  + NNBAST
      KPOTAO = KKINAO + NNBAST
      KWRK1  = KPOTAO + N2BAST
      LWRK1  = LWORK - KWRK1

      CALL DZERO(WORK,KWRK1-KH1AO)

      CALL RDONEL('ONEHAMIL',FOUND1,WORK(KH1AO), NNBAST)
      CALL RDONEL('KINETINT',FOUND2,WORK(KKINAO),NNBAST)

      REF_NA = 0.0D0
   
      K    = 0  
      DO 200 I=1,NBAST
        DO 100 J=1,I
C           WORK(KPOTAO+(I-1)*NBAST+J-1)=
C     &         WORK(KH1AO+K)-WORK(KKINAO+K)
C           WORK(KPOTAO+(J-1)*NBAST+I-1)=
C     &         WORK(KH1AO+K)-WORK(KKINAO+K)


          IF ( I .EQ. J ) THEN 
            REF_NA = REF_NA +  DMAT(I,J) *
     &        ( WORK(KH1AO+K)-WORK(KKINAO+K) )
          ELSE 
            REF_NA = REF_NA +  
     &        ( DMAT(I,J) + DMAT(J,I) ) *
     &        ( WORK(KH1AO+K)-WORK(KKINAO+K) )
          END IF 
          K=K+1
 100    CONTINUE
 200  CONTINUE

      CALL HEADER('Reference value for nuclear attraction',-1)

       WRITE (LUPRI,'(F30.16)') REF_NA

       IF ( .NOT. ( FOUND1 .AND. FOUND2 )) 
     &   WRITE (LUPRI,*) 'Error : Integrals not found !'

#endif

      RETURN
      END

      SUBROUTINE PRIRNAI(RNAI,KODBC2,NPRFAB,JMAX)
#include "implicit.h"
#include "priunit.h"       
      DIMENSION RNAI(NPRFAB,0:JMAX,0:JMAX,0:JMAX)
      INTEGER T,U,V

      CALL HEADER('RNAI',-1)

C      WRITE(LUPRI,'(/,7X,3(A,I5,A,3X),/)')
C     &            'OD batch number :',KODBC2,',',
C     &            'NPRFAB :', NPRFAB,',',
C     &            'JMAX :', JMAX,' '


      WRITE(LUPRI,'(7X,A)') 
     &            'I    T    U    V     RNAI(I,T,U,V)'

      DO I = 1, NPRFAB

       DO T = 0, JMAX
        DO U = 0, JMAX
         DO V = 0, JMAX
          IF ( ( T+U+V .LE. JMAX ) .AND. 
     &         ( ABS(RNAI(I,T,U,V)) .GT. 1.0D-10 ) ) THEN
           WRITE(LUPRI,'(3X,4I5,F16.8)') I,T,U,V,RNAI(I,T,U,V)
          END IF
         END DO
        END DO
       END DO

      END DO

      RETURN
      END


      SUBROUTINE PRIFTUV(FTUV,KODBC2,NPRFAB,JMAX)
#include "implicit.h"
#include "priunit.h"       
      DIMENSION FTUV(0:JMAX,0:JMAX,0:JMAX,NPRFAB)
      INTEGER T,U,V

      CALL HEADER('FTUV',-1)

      WRITE(LUPRI,'(/,7X,3(A,I5,A,3X),/)')
     &            'OD batch number :',KODBC2,',',
     &            'NPRFAB :', NPRFAB,',',
     &            'JMAX :', JMAX,' '


      WRITE(LUPRI,'(7X,A)') 
     &            'I    T    U    V     FTUV(T,U,V,I)'

      DO I = 1, NPRFAB

       DO T = 0, JMAX
        DO U = 0, JMAX
         DO V = 0, JMAX

          IF ( T+U+V .LE. JMAX ) THEN
            WRITE(LUPRI,'(3X,4I5,F16.8)') I,T,U,V,FTUV(T,U,V,I)
          ELSE
            IF ( FTUV(T,U,V,I) .GT. 0D0 ) THEN 
             WRITE(LUPRI,'(3X,4I5,F16.8,A)') I,T,U,V,FTUV(T,U,V,I),
     &                   ' ERROR - T+U+V > JMAX ! '
            END IF 
          END IF 
         END DO
        END DO
       END DO

      END DO

      RETURN
      END

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


      SUBROUTINE PRIIFTUV(IFTUV)
#include "implicit.h"
#include "priunit.h"

      DIMENSION IFTUV(*)

      WRITE(LUPRI,'(A,6I6)') 'IFTUV-beg.0: ', (IFTUV(I),I=1,6)

      RETURN 
      END


      SUBROUTINE PRIFULLIFTUV(IFTUVCL,IFTUV,NFTUVCL)
#include "implicit.h"
#include "priunit.h"

      DIMENSION IFTUVCL(10,*),IFTUV(*)

      CALL HEADER('PRIFULLIFTUV',-1)

      DO I=1,NFTUVCL
        NTFTUV     = IFTUVCL( 1,I)
        KIFTUV_REL = IFTUVCL( 9,I)
        DO J=1,NTFTUV
          WRITE(LUPRI,'(A,6I6)') 'IFTUV :', 
     &            I, J, (IFTUV(KIFTUV_REL+2*(J-1)+K-1),K=1,2)
        END DO
      END DO


      RETURN 
      END
