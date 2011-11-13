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
C
#if defined (VAR_SINGLE)
#include "single.h"
#endif
C
C  Definition of the FCK3AOB - ao batches structure...
C  

      SUBROUTINE SETAOBTCH(SPEXP,CONTRACTED,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      DIMENSION WORK(LWORK)
      LOGICAL SPEXP, CONTRACTED
      DIMENSION IORB(3), KHKT_REL(3)
#include "shells.h"
#include "symmet.h"
#include "nuclei.h"
#include "fck3aob.h"
#include "fck3mem.h"
#include "fck3rel.h"
#include "ibtfun.h"
      CALL QENTER('SETAOBTCH')
C
C     Calculate the memory requirements of the fck3aob common block arrays
C ( Note that EXPBTCH and CCFBTCH can be considerably reduced in size if
C   we utilize the fact that we usually have (and especially for large 
C   systems) many atoms sharing the same basis set )
C

      MEMAOB = MXPRIM        + MXPRIM*MXCONT
     &       + MXSHEL*MXCONT + 4*MXSHEL      + 2*MXAOBTCH + MXAOBTCH
     &       + MXAOBTCH      + MXAOBTCH      + MXAOBTCH   + 2*MXAOBTCH
     &       + MXAOBTCH      + MXAOBTCH      + MXAOBTCH   + 2*MXAOBTCH

C corresponding to arrays:
C              EXPBTCH       CCFBTCH       
C              INDXBTCH       IAOCLS          KNDXBTCH    NHKTBTCH    
C              NPRFBTCH      NCTFBTCH         KEXPBTCH    KCCFBTCH
C              NTYPBTCH      NCNTBTCH

C
C     Number of basis functions + max. NHKT
C     =====================================
C
      MAXAONHK = 1
      NBASE  = 0

      DO K = 0, MAXREP
cSimen Should change MAXAONHK to normal and auxilliary part, might 
cSimen also want NBASE and NBASEAUX
      DO I = 1, KMAX
         IF (IBTAND(K,ISTBAO(I)).EQ.0) NBASE = NBASE + KHKT(I)
         MAXAONHK = MAX( MAXAONHK,  NHKT(I) )  
      END DO
      DO I = KMAX+1, KMAXTOT
         MAXAONHK = MAX( MAXAONHK,  NHKT(I) )  
      END DO
      END DO
      MAXAONHK = MAXAONHK + NDER

      MAXAOKCK =     MAXAONHK * ( MAXAONHK + 1 ) / 2 
      MAXAOKHK = 2 * MAXAONHK - 1

C
C
C     AO batches
C     ==========
C

      EXPAOMIN = 1.0D10
      CCFAOMAX = 0.0D0

      MAXAOPRF = 1 
      MAXAOCTF = 1 
C     Set the offsets for the large and small component
      IORB(1) = 1
      IORB(2) = NLARGE+1
      IORB(3) = NLARGE+1
      IAOIND = 1
      NB = -1
      NC = -1
      ND = -1

      DO ISHELL = 1, KMAXTOT

c       The small component shells (LCLASS.eq.2) are done implicitly, using
c       kinetic balance on the large component. This generates the (modified)
c       solid harmonics from which we can construct all that is needed.
        KHKT_REL(1) = KHK_SC(NHKT(ISHELL),0)
        KHKT_REL(2) = KHK_SC(NHKT(ISHELL),1)
        KHKT_REL(3) = KHK_SC(NHKT(ISHELL),-1)
C       Do some error checking (LC should be spherical)
        IF (KHKT_REL(1).NE.KHKT(ISHELL)) 
     &  CALL QUIT ("Wrong choice of angular momentum options for FCK3")

        IF (.NOT.LCLASS(ISHELL).EQ.2) THEN
C         The modified solid harmonics follow after the regular ones
          IORB(3) = IORB(2) + KHKT_REL(2)
          IF (NB.NE.NBCH(ISHELL) .OR. NC.NE.NCENT(ISHELL) .OR.
     &        ND.NE.LCLASS(ISHELL) ) THEN
           NB = NBCH(ISHELL)
           NC = NCENT(ISHELL)
           ND = LCLASS(ISHELL)
           
           CALL AOBTCH1(IORB,IAOIND,NB,ISHELL,WORK,LWORK,IPRINT)
          END IF
C         Increase the orbital counter for the large components
          IORB(1) = IORB(1) + KHKT_REL(1)
C         Increase the orbital counter for the kinetically balanced small components
          IORB(2) = IORB(2) + KHKT_REL(2) + KHKT_REL(3)
       ENDIF 
      END DO


C     Check s-p exponents
C     ===================
C

      IF ( SPEXP ) CALL SETSPEXP(IPRINT)
C
C     AO classes
C     ==========
C
      CALL SETAOCLS(IPRINT)
C
C
C     If requested - internally decontract AO's  
C     =========================================
C (Only for the internal use in FCK3 and FMM. Doesn't 
C  work with general contractions )
C
      IF ( .NOT. CONTRACTED ) CALL DECONTRACTAO(WORK,LWORK,IPRINT)
C
C     Print section
C     =============
C

      IF (IPRINT .GT. 2) THEN
       CALL HEADER('AO classes',-1)

       WRITE(LUPRI,'(3X,A,I5)') 'NAOBCLS:', NAOBCLS
       DO ICLS = 1, NAOBCLS
        KSTAOB = IAOCLS(2,ICLS)
        NAOB   = IAOCLS(3,ICLS)
        WRITE(LUPRI,'(3X,A,5I5)') ' Class number :',ICLS,IAOCLS(1,ICLS),
     &           IAOCLS(2,ICLS),IAOCLS(3,ICLS)
        DO IAOB = KSTAOB, ( KSTAOB + NAOB - 1 )
          WRITE(LUPRI,'(5X,A,I5)') 'IAOB', IAOB 
        END DO
       END DO

       CALL HEADER('Common block FCK3AOBTCH in SETAOBTCH',-1)
       WRITE (LUPRI,'(2(3X,A,I8),3X,A,5I5)') 
     &             ' Number of classes NAOBCLS:',NAOBCLS,
     &             ' Number of batches NAOBTCH:',NAOBTCH,
     &             ' MAXAONHK, MAXAOKCK, MAXAOKHK, MAXAOPRF, MAXAOCTF,',
     &               MAXAONHK, MAXAOKCK, MAXAOKHK, MAXAOPRF, MAXAOCTF

       CALL HEADER(
     &' # IKB  NCNT  NHKT  NPRF  NCTF  NTYP KNDX1 KNDX2 INDAO1 INDAO2'
     & ,1)
       DO IKB = 1, 3
       DO I = 1, NAOBTCH
        WRITE (LUPRI,'(11I6)') I,IKB, 
     &   NCNTBTCH(I), NHKTBTCH(I), NPRFBTCH(I),
     &   NCTFBTCH(I),NTYPBTCH(I),KNDXBTCH(1,I,IKB),KNDXBTCH(2,I,IKB),
     &   INDXBTCH(KNDXBTCH(1,I,IKB),IKB),INDXBTCH(KNDXBTCH(2,I,IKB),IKB)
       END DO
       END DO
      END IF

      IF (IPRINT .GT. 5) THEN

       CALL HEADER('Exponents and contraction coefficients',-1)
       DO I = 1, NAOBTCH
        WRITE (LUPRI,'(/,1X,A,I2,A,2I3,/)')
     &   ' Exp. and cont. coef. for AO batch',I,
     &   ' of dimensions',NPRFBTCH(I),NCTFBTCH(I)
        KEXP = KEXPBTCH(I)
        KCCF = KCCFBTCH(1,I)
        NPRF = NPRFBTCH(I)
        NCTF = NCTFBTCH(I)
        NTYP = NTYPBTCH(I)

        DO J = 1, NPRF
         WRITE (LUPRI,'(2X,E8.3,7(2X,E8.2))') EXPBTCH(KEXP+J-1),
     &     ( CCFBTCH(KCCF+J-1+NPRF*(K-1)), K=1, NCTF)
        END DO 
 
        IF ( NTYP .EQ. 2 ) THEN 
         KCCF = KCCFBTCH(2,I)
         WRITE (LUPRI,'(/,2X,A,/)') 'p-function with same exponents:'
         DO J = 1, NPRF
          WRITE (LUPRI,'(2X,E8.3,7(2X,E8.2))') EXPBTCH(KEXP+J-1),
     &     ( CCFBTCH(KCCF+J-1+NPRF*(K-1)), K=1, NCTF)           
         END DO
        END IF
       END DO

C       CALL HEADER('Coordinates for AO batches',-1)
C       DO I = 1, NAOBTCH
C        WRITE (LUPRI,'(2X,I5,5X,3F12.6)') I, CORD(1,NCNTBTCH(I)),
C     &     CORD(2,NCNTBTCH(I)), CORD(3,NCNTBTCH(I))
C       END DO

       CALL HEADER('Orbital indices for AO batches',-1)
       DO I = 1, NAOBTCH
        WRITE (LUPRI,'(2X,I5,3X,(40I3))') I,
     &        (INDXBTCH(KNDXBTCH(1,I,1)+J-1,1),J=1,NCTFBTCH(I))
       END DO
      END IF
C
      CALL QEXIT('SETAOBTCH')
      RETURN
      END
C
C ---------------------------------------------------------------
C

      SUBROUTINE AOBTCH1(IORB,IAOIND,NB,ISHELL,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      DIMENSION IORB(3)
      DIMENSION WORK(LWORK)
#include "shells.h"
      CALL QENTER('AOBTCH1')

C
      NUCA  = NUCO(ISHELL)
      NRCA  = NRCO(ISHELL)
      JSTRA = JSTRT(ISHELL)
C
      KMAT   = 1
      KDONE  = KMAT   + NUCA*NRCA
      KDIMP  = KDONE  + NUCA
      KDIMC  = KDIMP  + NUCA
      KPNTC  = KDIMC  + NUCA
      KPNTP  = KPNTC  + NUCA*NRCA
      KLAST  = KPNTP  + NUCA*NUCA

      IF (KLAST .GT. LWORK) CALL STOPIT('AOBTCH1',' ',KLAST,LWORK)
C
      CALL AOBTCH2(IAOIND,NB,NMAT,WORK(KMAT),NUCA,NRCA,
     &             JSTRA,WORK(KDONE),
     &             WORK(KDIMP),WORK(KDIMC),WORK(KPNTC),WORK(KPNTP),
     &             ISHELL,IORB,IPRINT)
C
      CALL QEXIT('AOBTCH1')
      RETURN
      END

C
C -----------------------------------------------------
C

      SUBROUTINE AOBTCH2(IAOIND,NB,NMAT,MAT,NUCA,NRCA,JSTRA,IDONE,NDIMP,
     &                   NDIMC,IPNTC,IPNTP,ISHELL,IORB,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
      PARAMETER (THR = 1.0D-17)
C
      DIMENSION MAT(NUCA,NRCA), IDONE(NUCA), NDIMP(NUCA), NDIMC(NUCA),
     &          IPNTC(NUCA,NRCA),IPNTP(NUCA,NUCA)
      DIMENSION IORB(3),IADRO(3),KHKT_REL(3)
C
#include "symmet.h"
#include "fck3aob.h"
#include "shells.h"
#include "primit.h"
#include "nuclei.h"
C
      SAVE IADRP, IADRC, IADRO, NBTCH

C     First get the degeneracy of the (modified) solid harmonics
C     NB: further down IKB will denote 1:LC, 2:SC1, 3:SC2
C     with SC1 and SC2 the functions generated by 
C     kinetic balance upwards or downwards respectively
C
      KHKT_REL(1) = KHK_SC(NHKT(ISHELL),0)
      KHKT_REL(2) = KHK_SC(NHKT(ISHELL),1)
      KHKT_REL(3) = KHK_SC(NHKT(ISHELL),-1)
C
C     Construct matrix MAT containing 1 and 0 only, then fill
C     out all elements within each contraction matrix
C
      CALL IZERO(MAT,NUCA*NRCA)
C
      DO I = 1, NUCA
      DO J = 1, NRCA
         IF (ABS(PRICCF(JSTRA+I,J)).GT.THR) MAT(I,J) = 1
      END DO
      END DO
C
      DO I = 1, NRCA
      DO J = 1, NRCA
         IF (INPRD(MAT(1,I),MAT(1,J),NUCA) .GT. 0) THEN
            DO K = 1, NUCA
               MAT(K,I) = MAX(MAT(K,I),MAT(K,J))
               MAT(K,J) = MAX(MAT(K,I),MAT(K,J))
            END DO
         END IF
      END DO
      END DO
C
C     NMAT, NDIMP(), NDIMC(), IPNTP(), IPNTC()
C
      NMAT = 0
      CALL IZERO(IDONE,NUCA)
      DO I = 1, NRCA
         IF (IDONE(I).EQ.0) THEN
C
            NMAT = NMAT + 1
            DO J = I, NRCA
               IDONE(J) =  INPRD(MAT(1,I),MAT(1,J),NUCA)
            END DO
C
C           First non-zero primitive
C
            KCONT = I
            DO J = 1, NUCA
               IF (MAT(J,KCONT) .EQ. 1) THEN
                  KPRIM = J
                  GO TO 100
               END IF
            END DO
  100       CONTINUE
C
C           IPNTP and NDIMP
C
            IPRI = 0
            DO J = 1, NUCA
            IF (MAT(J,KCONT) .EQ. 1) THEN
               IPRI = IPRI + 1
               IPNTP(NMAT,IPRI) = JSTRA + J
            END IF
            END DO
            NDIMP(NMAT) = IPRI
C
C           IPNTC and NDIMC
C
            ICNT = 0
            DO J = 1, NRCA
            IF (MAT(KPRIM,J) .EQ. 1) THEN
               ICNT = ICNT + 1
               IPNTC(NMAT,ICNT) = J
            END IF
            END DO
            NDIMC(NMAT) = ICNT
C
         END IF
      END DO
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from AOBTCH1',-1)
         WRITE (LUPRI,'(1X,A,I3,A)')
     &      ' Input contraction matrix for block', NB,':'
         CALL OUTPUT(PRICCF(JSTRA+1,1),1,NUCA,1,NRCA,
     &               MXPRIM,NRCA,1,LUPRI)
         WRITE (LUPRI,'(/,1X,A,I3)') ' Number of AO batches:',NMAT
         DO I = 1, NMAT
            WRITE (LUPRI,'(1X,A,I3,A,2I3)')
     &      ' Dimensions (prim. - cont.) for batch',I,':',
     &        NDIMP(I),NDIMC(I)
            WRITE (LUPRI,'(1X,A,17I3/,(29X,17I3))')
     &      ' Contributing primitives:  ',(IPNTP(I,J),J = 1,NDIMP(I))
            WRITE (LUPRI,'(1X,A,17I3/,(29X,17I3))')
     &      ' Contributing contracted:  ',(IPNTC(I,J),J = 1,NDIMC(I))
         END DO
      END IF
C
C     Assign to COMMON /FCK3AOBTCH/
C
      IF (ISHELL .EQ. 1) THEN
         IADRP = 1
         IADRC = 1
         IADRO(1) = 1
         IADRO(2) = 1
         IADRO(3) = 1
         NBTCH = 0
      END IF
C
      DO IMAT = 1, NMAT
         IBATCH = NBTCH + IMAT

         NCNTBTCH(IBATCH)  = NCENT(ISHELL)
         NHKTBTCH(IBATCH)  = NHKT(ISHELL)

C
C        We first set type to 1 (check for type 2 is afterwards )
C
         NTYPBTCH(IBATCH)  = 1

         NCTFBTCH(IBATCH)  = NDIMC(IMAT)
         NPRFBTCH(IBATCH)  = NDIMP(IMAT)

         MAXAOPRF = MAX( MAXAOPRF, NDIMP(IMAT) )
         MAXAOCTF = MAX( MAXAOCTF, NDIMC(IMAT) )
         KCLSBTCH(IBATCH) = LCLASS(ISHELL)
C
C        exponents
C
         KEXPBTCH(IBATCH) = IADRP
         DO IPRI = 1, NDIMP(IMAT)
            EXPBTCH(IADRP) = PRIEXP(IPNTP(IMAT,IPRI))
            IF (.NOT.LCLASS(ISHELL).EQ.0)
     &         EXPAOMIN = MIN( EXPAOMIN, EXPBTCH( IADRP ) )
            IADRP = IADRP + 1
         END DO
C
C        contraction coefficients
C        ( second index we initialize to 1 since even if we do not 
C          need it it might be pass as the index in CCFBT into 
C          subroutines )
C          
C
         KCCFBTCH(1,IBATCH) = IADRC
         KCCFBTCH(2,IBATCH) = 1

         DO ICNT = 1, NDIMC(IMAT)
         DO IPRI = 1, NDIMP(IMAT)
            CCFBTCH(IADRC) = PRICCF(IPNTP(IMAT,IPRI),IPNTC(IMAT,ICNT))
            CCFAOMAX = MAX( CCFAOMAX, ABS( CCFBTCH(IADRC) ) )
            IADRC = IADRC + 1
         END DO
         END DO
C
C        Loop over contracted functions
C        NB : to be generalized (we assume as many s.c contracted as l.c. contracted)
C
         DO IKB = 1, 3
          KNDXBTCH(1,IBATCH,IKB) = IADRO(IKB)
          KNDXBTCH(2,IBATCH,IKB) = 1
          DO ICNT = 1, NDIMC(IMAT)
C            pointer to the start of the contracted functions for this class
             INDXBTCH(IADRO(IKB),IKB) = IORB(IKB) +
     &                             (IPNTC(IMAT,ICNT)-1)*KHKT_REL(IKB)
             IADRO(IKB) = IADRO(IKB) + 1
          END DO
         END DO
C
C        AO index
C
         KAOBTCH(1,IBATCH) = IAOIND
         KAOBTCH(2,IBATCH) = 1
         NHK    = NHKTBTCH(IBATCH)
         KHK    = NHK * 2 - 1
         IAOIND  = IAOIND  + KHK * NCTFBTCH(IBATCH)
      END DO
C
      NBTCH   = NBTCH + NMAT
      NAOBTCH = NBTCH
C
      RETURN
      END

C
C ---------------------------------------------------------------
C

      SUBROUTINE SETAOCLS(IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"

      ICLS = 1 
      KSTAOB = 1
 
      NAOB  = 0
      ICNT0 = NCNTBTCH(1)
      LCLASS0 = KCLSBTCH(1)


      DO IAOB = 1, NAOBTCH
       ICNT = NCNTBTCH(IAOB)
       LCLASS = KCLSBTCH(IAOB)

       IF ( (ICNT .NE. ICNT0) .OR. (LCLASS .NE. LCLASS0) )  THEN 
         
         IF ((LCLASS .LT. LCLASS0).AND.(LCLASS.NE.0)) THEN
           Write(LUPRI,*) 'LCLASS changed from ', LCLASS0, ' to ', 
     &                     LCLASS, '.'
           Write(LUPRI,*) 'Shells should be ordered according to '//
     &                    'LCLASS. LCLASS = {1,2,3,0} '//
     &                    '= {Large, Small, Auxilliary, 0-exponent'
           CALL QUIT('Shells should be ordered according to '//
     &               'LCLASS = {1,2,3,0}.')
         END IF

        IAOCLS(1,ICLS) = ICNT0
        IAOCLS(2,ICLS) = KSTAOB
        IAOCLS(3,ICLS) = NAOB
        IAOCLS(4,ICLS) = LCLASS0

        ICLS = ICLS + 1
        KSTAOB = IAOB
        NAOB   = 0
        LCLASS0 = LCLASS

        ICNT0  = ICNT
       END IF

       NAOB = NAOB + 1
      END DO

      IAOCLS(1,ICLS) = ICNT0
      IAOCLS(2,ICLS) = KSTAOB
      IAOCLS(3,ICLS) = NAOB
      IAOCLS(4,ICLS) = LCLASS0

      NAOBCLS = ICLS 

      RETURN
      END

C
C ---------------------------------------------------------------
C

      SUBROUTINE SETSPEXP(IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3thr.h"

      IAOBLASTS = 1
      IAOBLAST  = 0
      
      DO IAOB = 1, NAOBTCH
       NHKT = NHKTBTCH( IAOB )

       NTYP = 1
       IF ( NHKT .EQ. 2 ) THEN
        NHKTLASTS = NHKTBTCH( IAOBLASTS )
        DO WHILE ( ( NHKTLASTS .NE. 1 ) .AND.
     &             ( IAOBLASTS .LT. IAOBLAST ) )
         IAOBLASTS = IAOBLASTS + 1
         NHKTLASTS = NHKTBTCH( IAOBLASTS )
        END DO
        DO WHILE ( ( NTYP .EQ. 1 ) .AND. ( NHKTLASTS .EQ. 1 )
     &       .AND. ( IAOBLASTS .LE. IAOBLAST ) )
         IF ( NPRFBTCH( IAOB ) .EQ. NPRFBTCH( IAOBLASTS )) THEN
           CALL CHECKSPEXP( NTYP, IAOB, IAOBLASTS )
         END IF
         IAOBLASTS = IAOBLASTS + 1
         NHKTLASTS = NHKTBTCH( IAOBLASTS )
        END DO
       END IF

       IF ( NTYP .EQ. 2 ) THEN 
         NTYPBTCH( IAOBLASTS-1 ) = 2
         NHKTBTCH( IAOBLASTS-1 ) = 2
         KNDXBTCH( 2, IAOBLASTS-1, 1 ) = KNDXBTCH( 1, IAOB, 1 )
         KNDXBTCH( 2, IAOBLASTS-1, 2 ) = KNDXBTCH( 1, IAOB, 2 )
         KNDXBTCH( 2, IAOBLASTS-1, 3 ) = KNDXBTCH( 1, IAOB, 3 )
         KCCFBTCH( 2, IAOBLASTS-1 ) = KCCFBTCH( 1, IAOB )
         KCLSBTCH( IAOBLAST - 1 )   = KCLSBTCH(IAOB)
         KAOBTCH(  2, IAOBLASTS-1 ) = KAOBTCH(  1, IAOB )     
       ELSE
         IAOBLAST = IAOBLAST + 1
         IF ( IAOBLAST .LT. IAOB ) THEN 
          NTYPBTCH( IAOBLAST )    = 1
          NPRFBTCH( IAOBLAST )    = NPRFBTCH( IAOB )
          NCTFBTCH( IAOBLAST )    = NCTFBTCH( IAOB )
          NHKTBTCH( IAOBLAST )    = NHKTBTCH( IAOB )
          NCNTBTCH( IAOBLAST )    = NCNTBTCH( IAOB )
          KEXPBTCH( IAOBLAST )    = KEXPBTCH( IAOB )
          KNDXBTCH( 1, IAOBLAST, 1 ) = KNDXBTCH( 1, IAOB, 1 )
          KNDXBTCH( 1, IAOBLAST, 2 ) = KNDXBTCH( 1, IAOB, 2 )
          KNDXBTCH( 1, IAOBLAST, 3 ) = KNDXBTCH( 1, IAOB, 3 )
          KNDXBTCH( 2, IAOBLAST, 1 ) = 1
          KNDXBTCH( 2, IAOBLAST, 2 ) = 1
          KNDXBTCH( 2, IAOBLAST, 3 ) = 1
          KCCFBTCH( 1, IAOBLAST ) = KCCFBTCH( 1, IAOB )
          KCCFBTCH( 2, IAOBLAST ) = 1
          KCLSBTCH( IAOBLAST )    = KCLSBTCH(IAOB)
          KAOBTCH(  1, IAOBLAST ) = KAOBTCH(  1, IAOB )     
          KAOBTCH(  2, IAOBLAST ) = 1
         END IF
       END IF
         
      END DO

      NAOBTCH = IAOBLAST

      RETURN
      END

      SUBROUTINE CHECKSPEXP(NTYP, IAOB1, IAOB2)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3thr.h"

      NPRF1 = NPRFBTCH( IAOB1 )
      NCTF1 = NCTFBTCH( IAOB1 )
      NCNT1 = NCNTBTCH( IAOB1 )
      KEXP1 = KEXPBTCH( IAOB1 )

      NPRF2 = NPRFBTCH( IAOB2 )
      NCTF2 = NCTFBTCH( IAOB2 )
      NCNT2 = NCNTBTCH( IAOB2 )
      KEXP2 = KEXPBTCH( IAOB2 )

      NTYP = 1

      IF ( ( NPRF1 .EQ. NPRF2 ) .AND. ( NCTF1 .EQ. NCTF2 ) .AND.
     &     ( NCNT1 .EQ. NCNT2 ) ) THEN 

       EXPDIFSUM  = 0.0D0
       DO IPRF1=1,NPRF1
         EXPDIFSUM = EXPDIFSUM + 
     &     ABS(EXPBTCH(KEXP1+IPRF1-1)-EXPBTCH(KEXP2+IPRF1-1))
       END DO
      
       IF ( EXPDIFSUM .LT. THRSPEXP ) NTYP = 2

      END IF

      RETURN 
      END 
            
C
C -----------------------------------------------------
C

      SUBROUTINE GETAOB(NPRF,NCTF,NTYP,NHKT,JMAX,KHKT,KCKT,
     &                  KCSQ,KEXP,KCCF,IAOB )
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
      DIMENSION KCCF(2)

      NPRF  = NPRFBTCH( IAOB )
      NCTF  = NCTFBTCH( IAOB )
      NHKT  = NHKTBTCH( IAOB )
      NTYP  = NTYPBTCH( IAOB ) 
      KEXP  = KEXPBTCH( IAOB ) 
      KCCF(1) = KCCFBTCH( 1, IAOB )
      KCCF(2) = KCCFBTCH( 2, IAOB )

      JMAX  = NHKT - 1
      KHKT  = 2 * NHKT - 1 
      KCKT  = NHKT * ( NHKT + 1 ) / 2 

      KCSQ  = KSQADR(JMAX,0,0,0,0)

      RETURN
      END 


      SUBROUTINE DECONTRACTAO(WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "iratdef.h"
      DIMENSION WORK(LWORK)

      KKNDXTMP = 1 
      KNHKTTMP = KKNDXTMP + NAOBTCH*6/IRAT + 1
      KKEXPTMP = KNHKTTMP + NAOBTCH  /IRAT + 1
      KKCCFTMP = KKEXPTMP + NAOBTCH  /IRAT + 1
      KNTYPTMP = KKCCFTMP + NAOBTCH*2/IRAT + 1
      KLAST    = KNTYPTMP + NAOBTCH  /IRAT + 1

      IF ( KLAST .GT. LWORK ) THEN 
        CALL STOPIT('SETAOBTCH','DECONTRACTAO',KLAST,LWORK)
      END IF 

      CALL ICOPY(6*NAOBTCH,KNDXBTCH,1,WORK(KKNDXTMP),1) 
      CALL ICOPY(  NAOBTCH,NHKTBTCH,1,WORK(KNHKTTMP),1)
      CALL ICOPY(  NAOBTCH,KEXPBTCH,1,WORK(KKEXPTMP),1)
      CALL ICOPY(2*NAOBTCH,KCCFBTCH,1,WORK(KKCCFTMP),1)
      CALL ICOPY(  NAOBTCH,NTYPBTCH,1,WORK(KNTYPTMP),1)

      CALL DECONTRACTAO1(WORK(KKNDXTMP),WORK(KNHKTTMP),
     &                   WORK(KKEXPTMP),WORK(KKCCFTMP),
     &                   WORK(KNTYPTMP),IPRINT)

      RETURN
      END

      SUBROUTINE DECONTRACTAO1( KNDXTMP, NHKTTMP, KEXPTMP,
     &                          KCCFTMP, NTYPTMP, IPRINT )
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "iratdef.h"
      DIMENSION KNDXTMP(2,MXAOBTCH,3),NHKTTMP(MXAOBTCH),
     &          KEXPTMP(MXAOBTCH),KCCFTMP(2,MXAOBTCH),
     &          NTYPTMP(MXAOBTCH)

      IAOB = 1

      DO ICLS = 1, NAOBCLS
       KSTAOB = IAOB 

       KSTAOB_OLD = IAOCLS(2,ICLS)
       NAOB_OLD   = IAOCLS(3,ICLS)

       DO IAOB_OLD = KSTAOB_OLD, KSTAOB_OLD + NAOB_OLD - 1
        IF ( NCTFBTCH(IAOB_OLD) .GT. 1 ) THEN
         CALL QUIT('Uncontr. not implemented for general contractions.')
        END IF
        IF ( IAOB + NPRFBTCH(IAOB_OLD) .GT. MXAOBTCH ) THEN
         CALL QUIT('More uncontracted AO batches than MXAOBTCH !!!')
        END IF 

        DO I=1,NPRFBTCH(IAOB_OLD)
         KNDXBTCH(1,IAOB,1) = KNDXTMP(1,IAOB_OLD,1) 
         KNDXBTCH(2,IAOB,1) = KNDXTMP(2,IAOB_OLD,1) 
         NHKTBTCH(IAOB)   = NHKTTMP(IAOB_OLD)
         KEXPBTCH(IAOB)   = KEXPTMP(IAOB_OLD)   + I - 1
         KCCFBTCH(1,IAOB) = KCCFTMP(1,IAOB_OLD) + I - 1
         KCCFBTCH(2,IAOB) = KCCFTMP(2,IAOB_OLD) + I - 1
         NTYPBTCH(IAOB)   = NTYPTMP(IAOB_OLD)

         IAOB = IAOB + 1
        END DO
       END DO

       NAOB = IAOB - KSTAOB

       IAOCLS(2,ICLS) = KSTAOB
       IAOCLS(3,ICLS) = NAOB
      END DO

      NAOBTCH = IAOB - 1

      DO IAOB = 1, NAOBTCH 
       NPRFBTCH(IAOB) = 1
       NCTFBTCH(IAOB) = 1
      END DO

      RETURN
      END 


C     Routine for normalizing the deinsity-fitting basis wrt. the Coulomb
C     self-integrals (alpha|alpha)
      SUBROUTINE NRMFIT(GABALP,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      DIMENSION WORK(LWORK), GABALP(*)
#include "shells.h"
#include "symmet.h"
#include "nuclei.h"
#include "fck3aob.h"
#include "fck3mem.h"
#include "ibtfun.h"
      CALL QENTER('NRMFIT')

      IORB = 1
      IAOIND = 1
      NB = -1
      NC = -1
      ND = -1

      DO ISHELL = 1, KMAXTOT
        IF (LCLASS(ISHELL).EQ.3) THEN
          IF (NB.NE.NBCH(ISHELL) .OR. NC.NE.NCENT(ISHELL) .OR.
     &        ND.NE.LCLASS(ISHELL) ) THEN
           NB = NBCH(ISHELL)
           NC = NCENT(ISHELL)
           ND = LCLASS(ISHELL)
           
           CALL NRMFIT1(GABALP,IORB,IAOIND,NB,ISHELL,WORK,LWORK,IPRINT)
          END IF
          IORB = IORB + KHKT(ISHELL)
       ENDIF 
      END DO
      CALL QEXIT('NRMFIT')
      RETURN
      END
C
C ---------------------------------------------------------------
C

      SUBROUTINE NRMFIT1(GABALP,IORB,IAOIND,NB,ISHELL,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      DIMENSION GABALP(*), WORK(LWORK)
#include "shells.h"
      CALL QENTER('NRMFIT1')

C
      NUCA  = NUCO(ISHELL)
      NRCA  = NRCO(ISHELL)
      JSTRA = JSTRT(ISHELL)
C
      KMAT   = 1
      KDONE  = KMAT   + NUCA*NRCA
      KDIMP  = KDONE  + NUCA
      KDIMC  = KDIMP  + NUCA
      KPNTC  = KDIMC  + NUCA
      KPNTP  = KPNTC  + NUCA*NRCA
      KLAST  = KPNTP  + NUCA*NUCA

      IF (KLAST .GT. LWORK) CALL STOPIT('NRMFIT1',' ',KLAST,LWORK)
C
      CALL NRMFIT2(GABALP,IAOIND,NB,NMAT,WORK(KMAT),NUCA,NRCA,
     &             JSTRA,WORK(KDONE),
     &             WORK(KDIMP),WORK(KDIMC),WORK(KPNTC),WORK(KPNTP),
     &             ISHELL,IORB,IPRINT)
C
      CALL QEXIT('NRMFIT1')
      RETURN
      END

C
C -----------------------------------------------------
C

      SUBROUTINE NRMFIT2(GABALP,IAOIND,NB,NMAT,MAT,NUCA,NRCA,JSTRA,
     &                   IDONE,NDIMP,
     &                   NDIMC,IPNTC,IPNTP,ISHELL,IORB,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
      PARAMETER (THR = 1.0D-17)
C
#include "symmet.h"
#include "fck3aob.h"
#include "shells.h"
#include "primit.h"
#include "nuclei.h"
C
      DIMENSION MAT(NUCA,NRCA), IDONE(NUCA), NDIMP(NUCA), NDIMC(NUCA),
     &          GABALP(NBASISAUX),IPNTC(NUCA,NRCA),IPNTP(NUCA,NUCA)
C
      SAVE IADRP, IADRO, NBTCH
C
C aspg, 2006-06-27 making ftnchek happy
C
      NBTCH = 0
C
C     Construct matrix MAT containing 1 and 0 only, then fill
C     out all elements within each contraction matrix as in SETAOBTCH
C
      CALL IZERO(MAT,NUCA*NRCA)
C
      DO I = 1, NUCA
      DO J = 1, NRCA
         IF (ABS(PRICCF(JSTRA+I,J)).GT.THR) MAT(I,J) = 1
      END DO
      END DO
C
      DO I = 1, NRCA
      DO J = 1, NRCA
         IF (INPRD(MAT(1,I),MAT(1,J),NUCA) .GT. 0) THEN
            DO K = 1, NUCA
               MAT(K,I) = MAX(MAT(K,I),MAT(K,J))
               MAT(K,J) = MAX(MAT(K,I),MAT(K,J))
            END DO
         END IF
      END DO
      END DO
C
C     NMAT, NDIMP(), NDIMC(), IPNTP(), IPNTC()
C
      NMAT = 0
      CALL IZERO(IDONE,NUCA)
      DO I = 1, NRCA
         IF (IDONE(I).EQ.0) THEN
C
            NMAT = NMAT + 1
            DO J = I, NRCA
               IDONE(J) =  INPRD(MAT(1,I),MAT(1,J),NUCA)
            END DO
C
C           First non-zero primitive
C
            KCONT = I
            DO J = 1, NUCA
               IF (MAT(J,KCONT) .EQ. 1) THEN
                  KPRIM = J
                  GO TO 100
               END IF
            END DO
  100       CONTINUE
C
C           IPNTP and NDIMP
C
            IPRI = 0
            DO J = 1, NUCA
            IF (MAT(J,KCONT) .EQ. 1) THEN
               IPRI = IPRI + 1
               IPNTP(NMAT,IPRI) = JSTRA + J
            END IF
            END DO
            NDIMP(NMAT) = IPRI
C
C           IPNTC and NDIMC
C
            ICNT = 0
            DO J = 1, NRCA
            IF (MAT(KPRIM,J) .EQ. 1) THEN
               ICNT = ICNT + 1
               IPNTC(NMAT,ICNT) = J
            END IF
            END DO
            NDIMC(NMAT) = ICNT
C
         END IF
      END DO

C     Renormalization loop over blocks
      DO IMAT = 1, NMAT
C
C        AO index
C
         NHK    = NHKT(ISHELL)
         KHK    = NHK * 2 - 1

C        Renormalization
         DO ICNT = 1, NDIMC(IMAT)
         DO IPRI = 1, NDIMP(IMAT)
            PRICCF(IPNTP(IMAT,IPRI),IPNTC(IMAT,ICNT)) = 
     &          PRICCF(IPNTP(IMAT,IPRI),IPNTC(IMAT,ICNT))/GABALP(IAOIND)
         END DO
         END DO
         IAOIND  = IAOIND  + KHK * NDIMC(IMAT)

      END DO
C
      NBTCH   = NBTCH + NMAT
      NAOBTCH = NBTCH
C
      RETURN
      END

C     For debugging
      FUNCTION NACTIVE_BATCHES()
#include "implicit.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
      N=0
      DO I = 1, NAOBTCH 
        IF (ACTVAOBTCH(I)) N=N+1
      ENDDO
      NACTIVE_BATCHES = N
      RETURN
      END
C  /* Deck pick_act_batch */
      SUBROUTINE ACT_ALL_BATCHES
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <aovec.h>
#include <fck3aob.h>
C aspg, 2006-07-25
C putting outer do loop in line with dalton
C     DO I = 1, MXSHEL
      DO I = 1, NAOBTCH 
        ACTVAOBTCH(I) = .TRUE.
        DO J = 1, 4
          ACTVBTH(I,J) = .TRUE.
        ENDDO
      ENDDO
      RETURN
      END
