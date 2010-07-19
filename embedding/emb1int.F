!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!dirac_copyright_end

C VERSION : $Revision$
C DATE    : $Date$
C FILE    : emb1int.F
C
C
C  /* Deck oneprm */
      MODULE VC_GRID

      PRIVATE
      PUBLIC EMB_VC
      
C     Module for calculating the coulomb potential on a grid.
C     Contains:
C     EMB_VC    The interface for users. 
C     LOOP_PRIM Loops over primitive orbitals.
C     VC_1PRIM  The routine that calls for one primitive shell VC and VC_CART
C               to calculate the shell's contribution to the coulomb potential.
C     VC        Calculates the coulomb integrals in Hermite basis
C     VC_CART   Transforms the coulomb integrals to a cartesian basis.
C     
C
C     Calls the Dirac F77 routines: GETODC in abacus/her1odc.F
C     Based on the routines: HERNAI, CINT0 in abacus/her1car.F
C     The loop over the primitives is based on abacus/her1drv.F and abacus/her1sym.F
C
C     The include files that need to be set:
C     onecom.h
C     sdepre.h
C

      CONTAINS

      SUBROUTINE EMB_VC(DMAT,points,vc,IPRINT,WORK,LWORK)
C     Interface for the user. Calles LOOP_PRIM
C     calls VC_PRIM.
C     Needed to make contact with F77.
C     JHS DMAT is in MO basis. Transformations needed.

      IMPLICIT NONE
#include "priunit.h"
#include "maxaqn.h"
#include "ccom.h"

C     PASSED variables
      real(kind=8) :: points(:,:),vc(:)
      INTEGER :: IPRINT
      REAL(KIND=8) :: DMAT(*)
      REAL(KIND=8) :: WORK(*)
      INTEGER :: LWORK

C     JHS MAKE LATER OPTIONAL
      REAL(KIND=8) :: TOLOG, TOLS,FACTOR=1.0
      
      TOLS=THRS**2
      TOLOG=2*LOG(THRS)

      IF (IPRINT .GE. 5) CALL TITLER('Output from EMB_VC','*',103)

C     JHS Transform Density matrix to AO basis.

      CALL LOOP_PRIM(DMAT,points,vc,WORK,LWORK,
     &     IPRINT,TOLOG,TOLS)

      END SUBROUTINE

      SUBROUTINE LOOP_PRIM(DMAT,points,vc,WORK,LWORK,IPRINT,TOLOG,
     &                  TOLS)
C     Driver routine that loops over the primitive orbitals and
C     calls VC_PRIM.
      USE VECTOR_FUN, ONLY: NORM
      IMPLICIT NONE

#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#endif
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"

C     PASSED variables
      real(kind=8) :: points(:,:),vc(:)
      INTEGER :: IPRINT
      REAL(KIND=8) :: DMAT(NTBAS(0),NTBAS(0),NZ)
      REAL(KIND=8) :: WORK(*)
      REAL(KIND=8) :: TOLOG, TOLS,FACTOR=1.0
      INTEGER :: LWORK

C     Local Variables
      INTEGER :: ISHELA,ISHELB,ICA,ICB,NA,IREP,NB
      INTEGER :: MULTA,MULTB,NHKTAB,KAB,IDENB_TMP
      INTEGER :: ISYM
      INTEGER :: ipoint, i,j,k,l

C      REAL(KIND=8) :: CORA(3),CORB(3)
      REAL(KIND=8),ALLOCATABLE :: ADER(:,:,:)
C     Contains the coulomb integrals for one set of primitive orbitals.

C     JHS TMP STUFF FOR DEBUGGING
      REAL(KIND=8) :: gexp=0.0,vc_tmp
      REAL(KIND=8) :: dummy
      INTEGER(KIND=8) :: idummy

#include "cbisol.h"
#include "shells.h"
C#include "inforb.h"
#include "onecom.h"
#include "primit.h"
#include "symind.h"
#include "symmet.h"
#include "lmns.h"
#include "nuclei.h"

C     Initialization
      ISOFRA(:)=0


C     Loop over shells.
      IDENA=0
      DO ISHELA = 1,KMAX
        NHKTA = NHKT(ISHELA)
        KHKTA = KHKT(ISHELA)
        KCKTA = KCKT(ISHELA)
        ICA   = LCLASS(ISHELA)
        SPHRA = SPHR(ISHELA)
        CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
        NCENTA = NCENT(ISHELA)
        ICENTA = NUCNUM(NCENTA,1)
        MULA   = ISTBAO(ISHELA)
        MULTA  = MULT(MULA)
        NUCA   = NUCO(ISHELA)
        NUMCFA = NUMCF(ISHELA)
        JSTA   = JSTRT(ISHELA)
        CORAX  = CENT(ISHELA,1,1)
        CORAY  = CENT(ISHELA,2,1)
        CORAZ  = CENT(ISHELA,3,1)
C
C       Compute symmetry integral pointers for contributions
C       from this block.  Note that at present this assumes all
C       components from a shell are included.
C
        IDENB = 0
        ISOFRB(:) = 0
        INDFA(:,:) = -10 000 000
        DO NA = 1, KHKTA
          DO IREP = 0, MAXREP
          IF (IAND(MULA,IEOR(IREP,ISYMAO(NHKTA,NA)))==0) THEN
            ISOFRA(IREP+1)    = ISOFRA(IREP+1) + 1
            INDFA (IREP+1,NA) = ISOFRA(IREP+1)
          END IF
          END DO
        END DO
        DO ISHELB = 1,ISHELA
          LDIAG = ISHELA .EQ. ISHELB
          NHKTB = NHKT(ISHELB)
          KHKTB = KHKT(ISHELB)
          KCKTB = KCKT(ISHELB)
          ICB   = LCLASS(ISHELB)
          SPHRB = SPHR(ISHELB)
          CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
          NCENTB = NCENT(ISHELB)
          NHKTAB = NHKTA + NHKTB
          MULB   = ISTBAO(ISHELB)
          MULTB  = MULT(MULB)
          NUCB   = NUCO(ISHELB)
          NUMCFB = NUMCF(ISHELB)
          JSTB   = JSTRT(ISHELB)
C     JHS Later we need to multiply with the sign for the symmetry operations
          CORBX  = CENT(ISHELB,1,1)
          CORBY  = CENT(ISHELB,2,1)
          CORBZ  = CENT(ISHELB,3,1)
          KHKTAB = KHKTA*KHKTB
          KCKTAB = KCKTA*KCKTB
          MAB    = IOR(MULA,MULB)
          KAB    = IAND(MULA,MULB)
          HKAB   = FMULT(KAB)
C
          SPHRAB = SPHRA .OR. SPHRB
C     JHS Note that if SPHRAB is true we should correct the densitymatrix for this.
C
C        Compute symmetry integral pointers for contributions
C        from this block.  Note that at present this assumes all
C        components from a shell are included
C
          INDFB(:,:) = -10 000 000
          DO NB = 1, KHKTB
            DO IREP = 0, MAXREP
            IF (IAND(MULB,IEOR(IREP,ISYMAO(NHKTB,NB)))==0) THEN
              ISOFRB(IREP+1)    = ISOFRB(IREP+1) + 1
              INDFB (IREP+1,NB) = ISOFRB(IREP+1)
            END IF
            END DO
          END DO
C        jth - CLASS (which integral class: LL,LS,SL,SS)
C                                           0  1  1  2
C     JHS FOLLOWING STATEMENT NEEDS DEBUGGING FOR OTHER COMPILERS
         IF (ICA==1 .AND. ICB==1) THEN; INTCLASS=2
         ELSE IF (ICA==1.OR.ICB==1 ) THEN; INTCLASS=1
         ELSE; INTCLASS=0
         END IF

         IF (IPRINT>5) 
     &        WRITE(LUPRI,*) 'JHS ICA,ICB,INTCLASS', ICA,ICB,INTCLASS

         IF (IPRINT .GE. 05) WRITE (LUPRI, *) ISHELA, ISHELB
         IF (IPRINT .GE. 10) THEN
             WRITE (LUPRI,'(A,2I10)') ' NRCO',NRCO(ISHELA),NRCO(ISHELB)
             WRITE (LUPRI,'(A,2I10)') ' NHKT   ', NHKTA, NHKTB
             WRITE (LUPRI,'(A,2I10)') ' KHKT   ', KHKTA, KHKTB
             WRITE (LUPRI,'(A,2I10)') ' KCKT   ', KCKTA, KCKTB
             WRITE (LUPRI,'(A,2I10)') ' NCENT  ', NCENTA, NCENTB
             WRITE (LUPRI,'(A,2I10)') ' ISTBAO ', MULA, MULB
             WRITE (LUPRI,'(A,2I10)') ' MULT   ', MULTA, MULTB
             WRITE (LUPRI,'(A,2I10)') ' NUC    ', NUCA, NUCB
             WRITE (LUPRI,'(A,2I10)') ' NUMCF  ', NUMCFA, NUMCFB
             WRITE (LUPRI,'(A,2I10)') ' JST    ', JSTA, JSTB
             WRITE (LUPRI,'(A,2I10)') ' NBCH',NBCH(ISHELA),NBCH(ISHELB)
             WRITE (LUPRI,'(A,2I10)') ' NCENT  ', NCENTA,NCENTB
             WRITE (LUPRI,'(A,2F12.6)') ' CORAX    ', CORAX, CORBX
             WRITE (LUPRI,'(A,2F12.6)') ' CORAY    ', CORAY, CORBY
             WRITE (LUPRI,'(A,2F12.6)') ' CORAZ    ', CORAZ, CORBZ
             WRITE (LUPRI,'(A,2I10)') ' IDEN   ',IDENA,IDENB           
         END IF
C
C        Initialization for nuclear attraction integrals
C
         JMAX = NHKTAB - 2
         ISTEPU = JMAX + 1
         ISTEPV = ISTEPU*ISTEPU
         NAHGTF = ISTEPU*ISTEPV

         ALLOCATE(ADER(KCKTA,KCKTB,SIZE(points,2)))
C     Contains ader for one gridpoint in SO basis.


         IF (IPRINT>5) WRITE(LUPRI,*) 'JHS SYMMETRY,SYMOP',MAXOPR
C     JHS Now loop over symmetry dependend stuff.

C     JHS Maybe IDENB needs to be corrected by -KHKTB
         IDENB_TMP=IDENB-KHKTB
C     JHS IDENB_TMP variable necessary?
         DO ISYM=0,MAXOPR
            IF (IAND(ISYM,MAB)/=0) THEN
              IF (IAND(ISYM,MULB)==0) IDENB_TMP=IDENB+KHKTB
              CYCLE
            END IF
            IDENB_TMP=IDENB_TMP+KHKTB
            ICENTB=NUCNUM(NCENTB,ISYM+1)
            ONECEN=(ICENTA==ICENTB)
            SIGNBX=PT(IAND(ISYMAX(1,1),ISYM))
            SIGNBY=PT(IAND(ISYMAX(2,1),ISYM))
            SIGNBZ=PT(IAND(ISYMAX(3,1),ISYM))
            
C     JHS Correct the B coordinate for the symmetry operations.
C         Direct passing is not possible as CORA,CORB are used inside her1odc.F
            CORBX=CORBX*SIGNBX
            CORBY=CORBY*SIGNBY
            CORBZ=CORBZ*SIGNBZ
 
            IF (IPRINT>5) WRITE(LUPRI,*) 'JHS CORB',CORBX,CORBY,CORBZ
            
            CALL VC_SHELL(ADER,TOLOG,TOLS,points,IPRINT,WORK,
     &      LWORK)
            
C     JHS Transform to SO basis?

            IF (IPRINT>5) WRITE(LUPRI,*) 'JHS ADER(VC_SHELL)',ADER

C     JHS Start contraction with the density matrix.

C     JHS Density matrix is in AO basis. 
C     Contraction of ADER to AO basis done in CART_VC (in ODC coefficients)
C     Maybe move inside PRIM_VC
C     Similar to dirac/diraos.F:GETGAO and dirac/dirden.F:DENMAT
C
C     Contract ADER with density matrix, indices in ADER are swapped.
            IF (IPRINT>10) THEN
                 WRITE(LUPRI,*) 'DMAT Contraction.',ISHELA,ISHELB

                 WRITE(LUPRI,*) ISHELA,ISHELB,KCKTA,KCKTB,IDENA,
     &           IDENB_TMP,NZ,SPHRA,SPHRB            
            END IF

       DO ipoint=1,SIZE(vc)               
C     Contract with the right symmetry part of the density matrix.
C     This is given by INDFA and INDFB


C     No contraction over the non-real quaternionic part. Use index matrices of quaternions
               vc(ipoint)=vc(ipoint)+SUM(
     &              DMAT(idena+1:idena+KCKTA,idenb+1:idenb+KCKTB,1)*
     &              ADER(1:KCKTA,1:KCKTB,ipoint))

            END DO 
           

C     JHS Transform back to the original B coordinate.
            CORBX=CORBX/SIGNBX
            CORBY=CORBY/SIGNBY
            CORBZ=CORBZ/SIGNBZ
         END DO

         DEALLOCATE(ADER)


         IDENB = IDENB + KHKTB*MULTB
       END DO
       IDENA = IDENA + KHKTA*MULTA
      END DO
      

      DO ipoint=1,SIZE(vc)
        vc_tmp=0.0
        DO j=1,MXCENT          
          vc_tmp=vc_tmp-CHARGE(j)/(NORM(CORD(:,j)-
     &          points(1:3,ipoint)))
        END DO
        vc(ipoint)=vc(ipoint)-vc_tmp
      END DO

      END SUBROUTINE




      SUBROUTINE VC_SHELL(ADER,TOLOG,TOLS,points,IPRINT,WORK,
     &     LWORK)
C     Calculates the contribution for one primitive orbital set.
      USE VECTOR_FUN, ONLY: NORM
      IMPLICIT NONE

#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#endif
C     PASSED variables
      INTEGER :: IPRINT
      real(kind=8) :: points(:,:)
      REAL(KIND=8) :: ADER(KCKTA,KCKTB,SIZE(points,2))
      REAL(KIND=8) :: WORK(*)
      REAL(KIND=8) :: TOLOG, TOLS,FACTOR=1.0
      INTEGER :: LWORK


C     Local variables
      REAL(KIND=8),PARAMETER :: D0=0.0D0, D1=1.00 D00, D3=3.00 D00, 
     &                          THIRD=D1/D3, GNUTHR=1.D-4,D1P5=1.5D0
      REAL(KIND=8)::AHGTF(0:JMAX,0:JMAX,0:JMAX)
      REAL(KIND=8),ALLOCATABLE::ODC(:,:,:,:,:,:)

      LOGICAL :: KINODC
      REAL(KIND=8) :: DIFAB(3),DIFCP(3),CORP(3),ORIGIN(3),
     &     CORA(3),CORB(3)
      REAL(KIND=8) :: DISTAB,CONTA,EXPA,CONTB,EXPB,EXPP,EXPPI,
     &                EXPABQ,SAAB, ASAAB,SAAB13, EXPAPI,EXPBPI
      INTEGER :: JMAXA,JMAXB,JMAXD,JMAXT,JMAXM,ipoint,
     &           IPRIMA,IPRIMB,JPRIMA,JPRIMB
      INTEGER :: IDUMMY
#include "cbisol.h"
#include "onecom.h"
#include "ader.h"
#include "primit.h"
     
      IF (IPRINT .GE. 5) CALL TITLER('Output from VC_SHELL','*',103)

C     ALLOCATION
      JMAXD = 2
      JMAXA = NHKTA - 1
      JMAXB = NHKTB - 1
      JMAXT = JMAXA + JMAXB + JMAXD
      JMAXM = 0
      ALLOCATE(ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3))

C     Initialization
      ADER(:,:,:)=0.0
      CORA(1)=CORAX
      CORA(2)=CORAY
      CORA(3)=CORAZ
      CORB(1)=CORBX
      CORB(2)=CORBY
      CORB(3)=CORBZ

      DIFAB(:) = CORA(:) - CORB(:)
      DISTAB=NORM(DIFAB)

C
C     ****************************************
C     ***** Loop over primitive orbitals *****
C     ****************************************
C
      IF (IPRINT>10) 
     &     WRITE(LUPRI,*) 'JHS Start looping over prim orbitals'
      DO IPRIMA = 1,NUCA
         JPRIMA = JSTA + IPRIMA
         CONTA = PRICCF(JPRIMA,NUMCFA)
         EXPA = PRIEXP(JPRIMA)
         DO IPRIMB = 1,NUCB
            JPRIMB = JSTB + IPRIMB
            CONTB = PRICCF(JPRIMB,NUMCFB)
            EXPB = PRIEXP(JPRIMB)
            EXPP = EXPA + EXPB
            EXPPI = D1/EXPP
C
C           Calculate and test square root factor
C     JHS SAAB are the contraction coefficients from the primitive orbitals to the basisfunctions.
            EXPABQ = EXPA*EXPB*DISTAB*EXPPI
            SAAB = CONTA*CONTB*EXP(-EXPABQ)
            ASAAB = ABS(SAAB)
            IF (IPRINT>20) WRITE(LUPRI,*) 'JHS INTSCREENING',
     &           'EXPABQ=',EXPABQ,'TOLOG=',TOLOG,
     &           'ASAAB=',ASAAB,'TOLS=',TOLS
            IF (EXPABQ<TOLOG) CYCLE
            IF (ASAAB<TOLS) CYCLE
            SAAB13 = SIGN(ASAAB**THIRD,SAAB)
C
C           Calculate coordinates of product Gaussian P
C
            EXPAPI = EXPA*EXPPI
            EXPBPI = EXPB*EXPPI
            CORP(:)=EXPAPI*CORA(:) + EXPBPI*CORB(:)
C
C           *********************************************
C           ***** Overlap Distribution Coefficients *****
C           *********************************************
C
C           Expansion coefficients for undifferentiated orbitals
C
#ifdef PRG_DIRAC
C           jth: removed KINODC = .NOT.PROPTY
            KINODC = .FALSE.
#else
            KINODC = .TRUE.
#endif 
            IDUMMY = 0
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,.FALSE.,
     &           KINODC,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,WORK,LWORK,
     &           CORP(1),CORP(2),CORP(3),.TRUE.,.FALSE.,ORIGIN,IDUMMY)
C
C           **********************************************
C           ***** Calculation of Hermitian integrals *****
C           **********************************************
C
C           Nuclear attraction integrals
C
            DO ipoint=1,SIZE(points,2)
              DIFCP(:) = points(1:3,ipoint) - CORP(:)
              CALL VNUC(AHGTF,JMAX,EXPP,DIFCP,FACTOR,IPRINT)
C
C           **********************************************
C           ***** Calculation of Cartesian integrals *****
C           **********************************************
C             
              CALL CART_VC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &             ADER(1:KCKTA,1:KCKTB,ipoint:ipoint),AHGTF,IPRINT)
            END DO
         END DO
      END DO

      DEALLOCATE(ODC)

C
C     Print
C
      
      RETURN
      END SUBROUTINE
C  /* Deck hernai */




      SUBROUTINE VNUC(AHGTF,JMAX,PVAL,CP,FACINT,IPRINT)
      USE VECTOR_FUN, ONLY: NORM
      IMPLICIT NONE
C
C     This subroutine calculates the R integrals as defined by
C     McMurchie and Davidson in J. Comp. Phys. 26 (1978) 218.
C     The recursion formulas (4.6) - (4.8) are used.
C
C     Only slightly slower than the implementation HERNAI in abacus/her1car.F of TUH 
C
C     JHS 26-03-08
C     
#include "priunit.h"
#include "pi.h"
#include "maxaqn.h"

C     Arguments passed
      REAL(KIND=8) :: AHGTF(0:,0:,0:),CP(3)
      REAL(KIND=8) :: PVAL,FACINT
      INTEGER :: JMAX,IPRINT

C     Local Variables
      REAL(KIND=8) :: RJ(0:JMAX,0:JMAX,0:JMAX,0:JMAX)
      REAL(KIND=8) :: FACTOR
      INTEGER :: JVAL,IT,IU,IV

!     Temporary variables
#include "gamcom.h"

!     Initialisation:
      RJ=0.0

      IF (IPRINT>10) THEN
         WRITE (LUPRI,'(/A/)') ' <<<<< Output from VC >>>>> '
         WRITE (LUPRI,'(A,I10)')   ' JMAX   ', JMAX
         WRITE (LUPRI,'(A,F12.6)') ' PVAL   ', PVAL
         WRITE (LUPRI,'(A,F12.6)') ' FACINT ', FACINT
         WRITE (LUPRI,'(A,3F12.6)') ' CP ', CP
      END IF
C
C     *************************************
C     ***** Incomplete Gamma Function *****
C     *************************************
C
      WVAL = PVAL*NORM(CP)**2
      JMAX0 = JMAX
      CALL GAMFUN
      IF (IPRINT==20) THEN
         CALL HEADER ('FJW after GAMFUN',-1)
         WRITE (LUPRI,'(6F12.6)') (FJW(JVAL), JVAL = 0, JMAX)
      END IF
      
!     RJ(0,0,0)
      FACTOR = FACINT*2.0*PI/PVAL
      DO JVAL = 0, JMAX
         FJW(JVAL)      =   FACTOR*FJW(JVAL)
         FACTOR         = - 2*PVAL*FACTOR
         RJ(0,0,0,JVAL) =   FJW(JVAL)
      END DO

!     Contraction scheme of Davidson. 
!     Slicing is not possible as we need to perform the calculations in order.
      DO JVAL = JMAX, 1, -1
         RJ(0,0,1,JVAL-1)= -CP(3)*RJ(0,0,0,JVAL)
         DO IV=1,JMAX-JVAL
!           Having one loop and using array slices is a tiny bit slower.
            RJ(0,0,IV+1,JVAL-1)= -CP(3)*RJ(0,0,IV,JVAL)+
     &                             IV*RJ(0,0,IV-1,JVAL)
!     RJ(0,1,v) etc. is moved outside the loop for efficiency
            RJ(0,1,IV,JVAL-1)= -CP(2)*RJ(0,0,IV,JVAL)
            DO IU=1,JMAX-JVAL-IV
               RJ(0,IU+1,IV,JVAL-1)= -CP(2)*RJ(0,IU,IV,JVAL)+
     &                                 IU*RJ(0,IU-1,IV,JVAL)
               RJ(1,IU,IV,JVAL-1)= -CP(1)*RJ(0,IU,IV,JVAL)
               DO IT=1,JMAX-JVAL-IU-IV
!     USE TRANSFER(IT,real_8) etc. if this is platform dependend.
                  RJ(IT+1,IU,IV,JVAL-1)= -CP(1)*RJ(IT,IU,IV,JVAL)+
     &                                     IT*RJ(IT-1,IU,IV,JVAL)
               END DO
            END DO
         END DO
      END DO

      IF (IPRINT .GT. 20) THEN
         CALL HEADER ('Hermite integrals in HERNAI',-1)
         WRITE (LUPRI,'(6F12.6)') AHGTF(:,:,:)
      END IF
      
      AHGTF(:,:,:)=RJ(0:,0:,0:,0)

      END SUBROUTINE

C  /* Deck cint0 */
      SUBROUTINE CART_VC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                 ADER,AHGTF,IPRINT)
C     TUH
      IMPLICIT NONE
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
C     Variables passed to CART_VC
      REAL(KIND=8) :: AHGTF(0:,0:,0:), ADER(KCKTA,KCKTB),
     &                ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
      INTEGER :: JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,IPRINT

C     Local variables
      REAL(KIND=8) :: AINT,EV,EE,EEE
      INTEGER :: ICOMPA,LVALA,MVALA,NVALA,ICOMPB,LVALB,MVALB,NVALB,
     &           IV,IU,IT

#include "onecom.h"
#include "lmns.h"

      DO ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C     **************************************************
C     ***** CALCULATE NUCLEAR ATTRACTION INTEGRALS *****
C     **************************************************
C      
      AINT = 0.0
      DO IV = 0, NVALA + NVALB
         EV = ODC(NVALA,NVALB,IV,0,0,3)
         DO IU = 0, MVALA + MVALB
            EE = ODC(MVALA,MVALB,IU,0,0,2)*EV
            DO IT = 0, LVALA + LVALB
               EEE = ODC(LVALA,LVALB,IT,0,0,1)*EE
               AINT = AINT + EEE*AHGTF(IT,IU,IV)
            END DO
         END DO
      END DO
      ADER(ICOMPA,ICOMPB) = ADER(ICOMPA,ICOMPB) + AINT

      END DO
      END DO


      RETURN
      END SUBROUTINE



      SUBROUTINE SYM_CONTRACT(AO,SO,KB,MULA,MULB,NHKTA,NHKTB,
     &     KHKTA,KHKTB,HKAB,THRESH,IMAT0,IPRINT)
C
C     Take block of distinct AO integral (derivatives) and
C     generate symmetrized contributions to SO integral
C     (derivatives) for the totally symmetric case
C                                          880407  PRT
C     Modified tuh 880819
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
      PARAMETER (LUTEMP = 48, IREPO = 0)
      REAL(KIND=8) :: AO(:,:),SO(:,:)
#include "csym1.h"
#include "symmet.h"
#include "symind.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine SYM1S',-1)
         WRITE (LUPRI, '(A, I5)') ' KB ', KB
         WRITE (LUPRI, '(A,2I5)') ' NHKTA/B ', NHKTA, NHKTB
      END IF
C
C     Loop over all irreps in molecule
C
      DO IREP = 0, MAXREP
         INDOFF = NPARSU(IREP + 1)
C
C        Loop over AOs which are of symmetry IREP in stabilizer MULA
C
         DO NA = 1, KHKTA
         IF (IBTAND(MULA,IBTXOR(IREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
            DO NB = 1,KHKTB
            IF (IBTAND(MULB,IBTXOR(IREP,ISYMAO(NHKTB,NB))).EQ.0) THEN
C
C              Weight and parity factor
C
               FAC = HKAB*PT(IBTAND(KB,IBTXOR(IREP,ISYMAO(NHKTB,NB))))
C
C              Locate SO integrals to which AO's contribute
C
               INDA = INDFA(IREP + 1,NA)
               INDB = INDFB(IREP + 1,NB)
               INDM = MAX(INDA,INDB)
               RINT = FAC*AO(NA,NB)
               IND  = INDOFF + (INDM*(INDM - 3))/2 + INDA + INDB
               SO(INDA,INDB) = SO(INDA,INDB) + RINT
            END IF
            END DO
         END IF
         END DO
      END DO
      RETURN
      END SUBROUTINE


      END MODULE
