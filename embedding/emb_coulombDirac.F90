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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

! VERSION : $Revision: 13314 $
! DATE    : $Date: 2011-03-19 14:32:26 +0100 (Sat, 19 Mar 2011) $
! FILE    : emb1int.F
!
!
!  /* Deck oneprm */
      MODULE VC_GRID

      PRIVATE

!     REAL(KIND=8), save :: TOLOG, TOLS,FACTOR=1.0

      PUBLIC EMB_VC, ESP_EL, ESP_NUC
      
!     Module for calculating the coulomb potential on a grid.
!     Contains:
!     EMB_VC    The interface for users. 
!     LOOP_PRIM Loops over primitive orbitals.
!     VC_1PRIM  The routine that calls for one primitive shell VC and VC_CART
!               to calculate the shell's contribution to the coulomb potential.
!     VC        Calculates the coulomb integrals in Hermite basis
!     VC_CART   Transforms the coulomb integrals to a cartesian basis.
!     
!
!     Calls the Dirac F77 routines: GETODC in abacus/her1odc.F
!     Based on the routines: HERNAI, CINT0 in abacus/her1car.F
!     The loop over the primitives is based on abacus/her1drv.F and abacus/her1sym.F
!
!     The include files that need to be set:
!     onecom.h
!     sdepre.h
!

      CONTAINS

      SUBROUTINE EMB_VC(DMAT,points,vc,IPRINT,WORK,LWORK)
!     Interface for the user. Calles LOOP_PRIM
!     calls VC_PRIM.
!     Needed to make contact with F77.
!     JHS DMAT is in MO basis. Transformations needed.

      IMPLICIT NONE
#include "priunit.h"
#include "maxaqn.h"
#include "ccom.h"

!     PASSED variables
      real(kind=8) :: points(:,:),vc(:)
      INTEGER :: IPRINT
      REAL(KIND=8) :: DMAT(*)
      REAL(KIND=8) :: WORK(*)
      INTEGER :: LWORK

      
      IF (IPRINT .GE. 5) CALL TITLER('Output from EMB_VC','*',103)

!     Electrostatic potential: electronic contribution
      CALL ESP_EL(DMAT,points,vc,WORK,LWORK,IPRINT)

!     Electrostatic potential: nuclear contribution
      CALL ESP_NUC(points,vc,IPRINT)

      END SUBROUTINE

      SUBROUTINE ESP_EL(DMAT,points,vc,WORK,LWORK,IPRINT)
!     Driver routine that loops over the primitive orbitals and
!     calls VC_PRIM.
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
#include "ccom.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"

!     PASSED variables
      real(kind=8) :: points(:,:),vc(:)
      INTEGER :: IPRINT
      REAL(KIND=8) :: DMAT(NTBAS(0),NTBAS(0),NZ)
      REAL(KIND=8) :: WORK(*)
      INTEGER :: LWORK

!     Local Variables
      INTEGER :: ISHELA,ISHELB,ICA,ICB,NA,IREP,NB
      INTEGER :: MULTA,MULTB,NHKTAB,KAB,IDENB_TMP
      INTEGER :: ISYM
      INTEGER :: ipoint, i,j,k,l

!      REAL(KIND=8) :: CORA(3),CORB(3)
      REAL(KIND=8),ALLOCATABLE :: ADER(:,:,:)
!     Contains the coulomb integrals for one set of primitive orbitals.

!     JHS TMP STUFF FOR DEBUGGING
      REAL(KIND=8) :: gexp=0.0,vc_tmp
      REAL(KIND=8) :: dummy
      INTEGER(KIND=8) :: idummy
      REAL(KIND=8) :: TOLOG, TOLS,FACTOR=1.0

#include "cbisol.h"
#include "shells.h"
!#include "inforb.h"
#include "onecom.h"
#include "primit.h"
#include "symind.h"
#include "symmet.h"
#include "lmns.h"
#include "nuclei.h"

!     Initialization
      ISOFRA(:)=0
      TOLS=THRS**2
      TOLOG=2*LOG(THRS)

!     Loop over shells.
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
!
!       Compute symmetry integral pointers for contributions
!       from this block.  Note that at present this assumes all
!       components from a shell are included.
!
        IDENB = 0
        ISOFRB(:) = 0
        INDFA(:,:) = -10000000
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
!     JHS Later we need to multiply with the sign for the symmetry operations
          CORBX  = CENT(ISHELB,1,1)
          CORBY  = CENT(ISHELB,2,1)
          CORBZ  = CENT(ISHELB,3,1)
          KHKTAB = KHKTA*KHKTB
          KCKTAB = KCKTA*KCKTB
          MAB    = IOR(MULA,MULB)
          KAB    = IAND(MULA,MULB)
          HKAB   = FMULT(KAB)
!
          SPHRAB = SPHRA .OR. SPHRB
!     JHS Note that if SPHRAB is true we should correct the densitymatrix for this.
!
!        Compute symmetry integral pointers for contributions
!        from this block.  Note that at present this assumes all
!        components from a shell are included
!
          INDFB(:,:) = -10000000
          DO NB = 1, KHKTB
            DO IREP = 0, MAXREP
            IF (IAND(MULB,IEOR(IREP,ISYMAO(NHKTB,NB)))==0) THEN
              ISOFRB(IREP+1)    = ISOFRB(IREP+1) + 1
              INDFB (IREP+1,NB) = ISOFRB(IREP+1)
            END IF
            END DO
          END DO
!        jth - CLASS (which integral class: LL,LS,SL,SS)
!                                           0  1  1  2
!     JHS FOLLOWING STATEMENT NEEDS DEBUGGING FOR OTHER COMPILERS
         IF (ICA==1 .AND. ICB==1) THEN; INTCLASS=2
         ELSE IF (ICA==1.OR.ICB==1 ) THEN; INTCLASS=1
         ELSE; INTCLASS=0
         END IF

         IF (IPRINT>5)                                                  &
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
!
!        Initialization for nuclear attraction integrals
!
         JMAX = NHKTAB - 2
         ISTEPU = JMAX + 1
         ISTEPV = ISTEPU*ISTEPU
         NAHGTF = ISTEPU*ISTEPV

         ALLOCATE(ADER(KCKTA,KCKTB,SIZE(points,2)))
!        Contains ader for one gridpoint in SO basis.


         IF (IPRINT>5) WRITE(LUPRI,*) 'JHS SYMMETRY,SYMOP',MAXOPR
!     JHS Now loop over symmetry dependend stuff.

!     JHS Maybe IDENB needs to be corrected by -KHKTB
         IDENB_TMP=IDENB-KHKTB
!     JHS IDENB_TMP variable necessary?
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
            
!     JHS Correct the B coordinate for the symmetry operations.
!         Direct passing is not possible as CORA,CORB are used inside her1odc.F
            CORBX=CORBX*SIGNBX
            CORBY=CORBY*SIGNBY
            CORBZ=CORBZ*SIGNBZ
 
            IF (IPRINT>5) WRITE(LUPRI,*) 'JHS CORB',CORBX,CORBY,CORBZ
            
            CALL VC_SHELL(ADER,TOLOG,TOLS,points,IPRINT,WORK,LWORK)
            
!     JHS Transform to SO basis?

            IF (IPRINT>5) WRITE(LUPRI,*) 'JHS ADER(VC_SHELL)',ADER

!     JHS Start contraction with the density matrix.

!     JHS Density matrix is in AO basis. 
!     Contraction of ADER to AO basis done in CART_VC (in ODC coefficients)
!     Maybe move inside PRIM_VC
!     Similar to dirac/diraos.F:GETGAO and dirac/dirden.F:DENMAT
!
!     Contract ADER with density matrix, indices in ADER are swapped.
            IF (IPRINT>10) THEN
                 WRITE(LUPRI,*) 'DMAT Contraction.',ISHELA,ISHELB

                 WRITE(LUPRI,*) ISHELA,ISHELB,KCKTA,KCKTB,IDENA,        &
     &           IDENB_TMP,NZ,SPHRA,SPHRB            
            END IF

       DO ipoint=1,SIZE(vc)               
!     Contract with the right symmetry part of the density matrix.
!     This is given by INDFA and INDFB


!     No contraction over the non-real quaternionic part. Use index matrices of quaternions
               vc(ipoint)=vc(ipoint)+SUM(                               &
     &              DMAT(idena+1:idena+KCKTA,idenb+1:idenb+KCKTB,1)*    &
     &              ADER(1:KCKTA,1:KCKTB,ipoint))

       END DO 
           

!     JHS Transform back to the original B coordinate.
            CORBX=CORBX/SIGNBX
            CORBY=CORBY/SIGNBY
            CORBZ=CORBZ/SIGNBZ
         END DO

         DEALLOCATE(ADER)


         IDENB = IDENB + KHKTB*MULTB
       END DO
       IDENA = IDENA + KHKTA*MULTA
      END DO

      END SUBROUTINE
      
      SUBROUTINE ESP_NUC(points,vc,IPRINT)
!     Routine for calcuating electrostatic potential from nuclei in points
!     provided by array points
      USE VECTOR_FUN, ONLY: NORM
      IMPLICIT NONE
#include "priunit.h"
#include "mxcent.h"
!     PASSED variables
      real(kind=8) :: points(:,:),vc(:)
      INTEGER :: IPRINT
!     Local Variables
      INTEGER :: ipoint, i
      REAL(KIND=8) :: vc_tmp

#include "nuclei.h"

      DO ipoint=1,SIZE(vc)
        vc_tmp=0.0D0
        DO i=1,MXCENT          
          vc_tmp=vc_tmp-CHARGE(i)/(NORM(CORD(:,i)-                      &
     &          points(1:3,ipoint)))
        END DO
        vc(ipoint)=vc(ipoint)-vc_tmp
      END DO

      END SUBROUTINE


      SUBROUTINE VC_SHELL(ADER,TOLOG,TOLS,points,IPRINT,WORK,           &
     &     LWORK)
!     Calculates the contribution for one primitive orbital set.
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
!     PASSED variables
      INTEGER :: IPRINT
      real(kind=8) :: points(:,:)
      REAL(KIND=8) :: ADER(KCKTA,KCKTB,SIZE(points,2))
      REAL(KIND=8) :: WORK(*)
      INTEGER :: LWORK

      REAL(KIND=8) :: TOLOG, TOLS, factor=1.0


!     Local variables
      REAL(KIND=8),PARAMETER :: D0=0.0D0, D1=1.00D00, D3=3.00D00,       &
     &                          THIRD=D1/D3, GNUTHR=1.D-4,D1P5=1.5D0
      REAL(KIND=8)::AHGTF(0:JMAX,0:JMAX,0:JMAX)
      REAL(KIND=8),ALLOCATABLE::ODC(:,:,:,:,:,:)

      LOGICAL :: KINODC
      REAL(KIND=8) :: DIFAB(3),DIFCP(3),CORP(3),ORIGIN(3),              &
     &     CORA(3),CORB(3)
      REAL(KIND=8) :: DISTAB,CONTA,EXPA,CONTB,EXPB,EXPP,EXPPI,          &
     &                EXPABQ,SAAB, ASAAB,SAAB13, EXPAPI,EXPBPI
      INTEGER :: JMAXA,JMAXB,JMAXD,JMAXT,JMAXM,ipoint,                  &
     &           IPRIMA,IPRIMB,JPRIMA,JPRIMB
      INTEGER :: IDUMMY
#include "cbisol.h"
#include "onecom.h"
#include "ader.h"
#include "primit.h"
     
      IF (IPRINT .GE. 5) CALL TITLER('Output from VC_SHELL','*',103)

!     ALLOCATION
      JMAXD = 2
      JMAXA = NHKTA - 1
      JMAXB = NHKTB - 1
      JMAXT = JMAXA + JMAXB + JMAXD
      JMAXM = 0
      ALLOCATE(ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3))

!     Initialization
      ADER(:,:,:)=0.0
      CORA(1)=CORAX
      CORA(2)=CORAY
      CORA(3)=CORAZ
      CORB(1)=CORBX
      CORB(2)=CORBY
      CORB(3)=CORBZ

      DIFAB(:) = CORA(:) - CORB(:)
      DISTAB=NORM(DIFAB)

!
!     ****************************************
!     ***** Loop over primitive orbitals *****
!     ****************************************
!
      IF (IPRINT>10)                                                    &
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
!
!           Calculate and test square root factor
!     JHS SAAB are the contraction coefficients from the primitive orbitals to the basisfunctions.
            EXPABQ = EXPA*EXPB*DISTAB*EXPPI
            SAAB = CONTA*CONTB*EXP(-EXPABQ)
            ASAAB = ABS(SAAB)
            IF (IPRINT>20) WRITE(LUPRI,*) 'JHS INTSCREENING',           &
     &           'EXPABQ=',EXPABQ,'TOLOG=',TOLOG,                       &
     &           'ASAAB=',ASAAB,'TOLS=',TOLS
            IF (EXPABQ<TOLOG) CYCLE
            IF (ASAAB<TOLS) CYCLE
            SAAB13 = SIGN(ASAAB**THIRD,SAAB)
!
!           Calculate coordinates of product Gaussian P
!
            EXPAPI = EXPA*EXPPI
            EXPBPI = EXPB*EXPPI
            CORP(:)=EXPAPI*CORA(:) + EXPBPI*CORB(:)
!
!           *********************************************
!           ***** Overlap Distribution Coefficients *****
!           *********************************************
!
!           Expansion coefficients for undifferentiated orbitals
!
#ifdef PRG_DIRAC
!           jth: removed KINODC = .NOT.PROPTY
            KINODC = .FALSE.
#else
            KINODC = .TRUE.
#endif 
            IDUMMY = 0
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,.FALSE.,      &
     &           KINODC,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,WORK,LWORK,&
     &           CORP(1),CORP(2),CORP(3),.TRUE.,.FALSE.,ORIGIN,IDUMMY)
!
!           **********************************************
!           ***** Calculation of Hermitian integrals *****
!           **********************************************
!
!           Nuclear attraction integrals
!
            DO ipoint=1,SIZE(points,2)
              DIFCP(:) = points(1:3,ipoint) - CORP(:)
              CALL VNUC(AHGTF,JMAX,EXPP,DIFCP,FACTOR,IPRINT)
!
!           **********************************************
!           ***** Calculation of Cartesian integrals *****
!           **********************************************
!             
              CALL CART_VC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,           &
     &             ADER(1:KCKTA,1:KCKTB,ipoint:ipoint),AHGTF,IPRINT)
            END DO
         END DO
      END DO

      DEALLOCATE(ODC)

!
!     Print
!
      
      RETURN
      END SUBROUTINE
!  /* Deck hernai */




      SUBROUTINE VNUC(AHGTF,JMAX,PVAL,CP,FACINT,IPRINT)
      USE VECTOR_FUN, ONLY: NORM
      IMPLICIT NONE
!
!     This subroutine calculates the R integrals as defined by
!     McMurchie and Davidson in J. Comp. Phys. 26 (1978) 218.
!     The recursion formulas (4.6) - (4.8) are used.
!
!     Only slightly slower than the implementation HERNAI in abacus/her1car.F of TUH 
!
!     JHS 26-03-08
!     
#include "priunit.h"
#include "pi.h"
#include "maxaqn.h"

!     Arguments passed
      REAL(KIND=8) :: AHGTF(0:,0:,0:),CP(3)
      REAL(KIND=8) :: PVAL,FACINT
      INTEGER :: JMAX,IPRINT

!     Local Variables
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
!
!     *************************************
!     ***** Incomplete Gamma Function *****
!     *************************************
!
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
            RJ(0,0,IV+1,JVAL-1)= -CP(3)*RJ(0,0,IV,JVAL)+                &
     &                             IV*RJ(0,0,IV-1,JVAL)
!     RJ(0,1,v) etc. is moved outside the loop for efficiency
            RJ(0,1,IV,JVAL-1)= -CP(2)*RJ(0,0,IV,JVAL)
            DO IU=1,JMAX-JVAL-IV
               RJ(0,IU+1,IV,JVAL-1)= -CP(2)*RJ(0,IU,IV,JVAL)+           &
     &                                 IU*RJ(0,IU-1,IV,JVAL)
               RJ(1,IU,IV,JVAL-1)= -CP(1)*RJ(0,IU,IV,JVAL)
               DO IT=1,JMAX-JVAL-IU-IV
!     USE TRANSFER(IT,real_8) etc. if this is platform dependend.
                  RJ(IT+1,IU,IV,JVAL-1)= -CP(1)*RJ(IT,IU,IV,JVAL)+      &
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

!  /* Deck cint0 */
      SUBROUTINE CART_VC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,             &
     &                 ADER,AHGTF,IPRINT)
!     TUH
      IMPLICIT NONE
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
!     Variables passed to CART_VC
      REAL(KIND=8) :: AHGTF(0:,0:,0:), ADER(KCKTA,KCKTB),               &
     &                ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
      INTEGER :: JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,IPRINT

!     Local variables
      REAL(KIND=8) :: AINT,EV,EE,EEE
      INTEGER :: ICOMPA,LVALA,MVALA,NVALA,ICOMPB,LVALB,MVALB,NVALB,     &
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
!
!     **************************************************
!     ***** CALCULATE NUCLEAR ATTRACTION INTEGRALS *****
!     **************************************************
!      
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



      SUBROUTINE SYM_CONTRACT(AO,SO,KB,MULA,MULB,NHKTA,NHKTB,           &
     &     KHKTA,KHKTB,HKAB,THRESH,IMAT0,IPRINT)
!
!     Take block of distinct AO integral (derivatives) and
!     generate symmetrized contributions to SO integral
!     (derivatives) for the totally symmetric case
!                                          880407  PRT
!     Modified tuh 880819
!
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
!
!     Loop over all irreps in molecule
!
      DO IREP = 0, MAXREP
         INDOFF = NPARSU(IREP + 1)
!
!        Loop over AOs which are of symmetry IREP in stabilizer MULA
!
         DO NA = 1, KHKTA
         IF (IBTAND(MULA,IBTXOR(IREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
            DO NB = 1,KHKTB
            IF (IBTAND(MULB,IBTXOR(IREP,ISYMAO(NHKTB,NB))).EQ.0) THEN
!
!              Weight and parity factor
!
               FAC = HKAB*PT(IBTAND(KB,IBTXOR(IREP,ISYMAO(NHKTB,NB))))
!
!              Locate SO integrals to which AO's contribute
!
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
