      SUBROUTINE OPTIM_CC2_KRCC(CC_AMP,CCVEC1,CCVEC2,VEC1,VEC2,
     &                         CCTYPE,ECC,WORK,KFREE,LFREE)
*
* Optimize CC wave function
*
* Initial version, Sept 98
*
* Jeppe Olsen
*
* Summer of 99 : CCTYPE added so both CCSD and General CC can be realized
* Feb  2000    : ICCSOLVE added for DIIS etc
* June 2001    : Version with just three vectors
*
* Input :
* ========
*
* CC_AMP : Initial set of CC amplitiudes
* CCVEC6 : Imaginary Initial set of CC amplitiudes
*
* Output
* ========
*
* CC_AMP : Final set of amplitudes 
* CCVEC6 : Imaginary Final set of amplitudes
*
* Scratch
* ========
*
* VEC1, VEC2 : Two vectors for CI calculations (old CCVEC)
* CCVEC1, CCVEC2 : Two complete CC vectors
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
      CHARACTER*6 CCTYPE
#include "mxpdim.inc"
C#include "wrkspc.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cecore.inc"
#include "clunit.inc"
#include "ctcc.inc"
#include "cprnt.inc"
#include "cc_exc.inc"
#include "cands.inc"
#include "cstate.inc"
#include "cintfo.inc"
#include "cgas.inc"
#include "oper.inc"
C added
#include "cicisp.inc"
      COMMON/CNTDL_MAX/NTDL_MAX_ACT
*. Scratch 
      DIMENSION VEC1(*),VEC2(*), CCVEC1(*), CCVEC2(*)
*. Input and output
      DIMENSION CC_AMP(*),WORK(*)
* Test
      INTEGER INUMBER
*
      IDUMMY = 0
C     CALL MEMMAR(IDUMMY,IDUMMY,'MARK  ',IDUMMY,'CC_OPT')
*
      LBLK = -1
      NTDL_MAX_ACT = 0
*
      NTEST = 000 
      NTEST = MAX(NTEST,IPRCC)
      IF(NTEST.GE.5) THEN
        WRITE(6,'(A,A)') ' OPTIM_CC: Form of CC expansion ',CCTYPE
      END IF
*. Optimization of coupled cluster wavefunction defines symmetry as
      ICSM = IREFSM
      ISSM = IREFSM
      ITEX_SM = 1
*      
*. Coupled cluster flag
      ICC_EXC = 1
*
* Switch between real and complex algebra
* a multiplication factor to determine lenght of vector
      IF(IRECOM.EQ.1) THEN
        IMULTFAC = 1
      ELSE
        IMULTFAC = 2
      END IF
* ===========================================
*             Construct diagonal 
* ===========================================
*
C     IF(I_DO_NEWCCV.GE.1) THEN
C
C Inactive Fock Matrix only used in new cc vector function hence only contructed there. Lasse
C Comes in a lttle differently in CI-driven
* Obtain Inactive Fock matrix
*
      CALL FI_KRCC(WORK(KT_CC),WORK(KFI),ECC,1,ENVIRO,ISPINFREE)
C Imaginary part of inactive fock-matrix, At the moment only real.
*      CALL FI_REL(WORK(KT_CC+N1ELINT+N2ELINT),WORK(KFI+NTOOB ** 2)
*     &                 ECCCOM,1,ENVIRO,ISPINFREE)
*C      CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
*. New CC vector function uses particle hole form so 
*. use inactive Fock matrix and modified core-energy 
* May need to be changed. Lasse (look if these quantities are present)
        ECORE_H = ECORE_AB
        ECORE = ECORE_ORIG + ECORE_H
        WRITE(6,*) ' Updated core energy ',ECORE 
C Changed to a copy to KT_CC to save memory
Cori        CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
CNEW        CALL COPVEC(WORK(KFI),WORK(KT_CC),NINT1)
C     ELSE
C       ECORE_H = 0.0D0
C       ECORE_ORIG = 0.D0
C     END IF
C      print*,'real diagonal'
C      print*,(WORK(KT_CC+II),II=0,N1ELINT+N2ELINT-1)
C      print*,'imaginay diagonal'
C      print*,(WORK(KT_CC+II),II=N1ELINT+N2ELINT,2*(N1ELINT+N2ELINT))
*. And the diagonal 
      ZERO = 0.0D0
      CALL SETVEC(CCVEC1,ZERO,IMULTFAC*N_CC_AMP)
C
      IF(I_DO_NEWCCV.EQ.0) THEN
      CALL GENCC_F_DIAG_M_KRCC(WORK(KLSOBEX_CC),NSPOBEX_TP_CC,CCVEC1,
     &      1,WORK(KVEC1P),WORK(KVEC2P),MX_ST_TSOSO_MX,
     &      MX_ST_TSOSO_BLK_MX,N1ELINT,N2ELINT,WORK(KT_CC),
     &      WORK,KFREE,LFREE)
      ELSE
      CALL GENCC_F_DIAG_M_KRCC(WORK(KLSOBEX_CC),NSPOBEX_TP_CC,CCVEC1,
     &      1,KDUM,KDUM,MX_ST_TSOSO_MX,
     &      MX_ST_TSOSO_BLK_MX,N1ELINT,N2ELINT,WORK(KT_CC),
     &      WORK,KFREE,LFREE)
      END IF
*
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' OPTIM_CC : Diagonal '
        CALL WRTMAT(CCVEC1,1,2*N_CC_AMP,1,2*N_CC_AMP)
      END IF
*******
* Not in use. Lasse
*******
C     IF(I_DO_MASK_CC.EQ.1) THEN
*. Mask the excitations not arising from Mask det
C       VALUE = 1.0D+15
C       CALL MASK_CCVEC(WORK(KLSOBEX),NSPOBEX_TP_CC,CCVEC1,1,
C    &       MASK_SD,MSK_AEL,MSK_BEL,VALUE,MX_ST_TSOSO_BLK_MX)
C     END IF
C?    STOP ' Enforced stop after diagonal '
*. Save Diagonal on Disc 
*  Should probably also save the imaginary part, needed later
*  so just here a null vector.
      CALL VEC_TO_DISC(CCVEC1,IMULTFAC*N_CC_AMP,1,LBLK,LUDIA)
C reset CCVEC1
      ZERO = 0.0D0
      CALL SETVEC(CCVEC1,ZERO,IMULTFAC*N_CC_AMP)
C      print*,'n_cc_amp',n_cc_amp
C setup vector structure
C XISPSM(MXPCSM,MXPICI)
C
C
*
      MAX_DIIS_VEC=MXCIV_ORIG
      NVEC = 0
      NPVEC = -1
      IM_CONVERGED = 0
*. (Maxit is obtained from CRUN)
      DO ITER = 1, MAXIT
        IITER = ITER
        NTEST = 2
        IF(NTEST.GE.2) THEN
           WRITE(6,*)
           WRITE(6,*) ' =================================='
           WRITE(6,*) ' Information from iteration ', ITER
           WRITE(6,*) ' =================================='
           WRITE(6,*)
        END IF
*. Calculate CC_vector function for current set of parameters, 
*. result in CCVEC1
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Current CC amplitudes '
          CALL WRT_CC_VEC2(CC_AMP,LU,CCTYPE)
        END IF
C**************************************************************
C NOTE perhaps CCVEC6 should be added, Done
C ADD ECC1_IMG ?
C       CALL CC_VEC_FNC_REL(CC_AMP,CCVEC1,ECC1,VEC1,VEC2,
C    &                      CCTYPE,CCVEC2)
        call rmemchk               
C**************************************************************
*. Save on LUHC
C Save real and imagenary if complex algebra otherwise just real
        CALL VEC_TO_DISC(CCVEC1,IMULTFAC*N_CC_AMP,1,LBLK,LUHC)
C FIX FOR ECORE
        IF(ITER.EQ.1) THEN
          CALL GET_ECORE(ECORE2,I_DO_NEWCCV)
          print*,'ECC1',ECC1
        END IF
        IF(I_DO_NEWCCV.EQ.1) THEN
          ECC = ECC1 + ECORE2
        ELSE
          ECC = ECC1 + ECORE2
*          ECC = ECC1
        END IF
C        print*,'ECC,ECC1,ECORE2,iter',ECC,ECC1,ECORE2,iter
C        print*,'ECORE in optim',ECORE
C Norm of real + imagenary
        IDEBUG = 0
        IF(IDEBUG.GE.1) THEN
        print*
        do i=1,N_CC_AMP+1
          print*,'CCVEC1',CCVEC1(i)
        end do
        print*
        do i=1,N_CC_AMP
          print*,'CC_AMP',CC_AMP(i)
        end do
        print*
        do i=1,N_CC_AMP
          print*,'VEC1',VEC1(i)
        end do
        print*
        do i=1,N_CC_AMP
          print*,'VEC2',VEC2(i)
        end do
        print*
        do i=1,N_CC_AMP
          print*,'CCVEC2',CCVEC2(i)
        end do


        END IF
          XNORM = SQRT(INPROD(CCVEC1,CCVEC1,IMULTFAC*N_CC_AMP))
C        print*,'after XNORM',XNORM
        call rmemchk
****************************************
****************************************
**  Convergence threshold
**  Convergence threshold
****************************************
****************************************
Cori        CCCONV = 1.0D-6
CDEBUG        CCCONV = 1.0D-13
        CCCONV = 1.0D-10
****************************************
****************************************
        IF(XNORM.LE.CCCONV) IM_CONVERGED = 1
*
        IF(NTEST.GE.2) WRITE(6,'(A,I5,2F25.12)') 
     &  ' It, Energy, Vecnorm',ITER,ECC,XNORM
        IF(IM_CONVERGED.EQ.1) GOTO 1001
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' vector function : '
          CALL WRT_CC_VEC2(CCVEC1,LU,CCTYPE)
* (LU is not active)
        END IF
        IF(XNORM.GE.50) STOP 'NOT CONVERGING'
*
* =========================
* Obtain  correction vector 
* =========================
*
* Two methods p.t (Based on I_DO_ISBSPJA):
* 1 : Simple perturbation correction based in diagonal approximation to 
*     Jacobian
* 2 : Use Diagonal + the current set of CC vector functions to obtain 
*     the improved guess 
* Hardwirering
        I_DO_SBSPJA = 1
        IF(I_DO_SBSPJA.EQ.0) THEN
*
* Simple perturbation update
*. find change of CC coefficients as - Scaled(diag)-1 * CC_VEC_FNC on CCVEC1
*
          ZERO = 0.0D0
C              VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
C Both real and imagenary
C          print*,'before VEC_FROM_DISC'
          CALL VEC_FROM_DISC(CCVEC2,IMULTFAC*N_CC_AMP,1,LBLK,LUDIA)
C          print*,'after VEC_FROM_DISC'
          call rmemchk
C          print*,'after rmemchk'
C          do i=1,IMULTFAC*N_CC_AMP
C            print*,'CCVEC1,CCVEC2',CCVEC1(i),CCVEC2(i),i
C          end do
          CALL DIAVC2(CCVEC1,CCVEC1,CCVEC2,ZERO,IMULTFAC*N_CC_AMP)
C          print*,'after DIAVC2'
C          do i=1,IMULTFAC*N_CC_AMP
C            print*,'CCVEC1',CCVEC1(i),CCVEC2(i),i
C          end do
          call rmemchk
C          stop
C          print*,'after 2'
          ONEM = -1.0D0
          CALL SCALVE(CCVEC1,ONEM,IMULTFAC*N_CC_AMP)
        ELSE
*
* Use diagonal + subspace approximation 
*
* Add current correction vector and the corresponding Jacobian 
* vector to disc.
*
* The subspace approach is based upon assuming that the jacobian 
* is nearly constant so we can obtain information about from the 
* difference between two vector functions
*
* CC(X_{i+1}) - CC(X_{i}) = J(X_{i+1}) - CC(X_{i})
*
* So the approach is based upon a sequence of changes of 
* CC coefficients and CC vector functions
*
*. Place last set of differences of coefs and vector functions
*. on disc
          CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUSC41)
          IF(ITER.GT.1) THEN
*. Obtain X_{i+1} - X_{i} and save as vector ITER-1 in FILE LU_CCVECT
*. X_{i} is in LU_CCVECL
            CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUSC1)
            CALL SKPVCD(LU_CCVECT,ITER-2,CCVEC1,1,LBLK)
            CALL REWINO(LUSC1)
            CALL REWINO(LU_CCVECL)
            ONE = 1.0D0
            ONEM = -1.0D0
C            CALL VECSMD(CCVEC1,CCVEC2,ONE,ONEM,
C     &                  LUSC1,LU_CCVECL,LU_CCVECT,0,LBLK)
            CALL VECSMD(CC_AMP,CCVEC1,ONE,ONEM,
     &                  LUSC1,LU_CCVECL,LU_CCVECT,0,LBLK)
C          CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LU_CCVECL)
*. Obtain CCVF(X_{i+1})-CCVF(X_{i}) as vector ITER-1 in FILE LU_CCVECF
*. Previous CC vector function resides on LU_CCVECFL
            CALL SKPVCD(LU_CCVECF,ITER-2,CCVEC1,1,LBLK)
            CALL REWINO(LUHC)
            CALL REWINO(LU_CCVECFL)
            ONE = 1.0D0
            ONEM = -1.0D0
C            CALL VECSMD(CCVEC1,CCVEC2,ONE,ONEM,
C     &                  LUHC,LU_CCVECFL,LU_CCVECF,0,LBLK)
            CALL VECSMD(CCVEC1,CC_AMP,ONE,ONEM,
     &                  LUHC,LU_CCVECFL,LU_CCVECF,0,LBLK)
          END IF
*. The previous amplitudes and vector function can now 
*. be overwritten
*. X_{i+1}
C          IF(ITER.EQ.1) THEN
          CALL REWINO(LUSC41)
          CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUSC41)
          CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LU_CCVECL)
C          END IF
*. CCvec_{i+1}
C              COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
          CALL COPVCD(LUHC,LU_CCVECFL,CCVEC1,1,LBLK)
*. Multiply current CC vector function with approximate Jacobin 
*. to obtain new step
          NSBSPC_VEC = ITER-1
          MAXVEC = MIN(15,MAXIT)
C          print*,'NSBSPC_VEC ',NSBSPC_VEC
          CALL APRJAC_TV(NSBSPC_VEC,LU_CCVECFL,LUSC1,LU_CCVECT,
C     &                   LU_CCVECF,LUDIA,CCVEC1,CCVEC2,
     &                   LU_CCVECF,LUDIA,CCVEC1,CC_AMP,
     &                   WORK(KSBSPJA),N_CC_AMP,LUSC2,LUSC3,
     &                   MAXVEC)
*. The new correction vector is now residing in LUSC1,
*. Fetch and multiply with -1
          CALL VEC_FROM_DISC(CCVEC1,N_CC_AMP,1,LBLK,LUSC1)
          ONEM = -1.D0
          CALL SCALVE(CCVEC1,ONEM,N_CC_AMP)
          CALL REWINO(LUSC41)
          CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUSC41)
C*
C* Use diagonal + subspace approximation 
C*
C* Add current correction vector and the corresponding Jacobian 
C* vector to disc.
C*
C* The subspace approach is based upon assuming that the jacobian 
C* is nearly constant so we can obtain information about from the 
C* difference between two vector functions
C*
C* CC(X_{i+1}) - CC(X_{i}) = J(X_{i+1}) - CC(X_{i})
C*
C* So the approach is based upon a sequence of changes of 
C* CC coefficients and CC vector functions
C*
C*. Place last set of differences of coefs and vector functions
C*. on disc
C          print*,'iter',iter
C          IF(ITER.GT.1) THEN
C*. Obtain X_{i+1} - X_{i} and save as vector ITER-1 in FILE LU_CCVECT
C*. X_{i} is in LU_CCVECL
CC real
C            CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUSC1)
CC imagenary
C            CALL VEC_TO_DISC(CC_AMP(N_CC_AMP+1),N_CC_AMP,0,LBLK,LUSC1)
CC real part
C            CALL SKPVCD(LU_CCVECT,2*ITER-4,CCVEC1,1,LBLK)
C            CALL REWINO(LUSC1)
C            CALL REWINO(LU_CCVECL)
C            ONE = 1.0D0
C            ONEM = -1.0D0
C            CALL VECSMD(CCVEC1,CCVEC2,ONE,ONEM,
C     &                  LUSC1,LU_CCVECL,LU_CCVECT,0,LBLK)
CC imagenary part
CC input vector just dummys
CC No need to skip or rewind since files is in right place
CC            CALL SKPVCD(LU_CCVECT,2*ITER-5,CCVEC1,1,LBLK)
CC no rewind since they now should be at right place
CC            CALL REWINO(LUSC1)
CC            CALL REWINO(LU_CCVECL)
C            CALL VECSMD(CCVEC1,CCVEC2,ONE,ONEM,
C     &                  LUSC1,LU_CCVECL,LU_CCVECT,0,LBLK)
C*. Obtain CCVF(X_{i+1})-CCVF(X_{i}) as vector ITER-1 in FILE LU_CCVECF
C*. Previous CC vector function resides on LU_CCVECFL
C real
C            CALL SKPVCD(LU_CCVECF,2*ITER-4,CCVEC1,1,LBLK)
C            CALL REWINO(LUHC)
C            CALL REWINO(LU_CCVECFL)
C            ONE = 1.0D0
C            ONEM = -1.0D0
C            CALL VECSMD(CCVEC1,CCVEC2,ONE,ONEM,
C     &                  LUHC,LU_CCVECFL,LU_CCVECF,0,LBLK)
CC imagenary
CC Same here
CC            CALL SKPVCD(LU_CCVECF,2*ITER-5,CCVEC1,1,LBLK)
CC            CALL REWINO(LUHC)
CC            CALL REWINO(LU_CCVECFL)
C            CALL VECSMD(CCVEC1,CCVEC2,ONE,ONEM,
C     &                  LUHC,LU_CCVECFL,LU_CCVECF,0,LBLK)
C          END IF
C*. The previous amplitudes and vector function can now 
C*. be overwritten
C*. X_{i+1}
CC REAL
C          CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LU_CCVECL)
CC IMAGENARY
C          CALL VEC_TO_DISC(CC_AMP(N_CC_AMP+1),N_CC_AMP,0,LBLK,LU_CCVECL)
C*. CCvec_{i+1}
CC              COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C          CALL COPVCDC(LUHC,LU_CCVECFL,CCVEC1,1,LBLK)
C*. Multiply current CC vector function with approximate Jacobin 
C*. to obtain new step
C          NSBSPC_VEC = ITER-1
C          MAXVEC = MIN(15,MAXIT)
C          CALL APRJAC_TV_REL(NSBSPC_VEC,LU_CCVECFL,LUSC1,LU_CCVECT,
C     &                   LU_CCVECF,LUDIA,CCVEC1,CCVEC2,
C     &                   WORK(KSBSPJA),N_CC_AMP,LUSC2,LUSC3,
C     &                   MAXVEC)
C*. The new correction vector is now residing in LUSC1,
C*. Fetch and multiply with -1
CC REAL
C          CALL VEC_FROM_DISC(CCVEC1,N_CC_AMP,1,LBLK,LUSC1)
CC IMAGENARY
C          CALL VEC_FROM_DISC(CCVEC1(N_CC_AMP+1),N_CC_AMP,0,LBLK,LUSC1)
C          ONEM = -1.D0
C          CALL SCALVE(CCVEC1,ONEM,2*N_CC_AMP)
        END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Correction  vector '
          CALL WRT_CC_VEC2(CCVEC1,LU,CCTYPE)
* (LU is not active)
        END IF
*
* Hardwirering
        ICCSOLVE = 1
        IF(ICCSOLVE.EQ.1) THEN
*. Obtain new vector by adding correction to previous guess
          ONE = 1.0D0
* REAL AND IMG STORED AFTER EACH OTHER so just 2*N_CC_AMP
          CALL VECSUM(CC_AMP,CC_AMP,CCVEC1,ONE,ONE,IMULTFAC*N_CC_AMP)
        ELSE IF(ICCSOLVE.EQ.2) THEN 
*. DIIS
          IF(ITER.EQ.1) THEN
*. Two possibilities : 1 Norm of CC_AMP = 0 => Initial it
*                      2 Norm of CC_AMP.gt.0 => Some kind of restart
           XINI = INPROD(CC_AMP,CC_AMP,N_CC_AMP)
           IF(XINI.GT.0.0D0) THEN
*. CC_AMP is a nontrivial vector, save as first vector
            CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LU_CCVEC)
            NVEC = 1
           END IF
           IRESET = 0
          END IF
*         ^ End of check in first iteration
*. The DIIS routines works with vectors on DISC, 
*. copy current correction vector to LUSC3
          CALL VEC_TO_DISC(CCVEC1,N_CC_AMP,1,LBLK,LUSC3)
*. Add current vector to DIIS subspace, reset if req. 
          CALL ADDVEC_DIIS_SUBSPC(MAX_DIIS_VEC,NVEC,1,LU_CCVEC,
     &         LUSC3,LUSC1,LUSC2,CCVEC1,CCVEC2,NVECO,NVECD,
     &         IRESET)
*. And call DIIS accelerator to obtained improved guess on LUHC 
          IACTCOR = 1
*. Not pt taking advatage of calculated part of B matrix
          NVECP = 0
          NVEC = NVECO
          CALL DIIS_ACC(LU_CCVEC,NVEC,NVECP,WORK(KSCDIIS),
     &         WORK(KBMAT),LUHC,CCVEC1,CCVEC2,LUSC1,
     &         IACTCOR,IRESET)
*. Enforced reset ?
          I_DO_RESET = 0
          IF(I_DO_RESET.EQ.1.AND.NVEC.EQ.MAX_DIIS_VEC) THEN
            IRESET = 1
          END IF
*. And retrieve improved vector
C              VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
          CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUHC)
        END IF
*       ^ End of PERT/DIIS switch 
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Updated amplitudes '
          CALL WRT_CC_VEC2(CC_AMP,LU,CCTYPE)
* (LU is not active)
        END IF
*. Dump to LU_CCAMP
      OPEN(Unit=LU_CCAMP,File='CCAMP',Status='UNKNOWN',
     &     Form='UNFORMATTED')
        CALL REWINO(LU_CCAMP)
        I_FORMATTED = 0
C        print*,'2*N_CC_AMP ',2*N_CC_AMP
        IF(I_FORMATTED.EQ.1) THEN
          WRITE(LU_CCAMP,'(I9)') N_CC_AMP
          DO I = 1, IMULTFAC*N_CC_AMP
            WRITE(LU_CCAMP,'(E25.15)') CC_AMP(I)
          END DO
        ELSE
          WRITE(LU_CCAMP) N_CC_AMP
          WRITE(LU_CCAMP) (CC_AMP(I),I=1,IMULTFAC*N_CC_AMP)
        END IF
        CALL REWINO(LU_CCAMP)
      CLOSE(LU_CCAMP)
*
C      stop
      END DO
*     ^ End of loop over iterations
 1001 CONTINUE
*.    ^ Anchor for going out with convergence criterion fulfilled
      IF(IM_CONVERGED.EQ.0) THEN
        WRITE(6,'(A,I4,A)') 
     &  ' Convergence not obtained in ', IITER, ' iterations'
      ELSE
        WRITE(6,'(A,I4,A)') 
     &  ' Convergence obtained in ', IITER, ' iterations'
        print*,'Convergerged energy',ECC
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Final set of CC_amplitudes '
        CALL WRT_CC_VEC_REL(CC_AMP,6)
      END IF
*
      WRITE(6,*) ' Largest T(D,L) block in use = ', NTDL_MAX_ACT
C     CALL MEMMAR(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'CC_OPT')
      RETURN
      END





*
      SUBROUTINE GET_ECC(ECC)
C FOR DIRTY FIX
#include "implicit.inc"
#include "dcbdhf.h"
      ECC=DHFERG
      print*,'ECC in GET_ECC',ECC
      RETURN
      END







*
      SUBROUTINE GET_ECORE(ECORE2,I_DO_NEWCCV)
C FOR DIRTY FIX
#include "implicit.inc"
C#include "dcbdhf.h"
C#include "dcbgen.h"
C#include "frame.h"
#include "cecore.inc"
C      ECORE2 = DHFERG - ELERGY + ECORE
      LOGICAL LOGI
      CHARACTER TEXT*74
      DIMENSION IDIM(3,2)
      IF(I_DO_NEWCCV.EQ.0) THEN
      ECORE2 = ECORE 
      ELSE
      ECORE2 = E_TEMP
      END IF
      print*,'ECORE2',ECORE2
      RETURN
      END
*
      SUBROUTINE APRJAC_TV(NVEC,LUIN,LUOUT,LUVEC,LUJVEC,
     &           LUJDIA,VEC1,VEC2,SCR,N_CC_AMP,LUSCR,LUSCR2,
     &           MAXVEC)
*
* An approximate Jacobian is given in the form of
* a diagonal approximation, JACDIA, and
* NVEC vectors(in LUVEC)  and Jacobian times these vectors(LUJVEC).
*
* Find Inverse approximate Jacobian times the vector in luin to
* obtain vector in luout.
* Largest number vectors in subspave of approximate in Jacobian is MAXVEC
*
* Jeppe Olsen, March 2000
*
#include "implicit.inc"
*. Local vectors
      DIMENSION VEC1(*),VEC2(*)
*. Scratch space : Should atleast be length :
      DIMENSION SCR(*)
*.
      REAL*8 INPRDD
*
      NTEST = 000
      LBLK = -1
*.
*. Obtain overlap of input vectors
*
      KLS = 1
      KLFREE = KLS + NVEC**2
*
      KLS2 = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLEIGVEC = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLU = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLSCR = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLVEC = KLFREE
      KLFREE = KLFREE + NVEC
*
      KLVEC2 = KLFREE
      KLFREE = KLFREE + NVEC
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'APRJA_TV reporting to work '
        WRITE(6,*) '==========================='
        WRITE(6,*) ' LUVEC, LUJVEC ',LUVEC,LUJVEC
        WRITE(6,*) ' NVEC, MAXVEC = ', NVEC,MAXVEC
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Vectors on file LUVEC '
        CALL REWINO(LUVEC)
        DO IVEC = 1, NVEC
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
          CALL WRTVCD(VEC1,LUVEC,0,LBLK)
        END DO
        WRITE(6,*) ' Vectors in file LUJVEC '
        CALL REWINO(LUJVEC)
        DO IVEC = 1, NVEC
          CALL WRTVCD(VEC1,LUJVEC,0,LBLK)
        END DO
      END IF
*
      CALL REWINO(LUVEC)
      DO JVEC = 1, NVEC
        CALL REWINO(LUSCR)
        CALL SKPVCD(LUVEC,JVEC-1,VEC1,1,LBLK)
        CALL COPVCD(LUVEC,LUSCR,VEC1,0,LBLK)
C  INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
        SJJ = INPRDD(VEC1,VEC2,LUSCR,LUSCR,1,LBLK)
        JJ = JVEC*(JVEC+1)/2
        SCR(KLS-1+JJ) = SJJ
        DO IVEC = JVEC+1,NVEC
          IJ = IVEC*(IVEC-1)/2+JVEC
          CALL REWINO(LUSCR)
          SIJ = INPRDD(VEC1,VEC2,LUSCR,LUVEC,0,LBLK)
          SCR(KLS-1+IJ) = SIJ
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Overlap matrix '
        CALL PRSYM(SCR(KLS),NVEC)
      END IF
*
*
*. Check for linear dependencies by diagonalizing
*. Eliminate vectors ( start with the oldest) to obtain
*. a non-singular basis
*
*
*. Allowed ratio between largest and smallest eigenvalue
      XMXRATIO = 1.0D10
*. for the case where the are no vectors in subspace
      NVEC_EFF = 0
*
      DO NELI = 0, NVEC-1
       I1VEC = NELI + 1
       NVEC_EFF = NVEC - NELI
       LEN = NVEC_EFF*(NVEC_EFF+1)/2
*. Obtain overlap of last NVEC_EFF vectors
       DO I = NELI+1,NVEC
         DO J = NELI+1,I
          IJ = (I-NELI)*(I-NELI-1)/2 + J - NELI
          SCR(KLS2-1+IJ) = SCR(KLS-1+I*(I-1)/2  + J)
         END DO
       END DO
*. Diagonalize
C      CALL EIGEN(A,R,N,MV,MFKR)
       CALL EIGEN(SCR(KLS2),SCR(KLEIGVEC),NVEC_EFF,0,1)
*. Extract eigenvalues
       CALL COPDIA(SCR(KLS2),SCR(KLS2),NVEC_EFF,1)
*. check eigenvalue ratio and explicit singularty
       I_AM_FINE = 1
       IF(SCR(KLS2).LE.0) THEN
         I_AM_FINE = 0
       ELSE
         RATIO  = SCR(KLS2-1+NVEC_EFF)/SCR(KLS2-1+1)
         IF(RATIO.GT.XMXRATIO) I_AM_FINE = 0
       END IF
       IF(I_AM_FINE.EQ.1) GOTO 1001
      END DO
 1001 CONTINUE
C?    WRITE(6,*) ' Eigenvalues of nonsingular Metric  : '
C?    CALL WRTMAT(SCR(KLS2),1,NVEC_EFF,1,NVEC_EFF)
*     ^ End of loop over number of vectors to be truncated
*. We now know that the last NVEC_EFF vectors spans an
*. orthonormal basis
*
      NVEC_SUB = MIN(MAXVEC,NVEC_EFF)
      NVEC_SKIP = NVEC - NVEC_SUB
      I1VEC = NVEC-NVEC_SUB + 1
      IF(NVEC_SUB.EQ.0) GOTO 1002
*
C?    WRITE(6,*) ' APRJAC_TV : nvec_sub, nvec_skip ',
C?   &             NVEC_SUB, NVEC_SKIP
      CALL COP_SYMMAT(SCR(KLS),SCR(KLS2),I1VEC,NVEC_SUB)
*. Metric resides now in SCR(KLS2) in packed form
*
*. Obtain inverse of metrix
*
*.
      CALL TRIPAK(SCR(KLS),SCR(KLS2),2,NVEC_SUB,NVEC_SUB,1.0D0)
C          INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
C trying to switch to 2
      CALL INVERT_BY_DIAG2(SCR(KLS),SCR(KLEIGVEC),SCR(KLSCR),
     &                    SCR(KLVEC),NVEC_SUB)
*. And pack it
      CALL TRIPAK(SCR(KLEIGVEC),SCR(KLS),1,NVEC_SUB,NVEC_SUB,1.0D0)
      LEN = NVEC_SUB*(NVEC_SUB+1)/2
      CALL COPVEC(SCR(KLEIGVEC),SCR(KLS),LEN)
*. Inverse matrix resides now in SCR(KLS) in packed form

*
* We will use vectors I1VEC to NVEC to span the
* subspace used for
*
* The projection to the subspace is
*
* P = sum_{ij} x_i s_{ij}-1 x^t_j
*
* The approximate Jacobian is thus (J0 is diagonal part of Jacobian)
*
* J_apr = J0 (1-P) + Jex P
*       = J0 + (Jex-J0) P
*       = J0 + sum_ij (s_i - s0_i) S_{ij}-1 x^t_j (s_i = Jex x_i,s0_i = J0 x_i)
*       = J0  + R S^{-1} T^t
*
*
* where R is the collection of NVEC_SUB column vectors (s_i - s0_i)
* and T is the collection of NVEC_SUB column vectors x_i.
*
* the Sherman-Morrison formula gives the inverse of a rank-n updated matrix
*
* A' = A + R X T^t
*
* as
*
*   A'^{-1} = A^{-1} -  A^{-1}R U^{-1} T^t A^{-1}
*
* U = X^{-1} + T^t A^{-1} R
*
*
* In our case the U matrix become
* U_ij =  S - T^t_i A^{-1} R_j =
*          x_i^t J_0^{-1} J x_j
*
      CALL SKPVCD(LUJVEC,NVEC_SKIP,VEC1,1,LBLK)
      DO J = 1,NVEC_SUB
*. Diagonal times Jx_j
        ZERO = 0.0D0
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
        CALL REWINO(LUJDIA)
        CALL REWINO(LUSCR)
        CALL DMTVCD(VEC1,VEC2,LUJDIA,LUJVEC,LUSCR,ZERO,0,1,LBLK)
        CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
        DO I = 1, NVEC_SUB
*. x_i^t J0^{-1}Jx_j
          CALL REWINO(LUSCR)
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
          UIJ = INPRDD(VEC1,VEC2,LUVEC,LUSCR,0,LBLK)
          IJEXP = (J-1)*NVEC_SUB + I
          SCR(KLU-1+IJEXP) = UIJ
CERROR. and remember the overlap term
CERROR          IJSYM = I*(I-1)/2 + J
CERROR          SCR(KLU-1+ IJEXP) = SCR(KLU-1+IJEXP)
CERROR     &                      + SCR(KLS-1+IJSYM) - SCR(KLS2-1+IJSYM)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' U matrix '
        CALL WRTMAT(SCR(KLU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
      END IF
*. Inverse U matrix and keep in KLU
C          INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
C trying to switch to 2
      CALL INVERT_BY_DIAG2(SCR(KLU),SCR(KLEIGVEC),SCR(KLSCR),
     &                    SCR(KLVEC),NVEC_SUB)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverted U matrix '
        CALL WRTMAT(SCR(KLU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
      END IF
*
* And then J0{-1} R U{-1} T^t J0^{-1} VEC
*
*. 1 : J0^{-1} VEC
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,ZERO,1,1,LBLK)
*. 2 : T^t J0^{-1} VEC
      CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
      DO I = 1, NVEC_SUB
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
        CALL REWINO(LUSCR)
        SCR(KLVEC-1+I) = INPRDD(VEC1,VEC2,LUVEC,LUSCR,0,LBLK)
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T(t) J0(-1) Vecin '
        CALL WRTMAT(SCR(KLVEC),1,NVEC_SUB,1,NVEC_SUB)
      END IF
*. 3 :   VEC' = U{-1} T^t J0^{-1} VEC (Vec' is a vector in the subspace)
      CALL MATVCB(SCR(KLU),SCR(KLVEC),SCR(KLVEC2),NVEC_SUB,
     &            NVEC_SUB,0)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' U(-1) T(t) J0(-1) Vecin '
        CALL WRTMAT(SCR(KLVEC2),1,NVEC_SUB,1,NVEC_SUB)
      END IF
*  4 : J0^{-1} R U{-1} T^t J0^{-1} VEC
*      = J0^{-1} ( sum_j Jx_j vec_j)  - sum_j x_j vec_j
*
*. sum_j Jx_j vec_j on LUSCR
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '  LUSCR before call to MVCSMD'
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
      CALL REWINO(LUSCR)
      CALL REWINO(LUSCR2)
      CALL SKPVCD(LUJVEC,NVEC_SKIP,VEC1,1,LBLK)
      CALL MVCSMD(LUJVEC,SCR(KLVEC2),LUSCR,LUSCR2,VEC1,VEC2,
     &             NVEC_SUB,0,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '  (Jx)(U-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
*.  J0^{-1}  sum_j Jx_j vec_j on LUSCR2
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUSCR,LUSCR2,ZERO,1,1,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' J0(-1) (Jx)(U-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR2,1,LBLK)
      END IF
*.  sum_j x_j vec_j on LUSCR (LUOUT is used as scratch)
      CALL REWINO(LUSCR)
      CALL REWINO(LUOUT)
      CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
      CALL MVCSMD(LUVEC,SCR(KLVEC2),LUSCR,LUOUT,VEC1,VEC2,
     &             NVEC_SUB,0,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' xU(-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
*. and the synthesis in LUOUT
C  VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK
      ONE = 1.0D0
      ONEM = -1.0D0
      CALL VECSMD(VEC1,VEC2,ONE,ONEM,LUSCR2,LUSCR,LUOUT,1,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Total low rank correction'
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
      END IF
*. And the remainders
*. J0^{-1} Vecin on LUSCR
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,ZERO,1,1,LBLK)
*. (J0^{-1} -  J0^{-1} R U{-1} T^t J0^{-1} ) VECIN
      CALL VECSMD(VEC1,VEC2,ONE,ONEM,LUSCR,LUOUT,LUSCR2,1,LBLK)
*. Well the result ended up on LUSCR2, place it properly on LUOUT
C      COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
      CALL COPVCD(LUSCR2,LUOUT,VEC1,1,LBLK)
*
 1002 CONTINUE
      IF(NVEC_SUB.EQ.0) THEN
        ZERO = 0.0D0
        CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUOUT,ZERO,1,1,LBLK)
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input vector '
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUIN,1,LBLK)
        WRITE(6,*) ' Inverted approximate Jacobian times Vector '
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
      END IF
*
      RETURN
      END

