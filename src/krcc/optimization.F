      SUBROUTINE OPTIM_CC2_KRCC(CC_AMP,CCVEC1,
     &                         ECC,WORK,KFREE,LFREE)
*
* Optimize CC wave function
*
* Initial version, Sept 98
*
* Jeppe Olsen
*
* Input :
* ========
*
* CC_AMP : Initial set of CC amplitiudes
*
* Output
* ========
*
* CC_AMP : Final set of amplitudes 
*
* Scratch
* ========
*
* CCVEC1 : Complete CC vector
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
#include "mxpdim.inc"
#include "ipoist8.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cecore.inc"
#include "clunit.inc"
#include "ctcc.inc"
#include "ctccp.inc"
#include "cands.inc"
#include "cstate.inc"
#include "cintfo.inc"
#include "cgas.inc"
C added
#include "cicisp.inc"
*. Scratch 
      DIMENSION CCVEC1(*)
*. Input and output
      DIMENSION CC_AMP(*),WORK(*)
*
      IDUMMY = 0
      I_WRITE_INTERMEDIATES = 0
C     CALL MEMMAR(IDUMMY,IDUMMY,'MARK  ',IDUMMY,'CC_OPT')
*
      LBLK = -1
*
      NTEST = 000 
*. Optimization of coupled cluster wavefunction defines symmetry as
      ICSM = IREFSM
      ISSM = IREFSM
*      
*
* Switch between real and complex algebra
* a multiplication factor to determine lenght of vector
      IF(IRECOM.EQ.1) THEN
        IMULTFAC = 1
      ELSE
        IMULTFAC = 2
      END IF
      print*,'IRECOM',IRECOM
* ===========================================
*             Construct diagonal 
* ===========================================
*
* Obtain Inactive Fock matrix
*
      CALL FI_KRCC(WORK(KT_CC),WORK(KFI),ECC,1,ISPINFREE)
C Imaginary part of inactive fock-matrix, At the moment only real.
*      CALL FI_REL(WORK(KT_CC+N1ELINT+N2ELINT),WORK(KFI+NTOOB ** 2)
*     &                 ECCCOM,1,ISPINFREE)
*. New CC vector function uses particle hole form so 
*. use inactive Fock matrix and modified core-energy 
* May need to be changed. Lasse (look if these quantities are present)
        ECORE_H = ECORE_AB
        ECORE = ECORE_ORIG + ECORE_H
        WRITE(6,*) ' Updated core energy ',ECORE 
*. And the diagonal 
      ZERO = 0.0D0
      CALL SETVEC(CCVEC1,ZERO,IMULTFAC*N_CC_AMP)
C
C Remember the 1 is the symmetry
C
      CALL GENCC_F_DIAG_M_KRCC(WORK(KLSOBEX_CC),NSPOBEX_TPE,CCVEC1,
     &      1,KDUM,KDUM,MX_ST_TSOSO_MX,
     &      MX_ST_TSOSO_BLK_MX,N1ELINT,N2ELINT,WORK(KT_CC),
     &      WORK,KFREE,LFREE)
*
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' OPTIM_CC : Diagonal '
        CALL WRTMAT(CCVEC1,1,2*N_CC_AMP,1,2*N_CC_AMP)
      END IF
*. Save Diagonal on Disc 
*  Should probably also save the imaginary part, needed later
*  so just here a null vector.
      CALL VEC_TO_DISC(CCVEC1(2),IMULTFAC*N_CC_AMP-1,1,LBLK,LUDIA)
C reset CCVEC1
      ZERO = 0.0D0
      CALL SETVEC(CCVEC1,ZERO,IMULTFAC*N_CC_AMP)
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
          CALL WRT_CC_VEC2(CC_AMP,LU,
     &                     WORK,KFREE,LFREE)
        END IF
C**************************************************************
C NOTE perhaps CCVEC6 should be added, Done
C ADD ECC1_IMG ?
        CALL CC_VEC_FNC_KRCC(CC_AMP,CCVEC1,ECC1,ITER,
     &                       WORK,KFREE,LFREE)
C       stop 'after CC_VEC_FNC_KRCC'
        CALL MEMCHK_KRCC(WORK)
C**************************************************************
*. Save on LUHC
C Save real and imagenary if complex algebra otherwise just real
C         CALL FIRST_ELEMENT_LAST(CC_AMP,IMULTFAC*N_CC_AMP)
C         CALL FIRST_ELEMENT_LAST(CCVEC1,IMULTFAC*N_CC_AMP)
        CALL VEC_TO_DISC(CCVEC1(2),IMULTFAC*N_CC_AMP-1,1,LBLK,LUHC)
C FIX FOR ECORE
        IF(ITER.EQ.1) THEN
          CALL GET_ECC(ECORE2)
        END IF
          ECC = ECC1 + ECORE2
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
        do i=1,N_CC_AMP +1
          print*,'CC_AMP',CC_AMP(i)
        end do
        END IF
C Notice the identety is now number one and not the last element!
         XNORM = SQRT(INPROD(CCVEC1(2),CCVEC1(2),IMULTFAC*N_CC_AMP))
         print*,'after XNORM',XNORM
         CALL MEMCHK_KRCC(WORK)
****************************************
****************************************
**  Convergence threshold
**  Convergence threshold
****************************************
****************************************
****************************************
****************************************
        IF(XNORM.LE.CCCONV) IM_CONVERGED = 1
*
        IF(NTEST.GE.2) WRITE(6,'(A,I5,2F25.12)') 
     &  ' It, Energy, Vecnorm',ITER,ECC,XNORM
        IF(IM_CONVERGED.EQ.1) GOTO 1001
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' vector function : '
          CALL WRT_CC_VEC2(CCVEC1,LU,
     &                     WORK,KFREE,LFREE)
* (LU is not active)
        END IF
        IF(XNORM.GE.500) STOP 'NOT CONVERGING'
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
        CALL MEMGET('REAL',KDIA,IMULTFAC*N_CC_AMP,WORK,KFREE,LFREE)
*
* Simple perturbation update
*. find change of CC coefficients as - Scaled(diag)-1 * CC_VEC_FNC on CCVEC1
*
          ZERO = 0.0D0
C              VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
C Both real and imagenary
C          print*,'before VEC_FROM_DISC'
         CALL VEC_FROM_DISC(WORK(KDIA),IMULTFAC*N_CC_AMP-1,1,LBLK,LUDIA)
           print*,'after VEC_FROM_DISC'
          CALL MEMCHK_KRCC(WORK)
           print*,'after rmemchk'
C          do i=1,IMULTFAC*N_CC_AMP
C            print*,'CCVEC1,CCVEC2',CCVEC1(i),CCVEC2(i),i
C          end do
          CALL DIAVC2(CCVEC1(2),CCVEC1(2),WORK(KDIA),
     &                ZERO,IMULTFAC*N_CC_AMP-1)
           print*,'after DIAVC2'
C          CALL DUMDUM(N_CC_AMP,CCVEC1,WORK(KDIA))
C          do i=1,IMULTFAC*N_CC_AMP
C            print*,'CCVEC1',CCVEC1(i),CCVEC2(i),i
C          end do
          CALL MEMCHK_KRCC(WORK)
           print*,'after 2'
          ONEM = -1.0D0
          CALL SCALVE(CCVEC1(2),ONEM,IMULTFAC*N_CC_AMP-1) !see if correct
        CALL MEMREL('MAPPING',WORK,KDIA,KDIA,KFREE,LFREE)
        ELSE
* Silly Lasse edition due to laziness
* Reason I have the identity first Jeppe has it last
* and I don't feel like changing his optimizer
C         CALL FIRST_ELEMENT_LAST(CC_AMP,IMULTFAC*N_CC_AMP)
C         CALL FIRST_ELEMENT_LAST(CCVEC1,IMULTFAC*N_CC_AMP)
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
          CALL VEC_TO_DISC(CC_AMP(2),N_CC_AMP-1,1,LBLK,LUSC41)
          IF(ITER.GT.1) THEN
*. Obtain X_{i+1} - X_{i} and save as vector ITER-1 in FILE LU_CCVECT
*. X_{i} is in LU_CCVECL
            CALL VEC_TO_DISC(CC_AMP(2),N_CC_AMP-1,1,LBLK,LUSC1)
            CALL SKPVCD(LU_CCVECT,ITER-2,CCVEC1(2),1,LBLK)
            CALL REWINO(LUSC1)
            CALL REWINO(LU_CCVECL)
            ONE = 1.0D0
            ONEM = -1.0D0
C            CALL VECSMD(CCVEC1,CCVEC2,ONE,ONEM,
C     &                  LUSC1,LU_CCVECL,LU_CCVECT,0,LBLK)
            CALL VECSMD(CC_AMP(2),CCVEC1(2),ONE,ONEM,
     &                  LUSC1,LU_CCVECL,LU_CCVECT,0,LBLK)
C          CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LU_CCVECL)
*. Obtain CCVF(X_{i+1})-CCVF(X_{i}) as vector ITER-1 in FILE LU_CCVECF
*. Previous CC vector function resides on LU_CCVECFL
            CALL SKPVCD(LU_CCVECF,ITER-2,CCVEC1(2),1,LBLK)
            CALL REWINO(LUHC)
            CALL REWINO(LU_CCVECFL)
            ONE = 1.0D0
            ONEM = -1.0D0
C            CALL VECSMD(CCVEC1,CCVEC2,ONE,ONEM,
C     &                  LUHC,LU_CCVECFL,LU_CCVECF,0,LBLK)
            CALL VECSMD(CCVEC1(2),CC_AMP(2),ONE,ONEM,
     &                  LUHC,LU_CCVECFL,LU_CCVECF,0,LBLK)
          END IF
*. The previous amplitudes and vector function can now 
*. be overwritten
*. X_{i+1}
C          IF(ITER.EQ.1) THEN
          CALL REWINO(LUSC41)
          CALL VEC_FROM_DISC(CC_AMP(2),N_CC_AMP-1,1,LBLK,LUSC41)
          CALL VEC_TO_DISC(CC_AMP(2),N_CC_AMP-1,1,LBLK,LU_CCVECL)
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
     &                   LU_CCVECF,LUDIA,CCVEC1(2),CC_AMP(2),
     &                   WORK(KSBSPJA),N_CC_AMP-1,LUSC2,LUSC3,
     &                   MAXVEC)
*. The new correction vector is now residing in LUSC1,
*. Fetch and multiply with -1
          CALL VEC_FROM_DISC(CCVEC1(2),N_CC_AMP-1,1,LBLK,LUSC1)
          ONEM = -1.D0
          CALL SCALVE(CCVEC1(2),ONEM,N_CC_AMP-1)
          CALL REWINO(LUSC41)
          CALL VEC_FROM_DISC(CC_AMP(2),N_CC_AMP-1,1,LBLK,LUSC41)
* Lasse addition for laziness
C         CALL LAST_ELEMENT_FIRST(CC_AMP,IMULTFAC*N_CC_AMP)
C         CALL LAST_ELEMENT_FIRST(CCVEC1,IMULTFAC*N_CC_AMP)
        END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Correction  vector '
          CALL WRT_CC_VEC2(CCVEC1,LU,
     &                     WORK,KFREE,LFREE)
* (LU is not active)
        END IF
*
* Hardwirering
        ICCSOLVE = 1
        IF(ICCSOLVE.EQ.1) THEN
*. Obtain new vector by adding correction to previous guess
          ONE = 1.0D0
* REAL AND IMG STORED AFTER EACH OTHER so just 2*N_CC_AMP
          CALL VECSUM(CC_AMP(2),CC_AMP(2),CCVEC1(2),
     &                ONE,ONE,IMULTFAC*N_CC_AMP-1)
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
          CALL WRT_CC_VEC2(CC_AMP,LU,
     &                     WORK,KFREE,LFREE)
* (LU is not active)
        END IF
*. Dump to LU_CCAMP
      OPEN(Unit=LU_CCAMP,File='CCAMP',Status='UNKNOWN',
     &     Form='UNFORMATTED')
C    &     Form='FORMATTED')
        CALL REWINO(LU_CCAMP)
        I_FORMATTED = 2
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
        WRITE(6,*) 'Convergerged energy',ECC
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Final set of CC_amplitudes '
        CALL WRT_CC_VEC2(CC_AMP,6,
     &                   WORK,KFREE,LFREE)
      END IF
*
C     WRITE(6,*) ' Largest T(D,L) block in use = ', NTDL_MAX_ACT
*. Dump to LU_CCAMP
      OPEN(Unit=LU_CCAMP,File='CCAMP',Status='UNKNOWN',
     &     Form='UNFORMATTED')
C    &     Form='FORMATTED')
        CALL REWINO(LU_CCAMP)
        I_FORMATTED = 2
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
C     CALL MEMMAR(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'CC_OPT')
      RETURN
      END



      SUBROUTINE GET_ECC(ECC)
C FOR DIRTY FIX
#include "implicit.inc"
#include "dcbdhf.h"
#include "cecore.inc"
      ECC=DHFERG
      print*,'ECC in GET_ECC',ECC
      print*,'ECORE,E_TEMP',ECORE,E_TEMP
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
      CALL INVERT_BY_DIAG3(SCR(KLS),SCR(KLEIGVEC),SCR(KLSCR),
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
        CALL REWINO(LUJDIA)
        CALL REWINO(LUSCR)
        CALL DMTVCD(VEC1,VEC2,LUJDIA,LUJVEC,LUSCR,ZERO,0,1,LBLK)
        CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
        DO I = 1, NVEC_SUB
*. x_i^t J0^{-1}Jx_j
          CALL REWINO(LUSCR)
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
      CALL INVERT_BY_DIAG3(SCR(KLU),SCR(KLEIGVEC),SCR(KLSCR),
     &                    SCR(KLVEC),NVEC_SUB)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverted U matrix '
        CALL WRTMAT(SCR(KLU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
      END IF
*
* And then J0{-1} R U{-1} T^t J0^{-1} VEC
*
*. 1 : J0^{-1} VEC
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,ZERO,1,1,LBLK)
*. 2 : T^t J0^{-1} VEC
      CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
      DO I = 1, NVEC_SUB
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
*
      SUBROUTINE ADDVEC_DIIS_SUBSPC(MAXVEC,NVECIN,NVECADD,
     &           LUVEC,LUADD,LUSCR1,LUSCR2,VEC1,VEC2,NVECOUT,NVECDEL,
     &           IRESET)
*
* Add NVECADD vectors to DIIS subspace. If the number of vectors in subspace 
* before addition is MAXVEC, remove also the oldest correction
*
* LUVEC : Initial approximation to solution vector and 
*         NVECIN - 1 correction correction
*     
* IF IRESET = 1, then the LUVEC file is reset : current approximation
* is constructed on LUVEC and the NVECADD vectors are added
* Jeppe Olsen, Feb. 2000
*     
#include "implicit.inc"
      DIMENSION VEC1(*),VEC2(*)
*. Local scratch : Assuming max dim of subspace is 100
      PARAMETER (MXP_DIM_SUBSPACE = 100)
      DIMENSION SCR(MXP_DIM_SUBSPACE)
*
      IF(MAXVEC.GT.MXP_DIM_SUBSPACE) THEN
        WRITE(6,*) ' Potential problem for ADDVEC_DIIS_SUBSPC '
        WRITE(6,*) ' Max dim of subspace larger than program limit'
        WRITE(6,*) ' MAXVEC, MXP_DIM_SUBSPACE = ',
     &               MAXVEC, MXP_DIM_SUBSPACE
        WRITE(6,*) ' Increase  MXP_DIM_SUBSPACE or decrease MAXVEC'
        STOP       ' Increase  MXP_DIM_SUBSPACE or decrease MAXVEC'
      END IF
*
      LBLK = -1
*
      IF(IRESET.EQ.0) THEN
        NVECOUT = MIN(MAXVEC,NVECIN+NVECADD)
        NVECDEL = MAX(0,NVECIN+NVECADD-MAXVEC)
      ELSE
        NVECOUT = 1 + NVECADD
        NVECDEL = NVECIN - 1
      END IF
C?    WRITE(6,*) ' NVECIN, NVECADD,NVECOUT, NVECDEL',
C?   &             NVECIN, NVECADD,NVECOUT, NVECDEL
      IF(NVECDEL.EQ.0) THEN
*. Just add vectors
       CALL SKPVCD(LUVEC,NVECIN,VEC1,1,LBLK)
       CALL REWINO(LUADD)
       DO IVEC = 1, NVECADD
         CALL COPVCD(LUADD,LUVEC,VEC1,0,LBLK)
       END DO
      ELSE IF( NVECDEL.GT.0 ) THEN
*. Reset and add 
*. New first vector (initial approximation) 
*  should now by sum of first NVECDEL+1 vectors 
*. on LUVEC  
       NVECDELP1 = NVECDEL + 1
       CALL REWINO(LUVEC)
       ONE = 1.0D0
       CALL SETVEC(SCR,ONE,NVECDELP1)
C           MVCSMD(LUIN, FAC,LUOUT, LUSCR, VEC1,VEC2,NVEC,IREW,LBL
       CALL MVCSMD(LUVEC,SCR,LUSCR1,LUSCR2,VEC1,VEC2,NVECDELP1,1,LBLK)
*. Copy to LUSCR1 the vectors that should be kept after new vector 
       CALL SKPVCD(LUVEC,NVECDELP1,VEC1,1,LBLK)
       DO IVEC = NVECDELP1+1,NVECIN
         CALL COPVCD(LUVEC,LUSCR1,VEC1,0,LBLK)
       END DO
       CALL REWINO(LUSCR1)
       CALL REWINO(LUVEC)
       DO IVEC = 1, NVECIN-NVECDEL
         CALL COPVCD(LUSCR1,LUVEC,VEC1,0,LBLK)
       END DO
*. And new vectors from LUADD
       CALL REWINO(LUADD)
       DO IVEC = 1, NVECADD
         CALL COPVCD(LUADD,LUVEC,VEC1,0,LBLK)
       END DO
      END IF
*     ^ End of switch : NVECIN
      RETURN
      END
*
      SUBROUTINE DIIS_ACC(LUVEC_IN,NVEC,NPVEC,SCR,BMAT,LUVEC_OUT,
     &                    VEC1,VEC2,LUSCR,IACTCOR,IRESET)
*
* DIIS acceleration
*
* A set of corrections, e_k,  are residing  on LUVEC_IN, 
* combine these to an obtained improved approximation and 
* write improved eigenvector on LUVEC_OUT as first vectors
*
* IF IACTCOR.EQ.1 Then overwrite the last correction vector 
*                 with the actual difference between current 
*                 and previous approximations
*                 The previous approximation is assumed to be 
*                 the sum of the first NVEC-1 vectors
*     
* IRESET is advice given by this routine to control routine.
* IRESET = 0 suggest no reset
* IRESET = 1 suggest reset
*
* The routine minimizes the vector
*     
* sum_k w_k e_k with the constraint that sum_k w_k = 1
*
* Jeppe Olsen, Feb. 1999
*
#include "implicit.inc"
*. Scratch for holding two blocks of vectors
      DIMENSION VEC1(*),VEC2(*)
*. Scratch Atleast length 3*(NVEC+1)**2 + NVEC+1
      DIMENSION SCR(*)
*. Input/Output
      DIMENSION BMAT(*)
*
      REAL*8 INPRDD
*. Form of files
      LBLK = -1
*
      IRESET = 0
      NVECP1 = NVEC + 1
      KLINV  = 1
      KLSCR1 = KLINV  + NVECP1**2
      KLSCR2 = KLSCR1 + NVECP1**2
      KLW    = KLSCR2 + NVECP1**2
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to DIIS_ACC '
        WRITE(6,*) ' Number of vectors in subspace', NVEC
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vectors on file LUVEC_IN'
        CALL REWINO(LUVEC_IN)
        DO IVEC = 1, NVEC
          CALL WRTVCD(VEC1,LUVEC_IN,0,LBLK)
        END DO
      END IF
*. Currently no reuse of B matrix so
      IF(NPVEC.GT.0) THEN
        WRITE(6,*) ' DIIS_ACC : ERROR NPVEC .ne. 0 '
        STOP       ' DIIS_ACC : ERROR NPVEC .ne. 0 '
      END IF
*
*
*. Augment matrix B = <e_i!e_j>
*. It is assumed that the B-matrix previously has been constructed
*. for NPVEC vectors and is stored in BMAT in packed form. 
      DO IVEC = NPVEC+1, NVEC
*. Copy vector IVEC to LUSCR
        CALL SKPVCD(LUVEC_IN,IVEC-1,VEC1,1,LBLK)
        CALL REWINO(LUSCR)
        CALL COPVCD(LUVEC_IN,LUSCR,VEC1,0,LBLK)
        CALL REWINO(LUVEC_IN)
        DO JVEC = 1, IVEC
          CALL REWINO(LUSCR)
          X = INPRDD(VEC1,VEC2,LUSCR,LUVEC_IN,0,LBLK)
          BMAT(IVEC*(IVEC-1)/2+JVEC) = X
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Updated B = EE  matrix '
         CALL PRSYM(BMAT,NVEC)
      END IF
*
*. We will start by allowing individual weights for
*. all NVEC vectors. If the system is to loose -
*. as indicated by the first vectors getting a high 
*. weight, we will first combine the first two vectors 
*. Into a single vector, then the first three vectors in 
*. a single vector etc.
*
        I_DO_WCHECK = 1
        IF(I_DO_WCHECK.EQ.1) THEN
          NTRUNC_MX = 2
        ELSE
          NTRUNC_MX = 1
        END IF
        DO NTRUNC = 1, NTRUNC_MX
          NVEC_ACT   = NVEC-NTRUNC
          NVEC_ACTP1 = NVEC_ACT +1
*. Set up matrix A
*
*       ( B11  B12  B13 ... B1N  -1)
*       ( B21  B22  B23 ... B1N  -1)
*       .
*       .
*       ( BN1  BN2  BN3 ... BNN  -1)
*       (  -1  -1   -1  ... -1    0)
*
      CALL DIIS_AMAT(BMAT,SCR(KLSCR1),NVEC)
      CALL COPVEC(SCR(KLSCR1),BMAT,NVECP1*(NVECP1+1)/2)
*. Modifications connected with truncation of basis set
* B(1,1) = SUM(I=1,NTRUNC,J=1,NTRUNC) B(I,J)
      B11 = 0.0D0
      IJ = 0
      DO I = 1, NTRUNC
        DO J =1 , I-1
         IJ = IJ + 1
         B11 = B11 + 2.0D0*BMAT(IJ)
        END DO
        IJ = IJ + 1
        B11 = B11 + BMAT(IJ)
      END DO
      BMAT(1) = B11
*. B(NTRUNC,1) = 0
      IF(NTRUNC.GT.1) THEN
        BMAT(NTRUNC*(NTRUNC-1)/2+1) = 0.0D0
      END IF
*. B(I,1),I .gt. NTRUNC
      DO I = NTRUNC+1, NVEC
       BI1 = 0.0D0
*. B(I,1) = Sum(J=1,NTRUNC) B(I,J)
       DO J = 1, NTRUNC
         IJ = MAX(I,J)*(MAX(I,J)-1)/2 + MIN(I,J)
         BI1 = BI1 + BMAT(IJ)
       END DO
       BMAT(I*(I-1)/2+1) = BI1
      END DO
*
C     Zero all off diagonal elemnents for truncated vectors
      DO J = 2, NTRUNC
       DO I = J+1,NVECP1
         IJ = I*(I-1)/2+J
         BMAT(IJ) = 0.0D0
       END DO
      END DO
C?    WRITE(6,*) ' BMAT after truncation '
C?    CALL PRSYM(BMAT,NVECP1)
*
      CALL TRIPAK(SCR(KLINV),BMAT,2,NVECP1,NVECP1,1.0D0)
*. Invert 
      CALL INVERT_BY_DIAG3(SCR(KLINV),SCR(KLSCR1),SCR(KLSCR2),
     &                    SCR(KLW),NVECP1)
*. Multiply inverse matrix and vector (0,0, .. ,-1)
      CALL COPVEC(SCR(KLINV+NVEC*(NVEC+1)),SCR(KLW),NVEC)
*. Copy the weight from the first vector to vectors 2-NTRUNC 
      DO I = 2, NTRUNC
        SCR(KLW-1+I) = SCR(KLW-1+1)
      END DO
*
      ONEM = -1.0D0
      CALL SCALVE(SCR(KLW),ONEM,NVEC)
*
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Weight vector '
       CALL WRTMAT(SCR(KLW),1,NVEC,1,NVEC)
      END IF
*. Is this weight vector okay ?
*. Current criteria sum_i=1,NVEC-1 Abs(W_I) .le. ABS(W_NVEC)
      WPREV = 0.0D0
      DO I = 1, NVEC-1
       WPREV = WPREV + ABS(SCR(KLW-1+I))
      END DO
      FACTOR = 1.0D0
      WNVEC = ABS(SCR(KLW-1+NVEC))
      IF(WPREV.LE.FACTOR*WNVEC.OR.NVEC.LE.3.OR.I_DO_WCHECK.EQ.0)THEN
        I_AM_OKAY = 1
        IF(NTRUNC.EQ.1) IRESET = 0
      ELSE
        I_AM_OKAY = 0
        IF(NTRUNC.EQ.1) IRESET = 1
      END IF
      IF(NTEST.GE.1) WRITE(6,'(A,I3,2E13.6)')
     &' DIIS : nvec, w_prev, w_nvec',NVEC,WPREV,WNVEC
      IF(I_AM_OKAY.EQ.1) GOTO 1001
      END DO
*     ^ End of loop over truncation level
 1001 CONTINUE
*. Form new vector 
*. First we have only the corrections stored on LUVEC_IN, 
*. weights are given as som of sum of solution vectors, so 
*. change weights
*. sum w_k x_k = sum (n+1-k)w_k e_k
*
      DO I = 1, NVEC
        SUM = 0
        DO K = I, NVEC
          SUM = SUM + SCR(KLW-1+K)
        END DO
        SCR(KLINV-1+I) = SUM
      END DO
      IF(IACTCOR.EQ.0) THEN
        CALL MVCSMD(LUVEC_IN,SCR(KLINV),LUVEC_OUT,LUSCR,
     &              VEC1,VEC2,NVEC,1,LBLK)
      ELSE
*. Obtain correction vector
       DO IVEC = 1, NVEC - 1
         SCR(KLINV-1+IVEC) = SCR(KLINV-1+IVEC) -1.0D0
       END DO
       CALL MVCSMD(LUVEC_IN,SCR(KLINV),LUVEC_OUT,LUSCR,
     &              VEC1,VEC2,NVEC,1,LBLK)
*. And save on LUVEC_IN
       CALL SKPVCD(LUVEC_IN,NVEC-1,VEC1,1,LBLK)
       CALL REWINO(LUVEC_OUT)
       CALL COPVCD(LUVEC_OUT,LUVEC_IN,VEC1,0,LBLK)
*. And obtain new total vector on LUVEC_OUT
       ONE = 1.0D0
       CALL SETVEC(SCR(KLINV),ONE,NVEC)
       CALL MVCSMD(LUVEC_IN,SCR(KLINV),LUVEC_OUT,LUSCR,
     &              VEC1,VEC2,NVEC,1,LBLK)
      END IF
*     ^ End of IACTCOR switch 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' New approximation to solution vector'
        CALL WRTVCD(VEC1,LUVEC_OUT,1,LBLK)
      END IF
*
      RETURN
      END
*
      SUBROUTINE DIIS_AMAT(EE,A,NVEC)
*     
* Construct B matrix of DIIS method from error matrix
*
* Input :
* 
* NVEC : Number of vectors in DIIS subspace
* EE(I,J) =  <e_i!e_j>
*
* Output
*     
*   A =  
*
*       ( E11  E12  E13 ... E1N  -1)
*       ( E21  E22  E23 ... E1N  -1)
*       .
*       .
*       ( EN1  EN2  EN3 ... ENN  -1)
*       (  -1  -1   -1  ... -1    0)
*
* Input and output matrices are assumed in lower diag packed form
*
* Jeppe Olsen, Feb28, 2000
*
#include "implicit.inc"
*. Input
      DIMENSION EE(*)
*. Output
      DIMENSION A(*)
*
      LEN = NVEC*(NVEC+1)/2
      CALL COPVEC(EE,A,LEN)
      ONEM = -1.0D0
      CALL SETVEC(A(LEN+1),ONEM,NVEC)
      NVECP1 = NVEC+1
      A(NVECP1*(NVECP1+1)/2) = 0.0D0
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' A matrix from DIIS_AMAT'
        CALL PRSYM(A,NVECP1)
      END IF
*
      RETURN
      END
*
      SUBROUTINE INVERT_BY_DIAG4
      print*,'you should call INVERT_BY_DIAG2 instead'
      stop 'you should call INVERT_BY_DIAG2 instead'
      END 
*

      SUBROUTINE DUMDUM(N_CC_AMP,CCVEC1,CCVEC2)
#include "implicit.inc"
      DIMENSION CCVEC1(*),CCVEC2(*)
          do i=1,N_CC_AMP
            print*,'CCVEC1',CCVEC1(i),CCVEC2(i),i
          end do
      END
*
      SUBROUTINE FIRST_ELEMENT_LAST(VEC,NDIM)
#include "implicit.inc"
      DIMENSION VEC(NDIM)
*
      NTEST = 100
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Enter FIRST_ELEMENT_LAST'
        CALL WRTMAT(VEC,1,NDIM,1,NDIM)
      END IF
*
      A_FIRST_ELEMENT = VEC(1)
      CALL COPVEC(VEC(2),VEC,NDIM-1)
      VEC(NDIM) = A_FIRST_ELEMENT
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Exit FIRST_ELEMENT_LAST'
        CALL WRTMAT(VEC,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
*
      SUBROUTINE LAST_ELEMENT_FIRST(VEC,NDIM)
#include "implicit.inc"
      DIMENSION VEC(NDIM)
*
      NTEST = 100
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Enter LAST_ELEMENT_FIRST'
        CALL WRTMAT(VEC,1,NDIM,1,NDIM)
      END IF
*
      A_LAST_ELEMENT = VEC(NDIM)
      DO I=NDIM,2,-1
        VEC(I)= VEC(I-1)
      END DO
      VEC(1) = A_LAST_ELEMENT
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Exit LAST_ELEMENT_FIRST'
        CALL WRTMAT(VEC,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
*
      SUBROUTINE INVERT_BY_DIAG3(A,B,SCR,VEC,NDIM)
*
* Invert symmetric  - hopefully nonsingular - matrix A 
* by diagonalization
*
* Jeppe Olsen, Oct 97 to check INVMAT
*              March 00 : Scale initial matrix to obtain unit diagonal
*
#include "implicit.inc"
*. Input and output matrix
      DIMENSION A(*)       
*. Scratch matrices and vector
      DIMENSION B(*),SCR(*),VEC(*)
*
      NTEST = 00
*. Reform a to symmetric packed form
      CALL TRIPAK(A,SCR,1,NDIM,NDIM,1.0D0)
*. Extract diagonal
      CALL COPDIA(SCR,VEC,NDIM,1)
*
*.scale 
*
      DO I = 1, NDIM
*. Scaling vector
        IF(VEC(I).EQ.0.0D0) THEN
          VEC(I) = 1.0D0
        ELSE
          VEC(I) = 1.0D0/SQRT(ABS(VEC(I)))
        END IF
      END DO
*. Scale matrix
      IJ = 0
      DO I = 1, NDIM
        DO J = 1, I
          IJ = IJ + 1
          SCR(IJ) = SCR(IJ)*VEC(I)*VEC(J)
        END DO
      END DO
C     DO I = 1, NDIM
C       VEC(I) = 1.0D0/VEC(I)
C     END DO
*. Diagonalize
      CALL EIGEN(SCR,B,NDIM,0,1)
*. Scale eigenvectors
      DO IVEC = 1, NDIM 
        IOFF = 1 + (IVEC-1)*NDIM
        CALL VVTOV(B(IOFF),VEC,B(IOFF),NDIM)
      END DO
*. 
      CALL COPDIA(SCR,VEC,NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of scaled matrix : '
        CALL WRTMAT(VEC,NDIM,1,NDIM,1)
      END IF
*. Invert diagonal elements 
      DO I = 1, NDIM
       IF(ABS(VEC(I)).GT.1.0D-15) THEN
         VEC(I) = 1.0D0/VEC(I)
       ELSE
         VEC(I) = 0.0D0
         WRITE(6,*) ' Singular mode activated '
       END IF
      END DO
*. and obtain inverse matrix by transformation
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(A,B,VEC,NDIM,SCR)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverse matrix from INVERSE_BY_DIAG'
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END 
