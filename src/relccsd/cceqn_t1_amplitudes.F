      SUBROUTINE amplitude_equation_T1 (FVO,T1,T2,
     &        HO,HV,HOV,BUF1,BUF2,BUF3,NBUF3,S1)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Evaluate T1 or Lambda1 equations
C     Here we have to take care of the fact that certain integral classes
C     are now split over individual nodes. The corresponding contributions
C     which are linear in the integrals are added at the end of the T1
C     routines. Therefore we have to take care of calculating the
C     non-distributed contributions to the T1 amplitudes only on the
C     master node in order to avoid double-counting in the PARTS routine !
C     all direct VOVV/VVVV contributions and VOVV/VVVV contributions stemming
C     from intermediates have to be calculated on EACH node.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Reference----------------------------------------------
C
C     [VLD96] L. Visscher, T. J. Lee and K. G. Dyall,
C             J. Chem. Phys. 105, 8769 (1996), doi:10.1063/1.472655
C
C     The T1 equations are VLD96 Eq. (4). Terms are numbered in order
C     of appearance; implicit sums run over repeated indices:
C
C       1   F_I^A                       (added in S1S2_ADDFO)
C       2   -2 F_C^K T_K^A T_I^C
C       3   + H_C^A T_I^C               H_C^A = HV, Eq. (9)
C       4   - H_I^K T_K^A               H_I^K = HO, Eq. (10)
C       5a  + H_C^K T_IK^AC             H_C^K = HOV, Eq. (11)
C       5b  + H_C^K T_K^A T_I^C
C       6   + V_IC^AK T_K^C
C       7   + V_CD^AK tau_IK^CD  (C<D)  tau = Eq. (8)
C       8   - V_IC^KL tau_KL^AC  (K<L)
C
C     HO and HV enter with the orbital energies subtracted from their
C     diagonal (see HOINTM/HVINTM), so that the new amplitudes follow
C     from S1 by division by the energy denominators.
C
C---------------Calling variables--------------------------------------
C
      external t1eqn1
      external t1eqn2
     

      INTEGER NBUF3
      REAL*8 FVO(*)
      REAL*8 T1(*),T2(*)
      REAL*8 S1(*)
      REAL*8 HO(*),HV(*),HOV(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ihm.inc"
#include "complex.inc"
#include "ccpar.inc"
C

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface

C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ,USEDZ,RIGHT
      REAL*8 SDOT,DDOT
      INTEGER irp,k,k1,k2,m,n,off1,off2
C
C---------------Executable code--------------------------------------
C
C-------------------
C Initialize to zero
C-------------------
      CALL XCOPY (NDIMT1,A0,0,S1,1)
C--------------------------------
C Terms specific for T1 equations
C--------------------------------
                   
C ********************
C ** TASK SPLITTING **
C ********************
        IF (MYPROC.EQ.MASTER) THEN
C------------------------------------------------------------
C S(A,I) = S(A,I) + T(A,K) * [HOV(K,C) - 2 * F(K,C)] * T(C,I)
C VLD96 Eq. (4), terms 2 + 5b
C------------------------------------------------------------
         CALL XTIME(0,1,'--- T1EQNS T*[HOV - F]*T      ')
         CALL XCOPY (NFVO,HOV,1,BUF1,1)
         IF (CARITH) CALL CONJUGA (NFVO,BUF1,1)
         CALL XAXPY (NFVO,-A2,FVO,1,BUF1,1)
         CALL DELFCK ('VO','DD',BUF1)
         OFF1 = 1
         DO IRP = 1, NREP
            M = NV(IRP)
            N = NO(IRP)
            K1 = NO(IRP)
            K2 = NV(IRP)
            CALL XGEMM ('C','N',K1,N,K2,A1,BUF1(OFF1),K2,T1(OFF1),K2,
     &                   A0,BUF2,K1)
            CALL XGEMM ('N','N',M,N,K1,A1,T1(OFF1),M,BUF2,K1,
     &                   A1,S1(OFF1),M)
            OFF1 = OFF1 + M * N * RCW
         ENDDO
         CALL XTIME(0,2,'--- T1EQNS T*[HOV - F]*T      ')
        ENDIF
C ****************************************************
C ** END OF TASK SPLITTING. Parallel execution      **
C ****************************************************
C-----------------------------------------------------------------
C S(A,I) = S(A,I) + HOV(CK) * T2(A,C,I,K)
C VLD96 Eq. (4), term 5a
C-----------------------------------------------------------------
         CALL XTIME(0,1,'--- T1EQNS HOV*T2(A,C,I,K       ')
         CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,T2,BUF1)

         K    = IDIST(1,3,1)
         OFF1 = IDIST(2,3,1)
         OFF2 = IDIST(3,3,1)
         M    = MVO(1)
         CALL XGEMV ('N',M,K,A1,BUF1(OFF1),M,HOV(OFF2),1,A1,S1,1)
         CALL XTIME(0,2,'--- T1EQNS HOV*T2(A,C,I,K       ')
           

C ********************
C ** TASK SPLITTING **
C ********************

      IF (MYPROC.EQ.MASTER) THEN

      CALL XTIME(0,1,'--- T1EQNS   HV*T / T*HO              ')
C-----------------------------------------------
C S(A,I) = S(A,I) + HV(A,C) * T(C,I)
C VLD96 Eq. (4), term 3
C-----------------------------------------------
      CALL CNTRCT ('N','N',NV,NO,NV,A1,HV,T1,A1,S1,NREP)
C-----------------------------------------------
C S(A,I) = S(A,I) - T(A,K) * HO(K,I)
C VLD96 Eq. (4), term 4
C-----------------------------------------------
      CALL CNTRCT ('N','N',NV,NO,NO,-A1,T1,HO,A1,S1,NREP)
      CALL XTIME(0,2,'--- T1EQNS   HV*T / T*HO            ')

C--------------------------------------------
C S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C VLD96 Eq. (4), term 8 (sign: V(K,L,C,I) = -V_IC^KL)
C--------------------------------------------
      CALL XTIME(0,1,'--- T1EQNS VOOO*TAU                ')
      CALL GETTAU (T1,T2,BUF3)
      CALL T1EQN1 (BUF3,S1)
      CALL XTIME(0,2,'--- T1EQNS VOOO*TAU               ')

C ****************************************************
C ** END OF TASK SPLITTING. THE COMMON STUFF STARTS **
C ****************************************************
      ENDIF

C--------------------------------------------
C S(A,I) = S(A,I) + V(A,K,C,D) * TAU(C,D,I,K)   ! VOVV contributions !
C VLD96 Eq. (4), term 7
C--------------------------------------------
      CALL GETTAU (T1,T2,BUF1)
      CALL XTIME(0,1,'--- T1EQNS VOVV contribution  ')
      CALL T1EQN2 (BUF1,BUF3,NBUF3,S1)   ! VOVV contractions, term 7 !
      CALL XTIME(0,2,'--- T1EQNS VOVV contribution  ')

C--------------------------------------
C S(A,I) = S(A,I) - W(A,K,C,I) * T(C,K)    ! VOVO contributions
C VLD96 Eq. (4), term 6 (sign: W(A,K,C,I) = -V_IC^AK)
C--------------------------------------
      CALL XTIME(0,1,'--- T1EQNS VOVO * T(C,K)         ')
      CALL GETVOVO (BUF1)
      CALL DELINT ('VOVO','KDDK',BUF1,0,0,0)
      CALL SRT16 (NREP,MULTB,LTR,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            BUF1,BUF2)
      M    = MVO(1)
      K    = IDIST(1,3,1)
      OFF1 = IDIST(2,3,1)
      OFF2 = IDIST(3,3,1)
      CALL XGEMV ('N',M,K,-A1,BUF2(OFF1),M,T1(OFF2),1,A1,S1,1)
      CALL XTIME(0,2,'--- T1EQNS VOVO * T(C,K)        ')
            
      END
