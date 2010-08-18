      SUBROUTINE GENCC_F_DIAG_M_KRCC(ITSS_TP,NTSS_TP,CCDIA,ISM,
     &           VEC1,VEC2,MX_ST_TSOSO_MX,MX_ST_TSOSO_BLK_MX,
     &           N1ELINT,N2ELINT,T,
     &           WORK,KFREE,LFREE)
*
* Obtain diagonal <0![F,T]!0>
*
* Jeppe Olsen, Summer of 99,
*              Modified to allow for alpha- and beta- dependent F matrices,
*              Oct 2002 ( to see if this improves converegence for Open shell)
*
#include "implicit.inc"
      DIMENSION T(*)
#include "mxpdim.inc"
C#include "wrkspc.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "cintfo.inc"
#include "orbinp.inc"
*Added by Lasse
#include "cgas.inc"
#include "lucinp.inc"
#include "crun.inc"
#include "ipoist8.inc"
*
      DIMENSION VEC1(*),VEC2(*)
*
      DIMENSION WORK(*)
*
      IDUM = 0
      CALL MEMMAR(IDUM,IDUM,'MARK  ',IDUM,'GCC_FD')
      CALL MEMMAR(KFDIA_AL,NTOOB,'ADDL  ',2,'FDIA_A')
      CALL MEMMAR(KFDIA_BE,NTOOB,'ADDL  ',2,'FDIA_B')
      IRECOM = 1
      IF(IRECOM.NE.1) THEN
      CALL MEMMAR(KFDIA_AL_IM,NTOOB,'ADDL  ',2,'FDIA_A_IM')
      CALL MEMMAR(KFDIA_BE_IM,NTOOB,'ADDL  ',2,'FDIA_B_IM')
      END IF
*. Four blocks of string occupations
C
C     MX_ST_TSOSO_BLK_MX NEEDS TO BE CORRECTED FOR LARGER CC CALCULATIONS
C
      CALL MEMMAR(KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC1')
      CALL MEMMAR(KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC2')
      CALL MEMMAR(KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC3')
      CALL MEMMAR(KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC4')
*. Four arrays ( will be used to store sums of orbital energies)
      CALL MEMMAR(KLSTR1,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR1 ')
      CALL MEMMAR(KLSTR2,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR2 ')
      CALL MEMMAR(KLSTR3,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR3 ')
      CALL MEMMAR(KLSTR4,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR4 ')
      IF(IRECOM.NE.1) THEN
*. And for the imaginary part
      CALL MEMMAR(KLSTR5,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR5 ')
      CALL MEMMAR(KLSTR6,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR6 ')
      CALL MEMMAR(KLSTR7,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR7 ')
      CALL MEMMAR(KLSTR8,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR8 ')
      END IF
*
* 1 : Obtain Inactive Fock matrix
*
*. One-body density matrix
      LBLK = -1
C      CALL COPVCDC(LUC,LUSC1,VEC1,1,LBLK)
C      CALL WRTVCDC(VEC1,LUC,1,LBLK)
Cori      CALL DENSI2_REL(1,WORK(KRHO1),WORK(KRHO2),VEC1,VEC2,
Cori     &                LUSC1,LUC,EXPS2,0,WORK(KRHO1))
C      CALL SIGDEN_CTRL(WORK(KVEC1),WORK(KVEC2),LUC,LUHC,WORK(KT_CC),2)
C      CALL SIGDEN_CTRL(VEC1,VEC2,LUC,LUSC1,WORK(KT_CC),2)
      CALL REWINO(LUC)
      CALL REWINO(LUHC)
CTEMP      CALL SIGDEN_CTRL(VEC1,VEC2,LUC,LUSC1,WORK(KRHO1),2)
C      stop 'after sigden ctrl call'
* NOTE : THIS MAY BE USED LATER, LASSE 011205
*CTF:
*  Use sigden_ctrl
*  subroutine sigden_ctrl(C,HC,LUC,LUHC,T,ISIGDEN)
*     Feed for densities accordingly
*     ISIGDEN = 2   for densities
*CTF
*
*. And Fock matrix
C       do i=1,2
C        do j=1,NGAS
C         print*,'NGSSH,i,j',NGSSH(i,j),i,j
C        end do
C       end do
*****
* Notice harwirering of I_RES_OR_UNRES . Lasse
*****
C      DO I=1,N1ELINT+N2ELINT
C       print*,'T,I',T(I),I
C      END DO
      I_RES_OR_UNRES = 2
      IF(I_RES_OR_UNRES.EQ.1) THEN
*. Standard spinrestricted Fock matrices
C       CALL COPVEC(WORK(KINT1O),WORK(KFI),NINT1)
C       CALL FIFAM_REL(WORK(KFI))
C       CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
C       ECORE_H = 0.0D0
*. Extract diagonal of fock-matrix
C       CALL GT1DIS(WORK(KFDIA_AL),IREOTS,WORK(KPINT1),
C    &            WORK(KFI),ISMFTO,IBSO,NACOB)
C       CALL COPVEC(WORK(KFDIA_AL),WORK(KFDIA_BE),NACOB)
      ELSE
C NOTE ONLY DIAGONAL IS USED!!!
C THEREFORE ONLY THAT WILL BE GENERATED
        I_ORI = 0
        IF(I_ORI.EQ.0) THEN
          CALL FI_HS_DIA_KRCC(T,ECC,
     &               IBSO,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
     &               NGAS,NGSSH,IHPVGAS_AB,
     &               NTOOBS,NTOOB,IREOTS,IREOTS2,IDOH2,
     &               ENVIRO,ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
     &               NTOOB,WORK(KFDIA_AL),WORK(KFDIA_BE))
        call rmemchk
*
          IF(IRECOM.NE.1) THEN
            print*,'img dia'
            CALL FI_HS_DIA_KRCC(T(N1ELINT+N2ELINT+1),ECCCOM,
     &               IBSO,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
     &               NGAS,NGSSH,IHPVGAS_AB,
     &               NTOOBS,NTOOB,IREOTS,IREOTS2,IDOH2,
     &               ENVIRO,ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
     &               NTOOB,WORK(KFDIA_AL_IM),WORK(KFDIA_BE_IM))
          END IF
        ELSE
*. Spin unrestricted Fock matrices
C        CALL FI_HS_REL(WORK(KINT1O),WORK(KFI_AL),WORK(KFI_BE),ECORE_X,1)
C      Note from Lasse to self. Check if FI_AL + more is correct
C         CALL FI_HS_REL(T,ECC,IBSO,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
C    &               NGAS,NGSSH,IHPVGAS_AB,
C    &               NTOOBS,NTOOB,IREOTS,IREOTS2,IDOH2,
C    &               ENVIRO,ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
C    &               NTOOB,WORK(KFI_AL),WORK(KFI_BE))
C and now the imaginary part
C          CALL FI_HS_REL(T(N1ELINT+N2ELINT+1),ECCCOM,IBSO,NSMOB,
C         CALL FI_HS_REL(T,ECCCOM,IBSO,NSMOB,
C    &               ITPFSO,ITPFSO2,IPHGAS,
C    &               NGAS,NGSSH,IHPVGAS_AB,
C    &               NTOOBS,NTOOB,IREOTS,IREOTS2,IDOH2,
C    &               ENVIRO,ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
C    &               NTOOB,WORK(KFI_AL_IM),WORK(KFI_BE_IM))
*. Extract diagonal of fock-matrix
COri        CALL GT1DIS(WORK(KFDIA_AL),IREOTS,WORK(KPINT1),
COri     &            WORK(KFI_AL),ISMFTO,IBSO,NACOB)
COri        CALL GT1DIS(WORK(KFDIA_BE),IREOTS,WORK(KPINT1),
COri     &            WORK(KFI_BE),ISMFTO,IBSO,NACOB)
* Replaced with
C         print*,'before FOCK_DIA'
C         CALL FOCK_DIA(WORK(KFI_AL),WORK(KFDIA_AL),NACOB)
C         CALL FOCK_DIA(WORK(KFI_BE),WORK(KFDIA_BE),NACOB)
        END IF
* end replacement. Lasse
      END IF
C      print*,'real diagonal 1'
C      print*,(WORK(KFDIA_AL+II),II=0,NACOB)
C      print*,'real diagonal 2'
C      print*,(WORK(KFDIA_BE+II),II=0,NACOB)
C      print*,'img diagonal 1'
C      print*,(WORK(KFDIA_AL_IM+II),II=0,NACOB)
C      print*,'img diagonal 2'
c      print*,(WORK(KFDIA_BE_IM+II),II=0,NACOB)
*. And then the diagonal
C NOTE : THIS ROUTINE SHOULD PERHAPS BE CHANGED, LASSE
      call rmemchk
C      print*,'before GENCC_F_DIAG_REL'
C      do i=1,100
C        print*,'me'
C      end do
      IF(IRECOM.EQ.1) THEN
      CALL GENCC_F_DIAG_KRCC(ITSS_TP,NTSS_TP,CCDIA,ISM,N_CC_AMP,IRECOM,
     &     WORK(KFDIA_AL),WORK(KFDIA_BE),
     &     KDUM,KDUM,
     &     WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
     &     WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),
     &     WORK(KLSTR1),WORK(KLSTR2),WORK(KLSTR3),WORK(KLSTR4),
     &     KDUM,KDUM,KDUM,KDUM,
     &     WORK,KFREE,LFREE)
      ELSE
      CALL GENCC_F_DIAG_KRCC(ITSS_TP,NTSS_TP,CCDIA,ISM,N_CC_AMP,IRECOM,
     &     WORK(KFDIA_AL),WORK(KFDIA_BE),
     &     WORK(KFDIA_AL_IM),WORK(KFDIA_BE_IM),
     &     WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
     &     WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),
     &     WORK(KLSTR1),WORK(KLSTR2),WORK(KLSTR3),WORK(KLSTR4),
     &     WORK(KLSTR5),WORK(KLSTR6),WORK(KLSTR7),WORK(KLSTR8),
     &     WORK,KFREE,LFREE)
      END IF
C      CALL GENCC_F_DIAG_REL(ITSS_TP,NTSS_TP,CCDIA(NTSS_TP+1),ISM,
C     &     WORK(KFDIA_AL),WORK(KFDIA_BE),
C     &     WORK(KFI_AL_IM),WORK(KFI_BE_IM),
C     &     WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
C     &     WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),
C     &     WORK(KLSTR1),WORK(KLSTR2),WORK(KLSTR3),WORK(KLSTR4),
C     &     WORK(KLSTR5),WORK(KLSTR6),WORK(KLSTR7),WORK(KLSTR8))
*
C?    WRITE(6,*) ' Enforced stop in GENCC_F_DIAG_M '
C?    STOP '  Enforced stop in GENCC_F_DIAG_M'
      CALL MEMMAR(IDUM,IDUM,'FLUSM ',IDUM,'GCC_FD')
*
      RETURN
      END
*
C      SUBROUTINE FI_HS_REL(T,FI_AL,FI_BE,ECC,IDOH2)
      SUBROUTINE FI_HS_DIA_KRCC(T,ECC,IOBSM,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
     &                   NGAS,NGSSH,IHPVGAS_AB,
     &                   LOBSM,NORBT,IREOTS,IREOTS2,IDOH2,
     &                   ENVIRO,ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
     &                   NTOOB,FI_AL,FI_BE)
*
* Obtain core energy and alpha and beta part of inactive Fock-matrix
* / effective 1- electron operator for high spin openshell
* case
*
* Jeppe Olsen, July 3, 2002
*
* The formulae goes
*
* FI_ALPHA(I,J) = H(I,J)
*               + SUM(K) (RHO_ALPHA(K,K) + RHO_BETA(K,K)) G(IJKK)
*               - SUM(K)  RHO_ALPHA(K,K) * G(IKKJ)
*               - SUM(K)  RHO_ALPHA(K,K) * G(IK_K_J) in rel. case
* FI_ALPHA(I_,J) = H(I_,J)
*               + SUM(K) (RHO_ALPHA(K,K) + RHO_BETA(K,K)) G(I_JKK)
*               - SUM(K)  RHO_ALPHA(K,K) * G(I_KKJ)
*               - SUM(K)  RHO_ALPHA(K,K) * G(I_K_K_J) in rel. case
* + rest
*
* ECORE = sum(i) ((RHO_ALPHA(I,I) + RHO_BETA (I,I)) H(II)
*       + 1/2 sum(i,j) (RHO_ALPHA(I,I)*RHO_ALPHA(J,J)) ( G(IIJJ) - G(IJJI) )
*       + 1/2 sum(i,j) (RHO_BETA (I,I)*RHO_BETA (J,J)) ( G(IIJJ) - G(IJJI) )
*       +     sum(i,j) RHO_ALPHA(I,I)*RHO_BETA(J,J)*G(IIJJ)
*
* The alpha and beta parts of the density matrices are not accessed,
* instead is IHPVGAS_AB used to obtain the occupied spaces
*
* Lasse, July 2005
*
* Changed to a structure similar to that in fih_rel
*
*
#include "implicit.inc"
*. General input
#include "wrkspc.inc"
#include "mxpdim.inc"
#include "ipoist8.inc"
#include "ctcc.inc"
#include "glbbas.inc"
#include "gasstr.inc"
C#include "cgas.inc"
C#include "orbinp.inc"
C#include "lucinp.inc"
*. Standard 1+2 electron integrals
      DIMENSION T(*)
*. Output
      DIMENSION FI_AL(*), FI_BE(*)
*. Input
      CHARACTER*6 ENVIRO
      INTEGER IOBSM(*),LOBSM(*),IREOTS(*),IREOTS2(*)
      INTEGER NGSSH(MXPIRR,NGAS)
      INTEGER ITPFSO(*), ITPFSO2(*), IPHGAS(*), IKPERM(4)
      INTEGER NGSOB(MXPOBS,*), NGSOB2(MXPOBS,*)
      INTEGER NGSOBT(*), NGSOBT2(*), IHPVGAS_AB(MXPNGAS,2)
*
* Is matrix symmetrically packed or ???
      ISYMPACK = 1
      ZERO = 0.0D0
C      print*,'N1DIAINT',N1DIAINT
      CALL SETVEC(FI_AL,ZERO,N1DIAINT/2)
      CALL SETVEC(FI_BE,ZERO,N1DIAINT/2)
      call rmemchk
*
      NTEST = 0
C      print*,'NTOOBS,NTOOBS,NSMOB',NTOOBS,NTOOBS,NSMOB
C      print*,'NORBT',NORBT
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' FI_HS : Initial 1 e integral list '
        CALL APRBLM2(T,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
      END IF
*
* Core energy
*
      ECC = 0.0D0
      IJSM = 1
      IIOFF = 0
      NINT12 = N1ELINT + N2ELINT
      DO I=1,NORBT
        FI_AL(I) = 0.0D0
        FI_BE(I) = 0.0D0
      END DO
*
* GAS-SYMM ordering (GAS is outer loop!)
      IMARK = 1
      do IGAS=1,NGAS,1
        do ISMOB=1,NSMOB/2,1
* Determine index range for this symmetry/GAS block:
          NORBSM = 0
          do IGASS=1,IGAS-1,1
            do ISM=1,NSMOB/2,1
              NORBSM = NORBSM + NGSSH(ISM,IGASS)
            end do
          end do
          do JSM=1,ISMOB-1,1
            NORBSM = NORBSM + NGSSH(JSM,IGAS)
          end do
          ISTART = NORBSM + 1
          IEND = NORBSM + NGSSH(ISMOB,IGAS)
          if (NTEST.ge.50) then
            write(6,'(2X,A,5I4)') 'IGAS,ISMOB,NORBSM,ISTART,IEND',
     &                             IGAS,ISMOB,NORBSM,ISTART,IEND
            write(6,*) 'N1ELINT,NGSSH(ISMOB,IGAS)',
     &                  N1ELINT,NGSSH(ISMOB,IGAS)
          end if
C          print*,'IPHGAS(ITPFSO(ISTART)',IPHGAS(ITPFSO(ISTART))
C
          FACTOR_AL = 0.0D0
          FACTOR_BE = 0.0D0
          ICOUNT = ITPFSO(ISTART)
CTEMP          IF(IPHGAS(ITPFSO(ISTART)).EQ.2) then
            IUB = 1
C            Print*,'IMARK',IMARK
            call geth1dia_sigden_KRCC(IUB,FI_AL(IMARK),T,
     &                           WORK(KMDIASV),WORK(KPOSDIA),
     &                           N1DIAINT,NDIAINT,NGSSH(ISMOB,IGAS),
     &                           ISTART,IEND,IREOTS,NTEST)
C          print*,'FI_AL(IMARK),IMARK',FI_AL(IMARK),IMARK
C         IF(IPHGAS(ITPFSO(I)).EQ.2) ECC = ECC + 2*FI(II)
C            call copvec(WORK(KLH1),H1DIA_U(IOFFSD),NGSSH(ISMOB,IGAS))
C           print*,'NGSSH,ICOUNT',NGSSH(ISMOB,IGAS),ICOUNT
           IF(IHPVGAS_AB(ICOUNT,1).EQ.1) THEN
             FACTOR_AL = 1.0D0
             DO II=1,NGSSH(ISMOB,IGAS)
C               IF(IHPVGAS_AB(II,1).EQ.1) FACTOR_AL = 1.0D0
               ECC = ECC + FACTOR_AL * FI_AL(II)
C               FACTOR_AL = 0.0D0
C             print*,'ECC,IHPVGAS_AB(II,1)',ECC,IHPVGAS_AB(II,1)
             END DO
           END IF
            IUB = 2
            call geth1dia_sigden_KRCC(IUB,FI_BE(IMARK),T,
     &                           WORK(KMDIASV),WORK(KPOSDIA),
     &                           N1DIAINT,NDIAINT,NGSSH(ISMOB,IGAS),
     &                           ISTART,IEND,IREOTS2,NTEST)
C            call copvec(WORK(KLH1),H1DIA_B(IOFFSD),NGSSH(ISMOB,IGAS))
           IF(IHPVGAS_AB(ICOUNT,2).EQ.1) THEN
             FACTOR_BE = 1.0D0
             DO II=1,NGSSH(ISMOB,IGAS)
C               IF(IHPVGAS_AB(II,2).EQ.1) FACTOR_BE = 1.0D0
               ECC = ECC + FACTOR_BE * FI_BE(II)
C               FACTOR_BE = 0.0D0
             END DO
           END IF
CTEMP          end if
C          IOFFSD = IOFFSD + NGSSH(ISMOB,IGAS)
            IMARK = IMARK + NGSSH(ISMOB,IGAS)
        end do
      end do
C      print*,'ECC',ECC
      call rmemchk
*
* Testing
*. Number of 1-e integrals
C  NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
C        print*,'N1ELINT,NALLINT,N2ELINT',N1ELINT,NALLINT,N2ELINT
C        CALL COPVEC(WORK(KLH1),FI_AL,N1DIAINT)
C        CALL COPVEC(WORK(KLH1+N1DIAINT),FI_BE,N1DIAINT)
C        stop
C        CALL COPVEC(T(N1ELINT+N2ELINT+1),FI_AL_IM,N1ELINT/2)
C        CALL COPVEC(T(N1ELINT+N2ELINT+N1ELINT/2+1),FI_BE_IM,N1ELINT/2)
*
*
* 2-electron
*
* Coulomb case
        IKFAC = 1
        IKPERM(1) = 1
        IKPERM(2) = 2
        IKPERM(3) = 3
        IKPERM(4) = 4
        if (ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0) then
* Storage mode for integrals contains "-" sign for beta/barred indices.
          IKFAC = -1
* Use Kramers symmetry on second index pair.
          IKPERM(1) = 1
          IKPERM(2) = 2
          IKPERM(3) = 4
          IKPERM(4) = 3
        end if
C        print*,'ISPINFREE,ENVIRO,IKPERM(4)',ISPINFREE,ENVIRO,IKPERM(4)
        IGSOFF = 0
        JGSOFF = 0
        ISMOFF = 0
        JSMOFF = 0
C        IJ = 1 + N1DIAINT/2
        IJ = 1
        do IGS = 1,NGAS,1
C          print*,'NGASOCC,IGS',NGASOCC(IGS),IGS
          IF(NGASOCC(IGS).eq.1) THEN
          if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
          do JGS = 1,NGAS,1
C          print*,'NGASOCC(JGS),JGS',NGASOCC(JGS),JGS
          IF(NGASOCC(JGS).eq.1) THEN
            if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT(JGS-1)
            do ISM = NSMOB/2+1, NSMOB
              if (ISM.gt.NSMOB/2+1) ISMOFF = ISMOFF + NGSOB(ISM-1,IGS)
              do JSM = NSMOB/2+1, NSMOB
                if (JSM.gt.NSMOB/2+1) JSMOFF = JSMOFF + NGSOB(JSM-1,JGS)
                NI = NGSOB(ISM,IGS)
                NJ = NGSOB(JSM,JGS)
                IOFF = IGSOFF + ISMOFF
                JOFF = JGSOFF + JSMOFF
                if (NTEST.ge.10) then
                  write(6,*) 'IGS,JGS,ISM,JSM ',IGS,JGS,ISM,JSM
                  write(6,*) ' Index range :',1+IOFF,'...',NI+IOFF
                  write(6,*) ' Index range :',1+JOFF,'...',NJ+JOFF
                end if
C                IF(IGS.EQ.2.and.JGS.eq.2) then
C                  print*,'IHPVGAS_AB(ITPFSO(IOFF+1),1)',
C     &                    IHPVGAS_AB(ITPFSO(IOFF+1),1)
C                  print*,'IHPVGAS_AB(ITPFSO(IOFF+1),2)',
C     &                    IHPVGAS_AB(ITPFSO(IOFF+1),2)
C                end if
                IF(IHPVGAS_AB(ITPFSO(IOFF+1),1).EQ.1) THEN
C     &             IHPVGAS_AB(ITPFSO(JOFF+1),1).EQ.1) THEN
                do I = 1+IOFF,NI+IOFF,1
                  do J = 1+JOFF,NJ+JOFF,1
C                    print*,'I,J,NINT12,N1INACINT,NINACINT,NTEST',
C     &                      I,J,NINT12,N1INACINT,NINACINT,NTEST
*. Factor for Coulomb and exchange
C              FACTOR_C = 0.0D0
C              FACTOR_E = 0.0D0
C              IF(IHPVGAS_AB(IGAS,1).EQ.1.AND.IHPVGAS_AB(JGAS,1).EQ.1)
C     &        THEN
C               FACTOR_C = 0.5D0
C              END IF
C              IF(IHPVGAS_AB(IGAS,2).EQ.1.AND.IHPVGAS_AB(JGAS,2).EQ.1)
C     &        THEN
C               FACTOR_C = FACTOR_C + 0.5D0
C              END IF
C              IF(IHPVGAS_AB(IGAS,1).EQ.1.AND.IHPVGAS_AB(JGAS,2).EQ.1)
C     &        THEN
C               FACTOR_C = FACTOR_C + 1.0D0
C              END IF
C              IF(FACTOR_C.NE.0.0D0)
C     &        ECC = ECC + FACTOR_C*GTIJKL(IP,IP,JP,JP)
                    FACTOR_C = 0.0D0
                    FACTOR_CE = 0.0D0
                    IF(IHPVGAS_AB(ITPFSO(I),1).EQ.1.AND.
     &                 IHPVGAS_AB(ITPFSO(J),1).EQ.1) THEN
                      FACTOR_CE = 0.5D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),2).EQ.1.AND.
     &                 IHPVGAS_AB(ITPFSO(J),2).EQ.1) THEN
                      FACTOR_CE = FACTOR_CE + 0.5D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),1).EQ.1.AND.
     &                 IHPVGAS_AB(ITPFSO(J),2).EQ.1) THEN
                      FACTOR_CE = FACTOR_CE + 1.0D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),1).EQ.1) THEN
                      FACTOR_C = 1.0D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),2).EQ.1) THEN
                      FACTOR_C = FACTOR_C + 1.0D0
                    END IF
                    call geth2dia_sigden_KRCC(ECOU,T,
C     &                                   WORK(KT_CC),
     &                                   I,J,
     &                                   'COUL',WORK(KLIBSOBEX),NINT12,
     &                                   I2DIATP,IKFAC,IKPERM,
     &                                   IREOTS,IREOTS,
     &                                   WORK(KMDIASV),WORK(KPOSDIA),
     &                                   N1DIAINT,NDIAINT,NTEST)
                    ECC = ECC + FACTOR_CE*ECOU
C                    print*,'FACTOR_C,ECOU',FACTOR_CE,ECOU
                    IJ = J
C                    print*,'IJ',IJ
                    FI_AL(IJ) =
     &              FI_AL(IJ) + FACTOR_C*ECOU
                    FI_BE(IJ) =
     &              FI_BE(IJ) + FACTOR_C*ECOU
C      print*,'ECC,ECOU,FACTOR_C,FACTOR_C2,FACTOR_C3',
C     &        ECC,ECOU,FACTOR_C,FACTOR_C2,FACTOR_C3
                  end do
C                    IF(ECOU.NE.0.0D0) THEN
C                    IJ = IJ + 1
C                    END IF
C                    PRINT*,'IJ',IJ
                    ECOU = 0.0D0
                end do
                END IF
C                END IF
C                END IF
              end do
              JSMOFF = 0
            end do
            ISMOFF = 0
            END IF
          end do
          JGSOFF = 0
          ELSE IF(NGASOCC(IGS).eq.0.and.IGS.eq.1) THEN
            DO I=1,NTOOB
              DO J=1,NTOOB
                REPLACED = ZERO
              END DO
            END DO
          END IF
        end do
C       print*,'ECC',ECC
       call rmemchk
*
* Exchange
*
        IKFAC = 1
        IKPERM(1) = 1
        IKPERM(2) = 2
        IKPERM(3) = 3
        IKPERM(4) = 4
        if (ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0) then
          IKFAC = -1
* Use Kramers symmetry on second index pair.
          IKPERM(1) = 1
          IKPERM(2) = 2
          IKPERM(3) = 4
          IKPERM(4) = 3
        end if
        IGSOFF = 0
        JGSOFF = 0
        ISMOFF = 0
        JSMOFF = 0
        IJ = 1
        do IGS = 1,NGAS,1
          IF(NGASOCC(IGS).eq.1) THEN
          if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
          do JGS = 1,NGAS,1
            IF(NGASOCC(JGS).eq.1) THEN
            if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT(JGS-1)
            do ISM = NSMOB/2+1, NSMOB
              if (ISM.gt.NSMOB/2+1) ISMOFF = ISMOFF + NGSOB(ISM-1,IGS)
              do JSM = NSMOB/2+1, NSMOB
                if (JSM.gt.NSMOB/2+1) JSMOFF = JSMOFF + NGSOB(JSM-1,JGS)
                NI = NGSOB(ISM,IGS)
                NJ = NGSOB(JSM,JGS)
                IOFF = IGSOFF + ISMOFF
                JOFF = JGSOFF + JSMOFF
                if (NTEST.ge.10) then
                  write(6,*) 'IGS,JGS,ISM,JSM ',IGS,JGS,ISM,JSM
                  write(6,*) ' Index range :',1+IOFF,'...',NI+IOFF
                  write(6,*) ' Index range :',1+JOFF,'...',NJ+JOFF
                end if
                IF(IHPVGAS_AB(ITPFSO(IOFF+1),1).EQ.1) THEN
C     &             IHPVGAS_AB(ITPFSO(JOFF+1),1).EQ.1) THEN
                do I = 1+IOFF,NI+IOFF,1
                  do J = 1+JOFF,NJ+JOFF,1
                    FACTOR_E_AL = 0.0D0
                    FACTOR_E_BE = 0.0D0
                    FACTOR_E = 0.0D0
                    IF(IHPVGAS_AB(ITPFSO(I),1).EQ.1.AND.
     &                 IHPVGAS_AB(ITPFSO(J),1).EQ.1) THEN
                      FACTOR_E = -0.5D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),2).EQ.1.AND.
     &                 IHPVGAS_AB(ITPFSO(J),2).EQ.1) THEN
                      FACTOR_E = FACTOR_E - 0.5D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),1).EQ.1) THEN
                      FACTOR_E_AL = -1.0D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),2).EQ.1) THEN
                      FACTOR_E_BE = -1.0D0
                    END IF
                    call geth2dia_sigden_KRCC(EEXC,T,
C     &                                   WORK(KT_CC),
     &                                   I,J,
     &                                   'EXC1',WORK(KLIBSOBEX),NINT12,
     &                                   I2DIATP,IKFAC,IKPERM,
     &                                   IREOTS,IREOTS,
     &                                   WORK(KMDIASV),WORK(KPOSDIA),
     &                                   N1DIAINT,NDIAINT,NTEST)
C                    print*,'FACTOR_E',FACTOR_E
                    ECC = ECC + FACTOR_E*EEXC
                    IJ = J
C                    print*,'IJ',IJ
                    FI_AL(IJ) = FI_AL(IJ) + FACTOR_E_AL*EEXC
                    FI_BE(IJ) = FI_BE(IJ) + FACTOR_E_BE*EEXC
                    EEXC = 0.0D0
                  end do
C                    IF(EEXC.NE.0.0D0) THEN
C                    IJ = IJ + 1
C                    END IF
C                    PRINT*,'IJ',IJ
C                    EEXC = 0.0D0
                end do
                END IF
C                END IF
              end do
              JSMOFF = 0
            end do
            ISMOFF = 0
            END IF
          end do
          JGSOFF = 0
          ELSE IF(NGASOCC(IGS).eq.0.and.IGS.eq.1) THEN
            DO I=1,NTOOB
              DO J=1,NTOOB
                REPLACED = ZERO
              END DO
            END DO
          END IF
        end do
C       print*,'ECC',ECC
       call rmemchk
*
* Exchange 2
*
          IKFAC = 1
          IKPERM(1) = 1
          IKPERM(2) = 2
          IKPERM(3) = 3
          IKPERM(4) = 4
          if (ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0) then
            IKFAC = -1
            IKPERM(1) = 1
            IKPERM(2) = 2
            IKPERM(3) = 3
            IKPERM(4) = 4
          end if
          IGSOFF = 0
          JGSOFF = 0
          ISMOFF = 0
          JSMOFF = 0
          IJ = 1
          do IGS = 1,NGAS,1
            IF(NGASOCC(IGS).eq.1) THEN
            if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
            do JGS = 1,NGAS,1
              IF(NGASOCC(JGS).eq.1) THEN
              if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT2(JGS-1)
              do ISM = NSMOB/2+1, NSMOB
                if (ISM.gt.NSMOB/2+1) ISMOFF = ISMOFF + NGSOB(ISM-1,IGS)
                do JSM = NSMOB/2+1, NSMOB
                  if (JSM.gt.NSMOB/2+1)
     &                       JSMOFF = JSMOFF + NGSOB2(JSM-1,JGS)
                  NI = NGSOB(ISM,IGS)
                  NJ = NGSOB2(JSM,JGS)
                  IOFF = IGSOFF + ISMOFF
                  JOFF = JGSOFF + JSMOFF
                  if (NTEST.ge.10) then
                    write(6,*) 'IGS,JGS,ISM,JSM ',IGS,JGS,ISM,JSM
                    write(6,*) ' Index range :',1+IOFF,'...',NI+IOFF
                    write(6,*) ' Index range :',1+JOFF,'...',NJ+JOFF
                  end if
                IF(IHPVGAS_AB(ITPFSO(IOFF+1),1).EQ.1) THEN
C     &             IHPVGAS_AB(ITPFSO(JOFF+1),1).EQ.1) THEN
                  do I = 1+IOFF,NI+IOFF,1
                    do J = 1+JOFF,NJ+JOFF,1
                    FACTOR_E_AL = 0.0D0
                    FACTOR_E_BE = 0.0D0
                    FACTOR_E = 0.0D0
                    IF(IHPVGAS_AB(ITPFSO(I),1).EQ.1.AND.
     &                 IHPVGAS_AB(ITPFSO(J),1).EQ.1) THEN
                      FACTOR_E = -0.5D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),2).EQ.1.AND.
     &                 IHPVGAS_AB(ITPFSO(J),2).EQ.1) THEN
                      FACTOR_E = FACTOR_E - 0.5D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),1).EQ.1) THEN
                      FACTOR_E_AL = -1.0D0
                    END IF
                    IF(IHPVGAS_AB(ITPFSO(I),2).EQ.1) THEN
                      FACTOR_E_BE = -1.0D0
                    END IF
                      call geth2dia_sigden_KRCC(EEXC2,T,
C     &                                   WORK(KT_CC),
     &                                   I,J,
     &                                   'EXC2',WORK(KLIBSOBEX),NINT12,
     &                                   I2DIATP,IKFAC,IKPERM,
     &                                   IREOTS,IREOTS2,
     &                                   WORK(KMDIASV),WORK(KPOSDIA),
     &                                   N1DIAINT,NDIAINT,NTEST)
C                    print*,'EEXC2,FACTOR_E',EEXC2,FACTOR_E
                    ECC = ECC + FACTOR_E*EEXC2
                    IJ = J
C                    print*,'IJ',IJ
                    FI_AL(IJ) = FI_AL(IJ) + FACTOR_E_AL*EEXC2
                    FI_BE(IJ) = FI_BE(IJ) + FACTOR_E_BE*EEXC2
                    EEXC2 = 0.0D0
                    end do
C                  IJ = IJ + 1
C                  PRINT*,'IJ',IJ
C                    IF(EEXC2.NE.0.0D0) THEN
C                    IJ = IJ + 1
C                    END IF
C                    PRINT*,'IJ',IJ
C                    EEXC2 = 0.0D0
                  end do
C                  IJ = IJ + 1
C                  PRINT*,'IJ',IJ
                END IF
C                END IF
                end do
                JSMOFF = 0
              end do
              ISMOFF = 0
              END IF
            end do
            JGSOFF = 0
            ELSE IF(NGASOCC(IGS).eq.0.and.IGS.eq.1) THEN
              DO I=1,NTOOB
                DO J=1,NTOOB
                  REPLACED = ZERO
                END DO
              END DO
            ELSE IF(NGASOCC(IGS).eq.0.and.IGS.ne.1.and.NGASOCC(1).ne.0)
     &        THEN
              WRITE(6,*) 'WARNING FROM GET_DIA_2E'
              WRITE(6,*) 'CHECK YOUR 2-ELECTRON INTEGRALS'
              WRITE(6,*) 'FOR CI DIAGONAL'
              WRITE(6,*) 'LASSE 25-11-04'
            END IF
          end do
C       print*,'ECC',ECC
       call rmemchk
*
      IF(NTEST.GT.111110.AND.IDOH2.NE.0) THEN
       WRITE(6,*) ' FI_alpha in Symmetry blocked form '
       WRITE(6,*) ' =================================='
       WRITE(6,*)
       ISYM = 0
       CALL APRBLM2(FI_AL,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
       WRITE(6,*) ' FI_beta in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*)
       ISYM = 0
       CALL APRBLM2(FI_BE,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
      END IF
      print*,'ECC',ECC
*
      RETURN
      END
*
*       CALL PLACE_1_EL_INT(FI_BE,N_BE)
      SUBROUTINE PLACE_1_EL_INT(ARR,N)
#include "implicit.inc"
      DIMENSION ARR(N**2)
      IF (N.GT.1) THEN
        DO I=N,2,-1
          ARR(I**2) = ARR(I)
C          print*,'ARR(I)',ARR(I),I
          ARR(I) = 0.0D0
        END DO
       END IF
      RETURN
      END
*
      SUBROUTINE GENCC_F_DIAG_KRCC(ITSS_TP,NTSS_TP,CCDIA,
     &           ISM,N_CC_AMP,IRECOM,
     &           FDIA_AL,FDIA_BE,FDIA_AL_IM,FDIA_BE_IM,
     &           IOCC_CA,IOCC_CB,IOCC_AA,IOCC_AB,
     &           E_CA,E_CB,E_AA,E_AB,
     &           E_CA_IM,E_CB_IM,E_AA_IM,E_AB_IM,
     &           WORK,KFREE,LFREE)
*
* <0![F,T]!0)
*
* Jeppe Olsen, Summer of 99
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "cc_exc.inc"
#include "symm.inc"
*. Specific input
      INTEGER ITSS_TP(4*NGAS,NTSS_TP)
      DIMENSION FDIA_AL(*), FDIA_BE(*), FDIA_AL_IM(*), FDIA_BE_IM(*)
*. Output
      DIMENSION CCDIA(*)
*. Scratch
      INTEGER IOCC_CA(*),IOCC_CB(*),IOCC_AA(*),IOCC_AB(*)
      DIMENSION E_CA(*), E_CB(*), E_AA(*), E_AB(*)
      DIMENSION E_CA_IM(*), E_CB_IM(*), E_AA_IM(*), E_AB_IM(*)
*. Local scratch
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS)
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
*
      DIMENSION WORK(*)
*
      NTEST = 00
      IT = 0
C      print*,'NTSS_TP',NTSS_TP
C      II = 4*NGAS
C      DO I =1, II
C        DO J = 1, NTSS_TP
C          print*,'ITSS_TP',ITSS_TP(I,J),I,J
C        END DO
C      END DO
      DO ITSS = 1, NTSS_TP
*. Transform from occupations to groups
       CALL OCC_TO_GRP_CC_KRCC(ITSS_TP(1+0*NGAS,ITSS),IGRP_CA,1      )
       CALL OCC_TO_GRP_CC_KRCC(ITSS_TP(1+1*NGAS,ITSS),IGRP_CB,1      )
       CALL OCC_TO_GRP_CC_KRCC(ITSS_TP(1+2*NGAS,ITSS),IGRP_AA,1      )
       CALL OCC_TO_GRP_CC_KRCC(ITSS_TP(1+3*NGAS,ITSS),IGRP_AB,1      )
*
       NEL_CA = IELSUM(ITSS_TP(1+0*NGAS,ITSS),NGAS)
       NEL_CB = IELSUM(ITSS_TP(1+1*NGAS,ITSS),NGAS)
       NEL_AA = IELSUM(ITSS_TP(1+2*NGAS,ITSS),NGAS)
       NEL_AB = IELSUM(ITSS_TP(1+3*NGAS,ITSS),NGAS)
C       print*,'NEL_CA,NEL_CB,NEL_AA,NEL_AB',
C     &         NEL_CA,NEL_CB,NEL_AA,NEL_AB
*
* find these in gnside_rel (LUCIAREL) in sigden_r.F
*
*. Diagonal block ?
       CALL DIAG_EXC_KRCC(ITSS_TP(1+0*NGAS,ITSS),
     &                  ITSS_TP(1+1*NGAS,ITSS),
     &                  ITSS_TP(1+2*NGAS,ITSS),
     &                  ITSS_TP(1+3*NGAS,ITSS),NGAS,IDIAG)
       IF(MSCOMB_CC.EQ.0.OR.IDIAG.EQ.0) THEN
         IRESTRICT = 0
       ELSE
         IRESTRICT = 1
       END IF
       DO ISM_C = 1, NSMST
Cori     ISM_A = MULTD2H(ISM,ISM_C)
C         ISM_A = IDBGMULT(ISM_C,IADJSYM(ISM))
C         ISM_A = IDBGMULT(ISM,IADJSYM(ISM))
C         ISM_A = IDBGMULT(ISM,IADJSYM(ISM_C))
         ISM_A = IDBGMULT(ISM,ISM_C)
* ADJSYM(MULTD2H_DBG)
         DO ISM_CA = 1, NSMST
Cori       ISM_CB = MULTD2H(ISM_C,ISM_CA)
           ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
           DO ISM_AA = 1, NSMST
Cori        ISM_AB =  MULTD2H(ISM_A,ISM_AA)
            ISM_AB =  IDBGMULT(ISM_A,INVELM(ISM_AA))
*
            ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
            ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
            IF(IRESTRICT.EQ.1.AND.ISM_BETA.GT.ISM_ALPHA) GOTO 777
            IF(IRESTRICT.EQ.0.OR.ISM_ALPHA.GT.ISM_BETA) THEN
             IRESTRICT_LOOP = 0
            ELSE
             IRESTRICT_LOOP = 1
            END IF
C?          WRITE(6,'(A,4I5)') ' ISM_CA, ISM_CB, ISM_AA, ISM_AB',
C?   &                           ISM_CA, ISM_CB, ISM_AA, ISM_AB
*. obtain strings IUB added
            IUB = 1
            CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CA,NGAS,ISM_CA,NEL_CA,
     &           NSTR_CA,IOCC_CA, NORBT,0,IDUM,IDUM,
     &           WORK,KFREE,LFREE)
            IUB = 2
            CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CB,NGAS,ISM_CB,NEL_CB,
     &           NSTR_CB,IOCC_CB, NORBT,0,IDUM,IDUM,
     &           WORK,KFREE,LFREE)
            IUB = 1
            CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AA,NGAS,ISM_AA,NEL_AA,
     &           NSTR_AA,IOCC_AA, NORBT,0,IDUM,IDUM,
     &           WORK,KFREE,LFREE)
            IUB = 2
            CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AB,NGAS,ISM_AB,NEL_AB,
     &           NSTR_AB,IOCC_AB, NORBT,0,IDUM,IDUM,
     &           WORK,KFREE,LFREE)
C     GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
C    &                              NORBT,IDOREO,IZ,IREO)
*. Set up sums of F elements
             CALL SUM_FDIA_FOR_STR(E_CA,FDIA_AL,NSTR_CA,NEL_CA,IOCC_CA)
             CALL SUM_FDIA_FOR_STR(E_CB,FDIA_BE,NSTR_CB,NEL_CB,IOCC_CB)
             CALL SUM_FDIA_FOR_STR(E_AA,FDIA_AL,NSTR_AA,NEL_AA,IOCC_AA)
             CALL SUM_FDIA_FOR_STR(E_AB,FDIA_BE,NSTR_AB,NEL_AB,IOCC_AB)
*. And for the imaginary part added by Lasse
        IF(IRECOM.NE.1) THEN
        CALL SUM_FDIA_FOR_STR(E_CA_IM,FDIA_AL_IM,NSTR_CA,NEL_CA,IOCC_CA)
        CALL SUM_FDIA_FOR_STR(E_CB_IM,FDIA_BE_IM,NSTR_CB,NEL_CB,IOCC_CB)
        CALL SUM_FDIA_FOR_STR(E_AA_IM,FDIA_AL_IM,NSTR_AA,NEL_AA,IOCC_AA)
        CALL SUM_FDIA_FOR_STR(E_AB_IM,FDIA_BE_IM,NSTR_AB,NEL_AB,IOCC_AB)
        END IF
*. Loop over T elements as matrix T(I_CA, I_CB, IAA, I_AB)
             DO I_AB = 1, NSTR_AB
              IF(IRESTRICT_LOOP.EQ.1) THEN
                I_AA_MIN = I_AB
              ELSE
                I_AA_MIN = 1
              END IF
              DO I_AA = I_AA_MIN, NSTR_AA
               DO I_CB = 1, NSTR_CB
                IF(IRESTRICT_LOOP.EQ.1.AND.I_AB.EQ.I_AA) THEN
                 ICA_MIN = I_CB
                ELSE
                 ICA_MIN = 1
                END IF
                DO I_CA = ICA_MIN, NSTR_CA
                 IT = IT + 1
                 CCDIA(IT) =
     &           E_CA(I_CA) + E_CB(I_CB) - E_AA(I_AA) - E_AB(I_AB)
C                 print*,'E_CA(I_CA),E_CB(I_CB),E_AA(I_AA),E_AB(I_AB)'
C                 print*, E_CA(I_CA),E_CB(I_CB),E_AA(I_AA),E_AB(I_AB)
                 IF(IRECOM.NE.1) THEN
                 CCDIA(IT+N_CC_AMP) =
     &           E_CA_IM(I_CA)+E_CB_IM(I_CB)-E_AA_IM(I_AA)-E_AB_IM(I_AB)
                 END IF
                 YTEST = 1.D-10
*. A temp measure for eliminating some redundant operations
C                 IF(ABS(CCDIA(IT)).LE.YTEST) THEN
C                   CCDIA(IT) = 1.0D+10
C                 END IF
C?               WRITE(6,'(A,4I4)') ' I_AB, I_AA, I_CB, I_CA',
C?   &                                I_AB, I_AA, I_CB, I_CA
C?               WRITE(6,*) ' IT and Diag(IT)', IT,CCDIA(IT)
                END DO
               END DO
              END DO
             END DO
*            ^ End of loop over elements of block
  777       CONTINUE
            END DO
*           ^ End of loop over ISM_AA
         END DO
*        ^ End of loop over ISM_CA
       END DO
*      ^ End of loop over ISM_C
      END DO
*     ^ End of loop over ITSS
*
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' <0![F,T]!0> constructed '
        WRITE(6,*) ' Number of elements ', IT
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' The real diagonal   '
        CALL WRTMAT(CCDIA,1,IT,1,IT)
        IF(IRECOM.NE.1) THEN
        WRITE(6,*) ' The img diagonal   '
        CALL WRTMAT(CCDIA(IT+1),1,IT,1,IT)
        END IF
      END IF
*
      RETURN
      END
*
* Get diagonal one-electron integrals in SIGDEN formalism.
*
      subroutine geth1dia_sigden_KRCC(IUB,XOUT,XINT,
     &                           MDIASV,IPOSDIA,N1E,N1DIAINT,NUMINT,
     &                           ISTART,IEND,IREOTS,IPRNT)
*
#include "implicit.inc"
#include "ipoist8.inc"
*
      dimension IREOTS(*),MDIASV(4,*),IPOSDIA(*)
      dimension XINT(N1E),XOUT(*)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      IZERO = 0
      ZERO = 0.D0
      ICT = 1
C      do II=1,10
C        print*,'IPOSDIA(II),XINT(IPOSDIA(II))',
C     &  IPOSDIA(II),XINT(IPOSDIA(II))
C      end do
      do ISEA=ISTART,IEND,1
*
* Rename to SYMM-TYPE ordering on external list:
        ISEA_REO = IREOTS(ISEA)
C       do INT=1,N1E,1
        do INT=1,N1DIAINT,1
C         if (INDEX(IUB,INT).eq.ISEA_REO.and.INDEX(IUB+2,INT).eq.
C    &        ISEA_REO.and.INDEX(IUB+1,INT).eq.IZERO) then
          if (MDIASV(IUB,INT).eq.ISEA_REO.and.MDIASV(IUB+2,INT).eq.
     &        ISEA_REO.and.MDIASV(IUB+1,INT).eq.IZERO) then
            if (NTEST.ge.100) then
C             write(6,'(A,I5,A,I5,A,1F18.12)')
C    &                'Writing integral no. ',INT,'to ',ICT,
C    &                '  Integral ',XINT(INT)
              write(6,'(A,I5,A,I5,A,1F18.12)')
     &                'Writing integral no. ',IPOSDIA(INT),'to ',ICT,
     &                '  Integral ',XINT(IPOSDIA(INT))
C            print*,'XINT(IPOSDIA(INT)',XINT(IPOSDIA(INT))
            end if
C           XOUT(ICT) = XINT(INT)
            XOUT(ICT) = XINT(IPOSDIA(INT))
            ICT = ICT + 1
            goto 10
          end if
        end do
        write(6,*) 'No integral found in geth1dia_sigden.'
        write(6,*) 'Index pair ',ISEA_REO,ISEA_REO
        stop 'Quitting.'
  10    continue
      end do
      if ((ICT-1).ne.NUMINT) then
        write(6,*) 'Number of processed integrals: ',ICT-1
        write(6,*) 'Number of diagonal elements  : ',NUMINT
        write(6,*) 'Problem in geth1dia_sigden.'
        stop 'Quitting.'
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Get diagonal two-electron integrals in SIGDEN formalism.
*   Timo Fleig, August 2001
* Modified to account for integrals over true spinors, i.e.
* explicit reference to unbarred (+) and barred (-) indices (IKFAC).
* Revised routine.
*   Timo Fleig, February 2002
*
      subroutine geth2dia_sigden_KRCC(XOUT,XINT,ISEA,JSEA,TYPE,
     &                           IOFF2,N12,I2DIATP,IKFAC,IKPERM,
     &                           IREOTS,IREOTS2,
     &                           MDIASV,IPOSDIA,N1DIAINT,NDIAINT,IPRNT)
*
#include "implicit.inc"
#include "ipoist8.inc"
*
      character*4 TYPE
*
      dimension IAC(2),JAC(2),IOFF2(*),IKPERM(4),
     &          IREOTS(*),IREOTS2(*),MDIASV(4,*),IPOSDIA(*)
      dimension XINT(N12)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      IERROR = 1
      N2B = IOFF2(I2DIATP)
C     do INT = N2B,N12,1
      do INT = N1DIAINT+1,NDIAINT,1
        if (TYPE.eq.'COUL') then
C         I = INDEX(IKPERM(1),INT)
C         J = INDEX(IKPERM(2),INT)
C         K = INDEX(IKPERM(3),INT) * IKFAC
C         L = INDEX(IKPERM(4),INT) * IKFAC
          I = MDIASV(IKPERM(1),INT)
          J = MDIASV(IKPERM(2),INT)
          K = MDIASV(IKPERM(3),INT) * IKFAC
          L = MDIASV(IKPERM(4),INT) * IKFAC
        else if (TYPE.eq.'EXC1') then
C         I = INDEX(IKPERM(1),INT)
C         J = INDEX(IKPERM(4),INT) * IKFAC
C         K = INDEX(IKPERM(2),INT)
C         L = INDEX(IKPERM(3),INT) * IKFAC
          I = MDIASV(IKPERM(1),INT)
          J = MDIASV(IKPERM(4),INT) * IKFAC
          K = MDIASV(IKPERM(2),INT)
          L = MDIASV(IKPERM(3),INT) * IKFAC
        else if (TYPE.eq.'EXC2') then
C         I = INDEX(IKPERM(1),INT)
C         J = INDEX(IKPERM(4),INT)
C         K = INDEX(IKPERM(2),INT) * IKFAC
C         L = INDEX(IKPERM(3),INT) * IKFAC
          I = MDIASV(IKPERM(1),INT)
          J = MDIASV(IKPERM(4),INT)
          K = MDIASV(IKPERM(2),INT) * IKFAC
          L = MDIASV(IKPERM(3),INT) * IKFAC
        end if
*
* Rename: External ST ordering, internal TS ordering:
        if (I.eq.IREOTS(ISEA).and.
     &      J.eq.IREOTS(ISEA).and.
     &      K.eq.IREOTS2(JSEA).and.
     &      L.eq.IREOTS2(JSEA)) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
C            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'Diagonal int type ',TYPE,
     &             IPOSDIA(INT),(MDIASV(II,IPOSDIA(INT)),II=1,4),
     &             XINT(IPOSDIA(INT))
            write(6,*)
          end if
C         XOUT = XINT(INT)
          XOUT = XINT(IPOSDIA(INT))
          IERROR = 0
* and exit.
          goto 10
        end if
      end do
      if (IERROR.eq.1) then
        write(6,*) 'No integral found in geth2dia_sigden.'
C       write(6,*) 'N2B,N12 = ',N2B,N12
        write(6,*) 'Type of integral = ',TYPE
        write(6,*) 'Index pair ',IREOTS(ISEA),IREOTS2(JSEA)
        write(6,*) 'Indices ',I,J,K,L
        stop 'Quitting.'
      end if
*
  10  continue
      return
      end
*
