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
#include "ipoist8.inc"
      DIMENSION T(*)
#include "mxpdim.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "cintfo.inc"
#include "orbinp.inc"
*Added by Lasse
#include "strinp.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "crun.inc"
*
      DIMENSION VEC1(*),VEC2(*)
*
      DIMENSION WORK(*)
*
      IDUM = 0

      CALL MEMGET('REAL',KFDIA_AL,NTOOB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFDIA_BE,NTOOB,WORK,KFREE,LFREE)
      IRECOM = 1
      IF(IRECOM.NE.1) THEN
      CALL MEMGET('REAL',KFDIA_AL_IM,NTOOB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFDIA_BE_IM,NTOOB,WORK,KFREE,LFREE)
      END IF
*. Four blocks of string occupations
C
C     MX_ST_TSOSO_BLK_MX NEEDS TO BE CORRECTED FOR LARGER CC CALCULATIONS
C
      CALL MEMGET('INTE',KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
*. Four arrays ( will be used to store sums of orbital energies)
      CALL MEMGET('REAL',KLSTR1,MX_ST_TSOSO_MX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLSTR2,MX_ST_TSOSO_MX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLSTR3,MX_ST_TSOSO_MX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLSTR4,MX_ST_TSOSO_MX,WORK,KFREE,LFREE)
      IF(IRECOM.NE.1) THEN
*. And for the imaginary part
      CALL MEMGET('REAL',KLSTR5,MX_ST_TSOSO_MX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLSTR6,MX_ST_TSOSO_MX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLSTR7,MX_ST_TSOSO_MX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLSTR8,MX_ST_TSOSO_MX,WORK,KFREE,LFREE)
      END IF
*
* 1 : Obtain Inactive Fock matrix
*
*. One-body density matrix
      LBLK = -1
C      CALL COPVCDC(LUC,LUSC1,VEC1,1,LBLK)
C      CALL WRTVCDC(VEC1,LUC,1,LBLK)
C      CALL SIGDEN_CTRL(WORK(KVEC1),WORK(KVEC2),LUC,LUHC,WORK(KT_CC),2)
C      CALL SIGDEN_CTRL(VEC1,VEC2,LUC,LUSC1,WORK(KT_CC),2)
      CALL REWINO(LUC)
      CALL REWINO(LUHC)
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
*****
* Notice harwirering of I_RES_OR_UNRES . Lasse
*****
      I_RES_OR_UNRES = 2
      IF(I_RES_OR_UNRES.EQ.1) THEN
*. Standard spinrestricted Fock matrices
*. Extract diagonal of fock-matrix
      ELSE
C NOTE ONLY DIAGONAL IS USED!!!
C THEREFORE ONLY THAT WILL BE GENERATED
        I_ORI = 0
        IF(I_ORI.EQ.0) THEN
C much easier just to fetch the diagonal from inactive Fock matrix hence
          IOFF = 0
          CALL GET_DIA_KRCC(WORK(KFDIA_AL),WORK(KFDIA_BE),
     &                      WORK(KFI),NTOOB,IOFF)
C         CALL FI_HS_DIA_KRCC(T,ECC,
C    &               IBSO,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
C    &               NGAS,NGSSH,IHPVGAS_AB,
C    &               NTOOBS,NTOOB,IDOH2,
C    &               ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
C    &               NTOOB,WORK(KFDIA_AL),WORK(KFDIA_BE),NELEC)
*
          IF(IRECOM.NE.1) THEN
            print*,'img dia'
          IOFF = NTOOB
          CALL GET_DIA_KRCC(WORK(KFDIA_AL_IM),WORK(KFDIA_BE_IM),
     &                      WORK(KFI),NTOOB,IOFF)
C           CALL FI_HS_DIA_KRCC(T(N1ELINT+N2ELINT+1),ECCCOM,
C    &               IBSO,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
C    &               NGAS,NGSSH,IHPVGAS_AB,
C    &               NTOOBS,NTOOB,IDOH2,
C    &               ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
C    &               NTOOB,WORK(KFDIA_AL_IM),WORK(KFDIA_BE_IM),NELEC)
          END IF
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
C
C New addition by Lasse. Calculate the energy
C
      CALL FOCK_ENERGY(ESOMETHING,T,NTOOB,WORK(KFDIA_AL),WORK(KFDIA_BE))
*
C?    WRITE(6,*) ' Enforced stop in GENCC_F_DIAG_M '
C?    STOP '  Enforced stop in GENCC_F_DIAG_M'
      CALL MEMREL('F_DIAG',WORK,KFDIA_AL,KFDIA_AL,KFREE,LFREE)
*
      RETURN
      END
*
      SUBROUTINE FOCK_ENERGY(ECC,T,NTOOB,DIA_AL,DIA_BL)
*
* Should be made 
*
#include "implicit.inc"
      DIMENSION T(*),DIA_AL(NTOOB),DIA_BL(NTOOB)
      INTEGER IDUM(2)
*
      ECC = 0.0D0
*
      IDUM(1) = 1
      IDUM(2) = 4
      DO I=1,2
        ECC = T(IDUM(I)) + DIA_AL(I) + ECC
        print*,'T(IDUM(I)),DIA_AL(I),DIA_BL(I),I,ECC',
     &          T(IDUM(I)),DIA_AL(I),DIA_BL(I),I,ECC
      END DO
      print*,'ECC',ECC
*
      RETURN
      END
*
      SUBROUTINE GET_DIA_KRCC(DIA_AL,DIA_BE,
     &                      FOCK,NTOOB,IOFF)
*
* The easy way to extract the diagonal from the inactive fock matrix
* At the moment only double occupied spinors.
*
#include "implicit.inc"
*
      DIMENSION DIA_AL(NTOOB),DIA_BE(NTOOB),FOCK(NTOOB**2)
*
      DO I=1,NTOOB
        IFETCH = IOFF + NTOOB*(I-1) + I
C       print*,'IFETCH',IFETCH
        DIA_AL(I) = FOCK(IFETCH)
        DIA_BE(I) = FOCK(IFETCH)
        print*,'DIA_AL(I),I',DIA_AL(I),I
      END DO
*
      RETURN
      END
*
      SUBROUTINE FI_HS_DIA_KRCC(T,ECC,IOBSM,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
     &                   NGAS,NGSSH,IHPVGAS_AB,
C    &                   LOBSM,NORBT,IDOH2,
     &                   NTOOBS,NORBT,IDOH2,
     &                   ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
     &                   NTOOB,FI_AL,FI_BE,NELEC)
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
* Made the final changes so we could get rid og the arrays with
* diagonal integrals stored in them.
*
* Lasse, 2011
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "ctcc.inc"
#include "gasstr.inc"
*
*. Standard 1+2 electron integrals
      DIMENSION T(*)
*. Output
      DIMENSION FI_AL(*), FI_BE(*)
*. Input
C     INTEGER IOBSM(*),LOBSM(*)
      INTEGER IOBSM(*),NTOOBS(*)
      INTEGER NGSSH(MXPIRR,NGAS), NELEC(*)
      INTEGER ITPFSO(*), ITPFSO2(*), IPHGAS(*), IKPERM(4)
      INTEGER NGSOB(MXPOBS,*), NGSOB2(MXPOBS,*)
      INTEGER NGSOBT(*), NGSOBT2(*), IHPVGAS_AB(MXPNGAS,2)
*
* Is matrix symmetrically packed or ???
      ISYMPACK = 1
      ZERO = 0.0D0
      HALF = 0.5D0
      CALL SETVEC(FI_AL,ZERO,NORBT)
      CALL SETVEC(FI_BE,ZERO,NORBT)
*
      NTEST = 00
      print*,'ISYMPACK,NSMOB',ISYMPACK,NSMOB
      print*,'NORBT',NORBT
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
C IINT is dimension of FI, ISIDE is square root of IINT
      IINT = 0
      ISIDE = 0
      DO IGAS = 1,NGAS
        DO ISMOB = 1,NSMOB
          DO JGAS = 1,NGAS
            DO JSMOB = 1,NSMOB
              IINT = IINT + NGSSH(ISMOB,IGAS) * NGSSH(JSMOB,JGAS)
            END DO
          END DO
          IF(NGSSH(ISMOB,IGAS).NE.0) THEN
            ISIDE = ISIDE + NGSSH(ISMOB,IGAS)
          END IF
        END DO
      END DO
*
      DO I=1,NORBT
        FI_AL(I) = 0.0D0
        FI_BE(I) = 0.0D0
      END DO
*
* GAS-SYMM ordering (GAS is outer loop!)
      ICOUNT = 1
*
      do IGAS=1,NGAS,1
        do ISMOB=1,NSMOB/2,1
          ISTART = 0
* Determine index range for this symmetry/GAS block:
          ININT = NGSSH(ISMOB,IGAS) * NGSSH(ISMOB,IGAS)
          DO IGASS =1,IGAS-1,1
            DO ISM=1,NSMOB/2,1
              ISTART = ISTART + NGSSH(ISM,IGASS)
            END DO
          END DO
          ISTART = ISTART + 1
          IEND = ISTART + NGSSH(ISMOB,IGAS) - 1
          if (NTEST.ge.50) then
            write(6,*) 'IGAS,ISMOB,ISTART,IEND',
     &                  IGAS,ISMOB,ISTART,IEND
            write(6,*) 'ININT,NGSSH(ISMOB,IGAS)',
     &                  ININT,NGSSH(ISMOB,IGAS)
          end if
C         ICOUNT = ITPFSO(ISTART)
          print*,'ICOUNT',ICOUNT
          IUB = 1
          IF(ININT.GE.1) THEN
            IUB = 1
            CALL ONE_ELEC_DIA(IUB,FI_AL(ICOUNT),T,ISTART,IEND,ISIDE)
            IF(IHPVGAS_AB(ICOUNT,1).EQ.1) THEN
              DO II=1,NGSSH(ISMOB,IGAS)
                ECC = ECC + FI_AL(II)
              END DO
            END IF
            IUB = 2
            CALL ONE_ELEC_DIA(IUB,FI_BE(ICOUNT),T,ISTART,IEND,ISIDE)
            IF(IHPVGAS_AB(ICOUNT,2).EQ.1) THEN
              DO II=1,NGSSH(ISMOB,IGAS)
                ECC = ECC + FI_BE(II)
              END DO
            END IF
          END IF
          ICOUNT = ICOUNT + IEND - ISTART + 1
        END DO
      end do
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Printing Diagonal 1e-integrals: alpha'
        DO I = 1, ISIDE
          WRITE(6,*) FI_AL(I),I
        END DO
        WRITE(6,*) ' Printing Diagonal 1e-integrals: beta'
        DO I = 1, ISIDE
          WRITE(6,*) FI_BE(I),I
        END DO
      END IF
      print*,'grep for this'
      print*,'ECC',ECC
*
*
* 2-electron case. Now fetching all integrals in one loop.
* Fetching is done similar to GET_OPINT4_REL.
*
        ICOUNT = 1
        IGSOFF = 0
        JGSOFF = 0
        ISMOFF = 0
        JSMOFF = 0
        IJ = 1
        do IGS = 1,NGAS,1
          print*,'NGASOCC,IGS',NGASOCC(IGS),IGS
          IF(NGASOCC(IGS).eq.1) THEN
            if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
            do JGS = IGS,NGAS,1
          print*,'NGASOCC(JGS),JGS',NGASOCC(JGS),JGS
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
C
                IF(IHPVGAS_AB(ITPFSO(IOFF+1),1).EQ.1) THEN
C
                print*,'IOFF,,NI,JOFF,NJ',IOFF,NI,JOFF,NJ
                do I = 1+IOFF,NI+IOFF,1
                 print*,'I',I
                  do J = 1+JOFF,NJ+JOFF,1
                   print*,'J',J
C                  IF(J.GE.3) CYCLE
C
C Assuming a bottom up occupation of electrons
C IUB loop to only count occpied contributions for IAB
C                   DO IUB = 1,2,1
C                     DO IAB = 1,NELEC(IUB)
C                     print*,'I,J,IAB,IPOS',I,J,IAB,IPOS
C
C Here a call to get int with nalpha = 0 will fetch first term in sum
C and nalpha = 2 second part ( which is the one with the special exchange 2 integrals)
C
                      NOP = 4
                      XOUT = 0.0D0
                      XOUT2 = 0.0D0
                      XOUT3 = 0.0D0
                      XOUT4 = 0.0D0
                      IF(I.GT.J) THEN
                      II = I
                      JJ = J
                      ELSE
                      II = J
                      JJ = I
                      END IF
C                     IABA = IAB
C                     IABB = IAB
C
C Check if integral for nalpha = 0 or 4 is present and fetch it.
C
C                     DO III=1,100
C                       print*,'Lasse is not sure if this is right'
C                       print*,'please debug completely'
C                     END DO
                      print*,'II,JJ',II,JJ
                      IF(IHPVGAS_AB(ITPFSO(JOFF+1),1).EQ.1) THEN
                      INOINT = 0
                      NALPHA = 4
                      IF(I.EQ.J) THEN
                        INOINT = 1
                      END IF
                      IF(NOP.EQ.4.AND.INOINT.EQ.0) THEN
C                     CALL GET_INT2(XOUT,T,II,II,JJ,JJ,NALPHA)
C                     II = I
C                     JJ = J
C                     CALL GET_INT2(XOUT,T,II,JJ,II,JJ,NALPHA)
C                     II = I
C                     JJ = J
C                     CALL GET_INT2(XOUT,T,II,JJ,JJ,II,NALPHA)
C                     II = I
C                     JJ = J
C                     CALL GET_INT2(XOUT,T,JJ,II,II,JJ,NALPHA)
C                     II = I
C                     JJ = J
C                     CALL GET_INT2(XOUT,T,JJ,II,JJ,II,NALPHA)
C                     II = I
C                     JJ = J
                      CALL GET_INT2(XOUT,T,JJ,JJ,II,II,NALPHA)
                      END IF
                      END IF
                      print*,'XOUT',XOUT
C
C Fetch second term with nalpha = 2
C
C                     print*,'XOUT',XOUT
                      INOINT = 0
                      NALPHA = 2
                      II = I
                      JJ = J
                      IABA = IAB
                      IABB = IAB
C                     print*,'II,JJ',II,JJ
C                     CALL GET_INT2(XOUT2,T,II,II,JJ,JJ,NALPHA)
C                     CALL GET_INT2(XOUT2,T,JJ,JJ,II,II,NALPHA)
C                     CALL GET_INT2(XOUT2,T,JJ,II,JJ,II,NALPHA)

C                     CALL GET_INT2(XOUT2,T,II,II,JJ,JJ,NALPHA)
C                     II = I
C                     JJ = J
C                     CALL GET_INT2(XOUT2,T,II,JJ,II,JJ,NALPHA)
C                     II = I
C                     JJ = J
C                     CALL GET_INT2(XOUT2,T,II,JJ,JJ,II,NALPHA)
C                     II = I
C                     JJ = J
C                     CALL GET_INT2(XOUT2,T,JJ,II,II,JJ,NALPHA)
C                     II = I
C                     JJ = J
C                     CALL GET_INT2(XOUT2,T,JJ,II,JJ,II,NALPHA)
C                     II = I
C                     JJ = J
                      CALL GET_INT2(XOUT2,T,JJ,JJ,II,II,NALPHA)
                      print*,'XOUT2',XOUT2
C                     IF(II.NE.JJ) THEN
C                     print*,'II,JJ',II,JJ
C                       CALL GET_INT2(XOUT3,T,JJ,JJ,II,II,NALPHA)
C                       CALL GET_INT2(XOUT4,T,II,JJ,JJ,II,NALPHA)
C                       XOUT2 = XOUT2 + XOUT3 - XOUT4
C                       print*,'XOUT3,XOUT4',XOUT3,XOUT4
C                     END IF
C
C Add the integrals to inactive fock matrix
C
                      FI_AL(J) = (XOUT + XOUT2) + FI_AL(J)
                      FI_BE(J) = (XOUT + XOUT2) + FI_BE(J)
C                     END DO
C                   END DO
                  end do
                end do
                END IF
          ININT = NGSSH(ISMOB,IGAS) * NGSSH(JSMOB,JGAS)
              end do
              JSMOFF = 0
            end do
            ISMOFF = 0
            END IF
            end do
C           ^ JGAS
          JGSOFF = 0
          ELSE IF(NGASOCC(IGS).eq.0.and.IGS.eq.1) THEN
            DO I=1,NTOOB
              DO J=1,NTOOB
                REPLACED = ZERO
              END DO
            END DO
          END IF
C         ^ NGASOCC
        end do
C       ^ IGAS
       NTEST = 100
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Printing diagonal 2e-case : alpha '
         DO I = 1, ISIDE
           WRITE(6,*) FI_AL(I), I
         END DO
         WRITE(6,*) ' Printing diagonal 2e-case : beta '
         DO I = 1, ISIDE
           WRITE(6,*) FI_BE(I), I
         END DO
       END IF
*
*
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
      stop
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
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "orbinp.inc"
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
      NTEST = 000
      IT = 0
C     print*,'NTSS_TP',NTSS_TP
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
C      print*,'NEL_CA,NEL_CB,NEL_AA,NEL_AB',
C    &         NEL_CA,NEL_CB,NEL_AA,NEL_AB
*
* find these in gnside_rel (LUCIAREL) in sigden_r.F
*
*. Diagonal block ?
       CALL DIAG_EXC_KRCC(ITSS_TP(1+0*NGAS,ITSS),
     &                  ITSS_TP(1+1*NGAS,ITSS),
     &                  ITSS_TP(1+2*NGAS,ITSS),
     &                  ITSS_TP(1+3*NGAS,ITSS),NGAS,IDIAG)
C      IF(IDIAG.EQ.0) THEN
         IRESTRICT = 0
C      ELSE
C        IRESTRICT = 1
C      END IF
       DO ISM_C = 1, NSMST
Cori     ISM_A = MULTD2H(ISM,ISM_C)
C         ISM_A = IDBGMULT(ISM_C,IADJSYM(ISM))
C         ISM_A = IDBGMULT(ISM,IADJSYM(ISM))
C         ISM_A = IDBGMULT(ISM,IADJSYM(ISM_C))
C        ISM_A = IDBGMULT(ISM,INVELM(ISM_C))
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
C           WRITE(6,'(A,4I5)') ' ISM_CA, ISM_CB, ISM_AA, ISM_AB',
C    &                           ISM_CA, ISM_CB, ISM_AA, ISM_AB
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
C            print*,'NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA',
C    &               NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA
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
C                print*,'ITSS',ITSS
C                print*,'E_CA(I_CA),E_CB(I_CB),E_AA(I_AA),E_AB(I_AB)'
C                print*, E_CA(I_CA),E_CB(I_CB),E_AA(I_AA),E_AB(I_AB)
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
C
      SUBROUTINE ONE_ELEC_DIA(IUB,XOUT,XINT,ISTART,IEND,ISIDE)
C
C New routine to fetch one electron integrals from KT_CC
C This will fetch a block of a type and sort it to a matrix form.
C IUB not actively used but can be if there is no Kramers symmetry
C of our one particle functions.
C                                                         Lasse
C
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
*
      DIMENSION XINT(*),XOUT(ISIDE**2)
      INTEGER I_INDX_TO_GAS(2),ISUBLENGHT(NGAS+1),IORB(NGAS+1)
*
      NTEST = 100
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Start index = ',ISTART
        WRITE(6,*) ' End index = ',IEND
        WRITE(6,*) ' Dimension of side = ',ISIDE
      END IF
*
      IPOS = 0
      IORBCOUNT = 0
      ILENGHT = 0
      IPOS = 0
*
C Map first indix to orbital subspace
      DO IGAS = 1,NGAS
        IORBCOUNT = IORBCOUNT + NGSOBT(IGAS)
        IF(IEND.EQ.IORBCOUNT) THEN
          I_INDX_TO_GAS(1) = IGAS
          I_INDX_TO_GAS(2) = IGAS
        END IF
        ILENGHT = ILENGHT + IGAS
        ISUBLENGHT(IGAS+1) =  ILENGHT
        IORB(IGAS+1) = IORBCOUNT
      END DO
*
C Calculate relative TYPE
      ITYPE = I_INDX_TO_GAS(2) + NGAS*(I_INDX_TO_GAS(1) - 1)
C Absolute type
      IF(IUB.EQ.1) THEN
        ITYPE = ITYPE
      ELSE
        ITYPE = ITYPE + IRELTYPE(2) - 1 
      END IF
C Calculate INTEGRAL offset for absolute type TYPE
      IOFF = IOFFINTTYPE(ITYPE)
      ILL = IEND - ISTART + 1
C
      DO ISEA=1,ILL
C
C Calculate position of integral in array
C
        IPOS = IPOS + 1 
        IFETCH = IOFF + ISEA + (ISEA-1)*ILL -1 
C
C Calculate absolute offset off integral with help from relative
C
        print*,'IPOS,IFETCH',IPOS,IFETCH
        XOUT(IPOS) = XINT(IFETCH)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'DUMPING ONE ELECTRON INTEGRALS FETCHED'
        DO I = 1,ILL
          WRITE(6,*) XOUT(I)
        END DO
      END IF
*
      RETURN
      END
