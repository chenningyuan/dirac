********************************************************************
*                                                                  *
* LUCIA, by Jeppe Olsen,                                           *
* DIRAC adaptation to ARDUCCA by Timo Fleig, Lasse Sorensen        *
*                                                                  *
********************************************************************
*
      SUBROUTINE FI_KRCC(T,FIMAT,ECC,IDOH2,ISPINFREE)
*
*. Construct inactive fockmatrix + core-core interaction energy.
*. I.e. add contributions from all orbitals
*  that belong to hole orbital spaces ( as defined by IPHGAS).
*
* Note that this is a more general definition of the
* Inactive Fockmatrix than usually used.
*
*. On input FIMAT should be the inactive Fock matrix, in symmetry packed form
*
* If I_USE_SIMTRH = 0 input and output matrices are assumed lower half packed
*                 = 1 Input and output matrices are assumed complete blocks
*
* Jeppe Olsen
*
* Revision : Dec 97 : General hole spaces
*            aug 00 : I_USE_SIMTRH switch added
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION T(*),FIMAT(*)
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cgas.inc"
#include "strinp.inc"
*
*
      I_USE_SIMTRH = 0
      IF(I_USE_SIMTRH.EQ.0) THEN
        CALL FIH_KRCC(T,FIMAT,ECC,IBSO,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
     &               NGAS,NGSSH,IHPVGAS_AB,
     &               NTOOBS,NTOOB,IDOH2,
     &               ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
     &               NTOOB,NELEC)
      ELSE
Commented out by lasse 2010
C       CALL FIHA_KRCC(T,FIMAT,ECC,IBSO,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
C    &                NTOOBS,NTOOB,IDOH2,
C    &                ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
C    &                NTOOB)
      END IF
*
      RETURN
      END
*
C      SUBROUTINE FI_HS_REL(T,FI_AL,FI_BE,ECC,IDOH2)
      SUBROUTINE FIH_KRCC(T,FI,ECC,IOBSM,NSMOB,ITPFSO,ITPFSO2,IPHGAS,
     &                   NGAS,NGSSH,IHPVGAS_AB,
     &                   LOBSM,NORBT,IDOH2,
     &                   ISPINFREE,NGSOB,NGSOB2,NGSOBT,NGSOBT2,
     &                   NTOOB,NELEC)
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
*
* Changed to a structure similar to that in fih_rel
*
* Added more efficient integral fetching and fixed bug
* Work only for real groups at the moment since it is
* exploting the integral symmtry of these along with
* time reversal symmetry.
*                                  Lasse june 2008
*
* Non rel case
*
* F(n,m) = h(n,m) + SUM(i) ( 2g(n,m,i,i) - g(n,i,i,m))
*
* Rel case (with timereversal symmtry). Notice i is orbital/kramers pair index
*
* F(n,m) = h(n,m) + SUM(i) ( 2g(n,m,i,i) - g(n,i,i,m) - g(n,i_,i_,m))
*
* Since integrals are not stored like this we will rewrite this
*
* F(n,m) = h(n,m) + SUM(i) ( (g(n,m,i,i) - g(n,i,i,m)) + (g(n,m,i_,i_) - g(n,i_,i_,m)) )
*
* Where we have NALPHA = 0 for first term and NALPHA = 2 for second
*
* This in turn will be fetched like this
*
* F(n,m) = h(n,m) + SUM(2*i) ( 1/2 * (g(n,m,i,i) - g(n,i,i,m)) + (g(n,m,i_,i_) - g(n,i_,i_,m)) )
*
* since we will sum over spinors where we may not have same occupation for the Fermi vacuum
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "ctcc.inc"
#include "gasstr.inc"
*. Standard 1+2 electron integrals
      DIMENSION T(*)
*. Output
      DIMENSION FI(*)
*. Input
      INTEGER IOBSM(*),LOBSM(*)
      INTEGER NGSSH(MXPIRR,NGAS), NELEC(*)
      INTEGER ITPFSO(*), ITPFSO2(*), IPHGAS(*)
      INTEGER NGSOB(MXPOBS,*), NGSOB2(MXPOBS,*)
      INTEGER NGSOBT(*), NGSOBT2(*), IHPVGAS_AB(MXPNGAS,2)
*
      NTEST = 00
C      print*,'NTOOBS,NTOOBS,NSMOB',NTOOBS,NTOOBS,NSMOB
C      print*,'NORBT',NORBT
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' FI_HS : Initial 1 e integral list '
C       CALL APRBLM2(T,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
        CALL APRBLM2(T,NORBT,NORBT,NSMOB,ISYMPACK)
      END IF
*
      ZERO = 0.0D0
      HALF = 0.5D0
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
      CALL SETVEC(FI,ZERO,IINT)
*
* Core energy
*
      ECC = 0.0D0
      IJSM = 1
      IIOFF = 0
      NINT12 = N1ELINT + N2ELINT
*
* GAS-SYMM ordering (GAS is outer loop!)
      do IGAS=1,NGAS,1
        do JGAS=1,NGAS,1
          do ISMOB=1,NSMOB/2,1
            do JSMOB=1,NSMOB/2,1
              ISTART = 0
              JSTART = 0
* Determine index range for this symmetry/GAS block:
          ININT = NGSSH(ISMOB,IGAS) * NGSSH(JSMOB,JGAS)
          DO IGASS =1,IGAS-1,1
            DO ISM=1,NSMOB/2,1
              ISTART = ISTART + NGSSH(ISM,IGASS)
            END DO
          END DO
C          DO ISM = 1,ISMOB-1,1
C            ISTART = ISTART + NGSSH(ISM,IGAS)
C          END DO
          ISTART = ISTART + 1
          IEND = ISTART + NGSSH(ISMOB,IGAS) - 1
          DO JGASS =1,JGAS-1,1
            DO JSM=1,NSMOB/2,1
              JSTART = JSTART + NGSSH(JSM,JGASS)
C              print*,'NGSSH,JSM,JGASS',NGSSH(JSM,JGASS),JSM,JGASS
            END DO
          END DO
          JSTART = JSTART + 1
          JEND = JSTART + NGSSH(JSMOB,JGAS) - 1
          if (NTEST.ge.50) then
            write(6,*) 'IGAS,ISMOB,ISTART,IEND',
     &                  IGAS,ISMOB,ISTART,IEND
            write(6,*) 'ININT,NGSSH(ISMOB,IGAS)',
     &                  ININT,NGSSH(ISMOB,IGAS)
            write(6,*) 'JGAS,JSMOB,JSTART,JEND',
     &                  JGAS,JSMOB,JSTART,JEND
            write(6,*) 'ININT,NGSSH(JSMOB,JGAS)',
     &                  ININT,NGSSH(JSMOB,JGAS)
          end if
          ICOUNT = ITPFSO(ISTART)
            IUB = 1
            IF(ININT.GE.1) THEN
            CALL ONE_ELEC_INACTIV(IUB,FI,T,
     &                            ISTART,IEND,JSTART,JEND,ISIDE)
            END IF
            END DO
          END DO
        end do
      end do
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Printing entire one-electron inactive Fock-matrix '
        DO I = 1, ISIDE
          DO J = 1, ISIDE
            WRITE(6,*) FI(I + (J-1)*ISIDE) , I, J
          END DO
        END DO
      END IF
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
C          print*,'NGASOCC,IGS',NGASOCC(IGS),IGS
          IF(NGASOCC(IGS).eq.1) THEN
            if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
            do JGS = 1,NGAS,1
C          print*,'NGASOCC(JGS),JGS',NGASOCC(JGS),JGS
              IF(NGASOCC(JGS).eq.1) THEN
                if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT(JGS-1)
C                do KGS = 1,NGAS,1
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
C                IF(IHPVGAS_AB(ITPFSO(IOFF+1),1).EQ.1) THEN
C     &             IHPVGAS_AB(ITPFSO(JOFF+1),1).EQ.1) THEN
                do I = 1+IOFF,NI+IOFF,1
                  do J = 1+JOFF,NJ+JOFF,1
C
C Calculate position in matrix. Complete element calculated in one loop.
C
                    IPOS = I + (J-1)*ISIDE
C
C Assuming a bottom up occupation of electrons
C IUB loop to only count occpied contributions for IAB
                    DO IUB = 1,2,1
                      DO IAB = 1,NELEC(IUB)
C                     print*,'I,J,IAB,IPOS',I,J,IAB,IPOS
C
C Here a call to get int with nalpha = 0 will fetch first term in sum
C and nalpha = 2 second part ( which is the one with the special exchange 2 integrals)
C
                      NOP = 4
                      XOUT = 0.0D0
                      XOUT2 = 0.0D0
                      II = I
                      JJ = J
                      IABA = IAB
                      IABB = IAB
C
C Check if integral for nalpha = 0 or 4 is present and fetch it.
C
                      INOINT = 0
                      NALPHA = 4
                      IF(I.EQ.IAB.OR.
     &                   J.EQ.IAB) THEN
                        INOINT = 1
                      END IF
                      IF(NOP.EQ.4.AND.INOINT.EQ.0) THEN
                      CALL GET_INT2(XOUT,T,II,JJ,
     &                                 IABA,IABB,
     &                                 NALPHA)
                      END IF
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
                      CALL GET_INT2(XOUT2,T,II,JJ,
     &                                 IABA,IABB,
     &                                 NALPHA)
C
C Add the integrals to inactive fock matrix
C
                      FI(IPOS) = HALF * (XOUT + XOUT2) + FI(IPOS)
C                     print*,'XOUT,XOUT2,IPOS',XOUT,XOUT2,IPOS
                      END DO
                    END DO
                  end do
                end do
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
C
       NTEST =00
       IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Printing entire inactive Fock-matrix '
         DO I = 1, ISIDE
           DO J = 1, ISIDE
             WRITE(6,*) FI(I + (J-1)*ISIDE) , I, J
C            IF(I.GE.J) WRITE(6,*) FI(I + (J-1)*ISIDE) , I, J
           END DO
         END DO
       END IF
*
*
      RETURN
      END
C
***********************************************************************
C
      SUBROUTINE ONE_ELEC_INACTIV(IUB,XOUT,XINT,
     &                            ISTART,IEND,JSTART,JEND,ISIDE)
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
      NTEST = 00
*
      IPOS = 0
      IORBCOUNT = 0
      ILENGHT = 0
*
C Map first indix to orbital subspace
      DO IGAS = 1,NGAS
        IORBCOUNT = IORBCOUNT + NGSOBT(IGAS)
        IF(IEND.EQ.IORBCOUNT) THEN
          I_INDX_TO_GAS(1) = IGAS
        END IF
        IF(JEND.EQ.IORBCOUNT) THEN
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
        ITYPE = ITYPE + IRELTYPE(2)
      END IF
C Calculate INTEGRAL offset for absolute type TYPE
      IOFF = IOFFINTTYPE(ITYPE)
C      print*,'ITYPE',ITYPE
C      print*,'IOFF',IOFF
      do JSEA=JSTART,JEND,1
        do ISEA=ISTART,IEND,1
C
C Calculate position of integral in matrix
C
          IPOS = ISEA + (JSEA-1)*ISIDE
C
C Calculate absolute offset off integral with help from relative
C
          XOUT(IPOS) = XINT(IOFF)
C          print*,'IPOS',IPOS
C          print*,'XOUT(IPOS)',XOUT(IPOS)
          IOFF = IOFF + 1
        end do
      end do
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'DUMPING ONE ELECTRON INTEGRALS FETCHED'
        DO I = 1,ISIDE**2
          WRITE(6,*) XOUT(I)
        END DO
      END IF
*
      RETURN
      END
