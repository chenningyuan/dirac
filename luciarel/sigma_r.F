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

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ADVICE_SIGMA_REL(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,LADVICE)
*
* Advise sigma routine about best route to take
*
* ITERM : Term  to be studied :
*         =1 alpha-beta term
*         ....... ( to be continued )
*
* LADVICE : ADVICE given ( short, an integer !!)
*
* For ITERM = 1 :
*      LADVICE = 1 : Business as usual, no transpose of matrix
*                    (resolution on alpha strings, direct exc on beta)
*      LADVICE = 2 = Transpose matrices
*                    (resolution on beta strings, direct exc on alpha)
*
* Jeppe Olsen, Tirstrup Airport, Jan 12, 98
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "gasstr.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "crun.inc"
*. Specific input
      INTEGER IAOCC(*),IBOCC(*),JAOCC(*),JBOCC(*)
*. Local Scratch
       DIMENSION ITP(16),JTP(16),KTP(16),LTP(16)
*. Temp Fix : All orbitals spaces set to particle spaces
CTI Has been set in cgas.inc
CTI    INTEGER IPHGAS(MXPNGAS)
       DO IGAS = 1, NGAS
         IPHGAS(IGAS) = 1
       END DO
*
      NTEST = 00
      IF(ITERM.EQ.1) THEN
*.
*. sigma(i,Ka,Ib) = sum(i,kl)<Ib!Eb_kl!Jb>(ij!kl)C(j,Ka,Jb)
*
* Number of ops : Number of sx(kl) N_i*N_j_dimension of C(j,Ka,Jb)
*.No absolute calc of flops is made, only a relative measure
*
* Single excitations connecting the two types
*
C            SXTYP2_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP,IPHGAS)
        CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOCC,JAOCC,IPHGAS)
        CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOCC,JBOCC,IPHGAS)
C?      WRITE(6,*) 'NIJTYP, NKLTYP', NIJTYP,NKLTYP
*. P-h modifications ( I cannot predict these at the moment
        IF(NIJTYP.GE.1.AND.NKLTYP.GE.1) THEN
*
        IF((IPHGAS(ITP(1)).EQ.2.AND.IPHGAS(JTP(1)).EQ.2).OR.
     &     (IPHGAS(KTP(1)).EQ.2.AND.IPHGAS(LTP(1)).EQ.2)     ) THEN
           IPHMODI = 1
         ELSE
           IPHMODI = 0
         END IF
        ELSE
           IPHMODI = 0
        END IF
*
        IF(IPHMODI.EQ.1.OR.NIJTYP.NE.1.OR.NKLTYP.NE.1
     &     .OR.IADVICE.EQ.0) THEN
*. Several connections, i.e. the alpha or the beta blocks are identical,
*. or ph modifications
*. just continue
          LADVICE = 1
        ELSE
* =========================================
*.. Index for flops along C(j,Ka,Jb) route
* =========================================
*.Dim of C(j,Ka,Jb) relative to C(Ja,Jb)
*. going from Ja to  Ka reduces occ by one elec, changes dim by n/(N-n+1)
          XNJOB = FLOAT(NOBPT(JTP(1)))
          XNJEL = FLOAT(JAOCC(JTP(1)))
          XCJKAJB = XNJOB*XNJEL/(XNJOB-XNJEL+1)
*. Number of kl excitations per beta string :
          XNKLSX = FLOAT((NOBPT(KTP(1))-JBOCC(KTP(1)))*JBOCC(LTP(1)))
*. Number of ops (relative to dim of C)
          XNIOB = FLOAT(NOBPT(ITP(1)))
          XFLOPA = XCJKAJB*XNKLSX*XNIOB
* =========================================
*.. Index for flops along C(l,Ja,Kb) route
* =========================================
*.Dim of C(l,Ja,Kb) relative to C(Ja,Jb)
          XNLOB = FLOAT(NOBPT(LTP(1)))
          XNLEL = FLOAT(JBOCC(LTP(1)))
          XCLJAKB = XNLOB*XNLEL/(XNLOB-XNLEL+1)
*. Number of ij excitations per alpha string :
          XNIJSX = FLOAT((NOBPT(ITP(1))-JAOCC(ITP(1)))*JAOCC(JTP(1)))
*. Number of ops (relative to dim of C)
          XNKOB = FLOAT(NOBPT(KTP(1)))
          XFLOPB = XCLJAKB*XNIJSX*XNKOB
*. Switch to second route if atleast 20 percent less work
          IF(XFLOPB.LE.0.8*XFLOPA) THEN
            LADVICE = 2
          ELSE
            LADVICE = 1
          END IF
*. Well, an additional consideration :
* If the C block involes the smallest allowed number of elecs in hole space,
* and the annihilation is in hole space
* then we do the annihilation in the space with the smallest number of
* hole electrons.
          LHOLEA =0
          LHOLEB =0
          DO IGAS = 1, NGAS
            IF(IPHGAS(IGAS).EQ.2) THEN
              LHOLEA = LHOLEA + JAOCC(IGAS)
              LHOLEB = LHOLEB + JBOCC(IGAS)
            END IF
          END DO
*
          IF(LHOLEA+LHOLEB.EQ.MNHL.AND.
     &       (IPHGAS(JTP(1)).EQ.2.OR.IPHGAS(LTP(1)).EQ.2))  THEN
             IF(LHOLEA.EQ.LHOLEB) THEN
               LLADVICE = LADVICE
             ELSE IF(LHOLEA.LT.LHOLEB) THEN
               LLADVICE= 1
             ELSE
               LLADVICE = 2
             END IF
             IF(NTEST.GE.100.AND.LADVICE.NE.LLADVICE) THEN
               WRITE(6,*) ' Advice changed by hole considetions'
               WRITE(6,*) ' LADVICE, LLADVICE', LADVICE,LLADVICE
             END IF
             LADVICE = LLADVICE
          END IF
*
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' ADVICE active '
            WRITE(6,*) ' IAOCC JAOCC IBOCC JBOCC'
            CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
            WRITE(6,*) ' ITP JTP KTP LTP ',ITP(1),JTP(1),KTP(1),LTP(1)
            WRITE(6,*) ' XFLOPA,XFLOPB', XFLOPA,XFLOPB
            WRITE(6,*) ' ADVICE given : ', LADVICE
          END IF
        END IF
*       ^ End if several types/ph modi
      END IF
*     ^ End if ITERM test ( type of excitation)
C     WRITE(6,*) ' MEMCHECK at end of ADVICE'
C     CALL RMEMCHK
C     WRITE(6,*) ' MEMCHECK passed '
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ALG_ROUTERX_REL(ILOC,IROC,NOP,IOP_TYP,IOP_AC,IOP_REO,
     &                      SIGN)
*
* Decide route for calculating <ILOC! operator string !ROC>
* with smallest amount of operations/storage
*
* The operator string contains NOP elementary operators
* defined by
*
* IOP_TYP : Orbital type of each operator
* IOP_AC  : Creation/annihilation operator
*
*. Output :
* ==========
*
* IOP_REO : New to old order
* sign    : sign of full rank operator
*
* Method : use IPHGAS , and move creation of holes
*                       and annihilation of particles to the left
*
* Jeppe Olsen, October 1997
* version  of : dec 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
*. Specific input, I*OC : Contains number of electrons in each space
*  (not used !!!)
       INTEGER ILOC(NGAS),IROC(NGAS)
       INTEGER IOP_TYP(*), IOP_AC(*)
*. Output, new to old order :
       INTEGER IOP_REO(*)
*
*. The loops are assumed to be of the resolution type, so
*. we will insert a resolution of the identity after NOP/2
*  operators(we assume that NOP is even ) :
*
* <ILOC! IOPL !K><K! IOPR ! IROC >
*=<ILOC! IOPL !K><IROC! (IOPR)+ ! K >
*
*
* and the question is now how to order the operators.
* The operation count can in general be mininized by
*. According to above : Should operators stand left or
*                       right of the resolution K
*
      NTEST = 00
*. First time around : not very general
      NMOVE = 0
      IF(NOP.EQ.2) THEN
*
* One-Electron operator
*
        IF(IOP_AC(1).EQ.2.AND.IOP_AC(2).EQ.1) THEN
* a+1 a 2 : Move around only if both spaces are hole spaces
          NMOVE = 0
          IF(IPHGAS(IOP_TYP(1)).EQ.2) NMOVE = 1
          IF(IPHGAS(IOP_TYP(2)).EQ.2) NMOVE = NMOVE + 1
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) '  IOP_TYP(1), IOP_TYP(2), NMOVE ',
     &                    IOP_TYP(1), IOP_TYP(2), NMOVE
          END IF
        ELSE
          WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
          WRITE(6,*) ' NOP  = ', NOP
          WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
          Call Abend2( ' ALG_ROUTER : Path not implemented ' )
        END IF
*
        IF(NMOVE.EQ.2) THEN
*. Move the operators around in expansion
          IOP_REO(1) = 2
          IOP_REO(2) = 1
          SIGN = -1.0D0
        ELSE
*. No reorganization adviced
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          SIGN = 1.0D0
        END IF
      ELSE IF(NOP.EQ.4) THEN
*
* Two-electron operator
        IF(IOP_AC(1).EQ.2.AND.IOP_AC(2).EQ.2.AND.
     &     IOP_AC(3).EQ.1.AND.IOP_AC(4).EQ.1) THEN
* a+1 a+2 a3 a4
          IF(IPHGAS(IOP_TYP(1)).EQ.2) THEN
            IMOVE1 = 1
          ELSE
            IMOVE1 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(2)).EQ.2) THEN
            IMOVE2 = 1
          ELSE
            IMOVE2 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(3)).EQ.2) THEN
            IMOVE3 = 1
          ELSE
            IMOVE3 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(4)).EQ.2) THEN
            IMOVE4 = 1
          ELSE
            IMOVE4 = 0
          END IF
        ELSE
*. Not yet implemented
          WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
          WRITE(6,*) ' NOP  = ', NOP
          WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
          Call Abend2( ' ALG_ROUTER : Path not implemented ' )
        END IF
*. Number of left and right operators that would like to be moved
        IMOVEL = IMOVE1 + IMOVE2
        IMOVER = IMOVE3 + IMOVE4
        IMOVE = MIN(IMOVEL,IMOVER)
        IF(IMOVE.EQ.2) THEN
*. Well, we are all moving
          IOP_REO(1) = 3
          IOP_REO(2) = 4
          IOP_REO(3) = 1
          IOP_REO(4) = 2
          SIGN = 1.0D0
        ELSE IF(IMOVE.EQ.0) THEN
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          IOP_REO(3) = 3
          IOP_REO(4) = 4
          SIGN = 1.0D0
        ELSE IF(IMOVEL.EQ.1.AND.IMOVER.EQ.1) THEN
* a+ a a+ a
* 1 : Position the two operators to be moved as operators 2, 3
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          IOP_REO(3) = 3
          IOP_REO(4) = 4
          SIGN = 1.0D0
          IF(IMOVE1.EQ.1) THEN
            IOP_REO(1) = 2
            IOP_REO(2) = 1
            SIGN = -SIGN
          END IF
          IF(IMOVE4.EQ.1) THEN
           IOP_REO(3) = 4
           IOP_REO(4) = 3
           SIGN = - SIGN
          END IF
* 2 : and interchange operators 2 and 3
          ISAVE = IOP_REO(2)
          IOP_REO(2) = IOP_REO(3)
          IOP_REO(3) = ISAVE
          SIGN = - SIGN
        ELSE IF (IMOVEL.EQ.2.AND.IMOVER.EQ.1) THEN
* Final operator is a a+ a+ a
          IF(IMOVE4.EQ.1) THEN
            IOP_REO(3) = 4
            IOP_REO(4) = 3
            SIGN = -1.0D0
          ELSE
            IOP_REO(3) = 3
            IOP_REO(4) = 4
            SIGN = 1.0D0
          END IF
          IOP_REO(1) = IOP_REO(3)
          IOP_REO(2) = 1
          IOP_REO(3) = 2
          SIGN = SIGN
        ELSE IF( IMOVEL.EQ.1 .AND. IMOVER .EQ. 2 ) THEN
* Final operator is a+ a a a+
          IF( IMOVE1.EQ.1) THEN
            IOP_REO(1) = 2
            IOP_REO(2) = 1
            SIGN = -1.0D0
          ELSE
            IOP_REO(1) = 1
            IOP_REO(2) = 2
            SIGN = 1.0D0
          END IF
          IOP_REO(4) = IOP_REO(2)
          IOP_REO(2) = 3
          IOP_REO(3) = 4
          SIGN = SIGN
        END IF
*       ^ End of number of pairs to be moved around
      ELSE
*. Not one- or two- electron operator
        WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
        WRITE(6,*) ' NOP  = ', NOP
        WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
        Call Abend2( ' ALG_ROUTER : Path not implemented ' )
      END IF
*

      IF(NOP.EQ.4.AND.NTEST.GE.100) THEN
        WRITE(6,*) ' Information from ALG_ROUTER '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
     &  ' Input : anni- or crea- operators + orbital types'
        CALL IWRTMA(IOP_AC,1,NOP,1,NOP)
        CALL IWRTMA(IOP_TYP,1,NOP,1,NOP)
        WRITE(6,*)
        WRITE(6,*) ' Suggested output order '
        CALL IWRTMA(IOP_REO,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE C_TO_CKK(C,CKK,IWAY,NJA,NJB,NKA,NKB,NTA,NTB,
     &                    IKJA,XKJA,IKJB,XKJB)
*
* IWAY = 1 : C(Ja,Jb) => C(Ka,Kb,ITa,ITb)
* IWAY = 2 : C(Ka,KB,Ita,Itb) => C(Ja,Jb) + C(Ka,Kb,Ita,Itb)
*
* Jeppe Olsen, Summer of 99
*
* Version where matrix notation is eliminated
*
* Notice : K and T is in general batched. This is not visible here
*
*. Input
#include "implicit.inc"
#include "ipoist8.inc"
#if defined(VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif      
#include "parluci.h"
      INTEGER    IKJA(NTA*NKA), IKJB(NTB,NKB)             
      DIMENSION  XKJA(NTA*NKA), XKJB(NTB,NKB)             
*. Input/Output 
C     DIMENSION C(NJA,NJB),CKK(NKA,NKB,*)
      DIMENSION C(NJA*NJB),CKK(*)
csk      WRITE(LUWRT,*)'NJA,NJB, NKA, NKB',NJA,NJB, NKA, NKB
csk      WRITE(LUWRT,*)'NTA,NTB ',NTA,NTB
csk      WRITE(LUWRT,*)'vector IKJA'
csk      CALL IWRTMAMN(IKJA,1,NTA*NKA,1,NTA*NKA,6)
csk      WRITE(LUWRT,*)'vector IKJB', IKJB
csk      CALL IWRTMAMN(IKJB,NTB,NKB,NTB,NKB,6)
*
      DO ITB = 1, NTB
      DO ITA = 1, NTA
       ITAB = (ITB-1)*NTA+ITA-1
*
       KAKBITAB = ITAB*NKA*NKB
csk       WRITE(LUWRT,*)'ITAB',ITAB
       IF(IWAY.EQ.1) THEN
         DO KB = 1, NKB         
          JB = IKJB(ITB,KB)
          IF(JB.NE.0) THEN
csk           WRITE(LUWRT,*)'JB',JB
           SIGNB = XKJB(ITB,KB)
           JAJB0 = (JB-1)*NJA 
           ITAKA = ITA-NTA
           DO KA = 1, NKA           
            ITAKA = ITAKA + NTA
csk            WRITE(LUWRT,*)'ITAKA,NTA',ITAKA,NTA
            JA = IKJA(ITAKA)
csk            WRITE(LUWRT,*)'JA',JA
            JA = IKJA(ITAKA)
            KAKBITAB = KAKBITAB + 1
            IF(JA.NE.0) THEN
                JAJB = JAJB0 + JA
                SIGN = XKJA(ITAKA)*SIGNB
csk                WRITE(LUWRT,*)'KAKBITAB, JAJB',KAKBITAB, JAJB, JAJB0, JA
                CKK(KAKBITAB) = SIGN*C(JAJB)
            ELSE 
               CKK(KAKBITAB) = 0.0D0
            END IF
           END DO
          ELSE IF (JB.EQ.0) THEN
           DO KAKBITAB_EFF =  KAKBITAB+1,  KAKBITAB + NKA  
               CKK(KAKBITAB_EFF) = 0.0D0
           END DO
           KAKBITAB  =  KAKBITAB + NKA
          END IF
*         ^ End if JB .NE. 0
         END DO
       ELSE 
         DO KB = 1, NKB         
          JB = IKJB(ITB,KB)
          IF(JB.NE.0) THEN
           JAJB0 = (JB-1)*NJA
           SIGNB = XKJB(ITB,KB)
           ITAKA = ITA-NTA
           DO KA = 1, NKA           
            ITAKA = ITAKA + NTA
            JA = IKJA(ITAKA)
            KAKBITAB = KAKBITAB + 1
            IF(JA.NE.0) THEN
             JAJB = JAJB0 + JA
             SIGN = XKJA(ITAKA)*SIGNB
             C(JAJB) = C(JAJB) + SIGN*CKK(KAKBITAB)
            END IF
           END DO
          ELSE
            KAKBITAB =  KAKBITAB + NKA
          END IF
*         ^ End if JB = 0
         END DO
       END IF
*      ^ End of IWAY switch
      END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Output from C_TO_CKK, IWAY =  ', IWAY
        WRITE(6,*) ' ======================================='
        WRITE(6,*)
        WRITE(6,*) ' C(Ja,Jb) '
        CALL WRTMAT(C,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' C(Ka,Kb,ITa, Itb) '
        ITAB = 0 
        DO ITB = 1, NTB
        DO ITA = 1, NTA
          ITAB = ITAB + 1
          WRITE(6,*) ' C(Ka,Kb, Ita, Itb) for Ita,Itb =', Ita, Itb
          CALL WRTMAT(CKK(1+ITAB*NKA*NKB),NKA,NKB,NKA,NKB)
        END DO
        END DO
      END IF
*
      RETURN
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE COMPRS2LST_B(I1I,XI1I,N1,I2I,XI2I,N2,NKIN,NKOUT,
     &                        I1O,XI1O,I2O,XI2O)
*
* Two lists of excitations/annihilations/creations are given.
* Compress to common nonvanishing entries
*
* Jeppe Olsen, July 2000 from COMPRS2LST 
*
* Compared to COMPRS2LST : Order of arrys in I1 and I2 interchanged
*                          Output lists differs from input lists
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input lists
      DIMENSION I1I(N1,NKIN),XI1I(N1,NKIN)
      DIMENSION I2I(N2,NKIN),XI2I(N2,NKIN)
*. Output lists
      DIMENSION I1O(N1,NKIN),XI1O(N1,NKIN)
      DIMENSION I2O(N2,NKIN),XI2O(N2,NKIN)
*
      NKOUT = 0
      DO K = 1, NKIN
        I1ACT  = 0
        DO I = 1, N1
          IF(I1I(I,K).NE.0) I1ACT = 1
        END DO
        I2ACT = 0
        DO I = 1, N2
          IF(I2I(I,K).NE.0) I2ACT = 1
        END DO
        IF(I1ACT.EQ.1.AND.I2ACT.EQ.1) THEN
          NKOUT = NKOUT + 1
            DO I = 1, N1
               I1O(I,NKOUT) = I1I(I,K)
              XI1O(I,NKOUT) =XI1I(I,K)
            END DO
            DO I = 1, N2
               I2O(I,NKOUT) = I2I(I,K)
              XI2O(I,NKOUT) =XI2I(I,K)
            END DO
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
CTF          call dmk0one(IBSM,IBTP,JBSM,JBTP,
CTF                       IJBGRP,NIA,
CTF                       NGAS,IBOC,JBOC,
CTF                       SB,CB,
CTF                       ADSXA,STSTSX,
CTF                       NOBPTS,IOBPTS,
CTF                       MAXI,MAXK,
CTF                       SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
CTF                       NSMOB,NSMST,MOCAA,
CTF                       SCLFAC,ONEINTS,INTSPACE,
CTF                       IRI,IRICO,IPM1,IUB,
CTF                       NTEST)
*
*
      subroutine dmk0one(ISCSM,      ! sym of sigma columns
     &                   ISCTP,      ! type of sigma columns
     &                   ICCSM,      ! sym of C columns
     &                   ICCTP,      ! type of C columns
     &                   IGRP,       ! string group of columns
     &                   NROW,       ! number of rows in S and C block
     &                   NGAS,       ! number of active sets
     &                   ISEL,       ! occupation in each active set for S block
     &                   ICEL,       ! occupation in each active set for C block
     &                   SB,         ! output sigma block
     &                   CB,         ! input C block
     &                   ADSXA,      ! Sym(a+a) = (Sym(a+),Sym(a))
     &                   STSTSX,     ! Sym of !st>,sx!st'> => sym of sx
*                                      so <st|sx|st'> is totally symmetric
     &                   NOBPTS,     ! # orbitals per type/symmetry
     &                   IOBPTS,     ! offsets for orbs of type/symmetry
     &                   MAXI,MAXK,  ! max # spectator/inner resol. strings
*                                      treated simultaneously
     &                   SSCR,CSCR,
     &                   I1,XI1S,
     &                   I2,XI2S,
     &                   H,          ! list of transformed integrals
     &                   NSMOB,      ! number of orbital syms
     &                   NSMST,      ! number of string syms
     &                   MOC,        ! use MOC method (inst. of N-1 res.)
     &                   SCLFAC,     ! scaling factor (Jeppe)
     &                   ONEINTS,    ! transformed onel-ints from WORK(KPICI1)
     &                   INTSPACE,   ! length of this vector
     &                   IRI,        ! real/imag. part of sigma vector
     &                   IRICO,      ! real/imag. part of coefficients
     &                   IPM1,       ! sign of Delta Mk (0 here!)
     &                   IUB,        ! UU part (1) or BB (2); 0 if IPM1.ne.0
     &                   IPRNT)     ! Print flag (IPRSIG in input)
*
* One electron excitations on column strings
*. If IH2TRM .ne. 0 then the diagonal and one-electron
*  excitations arising from the two body operator are also
* included
*
*  Real/imaginary logics:
*    S(r) = +h(r) . C(r) - h(i) . C(i)
*    iS(i) = +ih(i) . C(r) - h(r) . iC(i)
*     |      |       |     |       |
*    IRI    SCALEF  IRICO  SCALEF  IRICO
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for one electron integrals
*
*  Jeppe Olsen, Winter of 1991
*
*  changed for relativistic double group calculation
*  Timo Fleig, Summer of 1997
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*  symmetry info for double group CI
#include "symm.inc"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
*  we don't have 2 sets here, as only the relevant one is passed
*  to the routine.
      INTEGER IOBPTS(MXPNGAS,*)
*.Specific Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      character*72 ENVIRO
*
      dimension ONEINTS(INTSPACE)
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),H(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
*.Local arrays
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS)
*
      call qenter('DMK0ONE')
* Type of single excitations that connects the two column strings
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.500)THEN
        WRITE(6,*)
        WRITE(6,*) ' ========================== '
        WRITE(6,*) ' Information from DMK0ONE '
        WRITE(6,*) ' ========================== '
        WRITE(6,*)
CTF     WRITE(6,*) ' DMK0ONE : MOC,IH2TRM ', MOC,IH2TRM
        WRITE(6,*) ' ISEL : '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' ICEL : '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*
      IPM1 = 0
*  (fortunately, +0 = -0)
*
*  The mapping of strings is only done the first time around, i.e.
*  in the first loop for the first Kbatch of strings.
*  After that, IFRST, JFRST (and KFRST, LFRST in dmk1one) are set
*  to zero. In adst_gas, the arc weights if the first index (I) is
*  .ne.0 and the occupation of strings if the second index (here J)
*  is .ne.0.
      IFRST = 1
      JFRST = 1
*
*  1) determine D(Ja,Ib,j+,i) = sum(Jb) <Ib|j_+ i_|Jb> . C(Ja,Jb)
*    a) find single excitations connecting spectator <Ib| = ISEL
*                                         and column |Jb> = ICEL
*
      CALL SXTYP_GAS_REL(NSXTP,JTP,ITP,NGAS,MXPNGAS,ISEL,ICEL)
*
*    b) determine S(j_+ i_) = S*(|Ib>) x S*(|Jb>)
*       i.e. we need adjoint symmetries of the strings
*
      ISCASM = IADJSYM(ISCSM)
      ICCASM = IADJSYM(ICCSM)
      JISM = IDBGMULT(ICCSM,IADJSYM(ISCSM))
      if (JISM.eq.0) goto 1000
*
*  loop over types of single excitations
*
      do 900 JITP = 1,NSXTP,1
         JTYP = JTP(JITP)
         ITYP = ITP(JITP)
         if (NTEST.ge.1000) write(6,*) ' ITYP JTYP ', ITYP,JTYP
         do 800 JSM = 1,NSMOB,1
*
*  loop over symmetries of creator J, resulting in symmetry of creator I (ISM)
*  (sym of creator = adjsym of annihilator)
*
            ISM = IDBGMULT(IADJSYM(JSM),JISM)
            ISM = IADJSYM(ISM)
*  New intermediate strings will be accessed, so
            KFRST = 1
            if (ISM.eq.0) goto 800
            if (NTEST.ge.1000) write(6,*) ' JSM ISM ', JSM,ISM
            NJORB = NOBPTS(JTYP,JSM)
            NIORB = NOBPTS(ITYP,ISM)
            if (NJORB.eq.0.or.NIORB.eq.0) goto 800
*
* Using N-1 resolution method (there is also MOC)
*
* loop over partitionings of N-1 strings
* start with negative values up to zero
* MAXK is given as maximum number of strings per K batch
* inserting batches of N-1 strings according to
*   sum(Kb) <Ib| j_+ |Kb><Kb| i_ |Jb>
*
            KBOT = 1 - MAXK
            KTOP = 0
 700        continue
               KBOT = KBOT + MAXK
               KTOP = KTOP + MAXK
*  offset number for GAS ITYP and symmetry ISM, annihilator
               IBORB = IOBPTS(ITYP,ISM)
               I12 = 1
               K12 = 1
               if (IFRST.eq.1) KFRST = 1
*  determine mapping <Kb| i_ |Jb>, meaning <Jb| i_+ |Kb>, which is real
*
*  Input
*    IBORB : offset number for first spinor of ITYP and ISM
*    NIORB : number of spinors in block ITYP,ISM
*    ICCTP : Type of column vectors Jb
*    ICCSM : Symmetry of column vectors Jb
*    IGRP  : String group of columns
*    KBOT  : First string from K batch
*    KTOP  : Last string from K batch
*    MAXK  : max. no. of K strings per batch
*    IFRST : 1 at start, then 0
*    KFRST : 1 at start, then 0
*    I12   : for local scratch
*    K12   : for local scratch
*    SCLFAC: multiplied with sign in ads1_gas; scaling factor
*
*  Output
*    NKBTC : no. of K strings in present batch
*            in subroutine: NK, no. of K strings checked
*    I1    : I1(KSTR,IORB) depending on ket string K from batch and
*                          annihilator (transferred to creator) IORB
*            value is either 0 (no coupling of Jb and Kb)
*                            +-1 (coupling, strings are connected by IORB)
*    XI1S  : +- a sign might turn up for successful coupling
*    KEND  : if .eq.0 not all K strings have been accessed
*    KKIACT: if some nonvanishing excitations were found, KACT is set to 1,
*            else it is zero.
*
               call ADST_GAS(IBORB,NIORB,ICCTP,ICCSM,IGRP,KBOT,KTOP,
     &                       I1,XI1S,MAXK,NKBTC,KEND,IFRST,KFRST,
     &                       I12,K12,KKIACT,SCLFAC)
*
               IFRST = 0
               KFRST = 0
*. For testing
               IFRST = 1
               KFRST = 1
*
*. Determine mapping <Ib| j_+ |Kb>
*. Set up I2(KSTR,JORB)
               JBORB = IOBPTS(JTYP,JSM)
               J12 = 2
               if (JFRST.eq.1) KFRST = 1
               ONE = 1.0D0
               call ADST_GAS(JBORB,NJORB,ISCTP,ISCSM,IGRP,KBOT,KTOP,
     &                       I2,XI2S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &                       J12,K12,KKJACT,ONE)
               JFRST = 0
               KFRST = 0
* for testing
               JFRST = 1
               KFRST = 1
*
*. loop over partitions of spectator strings
*
               NIPART = NROW/MAXI
               if (NIPART*MAXI.ne.NROW) NIPART = NIPART + 1
               do 701 IPART = 1,NIPART,1
                  IBOT = (IPART-1)*MAXI + 1
                  ITOP = MIN(IBOT+MAXI-1,NROW)
                  NIBTC = ITOP - IBOT + 1
*
*. obtain C(Ib,Kb,i_+) = sum(Jb) <Jb| i_+ |Kb> . C(Ib,Jb)
*  as an intermediate for D, as we still have a partitioning of
*  strings in K batches
*
                  do IIORB = 1,NIORB,1
                     ICGOFF = 1 + (IIORB-1)*NKBTC*NIBTC
                     call matcg(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                          NKBTC,I1(1,IIORB),XI1S(1,IIORB))
                  end do
*
*. obtain one-electron integrals (ISM,ITP,JSM,JTP)
*
                  if (NKBTC.ne.0) then
                     ZERO = 0.D0
                     call setvec(H,ZERO,NIORB*NJORB)
                     call geth1rel3(ONEINTS,H,ISM,ITYP,JSM,JTYP,
     &                       IRI,IRICO,IPM1,IUB,IRETURN,SCALEF,
     &                       INTSPACE,ENVIRO)
*
*  Some important comments:
*
                     if (NTEST.ge.1000) then
                        write(6,*) ' dmk0one H block '
                        if ((IRI.eq.1.and.IRICO.eq.1) .or.
     &                      (IRI.eq.2.and.IRICO.eq.2)) then
                           write(6,*) ' Real part'
                           call wrtmat(H,NIORB,NJORB,NIORB,NJORB)
                        else if ((IRI.eq.1.and.IRICO.eq.2) .or.
     &                           (IRI.eq.2.and.IRICO.eq.1)) then
                           write(6,*) 'Imaginary part'
                           call wrtmat(H,NIORB,NJORB,NIORB,NJORB)
                        else
                           write(6,*) 'Error in r/i structure'
                           Call Abend()
                        end if
                     end if
                  end if
*
*. multiply C(Ib,Kb,i_+) with h(j_,i_) to obtain Sscr(Ib,Kb,j_)
*
                  NIK = NIBTC*NKBTC
                  FACTORC = 0.0D0
                  FACTORAB = SCALEF
                  call matml7(SSCR,CSCR,H,NIK,NJORB,NIK,NIORB,NIORB,
     &                        NJORB,FACTORC,FACTORAB,0)
*. should be modified later, as we have real and imaginary parts of integrals
                  do JJORB = 1,NJORB,1
                     JBSOFF = 1 + (JJORB-1)*NKBTC*NIBTC
                     IF(IBOT.LE.0) THEN
                       WRITE(6,*) ' Start of rows (IBOT) le. 0', IBOT
                       Call Abend2( ' Problem in dmk0one .......' )
                     END IF
*
                     call matcas(SSCR(JBSOFF),SB,NIBTC,NROW,IBOT,NKBTC,
     &                           I2(1,JJORB),XI2S(1,JJORB))
                  end do
 701           continue
*              ^ end of loop over I batches
            if (KEND.eq.0) goto 700
*           ^ end of loop over K batches
 800     continue
*        ^ end of loop over index symmetries
 900  continue
*     ^ end of loop over types of single excitations
1000  continue
*     ^ JISM = 0 condition
2000  continue
      call qexit('DMK0ONE')
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine dmk1one(ISCSMA,ISCSMB,ISCTPA,ISCTPB,
     &                   ICCSMA,ICCSMB,ICCTPA,ICCTPB,
     &                   IAGRP,IBGRP,JAGRP,JBGRP,
     &                   NIROWA,NIROWB,
     &                   NJROWA,NJROWB,
     &                   NGAS,
     &                   ISELA,ISELB,ICELA,ICELB,
     &                   SB,CB,
     &                   ADSXA,STSTSX,
     &                   NOBPTS_U,NOBPTS_B,
     &                   IOBPTS_U,IOBPTS_B,MAXL,MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,
     &                   H,
     &                   NSMOB,NSMST,
     &                   ONEINTS,INTSPACE,SCLFAC,
     &                   IRI,IRICO,IPM1,IUB,IPRNT)
*
*  Info for alpha strings has been added. We need both info for
*  alpha and beta at the same time!!
*  This has to be taken care of in the calling routine,
*  which is: rssbcbrel
*  Be careful: The CA strings are assumed to be N-1 strings already.
*  But the CB strings are assumed to be N+1 strings.
*
* =====
* Input
* =====
*
* ISCSMA,ISCTPA : Symmetry and type of sigma columns for alpha type
* ICCSMA,ICCTPA : Symmetry and type of C     columns for alpha type
* IAGRP : String group of columns, unbarred strings
* IBGRP : String group of columns, barred strings
* NGAS : Number of active sets
* ISELA : Occupation in each active set for sigma block for alpha type
* ICELA : Occupation in each active set for C     block for alpha type
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST : Number of symmetries of orbitals and strings
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* IRI : real or imaginary part of sigma vector?
*       real :1
*       imaginary :2
* IRICO : real or imaginary part of coefficients as input?
*       real :1
*       imaginary :2
* IPM1 : either +1 (delta Mk=+1) or -1 (delta Mk=-1)
* IUB : should be set to 0 here (only relevant in dmk0one)
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for one electron integrals
* (The program will itself decide what type (real or imaginary) of
*   integrals will be needed and these fetched)
*
*  Jeppe Olsen, Winter of 1991
*
*  changed for relativistic double group calculation
*  Timo Fleig, Summer of 1997
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*  symmetry info for double group CI
#include "symm.inc"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS_U(MXPNGAS,*),NOBPTS_B(MXPNGAS,*)
      INTEGER IOBPTS_U(MXPNGAS,*),IOBPTS_B(MXPNGAS,*)
*.Specific Input
      INTEGER ISELA(NGAS),ICELA(NGAS),ISELB(NGAS),ICELB(NGAS)
      character*72 ENVIRO
*
* Input
      dimension ONEINTS(INTSPACE)
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*),H(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
      dimension SSCRH1(MXNSPIPDMGRP*MXNSPIPDMGRP)
*.Local arrays ( assume MPNGAS = 16 ) !!!
C      DIMENSION ITP(16),JTP(16)
* Type of single excitations that connects the two column strings
*
      call qenter('DMK1ONE')
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.500)THEN
        write(6,*)
        write(6,*) ' ========================= '
        write(6,*) ' Information from DMK1ONE '
        write(6,*) ' ========================= '
        write(6,*)
        write(6,*) ' IPM1 = ',IPM1
        write(6,*) ' ISELA : '
        CALL IWRTMA(ISELA,1,NGAS,1,NGAS)
        write(6,*) ' ICELA : '
        CALL IWRTMA(ICELA,1,NGAS,1,NGAS)
        write(6,*) ' ISELB : '
        call iwrtma(ISELB,1,NGAS,1,NGAS)
        write(6,*) ' ICELB : '
        call iwrtma(ICELB,1,NGAS,1,NGAS)
*
        write(6,*) ' NIROWA, NIROWB ', NIROWA, NIROWB
        write(6,*) ' NJROWA, NJROWB ', NJROWA, NJROWB
*
        if (NTEST.ge.2000) then
          WRITE(6,*)  ' Input C block '
          CALL WRTMAT(CB, NJROWA, NJROWB, NJROWA, NJROWB)
          WRITE(6,*)  ' Input S blocks '
          CALL WRTMAT(SB,NIROWA, NIROWB,NIROWA, NIROWB)
        end if
*
        WRITE(6,*) ' ISCSMA,ISCSMB ', ISCSMA,ISCSMB
        WRITE(6,*) ' ICCSMA,ICCSMB ', ICCSMA,ICCSMB
*
      END IF
*
      IFRST = 1
      JFRST = 1
*
*  1) Determine C(Ja,Ib,j) = sum(Jb) <Jb| j_+ |Ib> . C(Ja,Jb)
*     (if IPM1 = -1 then
*               C(Jb,Ia,j) = sum(Ja) <Ja| j+ |Ia> . C(Jb,Ja) )
*    a) find single excitations connecting spectator |Ib> = ISELB
*                                         and column <Jb| = ICELB
*       same for <Ia| i+ |Ja>
*               (<Ib| i_+ |Jb>)
*
* ===========================
* The algorithm for IPM1 = 1
* ===========================
* Sigma(Ia,Ib) = <Ia|i+|Ja><Ib|j_|Jb>C(Ja,Jb) H(i,j)
*
* Part 1 : C(Ja,Kb,j) = C(Ja,Ib,j) <Ib|j_|Jb>C(Ja,Jb)
*          (Proceed as normal 1-electron case)
*
* Part 2 : S(Ja,Ib,i) = sum(j) C(Ja,Ib,j) h(i,j)
*          (as normal 1-electron case )
*
* Part 3 : S(Ia,Ib) = S(Ia,Ib) + <Ia|i+|Ja>S(Ja,Ib,i)
*          Notice this is ROW scattering, to use our normal
*          routines we do the calculation as updates to transposed block
*           Loop over orbital i
*             transpose S(Ja,Ib,i) to S(Ib,Ja,i)
*             Update S(Ib,Ia) = S(Ib,Ia) + <Ia|i+|Ja>S(Ja,Ib,i)
*           End of loop over orbital I
*
*
* =============================
* The algorithm for IPM1 = -1
* =============================
*
* Sigma(Ia,Ib) = <Ia|j|Ja><Ib|i_+|Jb>C(Ja,Jb) h(j,i)
*
* Part 1 : C(         ) = <Ia|j|Ja>C(Ja,Jb)
*    Notice this is row scattering, we therefore proceed as
*    a) C(Jb,Ia,j) = <Ia|j|Ja> C(Jb,Ja)
*    b) C(Jb,Ia,j) => C(Ia,Jb,j)
*
* Part2 : S(Ia,Jb,i) = C(Ia,Jb,j) h(j,i)
*         (as normal 1-electron block)
*
* Part3 : S(Ia,Ib) = S(Ia,Ib) + <Ib|i_+|Jb>S(Ia,Jb,i)
*         (as normal 1-electron block )
*
      IF(IPM1.EQ.1) THEN
*. Transpose sigma block
        CALL trpmt3(SB,NIROWA,NIROWB,CSCR)
        CALL COPVEC(CSCR,SB,NIROWA*NIROWB)
      ELSE
*. Transpose C block
        CALL trpmt3(CB,NJROWA,NJROWB,CSCR)
        CALL COPVEC(CSCR,CB,NJROWA*NJROWB)
      END IF
*
      if (IPM1.eq.1) then
         call adtyp(NSXTPB,JTP,NGAS,ICELB,ISELB)
         call adtyp(NSXTPA,ITP,NGAS,ISELA,ICELA)
      else if (IPM1.eq.-1) then
         call adtyp(NSXTPA,JTP,NGAS,ICELA,ISELA)
         call adtyp(NSXTPB,ITP,NGAS,ISELB,ICELB)
      else
         write(6,*) ' IPM1 out of range '
         write(6,*) ' IPM1 = ',IPM1
         Call Abend2( ' Problem in dmk1one ..... ' )
      end if
*
      if (NTEST.ge.1000) then
         write(6,*) ' ITP,JTP ',ITP,JTP
      end if
*
      if (NSXTPB.eq.0 .or. NSXTPA.eq.0) goto 2000
*
* determine sym(J) for non-vanishing <Jb| j_+ |Ib>
*                                   (<Ja| j+ |Ia>)
* determine sym(I) for non-vanishing <Ia| i+ |Ja>
*                                   (<Ib| i_+ |Jb>)
      if (IPM1.eq.1) then
         JSM = IDBGMULT(ICCSMB,IADJSYM(ISCSMB))
         ISM = IDBGMULT(ISCSMA,IADJSYM(ICCSMA))
         NJORB = NOBPTS_B(JTP,JSM)
         NIORB = NOBPTS_U(ITP,ISM)
      else
         JSM = IDBGMULT(ICCSMA,IADJSYM(ISCSMA))
         ISM = IDBGMULT(ISCSMB,IADJSYM(ICCSMB))
         NJORB = NOBPTS_U(JTP,JSM)
         NIORB = NOBPTS_B(ITP,ISM)
      end if
      if (NTEST.ge.500) write(6,*) ' JSM ISM ', JSM,ISM
      if (NJORB.eq.0 .or. NIORB.eq.0) goto 2000
*
* start with batches of K strings according to sum(Kb) <Jb| j_+ |Kb><Kb|1|Ib>
*                                             (sum(Ka) <Ja| j+ |Ka><Ka|1|Ia>)
      KBOT = 1 - MAXK
      KTOP = 0
 700  continue
         KBOT = KBOT + MAXK
         KTOP = KTOP + MAXK
*
* determine mapping <Jb| j_+ |Kb>
*                  (<Ja| j+ |Ka>)
*
*  For internal affairs (unclear at the moment, to me at least)
*  (Bug fix from middle of march (2000?), set to 2 instead of 1)
         J12 = 2
         K12 = 2
*
         if (JFRST.eq.1) KFRST = 1
         if (IPM1.eq.1) then
            JBORB = IOBPTS_B(JTP,JSM)
            call adst_gas(JBORB,NJORB,ICCTPB,ICCSMB,JBGRP,KBOT,KTOP,
     &                    I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &                    J12,K12,KKJACT,SCLFAC)
         else
            JBORB = IOBPTS_U(JTP,JSM)
            call adst_gas(JBORB,NJORB,ICCTPA,ICCSMA,JAGRP,KBOT,KTOP,
     &                    I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &                    J12,K12,KKJACT,SCLFAC)
         end if
         if (NKBTC.le.0) then
           write(6,*) 'NKBTC after adst: ',NKBTC
           goto 2000
         end if
*
         if (NTEST.ge.1000) write(6,*) 'NKBTC = ',NKBTC
*
* finding out whether e.g. <Kb|Ib> is non-vanishing is done
* "behind the curtain".
*
         JFRST = 0
         KFRST = 0
* loop over batches of La, according to sum(La) <Ia| i+ |La><La|Ja>
*                                      (sum(Lb) <Ib| i_+ |Lb><Lb|Jb>)
         LBOT = 1 - MAXL
         LTOP = 0
 600     continue
            LBOT = LBOT + MAXL
            LTOP = LTOP + MAXL
*
*  determine mapping <Ia| i+ |La>
*                   (<Ib| i_+ |Lb>)
            I12 = 1
            L12 = 1
            if (IFRST.eq.1) LFRST = 1
            if (IPM1.eq.1) then
               IBORB = IOBPTS_U(ITP,ISM)
               call adst_gas(IBORB,NIORB,ISCTPA,ISCSMA,IAGRP,LBOT,LTOP,
     &                       I2,XI2S,MAXL,NLBTC,LEND,IFRST,LFRST,
     &                       I12,L12,LLIACT,SCLFAC)
            else
               IBORB = IOBPTS_B(ITP,ISM)
               call adst_gas(IBORB,NIORB,ISCTPB,ISCSMB,IBGRP,LBOT,LTOP,
     &                       I2,XI2S,MAXL,NLBTC,LEND,IFRST,LFRST,
     &                       I12,L12,LLIACT,SCLFAC)
            end if
            if (NLBTC.le.0) then
              write(6,*) 'NLBTC after adst: ',NLBTC
              goto 2000
            end if
*
            if (NTEST.ge.1000) write(6,*) 'NLBTC = ',NLBTC
*
*  obtain C(Ja,Kb,j) = sum(Jb) <Jb| j_+ |Kb><Kb|Ib> C(Ja,Jb)
*        (C(Jb,Ka,j) = sum(Ja) <Ja| j+ |Ka><Ka|Ia>)
*
            IF(LBOT.LE.0.OR.KBOT.LE.0) THEN
               WRITE(6,*) ' LBOT, KBOT = ',  LBOT, KBOT
               Call Abend2( ' Negative LBOT or KBOT in dmk1one ' )
            END IF
*
            if (IPM1.eq.1) then
               do JJORB = 1,NJORB,1
                  ICGOFF = 1 + (JJORB-1)*NKBTC*NLBTC
                  call matcg(CB,CSCR(ICGOFF),NJROWA,NLBTC,LBOT,NKBTC,
     &                       I1(1,JJORB),XI1S(1,JJORB))
               end do
            else
               do JJORB = 1,NJORB,1
                  ICGOFF = 1 + (JJORB-1)*NKBTC*NLBTC
                  call matcg(CB,SSCR(ICGOFF),NJROWB,NLBTC,LBOT,NKBTC,
     &                       I1(1,JJORB),XI1S(1,JJORB))
                  CALL COPVEC(SSCR(ICGOFF),CSCR(ICGOFF),NKBTC*NLBTC)
               end do
            end if
*
            if (NKBTC.ne.0.and.NLBTC.ne.0) then
               call geth1rel3(ONEINTS,H,ISM,ITP,JSM,JTP,
     &                        IRI,IRICO,IPM1,IUB,IRETURN,SCALEF,
     &                        INTSPACE,ENVIRO)
*
*  Transposition in BU case; structure is always adapted to "alpha-beta"
*  Transposition is required in both cases.
*  Reason is unclear at the moment!
c              if (IPM1.eq.1) then
                 call trpmt3(H,NIORB,NJORB,SSCRH1(1))
                 call copvec(SSCRH1(1),H,NIORB*NJORB)
c              end if
*
               if (NTEST.ge.500) then
               write(6,*) ' dmkp1one H block '
                  if ((IRI.eq.1.and.IRICO.eq.1).or.
     &                (IRI.eq.2.and.IRICO.eq.2)) then
                     write(6,*) ' Real integrals '
                     call wrtmat(H,NIORB,NJORB,NIORB,NJORB)
                  else if ((IRI.eq.1.and.IRICO.eq.2).or.
     &                     (IRI.eq.2.and.IRICO.eq.1)) then
                     write(6,*) ' Imaginary integrals '
                     call wrtmat(H,NIORB,NJORB,NIORB,NJORB)
                  else
                     write(6,*) 'Error in real/imaginary structure'
                     write(6,*) 'Stopping in dmk1one'
                     Call Abend()
                  end if
                  if (NTEST.ge.1000) then
                     write(6,*) 'Checking integral vector'
                     write(6,*) 'in dmk1one ...'
                     do I=1,NIORB*NJORB,1
                        write(6,*) 'H(',I,') = ',H(I)
                     end do
                  end if
               end if
            end if
*
            NLK = NLBTC*NKBTC
            FACTORC = 0.0D0
            FACTORAB = SCALEF
*
*  obtain S(Ja,Kb,i) = sum(j) C(Ja,Kb,j) h(j_,i)
*        (S(Jb,Ka,i) = sum(j) C(Jb,Ka,j) h(j,i_))
            call matml7(SSCR,CSCR,H,NLK,NIORB,NLK,NJORB,NJORB,
     &                  NIORB,FACTORC,FACTORAB,0)
      call memmar(0,0,'CHECK ',0,'MATML7')
*
*  transpose matrix S(La,Kb,i) --> S(Kb,La,i)
*  determine S(Ia,Kb) = S(Ia,Kb) + <Ia| i+ |La> S(Ja,Kb,i)
*           (S(Ib,Ka) = S(Ib,Ka) + <Ib| i_+ |Lb> S(Jb,Ka,i))
*
*  write transpose matrix to CSCR array which not in use at this point
            if (IPM1.eq.1) then
               do IIORB = 1,NIORB,1
                  KOFF = (IIORB-1)*NKBTC*NLBTC + 1
                  call trpmt3(SSCR(KOFF),NLBTC,NKBTC,CSCR(KOFF))
                  call matcas(CSCR(KOFF),SB,NKBTC,NIROWB,KBOT,NLBTC,
     &                        I2(1,IIORB),XI2S(1,IIORB))
               end do
            else
               IF(NKBTC.LE.0.OR.NLBTC.LE.0) THEN
                 WRITE(6,*) ' NKBTC, NLBTC = ', NKBTC, NLBTC
                 Call Abend2( ' Negative batch sizes ' )
               END IF
               do IIORB = 1,NIORB,1
                  KOFF = (IIORB-1)*NKBTC*NLBTC + 1
                  call trpmt3(SSCR(KOFF),NLBTC,NKBTC,CSCR(KOFF))
                  call matcas(CSCR(KOFF),SB,NKBTC,NIROWA,KBOT,NLBTC,
     &                        I2(1,IIORB),XI2S(1,IIORB))
C                      MATCAS(CIN,COUT,NROWI,NROWO,IROWO1,NGCOL,ISCA,SCASGN
               end do
            end if
*
*  final result is on array SB (updated sigma block)
         if (LEND.eq.0) goto 600
      if (KEND.eq.0) goto 700
2000  continue
*. Transpose back to obtain usual sigma block
*
      IF(IPM1.EQ.1) THEN
*. Transpose sigma block
        call trpmt3(SB,NIROWB,NIROWA,CSCR)
        call COPVEC(CSCR,SB,NIROWA*NIROWB)
      ELSE IF( IPM1.EQ.-1) THEN
*. Transpose C block
        call trpmt3(CB,NJROWB,NJROWA,CSCR)
        call COPVEC(CSCR,CB,NJROWB*NJROWA)
      END IF
*
      call qexit('DMK1ONE')
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
* 2 electron contributions with mk=1
*
* The generic form is
*
* sigma(Ia,Ib) = sum(ijkl) <Ia!a+ia!Ja> <Ib!a+jb a kb a lb!Jb>
*                (ij!kl) C(Ja,Jb)
* The current code goes as
*
* Loop over batches of Ib
*
* set up C(Ka,I,Jb) = sum(Ja) <Ia!a+ ia!Ja> C(Ja,Jb>
*
* Set up excitations <Ib!a+jb !Kb>, <Kb! a kb !Mb> <Mb! a lb !Jb>
*
      SUBROUTINE DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
     &                   NIB,IBSM,IBTP,IBGRP,IBOC,
     &                   NJA,JASM,JATP,JAGRP,JAOC,
     &                   NJB,JBSM,JBTP,JBGRP,JBOC,
     &          I_AC,JKL_AC,NGAS,SB,CB,
     &          NALPTS,NBEPTS,MAXK,
     &          SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &          XINT,NSMOB,NSMST,IUSEAB,
     &          CJRES,SIRES,SCLFAC,NTEST,
     &          IUSE_PH,IPHGAS,IBU,JBU,KBU,LBU,
     &          XIJILS)
*
* Note on symmetry handling : I assume that double group symmetry is
* in use, and use info in DBGSYM
*
* =====
* Input
* =====
*
* NIA(NIB,NJA,NJB)         : Number of IA strings
* IASM(IBSM,JASM,JBSM)     : Symmetry of IA strings
* IATP(IBTP,JATP,JBTP)     : Supergroup of strings IA (relative to start of type!)
* IAGRP(IBGRP,JAGRP,JBGRP) : Type (number of electrons) of IA strings
* IAOCC(IBGRP,JAOCC,JBOCC) : Occupation of  each GAS space for IA strings
*
* I_AC = 1, JKL_AC(1-3) = 2,2,1 => a ia a+jb a+kb alb
* I_AC = 2, JKL_AC(1-3) = 2,1,1 => a+ia a+jb a kb alb
*
* NALPTS : Number of alpha spinors per sym and type
* NALPTS : Number of beta  spinors per sym and type
* NGAS : Number of GAS orbital spaces
* CB   : Input C block
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
*       single excitations
* MAXK   : Largest number of inner resolution strings treated at simult.
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT  : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1997
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*. Symmetry info
#include "symm.inc"
*. General input
      INTEGER NALPTS(MXPNGAS,*),NBEPTS(MXPNGAS,*)
*
*.Input
      DIMENSION CB(*)
      INTEGER I_AC(1),JKL_AC(3)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*)
      DIMENSION CJRES(*),SIRES(*)
*
*.Local arrays
      PARAMETER (MXNX = 256)
      DIMENSION ITP(256),JKLTP(256,3)
      DIMENSION IOP_TYP(2),IOP_AC(2),IOP_REO(2)
*
      DIMENSION IJ_TYP(2),IJ_DIM(2),IJ_REO(2),IJ_AC(2),IJ_SYM(2)
      DIMENSION KL_TYP(2),KL_DIM(2),KL_REO(2),KL_AC(2),KL_SYM(2)
*
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
*
      INTEGER MBOC(20),KBOC(20),LBOC(20), KBSPGP(20),MBSPGP(20)

*. Arrays for reorganization
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
C     CHARACTER*6 DBGSYM
C     COMMON/DBGSYM/IDBGMULT(8,8),IPGMULT(4,4),
C    &              NDBGIRR,IMOSF_SP(8,2),ISPSF_MO(8,2),
C    &              IADJSYM(8),INVELM(8)
* Contents of DBGSYM
*    IDBGMULT : Group multiplication table  in double group
*    IPGMULT  : Group multiplication table  in point  group
*    NDBGIRR  : Number of irreps in double group
*    IMOSF_SP : Translation table from MO symmetry to spinor symmetry
*               IMOSF_SP(IOBSM,IAB) Gives the spinor symmetry
*               of an orbital of symmetry IOBSM times an alpha (IAB=1)
*               or and beta spin-function (IAB = 2 )
*    ISPSF_MO : Translation table from spinor symmetry to MO symmetry
*               ISPSF(ISPSM,IAB) gives the symmetry of an MO that
*               combined with an alpha(IAB=1) or beta(IAB=2)
*               spin function give a spinor with symmetry ISPSM
*
*   IADJSYM   : IADJSYM(ISPSYM) is the adjoint symmetry of ISPSYM
*   INVELM    : Inverse elements
*
* The irreps are organized in the usual way with the boson irreps first
*
      CALL QENTER('DMK1TWO')
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ==================== '
        WRITE(6,*) '  SIGMA_MK1 speaking '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
        WRITE(6,*) ' Occupation of groups (IA,IB,JA,JB) '
        WRITE(6,*)
        CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Groups defining each supergroup ( Group number for each gas space )
      CALL GET_SPGP_INF_REL(IATP,IAGRP,IASPGP)
      CALL GET_SPGP_INF_REL(JATP,JAGRP,JASPGP)
      CALL GET_SPGP_INF_REL(IBTP,IBGRP,IBSPGP)
      CALL GET_SPGP_INF_REL(JBTP,JBGRP,JBSPGP)
*. Symmetry of allowed excitations
      IISM = IDBGMULT(IADJSYM(IASM),JASM)
      IF(I_AC(1).EQ.1) THEN
        ISM = IADJSYM(IISM)
      ELSE
        ISM = IISM
      END IF
      JKLSM = IDBGMULT(IADJSYM(IBSM),JBSM)
      IF(ISM.EQ.0.OR.JKLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM ISM ',IASM,JASM,ISM
        write(6,*) ' IBSM JBSM JKLSM ',IBSM,JBSM,JKLSM
      END IF
*.Types of Excitations  that connects the two strings
C     GXFTP(NOP,IAC,ILOCC,IROCC,NGAS,NX,IX,MXNX)
      CALL GXFTP(1,I_AC,IAOC,JAOC,NGAS,NIX,ITP,MXNX)
      CALL GXFTP(3,JKL_AC,IBOC,JBOC,NGAS,NJKLX,JKLTP,MXNX)
      IF(NIX.EQ.0.OR.NJKLX.EQ.0) GOTO 9999
      DO 2001 IX = 1, NIX
*
        ITYP = ITP(IX)
        KAFRST = 1
        NI = NALPTS(ITYP,ISM)
*
*
*. Generate creation- or annihilation- mappings <Ia!a+ ia! Ja>
*. Ordered according to Ia so we set up tables <Ja! a ia !Ia>
        IF(I_AC(1).EQ.1) THEN
          IAC_CONJ = 2
        ELSE
          IAC_CONJ = 1
        END IF
*
        CALL ADAST_GAS_REL(ISM,ITYP,NGAS,JASPGP,JASM,
     &       I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFAC,IAC_CONJ)
*. No compression at the moment
        NKAEFF = NKASTR
*. Loop over batches of KA strings (KA strings = IA strings)
        NKABTC = NKAEFF/MAXK
        IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
*
        DO 1801 IKABTC = 1, NKABTC
           KABOT = (IKABTC-1)*MAXK + 1
           KATOP = MIN(KABOT+MAXK-1,NKAEFF)
           LKABTC = KATOP-KABOT+1
*. Obtain C(ka,i,JB) for Ka in batch
            DO II = 1, NI
              CALL GET_CKAJJB(CB,NI,NJA,CJRES,LKABTC,NJB,
     &             II,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*
            CALL SETVEC(SIRES,ZERO,NIB*LKABTC)
            FACS = 1.0D0
*
            DO 2000 JKL = 1, NJKLX
              JTYP = JKLTP(JKL,1)
              KTYP = JKLTP(JKL,2)
              LTYP = JKLTP(JKL,3)
              IF(NTEST.GE.100) THEN
                WRITE(6,*) ' JTYP, KTYP, LTYP', JTYP,KTYP, LTYP
              END IF
*
              DO 1940 JSM = 1, NSMOB
              DO 1930 KSM = 1, NSMOB
*
                IFIRST = 1
*
                IF(JKL_AC(1).EQ.1) THEN
                  JJSM = IADJSYM(JSM)
                ELSE
                  JJSM = JSM
                END IF
*
                IF(JKL_AC(2).EQ.1) THEN
                  KKSM = IADJSYM(KSM)
                ELSE
                  KKSM = KSM
                END IF
*
                JKSM = IDBGMULT(JJSM,KKSM)
                LLSM = IDBGMULT(INVELM(JKSM),JKLSM)
                IF(JKL_AC(3).EQ.1) THEN
                  LSM = IADJSYM(LLSM)
                ELSE
                  LSM = LLSM
                END IF
*
                IF(LSM.EQ.0) GOTO 1930
                NJ = NBEPTS(JTYP,JSM)
                NK = NBEPTS(KTYP,KSM)
                NL = NBEPTS(LTYP,LSM)
                IF(NJ.EQ.0.OR.NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
* we will set the following mappings up
*   <KB! a j (conj) !IB>
*   <MB! a k (conj) !KB>
*   <JB! a l (conj) !MB>, drivenm by the strings to the right
*. Find symmetry, group and type of  KB and MB
                KBSM = IDBGMULT(JJSM,IBSM)
                MBSM = IDBGMULT(KKSM,KBSM)
*
                CALL ICOPVE(IBOC,KBOC,NGAS)
                IF(JKL_AC(1).EQ.1) THEN
                  KBOC(JTYP) = KBOC(JTYP) + 1
                  KBTP = IBTP + 2
                  JCONJ = 2
                ELSE
                  KBOC(JTYP) = KBOC(JTYP) - 1
                  KBTP = IBTP - 2
                  JCONJ = 1
                END IF
*
                CALL ICOPVE(KBOC,MBOC,NGAS)
                IF(JKL_AC(2).EQ.1) THEN
                  MBOC(JTYP) = KBOC(JTYP) + 1
                  MBTP = KBTP + 2
                  KCONJ = 2
                ELSE
                  MBOC(JTYP) = KBOC(JTYP) - 1
                  MBTP = KBTP - 2
                  KCONJ = 1
                END IF
*
                IF(JKL_AC(3).EQ.1) THEN
                  LCONJ = 2
                ELSE
                  LCONJ = 1
                END IF
*   <KB! a j (conj) !IB>
                CALL ADAST_GAS_REL(JSM,JTYP,NGAS,KBSPGP,KBSM,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,JKL_AC(1))
                IF(NKBSTR.EQ.0) GOTO 1930
*   <MB! a k (conj) !KB>
                CALL ADAST_GAS_REL(KSM,KTYP,NGAS,MBSPGP,MBSM,
     &               I3,XI3S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,JKL_AC(2))
                IF(NKBSTR.EQ.0) GOTO 1930
*   <JB! a l (conj) !MB>
                CALL ADAST_GAS_REL(LSM,LTYP,NGAS,JBSPGP,JBSM,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,JKL_AC(3))
* Fetch Integrals as (i j |  k l )
*
                IXCHNG = 0
                ICOUL = 1
                CALL GETINT_DBG(XINT,ISM,ITYP,IBU,JSM,JTYP,JBU,
     &                           KSM,KTYP,KBU,LSM,LTYP,LBU,
     &                           IRI,IXCHNG,0,0,ICOUL)
*
* S(Ka,Ib) = sum(i,j,k,l,Jb)<Ib!a jb a kb a lb!Jb>C(Ka,i,Jb)*(ij!kl)
*
                CALL SKICKJ_MK1(SB,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,NKBSTR,NKBSTR,XINT,NI,NJ,NK,NL,
     &               NKBSTR,NKBSTR,I2,XI2S,I3,XI3S,I4,XI4S,IKORD,
     &               FACS,MXTSOB,XIJILS)
C               SKICKJ_MK1(SKII,CKJJ,NKA,NIB,NJB,NKB,NMB,XIJKL,
C    &                  NI,NJ,NK,NL,MAXI,MAXK,MAXM
C    &                  KBIB,XKBIB,MBKB,XMBKB,JBMB,XJBMB,IKORD,
C    &                  FACS,MXTSOB,XIJILS)

 1930         CONTINUE
*             ^ End of loop over KSM
 1940       CONTINUE
*          .^ End of loop over ISM
 2000       CONTINUE
*           ^ End of loop over KLTYP
 1801     CONTINUE
*.        ^End of loop over partitioning of alpha strings
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      CALL QEXIT('DMK1TWO')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE DXFSTR(IKIN1,IKIN2,IKST,IKSG,NIK,
     &                  KJIN1,KJIN2,KJST,KJSG,NKJ,
     &                  NI,NJ,NK,NL,
     &                  IOFF,JOFF,KOFF,LOFF,
     &                  ISTR,MXSXST,XIJKL,
     *                  IEXSTR,FACSTR,NEX,NTESTG)
*. Two list of single excitations occurs
*
*  IK : String I => String K (IKIN1,IKIN2,IKST,IKSG)
*  KJ : String K => String J (KJIN1,KJIN2,KJST,KJSG)
*
* Obtain effect of two-electron excitation part of two-electron operator
*
*
* g = 1/2sum(i,j,k,l) a+ia+ka l a j (ij!kl)
*   = sum(i>k,j>l)    a+ia+ka l a j ((ij!kl) -(il!jk))
*   = sum(i>k,j>l,j.ne.k)    a+ia j a+k a l ((ij!kl) -(il!jk))
*   + sum(i>k,k>l)    a+ia k a+k a l ((ik!kl) -(il!kk))
*   = sum(i>k,j>l,j.ne.k)    a+ia j a+k a l ((ij!kl) -(il!jk))
*   - sum(i>k,k>l)    a+ia+k a k a l ((ik!kl) -(il!kk))
*   = sum(i>k,j>l,j.ne.k)    a+ia j a+k a l ((ij!kl) -(il!jk))
*   - sum(i>k,k>j)    a+ia j a+k a k ((ik!kj) -(ij!kk))
*
* COLD :
*   = sum(i>k,j.ne.k, j.gt.l)
*     a+ia j a+k a l ((ij!kl) -(il!jk))
*   + sum(i>k, j>l)
*     a+ia j a+j a l ((ij!jl) -(il!jj))
*
* as
*   G !ISTR> = sum(i>k,j.ne.k,if(k.ne.l ) j.gt.l)
*     a+ia ji!KSTR><KSTR! a+k a l!ISTR> ((ij!kl) -(il!jk))
*
* Version only employing single annihilation tables
*
* =====
* Input
* =====
*
* IKIN1 : Index one of <K STR! SX ! ISTR > i.e k
* IKIN2 : Index two of <K STR! SX ! ISTR > i.e l
* IKST  : Corresponding K string
* IKSG :  Phase shift of excitation
* NIK  : Number of excitations per string I=> K
*
* KJIN1 : Index one of <J STR! SX ! KSTR > i.e i
* KJIN2 : Index two of <J STR! SX ! KSTR > i.e j
* KJST  : Corresponding J string
* KJSG :  Phase shift of excitation
* NKJ  : Number of excitations per string K => J
*
* NI,NJ,NK,NL : Number of orbitals in each group
* IOFF,JOFF,KOFF,LOFF : Offset for each group of orbitals.
*
*
* XIJKL : Integrals stored as (j,l,i,k)
* ======
* Output
* ======
*
*     NEX    : Number of excitations
*
*     IEXSTR : IEXSTR(K) is string number of excited string K
*              relative to start of TS block
*     FACSTR : FACSTR(K) is the factor that string INDSTR(K) is multiplied
*              with
*
*     Jeppe Olsen, March 1994, LUCIA version of DXFSTR
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      INTEGER IKIN1(MXSXST,*),IKIN2(MXSXST,*)
      INTEGER IKST(MXSXST,*),NIK(*)
      REAL*8  IKSG(MXSXST,*)
*
      INTEGER KJIN1(MXSXST,*),KJIN2(MXSXST,*)
      INTEGER KJST(MXSXST,*),NKJ(*)
      REAL*8  KJSG(MXSXST,*)
*
      DIMENSION XIJKL(NJ,NL,NI,NK)
*
*. Output
      INTEGER IEXSTR(*)
      DIMENSION FACSTR(*)
*
      CALL QENTER('DXFST')
*
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ****************************************'
        WRITE(6,*) ' DXFSTR at your command : ISTR = ',ISTR
        WRITE(6,*) ' ****************************************'
      END IF
*
      NEX = 0
C?    WRITE(6,*) ' Number of excitations from string ', NIK(ISTR)
      DO 1000 KLEX = 1, NIK(ISTR)
        K = IKIN1(KLEX,ISTR)
        L = IKIN2(KLEX,ISTR)
        KSTR = IKST(KLEX,ISTR)
        SIGN1 = IKSG(KLEX,ISTR)
C?      WRITE(6,*) ' KLEX K L ', KLEX,K,L
C?      WRITE(6,*) ' KSTR SIGN1 ',KSTR,SIGN1
        DO 900 IJEX = 1, NKJ(KSTR)
          I = KJIN1(IJEX,KSTR)
          J = KJIN2(IJEX,KSTR)
C?        WRITE(6,*) ' I J ', I,J
          ICASE = 0
          IF(I.GT.K .AND. J.GT.L .AND.
     &       I.NE.L.AND.I.NE.J.AND.K.NE.L.AND.K.NE.J) THEN
            ICASE = 1
          END IF
          IF(ICASE.EQ.0) GOTO 900
          JSTR = KJST(IJEX,KSTR)
C?        write(6,*) ' KSTR = ', KSTR
C?        write(6,'(A,5I3)')
C?   &   ' Allowed excitation i,j,k,l  JSTR ',I,J,K,L,JSTR
          IF(KJSG(IJEX,KSTR).EQ.1.0D0) THEN
            SIGN =   SIGN1
          ELSE
            SIGN = - SIGN1
          END IF
C?        write(6,*) ' jstr,kstr sign ', jstr,kstr,sign
*
          NEX = NEX + 1
          IEXSTR(NEX) = JSTR
          IF(ICASE.EQ.1) THEN
C    &      SIGN*XIJKL(J-JOFF+1,L-LOFF+1,I-IOFF+1,K-KOFF+1)
            FACSTR(NEX) =
     &      SIGN*(GTIJKL_REL(I,J,K,L)-GTIJKL_REL(I,L,K,J))
C    &      SIGN*XIJKL(J-JOFF+1,L-LOFF+1,I-IOFF+1,K-KOFF+1)
          END IF
C         IF(FORM1.NE.FACSTR(NEX) ) THEN
C           write(6,*) ' problemo in DXFSTR'
C           write(6,*) ' i j k l Form1 FACSTR '
C           write(6,*) i,j,k,l,FORM1,FACSTR(NEX)
C           write(6,*) ' ijkl and ilkj '
C           write(6,*) GTIJKL_REL(I,J,K,L),GTIJKL_REL(I,L,K,J)
C         END IF
  900   CONTINUE
 1000 CONTINUE
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 'Number of excitations obtained ', NEX
        WRITE(6,*)
        WRITE(6,*) ' iex      exc.string         factor '
        WRITE(6,*) ' =================================='
        DO IEX = 1, NEX
          WRITE(6,'(I4,4X,I8,6X,F13.5)')
     &    IEX,IEXSTR(IEX),FACSTR(IEX)
        END DO
      END IF
*
      CALL QEXIT ('DXFST')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GXFTP(NOP,IAC,ILOCC,IROCC,NGAS,NX,IX,MXNX)
*
* Occupation of left and right string given by ILOCC,IROCC
* These are connected with a string containing NOP operators.
*
* <ILOCC! IOP ! IROCC>
*
* Find possible gasspaces of the each operator
*
* Jeppe Olsen, Winter of 1998
*
*      IMPLICIT REAL*8(A-H,O-Z)
*
* ======
*. Input
* ======
*
* NOP : Number of operators
* IACI(I) : =1 => operator I is annihilation
*           =2 => operator I is creation
* ILOCC   : Occupation of left  string
* IROCC   : Occupation of right string
*
* ======
* Output
* ======
*
* NX : Number of connecting excitations
* IX : Connecting excitations
* MXNX : Max number of connecting excitations
*
*. Input
CTI  This array has to be passed from above or put on a common block
CTI      DIMENSION IOP(NOP)
      INTEGER ILOCC(NGAS),IROCC(NGAS)
#include "ipoist8.inc"
*. output
      DIMENSION IX(MXNX,NOP)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE K_TO_J_TOT(IKJ,XKJ,KSM,IK_B,IK_E,IT,ITSM,NTOP,
     &                      IT_B,IT_E,IM,XM,IBM,NK,LTOP,IZERO_MAP,IUBL)
*
* Obtain Total map !J> = T-Oper !K>
* T-oper : string of elementary operators 
* !K> String of sym KSM
*
* !K> are of sym KSM and are restricted to IK_B to IK_E ( within sym)
* T_operators are of sym ITSM and are restricted to IT_B to IT_E
*
* Map for each elementary operator is provided by 
*  IM with offset mat IBM
*
* Jeppe Olsen, Summer of 99
* Double group symmetry incorporated for relativistic calculations
* IUBL flag for unbarred/barred operator strings introduced.
*  Timo Fleig, July 2001, January 2002
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "cgas.inc"
#include "csm.inc"
*. Input
      INTEGER IT(LTOP,*),IM(*),
     &        IBM(MXNDGIRR,MXNDGIRR,LTOP),NK(MXPOBS,LTOP)
      DIMENSION XM(*)
*. Output
      INTEGER IKJ(NTOP,(IK_E-IK_B+1))
      DIMENSION XKJ(NTOP,(IK_E-IK_B+1))
*
      NTEST = 0000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Information from K_TO_J.... '
        WRITE(6,*) ' =========================== '
        WRITE(6,*) ' LTOP, NTOP = ', LTOP,NTOP
        WRITE(6,*) ' IK_B,IK_E,IT_B,IT_E = ', IK_B,IK_E,IT_B,IT_E
        WRITE(6,*) ' IUBL, KSM = ',IUBL,KSM
        write(6,*) 'Offset matrix IBM:'
        do L=1,LTOP,1
          write(6,*) 'Ltop value ',L
          call iwrtma(IBM(1,1,L),NIRR_DG,NIRR_DG,MXNDGIRR,MXNDGIRR)
        end do
      END IF
C?    WRITE(6,*) ' Dim of IKJ, XKJ ',NTOP,IK_E-IK_B+1
*
      THRESH = 1.D-30
      IZERO_MAP = 1
      IF(LTOP.GT.0) THEN
        DO KSTR = IK_B, IK_E
C?        WRITE(6,*) ' KSTR = ', KSTR
          DO ITOP = IT_B, IT_E
C?          WRITE(6,*) ' ITOP = ', ITOP
            KNSM = KSM
            SIGN = 1.0D0
            KNSTR = KSTR
C?          WRITE(6,'(A,10I4)') ' T-op : ', (IT(JOP,ITOP),JOP=1,LTOP)
            DO IOP = 1, LTOP
C?            WRITE(6,*) ' Info for IOP = ', IOP
              JOB_ABS = IT(IOP,ITOP)
              if (IUBL.eq.1) then
                JOB_SM = ISMFTO(JOB_ABS)
                JOB_TP = ITPFTO(JOB_ABS)
              else if (IUBL.eq.2) then
                JOB_SM = ISMFTO2(JOB_ABS)
                JOB_TP = ITPFTO2(JOB_ABS)
              end if
C?            WRITE(6,*) ' JOB_ABS, IBTS, JOB_SM ',
C?   &        JOB_ABS, JOB_TP, JOB_SM
              if (IUBL.eq.1) then
                JOB_REL = JOB_ABS - IOBPTS(JOB_TP,JOB_SM) + 1
              else if (IUBL.eq.2) then
                JOB_REL = JOB_ABS - IOBPTS2(JOB_TP,JOB_SM) + 1
              end if
              LK = NK(KNSM,IOP)
C?            WRITE(6,*) ' IBM,JOB_REL,LK,KNSTR',
C?   &        IBM(JOB_SM,KNSM,IOP),JOB_REL,LK,KNSTR
              IADR = IBM(JOB_SM,KNSM,IOP) -1 + (JOB_REL-1)*LK+KNSTR
C?            WRITE(6,*) ' IADR = ', IADR
              KNSTR = IM(IADR)
C?            WRITE(6,*) ' IADR, KNSTR =  ', IADR, KNSTR
              SIGN = SIGN*XM(IADR) 
              IF(KNSTR.EQ.0) GOTO 1001
              KNSM = IDBGMULT(KNSM,JOB_SM)
            END DO
*
 1001       CONTINUE
            IF(KNSTR.NE.0)  IZERO_MAP = 0
            IKJ(ITOP-IT_B+1,KSTR-IK_B+1) = KNSTR
            XKJ(ITOP-IT_B+1,KSTR-IK_B+1) = SIGN
C?          WRITE(6,*) ' K_TO_J Row and col : ', 
C?   &      ITOP-IT_B+1,KSTR-IK_B+1
          END DO
        END DO
      ELSE
*. No K-operators => Identity map
        IZERO_MAP = 0
        DO KSTR = IK_B,  IK_E
          IKJ(1,KSTR-IK_B+1) = KSTR
          XKJ(1,KSTR-IK_B+1) = 1.0D0
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Operators '
        print*,'IT_B is ',IT_B
        CALL IWRTMA(IT(1,IT_B),LTOP,IT_E-IT_B+1,LTOP,IT_E-IT_B+1)
        WRITE(6,*) ' Output from K_TO_J_TOT ' 
        WRITE(6,*) 'IKJ and XKJ arrays '
        LT = IT_E-IT_B+1
        LK = IK_E-IK_B+1
        CALL IWRTMA(IKJ,LT,LK,LT,LK)
        CALL WRTMAT(XKJ,LT,LK,LT,LK)
      END IF
*
      RETURN
      END
* 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MSAXPYN6(AX,A,X,NDIM,NVEC,INDEX)
*
* AX(I) = AX(I) + SUM(L=1,NVEC) A(L)*X(I,INDEX(L))
*
* Initial RISC/6000 version
*
* Jeppe Olsen
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION AX(*),X(NDIM,*)
      DIMENSION A(*) ,INDEX(*)
*
* level of unrolling of column loop
      NROL = 5
      NRES = MOD(NVEC,NROL)
*
*. part that cannot be column unrolled
*
      IF(NRES.EQ.0.OR.NVEC.EQ.0) THEN
*. Nothing to do !
      ELSE IF (NRES.EQ.1) THEN
        I1 = INDEX(1)
        A1 = A(1)
        DO 20 I = 1, NDIM
          AX(I) = AX(I) + A1*X(I,I1)
   20   CONTINUE
      ELSE IF (NRES.EQ.2) THEN
        I1 = INDEX(1)
        A1 = A(1)
        I2 = INDEX(2)
        A2 = A(2)
        DO 30 I = 1, NDIM
         AX(I) = AX(I) + A1*X(I,I1)+A2*X(I,I2)
   30   CONTINUE
      ELSE IF (NRES.EQ.3) THEN
        I1 = INDEX(1)
        A1 = A(1)
        I2 = INDEX(2)
        A2 = A(2)
        I3 = INDEX(3)
        A3 = A(3)
        DO 40 I = 1, NDIM
         AX(I) =
     &   AX(I) + A1*X(I,I1)+A2*X(I,I2)+A3*X(I,I3)
   40   CONTINUE
      ELSE IF (NRES.EQ.4) THEN
        I1 = INDEX(1)
        A1 = A(1)
        I2 = INDEX(2)
        A2 = A(2)
        I3 = INDEX(3)
        A3 = A(3)
        I4 = INDEX(4)
        A4 = A(4)
        DO 50 I = 1, NDIM
         AX(I) =
     &   AX(I) + A1*X(I,I1)+A2*X(I,I2)+A3*X(I,I3)+A4*X(I,I4)
   50   CONTINUE
      ELSE IF (NRES.GE.5) THEN
        WRITE(6,*) ' NRES to large in MSAXPY so STOP '
        Call Abend2( ' MSAXPY ' )
      END IF
*. unrolled column loops
      DO 285 L = 1+NRES,NVEC,NROL
        AL= A(L)
        IL= INDEX(L)
        AL1= A(L+1)
        IL1 = INDEX(L+1)
        AL2= A(L+2)
        IL2 = INDEX(L+2)
        AL3= A(L+3)
        IL3 = INDEX(L+3)
        AL4= A(L+4)
        IL4 = INDEX(L+4)
        DO 205 I = 1, NDIM
          AX(I) = AX(I) + AL*X(I,IL)+AL1*X(I,IL1)
     &          + AL2*X(I,IL2)+AL3*X(I,IL3)+AL4*X(I,IL4)
  205   CONTINUE
  285 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MV7_REL(C,HC,LUC,LUHC,INONREL)
*
* Outer routine for sigma vector generation
* GAS version !!!!
*
* Written in terms of RASG3/SBLOCK, May 1997
* Works only for ICISTR .ne. 1, i.e. pure memory version
* has been eliminated.
*
* Generalized to allow for relativistic sigma, T.F. + J.O., march 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma
#include "cands.inc"
*
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
*
      common/glbref/ISPC,ISM,ITERGLB
*
CTF  Variable INONREL added to list. It is passed from CMICDV
*
      CALL QENTER('MV7REL')
      IDUM = 0
      KDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'MV7REL')
*
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*) ' MV7_REL does not work for ICISTR = 1'
        WRITE(6,*) ' SWITCH to ICISTR = 2,3 or program'
        Call Abend2( ' MV7_REL does not work for ICISTR = 1' )
      END IF
*. Arrays for partitioning of sigma
* (*_MS2 arrays are in CSTATE )
      CALL Z_BLKFO_REL(IDC,NMS2VAL,
     &                 ISPC,ISM,KLSLBT,KLSLEBT,KLSI1BT,KLSIBT,KLSBLTP,
     &                 NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,
     &                 IBAT_MS2,IDUMMY,.FALSE.,0,0,IDUMMY)
*
      CALL RASSG3_REL(C,HC,NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                WORK(KLSI1BT),WORK(KLSIBT),LUC,LUHC)
C          RASSG3(CB,SB,LBATS,LEBATS,I1BATS,IBATS,
C    &     LUC,LUHC)
*. Eliminate local memory
      call memmar(KDUM ,IDUM,'FLUSM ',2,'MV7REL')
*
      CALL QEXIT('MV7REL')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MXRESC_REL(IAB,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                      NSMST,NSTFSMSPGP,NSTFSMSPGP2,
     &                      MXPNSMST,
     &                      NSMOB,MXPTOB,NTPOB,NTSOB,NTSOB2,
     &                      NTESTG,MXPKA,
     &                      NEL1234,
     &                      MXCJ,MXCIJA,
     &                      MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
*
* Find largest dimension of matrix C(Ka,Ib,J)
* Find largest dimension of matrix C(ij,Ka,Ib)
* Find largest dimension of matrix C(ij,Ia,Kb)
* Find largest dimension of matrix C(ij,Ka,Kb)
*
* Largest block of single excitations MXSXBL

*. Input
* IAB :allowed combination of alpha and beta supergroups
* IOCPTA : Number of first active alpha supergroup
* IOCPTB : Number of first active beta  supergroup
* NOCTPA : Number of active alpha supergroups
* NOCTPB : Number of active alpha supergroups

      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION IAB(NOCTPA,NOCTPB)
      DIMENSION NSTFSMSPGP(MXPNSMST,*),NSTFSMSPGP2(MXPNSMST,*)
      DIMENSION NTSOB(MXPTOB,NSMOB),NTSOB2(MXPTOB,NSMOB)
      DIMENSION NEL1234(MXPTOB,*)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.100) WRITE(6,*) ' MXRESC : MXPKA ', MXPKA
*
* matrix C(j,Ka,Ib)
*
*. Note : Only done for alpha-strings, problems when transposing
*         constructing C(J,Ia,Kb)
      MXCJ = 0
      MXADKBLK = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(j,Ka,Ib)'
        write(6,*)
      end if
      DO 100 IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO 200 IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXB = 0
            DO 210 ISM = 1, NSMST
              MXB = MAX(MXB,NSTFSMSPGP2(ISM,IBTPABS))
  210       CONTINUE
            IF(NTEST.GE.100) WRITE(6,*) ' MXB = ', MXB
            DO 300 IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP 
*  from IATP
              CALL NEWTYP_REL(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO 310 KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
  310           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                MXKAO = MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKA= MXPKA
                MXSOB = 0
                DO 320 ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
  320           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKAO)
                LCJBLK = MXSOB*MXKA*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
  300       CONTINUE
          END IF
  200   CONTINUE
  100 CONTINUE
*
* matrix C(j,Ia,Kb)
*
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(j,Ia,Kb)'
        write(6,*)
      end if
*
      DO IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXA = 0
            DO ISM = 1, NSMST
              MXA = MAX(MXA,NSTFSMSPGP(ISM,IATPABS))
            END DO
            IF(NTEST.GE.100) WRITE(6,*) ' MXA = ', MXA
            DO IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP 
*  from IATP
              CALL NEWTYP_REL(IBTPABS,1,IOBTP,1,KBTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KBTP ',IOBTP,KBTP
              IF(KBTP.GT.0) THEN
                MXKB = 0
                DO KSM = 1, NSMST
                  MXKB = MAX(MXKB,NSTFSMSPGP2(KSM,KBTP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                MXKBO = MXKB
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKB= MXPKA
                MXSOB = 0
                DO ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB2(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKBO)
C               LCJBLK = MXSOB*MXKB*MXB
                LCJBLK = MXSOB*MXKB*MXA
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXADKBLK,MXCJ ', MXADKBLK,MXCJ
      END IF
*
* matrix C(ij,Ka,Ib)
* both Ka and Ib blocked
*
      MXCIJA = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ka,Ib)'
        write(6,*)
      end if
*
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA -1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIB = 0
            DO  ISM = 1, NSMST
              MXIB = MAX(MXIB,NSTFSMSPGP2(ISM,IBTPABS))
            END DO
            IF(MXIB.GT.MXPKA) MXIB = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP 
*  from IATP
              CALL NEWTYP_REL(IATPABS,1,IOBTP,1,K1ATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1ATP ',IOBTP,K1ATP
              IF(K1ATP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP 
*  from K1ATP
                  CALL NEWTYP_REL(K1ATP,1,JOBTP,1,KATP)
                  IF(KATP.GT.0) THEN
                    MXKA = 0
                    DO KSM = 1, NSMST
                      MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                    IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &              MXKA= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKA*MXIB
                    MXCIJA = MAX(MXCIJA,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJA ', MXCIJA
      END IF
*
*
* matrix C(ij,Ia,Kb)
* both Ka and Ib blocked
*
      IF(NTEST.GE.100) WRITE(6,*) ' MXCIJB under development '
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ia,Kb)'
        write(6,*)
      end if
*
      MXCIJB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIA = 0
            DO  ISM = 1, NSMST
              MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS ))
            END DO
            IF(MXIA.GT.MXPKA) MXIA = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP 
*  from IBTP
              CALL NEWTYP_REL(IBTPABS,1,IOBTP,1,K1BTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1BTP ',IOBTP,K1BTP
              IF(K1BTP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB2(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP 
*  from K1ATP
                  CALL NEWTYP_REL(K1BTP,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP2(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB2(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXIA
                    MXCIJB = MAX(MXCIJB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJB ', MXCIJB
      END IF
*
*
* matrix C(ij,Ka,Kb)
* both Ka and Kb blocked
*
      MXCIJAB = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ka,Kb)'
        write(6,*)
      end if
*
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            DO  IOBTP = 1, NTPOB
*. type of Ka string obtained by removing one elec of type IOPBTP 
*  from IATP
              CALL NEWTYP_REL(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) MXKA= MXPKA
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP 
*  from IBTP
                  CALL NEWTYP_REL(IBTPABS,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP2(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB2(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXKA
                    MXCIJAB = MAX(MXCIJAB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
*
* Largest block of single excitations :
* Strings of given type and sym, orbitals of given type and sym
*
* Largest block of creations : a+i !kstring> where K string is
* obtained as single annihilations
      MXSXBL = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Largest block of creations a+i |kstring>'
        write(6,*) 'alpha strings'
        write(6,*)
      end if
*
*. For alpha strings :
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        MXIA = 0
        DO  ISM = 1, NSMST
          MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP_REL(IATPABS,1,JOBTP,1,KATP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KATP ',JOBTP,KATP
          IF(KATP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IATPABS))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP_REL(KATP,2,IOBTP,1,JATP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
*. For beta  strings :
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Largest block of creations a+i |kstring>'
        write(6,*) 'beta strings'
        write(6,*)
      end if
*
      DO  IBTP = 1, NOCTPB
        IBTPABS = IBTP + IOCTPB - 1
        MXIB = 0
        DO  ISM = 1, NSMST
          MXIB = MAX(MXIB,NSTFSMSPGP2(ISM,IBTPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP_REL(IBTPABS,1,JOBTP,1,KBTP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KBTP ',JOBTP,KBTP
          IF(KBTP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB2(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IBTP))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP_REL(KBTP,2,IOBTP,1,JBTP)
              IF(JBTP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB2(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIB)
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC: MXSXBL : ', MXSXBL
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE RASSG3_REL(CB,SB,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &                      LLUC,LLUHC)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen   Winter of 1991
*               May 1997 : Connected to SBLOCK
*
*               March 1998: Relativistic version
*                           real/imaginary loops,
*                           MK2 coupling
*               (Jeppe Olsen, with a little help from my friend)
*
*
*=========================================================
* SUMMARY OF THE VECTOR HANDLING (I/O) IN SIGMA ROUTINES:
*=========================================================
*
* rassg3:	Set current sigma batch to zero (no disc operation)
*                 CALL SETVEC(SB,ZERO,LS)
* rassg3:	Copy blocks with given MK2_C to LUSC37
*                 CALL COPNBLKD(LLUC,LUSC37,CB,NBL_C,0,LBLK)
* rassg3:       Call to sblock with given IRIS,MK2_S,IRIC,MK2_C
* rassg3:	Transfer S block to permanent storage
*                 CALL TODSC(SB(IOFF),LEN_T,-1,LLUHC)
*=========================================================
*
*=========================================================
* OVERVIEW OVER TRANSPOSITIONS OF S AND C BLOCKS
*=========================================================
*
* The blocks are assumed transposed at the beginning, meaning
* a global transposition is needed in sblocks/dblocks.
* General rule: For processing a beta block, a transposition is needed.
*
* sblocks:	Offdiagonal blocks, explicit transposing
*               (Always done with PSSIGN and PLSIGN = 0 => IDC = 1
*               CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
*               CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
*                - - - - - - - - - - - - - - - - - - - - - - -
* rssbcbrel:	  <beta || beta>      TRP(S)  TRP(C)
*                                     do matrix operations
*                                     TRP(S)  TRP(C)
*                - - - - - - - - - - - - - - - - - - - - - - -
*                 <alpha || alpha>    do matrix operations
*                - - - - - - - - - - - - - - - - - - - - - - -
*                 <alpha || beta>     TRP(S)  TRP(C)
*                                     TRP(S)
*                                     do matrix operations
*                                     TRP(S)
*                                     TRP(S)  TRP(C)
*                - - - - - - - - - - - - - - - - - - - - - - -
*                 <beta || alpha>             TRP(C)
*                                     do matrix operations
*                                             TRP(C)
*                                     TRP(S)  TRP(C)
*                - - - - - - - - - - - - - - - - - - - - - - -
* sblocks:	Transpose final sigma block back.
*
* =====
* Input
* =====
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Batches of sigma
      INTEGER LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
*.Scratch
      DIMENSION SB(*),CB(*)
*
#include "mxpdim.inc"
#include "cstate.inc"
#include "clunit.inc"
#include "gasstr.inc"
#include "cprnt.inc"
*
*.Definition of c and sigma spaces
#include "cands.inc"
*
* IA/BS/CTP : Type (number of elecs in A/B,S/C)
*. CSTATE
C     COMMON/CSTATE/IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS,
C    &              NROOT,IROOT(100),IDC,PSSIGN,PLSIGN,INTSEL,IALIGN,
C    &              ETHRES,CTHRES,EROOT(100),NMS2VAL,MS2VAL(100)

*
      CALL QENTER('RASSG3')
*
      NTESTL = 0005
      NTEST = MAX(NTESTL,IPRSIG)
*
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' RASSG3 speaking :'
        WRITE(6,*) ' ================='
        write(6,*)
        WRITE(6,*) ' RASSG3 : NBATS = ',NBATS
      END IF
*
      LBLK = -1
      CALL REWINE(LLUC,-1)
      CALL REWINE(LLUHC,-1)
*
*. Loop over vectors of sigma and C : MK2 space + Real/Imag
      DO IRIS = 1, 2
        JBATABS = 0
        DO IMK2_S = 1, NMS2VAL
*
          MK2_S = MS2VAL(IMK2_S)
          IASTP = IST_FOR_DT(1,IMK2_S)
          IBSTP = IST_FOR_DT(2,IMK2_S)
*. Loop over batches of sigma for given IRIS and IMK2_S
          DO JBATS = 1, NBAT_MS2(IMK2_S)
            JBATABS = JBATABS + 1
            if (NTEST.ge.500) then
               WRITE(6,*)
               WRITE(6,*)
     &         ' =============================================='
               WRITE(6,*)
     &         ' Start of construction of Sigma batch ', JBATABS
               WRITE(6,*)
     &         ' =============================================='
               WRITE(6,*)
            end if
*. Zero sigma batch
            LS = LEBATS(JBATABS)
            ZERO = 0.0D0
            CALL SETVEC(SB,ZERO,LS)
*
            CALL REWINE(LLUC,-1)
            DO IRIC = 1, 2
              IF(IRIC.EQ.2) THEN
*. Skiup EOV mark between real and imaginary part
                 CALL IFRMDS(IONEM,1,-1,LLUC)
              END IF
              DO IMK2_C = 1, NMS2VAL
*
                IF(IRIC.EQ.2.AND.IMK2_C.EQ.NMS2VAL) THEN
                  I_DO_TRNSP = 1
                ELSE
                  I_DO_TRNSP = 0
                END IF
*
                MK2_C = MS2VAL(IMK2_C)
                IACTP = IST_FOR_DT(1,IMK2_C)
                IBCTP = IST_FOR_DT(2,IMK2_C)
*. Copy blocks with given MK2_C to LUSC37
                CALL REWINE(LUSC37,-1)
                NBL_C = NBLK_MS2(IMK2_C)
                CALL COPNBLKD(LLUC,LUSC37,CB,NBL_C,0,LBLK)
C                    COPNBLKD(LUIN,LUOUT,SEGMNT,NBLOCK,IREW,LBLK)
*. End of vector mark
                CALL ITODS(-1,1,-1,LUSC37)
                if (NTEST.ge.50) then
                  WRITE(6,*)
     &            ' ********************************************* '
                  WRITE(6,*)
     &            ' *                                           * '
                  WRITE(6,*)
     &            ' * Start of contribution from given CISPACE  * '
                  WRITE(6,*)
     &            ' *                                           * '
                  WRITE(6,*)
     &            ' ********************************************* '
                  WRITE(6,*)
                  WRITE(6,'(4X,A,2I3)')
     &            ' MK2/IRI vector of S :  MK2_S, IRIS =', MK2_S,IRIS
                  WRITE(6,'(4X,A,2I3)')
     &            ' MK2/IRI vector of C :  MK2_C, IRIC =', MK2_C,IRIC
                  WRITE(6,*)
                end if
*
C               rewind LUSC37
C               CALL WRTVCD(CB,LUSC37,1,-1)
                IF(ABS(MK2_S-MK2_C).LE.4) THEN
*
*. Connection. The Hamiltonian is at most a two-particle operator!
*. Obtain contribution to batch sigma from all blocks with MK2_C
CM                WRITE(6,*) ' Memchk before call to SBLOCK '
CM                CALL RMEMCHK
                  IF(NTEST.GE.100) THEN
                   WRITE(6,*) ' Block info for S batch, JBATS:',JBATABS
                   WRITE(6,*) ' Number of blocks : ', LBATS(JBATABS)
                   WRITE(6,*) ' Offset of batch ', I1BATS(JBATABS)
                   WRITE(6,*) ' TTSSOOLL array for blocks '
                   LLX = LBATS(JBATABS)
                   CALL IWRTMA(IBATS(1,I1BATS(JBATABS)),8,LLX,8,LLX)
                  END IF
*
                  CALL SBLOCK_REL(LBATS(JBATABS),
     &                            IBATS(1,I1BATS(JBATABS)),
     &                            1,
     &                            CB,SB,LUSC37,0,0,0,0,
     &                            0,IRIS,IRIC,I_DO_TRNSP)
*
                  if (NTEST.ge.100) then
                    WRITE(6,*) ' Updated array containing sblocks '
                    CALL WRTMAT(SB,1,LS,1,LS)
                  end if
                END IF
*               ^End if interacting
              END DO
*             ^End of loop over IMK2_C
            END DO
*           ^End of loop over real/imaginary part of C
*. A batch of sigma has been constructed so save
*. Transfer S block to permanent storage
            NTEST = 10
            if (NTEST.ge.10) then
               WRITE(6,*)
               WRITE(6,*) ' **********************************'
               WRITE(6,*) ' *                                *'
               WRITE(6,*) ' * Array containing final sbatch  *'
               WRITE(6,*) ' *                                *'
               WRITE(6,*) ' **********************************'
               WRITE(6,*)
               CALL WRTMAT(SB,1,LS,1,LS)
                  stop
            end if
            KBATSEND = I1BATS(JBATABS) + LBATS(JBATABS)
            DO ISBLK = I1BATS(JBATABS),KBATSEND-1
              IOFF = IBATS(6,ISBLK)
              LEN_T  = IBATS(8,ISBLK)
C?            write(6,*) 'RASSG3 : IOFF, SB(IOFF)',IOFF,SB(IOFF)
              CALL ITODS(LEN_T,1,-1,LLUHC)
              CALL TODSC(SB(IOFF),LEN_T,-1,LLUHC)
            END DO
          END DO
*         ^End of batches for given MK2_S
        END DO
*       ^End of loop over IMK2_S
*. End of real/imag part of sigma => end of vector
        CALL ITODS(-1,1,-1,LLUHC)
      END DO
*     ^End of loop over real/imag part of sigma
*
      IF(NTEST.GE.01) THEN
        WRITE(6,*) ' Final S-vector on disc'
        CALL WRTVCDC(SB,LLUHC,1,2,-1)
      END IF
*
      CALL QEXIT('RASSG3')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
      SUBROUTINE RSBB2A_REL(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NSCOL,
     &                  NGAS,ISOC,ICOC,
     &                  SB,CB,
     &                  ADSXA,DXSTST,STSTDX,SXDXSX,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                  NSMOB,NSMST,NSMSX,NSMDX,
     &                  RIKSX,RKJSX,MXSXST,MXSXBL,IMOC,SCLFAC,NTESTG,
     &                  NSEL2E,ISEL2E,IRI,IBU,SCR)
*
* two electron excitations on column strings
*
* Updated for double group calculation (IRI_INT, IBU added)
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NSCOL : Number of columns in S block
* ISEL1(3) : Number of electrons in RAS1(3) for S block
* ICEL1(3) : Number of electrons in RAS1(3) for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* DXSTST : Sym of sx,!st> => sym of sx !st>
* STSTDX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* SXDXSX : Symmetry of SX1,SX1*SX2 => symmetry of SX2
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals, strings,
*                           single excitations, double excitations
* MAXI   : Largest number of 'spectator strings' treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
*
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),DXSTST(NSMDX,NSMST)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
*
      INTEGER ISEL2E(*)
C     INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      DIMENSION CB(*)
      INTEGER ISOC(NGAS),ICOC(NGAS)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),XINT(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
      DIMENSION RIKSX(MXSXBL,4),RKJSX(MXSXBL,4)
      dimension SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
*
*.Local arrays
      DIMENSION ITP(256),JTP(256),KTP(256),LTP(256)
*
      INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
CTF   PARAMETER(MXTSOB = 40)
*
#include "symm.inc"
*
*   All important symmetry variables not on pointers
*
C     character*6 DBGSYM
C     integer   NIRR_PN,DDIM,NIRR_DG,DOUGRP,
C    &          SYDI,SYMDIM,ITOTDM,
C    &          IDIMGRP,NSDGIRR,NSPGIRR,
C    &          N_SPINR_EXP,NGASPI,NGSOBPI,IMAXDIM,
C    &          ITRIA,ITRIADB,ITRIASUM,ITRIADBSUM
C     dimension IDBGMULT(8,8),IPGMULT(4,4),IMOSF_SP(8,2),
C    &          ISPSF_MO(8,2),IADJSYM(8),INVELM(8)
C     dimension SYDI(8),SYMDIM(4),
C    &          NGSOBPI(8,16),
C    &          ITRIA(8),ITRIADB(4),ITOTDM(8)
C     common /DGSYMM/NIRR_PN,DDIM,NIRR_DG,DOUGRP,
C    &          SYDI,SYMDIM,ITOTDM,
C    &          IDIMGRP,NSDGIRR,NSPGIRR,
C    &          N_SPINR_EXP,NGASPI,NGSOBPI,IMAXDIM,
C    &          ITRIA,ITRIADB,ITRIASUM,ITRIADBSUM,
C    &          IDBGMULT,IPGMULT,IMOSF_SP,ISPSF_MO,
C    &          IADJSYM,INVELM

*
* DDIM	 	: double dimension (of SYMDIM(x))
* NIRR_DG	: number of double group irreps
* DOUGRP	: double group
* NIRR_PN	: number of spatial irreps
* SYDI		: symmetry dimensions of spatial irreps
* SYMDIM	: two SYDI's added (depending on double group)
* IDIMGRP	: dimension group for int trafo (1 SYMDIM)
* NSDGIRR	: number of sub double group irreps
* NSPGIRR	: number of sub point group irreps
* N_SPINR_EXP   : number of spinor expansions
* NGASPI	: number of GAS spaces (in picasso)
* NGSOBPI	: matrix with number orbitals in GAS and symmetries
* IMAXDIM	: largest dimension of spatial irreps (number from ordering
*                    in Koster, Dimmock, Wheeler, Statz)
* ITRIA         : triangle lengths of point group blocks
* ITRIASUM      : sum of these
* ITRIADB       : triangle lengths of two connected blocks in double group
* ITOTDM        : total orbital space dimensions per spatial irrep
* IDBGMULT      : Double Group multiplication table
* IPGMULT       : Point group mult. table
* IMOSF_SP      : Translation table from MO symmetry to spinor symmetry
*                 IMOSF_SP(IOBSM,IAB) Gives the spinor symmetry
*                 of an orbital of symmetry IOBSM times an alpha (IAB=1)
*                 or and beta spin-function (IAB = 2 )
* ISPSF_MO      : Translation table from spinor symmetry to MO symmetry
*                 ISPSF(ISPSM,IAB) gives the symmetry of an MO that
*                 combined with an alpha(IAB=1) or beta(IAB=2)
*                 spin function give a spinor with symmetry ISPSM
*
* IADJSYM       : IADJSYM(ISPSYM) is the adjoint symmetry of ISPSYM
* INVELM        : Inverse elements
*
* The irreps are organized in the usual way with the boson irreps first
*
      DIMENSION IACAR(2),ITPAR(2)
      CALL QENTER('RS2A')
*
      IMOC = 0
      NTESTL = 0000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' RSBB2A speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' ISCSM and ICCSM = ', ISCSM,ICCSM
        WRITE(6,*) ' IBU = ', IBU
C?      WRITE(6,*) ' MXSXST MXSXBL = ', MXSXST,MXSXBL
        WRITE(6,*) ' RSBB2A : IMOC ', IMOC
        WRITE(6,*) ' ISOC and ICOC : '
        CALL IWRTMA(ISOC,1,NGAS,1,NGAS)
        CALL IWRTMA(ICOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' NTESTG and NSEL2E', NTESTG, NSEL2E
      END IF
*
      IFRST = 1
      JFRST = 1
*
*.Types of DX that connects the two strings
*
C     IDXSM = STSTDX(ISCSM,ICCSM)
      IDXSM = IDBGMULT(ISCSM,INVELM(ICCSM))
      IF(IDXSM.EQ.0)  GOTO 2001
*. Connecting double excitations
      CALL DXTYP_GAS(NDXTYP,ITP,JTP,KTP,LTP,NGAS,ISOC,ICOC)
      DO 2000 IDXTYP = 1, NDXTYP
        ITYP = ITP(IDXTYP)
        JTYP = JTP(IDXTYP)
        KTYP = KTP(IDXTYP)
        LTYP = LTP(IDXTYP)
        NIJKL1 = 0
        IF(ITYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(JTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(KTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(LTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(NIJKL1.EQ.0) CALL QENTER('BB2A0')
        IF(NIJKL1.EQ.1) CALL QENTER('BB2A1')
        IF(NIJKL1.EQ.2) CALL QENTER('BB2A2')
        IF(NIJKL1.EQ.3) CALL QENTER('BB2A3')
        IF(NIJKL1.EQ.4) CALL QENTER('BB2A4')
*. Integrals included ??
        IF(NSEL2E.NE.0) THEN
          IAMOKAY=0
          IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP) THEN
            DO JSEL2E = 1, NSEL2E
              IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY=1
            END DO
          END IF
          IF(IAMOKAY.EQ.0) GOTO 2000
        END IF
*. Largest possible number of orbital pairs
        MI = 0
        MJ = 0
        MK = 0
        ML = 0
        DO IOBSM = 1, NSMST
          MI = MAX(MI,NOBPTS(ITYP,IOBSM))
          MJ = MAX(MJ,NOBPTS(JTYP,IOBSM))
          MK = MAX(MK,NOBPTS(KTYP,IOBSM))
          ML = MAX(ML,NOBPTS(LTYP,IOBSM))
        END DO
        MXPAIR = MAX(MI*MK,MJ*ML)
C       WRITE(6,*) ' MXPAIR  = ', MXPAIR
*. Largest posssible
*. Symmetry of allowed Double excitation, loop over excitations
        DO 1950 IKOBSM = 1, NSMOB
C         JLOBSM = SXDXSX(IKOBSM,IDXSM)
          JLOBSM = IDBGMULT(IDXSM,INVELM(IKOBSM))
C         WRITE(6,*) ' IDXSM IKOBSM JLOBSM ',IDXSM,IKOBSM,JLOBSM
          IF(JLOBSM.EQ.0) GOTO 1950
*. types + symmetries defined => K strings are defined
          KFRST = 1
          K2FRST = 1
*
*. Number of batchs of symmetry pairs IK
*
          LENGTH = 0
          NIKBT = 0
          NBLK = 0
          NBLKT = 0
          DO ISM = 1, NSMOB
C           KSM = ADSXA(ISM,IKOBSM)
            KSM = IDBGMULT(IKOBSM,INVELM(ISM))
            NI = NOBPTS(ITYP,ISM)
            NK = NOBPTS(KTYP,KSM)
            IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
              NIK = NI*(NI+1)/2
C           ELSE IF(ITYP.GT.KTYP.OR.ITYP.EQ.KTYP.AND.ISM.GT.KSM) THEN
            ELSE IF( IOBPTS(ITYP,ISM).GT.IOBPTS(KTYP,KSM)) THEN
              NIK = NI*NK
            ELSE
              NIK = 0
            END IF
            IF(NIK.NE.0) THEN
              NBLKT = NBLKT + 1
              IF(LENGTH+NIK .GT. MXPAIR) THEN
*. The present batch is complete
                NIKBT = NIKBT+1
                IKBT(1,NIKBT) = NBLKT - NBLK
                IKBT(2,NIKBT) = NBLK
                IKBT(3,NIKBT) = LENGTH
                LENGTH = 0
                NBLK = 0
              END IF
              NBLK = NBLK + 1
              LENGTH = LENGTH + NIK
              IKSMBT(1,NBLKT) = ISM
              IKSMBT(2,NBLKT) = KSM
            END IF
          END DO
*. The last batch
          IF(NBLK.NE.0) THEN
            NIKBT = NIKBT+1
            IKBT(1,NIKBT) = NBLKT - NBLK + 1
            IKBT(2,NIKBT) = NBLK
            IKBT(3,NIKBT) = LENGTH
          END IF
*.
          IF(NTEST.GE.2000) THEN
            WRITE(6,*) ' ITYP, KTYP, IKOBSM,  NIKBT = ',
     &                   ITYP, KTYP, IKOBSM,  NIKBT
            WRITE(6,*) ' IKBT : Offset, number, length '
            DO JIKBT = 1, NIKBT
              WRITE(6,'(3i3)') (IKBT(II,JIKBT), II = 1, 3)
            END DO
            WRITE(6,*) ' IKSMBT '
            CALL IWRTMA(IKSMBT,2,NBLKT,2,8)
          END IF
*
*. Number of batchs of symmetry pairs JL
*
          LENGTH = 0
          NJLBT = 0
          NBLK = 0
          NBLKT = 0
          DO JANISM = 1, NSMOB
C           LSM = ADSXA(JSM,JLOBSM)
            LANISM = IDBGMULT(JLOBSM,INVELM(JANISM))
*. JANI and LANI are annihilation symmetries, creation symmetries
            JSM = IADJSYM(JANISM)
            LSM = IADJSYM(LANISM)
C           WRITE(6,*) 'JLOBSM,JANISM,LANISM, JSM,LSM',
C    &                 JLOBSM,JANISM,LANISM, JSM,LSM

*
            NJ = NOBPTS(JTYP,JSM)
            NL = NOBPTS(LTYP,LSM)
            IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
              NJL = NJ*(NJ+1)/2
C           ELSE IF(JTYP.GT.LTYP.OR.JTYP.EQ.LTYP.AND.JSM.GT.LSM) THEN
            ELSE IF (IOBPTS(JTYP,JSM).GT.IOBPTS(LTYP,LSM) ) THEN
              NJL = NJ*NL
            ELSE
              NJL = 0
            END IF
            IF(NJL.NE.0) THEN
              NBLKT = NBLKT + 1
              IF(LENGTH+NJL .GT. MXPAIR) THEN
*. The present batch is complete
                NJLBT = NJLBT+1
                JLBT(1,NJLBT) = NBLKT - NBLK
                JLBT(2,NJLBT) = NBLK
                JLBT(3,NJLBT) = LENGTH
                LENGTH = 0
                NBLK = 0
              END IF
              NBLK = NBLK + 1
              LENGTH = LENGTH + NJL
              JLSMBT(1,NBLKT) = JSM
              JLSMBT(2,NBLKT) = LSM
            END IF
          END DO
*. The last batch
          IF(NBLK.NE.0) THEN
            NJLBT = NJLBT+1
            JLBT(1,NJLBT) = NBLKT - NBLK + 1
            JLBT(2,NJLBT) = NBLK
            JLBT(3,NJLBT) = LENGTH
          END IF

*.
          IF(NTEST.GE.2000) THEN
            WRITE(6,*) ' JTYP, LTYP, JLOBSM,  NJLBT = ',
     &                   JTYP, LTYP, JLOBSM,  NJLBT
            WRITE(6,*) ' JLBT : Offset, number, length '
            DO JJLBT = 1, NJLBT
              WRITE(6,'(3i3)') (JLBT(II,JJLBT), II = 1, 3)
            END DO
            WRITE(6,*) ' JLSMBT '
            CALL IWRTMA(JLSMBT,2,NBLKT,2,8)
          END IF
*
*. Loop over batches of IK strings
          DO 1940 IKBTC = 1, NIKBT
          IF(NTEST.GE.1000) WRITE(6,*) ' IKBTC = ', IKBTC
*. Loop over batches of JL strings
          DO 1930 JLBTC = 1, NJLBT
*
CMOVED        KFRST = 1
CMOVED        K2FRST = 1
*
              IF(IMOC.EQ.0) THEN
*
* =========================================================================
*                    Use N-2 projection method
* =========================================================================
*
              IFIRST = 1
*. Loop over batches of I strings
              NPART = NROW/MAXI
              IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
              IF(NTEST.GE.2000)
     &        write(6,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
              DO 1801 IPART = 1, NPART
                IBOT = 1+(IPART-1)*MAXI
                ITOP = MIN(IBOT+MAXI-1,NROW)
                NIBTC = ITOP-IBOT+1
*.Loop over batches of intermediate strings
                KBOT = 1- MAXK
                KTOP = 0
 1800           CONTINUE
                  KBOT = KBOT + MAXK
                  KTOP = KTOP + MAXK
*
                  IONE = 1
                  JLBOFF = 1
                  NJLT = JLBT(3,JLBTC)
                  DO JLPAIR = 1, JLBT(2,JLBTC)
                    JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                    LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
                    NJ = NOBPTS(JTYP,JSM)
                    NL = NOBPTS(LTYP,LSM)
                    IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                      NJL = NJ*(NJ+1)/2
                      JLSM = 1
                    ELSE
                      NJL = NJ * NL
                      JLSM = 0
                    END IF
C                   WRITE(6,*) ' JSM LSM JLSM ' , JSM, LSM,JLSM

*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Obtain all double excitations from this group of K strings
CT                CALL QENTER('ADADS')
                  I12C = 1
                  K12C = 1
                  IONE = 1
                  CALL ADADST_GAS_REL(IONE,JSM,JTYP,NJ,
     &                            IONE,LSM,LTYP,NL,
     &                        ICCTP,ICCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                        JFRST,KFRST,I12C,K12C,SCLFAC)
                  JFRST = 0
                  KFRST = 0
*
CT                CALL QEXIT('ADADS')
C                 WRITE(6,*) ' NKBTC after ADADST_GAS',NKBTC
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  J = 0
                  L = 1
*
CT                CALL QENTER('MATCG')
                  DO  IJL = 1, NJL
                    CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                    I1JL = (L-1)*NJ+J
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                    JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                    IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                      ZERO = 0.0D0
                      CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                    ELSE
                      CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,
     &                            I1(1,I1JL),XI1S(1,I1JL))
                    END IF
                  END DO
CT                CALL QEXIT ('MATCG')
*
                  JLBOFF = JLBOFF + NJL
                  END DO
*. ( End of loop over jl in batch )
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
*.Obtain two electron integrals xint(ik,jl) = (ij!kl)-(il!kj)
                  IF(IFIRST.EQ.1) THEN
                    IXCHNG = 1
* Obtain integrals in ik batch
                    NIKT = IKBT(3,IKBTC)
                    NJLT = JLBT(3,JLBTC)
                    JLOFF = 1
                    DO JLPAIR = 1, JLBT(2,JLBTC)
                    IKOFF = 1
                    DO IKPAIR = 1, IKBT(2,IKBTC)
*
                      ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                      KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                      JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                      LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
*
                      IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                        IKSM = 1
                        NIK = NOBPTS(ITYP,ISM)*(NOBPTS(ITYP,ISM)+1)/2
                      ELSE
                        IKSM = 0
                        NIK = NOBPTS(ITYP,ISM)*NOBPTS(KTYP,KSM)
                      END IF
*
                      IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                        JLSM = 1
                        NJL = NOBPTS(JTYP,JSM)*(NOBPTS(JTYP,JSM)+1)/2
                      ELSE
                        JLSM = 0
                        NJL = NOBPTS(JTYP,JSM)*NOBPTS(LTYP,LSM)
                      END IF

                      ICOUL = 0
                      CALL GETINT_DBG(SCR,ITYP,ISM,IBU,JTYP,JSM,IBU,
     &                                KTYP,KSM,IBU,LTYP,LSM,IBU,
     &                                IRI,IXCHNG,IKSM,JLSM,ICOUL)

C                     CALL GETINT(SCR,ITYP,ISM,JTYP,JSM,KTYP,KSM,LTYP,
C    &                          LSM,IXCHNG,IKSM,JLSM,ICOUL)
                      DO JL = 1, NJL
                        CALL COPVEC(SCR((JL-1)*NIK+1),
     &                             XINT((JLOFF-1+JL-1)*NIKT+IKOFF),NIK)
                      END DO
                      IKOFF = IKOFF + NIK
                    END DO
                    JLOFF = JLOFF + NJL
                    END DO
                  END IF
                  IFIRST = 0
*.and now ,to the work
                  LIKB = NIBTC*NKBTC
                  IF(NTEST.GE.3000) THEN
                   WRITE(6,*) ' Integral block '
                   CALL WRTMAT(XINT,NIKT,NJLT,NIKT,NJLT)
                  END IF
                  IF(NTEST.GE.3000) THEN
                    WRITE(6,*) ' CSCR matrix '
                    CALL WRTMAT(CSCR,LIKB,NJLT,LIKB,NJLT)
                  END IF
*
                  FACTORC = 0.0D0
                  FACTORAB = 1.0D0
                  CALL MATML7(SSCR,CSCR,XINT,
     &                        LIKB,NIKT,LIKB,NJLT,NIKT,NJLT,
     &                        FACTORC,FACTORAB,2)
                  IF(NTEST.GE.3000) THEN
                    WRITE(6,*) ' SSCR matrix '
                    CALL WRTMAT(SSCR,LIKB,NIKT,LIKB,NIKT)
                    WRITE(6,*) ' Mem check '
                    CALL RMEMCHK
                  END IF
* ============================
* Loop over ik and scatter out
* ============================
*. Generate double excitations from K strings
*. I strings connected with K strings in batch <I!a+i a+k!K)
                  IF(NTEST.GE.3000) THEN
                    WRITE(6,*) ' S block second call to ADADST_GAS'
                    CALL WRTMAT(SB,NROW,NSCOL,NROW,NSCOL)
                  END IF
                  I12C = 2
*
                  IONE = 1
                  IKBOFF = 1
                  DO IKPAIR = 1, IKBT(2,IKBTC)
                    ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                    KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                    NI = NOBPTS(ITYP,ISM)
                    NK = NOBPTS(KTYP,KSM)
                    IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                      NIK = NI*(NI+1)/2
                      IKSM = 1
                    ELSE
                      NIK = NI * NK
                      IKSM = 0
                    END IF
CT                  CALL QENTER('ADADS')
                    IF(IFRST.EQ.1) KFRST = 1
                    ONE = 1.0D0
                    CALL ADADST_GAS_REL(IONE,ISM,ITYP,NI,
     &                              IONE,KSM,KTYP,NK,
     &                              ISCTP,ISCSM,IGRP,
     &                              KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                              IFRST,KFRST,I12C,K12C,ONE)
*
                    IFRST = 0
                    K2FRST = 0
CT                  CALL QEXIT ('ADADS')
*
CT                  CALL QENTER('MATCS')
                    I = 0
                    K = 1
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' S block before MATCA '
                      CALL WRTMAT(SB,NROW,NSCOL,NROW,NSCOL)
                    END IF
                    DO IK = 1, NIK
                      CALL NXTIJ(I,K,NI,NK,IKSM,NONEW)
                      IKOFF = (K-1)*NI + I
                      ISBOFF = 1+(IKBOFF-1+IK-1)*NIBTC*NKBTC
                      IF(IKSM.EQ.1.AND.I.EQ.k) THEN
* a+ i a+i gives trivially zero
                      ELSE
                        CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                              NKBTC,I1(1,IKOFF),XI1S(1,IKOFF))
                      END IF
                    END DO
CT                  CALL QEXIT ('MATCS')
                    IKBOFF = IKBOFF + NIK
                    IF(NTEST.GE.3000) THEN
                     WRITE(6,*) ' Updated S block'
                     CALL WRTMAT(SB,NROW,NSCOL,NROW,NSCOL)
                    END IF
                  END DO
*
                IF(KEND.EQ.0) GOTO 1800
*. End of loop over partitionings of resolution strings
 1801         CONTINUE
                ELSE IF (IMOC.EQ.1) THEN
*
* =========================================================================
*                    Use MOC method
* =========================================================================
*
* Fetch two electron integrals as (i,k,j,l)
                  IXCHNG = 1
                    ICOUL = 0
                    CALL GETINT_REL(XINT,
     &                              ITYP,ISM,JTYP,JSM,KTYP,KSM,LTYP,
     &                              LSM,IXCHNG,0,0,ICOUL)


*. Generate off diagonal I => K single excitations
              CALL SXSTR(ISCSM,ISCTP,IGRP,LSM,LTYP,
     &                   KSM,KTYP,MXSXST,RIKSX(1,1),RIKSX(1,2),
     &                   RIKSX(1,3),RIKSX(1,4),I1,1,NTEST)
*. Generate off-diagonal K => J single excitations
      IACAR(1) = 2
      IACAR(2) = 1
*. We are actually using conjugated operators so
      ITPAR(1) = LTYP
      ITPAR(2) = KTYP
*     ITPAR(1) = KTYP
*     ITPAR(2) = LTYP
CSK - 6th of October 2009 - commented out because of wrong number of
Carguments in calling list - thanks to Miro for pointig it out
CSK   CALL NEWTYP_REL(ISCTP,IACAR,ITPAR,2,KSTTYP)
C     NEWTYP_REL(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
C      SYMCOM(ITASK,IOBJ,I1,I2,I12)
      CALL SYMCOM_REL(3,0,KLSM,ISCSM,KSTRSM)
C?    WRITE(6,*) ' KSTRSM KSTTYP = ', KSTRSM,KSTTYP
*. Generate off diagonal K=> J single excitations
              CALL SXSTR(KSTRSM,KSTTYP,IGRP,JSM,JTYP,
     &                   ISM,ITYP,MXSXST,RKJSX(1,1),RKJSX(1,2),
     &                   RKJSX(1,3),RKJSX(1,4),XI1S,1,NTEST)
C?     write(6,*) ' memory check in after second SXSTR '
C?    write(6,*) ' MXSXST 3 ', MXSXST
C?     call memchk
*. Combine single excitations to double excitations
*. Loop over columns trings of sigma
       DO ISTR = 1, NSCOL
*. Double excitations out from DXFSTR
       CALL DXFSTR(RIKSX(1,1),RIKSX(1,2),RIKSX(1,3),RIKSX(1,4),
     &             I1,
     &             RKJSX(1,1),RKJSX(1,2),RKJSX(1,3),RKJSX(1,4),
     &             XI1S,
     &             NI,NJ,NK,NL,
     &             IOFF,JOFF,KOFF,LOFF,
     &             ISTR,MXSXST,XINT,
     &             I2,XI2S,NEX,NTEST)
C?     write(6,*) ' memory check in after DXSTR        '
C?     call memchk
C?    write(6,*) ' MXSXST 4 ', MXSXST
*. Msaxpy update
        CALL MSAXPYN6(SB((ISTR-1)*NROW+1),XI2S,CB,NROW,
     &                NEX,I2)
C      write(6,*) ' memory check in RSBB2A '
C      call memchk
        END DO
              END IF
 1930       CONTINUE
 1940     CONTINUE
 1950   CONTINUE
        IF(NIJKL1.EQ.0) CALL QEXIT ('BB2A0')
        IF(NIJKL1.EQ.1) CALL QEXIT ('BB2A1')
        IF(NIJKL1.EQ.2) CALL QEXIT ('BB2A2')
        IF(NIJKL1.EQ.3) CALL QEXIT ('BB2A3')
        IF(NIJKL1.EQ.4) CALL QEXIT ('BB2A4')
 2000 CONTINUE
*
 2001 CONTINUE
      CALL QEXIT('RS2A ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE RSBB2B_REL(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                  JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,NGAS,
     &                  IAOC,IBOC,JAOC,JBOC,
     &                  SB,CB,
     &                  ADSXA,STSTSX,
     &                  NOBPTS,IOBPTS,NOBPTS2,IOBPTS2,
     &                  ITSOB,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  XINT,
     &                  NSMOB,NSMST,NSMSX,NSMDX,IUSEAB,
     &                  ICJKAIB,CJRES,SIRES,S2,SCLFAC,NTESTG,
     &                  NSEL2E,ISEL2E,IRI,SCR)
*
* Combined alpha-beta double excitation
* contribution from given C block to given S block
*. If IUSAB only half the terms are constructed
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha  strings in sigma
* IBSM,IBTP : Symmetry and type of beta   strings in sigma
* JASM,JATP : Symmetry and type of alpha  strings in C
* JBSM,JBTP : Symmetry and type of beta   strings in C
* NIA,NIB : Number of alpha-(beta-) strings in sigma
* NJA,NJB : Number of alpha-(beta-) strings in C
* IAGRP : String group of alpha strings
* IBGRP : String group of beta strings
* IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
* IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
* JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
* JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
* CB   : Input C block
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
*       single excitations
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ICJKAIB =1 =>  construct C(Ka,Jb,j) and S(Ka,Ib,i) as intermediate
*                 matrices in order to reduce overhead
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + ICJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks
*
*. October 96 : New routines for accessing annihilation information
*.              Cleaned and shaved, only IROUTE = 3 option active
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),
     &        NOBPTS2(MXPNGAS,*),IOBPTS2(MXPNGAS,*),
     &        ITSOB(*)
*
      INTEGER ISEL2E(*)
C     INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION S2(*)
      dimension SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
*.Local arrays
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS),KTP(MXPNGAS),LTP(MXPNGAS)
      INTEGER IREO(4),IILEN(4)
*
CTF
* Added these nice things ...
#include "symm.inc"
*
      CALL QENTER('RS2B ')
      NTESTL = 0
      NTEST = MAX(NTESTL,NTESTG)
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' RSBB2B speaking '
        WRITE(6,*) ' =============== '
        WRITE(6,*)
        WRITE(6,*) ' Initial sigma block '
        CALL WRTMAT(SB,NIB,NIA,NIB,NIA)

      END IF

      IROUTE = 3
* IROUTE = 1 : Normal (i.e. old) route,
* IROUTE = 2 : New route with j first
* IROUTE = 3 : C(Ka,j,Jb)
*
*
*. Symmetry of allowed excitations
      IJSM = IDBGMULT(IADJSYM(IASM),JASM)
      KLSM = IDBGMULT(IADJSYM(IBSM),JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP_GAS_REL(NKLTYP,KTP,LTP,NGAS,MXPNGAS,IBOC,JBOC)
      CALL SXTYP_GAS_REL(NIJTYP,ITP,JTP,NGAS,MXPNGAS,IAOC,JAOC)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
*
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = IDBGMULT(ISM,IADJSYM(IJSM))
          IF(NTESTL.GE.2000)
     &    WRITE(6,*) ' Loop 1940, ISM, KSM ', ISM,JSM
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          if(ntest.ge.1500) write(6,*) ' ISM JSM ', ISM,JSM
          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Generate annihilation mappings for all K strings
*. a+j!ka> = +/-/0 * !Ja>
          CALL ADSTN_GAS_REL(JSM,JTYP,JATP,JASM,IAGRP,
     &                   I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFAC)
*. a+i!ka> = +/-/0 * !Ia>
          ONE = 1.0D0
          CALL ADSTN_GAS_REL(ISM,ITYP,IATP,IASM,IAGRP,
     &                   I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE)
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
C             COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
              CALL COMPRS2LST(I1,XI1S,NJ,I3,XI3S,NI,NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF

*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
C         IF(NKABTC.GT.1) THEN
C           WRITE(6,*) ' RSBB2B : NKABTC .GT. 1 '
C           WRITE(6,*) ' I am not prepared for this '
C           STOP 'RSBB2B :  NKABTC .GT. 1'
C         END IF
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, NJ
              CALL GET_CKAJJB(CB,NJ,NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*
            ZERO = 0.0D0
            CALL SETVEC(SIRES,ZERO,NIB*LKABTC*NI)
            FACS = 1.0D0
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
*. Should this group of excitations be included
              IF(NSEL2E.NE.0) THEN
               IAMOKAY=0
               IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP)THEN
                 DO JSEL2E = 1, NSEL2E
                   IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY = 1
                 END DO
               END IF
               IF(IAMOKAY.EQ.0) GOTO 2000
              END IF
*
              DO 1930 KSM = 1, NSMOB
                IFIRST = 1
                LSM = IDBGMULT(KSM,IADJSYM(KLSM))
                IF(NTEST.GE.2000)
     &          WRITE(6,*) ' Loop 1930, KSM, LSM = ',KSM,LSM
                IF(LSM.EQ.0) GOTO 1930
                KOFF = IOBPTS2(KTYP,KSM)
                LOFF = IOBPTS2(LTYP,LSM)
                NK = NOBPTS2(KTYP,KSM)
                NL = NOBPTS2(LTYP,LSM)
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
                CALL ADSTN_GAS_REL(LSM,LTYP,JBTP,JBSM,IBGRP,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADSTN_GAS_REL(KSM,KTYP,IBTP,IBSM,IBGRP,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*
* Fetch Integrals as (j i k l )
*
*. We want the integrals out in the order (J I ! K L).
*. However, for complex integrals we cannot just switch the
* indeces around, so obtain the correct integrals and
* reorder subsequently
*
* 1 : ( I J ! K_ L_)
                IXCHNG = 0
                ICOUL = 1
                CALL GETINT_DBG(SCR,ITYP,ISM,2,JTYP,JSM,2,KTYP,KSM,1,
     &           LTYP,LSM,1,IRI,IXCHNG,0,0,ICOUL)
*. Reorder to ( J I ! K L)
                 IREO(1) = 2
                 IREO(2) = 1
                 IREO(3) = 3
                 IREO(4) = 4
                 IILEN(1) = NI
                 IILEN(2) = NJ
                 IILEN(3) = NK
                 IILEN(4) = NL
C      REO4MT(XIN,XOUT,LIN,IREO,FACIN,FACOUT)
                 ONE = 1.0D0
                 ZERO = 0.0D0
                 CALL REO4MT(SCR,XINT,IILEN,IREO,ONE,ZERO)
*. 2 : (K_ J !L I_)
C                CALL GETINT_DBG(SCR,KTYP,KSM,1,JTYP,JSM,2,LTYP,LSM,1,
C    &                           ITYP,ISM,2,IRI,IXCHNG,0,0,ICOUL)
*. Reorder to (J I ! K L)
C                IREO(1) = 3
C                IREO(2) = 1
C                IREO(3) = 4
C                IREO(4) = 2
C                IILEN(1) = NK
C                IILEN(2) = NJ
C                IILEN(3) = NL
C                IILEN(4) = NI
C                CALL REO4MT(SCR,XINT,IILEN,IREO,ONE,ONE)
*
                 if (NTEST.ge.1000) then
                   WRITE(6,'(A,8I3)') ' Sym and type for ijkl ',
     &             ISM,ITYP,JSM,JTYP,KSM,KTYP,LSM,LTYP
                   WRITE(6,*) ' XINT as (J I ! K L ) '
                   CALL WRTMAT(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
                 end if
*
* S(Ka,j,Ib) = sum(k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)
*
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' C(Ka,j,JB) as C(Kaj,JB) '
                  CALL WRTMAT(CJRES,LKABTC*NJ,NJB,LKABTC*NJ,NJB)
                END IF
*
                CALL SKICKJ_REL(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,XINT,NI,NJ,NK,NL,NKBSTR,
     &               I4,XI4S,I2,XI2S,IKORD,
     &               FACS,IROUTE,MXTSOB)
                FACS = 1.0D0
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Updated S(Ka,i,IB) as S(Kai,IB) '
                  CALL WRTMAT(SIRES,LKABTC*NI,NIB,LKABTC*NI,NIB)
                END IF
*
*
 1930         CONTINUE
 2000       CONTINUE
*
*. Scatter out from s(Ka,Ib,i)
*
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Final S(Ka,i,IB) as S(Kai,IB) '
                  CALL WRTMAT(SIRES,LKABTC*NI,NIB,LKABTC*NI,NIB)
                END IF
            DO II = 1, NI
              CALL ADD_SKAIIB(SB,NI,NIA,SIRES,LKABTC,NIB,II,
     &             I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
 1801     CONTINUE
*. End of loop over partitioning of alpha strings
 1940   CONTINUE
 2001 CONTINUE
*
 9999 CONTINUE
*
*
      CALL QEXIT('RS2B ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE RSSBCBREL(ONEINTS,INTSPACE,
     &     IASM,IATP,IOCTPA,IAOC,IAGRP,NIA,
     &     IBSM,IBTP,IOCTPB,IBOC,IBGRP,NIB,
     &     JASM,JATP,JOCTPA,JAOC,JAGRP,NJA,
     &     JBSM,JBTP,JOCTPB,JBOC,JBGRP,NJB,
     &     NGAS,SB,CB,IDOH2,
     &     ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &     NOBPTS,IOBPTS,NOBPTS2,IOBPTS2,
     &     MXPNGAS,MAXI,MAXK,
     &     SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &     NSMOB,NSMST,NSMSX,NSMDX,
     &     MXPOBS,IDC,
     &     CJRES,SIRES,I3,XI3S,I4,XI4S,
     &     MXSXBL,MXSXST,IPRNT,ISPINFREE,
     &     IPTSPC,JPTSPC,IJOP,SCLFAC,
     &     IRI,IRICO,IELCI,
     &     MXTSOB,SOMESCR,XIJILS,ISMDST)
*
*     Contributions to sigma block (iasm iatp, ibsm ibtp ) from
*     C block (jasm jatp , jbsm, jbtp)
*
*
*     Cleaned and shaved by Jeppe + Timo, March19, 1998
*     =====
*     Input
*     =====
*
*     IASM,IATP : Symmetry and type of alpha strings in sigma
*     IBSM,IBTP : Symmetry and type of beta  strings in sigma
*     JASM,JATP : Symmetry and type of alpha strings in C
*     JBSM,JBTP : Symmetry and type of beta  strings in C
*     NGAS      : Number of active spaces in calculation
*     IAOC,IBOC : Number of electrons in each AS for sigma supergroups
*     JAOC,JBOC : Number of electrons in each AS for C     supergroups
*     NAEL : Number of alpha electrons
*     NBEL : Number of  beta electrons
*     IJAGRP    : IA and JA belongs to this group of strings
*     IJBGRP    : IB and JB belongs to this group of strings
*     CB : Input c block
*     IDOH2 : = 0 => no two electron operator
*     IDOH2 : = 1 =>    two electron operator
*     ADASX : sym of a+, a => sym of a+a
*     ADSXA : sym of a+, a+a => sym of a
*     SXSTST : Sym of sx,!st> => sym of sx !st>
*     STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*     is nonvanishing by symmetry
*     DXSTST : Sym of dx,!st> => sym of dx !st>
*     STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*     is nonvanishing by symmetry
*     NTSOB  : Number of orbitals per type and symmetry
*     IBTSOB : base for orbitals of given type and symmetry
*     IBORB  : Orbitals of given type and symmetry
*     MAXI   : Largest Number of ' spectator strings 'treated simultaneously
*     MAXK   : Largest number of inner resolution strings treated at simult.
*
CTF
CTF   We need separate calls for real/imaginary parts according to
CTF   IRI: r/i part of sigma block
CTF   IRICO: r/i part of C block
CTF   (the appropriate integral list is then determined automatically)
CTF
CTF   IUB is given directly before routines are called for update of S
CTF   same for IPM1
CTF
*
*     ======
*     Output
*     ======
*     SB : fresh sigma block
*
*     =======
*     Scratch
*     =======
*     SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*     largest number of orbital pairs of given symmetries and
*     types.
*     I1, XI1S   : at least MXSTSO : Largest number of strings of given
*     type and symmetry
*     I1, XI1S   : at least MXSTSO : Largest number of strings of given
*     type and symmetry
*     C2 : Must hold largest STT block of sigma or C
*
*     XINT : Scratch space for integrals.
*
*     Jeppe Olsen , Winter of 1991
*
*     modified for relativistic double group calculations
*     Timo Fleig + Jeppe Olsen, March 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      REAL*8 INPROD
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
*
CTF DELTA_MK2 should be an integer!!
      integer DELTA_MK2
      INTEGER JKL_AC(3)
*
*     Input
      dimension CB(*)
*     Input (from WORK(KPICI1) in ... sblock)
      dimension ONEINTS(INTSPACE)
      INTEGER ISMDST(*)
*
*     . Output
      DIMENSION SB(*)
*     . Scratch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*)
      DIMENSION C2(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION IBLOCK(8)
*
      call qenter('RSSBCBREL')
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      XNORM = INPROD(CB,CB,NJA*NJB)
C     IF(XNORM.EQ.0.0D0) GOTO 1001
*
      if (NTEST.ge.100) call rmemchk
*
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' ====================='
         WRITE(6,*) ' RSSBCBREL :  C block '
         WRITE(6,*) ' ====================='
         CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
         WRITE(6,*) ' ============================='
         WRITE(6,*) ' RSSBCBREL : Initial  S block '
         WRITE(6,*) ' ============================='
         CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
         WRITE(6,*) ' Overall scalefactor ',SCLFAC
      END IF
*
      if (NTEST.ge.2000) then
        write(6,*)
        write(6,*) 'Checking integrals in rssbcbrel '
        do INT=1,INTSPACE,1
          write(6,'(I7,1F18.14)') INT,ONEINTS(INT)
        end do
      end if
*
      IF(NTEST.GE.500) THEN
         WRITE(6,*) ' IAOC and IBOC '
         CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
         CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
         WRITE(6,*) ' JAOC and JBOC  : '
         CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
         CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
         WRITE(6,*) ' IASM IATP JASM JATP ', IASM,IATP,JASM,JATP
         WRITE(6,*) ' IBSM IBTP JBSM JBTP ', IBSM,IBTP,JBSM,JBTP
         WRITE(6,*) ' NJA NJB NIA NIB ', NJA,NJB,NIA,NIB
      END IF
*     . No use of alpha-beta symmetry at the moment
      IUSEAB = 0
CTF
* Should the corresponding Hamiltonian matrix block be
* calculated exactly or approximately
*       NSEL2E = NNSEL2E
*  This is how it could be done somewhen (by passing NNSEL2E down)
*. Operator specified by input (for now, all excitations included)
      NSEL2E = 0
CTF
*     . Number of alpha- and beta- electrons in strings
      NIAEL = IELSUM(IAOC,NGAS)
      NIBEL = IELSUM(IBOC,NGAS)
      NJAEL = IELSUM(JAOC,NGAS)
      NJBEL = IELSUM(JBOC,NGAS)
*
      NIEL = NIAEL + NIBEL
      NJEL = NJAEL + NJBEL
*
      IF(NIEL.NE.NJEL) THEN
         WRITE(6,*)
     &        ' Total number of electrons in IaIb and JaJb differs'
         WRITE(6,*) ' NIAEL NIBEL NJAEL NJBEL', NIAEL,NIBEL,NJAEL,NJBEL
         Call Abend2
     &   ('RSSBCBREL: Tot. number of els. in IaIb and JaJb differs')
      END IF
*     . Difference in unbarred (A) and barred (B) occupation
*     ( number in sigma - number in C)
      DELTA_MK2 = 2*(NIAEL - NJAEL)
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' DELTA_MK2 = ', DELTA_MK2
      END IF
*
      IF (DELTA_MK2.EQ.0) THEN
         if (NTEST.ge.100) then
            write(6,*) '.........................'
            write(6,*) 'Entering DMK2=0 part ...'
            write(6,*) '.........................'
         end if
*
*     =================================
*     Sigma (delta Mk = 0) contribution
*     =================================
*
*     Sigma bb(IA,IB) = sum(i.gt.k,j.gt.l)<IB!Eb(ij)Eb(kl)!JB>
*     * ((ij!kl)-(il!kj)) C(IA,JB)
*     + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)
*     .One electron part beta
         JACTIVE = 0
         IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.NIBEL.GE.1) THEN
            if (NTEST.ge.100) then
               write(6,*) '........................'
               write(6,*) 'I have barred electrons'
               write(6,*) '........................'
            end if
            JACTIVE = 1
            CALL TRPMT3(SB,NIB,NIA,C2)
            CALL COPVEC(C2,SB,NIA*NIB)
            CALL TRPMT3(CB,NJB,NJA,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
*
            IF(NTEST.GE.500) THEN
               WRITE(6,*) ' SB before dmk0one'
               call wrtmat(sb,nia,nib,nia,nib)
            END IF
            IF(NTEST.GE.101)
     &           WRITE(6,*) ' I am going to call dmk0one'
*
*     Call new relativistic routines
C?          WRITE(6,*) 'Do memory check ...before dmk0one'
            if (NTEST.ge.100) call rmemchk
            IUB = 2
            IPM1 = 0
CSKold      call dmk0one(IBSM,IBTP,JBSM,JBTP,
CSKold           IBGRP,NIA,
CSKold           NGAS,IBOC,JBOC,
CSKold           SB,CB,
CSKold           ADSXA,STSTSX,
CSKold           NOBPTS2,IOBPTS2,
CSKold           MAXI,MAXK,
CSKold           SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
CSKold           NSMOB,NSMST,MOCAA,
CSKold           SCLFAC,
CSKold           ONEINTS,INTSPACE,
CSKold           IRI,IRICO,IPM1,IUB,ISMDST,NSMST)
*
            IF(NTEST.GE.500) THEN
               WRITE(6,*) ' SB after dmk0one'
               call wrtmat(sb,nia,nib,nia,nib)
            END IF
*
CTF  First call for two electron routines ...
CTF  First call for two electron routines ...
CTF  First call for two electron routines ...
CTF  First call for two electron routines ...
*
            IF(IDOH2.NE.0.AND.NIBEL.GE.2.and.IELCI.ge.2) THEN
               if (NTEST.ge.100) then
                  write(6,*) '...........................'
                  write(6,*) 'I even have > 1 barred el.'
                  write(6,*) '...........................'
               end if
*     . Two electron part beta-beta
               IF(NTEST.GE.101) then
                  WRITE(6,*) ' I am going to call dmk0two'
                  write(6,*) ' (Work is actually done by RSBB2A)'
               end if
*
C?             WRITE(6,*) 'Do memory check ...before RSBB2A'
C?             call memchk
*
*     . barred excitations ( no they are not barred)
               IBU = 1
               IF(IRI.NE.IRICO) THEN
                  IRI_INT = 2
               ELSE
                  IRI_INT = 1
               END IF
               IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag*imag => real gives minus
                  SCLFAC2 = - SCLFAC
               ELSE
                  SCLFAC2 = SCLFAC
               END IF
               CALL RSBB2A_REL(IBSM,IBTP,JBSM,JBTP,IBGRP,NIA,NIB,
     &              NGAS,IBOC,JBOC,
     &              SB,CB,
     &              ADSXA,DXSTST,STSTDX,SXDXSX,
     &              NOBPTS2,IOBPTS2,ITSOB,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &              NSMOB,NSMST,NSMSX,NSMDX,
     &              CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC2,NTEST,
     &              NSEL2E,ISEL2E,IRI_INT,IBU,SOMESCR)
               IF(NTEST.GE.500) THEN
                  WRITE(6,*) ' SB after dmk0two'
                  call wrtmat(sb,nia,nib,nia,nib)
               END IF
*
            END IF
*           ^End if active 2 electron part
            CALL TRPMT3(SB,NIA,NIB,C2)
            CALL COPVEC(C2,SB,NIA*NIB)
            CALL TRPMT3(CB,NJA,NJB,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
         END IF
*        ^ End if IASM = JASM and IATP = JATP and atleast one electron
*
            IF(NIAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
               if (NTEST.ge.100) then
                  write(6,*) '.......................'
                  write(6,*) 'I have unbarred elecs.'
                  write(6,*) '.......................'
               end if
*
*     alpha single excitation
*
               IF(NTEST.GE.101)
     &            WRITE(6,*) 'I am going to call dmk0one (last time)'
*
               IUB = 1
               IPM1 = 0
C?             WRITE(6,*) 'Do memory check ...before dmk0one'
C?             call memchk
*
CSK - 6th of October 2009 - commented out because of wrong number of
Carguments in calling list - thanks to Miro for pointig it out
CSK            call dmk0one(IASM,IATP,JASM,JATP,
CSK  &              IAGRP,NIB,
CSK  &              NGAS,IAOC,JAOC,
CSK  &              SB,CB,
CSK  &              ADSXA,STSTSX,
CSK  &              NOBPTS,IOBPTS,
CSK  &              MAXI,MAXK,
CSK  &              SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
CSK  &              NSMOB,NSMST,MOCAA,
CSK  &              SCLFAC,
CSK  &              ONEINTS,INTSPACE,
CSK  &              IRI,IRICO,IPM1,IUB,ISMDST,NTEST)
*
               IF(NTEST.GE.500) THEN
                  WRITE(6,*) ' SB transposed  after dmk0one'
                  call wrtmat(SB,nib,nia,nib,nia)
               END IF
            end if
*           ^ End NIAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM
*
*     alpha double excitation
*
*
CTF  Second call for two electron routines ...
CTF  Second call for two electron routines ...
CTF  Second call for two electron routines ...
CTF  Second call for two electron routines ...
*
            if (IDOH2.ne.0.and.NIAEL.ge.2.and.IELCI.ge.2.AND.
     &          IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) then
               if (NTEST.ge.100) then
                WRITE(6,*) ' I am going to call RSBB2A (last time )'
                write(6,*) '..................................'
                write(6,*) 'meaning I have >1 unbarred elecs.'
                write(6,*) '..................................'
               end if
*
C?             WRITE(6,*) 'Do memory check ...before RSBB2A'
C?             call memchk
*
*     . unbarred excitations
               IBU = 2
               IF(IRI.NE.IRICO) THEN
                  IRI_INT = 2
               ELSE
                  IRI_INT = 1
               END IF
               IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
                  SCLFAC2 = - SCLFAC
               ELSE
                  SCLFAC2 = SCLFAC
               END IF
               CALL RSBB2A_REL(IASM,IATP,JASM,JATP,IAGRP,NIB,NIA,
     &              NGAS,IAOC,JAOC,
     &              SB,CB,
     &              ADSXA,DXSTST,STSTDX,SXDXSX,
     &              NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &              NSMOB,NSMST,NSMSX,NSMDX,
     &              CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC2,NTEST,
     &              NSEL2E,ISEL2E,IRI_INT,IBU,SOMESCR)
            end if
*           ^ End if 2-electron aa is active
*
*
            IF(NTEST.GE.500) THEN
               WRITE(6,*) ' SB after dmk0two'
               call wrtmat(sb,nia,nib,nia,nib)
            END IF
*
CE?      END IF
*        ^End of IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM
*
*     Alpha-beta excitation
*
CTF  Third call for two electron routines ...
CTF  Third call for two electron routines ...
CTF  Third call for two electron routines ...
CTF  Third call for two electron routines ...
*
         if (NIAEL.ge.1.and.NIBEL.ge.1.and.IELCI.ge.2) then
            if (NTEST.ge.100) then
               write(6,*) '.................................'
               write(6,*) 'I have barred AND unbarred elecs.'
               write(6,*) '.................................'
            end if
*
C?          WRITE(6,*) 'Do memory check ...before RSBB2B'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
            CALL RSBB2B_REL(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &           JASM,JATP,JBSM,JBTP,NJA,NJB,
     &           IAGRP,IBGRP,NGAS,
     &           IAOC,IBOC,JAOC,JBOC,
     &           SB,CB,
     &           ADSXA,STSTSX,
     &           NOBPTS,IOBPTS,NOBPTS2,IOBPTS2,
     &           ITSOB,MAXK,
     &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           XINT,NSMOB,NSMST,NSMSX,NSMDX,
     &           IUSEAB,ICJKAIB,CJRES,SIRES,C2,SCLFAC,NTEST,
     &           NSEL2E,ISEL2E,IRI_INT,SOMESCR)
CSK old          NSEL2E,ISEL2E,IRI_INT,SOMESCR,ISMDST)
         end if
      END IF
*     ^End of delta_mk2 = 0 part
*
      IF(NTEST.GE.500) THEN
         WRITE(6,*) ' SB after DMK2=0'
         call wrtmat(sb,nia,nib,nia,nib)
      END IF
*
      IF(DELTA_MK2 .EQ. 2.and.ISPINFREE.eq.0) THEN
         if (NTEST.ge.100) then
            write(6,*) '........................'
            write(6,*) 'Entering DMK2=2 part ...'
            write(6,*) '........................'
         end if
*
*     ===================================
*     Sigma (delta Mk = +1) contribution
*     ===================================
*
         IF(NTEST.GE.101)
     &        WRITE(6,*) ' I am going to call DMK1ONE'
*
         CALL TRPMT3(SB,NIB,NIA,C2)
         CALL COPVEC(C2,SB,NIA*NIB)
         CALL TRPMT3(CB,NJB,NJA,C2)
         CALL COPVEC(C2,CB,NJA*NJB)
*
         IPM1 = 1
         IUB = 0
         if (NTEST.ge.100) then
            write(6,*) '...................'
            write(6,*) 'dmk1one call with :'
            write(6,*) 'IPM1 =  ',IPM1
            write(6,*) 'IUB  =  ',IUB
            write(6,*) '...................'
         end if
*
C?       WRITE(6,*) 'Do memory check ...before dmk1one'
C?       call memchk
*
CSKold   call dmk1one(IASM,IBSM,IATP,IBTP,JASM,JBSM,
CSKold        JATP,JBTP,IAGRP,IBGRP,JAGRP,JBGRP,NIA,NIB,
CSKold        NJA,NJB,NGAS,IAOC,IBOC,JAOC,JBOC,
CSKold        SB,CB,
CSKold        ADSXA,STSTSX,NOBPTS,NOBPTS2,
CSKold        IOBPTS,IOBPTS2,MAXI,MAXK,
CSKold        SSCR,CSCR,I1,XI1S,I2,XI2S,
CSKold        XINT,NSMOB,NSMST,
CSKold        ONEINTS,INTSPACE,SCLFAC,
CSKold        IRI,IRICO,IPM1,IUB,ISMDST,NTEST)
*
C?       WRITE(6,*) 'Do memory check ...after dmk1one'
C?       call memchk
*
         if (NTEST.ge.500) then
            if (IPM1.eq.-1) write(6,*) 'MARKER-1'
            if (IPM1.eq.1) write(6,*) 'MARKER+1'
            write(6,*) 'Updated sigma after dmk1one'
            call wrtmat(SB,NIA,NIB,NIA,NIB)
         end if
*
*
CTF  Fourth call for two electron routines ...
CTF  Fourth call for two electron routines ...
CTF  Fourth call for two electron routines ...
CTF  Fourth call for two electron routines ...
*
         if (NIEL.ge.2.and.IELCI.ge.2) then
            if (NTEST.ge.100) then
               write(6,*) '.................'
               write(6,*) '>1 sigma electron'
               write(6,*) '.................'
            end if
*     . a+ ia a+ kb a lb a jb
            I_AC = 2
            JKL_AC(1) = 2
            JKL_AC(2) = 1
            JKL_AC(3) = 1
*
C?          WRITE(6,*) 'Do memory check ...before DMK1TWO'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
            IBU = 2
            JBU = 1
            KBU = 2
            LBU = 2
*
CTF
*  These routines are not tested yet!
*  In an initial relativistic CI run, all integrals
*  required in this part are = 0. This is because there
*  is no spin mixing of barred and unbarred spinors
*  (unbarred spinors are plain alpha, barred spinors
*   plain beta.) and therefore, spin integration makes
*  the integrals vanish.
CTF
*
c           CALL DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
c    &           NIB,IBSM,IBTP,IBGRP,IBOC,
c    &           NJA,JASM,JATP,JAGRP,JAOC,
c    &           NJB,JBSM,JBTP,JBGRP,JBOC,
c    &           I_AC,JKL_AC,NGAS,SB,CB,
c    &           NOBPTS,NOBPTS2,MAXK,
c    &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c    &           SCR,NSMOB,NSMST,IUSEAB,CJRES,SIRES,
c    &           SCLFAC2,NTEST,IUSE_PH,IPHGAS,IRI_INT,
c    &           IBU,JBU,KBU,LBU,
c    &           XIJILS)
*
C     DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
C     &                    NIB,IBSM,IBTP,IBGRP,IBOC,
C     &                    NJA,JASM,JATP,JAGRP,JAOC,
C     &                    NJB,JBSM,JBTP,JBGRP,JBOC,
C     &          I_AC,JKL_AC,NGAS,SB,CB,MXPNGAS,
C     &          NALPTS,NBEPTS,MAXK,
C     &          SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
C     &          XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
C     &          CJRES,SIRES,SCLFAC,NTEST,
C     &          IUSE_PH,IPHGAS,
C     &          XIJILS)
*
         end if
*     . Restore order
         CALL TRPMT3(SB,NIA,NIB,C2)
         CALL COPVEC(C2,SB,NIA*NIB)
         CALL TRPMT3(CB,NJA,NJB,C2)
         CALL COPVEC(C2,CB,NJA*NJB)
*
         if (NTEST.ge.500) then
            write(6,*) 'Updated sigma after dmk1two'
            call wrtmat(SB,NIB,NIA,NIB,NIA)
         end if
*
CTF  Fifth call for two electron routines ...
CTF  Fifth call for two electron routines ...
CTF  Fifth call for two electron routines ...
CTF  Fifth call for two electron routines ...
*
         if (NIEL.ge.2.and.IELCI.ge.2) then
*     . a ib a+  ka a+ la a ja term
            I_AC = 1
            JKL_AC(1) = 2
            JKL_AC(2) = 2
            JKL_AC(3) = 1
            IBU = 1
            JBU = 2
            KBU = 1
            LBU = 1
*     . Proceed with blocks in the form C(Jb,Ja) and S(Ib,Ia) ( as on entry  )
C?          WRITE(6,*) 'Do memory check ...before DMK1TWO'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
c           CALL DMK1TWO(NIB,IBSM,IBTP,IBGRP,IBOC,
c    &           NIA,IASM,IATP,IAGRP,IAOC,
c    &           NJB,JBSM,JBTP,JBGRP,JBOC,
c    &           NJA,JASM,JATP,JAGRP,JAOC,
c    &           I_AC,JKL_AC,NGAS,SB,CB,
c    &           NOBPTS2,NOBPTS,MAXK,
c    &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c    &           SCR,NSMOB,NSMST,IUSEAB,CJRES,SIRES,
c    &           SCLFAC2,NTEST,IUSE_PH,IPHGAS,IRI_INT,
c    &           IBU,JBU,KBU,LBU,XIJILS)
         end if
      END IF
*     ^End of DELTA_MK = 1 part
*
      IF(NTEST.GE.500) THEN
         WRITE(6,*) ' SB after DELTA_MK = +1'
         call wrtmat(sb,nia,nib,nia,nib)
      END IF
*
      IF(DELTA_MK2 .EQ. -2.and.ISPINFREE.eq.0) THEN
*
*     ===================================
*     Sigma (delta Mk = -1) contribution
*     ===================================
*
         IF(NTEST.GE.101)
     &        WRITE(6,*) ' I am going to call DMK1ONE'
*
         CALL TRPMT3(SB,NIB,NIA,C2)
         CALL COPVEC(C2,SB,NIA*NIB)
         CALL TRPMT3(CB,NJB,NJA,C2)
         CALL COPVEC(C2,CB,NJA*NJB)
*
         IPM1 =-1
         IUB = 0
         if (NTEST.ge.100) then
            write(6,*) '...................'
            write(6,*) 'dmk1one call with :'
            write(6,*) 'IPM1 =  ',IPM1
            write(6,*) 'IUB  =  ',IUB
            write(6,*) '...................'
         end if
*
C?       WRITE(6,*) 'Do memory check ...before dmk1one'
C?       call memchk
*
CSK - 6th of October 2009 - commented out because of wrong number of
Carguments in calling list - thanks to Miro for pointig it out
CSK      call dmk1one(IASM,IBSM,IATP,IBTP,JASM,JBSM,
CSK  &        JATP,JBTP,IAGRP,IBGRP,JAGRP,JBGRP,NIA,NIB,
CSK  &        NJA,NJB,NGAS,IAOC,IBOC,JAOC,JBOC,
CSK  &        SB,CB,
CSK  &        ADSXA,STSTSX,NOBPTS,NOBPTS2,
CSK  &        IOBPTS,IOBPTS2,MAXI,MAXK,
CSK  &        SSCR,CSCR,I1,XI1S,I2,XI2S,
CSK  &        XINT,NSMOB,NSMST,
CSK  &        ONEINTS,INTSPACE,SCLFAC,
CSK  &        IRI,IRICO,IPM1,IUB,ISMDST,NTEST)
*
         if (NTEST.ge.500) then
            if (IPM1.eq.-1) write(6,*) 'MARKER-1'
            if (IPM1.eq.1) write(6,*) 'MARKER+1'
            write(6,*) 'Updated sigma after dmk1one'
            call wrtmat(SB,NIA,NIB,NIA,NIB)
         end if
*
CTF  Sixth call for two electron routines ...
CTF  Sixth call for two electron routines ...
CTF  Sixth call for two electron routines ...
CTF  Sixth call for two electron routines ...
*
         if (NIEL.ge.2.and.IELCI.ge.2) then
*     . a  ia a+ jb a+ kb a lb
            IBU = 1
            JBU = 2
            KBU = 2
            LBU = 2
            I_AC = 1
            JKL_AC(1) = 2
            JKL_AC(2) = 2
            JKL_AC(3) = 1
C?          WRITE(6,*) 'Do memory check ...before DMK1TWO'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
c           CALL DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
c    &           NIB,IBSM,IBTP,IBGRP,IBOC,
c    &           NJA,JASM,JATP,JAGRP,JAOC,
c    &           NJB,JBSM,JBTP,JBGRP,JBOC,
c    &           I_AC,JKL_AC,NGAS,SB,CB,
c    &           NOBPTS,NOBPTS2,MAXK,
c    &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c    &           SCR,NSMOB,NSMST,IUSEAB,CJRES,SIRES,
c    &           SCLFAC2,NTEST,IUSE_PH,IPHGAS,IRI_INT,
c    &           IBU,JBU,KBU,LBU,XIJILS)
C     DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
C     &                    NIB,IBSM,IBTP,IBGRP,IBOC,
C     &                    NJA,JASM,JATP,JAGRP,JAOC,
C     &                    NJB,JBSM,JBTP,JBGRP,JBOC,
C     &          I_AC,JKL_AC,NGAS,SB,CB,MXPNGAS,
C     &          NALPTS,NBEPTS,MAXK,
C     &          SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
C     &          XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
C     &          CJRES,SIRES,SCLFAC,NTEST,
C     &          IUSE_PH,IPHGAS,XIJILS)
*
         end if
*     . Restore order
         CALL TRPMT3(SB,NIA,NIB,C2)
         CALL COPVEC(C2,SB,NIA*NIB)
         CALL TRPMT3(CB,NJA,NJB,C2)
         CALL COPVEC(C2,CB,NJA*NJB)
*
CTF  Seventh call for two electron routines ...
CTF  Seventh call for two electron routines ...
CTF  Seventh call for two electron routines ...
CTF  Seventh call for two electron routines ...
*
         if (NIEL.ge.2.and.IELCI.ge.2) then
*     . a+ ib a+  ja a ka a la term
            IBU  = 1
            JBU = 2
            KBU =  2
            LBU = 2
            I_AC = 2
            JKL_AC(1) = 2
            JKL_AC(2) = 1
            JKL_AC(3) = 1
*     . Proceed with blocks in the form C(Jb,Ja) and S(Ib,Ia) ( as on entry  )
C?          WRITE(6,*) 'Do memory check ...before DMK1TWO'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
c           CALL DMK1TWO(NIB,IBSM,IBTP,IBGRP,IBOC,
c    &           NIA,IASM,IATP,IAGRP,IAOC,
c    &           NJB,JBSM,JBTP,JBGRP,JBOC,
c    &           NJA,JASM,JATP,JAGRP,JAOC,
c    &           I_AC,JKL_AC,NGAS,SB,CB,
c    &           NOBPTS2,NOBPTS,MAXK,
c    &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c    &           SCR,NSMOB,NSMST,IUSEAB,CJRES,SIRES,
c    &           SCLFAC2,NTEST,IUSE_PH,IPHGAS,IRI_INT,
c    &           IBU,JBU,KBU,LBU,XIJILS)
         end if
      END IF
*     ^End of DELTA_MK = -1 part
*
      IF(NTEST.GE.500) THEN
         WRITE(6,*) ' SB after DELTA_MK = -1'
         call wrtmat(sb,nia,nib,nia,nib)
      END IF
*
CTF  Eighth call for two electron routines ...
CTF  Eighth call for two electron routines ...
CTF  Eighth call for two electron routines ...
CTF  Eighth call for two electron routines ...
*
      if (DELTA_MK2.eq.4.and.NIEL.ge.2.and.IELCI.ge.2.
     &                    and.ISPINFREE.eq.0) then
         if (NTEST.ge.100) then
            write(6,*) '........................'
            write(6,*) 'Entering DMK2=4 part ...'
            write(6,*) '........................'
         end if
*
*     ===================================
*     Sigma (delta Mk = +2) contribution
*     ===================================
*
*     +2 contribution:
C?       WRITE(6,*) 'Do memory check ...before dmk2two'
C?       call memchk
*
         IF(IRI.NE.IRICO) THEN
            IRI_INT = 2
         ELSE
            IRI_INT = 1
         END IF
         IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
            SCLFAC2 = - SCLFAC
         ELSE
            SCLFAC2 = SCLFAC
         END IF
CTF
         IIDMK2 = 2
*  See what's going on here in detail when routines are
*  tested ...
CTF
c        call dmk2two(IBSM,IBTP,NIB,IOCTPB,
c    &        IASM,IATP,NIA,IOCTPA,
c    &        JBSM,JBTP,NJB,JOCTPB,
c    &        JASM,JATP,NJA,JOCTPA,
c    &        NGAS,IBOC,IAOC,JBOC,JAOC,
c    &        SB,CB,
c    &        ADSXA,
c    &        NOBPTS,NOBPTS2,MAXI,MAXK,
c    &        SSCR,CSCR,I1,XI1S,I2,XI2S,
c    &        XINT,
c    &        NSMOB,NSMST,
c    &        SCLFAC2,NTEST,IRI_INT,IIDMK2,
c    &        SOMESCR)
      end if
*
*
*
CTF  Ninth call for two electron routines ...
CTF  Ninth call for two electron routines ...
CTF  Ninth call for two electron routines ...
CTF  Ninth call for two electron routines ...
*
      if (DELTA_MK2.eq.-4.and.NIEL.ge.2.and.IELCI.ge.2.
     &                    and.ISPINFREE.eq.0) then
         if (NTEST.ge.100) then
            write(6,*) '.........................'
            write(6,*) 'Entering DMK2=-4 part ...'
            write(6,*) '.........................'
         end if
*
*     ===================================
*     Sigma (delta Mk = -2) contribution
*     ===================================
*
C?       WRITE(6,*) 'Do memory check ...before dmk2two'
C?       call memchk
*
         IF(IRI.NE.IRICO) THEN
            IRI_INT = 2
         ELSE
            IRI_INT = 1
         END IF
         IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
            SCLFAC2 = - SCLFAC
         ELSE
            SCLFAC2 = SCLFAC
         END IF
         IIDMK2 = -2
c        call dmk2two(IASM,IATP,NIA,IOCTPA,
c    &        IBSM,IBTP,NIB,IOCTPB,
c    &        JASM,JATP,NJA,JOCTPA,
c    &        JBSM,JBTP,NJB,JOCTPB,
c    &        NGAS,IAOC,IBOC,JAOC,JBOC,
c    &        SB,CB,
c    &        ADSXA,
c    &        NOBPTS2,NOBPTS,MAXI,MAXK,
c    &        SSCR,CSCR,I1,XI1S,I2,XI2S,
c    &        XINT,
c    &        NSMOB,NSMST,
c    &        SCLFAC2,NTEST,IRI_INT,IIDMK2,
c    &        SOMESCR)
      end if
*
*
*
      IF(NTEST.GE.50) THEN
         WRITE(6,*) ' =========================='
         WRITE(6,*) ' RSSBCBREL : Final S block '
         WRITE(6,*) ' =========================='
         CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
      END IF
*
 1001 CONTINUE
      call qexit('RSSBCBREL')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SBLOCK_REL(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
     &                      LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,
     &                      IRIS,IRIC,I_DO_TRNSP)
*
* Generate a set of sigma blocks,
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order
* If LUCBLK.GT.0, the blocks of C corresponding to IBLOCK
* are stored on LUCBLK
*
* ICONSPA,ICONSPB  added October 1996
* ICBAT_RES, ICBAT_INI, IBBAT_END added august 1997
*
* If ICBAT_RES .eq.1 then it as assumed that only
* Cbatches ICBAT_INI to ICBAT_END are stored on  LUC
*
* Generalized to NOPAIR Hamiltonian, T.F + J.O., March21, 7.22
*                                    (not bad for a saturday)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
* IA/BS/CTP : Type (number of elecs in A/B,S/C)
#include "cands.inc"
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
*
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
*. Used : NSMOB
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
#include "cintfo.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      COMMON/HIDSCRR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*
CTF
* Scratch space for rsbb2a, rsbb2bn, rsbb2bvn, skickj
! miro: very large in size according to readelf printout
CSK   dimension SOMESCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
CSK   commented out: this LUCIAREL route is not active anymore!
      dimension XIJILS(MXTSOB)
*
      IDUM = 0
      KDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'SBLOCK')
*
      NTEST = 0
      IF(NTEST.GE.5)
     &WRITE(6,*) ' SBLOCK : ISSPC,ICSPC ', ISSPC,ICSPC
C?    WRITE(6,*) ' LUC in SBLOCK ', LUC
      IF(LUCBLK.GT.0) THEN
        CALL REWINE(LUCBLK,-1)
      END IF
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = IASTP
      IBTP = IBSTP
*
      JATP = IACTP
      JBTP = IBCTP
*. Number of supergroups
      NOCTPA_S = NOCTYP(IATP)
      NOCTPB_S = NOCTYP(IBTP)
*
      NOCTPA_C = NOCTYP(JATP)
      NOCTPB_C = NOCTYP(JBTP)
*. Offset for supergroups
      IOCTPA_S = IBSPGPFTP(IATP)
      IOCTPB_S = IBSPGPFTP(IBTP)
*
      IOCTPA_C = IBSPGPFTP(JATP)
      IOCTPB_C = IBSPGPFTP(JBTP)
*
      NAEL_S = NELEC(IATP)
      NBEL_S = NELEC(IBTP)
*
      NAEL_C = NELEC(JATP)
      NBEL_C = NELEC(JBTP)
*
*. connection matrices for supergroups, just st for for s-space
*  (only used when number of elecs are identical )
*
      call memmar(KCONSPA_S,NOCTPA_S**2,'ADDL  ',1,'CNSPAC')
      call memmar(KCONSPB_S,NOCTPB_S**2,'ADDL  ',1,'CNSPBC')
C     SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
*
      CALL SPGRPCON(IOCTPA_S,NOCTPA_S,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA_S),IPRCIX)
      CALL SPGRPCON(IOCTPB_S,NOCTPB_S,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB_S),IPRCIX)
*
*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      call memmar(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      call memmar(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM_REL(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MXSTBL = MXNSTR
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
      MAXIK = MAX(MAXI,MAXK)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
      MAXIJ = MXTSOB_L ** 2
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MXTSOB_L ** 4
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      call memmar(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*.Local scratch arrays for blocks of C and sigma
      IF(ICISTR.LE.2) THEN
        L0BLOCK = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        L0BLOCK = MXSOOB
      END IF
      LSCR1 = MAX(L0BLOCK,L2BLOCK)
Cold      LSCR1 = MAX(LSCR1,LCSBLK)
      if (IPRCIX.ge.2) then
        WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
      end if
*. Arrays giving allowed type combinations
      call memmar(KCIOIO,NOCTPA_C*NOCTPB_C,'ADDL  ',2,'CIOIO ')
      call memmar(KSIOIO,NOCTPA_S*NOCTPB_S,'ADDL  ',2,'SIOIO ')
*. sigma needed for MXRESC
C          IAIBCM_REL(ICISPC,IAIB)
      CALL IAIBCM_REL(ISSPC,IATP,IBTP,WORK(KSIOIO))
      CALL IAIBCM_REL(ICSPC,JATP,JBTP,WORK(KCIOIO))
*. Arrays giving block type
      call memmar(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. for SIgmav reflection - not active
CTF   KSVST = 1
*  Why the heck is WORK(1) used here??
*  This is damn dangerous. ISTRFL allocated statically in sblocks
*
*.scratch space for projected matrices and a CI block
      LSCR2 = LCJRES
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
C  memory for blocks already allocated, so
      KC2 = KVEC3
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*. Allocation of memory for annihilation/creation from supergroups
CTF setting this undefined variable to two:
      MXSXST = 2
      MXADKBLK = 2
*. The size of XI* should be obtained from MXRESCP, but for now ..
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0,
     &             MXSTBL*MXTSOB_L)
      call memmar(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      call memmar(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
*
      call memmar(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      call memmar(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
*
      call memmar(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      call memmar(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
*
      call memmar(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      call memmar(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
*. Type of blocks (all types set to 1 for now )
      CALL ZBLTP_REL(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP))
*.Some TTS arrays
      NTTS = MXNTTS
*. for partitioning of vector
      call memmar(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      call memmar(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      call memmar(KLI1BT,2*NTTS  ,'ADDL  ',1,'I1BTC ')
      call memmar(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*. For scaling for each TTS block
      call memmar(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')
*. Space for four blocks of string occupations and arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL_S,NAEL_C,NBEL_S,NBEL_C)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL_S,NBEL_S,NAEL_C,NBEL_C)) * NOCOB
      DO I1234 = 1, 2
      call memmar(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      call memmar(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
      call memmar(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      call memmar(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*
* Allocate memory for one array of symmetrydistributions
      CALL MEMMAR(KLISMDST,MAX_NSYMDIST,'ADDL  ',1,'ISMDST')
*
* Fed in directly from common block oper.inc.
* Set in readin !!
      IF(IHAM12.EQ.2) THEN
* One- and two-electron Hamiltonian active
        IDOH2 = 1
      ELSE
* Only one-electron Hamiltonian active
        IDOH2 = 0
      END IF
*
      CALL SBLOCKS_REL(WORK(KPICI1),
     &            NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KC2),
     &            WORK(KCIOIO),ISMOST(1,ICSM),WORK(KCBLTP),
     &            NACOB,WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &            WORK(KNSTSO(JATP)),WORK(KNSTSO2(JBTP)),
     &            NAEL_S,IATP,NBEL_S,IBTP,
     &            NAEL_C,JATP,NBEL_C,JBTP,
     &            IOCTPA_S,IOCTPB_S,NOCTPA_S,NOCTPB_S,
     &            IOCTPA_C,IOCTPB_C,NOCTPA_C,NOCTPB_C,
     &            NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,
     &            NOBPTS2,IOBPTS2,MXPNGAS,
     &            ITSOB,MAXIJ,MAXK,MAXI,INSCR,LSCR1,
     &            LSCR1,WORK(KINSCR),WORK(KCSCR),WORK(KSSCR),
     &            SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &            ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &            MXTSOB,SOMESCR,XIJILS,
     &            WORK(KI1),WORK(KXI1S),
     &            WORK(KI2),WORK(KXI2S),IDOH2,MXPOBS,
     &            PSSIGN,IPRSIG,LUC,ICJKAIB,WORK(KCJRES),
     &            WORK(KSIRES),WORK(KI3),WORK(KXI3S),
     &            WORK(KI4),WORK(KXI4S),MXSXST,MXSXBL,
     &            MOCAA,MOCAB,IAPR,
     &            WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &            IRESTRICT,WORK(KCONSPA_S),WORK(KCONSPB_S),
     &            WORK(KLSCLFAC),
     &            LUCBLK,IPERTOP,IH0INSPC,ISPINFREE,
CTI  &            WORK(KLH0SPC),
     &            ICBAT_RES,ICBAT_INI,ICBAT_END,
     &            IRELACI,IRIS,IRIC,I_DO_TRNSP,IELCI,
     &            WORK(KLISMDST) )
*
C?    write(6,*) 'Checking IDC  ',IDC
      IF(IDC.EQ.2) THEN
*. reform
        CALL RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,
     &             1,NSMST,NOCTPA_S,NOCTPB_S,
     &             WORK(KNSTSO(IATP)), WORK(KNSTSO2(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK,NSMST,NOCTPA_S,NOCTPB_S,
     &              WORK(KNSTSO(IATP_S)), WORK(KNSTSO2(IBTP_S)),
     &              IDC,1,NTEST)
      END IF
*
      IF(LUCBLK.GT.0) THEN
        CALL ITODS(-1,1,-1,LUCBLK)
      END IF
*. Eliminate local memory
      IDUM = 0
      call memmar(KDUM ,IDUM,'FLUSM ',2,'SBLOCK')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SBLOCKS_REL(ONEINTS,
     &                NSBLOCK,ISBLOCK,CB,SB,C2,
     &                ICOCOC,ICSMOS,ICBLTP,NACOB,NSSOA_S,NSSOB_S,
     &                NSSOA_C,NSSOB_C,
     &                NAEL_S,IAGRP_S,NBEL_S,IBGRP_S,
     &                NAEL_C,IAGRP_C,NBEL_C,IBGRP_C,
     &                IOCTPA_S,IOCTPB_S,NOCTPA_S,NOCTPB_S,
     &                IOCTPA_C,IOCTPB_C,NOCTPA_C,NOCTPB_C,
     &                NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,
     &                NOBPTS2,IOBPTS2,MXPNGAS,
     &                ITSOB,MAXIJ,MAXK,MAXI,LI,LC,LS,
     &                XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &                SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &                MXTSOB,SOMESCR,XIJILS,
     &                I1,XI1S,I2,XI2S,IDOH2,MXPOBS,
     &                PS,IPRNT,LUC,ICJKAIB,CJRES,SIRES,I3,XI3S,
     &                I4,XI4S,MXSXST,MXSXBL,MOCAA,MOCAB,IAPR,
     &                LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,IRESTRICT,
     &                ICONSPA,ICONSPB,SCLFAC,
     &                LUCBLK,IPERTOP,IH0INSPC,ISPINFREE,
CTF  &                IH0SPC,
     &                ICBAT_RES,ICBAT_INI,ICBAT_END,
     &                IRELACI,IRIS,IRIC,I_DO_TRNSP,IELCI,
     &                ISMDST)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen , Winter of 1991
*
* Modified for NO-PAIR DBG CI , T.F + J.O, March 21 98
*
* ISMDST added, Aug. 2003
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included
*   ISBLOCK(1,*) : alpha type of block
*   ISBLOCK(2,*) : beta type of block
*   ISBLOCK(3,*) : sym of alpha in block
*   ISBLOCK(4,*) : Offset of block
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* Note : Prescreen + Restrict features eliminated
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Specific input
      INTEGER ISBLOCK(8,*)
#include "picasso.inc"
      dimension ONEINTS(INTSPACE)
*.General input
      INTEGER ICOCOC(NOCTPA_C,NOCTPB_C)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA_S(NSMST,*),NSSOB_S(NSMST,*)
      INTEGER NSSOA_C(NSMST,*),NSSOB_C(NSMST,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST),ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NOBPTS2(MXPNGAS,*),IOBPTS2(MXPNGAS,*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA_S,NOCTPA_S),ICONSPB(NOCTPB_S,NOCTPB_S)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),I3(*),I4(*),
     &          XI1S(*),XI2S(*),XI3S(*),XI4S(*)
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      dimension ISTRFL(12)
      DIMENSION SCLFAC(*)
      INTEGER ISMDST(*)
*
      REAL * 8 INPROD
*
      CALL QENTER('SBLOCK_S')
C     WRITE(6,*) ' RMEMCHK at start of SBLOCKS'
C     CALL RMEMCHK
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' SBLOCKS speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(6,*) ' TTSS for each ACTIVE sigma block'
          DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0)
     &      WRITE(6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
          END DO
          WRITE(6,*) ' IDC PS ', IDC,PS
          WRITE(6,*) ' IDOH2 = ',IDOH2
*
          WRITE(6,*) ' NAEL_C, NBEL_C, NAEL_S,NBEL_S',
     &                 NAEL_C, NBEL_C, NAEL_S,NBEL_S
      END IF
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Initial C vector '
        CALL WRTVCD(CB,LUC,1,-1)
      END IF
      NTEST = 00
* ===========================
* 1 : Arrays for accessing C
* ============================
*. Find batches of C - strings
      ITTSS_ORD = 1
      CALL PART_CIV2(IDC,ICBLTP,NSSOA_C,NSSOB_C,NOCTPA_C,NOCTPB_C,NSMST,
     &               LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &               I1CBLOCK,ICBLOCK,0,ITTSS_ORD)
*. Find the active blocks on LUC, store info in SCLFAC
C     FIND_ACTIVE_BLOCKS(LUIN,LBLK,BLK_A,SEGMNT)
*     \/ Change to accomodate for ICBAT_RES
      CALL FIND_ACTIVE_BLOCKS(LUC,-1,SCLFAC,CB)
C     write(6,*) ' memchk after FIND_ACTIVE_BLOCKS'
C     CALL RMEMCHK
*
* Loop over batches over C blocks
      IF(IDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE
        MXEXC = 1
      END IF
C?    WRITE(6,*) ' IDOH2, MXEXC',IDOH2,MXEXC
      REWIND LUC
      IF(ICBAT_RES.EQ.1) THEN
        WRITE(6,*) ' Restricted set of C batches '
        WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
        JCBAT_INI = ICBAT_INI
        JCBAT_END = ICBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = NCBATCH
      END IF
*
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END
        if (NTEST.ge.1000) then
          write(6,*) 'Batch info from SBLOCKS:'
          write(6,*) 'JCBATCH = ',JCBATCH
          write(6,*) 'JCBAT_END = ',JCBAT_END
        end if
*
*. Read C blocks of this batch into core
*
        ICOFF = 1
        NJBLOCK = LCBLOCK(JCBATCH)
        DO JJCBLOCK = 1, NJBLOCK
          if (NTEST.ge.1000) then
            write(6,*) 'Block info from SBLOCKS:'
            write(6,*) 'JJCBLOCK = ',JJCBLOCK
            write(6,*) 'NJBLOCK = ',NJBLOCK
          end if
*
          JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
          JATP = ICBLOCK(1,JBLOCK)
          JBTP = ICBLOCK(2,JBLOCK)
          JASM = ICBLOCK(3,JBLOCK)
          JBSM = ICBLOCK(4,JBLOCK)
          JOFF = ICBLOCK(5,JBLOCK)
          ISCALE = 1
*!!!
CTF
*  Caution: C and ICOOSC are not initialized at the moment.
*           But we read CB in from disc, so they are currently
*           NOT used in GSTTBL. If LUC is zero, they will be
*           used, i.e. C is input array for CB !
CTF
          CALL GSTTBL(C,CB(JOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                NOCTPA_C,NOCTPB_C,NSSOA_C,NSSOB_C,PS,
     &                ICOOSC,IDC,
     &                PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK))
*
          if (NTEST.ge.1000) then
            NJATEST = NSSOA_C(JASM,JATP)
            NJBTEST = NSSOB_C(JBSM,JBTP)
            write(6,*) 'after GSTTBL in sblocks'
            write(6,*) 'JOFF is  ',JOFF
            call wrtmat(CB(JOFF),NJBTEST,NJATEST,
     &                  NJBTEST,NJATEST)
          end if
*
          SCLFAC(JBLOCK) = 1.0d0
        END DO
*
*. Loop over blocks of sigma and C in core and obtain  contribution from
*  given C block to given S block
*. Loop over C blocks
        DO 9000 ICBLK = I1CBLOCK(JCBATCH),
     &                I1CBLOCK(JCBATCH)-1+NJBLOCK
          if (NTEST.ge.1000) then
            write(6,*) 'Block info from SBLOCKS:'
            write(6,*) 'ICBLK = ',ICBLK
            write(6,*) 'I1CBLOCK(JCBATCH)-1+NJBLOCK = ',
     &                  I1CBLOCK(JCBATCH)-1+NJBLOCK
          end if
*
          JATP = ICBLOCK(1,ICBLK)
          JBTP = ICBLOCK(2,ICBLK)
          JASM = ICBLOCK(3,ICBLK)
          JBSM = ICBLOCK(4,ICBLK)
          ICOFF = ICBLOCK(5,ICBLK)
          NJA = NSSOA_C(JASM,JATP)
          NJB = NSSOB_C(JBSM,JBTP)
*
          IF(SCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
            CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                  LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
*. Start with transposed block
            DO 8765 IPERM = NPERM,1,-1
              if (NTEST.ge.1000) then
                write(6,*) 'Permutation info from SBLOCKS:'
                write(6,*) 'IPERM = ',IPERM
                write(6,*) 'NPERM = ',NPERM
              end if
*
              LLASM = LASM(IPERM)
              LLBSM = LBSM(IPERM)
              LLATP = LATP(IPERM)
              LLBTP = LBTP(IPERM)
              NLLA = NSSOA_C(LLASM,LLATP)
              NLLB = NSSOB_C(LLBSM,LLBTP)
*. The routines assumes on input that the blocks are transposed, so,
*. Initial tour, IPERM = 1 always corresponds to no transpose, so transpose!
              IF(IPERM.EQ.1) THEN
                IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM)
     &                   THEN
*. Diagonal blocks, Transposing corresponds to scaling
                  IF(PS.EQ.-1.0D0) THEN
                    CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
                  END IF
                ELSE
*. offdiagonal blocks, explicit transposing
                  CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
                  CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
                  if (LSGN(IPERM).EQ.-1) then
                    write(6,*) 'LCOL and LROW undefined in sblocks.'
                    Call Abend2( 'Quitting.' )
CTF                 CALL SCALVE(CB(ICOFF),-1.0D0,LROW*LCOL)
                  end if
                END IF
              END IF
*
              DO 10000 ISBLK = 1, NSBLOCK
                if (NTEST.ge.1000) then
                  write(6,*) 'Block info from SBLOCKS:'
                  write(6,*) 'ISBLK = ',ISBLK
                  write(6,*) 'NSBLOCK = ',NSBLOCK
                end if
*
                IF(ISBLOCK(1,ISBLK) .GT. 0 ) THEN
                  IATP = ISBLOCK(1,ISBLK)
                  IBTP = ISBLOCK(2,ISBLK)
                  IASM = ISBLOCK(3,ISBLK)
                  IBSM = ISBLOCK(4,ISBLK)
                  ISOFF = ISBLOCK(5,ISBLK)
                  NIA = NSSOA_S(IASM,IATP)
                  NIB = NSSOB_S(IBSM,IBTP)
*
                  IF(NIA*NIB.EQ.0) GOTO 10000
                  IF(NJA*NJB.EQ.0) GOTO 10000
*
                  CNORM = INPROD(CB(ICOFF),CB(ICOFF),NJA*NJB)
                  IF(CNORM.LT.1.0D-15) GOTO 10000
*
                  IF(NTEST.GE.100) THEN
                    WRITE(6,*) ' Next s block in batch : '
                    write(6,*) ' ISBLK IASM IBSM IATP IBTP'
                    write(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  END IF
*
                  IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &              ((LLBSM.GT.LLASM).OR.
     &              (LLASM.EQ.LLBSM).AND.(LLBTP.GT.LLATP)))
     &              GOTO 8764
*
                  IF(NTEST.GE.60) THEN
                    WRITE(6,*) ' RSSBCB will be called for '
                    WRITE(6,*) ' Sigma block : '
                    WRITE(6,*) ' ISOFF ', ISOFF
                    WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
                    WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                    WRITE(6,*) ' C     block : '
                    WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                    WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                    WRITE(6,*) ' ICOFF ', ICOFF
                    WRITE(6,*) ' Overall scale',SCLFAC(ICBLK)
                  END IF
*
                  XFAC = SCLFAC(ICBLK)
*
                  if (IRELACI.eq.0) then
*. Nonrelativistic Hamiltonian
C                   CALL RSSBCB2(IASM,IATP,IOCTPA_S,
C    &                 IBSM,IBTP,IOCTPB_S,
C    &                 LLASM,LLATP,LLBSM,LLBTP,NGAS,
C    &                 NELFSPGP(1,IATP+IOCTPA_S-1),
C    &                 NELFSPGP(1,IBTP+IOCTPB_S-1),
C    &                 NELFSPGP(1,LLATP+IOCTPA_S-1),
C    &                 NELFSPGP(1,LLBTP+IOCTPB_S-1),
C    &                 NAEL_S,NBEL_S,
C    &                 IAGRP_S,IBGRP_S,
C    &                 SB(ISOFF),CB(ICOFF),IDOH2,
C    &                 ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
C    &                 NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
C    &                 SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
C    &                 NSMOB,NSMST,NSMSX,NSMDX,
C    &                 NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
C    &                 ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
C    &                 MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,
C    &                 IPRNT,IPERTOP,IPTSPC,JPTSPC,IJOP,0,IDUM,XFAC,
C    &                 MXTSOB,SOMESCR,XIJILS)
                  else if (IRELACI.eq.1) then
*. No-pair relativistic Hamiltonian
                    IATPABS = IATP +IOCTPA_S-1
                    IBTPABS = IBTP +IOCTPB_S-1
                    JATPABS = JATP +IOCTPA_C-1
                    JBTPABS = JBTP +IOCTPB_C-1
*
                    NTEST = 200
                    if (NTEST.ge.200) then
                      write(6,*) 'String type offsets in sblocks:'
                      write(6,*) 'IATPABS = ',IATPABS
                      write(6,*) 'IBTPABS = ',IBTPABS
                      write(6,*) 'JATPABS = ',JATPABS
                      write(6,*) 'JBTPABS = ',JBTPABS
                      write(6,*) 'IATP,IOCTPA_S',IATP,IOCTPA_S
                      write(6,*) 'IBTP,IOCTPB_S',IBTP,IOCTPB_S
                      write(6,*) 'JATP,IOCTPA_C',JATP,IOCTPA_C
                      write(6,*) 'JBTP,IOCTPB_C',JBTP,IOCTPB_C
                    end if
*
                    call rssbcbrel(ONEINTS,INTSPACE,
     &         IASM,IATP,IOCTPA_S,NELFSPGP(1,IATPABS),IAGRP_S,NIA,
     &         IBSM,IBTP,IOCTPB_S,NELFSPGP(1,IBTPABS),IBGRP_S,NIB,
     &         JASM,JATP,IOCTPA_C,NELFSPGP(1,JATPABS),IAGRP_C,NJA,
     &         JBSM,JBTP,IOCTPB_C,NELFSPGP(1,JBTPABS),IBGRP_C,NJB,
     &                  NGAS,SB(ISOFF),CB(ICOFF),IDOH2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  NOBPTS,IOBPTS,NOBPTS2,IOBPTS2,
     &                  MXPNGAS,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &                  NSMOB,NSMST,NSMSX,NSMDX,
     &                  MXPOBS,IDC,
     &                  CJRES,SIRES,I3,XI3S,I4,XI4S,
     &                  MXSXBL,MXSXST,IPRNT,ISPINFREE,
     &                  IPTSPC,JPTSPC,IJOP,SCLFAC,
     &                  IRIS,IRIC,IELCI,
     &                  MXTSOB,SOMESCR,XIJILS,ISMDST)
C                   WRITE(6,*) 'before CB(19) = ', CB(19)
                  else
                    write(6,*) 'Something is wrong in SBLOCKS'
                    write(6,*) 'RELACI corrupted. IRELACI = ',IRELACI
                    Call Abend()
                  end if
*                 ^ End of switch relativistic/non-relativistic
 8764             CONTINUE        ! if great if block is true
                END IF
*               ^ end if ISBLOCK > 0
10000         CONTINUE
*.            ^ End of loop over sigma blocks
 8765       CONTINUE
*           ^ end loop start with transpose block
          END IF
*         ^ end  if SCLFAC > 0
 9000   CONTINUE
*.      ^ End of loop over C blocks in Batch
20000 CONTINUE
*.    ^End of loop over batches of C blocks
C     write(6,*) ' RMEMCHK after loop 20000'
C     CALL RMEMCHK
*
* Order
* The condition I_DO_TRNSP.EQ.1 is fulfilled if
*  IRIC.EQ.2.AND.IMK2_C.EQ.NMS2VAL
*  in rassg3. This means that the imaginary part (IRIC) of the
*  last MK2 value of the ket vectors (IMK2_C) for a given IRIS
*  and MK2_S is now being treated. After the sigma vector for
*  this partition is created, the final S block is reordered to
*  untransposed form (in the following).
*
      IF(I_DO_TRNSP.EQ.1) THEN
        DO  ISBLK = 1 , NSBLOCK
        IF(ISBLOCK(1,ISBLK).GT.0) THEN
          IATP = ISBLOCK(1,ISBLK)
          IBTP = ISBLOCK(2,ISBLK)
          IASM = ISBLOCK(3,ISBLK)
          IBSM = ISBLOCK(4,ISBLK)
          ISOFF  = ISBLOCK(5,ISBLK)
          ISOFFP = ISBLOCK(6,ISBLK)
          NIA = NSSOA_S(IASM,IATP)
          NIB = NSSOB_S(IBSM,IBTP)
*. Tranpose sigma block was obtained, transpose to obtain correct block
          CALL TRPMT3(SB(ISOFF),NIB,NIA,C2)
C?        WRITE(6,*) 'Transpose loop, NIA,NIB,ISOFF',NIA,NIB,ISOFF
          CALL COPVEC(C2,SB(ISOFF),NIA*NIB)
          IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP)
     &    CALL TRPAD3(SB(ISOFF),PS,NIA)
        END IF
        END DO
      END IF
*
      NTEST = 50
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' output blocks from SBLOCKS '
        CALL WRTTTS(SB,ISBLOCK,NSBLOCK,
     &              NSMST,NOCTPA_S,NOCTPB_S,NSSOA_S,NSSOB_S,1)
      END IF
      NTEST = 0
*
C?    WRITE(6,*) ' RMEMCHK at END of SBLOCKS'
C?    CALL RMEMCHK
      CALL QEXIT('SBLOCK_S')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SCTADRMT(ISCT,SSCT,AIN,AOUT,NROWIN,NROWOUT,NCOL,FACTOR)
*
* AOUT(ISCT(I),J) = AOUT(ISCT(I),J)+SSCT(I)*AIN(I,J)*FACTOR if ISCT(I) .NE.0
* No update                                 if ISCT(I) .EQ.0
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*.Input
      DIMENSION ISCT(*),SSCT(*),AIN(NROWIN,NCOL)
*.Output
      DIMENSION AOUT(NROWOUT,NCOL)
*
*. Currently programmed as stright loops,
* could be improved by processing batches of columns simultaneously
      DO J = 1, NCOL
        DO I =1, NROWIN
          IF(ISCT(I).NE.0) THEN
           AOUT(ISCT(I),J) = AOUT(ISCT(I),J) + SSCT(I)*FACTOR*AIN(I,J)
          END IF
        END DO
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SKICKJ_MK1(SKII,CKJJ,NKA,NIB,NJB,NKB,NMB,XIJKL,
     &                  NI,NJ,NK,NL,MAXI,MAXK,MAXM,
     &                  KBIB,XKBIB,MBKB,XMBKB,JBMB,XJBMB,IKORD,
     &                  FACS,MXTSOB,XIJILS)
*
*
* Calculate S(Ka,Ib) = FACS*S(Ka,Ib)
*          +SUM(i,j,k,l,Kb) <Ib!a  jb!Kb><Kb!a kb !Mb>*<Mb! a lb!Jb>
*           (ij!kl)*C(Ka,i,Jb)
*
* as
*
*          +SUM(i,j,k,l,Kb) <Kb!a  jb!Ib><Mb!a kb !Kb>*<Ib! a lb!Mb>
*           (ij!kl)*C(Ka,i,Jb)
*
* Jeppe Olsen, Winter of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION CKJJ(NKA,*)
      DIMENSION XIJKL(*)
*
      DIMENSION KBIB(MAXI,*),XKBIB(MAXI,*)
      DIMENSION MBKB(MAXK,*),XMBKB(MAXK,*)
      DIMENSION JBMB(MAXM,*),XJBMB(MAXM,*)
*
*. Input and output
      DIMENSION SKII(NKA,*)
*. Scratch
      DIMENSION XIJILS(MXTSOB)
      IF(NI.GT.MXTSOB.OR.NJ.GT.MXTSOB.OR.NK.GT.MXTSOB
     &   .OR.NL.GT.MXTSOB) THEN
         WRITE(6,*) ' SKICKJ : Too many orbs : NI > MXTSOB '
         WRITE(6,*) ' NI, MXTSOB ',MAX(NI,NJ,NK,NL),MXTSOB
         Call Abend2( ' Redim MXTSOB in SKICKJ' )
      END IF
*
      DO IB = 1, NIB
        DO J = 1, NJ
        IF(KBIB(IB,J).NE.0) THEN
          SIGNJ = XKBIB(IB,J)
          KB = KBIB(IB,J)
          DO K = 1, NK
          IF(MBKB(KB,K).NE.0) THEN
            MB = MBKB(KB,K)
            SIGNK = XMBKB(KB,K)
            DO L = 1, NL
            IF(JBMB(MB,L).NE.0) THEN
              JB = JBMB(MB,L)
              SIGNL = XJBMB(JB,L)
*. All connections have been made and
*. We have now a IB and Jb string, let's do it
C             ISOFF = (IB-1)*NKA + 1
C             ICOFF = (JB-1)*NI*NKA + 1
              ISOFF = IB
*. C(1,1,JB)
              ICOFF = (JB-1)*NI+ 1
*. Integrals are assumed stored as (ij!kl), offset to (1j!kl)
              INTOFF = ((L-1)*NK + K - 1 )*NI*NJ
     &               +  (J-1)*NI + 1
              CALL MATML7(SKII(1,ISOFF),CKJJ(1,ICOFF),XIJKL(INTOFF),
     &                    NKA,1,NKA,NI,NI,1,0,1.0D0,1.0D0)
            END IF
            END DO
          END IF
          END DO
        END IF
        END DO
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SKICKJ_REL(SKII,CKJJ,NKA,NIB,NJB,NKB,XIJKL,
     &     NI,NJ,NK,NL,MAXK,
     &     KBIB,XKBIB,KBJB,XKBJB,IKORD,
     &     FACS,IROUTE,MXTSOB)
*
*
* Calculate S(Ka,Ib,i) = FACS*S(Ka,Ib,i)
*          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
*
* Jeppe Olsen, Spring of 94
*
* : Note : Route 1 has retired, March 97
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION CKJJ(*)
      DIMENSION XIJKL(*)
      DIMENSION KBIB(MAXK,*),XKBIB(MAXK,*)
      DIMENSION KBJB(MAXK,*),XKBJB(MAXK,*)
*. Input and output
      DIMENSION SKII(*)
*. Scratch
C     COMMON/SOMESCR/SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
C     DIMENSION IBOFF(MXTSOB*MXTSOB),JBOFF(MXTSOB*MXTSOB)
C     DIMENSION SCR(NI*NL*NJ),IBOFF(NI),JBOFF(NL*NJ)
      DIMENSION XIJILS(MXTSOB)
      IF(NI.GT.MXTSOB.OR.NJ.GT.MXTSOB.OR.NK.GT.MXTSOB
     &     .OR.NL.GT.MXTSOB) THEN
         WRITE(6,*) ' SKICKJ : Too many orbs : NI > MXTSOB '
         WRITE(6,*) ' NI, MXTSOB ',MAX(NI,NJ,NK,NL),MXTSOB
         Call Abend2( ' Redim MXTSOB in SKICKJ' )
      END IF
*
C     CALL QENTER('SKICK')
      IF(IROUTE.EQ.3) THEN
* S(Ka,i,Ib) = S(Ka,i,Ib) + sum(j) (ji!kl) C(Ka,j,Jb)
         DO KB = 1, NKB
*. Number of nonvanishing connections from KB
            LL = 0
            KK = 0
            DO L = 1, NL
               IF(KBJB(KB,L).NE.0) LL = LL + 1
            END DO
            DO K = 1, NK
               IF(KBIB(KB,K).NE.0) KK = KK + 1
            END DO
*
            IF(KK.NE.0.AND.LL.NE.0) THEN
               DO K = 1, NK
                  IB = KBIB(KB,K)
                  IF(IB.NE.0) THEN
                     SGNK = XKBIB(KB,K)
                     DO L = 1, NL
                        JB = KBJB(KB,L)
                        IF(JB.NE.0) THEN
                           SGNL = XKBJB(KB,L)
                           FACTOR = SGNK*SGNL
*. We have now an IB and Jb string, let us do it
                           ISOFF = (IB-1)*NI*NKA + 1
                           ICOFF = (JB-1)*NJ*NKA + 1
C?                         write(6,*) 'ICOFF is now =  ',ICOFF
                           INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
                           IMAX = NI
*
                           IF(IKORD.NE.0) THEN
*. Restrict so (ij) .le. (kl)
                              IMAX  = K
                              JKINTOF = INTOF + (K-1)*NJ
C     CALL COPVEC(XIJKL(JKINTOF),XIJILS,NJ)
                              DO J = L,NL
                                 XIJILS(J) = XIJKL(JKINTOF-1+J)
                              END DO
                              XIJKL(JKINTOF-1+L)
     &                          = 0.5D0*XIJKL(JKINTOF-1+L)
                              DO J = L+1, NL
                                 XIJKL(JKINTOF-1+J) = 0.0D0
                              END DO
                           END IF
C     ONE = 1.0D0
                           CALL MATML7(SKII(ISOFF),CKJJ(ICOFF),
     &                          XIJKL(INTOF),NKA,IMAX,NKA,NJ,
     &                          NJ,IMAX,FACS,FACTOR ,0)
                           IF(IKORD.NE.0) THEN
                              DO J = L,NL
                                 XIJKL(JKINTOF-1+J) =  XIJILS(J)
                              END DO
C     CALL COPVEC(XIJILS,XIJKL(JKINTOF),NJ)
                           END IF
*
                        END IF
*     ^ JB.eq.0 or not
                     END DO
*                    ^ End of L to NL loop
*
                  END IF
*                 ^ IB.eq.0 or not
               END DO
*              ^ End of K to NK loop
            END IF
*           ^ KK, LL.eq.0 or not
         END DO
*        ^ End of loop over KB strings
*
      ELSE IF (IROUTE.EQ.1.OR.IROUTE.EQ.2) THEN

         WRITE(6,*) ' Sorry route 1and 2 have retired'
         Call Abend2( 'SKICKJ:Invalid route=1,2' )
      END IF
*. End of IROUTE branchning
*
C     CALL QEXIT('SKICK')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SXSTR(ISTRSM,ISTRTP,ISTRGP,
     &                  IOBSM,IOBTP,JOBSM,JOBTP,MXSXST,
     *                  ISXSTR,JSXSTR,IEXSTR,FACSTR,NEX,
     &                  IOFFDG,NTESTG)
*
* Obtain all excitations from a set of strings of given sym, type,
* and group that can be obtained by applying
* single excitations where each operator has a given sym and type
*
* IF IOFFDG .ne. 0 , only excitation a+i a j with i.ne.j are generated
*
*
* =====
* Input
* =====
*

*    ISTRSM,ISTTP,ISTRGP :  symmetry, type and group  of input string
*
*   IOBSM,IOBTP : symmetry and type of orbital i
*   JOBSM,JOBTP : symmetry and type of orbital j
*   MXSXST      : Max number of single excitations for given string
*
* ======
* Output
* ======
*     ISXSTR(IEX,ISTR) : I orbital indeces of SX
*     JSXSTR(IEX,ISTR) : J orbital indeces of SX
*     IEXSTR(IEX,ISTR) : Number of excited string ( relative to offset )
*     NEX(ISTR)        : Number of excitations
*     FACSTR(ISTR)     : Phase factor of excitation.
*
*     Jeppe Olsen, March 1994, LUCIA version
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. a few include blocks
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPOBS,MXPR4T),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),
C    &              NOBPTS(6+MXPR4T,MXPOBS),IOBPTS(6+MXPR4T,MXPOBS),
C    &              ITOOBS(MXPOBS),ITPFTO(MXPORB),ISMFSO(MXPORB)
C     COMMON/STRINP/NSTTYP,MNRS1(MXPSTT),MXRS1(MXPSTT),
C    &              MNRS3(MXPSTT),MXRS3(MXPSTT),NELEC(MXPSTT),
C    &              IZORR(MXPSTT),IAZTP,IBZTP,IARTP(3,10),IBRTP(3,10),
C    &              NZSTTP,NRSTTP,ISTTP(MXPSTT)
C     COMMON/STRBAS/KSTINF,KOCSTR(MXPSTT),KNSTSO(MXPSTT),KISTSO(MXPSTT),
C    &              KSTSTM(MXPSTT,2),KZ(MXPSTT),
C    &              KSTREO(MXPSTT),KSTSM(MXPSTT),KSTCL(MXPSTT),
C    &              KEL1(MXPSTT),KEL3(MXPSTT),KACTP(MXPSTT),
C    &              KCOBSM,KNIFSJ,KIFSJ,KIFSJO,KSTSTX
C    &             ,KNDMAP(MXPSTT),KNUMAP(MXPSTT)
*. A bit of local scratch
      DIMENSION IACAR(2),ITPAR(2)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ***************** '
        WRITE(6,*) ' Output from SXSTR '
        WRITE(6,*) ' ***************** '
      END IF
*
*. Number of orbitals and offsets
      NI   = NOBPTS(IOBTP,IOBSM)
      IOFF = IBTSOB(IOBTP,IOBSM)
*
      NJ   = NOBPTS(JOBTP,JOBSM)
      JOFF = IBTSOB(JOBTP,JOBSM)
*
*. Offset for input string
      ISTROF = IFRMR_REL(WORK,KISTSO(ISTRGP),
     &         (ISTRSM-1)*NOCTYP(ISTRGP)+ISTRTP)
*. Number of input strings
      NIST = IFRMR_REL(WORK,KNSTSO(ISTRGP),
     &         (ISTRSM-1)*NOCTYP(ISTRGP)+ISTRTP)
*. Group of output strings
      IACAR(1) = 2
      IACAR(2) = 1
*
      ITPAR(1) = IOBTP
      ITPAR(2) = JOBTP
CSKoldCALL NEWTYP_REL(ISTRGP,ISTRTP,IACAR,ITPAR,2,KSTRGP,KSTRTP)
*. Symmetry of output strings
      CALL SYMCOM_REL(3,0,IOBSM,JOBSM,IJSXSM)
      CALL SYMCOM_REL(3,0,ISTRSM,IJSXSM,KSTRSM)
*. Offset for output strings
      KSTROF = IFRMR_REL(WORK,KISTSO(ISTRGP),
     &         (KSTRSM-1)*NOCTYP(ISTRGP)+KSTRTP)
*
*. Type of mappings
*
*. Nel => Nel - 1 electrons
      IF(ISTAC(ISTRGP,1).NE.0.AND.ISTAC(ISTRGP,2).NE.0) THEN
*. full list
        IM1FUL = 1
        LM1    = NACOB
      ELSE
*. Truncated list
        IM1FUL = 0
        LM1 = NELEC(ISTRGP)
      END IF
*. Nel -1 => Nel mapping, things are not so simple anymore so
      IF(ISTAC(ISTRGP+1,1).NE.0.AND.ISTAC(ISTRGP+1,2).NE.0 ) THEN
        IP1FUL = 1
        LP1 = NACOB
      ELSE
*. Contains only creations, compact with indeces
        IP1FUL = 0
        LP1 = -1
      END IF
C?    write(6,*) ' SXSTR : ISTRGP,IP1FUL LP1' ,
C?   &                     ISTRGP,IP1FUL,LP1
      CALL  SXSTRS(ISTROF,NIST,KSTROF,
     &              NI,IOFF,NJ,JOFF,
     &              WORK(KSTSTM(ISTRGP,1)),WORK(KSTSTM(ISTRGP,2)),
     &              LM1,IM1FUL,
     &              WORK(KSTSTM(ISTRGP+1,1)),WORK(KSTSTM(ISTRGP+1,2)),
     &              LP1,IP1FUL,
     &              WORK(KSTSTMI(ISTRGP+1)), WORK(KSTSTMN(ISTRGP+1)),
     &              MXSXST,ISXSTR,JSXSTR,IEXSTR,FACSTR,NEX,IOFFDG,
     &              NTEST)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SXSTRS(ISTROF,NIST,KSTROF,
     &                   NI,IOFF,NJ,JOFF,
     &                   IAMAPO,IAMAPS,LAMAP,IAMPFL,
     &                   ICMAPO,ICMAPS,LCMAP,ICMPFL,
     &                   ICMPO,ICMPL,
     &                   MXSXST,ISXSTR,JSXSTR,IEXSTR,FACSTR,NEX,
     &                   IOFFDG,NTEST)
*
* Obtain single excitations from string ISTROF-ISTROF+NIST-1,
* Slave routine mastered by SXSTR
*
* ==================
*. Additional input ( compared to SXSTR)
* ==================
*
* ISTROF : Absolute number of first string to be excited from
* NIST   : Number of strings to be excited from
* KSTROF : Offset of strings in resulting type-symmetry block
* N*,*OFF,*=I,J : Number and offset for each orbital set
*
* IAMAPO : Annihilation mapping, orbital part
* IAMAPS : Annihilation mapping, string part
* LAMAP  : Row dimension of Annihilation map
* IAMPFL : Annihilation map complete ?
*
* ICMAPO : Creation     mapping, orbital part
* ICMAPS : Creation     mapping, string part
* ICMPFL : Creation map complete ?
* LCMAP  : Row dimension of Creation     map
*
* Jeppe Olsen, March 1994
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
*. Input
*
      INTEGER IAMAPO(LAMAP,*), IAMAPS(LAMAP,*)
C     INTEGER ICMAPO(LCMAP,*), ICMAPS(LCMAP,*)
      INTEGER ICMAPO(*),ICMAPS(*)
      INTEGER ICMPO(*),ICMPL(*)
*. Output
      INTEGER ISXSTR(MXSXST,*),JSXSTR(MXSXST,*)
      INTEGER IEXSTR(MXSXST,*)
      DIMENSION FACSTR(MXSXST,*)
      INTEGER NEX(*)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' SXFSTS : NTEST = ', NTEST
        WRITE(6,*) ' LCMAP LAMAP ', LCMAP,LAMAP
        WRITE(6,*) ' IAMPFL, ICMPFL ',IAMPFL, ICMPFL
        WRITE(6,*) ' MXSXST ', MXSXST
      END IF
*
      DO 1100 ISTR = ISTROF,ISTROF+NIST-1
C?      write(6,*) ' ISTR = ', ISTR
        LEX = 0
        DO 1002 JORB = JOFF,JOFF+NJ-1
C?        write(6,*) ' JORB = ', JORB
*
* =====================================
* 1 :  Remove orbital JORB from ISTR
* =====================================
*
          JOCC = 0
          IF(IAMPFL.EQ.1) THEN
*. Read from full map
            IF(IAMAPO(JORB,ISTR).EQ.-JORB) THEN
              JOCC = 1
              IF(IAMAPS(JORB,ISTR).GT.0) THEN
                JISTR = IAMAPS(JORB,ISTR)
                SJ = 1.0D0
               ELSE
                JISTR = -IAMAPS(JORB,ISTR)
                SJ = -1.0D0
              END IF
            END IF
          ELSE
*. Read from compact map
            DO JELEC = 1, LAMAP
              IF(IAMAPO(JELEC,ISTR).EQ.-JORB) THEN
                JOCC = 1
                IF(IAMAPS(JELEC,ISTR).GT.0) THEN
                   JISTR = IAMAPS(JELEC,ISTR)
                   SJ = 1.0D0
                ELSE
                   JISTR = -IAMAPS(JELEC,ISTR)
                   SJ = -1.0D0
                END IF
              END IF
            END DO
          END IF
C?        WRITE(6,*) ' JOCC = ', JOCC
          IF(JOCC.EQ.0) GOTO 1002
*
* ==================================
*. 2 : Add orbital I to string JISTR
* ==================================
*
          IF(ICMPFL.EQ.1) THEN
            JISTRO = (JISTR-1)*LCMAP+1
            NIEFF = LCMAP
          ELSE
            JISTRO = ICMPO(JISTR)
            NIEFF = ICMPL(JISTR)
          END IF
C?        WRITE(6,*) ' JISTRO NIEFF ', JISTRO,NIEFF
C         DO 1001 IORB = IOFF,IOFF+NI-1
          DO 1001 IIORB = 1, NIEFF
*
            IF(ICMPFL.EQ.1) THEN
              IORB = IIORB-1+IOFF
              IIIORB = IORB
            ELSE
              IORB = ABS(ICMAPO(JISTRO-1+IIORB))
              IIIORB = IIORB
              IF(IORB.LT.IOFF .OR. IORB .GT. IOFF + NI -1 )
     &        GOTO 1001
            END IF
*
C?          WRITE(6,*) ' IORB = ', IORB
            IF(IOFFDG.NE.0 .AND. IORB.EQ.JORB) GOTO 1001
            IJACT = 0
            IF(ICMAPO(JISTRO-1+IIIORB).EQ.+IORB) THEN
              IJACT  = 1
              IF(ICMAPS(JISTRO-1+IIIORB).GT.0) THEN
                IJISTR = ICMAPS(JISTRO-1+IIIORB)
                SIJ = SJ
              ELSE
                IJISTR = -ICMAPS(JISTRO-1+IIIORB)
                SIJ = -SJ
              END IF
            END IF
C?          WRITE(6,*) ' IIIORB IJISTR ', IIIORB, IJISTR
C?          WRITE(6,*) ' IJACT = ', IJACT
            IF(IJACT.EQ.0) GOTO 1001
*
*. A new excitation has been born, enlist it !
*
            LEX = LEX + 1
            ISXSTR(LEX,ISTR-ISTROF+1) = IORB
            JSXSTR(LEX,ISTR-ISTROF+1) = JORB
            IEXSTR(LEX,ISTR-ISTROF+1) = IJISTR-KSTROF+1
            FACSTR(LEX,ISTR-ISTROF+1) = SIJ
C?          WRITE(6,*) 'FACSTR = ',  FACSTR(LEX,ISTR-ISTROF+1)
 1001     CONTINUE
 1002   CONTINUE
        NEX(ISTR-ISTROF+1) = LEX
 1100 CONTINUE
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' **************** '
         WRITE(6,*) ' SXSTRS reporting '
         WRITE(6,*) ' **************** '
         WRITE(6,*)
         DO ISTR = ISTROF,ISTROF+NIST-1
           WRITE(6,*) ' excitations from string ',ISTR
           WRITE(6,*)
           WRITE(6,*) ' iorb jorb exc.string phase '
           WRITE(6,*) ' =========================='
           DO LEX = 1, NEX(ISTR-ISTROF+1)
             WRITE(6,'(2I4,I8,F8.3)')
     &       ISXSTR(LEX,ISTR-ISTROF+1),JSXSTR(LEX,ISTR-ISTROF+1),
     &       IEXSTR(LEX,ISTR-ISTROF+1),FACSTR(LEX,ISTR-ISTROF+1)
           END DO
           WRITE(6,*)
         END DO
      END IF
*
      RETURN
      END
