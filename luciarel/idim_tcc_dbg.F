      SUBROUTINE IDIM_TCC_DBG(ITSOSO_TP,NTSOSO_TP,ISYM,
     &           MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &           LTSOSO_TP,IBTSOSO_TP,IDIM_T,MSCOMB_CC,
     &           MX_TBLK_AS,IFTONE,IFTTWO,ILTTWO,N1ELINT,N2ELINT,IPRCC)
*
* Dimension of T operators in spin-orbital basis
*
* Largest number of strings of given sym in T(ICA,ICB,IAA,IAB)
* i.e. largest block of ICA, ICB,IAA,IAB of given sym 
*
* Size of block required to hold above blocks
*
* Jeppe Olsen, Summer of 99
*
* Double group version
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "csm.inc"
#include "symm.inc"
*. Specific input
      INTEGER ITSOSO_TP(4*NGAS,NTSOSO_TP)
*. Local scratch 
      INTEGER ICA_GRP(MXPNGAS), ICB_GRP(MXPNGAS)
      INTEGER IAA_GRP(MXPNGAS), IAB_GRP(MXPNGAS)
*.Output : Length of each type and offset for each type
      INTEGER LTSOSO_TP(NTSOSO_TP),IBTSOSO_TP(NTSOSO_TP)
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRCC)
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' IDIM_TCC : MSCOMB_CC, ISYM = ', MSCOMB_CC,ISYM
      END IF
*
      LENGTH = 0
      LENGTH2 = 0
      MX_ST_TSOSO = 0
      MX_ST_TSOSO_BLK = 0
      MX_TBLK = 0
      MX_SBSTR = 0
      MX_TBLK_AS = 0
*
      DO ITSS_TP = 1, NTSOSO_TP
       IF(NTEST.GE.100) WRITE(6,*) ' ITSS_TP = ', ITSS_TP 
*. Occupation to group translation 
       CALL OCC_TO_GRP(ITSOSO_TP(1+0*NGAS,ITSS_TP),ICA_GRP,1)
       CALL OCC_TO_GRP(ITSOSO_TP(1+1*NGAS,ITSS_TP),ICB_GRP,1)
       CALL OCC_TO_GRP(ITSOSO_TP(1+2*NGAS,ITSS_TP),IAA_GRP,1)
       CALL OCC_TO_GRP(ITSOSO_TP(1+3*NGAS,ITSS_TP),IAB_GRP,1)
*
       if (NTEST.ge.50) then
         write(6,*) ' Output from OCC_TO_GRP : '
         write(6,*) ' ICA_GRP : ',(ICA_GRP(I),I=1,NGAS,1)
         write(6,*) ' ICA_GRP : ',(ICB_GRP(I),I=1,NGAS,1)
         write(6,*) ' ICA_GRP : ',(IAA_GRP(I),I=1,NGAS,1)
         write(6,*) ' ICA_GRP : ',(IAB_GRP(I),I=1,NGAS,1)
       end if
*
       NEL_CA = IELSUM(ITSOSO_TP(1+0*NGAS,ITSS_TP),NGAS)
       NEL_CB = IELSUM(ITSOSO_TP(1+1*NGAS,ITSS_TP),NGAS)
       NEL_AA = IELSUM(ITSOSO_TP(1+2*NGAS,ITSS_TP),NGAS)
       NEL_AB = IELSUM(ITSOSO_TP(1+3*NGAS,ITSS_TP),NGAS)
*
       IBTSOSO_TP(ITSS_TP) = LENGTH + 1
       DO I_CR_SM = 1, NSMST
*. symmetry of annihilation strings
        I_AN_SM = IDBGMULT(ISYM,INVELM(I_CR_SM))
*. symmetry of creation strings corresponding to annihilation strings ...
        I_AN_SM = IADJSYM(I_AN_SM)
        DO I_CR_AL_SM = 1, NSMST
          DO II_AN_AL_SM = 1, NSMST
             I_AN_AL_SM = IADJSYM(II_AN_AL_SM)
*
             I_CR_BE_SM = IDBGMULT(INVELM(I_CR_AL_SM),I_CR_SM)
             I_AN_BE_SM = IDBGMULT(INVELM(I_AN_AL_SM),I_AN_SM)
*
             CALL NST_SPGRP_DBG(NGAS,ICA_GRP,I_CR_AL_SM,NSMST,
     &                          LEN_CA,NDIST_CA,1)
             CALL NST_SPGRP_DBG(NGAS,ICB_GRP,I_CR_BE_SM,NSMST,
     &                          LEN_CB,NDIST_CB,2)
             CALL NST_SPGRP_DBG(NGAS,IAA_GRP,I_AN_AL_SM,NSMST,
     &                          LEN_AA,NDIST_AA,1)
             CALL NST_SPGRP_DBG(NGAS,IAB_GRP,I_AN_BE_SM,NSMST,
     &                          LEN_AB,NDIST_AB,2)
*
             IF(NTEST.GE.1000) THEN
               WRITE(6,'(A,4I5)') ' sym of CA CB AA AB',
     &         I_CR_AL_SM, I_CR_BE_SM, I_AN_AL_SM, I_AN_BE_SM
               WRITE(6,'(A,4I5)') ' LEN_CA, LEN_CB, LEN_AA, LEN_AB ',
     &                      LEN_CA, LEN_CB, LEN_AA, LEN_AB
             END IF
*
             LENGTH = LENGTH + LEN_CA*LEN_CB*LEN_AA*LEN_AB
             LENGTH2 = LENGTH2 + LEN_CA*LEN_CB*LEN_AA*LEN_AB
             if (NTEST.ge.2000) then
               write(6,*) 'LENGTH ',LEN_CA*LEN_CB*LEN_AA*LEN_AB
             end if
*
             MX_ST_TSOSO = 
     &       MAX(MX_ST_TSOSO,LEN_CA,LEN_CB,LEN_AA,LEN_AB)
             MX_ST_TSOSO_BLK = 
     &       MAX(MX_ST_TSOSO_BLK,LEN_CA*NEL_CA,LEN_CB*NEL_CB,
     &                           LEN_AA*NEL_CA,LEN_AB*NEL_CB)
             MX_TBLK = MAX(MX_TBLK,LEN_CA*LEN_CB*LEN_AA*LEN_AB)
             MX_SBSTR = MAX(MX_SBSTR,LEN_CA,LEN_CB,LEN_AA,LEN_AB)
          END DO
        END DO
       END DO
      LTSOSO_TP(ITSS_TP) = LENGTH2
      MX_TBLK_AS = MAX(MX_TBLK_AS,LENGTH2)
      END DO
*
      IDIM_T = LENGTH
*
* Number of non-spin-flip one- and two electron integrals 
* (for diagonal):
      N1ELINT = IBTSOSO_TP(IFTTWO) - IBTSOSO_TP(IFTONE)
      N2ELINT = LENGTH - N1ELINT
*
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' Number of T-coefficients ', LENGTH
        WRITE(6,*) ' Largest symmetry block of T ', MX_TBLK
        WRITE(6,*) ' Largest block of T ', MX_TBLK_AS
        WRITE(6,*) ' Largest substring : ', MX_SBSTR
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) 
     &  ' Largest number of strings of given sym and type ',
     &               MX_ST_TSOSO
        WRITE(6,*) 
     & ' Largest block of strings of given sym and type ',
     &               MX_ST_TSOSO_BLK
        WRITE(6,*) ' Number of elements per block '
        CALL IWRTMA(LTSOSO_TP,1,NTSOSO_TP,1,NTSOSO_TP)
        WRITE(6,*) ' Offset for each block '
        CALL IWRTMA(IBTSOSO_TP,1,NTSOSO_TP,1,NTSOSO_TP)
      END IF
*
      RETURN
      END 
