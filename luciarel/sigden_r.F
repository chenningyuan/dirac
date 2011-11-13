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
      SUBROUTINE ANA_GENCC(T,ISM)
*
* Analyze T-coefficients of general coupled cluster wavefunction
*
* Jeppe Olsen, September 1999 ( modified GASANA)
*
*
#include "implicit.inc"
#include "ipoist8.inc"
*. specific input Coefficients and types
      DIMENSION T(*)
*
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "strbas.inc"
#include "cicisp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "ctcc.inc"
C     INCLUDE 'cc_exc.inc'
*
      CALL QENTER('ANAGC')
      call memmar(KLOFF,IDUMMY,'MARK  ',DUMMY,'ANAGCC')
*
*. Four blocks of string occupations
      call memmar(KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC1')
      call memmar(KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC2')
      call memmar(KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC3')
      call memmar(KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC4')
*. Number of terms to be printed 
      IF(IPRNCIV.EQ.0) THEN
        THRES = 0.02 
        MAXTRM = 200
      ELSE
        THRES = 0.0D0
*. Well atmost 100000 coefs - to save disk ..
        MAXTRM = 100000
      END IF
*Numbers and weight of spinorbital excitation type 
      LENGTH_SOX = NSPOBEX_TP*10
      call memmar(KNCPMT_SOX,LENGTH_SOX,'ADDL  ',1,'KNCPMT')
      call memmar(KWCPMT_SOX,LENGTH_SOX,'ADDL  ',2,'KWCPMT')
*. Number and weight of orbital excitation type
      LENGTH_OX = NOBEX_TP*10
      call memmar(KNCPMT_OX,LENGTH_OX,'ADDL  ',1,'KNCPMT')
      call memmar(KWCPMT_OX,LENGTH_OX,'ADDL  ',1,'KNCPMT')
      IUSLAB = 0
*
* FIXME
* !!! Define KLOBEX,KLSOX_TO_OX when debugging !!!
*
*. Occupation of strings of given sym and supergroup
      CALL ANA_GENCCS(T,WORK(KLSOBEX),WORK(KLOBEX),
     &            NSPOBEX_TP,NOBEX_TP,ISM,
     &            THRES,MAXTRM,
     &            WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
     &            WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),
     &            IUSLAB,
     &            IDUMMY,WORK(KNCPMT_SOX),WORK(KWCPMT_SOX),      
     &            WORK(KNCPMT_OX),WORK(KWCPMT_OX),      
     &            NTOOB,IPRNCIV,WORK(KLSOX_TO_OX))
   
      call memmar(KDUM,IDUM,'FLUSM',IDUM,'GASGCC')
      CALL QEXIT('ANAGC')
*
      RETURN
      END
      SUBROUTINE ANA_GENCCS(T,ISPOBEX_TP,IOBEX_TP,
     &                 NSPOBEX_TP,NOBEX_TP,ISM,
     &                 THRES,MAXTRM,
     &                 IOCC_CA, IOCC_CB, IOCC_AA, IOCC_AB,
     &                 IUSLAB,
     &                 IOBLAB,NCPMT_SOX,WCPMT_SOX,
     &                        NCPMT_OX, WCPMT_OX,
     &                 NORB,IPRNCIV,ISOX_TO_OX)
*
* Analyze T-CC vector :
*
*      1) Print atmost MAXTRM  operators with coefficients
*         larger than THRES
*
*      2) Number of coefficients in given range
*
*      3) Number of coefficients in given range for given 
*         excitation type, spin orbital types and orbital types
*
* Jeppe Olsen , September 1999 
*               July 2000, MSCOMB_CC added
*                                                  

*. If IUSLAB  differs from zero Character*6 array IOBLAB is used to identify
*  Orbitals
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "symm.inc"
#include "csm.inc"
#include "orbinp.inc"
*. Specific input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER IOBEX_TP(2*NGAS,NOBEX_TP)
      DIMENSION T(*)
      DIMENSION ISOX_TO_OX(NSPOBEX_TP)
*. Scratch
      INTEGER IOCC_CA(*),IOCC_CB(*),IOCC_AA(*),IOCC_AB(*)
*. Local scratch
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS) 
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
      CHARACTER*6 IOBLAB(*)
*. Output
      DIMENSION NCPMT_SOX(10,NSPOBEX_TP)                           
      DIMENSION WCPMT_SOX(10,NSPOBEX_TP)                          
      DIMENSION NCPMT_OX(10,NOBEX_TP)                           
      DIMENSION WCPMT_OX(10,NOBEX_TP)                          
*
      CALL ISETVC(NCPMT_SOX,0    ,10*NSPOBEX_TP)
      CALL SETVEC(WCPMT_SOX,0.0D0,10*NSPOBEX_TP)
      CALL ISETVC(NCPMT_OX,0    ,10*NOBEX_TP)
      CALL SETVEC(WCPMT_OX,0.0D0,10*NOBEX_TP)
*
* ===========================================================
*.1 : Printout of coefficients and largest occupation vectors       
*.2 : Group coefficients by type and size
* ===========================================================
*
      WRITE(6,*)
      WRITE(6,*) ' Operators are written as : '
      WRITE(6,*)
      WRITE(6,*)   ' Creation of alpha '
      WRITE(6,*)   ' Creation of beta '
      WRITE(6,*)   ' Annihilation of alpha '
      WRITE(6,*)   ' Annihilation of beta '
      WRITE(6,*)
      MINPRT = 0
      ITRM = 0
      ILOOP = 0
      NTVAR = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      TNORM = 0.0D0
      TTNORM = 0.0D0
2001  CONTINUE
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = 1.0D0
        XMIN = 1.0D0/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      IF(IPRNCIV.EQ.1) THEN
*. Print in one shot
       XMAX = 3006.1956
       XMIN =-3006.1956
      END IF
*
      WRITE(6,*)
      WRITE(6,'(A,E10.4,A,E10.4)')
     &'  Printout of coefficients in interval  ',XMIN,' to ',XMAX
      WRITE(6,'(A)')
     &'  ========================================================='
      WRITE(6,*)
*
      IT = 0
      IIT = 0
      DO ITSS = 1, NSPOBEX_TP
C?     WRITE(6,*) ' NSPOBEX_TP = ', NSPOBEX_TP
*. Transform from occupations to groups
       CALL OCC_TO_GRP(ISPOBEX_TP(1+0*NGAS,ITSS),IGRP_CA,1)
       CALL OCC_TO_GRP(ISPOBEX_TP(1+1*NGAS,ITSS),IGRP_CB,1)
       CALL OCC_TO_GRP(ISPOBEX_TP(1+2*NGAS,ITSS),IGRP_AA,1)
       CALL OCC_TO_GRP(ISPOBEX_TP(1+3*NGAS,ITSS),IGRP_AB,1)
*
       NEL_CA = IELSUM(ISPOBEX_TP(1+0*NGAS,ITSS),NGAS)
       NEL_CB = IELSUM(ISPOBEX_TP(1+1*NGAS,ITSS),NGAS)
       NEL_AA = IELSUM(ISPOBEX_TP(1+2*NGAS,ITSS),NGAS)
       NEL_AB = IELSUM(ISPOBEX_TP(1+3*NGAS,ITSS),NGAS)
*. Diagonal restricted type of spinorbital excitation ?
C      IF(MSCOMB_CC.EQ.1) THEN
C        CALL DIAG_EXC_CC(
C    &        ISPOBEX_TP(1+0*NGAS,ITSS),ISPOBEX_TP(1+1*NGAS,ITSS),
C    &        ISPOBEX_TP(1+2*NGAS,ITSS),ISPOBEX_TP(1+3*NGAS,ITSS),
C    &        NGAS,IDIAG)
C      ELSE
         IDIAG = 0
C      END IF
*
       DO ISM_C = 1, NSMST
        ISM_A = IDBGMULT(ISM,INVELM(ISM_C)) 
        DO ISM_CA = 1, NSMST
         ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
         DO ISM_AA = 1, NSMST
          ISM_AB =  IDBGMULT(ISM_A,INVELM(ISM_AA))
          ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
          ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
*
          IF(IDIAG.EQ.1.AND.ISM_ALPHA.LT.ISM_BETA) GOTO 777
          IF(IDIAG.EQ.1.AND.ISM_ALPHA.EQ.ISM_BETA) THEN 
           IRESTRICT_LOOP = 1
          ELSE
           IRESTRICT_LOOP = 0
          END IF
*. obtain strings
          IUB = 1
          CALL GETSTR2_TOTSM_SPGP_REL(IUB,IGRP_CA,NGAS,ISM_CA,NEL_CA,
     &                                NSTR_CA,
     &                                IOCC_CA,NORB,0,IDUM,IDUM)
          IUB = 2
          CALL GETSTR2_TOTSM_SPGP_REL(IUB,IGRP_CB,NGAS,ISM_CB,NEL_CB,
     &                                NSTR_CB,
     &                                IOCC_CB,NORB,0,IDUM,IDUM)
          IUB = 1
          CALL GETSTR2_TOTSM_SPGP_REL(IUB,IGRP_AA,NGAS,ISM_AA,NEL_AA,
     &                                NSTR_AA,
     &                                IOCC_AA,NORB,0,IDUM,IDUM)
          IUB = 2
          CALL GETSTR2_TOTSM_SPGP_REL(IUB,IGRP_AB,NGAS,ISM_AB,NEL_AB,
     &                                NSTR_AB,
     &                                IOCC_AB,NORB,0,IDUM,IDUM)
*. Loop over T elements as  matrix T(I_CA, I_CB, IAA, I_AB)
          DO I_AB = 1, NSTR_AB
           IF(IRESTRICT_LOOP.EQ.1) THEN
             I_AA_MIN = I_AB
           ELSE
             I_AA_MIN = 1
           END IF
           DO I_AA = I_AA_MIN, NSTR_AA
            DO I_CB = 1, NSTR_CB
             IF(IRESTRICT_LOOP.EQ.1.AND.I_AA.EQ.I_AB) THEN
               I_CA_MIN = I_CB
             ELSE
               I_CA_MIN = 1
             END IF
             DO I_CA = I_CA_MIN, NSTR_CA
              IT = IT + 1
C?            WRITE(6,*) ' IT, T(IT) = ', IT,T(IT)
*
              IF(ILOOP .EQ. 1 ) THEN 
                NTVAR = NTVAR + 1
                TNORM = TNORM + T(IT)**2
*. Classify element according to size 
                DO IPOT = 1, 10
                  IF(10.0D0 ** (-IPOT+1).GE.ABS(T(IT)).AND.
     &            ABS(T(IT)).GT. 10.0D0 ** ( - IPOT )      )THEN
                    IOEX = ISOX_TO_OX(ITSS)
                    NCPMT_OX(IPOT,IOEX)= NCPMT_OX(IPOT,IOEX)+1  
                    WCPMT_OX(IPOT,IOEX)= WCPMT_OX(IPOT,IOEX)+T(IT)**2
                  END IF
                END DO
*             ^ End of loop over powers of ten
              END IF
*             ^ end if we are in loop 1
              IF( XMAX .GE. ABS(T(IT)) .AND.
     &        ABS(T(IT)).GT. XMIN .AND. ITRM.LE.MAXTRM) THEN
                ITRM = ITRM + 1
                IIT = IIT + 1
                TTNORM = TTNORM + T(IT) ** 2
                WRITE(6,'(A)')
                WRITE(6,'(A)')
     &          '                 =================== '
                WRITE(6,*)

                WRITE(6,'(A,2I8,2X,E14.8)')
     &          '  Type, number, size of amplitude : ',
     &          ITSS, IT, T(IT)
                IF(IUSLAB.EQ.0) THEN
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_CA(IEL+(I_CA-1)*NEL_CA),IEL = 1, NEL_CA)
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_CB(IEL+(I_CB-1)*NEL_CB),IEL = 1, NEL_CB)
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_AA(IEL+(I_AA-1)*NEL_AA),IEL = 1, NEL_AA)
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_AB(IEL+(I_AB-1)*NEL_AB),IEL = 1, NEL_AB)
                END IF
              END IF
*             ^ End if could and should be printed
             END DO
*            ^ End of loop over alpha creation strings
            END DO
*           ^ End of loop over beta creation strings
           END DO
*          ^ End of loop over alpha annihilation 
          END DO 
*         ^ End of loop over beta annihilation 
  777    CONTINUE
         END DO
        END DO
       END DO
*      ^ End of loop over symmetry blocks
      END DO
*     ^ End of loop over over types of excitations
       IF(IIT .EQ. 0 ) WRITE(6,*) '   ( no coefficients )'
       IF( XMIN .GT. THRES .AND. ILOOP .LE. 20 ) GOTO 2001
*
       WRITE(6,'(A,E15.8)')
     & '  Norm of printed coefficients .. ', TTNORM
       WRITE(6,'(A,E15.8)')
     & '  Norm of all     coefficients .. ',  TNORM
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') '   Magnitude of T coefficients '
      WRITE(6,'(A)') '  =============================='
      WRITE(6,'(A)')
      WACC = 0.0D0
      NACC = 0
      DO 300 IPOT = 1, 10
        W = 0.0D0
        N = 0
        DO 290 ITSS = 1, NOBEX_TP
            N = N + NCPMT_OX(IPOT,ITSS)                    
            W = W + WCPMT_OX(IPOT,ITSS)                    
  290   CONTINUE
        WACC = WACC + W
        NACC = NACC + N
        WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &  '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
  300 CONTINUE
*
      WRITE(6,*) ' Number of coefficients less than  10-11',
     &           ' IS  ',NTVAR - NACC
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Magnitude of CI coefficients for each type of excitation '
      WRITE(6,'(A)') 
     & '  ========================================================='
      WRITE(6,'(A)')
      DO 400 ITSS   = 1, NOBEX_TP 
          N = 0
          DO 380 IPOT = 1, 10
            N = N + NCPMT_OX(IPOT,ITSS)                     
  380     CONTINUE
          IF(N .NE. 0 ) THEN
            WRITE(6,*) ' Orbital excitation type = ', ITSS
            WRITE(6,'(A,I9)')  
     &      '         Number of coefficients larger than 10-11 ', N
            WRITE(6,*)
            WACC = 0.0D0
            DO IPOT = 1, 10
              N =  NCPMT_OX(IPOT,ITSS)                    
              W =  WCPMT_OX(IPOT,ITSS)                    
              WACC = WACC + W
              WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &        '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
            END DO
          END IF 
  400 CONTINUE
*
*. Total weight and number of dets per excitation level
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Total weight and number of SD''s (> 10 ** -11 )  : '          
      WRITE(6,'(A)') 
     & '  ================================================='
      WRITE(6,'(A)')
      WRITE(6,*) '        N      Weight      Acc. Weight   Exc. type '
      WRITE(6,*) ' ==================================================='
      WACC = 0.0D0
      DO ITSS = 1, NOBEX_TP  
          N = 0
          W = 0.0D0
          DO IPOT = 1, 10
            N = N + NCPMT_OX(IPOT,ITSS)                   
            W = W + WCPMT_OX(IPOT,ITSS)                   
          END DO
          WACC = WACC + W
          IF(N .NE. 0 ) THEN
            WRITE(6,'(1X,I9,3X,E9.4,7X,E9.4,2X,I3)') 
     &      N,W,WACC,ITSS                                   
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
      SUBROUTINE CAAB_INCLUDED(IEXOP,INCLUDED_CAAB)
*
* Is excitation operator CAAB included in the 
* list of excitations
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
*. Specific input
      INTEGER IEXOP(NGAS,4)
*
      INCLUDED_CAAB = 1
* Creation alpha
      CALL OCC_INCLUDED(IEXOP(1,1),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
* Creation beta
      CALL OCC_INCLUDED(IEXOP(1,2),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
* Annihilation alpha
      CALL OCC_INCLUDED(IEXOP(1,3),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
* Annihilation beta
      CALL OCC_INCLUDED(IEXOP(1,4),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
*
#if defined LUCI_DEBUG
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' CAAB excitation operator '
       CALL WRT_SPOX_TP(IEXOP,1)
       WRITE(6,*) ' INCLUDED_CAAB = ', INCLUDED_CAAB
      END IF
#endif
*
      RETURN
      END

    
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE CCEX_OCC_OCC(INOCC,IOUTOCC,NGAS,IAC,ICCEXC,IZERO)
*
* A set of occupations INOCC is given
*
* Find occupation generated by applying coupled cluster operator ICCEXC
* in either creation or annihilation form
*
* Jeppe Olsen, Summer of 99
*
#include "implicit.inc"
#include "ipoist8.inc"
*. Input
      INTEGER INOCC(NGAS)
      INTEGER ICCEXC(NGAS)
*. Output
      INTEGER IOUTOCC(NGAS)
*
      IZERO = 0
      DO JGAS = 1, NGAS
        IF(IAC.EQ.1) THEN
          JDELTA = - ICCEXC(JGAS)
        ELSE
          JDELTA = ICCEXC(JGAS)
        END IF
        IOUTOCC(JGAS) = INOCC(JGAS) + JDELTA
        IF(IOUTOCC(JGAS).LT.0) IZERO = 1
      END DO
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        write(6,*) 'CCEX_OCC_OCC speaking:'
        WRITE(6,*) ' Input occupation '   
        CALL IWRTMA(INOCC,1,NGAS,1,NGAS)
        IF(IAC.EQ.1) THEN
          WRITE(6,*) ' String of annihilation operators'
        ELSE
          WRITE(6,*) ' String of creation operators'
        END IF
        CALL IWRTMA(ICCEXC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Output string '
        CALL IWRTMA(IOUTOCC,1,NGAS,1,NGAS)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE DIAG_EXC_CC(NCA,NCB,NAA,NAB,NGAS,IDIAG)
*
* Check of alpha and betaexcitation parts are identical
*
*. Output
* IDIAG = 0 : Not diagonal in alpha, beta
* IDIAG = 1 : Diag in alpha, beta
*
*
#include "implicit.inc"
#include "ipoist8.inc"
*. Input
      INTEGER NCA(NGAS),NCB(NGAS),NAA(NGAS),NAB(NGAS)
*
      IDIAG = 1 
      DO IGAS = 1, NGAS
        IF(NCA(IGAS).NE.NCB(IGAS)) IDIAG = 0
        IF(NAA(IGAS).NE.NAB(IGAS)) IDIAG = 0
      END DO
*
      NTEST = 00
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCA, NCB, NAA, NAB = '
        CALL IWRTMA(NCA,1,NGAS,1,NGAS)
        CALL IWRTMA(NCB,1,NGAS,1,NGAS)
        CALL IWRTMA(NAA,1,NGAS,1,NGAS)
        CALL IWRTMA(NAB,1,NGAS,1,NGAS)
        WRITE(6,*) ' IDIAG = ', IDIAG
      END IF
#endif
*
      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_HX_RELA(IFLAG,IHTYPE,NO_TTWO,NOBPT,NHX,
     &                       IFTONE,IFTTWO,ILTTWO,IFSFLT1,ILSFLT1,
     &                       IHX_LOC,IHX,IHINDX,
     &                       SIGN_NHX,I2DIATP,LABEL,IOP_REO,SIGN_OPREO,
     &                       IPRNT)
C***********************************************************************
C
C     General purpose: obtain number (IFLAG = 1) or number and types (IFLAG.NE.1) 
C     of excitations in relativistic Hamiltonian 
C
C     written by Jeppe Olsen and Timo Fleig - March/April 2001
C
C     last revisions:
C
C     Stefan Knecht Nov 2008 - fix one-electron operator matrix elements for 
C                              correct behavior under time-reversal
C                              symmetry. 
C
C     Stefan Knecht Dec 2009 - removed introduced "time-reversal
C                              correction". With my new and proper
C                              one-electron integral fetching we get the 
C                              correct integrals.
C============================================================================
C  
C     Hans Joergen Aa. Jensen and Stefan Knecht Mar 2008
C       NO_TTWO.gt.0 for only one-electron terms
C
C     Timo Fleig, June 2001
C       orbital and electron number limitations
C
C     JO and TF - Oct 2001
C       updated to exclude unwanted components
C       included LABEL specification of excitation type for DIRAC
C       (and other interfaces)
C============================================================================
C
C IHTYPE = 1 : Include only Delta MK = 0 terms
C IHTYPE = 2 : Include Delta Mk = 0 for two-electron part 
C              and complete one-electron part
C IHTYPE = 3 : Include all terms for relativisitic Hamiltonian
C              in a basis of non-relativistic orbitals
C IHTYPE = 4 : Include all terms for relativisitic Hamiltonian
C              in a basis of relativistic spinors
C
C LABEL :  1 for unbarred (alpha)
C         -1 for barred   (beta)
C         and refers to integral (not operator!) indexing !!
C
C IOP_REO : Reordering array for output density elements.
C           Relate definition of operator in this routine with
C           general output CA AA CB AB (or corresponding scheme
C           in DMK .ne. 0 cases). Not used in 1-el. case, but
C           due to NHX type counting need to be defined.
C
C SIGN_OPREO : Sign for this operator reordering
C          1.0 or -1.0 : Non-redundant density and sign
C          0.0         : Redundant density type
C
C
C An excitation is stored as : 
C (alpha creation, beta creation, alpha annihilation, beta annihilation)
C
C***********************************************************************
C
#include "implicit.inc"
#include "ipoist8.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "mxpdim.inc"
#include "cgas.inc"
C     INPUT
C     -----
      DIMENSION NOBPT(NGAS)
C     OUTPUT
C     ------
      INTEGER IHX(NGAS,4,*),IHINDX(4,*),LABEL(4,*),IOP_REO(4,*)
      DIMENSION SIGN_NHX(*), SIGN_OPREO(*)
C     local scratch
      INTEGER IHX_LOC(4*NGAS)
C
C
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
C
      ZERO = 0.D0
      ONE = 1.0D0
      DMONE = -1.0D0
C
C     WRITE(LUWRT,*) 'GET_HX_RELA: IHTYPE is',IHTYPE
C
C
      NHX = 0
C
C Determine maximum number of electrons per GAS first:
C FIXME  : Hardwired only one type of GASSPC calculation for the
C          moment. This would require other changes as well.
      MXNELGS(1) = IGSOCCX(1,2,1)
      do IGAS=2,NGAS,1
        MXNELGS(IGAS) = IGSOCCX(IGAS,2,1) - IGSOCCX(IGAS-1,1,1)
        MXNELGS(IGAS) = min(MXNELGS(IGAS),2*NOBPT(IGAS))
      end do
C
C  Delta MK = 0 terms
C
C Operator: a+i alpha a j alpha
C           where i belongs to IGAS and j belongs to JGAS
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          if (NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1) then
              CALL IZERO(IHX_LOC,4*NGAS)
              IHX_LOC(IGAS+(1-1)*NGAS) = 1
              IHX_LOC(JGAS+(3-1)*NGAS) = 1
*. Is this excitation type included ?
              CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN 
                WRITE(LUWRT,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(IHINDX(1,NHX),NHX) = 1
                  LABEL(IHINDX(2,NHX),NHX) = 1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
C                 WRITE(LUWRT,*) ' alpha alpha type ',NHX
                end if
              END IF
            end if
          end if
        END DO
      END DO
*
* operator: a+i beta a j beta
*           where i belongs to IGAS and j belongs to JGAS
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          if (NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1) then
              CALL IZERO(IHX_LOC,4*NGAS)
              IHX_LOC(IGAS+(2-1)*NGAS) = 1
              IHX_LOC(JGAS+(4-1)*NGAS) = 1
              CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN 
                WRITE(LUWRT,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
C                 WRITE(LUWRT,*) ' beta beta: NHX, SIGN_NHX = ',
C    &                             NHX,SIGN_NHX(NHX)
                  LABEL(IHINDX(1,NHX),NHX) = -1
                  LABEL(IHINDX(2,NHX),NHX) = -1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
C                 WRITE(LUWRT,*) ' beta beta type ',NHX
                end if
              END IF
            end if
          end if
        END DO
      END DO
*
* For general purposes: First type with 
*                       1-electron and 2-electron integrals:
      IFTONE = 1
C     if only one-elecron operator (NO_TTWO > 0): spin-flip 
C     integral types (one-electron integrals) are included as 
C     2e-integrals depending on the operator symmetry...
C     even for pure one-electron operators... - SK: 18 Sep 08
      IFTTWO = NHX + 1
C
      IF (NO_TTWO .gt. 0) GOTO 200
*
* operator:
*. a+i alpha  a+j alpha  a k alpha  a l alpha, ordered so i.ge.j, k.ge.l
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS 
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1))) then
                CALL IZERO(IHX_LOC,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS) + 1
                IHX_LOC(KGAS+(3-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS) + 1
C               IHX(IGAS,1,NHX) = 1
C               IHX(JGAS,1,NHX) = IHX(JGAS,1,NHX)+1
* means IGAS=JGAS, so the number of creators of this type is 2!
C               IHX(KGAS,3,NHX) = 1
C               IHX(LGAS,3,NHX) = IHX(LGAS,3,NHX)+1
                CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN 
                  WRITE(LUWRT,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  a | +a  a)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = 1
                    LABEL(IHINDX(3,NHX),NHX) = 1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 2
                    SIGN_OPREO(NHX) = -ONE
                  end if
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
* operator:
*. a+i beta a+j beta a k beta a l beta, ordered so i.ge.j, k.ge.l
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS 
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1))) then
                CALL IZERO(IHX_LOC,4*NGAS)
                IHX_LOC(IGAS+(2-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS) + 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS) + 1
                CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN 
                  WRITE(LUWRT,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  b | +b  b)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = -1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = -1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 2
                    SIGN_OPREO(NHX) = -ONE
                  end if
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
* First type of integrals for 2-electron diagonal
  200 I2DIATP = NHX + 1
      IF (NO_TTWO .GT. 0) GO TO 300
C     ... skip next set of two-electron terms
*
* operator:
*. a+i alpha  a+j beta  a k beta  a l alpha
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, NGAS 
              if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &            MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &            NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &            NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1) then
                CALL IZERO(IHX_LOC,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = 1
                CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN 
                  WRITE(LUWRT,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  a | +b  b)
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    SIGN_NHX(NHX) = -1.D0 * ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 2
                    IOP_REO(4,NHX) = 3
                    SIGN_OPREO(NHX) = -1.D0 * ONE
                  end if
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
      if (IHTYPE.ge.4) then
*
* operator:
*. a+i alpha  a+j beta  a k alpha  a l beta
*   Generates type (ub|bu) of integrals.
*
        do IGAS = 1,NGAS,1
          do JGAS = 1,NGAS,1
            do KGAS = 1,NGAS,1
              do LGAS = 1,NGAS,1
                if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &              MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1) then
                  CALL IZERO(IHX_LOC,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = 1
                  call caab_included(IHX_LOC,INCLUDED)
                  if (INCLUDED.eq.0) then
                    write(LUWRT,*) 'Hamiltonian type excluded '
                    call wrt_spox_tp(IHX_LOC,1)
                  else
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      call icopve(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  b | +b  a)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 3
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  end if
                end if
              end do
            end do
          end do
        end do
*
      end if
*
* Last type of non-spin-flip integrals
  300 ILTTWO = NHX
      IFSFLT1= ILTTWO + 1
*
* one-electron spin orbit
*
      IF (IHTYPE.GT.2) THEN
*
* operator: a+i alpha  j beta
*           where i belongs to IGAS and j belongs to JGAS
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &          NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
              CALL IZERO(IHX_LOC,4*NGAS)
              IHX_LOC(IGAS+(1-1)*NGAS) = 1
              IHX_LOC(JGAS+(4-1)*NGAS) = 1
              CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN 
                WRITE(LUWRT,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(1,NHX) = 1
                  LABEL(2,NHX) = -1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
C                 WRITE(LUWRT,*) ' alpha beta spin flip type ',NHX
                end if
              END IF
            END IF
          END DO
        END DO
*
* operator: a+i beta a j alpha 
*           where i belongs to IGAS and j belongs to JGAS
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &          NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
              CALL IZERO(IHX_LOC,4*NGAS)
              IHX_LOC(IGAS+(2-1)*NGAS) = 1
              IHX_LOC(JGAS+(3-1)*NGAS) = 1
              CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN 
                WRITE(LUWRT,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  SIGN_NHX(NHX) = ONE
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  LABEL(1,NHX) = -1
                  LABEL(2,NHX) = 1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
C                 WRITE(LUWRT,*) ' beta alpha spin flip type ',NHX
                end if
              END IF
            END IF
          END DO
*         quaternionic matrix groups only!
        END DO
      END IF
*     last one-electron spin-flip type integral
      ILSFLT1 = NHX
*
*     two-electron terms with Delta mk .ne. 0
*
      if (NO_TTWO .eq. 0) THEN
*         include two-electron terms when NO_TTWO .gt. 0
      if (IHTYPE.eq.3.or.IHTYPE.eq.5) THEN
*
**********************************************************
*  Spin-dependent spinor basis. E.g. DIRAC environment   *
*    (ub|uu), (bu|uu), (bu|bb), (ub|bb) integral classes *
**********************************************************
*
* operator:
*. a+i alpha a+j alpha a k alpha a l beta, i .ge. j
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS 
                if ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &               MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &               NOBPT(IGAS).ge.2.and.NOBPT(KGAS).ge.1.and.
     &               NOBPT(LGAS).ge.1).or.
     &               (IGAS.ne.JGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL IZERO(IHX_LOC,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS)+1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = 1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(LUWRT,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  a | b  a+)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = -1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = -1.D0 * ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
* operator:
*. a+i beta a+j alpha a k alpha a l alpha, k .ge. l
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS 
                if ((KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &               MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &               NOBPT(KGAS).ge.2.and.NOBPT(IGAS).ge.1.and.
     &               NOBPT(JGAS).ge.1).or.
     &               (KGAS.ne.LGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL IZERO(IHX_LOC,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS)+1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(LUWRT,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  a | a  a+)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = 1 
                      IOP_REO(1,NHX) = 3
                      IOP_REO(2,NHX) = 1
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
* operator:
*. a+i beta  a+j beta  a k beta  a l alpha , i .ge. j
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS 
                if ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &               MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &               NOBPT(IGAS).ge.2.and.NOBPT(KGAS).ge.1.and.
     &               NOBPT(LGAS).ge.1).or.
     &               (IGAS.ne.JGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL IZERO(IHX_LOC,4*NGAS)
                  IHX_LOC(IGAS+(2-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS)+1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = 1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(LUWRT,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  a | b  b+)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 4
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 2
C  ???                SIGN_NHX(NHX) = -1 * ONE
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = 1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
* operator:
*. a+i alpha  a+j beta  a k beta  a l beta , k .ge. l
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS 
                if ((KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &               MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &               NOBPT(KGAS).ge.2.and.NOBPT(IGAS).ge.1.and.
     &               NOBPT(JGAS).ge.1).or.
     &               (KGAS.ne.LGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL IZERO(IHX_LOC,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS)+1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(LUWRT,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  b | +b  b)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE  !!??
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = -1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
      end if
*     end if (IHTYPE.eq.3.or.IHTYPE.eq.5)
*
**********************************************************
*  Spin-dependent spinor basis. E.g. DIRAC environment   *
*      (bu|bu), (ub|ub) integral classes                 *
**********************************************************
*
      if (IHTYPE.eq.4) then
*
* operator:
*. a+i alpha  a+j alpha  a k beta  a l beta, i.ge.j, k.ge.l
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1))) then
                  CALL IZERO(IHX_LOC,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS)+1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS)+1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(LUWRT,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+a  b | +a  b )
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = 1 * ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = -1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 3
                      IOP_REO(4,NHX) = 2
                      SIGN_OPREO(NHX) = -ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
* operator:
*. a+i beta  a+j beta  a k alpha  a l alpha, i.ge.j, k.ge.l
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1))) then
                  CALL IZERO(IHX_LOC,4*NGAS)
                  IHX_LOC(IGAS+(2-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS)+1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS)+1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(LUWRT,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  (+b  a | +b  a)
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 2
                      IHINDX(4,NHX) = 4
                      SIGN_NHX(NHX) = 1 * ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = 1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 3
                      IOP_REO(4,NHX) = 2
                      SIGN_OPREO(NHX) = -ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
      end if
*     end if (IHTYPE .eq. 4)
*
      end if
*     end if (NO_TTWO .gt. 0)
*
#if defined LUCI_DEBUG
      if (IFLAG.ne.1.and.NTEST.ge.50) then
C     if(IFLAG.ne.1)then
        WRITE(LUWRT,'(/A/)') '  *** Information from  GET_HX_RELA ***'
        WRITE(LUWRT,'(/A,I10)') '  Number of excitation types ', NHX
        do I=1,NHX
          write(LUWRT,*) '  excitation type and sign SIGN_NHX ',
     &                     I,SIGN_NHX(I)
        end do
      end if
#endif
*
      if (IFLAG.ne.1) then
#if defined LUCI_DEBUG
        if (NTEST.ge.100) then
          write(LUWRT,'(/A)') '  Max. no of electrons per GAS:'
          call iwrtma(MXNELGS,1,NGAS,1,NGAS)
          WRITE(LUWRT,*) ' The excitation types '
          CALL WRT_SPOX_TP(IHX,NHX)
C              WRT_SPOX_TP(IEX_TP,NEX_TP)
        end if
#endif
      end if
*
      RETURN
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GNSIDE_REL(ISD,IRIS,IRIC,ICA,ICB,IAA,IAB,
     &           IAOC,IBOC,JAOC,JBOC,
     &           NIA,NIB,NJA,NJB,
     &           TB,T_BUFF,T_SCR_CC,IXTP,
     &           SB,CB,ISSM,ICSM,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           TBSUB,CJRES,SIRES,MAXLB,
     &           ICA_STR,ICB_STR,IAA_STR,IAB_STR,
     &           KJA,XKJA,KJB,XKJB,KIA,XKIA,KIB,XKIB,
     &           KJAD,XKJAD,KJBD,XKJBD,KIAD,XKIAD,KIBD,XKIBD,
     &           ITDIAG,IC_RESTRICT,IS_RESTRICT,
     &           SIGN_EXTERNAL,IPRS,IPRD,ISMDST
#if defined (VAR_MPI2)
     &          ,IOFF_BLK_C,IOFF_BLK_S,IACT_C,IACT_S
#endif
     &           )

*
* Sigma routine and density routine for general operator 
*
* sum(Ica,Icb,Iaa,Iab) T(Ica,Icb,Iaa,Iab) Ica+ Icb+ Iaa Iab
*
* contribution from C vector with occupation JAOC,IBOC
* to Sigma vector with occupation IAOC, IBOC
*
* All symmetryblocks of C and Sigma are treated 
*
* ITDIAG = 1 => Diagonal T-block 
*               Only symmetry blocks with 
*                  (ICASM,IAASM).GE.(ICBSM,IABSM)
*               are included, and blocks with 
*                  (ICASM,IAASM) = (ICBSM,IABSM)
*               are given in diagonal form 
*               (Iaa.ge.Iab, if Iaa=Iab, Ica .ge. Icb)
* ITABTRSNP = 1 => alpha-beta transposed block, input block is 
*               T(ICB,ICA,IAB,IAA)
* IX_RESTRICT = 1: BLocks in X are restricted to those with IASM.GE.IBSM
* IX_RESTRICT =-1: BLocks in X are restricted to those with IASM.LT.IBSM
*
* SIGN_EXTERNAL : A sign multiplying all contributions
*                 ( And he said : Mighty Lord, Give me a sign)
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
*
* Jeppe Olsen, August 1999    
* Timo Fleig, June 2001, modified for relativistic case
*             cont´d, 11th of July, Mallorca, around noon (I´m crazy)
*                     15th of July
*                     August 2001
*                     August 2003, ISMDST 
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
*. General input
#include "mxpdim.inc"
#include "gasstr.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "symm.inc"
#include "cgas.inc"
#include "ctcc.inc"
#include "parluci.h"
#include "krmc_shmem.h"
*     
*.Input
      DIMENSION CB(*)
      INTEGER NIA(*),NIB(*),NJA(*),NJB(*)
      INTEGER IAOC(*),IBOC(*),JAOC(*),JBOC(*)
      INTEGER ICA(*),ICB(*),IAA(*),IAB(*)
#if defined (VAR_MPI2)
      INTEGER IACT_C(*), IACT_S(*)
#endif
*.Input or Output
      DIMENSION SB(*)
      dimension TB(*), T_BUFF(*)
*.Local Scratch
      INTEGER KAOC(MXPNGAS),KBOC(MXPNGAS)
      INTEGER NIKAINTM(MXPOBS,MXPLCCOP),NIKBINTM(MXPOBS,MXPLCCOP)
      INTEGER NJKAINTM(MXPOBS,MXPLCCOP),NJKBINTM(MXPOBS,MXPLCCOP)
      INTEGER NKA(MXPOBS),NKB(MXPOBS)
*. 
      INTEGER NICA(MXPOBS),NICB(MXPOBS),NIAA(MXPOBS),NIAB(MXPOBS)
      INTEGER ICA_EXP(MXNOP),ICB_EXP(MXNOP),
     &        IAA_EXP(MXNOP),IAB_EXP(MXNOP)
      INTEGER ICAGP(MXPNGAS),ICBGP(MXPNGAS)
      INTEGER IAAGP(MXPNGAS),IABGP(MXPNGAS)
      integer KAGP(MXPNGAS),KBGP(MXPNGAS)
      integer JAGP(MXPNGAS),JBGP(MXPNGAS)
      INTEGER IAC_AR(MXNOP)
      INTEGER IB_CA(MXNDGIRR,MXNDGIRR,MXPLCCOP),
     &        IB_CB(MXNDGIRR,MXNDGIRR,MXPLCCOP)
      INTEGER IB_AA(MXNDGIRR,MXNDGIRR,MXPLCCOP),
     &        IB_AB(MXNDGIRR,MXNDGIRR,MXPLCCOP)
      INTEGER IB_T(MXNDGIRR,MXNDGIRR,MXNDGIRR)
      INTEGER IB_C(MXNDGIRR),IB_S(MXNDGIRR)
*. Scratch through input 
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION TBSUB(*), T_SCR_CC(*)
      DIMENSION ISMDST(*)
      INTEGER KJA(MAXLB),KJB(MAXLB),KIA(MAXLB),KIB(MAXLB)
      DIMENSION XKJA(MAXLB),XKJB(MAXLB),XKIA(MAXLB),XKIB(MAXLB)
      INTEGER KJAD(MAXLB),KJBD(MAXLB),KIAD(MAXLB),KIBD(MAXLB)
      DIMENSION XKJAD(MAXLB),XKJBD(MAXLB),XKIAD(MAXLB),XKIBD(MAXLB)
*. ... Must hold excitations for all intermediate strings of given sym 
*. and all orbitals of given type
      INTEGER ICA_STR(*),ICB_STR(*),IAA_STR(*),IAB_STR(*)
*     ... Must hold strings of given sym for T ops.
      DIMENSION CJRES(*),SIRES(*)
*
      COMMON/KKKDUMMY/LEN_C2,LEN_S2
      COMMON/MATMLST/XNFLOP,XNCALL,XLCROW,XLCCOL,XLCROWCOL,TMULT
      COMMON/CMXCJ/MXCJ
      INTEGER*8 ITOFF
*
      NTESTL = 0000 
      NTEST = MAX(NTESTL,max(IPRS,IPRD))
csk   NTEST = 100
*
      ISKIP = 0
#if defined (VAR_MPI2)
      LEN_C = LEN_TT_BLOCK_PAR(ICSM,NJA,NJB,NSMST,IACT_C,IOFF_BLK_C,0)
#else
      LEN_C = LEN_TT_BLOCK(ICSM,NJA,NJB,NSMST)
#endif
C     check norm
      XCNORM = DDOT(LEN_C,CB,1,CB,1)
      IF(XCNORM.EQ.0.0D0) GOTO 9999
*
CSK   if (NTEST.ge.100) write(LUWRT,*) 'Coupling with type ',IXTP
*
CSK   IF(NTEST.GE.500) THEN
C       WRITE(LUWRT,*) ' ================ '
C       WRITE(LUWRT,*) ' GNSIDE speaking '
C       WRITE(LUWRT,*) ' ================ '
C       write(LUWRT,*) ' Real / imaginary contribution :'
C       write(LUWRT,*) ' Left , Right ',IRIS,IRIC
C       WRITE(LUWRT,*)
CSK     WRITE(LUWRT,*) ' IAOC and IBOC in GENSIG'
CSK     CALL IWRTMAMN(IAOC,1,NGAS,1,NGAS,LUWRT)
CSK     CALL IWRTMAMN(IBOC,1,NGAS,1,NGAS,LUWRT)
CSK     WRITE(LUWRT,*)
CSK     WRITE(LUWRT,*) ' JAOC and JBOC in GENSIG'
CSK     CALL IWRTMAMN(JAOC,1,NGAS,1,NGAS,LUWRT)
CSK     CALL IWRTMAMN(JBOC,1,NGAS,1,NGAS,LUWRT)
CSK     WRITE(LUWRT,*)
CSK     WRITE(LUWRT,*) ' Type of operator '
CSK     CALL WRT_TP_GENOP(ICA,ICB,IAA,IAB,NGAS)
CSK     IF(ISD.EQ.1) THEN
CSK       WRITE(LUWRT,*) ' Sigma generation '
CSK     ELSE 
CSK       WRITE(LUWRT,*) ' Density generation '
CSK     END IF
*
CSK   WRITE(LUWRT,*) ' Input C vector '
#if defined (VAR_MPI2)
CSK   CALL WRTVH1_PAR(CB,ICSM,NJA,NJB,NSMST,0,IACT_C,IOFF_BLK_C,0)
#else
CSK   CALL WRTVH1(CB,ICSM,NJA,NJB,NSMST,0)
#endif
CSK   if (ISD.eq.2) then
CSK      WRITE(LUWRT,*) ' Input S vector '
#if defined (VAR_MPI2)
CSK      CALL WRTVH1_PAR(SB,ISSM,NIA,NIB,NSMST,0,IACT_S,IOFF_BLK_S,1)
#else
CSK      CALL WRTVH1(SB,ISSM,NIA,NIB,NSMST,0)
#endif
CSK    end if
CSK   END IF
C
C     symmetry of operator (Hamiltonian)
C     ----------------------------------
      IOPSM = IDBGMULT(INVELM(ICSM),ISSM)
C
C     total length of sigma block
#if defined (VAR_MPI2)
      LEN_S = LEN_TT_BLOCK_PAR(ISSM,NIA,NIB,NSMST,IACT_S,IOFF_BLK_S,1)
#else
      LEN_S = LEN_TT_BLOCK(ISSM,NIA,NIB,NSMST)
#endif
C     check S norm (only density runs)
      IF(ISD.EQ.2) THEN
       XSNORM = DDOT(LEN_S,SB,1,SB,1)
       IF(XSNORM.EQ.0.0D0) GOTO 9999
      END IF
#if defined (VAR_MPI2)
C
C     vector partitioning has been changed! we have to be careful
C     calculating the proper offset for each block ( sigma and c )!!!
C
C     adapted offset of c and sigma routines
      CALL IOFF_SYMBLK_MAT_PAR(NSMST,NJA,NJB,ICSM,IB_C,IC_RESTRICT,
     &                         IACT_C,IOFF_BLK_C,0)
      CALL IOFF_SYMBLK_MAT_PAR(NSMST,NIA,NIB,ISSM,IB_S,IS_RESTRICT,
     &                         IACT_S,IOFF_BLK_S,1)
C
#else
C     offset of c and sigma routines
      CALL IOFF_SYMBLK_MAT(NSMST,NJA,NJB,ICSM,IB_C,IC_RESTRICT)
      CALL IOFF_SYMBLK_MAT(NSMST,NIA,NIB,ISSM,IB_S,IS_RESTRICT)
#endif
C     types of strings in T in groupnotation
      CALL OCC_TO_GRP(ICA,ICAGP,1)
      CALL OCC_TO_GRP(ICB,ICBGP,1)
      CALL OCC_TO_GRP(IAA,IAAGP,1)
      CALL OCC_TO_GRP(IAB,IABGP,1)
C          OCC_TO_GRP(IOCC,IGRP,IWAY)
C     operators in T in expanded form
      CALL REF_OP(ICA,ICA_EXP,NCA_OP,NGAS,1,MXNOP)
      CALL REF_OP(ICB,ICB_EXP,NCB_OP,NGAS,1,MXNOP)
      CALL REF_OP(IAA,IAA_EXP,NAA_OP,NGAS,1,MXNOP)
      CALL REF_OP(IAB,IAB_EXP,NAB_OP,NGAS,1,MXNOP)
C          REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY)
C     sign for bringing operators and amplitudes into new order
      NPERM =  NAA_OP*NCB_OP 
     &       + NAA_OP*(NAA_OP-1)/2
     &       + NAB_OP*(NAB_OP-1)/2
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGNXXX  =  -1.0D0*SIGN_EXTERNAL
C       WRITE(LUWRT,*) ' new SIGNXXX for type IXTP =',IXTP
      ELSE
C       WRITE(LUWRT,*) ' old SIGNXXX for type IXTP =',IXTP
        SIGNXXX = SIGN_EXTERNAL
      END IF
C     type of Ka and Kb
      CALL CCEX_OCC_OCC(JAOC,KAOC,NGAS,1,IAA,IKA_ZERO)
      CALL CCEX_OCC_OCC(JBOC,KBOC,NGAS,1,IAB,IKB_ZERO)
C     check that annihilation has not acted upon empty GAS:
      IF(IKA_ZERO.EQ.0.AND.IKB_ZERO.EQ.0) THEN
C      get group of K mapped occupations
       call occ_to_grp(KAOC,KAGP,1)
       call occ_to_grp(KBOC,KBGP,1)
C      number of strings in this supergroup (new DBG routine)
       do KSM = 1,NSMST,1
         call nst_spgrp_dbg(NGAS,KAGP,KSM,NSMST,NKA(KSM),NKADST,1)
         call nst_spgrp_dbg(NGAS,KBGP,KSM,NSMST,NKB(KSM),NKBDST,2)
       end do
C      number of strings in excitation ops
       do ISYM = 1,NSMST,1
         call nst_spgrp_dbg(NGAS,ICAGP,ISYM,NSMST,NICA(ISYM),NCADST,1)
         call nst_spgrp_dbg(NGAS,ICBGP,ISYM,NSMST,NICB(ISYM),NCBDST,2)
C        all arrays refer to strings of creation operators of given occ.
C        so do not bother with difference between anni and crea here.
         call nst_spgrp_dbg(NGAS,IAAGP,ISYM,NSMST,NIAA(ISYM),NAADST,1)
         call nst_spgrp_dbg(NGAS,IABGP,ISYM,NSMST,NIAB(ISYM),NABDST,2)
       end do
C      offset for symmetryblocks of T
       CALL Z_TCC_OFF(IB_T,NICA,NICB,NIAA,NIAB,IOPSM,NSMST,ITDIAG)
C           Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST)
C      creation path from C/S strings to K strings
C      mappings from !I> to !K> are annihilations mappings 
       NMXOP_T = MAX(NCA_OP,NCB_OP,NAA_OP,NAB_OP)
       DO JOP = 1, NMXOP_T
         IAC_AR(JOP) = 1
       END DO
C      CA : KA => IA mappings in 1
C           MAP_EXSTR(IOP,IAC,NOP,IREFOC,IX,SX,NK,IB,SCLFAC,ISMDST,IAB)
       IAB_MAP = 1
       CALL MAP_EXSTR(ICA_EXP,IAC_AR,NCA_OP,IAOC,I1,XI1S,
     &                NIKAINTM,IB_CA,1.0D0,ISMDST,IAB_MAP)
C      CB : KB => IB     
       IAB_MAP = 2
       CALL MAP_EXSTR(ICB_EXP,IAC_AR,NCB_OP,IBOC,I2,XI2S,
     &                NIKBINTM,IB_CB,1.0D0,ISMDST,IAB_MAP)
C      AA : KA => JA
       IAB_MAP = 1
       CALL MAP_EXSTR(IAA_EXP,IAC_AR,NAA_OP,JAOC,I3,XI3S,
     &                NJKAINTM,IB_AA,1.0D0,ISMDST,IAB_MAP)
C      AB : KB => JB
       IAB_MAP = 2
       CALL MAP_EXSTR(IAB_EXP,IAC_AR,NAB_OP,JBOC,I4,XI4S,
     &                NJKBINTM,IB_AB,1.0D0,ISMDST,IAB_MAP)
C      loop over symmetry-blocks of C and Sigma
#if defined (VAR_MPI2)
C      initialize offsets for sigma- and c-block loops
       I_BLK_S_SCR = 0
       I_BLK_C_SCR = 0
       I_BLK_S_SCR = IOFF_BLK_S
       I_BLK_C_SCR = IOFF_BLK_C
#endif
       DO IASM = 1, NSMST 
C        parallel calculation: - sigma-blocks distributed
C                              - check of 'active' c-blocks enabled
#if defined (VAR_MPI2)
        IF( IACT_S( I_BLK_S_SCR ) .ne. MYPROC ) GOTO 3003
#endif
        IBSM = IDBGMULT(ISSM,INVELM(IASM))
#if defined (VAR_MPI2)
C
C        reset offset for c-block loop
         I_BLK_C_SCR = 0
         I_BLK_C_SCR = IOFF_BLK_C
#endif
        DO JASM = 1, NSMST
#if defined (VAR_MPI2)
        IF( IACT_C( I_BLK_C_SCR ) .eq. 0 ) GOTO 2002
#endif
         JBSM = IDBGMULT(ICSM,INVELM(JASM))
CSK      IF(NTEST.GE.1000) WRITE(LUWRT,'(A,4I4)') 
CSK  &   ' IASM, IBSM, JASM, JBSM ', IASM,IBSM,JASM,JBSM
         DO KASM = 1, NSMST
C         IAA_SM is symmetry of ANNIHILATOR string
          IAA_SM = IDBGMULT(KASM,INVELM(JASM))
          ICA_SM = IDBGMULT(IASM,INVELM(KASM))
C         to obtain the number of annnihilator strings of a given sym
C         we find the number of creator strings of conjugated symmetry
          LIAA_SM = NIAA(INVELM(IAA_SM))
          LICA_SM = NICA(ICA_SM)
          LKA_SM = NKA(KASM)
C
          DO KBSM = 1, NSMST
C
           IAB_SM = IDBGMULT(KBSM,INVELM(JBSM))
           ICB_SM = IDBGMULT(IBSM,INVELM(KBSM))
C
           LIAB_SM = NIAB(INVELM(IAB_SM))
           LICB_SM = NICB(ICB_SM)
C
           LKB_SM = NKB(KBSM)
C
           ICBOFF = IB_C(JASM)
           ISBOFF = IB_S(IASM)
C
           IF(LIAA_SM*LICA_SM*LIAB_SM*LICB_SM*LKA_SM*LKB_SM.NE.0)
     &     THEN
C
CSK          IF(NTEST.GE.700) THEN
CSK            WRITE(LUWRT,*) ' ICA_SM, ICB_SM, IAA_SM,IAB_SM =',
CSK  &                      ICA_SM, ICB_SM, IAA_SM,IAB_SM
CSK          END IF 
*
CSK          IF(NTEST.GE.1000) 
CSK  &       WRITE(LUWRT,'(A,2I4)') ' KASM,KBSM ', KASM, KBSM
CSK          IF(NTEST.GE.1000) 
CSK  &       WRITE(LUWRT,*) ' LKA_SM, LKB_SM', LKA_SM, LKB_SM
CSK          IF(NTEST.GE.1000) 
CSK  &       WRITE(LUWRT,'(A,4I3)')      
CSK  &       'LIAA_SM, LICA_SM,LIAB_SM,LICB_SM',
CSK  &        LIAA_SM, LICA_SM,LIAB_SM,LICB_SM
C
C           the operator strings of ICA, ICB, IAA, IAB 
            IUB = 1
            CALL GETSTR2_TOTSM_SPGP_REL(IUB,ICAGP,NGAS,ICA_SM,NCA_OP,
     &                                  LICA_SM2,ICA_STR,NORBT,0,0,0)
            IUB = 2
            CALL GETSTR2_TOTSM_SPGP_REL(IUB,ICBGP,NGAS,ICB_SM,NCB_OP,
     &                                  LICB_SM2,ICB_STR,NORBT,0,0,0)
            IUB = 1
*           we are generating strings of creation operators so 
            IAA_SM_ADJ = INVELM(IAA_SM)
            CALL GETSTR2_TOTSM_SPGP_REL(IUB,IAAGP,NGAS,IAA_SM_ADJ,
     &                                  NAA_OP,
     &                                  LIAA_SM2,IAA_STR,NORBT,0,0,0)
            IUB = 2
            IAB_SM_ADJ = INVELM(IAB_SM)
            CALL GETSTR2_TOTSM_SPGP_REL(IUB,IABGP,NGAS,IAB_SM_ADJ,
     &                                  NAB_OP,
     &                                  LIAB_SM2,IAB_STR,NORBT,0,0,0)
C          batching                   
           LEN_KA_BAT = MAXLB
           LEN_KB_BAT = MAXLB
C
           IF(LKA_SM.LT.MAXLB) THEN
             LEN_KA_BAT = LKA_SM
             LEN_KB_BAT = MIN(LKB_SM,MAXLB*MAXLB/LKA_SM)
           END IF
C
           IF(LKB_SM.LT.MAXLB) THEN
             LEN_KB_BAT = LKB_SM
             LEN_KA_BAT = MIN(LKA_SM,MAXLB*MAXLB/LKB_SM)
           END IF
C
           NKA_BT = LKA_SM/LEN_KA_BAT
           IF(NKA_BT*LEN_KA_BAT. LT. LKA_SM) NKA_BT = NKA_BT + 1
           NKB_BT = LKB_SM/LEN_KB_BAT
           IF(NKB_BT*LEN_KB_BAT. LT. LKB_SM) NKB_BT = NKB_BT + 1
C
           NCA_BT = LICA_SM/MAXLB
           IF(NCA_BT*MAXLB.LT.LICA_SM) NCA_BT = NCA_BT + 1
           NCB_BT = LICB_SM/MAXLB
           IF(NCB_BT*MAXLB.LT.LICB_SM) NCB_BT = NCB_BT + 1
           NAA_BT = LIAA_SM/MAXLB
           IF(NAA_BT*MAXLB.LT.LIAA_SM) NAA_BT = NAA_BT + 1
           NAB_BT = LIAB_SM/MAXLB
           IF(NAB_BT*MAXLB.LT.LIAB_SM) NAB_BT = NAB_BT + 1
C          loop over batches
           DO IKA_BT = 1, NKA_BT
           DO IKB_BT = 1, NKB_BT
            IKA_B = (IKA_BT-1)*LEN_KA_BAT + 1
            IKA_E = MIN(LKA_SM,IKA_B + LEN_KA_BAT - 1)
            NKA_B = IKA_E - IKA_B + 1
            NKA_B0 = NKA_B
            IKB_B = (IKB_BT-1)*LEN_KB_BAT + 1
            IKB_E = MIN(LKB_SM,IKB_B +  LEN_KB_BAT - 1)
            NKB_B = IKB_E - IKB_B + 1
            NKB_B0 = NKB_B
C           loop over annihilation batches
            DO IAA_BT = 1, NAA_BT
             IAA_B = (IAA_BT-1)*MAXLB + 1
             IAA_E = MIN(LIAA_SM,IAA_B + MAXLB -1)
             NAA_B = IAA_E-IAA_B+1
C            KA => JA mapping
             IUB_KJ = 1
             CALL K_TO_J_TOT(KJAD,XKJAD,KASM,IKA_B,IKA_E,
     &                       IAA_STR,IAA_SM,NAA_B,IAA_B,IAA_E,
     &                       I3,XI3S,IB_AA,NJKAINTM,NAA_OP,IZEROKAJA,
     &                       IUB_KJ)
            IF(IZEROKAJA.EQ.0) THEN
            DO IAB_BT = 1, NAB_BT
C            B(eginning) and E(nd) of each Batch
             IAB_B = (IAB_BT-1)*MAXLB + 1
             IAB_E = MIN(LIAB_SM,IAB_B + MAXLB -1)
             NAB_B = IAB_E-IAB_B+1
C            KB => JB mapping
             IUB_KJ = 2
             CALL K_TO_J_TOT(KJBD,XKJBD,KBSM,IKB_B,IKB_E,
     &                       IAB_STR,IAB_SM,NAB_B,IAB_B,IAB_E,
     &                       I4,XI4S,IB_AB,NJKBINTM,NAB_OP,IZEROKBJB,
     &                       IUB_KJ)
            IF(IZEROKBJB.EQ.0) THEN
C           loop over creation batches 
            DO ICA_BT = 1, NCA_BT
             ICA_B = (ICA_BT-1)*MAXLB + 1
             ICA_E = MIN(LICA_SM,ICA_B + MAXLB -1)
             NCA_B = ICA_E-ICA_B+1
C            KA => IA mapping
             IUB_KJ = 1
             CALL K_TO_J_TOT(KIAD,XKIAD,KASM,IKA_B,IKA_E,
     &                       ICA_STR,ICA_SM,NCA_B,ICA_B,ICA_E,
     &                       I1,XI1S,IB_CA,NIKAINTM,NCA_OP,IZEROKAIA,
     &                       IUB_KJ)
            DO ICB_BT = 1, NCB_BT
             ICB_B = (ICB_BT-1)*MAXLB + 1
             ICB_E = MIN(LICB_SM,ICB_B + MAXLB -1)
             NCB_B = ICB_E-ICB_B+1
C            KB => IB mapping
             IUB_KJ = 2
             CALL K_TO_J_TOT(KIBD,XKIBD,KBSM,IKB_B,IKB_E,
     &                       ICB_STR,ICB_SM,NCB_B,ICB_B,ICB_E,
     &                       I2,XI2S,IB_CB,NIKBINTM,NCB_OP,IZEROKBIB,
     &                       IUB_KJ)
             IF(IZEROKBJB.EQ.0.AND.IZEROKAJA.EQ.0.AND.
     &          IZEROKBIB.EQ.0.AND.IZEROKAIA.EQ.0     ) THEN
C            Compress Ka and Kb strings to active subset 
C            COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
C            Ka
             CALL COMPRS2LST_B(KIAD,XKIAD,NCA_B,KJAD,XKJAD,NAA_B,
     &                         NKA_B0,NKA_B,KIA,XKIA,KJA,XKJA)
C            Kb
             CALL COMPRS2LST_B(KIBD,XKIBD,NCB_B,KJBD,XKJBD,NAB_B,
     &                         NKB_B0,NKB_B,KIB,XKIB,KJB,XKJB)
             IM_ACTIVE = 1
C            offset and form of T-coefficients
             ITOFF = IB_T(ICA_SM,ICB_SM,INVELM(IAA_SM))
C            WRITE(LUWRT,*) ' output of IB_T for current block',ITOFF
             IF(ITOFF.LT.0) THEN
               ITOFF = -ITOFF
               ITRNSP = 1
             ELSE
               ITRNSP = 0
             END IF
             IF(ITDIAG.EQ.1.AND.ICA_SM.EQ.ICB_SM.AND.
     &          IAA_SM.EQ.IAB_SM) THEN
                IDIAG = 1
             ELSE
                IDIAG = 0
             END IF
C
             IF(ISD.EQ.1) THEN
               if ((IRIS.eq.1.and.IRIC.eq.2).or.
     &             (IRIS.eq.2.and.IRIC.eq.1)) then
                 if( SPLIT_IJKL )then
                   ITOFF = ITOFF + LEN_T_BUFF
                 end if
               end if
               CALL TCC_SUBBLK(TB(ITOFF),T_BUFF(ITOFF),T_SCR_CC(ITOFF),
     &                         TBSUB,1,
     &                         LICA_SM,NCA_B,ICA_B,LICB_SM,NCB_B,ICB_B,
     &                         LIAA_SM,NAA_B,IAA_B,LIAB_SM,NAB_B,IAB_B,
     &                         ITRNSP,IDIAG,NTEST)
C             WRITE(LUWRT,*) ' TBSUB for IRIC and IRIS',IRIC,IRIS
C             CALL WRTMATMN(TBSUB,1,NCA_B*NCB_B*NAA_B*NAB_B,1,
C    &                      NCA_B*NCB_B*NAA_B*NAB_B,LUWRT)
             ELSE
C             zero density subblock
              CALL DZERO(TBSUB,NCA_B*NCB_B*NAA_B*NAB_B)
             END IF
C  
C            C(Ka,Kb,Iaa,Iab) = sum(Ja,Jb) <Ka!Iaa!Ja><Kb!Iab!Jb> C(Ja,Jb)
C            -------------------------------------------------------------
             CALL C_TO_CKK(CB(ICBOFF),CJRES,1,
     &                     NJA(JASM),NJB(JBSM),
     &                     NKA_B,NKB_B,NAA_B,NAB_B,
     &                     KJA,XKJA,KJB,XKJB)
C
             IF(ISD.EQ.2) THEN
               CALL C_TO_CKK(SB(ISBOFF),SIRES,1,
     &                       NIA(IASM),NIB(IBSM),
     &                       NKA_B,NKB_B,NCA_B,NCB_B,
     &                       KIA,XKIA,KIB,XKIB)
             END IF
C
C            ISD == 1:
C            SIRES(Ka,Kb,Ica,Icb) = CJRES(Ka,Kb,Iaa,Iab)*T(Ica,Icb,Iaa,Iab)
C            ----------------------------------------------------------------
C            ISD == 2:
C            T(Ica,Icb,Iaa,Iab)   = SIRES(Ka,Kb,Ica,Icb)*CJRES(Ka,Kb,Ica,Icb)
C            ----------------------------------------------------------------
             LKAB = NKA_B*NKB_B
             LCAB = NCA_B*NCB_B 
             LAAB = NAA_B*NAB_B
             FACTORC = 0.0D0
             FACTORAB = 1.0D0*SIGNXXX
C            statistics
             XNCALL = XNCALL + 1
             T_INI = Second()
C
#if defined (VAR_MPI2)
             xxxsitime = MPI_WTIME()
#endif
             IF(ISD.EQ.1) THEN
C
C              call DGEMM - statistics
               XLCROW = XLCROW + LKAB
               XLCCOL = XLCCOL + LCAB
               XLCROWCOL = XLCROWCOL + LKAB*LCAB
               XNFLOP = XNFLOP + 2*LKAB*LCAB*LAAB
C              input
               LDA = MAX(1,LKAB)
               LDB = MAX(1,LCAB)
C              output
               LDC = MAX(1,LKAB)
C
               CALL DGEMM('N','T',LKAB,LCAB,LAAB,FACTORAB,CJRES,LDA,
     &                    TBSUB,LDB,FACTORC,SIRES,LDC)
             ELSE
C
C              call DGEMM - statistics
               XLCROW = XLCROW + LCAB
               XLCCOL = XLCCOL + LAAB
               XLCROWCOL = XLCROWCOL + LCAB*LAAB
               XNFLOP = XNFLOP + 2*LCAB*LAAB*LKAB
C              input
               LDA = MAX(1,LKAB)
               LDB = MAX(1,LKAB)
C              output
               LDC = MAX(1,LCAB)
C
               CALL DGEMM('T','N',LCAB,LAAB,LKAB,FACTORAB,SIRES,LDA,
     &                    CJRES,LDB,FACTORC,TBSUB,LDC)
C
             END IF
C
             T_END = Second()
             TMULT = TMULT + T_END - T_INI
C
#if defined (VAR_MPI2)
             xcomputesi = xcomputesi - xxxsitime + MPI_WTIME()
#endif
C            ISD == 1:
C            S(Ia,IB) = SI(Ia,Ib) + <Ka!O+_ca!Ia><Kb!O+cb!Ib>SIRES(Ka,Kb,Ica,Icb)
C            ----------------------------------------------------------------
C            ISD == 2:
C            just scatter out to complete density block
C            ----------------------------------------------------------------
             IF(ISD.EQ.1) THEN
               CALL C_TO_CKK(SB(ISBOFF),SIRES,
     &                       2,NIA(IASM),NIB(IBSM),
     &                       NKA_B,NKB_B,NCA_B,NCB_B,
     &                       KIA,XKIA,KIB,XKIB)
             ELSE IF(ITRNSP.EQ.0) THEN
               if ((IRIS.eq.1.and.IRIC.eq.2).or.
     &             (IRIS.eq.2.and.IRIC.eq.1)) then
                 if( SPLIT_IJKL )then
                   ITOFF = ITOFF + LEN_T_BUFF
                 end if
               end if
               CALL TCC_SUBBLK(TB(ITOFF),T_BUFF(ITOFF),T_SCR_CC(ITOFF),
     &                         TBSUB,2,
     &                         LICA_SM,NCA_B,ICA_B,LICB_SM,NCB_B,ICB_B,
     &                         LIAA_SM,NAA_B,IAA_B,LIAB_SM,NAB_B,IAB_B,
     &                         ITRNSP,IDIAG,NTEST)
             END IF
            END IF
*           ^ End if all maps are nonvanishings
            END DO
            END DO
*           ^ End of loop over batches of creation   operators
            END IF
*           ^ End if batch of beta annihilation was nontrivial
            END DO
            END IF
*           ^ End if batch of alpha annihilation was nontrivial
            END DO
*           ^ End of loop over batches of annihilation operators
           END DO
           END DO
*          ^ End of loop over batches of Ka,Kb
           END IF
*          ^ End if symmetry combinations have nonvanishing dimensions
         END DO
         END DO
*        ^ End of loop over KASM, KBSM
#if defined (VAR_MPI2)
 2002  CONTINUE
*      ^ skip c-block not included in active-block-list; only in parallel
*        calculations
       I_BLK_C_SCR = I_BLK_C_SCR + 1
#endif
        END DO
*       ^ End of loop over JASM
#if defined (VAR_MPI2)
 3003  CONTINUE
*      ^ skip s-block not included in node-list; only in parallel
*        calculations
       I_BLK_S_SCR = I_BLK_S_SCR + 1
#endif
       END DO
*      ^ End of loop over IASM
      END IF
*     ^ End if Ka and Kb are nontrivial strings
C
C     debug print section
C
CSK   IF(NTEST.GE.300) THEN 
C       IF(ISD.EQ.1) THEN
CSK       WRITE(LUWRT,*) ' Updated sigma block '
#if defined (VAR_MPI2)
CSK       CALL WRTVH1_PAR(SB,ISSM,NIA,NIB,NSMST,0,IACT_S,IOFF_BLK_S,1)
#else
CSK       CALL WRTVH1(SB,ISSM,NIA,NIB,NSMST,0)
#endif
C       ELSE IF (ISD.EQ.2) THEN
C         WRITE(LUWRT,*) ' Updated density block '
C         CALL WRT_TCC_BLK(TB,1,NICA,NICB,NIAA,NIAB,NSMST)
C              WRT_TCC_BLK(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
C       END IF
CSK   END IF
      NTEST = 000
 9999 RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE IDIM_TCC_DBG(ITSOSO_TP,NTSOSO_TP,ISYM,
     &                        MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &                        MSCOMB_CC,MX_TBLK_AS,IFTONE,IFTTWO,
     &                        ILTTWO,IFSFLT1,ILSFLT1,N1ELINT,IPRCC)
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
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "csm.inc"
#include "symm.inc"
#include "integrals_off.inc"
*. Specific input
      INTEGER ITSOSO_TP(4*NGAS,NTSOSO_TP)
*. Local scratch 
      INTEGER ICA_GRP(MXPNGAS), ICB_GRP(MXPNGAS)
      INTEGER IAA_GRP(MXPNGAS), IAB_GRP(MXPNGAS)
*.Output : offset for each type
csk   DIMENSION IBTSOSO_TP(NTSOSO_TP)
      INTEGER*8 LEN_SCR1, LEN_SCR2, LEN_SCR3, LEN_SCR4
      INTEGER*8 LEN_SCR5, LEN_SCR6, LEN_SCR7, LEN_SCR8
      INTEGER*8 LENGTH, LENGTH2, N2ELINT_I8
      INTEGER*8 MX_TBLK, MX_TBLK_AS, ITSS_TP_I8
*
      LEN_SCR1 = 0
      LEN_SCR2 = 0
      LEN_SCR3 = 0
      LEN_SCR4 = 0
      LEN_SCR5 = 0
      LEN_SCR6 = 0
      LEN_SCR7 = 0
      LEN_SCR8 = 0
      ITSS_TP_I8 = 0
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRCC)
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' IDIM_TCC : MSCOMB_CC, ISYM = ', MSCOMB_CC,ISYM
      END IF
#endif
*
      LENGTH = 0
      LENGTH2 = 0
      MX_ST_TSOSO = 0
      MX_ST_TSOSO_BLK = 0
      MX_TBLK = 0
      MX_SBSTR = 0
      MX_TBLK_AS = 0
      N1ELINT = 0
*
      DO ITSS_TP = 1, NTSOSO_TP
#if defined LUCI_DEBUG
       IF(NTEST.GE.100) WRITE(6,*) ' ITSS_TP = ', ITSS_TP 
#endif
*. Occupation to group translation 
       CALL OCC_TO_GRP(ITSOSO_TP(1+0*NGAS,ITSS_TP),ICA_GRP,1)
       CALL OCC_TO_GRP(ITSOSO_TP(1+1*NGAS,ITSS_TP),ICB_GRP,1)
       CALL OCC_TO_GRP(ITSOSO_TP(1+2*NGAS,ITSS_TP),IAA_GRP,1)
       CALL OCC_TO_GRP(ITSOSO_TP(1+3*NGAS,ITSS_TP),IAB_GRP,1)
*
#if defined LUCI_DEBUG
       if (NTEST.ge.50) then
         write(6,*) ' Output from OCC_TO_GRP : '
         write(6,*) ' ICA_GRP : ',(ICA_GRP(I),I=1,NGAS,1)
         write(6,*) ' ICA_GRP : ',(ICB_GRP(I),I=1,NGAS,1)
         write(6,*) ' ICA_GRP : ',(IAA_GRP(I),I=1,NGAS,1)
         write(6,*) ' ICA_GRP : ',(IAB_GRP(I),I=1,NGAS,1)
       end if
#endif
*
       NEL_CA = IELSUM(ITSOSO_TP(1+0*NGAS,ITSS_TP),NGAS)
       NEL_CB = IELSUM(ITSOSO_TP(1+1*NGAS,ITSS_TP),NGAS)
       NEL_AA = IELSUM(ITSOSO_TP(1+2*NGAS,ITSS_TP),NGAS)
       NEL_AB = IELSUM(ITSOSO_TP(1+3*NGAS,ITSS_TP),NGAS)
*
       LEN_SCR8 = LENGTH + 1
       ITSS_TP_I8 = ITSS_TP
       IBTSOSO_TP(ITSS_TP_I8) = LEN_SCR8
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
#if defined LUCI_DEBUG
             IF(NTEST.GE.1000) THEN
               WRITE(6,'(A,4I5)') ' sym of CA CB AA AB',
     &         I_CR_AL_SM, I_CR_BE_SM, I_AN_AL_SM, I_AN_BE_SM
               WRITE(6,'(A,4I5)') ' LEN_CA, LEN_CB, LEN_AA, LEN_AB ',
     &                      LEN_CA, LEN_CB, LEN_AA, LEN_AB
             END IF
#endif
*
             LEN_SCR3 = LEN_CA
             LEN_SCR4 = LEN_CB
             LEN_SCR5 = LEN_AA
             LEN_SCR6 = LEN_AB
C
             LEN_SCR1 = LEN_SCR3 * LEN_SCR4 
             LEN_SCR2 = LEN_SCR5 * LEN_SCR6
             LEN_SCR7 = LEN_SCR1 * LEN_SCR2
             LENGTH  = LENGTH  + LEN_SCR7
             LENGTH2 = LENGTH2 + LEN_SCR7
C            count one-electron integrals
             IF((ITSS_TP.lt.IFTTWO).or.(ITSS_TP.ge.IFSFLT1
     &           .and.ITSS_TP.le.ILSFLT1)) N1ELINT = N1ELINT + LEN_SCR7
#if defined LUCI_DEBUG
             if (NTEST.ge.2000) then
               write(6,*) 'LENGTH ', LEN_SCR1 * LEN_SCR2
             end if
#endif
*
             MX_ST_TSOSO = 
     &       MAX(MX_ST_TSOSO,LEN_CA,LEN_CB,LEN_AA,LEN_AB)
             MX_ST_TSOSO_BLK = 
     &       MAX(MX_ST_TSOSO_BLK,LEN_CA*NEL_CA,LEN_CB*NEL_CB,
     &                           LEN_AA*NEL_CA,LEN_AB*NEL_CB)
             MX_TBLK = MAX(MX_TBLK, LEN_SCR7)
             MX_SBSTR = MAX(MX_SBSTR,LEN_CA,LEN_CB,LEN_AA,LEN_AB)
          END DO
        END DO
       END DO
      MX_TBLK_AS = MAX( MX_TBLK_AS, LENGTH2 )
      END DO
C
      I_TOT_DIM_T = LENGTH
C
C     determine number of one-electron and two-electron integrals
C
CSK   N1ELINT = IBTSOSO_TP(IFTTWO)  - IBTSOSO_TP(IFTONE)
C
      NINT_2E    = LENGTH - N1ELINT
      
      N2ELINT_I8 = LENGTH - N1ELINT
*
CSK   NTEST = 5
#if defined LUCI_DEBUG
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' Number of 1e-integrals : ', N1ELINT
        WRITE(6,*) ' Number of 2e-integrals : ', NINT_2E
        WRITE(6,*) ' total number of integrals: ', NINT_2E + N1ELINT
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
        WRITE(6,*) ' Offset for each block '
         DO II = 1, NTSOSO_TP
           WRITE(6,*) 'IBTSOSO_TP(',II,') is ',IBTSOSO_TP(II)
         END DO
      END IF
#endif
      NTEST = 0000
*
      RETURN
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION ITDIANUM(ICA,ICB,IAA,IAB,NC,NA)
*
* Find adress of element T(ICA,ICB,IAA,IAB) in diagonal 
* block
* It is assumed that (ICA,IAA) .ge. (ICB,IAB) has been checked outside
*
* Jeppe Olsen, July 2000, HNIE
*
#include "implicit.inc"
*. T(1,1,Iab,Iab)
      I11IABIAB = ((IAB-1)*NA - IAB*(IAB-1)/2)*NC*NC+
     &             (IAB-1)*NC*(NC+1)/2
*. T(1,1,Iaa,Iab)
      IF(IAA.GT.IAB) THEN
       I11IAAIAB =  I11IABIAB + NC*(NC+1)/2 + (IAA-IAB-1)*NC*NC
      ELSE
       I11IAAIAB =  I11IABIAB 
      END IF
*. T(Ica,Icb,Iaa,Iab)
      IF(IAA.GT.IAB) THEN
       IADR = I11IAAIAB  + (ICB-1)*NC + ICA
      ELSE
       IADR = I11IAAIAB  + (ICB-1)*NC + ICA - ICB*(ICB-1)/2
      END IF
*
      ITDIANUM = IADR
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Address of diagonal T-element : '
        WRITE(6,*) ' ICA, ICB, IAA, IAB, NC,NA, IADR = ',
     &               ICA, ICB, IAA, IAB, NC,NA, IADR
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE LEN_GENOP_STR_MAP(NGENOP,IGENOP,NSPGRP,ISPGRP,NOBPT,
     &                             NGAS,MAXLEN)
*
* A set of general operators, IGENOP, and a set of supergroups, ISPGRP, 
* are defined. Obtain max.length of mappings from resolution strings 
* to istrings. 
*
* The resolution strings KSTR are assumed to be inserted between the 
* creation and annihilation operators 
*
*   <ISPGRP!ICREA!KSTR> <KSTR!IANNI!ISPGRP'>
* So the mappings from KSTR to ISTR are always creation mappings 
*
* Jeppe Olsen, July 2000 ( At summerschool, HNIE)
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
*
*. Specific input
      INTEGER IGENOP(2*NGAS,*), ISPGRP(MXPNGAS,*), NOBPT(NGAS)
*. Local scratch
      INTEGER KSPGRP(MXPNGAS), IGENOP_EXP(MXPLCCOP)
*
      MAXLEN = 0
      IONE = 1
      IMONE = -1
      DO ICA = 1, 2
C?     print*,'ICA ',ICA
      DO JSPGRP = 1, NSPGRP
C?     print*,'JSPGRP ',JSPGRP
*. ICREA/IANNI(dag)*!KSTR> = !ISTR>
       DO JGENOP = 1, NGENOP
C?      print*,'JGENOP ',JGENOP
*. Occupation of KSPGRP
         IF(ICA.EQ.1) THEN
           CALL IVCSUM(KSPGRP,ISPGRP(1,JSPGRP),IGENOP(1,JGENOP),
     &                 IONE,IMONE,NGAS)
         ELSE
           CALL IVCSUM(KSPGRP,ISPGRP(1,JSPGRP),IGENOP(NGAS+1,JGENOP),
     &                 IONE,IMONE,NGAS)
         END IF
*. Is KSPGRP a correct supergroup( all occ larger than zero) 
         I_AM_OKAY = 1
         DO IGAS = 1, NGAS
           IF(KSPGRP(IGAS).LT.0) I_AM_OKAY = 0
         END DO
         IF(I_AM_OKAY.EQ.1) THEN
*. IGENOP in expanded form
C  REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY)
          LEN = 0
          IF(ICA.EQ.1) THEN
            CALL REF_OP(IGENOP(1,JGENOP),IGENOP_EXP,NOP,
     &                  NGAS,1,MXPLCCOP)
          ELSE
            CALL REF_OP(IGENOP(NGAS+1,JGENOP),IGENOP_EXP,NOP,
     &                  NGAS,1,MXPLCCOP)
          END IF
          DO IOP = 1, NOP
C?         print*,'IOP ',IOP
            IOPTP = IGENOP_EXP(IOP)
*. Number of Kstrings
C                  NST_FOR_OCC(NELEC_PER_GAS,NORB_PER_GAS,NGAS)
            NSTR = NST_FOR_OCC(KSPGRP,NOBPT,NGAS)
            NORB = NOBPT(IOPTP)
            LEN = LEN + NSTR*NORB 
       if(LEN.gt.100000000) print*,'NSTR,NORB ',NSTR,NORB
*. Update KSPGRP 
            KSPGRP(IOPTP) = KSPGRP(IOPTP) + 1
          END DO
          MAXLEN = MAX(MAXLEN,LEN)
         END IF
       END DO
      END DO
      END DO
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Max length of KSTR => ISTR map = ', MAXLEN
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION LEN_TT_BLOCK(ISM,NIA,NIB,NSMST)
*
* Length of TT block with total sym ISM
*
* Jeppe Olsen, May 99 in Aarhus
* Modified for relativistic case, Timo Fleig, June 2001
*
#include "implicit.inc"
#include "mxpdim.inc"
*. Specific input
      INTEGER NIA(*), NIB(*)
*. General input
#include "symm.inc"
*
      LEN = 0
      DO IASM = 1, NSMST
        IBSM = IDBGMULT(INVELM(IASM),ISM)
        LEN = LEN + NIA(IASM)*NIB(IBSM)
      END DO
*
      LEN_TT_BLOCK = LEN
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' LEN_TT_BLOCK = ', LEN_TT_BLOCK 
      END IF
#endif
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION LEN_TT_BLOCK_PAR(ISM,NIA,NIB,NSMST,INCLUDE_LIST,IOFF1,
     &                          IVERS)
*
* Length of TT block with total sym ISM
*
* Jeppe Olsen, May 99 in Aarhus
* Modified for relativistic case, Timo Fleig, June 2001
* Modified for parallel runs, Stefan Knecht, May 2007
*
* IOFF1: absolute address value
* IVERS: 0: c-block version; 1: s-block version
*
#include "implicit.inc"
#include "mxpdim.inc"
*. Specific input
      INTEGER NIA(*), NIB(*), IOFF1, IOFF2, INCLUDE_LIST(*)
*. General input
#include "symm.inc"
#include "parluci.h"
*
      LEN   = 0
      IOFF2 = 0
      IOFF2 = IOFF1
CSK      WRITE(LUWRT,*) ' IOFF1 and IOFF2 are now',IOFF1,IOFF2
CSK      WRITE(LUWRT,*) ' IVERS is ',IVERS
CSK      WRITE(LUWRT,*) ' NSMST is ',NSMST
      
      IF( IVERS .eq. 0 ) THEN
        DO IASM = 1, NSMST
          IF( INCLUDE_LIST(IOFF2) .eq. 0 ) GOTO 100
          IBSM = IDBGMULT(INVELM(IASM),ISM)
          LEN = LEN + NIA(IASM)*NIB(IBSM)
 100      CONTINUE
          IOFF2 = IOFF2 + 1
        END DO
      ELSE
        DO IASM = 1, NSMST
          IF( INCLUDE_LIST(IOFF2) .ne. MYPROC ) GOTO 200
          IBSM = IDBGMULT(INVELM(IASM),ISM)
CSK          WRITE(LUWRT,*) ' IASM, IBSM',IASM,IBSM
          LEN = LEN + NIA(IASM)*NIB(IBSM)
 200      CONTINUE
          IOFF2 = IOFF2 + 1
        END DO
      END IF
*
      LEN_TT_BLOCK_PAR = LEN
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*) ' LEN_TT_BLOCK_PAR = ', LEN_TT_BLOCK_PAR
      END IF
#endif
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MAP_EXSTR(IOP,IAC,NOP,IREFOC,IX,SSX,NK,IBX,SCLFAC,
     &                     ISMDST,IAB)

* Information for mapping |Kstr> = IOP |Irefoc>
*
* Mapping for each creation/annihilation operator is given
*
*     Jeppe Olsen, Summer of 99
*     Timo Fleig, summer of 2001, mod. for rel. case
*                 IAB added for alpha (unbarred) and beta (barred)
*                 ISMDST added Aug. 2003
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "cgas.inc"
#include "csm.inc"
#include "symm.inc"
#include "ctcc.inc"
#include "parluci.h"
*. Input
      INTEGER IOP(NGAS),IAC(NGAS),IREFOC(NGAS)
*. Output
      INTEGER IX(*),IBX(MXNDGIRR,MXNDGIRR,NOP)
      INTEGER NK(MXPOBS,NOP)
      DIMENSION SSX(*)
*. Local Scratch
      INTEGER KOC(MXPNGAS),KGRP(MXPNGAS)
      INTEGER ISMDST(*)
*. IBX(ISM,JSM,JOP) : Start of map in JX for orbitals of symmetry ISM,
*.                    input strings of symmetry JSM for operator JOP
*. Notice at each level N the mapping is from !K(N)> to !K(N-1)> 
*  where !K(0)> is input type and !K(NOP)> is output type
*
C     MAP_EXSTR(IOP,IAC,NOP,IREFOC,IX,SSX,NK,IB,SCLFAC)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' MAP_EXSTR speaking' 
         WRITE(6,*) ' =================='
        WRITE(6,*) ' Reference occupation '
        CALL IWRTMA(IREFOC,1,NGAS,1,NGAS) 
        WRITE(6,*) ' NOP = ', NOP
        WRITE(6,*) ' Operator string, Orbspaces and C/A '
        CALL IWRTMA(IOP,1,NOP,1,NOP)
        CALL IWRTMA(IAC,1,NOP,1,NOP)
      END IF
#endif
*
      CALL ICOPY(NGAS,IREFOC,1,KOC,1)
      CALL OCC_TO_GRP(KOC,KGRP,1)
      IOFF = 1
#if defined LUCI_DEBUG
      IF(NTEST.GE.1000) WRITE(6,*) ' NOP = ', NOP
#endif
      DO JOP = NOP,1,-1
#if defined LUCI_DEBUG
       IF(NTEST.GE.1000) WRITE(6,*) ' JOP = ', JOP
#endif
*
       JGAS = IOP(JOP)
       JAC  = IAC(JOP)
*
       IF(JAC.EQ.1) THEN 
         JAC_ADJ = 2 
       ELSE 
         JAC_ADJ = 1 
         WRITE(6,*)' Map <I|a|K> not debugged for DBG symmetry'
         STOP' Map <I|a|K> not debugged for DBG symmetry'
       END IF
*. 
COLD
CTF Depending on the evaluation "direction" of the coupling coefficient,
CTF the symmetry handling changes: If |K> is on the right, the inverse
CTF element is correct (INVELM(JSMOB)), if |J> is on the right, the
CTF element is NOT inversed!
CTF
       DO JSMOB = 1, NSMOB
        DO JSMSPGP = 1, NSMST
*. Symmmetry of |Kgrp> = a(JSMOB) |JSMSPGP>
*. = symmetry obtained from |JSMSPGP> = a+(JSMOB) |Kgrp>
         KSMSPGP = IDBGMULT(INVELM(JSMOB),JSMSPGP)
*
         IBX(JSMOB,KSMSPGP,JOP) = IOFF
         CALL ADAST_GAS_REL(JSMOB,JGAS,NGAS,KGRP,JSMSPGP,
     &                  IX(IOFF),SSX(IOFF),NKSTR,IEND,
     &                  IFRST,KFRST,KACT,SCLFAC,JAC_ADJ,ISMDST,IAB)
C             ADAST_GAS_REL(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
C    &                  I1,XI1S,NKSTR,IEND,
C    &                  IFRST,KFRST,KACT,SCLFAC,IAC,ISMDST,IAB)
         NK(KSMSPGP,JOP) = NKSTR
         if (IAB.eq.1) then
           IOFF = IOFF + NOBPTS(JGAS,JSMOB)*NKSTR
         else if (IAB.eq.2) then
           IOFF = IOFF + NOBPTS2(JGAS,JSMOB)*NKSTR
         end if
         IF(IOFF-1.GT. MAXLEN_I1) THEN
           WRITE(6,*) ' MAP... MAXLEN_I1 too small '
           WRITE(6,*) ' IOFF and MAXLEN_I1 = ', IOFF, MAXLEN_I1
           STOP       ' MAP... MAXLEN_I1 too small '
         END IF
*
        END DO
       END DO
*. Updated Kgroup  
*. Info on !K(JOP)> 
       IF(JAC.EQ.1) THEN 
         JDELTA = -1 
       ELSE 
         JDELTA = +1
       END IF
       KOC(JGAS) = KOC(JGAS) + JDELTA
       CALL OCC_TO_GRP(KOC,KGRP,1)
#if defined LUCI_DEBUG
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' KOC and KGRP '
         CALL IWRTMA(KOC,1,NGAS,1,NGAS)
         CALL IWRTMA(KGRP,1,NGAS,1,NGAS)
       END IF
#endif
*
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from MAP_EXSTR '
        WRITE(6,*) ' Mappings from a+/a !K(J)> to !K(J-1)> '
*
        DO JOP = 1, NOP
         IF(IAC(JOP).EQ.1) THEN
           WRITE(6,*) 
     &     ' Operator ', JOP, ' is annihilation operator of type',
     &     IOP(NOP) 
         ELSE
           WRITE(6,*) 
     &     ' Operator ', JOP, ' is creation     operator of type',
     &     IOP(NOP) 
         END IF
         DO JSM = 1, NSMOB
           DO KNSM = 1, NSMST
             WRITE(6,*) ' Sym of operator and !K(N)> ',
     &       JSM,KNSM 
             LK = NK(KNSM,JOP)
             WRITE(6,*) ' Number of strings !K(N)> ',LK
             if (IAB.eq.1) then
               LJ = NOBPTS(IOP(JOP),JSM) 
             else if (IAB.eq.2) then
               LJ = NOBPTS2(IOP(JOP),JSM) 
             end if
             DO JJ = 1, LJ
              WRITE(6,*) ' Info for orbital ', JJ
              WRITE(6,*) ' Excited strings and sign '
              IOFF2 = IBX(JSM,KNSM,JOP) + (JJ-1)*LK
              WRITE(6,*) ' JOP JSM KNSM LK ',JOP,JSM,KNSM,LK
              WRITE(6,*) ' IOFF2 ', IOFF2
              CALL IWRTMA(IX(IOFF2),1,LK,1,LK)        
              CALL WRTMAT(SSX(IOFF2),1,LK,1,LK)        
             END DO
           END DO
         END DO
        END DO
      END IF
#endif
*
      RETURN
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Find number of elements in a given excitation class.
*
      integer*8 function nelm_cc(ITP,NXTP,NT_SCRATCH)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "integrals_off.inc"
*
      INTEGER*8 NT_SCRATCH
*
      if (ITP.lt.NXTP) then
        NELM_CC = IBTSOSO_TP(ITP+1) - IBTSOSO_TP(ITP)
      else if (ITP.eq.NXTP) then
        NELM_CC = NT_SCRATCH - IBTSOSO_TP(ITP) + 1
      else
        write(6,*) 'Requested excitation type ',ITP
        write(6,*) 'exceeds total number of types ',NXTP
        stop 'QUITTING in integer*8 function nelm_cc.'
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION NST_FOR_OCC(NELEC_PER_GAS,NORB_PER_GAS,NGAS)
*
* A supergroup is defined by NELEC_PER_GAS. Find the 
* Total number of strings of this supergroup
*
* Jeppe Olsen, July 2000 (HNIE)
*
#include "implicit.inc"
*
      INTEGER NELEC_PER_GAS(*),NORB_PER_GAS(*)
*
      NSTR = 1
      DO IGAS = 1, NGAS
       IF(NELEC_PER_GAS(IGAS).NE.0) THEN
        NSTR = NSTR*IBION(NORB_PER_GAS(IGAS),NELEC_PER_GAS(IGAS))
       END IF
      END DO
*
      NST_FOR_OCC = NSTR
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of electrons and orbitals per gasspace '
        CALL IWRTMA(NELEC_PER_GAS,1,NGAS,1,NGAS)
        CALL IWRTMA(NORB_PER_GAS,1,NGAS,1,NGAS)
        WRITE(6,*) ' Number of strings = ', NSTR
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE OCC_TO_GRP(IOCC,IGRP,IWAY)
*. Translate between occupation and group labels of supergroup
*
* IWAY = 1 => OCC to group
*      = 2 => Group to occ
*
* Number of strings in supergroup with given symmetry and number of 
* electrons in each space
*
* Jeppe Olsen, March 1999
*
#include "implicit.inc"
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
C     INTEGER JGRP(MXPNGAS)
*. Specific input/output
      INTEGER IOCC(NGAS),IGRP(NGAS)
*
      NTEST = 000
*
      IF(IWAY.EQ.1) THEN 
*. Occupation => Group number 
        DO IOBTP = 1, NGAS
          JJGRP = 0
          DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
            IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
          END DO
          IGRP(IOBTP) = JJGRP
*
#if defined LUCI_DEBUG
          IF(NTEST.ge.10) THEN
           WRITE(6,*) ' Group not included in list '
           WRITE(6,*) ' Input occupations : '
           CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
           write(6,*) 'Current orbital type :',IOBTP
           write(6,*) 'IBGPSTR array :'
           call iwrtma(IBGPSTR,1,NGAS,1,NGAS)
           write(6,*) 'NGPSTR array :'
           call iwrtma(NGPSTR,1,NGAS,1,NGAS)
           write(6,*) 'Group range ',
     &                 IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
           write(6,*) 'List of electron occupations: '
           do K=IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
             write(6,*) 'NELFGP(',K,') =', NELFGP(K)
           end do
          END IF
#endif
          if (JJGRP.eq.0) 
     &       CALL QUIT(' *** ERROR in OCC_TO_GRP: Group not included 
     &                   in list ***')
        END DO
*
      ELSE
*. Group => Occupation 
        DO IOBTP = 1, NGAS
          IOCC(IOBTP) = NELFGP(IGRP(IOBTP))
        END DO
      END IF
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Occupation and corresponding group array '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS) 
      END IF
#endif
*
      RETURN
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY,MXNOP)
* An operatorstring may be specifed as
*
* IOPGAS : Number of operators per GASspace
* IOP    : GASpace of each operator
*
* Transform between these two form
*
* Iway = 1 : IOPGAS => IOP
* Iway = 2 : IOP    => IOPGAS
*
* Jeppe Olsen, Summer of 99
*
#include "implicit.inc"
#include "ipoist8.inc"
*. Input/Output
      INTEGER IOPGAS(NGAS),IOP(MXNOP)
*
      IF(IWAY.EQ.1) THEN
        JOP = 0
        DO JGAS = 1, NGAS
          LJGAS = IOPGAS(JGAS)
          DO JJOP = 1, LJGAS
            JOP = JOP + 1
            IOP(JOP) = JGAS
          END DO
        END DO
        NOP = JOP
        if (NOP.gt.MXNOP) then
          write(6,*) 'NOP too large in REF_OP. NOP =',NOP
          write(6,*) 'Increase MXNOP or find bug !!'
        end if
      ELSE
        DO JGAS = 0, NGAS
          JOP = 0
          DO JJOP = 1, NOP
            IF(IOP(JJOP).EQ.JGAS) JOP = JOP +1 
          END DO
          IOPGAS(JGAS) = JOP
        END DO
      END IF
*
#if defined LUCI_DEBUG
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' IOPGAS => IOP '
        ELSE
         WRITE(6,*) ' IOP => IOPGAS '
        END IF
        WRITE(6,*) ' IOPGAS and IOP '
        CALL IWRTMA(IOPGAS,1,NGAS,1,NGAS)
        CALL IWRTMA(IOP,1,NOP,1,NOP)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SET_HOP_DBG(IHTYPE,NO_TTWO,ISYM_T,IPRNT)
*
*     Construct orbital excitation types for different choices
*     of Hamiltonian operator, double group version
*
*
*     IHTYPE = 1 : Include only Delta MK = 0 terms
*     IHTYPE = 2 : Include Delta MK = 0 for two-electron part
*                  and complete one-electron part
*     IHTYPE = 3 : Include all terms for general relativistic Hamiltonian
*
*     NO_TTWO .gt. 0 : do not include two-electron terms
*
*     added INTEGER*8 adaption: S. Knecht - Nov. 2007
*
#include "implicit.inc"
#include "ipoist8.inc"
*  
#include "dgroup.h"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "cprnt.inc"
#include "orbinp.inc"
C
C         OUTPUT
C
#include "ctcc.inc"
#include "integrals_off.inc"
#include "parluci.h"
C
C     specific input in order to handle large integral arrays
C     even on 32-bit systems
      INTEGER*8 MX_TBLK, MX_TBLK_AS, MSCOMB_T
C
C     number of spin-orbital excitation operators in H
C
      call memmar(KLSOBEX_TEMP,4*NGAS,'ADDS  ',1,'SPOBET')
      IFLAG = 1
      CALL GET_HX_RELA(IFLAG,IHTYPE,NO_TTWO,NOBPT,NHX,
     &                 IFTONE,IFTTWO,ILTTWO,IFSFLT1,ILSFLT1,
     &                 WORK(KLSOBEX_TEMP),
     &                 DUMMY,DUMMY,DUMMY,I2DIATP,DUMMY,DUMMY,DUMMY,
     &                 IPRNT)
      NSPOBEX_TP = NHX
C
C     actual spin-orbital excitations
C
      call memmar(KLSOBEX,4*NGAS*NSPOBEX_TP,'ADDS  ',1,'SPOBEX')
      call memmar(KLIHIND,     4*NSPOBEX_TP,'ADDS  ',1,'IHINDX')
      call memmar(KSIGNNHX,      NSPOBEX_TP,'ADDS  ',2,'SIGNDX')
      call memmar(KLABEXTP,    4*NSPOBEX_TP,'ADDS  ',1,'LABXTP')
      call memmar(KLOP_REO,    4*NSPOBEX_TP,'ADDS  ',1,'OP_REO')
      call memmar(KSIOPREO,      NSPOBEX_TP,'ADDS  ',2,'SIOPRE')
      IFLAG = 0
      CALL GET_HX_RELA(IFLAG,IHTYPE,NO_TTWO,NOBPT,NHX,
     &                 IFTONE,IFTTWO,ILTTWO,IFSFLT1,ILSFLT1,
     &                 WORK(KLSOBEX_TEMP),
     &                 WORK(KLSOBEX),WORK(KLIHIND),WORK(KSIGNNHX),
     &                 I2DIATP,WORK(KLABEXTP),WORK(KLOP_REO),
     &                 WORK(KSIOPREO),IPRNT)
C
C     Space for info in the different parts of Hamiltonian
C
C     KLSPOBEX_AC : Is this spinorbital excitation active
C
      call memmar(KLSPOBEX_AC,NSPOBEX_TP,'ADDS  ',1,'SPOBAC')
C
C     all spinorbital excitations are initially active
C
      IONE     = 1
      IZERO_I8 = 0
C
      IF( NSPOBEX_TP .gt. MXINT_TP ) THEN
          WRITE(LUWRT,*) ' STOP in subroutine SET_HOP_DBG !!!'
          WRITE(LUWRT,*) ' Number of SO-integral types       '
          WRITE(LUWRT,*) ' larger than MXINT_TP: ', NSPOBEX_TP, MXINT_TP
          WRITE(LUWRT,*) ' increase MXINT_TP (mxpdim.inc): to',
     &                     NSPOBEX_TP
          stop ' Quitting in subroutine SET_HOP_DBG'
      END IF
C
      CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TP)
      CALL ISETVC8_2(IBTSOSO_TP,IONE,NSPOBEX_TP)
C
C     dimensions for operator of symmetry ISYM_T (set in KRMC_LUCI2)
C
      MSCOMB_T = 0
C
      CALL IDIM_TCC_DBG(WORK(KLSOBEX),NSPOBEX_TP,ISYM_T,MX_ST_TSOSO,
     &                  MX_ST_TSOSO_BLK,MX_TBLK,MSCOMB_T,MX_TBLK_AS,
     &                  IFTONE,IFTTWO,ILTTWO,IFSFLT1,ILSFLT1,N1ELINT,
     &                  IPRNT)
C
C     output -- transfer information to common block
C
      LEN_T_VEC          = I_TOT_DIM_T
      NALLINT            = LEN_T_VEC * MIN(NZ,2)
      N2ELINT            = NINT_2E
      MX_TBLK_MX         = MX_TBLK
      MX_ST_TSOSO_MX     = MX_ST_TSOSO
      MX_ST_TSOSO_BLK_MX = MX_ST_TSOSO_BLK
C
C     print section
C
      IF( IPRNT .gt. 0)THEN
        WRITE(LUWRT,'(/A )')'  *** OUTPUT FROM SET_HOP_DBG *** ' 
        WRITE(LUWRT,'( A/)')'      _______________________     ' 
        WRITE(LUWRT,'(A,I20)')'  symmetry of T op. :',ISYM_T
        WRITE(LUWRT,'(A,I20)')'  LEN_T_VEC         :',LEN_T_VEC 
        WRITE(LUWRT,'(A,I20)')'  N2ELINT           :',N2ELINT
        WRITE(LUWRT,'(A,I20)')'  MX_TBLK_MX        :',MX_TBLK_MX
        WRITE(LUWRT,'(A,I20)')'  MX_ST_TSOSO_MX    :',MX_ST_TSOSO_MX
        WRITE(LUWRT,'(A,I20)')'  MX_ST_TSOSO_BLK_MX:',MX_ST_TSOSO_BLK_MX
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SIGDEN_CTRL(C,HC,LUC,LUHC,T,T_BUFF,ISIGDEN
#if defined (VAR_MPI2)
     &                   ,LUCLIST,LUHCLIST,IBLOCKL_C,NPARBLOCK_C,
     &                    IGROUPLIST,IPROCLIST,RCCTOS,
     &                    IT_TTPL,IT_TTOL
     &                   ,IBLOCKL_S,NPARBLOCK_S
#endif
     &                   )
*
* Outer routine for sigma vector generation
* in SIGDEN formalism
*
* Modified mv7 routine 
*  Timo Fleig, Halloween 2001  =:-()
*
#include "implicit.inc"
* some stuff relevant for parallel runs
#if defined (VAR_MPI2)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)     
      INTEGER RCCTOS(*)
#endif
#include "infpar.h"
#include "parluci.h"
#include "krciprop.h"
*
#include "ipoist8.inc"
#include "mxpdim.inc"
*
* ======
*.Output
* ======
*
#include "wrkspc.inc"
#include "cstate.inc"
#include "cands.inc"
#include "ctcc.inc"
*
      CALL QENTER('SIGDC')
      IDUM = 0
      CALL MEMMAR(KXDUM,IDUM,'MARK  ',IDUM,'SIGDC ')
C
C     arrays for partitioning of sigma + information for C vector 
C     *_MS2_* arrays are in CSTATE
C
C     DOSIGPROP controls output for NBLK_MS2_C; I_SET_L2BLOCK as ICSM 
C     carrier to z_blkfo_rel. note the minus sign!
      ICSM_SCR = 0
      IF( DOSIGPROP ) ICSM_SCR =  - ICSM
      call z_blkfo_rel(IDC,NMS2VAL,1,ISSM,
     &                 KLSLBT,KLSLEBT,KLSI1BT,KLSIBT,KXX_BLTP,
     &                 NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,
     &                 NBAT_MS2,IBAT_MS2,NBLK_MS2_C,DOSIGPROP,ICSM_SCR,
#if defined (VAR_MPI2)
     &                 2,NPARBLOCK_S
#else
     &                 0,IDUMMY
#endif
     &                 )
C
#if defined (VAR_MPI2)
C
C     ISCLFAC_GROUP array for c-blocks: 
C     = 0: block is not active
C     > 0: block is active
C
      CALL MEMMAR(KSCALLOC,NUM_BLOCKS2,'ADDL  ',1,'ICLLOC')
C
#endif
C
      call sigden_rel(C,HC,NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                WORK(KLSI1BT),WORK(KLSIBT),LUC,LUHC,T,T_BUFF,
     &                ISIGDEN
#if defined (VAR_MPI2)
     &               ,LUCLIST,LUHCLIST,IBLOCKL_C,NPARBLOCK_C,
     &                IGROUPLIST,IPROCLIST,RCCTOS,WORK(KSCALLOC),
     &                IT_TTPL,IT_TTOL
     &               ,IBLOCKL_S,NPARBLOCK_S
#endif
     &                 )
C
      CALL MEMMAR(KXDUM,IDUM,'FLUSM ',IDUM,'SIGDC ')
      CALL QEXIT('SIGDC')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SIGDEN_REL(C,HC,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &                      LLUC,LLUHC,T,T_BUFF,ISIGDEN
#if defined (VAR_MPI2)
     &                     ,LUCLIST,LUHCLIST,IBLOCKL_C,NPARBLOCK_C,
     &                      IGROUPLIST,IPROCLIST,RCCTOS,ISCLFAC_GROUP,
     &                      IT_TTPL,IT_TTOL
     &                     ,IBLOCKL_S,NPARBLOCK_S
#endif
     &                 )
C***********************************************************************
C
C     outer routine for sigma (ISIGDEN==1) and density (ISIGDEN==2) 
C     calculation for a general relativistic CI/CC code
C
C     for sigma T is the input set of integrals (CI) / cluster
C     amplitutes (CC)
C     for density T is the output set of density matrix elements
C
C     Jeppe Olsen, Summer of 99, March 2001
C 
C     Modifications for vector storage and sigma loops
C     Timo Fleig, Halloween 2001
C
C     revision and (parallel) generalization: Stefan Knecht - Oct 2008
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
#include "cands.inc"
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
#include "ctcc.inc"
#include "clunit.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      DIMENSION LUCLIST(*), IBLOCKL_C(*), NPARBLOCK_C(*), IGROUPLIST(*)
      DIMENSION IPROCLIST(*), IBLOCKL_S(*), NPARBLOCK_S(*)
      INTEGER   RCCTOS(*)
      DIMENSION ISCLFAC_GROUP(*), LUHCLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR1, IOFF_SCR2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR1_C, IOFF_SCR2_C
      CHARACTER SECTID*12,WALLTID*12,WALLTID2*12,WALLTID3*12
#endif /* if def VAR_MPI2 */
#include "parluci.h"
#include "krciprop.h"
#include "files_r.inc"
*
#include "dgroup.h"
      LOGICAL SBLK_READS2
*
*     batches of sigma
      DIMENSION LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
      DIMENSION HC(*),C(*)
*     scratch
      INTEGER IONE, INT_IOFF1, INT_IOFF2
      INTEGER INT_IOFF1_C, INT_IOFF2_C
*
      IONE        = 1
      INT_IOFF1   = 0
      INT_IOFF2   = 0
      ISI_CALC_BL = 0
      IBI_MULT_BL = 0
      xreadtimebi = 0.0D0
      xcomputesi  = 0.0D0
#if defined (VAR_MPI2)
      IOFF_SCR1 = 0
      IOFF_SCR2 = 0
      mem_cp_time = 0.0D0
#endif
*
      CALL QENTER('SIGDE')
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'SIGDEN')
*
      NTESTL = 0000
      if (ISIGDEN.eq.1) then
        NTEST = max(NTESTL,IPRSIG)
      else if (ISIGDEN.eq.2) then
        NTEST = max(NTESTL,IPRDEN)
      else
        write(LUWRT,*) ' Illegal value for ISIGDEN: ',ISIGDEN
        call quit('*** ERROR in SIGDEN_REL ***')
      end if
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.50) THEN 
        write(LUWRT,*)
        write(LUWRT,*)
        write(LUWRT,*) '++++++++++++++++++++++++++++++++++++++'
        write(LUWRT,*) '++++++++++++++++++++++++++++++++++++++'
        WRITE(LUWRT,*) '       SIGDEN_REL in action           '
        write(LUWRT,*) '++++++++++++++++++++++++++++++++++++++'
        write(LUWRT,*) '++++++++++++++++++++++++++++++++++++++'
        WRITE(LUWRT,*) '          ICSM, ISSM = ', ICSM, ISSM
      END IF
#endif
*
*     file LUSC37 is not needed in a MPI run
*
#if !defined (VAR_MPI2)
      Open(Unit=LUSC37,File=NLUSC37_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
#endif
*
*=====================================================
* Implementing time-reversal symmetry 
* in the even-particle case, according to
*   Cm_ = Cm*
*   Cm  = Cm_*
* No signs here, phase factors always 1. (In odd-particle
* case -1 for Nbeta odd, Nalpha even!)
*
*  Timo Fleig, 31.07.2001
*
* define symmetry, factors, and loop restrictions
      call evenodd(IEONELEC,NACTEL)
      if (IEONELEC.eq.1) then
        IKRA_SYM = 0
      else if (IEONELEC.eq.2) then
        IKRA_SYM = 1
      end if
      IKRA_SYM = 0         ! FOR THE MOMENT !
* 
*=====================================================
*
      LBLK = -1
*
*     LUCIAREL is a relativistic CI program. 
*     However, for quaternion real groups, we use real
*     algebra.
      if (NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
*=====================================================
#if defined (VAR_MPI2)
*
*     make an update of list ISCLFAC_GROUP using LUCLIST...
      BLOCKTIME = 0.0D0
      CALL UPDATE_LUC_LIST(ISCLFAC_GROUP,LUCLIST,RCCTOS,C,
     &                     NPARBLOCK_C,IBLOCKL_C,IGROUPLIST,
     &                     IPROCLIST,IRILP,BLOCKTIME)
*
*     debug printing
      NPTESTVAR = 00
      IF( NPTESTVAR .ge. 20 ) THEN
      IOFFSET_SCRATCH = 0
      WRITE(LUWRT,*) 'left vector = C-vector'
      DO I_COMPLEX = 1, IRILP
        IF( I_COMPLEX .eq. 1 ) THEN
          WRITE(LUWRT,*) ' real part'
        ELSE
          WRITE(LUWRT,*) ' imaginary part'
        END IF
        DO IBLK = 1, NUM_BLOCKS
*
          LENGTH_BLK = IBLOCKL_C(IBLK)
          IF(ISCLFAC_GROUP(((I_COMPLEX-1)*NUM_BLOCKS)+IBLK).ne.0)THEN
            CALL MPI_FILE_READ_AT(LLUC,IOFFSET_SCRATCH,C,LENGTH_BLK,
     &                            MPI_REAL8,ISTAT,IERR)
            WRITE(LUWRT,'(A,1X,I12,1X,A,1X,I6,1X,I6)') 'Read-in at',
     &            IOFFSET_SCRATCH,'for block', 
     &            (I_COMPLEX - 1 ) * NUM_BLOCKS  + IBLK,
     &            LENGTH_BLK
            CALL WRTMATMN(C,1,LENGTH_BLK,1,LENGTH_BLK,LUWRT)
          END IF
          IOFFSET_SCRATCH = IOFFSET_SCRATCH + LENGTH_BLK
*
*
          END DO
        END DO
      END IF
      NPTESTVAR = 00
      IOFFSET_SCRATCH = 0
*
*=====================================================
#else
*=====================================================

#if defined LUCI_DEBUG
      if (NTEST.ge.1000) then
        write(LUWRT,*) '--------------------------------------'
        write(LUWRT,*) ' Printing vectors in SIGDEN_REL.      '
        write(LUWRT,*) '--------------------------------------'
        write(LUWRT,*) '  Right vector: '
        call wrtvcdc(C,LLUC,1,IRILP,-1)
        if (ISIGDEN.eq.2) then
          write(LUWRT,*) '  Left vector: '
          call wrtvcdc(HC,LLUHC,1,IRILP,-1)
        end if
      end if
#endif
*
      CALL REWINE(LLUC,-1)
      CALL REWINE(LLUHC,-1)
*
#endif /* if def VAR_MPI2 */
*=====================================================
*
*
*     loop over vectors of sigma and C : MK2 space
*
      DO IRIS = 1,IRILP,1
#if defined (VAR_MPI2)
*       initialize writing (ISIGDEN == 1) resp. reading (ISIGDEN == 2)
        IOFF_SCR1 = 0
        IOFF_SCR2 = 0
        IOFF_SCR1 = MY_LU2_OFF
     &            + MY_VEC2_IOFF * ( JVEC_SF )
     &            + MY_VEC1_IOFF * ( IRIS -1 )
        INT_IOFF1  = 0
        INT_IOFF2  = 0
*
        INT_IOFF1  = 1 + MY_ACT_BLK2 * ( JVEC_SF )
     &                 + MY_ACT_BLK1 * ( IRIS -1 )
*
C?      WRITE(LUWRT,*) 'INT_IOFF1 in SIGDEN',INT_IOFF1
*
*
        NUM_BLK_ACT_S   = 0
        IOFFSET_SCRATCH = 0
*
#endif /* if def VAR_MPI2 */
*
#if defined LUCI_DEBUG
        if (NTEST.ge.50) then
          write(LUWRT,*)
          write(LUWRT,*) '-----------------------'
          write(LUWRT,*) ' SIGMA LOOP , IRIS = ',IRIS
          write(LUWRT,*) '-----------------------'
        end if
#endif
*
*       count the 'active' sigma blocks
* 
        NUM_SBLK_ACT_BT = 0
*
        JBATABS = 0
        DO IMK2_S = 1, NMS2VAL
          MK2_S = MS2VAL(IMK2_S)
#if defined LUCI_DEBUG
          if (NTEST.ge.50) then
            write(LUWRT,*) '- - - - - - - - - - - -'
            write(LUWRT,*) ' MK2_S LOOP , MK2_S = ',MK2_S
            write(LUWRT,*) '- - - - - - - - - - - -'
          end if
#endif
*
          IASTP = IST_FOR_DT(1,IMK2_S)
          IBSTP = IST_FOR_DT(2,IMK2_S)
          NOCTPA_S = NOCTYP(IASTP)
          NOCTPB_S = NOCTYP(IBSTP)
*
*         arrays giving allowed type combinations (of alpha and beta)
*         (integer!)
          call memmar(KSIOIO,NOCTPA_S*NOCTPB_S,'ADDL  ',1,'SIOIO ')
C         call memmar(KSIOIO,NOCTPA_S*NOCTPB_S,'ADDL  ',2,'SIOIO ')
          CALL IAIBCM_REL(ISSPC,IASTP,IBSTP,WORK(KSIOIO))
*
*         arrays giving block type (integer!)
          KSVST = 1
          call memmar(KSBLTP,NSMST,'ADDL  ',1,'SBLTP ')
C         call memmar(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
          CALL ZBLTP_REL(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP))
*
* Loop over batches of sigma for this MK2_S projection value
          do JBATS = 1,NBAT_MS2(IMK2_S)
*
            JBATABS = JBATABS + 1
*
#if defined (VAR_MPI2)
*
*           start various timings...
            starttime  = mpi_wtime()
            sbatchtime = 0.0D0
            writetime  = 0.0D0
            KBATSEND = I1BATS(JBATABS) + LBATS(JBATABS)
*
*           set new offsets for reading from ILU2
            IF( ISIGDEN .eq. 2 )THEN
              IOFF_SCR1 = IOFF_SCR1 + IOFFSET_SCRATCH
              INT_IOFF1 = INT_IOFF1 + NUM_BLK_ACT_S 
            END IF
            ICOMPUTE = 0
            LS = 0
            NUM_BLK_ACT_S   = 0
            IOFFSET_SCRATCH = 0
*
*           check whether a cpu needs to call sigden_rel2
            DO ISBLK = I1BATS(JBATABS),KBATSEND-1,1
              IF( NPARBLOCK_S(ISBLK) .eq. MYPROC ) THEN 
                ICOMPUTE = 1
                NUM_BLK_ACT_S = NUM_BLK_ACT_S + 1
                IOFFSET_SCRATCH = IOFFSET_SCRATCH + IBLOCKL_S(ISBLK)
              END IF
            END DO
CSK         WRITE(LUWRT,*) 'ICOMPUTE for batch JBATABS is',
CSK  &                      ICOMPUTE,JBATABS
            IF( ICOMPUTE .eq. 0 ) THEN
              starttimer = MPI_WTIME()
              GOTO 60
            ENDIF
#endif /* if def VAR_MPI2 */
C
            LS = LEBATS(JBATABS)
            CALL DZERO(HC,LS)
            if (ISIGDEN.eq.1) then
#if defined LUCI_DEBUG
              if (NTEST.ge.500) then
                 WRITE(LUWRT,*)
                 WRITE(LUWRT,*)
     &           ' =============================================='
                 WRITE(LUWRT,*)
     &           ' Start of construction of Sigma batch ', JBATABS
                 write(LUWRT,'(2X,A,I3)') ' for MK2_S value ',MK2_S
                 WRITE(LUWRT,*) 
     &           ' =============================================='
                 WRITE(LUWRT,*)
              end if
#endif
            else 
C             ... sigma blocks read from disk in density run?
C             this is important since we might skip complete s blocks 
C             because of MK2_DIFF > 4.
              SBLK_READS2 = .FALSE.
            end if
C
#if !defined (VAR_MPI2)
            call REWINE(LLUC,-1)
#endif
            DO IRIC = 1,IRILP,1
#if defined LUCI_DEBUG
            if (NTEST.ge.50) then
              write(LUWRT,*)
              write(LUWRT,*) '-----------------------'
              write(LUWRT,*) ' C     LOOP , IRIC = ',IRIC
              write(LUWRT,*) '-----------------------'
            end if
#endif
#if defined (VAR_MPI2)
*             initialize reading
              IOFF_SCR1_C = 0
              IOFF_SCR2_C = L_COMBI
              IOFF_SCR1_C = IOFF_SCR2_C * ( IRIC -1 )
              IOFF_SCR2_C = 0
*
              INT_IOFF1_C = 0
              INT_IOFF2_C = 0
*
              INT_IOFF1_C = 1 + MY_ACT_BLK_ALL * ( IRIC -1 )
              NUM_BLK = 0
*
#else
              IF(IRIC.EQ.2) THEN
*               skip EOV mark between real and imaginary part
                CALL IFRMDS(IONEM,1,-1,LLUC)
              END IF
#endif /* if def VAR_MPI2 */
*
              DO IMK2_C = 1, NMS2VAL
                MK2_C = MS2VAL(IMK2_C)
#if defined LUCI_DEBUG
                if (NTEST.ge.50) then
                  write(LUWRT,*) '- - - - - - - - - - - -'
                  write(LUWRT,*) ' MK2_C LOOP , MK2_C = ',MK2_C
                  write(LUWRT,*) '- - - - - - - - - - - -'
                end if
#endif
                MK2_DIFF = ABS(MK2_C - MK2_S)
                max_MK2_DIFF = 4
#if defined LUCI_DEBUG
                IF( NTEST .ge. 50 )THEN
                  IF (MK2_DIFF .gt. max_MK2_DIFF) THEN
                    write(LUWRT,*) 'MK2 skip, MK2_DIFF =',MK2_DIFF
                  END IF
                END IF
#endif
*               ... for any two-electron operator.
*                   For one-electron operators max_MK2_DIFF = 1 or 0
*                   (or zero, if NZ .lt. 4 and operator is not j- or
*                   k-imaginary) /HJAAJ Aug 2008. TODO TODO to define
*                   max_MK2_DIFF
*
* Kramers symmetry restriction
                if (IKRA_SYM.eq.1.and.MK2_C.lt.0) then
                  if (IRIC.eq.1) IKSFAC = 1
                  if (IRIC.eq.2) IKSFAC = -1
* copy relevant block of sigma vector to position and 
* multiply with IKSFAC.
* FIXME !!
                  stop 'KS not implemented yet.'
                else
                  IACTP = IST_FOR_DT(1,IMK2_C)
                  IBCTP = IST_FOR_DT(2,IMK2_C)
                  NBL_C = NBLK_MS2(IMK2_C)
C                 C and sigma may have different symmetry in general
                  IF( DOSIGPROP ) NBL_C = NBLK_MS2_C(IMK2_C)
#if defined (VAR_MPI2)
*                 set new offset
                  IOFF_SCR1_C = IOFF_SCR1_C + IOFF_SCR2_C
                  INT_IOFF1_C = INT_IOFF1_C + INT_IOFF2_C
CSK               WRITE(LUWRT,*) ' NEW OFFSET FOR C-VECTOR FILE',
CSK     &                          IOFF_SCR1_C, INT_IOFF1_C
#else
                  IF (MK2_DIFF .LE. max_MK2_DIFF) THEN
C                   copy blocks with given MK2_C to LUSC37
                    CALL REWINE(LUSC37,-1)
                    CALL COPNBLKD(LLUC,LUSC37,C,NBL_C,0,LBLK)
C                   end of vector mark
                    CALL ITODS(-1,1,-1,LUSC37)
                  ELSE
C                   skip blocks with given MK2_C
                    CALL COPNBLKD(LLUC,-1,C,NBL_C,0,LBLK)
                  END IF
#endif /* if def VAR_MPI2 */
                  IF (MK2_DIFF .LE. max_MK2_DIFF)
     &            call sigden_rel2(C,HC,JBATABS,LBATS,LEBATS,I1BATS,
     &                             IBATS,LUSC37,LLUHC,T,T_BUFF,
     &                             ISIGDEN,IASTP,IBSTP,IRIS,IRIC,
     &                             IMK2_C,IPRDEN,SBLK_READS2
#if defined (VAR_MPI2)
     &                            ,LUCLIST,IOFF_SCR1_C,INT_IOFF1_C,
     &                             ISCLFAC_GROUP,NPARBLOCK_C,IBLOCKL_C,
     &                             LUHCLIST,IOFF_SCR1,INT_IOFF1,LLUC,
     &                             IT_TTPL,IT_TTOL
     &                            ,IBLOCKL_S,NPARBLOCK_S
#endif
     &                             )
                end if
*               ^ End if Kramers restriction
#if defined (VAR_MPI2)
*               keep track of correct offset
                IOFF_SCR2_C = 0
                DO ICBLK = 1, NBL_C
                  IOFF_SCR2_C = IOFF_SCR2_C
     &                        + IBLOCKL_C( NUM_BLK + ICBLK )
                END DO
                NUM_BLK = NUM_BLK + NBL_C
                INT_IOFF2_C = NBL_C
#endif /* if def VAR_MPI2 */
              END DO
*             ^ End of loop over IMK2_C
            END DO
*           ^ End of loop over IRIC
*
 60       CONTINUE
            if(ISIGDEN.eq.1) then
#if defined LUCI_DEBUG
              if(NTEST.ge.100) then
                WRITE(LUWRT,'(/A)') ' ******************'
                WRITE(LUWRT,'( A)') ' *  final sbatch  *'
                WRITE(LUWRT,'(A/)') ' ******************'
                CALL WRTMATMN(HC,1,LEBATS(JBATABS),1,
     &                         LEBATS(JBATABS),LUWRT)
              end if
#endif
              KBATSEND = I1BATS(JBATABS) + LBATS(JBATABS)
#if defined (VAR_MPI2)
              starttimer = MPI_WTIME()
              IF( ICOMPUTE .eq. 0 ) GOTO 70
#endif /* if def VAR_MPI2 */
              do ISBLK = I1BATS(JBATABS),KBATSEND-1,1
#if defined (VAR_MPI2)
*               check if block is needed
                IF( NPARBLOCK_S(ISBLK) .ne. MYPROC ) GOTO 65
*
*               new offset
*
                IOFF_SCR1 = IOFF_SCR1 + IOFF_SCR2
                INT_IOFF1 = INT_IOFF1 + INT_IOFF2
*
                IOFF = IBATS(6,ISBLK)
                LEN_T = IBATS(8,ISBLK)
*               check norm of the block
                XXX = 0.0D0
                XXX = DDOT(LEN_T,HC(IOFF),1,HC(IOFF),1)
CSK             IF( XXX .eq. 0.0D0 ) THEN
CSK               WRITE(LUWRT,*) ' Norm of block', ISBLK,' is 0.0D0'
CSK             END IF
                IF( XXX .ne. 0.0D0 ) THEN
C                 mark block as nonzero
                  LUHCLIST(INT_IOFF1) = IONE
C                 transfer block to disc
!                 WRITE(LUWRT,*) 'I transfer to file block ISBLK',
!    &                            ISBLK
!                 CALL WRTMATMN(HC(IOFF),1,LEN_T,1,LEN_T,LUWRT)
!                 WRITE(LUWRT,*) 'offset,length, IRIS',
!    &                            IOFF_SCR1,LEN_T,IRIS
                  CALL MPI_FILE_WRITE_AT(LLUHC,IOFF_SCR1,HC(IOFF),LEN_T,
     &                                   MPI_REAL8,ISTAT,IERR)
                  ISI_CALC_BL = ISI_CALC_BL + 1
                END IF
*               keep track of correct offset
                IOFF_SCR2 = LEN_T
                INT_IOFF2 = 1
 65             CONTINUE
#else
                IOFF = IBATS(6,ISBLK)
                LEN_T = IBATS(8,ISBLK)
                call itods(LEN_T,1,-1,LLUHC)
                call todsc(HC(IOFF),LEN_T,-1,LLUHC)
                ISI_CALC_BL = ISI_CALC_BL + 1
#endif /* if def VAR_MPI2 */
              END DO
 70           CONTINUE
#if defined (VAR_MPI2)
              IF( ICOMPUTE .ne. 0 ) THEN
                writetime = writetime + MPI_WTIME() - starttimer
                sbatchtime = sbatchtime + MPI_WTIME() - starttime
CSK too much information?                WALLTID  = SECTID(writetime)
                WALLTID2 = SECTID(sbatchtime)
CSK too much information?                 WRITE(LUWRT,50000) WALLTID
                IF( .NOT. DOSIGPROP ) WRITE(LUWRT,60001) WALLTID2
CSK                WRITE(LUWRT,60000) WALLTID2, JBATABS
              END IF
#endif /* if def VAR_MPI2 */
            END IF
C           ^ ISIGDEN == 1
C
#if !defined (VAR_MPI2)
C           skip sigma blocks on disk
            IF( ISIGDEN .eq. 2 .and. .NOT. SBLK_READS2 )THEN
               ISKIP_S2_BLOCKS = 0
               do ISBLK = I1BATS(JBATABS),
     &                    I1BATS(JBATABS) + LBATS(JBATABS)-1
                  ISKIP_S2_BLOCKS = ISKIP_S2_BLOCKS + 1
               end do 
               CALL COPNBLKD(LLUHC,-1,C,ISKIP_S2_BLOCKS,0,LBLK)
            END IF
#endif /* if !def VAR_MPI2 */
          END DO
*         ^ End loop over S batches of given MK2_S value
        END DO
*       ^ End of loop over IMK2_S
*
#if !defined (VAR_MPI2)
        if (ISIGDEN.eq.1) then
* Write end of vector mark between real and imaginary part.
          call itods(-1,1,-1,LLUHC)
        else if (ISIGDEN.eq.2) then
* Skip end of vector mark between real and imaginary part.
          call ifrmds(IONEM,1,-1,LLUHC)
        end if
#endif /* if !def VAR_MPI2 */
*
      END DO
*     ^ End of loop over IRIS
      call memmar(KDUM ,IDUM,'FLUSM ',2,'SIGDEN')
*
      CALL QEXIT('SIGDE')
#if !defined (VAR_MPI2)
*     close scratch unit for Right vector
      close(unit=LUSC37,status='DELETE')
*
#else
*
      WALLTID3 = SECTID(BLOCKTIME)
      IF( .NOT. DOSIGPROP ) WRITE(LUWRT,80000) WALLTID3
      IF( SHARED_M .and. ISIGDEN .eq. 1 )THEN
        WALLTID3 = SECTID(mem_cp_time)
        WRITE(LUWRT,90000) WALLTID3
      END IF
C
      IF( TIMING .and. (.NOT. DOSIGPROP)) THEN
C
C       print statistics
C
        WRITE(LUWRT,'(/A)')
     &  '               H x b_i contraction run statistics  '
        WRITE(LUWRT,'(A/)')
     &  '              ____________________________________ '
        WRITE(LUWRT,'(2X,A,1X,I9)')
     &  ' number of s_i blocks calculated          : ',ISI_CALC_BL
        WRITE(LUWRT,'(2X,A,1X,I9)')
     &  ' total number of b_i read from disk       : ',IBI_MULT_BL
        WALLTID = SECTID(xreadtimebi)
        WRITE(LUWRT,'(2X,A,1X,A)')
     &  ' read time for b_i blocks from disk       : ', WALLTID
        WALLTID = SECTID(xcomputesi)
        WRITE(LUWRT,'(2X,A,1X,A/)')
     &  ' matrix x vector multiplication time      : ', WALLTID
C
      END IF
C
*
50000 FORMAT(' >>>  WALL TIME FOR WRITING SB TO DISK           :
     & ',A)
60000 FORMAT(' >>>  WALL TIME IN TOTAL FOR CURRENT S-BATCH     :
     & ',A,I4)
60001 FORMAT(' >>>  WALL TIME IN TOTAL FOR CURRENT S-BATCH     :
     & ',A)
80000 FORMAT(' >>>  WALL TIME FOR C-COEFFICIENT EXCHANGE       :
     & ',A)
90000 FORMAT(' >>>  WALL TIME FOR OBTAINING T-ELEMENTS         :
     & ',A)
*
#endif /* if !def VAR_MPI2 */
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SIGDEN_REL2(CB,SB,JBATABS,LBATS,LEBATS,I1BATS,IBATS,
     &                       LUC,LUHC,T,T_BUFF,ISIGDEN,IASTP,IBSTP,
     &                       IRIS,IRIC,IMK2_C,IPRNT,SBLK_READS2
#if defined (VAR_MPI2)
     &                      ,LUCLIST,IOFF_SCR1_C,INT_IOFF1_C,
     &                       ISCLFAC_GROUP,NPARBLOCK_C,IBLOCKL_C,
     &                       LUHCLIST,IOFF_SCR1,INT_IOFF1,LLUC,
     &                       IT_TTPL,IT_TTOL
     &                      ,IBLOCKL_S,NPARBLOCK_S
#endif
     &                       )
*
* First inner routine for general CC sigma/densi
*
* Complete S and C vectors
*
* Jeppe Olsen   Relativistic version, Winter of 2001
*
* =====
* Input
* =====
*
#include "implicit.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE), NPARBLOCK_C(*)
      DIMENSION LUCLIST(*), ISCLFAC_GROUP(*), IBLOCKL_C(*)
      DIMENSION LUHCLIST(*), IBLOCKL_S(*), NPARBLOCK_S(*) 
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR1, IOFF_SCR2 
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR_EXPL
      INTEGER INT_IOFF1, INT_IOFF2, INT_IOFF3, INT_IOFF_EXPL
#endif
#include "parluci.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "gasstr.inc"
#include "cstate.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "csm.inc"
*. Batches of sigma
      INTEGER LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
*.Scratch
      DIMENSION SB(*),CB(*)
      dimension T(*)
      LOGICAL SBLK_READS2
*
      NTESTL = 000
      NTEST = MAX(NTESTL,IPRNT)
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*) ' ================='
        WRITE(LUWRT,*) ' SIGDEN2 speaking :'
        WRITE(LUWRT,*) ' ================='
      END IF
#endif
#if defined (VAR_MPI2)
*     initialize scratch variables      
      IOFF_SCR2     = 0
      INT_IOFF3     = 0
      INT_IOFF_EXPL = 0
      IOFF_SCR_EXPL = 0
      ICOUNT_ACT    = 0
#endif
*
*      ISIGDEN = 2 : Read in batch of blocks of lhs vector
*      unique blocks are read in
*
#if defined (VAR_MPI2)
      INT_IOFF2 = 0
      INT_IOFF3 = I1BATS(JBATABS)
CSK      WRITE(LUWRT,*) 'NEW INT_IOFF3 in SIGDEN2 FOR JBATABS',INT_IOFF3,
CSK     &                JBATABS
      INT_IOFF_EXPL = INT_IOFF1
      IOFF_SCR_EXPL = IOFF_SCR1
#endif
*
      IF (ISIGDEN .eq. 2 )THEN
csk     if (IRIC .eq. 1 .and. IMK2_C .eq. 1 )then
        IF ( .NOT. SBLK_READS2 )THEN
C
          DO 100 ISBLK = I1BATS(JBATABS),
     &                   I1BATS(JBATABS)+LBATS(JBATABS)-1
#if defined (VAR_MPI2)
            IF( NPARBLOCK_S( ISBLK ) .ne. MYPROC ) GOTO 100
C           ... active CPU-block found
            ICOUNT_ACT = ICOUNT_ACT + 1
#endif
            ISATP = IBATS(1,ISBLK)
            ISBTP = IBATS(2,ISBLK)
            ISASM = IBATS(3,ISBLK)
            ISBSM = IBATS(4,ISBLK)
            ISOFF = IBATS(5,ISBLK)
            LEN2  = IBATS(7,ISBLK)
#if defined (VAR_MPI2)
C           CPU-block - test if non-zero
            IF( LUHCLIST( INT_IOFF_EXPL + ICOUNT_ACT - 1 ).eq.0) GOTO 80
            CALL GSTTBLD_PAR2(SB(ISOFF),LUHC,IOFF_SCR_EXPL,LEN2)
#else
            CALL GSTTBLD2(SB(ISOFF),LUHC)
#endif
*
C           NTEST = 100
C           IF ( LEN2 .gt. 0 .and. NTEST .ge. 100) THEN
C             WRITE(LUWRT,'(/A,2I6)') 
C    &        '  lhs vector block:',ISBLK, ISOFF
C             CALL WRTMATMN(SB(ISOFF),1,LEN2,1,LEN2,LUWRT)
C           END IF
C           NTEST = 000
  80        CONTINUE
#if defined (VAR_MPI2)
            IOFF_SCR2 = LEN2
            IOFF_SCR_EXPL = IOFF_SCR_EXPL + IOFF_SCR2
#endif
 100      CONTINUE
*         ^ end of loop over S-blocks in batch
          SBLK_READS2 = .TRUE.
        END IF
*       ^ only perform read if new S batch requested
      END IF
*     ^ only for density evaluation
*
*     obtain sigma/density for batch of blocks
      CALL SIGDEN_REL3(LBATS(JBATABS),IBATS(1,I1BATS(JBATABS)),1,
     &                 CB,SB,LUC,T,T_BUFF,ISIGDEN,IRIS,IRIC
#if defined (VAR_MPI2)
     &                ,LUCLIST,IOFF_SCR1_C,INT_IOFF1_C,
     &                 ISCLFAC_GROUP,NPARBLOCK_C,IBLOCKL_C,
     &                 LUHCLIST,INT_IOFF3,LLUC,
     &                 IT_TTPL,IT_TTOL
     &                ,IBLOCKL_S,NPARBLOCK_S
#endif
     &                 )
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SIGDEN_REL3(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,T,T_BUFF,
     &                       ISIGDEN,IRIS,IRIC
#if defined (VAR_MPI2)
     &                      ,LUCLIST,IOFF_SCR1_C,INT_IOFF1_C,
     &                       ISCLFAC_GROUP,NPARBLOCK_C,IBLOCKL_C,
     &                       LUHCLIST,INT_IOFF3,LLUC,
     &                       IT_TTPL,IT_TTOL
     &                      ,IBLOCKL_S,NPARBLOCK_S
#endif
     &                       )
*
* Generate a set of sigma blocks, 
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB
*
*============================================
* Proper allocation of resolution quantities
*    Timo Fleig, March 25, 2003
* Revised, more efficient.
*    Timo Fleig, March 31, 2003
*============================================
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
#include "cands.inc"
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
#include "ctcc.inc"
      COMMON/CMXCJ/MXCJ
*
      COMMON/HIDSCRR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
C
C     integer array may become very long
C
      INTEGER*8 INTSCR
C
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "cintfo.inc"
#include "parluci.h"
#if defined (VAR_MPI2)
#include "krmc_shmem.h"
#include "mxdim_mpi2.h"
      POINTER (MY_T_SCR_PTR, T_SCR_CC(IDUMMY_BUFF_SZ))
#endif
      DIMENSION CB(*), HCB(*)
C
C     set memory mark for local scratch
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'SIGD3X')
C
C     offsets for alpha and beta supergroups
      NOCTPA_C = NOCTYP(IACTP)
      NOCTPB_C = NOCTYP(IBCTP)
C
      NAEL_C = NELEC(IACTP)
      NBEL_C = NELEC(IBCTP)
C
C     arrays for storing NEL consecutive annihilations/creations
C     (resolution matrices)
C     number of alpha/beta excitation types
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &                      1,NAOBEX_TP,NBOBEX_TP,IDUMMY,IDUMMY)
C     and the alpha-and beta-excitations (integer arrays...)
      LENA = 2*NGAS*NAOBEX_TP
      LENB = 2*NGAS*NBOBEX_TP
      call memmar(KLAOBEX,LENA,'ADDL  ',1,'IAOBEX')
      call memmar(KLBOBEX,LENB,'ADDL  ',1,'IAOBEX')
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &                      0,NAOBEX_TP,NBOBEX_TP,
     &                      WORK(KLAOBEX),WORK(KLBOBEX))
*
*     for alpha excitations, C and S vector (might be flipped around!)
      IOCTPA_C = IBSPGPFTP(IACTP)
      CALL LEN_GENOP_STR_MAP(NAOBEX_TP,WORK(KLAOBEX),NOCTPA_C,
     &                       NELFSPGP(1,IOCTPA_C),NOBPT,NGAS,
     &                       MAXLENAC)
      IOCTPA_S = IBSPGPFTP(IASTP)
      NOCTPA_S = NOCTYP(IASTP)
      CALL LEN_GENOP_STR_MAP(NAOBEX_TP,WORK(KLAOBEX),NOCTPA_S,
     &                       NELFSPGP(1,IOCTPA_S),NOBPT,NGAS,
     &                       MAXLENAS)
*
*     for beta excitations
      IOCTPB_C = IBSPGPFTP(IBCTP)
      CALL LEN_GENOP_STR_MAP(NBOBEX_TP,WORK(KLBOBEX),NOCTPB_C,
     &                       NELFSPGP(1,IOCTPB_C),NOBPT,NGAS,
     &                       MAXLENBC)
      IOCTPB_S = IBSPGPFTP(IBSTP)
      NOCTPB_S = NOCTYP(IBSTP)
      CALL LEN_GENOP_STR_MAP(NBOBEX_TP,WORK(KLBOBEX),NOCTPB_S,
     &                       NELFSPGP(1,IOCTPB_S),NOBPT,NGAS,
     &                       MAXLENBS)
C
C     release pure local scratch
      IDUM = 0
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'SIGD3X')
C
      MAXLEN_I1 = MAX(MAXLENAC,MAXLENBC,MAXLENAS,MAXLENBS)
      LSCR3 = MAXLEN_I1
C
C     non-local scratch - used in called subroutines
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'SIGD3Y')
C
      call memmar(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      call memmar(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
C
      call memmar(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      call memmar(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
C
      call memmar(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      call memmar(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
C
      call memmar(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      call memmar(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
C
C     T-elements
      INTSCR = MXINKA**4
      INTSCR = MIN(INTSCR,MX_TBLK_MX)
      CALL MEMMAR_I8(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
C
C     arrays for storing occupations of T-operators
      call memmar(KTOCC1,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'TOCC1 ')
      call memmar(KTOCC2,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'TOCC2 ')
      call memmar(KTOCC3,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'TOCC3 ')
      call memmar(KTOCC4,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'TOCC4 ')
C
C     arrays for storing NEL consecutive annihilations/creations
      LSCR4 = MXINKA*MXINKA*MXINKA
      call memmar(KI1G  ,LSCR4,'ADDL  ',1,'I1G   ')
      call memmar(KXI1G,LSCR4,'ADDL  ',2,'XIG   ')
C
      call memmar(KI2G  ,LSCR4,'ADDL  ',1,'I2G   ')
      call memmar(KXI2G,LSCR4,'ADDL  ',2,'XI2G  ')
C
      call memmar(KI3G  ,LSCR4,'ADDL  ',1,'I3G   ')
      call memmar(KXI3G,LSCR4,'ADDL  ',2,'XI3G  ')
C
      call memmar(KI4G  ,LSCR4,'ADDL  ',1,'I4G   ')
      call memmar(KXI4G,LSCR4,'ADDL  ',2,'XI4G  ')
C
      call memmar(KI1GE  ,LSCR4,'ADDL  ',1,'I1GE  ')
      call memmar(KXI1GE,LSCR4,'ADDL  ',2,'XIGE  ')
C
      call memmar(KI2GE  ,LSCR4,'ADDL  ',1,'I2GE  ')
      call memmar(KXI2GE,LSCR4,'ADDL  ',2,'XI2GE ')
C
      call memmar(KI3GE  ,LSCR4,'ADDL  ',1,'I3GE  ')
      call memmar(KXI3GE,LSCR4,'ADDL  ',2,'XI3GE ')
C
      call memmar(KI4GE  ,LSCR4,'ADDL  ',1,'I4GE  ')
      call memmar(KXI4GE,LSCR4,'ADDL  ',2,'XI4GE ')
C     TTS arrays for partitioning of vector 
      NTTS = MXNTTS
      call memmar(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      call memmar(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      call memmar(KLI1BT,2*NTTS,'ADDL  ',1,'I1BTC ')
      call memmar(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
C
C     arrays giving allowed type combinations
      call memmar(KCIOIO,NOCTPA_C*NOCTPB_C,'ADDL  ',1,'CIOIO ')
      CALL IAIBCM_REL(ICSPC,IACTP,IBCTP,WORK(KCIOIO))
C
C     arrays giving block type
      call memmar(KCBLTP,NSMST,'ADDL  ',1,'CBLTP ')
      CALL ZBLTP_REL(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP))
C     use all symmetryblocks of given TT
      ITTSS_ORD = 2
C
C     local scratch arrays for blocks of C and sigma
      L0BLOCK = MXSOOB
      LSCR1 = MAX(L0BLOCK,LBLOCK)
#if defined (VAR_MPI2)
C
C     set proper C-block offset...
      ISTART_BLOCK = INT_IOFF1_C - 1
      IDEBUGPRNT = 00
      CALL PART_CIV_PAR3(IDC,WORK(KCBLTP),WORK(KNSTSO(IACTP)),
     &     WORK(KNSTSO2(IBCTP)),NOCTPA_C,NOCTPB_C,NSMST,LSCR1,
     &     WORK(KCIOIO),ISMOST(1,ICSM),NCBATCH,WORK(KLLBT),
     &     WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),0,ITTSS_ORD,
     &     ISCLFAC_GROUP,ISTART_BLOCK,0,IDEBUGPRNT)
C           use this line as line 2 for testing more batches ...
C     &     WORK(KNSTSO2(IBCTP)),NOCTPA_C,NOCTPB_C,NSMST,L0BLOCK,
#else
      CALL PART_CIV2_SPC(IDC,WORK(KCBLTP),WORK(KNSTSO(IACTP)),
     &     WORK(KNSTSO2(IBCTP)),NOCTPA_C,NOCTPB_C,NSMST,LSCR1,
     &     WORK(KCIOIO),ISMOST(1,ICSM),NCBATCH,WORK(KLLBT),
     &     WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),0,ITTSS_ORD)
#endif
C
ctf   CALL MEMMAR(KLISMDST,MAX_NSYMDIST,  'ADDL  ',1,'ISMDST')
      CALL MEMMAR(KLISMDST,MAX_NSYMDIST+1,  'ADDL  ',1,'ISMDST')
C
      LSIRES = MXINKA**4
      LCJRES = LCJRES_SAVE
C
C     allocation of resolution matrices: sigma + C
      call memmar(KSIRES,LSIRES,'ADDL  ',2,'KSIRES')
      call memmar(KCJRES,LCJRES,'ADDL  ',2,'KCJRES')
C
C     set offset for T_BUFF array
      LEN_T_BUFF_NZ = LEN_T_VEC
#if defined (VAR_MPI2)
C
C     dynamic allocation for scratch T block - MPI info if non-zero
      IF( IT_SHL .ge. 0 ) CALL GET_MEM_T_SCR_CC(0)
#endif
C
      CALL SIGDEN_REL4(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,
     &                 WORK(KNSTSO(IASTP)),WORK(KNSTSO2(IBSTP)),
     &                 WORK(KNSTSO(IACTP)),WORK(KNSTSO2(IBCTP)),
     &                 IASTP,IBSTP,IACTP,IBCTP,MXINKA,WORK(KINSCR),
     &                 WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
     &                 WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),
     &                 IPRSIG,IPRDEN,LUC,WORK(KCJRES),WORK(KSIRES),
     &                 WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),
     &                 WORK(KLIBT),ISSM,ICSM,ISIGDEN,NCBATCH,T,T_BUFF,
     &                 IDIAG,IRIS,IRIC,WORK(KLSOBEX),WORK(KLSPOBEX_AC),
     &                 WORK(KTOCC1),WORK(KTOCC2),WORK(KTOCC3),
     &                 WORK(KTOCC4),WORK(KI1G),WORK(KXI1G),WORK(KI2G),
     &                 WORK(KXI2G),WORK(KI3G),WORK(KXI3G),WORK(KI4G),
     &                 WORK(KXI4G),WORK(KI1GE),WORK(KXI1GE),WORK(KI2GE),
     &                 WORK(KXI2GE),WORK(KI3GE),WORK(KXI3GE),
     &                 WORK(KI4GE),WORK(KXI4GE),WORK(KSIGNNHX),
     &                 WORK(KLISMDST)
#if defined (VAR_MPI2)
     &                ,LUCLIST,IOFF_SCR1_C,INT_IOFF1_C,
     &                 ISCLFAC_GROUP,NPARBLOCK_C,IBLOCKL_C,
     &                 LUHCLIST,INT_IOFF3,ISTART_BLK,LLUC,T_SCR_CC,
     &                 IT_TTPL,IT_TTOL
     &                ,IBLOCKL_S,NPARBLOCK_S
#endif
     &     )
C
C     eliminate local memory
      IDUM = 0
      call memmar(KDUM ,IDUM,'FLUSM ',2,'SIGD3Y')
C
#if defined (VAR_MPI2)
      IF( IT_SHL .ge. 0 ) CALL REL_MEM_T_SCR(T_SCR_CC)
#endif
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SIGDEN_REL4(NSBLOCK,ISBLOCK,CB,SB,
     &           NSSOA_S,NSSOB_S,NSSOA_C,NSSOB_C,
     &           IASTP,IBSTP,IACTP,IBCTP,MAXK,
     &           XINT,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           IPRS,IPRD,LUC,CJRES,SIRES,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,
     &           ISSM,ICSM,ISIGDEN,NCBATCH,
     &           T,T_BUFF,
     &           IDIAG,IRIS,IRIC,ITSPOBEX_TP,ISPOBEX_AC,
     &           ITOCC1,ITOCC2,ITOCC3,ITOCC4,
     &           I1G,XI1G,I2G,XI2G,I3G,XI3G,I4G,XI4G,
     &           I1GE,XI1GE,I2GE,XI2GE,I3GE,XI3GE,I4GE,XI4GE,
     &           SIGN_NHX,ISMDST
#if defined (VAR_MPI2)
     &          ,LUCLIST,IOFF_SCR1_C,INT_IOFF1_C,
     &           ISCLFAC_GROUP,NPARBLOCK_C,IBLOCKL_C,
     &           LUHCLIST,INT_IOFF3,ISTART_BLK,LLUC,T_SCR_CC,
     &           IT_TTPL,IT_TTOL
     &          ,IBLOCKL_S,NPARBLOCK_S
#endif
     &           )
*
*
* Jeppe Olsen, Summer of 1999
*              July 2000, Combinations of T-amplitudes and 
*                         determinants introduced
*              September 2000, ISPOBEX_FRZ added
*
* Timo Fleig, 2001/2002
*                  Modifications for relativistic CI
*             Aug. 2002
*                  Real/Imag. factor fixed for S_r (Jeppe + Timo)
*             Aug. 2003
*                  ISMDST added 
*                                 
* Stefan Knecht, May - June 2007,
*                   parallel adaption
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included 
*
* NSSOA : Number of strings per type and symmetry for alpha strings
* NSSOB : Number of strings per type and symmetry for beta strings
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
*
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE), NPARBLOCK_S(*), IBLOCKL_S(*)
      DIMENSION LUCLIST(*), ISCLFAC_GROUP(*), LUHCLIST(*)
      DIMENSION NPARBLOCK_C(*), IBLOCKL_C(*), IT_TTPL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR_C, IOFF_SCR_C_EXPL
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR1_C, IOFF_SAVE
      INTEGER(KIND=MPI_OFFSET_KIND) IT_TTOL(*)
      INTEGER(KIND=MPI_ADDRESS_KIND) ICPL_I8, ITB_I8
#endif
#include "parluci.h"
C     general input
#include "mxpdim.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "gasstr.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "cstate.inc"
#include "strbas.inc"
#include "integrals_off.inc"
#include "ctcc.inc"
C     memory window
#include "krmc_shmem.h"
C     specific input
      INTEGER ISBLOCK(8,*)
      INTEGER LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
C     general input
      INTEGER NSSOA_S(NSMST,*),NSSOB_S(NSMST,*)
      INTEGER NSSOA_C(NSMST,*),NSSOB_C(NSMST,*)
      DIMENSION SIGN_NHX(NSPOBEX_TP)
C     scratch
      DIMENSION SB(*),CB(*)
      DIMENSION XINT(*)
      DIMENSION I1(*),I2(*),I3(*),I4(*)
      DIMENSION XI1S(*),XI2S(*),XI3S(*)
      DIMENSION ISMDST(*)
C
      DIMENSION CJRES(*),SIRES(*)
C
C     T-coefficients
C
      DIMENSION T(*), T_BUFF(*)
#if defined (VAR_MPI2)
      DIMENSION T_SCR_CC(*)
#endif

      dimension ITSPOBEX_TP(4*NGAS,*)
      INTEGER ISPOBEX_AC(*), ICPL
      INTEGER*8 ITB, ITB_BUFF, ITB_BUFF_ADD
C
      ONEM = -1.0D0
      IZERO = 0
      IST_C_BLK_CT = 0
      ICPL         = 0
      ICPL_I8      = 0
      ITB_I8       = 0
C     process id
      ID_ACCUM     = -1 
      ID_GET       = -1 
C
      NTESTL = 000
      NTEST = MAX(NTESTL,max(IPRS,IPRD))
C
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.800) THEN
        WRITE(LUWRT,*) ' ================'
        WRITE(LUWRT,*) ' SIGDEN4 speaking :'
        WRITE(LUWRT,*) ' ================'
        WRITE(LUWRT,*)  
        WRITE(LUWRT,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(LUWRT,*) ' TTSS and IOFF for each ACTIVE sigma block'
          DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0) 
     &      WRITE(LUWRT,'(10X,I6,4I3,I8)') 
     &      IBLOCK,(ISBLOCK(II,IBLOCK),II=1,5)
          END DO
      END IF
#endif
      NTEST = 000
*
#if !defined (VAR_MPI2)
#if defined LUCI_DEBUG
      IF(NTEST.GE.500) THEN
        WRITE(LUWRT,*) ' Initial C vector '
        CALL WRTVCD(CB,LUC,1,-1)
      END IF
#endif
*     rewind C file
      REWIND LUC
#else
      IOFF_SCR_C = 0
      IOFF_SCR_C_EXPL = 0
      IOFF_SCR_C_EXPL = IOFF_SCR1_C
      IOFF_SAVE = IOFF_SCR1_C
      start_timer = 0.0D0
CSK      WRITE(LUWRT,*) ' GOT NEW OFFSET FOR C-VECTOR FILE',
CSK     &                 IOFF_SAVE
CSK      WRITE(LUWRT,*) 'INT_IOFF3 in SIGDEN4',INT_IOFF3
CSK      WRITE(LUWRT,*) ' ISCLFAC_GROUP in SIGDEN4'
CSK      CALL IWRTMAMN(ISCLFAC_GROUP,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUWRT)
#endif
C
      IOCTPA = IBSPGPFTP(IASTP)
      IOCTPB = IBSPGPFTP(IBSTP)
      JOCTPA = IBSPGPFTP(IACTP)
      JOCTPB = IBSPGPFTP(IBCTP)
C
      NSTT_BLK = NSBLOCK/NSMST
C
C     loop over batches over C blocks      
C
      DO 20000 JCBATCH = 1, NCBATCH             
C
C      C blocks offset in current batch
C
        NJBLOCK = LCBLOCK(JCBATCH)
        NCTT_BLK = NJBLOCK/NSMST
        I1C = I1CBLOCK(JCBATCH)
CSK     WRITE(LUWRT,*) ' JCBATCH, NJBLOCK, NCTT_BLK, I1C',
CSK  &                   JCBATCH, NJBLOCK, NCTT_BLK, I1C
C
C       loop over TT blocks of sigma and C in batches and 
C       obtain contributions 
C
        DO 9000 ICTT_BLK = 1, NCTT_BLK
          IREAD = 0
C
C         first block of next TT block of C    
C
          ICBLK = I1C + (ICTT_BLK-1)*NSMST
          JATP  = ICBLOCK(1,ICBLK)
          JBTP  = ICBLOCK(2,ICBLK)
          ICOFF = ICBLOCK(5,ICBLK)
C
#if defined LUCI_DEBUG
          IF(NTEST.GE.500) THEN
            WRITE(LUWRT,*) ' Next block of C, ICBLK,JATP,JBTP,ICOFF ',
     &      ICBLK,JATP,JBTP
          END IF
#endif
C
#if defined (VAR_MPI2)
C
C         check if TT-blocks are 'active'! if not, 
C         skip all loops that will follow.
C
          ICOMPUTE_C = 0
          NCBLK = NSMST
          IC_OFF_BLK = INT_IOFF1_C
          IF( IRIC .eq. 1 ) THEN
            IST_C_BLK_CT = IC_OFF_BLK
          ELSE
            IST_C_BLK_CT = IC_OFF_BLK - NUM_BLOCKS
          END IF
          DO JBLK = ICBLK, ICBLK + NCBLK - 1
CSK         WRITE(LUWRT,*) ' Checking for c block',IC_OFF_BLK+JBLK-1
            IF( ISCLFAC_GROUP(IC_OFF_BLK+JBLK-1) .gt. 0 ) THEN
              ICOMPUTE_C = 1
CSK           WRITE(LUWRT,*) ' active c block found',IC_OFF_BLK+JBLK-1
            END IF
          END DO
          IF( ICOMPUTE_C .eq. 0 ) GOTO 8999
*
#endif
          IRATP = JATP + JOCTPA - 1
          IRBTP = JBTP + JOCTPB - 1
C
          DO 10000 ISTT_BLK = 1, NSTT_BLK
C
C           first block of next TT block of sigma
            ISBLK = (ISTT_BLK-1)*NSMST + 1
#if defined (VAR_MPI2)
C
C           check if TT-blocks need to be calculated! if not,
C           skip all loops that follow.
            ICOMPUTE_S = 0
            NSBLK = NSMST
            IS_OFF_BLK = 0
            IS_OFF_BLK = INT_IOFF3
            DO JJBLK = ISBLK, ISBLK + NSBLK - 1
             IF(NPARBLOCK_S(JJBLK+IS_OFF_BLK-1).eq.MYPROC) ICOMPUTE_S=1
            END DO
            IF( ICOMPUTE_S .eq. 0 ) GOTO 9999
#endif
            IATP = ISBLOCK(1,ISBLK)
            IBTP = ISBLOCK(2,ISBLK)
            ISOFF = ISBLOCK(5,ISBLK)
CSK         WRITE(LUWRT,*) ' Next block of S, ISBLK,IATP,IBTP,ISOFF ',
CSK  &                                        ISBLK,IATP,IBTP,ISOFF
            ILATP = IATP + IOCTPA - 1
            ILBTP = IBTP + IOCTPB - 1
C
C           loop over spinorbital-excitation types 
            ITB_BUFF = 1
C
            DO 8000 ITTP = 1, NSPOBEX_TP
              ICA_OFF = 1
              ICB_OFF = 1 +  NGAS
              IAA_OFF = 1 +2*NGAS
              IAB_OFF = 1 +3*NGAS
#if defined (VAR_MPI2)
              IF( SPLIT_IJKL ) THEN
C               determine T type for use in T list
                IXXTTP = ITTP
                IF( (IRIS .eq. 1 .and. IRIC .eq. 2)
     &             .or. ( IRIS .eq. 2 .and. IRIC .eq. 1) ) THEN
                     IXXTTP = ITTP + NSPOBEX_TP
                END IF
C               active T block? (complex case included)
                ID_GET = IT_TTPL( IXXTTP )
csk             WRITE(LUWRT,*) ' type IXXTTP and ID_GET',IXXTTP,ID_GET
                IF( ID_GET .le. 0 ) THEN
                  ITB_BUFF_ADD = 0
                  GOTO 8000
                ELSE
C                 length of T block to compute T_BUFF offset
                  ITB_BUFF_ADD = NELM_CC(ITTP,NSPOBEX_TP,LEN_T_BUFF_NZ)
                END IF
              END IF
#endif
C
C             check connections 
              CALL GXFSTR(NELFSPGP(1,ILATP),NELFSPGP(1,ILBTP),
     &                    NELFSPGP(1,IRATP),NELFSPGP(1,IRBTP),
     &                    ITSPOBEX_TP(ICA_OFF,ITTP), 
     &                    ITSPOBEX_TP(IAA_OFF,ITTP),
     &                    ITSPOBEX_TP(ICB_OFF,ITTP), 
     &                    ITSPOBEX_TP(IAB_OFF,ITTP),
     &                    NGAS,ICON)
csk           WRITE(LUWRT,*) ' type ITTP and ICON',ITTP,ICON
              IF( ICON .eq. 1 .and. ISPOBEX_AC(ITTP) .eq. 1 ) THEN
C
C                offset for T array
                 ITB = IBTSOSO_TP(ITTP)
                 IF( (IRIS .eq. 1 .and. IRIC .eq. 2)
     &                .or. ( IRIS .eq. 2 .and. IRIC .eq. 1) ) THEN
                      ITB = ITB + LEN_T_VEC
                 END IF
C
C                ISIGDEN ==1
C
C                "shared memory" MPI mode:
C                copy T-block from window to T_SCR_CC
#if defined (VAR_MPI2)
                 IF( IT_SHL .ge. 0 )THEN
C
C                  get ICPL8 elements of T
                   ICPL_I8     = NELM_CC(ITTP,NSPOBEX_TP,LEN_T_BUFF_NZ)
                   XNINT_TP_I8 = REAL( ICPL_I8 )
                   ICPL        = IGIVE_I_B(XNINT_TP_I8)
C                  determine T type for use in T list
                   IXXTTP = ITTP
                   IF( (IRIS .eq. 1 .and. IRIC .eq. 2)
     &                .or. ( IRIS .eq. 2 .and. IRIC .eq. 1) ) THEN
                      IXXTTP = ITTP + NSPOBEX_TP
                   END IF
C                  process id and window offset (complex case included)
                   ID_GET = IT_TTPL( IXXTTP )
C                  check for availability
                   IF( ID_GET .lt. 0) GOTO 8000
C                  window offset
                   ITB_I8 = IT_TTOL( IXXTTP )
C
                   IF( ISIGDEN .eq. 1 )THEN
csk                  write(luwrt,*) ' ID_GET, ITB_I8, IXXTTP, ICPL',
csk  &                                ID_GET, ITB_I8, IXXTTP, ICPL
                     start_timer = MPI_WTIME() 
                     IF( MYNEW_ID_SM .eq. ID_GET )THEN
                       CALL DCOPY(ICPL,T_BUFF(ITB_I8+1),1,T_SCR_CC,1)
                     ELSE
                       CALL MPIXGET(T_SCR_CC,ICPL,2,ID_GET,
     &                              ITB_I8,ICPL,2,MY_T_WIN)
                     END IF
                     mem_cp_time = mem_cp_time + 
     &                             MPI_WTIME() - start_timer
                   ELSE IF( ISIGDEN .eq. 2 )THEN
                     CALL DZERO8( T_SCR_CC, ICPL_I8)
                   END IF
C                  ^ ISIGDEN
                 END IF
C                ^ IT_SHL 
#endif
C
C                read C block in
                 IF(IREAD.EQ.0) THEN
C                  number of blocks
                   NCBLK = NSMST
C
C                  fetch c-blocks
#if defined (VAR_MPI2)
                   IC_OFF_BLK     = INT_IOFF1_C
                   IOFF_SCRATCH_C = 0
#endif
                   DO JBLK = ICBLK, ICBLK + NCBLK - 1
                     JASM = ICBLOCK(3,JBLK)
                     JBSM = ICBLOCK(4,JBLK)
                     JOFF = ICBLOCK(5,JBLK)
                     JLEN = ICBLOCK(7,JBLK)
#if defined (VAR_MPI2)
CSK                  WRITE(LUWRT,*) ' checking for IC_OFF_BLK+JBLK-1',
CSK  &                                IC_OFF_BLK+JBLK-1
CSK                  WRITE(LUWRT,*) ' IOFF_SCR_C is', IOFF_SCR_C
CSK                  WRITE(LUWRT,*) ' JLEN is for IST_C_BLK_CT + JBLK', 
CSK  &                                JLEN, IST_C_BLK_CT + JBLK - 1
C
C                    check for active blocks
                     IF(ISCLFAC_GROUP(IC_OFF_BLK+JBLK-1).eq.0) GOTO 8888
CSK                  WRITE(LUWRT,*) ' Read-in for JBLK',JBLK,
CSK  &                                IC_OFF_BLK+JBLK-1
C
C                    new offset in MPI-2 parallel file
                     IOFF_SCR_C_EXPL = IOFF_SAVE + IOFF_SCR_C 
CSK                  WRITE(LUWRT,*) ' Read-in at file offset',
CSK  &                                IOFF_SCR_C_EXPL
CSK                  WRITE(LUWRT,*) ' memory offset',
CSK  &                                JOFF
C
C                    read block in core
                     xxxbitime = MPI_WTIME()
                     CALL GSTTBLD_PAR2(CB(JOFF),LLUC,
     &                                 IOFF_SCR_C_EXPL,JLEN)
                     xreadtimebi = xreadtimebi-xxxbitime+MPI_WTIME()
#else
                     CALL GSTTBLD2(CB(JOFF),LUC)
#endif
C
                     IBI_MULT_BL = IBI_MULT_BL + 1
C                    real/imaginary scaling according to:
C                    S_r = (H_r . C_r) - (H_i . C_i)
                     if (IRIS.eq.1.and.IRIC.eq.2) then
                       CALL DSCAL(JLEN,ONEM,CB(JOFF),1)
                     end if
C                    NTEST = 200
C                    IF (NTEST.GE.200) THEN
C                      LEN = ICBLOCK(7,JBLK)
C                      IF( LEN .gt. 0 ) THEN
C                        WRITE(LUWRT,*) ' C block read in, '//
C    &                                  'JBLK, JOFF, IRIS, IRIC:',
C    &                   JBLK,JOFF,IRIS,IRIC
C                        LEN = ICBLOCK(7,JBLK)
C                        CALL WRTMATMN(CB(JOFF),1,LEN,1,LEN,LUWRT)
C                      END IF
C                    END IF
C                    NTEST = 000
#if defined (VAR_MPI2)
 8888                CONTINUE
                     IOFF_SCR_C = IOFF_SCR_C + JLEN
#endif
                   END DO
                   IREAD = 1
                 END IF
C                ^ End if block should be read in 
C
C                setting undefined variables to zero:
                 ICRESTRICT = 0
                 ISRESTRICT = 0
                 ITDIAG = 0
C                proper multiplication sign
                 SIGN_EXTERNAL = SIGN_NHX(ITTP)
C
                 CALL GNSIDE_REL(ISIGDEN,IRIS,IRIC,
     &                           ITSPOBEX_TP(ICA_OFF,ITTP),
     &                           ITSPOBEX_TP(ICB_OFF,ITTP),
     &                           ITSPOBEX_TP(IAA_OFF,ITTP),
     &                           ITSPOBEX_TP(IAB_OFF,ITTP),
     &                           NELFSPGP(1,ILATP),NELFSPGP(1,ILBTP),
     &                           NELFSPGP(1,IRATP),NELFSPGP(1,IRBTP),
     &                           NSSOA_S(1,IATP),NSSOB_S(1,IBTP),
     &                           NSSOA_C(1,JATP),NSSOB_C(1,JBTP),
     &                           T(ITB),T_BUFF(ITB_BUFF),T_SCR_CC,ITTP,
     &                           SB(ISOFF),CB(ICOFF),ISSM,ICSM,
     &                           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                           XINT,CJRES,SIRES,MAXK,
     &                           ITOCC1,ITOCC2,ITOCC3,ITOCC4,
     &                           I1G,XI1G,I2G,XI2G,I3G,XI3G,I4G,XI4G,
     &                           I1GE,XI1GE,I2GE,XI2GE,I3GE,XI3GE,I4GE,
     &                           XI4GE,ITDIAG,ICRESTRICT,ISRESTRICT,
     &                           SIGN_EXTERNAL,IPRS,IPRD,ISMDST
#if defined (VAR_MPI2)
     &                          ,INT_IOFF1_C+ICBLK-1,IS_OFF_BLK+ISBLK-1,
     &                           ISCLFAC_GROUP,NPARBLOCK_S
#endif
     &           )
C
C                ISIGDEN ==2
C
C                "shared memory" MPI mode:
C                all : T_SCR_CC --> MY_T_WIN 
C
                 IF( IT_SHL .ge. 0 )THEN
C
                   IF( ISIGDEN .eq. 2 )THEN
C
                     XXX = DDOT( ICPL, T_SCR_CC, 1, T_SCR_CC, 1 )
C                    update if necessary
                     IF( XXX .ne. 0.0D0 ) THEN
C                       process id and window offset
C                       (complex case included)
#if defined (VAR_MPI2)
C                       use ID_GET as ID_ACCUM
                        ID_ACCUM = ID_GET
C                       window offset determined above
                        CALL MPIXACCUM(T_SCR_CC,ICPL,2,ID_ACCUM,
     &                                 ITB_I8,ICPL,2,MY_T_WIN)
#endif
                     END IF
                   END IF
C                  ^ ISIGDEN
                 END IF
C                ^ IT_SHL
               END IF
C              ^ end if connection 
#if defined (VAR_MPI2)
            ITB_BUFF = ITB_BUFF + ITB_BUFF_ADD
#endif
C
 8000       CONTINUE
C           ^  end of loop over spinorbital-excitation types 
C
#if defined LUCI_DEBUG
            NTEST = 300
            IF(NTEST.GE.300) THEN
              WRITE(LUWRT,*) ' Updated sigma blocks ' 
#if defined (VAR_MPI2)
              CALL WRTVH1_PAR(SB(ISOFF),ISSM,NSSOA_S(1,IATP),
     &                        NSSOB_S(1,IBTP),NSMST,0,NPARBLOCK_S,
     &                        IS_OFF_BLK+ISBLK-1,1)
#else
              CALL WRTVH1(SB(ISOFF),ISSM,NSSOA_S(1,IATP),
     &                    NSSOB_S(1,IBTP),NSMST,0)
#endif
            END IF
            NTEST = 000
#endif /* debug mode */
C
 9999       CONTINUE
C           ^ skip TT-blocks for sigma if not included in nodelist
10000     CONTINUE
C         ^ End of loop over blocks of sigma 
#if defined (VAR_MPI2)
          IF(IREAD.EQ.0) THEN
C
C           these C blocks were not needed, skip them
            DO JBLK = ICBLK, ICBLK + NCBLK - 1
CSK            WRITE(LUWRT,*) 'IREAD == 0: skipping for block',
CSK  &                         IST_C_BLK_CT+JBLK - 1
              IOFF_SCR_C = IOFF_SCR_C + IBLOCKL_C(IST_C_BLK_CT+JBLK-1)
            END DO
C
          END IF 
C         ^ unused c-blocks
 8999     CONTINUE
C         ^ skip TT-c-blocks if not 'active'
          IF( ICOMPUTE_C .eq. 0 ) THEN
C           skip c-blocks on file
            DO JBLK = ICBLK, ICBLK + NCBLK - 1
CSK           WRITE(LUWRT,*) 'skipping for block',IST_C_BLK_CT+JBLK - 1
              IOFF_SCR_C = IOFF_SCR_C +IBLOCKL_C(IST_C_BLK_CT+JBLK-1)
            END DO
          END IF
C
#else
          IF(IREAD.EQ.0) THEN
C           these C blocks were not needed, skip them
            DO JBLK = 1, NSMST
              CALL IFRMDS(LBL,-1,1,LUC)
              CALL SKPRCD2(LBL,-1,LUC)
            END DO
          END IF
#endif
 9000   CONTINUE
C       ^ End of loop over TT C blocks in batch
20000 CONTINUE
C     ^End of loop over batches of C blocks
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE TCC_SUBBLK(TCC,TCC_BUFF,T_SCR_CC,TCC_SUB,IWAY,
     &                      NCA,NCA_SUB,ICA_B,NCB,NCB_SUB,ICB_B,
     &                      NAA,NAA_SUB,IAA_B,NAB,NAB_SUB,IAB_B,
     &                      IAB_TRNSP,IDIAG,IPRNT)
*
* Extract - or add- subbblock of a TCC block
*
* IWAY = 1 Obtain TCC_SUB from TCC
* IWAY = 2 Add TCC_SUB to TCC
*
* Jeppe Olsen, Summer of 99
*              Updated July 2000, IAB_TRNSP, IDIAG added
*
* IAB_TRNSP = 1 : Input block is TCC(ICB,ICA,IAB,IAA)
* IDIAG     = 1 : Input block is lower diagonal
*
* Stefan Knecht, Nov. 2007 
*                use of T_SCR_CC in MPI "shared memory" mode added
*                --> treated in same way as TCC_SUB
*
* Stefan Knecht, May  2008 
*                use of TCC_BUFF in split mode added
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "parluci.h"
*. Input and output
      DIMENSION TCC(NCA*NCB*NAA*NAB), TCC_BUFF(*), T_SCR_CC(*)
      DIMENSION TCC_SUB(NCA_SUB,NCB_SUB,NAA_SUB,NAB_SUB)
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
csk      NTEST = 500
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.500) THEN
        WRITE(LUWRT,*) ' Input to TCC_SUBBLK '
        WRITE(LUWRT,*) ' NCA, NCB, NAA, NAB', NCA,NCB,NAA,NAB
        WRITE(LUWRT,*) ' NCA_SUB, NCB_SUB, NAA_SUB, NAB_SUB',
     &               NCA_SUB, NCB_SUB, NAA_SUB, NAB_SUB
        WRITE(LUWRT,*) ' IAB_TRNSP, IDIAG = ', IAB_TRNSP, IDIAG
        WRITE(LUWRT,*) ' Iway = ', IWAY
        IF( SHARED_M )THEN
          WRITE(LUWRT,*) ' T_SCR_CC input BLOCK'
          CALL WRTMATMN(T_SCR_CC,NCA*NCB,NAA*NAB,NCA*NCB,NAA*NAB,LUWRT)
        ELSE IF(SPLIT_IJKL)THEN
          WRITE(LUWRT,*) ' TCC_BUFF input BLOCK'
          CALL WRTMATMN(TCC_BUFF,NCA*NCB,NAA*NAB,NCA*NCB,NAA*NAB,LUWRT)
        ELSE
          WRITE(LUWRT,*) ' TCC input BLOCK'
          CALL WRTMATMN(TCC,NCA*NCB,NAA*NAB,NCA*NCB,NAA*NAB,LUWRT)
        END IF
      END IF
#endif
*
      IF(IAB_TRNSP.EQ.0.AND.IDIAG.EQ.0) THEN

       DO ICA = 1, NCA_SUB
         DO ICB = 1, NCB_SUB
           DO IAA = 1, NAA_SUB
             DO IAB = 1, NAB_SUB
               IADR_FULL = (IAB+IAB_B-2)*NAA*NCB*NCA + 
     &                     (IAA+IAA_B-2)*NCB*NCA +
     &                     (ICB+ICB_B-2)*NCA +
     &                      ICA+ICA_B-1
C?             WRITE(6,*) ' ICA, ICB, IAA, IAB =', ICA,ICB,IAA,IAB
csk               WRITE(LUWRT,*) ' IADR_FULL 1 is',IADR_FULL
               IF( SHARED_M )THEN
                 IF(IWAY.EQ.1) THEN
                   TCC_SUB(ICA,ICB,IAA,IAB) =  T_SCR_CC(IADR_FULL)
                 ELSE
                   T_SCR_CC(IADR_FULL) = T_SCR_CC(IADR_FULL) + 
     &                                   TCC_SUB(ICA,ICB,IAA,IAB)
                 END IF
               ELSE IF( SPLIT_IJKL )THEN
                 IF(IWAY.EQ.1) THEN
                   TCC_SUB(ICA,ICB,IAA,IAB) =  TCC_BUFF(IADR_FULL)
                 ELSE
                   TCC_BUFF(IADR_FULL) = TCC_BUFF(IADR_FULL) + 
     &                                   TCC_SUB(ICA,ICB,IAA,IAB)
                 END IF
               ELSE 
                 IF(IWAY.EQ.1) THEN
                   TCC_SUB(ICA,ICB,IAA,IAB) =  TCC(IADR_FULL)
                 ELSE
                   TCC(IADR_FULL) = TCC(IADR_FULL) + 
     &                              TCC_SUB(ICA,ICB,IAA,IAB)
                 END IF
               END IF
C              ^ "shared memory" / split mode
             END DO
           END DO
         END DO
       END DO
      ELSE IF (IAB_TRNSP.EQ.1.AND.IDIAG.EQ.0) THEN
       DO ICA = 1, NCA_SUB
       DO ICB = 1, NCB_SUB
       DO IAA = 1, NAA_SUB
       DO IAB = 1, NAB_SUB
*. Obtain T(ICA,ICB,IAA,IAB) as T(ICB,ICA,IAB,IAA)
         IADR_FULL = (IAA+IAA_B-2)*NAB*NCA*NCB + (IAB+IAB_B-2)*NCA*NCB
     +             + (ICA+ICA_B-2)*NCB + ICB+ICB_B-1
csk               WRITE(LUWRT,*) ' IADR_FULL 2 is',IADR_FULL
C?       WRITE(6,*) ' ICA, ICB, IAA, IAB =', ICA,ICB,IAA,IAB
         IF( SHARED_M )THEN
           IF(IWAY.EQ.1) THEN
             TCC_SUB(ICA,ICB,IAA,IAB) =  T_SCR_CC(IADR_FULL)
           ELSE
             T_SCR_CC(IADR_FULL)  =  T_SCR_CC(IADR_FULL) + 
     &                               TCC_SUB(ICA,ICB,IAA,IAB)
           END IF
         ELSE IF (SPLIT_IJKL) THEN
           IF(IWAY.EQ.1) THEN
             TCC_SUB(ICA,ICB,IAA,IAB) =  TCC_BUFF(IADR_FULL)
           ELSE
             TCC_BUFF(IADR_FULL)  =  TCC_BUFF(IADR_FULL) + 
     &                               TCC_SUB(ICA,ICB,IAA,IAB)
           END IF
         ELSE
           IF(IWAY.EQ.1) THEN
             TCC_SUB(ICA,ICB,IAA,IAB) =  TCC(IADR_FULL)
           ELSE
             TCC(IADR_FULL)  =  TCC(IADR_FULL) + 
     &                          TCC_SUB(ICA,ICB,IAA,IAB)
           END IF
         END IF
C        ^ "shared memory" / split mode
       END DO
       END DO
       END DO
       END DO
      ELSE IF (IDIAG.EQ.1) THEN
       DO ICA = 1, NCA_SUB
       DO ICB = 1, NCB_SUB
       DO IAA = 1, NAA_SUB
       DO IAB = 1, NAB_SUB
         ICA_ABS = ICA + ICA_B -1
         ICB_ABS = ICB + ICB_B -1
         IAA_ABS = IAA + IAA_B -1
         IAB_ABS = IAB + IAB_B -1
*
         IF((IAA_ABS.GT.IAB_ABS).OR.
     &      (IAA_ABS.EQ.IAB_ABS.AND.ICA_ABS.GE.ICB_ABS)) THEN
            IADR_FULL = ITDIANUM(ICA_ABS,ICB_ABS,IAA_ABS,IAB_ABS,
     &                  NCA,NAA)
            ILOW = 1
         ELSE 
            IADR_FULL = ITDIANUM(ICB_ABS,ICA_ABS,IAB_ABS,IAA_ABS,
     &                  NCA,NAA)
            ILOW = 0
         END IF
csk               WRITE(LUWRT,*) ' IADR_FULL 3 is',IADR_FULL
*
         IF( SHARED_M )THEN
           IF(IWAY.EQ.1) THEN
             TCC_SUB(ICA,ICB,IAA,IAB) =  T_SCR_CC(IADR_FULL)
           ELSE IF(ILOW.EQ.1) THEN
             T_SCR_CC(IADR_FULL)  =  T_SCR_CC(IADR_FULL) + 
     &                               TCC_SUB(ICA,ICB,IAA,IAB)
           END IF
         ELSE IF (SPLIT_IJKL) THEN
           IF(IWAY.EQ.1) THEN
             TCC_SUB(ICA,ICB,IAA,IAB) =  TCC_BUFF(IADR_FULL)
           ELSE IF(ILOW.EQ.1) THEN
             TCC_BUFF(IADR_FULL)  =  TCC_BUFF(IADR_FULL) + 
     &                               TCC_SUB(ICA,ICB,IAA,IAB)
           END IF
         ELSE
           IF(IWAY.EQ.1) THEN
             TCC_SUB(ICA,ICB,IAA,IAB) =  TCC(IADR_FULL)
           ELSE IF(ILOW.EQ.1) THEN
             TCC(IADR_FULL)  =  TCC(IADR_FULL) + 
     &                          TCC_SUB(ICA,ICB,IAA,IAB)
           END IF
         END IF
C        ^ "shared memory" / split mode
       END DO
       END DO
       END DO
       END DO
      END IF
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.500) THEN 
        WRITE(LUWRT,*) ' TCC_SUBBLK : TCC_SUB '
        NCAB_SUB = NCA_SUB*NCB_SUB
        NAAB_SUB = NAA_SUB*NAB_SUB
        CALL WRTMATMN(TCC_SUB,NCAB_SUB,NAAB_SUB,NCAB_SUB,NAAB_SUB,LUWRT)
      END IF
#endif
*  
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE WRT_TCC_BLK(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
*
*. Write TCC block containing all symmetries, total sym of T is ITCC_SM
*
* Jeppe Olsen, summer of 99
*
*  Adapted to double group symmetry
*    Timo Fleig, Apr. 07, 2003
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "symm.inc"
*
*. Input
       DIMENSION TCC(*)
       INTEGER NCA(*), NCB(*), NAA(*), NAB(*)
*
       WRITE(6,*) ' Block of Computed Density Matrices '
       WRITE(6,*) '===================================='
       IOFF = 1
       DO ISM_C = 1, NSMST
         ISM_A = IDBGMULT(ISM_C,IADJSYM(ITCC_SM)) 
         DO ISM_CA = 1, NSMST
           ISM_CB = IDBGMULT(ISM_C,IADJSYM(ISM_CA))
           DO ISM_AA = 1, NSMST
             ISM_AB = IDBGMULT(ISM_A,IADJSYM(ISM_AA))
C?           WRITE(6,*) ' ISM_AB = ', ISM_AB
             LCA = NCA(ISM_CA)
             LCB = NCB(ISM_CB)
             LAA = NAA(ISM_AA)
             LAB = NAB(ISM_AB)
             LENGTH = LCA*LCB*LAA*LAB
*
             IF(LENGTH.NE.0) THEN
               WRITE(6,'(A,4I4)') ' Sym of CA, CB, AA, AB ',
     &         ISM_CA, ISM_CB, ISM_AA, ISM_AB  
               CALL WRTMAT(TCC(IOFF),1,LENGTH,1,LENGTH)
             END IF
             IOFF = IOFF + LENGTH
           END DO
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
      SUBROUTINE WRT_TP_GENOP(ICA,ICB,IAA,IAB,NGAS)
*
* Write occupation in each GASspace for general operator
*
* Jeppe Olsen, August 1999
*
#include "implicit.inc"
*
      INTEGER ICA(NGAS),ICB(NGAS),IAA(NGAS),IAB(NGAS)
*
      WRITE(6,*) ' Occupation of alpha-creation string '
      CALL IWRTMA(ICA,1,NGAS,1,NGAS)
      WRITE(6,*) ' Occupation of beta -creation string '
      CALL IWRTMA(ICB,1,NGAS,1,NGAS)
      WRITE(6,*) ' Occupation of alpha-annihilation string '
      CALL IWRTMA(IAA,1,NGAS,1,NGAS)
      WRITE(6,*) ' Occupation of beta -annihilation string '
      CALL IWRTMA(IAB,1,NGAS,1,NGAS)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE WRT_SPOX_TP(IEX_TP,NEX_TP)
*
* Print types of spin-orbital excitations
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
      INTEGER IEX_TP(4*NGAS,NEX_TP)
*
      WRITE(6,*)
      WRITE(6,*) ' ***************************************** '
      WRITE(6,*) ' Information about spinorbital excitations '
      WRITE(6,*) ' ***************************************** '
      WRITE(6,*)
*
      DO JEX_TP = 1, NEX_TP
        WRITE(6,*)
        WRITE(6,*) ' Included spinorbitalexcitation ', JEX_TP
        WRITE(6,'(A,16I4)')
     &  ' Creation of alpha     :',
     &  (IEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)')
     &  ' Creation of beta      :',
     &  (IEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)')
     &  ' Annihilation of alpha :',
     &  (IEX_TP(I+2*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)')
     &  ' Annihilation of beta  :',
     &  (IEX_TP(I+3*NGAS,JEX_TP),I=1,NGAS)
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
*
* Write one-electron integrals with symmetry IVSM
* ISYM = 1 => Only lower triangular matrix included
*
* Jeppe Olsen, Jan. 1999
* Timo Fleig, June 2001: relativistic case
*
#include "implicit.inc"
#include "mxpdim.inc"
*. General input
#include "symm.inc"
      INTEGER NRPSM(NSMOB),NCPSM(NSMOB)
*. Specific input
      DIMENSION H(*)
*
      IOFF = 1
        DO ISM = 1, NSMOB
          JSM = IDBGMULT(INVELM(ISM),IHSM)
          NI = NRPSM(ISM)
          NJ = NCPSM(JSM)
          IF(ISYM.EQ.0.OR.ISM.GT.JSM) THEN
*. Complete block
            WRITE(6,*) ' Block with symmetry ISM, JSM ',ISM,JSM
            CALL WRTMAT(H(IOFF),NI,NJ,NI,NJ)
            IOFF = IOFF + NI*NJ
          ELSE IF (ISYM.EQ.1.AND.ISM.EQ.JSM) THEN
            CALL PRSYM(H(IOFF),NI)
            IOFF = IOFF * NI*(NI+1)/2
          END IF
        END DO
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WRTVH1_PAR(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM,IACTL,IOFFLST,
     &                      IVERS)
*
* Write one-electron integrals with symmetry IVSM
* ISYM = 1 => Only lower triangular matrix included
*
* Jeppe Olsen, Jan. 1999
* Timo Fleig, June 2001: relativistic case
* Stefan Knecht, May 2007: parallel adaption
*
#include "implicit.inc"
#include "mxpdim.inc"
*. General input
#include "symm.inc"
#include "parluci.h"
      INTEGER NRPSM(NSMOB),NCPSM(NSMOB), IACTL(*)
*. Specific input
      DIMENSION H(*)
*
      IOFF = 1
      IRUN = 0
      IRUN = IOFFLST
      IF( IVERS .eq. 0 ) THEN
        DO ISM = 1, NSMOB
          IF( IACTL( IRUN ) .eq. 0 ) GOTO 100
          JSM = IDBGMULT(INVELM(ISM),IHSM)
          NI = NRPSM(ISM)
          NJ = NCPSM(JSM)
          IF(ISYM.EQ.0.OR.ISM.GT.JSM) THEN
*. Complete block
            WRITE(LUWRT,*) ' Block with symmetry ISM, JSM ',ISM,JSM
            CALL WRTMATMN(H(IOFF),NI,NJ,NI,NJ,LUWRT)
            IOFF = IOFF + NI*NJ
          ELSE IF (ISYM.EQ.1.AND.ISM.EQ.JSM) THEN
            CALL PRSYM(H(IOFF),NI)
            IOFF = IOFF * NI*(NI+1)/2
          END IF
 100      CONTINUE
          IRUN = IRUN + 1
        END DO
      ELSE
        DO ISM = 1, NSMOB
          IF( IACTL( IRUN ) .ne. MYPROC ) GOTO 200
          JSM = IDBGMULT(INVELM(ISM),IHSM)
          NI = NRPSM(ISM)
          NJ = NCPSM(JSM)
          IF(ISYM.EQ.0.OR.ISM.GT.JSM) THEN
*. Complete block
            WRITE(LUWRT,*) ' Block with symmetry ISM, JSM ',ISM,JSM
            CALL WRTMATMN(H(IOFF),NI,NJ,NI,NJ,LUWRT)
            IOFF = IOFF + NI*NJ
          ELSE IF (ISYM.EQ.1.AND.ISM.EQ.JSM) THEN
            CALL PRSYM(H(IOFF),NI)
            IOFF = IOFF * NI*(NI+1)/2
          END IF
 200      CONTINUE
          IRUN = IRUN + 1
        END DO
      END IF
*
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XSIGDEN_CTRL(VEC1,VEC2,T_X,COPY_NEW_C,NVC_LUC,IZ
#if defined (VAR_MPI2)
     &                        ,LU1LIST,LU2LIST,LUCLIST,LBATX,LEBATX,
     &                        I1BATX,IBATX,NBATX,IBLOCKL_X,
     &                        NPARBLOCK_X,IGROUPLIST,IPROCLIST,
     &                        RCCTOS_X,IXSYMPAIRS,ICSM,ISSM,IX_CT,
     &                        IX_ST,IOFF_LU1
#endif
     &                        )
C**********************************************************************
C
C     sigma vector calculation wrapper to be used in CI property runs.
C
C     written by S. Knecht - Nov 2008
C
C     last revision:
C
C**********************************************************************
#include "ipoist8.inc"
#include "krciprop.h"
#include "parluci.h"
#if defined (VAR_MPI2)
#include "mpif.h"
#include "infpar.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE), IBLOCKL_X(NUM_BLOCKS,*)
      DIMENSION NPARBLOCK_X(NUM_BLOCKS,*) 
      DIMENSION IXSYMPAIRS(MXPROPKRCI_SYM,*), LU1LIST(*), LU2LIST(*)
      DIMENSION LUCLIST(*), LBATX(*), LEBATX(*), I1BATX(*), IBATX(8,*)
      INTEGER RCCTOS_X(NUM_BLOCKS2,*)
#endif
      DIMENSION VEC1(*), VEC2(*), T_X(*)
#include "clunit.inc"
#include "dummy.h"
      LOGICAL COPY_NEW_C
C
C
      IF( .not. COPY_NEW_C ) GOTO 100
#if defined (VAR_MPI2)
C
C     erase old LUCLIST
      CALL IZERO(LUCLIST,IALL_LUC)
C
!     WRITE(LUWRT,*) ' IOFF_LU1 is',IOFF_LU1
      IF( IZ .eq. 1 ) THEN
C       real quaternion double group
        CALL COPVCD_PP_CC_B(ILU1,ILUC,VEC1,NBATX,
     &                      LBATX,LEBATX,I1BATX,IBATX,
     &                      MY_LU1_OFF,MY_LUC_OFF,
     &                      LU1LIST,LUCLIST,
     &                      IBLOCKL_X(1,IX_CT),IOFF_LU1)
      ELSE
C       complex quaternion double group
        CALL COPVCD_PP_CC_B_C(ILU1,ILUC,VEC1,NBATX,
     &                        LBATX,LEBATX,I1BATX,IBATX,
     &                        MY_LU1_OFF,MY_LUC_OFF,
     &                        LU1LIST,LUCLIST,
     &                        IBLOCKL_X(1,IX_CT),IOFF_LU1)
      END IF
      COPY_NEW_C = .FALSE.
C     find corresponding symmetry pair
 100  IXPAIR = 0
      IF( IXSYMPAIRS(ICSM,ISSM) .gt. 0 ) THEN 
        IXPAIR = IXSYMPAIRS(ICSM,ISSM)
      ELSE
        WRITE(LUWRT,'(/A,2I4)')' *** ERROR in XSIGDEN_CTRL. no pair'//
     &                        ' type found for ICSM and ISSM:',ICSM,ISSM
        CALL QUIT(' *** ERROR in XSIGDEN_CTRL: unknown pair type.')
      END IF
C     erase old LU2LIST
      CALL IZERO(LU2LIST,IALL_LU2)
#else
C     copy new C vector to scratch file
      CALL REWINE(LUSC61,-1)
      CALL SKPVCD(LUC,IZ*(NVC_LUC-1),VEC1,1,-1)
      CALL COPVCDC(LUC,LUSC61,VEC1,0,IZ,-1)
      COPY_NEW_C = .FALSE.
 100  CALL REWINE(LUSC61,-1)
      CALL REWINE(LUHC,-1)
#endif
C     set offset for sigma-file
      JVEC_SF = 0
C
C     start calculation: sigma = T_X x C
C     ----------------------------------
      CALL SIGDEN_CTRL(VEC1,VEC2,
#if defined (VAR_MPI2)
     &                 ILUC,ILU2,
#else
     &                 LUSC61,LUHC,
#endif
     &                 T_X,DUMMY,1
#if defined (VAR_MPI2)
     &                ,LUCLIST,LU2LIST,IBLOCKL_X(1,IX_CT),
     &                 NPARBLOCK_X(1,IX_CT),IGROUPLIST,IPROCLIST,
     &                 RCCTOS_X(1,IXPAIR),DUMMY,DUMMY,
     &                 IBLOCKL_X(1,IX_ST),NPARBLOCK_X(1,IX_ST)
#endif
     &                 )
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST,IDIAG)
*
* Offsets for symmetryblocks of TCC elements, sym of CA,CB,AA used 
* If Idiag.eq.1 only blocks (Icasm,Iaasm).ge.(Icbsm,Iabsm) are included
*
*
* Jeppe Olsen, Summer of 99
*              July 2000, HNIE : IDIAG added
* Timo Fleig, summer of 2001, modified for relativistic case
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "symm.inc"
*. Input
*. Number of strings of different types, assumed to be 
*. creations strings
      INTEGER NCA(*),NCB(*),NAA(*),NAB(*)
*. Output
      INTEGER IBT(MXNDGIRR,MXNDGIRR,MXNDGIRR)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Z_TCC_OFF speaking '
        WRITE(6,*) ' NSMST, IDIAG = ', NSMST, IDIAG
      END IF
#endif
*
      IOFF = 1
      DO ISM_C = 1, NSMST
*. Symmetry of annihilation part
        ISM_A = IDBGMULT(ITSYM,ISM_C)
        DO ISM_CA = 1, NSMST
          ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
          DO ISM_AA = 1, NSMST
            ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))
*
            ISM_ALPHA = (ISM_AA-1)*NSMST+ISM_CA
            ISM_BETA  = (ISM_AB-1)*NSMST+ISM_CB
*
            IF((IDIAG.EQ.0).OR.
     &         (IDIAG.EQ.1.AND.ISM_ALPHA.GT.ISM_BETA)) THEN
              IBT(ISM_CA,ISM_CB,ISM_AA) = IOFF
              IF(IDIAG.EQ.1.AND.ISM_ALPHA.GT.ISM_BETA) 
     &        IBT(ISM_CB,ISM_CA,ISM_AB) = -IOFF
              IOFF = IOFF + 
     &          NCA(ISM_CA)*NCB(ISM_CB)*NAA(ISM_AA)*NAB(ISM_AB)
*                 ^ NAA was number of creation strings per sym
            ELSE IF (IDIAG.EQ.1.AND.ISM_ALPHA.EQ.ISM_BETA) THEN
              IBT(ISM_CA,ISM_CB,ISM_AA) = IOFF
              LEN = NCA(ISM_CA)*NAA(ISM_AA)
              IOFF = IOFF + LEN*(LEN+1)/2
            END IF
#if defined LUCI_DEBUG
            if (NTEST.ge.100) then
              write(6,'(A,I2,A,I2,A,I2,A,I5)') 
     &                'IBT(',ISM_CA,',',ISM_CB,',',ISM_AA,') =',
     &                 IBT(ISM_CA,ISM_CB,ISM_AA)
            end if
#endif
*
          END DO
        END DO
      END DO
*
      RETURN
      END
