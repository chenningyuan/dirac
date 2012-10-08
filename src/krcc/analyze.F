      SUBROUTINE ANA_GENCC_CC_KRCC(T,ISM,
     &                             WORK,KFREE,LFREE)
*
* Analyze T-coefficients of general coupled cluster wavefunction
*
* Jeppe Olsen, September 1999 ( modified GASANA)
*
*
#include "implicit.inc"
*. specific input Coefficients and types
      DIMENSION T(*)
*
* =====
*.Input
* =====
*
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "ctcc.inc"
#include "ctccp.inc"
*
      DIMENSION WORK(*)
*
      CALL QENTER('ANAGC')
*
      NTEST = 0
      IPRNCIV = NTEST
*
*. Four blocks of string occupations
      CALL MEMGET('INTE',KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
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
C     LENGTH_SOX = NSPOBEX_TP_CC*10
      LENGTH_SOX = NSPOBEX_TPE*10
      CALL MEMGET('INTE',KNCPMT_SOX,LENGTH_SOX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KWCPMT_SOX,LENGTH_SOX,WORK,KFREE,LFREE)
*. Number and weight of orbital excitation type
      LENGTH_OX = (NOBEX_TP+1)*10
      CALL MEMGET('INTE',KNCPMT_OX,LENGTH_OX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KWCPMT_OX,LENGTH_OX,WORK,KFREE,LFREE)
      IUSLAB = 0
*
*. Occupation of strings of given sym and supergroup
      CALL ANA_GENCCS_CC_KRCC(T,WORK(KLSOBEX_CC),
     &            NSPOBEX_TPE,NOBEX_TP+1,ISM,
     &            THRES,MAXTRM,
     &            WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
     &            WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),
     &            IUSLAB,IDUMMY,
     &            WORK(KNCPMT_SOX),WORK(KWCPMT_SOX),
     &            WORK(KNCPMT_OX),WORK(KWCPMT_OX),
     &            NTOOB,IPRNCIV,WORK(KLSOX_TO_OX),
     &            WORK,KFREE,LFREE)

      CALL MEMREL('ANALYZE',WORK,KLSTR1_OCC,KLSTR1_OCC,KFREE,LFREE)
      CALL QEXIT('ANAGC')
*
      RETURN
      END
*
      SUBROUTINE ANA_GENCCS_CC_KRCC(T,ISPOBEX_TP,
     &                 NSPOBEX_TP,NOBEX_TP,ISM,
     &                 THRES,MAXTRM,
     &                 IOCC_CA, IOCC_CB, IOCC_AA, IOCC_AB,
     &                 IUSLAB,IOBLAB,
     &                 NCPMT_SOX,WCPMT_SOX,
     &                 NCPMT_OX, WCPMT_OX,
     &                 NORB,IPRNCIV,ISOX_TO_OX,
     &                 WORK,KFREE,LFREE)
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
*

*. If IUSLAB  differs from zero Character*6 array IOBLAB is used to identify
*  Orbitals
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "symm.inc"
*
      DIMENSION WORK(*)
*. Specific input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
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
C     DO I=10,25
C       T(I) = 0.8D0/I
C     END DO
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
        CALL WRT_SPOX_TP_CC_KRCC(ISPOBEX_TP(1,ITSS),1)
        WRITE(6,*) ' NSPOBEX_TP,ITSS = ', NSPOBEX_TP,ITSS
*. Transform from occupations to groups
       CALL OCC_TO_GRP_CC_KRCC(ISPOBEX_TP(1+0*NGAS,ITSS),IGRP_CA,1   )
       CALL OCC_TO_GRP_CC_KRCC(ISPOBEX_TP(1+1*NGAS,ITSS),IGRP_CB,1   )
       CALL OCC_TO_GRP_CC_KRCC(ISPOBEX_TP(1+2*NGAS,ITSS),IGRP_AA,1   )
       CALL OCC_TO_GRP_CC_KRCC(ISPOBEX_TP(1+3*NGAS,ITSS),IGRP_AB,1   )
*
       NEL_CA = IELSUM(ISPOBEX_TP(1+0*NGAS,ITSS),NGAS)
       NEL_CB = IELSUM(ISPOBEX_TP(1+1*NGAS,ITSS),NGAS)
       NEL_AA = IELSUM(ISPOBEX_TP(1+2*NGAS,ITSS),NGAS)
       NEL_AB = IELSUM(ISPOBEX_TP(1+3*NGAS,ITSS),NGAS)
*
      CALL MEMCHK_KRCC(WORK)
       DO ISM_C = 1, NSMST
        ISM_A =  IDBGMULT(ISM,INVELM(ISM_C))
        ISM_A = IADJSYM(ISM_A) !testing
        DO ISM_CA = 1, NSMST
         ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
         DO ISM_AAA = 1, NSMST !testing
          ISM_AA = IADJSYM(ISM_AAA) !testing
          ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))
C check ISM_ALPHA and ISM_BETA since these are CA A or B strings sym
          ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
          ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
*. obtain strings
          IUB = 1
          CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CA,NGAS,
     &         ISM_CA,NEL_CA,NSTR_CA,
     &         IOCC_CA, NORB,0,IDUM,IDUM,
     &         WORK,KFREE,LFREE)
          IUB = 2
          CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CB,NGAS,
     &         ISM_CB,NEL_CB,NSTR_CB,
     &         IOCC_CB, NORB,0,IDUM,IDUM,
     &         WORK,KFREE,LFREE)
          IUB = 1
          CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AA,NGAS,
     &         ISM_AA,NEL_AA,NSTR_AA,
     &         IOCC_AA, NORB,0,IDUM,IDUM,
     &         WORK,KFREE,LFREE)
          IUB = 2
          CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AB,NGAS,
     &         ISM_AB,NEL_AB,NSTR_AB,
     &         IOCC_AB, NORB,0,IDUM,IDUM,
     &         WORK,KFREE,LFREE)
*. Loop over T elements as  matrix T(I_CA, I_CB, IAA, I_AB)
          NSTR_TOT =NSTR_AB*NSTR_AA*NSTR_CB*NSTR_CA
          IF(NSTR_TOT.GE.1) THEN
C         print*,'NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA',
C    &            NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA
          END IF
          DO I_AB = 1, NSTR_AB
             I_AA_MIN = 1
           DO I_AA = I_AA_MIN, NSTR_AA
            DO I_CB = 1, NSTR_CB
               I_CA_MIN = 1
             DO I_CA = I_CA_MIN, NSTR_CA
              IT = IT + 1
              WRITE(6,*) ' IT, T(IT) = ', IT,T(IT)
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
C             IF( XMAX .GE. ABS(T(IT)) .AND.
C    &        ABS(T(IT)).GT. XMIN .AND. ITRM.LE.MAXTRM) THEN
                ITRM = ITRM + 1
                IIT = IIT + 1
                TTNORM = TTNORM + T(IT) ** 2
                WRITE(6,'(A)')
                WRITE(6,'(A)')
     &          '                 =================== '
                WRITE(6,*)

C      print*,'ISM_CA,ISM_CB,ISM_AA,ISM_AB',ISM_CA,ISM_CB,ISM_AA,ISM_AB
                WRITE(6,'(A,2I8,2X,E14.8)')
     &          '  Type, number, size of amplitude : ',
     &          ITSS, IT, T(IT)
C               IF(IUSLAB.EQ.0) THEN
                  WRITE(6,*) ' Creator alpha '
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_CA(IEL+(I_CA-1)*NEL_CA),IEL = 1, NEL_CA)
                  WRITE(6,*) ' Creator beta '
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_CB(IEL+(I_CB-1)*NEL_CB),IEL = 1, NEL_CB)
                  WRITE(6,*) ' Annihilator alpha '
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_AA(IEL+(I_AA-1)*NEL_AA),IEL = 1, NEL_AA)
                  WRITE(6,*) ' Annihilator beta '
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_AB(IEL+(I_AB-1)*NEL_AB),IEL = 1, NEL_AB)
C               END IF
C             END IF
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
        WRITE(6,'(A,I2,A,I2,3X,I9,X,E15.8,3X,E15.8)')
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
            WRITE(6,'(1X,I9,3X,E9.4,7X,E9.4,2X,I3))')
     &      N,W,WACC,ITSS
          END IF
      END DO
*
      RETURN
      END
