      SUBROUTINE CLASS_TRUNC_REL(NCLS,ICLS_L,CLS_CT,CLS_ET,CLS_C,CLS_E,
     &                       E_CONV,ICLS_A,N_TRN_CLS,E_TRUNC,W_TRUNC)
*
* Decide which classes of parameters that can be eliminated
*
* Jeppe Olsen, Jan 97
*              March '97 updated
*
*. Note in current version all energy contributions are
*. positive.
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ICLS_L(NCLS)
      DIMENSION CLS_CT(NCLS),CLS_ET(NCLS),CLS_C(NCLS),CLS_E(NCLS)
*. Output
      INTEGER ICLS_A(NCLS)
*. Giving  the truncated classes
*  Additional output is
*              N_TRN_CLS : Number of truncated classes
*              E_TRUNC   : Estimated error of eliminating these classes
*                          (compared to expansion with only largest
*                           coefficients included).
*                          E_TRUNC is thus the error that should be
*                          added to the estmated error arising from
*                          eliminating small terms
*
      WRITE(6,*)
      WRITE(6,*) ' Welcome to TRUNC_CLASS '
      WRITE(6,*) ' ======================='
CM    WRITE(6,*)
CM    WRITE(6,*) ' Information about classes '
CM    WRITE(6,*) ' =========================='
CM    WRITE(6,*)
CM    WRITE(6,*)
CM   & ' Class    Dimension      E          E(Trunc)      C',
CM   & '         C(Trunc)'
CM    WRITE(6,*)
CM   & ' =============================================================',
CM   & '========='
CM     DO JCLS = 1, NCLS
CM       IF(ABS(CLS_ET(JCLS)).GT.0.0D0.OR.ABS(CLS_CT(JCLS)).GT.0.0D0)
CM   &   THEN
CM         WRITE(6,'(1H ,2X,I4,1X,I10,3X,4(E12.6,1X))')
CM   &     JCLS,ICLS_L(JCLS),CLS_ET(JCLS),CLS_E(JCLS),
CM   &     CLS_CT(JCLS),CLS_C(JCLS)
CM       END IF
CM     END DO
*
       WRITE(6,*)
     & ' Required threshold for convergence of energy',E_CONV
*
*. First find total energy correction and
*. largest class contribution
*
      E_TOT = 0.0D0
      E_CLS_MAX= 0.0D0
      DO JCLS = 1, NCLS
       E_TOT = E_TOT + CLS_ET(JCLS)
       E_CLS_MAX = MAX( E_CLS_MAX,CLS_ET(JCLS))
      END DO
      WRITE(6,*) ' E_CLS_MAX = ', E_CLS_MAX
*
*. The truncation is done in two steps, first an overall
*  threshold for deleting classes is constructed,
*  and all classes with energy contributions less than this
*  threshold are excluded.
*. Then one takes a pass through remaining classes and
*  eliminates untill the quota is filled
*
*. Loop over thesholds, increase by a given amount (1.5) in each it
*. Until total energy threshold is passed
*
* Another strategy : start with E_CLS_MAX and decrease from there
      ILOOP = 0
      STEP = 1.2
 1000 CONTINUE
        ILOOP = ILOOP+1
        IF(ILOOP.EQ.1) THEN
C         E_THRES = E_CONV/NCLS
          E_THRES =  E_CLS_MAX
        ELSE
C         E_THRES = E_THRES*STEP
          E_THRES = E_THRES/STEP
        END IF
        E_ELI = 0.0D0
        N_ELI = 0
        DO JCLS = 1, NCLS
          IF(ABS(CLS_ET(JCLS)).LE.E_THRES) THEN
            E_ELI = E_ELI + CLS_ET(JCLS)
            N_ELI = N_ELI + 1
          END IF
        END DO
        IF(ILOOP.EQ.10000) THEN
          WRITE(6,*) ' Loop count exceeded 10000'
          WRITE(6,*) ' I am afraid I am in an infinite loop'
          WRITE(6,*) ' So I will stop '
          Call Abend2( 'CLASS_TRUNC: Iloop.eq.10000' )
        END If
      IF(E_ELI.GT. E_CONV) GOTO 1000
*. We stepped out when the energy of the eliminated classes
*  became to large, go back to previous threshold
COLD  E_THRES = E_THRES/STEP
      WRITE(6,*)
     & ' Overall Threshold for eliminating classes',E_THRES
*
      N_PAS_CLS = 0
      L_PAS_CLS = 0
      N_TRN_CLS = 0
      L_TRN_CLS = 0
      N_ACT_CLS = 0
      L_ACT_CLS = 0
*
      E_TRUNC = 0.0D0
      E_TRUNCT = 0.0D0
      W_TRUNCT= 0.0D0
      W_TRUNC = 0.0D0
*
*. Eliminate classes with energy contribution less than E_THRES
*
      IONE = 1
      CALL ISETVC(ICLS_A,IONE,NCLS)
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).NE.0.0D0.AND.CLS_ET(JCLS).LE. E_THRES) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        END IF
      END DO
      WRITE(6,*)
     &  ' Number of classes with contributions less than E_THRES',
     &   N_TRN_CLS
      WRITE(6,*) ' Energy contributions from these classes ',
     &   E_TRUNCT
* Eliminate remaining classes until thres hold is obtained
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).EQ.0.0D0) THEN
*. Passive class, no contribution before truncation
          N_PAS_CLS = N_PAS_CLS + 1
          L_PAS_CLS = L_PAS_CLS + ICLS_L(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1
     &  .AND.ABS(E_TRUNCT+CLS_ET(JCLS)).LT.E_CONV) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1) THEN
*. Class is active
          N_ACT_CLS = N_ACT_CLS + 1
          L_ACT_CLS = L_ACT_CLS + ICLS_L(JCLS)
           ICLS_A(JCLS) = 1
        END IF
      END DO
*. Correct for missing minus in first order correction
      E_TRUNC  = -E_TRUNC
      E_TRUNCT = -E_TRUNCT
*
C?    WRITE(6,*) ' ICLS_A giving active classes '
C?    WRITE(6,*) ' ============================='
C?    WRITE(6,*)
C?    CALL IWRTMA(ICLS_A,1,NCLS,1,NCLS)
      WRITE(6,*)
      WRITE(6,*) ' Estimated complete energy contribution of ',
     &            '  eliminated classes ', E_TRUNCT
      WRITE(6,*) ' Estimated truncated energy contribution of ',
     &            '  eliminated classes ', E_TRUNC
      WRITE(6,*) ' Truncation error in eliminated classes was ',
     &              E_TRUNC-E_TRUNCT
      WRITE(6,*) ' Estimated truncated weight contribution of ',
     &            '  eliminated classes ', W_TRUNC
      WRITE(6,*) ' Estimated energy contribution without trunc',
     &             -E_TOT
      WRITE(6,*) ' Energy contribution of of active classes ',
     &             -E_TOT-E_TRUNCT

      WRITE(6,*)
      WRITE(6,*) '  Class      Number     Dimension  '
      WRITE(6,*) ' ================================= '
      WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Passive  ', N_PAS_CLS,L_PAS_CLS
      WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Truncated', N_TRN_CLS,L_TRN_CLS
      WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Active   ', N_ACT_CLS,L_ACT_CLS
*
      WRITE(6,*)
      WRITE(6,*) ' Information about classes '
      WRITE(6,*) ' =========================='
      WRITE(6,*)
      WRITE(6,*)
     & ' Class    Dimension      E          E(Trunc)      C',
     & '         C(Trunc)   Active'
      WRITE(6,*)
     & ' =============================================================',
     & '================'
       DO JCLS = 1, NCLS
         IF(ABS(CLS_ET(JCLS)).GT.0.0D0.OR.ABS(CLS_CT(JCLS)).GT.0.0D0)
     &   THEN
           WRITE(6,'(1H ,2X,I4,1X,I10,3X,4(E12.6,1X),1X,I2)')
     &     JCLS,ICLS_L(JCLS),CLS_ET(JCLS),CLS_E(JCLS),
     &     CLS_CT(JCLS),CLS_C(JCLS),ICLS_A(JCLS)
         END IF
       END DO
      RETURN
      END
