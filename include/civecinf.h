C  CI vector partitioning information.
C  Used in KRMC - LUCIAREL communication.
C
C  Max. number of type-type-symmetry-symmetry blocks allowed:
C  (dynamically allocated in LUCIAREL)
      PARAMETER (MXNTTS_MC = 10000)
C
C  Max. no. of MK projection values allowed (same in LUCIAREL).
      PARAMETER (MXMKSTA = 100)
C
      DIMENSION IBATS_MC(8,MXNTTS_MC),LBATS_MC(MXNTTS_MC),
     &          I1BATS_MC(MXNTTS_MC),
     &          LEBATS_MC(MXNTTS_MC),NBAT_MS2_MC(MXMKSTA)
C
      COMMON/CIVCINF/IBATS_MC,LBATS_MC,I1BATS_MC,LEBATS_MC,NBAT_MS2_MC,
     &               NMK2VAL
