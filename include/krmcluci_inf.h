C  FILE : krmcluci_inf.h
C  Information used for KRMC and LUCIAREL interface
C  ------------------------------------------------
C
C  1) CI energy
C  2) CI vector partitioning information.
C  3) IRIQ for real/imaginary mismatch when running real groups
C     in LUCIAREL (preliminary implementation)
C  4) Some auxiliary variables
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
      COMMON/CIVCINF/ECI,
     &               IBATS_MC,LBATS_MC,I1BATS_MC,LEBATS_MC,NBAT_MS2_MC,
     &               NMK2VAL,
     &               IRIQ,
     &               NIRR_DG_LOC
