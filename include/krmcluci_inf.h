!  FILE : krmcluci_inf.h
!  Information used for KRMC and LUCIAREL interface
!  ------------------------------------------------
!
!  1) CI energy
!  2) CI vector partitioning information.
!  3) IRIQ for real/imaginary mismatch when running real groups
!     in LUCIAREL (preliminary implementation)
!  4) Some auxiliary variables
!
!  Max. number of type-type-symmetry-symmetry blocks allowed:
!  (dynamically allocated in LUCIAREL)
      INTEGER MXNTTS_MC
      PARAMETER (MXNTTS_MC = 10000)

!  Max. no. of MK projection values allowed (same in LUCIAREL).
      INTEGER MXMKSTA
      PARAMETER (MXMKSTA = 100)

      INTEGER IBATS_MC(8, MXNTTS_MC), LBATS_MC(MXNTTS_MC),              &
     &        I1BATS_MC(MXNTTS_MC), LEBATS_MC(MXNTTS_MC),               &
     &        NBAT_MS2_MC(MXMKSTA)

      REAL*8 ECI
      INTEGER NMK2VAL, IRIQ, NIRR_DG_LOC
      COMMON /CIVCINF/ ECI, IBATS_MC, LBATS_MC, I1BATS_MC, LEBATS_MC,   &
     &                 NBAT_MS2_MC, NMK2VAL, IRIQ, NIRR_DG_LOC
