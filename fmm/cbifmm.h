C
C variables required in a Multipole Method (FMM) run
C external to the main MM code
C (e.g. multipole moment generation order via HERMIT)
C
      LOGICAL         MMSKIP, MMCUT,
     &                MMSKIPDF, MMCUTDF
      INTEGER         IPRINT, LMAX, MMPRDF, LMAXDF
C
      COMMON /CBIFMM/ IPRINT, LMAX, MMSKIP, MMCUT
C
C     defaults:
C
      PARAMETER ( MMSKIPDF = .FALSE., MMCUTDF = .TRUE. )
      PARAMETER ( MMPRDF = 0, LMAXDF = 8 )

