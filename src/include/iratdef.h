!     IRAT  = (real word length) / (integer word length)
!     IRAT2 = (real word length) / (half-integer word length)
!             if available and used, otherwise IRAT2 = IRAT
!     LRAT  = (real word length) / (logical word length)
      INTEGER IRAT, IRAT2, LRAT
#if defined (INT_STAR8)
      PARAMETER (IRAT = 1, LRAT = 1)
#else
      PARAMETER (IRAT = 2, LRAT = 2)
#endif
#if defined (VAR_STAR2)
      PARAMETER (IRAT2 = 4)
#else
      PARAMETER (IRAT2 = IRAT)
#endif
