#include <linesep.h>
!     IRAT  = (real word length) / (integer word length)
!     IRAT2 = (real word length) / (half-integer word length)
!             if available and used, otherwise IRAT2 = IRAT
!     LRAT  = (real word length) / (logical word length)
#if defined (SYS_CRAY) || defined (SYS_T3D)
      INTEGER IRAT, IRAT2, LRAT
      PARAMETER (IRAT = 1, IRAT2 = 1, LRAT = 1)

#else
#if defined (VAR_STAR2)
      PARAMETER (IRAT = 2, IRAT2 = 4, LRAT = 2)

#else
      PARAMETER (IRAT = 2, IRAT2 = 2, LRAT = 2)

#endif
#endif
