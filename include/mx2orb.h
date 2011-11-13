!
! VERSION : $Revision: 4924 $
! DATE    : $Date: 2006-07-12 09:17:34 +0200 (Wed, 12 Jul 2006) $
! FILE    : maxorb.h
!
!     MAXORB = maximum number of orbitals
!
      INTEGER MAXORB
#if ! defined (LARGE_MOLECULES)
      PARAMETER (MAXORB = 1200)
#else
      PARAMETER (MAXORB = 3000)
#endif 
