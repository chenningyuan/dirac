!
! VERSION : $Revision$
! DATE    : $Date$
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
