!-----------------------------------------------------------------------------
!
! IPRCAP     - Print level in CAP related routines
! IPCAPR     - Pointer to real part of CAP operator
! IPCAPI     - Pointer to imaginary part of CAP operator
! IPCAP_S    - Pointer for overlap integrals from CAP (to be compared with DIRAC)
!
!-----------------------------------------------------------------------------

! ... flag for comparison between DIRAC-CAP overlap integrals
      LOGICAL CAPOVLP,NOOVLP
      COMMON / DCCAPOVL / CAPOVLP,NOOVLP

!     INTEGER*4  IPRCAP 
      INTEGER  IPRCAP, IPCAPR, IPCAPI, IPCAP_S
      COMMON / DCCAPI / IPRCAP, IPCAPR, IPCAPI, IPCAP_S

