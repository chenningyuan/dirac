!-----------------------------------------------------------------------------
!
! CAPOVLP - do only comparison between CAP and DIRAC overlap integrals
! NOOVLP - no comparison of CAP-DIRAC overlap integrals before calculation of CAP 
!
! DONRCAP - calculate only non-relativistic (Sajeev's original) version of the CAP
! otherwise calculate four-component relativistic CAP
!
! IPRCAP     - Print level in CAP related routines
! IPCAPR     - Pointer to real part of CAP operator
! IPCAPI     - Pointer to imaginary part of CAP operator
! IPCAP_S    - Pointer for overlap integrals from CAP (to be compared with DIRAC)
!
! IPCAP_DXR,IPCAP_DXI,IPCAP_DYR ... pointers for the real and imaginary parts of the Dirac Ham.components
!
! LUPEMAT    - lu number of the 'PEMAT' file with all CAP integrals
!
!-----------------------------------------------------------------------------

! ... flag for comparison between DIRAC-CAP overlap integrals
      LOGICAL CAPOVLP, NOOVLP, DONRCAP
      COMMON / DCCAPOVL / CAPOVLP, NOOVLP, DONRCAP

      INTEGER  IPRCAP, IPCAPR, IPCAPI,                                  &
     &         IPCAP_DR, IPCAP_DI, IPCAP_S, LUPEMAT,                    &
     &         IPCAP_DXR, IPCAP_DXI, IPCAP_DYR, IPCAP_DYI,              &
     &         IPCAP_DZR, IPCAP_DZI 

      COMMON / DCCAPI / IPRCAP, IPCAPR, IPCAPI,                         &
     &         IPCAP_DR, IPCAP_DI, IPCAP_S, LUPEMAT,                    &
     &         IPCAP_DXR, IPCAP_DXI, IPCAP_DYR, IPCAP_DYI,              &
     &         IPCAP_DZR, IPCAP_DZI
 
