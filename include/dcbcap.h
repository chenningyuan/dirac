!-----------------------------------------------------------------------------
!
! CAPOVLP - do only comparison between CAP and DIRAC overlap integrals
! NOOVLP - no comparison between CAP and DIRAC overlap integrals
!           before calculating true CAP elements 
!
! DONRCAP - calculate only non-relativistic (Sajeev's original) version of the CAP
!            otherwise calculate four-component relativistic CAP
!
! NOCAPRE - in MO transformation step do NOT add real part of the CAP operator
! NOCAPIM - in the moltra step do NOT add imaginary part of the CAP operator
!
! CAPIMRE - in the moltra step add imaginary part of the CAP operator as REAL
!
! IPRCAP     - Print level in CAP related routines
!
! IPCAPR     - Pointer to real part of the non-relativistic CAP operator
! IPCAPI     - Pointer to imaginary part of the non-relativistic CAP operator
!
! IPCAP_DR   - Pointer to real part of the relativistic CAP operator
! IPCAP_DI   - Pointer to imaginary part of the relativistic CAP operator
!
! IPCAP_S    - Pointer for overlap integrals from CAP (to be compared with DIRAC)
!
!
! LUPEMAT    - lu number of the 'PEMAT' file with all CAP integrals
!
!-----------------------------------------------------------------------------

! ... flag for comparison between DIRAC-CAP overlap integrals
      LOGICAL CAPOVLP, NOOVLP, DONRCAP, NOCAPRE, NOCAPIM, CAPIMRE
      COMMON / DCCAPOVL / CAPOVLP, NOOVLP, DONRCAP,                     &
     &        NOCAPRE, NOCAPIM, CAPIMRE

      INTEGER  IPRCAP, IPCAPR, IPCAPI,                                  &
     &         IPCAP_DR, IPCAP_DI, IPCAP_S, LUPEMAT

      COMMON / DCCAPI / IPRCAP, IPCAPR, IPCAPI,                         &
     &         IPCAP_DR, IPCAP_DI, IPCAP_S, LUPEMAT
 
