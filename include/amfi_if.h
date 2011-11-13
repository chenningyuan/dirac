!----------------------------------------------------------------------------
!
!  Common file for passing information from DIRAC to AMFI/RELSCF
! ---------------------------------------------------------------
!
!   IPR_AMFI  - print level variable for AMFI part of B.Schimm.
!
!  IPR_RELSCF - print level for relscf program interfaced to AMFI
!
!   MXXTRP_INP - maximum number of iterations for relscf
!
! ... other parameters can be added here as well ...
!
!----------------------------------------------------------------------------
      INTEGER IPR_AMFI,IPR_RELSCF,MXXTRP_INP

      COMMON /AMFI_RELSCF/ IPR_AMFI,IPR_RELSCF,MXXTRP_INP

