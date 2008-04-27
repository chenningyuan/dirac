!----------------------------------------------------------------------------
!
!
!   IPR_AMFI  - print level variable for AMFI part of B.Schimm.
!
!  IPR_RELSCF - print level for relscf program interfaced to AMFI
!
!
!----------------------------------------------------------------------------
      LOGICAL AMFI_VERBOSE,AMFI1_VERBOSE,AMFI2_VERBOSE,
     &                     AMFI3_VERBOSE,RELSCF_VERBOSE,DIRBSS_VERBOSE

      COMMON /AMFI_CONTROL/AMFI_VERBOSE,AMFI1_VERBOSE,AMFI2_VERBOSE,
     &                     AMFI3_VERBOSE,RELSCF_VERBOSE,DIRBSS_VERBOSE

      INTEGER IPR_AMFI,IPR_RELSCF

      COMMON /AMFI_RELSCF/ IPR_AMFI,IPR_RELSCF

