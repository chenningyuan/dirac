!
!TSAUE: Added GROUPS
!
!     SYMOP - Symmetry operations
!     ===========================
!
!         zyx
!     0   000    E   identity
!     1   001   Oyz  reflection in the yz-plane
!     2   010   Oxz  reflection in the xz-plane
!     3   011   C2z  rotation about the z-axis
!     4   100   Oxy  reflection in the xy-plane
!     5   101   C2y  rotation about the y-axis
!     6   110   C2x  rotation about the x-axis
!     7   111    i   inversion centre
!
      CHARACTER*3 GROUP, SYMOP(0:7), REP(0:7), GROUPS(0:7)
      COMMON /PGROUP/ GROUP, SYMOP, REP, GROUPS
