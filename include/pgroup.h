C
CTSAUE: Added GROUPS
C
C     SYMOP - Symmetry operations
C     ===========================
C
C         zyx
C     0   000    E   identity
C     1   001   Oyz  reflection in the yz-plane
C     2   010   Oxz  reflection in the xz-plane
C     3   011   C2z  rotation about the z-axis
C     4   100   Oxy  reflection in the xy-plane
C     5   101   C2y  rotation about the y-axis
C     6   110   C2x  rotation about the x-axis
C     7   111    i   inversion centre
C
      CHARACTER GROUP*3,SYMOP*3,REP*3,GROUPS*3
      COMMON/PGROUP/GROUP,SYMOP(0:7),REP(0:7),GROUPS(0:7)
