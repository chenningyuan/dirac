C     common block for asymptotic corrections
C     =======================================
C
C   1 LBIT_DFTAC =        1 * GRAC
C   2                 +   2 * GRAC is active
C   3                 +   4 * SAOP
C   4                 +   8 * SAOP is active
C   5                 +  16 * SAOP with response part
C   6                 +  32 * asymptotic part is LB94
C   7                 +  64 * asymptotic part is LBalpha
C   8                 + 128 * debug
C
      COMMON /CBLACR/ DFTAC_IP, DFTAC_THRESHOLD, HOMO_EIGENVALUE,       &
     &                GRAC_ALPHA, GRAC_BETA
C
      COMMON /CBLACI/ LBIT_DFTAC, ICOUNTER_DFTAC_DEBUG
