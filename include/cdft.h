C     common block for spin/current density response
C     ==============================================
C
C     LBIT(LBIT_CDFT,1) spin density contribution
C     LBIT(LBIT_CDFT,2) current density contribution
C     LBIT(LBIT_CDFT,3) collinear approximation
C     LBIT(LBIT_CDFT,4) spin density operator is \beta*\Sigma instead of \Sigma
C     LBIT(LBIT_CDFT,5) X\alpha type current density contribution
C
      COMMON /CBLCDI/ LBIT_CDFT
      COMMON /CBLCDR/ CDFT_ALPHA
