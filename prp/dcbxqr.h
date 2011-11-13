C dcbxqr.h - control common block for pamxqr.h and pamtpa.h
C **************************************************************
C ***********    Quadratic response  <<A;B,C>>    **************
C **************************************************************
C
C MAXQR  - max number of allowed QR functions
C MAXLQR - max number of operators specified
C MAXFQR - max number of frequencies specified
C MZYEE  - max length of electronic response vectors
C MZYEP  - max length of positronic response vectors
C
C Used to read in USER specified requests
C----------------------------------------
C LAQROP - list of A-operator pointers
C LBQROP - list of B-operator pointers
C LCQROP - list of C-operator pointers
C IPRXQR - general print level
C ITRXQR - max number of iterations 
C NBQRFR - number of B frequencies 
C BQRFR  - B frequencies
C NCQRFR - number of C frequencies 
C CQRFR  - C frequencies
C XQR_SKIPEE - do not include e-e rotations
C XQR_SKIPEP - do not include e-p rotations
C XQR_ALLCMB - calculate all response functions disregarding 
C              permutational symmetries (debug option) 
C XQR_LEVICI - Include only functions with 3 different cartesian
C              indices (for .VERDET)   
C XQR_INDSTR - list of orbitals to include in rotation vectors
C
C Used to store the necessary linear response equations
C------------------------------------------------------
C NQROP  - number of unique operator pointers which are to
C          to be solved linear response for
C LQROP  - list of operator pointers, and number of 
C          frequencies for that operator
C QRFREQ - list of frequencies for operators
C
C Used to define nontrivial QR response functions
C------------------------------------------------
C NQRHYP  - number of nontrivial QR functions to be computed
C LQRHYP  - list of A/B/C operator pointers
C QRFRHYP - list of A/B/C frequencies
C
C Used in TPA calculations (EXCIT is also used by EXCPRP)
C--------------------------------------------------------
C IEXCOFF - offset when creating direct access labels for excitations
C TPACNV  - number of states per irrep
C TPASIM  - 
C TPASTV  - 
C EXCIT   - list of excitation frequencies in each irrep
C
C Used in EXCPRP calculations
C ---------------------------
C LSTATE - list of left excited states
C RSTATE - list of right excited states
C
      INTEGER MAXLQR,MAXFQR,MAXQR,IEXCOFF,LUQRINFO,LUQREE,LUQREP
C
      PARAMETER ( MAXLQR=10,MAXFQR=15,MAXQR=70,IEXCOFF=90 )
      PARAMETER ( LUQRINFO=40,LUQREE=41,LUQREP=42 )
C
      LOGICAL XQRNRM,XQRDIH,XQR_SKIPEE,XQR_SKIPEP,XQR_ALLCMB,
     &     XQR_LEVICI
      INTEGER IPRXQR,ITRXQR,NBQRFR,NCQRFR,NQROP,MAXQRM,
     &     LAQROP,NAQROP,LBQROP,NBQROP,
     &     LCQROP,NCQROP,LQROP,LQRHYP,NQRHYP,INTXQR,ITRIQR,
     &     MZYEE,MZYEP,TPACNV,TPASIM,TPASTV,LSTATE,RSTATE
      REAL*8 THCQR,BQRFR,CQRFR,QRFREQ,QRFRHYP,CNVXQR,RESXQR,EXCIT
      CHARACTER XQR_INDSTR*72
C
      COMMON /XCBLQR/ XQRNRM,XQRDIH,XQR_SKIPEE,XQR_SKIPEP,XQR_ALLCMB,
     &     XQR_LEVICI
      COMMON /XCBRQR/ THCQR,BQRFR(MAXFQR),CQRFR(MAXFQR),
     &     QRFREQ(3*MAXFQR,3*MAXLQR),QRFRHYP(MAXQR,3),
     &     CNVXQR(2),RESXQR,EXCIT(MAXFQR,8)
      COMMON /XCBIQR/ IPRXQR,ITRXQR,NBQRFR,NCQRFR,NQROP,MAXQRM,
     &     LAQROP(MAXLQR),NAQROP,LBQROP(MAXLQR),NBQROP,
     &     LCQROP(MAXLQR),NCQROP,LQROP(3*MAXLQR,2),
     &     LQRHYP(MAXQR,3),NQRHYP,INTXQR,ITRIQR(2),MZYEE,MZYEP,
     &     TPACNV(8),TPASIM(8),TPASTV(8),LSTATE(8),RSTATE(8)
      COMMON /XCBCQR/XQR_INDSTR(3,2)
C
C End of dcbxqr.h
C
