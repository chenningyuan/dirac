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
C QREE   - include e-e rotations
C QREP   - include e-p rotations
C
C Used to store the necessary linear response equations
C------------------------------------------------------
C NQROP - number of unique operator pointers which are to
C         to be solved linear response for
C LQROP - list of operator pointers, and number of 
C         frequencies for that operator
C QRFREQ - list of frequencies for operators
C
C Used to define nontrivial QR response functions
C------------------------------------------------
C NQRHYP - number of nontrivial QR functions to be computed
C LQRHYP - list of A/B/C operator pointers
C QRFRHYP - list of A/B/C frequencies
C
*COMDECK XCBQR
      PARAMETER ( MAXLQR=10,MAXFQR=10,MAXQR=60 )
      PARAMETER ( LUQRINFO = 40, LUQREE = 41, LUQREP = 42 )
      LOGICAL XQRNRM,XQRDIH,QREE,QREP
      COMMON /XCBLQR/ XQRNRM,XQRDIH,QREE,QREP
      COMMON /XCBRQR/ THCQR,BQRFR(MAXFQR),CQRFR(MAXFQR),
     &     QRFREQ(3*MAXFQR,3*MAXLQR),QRFRHYP(MAXQR,3),
     &     CNVXQR(2),RESXQR
      COMMON /XCBIQR/ IPRXQR,ITRXQR,NBQRFR,NCQRFR,NQROP,MAXQRM,
     &     LAQROP(MAXLQR),NAQROP,LBQROP(MAXLQR),NBQROP,
     &     LCQROP(MAXLQR),NCQROP,LQROP(3*MAXLQR,2),
     &     LQRHYP(MAXQR,3),NQRHYP,INTXQR,ITRIQR(2),MZYEE,MZYEP
