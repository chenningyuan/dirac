C
C   See the description in fck3aob.F
C
C   (the size of CCFTCH can be reduced by storing only 
C    the unique contraction coefficients )
C
C      PARAMETER (MXAOBTCH=MXSHEL)
      PARAMETER (MXAOBTCH=MXSHEL*MXCONT)
C
      COMMON /FCK3AOBTCH/ EXPAOMIN, CCFAOMAX,
     &         EXPBTCH(MXPRIM),    CCFBTCH(MXPRIM*MXCONT),   
     &         INDXBTCH(MXSHEL*MXCONT), IAOCLS(3,MXSHEL),
     &         KNDXBTCH(2,MXAOBTCH),NHKTBTCH(MXAOBTCH),  
     &         NPRFBTCH(MXAOBTCH),  NCTFBTCH(MXAOBTCH),
     &         KEXPBTCH(MXAOBTCH),  KCCFBTCH(2,MXAOBTCH),
     &         NTYPBTCH(MXAOBTCH),  NCNTBTCH(MXAOBTCH),
     &         MAXAONHK, MAXAOKCK, MAXAOKHK, MAXAOPRF,
     &         MAXAOCTF, NAOBTCH, NAOBCLS, NBASE
