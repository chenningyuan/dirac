C
C   See the description in fck3aob.F
C
C   (the size of CCFTCH can be reduced by storing only 
C    the unique contraction coefficients )
C
C   EXPBTCH  - exponents
C   CCFBTCH  - contraction coefficients
C   INDXBTCH - shell indices for AO-batches
C   KNDXBTCH(1,IBATCH) - start adress for orbital of IBATCH
C   KNDXBTCH(2,IBATCH) - to be used for p with identical exponent
C   IAOCLS
C   NHKTBCH
C   NPRFBTCH
C   NCTFBTCH
C   KEXPBTCH
C   KCCFBCTH
C   NTYPBTCH
C   NCNTBTCH
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
