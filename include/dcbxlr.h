C
C *** Linear response  <<A,B>>  ***
C
C MAXLLR - max number of operators specified
C MAXFLR - max number of frequencies specified
C
C NLRAP/NLRBP(NBSYM,ITIM) - number of A/B-operators classified 
C      according to bosonsymmetries and time reversal symmetry
C LLRAPU/LLRBPU - pointer to list of A/B - operators (unsorted)
C LLRAPS/LLRBPS - pointer to list of A/B - operators sorted on symmetries
C JLRAP/JLRBP   - symmetry offsets in LLRAPS/LLRBPS
C
C IPRXLR - general print level
C ITRXLR - max number of iterations 
C MAXRM  - max dimension of reduced matrix
C THCLR  - convergence threshold
C NBFREQ - number of frequencies 
C BFREQ  - frequencies
C INDXLR is pointer to current B-operator in LLRAPU
C LSFG(1)=T   - call GMOLI in sigma-vector generation
C LSFG(2)=T   - call FMOLI in sigma-vector generation
C
*COMDECK XCBLR
      PARAMETER ( MAXLLR = 40, MAXFLR = 20 )
      LOGICAL ALLCMB,SKIPEP,QEDRUN,LSFG,XLRNRM,XLRDIH,TRIAB,
     &        COMPRS,ORBXLR,XSTERN,SKIPEE,XLRANA
      COMMON /XCBLLR/ALLCMB,SKIPEP,QEDRUN,LSFG(2),XLRNRM,XLRDIH,
     &               TRIAB,COMPRS,ORBXLR,XSTERN,SKIPEE,XLRANA
      COMMON /XCBRLR/THCLR,BFREQ(MAXFLR),THRCOM,CNVXLR(2),RESXLR
      COMMON /XCBILR/ IPRXLR,ITRXLR,NBFREQ,INDBP,MAXRM,MXLOAD,
     &     LLRAPU(MAXLLR),LLAPSU(MAXLLR),NLRAPT,NLRAP(8),JLRAP(8),
     &     LLRBPU(MAXLLR),LLBPSU(MAXLLR),NLRBPT,NLRBP(8),JLRBP(8),
     &     INTXLR,ITRILR(2),NCLS(2,3,2)
      CHARACTER XLR_INDSTR*72, XVCFIL*6
      COMMON/DCCXLR/XLR_INDSTR(3,2),XVCFIL
