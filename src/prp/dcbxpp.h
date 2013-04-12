C
C *** PP (Polarization propagator) excitation properties  ***
C
C MAXLPP - max number of A operators specified
C
C NPPAP(NBSYM) - number of A-operators classified 
C      according to boson symmetries
C LPPAPU - pointer to list of A - operators (unsorted)
C LPPAPS - pointer to list of A - operators sorted on symmetries
C JPPAP  - symmetry offsets in LPPAPS
C
C IPRXPP - general print level
C ITRXPP - max number of iterations 
C THCPP  - convergence threshold
C KEXCNV(8) - number of excitations to converge in each boson sym.
C KEXSIM(8) - number of simult. trial vectors in each boson sym.
C KEXSTV(8) - number of start   trial vectors in each boson sym.
C MAXEXC = max(KEXCNV(i), i = 1,nbsym)
C
C    dcbxpp.h -- Feb. 18, 1999, hjaaj
C
      PARAMETER ( MAXLPP = 40 )
      LOGICAL        XPP_SKIPEE,XPP_SKIPEP,XPP_LSFG   ,
     &        XPPNRM,ORBXPP,XPPANA,XPP_TRIPLET,XPP_E2CHEK
      COMMON /XCBLPP/XPP_SKIPEE,XPP_SKIPEP,XPP_LSFG(2),
     &        XPPNRM,ORBXPP,XPPANA,XPP_TRIPLET,XPP_E2CHEK
      COMMON /XCBRPP/THCPP,CNVXPP(2),RESXPP,XPPERG
      COMMON /XCBIPP/IPRXPP,ITRXPP,MAXEXC,
     &     KEXCNV(8), KEXSIM(8), KEXSTV(8),
     &     LPPAPU(MAXLPP),LPPAPS(MAXLPP),NPPAPT,NPPAP(8),JPPAP(8),
     &     INTXPP,ITRIPP(2),NCLS_XPP(2,3,2)
      CHARACTER XPP_INDSTR*72, XPPFIL*6
      COMMON/DCCXPP/XPP_INDSTR(3,2),XPPFIL
