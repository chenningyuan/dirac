C*** Reduced response equations
C
C    JSYMOP - point group symmetry of operator
C    JTIMOP - time reversal symmetry of operator: 1(+),2(-)
C    JOPSY  - pointer to fermion ircop which spans JSYMOP
C
C    IPRXRS - print level
C    MAXITR  - max number of iterations 
C    NFREQ  - number of frequencies
C    NEXSIM - number of simultaneous solutions
C    NEXSTV - number of start vectors
C
C    THCXRS - threshold for convergence of reduced response equations
C     
C    Configuration and orbital rotation parameters
C    are held in separate vectors
C
C    Trial vectors have the general structure:
C
C     (b_1)^T = [ Z  Y^*]     (b_2)^T  =  [Y   Z^*]
C
C    
      CHARACTER XRSLBL*8,INTTYP*9
      COMMON/ XRSCCB/ XRSLBL(3),INTTYP
      LOGICAL RSREST, STATIC, LINEQ, LSVCFG, TKNORM, DIAGHE,STERNH
      COMMON /XRSICB/ IPRXRS,JSYMOP,JTIMOP,JOPSY,MAXITR,ITMIC,NFREQ,
     &                NEXSIM,NEXSTV,NEXCNV,LOFFTY,IPX,IMX,
     &                NZVAR, NZVARQ ,NZCONF,NZCONFQ,NTYP,
     &                NZXOPT,NZXOPTQ,NZXOPE,NZXOPEQ,NZXOPP,NZXOPPQ,
     &                NREDM ,N2REDM ,NZRED, NCRED ,NERED  ,NPRED, 
     &                KCONV ,NCSIM, NESIM, NPSIM, NOSIM, NTSIM,
     &                IZRED,ICRED,IERED,IPRED,LSVCFG(2),INTXRS,
     &                INTFLG,INTDEF,INTBUF,ITRINT(2)
      COMMON /XRSRCB/ THCXRS,CNVINT(2),RESFAC
      COMMON /WRKFLG/ STATIC, RSREST, LINEQ, TKNORM, DIAGHE,STERNH
      COMMON /MEMXRS/ KBTOT,KBCI,KBOE,KBOP,KZCONF,KZXOPE,KZXOPP
      COMMON /XRSIOU/ LUBCI,LUBOE,LUBOP,LUSCI,LUSOE,LUSOP,LURSP,LUXVC,
     &                LURST,LUCYCL
C      COMMON /INFXRS/ THCXRS, JSYMOP, JTIMOP, MAXIT, MXZRED,
C     &                NZVAR,  NZYVAR, NZCONF, NZYCONF,
C     &                NZXOPE, NZYXOPE,NZXOPP, NZYXOPP,NZXOPT, NZYXOPT,
C     &                LZVAR,  LZYVAR, LZCONF, LZYCONF,
C     &                LZXOPE, LZYXOPE,LZXOPP, LZYXOPP,LZXOPT, LZYXOPT,
C     &                JZC, JZE, JZP, JYC, JYE, JYP,
C     &                NEXCNV, NEXSIM, NEXSTV, LOFFTY, NZP,
C     &                RSREST, STATIC
C      COMMON /WRKXRS/ NZRED, NZYRED, JEXSIM, KCONV,
C     &                MCSIM, MESIM, MPSIM, MOSIM, MTSIM
C      CHARACTER*4 GPQUAT, E2QUAT, S2QUAT, BQUAT, SQUAT, EVRQUAT
C      COMMON /CINFXRS/ GPQUAT, E2QUAT, S2QUAT, BQUAT, SQUAT, EVRQUAT
