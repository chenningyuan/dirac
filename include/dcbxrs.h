C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbxrs.h
C SOURCE  : $Source$
C
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
      CHARACTER XRSLBL*8,INTTYP*9,BVTYP*2
      COMMON/ XRSCCB/ XRSLBL(3),INTTYP,BVTYP(6)
      LOGICAL RSREST, STATIC, LINEQ, LSVCFG, TKNORM, DIAGHE,STERNH
      COMMON /XRSICB/ IPRXRS,JSYMOP,JTIMOP,JOPSY,MAXITR,ITMIC,NFREQ,
     &                NEXSIM,NEXSTV,NEXCNV,LOFFTY,IPX,IMX,
     &                NZVAR, NZVARQ ,NZCONF,NZCONFQ,NTYP,
     &                NZXOPT,NZXOPTQ,NZXOPE,NZXOPEQ,NZXOPP,NZXOPPQ,
     &                NVPAR ,IPEP,IPEM,IPPP,IPPM,IPCP,IPCM,
     &                NREDM ,N2REDM ,NZRED, NCRED ,NERED  ,NPRED, 
     &                KCONV ,NCSIM, NESIM, NPSIM, NOSIM, NTSIM,
     &                IZRED,ICRED,IERED,IPRED,LSVCFG(2),INTXRS,
     &                INTFLG,INTDEF,INTBUF,ITRINT(2),NFOCKMAT
      COMMON /XRSRCB/ THCXRS,CNVINT(2),RESFAC
      COMMON /WRKFLG/ STATIC, RSREST, LINEQ, TKNORM, DIAGHE,STERNH
      COMMON /MEMXRS/ KBTOT,KBCI,KBOE,KBOP,KZCONF,KZXOPE,KZXOPP
      COMMON /XRSIOU/ LUBCI,LUBOE,LUBOP,LUSCI,LUSOE,LUSOP,LURSP,LUXVC,
     &                LURST,LUCYCL
