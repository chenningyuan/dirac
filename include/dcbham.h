#if !defined(__CVERSION__)
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbham.h
C SOURCE  : $Source$
C
C*****************************************************************************
C  DCBHAM:  Hamiltonian used for generation of wavefunction
C
C     IPRHAM     - Print level
C     IPOVERLAP  - Pointer to overlap in DCBXPR
C     IPBETAMT   - Pointer to beta matrix in DCBXPR
C     IPMOLFLD   - Pointer to nuclear attraction in DCBXPR
C     IPKINERG   - Pointer to kinetic energy in DCBXPR
C     N1OPER     - Number of additional operators (finite field calculations)
C     IPR1OP     - Pointer to DCBXPR for additional operators
C     ONESYS     - No two-electron operators
C     FREEPJ     - Hamiltonian embedded in free particle projections
C     VEXTPJ     - Hamiltonian embedded in molecular field projections
C     PROJEC     - any projection, defined by a set of coefficients
C     NOSMLV     - No small nuclear attraction integrals
C     URKBAL     - Unrestricted kinetic balance
C     LEVYLE     - Levy-Leblond Hamiltonian
C     SSMTRC     - Scaling factor for SS overlap matrix
C     SPINFR     - Neglect spin-orbit coupling
C     ZORA       - Make ZORA approximation (ZORA: approx density, cheap)
C     ZORA4      - Make ZORA approximation (ZORA4: full density, expensive)
C     SMLV1C     - Neglect two-center blocks of VSS
C
      PARAMETER(MAX1OPER=10,MPRJREF=6)
      CHARACTER PRJFIL*6,VCPROJ*72
      COMMON/DCCHAM/PRJFIL(MPRJREF),VCPROJ(2,MPRJREF)
      COMMON/DCIHAM/IPRHAM,IPOVRLAP,IPBETAMT,IPMOLFLD,IPKINERG,
     &              N1OPER,IPR1OP(MAX1OPER),LSOLMX,IPRSOL,
     &              NPRJREF
      LOGICAL ONESYS,FREEPJ,VEXTPJ,NOSMLV,URKBAL,LEVYLE,SOLVEN,
     &        INERSI,INERSF,SPINFR,ZORA,ZORA4,ZORASC,SMLV1C,
     &        PROJEC,PROOWN,DFTHAM
      PARAMETER (NDCLHAM=17)
      COMMON/DCLHAM/ONESYS,FREEPJ,VEXTPJ,NOSMLV,URKBAL,LEVYLE,SOLVEN,
     &        INERSI,INERSF,SPINFR,ZORA,ZORA4,ZORASC,SMLV1C,PROJEC,
     &        PROOWN,DFTHAM
      COMMON/DCRHAM/SSMTRC,EPSOL,EPSTAT,RSOLAV,PRJTHR,HFXFAC
#else  /* __CVERSION__ */
#define MAX1OPER 10
#define MPRJREF 6
  extern struct common_dcrham {
    real ssmtrc,epsol,epstat,rsolav,prjthr,hfxfac;
  } dcrham_;

#endif /* __CVERSION__ */
