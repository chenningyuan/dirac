#if !defined(__CVERSION__)
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbham.h
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
C     QED        - use QED Hamiltonian (filled positronic levels)
C     ONESYS     - No two-electron operators
C     FREEPJ     - Hamiltonian embedded in free particle projections
C     VEXTPJ     - Hamiltonian embedded in molecular field projections
C     NOSMLV     - No small nuclear attraction integrals
C     URKBAL     - Unrestricted kinetic balance
C     LEVYLE     - Levy-Leblond Hamiltonian
C     SSMTRC     - Scaling factor for SS overlap matrix
C     SPINFR     - Neglect spin-orbit coupling
C     SPINFR2    - Neglect 2-electron spin-orbit coupling
C     ZORA       - Make ZORA approximation (ZORA: approx density, cheap)
C     ZORA4      - Make ZORA approximation (ZORA4: full density, expensive)
C     SMLV1C     - Neglect two-center blocks of VSS
C     ONECAP     - Do one-center approximation depending on model INTV1C
C                  Default is type 2 : Neglect two-center blocks of VSS
C                  ONECAP type 1 model. One-center (SS|SS) and (LL|SS).
C                  Correction made to one-electron Fock-matrix.
C                  ONECAP type 3 model. One-center (SS|SS) and (LL|SS).
C                  Correction made to two-electron Fock-matrix.
C     ONECOFF    - internal flag for turning off ONECAP
C
      PARAMETER(MAX1OPER=10)
      COMMON/DCIHAM/IPRHAM,IPOVRLAP,IPBETAMT,IPMOLFLD,IPKINERG,
     &              N1OPER,IPR1OP(MAX1OPER),LSOLMX,IPRSOL,
     &              ICTLV1C(2)
      EQUIVALENCE (INTV1C, ICTLV1C(1))
      LOGICAL QED,ONESYS,FREEPJ,VEXTPJ,NOSMLV,URKBAL,LEVYLE,SOLVEN,
     &        INERSI,INERSF,SPINFR,SPINFR2,ZORA,ZORA4,ZORASC,
     &        PROJEC,PROOWN,DFTHAM,SMLV1C,ONECAP,ONECOFF
      PARAMETER (NDCLHAM=17)
      COMMON/DCLHAM/
     &        QED,ONESYS,FREEPJ,VEXTPJ,NOSMLV,URKBAL,LEVYLE,SOLVEN,
     &        INERSI,INERSF,SPINFR,SPINFR2,ZORA,ZORA4,ZORASC,
     &        PROJEC,PROOWN,DFTHAM,SMLV1C,ONECAP,ONECOFF
      COMMON/DCRHAM/SSMTRC,EPSOL,EPSTAT,RSOLAV,HFXFAC,ONECNV
#else  /* __CVERSION__ */
#define MAX1OPER 10
  extern struct common_dcrham {
    real ssmtrc,epsol,epstat,rsolav,hfxfac;
  } dcrham_;

#endif /* __CVERSION__ */
