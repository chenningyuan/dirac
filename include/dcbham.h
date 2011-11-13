#if !defined(__CVERSION__)

! VERSION : $Revision: 1.18 $
! DATE    : $Date: 2004/05/26 13:06:54 $
! FILE    : dcbham.h
!
!*****************************************************************************
!  DCBHAM:  Hamiltonian used for generation of wavefunction
!
!     IPRHAM     - Print level
!     IPOVERLAP  - Pointer to overlap in DCBXPR
!     IPBETAMT   - Pointer to beta matrix in DCBXPR
!     IPMOLFLD   - Pointer to nuclear attraction in DCBXPR
!     IPKINERG   - Pointer to kinetic energy in DCBXPR
!     IPANGMOM   - Pointer to z-component of angular momentum in DCBXPR
!     IPSPNMOM   - Pointer to z-component of spin in DCBXPR
!     N1OPER     - Number of additional operators (finite field calculations)
!     IPR1OP     - Pointer to DCBXPR for additional operators
!     QED        - use QED Hamiltonian (filled positronic levels)
!     ONESYS     - No two-electron operators
!     FREEPJ     - Hamiltonian embedded in free particle projections
!     VEXTPJ     - Hamiltonian embedded in molecular field projections
!     NOSMLV     - No small nuclear attraction integrals
!     URKBAL     - Unrestricted kinetic balance
!     LEVYLE     - Levy-Leblond Hamiltonian
!     SUB_BL     - Diagonalize Hamiltonian blockwise (neglect couplings)
!     SSMTRC     - Scaling factor for SS overlap matrix
!     SPINFR     - Neglect spin-orbit coupling
!     SPINFR2    - Neglect 2-electron spin-orbit coupling
!     ZORA       - Make ZORA approximation (ZORA: approx density, cheap)
!     ZORA4      - Make ZORA approximation (ZORA4: full density, expensive)
!     SMLV1C     - Neglect two-center blocks of VSS
!     ONECAP     - Do one-center approximation depending on model INTV1C
!                  Default is type 2 : Neglect two-center blocks of VSS
!                  ONECAP type 1 model. One-center (SS|SS) and (LL|SS).
!                  Correction made to one-electron Fock-matrix.
!                  ONECAP type 3 model. One-center (SS|SS) and (LL|SS).
!                  Correction made to two-electron Fock-matrix.
!     ONECOFF    - internal flag for turning off ONECAP
!
      INTEGER MAX1OPER
      PARAMETER (MAX1OPER = 10)

      INTEGER IPRHAM, IPOVRLAP, IPBETAMT, IPMOLFLD, IPKINERG, IPANGMOM, &
     &        IPSPNMOM, N1OPER, IPR1OP(MAX1OPER), LSOLMX, IPRSOL,       &
     &        ICTLV1C(2)
      COMMON /DCIHAM/ IPRHAM, IPOVRLAP, IPBETAMT, IPMOLFLD, IPKINERG,   &
     &                IPANGMOM, IPSPNMOM, N1OPER, IPR1OP, LSOLMX,       &
     &                IPRSOL, ICTLV1C

      INTEGER INTV1C
      EQUIVALENCE (INTV1C, ICTLV1C (1))

      LOGICAL QED, ONESYS, FREEPJ, VEXTPJ, NOSMLV, URKBAL, LEVYLE,      &
     &        SUB_BL, SOLVEN, INERSI, INERSF, SPINFR, SPINFR2, ZORA,    &
     &        ZORA4, ZORASC, PROJEC, PROOWN, DFTHAM, SMLV1C, ONECAP,    &
     &        ONECOFF, NOSPIN
      INTEGER NDCLHAM
      PARAMETER (NDCLHAM = 18)

      COMMON /DCLHAM/ QED, ONESYS, FREEPJ, VEXTPJ, NOSMLV, URKBAL,      &
     &                LEVYLE, SUB_BL, SOLVEN, INERSI, INERSF, SPINFR,   &
     &                SPINFR2, ZORA, ZORA4, ZORASC, PROJEC, PROOWN,     &
     &                DFTHAM, SMLV1C, ONECAP, ONECOFF, NOSPIN

      REAL*8 SSMTRC, EPSOL, EPSTAT, RSOLAV, HFXFAC, ONECNV
      COMMON /DCRHAM/ SSMTRC, EPSOL, EPSTAT, RSOLAV, HFXFAC, ONECNV

#else  /* __CVERSION__ */
#define MAX1OPER 10
  extern struct common_dcrham {
    real ssmtrc,epsol,epstat,rsolav,hfxfac;
  } dcrham_;

#endif /* __CVERSION__ */
