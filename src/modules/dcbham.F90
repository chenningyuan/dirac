!-----------------------------------------
! Module with variables previously defined
! at dcbham.h include file (COMMON /DCIHAM/, /DCLHAM/, /DCRHAM/)
! HFXFAC and HFXMU are also used from C (src/cfun), hence bind(C).
! INTV1C was EQUIVALENCEd to ICTLV1C(1); only ICTLV1C(2) is used besides it,
! so ICTLV1C keeps just that element (a stale ICTLV1C(1) is out of bounds).
!-----------------------------------------
!*****************************************************************************
!  DCBHAM:  Hamiltonian used for generation of wavefunction
!     IPRHAM     - Print level
!     IPOVERLAP  - Pointer to the overlap in DCBXPR
!     IPBETAMT   - Pointer to the beta matrix in DCBXPR
!     IPESHIFT   - Energy shift
!     IPMOLFLD   - Pointer to the nuclear attraction in DCBXPR
!     IPKINERG   - Pointer to the kinetic energy in DCBXPR
!     IPSOPP     - Pointer to the ECP spin-orbit matrix, SOPP DCBXPR
!     IPANGMOM   - Pointer to the z-component of angular momentum in DCBXPR
!     IPSPNMOM   - Pointer to the z-component of spin in DCBXPR
!     IPSPNORB   - Pointer to the dot between spin and angular moment in DCBXPR
!     MC         - number of components (1=only large, 2=large+small)
!     N1OPER     - Number of additional operators (finite field calculations)
!     IPR1OP     - Pointer to DCBXPR for additional operators
!     QED        - use QED Hamiltonian (filled positronic levels)
!     ONESYS     - No two-electron operators
!     FREEPJ     - Hamiltonian embedded in free particle projections
!     VEXTPJ     - Hamiltonian embedded in molecular field projections
!     NOSMLV     - No small nuclear attraction integrals
!     URKBAL     - Unrestricted kinetic balance
!     MDIRAC     - modified Dirac Hamiltonian
!     LEVYLE     - 4-component non-relativistic Levy-Leblond Hamiltonian
!     NONRAL     - 2-component non-relativistic Pauli Hamiltonian
!     ECPCALC    - RECP calculation
!     SUB_BL     - Diagonalize Hamiltonian blockwise (neglect couplings)
!     SSMTRC     - Scaling factor for SS overlap matrix
!     SSMTRC_SAVE- Scaling factor for SS overlap matrix, temporary variable
!     SPINFR     - Neglect spin-orbit coupling
!     SPINFR2    - Neglect 2-electron spin-orbit coupling
!     X2CMMF     - Use molecular mean field approximation (^2 DC(G) (U_mol)) in correlated calcs.
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
!     BSS        - flag for generating matrix elements of
!                  Barysz-Sadlej-Snijders transformed two-component Hamiltonians
!     IBSS       - what kind of BSS Hamiltonian(1.,2. or infinite order; scalar of full
!                  relativity)
!     X2C        - flag for generating matrix elements of
!                  the exact two-component Hamiltonian using the X2C module (i.e. X2Cmod_x2c == .true. (module x2cmod_cfg))
!     X2C_aba    - flag for generating matrix elements of
!                  the exact two-component Hamiltonian using the X2C module (i.e. X2Cmod_x2c == .true. (module x2cmod_cfg))
!                  in the molecular fragment/atom-by-atom approach
!     IMFCH      - charge of the system to be spliited among mean-field summations of AMFI centers
!     NOPRTR     - do not perform the preliminary transformation to the ON basis in the 2c mode
!                  (only for infinite order Tronds approach), if .false., we do the preliminary
!                  transformation to the free-particle basis or to the Vext-particle ON basis
!     BLOCKD     - do the numerical block diagonalization (only for the infinite order)
!     DO2C4C     - run SCF first in 2c mode, then continue in the 4-component mode
!     START2C    - if .true., run preliminary 2com BSS-SCF before 4comp DC-SCF (used together with DO2C4C)
!     INI2C      - integer number specifying what two-component quantity will be used for (re)starting
!                  the DC-SCF from previous BSS-SCF
!     DO4C2C     - after finishing 4c DC-SCF do block diagonalize 4c Fock and continue
!                  at the 2c level
!     USEDF      - if true, use the (converged) Dirac-Fock operator to be block-diagonalized
!                  (instead of 1 electron Dirac operator)
!     CONT2C     - continue with 2-component (inf.order) SCF after DC-SCF
!     I2CHAM     - integer number (1,2...) specifynig the type of post-DC-SCF 2component Hamiltonian
!     RESTFCK    - flag to denote restarting DC-SCF from BSS-SCF, or vice versa, from 2c/4c Fock matrix in MO basis
!     I2COFK     - integer number specifying which 2comp. Hamiltonian (in AO basis) is to be read in the ONEFCK routine
!     TRMO4C2C   - flag meaning transformation of 4c MO into 2c form using the picture change transformation matrix
!     ALLAMFA    - flag about if all atoms are present for the AMFI contributions.
!     ISORDER_AMFI_x2c - integer value determining the order of AMFI contributions to the X2C Hamiltonian.
!     CMPEIG     - for IOTC/BSS: compare eigenvalues of 2c and 4c Hamiltonian
!     DKH1       - own flag for the 1st order DKH,  very rough approximation
!     DKH2       - own flag for the second order DKH, classic
!     GAUNT      - Include Gaunt interaction in the 2-e operator
!     CAP        - flag for complex absoption potential used only with the FSCC method
!     NOAMFI      - do not include AMFI contributions at all
!     QDOTS      - use quantum dots harmonic potential instead of the electron-nuclei attraction potential
!     FAKE2C     - do not compute the 2C 1e-Hamiltonian and take the 2C spinors
!                  from previous run or import from different QC package
!     LSCALE_DFT_GAUNT - If true the Gaunt contribution will be scaled together with exchange in DFT
!                  calculations. If false 100% of the Gaunt contribution (if turned on) will be
!                  included in DFT calculations.
!     NOSFMU     - in spinfree calculations, do not generate spinfree multiplicaton tables, that is
!                  constructing direct products of spin and spatial symmetries
!     TRS_MIX - if .true.,allow mixing T+ Hamiltonian with chosen T- perturbative operator
   module dcbham_mod
   use, intrinsic :: iso_c_binding, only: c_double
   implicit none
   private :: c_double
   integer, parameter :: MAX1OPER = 10
   real*8 :: SSMTRC, SSMTRC_SAVE, EPSOL, E_PE, EPSTAT, RSOLAV, HFXATT, ONECNV, ESHIFT
   real(c_double), bind(C, name='dcrham_hfxfac') :: HFXFAC
   real(c_double), bind(C, name='dcrham_hfxmu') :: HFXMU
   integer :: IPRHAM, IPOVRLAP, IPBETAMT, IPMOLFLD, IPKINERG, IPANGMOM, IPSPNMOM, IPSPNORB, &
              N1OPER, IPR1OP(MAX1OPER), LSOLMX, IPRSOL, INTV1C, ICTLV1C(2:2), IBSS, MC, IMFCH, &
              I2CHAM, INI2C, I2COFK, IPVEMB0, IPSOPP, ISORDER_AMFI_X2C, IPESHIFT
   logical :: QED, ONESYS, FREEPJ, VEXTPJ, NOSMLV, URKBAL, LEVYLE, SUB_BL, SOLVEN, INERSI, &
              INERSF, SPINFR, SPINFR2, ZORA, ZORA4, ZORASC, PROJEC, PROOWN, SMLV1C, ONECAP, &
              ONECOFF, NOSPIN, BSS, NOPRTR, BLOCKD, YREQ1, DO2C4C, DO4C2C, USEDF, RESTFCK, &
              CONT2C, START2C, GAUNT, ALLAMFA, RDSUSLLMOD, CAP, CMPEIG, DKH1, DKH2, NOAMFI, &
              QDOTS, LSCALE_DFT_GAUNT, NOSFMU, TRMO4C2C, X2CMMF, X2C, X2C_ABA, ECPCALC, &
              MDIRAC, TRS_MIX, WRITE_JZ_MATRIX, NONREL, WRITE_K_MATRIX, FAKE2C, DOEFFQED, &
              DOASHIFT
   end module dcbham_mod
