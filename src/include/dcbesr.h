!*****************************************************************************
!  dcbesr.h:  ESR properties
!
!  IPSIGMA(i) : pointer to the spin operator Sigma_i
!  LVL_CIESR(i)>0: include spin-polarization with
!               GAS0 with 0:LVL_CIESR(0) holes in NCSH_ESR(1:2) core orbitals
!               GAS1 with 0:LVL_CIESR(1) holes in NVSH_ESR(1:2) valence orbitals
!               GAS2 with open shells from KR-HF or KR-MCSCF (no occ. restriction)
!               GAS3 with 0:LVL_CIESR(2) electrons in NSSH_ESR(1:2) virtual orbitals
!
      INTEGER MAXESR, NPSIGMA
      PARAMETER (MAXESR = 90, NPSIGMA = 3)
      REAL*8  THR_CIESR, THR_CVEC_ESR, THRPCI_ESR
      INTEGER IPRESR,NESRP,LESRP(MAXESR),N_CIESR,IPSIGMA(NPSIGMA),      &
     &        MULTIPESR, MAX_CIESR_IT,                                  &
     &        NCSH_ESR(2),NVSH_ESR(2),NSSH_ESR(2),LVL_CIESR(0:2)
      COMMON /DCBESR/ THR_CIESR,THR_CVEC_ESR,THRPCI_ESR,                &
     &        IPRESR,NESRP,LESRP,N_CIESR,IPSIGMA, MULTIPESR,            &
     &        MAX_CIESR_IT,NCSH_ESR,NVSH_ESR,NSSH_ESR,LVL_CIESR
! --- end of dcbesr.h ---
