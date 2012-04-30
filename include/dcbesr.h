!*****************************************************************************
!  dcbesr.h:  ESR properties
!
!  IPSIGMA(i) : pointer to the spin operator Sigma_i
!  LVL_CIESR>0: include spin-polarization with
!               GAS1 with 0:LVL_CIESR holes in NCESR(2) core orbitals
!               GAS2 with open shell from DHF (no occ. restriction)
!               GAS3 with 0:LVL_CIESR elec.s in NVESR(2) virtual orbitals
!
      INTEGER MAXESR, NPSIGMA
      PARAMETER (MAXESR = 90, NPSIGMA = 3)
      REAL*8  THR_CIESR, THR_CVEC_ESR
      INTEGER IPRESR,NESRP,LESRP(MAXESR),N_CIESR,IPSIGMA(NPSIGMA),      &
     &        MULTIPESR, MAX_CIESR_IT,                                  &
     &        NCESR(2),NVESR(2),LVL_CIESR
      COMMON /DCBESR/ THR_CIESR,THR_CVEC_ESR,                           &
     &        IPRESR,NESRP,LESRP,N_CIESR,IPSIGMA,                       &
     &        MULTIPESR, MAX_CIESR_IT,NCESR,NVESR,LVL_CIESR
! --- end of dcbesr.h ---
