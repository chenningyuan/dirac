!*****************************************************************************
!  DCBESR:  ESR properties
!
!  IPSIGMA(i): pointer to the spin operator Sigma_i
!
      INTEGER MAXESR
      PARAMETER (MAXESR = 50)

      INTEGER IPRESR, NESRP, LESRP(MAXESR), NCIESR, IPSIGMA(3)
      COMMON /DCBESR/ IPRESR, NESRP, LESRP, NCIESR, IPSIGMA
