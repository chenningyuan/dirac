!
! *** Wave function module ***
!
!  DOSCF - either DODHF or DODFT true
!  DODHF - perform Dirac-Hartee-Fock calculation
!  DOMP2 - perform direct second-order Moeller-Plesset calculation
!  DODFT - perform Kohn-Sham DFT calculation
!  L1ORBM- orbital pre-modifications (reordering, rotations) 
!  L2ORBM- orbital post-modifications (reordering, phase adjustment) 
!  DORES - resolution of open-shell manifold with GOSCIP
!  DOMVO - make modified virtual orbitals
!  DOCCM - perform coupled-cluster calculation with RELCCSD
!  DOHSCC - perform HSMRCC calculation with external code
!  DOCOSCI - perform CI calculation with GOSCIP
!  DOCIM - perform CI calculation with DIRRCI/GOSCIP
!  DOLUCT- perform CI calculation with LUCITA
!  DOADC - compute ionization spectrum via Green's functions

      LOGICAL DOSCF, DODHF, DODFT, DOMP2, DOMVO, DORES, DOCCM, DOCIM,   &
     &        DOHSCC, DOLUCT, DOKRMC, DOLUCIAR,                         &
     &        L1ORBM, L2ORBM, PHCOEF, DOCOSCI, DOADC

      COMMON /DCLPSI/ DOSCF, DODHF, DODFT, DOMP2, DOMVO, DORES, DOCCM,  &
     &                DOCIM, DOHSCC, DOLUCT, DOKRMC, DOLUCIAR, L1ORBM,  &
     &                L2ORBM, PHCOEF, DOCOSCI, DOADC
