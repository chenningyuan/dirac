! *** Wave function module ***
!  DOSCF - either DODHF or DODFT true
!  DODHF - perform Dirac-Hartee-Fock calculation
!  DOMP2 - perform direct second-order Moeller-Plesset calculation
!  DODFT - perform Kohn-Sham DFT calculation
!  DORES - resolution of open-shell manifold with GOSCIP
!  DOMVO - make modified virtual orbitals
!  DOCCM - perform coupled-cluster calculation with RELCCSD
!  DOCIM - perform CI calculation with DIRRCI/GOSCIP
!  DOLUCT- perform CI calculation with LUCITA
      LOGICAL DOSCF, DODHF, DODFT, DOMP2, DOMVO, DORES, DOCCM, DOCIM,   &
     &        DOLUCT, DOKRMC, DOLUCIAR
      COMMON /DCLPSI/ DOSCF, DODHF, DODFT, DOMP2, DOMVO, DORES, DOCCM,  &
     &                DOCIM, DOLUCT, DOKRMC, DOLUCIAR
