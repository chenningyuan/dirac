C *** Wave function module ***
C  DOSCF - either DODHF or DODFT true
C  DODHF - perform Dirac-Hartee-Fock calculation
C  DOMP2 - perform direct second-order Moeller-Plesset calculation
C  DODFT - perform Kohn-Sham DFT calculation
C  DORES - resolution of open-shell manifold with GOSCIP
C  DOMVO - make modified virtual orbitals
C  DOCCM - perform coupled-cluster calculation with RELCCSD
C  DOCIM - perform CI calculation with DIRRCI/GOSCIP
C  DOLUCT- perform CI calculation with LUCITA
      LOGICAL       DOSCF,DODHF,DODFT, DOMP2, DOMVO,DORES,
     &              DOCCM,DOCIM,DOLUCT,DOKRMC,DOLUCIAR
      COMMON/DCLPSI/DOSCF,DODHF,DODFT, DOMP2, DOMVO,DORES,
     &              DOCCM,DOCIM,DOLUCT,DOKRMC,DOLUCIAR
