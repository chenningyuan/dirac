C*****************************************************************************
C  DCBEXP:  Expectation values
C
      PARAMETER(MAXEXP=350)
      LOGICAL ORBANA
      COMMON/DCBEXP/IPREXP,NEXPP,LEXPP(MAXEXP),ORBANA,
     &              IPDIPOLE(3),IPQUADRU(6),IPEFG(6*MXCENT)

