C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbexp.h
C
C*****************************************************************************
C  DCBEXP:  Expectation values
C
      PARAMETER(MAXEXP=350)
      LOGICAL ORBANA
      COMMON/DCBEXP/IPREXP,NEXPP,LEXPP(MAXEXP),
     &              IPDIPOLE(3),IPQUADRU(6),IPEFG(6*MXCENT),
     &              IPPVC(MXCENT),IPEFT(15*MXCENT),IPEFTNTL(15*MXCENT),
     &              IPEFN(3*MXCENT),ORBANA

