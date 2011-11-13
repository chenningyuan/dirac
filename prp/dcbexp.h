C
C VERSION : $Revision: 1.4 $
C DATE    : $Date: 2003/04/28 15:23:12 $
C FILE    : dcbexp.h
C
C*****************************************************************************
C  DCBEXP:  Expectation values
C
      PARAMETER(MAXEXP=350)
      LOGICAL ORBANA,PRJANA
      COMMON/DCBEXP/IPREXP,NEXPP,LEXPP(MAXEXP),
     &              IPDIPOLE(3),IPQUADRU(6),IPEFG(6*MXCENT),
     &              IPPVC(MXCENT),IPEFT(15*MXCENT),IPEFTNTL(15*MXCENT),
     &              IPEFN(3*MXCENT),ORBANA,PRJANA

