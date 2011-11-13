C
C VERSION : $Revision: 4258 $
C DATE    : $Date: 2005-02-10 20:19:14 +0100 (Thu, 10 Feb 2005) $
C FILE    : dcbexp.h
C
C*****************************************************************************
C  DCBEXP:  Expectation values
C
      PARAMETER(MAXEXP=350)
      LOGICAL ORBANA,PRJANA,PRPCAN
      COMMON/DCBEXP/IPREXP,NEXPP,LEXPP(MAXEXP),
     &              IPDIPOLE(3),IPQUADRU(6),IPEFG(6*MXCENT),
     &              IPPVC(MXCENT),IPEFT(15*MXCENT),IPEFTNTL(15*MXCENT),
     &              IPDSUSC(39),
     &              IPEFN(3*MXCENT),ORBANA,PRJANA,PRPCAN

