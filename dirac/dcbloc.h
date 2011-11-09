! ** COMMON block for Molecular orbital localization
!     IPRLOC - print level
!     PRJLOC - localization from projection analysis
!
!#include <dcblab.h>
      LOGICAL PRJLOC,LGFULL,LGDIAG,LGGRAD,LGCHCK
      INTEGER IPRLOC,ITRLOC
      REAL*8  THFULL,THDIAG,THGRAD
      CHARACTER*4 HESLOC
      COMMON/CBILOC/IPRLOC,ITRLOC,THFULL,THDIAG,THGRAD,HESLOC,PRJLOC,   &
     &              LGFULL,LGDIAG,LGGRAD,LGCHCK
