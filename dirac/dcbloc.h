! ** COMMON block for Molecular orbital localization
!     IPRLOC - print level
!     PRJLOC - localization from projection analysis
!
!#include <dcblab.h>
      LOGICAL PRJLOC,LGFULL,LGDIAG
      INTEGER IPRLOC,ITRLOC
      REAL*8  THFULL,THDIAG
      CHARACTER*4 HESLOC
      COMMON/CBILOC/IPRLOC,ITRLOC,THFULL,THDIAG,HESLOC,PRJLOC,          &
     &              LGFULL,LGDIAG
