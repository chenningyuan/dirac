! ** COMMON block for Molecular orbital localization
!     IPRLOC - print level
!     PRJLOC - localization from projection analysis
!
!#include <dcblab.h>
      LOGICAL PRJLOC,HESLOC,AUTO
      INTEGER IPRLOC,ITRLOC
      REAL*8  THRLOC,HESTHR
      COMMON/CBILOC/IPRLOC,AUTO,PRJLOC,HESLOC,THRLOC,HESTHR,ITRLOC
