MODULE RECP_CHK
CONTAINS

SUBROUTINE RECP_CHK_MAIN
  IMPLICIT NONE
#include "inc_print.h"


! check inconsistency in property options 
  CALL RECP_CHK_PRP

! check inconsistency in optimize options
  CALL RECP_CHK_OPT

END SUBROUTINE RECP_CHK_MAIN


SUBROUTINE RECP_CHK_PRP
! check inconsistency in property options 
  IMPLICIT NONE
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "inc_print.h"
  IF (DOPRP) THEN
     IF (SHIELD) THEN
        WRITE(RECP_OUT,'(10X,A)')'>>> RECP error' 
        WRITE(RECP_OUT,'(10X,A)')'NMR is not supported in RECP calculation'
        CALL QUIT('RECP_CHK_PRP: NMR is not supported in RECP') 
     ENDIF
  ENDIF
END SUBROUTINE RECP_CHK_PRP


SUBROUTINE RECP_CHK_OPT
! check inconsistency in optimize options
  IMPLICIT NONE
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbgrd.h"
#include "inc_print.h"

  IF (OPTIMI) THEN
!    Currently RECP only support numerical gradient in geo-opt
     IF (.NOT.DONGRD) THEN
        WRITE(RECP_OUT,'(10X,A)')'>>> RECP error' 
        WRITE(RECP_OUT,'(10X,A)')'Only numerical gradient in geometry optimization' 
        WRITE(RECP_OUT,'(10X,A)')'is supported for RECP calculation' 
        WRITE(RECP_OUT,'(10X,A)')'Use NUMGRA in optimization'
        CALL QUIT('RECP_CHK_OPT: Use NUMGRA in optimization') 
     ENDIF
  ENDIF
END SUBROUTINE RECP_CHK_OPT

END MODULE RECP_CHK
