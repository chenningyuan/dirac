MODULE RECP_OUTPUT

CONTAINS

! ---------------------------------------------------------------

SUBROUTINE OUTPUT_TITLE
!  OUTPUT : print title
   IMPLICIT NONE
#include "inc_print.h"

   IF (PRINT_LEVEL .GE. 5) THEN
      WRITE(RECP_OUT,*) '   Relativistic Effective Core Potential Integral routine'
      WRITE(RECP_OUT,*) '   This routine is based on ARGOS integral'
      WRITE(RECP_OUT,*) '   with the permission of R. M. Pitzer (Ohio State University)'
      WRITE(RECP_OUT,*) '   and maintained by Y. C. Park (KAIST, youngc@kaist.ac.kr)'
   ENDIF
END SUBROUTINE OUTPUT_TITLE

! ---------------------------------------------------------------

SUBROUTINE OUTPUT_LOC(RECP_LOC,RECP_STATUS)
!  OUTPUT : print location
   IMPLICIT NONE
#include "inc_print.h"
   CHARACTER(LEN=*) :: RECP_LOC
   CHARACTER*1  :: RECP_STATUS

!  FIXIT : Count number of RECP_LOC and print RECP_LOC(1:#RECP_LOC)

   IF (RECP_DBG .GE. 5) THEN
      IF(RECP_STATUS .EQ. 'E') WRITE(RECP_OUT,*) ' - RECP : Entering.... ', RECP_LOC
      IF(RECP_STATUS .EQ. 'X') WRITE(RECP_OUT,*) ' - RECP : Exiting..... ', RECP_LOC
      IF(RECP_STATUS .EQ. 'H') WRITE(RECP_OUT,*) ' - RECP : Here we are, ', RECP_LOC
   ENDIF
END SUBROUTINE OUTPUT_LOC

! ---------------------------------------------------------------

SUBROUTINE OUTPUT_ERROR(RECP_LOC,ERR_MESSAGE,ERR_VALUE,ERR_TYPE)
  IMPLICIT NONE
#include "inc_print.h"
  CHARACTER(LEN=*) :: RECP_LOC, ERR_MESSAGE
  INTEGER          :: ERR_VALUE, ERR_TYPE
  IF (ERR_TYPE.EQ.3) THEN    
     WRITE(RECP_OUT,*) RECP_LOC,':',ERR_MESSAGE,':',ERR_VALUE
     STOP
  ENDIF

END SUBROUTINE OUTPUT_ERROR

END MODULE RECP_OUTPUT
