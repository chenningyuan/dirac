!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

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
     IF (DSO) THEN
        WRITE(RECP_OUT,'(10X,A)')'>>> RECP error' 
        WRITE(RECP_OUT,'(10X,A)')'DSO is not supported in RECP calculation'
        CALL QUIT('RECP_CHK_PRP: DSO is not supported in RECP') 
     ELSEIF (MGRAD) THEN
        WRITE(RECP_OUT,'(10X,A)')'>>> RECP error' 
        WRITE(RECP_OUT,'(10X,A)')'Molecular gradient is not supported in RECP calculation'
        CALL QUIT('RECP_CHK_PRP: Molecular gradient is not supported in RECP') 
     ELSEIF (SHIELD) THEN
        WRITE(RECP_OUT,'(10X,A)')'>>> RECP error' 
        WRITE(RECP_OUT,'(10X,A)')'NMR is not supported in RECP calculation'
        CALL QUIT('RECP_CHK_PRP: NMR is not supported in RECP') 
     ELSEIF (PVC) THEN
        WRITE(RECP_OUT,'(10X,A)')'>>> RECP error' 
        WRITE(RECP_OUT,'(10X,A)')'PVC is not supported in RECP calculation'
        CALL QUIT('RECP_CHK_PRP: PVC is not supported in RECP') 
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
