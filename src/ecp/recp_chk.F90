!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
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

! linear symmetry is not working yet, switch it off
  CALL RECP_CHK_SYM

END SUBROUTINE RECP_CHK_MAIN


SUBROUTINE RECP_CHK_PRP
! check inconsistency in property options 
  use mxcent_mod
  IMPLICIT NONE
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
  use mxcent_mod
  IMPLICIT NONE
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

SUBROUTINE RECP_CHK_SYM
! switch off linear symmetry if this is active
  IMPLICIT NONE
#include "inc_print.h"
#include "dcbham.h"
#include "dgroup.h"
  IF (LINEAR.AND.ECPCALC) THEN
!       * Switch off linear symmetry
!       WRITE(RECP_OUT,'(10X,A)')'>>> RECP incompatibility: '
!       WRITE(RECP_OUT,'(10X,A)')'Linear symmetry is switched off' 
!       LINEAR = .FALSE.

!       * Keep linear symmetry and show warning (for auto-symmetry/Dinfh/Cinfv)
        WRITE(RECP_OUT,'(X,A,A)')'* Warning: Symmetries higher than D2h in RECP integral ', &
                                 'is not fully tested.'
        WRITE(RECP_OUT,'(X,A)')  '           The calculation results can be erroneous.'
        WRITE(RECP_OUT,'(X,A)')  '           Use auto-symmetry with care.'
  ENDIF
       
END SUBROUTINE RECP_CHK_SYM

END MODULE RECP_CHK
