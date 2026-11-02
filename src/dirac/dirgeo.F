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

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamgeo */
      SUBROUTINE PAMGEO()
C*****************************************************************************
C
C     Module for geometry optimization
C
C     Written by Joern Thyssen 1997/06/17
C     Last revision : 1997/06/17 jth
C
C*****************************************************************************
!      use dcbgrd_mod, only: DONGRD
!      use dcbpsi_mod, only: DOCCM, DOCIM, DOKRMC, DOLUCIAR, DOLUCT, DOMP2, DORES
!      use dcbham_mod, only: BSS, LEVYLE, SPINFR, ZORA
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCOOR
         use dirac_cfg
#include "implicit.h"
C
#include "dcbham.h"
#include "dcbpsi.h"
#include "dcbgrd.h"
C
      LOGICAL NOGRAD
C
      CALL QENTER('PAMGEO')
C
C     Determine whether analytical or numerical gradient is available
C
      NOGRAD= BSS.OR.ZORA.OR.(SPINFR .AND. .NOT.LEVYLE)
C
      IF(.NOT.DONGRD) THEN
C       We want the analytical gradient, check whether it is implemented
        DONGRD = NOGRAD
     &           .or. domp2
     &           .or. dores
     &           .or. doccm
     &           .or. docim
     &           .or. doluct
     &           .or. dokrmc
     &           .or. doluciar
     &           .or. dirac_cfg_pcm
        IF(DONGRD) THEN
          WRITE(LUPRI,'(/A/A)')
     &     ' * PAMGEO INFO: will use NUMERICAL MOLECULAR GRADIENT',
     &     '         (analytical gradient not implemented for this wave'
     *     //' function type).'
        ENDIF
      ENDIF
      CALL OPTMIN()
C
      CALL QEXIT('PAMGEO')
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck exedir */
      SUBROUTINE EXEDIR()
C*****************************************************************************
C
C     (c) 1997/07/09 by J. Thyssen
C
C     EXEcute DIRac; Based on PAMDRV with a few modifications.
C
C     Last revision: 1997/07/08 jth
C
C*****************************************************************************
!      use dcbgen_mod, only: DOANA, DOPRP, DOPSI, INPTES
!      use dcbdhf_mod, only: LUCYCL
      use priunit_mod, only: LUPRI
!      use dcbham_mod, only: IF, X2CMMF
      use dirac_cfg
#include "implicit.h"
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
C 
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
C
      CALL QENTER('EXEDIR')
C
C     Delete DFCYCL so we do not restart DIIS
C     also delete MCCRES to avoid restarting in relccsd.
C     
      CALL DELFILE(LUCYCL,'DFCYCL')
      CALL DELFILE(16,'MCCRES')
C
C     To avoid readin twice in the 0. iteration skip
C     the readin/setup.
C
C     Input and initialization

C     Always do wave function (implicit .WAVE FUNCTION parameter)
      DOPSI = .TRUE.

C
      IF ( FIRST ) THEN
         FIRST = .FALSE.
      ELSE
         IF (X2CMMF) THEN
!            we need to restart from the 4C vectors
#if defined SYS_WINDOWS
            call system('copy CHECKPOINT.4C CHECKPOINT.h5')
#else
            call system('cp CHECKPOINT.4C CHECKPOINT.h5')
#endif
         END IF
C
C***********************************************************************
C*****   I  N  P  U  T    S  E  C  T  I  O  N   ************************
C***********************************************************************
C
         CALL PAMINP()
         IF ( INPTES ) RETURN
C
C*****************************************************************************
C*****  S E T U P  -  M O D U L E  *******************************************
C*****************************************************************************
C
         CALL PAMSET()
      END IF
C
C
C******************************************************************************
C*****  G E T    W A V E    F U N C T I O N   *********************************
C******************************************************************************
C
C
C     Since there is no fun in doing a geometry optimization
C     without a wave function, PAMPSI is always called.
C
      CALL PAMPSI()
C
C
C*****************************************************************************
C*****  A N A L Y S I S    M O D U L E  **************************************
C*****************************************************************************
C
      IF (DOANA) CALL PAMANA()
C
C
C*****************************************************************************
C*****  R E S P O N S E    M O D U L E  **************************************
C*****************************************************************************
C
C
      IF (DOPRP) CALL PAMPRP()
C
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('EXEDIR')
      RETURN
      END
C  /* Deck rgetccsd */
      FUNCTION RGETCCSD()
      use priunit_mod, only: LUPRI
#include "implicit.h"
C
C     Get energy from RELCCSD
C
      CHARACTER CENERGY*50
#include "../relccsd/results.inc" 
#include "../relccsd/inpt.inc" 
      IF (DOCCSDT) THEN
         CENERGY = 'CCSD(T) energy from RELCCSD'
      ELSE IF (DOCCSD) THEN
         CENERGY = 'CCSD energy from RELCCSD'
      ELSE IF (DOMP2) THEN
         CENERGY = 'MP2 energy from RELCCSD'
      ELSE IF (DOFSPC) THEN
         CENERGY = '(IH-)Fock space CC energy for selected state'
      ELSE
         call quit('no energy from RELCCSD')
      END IF

      RGETCCSD            = ETOT
      if(dofspc) rgetccsd = efs(1)

      WRITE(LUPRI,'(2A,1P,E20.10)')
     &   CENERGY,' = ',rgetccsd
      RETURN
      END
C  /* Deck rgetccsd */
      FUNCTION RGETCI()
C***********************************************************************
C
C     Get energy from DIRRCI
C     Taking the lowest root....
C     
C     Written by T. Saue Dec 9 2004
C***********************************************************************
#include "implicit.h"
C
#include "../dirrci/param.inc"  ! For iterat.inc
#include "../dirrci/iterat.inc" 
      RGETCI = EIGEN2(1)
      RETURN
      END
