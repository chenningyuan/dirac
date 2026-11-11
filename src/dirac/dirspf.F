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

C     Contains the modifications to activate the spinfree option in 
C     the code
C
C
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck spfscf */
      SUBROUTINE SPFSCF(FMO,WORK,LWORK)
C***********************************************************************
C
C     Zero out the spin-dependent part of the Fock matrix.
C     On input:
C       FMO - Full Fock matrix in MO-basis
C     On output:
C       FMO - Spinfree Fock matrix in MO-basis
C
C     Note that the MO transformation matrix should be defined according 
C     to the SPINFR option (having boson symmetry adapted blocks)
C
C     Written by L.Visscher, august 1998
C
C***********************************************************************
!      use dcborb_mod, only: I2ORBT, NBORB, NORB
      use dgroup_mod, only: NFSYM, NZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I,J
C
#include "dcborb.h"
      REAL(dp), PARAMETER :: DP5 = 0.50D00, D0 = 0.00D00, D1 = 1.00D00,
     &                       DM1 = -1.00D00, D2 = 2.00D00
C
      DIMENSION FMO(*),WORK(LWORK)
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
C
C Zero out the coupling blocks between different boson irreps and the
C imaginary or q-imaginary blocks. With a 1-1 kinetic balance basis this
C will correspond to the spin-free formalism defined by Dyall.
C
      DO IFRP = 1, NFSYM
         NORB2 = NORB(IFRP)*NORB(IFRP)
C
C        Zero the imaginary or q-imaginary matrices completely
C        -----------------------------------------------------
C
         CALL DZERO(FMO(1+I2ORBT(IFRP)+NORB2),NORB2*(NZ-1))
C
C        Zero off-diagonal blocks between different boson irreps
C        -------------------------------------------------------
C
         NBRP = 4 / NZ
         IJ = I2ORBT(IFRP)
         JI = I2ORBT(IFRP)
         DO JZ = 2, 1, -1
          DO JSYM = 1, NBRP
           NORBJ = NBORB(JSYM,IFRP,JZ)
            DO J = 1, NORBJ
             JI = JI + 1
             DO IZ = 2, 1, -1
              DO ISYM = 1, NBRP
               NORBI = NBORB(ISYM,IFRP,IZ)
               DO I = 1, NORBI
                IJ = IJ + 1
                IF (ISYM.NE.JSYM) THEN
                    FMO(IJ) = D0
                    FMO(JI) = D0
                ENDIF
                JI = JI + NORB(IFRP)
               ENDDO
              ENDDO
             ENDDO
             JI = JI - NORB2
           ENDDO
          ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck spfao */
      SUBROUTINE SPFAO(FMO,FAO,WORK,LWORK)
C***********************************************************************
C
C     Zero out the spin-dependent part of the Fock matrix.
C     On input:
C       FAO - Full Fock matrix in AO-basis
C     On output:
C       FAO - Spinfree Fock matrix in AO-basis
C
C     Note that the MO transformation matrix should be defined according 
C     to the SPINFR option (having boson symmetry adapted blocks)
C
C     Written by L.Visscher, march 1999
C
C***********************************************************************
      use dcbbas_mod, only: I2BASX, I2TMT, N2TMT, NFBAS, NTBAS, NZT
      use priunit_mod, only: LUPRI
      use dcbgen_mod, only: LUTMAT
!      use dcborb_mod, only: I2ORBT, NORB
      use dgroup_mod, only: IPQTOQ, NFSYM, NZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
C
!#include "dcbdhf.h"
!#include "dcbham.h"
#include "dcborb.h"
      REAL(dp), PARAMETER :: DP5 = 0.50D00, D0 = 0.00D00, D1 = 1.00D00,
     &                       DM1 = -1.00D00, D2 = 2.00D00
      DIMENSION FMO(*),FAO(*),WORK(LWORK)
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
      CALL QENTER('SPFAO')
C
C     We have no information about the print level : make IPRINT zero
C
      IPRINT = 0

      CALL MEMGET('REAL',KBUF,N2TMT,WORK,KFREE,LFREE)
C
C     Get transformation matrix S(-1/2) to go to the orthogonal MO-basis
C
      CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','SPFAO')
      CALL READT(LUTMAT,N2TMT,WORK(KBUF))
C
C     Transform the Fock matrix to spinfree MO-basis 
C
      CALL MKMOFK2(FMO,FAO,WORK(KBUF),WORK(KFREE),LFREE)
      IF (IPRINT .GE. 3 ) THEN
       DO I = 1, NFSYM
         CALL HEADER('Fock matrix in orth. basis without spin-orbit',-1)
         CALL PRQMAT(FMO(1+I2ORBT(I)),NORB(I),NORB(I),
     &               NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDDO
      END IF
C
C     Get transformation matrix S(+1/2) to go back to AO-basis
C
      CALL READT(LUTMAT,N2TMT,WORK(KBUF))
      CLOSE(LUTMAT,STATUS='KEEP')
C
C     Tranform back to AO-basis
C
      DO I = 1,NFSYM
         IF(NORB(I).NE.0) CALL QTRANS('MOAO','S',D0,
     &            NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &            FAO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,
     &            IPQTOQ(1,0),
     &            FMO(1+I2ORBT(I)),NORB(I),NORB(I),NZ,
     &            IPQTOQ(1,0),
     &            WORK(KBUF+I2TMT(I)),NFBAS(I,0),NORB(I),
     &            NZT,IPQTOQ(1,0),
     &            WORK(KBUF+I2TMT(I)),NFBAS(I,0),NORB(I),
     &            NZT,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)
      ENDDO
C
      CALL MEMREL('SPFAO',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SPFAO')
C
      RETURN
      END
