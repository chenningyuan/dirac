C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C     Contains the modifications to activate the spinfree option in the code
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dirspf.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck spfscf */
      SUBROUTINE SPFSCF(FMO,WORK,LWORK)
C*****************************************************************************
C
C     Zero out the spin-dependent part of the Fock matrix.
C     On input:
C       FMO - Full Fock matrix in MO-basis
C     On output:
C       FMO - Spinfree Fock matrix in MO-basis
C
C     Note that the MO matrix should be defined according to the SPINFR option
C     in 
C
C     Written by L.Visscher, august 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      PARAMETER ( DP5 = 0.50D00 , D0 = 0.00D00, D1 = 1.00D00,
     &            DM1 = -1.00D00, D2 = 2.00D00)
C
      DIMENSION FMO(*),WORK(LWORK)
C
#include <memint.h>
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
