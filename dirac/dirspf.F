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
C     Note that the MO transformation matrix should be defined according 
C     to the SPINFR option (having boson symmetry adapted blocks)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck spfao */
      SUBROUTINE SPFAO (FMO,FAO,WORK,LWORK)
C*****************************************************************************
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
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcbham.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      PARAMETER ( DP5 = 0.50D00 , D0 = 0.00D00, D1 = 1.00D00,
     &            DM1 = -1.00D00, D2 = 2.00D00)
      DIMENSION FMO(*),FAO(*),WORK(LWORK)
C
#include <memint.h>
      CALL QENTER('SPFAO')
C
C     We have no information about the print level : make IPRINT zero
C
      IPRINT = 0
C
C     Get transformation matrix S(-1/2) to go to the orthogonal MO-basis
C
      CALL MEMGET('REAL',KTMAT,N2TMT,WORK,KFREE,LFREE)
      CALL OPNFIL(LUTMAT,'DFTMAT','OLD','PREDHF')
      CALL READT(LUTMAT,N2TMT,WORK(KTMAT))
      CLOSE(LUTMAT,STATUS='KEEP')
C
C     Transform the Fock matrix to MO-basis
C
      CALL MKMOFK2(FMO,FAO,WORK(KTMAT),WORK(KFREE),LFREE)
C
C     Make the Fock matrix spinfree
C
      IF (IPRINT .GE. 6 ) THEN
       DO I = 1, NFSYM
         CALL HEADER('Fock matrix in orth. basis with spin-orbit',-1)
         CALL PRQMAT(FMO(1+I2ORBT(I)),NORB(I),NORB(I),
     &               NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDDO
      END IF
C
      CALL SPFSCF(FMO,WORK(KFREE),LFREE)
C
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
      ISIZE = 0
      DO I = 1,NFSYM
         ISIZE = ISIZE + NZ*NORB(I)*NFBAS(I,0)
      END DO
      CALL MEMGET('REAL',KSAOMO,ISIZE,WORK,KFREE,LFREE)
      CALL MKSAOMO(WORK(KSAOMO),WORK(KTMAT),IPRINT,WORK(KFREE),LFREE)
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
     &            WORK(KSAOMO+I2TMT(I)),NFBAS(I,0),NORB(I),
     &            NZT,IPQTOQ(1,0),
     &            WORK(KSAOMO+I2TMT(I)),NFBAS(I,0),NORB(I),
     &            NZT,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)
      ENDDO
C
      CALL MEMREL('KTMAT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SPFAO')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MKSAOMO */
      SUBROUTINE MKSAOMO(SAOMO,TMAT,IPRINT,WORK,LWORK)
C****************************************************************************
C     PURPOSE:
C        Construct half transformed overlap matrix, Saomo = Sao C
C
C        Written by L. Visscher (based on MKSMOAO by J. Thyssen).
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcbham.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION SAOMO(*),TMAT(*),WORK(*)
C
#include <memint.h>
C
      CALL QENTER('MKSAOMO')
C
      ISIZE = 0
      DO I = 1,NFSYM
         ISIZE = ISIZE + NZ*NFBAS(I,0)*NORB(I)
      END DO
      CALL DZERO(SAOMO,ISIZE)
C
C     Get AO overlap matrix
C
      CALL MEMGET('REAL',KSAO,N2BBASX,WORK,KFREE,LFREE)
      CALL GTOVLX(WORK(KSAO),SSMTRC,WORK(KFREE),LFREE)
C
      IOFF = 1
      DO I = 1,NFSYM
         IF ( NORB(I) .GT. 0 ) THEN
            CALL QGEMM(NFBAS(I,0),NORB(I),NFBAS(I,0),D1,
     &         'N','N',IPQTOQ(1,0),WORK(KSAO+I2BASX(I,I)),
     &         NTBAS(0),NTBAS(0),1,
     &         'N','N',IPQTOQ(1,0),TMAT(I2TMT(I)+1),
     &         NFBAS(I,0),NORB(I),NZT,
     &         D0,IPQTOQ(1,0),SAOMO(IOFF),
     &         NFBAS(I,0),NORB(I),NZ)
C
            IF ( IPRINT .GE. 10 ) THEN
               CALL HEADER('PREDHF: Right index transf. overlap'//
     &            ' matrix Smomo',-1)
               CALL PRQMAT(SAOMO(IOFF),NFBAS(I,0),NORB(I),
     &            NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
               CALL FLSHFO(LUPRI)
            END IF
         END IF
      IOFF = IOFF + NZ*NFBAS(I,0)*NORB(I)
      END DO
C
      CALL MEMREL('SAOMO',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MKSAOMO')
      RETURN
      END
