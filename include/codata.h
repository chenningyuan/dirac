C
C  codata.h revised 1999/04/29 and combined with units.h
C
C
C  From
C  "The 1986 CODATA Recommended Values of the Fundamental Physical Constants"
C  E.R.Cohen and B.N.Taylor, J.Phys.Chem.Ref.Data, 17,1795 (1988)
C  (NB! typo for XTKAYS in Table 4 in Cohen & Taylor 1988)
C
C -- Fundamental constants, atomic units, Avogadro's constant
#include <alphac.h>
C
      PARAMETER ( XTANG = 0.529 177 249 D0, XFAMU  = 1822.88852 D0,
     &            ECHARGE = 1.602 177 33D-19, HBAR = 1.054 572 66D-34)
      PARAMETER ( XFMOL = 6.022 136 7 D23 )
C -- conversion from hartree (au) to:
      PARAMETER ( XTJ  = 4.359 748 2 D-18,  XTKAYS = 219 474.63067 D0,
     *            XTHZ = 6.579 683 899 9 D15, XTEV = 27.211 396 D0,
     &            XKJMOL = XTJ*XFMOL*1.D-3, XKCMOL = XKJMOL/4.184 D0,
     &            XTKJML = XKJMOL,          XTKCML = XKCMOL,
     *            XTNM = 1.D7/XTKAYS,       XAJOUL = 1.0 D18*XTJ)
C
C -- other
      PARAMETER ( AUTIME = HBAR/XTJ, XFSEC = AUTIME)
      PARAMETER ( XTKMML = 974.864 D0)
      PARAMETER ( TESLA=(XTANG*XTANG*ECHARGE/HBAR)*1.D-20 )
C
