      SUBROUTINE CAAB_INCLUDED(IEXOP,INCLUDED_CAAB)
*
* Is excitation operator CAAB included in the 
* list of excitations
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
*. Specific input
      INTEGER IEXOP(NGAS,4)
*
      INCLUDED_CAAB = 1
* Creation alpha
      CALL OCC_INCLUDED(IEXOP(1,1),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
* Creation beta
      CALL OCC_INCLUDED(IEXOP(1,2),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
* Annihilation alpha
      CALL OCC_INCLUDED(IEXOP(1,3),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
* Annihilation beta
      CALL OCC_INCLUDED(IEXOP(1,4),INCLUDED)
      IF(INCLUDED.EQ.0) INCLUDED_CAAB = 0
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' CAAB excitation operator '
       CALL WRT_SPOX_TP(IEXOP,1)
       WRITE(6,*) ' INCLUDED_CAAB = ', INCLUDED_CAAB
      END IF
*
      RETURN
      END

    
