      SUBROUTINE ZSYMDG(IPRNT)
*
* Number of symmetries for d2h
* Symmetry connecting arrays
*
*
#include "mxpdim.inc"
#include "csm.inc"
#include "symm.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
      NSMSX = NIRR_DG
      NSMDX = NIRR_DG
      NSMST = NIRR_DG
      NSMCI = NIRR_DG
      NSMXT = NIRR_DG
      ITSSX = 1
      ITSDX = 1
      ITSXT = 1

*
      IZERO = 0
      CALL ISETVC(ADASX,IZERO,MXPOBS**2)
      CALL ISETVC(ASXAD,IZERO,2*MXPOBS**2)
      CALL ISETVC(ADSXA,IZERO,2*MXPOBS**2)
*
      DO IREP = 1, NIRR_DG
        DO JREP = 1, NIRR_DG
*A+ A => SX
          ADASX(IREP,JREP) = IDBGMULT(IREP,IADJSYM(JREP))
        END DO
      END DO
*
      DO IREP = 1, NIRR_DG
        DO JREP = 1, NIRR_DG
          DO IJSM = 1, NIRR_DG
            IF(ADASX(IREP,JREP).EQ.IJSM) THEN
              ADSXA(IREP,IJSM) = JREP
              ASXAD(JREP,IJSM) = IREP
            END IF
          END DO
        END DO
      END DO

      CALL ICPMT2(IDBGMULT,SXSXDX,8,8,2*MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(IDBGMULT,SXDXSX,8,8,2*MXPOBS,4*MXPOBS,1)
*
      RETURN
      END
