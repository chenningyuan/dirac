      SUBROUTINE ZSYM2_REL(IPRNT)
*
* Symmetry connecting arrays
*
* ======
*. Input
* ======
*
*./LUCINP/
C     INTEGER PNTGRP,CITYP,EXTSPC
C     PARAMETER(MXPIRR = 20)
C     PARAMETER(MXPOBS = 20 )
C     PARAMETE
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C    &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30
#include "mxpdim.inc"
#include "lucinp.inc"
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONABR/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*./CSM
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
* ======
*.Output
* ======
*
*./CSMPRD/
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
**. ADASX,ASXAD,ADSXA
      CALL ISETVC(ADASX,0,MXPOBS**2)
      CALL ISETVC(ASXAD,0,2*MXPOBS**2)
      CALL ISETVC(ADSXA,0,2*MXPOBS**2)
*
      DO 100 ISM = 1, NSMOB
C       MLSM_REL(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM_REL(IML,IPARI,ISM,'OB',2)
        DO 90 JSM = 1, NSMOB
          CALL MLSM_REL(JML,JPARI,JSM,'OB',2)
*.a+ i a j symmetry
          IJML = IML - JML
          IF((IPARI.EQ.1.AND.JPARI.EQ.1).OR.
     &       (IPARI.EQ.2.AND.JPARI.EQ.2)) THEN
            IJPARI = 1
          ELSE
            IJPARI = 2
          END IF
          IJSM = (IJPARI-1)*NMLSX + IJML - MNMLSX + 1
          ADASX(ISM,JSM) = IJSM
          ASXAD(JSM,IJSM) = ISM
          ADSXA(ISM,IJSM) = JSM
   90   CONTINUE
  100 CONTINUE
*.SXSXDX,SXDXSX
      DO 200 ISX = 1, NSMSX
C       MLSM_REL(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM_REL(IML,IPARI,ISX,'SX',2)
        DO 190 JSX = 1, NSMSX
          CALL MLSM_REL(JML,JPARI,JSX,'SX',2)
          IF((IPARI.EQ.1.AND.JPARI.EQ.1).OR.
     &       (IPARI.EQ.2.AND.JPARI.EQ.2)) THEN
            IJPARI = 1
          ELSE
            IJPARI = 2
          END IF
          IJML = IML + JML
          IJSM = (IJPARI-1)*NMLDX+IJML - MNMLDX + 1
          SXSXDX(ISX,JSX) = IJSM
          SXDXSX(ISX,IJSM) = JSX
  190   CONTINUE
  200 CONTINUE
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' ADASX '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ADASX,NSMOB,NSMOB,MXPOBS,MXPOBS)
         WRITE(6,*) ' ASXAD '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ASXAD,NSMOB,NSMSX,MXPOBS,2*MXPOBS)
         WRITE(6,*) ' ADSXA '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ADSXA,NSMOB,NSMSX,MXPOBS,2*MXPOBS)
         WRITE(6,*) ' SXSXDX'
         WRITE(6,*) ' ======'
         CALL IWRTMA(SXSXDX,NSMSX,NSMSX,2*MXPOBS,2*MXPOBS)
         WRITE(6,*) ' SXDXSX'
         WRITE(6,*) ' ======'
         CALL IWRTMA(SXDXSX,NSMSX,NSMDX,2*MXPOBS,4*MXPOBS)
      END IF
*
      RETURN
      END
