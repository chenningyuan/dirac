C*** Labels ***
C    PLABEL(*,1), IATTR(*,1), IPLAB(*,1) for AO basis
C    PLABEL(*,2), IATTR(*,2), IPLAB(*,2) for SO basis
C    NPLAB(0:2): number of labels total/large comp/small comp
      PARAMETER(MAXLAB = 600)
      CHARACTER*12 PLABEL
      COMMON/LABNAM/PLABEL(MAXLAB,2)
      COMMON/LABELS/IATTR(MAXLAB,2),IPLAB(MXCORB,2),
     &              NPLAB(0:2)
