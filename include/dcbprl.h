C
C     PROPERTY LABELS:
C
C       PRPLBL  - property label
C       IPRREP  - boson irrep of scalar operator (0 - 7)
C       IPRLTYP - symmetry of scalar operator under matrix transpose (-1,0,1)
C
      PARAMETER ( MAXPRPLBL = 200 )
      CHARACTER   PRPLBL*8, PDOINT*4
      COMMON/XCBPRL/ NPRPLBL,IPRLREP(MAXPRPLBL),IPRLTYP(MAXPRPLBL)
      COMMON/XCBLBL/ PRPLBL(MAXPRPLBL),PDOINT(MAXPRPLBL)
