C     LMOORD: Reorder orbitals?
C     MXREORD: maximum number of reorderings
C     IREORD, IMOORD: orbital reordering
      PARAMETER ( MXREORD = 50 )
      LOGICAL LMOORD
      COMMON /CBIREO/ IMOORD(MXREORD*2,2),IREORD(2)
      COMMON /CBLREO/ LMOORD
