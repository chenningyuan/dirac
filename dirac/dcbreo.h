C     LMOORD: Reorder orbitals?
C     MXREORD: maximum number of reorderings
C     IREORD, IMOORD: orbital reordering
      PARAMETER ( MXREORD = 200 )
      LOGICAL LMOORD
      COMMON /CBIREO/ IMOORD(MXREORD,2),IREORD(2)
      COMMON /CBLREO/ LMOORD
C
C     LMOORF: Post-DHF Reorder orbitals?
C     IREORF, IMOORF: orbital reordering
      LOGICAL LMOORF
      COMMON /CBI2RE/ IMOORF(MXREORD,2),IREORF(2)
      COMMON /CBL2RE/ LMOORF
