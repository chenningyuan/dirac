C
C     i    IINFO      DINFO
C    ---  -------    -------
C     1    ITMAC      EMY
C     2    ITMIC      EACTIV
C     3    NREDL      EMCSCF
C     4               DEPRED
C     5               DEACT
C     6               RATIO
C     7               BETA
C     8               STPLNG
C     9               RTRUST
C    10               GAMMA
C    11               GCINRM ( GNORM(1) )
C    12               GEENRM ( GNORM(2) )
C    13               GEPNRM ( GNORM(3) )
C    14               GOBNRM ( GNORM(4) )
C    15               GTOTNRM( GNORM(5) )
C    16               TIMMAC
C    17               TIMMIC
C    18               TIMITR
C
      PARAMETER (LDINFO = 20, LIINFO = 20)
      COMMON /IITINFO/ IINFO(LIINFO)
      COMMON /RITINFO/ DINFO(LDINFO)
