C
C     i     IINFO      DINFO
C    ---   -------    -------
C     01    ITMAC      EMY
C     02    ITMIC      EACTIV
C     03    NZRED      EMCSCF
C     04    NCRED      DEPRED
C     05    NERED      DEACT
C     06    NPRED      RATIO
C     07               BETA
C     08               STPLNG
C     09               RTRUST
C     10               GAMMA
C     11               GCINRM ( GNORM(1) )
C     12               GEENRM ( GNORM(2) )
C     13               GEPNRM ( GNORM(3) )
C     14               GOBNRM ( GNORM(4) )
C     15               GTOTNRM( GNORM(5) )
C     16               TIMMAC
C     17               TIMMIC
C     18               TIMITR
C     19               CI STPLNG
C     20               E STPLNG
C     21               P STPLNG
C     22               LEVEL SHIFT
C
      PARAMETER (LDINFO = 22, LIINFO = 6)
      COMMON /IITINFO/ IINFO(LIINFO)
      COMMON /RITINFO/ DINFO(LDINFO)
