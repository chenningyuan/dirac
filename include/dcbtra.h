C*** Information for 4-index transformation program **
C
C     Length of common block : keep up to date !!
      LOGICAL TRA_ANTIS,TRA_MSOUT,NOPAIR,NO2IND,NO4IND,PRPTRA
      CHARACTER*72 TRA_INDSTR,TRA2_INDSTR,TRA4_INDSTR,TRA_CORSTR
      PARAMETER (NDCBTRI=269)
      COMMON/DCBTC/TRA_INDSTR(2),TRA2_INDSTR(2,2),TRA4_INDSTR(4,2),
     &             TRA_CORSTR(2)
      COMMON/DCBTRI/IPRTRA,ISTRTRA,ITRA_INTFLG,ITRA_INTFL2,ITRA_INTFL4,
     &              ISAME(4),ICLARR(4,4,16),IQCMPL(4)
      COMMON/DCBTRL/TRA_ANTIS,TRA_MSOUT,NOPAIR,NO2IND,NO4IND,PRPTRA
      COMMON/DCBTRR/THROUT
C     Units for MOLFDIR interface files
      PARAMETER (LUMLF1=75,LUMLF2=76,LUMLF3=76)
C     Units for scratch files (LUTRA1 is an offset (-1) for 2 files !)
      PARAMETER (LUTRA1=78,LUTRA2=81)
#include <comdis.h>

