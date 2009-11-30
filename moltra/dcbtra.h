C*** Information for 4-index transformation program **
C
C     Length of common block : keep up to date !!
      LOGICAL TRA_ANTIS,TRA_MSOUT,NOPAIR,NO2IND,NO4IND,PRPTRA,RCORBS
      LOGICAL FKRMC, NOMDCINT, MDCSCAT, TRA_ASCII, PRPSYA
c aspg, 20090424, for fine-sort of 1HT
      logical do_finesort,set_auto_finesort
c
      CHARACTER*72 TRA_INDSTR,TRA2_INDSTR,TRA4_INDSTR,TRA_CORSTR,
     &             TRA2_INDPRP
      PARAMETER (NDCBTRI=267,NDCBTRR=258,NDCBTRL=12)
      COMMON/DCBTC/TRA_INDSTR(2),TRA2_INDSTR(2,2),TRA4_INDSTR(4,2),
     &             TRA_CORSTR(2), TRA2_INDPRP(2)
      COMMON/DCBTRI/IPRTRA,ISTRAT,ITRA_INTFLG,ITRA_INTFL2,ITRA_INTFL4,
     &              NQQCLASS,ISAME(4),ICLARR(4,4,16),
     &              NFPCK12(2),NFPCK34(2),NFPCK12T,NFPCK34T,
     &              IFPCK12(2,2),IFPCK34(2,2),N2EFIL,NTKL,IPAR4BS
      COMMON/DCBTRL/TRA_ANTIS,TRA_MSOUT,NOPAIR,NO2IND,NO4IND,PRPTRA,
     &              RCORBS,FKRMC,NOMDCINT,MDCSCAT,TRA_ASCII,
     &              do_finesort,set_auto_finesort, PRPSYA
      COMMON/DCBTRR/THROUT,SCRTRA,RINVICLARR(4,4,16)
C     Unit for ASCII interface file
      PARAMETER (LUASCII=74)
C     Units for MOLFDIR interface files
      PARAMETER (LUMLF1=75,LUMLF2=76,LUMLF3=76)
C     Units for scratch files (LUTRA1 is an offset (-1) for 2 files !)
      PARAMETER (LUTRA1=78,LUTRA2=81)
#include <comdis.h>


