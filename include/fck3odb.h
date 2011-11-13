C
C   See the description in fck3odb.F
C
      LOGICAL DOSCREENCS, SYMDENS
      COMMON /FCK3GAB/  DOSCREENCS
      COMMON /FCK3ODBTCH/ SYMDENS, NODBCLS, NODBTCH,
     &                    NPRFA, NCTFA, NPRFB, NCTFB,
     &                    NCSQ1, NCSQ2, JTOP,NRTOP,
     &                    MAXODNHK, MAXODKCK, MAXODKHK,
     &                    MAXODPRF, MAXODCTF,
     &                    MAXFCLSPRF,NFTUVCLS,NFTUVBTCH

      COMMON /FCK3COOR12/ COOR1(3,0:7),COOR2(3,0:7),
     &                    R12(3,0:7,0:7),R122(0:7,0:7),
     &                    INDREP(2,0:7)
