      REAL*8          SPRNUC(3,3,MXCENT),SPRNUC2(3,3,MXCENT),
     &                SPRNUC3(3,3,MXCENT),
     &                ELSPRD(3,MXCOOR), ELSPRP(3,MXCOOR),
     &                SPRDNL(3,3,MXCENT), GTRANT(3,3,MXCENT)

      LOGICAL         LINDET

      COMMON /SPINRO/ SPRNUC,SPRNUC2,SPRNUC3,ELSPRD, ELSPRP, SPRDNL,
     &                GTRANT, LINDET
