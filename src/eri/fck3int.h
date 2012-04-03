C
C     All the coordinates of have eight possible locations due
C     to symmetry, and there can be four centers
C
      LOGICAL AUX1, AUX2, GAB, TCMPAB, TCMPCD, ODTR12,
     &        TKMPAB, TKMPCD
      LOGICAL SPHERA, SPHERB, SPHERC, SPHERD, SPHERAB, SPHERCD
      COMMON /FCK3INT/  AUX1, AUX2, GAB

      COMMON /FCK3COOR/ COOR(3,0:7,4),
     &                  DIST12(3,0:7,0:7,2), DIST122(0:7,0:7,2)


      COMMON /FCK3DIM/  NDIM(4), IABCLSSTART, IABCLSEND,
     &                  ICDCLSSTART, ICDCLSEND

      COMMON /FCK3SYM/  IPQXYZFCK3, IHCBITFCK3, ICCXYZFCK3

      COMMON /FCK3SHR/  SPHERA, SPHERB, SPHERC, SPHERD,
     &                  SPHERAB, SPHERCD, TCMPAB, TCMPCD, ODTR12,
     &                  TKMPAB, TKMPCD
