C *** Quaternion algebra ***
C
C  This COMMON block provides multiplication rules for
C  quaternion unit vectors 1,i,j and k
C
C  IQMULT(e_1,e_2,e_M) and IQPHAS(e_1,e_2,e_M) gives the
C  resulting unit vector and phase (+/- 1) of a quaternion
C  triple product (e_1)(e_M)(e_2) :
C
C   1 |  1  i  j  k    i |  1  i  j  k
C   ---------------    ---------------
C   1 |  1  i  j  k    1 |  i -1  k -j
C   i |  i -1  k -j    i | -1 -i -j -k
C   j |  j -k -1  i    j | -k -j  i  1
C   k |  k  j -i -1    k |  j -k -1  i
C
C   j |  1  i  j  k    k |  1  i  j  k
C   ---------------    ---------------
C   1 |  j -k -1  i    1 |  k  j -i -1
C   i |  k  j -i -1    i | -j  k  1 -i
C   j | -1 -i -j -k    j |  i -1  k -j
C   k | -i  1 -k  j    k | -1 -i -j -k
C
C   IQSIGN(e_1,FX,TX) gives the sign of the quaternion units under
C                     various operations:
C   FX = 1 : Standard quaternion:   q  = a  + bj
C   FX = 2 : Complex conjugation:   q* = a* - bj
C
C   TX = 1 : Standard quaternion:   q  = a  + bj
C   TX = 2 : i-transform        : -iqi = a  - bj
C   TX = 3 : j-transform        : -jqj = a* + b*j
C   TX = 4 : k-transform        : -kqk = a* - b*j
C  
      COMMON/QALGEBRA/IQMULT(4,4,4),IQPHASE(4,4,4),
     &                IQSIGN(4,2,4)
      CHARACTER QUNIT*1
      COMMON/QUNITS/QUNIT(4)



