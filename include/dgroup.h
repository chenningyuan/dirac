C
C     NZ - number of matrices (1,2,4)
C     NFSYM - number of fermion irreps
C     NBSYM - number of boson irreps
C     JSPINR(IM,IC,IFRP) - boson irrep for a given spinor position
C     JBTOF(IBRP,IC) - points to fermion irrep for given boson irrep/component
C     JFSYM(NZ,IFRP) - points to boson symmetries associated with given fermion
C                      ircop
C     IQDEF(4)
C     JQBAS(IBRP,IC) - points to quaternionic position
C     IPQTOQ(4,IBRP) - pointer from packed quaternion to quaternion unit
C     IQTOPQ(4,IBRP) - pointer from quaternion to packed quaternion unit
C     IHQMAT(4,IH, ) - matrix symmetry(symmetric/antisymmetric)
C     IRQMAT(4,IBRP)    - irrep of matrix
C
      CHARACTER FREP*3
      COMMON/DNAME/FREP(2)
      PARAMETER (NDGROUP=192)
      COMMON/DGROUP/NZ,NFSYM,NBSYM,IFDIRP(4,2),JBTOF(0:7,2),
     &              JMROI(4),JSPINR(4,2,2),JQBAS(0:7,2),JFSYM(4,2),
     &              IPQTOQ(4,0:7),IQTOPQ(4,0:7),IQPH(4,2),
     &              IHQMAT(4,-1:1),IRQMAT(4,0:7),IQDEF(4),ITOT
C
C *** 4-component operators : Matrix part ***
C
C     M = I_4,iA_z,iA_y,iA_x
C
C     JM4REP(0:7) - bosonirrep of matrix operator
C     JM4POS(0:7) - quaternionic position of matrix operator (1,2(i),3(j),4(k))
C     JM4FAS(0:7) - sign of matrix operator (+/- 1)
C
C     OPERATORS (18 types):
C
C     MCMP  - number of components for a given operator
C     JM4   - pointer to matrix operator for a given component
C     JCOM  - coefficient of component
C
      COMMON/ONEOPS/MCMP(18),JM4(3,18),JCOM(3,18),
     &              JM4REP(0:7),JM4POS(0:7),JM4FAS(0:7)

