C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ./include/dcbmp2.F
C
C  DCBMP2:  Parameters for MP2 calculation
C
C     MP2_NSTR(OCC/VIR, Fermion irrep)   - number of occ/virt pairs
C     MP2_INDSTR(OCC/VIR, Fermion irrep) - string of active pairs
C
      LOGICAL MP2_MSOUT,MP2_ANTIS,MP2NOP
      COMMON/DCLMP2/MP2_MSOUT,MP2_ANTIS,MP2NOP
      CHARACTER*72 MP2_INDSTR
      COMMON/DCCMP2/MP2_INDSTR(2,2)
      COMMON/DCBMP2/DMP2_VIRTHR,IPRMP2,ISTRMP2,MP2_INTFLG,
     &              MP2_NSTR(2,2)
      COMMON/DCRMP2/SCRMP2,EMP2
