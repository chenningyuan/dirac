C  DCBMP2:  Parameters for MP2 calculation
C
C     MP2_NSTR(OCC/VIR, Fermion irrep)   - number of occ/virt pairs
C     MP2_INDSTR(OCC/VIR, Fermion irrep) - string of active pairs
C
      LOGICAL MP2_MSOUT,MP2_ANTIS,NOPAIR
      COMMON/DCLMP2/MP2_MSOUT,MP2_ANTIS,NOPAIR
      CHARACTER*72 MP2_INDSTR
      COMMON/DCCMP2/MP2_INDSTR(2,2)
      COMMON/DCBMP2/DMP2_VIRTHR,IPRMP2,ISTRMP2,MP2_INTFLG,
     &       MP2_NSTR(2,2)
      COMMON/DCCMPE/EMP2

