C*****************************************************************************
C  DCBHAM:  Hamiltonian used for generation of wavefunction
C
C     IPRHAM     - Print level
C     IPOVERLAP  - Pointer to overlap in DCBXPR
C     IPBETAMT   - Pointer to beta matrix in DCBXPR
C     IPMOLFLD   - Pointer to nuclear attraction in DCBXPR
C     IPKINERG   - Pointer to kinetic energy in DCBXPR
C     N1OPER     - Number of additional operators (finite field calculations)
C     IPR1OP     - Pointer to DCBXPR for additional operators
C     ONESYS     - No two-electron operators
C     FREEPJ     - Hamiltonian embedded in free particle projections
C     VEXTPJ     - Hamiltonian embedded in molecular field projections
C     NOSMLV     - No small nuclear attraction integrals
C     URKBAL     - Unrestricted kinetic balance
C     LEVYLE     - Levy-Leblond Hamiltonian
C     SSMTRC     - Scaling factor for SS overlap matrix
C
      PARAMETER(MAX1OPER=10)
      COMMON/DCIHAM/IPRHAM,IPOVRLAP,IPBETAMT,IPMOLFLD,IPKINERG,
     &              N1OPER,IPR1OP(MAX1OPER)
      LOGICAL ONESYS,FREEPJ,VEXTPJ,NOSMLV,URKBAL,LEVYLE
      COMMON/DCLHAM/ONESYS,FREEPJ,VEXTPJ,NOSMLV,URKBAL,LEVYLE
      COMMON/DCRHAM/SSMTRC
