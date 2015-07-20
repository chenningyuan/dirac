 
 
!     Register the status of the independent parts of the code
!
!     Status given by integer : 
!     0 - Never asked for
!     1 - On menu but not yet called
!     2 - Running
!     3 - Completed, restart information available
!     4 - Completed, no restart information
!
!     The indices of the different tasks are :
!      1 - Integral sort up to VOVO integrals
!      2 - Integral sort of VOVV and VVVV integrals
!      3 - Construction of Fock matrix
!      4 - MP2  energy and T1 + T2 amplitude calculation
!      5 - CCSD energy and T1 + T2 amplitude calculation
!      6 - CCSD(T) energy calculation
!      7 - 30 reserved for future use
!
      PARAMETER (NSTATS=30)
      COMMON/STATUS/IMSTAT(NSTATS)
!
!     Gather results of the energy calculations
!     ETOT is the energy at the current highest level of theory
!     EFS is the energy of the lowest/selected eigenstate in the FSCC run (last sector)
!
      COMMON/ENERG/ECORE,ETOT,ESCF,EMP2,ECCSD,ECCSDIM,ET1,ET2,ET3,EFS(2)
      REAL*8 ECORE,ETOT,ESCF,EMP2,ECCSD,ECCSDIM,ET1,ET2,ET3,EFS
