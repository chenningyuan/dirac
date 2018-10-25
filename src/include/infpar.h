#ifdef COMMENT
! -- infpar.h --
!     my_MPI_INTEGER is used both in .c and .F routines in MPI calls
!        so we can handle "-i8" compilations on 32-bit machines,
!        using INT_STAR8 /Jan-2007 hjaaj
#endif
#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

#if defined(__CVERSION__)
#define nodtot numnod
extern struct common_diracinfpar {
    double totval;
#if defined (INT_STAR8)
    long iprpar, ntask, ncode, ndegdi, mynum, mytid;
    long numnod, icount, mparid, node, nodes, iparex;
    long max_omp_threads;
    long parher, slave, timing, debug, pario;
#else
    int  iprpar, ntask, ncode, ndegdi, mynum, mytid;
    int  numnod, icount, mparid, node, nodes, iparex;
    int  max_omp_threads;
    int  parher, slave, timing, debug, pario;
#endif 
} diracinfpar_;
#else /* __CVERSION__ */
!  -- File: infpar.h for DIRAC; specific information for parallel calculations
!     .. mpi variables (used to be in an infmpi.h /hjaaj Jan 2003):
!        NODTOT = NUMNOD = number of slaves = number of nodes - 1 (from MPI_Comm_rank)
!        MYTID  = mpi ID of this process (from MPI_Comm_rank)
      REAL*8 TOTWAL
      INTEGER IPRPAR, NTASK, NCODE, NDEGDI, MYNUM, MYTID, NODTOT
      INTEGER NUMNOD, ICOUNT,MPARID,NODE,NODES,IPAREX
      INTEGER MAX_OMP_THREADS ! for 
      LOGICAL PARHER, SLAVE, TIMING, DEBUG, PARIO
      COMMON /DIRACINFPAR/ TOTWAL                                        &
     &       ,IPRPAR, NTASK, NCODE, NDEGDI, MYNUM, MYTID, NODTOT         &
     &       ,ICOUNT, MPARID, NODE, NODES, IPAREX                        &
     &       ,MAX_OMP_THREADS                                            &
     &       ,PARHER, SLAVE, TIMING, DEBUG, PARIO

      EQUIVALENCE (NODTOT, NUMNOD)
!     ... for some obscure reason NODTOT in Dalton is called NUMNOD in Dirac.
!         This equivalence is important such that both e.g. abacus routines using NODTOT
!         and dirac routines using NUMNOD are working correctly together /hjaaj Dec 07

! -- end of infpar.h --
#endif
