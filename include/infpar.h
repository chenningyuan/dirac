#ifdef __CVERSION__
#define MAXNOD 100
struct common_infpar {
    double totval;
    int icount, kprint, mtids[MAXNOD], mparid, mytid;
    int node, nodes, numnod, ntask, ncode, parher, timing;
} infpar_;
#else /* __CVERSION__ */
C     .. mpi variables (used to be in an infmpi.h /hjaaj Jan 2003):
C        NUMNOD = number of slaves = number of nodes - 1 (from MPI_Comm_rank)
C        MYTID  = mpi ID of this process (from MPI_Comm_rank)
     PARAMETER (MAXNOD=100, MAXCL2=10000, MAXTSK=MXSHEL*(MXSHEL+1)/2)
      REAL*8 TOTWAL
      LOGICAL PARHER, TIMING, SLAVE, DEBUG
      COMMON /INFPAR/ TOTWAL,
     &                ICOUNT, IPRPAR, MTIDS(MAXNOD), MPARID, MYTID,
     &                NODE, NODES, NUMNOD, NTASK, NCODE, PARHER, TIMING,
     &                SLAVE, NODTOT, NDEGDI, DEBUG
#endif
