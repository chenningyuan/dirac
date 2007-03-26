#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

	  
#ifdef __CVERSION__
#define MAXNOD 200
struct common_infpar {
    double totval;
    int icount, kprint, mtids[MAXNOD], mparid, mytid;
    int node, nodes, numnod, ntask, ncode, parher, timing;
} infpar_;
#else /* __CVERSION__ */
!     my_MPI_INTEGER is used both in .c and .F routines in MPI calls
!        so we can handle "-i8" compilations on 32-bit machines,
!        using INT_STAR8 /Jan-2007 hjaaj
!     .. mpi variables (used to be in an infmpi.h /hjaaj Jan 2003):
!        NUMNOD = number of slaves = number of nodes - 1 (from MPI_Comm_rank)
!        MYTID  = mpi ID of this process (from MPI_Comm_rank)
      INTEGER MAXNOD, MAXCL2
      PARAMETER (MAXNOD = 200, MAXCL2 = 10000)
      REAL*8 TOTWAL
      LOGICAL PARHER, TIMING, SLAVE, DEBUG
      INTEGER ICOUNT, IPRPAR, MTIDS(MAXNOD), MPARID, MYTID, NODE, NODES,&
     &        NUMNOD, NTASK, NCODE, NODTOT, NDEGDI, MYNUM
      COMMON /INFPAR/ TOTWAL, ICOUNT, IPRPAR, MTIDS, MPARID, MYTID,     &
     &                NODE, NODES, NUMNOD, NTASK, NCODE, PARHER, TIMING,&
     &                SLAVE, NODTOT, NDEGDI, DEBUG, MYNUM
#if defined (VAR_MPI)
#if defined (INT_STAR8)
!     integer array ISTAT contains MPI_SOURCE information.
!     Proper use of ISTAT on 64-bit machines in 
!     combination with INT_STAR8 requires explicit declaration 
!     as INTEGER*4 /March-2007 sk 
      INTEGER*4 ISTAT
#endif
#endif

#endif
