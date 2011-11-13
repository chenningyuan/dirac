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
#define MAXNOD 200
#define MAXCL2 10000
#define nodtot numnod
extern struct common_diracinfpar {
    double totval;
#if defined (INT_STAR8)
    long iprpar, ntask, ncode, ndegdi, mynum, mytid;
    long numnod, icount, mtids[MAXNOD+1], mparid, node, nodes, iparex;
    long parher, slave, timing, debug, pario;
#else
    int  iprpar, ntask, ncode, ndegdi, mynum, mytid;
    int  numnod, icount, mtids[MAXNOD+1], mparid, node, nodes, iparex;
    int  parher, slave, timing, debug, pario;
#endif 
} diracinfpar_;
#else /* __CVERSION__ */
!  -- File: infpar.h for DIRAC; specific information for parallel calculations
!     .. mpi variables (used to be in an infmpi.h /hjaaj Jan 2003):
!        NODTOT = NUMNOD = number of slaves = number of nodes - 1 (from MPI_Comm_rank)
!        MYTID  = mpi ID of this process (from MPI_Comm_rank)
      INTEGER MAXNOD, MAXCL2
      PARAMETER (MAXNOD = 200, MAXCL2 = 10000)
      REAL*8 TOTWAL
!CSK      INTEGER IPRPAR, NTASK, NCODE, NDEGDI, MASTER, MYNUM, MYTID, NODTOT
      INTEGER IPRPAR, NTASK, NCODE, NDEGDI, MYNUM, MYTID, NODTOT
      INTEGER NUMNOD, ICOUNT,MTIDS(0:MAXNOD), MPARID,NODE,NODES,IPAREX
      LOGICAL PARHER, SLAVE, TIMING, DEBUG, PARIO
      CHARACTER*20   NODNAM(0:MAXNOD), MYNAME
!CSK      COMMON /DIRACINFPAR/ TOTWAL                                        &
!CSK     &       ,IPRPAR, NTASK, NCODE, NDEGDI, MASTER, MYNUM, MYTID, NODTOT &
      COMMON /DIRACINFPAR/ TOTWAL                                        &
     &       ,IPRPAR, NTASK, NCODE, NDEGDI, MYNUM, MYTID, NODTOT         &
     &       ,ICOUNT, MTIDS, MPARID, NODE, NODES, IPAREX                 &
     &       ,PARHER, SLAVE, TIMING, DEBUG, PARIO                        &
     &       ,NODNAM, MYNAME

      EQUIVALENCE (NODTOT, NUMNOD)
!C     ... for some obscure reason NODTOT in Dalton is called NUMNOD in Dirac.
!C         This equivalence is important such that both e.g. abacus routines using NODTOT
!C         and dirac routines using NUMNOD are working correctly together /hjaaj Dec 07

#if defined (INT_STAR8)
!     integer array ISTAT contains MPI_SOURCE information.
!     Proper use of ISTAT on 64-bit machines in 
!     combination with INT_STAR8 requires explicit declaration 
!     as INTEGER*4 /March-2007 sk - This is only true for 
!     32-bit MPI libraries.
!     64-bit MPI libraries require INTEGER ISTAT (which will then be 
!     INTEGER*8)
!
!     radovan: it is not ideal that istat type is not defined ifndef INT_STAR8
!              i do not want to fix it, because maybe there is a reason
      INTEGER   ISTAT
#endif

!C -- end of infpar.h --
#endif
