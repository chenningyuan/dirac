#ifdef __CVERSION__
#define MAXNOD 100
struct common_infpar {
    double totval;
    int icount, kprint, mtids[MAXNOD], mparid, node, nodes;
    int ntask, ncode, parher, timing;
} infpar_;
#else /* __CVERSION__ */
      PARAMETER (MAXNOD=100)
      REAL*8 TOTWAL
      LOGICAL PARHER, TIMING
      COMMON /INFPAR/ TOTWAL, ICOUNT,
     &                KPRINT, MTIDS(MAXNOD), MPARID, NODE, NODES,
     &                NTASK, NCODE, PARHER, TIMING
#endif
