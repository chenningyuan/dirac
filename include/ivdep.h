#if defined (SYS_CRAY) || defined (SYS_T3D)
CDIR$ IVDEP
#endif
#if defined (SYS_ALLIANT)
CVD$ NODEPCHK
#endif
#if defined (SYS_CONVEX)
C$DIR NO_RECURRENCE
#endif
#if defined (SYS_IBM)
Cvectorization note: ignore vector dependence
C*VDIR: IGNORE RECRDEPS
#endif
#if defined (SYS_SX)
*VDIR NODEP
#endif
#if !defined (SYS_CRAY) && !defined (SYS_ALLIANT) && !defined (SYS_IBM) && !defined (SYS_CONVEX) && !defined (SYS_T3D) && !defined (SYS_SX)
Cvectorization note: ignore vector dependence
#endif
