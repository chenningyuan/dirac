#if defined (SYS_CRAY) || defined (SYS_T3D)
!DIR$ IVDEP
#endif
#if defined (SYS_ALLIANT)
!VD$ NODEPCHK
#endif
#if defined (SYS_CONVEX)
!$DIR NO_RECURRENCE
#endif
#if defined (SYS_IBM)
!vectorization note: ignore vector dependence
!*VDIR: IGNORE RECRDEPS
#endif
#if defined (SYS_SX)
!VDIR NODEP
#endif
#if !defined (SYS_CRAY) && !defined (SYS_ALLIANT) && !defined (SYS_IBM) && !defined (SYS_CONVEX) && !defined (SYS_T3D) && !defined (SYS_SX)
!vectorization note: ignore vector dependence
#endif
