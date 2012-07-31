#if defined (SYS_ALLIANT)
!VD$ NOVECTOR
#endif
#if defined (SYS_CONVEX)
!$DIR SCALAR
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D)
!DIR$ NEXTSCALAR
#endif
#if defined (SYS_SX)
!VDIR NOVECTOR
#endif
