#if defined (SYS_ALLIANT)
CVD$ NOVECTOR
#endif
#if defined (SYS_CONVEX)
C$DIR SCALAR
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D)
CDIR$ NEXTSCALAR
#endif
#if defined (SYS_SX)
*VDIR NOVECTOR
#endif
