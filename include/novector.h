C next loop may be run in scalar or parallel mode,
C but not in vector mode (vector mode inefficient)
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
