#if defined (SYS_CRAY) || defined (SYS_T3D)
C CRAY: do not vectorize this subroutine
CDIR$ NOVECTOR
#endif
#if defined (SYS_SX)
*VDIR NOVECTOR
#endif
