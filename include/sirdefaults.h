#if !defined (SYS_VAX) && !defined (SYS_CRAY) && !defined (SYS_ALLIANT) && !defined (SYS_CONVEX) && !defined (SYS_AIX) && !defined (SYS_PARAGON) && !defined (SYS_DEC) && !defined (SYS_IRIX) && !defined (SYS_T3D)
      MWORK  = 300 000
#endif
#if defined (SYS_VAX)
      MWORK  = 96 000
C     ------   96 000 = 1 500 VAX memory pages (each 512 bytes)
#endif
#if defined (SYS_CRAY) || defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_T3D)
      MWORK  = 2 000 000
#endif
C     ------ size of working set (real memory, not virtual memory)
