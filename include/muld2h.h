C
C$Id: muld2h.h 981113 hjaaj
C
#if defined (SYS_VAX) || defined (SYS_ALLIANT) || defined (SYS_IBM) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_HPUX) || defined (SYS_SUN) || defined (SYS_LINUX) || defined (SYS_DARWIN)
      MULD2H(I,J) = IEOR(I-1,J-1) + 1
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D)
      MULD2H(I,J) = XOR(I-1,J-1) + 1
#endif
#if !defined (SYS_VAX) && !defined (SYS_ALLIANT) && !defined (SYS_IBM) && !defined (SYS_CONVEX) && !defined (SYS_CRAY) && !defined (SYS_AIX) && !defined (SYS_PARAGON) && !defined (SYS_DEC) && !defined (SYS_IRIX) && !defined (SYS_HPUX) && !defined (SYS_SUN) && !defined (SYS_T3D) && !defined (SYS_LINUX) && !defined (SYS_DARWIN)
      You must define muld2h in muld2h.h file for this computer.
#endif
