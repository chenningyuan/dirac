C
C$Id$
C
#if defined (SYS_VAX) || defined (SYS_ALLIANT) || defined (SYS_IBM) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_HPUX) || defined (SYS_SUN) || defined (SYS_LINUX) || defined (SYS_SX) || defined (SYS_DARWIN)
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = IEOR(I,J)
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = SHIFTL(I,J)
      IBTSHR(I,J) = SHIFTR(I,J)
      IBTXOR(I,J) = XOR(I,J)
#endif
#if !defined (SYS_VAX) && !defined (SYS_ALLIANT) && !defined (SYS_IBM) && !defined (SYS_CONVEX) && !defined (SYS_CRAY) && !defined (SYS_AIX) && !defined (SYS_PARAGON) && !defined (SYS_DEC) && !defined (SYS_IRIX) && !defined (SYS_HPUX) && !defined (SYS_SUN) && !defined (SYS_T3D) && !defined (SYS_LINUX) && !defined (SYS_SX) && !defined (SYS_DARWIN)
      You must define IBTFUN in comdeck file for this computer.
#endif
