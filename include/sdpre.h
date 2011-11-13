#if defined (SYS_CRAY) || defined (SYS_T3D)
      SDPRE(I) = FLOAT(I)
#else
      SDPRE(I) = DFLOAT(I)
#endif
