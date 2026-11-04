!*** Unpack four-indices ***
      integer ibit2,IBIT2M,IBIT4, IBIT4M
#if defined (INT_STAR8)
      PARAMETER (IBIT2=65536,IBIT2M=65535,IBIT4=4294967296,             &
     & IBIT4M=4294967295)
#else
      PARAMETER (IBIT2 = 256,IBIT2M = 255,IBIT4 = 65536,IBIT4M = 65535)
#endif
#if defined (INT_STAR8)
      INDGET(I,J) = IAND(ISHFT(I,-J),IBIT2M)
      LGET(I) = INDGET(I,0)
      KGET(I) = INDGET(I,16)
      JGET(I) = INDGET(I,32)
      IGET(I) = INDGET(I,48)
      IPACK(I,J,K,L) = IOR(IOR(IOR(ISHFT(L,0),                          &
     &ISHFT(K,16)),ISHFT(J,32)),ISHFT(I,48))
#else
      INDGET(I,J) = IAND(ISHFT(I,-J),IBIT2M)
      LGET(I) = INDGET(I,0)
      KGET(I) = INDGET(I,8)
      JGET(I) = INDGET(I,16)
      IGET(I) = INDGET(I,24)
      IPACK(I,J,K,L) = IOR(IOR(IOR(ISHFT(L,0),                          &
     &ISHFT(K,8)),ISHFT(J,16)),ISHFT(I,24))
#endif
