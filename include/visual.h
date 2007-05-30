C     visualization/integration common block
C     ======================================
C
C     LBIT(LBIT_VISUAL,1)
C     LBIT(LBIT_VISUAL,2) 2D
C     LBIT(LBIT_VISUAL,3) 3D
C     LBIT(LBIT_VISUAL,4)
C     LBIT(LBIT_VISUAL,5) skip +- (SKIPEP)
C     LBIT(LBIT_VISUAL,6) skip ++ (SKIPEE)
C     LBIT(LBIT_VISUAL,7) use orbital string
C     LBIT(LBIT_VISUAL,8) integrate
C
C     IQ_N           =  1 !    number density
C     IQ_M           =  2 !    magnetization
C     IQ_DIVM        =  3 !div(magnetization)
C     IQ_ROTM        =  4 !rot(magnetization)
C     IQ_MDIP        =  5 !    magnetic dipole density
C     IQ_J           =  6 !    current density
C     IQ_DIVJ        =  7 !div(current density)
C     IQ_ROTJ        =  8 !rot(current density)
C     IQ_JPARA       =  9 !    paramagnetic current density
C     IQ_DIVJPARA    = 10 !div(paramagnetic current density)
C     IQ_ROTJPARA    = 11 !rot(paramagnetic current density)
C     IQ_JDIA        = 12 !     diamagnetic current density
C     IQ_VPARA       = 13 !    paramagnetic velocity density
C     IQ_DIVVPARA    = 14 !div(paramagnetic velocity density)
C     IQ_ROTVPARA    = 15 !rot(paramagnetic velocity density) a.k.a. vorticity
C     IQ_ELF         = 16 !    ELF
C     IQ_K           = 17 !    kinetic energy density
C     IQ_H           = 18 !    helicity density
C     IQ_GAMMA5      = 19 !    \gamma_5 density
C     IQ_NUC_SPIN_IX = 20 !    give me a name ...
C     IQ_NUC_SPIN_JX = 21 !    give me a name ...
C     IQ_NUC_SPIN_IY = 22 !    give me a name ...
C     IQ_NUC_SPIN_JY = 23 !    give me a name ...
C     IQ_NUC_SPIN_IZ = 24 !    give me a name ...
C     IQ_NUC_SPIN_JZ = 25 !    give me a name ...
C     IQ_EDIP        = 26 !    electric dipole density
C
      INTEGER MAX_NR_XVECTORS
      PARAMETER (MAX_NR_XVECTORS = 25)
C
      CHARACTER       XVECTOR_FILE(MAX_NR_XVECTORS)*6,                  &
     &                ORBITAL_STRING(2)*8,                              &
     &                ORIENTATION(3)*2
      LOGICAL         DO_VISUAL
C
      COMMON /CBVISC/ XVECTOR_FILE,                                     &
     &                ORBITAL_STRING,                                   &
     &                ORIENTATION
      COMMON /CBVISL/ DO_VISUAL
      COMMON /CBVISI/ LBIT_VISUAL,                                      &
     &                NCUBE(0:3),                                       &
     &                NR_XVECTORS,                                      &
     &                IQ_XVECTOR(MAX_NR_XVECTORS),                      &
     &                IREP_XVECTOR(MAX_NR_XVECTORS),                    &
     &                IWHICH_XVECTOR(MAX_NR_XVECTORS),                  &
     &                IOPTION_XVECTOR(MAX_NR_XVECTORS)
      COMMON /CBVISR/ ADD_3D, COMMON_FACTOR
