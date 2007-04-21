C     visualization/integration common block
C     ======================================
C
C     LBIT(LBIT_VISUAL,1)            1D (not implemented)
C     LBIT(LBIT_VISUAL,2)            2D
C     LBIT(LBIT_VISUAL,3)            3D
C     LBIT(LBIT_VISUAL,4)            read solution vector
C     LBIT(LBIT_VISUAL,5)            skip +- (SKIPEP)
C     LBIT(LBIT_VISUAL,6)            skip ++ (SKIPEE)
C     LBIT(LBIT_VISUAL,7)            use orbital string
C     LBIT(LBIT_VISUAL,8)            integrate
C
C     LBIT(LBIT_VISUAL_QUANTITY,1)   density
C     LBIT(LBIT_VISUAL_QUANTITY,2)   magnetization
C     LBIT(LBIT_VISUAL_QUANTITY,3)   current density
C     LBIT(LBIT_VISUAL_QUANTITY,4)   ELF
C     LBIT(LBIT_VISUAL_QUANTITY,5)   r \times j
C     LBIT(LBIT_VISUAL_QUANTITY,6)   div(magnetization)
C     LBIT(LBIT_VISUAL_QUANTITY,7)   div(current density)
C     LBIT(LBIT_VISUAL_QUANTITY,8)   rot(magnetization)
C     LBIT(LBIT_VISUAL_QUANTITY,9)   rot(current density)
C     LBIT(LBIT_VISUAL_QUANTITY,20)  diamagnetic current density
C     LBIT(LBIT_VISUAL_QUANTITY,21)  kinetic energy density
C     LBIT(LBIT_VISUAL_QUANTITY,22)  helicity density
C     LBIT(LBIT_VISUAL_QUANTITY,23)  \gamma_5 density
C     LBIT(LBIT_VISUAL_QUANTITY,24)  paramagnetic current density
C     LBIT(LBIT_VISUAL_QUANTITY,25)  div(paramagnetic current density)
C     LBIT(LBIT_VISUAL_QUANTITY,26)  rot(paramagnetic current density)
C
      INTEGER NR_QUANTITIES
      PARAMETER (NR_QUANTITIES = 26)
C
      CHARACTER*6     XVECFILE
      CHARACTER*8     ORBSTRING1, ORBSTRING2
      CHARACTER       QUANTITY_NAME(NR_QUANTITIES)*30
      CHARACTER       ORIENTATION(3)*2
      LOGICAL         DO_VISUAL, QUANTITY_IS_A_VECTOR(NR_QUANTITIES)
C
      COMMON /CBVISC/ XVECFILE, ORBSTRING1, ORBSTRING2, QUANTITY_NAME,  &
     &                ORIENTATION
      COMMON /CBVISL/ DO_VISUAL, QUANTITY_IS_A_VECTOR
      COMMON /CBVISI/ NCUBE(0:3), NXVEC, IREP_VISUAL,                   &
     &                LBIT_VISUAL, LBIT_VISUAL_QUANTITY
      COMMON /CBVISR/ ADD_3D
