C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : gencon.h
C
C   - Parameter for maximum number of functions in a block:
C     MXCON_PARAM = NPRIM * l(l+1)/2 * NCONTR
C   - DIRCON true if contracted basis set is used
C
      PARAMETER (MXCON_PARAM = 300)
C     
      LOGICAL DIRCON(2)
      COMMON /GENCON/ DIRCON
