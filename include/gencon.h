C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : gencon.h
C SOURCE  : $Source$
C
C   - Parameter for maximum number of functions in a block:
C     MXCON_PARAM = NPRIM * l(l+1)/2 * NCONTR
C   - LARGE/SMALLCON logicals tell if basis set is general 
C     contracted or not.
C   - ILARGEC/ISMALLC give number of primitives (this data
C     gets lost because we have repetition of same exponents
C     in more then 1 block.
C   - DIRCON true if contracted basis set is used
C   - MOLFCON true if molfdir basis set is used
C
      PARAMETER (MXCON_PARAM = 300)
C     
      COMMON /GENCON/ LARGECON(MXCENT), SMALLCON(MXCENT),
     &                DIRCON,MOLFCON
      COMMON /GCDEFN/ ILARGEC(MXCENT,MXQN),ISMALLC(MXCENT,MXQN)
C
      LOGICAL LARGECON, SMALLCON, DIRCON, MOLFCON
