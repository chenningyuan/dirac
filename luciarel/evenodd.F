*  This (stupid) little subroutine merely checks whether an
*  integer is odd or even and passes this info to the calling
*  routine.
*
      subroutine evenodd(IEVOD,INTTEST)
      implicit real*8 (A-H,O-Z)
* add 1
      INTTESTP = INTTEST + 1
* put to real value
      A=INTTEST
      B=INTTESTP
* divide both values by 2
      FIRST = A/2
      SECOND = B/2
* cut off places right of point
      ICUT1 = FIRST
      ICUT2 = SECOND
* compare results
      if (ICUT1.eq.ICUT2) then
* INTTEST is even
         IEVOD = 2
      else
* INTTEST is odd
         IEVOD = 1
      end if
      return
      end
