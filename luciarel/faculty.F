*  Calculate faculty of integer
*
      subroutine faculty(FACUL,N)
      implicit none
      integer FACUL,N
      integer I
      FACUL=1
      do 10 I=1,N,1
         FACUL = I*FACUL
10    continue
      return
      end
