*
*   Routine(1) for error message
*
      subroutine hammer(LAST,IMAX)
      implicit none
      integer LAST,IMAX
      write(*,*) 'The dimensions have exceeded the maximum work space!'
      write(*,*) 'Dimensions are  ',LAST
      write(*,*) 'Maximum space for this part is ',IMAX
      write(*,*) 'Picasso will put down the brush.'
      write(*,*) 'Bye-bye.'
      Call Abend()
      end
