*
*   Routine for multiplying matrix (or vector) with scalar
*
      subroutine matfac(A,NROWS,NCOLUMNS,IFAC)
      implicit real*8 (A-H,O-Z)
      dimension A(NROWS,NCOLUMNS)
*
*  multiply
*
      do J=1,NROWS,1
         do I=1,NCOLUMNS,1
            A(I,J) = IFAC * A(I,J)
         end do
      end do
      return
      end
