*
*   Routine for (complex) transposing rectangular matrix A
*
      subroutine adjoint(B,A,NROWS,NCOLUMNS,IMRECO)
      implicit real*8 (A-H,O-Z)
      dimension A(NROWS,NCOLUMNS),B(NROWS,NCOLUMNS)
*
*   Write zeroes to matrix B
*
      do K=1,NCOLUMNS,1
         do L=1,NROWS,1
            B(L,K) = 0.
         end do
      end do
*
*  transpose
*
      do I=1,NCOLUMNS,1
         do J=1,NROWS,1
            B(I,J) = IMRECO * A(J,I)
         end do
      end do
      return
      end
