*  Unfold lower triangle square matrix to full size;
*  I suppose that the input matrix is a one-dim. field with
*  N(N/2 + 1) elements
*  1) write vector to matrix
*  2) unfold considering symmetry (e.g. hermitian or not)
*
      subroutine unfold(ARRAY,VECTOR,IDIM,ISYM,NTRIAEL,IFAC)
      implicit real*8 (A-H,O-Z)
      dimension ARRAY(IDIM,IDIM),VECTOR(NTRIAEL)
*
*  unfold
*
      ICOUNT = 0
      do N = 1,IDIM,1
         do M = 1,N,1
            ARRAY(N,M) = IFAC * VECTOR(ICOUNT+M)
         end do
         ICOUNT = ICOUNT + N
      end do
*
*  multiply with sign
*
      do I=1,IDIM,1
         do J=I+1,IDIM,1
            ARRAY(I,J) = ISYM * ARRAY(J,I)
         end do
      end do
      return
      end
