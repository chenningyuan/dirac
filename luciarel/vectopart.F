*  Unfold lower triangle square matrix to full size;
*  I suppose that the input matrix is a one-dim. field with
*  N(N/2 + 1) elements
*  1) write vector to matrix
*  2) unfold considering symmetry (e.g. hermitian or not)
*  output matrix id partition of a fully dimensioned matrix!
*
      subroutine vectopart(ARRAY,VECTOR,ITOTDIM,ISYM,NTRIAEL,
     &                     IFAC,IPARTDIM,IROWOFF,ICOLOFF)
      implicit real*8 (A-H,O-Z)
      dimension ARRAY(ITOTDIM,ITOTDIM),VECTOR(NTRIAEL)
*
      NTEST = 000
*
*  copy to array
*
      ICOUNT = 0
      do N = 1,IPARTDIM,1
         do M = 1,N,1
            ARRAY(N+IROWOFF-1,M+ICOLOFF-1)
     &            = IFAC * VECTOR(ICOUNT+M)
         end do
         ICOUNT = ICOUNT + N
      end do
*
*  multiply with sign and generate second entry (unfold)
*
      do I=1,IPARTDIM,1
         do J=I+1,IPARTDIM,1
            ARRAY(I+IROWOFF-1,J+ICOLOFF-1)
     &            = ISYM * ARRAY(J+IROWOFF-1,I+ICOLOFF-1)
         end do
      end do
*
      if (NTEST.ge.100) then
        write(6,*) 'Info from vectopart:'
        write(6,*) 'Dimension of partition is: ',IPARTDIM
        write(6,*) 'Row offset is ',IROWOFF
        write(6,*) 'Column offset is ',ICOLOFF
        call wrtmat(ARRAY,ITOTDIM,ITOTDIM,ITOTDIM,ITOTDIM)
      end if
*
      return
      end
