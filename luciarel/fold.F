*  Fold square matrix to lower tringle form
*  The output matrix is a one-dim. field with
*  N(N+1)/2 elements
*
      subroutine fold(VECTOR,ARRAY,IDIM,ISYM,NTRIAEL)
      implicit real*8 (A-H,O-Z)
      dimension ARRAY(IDIM,IDIM),VECTOR(NTRIAEL)
*
      NTEST = 0
      do M=1,IDIM,1
         J = 0
         do L=1,M,1
            J = J + L
         end do
         ICOUNT = J
         do N=M,IDIM,1
            VECTOR(ICOUNT) = ARRAY(N,M)
            if (N.eq.IDIM.and.M.eq.IDIM) goto 100
            ICOUNT = ICOUNT + N
         end do
      end do
100   continue
      if (ICOUNT.eq.NTRIAEL) then
         write(*,*) 'I have folded successfully.'
      end if
      if (NTEST.ne.0) then
        write(*,*) 'Vector is:'
        do I=1,NTRIAEL,1
           write(*,*) VECTOR(I)
        end do
      end if
      return
      end
