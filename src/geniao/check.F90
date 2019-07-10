module check

contains

 subroutine check_normalization(matrix,row,col)
       complex(8), parameter   :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
       integer, intent(in)     :: row,col
       complex(8), intent(in)  :: matrix(row,col)
       integer                 :: ierr,i,j
       complex(8)              :: sum_column

       ierr=0
       do i=1,col
         sum_column = ZERO
         do j=1,row
           sum_column = sum_column + conjg(matrix(j,i))*matrix(j,i)
         end do
         if (REAL(sum_column) .lt. 0.999D0 .or. AIMAG(sum_column) .gt. 1.0D-3) then
            ierr=1
         endif
       end do

       if (ierr.eq.1) then
        write(*,*)"**Warning** The matrix is not normalized."
       else if (ierr.eq.0) then
        write(*,*)"The matrix is normalized."
       endif
 end subroutine check_normalization

 subroutine check_SymOrth(C,S,row,col)
       complex(8), parameter   :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
       integer, intent(in)     :: row,col
       complex(8), intent(in)  :: C(row,col), S(row,row)
       complex(8)              :: check(col,col),Id(col,col)
       integer                 :: i,j

       Id = ZERO
       forall (i=1:col) Id(i,i)=ONE
       check = matmul(conjg(transpose(C)),matmul(S,C))
       write(*,*)"Deviation from orthogonality      :", &
                 norm2([norm2(real(check-Id)),norm2(aimag(check-Id))])
 end subroutine check_SymOrth

 subroutine deviation_nelec(rdm,n,nelec)
       complex(8), parameter   :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
       integer,    intent(in)  :: nelec,n
       complex(8), intent(in)  :: rdm(n,n)
       complex(8)              :: sum_diag
     
       sum_diag = ZERO
       do i=1,n
         sum_diag = sum_diag + rdm(i,i)
       enddo
       write(*,*)"Deviation in number of electrons  :", real(sum_diag - nelec)
 end subroutine deviation_nelec
end module
