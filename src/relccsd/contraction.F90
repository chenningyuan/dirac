  module contraction
  use wrapper_sorting
    implicit none

    public contraction_222
    public contraction_442   
    public contraction_422

     private

     contains

    subroutine contraction_222 (varray, dim_v, tarray, dim_t, warray, vtrans, ttrans, wtrans, &
           &        ext_v, ext_t, index_sum, dia_factor, dia_sum_factor)


! this routine will be modified soon....


!----------calling variables-----------------------------------
    real(8), intent(in) :: varray(*), tarray(*)      
    character(*), intent(in) :: vtrans, ttrans, wtrans
    complex(8), intent(in) :: dia_factor, dia_sum_factor               !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    integer, intent(in) :: dim_v, dim_t, ext_v(nrep), ext_t(nrep), index_sum(nrep)
    real(8), intent(out) :: warray(*)
!--------------------------------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"


    if (vtrans.eq.'C') then
    if (carith) call conjuga (dim_v,varray,1) 
    endif
    
    if (ttrans.eq.'C') then
    if (carith) call conjuga (dim_t,tarray,1) 
    endif
    call CNTRCT (vtrans,ttrans,ext_v,ext_t,index_sum,dia_factor,varray,tarray,dia_sum_factor,warray,NREP)

    end subroutine


   subroutine contraction_442 (left_order, right_order, im_order, rightTensor, imTensor, &
           &                   dia_factor, dia_sum_factor)
     
!----------calling variables-----------------------------------
    real(8), intent(in) :: dia_factor, dia_sum_factor   !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    character (*), dimension(4), intent(in) :: left_order, right_order
    character (*), dimension(2), intent(in) :: im_order
    real(8), intent(in)   :: rightTensor(:)   
    integer,allocatable   :: x_row(:),x_column(:)
    integer, allocatable  :: y_column(:)              
    real(8), intent(inout):: imTensor(*)
    character (1) :: transA
    real(8)  :: v_sort(10000)
    real(8)  :: buf4(100000)    
!--------------------------------------------------------------
#include "param.inc"
#include "symm.inc"

      allocate(y_column(nrep))
      allocate(x_column(nrep))
      allocate(x_row(nrep))

      y_column = no

      v_sort = 0.0d0
! 1) fetch the 2-e integrals and sort them for BLAS. 
      call v_for_blas (13,left_order,right_order,im_order, transA, v_sort, x_row, x_column)
! 2) the call of SRT22 routine will not be there. it will be replaced by general sorting routine for fixed arrays. 
!     likewise, BUF4 will also be replaced.   
      CALL SRT22 (NREP,MULTB,NVVT,NO,NO,NVVOT,LOVVOT,LLOVVT,rightTensor,BUF4)
! 3) contraction routine for symmetry-packed arrays.
      call CNTRCT (transA,'N',x_row,y_column,x_column,dia_factor*A1,v_sort, BUF4,dia_sum_factor*A1,imTensor,NREP) 
  
     deallocate(y_column)
     deallocate(x_column)
     deallocate(x_row)

    end subroutine

   subroutine contraction_422 (left_order, right_order, im_order, rightTensor, imTensor, &
           &                   dia_factor, dia_sum_factor, nrep)
     
!----------calling variables-----------------------------------
    character(*), dimension(4), intent(in)   :: left_order, right_order
    character(*), dimension(2), intent(in)   :: im_order
!    real(8), intent(in), optional :: left_array(:), right_array(:)   !when the left or right array is a fixed array then we will provide input.    
    real(8), intent(in) :: rightTensor(*)   
    integer, intent(in) :: nrep
    real(8), intent(in)        :: dia_factor, dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)     :: imTensor(*)
!--------------------------------------------------------------
    integer  :: x_row(nrep), x_column(nrep),i
    character (1)  :: trans
    real(8)  :: v_sort(10000)
#include "param.inc"
!  1) fetch the fixed array, properly sorted to contract with BLAS.

!  2) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS.

      v_sort = 0.0d0

      call v_for_blas (22,left_order,right_order,im_order, trans, v_sort, x_row, x_column)

!  3) determine the row and column index for the matrix-vector multiplication routine. 
!     I will like to fetch a flexible array   

      CALL XGEMV (trans,x_row,x_column,dia_factor*A1,v_sort,x_row,rightTensor,1,dia_sum_factor*A1,imTensor,1) 

    end subroutine


  end module
