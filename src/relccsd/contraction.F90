  module contraction
  use wrapper_sorting
    implicit none

    public contraction_222
    public contraction_442   
    public contraction_422
    public contraction_444
    public contraction_424
    public contraction_244


     private

     contains

   subroutine contraction_222 (left_order,right_order,im_order,leftTensor,rightTensor,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep)

!----------calling variables-----------------------------------
    character(2), dimension(4), intent(in)   :: left_order
    character(2), dimension(2), intent(in)   :: im_order, right_order
!    real(8), intent(in), optional :: left_array(:), right_array(:)   !when the left or right array is a fixed array then we will provide input.    
    real(8), intent(in) :: rightTensor(*),leftTensor(*)   
    integer, intent(in) :: nrep
    real(8), intent(in) :: dia_factor, dia_sum_factor             !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)  :: imTensor(*)
!--------------------------------------------------------------
    integer        :: x_row(nrep),x_column(nrep),y_row(nrep),y_column(nrep)
    character(1)   :: transA,transB
#include "param.inc"

    transA = 'N'  
    transB = 'N'  

!  1) fetch the fixed arrays, properly sorted to contract with BLAS.
   
   call onebody_for_blas (left_order,im_order,right_order,transA,x_row,x_column)

   call onebody_for_blas (right_order,left_order,im_order,transB,y_row,y_column)

  ! call CNTRCT (vtrans,ttrans,ext_v,ext_t,index_sum,dia_factor,varray,tarray,dia_sum_factor,warray,NREP)

   call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,leftTensor,rightTensor,dia_sum_factor*A1,imTensor,NREP) 
    end subroutine

   subroutine contraction_442 (left_order, right_order, im_order, rightTensor, imTensor, &
           &                   dia_factor, dia_sum_factor, nrep)
     
!----------calling variables-----------------------------------
    real(8), intent(in) :: dia_factor, dia_sum_factor   !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    character (*), dimension(4), intent(in) :: left_order, right_order
    character (*), dimension(2), intent(in) :: im_order
    real(8), intent(in)   :: rightTensor(:)   
    integer, intent(in)   :: nrep
    integer, allocatable  :: x_row(:),x_column(:)
    integer, allocatable  :: y_column(:),y_row(:)              
    real(8), intent(inout):: imTensor(*)
    integer               :: match_index  
    character(1)          :: transA,transB
    real(8)  :: v_sort(10000)
    real(8)  :: right_sorted(100000)    
!--------------------------------------------------------------
!     We need to eliminate the direct use of the kvvoot (and other) symmetry arrays here to be able to remove symm.inc (should be a function call that give the dimension arrays)

#include "param.inc"

      allocate(y_column(nrep))
      allocate(x_column(nrep))
      allocate(x_row(nrep))
      allocate(y_row(nrep))

      v_sort = 0.0d0
! 1) fetch the 2-e integrals and sort them for BLAS. 

     call decide_sorting_routine(13,left_order,right_order,im_order,transA,match_index)

     call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column)

!     call v_for_blas (13,left_order,right_order,im_order, transA, v_sort, x_row, x_column)

! 2) general sorting routine for fixed arrays. 

!     call decide_sorting_routine(31,right_order,left_order,im_order,transB,match_index)

!     call sorting_for_blas(match_index,.false.,transB,right_sorted,y_row,y_column,rightTensor)

!      there is a bit of problem in using the same subroutines because in align_array routine the call sequence is somewhat different. well, I can use a different call sequence for the arguments. need a bit of rechecking of the dimensionbs of the character arrays. ideally, they should be generalized for any number of characters. 

     call T_for_blas(31,right_order,left_order,im_order,transB,rightTensor,right_sorted,y_row,y_column)

! 3) contraction routine for symmetry-packed arrays.

      call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,v_sort,right_sorted,dia_sum_factor*A1,imTensor,NREP) 
  
     deallocate(y_column)
     deallocate(x_column)
     deallocate(x_row)
     deallocate(y_row)

    end subroutine

   subroutine contraction_422 (left_order, right_order, im_order, rightTensor, imTensor, &
           &                   dia_factor, dia_sum_factor, nrep)
     
!----------calling variables-----------------------------------
    character(2), dimension(4), intent(in)   :: left_order
    character(2), dimension(2), intent(in)   :: im_order, right_order
!    real(8), intent(in), optional :: left_array(:), right_array(:)   !when the left or right array is a fixed array then we will provide input.    
    real(8), intent(in) :: rightTensor(*)   
    integer, intent(in) :: nrep
    real(8), intent(in)     :: dia_factor, dia_sum_factor             !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)  :: imTensor(*)
!--------------------------------------------------------------
    integer        :: x_row(nrep), x_column(nrep)
    character(1)   :: transA
    integer        :: match_index 
    real(8)        :: v_sort(10000)
#include "param.inc"
!  1) fetch the fixed array, properly sorted to contract with BLAS.

!  2) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS.

     v_sort = 0.0d0

     transA = 'N'
     call decide_sorting_routine(22,left_order,right_order,im_order,transA,match_index)

     call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column)
      
!    call v_for_blas (22,left_order,right_order,im_order, trans, v_sort, x_row, x_column)

!  3) determine the row and column index for the matrix-vector multiplication routine. 
!     I will like to fetch a flexible array   

      CALL XGEMV (transA,x_row,x_column,dia_factor*A1,v_sort,x_row,rightTensor,1,dia_sum_factor*A1,imTensor,1) 

    end subroutine

   subroutine contraction_444 (left_order,right_order,im_order,rightTensor,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep)
    
!----------calling variables-----------------------------------
    character(*), dimension(4), intent(in)   :: left_order, right_order, im_order
!   real(8), intent(in), optional :: left_array(:), right_array(:)   !when the left or right array is a fixed array then we will provide input.    
    real(8), intent(in) :: rightTensor(:)   
    integer, intent(in) :: nrep
    real(8), intent(in) :: dia_factor, dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)     :: imTensor(*)
!--------------------------------------------------------------
   integer        :: x_row(nrep), x_column(nrep), y_row(nrep), y_column(nrep)
   character(1)   :: transA, transB
   integer        :: save_index
   real(8)        :: v_sort(10000)
   real(8)        :: right_sorted(10000)
#include "param.inc"

   v_sort = 0.0d0
   transA = 'N'
   transB = 'N'

!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

     call decide_sorting_routine(22,left_order,right_order,im_order,transA,save_index)

     call sorting_for_blas(save_index,.false.,transA,v_sort,x_row,x_column)

!     call v_for_blas (22,left_order,right_order,im_order, transA, v_sort, x_row, x_column)

!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.
 
     call T_for_blas(22,right_order,left_order,im_order,transB,rightTensor,right_sorted,y_row,y_column)

!  3) matrix-matrix multiplication routine runs for each irreps.

    call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,v_sort,right_sorted,dia_sum_factor*A1,imTensor,NREP) 
!    call CNTRCT (transA,transB,NOOT,NOOT,NVVT,dia_factor*A1,v_sort,right_sorted,dia_sum_factor*A1,imTensor,NREP)     
   end subroutine

   subroutine contraction_424 (left_order,right_order,im_order,rightTensor,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep,nv1)

! careful about nv1. soon replace it by common input index. 

!----------calling variables-----------------------------------
    character(*), dimension(4), intent(in)   :: left_order, right_order, im_order
!   real(8), intent(in), optional :: left_array(:), right_array(:)   !when the left or right array is a fixed array then we will provide input.    
    real(8), intent(in) :: rightTensor(*)   
    integer, intent(in) :: nrep,nv1
    real(8), intent(in) :: dia_factor, dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout) :: imTensor(*)
!--------------------------------------------------------------
   integer :: x_row(nrep), x_column(nrep), y_row(nrep), y_column(nrep)
   character(1) :: transA, transB
   integer :: match_index
   character(2) :: FreeIndx,sorted_opArray(3)
   real(8) :: v_sort(10000)
   real(8) :: im_temp(10000)
#include "param.inc"

   v_sort =  0.0d0
   transA = 'N'
   transB = 'N'
!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

  call decide_sorting_routine(31,left_order,right_order,im_order,transA,match_index,sorted_opArray)

  call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column)

!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.

   call onebody_for_blas (right_order,left_order,im_order,transB,y_row,y_column,FreeIndx)

!  3) matrix-matrix multiplication routine runs for each irreps.

   call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,v_sort,rightTensor,A0,im_temp,NREP) 

   call decide_sorting_routine(31,im_order,left_order,(/sorted_opArray,FreeIndx/),transA,match_index,inverse=.true.)

   call sorting_for_blas(match_index,.true.,transA,v_sort,x_row,x_column,im_temp)

   call xaxpy (nv1,A1,v_sort,1,imTensor,1)

  end subroutine

  end module
