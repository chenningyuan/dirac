  module contraction
  use wrapper_sorting
    implicit none

    public contraction_222
    public contraction_442   
    public contraction_242
    public contraction_422
    public contraction_444
    public contraction_424
    public contraction_244

     private

     contains

   subroutine contraction_222 (left_order,right_order,im_order,leftTensor,rightTensor,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep)

!----------calling variables-----------------------------------
    character(2), dimension(2), intent(in)   ::left_order,im_order,right_order
    real(8), intent(in)                      :: rightTensor(*),leftTensor(*)   
    integer, intent(in)                      :: nrep
    real(8), intent(in)                      :: dia_factor, dia_sum_factor             !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)                   :: imTensor(*)
!--------------------------------------------------------------
    integer        :: x_row(nrep),x_column(nrep),y_row(nrep),y_column(nrep)
    character(1)   :: transA,transB
    character(2)   :: freeindex_left(2)
    character(2)   :: freeindex_right(2)
    logical        :: array_equal = .true. 
    real*8,allocatable :: ImTensor_temp(:) 
    real*8,allocatable :: ImTensor_temp_temp(:) 

#include "param.inc"

    transA = 'N'  
    transB = 'N'  

!  1) fetch the fixed arrays, properly sorted to contract with BLAS.
   
   call onebody_for_blas (left_order,im_order,right_order,transA,x_row,x_column,freeindex_left)

   call onebody_for_blas (right_order,left_order,im_order,transB,y_row,y_column,freeindex_right)
   
   array_equal = all((/freeindex_left(1),freeindex_right(2)/)==im_order)

   if (array_equal) then

   call cntrct (transA,transB,x_row,y_column,x_column,dia_factor*A1,leftTensor,rightTensor,dia_sum_factor*A1,imTensor,NREP) 

   else

   allocate (ImTensor_temp(dot_product(x_row,y_column)))

   call cntrct (transA,transB,x_row,y_column,x_column,dia_factor*A1,leftTensor,rightTensor,A0,imTensor_temp,NREP) 

   allocate (ImTensor_temp_temp(dot_product(x_row,y_column)))

   call srt1c1 (nrep,x_row,y_column,imTensor_temp,imTensor_temp_temp)

   deallocate (ImTensor_temp)

   call xaxpy (dot_product(x_row,y_column),dia_sum_factor*A1,imTensor_temp_temp,1,imTensor,1) 

   deallocate (ImTensor_temp_temp)

   endif

    end subroutine

   subroutine contraction_442 (left_order,right_order,im_order,ImTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor,RightTensor)
     
!----------calling variables-----------------------------------
    real(8), intent(in) :: dia_factor, dia_sum_factor   !dia_factor = \alpha in DGEMM ; dia_sum_factor= \beta in DGEMM
    character(2),dimension(4), intent(in) :: left_order, right_order
    character(2),dimension(2), intent(in) :: im_order
    real(8),optional,intent(in)           :: RightTensor(:),LeftTensor(:)   
    integer,intent(in)    :: nrep
    integer,allocatable   :: x_row(:),x_column(:)
    integer,allocatable   :: y_column(:),y_row(:)              
    real(8),intent(inout) :: ImTensor(*)
    integer               :: match_index  
    character(1)          :: transA,transB
    real(8)               :: v_sort(10000)
    real(8)               :: right_sorted(10000)    

    character(2)          :: sorted_oparray_left(4)
    character(2)          :: sorted_oparray_right(4)
    logical               :: array_equal = .true. 
    real(8),allocatable   :: ImTensor_temp(:) 
    real(8),allocatable   :: ImTensor_temp_temp(:) 
!--------------------------------------------------------------
!     We need to eliminate the direct use of the kvvoot (and other) symmetry arrays here to be able to remove symm.inc (should be a function call that give the dimension arrays)

#include "param.inc"

      allocate(y_column(nrep))
      allocate(x_column(nrep))
      allocate(x_row(nrep))
      allocate(y_row(nrep))
 
      transA = 'N'
      transB = 'N'
      v_sort = 0.0d0
      right_sorted = 0.0d0

! 1) fetch 2-e integrals and sort them for BLAS.

    call decide_sorting_routine(13,left_order,right_order,im_order,transA,match_index,sorted_oparray_left)


     if (present(LeftTensor)) then 
      call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column,LeftTensor)
     else
      call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column)
     endif 

! 2) general sorting routine for fixed arrays. 

     call decide_sorting_routine(31,right_order,im_order,sorted_oparray_left,transB,match_index,sorted_oparray_right)


     array_equal = all((/sorted_oparray_left(1),sorted_oparray_right(4)/)==im_order)
  

     if (present(RightTensor)) then 
     call sorting_for_blas(match_index,.false.,transB,right_sorted,y_row,y_column,rightTensor)
     else
     call sorting_for_blas(match_index,.false.,transB,right_sorted,y_row,y_column)
     endif  

! 3) contraction routine for symmetry-packed arrays.

   if (array_equal) then
      call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,v_sort,right_sorted,dia_sum_factor*A1,imTensor,NREP) 

    else

   allocate (ImTensor_temp(dot_product(x_row,y_column)))
     call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,v_sort,right_sorted,A0,imTensor_temp,NREP) 

   allocate (ImTensor_temp_temp(dot_product(x_row,y_column)))

   call srt1c1 (nrep,x_row,y_column,imTensor_temp,imTensor_temp_temp)

   deallocate (ImTensor_temp)

   call xaxpy (dot_product(x_row,y_column),dia_sum_factor*A1,imTensor_temp_temp,1,imTensor,1) 

   deallocate (ImTensor_temp_temp)

   endif

     deallocate(y_column)
     deallocate(x_column)
     deallocate(x_row)
     deallocate(y_row)

    end subroutine

   subroutine contraction_422 (left_order,right_order,im_order,rightTensor,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor)
     
!----------calling variables-----------------------------------
    character(2),dimension(4),intent(in)   :: left_order
    character(2),dimension(2),intent(in)   :: im_order, right_order
    real(8),intent(in),optional :: LeftTensor(:)            !when the left array is a fixed array then we will provide input.    
    real(8),intent(in)     :: rightTensor(*)   
    integer,intent(in)     :: nrep
    real(8),intent(in)     :: dia_factor, dia_sum_factor  !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8),intent(inout)  :: imTensor(*)
!--------------------------------------------------------------
    integer        :: x_row(nrep), x_column(nrep)
    character(1)   :: transA
    integer        :: match_index 
    real(8)        :: v_sort(10000)
#include "param.inc"

!  1) fetch the array which is flexible in terms of indices. sort it to make suitable for BLAS.
!     also determine the row and column index for the matrix-vector multiplication routine.

     v_sort = 0.0d0

     transA = 'N'

     call decide_sorting_routine(22,left_order,right_order,im_order,transA,match_index)

     if (present(LeftTensor)) then

     call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column,presorted=LeftTensor,if_totsym=.true.)

     else

     call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column,if_totsym=.true.)

     endif


!  2) carry out the contraction. 

      call xgemv (transA,x_row,x_column,dia_factor*A1,v_sort,x_row,rightTensor,1,dia_sum_factor*A1,imTensor,1) 

    end subroutine

   subroutine contraction_242 (left_order,right_order,im_order,LeftTensor,ImTensor, &
           &                   dia_factor,dia_sum_factor,nrep,RightTensor)

!-------------------------------------------------------------
!    use of xgemv is different from dgemm because when I mention 'T', then the dimensionality of the incoming array is also changed.
!-------------------------------------------------------------

!----------calling variables-----------------------------------
    character(2), dimension(2), intent(in)   :: left_order,im_order
    character(2), dimension(4), intent(in)   :: right_order
    real(8), intent(in)                      :: LeftTensor(:)
    real(8), intent(in), optional            :: RightTensor(:)   
    integer, intent(in)                      :: nrep
    real(8), intent(in)                      :: dia_factor,dia_sum_factor             !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)                   :: ImTensor(*)
!--------------------------------------------------------------
    integer        :: y_row(nrep), y_column(nrep)
    character(1)   :: transA,transB
    integer        :: match_index,i 
    real(8)        :: right_sort(10000)
!---------------------------------------------------------------    
#include "param.inc"
#include "complex.inc"

     right_sort = 0.0d0
     transA = 'N'
!  1> fetch the integral array, properly sorted to contract with BLAS.

   call decide_sorting_routine(22,right_order,im_order,left_order,transA,match_index)

   if (present(RightTensor)) then
   call sorting_for_blas(match_index,.false.,transA,right_sort,y_row,y_column,RightTensor,if_totsym=.true.)
   else
   call sorting_for_blas(match_index,.false.,transA,right_sort,y_row,y_column,if_totsym=.true.)
   endif


!  2>  to do B = XA  type operation with xgemv, we are doing B^T = A^T*X^T then transpose(B^T). that implies following modifications.   

  if (transA == 'N') then
   transB = 'C'
  else
   transB = 'N'
   if (carith) call conjuga(y_row(1)*y_column(1),right_sort,1)
  endif

  if (transB == 'N') then
  if (carith) call conjuga(y_column(1),LeftTensor,1)
  else
  if (carith) call conjuga(y_row(1),LeftTensor,1)
  endif

   call xgemv (transB,y_row(1),y_column(1),dia_factor*A1,right_sort,y_row(1),LeftTensor,1,dia_sum_factor*A1,ImTensor,1) 

  if (transB == 'N') then
  if (carith) call conjuga(y_row(1),ImTensor,1)
  else
  if (carith) call conjuga(y_column(1),ImTensor,1)
  endif

    end subroutine

   subroutine contraction_444 (left_order,right_order,im_order,rightTensor,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor)
    
!----------calling variables-----------------------------------
    character(*), dimension(4), intent(in)   :: left_order, right_order, im_order
    real(8), intent(in), optional :: LeftTensor(:)   !when the left array is a fixed array then we will provide input.    
    real(8), intent(in) :: RightTensor(:)   
    integer, intent(in) :: nrep
    real(8), intent(in) :: dia_factor, dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)  :: imTensor(*)
!--------------------------------------------------------------
   integer        :: x_row(nrep),x_column(nrep),y_row(nrep),y_column(nrep)
   character(1)   :: transA, transB
   integer        :: save_index,i
   character(2)   :: sorted_opArray_left(4),sorted_opArray_right(4)
   real(8)        :: v_sort(10000)
   real(8)        :: right_sorted(10000)
   real(8),allocatable :: im_temp(:)
   real(8),allocatable :: im_temp_temp(:)

#include "param.inc"

   v_sort       = 0.0d0
   right_sorted = 0.0d0
   transA = 'N'
   transB = 'N'

!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

    call decide_sorting_routine(22,left_order,right_order,im_order,transA,save_index,sorted_opArray_left)

    if (present(LeftTensor)) then 
    call sorting_for_blas(save_index,.false.,transA,v_sort,x_row,x_column, presorted=LeftTensor,if_TotSym=.false.)
    else
    call sorting_for_blas(save_index,.false.,transA,v_sort,x_row,x_column,if_TotSym=.false.)
    endif

!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.

    call decide_sorting_routine(22,right_order,im_order,sorted_opArray_left,transB,save_index,sorted_opArray_right)

    call sorting_for_blas(save_index,.false.,transB,right_sorted,y_row,y_column,rightTensor,if_TotSym=.false.)

!  3) matrix-matrix multiplication routine runs for each irreps.

   if (dia_sum_factor==0.0d0) then

!  if (transB=='n') then
   allocate(im_temp(dot_product(x_row,y_column)))
!  else
!  allocate(im_temp(dot_product(x_row,y_row)))
!  endif

   im_temp = 0.0d0

   call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,v_sort,right_sorted,A0,im_temp,NREP) 

   v_sort = 0.0d0

   transA = 'N'

   call decide_sorting_routine(22,im_order,left_order,(/sorted_opArray_left(1:2),sorted_opArray_right(3:4)/), &
                       & transA,save_index,inverse=.true.)

   if (transA=='C') then

!  if (transB=='n') then
   allocate(im_temp_temp(dot_product(x_row,y_column)))
!  else
!  allocate(im_temp_temp(dot_product(x_row,y_row)))
!  endif

   im_temp_temp= 0.0d0

   call SRT1C1 (nrep,x_row,y_column,im_temp,im_temp_temp)

   call sorting_for_blas(save_index,.true.,transA,v_sort,x_row,x_column,im_temp_temp,if_TotSym=.false.)

     deallocate(im_temp_temp)
   else

   call sorting_for_blas(save_index,.true.,transA,v_sort,x_row,x_column,im_temp)

   endif

   deallocate(im_temp)
   call xaxpy (dot_product(x_row,x_column),A1,v_sort,1,imTensor,1)

   else
    call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,v_sort,right_sorted,dia_sum_factor*A1,imTensor,NREP) 
   endif

  end subroutine

   subroutine contraction_424 (left_order,right_order,im_order,ImTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor,RightTensor)

!----------calling variables-----------------------------------
    character(2), dimension(4), intent(in)   :: left_order,im_order
    character(2), dimension(2), intent(in)   :: right_order
    real(8), intent(in), optional :: LeftTensor(:), RightTensor(:) 
    integer, intent(in) :: nrep
    real(8), intent(in) :: dia_factor, dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout) :: ImTensor(*)
    integer :: x_row(nrep), x_column(nrep), y_row(nrep), y_column(nrep)
    character(1) :: transA, transB
    integer :: match_index,i
    character(2) :: FreeIndx(2),sorted_opArray(3)
    real(8) :: v_sort (10000)
    real(8),allocatable :: im_temp(:)
    real(8),allocatable :: im_temp_temp(:)

!   real(8) :: im_temp(10000)
!   real(8) :: im_temp_temp(10000)

#include "param.inc"
   v_sort = 0.0d0
   transA = 'N'
   transB = 'N'
!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

   call decide_sorting_routine(31,left_order,right_order,im_order,transA,match_index,sorted_opArray)

 if (present(LeftTensor)) then
 
   call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column,LeftTensor)
 
 else

   call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column)

 endif
!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.

   call onebody_for_blas (right_order,left_order,im_order,transB,y_row,y_column,FreeIndx)

!  3) matrix-matrix multiplication routine runs for each irreps.

   allocate(im_temp(dot_product(x_row,y_column)))
   im_temp = 0.0d0

  if (present(RightTensor)) call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,v_sort,RightTensor,A0,im_temp,NREP) 

   v_sort = 0.0d0
!  4) re-sort the arrays to align with the desired structure of the imtermediate.

   transA = 'N'

   call decide_sorting_routine(31,im_order,left_order,(/sorted_opArray,FreeIndx(2)/),transA,match_index,inverse=.true.)

   if (transA=='C') then

   allocate(im_temp_temp(dot_product(x_row,y_column)))
   im_temp_temp = 0.0d0
   call SRT1C1 (nrep,x_row,y_column,im_temp,im_temp_temp)

  deallocate(im_temp)
   call sorting_for_blas(match_index,.true.,transA,v_sort,x_row,x_column,im_temp_temp)

  deallocate(im_temp_temp)
   else

   call sorting_for_blas(match_index,.true.,transA,v_sort,x_row,x_column,im_temp)

  deallocate(im_temp)
   endif

   call xaxpy (dot_product(x_row,x_column),A1,v_sort,1,ImTensor,1)

  end subroutine

  subroutine contraction_244 (left_order,right_order,im_order,ImTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor,RightTensor)

!----------calling variables-----------------------------------
    character(2),dimension(4),intent(in) :: right_order,im_order
    character(2),dimension(2),intent(in) :: left_order
    real(8),intent(in),optional :: LeftTensor(:), RightTensor(:) 
    integer,intent(in) :: nrep
    real(8),intent(in) :: dia_factor, dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8),intent(inout) :: imTensor(*)
    integer :: x_row(nrep), x_column(nrep), y_row(nrep), y_column(nrep)
    character(1) :: transA, transB
    integer :: match_index,i
    character(2) :: FreeIndx(2),sorted_opArray(4)
    real(8) :: right_sort (10000)
    real(8),allocatable :: im_temp(:)
    real(8),allocatable :: im_temp_temp(:)

!   real(8) :: im_temp(10000)
!   real(8) :: im_temp_temp(10000)

#include "param.inc"
   right_sort = 0.0d0
   transA = 'N'
   transB = 'N'
!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

   call decide_sorting_routine(13,right_order,im_order,left_order,transB,match_index,sorted_opArray)

 if (present(RightTensor)) then
 
   call sorting_for_blas(match_index,.false.,transB,right_sort,y_row,y_column,RightTensor)
 
 else

   call sorting_for_blas(match_index,.false.,transB,right_sort,y_row,y_column)

 endif

!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.

   call onebody_for_blas (left_order,im_order,sorted_opArray,transA,x_row,x_column,FreeIndx)

!  3) matrix-matrix multiplication routine runs for each irreps.

!  if (transB=='n') then
   allocate(im_temp(dot_product(x_row,y_column)))
!  else
!  allocate(im_temp(dot_product(x_column,y_column)))
!  endif
   im_temp = 0.0d0

  if (present(LeftTensor)) call CNTRCT (transA,transB,x_row,y_column,x_column,dia_factor*A1,LeftTensor,right_sort,A0,im_temp,NREP) 

   right_sort = 0.0d0
!  4) re-sort the arrays to align with the desired structure of the imtermediate.

   transA = 'N'

   call decide_sorting_routine(13,im_order,left_order,(/FreeIndx(1),sorted_opArray(2:4)/),transA,match_index,inverse=.true.)

   if (transA=='C') then

!  if (transB=='n') then
  allocate(im_temp_temp(dot_product(x_row,y_column)))
!  else
!  allocate(im_temp_temp(dot_product(x_column,y_column)))
!  endif
   im_temp_temp = 0.0d0
   call SRT1C1 (nrep,x_row,y_column,im_temp,im_temp_temp)

 deallocate(im_temp)
   call sorting_for_blas(match_index,.true.,transA,right_sort,x_row,x_column,im_temp_temp)

 deallocate(im_temp_temp)
   else
   call sorting_for_blas(match_index,.true.,transA,right_sort,x_row,x_column,im_temp)
 deallocate(im_temp)
   endif
   call xaxpy (dot_product(x_row,x_column),A1,right_sort,1,ImTensor,1)

  end subroutine

  end module
