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

     external cntrct
     external srt1c1
     external srt1c1n

  !!   integer,parameter :: maxdim = 10000000

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
#include "complex.inc"

    transA = 'N'  
    transB = 'N'  

!  1) fetch the fixed arrays, properly sorted to contract with BLAS.
   
   call onebody_for_blas (left_order,im_order,right_order,transA,x_row,x_column,freeindex_left)

   call onebody_for_blas (right_order,left_order,im_order,transB,y_row,y_column,freeindex_right)
   
   array_equal = all((/freeindex_left(1),freeindex_right(2)/)==im_order)

   if (array_equal) then

   call cntrct (transA,transB,x_row,y_column,x_column,dia_factor*A1,leftTensor,rightTensor,dia_sum_factor*A1,imTensor,NREP) 

   else

   allocate (ImTensor_temp(dot_product(x_row,y_column)*rcw))

   ImTensor_temp = 0.0d0

   call cntrct (transA,transB,x_row,y_column,x_column,dia_factor*A1,leftTensor,rightTensor,A0,imTensor_temp,NREP) 

   allocate (ImTensor_temp_temp(dot_product(x_row,y_column)*rcw))

   call srt1c1n (nrep,x_row,y_column,imTensor_temp,imTensor_temp_temp)

!   deallocate (ImTensor_temp)

   call xaxpy (dot_product(x_row,y_column),dia_sum_factor*A1,imTensor_temp_temp,1,imTensor,1) 

   deallocate (ImTensor_temp)
   deallocate (ImTensor_temp_temp)

   endif

    end subroutine

   subroutine contraction_442 (left_order,right_order,im_order,ImTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor,RightTensor)
     
!----------calling variables-----------------------------------
    real(8), intent(in) ::  dia_sum_factor   !dia_factor = \alpha in DGEMM ; dia_sum_factor= \beta in DGEMM
    real(8), intent(in) ::  dia_factor        
    character(2),dimension(4), intent(in) :: left_order, right_order
    character(2),dimension(2), intent(in) :: im_order
    real(8),optional,intent(in)           :: RightTensor(:),LeftTensor(:)   
    integer,intent(in)    :: nrep
    integer   :: x_row(nrep),x_column(nrep)
    integer   :: y_column(nrep),y_row(nrep)              
    real(8),intent(inout) :: ImTensor(*)
    integer               :: match_index  
    character(1)          :: transA,transB
    real(8),allocatable   :: v_sort(:)
    real(8),allocatable   :: right_sorted(:)    
    character(2)          :: sorted_oparray_left(4)
    character(2)          :: sorted_oparray_right(4)
    logical               :: array_equal = .true. 
    integer               :: asym_factor_R,asym_factor_L,i 
    integer               :: maxdim1, maxdim2 
    real(8),allocatable   :: ImTensor_temp(:) 
    real(8),allocatable   :: ImTensor_temp_temp(:) 
!--------------------------------------------------------------

#include "param.inc"
#include "complex.inc"



      maxdim1 = sorted_dimension(left_order)
      maxdim2 = sorted_dimension(right_order)

  allocate(v_sort(maxdim1*rcw))  
  allocate(right_sorted(maxdim2*rcw))  

 
      transA = 'N'
      transB = 'N'
!!      v_sort = 0.0d0
!!      right_sorted = 0.0d0

! 1) fetch 2-e integrals and sort them for BLAS.

    call decide_sorting_routine(13,left_order,right_order,im_order,transA,match_index,asym_factor_L, &
    &   .true.,sorted_opArray=sorted_oparray_left)

!    write(*,*)'left array',(sorted_opArray_left(i),i=1,4),match_index,transA
     

     if (present(LeftTensor)) then 
      call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column,LeftTensor)
     else
      call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column)
     endif 

! 2) general sorting routine for fixed arrays. 

     call decide_sorting_routine(31,right_order,im_order,sorted_oparray_left,transB,match_index,asym_factor_R, &
    &         .false.,sorted_opArray=sorted_oparray_right)

!    write(*,*)'right array',(sorted_opArray_right(i),i=1,4)

     array_equal = all((/sorted_oparray_left(1),sorted_oparray_right(4)/)==im_order)

     if (present(RightTensor)) then 
     call sorting_for_blas(match_index,.false.,transB,right_sorted,y_row,y_column,rightTensor)
     else
     call sorting_for_blas(match_index,.false.,transB,right_sorted,y_row,y_column)
     endif  

! 3) contraction routine for symmetry-packed arrays.

   if (array_equal) then
      call cntrct (transA,transB,x_row,y_column,x_column,asym_factor_R*asym_factor_L*dia_factor*A1,v_sort, &
   &                right_sorted,dia_sum_factor*A1,imTensor,nrep) 

    else

   allocate (ImTensor_temp(dot_product(x_row,y_column)*rcw))
   ImTensor_temp = 0.0d0
     call cntrct (transA,transB,x_row,y_column,x_column,asym_factor_R*asym_factor_L*dia_factor*A1,v_sort, &
    &              right_sorted,A0,imTensor_temp,nrep) 

   allocate (ImTensor_temp_temp(dot_product(x_row,y_column)*rcw))

   call srt1c1n (nrep,x_row,y_column,imTensor_temp,imTensor_temp_temp)

   deallocate (ImTensor_temp)

   call xaxpy (dot_product(x_row,y_column),dia_sum_factor*A1,imTensor_temp_temp,1,imTensor,1) 

   deallocate (ImTensor_temp_temp)

   endif

  deallocate(v_sort)  
  deallocate(right_sorted)  

    end subroutine

   subroutine contraction_422 (left_order,right_order,im_order,rightTensor,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor)
     
!----------calling variables-----------------------------------
    character(2),dimension(4),intent(in)   :: left_order
    character(2),dimension(2),intent(in)   :: im_order, right_order
    real(8),intent(in),optional :: LeftTensor(:)            !when the left array is a fixed array then we will provide input.    
    real(8),intent(in)     :: rightTensor(*)   
    integer,intent(in)     :: nrep
    real(8),intent(in)     :: dia_sum_factor  !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8),intent(in)     :: dia_factor  !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8),intent(inout)  :: imTensor(*)
!--------------------------------------------------------------
    integer        :: x_row(nrep), x_column(nrep),y_row(nrep),y_column(nrep)
    character(1)   :: transA,transB
    integer        :: match_index,i 
    integer        :: asym_factor
    integer        :: maxdim, lda
    logical        :: array_equal
    real(8),allocatable        :: v_sort(:)
    real(8),allocatable        :: RightTensor_sort(:)
    character(2)   :: sorted_left(4),freeindx(2)
#include "param.inc"
#include "complex.inc"

!---------------known problem----------------------------------------------
! call to an XGEMV routine is different from DGEMM. we need to provide the untransposed dimensions even for the cases when
!transposition is necessary. way out is to use DGEMM in this case also. however, as of now, transA='C' case has not
!appeared.
!--------------------------------------------------------------------------
!  1> fetch the array which is flexible in terms of indices. sort it to make suitable for BLAS.
!     also determine the row and column index for the matrix-vector multiplication routine.

  maxdim = sorted_dimension(left_order)

  allocate(v_sort(maxdim*rcw))  

!!     v_sort = 0.0d0

     transA = 'N'

     call decide_sorting_routine(22,left_order,right_order,im_order,transA,match_index,asym_factor,.true., &
    &                        sorted_opArray=sorted_left)


     if (present(LeftTensor)) then

     call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column,presorted=LeftTensor,if_totsym=.true.)

     else

     call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column,if_totsym=.true.)

     endif

     call onebody_for_blas (right_order,sorted_left(3:4),im_order,transB,y_row,y_column) 

     array_equal = all((/sorted_left(1),sorted_left(2)/)==im_order)

     write(*,*)'array_equal, always T eXpected',array_equal,transA

         if (transA.eq.'n'.or.transA.eq.'N') then
            lda = x_row(1)
         else
            lda = x_column(1)
         endif

     if (transB=='T') then

      allocate(rightTensor_sort(x_column(1)*rcw))

      call srt1c1n (nrep,y_column,y_row,rightTensor,rightTensor_sort)

!!      if (carith) call conjuga(x_column(1),rightTensor,1)

      call xgemm (transA,'N',x_row(1),1,x_column(1),asym_factor*dia_factor*A1,v_sort, &
      &       lda,rightTensor_sort,x_column(1),dia_sum_factor*A1,imTensor,x_row(1))

!     call xgemv (transA,x_row(1),x_column(1),asym_factor*dia_factor*A1,v_sort,x_row(1),rightTensor_sort, &
!     &            1,dia_sum_factor*A1,imTensor,1) 

      deallocate(rightTensor_sort)
     
      else

      call xgemm (transA,'N',x_row(1),1,x_column(1),asym_factor*dia_factor*A1,v_sort, &
      &       lda,rightTensor,x_column(1),dia_sum_factor*A1,imTensor,x_row(1))

!      call xgemv (transA,x_row(1),x_column(1),asym_factor*dia_factor*A1,v_sort,x_row(1),rightTensor,1,dia_sum_factor*A1,imTensor,1) 
      endif

      deallocate(v_sort)
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
    real(8), intent(in)                      :: dia_sum_factor        !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(in)                      :: dia_factor            !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)                   :: ImTensor(*)
!--------------------------------------------------------------
    character(2)   :: freeindex_left(2)
    character(2)   :: sorted_oparray_right(4) 
    integer        :: y_row(nrep),y_column(nrep)
    integer        :: x_row(nrep),x_column(nrep)
    character(1)   :: transA,transB,trans1A
    integer        :: match_index,i 
    integer        :: asym_factor
    integer        :: maxdim
    integer        :: ldb 
    logical        :: array_equal
    real(8),allocatable        :: right_sort(:)
    real(8),allocatable        :: left_sort(:)
    real(8),allocatable        :: im_temp(:)
!---------------------------------------------------------------    
#include "param.inc"
#include "complex.inc"

    maxdim = sorted_dimension(right_order)

   allocate(right_sort(maxdim*rcw))  

!!     right_sort = 0.0d0
     transA = 'N'
     trans1A = 'N'

!**  1> fetch the integral array, properly sorted to contract with BLAS.

   call decide_sorting_routine(22,right_order,im_order,left_order,transA,match_index,asym_factor,.true.,&
   &                          sorted_opArray=sorted_oparray_right)


   if (present(RightTensor)) then
   call sorting_for_blas(match_index,.false.,transA,right_sort,y_row,y_column,presorted=RightTensor,if_totsym=.true.)
   else
   call sorting_for_blas(match_index,.false.,transA,right_sort,y_row,y_column,if_totsym=.true.)
   endif

!**  2> sort onebody integral properly to use in BLAS.

   call onebody_for_blas (left_order,im_order,sorted_oparray_right,trans1A,x_row,x_column)

   allocate (left_sort(y_row(1)*rcw))


   if (trans1A=='T') then


!!   call xcopy (y_row(1),LeftTensor,1,left_sort,1)  
   

!!   if (carith) call conjuga(y_row(1),left_sort,1)

    call srt1c1n (nrep,x_column,x_row,LeftTensor,left_sort)

   else
    
   call xcopy (y_row(1),LeftTensor,1,left_sort,1)  

   endif

     array_equal = all((/sorted_oparray_right(3),sorted_oparray_right(4)/)==im_order)
     
     write(*,*)'array_equal, always T eXpected',array_equal,transA

         if (transA.eq.'N'.or.transA.eq.'n') then
            ldb = y_row(1)
         else
            ldb = y_column(1)
         endif

   allocate(im_temp(y_column(1)*rcw))

   im_temp =0.0d0
!**  3>  to do B = XA  type operation with xgemv, we are doing B^T = A^T*X^T then transpose(B^T). that implies following modifications.   
   call xgemm ('N',transA,1,y_column(1),y_row(1),asym_factor*dia_factor*A1,left_sort, &
      &       1,right_sort,ldb,dia_sum_factor*A0,im_temp,1)


!!!   if (carith) call conjuga(y_column(1),im_temp,1)

   call xaxpy (y_column(1),A1,im_temp,1,imTensor,1)


!!if (transA == 'N') then
!! transB = 'C'
!!else
!! transB = 'N'
!! if (carith) call conjuga(y_row(1)*y_column(1),right_sort,1)
!!endif

!!if (transB == 'N') then
!!if (carith) call conjuga(y_column(1),left_sort,1)
!!else
!!if (carith) call conjuga(y_row(1),left_sort,1)
!!endif

!! call xgemv (transB,y_row(1),y_column(1),asym_factor*dia_factor*A1,right_sort,y_row(1),left_sort,1,dia_sum_factor*A0,im_temp,1) 

!!if (transB == 'N') then
!!if (carith) call conjuga(y_row(1),im_temp,1)
!! call xaxpy (y_row(1),A1,im_temp,1,imTensor,1)
!!else
!!if (carith) call conjuga(y_column(1),im_temp,1)
!! call xaxpy (y_column(1),A1,im_temp,1,imTensor,1)
!!endif

   deallocate(im_temp)

   deallocate (left_sort)
   deallocate (right_sort)

    end subroutine

   subroutine contraction_444 (left_order,right_order,im_order,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep,RightTensor,LeftTensor)
    
!----------calling variables-----------------------------------
    character(*), dimension(4), intent(in)   :: left_order, right_order, im_order
    real(8), intent(in), optional :: LeftTensor(:)   !when the left array is a fixed array then we will provide input.    
    real(8), intent(in), optional :: RightTensor(:)   
    integer, intent(in) :: nrep
    real(8), intent(in) :: dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(in) :: dia_factor   
    real(8), intent(inout)  :: imTensor(*)
!--------------------------------------------------------------
   integer        :: x_row(nrep),x_column(nrep),y_row(nrep),y_column(nrep)
   character(1)   :: transA, transB
   integer        :: save_index,i
   character(2)   :: sorted_opArray_left(4),sorted_opArray_right(4)
   real(8),allocatable  :: v_sort(:)
   real(8),allocatable  :: right_sorted(:)
   real(8),allocatable :: im_temp(:)
   real(8),allocatable :: im_temp_temp(:)
   integer             :: asym_factor_R,asym_factor_L
   integer             :: maxdim1, maxdim2

#include "param.inc"
#include "complex.inc"


  maxdim1 = max(sorted_dimension(left_order),sorted_dimension(im_order))
  maxdim2 = sorted_dimension(right_order)


  allocate(v_sort(maxdim1*rcw))  
  allocate(right_sorted(maxdim2*rcw))  

!!   v_sort       = 0.0d0
!!   right_sorted = 0.0d0
   transA = 'N'
   transB = 'N'

!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

    call decide_sorting_routine(22,left_order,right_order,im_order,transA,save_index,asym_factor_L,.true., &
    &  sorted_opArray=sorted_opArray_left)

    if (present(LeftTensor)) then 
    call sorting_for_blas(save_index,.false.,transA,v_sort,x_row,x_column, presorted=LeftTensor,if_TotSym=.false.)
    else
    call sorting_for_blas(save_index,.false.,transA,v_sort,x_row,x_column,if_TotSym=.false.)
    endif

!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.


    call decide_sorting_routine(22,right_order,im_order,sorted_opArray_left,transB,save_index,asym_factor_R,&
    &          .false.,sorted_opArray=sorted_opArray_right)

    if (present(RightTensor)) then 
    call sorting_for_blas(save_index,.false.,transB,right_sorted,y_row,y_column, presorted=RightTensor,if_TotSym=.false.)
    else
    call sorting_for_blas(save_index,.false.,transB,right_sorted,y_row,y_column,rightTensor,if_TotSym=.false.)
    endif

!  3) matrix-matrix multiplication routine runs for each irreps.

   if (dia_sum_factor==0.0d0) then

   allocate(im_temp(dot_product(x_row,y_column)*rcw))

   im_temp = 0.0d0

   call cntrct (transA,transB,x_row,y_column,x_column,asym_factor_L*asym_factor_R*dia_factor*A1,v_sort,right_sorted,A0,im_temp,NREP) 

 !!  v_sort = 0.0d0

   transA = 'N'

   call decide_sorting_routine(22,im_order,left_order,(/sorted_opArray_left(1:2),sorted_opArray_right(3:4)/), &
                       & transA,save_index,asym_factor_R,.false.,sorted_opArray=sorted_opArray_right,inverse=.true.)

   if (transA=='T') then

   allocate(im_temp_temp(dot_product(x_row,y_column)*rcw))

   call srt1c1n (nrep,x_row,y_column,im_temp,im_temp_temp)

   call sorting_for_blas(save_index,.true.,transA,v_sort,x_row,x_column,im_temp_temp,if_TotSym=.false.)

     deallocate(im_temp_temp)
   else

   call sorting_for_blas(save_index,.true.,transA,v_sort,x_row,x_column,im_temp,if_TotSym=.false.)

   endif

   deallocate(im_temp)
   call xaxpy (dot_product(x_row,x_column),A1*asym_factor_R,v_sort,1,imTensor,1)

   else
    call cntrct (transA,transB,x_row,y_column,x_column,asym_factor_R*asym_factor_L*dia_factor*A1,v_sort,right_sorted, &
   &              dia_sum_factor*A1,imTensor,nrep) 
   endif

   deallocate(v_sort)
   deallocate(right_sorted)

  end subroutine

   subroutine contraction_424 (left_order,right_order,im_order,ImTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor,RightTensor)

!----------calling variables-----------------------------------
    character(2), dimension(4), intent(in)   :: left_order,im_order
    character(2), dimension(2), intent(in)   :: right_order
    real(8), intent(in), optional :: LeftTensor(:), RightTensor(:) 
    integer, intent(in) :: nrep
    real(8), intent(in) :: dia_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(in) :: dia_sum_factor  !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout) :: ImTensor(*)
    integer :: x_row(nrep), x_column(nrep), y_row(nrep), y_column(nrep)
    character(1) :: transA, transB
    integer :: match_index,i
    character(2) :: FreeIndx(2),sorted_opArray(4)
    real(8),allocatable :: v_sort (:)
    integer :: asym_factor 
    integer :: maxdim 
    real(8),allocatable :: im_temp(:)
    real(8),allocatable :: im_temp_temp(:)

#include "param.inc"
#include "complex.inc"


  maxdim = max (sorted_dimension(left_order),sorted_dimension(im_order))

  allocate(v_sort(maxdim*rcw))  

!!   v_sort = 0.0d0
   transA = 'N'
   transB = 'N'
!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

   call decide_sorting_routine(31,left_order,right_order,im_order,transA,match_index,asym_factor, &
   &             .true.,sorted_opArray=sorted_opArray)

!  write(*,*)'right array',(sorted_opArray(i),i=1,4)


 if (present(LeftTensor)) then

   call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column,presorted=LeftTensor)
 
 else

   call sorting_for_blas(match_index,.false.,transA,v_sort,x_row,x_column)

 endif

!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.

   call onebody_for_blas (right_order,sorted_opArray,im_order,transB,y_row,y_column,FreeIndx)

!  3) matrix-matrix multiplication routine runs for each irreps.

   allocate(im_temp(dot_product(x_row,y_column)*rcw))
   im_temp = 0.0d0

  if (present(RightTensor)) call cntrct (transA,transB,x_row,y_column,x_column,asym_factor*dia_factor*A1,v_sort, &
   &                                      RightTensor,A0,im_temp,NREP) 

!!   v_sort = 0.0d0
!  4) re-sort the arrays to align with the desired structure of the intermediate.

   transA = 'N'

   call decide_sorting_routine(31,im_order,left_order,(/sorted_opArray(1:3),FreeIndx(2)/),transA,match_index, &
           &                    asym_factor,.false.,inverse=.true.,sorted_opArray=sorted_opArray)

!  write(*,*)'right array',(sorted_opArray(i),i=1,4)
   if (transA=='T') then

   allocate(im_temp_temp(dot_product(x_row,y_column)*rcw))

   call srt1c1n (nrep,x_row,y_column,im_temp,im_temp_temp)

!!      if (carith) call conjuga(dot_product(x_row,y_column),im_temp_temp,1)

   deallocate(im_temp)

   call sorting_for_blas(match_index,.true.,transA,v_sort,x_row,x_column,im_temp_temp,if_TotSym=.false.)

   deallocate(im_temp_temp)

   else

   call sorting_for_blas(match_index,.true.,transA,v_sort,x_row,x_column,im_temp,if_TotSym=.false.) 

   deallocate(im_temp)

   endif

   call xaxpy (dot_product(x_row,x_column),A1*asym_factor*dia_sum_factor,v_sort,1,ImTensor,1)

   deallocate(v_sort)

  end subroutine

  subroutine contraction_244 (left_order,right_order,im_order,ImTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor,RightTensor)

!----------calling variables-----------------------------------
    character(2),dimension(4),intent(in) :: right_order,im_order
    character(2),dimension(2),intent(in) :: left_order
    real(8),intent(in),optional :: LeftTensor(:), RightTensor(:) 
    integer,intent(in) :: nrep
    real(8),intent(in) ::  dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8),intent(in) :: dia_factor           !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8),intent(inout) :: imTensor(*)
    integer :: x_row(nrep), x_column(nrep), y_row(nrep), y_column(nrep)
    character(1) :: transA, transB
    integer :: match_index,i
    character(2) :: FreeIndx(2),sorted_opArray(4)
    real(8),allocatable :: right_sort (:)
    integer :: asym_factor
    integer :: maxdim 
    real(8),allocatable :: im_temp(:)
    real(8),allocatable :: im_temp_temp(:)

#include "param.inc"
#include "complex.inc"

  maxdim = max (sorted_dimension(right_order),sorted_dimension(im_order))

  allocate(right_sort(maxdim*rcw))  

!!   right_sort = 0.0d0
   transA = 'N'
   transB = 'N'
!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

   call decide_sorting_routine(13,right_order,im_order,left_order,transB,match_index,asym_factor, &
    &                    .true.,sorted_opArray=sorted_opArray)

 if (present(RightTensor)) then
 
   call sorting_for_blas(match_index,.false.,transB,right_sort,y_row,y_column,RightTensor)
 
 else

   call sorting_for_blas(match_index,.false.,transB,right_sort,y_row,y_column)

 endif

!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.

   call onebody_for_blas (left_order,im_order,sorted_opArray,transA,x_row,x_column,FreeIndx)

!  3) matrix-matrix multiplication routine runs for each irreps.

   allocate(im_temp(dot_product(x_row,y_column)*rcw))
   im_temp = 0.0d0

  if (present(LeftTensor)) call cntrct (transA,transB,x_row,y_column,y_row,asym_factor*dia_factor*A1,LeftTensor,&
  &        right_sort,A0,im_temp,nrep) 

!!   right_sort = 0.0d0
!  4) re-sort the arrays to align with the desired structure of the imtermediate.

   transA = 'N'

!!   write(*,*)'write sorted before final',match_index,(/FreeIndx(1),sorted_opArray(2:4)/)
   call decide_sorting_routine(13,im_order,left_order,(/FreeIndx(1),sorted_opArray(2:4)/),transA,match_index, &
   &                        asym_factor,.false.,inverse=.true.,sorted_opArray=sorted_opArray)

!!   write(*,*)'write sorted final',match_index,(sorted_oparray(i),i=1,4),transA

   if (transA=='T') then

  allocate(im_temp_temp(dot_product(x_row,y_column)*rcw))
   call srt1c1n (nrep,x_row,y_column,im_temp,im_temp_temp)

 deallocate(im_temp)
   call sorting_for_blas(match_index,.true.,transA,right_sort,x_row,x_column,im_temp_temp,if_TotSym=.false.)

 deallocate(im_temp_temp)
   else
   call sorting_for_blas(match_index,.true.,transA,right_sort,x_row,x_column,im_temp,if_TotSym=.false.)
 deallocate(im_temp)
   endif
   call xaxpy (dot_product(x_row,x_column),A1*asym_factor*dia_sum_factor,right_sort,1,ImTensor,1)

  deallocate(right_sort)

  end subroutine

 integer function dotprod_symmpack(row,column) 

#include "symm.inc"

  integer, intent(in) :: row(nrep),column(nrep)
  integer             :: array_tot(nrep)
  integer             :: irep,jrep 

   do irep = 1, nrep

       jrep = multb(irep+nrep,1+nrep,2)
       array_tot(irep+1) = array_tot(irep) + row(irep)*column(irep)  
       
   enddo
      
       dotprod_symmpack = array_tot(nrep+1) 

 end function

  end module
