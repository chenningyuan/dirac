  module contraction
  use wrapper_sorting
  use symmetry_offset
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

! 1) fetch 2-e integrals and sort them for BLAS.

    call decide_sorting_routine(13,left_order,right_order,im_order,transA,match_index,asym_factor_L, &
    &   .true.,sorted_opArray=sorted_oparray_left)

     if (present(LeftTensor)) then 
      call sorting_for_blas(13,left_order,match_index,.false.,transA,v_sort,x_row,x_column,LeftTensor)
     else
      call sorting_for_blas(13,left_order,match_index,.false.,transA,v_sort,x_row,x_column)
     endif 

! 2) general sorting routine for fixed arrays. 

     call decide_sorting_routine(31,right_order,im_order,sorted_oparray_left,transB,match_index,asym_factor_R, &
    &         .false.,sorted_opArray=sorted_oparray_right)


     array_equal = all((/sorted_oparray_left(1),sorted_oparray_right(4)/)==im_order)

     if (present(RightTensor)) then 
     call sorting_for_blas(31,right_order,match_index,.false.,transB,right_sorted,y_row,y_column,rightTensor)
     else
     call sorting_for_blas(31,right_order,match_index,.false.,transB,right_sorted,y_row,y_column)
     endif  

! 3) contraction routine for symmetry-packed arrays.

   if (array_equal) then
      call cntrct_new (transA,transB,x_row,y_column,x_column,asym_factor_R*asym_factor_L*dia_factor*A1,v_sort, &
   &                right_sorted,dia_sum_factor*A1,imTensor,1,1,.false.) 

    else

   allocate (ImTensor_temp(dot_product(x_row,y_column)*rcw))
   ImTensor_temp = 0.0d0
     call cntrct_new (transA,transB,x_row,y_column,x_column,asym_factor_R*asym_factor_L*dia_factor*A1,v_sort, &
    &              right_sorted,A0,imTensor_temp,1,1,.false.) 

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
    real(8),intent(in),optional :: LeftTensor(:)                
    real(8),intent(in)     :: rightTensor(*)   
    integer,intent(in)     :: nrep
    real(8),intent(in)     :: dia_sum_factor,dia_factor   !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
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

!-------------------------------------------------------------
!    use of xgemv is different from dgemm because when I mention 'T', then the dimensionality of the incoming array is also changed.
!-------------------------------------------------------------

!--------------------------------------------------------------------------
!  1> fetch the array which is flexible in terms of indices. sort it to make suitable for BLAS.
!     also determine the row and column index for the matrix-vector multiplication routine.

  maxdim = sorted_dimension(left_order)


  allocate(v_sort(maxdim*rcw))  

     transA = 'N'

     call decide_sorting_routine(22,left_order,right_order,im_order,transA,match_index,asym_factor,.true., &
    &                        sorted_opArray=sorted_left)


     if (present(LeftTensor)) then

     call sorting_for_blas(22,left_order,match_index,.false.,transA,v_sort,x_row,x_column,presorted=LeftTensor,if_totsym=.true.)

     else

     call sorting_for_blas(22,left_order,match_index,.false.,transA,v_sort,x_row,x_column,if_totsym=.true.)

     endif

     call onebody_for_blas (right_order,sorted_left(3:4),im_order,transB,y_row,y_column) 

     array_equal = all((/sorted_left(1),sorted_left(2)/)==im_order)

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

!!       call xgemv (transA,x_row(1),x_column(1),asym_factor*dia_factor*A1,v_sort,x_row(1),rightTensor,1,dia_sum_factor*A1,imTensor,1) 
      endif

      deallocate(v_sort)
    end subroutine

   subroutine contraction_242 (left_order,right_order,im_order,LeftTensor,ImTensor, &
           &                   dia_factor,dia_sum_factor,nrep,RightTensor)


!----------calling variables-----------------------------------
    character(2), dimension(2), intent(in)   :: left_order,im_order
    character(2), dimension(4), intent(in)   :: right_order
    real(8), intent(in)                      :: LeftTensor(:)
    real(8), intent(in), optional            :: RightTensor(:)   
    integer, intent(in)                      :: nrep
    real(8), intent(in)                      :: dia_sum_factor        
    real(8), intent(in)                      :: dia_factor            
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
!---------------------------------------------------------------    
#include "param.inc"
#include "complex.inc"

    maxdim = sorted_dimension(right_order)

   allocate(right_sort(maxdim*rcw))  

     transA = 'N'
     trans1A = 'N'

!**  1> fetch the integral array, properly sorted to contract with BLAS.

   call decide_sorting_routine(22,right_order,im_order,left_order,transA,match_index,asym_factor,.true.,&
   &                          sorted_opArray=sorted_oparray_right)


   if (present(RightTensor)) then
   call sorting_for_blas(22,right_order,match_index,.false.,transA,right_sort,y_row,y_column,presorted=RightTensor,if_totsym=.true.)
   else
   call sorting_for_blas(22,right_order,match_index,.false.,transA,right_sort,y_row,y_column,if_totsym=.true.)
   endif

!**  2> sort onebody integral properly to use in BLAS.

   call onebody_for_blas (left_order,im_order,sorted_oparray_right(1:2),trans1A,x_row,x_column)

     array_equal = all((/sorted_oparray_right(3),sorted_oparray_right(4)/)==im_order)

         if (transA.eq.'N'.or.transA.eq.'n') then
            ldb = y_row(1)
         else
            ldb = y_column(1)
         endif

!**  3>  to do B = XA  type operation with xgemv, we are doing B^T = A^T*X^T then transpose(B^T). that implies following modifications.   

   if (trans1A=='T') then

   allocate (left_sort(y_row(1)*rcw))

   call srt1c1n (nrep,x_column,x_row,LeftTensor,left_sort)

   call xgemm ('N',transA,1,y_column(1),y_row(1),asym_factor*dia_factor*A1,left_sort, &
      &       1,right_sort,ldb,dia_sum_factor*A1,imTensor,1)

   deallocate (left_sort)

   else
    
   call xgemm ('N',transA,1,y_column(1),y_row(1),asym_factor*dia_factor*A1,LeftTensor, &
      &       1,right_sort,ldb,dia_sum_factor*A1,imTensor,1)

   endif
   deallocate (right_sort)

    end subroutine

   subroutine contraction_444 (left_order,right_order,im_order,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep,RightTensor,LeftTensor)
    
!----------calling variables-----------------------------------
    character(*), dimension(4), intent(in)   :: left_order, right_order, im_order
    real(8), intent(in), optional :: LeftTensor(:)       
    real(8), intent(in), optional :: RightTensor(:)   
    integer, intent(in) :: nrep
    real(8), intent(in) :: dia_factor,dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)  :: imTensor(*)
!--------------------------------------------------------------
   integer        :: x_row(nrep),x_column(nrep),y_row(nrep),y_column(nrep)
   character(1)   :: transA, transB
   integer        :: save_index_l,save_index_r,save_index,i
   character(2)   :: sorted_opArray_left(4),sorted_opArray_right(4)
   real(8),allocatable :: v_sort(:)
   real(8),allocatable :: right_sorted(:)
   real(8),allocatable :: im_temp(:)
   real(8),allocatable :: im_temp_temp(:)
   integer             :: asym_factor_R,asym_factor_L
   integer             :: maxdim1, maxdim2
   integer             :: contraction_type

#include "param.inc"
#include "complex.inc"

  maxdim1 = max(sorted_dimension(left_order),sorted_dimension(im_order))
  maxdim2 = sorted_dimension(right_order)

  allocate(v_sort(maxdim1*rcw))  
  allocate(right_sorted(maxdim2*rcw))  

   transA = 'N'
   transB = 'N'

!  1) fetch/provide the left tensor and sort properly for BLAS. also extract row and column dimensions of the sorted tensor.

    call decide_sorting_routine(22,left_order,right_order,im_order,transA,save_index_l,asym_factor_L,.true., &
    &  sorted_opArray=sorted_opArray_left)

    if (present(LeftTensor)) then 
    call sorting_for_blas(22,left_order,save_index_l,.false.,transA,v_sort,x_row,x_column, presorted=LeftTensor,if_TotSym=.false.)
    else
    call sorting_for_blas(22,left_order,save_index_l,.false.,transA,v_sort,x_row,x_column,if_TotSym=.false.)
    endif

!  2) fetch/provide the right tensor and sort properly for BLAS. also extract row and column dimensions of the sorted tensor.

    call decide_sorting_routine(22,right_order,im_order,sorted_opArray_left,transB,save_index_r,asym_factor_R,&
    &          .false.,sorted_opArray=sorted_opArray_right)

    if (present(RightTensor)) then 
    call sorting_for_blas(22,right_order,save_index_r,.false.,transB,right_sorted,y_row,y_column,&
   &                  presorted=RightTensor,if_TotSym=.false.)
    else
    call sorting_for_blas(22,right_order,save_index_r,.false.,transB,right_sorted,y_row,y_column,if_TotSym=.false.)
    endif

!  3) Distinguish the contraction cases based on symmetry ordering of the sorted tensors.
!     This disntiction is required since some of them use multb(*,*,1) for final
!     packing and some others multb(*,*,2). As of now I am hard-coding them, until
!     I get a better solution.

    if ((save_index_l==2).and.(save_index_r==1)) then
      contraction_type = 2 
    elseif ((save_index_l==1).and.(save_index_r==2)) then
      contraction_type = 3
    elseif ((save_index_l==2).and.(save_index_r==2)) then
      contraction_type = 4
    else
      contraction_type = 1
    endif  

    write(*,*)'what is the contraction type?',contraction_type

!  4) matrix-matrix multiplication routine runs for each irreps.

   if (dia_sum_factor==0.0d0) then

    allocate(im_temp(dot_product(x_row,y_column)*rcw))

    im_temp = 0.0d0

   select case(contraction_type)

    case(1)

    call cntrct_new (transA,transB,x_row,y_column,x_column,asym_factor_L*asym_factor_R*dia_factor*A1,v_sort,&
    &    right_sorted,A0,im_temp,1,1,.true.) 

    case(2)

    call cntrct_1 (transA,transB,x_row,y_column,x_column,asym_factor_L*asym_factor_R*dia_factor*A1,v_sort, &
    &      right_sorted,A0,im_temp,1,1,.true.) 

    case(3)

    call cntrct_2 (transA,transB,x_row,y_column,x_column,asym_factor_L*asym_factor_R*dia_factor*A1,v_sort, &
    &              right_sorted,A0,im_temp,1,1,.true.) 
    case(4)

    call cntrct_3 (transA,transB,x_row,y_column,x_column,asym_factor_L*asym_factor_R*dia_factor*A1,v_sort, &
    &           right_sorted,A0,im_temp,1,1,.true.) 

  end select  

!  5) resort and antisymmetrize the product tensor. sum that up with the particular class of intermediate/residue/sigma etc. 

   transA = 'N'

   v_sort = 0.0d0

   call decide_sorting_routine(22,im_order,left_order,(/sorted_opArray_left(1:2),sorted_opArray_right(3:4)/), &
                       & transA,save_index,asym_factor_R,.false.,sorted_opArray=sorted_opArray_right,inverse=.true.)

   if (transA=='T') then
    write(*,*)'possibility of mistake',x_row(1),y_column(1),x_row(2),y_column(2)
   allocate(im_temp_temp(dot_product(x_row,y_column)*rcw))

   call srt1c1n (nrep,x_row,y_column,im_temp,im_temp_temp)

   call sorting_for_blas(22,im_order,save_index,.true.,transA,v_sort,x_row,x_column,im_temp_temp,if_TotSym=.false.)

     deallocate(im_temp_temp)
   else

   call sorting_for_blas(22,im_order,save_index,.true.,transA,v_sort,x_row,x_column,im_temp,if_TotSym=.false.)

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
    real(8), intent(in) :: dia_factor,dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout) :: ImTensor(*)
    integer      :: x_row(nrep),x_column(nrep),y_row(nrep),y_column(nrep)
    character(1) :: transA, transB
    integer      :: match_index,i
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

   transA = 'N'
   transB = 'N'


!  1) fetch/provide the left tensor and sort properly for BLAS. also extract row and column dimensions of the sorted tensor.

   call decide_sorting_routine(31,left_order,right_order,im_order,transA,match_index,asym_factor, &
   &             .true.,sorted_opArray=sorted_opArray)


 if (present(LeftTensor)) then

   call sorting_for_blas(31,left_order,match_index,.false.,transA,v_sort,x_row,x_column,presorted=LeftTensor)
 
 else

   call sorting_for_blas(31,left_order,match_index,.false.,transA,v_sort,x_row,x_column)

 endif

   call onebody_for_blas (right_order,sorted_opArray,im_order,transB,y_row,y_column,FreeIndx)

!  2) matrix-matrix multiplication routine runs for each irreps.

   allocate(im_temp(dot_product(x_row,y_column)*rcw))
   im_temp = 0.0d0

    if (present(RightTensor)) call cntrct_new (transA,transB,x_row,y_column,x_column,asym_factor*dia_factor*A1,v_sort, &
    &                                      RightTensor,A0,im_temp,1,1,.false.)

!  3) re-sort the arrays to align with the desired structure of the intermediate.

   transA = 'N'

   call decide_sorting_routine(31,im_order,left_order,(/sorted_opArray(1:3),FreeIndx(2)/),transA,match_index, &
           &                    asym_factor,.false.,inverse=.true.,sorted_opArray=sorted_opArray)

   if (transA=='T') then

   allocate(im_temp_temp(dot_product(x_row,y_column)*rcw))

   call srt1c1n (nrep,x_row,y_column,im_temp,im_temp_temp)

!!      if (carith) call conjuga(dot_product(x_row,y_column),im_temp_temp,1)

   deallocate(im_temp)

   call sorting_for_blas(31,im_order,match_index,.true.,transA,v_sort,x_row,x_column,im_temp_temp,if_TotSym=.false.)

   deallocate(im_temp_temp)

   else

   call sorting_for_blas(31,im_order,match_index,.true.,transA,v_sort,x_row,x_column,im_temp,if_TotSym=.false.) 

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

   transA = 'N'
   transB = 'N'
!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

   call decide_sorting_routine(13,right_order,im_order,left_order,transB,match_index,asym_factor, &
    &                    .true.,sorted_opArray=sorted_opArray)

 if (present(RightTensor)) then
 
   call sorting_for_blas(13,right_order,match_index,.false.,transB,right_sort,y_row,y_column,RightTensor)
 
 else

   call sorting_for_blas(13,right_order,match_index,.false.,transB,right_sort,y_row,y_column)

 endif

!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.

   call onebody_for_blas (left_order,im_order,sorted_opArray,transA,x_row,x_column,FreeIndx)

!  3) matrix-matrix multiplication routine runs over each irreps.

   allocate(im_temp(dot_product(x_row,y_column)*rcw))
   im_temp = 0.0d0

! if (present(LeftTensor)) call cntrct (transA,transB,x_row,y_column,y_row,asym_factor*dia_factor*A1,LeftTensor,&
! &        right_sort,A0,im_temp,nrep) 


  if (present(LeftTensor)) call cntrct_new (transA,transB,x_row,y_column,y_row,asym_factor*dia_factor*A1,LeftTensor,&
  &        right_sort,A0,im_temp,1,1,.false.) 

!!   right_sort = 0.0d0
!  4) re-sort the arrays to align with the desired structure of the imtermediate.

   transA = 'N'

   call decide_sorting_routine(13,im_order,left_order,(/FreeIndx(1),sorted_opArray(2:4)/),transA,match_index, &
   &                        asym_factor,.false.,inverse=.true.,sorted_opArray=sorted_opArray)

   if (transA=='T') then

  allocate(im_temp_temp(dot_product(x_row,y_column)*rcw))
   call srt1c1n (nrep,x_row,y_column,im_temp,im_temp_temp)

 deallocate(im_temp)
   call sorting_for_blas(13,im_order,match_index,.true.,transA,right_sort,x_row,x_column,im_temp_temp,if_TotSym=.false.)

 deallocate(im_temp_temp)
   else
   call sorting_for_blas(13,im_order,match_index,.true.,transA,right_sort,x_row,x_column,im_temp,if_TotSym=.false.)
 deallocate(im_temp)
   endif
   call xaxpy (dot_product(x_row,x_column),A1*asym_factor*dia_sum_factor,right_sort,1,ImTensor,1)

  deallocate(right_sort)

  end subroutine

!===========================================================================================

      subroutine cntrct_new (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine is
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*),c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim(nrep),ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)

      shift = 0
      if(shift_for_boson) shift=nrep 

       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif
       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif
        call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)

      do krp = 1, nrep
         mrp = multb(krp+shift,irrepa+nrep,1)  
         nrp = multb(krp+shift,irrepb+nrep,1)  

         m = mdim(mrp)
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
          off1 = 1 + e%twoDirac(mrp,krp)*rcw
         else 
          off1 = 1 + e%twoDirac(krp,mrp)*rcw
         endif

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twoDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twoDirac(nrp,krp)*rcw
         endif

         off3 = 1 + g%twoDirac(mrp,nrp)*rcw

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
      enddo

      call dealloc_array(e,f,g)
      end subroutine

      subroutine cntrct_1 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine is
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*),c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim(nrep),ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)
      shift = 0
      if(shift_for_boson) shift=nrep 

       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif
       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif
        call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)

      do krp = 1, nrep
         mrp = multb(krp+shift,irrepa+nrep,1)  
         nrp = multb(krp+shift,irrepb+nrep,2)  

         m = mdim(mrp)
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
          off1 = 1 + e%twoDirac(mrp,krp)*rcw
         else 
          off1 = 1 + e%twoDirac(krp,mrp)*rcw
         endif

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twononDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twononDirac(nrp,krp)*rcw
         endif

         off3 = 1 + g%twononDirac(mrp,nrp)*rcw

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
      enddo
     
      call dealloc_array(e,f,g)

      end subroutine

      subroutine cntrct_2 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine is
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*),c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim(nrep),ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)
      shift = 0
      if(shift_for_boson) shift=nrep 

       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif
       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif
        call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)

      do krp = 1, nrep
         mrp = multb(krp+shift,irrepa+nrep,2)  
         nrp = multb(krp+shift,irrepb+nrep,1)  

         m = mdim(mrp)
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
          off1 = 1 + e%twononDirac(mrp,krp)*rcw
         else 
          off1 = 1 + e%twononDirac(krp,mrp)*rcw
         endif

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twoDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twoDirac(nrp,krp)*rcw
         endif

         off3 = 1 + g%twononDirac(mrp,nrp)*rcw
!          off3 = 1 + g%twoDirac(mrp,nrp)*rcw

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
      enddo

      call dealloc_array(e,f,g)
      end subroutine

      subroutine cntrct_3 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine is
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*),c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim(nrep),ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)
      shift = 0
      if(shift_for_boson) shift=nrep 

       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif
       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif
        call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)

      do krp = 1, nrep
         mrp = multb(krp+shift,irrepa+nrep,2)  
         nrp = multb(krp+shift,irrepb+nrep,2)  

         m = mdim(mrp)
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
          off1 = 1 + e%twononDirac(mrp,krp)*rcw
         else 
          off1 = 1 + e%twononDirac(krp,mrp)*rcw
         endif

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twononDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twononDirac(nrp,krp)*rcw
         endif

         off3 = 1 + g%twononDirac(mrp,nrp)*rcw


         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
      enddo

      call dealloc_array(e,f,g)
      end subroutine

  end module
