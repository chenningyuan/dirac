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

   subroutine contraction_442 (varray, tarray, warray, vorder, torder, worder, &
           &        ext_v, ext_t, index_sum, dia_factor, dia_sum_factor,buf2,buf1)

     
!----------calling variables-----------------------------------
    real(8), intent(in) :: varray(*), tarray(*)      
    complex(8), intent(in) :: dia_factor, dia_sum_factor               !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
!   character(*), intent(in) :: vorder, torder, worder
    character (*), dimension(4), intent(in) :: vorder, torder
    character (*), dimension(2), intent(in) :: worder
    real(8), intent(in) :: buf2(*), buf1(*)      
    integer, intent(in) :: ext_v(nrep), ext_t(nrep), index_sum(nrep)
    real(8), intent(out) :: warray(*)
!--------------------------------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
   
   
      if (worder(1).eq.'h1'.and.worder(2).eq.'h2') then
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,varray,BUF2)
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,tarray,BUF1)
      endif
 
      if (worder(1).eq.'p1'.and.worder(2).eq.'p2') then
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT, &
     &             varray,BUF2)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT, &
     &             tarray,BUF1)
      endif

     CALL CNTRCT ('C','N',ext_v,ext_t,index_sum,dia_factor,BUF2,BUF1,dia_sum_factor,warray,NREP) 
 
  
    end subroutine


   subroutine contraction_422 (left_order, right_order, im_order, right_array, im_array, &
           &                   dia_factor, dia_sum_factor, nrep)

     
!----------calling variables-----------------------------------
    character(*), dimension(4), intent(in)   :: left_order, right_order
    character(*), dimension(2), intent(in)   :: im_order
!    real(8), intent(in), optional :: left_array(:), right_array(:)   !when the left or right array is a fixed array then we will provide input.    
    real(8), intent(in) :: right_array(*)   
    integer, intent(in) :: nrep
    real(8), intent(in)        :: dia_factor, dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)     :: im_array(*)
!--------------------------------------------------------------
    integer  :: x_row(nrep), x_column(nrep)
    character (1)  :: trans
    real(8)  :: v_sort(100000)
#include "param.inc"
!  1) fetch the fixed array, properly sorted to contract with BLAS.

!  2) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS.

      call v_for_blas (left_order,right_order,im_order, trans, v_sort, x_row, x_column)

!  3) determine the row and column index for the matrix-vector multiplication routine. 
!     I will like to fetch a flexible array   

      CALL XGEMV (trans,x_row,x_column,dia_factor*A1,v_sort,x_row,right_array,1,dia_sum_factor*A1,im_array,1) 
  
    end subroutine


  end module
