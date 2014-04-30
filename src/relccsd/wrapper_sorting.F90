  module wrapper_sorting

    implicit none

    public sort_22_unrestricted
!    public sort_13
!    public sort_31   
!    public contraction_422

     private


     contains


    subroutine sort_22_unrestricted(vorder,v_dummy_order,torder, worder)

!----------calling variables-----------------------------------
    character (*), dimension(4), intent(in) :: vorder, v_dummy_order
    character (*), dimension(*), intent(in) :: worder, torder
!--------------------------------------------------------------

! ------- local variables --------
    integer :: reorder_v(4)
    integer :: i,j,k,l,dim_t,temp  
    logical :: t_transpose
! --------------------------------

!    SRT1SS4 = 1324
!    SRT1L1  = 2134
!    SRT1C1  = 3412
!    SRT16   = 1432
!    SRT1R1  = 1243  


   type sorting_list
   integer, dimension(5) :: reorder_index = (/1234, 1324, 2134, 3412, 1432, 1243/)
   integer :: number  
   endtype sorting_list 

   type(sorting_list) ::  srt

       srt%reorder_index(1) =  SRT1SS4   
       srt%reorder_index(2) =  SRT1L1    
       srt%reorder_index(3) =  SRT1C1    
       srt%reorder_index(4) =  SRT16     
       srt%reorder_index(5) =  SRT1R1    


     t_transpose = (torder(1).ne.vorder(reorder_v(3)))


       k = 1 
       l = 3
     do i = 1, 4       ! loop over vorder 
       do j = 1, 2     ! loop over worder 
          
            if (vorder(i).eq.worder(j)) then
               reorder_v(k) = i
               k = k+1
            else 
               reorder_v(l) = i
               l = l+1 
            endif
           
       enddo
     enddo  

    dim_t = size(torder,1)
    
   if (dim_t.eq.2) then
      if (t_transpose) then
       temp = reorder_v(3)
       reorder_v(3)=reorder_v(4)
       reorder_v(4)=temp
      endif
!   else 
!   if (dim_t.eq.4) then
!      
!   do i = 1, 4       ! loop over vorder 
!      do j = 1, 2     ! loop over worder 

!           if (vorder(i).eq.worder(j)) then
!              reorder_v(k) = i
!              k = k+1
!           else
!              reorder_v(l) = i
!              l = l+1
!           endif

!      enddo
!    enddo
!


!   endif
   endif  

! do i = 1,4

!  if V_DUMMY_ORDEr(i).ne.

!  
!  type(sorting_list)::a
!  a%number = 1000*reorder_v(1) + 100*reorder_v(2) + 10*reorder_v(3) + reorder_v(4)  

!  do 
!  
!  if (a%number .eq.  
     
     

    end subroutine
 
    

 
  end module   
