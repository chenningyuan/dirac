  module wrapper_sorting

    implicit none

    public sort_22_unrestricted
!    public sort_13
!    public sort_31   
!    public contraction_422

     private
     shuffle_


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

  type(sorting_list)::a
  a%number = 1000*reorder_v(1) + 100*reorder_v(2) + 10*reorder_v(3) + reorder_v(4)  

   do 
    if (a%number == reorder_index(i)) then
    save_index = i
    no_match = .FALSE.
    exit
     i = i+1
    if (i > 5) no_match = .TRUE.
    exit
   enddo 

  if (.no_match) then
    
  do 

   call shuffle_external(reorder_v, no_match)
   
   if (.not.no_match) exit
   
   call shuffle_internal(reorder_v, no_match)

   if (.not.no_match) exit
 
   call shuffle_external_internal(reorder_v, no_match)

   if (.not.no_match) exit

   write (*,*) "no matching sorting subroutine available"

  enddo 
  
  if (.not.no_match)

!call srt%a 
   call srt%reorder_index(save_index) 



   end subroutine

   subroutine shuffle_external(reorder_array, no_matchi, save_index)

! ------- variables --------
    integer, intent(inout) :: reorder_array(4)
    integer :: temp  
    logical, intent(inout) :: no_match
! --------------------------------
   temp = reorder_array(1)
   reorder_array(1) = reorder_array(2)
   reorder_array(2) = temp
   
   a%number = 1000*reorder_array(1) + 100*reorder_array(2) + 10*reorder_array(3) + reorder_array(4)
   
   call match_array(no_match)   

   end subroutine

  subroutine match_array(no_match, save_index)

!-----------variable-------------
  logical, intent(out) :: no_match
  integer, intent(out) :: save_index 
!--------------------------------
   do 
    if (a%number == reorder_index(i)) then
    save_index = i
    no_match = .FALSE.
    exit
     i = i+1
    if (i > 5) no_match = .TRUE.
    exit
   enddo 

  end subroutine

  subroutine

  end module   
