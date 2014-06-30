  module wrapper_sorting

    implicit none

    public sort_22_unrestricted
!    public sort_13
!    public sort_31   
    public shuffle_indices
    public match_array
    public v_for_blas


    private

     type sorting_list
       character(20)::unrestricted, restricted(3) 
       integer :: num 
       integer :: trans 
     endtype sorting_list 

    CONTAINS

    subroutine sort_22_unrestricted(vorder,torder, worder,vtrans)

!----------calling variables-----------------------------------
    character (*), dimension(4), intent(in) :: vorder
    character (*), dimension(:), intent(in) :: worder, torder
    character (*), intent(out) :: vtrans
!--------------------------------------------------------------

! ------- local variables --------
    integer      :: reorder_v(4)
    integer      :: i,j,k,l,dim_t,save_index  
    logical      :: t_transpose

! --------------------------------

      call v_for_blas (vorder, worder, torder, vtrans)

!      k = 1 
!      l = 3
!    do i = 1, 4       ! loop over vorder 
!      do j = 1, 2     ! loop over worder 
!         
!           if (vorder(i).eq.worder(j)) then
!              reorder_v(k) = i
!              k = k+1
!           else 
!              reorder_v(l) = i
!              l = l+1 
!           endif
!         
!      enddo
!    enddo  

!   dim_t = size(torder,1)
!   
!  if (dim_t.eq.2) then
!     if (t_transpose) then
!      temp = reorder_v(3)
!      reorder_v(3)=reorder_v(4)
!      reorder_v(4)=temp
!     endif
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
!  endif  

! do i = 1,4

!  if V_DUMMY_ORDEr(i).ne.

!  type(sorting_list)::a

!  a%num = 1000*reorder_v(1) + 100*reorder_v(2) + 10*reorder_v(3) + reorder_v(4)  

!  call match_array (a%num, no_match, save_index)  

! if (.no_match) then
!   
! do 

!  call shuffle(external, reorder_v, no_match, save_index)
!  
!  if (.not.no_match) exit
!  
!  call shuffle(internal,reorder_v, no_match,save_index)

!  if (.not.no_match) exit
!
!  call shuffle(external_and_internal,reorder_v, no_match,save_index)

!  if (.not.no_match) exit

!  write (*,*) "no matching sorting subroutine available"

! enddo 
! 
! if (.not.no_match)

!call srt%a 
!   call srt%reorder_index(save_index) 
!   write(*,*) srt%reorder_index(save_index) 

   end subroutine

   subroutine shuffle_indices(arraytype, shuffle, reorder_array, no_match, save_index, vtrans)

! ------- variables --------
    integer, intent(inout)   :: reorder_array(4)
    logical, intent(inout)   :: no_match
    character(*), intent(in) :: shuffle, arraytype
    character(*), intent(out):: vtrans
    integer                  :: temp  
    integer, intent(out)     :: save_index
    type(sorting_list)       :: a
! --------------------------

  select case(shuffle)
    
    case('external')

     temp = reorder_array(1)
     reorder_array(1) = reorder_array(2)
     reorder_array(2) = temp
   
     a%num = 1000*reorder_array(1) + 100*reorder_array(2) + 10*reorder_array(3) + reorder_array(4)

     a%trans = 1000*reorder_array(3) + 100*reorder_array(4) + 10*reorder_array(1) + reorder_array(2)  

    case('internal')

     temp = reorder_array(3)
     reorder_array(3) = reorder_array(4)
     reorder_array(4) = temp
   
     a%num = 1000*reorder_array(1) + 100*reorder_array(2) + 10*reorder_array(3) + reorder_array(4)

     a%trans = 1000*reorder_array(3) + 100*reorder_array(4) + 10*reorder_array(1) + reorder_array(2)  
    case('external_and_internal') 

     temp = reorder_array(1)
     reorder_array(1) = reorder_array(2)
     reorder_array(2) = temp

     temp = reorder_array(3)
     reorder_array(3) = reorder_array(4)
     reorder_array(4) = temp

     a%num = 1000*reorder_array(1) + 100*reorder_array(2) + 10*reorder_array(3) + reorder_array(4)
     a%trans = 1000*reorder_array(3) + 100*reorder_array(4) + 10*reorder_array(1) + reorder_array(2)  
  end select
   
   call match_array (arraytype, a%num, a%trans, no_match, save_index, vtrans)  

   end subroutine

  subroutine match_array (arraytype, permute, permute_transpose, no_match, save_index, vtrans) 

!-----------variable-------------
   integer, intent(in) :: permute, permute_transpose
   logical, intent(inout):: no_match
   integer, intent(out):: save_index
   character(*), intent(out) :: vtrans
   character(*), intent(in)  :: arraytype   
!--------------------------------
!---------local variable---------
  integer, allocatable :: A(:)
!  integer, dimension(6) :: B = (/1234, 1324, 2134, 3412, 1432, 1243/)
!  integer, dimension(6) :: C = (/1234, 1324, 2134, 3412, 1432, 1243/)
  integer :: i, dim_A
!--------------------------------

   vtrans='N'

  select case(arraytype)

  case('unrestricted')

  allocate(A(6))

  A = (/1234, 1324, 2134, 3412, 1432, 1243/)

  case('restricted_braket')

  allocate(A(2))
  A = (/1234,1324/)
 
  case('restricted_bra')

  allocate(A(2))
  A = (/1234,1432/)
 
  case('restricted_ket')

  allocate(A(4))
  A = (/1234,1324,1432,1342/)

  end select

  dim_A = size(A)

   i = 1
   do 
   if (permute == A(i)) then
   save_index = i
   no_match = .FALSE.
   exit
   endif
     i = i+1
   if (i > dim_A) exit
   enddo 

   if (no_match) then
   i = 1
   do 
    if (permute_transpose == A(i)) then 
    save_index = i
    no_match = .FALSE.
    vtrans ='C'  
    exit
    endif
     i = i+1
    if (i > dim_A) exit
   enddo 
   endif

  deallocate(A)

  end subroutine

  subroutine v_for_blas (array_left, worder, torder, vtrans)

!--------------calling variables---------------
    character (*), dimension(4), intent(in) :: array_left
    character (*), dimension(:), intent(in) :: worder, torder
    character (*), intent(out)              :: vtrans
!----------------------------------------------

!-------local variables--------
    character(2), dimension(4) :: array_exchange
    logical      :: v_direct, v_exchange
    integer      :: reorder_v(4)
    integer      :: i,j,k,l,dim_t,save_index,temp 
    logical      :: no_match, t_transpose  
    character(2) :: p,q,r,s
    character(20):: arraytype
    type(sorting_list) :: a

    type(sorting_list) :: srt(0:5)
       srt(1)%unrestricted  = 'SRT1SS4'   
       srt(2)%unrestricted  = 'SRT1L1'    
       srt(3)%unrestricted  = 'SRT1C1'    
       srt(4)%unrestricted  = 'SRT16'     
       srt(5)%unrestricted  = 'SRT1R1'
       srt(1)%restricted(1) = 'SRT1TT4'
!      srt(1)%restricted(2) = 'SRT1LS1'
       srt(1)%restricted(2) = 'SRT36'
       srt(1)%restricted(3) = 'SRT22'
       srt(2)%restricted(3) = 'SRT26'
       srt(3)%restricted(3) = 'SRT20D'
     
!------------------------------

     arraytype = 'unrestricted'

     p = array_left(1)
     q = array_left(2)
     r = array_left(3)
     s = array_left(4)

    if (p(1:1) == q(1:1)) then
      if (r(1:1) == s(1:1)) then
       arraytype = 'restricted_braket'
      else
       arraytype = 'restricted_bra'
      endif 
      elseif (r(1:1)==s(1:1)) then
       arraytype = 'restricted_ket'
    endif

    dim_t = size(torder)
    v_direct = .TRUE.
    v_exchange = .FALSE.
    no_match = .TRUE.

    if (v_direct) then
       
       k = 1
       l = 3
     do i = 1, 4       ! loop over vorder 
       do j = 1, 2     ! loop over worder 

            if (array_left(i).eq.worder(j)) then
               reorder_v(k) = i
               k = k+1
               exit
               elseif (array_left(i).eq.torder(j)) then !CAUTION. this is not completely general. I have to change it.
               reorder_v(l) = i
               l = l+1
            endif
       
       enddo
     enddo

     t_transpose = (torder(1)/=array_left(reorder_v(3)))

    if (dim_t == 2) then
      if (t_transpose) then
       temp = reorder_v(3)
       reorder_v(3)=reorder_v(4)
       reorder_v(4)=temp
      endif
    endif

   a%num = 1000*reorder_v(1) + 100*reorder_v(2) + 10*reorder_v(3) + reorder_v(4)

   a%trans = 1000*reorder_v(3) + 100*reorder_v(4) + 10*reorder_v(1) + reorder_v(2)  


   call match_array (arraytype, a%num, a%trans, no_match, save_index, vtrans)

   if (no_match) call shuffle_indices (arraytype,'external', reorder_v, no_match, save_index,vtrans)

!   call shuffle(internal,reorder_v, no_match,save_index)

!   if (.not.no_match) exit

!   call shuffle(external_and_internal,reorder_v, no_match,save_index)

!   if (.not.no_match) exit

     if (no_match) then
      write (*,*) "no matching sorting subroutine available for direct mode of contraction."
      v_exchange = .TRUE.
     endif
  
   endif 

   if (v_exchange) then

    array_exchange = array_left

    call swap_character (array_exchange,3,4) 

    arraytype = 'unrestricted'

     p = array_left(1)
     q = array_left(2)
     r = array_left(3)
     s = array_left(4)

    if (p(1:1) == q(1:1)) then
      if (r(1:1) == s(1:1)) then
       arraytype = 'restricted_braket'
      else
       arraytype = 'restricted_bra'
      endif 
      elseif (r(1:1)==s(1:1)) then
       arraytype = 'restricted_ket'
    endif

   write(*,*)'contraction will be done in exchange mode. It suggests we will permute only contracted indices....'

       k = 1
       l = 3
     do i = 1, 4       ! loop over vorder 
       do j = 1, 2     ! loop over worder 

            if (array_exchange(i).eq.worder(j)) then
               reorder_v(k) = i
               k = k+1
               exit
               elseif (array_exchange(i).eq.torder(j)) then !this is not completely general. I have to change it.
               reorder_v(l) = i
               l = l+1
            endif
       
       enddo
     enddo

     t_transpose = (torder(1)/=array_exchange(reorder_v(3)))

    if (dim_t == 2) then
      if (t_transpose) then
       temp = reorder_v(3)
       reorder_v(3)=reorder_v(4)
       reorder_v(4)=temp
      endif
    endif

     a%num = 1000*reorder_v(1) + 100*reorder_v(2) + 10*reorder_v(3) + reorder_v(4)

     a%trans = 1000*reorder_v(3) + 100*reorder_v(4) + 10*reorder_v(1) + reorder_v(2)  

     call match_array (arraytype,a%num, a%trans, no_match, save_index,vtrans)

    if (no_match) call shuffle_indices (arraytype,'external', reorder_v, no_match, save_index, vtrans)

    endif

     if (no_match) then
      write (*,*) "no matching sorting subroutine available. have to do something different...."
     endif

    if (.not.no_match) then

     save_index = save_index - 1

     if (save_index > 0) then

     if (arraytype == 'unrestricted') then
 !call srt%a 
     write(*,*) srt(save_index)%unrestricted
     endif

     if (arraytype == 'restricted_braket') then
 !call srt%a 
 !   call srt%reorder_index(save_index) 
     write(*,*) srt(save_index)%restricted(1)
     endif

     if (arraytype == 'restricted_bra') then
 !call srt%a 
 !   call srt%reorder_index(save_index) 
     write(*,*) srt(save_index)%restricted(2)
     endif

     if (arraytype == 'restricted_ket') then
 !call srt%a 
 !   call srt%reorder_index(save_index) 
     write(*,*) srt(save_index)%restricted(3)
     endif

     else
     write(*,*)'sorting is not needed'
     endif
   endif
 
  end subroutine
  subroutine swap_character(array,k,l)

  character(*),intent(inout) :: array(*)
  integer, intent(in)        :: k,l  
  character(2)               :: temp

  temp = array(k)

  array(k) = array(l)
  array(l) = temp

  end subroutine

  end module   
