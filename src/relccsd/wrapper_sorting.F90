  module wrapper_sorting

    implicit none

    public shuffle_indices
    public match_array
    public v_for_blas
    public T_for_blas
    public onebody_for_blas
    public sorting_for_blas
    public decide_sorting_routine 



    private

    integer      :: sortType
    logical      :: do_direct, do_exchange
    logical      :: t_direct, t_exchange
    character(20):: arraytype
    character(2),pointer :: array_original(:)
    character(2)         :: exchange(4)
    integer              :: reorder_left(4) 
    type sorting_list
      integer :: num 
      integer :: trans
!     contains
      procedure(),pointer, nopass :: unrestricted !, restricted 
      procedure(),pointer, nopass :: restricted_braket, restricted_bra, restricted_ket 
    endtype sorting_list 

    type indices
      integer, dimension(:,:),pointer :: tuple
      integer, dimension(:),pointer   :: tuple_combInd, orbital_type, full_array 
    endtype indices

    type, extends(indices) :: middle_indices
       type(indices) :: triangular
    endtype

   type, extends(middle_indices) :: more_indices
      type(middle_indices) :: pqr,qrs
   end type

    type bit_processing
       integer :: bra,ket,braKet,I,J,K,L,ijk,jkl,ijkl
       integer :: I_reorder,J_reorder,K_reorder,L_reorder,ijkl_reorder
       integer :: bra_reorder,ket_reorder,ijk_reorder,jkl_reorder
    endtype bit_processing 

    type presorting
       real*8, allocatable :: presorted(:) 
    endtype
   external SRT1SS4
   external SRT1L1
   external SRT1C1
   external SRT16
   external SRT1R1
   external SRT1TT4
   external SRT36 
   external SRT22
   external SRT26
   external SRT20D
   external SRT1T3
   external SRT1T2
   external SRT1S2
   external SRT19
   external SRT1S3
   external SRT7
   external SRT1TS4

    CONTAINS
    
   subroutine shuffle_indices(reorder_array,no_match,save_index,trans)

! ------- variables --------
    integer,intent(inout),target :: reorder_array(4)
    logical,intent(inout)        :: no_match
    character(1),intent(inout)   :: trans
    integer,pointer              :: shifted_array(:)
    integer,allocatable          :: P(:,:)
    integer                      :: temp,i  
    integer, intent(out)         :: save_index
    type(sorting_list)           :: a
! --------------------------

  select case(sortType)
    
    case(13)

     shifted_array => reorder_array(2:4)

     temp = reorder_array(1)
!     reorder_array(1) = reorder_array(2)
!     reorder_array(2) = temp

!    a%num = toDecimal(reorder_array)
!    a%trans = toDecimal(reorder_array,.true.)    

    allocate (P(6,3))

    call permutate(shifted_array,P)

    do i = 1, size(P,1)
      a%num = toDecimal((/temp,P(i,:)/))
      a%trans = toDecimal((/temp,P(i,:)/),.true.)    

      call match_array (a%num, no_match, save_index, trans)  
      if (no_match) call match_array_transpose (a%trans, no_match, save_index, trans)  

      if (.not.no_match) then
        reorder_array(2:4)=P(i,:)
        exit
      endif

    enddo
    deallocate(P)

     case(31) 

     shifted_array => reorder_array(1:3)

     temp = reorder_array(4)

    allocate (P(6,3))

    call permutate(shifted_array,P)

    do i = 1, size(P,1)
!    write(*,*)P(i,:)
    a%num = toDecimal((/P(i,:),temp/))
    a%trans = toDecimal((/P(i,:),temp/),.true.)    

    call match_array (a%num, no_match, save_index, trans)  
    if (no_match) call match_array_transpose (a%trans, no_match, save_index, trans)  
    if (.not.no_match) then
    reorder_array(1:3)=P(i,:)
    exit
    endif
    enddo
    deallocate(P)

   case (22) 

     temp = reorder_array(1)
     reorder_array(1) = reorder_array(2)
     reorder_array(2) = temp

    a%num = toDecimal(reorder_array)
    a%trans = toDecimal(reorder_array,.true.)    

   call match_array (a%num,no_match,save_index,trans)  
   if (no_match) call match_array_transpose (a%trans,no_match,save_index,trans)  
  end select
   

   end subroutine

  subroutine match_array (permute,no_match,save_index,trans) 

!-----------calling variable-------------
   integer, intent(in)        :: permute 
   logical, intent(inout)     :: no_match
   integer, intent(out)       :: save_index
   character(1), intent(inout):: trans
!--------------------------------
!---------local variable---------
   integer, allocatable       :: A(:)
   integer                    :: i,dim_A,upbound_A
!--------------------------------

  select case(arraytype)

  case('unrestricted')

  select case(sortType)

  case(22) 

  allocate(A(6))

  A = (/1234, 1324, 2134, 3412, 1432, 1243/)

  case(31)

  allocate(A(2))

  A = (/1234,1243/)

  case(13)

  allocate(A(1))

  A = (/1234/)

  end select

  case('restricted_braket')

  select case(sortType) 

  case(22)
  allocate(A(2))
  A = (/1234,1324/)

  case(31)
 
  allocate(A(3))

  A= (/1234,4123,1342/)
  
  case(13) 
  allocate(A(1))

  A= (/1234/)

  end select
 
  case('restricted_bra')

  select case(sortType) 

  case(22)
  allocate(A(2))
  A = (/1234,1432,1324/)

  case(31)
  allocate(A(2))
  A = (/1342,1243/)

  case(13)
  allocate(A(1))
  A = (/1234/)

  end select
 
  case('restricted_ket')

  select case(sortType)

  case(22)

  allocate(A(4))
  A = (/1234,1324,1432,1342/)

  case(31)

  allocate(A(2))
  A = (/1234,4123/)

  end select

  end select

  dim_A = size(A)
  upbound_A = dim_A

   i = 1

   do 
   if (permute == A(i)) then
   save_index = i
   no_match = .FALSE.
   exit
   endif
     i = i+1
   if (i > upbound_A) exit
   enddo 

  deallocate(A)

  end subroutine

  subroutine match_array_transpose (permute_transpose, no_match, save_index, trans) 

!-----------variable-------------
   integer,intent(in) :: permute_transpose
   logical,intent(inout):: no_match
   integer,intent(out):: save_index
   character(1),intent(inout) :: trans
!--------------------------------
!---------local variable---------
  integer, allocatable :: A(:)
  integer :: i, dim_A, upbound_A
  integer :: sort_type
!--------------------------------

  if (sortType == 31) sort_type = 13
  if (sortType == 13) sort_type = 31

 select case(arraytype)

  case('unrestricted')

  select case(sort_type)

  case(22) 

  allocate(A(6))

  A = (/1234, 1324, 2134, 3412, 1432, 1243/)

  case(31)

  allocate(A(2))

  A = (/1234,1243/)

  case(13)

  allocate(A(1))

  A = (/1234/)

  end select

  case('restricted_braket')

  select case(sort_type) 

  case(22)
  allocate(A(2))
  A = (/1234,1324/)

  case(31)
 
  allocate(A(3))

  A= (/1234,4123,1342/)
  
  case(13) 
  allocate(A(1))

  A= (/1234/)

  end select
 
  case('restricted_bra')

  select case(sort_type) 

  case(22)
  allocate(A(2))
  A = (/1234,1432,1324/)

  case(31)
  allocate(A(1))
  A = (/1342/)

  case(13)
  allocate(A(1))
  A = (/1234/)

  end select
 
  case('restricted_ket')

  select case(sort_type)

  case(22)

  allocate(A(4))
  A = (/1234,1324,1432,1342/)

  case(31)

  allocate(A(2))
  A = (/1234,4123/)

  end select

  end select

  dim_A = size(A)
  upbound_A = dim_A

   i = 1
   do 
    if (permute_transpose == A(i)) then 
      save_index = i
      no_match = .FALSE.
      trans ='C'  
      exit
    endif
     i = i+1
    if (i > upbound_A) exit
   enddo 

  deallocate(A)

  end subroutine

  subroutine v_for_blas (sort_type,opArray_left,opArray_right,opArray_target,if_trans,v_sorted, &
                      &  row,column,sorted_opArray)

!--------------calling variables---------------
    character(*),intent(in),target        :: opArray_left(:)
    character(*),dimension(:),intent(in)  :: opArray_target, opArray_right
    integer,intent(in)                    :: sort_type
    character(*),intent(out)              :: if_trans
    integer,intent(out)                   :: row(:), column(:)
    real(8),intent(out)                   :: v_sorted(:) !think a bit more carefully to replace it by dynamic allocation. not easy.
    character(2),intent(out),optional     :: sorted_opArray(size(opArray_left))
!----------------------------------------------

!-------local variables--------
    character(2) :: array_exchange(size(opArray_left))
    integer      :: reorder_v(4)
    integer      :: i,j,k,l,save_index,temp 
    logical      :: no_match, t_transpose  
    character(2) :: p,q,r,s
    integer      :: sort_type_dummy
#include "symm.inc"
#include "param.inc"

    type(sorting_list) :: a
    type(middle_indices)   :: c(0:5)  
    type(more_indices)     :: m(0:15)  
    type(bit_processing)   :: d 
    type(presorting)       :: V 

    type(sorting_list) :: srt_22(0:5)
    type(sorting_list) :: srt_31(0:5)
    type(sorting_list):: srt_13(0:5)
! (2,2)->(2,2) type of sorting routines.         
       srt_22(1)%unrestricted => SRT1SS4   
       srt_22(2)%unrestricted => SRT1L1    
       srt_22(3)%unrestricted => SRT1C1    
       srt_22(4)%unrestricted => SRT16     
       srt_22(5)%unrestricted => SRT1R1

       srt_22(1)%restricted_braket =>  SRT1TT4
!      srt_22(1)%restricted_bra =  SRT1LS1
       srt_22(1)%restricted_bra => SRT36
       srt_22(1)%restricted_ket => SRT1ST4
       srt_22(2)%restricted_ket => SRT26
       srt_22(3)%restricted_ket => SRT20D

! (2,2)->(3,1) type of sorting routines.
       srt_31(0)%unrestricted => SRT19
       srt_31(1)%unrestricted => SRT1S3
       srt_31(0)%restricted_ket => SRT1T3
       srt_31(1)%restricted_ket => SRT22
       srt_31(0)%restricted_braket => SRT1T3
       srt_31(1)%restricted_braket => SRT22
       srt_31(2)%restricted_braket => SRT7
       srt_31(0)%restricted_bra => SRT7
       srt_31(1)%restricted_bra => SRT19
    
!(2,2)->(1,3) type of sorting routines.
       srt_13(0)%unrestricted => SRT1S2
       srt_13(0)%restricted_bra => SRT1T2
       srt_13(0)%restricted_braket => SRT1T2

!------------------------------

     sortType = sort_type       !very dirty trick. will soon amend it.

     array_original => opArray_left 

     arraytype = 'unrestricted'

     p = opArray_left(1)
     q = opArray_left(2)
     r = opArray_left(3)
     s = opArray_left(4)

    if (p(1:1) == q(1:1)) then
      if (r(1:1) == s(1:1)) then
       arraytype = 'restricted_braket'
      else
       arraytype = 'restricted_bra'
      endif 
      elseif (r(1:1)==s(1:1)) then
       arraytype = 'restricted_ket'
    endif
    if_trans   = 'N'
    do_direct  = .TRUE.
    do_exchange= .FALSE.
    no_match   = .TRUE.

    if (do_direct) then

    call align_array(opArray_left, opArray_target, opArray_right, reorder_v)

    !reorder_v=0

    !  k = 1
    !  l=3
    !do i = 1, 4       ! loop over vorder 
    !  do j = 1, dim_w     ! loop over opArray_target 

    !       if (opArray_left(i).eq.opArray_target(j)) then
    !          reorder_v(k) = i
    !          k = k+1
    !          exit
!   !          elseif (opArray_left(i).eq.opArray_right(j)) then !CAUTION. this is not completely general. I have to change it.
!   !          reorder_v(l) = i
!   !          l = l+1
    !       endif
    !  
    !  enddo
    !enddo
    !write(*,*)'value of k',k 
    !do i =1, 4      !loop over vorder
    !  do j = 1, dim_t     ! loop over opArray_right 

    !       if (opArray_left(i).eq.opArray_right(j)) then
    !          reorder_v(k) = i
    !          k = k+1
    !          exit
! ! !          elseif (opArray_left(i).eq.opArray_right(j)) then !CAUTION. this is not completely general. I have to change it.
! ! !          reorder_v(l) = i
! ! !          l = l+1
    !       endif
    !  
    !  enddo
    !enddo

  
    
!!!findloc feature is not available for intel compilers. so I can't use it now.

    a%trans = toDecimal(reorder_v,.true.)    
    a%num = toDecimal(reorder_v)    

   call match_array (a%num, no_match, save_index, if_trans)

   if (no_match) call match_array_transpose (a%trans, no_match, save_index, if_trans)  

   if (no_match) call shuffle_indices (reorder_v, no_match, save_index,if_trans)

!   call shuffle(internal,reorder_v, no_match,save_index)

!   if (.not.no_match) exit

!   call shuffle(external_and_internal,reorder_v, no_match,save_index)

!   if (.not.no_match) exit

     if (no_match) then
      write (*,*) "no matching sorting subroutine available for direct mode of contraction."
      do_exchange = .TRUE.
     endif
  
   endif 

   if (do_exchange) then

    array_exchange = opArray_left

    call swap_character(array_exchange,3,4) !it is possible to impart more flexibility.


    arraytype = 'unrestricted'

     p = opArray_left(1)
     q = opArray_left(2)
     r = opArray_left(3)
     s = opArray_left(4)

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

!      k = 1
!      l = 3
!    do i = 1, 4       ! loop over vorder 
!      do j = 1, 2     ! loop over opArray_target 
!           if (array_exchange(i).eq.opArray_target(j)) then
!              reorder_v(k) = i
!              k = k+1
!              exit
!              elseif (array_exchange(i).eq.opArray_right(j)) then !this is not completely general. I have to change it.
!              reorder_v(l) = i
!              l = l+1
!           endif
!      
!      enddo
!    enddo


     call align_array(array_exchange, opArray_target, opArray_right, reorder_v)

    a%num   = toDecimal(reorder_v)
    a%trans = toDecimal(reorder_v,.true.)    

    call match_array (a%num, no_match, save_index,if_trans)
      
    if (no_match) call match_array_transpose (a%trans, no_match, save_index, if_trans)  

    if (no_match) call shuffle_indices (reorder_v, no_match, save_index, if_trans)

    endif

     if (no_match) then
      write (*,*) "no matching sorting subroutine available. have to do something different...."
     endif

    if (.not.no_match) then

     save_index = save_index - 1

     call integral_fetching_incore(V) 

     if (if_trans == 'C') then
       if (sortType == 31) sort_type_dummy = 13
       if (sortType == 13) sort_type_dummy = 31
     endif

    if (sortType == 22) then

     if (save_index > 0) then
      select case(arraytype)
        
         case('unresticted')

         case('restricted_braket') 

         case('restricted_bra')

        if (present(sorted_opArray)) then 
          if (.NOT.do_exchange) then
            sorted_opArray = get_operator_array(opArray_left, a%num)
          else
            sorted_opArray = get_operator_array(array_exchange, a%num)
          endif
        endif
          call input_index(c,d)

           row = c(d%bra_reorder)%tuple_combInd
           column = c(d%ket_reorder)%tuple_combInd     

          call srt_22(save_index)%restricted_bra(nrep, multb,.true.,.false.,c(d%I)%orbital_type, c(d%J)%orbital_type,&
                  & c(d%K)%orbital_type, c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braKet)%full_array,  &
                  & c(d%ket_reorder)%tuple,c(d%bra_reorder)%tuple,V%presorted,v_sorted)  

         case('restricted_ket')

       end select
     else
      write(*,*)'==sorting is not needed=='
       v_sorted = V%presorted  !dirtiest solution. will remove it soon.  

       call input_index(c,d)

       row = c(d%bra)%triangular%tuple_combInd
       column = c(d%ket)%triangular%tuple_combInd     

!      call move_alloc(V%presorted,v_sorted)
     endif

     elseif (sort_type_dummy == 31) then
      select case(arraytype)
        
         case('unresticted')

         case('restricted_braket')

          call input_index(c,d,m)

          column = m(d%jkl)%qrs%triangular%tuple_combInd
          row = c(d%j)%orbital_type  

          call srt_31(save_index)%restricted_braket(nrep,multb,.false.,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type, &
                     & c(d%L)%orbital_type,m(d%jkl)%qrs%triangular%tuple_combInd, & 
                     & m(d%ijkl)%qrs%full_array,m(d%jkl)%qrs%triangular%tuple,v%presorted,v_sorted,c(d%I)%orbital_type, &
                     & c(d%J)%orbital_type,c(d%ket)%triangular%tuple_combInd)

         case('restricted_bra')

         case('restricted_ket')

      end select

     elseif (sort_type_dummy == 13) then
      select case(arraytype)
        
         case('unresticted')

         case('restricted_braket')

          call input_index(c,d,m)

         case('restricted_bra')

         case('restricted_ket')

      end select
   endif
     deallocate(V%presorted)
 endif
 
  end subroutine


  subroutine decide_sorting_routine(sort_type,opArray_left,opArray_right,opArray_target,if_trans, &
                                &    save_index,sorted_opArray,inverse)

!--------------calling variables---------------
    character(2), intent(in),target        :: opArray_left(:)
    character(2), intent(in)               :: opArray_target(:),opArray_right(:)
    integer,      intent(in)               :: sort_type
    logical,      intent(in), optional     :: inverse  
    character(1), intent(out)              :: if_trans
    integer,      intent(out)              :: save_index 
    character(2), intent(out),optional     :: sorted_opArray(size(opArray_left))
!----------------------------------------------

!-------local variables--------
    character(2) :: array_exchange(size(opArray_left))
    integer      :: reorder_v(size(opArray_left))
    integer      :: i,j,k,l 
    logical      :: no_match, t_transpose  
    character(2) :: p,q,r,s

    type(sorting_list)     :: a
!------------------------------
     sortType = sort_type       !very dirty trick. will soon amend it.

     array_original => opArray_left 

     arraytype = 'unrestricted'

     p = opArray_left(1)
     q = opArray_left(2)
     r = opArray_left(3)
     s = opArray_left(4)

    if (p(1:1) == q(1:1)) then
      if (r(1:1) == s(1:1)) then
       arraytype = 'restricted_braket'
      else
       arraytype = 'restricted_bra'
      endif 
      elseif (r(1:1)==s(1:1)) then
       arraytype = 'restricted_ket'
    endif
    
    do_direct = .true.
    do_exchange = .false.
    no_match = .true.

   if (present(inverse)) then 

    call align_array(array1=opArray_left,array2=opArray_target,reorder_array=reorder_v)

!findloc feature is not available for intel compilers. so I can't use it now.

    a%num = toDecimal(reorder_v)    

    call match_array (a%num,no_match,save_index,if_trans)
   
    endif

    if (do_direct) then

    call align_array(opArray_left,opArray_target,opArray_right,reorder_v)

!findloc feature is not available for intel compilers. so I can't use it now.

    a%trans = toDecimal(reorder_v,.true.)    
    a%num = toDecimal(reorder_v)    

    call match_array (a%num,no_match,save_index,if_trans)

    if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  

    if (no_match) call shuffle_indices (reorder_v,no_match,save_index,if_trans)

     if (no_match) then
      write (*,*) "====no matching sorting subroutine available for direct mode of contraction.====="
      do_exchange = .true.
     endif

    endif 

   if (do_exchange) then

      write(*,*)'=====contraction will be done in exchange mode====='

       array_exchange = opArray_left

       if (arraytype == 'restricted_ket') call swap_character(array_exchange,3,4)   
       if (arraytype == 'restricted_bra') call swap_character(array_exchange,1,2)   

       call align_array(array_exchange,opArray_target,opArray_right,reorder_v)

       a%num   = toDecimal(reorder_v)
       a%trans = toDecimal(reorder_v,.true.)    

       call match_array (a%num,no_match,save_index,if_trans)
      
       if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  

       if (no_match) call shuffle_indices (reorder_v,no_match,save_index,if_trans)


    if (arraytype == 'restricted_braket') then
       k = 1
       i = 1 
     do

       call swap_character(array_exchange,i,i+1)   

       call align_array(array_exchange,opArray_target,opArray_right,reorder_v)

       a%num   = toDecimal(reorder_v)
       a%trans = toDecimal(reorder_v,.true.)    

       call match_array (a%num,no_match,save_index,if_trans)
      
       if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  
    
       if (.not.no_match) exit
       if (k>2)           exit
       i = i + size(opArray_left)/2 
       array_exchange = opArray_left
       k = k+1

     enddo

    endif
   endif

    exchange = array_exchange

    if (no_match) then
      call quit('No Sorting subroutine available. code stops here.')
    endif

    if (present(sorted_opArray)) then 
      if (.not.do_exchange) then
        sorted_opArray = get_operator_array(opArray_left,toDecimal(reorder_v))
      else
        sorted_opArray = get_operator_array(array_exchange,toDecimal(reorder_v))
      endif
    endif

  reorder_left = reorder_v 

  end subroutine

  subroutine sorting_for_blas(save_index,do_inverse,if_trans,sorted,row,column,presorted) 

!--------------calling variables---------------
    real(8), intent(in),optional           :: presorted(:)
    integer, intent(inout)                 :: save_index
    logical, intent(in)                    :: do_inverse  
    character(1), intent(inout)            :: if_trans
    integer, intent(out)                   :: row(:), column(:)
    real(8), intent(out)                   :: sorted(:) !think a bit more carefully to replace it by dynamic allocation. not easy.
!----------------------------------------------

!-------local variables------------------------
!   character(2) :: array_exchange(4)
    integer      :: reorder_v(4)
    integer      :: i,j,k,l 
    logical      :: t_transpose  
    character(2) :: p,q,r,s
    integer      :: sort_type_dummy
#include "symm.inc"
#include "param.inc"

    type(sorting_list)     :: a
    type(middle_indices)   :: c(0:5)  
    type(more_indices)     :: m(0:15)  
    type(bit_processing)   :: d 
    type(presorting)       :: V 

    type(sorting_list) :: srt_22(0:5)
    type(sorting_list) :: srt_31(0:5)
    type(sorting_list) :: srt_13(0:5)
! (2,2)->(2,2) type of sorting routines.         
       srt_22(1)%unrestricted => SRT1SS4   
       srt_22(2)%unrestricted => SRT1L1    
       srt_22(3)%unrestricted => SRT1C1    
       srt_22(4)%unrestricted => SRT16     
       srt_22(5)%unrestricted => SRT1R1

       srt_22(1)%restricted_braket =>  SRT1TT4
!      srt_22(1)%restricted_bra =  SRT1LS1
       srt_22(1)%restricted_bra => SRT36
       srt_22(1)%restricted_ket => SRT1ST4
       srt_22(2)%restricted_ket => SRT26
       srt_22(3)%restricted_ket => SRT20D

! (2,2)->(3,1) type of sorting routines.
       srt_31(0)%unrestricted => SRT19
       srt_31(1)%unrestricted => SRT1S3
       srt_31(0)%restricted_ket => SRT1T3
       srt_31(1)%restricted_ket => SRT22
       srt_31(0)%restricted_braket => SRT1T3
       srt_31(1)%restricted_braket => SRT22
       srt_31(2)%restricted_braket => SRT7
       srt_31(0)%restricted_bra => SRT7
       srt_31(1)%restricted_bra => SRT19
    
!(2,2)->(1,3) type of sorting routines.
       srt_13(0)%unrestricted => SRT1S2
       srt_13(0)%restricted_bra => SRT1T2
       srt_13(0)%restricted_braket => SRT1T2

!------------------------------

     save_index = save_index - 1

     if (.not.present(presorted)) call integral_fetching_incore(V) 
     if (present(presorted)) v%presorted = presorted 

     sort_type_dummy = sortType

     if (if_trans == 'C') then
       if (sortType == 31) sort_type_dummy = 13
       if (sortType == 13) sort_type_dummy = 31
     endif

    if (sortType == 22) then

     if (save_index > 0) then
      select case(arraytype)
        
         case('unresticted')

         case('restricted_braket') 

         case('restricted_bra')

          call input_index(c,d)

           row = c(d%bra_reorder)%tuple_combInd
           column = c(d%ket_reorder)%tuple_combInd     

          call srt_22(save_index)%restricted_bra(nrep, multb,.true.,do_inverse,c(d%I)%orbital_type, c(d%J)%orbital_type,&
                  & c(d%K)%orbital_type, c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braKet)%full_array,  &
                  & c(d%ket_reorder)%tuple,c(d%bra_reorder)%tuple,V%presorted,sorted)  

         case('restricted_ket')

       end select
     else
      write(*,*)'====sorting is not needed===='
     if(present(presorted))then
     call xcopy(ndimt2,presorted,1,sorted,1) 
     else
     sorted = V%presorted  !dirtiest solution. will remove it soon.  
     endif
       call input_index(c,d)

       row = c(d%bra)%triangular%tuple_combInd
       column = c(d%ket)%triangular%tuple_combInd     

!      call move_alloc(V%presorted,sorted)
     endif

     elseif (sort_type_dummy == 31) then
      select case(arraytype)
        
         case('unresticted')

         case('restricted_braket')

          call input_index(c,d,m)

       if (save_index == 2) then

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%qrs%triangular%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%i)%orbital_type  
          column = m(d%jkl)%qrs%triangular%tuple_combInd
       endif

    call srt_31(save_index)%restricted_braket(nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type, &
                & c(d%L)%orbital_type,m(d%jkl)%qrs%triangular%tuple_combInd, & 
                & m(d%ijkl)%qrs%full_array,m(d%jkl)%qrs%triangular%tuple,v%presorted,sorted,c(d%I)%orbital_type, &
                & c(d%J)%orbital_type,c(d%ket)%triangular%tuple_combInd)
       elseif (save_index == 1) then

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%pqr%triangular%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%pqr%triangular%tuple_combInd
       endif

     call srt_31(save_index)%restricted_braket(nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type,&
              & c(d%L)%orbital_type,m(d%ijk_reorder)%pqr%triangular%tuple_combInd,m(d%ijkl)%pqr%full_array, &
              & m(d%ijk_reorder)%qrs%triangular%tuple,v%presorted,sorted,c(d%I)%orbital_type, &
              & c(d%J)%orbital_type,c(d%ket)%triangular%tuple_combInd)

       else

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%pqr%triangular%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%pqr%triangular%tuple_combInd
       endif

     call srt_31(save_index)%restricted_braket(nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type,&
              & c(d%L)%orbital_type,m(d%ijk_reorder)%pqr%triangular%tuple_combInd,m(d%ijkl_reorder)%pqr%full_array, &
              & m(d%ijk_reorder)%pqr%triangular%tuple,v%presorted,sorted,c(d%I)%orbital_type, &
              & c(d%J)%orbital_type,c(d%ket)%triangular%tuple_combInd)

       endif

        case('restricted_bra')

          call input_index(c,d,m)

          column = c(d%l_reorder)%orbital_type
          row = m(d%ijk_reorder)%pqr%triangular%tuple_combInd  

!         CALL SRT19 (NREP,MULTB,LFA,NOOT,NV,NO,NOOOT,KOOOVT,KKOOOT,BUF5,v_sort)

        call srt_31(save_index)%restricted_bra (nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type, &
             & c(d%L)%orbital_type,m(d%ijk_reorder)%pqr%triangular%tuple_combInd,m(d%ijkl_reorder)%pqr%full_array, &
             & m(d%ijk_reorder)%pqr%triangular%tuple,v%presorted,sorted) 

         case('restricted_ket')

      end select

     elseif (sort_type_dummy == 13) then
      select case(arraytype)
        
         case('unresticted')

         case('restricted_braket')

          call input_index(c,d,m)

         case('restricted_bra')

         case('restricted_ket')

      end select
   endif
     deallocate(V%presorted)
 
  end subroutine

  subroutine T_for_blas(sort_type,array_right,opArray_left,array_im,Ttrans,t2,right_sorted,row,column)

!--------------calling variables---------------
  character(2), dimension(4), intent(in),target :: array_right
  character(2), dimension(:), intent(in) :: opArray_left, array_im
  character(1), intent(out)              :: Ttrans
  integer, intent(in)                     :: sort_type
  integer, intent(out)                    :: row(:), column(:)
  real*8, intent(in)                      :: t2(:)
  real*8, intent(out)                     :: right_sorted(:) 
!----------------------------------------------

!-------local variables--------
  character(2), dimension(4) :: array_exchange
! logical      :: t_direct, t_exchange
  integer      :: reorder_t(4)
  integer      :: i,j,k,l,save_index,temp 
  logical      :: no_match, t_transpose  
  character(2) :: p,q,r,s
  integer      :: sort_type_dummy
#include "symm.inc"
#include "param.inc"
  type(sorting_list) :: a

  type(middle_indices)   :: c(0:5)  
  type(more_indices)     :: m(0:15)  
  type(bit_processing)   :: d 

  type(sorting_list) :: srt_22(0:5)
  type(sorting_list) :: srt_31(0:5)
  type(sorting_list) :: srt_13(0:5)

! (2,2)->(2,2) type of sorting routines.         

  srt_22(1)%restricted_braket =>  SRT1TT4

! (2,2)->(3,1) type of sorting routines.

  srt_31(0)%restricted_braket => SRT1T3
  srt_31(1)%restricted_braket => SRT22
  srt_31(2)%restricted_braket => SRT7

! (2,2)->(1,3) type of sorting routines.

   srt_13(0)%restricted_braket => SRT1T2

!----------------------------------------------
   sortType = sort_type

   Ttrans = 'N'

   array_original => array_right 
!   arraytype = 'restricted_braket'

   t_direct = .TRUE.
   t_exchange = .FALSE.
   no_match = .TRUE.

!1) align the array according to the already ordered v-array.

   if (t_direct) then

   call align_array(array_right, opArray_left, array_im, reorder_t)

   a%num = toDecimal(reorder_t)    
   a%trans = toDecimal(reorder_t,.true.)    

!2) check whether it is possible to match with available sorting routine. both original and transpose array will work in this purpose. we can shuffle the contractable indices as well the uncontractable ones. 

   call match_array (a%num, no_match, save_index, Ttrans)
   if (no_match) call match_array_transpose (a%trans, no_match, save_index, Ttrans) 
   if (no_match) call shuffle_indices (reorder_t, no_match, save_index,Ttrans)

!3) if we don't find a suitable subroutine, we will generate all the exchange partners of original array and try to match.  

     if (no_match) then
      write (*,*) "===no matching sorting subroutine available for direct mode of contraction==="
      t_exchange = .TRUE.
     endif
   endif

   if (t_exchange) then

    array_exchange = array_right

    call swap_character(array_exchange,3,4)

    call align_array(array_exchange,opArray_left,array_im,reorder_t)

    a%num   = toDecimal(reorder_t)
    a%trans = toDecimal(reorder_t,.true.)    

   call match_array (a%num, no_match, save_index, Ttrans)
   if (no_match) call match_array_transpose (a%trans, no_match, save_index, Ttrans) 
   if (no_match) call shuffle_indices (reorder_t, no_match, save_index,Ttrans)

   endif

   if (no_match) then
      write (*,*) "=====no matching sorting subroutine available. have to do something different====="
   endif


!4) call the identified sorting routine.  

   if (.not.no_match) then

     save_index = save_index - 1

      sort_type_dummy = sortType

     if (Ttrans == 'C') then
       if (sortType == 31) sort_type_dummy = 13
       if (sortType == 13) sort_type_dummy = 31
     endif

    if (sortType == 22) then

     if (save_index > 0) then

     else
      write(*,*)'sorting is not needed'
     call xcopy(ndimt2,t2,1,right_sorted,1)  !dirtiest solution. will remove it soon.  
    
     call input_index(a=c,b=d,reorder_v=reorder_t)

     row    = c(d%bra)%triangular%tuple_combInd
     column = c(d%ket)%triangular%tuple_combInd


!     call move_alloc(t2,t_sorted)
    endif

     elseif (sort_type_dummy == 31) then
             
     call input_index(a=c,b=d,k=m,reorder_v=reorder_t)

     row = m(d%ijk)%pqr%triangular%tuple_combInd
     column = c(d%k)%orbital_type  

!    call srt_31(save_index)%restricted_braket(nrep,multb,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type, &
!               & c(d%L)%orbital_type,m(d%ijk)%pqr%triangular%tuple_combInd,m(d%ijkl)%pqr%full_array, &
!               & m(d%ijk)%qrs%triangular%tuple,t2,right_sorted)

     call srt_31(save_index)%restricted_braket(nrep,multb,.false.,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type, &
                & c(d%L)%orbital_type, m(d%ijk)%pqr%triangular%tuple_combInd,m(d%ijkl)%pqr%full_array, &           
                & m(d%ijk)%qrs%triangular%tuple,t2,right_sorted,c(d%I)%orbital_type, &
                & c(d%J)%orbital_type,c(d%ket)%triangular%tuple_combInd)

! at this point there is a conflict of indices depending on whether the subroutine stems from I>J or K>L case. right now the implementation is for K>L case. I have to extend the case for I>J case too. easiest way to do this is to use a if statement.   

     elseif (sort_type_dummy == 13) then

    endif
   endif
 end subroutine

  subroutine onebody_for_blas (array1,array2,array3,trans_1b,row,column,FreeIndx) 

!------------------------------
   character(*),dimension(:),intent(in)  :: array1, array2, array3 ! follows the same instruction as it is in the align_array.
   character(*),intent(out)        :: trans_1b
   integer :: reorder_array(size(array1))      
   integer,intent(out)             :: row(:), column(:) 
   type(middle_indices)             :: a(0:2)  
   integer,dimension(2) :: bit_original
   character(2),optional:: FreeIndx
   integer              :: i
!------------------------------
#include "symm.inc"

  call align_array(array1,array2,array3,reorder_array)   

  if ((toDecimal(reorder_array)) == 21) then
  trans_1b = 'T'
  else
  trans_1b = 'N'
  endif
  
   a(0)%orbital_type => NO
   a(1)%orbital_type => NV

   bit_original = bitString (array1, (/1,2/)) 

   row = a(bit_original(1))%orbital_type
   column = a(bit_original(2))%orbital_type

  if (present(FreeIndx)) then

   FreeIndx = array1(reorder_array(2))

  endif 

  end subroutine 


  subroutine swap_character(array,k,l)

    character(*),intent(inout) :: array(:)
    integer, intent(in)        :: k,l  
    character(2)               :: temp

    temp = array(k)

    array(k) = array(l)
    array(l) = temp

  end subroutine

  subroutine align_array(array1, array2, array3, reorder_array)

!-------------------------------------------------------------
    character(*),intent(in) :: array1(:), array2(:) ! array1 is the array to be reordered
    character(*),intent(in),optional :: array3(:) 
    integer, dimension(size(array1)),intent(out) :: reorder_array     ! array2 is the array from where we will extract the first indices in the reordered 
                                                                       ! array i.e, for left array it is always intermediate and for right always    
    integer                                  :: k,i                    ! left array. array3 is the array to extract rest of the indices e.g, for left array
                                                                       ! always right and for right array always intermediate.
!-------------------------------------------------------------
     k = 1
     do i = 1, size(array2)       ! loop over vorder 
        if (findloc(array1,array2(i)) > 0) then 
          reorder_array(k) = findloc(array1, array2(i))
          k=k+1
        endif
     enddo

    if (present(array3)) then 
     do i = 1, size(array3)       ! loop over vorder 
        if (findloc(array1,array3(i)) > 0) then 
          reorder_array(k) = findloc(array1,array3(i))
          k=k+1
        endif
     enddo
    endif
    
  end subroutine


! subroutine integral_fetching_out_of_core()

! bit_original = bitString (array_original, (/1,2,3,4/)) 

!  value = 10
!  if (popcnt(bit_original)==2) then
!   value = ior(iand(bit_original(1),bit_original(2)),iand(bit_original(3),bit_original(4))) 
!  endif



!  select case(integrals) 

!   case (popcnt(bit_original)==0)
!     allocate(Voooo(NV1))
!     call getOOOO(Voooo)
!   case (popcnt(bit_original)==1)
!     allocate(Vvooo(NV2))
!     call getVOOO(Vvooo)
!   case (value==1)
!     allocate(VVVOO(NV3))
!     call getVVOO(VVVOO) 
!   case (value==0)
!     allocate(VVOVO(NV4))
!     call getVOVO(VVOVO)

!  end select

! end subroutine


  subroutine integral_fetching_incore(G)
  implicit none
!-------------------------------------------------------------
  integer :: value
  integer,dimension(4):: bit_original, bit_available
  real*8, allocatable :: VOOOO(:),VVOOO(:),VVOVO(:),VVVOO(:) 
  integer :: AllocateStatus 
  type(bit_processing):: b,c 
  type(presorting), intent(out) :: G
!-------------------------------------------------------------
#include "symm.inc"

    bit_original = bitString (array_original, (/1,2,3,4/)) 

    b%bra = ior(bit_original(1),bit_original(2)) + iand(bit_original(1),bit_original(2))  
    b%ket = ior(bit_original(3),bit_original(4)) + iand(bit_original(3),bit_original(4)) 
    value = popcnt(bin2dec(bit_original))+iand(ior(bit_original(1),bit_original(2)),ior(bit_original(3),bit_original(4))) 

!here you can think of keeping all the VVVO integrals in fast memory. 

   select case(value) 

    case (0)
      allocate(Voooo(NV1))
      call getOOOO(Voooo)
    case (1)
      allocate(Vvooo(NV2))
      call getVOOO(Vvooo)

      bit_available = bitString((/"p","o","o","o"/),(/1,2,3,4/))
     
      c%bra = ior(bit_available(1),bit_available(2)) + iand(bit_available(1),bit_available(2))  
      c%ket = ior(bit_available(3),bit_available(4)) + iand(bit_available(3),bit_available(4)) 
      if (c%bra == b%ket)then

         allocate(G%presorted(NV2))

       call SRT1C1 (NREP,NVO,NOOT,VVOOO,G%presorted) 

       deallocate(VVOOO) 

      endif
!      if (bin2dec((/bit_available(2),bit_available(1)/))== bin2dec((/bit_original(1),bit_original(2)/))) call
    case (2)

      ALLOCATE (VVVOO(NV3),STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

      call getVVOO(VVVOO) 

      bit_available = bitString((/"p","p","o","o"/),(/1,2,3,4/))
     
      c%bra = ior(bit_available(1),bit_available(2)) + iand(bit_available(1),bit_available(2))  
      c%ket = ior(bit_available(3),bit_available(4)) + iand(bit_available(3),bit_available(4)) 

      if (c%bra == b%ket)then

      allocate(G%presorted(NV3),STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      G%presorted = 0.0d0 

      call SRT1C1 (NREP,NVVT,NOOT,VVVOO,G%presorted) 

      deallocate(VVVOO)
 
      endif
    case (3)
      allocate(VVOVO(NV4))
 !     allocate(G%presorted(NV4))
      call getVOVO(VVOVO)

!      v_presorted => VVOVO
    case default
      write(*,*)'no suitable array to fetch'
   end select

  if (allocated(VOOOO)) call move_alloc(VOOOO,G%presorted)
  if (allocated(VVOOO)) call move_alloc(VVOOO,G%presorted)
  if (allocated(VVOVO)) call move_alloc(VVOVO,G%presorted)
  if (allocated(VVVOO)) call move_alloc(VVVOO,G%presorted)

  end subroutine

  subroutine input_index(a,b,k,reorder_v)

!-------------------------------
 integer,dimension(:),allocatable,target :: bit_reordered, bit_original
 integer,dimension(4),intent(in),optional:: reorder_v 
! logical, intent(in)                :: do_exchange
! character(2), intent(in)           :: array_exchange(4)
 type(middle_indices),intent(out)     :: a(0:5)  
 type(bit_processing),intent(out)   :: b 
 type(more_indices),intent(out),optional :: k(0:15)
 integer :: i
 integer, pointer :: bitpqr(:),bitqrs(:),bitpqr_reorder(:),bitqrs_reorder(:)
!-------------------------------
#include "symm.inc"
     a(0)%tuple => JJOO         
     a(1)%tuple => JJVO
     a(2)%tuple => JJVV
     a(0)%tuple_combInd => MOO
     a(1)%tuple_combInd => MVO
     a(2)%tuple_combInd => MVV

     a(0)%triangular%tuple_combInd => NOOT
     a(2)%triangular%tuple_combInd => NVVT

     a(0)%full_array => JOOOO 
     a(1)%full_array => JVOOO
     a(2)%full_array => JVOVO
     a(3)%full_array => JVVOO  
!    a(4)%full_array => JVVVO
!    a(5)%full_array => 
     a(0)%orbital_type => NO
     a(1)%orbital_type => NV

  if (present(k)) then

     k(0)%pqr%triangular%tuple => KKOOOT
     k(1)%pqr%triangular%tuple => KKOOVT
     k(2)%pqr%triangular%tuple => KKVVOT
     k(3)%pqr%triangular%tuple => KKVVVT
 
     k(0)%pqr%triangular%tuple_combInd => NOOOT
     k(1)%pqr%triangular%tuple_combInd => NOOVT
     k(2)%pqr%triangular%tuple_combInd => NVVOT
     k(3)%pqr%triangular%tuple_combInd => NVVVT

     k(0)%qrs%triangular%tuple => LLOOOT
     k(1)%qrs%triangular%tuple => LLVOOT
     k(2)%qrs%triangular%tuple => LLOVVT
     k(3)%qrs%triangular%tuple => LLVVVT

!    k(0)%qrs%triangular%tuple_combInd => NOOOT2
     k(1)%qrs%triangular%tuple_combInd => NVOOT
     k(2)%qrs%triangular%tuple_combInd => NOVVT
!    k(3)%qrs%triangular%tuple_combInd => NVVVT2

     k(0)%pqr%full_array => KOOOOT
     k(1)%pqr%full_array => KOOOVT
     k(2)%pqr%full_array => KOOVOT
     k(3)%pqr%full_array => KOOVVT
     k(12)%pqr%full_array=> KVVOOT
!    k(14)%pqr%full_array=> KVVVOT  !this is commented because KVVVOT is int*8 type. plan is to cast it as int*8 from int*4 

     k(3)%qrs%full_array => LOOVVT

  endif

  if (present(reorder_v)) reorder_left = reorder_v  

  allocate (bit_reordered(size(array_original)))
  allocate (bit_original(size(array_original)))

  if (do_exchange) then
  bit_reordered = bitString (exchange,reorder_left)
  bit_original  = bitString (exchange, (/1,2,3,4/)) 
  else
  bit_reordered  = bitString (array_original,reorder_left)
  bit_original  = bitString (array_original, (/1,2,3,4/)) 
  endif

! unique indexing for pair indices.  

    b%bra_reorder = ior(bit_reordered(1),bit_reordered(2)) + iand(bit_reordered(1),bit_reordered(2))  
    b%ket_reorder = ior(bit_reordered(3),bit_reordered(4)) + iand(bit_reordered(3),bit_reordered(4)) 
    b%braKet = b%bra_reorder + b%ket_reorder + &
            &  ior(iand(bit_reordered(1),bit_reordered(2)),iand(bit_reordered(3),bit_reordered(4)))  

    b%bra = ior(bit_original(1),bit_original(2)) + iand(bit_original(1),bit_original(2))  
    b%ket = ior(bit_original(3),bit_original(4)) + iand(bit_original(3),bit_original(4)) 

! unique indexing for combined triple indices.

    bitpqr=>bit_original(1:3)
    bitqrs=>bit_original(2:4)

    b%ijk = popcnt(bin2dec(bitpqr)) 
    b%jkl = popcnt(bin2dec(bitqrs))

    bitpqr_reorder=>bit_reordered(1:3)  
    bitqrs_reorder=>bit_reordered(2:4)  

    b%ijk_reorder = popcnt(bin2dec(bitpqr_reorder)) 
    b%jkl_reorder = popcnt(bin2dec(bitqrs_reorder))

! unique indexing for individual label of an array.

    b%I = bit_original(1) 
    b%J = bit_original(2)
    b%K = bit_original(3)
    b%L = bit_original(4)

    b%I_reorder = bit_reordered(1) 
    b%J_reorder = bit_reordered(2)
    b%K_reorder = bit_reordered(3)
    b%L_reorder = bit_reordered(4)
! unique indexing for full array when triangular index tuple are being used.
    b%ijkl_reorder = bin2dec(bit_reordered)  
    b%ijkl = bin2dec(bit_original)  

  deallocate (bit_reordered)
  deallocate (bit_original)
  end subroutine

  integer function findloc(array,value)

  character(*), dimension(:),intent(in) :: array
  character(2),intent(in) :: value
  integer      :: i  

     do i = 1, size(array)       ! loop over vorder 
        if (array(i)==value) then
          findloc = i
          exit
          else
          findloc = 0
        endif
     enddo
   end function

  function bitString (generic_array, reorder_index_array)

  character(*), dimension(:),intent(in) :: generic_array
  integer,      dimension(:),intent(in) :: reorder_index_array
  integer      :: i,k
  integer, dimension(size(generic_array)) ::bitString 
  character(2), allocatable :: temp_array(:)

  allocate(temp_array(size(generic_array)))     

     do i = 1, size(generic_array)
        temp_array(i)=generic_array(reorder_index_array(i))
        k = ichar(temp_array(i)(1:1))-ichar('o')
        bitString(i)=k
     enddo

  deallocate(temp_array)

  end function   

  integer function toDecimal(index_array, trans)
   
  integer, dimension(:),intent(in),target :: index_array
  logical, optional, intent(in) :: trans
  integer :: scratch(size(index_array))
  integer,pointer :: scratch_bra(:)
  integer,pointer :: scratch_ket(:)
  logical                       :: do_trans 
  integer                       :: temp, i
  
  scratch = index_array
  if (present(trans)) then
   do_trans = trans
  else
   do_trans = .false.
  endif
 
  if (do_trans) then

! it should not be case dependent. if possible try to modify it.

  select case (sortType)
     case (22)
     scratch_bra => index_array(3:4) 
     scratch_ket => index_array(1:2) 
     case (31)
     scratch_bra => index_array(4:4)
     scratch_ket => index_array(1:3)
     case (13)
     scratch_bra => index_array(2:4)
     scratch_ket => index_array(1:1)
  end select  

     scratch = (/scratch_bra,scratch_ket/)  
  
  endif   

    toDecimal=0

    do i = size(index_array), 1 , -1  
      toDecimal = toDecimal + scratch(i)*10**(size(index_array)-i)
    enddo 

  end function

  integer function bin2dec(bitstring)

    integer, intent(in) :: bitstring(:)
    integer             :: i 

    bin2dec=0
 
    do i = 1, size(bitstring)
     bin2dec = bin2dec  + bitstring(i)*2**(size(bitstring)-i)
    enddo

  end function

  recursive subroutine permutate(E, P) 
    integer, intent(in)  :: E(:)       ! array of objects 
    integer, intent(out) :: P(:,:)     ! permutations of E 
    integer  :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1) 
       N = size(E); Nfac = size(P,1); 
     do i=1,N                           ! cases with E(i) in front 
          if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S) 
         forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/) 
     end do 
  end subroutine permutate 

  function get_operator_array(character_array,integer_for_array)
  integer,intent(in) :: integer_for_array
  character(*),intent(in) :: character_array(:)  
  character(2)             :: get_operator_array(size(character_array))  
  integer :: y, temp
  integer :: x, i

  temp = integer_for_array
  i = 1

  do  
     y = temp/10**(size(character_array)-i)  
     x = mod(temp, 10**(size(character_array)-i))
     temp = x
     get_operator_array(i) = character_array(y)
     i = i+1
     if (i > size(character_array)) exit
  enddo 
  end function

  character(len=10) function seqstring(seqnum) 
    integer, intent(in) :: seqnum
    character(10)       :: temp
    write (temp,'(I0)') seqnum
    seqstring = trim(temp)
  end function

  end module   
