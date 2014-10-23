  module wrapper_sorting

    implicit none

    public shuffle_indices
    public match_array
    public v_for_blas

    private

    integer    :: sortType
    logical      :: v_direct
    character(2),pointer :: array_original(:)

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
       integer :: bra, ket, braKet, I, J, K, L, ijk, jkl, ijkl, bra_reorder, ket_reorder
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
    
   subroutine shuffle_indices(arraytype, shuffle, reorder_array, no_match, save_index, vtrans)

! ------- variables --------
    integer, intent(inout),target   :: reorder_array(4)
    logical, intent(inout)   :: no_match
    character(*), intent(in) :: shuffle, arraytype
    character(*), intent(out):: vtrans
    integer,pointer :: shifted_array(:)
    integer                  :: temp  
    integer, intent(out)     :: save_index
    type(sorting_list)       :: a
! --------------------------

  select case(shuffle)
    
    case('external')

     shifted_array => reorder_array(2:4)

     temp = reorder_array(1)
!     reorder_array(1) = reorder_array(2)
!     reorder_array(2) = temp

!    a%num = toDecimal(reorder_array)
!    a%trans = toDecimal(reorder_array,.true.)    

    a%num = toDecimal((/temp,cshift(shifted_array,shift=-1)/))
    a%trans = toDecimal((/temp,cshift(shifted_array,shift=-1)/),.true.)    

   if (sortType ==22) then

     temp = reorder_array(1)
     reorder_array(1) = reorder_array(2)
     reorder_array(2) = temp

    a%num = toDecimal(reorder_array)
    a%trans = toDecimal(reorder_array,.true.)    

   endif

    case('internal')
    
     temp = reorder_array(3)
     reorder_array(3) = reorder_array(4)
     reorder_array(4) = temp
   
    a%num = toDecimal(reorder_array)
    a%trans = toDecimal(reorder_array,.true.)    
    case('external_and_internal') 

     temp = reorder_array(1)
     reorder_array(1) = reorder_array(2)
     reorder_array(2) = temp

     temp = reorder_array(3)
     reorder_array(3) = reorder_array(4)
     reorder_array(4) = temp

    a%num = toDecimal(reorder_array)
    a%trans = toDecimal(reorder_array,.true.)    
  end select
   
   call match_array (arraytype, a%num, no_match, save_index, vtrans)  
   if (no_match) call match_array_transpose (arraytype, a%trans, no_match, save_index, vtrans)  

   end subroutine

  subroutine match_array (arraytype, permute, no_match, save_index, vtrans) 

!-----------variable-------------
   integer, intent(in) :: permute 
   logical, intent(inout):: no_match
   integer, intent(out):: save_index
   character(*), intent(out) :: vtrans
   character(*), intent(in)  :: arraytype   
!--------------------------------
!---------local variable---------
   integer, allocatable :: A(:)
   integer :: i, dim_A, upbound_A
!--------------------------------

   vtrans='N'

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
  allocate(A(1))
  A = (/1342/)

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

  subroutine match_array_transpose (arraytype, permute_transpose, no_match, save_index, vtrans) 

!-----------variable-------------
   integer, intent(in) :: permute_transpose
   logical, intent(inout):: no_match
   integer, intent(out):: save_index
   character(*), intent(out) :: vtrans
   character(*), intent(in)  :: arraytype   
!--------------------------------
!---------local variable---------
  integer, allocatable :: A(:)
  integer :: i, dim_A, upbound_A
  integer :: sort_type
!--------------------------------

   vtrans='N'

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
    vtrans ='C'  
    exit
    endif
     i = i+1
    if (i > upbound_A) exit
   enddo 

  deallocate(A)

  end subroutine

  subroutine v_for_blas (sort_type, array_left, torder, worder, vtrans, v_sorted, row, column)

!--------------calling variables---------------
    character (*), dimension(4), intent(in),target :: array_left
    character (*), dimension(:), intent(in) :: worder, torder
    character (*), intent(out)              :: vtrans
    integer, intent(in)                     :: sort_type
    integer, intent(out)                    :: row(:), column(:)
    real*8, intent(out)                     :: v_sorted(:) !think a bit more carefully to replace it by dynamic allocation. not easy.
!----------------------------------------------

!-------local variables--------
    character(2), dimension(4) :: array_exchange
    logical      :: v_direct, v_exchange
    integer      :: reorder_v(4)
    integer      :: i,j,k,l,dim_t,dim_w,save_index,temp 
    logical      :: no_match, t_transpose  
    character(2) :: p,q,r,s
    character(20):: arraytype
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
       srt_22(1)%unrestricted => SRT1SS4   
       srt_22(2)%unrestricted => SRT1L1    
       srt_22(3)%unrestricted => SRT1C1    
       srt_22(4)%unrestricted => SRT16     
       srt_22(5)%unrestricted => SRT1R1

       srt_22(1)%restricted_braket =>  SRT1TT4
!      srt_22(1)%restricted(2) =  SRT1LS1
       srt_22(1)%restricted_bra => SRT36
       srt_22(2)%restricted_bra => SRT1TS4
       srt_22(1)%restricted_ket => SRT22
       srt_22(2)%restricted_ket => SRT26
       srt_22(3)%restricted_ket => SRT20D

! next sorting routines are (2,2)-> (3,1) type
       srt_31(0)%unrestricted => SRT19
       srt_31(1)%unrestricted => SRT1S3
       srt_31(0)%restricted_ket => SRT1T3
       srt_31(1)%restricted_ket => SRT22
       srt_31(0)%restricted_braket => SRT1T3
       srt_31(1)%restricted_braket => SRT22
       srt_31(2)%restricted_braket => SRT7
       srt_31(0)%restricted_bra => SRT7
    
! next sorting routines are (2,2)-> (1,3) type
       srt_13(0)%unrestricted => SRT1S2
       srt_13(0)%restricted_bra => SRT1T2
       srt_13(0)%restricted_braket => SRT1T2

!------------------------------

     sortType = sort_type       !very dirty trick. will soon amend it.

     array_original => array_left 

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
    dim_w = size(worder)
    
    v_direct = .TRUE.
    v_exchange = .FALSE.
    no_match = .TRUE.

    if (v_direct) then

    call align_array(array_left, worder, torder, reorder_v)


    !reorder_v=0

    !  k = 1
    !  l=3
    !do i = 1, 4       ! loop over vorder 
    !  do j = 1, dim_w     ! loop over worder 

    !       if (array_left(i).eq.worder(j)) then
    !          reorder_v(k) = i
    !          k = k+1
    !          exit
!   !          elseif (array_left(i).eq.torder(j)) then !CAUTION. this is not completely general. I have to change it.
!   !          reorder_v(l) = i
!   !          l = l+1
    !       endif
    !  
    !  enddo
    !enddo
    !write(*,*)'value of k',k 
    !do i =1, 4      !loop over vorder
    !  do j = 1, dim_t     ! loop over torder 

    !       if (array_left(i).eq.torder(j)) then
    !          reorder_v(k) = i
    !          k = k+1
    !          exit
! ! !          elseif (array_left(i).eq.torder(j)) then !CAUTION. this is not completely general. I have to change it.
! ! !          reorder_v(l) = i
! ! !          l = l+1
    !       endif
    !  
    !  enddo
    !enddo


    
!!!findloc feature is not available for intel compilers. so I can't use it now.

    a%trans = toDecimal(reorder_v,.true.)    
    a%num = toDecimal(reorder_v)    

   call match_array (arraytype, a%num, no_match, save_index, vtrans)

   if (no_match) call match_array_transpose (arraytype, a%trans, no_match, save_index, vtrans)  

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

    call swap_character(array_exchange,3,4) !it is possible to impart more flexibility.

    write(*,*)'array exchange',(array_exchange(i),i=1,4) 

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

!      k = 1
!      l = 3
!    do i = 1, 4       ! loop over vorder 
!      do j = 1, 2     ! loop over worder 
!           if (array_exchange(i).eq.worder(j)) then
!              reorder_v(k) = i
!              k = k+1
!              exit
!              elseif (array_exchange(i).eq.torder(j)) then !this is not completely general. I have to change it.
!              reorder_v(l) = i
!              l = l+1
!           endif
!      
!      enddo
!    enddo


     call align_array(array_exchange, worder, torder, reorder_v)

    a%num   = toDecimal(reorder_v)
    a%trans = toDecimal(reorder_v,.true.)    

    call match_array (arraytype,a%num, no_match, save_index,vtrans)
      
    if (no_match) call match_array_transpose (arraytype, a%trans, no_match, save_index, vtrans)  

    if (no_match) call shuffle_indices (arraytype,'external', reorder_v, no_match, save_index, vtrans)

    endif

     if (no_match) then
      write (*,*) "no matching sorting subroutine available. have to do something different...."
     endif

    if (.not.no_match) then

     save_index = save_index - 1

     call integral_fetching_incore(V) 

     if (vtrans == 'C') then
       if (sortType == 31) sort_type_dummy = 13
       if (sortType == 13) sort_type_dummy = 31
     endif

    if (sortType == 22) then

     if (save_index > 0) then
      select case(arraytype)
        
         case('unresticted')

         case('restricted_braket') 

         case('restricted_bra')

          call input_index(array_exchange,v_exchange,reorder_v,c,d)

           row = c(d%bra_reorder)%tuple_combInd
           column = c(d%ket_reorder)%tuple_combInd     

          call srt_22(save_index)%restricted_bra(nrep, multb, ltr, lfa, c(d%I)%orbital_type, c(d%J)%orbital_type,&
                  & c(d%K)%orbital_type, c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braKet)%full_array,  &
                  & c(d%ket_reorder)%tuple,c(d%bra_reorder)%tuple,V%presorted,v_sorted)  
    
         write(*,*)'successfully called the subroutine'

         case('restricted_ket')

       end select
     else
      write(*,*)'sorting is not needed'
       v_sorted = V%presorted  !dirtiest solution. will remove it soon.  
!      call move_alloc(V%presorted,v_sorted)
     endif

     elseif (sort_type_dummy == 31) then
      select case(arraytype)
        
         case('unresticted')

         case('restricted_braket')

          call input_index(array_exchange,v_exchange,reorder_v,c,d,m)

          column = m(d%jkl)%qrs%triangular%tuple_combInd
          row = c(d%j)%orbital_type  

          call srt_31(save_index)%restricted_braket(nrep,multb,c(d%bra)%triangular%tuple_combInd,c(d%I)%orbital_type, &
                     & c(d%J)%orbital_type,c(d%ket)%triangular%tuple_combInd,m(d%jkl)%qrs%triangular%tuple_combInd, & 
                     & m(d%ijkl)%qrs%full_array,m(d%jkl)%qrs%triangular%tuple,v%presorted,v_sorted)

         case('restricted_bra')

         case('restricted_ket')

      end select

     elseif (sort_type_dummy == 13) then
      select case(arraytype)
        
         case('unresticted')

         case('restricted_braket')

          call input_index(array_exchange,v_exchange,reorder_v,c,d,m)

         case('restricted_bra')

         case('restricted_ket')

      end select
   endif
     deallocate(V%presorted)
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
    character (*), dimension(:), intent(in)  :: array1, array2, array3 ! array1 is the array to be reordered
    integer, dimension(size(array1)), intent(out) :: reorder_array     ! array2 is the array from where we will extract the first indices in the reordered 
                                                                       ! array i.e, for left array it is always intermediate and for right always    
    integer                                  :: k,i                    ! left array. array3 is the array to extract rest of the indices e.g, for left array
                                                                       ! always right and for right array always intermediate.
!-------------------------------------------------------------
     k = 1
     do i = 1, size(array2)       ! loop over vorder 
        if (findloc(array1, array2(i)) > 0) then 
          reorder_array(k) = findloc(array1, array2(i))
          k=k+1
        endif
     enddo

     do i = 1, size(array3)       ! loop over vorder 
        if (findloc(array1, array3(i)) > 0) then 
          reorder_array(k) = findloc(array1, array3(i))
          k=k+1
        endif
     enddo

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
  integer, dimension(4) :: bit_original, bit_available
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

  subroutine input_index(array_exchange,v_exchange,reorder_v,a,b,k)

!-------------------------------
 integer, dimension(:), allocatable,target :: bit_reordered, bit_original
 integer, dimension(4),intent(in)   :: reorder_v 
 logical, intent(in)                :: v_exchange
 character(2), intent(in)           :: array_exchange(4)
 type(middle_indices),intent(out)     :: a(0:5)  
 type(bit_processing),intent(out)   :: b 
 type(more_indices),intent(out),optional :: k(0:15)
 integer :: i
 integer, pointer :: bitpqr(:), bitqrs(:)
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

  allocate (bit_reordered(size(array_original)))
  allocate (bit_original(size(array_original)))

  if (v_exchange) then
  bit_reordered = bitString (array_exchange, reorder_v)
  bit_original  = bitString (array_exchange, (/1,2,3,4/)) 
  else
  bit_reordered  = bitString (array_original, reorder_v)
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

! unique indexing for individual label of an array.

    b%I = bit_original(1) 
    b%J = bit_original(2)
    b%K = bit_original(3)
    b%L = bit_original(4)

! unique indexing for full array when triangular index tuple are being used.
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

  end module   
