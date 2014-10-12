  module wrapper_sorting

    implicit none

    public sort_22_unrestricted
!    public sort_13
!    public sort_31   
    public shuffle_indices
    public match_array
    public v_for_blas


    private
  
!    character(2) :: array_exchange(4)
    logical      :: v_direct
!    integer      :: reorder_v(4)
    character(2),pointer :: array_original(:)

    type sorting_list
      integer :: num 
      integer :: trans
!     contains
      procedure(),pointer, nopass :: unrestricted !, restricted 
      procedure(),pointer, nopass :: restricted_braket, restricted_bra, restricted_ket 
    endtype sorting_list 

    type indices
      integer,dimension(:,:),pointer :: pair
      integer, dimension(:),pointer  :: pair_combInd, orbital_type, full_array  
    endtype indices

    type bit_processing
       integer :: bra, ket, braKet, I, J, K, L
    endtype bit_processing 

    type presorting
       real*8, allocatable :: presorted(:) 
    endtype
<<<<<<< HEAD
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

!! abstract interface
!! 
!!   subroutine SRT1SS4
!!   end subroutine SRT1SS4

!!   subroutine SRT1L1
!!   end subroutine SRT1L1

!!    subroutine SRT1C1
!!    end subroutine SRT1C1

!!   subroutine SRT16
!!   end subroutine SRT16

!!   subroutine SRT1R1
!!   end subroutine SRT1R1

!!   subroutine SRT1TT4
!!   end subroutine SRT1TT4

!!   subroutine SRT36 
!!   end subroutine SRT36 

!!   subroutine SRT22
!!   end subroutine SRT22

!!   subroutine SRT26
!!   end subroutine SRT26

!!   subroutine SRT20D
!!   end subroutine SRT20D

!!  end interface  
=======


   abstract interface
   
     subroutine SRT1SS4
     end subroutine SRT1SS4

     subroutine SRT1L1
     end subroutine SRT1L1

!     subroutine SRT1C1
!     end subroutine SRT1C1

     subroutine SRT16
     end subroutine SRT16

     subroutine SRT1R1
     end subroutine SRT1R1

     subroutine SRT1TT4
     end subroutine SRT1TT4

     subroutine SRT36 
     end subroutine SRT36 

     subroutine SRT22
     end subroutine SRT22

     subroutine SRT26
     end subroutine SRT26

     subroutine SRT20D
     end subroutine SRT20D

    end interface  
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0

    CONTAINS
    

    subroutine sort_22_unrestricted(vorder,torder, worder,vtrans)

!----------calling variables-----------------------------------
    character (2), dimension(4), intent(in) :: vorder
    character (*), dimension(:), intent(in) :: worder, torder
    character (*), intent(out) :: vtrans
!--------------------------------------------------------------

! ------- local variables --------
    integer      :: reorder_v(4)
    integer      :: i,j,k,l,dim_t,save_index  
    logical      :: t_transpose
    real*8       :: v_final(10000) 

! --------------------------------
   

!   call integral_fetching_incore(V) 
<<<<<<< HEAD
  ! call v_for_blas (vorder, worder, torder, vtrans, v_final)
=======
   call v_for_blas (vorder, worder, torder, vtrans, v_final)
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0

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
   
    a%num = num(reorder_array)
    a%trans = num(reorder_array,.true.)    

    case('internal')
    
     temp = reorder_array(3)
     reorder_array(3) = reorder_array(4)
     reorder_array(4) = temp
   
    a%num = num(reorder_array)
    a%trans = num(reorder_array,.true.)    
    case('external_and_internal') 

     temp = reorder_array(1)
     reorder_array(1) = reorder_array(2)
     reorder_array(2) = temp

     temp = reorder_array(3)
     reorder_array(3) = reorder_array(4)
     reorder_array(4) = temp

    a%num = num(reorder_array)
    a%trans = num(reorder_array,.true.)    
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
  integer :: i, dim_A, upbound_A
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
    if (i > upbound_A) exit
   enddo 
   endif

  deallocate(A)

  end subroutine

<<<<<<< HEAD
  subroutine v_for_blas (array_left, torder, worder, vtrans, v_sorted, row, column)

!--------------calling variables---------------
    character (*), dimension(4), intent(in),target :: array_left
=======
  subroutine v_for_blas (array_left, worder, torder, vtrans, v_sorted)

!--------------calling variables---------------
    character (2), dimension(4), intent(in),target :: array_left
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0
    character (*), dimension(:), intent(in) :: worder, torder
!    real*8, dimension(:), intent(in) :: v_presorted
    character (*), intent(out)              :: vtrans
    integer, intent(out)                    :: row(:), column(:)
!----------------------------------------------

!-------local variables--------
    character(2), dimension(4) :: array_exchange
!    real*8, allocatable, dimension(:)     :: v_sorted
    real*8, intent(out)       :: v_sorted(*) !think a bit more carefully to replace it by dynamic allocation. not easy.
    logical      :: v_direct, v_exchange
    integer      :: reorder_v(4)
    integer      :: i,j,k,l,dim_t,dim_w,save_index,temp 
    logical      :: no_match, t_transpose  
    character(2) :: p,q,r,s
    character(20):: arraytype
    procedure(swap_character),pointer:: xxx 

#include "symm.inc"
#include "param.inc"

    type(sorting_list) :: a
    type(indices)          :: c(0:5)  
    type(bit_processing)   :: d 
    type(presorting)       :: V 

    type(sorting_list) :: srt(0:5)
       srt(1)%unrestricted =>  SRT1SS4   
       srt(2)%unrestricted =>  SRT1L1    
<<<<<<< HEAD
       srt(3)%unrestricted =>  SRT1C1    
=======
!       srt(3)%unrestricted =>  SRT1C1    
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0
       srt(4)%unrestricted =>  SRT16     
       srt(5)%unrestricted =>  SRT1R1
       srt(1)%restricted_braket   =>  SRT1TT4
!      srt(1)%restricted(2) =  SRT1LS1
<<<<<<< HEAD
       srt(1)%restricted_bra => SRT36
       srt(1)%restricted_ket => SRT22
       srt(2)%restricted_ket => SRT26
       srt(3)%restricted_ket => SRT20D
=======
       srt(2)%restricted_bra   =>  SRT36
       srt(3)%restricted_ket   =>  SRT22
       srt(4)%restricted_ket   =>  SRT26
       srt(5)%restricted_ket   =>  SRT20D
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0
!------------------------------
!      allocate(v_sorted(size(v_presorted))) 

      array_original => array_left 

       xxx => swap_character

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


    !write (*,*)'first reorder_v, probably wrong',(reorder_v(i),i=1,4)

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


!    t_transpose = (torder(1)/=array_left(reorder_v(3)))

    a%trans = num(reorder_v,.true.)    
    a%num = num(reorder_v)    


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

    call xxx(array_exchange,3,4)

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
<<<<<<< HEAD

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


!   t_transpose = (torder(1)/=array_exchange(reorder_v(3)))

=======

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


!   t_transpose = (torder(1)/=array_exchange(reorder_v(3)))

>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0
!   if (dim_t == 2) then
!     if (t_transpose) then
!      temp = reorder_v(3)
!      reorder_v(3)=reorder_v(4)
!      reorder_v(4)=temp
!     endif
!   endif

    a%num   = num(reorder_v)
    a%trans = num(reorder_v,.true.)    

<<<<<<< HEAD
     write(*,*)'reorder v', (reorder_v(i),i=1,4)
=======
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0

    call match_array (arraytype,a%num, a%trans, no_match, save_index,vtrans)


    if (no_match) call shuffle_indices (arraytype,'external', reorder_v, no_match, save_index, vtrans)

    endif

     if (no_match) then
      write (*,*) "no matching sorting subroutine available. have to do something different...."
     endif

    if (.not.no_match) then

     save_index = save_index - 1


     call integral_fetching_incore(V) 

<<<<<<< HEAD
=======
    write(*,*)'going for hunting2'
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0
     if (save_index > 0) then

     if (arraytype == 'unrestricted') then
 !call srt%a 
!     write(*,*) srt(save_index)%unrestricted
     endif

     if (arraytype == 'restricted_braket') then
 !call srt%a 
 !   call srt%reorder_index(save_index) 
!     write(*,*) srt(save_index)%restricted_braket
     endif

     if (arraytype == 'restricted_bra') then
     
<<<<<<< HEAD
  
     write(*,*)'reorder v', (reorder_v(i),i=1,4)
 
     call input_index(array_exchange,v_exchange,reorder_v,c,d)


    row = c(d%bra)%pair_combInd
    column = c(d%ket)%pair_combInd     

     write(*,*)'indices',d%I, d%J,&
        &   d%k, d%L,d%bra,d%braKet,d%ket,d%bra

     call srt(save_index)%restricted_bra(nrep, multb, ltr, lfa, c(d%I)%orbital_type, c(d%J)%orbital_type,&
        &   c(d%k)%orbital_type, c(d%L)%orbital_type,c(d%bra)%pair_combInd,c(d%braKet)%full_array,c(d%ket)%pair,&
        &   c(d%bra)%pair,V%presorted,v_sorted)  
    
    write(*,*)'successfully called the subroutine'
 
     deallocate(V%presorted)
 
=======
    write(*,*)'going for hunting..',(reorder_v(i),i=1,4),v_exchange 

   
     call input_index(array_exchange,v_exchange,reorder_v,c,d)

  write(*,*)'value of some indices',c(d%bra)%pair_combInd(1)
!     call srt(save_index)%restricted_bra(nrep, multb, ltr, lfa, a(b%I)%orbital_type, a(b%J)%orbital_type,&
!        &   c(d%k)%orbital_type, c(d%L)%orbital_type,c(d%bra)%pair_combInd,c(d%braKet)%full_array,c(d%ket)%pair,&
!        &   c(d%bra)%pair,V%presorted,v_sorted)  
    
      write(*,*)'successfully called the subroutine'
 
     deallocate(V%presorted) 
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0
     endif

     if (arraytype == 'restricted_ket') then
 !call srt%a 
 !call srt%reorder_index(save_index) 
 !    write(*,*) srt(save_index)%restricted_ket
     endif

     else
     write(*,*)'sorting is not needed'
     endif
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
    integer, dimension(size(array1)), intent(out) :: reorder_array          ! array2 is the array from where we will extract the first indices in the reordered array i.e, for left array it is always intermediate and for right always    
    integer                                  :: k,i         ! left array. array3 is the array to extract rest of the indices e.g, for left array always right and for right array always intermediate.
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

!  real*8, allocatable :: VOOO_presorted(:),VOVO_presorted(:),VVOO_presorted(:)
  type(bit_processing):: b,c 
  type(presorting), intent(out) :: G
  !logical :: integrals
!-------------------------------------------------------------
#include "symm.inc"

    bit_original = bitString (array_original, (/1,2,3,4/)) 

    b%bra = ior(bit_original(1),bit_original(2)) + iand(bit_original(1),bit_original(2))  
    b%ket = ior(bit_original(3),bit_original(4)) + iand(bit_original(3),bit_original(4)) 
    value = popcnt(bin2dec(bit_original))+ior(iand(bit_original(1),bit_original(2)),iand(bit_original(3),bit_original(4))) 

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

!       select type(sh)
!         type is(presorting)
         allocate(G%presorted(NV2))
!       end select

<<<<<<< HEAD

       call SRT1C1 (NREP,NVO,NOOT,VVOOO,G%presorted) 
=======
       call SRT1C1 (NREP,MVO,MOO,VVOOO,G%presorted) 
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0

       deallocate(VVOOO) 

!       call SRT1C1(NREP,MVO,MOO,VVOOO,VOOO_presorted) 
      endif
!      if (bin2dec((/bit_available(2),bit_available(1)/))== bin2dec((/bit_original(1),bit_original(2)/))) call
    case (2)
      allocate(VVVOO(NV3))
 !     allocate(G%presorted(NV3))
      call getVVOO(VVVOO) 
      
!      v_presorted => VVVOO 
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


  subroutine input_index(array_exchange,v_exchange,reorder_v,a,b)

!-------------------------------
 integer, dimension(:), allocatable :: bit_reordered, bit_original
 integer, dimension(4),intent(in)   :: reorder_v 
 logical, intent(in)                :: v_exchange
 character(2), intent(in)           :: array_exchange(4)
 type(indices),intent(out)          :: a(0:5)  
 type(bit_processing),intent(out)   :: b  
 integer :: i
!-------------------------------
#include "symm.inc"
     a(0)%pair => JJOO         
     a(1)%pair => JJVO
     a(2)%pair => JJVV
     a(0)%pair_combInd => MOO
     a(1)%pair_combInd => MVO
     a(2)%pair_combInd => MVV
     a(0)%full_array => JOOOO 
     a(1)%full_array => JVOOO
     a(2)%full_array => JVOVO
     a(3)%full_array => JVVOO  
!    a(4)%full_array => JVVVO
!    a(5)%full_array => 
     a(0)%orbital_type => NO
     a(1)%orbital_type => NV
  allocate (bit_reordered(size(array_original)))
  allocate (bit_original(size(array_original)))

<<<<<<< HEAD
  if (v_exchange) then
  bit_reordered = bitString (array_exchange, reorder_v)
  bit_original  = bitString (array_exchange, (/1,2,3,4/)) 
  else
  bit_reordered  = bitString (array_original, reorder_v)
  bit_original  = bitString (array_original, (/1,2,3,4/)) 
  endif

=======
  write(*,*)'exchange',v_exchange
  write(*,*)'reorder',(array_exchange(i),i=1,4)

  if (v_exchange) then
  bit_reordered = bitString (array_exchange, reorder_v)
  else
  bit_reordered  = bitString (array_original, reorder_v)
  endif
  bit_original  = bitString (array_original, (/1,2,3,4/)) 


  write(*,*)'bit reordered', bit_reordered(4)
>>>>>>> 6da04c9f1b49f5425a4f6bb72c613c5e184d5ad0

! be careful. this is not general. I have to think of different algorithm for (2,2)->(1,3) or (3,1) case.

    b%bra = ior(bit_reordered(1),bit_reordered(2)) + iand(bit_reordered(1),bit_reordered(2))  
    b%ket = ior(bit_reordered(3),bit_reordered(4)) + iand(bit_reordered(3),bit_reordered(4)) 
    b%braKet = b%bra + b%ket + ior(iand(bit_reordered(1),bit_reordered(2)),iand(bit_reordered(3),bit_reordered(4)))  

    b%I = bit_original(1) 
    b%J = bit_original(2)
    b%K = bit_original(3)
    b%L = bit_original(4)

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

  integer function num(index_array, trans)
   
  integer, dimension(:),intent(in) :: index_array
  integer :: scratch(size(index_array))
  logical, optional, intent(in) :: trans
  logical                       :: do_trans 
  integer                       :: temp, i 

  scratch = index_array
  if (present(trans)) then
   do_trans = trans
  else
   do_trans = .false.
  endif
 
  if (do_trans) then
   do i = 1, size(index_array)/2
     temp = scratch(i)
     scratch(i) = scratch(size(index_array)/2 + i) 
     scratch(size(index_array)/2 + i) = temp
   enddo  
  endif   

    do i = size(index_array), 1 , -1
      num = num + scratch(i)*10**(size(index_array)-i)
    enddo 

  end function

    integer function bin2dec(bitstring)

    integer, intent(in) :: bitstring(:)
    integer             :: value,i 

    value=0
 
    do i = 1, size(bitstring)
     bin2dec = value + bitstring(i)*2**(i-1)
    enddo

    end function

  end module   
