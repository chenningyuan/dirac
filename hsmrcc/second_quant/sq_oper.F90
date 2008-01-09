Module SecQuantOperator

type, public :: SQOper
 !private
   !integer(kind=4) :: p = 0
   integer :: p = 0
   logical :: is_creation_oper = .false.
   !type(SQOper), pointer :: next => null()
end type SQOper
 
 contains

function make_sqo(i,is_creat) result (this)
 !integer(kind=4), intent(in) :: i
 integer, intent(in) :: i
 logical, optional, intent(in) :: is_creat
 type(SQOper) :: this

! write(*,'(2x,a,i3,a,l)')   &
!  'make_sqo: entering i=',i,' entering is_creat=',is_creat
  
 if (i>0) then
  this%p = i
  if ( present(is_creat) ) then
    this%is_creation_oper = is_creat
  else
    ! ... by default ...
    this%is_creation_oper = .FALSE.
  endif
 else
   stop 'make_sqo: bad entering parameter i !'
 endif
end function make_sqo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_sqo(this)
!----------------------------------------------------------------------------
!
! print sec.quantization operator info
!
!----------------------------------------------------------------------------
type(SQOper), intent(in) :: this
write(*,'(2x,a,i3,$)') ' sec.quant.oper. - index no',this%p
if (this%is_creation_oper) then
  write(*,*) ',creation oper.'
else 
  write(*,*) ',anihilation oper.'
endif
end subroutine print_sqo

End Module SecQuantOperator

