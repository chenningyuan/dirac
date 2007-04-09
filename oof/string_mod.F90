  module mystring

 type string
  private
   character*1, dimension (:), pointer :: strp  => NULL()
   logical :: is_alloc = .false.
 end type string

 contains 

subroutine make_string(self, text)
 type(string), dimension(:), intent(inout) :: self
 character(*), intent(in) :: text

 integer :: i, alloc_stat

 if (.not.self(1)%is_alloc) then
   allocate( self(1)%strp(len(text)), stat = alloc_stat )
   if ( alloc_stat > 0 ) then
     write(*,*) 'alloc_stat erro !'
     stop 
   endif

  ! copy the text ....
  do i = 1, len(text)
    self(1)%strp(i) = text(i:i)
  enddo
   self(1)%is_alloc=.true.
 else
  write(*,*) 'string is allocated: ',self(1)%strp
 endif

 end subroutine make_string

 subroutine cancel_string(self)
 type(string), dimension(:), intent(inout) :: self
 integer :: dealloc_stat

     deallocate(self(1)%strp, stat = dealloc_stat)
     if (dealloc_stat > 0) then
       write(*,*) 'dellocation error !'; stop
     endif

    nullify( self(1)%strp )
    self(1)%is_alloc=.false.

 end subroutine cancel_string

 subroutine print_string(self)
 type(string), dimension(:), intent(in) :: self
   if (self(1)%is_alloc) then
    print *,(self(1)%strp)
   endif
 end subroutine print_string

  end module mystring

