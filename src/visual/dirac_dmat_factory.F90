module dirac_dmat_factory
!
!   gosia: this is a temporary solution for handling density matrices;
!   it is now used in VISUAL and FDE;
!   also see https://gitlab.com/dirac/dirac-private/-/issues/550
!
  implicit none
  private

  public dmat_t

  type dmat_t
    private
    integer, public:: mat_dim = 0
    integer, public:: mat_nz = 0
    real(8), public, pointer:: elements(:)
    real(8), public, pointer:: elements_q(:,:,:)
  contains
    private
    procedure, pass, public:: alloc => allocate_elements
    procedure, pass, public:: alloc_q => allocate_elements_q
    procedure, pass, public:: associated_elements
    procedure, pass, public:: associated_elements_q
    procedure, pass, public:: scale_elements
    procedure, pass, public:: get_dmat_sector
    procedure, pass, public:: erase   => erase_dmat
    procedure, pass, public:: display => display_dmat
    procedure, pass, public:: reshape_to_q => reshape_dmat_to_q
  end type dmat_t

!  type dmat_t
!    private
!    integer :: mat_dim = 0
!    integer :: mat_nz = 0
!    real(8), pointer :: elements(:) ! or allocatable??
!    real(8), pointer :: elements_q(:,:,:)
!  contains
!    private
!    procedure, pass, public:: is_allocated
!    procedure, pass, public:: alloc => allocate_elements
!    procedure, pass, public:: mo_to_ao => transform_from_mo_to_ao_c1
!    !procedure, pass, public:: insert_symmetry_info
!    !procedure, pass, public:: set_dimensions
!    !procedure, pass, public:: alloc_q => allocate_elements_q
!    !procedure, pass, public:: associated_elements_q
!    procedure, pass, public:: scale_elements
!    !procedure, pass, public:: get_dmat_sector
!    procedure, pass, public:: erase   => erase_dmat
!    procedure, pass, public:: display => display_dmat
!    !procedure, pass, public:: reshape_to_q => reshape_dmat_to_q
!  end type dmat_t

contains

  subroutine allocate_elements(self)
    class(dmat_t), intent(inout):: self
    integer :: n

    if (self%associated_elements()) return

    n = self%mat_dim*self%mat_dim*self%mat_nz
    nullify(self%elements)
    allocate(self%elements(n))
    self%elements(1:n) = 0.0d0
  end subroutine

  subroutine allocate_elements_q(self)
    class(dmat_t), intent(inout):: self

    if (self%associated_elements_q()) return

    nullify(self%elements_q)
    allocate(self%elements_q(self%mat_dim,self%mat_dim,self%mat_nz))
    self%elements_q(1:self%mat_dim,1:self%mat_dim,1:self%mat_nz) = 0.0d0
  end subroutine

  subroutine scale_elements(self, f, range_start, range_end)
    class(dmat_t), intent(inout):: self
    real(8), intent(in):: f
    integer, intent(in), optional:: range_start, range_end
    integer:: s,e

    if (.not. self%associated_elements()) return

    if (present(range_start)) then
      s = range_start
    else
      s = 1
    end if

    if (present(range_end)) then
      e = range_end
    else
      e = size(self%elements)
    end if

    self%elements(s:e) = f * self%elements(s:e)
  end subroutine


  subroutine get_dmat_sector(self, sector)
    class(dmat_t), intent(inout):: self
    character(2), intent(in):: sector
    integer:: s1,s2,e1,e2

    if (.not. self%associated_elements()) return

    if (.not. self%associated_elements_q()) then
      call self%alloc_q()
      call self%reshape_to_q()
    end if

    !select case (sector)
    !  case ('ii')
    !    s1 = 1
    !    e1 = s1 + 
    !    s2 = 
    !    e2 = 
    !    self%elements_q(s1:e1,s2:e2,1:self%mat_nz) = 0.0d0
    !  case ('ip')
    !  case default
    !    if (present(range_start)) then
    !    self%elements_q(s1:e1,s2:e2,1:self%mat_nz) = 0.0d0
    !end select
  end subroutine



  pure function associated_elements(self)
    class(dmat_t), intent(in):: self
    logical:: associated_elements
    associated_elements = associated(self%elements)
  end function

  pure function associated_elements_q(self)
    class(dmat_t), intent(in):: self
    logical:: associated_elements_q
    associated_elements_q = associated(self%elements_q)
  end function

  subroutine erase_dmat(self)
    class(dmat_t), intent(out):: self
  end subroutine

  subroutine reshape_dmat_to_q(self)
    class(dmat_t), intent(inout):: self
    integer :: ndim, nz

    if (.not. self%associated_elements() .and.&
        .not. self%associated_elements_q()) return

    ndim = self%mat_dim
    nz   = self%mat_nz

    self%elements_q  = reshape(self%elements,(/ndim,ndim,nz/))

  end subroutine

  subroutine display_dmat(self,node_id)
    class(dmat_t), intent(in):: self
    integer, optional, intent(in):: node_id
    integer:: i, j, k
    write(*, *)
    if (present(node_id)) then
      write(*, *) ' DMAT on node : ', node_id
    else
      write(*, *) ' DMAT: '
    end if
    write(*, *) ' mat_dim = ', self%mat_dim
    write(*, *) ' mat_nz  = ', self%mat_nz
    if (self%associated_elements()) then
      do i = 1, self%mat_dim*self%mat_dim*self%mat_nz
        write(*, '(a,i6,es16.8e2)') 'elements(i) = ', i, self%elements(i)
      end do
    end if
    if (self%associated_elements_q()) then
      do k = 1, self%mat_nz
        do j = 1, self%mat_dim
          do i = 1, self%mat_dim
            write(*, '(a,3i6,es16.8e2)') 'elements_q(i,j,k) = ', i, j, k, self%elements_q(i,j,k)
          end do
        end do
      end do
    end if
    write(*, *)
  end subroutine



! NEW VERSION:
!! ------------------------------------------------------------------------------
!  subroutine allocate_elements(self)
!!OK
!    class(dmat_t), intent(inout):: self
!    integer :: ndim
!
!    if (self%is_allocated()) return
!
!    ndim = self%mat_dim*self%mat_dim*self%mat_nz
!    nullify(self%elements)
!    allocate(self%elements(ndim))
!    self%elements(1:ndim) = 0.0d0
!  end subroutine
!
!! ------------------------------------------------------------------------------
!  pure function is_allocated(self)
!!OK
!    class(dmat_t), intent(in):: self
!    logical:: is_allocated
!    is_allocated = associated(self%elements)
!  end function
!
!
!! ------------------------------------------------------------------------------
!  subroutine transform_from_mo_to_ao_c1(cmo, mat_mo, mat_ao)
!!OK
!
!    real(8), intent(in):: cmo
!    type(dmat_t), intent(inout):: dmat_ao
!    type(dmat_t), intent(inout):: dmat_mo
!    integer:: iprint
!    integer:: nao, nmo, nz
!!
!!   dmat_ao(mu,nu)=cmo(mu,p)*dmat_mo(p,q)*cmo^{-1}(q,nu)
!!   gosia fixme: no symmetry
!! ------------------------------------------------------------------------------
!
!    if (.not. mat_ao%is_allocated()) stop
!    if (.not. mat_mo%is_allocated()) stop
!    if (mat_mo%mat_nz /= mat_ao%mat_nz) stop
!
!    iprint = 0
!
!    nao  = mat_ao%mat_dim
!    nmo  = mat_mo%mat_dim
!    nz   = mat_ao%mat_nz
!
!    call qtrans90('MOAO','S',0.D0,nao,nao,nmo,nmo,&
!                  mat_ao%elements,nao,nao,nz,(/1,2,3,4/),& 
!                  mat_mo%elements,nmo,nmo,nz,(/1,2,3,4/),& 
!                  cmo,nao,nmo,nz,(/1,2,3,4/), &
!                  cmo,nao,nmo,nz,(/1,2,3,4/), &
!                  iprint)
!
!  end subroutine
!
!
!! ------------------------------------------------------------------------------
!  pure function associated_elements_q(self)
!    class(dmat_t), intent(in):: self
!    logical:: associated_elements_q
!    associated_elements_q = associated(self%elements_q)
!  end function
!
!! ------------------------------------------------------------------------------
!  subroutine erase_dmat(self)
!!OK
!    class(dmat_t), intent(out):: self
!  end subroutine
!
!
!  subroutine reshape_dmat_to_q(self)
!    class(dmat_t), intent(inout):: self
!    integer :: ndim, nz
!
!    if (.not. self%associated_elements() .and.&
!        .not. self%associated_elements_q()) return
!
!    ndim = self%mat_dim
!    nz   = self%mat_nz
!
!    self%elements_q  = reshape(self%elements,(/ndim,ndim,nz/))
!
!  end subroutine
!
!
!  subroutine allocate_elements_q(self)
!    class(dmat_t), intent(inout):: self
!
!    if (self%associated_elements_q()) return
!
!    nullify(self%elements_q)
!    allocate(self%elements_q(self%mat_dim,self%mat_dim,self%mat_nz))
!    self%elements_q(1:self%mat_dim,1:self%mat_dim,1:self%mat_nz) = 0.0d0
!  end subroutine
!
!! ------------------------------------------------------------------------------
!  subroutine scale_elements(self, f, range_start, range_end)
!!OK
!    class(dmat_t), intent(inout):: self
!    real(8), intent(in):: f
!    integer, intent(in), optional:: range_start, range_end
!    integer:: s,e
!
!    if (.not. self%is_allocated()) return
!
!    if (present(range_start)) then
!      s = range_start
!    else
!      s = 1
!    end if
!
!    if (present(range_end)) then
!      e = range_end
!    else
!      e = size(self%elements)
!    end if
!
!    self%elements(s:e) = f * self%elements(s:e)
!  end subroutine
!
!
!  subroutine get_dmat_sector(self, sector)
!    class(dmat_t), intent(inout):: self
!    character(2), intent(in):: sector
!    integer:: s1,s2,e1,e2
!
!    if (.not. self%associated_elements_q()) return
!
!    select case (sector)
!      case ('ii')
!        self%elements_q(s1:e1,s2:e2,1:self%mat_nz) = 0.0d0
!    end select
!  end subroutine
!
!  subroutine display_dmat(self,node_id)
!    class(dmat_t), intent(in):: self
!    integer, optional, intent(in):: node_id
!    integer:: i, j, k
!    write(*, *)
!    if (present(node_id)) then
!      write(*, *) ' DMAT on node : ', node_id
!    else
!      write(*, *) ' DMAT: '
!    end if
!    write(*, *) ' mat_dim = ', self%mat_dim
!    write(*, *) ' mat_nz  = ', self%mat_nz
!    if (self%is_allocated()) then
!      do i = 1, self%mat_dim*self%mat_dim*self%mat_nz
!        write(*, '(a,i6,es16.8e2)') 'elements(i) = ', i, self%elements(i)
!      end do
!    end if
!    !if (self%associated_elements_q()) then
!    !  do k = 1, self%mat_nz
!    !    do j = 1, self%mat_dim
!    !      do i = 1, self%mat_dim
!    !        write(*, '(a,3i6,es16.8e2)') 'elements_q(i,j,k) = ', i, j, k, self%elements_q(i,j,k)
!    !      end do
!    !    end do
!    !  end do
!    !end if
!    write(*, *)
!  end subroutine


end module
