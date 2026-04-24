module exacorr_checkpoint_mo

  use checkpoint_object, only: checkpoint_object_t
  use exacorr_mo, only: cmo, alloc_mo, copy_mo, dealloc_mo, &
      convert_mo_to_complex, convert_mo_to_quaternion

  implicit none

  private

  public :: checkpoint_mo_t

  type, extends(checkpoint_object_t) :: checkpoint_mo_t
     type(cmo) :: mo
     logical   :: coeff_loaded = .false.
  contains
     procedure :: init => checkpoint_mo_init
     procedure :: read_from_checkpoint => checkpoint_mo_read
     procedure :: write_to_checkpoint => checkpoint_mo_write
     procedure :: mark_all_dirty => checkpoint_mo_mark_all_dirty
  end type checkpoint_mo_t

contains

  subroutine checkpoint_mo_init(self, root)

    class(checkpoint_mo_t), intent(inout) :: self
    character(len=*), intent(in)          :: root

    call self%set_root(root)
    call register_mobasis_fields(self)

  end subroutine checkpoint_mo_init

  subroutine checkpoint_mo_read(self, read_coeff)

    class(checkpoint_mo_t), intent(inout) :: self
    logical, intent(in), optional         :: read_coeff

    integer                              :: n_basis(2), n_mo(2), n_po(2)
    integer                              :: nz, n_fsym
    integer                              :: nao, nmo, npo, nto
    integer                              :: ncoeff, i
    real(8), allocatable                 :: coeff(:), eig(:), coeff_q(:,:,:)
    integer, allocatable                 :: symmetry(:), index(:)
    logical                              :: c1_info, need_coeff

    call self%read_integer_array('n_basis', n_basis)
    call self%read_integer_array('n_mo', n_mo)
    call self%read_integer_array('n_po', n_po)
    call self%read_integer('nz', nz)
    call self%read_integer('n_fsym', n_fsym)

    nao = sum(n_basis(1:n_fsym))
    nmo = sum(n_mo(1:n_fsym))
    npo = sum(n_po(1:n_fsym))
    nto = npo + nmo

    allocate(eig(nto))
    allocate(symmetry(nto))
    allocate(index(nmo))

    need_coeff = .true.
    if (present(read_coeff)) need_coeff = read_coeff

    call self%query_field('eigenvalues_C1', c1_info)
    if (c1_info) then
       call self%read_real_array('eigenvalues_C1', eig)
    else
       call self%read_real_array('eigenvalues', eig)
    end if
    call self%read_integer_array('symmetry', symmetry)

    do i = 1, nmo
       index(i) = i
    end do

    if (self%mo%allocated) call dealloc_mo(self%mo)
    self%mo%algebra = 'c'
    call alloc_mo(self%mo, nao, nmo)
    self%mo%nz = nz
    self%mo%index = index
    self%mo%energy = eig(npo+1:nto)
    self%mo%symmetry = symmetry(npo+1:nto)
    self%coeff_loaded = .false.

    if (need_coeff) then
       ncoeff = nao * nto * nz
       allocate(coeff(ncoeff))
       allocate(coeff_q(nao,nto,nz))

       call self%query_field('orbitals_C1', c1_info)
       if (c1_info) then
          call self%read_real_array('orbitals_C1', coeff)
       else
          call self%read_real_array('orbitals', coeff)
       end if

       deallocate(self%mo%coeff)
       self%mo%algebra = 'q'
       allocate(self%mo%coeff_q(self%mo%nao, self%mo%nmo, nz))
       coeff_q = reshape(coeff, (/nao, nto, nz/))
       self%mo%coeff_q = coeff_q(:,npo+1:nto,:)
       call convert_mo_to_complex(self%mo)
       self%coeff_loaded = .true.

       deallocate(coeff, coeff_q)
    end if

    deallocate(eig, symmetry, index)

  end subroutine checkpoint_mo_read

  subroutine checkpoint_mo_write(self, mobasis_id)

    class(checkpoint_mo_t), intent(inout) :: self
    integer, intent(in)                   :: mobasis_id

    type(cmo)                             :: qorbital
    integer                               :: ncoeff

    call copy_mo(self%mo, qorbital)
    call convert_mo_to_quaternion(qorbital)

    ncoeff = qorbital%nao * qorbital%nmo * qorbital%nz

    call self%write_integer('mobasis_id', mobasis_id)
    call self%write_real_array('eigenvalues', qorbital%energy)
    call self%write_integer('nz', qorbital%nz)
    call self%write_integer('n_fsym', 1)
    call self%write_integer('n_basis', qorbital%nao)
    call self%write_integer('n_mo', qorbital%nmo)
    call self%write_integer('n_po', 0)
    call self%write_integer_array('symmetry', qorbital%symmetry)
    call self%write_real_array('orbitals', reshape(qorbital%coeff_q, (/ncoeff/)))

    call dealloc_mo(qorbital)

  end subroutine checkpoint_mo_write

  subroutine checkpoint_mo_mark_all_dirty(self)

    class(checkpoint_mo_t), intent(inout) :: self
    integer                               :: i

    if (.not.allocated(self%fields)) return
    do i = 1, size(self%fields)
       call self%mark_dirty(self%fields(i)%name)
    end do

  end subroutine checkpoint_mo_mark_all_dirty

  subroutine register_mobasis_fields(self)

    class(checkpoint_mo_t), intent(inout) :: self

    call self%register_field('mobasis_id', 'integer', 'single', required=.true.)
    call self%register_field('n_basis', 'integer', 'single', required=.true.)
    call self%register_field('n_mo', 'integer', 'single', required=.true.)
    call self%register_field('n_po', 'integer', 'single', required=.false.)
    call self%register_field('nz', 'integer', 'single', required=.true.)
    call self%register_field('n_fsym', 'integer', 'single', required=.true.)
    call self%register_field('orbitals', 'real', 'array', required=.true.)
    call self%register_field('eigenvalues', 'real', 'array', required=.true.)
    call self%register_field('symmetry', 'integer', 'array', required=.false.)
    call self%register_field('orbitals_C1', 'real', 'array', required=.false.)
    call self%register_field('eigenvalues_C1', 'real', 'array', required=.false.)

  end subroutine register_mobasis_fields

end module exacorr_checkpoint_mo
