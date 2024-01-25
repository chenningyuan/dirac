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
    integer, public:: ih_sym = 0 ! hermiticity
    integer, public:: tr_sym = 0 ! time-reversal symmetry
    integer, public:: irep   = 0 ! irep
    real(8), public, allocatable:: elements(:)
    real(8), public, allocatable:: elements_q(:,:,:)
  contains
    private
    procedure, pass, public:: alloc => allocate_elements
    procedure, pass, public:: alloc_q => allocate_elements_q
    procedure, pass, public:: allocated_elements
    procedure, pass, public:: allocated_elements_q
    procedure, pass, public:: scale_elements
    procedure, pass, public:: get_dmat_q_mo_sector
    procedure, pass, public:: erase   => erase_dmat
    procedure, pass, public:: display => display_dmat
    procedure, pass, public:: reshape_to_q => reshape_dmat_to_q
    procedure, pass, public:: reshape_from_q => reshape_dmat_from_q
    procedure, pass, public:: assign_ih
    procedure, pass, public:: assign_tr
    procedure, pass, public:: assign_irep
    procedure, pass, public:: symmetrize
  end type dmat_t


contains

! ------------------------------------------------------------------------------
  subroutine allocate_elements(self)
    class(dmat_t), intent(inout):: self
    integer :: n

    if (self%allocated_elements()) return

    n = self%mat_dim*self%mat_dim*self%mat_nz
    allocate(self%elements(n))
    self%elements(1:n) = 0.0d0
  end subroutine

! ------------------------------------------------------------------------------
  subroutine allocate_elements_q(self)
    class(dmat_t), intent(inout):: self

    if (self%allocated_elements_q()) return

    allocate(self%elements_q(self%mat_dim,self%mat_dim,self%mat_nz))
    self%elements_q(1:self%mat_dim,1:self%mat_dim,1:self%mat_nz) = 0.0d0
  end subroutine

! ------------------------------------------------------------------------------
  subroutine scale_elements(self, f, range_start, range_end)
    class(dmat_t), intent(inout):: self
    real(8), intent(in):: f
    integer, intent(in), optional:: range_start, range_end
    integer:: s,e

    if (.not. self%allocated_elements()) return

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


! ------------------------------------------------------------------------------
  subroutine assign_ih(self, ih)
    class(dmat_t), intent(inout):: self
    integer, intent(in):: ih
    if ((ih .eq. 1) .or. (ih .eq. -1)) then
      self%ih_sym = ih
    else
      print *, 'Illegal ih. stop'
      stop
    end if
  end subroutine

! ------------------------------------------------------------------------------
  subroutine assign_tr(self, tr)
    class(dmat_t), intent(inout):: self
    integer, intent(in):: tr
    if ((tr .eq. 1) .or. (tr .eq. -1)) then
      self%tr_sym = tr
    else
      print *, 'Illegal tr. stop'
      stop
    end if
  end subroutine

! ------------------------------------------------------------------------------
  subroutine assign_irep(self, irep)
    class(dmat_t), intent(inout):: self
    integer, intent(in):: irep
    if ((irep < 0) .or. (irep > 7)) then
      print *, 'Illegal irep. stop'
      stop
    else
      self%irep = irep
    end if
  end subroutine

! ------------------------------------------------------------------------------
  subroutine symmetrize(self)
    class(dmat_t), intent(inout):: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine get_dmat_q_mo_sector(self, sector, ioffs, nmos, nfsym)
    class(dmat_t), intent(inout):: self
    character(3), intent(in):: sector
    integer, intent(in):: ioffs(6), nmos(2), nfsym
    integer:: ns(2), pi(2), ps(2), n_mo(2)
    integer:: i,j,isym,iz
!
!   return dmat with all elements zeroed-out
!   except the ones defined by sector
!

    if (.not. self%allocated_elements()) return

    if (.not. self%allocated_elements_q()) then
      call self%alloc_q()
      call self%reshape_to_q()
    end if

    ns(1)   = ioffs(1)
    pi(1)   = ioffs(2) 
    ps(1)   = ioffs(3) 
    n_mo(1) = nmos(1)
    if (nfsym == 2) then
      ns(2)   = ioffs(4)
      pi(2)   = ioffs(5) 
      ps(2)   = ioffs(6) 
      n_mo(1) = nmos(2)
    end if

!   todo: consistency check

    select case (sector)
      case ('ii ')
      ! get inactive-inactive block

        do iz = 1, self%mat_nz
          do isym = 1, nfsym

            do j = ns(isym), n_mo(isym)
              do i = ns(isym), pi(isym) - 1
                self%elements_q(i, j, iz) = 0.0d0
              end do
              do i = ps(isym), n_mo(isym)
                self%elements_q(i, j, iz) = 0.0d0
              end do
            end do

            do j = ns(isym), pi(isym) - 1
              do i = pi(isym), ps(isym) - 1
                self%elements_q(i, j, iz) = 0.0d0
              end do
            end do

            do j = ps(isym), n_mo(isym)
              do i = pi(isym), ps(isym) - 1
                self%elements_q(i, j, iz) = 0.0d0
              end do
            end do
          end do
        end do

      case ('si ')
      ! get secondary-inactive and inactive-secondary blocks
      ! (both: positive and negative energy secondary blocks)

        do iz = 1, self%mat_nz
           do isym = 1, nfsym
       
              do j = pi(isym), ps(isym) - 1
                 do i = pi(isym), ps(isym) - 1
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
       
              do j = ns(isym), pi(isym) - 1
                 do i = ns(isym), pi(isym) - 1
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
                 do i = ps(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
              do j = ps(isym), n_mo(isym)
                 do i = ns(isym), pi(isym) - 1
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
                 do i = ps(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
       
           end do
        end do

      case ('si+')
      ! get secondary-inactive and inactive-secondary blocks;
      ! only positive energy secondary

        do iz = 1, self%mat_nz
           do isym = 1, nfsym
       
              do j = pi(isym), ps(isym) - 1
                 do i = pi(isym), ps(isym) - 1
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
                 do i = ns(isym), pi(isym) - 1
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
       
              do j = ns(isym), pi(isym) - 1
                 do i = ns(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
              do j = ps(isym), n_mo(isym)
                 do i = ns(isym), pi(isym) - 1
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
                 do i = ps(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
       
           end do
        end do

      case ('si-')
      ! get secondary-inactive and inactive-secondary blocks
      ! (both: positive and negative energy secondary blocks)

        do iz = 1, self%mat_nz
           do isym = 1, nfsym
       
              do j = pi(isym), ps(isym) - 1
                 do i = pi(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
                 do i = ps(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
       
              do j = ns(isym), pi(isym) - 1
                 do i = ns(isym), pi(isym) - 1
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
                 do i = ps(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
              do j = ps(isym), n_mo(isym)
                 do i = ns(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
       
           end do
        end do

      case ('pi ')
      ! get general-inactive and inactive-general blocks
      ! (both: positive and negative energy general blocks)

        do iz = 1, self%mat_nz
           do isym = 1, nfsym
       
              do j = ns(isym), pi(isym) - 1
                 do i = ns(isym), pi(isym) - 1
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
                 do i = ps(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
              do j = ps(isym), n_mo(isym)
                 do i = ns(isym), pi(isym) - 1
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
                 do i = ps(isym), n_mo(isym)
                    self%elements_q(i, j, iz) = 0.0d0
                 end do
              end do
       
           end do
        end do

      case default
        stop
    end select

    call self%reshape_from_q()

  end subroutine



! ------------------------------------------------------------------------------
  pure function allocated_elements(self)
    class(dmat_t), intent(in):: self
    logical:: allocated_elements
    allocated_elements = allocated(self%elements)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_elements_q(self)
    class(dmat_t), intent(in):: self
    logical:: allocated_elements_q
    allocated_elements_q = allocated(self%elements_q)
  end function

! ------------------------------------------------------------------------------
  subroutine erase_dmat(self)
    class(dmat_t), intent(out):: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine reshape_dmat_to_q(self)
    class(dmat_t), intent(inout):: self
    integer :: ndim, nz

    if (.not. self%allocated_elements() .and.&
        .not. self%allocated_elements_q()) return

    ndim = self%mat_dim
    nz   = self%mat_nz

    self%elements_q  = reshape(self%elements,(/ndim,ndim,nz/))

  end subroutine

! ------------------------------------------------------------------------------
  subroutine reshape_dmat_from_q(self)
    class(dmat_t), intent(inout):: self
    integer :: ndim, nz

    if (.not. self%allocated_elements() .and.&
        .not. self%allocated_elements_q()) return

    ndim = self%mat_dim
    nz   = self%mat_nz

    self%elements  = reshape(self%elements_q,(/ndim*ndim*nz/))

  end subroutine

! ------------------------------------------------------------------------------
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
    if (self%allocated_elements()) then
      do i = 1, self%mat_dim*self%mat_dim*self%mat_nz
        write(*, '(a,i6,es16.8e2)') 'elements(i) = ', i, self%elements(i)
      end do
    end if
    if (self%allocated_elements_q()) then
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


end module
