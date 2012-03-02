module openrsp_interface_response

   use dirac_interface
   use memory_allocator
   use matrix_defop
   use openrsp_cfg

   implicit none

   public response_solver

   private

#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dgroup.h"

contains

  subroutine set_orbital_rotation_indices(irep,       &
                                          length,     &
                                          from_block, &
                                          to_block,   &
                                          orbital_rotation_indices)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: irep
    integer, intent(in)  :: length
    integer, intent(in)  :: from_block
    integer, intent(in)  :: to_block
    integer, intent(out) :: orbital_rotation_indices(2, length)
!   ----------------------------------------------------------------------------
    integer, parameter   :: ns = 1
    integer, parameter   :: pi = 2
    integer, parameter   :: ps = 3

    integer, parameter   :: g  = 1
    integer, parameter   :: u  = 2

    logical              :: gu_combine(2, 2)

    integer              :: nr_mo(2, 3)       = 0
    integer              :: range_mo(2, 3, 2) = 0

    integer              :: i, s, k
    integer              :: m, m1, m2
    integer              :: gu, gu1, gu2
    integer              :: block, block1, block2

    logical              :: debug_me = .false.
!   ----------------------------------------------------------------------------

    nr_mo(g, ns) = npsh(g)
    nr_mo(g, pi) = nish(g)
    nr_mo(g, ps) = nesh(g) - nish(g)

    nr_mo(u, ns) = npsh(u)
    nr_mo(u, pi) = nish(u)
    nr_mo(u, ps) = nesh(u) - nish(u)

    k = 0
    do gu = 1, 2
      do block = 1, 3
        do m = 1, nr_mo(gu, block)
          k = k + 1
          if (k > range_mo(gu, block, 1)) range_mo(gu, block, 1) = k
          range_mo(gu, block, 2) = k
        end do
      end do
    end do

    gu_combine = .false.
    if (jbtof(irep, 1) == 2) then
!     ungerade perturbation
      gu_combine(1, 2) = .true.
      gu_combine(2, 1) = .true.
    else
!     gerade perturbation
      gu_combine(1, 1) = .true.
      gu_combine(2, 2) = .true.
    end if

    k = 0
    i = 0
    do gu1 = 1, 2
      do block1 = 1, 3
        do m1 = 1, nr_mo(gu1, block1)
          i = i + 1
          s = 0
          do gu2 = 1, 2
            do block2 = 1, 3
              if (gu_combine(gu1, gu2)         &
                  .and. (block1 == from_block) &
                  .and. (block2 == to_block)) then
                do m2 = 1, nr_mo(gu2, block2)
                  k = k + 1
                  s = s + 1
                  orbital_rotation_indices(1, k) = i
                  orbital_rotation_indices(2, k) = s
                end do
              else
                s = s + nr_mo(gu2, block2)
              end if
            end do
          end do
        end do
      end do
    end do

    if (debug_me) then
      call header('debug orbital_rotation_indices', -1)
      do k = 1, length
        write(*, *) k, orbital_rotation_indices(1, k), &
                 '->', orbital_rotation_indices(2, k)
      end do
    end if

  end subroutine

  subroutine response_solver(freq, RHS, Dp)

!   1. RHS (AO) -> RHS (MO)
!   2. RHS (MO) -> property gradient
!   3. solve response equation, in:  property gradient
!                               out: response vector
!   4. response vector -> Wp (MO)
!   5. Wp (MO)         -> Dp (AO)

!   ----------------------------------------------------------------------------
    real(8),      intent(in)    :: freq
    type(matrix), intent(inout) :: RHS
    type(matrix), intent(inout) :: Dp
!   ----------------------------------------------------------------------------
    type(matrix)                :: RHS_mo, Wp
    type(matrix)                :: C, Cig, Ciu, Csg, Csu

!                                  for sf calcs
    integer,      allocatable   :: ibeig(:)
    real(8),      allocatable   :: mo_coef(:)
    integer                     :: i1, i2, irep

    integer,      allocatable   :: ibtyp(:)
    integer,      allocatable   :: ibtyp_pointer_pp(:)
    integer,      allocatable   :: ibtyp_pointer_pn(:)

    real(8),      allocatable   :: eigval(:)
    real(8),      allocatable   :: eigvec(:)
    real(8),      allocatable   :: convergence(:)

    real(8),      allocatable   :: response_vector_pph(:, :, :)
    real(8),      allocatable   :: response_vector_ppa(:, :, :)
    real(8),      allocatable   :: response_vector_pnh(:, :, :)
    real(8),      allocatable   :: response_vector_pna(:, :, :)

    real(8),      pointer       :: work(:)
    integer                     :: lwork, kfree, lfree

    real(8),      allocatable   :: prop_gradient_pp(:, :)
    real(8),      allocatable   :: prop_gradient_pn(:, :)

    integer,      allocatable   :: orbital_rotation_indices_pp(:, :)
    integer,      allocatable   :: orbital_rotation_indices_pn(:, :)

    integer                     :: length_pp
    integer                     :: length_pn

    integer                     :: iz, k, i, s, off
    integer                     :: nr_freq

    logical                     :: debug_me = .false.

    real(8),      allocatable   :: cmo_from_file(:)
!   ----------------------------------------------------------------------------

!   fixme: nssh not always properly set, change to nesh minus nish

    if (jbtof(RHS%irep, 1) == 2) then
!     ungerade perturbation
      length_pp = nish(1)*nssh(2) + nish(2)*nssh(1)
      length_pn = nish(1)*npsh(2) + nish(2)*npsh(1)
    else
!     gerade perturbation
      length_pp = nish(1)*nssh(1) + nish(2)*nssh(2)
      length_pn = nish(1)*npsh(1) + nish(2)*npsh(2)
    end if

    nzconf  = 0
    nzxope  = 0
    nzxopp  = 0

    if (include_pp_rotations) nzxope = length_pp
    if (include_pn_rotations) nzxopp = length_pn

    nzxopt  = nzxope + nzxopp

    nzconfq = nzconf*nz
    nzxopeq = nzxope*nz
    nzxoppq = nzxopp*nz
    nzxoptq = nzxopt*nz

    nzvar   = nzxopt
    nzvarq  = nzxoptq

    call alloc(orbital_rotation_indices_pp, 2, length_pp)
    if (length_pp > 0) then
       call set_orbital_rotation_indices(RHS%irep,  &
                                         length_pp, &
                                         2,         &
                                         3,         &
                                         orbital_rotation_indices_pp)
    end if

    call alloc(orbital_rotation_indices_pn, 2, length_pn)
    if (length_pn > 0) then
       call set_orbital_rotation_indices(RHS%irep,  &
                                         length_pn, &
                                         2,         &
                                         1,         &
                                         orbital_rotation_indices_pn)
    end if

    if (openrsp_spinfree .and. (RHS%tr_sym == 1)) then

       call alloc(mo_coef, ncmotq)
       call alloc(ibeig, ntbas(0))

       call read_mo_coef(mo_coef)
       call read_ibeig(ibeig)

       call init_mat(RHS_mo, norbt, norbt)
       RHS_mo%elms = 0.0d0
       irep = RHS%irep

       lwork = len_f77_work
       call di_select_wrk(work, lwork)

       do i1 = 1, nfsym
          i2 = mod(i1 + jbtof(irep, 1), 2) + 1
          if (nfbas(i1, 0) > 0) then
             if (nfbas(i2, 0) > 0) then
                call qbtrans(irep,                                           &
                             'AOMO',                                         &
                             'S',                                            &
                             0.0d0,                                          &
!                            -------------------------------------------------
                             nfbas(i1, 0), nfbas(i2, 0), norb(i1), norb(i2), &
!                            -------------------------------------------------
                             RHS%elms(i2basx(i1, i2) + 1),                   &
                             ntbas(0), ntbas(0), nz, ipqtoq(1, irep),        &
!                            -------------------------------------------------
                             RHS_mo%elms(i2orbx(i1, i2) + 1),                &
                             norbt, norbt, nz, ipqtoq(1, irep),              &
!                            -------------------------------------------------
                             mo_coef(icmoq(i1) + 1),                         &
                             nfbas(i1, 0), norb(i1), nz, ipqtoq(1, 0),       &
                             ibeig(iorb(i1) + 1),                            &
!                            -------------------------------------------------
                             mo_coef(icmoq(i2) + 1),                         &
                             nfbas(i2, 0), norb(i2), nz, ipqtoq(1, 0),       &
                             ibeig(iorb(i2) + 1),                            &
!                            -------------------------------------------------
                             work, lwork, 0)
             end if
          end if
       end do

       call di_deselect_wrk(work, lwork)
       call dealloc(mo_coef)
       call dealloc(ibeig)

    else

       allocate(cmo_from_file(n2bbasxq))
       call read_mo_coef(cmo_from_file)
       call get_C(C, cmo_from_file, i=1.0d0, s=1.0d0, g=1.0d0, u=1.0d0)
       deallocate(cmo_from_file)
       RHS_mo = dag(C)*(RHS*C)
       C = 0

    end if

    RHS_mo%ih_sym = RHS%ih_sym
    RHS_mo%tr_sym = RHS%tr_sym
    RHS_mo%irep   = RHS%irep

    if (debug_me) then
      call print_mat(RHS,    label = 'debug RHS in AO basis')
      call print_mat(RHS_mo, label = 'debug RHS in MO basis')
    end if


!   get positive -> positive property gradient
!   ==========================================

    if (include_pp_rotations) then

       call alloc(prop_gradient_pp, length_pp, nz)

       do iz = 1, nz
          off = (iz - 1)*norbt*norbt
          do k = 1, length_pp
             i = orbital_rotation_indices_pp(1, k)
             s = orbital_rotation_indices_pp(2, k)
             prop_gradient_pp(k, iz) = -2.0d0*RHS_mo%elms(off + (i - 1)*norbt + s)
          end do
       end do

       if (debug_me) then
          call header('debug prop_gradient_pp', -1)
          call print_q_vector(prop_gradient_pp, length_pp)
       end if

    end if


!   get positive -> negative property gradient
!   ==========================================

    if (include_pn_rotations) then

       call alloc(prop_gradient_pn, length_pn, nz)

       do iz = 1, nz
          off = (iz - 1)*norbt*norbt
          do k = 1, length_pn
             i = orbital_rotation_indices_pn(1, k)
             s = orbital_rotation_indices_pn(2, k)
             prop_gradient_pn(k, iz) = -2.0d0*RHS_mo%elms(off + (i - 1)*norbt + s)
          end do
       end do

       if (debug_me) then
          call header('debug prop_gradient_pn', -1)
          call print_q_vector(prop_gradient_pn, length_pn)
       end if

    end if


!   inherit shape
    Wp     = RHS_mo
    RHS_mo = 0

!   fixme:
!   call with proper limit
!   at the moment i don't know the limit
!   simply give entire work and hope for the best
    lwork = len_f77_work
    call di_select_wrk(work, lwork)

    kfree = 1
    lfree = len_f77_work
    call memget('INTE', kzxope, 2*length_pp, work, kfree, lfree)
    call memget('INTE', kzxopp, 2*length_pn, work, kfree, lfree)

    if (include_pn_rotations) then
       call fill_work(length_pn, orbital_rotation_indices_pn, work(kzxope))
    end if
    if (include_pp_rotations) then
       call fill_work(length_pp, orbital_rotation_indices_pp, work(kzxopp))
    end if


!   parameters
!   ==========

    call setrsp

!   .URKBAL response with frequency=0.0 but static=.false. converges much worse
!   than frequency=0.0 with static=.true.
!   it is important to correctly set static
    static = .true.
    if (dabs(freq) > tiny(0.0d0)) then
       static = .false.
    end if

    lineq     = .true.
    lsvcfg(1) = .true.
    lsvcfg(2) = .true.
    tknorm    = .true.
    diaghe    = .true.
    iprxrs    = 0
    thcxrs    = openrsp_cfg_threshold_response
    resfac    = 1.0d3
    maxitr    = 150
    nr_freq   = 1
    maxsim    = -1
    nredm     = 400
    n2redm    = nredm*nredm
    loffty    = 0
    uncoup    = .false.
    imfreq    = .false.
    triplet   = .false.
    cnvint(1) = 1.0d20
    cnvint(2) = 1.0d20
    itrint(1) = 1
    itrint(2) = 1
    intdef    = integral_flag
    sternh    = .false.
    sternc    = .false.
    e2chek    = .false.

!   nr_freq   = size(freq)
    jsymop    = RHS%irep + 1
    jtimop    = RHS%tr_sym
    jopsy     = jbtof(RHS%irep, 1)
    nfreq     = nr_freq
    nexsim    = nr_freq
    nexstv    = nr_freq
    nexcnv    = nr_freq
    ncred     = 0
    nered     = 0
    npred     = 0
    nzred     = ncred + nered + npred


!   solve response equation in the reduced space spanned by the trial vectors
!   =========================================================================

    call alloc(ibtyp,            2*nredm)
    call alloc(ibtyp_pointer_pp, nredm)
    call alloc(ibtyp_pointer_pn, nredm)
    call alloc(convergence,      nr_freq)
    call alloc(eigval,           nr_freq)
    call alloc(eigvec,           nredm*nr_freq)

    eigval = freq

    call xrsctl((/0.0d0/),        &
                prop_gradient_pp, &
                prop_gradient_pn, &
                ibtyp,            &
                (/0/),            &
                ibtyp_pointer_pp, &
                ibtyp_pointer_pn, &
                convergence,      &
                eigval,           &
                eigvec,           &
                work, kfree, lfree)

    call dealloc(convergence)
    call dealloc(eigval)
    call di_deselect_wrk(work, lwork)


!   construct response vector
!   =========================

    if (include_pp_rotations) then

      call alloc(response_vector_pph, length_pp, nz, nr_freq)
      call alloc(response_vector_ppa, length_pp, nz, nr_freq)

      response_vector_pph = 0.0d0
      response_vector_ppa = 0.0d0

      call construct_response_vector('pp',                &
                                     1,                   &
                                     response_vector_pph, &
                                     length_pp,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pp,    &
                                     eigvec)

      call construct_response_vector('pp',                &
                                     -1,                  &
                                     response_vector_ppa, &
                                     length_pp,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pp,    &
                                     eigvec)

      if (debug_me) then
        call header('debug response_vector_pph', -1)
        call prbvec(6, response_vector_pph, nr_freq, length_pp)
        call header('debug response_vector_ppa', -1)
        call prbvec(6, response_vector_ppa, nr_freq, length_pp)
      end if

    end if

    if (include_pn_rotations) then

      call alloc(response_vector_pnh, length_pn, nz, nr_freq)
      call alloc(response_vector_pna, length_pn, nz, nr_freq)

      response_vector_pnh = 0.0d0
      response_vector_pna = 0.0d0

      call construct_response_vector('pn',                &
                                     1,                   &
                                     response_vector_pnh, &
                                     length_pn,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pn,    &
                                     eigvec)

      call construct_response_vector('pn',                &
                                     -1,                  &
                                     response_vector_pna, &
                                     length_pn,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pn,    &
                                     eigvec)

      if (debug_me) then
        call header('debug response_vector_pnh', -1)
        call prbvec(6, response_vector_pnh, nr_freq, length_pn)
        call header('debug response_vector_pna', -1)
        call prbvec(6, response_vector_pna, nr_freq, length_pn)
      end if

    end if

    call dealloc(ibtyp)
    call dealloc(ibtyp_pointer_pp)
    call dealloc(ibtyp_pointer_pn)
    call dealloc(eigvec)


!   scatter response vectors
!   ========================

    Wp%elms = 0.0d0

    if (include_pp_rotations) then
      call scatter_vector(length_pp,                   &
                          orbital_rotation_indices_pp, &
                          1.0d0,                       &
                          response_vector_pph,         &
                          Wp%elms,                     &
                          Wp%irep)
      call scatter_vector(length_pp,                   &
                          orbital_rotation_indices_pp, &
                         -1.0d0,                       &
                          response_vector_ppa,         &
                          Wp%elms,                     &
                          Wp%irep)

      call dealloc(orbital_rotation_indices_pp)
      call dealloc(response_vector_pph)
      call dealloc(response_vector_ppa)
    end if

    if (include_pn_rotations) then
      call scatter_vector(length_pn,                   &
                          orbital_rotation_indices_pn, &
                          1.0d0,                       &
                          response_vector_pnh,         &
                          Wp%elms,                     &
                          Wp%irep)
      call scatter_vector(length_pn,                   &
                          orbital_rotation_indices_pn, &
                         -1.0d0,                       &
                          response_vector_pna,         &
                          Wp%elms,                     &
                          Wp%irep)

      call dealloc(orbital_rotation_indices_pn)
      call dealloc(response_vector_pnh)
      call dealloc(response_vector_pna)
    end if


!   get coefficients
!   ================

    allocate(cmo_from_file(n2bbasxq))
    call read_mo_coef(cmo_from_file)
    call get_C(Cig, cmo_from_file, i=1.0d0, s=0.0d0, g=1.0d0, u=0.0d0)
    call get_C(Csg, cmo_from_file, i=0.0d0, s=1.0d0, g=1.0d0, u=0.0d0)
    if (nfsym == 2) then
      call get_C(Ciu, cmo_from_file, i=1.0d0, s=0.0d0, g=0.0d0, u=1.0d0)
      call get_C(Csu, cmo_from_file, i=0.0d0, s=1.0d0, g=0.0d0, u=1.0d0)
    end if
    deallocate(cmo_from_file)


!   construct perturbed AO density matrix
!   =====================================

    if (nfsym == 2) then
      if (jbtof(Wp%irep, 1) == 2) then
!       ungerade perturbation
        Dp = (Cig*(Wp*dag(Csu))) &
           + (Ciu*(Wp*dag(Csg))) &
           - (Csg*(Wp*dag(Ciu))) &
           - (Csu*(Wp*dag(Cig)))
      else
!       gerade perturbation
        Dp = (Cig*(Wp*dag(Csg))) &
           + (Ciu*(Wp*dag(Csu))) &
           - (Csg*(Wp*dag(Cig))) &
           - (Csu*(Wp*dag(Ciu)))
      end if
      Ciu = 0
      Csu = 0
    else
      Dp = (Cig*(Wp*dag(Csg))) &
         - (Csg*(Wp*dag(Cig)))
    end if
    Cig = 0
    Csg = 0

    Dp%ih_sym = Wp%ih_sym
    Dp%tr_sym = Wp%tr_sym
    Dp%irep   = Wp%irep

    if (debug_me) then
       call print_mat(Wp, label = 'debug Wp')
    end if

    Wp = 0

    if (debug_me) then
       call print_mat(Dp, label = 'debug Dp')
    end if

  end subroutine

  subroutine construct_response_vector(vector_type,     &
                                       ih,              &
                                       response_vector, &
                                       length,          &
                                       nr_freq,         &
                                       ibtyp,           &
                                       ibtyp_pointer,   &
                                       eigvec)

!   ----------------------------------------------------------------------------
    character(*), intent(in)  :: vector_type
    integer,      intent(in)  :: ih
    real(8),      intent(out) :: response_vector(*)
    integer,      intent(in)  :: length
    integer,      intent(in)  :: nr_freq
    integer,      intent(in)  :: ibtyp(*)
    integer,      intent(in)  :: ibtyp_pointer(*)
    real(8),      intent(in)  :: eigvec(*)
!   ----------------------------------------------------------------------------
    integer,      allocatable :: ivecs(:)
    real(8),      allocatable :: buffer(:)
    character(6)              :: file_name
    integer                   :: file_unit, nbtyp
!   ----------------------------------------------------------------------------

    select case (vector_type)
      case ('pp')
        nbtyp     = 1
        file_unit = luboe
        file_name = 'PAMBOE'
      case ('pn')
        nbtyp     = 2
        file_unit = lubop
        file_name = 'PAMBOP'
      case default
        call quit('error in construct_response_vector: unknown vector type')
    end select

    call alloc(ivecs,  nr_freq)
    call alloc(buffer, length*nz)

    open(file_unit,              &
         file   = file_name,     &
         form   = 'unformatted', &
         access = 'direct',      &
         recl   = 8*length*nz,   &
         status = 'old')

    call xrsxv1(ih,              &
                nbtyp,           &
                response_vector, &
                eigvec,          &
                nr_freq,         &
                ibtyp,           &
                ibtyp_pointer,   &
                ivecs,           &
                buffer)

    close(file_unit, status = 'keep')

    call dealloc(ivecs)
    call dealloc(buffer)

  end subroutine

  subroutine print_vector(vector, vector_dimension)

    real(8), intent(in) :: vector(:)
    integer,       intent(in) :: vector_dimension

    call prqmat(vector,           &
                vector_dimension, &
                1,                &
                vector_dimension, &
                1,                &
                1,                &
                ipqtoq(1, 0),     &
                6)

  end subroutine

  subroutine print_q_vector(vector, vector_dimension)

    real(8), intent(in) :: vector(*)
    integer,       intent(in) :: vector_dimension

    call prqmat(vector,           &
                vector_dimension, &
                1,                &
                vector_dimension, &
                1,                &
                nz,               &
                ipqtoq(1, 0),     &
                6)

  end subroutine

  subroutine scatter_vector(length,                   &
                            orbital_rotation_indices, &
                            h,                        &
                            response_vector,          &
                            matrix,                   &
                            irep)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: length
    integer, intent(in)  :: orbital_rotation_indices(2, length)
    real(8), intent(in)  :: h
    real(8), intent(in)  :: response_vector(length, nz)
    real(8), intent(out) :: matrix(norbt, norbt, nz)
    integer, intent(in)  :: irep
!   ----------------------------------------------------------------------------
    integer              :: i, s, is, iz
    real(8)              :: f
!   ----------------------------------------------------------------------------

    do is = 1, length

      i = orbital_rotation_indices(1, is)
      s = orbital_rotation_indices(2, is)

      do iz = 1, nz

        if (ipqtoq(iz, irep) > 1) then
          f = -1.0
        else
          f =  1.0
        end if

!              row
!              ¦
!              ¦  column
!              ¦  ¦
        matrix(s, i, iz) = matrix(s, i, iz) &
                         +     response_vector(is, iz)
        matrix(i, s, iz) = matrix(i, s, iz) &
                         - f*h*response_vector(is, iz)
      end do
    end do

  end subroutine



end module

subroutine fill_work(length, orbital_rotation_vector, jxop)

   implicit none
 
   integer :: length
   integer :: orbital_rotation_vector(2, *)
   integer :: jxop(2, *)
   integer :: i
 
   do i = 1, length
      jxop(1, i) = orbital_rotation_vector(1, i)
      jxop(2, i) = orbital_rotation_vector(2, i)
   end do

end subroutine
