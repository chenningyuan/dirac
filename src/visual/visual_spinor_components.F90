!
! Experimental support code for writing individual four-component spinor
! components from the DIRAC VISUAL pipeline.
! Added by Aleksandr Zaichenko,2026

module visual_spinor_components
  use visual_datatypes, only: cmo_t, grid_function_t, grid_t, mol_t
  use checkpoint
  implicit none
  private

  public :: visual_spinor_components_handle_grid_function
  public :: visual_spinor_components_evaluate_grid_function
  public :: visual_spinor_components_report_selection
  public :: visual_spinor_components_locate_mo
  public :: visual_spinor_components_extract_coefficients
  public :: visual_spinor_components_split_quaternion

  type :: spinor_c1_basis_t
   integer :: n_ao = 0
   integer :: n_shells = 0
   integer :: angular = 0

   integer, allocatable :: orbmom(:)
   integer, allocatable :: n_prim(:)
   integer, allocatable :: n_cont(:)

   real(8), allocatable :: centers(:,:)       ! (3,n_shells)
   real(8), allocatable :: exponents(:)
   real(8), allocatable :: contractions(:)

   integer, allocatable :: shell_exp_offset(:)
   integer, allocatable :: shell_coeff_offset(:)

   integer, allocatable :: ao_shell(:)
   integer, allocatable :: ao_cont(:)
   integer, allocatable :: ao_lx(:)
   integer, allocatable :: ao_ly(:)
   integer, allocatable :: ao_lz(:)
  end type spinor_c1_basis_t


contains

  ! ---------------------------------------------------------------------------
  subroutine visual_spinor_components_handle_grid_function(cmo, grid_function, print_unit)
    !
    ! Entry point called from VISUAL for every grid function.
    !
    ! At this stage this is only a diagnostic report.  In the next stage this
    ! routine will dispatch the actual spinor-component evaluation/writer.
    !
    type(cmo_t), intent(in) :: cmo
    type(grid_function_t), intent(in) :: grid_function
    integer, intent(in), optional :: print_unit

    integer :: lu

    lu = 6
    if (present(print_unit)) lu = print_unit

    if (.not. allocated(grid_function%name)) return
    if (trim(grid_function%name) /= 'spinor_components') return

    call visual_spinor_components_report_selection(cmo, grid_function%spinor_index, lu)

    if (allocated(grid_function%spinor_phase)) then
       write(lu,'(a,a)') '   requested phase mode       : ', trim(grid_function%spinor_phase)
    else
       write(lu,'(a)') '   requested phase mode       : none'
    end if

  end subroutine visual_spinor_components_handle_grid_function


  ! ---------------------------------------------------------------------------
  subroutine visual_spinor_components_evaluate_grid_function(cmo, grid_function, grid, mol, print_unit)
    !
    ! Permanent VISUAL-driver entry point for name='spinor_components'.
    !
    ! At this stage this routine is intentionally a no-op evaluator: it owns
    ! the spinor_components branch in visual_driver, so this grid function is
    ! not passed to the density/scalar evaluator.  The actual grid evaluation
    ! and cube writer will be implemented here in the next patch.
    !
    ! This is preferable to a bare cycle in visual_driver because the control
    ! flow is already final:
    !
    !   if spinor_components:
    !      call visual_spinor_components_evaluate_grid_function(...)
    !   else:
    !      call get_quantity_in_point_refactored(...)
    !
    type(cmo_t), intent(in) :: cmo
    type(grid_function_t), intent(in) :: grid_function
    type(grid_t), intent(in) :: grid
    type(mol_t), intent(in) :: mol
    type(spinor_c1_basis_t) :: basis_l
    type(spinor_c1_basis_t) :: basis_s
    integer, intent(in), optional :: print_unit

    integer :: lu
    integer :: ifsym
    integer :: imo_local
    integer :: n_lao
    integer :: n_sao
    integer :: n_lao_total
    integer :: n_sao_total
    integer :: ic1
    logical :: ok

    integer :: maxopr
    logical :: sym_info

    real(8), allocatable :: c_ao_q_c1(:,:)
    real(8), allocatable :: c_ao_q(:,:)
    real(8), allocatable :: ao_c1_test(:)
    real(8) :: xyz_test(3)
    
    complex(8), allocatable :: c_l_alpha(:)
    complex(8), allocatable :: c_l_beta(:)
    complex(8), allocatable :: c_s_alpha(:)
    complex(8), allocatable :: c_s_beta(:)

    lu = 6
    if (present(print_unit)) lu = print_unit

    if (.not. allocated(grid_function%name)) return
    if (trim(grid_function%name) /= 'spinor_components') return
    write(lu,'(a)') ' VISUAL spinor-components grid-function handler'
    write(lu,'(a)') '   status                     : selected; extracting spinor coefficients'
    
    call visual_spinor_components_report_selection(cmo, grid_function%spinor_index, lu)
    
    if (allocated(grid_function%spinor_phase)) then
       write(lu,'(a,a)') '   requested phase mode       : ', trim(grid_function%spinor_phase)
    else
       write(lu,'(a)') '   requested phase mode       : none'
    end if
    
    call visual_spinor_components_locate_mo(cmo, grid_function%spinor_index, &
         ifsym, imo_local, ok)
    
    if (.not. ok) then
       write(lu,'(a)') '   status                    : invalid spinor index; skip'
       return
    end if
    
    call visual_spinor_components_get_lcsc_dims(cmo, n_lao, n_sao)
    
    write(lu,'(a,i0)') '   LC AOs per fsym            : ', n_lao
    write(lu,'(a,i0)') '   SC AOs per fsym            : ', n_sao
    write(lu,'(a,i0)') '   LC+SC                      : ', n_lao+n_sao
    write(lu,'(a,i0)') '   cmo%n_ao(ifsym)            : ', cmo%n_ao(ifsym)

    n_lao_total = n_lao * cmo%n_fsym
    n_sao_total = n_sao * cmo%n_fsym
   
    maxopr = 0
    call checkpoint_query('/result/symmetry', sym_info)
    
    if (sym_info) then
       call checkpoint_read('/result/symmetry/maxopr', idata=maxopr)
    end if
    
    if (maxopr == 0) then
       call visual_spinor_components_read_c1_basis('/input/aobasis/1/', basis_l)
       call visual_spinor_components_read_c1_basis('/input/aobasis/2/', basis_s)
    else
       call visual_spinor_components_read_c1_basis('/input/aobasis/3/', basis_l)
       call visual_spinor_components_read_c1_basis('/input/aobasis/4/', basis_s)
    end if    
    
    write(lu,'(a,i0)') '   C1 LC basis n_ao           : ', basis_l%n_ao
    write(lu,'(a,i0)') '   C1 SC basis n_ao           : ', basis_s%n_ao
    write(lu,'(a,i0)') '   C1 LC basis n_shells       : ', basis_l%n_shells
    write(lu,'(a,i0)') '   C1 SC basis n_shells       : ', basis_s%n_shells
    write(lu,'(a,i0)') '   C1 LC angular              : ', basis_l%angular
    write(lu,'(a,i0)') '   C1 SC angular              : ', basis_s%angular
    write(lu,'(a,i0)') '   C1 LC generated AO labels  : ', size(basis_l%ao_shell)
    write(lu,'(a,i0)') '   C1 SC generated AO labels  : ', size(basis_s%ao_shell)

    if (basis_l%n_ao /= n_lao_total) then
       write(lu,'(a)') '   ERROR: LC C1 basis n_ao does not match n_lao_total'
       return
    end if
    
    if (basis_s%n_ao /= n_sao_total) then
       write(lu,'(a)') '   ERROR: SC C1 basis n_ao does not match n_sao_total'
       return
    end if


    write(lu,'(a,i0)') '   LC AOs total C1            : ', n_lao_total
    write(lu,'(a,i0)') '   SC AOs total C1            : ', n_sao_total

    if (maxopr == 0) then
       !
       ! Native NOSYM/C1 calculation: cmo%coeffs is already in the
       ! required C1 AO ordering.
       !
       call visual_spinor_components_extract_coefficients(cmo, ifsym, imo_local, &
            c_ao_q_c1)
       ic1 = grid_function%spinor_index
    else
       !
       ! Symmetry-adapted calculation: use the separately stored
       ! C1-expanded orbitals.
       !
       call visual_spinor_components_extract_c1_coefficients(cmo, &
            grid_function%spinor_index, n_lao_total, n_sao_total, &
            c_ao_q_c1, ic1)
    end if

    write(lu,'(a,i0)') '   mapped C1 spinor index     : ', ic1
    write(lu,'(a,i0)') '   C1 AO dimension            : ', size(c_ao_q_c1,1)
    write(lu,'(a,i0)') '   C1 quaternion components   : ', size(c_ao_q_c1,2)
    
    write(lu,'(a,es16.8)') '   C1 coeff q1 min            : ', minval(c_ao_q_c1(:,1))
    write(lu,'(a,es16.8)') '   C1 coeff q1 max            : ', maxval(c_ao_q_c1(:,1))
    write(lu,'(a,es16.8)') '   C1 coeff q2 min            : ', minval(c_ao_q_c1(:,2))
    write(lu,'(a,es16.8)') '   C1 coeff q2 max            : ', maxval(c_ao_q_c1(:,2))
    write(lu,'(a,es16.8)') '   C1 coeff q3 min            : ', minval(c_ao_q_c1(:,3))
    write(lu,'(a,es16.8)') '   C1 coeff q3 max            : ', maxval(c_ao_q_c1(:,3))
    write(lu,'(a,es16.8)') '   C1 coeff q4 min            : ', minval(c_ao_q_c1(:,4))
    write(lu,'(a,es16.8)') '   C1 coeff q4 max            : ', maxval(c_ao_q_c1(:,4))
    
    call visual_spinor_components_split_quaternion(c_ao_q_c1, n_lao_total, n_sao_total, &
         c_l_alpha, c_l_beta, c_s_alpha, c_s_beta)
    
    write(lu,'(a,es16.8)') '   max |L_alpha| C1 coeff     : ', maxval(abs(c_l_alpha))
    write(lu,'(a,es16.8)') '   max |L_beta| C1 coeff      : ', maxval(abs(c_l_beta))
    write(lu,'(a,es16.8)') '   max |S_alpha| C1 coeff     : ', maxval(abs(c_s_alpha))
    write(lu,'(a,es16.8)') '   max |S_beta| C1 coeff      : ', maxval(abs(c_s_beta))

    allocate(ao_c1_test(basis_l%n_ao + basis_s%n_ao))

    xyz_test(1) = 0.0d0
    xyz_test(2) = 0.0d0
    xyz_test(3) = 0.0d0
    
    call visual_spinor_components_eval_c1_ao_at_point(basis_l, basis_s, xyz_test, ao_c1_test)
    
    write(lu,'(a,3f14.8)') '   AO test point              :', xyz_test
    write(lu,'(a,es16.8)') '   max |AO_C1_L|              : ', maxval(abs(ao_c1_test(1:basis_l%n_ao)))
    write(lu,'(a,es16.8)') '   max |AO_C1_S|              : ', maxval(abs(ao_c1_test(basis_l%n_ao+1:)))
    write(lu,'(a,es16.8)') '   first AO_C1_L              : ', ao_c1_test(1)
    write(lu,'(a,es16.8)') '   first AO_C1_S              : ', ao_c1_test(basis_l%n_ao+1)
    
    deallocate(ao_c1_test)
    call visual_spinor_components_write_8cubes_c1(grid_function, grid, mol, &
         basis_l, basis_s, c_l_alpha, c_l_beta, c_s_alpha, c_s_beta, lu)    

    write(lu,'(a)') '   status                     : C1 coefficient extraction and cube writing finished'

  end subroutine visual_spinor_components_evaluate_grid_function


  ! ---------------------------------------------------------------------------
  subroutine visual_spinor_components_report_selection(cmo, imo_global, print_unit)
    !
    ! Report how a user-visible/native VISUAL MO index maps to the internal
    ! fermion-symmetry block and local MO index.
    !
    ! This deliberately uses cmo%n_mo and the native VISUAL ordering, not
    ! eigenvalues_C1/orbitals_C1 ordering.
    !
    type(cmo_t), intent(in) :: cmo
    integer, intent(in) :: imo_global
    integer, intent(in), optional :: print_unit

    integer :: lu
    integer :: ifsym
    integer :: imo_local
    logical :: ok

    lu = 6
    if (present(print_unit)) lu = print_unit

    call visual_spinor_components_locate_mo(cmo, imo_global, ifsym, imo_local, ok)

    write(lu,'(a)') ' VISUAL spinor-components selection'
    write(lu,'(a,i0)') '   requested global MO index  : ', imo_global
    write(lu,'(a,i0)') '   n_fsym                     : ', cmo%n_fsym
    write(lu,'(a,2(i0,1x))') '   n_mo                       : ', cmo%n_mo
    write(lu,'(a,2(i0,1x))') '   n_ao                       : ', cmo%n_ao
    write(lu,'(a,i0)') '   nz                         : ', cmo%nz
    write(lu,'(a,*(1x,i0))') '   n_po                       :', cmo%n_po

    if (.not. ok) then
       write(lu,'(a)') '   status                    : invalid MO index'
       return
    end if

    if (imo_local <= cmo%n_po(ifsym)) then
       write(lu,'(a)') '   WARNING: selected spinor is in negative-energy sector'
    else
       write(lu,'(a)') '   selected sector            : electronic'
    end if

    write(lu,'(a,i0)') '   fermion symmetry block     : ', ifsym
    write(lu,'(a,i0)') '   local MO index in block    : ', imo_local

    if (allocated(cmo%eigenvalues)) then
       if (imo_global >= 1 .and. imo_global <= size(cmo%eigenvalues)) then
          write(lu,'(a,es24.16)') '   energy                     : ', cmo%eigenvalues(imo_global)
       end if
    end if

    if (allocated(cmo%occupations)) then
       if (imo_global >= 1 .and. imo_global <= size(cmo%occupations)) then
          write(lu,'(a,es24.16)') '   occupation                 : ', cmo%occupations(imo_global)
       end if
    end if

    if (allocated(cmo%symmetry)) then
       if (imo_global >= 1 .and. imo_global <= size(cmo%symmetry)) then
          write(lu,'(a,i0)') '   symmetry label             : ', cmo%symmetry(imo_global)
       end if
    end if

    if (allocated(cmo%shell_id)) then
       if (imo_global >= 1 .and. imo_global <= size(cmo%shell_id)) then
          write(lu,'(a,i0)') '   shell id                   : ', cmo%shell_id(imo_global)
       end if
    end if

  end subroutine visual_spinor_components_report_selection


  ! ---------------------------------------------------------------------------
  subroutine visual_spinor_components_locate_mo(cmo, imo_global, ifsym, imo_local, ok)
    !
    ! Convert a VISUAL/native global MO index into:
    !
    !   ifsym     = fermion symmetry block, 1..cmo%n_fsym
    !   imo_local = local MO index inside that block, 1..cmo%n_mo(ifsym)
    !
    ! This is intentionally based only on cmo%n_mo.
    !
    type(cmo_t), intent(in)  :: cmo
    integer, intent(in)  :: imo_global
    integer, intent(out) :: ifsym
    integer, intent(out) :: imo_local
    logical, intent(out) :: ok

    integer :: i
    integer :: first
    integer :: last

    ok = .false.
    ifsym = 0
    imo_local = 0

    if (imo_global < 1) return

    first = 1
    do i = 1, cmo%n_fsym
       last = first + cmo%n_mo(i) - 1

       if (imo_global >= first .and. imo_global <= last) then
          ifsym = i
          imo_local = imo_global - first + 1
          ok = .true.
          return
       end if

       first = last + 1
    end do

  end subroutine visual_spinor_components_locate_mo

  subroutine visual_spinor_components_get_lcsc_dims(cmo, n_lao, n_sao)
    use checkpoint
    implicit none
  
    type(cmo_t), intent(in) :: cmo
    integer, intent(out) :: n_lao
    integer, intent(out) :: n_sao
  
    integer :: maxopr
    integer :: inversion
    integer :: n_lao_total
    integer :: n_sao_total
    logical :: sym_info
    character(len=:), allocatable :: group_lc
    character(len=:), allocatable :: group_sc
  
    maxopr = 0
    inversion = 1
  
    call checkpoint_query('/result/symmetry', sym_info)
  
    if (sym_info) then
       call checkpoint_read('/result/symmetry/maxopr', idata=maxopr)
       call checkpoint_read('/result/symmetry/inversion', idata=inversion)
    end if
  
    if (maxopr == 0) then
       group_lc = '/input/aobasis/1/'
       group_sc = '/input/aobasis/2/'
    else
       group_lc = '/input/aobasis/3/'
       group_sc = '/input/aobasis/4/'
    end if
  
    call checkpoint_read(group_lc//'n_ao', idata=n_lao_total)
    call checkpoint_read(group_sc//'n_ao', idata=n_sao_total)
  
    if (mod(n_lao_total, cmo%n_fsym) /= 0) then
       stop 'LC AO dimension is not divisible by n_fsym'
    end if
  
    if (mod(n_sao_total, cmo%n_fsym) /= 0) then
       stop 'SC AO dimension is not divisible by n_fsym'
    end if
  
    n_lao = n_lao_total / cmo%n_fsym
    n_sao = n_sao_total / cmo%n_fsym
  
  end subroutine visual_spinor_components_get_lcsc_dims

  subroutine visual_spinor_components_extract_coefficients(cmo, ifsym, imo_local, c_ao_q)
    !
    ! Extract one native VISUAL/MO spinor from cmo%coeffs.
    !
    ! Assumed local layout inside each fermion-symmetry block:
    !
    !   coeffs(iao, imo, iq)
    !
    ! with iao fastest.
    !
    type(cmo_t), intent(in) :: cmo
    integer, intent(in) :: ifsym
    integer, intent(in) :: imo_local
    real(8), allocatable, intent(out) :: c_ao_q(:,:)
  
    integer :: nao
    integer :: nmo
    integer :: nz
    integer :: iao
    integer :: iq
    integer :: isym
    integer :: idx0
    integer :: idx
  
    if (ifsym < 1 .or. ifsym > cmo%n_fsym) then
       stop 'visual_spinor_components_extract_coefficients: invalid ifsym'
    end if
  
    nao = cmo%n_ao(ifsym)
    nmo = cmo%n_mo(ifsym)
    nz  = cmo%nz
  
    if (imo_local < 1 .or. imo_local > nmo) then
       stop 'visual_spinor_components_extract_coefficients: invalid local MO index'
    end if
  
    allocate(c_ao_q(nao,nz))
    c_ao_q = 0.0d0
  
    !
    ! Native VISUAL coefficient blocks are stored consecutively by fermion symmetry:
    !
    !   block 1: n_ao(1) * n_mo(1) * nz
    !   block 2: n_ao(2) * n_mo(2) * nz
    !   ...
    !
    idx0 = 1
    do isym = 1, ifsym-1
       idx0 = idx0 + cmo%n_ao(isym) * cmo%n_mo(isym) * nz
    end do
  
    do iq = 1, nz
       do iao = 1, nao
  
          idx = idx0                                           &
              + (iq-1)*nao*nmo                                 &
              + (imo_local-1)*nao                               &
              + (iao-1)
  
          c_ao_q(iao,iq) = cmo%coeffs(idx)
  
       end do
    end do
  
  end subroutine visual_spinor_components_extract_coefficients

  subroutine visual_spinor_components_extract_c1_coefficients(cmo, spinor_index, &
       n_lao_total, n_sao_total, c_ao_q_c1, ic1_out)
    !
    ! Extract one spinor from CHECKPOINT mobasis/orbitals_C1.
    !
    ! The input spinor_index is the native compact VISUAL/MO global index.
    ! It is mapped to the C1 index by
    !
    !   ic1 = (imo_local - 1)*n_fsym + ifsym
    !
    ! where (ifsym, imo_local) are obtained from cmo%n_mo.
    !
    ! For C1 orbitals we assume the layout:
    !
    !   orbitals_C1(iao, imo_c1, iq)
    !
    ! with iao fastest.
    !
    use checkpoint
    implicit none
  
    type(cmo_t), intent(in) :: cmo
    integer, intent(in) :: spinor_index
    integer, intent(in) :: n_lao_total
    integer, intent(in) :: n_sao_total
    real(8), allocatable, intent(out) :: c_ao_q_c1(:,:)
    integer, intent(out) :: ic1_out
  
    integer :: ifsym
    integer :: imo_local
    integer :: ic1
    integer :: nmo_c1
    integer :: nao_c1
    integer :: nz_c1
    integer :: ncoef_c1
    integer :: iao
    integer :: iq
    integer :: idx
    real(8) :: delta_e
    logical :: ok
    logical :: has_c1
    real(8), allocatable :: coeff_c1(:)
  
    call visual_spinor_components_locate_mo(cmo, spinor_index, ifsym, imo_local, ok)
  
    if (.not. ok) then
       stop 'visual_spinor_components_extract_c1_coefficients: invalid compact spinor index'
    end if
  
    !
    ! Compact VISUAL ordering:
    ! maping to C1 spinor
    !
    call visual_spinor_components_find_c1_by_energy(cmo, spinor_index, ic1, delta_e)
    ic1_out = ic1
  
    nao_c1 = n_lao_total + n_sao_total
    nmo_c1 = sum(cmo%n_mo)
    nz_c1  = 4
  
    ncoef_c1 = nao_c1 * nmo_c1 * nz_c1
  
    call checkpoint_query('/result/wavefunctions/scf/mobasis/orbitals_C1', has_c1)
  
    if (.not. has_c1) then
       stop 'visual_spinor_components_extract_c1_coefficients: orbitals_C1 not found on CHECKPOINT'
    end if
  
    if (ic1 < 1 .or. ic1 > nmo_c1) then
       stop 'visual_spinor_components_extract_c1_coefficients: C1 spinor index out of range'
    end if
  
    allocate(coeff_c1(ncoef_c1))
    coeff_c1 = 0.0d0
  
    call checkpoint_read('/result/wavefunctions/scf/mobasis/orbitals_C1', rdata=coeff_c1)
  
    allocate(c_ao_q_c1(nao_c1,nz_c1))
    c_ao_q_c1 = 0.0d0
    do iq = 1, nz_c1
       do iao = 1, nao_c1
    
          !
          ! Correct orbitals_C1 layout:
          !
          idx = (iq-1)*nao_c1*nmo_c1 + &
                (ic1-1)*nao_c1 +       &
                iao
          c_ao_q_c1(iao,iq) = coeff_c1(idx)
    
       end do
    end do

    deallocate(coeff_c1)
  
  end subroutine visual_spinor_components_extract_c1_coefficients


  subroutine visual_spinor_components_split_quaternion(c_ao_q, n_lao, n_sao, &
       c_l_alpha, c_l_beta, c_s_alpha, c_s_beta)
    !
    ! Split AO quaternion coefficients into four complex spinor components.
    !
    ! Current convention, matching the external proof-of-concept:
    !
    !   alpha = q0 + i q1
    !   beta  = -q2 + i q3
    !
    ! This convention should later be cross-checked against DIRAC quaternion
    ! internals, but this is the working convention from the H2+ tests.
    !
    real(8), intent(in) :: c_ao_q(:,:)
    integer, intent(in) :: n_lao
    integer, intent(in) :: n_sao
  
    complex(8), allocatable, intent(out) :: c_l_alpha(:)
    complex(8), allocatable, intent(out) :: c_l_beta(:)
    complex(8), allocatable, intent(out) :: c_s_alpha(:)
    complex(8), allocatable, intent(out) :: c_s_beta(:)
  
    integer :: mu
    integer :: off_s
    integer :: nz
  
    nz = size(c_ao_q,2)
  
    if (nz /= 4) then
       stop 'visual_spinor_components_split_quaternion: currently requires nz=4'
    end if
  
    if (size(c_ao_q,1) /= n_lao + n_sao) then
       stop 'visual_spinor_components_split_quaternion: inconsistent LC/SC dimensions'
    end if
  
    allocate(c_l_alpha(n_lao))
    allocate(c_l_beta(n_lao))
    allocate(c_s_alpha(n_sao))
    allocate(c_s_beta(n_sao))
  
    do mu = 1, n_lao
       c_l_alpha(mu) = dcmplx(c_ao_q(mu,1),  c_ao_q(mu,2))
       c_l_beta(mu)  = dcmplx(-c_ao_q(mu,3), c_ao_q(mu,4))
    end do
  
    off_s = n_lao
  
    do mu = 1, n_sao
       c_s_alpha(mu) = dcmplx(c_ao_q(off_s+mu,1),  c_ao_q(off_s+mu,2))
       c_s_beta(mu)  = dcmplx(-c_ao_q(off_s+mu,3), c_ao_q(off_s+mu,4))
    end do
  
  end subroutine visual_spinor_components_split_quaternion

  subroutine visual_spinor_components_find_c1_by_energy(cmo, spinor_index, ic1, delta_e)
    !
    ! Map native compact VISUAL spinor index to C1 spinor index
    ! by matching eigenvalues_C1 to cmo%eigenvalues(spinor_index).
    !
    use checkpoint
    implicit none
  
    type(cmo_t), intent(in) :: cmo
    integer, intent(in) :: spinor_index
    integer, intent(out) :: ic1
    real(8), intent(out) :: delta_e
  
    integer :: nmo_c1
    integer :: i
    logical :: has_eval_c1
    real(8) :: e_target
    real(8) :: de
    real(8), allocatable :: eval_c1(:)
  
    nmo_c1 = sum(cmo%n_mo)
  
    if (spinor_index < 1 .or. spinor_index > size(cmo%eigenvalues)) then
       stop 'visual_spinor_components_find_c1_by_energy: spinor_index out of range'
    end if
 
    !
    ! In a native C1/NOSYM calculation the compact and C1 orbital
    ! orderings are identical.  No eigenvalues_C1 mapping is needed.
    !
    if (cmo%n_fsym == 1) then
       ic1 = spinor_index
       delta_e = 0.0d0
       return
    end if

    e_target = cmo%eigenvalues(spinor_index)
  
    call checkpoint_query('/result/wavefunctions/scf/mobasis/eigenvalues_C1', has_eval_c1)
  
    if (.not. has_eval_c1) then
       stop 'visual_spinor_components_find_c1_by_energy: eigenvalues_C1 not found'
    end if
  
    allocate(eval_c1(nmo_c1))
    eval_c1 = 0.0d0
  
    call checkpoint_read('/result/wavefunctions/scf/mobasis/eigenvalues_C1', rdata=eval_c1)
  
    ic1 = 1
    delta_e = abs(eval_c1(1) - e_target)
  
    do i = 2, nmo_c1
       de = abs(eval_c1(i) - e_target)
       if (de < delta_e) then
          delta_e = de
          ic1 = i
       end if
    end do
  
    deallocate(eval_c1)
  
  end subroutine visual_spinor_components_find_c1_by_energy

  subroutine visual_spinor_components_open_cube_files(prefix, units)
    character(len=*), intent(in) :: prefix
    integer, intent(out) :: units(8)
  
    open(newunit=units(1), file=trim(prefix)//'_L_alpha_re.cube', status='replace', action='write')
    open(newunit=units(2), file=trim(prefix)//'_L_alpha_im.cube', status='replace', action='write')
    open(newunit=units(3), file=trim(prefix)//'_L_beta_re.cube',  status='replace', action='write')
    open(newunit=units(4), file=trim(prefix)//'_L_beta_im.cube',  status='replace', action='write')
    open(newunit=units(5), file=trim(prefix)//'_S_alpha_re.cube', status='replace', action='write')
    open(newunit=units(6), file=trim(prefix)//'_S_alpha_im.cube', status='replace', action='write')
    open(newunit=units(7), file=trim(prefix)//'_S_beta_re.cube',  status='replace', action='write')
    open(newunit=units(8), file=trim(prefix)//'_S_beta_im.cube',  status='replace', action='write')
  
  end subroutine visual_spinor_components_open_cube_files
  
  subroutine visual_spinor_components_write_cube_headers(units, prefix, grid, mol)
    integer, intent(in) :: units(8)
    character(len=*), intent(in) :: prefix
    type(grid_t), intent(in) :: grid
    type(mol_t), intent(in) :: mol
  
    integer :: iu
    integer :: iatom
    integer :: iz
    integer :: nat
    real(8) :: q
  
    nat = mol%n_atoms
  
    do iu = 1, 8
  
       write(units(iu),'(a)') 'DIRAC VISUAL spinor component cube'
       write(units(iu),'(a)') trim(prefix)
  
       write(units(iu),'(i5,3f12.6)') nat, grid%origin(1), grid%origin(2), grid%origin(3)
  
       write(units(iu),'(i5,3f12.6)') grid%nr_points(1), grid%spacing(1), 0.0d0, 0.0d0
       write(units(iu),'(i5,3f12.6)') grid%nr_points(2), 0.0d0, grid%spacing(2), 0.0d0
       write(units(iu),'(i5,3f12.6)') grid%nr_points(3), 0.0d0, 0.0d0, grid%spacing(3)
  
       do iatom = 1, nat
          iz = nint(mol%nuc_charge(iatom))
          q  = mol%nuc_charge(iatom)
          write(units(iu),'(i5,f12.6,3f12.6)') iz, q, &
               mol%geometry(1,iatom), mol%geometry(2,iatom), mol%geometry(3,iatom)
       end do
  
    end do
  
  end subroutine visual_spinor_components_write_cube_headers
  
  subroutine visual_spinor_components_close_cube_files(units)
    integer, intent(in) :: units(8)
    integer :: i
  
    do i = 1, 8
       close(units(i))
    end do
  
  end subroutine visual_spinor_components_close_cube_files
  
  subroutine visual_spinor_components_get_phase(grid_function, c_l_alpha, phase)
    type(grid_function_t), intent(in) :: grid_function
    complex(8), intent(in) :: c_l_alpha(:)
    complex(8), intent(out) :: phase
  
    integer :: i
    integer :: imax
    real(8) :: amax
  
    phase = dcmplx(1.0d0, 0.0d0)
  
    if (.not. allocated(grid_function%spinor_phase)) return
  
    select case (trim(grid_function%spinor_phase))
  
    case ('none')
       phase = dcmplx(1.0d0, 0.0d0)
  
    case ('i')
       phase = dcmplx(0.0d0, 1.0d0)
  
    case ('-i')
       phase = dcmplx(0.0d0, -1.0d0)
  
    case ('auto_lalpha')
       imax = 1
       amax = abs(c_l_alpha(1))
  
       do i = 2, size(c_l_alpha)
          if (abs(c_l_alpha(i)) > amax) then
             amax = abs(c_l_alpha(i))
             imax = i
          end if
       end do
  
       if (amax > 0.0d0) then
          phase = conjg(c_l_alpha(imax)) / abs(c_l_alpha(imax))
       else
          phase = dcmplx(1.0d0, 0.0d0)
       end if
  
    case default
       phase = dcmplx(1.0d0, 0.0d0)
  
    end select
  
  end subroutine visual_spinor_components_get_phase
  
  subroutine visual_spinor_components_read_c1_basis(group, basis)
    use checkpoint
    implicit none
  
    character(len=*), intent(in) :: group
    type(spinor_c1_basis_t), intent(out) :: basis
  
    integer :: n_exp_total
    integer :: n_coeff_total
    integer :: ishell
    integer :: i
    integer :: off
    real(8), allocatable :: center_flat(:)
  
    call checkpoint_read(group//'n_ao', idata=basis%n_ao)
    call checkpoint_read(group//'n_shells', idata=basis%n_shells)
    call checkpoint_read(group//'angular', idata=basis%angular)
  
    allocate(basis%orbmom(basis%n_shells))
    allocate(basis%n_prim(basis%n_shells))
    allocate(basis%n_cont(basis%n_shells))
    allocate(basis%centers(3,basis%n_shells))
    allocate(basis%shell_exp_offset(basis%n_shells))
    allocate(basis%shell_coeff_offset(basis%n_shells))
  
    call checkpoint_read(group//'orbmom', idata=basis%orbmom)
    call checkpoint_read(group//'n_prim', idata=basis%n_prim)
    call checkpoint_read(group//'n_cont', idata=basis%n_cont)
  
    allocate(center_flat(3*basis%n_shells))
    center_flat = 0.0d0
  
    call checkpoint_read(group//'center', rdata=center_flat)
  
    do ishell = 1, basis%n_shells
       do i = 1, 3
          basis%centers(i,ishell) = center_flat(3*(ishell-1) + i)
       end do
    end do
  
    deallocate(center_flat)
  
    n_exp_total = sum(basis%n_prim)
    n_coeff_total = sum(basis%n_prim * basis%n_cont)
  
    allocate(basis%exponents(n_exp_total))
    allocate(basis%contractions(n_coeff_total))
  
    call checkpoint_read(group//'exponents', rdata=basis%exponents)
    call checkpoint_read(group//'contractions', rdata=basis%contractions)
  
    off = 1
    do ishell = 1, basis%n_shells
       basis%shell_exp_offset(ishell) = off
       off = off + basis%n_prim(ishell)
    end do
  
    off = 1
    do ishell = 1, basis%n_shells
       basis%shell_coeff_offset(ishell) = off
       off = off + basis%n_prim(ishell)*basis%n_cont(ishell)
    end do
  
    call visual_spinor_components_build_cartesian_ao_labels(basis)
  
  end subroutine visual_spinor_components_read_c1_basis

  subroutine visual_spinor_components_build_cartesian_ao_labels(basis)
    implicit none
  
    type(spinor_c1_basis_t), intent(inout) :: basis
  
    integer :: ishell
    integer :: ic
    integer :: iao
    integer :: l
  
    allocate(basis%ao_shell(basis%n_ao))
    allocate(basis%ao_cont(basis%n_ao))
    allocate(basis%ao_lx(basis%n_ao))
    allocate(basis%ao_ly(basis%n_ao))
    allocate(basis%ao_lz(basis%n_ao))
  
    iao = 0
  
    do ishell = 1, basis%n_shells
  
       ! Same convention as the external Python script.
       l = basis%orbmom(ishell) - 1
  
       do ic = 1, basis%n_cont(ishell)
  
          select case (l)
  
          case (0)
             call add_cart_ao(0,0,0)
  
          case (1)
             call add_cart_ao(1,0,0)
             call add_cart_ao(0,1,0)
             call add_cart_ao(0,0,1)
  
          case (2)
             call add_cart_ao(2,0,0)
             call add_cart_ao(0,2,0)
             call add_cart_ao(0,0,2)
             call add_cart_ao(1,1,0)
             call add_cart_ao(1,0,1)
             call add_cart_ao(0,1,1)
  
          case (3)
             call add_cart_ao(3,0,0)
             call add_cart_ao(0,3,0)
             call add_cart_ao(0,0,3)
             call add_cart_ao(2,1,0)
             call add_cart_ao(2,0,1)
             call add_cart_ao(1,2,0)
             call add_cart_ao(0,2,1)
             call add_cart_ao(1,0,2)
             call add_cart_ao(0,1,2)
             call add_cart_ao(1,1,1)
  
          case default
             stop 'visual_spinor_components_build_cartesian_ao_labels: l > 3 not implemented'
          end select
  
       end do
  
    end do
  
    if (iao /= basis%n_ao) then
       stop 'visual_spinor_components_build_cartesian_ao_labels: generated AO count != n_ao'
    end if
  
  contains
  
    subroutine add_cart_ao(lx, ly, lz)
      integer, intent(in) :: lx
      integer, intent(in) :: ly
      integer, intent(in) :: lz
  
      iao = iao + 1
  
      if (iao > basis%n_ao) then
         stop 'visual_spinor_components_build_cartesian_ao_labels: too many generated AOs'
      end if
  
      basis%ao_shell(iao) = ishell
      basis%ao_cont(iao)  = ic
      basis%ao_lx(iao)    = lx
      basis%ao_ly(iao)    = ly
      basis%ao_lz(iao)    = lz
  
    end subroutine add_cart_ao
  
  end subroutine visual_spinor_components_build_cartesian_ao_labels

  subroutine visual_spinor_components_eval_cart_basis_at_point(basis, xyz, values)
    implicit none
  
    type(spinor_c1_basis_t), intent(in) :: basis
    real(8), intent(in) :: xyz(3)
    real(8), intent(out) :: values(:)
  
    integer :: iao
    integer :: ishell
    integer :: ic
    integer :: ip
    integer :: nprim
    integer :: exp_off
    integer :: coeff_off
    integer :: lx
    integer :: ly
    integer :: lz
  
    real(8) :: dx
    real(8) :: dy
    real(8) :: dz
    real(8) :: r2
    real(8) :: alpha
    real(8) :: coeff
    real(8) :: poly
    real(8) :: val
  
    if (size(values) /= basis%n_ao) then
       stop 'visual_spinor_components_eval_cart_basis_at_point: wrong values size'
    end if
  
    values = 0.0d0
  
    do iao = 1, basis%n_ao
  
       ishell = basis%ao_shell(iao)
       ic     = basis%ao_cont(iao)
  
       lx = basis%ao_lx(iao)
       ly = basis%ao_ly(iao)
       lz = basis%ao_lz(iao)
  
       dx = xyz(1) - basis%centers(1,ishell)
       dy = xyz(2) - basis%centers(2,ishell)
       dz = xyz(3) - basis%centers(3,ishell)
  
       r2 = dx*dx + dy*dy + dz*dz
  
       poly = 1.0d0
       if (lx > 0) poly = poly * dx**lx
       if (ly > 0) poly = poly * dy**ly
       if (lz > 0) poly = poly * dz**lz
  
       nprim     = basis%n_prim(ishell)
       exp_off   = basis%shell_exp_offset(ishell)
       coeff_off = basis%shell_coeff_offset(ishell) + (ic-1)*nprim
  
       val = 0.0d0
  
       do ip = 1, nprim
          alpha = basis%exponents(exp_off + ip - 1)
          coeff = basis%contractions(coeff_off + ip - 1)
          val = val + coeff * poly * exp(-alpha*r2)
       end do
  
       values(iao) = val
  
    end do
  
  end subroutine visual_spinor_components_eval_cart_basis_at_point

  subroutine visual_spinor_components_eval_c1_ao_at_point(basis_l, basis_s, xyz, ao_c1)
    implicit none
  
    type(spinor_c1_basis_t), intent(in) :: basis_l
    type(spinor_c1_basis_t), intent(in) :: basis_s
    real(8), intent(in) :: xyz(3)
    real(8), intent(out) :: ao_c1(:)
  
    integer :: n_l
  
    n_l = basis_l%n_ao
  
    if (size(ao_c1) /= basis_l%n_ao + basis_s%n_ao) then
       stop 'visual_spinor_components_eval_c1_ao_at_point: wrong ao_c1 size'
    end if
  
    call visual_spinor_components_eval_cart_basis_at_point(basis_l, xyz, & 
            ao_c1(1:n_l))
    call visual_spinor_components_eval_cart_basis_at_point(basis_s, xyz, &
            ao_c1(n_l+1:n_l+basis_s%n_ao))
  
  end subroutine visual_spinor_components_eval_c1_ao_at_point

  subroutine visual_spinor_components_write_8cubes_c1(grid_function, grid, mol, &
       basis_l, basis_s, c_l_alpha, c_l_beta, c_s_alpha, c_s_beta, lu)
  
    type(grid_function_t), intent(in) :: grid_function
    type(grid_t), intent(in) :: grid
    type(mol_t), intent(in) :: mol
    type(spinor_c1_basis_t), intent(in) :: basis_l
    type(spinor_c1_basis_t), intent(in) :: basis_s
  
    complex(8), intent(in) :: c_l_alpha(:)
    complex(8), intent(in) :: c_l_beta(:)
    complex(8), intent(in) :: c_s_alpha(:)
    complex(8), intent(in) :: c_s_beta(:)
  
    integer, intent(in) :: lu
  
    integer :: units(8)
    integer :: ip
    integer :: iao
    integer :: iv
  
    real(8), allocatable :: ao_c1(:)
    real(8) :: xyz(3)
    real(8) :: vals(8)
  
    character(len=256) :: prefix
  
    complex(8) :: phase
    complex(8) :: l_alpha_val
    complex(8) :: l_beta_val
    complex(8) :: s_alpha_val
    complex(8) :: s_beta_val
  
    write(prefix,'(a,i0)') 'spinor_', grid_function%spinor_index
  
    allocate(ao_c1(basis_l%n_ao + basis_s%n_ao))
  
    call visual_spinor_components_get_phase(grid_function, c_l_alpha, phase)
  
    write(lu,'(a,2f14.8)') '   global phase Re/Im         :', dble(phase), aimag(phase)
  
    call visual_spinor_components_open_cube_files(trim(prefix), units)
    call visual_spinor_components_write_cube_headers(units, trim(prefix), grid, mol)
  
    do ip = 1, grid%nr_points_total
  
       xyz(1) = grid%points(1,ip)
       xyz(2) = grid%points(2,ip)
       xyz(3) = grid%points(3,ip)
  
       call visual_spinor_components_eval_c1_ao_at_point(basis_l, basis_s, xyz, ao_c1)
  
       l_alpha_val = dcmplx(0.0d0, 0.0d0)
       l_beta_val  = dcmplx(0.0d0, 0.0d0)
       s_alpha_val = dcmplx(0.0d0, 0.0d0)
       s_beta_val  = dcmplx(0.0d0, 0.0d0)
  
       do iao = 1, basis_l%n_ao
          l_alpha_val = l_alpha_val + c_l_alpha(iao) * ao_c1(iao)
          l_beta_val  = l_beta_val  + c_l_beta(iao)  * ao_c1(iao)
       end do
  
       do iao = 1, basis_s%n_ao
          s_alpha_val = s_alpha_val + c_s_alpha(iao) * ao_c1(basis_l%n_ao + iao)
          s_beta_val  = s_beta_val  + c_s_beta(iao)  * ao_c1(basis_l%n_ao + iao)
       end do
  
       l_alpha_val = phase * l_alpha_val
       l_beta_val  = phase * l_beta_val
       s_alpha_val = phase * s_alpha_val
       s_beta_val  = phase * s_beta_val
  
       vals(1) = dble(l_alpha_val)
       vals(2) = aimag(l_alpha_val)
       vals(3) = dble(l_beta_val)
       vals(4) = aimag(l_beta_val)
       vals(5) = dble(s_alpha_val)
       vals(6) = aimag(s_alpha_val)
       vals(7) = dble(s_beta_val)
       vals(8) = aimag(s_beta_val)
  
       do iv = 1, 8
          write(units(iv),'(1x,1pe13.5)', advance='no') vals(iv)
          if (mod(ip,6) == 0) write(units(iv),*)
       end do
  
    end do
  
    do iv = 1, 8
       write(units(iv),*)
    end do
  
    call visual_spinor_components_close_cube_files(units)
  
    deallocate(ao_c1)
  
    write(lu,'(a)') '   status                     : wrote 8 C1 spinor component cube files'
  
  end subroutine visual_spinor_components_write_8cubes_c1

end module visual_spinor_components
