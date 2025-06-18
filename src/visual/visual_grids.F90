module visual_grids

  use visual_datatypes, only: grid_t, mol_t
  use visual_io
  use visual_utils

  implicit none

  public get_grid_info
  public get_grid_points

  interface generate_grid_rectilinear
    module procedure generate_grid_rectilinear
  end interface

  private

contains

! ------------------------------------------------------------------------------
  subroutine get_grid_info(grid_id, grid, mol, dmat)
use dirac_dmat_factory
use num_grid_gen
    integer, intent(in):: grid_id
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    type(dmat_t), intent(in):: dmat
    character(len = 150):: errmsg
    logical:: tobe
!   get additional information about the grid,
!   (which was not set up via user input in 
!   input/input_reader_sections.F90)

    if (grid%status == 0) return

    if (grid%input  ==  1 .or. grid%input  ==  -1) then
!
!     try reading grid information from file
!
      if (grid%input  ==  -1) then
!       try reading grid information from "numerical_grid",
!       if not exists, then first regenerate it
        inquire(file='numerical_grid',exist=tobe)
        if (.not. tobe) then
          write(*, *) 'VISUAL: numerical_grid will be regenerated'
          call generate_num_grid(dmat%elements)
        end if
      end if
      call get_grid_info_from_file(grid)
      if (grid%id == 0) grid%id = grid_id

!
!     in other cases - set up what is needed to generate grid points
!
    else if (grid%input==0 .or. grid%input==5 .or. grid%input==6) then
      ! 3d grid

      if (grid_id /= grid%id) then
        write (errmsg, *) 'ERROR: problem with id in get_grid_info:',&
                          ' grid%id=', grid%id,&
                          ' grid_id=', grid_id
        call visual_quit(errmsg)
      end if

      if (grid%typ == 1) then
        ! rectilinear uniform 3d grid
        call get_3d_grid_essentials(grid, mol)

        ! center grid around molecule,
        ! modify grid extents and grid%origin (but do not move molecule)
        ! todo - add a keyword + an option to control grid extents
        call center_mol_in_3d_grid(grid)

      else if (grid%typ == 2) then
        ! radial 3d grid
        call get_3d_radial_grid_essentials(grid)
      end if

    else if (grid%input  ==  3) then
      ! line
      call get_line_grid_essentials(grid)

    end if

    ! we should now have all information to work on the grid
    grid%status = 0

  end subroutine get_grid_info


! ------------------------------------------------------------------------------
  subroutine get_grid_points(grid)
    type(grid_t), intent(inout):: grid

    if (grid%status /= 0) return
    if (.not. grid%allocated_points()) return
    if (grid%requires_weights() .and..not. grid%allocated_weights()) return
    if (grid%input == 2) return ! grid points should already be read from input

    if (grid%input == 1 .or. grid%input == -1) then
      call get_grid_points_from_file(grid)
    else
      call generate_grid_points(grid)
    end if
  end subroutine


! ------------------------------------------------------------------------------
  subroutine generate_grid_radial(grid)
use num_grid_cfg
    type(grid_t), intent(inout):: grid
    real(8),allocatable:: tmp(:,:)
    real(8):: rad, radial_spacing
    integer:: nleb, i, m, n, istep, ipoint, grid_angint
    character(len = 150):: errmsg

    ! Get angular points and weights
    allocate(tmp(grid%radial_nr_points_angular, 4))
    tmp=0.0d0
    CALL LEBEDEV(grid%radial_lebedev_l,&
                 tmp(1,1),tmp(1,2),tmp(1,3),tmp(1, 4),n)
    CALL DSCAL(grid%radial_nr_points_angular, 4*ACOS(-1.0D0), tmp(1,4), 1)

    call generate_grid_radial_points_weights(grid,tmp)

    deallocate(tmp)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine generate_grid_radial_points_weights(grid,tmp)
    type(grid_t), intent(inout):: grid
    real(8),intent(in):: tmp(grid%radial_nr_points_angular,*)
    real(8):: rad, radial_spacing
    integer:: i, istep, ipoint

    ! Loop over radial points
    if (grid%radial_nr_steps > 1) then
      radial_spacing = grid%radial_length/(grid%radial_nr_steps-1)
    else
      radial_spacing = grid%radial_length
    end if

    ipoint = 1
    do istep = 0, grid%radial_nr_steps-1
      rad = istep*radial_spacing
      do i = 1, grid%radial_nr_points_angular
        grid%points(1:3,ipoint) = grid%origin(1:3) + rad*tmp(i, 1:3)
        ! grid%weights will be used for angular integration:
        if (grid%allocated_weights()) then
          grid%weights(ipoint) = tmp(i, 4)*rad*rad
        end if
        ipoint = ipoint + 1
      end do
      if (grid%allocated_radial_coord()) then
        grid%radial_coord(istep+1) = rad
      end if
    end do

  end subroutine


! ------------------------------------------------------------------------------
  subroutine generate_grid_rectilinear(grid)
! rectilinear equidistant
    type(grid_t), intent(inout):: grid
    integer:: i,j,k, ipoint, N(grid%get_ndim())
    integer:: ix,iy 
    real(8):: v(3), vx(3), vy(3), vz(3), p(3), p_off(3)
    real(8):: vx_norm, vy_norm, vz_norm, s, f
    real(8),allocatable::x(:),w(:)
    character(len = 150):: errmsg

    if (grid%typ /= 1) return

    if (grid%input == 0 .or. grid%input == 5) then
!     3d
      ipoint = 1
      do i = 1, grid%nr_points(1)
        v(1) = grid%origin(1) + (i - 1)*grid%spacing(1)
        do j = 1, grid%nr_points(2)
          v(2) = grid%origin(2) + (j - 1)*grid%spacing(2)
          do k = 1, grid%nr_points(3)
            v(3) = grid%origin(3) + (k - 1)*grid%spacing(3)
            grid%points(1:3, ipoint) = v(1:3)
            ipoint = ipoint + 1
          end do
        end do
      end do
    else if (grid%input == 4) then
!     2d
      vx(1:3) = (grid%extents(1,1:3) - grid%extents(0,1:3))/(grid%nr_points(1)-1)
      vy(1:3) = (grid%extents(2,1:3) - grid%extents(0,1:3))/(grid%nr_points(2)-1)
      vz(1) = vx(2)*vy(3) - vx(3)*vy(2)
      vz(2) = vx(3)*vy(1) - vx(1)*vy(3)
      vz(3) = vx(1)*vy(2) - vx(2)*vy(1)
      vx_norm = dsqrt(vx(1)*vx(1) + vx(2)*vx(2) + vx(3)*vx(3))
      vy_norm = dsqrt(vy(1)*vy(1) + vy(2)*vy(2) + vy(3)*vy(3))
      vz_norm = dsqrt(vz(1)*vz(1) + vz(2)*vz(2) + vz(3)*vz(3))
      ipoint = 1
      do i = 1, grid%nr_points(1)
        do j = 1, grid%nr_points(2)
          grid%points(1:3, ipoint) = grid%origin(1:3) + (i-1)*vx(1:3) + (j-1)*vy
          ipoint = ipoint + 1
        end do
      end do
    else if (grid%input == -4) then
    ! check that angle between 1-2 and 1-3 is 90 degs
      s = vx(1)*vy(1) + vx(2)*vy(2) + vx(3)*vy(3)
      if (s > tiny(0.0d0)) then
        call visual_quit('2d integration: lines connecting points 1-2 and 1-3'//&
                         ' do not enclose 90 degrees')
      end if
      vx(1:3) = grid%extents(1,1:3) - grid%extents(0,1:3)
      vy(1:3) = grid%extents(2,1:3) - grid%extents(0,1:3)
      vx_norm = dsqrt(vx(1)*vx(1) + vx(2)*vx(2) + vx(3)*vx(3))
      vy_norm = dsqrt(vy(1)*vy(1) + vy(2)*vy(2) + vy(3)*vy(3))
      allocate(x(grid%integration_order_2d))
      allocate(w(grid%integration_order_2d))
      ! nodes, weights in [-1,1]
      call lobatto_1d_grid(grid%integration_order_2d, x, w)
      ! nodes, weights in an interval set by grid%extents
      ipoint = 0
      f = 4.0d0*grid%nr_points(1)*grid%nr_points(2)
      do ix = 0, grid%nr_points(1)-1
         do iy = 0, grid%nr_points(2)-1
            p_off = grid%origin + ix*vx/grid%nr_points(1) + iy*vy/grid%nr_points(2)
            do i = 1, grid%integration_order_2d
               do j = 1, grid%integration_order_2d
                  ipoint = ipoint+1
                  p = p_off + 0.5d0*(x(i) + 1.0d0)*vx/grid%nr_points(1)  &
                            + 0.5d0*(x(j) + 1.0d0)*vy/grid%nr_points(2)
                  grid%points(1:3, ipoint) = p(1:3)
                  grid%weights(ipoint) = &
                  w(i)*w(j)*vx_norm*vy_norm/f
               end do
            end do
         end do
      end do
      deallocate(x)
      deallocate(w)
    else if (grid%input == 3) then
!     line
      ipoint = 1
      do i = 1, grid%nr_points(1)
        v(1:3) = grid%origin(1:3) + (i - 1)*grid%spacing(1:3)
        grid%points(1:3, ipoint) = v(1:3)
        ipoint = ipoint + 1
      end do
    else
      write (errmsg, *) 'ERROR: unsupported grid input '
      call visual_quit(errmsg)
    end if

  end subroutine

! ------------------------------------------------------------------------------
  subroutine grid_file_txt_write(file_unit, float_format, nr_points, x, y, z)
     integer, intent(in):: file_unit
     character(*), intent(in):: float_format
     integer, intent(in):: nr_points
     real(8), intent(in):: x(:), y(:), z(:)
     character(len = 8):: xc, yc, zc
     integer             :: i

     do i = 1, nr_points
        write(file_unit, float_format) x(i), y(i), z(i)
     enddo
  end subroutine


! ------------------------------------------------------------------------------
  subroutine grid_file_txt_read(file_unit, float_format, nr_points, x, y, z)
     real(8), intent(inout), dimension(:):: x, y, z
     integer, intent(in)  :: file_unit
     character(*), intent(in):: float_format
     integer, intent(in)  :: nr_points
     integer              :: i

     rewind(file_unit)

     i = 1
     do
        read(file_unit, *, end = 98) x(i), y(i), z(i)
        i = i+1
     end do
 98  continue

  end subroutine


! ------------------------------------------------------------------------------
  subroutine grid_origin_extents_from_mol(mol, grid)
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    real(8):: extents(0:1, 3)
    real(8):: atoms_min(1:3), atoms_max(1:3)
    integer::i, j

    atoms_min(1:3) = mol%geometry(1:3, 1)
    atoms_max(1:3) = mol%geometry(1:3, 1)

    do i = 1, mol%n_atoms
      do j = 1, 3
        atoms_min(j) = min(mol%geometry(j, i), atoms_min(j))
        atoms_max(j) = max(mol%geometry(j, i), atoms_max(j))
      end do
    end do

    ! add space around
    do i = 1, 3
      atoms_min(i) = atoms_min(i) - grid%margin
      atoms_max(i) = atoms_max(i) + grid%margin
      extents(0, i) = atoms_min(i)
      extents(1, i) = atoms_max(i)
    end do

    ! update grid info
    grid%origin(1:3) = extents(0, 1:3)
    grid%origin_known = .true.
    grid%extents(0, 1:3) = extents(0, 1:3)
    grid%extents(1, 1:3) = extents(1, 1:3)
    grid%extents_known = .true.

  end subroutine grid_origin_extents_from_mol


! ------------------------------------------------------------------------------
  subroutine grid_3d_npoints_from_spacing(grid,i)
    type(grid_t), intent(inout):: grid
    integer,intent(in):: i

    if (.not.grid%extents_known) return
    if (.not.grid%spacing_known(i)) return

    grid%nr_points(i) = int((grid%extents(1,i) - grid%extents(0,i))/grid%spacing(i)) + 1
    grid%nr_points_known(i) = .true.

  end subroutine


! ------------------------------------------------------------------------------
  subroutine grid_3d_spacing_from_npoints(grid,i)
    type(grid_t), intent(inout):: grid
    integer,intent(in):: i
    character(len = 150):: errmsg

    grid%spacing(i) = 0.0d0
    if (grid%nr_points_known(i)) then
      if (grid%nr_points(i) > 1) then
        grid%spacing(i) = (grid%extents(1,i) - grid%extents(0, i))/(grid%nr_points(i) - 1)
      else
        grid%spacing(i) = grid%extents(1,i) - grid%extents(0,i)
      end if
      grid%spacing_known(i) = .true.
    else
      write (errmsg, *) 'ERROR in grid_3d_spacing_from_npoints'
      call visual_quit(errmsg)
    end if
  end subroutine

! ------------------------------------------------------------------------------
  subroutine grid_line_spacing_from_npoints(grid)
    type(grid_t), intent(inout):: grid
    integer:: i
    character(len = 150):: errmsg

    do i = 1, grid%get_ndim()
      if (grid%spacing_known(i)) return
      grid%spacing(i) = 0.0d0
      if (grid%nr_points_total_known) then
        if (grid%nr_points_total > 1) then
          grid%spacing(i) = (grid%extents(1,i) - grid%extents(0, i))/(grid%nr_points_total - 1)
        else
          grid%spacing(i) = grid%extents(1,i) - grid%extents(0,i)
        end if
        grid%spacing_known(i) = .true.
      else
        write (errmsg, *) 'ERROR in grid_3d_spacing_from_npoints'
        call visual_quit(errmsg)
      end if
    end do
  end subroutine grid_line_spacing_from_npoints


! ------------------------------------------------------------------------------
  subroutine center_mol_in_3d_grid(grid)
#include "orgcom.h"
    type(grid_t), intent(inout):: grid
    real(8):: extents(0:1, 3)
    real(8):: cm_0, cm_1
    integer:: i

    do i = 1, 3
      ! get actual grid extents
      extents(0, i) = grid%origin(i)
      extents(1, i) = grid%origin(i) + (grid%nr_points(i) - 1)*grid%spacing(i)
      cm_0 = 0.5d0*(extents(1, i) + extents(0, i))
      !cm_1 = cmxyz(i)  ! this is center-of-mass from "orgcom.h" in a.u.
      cm_1 = gagorg(i) ! this is a "magnetic gauge origin" (from "orgcom.h" in a.u.), was used in old DIRAC-VISUAL code 
      extents(0, i) = extents(0, i) + (cm_1 - cm_0)
      extents(1, i) = extents(1, i) + (cm_1 - cm_0)
    end do
    grid%origin(1:3) = extents(0, 1:3)
    grid%origin_known = .true.
  end subroutine center_mol_in_3d_grid


! ------------------------------------------------------------------------------
  subroutine grid_get_nr_points_total(grid)
    type(grid_t), intent(inout):: grid
    integer:: i

    if (grid%nr_points_total_known) return
    grid%nr_points_total = 1
    do i = 1, 3
      if (grid%nr_points(i) > 0) then
        grid%nr_points_total = grid%nr_points_total * grid%nr_points(i)
      end if
    end do
    grid%nr_points_total_known = .true.
  end subroutine grid_get_nr_points_total


! ------------------------------------------------------------------------------
  subroutine get_3d_grid_essentials(grid, mol)
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    integer:: i

    if (grid%typ /= 1) return

    call grid_origin_extents_from_mol(mol, grid)

    do i = 1, grid%get_ndim()
      if (.not. grid%spacing_known(i)) call grid_3d_spacing_from_npoints(grid,i)
      if (.not. grid%nr_points_known(i)) call grid_3d_npoints_from_spacing(grid,i)
    end do

    if (.not. grid%nr_points_total_known) call grid_get_nr_points_total(grid)
  end subroutine get_3d_grid_essentials


! ------------------------------------------------------------------------------
  subroutine get_3d_radial_grid_essentials(grid)
use num_grid_cfg
    type(grid_t), intent(inout):: grid
    integer:: m, nleb

    if (grid%typ /= 2) return

    m = NLEB(num_grid_cfg_angint)

    write(*, '(/,4x, a)')  '** Radial distribution **'
    write(*, '(6x, A, I3)')  '   Lebedev angular grid L = ',num_grid_cfg_angint
    write(*, '(6x, A, I10)') '   Angular points:   ',m

    grid%radial_lebedev_l = num_grid_cfg_angint
    grid%radial_nr_points_angular = m
    if (.not. grid%nr_points_total_known) then
      grid%nr_points_total = m*grid%radial_nr_steps
      grid%nr_points_total_known = .true.
    end if

  end subroutine get_3d_radial_grid_essentials


! ------------------------------------------------------------------------------
  subroutine get_line_grid_essentials(grid)
    type(grid_t), intent(inout):: grid
    real(8):: origin(3)
    real(8):: extents(0:1, 3)
    integer:: i

    if (grid%typ /= 1) return
    call grid_line_spacing_from_npoints(grid)

  end subroutine get_line_grid_essentials


! ------------------------------------------------------------------------------
  subroutine generate_grid_points(grid)
    type(grid_t), intent(inout):: grid
    integer :: i
    character(len = 150):: errmsg

    select case (grid%typ)
      case (1)
        call generate_grid_rectilinear(grid)
      case (2)
        call generate_grid_radial(grid)
      case default
        write (errmsg, *) 'ERROR (generate_grid_points): unsupported grid type.'
        call visual_quit(errmsg)
    end select
  end subroutine

! ---- copied from old VISUAL:

  subroutine legendre_polynomial(p, n, nmax, x)
    integer,intent(in)::n,nmax
    real(8),intent(in)::x
    real(8),intent(out)::p(0:nmax)
    integer::i
    p(0) = 1
    p(1) = x
    do i = 1, n-1
      p(i+1) = (2*i + 1)*x*p(i) - i*p(i-1)
      p(i+1) = p(i+1)/(i+1)
    enddo
  end subroutine

  subroutine legendre_polynomial_derivatives(p, dp, n, nmax, x)
    integer,intent(in)::n, nmax
    real(8),intent(in)::x
    real(8),intent(out)::p(0:nmax),dp(nmax)
    integer::i
    do i = 1, n
      dp(i) = (i*x*p(i) - i*p(i-1))/(x*x - 1.0d0)
    enddo
  end subroutine

  subroutine lobatto_1d_grid(n, x, w)
    integer,intent(in)::n
    real(8),intent(out)::w(*)
    real(8),intent(out)::x(*)
    integer,parameter::nmax = 1000
    real(8),parameter::threshold = 1.0d-12
    real(8)::p(0:nmax),dp(nmax)
    real(8)::x_old,y_old,x_new,y_new,step
    integer::i

    x_old = -1.0d0+threshold
    call legendre_polynomial(p, n, nmax, x_old)
    call legendre_polynomial_derivatives(p, dp, n, nmax, x_old)
    do i = 2, n-1
      step  = 1.0d-4
      do while(dabs(step) .gt. threshold)
        y_old = dp(n-1)
        x_new = x_old+step
        call legendre_polynomial(p, n, nmax, x_new)
        call legendre_polynomial_derivatives(p, dp, n, nmax, x_new)
        y_new = dp(n-1)
        if((y_old/y_new) .lt. 0.0d0) then
          step = -step/10.0d0
        endif
        x_old = x_new
      enddo
      x(i) = x_new
      w(i) = n*(n-1)*p(n-1)*p(n-1)
      w(i) = 2.0d0/w(i)
    enddo
!   end points
    x(1) = -1.0d0
    x(n) =  1.0d0
    w(1) = 2.0d0/(n*(n-1.0d0))
    w(n) = w(1)
  end subroutine

  subroutine legendre_1d_grid(n, x, w)
    integer,intent(in)::n
    real(8),intent(inout)::w(*)
    real(8),intent(inout)::x(*)
    integer,parameter::nmax = 1000
    real(8),parameter::threshold = 1.0d-12
    real(8)::p(0:nmax),dp(nmax)
    real(8)::x_old,y_old,x_new,y_new,step
    integer::i
    x_old = -1.0d0
    call legendre_polynomial(p, n, nmax, x_old)
    do i = 1, n
      step  = 1.0d-4
      do while(dabs(step) .gt. threshold)
        y_old = p(n)
        x_new = x_old+step
        call legendre_polynomial(p, n, nmax, x_new)
        y_new = p(n)
        if((y_old/y_new) .lt. 0.0d0) then
          step = -step/10.0d0
        endif
        x_old = x_new
      enddo
      x(i) = x_new
      call legendre_polynomial_derivatives(p, dp, n, nmax, x(i))
      w(i) = (1.0d0-x(i)*x(i))*dp(n)*dp(n)
      w(i) = 2.0d0/w(i)
    enddo
  end subroutine

end module
