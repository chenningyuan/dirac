module dft_density

!gosia:
!!!!!#define debug_me

  use dft_cb_interface

  implicit none

  public get_n
  public get_gn
  public get_gn_nonhermitian
  public get_s
  public get_gs
  public get_gs_nonhermitian
  public lda_distribute
  public lda_london_distribute !added
  public gga_distribute
  public lda_distribute_qi
  public gga_distribute_qi
  public gga_sym

  private

  integer :: jxyz_start = 1 !noncollinear
! integer :: jxyz_start = 3 !collinear

contains

  subroutine gga_sym(nr_mat, &
                     mat,    &
                     isym,   &
                     ih)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: nr_mat
    real(8), intent(inout) :: mat(nr_ao_total, nr_ao_total, nr_qu, nr_mat)
    integer, intent(in)    :: isym(nr_mat)
    integer, intent(in)    :: ih(nr_mat)
!   ----------------------------------------------------------------------------
    integer                :: irep, iblock, imat, jxyz
    integer                :: i, j, ii, jj
    integer                :: st1, st2, nr1, nr2
    logical                :: skip(0:3)
    real(8)                :: average
!   ----------------------------------------------------------------------------

    do imat = 1, nr_mat

      skip = .false.
      if (ih(imat) ==  1) skip(1:3) = .true.
      if (ih(imat) == -1) skip(0:0) = .true.

      irep = isym(imat) - 1

      do jxyz = 0, 3
        if (.not. skip(jxyz)) then

          do iblock = 1, nr_density_blocks(jxyz, irep)
          
            nr1 = nr_ao(1, iblock, jxyz, irep)
            st1 = st_ao(1, iblock, jxyz, irep)
            nr2 = nr_ao(2, iblock, jxyz, irep)
            st2 = st_ao(2, iblock, jxyz, irep)

            if (.not. (st2 > st1)) then

              do i = 1, nr1
                do j = 1, nr2

                  ii = st1 + i - 1
                  jj = st2 + j - 1

                  average = mat(ii, jj, iq(jxyz), imat) &
                          + mat(jj, ii, iq(jxyz), imat)
                  average = 0.5d0*average

                  mat(ii, jj, iq(jxyz), imat) = average
                  mat(jj, ii, iq(jxyz), imat) = average

                end do
              end do
            end if
          end do
        end if
      end do
    end do

  end subroutine

  subroutine get_n(n, irep, dmat, buffer, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: n
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    n = 0.0d0

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   nr_ao_total,              &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        n = n + ddot(nr1, ao(st1), 1, buffer, 1)

      end do
    end do

  end subroutine

  subroutine get_s(s, irep, dmat, buffer, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: s(3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    s = 0.0d0

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   nr_ao_total,              &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        s(jxyz) = s(jxyz) + ddot(nr1, ao(st1), 1, buffer, 1)

      end do
    end do

  end subroutine

  subroutine get_gn(gn,     &
                    irep,   &
                    dmat,   &
                    buffer, &
                    ao,     &
                    ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
    real(8), intent(in)  :: ao1(nr_ao_total, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    gn = 0.0d0

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   nr_ao_total,              &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        do ixyz = 1, 3
          gn(ixyz) = gn(ixyz) + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)
        end do

      end do
    end do

    gn = 2.0d0*gn

  end subroutine

  subroutine get_gn_nonhermitian(gn,     &
                                 irep,   &
                                 dmat,   &
                                 buffer, &
                                 ao,     &
                                 ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
    real(8), intent(in)  :: ao1(nr_ao_total, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    gn = 0.0d0

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do ixyz = 1, 3

          call dgemv('N',                      &
                     nr1,                      &
                     nr2,                      &
                     1.0d0,                    &
                     dmat(st1, st2, iq(jxyz)), &
                     nr_ao_total,              &
                     ao(st2),                  &
                     1,                        &
                     0.0d0,                    &
                     buffer,                   &
                     1)

          gn(ixyz) = gn(ixyz) + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)

          call dgemv('N',                      &
                     nr1,                      &
                     nr2,                      &
                     1.0d0,                    &
                     dmat(st1, st2, iq(jxyz)), &
                     nr_ao_total,              &
                     ao1(st2, ixyz),           &
                     1,                        &
                     0.0d0,                    &
                     buffer,                   &
                     1)

          gn(ixyz) = gn(ixyz) + ddot(nr1, ao(st1), 1, buffer, 1)

        end do
      end do
    end do

  end subroutine

  subroutine get_gs(gs,     &
                    irep,   &
                    dmat,   &
                    buffer, &
                    ao,     &
                    ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gs(3, 3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
    real(8), intent(in)  :: ao1(nr_ao_total, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    gs = 0.0d0

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   nr_ao_total,              &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        do ixyz = 1, 3
          gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                         + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)
        end do

      end do
    end do

    gs = 2.0d0*gs

  end subroutine

  subroutine get_gs_nonhermitian(gs,     &
                                 irep,   &
                                 dmat,   &
                                 buffer, &
                                 ao,     &
                                 ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gs(3, 3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
    real(8), intent(in)  :: ao1(nr_ao_total, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    gs = 0.0d0

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do ixyz = 1, 3

          call dgemv('N',                      &
                     nr1,                      &
                     nr2,                      &
                     1.0d0,                    &
                     dmat(st1, st2, iq(jxyz)), &
                     nr_ao_total,              &
                     ao(st2),                  &
                     1,                        &
                     0.0d0,                    &
                     buffer,                   &
                     1)

          gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                         + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)

          call dgemv('N',                      &
                     nr1,                      &
                     nr2,                      &
                     1.0d0,                    &
                     dmat(st1, st2, iq(jxyz)), &
                     nr_ao_total,              &
                     ao1(st2, ixyz),           &
                     1,                        &
                     0.0d0,                    &
                     buffer,                   &
                     1)

          gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                         + ddot(nr1, ao(st1), 1, buffer, 1)

        end do
      end do
    end do

  end subroutine

  subroutine lda_distribute(s, irep, fmat, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(nr_ao_total)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
    real(8)                :: f
!   ----------------------------------------------------------------------------
#ifdef debug_me
   write(*, *) 'debug lda_distribute: nr1, st1, nr2, st2, &
  &             ii, ao(ii), f, fmat: '      
#endif

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do i = 1, nr2

          ii = st2 + i - 1

          f = s*ao(ii)

          call daxpy(nr1,                     &
                     f,                       &
                     ao(st1),                 &
                     1,                       &
                     fmat(st1, ii, iq(jxyz)), &
                     1)
#ifdef debug_me
          write(*, *) nr1, st1, nr2, st2, ii, ao(ii), &
                      f, fmat(st1, ii, iq(jxyz))
#endif
        end do
      end do
    end do

  end subroutine

 
! based on lda_distribute
! it calculates the exchange-correlation contribution to KS matrix
! differentiated wrt magnetic field
  subroutine lda_london_distribute(s, irep, fmat,       &
                                   ao,                  &
                                   ao_b1,               &
                                   cx, cy, cz,          &
                                   origx, origy, origz)   

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(nr_ao_total)
    real(8), intent(in)    :: ao_b1(nr_ao_total*3)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz, jj, j
    integer                :: a, b, c
    real(8)                :: f, f_left, f_right
    real(8)                :: cx, cy, cz, origx, origy, origz
    real(8)                :: r_current(3), r_origin(3), r_co(3)
    real(8)                :: aob1(nr_ao_total, 3)
!   ----------------------------------------------------------------------------

! rename:
    r_current(1) = cx
    r_current(2) = cy
    r_current(3) = cz
    r_origin( 1) = origx
    r_origin( 2) = origy
    r_origin( 3) = origz
    aob1(:,   1) = ao_b1(1:nr_ao_total)                 !d/dBx
    aob1(:,   2) = ao_b1(nr_ao_total+1:2*nr_ao_total)   !d/dBy
    aob1(:,   3) = ao_b1(2*nr_ao_total+1:3*nr_ao_total) !d/dBz 
    r_co(     1) = r_current(1) - r_origin(1) 
    r_co(     2) = r_current(2) - r_origin(2)
    r_co(     3) = r_current(3) - r_origin(3)

!#ifdef debug_me
!    !write(*, *) 'gauge origin: ', r_origin(1), r_origin(2), r_origin(3)
!    write(*, *) 'debug lda_london_distribute.1.'
!    write(*, *) 'grid point: ', r_current(1), r_current(2), r_current(3)
!#endif
    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)
        
        do i = 1, nr2
          do j = 1, nr2 ! (will be removed)

            ii = st2 + i - 1
            jj = st2 + j - 1 ! (will be removed)
!#ifdef debug_me
!            write(*, *) 'debug lda_london_distribute.2.    &
!          &             jj, ao(jj), aob1 ..x, ..y, ..z: '
!            write(*, *) jj, ao(jj), aob1(jj, 1),           &
!                        aob1(jj, 2), aob1(jj, 3)
!#endif

! if we want to do it in similar way to dalton, then:

            do a = 1, 3
              select case (a)
                case (1)
                  b = 2; c = 3
                case (2)
                  b = 3; c = 1
                case (3)
                  b = 1; c = 2
              end select

              f_left  = (r_co(b)*aob1(ii, c) - r_co(c)*aob1(ii, b))*ao(jj)
              f_right = (r_co(b)*aob1(jj, c) - r_co(c)*aob1(jj, b))*ao(ii)
              f = (f_left - f_right)
              !
              !write(*, *) 'f_left, f_right, f: ', f_left, f_right, f
              !f_left  = (r_co(b)*aob1(ii, c) - r_co(c)*aob1(ii, b))*s*ao(jj)
              !f_right = ao(ii)*(r_co(b)*s*aob1(jj, c) - r_co(c)*s*aob1(jj, b))
              !f = f_left - f_right
!#ifdef debug_me              
!              write(*, *) 'debug lda_london_distribute.3.       &
!           &   f_left, f_right, f: '
!              write(*, *) f_left, f_right, f
!#endif

! now it is done the same as in dalton...
! but one do-loop and daxpy looks much much better 
! and it can be done (as is in lda_distribute)
! to rewrite later...

            fmat(st1, jj, iq(jxyz)) = fmat(st1, jj, iq(jxyz)) + s*f
!#ifdef debug_me
!            write(*, *) 'debug lda_london_distribute.4.         &
!          & fmat(st1, jj, iq(jxyz))' 
!            write(*, *) fmat(st1, jj, iq(jxyz))
!#endif            
            end do !a = 1, 3

! still some doubts:
! 1. we're missing imaginary unit in multiplication above(and so is dalton)
!    
! 2. in general ao can be complex and here we want to calculate
!    (ao)^{\dag} * s * ao,
!    so ao^{\dag} for complex ao looks to be more complicated
!    (f.ex. zdotc in blas)
!    so if we do not do this conjugation, probably it is incorporated 
!    in symmetries (these iq(jxyz)) or maybe not important here...

          end do ! =1, nr2 (will be removed)
        end do
      end do
    end do

  end subroutine

  subroutine gga_distribute(s, v, irep, fmat, ao, ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    real(8), intent(in)    :: v(3)
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(nr_ao_total)
    real(8), intent(in)    :: ao1(nr_ao_total, 3)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
    real(8)                :: f
!   ----------------------------------------------------------------------------

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do i = 1, nr2

          ii = st2 + i - 1

          f =       s   *ao (ii)    &
            + 2.0d0*v(1)*ao1(ii, 1) &
            + 2.0d0*v(2)*ao1(ii, 2) &
            + 2.0d0*v(3)*ao1(ii, 3)

          call daxpy(nr1,                     &
                     f,                       &
                     ao(st1),                 &
                     1,                       &
                     fmat(st1, ii, iq(jxyz)), &
                     1)

        end do
      end do
    end do

  end subroutine

  subroutine lda_distribute_qi(s, irep, fmat, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s(3)
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(nr_ao_total)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
    real(8)                :: f
!   ----------------------------------------------------------------------------

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do i = 1, nr2

          ii = st2 + i - 1

          f = s(jxyz)*ao(ii)

          call daxpy(nr1,                     &
                     f,                       &
                     ao(st1),                 &
                     1,                       &
                     fmat(st1, ii, iq(jxyz)), &
                     1)

        end do
      end do
    end do

  end subroutine

  subroutine gga_distribute_qi(s, v, irep, fmat, ao, ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s(3)
    real(8), intent(in)    :: v(3, 3)
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(nr_ao_total)
    real(8), intent(in)    :: ao1(nr_ao_total, 3)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
    real(8)                :: f
!   ----------------------------------------------------------------------------

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do i = 1, nr2

          ii = st2 + i - 1

          f =       s(   jxyz)*ao (ii)    &
            + 2.0d0*v(1, jxyz)*ao1(ii, 1) &
            + 2.0d0*v(2, jxyz)*ao1(ii, 2) &
            + 2.0d0*v(3, jxyz)*ao1(ii, 3)

          call daxpy(nr1,                     &
                     f,                       &
                     ao(st1),                 &
                     1,                       &
                     fmat(st1, ii, iq(jxyz)), &
                     1)

        end do
      end do
    end do

  end subroutine

end module
