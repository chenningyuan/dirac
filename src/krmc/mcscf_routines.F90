module mcscf_routines

public :: raddfq90, m2dnz390

private :: get_fq

contains

  subroutine raddfq90(fq,h2xy,h2yx,pv,nxy,ix,iy,irepij,iprint, &
                      nz,ipqtoq,lupri,norbt,nasht,nnashx,key)
  !***********************************************************************
  !
  !     Add up the infinite number of terms that constitutes the FQ fock
  !     matrix.
  !
  !     For proper documentation see Ph.D. thesis of Joern Thyssen where
  !     all the formulae below should be documented :-)
  !
  !     Formulas (B.40) - (B.43)
  !
  !     Input:
  !        IREPIJ: symmetry of integrals (currently unused)
  !        H2XY, H2YX: (p,v) and (v,p)^T integrals for compound index XY
  !        NXY   : the XY compound index
  !        IX,IY : x and y
  !        PV    : 2-particle density matrix in (NZ,3) format, see appendix B3
  !
  !     Output:
  !        FQ    : the aux. Fock matrix
  !
  !     Written by J. Thyssen - Nov 16 2000
  !
  !***********************************************************************
  use unit_test_generator
  implicit none

#include "unit-test.h"

  ! External variables
  integer, intent(in), optional :: key
  integer, intent(in) :: nxy, ix, iy, irepij, iprint
  integer, intent(in) :: nz, ipqtoq(4,0:7), lupri, norbt, nasht, nnashx
  real*8,  intent(in) :: h2xy(norbt,nasht,nz,3), h2yx(norbt,nasht,nz,3)
  real*8,  intent(in) :: pv(nasht,nasht,nnashx,nz,3)
  real*8,  intent(inout) :: fq(norbt,nasht,nz)

  ! Internal variables
  logical :: create_outputs
  integer :: i, j, k, izz, npvcol, ier

  real*8, parameter   :: d1 = 1.0d0, dm1 = -1.0d0

  !------------------------------------------------------------------------------!
  !           Write to the file input information for unit testing               !
  !------------------------------------------------------------------------------!
  create_outputs = .false.

#ifdef CREATE_UNIT_TESTS
  call set_suffix(nxy)
  call set_suffix(ix)
  call set_suffix(iy)
  call add_input_file ('raddfq90.inp')
  call set_driver_routine('raddfq90')
  call add_module('mcscf_routines')

  call add_opt_variable('key   ','1','integer')

  call add_inp_variable(nxy   ,'nxy   ')
  call add_inp_variable(ix    ,'ix    ')
  call add_inp_variable(iy    ,'iy    ')
  call add_inp_variable(irepij,'irepij')
  call add_inp_variable(iprint,'iprint')
  call add_inp_variable(nz    ,'nz    ')
  call add_inp_variable(lupri ,'lupri ')
  call add_inp_variable(norbt ,'norbt ')
  call add_inp_variable(nasht ,'nasht ')
  call add_inp_variable(nnashx,'nnashx')
  call add_inp_variable(ipqtoq,'ipqtoq')
  call add_inp_variable(h2xy  ,'h2xy  ')
  call add_inp_variable(h2yx  ,'h2yx  ')
  call add_inp_variable(pv    ,'pv    ')

  call add_inp_variable(fq    ,'fq    ')
  call add_out_variable(fq    ,'fq    ',compare='file-file',file='raddfq90.out')

  call generate_test()

  create_outputs = .true.
#endif

  !------------------------------------------------------------------------------!
  !                               Print section                                  !
  !------------------------------------------------------------------------------!
  if(iprint .ge. 20)then
    call header('Output from raddfq90',-1)
    write(lupri,'(a,i3)') 'nxy       = ',nxy
    write(lupri,'(a,i3)') 'ix        = ',ix
    write(lupri,'(a,i3)') 'iy        = ',iy
    write(lupri,'(a,i3)') 'irepij    = ',irepij
    write(lupri,'(/a)') ' fq on entry'
    call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
    npvcol = nasht * nnashx
    do i = 1, 3
      write(lupri,'(a,i2,a)') ' pv density matrix',i,' of 3'
      call prqmat(pv(1,1,nxy,1,i),nasht,nasht,nasht,npvcol,nz,ipqtoq(1,0),lupri)
    end do
    do i = 1, 3
      write(lupri,'(a,i2,a)') ' h2xy matrix',i,' of 3'
      call prqmat(h2xy(1,1,1,i),norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
    end do
    if (ix /= iy) then
      do i = 1, 3
        write(lupri,'(a,i2,a)') ' h2yx matrix',i,' of 3'
        call prqmat(h2yx(1,1,1,i),norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
      end do
    end if
  end if

  !------------------------------------------------------------------------------!
  !                          *** Real part of FQ ***                             !
  !------------------------------------------------------------------------------!
  ! (B.40): each "term" is a line of the formula

  !--------------------!
  ! Term 1 (real part) !
  !--------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;1,1}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,1,1,1,1)

  if(ix /= iy)then
    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;1,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,d1,1,1,1,1,1)
  endif

  !--------------------!
  ! Term 2 (real part) !
  !--------------------!

  if(nz >= 2)then

    ! - \sum_{x >= y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;2,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,1,2,1,1)

    if(ix /= iy)then
      ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;2,1}
      call get_fq('n','n',h2yx,pv,fq,nxy,dm1,2,1,2,1,1)
    endif

  endif

  !--------------------------!
  ! Term 3 and 4 (real part) !
  !--------------------------!

  if(nz == 4 .and. ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;3,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,1,3,1,1)

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;4,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,4,1,4,1,1)

  endif

  !--------------------!
  ! Term 5 (real part) !
  !--------------------!

  if(ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,2} P_{qv,xy;1,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,2,1,2,1)

  endif

  !--------------------!
  ! Term 6 (real part) !
  !--------------------!

  if(nz >= 2 .and. ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,2} P_{qv,xy;2,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,2,2,2,1)

  endif

  if(nz == 4)then

  !--------------------!
  ! Term 7 (real part) !
  !--------------------!

    ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{3,2} P_{qv,xy;3,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,2,3,2,1)

    if(ix /= iy)then

      ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,3} P_{vq,xy;3,3}
      call get_fq('n','n',h2yx,pv,fq,nxy,d1,3,3,3,3,1)

    endif

  !--------------------!
  ! Term 8 (real part) !
  !--------------------!

    ! - \sum_{x >= y} \sum_v \Twoint{pv|xy}_{4,2} P_{qv,xy;4,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,4,2,4,2,1)

    if(ix /= iy)then

      ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,3} P_{vq,xy;4,3}
      call get_fq('n','n',h2yx,pv,fq,nxy,dm1,4,3,4,3,1)

    endif

  endif

  !--------------------!
  ! Term 9 (real part) !
  !--------------------!

  if(ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;1,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,3,1,3,1)

  endif

  !---------------------!
  ! Term 10 (real part) !
  !---------------------!

  if(nz >= 2 .and. ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;2,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,3,2,3,1)

  endif

  if(nz == 4 .and. ix /= iy)then

  !---------------------!
  ! Term 11 (real part) !
  !---------------------!

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,1} P_{vq,xy;3,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,d1,3,1,3,1,1)

  !---------------------!
  ! Term 12 (real part) !
  !---------------------!

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,1} P_{vq,xy;4,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,4,1,4,1,1)

  endif

  ! Skip remaining 36 terms if NZ = 1.
  if(nz > 1)then ! goto 100

  !------------------------------------------------------------------------------!
  !                         *** i-imag part of FQ ***                            !
  !------------------------------------------------------------------------------!
  ! (B.41)

  !----------------------!
  ! Term 1 (i-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;2,1}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,1,2,1,2)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;2,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,1,1,2,1,2)

  endif

  !----------------------!
  ! Term 2 (i-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;1,1}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,1,1,1,2)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;1,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,2,1,1,1,2)

  endif

  !----------------------------!
  ! Term 3 and 4 (i-imag part) !
  !----------------------------!

  if(nz == 4 .and. ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;4,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,1,4,1,2)

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;3,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,4,1,3,1,2)

  endif

  if(ix /= iy)then

  !----------------------!
  ! Term 5 (i-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,2} P_{qv,xy;2,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,2,2,2,2)

  !----------------------!
  ! Term 6 (i-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,2} P_{qv,xy;1,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,2,1,2,2)

  endif

  if(nz == 4)then

  !----------------------!
  ! Term 7 (i-imag part) !
  !----------------------!

    ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{3,2} P_{qv,xy;4,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,2,4,2,2)

    if(ix /= iy)then

      ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,3} P_{vq,xy;4,3}
      call get_fq('n','n',h2yx,pv,fq,nxy,dm1,3,3,4,3,2)

    endif

  !----------------------!
  ! Term 8 (i-imag part) !
  !----------------------!

    ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{4,2} P_{qv,xy;3,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,4,2,3,2,2)

    if(ix /= iy)then

      ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,3} P_{vq,xy;3,3}
      call get_fq('n','n',h2yx,pv,fq,nxy,dm1,4,3,3,3,2)

    endif

  endif

  if(ix /= iy)then

  !----------------------!
  ! Term 9 (i-imag part) !
  !----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;2,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,1,3,2,3,2)

  !-----------------------!
  ! Term 10 (i-imag part) !
  !-----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;1,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,3,1,3,2)

  endif

  if(nz == 4 .and. ix /= iy)then

  !-----------------------!
  ! Term 11 (i-imag part) !
  !-----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,1} P_{vq,xy;4,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,3,1,4,1,2)

  !-----------------------!
  ! Term 12 (i-imag part) !
  !-----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,1} P_{vq,xy;3,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,4,1,3,1,2)

  endif

  ! Skip remaining 24 terms if NZ = 2.
  if(nz > 2)then ! goto 101

  !------------------------------------------------------------------------------!
  !                         *** j-imag part of FQ ***                            !
  !------------------------------------------------------------------------------!
  ! (B.42)

  !----------------------!
  ! Term 1 (j-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;3,3}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,1,3,3,3)

  if(ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;3,2}
    call get_fq('n','n',h2yx,pv,fq,nxy,d1,1,1,3,2,3)

  endif

  !----------------------!
  ! Term 2 (j-imag part) !
  !----------------------!

  ! - \sum_{x >= y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;4,3}
  call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,1,4,3,3)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;4,2}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,2,1,4,2,3)

  endif

  !----------------------!
  ! Term 3 (j-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{3,2} P_{vq,xy;1,1}
  call get_fq('n','n',h2xy,pv,fq,nxy,d1,3,2,1,1,3)

  if(ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,3} P_{qv,xy;1,1}
    call get_fq('n','t',h2yx,pv,fq,nxy,d1,3,3,1,1,3)

  endif

  !----------------------!
  ! Term 4 (j-imag part) !
  !----------------------!

  ! - \sum_{x >= y} \sum_v \Twoint{pv|xy}_{4,2} P_{vq,xy;2,1}
  call get_fq('n','n',h2xy,pv,fq,nxy,dm1,4,2,2,1,3)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,3} P_{qv,xy;2,1}
    call get_fq('n','t',h2yx,pv,fq,nxy,dm1,4,3,2,1,3)

  endif

  if(ix /= iy)then

  !----------------------!
  ! Term 5 (j-imag part) !
  !----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;3,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,1,3,3,1,3)

  !----------------------!
  ! Term 6 (j-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;4,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,3,4,1,3)

  !----------------------!
  ! Term 7 (j-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,2} P_{vq,xy;3,1}
    call get_fq('n','n',h2xy,pv,fq,nxy,d1,1,2,3,1,3)

  !----------------------!
  ! Term 8 (j-imag part) !
  !----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,2} P_{vq,xy;4,1}
    call get_fq('n','n',h2xy,pv,fq,nxy,dm1,2,2,4,1,3)

  !----------------------!
  ! Term 9 (j-imag part) !
  !----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,1} P_{qv,xy;1,2}
    call get_fq('n','t',h2yx,pv,fq,nxy,dm1,3,1,1,2,3)

  !-----------------------!
  ! Term 10 (j-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,1} P_{qv,xy;2,2}
    call get_fq('n','t',h2yx,pv,fq,nxy,d1,4,1,2,2,3)

  !-----------------------!
  ! Term 11 (j-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;1,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,1,1,3,3)

  !-----------------------!
  ! Term 12 (j-imag part) !
  !-----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;2,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,4,1,2,3,3)

  endif

  !------------------------------------------------------------------------------!
  !                         *** k-imag part of FQ ***                            !
  !------------------------------------------------------------------------------!
  ! (B.43)

  !----------------------!
  ! Term 1 (k-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;4,3}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,1,4,3,4)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;4,2}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,1,1,4,2,4)

  endif

  !----------------------!
  ! Term 2 (k-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;3,3}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,1,3,3,4)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;3,2}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,2,1,3,2,4)

  endif

  !----------------------!
  ! Term 3 (k-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{3,2} P_{vq,xy;2,1}
  call get_fq('n','n',h2xy,pv,fq,nxy,d1,3,2,2,1,4)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,3} P_{qv,xy;2,1}
    call get_fq('n','t',h2yx,pv,fq,nxy,dm1,3,3,2,1,4)

  endif

  !----------------------!
  ! Term 4 (k-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{4,2} P_{vq,xy;1,1}
  call get_fq('n','n',h2xy,pv,fq,nxy,d1,4,2,1,1,4)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,3} P_{qv,xy;1,1}
    call get_fq('n','t',h2yx,pv,fq,nxy,dm1,4,3,1,1,4)

  endif

  if(ix /= iy)then

  !----------------------!
  ! Term 5 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;4,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,3,4,1,4)

  !----------------------!
  ! Term 6 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;3,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,3,3,1,4)

  !----------------------!
  ! Term 7 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,2} P_{vq,xy;4,1}
    call get_fq('n','n',h2xy,pv,fq,nxy,d1,1,2,4,1,4)

  !----------------------!
  ! Term 8 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,2} P_{vq,xy;3,1}
    call get_fq('n','n',h2xy,pv,fq,nxy,d1,2,2,3,1,4)

  !----------------------!
  ! Term 9 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,1} P_{qv,xy;2,2}
    call get_fq('n','t',h2yx,pv,fq,nxy,d1,3,1,2,2,4)

  !-----------------------!
  ! Term 10 (k-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,1} P_{qv,xy;1,2}
    call get_fq('n','t',h2yx,pv,fq,nxy,d1,4,1,1,2,4)

  !-----------------------!
  ! Term 11 (k-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;2,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,1,2,3,4)

  !-----------------------!
  ! Term 12 (k-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;1,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,4,1,1,3,4)

  endif

  endif ! 101
  endif ! 100

  !----------------!
  ! Output section !
  !----------------!

  if (iprint .ge. 20) then
    write(lupri,'(a)') 'fq on exit'
    call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  end if

  !------------------------------------------------------------------------------!
  !           Write to the file ouput information for unit testing               !
  !------------------------------------------------------------------------------!
  if( present(key) .or. create_outputs )then

    call set_suffix(nxy)
    call set_suffix(ix)
    call set_suffix(iy)
    call add_output_file('raddfq90.out')

    ! open file
    open(60,file=get_output_name(1),status='unknown',form='formatted', &
                                    access='sequential',iostat=ier)
    if(ier /= 0)then
      write(6,*) 'raddfq90: ERROR while creating raddfq90.out file'
      stop
    endif

    rewind 60

    ! write output data to the file
    do i=1,nz
    do j=1,nasht
    do k=1,norbt
    write(60,*) fq(k,j,i)
    enddo
    enddo
    enddo

    ! close file
    close(60,status='keep')

    call unit_test_generator_destructor

  endif

  end subroutine raddfq90

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_fq(ch1,ch2,h2,pv,fq,nxy,alpha,izh,ich,izp,icp,izf)
  implicit none

  ! External variables
  character, intent( in) :: ch1*1, ch2*1
  integer,   intent( in) :: nxy, izh, ich, izp, icp, izf
  real*8,    intent( in) :: alpha
  real*8,    intent( in) :: h2(:,:,:,:), pv(:,:,:,:,:)
  real*8,    intent(out) :: fq(:,:,:)

  ! Internal variables
  integer :: norbt, nasht

  norbt = size(h2,1)
  nasht = size(h2,2)

  call dgemm(ch1,ch2,norbt,nasht,nasht,alpha, &
       h2(1,1,izh,ich),norbt,                 &
       pv(1,1,nxy,izp,icp),nasht,             &
       1.0d0,fq(1,1,izf),norbt)

  end subroutine get_fq

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE M2DNZ390(TSYMM,KSYMM,PSYMM,VMUUUU,VDUUUU,IPRINT, &
                      IPQTOQ,NASH,NASHT,NNASHX,NZ,NZINCI,NM5, &
                      NBSYM,NFSYM,LUPRI,SPINFR,LEVYLE)
!***********************************************************************
!
!     Transform a UUUU matrix in Molfdir format to a UUUU matrix
!     in Dirac (NZ,3) format.
!
!     Input:
!        TSYMM: transition density matrix symmetrization:
!           do <0|...|B> + <B|...0> symmetrization
!        KSYMM: do Kramers symmetrization (i.e., calculate P++)
!        PSYMM: particle symmetrization: pqrs + rspq
!
!     Written by J. Thyssen - Jul 21 2000
!     Last revision :
!
!***********************************************************************
      USE UNIT_TEST_GENERATOR
      IMPLICIT NONE

#include "unit-test.h"

      LOGICAL, INTENT( IN) :: PSYMM, KSYMM, TSYMM, SPINFR, LEVYLE
      INTEGER, INTENT( IN) :: NASHT, NNASHX, NZ, NZINCI, NM5, NBSYM, NFSYM, IPRINT, LUPRI
      INTEGER, INTENT( IN) :: NASH(2), IPQTOQ(4,0:7)

      REAL*8,  INTENT( IN) :: VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,NM5)
      REAL*8,  INTENT(OUT) :: VDUUUU(NASHT,NASHT,NNASHX,NZ,3)

      INTEGER :: I, J, K, L, M, IJ, II, JJ, KK, LL, ITR, JTR, KTR, LTR
      INTEGER :: IER, NCLASS, IFSYM, IZ, KR(-300:300)
      REAL*8  :: THR_PRINT
      INTEGER, SAVE :: ICOUNT = 0

  !------------------------------------------------------------------------------!
  !           Write to the file input information for unit testing               !
  !------------------------------------------------------------------------------!
      icount = icount + 1
      write(6,*) 'M2DNZ390',psymm,ksymm,tsymm,icount

#ifdef CREATE_UNIT_TESTS
      call set_suffix(psymm)
      call set_suffix(ksymm)
      call set_suffix(tsymm)
      call set_suffix(icount)
      call set_driver_routine('m2dnz390')
      call add_module('mcscf_routines')

      call add_inp_variable(psymm ,'psymm ')
      call add_inp_variable(ksymm ,'ksymm ')
      call add_inp_variable(tsymm ,'tsymm ')
      call add_inp_variable(spinfr,'spinfr')
      call add_inp_variable(levyle,'levyle')
      call add_inp_variable(nasht ,'nasht ')
      call add_inp_variable(nnashx,'nnashx')
      call add_inp_variable(nz    ,'nz    ')
      call add_inp_variable(nzinci,'nzinci')
      call add_inp_variable(nm5   ,'nm5   ')
      call add_inp_variable(nbsym ,'nbsym ')
      call add_inp_variable(nfsym ,'nfsym ')
      call add_inp_variable(iprint,'iprint')
      call add_inp_variable(lupri ,'lupri ')
      call add_inp_variable(nash  ,'nash  ')
      call add_inp_variable(ipqtoq,'ipqtoq')
      call add_inp_variable(vmuuuu,'vmuuuu')
#endif

  !------------------------------------------------------------------------------!
      NCLASS = NZ * NZ * NBSYM / NFSYM
      CALL DZERO(VDUUUU,NASHT * NASHT * NNASHX * NZ * 3)

      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO

 9001 FORMAT(4I5,5X,2F20.10)

      IF (IPRINT .GE. 30) THEN
         THR_PRINT = 1.0D-12
         WRITE(LUPRI,'(/A,3(1X,L1)//A)') &
              ' (M2DNZ3) TSYMM, KSYMM, PSYMM = ',TSYMM,KSYMM,PSYMM, &
              ' (M2DNZ3) Input VUUUU matrix in Molfdir format'
         DO L = 1,2*NASHT
            DO K = 1,2*NASHT
               DO J = 1,2*NASHT
                  DO I = 1,2*NASHT
                     IF ( (ABS(VMUUUU(I,J,K,L,1)) .GT. THR_PRINT   &
                      .OR. ABS(VMUUUU(I,J,K,L,2)) .GT. THR_PRINT)) &
                          WRITE(LUPRI,9001) &
                          I,J,K,L,(VMUUUU(I,J,K,L,IZ),IZ=1,NZINCI)
                  END DO
               END DO
            END DO
         END DO
      END IF


      IJ = 0
      DO II = 1, NASHT
         DO JJ = 1, II
            IJ = IJ + 1
            DO LL = 1, NASHT
               DO KK = 1, NASHT
                  I = KR(II)
                  J = KR(JJ)
                  K = KR(KK)
                  L = KR(LL)
                  ITR = KR(-II)
                  JTR = KR(-JJ)
                  KTR = KR(-KK)
                  LTR = KR(-LL)

!                 Calculate integral classes.

                  IF (TSYMM) THEN
                     IF (PSYMM) THEN
                        IF (KSYMM) THEN

!                          (1) Symm. transition density matrix
!                          (2) Do particle symmetrization
!                          (3) Do Kramer's symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) + VMUUUU(L,K,J,I,1) &
                                + VMUUUU(I,J,K,L,1) + VMUUUU(J,I,L,K,1) &
                                + VMUUUU(LTR,KTR,I,J,1) &
                                + VMUUUU(KTR,LTR,J,I,1) &
                                + VMUUUU(I,J,LTR,KTR,1) &
                                + VMUUUU(J,I,KTR,LTR,1) &
                                + VMUUUU(K,L,JTR,ITR,1) &
                                + VMUUUU(L,K,ITR,JTR,1) &
                                + VMUUUU(JTR,ITR,K,L,1) &
                                + VMUUUU(ITR,JTR,L,K,1) &
                                + VMUUUU(LTR,KTR,JTR,ITR,1) &
                                + VMUUUU(KTR,LTR,ITR,JTR,1) &
                                + VMUUUU(JTR,ITR,LTR,KTR,1) &
                                + VMUUUU(ITR,JTR,KTR,LTR,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(I,JTR,K,LTR,1) &
                                + VMUUUU(LTR,K,JTR,I,1) &
                                + VMUUUU(JTR,I,LTR,K,1) &
                                - VMUUUU(L,KTR,I,JTR,1) &
                                - VMUUUU(I,JTR,L,KTR,1) &
                                - VMUUUU(KTR,L,JTR,I,1) &
                                - VMUUUU(JTR,I,KTR,L,1) &
                                - VMUUUU(K,LTR,J,ITR,1) &
                                - VMUUUU(J,ITR,K,LTR,1) &
                                - VMUUUU(LTR,K,ITR,J,1) &
                                - VMUUUU(ITR,J,LTR,K,1) &
                                + VMUUUU(L,KTR,J,ITR,1) &
                                + VMUUUU(J,ITR,L,KTR,1) &
                                + VMUUUU(KTR,L,ITR,J,1) &
                                + VMUUUU(ITR,J,KTR,L,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(I,JTR,KTR,L,1) &
                                + VMUUUU(L,KTR,JTR,I,1) &
                                + VMUUUU(JTR,I,L,KTR,1) &
                                - VMUUUU(LTR,K,I,JTR,1) &
                                - VMUUUU(I,JTR,LTR,K,1) &
                                - VMUUUU(K,LTR,JTR,I,1) &
                                - VMUUUU(JTR,I,K,LTR,1) &
                                - VMUUUU(KTR,L,J,ITR,1) &
                                - VMUUUU(J,ITR,KTR,L,1) &
                                - VMUUUU(L,KTR,ITR,J,1) &
                                - VMUUUU(ITR,J,L,KTR,1) &
                                + VMUUUU(LTR,K,J,ITR,1) &
                                + VMUUUU(J,ITR,LTR,K,1) &
                                + VMUUUU(K,LTR,ITR,J,1) &
                                + VMUUUU(ITR,J,K,LTR,1)
                        ELSE

!                          (1) Symm. transition density matrix
!                          (2) Do particle symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) + VMUUUU(L,K,J,I,1) &
                                + VMUUUU(I,J,K,L,1) + VMUUUU(J,I,L,K,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(LTR,K,JTR,I,1) &
                                + VMUUUU(I,JTR,K,LTR,1) &
                                + VMUUUU(JTR,I,LTR,K,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(L,KTR,JTR,I,1) &
                                + VMUUUU(I,JTR,KTR,L,1) &
                                + VMUUUU(JTR,I,L,KTR,1)
                        END IF
                     ELSE
                        IF (KSYMM) THEN

!                          (1) Symm. transition density matrix
!                          (2) Kramers symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) &
                                + VMUUUU(L,K,J,I,1) &
                                + VMUUUU(LTR,KTR,I,J,1) &
                                + VMUUUU(KTR,LTR,J,I,1) &
                                + VMUUUU(K,L,JTR,ITR,1) &
                                + VMUUUU(L,K,ITR,JTR,1) &
                                + VMUUUU(LTR,KTR,JTR,ITR,1) &
                                + VMUUUU(KTR,LTR,ITR,JTR,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(LTR,K,JTR,I,1) &
                                - VMUUUU(L,KTR,I,JTR,1) &
                                - VMUUUU(KTR,L,JTR,I,1) &
                                - VMUUUU(K,LTR,J,ITR,1) &
                                - VMUUUU(LTR,K,ITR,J,1) &
                                + VMUUUU(L,KTR,J,ITR,1) &
                                + VMUUUU(KTR,L,ITR,J,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(L,KTR,JTR,I,1) &
                                - VMUUUU(LTR,K,I,JTR,1) &
                                - VMUUUU(K,LTR,JTR,I,1) &
                                - VMUUUU(KTR,L,J,ITR,1) &
                                - VMUUUU(L,KTR,ITR,J,1) &
                                + VMUUUU(LTR,K,J,ITR,1) &
                                + VMUUUU(K,LTR,ITR,J,1)
                        ELSE

!                          (1) Symm. transition density matrix

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) + VMUUUU(L,K,J,I,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(LTR,K,JTR,I,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(L,KTR,JTR,I,1)
                        END IF
                     END IF
                  ELSE
                     IF (PSYMM) THEN
                        IF (KSYMM) THEN

!                          (1) Particle symmetrization
!                          (2) Kramers symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) &
                                + VMUUUU(I,J,K,L,1) &
                                + VMUUUU(LTR,KTR,I,J,1) &
                                + VMUUUU(I,J,LTR,KTR,1) &
                                + VMUUUU(K,L,JTR,ITR,1) &
                                + VMUUUU(JTR,ITR,K,L,1) &
                                + VMUUUU(LTR,KTR,JTR,ITR,1) &
                                + VMUUUU(JTR,ITR,LTR,KTR,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(I,JTR,K,LTR,1) &
                                - VMUUUU(L,KTR,I,JTR,1) &
                                - VMUUUU(I,JTR,L,KTR,1) &
                                - VMUUUU(K,LTR,J,ITR,1) &
                                - VMUUUU(J,ITR,K,LTR,1) &
                                + VMUUUU(L,KTR,J,ITR,1) &
                                + VMUUUU(J,ITR,L,KTR,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(I,JTR,KTR,L,1) &
                                - VMUUUU(LTR,K,I,JTR,1) &
                                - VMUUUU(I,JTR,LTR,K,1) &
                                - VMUUUU(KTR,L,J,ITR,1) &
                                - VMUUUU(J,ITR,KTR,L,1) &
                                + VMUUUU(LTR,K,J,ITR,1) &
                                + VMUUUU(J,ITR,LTR,K,1)
                        ELSE

!                          (1) Particle symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) &
                                + VMUUUU(I,J,K,L,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(I,JTR,K,LTR,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(I,JTR,KTR,L,1)
                        END IF
                     ELSE
                        IF (KSYMM) THEN

!                          (1) Kramers symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) &
                                + VMUUUU(LTR,KTR,I,J,1) &
                                + VMUUUU(K,L,JTR,ITR,1) &
                                + VMUUUU(LTR,KTR,JTR,ITR,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                - VMUUUU(L,KTR,I,JTR,1) &
                                - VMUUUU(K,LTR,J,ITR,1) &
                                + VMUUUU(L,KTR,J,ITR,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                - VMUUUU(LTR,K,I,JTR,1) &
                                - VMUUUU(KTR,L,J,ITR,1) &
                                + VMUUUU(LTR,K,J,ITR,1)
                        ELSE

!                          No symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = VMUUUU(K,L,I,J,1)
                           VDUUUU(KK,LL,IJ,1,2) = VMUUUU(K,LTR,I,JTR,1)
                           VDUUUU(KK,LL,IJ,1,3) = VMUUUU(KTR,L,I,JTR,1)
                        END IF
                     END IF
                  END IF

!                 Complex point groups
!                 --------------------

                  IF (NZ .GE. 2) THEN
                     IF (TSYMM) THEN
                        IF (PSYMM) THEN
                           IF (KSYMM) THEN

!                             (1) Symm. transition density matrix
!                             (2) Do particle symmetrization
!                             (3) Do Kramer's symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   + VMUUUU(I,J,K,L,2) &
                                   - VMUUUU(L,K,J,I,2) &
                                   - VMUUUU(J,I,L,K,2) &
                                   + VMUUUU(LTR,KTR,I,J,2) &
                                   + VMUUUU(I,J,LTR,KTR,2) &
                                   - VMUUUU(KTR,LTR,J,I,2) &
                                   - VMUUUU(J,I,KTR,LTR,2) &
                                   + VMUUUU(K,L,JTR,ITR,2) &
                                   + VMUUUU(JTR,ITR,K,L,2) &
                                   - VMUUUU(L,K,ITR,JTR,2) &
                                   - VMUUUU(ITR,JTR,L,K,2) &
                                   + VMUUUU(LTR,KTR,JTR,ITR,2) &
                                   + VMUUUU(JTR,ITR,LTR,KTR,2) &
                                   - VMUUUU(KTR,LTR,ITR,JTR,2) &
                                   - VMUUUU(ITR,JTR,KTR,LTR,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   + VMUUUU(I,JTR,K,LTR,2) &
                                   - VMUUUU(LTR,K,JTR,I,2) &
                                   - VMUUUU(JTR,I,LTR,K,2) &
                                   - VMUUUU(L,KTR,I,JTR,2) &
                                   - VMUUUU(I,JTR,L,KTR,2) &
                                   + VMUUUU(KTR,L,JTR,I,2) &
                                   + VMUUUU(JTR,I,KTR,L,2) &
                                   - VMUUUU(K,LTR,J,ITR,2) &
                                   - VMUUUU(J,ITR,K,LTR,2) &
                                   + VMUUUU(LTR,K,ITR,J,2) &
                                   + VMUUUU(ITR,J,LTR,K,2) &
                                   + VMUUUU(L,KTR,J,ITR,2) &
                                   + VMUUUU(J,ITR,L,KTR,2) &
                                   - VMUUUU(KTR,L,ITR,J,2) &
                                   - VMUUUU(ITR,J,KTR,L,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   + VMUUUU(I,JTR,KTR,L,2) &
                                   - VMUUUU(L,KTR,JTR,I,2) &
                                   - VMUUUU(JTR,I,L,KTR,2) &
                                   - VMUUUU(LTR,K,I,JTR,2) &
                                   - VMUUUU(I,JTR,LTR,K,2) &
                                   + VMUUUU(K,LTR,JTR,I,2) &
                                   + VMUUUU(JTR,I,K,LTR,2) &
                                   - VMUUUU(KTR,L,J,ITR,2) &
                                   - VMUUUU(J,ITR,KTR,L,2) &
                                   + VMUUUU(L,KTR,ITR,J,2) &
                                   + VMUUUU(ITR,J,L,KTR,2) &
                                   + VMUUUU(LTR,K,J,ITR,2) &
                                   + VMUUUU(J,ITR,LTR,K,2) &
                                   - VMUUUU(K,LTR,ITR,J,2) &
                                   - VMUUUU(ITR,J,K,LTR,2)
                           ELSE

!                            (1) Symm. transition density matrix
!                            (2) Do particle symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   - VMUUUU(L,K,J,I,2) &
                                   + VMUUUU(I,J,K,L,2) &
                                   - VMUUUU(J,I,L,K,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   - VMUUUU(LTR,K,JTR,I,2) &
                                   + VMUUUU(I,JTR,K,LTR,2) &
                                   - VMUUUU(JTR,I,LTR,K,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   - VMUUUU(L,KTR,JTR,I,2) &
                                   + VMUUUU(I,JTR,KTR,L,2) &
                                   - VMUUUU(JTR,I,L,KTR,2)
                           END IF
                        ELSE
                           IF (KSYMM) THEN

!                            (1) Symm. transition density matrix
!                            (2) Kramers symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   - VMUUUU(L,K,J,I,2) &
                                   + VMUUUU(LTR,KTR,I,J,2) &
                                   - VMUUUU(KTR,LTR,J,I,2) &
                                   + VMUUUU(K,L,JTR,ITR,2) &
                                   - VMUUUU(L,K,ITR,JTR,2) &
                                   + VMUUUU(LTR,KTR,JTR,ITR,2) &
                                   - VMUUUU(KTR,LTR,ITR,JTR,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   - VMUUUU(LTR,K,JTR,I,2) &
                                   - VMUUUU(L,KTR,I,JTR,2) &
                                   + VMUUUU(KTR,L,JTR,I,2) &
                                   - VMUUUU(K,LTR,J,ITR,2) &
                                   + VMUUUU(LTR,K,ITR,J,2) &
                                   + VMUUUU(L,KTR,J,ITR,2) &
                                   - VMUUUU(KTR,L,ITR,J,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   - VMUUUU(L,KTR,JTR,I,2) &
                                   - VMUUUU(LTR,K,I,JTR,2) &
                                   + VMUUUU(K,LTR,JTR,I,2) &
                                   - VMUUUU(KTR,L,J,ITR,2) &
                                   + VMUUUU(L,KTR,ITR,J,2) &
                                   + VMUUUU(LTR,K,J,ITR,2) &
                                   - VMUUUU(K,LTR,ITR,J,2)
                           ELSE

!                            (1) Symm. transition density matrix

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   - VMUUUU(L,K,J,I,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   - VMUUUU(LTR,K,JTR,I,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   - VMUUUU(L,KTR,JTR,I,2)
                           END IF
                        END IF
                     ELSE
                        IF (PSYMM) THEN
                           IF (KSYMM) THEN

!                            (1) Particle symmetrization
!                            (2) Kramers symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   + VMUUUU(I,J,K,L,2) &
                                   + VMUUUU(LTR,KTR,I,J,2) &
                                   + VMUUUU(I,J,LTR,KTR,2) &
                                   + VMUUUU(K,L,JTR,ITR,2) &
                                   + VMUUUU(JTR,ITR,K,L,2) &
                                   + VMUUUU(LTR,KTR,JTR,ITR,2) &
                                   + VMUUUU(JTR,ITR,LTR,KTR,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   + VMUUUU(I,JTR,K,LTR,2) &
                                   - VMUUUU(L,KTR,I,JTR,2) &
                                   - VMUUUU(I,JTR,L,KTR,2) &
                                   - VMUUUU(K,LTR,J,ITR,2) &
                                   - VMUUUU(J,ITR,K,LTR,2) &
                                   + VMUUUU(L,KTR,J,ITR,2) &
                                   + VMUUUU(J,ITR,L,KTR,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   + VMUUUU(I,JTR,KTR,L,2) &
                                   - VMUUUU(LTR,K,I,JTR,2) &
                                   - VMUUUU(I,JTR,LTR,K,2) &
                                   - VMUUUU(KTR,L,J,ITR,2) &
                                   - VMUUUU(J,ITR,KTR,L,2) &
                                   + VMUUUU(LTR,K,J,ITR,2) &
                                   + VMUUUU(J,ITR,LTR,K,2)
                           ELSE

!                            (1) Particle symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   + VMUUUU(I,J,K,L,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   + VMUUUU(I,JTR,K,LTR,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   + VMUUUU(I,JTR,KTR,L,2)
                           END IF
                        ELSE
                           IF (KSYMM) THEN

!                            (1) Kramers symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   + VMUUUU(LTR,KTR,I,J,2) &
                                   + VMUUUU(K,L,JTR,ITR,2) &
                                   + VMUUUU(LTR,KTR,JTR,ITR,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   - VMUUUU(L,KTR,I,JTR,2) &
                                   - VMUUUU(K,LTR,J,ITR,2) &
                                   + VMUUUU(L,KTR,J,ITR,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   - VMUUUU(LTR,K,I,JTR,2) &
                                   - VMUUUU(KTR,L,J,ITR,2) &
                                   + VMUUUU(LTR,K,J,ITR,2)
                           ELSE

!                             No symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2)
                           END IF
                        END IF
                     END IF
                  END IF

          if(nz == 4)then  ! goto nz4

          if(tsymm)then
            if(psymm)then
              if(ksymm)then

!               (1) Symm. transition density matrix
!               (2) Particle symmetrization
!               (3) Kramer's symmetrization
!                   Order of the implementation: (3), (1), (2)

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) - vmuuuu(ltr,ktr,j  ,itr,1) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,1) - vmuuuu(k  ,l  ,j  ,itr,1) &
                                     + vmuuuu(l  ,k  ,jtr,i  ,1) - vmuuuu(ktr,ltr,itr,j  ,1) &
                                     + vmuuuu(ktr,ltr,jtr,i  ,1) - vmuuuu(l  ,k  ,itr,j  ,1) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,1) - vmuuuu(j  ,itr,ltr,ktr,1) &
                                     + vmuuuu(i  ,jtr,ltr,ktr,1) - vmuuuu(j  ,itr,k  ,l  ,1) &
                                     + vmuuuu(jtr,i  ,l  ,k  ,1) - vmuuuu(itr,j  ,ktr,ltr,1) &
                                     + vmuuuu(jtr,i  ,ktr,ltr,1) - vmuuuu(itr,j  ,l  ,k  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(ltr,ktr,j  ,itr,2) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,2) - vmuuuu(k  ,l  ,j  ,itr,2) &
                                     - vmuuuu(l  ,k  ,jtr,i  ,2) + vmuuuu(ktr,ltr,itr,j  ,2) &
                                     - vmuuuu(ktr,ltr,jtr,i  ,2) + vmuuuu(l  ,k  ,itr,j  ,2) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,2) - vmuuuu(j  ,itr,ltr,ktr,2) &
                                     + vmuuuu(i  ,jtr,ltr,ktr,2) - vmuuuu(j  ,itr,k  ,l  ,2) &
                                     - vmuuuu(jtr,i  ,l  ,k  ,2) + vmuuuu(itr,j  ,ktr,ltr,2) &
                                     - vmuuuu(jtr,i  ,ktr,ltr,2) + vmuuuu(itr,j  ,l  ,k  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) - vmuuuu(l  ,ktr,jtr,itr,1) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,1) + vmuuuu(k  ,ltr,jtr,itr,1) &
                                     + vmuuuu(ltr,k  ,j  ,i  ,1) - vmuuuu(ktr,l  ,itr,jtr,1) &
                                     - vmuuuu(ktr,l  ,j  ,i  ,1) + vmuuuu(ltr,k  ,itr,jtr,1) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,1) - vmuuuu(jtr,itr,l  ,ktr,1) &
                                     - vmuuuu(i  ,j  ,l  ,ktr,1) + vmuuuu(jtr,itr,k  ,ltr,1) &
                                     + vmuuuu(j  ,i  ,ltr,k  ,1) - vmuuuu(itr,jtr,ktr,l  ,1) &
                                     - vmuuuu(j  ,i  ,ktr,l  ,1) + vmuuuu(itr,jtr,ltr,k  ,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(l  ,ktr,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,2) + vmuuuu(k  ,ltr,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,j  ,i  ,2) + vmuuuu(ktr,l  ,itr,jtr,2) &
                                     + vmuuuu(ktr,l  ,j  ,i  ,2) - vmuuuu(ltr,k  ,itr,jtr,2) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,2) - vmuuuu(jtr,itr,l  ,ktr,2) &
                                     - vmuuuu(i  ,j  ,l  ,ktr,2) + vmuuuu(jtr,itr,k  ,ltr,2) &
                                     - vmuuuu(j  ,i  ,ltr,k  ,2) + vmuuuu(itr,jtr,ktr,l  ,2) &
                                     + vmuuuu(j  ,i  ,ktr,l  ,2) - vmuuuu(itr,jtr,ltr,k  ,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) - vmuuuu(ltr,k  ,jtr,itr,1) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,1) + vmuuuu(ktr,l  ,jtr,itr,1) &
                                     + vmuuuu(l  ,ktr,j  ,i  ,1) - vmuuuu(k  ,ltr,itr,jtr,1) &
                                     - vmuuuu(k  ,ltr,j  ,i  ,1) + vmuuuu(l  ,ktr,itr,jtr,1) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,1) - vmuuuu(jtr,itr,ltr,k  ,1) &
                                     - vmuuuu(i  ,j  ,ltr,k  ,1) + vmuuuu(jtr,itr,ktr,l  ,1) &
                                     + vmuuuu(j  ,i  ,l  ,ktr,1) - vmuuuu(itr,jtr,k  ,ltr,1) &
                                     - vmuuuu(j  ,i  ,k  ,ltr,1) + vmuuuu(itr,jtr,l  ,ktr,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(ltr,k  ,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,2) + vmuuuu(ktr,l  ,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,j  ,i  ,2) + vmuuuu(k  ,ltr,itr,jtr,2) &
                                     + vmuuuu(k  ,ltr,j  ,i  ,2) - vmuuuu(l  ,ktr,itr,jtr,2) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,2) - vmuuuu(jtr,itr,ltr,k  ,2) &
                                     - vmuuuu(i  ,j  ,ltr,k  ,2) + vmuuuu(jtr,itr,ktr,l  ,2) &
                                     - vmuuuu(j  ,i  ,l  ,ktr,2) + vmuuuu(itr,jtr,k  ,ltr,2) &
                                     + vmuuuu(j  ,i  ,k  ,ltr,2) - vmuuuu(itr,jtr,l  ,ktr,2)
              else

!               (1) Symm. transition density matrix
!               (2) Particle symmetrization
!                   Order of the implementation: (1), (2)

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) + vmuuuu(l  ,k  ,jtr,i  ,1) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,1) + vmuuuu(jtr,i  ,l  ,k  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(l  ,k  ,jtr,i  ,2) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,2) - vmuuuu(jtr,i  ,l  ,k  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) + vmuuuu(ltr,k  ,j  ,i  ,1) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,1) + vmuuuu(j  ,i  ,ltr,k  ,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(ltr,k  ,j  ,i  ,2) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,2) - vmuuuu(j  ,i  ,ltr,k  ,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) + vmuuuu(l  ,ktr,j  ,i  ,1) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,1) + vmuuuu(j  ,i  ,l  ,ktr,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(l  ,ktr,j  ,i  ,2) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,2) - vmuuuu(j  ,i  ,l  ,ktr,2)
              endif
            else
              if(ksymm)then

!               (1) Symm. transition density matrix
!               (3) Kramers symmetrization
!                   Order of the implementation: (3), (1)

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) - vmuuuu(ltr,ktr,j  ,itr,1) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,1) - vmuuuu(k  ,l  ,j  ,itr,1) &
                                     + vmuuuu(l  ,k  ,jtr,i  ,1) - vmuuuu(ktr,ltr,itr,j  ,1) &
                                     + vmuuuu(ktr,ltr,jtr,i  ,1) - vmuuuu(l  ,k  ,itr,j  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(ltr,ktr,j  ,itr,2) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,2) - vmuuuu(k  ,l  ,j  ,itr,2) &
                                     - vmuuuu(l  ,k  ,jtr,i  ,2) + vmuuuu(ktr,ltr,itr,j  ,2) &
                                     - vmuuuu(ktr,ltr,jtr,i  ,2) + vmuuuu(l  ,k  ,itr,j  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) - vmuuuu(l  ,ktr,jtr,itr,1) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,1) + vmuuuu(k  ,ltr,jtr,itr,1) &
                                     + vmuuuu(ltr,k  ,j  ,i  ,1) - vmuuuu(ktr,l  ,itr,jtr,1) &
                                     - vmuuuu(ktr,l  ,j  ,i  ,1) + vmuuuu(ltr,k  ,itr,jtr,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(l  ,ktr,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,2) + vmuuuu(k  ,ltr,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,j  ,i  ,2) + vmuuuu(ktr,l  ,itr,jtr,2) &
                                     + vmuuuu(ktr,l  ,j  ,i  ,2) - vmuuuu(ltr,k  ,itr,jtr,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) - vmuuuu(ltr,k  ,jtr,itr,1) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,1) + vmuuuu(ktr,l  ,jtr,itr,1) &
                                     + vmuuuu(l  ,ktr,j  ,i  ,1) - vmuuuu(k  ,ltr,itr,jtr,1) &
                                     - vmuuuu(k  ,ltr,j  ,i  ,1) + vmuuuu(l  ,ktr,itr,jtr,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(ltr,k  ,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,2) + vmuuuu(ktr,l  ,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,j  ,i  ,2) + vmuuuu(k  ,ltr,itr,jtr,2) &
                                     + vmuuuu(k  ,ltr,j  ,i  ,2) - vmuuuu(l  ,ktr,itr,jtr,2)
              else

!               (1) Symm. transition density matrix

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) + vmuuuu(l  ,k  ,jtr,i  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(l  ,k  ,jtr,i  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) + vmuuuu(ltr,k  ,j  ,i  ,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(ltr,k  ,j  ,i  ,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) + vmuuuu(l  ,ktr,j  ,i  ,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(l  ,ktr,j  ,i  ,2)
              endif
            endif
          else
            if(psymm)then
              if(ksymm)then

!               (2) Particle symmetrization
!               (3) Kramers symmetrization
!                   Order of the implementation: (3), (2)

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) - vmuuuu(ltr,ktr,j  ,itr,1) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,1) - vmuuuu(k  ,l  ,j  ,itr,1) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,1) - vmuuuu(j  ,itr,ltr,ktr,1) &
                                     + vmuuuu(i  ,jtr,ltr,ktr,1) - vmuuuu(j  ,itr,k  ,l  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(ltr,ktr,j  ,itr,2) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,2) - vmuuuu(k  ,l  ,j  ,itr,2) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,2) - vmuuuu(j  ,itr,ltr,ktr,2) &
                                     + vmuuuu(i  ,jtr,ltr,ktr,2) - vmuuuu(j  ,itr,k  ,l  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) - vmuuuu(l  ,ktr,jtr,itr,1) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,1) + vmuuuu(k  ,ltr,jtr,itr,1) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,1) - vmuuuu(jtr,itr,l  ,ktr,1) &
                                     - vmuuuu(i  ,j  ,l  ,ktr,1) + vmuuuu(jtr,itr,k  ,ltr,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(l  ,ktr,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,2) + vmuuuu(k  ,ltr,jtr,itr,2) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,2) - vmuuuu(jtr,itr,l  ,ktr,2) &
                                     - vmuuuu(i  ,j  ,l  ,ktr,2) + vmuuuu(jtr,itr,k  ,ltr,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) - vmuuuu(ltr,k  ,jtr,itr,1) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,1) + vmuuuu(ktr,l  ,jtr,itr,1) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,1) - vmuuuu(jtr,itr,ltr,k  ,1) &
                                     - vmuuuu(i  ,j  ,ltr,k  ,1) + vmuuuu(jtr,itr,ktr,l  ,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(ltr,k  ,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,2) + vmuuuu(ktr,l  ,jtr,itr,2) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,2) - vmuuuu(jtr,itr,ltr,k  ,2) &
                                     - vmuuuu(i  ,j  ,ltr,k  ,2) + vmuuuu(jtr,itr,ktr,l  ,2)
              else

!               (2) Particle symmetrization

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) + vmuuuu(i  ,jtr,k  ,l  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) + vmuuuu(i  ,jtr,k  ,l  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) + vmuuuu(i  ,j  ,k  ,ltr,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) + vmuuuu(i  ,j  ,k  ,ltr,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) + vmuuuu(i  ,j  ,ktr,l  ,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) + vmuuuu(i  ,j  ,ktr,l  ,2)
              endif
            else
              if(ksymm)then

!               (3) Kramers symmetrization

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) - vmuuuu(ltr,ktr,j  ,itr,1) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,1) - vmuuuu(k  ,l  ,j  ,itr,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(ltr,ktr,j  ,itr,2) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,2) - vmuuuu(k  ,l  ,j  ,itr,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) - vmuuuu(l  ,ktr,jtr,itr,1) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,1) + vmuuuu(k  ,ltr,jtr,itr,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(l  ,ktr,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,2) + vmuuuu(k  ,ltr,jtr,itr,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) - vmuuuu(ltr,k  ,jtr,itr,1) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,1) + vmuuuu(ktr,l  ,jtr,itr,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(ltr,k  ,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,2) + vmuuuu(ktr,l  ,jtr,itr,2)
              else

!               No symmetrization

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2)
              endif
            endif
          endif

          endif  ! label nz4

! END OF LOOPS FOR UUUU matrix transformation
        enddo
      enddo
    enddo
  enddo

!     remove integrals that should be zero in spinfree calculations but
!     are not due to numerical noise
      if(spinfr.or.levyle)then
        do i = 2,3
          call dzero(VDUUUU(1,1,1,1,i),NASHT*NASHT*NNASHX*NZ)
        end do
      end if

      IF (IPRINT .GE. 30) THEN
         CALL HEADER('VUUUU matrix in Dirac (NZ,3) format',-1)
         CALL PRDNZ3(VDUUUU,NASHT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
      END IF

  !------------------------------------------------------------------------------!
  !           Write to the file ouput information for unit testing               !
  !------------------------------------------------------------------------------!
#ifdef CREATE_UNIT_TESTS
    call add_out_variable(vduuuu,'vduuuu')
    call generate_test()
#endif

  end subroutine m2dnz390

end module mcscf_routines
