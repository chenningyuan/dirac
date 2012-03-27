module mcscf_routines

public raddfq90

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
  real*8, allocatable :: tmpmat(:,:,:)

  !------------------------------------------------------------------------------!
  !           Write to the file input information for unit testing               !
  !------------------------------------------------------------------------------!
  create_outputs = .false.

#ifdef CREATE_UNIT_TESTS
  call add_input_file ('raddfq90.inp_' // i2c(nxy) //'_'// i2c(ix) //'_'// i2c(iy) )
  call add_output_file('raddfq90.out_' // i2c(nxy) //'_'// i2c(ix) //'_'// i2c(iy) )
  call set_driver_routine('raddfq90')
  call add_module('mcscf_routines')
  call set_test_name(nxy)
  call set_test_name(ix)
  call set_test_name(iy)

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

  call generate_test
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
    if (ix .ne. iy) then
      do i = 1, 3
        write(lupri,'(a,i2,a)') ' h2yx matrix',i,' of 3'
        call prqmat(h2yx(1,1,1,i),norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
      end do
    end if
  end if

!
!   ------------- START TEST -------------
!
  call header('Output from raddfq90',-1)
  write(lupri,'(a,3i3)') ' nxy, ix, iy  = ',nxy,ix,iy
  write(lupri,'(a,3i3)') ' norbt, nasht = ',norbt,nasht
  do izz=1,nz
    do i=1,3
      write(lupri,'(/a,i1,a,i1,$a)')' pv matrix(',izz,',',i,')'
      call pmattmp(lupri,pv(1,1,nxy,izz,i),nasht,nasht)
    enddo
  enddo
  allocate(tmpmat(norbt,nasht,3))
  ! h2xy
  if (nz==1)then
    do k=1,3
    do j=1,nasht
    do i=1,norbt
      tmpmat(i,j,k) = h2xy(i,j,1,k)
    enddo
    enddo
    enddo
  elseif(nz==2)then
    do k=1,3
      if(k/=3)then
        do j=1,nasht
        do i=1,norbt
          tmpmat(i,j,k) = h2xy(i,j,1,k) + h2xy(i,j,2,k)
        enddo
        enddo
      else
        do j=1,nasht
        do i=1,norbt
          tmpmat(i,j,k) = h2xy(i,j,1,k) - h2xy(i,j,2,k)
        enddo
        enddo
      endif
    enddo
  endif
  do i=1,3
    write(lupri,'(/a,i1,a,i1,$a)')' h2xy matrix(',nz,',',i,')'
    call pmattmp(lupri,tmpmat(1,1,i),norbt,nasht)
  enddo
  ! h2yx
  if (nz==1)then
    do k=1,3
      do j=1,nasht
      do i=1,norbt
        tmpmat(i,j,k) = h2yx(i,j,1,k)
      enddo
      enddo
      do izz=1,nz
        write(lupri,'(/a,i1,a,i1,$a)')' h2yx matrix(',izz,',',k,')'
        call pmattmp(lupri,tmpmat(1,1,k),norbt,nasht)
      enddo
    enddo
  elseif(nz==2)then
!    do j=1,nasht
!    do i=1,norbt
!      tmpmat(i,j,1) = h2yx(i,j,1,1) - h2yx(i,j,2,1)
!    enddo
!    enddo
!    i = 1
!    write(lupri,'(/a,i1,a,i1,$a)')' h2yx matrix(',nz,',',i,')'
!    call pmattmp(lupri,tmpmat(1,1,1),norbt,nasht)
    do i=1,3
      do izz=1,nz
        write(lupri,'(/a,i1,a,i1,$a)')' h2yx matrix(',izz,',',i,')'
        call pmattmp(lupri,h2yx(1,1,izz,i),norbt,nasht)
      enddo
    enddo
  endif
  deallocate(tmpmat)
!     do i=1,3
!       do izz=1,nz
!         write(lupri,'(/a,i1,a,i1,$a)')' h2xy matrix(',izz,',',i,')'
!         call pmattmp(lupri,h2xy(1,1,izz,i),norbt,nasht)
!       enddo
!     enddo
!     do izz=1,nz
!       write(lupri,'(/a,i1,$a)')' fq matrix(',izz,')'
!       call pmattmp(lupri,fq(1,1,izz),norbt,nasht)
!     enddo
!
!   ------------- END TEST -------------
!

! ***********************
! *** Real part of FQ ***
! ***********************
! (B.40): each "term" is a line of the formula

! ------------------
! Term 1 (real part)
! ------------------

! \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;1,1}

  call dgemm('n','t',norbt,nasht,nasht,d1, &
       h2xy(1,1,1,1),norbt,                &
       pv(1,1,nxy,1,1),nasht,              &
       d1,fq(1,1,1),norbt)

  if (ix .ne. iy) then

!    \sum_{x>y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;1,1}

     call dgemm('n','n',norbt,nasht,nasht,d1,&
          h2yx(1,1,1,1),norbt,&
          pv(1,1,nxy,1,1),nasht,&
          d1,fq(1,1,1),norbt)

  end if

! ------------------
! Term 2 (real part)
! ------------------

  if (nz .ge. 2) then

!    - \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;2,1}

     call dgemm('n','t',norbt,nasht,nasht,dm1,&
          h2xy(1,1,2,1),norbt,&
          pv(1,1,nxy,2,1),nasht,&
          d1,fq(1,1,1),norbt)

     if (ix .ne. iy) then

!       - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;2,1}

        call dgemm('n','n',norbt,nasht,nasht,dm1,&
             h2yx(1,1,2,1),norbt,&
             pv(1,1,nxy,2,1),nasht,&
             d1,fq(1,1,1),norbt)

     end if

  end if

! ------------------------
! Term 3 and 4 (real part)
! ------------------------

  if (nz .ge. 4 .and. ix .ne. iy) then

!    \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;3,1}

     call dgemm('n','t',norbt,nasht,nasht,d1,&
          h2xy(1,1,3,1),norbt,&
          pv(1,1,nxy,3,1),nasht,&
          d1,fq(1,1,1),norbt)

!    - \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;4,1}

     call dgemm('n','t',norbt,nasht,nasht,dm1,&
          h2xy(1,1,4,1),norbt,&
          pv(1,1,nxy,4,1),nasht,&
          d1,fq(1,1,1),norbt)

  end if

! ------------------
! Term 5 (real part)
! ------------------

  if (ix .ne. iy) then

!    \sum_{x>y} \sum_v \Twoint{pv|xy}_{1,2} P_{qv,xy;1,2}

     call dgemm('n','t',norbt,nasht,nasht,d1,&
          h2xy(1,1,1,2),norbt,&
          pv(1,1,nxy,1,2),nasht,&
          d1,fq(1,1,1),norbt)

  end if

! ------------------
! Term 6 (real part)
! ------------------

  if (nz .ge. 2 .and. ix .ne. iy) then

!    - \sum_{x>y} \sum_v \Twoint{pv|xy}_{2,2} P_{qv,xy;2,2}

     call dgemm('n','t',norbt,nasht,nasht,dm1, &
          h2xy(1,1,2,2),norbt,                 &
          pv(1,1,nxy,2,2),nasht,               &
          d1,fq(1,1,1),norbt)

  end if

  if (nz .ge. 4) then

!    ------------------
!    Term 7 (real part)
!    ------------------

!    \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{3,2} P_{qv,xy;3,2}

     call dgemm('n','t',norbt,nasht,nasht,d1, &
          h2xy(1,1,3,2),norbt,                &
          pv(1,1,nxy,3,2),nasht,              &
          d1,fq(1,1,1),norbt)

     if (ix .ne. iy) then

!       \sum_{x>y} \sum_v \Twoint{vp|xy}_{3,3} P_{vq,xy;3,3}

        call dgemm('n','n',norbt,nasht,nasht,d1, &
             h2yx(1,1,3,3),norbt,                &
             pv(1,1,nxy,3,3),nasht,              &
             d1,fq(1,1,1),norbt)

     end if

!    ------------------
!    Term 8 (real part)
!    ------------------

!    - \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{4,2} P_{qv,xy;4,2}

     call dgemm('n','t',norbt,nasht,nasht,dm1, &
          h2xy(1,1,4,2),norbt,                 &
          pv(1,1,nxy,4,2),nasht,               &
          d1,fq(1,1,1),norbt)

     if (ix .ne. iy) then

!       - \sum_{x>y} \sum_v \Twoint{vp|xy}_{4,3} P_{vq,xy;4,3}

        call dgemm('n','n',norbt,nasht,nasht,dm1, &
             h2yx(1,1,4,3),norbt,                 &
             pv(1,1,nxy,4,3),nasht,               &
             d1,fq(1,1,1),norbt)

     end if

  end if

! ------------------
! Term 9 (real part)
! ------------------

  if (ix .ne. iy) then

!    \sum_{x>y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;1,3}

     call dgemm('n','t',norbt,nasht,nasht,d1, &
          h2xy(1,1,1,3),norbt,                &
          pv(1,1,nxy,1,3),nasht,              &
          d1,fq(1,1,1),norbt)

  end if

! ------------------
! Term 10 (real part)
! ------------------

  if (nz .ge. 2 .and. ix .ne. iy) then

!    - \sum_{x>y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;2,3}

     call dgemm('n','t',norbt,nasht,nasht,dm1, &
          h2xy(1,1,2,3),norbt,                 &
          pv(1,1,nxy,2,3),nasht,               &
          d1,fq(1,1,1),norbt)

  end if

  if (nz .ge. 4 .and. ix .ne. iy) then

!    ------------------
!    Term 11 (real part)
!    ------------------

!    \sum_{x>y} \sum_v \Twoint{vp|xy}_{3,1} P_{vq,xy;3,1}

     call dgemm('n','n',norbt,nasht,nasht,d1, &
          h2yx(1,1,3,1),norbt,                &
          pv(1,1,nxy,3,1),nasht,              &
          d1,fq(1,1,1),norbt)

!    ------------------
!    Term 12 (real part)
!    ------------------

!    - \sum_{x>y} \sum_v \Twoint{vp|xy}_{4,1} P_{vq,xy;4,1}

     call dgemm('n','n',norbt,nasht,nasht,dm1, &
          h2yx(1,1,4,1),norbt,                 &
          pv(1,1,nxy,4,1),nasht,               &
          d1,fq(1,1,1),norbt)

  end if

! Skip remaining 36 terms if NZ = 1.

  if (nz .eq. 1) goto 100

! *************************
! *** i-imag part of FQ ***
! *************************
! (B.41)

! --------------------
! Term 1 (i-imag part)
! --------------------

! \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;2,1}

  call dgemm('n','t',norbt,nasht,nasht,d1, &
       h2xy(1,1,1,1),norbt,                &
       pv(1,1,nxy,2,1),nasht,              &
       d1,fq(1,1,2),norbt)

  if (ix .ne. iy) then

!    - \sum_{x>y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;2,1}

     call dgemm('n','n',norbt,nasht,nasht,dm1, &
          h2yx(1,1,1,1),norbt,                 &
          pv(1,1,nxy,2,1),nasht,               &
          d1,fq(1,1,2),norbt)

  end if

! --------------------
! Term 2 (i-imag part)
! --------------------


! + \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;1,1}

  call dgemm('n','t',norbt,nasht,nasht,d1, &
       h2xy(1,1,2,1),norbt,                &
       pv(1,1,nxy,1,1),nasht,              &
       d1,fq(1,1,2),norbt)

  if (ix .ne. iy) then

!    - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;1,1}

     call dgemm('n','n',norbt,nasht,nasht,dm1, &
          h2yx(1,1,2,1),norbt,                 &
          pv(1,1,nxy,1,1),nasht,               &
          d1,fq(1,1,2),norbt)

  end if

! --------------------------
! Term 3 and 4 (i-imag part)
! --------------------------

  if (nz .ge. 4 .and. ix .ne. iy) then

!    \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;4,1}

     call dgemm('n','t',norbt,nasht,nasht,d1, &
          h2xy(1,1,3,1),norbt,                &
          pv(1,1,nxy,4,1),nasht,              &
          d1,fq(1,1,2),norbt)

!    \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;3,1}

     call dgemm('n','t',norbt,nasht,nasht,d1, &
          h2xy(1,1,4,1),norbt,                &
          pv(1,1,nxy,3,1),nasht,              &
          d1,fq(1,1,2),norbt)

  end if

  if (ix .ne. iy) then

!    ------------------
!    Term 5 (i-imag part)
!    ------------------

!    \sum_{x>y} \sum_v \Twoint{pv|xy}_{1,2} P_{qv,xy;2,2}

     call dgemm('n','t',norbt,nasht,nasht,d1, &
          h2xy(1,1,1,2),norbt,                &
          pv(1,1,nxy,2,2),nasht,              &
          d1,fq(1,1,2),norbt)

!    --------------------
!    Term 6 (i-imag part)
!    --------------------

!    \sum_{x>y} \sum_v \Twoint{pv|xy}_{2,2} P_{qv,xy;1,2}

     call dgemm('n','t',norbt,nasht,nasht,d1, &
          h2xy(1,1,2,2),norbt,                &
          pv(1,1,nxy,1,2),nasht,              &
          d1,fq(1,1,2),norbt)

  end if

  if (nz .ge. 4) then

!    --------------------
!    Term 7 (i-imag part)
!    --------------------

!    \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{3,2} P_{qv,xy;4,2}

     call dgemm('n','t',norbt,nasht,nasht,d1, &
          h2xy(1,1,3,2),norbt,                &
          pv(1,1,nxy,4,2),nasht,              &
          d1,fq(1,1,2),norbt)

     if (ix .ne. iy) then

!       - \sum_{x>y} \sum_v \Twoint{vp|xy}_{3,3} P_{vq,xy;4,3}

        call dgemm('n','n',norbt,nasht,nasht,dm1, &
             h2yx(1,1,3,3),norbt,                 &
             pv(1,1,nxy,4,3),nasht,               &
             d1,fq(1,1,2),norbt)

     end if

!    --------------------
!    Term 8 (i-imag part)
!    --------------------

!    \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{4,2} P_{qv,xy;3,2}

     call dgemm('n','t',norbt,nasht,nasht,d1, &
          h2xy(1,1,4,2),norbt,                &
          pv(1,1,nxy,3,2),nasht,              &
          d1,fq(1,1,2),norbt)

     if (ix .ne. iy) then

!       - \sum_{x>y} \sum_v \Twoint{vp|xy}_{4,3} P_{vq,xy;3,3}

        call dgemm('n','n',norbt,nasht,nasht,dm1, &
             h2yx(1,1,4,3),norbt,                 &
             pv(1,1,nxy,3,3),nasht,               &
             d1,fq(1,1,2),norbt)

     end if

  end if

  if (ix .ne. iy) then

!    --------------------
!    Term 9 (i-imag part)
!    --------------------

!    - \sum_{x>y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;2,3}

     call dgemm('n','t',norbt,nasht,nasht,dm1, &
          h2xy(1,1,1,3),norbt,                 &
          pv(1,1,nxy,2,3),nasht,               &
          d1,fq(1,1,2),norbt)

!    ---------------------
!    Term 10 (i-imag part)
!    ---------------------

!    - \sum_{x>y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;1,3}

     call dgemm('n','t',norbt,nasht,nasht,dm1, &
          h2xy(1,1,2,3),norbt,                 &
          pv(1,1,nxy,1,3),nasht,               &
          d1,fq(1,1,2),norbt)

  end if

  if (nz .ge. 4 .and. ix .ne. iy) then

!    ---------------------
!    Term 11 (i-imag part)
!    ---------------------

!    \sum_{x>y} \sum_v \Twoint{vp|xy}_{3,1} P_{vq,xy;4,1}

     call dgemm('n','n',norbt,nasht,nasht,d1, &
          h2yx(1,1,3,1),norbt,                &
          pv(1,1,nxy,4,1),nasht,              &
          d1,fq(1,1,2),norbt)

!    ---------------------
!    Term 12 (i-imag part)
!    ---------------------

!    - \sum_{x>y} \sum_v \Twoint{vp|xy}_{4,1} P_{vq,xy;3,1}

     call dgemm('n','n',norbt,nasht,nasht,dm1, &
          h2yx(1,1,4,1),norbt,                 &
          pv(1,1,nxy,3,1),nasht,               &
          d1,fq(1,1,2),norbt)

  end if

! Skip remaining 24 terms if NZ = 2.

  if (nz .eq. 2) goto 100

  call quit('*** ERROR in raddfq90 *** nz = 4 not implemented')

! --------------
! Output section
! --------------

 100  continue
  if (iprint .ge. 20) then
    write(lupri,'(a)') 'fq on exit'
    call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  end if

  !------------------------------------------------------------------------------!
  !           Write to the file ouput information for unit testing               !
  !------------------------------------------------------------------------------!
  if( present(key) .or. create_outputs )then
    ! open file
    open(60,file='raddfq90.out_' // i2c(nxy) //'_'// i2c(ix) //'_'// i2c(iy), &
         status='unknown',form='formatted',access='sequential',iostat=ier)
    if(ier /= 0)then
      write(6,*) 'raddfq90: ERROR while creating raddfq90.uout file'
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
  endif

  end subroutine raddfq90

end module mcscf_routines
