module utils

contains

   subroutine write_dffiles(filename_unformatted,filename_formatted,ao_basis,coeff,nao,nmo)

       use exacorr_datatypes

       implicit none

       integer, parameter :: nsym = 1, iunit = 98
       character(len=74) :: text="IAO coefficients"
       character(len=*) :: filename_unformatted, filename_formatted
       integer :: idim(3,nsym), i, j
       integer, intent(in) :: nao, nmo
       complex(8), intent(in)  :: coeff(nao,nmo,2)
       type(basis_set_info_t), intent(in) :: ao_basis
       real(8)  :: coeff_quaternion(nao,int(nmo/2),4)
       real(8) :: toterg
       real(8) :: eig(nmo)
       integer :: ibeig(nmo)

       toterg= 0.D0
       eig = 0.D0
       ibeig = 0

       idim(3,1) = nao
       idim(1,1) = 0
       idim(2,1) = int(nmo/2)

       ! I think that I need to have this format for the DFCOEF file ?
       call convert_complex_to_quaternion_format(coeff,coeff_quaternion,nao,nmo)

       open (UNIT = iunit,file=filename_unformatted,FORM='UNFORMATTED',access='SEQUENTIAL')
       write (iunit) text,nsym,((idim(i,j),i=1,3),j=1,nsym),toterg
       write (iunit) coeff_quaternion
       write (iunit) eig
       write (iunit) ibeig
     
       call write_basis_dfcoef(ao_basis,"MOLECULE.dfcoef",iunit)
       close(iunit)

       open (UNIT = iunit,file=filename_formatted,FORM='FORMATTED',access='SEQUENTIAL')
       write(iunit,'(A74)') text
       write(iunit,'(7(X,I0))') nsym,((idim(i,j),i = 1,3),j=1,nsym)
       write(iunit,'(E24.16)') toterg
       write(iunit,'(6F22.16)') coeff_quaternion
       write(iunit,'(6E22.12)') eig
       write(iunit,'(66(X,I0))') ibeig
       close(iunit)

   end subroutine write_dffiles

   subroutine convert_complex_to_quaternion_format(scoeff,coeff,nao,nmo)

      implicit none
      integer,    intent(in   ) :: nao, nmo
      complex(8), intent(in   ) :: scoeff(nao,nmo,2)
      real(8),    intent(  out) :: coeff(nao,int(nmo/2),4)
      integer                   :: i, j, spin

      do spin = 1, 2
        if (spin == 1) then
            do j = 1, int(nmo/2)
                do i = 1, nao
                coeff(i,j,1) = real(scoeff(i,2*j-1,spin))
                coeff(i,j,2) = aimag(scoeff(i,2*j-1,spin))
                coeff(i,j,3) = real(scoeff(i,2*j,spin))
                coeff(i,j,4) = aimag(scoeff(i,2*j,spin))
                end do
            end do
!        else if (spin == 2) then
!            do j = 1, int(nmo/2)
!                do i = 1, nao
!                coeff(i,j,3) = - real(scoeff(i,2*j-1,spin))
!                coeff(i,j,4) = aimag(scoeff(i,2*j-1,spin))
!                coeff(i,j,1) = real(scoeff(i,2*j,spin))
!                coeff(i,j,2) = - aimag(scoeff(i,2*j,spin))
!                end do
!            end do
        end if
      end do

   end subroutine convert_complex_to_quaternion_format

   subroutine diag_matrix(A,n,eigval)
        implicit none
        integer,    intent(in)       :: n
        integer,    parameter        :: lwmax=10000
        integer                      :: info, lwork, liwork, lrwork
        integer,    dimension(lwmax) :: iwork
        complex(8), dimension(lwmax) :: work
        real(8),    dimension(lwmax) :: rwork
        complex(8), intent(inout)    :: A(n,n)
        real(8),    intent(  out)    :: eigval(n)
!
!       Query the optimal workspace.
!
        lwork = -1
        liwork = -1
        lrwork = -1

        call zheevd( 'V', 'L', n, A, n, eigval, work, lwork, rwork, &
                    lrwork, iwork, liwork, info )
        lwork = min(lwmax,int(work(1)))
        lrwork = min(lwmax,int(rwork(1)))
        liwork = min(lwmax,iwork(1))
!
!       Solve eigenproblem.
!
        call zheevd( 'V', 'L', n, A, n, eigval, work, lwork, rwork, &
                    lrwork, iwork, liwork, info )
!
!       Check for convergence.
!
        if( info.gt.0 ) then
           write(*,*)'The algorithm failed to compute eigenvalues.'
           stop
        end if

   end subroutine diag_matrix

end module utils
