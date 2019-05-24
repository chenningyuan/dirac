module utils

contains

   subroutine write_dffiles(filename_unformatted,filename_formatted,coeff,basis1,basis2)

       implicit none

       integer, parameter :: nsym = 1, iunit = 98, lucoef = 97
       character(len=74) :: text="IAO coefficients"
       character(10) :: filename_unformatted, filename_formatted
       integer :: idim(3,nsym), i,j
       integer, intent(in)     :: basis1, basis2
       complex(8), intent(in)  :: coeff(basis1,basis2,2)
       real(8)  :: coeff_quaternion(basis1,int(basis2/2),4)
       real(8) :: toterg
       real(8) :: eig(basis2)
       integer :: ibeig(basis2)

       toterg= 0.D0
       eig = 0.D0
       ibeig = 0

       idim(3,1) = basis1
       idim(1,1) = 0
       idim(2,1) = int(basis2/2)

       ! I think that I need to have this format for the DFCOEF file ?
       call convert_complex_to_quaternion_format(coeff,coeff_quaternion,basis1,basis2)

       open (UNIT = iunit,file=filename_unformatted,FORM='UNFORMATTED',access='SEQUENTIAL')
       write (iunit) text,nsym,((idim(i,j),i=1,3),j=1,nsym),toterg
       write (iunit) coeff_quaternion
       write (iunit) eig
       write (iunit) ibeig
       close (iunit)

       open (UNIT = iunit,file=filename_formatted,FORM='FORMATTED',access='SEQUENTIAL')
       write(iunit,'(A74)') text
       write(iunit,'(7(X,I0))') nsym,((idim(i,j),i = 1,3),j=1,nsym)
       write(iunit,'(E24.16)') toterg
       write(iunit,'(6F22.16)') coeff_quaternion
       write(iunit,'(6E22.12)') eig
       write(iunit,'(66(X,I0))') ibeig


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

end module utils
