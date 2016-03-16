program test_davidson_lille 

   use cdavidson_lille

   implicit none

   double precision, allocatable :: m_real_packed(:)
   double precision, allocatable :: m_real_full(:,:)
   integer :: mdim

   mdim = 5

   allocate(m_real_packed(mdim*(mdim-1)/2))
   call setup_real_symmetric_matrix_packed(m_real_packed, mdim)
   print *, "test, packed storage matrix"
   call print_real_matrix_packed(m_real_packed, mdim)

   allocate(m_real_full(mdim,mdim))
   call setup_real_symmetric_matrix_full(m_real_full, mdim)
   print *, "test, full matrix"
   call print_real_matrix_full(m_real_full, mdim)

contains

   complex function delta_complex(i,j)
      integer :: i, j
      if (i .eq. j) then
         delta_complex = (1.0d0, 0.0d0)
      else
         delta_complex = (0.0d0, 0.0d0)
      endif
   end function

   double precision function delta_real(i,j)
      integer :: i, j
      if (i .eq. j) then
         delta_real = 1.0d0
      else
         delta_real = 0.0d0
      endif
   end function

   subroutine setup_real_symmetric_matrix_packed(m, mdim)
      integer :: i, j, k, ij
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:)

      k = mdim/2

      ij = 0
      m = 0.0d0
      do i = 1, mdim
         do j = i, mdim
            ij = (i-1)*mdim + j 
            m(ij) = i*delta_real(i,j) + (i - j - k*k)
         enddo
      enddo
   end subroutine
       
   subroutine setup_real_symmetric_matrix_full(m, mdim)
      integer :: i, j, k
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:,:)

      k = mdim/2

      m = 0.0d0
      do i = 1, mdim
         do j = i, mdim
            m(i,j) = i*delta_real(i,j) + (i - j - k*k)
            m(j,i) = m(i,j)
         enddo
      enddo
   end subroutine

   subroutine print_real_matrix_packed(m, mdim)
      integer :: i, j, ij
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:)

      do i = 1, mdim
         do j = 1, mdim
            ij = (i-1)*mdim + j 
            print "(4x,1f14.6,\)", m(ij)
         enddo
         print *, ""
      enddo
   end subroutine

   subroutine print_real_matrix_full(m, mdim)
      integer :: i, j
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:,:)

      do i = 1, mdim
         do j = 1, mdim
            print "(4x,1f14.6,\)", m(i,j)
         enddo
         print *, ""
      enddo

   end subroutine

!   call eiprmn(eso,hd,vr,vi,n,mvec,ier,vecfil)

   subroutine davidson_driver_packed(m, mdim)
      integer :: i, j
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:)
      double precision, allocatable :: trial_vectors(:,:)

!     call cdeigen(eso,hd,trial_vectors,mdim,mvec,ier,vecfil) 

   end subroutine

!  subroutine sigma_vectors_builder(sigma_vectors, m, mdim)

!  end subroutine


end program test_davidson_lille 
