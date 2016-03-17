program test_davidson_lille 

   use cdavidson_lille

   implicit none

! eigenvectors, eigenvalues 
   double precision, allocatable :: eval_real(:)
   double precision, allocatable :: evec_real(:,:)
   double precision, allocatable :: evec_real_l(:,:)
   double precision, allocatable :: evec_real_r(:,:)

   double precision, allocatable :: eval_r(:)
   double precision, allocatable :: eval_i(:)
   complex (kind=8), allocatable :: eval_complex(:)
   complex (kind=8), allocatable :: evec_complex(:,:)
   complex (kind=8), allocatable :: evec_complex_r(:,:)
   complex (kind=8), allocatable :: evec_complex_l(:,:)

! matrices to diagolize
   double precision, allocatable :: m_real_symmetric_packed(:)
   double precision, allocatable :: m_real_symmetric_full(:,:)
   double precision, allocatable :: m_real_nonsymmetric_full(:,:)

   complex (kind=8), allocatable :: m_complex_hermitian_full(:,:)
   complex (kind=8), allocatable :: m_complex_nonsymmetric_full(:,:)
!
   logical :: diagonalize_real_symm_packed  = .false.
   logical :: diagonalize_real_symm_full    = .false.
   logical :: diagonalize_complex_herm_full = .true.
   logical :: diagonalize_real_nonsymm_full = .false.
   logical :: diagonalize_complex_nonsymm_full = .false.
!
   integer :: mdim

   mdim = 5

! allocating eigenvectors and eigenvalues storage
   allocate(eval_real(mdim))
   allocate(eval_r(mdim))
   allocate(eval_i(mdim))
   allocate(eval_complex(mdim))
   eval_real = 0.0d0 
   eval_r = 0.0d0 
   eval_i = 0.0d0 
   eval_complex = cmplx(0.0d0, 0.0d0)

   allocate(evec_real(mdim,mdim))
   allocate(evec_real_l(mdim,mdim))
   allocate(evec_real_r(mdim,mdim))

   evec_real = 0.0d0 
   evec_real_l = 0.0d0 
   evec_real_r = 0.0d0 

   allocate(evec_complex(mdim,mdim))
   allocate(evec_complex_l(mdim,mdim))
   allocate(evec_complex_r(mdim,mdim))
   evec_complex = cmplx(0.0d0, 0.0d0)
   evec_complex_l = cmplx(0.0d0, 0.0d0)
   evec_complex_r = cmplx(0.0d0, 0.0d0)

! symmetric / hermitian matrices

   allocate(m_real_symmetric_packed(mdim*(mdim-1)/2))

   call setup_real_symmetric_matrix_packed(m_real_symmetric_packed, mdim)
   print *, "test, packed storage real (symmetric) matrix"
   call print_real_matrix_packed(m_real_symmetric_packed, mdim)

   if (diagonalize_real_symm_packed) then
      call diagonalize_real_symmetric_matrix_packed(m_real_symmetric_packed,eval_real,evec_real)
      print *, "eigenvalues of matrix real symmetric matrix (packed)"
      call print_real_eigenvalues(eval_real, mdim)
   endif

   allocate(m_real_symmetric_full(mdim,mdim))
   call setup_real_symmetric_matrix_full(m_real_symmetric_full, mdim)
   print *, "test, full real symmetric matrix"
   call print_real_matrix_full(m_real_symmetric_full, mdim)

   if (diagonalize_real_symm_full) then
      evec_real = m_real_symmetric_full
!     call diagonalize_real_symmetric_matrix_full(evec_real,eval_real)
      call diagonalize_real_general_matrix_full(m_real_symmetric_full,eval_r,eval_i,evec_real_r,evec_real_l)
      print *, "eigenvalues of real symmetric matrix (full)"
!     call print_real_eigenvalues(eval_real, mdim)
      call print_real_eigenvalues(eval_r, mdim)
   endif

! non-symmetric matrix

   allocate(m_real_nonsymmetric_full(mdim,mdim))
   call setup_real_nonsymmetric_matrix_full(m_real_nonsymmetric_full, mdim)
   print *, "test, full real non-symmetric matrix"
   call print_real_matrix_full(m_real_nonsymmetric_full, mdim)

!  call diagonalize_real_symmetric_matrix(m_real_nonsymmetric_full,mdim)
   

! complex matrices

   print *, ""
   print *, ""

   allocate(m_complex_hermitian_full(mdim,mdim))
   call setup_complex_hermitian_matrix_full(m_complex_hermitian_full, mdim)
   print *, "test, full complex hermitian matrix"
   call print_complex_matrix_full(m_complex_hermitian_full, mdim)

   if (diagonalize_complex_herm_full) then
      call diagonalize_complex_general_matrix_full(m_complex_hermitian_full,eval_complex,evec_complex_r,evec_complex_l)
      print *, "eigenvalues of hermitian matrix (full)"
!     call print_complex_eigenvalues(eval_complex, mdim)
   endif


   allocate(m_complex_nonsymmetric_full(mdim,mdim))
   call setup_complex_nonsymmetric_matrix_full(m_complex_nonsymmetric_full, mdim)
   print *, "test, full complex non-symmetric matrix"
   call print_complex_matrix_full(m_complex_nonsymmetric_full, mdim)

contains

   complex function delta_complex(i,j)
      integer :: i, j
      if (i .eq. j) then
         delta_complex = cmplx(1.0d0, 0.0d0)
      else
         delta_complex = cmplx(0.0d0, 0.0d0)
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
            print *, i, j, ij, m(ij)
         enddo
      enddo
   end subroutine
       
   subroutine setup_real_symmetric_matrix_full(m, mdim)
      integer :: i, j, k
      integer :: ii, jj
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:,:)

      k = mdim/2

      m = 0.0d0
      do i = 1, mdim
         ii = i - 1
         do j = i, mdim
            jj = j - 1

            m(i,j) = ii*delta_real(ii,jj) + (ii - jj - k*k)
            m(j,i) = m(i,j)
         enddo
      enddo
   end subroutine

   subroutine setup_complex_hermitian_matrix_full(m, mdim)
      integer :: i, j, k
      integer :: ii, jj
      integer, intent(in) :: mdim
      real (kind=8) :: real_part, imag_part
      complex (kind=8), intent(inout) :: m(:,:)

      k = mdim/2

      m = cmplx(0.0d0, 0.0d0)
      do i = 1, mdim
         ii = i - 1
         do j = i, mdim
            jj = j - 1

            real_part = ii*delta_real(i,j) - (ii - jj - k*k)
            imag_part = ii*delta_real(i,j) - (ii - jj)
            m(i,j) = cmplx(real_part, imag_part)
            m(j,i) = conjg(m(i,j))
         enddo
      enddo
   end subroutine

  subroutine setup_real_nonsymmetric_matrix_full(m, mdim)
      integer :: i, j, k
      integer :: ii, jj
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:,:)
    
      k = mdim/2 
            
      m = 0.0d0
      do i = 1, mdim 
         ii = i - 1
         do j = 1, mdim 
            jj = j - 1
            if (jj .le. k) then
               m(i,j) = ii*delta_real(i,j) - (ii - jj - k*k)
            else if ((jj .le. mdim).and.(jj .ge. (k+1))) then
               m(i,j) = ii*delta_real(i,j) + (ii - jj - k*k)
            endif
         enddo
      enddo
   end subroutine

  subroutine setup_complex_nonsymmetric_matrix_full(m, mdim)
      integer :: i, j, k
      integer :: ii, jj
      integer, intent(in) :: mdim
      real (kind=8) :: real_part, imag_part
      complex (kind=8), intent(inout) :: m(:,:)
    
      k = mdim/2 
            
      m = cmplx(0.0d0,  0.0d0)
      do i = 1, mdim 
         ii = i - 1
         do j = 1, mdim 
            jj = j - 1

            if (jj .le. k) then
               real_part = ii*delta_real(i,j) - (ii - jj - k*k)
               imag_part = ii*delta_real(i,j) - (ii - jj)
               m(i,j) = cmplx(real_part, imag_part)
            else if ((jj .le. mdim).and.(jj .ge. (k+1))) then
               real_part = i*delta_real(i,j) + (ii - jj - k*k)
               imag_part = i*delta_real(i,j) - (ii - jj)
               m(i,j) = cmplx(real_part, imag_part)
            endif
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

   subroutine print_real_eigenvalues(e, n)
      integer :: i, n
      double precision, intent(in) :: e(:)

      do i = 1, n
         print "(4x,1f14.6,\)", e(i)
      enddo
      print *, ""
   end subroutine

   subroutine print_complex_matrix_full(m, mdim)
      integer :: i, j
      integer, intent(in) :: mdim
      complex (kind=8), intent(inout) :: m(:,:)
      character(4) :: imag_sign = "   +"

      do i = 1, mdim
         do j = 1, mdim
            if (aimag(m(i,j)) .ge. 0.0d0) then
                imag_sign = "   +" 
            else
                imag_sign = "   -" 
            endif

            print "(4x,1f14.6,a,1f14.6,a,\)", real(m(i,j)),imag_sign,abs(aimag(m(i,j)))," i"
         enddo
         print *, ""
      enddo

   end subroutine


   subroutine diagonalize_real_symmetric_matrix_packed(m,eval,evec)

      use lapack95, only: spev

      character :: jobz = 'V'
      character :: uplo = 'U'
      integer   :: info
      double precision, intent(inout) :: m(:)  
      double precision, intent(out) :: eval(:)
      double precision, intent(out) :: evec(:,:)

      call spev(m, eval, uplo, evec, info)
      if (info .ne. 0) then
         print *, 'error in diagonalize_real_symmetric_matrix_packed, ',info
         call exit(info)
      endif
   end subroutine

   subroutine diagonalize_real_symmetric_matrix_full(m,eval)

      use lapack95, only: syev

!     character :: jobz = 'V'
!     character :: uplo = 'U'
!     integer   :: info
      double precision, intent(inout) :: m(:,:)  
      double precision, intent(out) :: eval(:)

!     call syev(m, eval)
!     call syev(m, eval, jobz, uplo, info)
!     if (info .ne. 0) then
!        print *, 'error in diagonalize_real_symmetric_matrix_full, ',info
!        call exit(info)
!     endif
   end subroutine


   subroutine diagonalize_real_general_matrix_full(m,eval_r,eval_i,evec_r,evec_l)

      use lapack95, only: geev

      integer   :: info
      double precision, intent(inout) :: m(:,:)
      double precision, intent(out) :: eval_r(:), eval_i(:)
      double precision, intent(out) :: evec_r(:,:), evec_l(:,:)

      call geev(m, eval_r, eval_i, evec_l, evec_r, info)
      print *, m
      print *,eval_r
      print *,eval_i
      if (info .ne. 0) then
         print *, 'error in diagonalize_real_general_matrix, ',info
         call exit(info)
      endif

   end subroutine 

   subroutine diagonalize_complex_general_matrix_full(m,eval,evec_r,evec_l)

      use lapack95, only: geev

      integer   :: info
      complex (kind=8), intent(inout) :: m(:,:)
      complex (kind=8), intent(out) :: eval(:)
      complex (kind=8), intent(out) :: evec_r(:,:), evec_l(:,:)

      call geev(m, eval, evec_l, evec_r, info)

      print *, m
      print *,eval

      if (info .ne. 0) then
         print *, 'error in diagonalize_complex_general_matrix, ',info
         call exit(info)
      endif

   end subroutine


   subroutine diagonalize_complex_hermitian_matrix_packed(m,eval,evec)

      use lapack95, only: hpev

      character :: jobz = 'V'
      character :: uplo = 'U'
      integer   :: info
      complex (kind=8), intent(inout) :: m(:)  
      double precision, intent(out) :: eval(:)
      complex (kind=8), intent(out) :: evec(:,:)

      call hpev(m, eval, uplo, evec, info)
      if (info .ne. 0) then
         print *, 'error in diagonalize_complex_hermitian_matrix_packed, ',info
         call exit(info)
      endif
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
