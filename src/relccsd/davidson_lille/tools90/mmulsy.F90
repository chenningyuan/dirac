module mmulsy
!
!*******************************************************************************
!     ----- this module contains subroutine to compute the matrix -----
!     ----- product op(A)*B = C where the resulting matrix C is   -----
!     ----- symmetric or hermitian and stored in upper triangular -----
!     ----- packed format                                         -----
!-------------------------------------------------------------------------------
!     usage:
!     -----
!     call mmsy(a,b,c,nr,nc)
!
!  * A,B are (nr,nc) matrices (INPUT)
!  * C  is a matrix of rank nc stored in upper triangular packed format (OUTPUT)
!    if the matrices are real, op(A) = transpose(A)
!    if the matrices are complex, op(A) = conjg(transpose(A))
!  * NR = number of rows of A and B (OPTIONAL)
!  * NC = number of columns of A and B and rank of C (OPTIONAL)
!
!===============================================================================
!     written by: J.-P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
   use kindvar
   use mconst
   use mdebug
   use f95_precision, only: wp => dp
   use blas95, only: dot,gemm
   implicit none
!
!     ----- parameters -----
!  
      integer (kind=K_I4), parameter :: SH=5
      integer (kind=K_I4), parameter :: NB=2**SH
   interface mmsy
      module procedure mmsy_r8,mmsy_c81,mmsy_c82
   end interface mmsy
   contains
   subroutine mmsy_r8(a,b,c,nr,nc)
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!===============================================================================
!     written by: J.-P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
!     ----- interfaces -----
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in), optional :: nr,nc
      real (kind=K_R8), intent(in), dimension (:,:) :: a,b
      real (kind=K_R8), intent(out), dimension (:) :: c
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: mr,mc
!
!     ----- allocatable variables -----
!
!     ----- formats -----
!
!
!     +++++
!     start
!     +++++
!
      if(present(nr)) then
         mr = nr
      else
         mr = size(a,1)
      endif
      if(present(nc)) then
         mc = nc
      else
         mc = size(a,2)
      endif
      call vhvsr(a(1:mr,1:mc),b(1:mr,1:mc),c,mc)
      return
   end subroutine mmsy_r8
   subroutine mmsy_c81(a,b,c,nr,nc)
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!===============================================================================
!     written by: J.-P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
!     ----- interfaces -----
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in), optional :: nr,nc
      complex (kind=K_C8), intent(in), dimension (:,:) :: a,b
      complex (kind=K_C8), intent(out), dimension (:) :: c
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: mr,mc
!
!     ----- allocatable variables -----
!
!     ----- formats -----
!
!
!     +++++
!     start
!     +++++
!
      if(present(nr)) then
         mr = nr
      else
         mr = size(a,1)
      endif
      if(present(nc)) then
         mc = nc
      else
         mc = size(a,2)
      endif
      call vhvsc(a(1:mr,1:mc),b(1:mr,1:mc),c,mc)
      return
   end subroutine mmsy_c81
   subroutine mmsy_c82(ar,ai,br,bi,c,nr,nc)
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!===============================================================================
!     written by: J.-P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
!     ----- interfaces -----
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in), optional :: nr,nc
      real (kind=K_R8), intent(in), dimension (:,:) :: ar,ai,br,bi
      complex (kind=K_C8), intent(out), dimension (:) :: c
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: mr,mc
!
!     ----- allocatable variables -----
!
!     ----- formats -----
!
!
!     +++++
!     start
!     +++++
!
      if(present(nr)) then
         mr = nr
      else
         mr = size(ar,1)
      endif
      if(present(nc)) then
         mc = nc
      else
         mc = size(ar,2)
      endif
!     CALL DBGPRT('C82: MR MC=',(/ MR,MC /),2)
      call vhvsc2(ar(1:mr,1:mc),ai(1:mr,1:mc),br(1:mr,1:mc),  &
                  bi(1:mr,1:mc),c,mc)
      return
   end subroutine mmsy_c82
   subroutine vhv3c(vhv,hvc,vc,n,nprev,ncor)
!
!***********************************************************************
!     ----- compute the <v|H|v> matrix -----
!=======================================================================
!       VHV (C8) = lower triangle <v|H|v> (OUTPUT)
!        HV (C8) = H*|v> matrix (INPUT)
!         V (C8) = vectors v(n,nc) (INPUT)
!         N (I4) = row dimension of VR and HVR (INPUT)
!     NPREV (I4) = number of vectors |v> allready considered (INPUT)
!      NCOR (I4) = number of corrections vectors |v> to be added (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer, intent(in) :: n,nprev,ncor
      complex (kind=K_C8), intent(inout), dimension (:) :: vhv
      complex (kind=K_C8), intent(in), dimension (:,:) :: hvc,vc
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: j,ij,ij1,ij2,n1,n2,ios,i,k,k2,i2,i3,kend,iend,k3
      complex (kind=K_C8) :: c11,c21,c22,c12
      complex (kind=K_C8) :: a11,a21,a22,a12
      complex (kind=K_C8) :: b11,b21,b22,b12
!
!     ----- allocatable variables -----
!
      complex (kind=K_C8), allocatable, dimension (:,:) :: temp
!
!     +++++
!     start
!     +++++
!
      n1 = nprev+1
      n2 = nprev+ncor
      ij1 = (nprev*(nprev+1))/2
      if(nprev.gt.0) then
         allocate(temp(nprev,ncor),stat=ios)
         if(ios.ne.0) then
            write(6,*) 'alloc error in VHV3C',ios
         endif
         call gemm(vc(:,1:nprev),hvc(:,n1:n2),temp,'C','N',CONE,CZERO)
         ij = ij1
         do i=1,ncor
            vhv(ij+1:ij+nprev) = temp(1:nprev,i)
            ij = ij+i+nprev
         enddo
         deallocate(temp)
      endif
      i2 = iand(ncor,1)
      k2 = iand(n,1)
      kend = n-k2
      iend = n2-i2
      k3 = kend+1
      ij = ij1
      do i=n1,iend,2
         ij2 = ij+i
         do j=n1,i-2,2
            c11 = CZERO
            c12 = CZERO
            c21 = CZERO
            c22 = CZERO
            do k=1,kend,2
               a11 = conjg(vc(k,j))
               a12 = conjg(vc(k,j+1))
               b11 = hvc(k,i)
               b12 = hvc(k,i+1)
               c11 = c11+(a11*b11)
               c21 = c21+(a11*b12)
               c12 = c12+(a12*b11)
               c22 = c22+(a12*b12)
               a21 = conjg(vc(k+1,j))
               a22 = conjg(vc(k+1,j+1))
               b21 = hvc(k+1,i)
               b22 = hvc(k+1,i+1)
               c11 = c11+(a21*b21)
               c21 = c21+(a21*b22)
               c12 = c12+(a22*b21)
               c22 = c22+(a22*b22)
            enddo
            do k=k3,n
               a11 = conjg(vc(k,j))
               a12 = conjg(vc(k,j+1))
               b11 = hvc(k,i)
               b12 = hvc(k,i+1)
               c11 = a11*b11+c11
               c12 = a12*b11+c12
               c21 = a11*b12+c21
               c22 = a12*b12+c22
            enddo
            vhv(ij+j) = c11
            vhv(ij+j+1) = c12
            vhv(ij2+j) = c21
            vhv(ij2+j+1) = c22
         enddo
         j = i
         c11 = CZERO
         c21 = CZERO
         c22 = CZERO
         do k=1,kend,2
            a11 = conjg(vc(k,j))
            a21 = conjg(vc(k+1,j))
            a12 = conjg(vc(k,j+1))
            a22 = conjg(vc(k+1,j+1))
            b11 = hvc(k,i)
            b21 = hvc(k+1,i)
            b12 = hvc(k,i+1)
            b22 = hvc(k+1,i+1)
            c11 = c11+(a11*b11)+(a21*b21)
            c21 = c21+(a11*b12)+(a21*b22)
            c22 = c22+(a12*b12)+(a22*b22)
         enddo
         do k=k3,n
            a11 = conjg(vc(k,j))
            a12 = conjg(vc(k,j+1))
            b11 = hvc(k,i)
            b12 = hvc(k,i+1)
            c11 = a11*b11+c11
            c21 = a11*b12+c21
            c22 = a12*b12+c22
         enddo
         vhv(ij+j) = c11
         vhv(ij2+j) = c21
         vhv(ij2+j+1) = c22
         ij = ij2+i+1
      enddo
!
!     ----- clean up I -----
!
      i3 = iend+1
      do i=i3,n2
         ij2 = ij+i
         do j=n1,i-2,2
            c11 = CZERO
            c12 = CZERO
            c21 = CZERO
            c22 = CZERO
            do k=1,kend,2
               a11 = conjg(vc(k,j))
               a12 = conjg(vc(k,j+1))
               b11 = hvc(k,i)
!              b12 = hvc(k,i+1)
               c11 = c11+(a11*b11)
!              c21 = c21+(a11*b12)
               c12 = c12+(a12*b11)
!              c22 = c22+(a12*b12)
               a21 = conjg(vc(k+1,j))
               a22 = conjg(vc(k+1,j+1))
               b21 = hvc(k+1,i)
!              b22 = hvc(k+1,i+1)
               c11 = c11+(a21*b21)
!              c21 = c21+(a21*b22)
               c12 = c12+(a22*b21)
!              c22 = c22+(a22*b22)
            enddo
            do k=k3,n
               a11 = conjg(vc(k,j))
               a12 = conjg(vc(k,j+1))
               b11 = hvc(k,i)
!              b12 = hvc(k,i+1)
               c11 = a11*b11+c11
               c12 = a12*b11+c12
!              c21 = a11*b12+c21
!              c22 = a12*b12+c22
            enddo
            vhv(ij+j) = c11
            vhv(ij+j+1) = c12
!           vhv(ij2+j) = c21
!           vhv(ij2+j+1) = c22
         enddo
!        if(i.eq.n2) cycle
         j = i
         c11 = CZERO
!        c21 = CZERO
!        c22 = CZERO
         do k=1,kend,2
            a11 = conjg(vc(k,j))
            a21 = conjg(vc(k+1,j))
!           a12 = conjg(vc(k,j+1))
!           a22 = conjg(vc(k+1,j+1))
            b11 = hvc(k,i)
            b21 = hvc(k+1,i)
!           b12 = hvc(k,i+1)
!           b22 = hvc(k+1,i+1)
            c11 = c11+(a11*b11)+(a21*b21)
!           c21 = c21+(a11*b12)+(a21*b22)
!           c22 = c22+(a12*b12)+(a22*b22)
         enddo
         do k=k3,n
            a11 = conjg(vc(k,j))
!           a12 = conjg(vc(k,j+1))
            b11 = hvc(k,i)
!           b12 = hvc(k,i+1)
            c11 = a11*b11+c11
!           c21 = a11*b12+c21
!           c22 = a12*b12+c22
         enddo
         vhv(ij+j) = c11
!        vhv(ij2+j) = c21
!        vhv(ij2+j+1) = c22
         ij = ij2+i+1
      enddo
      return
   end subroutine vhv3c
   subroutine vhvsc(vc,hvc,vhv,nc)
!
!***********************************************************************
!     ----- compute the <v|H|v> matrix -----
!=======================================================================
!       VHV (C8) = lower triangle <v|H|v> (OUTPUT)
!        HV (C8) = H*|v> matrix (INPUT)
!         V (C8) = vectors v(n,nc) (INPUT)
!        NC (I4) = number of vectors |v> (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
!     use kindvar
!     use mconst
!     use mdebug
!     use f95_precision, only: wp => dp
!     use blas95, only: dot,gemm
      implicit none
!
!     ----- parameters -----
!  
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: nc
      complex (kind=K_C8), intent(out), dimension (:) :: vhv
      complex (kind=K_C8), intent(in), dimension (:,:) :: hvc,vc
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: i,j,ij,ii,jj,ij0,iend,iu,m,ib
      complex (kind=K_C8), dimension (NB,NB) :: temp
      complex (kind=K_C8), dimension (NB*(NB+1)/2) :: tempc
!
!     +++++
!     start
!     +++++
!
      ij0 = 0
      iend = NB*ishft(nc,-SH)
      m = size(vc,1)
      do ii=1,iend,NB
         do jj=1,ii-1,NB
            ij = ij0
            call gemm(vc(:,jj:jj+(NB-1)),hvc(:,ii:ii+(NB-1)),temp,  &
                      'C','N',CONE,CZERO)
            do i=1,NB
               vhv((ij+jj):(ij+jj)+(NB-1)) = temp(1:NB,i)
               ij = ij+(ii-1)+i
            enddo
         enddo
!
!     ----- diagonal bloc -----
!
         call vhv3c(tempc,hvc(:,ii:ii+(NB-1)),vc(:,ii:ii+(NB-1)),m,0,NB)
         ij = ij0
         iu = 1
         do i=1,NB
            do j=1,i
               vhv((ij+jj-1)+j) = tempc(iu)
               iu = iu+1
            enddo
            ij = ij+(ii-1)+i
         enddo
         ij0 = ij
      enddo
!
!      ----- clean up II bloc -----
!
      ib = nc-iend
      if(ib.eq.0) return
      ij0 = (ii*(ii-1))/2
      do jj=1,ii-1,NB
         call gemm(vc(:,jj:jj+(NB-1)),hvc(:,ii:nc),temp(:,1:(nc+1)-ii),  &
                   'C','N',CONE,CZERO)
         do i=1,ib
            vhv((ij+jj):(ij+jj)+(NB-1)) = temp(1:NB,i)
            ij = ij+(ii-1)+i
         enddo
         ij = ij0
      enddo
!
!     ----- diagonal bloc -----
!
      call vhv3c(tempc,hvc(:,ii:nc),vc(:,ii:nc),m,0,ib)
      ij = ij0
      iu = 1
      do i=1,ib
         do j=1,i
            vhv((ij+jj-1)+j) = tempc(iu)
            iu = iu+1
         enddo
         ij = ij+(ii-1)+i
      enddo
      return
   end subroutine vhvsc
   subroutine vhv3r(vhv,hvr,vr,n,nprev,ncor)
!
!***********************************************************************
!     ----- compute the <v|H|v> matrix -----
!=======================================================================
!       VHV (R8) = lower triangle <v|H|v> (OUTPUT)
!        HV (R8) = H*|v> matrix (INPUT)
!         V (R8) = vectors v(n,nc) (INPUT)
!         N (I4) = row dimension of VR and HVR (INPUT)
!     NPREV (I4) = number of vectors |v> allready considered (INPUT)
!      NCOR (I4) = number of corrections vectors |v> to be added (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
!     use f95_precision, only: wp => dp
!     use blas95, only: gemm
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer, intent(in) :: n,nprev,ncor
      real (kind=K_R8), intent(inout), dimension (:) :: vhv
      double precision, intent(in), dimension (:,:) :: hvr,vr
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: j,ij,ij1,ij2,n1,n2,ios,i,k,k2,i2,i3,kend,iend,k3
      double precision :: c11,c21,c22,c12
      double precision :: a11,a21,a22,a12
      double precision :: b11,b21,b22,b12
!
!     ----- allocatable variables -----
!
      double precision, allocatable, dimension (:,:) :: temp
!
!     +++++
!     start
!     +++++
!
      n1 = nprev+1
      n2 = nprev+ncor
      ij1 = (nprev*(nprev+1))/2
      if(nprev.gt.0) then
         allocate(temp(nprev,ncor),stat=ios)
         if(ios.ne.0) then
            write(6,*) 'alloc error in VHV3R',ios
         endif
         call gemm(vr(:,1:nprev),hvr(:,n1:n2),temp,'T','N',ONE,ZERO)
         ij = ij1
         do i=1,ncor
            vhv(ij+1:ij+nprev) = temp(1:nprev,i)
            ij = ij+i+nprev
         enddo
         deallocate(temp)
      endif
      i2 = iand(ncor,1)
      k2 = iand(n,1)
      kend = n-k2
      iend = n2-i2
      k3 = kend+1
      ij = ij1
      do i=n1,iend,2
         ij2 = ij+i
         do j=n1,i-2,2
            c11 = ZERO
            c12 = ZERO
            c21 = ZERO
            c22 = ZERO
            do k=1,kend,2
               a11 = vr(k,j)
               a12 = vr(k,j+1)
               b11 = hvr(k,i)
               b12 = hvr(k,i+1)
               c11 = c11+(a11*b11)
               c21 = c21+(a11*b12)
               c12 = c12+(a12*b11)
               c22 = c22+(a12*b12)
               a21 = vr(k+1,j)
               a22 = vr(k+1,j+1)
               b21 = hvr(k+1,i)
               b22 = hvr(k+1,i+1)
               c11 = c11+(a21*b21)
               c21 = c21+(a21*b22)
               c12 = c12+(a22*b21)
               c22 = c22+(a22*b22)
            enddo
            do k=k3,n
               a11 = vr(k,j)
               a12 = vr(k,j+1)
               b11 = hvr(k,i)
               b12 = hvr(k,i+1)
               c11 = a11*b11+c11
               c12 = a12*b11+c12
               c21 = a11*b12+c21
               c22 = a12*b12+c22
            enddo
            vhv(ij+j) = c11
            vhv(ij+j+1) = c12
            vhv(ij2+j) = c21
            vhv(ij2+j+1) = c22
         enddo
         j = i
         c11 = ZERO
         c21 = ZERO
         c22 = ZERO
         do k=1,kend,2
            a11 = vr(k,j)
            a21 = vr(k+1,j)
            a12 = vr(k,j+1)
            a22 = vr(k+1,j+1)
            b11 = hvr(k,i)
            b21 = hvr(k+1,i)
            b12 = hvr(k,i+1)
            b22 = hvr(k+1,i+1)
            c11 = c11+(a11*b11)+(a21*b21)
            c21 = c21+(a11*b12)+(a21*b22)
            c22 = c22+(a12*b12)+(a22*b22)
         enddo
         do k=k3,n
            a11 = vr(k,j)
            a12 = vr(k,j+1)
            b11 = hvr(k,i)
            b12 = hvr(k,i+1)
            c11 = a11*b11+c11
            c21 = a11*b12+c21
            c22 = a12*b12+c22
         enddo
         vhv(ij+j) = c11
         vhv(ij2+j) = c21
         vhv(ij2+j+1) = c22
         ij = ij2+i+1
      enddo
!
!     ----- clean up I -----
!
      i3 = iend+1
      do i=i3,n2
         ij2 = ij+i
         do j=n1,i-2,2
            c11 = ZERO
            c12 = ZERO
            c21 = ZERO
            c22 = ZERO
            do k=1,kend,2
               a11 = vr(k,j)
               a12 = vr(k,j+1)
               b11 = hvr(k,i)
!              b12 = hvr(k,i+1)
               c11 = c11+(a11*b11)
!              c21 = c21+(a11*b12)
               c12 = c12+(a12*b11)
!              c22 = c22+(a12*b12)
               a21 = vr(k+1,j)
               a22 = vr(k+1,j+1)
               b21 = hvr(k+1,i)
!              b22 = hvr(k+1,i+1)
               c11 = c11+(a21*b21)
!              c21 = c21+(a21*b22)
               c12 = c12+(a22*b21)
!              c22 = c22+(a22*b22)
            enddo
            do k=k3,n
               a11 = vr(k,j)
               a12 = vr(k,j+1)
               b11 = hvr(k,i)
!              b12 = hvr(k,i+1)
               c11 = a11*b11+c11
               c12 = a12*b11+c12
!              c21 = a11*b12+c21
!              c22 = a12*b12+c22
            enddo
            vhv(ij+j) = c11
            vhv(ij+j+1) = c12
!           vhv(ij2+j) = c21
!           vhv(ij2+j+1) = c22
         enddo
!        if(i.eq.n2) cycle
         j = i
         c11 = ZERO
!        c21 = ZERO
!        c22 = ZERO
         do k=1,kend,2
            a11 = vr(k,j)
            a21 = vr(k+1,j)
!           a12 = vr(k,j+1)
!           a22 = vr(k+1,j+1)
            b11 = hvr(k,i)
            b21 = hvr(k+1,i)
!           b12 = hvr(k,i+1)
!           b22 = hvr(k+1,i+1)
            c11 = c11+(a11*b11)+(a21*b21)
!           c21 = c21+(a11*b12)+(a21*b22)
!           c22 = c22+(a12*b12)+(a22*b22)
         enddo
         do k=k3,n
            a11 = vr(k,j)
!           a12 = vr(k,j+1)
            b11 = hvr(k,i)
!           b12 = hvr(k,i+1)
            c11 = a11*b11+c11
!           c21 = a11*b12+c21
!           c22 = a12*b12+c22
         enddo
         vhv(ij+j) = c11
!        vhv(ij2+j) = c21
!        vhv(ij2+j+1) = c22
         ij = ij2+i+1
      enddo
      return
   end subroutine vhv3r
   subroutine vhvsr(vr,hvr,vhv,nc)
!
!***********************************************************************
!     ----- compute the <v|H|v> matrix -----
!=======================================================================
!       VHV (R8) = lower triangle <v|H|v> (OUTPUT)
!        HV (R8) = H*|v> matrix (INPUT)
!         V (R8) = vectors v(n,nc) (INPUT)
!        NC (I4) = number of vectors |v> (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- parameters -----
!  
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: nc
      double precision, intent(out), dimension (:) :: vhv
      double precision, intent(in), dimension (:,:) :: hvr,vr
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: i,j,ij,ii,jj,ij0,iend,iu,m,ib
      double precision, dimension (NB,NB) :: temp
      double precision, dimension (NB*(NB+1)/2) :: tmps
!
!     +++++
!     start
!     +++++
!
      ij0 = 0
      iend = NB*ishft(nc,-SH)
      m = size(vr,1)
      do ii=1,iend,NB
         do jj=1,ii-1,NB
            ij = ij0
            call gemm(vr(:,jj:jj+(NB-1)),hvr(:,ii:ii+(NB-1)),temp,  &
                      'T','N',ONE,ZERO)
            do i=1,NB
               vhv((ij+jj):(ij+jj)+(NB-1)) = temp(1:NB,i)
               ij = ij+(ii-1)+i
            enddo
         enddo
!
!     ----- diagonal bloc -----
!
         call vhv3r(tmps,hvr(:,ii:ii+(NB-1)),vr(:,ii:ii+(NB-1)),m,0,NB)
         ij = ij0
         iu = 1
         do i=1,NB
            do j=1,i
               vhv((ij+jj-1)+j) = tmps(iu)
               iu = iu+1
            enddo
            ij = ij+(ii-1)+i
         enddo
         ij0 = ij
      enddo
!
!      ----- clean up II bloc -----
!
      ib = nc-iend
      if(ib.eq.0) return
      ij0 = (ii*(ii-1))/2
      do jj=1,ii-1,NB
         call gemm(vr(:,jj:jj+(NB-1)),hvr(:,ii:nc),temp(:,1:(nc+1-ii)),  &
                   'T','N',ONE,ZERO)
!        call dgemm('T','N',NB,m,nc+1-ii,ONE,vr(1,jj),m,hvr(1,ii),m,ZERO,temp,NB)
         do i=1,ib
            vhv((ij+jj):(ij+jj)+(NB-1)) = temp(1:NB,i)
            ij = ij+(ii-1)+i
         enddo
         ij = ij0
      enddo
!
!     ----- diagonal bloc -----
!
      call vhv3r(tmps,hvr(:,ii:nc),vr(:,ii:nc),m,0,ib)
      ij = ij0
      iu = 1
      do i=1,ib
         do j=1,i
            vhv((ij+jj-1)+j) = tmps(iu)
            iu = iu+1
         enddo
         ij = ij+(ii-1)+i
      enddo
      return
   end subroutine vhvsr
   subroutine vhvsc2(vr,vi,hvr,hvi,vhv,nc)
!
!***********************************************************************
!     ----- compute the <v|H|v> matrix -----
!=======================================================================
!       VHV (C8) = lower triangle <v|H|v> (OUTPUT)
!        HV (C8) = H*|v> matrix (INPUT)
!         V (C8) = vectors v(n,nc) (INPUT)
!        NC (I4) = number of vectors |v> (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
!     use f95_precision, only: wp => dp
!     use blas95, only: dot,gemm
      implicit none
!
!     ----- parameters -----
!  
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: nc
      real (kind=K_R8), intent(in), dimension (:,:) :: vr,vi,hvr,hvi
      complex (kind=K_C8), intent(out), dimension (:) :: vhv
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: i,j,ij,ii,jj,ij0,iend,iu,m,ib
      real (kind=K_R8), dimension (NB,NB) :: tempr,tempi
      complex (kind=K_C8), dimension (NB*(NB+1)/2) :: tempc
!
!     +++++
!     start
!     +++++
!
      ij0 = 0
      iend = NB*ishft(nc,-SH)
      m = size(vr,1)
      do ii=1,iend,NB
         do jj=1,ii-1,NB
            ij = ij0
            call gemm(vr(:,jj:jj+(NB-1)),hvr(:,ii:ii+(NB-1)),tempr,  &
                      'T','N',ONE,ZERO)
            call gemm(vi(:,jj:jj+(NB-1)),hvi(:,ii:ii+(NB-1)),tempr,  &
                      'T','N',ONE,ONE)
            call gemm(vr(:,jj:jj+(NB-1)),hvi(:,ii:ii+(NB-1)),tempi,  &
                      'T','N',ONE,ZERO)
            call gemm(vi(:,jj:jj+(NB-1)),hvr(:,ii:ii+(NB-1)),tempi,  &
                      'T','N',-ONE,ONE)
            do i=1,NB
               vhv((ij+jj):(ij+jj)+(NB-1)) = cmplx(tempr(1:NB,i),  &
                                                   tempi(1:NB,i),kind=K_C8)
               ij = ij+(ii-1)+i
            enddo
         enddo
!
!     ----- diagonal bloc -----
!
         call vhv3c2(tempc,hvr(:,ii:ii+(NB-1)),hvi(:,ii:ii+(NB-1)),  &
                          vr(:,ii:ii+(NB-1)),vi(:,ii:ii+(NB-1)),m,0,NB)
         ij = ij0
         iu = 1
         do i=1,NB
            do j=1,i
               vhv((ij+jj-1)+j) = tempc(iu)
               iu = iu+1
            enddo
            ij = ij+(ii-1)+i
         enddo
         ij0 = ij
      enddo
!
!      ----- clean up II bloc -----
!
      ib = nc-iend
      if(ib.eq.0) return
      ij0 = (ii*(ii-1))/2
      do jj=1,ii-1,NB
         call gemm(vr(:,jj:jj+(NB-1)),hvr(:,ii:nc),tempr(:,1:ib),'T','N',ONE,ZERO)
         call gemm(vi(:,jj:jj+(NB-1)),hvi(:,ii:nc),tempr(:,1:ib),'T','N',ONE,ONE)
         call gemm(vr(:,jj:jj+(NB-1)),hvi(:,ii:nc),tempi(:,1:ib),'T','N',ONE,ZERO)
         call gemm(vi(:,jj:jj+(NB-1)),hvr(:,ii:nc),tempi(:,1:ib),'T','N',-ONE,ONE)
         do i=1,ib
            vhv((ij+jj):(ij+jj)+(NB-1)) = cmplx(tempr(1:NB,i),  &
                                                tempi(1:NB,i),kind=K_C8)
            ij = ij+(ii-1)+i
         enddo
         ij = ij0
      enddo
!
!     ----- diagonal bloc -----
!
      call vhv3c2(tempc,hvr(:,ii:nc),hvi(:,ii:nc),vr(:,ii:nc),vi(:,ii:nc),  &
                  m,0,ib)
      ij = ij0
      iu = 1
      do i=1,ib
         do j=1,i
            vhv((ij+jj-1)+j) = tempc(iu)
            iu = iu+1
         enddo
         ij = ij+(ii-1)+i
      enddo
!     CALL DBGPRT('VHVSC2: CLEAN UP FINI')
      return
   end subroutine vhvsc2
   subroutine vhv3c2(vhv,hvr,hvi,vr,vi,n,nprev,ncor)
!
!***********************************************************************
!             ----- compute the <v|H|v> matrix -----
!=======================================================================
!       VHV (C8) = lower triangle <v|H|v> (OUTPUT)
!    HV(R,I)(R8) = H*|v> matrix (INPUT)
!     V(R,I)(R8) = vectors v(n,nc) (INPUT)
!         N (I4) = dimension of hvr, hvi, vr, vi (INPUT)
!     NPREV (I4) = dimension of the original space (INPUT)
!      NCOR (I4) = number of new vectors |v> (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer, intent(in) :: n,nprev,ncor
      complex (kind=K_C8), intent(inout), dimension (:) :: vhv
      double precision, intent(in), dimension (:,:) :: hvr,vr,hvi,vi
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: j,ij,ij1,ij2,n1,n2,ios,i,k,k2,i2,kend,k3,iend,i3
      double precision :: dumr11,dumi11,dumr12,dumi12,dumr21,dumi21, &
                          dumr22,dumi22,  &
                          ar11,ai11,ar12,ai12,ar21,ai21,ar22,ai22,  &
                          br11,bi11,br12,bi12,br21,bi21,br22,bi22
!
!     ----- allocatable variables -----
!
      double precision, allocatable, dimension (:,:) :: tempr,tempi
!
!     +++++
!     start
!     +++++
!
      n1 = nprev+1
      n2 = nprev+ncor
      ij1 = (nprev*(nprev+1))/2
!     CALL DBGPRT('   ---VHV3C2 ENTREE: N1 N2 IJ1=',(/ N1,N2,IJ1 /),3)
      if(nprev.gt.0) then
         allocate(tempr(nprev,ncor),tempi(nprev,ncor),stat=ios)
         if(ios.ne.0) then
            write(6,*) 'alloc error in VHV4',ios
         endif
         call gemm(vr(:,1:nprev),hvr(:,n1:n2),tempr,'T','N',ONE,ZERO)
         call gemm(vi(:,1:nprev),hvi(:,n1:n2),tempr,'T','N',ONE,ONE)
         call gemm(vi(:,1:nprev),hvr(:,n1:n2),tempi,'T','N',-ONE,ZERO)
         call gemm(vr(:,1:nprev),hvi(:,n1:n2),tempi,'T','N',ONE,ONE)
         ij = ij1
         do i=1,ncor
            vhv(ij+1:ij+nprev) = dcmplx(tempr(1:nprev,i),tempi(1:nprev,i))
            ij = ij+i+nprev
         enddo
         deallocate(tempi,tempr)
      endif
      i2 = iand(ncor,1)
      k2 = iand(n,1)
      kend = n-k2
      iend = n2-i2
      k3 = kend+1
      ij = ij1
!     CALL DBGPRT('   ---VHV3C2 APRES NPREV:I2 K2 K3=',(/ I2,K2,K3 /),3)
      do i=n1,iend,2
         ij2 = ij+i
         do j=n1,i-2,2
            dumr11 = ZERO
            dumi11 = ZERO
            dumr21 = ZERO
            dumi21 = ZERO
            dumr12 = ZERO
            dumi12 = ZERO
            dumr22 = ZERO
            dumi22 = ZERO
            do k=1,kend,2
               ar11 = vr(k,j)
               ar21 = vr(k+1,j)
               br11 = hvr(k,i)
               br21 = hvr(k+1,i)
               bi11 = hvi(k,i)
               bi21 = hvi(k+1,i)
               dumr11 = dumr11+(ar11*br11)+  &
                               (ar21*br21)
               dumi11 = dumi11+(ar11*bi11)+  &
                               (ar21*bi21)
               ai11 = vi(k,j)
               ai21 = vi(k+1,j)
               dumr11 = dumr11+(ai11*bi11)+  &
                               (ai21*bi21)
               dumi11 = dumi11-(ai11*br11)-  &
                               (ai21*br21)
               br12 = hvr(k,i+1)
               br22 = hvr(k+1,i+1)
               bi12 = hvi(k,i+1)
               bi22 = hvi(k+1,i+1)
               dumr21 = dumr21+(ar11*br12)+  &
                               (ar21*br22)
               dumi21 = dumi21+(ar11*bi12)+  &
                               (ar21*bi22)
               dumr21 = dumr21+(ai11*bi12)+  &
                               (ai21*bi22)
               dumi21 = dumi21-(ai11*br12)-  &
                               (ai21*br22)
               ar12 = vr(k,j+1)
               ar22 = vr(k+1,j+1)
               dumr12 = dumr12+(ar12*br11)+  &
                               (ar22*br21)
               dumi12 = dumi12+(ar12*bi11)+  &
                               (ar22*bi21)
               ai12 = vi(k,j+1)
               ai22 = vi(k+1,j+1)
               dumr12 = dumr12+(ai12*bi11)+  &
                               (ai22*bi21)
               dumi12 = dumi12-(ai12*br11)-  &
                               (ai22*br21)
               dumr22 = dumr22+(ar12*br12)+  &
                               (ar22*br22)
               dumi22 = dumi22+(ar12*bi12)+  &
                               (ar22*bi22)
               dumr22 = dumr22+(ai12*bi12)+  &
                               (ai22*bi22)
               dumi22 = dumi22-(ai12*br12)-  &
                               (ai22*br22)
            enddo
            do k=k3,n
               ar11 = vr(k,j)
!              ar21 = vr(k+1,j)
               br11 = hvr(k,i)
               bi11 = hvi(k,i)
!              br21 = hvr(k+1,i)
!              bi21 = hvi(k+1,i)
               dumr11 = dumr11+(ar11*br11)
               dumi11 = dumi11+(ar11*bi11)
               ai11 = vi(k,j)
               br12 = hvr(k,i+1)
               bi12 = hvi(k,i+1)
               dumr11 = dumr11+(ai11*bi11)
               dumi11 = dumi11-(ai11*br11)
               dumr21 = dumr21+(ar11*br12)
               dumi21 = dumi21+(ar11*bi12)
               dumr21 = dumr21+(ai11*bi12)
               dumi21 = dumi21-(ai11*br12)
               ar12 = vr(k,j+1)
               ai12 = vi(k,j+1)
               dumr12 = dumr12+(ar12*br11)
               dumi12 = dumi12+(ar12*bi11)
               dumr12 = dumr12+(ai12*bi11)
               dumi12 = dumi12-(ai12*br11)
               dumr22 = dumr22+(ar12*br12)
               dumi22 = dumi22+(ar12*bi12)
               dumr22 = dumr22+(ai12*bi12)
               dumi22 = dumi22-(ai12*br12)
            enddo
            vhv(ij+j) = dcmplx(dumr11,dumi11)
            vhv(ij+j+1) = dcmplx(dumr12,dumi12)
            vhv(ij2+j) = dcmplx(dumr21,dumi21)
            vhv(ij2+j+1) = dcmplx(dumr22,dumi22)
         enddo
         j = i
         dumr11 = ZERO
         dumi11 = ZERO
         dumr21 = ZERO
         dumi21 = ZERO
         dumr22 = ZERO
         dumi22 = ZERO
         do k=1,kend,2
            ar11 = vr(k,j)
            ar21 = vr(k+1,j)
            br11 = hvr(k,i)
            br21 = hvr(k+1,i)
            bi11 = hvi(k,i)
            bi21 = hvi(k+1,i)
            dumr11 = dumr11+(ar11*br11)+  &
                            (ar21*br21)
            dumi11 = dumi11+(ar11*bi11)+  &
                            (ar21*bi21)
            ai11 = vi(k,j)
            ai21 = vi(k+1,j)
            dumr11 = dumr11+(ai11*bi11)+  &
                            (ai21*bi21)
            dumi11 = dumi11-(ai11*br11)-  &
                            (ai21*br21)
            br12 = hvr(k,i+1)
            br22 = hvr(k+1,i+1)
            bi12 = hvi(k,i+1)
            bi22 = hvi(k+1,i+1)
            dumr21 = dumr21+(ar11*br12)+  &
                            (ar21*br22)
            dumi21 = dumi21+(ar11*bi12)+  &
                            (ar21*bi22)
            dumr21 = dumr21+(ai11*bi12)+  &
                            (ai21*bi22)
            dumi21 = dumi21-(ai11*br12)-  &
                            (ai21*br22)
            ar12 = vr(k,j+1)
            ar22 = vr(k+1,j+1)
            dumr22 = dumr22+(ar12*br12)+  &
                            (ar22*br22)
            dumi22 = dumi22+(ar12*bi12)+  &
                            (ar22*bi22)
            ai12 = vi(k,j+1)
            ai22 = vi(k+1,j+1)
            dumr22 =dumr22+(ai12*bi12)+  &
                            (ai22*bi22)
            dumi22 = dumi22-(ai12*br12)-  &
                            (ai22*br22)
         enddo
         do k=k3,n
            ar11 = vr(k,j)
            br11 = hvr(k,i)
            bi11 = hvi(k,i)
            dumr11 = dumr11+(ar11*br11)
            dumi11 = dumi11+(ar11*bi11)
            ai11 = vi(k,j)
!           ai21 = vi(k+1,j)
            dumr11 = dumr11+(ai11*bi11)
            dumi11 = dumi11-(ai11*br11)
            br12 = hvr(k,i+1)
            bi12 = hvi(k,i+1)
            dumr21 = dumr21+(ar11*br12)
            dumi21 = dumi21+(ar11*bi12)
            dumr21 = dumr21+(ai11*bi12)
            dumi21 = dumi21-(ai11*br12)
            ar12 = vr(k,j+1)
            ai12 = vi(k,j+1)
            dumr22 = dumr22+(ar12*br12)
            dumi22 = dumi22+(ar12*bi12)
            dumr22 = dumr22+(ai12*bi12)
            dumi22 = dumi22-(ai12*br12)
         enddo
         vhv(ij+j) = dcmplx(dumr11,dumi11)
         vhv(ij2+j) = dcmplx(dumr21,dumi21)
         vhv(ij2+j+1) = dcmplx(dumr22,dumi22)
         ij = ij2+i+1
      enddo
!
!     ----- clean up I -----
!
!     CALL DBGPRT('   ---VHV3C2 CLEAN UP')
      i3 = iend+1
      do i=i3,n2
         ij2 = ij+i
         do j=n1,i-2,2
            dumr11 = ZERO
            dumi11 = ZERO
!           dumr21 = ZERO
!           dumi21 = ZERO
            dumr12 = ZERO
            dumi12 = ZERO
!           dumr22 = ZERO
!           dumi22 = ZERO
            do k=1,kend,2
               ar11 = vr(k,j)
               ar21 = vr(k+1,j)
               br11 = hvr(k,i)
               bi11 = hvi(k,i)
               br21 = hvr(k+1,i)
               bi21 = hvi(k+1,i)
               dumr11 = dumr11+(ar11*br11)+  &
                               (ar21*br21)
               dumi11 = dumi11+(ar11*bi11)+  &
                               (ar21*bi21)
               ai11 = vi(k,j)
               ai21 = vi(k+1,j)
               dumr11 = dumr11+(ai11*bi11)+  &
                               (ai21*bi21)
               dumi11 = dumi11-(ai11*br11)-  &
                               (ai21*br21)
!              dumr21 = dumr21+(ar11*br12)+  &
!                              (ar21*br22)
!              dumi21 = dumi21+(ar11*bi12)+  &
!                              (ar21*bi22)
!              dumr21 = dumr21+(ai11*bi12)+  &
!                              (ai21*bi22)
!              dumi21 = dumi21-(ai11*br12)-  &
!                              (ai21*br22)
               ar12 = vr(k,j+1)
!              br12 = hvr(k,i+1)
               ar22 = vr(k+1,j+1)
!              br22 = hvr(k+1,i+1)
               dumr12 = dumr12+(ar12*br11)+  &
                               (ar22*br21)
               dumi12 = dumi12+(ar12*bi11)+  &
                               (ar22*bi21)
!              dumr22 = dumr22+(ar12*br12)+  &
!                              (ar22*br22)
               ai12 = vi(k,j+1)
               ai22 = vi(k+1,j+1)
!              bi12 = hvi(k,i+1)
!              bi22 = hvi(k+1,i+1)
               dumr12 = dumr12+(ai12*bi11)+  &
                               (ai22*bi21)
               dumi12 = dumi12-(ai12*br11)-  &
                               (ai22*br21)
!              dumi22 = dumi22+(ar12*bi12)+  &
!                              (ar22*bi22)
!              dumr22 = dumr22+(ai12*bi12)+  &
!                              (ai22*bi22)
!              dumi22 = dumi22-(ai12*br12)-  &
!                              (ai22*br22)
            enddo
            do k=k3,n
               ar11 = vr(k,j)
               ai11 = vi(k,j)
               br11 = hvr(k,i)
               bi11 = hvi(k,i)
               dumr11 = dumr11+(ar11*br11)
               dumi11 = dumi11+(ar11*bi11)
               dumr11 = dumr11+(ai11*bi11)
               dumi11 = dumi11-(ai11*br11)
!              br12 = hvr(k,i+1)
!              bi12 = hvi(k,i+1)
!              dumr21 = dumr21+(ar11*br12)
!              dumi21 = dumi21+(ar11*bi12)
!              dumr21 = dumr21+(ai11*bi12)
!              dumi21 = dumi21-(ai11*br12)
               ar12 = vr(k,j+1)
               ai12 = vi(k,j+1)
               dumr12 = dumr12+(ar12*br11)
               dumi12 = dumi12+(ar12*bi11)
               dumr12 = dumr12+(ai12*bi11)
               dumi12 = dumi12-(ai12*br11)
!              dumr22 = dumr22+(ar12*br12)
!              dumi22 = dumi22+(ar12*bi12)
!              dumr22 = dumr22+(ai12*bi12)
!              dumi22 = dumi22-(ai12*br12)
            enddo
            vhv(ij+j) = dcmplx(dumr11,dumi11)
            vhv(ij+j+1) = dcmplx(dumr12,dumi12)
!           vhv(ij2+j) = dcmplx(dumr21,dumi21)
!           vhv(ij2+j+1) = dcmplx(dumr22,dumi22)
         enddo
!        if(i.eq.n2) cycle
         j = i
         dumr11 = CZERO
!        dumr21 = CZERO
!        dumr22 = CZERO
         dumi11 = CZERO
!        dumi21 = CZERO
!        dumi22 = CZERO
         do k=1,kend,2
            ar11 = vr(k,j)
            ar21 = vr(k+1,j)
!           ar12 = vr(k,j+1)
!           ar22 = vr(k+1,j+1)
            br11 = hvr(k,i)
            br21 = hvr(k+1,i)
!           br12 = hvr(k,i+1)
!           br22 = hvr(k+1,i+1)
            dumr11 = dumr11+(ar11*br11)+(ar21*br21)
!           dumr21 = dumr21+(ar11*br12)+(ar21*br22)
!           dumr22 = dumr22+(ar12*br12)+(ar22*br22)
            bi11 = hvi(k,i)
            bi21 = hvi(k+1,i)
!           bi12 = hvi(k,i+1)
!           bi22 = hvi(k+1,i+1)
            dumi11 = dumi11+(ar11*bi11)+(ar21*bi21)
!           dumi21 = dumi21+(ar11*bi12)+(ar21*bi22)
!           dumi22 = dumi22+(ar12*bi12)+(ar22*bi22)
            ai11 = vi(k,j)
            ai21 = vi(k+1,j)
!           ai12 = vi(k,j+1)
!           ai22 = vi(k+1,j+1)
            dumr11 = dumr11+(ai11*bi11)+(ai21*bi21)
!           dumr21 = dumr21+(ai11*bi12)+(ai21*bi22)
            dumi11 = dumi11-(ai11*br11)-(ai21*br21)
!           dumi21 = dumi21-(ai11*br12)-(ai21*br22)
!           dumr22 = dumr22+(ai12*bi12)+(ai22*bi22)
!           dumi22 = dumi22-(ai12*br12)-(ai22*br22)
         enddo
         do k=k3,n
            ar11 = vr(k,j)
            br11 = hvr(k,i)
            bi11 = hvi(k,i)
            dumr11 = ar11*br11+dumr11
            dumi11 = ar11*bi11+dumi11
            ai11 = vi(k,j)
            dumr11 = ai11*bi11+dumr11
            dumi11 = dumi11-ai11*br11
         enddo
         vhv(ij+j) = dcmplx(dumr11,dumi11)
!        vhv(ij2+j) = dcmplx(dumr21,dumi21)
!        vhv(ij2+j+1) = dcmplx(dumr22,dumi22)
         ij = ij2+i+1
      enddo
!     CALL DBGPRT('   ---VHV3C2 CLEAN UP FINI')
      return
   end subroutine vhv3c2
end module mmulsy
