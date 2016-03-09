module matpri
!
!***********************************************************************
!          ----- generic module for printing square matrices -----
!          ----- and packed triangular matrices              -----
!=======================================================================
!     contains:
!       subroutine carmat(a,n,m)  generic subroutine
!       subroutine matpr_i4(a,n,m)  subroutine for I*4 matrices
!       subroutine matpr_i2(a,n,m)  subroutine for I*2 matrices
!       subroutine matpr_i8(a,n,m)  subroutine for I*8 matrices
!       subroutine matpr_r4(a,n,m)  subroutine for R*4 matrices
!       subroutine matpr_r8(a,n,m)  subroutine for R*8 matrices
!       subroutine matpr_rv8(a,n,m)  subroutine for R*8 matrices
!       subroutine matpr_c4(a,n,m)  subroutine for C*8 matrices
!       subroutine matpr_c8(a,n,m)  subroutine for C*16 matrices
!       subroutine matpr_ri8(a,b,n,m)  subroutine for C*16 matrices
!                                        a+I*b
!       subroutine matpr_l4(a,n,m)  subroutine for L*4 matrices
!       subroutine carmat(e,a,n,m)  generic subroutine
!       subroutine trimat(a,n)        generic routine
!       subroutine matpre_r8(e,a,n,m)  subroutine for R*8 matrices
!                                          with R*8 eigenvalues
!       subroutine matpre_c8(e,a,n,m)  subroutine for C*8 matrices
!                                          with R*8 eigenvalues
!       subroutine tripr_i4(a,n)  subroutine for symmetric I*4 matrices
!       subroutine tripr_i2(a,n)  subroutine for symmetric I*2 matrices
!       subroutine tripr_i8(a,n)  subroutine for symmetric I*8 matrices
!       subroutine tripr_r4(a,n)  subroutine for symmetric R*4 matrices
!       subroutine tripr_r8(a,n)  subroutine for symmetric R*8 matrices
!       subroutine tripr_c4(a,n)  subroutine for symmetric C*4 matrices
!       subroutine tripr_c8(a,n)  subroutine for symmetric C*8 matrices
!
!       makfo4 is a private subroutine to set the format for integer data
!       makfo8 is a private subroutine to set the format for integer data
!-----------------------------------------------------------------------
!   usage:
!   ------
!   call carmat(a,n,m)  to print a matrix (dimension a(:,:) or a(:))
!   call carmat(reshape(a,(/1,n/),a,1,1,n)  to print a vector linewise
!   call carmat(reshape(a,(/n,1/),a,n,n,1)  to print a vector columnwise
!
!      A(??) = matrix to be printed (INPUT)
!    LDA(I4) = leading dimension of A (i.e. dimension a(:,:) (INPUT)
!      M(I4) = number of rows to be printed  (INPUT)
!      N(I4) = number of columns to be printed  (INPUT)
!
!   call carmat(a,n,m)  to print a matrix (dimension a(:,:))
!                           with eigenvalues
!      E(R8) = eigenvalues (INPUT)
!      A(R8) = matrix to be printed (INPUT)
!    LDA(I4) = leading dimension of A (i.e. dimension a(:,:) (INPUT)
!      M(I4) = number of rows to be printed  (INPUT)
!      N(I4) = number of columns to be printed  (INPUT)
!
!   call carmat(a,n)  to print a symmetrc (hermitian) square matrix
!                     only the lower triangle is stored: a(ij) = a(i,j)
!                     with ij = i*(i-1)/2+j , i.GE. j
!      A(R8) = matrix to be printed (INPUT)
!      N(I4) = number of rows/columns to be printed  (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
use kindvar
interface carmat
   module procedure matpr_r8, matpr_r4, matpr_i4, matpr_i2, matpr_i8,   &
                    matpr_c4, matpr_c8, tripr_i4, tripr_i2, tripr_i8,   &
                    tripr_r4, tripr_r8, tripr_c4, tripr_c8, matpr_l4,   &
                    matpre_r8, matpre_c8, matpr_ri8, matpr_rv8
end interface carmat
private :: makfo4,makfo8
contains
subroutine matpr_r8(a,mm,nn)
!
!***********************************************************************
!     ----- impression d'une matrice R*8 ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   double precision, intent(in), dimension (:,:) ::  a
!
!     ----- local variables -----
!
   integer :: i,j,k,m
!
 9999 format(//,6x,12(4x,i3,3x))
 9998 format(1x,i3,2x,12(f10.6))
 9997 format(/)
!
!     +++++
!     start
!     +++++
!
   do m=1,nn,12
      k = min(m+11,nn)
      write(6,9999) (j,j=m,k)
      write(6,9997)
      do i=1,mm
         write(6,9998) i,(a(i,j),j=m,k)
      enddo
   enddo
   return
end subroutine matpr_r8
subroutine matpr_rv8(a,mm,nn)
!
!***********************************************************************
!    ----- impression d'une matrice R*8 donnee sous forme vecteur ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   double precision, intent(in), dimension (:) ::  a
!
!     ----- local variables -----
!
   integer :: i,j,k,m,lda
!
 9999 format(//,6x,12(4x,i3,3x))
 9998 format(1x,i3,2x,12(f10.6))
 9997 format(/)
!
!     +++++
!     start
!     +++++
!
   lda = size(a,1)
   do m=1,nn,12
      k = min(m+11,nn)
      write(6,9999) (j,j=m,k)
      write(6,9997)
      do i=1,mm
         write(6,9998) i,(a(i+(j-1)*lda),j=m,k)
      enddo
   enddo
   return
end subroutine matpr_rv8
subroutine matpr_r4(a,mm,nn)
!
!***********************************************************************
!     ----- impression d'une matrice R*4 ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   real, intent(in), dimension (:,:) ::  a
!
!     ----- local variables -----
!
   integer :: i,j,k,m
!
 9999 format(//,6x,12(4x,i3,3x))
 9998 format(1x,i3,2x,12(f10.6))
 9997 format(/)
!
!     +++++
!     start
!     +++++
!
   do m=1,nn,12
      k = min(m+11,nn)
      write(6,9999) (j,j=m,k)
      write(6,9997)
      do i=1,mm
         write(6,9998) i,(a(i,j),j=m,k)
      enddo
   enddo
   return
end subroutine matpr_r4
subroutine matpr_i4(a,mm,nn)
!
!***********************************************************************
!     ----- impression d'une matrice R*4 ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   integer (kind=K_I4), intent(in), dimension (:,:) ::  a
!
!     ----- local variables -----
!
   character (len=95) :: f1,f2,f3
   character (len=15) :: string
   integer :: i,j,k,m,ncol
   integer (kind=K_I4) :: maxt
!
!     +++++
!     start
!     +++++
!
   maxt = maxval(abs(a(1:mm,1:nn)))
   maxt = max(maxt,1)
   call makfo4(f1,f2,f3,ncol,string,max(mm,nn),maxt)
   do m=1,nn,ncol
      k = min(m+ncol-1,nn)
      write(6,fmt=f1) (j,j=m,k)
      write(6,fmt=f2) (trim(string),j=m,k)
      do i=1,mm
         write(6,fmt=f3) i,(a(i,j),j=m,k)
      enddo
   enddo
   return
end subroutine matpr_i4
subroutine matpr_i2(a,mm,nn)
!
!***********************************************************************
!     ----- impression d'une matrice I*2 ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: mm,nn
   integer (kind=K_I2), intent(in), dimension (:,:) :: a
!
!     ----- local variables -----
!
   character (len=95) :: f1,f2,f3
   character (len=15) :: string
   integer :: i,j,k,m,ncol
   integer (kind=K_I4) :: maxt
!
!     +++++
!     start
!     +++++
!
   maxt = maxval(abs(a(1:mm,1:nn)))
   call makfo4(f1,f2,f3,ncol,string,max(mm,nn),maxt)
   do m=1,nn,ncol
      k = min(m+ncol-1,nn)
      write(6,fmt=f1) (j,j=m,k)
      write(6,fmt=f2) (trim(string),j=m,k)
      do i=1,mm
         write(6,fmt=f3) i,(a(i,j),j=m,k)
      enddo
   enddo
   return
end subroutine matpr_i2
subroutine matpr_i8(a,mm,nn)
!
!***********************************************************************
!     ----- impression d'une matrice I*2 ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: mm,nn
   integer (kind=K_I8), intent(in), dimension (:,:) :: a
!
!     ----- local variables -----
!
   character (len=95) :: f1,f2,f3
   character (len=15) :: string
   integer :: i,j,k,m,ncol
   integer (kind=K_I8) :: maxt
!
!     +++++
!     start
!     +++++
!
   maxt = maxval(abs(a(1:mm,1:nn)))
   call makfo8(f1,f2,f3,ncol,string,max(mm,nn),maxt)
   do m=1,nn,ncol
      k = min(m+ncol-1,nn)
      write(6,fmt=f1) (j,j=m,k)
      write(6,fmt=f2) (trim(string),j=m,k)
      do i=1,mm
         write(6,fmt=f3) i,(a(i,j),j=m,k)
      enddo
   enddo
   return
end subroutine matpr_i8
subroutine matpr_c4(a,mm,nn)
!
!**********************************************************************
!     ----- Print a double complex table -----
!**********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   complex (kind=K_C4), intent(in), dimension (:,:) :: a
!
!     ----- local variables -----
!
   integer :: i,j,k,m
!
 9999 format(i3,4('  (',f12.6,',',f12.6,')',:))
 9998 format(3x,4('   ',12x,i3,9x,' ',:))
 9997 format(/)
!
!     +++++
!     start
!     +++++
!
      do m=1,nn,4
         k = min(m+3,nn)
         write(6,9998) (j,j=m,k)
         write(6,9997)
         do i=1,mm
            write(6,9999) i,(a(i,j),j=m,k)
         enddo
      enddo
      return
end subroutine matpr_c4
subroutine matpr_c8(a,mm,nn)
!
!**********************************************************************
!     ----- Print a double complex table -----
!**********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   complex (kind=K_C8), intent(in), dimension (:,:) :: a
!
!     ----- local variables -----
!
   integer :: i,j,k,m
!
 9999 format(i5,4(' (',f12.6,',',f12.6,')',:))
 9998 format(5x,4('  ',12x,i3,9x,' ',:))
 9997 format(/)
!
!     +++++
!     start
!     +++++
!
      do m=1,nn,4
         k = min(m+3,nn)
         write(6,9998) (j,j=m,k)
         write(6,9997)
         do i=1,mm
            write(6,9999) i,(a(i,j),j=m,k)
         enddo
      enddo
      return
end subroutine matpr_c8
subroutine matpr_ri8(a,b,mm,nn)
!
!**********************************************************************
!     ----- Print a double complex table -----
!**********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   double precision, intent(in), dimension (:,:) :: a,b
!
!     ----- local variables -----
!
   integer :: i,j,k,m
!
 9999 format(i3,4('  (',f12.6,',',f12.6,')',:))
 9998 format(3x,4('   ',12x,i3,9x,' ',:))
 9997 format(/)
!
!     +++++
!     start
!     +++++
!
      do m=1,nn,4
         k = min(m+3,nn)
         write(6,9998) (j,j=m,k)
         write(6,9997)
         do i=1,mm
            write(6,9999) i,(a(i,j),b(i,j),j=m,k)
         enddo
      enddo
      return
end subroutine matpr_ri8
subroutine matpr_l4(a,mm,nn)
!
!**********************************************************************
!     ----- Print a double complex table -----
!**********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   logical, intent(in), dimension (:,:) :: a
!
!     ----- local variables -----
!
   integer :: i,j,k,m
!
 9999 format(1x,i4,1x,30(2x,l1,1x))
 9998 format(//,6x,30(i4))
 9997 format(/)
!
!     +++++
!     start
!     +++++
!
      do m=1,nn,30
         k = min(m+29,nn)
         write(6,9998) (j,j=m,k)
         write(6,9997)
         do i=1,mm
            write(6,9999) i,(a(i,j),j=m,k)
         enddo
      enddo
      return
end subroutine matpr_l4
subroutine tripr_i4(a,n)
!
!***********************************************************************
!     ----- Print an integer*4 triangular matrix -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: n
   integer (kind=K_I4), intent(in), dimension (:) ::  a
!
!     ----- local variables -----
!
   character (len=95) :: f1,f2,f3
   character (len=15) :: string
   integer :: m,k,i1,j,i,klim,ii,ncol
   integer (kind=K_I4) :: maxt
!
!     +++++
!     start
!     +++++
!
   maxt = maxval(a(1:n*(n+1)/2))
   call makfo4(f1,f2,f3,ncol,string,n,maxt)
   do m=1,n,ncol
      k = min(m+ncol-1,n)
      i1 = min(m,n)
      write(6,f1) (j,j=m,k)
      write(6,f2) (string,j=m,k)
      do i=i1,n
         klim = min(k,i)
         ii = (i*(i-1))/2
         write(6,f3) i,(a(ii+j),j=m,klim)
      enddo
   enddo
   return
end subroutine tripr_i4
subroutine tripr_i2(a,n)
!
!***********************************************************************
!     ----- Print an integer*2 triangular matrix -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: n
   integer (kind=K_I2), intent(in), dimension (:) ::  a
!
!     ----- local variables -----
!
   character (len=95) :: f1,f2,f3
   character (len=15) :: string
   integer :: m,k,i1,j,i,klim,ii,ncol
   integer (kind=K_I4) :: maxt
!
!     +++++
!     start
!     +++++
!
   maxt = maxval(a(1:n*(n+1)/2))
   call makfo4(f1,f2,f3,ncol,string,n,maxt)
   do m=1,n,ncol
      k = min(m+ncol-1,n)
      i1 = min(m,n)
      write(6,f1) (j,j=m,k)
      write(6,f2) (string,j=m,k)
      do i=i1,n
         klim = min(k,i)
         ii = (i*(i-1))/2
         write(6,f3) i,(a(ii+j),j=m,klim)
      enddo
   enddo
   return
end subroutine tripr_i2
subroutine tripr_i8(a,n)
!
!***********************************************************************
!     ----- Print an integer*4 triangular matrix -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: n
   integer (kind=K_I8), intent(in), dimension (:) ::  a
!
!     ----- local variables -----
!
   character (len=95) :: f1,f2,f3
   character (len=15) :: string
   integer :: m,k,i1,j,i,klim,ii,ncol
   integer (kind=K_I8) :: maxt
!
!     +++++
!     start
!     +++++
!
   maxt = maxval(a(1:n*(n+1)/2))
   call makfo8(f1,f2,f3,ncol,string,n,maxt)
   do m=1,n,ncol
      k = min(m+ncol-1,n)
      i1 = min(m,n)
      write(6,f1) (j,j=m,k)
      write(6,f2) (string,j=m,k)
      do i=i1,n
         klim = min(k,i)
         ii = (i*(i-1))/2
         write(6,f3) i,(a(ii+j),j=m,klim)
      enddo
   enddo
   return
end subroutine tripr_i8
subroutine tripr_r4(a,n)
!
!***********************************************************************
!     ----- Print a complex triangular matrix -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: n
   real, intent(in), dimension (:) ::  a
!
!     ----- local variables -----
!
   integer :: m,k,i1,j,i,klim,ii
!
!     ----- formats -----
!
   20 format(//,4x,10(6x,i3,3x))
   30 format(i4,10(f12.6))
   36 format(/)
!
!     +++++
!     start
!     +++++
!
   do m=1,n,10
      k = min(m+9,n)
      i1 = min(m,n)
      write(6,20) (j,j=m,k)
      write(6,36)
      do i=i1,n
         klim = min(k,i)
         ii = (i*(i-1))/2
         write(6,30) i,(a(ii+j),j=m,klim)
      enddo
   enddo
   return
end subroutine tripr_r4
subroutine tripr_r8(a,n)
!
!***********************************************************************
!     ----- Print a double precision triangular matrix -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: n
   double precision, intent(in), dimension (:) ::  a
!
!     ----- local variables -----
!
   integer :: m,k,i1,j,i,klim,ii
!
!     ----- formats -----
!
   20 format(//,4x,10(6x,i3,3x))
   30 format(i4,10(f12.6))
   36 format(/)
!
!     +++++
!     start
!     +++++
!
   do m=1,n,10
      k = min(m+9,n)
      i1 = min(m,n)
      write(6,20) (j,j=m,k)
      write(6,36)
      do i=i1,n
         klim = min(k,i)
         ii = (i*(i-1))/2
         write(6,30) i,(a(ii+j),j=m,klim)
      enddo
   enddo
   return
end subroutine tripr_r8
subroutine tripr_c4(a,n)
!
!***********************************************************************
!     ----- Print a complex triangular matrix -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: n
   complex, intent(in), dimension (:) ::  a
!
!     ----- local variables -----
!
   integer :: m,k,i1,j,i,klim,ii
!
!     ----- formats -----
!
   20 format(//,9x,4(12x,i3,12x))
   30 format(i3,4('  (',f12.6,',',f12.6,')',:))
   36 format(/)
!
!     +++++
!     start
!     +++++
!
   do m=1,n,4
      k = min(m+3,n)
      i1 = min(m,n)
      write(6,20) (j,j=m,k)
      write(6,36)
      do i=i1,n
         klim = min(k,i)
         ii = (i*(i-1))/2
         write(6,30) i,(a(ii+j),j=m,klim)
      enddo
   enddo
   return
end subroutine tripr_c4
subroutine tripr_c8(a,n)
!
!***********************************************************************
!     ----- Print a complex triangular matrix -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: n
   complex (kind=K_C8), intent(in), dimension (:) ::  a
!
!     ----- local variables -----
!
   integer :: m,k,i1,j,i,klim,ii
!
!     ----- formats -----
!
   20 format(//,9x,4(12x,i3,12x))
   30 format(i3,4('  (',f12.6,',',f12.6,')',:))
   36 format(/)
!
!     +++++
!     start
!     +++++
!
   do m=1,n,4
      k = min(m+3,n)
      i1 = min(m,n)
      write(6,20) (j,j=m,k)
      write(6,36)
      do i=i1,n
         klim = min(k,i)
         ii = (i*(i-1))/2
         write(6,30) i,(a(ii+j),j=m,klim)
      enddo
   enddo
   return
end subroutine tripr_c8
subroutine matpre_r8(e,a,mm,nn)
!
!***********************************************************************
!     ----- impression d'une matrice et valeurs propres ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   double precision, intent(in),dimension (:) ::  e
   double precision, intent(in),dimension (:,:) ::  a
!
!     ----- local variables -----
!
   integer :: i,j,k,m
!
!     +++++
!     start
!     +++++
!
 9999 format(//,6x,12(4x,i3,3x))
 9998 format(1x,i3,2x,12(f10.6))
 9997 format(/)
 9996 format(1x,3x,2x,12(f10.6))
!
   do m=1,nn,12
      k=m+11
      if(k.gt.nn) k=nn
      write(6,9999) (j,j=m,k)
      write(6,9997)
      write(6,9996) (e(j),j=m,k)
      write(6,9997)
      do i=1,mm
         write(6,9998) i,(a(i,j),j=m,k)
      enddo
   enddo
   return
end subroutine matpre_r8
subroutine matpre_c8(e,a,mm,nn)
!
!***********************************************************************
!     ----- impression d'une matrice et valeurs propres ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) ::  mm,nn
   double precision, intent(in),dimension (:) ::  e
   double complex, intent(in),dimension (:,:) ::  a
!
!     ----- local variables -----
!
   integer :: i,j,k,m
!
!     ----- formats -----
!
 9999 format(6x,4(14x,i3,11x))
 9998 format(1x,i5,4(' (',f12.6,',',f12.6,')':))
 9997 format(/)
 9996 format(1x,3x,2x,4(10x,f10.6,8x))
!
!     +++++
!     start
!     +++++
!
   do m=1,nn,4
      k=m+3
      if(k.gt.nn) k=nn
      write(6,9997)
      write(6,9999) (j,j=m,k)
      write(6,9997)
      write(6,9996) (e(j),j=m,k)
      write(6,9997)
      do i=1,mm
         write(6,9998) i,(a(i,j),j=m,k)
      enddo
   enddo
   return
end subroutine matpre_c8
subroutine makfo4(f1,f2,f3,ncol,string,mm,maxa)
!
!***********************************************************************
!     ----- impression d'une matrice R*4 ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   character (len=*), intent(out) :: f1,f2,f3,string
   integer, intent(in) :: mm
   integer (kind=K_I4), intent(in) :: maxa
   integer, intent(out) :: ncol
!
!     ----- local variables -----
!
   character (len=30) :: s1,s2,s3
   integer :: i,maxl,maxt,ic,nx1,nx2
!
!     +++++
!     start
!     +++++
!
   maxl = 1+int(log10(dble(mm)))
   maxt = 1+int(log10(dble(maxa)))
   ncol = (130-maxl)/(max(maxt,maxl)+2)
   nx1 = 1+maxl/10
   nx2 = 2+maxt/10
!
!     ----- format for column headers -----
!
   write(s1,*) "(1x,i",nx2,",1x:)"
   ic = 0
   do i=1,len_trim(s1)
      if(s1(i:i).ne.' ') then
         ic = ic+1
         s1(ic:ic) = s1(i:i)
      endif
   enddo
   s1(ic+1:) = ' '
   write(f1,*) "(//,",2+nx1,"x,'u|',",ncol,trim(s1),")"
   ic = 0
   do i=1,len_trim(f1)
      if(f1(i:i).ne.' ') then
         ic = ic+1
         if(f1(i:i).eq.'u') then
            f1(ic:ic) = ' '
         else
            f1(ic:ic) = f1(i:i)
         endif
      endif
   enddo
   f1(ic+1:) = ' '
!
!     ----- format for separator -----
!
   write(s2,*) "(a",nx2+2,":)"
   write(f2,*) "(1x,'",("-",i=1,nx1+2),"','|',",ncol,trim(s2),")"
   string(1:nx2+2) = repeat('-',nx2+2)
   ic = 0
   do i=1,len_trim(f2)
      if(f2(i:i).ne.' ') then
         ic = ic+1
         if(f2(i:i).eq.'u') then
            f2(ic:ic) = ' '
         else
            f2(ic:ic) = f2(i:i)
         endif
      endif
   enddo
   f2(ic+1:) = ' '
!
!     ----- format for data -----
!
   write(s3,*) "(1x,i",nx2+1," :)"
   write(f3,*) '(1x,i',nx1+1,",'u|',",ncol,trim(s3),")"
   ic = 0
   do i=1,len_trim(f3)
      if(f3(i:i).ne.' ') then
         ic = ic+1
         if(f3(i:i).eq.'u') then
            f3(ic:ic) = ' '
         else
            f3(ic:ic) = f3(i:i)
         endif
      endif
   enddo
   f3(ic+1:) = ' '
   return
end subroutine makfo4
subroutine makfo8(f1,f2,f3,ncol,string,mm,maxa)
!
!***********************************************************************
!     ----- impression d'une matrice R*4 ------
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   character (len=*), intent(out) :: f1,f2,f3,string
   integer, intent(in) :: mm
   integer (kind=K_I8), intent(in) :: maxa
   integer, intent(out) :: ncol
!
!     ----- local variables -----
!
   character (len=30) :: s1,s2,s3
   integer :: i,ic,nx1,nx2
   integer (kind=K_I8) :: maxl,maxt
!
!     +++++
!     start
!     +++++
!
   maxl = 1+int(log10(dble(mm)))
   maxt = 1+int(log10(dble(maxa)))
   ncol = (130-maxl)/(max(maxt,maxl)+2)
   nx1 = 1+maxl/10
   nx2 = 2+maxt/10
!
!     ----- format for column headers -----
!
   write(s1,*) "(1x,i",nx2,",1x:)"
   ic = 0
   do i=1,len_trim(s1)
      if(s1(i:i).ne.' ') then
         ic = ic+1
         s1(ic:ic) = s1(i:i)
      endif
   enddo
   s1(ic+1:) = ' '
   write(f1,*) "(//,",2+nx1,"x,'u|',",ncol,trim(s1),")"
   ic = 0
   do i=1,len_trim(f1)
      if(f1(i:i).ne.' ') then
         ic = ic+1
         if(f1(i:i).eq.'u') then
            f1(ic:ic) = ' '
         else
            f1(ic:ic) = f1(i:i)
         endif
      endif
   enddo
   f1(ic+1:) = ' '
!
!     ----- format for separator -----
!
   write(s2,*) "(a",nx2+2,":)"
   write(f2,*) "(1x,'",("-",i=1,nx1+2),"','|',",ncol,trim(s2),")"
   string(1:nx2+2) = repeat('-',nx2+2)
   ic = 0
   do i=1,len_trim(f2)
      if(f2(i:i).ne.' ') then
         ic = ic+1
         if(f2(i:i).eq.'u') then
            f2(ic:ic) = ' '
         else
            f2(ic:ic) = f2(i:i)
         endif
      endif
   enddo
   f2(ic+1:) = ' '
!
!     ----- format for data -----
!
   write(s3,*) "(1x,i",nx2+1," :)"
   write(f3,*) '(1x,i',nx1+1,",'u|',",ncol,trim(s3),")"
   ic = 0
   do i=1,len_trim(f3)
      if(f3(i:i).ne.' ') then
         ic = ic+1
         if(f3(i:i).eq.'u') then
            f3(ic:ic) = ' '
         else
            f3(ic:ic) = f3(i:i)
         endif
      endif
   enddo
   f3(ic+1:) = ' '
   return
end subroutine makfo8
end module matpri
