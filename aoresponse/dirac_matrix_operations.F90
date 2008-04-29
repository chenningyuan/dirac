module dirac_matrix_operations_m

! author:           Andreas Johan Thorvaldsen
! porting to DIRAC: Radovan Bast




!Manual allocation and deallocation of matrices is cumbersome to program.
!This module does book-keeping of the allocated data and their corresponding
!references (aliases), and performs deallocation once all referrers have
!come out of scope (on the stack) or been overwritten.
!
!Arithmetic operators are: A=B (equals), A+B (plus), A-B (minus),-A (negative),
!A*B (matmul), r*A (scale by real*8), i*A (scale by integer), A/r (divide by real*8),
!A/i (divide by integer), A**B (matrix dot product), trps(A) (transpose),
!A=0 (equals zero - deallocate)
!
!This file must be compiled without inlining

    use kinds_m
    use memory_allocator
    implicit none

#include "dgroup.h"

   TYPE Matrix
      !number of rows, number of columns
!ajt mark uninitialized with -1      INTEGER :: nrow, ncol
      INTEGER :: nrow=-1, ncol
      !room for various flags/information for the matrix
      integer, dimension(:),pointer :: idata
      !pointer to dob.prec. matrix element storage
      real(kind_dp), pointer   :: elms(:)
      !pointer to dob.prec. matrix element storage for beta part
      real(kind_dp), pointer   :: elmsb(:)
      !pointer to complex matrix element storage
      complex(kind_dp), pointer :: celms(:)
      !pointer to complex matrix element storage for beta part
      complex(kind_dp), pointer :: celmsb(:)
      !pointer to storage scheme 1 for sparse matrices
!     type(SMAT1),pointer :: selm1(:)
      !pointer to storage scheme 2 for sparse matrices
      !type(SMAT2),pointer :: selm2
      !If only blocks of the matrix are to be stored:
      !array of pointers to the matrix-blocks
      type(Matrix),pointer :: block(:)
      !the indexes hvere the blocks start
      integer,dimension(:,:), pointer :: blockpos!(2,:)
      !room for any auxaliary information
      integer, pointer     :: iaux(:)
      real(kind_dp), pointer :: raux(:)
   END TYPE Matrix



    public Matrix,print,trps,norm,iszero,isdef,MatrixStats
    public mat_init, mat_free, printmat
    public assignment(=),operator(+),operator(-)
    public operator(*),operator(/),operator(**),trab
    private !rest is private

    interface assignment(=)
        module procedure MatEqMat,MatEqZero
    end interface

    interface operator(+)
        module procedure MatPlusMat
    end interface

    interface operator(-)
        module procedure MatMinusMat,MatNeg
    end interface

    interface operator(*)
        module procedure RealTimesMat,IntTimesMat,MatMulMat
    end interface

    interface operator(/)
        module procedure MatDivReal,MatDivInt
    end interface

    interface operator(**)
        module procedure MatDotMat
    end interface

    interface trab
        module procedure TrMatMulMat
    end interface

    interface trps
        module procedure TrpsMat
    end interface

    interface norm
        module procedure NormMat
    end interface

    interface iszero
        module procedure IsZeroMat
    end interface

    interface isdef
        module procedure IsDefMat
    end interface

    interface print
        module procedure PrintMat
    end interface

    !this comes from 'use matrix_module'
    ! type Matrix !reference to matrix
    !     integer:: r=-1,c !rows,columns. r=-1 means undef
    !     real*8,pointer:: d(:) !data, if defined and non-zero
    !     !real*8:: f scale factor
    !     !character*8:: k kind, N/T: transposed, G/S/A: gen,sym,anti,
    !                      !C/O/S/T: closed-,open-shell,singlet,triplet,
    !                      !F/S: full,sparse, L/D: local,distributed
    ! end type

    !no array-of-pointer in Fortran, so use array
    !of derived type containing pointer instead
    type refp
        type(Matrix),pointer:: m
    end type

    type datp
        real*8,pointer:: d(:)
    end type

    !number of references, current book-keeping buffer size
    integer:: nref=0,maxref=0
    !number of matrix datas allocated, current book-keeping buffer size
    integer:: ndat=0,maxdat=0
    !growth direction of the stack, -1 or +1
    integer:: stkdir !different for i686,Itanium,XEON, fex.

    !pointers to all type(Matrix) referring to a matrix data (%elms)
    type(refp),allocatable:: ref(:) !1:numref<maxref
    !one pointer to each of the allocated %elms(:)
    type(datp),allocatable:: dat(:) !1:numdat<maxdat
    !index in dat referred to by each ref
    integer,allocatable:: idx(:)
    !number of references to each dat
    integer,allocatable:: cnt(:)


contains

#ifndef VAR_IFC
#define baddress loc
#endif /* VAR_IFC */

  subroutine mat_init(a,nrow,ncol)
     implicit none
     type(matrix) :: a
     integer, intent(in) :: nrow, ncol
     integer :: nsize

     nullify(a%elms)
     a%nrow = nrow
     a%ncol = ncol
     nsize = a%nrow*a%ncol
     call alloc(a%elms, nsize)
  end subroutine

  subroutine mat_free(a)
     implicit none
     TYPE(Matrix) :: a
     integer :: nsize
     if (.not.ASSOCIATED(a%elms)) then
       print*,'memory previously released!!'
       STOP 'Error in mat_dense_free - memory previously released'
     endif
     nsize = SIZE(a%elms)
     call dealloc(a%elms)
     NULLIFY(a%elms)
  end subroutine

    subroutine allocif()
    !if needed allocate/extend the tables ref,dat,idx,cnt
        type(refp):: tref(nref); type(datp):: tdat(maxdat)
        integer:: tidx(nref),tcnt(maxdat),i
        if (nref<maxref-1.and.ndat<maxdat-1) return
        if (nref>=maxref-1) then
            if (maxref==0) call getdir() !only once
            if (maxref/=0) tref = ref(1:nref)
            if (maxref/=0) tidx = idx(1:nref)
            if (maxref/=0) deallocate(ref,idx)
            allocate(ref(max(64,2*maxref)))
            allocate(idx(max(64,2*maxref)))
            ref(1:nref) = tref
            idx(1:nref) = tidx
            maxref = max(64,2*maxref)
        endif
        if (ndat>=maxdat-1) then
            if (maxdat/=0) tdat = dat(1:maxdat)
            if (maxdat/=0) tcnt = cnt(1:maxdat)
            if (maxdat/=0) deallocate(dat,cnt)
            allocate(dat(max(8,2*maxdat)))
            allocate(cnt(max(8,2*maxdat))); cnt=0
            dat(1:maxdat) = tdat
            cnt(1:maxdat) = tcnt
            maxdat = max(8,2*maxdat)
        endif
    contains
        subroutine getdir()
            integer:: j
            stkdir = merge(-1,1,baddress(j)<baddress(i))
        end subroutine
    end subroutine


    subroutine delabove(k)
        integer:: k,i,j
        call allocif()
        do i=nref,1,-1
            if ((baddress(ref(i)%m)-baddress(k))*stkdir<0) exit
            j = idx(i)
            cnt(j) = cnt(j)-1
            if (cnt(j)==0) then
                deallocate(dat(j)%d)
                ndat = ndat-1
            endif
        enddo
        nref = i
    end subroutine


    subroutine insref(A)
        type(Matrix),intent(in),target:: A
        integer:: i,j,k; logical:: found,zero,same
        !see if any reference from A is recorded
        do i=nref,0,-1
            if (i==0) exit
            if ((baddress(ref(i)%m)-baddress(A))*stkdir<0) exit
        enddo; j=i+1; k=idx(j)
        found  = (j/=nref+1.and.associated(ref(j)%m,A))
        zero = (.not.associated(A%elms).or.abs(baddress(A%elms)-baddress(i))<2**28)
        same = .false.; if (found) same = associated(dat(k)%d,A%elms)
        !if A is recorded, but now is zero or A%elms is on the stack, remove the previous record
        if (found.and.zero) then
             ref(j:nref-1) = ref(j+1:nref)
             idx(j:nref-1) = idx(j+1:nref)
             nref = nref-1
        endif
        !if A was recorded as referring to some other data, remove that
        if (found.and..not.same) then
            cnt(k) = cnt(k)-1
            if (cnt(k)==0) deallocate(dat(k)%d)
            if (cnt(k)==0) ndat = ndat-1
        endif
        !see if A%elms is already recorded as some dat(k)%d
        do i=maxdat,1,-1
            if (cnt(i)==0) then
                k = i; cycle
            elseif (associated(A%elms,dat(i)%d)) then
                k = i; exit
            endif
        enddo
        !if A is not already recorded, make space for it
        if (.not.found.and..not.zero) then
            ref(j+1:nref+1) = ref(j:nref)
            idx(j+1:nref+1) = idx(j:nref)
            nref = nref+1
        endif
        !insert A
        if (.not.zero) then
            ref(j)%m => A
            idx(j) = k
            dat(k)%d => A%elms
            if (cnt(k)==0) ndat = ndat+1
            if (.not.same) cnt(k) = cnt(k)+1
        endif
    end subroutine

    subroutine delref(A)
    !if A is not recorded, deallocate A%elms
    !if A is recorded as referring to dat(k)%d, but refers to A%elms - deallocate(A%elms)
    !if A is recorded, referring to recorded dat(k)%d - decrement refcount for dat(k)%d
        type(Matrix),intent(inout),target:: A
        integer:: i,j,k
        !find the reference, if exists
        do j=nref,1,-1
            if ((baddress(ref(j)%m)-baddress(A))*stkdir<=0) exit
        enddo; k=idx(j)
        !if A%elms is not the recorded dat(k)%d, see if A%elms is referenced by some
        !other matrix. If not, deallocate A%elms
        if (A%nrow/=-1.and.associated(A%elms).and..not. &
            (associated(ref(j)%m,A).and.associated(dat(k)%d,A%elms))) then
            do i=maxdat,0,-1
                if (cnt(i)==0) cycle
                if (associated(A%elms,dat(i)%d)) exit
                if (i==1) deallocate(A%elms)
            enddo
        endif
        !if a reference from A was found, remove it
        if (associated(ref(j)%m,A)) then
            ref(j:nref-1) = ref(j+1:nref)
            idx(j:nref-1) = idx(j+1:nref)
            nref = nref-1
            cnt(k) = cnt(k)-1
            if (cnt(k)==0) then
                deallocate(dat(k)%d)
                ndat = ndat-1
            endif
        endif
        A%nrow=-1; nullify(A%elms)
    end subroutine

    subroutine MatEqMat(A,B)
        type(Matrix),intent(inout):: A
        type(Matrix),intent(in):: B; integer:: i
        call delabove(i)
        if (B%nrow==-1) stop 'Matrix A=B : B undefined'
        if (baddress(A)/=baddress(B)) call delref(A)
        A%nrow = B%nrow
        A%ncol = B%ncol
        A%elms => B%elms
        call insref(A)
    end subroutine

    subroutine MatEqZero(A,z)
        type(Matrix),intent(inout):: A
        integer, intent(in) :: z; integer:: i,j
        call delabove(i)
        if (A%nrow==-1) stop 'Matrix A=0 : A undefined'
        if (z/=0) stop 'Matrix A=k : k must be zero (clear)'
        i=A%nrow; j=A%ncol
        call delref(A)
        A%nrow=i; A%ncol=j !A%elms is nullified, meaning zero
    end subroutine

    function MatPlusMat(A,B) result(C)
        type(Matrix),intent(in):: A,B
        type(Matrix):: C; integer:: i
        call delabove(i)
        if (A%nrow==-1) stop 'Matrix A+B : A undefined'
        if (B%nrow==-1) stop 'Matrix A+B : B undefined'
        if (A%nrow/=B%nrow.or.A%ncol/=B%ncol) stop 'Matrix A+B : different shapes'
        C%nrow=A%nrow; C%ncol=A%ncol; nullify(C%elms)
        if (associated(A%elms).and.associated(B%elms)) then
            allocate(C%elms(C%nrow*C%ncol))
            C%elms = A%elms+B%elms
        elseif (associated(A%elms)) then
            C%elms => A%elms
        elseif (associated(B%elms)) then
            C%elms => B%elms
        endif
        call insref(C)
    end function

    function MatMinusMat(A,B) result(C)
        type(Matrix),intent(in):: A,B
        type(Matrix):: C
        C = A + (-B)
    end function

    function RealTimesMat(r,A) result(B)
        type(Matrix),intent(in):: A
        real*8,intent(in):: r; type(Matrix):: B
        integer:: i
        call delabove(i)
        if (A%nrow==-1) stop 'Matrix r*A : A undefined'
        B%nrow=A%nrow; B%ncol=A%ncol; nullify(B%elms)
        if (r/=0.and.associated(A%elms)) then
            allocate(B%elms(B%nrow*B%ncol))
            B%elms = r*A%elms
        endif
        call insref(B)
    end function

    function MatNeg(A) result(B)
        type(Matrix),intent(in):: A
        type(Matrix):: B
        B = (-1d0)*A
    end function

    function IntTimesMat(i,A) result(B)
        type(Matrix),intent(in):: A
        integer,intent(in):: i; type(Matrix):: B
        B = (1d0*i)*A
    end function

    function MatDivReal(A,r) result(B)
        type(Matrix),intent(in):: A
        real*8,intent(in):: r; type(Matrix):: B
        B = (1d0/r)*A
    end function

    function MatDivInt(A,i) result(B)
        type(Matrix),intent(in):: A
        integer,intent(in):: i; type(Matrix):: B
        B = (1d0/i)*A
    end function

    function MatMulMat(A,B) result(C)
        type(Matrix),intent(in):: A,B
        type(Matrix):: C
        integer:: i,k,l,m,n
        call delabove(i)
        if (A%nrow==-1) stop 'Matrix A*B : A undefined'
        if (B%nrow==-1) stop 'Matrix A*B : B undefined'
        if (A%ncol/=B%nrow) stop 'Matrix A*B : incompatible shapes'
        C%nrow=A%nrow; C%ncol=B%ncol
        if (associated(A%elms).and.associated(B%elms)) then
            allocate(C%elms(C%nrow*C%ncol))
            C%elms = reshape(matmul(reshape(A%elms(1:A%nrow*A%ncol),(/A%nrow,A%ncol/)), &
                     reshape(B%elms(1:B%nrow*B%ncol),(/B%nrow,B%ncol/))),(/C%nrow*C%ncol/))
        endif
        call insref(C)
    end function

    function MatDotMat(A,B) result(r)
        type(Matrix),intent(in):: A,B; real*8:: r
        if (A%nrow==-1) stop 'Matrix A**B : A undefined'
        if (B%nrow==-1) stop 'Matrix A**B : B undefined'
        if (A%nrow/=B%nrow) stop 'Matrix A**B : different number of rows'
        if (A%ncol/=B%ncol) stop 'Matrix A**B : different number of cols'
        r = 0d0
        if (.not.(associated(A%elms).and.associated(B%elms))) return
        r = sum(A%elms*B%elms)
    end function

    function TrMatMulMat(A,B) result(r)
        type(Matrix),intent(in):: A,B; real*8:: r
        if (A%nrow==-1) stop 'Matrix trab(A,B) : A undefined'
        if (B%nrow==-1) stop 'Matrix trab(A,B) : B undefined'
        if (A%nrow/=B%ncol) stop 'Matrix trab(A,B) : not compatible shapes'
        if (A%ncol/=B%nrow) stop 'Matrix trab(A,B) : not compatible shapes'
        r = 0d0
        if (.not.(associated(A%elms).and.associated(B%elms))) return
        r = sum(transpose(reshape(A%elms,(/A%nrow,A%ncol/)))*reshape(B%elms,(/B%nrow,B%ncol/)))
    end function

    function TrpsMat(A) result(B)
        type(Matrix),intent(in):: A; type(Matrix):: B
        integer:: i; call delabove(i)
        if (A%nrow==-1) stop 'Matrix trps(A) : A undefined'
        B%nrow=A%ncol; B%ncol=A%nrow
        if (associated(A%elms)) then
            allocate(B%elms(B%nrow*B%ncol))
            B%elms = reshape(transpose(reshape(A%elms,(/A%nrow,A%ncol/))),(/B%nrow*B%ncol/))
        endif
        call insref(B)
    end function

    function NormMat(A) result(r)
        type(Matrix),intent(in):: A; real*8:: r
        if (A%nrow==-1) stop 'Matrix A**B : A undefined'
        r = 0d0
        if (.not.associated(A%elms)) return
        r = sqrt(sum(A%elms**2))
    end function

    function IsZeroMat(A) result(z)
        type(Matrix),intent(in):: A; logical:: z
        integer:: i; call delabove(i)
        if (A%nrow==-1) stop 'Matrix iszero(A) : A undefined'
        z = (A%nrow/=-1.and..not.associated(A%elms))
    end function

    function IsDefMat(A) result(d)
        type(Matrix),intent(in):: A; logical:: d
        integer:: i; call delabove(i)
        d = (A%nrow/=-1)
    end function

    subroutine PrintMat(A,w)
        type(Matrix),intent(in):: A; integer:: i,ww,d
        integer,optional:: w
#include "priunit.h"
        call delabove(i)
        if (A%nrow==-1) stop 'Matrix print(A) : A undefined'
        if (associated(A%elms)) then

    call prqmat(A%elms, &
                A%nrow, &
                A%ncol, &
                A%nrow, &
                A%ncol, &
                nz,               &
                ipqtoq(1, 0),     &
                lupri)
            d = max(1,ceiling(log10(max(maxval(A%elms),-10*minval(A%elms)))))
            ww = 8; if (present(w)) ww = max(w,d+2); d = ww-d-2
#ifdef VAR_IFC
            print '(<A%nrow>(<A%ncol>f<ww>.<d>/))',(A%elms(i:i+A%nrow*(A%ncol-1):A%nrow),i=1,A%nrow)
#endif /* VAR_IFC */
        else !zero
            ww = 8; if (present(w)) ww = max(w,1+2); d = ww-1-2
#ifdef VAR_IFC
            print '(<A%nrow>(<A%ncol>f<ww>.<d>/))',(0d0,i=1,A%nrow*A%ncol)
#endif /* VAR_IFC */
        endif
    end subroutine

    subroutine MatrixStats()
        integer:: i,j
        call delabove(i)
        print '("ndat="i2" maxdat="i2" nref="i2" maxref="i2)', &
                ndat,maxdat,nref,maxref
        i = maxval(idx(1:nref))
        if (nref==0.or.ndat==0) return
#ifdef VAR_IFC
        print '("ref="<nref>i13/"dat="<i>i13/"idx="<nref>i2"  cnt="<i>i2/)', &
                (baddress(ref(j)%m),j=1,nref),(baddress(dat(j)%d),j=1,i), &
                idx(1:nref),cnt(1:i)
#endif /* VAR_IFC */
    end subroutine

end module dirac_matrix_operations_m

! program prog
!     use TypeMatrix
!     type(Matrix):: A
!     integer,parameter:: n=5
!     call subr()
!     A = 0            ; print '("A= "i12)',baddress(A)              ; call print(A); call MatrixStats()
! contains
!     subroutine subr()
!         type(Matrix):: B,W,X,Y,Z
!         integer:: i,j,k; integer,parameter:: n=5
!         real*8,pointer:: d(:),e(:),f(:)
!         allocate(d(n*n),e(n*n),f(n*n))
!         e = (/((abs(i-j)*1d0,j=1,n),i=1,n)/)
!         d = (/(1d0,(0d0,j=1,n),i=1,n-1),1d0/)
!         f = (/((i-j,j=1,n),i=1,n)/)
!         A = Matrix(n,n,d); print '("A="2i12)',baddress(A),baddress(A%d); call print(A); call MatrixStats()
!         B = Matrix(n,n,e); print '("B="2i12)',baddress(B),baddress(B%d); call print(B); call MatrixStats()
!         W = Matrix(n,n,f); print '("W="2i12)',baddress(W),baddress(W%d); call print(W); call MatrixStats()
!         X = A            ; print '("X="2i12)',baddress(X),baddress(X%d); call print(X); call MatrixStats()
!         Y = B            ; print '("Y="2i12)',baddress(Y),baddress(Y%d); call print(Y); call MatrixStats()
!         Z = W            ; print '("Z="2i12)',baddress(Z),baddress(Z%d); call print(Z); call MatrixStats()
!         Z = func(A,B)    ; print '("Z="2i12)',baddress(Z),baddress(Z%d); call print(Z); call MatrixStats()
!         A = 0            ; print '("A= "i12)',baddress(A)              ; call print(A); call MatrixStats()
!         B = 0            ; print '("B= "i12)',baddress(B)              ; call print(B); call MatrixStats()
!         W = 0            ; print '("W= "i12)',baddress(W)              ; call print(W); call MatrixStats()
!         X = 0            ; print '("X= "i12)',baddress(X)              ; call print(X); call MatrixStats()
!         Y = 0            ; print '("Y= "i12)',baddress(Y)              ; call print(Y); call MatrixStats()
!         Z = 0            ; print '("Z= "i12)',baddress(Z)              ; call print(Z); call MatrixStats()
!     end subroutine
!     function func(A,B) result(C)
!         type(Matrix),intent(in):: A,B
!         type(Matrix):: C,D
!         D = A*B; print '("D="2i12)',baddress(D),baddress(D%d); call print(D); call MatrixStats()
!         C = D*D; print '("C="2i12)',baddress(C),baddress(C%d); call print(C); call MatrixStats()
!     end function
! end program
