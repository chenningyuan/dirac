!     Test program for lapack symmetric eigensolver DSYEVR
!     By ulfek 2008


      program test_dsyevr
      implicit none
      integer i, j, k, N, lwork, liwork, info, idum, M
!     Store matrices in contigous memory blocks with (:,:) ?
      double precision, allocatable :: A(:,:),vec(:,:),eig(:)
      double precision, allocatable :: work(:)
      double precision ddum, workdum, maxerr, maxvecerr, pi
      integer, allocatable :: iwork(:), isuppz(:)
      logical issorted
      N = 234
      pi = 2*acos(0.0D0)
      print *,'Testing DSYEVR diagonalization of',N,'**2 matrix'
      print *,'Using',sizeof(N)*8,' bit integers'
      print *,'Allocating matrix and vectors ',
     &     ((2*N**2 + 3*N)*sizeof(ddum))/2**20,' MB'
      call flush(6)
      allocate( A(N,N), vec(N,N), eig(N), isuppz(2*N) )
!     set A = 0.5 on the upper and lower second diagonals
      do i=1,N
         do j=1,N
            A(i,j) = 0
         enddo
      enddo
      do i=1,N-1
         A(i,i+1) = 0.5
         A(i+1,i) = 0.5
      enddo

      CALL DSYEVR('V','A','U',N,A,N,DDUM,DDUM,IDUM,IDUM,0.0D0,
     &     IDUM,EIG,VEC,N,IDUM,WORKDUM,-1,liwork,-1,INFO)
      if (info.ne.0) then
         print *, 'DSYEVR memory query failed:',info
         stop
      endif
      lwork = workdum
      print *,'Allocating work buffers',(lwork*sizeof(ddum)+
     &     liwork*sizeof(idum))/2**10,' kB'
      call flush(6)
      allocate( work(lwork), iwork(liwork) )
      print *,'Diagonalizing..'
      call flush(6)
      CALL DSYEVR('V','A','U',N,A,N,DDUM,DDUM,IDUM,IDUM,0.0D0,
     &     M,EIG,VEC,N,isuppz,work,lwork, iwork,liwork,info)
      if (info.ne.0.or.M.ne.N) then
         print *,'Diagonalization failed, info = ',info,' M =',M
         stop
      endif
      issorted = .true.
      do i=1,N-1
         if (eig(i).gt.eig(i+1)) then
            issorted = .false.
         endif
      enddo
      print *,'Eigenvalues sorted ascending:',issorted
      if (issorted) then
!     Check eigenvalues and vectors against analytical solution
!     Currently only implemented for ascending sorting
         maxerr = 0
         maxvecerr = 0
         do i=1,N
            j = (N-i+1)
            maxerr = max(maxerr,abs(eig(i) - cos((j*pi)/(N+1))))
!     Discard the phase of vectors with abs
            do k=1,N
               maxvecerr = max(maxvecerr, abs(abs(vec(k,j)) -
     &              abs(sqrt(2.0D0/(N+1))*sin((j*k*pi)/(N+1)))))
            enddo
         enddo
         print *,'Maximum absolute eigenvalue error:',maxerr
         print *,'Maximum absolute vector element error:',maxvecerr
      else
         print *,'Eigenvalues not sorted so not testing them'
         print *,'This is not an error in itself'
      endif
      end
