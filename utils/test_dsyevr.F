!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

!     Test program for lapack symmetric eigensolver DSYEVR
!     By ulfek 2008
      program test_dsyevr
      implicit none
      integer i, j, k, N, N2, N3, lwork, liwork, info, idum, M,
     &     idum2,idum3,idum4
!     Store matrices in contigous memory blocks with (:,:) ?
      double precision, allocatable :: A(:,:),vec(:,:),eig(:)
      double precision, allocatable :: work(:)
      double precision ddum, workdum, maxerr, maxvecerr, pi, verr,
     &     ddum2
      integer, allocatable :: iwork(:), isuppz(:)
      character(len=32) :: arg
      logical issorted
      if (iargc().gt.0) then
         call getarg(1,arg)
         read (arg,*), N
      else
         N = 19
      endif
      pi = 2*acos(0.0D0)
      print *,'Testing DSYEVR diagonalization of',N,'**2 matrix'
#ifdef FORTRAN2003
      print *,'Using',sizeof(N)*8,' bit integers'
      print *,'Allocating matrix and vectors ',
     &     ((2*N**2 + 3*N)*sizeof(ddum))/2**10,' kB'
#endif
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
#ifdef DEBUG
      print *,'Matrix to diagonalize:'
      do i=1,N
         write (6,'(100F4.1)') (A(i,j),j=1,N)
      enddo
#endif
      N2 = N
      N3 = N
      CALL DSYEVR('V','A','U',N,A,N2,DDUM,DDUM2,IDUM,IDUM2,0.0D0,
     &     IDUM3,EIG,VEC,N3,IDUM4,WORKDUM,-1,liwork,-1,INFO)
      if (info.ne.0) then
         print *, 'DSYEVR memory query failed:',info
         stop
      endif
      lwork = workdum
#ifdef FORTRAN2003
      print *,'Allocating work buffers',(lwork*sizeof(ddum)+
     &     liwork*sizeof(idum))/2**10,' kB'
#endif
      call flush(6)
      allocate( work(lwork), iwork(liwork) )
      print *,'Diagonalizing..'
      call flush(6)
      CALL DSYEVR('V','A','U',N,A,N2,DDUM,DDUM2,IDUM,IDUM2,0.0D0,
     &     M,EIG,VEC,N3,isuppz,work,lwork, iwork,liwork,info)
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
#ifdef DEBUG
            print *,'vector ',i,' eigval:',eig(i)
            print *,'  <Calculated>             <Exact>'
#endif
            do k=1,N
!     Discard the phase of vectors with abs
               verr = abs(eig(i)*(abs(vec(k,i)) -
     &              abs(sqrt(2.0D0/(N+1))*sin((j*k*pi)/(N+1)))))
               maxvecerr = max(maxvecerr, verr)
#ifdef DEBUG
               if (verr.gt.N*1.0D-14) then
                  print *,vec(k,i),sqrt(2.0D0/(N+1))*sin((j*k*pi)/(N+1))
     &                 ,'<<<<< WARNING'
               else
                  print *,vec(k,i),sqrt(2.0D0/(N+1))*sin((j*k*pi)/(N+1))
               endif
#endif
            enddo
#ifdef DEBUG
            print *,'maxvecerr so far:',maxvecerr
#endif
         enddo
         print *,'Maximum absolute eigenvalue error:',maxerr
         print *,'Maximum absolute vector element error:',maxvecerr
      else
         print *,'Eigenvalues not sorted so not testing them,'//
     &        ' this can be ok.'
      endif
      end
