subroutine calchv(hvr,hvi,vr,vi,n,nvec,hr,hi,icol,hbuf,delta,d,mdress,  &
                   heff,work,sfile,elfile,elstar,sostar)
!
!***********************************************************************
!     ----- Compute H*|v> -----
!=======================================================================
!        HV (R8) = array containing all the H*v vectors.
!                        Input = contains NVEC vectors
!                        Output = contains H*v in column number NVEC+1
!         V (R8) = V vector (INPUT)
!      NVEC (I4) = number of H*v already calculated (INPUT)
!         H (R8) = lower triangle of the Hamiltonian (INPUT)
!      ICOL (I4) = position of the non-zero elements of the Hamiltonian
!                  matrix
!      HBUF (I4) = dimension of the buffers in H (INPUT)
!     DELTA (R8) = array containing the effective Hamiltonian
!         D (R8) = array containing the eigenvectors of the diagonalisation
!                  of H, without SO and shift (INPUT)
!    MDRESS (I4) = number of vectors to be dressed for the effective
!                  Hamiltonian
!      HEFF (L4) = flag for the calculation of the effective Hamiltonian
!     SFILE (I4) = unit number of the SO part of H (INPUT)
!    ELFILE (I4) = unit number of the non relativistic part of H (INPUT)
!     ELSTAR(I4) = array containing the 1st column of the first non-zero
!                  element of line I for the electrostatic part (INPUT)
!     SOSTAR(I4) = array containing the 1st column of the first non-zero
!                  element of line I for the S-O part (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   use kindvar
   use mconst
   use mdebug
   use f95_precision, only: wp => dp
   use blas95, only: gemm
   implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
   integer (kind=K_I4), intent(in) :: n,nvec,hbuf,mdress,sfile,elfile
   logical (kind=K_L4), intent(in) :: heff
   integer (kind=K_I4), intent(inout), dimension (hbuf) :: icol
   integer (kind=K_I4), intent(inout), dimension (:) :: elstar,sostar
   double precision, intent(in), dimension (:) :: delta
   double precision, intent(in), dimension (:,:) :: d
   double precision, intent(in), dimension (:,:) :: vr,vi
   double precision, intent(inout), dimension (hbuf) :: hr,hi
   double precision, intent(inout), dimension (:,:) :: work,hvr,hvi
!
!     ----- local variables -----
!
   integer (kind=K_I4) :: i,j,i1,nbj,ib1,ib2,nb2,elbuf
!
!     ----- statement functions -----
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
   call subent('CALCHV')
#endif
   hvr(1:n,1:nvec) = ZERO
   hvi(1:n,1:nvec) = ZERO
!
!     ------------------------------
!     ----- electrostatic part -----
!     ------------------------------
!
   rewind(elfile)
   read(elfile) hr,icol
!  i1 = elstar(1)
   elbuf = 0
   do i=1,n
      nbj = elstar(i+1)-elstar(i)
      ib1 = elbuf+1
      ib2 = elbuf+nbj
      if(ib2.le.hbuf) then
!
!     ----- the whole I-block is within the buffer -----
!
         nb2 = nbj
!
!     ----- HV(i,k) = H(i,j)*V(j,k) -----
!
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vr(icol(ib1:ib2),1:nvec),hvr(i:i,1:nvec),'N','N',ONE,ONE)
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vi(icol(ib1:ib2),1:nvec),hvi(i:i,1:nvec),'N','N',  &
                          ONE,ONE)
!
!     ----- HV(j,k) = H(j,i)*V(i,k) = H(i,j)*V(i,k) -----
!     ----- don't count H(i,i) twice -----
!
         if(nb2.gt.1) then
!           call gemm(reshape(hr(ib1:ib2),(/ nb2-1,1 /)),  &
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),  &
                      vr(i:i,1:nvec),work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
            hvr(icol(ib1:ib2-1),1:nvec) = hvr(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),  &
                      vi(i:i,1:nvec),work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
            hvi(icol(ib1:ib2-1),1:nvec) = hvi(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
         endif
!          ----- end if  skip diagonal term -----
         if(ib2.eq.hbuf) then
            read(elfile) hr,icol
!           nbj = nbj-nb2
            elbuf = 0
            cycle
         endif
      else
!
!     ----- the I-block spreads over two buffers -----
!     ----- 1) the part in the 1st buffer        -----
!
    1    continue
         ib2 = hbuf
         nb2 = ib2-(ib1-1)
!
!     ----- HV(i,k) = H(i,j)*V(j,k) -----
!
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vr(icol(ib1:ib2),1:nvec),hvr(i:i,1:nvec),'N','N',ONE,ONE)
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vi(icol(ib1:ib2),1:nvec),hvi(i:i,1:nvec),  &
                          'N','N',ONE,ONE)
!
!     ----- HV(j,k) = H(j,i)*V(i,k) = H(i,j)*V(i,k) -----
!     ----- don't count H(i,i) twice -----
!
         if(icol(ib2).lt.i) then
            call gemm(reshape(hr(ib1:ib2),(/ nb2,1 /)),  &
                      vr(i:i,1:nvec),work(1:nb2,1:nvec),'N','N',ONE,ZERO)
            hvr(icol(ib1:ib2),1:nvec) = hvr(icol(ib1:ib2),1:nvec)+  &
                                        work(1:nb2,1:nvec)
            call gemm(reshape(hr(ib1:ib2),(/ nb2,1 /)),vi(i:i,1:nvec),  &
                      work(1:nb2,1:nvec),'N','N',ONE,ZERO)
            hvi(icol(ib1:ib2),1:nvec) = hvi(icol(ib1:ib2),1:nvec)+  &
                                        work(1:nb2,1:nvec)
         else if(icol(ib2).eq.i) then
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),vr(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
            hvr(icol(ib1:ib2-1),1:nvec) = hvr(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),vi(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
            hvi(icol(ib1:ib2-1),1:nvec) = hvi(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
         endif
!          ----- end if  skip diagonal term -----
         read(elfile) hr,icol
         nbj = nbj-nb2
         elbuf = 0
!
!     ----- the I-block spreads over two buffers -----
!     ----- 2) the part in the 2nd buffer        -----
!
         ib1 = elbuf+1
         ib2 = nbj
         if(ib2.ge.hbuf) goto 1
         nb2 = ib2-(ib1-1)
!
!     ----- HV(i,k) = H(i,j)*V(j,k) -----
!
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vr(icol(ib1:ib2),1:nvec),hvr(i:i,1:nvec),'N','N',ONE,ONE)
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vi(icol(ib1:ib2),1:nvec),hvi(i:i,1:nvec),'N','N',ONE,ONE)
!
!     ----- HV(j,k) = H(j,i)*V(i,k) = H(i,j)*V(i,k) -----
!     ----- don't count H(i,i) twice -----
!
         if(nb2.gt.1) then
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),vr(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
            hvr(icol(ib1:ib2-1),1:nvec) = hvr(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),vi(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
            hvi(icol(ib1:ib2-1),1:nvec) = hvi(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
         endif
!          ----- end if  skip diagonal term -----
      endif
!          ----- end if  cut buffer ? -----
      elbuf = elbuf+nbj
   enddo
!
!     -------------------
!     ----- SO part -----
!     -------------------
!
   rewind sfile
   read(sfile) hr,hi,icol
!  i1 = sostar(1)
   elbuf = 0
   do i=1,n
      nbj = sostar(i+1)-sostar(i)
      ib1 = elbuf+1
      ib2 = elbuf+nbj
      if(ib2.le.hbuf) then
!
!     ----- the whole I-block is within the buffer -----
!
         nb2 = nbj
!
!     ----- HV(i,k) = H(i,j)*V(j,k) -----
!
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vr(icol(ib1:ib2),1:nvec),hvr(i:i,1:nvec),'N','N',ONE,ONE)
         call gemm(reshape(hi(ib1:ib2),(/ 1,nb2 /)),  &
                   vi(icol(ib1:ib2),1:nvec),hvr(i:i,1:nvec),  &
                          'N','N',-ONE,ONE)
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vi(icol(ib1:ib2),1:nvec),hvi(i:i,1:nvec),  &
                          'N','N',ONE,ONE)
         call gemm(reshape(hi(ib1:ib2),(/ 1,nb2 /)),  &
                          vr(icol(ib1:ib2),1:nvec),hvi(i:i,1:nvec),  &
                          'N','N',ONE,ONE)
!
!     ----- HV(j,k) = H(j,i)*V(i,k) = conjg(H(i,j))*V(i,k) -----
!     ----- don't count H(i,i) twice -----
!
         if(nb2.gt.1) then
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),  &
                      vr(i:i,1:nvec),work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
!           hvr(icol(ib1:ib2-1),1:nvec) = hvr(icol(ib1:ib2-1),1:nvec)+  &
!                                         work(1:nb2-1,1:nvec)
            call gemm(reshape(hi(ib1:ib2-1),(/ nb2-1,1 /)),  &
                      vi(i:i,1:nvec),work(1:nb2-1,1:nvec),'N','N',ONE,ONE)
!                           ONE,ZERO)
            hvr(icol(ib1:ib2-1),1:nvec) = hvr(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),  &
                      vi(i:i,1:nvec),work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
!           hvi(icol(ib1:ib2-1),1:nvec) = hvi(icol(ib1:ib2-1),1:nvec)+  &
!                                         work(1:nb2-1,1:nvec)
            call gemm(reshape(hi(ib1:ib2-1),(/ nb2-1,1 /)),  &
                      vr(i:i,1:nvec),work(1:nb2-1,1:nvec),'N','N',-ONE,ONE)
!                            ONE,ZERO)
!           hvi(icol(ib1:ib2-1),1:nvec) = hvi(icol(ib1:ib2-1),1:nvec)-  &
            hvi(icol(ib1:ib2-1),1:nvec) = hvi(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
         endif
         if(ib2.eq.hbuf) then
            read(sfile) hr,hi,icol
!:          nbj = nbj-nb2
            elbuf = 0
            cycle
         endif
      else
!
!     ----- the I-block spreads over two buffers -----
!     ----- 1) the part in the 1st buffer        -----
!
    2    continue
         ib2 = hbuf
         nb2 = ib2-(ib1-1)
!
!     ----- HV(i,k) = H(i,j)*V(j,k) -----
!
      call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                vr(icol(ib1:ib2),1:nvec),hvr(i:i,1:nvec),'N','N',ONE,ONE)
      call gemm(reshape(hi(ib1:ib2),(/ 1,nb2 /)),  &
                vi(icol(ib1:ib2),1:nvec),hvr(i:i,1:nvec),'N','N',-ONE,ONE)
      call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                vi(icol(ib1:ib2),1:nvec),hvi(i:i,1:nvec),'N','N',ONE,ONE)
      call gemm(reshape(hi(ib1:ib2),(/ 1,nb2 /)),  &
                vr(icol(ib1:ib2),1:nvec),hvi(i:i,1:nvec),'N','N',ONE,ONE)
!
!     ----- HV(j,k) = H(j,i)*V(i,k) = conjg(H(i,j))*V(i,k) -----
!     ----- don't count H(i,i) twice -----
!
         if(icol(ib2).lt.i) then
            call gemm(reshape(hr(ib1:ib2),(/ nb2,1 /)),vr(i:i,1:nvec),  &
                      work(1:nb2,1:nvec),'N','N',ONE,ZERO)
!           hvr(icol(ib1:ib2),1:nvec) = hvr(icol(ib1:ib2),1:nvec)+  &
!                                       work(1:nb2,1:nvec)
            call gemm(reshape(hi(ib1:ib2),(/ nb2,1 /)),vi(i:i,1:nvec),  &
                      work(1:nb2,1:nvec),'N','N',ONE,ONE)
!                     work(1:nb2,1:nvec),'N','N',ONE,ZERO)
            hvr(icol(ib1:ib2),1:nvec) = hvr(icol(ib1:ib2),1:nvec)+  &
                                        work(1:nb2,1:nvec)
            call gemm(reshape(hr(ib1:ib2),(/ nb2,1 /)),vi(i:i,1:nvec),  &
                      work(1:nb2,1:nvec),'N','N',ONE,ZERO)
!           hvi(icol(ib1:ib2),1:nvec) = hvi(icol(ib1:ib2),1:nvec)+  &
!                                       work(1:nb2,1:nvec)
            call gemm(reshape(hi(ib1:ib2),(/ nb2,1 /)),vr(i:i,1:nvec),  &
                      work(1:nb2,1:nvec),'N','N',-ONE,ONE)
!                     work(1:nb2,1:nvec),'N','N',ONE,ZERO)
!           hvi(icol(ib1:ib2),1:nvec) = hvi(icol(ib1:ib2),1:nvec)-  &
            hvi(icol(ib1:ib2),1:nvec) = hvi(icol(ib1:ib2),1:nvec)+  &
                                        work(1:nb2,1:nvec)
         else if(icol(ib2).eq.i) then
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),vr(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
!           hvr(icol(ib1:ib2),1:nvec) = hvr(icol(ib1:ib2),1:nvec)+  &
!                                   work(1:nb2-1,1:nvec)
            call gemm(reshape(hi(ib1:ib2-1),(/ nb2-1,1 /)),vi(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ONE)
!                        work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
            hvr(icol(ib1:ib2-1),1:nvec) = hvr(icol(ib1:ib2-1),1:nvec)+  &
                                        work(1:nb2-1,1:nvec)
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),vi(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
!           hvi(icol(ib1:ib2),1:nvec) = hvi(icol(ib1:ib2),1:nvec)+  &
!                                       work(1:nb2-1,1:nvec)
            call gemm(reshape(hi(ib1:ib2-1),(/ nb2-1,1 /)),vr(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',-ONE,ONE)
!                        work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
!           hvi(icol(ib1:ib2),1:nvec) = hvi(icol(ib1:ib2),1:nvec)-  &
            hvi(icol(ib1:ib2-1),1:nvec) = hvi(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
         endif
         read(sfile) hr,hi,icol
         nbj = nbj-nb2
         elbuf = 0
!
!     ----- the I-block spreads over two buffers -----
!     ----- 2) the part in the 2nd buffer        -----
!
         ib1 = elbuf+1
         ib2 = nbj
         if(ib2.ge.hbuf) goto 2
         nb2 = ib2-(ib1-1)
!
!     ----- HV(i,k) = H(i,j)*V(j,k) -----
!
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vr(icol(ib1:ib2),1:nvec),hvr(i:i,1:nvec),'N','N',ONE,ONE)
         call gemm(reshape(hi(ib1:ib2),(/ 1,nb2 /)),  &
                   vi(icol(ib1:ib2),1:nvec),hvr(i:i,1:nvec),'N','N',-ONE,ONE)
         call gemm(reshape(hr(ib1:ib2),(/ 1,nb2 /)),  &
                   vi(icol(ib1:ib2),1:nvec),hvi(i:i,1:nvec),'N','N',ONE,ONE)
         call gemm(reshape(hi(ib1:ib2),(/ 1,nb2 /)),  &
                   vr(icol(ib1:ib2),1:nvec),hvi(i:i,1:nvec),'N','N',ONE,ONE)
!
!     ----- HV(j,k) = H(j,i)*V(i,k) = conjg(H(i,j))*V(i,k) -----
!     ----- don't count H(i,i) twice -----
!
         if(nb2.gt.1) then
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),vr(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
!           hvr(icol(ib1:ib2),1:nvec) = hvr(icol(ib1:ib2),1:nvec)+  &
!                                       work(1:nb2-1,1:nvec)
            call gemm(reshape(hi(ib1:ib2-1),(/ nb2-1,1 /)),vi(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ONE)
!                     work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
            hvr(icol(ib1:ib2-1),1:nvec) = hvr(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),vi(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
!           hvi(icol(ib1:ib2),1:nvec) = hvi(icol(ib1:ib2),1:nvec)+  &
!                                       work(1:nb2-1,1:nvec)
            call gemm(reshape(hi(ib1:ib2-1),(/ nb2-1,1 /)),vr(i:i,1:nvec),  &
                      work(1:nb2-1,1:nvec),'N','N',-ONE,ONE)
!                     work(1:nb2-1,1:nvec),'N','N',ONE,ZERO)
!           hvi(icol(ib1:ib2),1:nvec) = hvi(icol(ib1:ib2),1:nvec)-  &
            hvi(icol(ib1:ib2-1),1:nvec) = hvi(icol(ib1:ib2-1),1:nvec)+  &
                                          work(1:nb2-1,1:nvec)
         endif
      endif
      elbuf = elbuf+nbj
   enddo
!
!     ------------------------
!     ----- Energy shift -----
!     ------------------------
!
   if(heff) then
      work(1:mdress,1:nvec) = ZERO
      call gemm(d(1:n,1:mdress),vr(1:n,1:nvec),work(1:mdress,1:nvec),  &
                'T','N',ONE,ZERO)
      do j=1,nvec
         work(1:mdress,j) = work(1:mdress,j)*delta(1:mdress)
      enddo
      call gemm(d(1:n,1:mdress),work(1:mdress,1:nvec),hvr(1:n,1:nvec), &
                'N','N',ONE,ONE)
      work(1:mdress,1:nvec) = ZERO
      call gemm(d(1:n,1:mdress),vi(1:n,1:nvec),work(1:mdress,1:nvec),  &
                'T','N',ONE,ZERO)
      do j=1,nvec
         work(1:mdress,j) = work(1:mdress,j)*delta(1:mdress)
      enddo
      call gemm(d(1:n,1:mdress),work(1:mdress,1:nvec),hvi(1:n,1:nvec), &
                'N','N',ONE,ONE)
   endif
#ifdef _MYDEBUG_
   call subret('CALCHV',3)
#endif
   return
end subroutine calchv
