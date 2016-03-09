subroutine dclchv(hv,vr,n,nvec,hr,icol,hbuf,ifile,elstar,work)
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
!     ELSTAR(I4) = array containing the 1st column of the first non-zero
!                  element of line I for the electrostatic part (INPUT)
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
#ifdef _MYDEBUG_
   use mdebug
#endif
   use f95_precision, only: wp => dp
   use blas95, only: gemm
   implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
   integer (kind=K_I4), intent(in) :: n,nvec,hbuf,ifile
   integer (kind=K_I4), intent(inout), dimension (hbuf) :: icol
   integer (kind=K_I4), intent(inout), dimension (:) :: elstar
   double precision, intent(in), dimension (:,:) :: vr
   double precision, intent(out), dimension (:,:) :: hv,work
   double precision, intent(inout), dimension (hbuf) :: hr
!
!     ----- local variables -----
!
   integer (kind=K_I4) :: i,i1,nbj,nb2,elbuf,ib1,ib2
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
   call subent('DCLCHV')
#endif
   hv(1:n,1:nvec) = ZERO
   rewind(ifile)
   read(ifile) hr,icol
   i1 = elstar(1)
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
         call gemm(reshape(hr(ib1:ib2), (/ 1,nb2 /)),  &
                   vr(icol(ib1:ib2),1:nvec),hv(i:i,1:nvec),'N','N',  &
                   ONE,ONE)
!#          do j=ib1,ib2
!#             hv(i,1:nvec) = hr(j)*vr(icol(j),1:nvec)+hv(i,1:nvec)
!#          enddo
!
!     ----- HV(j,k) = H(j,i)*V(i,k) = H(i,j)*V(i,k) -----
!     ----- don't count H(i,i) twice -----
!
         if(nb2.gt.1) then
           call gemm(reshape(hr(ib1:ib2-1), (/ nb2-1,1 /)),  &
                            vr(i:i,1:nvec),work(1:nb2-1,1:nvec),'N','N',  &
                            ONE,ZERO)
           hv(icol(ib1:ib2-1),1:nvec) = hv(icol(ib1:ib2-1),1:nvec)+  &
                                        work(1:nb2-1,1:nvec)
!#           do j=ib1,ib2-1
!#              hv(icol(j),1:nvec) = hr(j)*vr(i,1:nvec)+hv(icol(j),1:nvec)
!#           enddo
         endif
         if(ib2.eq.hbuf) then
            read(ifile) hr,icol
            elbuf = 0
!           i1 = elstar(i+1)
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
         call gemm(reshape(hr(ib1:ib2), (/ 1,nb2 /)),  &
                   vr(icol(ib1:ib2),1:nvec),hv(i:i,1:nvec),'N','N',  &
                   ONE,ONE)
!#         do j=ib1,ib2
!#            hv(i,1:nvec) = hr(j)*vr(icol(j),1:nvec)+hv(i,1:nvec)
!#         enddo
!
!     ----- HV(j,k) = H(j,i)*V(i,k) = H(i,j)*V(i,k) -----
!     ----- don't count H(i,i) twice -----
!
         if(icol(ib2).lt.i) then
          call gemm(reshape(hr(ib1:ib2),(/ nb2,1 /)),  &
                    vr(i:i,1:nvec),work(1:nb2,1:nvec),'N','N',  &
                    ONE,ZERO)
           hv(icol(ib1:ib2),1:nvec) = hv(icol(ib1:ib2),1:nvec)+  &
                                        work(1:nb2,1:nvec)
!#           do j=ib1,ib2
!#              hv(icol(j),1:nvec) = hr(j)*vr(i,1:nvec)+hv(icol(j),1:nvec)
!#           enddo
         else if(icol(ib2).eq.i) then
           call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),  &
                     vr(i:i,1:nvec),work(1:nb2-1,1:nvec),'N','N',  &
                     ONE,ZERO)
           hv(icol(ib1:ib2-1),1:nvec) = hv(icol(ib1:ib2-1),1:nvec)+  &
                                        work(1:nb2-1,1:nvec)
!#           do j=ib1,ib2-1
!#              hv(icol(j),1:nvec) = hr(j)*vr(i,1:nvec)+hv(icol(j),1:nvec)
!#           enddo
         endif
         read(ifile) hr,icol
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
        call gemm(reshape(hr(ib1:ib2), (/ 1,nb2 /)),  &
                  vr(icol(ib1:ib2),1:nvec),hv(i:i,1:nvec),'N','N',  &
                  ONE,ONE)
!#         do j=ib1,ib2
!#            hv(i,1:nvec) = hr(j)*vr(icol(j),1:nvec)+hv(i,1:nvec)
!#         enddo
!
!     ----- HV(j,k) = H(j,i)*V(i,k) = H(i,j)*V(i,k) -----
!     ----- don't count H(i,i) twice -----
!
         if(nb2.gt.1) then
            call gemm(reshape(hr(ib1:ib2-1),(/ nb2-1,1 /)),  &
                      vr(i:i,1:nvec),work(1:nb2-1,1:nvec),'N','N',  &
                      ONE,ZERO)
            hv(icol(ib1:ib2-1),1:nvec) = hv(icol(ib1:ib2-1),1:nvec)+  &
                                         work(1:nb2-1,1:nvec)
!#           do j=ib1,ib2-1
!#              hv(icol(j),1:nvec) = hr(j)*vr(i,1:nvec)+hv(icol(j),1:nvec)
!#           enddo
         endif
      endif
      elbuf = elbuf+nbj
   enddo
#ifdef _MYDEBUG_
   call subret('DCLCHV',2)
#endif
   return
end subroutine dclchv

