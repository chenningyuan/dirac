subroutine hred(vhv,hr,isel,icol,n,hbuf,elfile,elstar)
!
!***********************************************************************
!     ----- Get the relevant Hamiltonian submatrix on the space -----
!     ----- of selected determinants                            -----
!=======================================================================
!       VHV (R8) = lower triangle of the Hamiltonian (OUTPUT)
!        HR (R8) = buffer to read the Hamiltonian (SCRATCH)
!      ISEL (I4) = indix of the selected determinants (INPUT)
!                 isel(i) = 0 if det. |i> not selected
!                 isel(i) = j if det. |i> is the jth in the selected space
!      ICOL (I4) = position of the non-zero elements of the Hamiltonian
!                  matrix
!         N (I4) = size of the full determinant space (INPUT)
!      HBUF (I4) = dimension of the buffers in H (INPUT)
!    ELFILE (I4) = unit number of the non relativistic part of H (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   use mconst
#ifdef _MYDEBUG_
   use mdebug
#endif
   implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
   integer, intent(in) :: n,hbuf,elfile
   integer, intent(in), dimension (:) :: isel
   integer, intent(inout), dimension (hbuf) :: icol
   integer, intent(inout), dimension (:) :: elstar
   double precision, intent(out), dimension (:) :: vhv
   double precision, intent(inout), dimension (hbuf) :: hr
!
!     ----- local variables -----
!
   integer :: i,j,k,i1,i2,j1,j2,elbuf,nbj,nb2,ib1,ib2
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
   call subent('HRED')
#endif
   vhv = ZERO
   rewind(elfile)
   read(elfile) hr,icol
   i1 = elstar(1)
   elbuf = 0
   do i=1,n
      i2 = elstar(i+1)
      nbj = i2-i1
 10   continue
      ib1 = elbuf+1
      ib2 = elbuf+nbj
      if(ib2.le.hbuf) then
         if(isel(i).ne.0) then
            do k=ib1,ib2
               j = icol(k)
               if(isel(j).ne.0) then
                  j1 = max(isel(i),isel(j))
                  j2 = min(isel(i),isel(j))
                  vhv(j1*(j1-1)/2+j2) = hr(k)
               endif
            enddo
         endif
         if(ib2.eq.hbuf) then
            nb2 = ib2-(ib1-1)
            read(elfile) hr,icol
            nbj = nbj-nb2
            elbuf = 0
            i1 = elstar(i+1)
            cycle
         endif
      else
         ib2 = hbuf
         nb2 = ib2-(ib1-1)
         if(isel(i).ne.0) then
            do k=ib1,ib2
               j = icol(k)
               if(isel(j).ne.0) then
                  j1 = max(isel(i),isel(j))
                  j2 = min(isel(i),isel(j))
                  vhv(j1*(j1-1)/2+j2) = hr(k)
               endif
            enddo
         endif
         read(elfile) hr,icol
         nbj = nbj-nb2
         elbuf = 0
         goto 10
      endif
      elbuf = elbuf+nbj
      i1 = elstar(i+1)
   enddo
#ifdef _MYDEBUG_
   call subret('HRED',2)
#endif
   return
end subroutine hred
