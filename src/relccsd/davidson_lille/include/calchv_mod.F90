module calchv_mod
   interface calchv
      subroutine calchv(hvr,hvi,vr,vi,n,nvec,hr,hi,icol,hbuf,delta,d,mdress,  &
                        heff,work,sfile,elfile,elstar,sostar)
         use kindvar
         integer, intent(in) :: n,nvec,hbuf,mdress,sfile,elfile
         logical, intent(in) :: heff
         integer, intent(inout), dimension (hbuf) :: icol
         integer, intent(inout), dimension (:) :: elstar,sostar
         double precision, intent(in), dimension (:) :: delta
         double precision, intent(in), dimension (:,:) :: d
         double precision, intent(in), dimension (:,:) :: vr,vi
         double precision, intent(inout), dimension (hbuf) :: hr,hi
         double precision, intent(inout), dimension (:,:) :: work,hvr,hvi
      end subroutine calchv
   end interface calchv
end module calchv_mod
