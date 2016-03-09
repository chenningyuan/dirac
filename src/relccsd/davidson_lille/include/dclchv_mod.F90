module dclchv_mod
   interface dclchv
      subroutine dclchv(hv,vr,n,nvec,hr,icol,hbuf,ifile,elstar,work)
         use kindvar
         integer, intent(in) :: n,nvec,hbuf,ifile
         integer, intent(inout), dimension (hbuf) :: icol
         integer, intent(inout), dimension (:) :: elstar
         double precision, intent(in), dimension (:,:) :: vr
         double precision, intent(out), dimension (:,:) :: hv,work
         double precision, intent(inout), dimension (hbuf) :: hr
      end subroutine dclchv
   end interface dclchv
end module dclchv_mod
