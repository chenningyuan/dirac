module hred_mod
   interface hred
      subroutine hred(vhv,hr,isel,icol,n,hbuf,elfile,elstar)
         integer, intent(in) :: n,hbuf,elfile
         integer, intent(in), dimension (:) :: isel
         integer, intent(inout), dimension (hbuf) :: icol
         integer, intent(inout), dimension (:) :: elstar
         double precision, intent(out), dimension (:) :: vhv
         double precision, intent(inout), dimension (hbuf) :: hr
      end subroutine hred
   end interface hred
end module hred_mod
