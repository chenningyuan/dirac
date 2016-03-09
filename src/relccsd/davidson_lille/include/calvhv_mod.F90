module calvhv_mod
   interface calvhv
      subroutine calvhv(vhv,hvr,hvi,vr,vi,ndet,nvec)
         use kindvar
         integer, intent(in) :: ndet,nvec
         complex (kind=K_C8), intent(out), dimension (:) :: vhv
         double precision, intent(in), dimension (:,:) :: hvr,hvi,vr,vi
      end subroutine calvhv
   end interface calvhv
end module calvhv_mod
