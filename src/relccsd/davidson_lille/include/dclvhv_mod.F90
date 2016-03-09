module dclvhv_mod
   interface dclvhv
      subroutine dclvhv(vhv,hv,vr,ndet,nvec)
         integer :: ndet,nvec
         double precision, intent(inout), dimension (:) :: vhv
         double precision, intent(in), dimension (:,:) :: hv,vr
      end subroutine dclvhv
   end interface dclvhv
end module dclvhv_mod
