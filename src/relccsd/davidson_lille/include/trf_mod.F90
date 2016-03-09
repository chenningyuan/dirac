module trf_mod
   interface trf
      subroutine trf(cv,vr,vi,ndet,nvec,win,wout,wkdim)
         use kindvar
         integer, intent(in) :: ndet,nvec,wkdim
         double precision, intent(inout), dimension (:,:) :: vr,vi
         complex(kind=K_C8), intent(in), dimension (:,:) :: cv
         complex(kind=K_C8), intent(inout), dimension (:,:) :: win,wout
      end subroutine trf
   end interface trf
end module trf_mod
