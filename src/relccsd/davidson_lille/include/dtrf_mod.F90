module dtrf_mod
   interface dtrf
      subroutine dtrf(cv,mv,pr,np,ndet,nvec,wr,wkdim)
         integer, intent(in) :: mv,np,ndet,nvec,wkdim
         double precision, intent(in), dimension (:,:) :: cv
         double precision, intent(inout), dimension (:,:) :: pr,wr
      end subroutine dtrf
   end interface dtrf
end module dtrf_mod
