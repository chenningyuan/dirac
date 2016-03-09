module savvls_mod
   interface savvls
      subroutine savvls(v,eps,ndet,mvec,iunit,labsym)
         use kindvar
         character (len=*), optional, dimension (:) :: labsym
         integer, intent(in) :: ndet,mvec,iunit
         double precision, intent(in), dimension (:,:) :: v
         double precision, intent(in), dimension (:) :: eps
      end subroutine savvls
   end interface savvls
end module savvls_mod
