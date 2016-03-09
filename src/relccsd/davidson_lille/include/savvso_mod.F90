module savvso_mod
   interface savvso
      subroutine savvso(vr,vi,eps,omega,ndet,mvec,iunit,labsym)
         use kindvar
         character (len=*), optional, dimension (:) :: labsym
         integer, intent(in) :: ndet,mvec,iunit
         double precision, intent(in), dimension (:,:) :: vr,vi
         double precision, intent(in), dimension (:) :: eps,omega
      end subroutine savvso
   end interface savvso
end module savvso_mod
