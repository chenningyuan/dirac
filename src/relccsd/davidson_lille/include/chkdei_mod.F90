module chkdei_mod
   interface chkdei
      function chkdei(c,e,m,n) result(value)
         integer, intent(in) :: m,n
         double precision, intent(inout), dimension (:) :: e
         double precision, intent(in), dimension (:,:) :: c
         logical :: value
         double precision, dimension(:), allocatable :: tmp
      end function chkdei
   end interface chkdei
end module chkdei_mod
