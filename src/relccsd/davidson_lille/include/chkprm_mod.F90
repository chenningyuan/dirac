module chkprm_mod
   interface chkprm
      function chkprm(c,e,m,n) result(value)
         use kindvar
         integer, intent(in) :: m,n
         double precision, intent(inout), dimension (:) :: e
         complex (kind=K_C8), intent(inout), dimension (:,:) :: c
         logical :: value
         complex (kind=K_C8),dimension(:), allocatable :: tmp
      end function chkprm
   end interface chkprm
end module chkprm_mod
