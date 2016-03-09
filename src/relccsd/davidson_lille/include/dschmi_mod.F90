module dschmi_mod
   interface dschmi
      subroutine dschmi(pr,n,ivec,ni,stol)
         use kindvar
         integer, intent(in) :: n,ivec
         integer, intent(inout) :: ni
         double precision, intent(in) :: stol
         double precision, intent(inout), dimension (:,:) :: pr
      end subroutine dschmi
   end interface dschmi
end module dschmi_mod
