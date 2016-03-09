module cschmi_mod
   interface cschmi
      subroutine cschmi(pr,pi,n,ivec,ni,stol)
         integer, intent(in) :: n,ivec
         integer, intent(inout) :: ni
         double precision, intent(in) :: stol
         double precision, intent(inout), dimension (:,:) :: pr,pi
      end subroutine cschmi
   end interface cschmi
end module cschmi_mod
