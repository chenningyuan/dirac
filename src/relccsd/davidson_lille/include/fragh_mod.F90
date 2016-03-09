module fragh_mod
   interface fragh
      subroutine fragh(e,mvec,mvmax,nseg,nfrag)
         integer, intent(in) :: mvec
         integer, intent(inout) :: mvmax
         integer, intent(out) :: nfrag
         integer, intent(out), dimension (:) :: nseg
         double precision, intent(in), dimension (:) :: e
      end subroutine fragh
   end interface fragh
end module fragh_mod
