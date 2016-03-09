module fiflsh_mod
   interface fiflsh
      subroutine fiflsh(nft)
         use kindvar
         integer, intent(in) :: nft
      end subroutine fiflsh
   end interface fiflsh
end module fiflsh_mod
