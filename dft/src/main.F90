program main

   implicit none

   call read_input()

contains

   subroutine read_input()

      use keyword
      use character_processing
      use dirac_cfg
      use interface_input

!     --------------------------------------------------------------------------
      character(kw_length) :: word
      character(kw_length) :: kw_section
!     --------------------------------------------------------------------------

!     this is to catch keywords that appear before some section starts
      kw_section = '       '

      open(unit_in, file='xc_input')
      rewind unit_in

      do while (.true.)

         read(unit_in, '(a7)', end=1) word

         if (word == '       ') then
!           blank line
         else
            select case (word(1:1))
               case ('!', '#')
!                 comment

               case ('*')
!                 section
                  kw_section = uppercase(word)

               case default
!                 keyword
                  select case (kw_section)
                     case ('*DFT   ')
                        call read_input_dft(word, kw_section)
                     case ('*END OF')
                        go to 1
                     case default
                  end select
            end select
         end if

      end do

1     close(unit_in)

   end subroutine

end program
