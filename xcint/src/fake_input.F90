module fake_input

   use keyword
   use character_processing
   use dirac_cfg
   use interface_input
   use interface_functional_read
   use xc_derv

   implicit none

   public read_input

   private

contains

   subroutine read_input()

!     --------------------------------------------------------------------------
      character(kw_length) :: word
      character(kw_length) :: kw_section
      character(80)        :: line
      real(8)              :: hfx_out, mu_out, beta_out
!     --------------------------------------------------------------------------

#ifdef VAR_DALTON
      line = 'lda       '// &
             '          '// &
             '          '// &
             '          '// &
             '          '// &
             '          '// &
             '          '// &
             '          '
      call parse_functional(line, xc_fun, hfx_out, mu_out, beta_out, .true.)
      fun_is_automatic = .true.
#else

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

1     rewind unit_in
      do while (.true.)
         read(unit_in, '(a80)', end=2) line

         if (line(1:7) == '.DFTAUT') then
            read(unit_in, '(a80)') line
            call parse_functional(line, xc_fun, hfx_out, mu_out, beta_out, .true.)
         end if
      end do

2     close(unit_in)
#endif

   end subroutine

end module
