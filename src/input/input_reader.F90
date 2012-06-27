
! primitive input reader
! (c) Radovan Bast and Stefan Knecht
! licensed under the GNU Lesser General Public License

module input_reader

   use keyword
   use character_processing
   use interface_input
   use input_reader_sections

   implicit none

   public read_input
   public move_to_next_star

   private

contains

   subroutine read_input()

!    radovan: - this routine should not include any common blocks
!               this should be done in called subroutines
!             - please keep ifndef traps inside subroutines for better readability
!               of this subroutine

!     --------------------------------------------------------------------------
      character(kw_length) :: word
      character(kw_length) :: kw_section
!     --------------------------------------------------------------------------

!     this is to catch keywords that appear before some section starts
      kw_section = '       '

      open(unit_in, file='DIRAC.INP')
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

!                    case ('**METHO', '**WAVE ')
!                       call read_input_method(word, kw_section)

                     case ('*DFT   ')
                        call read_input_dft(word, kw_section)

                     case ('**GRID ')
                        call read_input_grid(word, kw_section)

!                    *VISUAL is for backward compatibility
!                    *VISUAL used to be under **ANALYZE
                     case ('**VISUA', '*VISUAL')
                        call read_input_visual(word, kw_section)

                     case ('*OPENRS')
                        call read_input_openrsp(word, kw_section)

                     case ('*X2C')
                        call read_input_x2c(word, kw_section)

                     case ('**RELCC')
                        call read_input_relcc(word, kw_section)

                     case ('*CCENER')
                        call read_input_ccener(word, kw_section)

                     case ('*CCFOPR')
                        call read_input_ccfopr(word, kw_section)

                     case ('*CCSORT')
                        call read_input_ccsort(word, kw_section)

                     case ('*CCFSPC')
                        call read_input_ccfspc(word, kw_section)

                     case ('*CCIH  ')
                        call read_input_ccih(word, kw_section)
                     
                     case ('*PCM   ')
                        call read_input_pcm(word, kw_section)

                     case ('*END OF')
                        go to 1

!                    case ('       ')
!                       call quit('keyword '//word//' without section')

                     case default
!                       activate as soon as everything is merged to here
!                       call quit('section '//kw_section//' not recognized')

                  end select

            end select
         end if

      end do

1     close(unit_in)

   end subroutine

   subroutine move_to_next_star(word_io)

!     --------------------------------------------------------------------------
      character(kw_length), intent(inout) :: word_io
!     --------------------------------------------------------------------------
      character(kw_length)                :: word
!     --------------------------------------------------------------------------

      do while (.true.)
         read(unit_in, '(a7)') word
         if (word(1:1) == '*') then
            word_io = word
            return
         end if
      end do

   end subroutine

end module
