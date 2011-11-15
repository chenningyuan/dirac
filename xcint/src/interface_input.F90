module interface_input

   use keyword
   use character_processing

   implicit none

   public read_input_dft

   private

contains

   subroutine read_input_dft(word, kw_section)

      use dft_cfg
      use sigma_prefactor_setting

!     ---------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     ---------------------------------------------------------------------------
      character(80)                    :: line
      character(4)                     :: first_4
      integer                          :: i, n
!     ---------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.TINYDE')) then
         call kw_read(word, dft_cfg_tinydens)
      end if

      if (kw_matches(word, '.SCREEN')) then
         call kw_read(word, dft_cfg_screening)
      end if

      if (kw_matches(word, '.GRAC  ')) then
!        gradient regulated asymptotic correction
         dft_cfg_grac = .true.
         read(unit_in, '(a4)') first_4
         if      (lowercase(first_4) == 'lb94') then
            dft_cfg_asymptote_is_lb94 = .true.
         else if (lowercase(first_4) == 'lbal') then
            dft_cfg_asymptote_is_lbalpha = .true.
         else
            print *, 'keyword following .GRAC not recognized'
            stop
         end if
         call kw_read(word,               &
                      dft_cfg_grac_alpha, &
                      dft_cfg_grac_beta,  &
                      dft_cfg_ac_ip,      &
                      dft_cfg_ac_threshold)
      end if

      if (kw_matches(word, '.SAOP  ')) then
!        statistical averaging of (model) orbital potentials
         dft_cfg_saop = .true.
         read(unit_in, '(a4)') first_4
         if      (lowercase(first_4) == 'lb94') then
            dft_cfg_asymptote_is_lb94 = .true.
         else if (lowercase(first_4) == 'lbal') then
            dft_cfg_asymptote_is_lbalpha = .true.
         else
            print *, 'keyword following .SAOP not recognized'
            stop
         end if
      end if

      if (kw_matches(word, '.SAOP! ')) then
!        "original" SAOP as in JCP 112, 1344 (2000).
         dft_cfg_saop                    = .true.
         dft_cfg_saop_with_response_part = .true.
         dft_cfg_asymptote_is_lbalpha    = .true.
         dft_cfg_alda_hs                 = .true.
         dft_cfg_alda_ha                 = .true.
         call setaldahs(.true.) !activate alda for the h+ part
         call setaldaha(.true.) !activate alda for the h- part
      end if

      if (kw_matches(word, '.NOSDFT')) then
         dft_cfg_no_sdft = .true.
      end if

      if (kw_matches(word, '.COLLIN')) then
         dft_cfg_sdft_collinear = .true.
      end if

      if (kw_matches(word, '.BETASI')) then
         sigma_prefactor = -1.0d0
      end if

      if (kw_matches(word, '.ALDA  ')) then
         dft_cfg_alda_hs = .true.
         dft_cfg_alda_ha = .true.
         call setaldahs(.true.) !activate alda for the h+ part
         call setaldaha(.true.) !activate alda for the h- part
      end if

      if (kw_matches(word, '.XALDA ')) then
!        use        (1-f)*S + VWN + f*HFx
!        instead of       S + VWN
         dft_cfg_xalda   = .true.
         dft_cfg_alda_hs = .true.
         dft_cfg_alda_ha = .true.
         call setxalda(.true.)
         call setaldahs(.true.) !activate alda for the h+ part
         call setaldaha(.true.) !activate alda for the h- part
      end if

      if (kw_matches(word, '.ALDA+ ')) then
!        use alda only for the      hermitian part
         dft_cfg_alda_hs = .true.
         call setaldahs(.true.) !activate alda for the h+ part
      end if

      if (kw_matches(word, '.ALDA- ')) then
!        use alda only for the anti-hermitian part
         dft_cfg_alda_ha = .true.
         call setaldaha(.true.) !activate alda for the h- part
      end if

      if (kw_matches(word, '.XALDA+')) then
!        use xalda only for the      hermitian part
         dft_cfg_xalda   = .true.
         dft_cfg_alda_hs = .true.
         call setxalda(.true.)
         call setaldahs(.true.) !activate alda for the h+ part
      end if

      if (kw_matches(word, '.XALDA-')) then
!        use xalda only for the anti-hermitian part
         dft_cfg_xalda   = .true.
         dft_cfg_alda_ha = .true.
         call setxalda(.true.)
         call setaldaha(.true.) !activate alda for the h- part
      end if

      if (kw_matches(word, '.GAUTSC')) then
         call set_lscale_dft_gaunt()
      end if

!radovan: not documented on wiki, code not perfectly tested
      if (kw_matches(word, '.BLOCKE')) then
         dft_cfg_blocked = .true.
      end if

      if (kw_matches(word, '.POINTW')) then
         dft_cfg_pointwise = .true.
      end if

!radovan: not documented on wiki, code not perfectly tested
      if (kw_matches(word, '.OVLDIA')) then
         dft_cfg_overlap_diagnostic = .true.
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

end module
