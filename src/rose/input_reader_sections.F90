   subroutine read_input_rose(word, kw_section)

      use rose_input_data
      use input_reader
      use jagged_array

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      integer           :: i_frag, n_orbitals
      integer, allocatable :: orbital_indices(:)

      call reset_available_kw_list()

      if (kw_matches(word, '.VERSION')) then
         call kw_read(word, version)
      end if

      if (kw_matches(word, '.CHARGE')) then
         call kw_read(word, charge)
      end if

      if (kw_matches(word, '.EXPONENT')) then
         call kw_read(word, expnt)
      end if

      if (kw_matches(word, '.UNRESTRICTED')) then
         restricted = .false.
      end if

      if (kw_matches(word, '.SPINORS')) then
         spatial_orbitals = .false.
      end if

      if (kw_matches(word, '.NFRAGMENTS')) then
         call kw_read(word, n_fragments)
         ! The following arrays can be tweaked via input, so allocate them
         allocate(frag_ncor(n_fragments))
         allocate(frag_nval(n_fragments))
         allocate(frag_bias(n_fragments))
         frag_ncor = -1  ! signal that this is not set by the user
         frag_nval = -1  ! signal that this is not set by the user
         frag_bias = 0.0 ! by default we do not bias (except the first fragment, but this is done later)
      end if

      if (kw_matches(word, '.FILE_FORMAT')) then
         call kw_read(word, file_for_MOs)
         file_for_molfrag = file_for_MOs
      end if

      if (kw_matches(word, '.FRAG_FORMAT')) then
         call kw_read(word, file_for_molfrag)
      end if

      if (kw_matches(word, '.TEST  ')) then
         test = .true.
      end if

      if (kw_matches(word, '.INCLUDE_CORE')) then
         include_core = .true.
      end if

      if (kw_matches(word, '.ADDITIONAL_VIRTUALS_CUTOFF')) then
         call kw_read(word, additional_virtuals_cutoff)
      end if

      if (kw_matches(word, '.FRAG_THRESHOLD')) then
         call kw_read(word, fragment_cutoff)
      end if

      if (kw_matches(word, '.FRAG_VALENCE')) then
         call kw_read(word, i_frag)
         if (i_frag > 0 .and. i_frag <= n_fragments) then
            call kw_read(word, frag_nval(i_frag))
         else
            print*," You want to set the valence for fragment ",i_frag," but the number of fragments is currently ",n_fragments
            stop "invalid input"
         end if
      end if

      if (kw_matches(word, '.FRAG_CORE')) then
         call kw_read(word, i_frag)
         if (i_frag > 0 .and. i_frag <= n_fragments) then
            call kw_read(word, frag_ncor(i_frag))
         else
            print*," You want to set the core for fragment ",i_frag," but the number of fragments is currently ",n_fragments
            stop "invalid input"
         end if
      end if

      if (kw_matches(word, '.FRAG_BIAS')) then
         call kw_read(word, i_frag)
         if (i_frag > 0 .and. i_frag <= n_fragments) then
            read(get_file_unit(), *) frag_bias(i_frag)
         else
            print*," You want to bias ",i_frag," but the number of fragments is currently ",n_fragments
            stop "invalid input"
         end if
      end if

      if (kw_matches(word, '.FRAG_AVAS')) then
         call kw_read(word, i_frag)
         if (i_frag > 0 .and. i_frag <= n_fragments) then
            ! if not yet allocated array_type will be blank and we allocate 
            if (avas_orbitals%array_type == ' ') call NewJA(avas_orbitals,n_fragments,'i')
            read(get_file_unit(), *) n_orbitals
            allocate(orbital_indices(n_orbitals))
            read(get_file_unit(), *) orbital_indices
            call JASetColumn(avas_orbitals,i_frag,orbital_indices)
         else
            print*," You want to get avas orbitals for ",i_frag," but the number of fragments is currently ",n_fragments
            stop "invalid input"
         end if
      end if

      if (kw_matches(word, '.AVAS_THRESHOLD')) then
         call kw_read(word, avas_threshold)
      end if

      call check_whether_kw_found(word, kw_section)

  end subroutine

  subroutine rose_check_input

      use rose_input_data

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !           CHECK IF FILE TYPE IS SUPPORTED      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (((trim(file_for_MOs) .eq. "adf") .and. (trim(file_for_molfrag) .ne. "adf")) .or. &
            ((trim(file_for_MOs) .ne. "adf") .and. (trim(file_for_molfrag) .eq. "adf"))) then
           error stop 'ADF has to be used for both molecule and fragments, because ADF works with Slater type orbitals'
        endif

        if (trim(file_for_MOs) .ne. "dirac" .and. &
            trim(file_for_MOs) .ne. "dfcoef" .and. &
            trim(file_for_MOs) .ne. "h5" .and. &
            trim(file_for_MOs) .ne. "fchk" .and. &
            trim(file_for_MOs) .ne. "psi4" .and. &
            trim(file_for_MOs) .ne. "pyscf" .and. &
            trim(file_for_MOs) .ne. "ams" .and. &
            trim(file_for_MOs) .ne. "adf") then
            print*,'Wrong input for the file containing the MOs:',trim(file_for_MOs)
            error stop
        elseif (trim(file_for_MOs) .eq. "psi4" .or. file_for_MOs .eq. "pyscf") then
           print*," Psi4 and pyscf are now interfaced via .h5 files, please update your input"
           stop 'old psi4 and pyscf interfaces are no longer supported'
        elseif (trim(file_for_MOs) .eq. "dfcoef") then
           print*," Legacy format dfcoef no longer supported, please update your DIRAC installation"
           stop 'dfcoef no longer supported'
        endif
        if (trim(file_for_molfrag) .ne. "dirac" .and. &
            trim(file_for_molfrag) .ne. "dfcoef" .and. &
            trim(file_for_molfrag) .ne. "h5" .and. &
            trim(file_for_molfrag) .ne. "fchk" .and. &
            trim(file_for_molfrag) .ne. "psi4" .and. &
            trim(file_for_molfrag) .ne. "pyscf" .and. &
            trim(file_for_molfrag) .ne. "ams" .and. &
            trim(file_for_molfrag) .ne. "adf") then
            print*,'Wrong input for the file containing the MOs of molecular fragments:',trim(file_for_molfrag)
            error stop
        elseif (trim(file_for_molfrag) .eq. "psi4" .or. file_for_MOs .eq. "pyscf") then
           print*," Psi4 and pyscf are now interfaced via .h5 files, please update your input"
           error stop 'old psi4 and pyscf interfaces are no longer supported'
        elseif (trim(file_for_molfrag) .eq. "dfcoef") then
           print*," Legacy format dfcoef no longer supported, please update your DIRAC installation"
           stop 'dfcoef no longer supported'
        endif

   end subroutine
