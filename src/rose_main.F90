   program rose

        use rose_utils
        use rose_input_data
        use rose_global

        implicit none

        logical :: input_found=.false.
 
        write(luout,*)
        write(luout,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(luout,*) "@                                                @"
        write(luout,*) "@       Reduction of Orbital Space Extent        @"
        write(luout,*) "@                     (ROSE)                     @"
        write(luout,*) "@                                                @"
        write(luout,*) "@ By Bruno Senjean, Souloke Sen, Michal Repisky, @"
        write(luout,*) "@        Gerald Knizia and Lucas Visscher.       @"
        write(luout,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(luout,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           READ INPUT FILE AND OPTIONS         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call get_free_fileunit(luinp)
        call read_menu_input('INPUT_GENIBO', luinp, 'ALL', input_found)
        if (.not.input_found) error stop 'INPUT_GENIBO not found'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           CHECK IF OPTIONS ARE CORRECT         !
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

        call generate_ibo

!!!!!!!!!!!!!!!!!!!!
!       END        !
!!!!!!!!!!!!!!!!!!!!
        write(luout,*)
        write(luout,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(luout,*) "@                                                @"
        write(luout,*) "@       Reduction of Orbital Space Extent        @"
        write(luout,*) "@                                                @"
        write(luout,*) "@                      END                       @"
        write(luout,*) "@                                                @"
        write(luout,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

   end
