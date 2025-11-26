!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end
#if defined(HAS_PELIB)
module pelib_dirac_interface
    use pelib_precision, only : rp, ip, lp

    implicit none

    private
     
    public :: use_pelib, pelib_ifc_gspol
    public :: pelib_ifc_activate, pelib_ifc_deactivate
    public :: pelib_ifc_init, pelib_ifc_finalize
    public :: pelib_ifc_input_reader
    public :: pelib_ifc_localfield, pelib_ifc_do_lf
    public :: pelib_ifc_fock
    public :: pelib_ifc_energy
    public :: pelib_ifc_response
#if defined(VAR_MPI)
    public :: pelib_ifc_slave
#endif

    contains

    logical function use_pelib()
        use pelib, only: pelib_enabled
        implicit none
        if (pelib_enabled) then
            use_pelib = .true.
        else
            use_pelib = .false.
        end if
    end function use_pelib

    subroutine pelib_ifc_activate()
        use pelib, only: pelib_enabled
        call qenter('pelib_ifc_activate')
        if (use_pelib()) call quit('PElib already active')
        pelib_enabled = .true.
        call qexit('pelib_ifc_activate')
    end subroutine pelib_ifc_activate

    subroutine pelib_ifc_deactivate()
        use pelib, only: pelib_enabled
        call qenter('pelib_ifc_deactivate')
        if (.not. use_pelib()) call quit('PElib already deactivated')
        pelib_enabled = .false.
        call qexit('pelib_ifc_deactivate')
    end subroutine pelib_ifc_deactivate

    logical function pelib_ifc_gspol()
        use pelib_options, only: pelib_gspol
        implicit none
        if (pelib_gspol) then
            pelib_ifc_gspol = .true.
        else
            pelib_ifc_gspol = .false.
        end if
    end function pelib_ifc_gspol

    logical function pelib_ifc_do_lf()
        use pelib_options, only: pelib_lf
        implicit none
        if (pelib_lf) then
            pelib_ifc_do_lf = .true.
        else
            pelib_ifc_do_lf = .false.
        end if
    end function pelib_ifc_do_lf

    subroutine pelib_ifc_input_reader(word, use_pde)
        !use pelib_precision
        use pelib_options
        use pelib_constants
        use pelib_utils, only: chcase
        use pelib_cavity_generators, only: ntsatm
#include "priunit.h"
        character(len=*), intent(inout) :: word
        character(len=80) :: option
        character(len=2) :: auoraa
        integer :: i, j
        real(rp), dimension(2) :: temp
        logical, intent(inout) :: use_pde

        call qenter('pelib_ifc_input_reader')
  
        potfile = 'POTENTIAL.INP'
        h5pdefile = 'standard.h5'
        pelib_skipmul = .true.
  
        do
            read(lucmd, '(a80)') option
            call chcase(option)
  
            ! Read potential (optionally from potfile)
            if (trim(option(2:7)) == 'POTENT') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, '(a240)') potfile
                end if
            ! direct solver for induced moments
            else if (trim(option(2:7)) == 'DIRECT') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, '(a80)') option
                    call chcase(option)
                    if (trim(option(1:7)) == 'NOCHOL') then
                        chol = .false.
                    end if
                end if
                pelib_iter = .false.
            ! iterative solver for induced moments (default)
            else if (trim(option(2:7)) == 'ITERAT') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) thriter
                end if
                pelib_iter = .true.
            else if (trim(option(2:7)) == 'DIIS T') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) thrdiis
                end if
                thriter = thrdiis
            ! use reduced threshold in iterative induced moments solver
            else if (trim(option(2:7)) == 'REDTHR') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) redlvl
                end if
                pelib_redthr = .true.
            ! handling sites near quantum-classical border
            else if (trim(option(2:7)) == 'BORDER') then
               read(lucmd, '(a80)') option
               backspace(lucmd)
               if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                 & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                   read(lucmd, '(a)', advance='no') border_type
                   backspace(lucmd)
                   call chcase(border_type)
                   if ((trim(border_type) /= 'REMOVE') .and.&
                     & (trim(border_type) /= 'REDIST') .and.&
                     & (trim(border_type) /= 'CHGRED')) then
                       call quit('unknown handling of border sites')
                   else if (trim(border_type) == 'REMOVE') then
                       read(lucmd, *) border_type, Rmin, auoraa
                   else if (trim(border_type) == 'REDIST') then
                       read(lucmd, *) border_type, redist_order, Rmin, auoraa, nredist
                       if ((nredist > 3) .or. (nredist < 1)) then
                           call quit('parameters can only be distributed to&
                                & minimum one site and maximum three sites')
                       end if
                   else if (trim(border_type) == 'CHGRED') then
                       read(lucmd, *) border_type, Rmin, auoraa
                   else
                       call quit('unrecognized input in .BORDER option')
                   end if
                   call chcase(auoraa)
                   if (trim(auoraa) == 'AA') Rmin = Rmin * aa2bohr
               end if
               pelib_border = .true.
          ! damp electric field from induced multipoles
          else if (trim(option(2:7)) == 'DAMP I') then
              read(lucmd, '(a80)') option
              backspace(lucmd)
              if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                  read(lucmd, *) ind_damp
              end if
              pelib_ind_damp = .true.
            ! damp electric field from permanent multipoles
            else if (trim(option(2:7)) == 'DAMP M') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) mul_damp
                end if
                pelib_mul_damp = .true.
            ! damp electric field from core region
            else if (trim(option(2:7)) == 'DAMP C') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) core_damp
                    ! attempt to optionally read a custom specification of the polarizabilities
                    read(lucmd, '(a80)') option
                    backspace(lucmd)
                    if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                       & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                        read(lucmd, *) j
                        allocate(core_alphas(6,j))
                        core_alphas = 0.0_rp
                        do i = 1, j
                            read(lucmd, *) core_alphas(1,i)
                            core_alphas(4,i) = core_alphas(1,i)
                            core_alphas(6,i) = core_alphas(1,i)
                        enddo
                    end if
                end if
                pelib_core_damp = .true.
            ! damp electric fields using AMOEABA-style Thole damping
            else if (trim(option(2:7)) == 'DAMP A') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) amoeba_damp
                end if
                pelib_amoeba_damp = .true.
            ! the old deprecated DAMP option
            else if (trim(option(2:7)) == 'DAMP') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) ind_damp
                end if
                pelib_ind_damp = .true.
                write(luout, *) 'WARNING: the .DAMP option is deprecated, please use .DAMP INDUCED'
               ! write(luerr, *) 'WARNING: the .DAMP option is deprecated, please use .DAMP INDUCED'
            ! neglect dynamic response from environment
            else if (trim(option(2:7)) == 'GSPOL') then
                pelib_gspol = .true.
            ! neglect many-body interactions
            else if (trim(option(2:7)) == 'NOMB') then
                pelib_nomb = .true.
            ! Use existing files for restart
            else if (trim(option(2:7)) == 'RESTAR') then
                pelib_restart = .true.
            ! calculate intermolecular two-electron integrals
            else if (trim(option(2:7)) == 'TWOINT') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, '(a240)') h5pdefile
                end if
                pelib_twoints = .true.
            ! activate calculation of london contributions for PDE
            else if (trim(option(2:7)) == 'TWOLAO') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, '(a240)') h5pdefile
                end if
                pelib_lao_twoints = .true.
            ! save density matrix
            else if (trim(option(2:7)) == 'SAVE D') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, '(a240)') h5pdefile
                end if
                pelib_savden = .true.
            ! disable exchange repulsion
            else if (trim(option(2:7)) == 'NO REP') then
                pelib_repuls = .false.
            ! electrostatics and exchange repulsion from fragment densities
            else if (trim(option(2:7)) == 'PDE') then
#if !defined(VAR_PDE)
                call quit('Fatal error: .PDE not included in this build of Dirac')
#endif
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, '(a240)') h5pdefile
                end if
                pelib_fd = .true.
                use_pde = .true. 
            ! skip electrostatics from fragment densities
            else if (trim(option(2:7)) == 'NO FD') then
                pelib_fdes = .false.
            ! activate legacy PE field evaluation code (OLDFIELD)
            else if (trim(option(2:7)) == 'OLDFIE') then
                pelib_old_field = .true.
            ! activate Fast Multipole Method for evaluation of multipolar fields
            else if (trim(option(2:7)) == 'FMM') then
                pelib_fmm = .true.
            else if (trim(option(2:7)) == 'THETA') then
                read(lucmd, *) fmm_theta
            else if (trim(option(2:7)) == 'NCRIT') then
                read(lucmd, *) fmm_ncrit
            else if (trim(option(2:7)) == 'EXPANS') then
                read(lucmd, *) fmm_expansion_order
            ! options for (approximate) qmmm interaction 
            else if (trim(option(2:7)) == 'INT_TH') then ! INT_THETA
                read(lucmd, *) qmmm_interaction_theta
            else if (trim(option(2:7)) == 'INT_SC') then ! INT_SCHEME
                read(lucmd, *) qmmm_interaction_scheme
            else if (trim(option(2:7)) == 'INT_RC') then ! INT_RCUT
                read(lucmd, *) qmmm_interaction_rcut
            else if (trim(option(2:7)) == 'INT_RF') then ! INT_RFIT
                read(lucmd, *) qmmm_interaction_rfit
            else if (trim(option(2:7)) == 'INT_OR') then ! INT_ORDER
                read(lucmd, *) qmmm_interaction_order
            else if (trim(option(2:7)) == 'ESPF_S') then ! ESPF_SOLVER
                read(lucmd, *) ESPF_solver
            else if (trim(option(2:7)) == 'ESPF_G') then ! ESPF_GRID_TYPE
                read(lucmd, *) ESPF_grid_type
            else if (trim(option(2:7)) == 'SVDTOL') then ! SVDTOL
                read(lucmd, *) svdtol
            else if (trim(option(2:7)) == 'NF_INC') then ! NF_INCORE
                incore_nearfield = .true.
            else if (trim(option(2:7)) == 'NONF_I') then ! NONF_INCORE
                incore_nearfield = .false.
            ! request calculation of effective dipole integrals
            else if (trim(option(2:7)) == 'EEF') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) ncrds
                    allocate(crds(3,ncrds))
                    do i = 1, ncrds
                        read(lucmd, *) (crds(j,i), j = 1, 3)
                    end do
                end if
                pelib_lf = .true.
            ! provide LJ parameters for the QM region
            else if (trim(option(2:7)) == 'LJ') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) qmLJsites
                    allocate(qmLJs(2,qmLJsites))
                    do i = 1, qmLJsites
                        read(lucmd, *) (temp(j), j = 1, 2)
                        qmLJs(1,i) = temp(1) * 2.0_rp
                        qmLJs(2,i) = temp(2)
                    end do
                    lvdw = .true.
                 end if
            ! skip QM calculations, i.e. go directly into PE library
            else if (trim(option(2:7)) == 'SKIPQM') then
                pelib_skipqm = .true.
            ! Write cube files
            else if (trim(option(2:7)) == 'CUBE') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                  & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    do
                        read(lucmd, '(a80)') option
                        call chcase(option)
                        if (trim(option(1:7)) == 'COARSE') then
                            xgrid = 3
                            ygrid = 3
                            zgrid = 3
                        else if (trim(option(1:7)) == 'MEDIUM') then
                            xgrid = 6
                            ygrid = 6
                            zgrid = 6
                        else if (trim(option(1:7)) == 'FINE') then
                            xgrid = 12
                            ygrid = 12
                            zgrid = 12
                        else if (trim(option(1:7)) == 'GRID') then
                            read(lucmd, *) xsize, xgrid, ysize, ygrid, zsize, zgrid
                        else if (trim(option(1:7)) == 'FIELD') then
                            cube_field = .true.
                        else if (option(1:1) == '.' .or. option(1:1) == '*') then
                            backspace(lucmd)
                            exit
                        else if (option(1:1) == '!' .or. option(1:1) == '#') then
                            cycle
                        else
                            call quit('unknown input under .CUBE option')
                        end if
                    end do
                end if
                pelib_cube = .true.
            ! evaluate molecular electrostatic potential
            else if (trim(option(2:7)) == 'MEP') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    do
                        read(lucmd, '(a80)') option
                        call chcase(option)
                        if (trim(option(1:7)) == 'COARSE') then
                            xgrid = 3
                            ygrid = 3
                            zgrid = 3
                        else if (trim(option(1:7)) == 'MEDIUM') then
                            xgrid = 6
                            ygrid = 6
                            zgrid = 6
                        else if (trim(option(1:7)) == 'FINE') then
                            xgrid = 12
                            ygrid = 12
                            zgrid = 12
                        else if (trim(option(1:7)) == 'GRID') then
                            read(lucmd, *) xsize, xgrid, ysize, ygrid, zsize, zgrid
                        else if (trim(option(1:7)) == 'FIELD') then
                            mep_field = .true.
                        else if (trim(option(1:7)) == 'EXTFLD') then
                            read(lucmd, *) (extfld(i), i = 1, 3)
                            mep_extfld = .true.
                        else if (trim(option(1:7)) == 'LOCFLD') then
                            read(lucmd, *) lf_component
                            mep_lf = .true.
                        else if (trim(option(1:7)) == 'SKIPQM') then
                            mep_qmcube = .false.
                        else if (trim(option(1:7)) == 'SKIPMUL') then
                            mep_mulcube = .false.
                        else if (trim(option(1:7)) == 'INPUT') then
                            mep_input = .true.
                            read(lucmd, '(a240)') h5gridfile
                        else if (option(1:1) == '.' .or. option(1:1) == '*') then
                            backspace(lucmd)
                            exit
                        else if (option(1:1) == '!' .or. option(1:1) == '#') then
                            cycle
                        else
                            call quit('unknown option present in .MEP section.')
                        end if
                    end do
                end if
                pelib_mep = .true.
            ! continuum solvation calculation
            else if (trim(option(2:7)) == 'SOLVAT' .or. trim(option(2:7)) == 'FIXSOL') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) solvent
                else
                    solvent = 'H2O'
                end if
                pelib_sol = .true.
                pelib_diis = .true.
                pelib_fixsol = .true.
                pelib_polar = .true.
                chol = .false.
            ! equilibrium solvation (non-equilibrium is default)
            else if (trim(option(2:7)) == 'EQSOL') then
                pelib_noneq = .false.
            ! set the static dielectric constant
            else if (trim(option(2:7)) == 'EPS') then
                read(lucmd, *) eps
            ! set the optical dielectric constant
            else if (trim(option(2:7)) == 'EPSINF') then
                read(lucmd, *) epsinf
            ! Turn off interaction between induced dipoles and surface charges
            else if (trim(option(2:7)) == 'NOMUQ') then
                pelib_nomuq = .true.
            else if (trim(option(2:7)) == 'NOVMU') then
                pelib_novmu = .true.
            ! specify surface file
            else if (trim(option(2:7)) == 'SURFAC') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, '(a240)') surfile
                end if
                pelib_read_surf = .true.
            ! FixSol solvation using FIXPVA2 tesselation (optional: number of tessera per atom)
            else if (trim(option(2:7)) == 'NTESS') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) ntsatm
                    if ((ntsatm .ne. 60) .and. (ntsatm .ne. 240) .and. (ntsatm .ne. 960)) then
                       call quit('number of tessera per atom must be equal to 60, 240 or 960')
                    end if
                end if
            ! apply external electric field
            else if (trim(option(2:7)) == 'FIELD') then
                pelib_field = .true.
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) (efields(i), i = 1, 3)
                end if
            ! verbose output
            else if (trim(option(2:7)) == 'VERBOS') then
                pelib_verbose = .true.
            ! debug output
            else if (trim(option(2:7)) == 'DEBUG') then
                pelib_debug = .true.
                pelib_verbose = .true.
            ! isotropic polarizabilities
            else if (trim(option(2:7)) == 'ISOPOL') then
                pelib_isopol = .true.
            ! zero out the polarizabilities
            else if (trim(option(2:7)) == 'ZEROPO') then
                pelib_zeropol = .true.
            ! zero out higher-order multipoles
            else if (trim(option(2:7)) == 'ZEROMU') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(lucmd, *) zeromul_order
                    if (zeromul_order < 0) then
                        call quit('ZEROMUL order cannot be negative')
                    end if
                end if
                pelib_zeromul = .true.
            ! skip calculation of multipole-multipole interaction energy
            else if (trim(option(2:7)) == 'SKIPMU') then
                pelib_skipmul = .true.
            else if (trim(option(2:7)) == 'NOSKIP') then
                pelib_skipmul = .false.
            else if (trim(option(2:7)) == 'GAUGE') then
                read(lucmd, '(a80)') option
                backspace(lucmd)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    allocate(gauge_input(3))
                    gauge_input = 0.0_rp
                    do i = 1,3
                        read(lucmd, *) gauge_input(i)
                    enddo
                end if
                pelib_gauge = .true.
            else if (option(1:1) == '*') then
                word = option
                exit
            else if (option(1:1) == '!' .or. option(1:1) == '#') then
                cycle
            else
                write(luout, *) 'unknown option:', option
                call quit('unknown option in *PELIB section')
            end if
        end do
  
        ! check options
        if (pelib_diis .and. .not. pelib_iter) then ! Assume that you want to use a direct solver for FixSol
            pelib_diis = .false.
        end if
  
        if (pelib_sol .and. pelib_iter .and. .not. pelib_fixsol) then
            call quit('.SOLV without .FIXSOL requires .DIRECT')
        end if
  
        if (pelib_mep .and. pelib_cube) then
            call quit('.MEP and .CUBE are not compatible.')
        end if
  
        call qexit('pelib_ifc_input_reader')
  end subroutine pelib_ifc_input_reader

    subroutine pelib_ifc_init()
        use pelib, only: pelib_init

        use mxcent_mod
        implicit none
#include "priunit.h"
#include "nuclei.h"
        call qenter('pelib_ifc_init')
        call pelib_init(lupri, cord(1:3,1:natoms), charge(1:natoms))
        call qexit('pelib_ifc_init')
    end subroutine pelib_ifc_init

    subroutine pelib_ifc_finalize()
        use pelib, only: pelib_finalize
        implicit none
        call qenter('pelib_ifc_finalize')
        call pelib_finalize()
        call qexit('pelib_ifc_finalize')
    end subroutine pelib_ifc_finalize

    subroutine pelib_ifc_fock(dmat, fmat, tot_nrg, el_nrg)
        use pelib, only: pelib_master
        implicit none
#include "priunit.h"
#include "dcbbas.h"
#include "dcbham.h"
        ! Input/output.
        real(rp), intent(in   ) :: dmat(n2bbasx)
        real(rp), intent(inout) :: fmat(n2bbasx)
        real(rp), intent(out  ) :: tot_nrg, el_nrg
        ! Local variables.
        real(rp)              :: nrg(1)
        real(rp), allocatable :: dspmat(:)
        real(rp), allocatable :: fspmat(:)
        integer(ip)           :: i, j, k, l
        ! External.
        real(rp), external :: ddot
        
        call qenter('pelib_ifc_fock')

        allocate(dspmat(nnbbasx))
        allocate(fspmat(nnbbasx))
        dspmat = 0.0_rp
        fspmat = 0.0_rp

        call build_density_matrix(dspmat, dmat, nnbbasx, n2bbasx, ntbas, levyle)

#if defined(VAR_MPI)
        call pelib_ifc_start_slaves(1)
#endif
        ! fspmat is calculated in pelib_master.
        
        call pelib_master(runtype = 'full_fock', &
                          triang  = .true.,      &
                          ndim    = ntbas(0),    &
                          nmats   = 1,           & 
                          denmats = dspmat,      &
                          fckmats = fspmat,      &
                          expvals = nrg)
        tot_nrg = nrg(1)
        el_nrg  = ddot(nnbbasx, dspmat, 1, fspmat, 1)

        call add_to_fock_matrix(fmat, fspmat, nnbbasx, n2bbasx, ntbas, levyle)

        deallocate(dspmat)
        deallocate(fspmat)

        call qexit('pelib_ifc_fock')

    end subroutine pelib_ifc_fock

    subroutine pelib_ifc_energy(dmat)
        use pelib, only: pelib_master
        implicit none
#include "dcbbas.h"
#include "dcbham.h"
        ! Input/output.
        real(rp), dimension(n2bbasx), intent(in) :: dmat(n2bbasx)
        ! Local variables.
        real(rp), dimension(1) :: nrg
        real(rp), allocatable  :: dspmat(:)
        integer(ip)            :: i, j, k, l
      
        call qenter('pelib_ifc_energy')

        allocate(dspmat(nnbbasx))
        dspmat = 0.0_rp

        call build_density_matrix(dspmat, dmat, nnbbasx, n2bbasx, ntbas, levyle)

#if defined(VAR_MPI)
        call pelib_ifc_start_slaves(2)
#endif
        call pelib_master(runtype = 'print_energy', &
                          triang  =  .true.,        &
                          ndim    =  ntbas(0),      &
                          nmats   = 1,              &
                          denmats = dspmat) 

        deallocate(dspmat)

        call qexit('pelib_ifc_energy')  
  
    end subroutine pelib_ifc_energy

    subroutine pelib_ifc_response(fao, dens1, ndmat, dmat_ih_sym)
        use pelib, only: pelib_master
        use pelib_options, only: pelib_gspol
!       Written by : Erik Donovan Hedegrd (edh)
!       Modified by: Ernst Dennis Larsson (edl)
!       Purpose    : Calculate dynamical pelib response part and add it to the
!                    one index transformed fock matrix
        implicit none
#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
        ! Input/output.
        ! Note: dens1 set to dim. n2bbasxq, but is passed from routine GMOLITX1 with N2ORBXQ (for reasons of memory re-use) -
        integer(ip), intent(in   ) :: ndmat              ! number of b-vectors
        integer(ip), intent(in   ) :: dmat_ih_sym(ndmat) ! symmetry wrt. hermicity
        real(rp)   , intent(in   ) :: dens1(n2bbasxq, *) ! one-index-transformed density matrix
        real(rp)   , intent(inout) :: fao(n2bbasxq, *)   ! total one-index transformed fock matrix
        ! Local variables.
        real(rp), allocatable  :: dlmat_triangle(:)
        real(rp), allocatable  :: flmat_triangle(:)
        integer(ip)            :: triangle_index
        integer(ip)            :: matrix_id
        integer(ip), parameter :: skip_id = -1

        call qenter('pelib_ifc_response')

        write(lupri,*)'*** Calculate dynamical PE response'
        allocate(dlmat_triangle(nnbbasx*ndmat))
        allocate(flmat_triangle(nnbbasx*ndmat))
        dlmat_triangle = 0.0_rp
        flmat_triangle = 0.0_rp
      
        ! edh: We only construct dlmat_triangle for the symmetric real part of dens1
        ! ihrmop(l)=1 of dens1 (i.e. transformed density matrix)
        triangle_index = 1
        do matrix_id = 1, ndmat
            if (dmat_ih_sym(matrix_id) /= skip_id) then
                call build_density_matrix(dlmat_triangle,      & 
                                          dens1(:, matrix_id), &
                                          nnbbasx,             &
                                          n2bbasxq,            &
                                          ntbas,               &
                                          levyle,              &
                                          triangle_index)
            end if
        end do
#if defined(VAR_MPI)
        if (.not. pelib_gspol) then
            ! handle early return of pelib_master in case of pelib_gspol
            ! so the mpi calls don't get stuck
            call pelib_ifc_start_slaves(3)
        end if
#endif
        call pelib_master(runtype = 'dynamic_response', &
                          triang  =  .true.,            &
                          ndim    =  ntbas(0),          &
                          nmats   = ndmat,              &
                          denmats = dlmat_triangle,     &
                          fckmats = flmat_triangle)
        ! F = F + M
        triangle_index = 1
        do matrix_id = 1, ndmat
            if (dmat_ih_sym(matrix_id) /= skip_id) then
                call add_to_fock_matrix(fao(:, matrix_id), &
                                        flmat_triangle,    &
                                        nnbbasx,           &
                                        n2bbasxq,          &
                                        ntbas,             &
                                        levyle,            &
                                        triangle_index)
            end if
        end do

        deallocate(dlmat_triangle)
        deallocate(flmat_triangle)

        call qexit('pelib_ifc_response')

    end subroutine pelib_ifc_response

    subroutine build_density_matrix(dspmat, dmat, nnbbasx, n2bbasx, ntbas, levyle, kstart)
        implicit none
        ! Input/output.
        integer(ip), optional, intent(inout) :: kstart
        logical(lp),           intent(in   ) :: levyle
        integer(ip),           intent(in   ) :: ntbas(0:2)
        integer(ip),           intent(in   ) :: n2bbasx
        integer(ip),           intent(in   ) :: nnbbasx
        real(rp)   ,           intent(in   ) :: dmat(n2bbasx)
        real(rp)   ,           intent(inout) :: dspmat(:)
        ! Local variables.
        integer(ip) :: i, j, k
        integer(ip) :: id
        integer(ip) :: total_orbitals
        integer(ip) :: large_orbitals
        integer(ip) :: small_orbitals

        if (present(kstart)) then
            k = kstart
        else
            k = 1
        end if
    
        total_orbitals = ntbas(0)
        large_orbitals = ntbas(1)
        LLblock: do i = 1, large_orbitals
            do j = 1, i
                id = (i - 1) * total_orbitals + j
                if (i == j) then
                    dspmat(k) = 2.0_rp * dmat(id)
                else
                    dspmat(k) = 4.0_rp * dmat(id)
                end if
                k = k + 1
            end do
        end do LLblock

        SSBlock: if (.not. levyle) then
            small_orbitals = large_orbitals + 1
            do i = small_orbitals, total_orbitals
                k = k + large_orbitals
                do j = small_orbitals, i                    
                    id = (i - 1) * total_orbitals + j
                    if (i == j) then
                        dspmat(k) = 2.0_rp * dmat(id)
                    else
                        dspmat(k) = 4.0_rp * dmat(id)
                    end if
                    k = k + 1
                end do
            end do
        end if SSBlock

        if (present(kstart)) kstart = kstart + nnbbasx

    end subroutine build_density_matrix

    subroutine add_to_fock_matrix(fmat, fspmat, nnbbasx, n2bbasx, ntbas, levyle, kstart)
        implicit none
        ! Input/output.
        integer(ip), optional, intent(inout) :: kstart
        logical(lp),           intent(in   ) :: levyle
        integer(ip),           intent(in   ) :: ntbas(0:2)
        integer(ip),           intent(in   ) :: n2bbasx
        integer(ip),           intent(in   ) :: nnbbasx
        real(rp)   ,           intent(in   ) :: fspmat(:)
        real(rp)   ,           intent(inout) :: fmat(n2bbasx)
        ! Local variables.
        integer(ip) :: i, j, k
        integer(ip) :: i_id
        integer(ip) :: j_id
        integer(ip) :: total_orbitals
        integer(ip) :: large_orbitals
        integer(ip) :: small_orbitals

        if (present(kstart)) then
            k = kstart
        else
            k = 1
        end if

        total_orbitals = ntbas(0)
        large_orbitals = ntbas(1)
        LLblock: do i = 1, large_orbitals
            do j = 1, i
                i_id = (j - 1) * total_orbitals + i
                j_id = (i - 1) * total_orbitals + j
                if (i == j) then
                    fmat(j_id) = fmat(j_id) + fspmat(k)
                else
                    fmat(j_id) = fmat(j_id) + fspmat(k)
                    fmat(i_id) = fmat(i_id) + fspmat(k)
                end if
                k = k + 1
            end do
        end do LLblock

        SSBlock: if (.not. levyle) then
            small_orbitals = large_orbitals + 1
            do i = small_orbitals, total_orbitals
                k = k + large_orbitals
                do j = small_orbitals, i
                    i_id = (j - 1) * total_orbitals + i
                    j_id = (i - 1) * total_orbitals + j
                    if (i == j) then
                        fmat(j_id) = fmat(j_id) + fspmat(k)
                    else
                        fmat(j_id) = fmat(j_id) + fspmat(k)
                        fmat(i_id) = fmat(i_id) + fspmat(k)
                    end if
                    k = k + 1
                end do
            end do
        end if SSBlock

        if (present(kstart)) kstart = kstart + nnbbasx

    end subroutine add_to_fock_matrix

   subroutine pelib_ifc_london(fckmats)
! edh modify for DIRAC
      use pelib, only: pelib_master
      implicit none
#include "inforb.h"
      real*8, dimension(3*n2basx), intent(out) :: fckmats
      integer :: i, j, k, l, m
      real*8, dimension(:), allocatable :: fckmats_packed
      call qenter('pelib_ifc_london')
      !if (.not. use_pelib()) call quit('PElib not active')
      allocate(fckmats_packed(3*nnbasx))
#if defined(VAR_MPI)
      call pelib_ifc_start_slaves(4)
#endif
      call pelib_master('magnetic_gradient', &
                        triang  = .true.,    &
                        ndim    = nbast,     &
                        fckmats = fckmats_packed)
      do i = 1, 3
         j = (i - 1) * nnbasx + 1
         k = i * nnbasx
         l = (i - 1) * n2basx + 1
         m = i * n2basx
         call daptge(nbas, fckmats_packed(j:k), fckmats(l:m))
      end do
      deallocate(fckmats_packed)
      call qexit('pelib_ifc_london')
   end subroutine pelib_ifc_london

   subroutine pelib_ifc_localfield(soint, nbast, nelmnt, noptyp)
      use pelib, only: pelib_master
      implicit none
#include "priunit.h"
      integer :: nbast, nelmnt, noptyp
      real*8, dimension(nelmnt,noptyp), intent(inout) :: soint
      real*8, dimension(:), allocatable :: intmats
      integer :: i
      call qenter('pelib_ifc_localfield')
      if (.not. use_pelib()) call quit('PElib not active')
      write(lupri,'(/A/A)') 'Effective external field (EEF) correction'//            &
                            ' is applied to dipole integrals.',                      & 
                            'Note: This correction alters *all* properties (incl.'// &
                            'static properties) that depend on dipole integrals.'         
      if (noptyp /= 3) then
         call quit('ERROR: inconsistent dimension')
      end if
      allocate(intmats(noptyp*nelmnt))
#if defined(VAR_MPI)
      call pelib_ifc_start_slaves(8)
#endif
      call pelib_master(runtype = 'effdipole', &
                        triang  = .true.,      &
                        ndim    = nbast,       &
                        fckmats = intmats)
      do i = 1, 3
         soint(:,i) = soint(:,i) + intmats((i-1)*nelmnt+1:i*nelmnt)
      end do
      call qexit('pelib_ifc_localfield')
   end subroutine pelib_ifc_localfield

#if defined(VAR_MPI)
   subroutine pelib_ifc_slave()
      use interface_to_mpi
      use pelib, only: pelib_slave
      implicit none
#include "infpar.h"
      integer :: runtype
      call qenter('pelib_ifc_slave')
      CALL interface_MPI_BCAST(runtype, 1, MPARID, global_communicator)
      if (runtype == 1) then
          call pelib_slave('full_fock')
      else if (runtype == 2) then
          call pelib_slave('print_energy')
      else if (runtype == 3) then
          call pelib_slave('dynamic_response')
      else if (runtype == 4) then
          call pelib_slave('magnetic_gradient')
      else if (runtype == 5) then
          call pelib_slave('molecular_gradient')
      else if (runtype == 8) then
          call pelib_slave('effdipole')
      end if
      call qexit('pelib_ifc_slave')
   end subroutine pelib_ifc_slave
#endif

#if defined(VAR_MPI)
   subroutine pelib_ifc_start_slaves(runtyp)
      use interface_to_mpi
      implicit none
      integer :: runtyp
#include "infpar.h"
#include "dirac_partask.h"
      call qenter('pelib_ifc_start_slaves')
      call dirac_parctl(PELIB_PAR)
      CALL interface_MPI_BCAST(RUNTYP, 1, MPARID, global_communicator)
      call dirac_parctl(RELEASE_NODES)
      call qexit('pelib_ifc_start_slaves')
   end subroutine pelib_ifc_start_slaves
#endif

!
! edh: The remaining is left-over from dalton: quadratic and cubic response,
! PE-MCSCF (linear and quadratic response + energy).
!      subroutine pe_grad(cref, cmo, cindx, dv, grd, energy, wrk, nwrk)
!!
!!     Written by Erik Donovan Hedegrd (edh) and Jgvan Magnus H. Olsen
!!                based on PCMGRAD
!!
!!     Purpose:  calculate (MCSCF) energy and gradient contribution
!!               from a PE potential using the PE library
!!
!!     Output:
!!     grd       MCSCF gradient with PE contribution added
!!     energy    total PE energy
!!
!! Used from common blocks:
!!   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
!!   INFORB: NNASHX, NNBASX, NNORBX, etc.
!!   INFTAP: LUIT2
!!
!      implicit none
!
!#include "priunit.h"
!#include "infvar.h"
!#include "inforb.h"
!#include "inftap.h"
!
!      integer :: nwrk
!      real*8 :: energy
!      real*8, dimension(*) :: cref, cmo, cindx, dv, grd
!      real*8, dimension(nwrk) :: wrk
!      character*8 :: star8 = '********'
!      character*8 :: solvdi = 'SOLVDIAG'
!      character*8 :: eodata = 'EODATA  '
!
!      logical :: fndlab
!      integer :: nc4, nw4, i
!      real*8 :: solelm, ddot
!      real*8 :: tmo, tac, test
!      real*8, dimension(:), allocatable :: fckmo, fckac
!      real*8, dimension(:), allocatable :: pegrd, diape
!
!      call qenter('pe_grad')
!
!      allocate(fckmo(nnorbx))
!      call pe_fckmo(cmo, fckmo, energy, dv, wrk, nwrk)
!
!      allocate(fckac(nnashx))
!      if (nasht > 0) call getac2(fckmo, fckac)
!
!      tmo = solelm(dv, fckac, fckmo, tac)
!
!      allocate(pegrd(nvarh))
!      pegrd = 0.0d0
!      if (nconf > 1) then
!         ! edh: SOLGC calc. < u | Fg | 0 > + < 0 | Fg | 0 > c_u
!        call solgc(cref, fckac, tac, pegrd, cindx, wrk, nwrk)
!      end if
!      if (nwopt > 0) then
!        ! edh: SOLGO calc. 2 < 0 | [Ers, Fg] | 0 >
!        call solgo(2.0d0, dv, fckmo, pegrd(1+nconf:nvarh))
!      end if
!
!      allocate(diape(nvar))
!      diape = 0.0d0
!      call soldia(tac, fckac, cindx, fckmo, dv, diape, wrk, nwrk)
!      diape = - diape
!      deallocate(fckmo, fckac)
!
!       !--------------- Orthogonality test ----------------
!       test = ddot(nconf, cref, 1, pegrd, 1)
!       if (abs(test) > 1.0d-8) then
!          nwarn = nwarn + 1
!          write(lupri,*) ' >>> PE GRADIENT WARNING <<< '
!          write(lupri,*) ' < CREF | GRAD > =', test
!       end if
!       ! --------------------------------------------------
!       ! Add PE gradient contribution to MCSCF gradient
!       call daxpy(nvarh, 1.0d0, pegrd, 1, grd, 1)
!       deallocate(pegrd)
!       if (luit2 > 0) then
!          nc4 = max(nconf, 4)
!          nw4 = max(nwopt, 4)
!          rewind luit2
!          if (fndlab(eodata,luit2)) backspace luit2
!          write(luit2) star8, star8, star8, solvdi
!          if (nconf > 1) call writt(luit2, nc4, diape)
!          write(luit2) star8, star8, star8, eodata
!       end if
!       call qexit('pe_grad')
!       end subroutine pe_grad
!
!!------------------------------------------------------------------------------
!
!      subroutine pe_fckmo(cmo, fckmo, energy, dv, wrk, nwrk)
!!
!!     Written by Erik Donovan Hedegrd (edh) and Jgvan Magnus H. Olsen
!!
!!     Purpose:  Get PE operator in MO basis using PE library
!!
!!     Output:
!!     fckmo     PE operator in MO basis
!!     energy    PE energy
!!
!!     Used from common blocks:
!!     INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
!!     INFORB: NNASHX, NNBASX, NNORBX, etc.
!!     INFTAP: LUIT2
!!
!
!      use polarizable_embedding, only: pe_master
!
!      implicit none
!
!#include "priunit.h"
!#include "infvar.h"
!#include "inforb.h"
!#include "inftap.h"
!
!      integer :: nwrk
!      real*8, dimension(1) :: energy
!      real*8, dimension(*) :: fckmo
!      real*8, dimension(*) :: cmo, dv
!      real*8, dimension(nwrk) :: wrk
!
!      integer :: i
!      real*8, dimension(:), allocatable :: dcao, dvao, fdtao, fckao
!
!      call qenter('pe_fckmo')
!
!      allocate(dcao(n2basx), dvao(n2basx))
!      call fckden((nisht>0), (nasht>0), dcao, dvao, cmo, dv, wrk, nwrk)
!      if (nisht==0) dcao = 0.0d0
!      if (nasht > 0) dcao = dcao + dvao
!      deallocate(dvao)
!
!      allocate(fdtao(nnbasx))
!      call dgefsp(nbast, dcao, fdtao)
!      deallocate(dcao)
!
!      allocate(fckao(nnbasx))
!      call pe_master(runtype='fock', denmats=fdtao, fckmats=fckao,
!     &               nmats=1, energies=energy)
!      deallocate(fdtao)
!
!      call uthu(fckao, fckmo, cmo, wrk, nbast, norbt)
!      deallocate(fckao)
!
!      call qexit('pe_fckmo')
!
!      end subroutine pe_fckmo
!
!!------------------------------------------------------------------------------
!
!      subroutine pe_lin(ncsim, nosim, bcvecs, bovecs, cref, cmo, cindx,
!     &                  dv, dtv, scvecs, sovecs, orblin, wrk, nwrk)
!!
!! Written by Erik Donovan Hedegrd and Jgvan Magnus H. Olsen
!!            after original code by  Hans Joergen Aa. Jensen
!!
!! Common driver for pe_lnc and pe_lno
!!
!!   Used from common blocks:
!!   INFLIN : NWOPPT,NVARPT
!
!      implicit none
!
!#include "priunit.h"
!#include "inflin.h"
!#include "infvar.h"
!#include "inforb.h"
!
!      logical :: orblin
!      integer :: ncsim, nosim, nwrk
!      real*8, dimension(*) :: bcvecs, bovecs, scvecs, sovecs
!      real*8, dimension(*) :: cmo, cref, cindx, dv, dtv
!      real*8, dimension(nwrk) :: wrk
!
!      integer :: i, nso
!
!      call qenter('pe_lin')
!
!      if (ncsim > 0) then
!        call pe_lnc(ncsim, bcvecs, cref, cmo, cindx, dv, dtv, scvecs,
!     &              wrk, nwrk)
!      end if
!
!      if (nosim > 0) then
!        if (.not. orblin) then
!            nso = nvarpt
!        else
!            nso = nwoppt
!        end if
!        call pe_lno(nosim, bovecs, cref, cmo, cindx, dv, sovecs, nso,
!     &              wrk, nwrk)
!      end if
!
!      call qexit('pe_lin')
!
!      end subroutine pe_lin
!
!!------------------------------------------------------------------------------
!
!      subroutine pe_lnc(ncsim, bcvecs, cref, cmo, cindx, dv, dtv,
!     &                  scvecs, wrk, nwrk)
!!
!!  Written by Erik Donovan Hedegaard and Jgvan Magnus H. Olsen
!!             after original routine by Hans Joergen Aa. Jensen
!!
!!  Purpose:  Calculate Hessian contribution from a polarizable
!!            embedding potantial to a csf trial vector.
!!
!!
!!  Used from common blocks:
!!    INFORB : NNASHX, NNORBX, NNBASX, etc.
!!    INFVAR : NWOPH
!!    INFLIN : NCONST, NVARPT, NWOPPT
!!
!
!      use pe_variables, only: pe_polar
!      use polarizable_embedding, only: pe_master
!
!      implicit none
!
!#include "priunit.h"
!#include "dummy.h"
!#include "inforb.h"
!#include "infvar.h"
!#include "inflin.h"
!#include "infdim.h"
!
!      integer :: ncsim, nwrk
!      real*8, dimension(*) :: bcvecs, cref, cmo, cindx, dv
!      real*8, dimension(nnashx,*) :: dtv
!      real*8, dimension(nvarpt,*) :: scvecs
!      real*8, dimension(nwrk) :: wrk
!
!      logical :: fndlab
!      integer :: i, j, jscvec, mwoph
!      real*8 :: tfxc, tfyc, tfycac, energy, solelm
!      real*8, dimension(:), allocatable :: udtvao, fdtvaos, fxcaos
!      real*8, dimension(:), allocatable :: tfxcacs, fyc, fycac
!      real*8, dimension(:,:), allocatable :: fxcs, fxcacs
!
!      call qenter('pe_lnc')
!
!      allocate(fxcs(nnorbx,ncsim))
!      allocate(fxcacs(nnashx,ncsim))
!      allocate(tfxcacs(ncsim))
!      if (pe_polar) then
!        allocate(udtvao(n2basx))
!        allocate(fdtvaos(ncsim*nnbasx))
!        do i = 1, ncsim
!            j = (i - 1) * nnbasx + 1
!            call fckden(.false., .true., dummy, udtvao, cmo,
!     &                  dtv(:,i), wrk, nwrk)
!            call dgefsp(nbast, udtvao, fdtvaos(j))
!        end do
!        deallocate(udtvao)
!
!        !---------------------------------------
!        ! Fxc = -R<0|Fe|B>Fe in fxcaos
!        !--------------------------------------
!        allocate(fxcaos(ncsim*nnbasx))
!        call pe_master(runtype='response', denmats=fdtvaos,
!     &                 fckmats=fxcaos, nmats=ncsim)
!        deallocate(fdtvaos)
!
!        do i = 1, ncsim
!            j = (i - 1) * nnbasx + 1
!            call uthu(fxcaos(j), fxcs(:,i), cmo, wrk, nbast, norbt)
!            if (nasht > 0) call getac2(fxcs(:,i), fxcacs(:,i))
!            tfxc = solelm(dv, fxcacs(:,i), fxcs(:,i), tfxcacs(i))
!        end do
!        deallocate(fxcaos)
!      end if
!
!      !---------------------------------------
!      ! Fg = Vmul -R<0|Fe|0>Fe in fyc
!      !---------------------------------------
!      allocate(fyc(nnorbx), fycac(nnashx))
!      call pe_fckmo(cmo, fyc, energy, dv, wrk, nwrk)
!      if (nasht > 0) call getac2(fyc, fycac)
!      tfyc = solelm(dv, fycac, fyc, tfycac)
!
!!    ...CSF part of sigma vectors
!      call solsc(ncsim, 0, bcvecs, cref, scvecs, fxcacs, fycac, tfxcacs,
!     &          tfycac, cindx, wrk, nwrk)
!      deallocate(fxcacs, fycac, tfxcacs)
!
!      if (nwoppt > 0) then
!        mwoph = nwoph
!        nwoph = nwoppt
!        jscvec = 1 + nconst
!        do i = 1, ncsim
!            if (pe_polar) then
!                call solgo(2.0d0, dv, fxcs(:,i), scvecs(jscvec,i))
!            end if
!            call solgo(0.0d0, dtv(:,i), fyc, scvecs(jscvec,i))
!        end do
!        nwoph = mwoph
!      end if
!      deallocate(fyc, fxcs)
!
!      call qexit('pe_lnc')
!
!      end subroutine pe_lnc
!
!!------------------------------------------------------------------------------
!
!      subroutine pe_lno(nosim, bovecs, cref, cmo, cindx, dv,
!     &                  sovecs, nso, wrk, nwrk)
!!
!!  Written by Erik Donovan Hedegaard and Jgvan Magnus H. Olsen
!!             after original code by Hans Jorgen Aa. Jensen
!!
!!  Purpose:  Calculate Hessian contribution from a
!!            PE potential to an orbital trial vector.
!!
!!  NSVEC     may be NVAR or NWOPT, dependent on LINTRN
!!
!!
!!  Used from common blocks:
!!    INFORB : NNASHX, NNORBX, NNBASX, etc.
!!    INFVAR : JWOP
!!    INFLIN : NWOPPT, NVARPT, NCONST, NCONRF
!!
!
!      use pe_variables, only: pe_polar
!      use polarizable_embedding, only: pe_master
!
!      implicit none
!
!#include "priunit.h"
!#include "dummy.h"
!#include "inforb.h"
!#include "infvar.h"
!#include "inflin.h"
!
!      integer :: nosim, nso, nwrk
!      real*8, dimension(*) :: cref, cmo, dv, cindx
!      real*8, dimension(nwrk) :: wrk
!      real*8, dimension(nwoppt,*) :: bovecs
!      real*8, dimension(nso,*) :: sovecs
!
!      integer :: i, j, jsovec, mwoph, ncolim
!      logical :: fulhes, fndlab
!      real*8 :: solelm
!      real*8 :: txyo, energy
!      real*8, dimension(:), allocatable :: txyoacs
!      real*8, dimension(:), allocatable :: ubodcao, ubodvao
!      real*8, dimension(:), allocatable :: bodtaos, fxoaos
!      real*8, dimension(:), allocatable :: fckmo, fyo, ufyo
!      real*8, dimension(:,:), allocatable :: ubovecs, fxos
!      real*8, dimension(:,:), allocatable :: fxyos, fxyoacs
!
!
!      call qenter('pe_lno')
!
!      allocate(ubovecs(n2orbx,nosim))
!      if (nosim > 0) then
!          do i = 1, nosim
!            call upkwop(nwoppt, jwop, bovecs(:,i), ubovecs(:,i))
!          end do
!      end if
!
!      !---------------------------------------------------
!      ! 1. Calculation of Fxo = R*<0|Fe(k)|O>Fe
!      !    Store in fxos
!      !---------------------------------------------------
!      if (pe_polar) then
!        allocate(ubodcao(n2basx), ubodvao(n2basx))
!        allocate(bodtaos(nosim*nnbasx))
!        do i = 1, nosim
!            j = (i - 1) * nnbasx + 1
!            call tr1den(cmo, ubovecs(:,i), dv, ubodcao, ubodvao,
!     &                  wrk, nwrk)
!            if (nasht > 0) ubodcao = ubodcao + ubodvao
!            call dgefsp(nbast, ubodcao, bodtaos(j))
!        end do
!        deallocate(ubodcao, ubodvao)
!
!        allocate(fxoaos(nosim*nnbasx))
!        call pe_master(runtype='response', denmats=bodtaos,
!     &                 fckmats=fxoaos, nmats=nosim)
!        deallocate(bodtaos)
!
!
!        allocate(fxos(nnorbx,nosim))
!        do i = 1, nosim
!            j = (i - 1) * nnbasx + 1
!            call uthu(fxoaos(j), fxos(:,i), cmo, wrk, nbast, norbt)
!        end do
!        deallocate(fxoaos)
!      end if
!
!      !---------------------------------------------------
!      ! 2. Calculation of Fyo = V(k) + R<0|F|0>Fe(k)
!      !    Store in fyos
!      !---------------------------------------------------
!      allocate(fckmo(nnorbx))
!      call pe_fckmo(cmo, fckmo, energy, dv, wrk, nwrk)
!
!      allocate(fyo(n2orbx))
!      call dsptsi(norbt, fckmo, fyo)
!      deallocate(fckmo)
!
!      allocate(ufyo(n2orbx), txyoacs(nosim))
!      allocate(fxyos(nnorbx,nosim), fxyoacs(nnashx,nosim))
!      do i = 1, nosim
!        ufyo = 0.0d0
!        call tr1uh1(ubovecs(:,i), fyo, ufyo, 1)
!        call dgetsp(norbt, ufyo, fxyos(:,i))
!        if (pe_polar) then
!            call daxpy(nnorbx, 1.0d0, fxos(:,i), 1, fxyos(:,i), 1)
!        end if
!        if (nasht > 0) then
!            call getac2(fxyos(:,i), fxyoacs(:,i))
!        end if
!        txyo = solelm(dv, fxyoacs(:,i), fxyos(:,i), txyoacs(i))
!      end do
!      !---------------------------------------------------
!      ! 3.   /        <0[Epq,Fxo + Fyo]|0>      \  orbital part
!      !      \ 2<0|Fyo + Fxo|mu> - <0|Fyo|0>*c0 /  CSF part
!      !---------------------------------------------------
!      !     ... CSF part of sigma vectors
!      if (lsymrf == lsymst) then
!        ncolim = 1
!      else
!        ncolim = 0
!      end if
!
!      ! Determine if full Hessian or only orbital Hessian
!      fulhes = (nso == nvarpt)
!      if (fulhes) then
!        jsovec = 1 + nconst
!      else
!        jsovec = 1
!      end if
!
!      if (fulhes .and. (nconst > ncolim)) then
!        call solsc(0, nosim, dummy, cref, sovecs, fxyoacs, dummy,
!     &             txyoacs, dummy, cindx, wrk, nwrk)
!      end if
!
!      ! ... orbital part of sigma vectors
!      mwoph = nwoph
!      nwoph = nwoppt
!      ! ... tell SOLGO only to use the NWOPPT first JWOP entries
!      do i = 1, nosim
!        call solgo(2.0d0, dv, fxyos(:,i), sovecs(jsovec,i))
!      end do
!      nwoph = mwoph
!
!      call qexit('pe_lno')
!
!      end subroutine pe_lno
!
!!------------------------------------------------------------------------------
!
!      subroutine pe_rsplin(ncsim, nosim, bcvecs, bovecs, cref, cmo,
!     &                     cindx, udv, dv, udvtr, dvtr, dtv, dtvtr,
!     &                     scvecs, sovecs, wrk, nwrk)
!
!        implicit none
!
!#include "priunit.h"
!#include "dummy.h"
!#include "wrkrsp.h"
!#include "infrsp.h"
!#include "inftap.h"
!
!        integer :: ncsim, nosim, nwrk
!        real*8, dimension(*) :: bcvecs, bovecs
!        real*8, dimension(*) :: cref, cmo, cindx, udv, dv
!        real*8, dimension(*) :: udvtr, dvtr, dtv, dtvtr
!        real*8, dimension(*) :: scvecs, sovecs
!        real*8, dimension(nwrk) :: wrk
!
!        call qenter('pe_rsplin')
!
!!       if (.not. tdhf) then
!!            write(lupri,*) 'ERROR: PE-MCSCF response not implemented.'
!!           stop 'ERROR: PE-MCSCF response not implemented.'
!!       end if
!
!        if (ncsim > 0 .and. .not. soppa) then
!            call pe_rsplnc(ncsim, bcvecs, cref, cmo, cindx, udv, dv,
!     &                     udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)
!        end if
!
!        if (nosim > 0) then
!            call pe_rsplno(ncsim, nosim, bovecs, cref, cmo, cindx,
!     &                     udv, dv, udvtr, dvtr, sovecs, wrk, nwrk)
!        end if
!
!        call qexit('pe_rsplin')
!
!      end subroutine pe_rsplin
!
!!------------------------------------------------------------------------------
!!
!       subroutine pe_rsplnc(ncsim, bcvecs, cref, cmo, cindx, udv, dv,
!     &                      udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)
!
!         use pe_variables, only: pe_polar, pe_debug
!         use polarizable_embedding, only: pe_master
!
!         implicit none
!!
!#include "priunit.h"
!#include "dummy.h"
!#include "infrsp.h"
!#include "inftap.h"
!#include "wrkrsp.h"
!#include "inforb.h"
!#include "qrinf.h"
!#include "infvar.h"
!!
!! Used from common: n2ashx, kzyvar, kzconf, ncref, maxwop
!      integer :: i, j
!      integer :: ncsim, nwrk
!
!      real*8, dimension(*) :: bcvecs, cref, cmo, cindx, udv, dv
!      real*8, dimension(*) :: udvtr, dvtr
!      real*8, dimension(n2ashx,*) ::  dtv
!      real*8, dimension(n2ashx,*) :: dtvtr
!      real*8, dimension(kzyvar,*) :: scvecs
!      real*8, dimension(nwrk) :: wrk
!
!      real*8 :: ovlap, solelm, tfpeac, tfpe
!      real*8, dimension(:,:), allocatable :: udtv!, udtvao
!      real*8, dimension(:,:), allocatable :: fxcs, fxcacs
!      real*8, dimension(:), allocatable :: udtvao, fuxcs
!      real*8, dimension(:), allocatable :: fdtvaos, fxcaos
!      real*8, dimension(:), allocatable :: tfxc, tfxcacs
!      real*8, dimension(:), allocatable :: fpe, fupe, fpeac
!
!      logical :: lexist, lopen
!      logical :: fndlab
!      logical :: tdm, norho2
!
!      lopen = .false.
!      tdm = .true.
!      norho2 = .true.
!
!
!      call qenter('pe_rsplnc')
!
!      allocate(fxcs(nnorbx,ncsim))
!      allocate(fuxcs(n2orbx))
!      allocate(fxcacs(nnashx,ncsim))
!      allocate(tfxc(ncsim))
!      allocate(tfxcacs(ncsim))
!      fxcs    = 0.0d0
!      fuxcs   = 0.0d0
!      fxcacs  = 0.0d0
!      tfxc    = 0.0d0
!      tfxcacs = 0.0d0
!
!      !-----------------------------------
!      ! Fxc = R*(<0(L)|Fe|0> + <0|Fe|0(R)>)Fe
!      !-----------------------------------
!      if (pe_polar) then
!         call getref(cref, ncref)
!         ! ...Construct <0(L)|...|0> + <0|...|0(R)>
!         allocate(udtv(n2ashx,ncsim))
!         udtv = 0.0d0
!         call rsptdm(ncsim,irefsy,ksymst,ncref,kzconf,cref,
!     &                 bcvecs,udtv,dummy,0,0,.true.,.true.,
!     &                 cindx,wrk,1,nwrk)
!         udtv = -1.0d0*udtv
!!     SUBROUTINE RSPTDM(NCSIM,ILRESY,IRSYM,NCLREF,NCRDIM,CLREF,
!!    *                 CR, RHO1,RHO2, ISPIN1,ISPIN2,TDM,NORHO2,
!!    *                 XINDX,WORK,KFREE,LFREE)
!
!         if ( ncsim > 0 ) then
!            allocate(fdtvaos(nnbasx*ncsim))
!            fdtvaos = 0.0d0
!            allocate(udtvao(n2basx))
!            udtvao = 0.0d0
!            do i = 1, ncsim
!               j = (i - 1) * nnbasx + 1
!                call fckden2(.false.,.true., dummy, udtvao, cmo,
!     &                       udtv(:,i), wrk, nwrk)
!                call dgefsp(nbast, udtvao, fdtvaos(j))
!!                write(lupri,*) 'udtvao matrix'
!!                call output(udtvao,1,nbast,1,nbast,nbast,nbast,1,lupri)
!            end do
!            deallocate(udtv,udtvao)
!         end if
!
!         allocate(fxcaos(ncsim*nnbasx))
!         fxcaos = 0.0d0
!         call pe_master(runtype='response', denmats=fdtvaos,
!     &                 fckmats=fxcaos, nmats=ncsim)
!         deallocate(fdtvaos)
!
!         do i = 1, ncsim
!            j = (i - 1) * nnbasx + 1
!            call uthu(fxcaos(j), fxcs(:,i), cmo, wrk, nbast, norbt)
!            if (nasht > 0) call getac2(fxcs(:,i), fxcacs(:,i))
!!              if (trplet) then
!!                    tfxc = solelm(dvtr, fxcacs(:,i), fxcs(:,i), tfxcac)
!!                    tfxc = tfxcac
!!            else
!            tfxc = solelm(dv, fxcacs(:,i), fxcs(:,i), tfxcacs(i))
!!            end if
!         end do
!
!        deallocate(fxcaos)
!
!      end if
!
!      !---------------------------------------
!      ! Fg = V - <0|F|0>Fe -unpack into fupe
!      !---------------------------------------
!      if (.not. tdhf) then
!          allocate(fpe(nnorbx))
!          fpe = 0.0d0
!          if (lusifc <= 0) then
!              call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',
!     &                    idummy, .false.)
!              lopen = .true.
!          end if
!          rewind(lusifc)
!          call mollab('PETMAT  ', lusifc, lupri)
!          call readt(lusifc, nnorbx, fpe)
!          if (lopen) call gpclose(lusifc, 'KEEP')
!          allocate(fupe(n2orbx), fpeac(nnashx))
!          fupe = 0.0d0
!          fpeac = 0.0d0
!          call dsptsi(norbt, fpe, fupe)
!          if (nasht > 0) call getac2(fpe, fpeac)
!              tfpe = solelm(dv, fpeac, fpe, tfpeac)
!          deallocate(fpe)
!      end if
!
!C    ***************************************
!
!C     Calculate Fxc(Rxc) and Fg(Ryc) contributions to SCVECS(NVAR,NCSIM)
!C     =================================================================
!
!C     ... CSF part of sigma vectors
!C
!      if (pe_debug ) then
!         write(lupri,*)' Linear transformed configuration vector'
!         write(lupri,*)' **** Before slvsc in pe_rsplnc **** '
!         call output(scvecs,1,kzyvar,1,ncsim,kzyvar,ncsim,1,lupri)
!      endif
!
!      call slvsc(ncsim, 0, nnashx, bcvecs, cref, scvecs, fxcacs,
!     &           fpeac, tfxcacs, tfpeac, cindx, wrk, nwrk)
!      deallocate(fxcacs, tfxcacs, fpeac)
!
!      if (pe_debug) then
!          write(lupri,*)' Linear transformed configuration vector'
!          write(lupri,*)' **** After slvsc in pe_rsplnc **** '
!          call output(scvecs,1,kzyvar,1,ncsim,kzyvar,ncsim,1,lupri)
!      end if
!C
!C     ... orbital part of sigma vector(s)
!C
!      if (kzwopt .gt. 0) then
!          do i = 1,ncsim
!             fuxcs = 0.0d0
!             call dsptsi(norbt,fxcs(:,i), fuxcs)
!!            if (trplet) then
!!               call slvsor(.true.,.false.,1,udvtr,
!!    &                      scvecs(1,i),fuxcs)
!!            else
!             call slvsor(.true.,.true., 1, udv, scvecs(1,i), fuxcs)
!!            end if
!             if (pe_debug) then
!                 write(lupri,*)' **** After slvsor in pe_rsplnc **** '
!                 write(lupri,*)
!     &           ' Orbital part of linear transformed conf vec no',i
!                 write(lupri,*)' Txc contribution'
!                 call output(scvecs(1,i),1,kzyvar,1,1,kzyvar,1,1,lupri)
!             end if
!
!!            if (trplet) then
!!               call slvsor(.false.,.false., 1, dtvtr(1,i),
!!    *                       scvecs(1,i),fupe)
!!            else
!             call slvsor(.false.,.false.,1, dtv(1,i), scvecs(1,i),fupe)
!!            end if
!             if (pe_debug ) then
!                 write(lupri,*)
!     &           ' Orbital part of linear transformed conf vec no',i
!                 write(lupri,*)' Tg contribution'
!                 call output(scvecs(1,i),1,kzyvar,1,1,kzyvar,1,1,lupri)
!             end if
!         end do
!         deallocate(fupe, fuxcs)
!
!         if (pe_debug ) then
!               write(lupri,*)' linear transformed configuration vector'
!               write(lupri,*)' **** after slvsor  in pe_rsplnc **** '
!               call output(scvecs,1,kzyvar,1,ncsim,kzyvar,ncsim,1,lupri)
!         end if
!      end if
!
!      if (ncref .ne. kzconf) call quit('pe_rsplnc: ncref .ne. kzconf')
!
!      call qexit('pe_rsplnc')
!
!      end subroutine pe_rsplnc
!!
!!------------------------------------------------------------------------------
!
!      subroutine pe_rsplno(ncsim, nosim, bovecs, cref, cmo, cindx,
!     &                     udv, dv, udvtr, dvtr, sovecs, wrk, nwrk)
!
!        use pe_variables, only: pe_polar, pe_gspol
!        use polarizable_embedding, only: pe_master
!
!        implicit none
!
!#include "priunit.h"
!#include "dummy.h"
!#include "wrkrsp.h"
!#include "inforb.h"
!#include "infrsp.h"
!#include "inftap.h"
!
!        integer :: nosim, ncsim, nwrk
!        real*8, dimension(*) :: bovecs
!        real*8, dimension(kzyvar,*) :: sovecs
!        real*8, dimension(*) :: cref, cmo, cindx, udv, dv, udvtr, dvtr
!        real*8, dimension(nwrk) :: wrk
!
!        integer :: i, j
!        real*8 :: txyo
!        real*8 :: ddot, slvqlm
!        real*8, dimension(:), allocatable :: dcao, dvao
!!        real*8, dimension(:), allocatable :: dcaotr, dvaotr
!        real*8, dimension(:), allocatable :: daos, fckaos
!        real*8, dimension(:), allocatable :: daotrs
!!        real*8, dimension(:), allocatable :: daotrs, fckaotrs
!        real*8, dimension(:), allocatable :: evec
!        real*8, dimension(:,:), allocatable :: ubovecs, evecs, eacs
!!        real*8, dimension(:,:), allocatable :: evectrs, eactrs
!        real*8, dimension(:), allocatable :: fpemo,fupemo
!        real*8, dimension(:), allocatable :: txyoacs
!        real*8, dimension(:), allocatable :: ovlp
!        logical :: lexist, lopen
!
!        ! return if no polarization and not MCSCF
!        if (tdhf .and. .not. pe_polar) then
!            return
!        ! no polarization for triplet excitations in closed shell SCF
!        else if ((nasht == 0) .and. trplet) then
!            return
!        ! ground state polarization approximation
!        else if (pe_gspol) then
!            return
!        ! triplet response for open shell systems not ready yet
!        else if ((nasht > 0) .and. trplet) then
!            call quit('ERROR: triplet operators for open shell'//
!     &                ' systems not implemented')
!        end if
!
!        lopen = .false.
!
!        call qenter('pe_rsplno')
!
!        if (.not. tdhf) then
!            ! Read Fg = V - <0|F|0>Fe from file
!            allocate(fpemo(nnorbx))
!            if (lusifc <= 0) then
!                call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',
!     &                      idummy, .false.)
!                lopen = .true.
!            end if
!            rewind(lusifc)
!            call mollab('PETMAT  ', lusifc, lupri)
!            call readt(lusifc, nnorbx, fpemo)
!            if (lopen) call gpclose(lusifc, 'KEEP')
!            allocate(fupemo(n2orbx))
!            call dsptsi(norbt, fpemo, fupemo)
!            deallocate(fpemo)
!        end if
!
!        allocate(ubovecs(n2orbx,nosim))
!        call rspzym(nosim, bovecs, ubovecs)
!
!        ubovecs = - ubovecs
!
!        allocate(dcao(n2basx), dvao(n2basx), daos(nosim*nnbasx))
!!        if (trplet) then
!!            allocate(dcaotr(n2basx), dvaotr(n2basx),
!!     &               daotrs(nosim*nnbasx))
!!        end if
!        !---------------------------------------
!        ! Calculate Fxo = <0|Fe(k)|0>Fe
!        !---------------------------------------
!        do i = 1, nosim
!            j = (i - 1) * nnbasx + 1
!            call deq27(cmo, ubovecs(:,i), udv, dcao, dvao,
!     &                 wrk, nwrk)
!!            edh: Does deq27 assume symmetric matrix? is this
!!            reasonable?
!!            if (trplet) then
!!                call deq27(cmo, ubovecs(:,i), udvtr, dcaotr, dvaotr,
!!    &                     wrk, nwrk)
!!            end if
!            if (nasht > 0) then
!                dcao = dcao + dvao
!!                if (trplet) then
!!                    dcaotr = dcaotr + dvaotr
!!                end if
!            end if
!            call dgefsp(nbast, dcao, daos(j))
!!            if (trplet) then
!!                call dgefsp(nbast, dcaotr, daotrs(j))
!!            end if
!        end do
!        deallocate(dcao, dvao)
!!        if (trplet) then
!!            deallocate(dcaotr, dvaotr)
!!        end if
!
!        allocate(fckaos(nosim*nnbasx))
!        call pe_master(runtype='response', denmats=daos,
!     &                 fckmats=fckaos, nmats=nosim)
!        deallocate(daos)
!!        if (trplet) then
!!            allocate(fckaotrs(nosim*nnbasx))
!!            call pe_master(runtype='response', denmats=daotrs,
!!     &                     fckmats=fckaotrs, nmats=nosim)
!!            deallocate(daotrs)
!!        end if
!
!        allocate(evec(nnorbx))
!        allocate(evecs(n2orbx,nosim))
!        evecs = 0.0d0
!        if (.not. tdhf) then
!            allocate(eacs(n2ashx,nosim))
!            allocate(txyoacs(nosim))
!            eacs = 0.0d0
!            txyoacs = 0.0d0
!        end if
!!        if (trplet) then
!!            allocate(evectrs(n2orbx,nosim), eactrs(n2ashx,nosim))
!!        end if
!
!        do i = 1, nosim
!            j = (i - 1) * nnbasx + 1
!            call uthu(fckaos(j), evec, cmo, wrk, nbast, norbt)
!            call dsptsi(norbt, evec, evecs(:,i))
!!            if (trplet) then
!!                uthu(fckaotrs(j), evectrs, cmo, wrk, nbast, norbt)
!!                call dsptsi(norbt, evectrs, evecstrs(:,i))
!!            end if
!
!            !---------------------------------------
!            ! Fyo = V(k) - <0|F|0>Fe(k)
!            !---------------------------------------
!            if (.not. tdhf) then
!                call onexh1(ubovecs(:,i), fupemo, evecs(:,i))
!                call getacq(evecs(:,i), eacs(:,i))
!                txyo = slvqlm(udv, eacs(:,i), evecs(:,i), txyoacs(i))
!!                if (trplet) then
!!                    call getacq(evecstrs, eacstrs)
!!                    txyot  = slvqlm(udvtr, eacstrs, evecstrs, fyoat(i))
!!                end if
!            end if
!!            if (trplet) then
!!                call uthu(fckaotrs(j), evec, cmo, wrk, nbast, norbt)
!!                call dsptsi(norbt, evec, evectrs(:,i))
!!            end if
!!            if (nasht > 0) then
!!                call getacq(evecs(:,i), eacs(:,i))
!!                if (trplet) then
!!                    call getacq(evectrs(:,i), eactrs(:,i))
!!                end if
!!            end if
!!            tr = solelm(dv, fxyoacs(:,i), fxyos(:,i), txyoacs(i))
!        end do
!
!        deallocate(evec)
!        if (.not. tdhf) then
!            deallocate(fupemo)
!            call slvsc(0, nosim, n2ashx, dummy, cref, sovecs, eacs,
!     &                 dummy, txyoacs, dummy, cindx, wrk, nwrk)
!            deallocate(eacs)
!            deallocate(txyoacs)
!        end if
!
!!        if (trplet) then
!!            call slvsor(.true., .false., nosim, udvtr, sovecs, evectrs)
!!            call slvsor(.true., .true., nosim, udv, sovecs, evecs)
!!        else
!            call slvsor(.true., .true., nosim, udv, sovecs, evecs)
!!        end if
!
!        deallocate(evecs)
!!        if (trplet) then
!!            deallocate(evectrs, eactrs)
!!        end if
!
!!        write(lupri,*)' Linear transformed orbital vector'
!!        call output(sovecs,1,kzyvar,1,nosim,kzyvar,nosim,1,lupri)
!
!        call qexit('pe_rsplno')
!
!      end subroutine pe_rsplno
!
!!------------------------------------------------------------------------------
!
!
!      subroutine pe_rspmcqr(vecb, vecc, etrs, xindx, zymb, zymc,
!     &                      den1, udv, wrk, lfree, kzyva, kzyvb, kzyvc,
!     &                      isyma, isymb, isymc, cmo, mjwop)
!
!         use pe_variables, only: pe_polar
!         use polarizable_embedding, only: pe_master
!
!         implicit none
!
!#include "inforb.h"
!#include "infvar.h"
!#include "infdim.h"
!#include "qrinf.h"
!#include "priunit.h"
!#include "dummy.h"
!#include "inftap.h"
!#include "infrsp.h"
!#include "wrkrsp.h"
!
!         integer :: kzyva, kzyvb, kzyvc
!         integer :: isyma, isymb, isymc, isymbc
!         integer :: lfree
!         integer :: ilsym, irsym, ncl, ncr, kzvarl, kzvarr
!         integer :: isymdn, isymst
!         integer :: kcref, lorb, lcon, nzyvec, nzcvec
!         integer :: iprone, nzconf, nzvar
!         integer :: n2ash
!
!         real*8 :: ovlap
!         real*8 :: fact
!         real*8 :: ddot
!
!         real*8, dimension(*) :: wrk
!         real*8, dimension(1) :: tmpwrk
!         real*8, dimension(*) :: cmo, xindx
!
!         real*8, dimension(kzyva) :: etrs
!         real*8, dimension(kzyvb) :: vecb
!         real*8, dimension(kzyvc) :: vecc
!
!         real*8, dimension(norbt,norbt) :: zymb, zymc
!         real*8, dimension(nashdi,nashdi) :: udv, den1
!         real*8, dimension(nnashx) :: dv
!
!         integer, dimension(2,maxwop,8) :: mjwop
!
!         real*8, dimension(:), allocatable :: fpe
!         real*8, dimension(:), allocatable :: cref
!         real*8, dimension(:), allocatable :: dcaos, fcaos
!         real*8, dimension(:), allocatable :: udtv, udtvao
!         real*8, dimension(:), allocatable :: dvaao, dvbao, dvatr
!         real*8, dimension(:), allocatable :: udcao, udvao
!         real*8, dimension(:), allocatable :: udcmo, udvmo
!         real*8, dimension(:), allocatable :: fcmo
!
!         real*8, dimension(:,:), allocatable :: dva, dvb
!         real*8, dimension(:,:), allocatable :: fupe
!         real*8, dimension(:,:), allocatable :: fxpeb,fxpec, fx2pe
!         real*8, dimension(:,:), allocatable :: fxo1k, fxc1s
!         real*8, dimension(:,:), allocatable :: fxo2k, fxc2s
!         real*8, dimension(:,:), allocatable :: fcas2_1, fcas2_2
!         real*8, dimension(:,:), allocatable :: fcas3_1, fcas3_2
!         real*8, dimension(:,:), allocatable :: fxo, fxo1k2k, fxo2k1k
!         real*8, dimension(:,:), allocatable :: fxc1s2s
!         real*8, dimension(:,:), allocatable :: fxc1s2k, fxc1k2s
!
!         logical :: lexist, lopen
!         logical :: fndlab
!
!         lopen = .false.
!
!         call qenter('pe_rspmcqr')
!
!         call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
!         call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)
!
!!         CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
!!         CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)
!
!         !-----------------------------------------------------------
!         ! Get Fg = Vmul - R*<0|F|>Fe from file
!         !-----------------------------------------------------------
!         if (.not. tdhf) then
!            allocate(fpe(nnorbx))
!            if (lusifc <= 0) then
!                call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',
!     &       idummy, .false.)
!                    lopen = .true.
!            end if
!            rewind(lusifc)
!            call mollab('PETMAT  ', lusifc, lupri)
!            call readt(lusifc, nnorbx, fpe)
!            if (lopen) call gpclose(lusifc, 'KEEP')
!            allocate(fupe(norbt,norbt))
!            call dsptsi(norbt, fpe, fupe)
!            deallocate(fpe)
!         end if
!         !-----------------------------------------------------------
!         ! Density Factory ...
!         !-----------------------------------------------------------
!
!         allocate(cref(mzconf(1)))
!         call getref(cref, mzconf(1))
!
!         if (pe_polar) then
!            allocate(udcao(n2basx))
!            allocate(udvao(n2basx))
!            if (.not. tdhf) then
!               allocate(dcaos(10*nnbasx))
!               else
!               allocate(dcaos(4*nnbasx))
!            end if
!            dcaos = 0.0d0
!
!            !  DTX = D_pq(k1) = <0|[k1,Epq]|0>
!            allocate(udcmo(n2orbx),udvmo(n2orbx))
!            udcmo = 0.0d0
!            udvmo = 0.0d0
!            call deq27mo(isymb, zymb, udv, udcmo, udvmo,
!     &                   wrk, lfree)
!            if (nasht > 0) then
!               udcmo = udcmo + udvmo
!            end if
!            udcao = 0.0d0
!            call motoao(udcmo,udcao,cmo,isymb,wrk,lfree)
!            call dgefsp(nbast, udcao, dcaos(1:nnbasx))
!            ! needed to fit with HF code
!            dcaos(1:nnbasx) = 0.5d0*dcaos(1:nnbasx)
!
!            !  DT2X = D_pq(k2,k1) = <0|[k2,[k1,Epq]|0>
!            udvmo = 0.0d0
!            call oitd1(isymc,zymc,udcmo,udvmo,isymb)
!            ! DT2X in udvmo (re-used to save memory)
!            udcao = 0.0d0
!            isymbc = muld2h(isymb,isymc)
!            call motoao(udvmo,udcao,cmo,isymbc,wrk,lfree)
!            call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))
!            ! needed to fit with HF code
!            dcaos(nnbasx+1:2*nnbasx) = 0.5d0*dcaos(nnbasx+1:2*nnbasx)
!
!            !  DTX = D_pq(k2) = <0|[k2,Epq]|0>
!            udcmo = 0.0d0
!            udvmo = 0.0d0
!            call deq27mo(isymc, zymc, udv, udcmo, udvmo,
!     &                   wrk, lfree)
!            if (nasht > 0) then
!               udcmo = udcmo + udvmo
!            end if
!            udcao = 0.0d0
!            call motoao(udcmo,udcao,cmo,isymb,wrk,lfree)
!            call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))
!            ! needed to fit with HF code
!            dcaos(2*nnbasx+1:3*nnbasx) =
!     &      0.5d0*dcaos(2*nnbasx+1:3*nnbasx)
!
!            !  DT2X = D_pq(k1,k2) = <0|[k1,[k2,Epq]|0>
!            udvmo = 0.0d0
!            call oitd1(isymb,zymb,udcmo,udvmo,isymc)
!            ! DT2X in udvmo (re-used to save memory)
!            udcao = 0.0d0
!            isymbc = muld2h(isymc,isymb)
!            call motoao(udvmo,udcao,cmo,isymbc,wrk,lfree)
!            call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))
!            ! needed to fit with HF code
!            dcaos(3*nnbasx+1:4*nnbasx) =
!     &      0.5d0*dcaos(3*nnbasx+1:4*nnbasx)
!            deallocate(udcmo,udvmo)
!
!            if (tdhf) then
!            write(lupri,*) 'PE-DFT or HF QR detected: Skipping CI dens.'
!            end if
!            if (.not. tdhf ) then
!            write(lupri,*) 'PE-MCSCF QR detected: Constructing CI dens.'
!
!               ! Construct the density matrix <02L|..|0> + <0|..|02R>
!               ilsym  = irefsy
!               irsym  = muld2h(irefsy,isymc)
!               ncl    = mzconf(1)
!               ncr    = mzconf(isymc)
!               kzvarl = mzconf(1)
!               kzvarr = mzyvar(isymc)
!
!               den1 = 0.0d0 ! edh: This is equal to udtv later...
!               allocate(udtv(n2ashx), udtvao(n2basx))
!
!               udtv = 0.0d0
!               udtvao = 0.0d0
!               call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
!     &                     cref, vecc, ovlap, udtv, dummy, 0 ,0, .true.,
!     &                    .true., xindx, wrk, 1, lfree, .true.)
!               call fckden2(.false.,.true., dummy, udtvao, cmo,
!     &                   udtv, wrk, lfree)
!               call dgefsp(nbast, udtvao, dcaos(4*nnbasx+1:5*nnbasx))
!               dcaos(4*nnbasx+1:5*nnbasx) =
!     &         1.0d0*dcaos(4*nnbasx+1:5*nnbasx)
!
!               ! Construct the density matrix <01L|..|0> + <0|..|01R>
!               ilsym  = irefsy
!               irsym  = muld2h(irefsy,isymb)
!               ncl    = mzconf(1)
!               ncr    = mzconf(isymb)
!               kzvarl = mzconf(1)
!               kzvarr = mzyvar(isymb)
!
!               udtv = 0.0d0
!               udtvao = 0.0d0
!               call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
!     &                     cref, vecb, ovlap, udtv, dummy, 0 ,0, .true.,
!     &                     .true., xindx, wrk, 1, lfree, .true.)
!               call fckden2(.false.,.true., dummy, udtvao, cmo,
!     &                      udtv, wrk, lfree)
!               call dgefsp(nbast, udtvao, dcaos(5*nnbasx+1:6*nnbasx))
!                dcaos(5*nnbasx+1:6*nnbasx) =
!     &         1.0d0*dcaos(5*nnbasx+1:6*nnbasx)
!
!               if (mzconf(isymb) .gt. 0 .and. mzconf(isymc) .gt. 0) then
!
!                  ! Construct <01L|..|02R> + <02L|..|01R> density
!                  ilsym  = muld2h(irefsy,isymb)
!                  irsym  = muld2h(irefsy,isymc)
!                  ncl    = mzconf(isymb)
!                  ncr    = mzconf(isymc)
!                  kzvarl = mzyvar(isymb)
!                  kzvarr = mzyvar(isymc)
!                  isymdn = muld2h(ilsym,irsym)
!
!                  udtv = 0.0d0
!                  udtvao = 0.0d0
!                  call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
!     &                        vecb, vecc, ovlap, udtv, dummy, 0 ,0,
!     &                       .true., .true., xindx, wrk, 1, lfree,
!     &                       .false.)
!                  call fckden2(.false.,.true., dummy, udtvao, cmo,
!     &                       udtv, wrk, lfree)
!                  call dgefsp(nbast, udtvao, dcaos(6*nnbasx+1:7*nnbasx))
!               end if
!               dcaos(6*nnbasx+1:7*nnbasx) =
!     &         1.0d0*dcaos(6*nnbasx+1:7*nnbasx)
!
!               ! D_pq = <0|Epq|0>
!               udcao = 0.0d0
!               udvao = 0.0d0
!               call dgefsp(nasht, udv, dv)
!               call fckden((nisht>0), (nasht>0), udcao, udvao,
!     &                      cmo, dv, wrk, lfree)
!               if (nisht==0) udcao = 0.0d0
!               udcao = udcao + udvao
!               call dgefsp(nbast, udcao, dcaos(7*nnbasx+1:8*nnbasx))
!               dcaos(7*nnbasx+1:8*nnbasx) =
!     &         0.5d0*dcaos(7*nnbasx+1:8*nnbasx)
!
!               ! D_pq(S1,k2) = <01L|[k2,Epq]|0> + <0|[k2,Epq]|01R>
!               allocate(dva(norbt,nasht), dvb(norbt,nasht))
!               allocate(dvaao(n2basx), dvbao(n2basx), dvatr(n2basx))
!               dva    = 0.0d0
!               dvb    = 0.0d0
!               dvaao  = 0.0d0
!               dvbao  = 0.0d0
!               dvatr  = 0.0d0
!               udtvao = 0.0d0
!               call rsptr1(1, udv, zymb, dva, dvb)
!               call fckden2(.false.,.true., dummy, dvaao, cmo,
!     &                      dva, wrk, lfree)
!               call fckden2(.false.,.true., dummy, dvbao, cmo,
!     &                      dvb, wrk, lfree)
!               call mtrsp(nbast, nbast, dvaao, nbast, dvatr, nbast)
!               udtvao = dvbao - dvatr
!               call dgefsp(nbast, udtvao, dcaos(8*nnbasx+1:9*nnbasx))
!               dcaos(8*nnbasx+1:9*nnbasx) =
!     &         2.0d0*dcaos(8*nnbasx+1:9*nnbasx)
!
!               ! D_pq(k1,S2) = <02L|[k1,Epq]|0> + <0|[k1,Epq]|02R>
!               dva    = 0.0d0
!               dvb    = 0.0d0
!               dvaao  = 0.0d0
!               dvbao  = 0.0d0
!               dvatr  = 0.0d0
!               udtvao = 0.0d0
!               call rsptr1(1, udv, zymc, dva, dvb)
!               call fckden2(.false.,.true., dummy, dvaao, cmo,
!     &                      dva, wrk, lfree)
!               call fckden2(.false.,.true., dummy, dvbao, cmo,
!     &                      dvb, wrk, lfree)
!               call mtrsp(nbast, nbast, dvaao, nbast, dvatr, nbast)
!               udtvao = dvbao - dvatr
!               call dgefsp(nbast, udtvao, dcaos(9*nnbasx+1:10*nnbasx))
!               dcaos(9*nnbasx+1:10*nnbasx) =
!     &         2.0d0*dcaos(9*nnbasx+1:10*nnbasx)
!
!               deallocate(dva, dvb, dvaao, dvbao, dvatr)
!               deallocate(udtv,udtvao)
!
!            end if
!            deallocate(udcao, udvao)
!            !-----------------------------------------------------------
!            ! Calculate PE response operators in AO basis
!            !-----------------------------------------------------------
!
!            if (.not. tdhf) then
!               allocate(fcaos(10*nnbasx))
!               fcaos = 0.0d0
!               call pe_master(runtype='response',
!     &                        denmats=dcaos, fckmats=fcaos,
!     &                        nmats=10)
!            else
!               allocate(fcaos(4*nnbasx))
!               fcaos = 0.0d0
!               call pe_master(runtype='response',
!     &                        denmats=dcaos, fckmats=fcaos,
!     &                        nmats=4)
!            end if
!            deallocate(dcaos)
!
!         end if
!
!         if ( .not. tdhf ) then
!            !-----------------------------------------------------------
!            !case 1
!            !-----------------------------------------------------------
!            if ( mzconf(isymb) .eq. 0 .or. mzconf(isymc) .eq. 0 ) return
!
!            !/   <01L| [qj,TB] |02R>  + <02L| [qj,TB] |01R>  \
!            !|                       0                       |
!            !|   <01L| [qj+,TB] |02R> + <02L| [qj+,TB] |01R> |
!            !\                       0                       /
!
!            ! ionstruct <01L|..|02R> + <02L|..|01R> density
!            ilsym  = muld2h(irefsy,isymb)
!            irsym  = muld2h(irefsy,isymc)
!            ncl    = mzconf(isymb)
!            ncr    = mzconf(isymc)
!            kzvarl = mzyvar(isymb)
!            kzvarr = mzyvar(isymc)
!
!            den1 = 0.0d0
!            call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
!     &                  vecb, vecc, ovlap, den1, dummy, 0, 0, .true.,
!     &                  .true., xindx, wrk, 1, lfree, .false.)
!
!            ! Make the gradient
!            isymdn = muld2h(ilsym,irsym)
!
!            if ( mzwopt(isyma) .gt. 0 ) then
!               call orbsx(1, isyma, kzyva, etrs, fupe, ovlap,
!     &                    isymdn, den1, mjwop, 1, lfree)
!           end if
!        end if
!        !-----------------------------------------------------------
!        !case 2
!        !-----------------------------------------------------------
!
!        if (pe_polar) then
!
!           allocate(fcmo(nnorbx))
!           allocate(fxo1k(norbt,norbt))
!           fxo1k = 0.0d0
!           if (.not. tdhf) then
!              allocate(fxc1s(norbt,norbt))
!              fxc1s = 0.0d0
!           end if
!
!           ! Fxo = R*<0|Fe(1k)|0>Fe
!           fcmo = 0.0d0
!           call uthu(2.0d0*fcaos(1:nnbasx), fcmo, cmo,
!     &               wrk, nbast, norbt)
!           call dsptsi(norbt, fcmo, fxo1k)
!
!           if (.not. tdhf) then
!              ! Fxc(1S) = ( R*<01lE|0>+<0|E01R> )Fe
!              ! edh: Should it be 1.0d0 or 2.0d0 ???
!              fcmo = 0.0d0
!              call uthu(1.0d0*fcaos(4*nnbasx+1:5*nnbasx), fcmo, cmo,
!     &                  wrk, nbast, norbt)
!              call dsptsi(norbt, fcmo, fxc1s)
!           end if
!
!           ! fcas2_1 = Fa[1](1k)
!           allocate(fcas2_1(norbt,norbt))
!           fcas2_1 = 0.0d0
!           if (.not. tdhf) then
!              fcas2_1 = fxo1k + fxc1s
!              deallocate(fxc1s)
!           else
!              fcas2_1 = fxo1k
!           end if
!           deallocate(fxo1k)
!
!           if (.not. tdhf) then
!              if (mzconf(isymc) .le. 0) return
!
!             !/   0    \
!             !| Sj(2)  | * <0| Fa[1](1k) |0>
!             !|   0    |
!             !\ Sj(2)* /
!
!              if (isyma .eq. isymc) then
!                  ovlap = 1.0d0
!                  call melone(fcas2_1, 1, udv, ovlap, fact,
!     &                        200,'fact for Fxo(1k) + Fxc(1S) ')
!                  nzconf = mzconf(isyma)
!                  nzvar  = mzvar(isyma)
!                  call daxpy(nzconf, fact, vecc, 1, etrs, 1)
!                  call daxpy(nzconf,fact,
!     &                         vecc(nzvar+1), 1, etrs(nzvar+1), 1)
!              end if
!           end if
!
!           allocate(fxo2k(norbt,norbt))
!           fxo2k = 0.0d0
!           if (.not. tdhf) then
!              allocate(fxc2s(norbt,norbt))
!              fxc2s = 0.0d0
!           end if
!
!           ! Fxo(2k) = R*<0|[2k,Epq]|0>Fe
!           fcmo = 0.0d0
!           call uthu(2.0d0*fcaos(2*nnbasx+1:3*nnbasx),
!     &                fcmo, cmo, wrk, nbast, norbt)
!           call dsptsi(norbt, fcmo, fxo2k)
!
!           if (.not. tdhf) then
!              ! Fxc(1S) = ( R*<01lE|0>+<0|E01R> )Fe
!              ! edh: Should it be 1.0d0 or 2.0d0 ???
!              fcmo = 0.0d0
!              call uthu(1.0d0*fcaos(5*nnbasx+1:6*nnbasx), fcmo, cmo,
!     &                   wrk, nbast, norbt)
!               call dsptsi(norbt, fcmo, fxc2s)
!           end if
!
!           ! fcas2_2 = Fa[1](2k)
!           allocate(fcas2_2(norbt,norbt))
!           fcas2_2 = 0.0d0
!           if (.not. tdhf) then
!              fcas2_2 = fxo2k + fxc2s
!              deallocate(fxc2s)
!           else
!              fcas2_2 = fxo2k
!           end if
!           deallocate(fxo2k)
!
!           if (.not. tdhf) then
!              if (mzconf(isymb) .le. 0) return
!
!             !/   0    \
!             !| Sj(1)  | * <0| Fa[1](2k) |0>
!             !|   0    |
!             !\ Sj(1)* /
!
!              if (isyma .eq. isymb) then
!                  ovlap = 1.0d0
!                  call melone(fcas2_2, 1, udv, ovlap, fact,
!     &                        200,'fact for Fxo(1k) + Fxc(1S) ')
!                    nzconf = mzconf(isyma)
!                    nzvar  = mzvar(isyma)
!                    call daxpy(nzconf, fact, vecb, 1, etrs, 1)
!                    call daxpy(nzconf,fact,
!     &                         vecb(nzvar+1), 1, etrs(nzvar+1), 1)
!              end if
!           end if
!        end if
!!       !-----------------------------------------------------------
!!       !case 3
!!       !-----------------------------------------------------------
!        allocate(fxpeb(norbt,norbt))
!        fxpeb = 0.0d0
!        call oith1(isymb,zymb,fupe,fxpeb,1)
!
!        ! fcas3_1 = Fg(k1) + F[a1]
!        allocate(fcas3_1(norbt,norbt))
!        if (pe_polar) then
!            fcas3_1 = fxpeb + fcas2_1
!        else
!            fcas3_1 = fxpeb
!        end if
!
!        if (.not. tdhf) then
!
!           if (mzconf(isymc) .le. 0) return
!
!          !/   <0| [qj,TD1] |02R>  + <02L| [qj,TD1] |0>  \
!          !|   <j| TD1 |02R>                             |
!          !|   <0| [qj+,TD1] |02R> + <02L| [qj+,TD1] |0> |
!          !\  -<02L| TD1 |j>                             /
!
!          ! 1a. Construct the density matrix <02L|..|0> + <0|..|02R>
!           ilsym  = irefsy
!           irsym  = muld2h(irefsy,isymc)
!           ncl    = mzconf(1)
!           ncr    = mzconf(isymc)
!           kzvarl = mzconf(1)
!           kzvarr = mzyvar(isymc)
!
!           den1 = 0.0d0
!           call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
!     &                 cref, vecc, ovlap, den1, dummy, 0 ,0, .true.,
!     &                 .true., xindx, wrk, 1, lfree, .true.)
!
!!          1b. Make the gradient
!           isymdn = muld2h(ilsym,irsym)
!           isymst = muld2h(isyma,irefsy)
!           if ( isymst .eq. irefsy ) then
!              lcon = ( mzconf(isyma) .gt. 1 )
!           else
!              lcon = ( mzconf(isyma) .gt. 0 )
!           end if
!           lorb    = ( mzwopt(isyma) .gt. 0 )
!           nzyvec = mzyvar(isymc)
!           nzcvec = mzconf(isymc)
!
!           call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymc, etrs,
!     &               vecc, nzyvec, nzcvec, ovlap, isymdn, den1, fcas3_1,
!     &               xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
!        end if
!        deallocate(fcas3_1)
!
!        allocate(fxpec(norbt,norbt))
!        fxpec = 0.0d0
!        call oith1(isymc,zymc,fupe,fxpec,1)
!
!        ! fcas3_2 = Fg(2k) + Fa[1]
!        allocate(fcas3_2(norbt,norbt))
!        if (pe_polar) then
!            fcas3_2 = fxpec + fcas2_2
!        else
!            fcas3_2 = fxpec
!        end if
!
!        if (.not. tdhf) then
!
!           if (mzconf(isymb) .le. 0) return
!
!          !/   <0| [qj,TD2] |01R>  + <01L| [qj,TD2] |0>  \
!          !|   <j| TD2 |01R>                             |
!          !|   <0| [qj+,TD2] |01R> + <01L| [qj+,TD2] |0> |
!          !\  -<01L| TD2 |j>                             /
!
!          ! 2a. Construct the density matrix <01L|..|0> + <0|..|01R>
!           ilsym  = irefsy
!           irsym  = muld2h(irefsy,isymb)
!           ncl    = mzconf(1)
!           ncr    = mzconf(isymb)
!           kzvarl = mzconf(1)
!           kzvarr = mzyvar(isymb)
!
!           den1 = 0.0d0
!           call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
!     &               cref, vecb, ovlap, den1, dummy, 0 ,0, .true.,
!     &               .true., xindx, wrk, 1, lfree, .true.)
!
!          ! 2b. Make the gradient
!           isymdn = muld2h(ilsym,irsym)
!           isymst = muld2h(isyma,irefsy)
!           if ( isymst .eq. irefsy ) then
!              lcon = ( mzconf(isyma) .gt. 1 )
!           else
!              lcon = ( mzconf(isyma) .gt. 0 )
!           end if
!           lorb    = ( mzwopt(isyma) .gt. 0 )
!           nzyvec = mzyvar(isymb)
!           nzcvec = mzconf(isymb)
!
!           call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymb, etrs,
!     &               vecb, nzyvec, nzcvec, ovlap, isymdn, den1, fcas3_2,
!     &               xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
!        end if
!        deallocate(fcas3_2)
!!       !-----------------------------------------------------------
!!       !case 4
!!       !-----------------------------------------------------------
!
!!       !fx2pe = 0.5*Fg(1k,2k) + 0.5*Fg(2k,1k) + ...
!        allocate(fx2pe(norbt,norbt))
!        fx2pe = 0.0d0
!        if (.not. tdhf) then
!           call oith1(isymc, zymc, fxpeb, fx2pe, isymb)
!           call oith1(isymb, zymb, fxpec, fx2pe, isymc)
!           deallocate(fxpeb,fxpec)
!        end if
!        fx2pe = 0.5d0 * fx2pe
!
!        ! ... + fcas2_1(2k) + fcas2_2(k1)
!        if (pe_polar) then
!            call oith1(isymc,zymc,fcas2_1, fx2pe,  isymb)
!            call oith1(isymb,zymb,fcas2_2, fx2pe, isymc)
!            deallocate(fcas2_1,fcas2_2)
!
!            if (.not. tdhf ) then
!            ! + ( S(1)S*(2) + S(2)S*(1) ) * fxo  + ...
!               if ((isymb. eq. isymc) .and. (mzconf(isymb) .gt. 0)) then
!                  allocate(fxo(norbt,norbt))
!                  fcmo = 0.0d0
!                  fxo = 0.0d0
!                 call uthu(0.25d0*fcaos(7*nnbasx+1:8*nnbasx), fcmo, cmo,
!     &                      wrk, nbast, norbt)
!                  call dsptsi(norbt, fcmo, fxo)
!                  nzconf = mzconf(isymb)
!                  nzvar  = mzvar(isymb)
!                  fact   = ddot(nzconf, vecb, 1, vecc(nzvar+1), 1) +
!     &                     ddot(nzconf, vecc, 1, vecb(nzvar+1), 1)
!                           call daxpy(n2orbx, fact, fxo, 1, fx2pe, 1)
!                  deallocate(fxo)
!               end if
!            end if
!
!            allocate(fxo1k2k(norbt,norbt), fxo2k1k(norbt,norbt))
!            fxo1k2k = 0.0d0
!            fxo2k1k = 0.0d0
!            if (.not. tdhf) then
!               allocate(fxc1s2s(norbt,norbt))
!               allocate(fxc1k2s(norbt,norbt), fxc1s2k(norbt,norbt))
!               fxc1s2s = 0.0d0
!               fxc1s2k = 0.0d0
!               fxc1k2s = 0.0d0
!            end if
!
!            ! ... + fxo(1k,2k) + fxo(2k,1k)
!            if (mzwopt(isymb).gt.0 .and. mzwopt(isymc).gt.0) then
!               !fxo(1k,2k) = <0|Fe(1k,2k)|0>Fe
!               fcmo = 0.0d0
!               call uthu(1.0d0*fcaos(nnbasx+1:2*nnbasx),
!     &                   fcmo, cmo, wrk, nbast, norbt)
!               call dsptsi(norbt, fcmo, fxo1k2k)
!
!               !fxo(2k,1k) = R*<0|Fe(2k,1k)|0>Fe
!               fcmo    = 0.0d0
!               call uthu(1.0d0*fcaos(3*nnbasx+1:4*nnbasx),
!     &                   fcmo, cmo, wrk, nbast, norbt)
!               call dsptsi(norbt, fcmo, fxo2k1k)
!            end if
!
!            ! + ... fxc(1s2s) + fxc(2s1s) + 2fxc(1s2k) + 2fxc(1k2S)
!            if (.not. tdhf) then
!               ! fxc(1s2s) = R* ( <01L|..|02R> + <02L|..|01R> )Fe
!               fcmo    = 0.0d0
!               call uthu(0.25d0*fcaos(6*nnbasx+1:7*nnbasx), fcmo, cmo,
!     &                   wrk, nbast, norbt)
!               call dsptsi(norbt, fcmo, fxc1s2s)
!
!               ! + ... fxc(1s2k) + fxc(1k2S)
!               ! fxc(1s2k) = ( <01L|[k2,Epq]|0> + <0|[k2,Epq]|01R> )Fe
!               call uthu(2.0d0*fcaos(8*nnbasx+1:9*nnbasx), fcmo, cmo,
!     &                   wrk, nbast, norbt)
!               call dsptsi(norbt, fcmo, fxc1s2k)
!
!               ! fxc(1k2s) = ( <02L|[k1,Epq]|0> + <0|[k1,Epq]|02R> )Fe
!               call uthu(2.0d0*fcaos(9*nnbasx+1:10*nnbasx), fcmo, cmo,
!     &                   wrk, nbast, norbt)
!               call dsptsi(norbt, fcmo, fxc1k2s)
!            end if
!            deallocate(fcmo)
!
!            if (.not. tdhf) then
!               fx2pe = fx2pe + fxo1k2k + fxo2k1k
!     &               + fxc1s2s + fxc1s2k + fxc1k2s
!               deallocate(fxc1s2s, fxc1s2k, fxc1k2s)
!            else
!               fx2pe = fx2pe + fxo1k2k + fxo2k1k
!            end if
!            deallocate(fxo1k2k, fxo2k1k)
!
!       end if
!
!       !/ <0| [qj ,TE] |0> \
!       !| <j| TE |0>       |
!       !| <0| [qj+,TE] |0> |
!       !\ -<0| TE |j>      /
!
!        isymdn = 1
!        ovlap  = 1.0d0
!        isymst = muld2h(isyma, irefsy)
!        if ( isymst .eq. irefsy ) then
!           lcon = ( mzconf(isyma) .gt. 1 )
!        else
!           lcon = ( mzconf(isyma) .gt. 0 )
!        end if
!        lorb   = ( mzwopt(isyma) .gt. 0 )
!        nzyvec = mzconf(1)
!        nzcvec = mzconf(1)
!
!        call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,
!     &              cref, nzyvec, nzcvec, ovlap, isymdn, udv, fx2pe,
!     &              xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
!        deallocate(cref)
!        deallocate(fupe)
!
!        call qexit('pe_rspmcqr')
!
!        end subroutine pe_rspmcqr
!
!!------------------------------------------------------------------------------
!
!      subroutine pe_rspqro(vecb, vecc, etrs, xindx, zymb, zymc,
!     &                 udv, wrk, nwrk, kzyva, kzyvb, kzyvc,
!     &                 isyma, isymb, isymc, cmo, mjwop)
!
!         use polarizable_embedding, only: pe_master
!
!         implicit none
!
!#include "inforb.h"
!#include "infvar.h"
!#include "infrsp.h"
!#include "infdim.h"
!#include "qrinf.h"
!
!         integer :: kzyva, kzyvb, kzyvc
!         integer :: isyma, isymb, isymc
!         integer :: nwrk
!         real*8, dimension(nwrk) :: wrk
!         real*8, dimension(kzyva) :: etrs
!         real*8, dimension(kzyvb) :: vecb
!         real*8, dimension(kzyvc) :: vecc
!         real*8, dimension(ncmot) :: cmo
!         real*8, dimension(norbt,norbt) :: zymb, zymc
!         real*8, dimension(nashdi,nashdi) :: udv
!         real*8, dimension(lcindx) :: xindx
!         integer, dimension(2,maxwop,8) :: mjwop
!
!         integer :: i, j, k
!         integer :: idum = 1
!         real*8, dimension(:), allocatable :: udcao, ufcmo
!         real*8, dimension(:), allocatable :: dcaos, fcaos
!         real*8, dimension(:), allocatable :: fcmo
!
!         call qenter('pe_rspqro')
!
!         call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
!         call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)
!
!         allocate(udcao(n2basx))
!         allocate(ufcmo(n2orbx))
!         allocate(dcaos(4*nnbasx))
!         dcaos = 0.0d0
!
!         !  D(1k)
!         udcao = 0.0d0
!         call cdens1(isymb, cmo, zymb, udcao, wrk, nwrk)
!         call dgefsp(nbast, udcao, dcaos(1:nnbasx))
!
!         ! D(1k,2k)
!         udcao = 0.0d0
!         call cdens2(isymb, isymc, cmo, zymb, zymc, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))
!
!         !  D(2k)
!         udcao = 0.0d0
!         call cdens1(isymc, cmo, zymc, udcao, wrk, nwrk)
!         call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))
!
!         !  D(2k,1k)
!         udcao = 0.0d0
!         call cdens2(isymc, isymb, cmo, zymc, zymb, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))
!
!         deallocate(udcao)
!
!         allocate(fcaos(4*nnbasx))
!         fcaos = 0.0d0
!         call pe_master(runtype='response', denmats=dcaos,
!     &                  fckmats=fcaos, nmats=4)
!         deallocate(dcaos)
!
!         allocate(fcmo(nnorbx))
!         ufcmo = 0.0d0
!
!         i = 1
!         j = nnbasx
!         call uthu(2.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymc, zymc, wrk(1:n2orbx), ufcmo, isyma)
!
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         ufcmo = ufcmo + wrk(1:n2orbx)
!
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(2.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymb, zymb, wrk(1:n2orbx), ufcmo, isyma)
!
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         ufcmo = ufcmo + wrk(1:n2orbx)
!
!         call rsp1gr(1, kzyva, idum, 0, isyma, 0, 1, etrs,
!     &               wrk, idum, idum, 1.0d0, 1, udv, ufcmo, xindx,
!     &               mjwop, wrk, nwrk, .true., .false., .false.)
!
!         deallocate(fcaos, fcmo, ufcmo)
!
!        call qexit('pe_rspqro')
!
!      end subroutine pe_rspqro
!
!!------------------------------------------------------------------------------
!
!      subroutine pe_rspcro(vecb, vecc, vecd, etrs, xindx, zymb, zymc,
!     &                     zymd, udv, wrk, nwrk, kzyva, kzyvb, kzyvc,
!     &                     kzyvd, isyma, isymb, isymc, isymd, cmo,mjwop)
!
!         use polarizable_embedding, only: pe_master
!
!         implicit none
!
!#include "inforb.h"
!#include "infvar.h"
!#include "infrsp.h"
!#include "infdim.h"
!#include "qrinf.h"
!
!         integer :: kzyva, kzyvb, kzyvc, kzyvd
!         integer :: isyma, isymb, isymc, isymd
!         integer :: nwrk
!         real*8, dimension(nwrk) :: wrk
!         real*8, dimension(kzyva) :: etrs
!         real*8, dimension(kzyvb) :: vecb
!         real*8, dimension(kzyvc) :: vecc
!         real*8, dimension(kzyvd) :: vecd
!         real*8, dimension(ncmot) :: cmo
!         real*8, dimension(norbt,norbt) :: zymb, zymc, zymd
!         real*8, dimension(nashdi,nashdi) :: udv
!         real*8, dimension(lcindx) :: xindx
!         integer, dimension(2,maxwop,8) :: mjwop
!
!         integer :: i, j, k
!         integer :: idum = 1
!         real*8, dimension(:), allocatable :: udcao, ufcmo
!         real*8, dimension(:), allocatable :: dcaos, fcaos
!         real*8, dimension(:), allocatable :: fcmo
!
!         call qenter('pe_rspcro')
!
!         call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
!         call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)
!         call gtzymt(1, vecd, kzyvd, isymd, zymd, mjwop)
!
!         allocate(udcao(n2basx))
!         allocate(ufcmo(n2orbx))
!         allocate(dcaos(15*nnbasx))
!         dcaos = 0.0d0
!
!         udcao = 0.0d0
!         call cdens1(isymb, cmo, zymb, udcao, wrk, nwrk)
!         call dgefsp(nbast, udcao, dcaos(1:nnbasx))
!         udcao = 0.0d0
!         call cdens1(isymc, cmo, zymc, udcao, wrk, nwrk)
!         call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))
!         udcao = 0.0d0
!         call cdens1(isymd, cmo, zymd, udcao, wrk, nwrk)
!         call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))
!
!         udcao = 0.0d0
!         call cdens2(isymb, isymc, cmo, zymb, zymc, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))
!         udcao = 0.0d0
!         call cdens2(isymc, isymb, cmo, zymc, zymb, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(4*nnbasx+1:5*nnbasx))
!         udcao = 0.0d0
!         call cdens2(isymb, isymd, cmo, zymb, zymd, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(5*nnbasx+1:6*nnbasx))
!         udcao = 0.0d0
!         call cdens2(isymd, isymb, cmo, zymd, zymb, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(6*nnbasx+1:7*nnbasx))
!         udcao = 0.0d0
!         call cdens2(isymc, isymd, cmo, zymc, zymd, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(7*nnbasx+1:8*nnbasx))
!         udcao = 0.0d0
!         call cdens2(isymd, isymc, cmo, zymd, zymc, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(8*nnbasx+1:9*nnbasx))
!
!         udcao = 0.0d0
!         call cdens3(isymb, isymc, isymd, cmo, zymb, zymc, zymd, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(9*nnbasx+1:10*nnbasx))
!         udcao = 0.0d0
!         call cdens3(isymd, isymb, isymc, cmo, zymd, zymb, zymc, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(10*nnbasx+1:11*nnbasx))
!         udcao = 0.0d0
!         call cdens3(isymc, isymd, isymb, cmo, zymc, zymd, zymb, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(11*nnbasx+1:12*nnbasx))
!         udcao = 0.0d0
!         call cdens3(isymb, isymd, isymc, cmo, zymb, zymd, zymc, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(12*nnbasx+1:13*nnbasx))
!         udcao = 0.0d0
!         call cdens3(isymc, isymb, isymd, cmo, zymc, zymb, zymd, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(13*nnbasx+1:14*nnbasx))
!         udcao = 0.0d0
!         call cdens3(isymd, isymc, isymb, cmo, zymd, zymc, zymb, udcao,
!     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
!         call dgefsp(nbast, udcao, dcaos(14*nnbasx+1:15*nnbasx))
!
!         deallocate(udcao)
!
!         allocate(fcaos(15*nnbasx))
!         fcaos = 0.0d0
!         call pe_master(runtype='response', denmats=dcaos,
!     &                  fckmats=fcaos, nmats=15)
!         deallocate(dcaos)
!
!         allocate(fcmo(nnorbx))
!         ufcmo = 0.0d0
!
!         i = 1
!         j = nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:2*n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymc, zymc, wrk(1:n2orbx),
!     &              wrk(n2orbx+1:2*n2orbx), isyma)
!         call oith1(isymd, zymd, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:2*n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymb, zymb, wrk(1:n2orbx),
!     &              wrk(n2orbx+1:2*n2orbx), isyma)
!         call oith1(isymd, zymd, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:2*n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymb, zymb, wrk(1:n2orbx),
!     &              wrk(n2orbx+1:2*n2orbx), isyma)
!         call oith1(isymc, zymc, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
!         i = 1
!         j = nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:2*n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymd, zymd, wrk(1:n2orbx),
!     &              wrk(n2orbx+1:2*n2orbx), isyma)
!         call oith1(isymc, zymc, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:2*n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymd, zymd, wrk(1:n2orbx),
!     &              wrk(n2orbx+1:2*n2orbx), isyma)
!         call oith1(isymb, zymb, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:2*n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymc, zymc, wrk(1:n2orbx),
!     &              wrk(n2orbx+1:2*n2orbx), isyma)
!         call oith1(isymb, zymb, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
!
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymd, zymd, wrk(1:n2orbx), ufcmo, isyma)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymd, zymd, wrk(1:n2orbx), ufcmo, isyma)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymc, zymc, wrk(1:n2orbx), ufcmo, isyma)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymc, zymc, wrk(1:n2orbx), ufcmo, isyma)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymb, zymb, wrk(1:n2orbx), ufcmo, isyma)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         call oith1(isymb, zymb, wrk(1:n2orbx), ufcmo, isyma)
!
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         ufcmo = ufcmo + wrk(1:n2orbx)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         ufcmo = ufcmo + wrk(1:n2orbx)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         ufcmo = ufcmo + wrk(1:n2orbx)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         ufcmo = ufcmo + wrk(1:n2orbx)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         ufcmo = ufcmo + wrk(1:n2orbx)
!         i = i + nnbasx
!         j = j + nnbasx
!         call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
!         wrk(1:n2orbx) = 0.0d0
!         call dsptsi(norbt, fcmo, wrk(1:n2orbx))
!         ufcmo = ufcmo + wrk(1:n2orbx)
!
!         call rsp1gr(1, kzyva, idum, 0, isyma, 0, 1, etrs,
!     &               wrk, idum, idum, 1.0d0, 1, udv, ufcmo, xindx,
!     &               mjwop, wrk, nwrk, .true., .false., .false.)
!
!         deallocate(fcaos, fcmo, ufcmo)
!
!        call qexit('pe_rspcro')
!
!      end subroutine pe_rspcro
!
!!------------------------------------------------------------------------------
!
!      subroutine dalton_cavity(all_centers, all_charges, ncents,
!     &                         wrk, nwrk)
!
!         use mxcent_mod
!!        implicit none
!
!#include "implicit.h"
!#include "dummy.h"
!#include "priunit.h"
!#include "pcmdef.h"
!#include "pcm.h"
!#include "pcmlog.h"
!
!        integer :: ncents, nwrk
!        real*8, dimension(ncents) :: all_charges
!        real*8, dimension(3,ncents) :: all_centers
!        real*8, dimension(nwrk) :: wrk
!
!        integer :: i
!        integer :: ntess, lusurf, lucav
!        integer :: bck_nucind, bck_nucdep, bck_natoms
!        integer, dimension(:), allocatable :: isphere
!        real*8, dimension(:), allocatable :: areas, rvdw
!        real*8, dimension(:,:), allocatable :: tess_coords
!
!        call qenter('dalton_cavity')
!
!C       van der Waals radii taken from "The Elements", 2nd edition,
!C       John Emsley, Clarendon Press, Oxford, 1991.  Unknown values
!C       are simply set to D0, rather than trying to guess them.
!C
!C       A.Bondi, J.Phys.Chem. 68: 441-451(1964) gives alternate
!C       values, and a few transition metals.
!        allocate(rvdw(99))
!        rvdw = (/ 1.20d0, 1.22d0, 0.00d0, 0.00d0, 2.08d0, 1.85d0,
!     &            1.54d0, 1.40d0, 1.35d0, 1.60d0, 2.31d0, 0.00d0,
!     &            2.05d0, 2.00d0, 1.90d0, 1.85d0, 1.81d0, 1.91d0,
!     &            2.31d0, 13*0.0d0, 2.00d0, 2.00d0, 1.95d0, 1.98d0,
!     &            2.44d0, 13*0.0d0, 2.20d0, 2.20d0, 2.15d0, 0.00d0,
!     &            2.62d0, 27*0.0d0, 2.40d0, 16*0.0d0 /)
!C     override the above table with U. Pisa's experience
!C     as to what works best for singly bonded C,N,O
!        rvdw(6) = 1.70d0
!        rvdw(7) = 1.60d0
!        rvdw(8) = 1.50d0
!
!        if (ncents > mxsp) then
!            stop 'ERROR: ncents > mxsp'
!        end if
!
!        nesfp = ncents
!        nesf = nesfp
!
!        do i = 1, nesfp
!            xe(i) = all_centers(1,i)
!            ye(i) = all_centers(2,i)
!            ze(i) = all_centers(3,i)
!            if (nint(all_charges(i)) > 99) then
!                stop 'ERROR: charge is too big'
!            end if
!            rin(i) = rvdw(nint(all_charges(i)))
!            if (rin(i) == 0.0d0) then
!                stop 'ERROR: no vdw radius defined'
!            end if
!        end do
!
!        bck_nucind = nucind
!        bck_nucdep = nucdep
!        bck_natoms = natoms
!
!        nucind = ncents
!        nucdep = ncents
!        natoms = ncents
!
!        omega = 40.0d0
!        fro = 0.7d0
!        ret = 0.2d0
!        nrwcav = 2
!        areats = 0.4d0
!        oldcen = .false.
!        iprpcm = 0
!        dr = 1.0d-4
!        alpha = 1.2d0
!        rsolv = 1.385d0
!
!        call pedram(wrk, nwrk)
!
!        nucind = bck_nucind
!        nucdep = bck_nucdep
!        natoms = bck_natoms
!
!        lucav = -999
!        call gpopen(lucav, 'CAVDATA', 'OLD', 'SEQUENTIAL',
!     &              'UNFORMATTED', idummy, .FALSE.)
!        read(lucav) ntess
!        allocate(isphere(ntess))
!        allocate(tess_coords(3,ntess))
!        allocate(areas(ntess))
!        read(lucav) isphere
!        read(lucav) tess_coords(1,:)
!        read(lucav) tess_coords(2,:)
!        read(lucav) tess_coords(3,:)
!        read(lucav) areas
!        call gpclose(lucav, 'KEEP')
!
!        lusurf = -999
!        call gpopen(lusurf, 'SURFACE.INP', 'NEW', 'SEQUENTIAL',
!     &              'FORMATTED', idummy, .FALSE.)
!
!        write(lusurf,'(i6)') ntess
!        write(lusurf,'(a)') 'AA'
!        do i = 1, ntess
!            write(lusurf,'(4f12.6)') tess_coords(:,i), areas(i)
!        end do
!
!        call gpclose(lusurf,'KEEP')
!
!        call qexit('dalton_cavity')
!
!      end subroutine dalton_cavity
!
!!------------------------------------------------------------------------------
!
!      subroutine Tk_lao_integrals(Tk_ints, nints, ncomps, coord)
!
!      use pe_work
!      use mxcent_mod
!#include "implicit.h"
!#include "dummy.h"
!#include "qm3.h"
!#include "orgcom.h"
!
!      logical :: trimat
!      character(len=8), dimension(9*mxcent) :: labint
!      integer, dimension(9*mxcent) :: intrep, intadr
!
!      integer, intent(in) :: nints, ncomps
!      real(8), dimension(3), intent(in) :: coord
!      real(8), dimension(3*ncomps*nints), intent(out) :: Tk_ints
!
!      integer :: i, j, k, l, m, n
!      real(8) :: dd
!      real(8), dimension(3) :: backup
!      character(len=7) :: inttype
!      integer :: nwrk
!      real(8), dimension(:), allocatable :: temp
!
!      call qenter('Tk_lao_integrals')
!
!      dd = 0.001d0
!      k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * real(ncomps)) - 1.0d0)) - 1
!
!      backup = diporg
!      diporg = coord
!      runqm3 = .true.
!
!      ! magnetic properties are assymetric
!      trimat = .false.
!
!      if (k == 0) then
!          inttype = 'PCMBSOL'
!      else if (k == 1) then
!          inttype = 'EFIELB1'
!      else if (k == 2) then
!        !inttype = 'EFIELB1'
!        continue
!      else if (k >= 3) then
!          stop 'ERROR: electric field hessian and higher order integrals
!     & not implemented for london orbitals.'
!!        if (nwrk < 24 * nints) then
!!            print *, 'Not enough work space for T^(3) integrals!'
!!        end if
!!        inttype = 'ELFGRDC'
!!    else
!!        stop 'wrong order specified or not implemented'
!      end if
!
!      if (.not. associated(pewrk)) then
!        call quit('ERROR: pewrk not allocated')
!      endif
!
!      nwrk = size(pewrk)
!
!      Tk_ints = 0.0d0
!
!      if (k == 0) then
!          n = ncomps * 3
!          call get1in(Tk_ints(1), inttype, n, pewrk(1), nwrk,
!     &                labint, intrep, intadr, 0, .false., 0, trimat,
!     &                dummy, .false., dummy, 1)
!          Tk_ints = - Tk_ints
!      else if (k == 1) then
!!          n = 3
!!          m = nints*n
!!      ! numerical approach
!!      ! x first
!!          do j=1,3
!!          diporg(j) = coord(j) + dd
!!          call get1in(temp(1:m), 'PCMBSOL', n, pewrk(1), nwrk, labint,
!!     &                intrep, intadr, 0, .false., 0, trimat, dummy,
!!     &                .false., dummy, 1)
!!          diporg(j) = coord(j) - dd
!!          call get1in(temp(m+1:2*m), 'PCMBSOL', n, pewrk(1),
!!     &    nwrk, labint,
!!     &                intrep, intadr, 0, .false., 0, trimat, dummy,
!!     &                .false., dummy, 1)
!!         TK_ints((j-1)*m+1:j*m) = (temp(1:m) - temp(m+1:m*2))/(2.0d0*dd)
!!          diporg = coord
!!          enddo
!
!      ! analytical approach
!          n = 9
!          call get1in(Tk_ints(1), inttype, n, pewrk(1), nwrk,
!     &                labint, intrep, intadr, 0, .false., 0, trimat,
!     &                dummy, .false., dummy, 1)
!      else if (k == 2) then
!          allocate(temp(2*9*nints))
!          n = 9
!          m = nints * n
!          ! xx, xy, xz
!          diporg = coord
!          diporg(1) = coord(1) + dd
!          call get1in(temp(1:m), 'EFIELB1', n, pewrk(1), nwrk,
!     &                labint, intrep, intadr, 0, .false., 0, trimat,
!     &                dummy, .false., dummy, 1)
!          diporg(1) = coord(1) - dd
!          call get1in(temp(1+m:2*m), 'EFIELB1', n, pewrk(1), nwrk,
!     &                labint, intrep, intadr, 0, .false., 0, trimat,
!     &                dummy, .false., dummy, 1)
!          Tk_ints(1:m) = temp(1:m) - temp(m+1:2*m)
!
!          ! yy, yz
!          diporg = coord
!          diporg(2) = coord(2) + dd
!          call get1in(temp(1:m), 'EFIELB1', n, pewrk(1), nwrk,
!     &                labint, intrep, intadr, 0, .false., 0, trimat,
!     &                dummy, .false., dummy, 1)
!          diporg(2) = coord(2) - dd
!          call get1in(temp(1+m:2*m), 'EFIELB1', n, pewrk(1), nwrk,
!     &                labint, intrep, intadr, 0, .false., 0, trimat,
!     &                dummy, .false., dummy, 1)
!          i = 3 * nints
!          Tk_ints(m+1:m+2*i+1) = temp(i+1:m) - temp(m+i+1:2*m)
!
!          ! zz
!          diporg = coord
!          diporg(3) = coord(3) + dd
!          call get1in(temp(1:m), 'EFIELB1', n, pewrk(1), nwrk,
!     &                labint, intrep, intadr, 0, .false., 0, trimat,
!     &                dummy, .false., dummy, 1)
!          diporg(3) = coord(3) - dd
!          call get1in(temp(1+m:2*m), 'EFIELB1', n, pewrk(1), nwrk,
!     &                labint, intrep, intadr, 0, .false., 0, trimat,
!     &                dummy, .false., dummy, 1)
!
!          i = 6 * nints
!          Tk_ints(m+i+1:2*m) = temp(i+1:m) - temp(m+i+1:2*m)
!          ! suspecting sign from magnus
!          Tk_ints = - Tk_ints / (2.0d0 * dd)
!          deallocate(temp)
!!      else if (k == 3) then
!!          n = 6
!!          m = nints
!!          i = n * nints
!!          l = 0
!!          do j = 1, 3
!!              diporg(j) = diporg(j) + 0.01d0
!!              call get1in(temp(l*i+1), inttype, n, pewrk(1), nwrk,
!!     &                    labint, intrep, intadr, 0, .false., 0, trimat,
!!     &                    dummy, .false., dummy, 1)
!!              diporg(j) = diporg(j) - 2.0d0 * 0.01d0
!!              call get1in(temp(j*i+1), inttype, n, pewrk(1),
!!     &                    nwrk, labint, intrep, intadr, 0, .false., 0,
!!     &                    trimat, dummy, .false., dummy, 1)
!!              diporg(j) = coord(j)
!!              temp(l*i+1:j*i) = (temp(l*i+1:j*i) - temp(j*i+1:(j+1)*i))
!!     &                                           / (2.0d0 * 0.01d0)
!!              l = l + 1
!!          end do
!!          Tk_ints(1:m) = temp(1:m)
!!          Tk_ints(m+1:2*m) = temp(i+1:i+1+m)
!!          Tk_ints(2*m+1:3*m) = temp(2*i+1:2*i+1+m)
!!          Tk_ints(3*m+1:4*m) = temp(i+1+m:i+1+2*m)
!!          Tk_ints(4*m+1:5*m) = temp(2*i+1+m:2*i+1+2*m)
!!          Tk_ints(5*m+1:6*m) = temp(2*i+1+2*m:2*i+1+3*m)
!!          Tk_ints(6*m+1:7*m) = temp(i+1+3*m:i+1+4*m)
!!          Tk_ints(7*m+1:8*m) = temp(2*i+1+3*m:2*i+1+4*m)
!!          Tk_ints(8*m+1:9*m) = temp(2*i+1+4*m:2*i+1+5*m)
!!          Tk_ints(9*m+1:10*m) = temp(2*i+1+5*m:2*i+1+6*m)
!!          Tk_ints = -1.0d0 * Tk_ints
!      end if
!
!      runqm3 = .false.
!      diporg = backup
!
!      call qexit('Tk_lao_integrals')
!
!      end subroutine Tk_lao_integrals

end module pelib_dirac_interface

#endif

