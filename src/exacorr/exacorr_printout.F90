module exacorr_printout

  implicit none
  private

  public print_exacorr_logo
  public print_exaoutput
  public print_exainput
  public print_orbital_energy
  public print_iteration
  public print_iterative_wavefunction_choice
  public print_mp2fvno_occupation_number_details
 
  contains

!  Print out routines used in exacorr

    subroutine print_exacorr_logo()
    ! Prints logo and authors of exacorr

      print*, ""
      print*, ""
      print*, "                   *****************************************************"
      print*, "                   ***     Entering the ExaCorr module in DIRAC      ***"
      print*, "                   ***                                               ***"
      print*, "                   ***  Authors (in alphabetical order) :            ***"
      print*, "                   ***                   - Jan Brandejs              ***"
      print*, "                   ***                   - Chima Chibueze            ***"
      print*, "                   ***                   - Gabriele Fabbro           ***"
      print*, "                   ***                   - Loic Halbert              ***"
      print*, "                   ***                   - Thomas Lejeune            ***"
      print*, "                   ***                   - Anastasios Papadopoulos   ***"
      print*, "                   ***                   - Johann Pototschnig        ***"
      print*, "                   ***                   - Michal Repisky            ***"
      print*, "                   ***                   - Andre Severo Pereira Gomes***"
      print*, "                   ***                   - Lucas Visscher            ***"
      print*, "                   ***                   - Xiang Yuan                ***"
      print*, "                   ***                                               ***"
      print*, "                   ***  Features:        - CCD/CC2/CCSD              ***"
      print*, "                   ***                   - Higher-order CC methods   ***"
      print*, "                   ***                   - Frozen virtual NOs        ***"
      print*, "                   ***                   - 1RDM CCD/CCSD, 2RDM CCSD  ***"
      print*, "                   ***                   - CIS and CIS(D)            ***"
      print*, "                   ***                   - EOM-CCSD                  ***"
      print*, "                   ***                   - LR and QR CCSD            ***"
      print*, "                   ***                   - ReSpect interface         ***"
      print*, "                   ***                                               ***"
      print*, "                   ***  Relativistic Coupled Cluster code using      ***"
      print*, "                   ***  CPU and GPU-accelerated tensor libraries :   ***"
      print*, "                   ***                                               ***" 
      print*, "                   ***  - TALSH and ExaTensor libraries developed by ***"
      print*, "                   ***  Dmitry Lyakh (see doi: 10.1002/qua.25926)    ***"
      print*, "                   ***  - Backends supported by the Tensor Agebra    ***"
      print*, "                   ***  Processing Primitives (TAPP) API (see TAPP   ***"
      print*, "                   ***  whitepaper https://arxiv.org/abs/2601.07827) ***"
      print*, "                   ***                                               ***"
      print*, "                   *****************************************************"
      print*, ""
      print*, "References"
      print*, ""
      print*, " - Ground-state CCD/CC2/CCSD/CCSD(T) energies, CCD/CCSD expectation values:"
      print*, "    JV Pototschnig et al., J. Chem. Theory Comput. 17, 5509 (2021)." 
      print*, "    doi: 10.1021/acs.jctc.1c00260"
      print*, ""
      print*, " - MP2 virtual frozen natural orbitals, generation of CC natural orbitals:"
      print*, "    X Yuan, L Visscher, ASP Gomes, J. Chem. Phys. 156, 224108 (2022). "
      print*, "    doi: 10.1063/5.0087243"
      print*, ""
      print*, " - Virtual frozen natural orbitals for response properties:"
      print*, "    M Le, X Yuan, ASP Gomes, in preparation"
      print*, ""
      print*, " - Equation of motion coupled cluster (EOM-CCSD) for IP, EE and EA:"
      print*, "    X Yuan, PhD thesis, Vrije Universiteit Amsterdam and Universite de Lille"
      print*, "    doi: 10.5463/thesis.505"
      print*, ""
      print*, " - (EOM-)CCSD 1RDMs and CCSD 2RDM"
      print*, "    L Halbert, ASP Gomes, in preparation"
      print*, ""
      print*, " - CIS and CIS(D) excitation energies and transition moments"
      print*, "    X Yuan, T Lejeune, ASP Gomes, in preparation"
      print*, ""
      print*, " - Orbital-unrelaxed linear response (EOM-)CCSD theory:"
      print*, "    X Yuan et al., J. Chem. Theory Comput. (2024). "
      print*, "    doi: 10.1021/acs.jctc.3c00812"
      print*, ""
      print*, " - Orbital-unrelaxed quadratic response (EOM-)CCSD theory:"
      print*, "    X Yuan, L Halbert, L Visscher, ASP Gomes, J. Chem. Theory Comput. 19, 9248 (2023). "
      print*, "    doi: 10.1021/acs.jctc.3c01011"
      print*, "    X Yuan, L Halbert, L Visscher, ASP Gomes, J. Phys. Chem. A 129, 11695 (2025). "
      print*, "    doi: 10.1021/acs.jpca.5c04475"
      print*, ""
      print*, " - Higher-order coupled cluster methods"
      print*, "  - code generated with the Tenpi framework (https://gitlab.com/janbrandejs/tenpi)"
      print*, "    - ground-state energies (CCSDT, CCSDTQ) " 
      print*, "      J Brandejs, J Pototschnig, T Saue, J. Chem. Theory Comput. 21, 7320 (2025). "
      print*, "      doi: 10.1021/acs.jctc.5c00219"
      print*, "    - ground-state expectation values (CCSDT, CCSDTQ) " 
      print*, "      G Fabbro, J Brandejs, T Saue, J. Phys. Chem. A 129, 6942 (2025). "
      print*, "      doi: 10.1021/acs.jpca.5c02844"
      print*, "  - code generated with the mb-autogen framework (doi: 10.5281/zenodo.8094429)"
      print*, "    - ground-state energies (CCSDT/-1a/-1b/-2/-3, CCSDTQ, CCSDTQP) " 
      print*, "      L Halbert, ASP Gomes, in preparation"
      print*, ""

    end subroutine print_exacorr_logo

    subroutine print_exaoutput(talsh, iterative_scheme, perturbative_scheme,scf_energy, mp2_energy, cc_energy, t1diag, pert_energy)

      logical, intent(in)           :: talsh
      integer, intent(in)           :: iterative_scheme
      integer, intent(in)           :: perturbative_scheme
      real(8), intent(in)           :: scf_energy, mp2_energy, cc_energy
      real(8), intent(in), optional :: t1diag
      real(8), intent(in), optional :: pert_energy(3)
      character(len=23) :: label_iterative

!     character(len=12) :: wf_label
!      
!     call print_iterative_wavefunction_choice(iterative_scheme, wf_label)
      write(label_iterative,'(A)') "                    "
      if (iterative_scheme.eq.20) write(label_iterative,'(A)') "CCD energy :           "
      if (iterative_scheme.eq.21) write(label_iterative,'(A)') "CCSD energy :          "
      if (iterative_scheme.eq.22) write(label_iterative,'(A)') "CC2 energy :           "

      if (iterative_scheme.eq.30) write(label_iterative,'(A)') "CCSDT energy :         "
      if (iterative_scheme.eq.31) write(label_iterative,'(A)') "CCSDT-1a energy :      "
      if (iterative_scheme.eq.32) write(label_iterative,'(A)') "CCSDT-1b energy :      "
      if (iterative_scheme.eq.33) write(label_iterative,'(A)') "CCSDT-2 energy :       "
      if (iterative_scheme.eq.34) write(label_iterative,'(A)') "CCSDT-3 energy :       "
      if (iterative_scheme.eq.35) write(label_iterative,'(A)') "CC3 energy :           "

      if (iterative_scheme.eq.40) write(label_iterative,'(A)') "CCSDTQ energy :        "
      if (iterative_scheme.eq.41) write(label_iterative,'(A)') "CCSDTQ-1a energy :     "
      if (iterative_scheme.eq.42) write(label_iterative,'(A)') "CCSDTQ-1b energy :     "

      if (iterative_scheme.eq.50) write(label_iterative,'(A)') "CCSDTQP energy :       "

      write(*,*) ""
      write(*,'(A)') "-----------------------------------------------"
      write(*,'(A)') "     -     Final results from EXACORR     -    "
      write(*,'(A)') "-----------------------------------------------"
      if (talsh) then 
        write(*,'(A)') " Using the single-node version "
      else 
        write(*,'(A)') " Using the multi-node version (MPI) "
      end if
      write(*,*) ""
      write (*,'(A,F20.15)') "  MP2 energy :                      ",mp2_energy
      write (*,'(A31,F25.15)') "Final "//label_iterative,cc_energy
      if (present(pert_energy).and.(perturbative_scheme.eq.30)) then 
            write (*,'(A,F20.15)') "  4th order triples correction :    ",pert_energy(1)
            write (*,'(A,F20.15)') "  5th order triples (T) correction :",pert_energy(2)
            write (*,'(A,F20.15)') "  5th order triples -T  correction :",pert_energy(3)
            write (*,'(A,F20.15)') "  CCSD+T  correlation energy :      ",cc_energy+pert_energy(1)
            write (*,'(A,F20.15)') "  CCSD(T) correlation energy :      ",cc_energy+pert_energy(1)+pert_energy(2)
            write (*,'(A,F20.15)') "  CCSD-T  correlation energy :      ",cc_energy+pert_energy(1)+pert_energy(3)
      end if
      if (present(t1diag)) then
         write (*,'(A,F20.15)') "  T1 diagnostic      :              ",t1diag
      end if
      write(*,*) ""
      write (*,'(A,F25.15)') "  SCF energy :                 ",scf_energy
      write (*,'(A,F25.15)') "  Total MP2 energy :           ",scf_energy+mp2_energy
      write (*,'(A31,F25.15)') "Total "//label_iterative,(scf_energy+cc_energy)

      if (present(pert_energy).and.(perturbative_scheme.eq.30)) then
         write (*,'(A,F25.15)') "  Total CCSD+T  energy :       ",scf_energy+cc_energy+pert_energy(1)
         write (*,'(A,F25.15)') "  Total CCSD(T) energy :       ",scf_energy+cc_energy+pert_energy(1)+pert_energy(2)
         write (*,'(A,F25.15)') "  Total CCSD-T  energy :       ",scf_energy+cc_energy+pert_energy(1)+pert_energy(3)
      end if
      write(*,*) ""
      write(*,'(A)') "-----------------------------------------------"
      write(*,*) ""

    end subroutine print_exaoutput

    subroutine print_iteration(iteration,convergence,energy,print_level)
    ! Formatted print for iterations

      ! input variables
      integer, intent(in)           :: iteration
      real(8), intent(in)           :: convergence
      real(8), intent(in), optional :: energy
      integer, intent(in), optional :: print_level

      if (iteration.eq.1 .and. present(energy)) then
        write(*,*) ""
        write(*,*) "-----------------------------------------------"
        write(*,*) "  It.     Energy                  Convergence  "
        write(*,*) "-----------------------------------------------"
        write(*,'(A3,I2,A4,F20.15,A3,E10.1)') "   ", iteration, "    ", energy, "   ", convergence
      else if (iteration.eq.1) then
        write(*,*) ""
        write(*,*) "------------------------"
        write(*,*) "  It.      Convergence  "
        write(*,*) "------------------------"
        write(*,'(A3,I2,A4,E10.1)') "   ", iteration, "    ", convergence
      else if (iteration.ne.1 .and. present(energy)) then
        write(*,'(A3,I2,A4,F20.15,A3,E10.1)') "   ", iteration, "    ", energy, "   ", convergence
      else if (iteration.ne.1) then
        write(*,'(A3,I2,A4,E10.1)') "   ", iteration, "    ", convergence
      end if
      
      if (present(print_level)) then
        if (print_level.gt.4) then
          call print_date(' CC iteration done')
        end if
      end if

    end subroutine print_iteration

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
    subroutine write_tensor(exa_tensor, i_dims)

      use exatensor
      use talsh
      use tensor_algebra

      type(tens_rcrsv_t)              :: exa_tensor
      integer(INTL), dimension(2)     :: i_dims
      integer(INTD), dimension(2)     :: h_dims
      type(talsh_tens_t)              :: talsh_tensor
      complex(8), pointer, contiguous :: talsh_tens(:,:)
      type(C_PTR)                     :: body_p
      integer                         :: i,j
      integer(INTD)                   :: ierr
      integer(C_SIZE_T)               :: buf_size=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
      integer(C_INT)                  :: host_arg_max

      ierr=talsh_init(buf_size,host_arg_max)
      call print_date('Initialized talsh library')
      write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

      h_dims=i_dims
      write(*,*) '+++ write to file'
      ierr=talsh_tensor_construct(talsh_tensor,C8,h_dims,init_val=(0.D0,0.D0))
      write(*,*) 'construct err:', ierr
      ierr=exatns_tensor_get_slice(exa_tensor,talsh_tensor)
      write(*,*) 'get_slice err:', ierr
      ierr=talsh_tensor_get_body_access(talsh_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p,talsh_tens,h_dims) 
      do i=1,h_dims(2)
          write(*,*) (talsh_tens(i,j), j=1,h_dims(1))
      enddo
      ierr=talsh_tensor_destruct(talsh_tensor)

      ierr = talsh_shutdown()
      
    end subroutine write_tensor
#endif

    subroutine print_orbital_energy(eps_occ,nocc,eps_vir,nvir,info)
      implicit none

      real(8), intent(in)   :: eps_vir(:),eps_occ(:)
      integer, intent(in)   :: nocc,nvir,info

      integer :: i
      
      if (info==1) then
        write(*,*) ' --- shifted active occupied spinors --- '
      else
        write(*,*) ' --- recomputed active occupied spinors --- '
      end if 
      write(*,*) '   #             E   '
      do i=1,nocc
        write(*,'(I8,E18.9)') i, eps_occ(i)
      end do

      if (info==1) then
        write(*,*) ' --- shifted active virtual spinors --- '
      else
        write(*,*) ' --- recomputed active virtual spinors --- '
      end if
      write(*,*) '   #             E   '
      do i=1,nvir
        write(*,'(I8,E18.9)') i, eps_vir(i)
      end do

    end subroutine print_orbital_energy

    subroutine print_iterative_wavefunction_choice(iterative_scheme, wf_label)
      integer, intent(in) :: iterative_scheme
      character(12), optional, intent(inout) :: wf_label
      character(12) :: wavefunction

      if (iterative_scheme .eq. 10) wavefunction = "CCS" 

      if (iterative_scheme .eq. 20) wavefunction = "CCD" 
      if (iterative_scheme .eq. 21) wavefunction = "CCSD" 
      if (iterative_scheme .eq. 22) wavefunction = "CC2" 

      if (iterative_scheme .eq. 30) wavefunction = "CCSDT" 
      if (iterative_scheme .eq. 31) wavefunction = "CCSDT-1a" 
      if (iterative_scheme .eq. 32) wavefunction = "CCSDT-1b" 
      if (iterative_scheme .eq. 33) wavefunction = "CCSDT-2" 
      if (iterative_scheme .eq. 34) wavefunction = "CCSDT-3" 
      if (iterative_scheme .eq. 35) wavefunction = "CC3" 

      if (iterative_scheme .eq. 40) wavefunction = "CCSDTQ" 
      if (iterative_scheme .eq. 41) wavefunction = "CCSDTQ-1a" 
      if (iterative_scheme .eq. 42) wavefunction = "CCSDTQ-1b" 
      if (iterative_scheme .eq. 43) wavefunction = "CCSDTQ-2" 
      if (iterative_scheme .eq. 44) wavefunction = "CCSDTQ-3" 
      if (iterative_scheme .eq. 45) wavefunction = "CC4" 

      if (iterative_scheme .eq. 50) wavefunction = "CCSDTQP" 
      if (iterative_scheme .eq. 51) wavefunction = "CCSDTQP-1a" 
      if (iterative_scheme .eq. 52) wavefunction = "CCSDTQP-1b" 
      if (iterative_scheme .eq. 53) wavefunction = "CCSDTQP-2" 
      if (iterative_scheme .eq. 54) wavefunction = "CCSDTQP-3" 
      if (iterative_scheme .eq. 55) wavefunction = "CC5" 

      if (iterative_scheme .eq. 99) wavefunction = "CCSDTQP-xyz " 

      if (.not.present(wf_label)) then
         write(*,'(A)') ""
         write(*,'(1X,A,I4,A,I4,A)') "Iterations with CC(",iterative_scheme,") model: "//wavefunction
      else
         write(wf_label,'(A)') wavefunction
      end if   
   end subroutine

   subroutine print_exainput(exa_input,nao)

      use exacorr_cfg,   only : exacc_input

      type(exacc_input), intent(in) :: exa_input
      integer, intent(in) :: nao
      integer             :: i

      write(*,*) ""
      write(*,'(A)') "-----------------------------------------------"
      write(*,'(A)') "      Characteristics of this calculation      "
      write(*,'(A)') "-----------------------------------------------"
      if (exa_input%auto_kernels) then
         if (exa_input%auto_scheme.eq.2) then
            write(*,'(A,L1)') " - Using kernels generated with TENPI"
         elseif(exa_input%auto_scheme.eq.1) then
            write(*,'(A,L1)') " - Using kernels generated with mbauto"
          elseif(exa_input%auto_scheme.eq.3) then
            write(*,'(A,L1)') " - Using TENPI kernels with TAPP reference implementation"
         endif
      endif
      write(*,'(A,L1)') " - Solve CCD equations :                ", exa_input%ccd
      write(*,'(A,L1)') " - Solve Lambda equations :             ", exa_input%lambda
      write(*,'(A,L1)') " - Compute perturbative corrections :   ", exa_input%do_perturbative
      write(*,'(A,I5)') " - Integral transformation scheme : ", exa_input%moint_scheme
      write(*,'(A,I5)') " - Maximum number of cycles :       ", exa_input%ncycles
      write(*,'(A,ES10.1)') " - Target precision of amplitudes :  ", exa_input%t_econv
      if (exa_input%moint_scheme.eq.42) then
        write(*,'(A,ES10.1)') " - Cholesky decomposition threshold :",exa_input%t_cholesky
      end if
      if (exa_input%level_shift.gt.1.0D-14) then
        write(*,'(A,ES10.1)') " - Level Shift :                      ",exa_input%level_shift
        write(*,'(A)')        " WARNING: MP2 energies are wrong due to level shift"
      end if
      if (exa_input%nff(1).gt.0) then
        write(*,'(A,L1)') " - Adding finite fields                 "
        do i=1,exa_input%nff(1)
          write(*,'(A,A8)') " - - Property :                         ",exa_input%ff_names(i)
          write(*,*) "- - Field strengths :                  ",exa_input%ff(i,1:exa_input%nff(2))
        end do
      end if
      write(*,*) ""
      write(*,'(A,I5)') " - Number of atomic orbitals :      ", nao
      write(*,'(A,I5)') " - Number of occupied spinors :     ", exa_input%nocc
      write(*,'(A,I5)') " - Number of virtual spinors :      ", exa_input%nvir
      write(*,*) ""
      if (exa_input%talsh) then
        write(*,'(A)') " Using the single-node version "
        write(*,'(A,I5,A)') " - TALSH buffer size:               ",exa_input%talsh_buff," GB"
      else
        write(*,'(A)') " Using the multi-node version (MPI) "
        write(*,'(A,I5)') " - ExaTensor block size :           ", exa_input%exa_blocksize
        write(*,'(A)') "!! ATTENTION !! Properties require TALSH, density matrix must be smaller than: "
        write(*,'(A,I5,A)') " - TALSH buffer size:              ",exa_input%talsh_buff," GB"
      end if
      write(*,*) ""
      write(*,'(A,I5)') " - Print Level :                    ", exa_input%print_level
      write(*,*) ""
      write(*,'(A)') " - Estimating memory requirements"
      write(*,'(A,I9,A)') " -- nao  ( scheme 1-3  ) : ",3*16*((nao*nao)/1000)**2/1000," GB"
      write(*,'(A,I9,A)') " -- nvir ( scheme 4/42 ) : ",3*16*((exa_input%nvir*exa_input%nvir)/1000)**2/1000," GB"
      write(*,'(A)') "-----------------------------------------------"
      write(*,*) ""

    end subroutine print_exainput

   subroutine print_mp2fvno_occupation_number_details(threshold, N_reduced, N_OCCUPATION_reduced, Sum_N_OCCUPATION_total)
      integer, intent(in)                            :: N_reduced
      real(kind=8), intent(in), dimension(N_reduced) :: N_OCCUPATION_reduced
      real(kind=8), intent(in)                       :: threshold, Sum_N_OCCUPATION_total

      real(kind=8) :: Sum_N_OCCUPATION_reduced, Perc_OCCUPATION_reduced
      integer :: i

      Sum_N_OCCUPATION_reduced = 0.0
      do i=1, N_reduced
            Sum_N_OCCUPATION_reduced = Sum_N_OCCUPATION_reduced + N_OCCUPATION_reduced(i)
      end do

      Perc_OCCUPATION_reduced = Sum_N_OCCUPATION_reduced / Sum_N_OCCUPATION_total
      Perc_OCCUPATION_reduced = Perc_OCCUPATION_reduced * 100
      write(*,*) ""
      write(*,"(A50,2X,D10.4)") 'MP2 occupation number (o.n.) threshold (t) : ', threshold
      write(*,"(A50,2X,I10)")    'Number of MP2 FVNOs retained (o.n. > t)    : ', N_reduced
      write(*,*) ""
      write(*,"(A50,2X,D10.4)") 'Sum of MP2 FVNOs o.n., retained NOs (a)    : ', Sum_N_OCCUPATION_reduced
      write(*,"(A50,2X,D10.4)") 'Sum of MP2 FVNOs o.n., all NOs (b)         : ', Sum_N_OCCUPATION_total
      write(*,"(A50,2X,D10.4)") 'Percentage of occupation retained (100*a/b): ', Perc_OCCUPATION_reduced
      write(*,*) ""
   endsubroutine print_mp2fvno_occupation_number_details

end module exacorr_printout
