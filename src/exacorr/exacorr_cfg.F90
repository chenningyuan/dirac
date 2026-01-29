module exacorr_cfg

   use exacorr_tensor_interface, only: ETI_DEFAULT_EXECUTOR
  implicit none

!  Collection of everything that can be set via input into one type
!  LV2025: probably better refactor and make a hierarchy of types

   integer, parameter, private:: REALD=8 ! default real for this module

   type, public :: eom_config_t
           integer             :: nroots
           character(len=8)    :: flavor
           logical             :: do_cvs_prj
           logical             :: do_expval
           logical             :: do_linrsp
           logical             :: do_quarsp
           logical             :: do_cubrsp
           integer             :: rhs_expval_states
           integer             :: lhs_expval_states
           character(len=8)    :: properties
           complex(8), pointer :: rsp_pert_freqs(:)
   end type eom_config_t

   type, public :: exacc_input
            logical                 :: talsh = .true.           ! if false, use ExaTensor
            logical                 :: ccd = .false.            ! if false=>CCSD, if true=>CCD
            logical                 :: cc2 = .false.            !L if true CC2

            logical                 :: do_cis = .false.         ! enables CIS, and branch off away from the CC driver 
            logical                 :: do_cis_d = .false.       ! enables CIS(D), and branch off away from the CC driver
            logical                 :: do_cisrt = .false.

            logical                 :: lambda = .false.
            logical                 :: beta_occ = .false.
            logical                 :: beta_vir = .false.
            logical                 :: do_iterative   = .true.  ! this provides a way to switch off iterative methods, e.g. if we want to do only MP2 
            logical                 :: do_perturbative = .true. ! this generalizes the do_triples
            logical                 :: do_triples = .true.
            logical                 :: tripl_lap = .false.
            logical                 :: do_no = .false.          ! if true, generate natural orbital
            logical                 :: do_neo = .false.          ! if true, generate natural effective orbital
            logical                 :: eom_lr = .false.

            logical                 :: do_1rdm_gs=.false.
            logical                 :: do_1rdm_lr=.false.
            logical                 :: do_2rdm_gs=.false.

            logical                 :: auto_kernels = .false.   ! if true, use automatically generated kernels for t/l equations
            integer                 :: auto_scheme  = 2         ! 1: autogen;  2: tenpi; 3: tapp reference implementation

            integer                 :: exa_blocksize = 75
            integer                 :: talsh_buff = 50          ! (GB) total memory available for talsh
            integer                 :: executor = ETI_DEFAULT_EXECUTOR ! executor for tensor contractions (tapp, talsh, exatensor, ...)

            integer                 :: moint_scheme = 4

            integer                 :: print_occupation_number = 0 ! print all occupation number
            integer                 :: mp2fno_n_prop = 1
            integer, allocatable    :: mp2fno_prop(:)
            real(REALD), allocatable :: mp2fno_prop_weights(:)
            logical                 :: mp2fno_weighted_sum = .true.
            integer                 :: mp2no = 0                ! > 0 mp2 natural orbital, < 0 mp2lap
                                                                !  -1: MP2LAP
                                                                !   1: MP2N0
                                                                !   2: MP2NEO
                                                                !   3: MP2FNO++
                                                                !  10: MP2_LAP_NO
            complex(8)              :: mp2fno_prop_omega = (0.D0, 0.D0)       !  Frequency of the perturbation for FNO++

            integer                 :: ncycles = 30
            integer                 :: print_level = 0
            integer                 :: nkr_occ, nkr_vir         ! number of Kramer pairs
            integer                 :: nocc, nvir               ! number of spinors
            integer                 :: iterative_level = 2      ! at which excitation level iterative amplitudes/lagrange mult. will be calculated: 1, 2, 3, 4...
            integer                 :: iterative_scheme = 21    ! scheme/variant for the iterative determination of amplitudes 
                                                                ! 10: CCS (not currently implemented, but logically should follow as this)
                                                                ! 20: CCD      ; 21: CCSD       ; 22: CC2        ; 
                                                                ! 30: CCSDT    ; 31: CCSDT-1a   ; 32: CCSDT-1b   ; 33: CCSDT-2   ; 34: CCSDT-3   ; 35:  CC3
                                                                ! 40: CCSDTQ   ; 41: CCSDTQ-1a  ; 42: CCSDTQ-1b  ; 43: CCSDTQ-2  ; 44: CCSDTQ-3  ; 45:  CC4
                                                                ! 50: CCSDTQP  ; 51: CCSDTQP-1a ; 52: CCSDTQP-1b ; 53: CCSDTQP-2 ; 54: CCSDTQP-3 ; 55:  CC5
            integer                 :: perturbative_level = 3   ! at which excitation level pertubative corrections will be calculated: 2, 3, 4...
            integer                 :: perturbative_scheme = 30 ! scheme/variant for the iterative determination of energies/amplitudes/... 
                                                                ! 20: MP2 (MP2 is always performed)
                                                                ! 30: CCSD(T)

!  aspg, TODO: harmonize this tripl_block with the more general setup of excitation rank
            integer                 :: tripl_block(1:2) = -1    ! for running triples
            integer                 :: nff(1:2) = 0             ! num. prop, num. field strengths
            integer                 :: n_laplace = 10           ! number of laplace points to use

! aspg, response-related
            integer                 :: N_eom_trial = 1          ! number of trial vectors in CC linear response
            integer                 :: wf_type = 2              ! wave-function type
            integer                 :: N_Prop_A, N_Prop_B    ! number of property
            integer                 :: N_Prop_C, N_Prop_D    ! number of property
            integer                 :: N_omega_lr = 1               ! number of frequency in cclr
            integer                 :: N_omega_qr1 = 1               ! number of frequency1 in ccqr
            integer                 :: N_omega_qr2 = 1               ! number of frequency2 in ccqr
            integer                 :: N_Time_anti = 0          ! Number of time-antisymmetric operators
            logical                 :: ccrsp_symmetric = .true. ! Use symmetric formulation in CC response

            integer, allocatable    :: mokr_occ(:), mokr_vir(:) ! list of active Kramer pairs
            integer, allocatable    :: mo_occ(:), mo_vir(:)     ! list of active spinors
            integer, allocatable    :: Prop_A(:), Prop_B(:), Prop_C(:), Prop_D(:) ! index of Property

            INTEGER                 ::  TPA_state             ! index of final excited state in eom-tpa
            real(REALD), allocatable :: w_lr_re(:)            ! real frequency in CC linear response
            real(REALD), allocatable :: w_lr_im(:)            ! imaginary frequency in CC response

            real(REALD), allocatable :: w_qr1_re(:)            ! 1st real frequency in CC quadratic response
            real(REALD), allocatable :: w_qr1_im(:)            ! 1st imaginary frequency in CC quadratic response
            real(REALD), allocatable :: w_qr2_re(:)            ! 2nd real frequency in CC quadratic response
            real(REALD), allocatable :: w_qr2_im(:)            ! 2nd imaginary frequency in CC quadratic response

! NEO-related inputs
            integer, allocatable    :: neo_cis_states(:)       ! the CIS states  used to construct NEO orbitals
                                                               ! in the form of an array to allow granularity on which states to use
            integer                 :: neo_cis_nstates = 1     ! the number of CIS states used to construce the NEO orbitals

! real-time propagation input parameters (currently only cis-based)
! cis(d)-related input
            integer                 :: cis_d_pmin = 1
            integer                 :: cis_d_pmax = 5
! propagation and field 
            integer                 :: rt_Nsteps = 300
            integer                 :: rt_fieldtype = 1
            real(REALD)             :: rt_dt = 1.0D-2         

! other thresholds
            real(REALD)             :: t_econv = 1.0D-9         ! convergence criterium
            real(REALD)             :: level_shift = 0.D0       ! level shift
            real(REALD)             :: t_cholesky = 1.0D-9      ! threshold for cholesky / mp2lap
            real(REALD)             :: t_mp2no = 1.0D-5         ! threshold for mp2 natural orbital occupation number


            real(REALD), allocatable:: rew_energies(:,:) ! restricted excitation window, in terms of energy thresholds (a.u.)

! this just signals we should handle *EOMCC inputs under **EXACC
            logical                 :: eom = .false.
            integer                 :: current_eom_input = -1
        !     logical                 :: eom_lr = .false.
            logical                 :: eom_qr = .false.
            logical                 :: qr_tpa = .false.

! handle one or more *EOMCC blocks 
            type(eom_config_t), allocatable :: eom_input(:)

            complex(8), allocatable   :: ff(:,:)                  ! field strengths
            character(8), allocatable :: ff_names(:)            ! property labels

  end type exacc_input

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
  public exacc_sync_cw
#endif
  public exacc_deallocate_cw

  integer, parameter :: LENGTH=100
  character(LENGTH), public :: string_occupied, string_virtual, string_occ_beta, string_vir_beta
  type(exacc_input), public :: exa_input

contains

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))

  subroutine exacc_sync_cw(exa_MPI_rank,my_MPI_rank)
!
! please add bcast here if a new variable has been added to the list above
!
  use interface_to_mpi

  integer :: exa_MPI_rank,my_MPI_rank,i,j

  call interface_mpi_bcast_l0(exa_input%ccd,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%cc2,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%eom,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%lambda,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%beta_occ,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%beta_vir,1,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast_l0(exa_input%do_cis,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%do_cis_d,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%do_cisrt,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%do_neo,1,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast_l0(exa_input%do_perturbative,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%do_iterative,1,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast_l0(exa_input%tripl_lap,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%do_triples,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%do_1rdm_gs,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%do_1rdm_lr,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%do_2rdm_gs,1,exa_MPI_RANK,global_communicator)


  call interface_mpi_bcast_l0(exa_input%auto_kernels,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%ccrsp_symmetric,1,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast(exa_input%exa_blocksize,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%talsh_buff,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%moint_scheme,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%mp2no,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%ncycles,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%print_level,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%tripl_block,2,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%n_laplace,1,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast(exa_input%neo_cis_nstates,1,exa_MPI_RANK,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%neo_cis_states(exa_input%neo_cis_nstates))
  call interface_mpi_bcast(exa_input%neo_cis_states,1,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast(exa_input%cis_d_pmin,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%cis_d_pmax,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%rt_Nsteps,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%rt_fieldtype,1,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast(exa_input%iterative_level,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%iterative_scheme,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%perturbative_level,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%perturbative_scheme,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%auto_scheme,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%nkr_occ,1,exa_MPI_rank,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%mokr_occ(exa_input%nkr_occ))
  call interface_mpi_bcast(exa_input%nkr_vir,1,exa_MPI_rank,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%mokr_vir(exa_input%nkr_vir))
  call interface_mpi_bcast(exa_input%nocc,1,exa_MPI_rank,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%mo_occ(exa_input%nocc))
  call interface_mpi_bcast(exa_input%nvir,1,exa_MPI_rank,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%mo_vir(exa_input%nvir))

  call interface_mpi_bcast(exa_input%mokr_occ,exa_input%nkr_occ,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%mokr_vir,exa_input%nkr_vir,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%mo_occ,exa_input%nocc,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%mo_vir,exa_input%nvir,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast(exa_input%t_econv,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%t_mp2no,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%t_rt_dt,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%level_shift,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%t_cholesky,1,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast(exa_input%t_rt_dt,1,exa_MPI_RANK,global_communicator)

! aspg, the information on response could perhaps be regrouped into one datatype
  call interface_mpi_bcast(exa_input%N_eom_trial,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%wf_type,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%N_Prop_A,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%N_Prop_B,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%N_Prop_C,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%N_Prop_D,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%N_omega_lr,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%N_omega_qr1,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%N_omega_qr2,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%N_Time_anti,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%TPA_state,1,exa_MPI_RANK,global_communicator)

  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%Prop_A(exa_input%N_Prop_A))
  call interface_mpi_bcast(exa_input%Prop_A,1,exa_MPI_RANK,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%Prop_B(exa_input%N_Prop_B))
  call interface_mpi_bcast(exa_input%Prop_B,1,exa_MPI_RANK,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%Prop_C(exa_input%N_Prop_C))
  call interface_mpi_bcast(exa_input%Prop_C,1,exa_MPI_RANK,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%Prop_D(exa_input%N_Prop_D))
  call interface_mpi_bcast(exa_input%Prop_D,1,exa_MPI_RANK,global_communicator)

  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%w_lr_re(exa_input%N_omega_lr))
  call interface_mpi_bcast(exa_input%w_lr_re,1,exa_MPI_RANK,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%w_lr_im(exa_input%N_omega_lr))
  call interface_mpi_bcast(exa_input%w_lr_im,1,exa_MPI_RANK,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%w_qr1_re(exa_input%N_omega_qr1))
  call interface_mpi_bcast(exa_input%w_qr1_re,1,exa_MPI_RANK,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%w_qr1_im(exa_input%N_omega_qr1))
  call interface_mpi_bcast(exa_input%w_qr1_im,1,exa_MPI_RANK,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%w_qr2_re(exa_input%N_omega_qr2))
  call interface_mpi_bcast(exa_input%w_qr2_re,1,exa_MPI_RANK,global_communicator)
  if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%w_qr2_im(exa_input%N_omega_qr2))
  call interface_mpi_bcast(exa_input%w_qr2_im,1,exa_MPI_RANK,global_communicator)

  if (allocated(exa_input%eom_input)) then
     do i = 1, size(exa_input%eom_input,1)
        call interface_mpi_bcast(exa_input%eom_input(i)%nroots,1,exa_MPI_RANK,global_communicator)
        call interface_mpi_bcast(exa_input%eom_input(i)%flavor,1,exa_MPI_RANK,global_communicator)
        call interface_mpi_bcast_l0(exa_input%eom_input(i)%do_cvs_prj,1,exa_MPI_RANK,global_communicator)
        call interface_mpi_bcast_l0(exa_input%eom_input(i)%do_expval,1,exa_MPI_RANK,global_communicator)
        call interface_mpi_bcast_l0(exa_input%eom_input(i)%do_linrsp,1,exa_MPI_RANK,global_communicator)
        call interface_mpi_bcast_l0(exa_input%eom_input(i)%do_quarsp,1,exa_MPI_RANK,global_communicator)
        call interface_mpi_bcast_l0(exa_input%eom_input(i)%do_quarsp,1,exa_MPI_RANK,global_communicator)
        call interface_mpi_bcast(exa_input%eom_input(i)%rhs_expval_states,1,exa_MPI_RANK,global_communicator)
        call interface_mpi_bcast(exa_input%eom_input(i)%lhs_expval_states,1,exa_MPI_RANK,global_communicator)
        call interface_mpi_bcast(exa_input%eom_input(i)%properties,1,exa_MPI_RANK,global_communicator)
!       do j = 1, size(exa_input%eom_input(i)%rsp_pert_freqs)
! this requires extension of bcast to complex types, if the frequencies are complex numbers
!          call interface_mpi_bcast(exa_input%eom_input(i)%rsp_pert_freqs(j),1,exa_MPI_RANK,global_communicator)
!       end do
     end do
  end if

  call interface_mpi_bcast(exa_input%nff,2,exa_MPI_RANK,global_communicator)
  if (exa_input%nff(1)>0) then
    if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%ff_names(exa_input%nff(1)))
    call interface_mpi_bcast(exa_input%ff_names,len(exa_input%ff_names(1))*exa_input%nff(1), &
                           exa_MPI_RANK,global_communicator)
    if (my_MPI_rank /= exa_MPI_rank) allocate(exa_input%ff(exa_input%nff(1),exa_input%nff(2)))
    call interface_mpi_bcast(exa_input%ff,exa_input%nff(1)*exa_input%nff(2)*2,exa_MPI_RANK,global_communicator)
  end if

  end subroutine exacc_sync_cw
#endif


  subroutine exacc_deallocate_cw()
!
! please add bcast here if a new variable has been added to the list above
!

  deallocate (exa_input%mokr_occ)
  deallocate (exa_input%mokr_vir)
  deallocate (exa_input%mo_occ  )
  deallocate (exa_input%mo_vir  )

  if (exa_input%nff(1)>0) then
    deallocate (exa_input%ff_names)
    deallocate (exa_input%ff)
  end if

  end subroutine exacc_deallocate_cw

end module
