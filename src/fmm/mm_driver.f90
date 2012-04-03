! mm_driver:
! General module for driving multipole method calculations. 
! 2002 (c) Mark Watson(*), Pawel Salek
! 
! Notes: traditional dalton's WORK memory is used for some of the
! larger allocations before memory allocation in dalton is "fixed".
!
#if defined(SYS_AIX) || defined(VAR_PGF77)
   ! for some strange reason, this fails on AIX.
   #define ALLOCATION_WITHOUT_WORK 1
#endif

MODULE mm_main_driver

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_get_E_via_raw_potentials,      &
             mm_get_J_via_raw_potentials,      &
             mm_get_int_type

   REAL(REALK), POINTER, SAVE :: Vff(:,:)
   TYPE(raw_mm_data),    SAVE :: LHS_mms, RHS_mms

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_driver(scheme,WORK,LWORK)

      USE density_fitting
      USE mm_interface_mod,   ONLY: mm_get_raw_data
      USE mm_aux_qlm_builder, ONLY: mm_init_qlm_builder,   &
                                    mm_get_aux_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK),      INTENT(IN) :: LWORK
      REAL(REALK),     INTENT(INOUT) :: WORK(:) !"input value of LWORK" words

      CALL mm_nullify_local_pointers
      ! get multipole data straight from interface file
      CALL mm_get_raw_data(scheme,LHS_mms,RHS_mms)

      ! now ensure appropriate prefactors and normalisation
      CALL mm_init_qlm_builder
      CALL mm_get_aux_qlm(scheme,LHS_mms,RHS_mms)

      ! allocate the far field potential
      CALL mm_allocate_Vff(scheme,WORK,LWORK)

   END SUBROUTINE mm_init_driver

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_driver

      USE mm_aux_qlm_builder, ONLY: mm_free_qlm_builder

      IMPLICIT NONE
      CALL mm_free_qlm_builder
      IF (.NOT.ASSOCIATED(Vff)) CALL QUIT('Vff should be allocated by now!')
#ifdef ALLOCATION_WITHOUT_WORK
      DEALLOCATE(Vff)
#endif
      CALL mm_nullify_local_pointers

   END SUBROUTINE mm_free_driver

!-------------------------------------------------------------------------------

   SUBROUTINE mm_nullify_local_pointers

      IMPLICIT NONE
      NULLIFY (Vff)
      NULLIFY (LHS_mms%paras, LHS_mms%dens, LHS_mms%batch_map)
      NULLIFY (LHS_mms%qlm, LHS_mms%qlm_W, LHS_mms%qlm_T)
      NULLIFY (RHS_mms%paras, RHS_mms%dens, RHS_mms%batch_map)
      NULLIFY (RHS_mms%qlm, RHS_mms%qlm_W, RHS_mms%qlm_T)

   END SUBROUTINE mm_nullify_local_pointers

!-------------------------------------------------------------------------------

   SUBROUTINE mm_allocate_Vff(scheme,WORK,LWORK)

      USE mm_memory_manager_mod, ONLY: mm_allocate

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK),   INTENT(IN)    :: LWORK
      REAL(REALK),     INTENT(INOUT) :: WORK(:) ! input value of LWORK 

      INTEGER(INTK) :: i, lm_dim, mms_dim, max_Lmin

      IF (.NOT.ASSOCIATED(LHS_mms%paras)) STOP 'mms ptrs not set in mm_driver!'
      IF (ASSOCIATED(Vff)) CALL QUIT('Vff should NOT be allocated already!')

      mms_dim = SIZE(LHS_mms%paras)
      IF (scheme%dynamic_LMAX_on) THEN
         max_Lmin = 0
         DO i = 1, mms_dim
            max_Lmin = MAX(max_Lmin, LHS_mms%paras(i)%Lmin)
         END DO
         max_Lmin = MIN((max_Lmin+scheme%LEXTRA), scheme%raw_LMAX)
         lm_dim = (1+ max_Lmin)**2
      ELSE
         ! We assume we are generating a RAW potential with limited l,m
         ! ie. we shouldn't use this array for translated, BOXED potentials
         lm_dim = (1+ scheme%raw_LMAX)**2
      END IF

#ifdef ALLOCATION_WITHOUT_WORK
      CALL mm_allocate(MEM_RAW_VFF,Vff,lm_dim,mms_dim)
#else
      CALL set_v(Vff,lm_dim, mms_dim, WORK, LWORK/lm_dim)
      !LWORK = LWORK - lm_dim*mms_dim
#endif
      ! must zero as Vff is build additively in general
      Vff(:,:) = zero

   CONTAINS

      SUBROUTINE set_v(v,NROWS,NCOLS,W,NWORK)

         IMPLICIT NONE
         INTEGER(INTK), INTENT(IN) :: NROWS, NCOLS, NWORK
         REAL(REALK), POINTER      :: v(:,:)
         REAL(REALK), TARGET       :: W(NROWS,NWORK)

         IF(NWORK.LT.NCOLS) THEN
            WRITE(LUPRI,*) 'WORK space exceded for Vff!'
            WRITE(LUPRI,*) 'Require', NROWS*NCOLS, 'words'
            WRITE(LUPRI,*) 'Avaliable words =', NWORK*NROWS
            CALL QUIT('Work space exceded!')
         END IF
         v => w(:,1:NCOLS)

      END SUBROUTINE set_v

   END SUBROUTINE mm_allocate_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_E_via_raw_potentials(scheme,energy,text,WORK,LWORK)

      USE mm_Vff_driver,    ONLY: mm_get_raw_Vff
      USE mm_Vff_processor, ONLY: mm_get_E_from_Vff, mm_get_E_from_pkd_Vff
      USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(INOUT) :: scheme
      REAL(REALK),        INTENT(OUT)   :: energy
      CHARACTER(*),       INTENT(OUT)   :: text
      INTEGER(INTK),      INTENT(IN)    :: LWORK
      REAL(REALK),        INTENT(INOUT) :: WORK(:) ! "input value of LWORK"

      ! we always have the density on both sides when getting energies
      scheme%LHS_dens = .TRUE.
      scheme%RHS_dens = .TRUE.

      ! set up T-pair generator algorithm
      IF (scheme%LHS_mm_range == scheme%RHS_mm_range) THEN
         IF (.NOT.scheme%T_searcher(DO_FQ)%all_square)  & 
            scheme%T_searcher(DO_FQ)%shape = SHAPE_TRIANGULAR
         IF (.NOT.scheme%T_searcher(DO_NN)%all_square)  & 
            scheme%T_searcher(DO_NN)%shape = SHAPE_TRIANGULAR
      END IF

      ! now build energy additively
      energy = zero
      ! get moments
      CALL mm_init_driver(scheme,WORK,LWORK)
      ! get potential
      CALL mm_get_raw_Vff(scheme,LHS_mms%paras,RHS_mms,Vff)
      ! get energy
      IF (scheme%pack_LHS) THEN
         CALL mm_get_E_from_pkd_Vff(scheme,LHS_mms,Vff,energy,text)
      ELSE
         CALL mm_get_E_from_Vff(scheme,LHS_mms,Vff,energy,text)
      END IF
      CALL mm_free_driver
     
   END SUBROUTINE mm_get_E_via_raw_potentials

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_J_via_raw_potentials(scheme,J_matrix,energy,txt,WORK,LWRK)

      USE mm_Vff_driver,    ONLY: mm_get_raw_Vff
      USE mm_Vff_processor, ONLY: mm_get_J_from_Vff,       &
                                  mm_get_J_from_pkd_Vff,   &
                                  mm_get_E_from_pkd_Vff,   &
                                  mm_get_E_from_Vff
      USE mm_qlm_processor, ONLY: mm_factor_in_dens
      USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(INOUT) :: scheme
      REAL(REALK),        INTENT(OUT)   :: J_matrix(:,:)
      REAL(REALK),        INTENT(OUT)   :: energy
      CHARACTER(*),       INTENT(OUT)   :: txt
      INTEGER(INTK),      INTENT(IN)    :: LWRK
      REAL(REALK),        INTENT(INOUT) :: WORK(:)

      ! we only have the density on RHS when getting J-matrix
      ! via far-field potential
      scheme%LHS_dens = .FALSE.
      scheme%RHS_dens = .TRUE.
      ! get moments
      CALL mm_init_driver(scheme,WORK,LWRK)

      ! get potential
      CALL mm_get_raw_Vff(scheme,LHS_mms%paras,RHS_mms,Vff)

      ! get J-matrix
      J_matrix = zero
      energy = zero
      IF (scheme%pack_LHS) THEN
         CALL mm_get_J_from_pkd_Vff(scheme,LHS_mms,Vff,J_matrix)
         ! get energy after factoring in density to LHS
         CALL mm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
         CALL mm_get_E_from_pkd_Vff(scheme,LHS_mms,Vff,energy,txt)
      ELSE
         CALL mm_get_J_from_Vff(scheme,LHS_mms,Vff,J_matrix)
         ! get energy after factoring in density to LHS
         CALL mm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
         CALL mm_get_E_from_Vff(scheme,LHS_mms,Vff,energy,txt)
      END IF

      CALL mm_free_driver

   END SUBROUTINE mm_get_J_via_raw_potentials

!-------------------------------------------------------------------------------
! When not combining moments, it is possible to save a factor of two in the
! contractions by avoiding double-counting.  This is only possible when we
! do not build the potential (which mixes RHS contributions).
! i.e. used in NN phase to obtain factor of 2 speed-up.
!
!   SUBROUTINE mm_get_J_directly(scheme,J_matrix,WORK,LWORK)
!
!      USE mm_Vff_driver,      ONLY: mm_get_raw_Vff
!      USE mm_Vff_processor,   ONLY: mm_get_J_from_Vff
!      USE mm_J_driver,        ONLY: mm_drive_NN_J_builder
!      USE mm_qlm_processor,   ONLY: mm_factor_in_dens
!
!      IMPLICIT NONE
!      TYPE(scheme_paras), INTENT(INOUT) :: scheme
!      REAL(REALK),        INTENT(OUT)   :: J_matrix(:,:)
!      INTEGER(INTK),      INTENT(IN)    :: LWORK
!      REAL(REALK),        INTENT(INOUT) :: WORK(:) ! input value of LWORK
!
!      REAL(REALK), POINTER :: dens_ptr(:)
!
!      ! set up T-pair generator algorithm
!      IF (scheme%LHS_mm_range == scheme%RHS_mm_range) THEN
!         IF (.NOT.scheme%T_searcher(DO_FQ)%all_square)  & 
!            scheme%T_searcher(DO_FQ)%shape = SHAPE_TRIANGULAR
!         IF (.NOT.scheme%T_searcher(DO_NN)%all_square)  & 
!            scheme%T_searcher(DO_NN)%shape = SHAPE_TRIANGULAR
!      END IF
!
!      J_matrix = zero
!
!      ! get raw moments for nearest neighbour phase without densities
!      scheme%LHS_dens = .FALSE.
!      scheme%RHS_dens = .FALSE.
!      CALL mm_init_driver(scheme,WORK,LWORK)
!
!      ! get NN part (or total if FQ algorithm) without going via Vff
!      IF (scheme%LHS_mm_range /= scheme%RHS_mm_range) THEN
!         CALL QUIT('ranges must be equal to get_J_directly')
!      END IF
!      ! set up density pointer for J-matrix build
!      dens_ptr => RHS_mms%dens(:)  ! LHS .eq. RHS
!      CALL mm_drive_NN_J_builder(scheme,LHS_mms,RHS_mms,dens_ptr,J_matrix)
!
!      ! reset moments for far-field phase
!      ! we add the density to RHS when getting J-matrix via far-field potential
!!FIMXE: add checks because this wont work if we compact RHS first
!      CALL mm_factor_in_dens(RHS_mms%dens,RHS_mms%qlm_T)
!      CALL mm_factor_in_dens(RHS_mms%dens,RHS_mms%qlm_W)
!
!      ! get potential for FF work
!      IF (scheme%algorithm /= DO_FQ) THEN
!         CALL mm_get_raw_Vff(scheme,LHS_mms%paras,RHS_mms,Vff,skip_NN=.TRUE.)
!         CALL mm_get_J_from_Vff(scheme,LHS_mms,Vff,J_matrix)
!      END IF
!
!      CALL mm_free_driver
!
!   END SUBROUTINE mm_get_J_directly
!
!-------------------------------------------------------------------------------

   INTEGER FUNCTION mm_get_int_type(n_el)

      IMPLICIT NONE
      CHARACTER(6), INTENT(IN) :: n_el

      SELECT CASE(n_el)
      CASE('ONE_EL')
         mm_get_int_type = INTTYPE_ONE_EL
      CASE('TWO_EL')
         mm_get_int_type = INTTYPE_TWO_EL
      CASE('FULL_J')
         mm_get_int_type = INTTYPE_FULL_J
      CASE('NUC_AT')
         mm_get_int_type = INTTYPE_NUC_AT
      CASE DEFAULT
         mm_get_int_type = -1
      END SELECT

   END FUNCTION mm_get_int_type

!-------------------------------------------------------------------------------

END MODULE mm_main_driver

!===============================================================================

SUBROUTINE mm_get_energy(ERI_energy,n_el,WORK,LWRK)

   USE mm_global_paras_mod
   USE mm_stats_mod
   USE mm_memory_manager_mod, ONLY: mm_init_mem_man, mm_close_mem_man
   USE mm_scheme_builder, ONLY: mm_get_scheme
   USE mm_interface_mod,  ONLY: mm_connect_interface, mm_disconnect_interface
   USE mm_main_driver,    ONLY: mm_get_E_via_raw_potentials, mm_get_int_type

   IMPLICIT NONE
   REAL(REALK),   INTENT(IN)    :: ERI_energy
   CHARACTER(6),  INTENT(IN)    :: n_el
   INTEGER(INTK), INTENT(IN)    :: LWRK
   REAL(REALK),   INTENT(INOUT) :: WORK(LWRK)

   REAL(REALK)   :: energy
   CHARACTER(36) :: E_text
   LOGICAL       :: A, B
   INTEGER(INTK) :: INUM, KWRK, integral_type
   TYPE(scheme_paras), POINTER :: scheme

   KWRK = 1

   CALL mm_init_mem_man
   CALL mm_get_scheme(scheme)
   CALL mm_init_stats(scheme,GET_ENERGY)

   CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)

   integral_type = mm_get_int_type(n_el)
   SELECT CASE (integral_type)
   CASE (INTTYPE_ONE_EL)
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = NUCLEAR_ONLY 
   CASE (INTTYPE_TWO_EL)
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = ELECTRONIC_ONLY 
   CASE (INTTYPE_NUC_AT)
      scheme%LHS_mm_range = NUCLEAR_ONLY
      scheme%RHS_mm_range = NUCLEAR_ONLY 
   CASE DEFAULT
      CALL QUIT ('can only do n-n e-e or e-n energy calculation!')
   END SELECT

   CALL mm_get_E_via_raw_potentials(scheme,energy,E_text,              &
                                    WORK(KWRK:LWRK),LWRK-KWRK+1)

   !PRINT '(X,A," = ",F22.12)', E_text, energy
   write(LUPRI,'(X,A," = ",F22.12)') E_text, energy

   A = (scheme%LHS_mm_range == ELECTRONIC_ONLY)
   B = (scheme%RHS_mm_range == ELECTRONIC_ONLY)
   IF (A.AND.B) THEN
     !PRINT '(" error wrt ERI e-e energy   =",E11.3)', ABS(energy-ERI_energy)
     write(LUPRI,'(" error wrt ERI e-e energy =",E11.3)')ABS(energy-ERI_energy)
   END IF

   CALL mm_print_stats
   CALL mm_disconnect_interface
   CALL mm_close_mem_man

END SUBROUTINE mm_get_energy

!-------------------------------------------------------------------------------

SUBROUTINE mm_get_J_matrix(n_el,Fock_matrix,Nfock1,Nfock2,WORK,LWRK)

   USE mm_global_paras_mod
   USE mm_stats_mod
   USE mm_memory_manager_mod, ONLY: mm_init_mem_man, mm_close_mem_man
   USE mm_scheme_builder, ONLY: mm_get_scheme
   USE mm_interface_mod,  ONLY: mm_connect_interface, mm_disconnect_interface
   USE mm_main_driver,    ONLY: mm_get_J_via_raw_potentials, mm_get_int_type

   IMPLICIT NONE
   CHARACTER(6),  INTENT(IN)    :: n_el
   REAL(REALK),   INTENT(INOUT) :: Fock_matrix(Nfock1,Nfock2) 
   INTEGER(INTK), INTENT(IN)    :: LWRK, Nfock1,Nfock2
   REAL(REALK),   INTENT(INOUT) :: WORK(LWRK)

!FIXME: dont want to store a separate J_matrix if we can pass Fock_matrix
   REAL(REALK), ALLOCATABLE :: J_matrix(:,:)
   REAL(REALK)        :: energy
   INTEGER(INTK)      :: rank, ierr, integral_type, KWRK
   CHARACTER(36)      :: E_text
   LOGICAL            :: branch_free_flag
   TYPE(scheme_paras), POINTER :: scheme

#if defined(VAR_MPI)
   INCLUDE 'mpif.h'
#endif
   DOUBLE PRECISION SECOND, T1

   CALL mm_init_stats(scheme,GET_JMAT)
!   stat_NBAST = NBAST

   T1 = SECOND()
   CALL mm_init_mem_man
   integral_type = mm_get_int_type(n_el)
   stat_1el_2el_or_full = mm_get_int_type(n_el)

#if defined(VAR_MPI)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,IERR)
   CALL mm_getj_wake_up_slaves()
   IF(rank==0) CALL mm_get_scheme(scheme)
   CALL mm_getj_sync_slaves
#else
   CALL mm_get_scheme(scheme)
#endif
   KWRK = 1
   CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)

   branch_free_flag = scheme%branch_free
   ALLOCATE(J_matrix(Nfock1,Nfock2))
   SELECT CASE (integral_type)
   CASE (INTTYPE_ONE_EL)
      ! always run 1-el part without "branch-free" algorithm
      scheme%branch_free = .FALSE.
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = NUCLEAR_ONLY
      CALL mm_get_J_via_raw_potentials(scheme,J_matrix,energy,E_text,  &
                                       WORK(KWRK:LWRK),LWRK-KWRK+1)
      IF ( MM_STANDALONE ) THEN
         WRITE(LUPRI,'(X,A," = ",E20.12)') E_text, energy
      END IF
      !WRITE(LUPRI,'(A,E20.12)') ' MAX(ABS(J_matrix)) =',MAXVAL(ABS(J_matrix))
   CASE (INTTYPE_TWO_EL)
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = ELECTRONIC_ONLY
      CALL mm_get_J_via_raw_potentials(scheme,J_matrix,energy,E_text,  &
                                       WORK(KWRK:LWRK),LWRK-KWRK+1)
      !WRITE(LUPRI,'(X,A," = ",E20.12)') E_text, energy
      !WRITE(LUPRI,'(A,E20.12)') ' MAX(ABS(J_matrix)) =',MAXVAL(ABS(J_matrix))
   CASE (INTTYPE_FULL_J)
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = ALL_MOMENTS
      CALL mm_get_J_via_raw_potentials(scheme,J_matrix,energy,E_text,  &
                                       WORK(KWRK:LWRK),LWRK-KWRK+1)
      WRITE(LUPRI,'(X,A," = ",E20.12)') E_text, energy
   CASE DEFAULT
      CALL QUIT ('require 1, 2, or full J_matrix build!')
   END SELECT
   scheme%branch_free = branch_free_flag

#ifdef VAR_MPI
   CALL mm_getj_collect(SECOND()-T1)
#endif

   Fock_matrix = Fock_matrix + J_matrix
   DEALLOCATE(J_matrix)

   !write(6,*) "J_matrix from MM:"
   !CALL OUTPUT(J_matrix,1,Nfock1,1,Nfock2,Nfock1,Nfock2,1,6)

   CALL mm_print_stats
   CALL mm_disconnect_interface
   CALL mm_close_mem_man

END SUBROUTINE mm_get_J_matrix

!-------------------------------------------------------------------------------

#ifdef VAR_MPI
SUBROUTINE mm_getj_sync_slaves
   ! can be split into two parts: syncing scheme, and other parts.
   IMPLICIT NONE
   INCLUDE 'mpif.h'
   INTEGER ierr
   INTEGER, PARAMETER :: master = 0

#ifdef FIX_THIS_CODE
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNUM, IERR)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NNODES,IERR)
   CALL MPI_Bcast(integral_type,1,MPI_INTEGER,master,MPI_COMM_WORLD,IERR)
   CALL bcast_scheme_paras(master)
   CALL bcast_T_searcher_type(master,scheme%T_searcher)
   CALL bcast_T_contract_schm(master,scheme%T_con)
   CALL bcast_W_contract_schm(master,scheme%W_con)
#else
   call quit('fix code in mm_getj_sync_slaves')
#endif
END SUBROUTINE mm_getj_sync_slaves

! note that broadcasting system size is redundant, it is not an
! input PARAMETER.
SUBROUTINE bcast_scheme_paras(master)
   IMPLICIT NONE
   INCLUDE 'mpif.h'
   INTEGER, INTENT(IN) :: master
   INTEGER intEx, logEx, doubleEx, ierr, theType
   INTEGER, DIMENSION(4) :: bllen, types, disps
#ifdef FIX_THIS_CODE
   CALL MPI_TYPE_EXTENT(MPI_INTEGER,intex,ierr)
   CALL MPI_TYPE_EXTENT(MPI_LOGICAL,logex,ierr)
   CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,doubleEX,ierr)
   
   types = (/ MPI_INTEGER, MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION/)
   bllen = (/ 2,           4,          6,          4                   /)
   disps(1) = 0; 
   disps(2) = disps(1) + bllen(1)*intEx
   disps(3) = disps(2) + bllen(2)*logEx
   disps(4) = disps(3) + bllen(3)*intEx
   CALL MPI_Type_Struct(SIZE(bllen),bllen,disps,types,theType, ierr)
   CALL MPI_TYPE_COMMIT(theType,ierr)
   CALL MPI_Bcast(scheme,1,theType,master,MPI_COMM_WORLD,IERR)
   CALL MPI_TYPE_FREE(theType,ierr)
#else
   call quit('fix code in bcast_scheme_paras')
#endif
END SUBROUTINE bcast_scheme_paras

SUBROUTINE bcast_T_searcher_type(master,searcher)
   IMPLICIT NONE
   INCLUDE 'mpif.h'
   INTEGER, INTENT(IN)   :: master
#ifdef FIX_THIS_CODE
   TYPE(T_searcher_type) :: searcher(5)
#else
   INTEGER :: searcher
#endif

   INTEGER intEx, logEx, ierr, theType
   INTEGER, DIMENSION(2) :: bllen, types, disps

#ifdef FIX_THIS_CODE
   CALL MPI_TYPE_EXTENT(MPI_INTEGER,intex,ierr)
   CALL MPI_TYPE_EXTENT(MPI_LOGICAL,logex,ierr)
   
   types = (/ MPI_INTEGER, MPI_LOGICAL /)
   bllen = (/ 2,           1           /)
   disps(1) = 0; 
   disps(2) = disps(1) + bllen(1)*intEx
   CALL MPI_Type_Struct(SIZE(bllen),bllen,disps,types,theType, IERR)
   CALL MPI_TYPE_COMMIT(theType,IERR)
   CALL MPI_Bcast(searcher,5,theType,master,MPI_COMM_WORLD,IERR)
   CALL MPI_TYPE_FREE(theType,ierr)
#else
   call quit('fix code in bcast_T_searcher_type')
#endif
END SUBROUTINE bcast_T_searcher_type

SUBROUTINE bcast_T_contract_schm(master,contractor)
   IMPLICIT NONE
   INTEGER, INTENT(IN)                  :: master
#ifdef FIX_THIS_CODE
   TYPE(T_contract_schm), INTENT(INOUT) :: contractor
#else
   integer ::                              contractor
#endif
   INCLUDE 'mpif.h'
   INTEGER ierr
#ifdef FIX_THIS_CODE
   CALL MPI_Bcast(contractor,6,MPI_INTEGER,master,MPI_COMM_WORLD,IERR)
#else
   call quit('fix code in bcast_T_contract_schm')
#endif
END SUBROUTINE bcast_T_contract_schm

SUBROUTINE bcast_W_contract_schm(master,w_con)
   IMPLICIT NONE
   INCLUDE 'mpif.h'
   INTEGER, INTENT(In)                  :: master
#ifdef FIX_THIS_CODE
   TYPE(W_contract_schm), INTENT(INOUT) :: w_con
   integer ierr
   CALL MPI_Bcast(w_con,3,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
#else
   INTEGER::                               w_con
   call quit('fix code in bcast_W_contract_schm')
#endif
END SUBROUTINE bcast_W_contract_schm

SUBROUTINE mm_getj_collect(tm)
   IMPLICIT NONE
   INCLUDE 'mpif.h'
   INCLUDE 'dcbbas.h'
   INTEGER IERR, MASTER,I
   REAL(8) :: Work(ntbas(0)*ntbas(0))
   DOUBLE PRECISION, INTENT(IN) :: tm
#ifdef FIX_THIS_CODE
   DOUBLE PRECISION  :: tms(NNODES)
   MASTER = 0
   CALL DCOPY(NBAST*NBAST,J_matrix, 1, Work, 1)
   CALL MPI_Reduce(work,J_matrix,NBAST*NBAST,MPI_DOUBLE_PRECISION,&
        & MPI_SUM, MASTER, MPI_COMM_WORLD, IERR)
   CALL MPI_GATHER(tm,1,MPI_DOUBLE_PRECISION,&
        &          tms,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,IERR)
   WRITE(LUPRI,"(I3,'Node times:',/,20(10F7.2,/))") MYNUM,&
        & (tms(i),i=1,NNODES)
#else
   call quit('fix code in mm_getj_collect')
#endif
END SUBROUTINE mm_getj_collect
#endif /* ifdef VAR_MPI */
