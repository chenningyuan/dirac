! Notes: traditional dalton's WORK memory is used for some of the
! larger allocations before memory allocation in dalton is "fixed".
!
#if defined(SYS_AIX) || defined(VAR_PGF77)
   ! for some strange reason, this fails on AIX.
   #define ALLOCATION_WITHOUT_WORK 1
#endif

MODULE mm_interface_mod

   USE mm_global_paras_mod
   USE mm_stats_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_connect_interface,        &
             mm_disconnect_interface,     &
             mm_get_raw_data,             &
             mm_check_moms

   ! we use pointers here so we can allocate outside of the module
   TYPE(raw_mm_paras), POINTER, SAVE :: raw_paras(:)
   TYPE(J_index_type), POINTER, SAVE :: J_indices(:)
   REAL(REALK),        POINTER, SAVE :: dens(:)
   REAL(REALK),        POINTER, SAVE :: raw_qlm(:,:)
   ! Number of multipole moments generated externally
   TYPE(mm_counters), SAVE :: n_mms 

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_connect_interface(scheme,WORK,KWRK,LWORK)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(INOUT) :: scheme
      INTEGER(INTK),      INTENT(IN)    :: LWORK
      INTEGER(INTK),      INTENT(INOUT) :: KWRK
      REAL(REALK),        INTENT(INOUT) :: WORK(LWORK)
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      INTEGER ierr
#endif

      NULLIFY (dens, raw_qlm, raw_paras, J_indices)
      IF ( .NOT. MM_STANDALONE ) CALL mm_GTUNIT(LUPRI)
      CALL mm_init_moments(scheme%raw_LMAX,WORK,KWRK,LWORK)
#if defined (VAR_MPI)
      if(mynum==0) CALL mm_read_in_raw_data
      CALL mm_raw_data_sync
#else
      CALL mm_read_in_raw_data
#endif
!      CALL mm_check_moms(scheme%raw_LMAX,raw_qlm)
      CALL mm_translate_centres(scheme%system_size)

   END SUBROUTINE mm_connect_interface

!-------------------------------------------------------------------------------

   SUBROUTINE mm_disconnect_interface

      IMPLICIT NONE
      n_mms%elec = 0
      n_mms%nuc  = 0
      n_mms%tot  = 0
      CALL mm_free_moments

   END SUBROUTINE mm_disconnect_interface

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_moments(LMAX,WORK,KWRK,LWORK)
 
      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)    :: LMAX, LWORK
      INTEGER(INTK), INTENT(INOUT) :: KWRK
      REAL(REALK),   INTENT(INOUT) :: WORK(LWORK)
      integer*4 :: itemp4_lwlm


#ifdef VAR_MPI
      INCLUDE 'mpif.h'
      INTEGER ierr
      if(mynum==0) CALL mm_get_n_mms_from_file(LMAX)
      CALL MPI_Bcast(n_mms%nuc, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(n_mms%elec,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#else
      CALL mm_get_n_mms_from_file(LMAX)
#endif
      !CALL mm_allocate_mms_arrays(LMAX,WORK,LWORK/(LMAX+1)**2)
      itemp4_lwlm=LWORK/(LMAX+1)**2
      CALL mm_allocate_mms_arrays(LMAX,WORK,itemp4_lwlm)
      KWRK = KWRK + ((LMAX+1)**2)*n_mms%tot
   END SUBROUTINE mm_init_moments

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_n_mms_from_file(LMAX_in)
     USE density_fitting
      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: LMAX_in
      INTEGER(INTK) :: LMAX, NBAST, ID, LUINTM
      REAL(REALK)   :: RD

      LUINTM = 0
      CALL GPOPEN(LUINTM,INPUT_FILE2,'UNKNOWN',' ','UNFORMATTED',ID,.FALSE.)
      REWIND (LUINTM)
      READ (LUINTM) LMAX, NBAST, n_mms%elec, n_mms%nuc, ID, &
                    fit%NAUXMOM, fit%LHS_AUX, fit%RHS_AUX

      fit%AUX = fit%LHS_AUX .OR. fit%RHS_AUX
      n_mms%tot = n_mms%elec + n_mms%nuc 
      stat_raw_moms = n_mms%tot
      stat_nuc_moms = n_mms%nuc

      IF (LMAX /= LMAX_in) CALL QUIT('LMAX inconsistency in MM interface!')
      CALL GPCLOSE(LUINTM,'KEEP')
   END SUBROUTINE mm_get_n_mms_from_file

!-------------------------------------------------------------------------------

   SUBROUTINE mm_allocate_mms_arrays(LMAX,WORK,maxmom)

      USE mm_memory_manager_mod, ONLY: mm_allocate

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: LMAX, maxmom
      REAL(REALK),   TARGET     :: WORK((LMAX+1)**2,maxmom)
      INTEGER(INTK) :: i

      CALL mm_allocate(MEM_RAW_QLM,raw_paras,n_mms%tot)
      ALLOCATE( J_indices(n_mms%tot) )
      ! initialise parameters
      DO i = 1, n_mms%tot 
         raw_paras(i)%cntr = zero
         raw_paras(i)%ext = zero
         raw_paras(i)%id = 0
         raw_paras(i)%batch = 0
         raw_paras(i)%Lmin = 0
         raw_paras(i)%map_up = 0
         raw_paras(i)%box = 0 
         raw_paras(i)%bra = 0
         raw_paras(i)%box_cntr = zero
         J_indices(i)%i_indx = 0
         J_indices(i)%j_indx = 0
      END DO

      CALL mm_allocate(MEM_RAW_QLM,dens,n_mms%tot)
      ! use work for the stuff below:
#ifdef ALLOCATION_WITHOUT_WORK
      CALL mm_allocate(MEM_RAW_QLM,raw_qlm,(LMAX+1)**2,n_mms%tot)
#else
      !WRITE(LUPRI,*)'WORK Memory needed for ',n_mms%tot,' moments.'
      !WRITE(LUPRI,*) maxmom,' =',maxmom*((LMAX+1)**2),' available.'
      IF(maxmom<n_mms%tot) THEN
         WRITE(LUPRI,*)'WORK Memory needed for ',n_mms%tot,' moments.'
         WRITE(LUPRI,*)'require ',n_mms%tot*(LMAX+1)**2,' words.'
         CALL QUIT('Not enough memory for moments')
      END IF
      raw_qlm => WORK(1:(LMAX+1)**2,1:n_mms%tot)
#endif

      raw_qlm(:,:) = zero  ! crucial since only non-zero written explicitly

   END SUBROUTINE mm_allocate_mms_arrays

!------------------------------------------------------------------------------

   SUBROUTINE mm_free_moments

      IMPLICIT NONE

      IF (.NOT.ASSOCIATED(dens))      STOP "dens should be allocated!"
      IF (.NOT.ASSOCIATED(raw_qlm))   STOP "raw_qlm should be allocated!"
      IF (.NOT.ASSOCIATED(raw_paras)) STOP "raw_paras should be allocated!"
      IF (.NOT.ASSOCIATED(J_indices)) STOP "J_indices should be allocated!"
#ifdef ALLOCATION_WITHOUT_WORK
      IF (ASSOCIATED(raw_qlm)) DEALLOCATE(raw_qlm)
#endif
      IF (ASSOCIATED(dens)) DEALLOCATE(dens)
      IF (ASSOCIATED(raw_paras)) DEALLOCATE(raw_paras)
      IF (ASSOCIATED(J_indices)) DEALLOCATE(J_indices)
      NULLIFY (raw_paras, dens, raw_qlm, J_indices)

   END SUBROUTINE mm_free_moments

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_raw_data(scheme,LHS_mms,RHS_mms)
     USE mm_global_paras_mod
     USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(OUT) :: LHS_mms, RHS_mms
      INTEGER(INTK) :: lo,hi

      ! get LHS range multipole data
      lo = 1
      hi = n_mms%tot
      IF (scheme%LHS_mm_range == NUCLEAR_ONLY)    lo = 1+n_mms%elec
!orig      IF (scheme%LHS_mm_range == ELECTRONIC_ONLY) hi = n_mms%elec
      IF (scheme%LHS_mm_range == ELECTRONIC_ONLY) hi = n_mms%elec - fit%NAUXMOM
      IF (fit%LHS_AUX) THEN
         lo = n_mms%elec - fit%NAUXMOM + 1
         hi = n_mms%elec
      END IF
      LHS_mms%qlm       => raw_qlm(:,lo:hi)
      LHS_mms%dens      => dens(lo:hi)
      LHS_mms%paras     => raw_paras(lo:hi)
      LHS_mms%J_indices => J_indices(lo:hi)
!CAREFUL      stat_raw_moms_LHS = (hi-lo+1)
      ! note that now paras:moments mapping is simply by array position

      ! get RHS range multipole data
      lo = 1
      hi = n_mms%tot
      IF (scheme%RHS_mm_range == NUCLEAR_ONLY)    lo = 1+n_mms%elec
!orig      IF (scheme%RHS_mm_range == ELECTRONIC_ONLY) hi = n_mms%elec
      IF (scheme%RHS_mm_range == ELECTRONIC_ONLY) hi = n_mms%elec - fit%NAUXMOM
      IF (fit%RHS_AUX) THEN
         lo = n_mms%elec - fit%NAUXMOM + 1
         hi = n_mms%elec
      END IF
      RHS_mms%qlm       => raw_qlm(:,lo:hi)
      RHS_mms%dens      => dens(lo:hi)
      RHS_mms%paras     => raw_paras(lo:hi)
      RHS_mms%J_indices => J_indices(lo:hi)

!cSimen
!      CALL HEADER('LHS_mms%qlm in mm-driver',-1)
!      CALL OUTPUT(LHS_mms%qlm,1,121,1,hi-lo+1,121,hi-lo+1,1,2)
!      CALL HEADER('LHS_mms%dens in mm-driver',-1)
!      CALL OUTPUT(LHS_mms%dens,1,1,1,hi-lo+1,1,hi-lo+1,1,2)
!      CALL HEADER('LHS_mms%paras in mm-driver',-1)
!      CALL OUTPUT(LHS_mms%paras,1,1,1,hi-lo+1,1,hi-lo+1,1,2)
!      CALL HEADER('LHS_mms%J_indices in mm-driver',-1)
!      CALL OUTPUT(LHS_mms%J_indices,1,1,1,hi-lo+1,1,hi-lo+1,1,2)
!     CALL OUTPUT(RHS_mms%qlm,1,121,1,hi-lo+1,121,hi-lo+1,1,6)
!      CALL HEADER('CFIT in mm-driver',-1)
!      CALL OUTPUT(RHS_mms%dens,1,hi-lo+1,1,1,1,hi-lo+1,1,2)
!     CALL OUTPUT(RHS_mms%paras,1,1,1,hi-lo+1,1,hi-lo+1,1,6)
!     CALL OUTPUT(RHS_mms%J_indices,1,1,1,hi-lo+1,1,hi-lo+1,1,6)
!cSimen


!CAREFUL      stat_raw_moms_RHS = (hi-lo+1)
      ! note that now paras:moments mapping is simply by array position

   END SUBROUTINE mm_get_raw_data

!-------------------------------------------------------------------------------
! Read in multipole moment data from interface file.
! In all this MM code we assume the order of moments is:
!   (0),(-1,0,1),(-2,-1,0,1,2)...

   SUBROUTINE mm_read_in_raw_data

      IMPLICIT NONE

      REAL(REALK)   :: EXTP,PX,PY,PZ, SPH,DNS
      INTEGER(INTK) :: L,M, A,B, LMINA,LMINB, IBTCH, I, p
      INTEGER(INTK) :: ID, LUINTM

      LUINTM = 0
      CALL GPOPEN(LUINTM,INPUT_FILE1,'UNKNOWN',' ','UNFORMATTED',ID,.FALSE.)
!      OPEN(UNIT=LUINTM,FILE=INPUT_FILE,FORM='UNFORMATTED')

      ! first read electronic multipole moments into core
      REWIND (LUINTM)
      readloop: DO
         READ (LUINTM) L,M, A,B, LMINA,LMINB,IBTCH,I,EXTP,PX,PY,PZ,SPH,DNS
         IF (L .EQ. -1) EXIT readloop
         IF (I > SIZE(raw_qlm,2)) CALL QUIT ('interface file error 1')
         IF (((L+1)**2) > SIZE(raw_qlm,1)) CALL QUIT ('interface file error 2')
          ! indices to map moments to orbitals in J-matrix
         J_indices(I)%i_indx = A
         J_indices(I)%j_indx = B
          ! index to batches of unique centre and extent
         raw_paras(I)%batch  = IBTCH
          ! minimum order for exact treatment (if uncontracted basis set)
         raw_paras(I)%Lmin = LMINA + LMINB
         raw_paras(I)%cntr = (/ PX, PY, PZ /)
         raw_paras(I)%ext  = EXTP
          ! see defn p.424 "Electronic Structure Theory", Helgaker et al
         p = L*(L+1) +M +1 
         dens(I) = DNS
          ! components (l,m) of MM without density factorised in 
         raw_qlm(p,I) = SPH
      END DO readloop

      ! next read nuclei data: charge, and location
      DO I = 1, n_mms%nuc 
         p = n_mms%elec + I
         READ (LUINTM) raw_qlm(1,p), raw_paras(p)%cntr(:)
         raw_paras(p)%batch = raw_paras(n_mms%elec)%batch + I
      END DO
      raw_paras((n_mms%elec+1):)%Lmin = 0       ! spherical charge
      raw_paras((n_mms%elec+1):)%ext = zero     ! point charges
      raw_paras((n_mms%elec+1):)%ext = zero     ! point charges
      J_indices((n_mms%elec+1):)%i_indx = 0     ! not relevant
      J_indices((n_mms%elec+1):)%j_indx = 0     ! not relevant

      CALL GPCLOSE(LUINTM,'KEEP')
!      CLOSE(UNIT=LUINTM,STATUS='KEEP')
      dens((n_mms%elec+1):) = one

   END SUBROUTINE mm_read_in_raw_data

!-------------------------------------------------------------------------------
#ifdef VAR_MPI
   ! broadcast MM input data to all the involved nodes.
   SUBROUTINE mm_raw_data_sync

      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER master,ierr,RAW_PARAS_TYPE

      master = 0
      ! pass raw_paras
      CALL get_raw_paras_type(RAW_PARAS_TYPE)
      CALL MPI_Bcast(raw_paras, SIZE(raw_paras), RAW_PARAS_TYPE,&
           &         master,MPI_COMM_WORLD,IERR)
      CALL mpi_type_free(RAW_PARAS_TYPE)

      ! pass density
      CALL MPI_Bcast(dens,SIZE(dens),MPI_DOUBLE_PRECISION,&
           &         master,MPI_COMM_WORLD,IERR)

      ! pass nuclear moment data
      CALL MPI_Bcast(raw_qlm,SIZE(raw_qlm),MPI_DOUBLE_PRECISION,&
           &         master,MPI_COMM_WORLD,IERR)
   END SUBROUTINE mm_raw_data_sync

      ! fixme: This definition should go along TYPE(raw_mm_paras)
   SUBROUTINE get_raw_paras_type(tp)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER, INTENT(OUT) :: tp

      INTEGER intEx, doubleEx, ierr
      INTEGER, DIMENSION(3) :: bllen, types, disps
      CALL MPI_TYPE_EXTENT(MPI_INTEGER,intex,ierr)
      CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,doubleEX,ierr)
      bllen = (/ 4, 9, 3 /)
      types = (/ MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_DOUBLE_PRECISION /)
      disps(1) = 0; 
      disps(2) = disps(1) + bllen(1)*doubleEx
      disps(3) = disps(2) + bllen(2)*intEx

      CALL MPI_Type_Struct(SIZE(bllen),bllen,disps,types,tp, ierr)
      CALL MPI_TYPE_COMMIT(tp,ierr)
   END SUBROUTINE get_raw_paras_type
#endif   
!-------------------------------------------------------------------------------
! Routine to translate co-ordinates so all centres are at x,y,z >= 0
! Hence all box indices will be positive integers
!FIXME; is there an issue when we have ZERO SYSTEM SIZE???

   SUBROUTINE mm_translate_centres(system_size)

      IMPLICIT NONE
      REAL(REALK), INTENT(OUT) :: system_size 

      REAL(REALK)   :: sys_min(3), sys_max(3)
      INTEGER(INTK) :: i

      ! determine system size
      sys_min = raw_paras(1)%cntr
      sys_max = raw_paras(1)%cntr
      DO i = 1, SIZE(raw_paras)
         sys_min(:) = MIN(sys_min(:),raw_paras(i)%cntr(:))
         sys_max(:) = MAX(sys_max(:),raw_paras(i)%cntr(:))
      END DO
      system_size = MAXVAL(sys_max - sys_min)

      DO i = 1, SIZE(raw_paras)
         raw_paras(i)%cntr = raw_paras(i)%cntr - sys_min
      END DO

   END SUBROUTINE mm_translate_centres

!-------------------------------------------------------------------------------

   SUBROUTINE mm_check_moms(LMAX,qlm)

      IMPLICIT NONE

      REAL(REALK), INTENT(IN) :: qlm(:,:)
      INTEGER(INTK), INTENT(IN) :: LMAX
      INTEGER(INTK) :: i,j, lm
      REAL(REALK) :: tmp

      write(lupri,*) "length of L,M array  =", SIZE(qlm,1)
      write(lupri,*) "length of mom array  =", SIZE(qlm,2)
      write(lupri,*) "length of dens array =", SIZE(dens)
      write(lupri,*) "number of electronic =", n_mms%elec

      IF ((1+LMAX)**2 /= SIZE(qlm,1)) CALL QUIT('error in mm_check_moms')
      lm = 0
      DO i = 0, LMAX
         DO j = -i, i
            lm = lm+1
            tmp = (MAXVAL(ABS(qlm(lm,:))))
            IF ( tmp > zero ) write(*,'(2I4,2E25.15)')   &
                    i,j, MAXVAL(ABS(qlm(lm,:))), dens(lm)
         END DO
      END DO
      write(*,*)

!      DO i = 1, SIZE(qlm,2)
!         print *, i
!         print '(E25.15)', qlm(:,i)
!      END DO
!      write(*,*)

   END SUBROUTINE mm_check_moms

!-------------------------------------------------------------------------------

END MODULE mm_interface_mod

!===============================================================================

MODULE mm_qlm_processor

   USE mm_global_paras_mod
   USE mm_stats_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_renormalise_qlm,      &
             mm_factor_in_dens,       &
             mm_get_T_sym_qlm,        &
             mm_pack_raw_moments,     &
             mm_pack_raw_parameters,  &
             mm_analyse_qlm_paras
            
CONTAINS

!-------------------------------------------------------------------------------
   
   SUBROUTINE mm_renormalise_qlm(LMAX,qlm)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)    :: LMAX
      REAL(REALK),   INTENT(INOUT) :: qlm(:,:)

      INTEGER(INTK) :: i,L,M,p,pp
      REAL(REALK)   :: pref

      ! prefactor to symmetrize T-matrix
      DO i = 1, SIZE(qlm,2)
         DO L = 0, LMAX
            pp = L*(L+1) +1
            DO M = -L, -1
               pref = -one/(SQRT(two*FACTORIAL(l-m)*FACTORIAL(l+m)))
               p = pp+M
               qlm(p,i) = pref*qlm(p,i)
            END DO
            pref = one/FACTORIAL(L)
            p = pp  ! M=0
            qlm(p,i) = pref*qlm(p,i)
            DO M = 1, L
               pref = ((-1)**m)/SQRT(two*FACTORIAL(l-m)*FACTORIAL(l+m))
               p = pp+M
               qlm(p,i) = pref*qlm(p,i)
            END DO
         END DO
      END DO

   CONTAINS

      REAL(REALK) FUNCTION FACTORIAL(n)
         IMPLICIT NONE
         INTEGER(INTK), INTENT(IN) :: n 
         INTEGER(INTK) :: i
         FACTORIAL = 1 
         DO i = n, 2, -1 
            FACTORIAL = FACTORIAL*i
         END DO
      END FUNCTION FACTORIAL

   END SUBROUTINE mm_renormalise_qlm

!-------------------------------------------------------------------------------

   SUBROUTINE mm_factor_in_dens(dens,qlm)

      IMPLICIT NONE
      REAL(REALK),   INTENT(IN)    :: dens(:)
      REAL(REALK),   INTENT(INOUT) :: qlm(:,:)

      INTEGER(INTK) :: i

!FIXME: nicer way to do this with SPREAD??
      DO i = 1, SIZE(qlm,2)
         qlm(:,i) = qlm(:,i)*dens(i)
      END DO

   END SUBROUTINE mm_factor_in_dens

!-------------------------------------------------------------------------------
! Prefactorising moments to symmetrize modified T-matrix

   SUBROUTINE mm_get_T_sym_qlm(LMAX,qlm_in,qlm_out)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: qlm_in(:,:)
      REAL(REALK),   INTENT(OUT) :: qlm_out(:,:)

      INTEGER(INTK) :: i,L,u, hi,lo
      REAL(REALK)   :: pref

      DO i = 1, SIZE(qlm_in,2)
         DO L = 0, LMAX
            u = L*(L+1) +1     ! m=0
            hi = u+L
            lo = u-L
            pref = two*((-1)**L)
            qlm_out(lo:hi,i) = pref*qlm_in(lo:hi,i)
            qlm_out(u,i) = half*pref*qlm_in(u,i)
         END DO
      END DO

   END SUBROUTINE mm_get_T_sym_qlm

!-------------------------------------------------------------------------------
! get number of unique batches;
! also check raw_paras are sorted by batch (expected from input file);
! do this by checking the batch ID is always increasing

   SUBROUTINE get_nbatch(paras,nbatch)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)  :: paras(:)
      INTEGER(INTK),      INTENT(OUT) :: nbatch
      INTEGER(INTK) :: i, ndim

      ndim = SIZE(paras)
      nbatch = 1
      DO i = 2, ndim
         IF ( paras(i)%batch < paras(i-1)%batch ) THEN
            CALL QUIT('batches input from source file not sorted!')
         END IF
         ! some batch numbers are skipped so need to take care of this
         IF ( paras(i)%batch /= paras(i-1)%batch ) THEN
            nbatch = nbatch + 1
         END IF
      END DO

   END SUBROUTINE get_nbatch

!-------------------------------------------------------------------------------

   SUBROUTINE get_pkd_data(add_dens,raw_mms,pkd_paras,pkd_qlm)

      IMPLICIT NONE
      LOGICAL,            INTENT(IN)  :: add_dens
      TYPE(raw_mm_data),  INTENT(IN)  :: raw_mms
      TYPE(raw_mm_paras), INTENT(OUT) :: pkd_paras(:) 
      REAL(REALK),        INTENT(OUT) :: pkd_qlm(:,:)
      INTEGER(INTK) :: i,j, last_batch

      j = 0
      last_batch = -1
      DO i = 1, SIZE(raw_mms%paras)
         IF ( raw_mms%paras(i)%batch == last_batch ) THEN
            ! element same batch as previous
            IF (add_dens) THEN
               pkd_qlm(:,j) = pkd_qlm(:,j) + raw_mms%qlm(:,i)*raw_mms%dens(i)
            ELSE
               pkd_qlm(:,j) = pkd_qlm(:,j) + raw_mms%qlm(:,i)
            END IF
            ! need to store most conservative Lmin (order of exact expansion)
            pkd_paras(j)%Lmin = MAX(pkd_paras(j)%Lmin, raw_mms%paras(i)%Lmin) 
         ELSE
            ! element in new batch
            j = j + 1
            pkd_paras(j) = raw_mms%paras(i)
            pkd_paras(j)%id = j
            IF (add_dens) THEN
               pkd_qlm(:,j) = raw_mms%qlm(:,i)*raw_mms%dens(i)
            ELSE
               pkd_qlm(:,j) = raw_mms%qlm(:,i)
            END IF
         END IF
         last_batch = raw_mms%paras(i)%batch
      END DO

   END SUBROUTINE get_pkd_data

!-------------------------------------------------------------------------------

   SUBROUTINE get_screened_nmom(screen_thr,qlm,skip,nskip)

      IMPLICIT NONE
      REAL(REALK),   INTENT(IN)  :: screen_thr
      REAL(REALK),   INTENT(IN)  :: qlm(:,:)
      LOGICAL,       INTENT(OUT) :: skip(:)
      INTEGER(INTK), INTENT(OUT) :: nskip
      INTEGER(INTK) :: i,j

      nskip = 0
      batches: DO i = 1, SIZE(qlm,2)
         skip(i) = .TRUE.
         nskip = nskip + 1
         lm_loop: DO j = 1, SIZE(qlm,1)
            IF ( ABS(qlm(j,i)) > screen_thr ) THEN
               skip(i) = .FALSE.
               nskip = nskip - 1
               EXIT lm_loop
            END IF
         END DO lm_loop
      END DO batches

   END SUBROUTINE get_screened_nmom

!-------------------------------------------------------------------------------
! here we squeeze all the significant batches of moments to be sequential
! at the top of the array, with all the insignificant moments overwritten;
! we use skip(:) as a logical mask to direct the skipping.

   SUBROUTINE squeeze_significant_batches(skip,paras,qlm)

      IMPLICIT NONE
      LOGICAL,            INTENT(IN)    :: skip(:)
      TYPE(raw_mm_paras), INTENT(INOUT) :: paras(:) 
      REAL(REALK),        INTENT(INOUT) :: qlm(:,:)
      INTEGER(INTK) :: i, j

      IF (SIZE(paras) /= SIZE(qlm,2)) STOP 'paras and qlm should be same size!'
      IF (SIZE(paras) /= SIZE(skip)) STOP 'paras and skip should be same size!'

      j = 0
      DO i = 1, SIZE(paras)
         IF (skip(i)) CYCLE
         j = j + 1
         paras(j) = paras(i)
         qlm(:,j) = qlm(:,i)
      END DO

   END SUBROUTINE squeeze_significant_batches

!-------------------------------------------------------------------------------
! routine to drive the packing of a set of raw moments by batches,
! where members of a batch share a common centre and extent,
! including density factoring and screening if requested.
! note that no record is kept of the packing for later "unpacking".

   SUBROUTINE mm_pack_raw_moments(raw_mms,dens,dens_thr,pkd_paras,pkd_qlm,ndim)

      USE mm_memory_manager_mod, ONLY: mm_allocate

      IMPLICIT NONE
      TYPE(raw_mm_data),  INTENT(IN)  :: raw_mms
      LOGICAL,            INTENT(IN)  :: dens
      REAL(REALK),        INTENT(IN)  :: dens_thr
      TYPE(raw_mm_paras), POINTER     :: pkd_paras(:) 
      REAL(REALK),        POINTER     :: pkd_qlm(:,:)
      INTEGER(INTK),      INTENT(OUT) :: ndim
      INTEGER(INTK) :: nbatch, nskip
      LOGICAL, ALLOCATABLE :: skip(:)

      ! get number of unique batches;
      CALL get_nbatch(raw_mms%paras,nbatch)

      ! now build packed data by summing raw data in the same batch
      CALL mm_allocate(MEM_RAW_QLM,pkd_paras,nbatch)
!miro-------------------------------------------------------------------------
!      CALL mm_allocate(MEM_RAW_QLM,pkd_qlm,SIZE(raw_mms%qlm,1),nbatch)
!   Error: There is no specific subroutine for the generic 'mm_allocate' at (1)
      call quit('error in compilation, discarded !')

      CALL get_pkd_data(dens,raw_mms,pkd_paras,pkd_qlm)
      ndim = nbatch

      IF ( dens .AND. (dens_thr > DENS_SCREEN_CUT) ) THEN
         !WRITE(LUPRI,'(A,E9.1)') ' MM RHS density screening on:', dens_thr
         ! perform density-based screening;
         ! note we do not reallocate whole array, but push all the
         ! significant terms to the top, and only this array section
         ! of significant moments is then pointed to.
         ALLOCATE(skip(nbatch))
         CALL get_screened_nmom(dens_thr,pkd_qlm,skip,nskip)
         ndim = nbatch - nskip  ! number of significant moments
         CALL squeeze_significant_batches(skip,pkd_paras,pkd_qlm)
         DEALLOCATE(skip)
      END IF

      ! this is just for statistics
      stat_pkd_moms_RHS = nbatch
      stat_screened_moms_RHS = ndim

   END SUBROUTINE mm_pack_raw_moments

!-------------------------------------------------------------------------------
! routine to drive the packing of a set of raw mm parameters by batches,
! where members of a batch share a common centre and extent,
! including the build of a mapping between the packed paras and the
! original raw paras.  This map can be used for later "unpacking".

   SUBROUTINE mm_pack_raw_parameters(raw_paras,pkd_paras,batch_map)

      USE mm_memory_manager_mod, ONLY: mm_allocate

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN) :: raw_paras(:) 
      TYPE(raw_mm_paras), POINTER    :: pkd_paras(:) 
      TYPE(id_list),      POINTER    :: batch_map(:) 
      INTEGER(INTK) :: i,j, nbatch, last_batch

      ! get number of unique batches;
      CALL get_nbatch(raw_paras,nbatch)
      stat_pkd_moms_LHS = nbatch

      ! initialise packed paras and batch map
      CALL mm_allocate(MEM_RAW_QLM,pkd_paras,nbatch)
      ALLOCATE(batch_map(nbatch))
      DO i = 1, nbatch 
         batch_map(i)%occ = 0
         NULLIFY( batch_map(i)%head )
      END DO

      ! now build packed paras by compressing raw paras in same batch
      j = 0
      last_batch = -1
      DO i = 1, SIZE(raw_paras)
         IF ( raw_paras(i)%batch == last_batch ) THEN
            ! add raw parameter mapping to existing linked-list for this batch
            CALL add_batch_item(batch_map(j),i)
            ! need to store most conservative Lmin (order of exact expansion)
            pkd_paras(j)%Lmin = MAX(pkd_paras(j)%Lmin, raw_paras(i)%Lmin) 
         ELSE
            ! element in new batch
            j = j + 1
            pkd_paras(j) = raw_paras(i)
            pkd_paras(j)%id = j
            ! linked-list for this batch is empty, so start one
            batch_map(j)%occ = 1
            ALLOCATE( batch_map(j)%head )
            batch_map(j)%head%id = i           ! map back to raw paras here
            NULLIFY( batch_map(j)%head%next )  ! rest of list is empty
         END IF
         last_batch = raw_paras(i)%batch
      END DO

   CONTAINS

      SUBROUTINE add_batch_item(batch_list,raw_id)

         IMPLICIT NONE
         TYPE(id_list), INTENT(INOUT) :: batch_list
         INTEGER(INTK), INTENT(IN)    :: raw_id
         TYPE(id_node), POINTER :: new_node

         batch_list%occ = batch_list%occ + 1
         ALLOCATE( new_node )
         new_node%id = raw_id
         IF (ASSOCIATED(batch_list%head%next)) THEN
            ! more than one entry in list (including head)
            ! so point new_node to old second entry
            new_node%next => batch_list%head%next
            ! point head to new_node
            NULLIFY(batch_list%head%next)
            batch_list%head%next => new_node
         ELSE
            ! only head so far; make new_node our second entry
            batch_list%head%next => new_node
            NULLIFY(new_node%next)   ! end of list
         END IF

      END SUBROUTINE add_batch_item

   END SUBROUTINE mm_pack_raw_parameters

!-------------------------------------------------------------------------------

   SUBROUTINE mm_analyse_qlm_paras(paras)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN) :: paras(:) 
      REAL(REALK)   :: centre(3)
      INTEGER(INTK) :: i,j, ncen

      ncen = 1
      DO i = 2, SIZE(paras) 
         ncen = ncen +1
         centre = paras(i)%cntr
         jloop: DO j = i-1, 1, -1
            IF ( centre(1) == paras(j)%cntr(1) ) THEN
               IF ( centre(2) == paras(j)%cntr(2) ) THEN
                  IF ( centre(3) == paras(j)%cntr(3) ) THEN
                     ncen = ncen -1
                     EXIT jloop
                  END IF
               END IF
            END IF
         END DO jloop
      END DO
      WRITE(LUPRI,'(A,I8)') ' distinct centres:', ncen

   END SUBROUTINE mm_analyse_qlm_paras

!-------------------------------------------------------------------------------

END MODULE mm_qlm_processor

!===============================================================================

MODULE mm_aux_qlm_builder

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_init_qlm_builder,    &
             mm_free_qlm_builder,    &
             mm_get_aux_qlm
            
   ! we use pointers here so we can allocate outside of the module
   TYPE(raw_mm_paras), POINTER, SAVE :: LHS_paras(:) 
   TYPE(id_list),      POINTER, SAVE :: LHS_batch_map(:) 
   TYPE(raw_mm_paras), POINTER, SAVE :: RHS_paras(:) 
   REAL(REALK),        POINTER, SAVE :: RHS_qlm_T(:,:)
   REAL(REALK),        POINTER, SAVE :: RHS_qlm_W(:,:)

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_qlm_builder

      IMPLICIT NONE

      NULLIFY(LHS_paras, LHS_batch_map)
      NULLIFY(RHS_paras, RHS_qlm_T, RHS_qlm_W)

   END SUBROUTINE mm_init_qlm_builder

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_qlm_builder

      IMPLICIT NONE
      INTEGER(INTK) :: i

      IF (ASSOCIATED(LHS_paras)) DEALLOCATE(LHS_paras)
      IF (ASSOCIATED(LHS_paras)) DEALLOCATE(LHS_batch_map)
      IF (ASSOCIATED(RHS_paras)) DEALLOCATE(RHS_paras)
      IF (ASSOCIATED(RHS_qlm_T)) DEALLOCATE(RHS_qlm_T)
      IF (ASSOCIATED(RHS_qlm_W)) DEALLOCATE(RHS_qlm_W)

      IF (ASSOCIATED(LHS_batch_map)) THEN
         DO i = 1, SIZE(LHS_batch_map)
            CALL free_batch_map(LHS_batch_map(i)%head)
         END DO
      END IF

      NULLIFY(LHS_paras, LHS_batch_map)
      NULLIFY(RHS_paras, RHS_qlm_T, RHS_qlm_W)

   CONTAINS 

      RECURSIVE SUBROUTINE free_batch_map(node)

         IMPLICIT NONE
         TYPE(id_node), POINTER :: node

         IF (ASSOCIATED(node%next)) THEN
            CALL free_batch_map(node%next)
         END IF
         DEALLOCATE(node)
         NULLIFY(node)

      END SUBROUTINE free_batch_map

   END SUBROUTINE mm_free_qlm_builder

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_aux_qlm(scheme,LHS_mms,RHS_mms)
      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: LHS_mms, RHS_mms

       ! we use the SCALED solid harmonic formulation throughout
      IF (USE_UNSCALED_HARMONICS) CALL QUIT ('unscaled harmonics not coded!')
      CALL get_normalised_qlm(scheme,LHS_mms%qlm,RHS_mms%qlm)

       ! get RHS parameter mappings and preconditioned W and T matrices first
      CALL get_RHS_data(scheme,RHS_mms)
       ! get LHS parameter mappings and preconditioned W and T matrices next;
       ! note we get LHS after RHS because LHS moments may
       ! alter the raw moments subsequently used by RHS
      CALL get_LHS_data(scheme,LHS_mms)

      ! all relevant data should now be held in %qlm_T and %qlm_W
      NULLIFY (LHS_mms%qlm, RHS_mms%qlm)

   END SUBROUTINE mm_get_aux_qlm

!-------------------------------------------------------------------------------

   SUBROUTINE get_normalised_qlm(scheme,LHS_qlm,RHS_qlm)
      USE density_fitting
      USE mm_qlm_processor, ONLY: mm_renormalise_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      REAL(REALK),        INTENT(INOUT) :: LHS_qlm(:,:), RHS_qlm(:,:)
      INTEGER(INTK) :: LMAX, i, qlm_dim

      LMAX = scheme%raw_LMAX

      IF (scheme%LHS_mm_range == ALL_MOMENTS) THEN
         ! all data renormalised via LHS
         CALL mm_renormalise_qlm(LMAX,LHS_qlm) 
      ELSE IF (scheme%RHS_mm_range == ALL_MOMENTS) THEN
         ! all data renormalised via RHS
         CALL mm_renormalise_qlm(LMAX,RHS_qlm) 
      ELSE
         CALL mm_renormalise_qlm(LMAX,LHS_qlm) 
!CAREFUL affects the results if renormalized twice!!!!!!!!!!!
         IF ((scheme%RHS_mm_range /= scheme%LHS_mm_range) .OR. fit%AUX) THEN 
!         IF (scheme%RHS_mm_range /= scheme%LHS_mm_range) THEN 
            ! LHS and RHS are disjoint (must avoid scaling some twice!)
            CALL mm_renormalise_qlm(LMAX,RHS_qlm) 
         END IF
      END IF

   END SUBROUTINE get_normalised_qlm

!-------------------------------------------------------------------------------

   SUBROUTINE get_LHS_data(scheme,LHS_mms)

      USE mm_memory_manager_mod, ONLY: mm_allocate
      USE mm_qlm_processor,  ONLY: mm_factor_in_dens,      &
                                   mm_pack_raw_parameters
      USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: LHS_mms
      INTEGER(INTK) :: i, qlm_dim

       ! get LHS parameters from raw data
      IF (scheme%pack_LHS) THEN
         CALL mm_pack_raw_parameters(LHS_mms%paras,LHS_paras,LHS_batch_map)
         qlm_dim = SIZE(LHS_paras)
         LHS_mms%batch_map => LHS_batch_map(:)
      ELSE
         qlm_dim = SIZE(LHS_mms%qlm,2)
         ALLOCATE(LHS_paras(qlm_dim))
!         CALL mm_allocate(MEM_RAW_QLM,LHS_paras,qlm_dim)
         LHS_paras = LHS_mms%paras(:)
      END IF

       ! set up (packed) LHS parameters and paras:moment mappings
      ! the compact expression is:
      ! LHS_paras(:)%id = (/ (i,i=1,qlm_dim) /)
      ! but intel fortran compiler 7.1-generated code may deadlock
      ! for certain vector lengths (in particular, for insulin molecule).
      ! let's play safe here.
      DO i = 1, qlm_dim
         LHS_paras(i)%id = i
      END DO
      LHS_mms%paras => LHS_paras(:)

       ! LHS T-matrix
      SELECT CASE (scheme%T_con%LHS_mm_type)
      CASE (USE_RAW_QLM)
         LHS_mms%qlm_T => LHS_mms%qlm(:,:)
      CASE DEFAULT
         CALL QUIT('cannot reconcile LHS_mm_type') 
      END SELECT

       ! factorise in density if required
      IF ( scheme%LHS_dens ) THEN
         CALL mm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
         NULLIFY(LHS_mms%dens)  ! should no-longer be used
      END IF

       ! LHS W-matrix
      NULLIFY(LHS_mms%qlm_W) ! not actually used

   END SUBROUTINE get_LHS_data

!-------------------------------------------------------------------------------

   SUBROUTINE get_RHS_data(scheme,RHS_mms)

      USE mm_memory_manager_mod, ONLY: mm_allocate
      USE mm_qlm_processor,  ONLY: mm_get_T_sym_qlm,     &
                                   mm_pack_raw_moments,  &
                                   mm_factor_in_dens,    &
                                   mm_analyse_qlm_paras
      USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: RHS_mms
      INTEGER(INTK) :: LMAX, i, qlm_dim
      LOGICAL :: dens ! flags if density is factored in
      REAL(REALK) :: thr

      LMAX = scheme%raw_LMAX

       ! get RHS parameters and moments from raw data
      IF ( scheme%pack_RHS ) THEN
         ! we now pack all moments belonging to
         ! the same batch (same centre,extent);
         ! also factor in density and perform
         ! density-based screening if flagged
         dens = scheme%RHS_dens 
         thr = scheme%dens_screen_thr
         CALL mm_pack_raw_moments(RHS_mms,dens,thr,RHS_paras,RHS_qlm_W,qlm_dim)
      ELSE
         ! use unpacked moments and parameters (as above for LHS)
         qlm_dim = SIZE(RHS_mms%qlm,2)
         CALL mm_allocate(MEM_RAW_QLM,RHS_paras,qlm_dim)
         RHS_paras = RHS_mms%paras(:)
         ! RHS W-matrix
 

         !CALL mm_allocate(MEM_RAW_QLM,RHS_qlm_W,(1+LMAX)**2,qlm_dim)
         call quit('miro: discarded !')

         RHS_qlm_W = RHS_mms%qlm(:,:)
         IF ( scheme%RHS_dens ) THEN
            CALL mm_factor_in_dens(RHS_mms%dens,RHS_qlm_W)
            NULLIFY(RHS_mms%dens)  ! should no-longer be used
         END IF
      END IF

      ! note we must take care to only point to the required
      ! array section since screening may mean that the array is
      ! larger than the number of significant moments we want 
      
       ! set up (packed) RHS parameters and paras:moment mappings
      DO i = 1, qlm_dim
         RHS_paras(i)%id = i
      END DO
      RHS_mms%paras => RHS_paras(:qlm_dim)
      RHS_mms%qlm_W => RHS_qlm_W(:,:qlm_dim)

       ! precondition RHS T-matrix
      SELECT CASE (scheme%T_con%RHS_mm_type)
      CASE (USE_RAW_QLM)
         RHS_mms%qlm_T => RHS_qlm_W(:,:qlm_dim)
      CASE (USE_T_SYM_QLM)
        ! compilation error  There is no specific subroutine for the generic 'mm_allocate' 
         !CALL mm_allocate(MEM_RAW_QLM,RHS_qlm_T,(1+LMAX)**2,qlm_dim)
          call quit('miro: discarded!')
          ! build %qlm_T by rescaling significant %qlm_W
         CALL mm_get_T_sym_qlm(LMAX,RHS_mms%qlm_W,RHS_qlm_T)
         RHS_mms%qlm_T => RHS_qlm_T(:,:)
      CASE DEFAULT
         CALL QUIT('cannot reconcile RHS_mm_type') 
      END SELECT

!      CALL mm_analyse_qlm_paras(RHS_mms%paras)

   END SUBROUTINE get_RHS_data

!-------------------------------------------------------------------------------

END MODULE mm_aux_qlm_builder

