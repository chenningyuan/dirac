module aoresponse_main_m

  use kinds_m
  use scf_density_class

  implicit none

! type gao_data_t
!   integer(kind_i4)              :: max_derivative
!   integer(kind_i4), allocatable :: gao_center(:)
!   real(kind_dp),    allocatable :: gao(:)
!   real(kind_dp),    allocatable :: gao1(:)
!   real(kind_dp),    allocatable :: gao2(:)
!   real(kind_dp),    allocatable :: buffer(:)
! end type gao_data_t

  public open_response_driver, get_ao_pmat, solve_perturbed_ao_dmat

  private ! everything below

contains














! subroutine open_response_driver(work, kfree, lfree, ao_dmat_fo)
  subroutine open_response_driver(work, kfree, lfree)
#include "dcbbas.h"
#include "dgroup.h"
    real(kind_dp), intent(inout)         :: work(lfree)
    integer,       intent(inout)         :: kfree
    integer,       intent(inout)         :: lfree
!   real(kind_dp), intent(in)            :: ao_dmat_fo(ntbas(0)*ntbas(0)*nz)

    real(kind_dp),           allocatable :: pmat(:, :, :)
    real(kind_dp),    allocatable        :: ao_dmat_zo(:, :, :)
    real(kind_dp),    allocatable        :: mo_coef(:, :, :)
    real(kind_dp)                        :: ddot
    real(kind_dp)                        :: dip
    integer                              :: i

  

      call init_mo_coef(mo_coef)
      call get_mo_coef(mo_coef)

      call init_ao_dmat(ao_dmat_zo)
      call get_ao_dmat_zo(ao_dmat_zo, mo_coef)

      call del_mo_coef(mo_coef)



!   write(*, *) 'calling PRPXL1_clone'
!   call PRPXL1_clone(WORK,KFREE,LFREE)
 


    do i = 1, 10
       call init_ao_dmat(pmat)
       call get_ao_pmat(i, pmat, work, kfree, lfree)
       dip = ddot(ntbas(0)*ntbas(0)*nz, pmat, 1, ao_dmat_zo, 1)
!      dip = ddot(ntbas(0)*ntbas(0)*nz, pmat, 1, ao_dmat_fo, 1)
!      call print_square_q_matrix(pmat, ntbas(0))
       write(*, *) 'blub', i, dip

       call del_ao_dmat(pmat)
    end do

    call del_ao_dmat(ao_dmat_zo)




  end subroutine open_response_driver

  subroutine get_ao_pmat(property_index, pmat, work, kfree, lfree)
    integer(kind_i4), intent(in)    :: property_index
    real(kind_dp),    intent(out)   :: pmat(:, :, :)
    real(kind_dp),    intent(inout) :: work(lfree)
    integer,          intent(inout) :: kfree
    integer,          intent(inout) :: lfree

    logical                         :: file_exists, is_square
    real(kind_dp)                   :: dummy_dp
    integer(kind_i4)                :: file_unit, dummy_i4, idimension
    real(kind_dp),      allocatable :: ptri(:)
    real(kind_dp),      allocatable :: op1int(:)
    logical,            allocatable :: first(:)

#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbxpr.h"

    is_square = .false.
    if (iprptim(property_index) == 0) then
      is_square = .true.
    end if

    if (is_square) then
      idimension = n2bbasx
    else
      idimension = nnbbasx
    end if

    allocate(ptri(nz*idimension))
    allocate(op1int(idimension))
    allocate(first(nz))

    file_unit = LU1INT

    file_exists = .false.
    inquire(file = 'AOPROPER', exist = file_exists)
    if (.not. file_exists) then
      call quit('file AOPROPER not found')
    end if

    open(file_unit,          &
         file = 'AOPROPER',  &
         status = 'unknown', &
         form = 'unformatted')

    rewind file_unit

    call prpex2(property_index,dummy_dp, &
                .false.,idimension,      &
                ptri,pmat,op1int,first,  &
                work,kfree,lfree,0)

    close(file_unit, status = 'keep')

    deallocate(ptri)
    deallocate(op1int)
    deallocate(first)
  end subroutine get_ao_pmat












  subroutine solve_perturbed_ao_dmat(property_index, ao_dmat, work, kfree, lfree)
    integer(kind_i4), intent(in)    :: property_index
!   real(kind_dp),    intent(out)   :: ao_dmat(:, :, :)
    real(kind_dp)                   :: ao_dmat
    real(kind_dp),    intent(inout) :: work(lfree)
    integer,          intent(inout) :: kfree
    integer,          intent(inout) :: lfree


#include "dcbxpr.h"
#include "dcbxrs.h"
#include "../prp/dcbxqr.h"
#include "dgroup.h"

    integer :: IOP
    integer :: ksave
    integer :: ifreq, factor
    integer :: half_length, max_full_length

    real(kind_dp),    allocatable :: eigenvalues_reduced(:), eigenvectors_reduced(:), convergence(:)
    real(kind_dp),    allocatable :: prop_gradient_c(:), prop_gradient_pp(:), prop_gradient_pn(:)
    integer(kind_i4), allocatable :: ibtyp(:), ibtyp_pointer_c(:), ibtyp_pointer_pp(:), ibtyp_pointer_pn(:)

    real(kind_dp),    allocatable :: response_vector_pp(:), response_vector_pn(:)


    real(kind_dp)                        :: frequencies(1)

      LOGICAL EXST,ZERO
      CHARACTER TTYP(-1:1)*2
      SAVE TTYP
      DATA TTYP /'T-','T0','T+'/
    frequencies(1) = 0.0d0
 



!     initialize symmetry independent variables
!     =========================================

      LINEQ = .TRUE.
      FKRMC = .FALSE.
      LSVCFG(1) = .TRUE.
      LSVCFG(2) = .TRUE.
      TKNORM = XQRNRM
      DIAGHE = XQRDIH
      IPRXRS = 0
      THCXRS = THCQR
      RESFAC = RESXQR
      MAXITR = ITRXQR
      NREDM  = MAXQRM
      N2REDM = MAXQRM*MAXQRM
      LOFFTY = 0
      CNVINT(1) = CNVXQR(1)
      CNVINT(2) = CNVXQR(2)
      ITRINT(1) = ITRIQR(1)
      ITRINT(2) = ITRIQR(2)
      INTDEF = INTXQR
      UNCOUP = .FALSE.
      STERNH = .FALSE.








      iop = 1


!     initialize configurational and orbital parameters
!     =================================================
  

         KSAVE  = KFREE

         JSYMOP = IPRPSYM(property_index)
         JTIMOP = IPRPTIM(property_index)
         JOPSY  = JBTOF(JSYMOP-1,1)
         NFREQ  = LQROP(IOP,2)
         NEXSIM = NFREQ
         NEXSTV = NFREQ
         NEXCNV = NFREQ
         NCRED  = 0
         NERED  = 0
         NPRED  = 0
         NZRED  = NCRED+NERED+NPRED
         STATIC = (NFREQ.EQ.1 .AND. QRFREQ(1,IOP).EQ. 0.0d0)

         CALL XRSPAR(XQR_INDSTR,XQR_SKIPEE,XQR_SKIPEP, &
                     WORK,KFREE,LFREE, 0)





!   get property gradient
!   =====================

    allocate(prop_gradient_c(nzconfq))
    allocate(prop_gradient_pp(nzxopeq))
    allocate(prop_gradient_pn(nzxoppq))


         CALL PAMGRD(property_index,prop_gradient_c, prop_gradient_pp, prop_gradient_pn, &
              JOPSY,WORK(KZCONF),WORK(KZXOPE),WORK(KZXOPP),      &
              NZCONF,NZXOPE,NZXOPP,WORK,KFREE,LFREE, 0)

!   TODO: check the norm of the property gradient
!         and if below threshold then skip things


!   solve response equation in the reduced space spanned by the trial vectors
!   =========================================================================

    allocate(ibtyp(2*NREDM))
    allocate(ibtyp_pointer_c(NREDM))
    allocate(ibtyp_pointer_pp(NREDM))
    allocate(ibtyp_pointer_pn(NREDM))

    allocate(eigenvalues_reduced(nexsim))
    allocate(eigenvectors_reduced(NREDM*NEXSIM))
    allocate(convergence(nfreq))

    call dcopy(nfreq, frequencies, 1, eigenvalues_reduced, 1)





         CALL XRSCTL(prop_gradient_c, prop_gradient_pp, prop_gradient_pn, &
              ibtyp, ibtyp_pointer_c, ibtyp_pointer_pp, ibtyp_pointer_pn, &
              convergence,eigenvalues_reduced,eigenvectors_reduced,   &
              WORK,KFREE,LFREE)

    deallocate(prop_gradient_c)
    deallocate(prop_gradient_pp)
    deallocate(prop_gradient_pn)

    deallocate(ibtyp_pointer_c)

    deallocate(convergence)






!   construct the pp and pn response vectors
!   ========================================

    nesim = 0
    npsim = 0
    ncsim = 0

    if (static) then
      factor = 1
    else
      factor = 2
    end if

    allocate(response_vector_pp(mzyee*nz))
    allocate(response_vector_pn(mzyep*nz))

    if (nzxope > 0) then

      nesim           = nfreq*factor
      half_length     = nzxope
      max_full_length = mzyee

      call construct_response_vector('opp', property_index,                     &
                                     response_vector_pp,                        &
                                     half_length, max_full_length,              &
                                     ibtyp, ibtyp_pointer_pp,                   &
                                     eigenvalues_reduced, eigenvectors_reduced, &
                                     work, kfree, lfree)
    end if

    if (nzxopp > 0) then

      npsim           = nfreq*factor
      half_length     = nzxopp
      max_full_length = mzyep

      call construct_response_vector('opn', property_index,                     &
                                     response_vector_pn,                        &
                                     half_length, max_full_length,              &
                                     ibtyp, ibtyp_pointer_pn,                   &
                                     eigenvalues_reduced, eigenvectors_reduced, &
                                     work, kfree, lfree)
    end if

    deallocate(response_vector_pp)
    deallocate(response_vector_pn)

    deallocate(ibtyp)
    deallocate(ibtyp_pointer_pp)
    deallocate(ibtyp_pointer_pn)
  
    deallocate(eigenvalues_reduced)
    deallocate(eigenvectors_reduced)


!   construct perturbed AO density matrices
!   =======================================




!     CALL SETGMT(ISYMOP,IHRMOP,IFCKOP,IBCVC,IBEVC,IBPVC,IBTYP,IPRXRS)
!
!     Calculate AO density matrices
!     -----------------------------
!
!     Note: reuse memory and save in FMO:
!        DVT(1..NCSIM), [DCX, DVX](1..NESIM), [DCX, DVX](1..NPSIM)
!     Only the necessary matrices are saved.
!
!     CALL DMOLITX(FMO,BEVEC,BPVEC,CMO,DVT,DV,IHRMOP,NDMAT,
!    &     WORK,KFREE,LFREE)





  end subroutine solve_perturbed_ao_dmat

  subroutine construct_response_vector(vector_type, property_index,               &
                                       response_vector_full,                      &
                                       half_length, max_full_length,              &
                                       ibtyp, ibtyp_pointer,                      &
                                       eigenvalues_reduced, eigenvectors_reduced, &
                                       work, kfree, lfree)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "dcbxrs.h"
#include "../prp/dcbxqr.h"
#include "dcbibn.h"
#include "dcborb.h"
#include "dgroup.h"

#include "dcbxpr.h"
    character(*), intent(in) :: vector_type
    real(kind_dp),    intent(inout) :: work(lfree)
    integer,          intent(inout) :: kfree
    integer,          intent(inout) :: lfree
!
    real(kind_dp),    allocatable :: response_vector_z(:, :, :)
    real(kind_dp),    allocatable :: response_vector_y(:, :, :)
    real(kind_dp)                 :: response_vector_full(:)
    real(kind_dp),    allocatable :: buffer(:)
!
    integer :: property_index, IBTYP,ibtyp_pointer, nbtyp, nbred, lubvec, lufile
    integer :: factor, half_length, max_full_length
    integer :: KVECZ, KVECY, KRSPVC, Kbuffer, KIVECS, nzypar, ifr, iz, i, irec
      PARAMETER ( D1=1.0D0, DM1=-1.0D0, DH=0.5D0 )
      CHARACTER FILNAM*6,FILRSP*8,TYP*1
      DIMENSION IBTYP(2,*), ibtyp_pointer(*), &
                eigenvalues_reduced(*),eigenvectors_reduced(NREDM,*)





    allocate(response_vector_z(half_length, nz, nfreq))
    allocate(response_vector_y(half_length, nz, nfreq))
    allocate(buffer(half_length*nz))



















    select case (vector_type)
      case ('opp')
        nbtyp = 1
        NBRED  = NERED
        LUBVEC = LUBOE
        FILNAM = 'PAMBOE'
        TYP    = 'E'
        LUFILE = LUQREE
        FILRSP = 'QRVEC.EE'
      case ('opn')
        nbtyp = 2
        NBRED  = NPRED
        LUBVEC = LUBOP
        FILNAM = 'PAMBOP'
        TYP    = 'P'
        LUFILE = LUQREP
        FILRSP = 'QRVEC.EP'
      case default
        call quit('error in construct_response_vector_sub, unknown vector type')
    end select

!
!     Construct response solution vectors by summing the trial
!     vectors multiplied by the eigenvectors in the reduced
!     system. Note that the trial vectors are of the form 
!     (Z+Y) and (Z-Y).
!
    write(*, *) 'blub 1'
      OPEN(LUBVEC,FILE=FILNAM,FORM='UNFORMATTED', &
           ACCESS='DIRECT',RECL=8*half_length*NZ,STATUS='OLD')
      CALL DZERO(response_vector_z,half_length*NZ*NFREQ)
      CALL DZERO(response_vector_y,half_length*NZ*NFREQ)
      CALL MEMGET('INTE',KIVECS,NFREQ,WORK,KFREE,LFREE)
      CALL XRSXV1(1,NBTYP,response_vector_z,eigenvectors_reduced,NFREQ,IBTYP,ibtyp_pointer,     WORK(KIVECS),buffer)
      CALL XRSXV1(-1,NBTYP,response_vector_y,eigenvectors_reduced,NFREQ,IBTYP,ibtyp_pointer,   WORK(KIVECS),buffer)
      CALL MEMREL('LRVEC1',WORK,KIVECS,KIVECS,KFREE,LFREE)
      CLOSE(LUBVEC,STATUS='DELETE')
!
!     Write response vector to file: / Z \
!                                    \ Y*/
!
      OPEN(LUFILE,FILE=FILRSP,FORM='UNFORMATTED', &
           ACCESS='DIRECT',RECL=8*max_full_length*NZ,STATUS='UNKNOWN')
!
!     Zeroing response vector
!
      CALL DZERO(response_vector_full,max_full_length*NZ)

      NZYPAR = 2*half_length

      DO IFR=1,NFREQ

         DO IZ=1,NZ
            CALL DCOPY(half_length,response_vector_z(1,IZ,IFR),1,    response_vector_full(1+(IZ-1)*NZYPAR),1)
            CALL DCOPY(half_length,response_vector_z(1,IZ,IFR),1,     response_vector_full(1+half_length+(IZ-1)*NZYPAR),1)
            CALL DAXPY(half_length,D1,response_vector_y(1,IZ,IFR),1,        response_vector_full(1+(IZ-1)*NZYPAR),1)
            CALL DAXPY(half_length,DM1,response_vector_y(1,IZ,IFR),1,     response_vector_full(1+half_length+(IZ-1)*NZYPAR),1)
            IF (IPQTOQ(IZ,JSYMOP-1).GT.1)  &
                 CALL DSCAL(half_length,DM1,     response_vector_full(1+half_length+(IZ-1)*NZYPAR),1)
         END DO
         CALL DSCAL(max_full_length*NZ,DH,response_vector_full,1)
         CALL INDQR(property_index,eigenvalues_reduced(IFR),IREC)
         CALL WRTDAC(LUFILE,max_full_length*NZ,response_vector_full,IREC)
!
!     Write info to QRINFO file to be used on restart
!
         DNORM=DNORM2(max_full_length*NZ,response_vector_full,1)
!        CALL WRTQR(PRPNAM(property_index),TYP,eigenvalues_reduced(IFR),IREC,DNORM)

      END DO
      CLOSE(LUFILE,STATUS='KEEP')
!
!     Print section
!
!     IF(IPRXQR.GE.10) THEN
         WRITE(LUPRI,'(/2A,/A,5(/5F12.8))')  &
              'blub Writing response vectors with label: ', &
              PRPNAM(property_index),' at frequencies: ', &
              (eigenvalues_reduced(I), I=1,NFREQ)
         WRITE(LUPRI,'(2A)') ' to file: ',FILRSP
!
         CALL HEADER('X(+) vectors in LRVEC1',-1)
         CALL PRBVEC(LUPRI,response_vector_z,NFREQ,half_length)
         CALL HEADER('X(-) vectors in LRVEC1',-1)
         CALL PRBVEC(LUPRI,response_vector_y,NFREQ,half_length)
!     ENDIF








    deallocate(response_vector_z)
    deallocate(response_vector_y)
    deallocate(buffer)

  end subroutine construct_response_vector

end module aoresponse_main_m
