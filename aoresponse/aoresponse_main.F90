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
    integer :: ifreq

    real(kind_dp),    allocatable :: eigenvalues_reduced(:), eigenvectors_reduced(:), convergence(:)
    real(kind_dp),    allocatable :: prop_gradient_c(:), prop_gradient_pp(:), prop_gradient_pn(:)
    integer(kind_i4), allocatable :: ibtyp(:), ibtyp_pointer_c(:), ibtyp_pointer_pp(:), ibtyp_pointer_pn(:)


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








!   construct the response vector from the solution in the reduced space
!   ====================================================================

    call lrvec(property_index, ibtyp, ibtyp_pointer_c, ibtyp_pointer_pp, ibtyp_pointer_pn, &
               eigenvalues_reduced, eigenvectors_reduced, work(kfree), lfree)

    deallocate(ibtyp)
    deallocate(ibtyp_pointer_c)
    deallocate(ibtyp_pointer_pp)
    deallocate(ibtyp_pointer_pn)
  
    deallocate(eigenvalues_reduced)
    deallocate(eigenvectors_reduced)
    deallocate(convergence)

  end subroutine solve_perturbed_ao_dmat

end module aoresponse_main_m
