module fde_eval_dirac

   use fde_cfg
   use fde_data
   use vc_grid
   use density_eval
   use ao_eval
   use basis_info

   
   implicit none

   public fde_get_dmat_dirac
   public fde_get_ngn_dirac
   public fde_get_vc_dirac
   public fde_free_dmat_dirac

   integer, save :: nbasis, nbasis_orb, nbasis_aux
   integer, save :: dfcoef_unit

   logical, save :: is_dirac_initialized = .false.

   contains

      subroutine dirac_data_interface 
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#ifdef MOD_DNF
#include "densfit.h"
#endif
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"

         nbasis_orb = ntbas(0) 
!         nbasis_aux = nbasisaux
!#ifdef MOD_DNF
!         if (DF_FITVXC) nbasis = MAX(nbasis_orb, nbasis_aux)
!#else
         nbasis = nbasis_orb
!#endif
         dfcoef_unit = LUCOEF

      end subroutine dirac_data_interface 

      subroutine fde_get_dmat_dirac(dmat) 
         integer :: ierr, dummy
         real(kind=8), allocatable, intent(out) :: dmat(:,:,:)
         real(kind=8), allocatable :: cmo(:,:,:)
         real(kind=8), allocatable :: work(:) 

         if (.not.is_dirac_initialized) call dirac_data_interface

         allocate(dmat(nbasis_orb,nbasis_orb,4))

         select case(fde_export_level)
         case('MP2')
            write(fde_log_unit,*) 'Reading CC_DENSITY'
            call GET_CC_DENSITY(fde_export_level,'AO',DMAT,ierr)
            if (ierr.ne.0) then
               write(fde_log_unit,*) 'Error in reading density matrix'
               call quit('error readin density matrix')
            end if
         case default
!     HF density used
            allocate(cmo(nbasis_orb,nbasis_orb,4))
            allocate(work(nbasis_orb*nbasis_orb*32))

            CMO=0.0
            work=0

            write(fde_log_unit,*) 'Reading DFCOEF'
            CALL OPNFIL(dfcoef_unit,'DFCOEF','UNKNOWN','EMBDRV')
            CALL REACMO(dfcoef_unit,'DFCOEF',CMO,dummy,dummy,dummy,3,WORK,1,1)
            deallocate(work)
            if (fde_print_level.gt.20) then
               WRITE(fde_log_unit,*) 'CMO', CMO
               WRITE(fde_log_unit,*) 'DMAT OLD', DMAT
            end if
!     SCF construct density matrix
            CALL DENMAT(DMAT,CMO,fde_print_level)
            deallocate(cmo)
         end select
!     Correct DMAT for occupationnumbers, both for scf and mp2...
         dmat=2*dmat
      end subroutine fde_get_dmat_dirac

      subroutine fde_free_dmat_dirac(dmat) 
         real(kind=8), allocatable:: dmat(:,:,:)
         deallocate(dmat)
      end subroutine fde_free_dmat_dirac



      subroutine fde_get_ngn_dirac(dmat)
         real(kind=8) :: dmat(nbasis_orb,nbasis_orb,4)
         real(kind=8) :: x,y,z, n, gn(3)
         real(kind=8), allocatable :: ao(:)
         real(kind=8), allocatable :: buffer(:)
         integer :: mat_dim
         integer :: i

         mat_dim = nbasis_orb

         allocate(ao    (total_ao_vector_length))
         allocate(buffer(total_ao_vector_length))

         call ao_eval_init(0, 0, .false.)

         do i = 1, fde_grid_np

            call get_fde_grid_coords(i,x,y,z)
            call get_ao(x, y, z, ao, buffer)
            call get_gn(n, gn, 0, mat_dim, dmat, buffer, ao, ao(ao_start(1, 0):ao_end(1, 0)))

!           call fde_data_set_active_ngn(i,n,gn(1),gn(2),gn(3))
            fde_active_n(i) = n
            fde_active_gn(1,i) = gn(1)
            fde_active_gn(2,i) = gn(2)
            fde_active_gn(3,i) = gn(3)

         end do

         deallocate(ao)
         deallocate(buffer)

      end subroutine fde_get_ngn_dirac 

!! calculate vc, density and gradient, will be stored at the fde_data module
!       call fde_calculate_active_n_gn(dmat)
!       call fde_calculate_active_elpot(dmat)
      subroutine fde_get_vc_dirac(dmat)
         real(kind=8) :: dmat(nbasis_orb,nbasis_orb,4)
         real(kind=8), allocatable :: work(:)
         real(kind=8), pointer :: gp(:,:)
         integer :: lwork, wsize
      
         wsize = 1000000
         lwork = 1

         call alloc(work,wsize)

         call set_pointer_fde_grid(gp)

         call emb_vc(dmat,             &
                     gp,               &
                     fde_active_elpot, &
                     fde_print_level,  &
                     work,lwork)

         call dealloc(work)

      end subroutine fde_get_vc_dirac

end module fde_eval_dirac
