module fde_eval_dirac

   use fde_cfg
   use fde_data
   use vc_grid
   use density_eval
   use ao_eval
   use interface_ao

   
   implicit none

   public fde_get_dmat_dirac
   public fde_get_ngn_dirac
   public fde_get_vc_dirac
   public fde_free_dmat_dirac

   integer, save :: nba, nba_orb, nba_aux
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

         nba_orb = ntbas(0) 
!         nba_aux = nbasisaux
!#ifdef MOD_DNF
!         if (DF_FITVXC) nba = MAX(nba_orb, nba_aux)
!#else
         nba = nba_orb
!#endif
         dfcoef_unit = LUCOEF

      end subroutine dirac_data_interface 

      subroutine fde_get_dmat_dirac(dmat) 
         integer :: ierr, dummy
         real(kind=8), allocatable, intent(out) :: dmat(:,:,:)
         real(kind=8), allocatable :: cmo(:,:,:)
         real(kind=8), allocatable :: work(:) 

         if (.not.is_dirac_initialized) call dirac_data_interface

         allocate(dmat(nba_orb,nba_orb,4))

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
            allocate(cmo(nba_orb,nba_orb,4))
            allocate(work(nba_orb*nba_orb*32))

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
         real(kind=8) :: dmat(nba_orb,nba_orb,4)
         real(kind=8) :: n(1), gn(3)
         real(kind=8), allocatable :: ao(:,:), buffer(:,:)
         real(kind=8), allocatable :: x(:), y(:), z(:)
         integer :: mat_dim, block_length
         integer :: i

         mat_dim = nba_orb
         block_length = 1

         call alloc(ao,block_length,nr_ao_slices*nr_ao_cartesian)
         ao = 0.0d0
         call alloc(buffer,block_length,nr_ao_slices*nr_ao_cartesian)
         buffer = 0.0d0

         call alloc(x,block_length)
         call alloc(y,block_length)
         call alloc(z,block_length)

         call ao_eval_init(0, 0, .false.)

         do i = 1, fde_grid_np, block_length

            call get_fde_grid_coords(i,block_length,x,y,z)
            call get_ao(block_length,x, y, z, ao, buffer)
            call get_gn(n(1), gn, 0, mat_dim, dmat, buffer, ao)

!           call fde_data_set_active_ngn(i,n(1),gn(1),gn(2),gn(3))
            fde_active_n(i) = n(1)
            fde_active_gn(1,i) = gn(1)
            fde_active_gn(2,i) = gn(2)
            fde_active_gn(3,i) = gn(3)

         end do

         deallocate(ao)
         deallocate(buffer)
         call dealloc(ao)
         call dealloc(buffer)
         call dealloc(x)
         call dealloc(y)
         call dealloc(z)

      end subroutine fde_get_ngn_dirac 

!! calculate vc, density and gradient, will be stored at the fde_data module
!       call fde_calculate_active_n_gn(dmat)
!       call fde_calculate_active_elpot(dmat)
      subroutine fde_get_vc_dirac(dmat)
         real(kind=8) :: dmat(nba_orb,nba_orb,4)
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
