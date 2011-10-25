module fde_eval_dirac

   use memory_allocator
   use fde_types
   use fde_cfg
   use fde_data
   use fde_io
   use fde_integrate
   use vc_grid
   use density_eval
   use ao_eval
   use interface_ao
   use interface_mo_specific
   use xc_max_block_length
   use vc_grid

   implicit none

   public dirac_data_interface

#ifdef MOD_UNRELEASED
   public fde_set_nz
   public fde_get_nz

   public fde_calculate_vemb_dirac
   public fde_calculate_active_elpot
   public fde_get_active_density
   public fde_get_density_matrix
   public fde_old_density_matrix_generator
#endif
   
   private 

   integer, save :: nba, nba_orb, nba_aux
   integer, save :: dfcoef_unit
   integer, save :: nz = 4

   logical, save :: is_dirac_initialized = .false.

   contains

#ifdef MOD_UNRELEASED
!   ----------------------------------------------------------------------------
      subroutine fde_set_nz(value)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: value

         nz = value
      end subroutine fde_set_nz


!   ----------------------------------------------------------------------------
      subroutine fde_get_nz(value)
!   ----------------------------------------------------------------------------
        integer, intent(out) :: value

         value = nz
      end subroutine fde_get_nz


!   ----------------------------------------------------------------------------
      subroutine fde_calculate_active_elpot(dmat,vc)
!   ----------------------------------------------------------------------------
         real(kind=8) :: dmat(:)
         real(kind=8), allocatable :: work(:)
         real(kind=8) :: vc(:)
         integer      :: work_dim=1000000
         integer      :: iprint = 0
         
         call alloc(work,work_dim, id="scratch space emb_vc")
         write (*,*) 'Generating the Coulomb potential ...'
         call emb_vc(dmat,fde_grid,vc,iprint,WORK,work_dim) 
         write (*,*) 'done!'
         call dealloc(work)
                           
      end subroutine fde_calculate_active_elpot  
      

!   ----------------------------------------------------------------------------
      subroutine fde_calculate_vemb_dirac(mat_dim,dmat,vmat)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: mat_dim
         real(kind=8), target, intent(in) :: dmat(*)
         real(kind=8), target  :: vmat(*)

         call integrate_fde(fde_mat_dim   = mat_dim, &
                            fde_nz        = nz,      &
                            fde_dmat_0    = dmat,    &
                            fde_nr_dmat   = 0,       &
                            fde_nr_fmat   = 1,       &
                            fde_fmat      = vmat,    &
                            fde_do_potential = .true.)

      end subroutine fde_calculate_vemb_dirac
#endif
      
      
!   ----------------------------------------------------------------------------
      subroutine dirac_data_interface 
!   ----------------------------------------------------------------------------
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
         nba = nba_orb
         dfcoef_unit = LUCOEF

      end subroutine dirac_data_interface 

#ifdef MOD_UNRELEASED
!   ----------------------------------------------------------------------------
   subroutine fde_get_active_density(dmat, mat_dim, get_gradient)
!   ----------------------------------------------------------------------------
      real(kind=8), intent(in) :: dmat(*)
      integer, intent(in)      :: mat_dim
      logical, intent(in) :: get_gradient
 
      call quit('write me!')
   end subroutine fde_get_active_density
   

!   ----------------------------------------------------------------------------
      subroutine fde_get_density_matrix(dmat,mat_dim,level)
!   ----------------------------------------------------------------------------
         character(len=3), intent(in) :: level
         real(kind=8), allocatable, intent(out) :: dmat(:)
         integer, intent(out) :: mat_dim
         call quit('write me!')
      end subroutine fde_get_density_matrix


!   ----------------------------------------------------------------------------
      subroutine fde_old_density_matrix_generator(dmat)
!   ----------------------------------------------------------------------------
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "mxcent.h"
#include "nuclei.h"
!     Passed variables
      type(fde_files)  :: ftmp
      type(fde_export) :: etmp
      integer :: file_unit
      character(len=60) :: file_name
      character(len=4) :: level
      real(kind=8), allocatable, intent(out) :: dmat(:)
      real(kind=8) :: test_int_density

      real(kind=8), allocatable :: cmo(:)
      REAL(kind=8), pointer :: points(:,:)
      real(8), allocatable :: grad_rho(:, :), hes_rho(:, :), vc(:), rho(:, :)

      INTEGER ::i,ierr
      INTEGER :: NDIMENSION
!     JHS THROW AWAY IN THE END?
      INTEGER :: NDER=0, DOLND=0

      REAL(KIND=8),ALLOCATABLE :: GAO(:), GAO1(:),GAO2(:),BUF(:),NCNT(:)
      REAL(KIND=8),ALLOCATABLE :: GAB1(:)

      REAL(KIND=8)::dummy,TOLS,TOLOG
      REAL(KIND=8)::gexp=0.0,ncentc=1.0,factor=1.0
      
      NDIMENSION=NTBAS(0)

      call alloc(gao ,   ndimension)
      call alloc(gao1, 3*ndimension)
      call alloc(gao2, 6*ndimension)
      call alloc(gab1, 3*ndimension)
      call alloc(ncnt,   ndimension)
      call alloc(buf , 4*ndimension)

!     Initialize the matrices to zero.

      GAO  = 0.0
      GAO1 = 0.0
      GAO2 = 0.0
      BUF  = 0.0
      NCNT = 0.0
      GAB1 = 0.0

      call alloc(dmat, n2bbasxq)
      dmat = 0.0d0
      call alloc(cmo, ncmotq)
      cmo = 0.0d0

      call fde_get_export_info(etmp)
      level = etmp%ex_level
      
      select case(level)
      case('MP2')
         write(lupri,*) 'Reading CC_DENSITY'
         call GET_CC_DENSITY(level,'AO',DMAT,ierr)
         if (ierr.ne.0) then
            write(lupri,*) 'Error in reading density matrix'
            call quit('error readin density matrix')
         end if
      case default
!     HF density used
         write(lupri,*) 'Reading DFCOEF'
         CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','EMBDRV')
         call reacmo_no_work(lucoef,    &
                             'DFCOEF',  &
                             cmo,       &
                             (/0.0d0/), &
                             (/0/),     &
                             dummy,     &
                             2)
!     SCF construct density matrix
      CALL DENMAT(DMAT,CMO,0)

!     TMP scaling
      end select
!     Correct DMAT for occupationnumbers, both for scf and mp2...
      dmat = 2.0d0*dmat

!     GET an id for LUGRIDOUT.
      call fde_get_files_info(ftmp)
      
      file_unit = ftmp%export%unit
      file_name = ftmp%export%name
      call fde_open_file(file_name,file_unit)

      WRITE(LUPRI,*) 'Reading grid.'
      nullify(points)
      CALL READ_GRID(file_unit,points)
      call fde_close_file(file_unit)
      WRITE(LUPRI,*) 'Grid read from file.'
      write(lupri,*) 'Number of gridpoints:', size(points,2)

      call alloc(rho,     1,int(SIZE(points,2)))
      call alloc(grad_rho,3,int(SIZE(points,2)))
      call alloc(hes_rho, 6,int(SIZE(points,2)))
      call alloc(vc,        int(SIZE(points,2)))

      rho      = 0.0d0
      grad_rho = 0.0d0
      hes_rho  = 0.0d0
      vc       = 0.0d0

      test_int_density = 0.0d0

      WRITE(LUPRI,*) 'Calculating the density and its derivatives.'
      DO i=1,size(points,2)
!     JHS Get the values of the basisfunctions in every gridpoint
        CALL GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,  &
                    points(1,i),points(2,i),points(3,i), &
                    BUF,NDIMENSION,2,DOLND,0)
!     JHS calculate the SCF density and its derivatives. 
        CALL GETRHO(rho(1,i),0,GAO,DMAT,BUF)
        CALL GETGRHO(grad_rho(1:3,i),0,GAO,GAO1,DMAT,BUF)
        CALL GETLRHO(hes_rho(1,i),hes_rho(2,i),hes_rho(3,i), &
             hes_rho(4,i),hes_rho(5,i),hes_rho(6,i), &
             GAO,GAO1,GAO2,DMAT,BUF)
!     JHS Calculate the coulomb potential.
         test_int_density = test_int_density + rho(1,i)*points(4,i)

         fde_active_n(i) = rho(1,i)
         fde_active_gn(i,1) = grad_rho(1,i)
         fde_active_gn(i,2) = grad_rho(2,i)
         fde_active_gn(i,3) = grad_rho(3,i)

         gf_active%n(i)    = rho(1,i)
         gf_active%gn(i,1) = grad_rho(1,i) 
         gf_active%gn(i,2) = grad_rho(2,i) 
         gf_active%gn(i,3) = grad_rho(3,i) 
      END DO

      WRITE(LUPRI,*) ' Integrated density adds up to ',test_int_density,&
     &               ' electrons.'

      call dealloc (cmo)
      call dealloc (gao)
      call dealloc (gao1)
      call dealloc (gao2)
      call dealloc (gab1)
      call dealloc (ncnt) 
      call dealloc (rho)
      call dealloc (grad_rho)
      call dealloc (hes_rho)
      call dealloc (vc)
      call deallocp (points)     
      end  subroutine fde_old_density_matrix_generator

#endif /* MOD_UNRELEASED */

end module fde_eval_dirac
