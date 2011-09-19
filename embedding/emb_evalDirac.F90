module fde_eval_dirac

   use memory_allocator
   use fde_types
   use fde_cfg
   use fde_data
   use fde_integrate
   use vc_grid
   use density_eval
   use ao_eval
   use interface_ao
   use interface_mo_specific
   use xc_max_block_length
   use vc_grid

   implicit none

   public fde_set_nz
   public fde_get_nz

   public fde_calculate_vemb_dirac
   public fde_calculate_active_elpot
   public fde_get_active_density
   public fde_get_density_matrix
   public dirac_data_interface

   private 

   integer, save :: nba, nba_orb, nba_aux
   integer, save :: dfcoef_unit
   integer, save :: nz = 4

   logical, save :: is_dirac_initialized = .false.

   contains


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
      subroutine fde_calculate_active_elpot(dmat)
!   ----------------------------------------------------------------------------
         real(kind=8) :: dmat(:)
         real(kind=8), allocatable :: work(:)
         integer      :: work_dim=1000000
         integer      :: iprint = 0
         
         call alloc(work,work_dim, id="scratch space emb_vc")
         write (*,*) 'Generating the Coulomb potential ...'
         call emb_vc(dmat,fde_grid,fde_active_elpot,iprint,WORK,work_dim) 
         write (*,*) 'done!'
         call dealloc(work)
         
      end subroutine fde_calculate_active_elpot  
      

!   ----------------------------------------------------------------------------
      subroutine fde_calculate_vemb_dirac(mat_dim,dmat,vmat)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: mat_dim
         real(kind=8), target, intent(in) :: dmat(*)
         real(kind=8), target  :: vmat(*)

         write (*,*) 'debug, inside fde_calculate_vemb_dirac'
!        write (*,*) '       before fde_get_ngn_dirac'
!        call fde_get_ngn_dirac(dmat,mat_dim)
!        write (*,*) '       after  fde_get_ngn_dirac'
         write (*,*) '       before integrate_fde'
         call integrate_fde(fde_mat_dim   = mat_dim, &
                            fde_nz        = nz,      &
                            fde_dmat_0    = dmat,    &
                            fde_nr_dmat   = 0,       &
                            fde_nr_fmat   = 1,       &
                            fde_fmat      = vmat,    &
                            fde_do_potential = .true.)
         write (*,*) '       after  integrate_fde'

      end subroutine fde_calculate_vemb_dirac
      
      
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


!   ----------------------------------------------------------------------------
   subroutine fde_get_active_density(dmat, mat_dim, get_gradient)
!   ----------------------------------------------------------------------------
      real(kind=8), intent(in) :: dmat(*)
      integer, intent(in)      :: mat_dim
      logical, intent(in) :: get_gradient
 
      integer :: ipoint
      real(8) :: nelec, x(1), y(1), z(1), w(1)
      integer :: nr_ao
 
#ifndef NEWSCHEME    
      real(kind=8) :: n_0(max_block_length)
      real(kind=8) :: gn_0_pointwise(3)
      real(8), allocatable :: ao(:, :)
      real(8), allocatable :: buffer(:, :)
      integer :: i
      integer :: max_ao_g_order, max_ao_m_order

      max_ao_g_order = 0
      if (get_gradient) then
        max_ao_g_order = 1
      end if
     
      max_ao_m_order = 0
     
      call ao_eval_init(max_ao_g_order, &
                        max_ao_m_order, &
                        .true.)
     
                         
      nr_ao = nr_ao_slices*nr_ao_cartesian
      
      call alloc(ao,     1, nr_ao, id="ao")
      call alloc(buffer, 1, nr_ao, id="buffer")
#else
      INTEGER :: NDER=0, DOLND=0
      integer :: pl
      real(kind=8), allocatable :: GAO(:), GAO1(:),GAO2(:)
      real(kind=8), allocatable :: GAB1(:),BUF(:),NCNT(:)
      
      nr_ao = nba
      write (*,*) 'debug, nr_ao is',nr_ao
      
      call fde_get_print_level(pl)
      call alloc(gao ,   nr_ao, id="gao")
      call alloc(gao1, 3*nr_ao, id="gao1")
      call alloc(gao2, 6*nr_ao, id="gao2")
      call alloc(gab1, 3*nr_ao, id="gab1")
      call alloc(ncnt,   nr_ao, id="ncnt")
      call alloc(buf , 4*nr_ao, id="buf")
      
#endif
      nelec = 0.0
      do ipoint = 1, fde_grid_np

         x(1) = fde_grid(ipoint,1)
         y(1) = fde_grid(ipoint,2)
         z(1) = fde_grid(ipoint,3)
         w(1) = fde_grid(ipoint,4)

#ifndef NEWSCHEME         
         ao     = 0.0d0
         buffer = 0.0d0
!         time_ao_start = second()
         call get_ao(1,    &
                     x(1), &
                     y(1), &
                     z(1), &
                     ao,   &
                     buffer)
                     
!         time_ao = time_ao + second() - time_ao_start
!         time_density_start = second()
         if (need_ao_order(1, 0)) then
!        density and density gradient
            call get_gn(n_0(1), gn_0_pointwise, 0, mat_dim, dmat, buffer, ao)
         else
!        density
            call get_n(n_0(1), 0, mat_dim, dmat, buffer, ao)
!           for mcscf/srdft also evaluate spin density here
!          call get_s(spin_density, 0, mat_dim, dmat_0, buffer, ao)
         end if
         
         if ( n_0(1) > fde_cfg_screening ) then
            fde_active_n(ipoint) = n_0(1)
            if (get_gradient) then
               fde_active_gn(ipoint,1) = gn_0_pointwise(1)
               fde_active_gn(ipoint,2) = gn_0_pointwise(2)
               fde_active_gn(ipoint,3) = gn_0_pointwise(3)
            endif
         else
            fde_active_n(ipoint)    = 0.0
            fde_active_gn(ipoint,1) = 0.0 
            fde_active_gn(ipoint,2) = 0.0 
            fde_active_gn(ipoint,3) = 0.0 
         endif

#else
      GAO  = 0.0
      GAO1 = 0.0
      GAO2 = 0.0
      BUF  = 0.0
      NCNT = 0.0
      GAB1 = 0.0

        CALL GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,  &
                    x(1),y(1),z(1), &
                    BUF,nr_ao,2,DOLND,pl)
        CALL GETRHO(fde_active_n(ipoint),0,GAO,DMAT,BUF)
        CALL GETGRHO(fde_active_gn(ipoint,1:3),0,GAO,GAO1,DMAT,BUF)
#endif
         write (*,*) 'i:',ipoint,' (',x(1),',',y(1),',',z(1),') = ', &
                     fde_active_n(ipoint),'  w*n =',w(1) * fde_active_n(ipoint)
         nelec = nelec + w(1) * fde_active_n(ipoint)

!         time_density = time_density + second() - time_density_start

      enddo
      write (*,*) 'debug, nelec is ',nelec

#ifndef NEWSCHEME      
      call dealloc(ao)
      call dealloc(buffer)
#else
      call dealloc(gao )
      call dealloc(gao1)
      call dealloc(gao2)
      call dealloc(gab1)
      call dealloc(ncnt)
      call dealloc(buf )

#endif      
   end subroutine fde_get_active_density
   

!   ----------------------------------------------------------------------------
      subroutine fde_get_density_matrix(dmat,mat_dim,level)
!   ----------------------------------------------------------------------------
#include "implicit.h"
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
         character(len=3), intent(in) :: level
         real(kind=8), allocatable, intent(out) :: dmat(:)
         integer, intent(out) :: mat_dim
         real(kind=8), allocatable :: cmo(:)
         integer :: ierr

         call alloc(dmat, n2bbasxq)
         call alloc(cmo, ncmotq)
         mat_dim = n2bbasxq

         dmat = 0.0d0
         cmo  = 0.0d0

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
               if (fde_print_level.gt.20) then
                  WRITE(LUPRI,*) 'CMO', CMO
                  WRITE(LUPRI,*) 'DMAT OLD', DMAT
               end if
!     SCF construct density matrix
         call denmat(dmat,cmo,fde_print_level)
         call dealloc(cmo)

         end select
!     Correct DMAT for occupationnumbers, both for scf and mp2...
         dmat = 2.0d0*dmat
         
         write (*,*) 'debug, in fde_get_dm, dim is',mat_dim
      end subroutine fde_get_density_matrix


end module fde_eval_dirac
