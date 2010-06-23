module visual_london

  use memory_allocator
  use visual_cfg
  use ao_eval
  use basis_info
  use matrix_defop
  use matrix_genop
  use matrix_operations
  use dirac_interface


  implicit none

  public visual_london_init
  public setup_fo_conmat_ao
  public get_quantity_lao_in_point
  public create_rkl


  save

  private

#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "shells.h"
#include "orgcom.h"



  integer, parameter :: ncomp = 3
  integer, parameter :: print_level = 0



contains

  subroutine visual_london_init(ivector)
!   --------------------------------------------------------------------
    integer, intent(in) :: ivector

    integer :: key_l
    logical :: file_exists

!   initialize visual_london and print messages
!   --------------------------------------------------------------------

!   london can be done only if iq_j, iq_nucdip, iq_mdip
!   ---------------------------------------------------
    if ((visual_cfg_iq_xvector(ivector) .ne. iq_j)       .and.   &
        (visual_cfg_iq_xvector(ivector) .ne. iq_nucdipx) .and.   &
        (visual_cfg_iq_xvector(ivector) .ne. iq_nucdipy) .and.   &
        (visual_cfg_iq_xvector(ivector) .ne. iq_nucdipz) .and.   &
        (visual_cfg_iq_xvector(ivector) .ne. iq_mdipx)   .and.   &
        (visual_cfg_iq_xvector(ivector) .ne. iq_mdipy)   .and.   &
        (visual_cfg_iq_xvector(ivector) .ne. iq_mdipz)) then

      call quit('visual_london_init: london option only for iq_j, iq_nucdip and iq_mdip')

    end if


!   file_with_connection_matrices  - AOPROPER
!   -----------------------------------------
    if (visual_cfg_conmat_file(ivector) .ne. 'AOPROPER') then
      call quit('error in visual: property file has to be AOPROPER')
    end if

    inquire(file = visual_cfg_conmat_file(ivector), exist = file_exists)

    if (.not. file_exists) then
      call quit('error in visual: no property file found')
    end if

!   which_component   - x/y/z
!   -------------------------
    if (visual_cfg_lao_chcomp(ivector) == 'X') then
      visual_cfg_lao_icomp(ivector) = 1
    else if (visual_cfg_lao_chcomp(ivector) == 'Y') then
      visual_cfg_lao_icomp(ivector) = 2
    else if (visual_cfg_lao_chcomp(ivector) == 'Z') then
      visual_cfg_lao_icomp(ivector) = 3
    else
      call quit('error in visual_london: icomp != X,Y,Z')
    end if

!   which_connection  - symmetric/natural/none
!   ------------------------------------------
    key_l = len_trim(visual_cfg_connection_key(ivector))
    if (visual_cfg_connection_key(ivector)(1:key_l) == 'symmetric') then
      visual_cfg_connection(ivector) = 2
      write(lupri, *) 'visual_london: symmetric connection used'
    else if (visual_cfg_connection_key(ivector)(1:key_l) == 'natural') then
      visual_cfg_connection(ivector) = 1
      write(lupri, *) 'visual_london: natural connection used'
    else if (visual_cfg_connection_key(ivector)(1:key_l) == 'none') then
      visual_cfg_connection(ivector) = 0
      write(lupri, *) 'visual_london: <none> connection used'
    else
      write(lupri, *) 'no connection is set, the default is <natural>'
    end if

!   additional options
!   ------------------
! which later will be put as optional in the input

    if (visual_cfg_connection(ivector) == 0) then
      if (visual_cfg_none_integrals(ivector) == 'none_symmetric') then
        write(lupri, *) 'symmetric connection integrals taken'
      else if (visual_cfg_none_integrals(ivector) == 'none_natural') then
        write(lupri, *) 'natural connection integrals taken'
      else
        write(lupri, *) 'no choice'
      end if
    end if


  end subroutine visual_london_init



  subroutine create_rkl(mat_dim, r_kl)
!   --------------------------------------------------------------------
    integer, intent(in)  :: mat_dim
    real(8) :: r_kl(mat_dim, mat_dim, 3)

    integer :: i, j, k, l, ixyz, ishell, isym
    real(8) :: px, py, pz
    real(8), allocatable :: p(:, :)

    logical :: debug_me

!   it is the same as xc_london_c1_init
!   but additionally returns r_kl array
!   r_kl array contains (R_k - R_l) where k,l - centers of ao's
!   --------------------------------------------------------------------

    call alloc(p, mat_dim, 3)

    i = 1
    do ishell = 1, kmax
      do isym = 0, maxopr
        if (iand(isym, istbao(ishell)) == 0) then

          px = pt(iand(isymax(1, 1), isym))*cent(ishell, 1, 1)
          py = pt(iand(isymax(2, 1), isym))*cent(ishell, 2, 1)
          pz = pt(iand(isymax(3, 1), isym))*cent(ishell, 3, 1)

          do j = i, i + spherical_deg(ishell) - 1
            p(j, 1) = px
            p(j, 2) = py
            p(j, 3) = pz
          end do

          i = i + spherical_deg(ishell)
        end if
      end do
    end do

    do ixyz = 1, 3
      do k = 1, mat_dim
        do l = 1, mat_dim
          r_kl(k, l, ixyz) = p(k, ixyz) - p(l, ixyz)
        end do
      end do
    end do

    call dealloc(p)


!    debug_me = .true.

    if (debug_me) then

      write(lupri, *) 'create_rkl: debug r_kl'
      do ixyz = 1, 3
        do k = 1, mat_dim
          do l = 1, mat_dim
            if (r_kl(k, l, ixyz) .ne. 0.0d0) then
              write(lupri, *) 'ixyz, k, l, r_kl(k, l, ixyz) = ',&
                               ixyz, k, l, r_kl(k, l, ixyz)
            end if
          end do
        end do
      end do

    end if

  end subroutine create_rkl



  subroutine setup_fo_conmat_ao(D_con_ao, ivector)
!   --------------------------------------------------------------------
    type(matrix)         :: D_con_ao
    integer, intent(in)  :: ivector

    integer    :: icomp, file_unit, idimens, print_level, idoint(4)
    integer    :: im4, iz, ibrp, p_typ
    character  :: p_name*10, p_int*4, rtnlbl(2)*8
    logical    :: debug_me, file_exists, ldoint(4)
    real(8)    :: pfac, factor
    real(8), allocatable  :: buffer(:)

    integer :: p_index

!   -------------------------------------------------------------------

    debug_me = .false.

!   connection matrices info
!   ========================

    if (visual_cfg_connection(ivector) == 2) then

      write(lupri, *) 'symmetric connection is used in visual_london'
      p_name  = 'dS/dB'//visual_cfg_lao_chcomp(ivector)
      pfac    = -0.5d0
      p_int   = '+00+'
      idimens = nnbbasx

    else if (visual_cfg_connection(ivector) == 1) then

      write(lupri, *) 'natural connection is used in visual_london'
      p_name  = 'd|S>/dB'//visual_cfg_lao_chcomp(ivector)
      pfac    = -1.0d0
      p_int   = '+00+'
      idimens = n2bbasx

    else  if (visual_cfg_connection(ivector) == 0) then

!   visual_cfg_connection(ivector) == 0
!   gosia: if <none> connection then anyway read 'dS/dB' integrals
!   as it is easier...

      write(lupri, *) '<none> connection is used in visual_london, dS/dB integrals read'
      p_name  = 'dS/dB'//visual_cfg_lao_chcomp(ivector)
      pfac    = -0.5d0
      p_int   = '+00+'
      idimens = nnbbasx

    end if

    if(levyle .or. (zora .and. .not. zora4)) then
      p_int(4:4) = '0'
    end if

!   read connection matrix in AO basis into buffer
!   ==============================================

    file_unit = lu1int

    open(file_unit,            &
         file = visual_cfg_conmat_file(ivector), &
         status = 'OLD',       &
         form = 'UNFORMATTED', &
         access = 'SEQUENTIAL')

    call alloc(buffer, idimens)

    call prprea(file_unit,   &
                p_name(1:len_trim(p_name)), &
                rtnlbl,      &
                p_int,       &
                buffer,      &
                idimens,     &
                print_level)

    close(file_unit, status='keep')

    if (debug_me) then
!     debug after reading integrals (to compare with PRPMAO for example)
      write(lupri, *) 'setup_fo_conmat_ao: buffer after calling prprea'
      write(lupri, *) 'p_name, pfac, p_int, idimens = ', &
                       p_name, pfac, p_int, idimens
      if (visual_cfg_connection(ivector) == 1 .or.       &
          visual_cfg_connection(ivector) == 0) then
!       natural connection - square matrix:
        call output(buffer, 1, ntbas(0), 1, ntbas(0), &
                    ntbas(0), ntbas(0), 1, lupri)
      else
!       symmetric connection - lower triangle:
        call outpak(buffer, ntbas(0), 1, lupri)
      end if
    end if


!gosia fixme: for 4-component without symmetry
! to be fixed later (as in prpmsao for example)

!   scale matrices by pfac:
!   =======================
    factor = pfac
    call dscal(idimens, factor, buffer, 1)


!   create the full matrix if symmetric connection
!   ==============================================

    !if (visual_cfg_connection(ivector) == 2) then
    if ((visual_cfg_connection(ivector) == 2) .or. &
        (visual_cfg_connection(ivector) == 0)) then
      call daptge(ntbas(0), buffer, D_con_ao%elms)
    end if


!   if symmetric connection resort
!   ==============================
!    if (visual_cfg_connection(ivector) == 2) then
    if ((visual_cfg_connection(ivector) == 2) .or. &
        (visual_cfg_connection(ivector) == 0)) then
      call butobs_no_work(D_con_ao%elms, 1)
    else
      call dcopy(idimens, buffer, 1, D_con_ao%elms, 1)
    end if

    call dealloc(buffer)


!   assign symmetries 
!   =================
!(now hardcoded, later on - call prpsym? gosia fixme):
    D_con_ao%ih_sym = 1
    D_con_ao%tr_sym = -1
    D_con_ao%irep   = D_con_ao%ih_sym - 1


    if (debug_me) then
!     should be the same as connection matrix in SA-AO basis
!     which can be printed from e.g. PRPMAO when running response
      call print_mat(D_con_ao, label = 'setup_fo_conmat_ao: debug D_con_ao in SA-AO basis')
    end if


  end subroutine setup_fo_conmat_ao



  subroutine get_quantity_lao_in_point(nr_dmat, D,         &
                                       nr_dmat_con, D_con, &
                                       ao,                 &
                                       px, py, pz,         &
                                       r_kl,   &
                                       buffer,             &
                                       quantity)
!   --------------------------------------------------------------------
    type(matrix)             :: D(nr_dmat), D_con(nr_dmat_con)
    real(8)                  :: ao(*), buffer(*)
    integer, intent(in)      :: nr_dmat, nr_dmat_con
    real(8), intent(in)      :: px, py, pz
    !real(8), intent(in)      :: r_kl(ntbas(0), ntbas(0), *)
    real(8)      :: r_kl(ntbas(0), ntbas(0), *)
    real(8), intent(out)     :: quantity(0:3)

    real(8) :: quantity_ixvector(0:3), vector_a(3), vector_b(3)
    real(8) :: rx, ry, rz, rsq, f, factor
    integer :: ixvector, i, j, ixyz, imat, imat_con
    logical :: nodirect, norespon, noreorth, noreortres, nooverlap
    logical :: debug_me

!   visualize j^B and property densities based on j^B with lao:
!   1) j^B with lao:             iq_j = 5
!   2) shieldings with lao:      iq_nucdip x/y/z = 16/17/18
!   3) magnetizability with lao: iq_mdip x/y/z   = 22/23/24
!   --------------------------------------------------------------------

    call get_ao(px, py, pz, ao, buffer)

    call dzero(quantity, 4)

    do ixvector = 1, visual_cfg_nr_xvectors

      call dzero(quantity_ixvector, 4)

!     j^B with LAO
!     ============
      if(visual_cfg_iq_xvector(ixvector) .eq. iq_j) then

        vector_a = 0.0d0

        call get_jblao(ixvector,            &
                       nr_dmat, D,          &
                       nr_dmat_con, D_con,  &
                       ao, &
                       buffer,          &
                       px, py, pz, r_kl,    &
                       vector_a)

        f = 1.0d0
        call daxpy(3, f, vector_a, 1, quantity_ixvector(1), 1)

      end if


!     shielding density with LAO
!     ==========================

      if((visual_cfg_iq_xvector(ixvector) .eq. iq_nucdipx) .or. &
         (visual_cfg_iq_xvector(ixvector) .eq. iq_nucdipy) .or. &
         (visual_cfg_iq_xvector(ixvector) .eq. iq_nucdipz)) then

!       get j^B with lao to vector_a

        call get_jblao(ixvector,            &
                       nr_dmat, D,          &
                       nr_dmat_con, D_con,  &
                       ao, buffer,          &
                       px, py, pz, r_kl,    &
                       vector_a)

        rx = px - cord(1, visual_cfg_ioption_xvector(ixvector))
        ry = py - cord(2, visual_cfg_ioption_xvector(ixvector))
        rz = pz - cord(3, visual_cfg_ioption_xvector(ixvector))

        vector_b(1) = ry*vector_a(3) - rz*vector_a(2)
        vector_b(2) = rz*vector_a(1) - rx*vector_a(3)
        vector_b(3) = rx*vector_a(2) - ry*vector_a(1)

        rsq = rx*rx + ry*ry + rz*rz
        f   = -1.0d0/(rsq**1.5d0)

        if(visual_cfg_iq_xvector(ixvector) .eq. iq_nucdipx) &
          call daxpy(1, f, vector_b(1), 1, quantity_ixvector(0), 1)
        if(visual_cfg_iq_xvector(ixvector) .eq. iq_nucdipy) &
          call daxpy(1, f, vector_b(2), 1, quantity_ixvector(0), 1)
        if(visual_cfg_iq_xvector(ixvector) .eq. iq_nucdipz) &
          call daxpy(1, f, vector_b(3), 1, quantity_ixvector(0), 1)

      end if


!     magnetizability tensor with LAO
!     ===============================
!     gosia fixme it won't work now
      if((visual_cfg_iq_xvector(ixvector) .eq. iq_mdipx) .or.   &
         (visual_cfg_iq_xvector(ixvector) .eq. iq_mdipy) .or.   &
         (visual_cfg_iq_xvector(ixvector) .eq. iq_mdipz)) then

!       get j^B with lao to vector_a

        vector_a = 0.0d0

        call get_jblao(ixvector,            &
                       nr_dmat, D,          &
                       nr_dmat_con, D_con,  &
                       ao, buffer,          &
                       px, py, pz, r_kl,    &
                       vector_a)


        rx = px - visual_cfg_gauge_origin(1)
        ry = py - visual_cfg_gauge_origin(2)
        rz = pz - visual_cfg_gauge_origin(3)

        vector_b(1) = ry*vector_a(3) - rz*vector_a(2)
        vector_b(2) = rz*vector_a(1) - rx*vector_a(3)
        vector_b(3) = rx*vector_a(2) - ry*vector_a(1)

        f = -0.5d0

        if(visual_cfg_iq_xvector(ixvector) .eq. iq_mdipx)            &
          call daxpy(1, f, vector_b(1), 1, quantity_ixvector(0), 1)
        if(visual_cfg_iq_xvector(ixvector) .eq. iq_mdipy)            &
          call daxpy(1, f, vector_b(2), 1, quantity_ixvector(0), 1)
        if(visual_cfg_iq_xvector(ixvector) .eq. iq_mdipz)            &
          call daxpy(1, f, vector_b(3), 1, quantity_ixvector(0), 1)
      
      end if


!     multiply by a factor given by user and add
!     ==========================================
      factor = visual_cfg_common_factor
      call daxpy(4, factor, quantity_ixvector(0), 1, quantity(0), 1)

    end do


  end subroutine get_quantity_lao_in_point
 


  subroutine get_jblao(ixvector,            &
                       nr_dmat, D,          &
                       nr_dmat_con, D_con,  &
                       ao, buffer,          &
                       px, py, pz, &
                       r_kl,    &
                       vector_a)
!   --------------------------------------------------------------------
    type(matrix)             :: D(nr_dmat), D_con(nr_dmat_con)
    real(8), intent(in)      :: ao(*), buffer(*)
    integer, intent(in)      :: ixvector, nr_dmat, nr_dmat_con
    real(8), intent(in)      :: px, py, pz
    !real(8), intent(in)      :: r_kl(ntbas(0), ntbas(0), *)
    real(8)      :: r_kl(ntbas(0), ntbas(0), 3)
    real(8), intent(out)     :: vector_a(3)

    real(8) :: vector_b(3)
    real(8) :: factor
    integer :: icomp, imat, imat_con
    logical :: nodirect, noreorth, norespon, noreortres, nooverlap 
    logical :: debug_me

    integer :: i, j, ixyz
    real(8) :: rx, ry, rz, rsq, f
!   --------------------------------------------------------------------

    vector_a = 0.0d0

    debug_me = .false.

    imat     = vector_to_matrix(ixvector)
    imat_con = vector_to_matrix_con(ixvector)
    
!   which contributions we want to exclude (for testing)
!   ====================================================
    nodirect = .false.
    noreorth = .false.   ! <none> connection -> no overlap 
    norespon = .false.   ! <none> connection -> no (T^B-\kappa^B)


!   "direct" LAO contribution
!   =========================
    if (.not. nodirect) then
    
      vector_b = 0.0d0

!      call get_jblao_direct_temp(D(1), &
!                                 vector_b, &
!                                 ao(ao_start(0, 0):ao_end(0, 0)), &
!                                 ao(ao_start(0, 1):ao_end(0, 1)), &
!                                 (/px, py, pz/), &
!                                 ixvector, buffer)


      call get_jblao_direct(D(1),           &
                            vector_b,       &
                            ao(ao_start(0, 0):ao_end(0, 0)),  &
                            (/px, py, pz/), &
                            r_kl,           &
                            ixvector,       &
                            buffer)
    
      f = - 1.0d0
      call daxpy(3, f, vector_b, 1, vector_a, 1)

      if (debug_me) then
        write(lupri, *) 'direct contribution to j^B with lao is included'
        write(lupri, *) 'j^b with lao after direct = ', vector_b
      end if

    
    end if
    
    
!     reorthonormalization contribution
!     =================================
      if (.not. noreorth) then
      
        vector_b = 0.0d0

        call get_j(.true.,                 &
                   vector_b,               &
                   !visual_cfg_irep_conmat(ixvector), &
                   0, &
                   ao(ao_start(0, 0):ao_end(0, 0)),   &
                   D_con(imat_con)%elms,   &
                   buffer)
      
        f = -1.0d0
        call daxpy(3, f, vector_b, 1, vector_a, 1)

        if (debug_me) then
          if (visual_cfg_connection(ixvector) == 0) then
            write(lupri, *) 'contribution from overlap to j^b with lao is included'
            write(lupri, *) 'j^b with lao only overlap = ', vector_b
          else
            write(lupri, *) 'reorthonormalization contribution to j^b with lao is included'
            write(lupri, *) 'j^b with lao only reortho = ', vector_b
          end if
        end if

      end if
      
!     "response" contribution
!     =======================
      if (.not. norespon) then
      
        vector_b = 0.0d0
      
        call get_j(.true.,                 &
                   vector_b,               &
                   !visual_cfg_irep_xvector(ixvector), &
                   0, &
                   ao(ao_start(0, 0):ao_end(0, 0)),   &
                   D(imat)%elms,           &
                   buffer)
      
        f = - 1.0d0
        call daxpy(3, f, vector_b, 1, vector_a, 1)
 
        if (debug_me) then
          if (visual_cfg_connection(ixvector) == 0) then
            write(lupri, *) 'contribution from (T^B-\kappa^B) to j^b with lao is included'
            write(lupri, *) 'j^b with lao only (T^B-\kappa^B) = ', vector_b
          else
            write(lupri, *) 'response contribution to j^b with lao is included'
            write(lupri, *) 'j^b with lao only response = ', vector_b
          end if
        end if
     
      end if
    

  end subroutine get_jblao



  subroutine get_jblao_direct(D_0, vector, ao, r, &
                              r_kl, &
                              ixvector, buffer)
!   ---------------------------------------------------------------------------------
    type(matrix)             :: D_0
    integer, intent(in)      :: ixvector
    real(8), intent(out)     :: vector(3)
    real(8), intent(in)      :: r(3)
    real(8)      :: r_kl(ntbas(0), ntbas(0), *)
    real(8)                  :: ao(*), buffer(*)

    integer                  :: k, l, b, c, iz, ioff
    real(8), allocatable     :: D_r(:, :, :)
    real(8)                  :: f, rkl_cross_r

    logical :: debug_me
    integer                  :: ixyz, i, j
!   ----------------------------------------------------------------------------

    debug_me = .false.

    vector = 0.0d0

    call alloc(D_r, ntbas(0), ntbas(0), nz)
    call dzero(D_r, ntbas(0)*ntbas(0)*nz)


!   build 'modified' density matrix,
!   D_r = - 0.5d0*(R_mn \times r)_icomp D_0
!   -----------------------------------------

    if (visual_cfg_lao_icomp(ixvector) == 1) then
      b = 2; c = 3
    else if (visual_cfg_lao_icomp(ixvector) == 2) then
      b = 3; c = 1
    else if (visual_cfg_lao_icomp(ixvector) == 3) then
      b = 1; c = 2
    end if

    do k = 1, ntbas(0)
      do l = 1, k - 1

        rkl_cross_r = r_kl(k, l, b)*r(c) &
                    - r_kl(k, l, c)*r(b)


        if (debug_me) then
          write(lupri, *) 'blabla : k, l, r_kl(k, l, b), r_kl(k, l, c) = ', &
                                    k, l, r_kl(k, l, b), r_kl(k, l, c)
        end if


        do iz = 2, 4
          ioff = (iz - 1)*ntbas(0)*ntbas(0) + ntbas(0)*(k-1) + l
          f = 0.5d0*rkl_cross_r*D_0%elms(ioff)

!         D matrices are scaled by 2.0d0 in visual.F90
!         D_0 for iz = 2, 4 is antisymmetric
!         and now it is multiplied by entisymmetric (R_mn \times r)
!         so finally we build a symmetric matrix:
          D_r(k, l, iz) = + 0.50d0*f
          D_r(l, k, iz) = + 0.50d0*f
        end do

      end do
    end do

    if (debug_me) then
      write(*, *) 'bla d_0'
      call prqmat(D_0%elms, ntbas(0),ntbas(0), ntbas(0), ntbas(0), 4, (/1,2,3,4/), lupri)
      write(*, *) 'bla d_r'
      call prqmat(D_r, ntbas(0),ntbas(0), ntbas(0), ntbas(0), 4, (/1,2,3,4/), lupri)
    end if

!   multiply D_r(k, l)*(ao(k)*(-c \alpha)*ao(l))
!   --------------------------------------------

    call get_j(.true., &
               vector, &
               0,      &
               ao,     &
               D_r,    &
               buffer)

    call dealloc(D_r)


  end subroutine get_jblao_direct




  subroutine get_p_index(file_unit, prop_name, prop_index)
!   -------------------------------------------------------------------
    integer, intent(in)       :: file_unit
    character(*), intent(in)  :: prop_name
    integer, intent(out)      :: prop_index

    integer :: i, j
    character :: t(4)*8, line*72
    logical :: prop_on_file, end_of_file, in_file

!   --------------------------------------------------------------------

    rewind file_unit
    end_of_file  = .false.
    prop_on_file = .false.
    i = 0

    do while (.not. end_of_file)
      read(file_unit, end = 10, err = 20) t
      i = i + 1
      if (t(1) == '********') then
        j = j + 1
        write(lupri, *) 'in file: i, t(4) = ', i, t(4)
        if (t(4) == prop_name(1:len_trim(prop_name))) then
          prop_index = i
          prop_on_file = .true.
        end if
      end if
    end do

 10 end_of_file  = .true.
 20 write(lupri, *) 'error in reading record i = ', i


    if (.not. prop_on_file) then
      write(lupri, *) 'error in reading AOPROPER: no property ', prop_name, &
                      ' on file'
    end if


  end subroutine get_p_index



  subroutine get_jblao_direct_temp(D_0, vector, ao, aob1, r, &
                              ixvector, buffer)
!   ---------------------------------------------------------------------------------
    type(matrix)             :: D_0
    integer, intent(in)      :: ixvector
    real(8), intent(out)     :: vector(3)
    real(8), intent(in)      :: r(3)
    real(8)                  :: ao(*), aob1(ntbas(0), *), buffer(*)

    real(8) :: a, b, f, vector_a(3)
    integer :: k, l, m, n, ixyz, j, ibcomp, irep, i2, i3
    integer :: iblock, nr1, nr2, st1, st2, ii, jj
    real(8) :: f1(ntbas(0)), f2(ntbas(0))

    logical :: debug_me
!   ----------------------------------------------------------------------------

    ibcomp = visual_cfg_lao_icomp(ixvector)
    irep   = visual_cfg_irep_xvector(ixvector)
    vector = 0.0d0

    if (ibcomp == 1) i2 = 2; i3 = 3
    if (ibcomp == 2) i2 = 3; i3 = 1
    if (ibcomp == 3) i2 = 1; i3 = 2

    do jj = 1, ntbas(0)
      f1(jj) = r(i3)*aob1(jj, i2) - r(i2)*aob1(jj, i3)
    end do

    a = -cval
    b = -cval
    k = 1
    l = 2
    m = 2
    n = 1
    j = 1

!   left side: 0.5d0*r_M*ao_left*ao_right
!   =========
    vector_a = 0.0d0

    do ixyz = 1, 3
      call get_density_sub(.true., &
                           .true., &
                           vector_a(ixyz), &
                           a, b, k, l, m, n, &
                           ixyz, j, &
                           0, &
                           !aob1(1, ibcomp), &
                           f1, &
                           ao, &
                           D_0%elms, &
                           buffer)

    end do

!    vector(1) = vector(1) + (r(3)*vector_a(2) - r(2)*vector_a(3))
!    vector(2) = vector(2) + (r(1)*vector_a(3) - r(3)*vector_a(1))
!    vector(3) = vector(3) + (r(2)*vector_a(1) - r(1)*vector_a(2))

    f = 1.0d0
    call daxpy(3, f, vector_a, 1, vector, 1)

!   right side 
!   ==========
    vector_a = 0.0d0

    do ixyz = 1, 3
      call get_density_sub(.true., &
                           .true., &
                           vector_a(ixyz), &
                           a, b, k, l, m, n, &
                           ixyz, j, &
                           0, &
                           ao, &
                           f1, &
                           !aob1(1, ibcomp), &
                           D_0%elms, &
                           buffer)
    
    end do

!    vector(1) = vector(1) - (r(3)*vector_a(2) - r(2)*vector_a(3))
!    vector(2) = vector(2) - (r(1)*vector_a(3) - r(3)*vector_a(1))
!    vector(3) = vector(3) - (r(2)*vector_a(1) - r(1)*vector_a(2))


   f =  1.0d0
   call daxpy(3, f, vector_a, 1, vector, 1)



  end subroutine get_jblao_direct_temp







end module
