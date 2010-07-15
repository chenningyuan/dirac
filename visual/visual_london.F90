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
  public visual_london_alloc
  public visual_london_dealloc


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


  real(8), allocatable     :: D_r(:, :, :)
  real(8), allocatable     :: r_kl(:, :, :)
  type(matrix)             :: D_help



contains

  subroutine visual_london_init(ivector)
!   --------------------------------------------------------------------
    integer, intent(in) :: ivector

    integer :: key_l
    logical :: file_exists
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
!   gosia: i don't know if it is the best, maybe it will be problematic
!   with some compilers...
    if (btest(visual_cfg_lao_contributions(ivector), 0)) &
      write(lupri, *) 'visual_london: direct contribution included'
    if (btest(visual_cfg_lao_contributions(ivector), 1)) then
      if (visual_cfg_connection(ivector) == 0) then
        write(lupri, *) 'visual_london: (response-reorth)_ai included'
      else
        write(lupri, *) 'visual_london: response contribution included'
      end if
    end if
    if (btest(visual_cfg_lao_contributions(ivector), 2)) then
      if (visual_cfg_connection(ivector) == 0) then
        write(lupri, *) 'visual_london: overlap_ij included'
      else
        write(lupri, *) 'visual_london: reorthonorm. contrib. included'
      end if
    end if


  end subroutine visual_london_init




  subroutine visual_london_alloc()
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------
    if (.not. allocated(D_r)) call alloc(D_r, ntbas(0), ntbas(0), nz)
    call mat_init(D_help, ntbas(0), ntbas(0))
!    if (.not. visual_cfg_london_with_sym) then
      if (.not. allocated(r_kl)) call alloc(r_kl, ntbas(0), ntbas(0), 3)
!    end if
  end subroutine visual_london_alloc



  subroutine visual_london_dealloc()
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------
    if (allocated(D_r)) call dealloc(D_r)
    D_help = 0
    if (allocated(r_kl)) call dealloc(r_kl)
  end subroutine visual_london_dealloc



  subroutine create_rkl()
!   --------------------------------------------------------------------
    integer  :: mat_dim

    integer :: i, j, k, l, ixyz, ishell, isym
    real(8) :: px, py, pz
    real(8), allocatable :: p(:, :)

    logical :: debug_me

!   it is the same as xc_london_c1_init
!   --------------------------------------------------------------------
    mat_dim = ntbas(0)

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


    debug_me = .false.

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
    print_level = 0

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
!     gosia: if <none> connection then anyway read 'dS/dB' integrals

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

!   i do this to avoid 1 dcopy:
    if (visual_cfg_connection(ivector) == 1) then
      call prprea(file_unit,   &
                  p_name(1:len_trim(p_name)), &
                  rtnlbl,      &
                  p_int,       &
                  D_con_ao%elms,      &
                  idimens,     &
                  print_level)
    else
      call prprea(file_unit,   &
                  p_name(1:len_trim(p_name)), &
                  rtnlbl,      &
                  p_int,       &
                  buffer,      &
                  idimens,     &
                  print_level)
    end if

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


!   scale matrices by pfac:
!   =======================
    factor = pfac
    if (visual_cfg_connection(ivector) == 1) then
      call mat_scal(factor, D_con_ao)
    else
      call dscal(idimens, factor, buffer, 1)
    end if


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
!    else
!      call dcopy(idimens, buffer, 1, D_con_ao%elms, 1)
    end if

    call dealloc(buffer)


!   assign symmetries 
!   =================
!fixme
    !D_con_ao%ih_sym = 1
    D_con_ao%ih_sym = visual_cfg_irep_xvector(ivector) + 1
    D_con_ao%tr_sym = -1
    D_con_ao%irep   = D_con_ao%ih_sym - 1

    !call symmetry_info(p_name(1:len_trim(p_name))

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
                                       buffer,             &
                                       quantity)
!   --------------------------------------------------------------------
    type(matrix)             :: D(nr_dmat), D_con(nr_dmat_con)
    real(8)                  :: ao(*), buffer(*)
    integer, intent(in)      :: nr_dmat, nr_dmat_con
    real(8), intent(in)      :: px, py, pz
    real(8), intent(out)     :: quantity(0:3)

    real(8) :: quantity_ixvector(0:3), vector_a(3), vector_b(3)
    real(8) :: rx, ry, rz, rsq, f, factor
    integer :: ixvector

    integer :: i
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
                       ao,                  &
                       buffer,              &
                       px, py, pz,          &
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
        vector_a = 0.0d0

        call get_jblao(ixvector,            &
                       nr_dmat, D,          &
                       nr_dmat_con, D_con,  &
                       ao, buffer,          &
                       px, py, pz,          &
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
                       px, py, pz,     &
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
                       px, py, pz,          &
                       jblao)
!   --------------------------------------------------------------------
    type(matrix)             :: D(nr_dmat), D_con(nr_dmat_con)
    real(8), intent(in)      :: ao(*), buffer(*)
    integer, intent(in)      :: ixvector, nr_dmat, nr_dmat_con
    real(8), intent(in)      :: px, py, pz
    real(8), intent(out)     :: jblao(3)

    integer :: imat, imat_con
    logical :: nodirect, noreorth, norespon, debug_me

    integer :: i, j, ixyz
    real(8) :: rx, ry, rz, rsq, f

    real(8) :: vector_a(3)
!   --------------------------------------------------------------------

    D_help%elms   = 0.0d0

    debug_me = .false.

    imat     = vector_to_matrix(ixvector)
    imat_con = vector_to_matrix_con(ixvector)
    

    if (btest(visual_cfg_lao_contributions(ixvector), 0)) then

      !vector_a = 0.0d0

      if (visual_cfg_london_with_sym) then

        !call get_jblao_direct_temp(D(1), &
        !                           vector_a, & 
        !                           ao(ao_start(0,0):ao_end(0,0)), &
        !                           ao(ao_start(0,1):ao_end(0,1)), &
        !                           (/px, py, pz/), &
        !                           ixvector, &
        !                           buffer)


        call test_direct(ixvector, &
                         1, &
                         D(1)%elms, &
                         ntbas(0), &
                         ao(ao_start(0,0):ao_end(0,0)), &
                         ao(ao_start(0,1):ao_end(0,1)), &
                         (/px, py, pz/), &
                         vector_a)

      else

        call get_jblao_direct(D(1),           &
                              ao(ao_start(0, 0):ao_end(0, 0)),  &
                              (/px, py, pz/), &
                              ixvector,       &
                              buffer)
        
        f = -1.0d0
        call daxpy(ntbas(0)*ntbas(0)*nz, f, D_r, 1, D_help%elms, 1)
      end if

    end if

    if (btest(visual_cfg_lao_contributions(ixvector), 1)) then
      f = -1.0d0
      call mat_daxpy(f, D(imat), D_help)
    end if

    if (btest(visual_cfg_lao_contributions(ixvector), 2)) then
      f = -1.0d0
      call mat_daxpy(f, D_con(imat_con), D_help)
    end if

!gosia: fixme: visual_cfg_irep_xvector jest zle!!
    call get_j(.true.,                          &
               jblao,                           &
               visual_cfg_irep_xvector(ixvector), &
               !0,                               &
               ao(ao_start(0, 0):ao_end(0, 0)), &
               D_help%elms,                     &
               buffer)


    if (visual_cfg_london_with_sym) then
      f = -1.0d0
      call daxpy(3, f, vector_a, 1, jblao, 1)
    end if

  end subroutine get_jblao



  subroutine get_jblao_direct(D_0, ao, r, ixvector, buffer)
!   ---------------------------------------------------------------------------------
    type(matrix)         :: D_0
    integer, intent(in)  :: ixvector
    real(8), intent(in)  :: r(3)
    real(8)              :: ao(*), buffer(*)

    integer              :: k, l, b, c, iz, ioff
    integer              :: ixyz, i, j
    real(8)              :: f, rkl_cross_r
    logical              :: debug_me
!   ----------------------------------------------------------------------------

    debug_me = .false.

!   build 'modified' density matrix,
!   D_r = - 0.5d0*(R_mn \times r)_icomp*D_0
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

        do iz = 2, 4
          ioff = (iz - 1)*ntbas(0)*ntbas(0) + ntbas(0)*(k-1) + l
          f = 0.5d0*rkl_cross_r*D_0%elms(ioff)

!         D matrices are scaled by 2.0d0 in visual.F90
!         D_0 for iz = 2, 4 is antisymmetric
!         and now it is multiplied by antisymmetric (R_mn \times r)
!         so finally we build a symmetric matrix:
          D_r(k, l, iz) = + 0.50d0*f
          D_r(l, k, iz) = + 0.50d0*f
        end do

      end do
    end do

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

    !do while (.not. end_of_file)
    do while (.true.)
      !read(file_unit, end = 10, err = 20) t
      read(file_unit, err = 20) t
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



  subroutine get_jblao_direct_temp(D_0, vector, ao, aob1, r, ixvector, buffer)
!   ---------------------------------------------------------------------------------
    type(matrix)             :: D_0
    integer, intent(in)      :: ixvector
    real(8), intent(out)     :: vector(3)
    real(8), intent(in)      :: r(3)
    real(8)                  :: ao(ntbas(0)), aob1(ntbas(0), *), buffer(*)

    real(8) :: a, b, f, vector_a(3), vector_b(3)
    integer :: k, l, m, n, ixyz, j, ibcomp, irep, i2, i3

    logical :: debug_me
    integer :: i, isym
    real(8) :: p1, p2, p3
!   ----------------------------------------------------------------------------

    ibcomp = visual_cfg_lao_icomp(ixvector)
    irep   = visual_cfg_irep_xvector(ixvector)
    vector = 0.0d0
    vector_a = 0.0d0
    vector_b = 0.0d0

    if (ibcomp == 1) i2 = 2; i3 = 3
    if (ibcomp == 2) i2 = 3; i3 = 1
    if (ibcomp == 3) i2 = 1; i3 = 2

    a = -cval
    b = -cval
    k = 1
    l = 2
    m = 2
    n = 1
    j = 1


!   (R \times r)_i ao = (R*ao)_j * r_k - (R*ao)_k * r_j
    do ixyz = 1, 3
!     (R*ao)_j
      call get_density_sub(.true.,         &
                           .true.,         &
                           vector_a(ixyz), &
                           a, b,           &
                           k, l, m, n,     &
                           ixyz, j,        &
                           irep, &
                           !0,              &
                           aob1(1, i2),    &
                           ao,             &
                           D_0%elms,       &
                           buffer)

!     (R*ao)_k
      call get_density_sub(.true.,         &
                           .true.,         &
                           vector_b(ixyz), &
                           a, b,           &
                           k, l, m, n,     &
                           ixyz, j,        &
                           irep, &
                           !0,              &
                           aob1(1, i3),    &
                           ao,             &
                           D_0%elms,       &
                           buffer)

      !write(*, *) 'OK: ixyz, vector_a, vector_b', ixyz, vector_a(ixyz), vector_b(ixyz)
    end do

    vector(1) = vector(1) + (r(i3)*vector_a(1) - r(i2)*vector_a(1))
    vector(2) = vector(2) + (r(i3)*vector_a(2) - r(i2)*vector_a(2))
    vector(3) = vector(3) + (r(i3)*vector_a(3) - r(i2)*vector_a(3))





  end subroutine get_jblao_direct_temp







  subroutine test_direct(ivector, isym, D_0, mat_dim, ao, ao_b1, r_co, vector)
!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ivector
    integer, intent(in)    :: isym
    integer, intent(in)    :: mat_dim
    real(8), intent(in)    :: D_0(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
    real(8), intent(in)    :: ao_b1(mat_dim, *)
    real(8), intent(in)    :: r_co(3)
    real(8), intent(out)   :: vector(3)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz, jj, j
    integer                :: a, b, c, ixyz, iz, ioff
    real(8)                :: f(3)
    real(8)                :: vector_a(3)
!   ----------------------------------------------------------------------------

    vector_a = 0.0d0

!   calculate for one component of B-field at a time
    a = visual_cfg_lao_icomp(ivector)

    if (a == 1) then
      b = 2; c = 3
    else if (a == 2) then
      b = 3; c = 1
    else if (a == 3) then
      b = 1; c = 2
    end if
    
    write(*, *) 'bla: a, isym, r_co = ', a, isym, r_co(1), r_co(2), r_co(3)

!   3 component of j-vector
    do ixyz = 1, 3
      do iblock = 1, ls_nr_density_blocks(ixyz, isym - 1)
      !do iblock = 1, ls_nr_density_blocks(0, isym - 1)

        nr1 = ls_nr_ao(1, iblock, ixyz, isym - 1)
        st1 = ls_st_ao(1, iblock, ixyz, isym - 1)
        nr2 = ls_nr_ao(2, iblock, ixyz, isym - 1)
        st2 = ls_st_ao(2, iblock, ixyz, isym - 1)

        !nr1 = ls_nr_ao(1, iblock, 0, isym - 1)
        !st1 = ls_st_ao(1, iblock, 0, isym - 1)
        !nr2 = ls_nr_ao(2, iblock, 0, isym - 1)
        !st2 = ls_st_ao(2, iblock, 0, isym - 1)

        do ii = st1, st1 + nr1 - 1
          do jj = st2, st2 + nr2 - 1

            f(ixyz) = (r_co(c)*ao_b1(jj, b) - r_co(b)*ao_b1(jj, c))*ao(ii) &
                    - (r_co(c)*ao_b1(ii, b) - r_co(b)*ao_b1(ii, c))*ao(jj)

            vector_a(ixyz) = vector_a(ixyz) + f(ixyz)*D_0(ii, jj, ls_iq(ixyz))

            write(*, *) 'bla ixyz, ls_iq(ixyz), f(ixyz), vector_a(ixyz)', &
                             ixyz, ls_iq(ixyz), f(ixyz), vector_a(ixyz)
            if (ls_iq(ixyz) == 0) call quit('ls_iq(ixyz) = 0')

          end do
        end do

      end do
      vector(ixyz) = 0.5d0*cval*vector_a(ixyz)
      !write(*, *) 'bla ixyz, vector_a, vector', &
      !                 ixyz, vector_a(ixyz), vector(ixyz)
    end do

  end subroutine


end module
