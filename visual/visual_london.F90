!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://wiki.chem.vu.nl/dirac
!dirac_copyright_end


module visual_london

  use memory_allocator
  use visual_cfg
  use ao_eval
  use density_eval
  use basis_info
  use matrix_defop
  use dirac_interface


  implicit none

  public visual_london_init
  public setup_fo_conmat_ao
  public get_quantity_lao_in_point
!  public visual_london_alloc
!  public visual_london_dealloc


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

#include "dcbprl.h"


! information about property integrals read from AOPROPER file:

  character(16) ::   property_name
  integer       ::   property_type
  integer       ::   number_of_components

  real(8)       ::   property_factor(3)
  integer       ::   property_label(3)
  character(4)  ::   property_do_int(3)
  integer       ::   property_lrep(3)
  integer       ::   property_ltyp(3)
  integer       ::   property_jm4(3)
  integer       ::   property_jm4pos(3)
  integer       ::   property_iqmult(3)
  integer       ::   property_iqtopq(3)

  integer       ::   property_ihqmat(4)
  integer       ::   number_nz


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




!  subroutine visual_london_alloc()
!
!
!  end subroutine visual_london_alloc
!
!
!
!  subroutine visual_london_dealloc()
!
!
!  end subroutine visual_london_dealloc



  subroutine setup_fo_conmat_ao(D_con_ao, ivector)
!   --------------------------------------------------------------------
    type(matrix)         :: D_con_ao
    integer, intent(in)  :: ivector

    integer    :: icomp, file_unit, idimens, print_level, idoint(4)
    integer    :: im4, iz, ibrp, p_typ
    character  :: p_name*8, p_int*4, rtnlbl(2)*8
    logical    :: debug_me, file_exists, ldoint(4)
    real(8)    :: pfac, factor
    real(8), allocatable  :: buffer(:)

    integer :: p_index
    logical :: doint(4), first(4)
    integer :: i
    integer :: iq, iqp, ioff, ioff2

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

    do i = 1, 4
      if (p_int(i:i) == '0') then
        doint(i) = .false.
      else
        doint(i) = .true.
      end if
    end do

!   get info about symmetries from file (for selected p_name)
!   =========================================================
    call D_con_symmetry_info(ivector, p_name(1:len_trim(p_name)))


!   read property integrals into buffer
!   ===================================
    file_unit = lu1int

    open(file_unit,            &
         file = visual_cfg_conmat_file(ivector), &
         status = 'OLD',       &
         form = 'UNFORMATTED', &
         access = 'SEQUENTIAL')

    call alloc(buffer, idimens)

!!   i do this to avoid 1 dcopy:
!    if (visual_cfg_connection(ivector) == 1) then
!      call prprea(file_unit,   &
!                  p_name(1:len_trim(p_name)), &
!                  rtnlbl,      &
!                  p_int,       &
!                  D_con_ao%elms,      &
!                  idimens,     &
!                  print_level)
!    else
      call prprea(file_unit,   &
                  p_name(1:len_trim(p_name)), &
                  rtnlbl,      &
                  p_int,       &
                  buffer,      &
                  idimens,     &
                  print_level)
!    end if

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


!   scale integrals by pfac
!   =======================
    factor = pfac
    call dscal(idimens, factor, buffer, 1)


!   insert phases
!   =============
    if (nz .lt. 4) then
      if (visual_cfg_connection(ivector) == 1) then
        call qmphase('S', property_jm4pos, property_lrep, buffer, doint)
      else
        call qmphase('L', property_jm4pos, property_lrep, buffer, doint)
      end if
    end if


    if (debug_me) then
      write(lupri, *) 'phase-transformed integrals'
      if (visual_cfg_connection(ivector) == 1) then
        call output(buffer, 1, ntbas(0), 1, ntbas(0), ntbas(0), ntbas(0), 1,lupri)
      else
        call outpak(buffer, ntbas(0), 1, lupri)
      end if
    end if


!   gosia i'll fix it later
    do i = 1, 4 
      first(i) = .false.
    end do

    do i = 1, nz
      first(i) = .true.
    end do


    if ((visual_cfg_connection(ivector) == 2) .or. &
        (visual_cfg_connection(ivector) == 0)) then

      ioff = ntbas(0)*ntbas(0)*(property_iqtopq(1) - 1) + 1
      if (first(property_iqtopq(1))) then
        call dcopy(idimens, buffer, 1, D_con_ao%elms(ioff), 1)
        first(property_iqtopq(1)) = .false.
      end if


      do i = 1, nz
        ioff2 = ntbas(0)*ntbas(0)*(i - 1) + 1
        if (.not. first(i)) then
          call dzero(buffer, idimens)
          call dcopy(idimens, D_con_ao%elms(i), 1, buffer, 1)
          if (property_ihqmat(i) == 2) then
            call daptge(ntbas(0), buffer, D_con_ao%elms(ioff2))
          else if (property_ihqmat(i) == 1) then
            call dsptsi(ntbas(0), buffer, D_con_ao%elms(ioff2))
          end if
          call butobs_no_work(D_con_ao%elms(ioff2), 1)
        else
          call dzero(D_con_ao%elms(ioff2), n2bbasx)
        end if
      end do 

    else

      ioff = ntbas(0)*ntbas(0)*(property_iqtopq(1) - 1) + 1
      call dcopy(idimens, buffer, 1, D_con_ao%elms(ioff), 1)
      call butobs_no_work(D_con_ao%elms, nz)

    end if



    call dealloc(buffer)


!   assign symmetries 
!   =================
    D_con_ao%ih_sym = visual_cfg_irep_conmat(ivector) + 1
    D_con_ao%tr_sym = visual_cfg_itim_conmat(ivector)
    D_con_ao%irep   = visual_cfg_irep_conmat(ivector)


    if (debug_me) then
!     should be the same as connection matrix in SA-AO basis
!     which can be printed from e.g. PRPMAO when running response
      call print_mat(D_con_ao, label = 'setup_fo_conmat_ao: debug D_con_ao in SA-AO basis')
    end if


  end subroutine setup_fo_conmat_ao


  subroutine D_con_symmetry_info(ivector, p_name)
!   --------------------------------------------------------------------
    integer, intent(in)       :: ivector
    character(*), intent(in)  :: p_name

    integer       :: isym, itim, indxpr, i
    character(8)  :: line
    logical :: in_file
    logical :: debug_me
!   --------------------------------------------------------------------

    debug_me = .true.

    if (debug_me) then
      write(lupri, *) 'read symmetry info for connection matrices'
      write(lupri, *) 'entering ivector, p_name, len_trim(p_name) = ', &
                                ivector, p_name, len_trim(p_name)
    end if


    open(11, file='property_info',form='formatted',access='sequential')
    rewind(11)
    in_file = .true.

    do while (in_file)

      read(11, *, end = 99) line
      write(*, *) 'line(1:len_trim(p_name)', line(1:len_trim(p_name))

      !if (p_name == line(1:len_trim(p_name))) then
      if (line(1:2) == p_name(1:2)) then

        property_name = p_name(1:len_trim(p_name))
        read(11, *, end = 99) property_type
        read(11, *, end = 99) isym
        read(11, *, end = 99) itim

!
        read(11, *, end = 99) number_nz
        do i = 1, number_nz
          read(11, *, end = 99) property_ihqmat(i)
        end do
!
        read(11, *, end = 99) number_of_components
        do i = 1, number_of_components
          read(11, *, end = 99) property_factor(i)
          read(11, *, end = 99) property_label(i)
          read(11, *, end = 99) property_do_int(i)
          read(11, *, end = 99) property_lrep(i)
          read(11, *, end = 99) property_ltyp(i)
          read(11, *, end = 99) property_jm4(i)
          read(11, *, end = 99) property_jm4pos(i)
          read(11, *, end = 99) property_iqmult(i)
          read(11, *, end = 99) property_iqtopq(i)
        end do

        visual_cfg_irep_conmat(ivector) = isym - 1
        visual_cfg_itim_conmat(ivector) = itim    

      end if

    end do

 99 in_file = .false.

    close(11)

    debug_me = .false.
    if (debug_me) then

      write(lupri, *) 'gosia output from read_property_info'

      write(lupri, *) 'property_name        = ', property_name
      write(lupri, *) 'property_type        = ', property_type
      write(lupri, *) 'visual_cfg_irep_conmat(ivector) = ', visual_cfg_irep_conmat(ivector)
      write(lupri, *) 'isym                 = ', isym
      write(lupri, *) 'visual_cfg_itim_conmat(ivector) = ', visual_cfg_itim_conmat(ivector)
      write(lupri, *) 'itim                 = ', itim

      write(lupri, *) 'number_nz            =', number_nz
      do i = 1, number_nz
        write(lupri, *) 'i, property_ihqmat(i) = ', i, property_ihqmat(i)
      end do
!

      write(lupri, *) 'number_of_components = ', number_of_components

      do i = 1, number_of_components
       write(lupri, *) '                i   = ', i 
       write(lupri, *) 'property_factor(i)  = ', property_factor(i)
       write(lupri, *) 'property_label(i)   = ', property_label(i)
       write(lupri, *) 'property_do_int(i)  = ', property_do_int(i)
       write(lupri, *) 'property_lrep(i)    = ', property_lrep(i)
       write(lupri, *) 'property_ltyp(i)    = ', property_ltyp(i)
       write(lupri, *) 'property_jm4(i)     = ', property_jm4(i)
       write(lupri, *) 'property_jm4pos(i)  = ', property_jm4pos(i)
       write(lupri, *) 'property_iqmult(i)  = ', property_iqmult(i)
       write(lupri, *) 'property_iqtopq(i)  = ', property_iqtopq(i)
      end do

    end if


  end subroutine D_con_symmetry_info



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

    jblao = 0.0d0

    debug_me = .false.

    imat     = vector_to_matrix(ixvector)
    imat_con = vector_to_matrix_con(ixvector)

!   'direct' contribution
    if (btest(visual_cfg_lao_contributions(ixvector), 0)) then

      vector_a = 0.0d0

      call get_jblao_direct(ixvector, &
                            1, &
                            D(1)%elms, &
                            ntbas(0), &
                            ao(ao_start(0,0):ao_end(0,0)), &
                            ao(ao_start(0,1):ao_end(0,1)), &
                            (/px, py, pz/), &
                            vector_a)

      f = -1.0d0
      call daxpy(3, f, vector_a, 1, jblao, 1)

    end if

!   'response' contribution
    if (btest(visual_cfg_lao_contributions(ixvector), 1)) then

      vector_a = 0.0d0

           call get_j(vector_a,                             &
                      visual_cfg_irep_xvector(ixvector), &
                      D(imat)%nrow,                       &
                      D(imat)%elms,                       &
                      buffer,                            &
                      ao)


      f = -1.0d0
      call daxpy(3, f, vector_a, 1, jblao, 1)

    end if

!   'reorthonormalization' contribution
    if (btest(visual_cfg_lao_contributions(ixvector), 2)) then

      vector_a = 0.0d0

      call get_j(vector_a,                          &
                 visual_cfg_irep_conmat(ixvector), &
                 D_con(imat_con)%nrow,                       &
                 D_con(imat_con)%elms,                       &
                 buffer, &
                 ao(ao_start(0, 0):ao_end(0, 0)))

      f = -1.0d0
      call daxpy(3, f, vector_a, 1, jblao, 1)

    end if


  end subroutine get_jblao



  subroutine get_jblao_direct(ivector, isym, D_0, mat_dim, ao, ao_b1, r_co, vector)
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

    a = visual_cfg_lao_icomp(ivector)

    if (a == 1) then
      b = 2; c = 3
    else if (a == 2) then
      b = 3; c = 1
    else if (a == 3) then
      b = 1; c = 2
    end if
    
    do ixyz = 1, 3
      do iblock = 1, lssl_nr_density_blocks(ixyz, isym - 1)

        nr1 = lssl_nr_ao(1, iblock, ixyz, isym - 1)
        st1 = lssl_st_ao(1, iblock, ixyz, isym - 1)
        nr2 = lssl_nr_ao(2, iblock, ixyz, isym - 1)
        st2 = lssl_st_ao(2, iblock, ixyz, isym - 1)

        do ii = st1, st1 + nr1 - 1
          do jj = st2, st2 + nr2 - 1

            f(ixyz) = (r_co(c)*ao_b1(jj, b) - r_co(b)*ao_b1(jj, c))*ao(ii) &
                    - (r_co(c)*ao_b1(ii, b) - r_co(b)*ao_b1(ii, c))*ao(jj)

            vector_a(ixyz) = vector_a(ixyz) + f(ixyz)*D_0(ii, jj, iq(ixyz))

          end do
        end do

      end do
      vector(ixyz) = 0.5d0*cval*vector_a(ixyz)
    end do

  end subroutine get_jblao_direct



end module
