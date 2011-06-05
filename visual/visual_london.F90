
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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end


module visual_london

#ifdef MOD_UNRELEASED
  use memory_allocator
  use interface_ao
  use visual_cfg
  use ao_eval
  use density_eval
  use matrix_defop
  use dirac_interface

  implicit none

  public setup_fo_conmat_ao
  public get_jblao
  public get_jblao_ll

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

    if (visual_cfg_london_symmetric) then

      write(lupri, *) 'symmetric connection is used in visual_london'
      p_name  = 'dS/dB'//visual_cfg_london_component
      pfac    = -0.5d0
      p_int   = '+00+'
      idimens = nnbbasx

    else if (visual_cfg_london_natural) then

      write(lupri, *) 'natural connection is used in visual_london'
      p_name  = 'd|S>/dB'//visual_cfg_london_component
      pfac    = -1.0d0
      p_int   = '+00+'
      idimens = n2bbasx

    else if (visual_cfg_london_none) then
!     gosia: if <none> connection then anyway read 'dS/dB' integrals

      write(lupri, *) '<none> connection is used in visual_london, dS/dB integrals read'
      p_name  = 'dS/dB'//visual_cfg_london_component
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
         file = 'AOPROPER',    &
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
      if (visual_cfg_london_natural .or.       &
          visual_cfg_london_none) then
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
      if (visual_cfg_london_natural) then
        call qmphase('S', property_jm4pos, property_lrep, buffer, doint)
      else
        call qmphase('L', property_jm4pos, property_lrep, buffer, doint)
      end if
    end if


    if (debug_me) then
      write(lupri, *) 'phase-transformed integrals'
      if (visual_cfg_london_natural) then
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


    if (visual_cfg_london_symmetric .or. &
        visual_cfg_london_none) then

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
    D_con_ao%tr_sym = visual_cfg_itim_conmat
    D_con_ao%irep   = visual_cfg_irep_conmat


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
!      write(*, *) 'line(1:len_trim(p_name)', line(1:len_trim(p_name))

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

        visual_cfg_irep_conmat = isym - 1
        visual_cfg_itim_conmat = itim    

      end if

    end do

 99 in_file = .false.

    close(11)

    debug_me = .false.
    if (debug_me) then

      write(lupri, *) 'gosia output from read_property_info'

      write(lupri, *) 'property_name        = ', property_name
      write(lupri, *) 'property_type        = ', property_type
      write(lupri, *) 'visual_cfg_irep_conmat(ivector) = ', visual_cfg_irep_conmat
      write(lupri, *) 'isym                 = ', isym
      write(lupri, *) 'visual_cfg_itim_conmat(ivector) = ', visual_cfg_itim_conmat
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

   subroutine get_jblao(D_0, mat_dim, ao, r, j)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: D_0(mat_dim, mat_dim, *)
      real(8), intent(in)  :: ao(*)
      real(8), intent(in)  :: r(3)
      real(8), intent(out) :: j(3)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, k, l, b, c, ixyz
      real(8)              :: f
!     --------------------------------------------------------------------------

      select case (visual_cfg_london_component)
         case ('X')
            b = 2
            c = 3
         case ('Y')
            b = 3
            c = 1
         case ('Z')
            b = 1
            c = 2
      end select
    
      j = 0.0d0

      do ixyz = 1, 3
         do iblock = 1, lssl_nr_density_blocks(ixyz, 0)
            nr1 = lssl_nr_ao(1, iblock, ixyz, 0)
            st1 = lssl_st_ao(1, iblock, ixyz, 0)
            nr2 = lssl_nr_ao(2, iblock, ixyz, 0)
            st2 = lssl_st_ao(2, iblock, ixyz, 0)
            do k = st1, st1 + nr1 - 1
               do l = st2, st2 + nr2 - 1
                  f = (r(c)*ao(ao_off_g0_m1(0, b) + l) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + l))*ao(k) &
                    - (r(c)*ao(ao_off_g0_m1(0, b) + k) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + k))*ao(l)
                  
                  j(ixyz) = j(ixyz) + f*D_0(l, k, iq(ixyz))
               end do
            end do
         end do
      end do

      j = -0.5d0*cval*j

   end subroutine

   subroutine get_jblao_ll(D_0,     &
                           mat_dim, &
                           ao,      &
                           r,       &
                           j)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: D_0(mat_dim, mat_dim, *)
      real(8), intent(in)  :: ao(*)
      real(8), intent(in)  :: r(3)
      real(8), intent(out) :: j(3)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, k, l, b, c, ixyz
      real(8)              :: f
!     --------------------------------------------------------------------------

      select case (visual_cfg_london_component)
         case ('X')
            b = 2
            c = 3
         case ('Y')
            b = 3
            c = 1
         case ('Z')
            b = 1
            c = 2
      end select
    
      j = 0.0d0

      do iblock = 1, llss_nr_density_blocks(0, 0)
         nr1 = llss_nr_ao(1, iblock, 0, 0)
         st1 = llss_st_ao(1, iblock, 0, 0)
         nr2 = llss_nr_ao(2, iblock, 0, 0)
         st2 = llss_st_ao(2, iblock, 0, 0)
         do k = st1, st1 + nr1 - 1
            do l = st2, st2 + nr2 - 1
               do ixyz = 1, 3
                  f = (r(c)*ao(ao_off_g0_m1(0, b) + l) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + l))*ao(ao_off_g1_m0(ixyz, 0) + k) &
                    + (r(c)*ao(ao_off_g0_m1(0, b) + k) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + k))*ao(ao_off_g1_m0(ixyz, 0) + l) &
                    - (r(c)*ao(ao_off_g1_m1(ixyz, b) + l) &
                     - r(b)*ao(ao_off_g1_m1(ixyz, c) + l))*ao(k) &
                    - (r(c)*ao(ao_off_g1_m1(ixyz, b) + k) &
                     - r(b)*ao(ao_off_g1_m1(ixyz, c) + k))*ao(l)

                  j(ixyz) = j(ixyz) + f*D_0(l, k, iq(0))
               end do
            end do
         end do
      end do

      j = -0.5d0*j

   end subroutine

#endif /* #ifdef MOD_UNRELEASED */
end module
