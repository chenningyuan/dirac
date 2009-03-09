module aoresponse_input_m

  use kinds_m
  use memory_allocator

  implicit none

  integer, parameter            :: max_nr_op = 137

  character(16),           save :: op_name_to_aoproper(max_nr_op)

  logical,                 save :: aoresponse_skip_pp
  logical,                 save :: aoresponse_skip_pn

  real(kind_dp),           save :: aoresponse_threshold_response
  real(kind_dp),           save :: aoresponse_threshold_rhs
  real(kind_dp),           save :: aoresponse_freq(4)
  real(kind_dp),           save :: aoresponse_freq_process

  logical,                 save :: aoresponse_alpha
  logical,                 save :: aoresponse_beta
  logical,                 save :: aoresponse_beta_static
  logical,                 save :: aoresponse_beta_2np1
  logical,                 save :: aoresponse_gamma

  logical,                 save :: aoresponse_gamma_static
  logical,                 save :: aoresponse_dc_shg
  logical,                 save :: aoresponse_thg
  logical,                 save :: aoresponse_idri
  logical,                 save :: aoresponse_dc_kerr

  logical,                 save :: aoresponse_integrate
  logical,                 save :: aoresponse_visual
  logical,                 save :: aoresponse_test_xc
  logical,                 save :: aoresponse_test_geo
  logical,                 save :: aoresponse_test_london

contains

  function find_op_index(op_name) result(op_index)

    character(*), intent(in) :: op_name
    integer                  :: op_index, i

    op_index = -1

    do i = 1, max_nr_op
      if (op_name_to_aoproper(i) == op_name) then
        op_index = i
      end if
    end do

    if (op_index == -1) then
      call quit('find_op_index: index for '//op_name//' not found')
    end if

  end function

  subroutine add_op_by_hand(op_name,    &
                            typ,        &
                            pg_sym,     &
                            tr_sym,     &
                            prefactors, &
                            pointers,   &
                            blocks)

!   ----------------------------------------------------------------------------
    character(16), intent(in)  :: op_name
    integer,       intent(in)  :: typ
    integer,       intent(in)  :: pg_sym
    integer,       intent(in)  :: tr_sym
    real(8),       intent(in)  :: prefactors(3)
    integer,       intent(in)  :: pointers(3)
    character(4),  intent(in)  :: blocks
!   ----------------------------------------------------------------------------
    integer                    :: i
!   ----------------------------------------------------------------------------

#include "dcbxpr.h"
#include "dcbprl.h"

    i = nprps + 1

    prpnam(i)  = op_name
    iprptyp(i) = typ
    iprpsym(i) = pg_sym
    iprptim(i) = tr_sym

    facprp(1:3,  i) = prefactors(1:3)
    iprplbl(1:3, i) = pointers(1:3)

    nprps = i

    op_name_to_aoproper(i) = op_name

    i = nprplbl + 1

    prplbl(i)  = op_name
    pdoint(i)  = blocks
    iprlrep(i) = pg_sym - 1
    iprltyp(i) = tr_sym

    nprplbl = i

  end subroutine

  subroutine add_op(op_name,       &
                    typ,           &
                    nr_components, &
                    prefactors,    &
                    labels)

!   ----------------------------------------------------------------------------
    character(7), intent(in)  :: op_name
    integer,      intent(in)  :: typ
    integer,      intent(in)  :: nr_components
    real(8),      intent(in)  :: prefactors(nr_components)
    character(8), intent(in)  :: labels(nr_components)
!   ----------------------------------------------------------------------------
    integer                   :: i
    integer                   :: pg_sym
    integer                   :: tr_sym
    real(8),      allocatable :: work(:)
    integer                   :: lwork
!   ----------------------------------------------------------------------------

    lwork = 200000
    call alloc(work, lwork)

    call xprind(op_name//'         ', &                                
                typ,                  &              
                nr_components,        &                           
                prefactors,           &                                
                labels,               &                                
                i,                    &            
                pg_sym,               &                 
                tr_sym,               &                     
                work,                 &               
                lwork,                &                
                0)

    call dealloc(work)

    op_name_to_aoproper(i) = op_name

  end subroutine

end module
