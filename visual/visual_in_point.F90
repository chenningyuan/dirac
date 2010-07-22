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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!dirac_copyright_end

module visual_in_point

  use memory_allocator
  use visual_cfg
  use codata
  use ao_eval
  use density_eval
  use basis_info
  use matrix_defop
  !use visual_london_in_point
  use visual_london

  implicit none

  public get_quantity_in_point

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



contains

      SUBROUTINE GET_QUANTITY_IN_POINT(nr_dmat, D, &
     &                                 ao, &
     &                                 PX,PY,PZ, &
     &                                 buffer, &
     &                                 QUANTITY_ADDMEUP)

      integer, intent(in) :: nr_dmat
      type(matrix)        :: D(nr_dmat)
      real(8) ::          AO(*)
      real(8) ::          buffer(*)
      real(8) ::          QUANTITY(0:3)
      real(8) ::          QUANTITY_ADDMEUP(0:3)
      real(8) ::          VECTOR_A(3)
      real(8) ::          VECTOR_B(3)
      real(8) :: factor, f, elf_chi, tau
      real(8) :: density, density_gradient(3), density_gradient_squared
      real(8) :: RSQ
      real(8) :: rxx, rxy, rxz, ryy, ryz, rzz
      real(8) :: px, py, pz, ddot, scalar
      real(8) :: rx, ry, rz, t, b, distance_to_origin
      integer :: ic, iq, IXVECTOR, ixyz, imat

      call get_ao(px, py, pz, ao, buffer)

      distance_to_origin = dsqrt(px*px + py*py + pz*pz)

      CALL DZERO(QUANTITY_ADDMEUP,4)

      DO IXVECTOR = 1, visual_cfg_NR_XVECTORS
        imat = vector_to_matrix(ixvector)

!       zero out quantity in point
        quantity = 0.0d0


!       electron density
!       ================

        if (visual_cfg_iq_xvector(ixvector) == iq_n) then
           call get_n(density,                           &
                      visual_cfg_irep_xvector(ixvector), &
                      D(imat)%nrow,                      &
                      D(imat)%elms,                      &
                      buffer,                            &
                      ao)

           quantity(0) = quantity(0) + density
        end if


!       electron density gradient
!       =========================

        if (visual_cfg_iq_xvector(ixvector) == iq_gradient) then
           call get_gn(density,                           &
                       density_gradient,                  &
                       visual_cfg_irep_xvector(ixvector), &
                       D(imat)%nrow,                      &
                       D(imat)%elms,                      &
                       buffer,                            &
                       ao,                                &
                       ao(ao_start(1, 0):ao_end(1, 0)))

           quantity(1) = quantity(1) + density_gradient(1)
           quantity(2) = quantity(2) + density_gradient(2)
           quantity(3) = quantity(3) + density_gradient(3)
        end if


!       elf: electron localization function
!       ===================================

        if (visual_cfg_iq_xvector(ixvector) == iq_elf) then
           call get_tau(tau,                               &
                        visual_cfg_irep_xvector(ixvector), &
                        D(imat)%nrow,                      &
                        D(imat)%elms,                      &
                        buffer,                            &
                        ao(ao_start(1, 0):ao_end(1, 0)))
           call get_gn(density,                           &
                       density_gradient,                  &
                       visual_cfg_irep_xvector(ixvector), &
                       D(imat)%nrow,                      &
                       D(imat)%elms,                      &
                       buffer,                            &
                       ao,                                &
                       ao(ao_start(1, 0):ao_end(1, 0)))
          
           density_gradient_squared = density_gradient(1)*density_gradient(1) &
                                    + density_gradient(2)*density_gradient(2) &
                                    + density_gradient(3)*density_gradient(3)
          
           b = density_gradient_squared/(8.0d0*density)
          
!              2.87123400018819d0 = 0.3d0*((3.0d0*pi*pi)**(2.0d0/3.0d0))
           f = 2.87123400018819d0*(density**(5.0d0/3.0d0))
          
           elf_chi = (tau - b)/f
          
           quantity(0) = quantity(0) + 1.0d0/(1.0d0 + elf_chi*elf_chi)
        end if


!       gamma5 density: \psi^\dag \gamma^5 \psi
!       =======================================

        if (visual_cfg_iq_xvector(ixvector) == iq_gamma5) then
           call get_gamma5(density,                           &
                           visual_cfg_irep_xvector(ixvector), &
                           D(imat)%nrow,                      &
                           D(imat)%elms,                      &
                           buffer,                            &
                           ao)

           quantity(0) = quantity(0) + density
        end if


!       spin density: \psi^\dag \Sigma \psi
!       ===================================

        if (visual_cfg_iq_xvector(ixvector) == iq_m) then
           call get_s(vector_a,                          &
                      visual_cfg_irep_xvector(ixvector), &
                      D(imat)%nrow,                      &
                      D(imat)%elms,                      &
                      buffer,                            &
                      ao)

           quantity(1) = quantity(1) + vector_a(1)
           quantity(2) = quantity(2) + vector_a(2)
           quantity(3) = quantity(3) + vector_a(3)
        end if


!       nr diamagnetic current density (hardcoded for B along z)
!       ========================================================

        if (visual_cfg_iq_xvector(ixvector) == iq_jdia) then
           call get_n(density,                           &
                      visual_cfg_irep_xvector(ixvector), &
                      D(imat)%nrow,                      &
                      D(imat)%elms,                      &
                      buffer,                            &
                      ao)

           rx = px - visual_cfg_gauge_origin(1)
           ry = py - visual_cfg_gauge_origin(2)
           rz = pz - visual_cfg_gauge_origin(3)
          
           vector_a(1) = -0.5d0*ry*density
           vector_a(2) =  0.5d0*rx*density
           vector_a(3) =  0.0d0

           quantity(1) = quantity(1) + vector_a(1)
           quantity(2) = quantity(2) + vector_a(2)
           quantity(3) = quantity(3) + vector_a(3)
        end if


!       current density: -ec \psi^\dag \alpha \psi
!       ==========================================

        if (visual_cfg_iq_xvector(ixvector) == iq_j) then
           call get_j(vector_a,                          &
                      visual_cfg_irep_xvector(ixvector), &
                      D(imat)%nrow,                      &
                      D(imat)%elms,                      &
                      buffer,                            &
                      ao)

           quantity(1) = quantity(1) + vector_a(1)
           quantity(2) = quantity(2) + vector_a(2)
           quantity(3) = quantity(3) + vector_a(3)
        end if




! Non-relativistic paramagnetic current density (orbital part):  j_para
        IQ = IQ_JPARA
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_JPARA(.TRUE., &
     &                   VECTOR_A, &
     &                   visual_cfg_irep_xvector(IXVECTOR), &
     &                   AO, &
     &                   ao(ao_start(1, 0):ao_end(1, 0)), &
     &                   D(imat)%elms, &
     &                   buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

! Divergence of spin density:    \nabla \cdot  \psi^\dag \Sigma \psi
        IQ = IQ_DIVM
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVM(.TRUE., &
     &                  SCALAR, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(1,1.0d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

! Divergence of current density: -ec \nabla \cdot  \psi^\dag \alpha \psi
        IQ = IQ_DIVJ
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVJ(.TRUE., &
     &                  SCALAR, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(1,1.0d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

!                  \psi^\dag (\Sigma \cdot p) \psi
        IQ = IQ_H
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVM(.TRUE., &
     &                  SCALAR, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(1,0.25d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

!                c \psi^\dag (\alpha \cdot p) \psi
        IQ = IQ_K
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_DIVJ(.TRUE., &
     &                  SCALAR, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(1,0.5d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

! Curl of spin density:     \nabla \times \psi^\dag \Sigma \psi
        IQ = IQ_ROTM
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_ROTM(.TRUE., &
     &                  VECTOR_A, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

! Curl of current density: -ec \nabla \times \psi^\dag \alpha \psi
        IQ = IQ_ROTJ
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_ROTJ(.TRUE., &
     &                  VECTOR_A, &
     &                  visual_cfg_irep_xvector(IXVECTOR), &
     &                  AO, &
     &                  ao(ao_start(1, 0):ao_end(1, 0)), &
     &                  D(imat)%elms, &
     &                  buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

! Curl of non-relativistic paramagnetic current density:    \nabla \times j_para
        IQ = IQ_ROTJPARA
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_ROTJPARA(.TRUE., &
     &                      VECTOR_A, &
     &                      visual_cfg_irep_xvector(IXVECTOR), &
     &                      AO, &
     &                      ao(ao_start(1, 0):ao_end(1, 0)), &
     &                      D(imat)%elms, &
     &                      buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF


! IQ_EDIPX       = 19                   r_x n
! IQ_EDIPY       = 20                   r_y n
! IQ_EDIPZ       = 21                   r_z n

        IF((visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPX) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPY) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPZ)) THEN

           call get_n(density,                           &
                      visual_cfg_irep_xvector(ixvector), &
                      D(imat)%nrow,                      &
                      D(imat)%elms,                      &
                      buffer,                            &
                      ao)

          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPX) &
     &      CALL DAXPY(1,1.0d0,DENSITY*PX,1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPY) &
     &      CALL DAXPY(1,1.0d0,DENSITY*PY,1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPZ) &
     &      CALL DAXPY(1,1.0d0,DENSITY*PZ,1,QUANTITY(0),1)
        ENDIF

! IQ_MDIPX       = 22             -1/2 (r_gauge \times j)_x
! IQ_MDIPY       = 23             -1/2 (r_gauge \times j)_y
! IQ_MDIPZ       = 24             -1/2 (r_gauge \times j)_z






!       responding operator is nuclear magnetic dipole
!       ==============================================

        if ((visual_cfg_iq_xvector(ixvector) == iq_nucdipx) .or. &
            (visual_cfg_iq_xvector(ixvector) == iq_nucdipy) .or. &
            (visual_cfg_iq_xvector(ixvector) == iq_nucdipz)) then

           call get_j(vector_a,                          &
                      visual_cfg_irep_xvector(ixvector), &
                      D(imat)%nrow,                      &
                      D(imat)%elms,                      &
                      buffer,                            &
                      ao)

           rx = px - cord(1, visual_cfg_ioption_xvector(ixvector))
           ry = py - cord(2, visual_cfg_ioption_xvector(ixvector))
           rz = pz - cord(3, visual_cfg_ioption_xvector(ixvector))
           
           vector_b(1) = ry*vector_a(3) - rz*vector_a(2)
           vector_b(2) = rz*vector_a(1) - rx*vector_a(3)
           vector_b(3) = rx*vector_a(2) - ry*vector_a(1)
           
           rsq = rx*rx + ry*ry + rz*rz
           f = -1.0d0/(rsq**1.5d0)
           
           if (visual_cfg_iq_xvector(ixvector) == iq_nucdipx) then
              quantity(0) = quantity(0) + f*vector_b(1)
           end if
           if (visual_cfg_iq_xvector(ixvector) == iq_nucdipy) then
              quantity(0) = quantity(0) + f*vector_b(2)
           end if
           if (visual_cfg_iq_xvector(ixvector) == iq_nucdipz) then
              quantity(0) = quantity(0) + f*vector_b(3)
           end if
        end if


!       responding operator is magnetic dipole
!       ======================================

        if ((visual_cfg_iq_xvector(ixvector) == iq_mdipx) .or. &
            (visual_cfg_iq_xvector(ixvector) == iq_mdipy) .or. &
            (visual_cfg_iq_xvector(ixvector) == iq_mdipz)) then

           call get_j(vector_a,                          &
                      visual_cfg_irep_xvector(ixvector), &
                      D(imat)%nrow,                      &
                      D(imat)%elms,                      &
                      buffer,                            &
                      ao)

           rx = px - visual_cfg_gauge_origin(1)
           ry = py - visual_cfg_gauge_origin(2)
           rz = pz - visual_cfg_gauge_origin(3)
           
           vector_b(1) = ry*vector_a(3) - rz*vector_a(2)
           vector_b(2) = rz*vector_a(1) - rx*vector_a(3)
           vector_b(3) = rx*vector_a(2) - ry*vector_a(1)
           
           f = -0.5d0
           
           if (visual_cfg_iq_xvector(ixvector) == iq_mdipx) then
              quantity(0) = quantity(0) + f*vector_b(1)
           end if
           if (visual_cfg_iq_xvector(ixvector) == iq_mdipy) then
              quantity(0) = quantity(0) + f*vector_b(2)
           end if
           if (visual_cfg_iq_xvector(ixvector) == iq_mdipz) then
              quantity(0) = quantity(0) + f*vector_b(3)
           end if
        end if


        f = visual_cfg_common_factor
        if (visual_cfg_radial) then
           f = f*(distance_to_origin**2.0d0)
        end if
        quantity_addmeup(0) = quantity_addmeup(0) + f*quantity(0)
        quantity_addmeup(1) = quantity_addmeup(1) + f*quantity(1)
        quantity_addmeup(2) = quantity_addmeup(2) + f*quantity(2)
        quantity_addmeup(3) = quantity_addmeup(3) + f*quantity(3)

      ENDDO

      END subroutine



end module
