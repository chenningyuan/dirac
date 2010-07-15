module visual_in_point

  use memory_allocator
  use visual_cfg
  use codata
  use ao_eval
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
      real(8) :: factor, fermi_kin, elf_chi, GRADIENT_SQUARED, density
      real(8) :: RSQ
      real(8) :: rxx, rxy, rxz, ryy, ryz, rzz
      real(8) :: px, py, pz, ddot, scalar
      real(8) :: rx, ry, rz, TF_KIN, BOSON_KIN, distance_to_origin
      integer :: ic, iq, IXVECTOR, ixyz, imat

      call get_ao(px, py, pz, ao, buffer)

      distance_to_origin = dsqrt(px*px + py*py + pz*pz)

      CALL DZERO(QUANTITY_ADDMEUP,4)

      DO IXVECTOR = 1, visual_cfg_NR_XVECTORS
        imat = vector_to_matrix(ixvector)

!       IF(visual_cfg_XVECTOR_FILE(IXVECTOR) .EQ. 'ZORDER') &

        CALL DZERO(QUANTITY,4)
!       Probability density:           \psi^\dag        \psi = n
        IQ = IQ_N
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          DENSITY = 0.0d0

          DO IC = 1,MC
            IF(IC .EQ. 1) THEN
              IF(.NOT. visual_cfg_IOPTION_XVECTOR(IXVECTOR) .EQ. 2) THEN
                CALL GETRH1(DENSITY,1,visual_cfg_irep_xvector(IXVECTOR), &
     &                      AO,D(imat)%elms,buffer)
              ENDIF
            ENDIF
            IF(IC .EQ. 2) THEN
              IF(.NOT. visual_cfg_IOPTION_XVECTOR(IXVECTOR) .EQ. 1) THEN
                CALL GETRH1(DENSITY,2,visual_cfg_irep_xvector(IXVECTOR), &
     &                      AO,D(imat)%elms,buffer)
              ENDIF
            ENDIF
          ENDDO

          CALL DAXPY(1,1.0d0,DENSITY,1,QUANTITY(0),1)
        ENDIF
! Gradient of probability density:  \nabla n 
        IQ = IQ_GRADIENT
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GETGRHO(VECTOR_A, &
     &                 visual_cfg_irep_xvector(IXVECTOR), &
     &                 AO, &
     &                 ao(ao_start(1, 0):ao_end(1, 0)), &
     &                 D(imat)%elms, &
     &                 buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF
! Laplacian of probability density: \nabla^2 n
        IQ = IQ_LAPLACIAN
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
!         CALL GETLRHO(RXX, &
!    &                 RXY, &
!    &                 RXZ, &
!    &                 RYY, &
!    &                 RYZ, &
!    &                 RZZ, &
!    &                 GAO, &
!    &                 GAO1, &
!    &                 GAO2, &
!    &                 D(imat)%elms, &
!    &                 buffer)

!         CALL DAXPY(1,1.0d0,RXX + RYY + RZZ,1,QUANTITY(0),1)
        ENDIF

! ELF: Electron Localization Function
        IQ = IQ_ELF
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          FERMI_KIN = 0.0d0
          DO IXYZ = 1,3
            DO IC = 1,MC
              SCALAR = 0.0d0

              CALL GETRH1(SCALAR, &
     &                    IC, &
     &                    visual_cfg_irep_xvector(IXVECTOR), &
     &                    AO(ao_start(1, 0)+(IXYZ-1)*NTBAS(0)), &
     &                    D(imat)%elms, &
     &                    buffer)

              FERMI_KIN = FERMI_KIN + SCALAR
            ENDDO
          ENDDO

          FERMI_KIN = FERMI_KIN/2.0d0

          CALL GETRHO(DENSITY, &
     &                visual_cfg_irep_xvector(IXVECTOR), &
     &                AO, &
     &                D(imat)%elms, &
     &                buffer)

          CALL GETGRHO(VECTOR_A, &
     &                 visual_cfg_irep_xvector(IXVECTOR), &
     &                 AO, &
     &                 ao(ao_start(1, 0):ao_end(1, 0)), &
     &                 D(imat)%elms, &
     &                 buffer)

          GRADIENT_SQUARED = DDOT(3,VECTOR_A,1,VECTOR_A,1)

          BOSON_KIN = GRADIENT_SQUARED/(8.0D0*DENSITY)

          TF_KIN = 0.3D0*((3.0d0*PI*PI)**(2.0d0/3.0d0))*(DENSITY**(5.0d0/3.0d0))

          ELF_CHI = (FERMI_KIN - BOSON_KIN)/TF_KIN

          SCALAR = 1.0d0/(1.0d0 + ELF_CHI*ELF_CHI)

          CALL DAXPY(1,1.0d0,SCALAR,1,QUANTITY(0),1)
        ENDIF

! j_dia (hardcoded for B along z)
        IQ = IQ_JDIA
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GETRHO(DENSITY, &
     &                visual_cfg_irep_xvector(IXVECTOR), &
     &                AO, &
     &                D(imat)%elms, &
     &                buffer)

          RX = PX - visual_cfg_GAUGE_ORIGIN(1)
          RY = PY - visual_cfg_GAUGE_ORIGIN(2)
          RZ = PZ - visual_cfg_GAUGE_ORIGIN(3)
          VECTOR_A(1) = -0.5d0*RY*DENSITY
          VECTOR_A(2) =  0.5d0*RX*DENSITY
          VECTOR_A(3) =  0.0d0

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

! Spin density: \psi^\dag \Sigma \psi
        IQ = IQ_M
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_M(.TRUE., &
     &               VECTOR_A, &
     &               visual_cfg_irep_xvector(IXVECTOR), &
     &               AO, &
     &               D(imat)%elms, &
     &               buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

! Current density: -ec               \psi^\dag \alpha \psi = j
        IQ = IQ_J
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_J(.TRUE., &
     &               VECTOR_A, &
     &               visual_cfg_irep_xvector(IXVECTOR), &
     &               AO, &
     &               D(imat)%elms, &
     &               buffer)

          CALL DAXPY(3,1.0d0,VECTOR_A,1,QUANTITY(1),1)
        ENDIF

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

!                  \psi^\dag  \gamma^5        \psi
        IQ = IQ_GAMMA5
        IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ) THEN
          CALL GET_GAMMA5(.TRUE., &
     &                    SCALAR, &
     &                    visual_cfg_irep_xvector(IXVECTOR), &
     &                    AO, &
     &                    D(imat)%elms, &
     &                    buffer)

          CALL DAXPY(1,1.0d0,SCALAR,1,QUANTITY(0),1)
        ENDIF
! IQ_NUCDIPX     = 16       -1/(r_K^3) (r_K \times j)_x
! IQ_NUCDIPY     = 17       -1/(r_K^3) (r_K \times j)_y
! IQ_NUCDIPZ     = 18       -1/(r_K^3) (r_K \times j)_z

        IF((visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPX) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPY) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPZ)) THEN

          CALL GET_J(.TRUE., &
     &               VECTOR_A, &
     &               visual_cfg_irep_xvector(IXVECTOR), &
     &               AO, &
     &               D(imat)%elms, &
     &               buffer)

          RX = PX - CORD(1,visual_cfg_IOPTION_XVECTOR(IXVECTOR))
          RY = PY - CORD(2,visual_cfg_IOPTION_XVECTOR(IXVECTOR))
          RZ = PZ - CORD(3,visual_cfg_IOPTION_XVECTOR(IXVECTOR))

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          RSQ = RX*RX + RY*RY + RZ*RZ
          FACTOR = -1.0d0/(RSQ**1.5d0)

          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPX) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(1),1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPY) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(2),1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_NUCDIPZ) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(3),1,QUANTITY(0),1)
        ENDIF

! IQ_EDIPX       = 19                   r_x n
! IQ_EDIPY       = 20                   r_y n
! IQ_EDIPZ       = 21                   r_z n

        IF((visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPX) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPY) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_EDIPZ)) THEN

          CALL GETRHO(DENSITY, &
     &                visual_cfg_irep_xvector(IXVECTOR), &
     &                AO, &
     &                D(imat)%elms, &
     &                buffer)

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

        IF((visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPX) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPY) .OR. &
     &     (visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPZ)) THEN

          CALL GET_J(.TRUE., &
     &               VECTOR_A, &
     &               visual_cfg_irep_xvector(IXVECTOR), &
     &               AO, &
     &               D(imat)%elms, &
     &               buffer)

          RX = PX - visual_cfg_GAUGE_ORIGIN(1)
          RY = PY - visual_cfg_GAUGE_ORIGIN(2)
          RZ = PZ - visual_cfg_GAUGE_ORIGIN(3)

          VECTOR_B(1) = RY*VECTOR_A(3) - RZ*VECTOR_A(2)
          VECTOR_B(2) = RZ*VECTOR_A(1) - RX*VECTOR_A(3)
          VECTOR_B(3) = RX*VECTOR_A(2) - RY*VECTOR_A(1)

          FACTOR = -0.5d0

          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPX) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(1),1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPY) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(2),1,QUANTITY(0),1)
          IF(visual_cfg_iq_xvector(IXVECTOR) .EQ. IQ_MDIPZ) &
     &      CALL DAXPY(1,FACTOR,VECTOR_B(3),1,QUANTITY(0),1)
        ENDIF

        factor = visual_cfg_common_factor
        if (visual_cfg_radial) then
          factor = factor*(distance_to_origin**2.0d0)
        end if
        call daxpy(4, factor, quantity(0), 1, quantity_addmeup(0), 1)

      ENDDO

      RETURN
      END subroutine



end module
