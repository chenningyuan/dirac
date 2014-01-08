      module qcorr_cfg

!     variable block for Q-correction module for LUCI* CI modules. The naming convention is as follows:
!     _qcorr: input parameters
!     _q    : internal variables transferred from each CI module in the interface routine
!     written by S. Heinen and S. Knecht, ETH Zuerich winter 2013/14.

      implicit none

! logical block
      logical, public                  :: luci_cfg_qcorr = .false.
! integer block
      integer, public                  :: nref_qcorr  = -1
      integer, public                  :: print_qcorr = -1
      integer, public                  :: NGAS_q      =  0
      integer, public                  :: MXNSTR_q    =  0
      integer, public                  :: NMS2VAL_q   =  0
      integer, public                  :: NORB_q      =  0
      integer, public                  :: MXPNGAS_q   =  0
      integer, public                  :: NSMST_q     =  0
      integer, public,   allocatable   :: IGSOCC_q(:,:)
      integer, public,   allocatable   :: IBSPGPFTP_q(:)
      integer, public,   allocatable   :: NELFSPGP_q(:,:)
      integer*8, public, allocatable   :: KNSTSO_q(:)
      integer*8, public, allocatable   :: KNSTSO2_q(:)
      integer, public,   allocatable   :: MS2VAL_q(:)
      integer, public,   allocatable   :: NBLK_MS2_q(:)
      integer, public,   allocatable   :: NELEC_q(:)
      integer, public,   allocatable   :: NOCTYP_q(:)
      integer, public,   allocatable   :: IST_FOR_DT_q(:,:)

! real*8 block
      real(8), allocatable             :: ref_wavefunction_coeff_qcorr(:)
      real(8), public                  :: NACTEL_q    = 0

! character block
      character (len=200), allocatable :: ref_wavefunction_qcorr(:,:)

      end module 
