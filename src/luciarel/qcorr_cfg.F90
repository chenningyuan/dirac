      module qcorr_cfg

      implicit none

! integer block
      integer, public              :: NGAS_q      = 0
      integer, public              :: IBSPGPFTP_q = 0
      integer, public              :: MXNSTR_q    = 0
      integer, public              :: NMS2VAL_q   = 0
      integer, public, allocatable :: IGSOCC_q(:)
      integer, public, allocatable :: NELFSPGP_q(:,:)
      integer, public, allocatable :: KNSTSO_q(:)
      integer, public, allocatable :: KNSTSO2_q(:)
      integer, public, allocatable :: MS2VAL_q(:)
      integer, public, allocatable :: NBLK_MS2_q(:)
      integer, public, allocatable :: NELEC_q(:)
      integer, public, allocatable :: NOCTYP_q(:)
      integer, public, allocatable :: IST_FOR_DT_q(:,:)

      end module 
