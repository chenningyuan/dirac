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

module visual_cfg

  implicit none

  save

  integer, parameter, public :: max_nr_xvectors       = 6
  integer, parameter, public :: max_nr_points_in_list = 100

  character(6),  public :: visual_cfg_xvector_file(max_nr_xvectors)
  character(8),  public :: visual_cfg_conmat_file(max_nr_xvectors)
  character(9),  public :: visual_cfg_connection_key(max_nr_xvectors)
  character(1),  public :: visual_cfg_lao_chcomp(max_nr_xvectors)
  character(72), public :: visual_cfg_orbital_string(2)

  logical,       public :: visual_cfg_list               = .false.
  logical,       public :: visual_cfg_line               = .false.
  logical,       public :: visual_cfg_2d                 = .false.
  logical,       public :: visual_cfg_3d                 = .false.
  logical,       public :: visual_cfg_use_orbital_string = .false.
  logical,       public :: visual_cfg_integrate          = .false.
  logical,       public :: visual_cfg_radial             = .false.
  logical,       public :: visual_cfg_force_small_ao     = .false.
  logical,       public :: visual_cfg_london             = .false.
  logical,       public :: visual_cfg_london_with_sym    = .true.

  integer,       public :: visual_cfg_ncube(0:3)                         = 80
  integer,       public :: visual_cfg_nr_points_in_list                  = 0
  integer,       public :: visual_cfg_nr_xvectors                        = 0
  integer,       public :: visual_cfg_iq_xvector(max_nr_xvectors)        = 0
  integer,       public :: visual_cfg_irep_xvector(max_nr_xvectors)      = 0
  integer,       public :: visual_cfg_irep_conmat(max_nr_xvectors)       = 0
  integer,       public :: visual_cfg_itim_conmat(max_nr_xvectors)       = 0
  integer,       public :: visual_cfg_iwhich_xvector(max_nr_xvectors)    = 1
  integer,       public :: visual_cfg_ioption_xvector(max_nr_xvectors)   = 0
  integer,       public :: visual_cfg_line_nr_steps                      = 0
  integer,       public :: visual_cfg_orientation                        = 1
  integer,       public :: visual_cfg_lao_icomp(max_nr_xvectors)         = 0
  integer,       public :: visual_cfg_connection(max_nr_xvectors)        = 1
  integer,       public :: visual_cfg_lao_contributions(max_nr_xvectors) = 7 !all contributions by default

  real(8),       public :: visual_cfg_add_3d                             = 4.0d0
  real(8),       public :: visual_cfg_common_factor                      = 1.0d0
  real(8),       public :: visual_cfg_xyz_list(max_nr_points_in_list, 3) = 0.0d0
  real(8),       public :: visual_cfg_origin_2d(3)                       = 0.0d0
  real(8),       public :: visual_cfg_gauge_origin(3)                    = 0.0d0
  real(8),       public :: visual_cfg_add_2d                             = 4.0d0
  real(8),       public :: visual_cfg_step_2d                            = 0.1d0
  real(8),       public :: visual_cfg_line_from(3)                       = 0.0d0
  real(8),       public :: visual_cfg_line_to(3)                         = 0.0d0


  integer,       public :: vector_to_matrix(max_nr_xvectors) = 0
  integer,       public :: vector_to_matrix_con(max_nr_xvectors) = 0

  integer, parameter, public :: IQ_N           =  1 !                  \psi^\dag        \psi = n
  integer, parameter, public :: IQ_M           =  2 !                  \psi^\dag \Sigma \psi
  integer, parameter, public :: IQ_DIVM        =  3 !    \nabla \cdot  \psi^\dag \Sigma \psi
  integer, parameter, public :: IQ_ROTM        =  4 !    \nabla \times \psi^\dag \Sigma \psi
  integer, parameter, public :: IQ_J           =  5 !-ec               \psi^\dag \alpha \psi = j
  integer, parameter, public :: IQ_DIVJ        =  6 !-ec \nabla \cdot  \psi^\dag \alpha \psi
  integer, parameter, public :: IQ_ROTJ        =  7 !-ec \nabla \times \psi^\dag \alpha \psi
  integer, parameter, public :: IQ_JPARA       =  8 !                  j_para
  integer, parameter, public :: IQ_DIVJPARA    =  9 !    \nabla \cdot  j_para
  integer, parameter, public :: IQ_ROTJPARA    = 10 !    \nabla \times j_para
  integer, parameter, public :: IQ_JDIA        = 11 !                  j_dia (hardcoded for B along z)
  integer, parameter, public :: IQ_ELF         = 12 !                  ELF
  integer, parameter, public :: IQ_K           = 13 !                c \psi^\dag (\alpha \cdot p) \psi
  integer, parameter, public :: IQ_H           = 14 !                  \psi^\dag (\Sigma \cdot p) \psi
  integer, parameter, public :: IQ_GAMMA5      = 15 !                  \psi^\dag  \gamma^5        \psi

! WARNING: the factor 1/c^2 is missing on purpose for compatibility with the
!          operators defined in prp/pamprp.F
  integer, parameter, public :: IQ_NUCDIPX     = 16 !      -1/(r_K^3) (r_K \times j)_x
  integer, parameter, public :: IQ_NUCDIPY     = 17 !      -1/(r_K^3) (r_K \times j)_y
  integer, parameter, public :: IQ_NUCDIPZ     = 18 !      -1/(r_K^3) (r_K \times j)_z

  integer, parameter, public :: IQ_EDIPX       = 19 !                  r_x n
  integer, parameter, public :: IQ_EDIPY       = 20 !                  r_y n
  integer, parameter, public :: IQ_EDIPZ       = 21 !                  r_z n
  integer, parameter, public :: IQ_MDIPX       = 22 !            -1/2 (r_gauge \times j)_x
  integer, parameter, public :: IQ_MDIPY       = 23 !            -1/2 (r_gauge \times j)_y
  integer, parameter, public :: IQ_MDIPZ       = 24 !            -1/2 (r_gauge \times j)_z

  integer, parameter, public :: IQ_GRADIENT    = 25 ! \nabla   n 
  integer, parameter, public :: IQ_LAPLACIAN   = 26 ! \nabla^2 n



  private

end module
