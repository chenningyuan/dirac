module openrsp_main

  use prop_test
  use prop_contribs
  use matrix_defop
  use matrix_stat
  use memory_bouncer
  use dirac_interface
  use openrsp_cfg
  use vib_prop
  use birefring

  implicit none

contains

  subroutine openrsp_driver(work, lwork)

#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    type(prop_molcfg)             :: mol
    type(matrix)  :: S, C_i, D, F
    complex(8) :: cfg_rsp_freq(10)
    complex(8) :: w(10)
    real(8) :: sum_up
    integer       :: i, j, n, k, l, number_atoms, num_coord, ierr
    integer, intent(in)    :: lwork
    real(8), intent(inout) :: work(lwork)
    complex(8), allocatable :: tsr(:)
    complex(8), allocatable :: tsr2(:)
    real(8), parameter :: zero = 0.0d0
    complex(8), parameter :: cplx_one = (1.0d0, 0.0d0)
    real(8), allocatable :: cmo_from_file(:)

!   initialize fortran 77 work memory interface
    call di_set_f77_work(work, lwork)
    
!   for debugging, remove later
!   call print_prp_common_blocks()

    cfg_rsp_freq = 0.0d0
    do i = 1, openrsp_cfg_nr_freq 
!      dirac frequencies are real, convert to complex
       cfg_rsp_freq(i) = dcmplx(openrsp_cfg_freq(i), 0.0d0)
    end do

    call get_S(S)
    call get_natoms(number_atoms)
    mol = prop_molcfg(0d0*S, number_atoms, lupri)

    allocate(cmo_from_file(n2bbasxq))
    call read_mo_coef(cmo_from_file)
    call get_C(C_i, cmo_from_file, i=1.0d0, s=0.0d0, g=1.0d0, u=1.0d0)
    deallocate(cmo_from_file)
    D = C_i*dag(C_i)

    call get_F(D, F)

!   call test_matrix_defop(S, D, F, C_i)
    C_i = 0

    write(*, *) 'electron count from S**D:', dot(S, D)

    if (openrsp_cfg_quadrupole) then
      allocate(tsr(6))
      call elec_quadrupole(mol,S, D, F, tsr)
      deallocate(tsr)
    end if

!radovan: this needs an if check
!   allocate(tsr(3))
!   call elec_dipmom(S, D, F, tsr)
!   deallocate(tsr)

    if (openrsp_cfg_alpha) then
      n = 2
      allocate(tsr(3**n))
      call elec_polariz(mol, S, D, F, cfg_rsp_freq(1), tsr)
      deallocate(tsr)
    end if

    if (openrsp_cfg_beta) then
      n = 3
      allocate(tsr(3**n))
      call elec_hypolar(mol, S, D, F, (/cfg_rsp_freq(1:n-1), -sum(cfg_rsp_freq(1:n-1))/), tsr)
      deallocate(tsr)
    end if

    if (openrsp_cfg_beta_2np1) then
      n = 3
      allocate(tsr(3**n))
      call alt_elec_hypol(mol, S, D, F, (/cfg_rsp_freq(1:n-1), -sum(cfg_rsp_freq(1:n-1))/), tsr)
      deallocate(tsr)
    end if

    if (openrsp_cfg_gamma) then
      n = 4
      allocate(tsr(3**n))
      call elec_sechyp(mol, S, D, F, (/cfg_rsp_freq(1:n-1), -sum(cfg_rsp_freq(1:n-1))/), tsr)
      deallocate(tsr)
    end if

    if (openrsp_cfg_magnetizability) then
      allocate(tsr(3*3))
      allocate(tsr2(3*3))
      call magnetiz(mol, S, D, F, cfg_rsp_freq(1), tsr, tsr2)
      deallocate(tsr)
      deallocate(tsr2)
    end if

    if (openrsp_cfg_dg_energy) then
!     call test_grcont(D%nrow, natoms_ifc(), .false., .true.)
!     call debug_1grad(D, F, natoms_ifc())
!     call debug_2grad(D,    natoms_ifc())

!     call dg_energy(S, D, F, 3*nr_atoms)
      call prop_test_gradient(mol,S, D, F, 3*nr_atoms)
    end if

        if (openrsp_cfg_roa) then
           call test_available_memory(102, D%nrow, nz)

!          print '(//a,f11.8/)','linscf calls roa_pol_Gpri_Aten_grads, freq =',cfg_rsp_freq(1)
           call get_natoms(number_atoms)
           l = 3*number_atoms
           allocate(tsr(3*3+3*3+6*3+l*3*3+l*3*3+l*6*3)) !ajt: hack
           call roa_pol_Gpri_Aten_grads(mol,S,D,F,l,cfg_rsp_freq(1)*(1d0,0d0),tsr(1:9), &
                                        tsr(10:18),tsr(19:36),tsr(37:36+9*l), &
                                        tsr(37+9*l:36+18*l),tsr(37+18*l:36+36*l))
           call VIBCTL_ifc(l,cfg_rsp_freq(1), &
                           (/((dreal(tsr(1 +i+3*j)),j=0,2),i=0,2)/), &
                           (/((dreal(tsr(10+i+3*j)),j=0,2),i=0,2)/), &
                           (/((dreal(tsr(19+i+6*j)),j=0,5),i=0,2)/), &
                           (/(0d0,i=1,3*l)/), &
                           (/(((dreal(tsr(37+     i+l*j+l*3*k)),k=0,2),j=0,2),i=0,l-1)/), &
                           (/(((dreal(tsr(37+ 9*l+i+l*j+l*3*k)),k=0,2),j=0,2),i=0,l-1)/), &
                           (/(((dreal(tsr(37+18*l+i+l*j+l*6*k)),k=0,2),j=0,5),i=0,l-1)/), &
                           work,lwork)
           deallocate(tsr)
        endif
        if (openrsp_cfg_cars) then
!          print '(//a,f11.8/)','linscf calls cars_pol_shyp_polgra, freq =',cfg_rsp_freq(1)
           call get_natoms(number_atoms)
           l = 3*number_atoms
           allocate(tsr(3+3*3+3*3*3*3+l*3*3)) !ajt: hack
           call cars_pol_shyp_polgra(mol,S,D,F,l,cfg_rsp_freq(1)*(1d0,0d0),tsr(1:3),tsr(4:12), &
                                     tsr(13:93),tsr(94:93+9*l))
           call print_shypol(dreal(cfg_rsp_freq(1))*(/-1,1,-1,1/),dreal(tsr(1:3)), &
                             (/dreal(tsr(4:12)),(0d0,i=1,9*3)/),(/(0d0,i=1,27*6)/), &
                             dreal(tsr(13:93)),lupri)
           call vib_ana_polari(mol, dreal(cfg_rsp_freq(1))*(/-1,1,-1,1/),dreal(tsr(1:3)),l,(/(0d0,i=1,3*l*4)/), &
                               (/((dreal(tsr(94+i+l*j)),j=0,8),i=0,l-1), (0d0,i=1,9*l), &
                                 ((dreal(tsr(94+i+l*j)),j=0,8),i=0,l-1), &
                                 ((dreal(tsr(94+i+l*j)),j=0,8),i=0,l-1), (0d0,i=1,9*l), &
                                 ((dreal(tsr(94+i+l*j)),j=0,8),i=0,l-1)/), &
                                (/(0d0,i=1,27*l*4)/),lupri)
           deallocate(tsr)
        endif

        ! calculates vibrational hyperpolarizability
        if (openrsp_cfg_pv_beta) then
           call get_natoms( number_atoms )
           num_coord = 3*number_atoms
           allocate( tsr(3+3*3*3+3*3*3+num_coord*3*3+num_coord*3*3*3), stat=ierr )
           if ( ierr /= 0 ) call quit( 'Failed to allocate tsr!')
           ! static case
           call vibhyp_hyp_dipgra_polgra( mol, S, D, F, num_coord, cplx_one* &
                                          (/zero,zero,zero/),              &
                                          tsr(1:3), tsr(4:30), tsr(31:57), &
                                          tsr(58:57+9*num_coord),          &
                                          tsr(58+9*num_coord:57+36*num_coord) )
           call print_shypol( (/zero,zero,zero,zero/),               &
                              dreal(tsr(1:3)),                       &
                              (/dreal(tsr(4:30)),(zero,i=1,9)/),     &
                              (/dreal(tsr(31:57)),(zero,i=1,27*5)/), &
                              (/(zero,i=1,81)/), lupri )
!          do iw = 0, size( this_info%real_freqs )-2, 2
             call vib_ana_polari( mol,   &
                                  (/zero,   &
                                    zero,   &
                                    zero,   &
                                    zero,   &
                                    zero/), &
                                  dreal(tsr(1:3)), num_coord,                                  &
!      results from \fn(vibhyp_hyp_dipgra_polgra) is (ng,3,3)
!      to \fn(vib_ana_polari) should be (3,ng,3)
                                  (/(((dreal(tsr(58+i+num_coord*j+num_coord*3*k)),             &
                                       j=0,2),i=0,num_coord-1),k=0,2),                         &
                                    (zero,i=1,num_coord*3*1)/),                                &
!      results from \fn(vibhyp_hyp_dipgra_polgra) is (ng,3,3,3)
!      to \fn(vib_ana_polari) should be (3,3,ng,3)
                                  (/(((dreal(tsr(58+num_coord*9+i+num_coord*j+num_coord*9*k)), &
                                       j=0,8),i=0,num_coord-1),k=0,2),                         &
                                    (zero,i=1,num_coord*9*3)/),                                &
                                  (/(zero,i=1,num_coord*27*4)/), lupri )
!          end do
           deallocate( tsr )
         end if

        if (openrsp_cfg_pv_gamma) then
      call get_natoms( number_atoms )
      num_coord = 3*number_atoms
      allocate( tsr(3+3*3*3*3+(3*4+3*3*6+3*3*3*4)*num_coord), stat=ierr )
      if ( ierr /= 0 ) call quit( 'Failed to allocate tsr!')
!     do iw = 0, size( this_info%real_freqs )-4, 4
           call vibshyp_shyp_dipg_polg_hypg( mol, S, D, F, num_coord, cplx_one* &
                                          (/zero,zero,zero,zero/),              &
                                          tsr(1:3), tsr(4:84),                       &
                                          tsr(85:84+12*num_coord),                   &
                                          tsr(85+12*num_coord:84+66*num_coord),      &
                                          tsr(85+66*num_coord:84+174*num_coord) )
        !ajt FIXME if freq complex, the nonzero imaginary part is never printed
           call print_shypol( (/zero,zero,zero,zero/),               &
                           dreal(tsr(1:3)),                 &
                           (/(zero,i=1,3*3*4)/),            &
                           (/(zero,i=1,3*3*3*6)/),          &
                           (/dreal(tsr(4:84))/), lupri )
        !ajt FIXME reverse index ordering in vib_ana_polari
             call vib_ana_polari( mol,   &
                                  (/zero,   &
                                    zero,   &
                                    zero,   &
                                    zero/), &
                             dreal(tsr(1:3)), num_coord,                     &
                             (/(((dreal(tsr(85+i+num_coord*(j+3*k))),        &
                                  j=0,3-1),i=0,num_coord-1),k=0,4-1)/),      &
                             (/(((dreal(tsr(85+i+num_coord*(12+j+3*3*k))),   &
                                  j=0,3*3-1),i=0,num_coord-1),k=0,6-1)/),    &
                             (/(((dreal(tsr(85+i+num_coord*(66+j+3*3*3*k))), &
                                  j=0,3*3*3-1),i=0,num_coord-1),k=0,4-1)/),  &
                             lupri )
!     end do
      deallocate( tsr )
    end if



        if (openrsp_cfg_efgb) then
           do i=1,1
!             print '(//a,f11.8/)','linscf calls efgb_Jpri_Bten_Bcal w=',cfg_rsp_freq(i)
              allocate(tsr(3+6+9+9+9+9+18+27+27+27+27+54+54)) !ajt: This tsr() thing is a hack.
                  !Response functions should perhaps eventually be stored in a dedicated module
              call efgb_Jpri_Bten_Bcal(mol,S,D,F,                           &
                      & cfg_rsp_freq(i)*(1d0,0d0)*(/1,-1,0/),                  &
                      & tsr(1:3),tsr(4:9),tsr(10:18),tsr(19:27),tsr(28:36),    &
                      & tsr(37:45),tsr(46:63),tsr(64:90),tsr(91:117),          &
                      & tsr(118:144),tsr(145:171),tsr(172:225),tsr(226:279))
              call efgb_output(dreal(cfg_rsp_freq(i)),                                 &
                      & dreal(tsr(1:3)),dreal(tsr(4:9)),dreal(tsr(10:18)),      &
                      & (/((dreal(tsr(19+i+3*j)),j=0,2),i=0,2)/),               &
                      & (/((dreal(tsr(28+i+3*j)),j=0,2),i=0,2)/),               &
                      & (/((dreal(tsr(37+i+3*j)),j=0,2),i=0,2)/),               &
                      & (/((dreal(tsr(46+i+6*j)),j=0,2),i=0,5)/),               &
                      & (/(((dreal(tsr(64 +i+3*j+9*k)),k=0,2),j=0,2),i=0,2)/),  &
                      & (/(((dreal(tsr(91 +i+3*j+9*k)),k=0,2),j=0,2),i=0,2)/),  &
                      & (/(((dreal(tsr(118+i+3*j+9*k)),k=0,2),j=0,2),i=0,2)/),  &
                      & (/(((dreal(tsr(145+i+3*j+9*k)),k=0,2),j=0,2),i=0,2)/),  &
                      & (/(((dreal(tsr(172+i+6*j+18*k)),k=0,2),j=0,2),i=0,5)/), &
                      & dreal(tsr(226:279)),lupri)
              deallocate(tsr)
           enddo
        endif

        if (openrsp_cfg_cme) then
           allocate(tsr(3+6+3*3+3*3+3*3+3*3+6*3+3*3+3*3*3+3*3*3+6*3*3 &
                        +6*3*3+3*3*3*3+3*3*3*3+6*3*3*3+6*3*3*3)) !ajt: hack
        !  do i=1,cfg_rsp_nfreq
           do i=1,1
!             print '(//a,f11.8/)','linscf calls cme_jones_eta_apri, freq =',cfg_rsp_freq(i)
              !ajt fixme Apri requires 6 extra equations and some more contractions,
              !          but is not used for Cotton-Mouton
              call cme_jones_eta_apri(mol,S,D,F,cfg_rsp_freq(i)*(1d0,0d0)*(/0,0,-1,1/), &
                             tsr(1:3),tsr(4:9), &
                             tsr(10:18),tsr(19:27),tsr(28:36),tsr(37:45),tsr(46:63),       &
                             tsr(64:72),tsr(73:99),tsr(100:126),tsr(127:180),tsr(181:234), &
                             tsr(235:315),tsr(316:396),tsr(397:558),tsr(559:720))
              call cme_output(dreal(cfg_rsp_freq(i)),dreal(tsr(64:72)), &
                              dreal(tsr(10:18)),dreal(tsr(19:27)), &
                              (/((((dreal(tsr(235+i+3*j +9*k+27*l)),l=0,2),k=0,2),j=0,2),i=0,2)/), &
                              (/((((dreal(tsr(316+i+3*j +9*k+27*l)),l=0,2),k=0,2),j=0,2),i=0,2)/), &
                              lupri)
           enddo
           deallocate(tsr)
        endif

        if (openrsp_cfg_jones) then
           allocate(tsr(3+6+3*3+3*3+3*3+3*3+6*3+3*3+3*3*3+3*3*3+6*3*3+6*3*3 &
                        +3*3*3*3+3*3*3*3+6*3*3*3+6*3*3*3))
        !  do i=1,cfg_rsp_nfreq
           do i=1,1
!             print '(//a,4f11.8/)','linscf calls cme_jones_eta_api, freq =',cfg_rsp_freq(i)
              call cme_jones_eta_apri(mol,S,D,F,cfg_rsp_freq(i)*(1d0,0d0)*(/1,0,-1,0/), &
                             tsr(1:3),tsr(4:9),tsr(10:18),tsr(19:27),tsr(28:36),    &
                             tsr(37:45),tsr(46:63),tsr(64:72),tsr(73:99),           &
                             tsr(100:126),tsr(127:180),tsr(181:234),                &
                             tsr(235:315),tsr(316:396),tsr(397:558),tsr(559:720))
              !no-London output (.false.)
              call jones_output(cfg_rsp_freq(i)*(1d0,0d0),.false.,tsr(1:3), &
                                (/(((tsr(73 +i+3*j +9*k),k=0,2),j=0,2),i=0,2)/), &
                                (/(((tsr(127+i+6*j+18*k),k=0,2),j=0,2),i=0,5)/), &
                                (/((((tsr(235+i+3*j +9*k+27*l),l=0,2),k=0,2),j=0,2),i=0,2)/), &
                                (/((((tsr(397+i+6*j+18*k+54*l),l=0,2),k=0,2),j=0,2),i=0,5)/), &
                                lupri)
              !London output (.true.)
              call jones_output(cfg_rsp_freq(i)*(1d0,0d0),.true.,tsr(1:3), &
                                (/(((tsr(100+i+3*j +9*k),k=0,2),j=0,2),i=0,2)/), &
                                (/(((tsr(181+i+6*j+18*k),k=0,2),j=0,2),i=0,5)/), &
                                (/((((tsr(316+i+3*j +9*k+27*l),l=0,2),k=0,2),j=0,2),i=0,2)/), &
                                (/((((tsr(559+i+6*j+18*k+54*l),l=0,2),k=0,2),j=0,2),i=0,5)/), &
                                lupri)
           enddo
           deallocate(tsr)
        endif



    if (openrsp_cfg_gamma_lin_dc_kerr) then

      w    = 0.0d0
      w(1) = openrsp_cfg_freq_process
      w(2) = 0.0d0
      w(3) = 0.0d0

      sum_up =          8.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 1, 1)
      sum_up = sum_up + 8.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 3, 3, 1)
      sum_up = sum_up + 2.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 3, 3)
      sum_up = sum_up + 2.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 3, 3, 1, 1)
      sum_up = sum_up + 3.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 3, 3, 3, 3)
      sum_up = sum_up/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged gamma dc Kerr',       &
            real(-sum(w)),                    &
            real(w(1)),                       &
            real(w(2)),                       &
            real(w(3)),                       &
            sum_up

    end if

    if (openrsp_cfg_gamma_lin_dc_shg) then

      w    = 0.0d0
      w(1) = openrsp_cfg_freq_process
      w(2) = openrsp_cfg_freq_process
      w(3) = 0.0d0

      sum_up =          8.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 1, 1)
      sum_up = sum_up + 2.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 3, 3, 1)
      sum_up = sum_up + 2.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 3, 1, 1, 3)
      sum_up = sum_up + 4.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 3, 3)
      sum_up = sum_up + 4.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 3, 3, 1, 1)
      sum_up = sum_up + 3.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 3, 3, 3, 3)
      sum_up = sum_up/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged gamma dc SHG',        &
            real(-sum(w)),                    &
            real(w(1)),                       &
            real(w(2)),                       &
            real(w(3)),                       &
            sum_up

    end if

    if (openrsp_cfg_gamma_lin_thg) then

      w    = 0.0d0
      w(1) = openrsp_cfg_freq_process
      w(2) = openrsp_cfg_freq_process
      w(3) = openrsp_cfg_freq_process

      sum_up =          8.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 1, 1)
      sum_up = sum_up + 6.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 3, 3)
      sum_up = sum_up + 6.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 3, 3, 1, 1)
      sum_up = sum_up + 3.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 3, 3, 3, 3)
      sum_up = sum_up/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged gamma THG',           &
            real(-sum(w)),                    &
            real(w(1)),                       &
            real(w(2)),                       &
            real(w(3)),                       &
            sum_up

    end if

    if (openrsp_cfg_gamma_lin_idri) then

      w    =  0.0d0
      w(1) =  openrsp_cfg_freq_process
      w(2) = -openrsp_cfg_freq_process
      w(3) =  openrsp_cfg_freq_process

      sum_up =          8.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 1, 1)
      sum_up = sum_up + 8.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 3, 3)
      sum_up = sum_up + 4.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 3, 1, 3)
      sum_up = sum_up + 3.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 3, 3, 3, 3)
      sum_up = sum_up/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged gamma IDRI',          &
            real(-sum(w)),                    &
            real(w(1)),                       &
            real(w(2)),                       &
            real(w(3)),                       &
            sum_up

    end if

    if (openrsp_cfg_gamma_lin_static) then

      w = 0.0d0

      sum_up =           8.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 1, 1)
      sum_up = sum_up + 12.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 1, 1, 3, 3)
      sum_up = sum_up +  3.0d0*gamma_element(mol,S, D, F, (/w(1:3)/), 3, 3, 3, 3)
      sum_up = sum_up/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged static gamma',        &
            real(-sum(w)),                    &
            real(w(1)),                       &
            real(w(2)),                       &
            real(w(3)),                       &
            sum_up

    end if

    if (openrsp_cfg_test_delta) then

!     this is 4hg (fourth harmonic generation)
      w    =  0.0d0
      w(1) =  openrsp_cfg_freq_process
      w(2) =  openrsp_cfg_freq_process
      w(3) =  openrsp_cfg_freq_process
      w(4) =  openrsp_cfg_freq_process

!     for the moment do not sum up
!     until i have the final isotropic average formula
!     just print tensor elements
      sum_up = delta_element(mol,S, D, F, (/w(1:4)/), 3, 1, 1, 1, 1)
      sum_up = delta_element(mol,S, D, F, (/w(1:4)/), 3, 3, 1, 1, 1)
      sum_up = delta_element(mol,S, D, F, (/w(1:4)/), 3, 3, 3, 1, 1)
      sum_up = delta_element(mol,S, D, F, (/w(1:4)/), 3, 3, 3, 3, 1)
      sum_up = delta_element(mol,S, D, F, (/w(1:4)/), 3, 3, 3, 3, 3)

    end if

#ifdef MERGE_WITH_LINSCA
    if (openrsp_cfg_zzz) call openrsp_zzz(openrsp_cfg_zzz_order, S, D, F)
#endif /* ifdef MERGE_WITH_LINSCA */

    S = 0
    D = 0
    F = 0

  end subroutine

end module
