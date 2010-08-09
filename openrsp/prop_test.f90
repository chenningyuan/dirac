module prop_test
!ajt/radovan: Response-related testing routines and some calculations

   use decompMod, only: DecompItem
   use matrix_defop
   use prop_contribs
   use rsp_equations
#ifdef VAR_LINSCA
   use scf_config, only: lupri
#endif
#ifdef PRG_DIRAC
   use dirac_interface
   use openrsp_output
#endif

   implicit none

   public elec_dipmom
   public elec_polariz
   public elec_hypolar
   public alt_elec_hypol
   public elec_sechyp
   public alt_elec_sechyp
   public alt2_elec_sechyp
   public magnetiz
   public elec_quadrupole
   public prop_test_gradient
#ifdef PRG_DIRAC
!  this should also work with Dalton
!  but for this need to adapt output
   public element_ooff
   public gamma_element
   public delta_element
#endif

   private

   !physical constants for conversion to non-atomic units
   real(8), parameter:: cm1 = 1/219474.631371d0, & !1 centimeter-to-minus-one in au
                        cvl = 137.03599907d0,    & !c-the-velocity-of-light in au
                        nm  = 10/0.52917706d0,   & !1 nanometer in au
                        pi  = 3.14159265358979323846D0 !acos(-1.d0)

   !field component lables for printing
   character*2:: fc(3) = (/'Fx','Fy','Fz'/), &
                 bc(3) = (/'Bx','By','Bz'/)

#ifndef VAR_LINSCA
! DALTON_AO_RSP and PRG_DIRAC take lupri from common block
#include "priunit.h"
#endif

contains

   subroutine prop_test_gradient(S, D, F, ng)
      type(matrix), intent(in) :: S, D, F
      integer,      intent(in) :: ng !number of geometrical coordinates
      complex(8)   :: Eg(ng) !NB: zero first
      type(matrix) :: DFD
#ifdef PRG_DIRAC
      integer :: iatom, ix, iy, iz
#endif

      ! contract first-order geometry-differentiated integrals
      ! with unperturbed density D and energy-weighted density DFD
      DFD = D*F*D
      Eg = 0 !zero first, since prop_one/twoave are incremental
      call prop_oneave(S, (/'GEO'/), (/D/), (/ng/), Eg, DFD=(/DFD/))
      call prop_twoave((/'GEO'/), (/D/), (/ng/), Eg)
      ! print
#ifndef PRG_DIRAC
      call print_tensor((/ng/), Eg, 'gradient = Eg = E0g - Sg DFD')
#endif
      ! free DFD
      DFD = 0

#ifdef PRG_DIRAC
!     this is needed in DIRAC tests
      call header('molecular gradient', -1)
  
      do iatom = 1, nr_atoms
        ix = 3*(iatom - 1) + 1
        iy = 3*(iatom - 1) + 2
        iz = 3*(iatom - 1) + 3
        write(*, '(a6, 3e26.16)') sym_dep_nuc_name(iatom), real(Eg(ix)), real(Eg(iy)), real(Eg(iz))
      end do
  
      call prsymb(unit_output, '-', 46, 0)
#endif

   end subroutine

   subroutine prop_test_Gprime_Df(decomp, S, D, F, idum, freq) !Gprime_Df(-freq)
      type(decompItem), intent(inout) :: decomp
      type(matrix), intent(in) :: S, D, F
      integer,      intent(in) :: idum    !not used
      complex(8),   intent(in) :: freq(:) !frequency
      complex(8)   :: Efb(3,3) = 0
      type(matrix) :: Df(3), Ff(3), DFDf(3)
      integer      :: i, j
      ! the following contribution is zero
      ! call prop_oneave(S, (/'EL ','MAG'/), (/D/), shape(Efb), Efb)
      ! solve EL response equations for the negative frequency
      call pert_dens(decomp, S, (/'EL'/), int(shape(Df)), (/D/), (/F/), &
                     Df, Ff, freq=-freq(1:1))
      ! compute the EL-perturbed energy-weighted density matrices,
      ! to be contracted against MAG-perturbed overlap integrals
      do i = 1, 3
         DFDf(i) = Df(i)*(F - freq(1)/2*S)*D + D*Ff(i)*D &
                     + D*(F + freq(1)/2*S)*Df(i)
      end do
      ! contract Df and DFDf with MAG-perturbed integrals
      call prop_oneave(S, (/'MAG'/), (/Df/), int(shape(Efb)), Efb, &
                       perm=(/2,1/), freq=freq(1:1), DFD=DFDf)
      call prop_twoave((/'MAG'/), (/D,Df/), int(shape(Efb)), Efb, perm=(/2,1/))
      ! print
      call print_tensor(int(shape(Efb)), Efb, 'G-prime = Efb = E1b Df - Sb DFDf')
      ! free Df(:), Ff(:), DFDf(:)
      Df=0; Ff=0; DFDf=0
   end subroutine


   subroutine prop_test_Gprime_Db(decomp, S, D, F, idum, freq) !Gprime_Db(+freq)
      type(decompItem), intent(inout) :: decomp
      type(matrix), intent(in) :: S, D, F
      integer,      intent(in) :: idum    !not used
      complex(8),   intent(in) :: freq(:) !frequency
      complex(8)   :: Efb(3,3) = 0
      type(matrix) :: Db(3), Fb(3)
      integer      :: i
      ! the following contribution is zero
      call prop_oneave(S, (/'EL ','MAG'/), (/D/), int(shape(Efb)), Efb)
      ! solve EL response equations for the negative frequency
      call pert_dens(decomp, S, (/'MAG'/), int(shape(Db)), (/D/), (/F/), &
                     Db, Fb, freq=freq(1:1))
      ! contract Db with -dipole integrals
      call prop_oneave(S, (/'EL'/), (/Db/), int(shape(Efb)), Efb)
      ! print
      call print_tensor(int(shape(Efb)), Efb, 'G-prime = Efb = E1f Db')
      ! free
      Db = 0;  Fb = 0
   end subroutine


   subroutine elec_dipmom(S, D, F, Ef)
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(out) :: Ef(3)
      Ef = 0
      call prop_oneave(S, (/'EL'/), (/D/), (/3/), Ef)
      call print_tensor((/3/), -Ef, 'mu = -Ef')
      call print_tensor((/3/), -Ef, 'mu = -Ef', unit=lupri)
   end subroutine


   subroutine elec_quadrupole(S, D, F, t)
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(out) :: t(6)
      t = 0.0d0
      call prop_oneave(S, (/'ELGR'/), (/D/), (/6/), t)
      call print_tensor((/6/), -t, 'quadrupole = -t')
   end subroutine


   subroutine elec_polariz(decomp, S, D, F, freq, Eff)
      type(decompItem), intent(inout) :: decomp
      type(matrix),      intent(in)  :: S, D, F
      complex(8),        intent(in)  :: freq
      complex(8),        intent(out) :: Eff(3,3)
      type(matrix) :: Df(3), Ff(3)
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq/))
      Eff = 0 !zero, as prop_oneave works incrementally
      call prop_oneave(S, (/'EL'/), (/Df/), (/3,3/), Eff)
      call print_tensor((/3,3/), -Eff, 'alpha = -Eff')
      call print_tensor((/3,3/), -Eff, 'alpha = -Eff', unit=decomp%lupri)
      !delete response matrices
      Df = 0; Ff = 0
      !change sign, since polarizability is minus quasi-energy derivative
      Eff = -Eff
   end subroutine


   subroutine elec_hypolar(decomp, S, D, F, freq, Efff)
      type(decompItem), intent(inout) :: decomp
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq(3)
      complex(8),   intent(out) :: Efff(3,3,3)
      type(matrix) :: Df(3), De(3), Dfe(3,3)
      type(matrix) :: Ff(3), Fe(3), Ffe(3,3)
      integer      :: i, j
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/elec_hypolar: sum(freq) should be zero!')
      !solve equations
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=freq(2:2))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=freq(3:3))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,Df,De/), (/F,Ff,Fe/), &
                     Dfe, Ffe, freq=(/freq(2),freq(3)/))
      !just one energy contribution, no reort or multiplier contributions
      Efff=0
      call prop_oneave(S, (/'EL'/), (/Dfe/), (/3,3,3/), Efff)
      !print
      call print_tensor((/3,3,3/), -Efff, 'beta = -Efff', freq)
      call print_tensor((/3,3,3/), -Efff, 'beta = -Efff', freq, unit=decomp%lupri)
      !free
      Df=0; Ff=0; De=0; Fe=0; Dfe=0; Ffe=0
      !change sign to get beta
      Efff = -Efff
   end subroutine


   subroutine alt_elec_hypol(decomp, S, D, F, freq, Efff)
      type(decompItem), intent(inout) :: decomp
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq(3)
      complex(8),   intent(out) :: Efff(3,3,3)
      type(matrix) :: De(3), Df(3), Dg(3), DeSD(3), FDSfg
      type(matrix) :: Fe(3), Ff(3), Fg(3), FeDS(3), DSDfg
      character(4) :: UNP(0)
      integer      :: i, j, k
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/alt_elec_hypol: sum(w) should be zero!')
      !solve equations
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=freq(1:1))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=freq(2:2))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Dg, Fg, freq=freq(3:3))
      !no energy contributions for HF, one for KSDFT
      Efff = 0
      call prop_twoave(UNP, (/D,De,Df,Dg,(0d0*D,i=1,54)/), (/3,3,3/), Efff)
      !gradient Lagrange multiplier contribution -tr (DeSD-) (FDS-)fg
      do i = 1, 3
         DeSD(i) = De(i)*S*D - D*S*De(i)
      end do
      do k = 1, 3
         do j = 1, 3
            FDSfg = (Fe(k)*Df(j) + Ff(j)*De(k)) * S &
              - S * (De(k)*Ff(j) + Df(j)*Fe(k))
            do i = 1, 3
               Efff(i,j,k) = Efff(i,j,k) - tr(DeSD(i),FDSfg)
            end do; FDSfg=0
         end do
      end do; DeSD=0
      !idempotency Lagrange multiplier contribution -tr (FeDS+) (DSD)fg
      do i = 1, 3
         FeDS(i) = Fe(i)*D*S + S*D*Fe(i) - Fe(i)
      end do
      do k = 1, 3
         do j = 1, 3
            DSDfg = De(k)*S*Df(j) + Df(j)*S*De(k)
            do i = 1, 3
               Efff(i,j,k) = Efff(i,j,k) - tr(FeDS(i),DSDfg)
            end do; DSDfg=0
         end do
      end do; FeDS=0
      !print
      call print_tensor((/3,3,3/), -Efff, 'beta = -Efff', freq)
      call print_tensor((/3,3,3/), -Efff, 'beta = -Efff', freq, unit=decomp%lupri)
      !free
      De=0; Fe=0; Df=0; Ff=0; Dg=0; Fg=0
      !change sign for beta
      Efff = -Efff
   end subroutine


   subroutine elec_sechyp(decomp, S, D, F, freq, Effff)
      type(decompItem), intent(inout) :: decomp
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq(4)
      complex(8),   intent(out) :: Effff(3,3,3,3)
      type(matrix) :: De(3), Df(3), Dg(3), Def(3,3),   &
                      Fe(3), Ff(3), Fg(3), Fef(3,3),   &
                      Deg(3,3), Dfg(3,3), Defg(3,3,3), &
                      Feg(3,3), Ffg(3,3), Fefg(3,3,3)
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/elec_sechyp: sum(freq) should be zero!')
      !solve equations
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=(/freq(2)/))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq(3)/))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Dg, Fg, freq=(/freq(4)/))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,De,Df/), (/F,Fe,Ff/), &
                     Def, Fef, freq=(/freq(2),freq(3)/))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,De,Dg/), (/F,Fe,Fg/), &
                     Deg, Feg, freq=(/freq(2),freq(4)/))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,Df,Dg/), (/F,Ff,Fg/), &
                     Dfg, Ffg, freq=(/freq(3),freq(4)/))
      call pert_dens(decomp, S, (/'EL','EL','EL'/), (/3,3,3/), &
                     (/D,De,Df,Dg,Def,Deg,Dfg/), &
                     (/F,Fe,Ff,Fg,Fef,Feg,Ffg/), &
                     Defg, Fefg, freq=(/freq(2),freq(3),freq(4)/))
      !just one property contribution
      Effff=0
      call prop_oneave(S, (/'EL'/), (/Defg/), (/3,3,3,3/), Effff)
      !print
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq)
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq, unit=decomp%lupri)
      !free
      De=0; Df=0; Dg=0; Def=0; Deg=0; Dfg=0; Defg=0
      Fe=0; Ff=0; Fg=0; Fef=0; Feg=0; Ffg=0; Fefg=0
      !change sign to get gamma
      Effff = -Effff
   end subroutine

   subroutine alt_elec_sechyp(decomp, S, D, F, freq, Effff)
      type(decompItem), intent(inout) :: decomp
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq(4)
      complex(8),   intent(out) :: Effff(3,3,3,3)
      type(matrix) :: De(3), Df(3), Dg(3), Dh(3), DeSD(3), &
                      Dfg(3,3), Dfh(3,3), Dgh(3,3), FDSfgh, &
                      Fe(3), Ff(3), Fg(3), Fh(3), FeDS(3), &
                      Ffg(3,3), Ffh(3,3), Fgh(3,3), DSDfgh
      character(4) :: UNP(0)
      integer      :: i, j, k, l
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/alt_elec_sechyp: sum(w) should be zero!')
      !solve equations
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=(/freq(1)/))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq(2)/))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Dg, Fg, freq=(/freq(3)/))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Dh, Fh, freq=(/freq(4)/))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,Df,Dg/), (/F,Ff,Fg/), &
                     Dfg, Ffg, freq=(/freq(2),freq(3)/))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,Df,Dh/), (/F,Ff,Fh/), &
                     Dfh, Ffh, freq=(/freq(2),freq(4)/))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,Dg,Dh/), (/F,Fg,Fh/), &
                     Dgh, Fgh, freq=(/freq(3),freq(4)/))
      ! one energy contribution
      Effff=0
      call prop_twoave(UNP, (/D,De,Df,Dg,Dh, &
                              (0d0*D,i=1,18),Dfg,(0d0*D,i=1,9), &
                              Dfh,Dgh,(0d0*D,i=1,189)/), &
                       int(shape(Effff)), Effff)
      ! gradient Lagrange multiplier contribution -tr (DeSD-) (FDS-)fgh
      do i = 1, 3
         DeSD(i) = De(i)*S*D - D*S*De(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               FDSfgh = (Ff(j)*Dgh(k,l) + Fg(k)*Dfh(j,l) &
                      +  Fh(l)*Dfg(j,k) + Ffg(j,k)*Dh(l) &
                      +  Ffh(j,l)*Dg(k) + Fgh(k,l)*Df(j)) * S &
                  - S * (Df(j)*Fgh(k,l) + Dg(k)*Ffh(j,l) &
                      +  Dh(l)*Ffg(j,k) + Dfg(j,k)*Fh(l) &
                      +  Dfh(j,l)*Fg(k) + Dgh(k,l)*Ff(j))
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(DeSD(i),FDSfgh)
               end do; FDSfgh=0
            end do
         end do
      end do; DeSD=0
      ! idempotency Lagrange multiplier contribution -tr (FhDS+) (DSD)efg
      do i = 1, 3
         FeDS(i) = Fe(i)*D*S + S*D*Fe(i) - Fe(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               DSDfgh = Df(j)*S*Dgh(k,l) + Dg(k)*S*Dfh(j,l) &
                      + Dh(l)*S*Dfg(j,k) + Dfg(j,k)*S*Dh(l) &
                      + Dfh(j,l)*S*Dg(k) + Dgh(k,l)*S*Df(j)
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(FeDS(i),DSDfgh)
               end do; DSDfgh=0
            end do
         end do
      end do; FeDS=0
      ! print
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq)
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq, unit=decomp%lupri)
      !free
      De=0; Df=0; Dg=0; Dh=0; Dfg=0; Dfh=0; Dgh=0
      Fe=0; Ff=0; Fg=0; Fh=0; Ffg=0; Ffh=0; Fgh=0
      !change sign to get gamma
      Effff = -Effff
   end subroutine


   subroutine alt2_elec_sechyp(decomp, S, D, F, freq, Effff)
      type(decompItem), intent(inout) :: decomp
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq(4)
      complex(8),   intent(out) :: Effff(3,3,3,3)
      type(matrix) :: De(3), Df(3), Dg(3), Dh(3), DSD,      &
                      Fe(3), Ff(3), Fg(3), Fh(3), FDS,      &
                      Def(3,3), Deg(3,3), Deh(3,3),         &
                      Fef(3,3), Feg(3,3), Feh(3,3),         &
                      DeSDf(3,3), FDSgh, FeDSf(3,3), DSDgh, &
                      DeSDg(3,3), FDSfh, FeDSg(3,3), DSDfh, &
                      DeSDh(3,3), FDSfg, FeDSh(3,3), DSDfg
      integer      :: i, j, k, l
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/alt_elec_sechyp: sum(freq) should be zero!')
      !solve equations
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=(/freq(1)/))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq(2)/))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Dg, Fg, freq=(/freq(3)/))
      call pert_dens(decomp, S, (/'EL'/), (/3/), (/D/), (/F/), Dh, Fh, freq=(/freq(4)/))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,De,Df/), (/F,Fe,Ff/), &
                     Def, Fef, freq=(/freq(1),freq(2)/))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,De,Dg/), (/F,Fe,Fg/), &
                     Deg, Feg, freq=(/freq(1),freq(3)/))
      call pert_dens(decomp, S, (/'EL','EL'/), (/3,3/), (/D,De,Dh/), (/F,Fe,Fh/), &
                     Deh, Feh, freq=(/freq(1),freq(4)/))
      !no energy contributions
      Effff=0
      !gradient Lagrange multiplier contribution -tr (DeSD-)f (FDS-)gh
      do j = 1, 3
         do i = 1, 3
            DeSDf(i,j) = Def(i,j)*S*D  - D*S*Def(i,j)  &
                       + De(i)*S*Df(j) - Df(j)*S*De(i)
         end do
      end do
      do l = 1, 3
         do k = 1, 3
            FDSgh = (Fg(k)*Dh(l) + Fh(l)*Dg(k)) * S &
              - S * (Dh(l)*Fg(k) + Dg(k)*Fh(l))
            do j = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(DeSDf(i,j),FDSgh)
               end do
            end do; FDSgh=0
         end do
      end do; DeSDf=0
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq)
      !idempotency Lagrange multiplier contribution -tr (FeDS+)f (DSD)gh
      do j = 1, 3
         do i = 1, 3
            FeDSf(i,j) = (Fef(i,j)*D + Fe(i)*Df(j)) * S &
                   + S * (D*Fef(i,j) + Df(j)*Fe(i)) - Fef(i,j)
         end do
      end do
      do l = 1, 3
         do k = 1, 3
            DSDgh = Dg(k)*S*Dh(l) + Dh(l)*S*Dg(k)
            do j = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(FeDSf(i,j),DSDgh)
               end do
            end do; DSDgh=0
         end do
      end do; FeDSf=0
      !gradient Lagrange multiplier contribution -tr (DeSD-)g (FDS-)fh
      do k = 1, 3
         do i = 1, 3
            DeSDg(i,k) = Deg(i,k)*S*D  - D*S*Deg(i,k)  &
                       + De(i)*S*Dg(k) - Dg(k)*S*De(i)
         end do
      end do
      do l = 1, 3
         do j = 1, 3
            FDSfh = (Ff(j)*Dh(l) + Fh(l)*Df(j)) * S &
              - S * (Df(j)*Fh(l) + Dh(l)*Ff(j))
            do k = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(DeSDg(i,k),FDSfh)
               end do
            end do; FDSfh=0
         end do
      end do; DeSDg=0
      !idempotency Lagrange multiplier contribution -tr (FeDS+)g (DSD)fh
      do k = 1, 3
         do i = 1, 3
            FeDSg(i,k) = (Feg(i,k)*D + Fe(i)*Dg(k)) * S &
                   + S * (D*Feg(i,k) + Dg(k)*Fe(i)) - Feg(i,k)
         end do
      end do
      do l = 1, 3
         do j = 1, 3
            DSDfh = Df(j)*S*Dh(l) + Dh(l)*S*Df(j)
            do k = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(FeDSg(i,k),DSDfh)
               end do
            end do; DSDfh=0
         end do
      end do; FeDSg=0
      !gradient Lagrange multiplier contribution -tr (DeSD-)h (FDS-)fg
      do l = 1, 3
         do i = 1, 3
            DeSDh(i,l) = Deh(i,l)*S*D  - D*S*Deh(i,l)  &
                       + De(i)*S*Dh(l) - Dh(l)*S*De(i)
         end do
      end do
      do k = 1, 3
         do j = 1, 3
            FDSfg = (Ff(j)*Dg(k) + Fg(k)*Df(j)) * S &
              - S * (Df(j)*Fg(k) + Dg(k)*Ff(j))
            do l = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(DeSDh(i,l),FDSfg)
               end do
            end do; FDSfg=0
         end do
      end do; DeSDh=0
      !idempotency Lagrange multiplier contribution -tr (FeDS+)h (DSD)fg
      do l = 1, 3
         do i = 1, 3
            FeDSh(i,l) = (Feh(i,l)*D + Fe(i)*Dh(l)) * S &
                   + S * (D*Feh(i,l) + Dh(l)*Fe(i)) - Feh(i,l)
         end do
      end do
      do k = 1, 3
         do j = 1, 3
            DSDfg = Df(j)*S*Dg(k) + Dg(k)*S*Df(j)
            do l = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(FeDSh(i,l),DSDfg)
               end do
            end do; DSDfg=0
         end do
      end do; FeDSh=0
      !print
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq)
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq, unit=decomp%lupri)
      !free
      De=0; Df=0; Dg=0; Dh=0; Def=0; Deg=0; Deh=0
      Fe=0; Ff=0; Fg=0; Dh=0; Fef=0; Feg=0; Feh=0
      !change sign to get gamma
      Effff = -Effff
   end subroutine


   subroutine magnetiz(decomp, S, D, F, freq, Eoo, Ebb)
   !no-Lon and London magnetizability
      type(decompItem), intent(inout) :: decomp
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq
      complex(8),   intent(out) :: Eoo(3,3), Ebb(3,3)
      type(matrix) :: Db(3), Fb(3), DFD, DFDb(3)
      integer      :: i, j
      !no-London
      call pert_dens(decomp, S, (/'MAGO'/), (/3/), (/D/), (/F/), Db, Fb, freq=(/freq/))
      Eoo = 0
      call prop_oneave(S, (/'MAGO','MAGO'/), (/D/), (/3,3/), Eoo)
      !call print_tensor( (/3,3/), Eoo, 'E0oo'); Eoo=0
      call prop_oneave(S, (/'MAGO'/), (/Db/), (/3,3/), Eoo)
      !call print_tensor( (/3,3/), Eoo, 'E1oDo'); Eoo=0
      call print_tensor((/3,3/), Eoo, 'no-London Magnetizability = Eoo', (/-freq,freq/))
      Db=0; Fb=0
#ifdef PRG_DIRAC
      Ebb = 0
#else
      !London diamag
      call pert_dens(decomp, S, (/'MAG'/), (/3/), (/D/), (/F/), Db, Fb, freq=(/freq/))
      Ebb = 0
      call prop_twoave((/'MAG','MAG'/), (/D/), (/3,3/), Ebb)
      DFD = D*F*D
      call prop_oneave(S, (/'MAG','MAG'/), (/D/), (/3,3/), Ebb, &
                       DFD=(/DFD/), freq=(/-freq,freq/))
      DFD = 0
      !call print_tensor( (/3,3/), Ebb, 'E0bb-i/2TbbD-SbbW'); Ebb=0
      !London paramag
      call prop_twoave((/'MAG'/), (/D,Db/), (/3,3/), Ebb)
      do j = 1, 3
         DFDb(j) = Db(j)*(F+(freq/2)*S)*D + D*Fb(j)*D &
                 +     D*(F-(freq/2)*S)*Db(j)
      end do
      call prop_oneave(S, (/'MAG'/), (/Db/), (/3,3/), Ebb, &
                       DFD=DFDb, freq=(/-freq/))
      Db=0; Fb=0; DFDb=0
      !call print_tensor( (/3,3/), Ebb, 'E1bDb-i/2TbDb-SbDFDb'); Ebb=0
      call print_tensor((/3,3/), Ebb, 'London Magnetizability = Ebb', (/-freq,freq/))
#endif /* ifdef PRG_DIRAC */
   end subroutine


  subroutine test_grcont(n, nr_atoms, test_london, test_geo)

    integer              :: n, nr_atoms
    logical              :: test_london, test_geo
!   ----------------------------------------------------------------------------
    real(8), allocatable :: work(:)
    integer              :: lwork

    integer              :: i, j, k, l, m
#ifdef PRG_DIRAC
    real(8)              :: D(n, n, 8)
#else                   
    real(8)              :: D(n, n, 2)
#endif /* PRG_DIRAC */
    real(8)              :: F(n, n, 3), g(3, nr_atoms)

#ifdef PRG_DIRAC
#include "dcbbas.h"
#endif /* PRG_DIRAC */

    lwork = 10000*n + 50*n*n
    allocate(work(lwork))

    if (test_london) then

      print '(a)', 'London 2e integrals'
      print '(a, 12x, a2, 18x, a2, 18x, a2)', '  i  j  k  l', 'Bx', 'By', 'Bz'

!     do l = 1, n
!       do k = 1, n
      do l = 1, 1
        do k = 1, 1

          D(:, :, 1) = 0.0d0
          D(k, l, 1) = 1.0d0
          F          = 0.0d0

          call grcont(work, lwork, F, n*n*3, &
                      .false., .true., 1, 0, .false., .true., D(:, :, 1), 1)

          do j = 1, n
            do i = 1, n
              write(0, '(4i3, 3f20.12)') i, j, k, l, F(i, j, :)
            end do
          end do
        end do
      end do
    end if

    if (test_geo) then

      print '(a)', 'Geometry 2e integrals'
      print '(a, 12x, a2, 18x, a2, 18x, a2)', '  i  j  k  l  a', 'Rx', 'Ry', 'Rz'

      D = 0.0d0

!     do l = 1, n
!       do k = 1, n
      do l = 1, 1
        do k = 1, 1
#ifdef PRG_DIRAC
          do j = 1, ntbas(1)
            do i = 1, ntbas(1)
#else
          do j = 1, n
            do i = 1, n
#endif /* PRG_DIRAC */

              D(k, l, 1) = 1.0d0
#ifdef PRG_DIRAC
              D(i, j, 5) = 1.0d0
#else
              D(i, j, 2) = 1.0d0
#endif /* PRG_DIRAC */
              g          = 0.0d0

              call grcont(work, lwork, g, 3*nr_atoms, &
                          .true., .false., 1, 0, .true., .false., D, 2)

              D(k, l, 1) = 0.0d0
#ifdef PRG_DIRAC
              D(i, j, 5) = 0.0d0
#else
              D(i, j, 2) = 0.0d0
#endif /* PRG_DIRAC */

              do m = 1, nr_atoms
                write(0, '(5i3, 3f20.12)') i, j, k, l, m, g(:, m)
              end do
            end do
          end do
        end do
      end do
    end if

    deallocate(work)

  end subroutine


   subroutine print_tensor(dims, tensor, title, freqs, colwidth, unit)
   !ajt This is used for printing response function tensors
   !TODO: Add comp. lables argument. Make space between blocks of rows
      integer,                intent(in) :: dims(:)
      complex(8),             intent(in) :: tensor(*) !(product(dims))
      character(*), optional, intent(in) :: title
      integer,      optional, intent(in) :: colwidth, unit
      complex(8),   optional, intent(in) :: freqs(:)
      integer     :: uni, siz, dec, cwid, i
      character*9 :: fmt
      uni = 6
      if (present(unit)) uni = unit
      if (present(title)) then
         if (present(freqs)) then
            fmt = '(a,nf6.3)'
            write (fmt(4:4),'(i1)') size(freqs)
            write (uni,fmt) title // ' w:', dreal(freqs)
         else
            write (uni,'(a)') title
         end if
      end if
      siz = product(dims)
      dec = max(2, 1+ceiling(log10(maxval(abs(tensor(1:siz))))))
      !dec = max(dec,ceiling(log10( maxval(tensor(1:siz)))))
      cwid = 15
      if (present(colwidth)) cwid = max(colwidth,dec+2)
      dec = cwid - dec - 1
      fmt = '(fww.dd)'
      write (fmt(3:7),'(i2,a1,i2)') cwid, '.', dec
      if (any(dimag(tensor(1:siz)) /= 0)) write (uni,'(a)') '    (real part)'
      call subr(dims, dreal(tensor(1:siz)), uni, cwid, fmt)
      if (any(dimag(tensor(1:siz)) /= 0)) then
         write (uni,'(a)') '    (imaginary part)'
         call subr(dims, dimag(tensor(1:siz)), uni, cwid, fmt)
      end if
      write (uni,'()') !final blank line
   contains
      subroutine subr(dims, tnz, unit, cwid, fmt)
         integer,     intent(in) :: dims(:), unit, cwid
         real(8),     intent(in) :: tnz(*)
         character*8, intent(in) :: fmt
         integer :: i, j, l
         character*(dims(1)*(cwid+1)) :: line
         do i = 1, product(dims), dims(1)
            l = 1
            do j = 1, dims(1)
               line(l:l) = ' '
               write (line(l+1:l+cwid), fmt) tnz(i+j-1)
               l = l + 1 + cwid
            end do
            write (unit,'(a)') line
         end do
      end subroutine
   end subroutine


!radovan: please do not leave commented-out code
!         (without describing the reason why it is commented out)
!         in the repository

!         function RandomReal() result(r)
!         !Random real numbers with Gaussian distribution
!         !The generator is seeded on the first call
!         !Two reals are generated each time, and returned in consecutive calls
!         !Random_Seed and Random_Number are provided by Intel and GNU compilers
!         !To get different numbers, change seed.
!             real(8),parameter :: pi=3.14159265358979323846d00
!             real(8) :: r,a,x,y
!             real(8),save:: rnext
!             logical,save:: first=.true.,flip=.true.
!             if (first) then !reset fortran's random generator
!                 call Random_Seed(put=(/628468422,281827387/))
!                 first = .false.
!             end if
!             r = rnext
!             flip = .not.flip
!             if (flip) return !return sav
!             call Random_Number(x)
!             a = 2*Pi*x !angle in [0,2pi]
!             call Random_Number(x)
!             y = sqrt(-log(x))
!             rnext = cos(a)*y
!             r = sin(a)*y
!         end function

#ifdef PRG_DIRAC
   function element_ooff(S, D, F, w, h, i, j, k)

!     --------------------------------------------------------------------------
      real(8)                     :: element_ooff
      type(matrix), intent(inout) :: S, D, F
      integer,      intent(in)    :: h, i, j, k
      complex(8),   intent(in)    :: w(3)
!     --------------------------------------------------------------------------
      type(matrix)                :: D_i(3)
      type(matrix)                :: D_j(3)
      type(matrix)                :: D_k(3)
      type(matrix)                :: D_ij(3, 3)
      type(matrix)                :: D_ik(3, 3)
      type(matrix)                :: D_jk(3, 3)
      type(matrix)                :: D_ijk(3, 3, 3)
      type(matrix)                :: F_i(3)
      type(matrix)                :: F_j(3)
      type(matrix)                :: F_k(3)
      type(matrix)                :: F_ij(3, 3)
      type(matrix)                :: F_ik(3, 3)
      type(matrix)                :: F_jk(3, 3)
      type(matrix)                :: F_ijk(3, 3, 3)
      complex(8)                  :: e(3, 3, 3, 3)
      integer                     :: c1, c2, c3
      type(decompitem)            :: decomp
!     --------------------------------------------------------------------------

      S = 1.0d0*S
      D = 1.0d0*D
      F = 1.0d0*F
     
      do c1 = 1, 3
         D_i(c1) = tiny(0.0d0)*D
         D_j(c1) = tiny(0.0d0)*D
         D_k(c1) = tiny(0.0d0)*D
         F_i(c1) = tiny(0.0d0)*D
         F_j(c1) = tiny(0.0d0)*D
         F_k(c1) = tiny(0.0d0)*D
         do c2 = 1, 3
            D_ij(c1, c2) = tiny(0.0d0)*D
            D_ik(c1, c2) = tiny(0.0d0)*D
            D_jk(c1, c2) = tiny(0.0d0)*D
            F_ij(c1, c2) = tiny(0.0d0)*D
            F_ik(c1, c2) = tiny(0.0d0)*D
            F_jk(c1, c2) = tiny(0.0d0)*D
            do c3 = 1, 3
               D_ijk(c1, c2, c3) = tiny(0.0d0)*D
               F_ijk(c1, c2, c3) = tiny(0.0d0)*D
            end do
         end do
      end do

      call pert_dens(decomp, S, (/'MAGO'/), (/1/), (/D/), (/F/), D_i(i), F_i(i), comp=(/i/), freq=(/w(1)/))
      call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_j(j), F_j(j), comp=(/j/), freq=(/w(2)/))
      call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_k(k), F_k(k), comp=(/k/), freq=(/w(3)/))

      call pert_dens(decomp,             &
                     S,                  &
                   (/'MAGO', 'EL  '/),       &
                   (/1, 1/),             &
                   (/                    &
                     D,                  &
                     D_i(i),             &
                     D_j(j)              &
                    /),                  &
                   (/                    &
                     F,                  &
                     F_i(i),             &
                     F_j(j)              &
                    /),                  &
                     D_ij(i, j),         &
                     F_ij(i, j),         &
                     comp=(/i, j/),      &
                     freq=(/w(1), w(2)/))
  
      call pert_dens(decomp,             &
                     S,                  &
                   (/'MAGO', 'EL  '/),       &
                   (/1, 1/),             &
                   (/                    &
                     D,                  &
                     D_i(i),             &
                     D_k(k)              &
                    /),                  &
                   (/                    &
                     F,                  &
                     F_i(i),             &
                     F_k(k)              &
                    /),                  &
                     D_ik(i, k),         &
                     F_ik(i, k),         &
                     comp=(/i, k/),      &
                     freq=(/w(1), w(3)/))
  
      call pert_dens(decomp,             &
                     S,                  &
                   (/'EL', 'EL'/),       &
                   (/1, 1/),             &
                   (/                    &
                     D,                  &
                     D_j(j),             &
                     D_k(k)              &
                    /),                  &
                   (/                    &
                     F,                  &
                     F_j(j),             &
                     F_k(k)              &
                    /),                  &
                     D_jk(j, k),         &
                     F_jk(j, k),         &
                     comp=(/j, k/),      &
                     freq=(/w(2), w(3)/))
  
      call pert_dens(decomp,             &
                     S,                  &
                   (/'MAGO', 'EL  ', 'EL  '/), &
                   (/1, 1, 1/),          &
                   (/                    &
                     D,                  &
                     D_i(i),             &
                     D_j(j),             &
                     D_k(k),             &
                     D_ij(i, j),         &
                     D_ik(i, k),         &
                     D_jk(j, k)          &
                    /),                  &
                   (/                    &
                     F,                  &
                     F_i(i),             &
                     F_j(j),             &
                     F_k(k),             &
                     F_ij(i, j),         &
                     F_ik(i, k),         &
                     F_jk(j, k)          &
                    /),                  &
                     D_ijk(i, j, k),     &
                     F_ijk(i, j, k),     &
                     comp=(/i, j, k/),   &
                     freq=(/w(1), w(2), w(3)/))

      e = 0.0d0
      call prop_oneave(S, (/'MAGO'/), (/D_ijk/), (/3, 3, 3, 3/), e)
  
      D_i   = 0
      D_j   = 0
      D_k   = 0
      D_ij  = 0
      D_ik  = 0
      D_jk  = 0
      D_ijk = 0
      F_i   = 0
      F_j   = 0
      F_k   = 0
      F_ij  = 0
      F_ik  = 0
      F_jk  = 0
      F_ijk = 0
  
      element_ooff = -e(h, i, j, k)
      print '("@ ooff (w = "4f9.4")" a10 a2 a2 a2 f26.10)', real(-sum(w)), &
                                                            real(w(1)),    &
                                                            real(w(2)),    &
                                                            real(w(3)),    &
                                                            fc(h),         &
                                                            fc(i),         &
                                                            fc(j),         &
                                                            fc(k),         &
                                                            element_ooff

   end function


  function gamma_element(S, D, F, w, h, i, j, k)

!   ----------------------------------------------------------------------------
    real(8)                     :: gamma_element
    type(matrix), intent(inout) :: S, D, F
    integer,      intent(in)    :: h, i, j, k
    complex(8),   intent(in)    :: w(3)
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3)
    type(matrix)                :: D_j(3)
    type(matrix)                :: D_k(3)
    type(matrix)                :: D_ij(3, 3)
    type(matrix)                :: D_ik(3, 3)
    type(matrix)                :: D_jk(3, 3)
    type(matrix)                :: D_ijk(3, 3, 3)
    type(matrix)                :: F_i(3)
    type(matrix)                :: F_j(3)
    type(matrix)                :: F_k(3)
    type(matrix)                :: F_ij(3, 3)
    type(matrix)                :: F_ik(3, 3)
    type(matrix)                :: F_jk(3, 3)
    type(matrix)                :: F_ijk(3, 3, 3)
    complex(8)                  :: e(3, 3, 3, 3)
    integer                     :: c1, c2, c3
    type(decompitem)            :: decomp
!   ----------------------------------------------------------------------------

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    do c1 = 1, 3
      D_i(c1) = tiny(0.0d0)*D
      D_j(c1) = tiny(0.0d0)*D
      D_k(c1) = tiny(0.0d0)*D
      F_i(c1) = tiny(0.0d0)*D
      F_j(c1) = tiny(0.0d0)*D
      F_k(c1) = tiny(0.0d0)*D
      do c2 = 1, 3
        D_ij(c1, c2) = tiny(0.0d0)*D
        D_ik(c1, c2) = tiny(0.0d0)*D
        D_jk(c1, c2) = tiny(0.0d0)*D
        F_ij(c1, c2) = tiny(0.0d0)*D
        F_ik(c1, c2) = tiny(0.0d0)*D
        F_jk(c1, c2) = tiny(0.0d0)*D
        do c3 = 1, 3
          D_ijk(c1, c2, c3) = tiny(0.0d0)*D
          F_ijk(c1, c2, c3) = tiny(0.0d0)*D
        end do
      end do
    end do

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_i(i), F_i(i), comp=(/i/), freq=(/w(1)/))
    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_j(j), F_j(j), comp=(/j/), freq=(/w(2)/))
    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_k(k), F_k(k), comp=(/k/), freq=(/w(3)/))

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j)              &
                  /),                  &
                   D_ij(i, j),         &
                   F_ij(i, j),         &
                   comp=(/i, j/),      &
                   freq=(/w(1), w(2)/))

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_k(k)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_k(k)              &
                  /),                  &
                   D_ik(i, k),         &
                   F_ik(i, k),         &
                   comp=(/i, k/),      &
                   freq=(/w(1), w(3)/))

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_j(j),             &
                   D_k(k)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_j(j),             &
                   F_k(k)              &
                  /),                  &
                   D_jk(j, k),         &
                   F_jk(j, k),         &
                   comp=(/j, k/),      &
                   freq=(/w(2), w(3)/))

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j),             &
                   D_k(k),             &
                   D_ij(i, j),         &
                   D_ik(i, k),         &
                   D_jk(j, k)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j),             &
                   F_k(k),             &
                   F_ij(i, j),         &
                   F_ik(i, k),         &
                   F_jk(j, k)          &
                  /),                  &
                   D_ijk(i, j, k),     &
                   F_ijk(i, j, k),     &
                   comp=(/i, j, k/),   &
                   freq=(/w(1), w(2), w(3)/))

    e = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ijk/), (/3, 3, 3, 3/), e)

    D_i   = 0
    D_j   = 0
    D_k   = 0
    D_ij  = 0
    D_ik  = 0
    D_jk  = 0
    D_ijk = 0
    F_i   = 0
    F_j   = 0
    F_k   = 0
    F_ij  = 0
    F_ik  = 0
    F_jk  = 0
    F_ijk = 0

    gamma_element = -e(h, i, j, k)
    print '("@ gamma (w = "4f9.4")" a10 a2 a2 a2 f26.10)', real(-sum(w)), &
                                                           real(w(1)),    &
                                                           real(w(2)),    &
                                                           real(w(3)),    &
                                                           fc(h),         &
                                                           fc(i),         &
                                                           fc(j),         &
                                                           fc(k),         &
                                                           gamma_element

  end function

  function delta_element(S, D, F, w, h, i, j, k, l)

!   ----------------------------------------------------------------------------
    real(8)                     :: delta_element
    type(matrix), intent(inout) :: S, D, F
    integer,      intent(in)    :: h, i, j, k, l
    complex(8),   intent(in)    :: w(4)
!   ----------------------------------------------------------------------------
!   pascal triangle: 1 0th order matrix
!                    4 1st order matrices
!                    6 2nd
!                    4 3rd
!                    1 4th
    type(matrix)                :: D_i(3)
    type(matrix)                :: D_j(3)
    type(matrix)                :: D_k(3)
    type(matrix)                :: D_l(3)
    type(matrix)                :: D_ij(3, 3)
    type(matrix)                :: D_ik(3, 3)
    type(matrix)                :: D_il(3, 3)
    type(matrix)                :: D_jk(3, 3)
    type(matrix)                :: D_jl(3, 3)
    type(matrix)                :: D_kl(3, 3)
    type(matrix)                :: D_ijk(3, 3, 3)
    type(matrix)                :: D_ijl(3, 3, 3)
    type(matrix)                :: D_ikl(3, 3, 3)
    type(matrix)                :: D_jkl(3, 3, 3)
    type(matrix)                :: D_ijkl(3, 3, 3, 3)

    type(matrix)                :: F_i(3)
    type(matrix)                :: F_j(3)
    type(matrix)                :: F_k(3)
    type(matrix)                :: F_l(3)
    type(matrix)                :: F_ij(3, 3)
    type(matrix)                :: F_ik(3, 3)
    type(matrix)                :: F_il(3, 3)
    type(matrix)                :: F_jk(3, 3)
    type(matrix)                :: F_jl(3, 3)
    type(matrix)                :: F_kl(3, 3)
    type(matrix)                :: F_ijk(3, 3, 3)
    type(matrix)                :: F_ijl(3, 3, 3)
    type(matrix)                :: F_ikl(3, 3, 3)
    type(matrix)                :: F_jkl(3, 3, 3)
    type(matrix)                :: F_ijkl(3, 3, 3, 3)

    complex(8)                  :: e1(            3)
    complex(8)                  :: e2(         3, 3)
    complex(8)                  :: e3(      3, 3, 3)
    complex(8)                  :: e4(   3, 3, 3, 3)
    complex(8)                  :: e5(3, 3, 3, 3, 3)

    integer                     :: c1, c2, c3, c4
    type(decompitem)            :: decomp
!   ----------------------------------------------------------------------------

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    do c1 = 1, 3
      D_i(c1) = tiny(0.0d0)*D
      D_j(c1) = tiny(0.0d0)*D
      D_k(c1) = tiny(0.0d0)*D
      D_l(c1) = tiny(0.0d0)*D
      F_i(c1) = tiny(0.0d0)*D
      F_j(c1) = tiny(0.0d0)*D
      F_k(c1) = tiny(0.0d0)*D
      F_l(c1) = tiny(0.0d0)*D
      do c2 = 1, 3
        D_ij(c1, c2) = tiny(0.0d0)*D
        D_ik(c1, c2) = tiny(0.0d0)*D
        D_il(c1, c2) = tiny(0.0d0)*D
        D_jk(c1, c2) = tiny(0.0d0)*D
        D_jl(c1, c2) = tiny(0.0d0)*D
        D_kl(c1, c2) = tiny(0.0d0)*D
        F_ij(c1, c2) = tiny(0.0d0)*D
        F_ik(c1, c2) = tiny(0.0d0)*D
        F_il(c1, c2) = tiny(0.0d0)*D
        F_jk(c1, c2) = tiny(0.0d0)*D
        F_jl(c1, c2) = tiny(0.0d0)*D
        F_kl(c1, c2) = tiny(0.0d0)*D
        do c3 = 1, 3
          D_ijk(c1, c2, c3) = tiny(0.0d0)*D
          D_ijl(c1, c2, c3) = tiny(0.0d0)*D
          D_ikl(c1, c2, c3) = tiny(0.0d0)*D
          D_jkl(c1, c2, c3) = tiny(0.0d0)*D
          F_ijk(c1, c2, c3) = tiny(0.0d0)*D
          F_ijl(c1, c2, c3) = tiny(0.0d0)*D
          F_ikl(c1, c2, c3) = tiny(0.0d0)*D
          F_jkl(c1, c2, c3) = tiny(0.0d0)*D
          do c4 = 1, 3
            D_ijkl(c1, c2, c3, c4) = tiny(0.0d0)*D
            F_ijkl(c1, c2, c3, c4) = tiny(0.0d0)*D
          end do
        end do
      end do
    end do

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_i(i), F_i(i), comp=(/i/), freq=(/w(1)/))
    e2 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_i/), int(shape(e2)), e2)
    call print_polarizability(e2, (/w(1)/), i)

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_j(j), F_j(j), comp=(/j/), freq=(/w(2)/))
    e2 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_j/), int(shape(e2)), e2)
    call print_polarizability(e2, (/w(2)/), j)

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_k(k), F_k(k), comp=(/k/), freq=(/w(3)/))
    e2 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_k/), int(shape(e2)), e2)
    call print_polarizability(e2, (/w(3)/), k)

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_l(l), F_l(l), comp=(/l/), freq=(/w(4)/))
    e2 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_l/), int(shape(e2)), e2)
    call print_polarizability(e2, (/w(4)/), l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j)              &
                  /),                  &
                   D_ij(i, j),         &
                   F_ij(i, j),         &
                   comp=(/i, j/),      &
                   freq=(/w(1), w(2)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ij/), int(shape(e3)), e3)
    call print_polarizability(e3, (/w(1), w(2)/), i, j)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_k(k)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_k(k)              &
                  /),                  &
                   D_ik(i, k),         &
                   F_ik(i, k),         &
                   comp=(/i, k/),      &
                   freq=(/w(1), w(3)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ik/), int(shape(e3)), e3)
    call print_polarizability(e3, (/w(1), w(3)/), i, k)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_l(l)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_l(l)              &
                  /),                  &
                   D_il(i, l),         &
                   F_il(i, l),         &
                   comp=(/i, l/),      &
                   freq=(/w(1), w(4)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_il/), int(shape(e3)), e3)
    call print_polarizability(e3, (/w(1), w(4)/), i, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_j(j),             &
                   D_k(k)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_j(j),             &
                   F_k(k)              &
                  /),                  &
                   D_jk(j, k),         &
                   F_jk(j, k),         &
                   comp=(/j, k/),      &
                   freq=(/w(2), w(3)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_jk/), int(shape(e3)), e3)
    call print_polarizability(e3, (/w(2), w(3)/), j, k)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_j(j),             &
                   D_l(l)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_j(j),             &
                   F_l(l)              &
                  /),                  &
                   D_jl(j, l),         &
                   F_jl(j, l),         &
                   comp=(/j, l/),      &
                   freq=(/w(2), w(4)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_jl/), int(shape(e3)), e3)
    call print_polarizability(e3, (/w(2), w(4)/), j, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_k(k),             &
                   D_l(l)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_k(k),             &
                   F_l(l)              &
                  /),                  &
                   D_kl(k, l),         &
                   F_kl(k, l),         &
                   comp=(/k, l/),      &
                   freq=(/w(3), w(4)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_kl/), int(shape(e3)), e3)
    call print_polarizability(e3, (/w(3), w(4)/), k, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j),             &
                   D_k(k),             &
                   D_ij(i, j),         &
                   D_ik(i, k),         &
                   D_jk(j, k)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j),             &
                   F_k(k),             &
                   F_ij(i, j),         &
                   F_ik(i, k),         &
                   F_jk(j, k)          &
                  /),                  &
                   D_ijk(i, j, k),     &
                   F_ijk(i, j, k),     &
                   comp=(/i, j, k/),   &
                   freq=(/w(1), w(2), w(3)/))
    e4 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ijk/), int(shape(e4)), e4)
    call print_polarizability(e4, (/w(1), w(2), w(3)/), i, j, k)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j),             &
                   D_l(l),             &
                   D_ij(i, j),         &
                   D_il(i, l),         &
                   D_jl(j, l)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j),             &
                   F_l(l),             &
                   F_ij(i, j),         &
                   F_il(i, l),         &
                   F_jl(j, l)          &
                  /),                  &
                   D_ijl(i, j, l),     &
                   F_ijl(i, j, l),     &
                   comp=(/i, j, l/),   &
                   freq=(/w(1), w(2), w(4)/))
    e4 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ijl/), int(shape(e4)), e4)
    call print_polarizability(e4, (/w(1), w(2), w(4)/), i, j, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_k(k),             &
                   D_l(l),             &
                   D_ik(i, k),         &
                   D_il(i, l),         &
                   D_kl(k, l)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_k(k),             &
                   F_l(l),             &
                   F_ik(i, k),         &
                   F_il(i, l),         &
                   F_kl(k, l)          &
                  /),                  &
                   D_ikl(i, k, l),     &
                   F_ikl(i, k, l),     &
                   comp=(/i, k, l/),   &
                   freq=(/w(1), w(3), w(4)/))
    e4 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ikl/), int(shape(e4)), e4)
    call print_polarizability(e4, (/w(1), w(3), w(4)/), i, k, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_j(j),             &
                   D_k(k),             &
                   D_l(l),             &
                   D_jk(j, k),         &
                   D_jl(j, l),         &
                   D_kl(k, l)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_j(j),             &
                   F_k(k),             &
                   F_l(l),             &
                   F_jk(j, k),         &
                   F_jl(j, l),         &
                   F_kl(k, l)          &
                  /),                  &
                   D_jkl(j, k, l),     &
                   F_jkl(j, k, l),     &
                   comp=(/j, k, l/),   &
                   freq=(/w(2), w(3), w(4)/))
    e4 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_jkl/), int(shape(e4)), e4)
    call print_polarizability(e4, (/w(2), w(3), w(4)/), j, k, l)

    call pert_dens(decomp,                   &
                   S,                        &
                 (/'EL', 'EL', 'EL', 'EL'/), &
                 (/1, 1, 1, 1/),             &
                 (/                          &
                   D,                        &
                   D_i(i),                   &        
                   D_j(j),                   &        
                   D_k(k),                   &        
                   D_l(l),                   &        
                   D_ij(i, j),               &           
                   D_ik(i, k),               &           
                   D_il(i, l),               &           
                   D_jk(j, k),               &           
                   D_jl(j, l),               &           
                   D_kl(k, l),               &
                   D_ijk(i, j, k),           &   
                   D_ijl(i, j, l),           &   
                   D_ikl(i, k, l),           &   
                   D_jkl(j, k, l)            &   
                  /),                        &
                 (/                          &
                   F,                        &
                   F_i(i),                   &        
                   F_j(j),                   &        
                   F_k(k),                   &        
                   F_l(l),                   &        
                   F_ij(i, j),               &           
                   F_ik(i, k),               &           
                   F_il(i, l),               &           
                   F_jk(j, k),               &           
                   F_jl(j, l),               &           
                   F_kl(k, l),               &
                   F_ijk(i, j, k),           &   
                   F_ijl(i, j, l),           &   
                   F_ikl(i, k, l),           &   
                   F_jkl(j, k, l)            &   
                  /),                        &
                   D_ijkl(i, j, k, l),       &
                   F_ijkl(i, j, k, l),       &
                   comp=(/i, j, k, l/),      &
                   freq=(/w(1), w(2), w(3), w(4)/))
    e5 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ijkl/), int(shape(e5)), e5)
    call print_polarizability(e5, (/w(1), w(2), w(3), w(4)/), i, j, k, l)

    delta_element = -e5(h, i, j, k, l)

    D_i    = 0
    D_j    = 0
    D_k    = 0
    D_l    = 0
    D_ij   = 0
    D_ik   = 0
    D_il   = 0
    D_jk   = 0
    D_jl   = 0
    D_kl   = 0
    D_ijk  = 0
    D_ijl  = 0
    D_ikl  = 0
    D_jkl  = 0
    D_ijkl = 0

    F_i    = 0
    F_j    = 0
    F_k    = 0
    F_l    = 0
    F_ij   = 0
    F_ik   = 0
    F_il   = 0
    F_jk   = 0
    F_jl   = 0
    F_kl   = 0
    F_ijk  = 0
    F_ijl  = 0
    F_ikl  = 0
    F_jkl  = 0
    F_ijkl = 0

  end function
#endif

end module
