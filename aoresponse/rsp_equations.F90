! Copyright 2009 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

module rsp_equations
! calculates response densities / perturbed densities,
! using integrals from module prop_integrals
!
! ajt Detailed description of contents to be added here ;-)

   use matrix_defop
   use prop_integrals

   implicit none

   public pert_dens      !perturbed (or response-) densities and Fock-matrices
   public pert_fock      !contract perturbed Fock-matrices
   public pert_scf_eq    !contract perturbed scf equation residuals
   public excit_energies !known excitation energies
   public rsp_equations_debug

   private

   ! array for storing known excitation energies. When a frequency is encountered which
   ! matches a known excitation energy, it is assumed that a residue is wanted, and the
   ! limit is taken (dropping non-singular terms). (ajt Still not implemented)
   real(8),allocatable :: excit_energies(:)

   !turn on or off debugging in this file
   logical :: rsp_equations_debug = .false.

contains

   subroutine pert_fock(p, D, F, c)
      character(*),      intent(in)    :: p(:) !perturbation lables p(n)
      type(matrix),      intent(in)    :: D(:) !un- and perturbed density matrices
      type(matrix),      intent(inout) :: F(*) !output, perturbed Fock matrices
                                            !(*) to permit re-ranking
      integer, optional, intent(in) :: c(:) !range within each p, c(2*n) (lower,upper indices)
      integer :: n, cc(size(p)*2), b(size(p)), i, pb, pb1, ne
      call prop_check('pert_fock', D(1), p, c=c, np=n, cc=cc, ne=ne)
      if (size(D) /= product( cc(n+1:n+n) - cc(1:n) + 2 )) &
          call quit('pert_fock: size(D) should be product 1+rank p')
      b = cc(n+1:n+n) - cc(1:n) + 1
      pb = product(b)
      pb1 = product(b+1)
      !zero the Fock matrices
      do i = 1, pb
         F(i) = 0.0d0*D(1)
      end do
      if (n==0) then
          call quit('pert_fock: unperturbed Fock matrix requested')
      else if (all(p=='EL  ' .or. p=='MAGO')) then
         !lowest order contribution
         if (n==1 .or. n==2 .and. all(p=='MAGO')) &
             call prop_fock(p, D(1:1), F, c=c)
          !highest order contribution
         call prop_fock(p(1:0), D, F, b=b) !p(1:0) is empty
      else if (n==1) then
         call prop_fock(p, (/D(1)/), F)
         call prop_fock(p(1:0), (/ D(1), D(pb1-pb+1:pb1) /), F)
      else if (n==2) then
         call prop_fock(p, (/D(1)/), F)
         call prop_fock(p(2:2), (/ D(1), D(1+1:1+b(1)) /), F)
         call prop_fock(p(1:1), (/ D(1), D( 1+b(1)+1 : 1+b(1)+b(2) ) /), F, o=(/2/) )
         call prop_fock(p(1:0), (/ D(1), D(1+1:1+b(1)), D( 1+b(1)+1 : 1+b(1)+b(2) ), &
                        D( pb1-pb+1 : pb1 ) /), F, b=(/b(1),b(2)/) ) !only works for HF (now also DFT)
      else !ajt: this needs to be done, at some point
         call quit('pert_fock: general case not implemented')
      end if
   end subroutine

   subroutine pert_scf_eq(p, S0, D, F, DSD, FDS, c, w)
   !p'th derivative of idempotency DSD-D and SCF equation FDS-SDF
      character(*),    intent(in) :: p(:) !perturbation lables p(n)
      type(matrix),    intent(in) :: S0, D(:) !overlap, un- and perturbed density matrices
      type(matrix),    intent(in) :: F(:) !perturbed Fock matrices
      type(matrix), intent(inout) :: DSD(*), FDS(*) !(*) to permit re-ranking
      integer,    optional, intent(in) :: c(:) !range within each p, c(2*n) (lower,upper indices)
      complex(8), optional, intent(in) :: w(:) !signed frequencies of each p
      integer      :: n, cc(2*size(p)), b(size(p)), op(size(p)), ne
      integer      :: i, j, k, i0, j0, k0, ij, ik, jk, ij0, ik0, jk0, ijk, ijk0, ji, jj
      complex(8)   :: ww(size(p)+1)
      type(matrix) :: Sp(3), Tp(3), SDF(1), X
      integer      :: imax
      call prop_check('pert_scf_eq', S0, p, c=c, w=w, np=n, op=op, cc=cc, ww=ww, ne=ne)
      ww(n+1) = 0
      if (any( op /= (/ (i,i=1,n) /) )) &
         call quit('PertSCFeq: General case not implemented')
      b = cc(n+1:n+n) - cc(1:n) + 1
      if (size(D) /= product(b+1)) &
         call quit('pert_scf_eq: size(D) /= product(b+1)' &
                   // ' wrong number of density matrices')
      if (size(F) /= product(b+1)) &
         call quit('pert_scf_eq: size(F)/=product(b+1)' &
                   // ' wrong number of Fock matrices')
      !ajt fixme: This should be written optimized an for general orders

!#define AOR_ONLY_INDIVIDUAL_COMPONENTS
#ifdef AOR_ONLY_INDIVIDUAL_COMPONENTS
!radovan: this is useful for high order rsp functions where the combinatorics
!         becomes painful

        if (all(p == 'EL' .or. p == 'MAGO')) then

          X = -sum(ww/2.0d0)*S0

          imax = size(D)

!radovan: this has to be done in this order
!         otherwise gamma is always zero in d2h
!         no idea why
          DSD(1) = D(imax)*S0*D(1) + D(1)*S0*D(imax) - D(imax)
          do i = 2, imax - 1
            j = imax - i + 1
            DSD(1) = DSD(1) + D(j)*S0*D(i)
          end do

          FDS(1) = (F(1) - X)*D(imax)
          SDF(1) =            D(imax)*(F(1) + X)
          do i = 1, imax
            j = imax - i + 1
            FDS(1) = FDS(1) + F(i)*D(j)
            SDF(1) = SDF(1) + D(i)*F(j)
          end do

          FDS(1) = FDS(1)*S0 - S0*SDF(1)

          SDF(1) = 0
          X      = 0
#else

      if (n==1 .and. all( p=='EL' .or. p=='MAGO' )) then
         do i = 1, b(1)
            DSD(i) = D(1+i)*S0*D(1) + D(1)*S0*D(1+i) - D(1+i)
            FDS(i) = F(1+i)*D(1)*S0 - S0*D(1)*F(1+i) &
                   + F(1)*D(1+i)*S0 - S0*D(1+i)*F(1) - ww(1)*S0*D(1+i)*S0
         end do
      else if (n==2 .and. all( p=='EL' .or. p=='MAGO' )) then
         do i0 = 1, b(1)
            do j0 = 1, b(2)
               ij0 = i0+b(1)*(j0-1)
               i = 1+i0
               j = 1+b(1)+j0
               ij = 1+b(1)+b(2)+ij0
               DSD(ij0) = D(ij)*S0*D(1)+D(1)*S0*D(ij)-D(ij)+D(j)*S0*D(i)+D(i)*S0*D(j)
               FDS(ij0) = (F(ij)*D(1)+F(j)*D(i)+F(i)*D(j)+(F(1)-sum(ww/2)*S0)*D(ij))*S0 &
                        - S0*(D(1)*F(ij)+D(i)*F(j)+D(j)*F(i)+D(ij)*(F(1)+sum(ww/2)*S0))
            end do
         end do
      else if (n==3 .and. all( p=='EL' .or. p=='MAGO' )) then
         do i0 = 1, b(1)
            do j0 = 1, b(2)
               do k0 = 1,b(3)
                  ij0 = i0+b(1)*(j0-1)
                  ik0 = i0+b(1)*(k0-1)
                  jk0 = j0+b(2)*(k0-1)
                  ijk0 = i0+b(1)*(j0-1)+b(1)*b(2)*(k0-1)
                  i = 1+i0
                  j = 1+b(1)+j0
                  k = 1+b(1)+b(2)+k0
                  ij = 1+b(1)+b(2)+b(3)+ij0
                  ik = 1+b(1)+b(2)+b(3)+b(1)*b(2)+ik0
                  jk = 1+b(1)+b(2)+b(3)+b(1)*b(2)+b(1)*b(3)+jk0
                  ijk = 1+b(1)+b(2)+b(3)+b(1)*b(2)+b(1)*b(3)+b(2)*b(3)+ijk0
                  DSD(ijk0) = D(ijk)*S0*D(1) + D(1)*S0*D(ijk) - D(ijk) &
                            + D(jk)*S0*D(i) + D(i)*S0*D(jk) &
                            + D(ik)*S0*D(j) + D(j)*S0*D(ik) &
                            + D(ij)*S0*D(k) + D(k)*S0*D(ij)
                  FDS(ijk0) = ( F(ijk)*D(1) + F(jk)*D(i) + F(ik)*D(j)    &
                               + F(ij)*D(k) + F(k)*D(ij) + F(j)*D(ik)    &
                               + F(i)*D(jk) + (F(1)-sum(ww/2)*S0)*D(ijk) ) * S0 &
                            - S0 * ( D(1)*F(ijk) + D(i)*F(jk) + D(j)*F(ik)   &
                                   + D(k)*F(ij) + D(ij)*F(k) + D(ik)*F(j)    &
                                   + D(jk)*F(i) + D(ijk)*(F(1)+sum(ww/2)*S0) )
               end do
            end do
         end do
#endif /* AOR_ONLY_INDIVIDUAL_COMPONENTS */
      elseif (n==1) then
         do i = 1, b(1)
            call prop_overlap(p, S0, Sp, c=(/cc(1)+i-1/) )
            call prop_overlap(p, S0, Tp, c=(/cc(1)+i-1/), w=(/ww(1),(0d0,0d0)/) )
            DSD(i) = D(1+i)*S0*D(1) + D(1)*S0*D(1+i) - D(1+i) + D(1)*Sp(1)*D(1)
            FDS(i) = (F(1+i) - 1/2d0*Tp(1)) * D(1)*S0                 &
                   - S0*D(1) * (F(1+i) - 1/2d0*Tp(1))                 &
                   + (F(1) - (ww(1)/2)*S0) * (D(1+i)*S0 + D(1)*Sp(1)) &
                   - (Sp(1)*D(1) + S0*D(1+i)) * (F(1) + (ww(1)/2)*S0)
            Sp(1) = 0
            Tp(1) = 0
         end do
      elseif (n==2) then
         do j = 0, b(2)-1
            do i = 1, b(1)
               jj = 1+b(1)+j+1
               ij = i+b(1)*j
               ji = 1+b(1)+b(2)+ij
               call prop_overlap(p(1:1), S0, Sp(1), c=(/cc(1)+i-1/) )
               call prop_overlap(p(1:1), S0, Tp(1), c=(/cc(1)+i-1/), w=(/-ww(1)/2,(0d0,0d0)/) )
               call prop_overlap(p(2:2), S0, Sp(2), c=(/cc(2)+j/) )
               call prop_overlap(p(2:2), S0, Tp(2), c=(/cc(2)+j/), w=(/-ww(2)/2,(0d0,0d0)/) )
               call prop_overlap(p(1:2), S0, Sp(3), c=(/cc(1)+i-1,cc(2)+j/) )
               call prop_overlap(p(1:2), S0, Tp(3), c=(/cc(1)+i-1,cc(2)+j/), &
                                 w = (/ -ww(1)/2, -ww(2)/2, (0d0,0d0) /) )
               DSD(ij) = D(1)*Sp(3)*D(1) + D(1+i)*Sp(2)*D(1) + D(1)*Sp(2)*D(1+i) &
                       + (D(1+i)*S0 + D(1)*Sp(1)) * D(jj) &
                       + D(jj) * (S0*D(1+i) + Sp(1)*D(1))
               if (.not.iszero(D(ji))) &
                  DSD(ij) = DSD(ij) + D(ji)*S0*D(1) + D(1)*S0*D(ji) - D(ji)
               FDS(ij) = (F(ji) + Tp(3) - (ww(1)+ww(2))/2 * Sp(3)) * D(1)*S0          &
                       - S0*D(1) * (F(ji) + Tp(3) + (ww(1)+ww(2))/2 * Sp(3))          &
                       + F(1)*(D(1)*Sp(3)) - (Sp(3)*D(1))*F(1)                        &
                       + (F(jj) + Tp(2) - ww(2)/2 * Sp(2)) * (D(1+i)*S0 + D(1)*Sp(1)) &
                       - (S0*D(1+i) + Sp(1)*D(1)) * (F(jj) + Tp(2) + ww(2)/2 * Sp(2)) &
                       + (F(1+i) + Tp(1) - ww(1)/2 * Sp(1)) * D(jj)*S0                &
                       - S0*D(jj) * (F(1+i) + Tp(1) + ww(1)/2 * Sp(1))
               if (.not.iszero(D(ji))) &
                  FDS(ij) = FDS(ij) + (F(1) - (ww(1)+ww(2))/2 * S0) * D(ji)*S0 &
                                     - S0*D(ji) * (F(1) + (ww(1)+ww(2))/2 * S0)
               if (.not.iszero(Sp(2))) &
                  FDS(ij) = FDS(ij) &
                       + ( (F(1+i) + Tp(1) - ww(1)/2 * Sp(1)) * D(1)         &
                         + (F(1) - ww(1) * S0) * D(1+i) ) * Sp(2)            &
                       - Sp(2) * ( D(1) * (F(1+i) + Tp(1) + ww(1)/2 * Sp(1)) &
                                 + D(1+i) * (F(1) + ww(1) *S0) )
               if (.not.iszero(Sp(1))) &
                  FDS(ij) = FDS(ij) &
                       + (F(1)-ww(2)*S0) * D(jj)*Sp(1) &
                       - Sp(1)*D(jj) * (F(1)+ww(2)*S0)
               Sp = 0
               Tp = 0
            end do
         end do
      else
         stop 'PertSCFeq: General case not implemented yet'
      end if
   end subroutine

   subroutine pert_dens(p, S0, D, F, DD, FF, c, w)
   !p'th derivative of idempotency DSD-D and SCF equation FDS-SDF
      character(*),    intent(in) :: p(:) !perturbation lables p(n)
      type(matrix),    intent(in) :: S0, D(:), F(:) !unperturbed overlap matrix,
                                !lower-order density and Fock matrices
      type(matrix), intent(inout) :: DD(*), FF(*) !highest order perturbed density and Fock
                                !(*) is to permit variable ranks - size cannot be checked
      integer,    optional, intent(in) :: c(:) !range within each p, c(2*n) (lower,upper indices)
      complex(8), optional, intent(in) :: w(:) !signed frequency of each p
      integer    :: n, cc(2*size(p)), b(size(p)), op(size(p)), i, j, pb, pb1, ne
      complex(8) :: ww(size(p)+1)
      call prop_check('pert_dens', S0, p, c=c, w=w, np=n, op=op, cc=cc, ww=ww, ne=ne)
      ww(n+1) = 0
      b = cc(n+1:n+n) - cc(1:n) + 1
      pb = product(b)
      pb1 = product(b+1)
      if (size(D) /= pb1-pb) &
          call quit('pert_dens: D must have size product b+1 minus product b')
      if (size(F) /= pb1-pb) &
          call quit('pert_dens: F must have size product b+1 minus product b')
      call subr()
   contains
      subroutine subr() !called once pb, size of DSD and FDS is known
         type(matrix) :: DSD(pb), FDS(pb)
         do i = 1, pb
            DD(i) = 0d0*S0
            FF(i) = 0d0*S0
         end do
         if (all(p=='EL  ' .or. p=='MAGO') .or. n==1 .or. n==2) then
            call pert_scf_eq(p, S0, (/D,DD(1:pb)/), (/F,FF(1:pb)/), DSD, FDS, c=c, w=ww(1:n))
            call pert_fock(p, (/D,DD(1:pb)/), FF(1:pb), c=c)
            do i = 1, pb
               call solve_scf_eq(sum(ww(1:n)), S0, D(1), F(1), DSD(i), FDS(i), DD(i), FF(i))
            end do
            print* !ajt newline after prints from solve_scf_eq
         else
            stop 'pert_dens: General case not not implemented yet'
         end if
      end subroutine
   end subroutine


   subroutine solve_scf_eq(w, S, D, F, DSD, FDS, Dp, Fp)
   !          in:          --------------------
   !          out:                               ------
   ! ajt: should at some point also return the final residuals,
   ! at least FDS, for use in Sellers-type higher-precision formulas
   ! Bin Gao: add the interface for Dalton main branch
   !          use -DOLD_DALTON

#ifdef PRG_DIRAC
      use dirac_interface
      use aor_cfg
#else
#if defined (OLD_DALTON)
      use mosolver
#else
      use newsolver, only: rsp2_init, rsp_solver2
      use complexsolver, only: rsp_complex_init, rsp_complex_solver
#endif /* OLD_DALTON */
      use scf_config, only: cfg_rsp_complex, cfg_rsp_gamma
#endif /* PRG_DIRAC */

      complex(8),   intent(in)    :: w
      type(matrix), intent(in)    :: S, D, F
      type(matrix), intent(inout) :: DSD, FDS, Dp, Fp
      !----------------------------------------------------------------------------
      real(8)        :: nrm, ww(1), gamsav
      type(matrix)   :: DS, RHS(3), SOL(4) !(3) for tot-real-imag parts
#ifdef PRG_DIRAC               
      type(matrix)   :: R, RHSdag
#endif /* PRG_DIRAC */
      character*4    :: nop(0)
      logical, save  :: first = .true., firstcplx = .true.
      logical        :: hasimag

      DS = D*S
      Dp = DSD - DSD*dag(DS) - DS*DSD
      RHS(1) = FDS*DS - dag(DS)*FDS + Fp
      FDS = 0

      call prop_fock(nop, (/D, Dp/), RHS)

      RHS(1) = -RHS(1)*DS + dag(DS)*RHS(1) !-(..) because RHS
      call init_mat(SOL(1), RHS(1))

#ifndef PRG_DIRAC
      nrm = norm(RHS(1))
      print *, 'before response solver: norm(RHS) = ', nrm

      if (nrm < 1d-10) then
          print *, '=> skipping this equation'
          SOL(1) = 0d0*SOL(1)
      else if (RHS(1)%complex .or. cfg_rsp_complex .or. imag(w) /= 0) then
          hasimag = RHS(1)%complex
          if (.not.hasimag) then
              call init_mat(RHS(2), RHS(1), alias=.true.)
              call init_mat(RHS(1), reset=.true.)
              RHS(1) = (1d0 + tiny(1d0)*(0d0,1d0)) * RHS(2)
              RHS(2) = 0 !free
              SOL(1) = 0 !free
              call init_mat(SOL(1), RHS(1))
          end if
          !create aliases for the real and imaginary parts of RHS and SOL
          call init_mat(RHS(2), re(RHS(1)), alias=.true.)
          call init_mat(RHS(3), im(RHS(1)), alias=.true.)
          call init_mat(SOL(2), re(SOL(1)), alias=.true.)
          call init_mat(SOL(3), im(SOL(1)), alias=.true.)
          !on the first run, initialize the solvers
          if (firstcplx) then
              call rsp2_init(2, 2, 2, 1, 2)
              call rsp_complex_init(1, 1, 1, 1, 1)
              firstcplx = .false.
          end if
          !configure the frequency and damping factor
          ww(1) = -real(w)
          gamsav = cfg_rsp_gamma !save, for restore after solver
          if (cfg_rsp_complex) then
              cfg_rsp_gamma = cfg_rsp_gamma - imag(w)
          else
              cfg_rsp_gamma = -imag(w)
          end if
          !the solver currently doesn't seem to work with zero gamma
          if (abs(cfg_rsp_gamma) < 1d-6) &
              cfg_rsp_gamma = -2d-6
          !solve for the real part of the RHS
          call rsp_complex_solver(F, D, S, 1, RHS(2:2), ww(1:1), 1, &
                                  SOL(2:2), SOL(3:3), .false.)
          !clear re and im aliases of SOL
          call init_mat(SOL(2), reset=.true.)
          call init_mat(SOL(3), reset=.true.)
          !imaginary part of RHS must be solved for separately
          if (hasimag) then
              call init_mat(SOL(4), SOL(1))
              call init_mat(SOL(2), re(SOL(4)), alias=.true.)
              call init_mat(SOL(3), im(SOL(4)), alias=.true.)
              call rsp_complex_solver(F, D, S, 1, RHS(3:3), ww(1:1), 1, &
                                      SOL(2:2), SOL(3:3), .false.)
              SOL(1) = SOL(1) + (0d0,1d0)*SOL(4)
              call init_mat(SOL(2), reset=.true.)
              call init_mat(SOL(3), reset=.true.)
              SOL(4) = 0
          end if
          !clear re and im aliases of RHS
          call init_mat(RHS(2), reset=.true.)
          call init_mat(RHS(3), reset=.true.)
          !restore gamma to avoid breaking other code
          cfg_rsp_gamma = gamsav
      else
          if (first) then
              call rsp2_init(1, 1, 1, 1, 1)
              first = .false.
          end if
          ww(1) = -w !rsp_solver2 seems to solve for negative frequency
          call rsp_solver2(D, S, F, .true., 1, RHS(1:1), ww(1:1), SOL(1:1))
      end if

      Dp = Dp + 2d0*(SOL(1)*dag(DS) - DS*SOL(1))
#else

      RHSdag = dag(RHS(1))

!     isolate h+ part
      R = RHS(1) - RHSdag
      R = -0.5d0*R
      if (rsp_equations_debug) then
        call print_mat(R, label = 'debug solve_scf_eq: R h+')
      end if
!     consider only if significant
      if (real(dot(R, R)) > aor_cfg_threshold_rhs) then
        R%tr_sym = 1
        call response_solver(real(w), R, SOL(1))
        Dp = Dp + SOL(1)
      end if
      R = 0

!     isolate h- part
      R = RHS(1) + RHSdag
      R = -0.5d0*R
      if (rsp_equations_debug) then
        call print_mat(R, label = 'debug solve_scf_eq: R h-')
      end if
!     consider only if significant
      if (real(dot(R, R)) > aor_cfg_threshold_rhs) then
        R%tr_sym = -1
        call response_solver(real(w), R, SOL(1))
        Dp = Dp + SOL(1)
      end if
      R = 0

      RHSdag = 0

#endif /* ifndef PRG_DIRAC */

      DSD    = 0
      DS     = 0
      RHS(1) = 0
      SOL(1) = 0

      !add contribution G(Dp) to Fp. Since prop_fock takes
      !array argument, move to RHS(1) before calling, and then back
      call init_mat(RHS(1), Fp, alias=.true.)
      call prop_fock(nop, (/D, Dp/), RHS(1:1))
      call init_mat(Fp, RHS(1), alias=.true.)
      call init_mat(RHS(1), reset=.true.)

   end subroutine

end module
