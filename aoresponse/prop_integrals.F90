! Copyright 2009 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

module prop_integrals
! routines for property contractions and matrices

   use matrix_defop
   use character_processing
#ifndef PRG_DIRAC
   use dal_interface, only:              &
           di_read_operator_int,         &
           di_GET_GbDs,                  &
           di_get_sigma_xc_cont,         &
           di_get_T_xc_cont,             &
           di_get_MagDeriv_FxD_DFT,      &
           di_get_MagDeriv_GxD_DFT,      &
           di_get_geomDeriv_molgrad_DFT, &
           di_get_geomDeriv_FxD_DFT,     &
           di_get_geomDeriv_GxD_DFT,     &
           di_select_wrk,                &
           di_deselect_wrk
                   
   use matrix_operations, only: mat_to_full, mat_set_from_full
#else
    use dirac_interface
    use dft_cfg
    use dft_temp
#endif /* PRG_DIRAC */

   implicit none
   public prop_energy
   public prop_fock
   public prop_overlap
   public prop_reort
   public prop_shape
   public prop_imag
   public prop_check
   public prop_integrals_debug

   private

   type ptb !private: 'perturbation'
      character*4  :: code !four-letter abbreviation
      character*64 :: name !long name
      integer      :: rank !number of components (when known, 0 otw.)
      logical      :: imag !imaginary (1st ord.) perturbed integrals
      logical      :: bas  !basis dependent (sa. GEO and MF)
      logical      :: lin  !one-electron operator linear in field strength (EL)
      logical      :: quad !one-electron operator quadratic in field strength (MAGO)
   end type

   !starting with two fields
   logical,parameter:: T=.true.,F=.false.
   type(ptb) :: tab(6) = &                                 !rk im ba ln qu
       (/ ptb('EL'  ,'Electric field'                      , 3, F, F, T, F), &
          ptb('VEL' ,'Velocity'                            , 3, T, F, T, F), &
          ptb('MAGO','Magnetic field w/o. London orbitals' , 3, T, F, F, T), &
          ptb('MAG' ,'Magnetic field with London orbitals' , 3, T, T, F, F), &
          ptb('ELGR','Electric field gradient'             , 6, F, F, T, F), &
          ptb('GEO' ,'Nuclear coordinates'                 , 0, F, T, F, F) /) !0=uninitialized

   integer:: nbas=0, natom=0 !unitialized

   character(1), parameter :: xyz(3) = (/'X','Y','Z'/)

#ifndef PRG_DIRAC
  logical, external :: do_dft
#endif /* PRG_DIRAC */

   logical :: prop_integrals_debug = .false.

contains


   subroutine prop_energy(p, D, e, c, b, o)
   ! Front for energy_core below, filling in optionals, checking sizes, and making permutations
   ! NB: The result of the contraction is ADDED to 'e' - remember to zero before first call
      character(*), intent(in)    :: p(:) !p(n) perturbation lables
      type(matrix), intent(in)    :: D(:) !density matrix expansion to contract integrals with
      complex(8),   intent(inout) :: e(*) !result is added to e(*).
                   !(*) to permit re-ranking. size(e) is not checked
      integer, optional, intent(in) :: &
               c(:), & !c(n) component indices, or c(n+n) components ranges of each p
               b(:), & !b(m) ranks of D(:) expansion: size(D)=product(b+1)
               o(:) !o(m) indices in e of each rank in D; default 1,2..m
      integer :: n, m, ne, op(size(p)), cc(size(p)*2), i
      logical :: zerg,imag
      !process and check arguments
      call prop_check('prop_energy', D(1), p, c=c, nd=size(D), b=b, o=o,    &
                      np=n, op=op, cc=cc, m=m, ne=ne, zerg=zerg, imag=imag)
      if (zerg .or. m==0.and.imag) return !the result is zero
      call subr() !since bb(m) and oe(ne) and ee(ne)
   contains
      subroutine subr()
         integer    :: bb(m), oe(ne)
         complex(8) :: ee(ne)
         call prop_order(n, op, cc, m, size(D), ne, bb, oe, b=b, o=o)
         !call core function erg
         ee = 0
         call energy_core(n, m, (/ (p(op(i)), i=1,n) /),                 &
                          (/ (cc(op(i)), i=1,n), (cc(n+op(i)), i=1,n) /), &
                          D, bb, ee)
         !add the reordered output to e
         e(1:ne) = e(1:ne) + (/ (ee(oe(i)), i=1,ne) /)
      end subroutine
   end subroutine


   subroutine prop_fock(p, D, F, c, b, o)
   ! Front for fock_core below, checking dimensions and making permutations
   ! Note that the result of the contraction is ADDED to F, ie. remeber to zero F before calling
      character(*), intent(in)    :: p(:) !p(n) perturbation lables
      type(matrix), intent(in)    :: D(:) !density matrix expansion to contract integrals with
      type(matrix), intent(inout) :: F(*) !result is added to F(*).
                  !F(*) to permit re-ranking, so size not checked
      integer, optional, intent(in) :: &
               c(:), & !c(n) component indices, or c(n+n) components ranges of each p
               b(:), & !b(m) ranks of D(:) expansion: size(D)=product(b+1)
               o(:)    !o(m) indices in e of each rank in D; default 1,2..m
      integer :: n, m, nf, op(size(p)), cc(2*size(p)), i
      logical :: zfck
      call prop_check('prop_fock',D(1),p, c=c,nd=size(D),b=b,o=o &
                     ,np=n,op=op,cc=cc,m=m,ne=nf,zfck=zfck)
      if (zfck) return !the result is zero
      call subr() !since bb(m) and oe(ne) and FF(ne)
   contains
      subroutine subr()
         integer      :: bb(m), of(nf)
         type(matrix) :: FF(nf)
         call prop_order(n, op, cc, m, size(D), nf, bb, of, b=b, o=o)
         !order F into FF, call core function fock_core, then reorder back into F
         do i = 1, nf
            FF(of(i)) = F(i)
            F(i) = 0
         end do
         call fock_core(n, m, (/ (p(op(i)), i=1,n) /),                  &
                        (/ (cc(op(i)), i=1,n), (cc(n+op(i)), i=1,n) /), &
                        D, bb, FF)
         do i = 1, nf
            F(i) = FF(of(i))
            FF(of(i)) = 0
         end do

      end subroutine
   end subroutine


   subroutine prop_overlap(p, S0, S, c, w)
   !Front for 'overlap_core' below, checking dimensions. Argument S0 is for context only
   !S will on exit contain the array of perturbed overlap matrices (which may all be zero (0*S0))
      character(*), intent(in)    :: p(:) !p(n) perturbation lables
      type(matrix), intent(in)    :: S0 !overlap matrix, for context (shape, basis set etc.)
      type(matrix), intent(inout) :: S(*) !resulting perturbed overlap matrices.
                           !intent(inout) is dropped to permit re-ranking: (/S(:,:)/)
      integer, intent(in), optional :: &
               c(:) !c(n) component indices, or c(n+n) components ranges of each p
      complex(8), intent(in), optional :: w(:)    !w(n+1) frequences for each p
      integer    :: n, op(size(p)), cc(2*size(p)), ns, i
      complex(8) :: ww(size(p)+1)
      logical    :: zovl
      !process and check arguments
      call prop_check('prop_overlap', S0, p, c=c, w=w, &
                      np=n, op=op, cc=cc, ww=ww, ne=ns, zovl=zovl)
      !delete any previous S(:)
      do i = 1, ns
         S(i) = 0d0*S0
      end do
      if (zovl) return !the result is zero
      call subr() !since bb(m) and oe(ne) and FF(ne)
   contains
      subroutine subr()
         integer      :: os(ns)
         type(matrix) :: SS(ns)
         call prop_order(n, op, cc, 0, 0, ns, cc(1:0), os)
         !call core function ovl, then reorder back into S
         do i = 1, ns
            SS(i) = 0d0*S0
         end do
         call overlap_core(n, (/ (p(op(i)),i=1,n) /), &
                           (/ (cc(op(i)), i=1,n), (cc(n+op(i)), i=1,n) /), &
                           (/ (ww(op(i)), i=1,n), ww(n+1) /), S0, SS)
         S(1:ns) = (/ (SS(os(i)),i=1,ns) /) !reorder (aliasing)
      end subroutine
   end subroutine


   subroutine prop_reort(p, c, S0, D, DFD, b, e, w)
   !front for 'reort_core' below, checking dimensions and taking care of permutations (not yet, actually)
      character(*) :: p(:) !p(n) perturbation lables
      integer      :: n, c(:) !c(n*2) component range for each p
      integer      :: m, b(:) !b(m) ranks of D and DFD
      type(matrix) :: S0 !uneperturbed overlap matrix
      type(matrix) :: D(*) !perturbed density matrices, size(D)=product(b)
      type(matrix) :: DFD(*) !energy-weighted density matrices, size(DFD)=product(b)
      complex(8), intent(inout) :: e(*) !property contributions are ADDED to e(*)
                                       !size(e)=product(b)*product(c(size(p)+1:)-c(:size(p)))
      complex(8), optional :: w(:) !frequencies w(n+1)
      integer    :: i, j, k, ip(size(p))
      complex(8) :: ww(size(p)+1)
      n = size(p)
      m = size(b)
      if (size(c) /= n*2) &
         stop 'prop_reort: size(c)/=2*size(p)'
      if (present(w)) then
         if (size(w) /= n+1) stop 'prop_reort: wrong number of frequencies, size(w)/=size(p)+1'
         ww = w
      end if
      call init(S0)
   end subroutine


   subroutine energy_core(n, m, p, c, D, b, e)
      integer      :: n, m, c(n+n), b(m)
      character(*) :: p(n) !perturbation lables
      type(matrix) :: D(*) !density matrix expansion to contract over, size(D)=product(b+1)
      complex(8), intent(out) :: e(*) !property contributions are ADDED to e
                                 !size(e)=product(b)*product(c(n+1:n+n)-c(1:n)+1)
      integer      :: i, j, k, l, ii, jj, kk, ll, i1, j1, k1, l1, pb, pb1
      type(matrix) :: A(6), T
      real(8)      :: r(3*natom*3) !maximum scratch needed
      do i = 1, size(A)
          A(i) = 0d0*D(1)
      end do
      pb = product(b)
      pb1 = product(b+1)
      if (n==0) then
         if (m==0) stop 'prop_energy: unperturbed energy requested'
         !The highest order Ds multiply the unperturbed F, which is unknown here,
         !so insist that no highest-order Ds are present
         if (.not.all((/(iszero(D(i)),i=pb1-pb+1,pb1)/))) &
            stop 'prop_energy: unperturbed Fock matrix needed'
         if (m==2) then !stopping at m=2, for now
            call init_mat(A(1), A(1)) !de_GET_GbDs expects initialized
            if (do_dft()) call init_mat(A(2), A(2)) !so does di_get_sigma..
            do j = 0, b(2)-1
               call di_GET_GbDs(D(1+b(1)+1+j), A(1))

               if (do_dft()) then
!radovan please leave here for debug:
!  call xc_contribution(2, D(1), (/D(1+b(1)+1+j)/), (/A(1)/))
                  call di_get_sigma_xc_cont(D(1+b(1)+1+j), D(1), D(1), A(2))
                  A(1) = A(1) + A(2)

               end if
               do i = 1, b(1)
                  e(i+b(1)*j) = tr(A(1), D(1+i)) !factor 2 hidden in tr
               end do
            end do
         else if (m>2) then
            call quit('prop_energy: m>2 not implemented at the moment')
         end if
      else if (n==1 .and. p(1)=='EL') then
         if (m==0) call DIPNUC_ifc(natom,r(1:3),(/(0d0,i=1,9*natom)/))
         if (m==0) e(1:c(2)-c(1)+1) = e(1:c(2)-c(1)+1) - r(c(1):c(2)) !sign change
         do j = 0, c(2)-c(1)
            jj = c(1)+j
            call oneint(xyz(jj) // 'DIPLEN ',A(1))
            do i = 1, pb
               ii = pb1-pb+i
               e(i+pb*j) = e(i+pb*j) + dot(A(1),D(ii)) !factor 2 hidden in dot
            end do
         end do
      else if (n==1 .and. p(1)=='VEL') then
         do j = 0, c(2)-c(1)
            jj = c(1)+j
            call oneint(xyz(jj) // 'DIPVEL ', A(1))
            do i = 1, pb
               ii = pb1-pb+i
               e(i+pb*j) = e(i+pb*j) + tr(A(1), D(ii)) !factor 2 hidden in tr
            end do
         end do
      else if (n==1 .and. p(1)=='MAGO') then
         do j = 0, c(2)-c(1)
            jj = c(1)+j
            call oneint(xyz(jj) // 'ANGMOM ', A(1))
            do i = 1, pb
               ii = pb1-pb+i
#ifndef PRG_DIRAC
               e(i+pb*j) = tr(A(1), D(ii))/2.0d0 !half, factor 2 hidden in tr
#else
               e(i+pb*j) = tr(A(1), D(ii))
#endif
            end do
         end do
      else if (n==1 .and. p(1)=='MAG' .and. m/=0) then !no m==0 because because MAG imag
         !one-electron integral contribution
         do k = 0, c(2)-c(1)
            kk = c(1)+k
            call oneint('dh/dB' // xyz(kk) // '  ', A(1))
            do i = 1, pb
               ii = pb1-pb+i
               e(i+pb*k) = e(i+pb*k) + tr(A(1),D(ii))
            end do
         end do
         !two electron integral contributions
         if (.not.all((/ (iszero(D(i)), i = pb1-pb+1, pb1) /))) then
            call twofck('M', D(1:1), A(1:3))
            if (do_dft()) then
               do i = 1, 3
                  call init_mat(A(3+i), A(3+i))
               end do
               call di_get_MagDeriv_FxD_DFT(A(4:6), D(1))
               do i = 1, 3
                  A(i) = A(i) + A(3+i)
                  A(3+i) = 0
               end do
            end if
            do i = 1, pb
               ii = pb1-pb+i
               do k = 0, c(2)-c(1)
                  kk = c(1)+k
                  e(i+pb*k) = e(i+pb*k) + tr(A(kk), D(ii))
               end do
            end do
         end if
         if (m==1) then
             !nothing more
         else if (m==2) then
            do i = 1, b(1)
               if (iszero(D(1+i))) cycle
               call twofck('M', D(1+i:1+i), A(1:3))
               if (do_dft()) then
                  do k = 4, 6 !since di_get_MagDeriv expects initialized
                     if (iszero(A(k))) call init_mat(A(k), A(k))
                  end do
                  call di_get_MagDeriv_GxD_DFT(D(1), D(1+i), A(4:6))
                  do k = 1, 3
                     A(k) = A(k) - A(3+k)
                  end do !sign on KS contrib
               end if
               do j = 0, b(2)-1
                  jj = 1+b(1)+1+j
                  do k = 0, c(2)-c(1)
                     kk = c(1)+k
                     e(i+b(1)*j+pb*k) = e(i+b(1)*j+pb*k) + tr(A(kk), D(jj))
                  end do
               end do
            end do
         else
            call quit('prop_energy: MAG, m>2 not implemented yet!')
         end if
      else if (n==1 .and. p(1)=='ELGR') then
         if (m==0) call QDRNUC_ifc(r(1:6))
         if (m==0) e(1:c(2)-c(1)+1) = e(1:c(2)-c(1)+1) - r(c(1):c(2)) !sign change
         do k = 0, c(2)-c(1)
            kk = c(1)+k
            k1 = kk-1
            if (kk>1) k1 = k1+2
            if (kk>3) k1 = k1+1
            call oneint(xyz(1+mod(k1,3))//xyz(1+k1/3)//'THETA ', A(1))
            do i = 1, pb
               ii = pb1-pb+i
               e(i+pb*k) = e(i+pb*k) + tr(A(1), D(ii))
            end do
         end do
      else if (n==1 .and. p(1)=='GEO') then
         if (m==0) then !nuclear repulsion contribution
            call GRADNN_ifc(natom, r(1:3*natom))
            do i = 1, c(2)-c(1)+1
               ii = c(1)+i-1
               e(i) = e(i) + r(ii)
            end do
         end if
         !one-electron integral contribution
         do k = 0, c(2)-c(1)
            kk = c(1)+k

#ifndef PRG_DIRAC
            call oneint('1DHAM'// &
                        prefix_zeros(kk, 3), A(1))
#else
            call oneint('G1N'// &
                        prefix_zeros(kk, 3), A(1))
            T = 1.0d0*A(1)
            call oneint('G1B'// &
                        prefix_zeros(kk, 3), T)
            A(1) = A(1) + T
            call oneint('G1KX'// &
                        prefix_zeros(kk, 3), T)
            call daxpy_iz_jz(1.0d0, T, 1, A(1), 4)
            call oneint('G1KY'// &
                        prefix_zeros(kk, 3), T)
            call daxpy_iz_jz(1.0d0, T, 1, A(1), 3)
            call oneint('G1KZ'// &
                        prefix_zeros(kk, 3), T)
            call daxpy_iz_jz(1.0d0, T, 1, A(1), 2)
#endif /* PRG_DIRAC */

            do i = 1, pb
               ii = pb1-pb+i
               e(i+pb*k) = e(i+pb*k) + tr(A(1), D(ii))
            end do
         end do
         !two-electron integral contribution
         do i = 1, pb
            ii = pb1-pb+i
            call twoctr('G', D(1), D(ii), r(1:3*natom))
            r(1:3*natom) = r(1:3*natom)*8 !factor 8 on these integrals
            if (m == 0) then
!             another factor for molecular gradient
              r(1:3*natom) = 0.5d0*r(1:3*natom)
            end if
            if (do_dft()) then
               if (m==0) call di_get_geomDeriv_molgrad_DFT( r( 3*natom+1 : 3*natom*2 ), &
                                                           natom, D(1))
               if (m/=0) call di_get_geomDeriv_FxD_DFT( r( 3*natom+1 : 3*natom*2 ), &
                                                       natom, D(1), D(ii))
               r(1:3*natom) = r(1:3*natom) + r( 3*natom+1 : 3*natom*2 )
            end if
            do k = 0, c(2)-c(1)
               kk = c(1)+k
               e(i+pb*k) = e(i+pb*k) + r(kk)
            end do
         end do
         if (m==0 .or. m==1) then !nothing more
         else if (m==2) then !integrals over products of first order densities
            do j = 0, b(2)-1
               jj = 1+b(1)+1+j
               do i = 1, b(1)
                  call twoctr('G', D(1+i), D(jj), r(1:3*natom))
                  r(1:3*natom) = r(1:3*natom)*8 !factor 8 on these integrals
                  if (do_dft()) then
                     call di_get_geomDeriv_GxD_DFT( r( 3*natom+1 : 3*natom*2), &
                                                   natom, D(1), D(1+i), D(jj))
                     r(1:3*natom) = r(1:3*natom) + r( 3*natom+1 : 3*natom*2 )
                  end if
                  do k = 0, c(2)-c(1)
                     kk = c(1)+k
                     e(i+b(1)*j+pb*k) = e(i+b(1)*j+pb*k) + r(kk) !factor 8 here
                  end do
               end do
            end do
         else
            call quit('prop_energy: GEO, m>2 not implemented')
         end if
         if (do_dft()) print* !after all the "...integrated to nn electrons..." prints
      else if (n==2 .and. all(p==(/'MAGO','MAGO'/))) then

#ifdef PRG_DIRAC
         if (.not. diamagnetic_via_pn) then
#endif
         do k = 0, c(4)-c(2)
            kk = c(2)+k
            do j = 0, c(3)-c(1)
               jj = c(1)+j
               call oneint(xyz(min(jj,kk)) // xyz(max(jj,kk)) &
                           // 'SUSCGO', A(1))
               do i = 1, pb
                  ii = pb1-pb+i
                  e(i+pb*j+pb*3*k) = tr(A(1), D(ii)) !factor 2 hidden in tr
               end do
            end do
         end do
#ifdef PRG_DIRAC
         endif
#endif

      else if (n==2 .and. all(p==(/'EL ','MAG'/)) .and. m/=0) then !no m==0 because MAG imag
         do l = 0, c(4)-c(2)
            ll = c(2)+l
            do k = 0, c(3)-c(1)
               kk = c(1)+k
               call oneint(xyz(kk) // '-CM1 ' // xyz(ll) // ' ', A(1))
               kk = k+l*(c(3)-c(1)+1)
               do i = 1, pb
                  ii = pb1-pb+i
                  e(i+pb*kk) = e(i+pb*kk) + tr(A(1), D(ii))
               end do
            end do
         end do
      else if (n==2 .and. all(p==(/'MAG','MAG'/))) then
         !one-electron integral contribution
         do l = 0, c(4)-c(2)
            ll = c(2)+l
            do k = 0, c(3)-c(1)
               kk = c(1)+k
               call oneint(xyz(min(kk,ll)) // xyz(max(kk,ll)) &
                           // 'dh/dB2', A(1))
               kk = k+l*(c(3)-c(1)+1)
               do i = 1, pb
                  ii = pb1-pb+i
                  e(i+pb*kk) = e(i+pb*kk) + tr(A(1), D(ii))
               end do
            end do
         end do
         if (m==0 .or. .not.all((/ (iszero(D(i)), i = pb1-pb+1, pb1) /))) then
            call twofck('MM', D(1:1), A(1:6))
            do i = 1, pb
               ii = pb1-pb+i
               do l = 0, c(4)-c(2)
                  ll = c(2)+l
                  do k = 0, c(3)-c(1)
                     kk = c(1)+k
                     kk = min(kk,ll) + max(kk,ll) * (max(kk,ll)-1) / 2
                     k1 = k+(c(3)-c(1)+1)*l
                     e(i+pb*k1) = e(i+pb*k1) + tr(A(kk), D(ii))/merge(2,1,m==0)
                  end do
               end do
            end do
         end if
         if (m==0 .or. m==1) then
             !nothing more
         else if (m==2) then
            do i = 1, b(1)
               if (iszero(D(1+i))) cycle
               call twofck('MM', D(1+i:1+i), A(1:6))
               do j = 0, b(2)-1
                  jj = 1+b(1)+1+j
                  do l = 0, c(4)-c(2)
                     ll = c(2)+l
                     do k = 0, c(3)-c(1)
                        kk = c(1)+k
                        kk = min(kk,ll) + max(kk,ll) * (max(kk,ll)-1) / 2
                        k1 = k+(c(4)-c(2)+1)*l
                        e(i+b(1)*j+pb*k1) = e(i+b(1)*j+pb*k1) + tr(A(kk), D(jj))
                     end do
                  end do
               end do
            end do
         else
            call quit('prop_energy: MAG-MAG, m>2 not implemented yet!')
         end if
      else if (n==2 .and. all(p==(/'MAG ','ELGR'/))) then
         do k = 0, c(4)-c(2) !loop over ELGR indices
            kk = c(2)+k
            k1 = merge(kk+3, merge(kk+2, kk, kk>1), kk>3) - 1
            do j = 0, c(3)-c(1)
               jj = c(1)+j !loop over MAG indices
               call oneint(xyz(1+mod(k1,3)) // xyz(1+k1/3) &
                           // '-QDB ' // xyz(jj), A(1))
               kk = j+k*(c(3)-c(1)+1)
               do i = 1, pb
                  ii = pb1-pb+i !loop over highest-order Ds
                  e(i+pb*kk) = e(i+pb*kk) + tr(A(1), D(ii))
               end do
            end do
         end do
         if (c(2)==1 .and. c(4)==6) then !remove trace
            i1 = pb*(c(3)-c(1)+1)
            do i = 1, i1
               r(1) = e(i+i1*0)+e(i+i1*2)+e(i+i1*5) !trace
               e(i+i1*0) = e(i+i1*0) - r(1)/3
               e(i+i1*2) = e(i+i1*2) - r(1)/3
               e(i+i1*5) = e(i+i1*5) - r(1)/3
            end do
         end if
      else if (n==2 .and. all(p==(/'EL ','GEO'/))) then
         if (m==0) then !nuclear repulsion contribution
            call DIPNUC_ifc(natom, (/0d0,0d0,0d0/), r(1:3*3*natom))
            do j = 0, c(4)-c(2)
               jj = c(2)+j
               do i = 1, c(3)-c(1)+1
                  ii = c(1)+i-1
                  i1 = i+(c(3)-c(1)+1)*j
                  e(i1) = e(i1) - r(ii+3*(jj-1))
               end do
            end do
         end if
         !one-electron integral contribution
         do k = 0, c(4)-c(2)
            kk = c(2)+k
            do j = 0, c(3)-c(1)
               jj = c(1)+j
               call oneint(prefix_zeros(kk, 2) &
                           // ' DPG ' // xyz(jj), A(1))
               do i = 1, pb
                  ii = pb1-pb+i
                  i1 = i+pb*(j+(c(3)-c(1)+1)*k)
                  e(i1) = e(i1) - tr(A(1), D(ii)) !factor -1
               end do
            end do
         end do
      else if (n==2 .and. all(p==(/'MAGO','GEO '/)) .and. m/=0) then
         do k = 0, c(4)-c(2)
            kk = c(2)+k
            do j = 0, c(3)-c(1)
               jj = c(1)+j
               call oneint(prefix_zeros(kk, 3) &
                           // 'AMDR' // xyz(jj), A(1))
               do i = 1, pb
                  ii = pb1-pb+i
                  i1 = i+pb*(j+(c(3)-c(1)+1)*k)
                  e(i1) = e(i1) + 0.5d0*tr(A(1), D(ii))
               end do
            end do
         end do
      else if (n==2 .and. all(p==(/'ELGR','GEO '/))) then
         if (m==0) then !nuclear repulsion
             call quit('prop_energy: ELGR,GEO - nucrep conrib not implemented')
         end if
         !one-electron integral contribution
         do k = 0, c(4)-c(2)
            kk = c(2)+k
            do j = 0, c(3)-c(1)
               jj = c(1)+j
               j1 = jj-1
               if (jj>1) j1 = j1+2
               if (jj>3) j1 = j1+1
               call oneint(prefix_zeros(kk, 2) // 'QDG ' &
                           // xyz(1+mod(j1,3)) // xyz(1+j1/3), A(1))
               do i = 1, pb
                  ii = pb1-pb+i
                  i1 = i+pb*(j+(c(3)-c(1)+1)*k)
                  e(i1) = e(i1) - 3/2d0*tr(A(1), D(ii)) !factor -3/2
               end do
            end do
            if (c(1)==1 .and. c(3)==6) then !remove trace
               do i = 1, pb
                  i1 = i+pb*(c(3)-c(1)+1)*k
                  r(1) = e(i1+pb*0) + e(i1+pb*2) + e(i1+pb*5) !trace
                  e(i1+pb*0) = e(i1+pb*0) - r(1)/3
                  e(i1+pb*2) = e(i1+pb*2) - r(1)/3
                  e(i1+pb*5) = e(i1+pb*5) - r(1)/3
               end do
            end if
         end do
      else if (n==3 .and. all(p==(/'EL ','MAG','MAG'/)) .and. (m==0.or.m==1)) then
         do l = 0, c(6)-c(3)
            ll = c(3)+l
            do k = 0, c(5)-c(2)
               kk = c(2)+k
               do j = 0, c(4)-c(1)
                  jj = c(1)+j
                  call oneint(xyz(jj) // '-CM2' // xyz(min(kk,ll)) &
                              // xyz(max(kk,ll)) // ' ', A(1))
                  jj = j+(c(4)-c(1)+1)*(k+(c(5)-c(2)+1)*l)
                  do i = 1, pb
                     ii = pb1-pb+i
                     e(i+pb*jj) = e(i+pb*jj) + tr(A(1), D(ii))
                  end do
               end do
            end do
         end do
      else
         print *,'prop_energy: no integrals for these perturbations: ', &
                 (p(i)//' ', i = 1, n)
         call quit('prop_energy: no such integrals')
      end if
      A(:) = 0 !delete scratch matrices
   end subroutine


   subroutine fock_core(n, m, p, c, D, b, F)

#ifdef PRG_DIRAC
      use dft_cfg
      use dft_main
#include "dgroup.h"
#include "dcbbas.h"
#endif /* PRG_DIRAC */

      integer,      intent(in)    :: n,m,c(n+n),b(m)
      character(*), intent(in)    :: p(n) !perturbation lables
#ifdef PRG_DIRAC
      type(matrix):: D(product(b+1))
#else
      type(matrix), intent(in)    :: D(*) !density matrix expansion to contract over,
                                          !size(D)=product(b+1)
#endif
      type(matrix), intent(inout) :: F(*) !perturbed Fock matrix contributions are ADDED to F(*)
                                          !size(F)=product(b)*product(c(n+1:n+n)-c(1:n)+1)
      integer      :: i, j, k, ii, jj, kk, i1, pb, pb1
      type(matrix) :: A(9)

#ifdef PRG_DIRAC
      real(8) :: z(ntbas(0)*ntbas(0)*nz)
      integer :: ones(50)     
      integer :: nr_dmat
      type(dft_options_t) :: dft_options
      ones = 1
#endif /* PRG_DIRAC */

      do i = 1, size(A)
         A(i) = 0d0*D(1)
      end do
      pb = product(b)
      pb1 = product(b+1)
      if (n==0) then
         if (m==0) stop 'prop_fock: unperturbed Fock matrix requested'
         call init_mat(A(1),A(1)) !di_GET_GbDs expects initalized
         do i = 1, pb
            ii = pb1-pb+i
            if (iszero(D(ii))) cycle

#ifdef PRG_DIRAC
!radovan: these many dftexc calls will be simplified soon to only one
            call di_GET_GbDs(D(ii), A(1))

            F(i) = F(i) + A(1)

          end do

            if (do_dft()) then !DFT contribution

            z       = 0.0d0
            nr_dmat = size(D)

            do i = 1, nr_dmat
              if (iszero(D(i))) then
                D(i) = tiny(0.0d0)*D(1)
              end if
            end do

            call reset_dft_options(dft_options)
            dft_options%aor_xc_order = m + 1

            if (m == 1) then
              call dftexc(2.0d0*D(1)%elms, &
                          (/               &
                            D(2)%elms      &
                           /),             &
                          z,               &
                          1,               &
                          1,               &
                          (/ones(1:1)/),   &
                          (/ones(1:1)/),   &
                          dft_options,     &
                          0)
            end if

            if (m == 2) then
              call dftexc(2.0d0*D(1)%elms, &
                          (/               &
                            D(2)%elms,     &
                            D(3)%elms,     &
                            D(4)%elms      &
                           /),             &
                          z,               &
                          3,               &
                          1,               &
                          (/ones(1:3)/),   &
                          (/ones(1:3)/),   &
                          dft_options,     &
                          0)
            end if

            if (m == 3) then
!             restriction: 
!             all perturbation operators are same
!             so i pass only the representative density matrices
              call dftexc(2.0d0*D(1)%elms,   &
                          (/                 &
                            D(3+1)%elms,     &
                            D(6+1)%elms      &
                           /),               &
                          z,                 &
                          m - 1,             &
                          1,                 &
                          (/ones(1:m - 1)/), &
                          (/ones(1:m - 1)/), &
                          dft_options,       &
                          0)
            end if

            if (m == 4) then
!             restriction: 
!             all perturbation operators are same
!             so i pass only the representative density matrices
              call dftexc(2.0d0*D(1)%elms,   &
                          (/                 &
                            D( 4+1)%elms,    &
                            D(10+1)%elms,    &
                            D(14+1)%elms     &
                           /),               &
                          z,                 &
                          m - 1,             &
                          1,                 &
                          (/ones(1:m - 1)/), &
                          (/ones(1:m - 1)/), &
                          dft_options,       &
                          0)
            end if

            if (m == 5) then
!             restriction: 
!             all perturbation operators are same
!             so i pass only the representative density matrices
              call dftexc(2.0d0*D(1)%elms,   &
                          (/                 &
                            D( 5+1)%elms,    &
                            D(15+1)%elms,    &
                            D(25+1)%elms,    &
                            D(30+1)%elms     &
                           /),               &
                          z,                 &
                          m - 1,             &
                          1,                 &
                          (/ones(1:m - 1)/), &
                          (/ones(1:m - 1)/), &
                          dft_options,       &
                          0)
            end if

            if (m == 6) then
!             restriction: 
!             all perturbation operators are same
!             so i pass only the representative density matrices
              call dftexc(2.0d0*D(1)%elms,   &
                          (/                 &
                            D( 6+1)%elms,    &
                            D(21+1)%elms,    &
                            D(41+1)%elms,    &
                            D(56+1)%elms,    &
                            D(62+1)%elms     &
                           /),               &
                          z,                 &
                          m - 1,             &
                          1,                 &
                          (/ones(1:m - 1)/), &
                          (/ones(1:m - 1)/), &
                          dft_options,       &
                          0)
            end if

            if (m == 7) then
!             restriction: 
!             all perturbation operators are same
!             so i pass only the representative density matrices
              call dftexc(2.0d0*D(1)%elms,   &
                          (/                 &
                            D( 7+1)%elms,    &
                            D(28+1)%elms,    &
                            D(63+1)%elms,    &
                            D(98+1)%elms,    &
                            D(119+1)%elms,   &
                            D(126+1)%elms    &
                           /),               &
                          z,                 &
                          m - 1,             &
                          1,                 &
                          (/ones(1:m - 1)/), &
                          (/ones(1:m - 1)/), &
                          dft_options,       &
                          0)
            end if

              A(1) = tiny(0.0d0)*D(1)
              call daxpy(ntbas(0)*ntbas(0)*nz, 4.0d0, z, 1, A(1)%elms, 1)
              F(1) = F(1) + A(1)

            end if




#else /* ifdef PRG_DIRAC */
            call twofck('', D(ii:ii), A(1:1))
            F(i) = F(i) + A(1)


            if (do_dft()) then !DFT contribution
               call di_get_sigma_xc_cont(D(ii), D(1), D(1), A(1))
               F(i) = F(i) + A(1)
            end if
         end do
         if (do_dft().and.m==1) then
            !nothing more
         else if (do_dft().and.m==2) then
            do j = 0, b(2)-1
               jj = 1+b(1)+1+j
               if (iszero(D(jj))) cycle
               do i = 1, b(1)
                  ii = 1+i
                  if (iszero(D(ii))) cycle
                  call di_get_T_xc_cont(D(1), D(ii), D(jj), A(1))
                  F(i+b(1)*j) = F(i+b(1)*j) + A(1)
               end do
            end do
         else if (do_dft()) then
            call quit('prop_fock: DFT m==3 not omplemented')
         end if


#endif /* ifdef PRG_DIRAC */


      else if (n==1 .and. p(1)=='EL') then !only m==0
         do i = 1, c(2)-c(1)+1
            ii = c(1)+i-1
            call oneint(xyz(ii) // 'DIPLEN ', A(1))
            F(i) = F(i) + A(1)
         end do
      else if (n==1 .and. p(1)=='MAGO' .and. m==0) then !only m==0
         do i = 1, c(2)-c(1)+1
            ii = c(1)+i-1
            call oneint(xyz(ii) // 'ANGMOM ', A(1))
#ifndef PRG_DIRAC
            F(i) = F(i) + 0.5d0*A(1)
#else
            F(i) = F(i) +       A(1)
#endif
         end do
      else if (n==1 .and. p(1)=='MAG') then
         if (m==0) then
            do i = 1, c(2)-c(1)+1
               ii = c(1)+i-1
               call oneint('dh/dB' // xyz(ii) // '  ', A(1))
               F(i) = F(i) + A(1)
            end do
         end if
         do i = 1, pb
            ii = pb1-pb+i
            if (iszero(D(ii))) cycle
            call twofck('M', D(ii:ii), A(1:3))
            do k = 0, c(2)-c(1)
               kk = c(1)+k
               F(i+pb*k) = F(i+pb*k) + A(kk)
            end do
            if (do_dft()) then
               if (m==0) call di_get_MagDeriv_FxD_DFT(A(1:3), D(1))
               if (m/=0) call di_get_MagDeriv_GxD_DFT(D(1), D(ii), A(1:3))
               do k = 0, c(2)-c(1)
                  kk = c(1)+k
                  F(i+pb*k) = F(i+pb*k) + A(kk)
               end do
            end if
         end do
      else if (n==1 .and. p(1)=='ELGR' .and. m==0) then !only m==0
         do i = 1, c(2)-c(1)+1
            ii = c(1)+i-1
            i1 = ii-1
            if (ii>1) i1 = i1+2
            if (ii>3) i1 = i1+1
            call oneint(xyz(1+mod(i1,3)) // xyz(1+i1/3) // 'THETA ', A(1))
            F(i) = F(i) + A(1)
         end do
      else if (n==1 .and. p(1)=='GEO') then
         if (m==0) then
            do i = 1, c(2)-c(1)+1
               ii = c(1)+i-1
#ifndef PRG_DIRAC
               call oneint('1DHAM'// &
                           prefix_zeros(ii, 3), A(1))
#else
               call quit('implement me')
#endif
               F(i) = F(i) + A(1)
            end do
         end if
         do i = 1, 3
            call init_mat(A(i), D(1))
         end do
         do i = 1, pb
            ii = pb1-pb+i
            if (iszero(D(ii))) cycle
            do k = 0, c(2)-c(1)
               kk = c(1)+k
               if (k==0 .or. mod(kk-1,3) == 0) &
                  call twofck('G', D(ii:ii), A(1:3), a=1+(kk-1)/3)
               F(i+pb*k) = F(i+pb*k) + A(1+mod(kk-1,3))
            end do
            if (do_dft()) call quit('prop_fock: GEO with DFT not implemented')
         end do
      else if (n==2 .and. all(p==(/'EL ','MAG'/)) .and. m==0) then
         do j = 0,c(4)-c(2)
            jj = c(2)+j;
            do i = 1, c(3)-c(1)+1
               ii = c(1)+i-1
               kk = i+(c(3)-c(1)+1)*j
               call oneint(xyz(ii) // '-CM1 ' // xyz(jj) // ' ', A(1))
               F(kk) = F(kk) + A(1)
            end do
         end do
      else if (n==2 .and. all(p==(/'MAGO','MAGO'/))) then !only m==0

#ifdef PRG_DIRAC
         if (.not. diamagnetic_via_pn) then
#endif
         do j = 0, c(4)-c(2)
            jj = c(2)+j
            do i = 1, c(3)-c(1)+1
               ii = c(1)+i-1
               kk = i+(c(3)-c(1)+1)*j
               call oneint(xyz(min(ii,jj)) // xyz(max(ii,jj)) // 'SUSCGO', A(1))
               F(kk) = F(kk) + A(1)
            end do
         end do
#ifdef PRG_DIRAC
         endif
#endif

      else if (n==2 .and. all(p==(/'MAG','MAG'/))) then
         !call twofck('MM',D(1),(/(D(1),i=1,9)/))
         call quit('prop_fock: MAG-MAG has a bug and is disabled')
      else
         print *,'prop_fock: no integrals for these perturbations:', &
                    (' ' // p(i), i = 1, n)
         call quit('prop_fock: no such integrals')
      end if
      A(:) = 0
   end subroutine


   subroutine overlap_core(n,p,c,w,S0,S)
   !frequency generalized perturbed overlap matrices
      integer,      intent(in)    :: n, c(n+n)
      character(*), intent(in)    :: p(n) !perturbation lables
      complex(8),   intent(in)    :: w(n+1) !frequencies
      type(matrix), intent(in)    :: S0 !unperturbed overlap matrix, for context
      type(matrix), intent(inout) :: S(*) !result, perturbed overlap matrices
                                         !size(S)=product(c(n+1:n+n)-c(1:n)+1)
      integer      :: i, j, k, ii, jj, kk
      type(matrix) :: A
      A = 0d0*S0 !scratch matrix
      if (n==0) then
         stop 'prop_overlap: unperturbed overlap matrix requested'
      else if (n==1 .and. p(1)=='MAG') then
         do i = 1, c(2)-c(1)+1
            ii = c(1)+i-1
            if (w(1)==0 .and. w(2)==1) then
               call oneint('dS/dB' // xyz(ii) // '  ', S(i))
            else
               call oneint('d<S|/dB' // xyz(ii), S(i))
               call oneint('d|S>/dB' // xyz(ii), A)
               S(i) = (w(2)-w(1)) * S(i) + (w(2)+w(1)) * A
            end if
         end do
      else if (n==1.and.p(1)=='GEO') then
         if (.not.all(w==(/0d0,1d0/))) &
            call quit('PropInt.prop_overlap: unexpected w')
         do i = 1, c(2)-c(1)+1
            ii = c(1)+i-1
#ifndef PRG_DIRAC
            call oneint('1DOVL'// &
                        prefix_zeros(ii, 3), S(i))
#else
            call oneint('G1O'// &
                        prefix_zeros(ii, 3), S(i))
#endif /* PRG_DIRAC */
            S(i) = -S(i) !notice minus sign
         end do
      else if (n==2.and.p(1)=='MAG'.and.p(2)=='MAG') then
         do j = 0, c(4)-c(2)
            jj = c(2)+j !loop over requested second components
            do i = 1, c(3)-c(1)+1
               ii = c(1)+i-1 !loop over first components
               if (w(1)==0 .and. w(2)==0 .and. w(3)==1) then
                  call oneint('dS/dB2' // xyz(min(ii,jj)) // &
                              xyz(max(ii,jj)), S(i+(c(3)-c(1)+1)*j))
               else
                  !only load >> XX XY YY XZ YZ ZZ, use dag for <<
                  call oneint('>>S/B2' // xyz(min(ii,jj)) // xyz(max(ii,jj)), A)
                  S(i+(c(3)-c(1)+1)*j) = (+w(1)+w(2)+w(3)) * A &
                                       + (-w(1)-w(2)+w(3)) * dag(A)
                  !only load <> XX XY YY XZ YZ ZZ, use dag for the opposite
                  call oneint('<>S/B2' // xyz(min(ii,jj)) // xyz(max(ii,jj)), A)
                  S(i+(c(3)-c(1)+1)*j) = S(i+(c(3)-c(1)+1)*j) &
                          + merge( -w(1)+w(2)+w(3), +w(1)-w(2)+w(3), ii<jj) * A &
                          + merge( +w(1)-w(2)+w(3), -w(1)+w(2)+w(3), ii<jj) * dag(A)
               end if
            end do
         end do
      else
         print *,'prop_overlap: no integrals for these perturbations:', &
                 (' ' // p(i), i = 1, n)
         call quit('prop_overlap: no integrals for these perturbations')
      end if
      A = 0 !delete scratch matrix
   end subroutine


   subroutine init(SD)
   !ajt: this should be removed, somehow
      type(matrix) :: SD
      if (SD%nrow==nbas) return !stub
      nbas = SD%nrow
      call get_natoms(natom)
      tab(idx('GEO'))%rank = 3*natom
      !AllPrt(idx('NS'))%rank = depends on which isotopes
   end subroutine

   subroutine oneint(l, A)
   
!     radovan: frontend to routine which reads AO integrals
!              (typically from AOPROPER)
!              with the label l into matrix A
!              please don't use fixed length of l, 
!              rather (*) to be more flexible
    
      character(*), intent(in)    :: l
      type(matrix), intent(inout) :: A
    
      if (iszero(A)) call init_mat(A, A)
      call di_read_operator_int(l, A)
   
   end subroutine


   subroutine twofck(what, D, F, a, b)
   !what='G' : 3 geometry-derivative Fock matrices; atom a
   !what='M' : 3 magnetic-derivative Fock matrices
   !what='GG': 9 2nd geom.-derivative Fock matrices; atoms a,b
   !what='MM': 9 2nd magn.-derivative Fock matrices
   !ajt apr09 Generalized to complex D
      character(*),      intent(in) :: what
      type(matrix),      intent(in) :: D(:)
      type(matrix),   intent(inout) :: F(:)
      integer, optional, intent(in) :: a, b
      real(8),pointer :: wrk(:)
      type(matrix)    :: dupF
      integer         :: i, j, k, nd, nf, n2, lwrk
      nd = size(D)
      nf = size(F)
      n2 = nbas**2
      if (len(what) == 0) then
         if (D(1)%complex) then
            if (.not.iszero(F(1)) .and. .not.F(1)%complex) F(1) = 0
            if (iszero(F(1))) call init_mat(F(1), F(1), complex=.true.)
            call init_mat(dupF, re(F(1)), alias=.true.)
            call di_GET_GbDs(re(D(1)), dupF)
            call init_mat(dupF, im(F(1)), alias=.true.)
            call di_GET_GbDs(im(D(1)), dupF)
         else
            if (iszero(F(1))) call init_mat(F(1), F(1))
            call di_GET_GbDs(D(1), F(1))
         end if
         return
      else if (what == 'M' .and. nf == 3*nd) then
         lwrk = 50*n2 + 10000*nbas + 5000000
      else if (what == 'G' .and. nf == 3*nd) then
         if (.not.present(a)) &
            call quit("prop_integ-twofck: Geometry Fock requires atom 'a' specified")
         lwrk = 50*n2 + 10000*nbas + 5000000
      else if (what == 'MM' .and. nf == 6*nd) then
         lwrk = 50*n2 + 10000*nbas + 5000000
      else if (what == 'GG' .and. nf == 3*3*nd) then
         if (.not.present(a) .or. .not.present(b)) &
            call quit("prop_integ-twock: 2nd ord. geom Fock reqires atoms 'a','b' specified")
         lwrk = (nd+2*nf+50)*n2 + 10000*nbas + 5000000
      else
         call quit('prop_integ-twofck: wrong size(F) for what='//what)
      end if
      call di_select_wrk(wrk, lwrk)
      do i = 1, nd
         wrk(1+n2*(i-1) : n2*i) = D(i)%elms
      end do
      !ajt Redundant: wrk(nbas**2+1:nbas**2*(nf+1)) = 0
#ifndef PRG_DIRAC
      call GRCONT(wrk( 1+n2*nd+n2*nf : lwrk ), (lwrk-n2*nd-n2*nf),   &
                  wrk( 1+n2*nd : n2*nd+n2*nf ), n2*nf,               &
                  (what(1:1) == 'G'), (what(1:1) == 'M'), len(what), &
                  0, .false., .true., wrk(1:n2*nd), nd)
#else
      call quit('debug grcont inside twofck')
#endif /* ifndef PRG_DIRAC */
      do i = 1, nf
         if (iszero(F(i))) call init_mat(F(i), F(i))
         F(i)%elms = wrk( 1+n2*(nd+i-1) : n2*(nd+i) )
      end do
      call di_deselect_wrk(wrk, lwrk)
   end subroutine


   subroutine twoctr(what, Da, Db, E)
   !gm=T -> 3*NATOMS geometriske
   !gm=F -> tre magnetiske
      character(*), intent(in)    :: what
      type(matrix), intent(in)    :: Da, Db
      real(8),      intent(inout) :: E(:)
      real(8), pointer :: wrk(:)
      integer          :: lwrk, nr_qu = 1
#ifdef PRG_DIRAC
!radovan: it's not ideal to have explicit matrix dimensions here
!         DIRAC: nbas*nbas*4
      nr_qu = 4 !nr of quaternion units
#endif /* ifdef PRG_DIRAC */
      if (nbas==0) stop 'prop_integrals/twoctr: nbas not initialized'
      !ajt MagSus contraction doesn't work, so I disabled this. Use twofck instead
      if (what=='MM' .or. what=='MM') &
         call quit("prop_integ-twoctr: 'M' or 'MM' not yet implemented")
      if     (what=='G' .and. size(e) == 3*natom) then
         lwrk = 50*nr_qu*nbas**2+10000*nr_qu*nbas+5000000
      else if (what=='M' .and. size(e)==3) then
         lwrk = 50*nr_qu*nbas**2+10000*nr_qu*nbas+5000000
      else if (what=='GG' .and. size(e) == 3*3*natom**2) then
         lwrk = 50*nr_qu*nbas**2+10000*nr_qu*nbas+5000000
      else if (what=='MM' .and. size(e)==3*3) then
         lwrk = 50*nr_qu*nbas**2+10000*nr_qu*nbas+5000000
      else
          call quit('prop_integrals/twoctr: wrong size(E) for what=' // what)
      end if

#ifdef PRG_DIRAC
!     try with the entire work array aoresponse got
!     di_select_wrk will allocate more if not enough
      lwrk = len_f77_work
#endif /* ifdef PRG_DIRAC */

      call di_select_wrk(wrk, lwrk)

      e = 0.0d0

      call grcont(wrk, lwrk, &
                  e, size(e), what(1:1)=='G', what(1:1)=='M', len(what), &
                  0, .true., .false., (/Da%elms, Db%elms/), 2)

      call di_deselect_wrk(wrk, lwrk)

   end subroutine


   function idx(p)
      character(*) :: p
      integer      :: idx
      do idx = 1, size(tab)
          if (tab(idx)%code == p) return
      end do
      call quit('Perturbation not found: ' // p)
   end function


   function prop_imag(p)
      character(*) :: p
      logical      :: prop_imag
      prop_imag = tab(idx(p))%imag
   end function


!    function prop_rank(p)
!       character(*) :: p
!       integer      :: prop_rank
!       prop_rank = tab(idx(p))%rank
!       if (prop_rank == 0) stop 'prop_rank: rank not initialized'
!    end function


   function prop_shape(p)
      character(*), intent(in) :: p(:)
      integer                  :: prop_shape(size(p)), i
      prop_shape = (/ ( tab(idx(p(i)))%rank, i=1,size(p) ) /)
   end function


   subroutine prop_check(who, S0, p, c, w, nd, b, o,                  &
                         np, op, cc, ww, m, ne, zerg, zfck, zovl, imag)
   ! Process, p,c,w,nd,b finding the order p(op(k)) which sorts p by idx,
   ! fill in from c the component ranges cc(1:n+n), and the frequencies ww(1:n+1).
   ! Determine m from nd and b
      character(*),intent(in) :: who  !who called, prepended to error messages
      character(*),intent(in) :: p(:) !perturbation lables
      type(matrix),intent(in) :: S0   !or D0, matrix needed as reference to
             !molecule, to be able to tell the size of 'GEO' (3*natom), etc.
      integer,    optional,intent(in)  :: c(:), nd, b(:), o(:)
      complex(8), optional,intent(in)  :: w(:)
      integer,    optional,intent(out) :: np, op(size(p)), cc(2*size(p)), m, ne
      complex(8), optional,intent(out) :: ww(size(p)+1)
      logical,    optional,intent(out) :: zerg, zfck, zovl, imag
      integer :: n, ip(size(p)), i, j, k
      logical :: one,zer
      call init(S0)
      n = size(p)
      if (present(np)) np = n
      ip = (/ (idx(p(i)), i=1,n) /)
      !find out if the perturbed integrals are zero,
      !and if it's a one-electron operator
      j = count((/ (tab(ip(i))%lin , i=1,n) /))
      k = count((/ (tab(ip(i))%quad, i=1,n) /))
      one = (j>0 .or. k>0)
      zer = (j>1 .or. k>2 .or. (j==1 .and. k>=1))
      if (j==0 .and. k==2) then
         j=0
         do i=1,n
            if (tab(i)%quad) then
               if (j/=0) zer = (ip(i)/=ip(j))
               j=i
            end if
         end do
      end if
      if (present(zerg)) zerg = zer
      if (present(zfck)) zfck = zer
      !check whether the perturbed overlap integrals are zero
      if (present(zovl)) &
         zovl = any((/ (.not.tab(ip(i))%bas, i=1,n) /))
      !the integrals are imaginary when
      if (present(imag)) & 
         imag = (mod(count((/ (tab(ip(i))%imag, i=1,n) /)),2) == 1)
      !process optional c to produce cc
      if (present(c)) then
         if (size(c) /= n .and. size(c) /= n+n) &
            call quit(who // ': c must be c(n) or c(n+n), with p(n)')
         if (any((/ (1>c(i) .or. c(i)>tab(ip(i))%rank, i=1,n) /))) &
            call quit(who // ': c must be 1 <= c(1:n) <= p%comp')
         if (present(cc)) cc = (/ c(1:n), c(1:n) /)
         if (size(c) == n+n) then
            if (any((/ (c(i)>c(n+i) .or. c(n+i)>tab(ip(i))%rank, i=1,n) /))) &
               call quit(who // ': c must be c(1:n) <= c(n+1:n+n) <= p%comp')
            if (present(cc)) cc = c
         end if
         if (present(ne) .and. size(c) == n  ) ne = 1
         if (present(ne) .and. size(c) == n+n) ne = product( c(n+1:n+n) - c(1:n) + 1 )
      else !ranges 1:p%comp
         if (present(cc)) cc = (/ (1, i=1,n), (tab(ip(i))%rank, i=1,n) /)
         if (present(ne)) ne = product((/ (tab(ip(i))%rank, i=1,n) /))
      end if
      !process optional w, produce ww
      if (present(w)) then
         if (size(w) /= n .and. size(w) /= n+1) &
            call quit(who // ': w must be w(n) or w(n+1), with p(n)')
         if (present(ww)) ww(1:n+1) = (/w(1:n),(0d0,0d0)/)
         if (present(ww).and.size(w)==n+1) ww(n+1) = w(n+1)
      else if (present(ww)) then
         ww(1:n) = 0
         ww(n+1) = 1
      end if
      !find order of p(i) according to idx
      if (present(op)) then
         do i = 1, n
            op(i:i) = minloc(ip)
            ip(op(i)) = huge(1) !not found again by minloc
         end do
      end if
      !if nd and b, check product(b+1)=nd, set m
      if (present(nd).and.present(m)) then
         if (present(b)) then
            if (product(b+1)/=nd) call quit(who // &
                    ': wrong number of density matrices: product(b+1)/=nd')
            m = size(b)
            if (present(ne)) ne = ne*product(b)
         else
            m = merge(0,1,nd==1)
            if (present(ne)) ne = ne*merge(1,nd-1,nd==1)
         end if
         if (present(o)) then
            if (size(o)/=m) call quit(who // ': o has wrong size')
            if (any(1>o .or. o>m+n)) call quit(who // ': o out of range')
            if (any((/ ( (o(i)==o(j), i=j+1,m), j=1,m-1) /))) &
               call quit(who // ': repeated element in o')
         end if
         !knowing m, Fock contributions may be zero
         if (one .and. present(zfck) .and. m>0) zfck = .true.
      else if (present(b) .or. present(o)) then
         call quit(who // ': if b or o present, so must nd and m be')
      end if
      !is the energy contribution
   end subroutine


   subroutine prop_order(n, op, cc, m, nd, ne, bb, oe, b, o)
      integer,           intent(in)  :: n, op(n), cc(n+n), m, nd, ne
      integer,           intent(out) :: bb(m), oe(ne)
      integer, optional, intent(in)  :: b(:), o(:) !b(m),o(m) will break ifort
      integer :: kk(m+n,3), k(3), opi(n), i, j, l
      !cumulative product of dimensions (pitches)
      bb = nd-1
      if (present(b)) bb = b
      kk(1:m+n,1) = (/bb, (cc(n+op(i)) - cc(op(i)) + 1, i=1,n) /)
      kk(1,2)=1
      do i = 2, m+n
         kk(i,2) = kk(i-1,2) * kk(i-1,1)
      end do
      !reorder dimensions and pitches into output order
      do i = 1, n !inverse permutation
         opi(op(i)) = i
      end do
      kk(m+1:m+n,1) = (/ (kk(m+opi(i),1), i=1,n) /)
      kk(m+1:m+n,2) = (/ (kk(m+opi(i),2), i=1,n) /)
      if (present(o)) then
         kk(1:m,3) = o(1:m)
         do i = m, 1, -1
            k(1:1) = maxloc(kk(1:i,3))
            j = k(1)
            k = kk(j,:)
            kk(j:k(3)-1,:) = kk(j+1:k(3),:)
            kk(k(3),:) = k
         end do
      end if
      !construct oe from dimensions kk(:,1) and pitches kk(:,2)
      oe(1) = 1; j = 1
      do i = 1, m+n
         oe(j+1:j*kk(i,1)) = (/ (oe(1:j)+l*kk(i,2), l=1, kk(i,1)-1) /)
         j = j*kk(i,1)
      end do
   end subroutine


end module
