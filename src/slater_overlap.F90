module slater_overlap

! Overlap integrals between Slater-type orbitals (STOs), independent of any
! external Slater-basis code.
!
! -----------------------------------------------------------------------------
! Definitions
! -----------------------------------------------------------------------------
!
! A normalized primitive STO is taken in the usual form (Slater, Phys. Rev. 36,
! 57 (1930)),
!
!    chi_{nlm}(r) = N_n(zeta) r^{n-1} exp(-zeta r) Z_{lm}(theta,phi),
!    N_n(zeta)    = (2 zeta)^{n+1/2} / sqrt((2n)!),
!
! with Z_{lm} the real spherical harmonics normalized to unity on the sphere.
! The components of a shell are ordered m = -l, ..., 0, ..., +l, which is the
! ordering ROSE uses elsewhere for spherical functions.
!
! -----------------------------------------------------------------------------
! Method
! -----------------------------------------------------------------------------
!
! Two-center integrals over exponential functions are classically evaluated in
! prolate spheroidal (confocal elliptic) coordinates, in which the exponential
! factors separate:
!
!    xi  = (r_A + r_B)/R  in [1,infinity),   eta = (r_A - r_B)/R  in [-1,1],
!    dV  = (R^3/8) (xi^2 - eta^2) dxi deta dphi,
!    zeta_A r_A + zeta_B r_B = a xi + b eta,
!    a = (zeta_A + zeta_B) R/2,   b = (zeta_A - zeta_B) R/2.
!
! See R. S. Mulliken, C. A. Rieke, D. Orloff and H. Orloff, J. Chem. Phys. 17,
! 1248 (1949), and C. C. J. Roothaan, J. Chem. Phys. 19, 1445 (1951), who reduce
! the integral to the auxiliary functions A_k(a) = int_1^inf xi^k e^{-a xi} dxi
! and B_k(b) = int_{-1}^{1} eta^k e^{-b eta} deta.
!
! This implementation uses the same coordinates but evaluates the integral by
! quadrature rather than through A_k and B_k. The reason is numerical: the
! classical reduction produces alternating sums of A_k and B_k whose terms grow
! with k while the result does not, so it loses digits catastrophically for
! nearly equal exponents (b -> 0) and for higher n. That failure is well
! documented; see for instance the stability discussions in E. O. Steinborn and
! E. Filter, Theor. Chim. Acta 38, 273 (1975), and I. I. Guseinov, J. Phys. B 3,
! 1399 (1970). Positive-weight quadrature of the same integrand has no such
! cancellation.
!
! The quadrature is exact, not merely convergent, in two of the three variables:
!
!   xi   After the substitution xi = 1 + t/a the non-exponential part of the
!        integrand is a polynomial in t of degree at most n_a + n_b, so
!        Gauss-Laguerre with ceil((n_a+n_b+1)/2) nodes integrates it exactly.
!
!   phi  At fixed xi and eta the integrand is a trigonometric polynomial of
!        degree at most l_a + l_b in the azimuth about the bond axis, because
!        each solid harmonic is a homogeneous polynomial of its own degree in
!        the Cartesian components. The trapezoidal rule with n_phi > l_a + l_b
!        points is exact for such a function.
!
!   eta  The integrand is a polynomial times exp(-b eta), which is analytic on
!        [-1,1]. Gauss-Legendre converges geometrically; the node count is
!        chosen from the polynomial degree and from |b|.
!
! Because the rule is a tensor product, the order of summation is immaterial:
! summing over phi first leaves exactly the azimuthal integral, which contains
! only even powers of the cylindrical radius and is therefore the polynomial in
! xi and eta that the Laguerre rule integrates exactly.
!
! No rotation of the harmonics to the bond frame is needed. The classical
! treatments rotate so that only equal |m| pairs survive; here the quadrature
! nodes are converted to global Cartesian coordinates and the harmonics are
! evaluated there, so the rotation is performed implicitly and exactly by the
! azimuthal rule. This removes the need for real-harmonic rotation matrices
! (Ivanic and Ruedenberg, J. Phys. Chem. 100, 6342 (1996)) altogether.
!
! The one-center case is analytic and is handled separately.
!
! -----------------------------------------------------------------------------
! Conventions, and how they relate to the rest of ROSE
! -----------------------------------------------------------------------------
!
! The components of a shell are ordered m = -l..l, as in get_spherical_projection
! in import_fragments.F90. The harmonics here were checked against the
! Schlegel-Frisch coefficients that rcsp in that file produces (H. B. Schlegel
! and M. J. Frisch, Int. J. Quantum Chem. 54, 83 (1995), Eq. 15): for every l
! and m and at every test point the two agree in sign and in ordering and differ
! only by the single positive factor
!
!     S_here(l,m) = sqrt((2l+1)!!/(4 pi)) * S_rcsp(l,m),
!
! which arises from the unnormalised Cartesian Gaussian convention that rcsp
! carries and has nothing to do with the harmonics themselves. The functions
! returned here are normalised, so the diagonal of an overlap matrix built from
! them is unity.
!
! -----------------------------------------------------------------------------
! Thread safety
! -----------------------------------------------------------------------------
!
! The module caches the last Gauss-Laguerre and Gauss-Legendre rule and a table
! of harmonic normalisations. The table is written once and then only read; the
! quadrature caches are rewritten whenever the node count changes. Calling the
! overlap routines from several threads therefore requires either that the
! caches be made thread private or that sto_overlap_warmup be called first and
! the node counts be held fixed. Nothing else in the module carries state.
!
! -----------------------------------------------------------------------------

  use rose_ao, only: slater_func_info_t

  implicit none
  private

  public :: sto_norm
  public :: sto_primitive_overlap
  public :: sto_shell_overlap
  public :: sto_overlap_matrix
  public :: real_solid_harmonics
  public :: cartesian_monomials
  public :: cart_angular_norm
  public :: sto_ncomponents
  public :: sto_overlap_set_accuracy
  public :: sto_overlap_selftest
  public :: sto_overlap_warmup

  integer, parameter :: wp = 8

  ! Contracted shells are described by slater_func_info_t of rose_ao, the Slater
  ! counterpart of gaussian_func_info_t. Note that its orb_momentum is l+1, as it
  ! is for the Gaussians, so the routines below subtract one on the way in.
  !
  ! Two angular conventions are supported, selected by the optional argument
  ! angular that follows the same code as basis_set_info_t%basis_angular:
  !
  !   angular = 2 (spherical, the default)  the 2l+1 real spherical harmonics,
  !              ordered m = -l..l.  This is the DIRAC and pyscf convention.
  !
  !   angular = 1 (Cartesian)  the (l+1)(l+2)/2 monomials x^kx y^ky z^kz with
  !              kx+ky+kz = l, ordered with kx descending outermost and ky
  !              descending within it, which is the order ROSE already uses for
  !              Cartesian Gaussians in get_spherical_projection and the order
  !              ADF uses for its Slater functions.
  !
  ! The two differ only in the angular factor. Both are homogeneous polynomials
  ! of degree l normalised to unit mean square on the unit sphere, so the radial
  ! factor r^{n-1-l}, the normalisation sto_norm and every node count of the
  ! quadrature are the same for either choice; see cartesian_monomials.

  ! Accuracy controls. The defaults give close to full double precision for
  ! the n and l values that occur in Slater basis sets; sto_overlap_set_accuracy
  ! allows them to be tightened for testing.
  integer  :: extra_eta_nodes = 16        ! added to the estimate from |b|
  integer  :: max_eta_nodes   = 96
  real(wp) :: screen_exponent = 60.0_wp   ! skip a pair when min(za,zb)*R exceeds this

  real(wp), parameter :: pi       = 3.141592653589793238462643_wp
  real(wp), parameter :: same_tol = 1.0e-10_wp  ! centers closer than this are one center

  ! Cached quadrature rules. Node counts repeat constantly over a basis, and
  ! generating a rule costs more than using it, so the last rule of each kind
  ! is kept. Not thread safe as written; see the note in the header.
  ! Normalization factors of the real spherical harmonics, filled on first use.
  integer, parameter    :: lmax_table = 10
  logical               :: harm_ready = .false.
  real(wp)              :: harm_norm(0:lmax_table,0:lmax_table)

  integer               :: cached_nlag = 0, cached_nleg = 0
  real(wp), allocatable :: cached_xlag(:), cached_wlag(:)
  real(wp), allocatable :: cached_xleg(:), cached_wleg(:)

contains

! =============================================================================
!  Public entry points
! =============================================================================

  subroutine sto_overlap_set_accuracy (extra_eta, max_eta, screen)
     integer,  intent(in), optional :: extra_eta, max_eta
     real(wp), intent(in), optional :: screen
     if (present(extra_eta)) extra_eta_nodes = extra_eta
     if (present(max_eta))   max_eta_nodes   = max_eta
     if (present(screen))    screen_exponent = screen
  end subroutine sto_overlap_set_accuracy


  subroutine sto_overlap_warmup (nmax,lmax)
  ! Fill the normalisation table and the quadrature caches ahead of any
  ! threaded region. nmax and lmax are the largest principal and angular
  ! quantum numbers that will occur.
     integer, intent(in) :: nmax, lmax
     integer  :: nc
     real(wp), allocatable :: blk(:,:)
     nc = max(2*lmax+1,(lmax+1)*(lmax+2)/2)
     allocate (blk(nc,nc))
     call sto_primitive_overlap (nmax,lmax,1.0_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                 nmax,lmax,1.0_wp,[0.0_wp,0.0_wp,1.0_wp], &
                                 blk(1:2*lmax+1,1:2*lmax+1), angular=2)
     nc = (lmax+1)*(lmax+2)/2
     call sto_primitive_overlap (nmax,lmax,1.0_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                 nmax,lmax,1.0_wp,[0.0_wp,0.0_wp,1.0_wp], &
                                 blk(1:nc,1:nc), angular=1)
     deallocate (blk)
  end subroutine sto_overlap_warmup


  real(wp) function sto_norm (n,zeta) result (anorm)
  ! Normalization constant of a primitive STO: (2 zeta)^{n+1/2} / sqrt((2n)!).
  ! Formed in logarithms so that it survives large n and large zeta.
     integer,  intent(in) :: n
     real(wp), intent(in) :: zeta
     integer  :: i
     real(wp) :: logfac
     logfac = 0.0_wp
     do i = 2, 2*n
        logfac = logfac + log(real(i,wp))
     end do
     anorm = exp((real(n,wp)+0.5_wp)*log(2.0_wp*zeta) - 0.5_wp*logfac)
  end function sto_norm


  subroutine sto_overlap_matrix (shells_bra,shells_ket,s,angular)
  ! Overlap between two lists of contracted Slater shells. The row and column
  ! ordering follows the shell order, with the components of a shell ordered
  ! as described in the header for the angular convention in force.
     type(slater_func_info_t), intent(in) :: shells_bra(:), shells_ket(:)
     real(wp),          intent(out) :: s(:,:)
     integer, intent(in), optional  :: angular

     real(wp), allocatable :: block(:,:)
     integer :: ish, jsh, ioff, joff, ni, nj, ang

     ang = 2
     if (present(angular)) ang = angular

     s = 0.0_wp
     ioff = 0
     do ish = 1, size(shells_bra)
        ni = sto_ncomponents(shells_bra(ish)%orb_momentum-1,ang)
        joff = 0
        do jsh = 1, size(shells_ket)
           nj = sto_ncomponents(shells_ket(jsh)%orb_momentum-1,ang)
           allocate (block(ni,nj))
           call sto_shell_overlap (shells_bra(ish),shells_ket(jsh),block,ang)
           s(ioff+1:ioff+ni,joff+1:joff+nj) = block
           deallocate (block)
           joff = joff + nj
        end do
        ioff = ioff + ni
     end do
  end subroutine sto_overlap_matrix


  subroutine sto_shell_overlap (sa,sb,block,angular)
  ! Overlap block between two contracted Slater shells, with both index ranges
  ! ordered as described in the header for the angular convention in force.
     type(slater_func_info_t), intent(in) :: sa, sb
     real(wp),                 intent(out) :: block(:,:)
     integer, intent(in), optional         :: angular

     real(wp), allocatable :: prim(:,:)
     integer  :: ip, jp, la, lb, ang

     ang = 2
     if (present(angular)) ang = angular

     block = 0.0_wp
     la = sa%orb_momentum - 1
     lb = sb%orb_momentum - 1
     allocate (prim(sto_ncomponents(la,ang),sto_ncomponents(lb,ang)))
     do ip = 1, sa%n_primitives
        do jp = 1, sb%n_primitives
           call sto_primitive_overlap (sa%n(ip),la,sa%exponent(ip),sa%coord, &
                                       sb%n(jp),lb,sb%exponent(jp),sb%coord, prim, ang)
           block = block + sa%coefficient(ip)*sb%coefficient(jp)*prim
        end do
     end do
     deallocate (prim)
  end subroutine sto_shell_overlap


  subroutine sto_primitive_overlap (na,la,za,ca, nb,lb,zb,cb, block, angular)
  ! Overlap between two normalized primitive Slater shells. On return
  ! block(ia,ib) is the overlap of component ia of the first with component ib
  ! of the second, the components being ordered as described in the header.
     integer,  intent(in)  :: na, la, nb, lb
     real(wp), intent(in)  :: za, zb, ca(3), cb(3)
     real(wp), intent(out) :: block(:,:)
     integer, intent(in), optional :: angular

     real(wp) :: d(3), rab
     integer  :: m, ang

     ang = 2
     if (present(angular)) ang = angular

     block = 0.0_wp
     if (na < la+1 .or. nb < lb+1) then
        error stop 'sto_primitive_overlap: principal quantum number below l+1'
     end if
     if (za <= 0.0_wp .or. zb <= 0.0_wp) then
        error stop 'sto_primitive_overlap: non-positive Slater exponent'
     end if

     d   = cb - ca
     rab = sqrt(dot_product(d,d))

     if (rab < same_tol) then
        ! One center: the radial integral is elementary and the angular one is
        ! done in closed form. For the spherical harmonics it is a Kronecker
        ! delta; for the Cartesian monomials it is neither diagonal nor zero
        ! for la /= lb, <s|x^2> being the simplest counterexample.
        if (ang == 1) then
           call one_center_cartesian (la,lb,block)
           block = block * sto_norm(na,za)*sto_norm(nb,zb) &
                         * gamma_integral(na+nb,za+zb)
        else
           if (la /= lb) return
           do m = 1, 2*la+1
              block(m,m) = sto_norm(na,za)*sto_norm(nb,zb) &
                         * gamma_integral(na+nb,za+zb)
           end do
        end if
     else
        call two_center_block (na,la,za, nb,lb,zb, d,rab, block, ang)
     end if
  end subroutine sto_primitive_overlap

! =============================================================================
!  The two-center integral
! =============================================================================

  subroutine two_center_block (na,la,za, nb,lb,zb, d,rab, block, ang)
     integer,  intent(in)  :: na, la, nb, lb, ang
     real(wp), intent(in)  :: za, zb, d(3), rab
     real(wp), intent(out) :: block(:,:)

     real(wp) :: ehat(3), uhat(3), vhat(3)
     real(wp) :: a, b, halfR, rhalf_xi, eta, rho, rho2, zloc, expfac, vol
     real(wp) :: ra, rb, wgt, radial, cphi, sphi
     real(wp) :: pa(3), pb(3)
     real(wp), allocatable :: tx(:), tw(:), ex(:), ew(:)
     real(wp), allocatable :: ya(:), yb(:), cp(:), sp(:)
     integer  :: nxi, neta, nphi, it, ie, ip, ia, ib
     integer  :: pwa, pwb, nca, ncb

     block = 0.0_wp

     a = 0.5_wp*(za+zb)*rab
     b = 0.5_wp*(za-zb)*rab

     ! The integrand is bounded by exp(-(a-|b|)) = exp(-min(za,zb)*R) times a
     ! polynomial, so a well separated pair with tight exponents contributes
     ! nothing representable.
     if (min(za,zb)*rab > screen_exponent) return

     call local_frame (d,rab,ehat,uhat,vhat)

     nxi  = (na+nb)/2 + 2
     nphi = la + lb + 2
     neta = max((na+nb)/2 + 2, int(abs(b)) + extra_eta_nodes)
     neta = min(neta,max_eta_nodes)

     call cached_laguerre (nxi,tx,tw)
     call cached_legendre (neta,ex,ew)

     halfR = 0.5_wp*rab
     pwa   = na - 1 - la
     pwb   = nb - 1 - lb

     ! The angular factors are homogeneous of degree la and lb whichever
     ! convention is in force, so nxi, neta and nphi above are unaffected.
     nca = sto_ncomponents(la,ang)
     ncb = sto_ncomponents(lb,ang)
     allocate (ya(nca), yb(ncb), cp(nphi), sp(nphi))
     do ip = 1, nphi
        cp(ip) = cos(2.0_wp*pi*real(ip-1,wp)/real(nphi,wp))
        sp(ip) = sin(2.0_wp*pi*real(ip-1,wp)/real(nphi,wp))
     end do

     do it = 1, nxi
        ! xi = 1 + t/a, written so that R*xi/2 stays a physical length and no
        ! large power of xi is ever formed.
        rhalf_xi = halfR + tx(it)/(za+zb)
        do ie = 1, neta
           eta    = ex(ie)
           ra     = rhalf_xi + halfR*eta
           rb     = rhalf_xi - halfR*eta
           rho2   = (rhalf_xi - halfR)*(rhalf_xi + halfR)*(1.0_wp-eta)*(1.0_wp+eta)
           if (rho2 < 0.0_wp) rho2 = 0.0_wp
           rho    = sqrt(rho2)
           zloc   = halfR + eta*rhalf_xi
           expfac = exp(-(a + b*eta))
           ! (R^3/8)(xi^2-eta^2)/a, in physical units
           vol    = (rab/(2.0_wp*a))*(rhalf_xi**2 - (halfR*eta)**2)
           radial = tw(it)*ew(ie)*expfac*vol*(ra**pwa)*(rb**pwb)
           if (radial == 0.0_wp) cycle
           radial = radial*(2.0_wp*pi/real(nphi,wp))
           do ip = 1, nphi
              cphi = cp(ip)
              sphi = sp(ip)
              ! global position relative to each center
              pa = rho*(cphi*uhat + sphi*vhat) + zloc*ehat
              pb = pa - rab*ehat
              if (ang == 1) then
                 call cartesian_monomials (la,pa,ya)
                 call cartesian_monomials (lb,pb,yb)
              else
                 call real_solid_harmonics (la,pa,ya)
                 call real_solid_harmonics (lb,pb,yb)
              end if
              do ib = 1, ncb
                 wgt = radial*yb(ib)
                 do ia = 1, nca
                    block(ia,ib) = block(ia,ib) + wgt*ya(ia)
                 end do
              end do
           end do
        end do
     end do

     block = block * sto_norm(na,za) * sto_norm(nb,zb)

     deallocate (ya,yb,cp,sp,tx,tw,ex,ew)
  end subroutine two_center_block


  subroutine cached_laguerre (n,x,w)
     integer,               intent(in)  :: n
     real(wp), allocatable, intent(out) :: x(:), w(:)
     if (n /= cached_nlag) then
        if (allocated(cached_xlag)) deallocate (cached_xlag,cached_wlag)
        call gauss_laguerre (n,cached_xlag,cached_wlag)
        cached_nlag = n
     end if
     allocate (x(n),w(n))
     x = cached_xlag
     w = cached_wlag
  end subroutine cached_laguerre


  subroutine cached_legendre (n,x,w)
     integer,               intent(in)  :: n
     real(wp), allocatable, intent(out) :: x(:), w(:)
     if (n /= cached_nleg) then
        if (allocated(cached_xleg)) deallocate (cached_xleg,cached_wleg)
        call gauss_legendre (n,cached_xleg,cached_wleg)
        cached_nleg = n
     end if
     allocate (x(n),w(n))
     x = cached_xleg
     w = cached_wleg
  end subroutine cached_legendre


  subroutine local_frame (d,rab,ehat,uhat,vhat)
  ! Orthonormal frame with ehat along the internuclear vector. The choice of
  ! uhat and vhat is arbitrary: the azimuthal integration is exact, so the
  ! result does not depend on it.
     real(wp), intent(in)  :: d(3), rab
     real(wp), intent(out) :: ehat(3), uhat(3), vhat(3)
     real(wp) :: t(3), tn

     ehat = d/rab
     if (abs(ehat(1)) < 0.9_wp) then
        t = [1.0_wp,0.0_wp,0.0_wp]
     else
        t = [0.0_wp,1.0_wp,0.0_wp]
     end if
     uhat = t - dot_product(t,ehat)*ehat
     tn   = sqrt(dot_product(uhat,uhat))
     uhat = uhat/tn
     vhat = [ ehat(2)*uhat(3)-ehat(3)*uhat(2), &
              ehat(3)*uhat(1)-ehat(1)*uhat(3), &
              ehat(1)*uhat(2)-ehat(2)*uhat(1) ]
  end subroutine local_frame


  real(wp) function gamma_integral (k,p) result (g)
  ! int_0^infinity r^k exp(-p r) dr = k! / p^{k+1}, in logarithms.
     integer,  intent(in) :: k
     real(wp), intent(in) :: p
     integer  :: i
     real(wp) :: logfac
     logfac = 0.0_wp
     do i = 2, k
        logfac = logfac + log(real(i,wp))
     end do
     g = exp(logfac - real(k+1,wp)*log(p))
  end function gamma_integral

! =============================================================================
!  Cartesian angular factors
! =============================================================================
!
! ADF stores its Slater functions as Cartesian powers,
!
!     chi = N x^kx y^ky z^kz r^kr exp(-alpha r),
!
! with kx+ky+kz = l and kr = n-1-l, so the angular factor is a monomial rather
! than a spherical harmonic and a d shell has six components, not five. The
! routines below supply that factor in the same normalisation as
! real_solid_harmonics, namely unit mean square over the unit sphere, which is
! what makes the rest of the module indifferent to the choice: the radial
! factor and sto_norm are identical, and the quadrature's exactness argument
! rests only on the angular factor being a polynomial of degree l.

  integer function sto_ncomponents (l,angular) result (nc)
  ! Number of components of a shell of angular momentum l: the Cartesian
  ! monomials of degree l, or the 2l+1 real spherical harmonics.
     integer, intent(in) :: l, angular
     if (angular == 1) then
        nc = (l+1)*(l+2)/2
     else
        nc = 2*l + 1
     end if
  end function sto_ncomponents


  subroutine cart_powers (l,i,kx,ky,kz)
  ! Cartesian powers of component i of a shell of degree l, in ROSE's order:
  ! kx descending outermost, ky descending within it. For l = 2 this gives
  ! xx, xy, xz, yy, yz, zz, which is also the order ADF writes.
     integer, intent(in)  :: l, i
     integer, intent(out) :: kx, ky, kz
     integer :: jx, jy, j
     j = 0
     do jx = l, 0, -1
        do jy = l-jx, 0, -1
           j = j + 1
           if (j == i) then
              kx = jx ; ky = jy ; kz = l - jx - jy
              return
           end if
        end do
     end do
     error stop 'cart_powers: component index outside the shell'
  end subroutine cart_powers


  real(wp) function dfact (n) result (f)
  ! Double factorial n!!, with the usual convention (-1)!! = 0!! = 1.
     integer, intent(in) :: n
     integer :: i
     f = 1.0_wp
     do i = n, 2, -2
        f = f*real(i,wp)
     end do
  end function dfact


  real(wp) function cart_angular_norm (kx,ky,kz) result (c)
  ! Normalisation of the angular factor (x/r)^kx (y/r)^ky (z/r)^kz, fixed by
  ! requiring unit mean square over the unit sphere:
  !
  !     int (x/r)^2kx (y/r)^2ky (z/r)^2kz dOmega
  !            = 4 pi (2kx-1)!!(2ky-1)!!(2kz-1)!! / (2l+1)!!
  !
  ! so that  c = sqrt( (2l+1)!! / (4 pi (2kx-1)!!(2ky-1)!!(2kz-1)!!) ).
  ! The full normalisation of a Cartesian Slater primitive is then
  ! sto_norm(n,zeta)*cart_angular_norm(kx,ky,kz), the radial part being the
  ! same as in the spherical case. For l = 2 this reproduces ADF's bnorm,
  ! including the factor sqrt(3) between xx and xy; see t_cartesian.
     integer, intent(in) :: kx, ky, kz
     integer :: l
     l = kx + ky + kz
     c = sqrt( dfact(2*l+1) &
             / (4.0_wp*pi*dfact(2*kx-1)*dfact(2*ky-1)*dfact(2*kz-1)) )
  end function cart_angular_norm


  subroutine cartesian_monomials (l,r,y)
  ! Normalised Cartesian angular factors of degree l evaluated at r, i.e.
  ! y(i) = c_i x^kx y^ky z^kz with the r^l folded in, ordered by cart_powers.
  ! The Cartesian counterpart of real_solid_harmonics, and interchangeable
  ! with it everywhere in this module.
     integer,  intent(in)  :: l
     real(wp), intent(in)  :: r(3)
     real(wp), intent(out) :: y(:)
     integer :: kx, ky, kz, i
     i = 0
     do kx = l, 0, -1
        do ky = l-kx, 0, -1
           kz = l - kx - ky
           i  = i + 1
           y(i) = cart_angular_norm(kx,ky,kz) &
                * r(1)**kx * r(2)**ky * r(3)**kz
        end do
     end do
  end subroutine cartesian_monomials


  subroutine one_center_cartesian (la,lb,block)
  ! Angular part of the one-center overlap of two normalised Cartesian angular
  ! factors. With K = ka + kb componentwise,
  !
  !     int Y_a Y_b dOmega = c_a c_b 4 pi (Kx-1)!!(Ky-1)!!(Kz-1)!! / (K+1)!!
  !
  ! when Kx, Ky and Kz are all even, and zero otherwise. This is not diagonal,
  ! and it does not vanish for la /= lb.
     integer,  intent(in)  :: la, lb
     real(wp), intent(out) :: block(:,:)
     integer :: ia, ib, ax, ay, az, bx, by, bz, kx, ky, kz
     block = 0.0_wp
     do ia = 1, (la+1)*(la+2)/2
        call cart_powers (la,ia,ax,ay,az)
        do ib = 1, (lb+1)*(lb+2)/2
           call cart_powers (lb,ib,bx,by,bz)
           kx = ax + bx ; ky = ay + by ; kz = az + bz
           if (mod(kx,2) /= 0 .or. mod(ky,2) /= 0 .or. mod(kz,2) /= 0) cycle
           block(ia,ib) = cart_angular_norm(ax,ay,az)*cart_angular_norm(bx,by,bz) &
                        * 4.0_wp*pi*dfact(kx-1)*dfact(ky-1)*dfact(kz-1) &
                        / dfact(kx+ky+kz+1)
        end do
     end do
  end subroutine one_center_cartesian


! =============================================================================
!  Real solid harmonics
! =============================================================================

  subroutine real_solid_harmonics (l,r,y)
  ! Real solid harmonics S_lm(x,y,z) = r^l Z_lm(theta,phi) for all m of one l,
  ! ordered m = -l..l, with Z_lm the real spherical harmonics normalized to
  ! unity on the unit sphere.
  !
  ! Built from the standard associated Legendre recurrences in the form that
  ! avoids the explicit sin(theta)^m factor, so that the result is polynomial
  ! in the Cartesian components and nothing is divided by r.
     integer,  intent(in)  :: l
     real(wp), intent(in)  :: r(3)
     real(wp), intent(out) :: y(:)

     real(wp) :: x, yy, z, r2, c, s, ctmp, stmp
     real(wp) :: p(0:l)
     integer  :: m
     real(wp) :: fac

     if (.not. harm_ready) call build_harm_norm
     if (l > lmax_table) error stop 'real_solid_harmonics: l above the tabulated maximum'

     x = r(1) ; yy = r(2) ; z = r(3)
     r2 = x*x + yy*yy + z*z

     ! p(m) = r^l * P_l^m(z/r) / (x^2+y^2)^{m/2}, i.e. the associated Legendre
     ! function with the (1-u^2)^{m/2} factor stripped and the r^l folded in.
     ! Recurrences, with u = z/r and w^2 = x^2+y^2:
     !   q_m^m = (2m-1)!! * (-1)^m        (Condon-Shortley included below)
     !   q_{m+1}^m = (2m+1) z q_m^m
     !   (l-m) q_l^m = (2l-1) z q_{l-1}^m - (l+m-1) r^2 q_{l-2}^m
     call legendre_stripped (l,z,r2,p)

     ! c = Re (x + i y)^m, s = Im (x + i y)^m, built by recursion
     y(l+1) = harm_norm(l,0)*p(0)

     c = 1.0_wp ; s = 0.0_wp
     do m = 1, l
        ctmp = c*x - s*yy
        stmp = c*yy + s*x
        c = ctmp ; s = stmp
        fac = harm_norm(l,m)*p(m)
        y(l+1+m) = fac*c
        y(l+1-m) = fac*s
     end do

  end subroutine real_solid_harmonics


  subroutine build_harm_norm
  ! harm_norm(l,m) = sqrt((2l+1)/4pi) sqrt((l-m)!/(l+m)!) and, for m > 0, an
  ! extra sqrt(2) from the real combination.
     integer  :: l, m, k
     real(wp) :: a, f
     harm_norm = 0.0_wp
     do l = 0, lmax_table
        a = sqrt(real(2*l+1,wp)/(4.0_wp*pi))
        harm_norm(l,0) = a
        do m = 1, l
           f = 0.0_wp
           do k = l-m+1, l+m
              f = f + log(real(k,wp))
           end do
           harm_norm(l,m) = a*sqrt(2.0_wp)*exp(-0.5_wp*f)
        end do
     end do
     harm_ready = .true.
  end subroutine build_harm_norm


  subroutine legendre_stripped (l,z,r2,p)
  ! p(m) = r^l P_l^m(z/r) / (x^2+y^2)^{m/2} for m = 0..l, with the
  ! Condon-Shortley phase removed, i.e. p(m) >= 0 for z > 0.
     integer,  intent(in)  :: l
     real(wp), intent(in)  :: z, r2
     real(wp), intent(out) :: p(0:l)

     real(wp) :: pmm, pmmp1, pll
     integer  :: m, ll

     p = 0.0_wp
     do m = 0, l
        ! q_m^m with the Condon-Shortley phase removed: (2m-1)!!
        pmm = 1.0_wp
        do ll = 1, 2*m-1, 2
           pmm = pmm*real(ll,wp)
        end do
        if (m == l) then
           p(m) = pmm
           cycle
        end if
        pmmp1 = real(2*m+1,wp)*z*pmm
        if (m+1 == l) then
           p(m) = pmmp1
           cycle
        end if
        do ll = m+2, l
           pll   = (real(2*ll-1,wp)*z*pmmp1 - real(ll+m-1,wp)*r2*pmm)/real(ll-m,wp)
           pmm   = pmmp1
           pmmp1 = pll
        end do
        p(m) = pmmp1
     end do
  end subroutine legendre_stripped

! =============================================================================
!  Quadrature rules
! =============================================================================

  subroutine gauss_legendre (n,x,w)
  ! Nodes and weights of the n-point Gauss-Legendre rule on [-1,1], by Newton
  ! iteration on the Legendre polynomial evaluated through its three-term
  ! recurrence.
     integer,               intent(in)  :: n
     real(wp), allocatable, intent(out) :: x(:), w(:)
     integer  :: i, j, it
     real(wp) :: z, z1, p0, p1, p2, pp

     allocate (x(n),w(n))
     do i = 1, (n+1)/2
        z = cos(pi*(real(i,wp)-0.25_wp)/(real(n,wp)+0.5_wp))
        do it = 1, 100
           p1 = 1.0_wp ; p2 = 0.0_wp
           do j = 1, n
              p0 = p2 ; p2 = p1
              p1 = (real(2*j-1,wp)*z*p2 - real(j-1,wp)*p0)/real(j,wp)
           end do
           pp = real(n,wp)*(z*p1 - p2)/(z*z - 1.0_wp)
           z1 = z
           z  = z1 - p1/pp
           if (abs(z-z1) <= 1.0e-15_wp) exit
        end do
        x(i)     = -z
        x(n+1-i) =  z
        w(i)     = 2.0_wp/((1.0_wp - z*z)*pp*pp)
        w(n+1-i) = w(i)
     end do
  end subroutine gauss_legendre


  subroutine gauss_laguerre (n,x,w)
  ! Nodes and weights of the n-point Gauss-Laguerre rule for the weight
  ! exp(-t) on [0,infinity), by Newton iteration on L_n evaluated through its
  ! three-term recurrence. The starting values are the customary ones for the
  ! first two roots and linear extrapolation thereafter.
     integer,               intent(in)  :: n
     real(wp), allocatable, intent(out) :: x(:), w(:)
     integer  :: i, j, it
     real(wp) :: z, z1, p0, p1, p2, pp, rn

     allocate (x(n),w(n))
     rn = real(n,wp)
     z  = 0.0_wp
     do i = 1, n
        if (i == 1) then
           z = 3.0_wp/(1.0_wp + 2.4_wp*rn)
        else if (i == 2) then
           z = z + 15.0_wp/(1.0_wp + 2.6_wp*rn)
        else
           z = z + (1.0_wp + 2.55_wp*real(i-2,wp))/(1.9_wp*real(i-2,wp))*(z - x(i-2))
        end if
        do it = 1, 200
           p1 = 1.0_wp ; p2 = 0.0_wp
           do j = 1, n
              p0 = p2 ; p2 = p1
              p1 = ((real(2*j-1,wp) - z)*p2 - real(j-1,wp)*p0)/real(j,wp)
           end do
           pp = (rn*p1 - rn*p2)/z
           z1 = z
           z  = z1 - p1/pp
           if (abs(z-z1) <= 1.0e-14_wp*abs(z)) exit
        end do
        x(i) = z
        w(i) = -1.0_wp/(pp*rn*p2)
     end do
  end subroutine gauss_laguerre


! =============================================================================
!  Self-test
! =============================================================================

  subroutine sto_overlap_selftest (nfail,verbose)
  ! Checks the module against results that are known independently of it:
  ! orthonormality of the harmonics, the analytic one-center integral, the
  ! classical closed form for 1s-1s, the vanishing of off-diagonal m on the
  ! axis, invariance of the block norm under rotation of the pair, bra-ket
  ! symmetry, continuity into the one-center limit, and a direct quadrature in
  ! spherical coordinates that shares none of the machinery under test.
     integer, intent(out)          :: nfail
     logical, intent(in), optional :: verbose

     logical :: loud

     loud = .true.
     if (present(verbose)) loud = verbose
     nfail = 0

     call t_harmonics        (nfail,loud)
     call t_one_center       (nfail,loud)
     call t_1s1s             (nfail,loud)
     call t_axis_structure   (nfail,loud)
     call t_rotation         (nfail,loud)
     call t_symmetry         (nfail,loud)
     call t_small_r          (nfail,loud)
     call t_brute_force      (nfail,loud)
     call t_shells           (nfail,loud)
     call t_cartesian        (nfail,loud)

     if (loud) then
        write(*,*)
        if (nfail == 0) then
           write(*,'(A)') ' slater_overlap: all checks passed'
        else
           write(*,'(A,I0,A)') ' slater_overlap: ',nfail,' check(s) FAILED'
        end if
     end if
  end subroutine sto_overlap_selftest


  subroutine t_shells (nfail,loud)
  ! The contracted-shell and matrix drivers: a one-primitive shell must
  ! reproduce the primitive, a normalised contraction must give unit diagonal,
  ! and a basis overlap matrix must be symmetric and positive definite.
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     type(slater_func_info_t) :: sh(3)
     real(wp) :: prim(3,3), blk(3,3), err
     real(wp), allocatable :: smat(:,:), ev(:)
     integer  :: i, j, n

     err = 0.0_wp

     ! single primitive shell versus the primitive routine
     sh(1)%orb_momentum = 2 ; sh(1)%n_primitives = 1 ; sh(1)%coord = [0.0_wp,0.0_wp,0.0_wp]
     allocate (sh(1)%n(1),sh(1)%exponent(1),sh(1)%coefficient(1))
     sh(1)%n = 2 ; sh(1)%exponent = 1.5_wp ; sh(1)%coefficient = 1.0_wp
     sh(2)%orb_momentum = 2 ; sh(2)%n_primitives = 1 ; sh(2)%coord = [0.0_wp,0.4_wp,1.3_wp]
     allocate (sh(2)%n(1),sh(2)%exponent(1),sh(2)%coefficient(1))
     sh(2)%n = 3 ; sh(2)%exponent = 1.1_wp ; sh(2)%coefficient = 1.0_wp
     call sto_shell_overlap (sh(1),sh(2),blk)
     call sto_primitive_overlap (2,1,1.5_wp,sh(1)%coord, 3,1,1.1_wp,sh(2)%coord, prim)
     err = max(err,maxval(abs(blk-prim)))

     ! a two-primitive shell, normalised through its own self-overlap
     sh(3)%orb_momentum = 1 ; sh(3)%n_primitives = 2 ; sh(3)%coord = [0.7_wp,0.0_wp,-0.2_wp]
     allocate (sh(3)%n(2),sh(3)%exponent(2),sh(3)%coefficient(2))
     sh(3)%n = [1,2] ; sh(3)%exponent = [1.8_wp,1.0_wp] ; sh(3)%coefficient = [0.6_wp,0.5_wp]
     call sto_shell_overlap (sh(3),sh(3),blk(1:1,1:1))
     sh(3)%coefficient = sh(3)%coefficient/sqrt(blk(1,1))
     call sto_shell_overlap (sh(3),sh(3),blk(1:1,1:1))
     err = max(err,abs(blk(1,1)-1.0_wp))

     ! the matrix driver
     n = 0
     do i = 1, 3
        n = n + 2*(sh(i)%orb_momentum-1) + 1
     end do
     allocate (smat(n,n),ev(n))
     call sto_overlap_matrix (sh,sh,smat)
     do i = 1, n
        do j = 1, n
           err = max(err,abs(smat(i,j)-smat(j,i)))
        end do
     end do
     call small_eigenvalues (smat,ev)
     if (minval(ev) <= 0.0_wp) err = max(err,1.0_wp)
     deallocate (smat,ev)
     do i = 1, 3
        deallocate (sh(i)%n,sh(i)%exponent,sh(i)%coefficient)
     end do

     call t_report ('contracted shells and matrix',err,1.0e-13_wp,nfail,loud)
  end subroutine t_shells


  subroutine small_eigenvalues (a,ev)
  ! Jacobi eigenvalues of a small symmetric matrix, used only by the self-test
  ! so that the module stays free of LAPACK.
     real(wp), intent(in)  :: a(:,:)
     real(wp), intent(out) :: ev(:)
     real(wp), allocatable :: m(:,:)
     real(wp) :: off, theta, t, c, s, aip, aiq
     integer  :: n, p, q, i, sweep
     n = size(a,1)
     allocate (m(n,n))
     m = a
     do sweep = 1, 100
        off = 0.0_wp
        do p = 1, n-1
           do q = p+1, n
              off = off + m(p,q)**2
           end do
        end do
        if (off < 1.0e-30_wp) exit
        do p = 1, n-1
           do q = p+1, n
              if (abs(m(p,q)) < 1.0e-18_wp) cycle
              theta = (m(q,q)-m(p,p))/(2.0_wp*m(p,q))
              t = sign(1.0_wp,theta)/(abs(theta)+sqrt(theta*theta+1.0_wp))
              c = 1.0_wp/sqrt(t*t+1.0_wp)
              s = t*c
              do i = 1, n
                 aip = m(i,p) ; aiq = m(i,q)
                 m(i,p) = c*aip - s*aiq
                 m(i,q) = s*aip + c*aiq
              end do
              do i = 1, n
                 aip = m(p,i) ; aiq = m(q,i)
                 m(p,i) = c*aip - s*aiq
                 m(q,i) = s*aip + c*aiq
              end do
           end do
        end do
     end do
     do i = 1, n
        ev(i) = m(i,i)
     end do
     deallocate (m)
  end subroutine small_eigenvalues


  subroutine t_report (name,err,tol,nfail,loud)
     character(*), intent(in)    :: name
     real(wp),     intent(in)    :: err, tol
     integer,      intent(inout) :: nfail
     logical,      intent(in)    :: loud
     if (err > tol) nfail = nfail + 1
     if (loud) then
        if (err <= tol) then
           write(*,'(A,A34,A,ES10.2,A,ES8.1,A)') ' ok   ',name,' err=',err,'  (tol ',tol,')'
        else
           write(*,'(A,A34,A,ES10.2,A,ES8.1,A)') ' FAIL ',name,' err=',err,'  (tol ',tol,')'
        end if
     end if
  end subroutine t_report


  subroutine t_harmonics (nfail,loud)
  ! Gauss-Legendre in cos(theta) with the trapezoidal rule in phi is exact for
  ! products of two spherical harmonics, so this isolates the harmonics.
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     integer  :: l, i, k, ma, mb, nu, np
     real(wp) :: ph, r(3), st, err, v
     real(wp), allocatable :: y(:), acc(:,:), u(:), wu(:)
     err = 0.0_wp
     do l = 0, 5
        nu = l + 2
        np = 2*l + 2
        call gauss_legendre (nu,u,wu)
        allocate (y(2*l+1), acc(2*l+1,2*l+1))
        acc = 0.0_wp
        do i = 1, nu
           st = sqrt(max(0.0_wp,1.0_wp-u(i)**2))
           do k = 1, np
              ph = 2.0_wp*pi*real(k-1,wp)/real(np,wp)
              r  = [st*cos(ph), st*sin(ph), u(i)]
              call real_solid_harmonics (l,r,y)
              do mb = 1, 2*l+1
                 do ma = 1, 2*l+1
                    acc(ma,mb) = acc(ma,mb) + wu(i)*(2.0_wp*pi/real(np,wp))*y(ma)*y(mb)
                 end do
              end do
           end do
        end do
        do mb = 1, 2*l+1
           do ma = 1, 2*l+1
              v = acc(ma,mb)
              if (ma == mb) v = v - 1.0_wp
              err = max(err,abs(v))
           end do
        end do
        deallocate (y,acc,u,wu)
     end do
     call t_report ('real harmonics orthonormal',err,1.0e-13_wp,nfail,loud)
  end subroutine t_harmonics


  subroutine t_one_center (nfail,loud)
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     real(wp) :: blk(9,9), err, c(3)
     integer  :: n, l, m
     c   = [0.1_wp,-0.3_wp,0.7_wp]
     err = 0.0_wp
     do n = 1, 6
        do l = 0, min(n-1,4)
           call sto_primitive_overlap (n,l,1.7_wp,c, n,l,1.7_wp,c, blk(1:2*l+1,1:2*l+1))
           do m = 1, 2*l+1
              err = max(err,abs(blk(m,m)-1.0_wp))
           end do
        end do
     end do
     call t_report ('one-center normalization',err,1.0e-13_wp,nfail,loud)
  end subroutine t_one_center


  subroutine t_1s1s (nfail,loud)
  ! S(1s,1s) = exp(-p) (1 + p + p^2/3) with p = zeta R, for equal exponents.
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     real(wp) :: blk(1,1), z, rr, p, ref, err
     integer  :: i
     err = 0.0_wp
     z   = 1.3_wp
     do i = 1, 16
        rr  = 0.25_wp*real(i,wp)
        p   = z*rr
        ref = exp(-p)*(1.0_wp + p + p*p/3.0_wp)
        call sto_primitive_overlap (1,0,z,[0.0_wp,0.0_wp,0.0_wp], &
                                    1,0,z,[0.0_wp,0.0_wp,rr], blk)
        err = max(err,abs(blk(1,1)-ref))
     end do
     call t_report ('1s-1s closed form',err,1.0e-13_wp,nfail,loud)
  end subroutine t_1s1s


  subroutine t_axis_structure (nfail,loud)
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     real(wp) :: blk(9,9), err
     integer  :: ia, ib, la, lb
     err = 0.0_wp
     do la = 0, 3
        do lb = 0, 3
           call sto_primitive_overlap (la+1,la,1.4_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                       lb+2,lb,0.9_wp,[0.0_wp,0.0_wp,1.8_wp], &
                                       blk(1:2*la+1,1:2*lb+1))
           do ib = 1, 2*lb+1
              do ia = 1, 2*la+1
                 if (ia-la-1 /= ib-lb-1) err = max(err,abs(blk(ia,ib)))
              end do
           end do
        end do
     end do
     call t_report ('m-diagonal on the axis',err,1.0e-14_wp,nfail,loud)
  end subroutine t_axis_structure


  subroutine t_rotation (nfail,loud)
  ! The block transforms as D_a S D_b^T with D orthogonal, so its Frobenius
  ! norm does not depend on the orientation of the pair.
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     real(wp) :: blk(9,9), c2(3), fn, fn0, err, al, be
     integer  :: k, la, lb
     err = 0.0_wp
     do la = 0, 3
        do lb = 0, 3
           fn0 = -1.0_wp
           do k = 0, 6
              al = 0.37_wp*real(k,wp)
              be = 0.21_wp*real(k,wp) + 0.4_wp
              c2 = 2.1_wp*[sin(be)*cos(al), sin(be)*sin(al), cos(be)]
              call sto_primitive_overlap (la+2,la,1.1_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                          lb+1,lb,1.6_wp,c2, blk(1:2*la+1,1:2*lb+1))
              fn = sqrt(sum(blk(1:2*la+1,1:2*lb+1)**2))
              if (fn0 < 0.0_wp) then
                 fn0 = fn
              else
                 err = max(err,abs(fn-fn0))
              end if
           end do
        end do
     end do
     call t_report ('rotational invariance',err,1.0e-13_wp,nfail,loud)
  end subroutine t_rotation


  subroutine t_symmetry (nfail,loud)
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     real(wp) :: b1(9,9), b2(9,9), err, c1(3), c2(3)
     integer  :: la, lb, ia, ib
     c1  = [0.2_wp,0.1_wp,-0.4_wp]
     c2  = [-1.0_wp,0.8_wp,1.3_wp]
     err = 0.0_wp
     do la = 0, 3
        do lb = 0, 3
           call sto_primitive_overlap (la+2,la,1.2_wp,c1, lb+1,lb,2.0_wp,c2, b1(1:2*la+1,1:2*lb+1))
           call sto_primitive_overlap (lb+1,lb,2.0_wp,c2, la+2,la,1.2_wp,c1, b2(1:2*lb+1,1:2*la+1))
           do ib = 1, 2*lb+1
              do ia = 1, 2*la+1
                 err = max(err,abs(b1(ia,ib)-b2(ib,ia)))
              end do
           end do
        end do
     end do
     call t_report ('bra-ket symmetry',err,1.0e-14_wp,nfail,loud)
  end subroutine t_symmetry


  subroutine t_small_r (nfail,loud)
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     real(wp) :: blk(3,3), err, rr
     integer  :: i
     err = 0.0_wp
     do i = 1, 6
        rr = 1.0e-3_wp/real(2**i,wp)
        call sto_primitive_overlap (2,1,1.5_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                    2,1,1.5_wp,[0.0_wp,0.0_wp,rr], blk)
        err = max(err,abs(blk(1,1)-1.0_wp),abs(blk(2,2)-1.0_wp),abs(blk(3,3)-1.0_wp))
     end do
     call t_report ('small-R continuity',err,1.0e-5_wp,nfail,loud)
  end subroutine t_small_r


  subroutine t_brute_force (nfail,loud)
  ! Direct quadrature in spherical coordinates about the first center. It uses
  ! none of the spheroidal machinery, so it is an independent check of the
  ! whole construction. It converges slowly because of the cusp at the second
  ! center, hence the loose tolerance.
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     real(wp) :: blk(7,7), ref, err, c2(3)
     integer  :: la, lb, na, nb, ia, ib
     err = 0.0_wp
     c2  = [0.0_wp,0.0_wp,1.7_wp]
     do la = 0, 2
        do lb = 0, 2
           na = la + 1
           nb = lb + 2
           call sto_primitive_overlap (na,la,1.3_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                       nb,lb,0.95_wp,c2, blk(1:2*la+1,1:2*lb+1))
           do ib = 1, 2*lb+1
              do ia = 1, 2*la+1
                 call t_brute (na,la,ia,1.3_wp, nb,lb,ib,0.95_wp, c2, 2, ref)
                 err = max(err,abs(blk(ia,ib)-ref))
              end do
           end do
        end do
     end do
     ! The difference is second order in the brute-force step and has been
     ! verified to fall by four on each halving, so it is that rule's error.
     call t_report ('vs independent quadrature',err,1.0e-4_wp,nfail,loud)
  end subroutine t_brute_force


  subroutine t_brute (na,la,ia,za, nb,lb,ib,zb, c2, ang, s)
  ! ia and ib are component indices in the convention ang.
     integer,  intent(in)  :: na, la, ia, nb, lb, ib, ang
     real(wp), intent(in)  :: za, zb, c2(3)
     real(wp), intent(out) :: s
     integer, parameter :: nr = 240, nu = 240, np = 24
     integer  :: i, j, k
     real(wp) :: r, ph, pa(3), pb(3), rb, va, vb, st, rmax, dr, du, w
     real(wp), allocatable :: ya(:), yb(:), uu(:)
     allocate (ya(sto_ncomponents(la,ang)), yb(sto_ncomponents(lb,ang)))
     rmax = 14.0_wp
     dr   = rmax/real(nr,wp)
     du   = 2.0_wp/real(nu,wp)
     allocate (uu(nu))
     do j = 1, nu
        uu(j) = -1.0_wp + du*(real(j,wp)-0.5_wp)
     end do
     s = 0.0_wp
     do i = 1, nr
        r = dr*(real(i,wp)-0.5_wp)
        do j = 1, nu
           st = sqrt(max(0.0_wp,1.0_wp-uu(j)**2))
           do k = 1, np
              ph = 2.0_wp*pi*(real(k,wp)-0.5_wp)/real(np,wp)
              pa = r*[st*cos(ph), st*sin(ph), uu(j)]
              pb = pa - c2
              rb = sqrt(dot_product(pb,pb))
              if (rb < 1.0e-12_wp) cycle
              if (ang == 1) then
                 call cartesian_monomials (la,pa,ya)
                 call cartesian_monomials (lb,pb,yb)
              else
                 call real_solid_harmonics (la,pa,ya)
                 call real_solid_harmonics (lb,pb,yb)
              end if
              va = r**(na-1-la)*exp(-za*r)*ya(ia)
              vb = rb**(nb-1-lb)*exp(-zb*rb)*yb(ib)
              w  = r*r*dr*du*(2.0_wp*pi/real(np,wp))
              s  = s + w*va*vb
           end do
        end do
     end do
     s = s*sto_norm(na,za)*sto_norm(nb,zb)
     deallocate (ya,yb,uu)
  end subroutine t_brute

  subroutine t_cartesian (nfail,loud)
  ! The Cartesian convention. Five independent checks: the normalisation
  ! against the constants ADF itself prints, the exactness of that
  ! normalisation through the one-center self-overlap, agreement with the
  ! already validated spherical path wherever the two bases coincide (l <= 1),
  ! bra-ket symmetry, and a direct quadrature for a d shell, where the six
  ! Cartesian components have no spherical counterpart.
     integer, intent(inout) :: nfail
     logical, intent(in)    :: loud
     real(wp) :: err, blk(10,10), bsph(10,10), bt(10,10), ref, rr
     real(wp) :: nxx, nxy, c2(3)
     integer  :: n, l, i, j, nc, la, lb, na, nb, ia, ib
     integer, parameter :: pperm(3) = [3,1,2]   ! Cartesian x,y,z -> spherical m = +1,-1,0

     ! 1. ADF prints bnorm = 3.0091 for xx and 5.2120 for xy in a 3d shell with
     !    zeta = 2, a ratio of exactly sqrt(3). Both must come out of sto_norm
     !    times cart_angular_norm with nothing else supplied. The tolerance is
     !    set by the precision to which ADF prints them.
     nxx = sto_norm(3,2.0_wp)*cart_angular_norm(2,0,0)
     nxy = sto_norm(3,2.0_wp)*cart_angular_norm(1,1,0)
     err = max(abs(nxx-3.0090_wp),abs(nxy-5.2118_wp))
     call t_report ('Cartesian norm vs ADF bnorm',err,1.0e-4_wp,nfail,loud)
     err = abs(nxy/nxx - sqrt(3.0_wp))
     call t_report ('Cartesian norm xx:xy ratio',err,1.0e-14_wp,nfail,loud)

     ! 2. One-center self-overlap: unit diagonal for every component, which
     !    tests cart_angular_norm and one_center_cartesian together.
     err = 0.0_wp
     do n = 1, 6
        do l = 0, min(n-1,3)
           nc = sto_ncomponents(l,1)
           call sto_primitive_overlap (n,l,1.7_wp,[0.1_wp,-0.3_wp,0.7_wp], &
                                       n,l,1.7_wp,[0.1_wp,-0.3_wp,0.7_wp], &
                                       blk(1:nc,1:nc), angular=1)
           do i = 1, nc
              err = max(err,abs(blk(i,i)-1.0_wp))
           end do
        end do
     end do
     call t_report ('Cartesian one-center norm',err,1.0e-13_wp,nfail,loud)

     ! 3. For l <= 1 the two bases are the same functions in a different order,
     !    so the Cartesian path must reproduce the spherical one exactly.
     err = 0.0_wp
     c2  = [0.6_wp,-0.4_wp,1.3_wp]
     do la = 0, 1
        do lb = 0, 1
           na = la + 2 ; nb = lb + 1
           call sto_primitive_overlap (na,la,1.3_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                       nb,lb,0.95_wp,c2, blk(1:2*la+1,1:2*lb+1), angular=1)
           call sto_primitive_overlap (na,la,1.3_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                       nb,lb,0.95_wp,c2, bsph(1:2*la+1,1:2*lb+1), angular=2)
           do ib = 1, 2*lb+1
              do ia = 1, 2*la+1
                 i = ia ; j = ib
                 if (la == 1) i = pperm(ia)
                 if (lb == 1) j = pperm(ib)
                 err = max(err,abs(blk(ia,ib)-bsph(i,j)))
              end do
           end do
        end do
     end do
     call t_report ('Cartesian vs spherical, l<=1',err,1.0e-13_wp,nfail,loud)

     ! 4. Bra-ket symmetry of a d-p block.
     call sto_primitive_overlap (3,2,1.4_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                 2,1,0.9_wp,c2, blk(1:6,1:3), angular=1)
     call sto_primitive_overlap (2,1,0.9_wp,c2, &
                                 3,2,1.4_wp,[0.0_wp,0.0_wp,0.0_wp], bt(1:3,1:6), angular=1)
     err = 0.0_wp
     do i = 1, 6
        do j = 1, 3
           err = max(err,abs(blk(i,j)-bt(j,i)))
        end do
     end do
     call t_report ('Cartesian bra-ket symmetry',err,1.0e-14_wp,nfail,loud)

     ! 5. Continuity into the one-center limit, which compares the quadrature
     !    with the closed form of one_center_cartesian including its
     !    off-diagonal elements.
     err = 0.0_wp
     call sto_primitive_overlap (3,2,1.5_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                 3,2,1.5_wp,[0.0_wp,0.0_wp,0.0_wp], bsph(1:6,1:6), angular=1)
     do i = 1, 5
        rr = 1.0e-3_wp/real(2**i,wp)
        call sto_primitive_overlap (3,2,1.5_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                    3,2,1.5_wp,[0.0_wp,0.0_wp,rr], blk(1:6,1:6), angular=1)
        err = max(err,maxval(abs(blk(1:6,1:6)-bsph(1:6,1:6))))
     end do
     ! The same for la /= lb, where the spherical one-center overlap vanishes
     ! but the Cartesian one does not: <s|xx> = <s|yy> = <s|zz> /= 0.
     call sto_primitive_overlap (1,0,1.5_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                 3,2,1.5_wp,[0.0_wp,0.0_wp,0.0_wp], bsph(1:1,1:6), angular=1)
     if (maxval(abs(bsph(1,1:6))) == 0.0_wp) err = max(err,1.0_wp)
     do i = 1, 5
        rr = 1.0e-3_wp/real(2**i,wp)
        call sto_primitive_overlap (1,0,1.5_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                    3,2,1.5_wp,[0.0_wp,0.0_wp,rr], blk(1:1,1:6), angular=1)
        err = max(err,maxval(abs(blk(1,1:6)-bsph(1,1:6))))
     end do
     call t_report ('Cartesian small-R continuity',err,1.0e-5_wp,nfail,loud)

     ! 6. A d shell against the independent spherical-coordinate quadrature.
     err = 0.0_wp
     c2  = [0.0_wp,0.0_wp,1.7_wp]
     call sto_primitive_overlap (3,2,1.3_wp,[0.0_wp,0.0_wp,0.0_wp], &
                                 2,1,0.95_wp,c2, blk(1:6,1:3), angular=1)
     do ib = 1, 3
        do ia = 1, 6
           call t_brute (3,2,ia,1.3_wp, 2,1,ib,0.95_wp, c2, 1, ref)
           err = max(err,abs(blk(ia,ib)-ref))
        end do
     end do
     call t_report ('Cartesian vs indep. quadrature',err,1.0e-4_wp,nfail,loud)

  end subroutine t_cartesian


end module slater_overlap
