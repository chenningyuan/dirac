module rose_utils

   use rose_global

contains

   subroutine get_free_fileunit(unit)
    ! idea from http://fortranwiki.org, will become obsolete once fortran2008 is supported everywhere
    ! gets free file unit
      integer, intent(out) :: unit
      integer, parameter   :: LUN_MIN=100, LUN_MAX=1000
      logical :: opened
      integer :: lun
      unit=-1
      do lun=LUN_MIN,LUN_MAX
        inquire(unit=lun,opened=opened)
        if (.not. opened) then
          unit=lun
          exit
        end if
      end do
    end subroutine get_free_fileunit

    subroutine print_matrix(M)
       implicit none
       integer             :: rcq,i,j
       real(8), intent(in) :: M(:,:,:)

       rcq = size(M,3)
       if (rcq .eq. 1) then
          do i=1,size(M,1)
             write(luout,555) (M(i,j,1),j=1,size(M,2))
          enddo
       else if (rcq .eq. 2) then
          do i=1,size(M,1)
             write(luout,556) (M(i,j,1),M(i,j,2),j=1,size(M,2))
          enddo
       else if (rcq .eq. 4) then
          do i=1,size(M,1)
             write(luout,557) (M(i,j,1),M(i,j,2),M(i,j,3),M(i,j,4),j=1,size(M,2))
          enddo
       else
          print*,"Wrong algebra to print matrix."
          stop
       endif

     555 FORMAT( 100G10.2)
     556 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,')') )
     557 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,',',G10.2,',',G10.2,')') )

     end subroutine print_matrix

     subroutine SCF_density(occMOcoeff,restricted,res)
       ! get the scf density as a vector of dimension (N(N-1)/2 + N) (i.e. lower triangular matrix)
       use linear_algebra

       implicit none
       real(8), intent(inout) :: res(:,:)
       real(8), intent(in)    :: occMOcoeff(:,:,:)
       logical, intent(in)    :: restricted
       real(8), allocatable   :: density_matrix(:,:,:), occMOcoeff_bar(:,:,:)
       integer                :: i,j,k,ij,nocc,nao,rcq

       nao  = size(occMOcoeff,1)
       nocc = size(occMOcoeff,2)
       rcq  = size(occMOcoeff,3)

       allocate(density_matrix(nao,nao,rcq))
       allocate(occMOcoeff_bar(nocc,nao,rcq))
       call transpose_conjg(occMOcoeff,occMOcoeff_bar)
       call matrix_mul(occMOcoeff,occMOcoeff_bar,density_matrix)
       
       do i=1,nao
        do j=1,i
           ij=i*(i-1)/2+i
           if (restricted) then
              res(ij,:) = 2*density_matrix(i,j,:)
           else
              res(ij,:) = density_matrix(i,j,:)
           end if
           if (res(ij,1) .lt. 1.D-12) res(ij,1) = 0.D0
        enddo
       enddo
         
       deallocate(occMOcoeff_bar,density_matrix)

     end subroutine SCF_density

     function normalization_constant_psi4(l,expo) result(norm_cons)
      ! Compute the normalization constant of cartesian GTOs. (except for DIRAC interface where this is already accounted)
      implicit none
      integer,intent(in) :: l
      real(8),intent(in) :: expo
      real(8)            :: norm_cons
      real(8)            :: pi
      real(8)            :: double_factorial, factorial

      pi = 4.D0*datan(1.D0)

      ! DIRAC NORMALIZATION FACTOR:
      norm_cons = (2.D0*expo/pi)**(0.75D0)*dsqrt(2.D0**l)*dsqrt(2.D0*expo)**l

      ! GTO in cartesian coordinate we have G_{αijk} = N_{αijk} x_i y_j z_k exp[−αr^2] where alpha = exponent.
      ! such that
      ! norm_cons = (2.D0*expo/pi)**(0.75D0)*dsqrt(2.D0**l/Fijk)*sqrt(2.D0*expo)**l
      ! and
      ! Fijk = double_factorial(2*i-1)*double_factorial(2*j-1)*double_factorial(2*k-1)
      ! We ignore Fijk because it is taking into account afterwards in the get_ao_overlap subroutine. (InteRest library).
      ! In the HERMIT integral code a single normalization is chosen for each shell by ignoring Fijk
      ! meaning that Cartesian Gaussians are normalized to ⟨Gαijk|Gαijk⟩=Fijk. So we also ignore it in ROSE.

      ! From PSI4: For Cartesian functions we need to add a basis function normalization constant of
      !      _______________________________
      !     / (2lx-1)!! (2ly-1)!! (2lz-1)!!
      !    /  -----------------------------
      !  \/             (2l-1)!!
      !
      ! which is omitted in the CCA standard, adopted by Psi4.

      ! So in this normalization factor one can recognize Fijk. But the sqrt(1/double_factorial(2l-1)) is missing.
      ! Therefore, we add it to norm_cons, such that the normalization constant will be different than in DIRAC.
      ! But it will compensate with the MO coefficients.

      norm_cons = norm_cons*dsqrt(1.0/double_factorial(2*l-1))
      ! double_factorial(2l-1) is equal to 1 for s and p functions, 3 for d functions, 15 for f functions and 105 for g functions.

     end function normalization_constant_psi4

     function normalization_constant_pyscf(l,expo) result(norm_cons)
      ! Compute the normalization constant of cartesian GTOs.
      implicit none
      integer,intent(in) :: l
      real(8),intent(in) :: expo
      real(8)            :: norm_cons
      real(8)            :: pi
      real(8)            :: double_factorial, factorial
          
      pi = 4.D0*datan(1.D0)

      ! DIRAC NORMALIZATION FACTOR:
      !norm_cons = (2.D0*expo/pi)**(0.75D0)*dsqrt(2.D0**l)*dsqrt(2.D0*expo)**l

      ! In PYSCF : https://sunqm.github.io/pyscf/gto.html#pyscf.gto.mole.Mole.gto_norm
      ! However, applying it does not give correct results...
      ! norm_cons = dsqrt(2**(2.D0*l+3)*factorial(l+1)*(2*expo)**(l+1.5D0)/(factorial(2*l+2)*dsqrt(pi)))
      ! Maybe this is not for Cartesian but Spherical GTOs ?

      ! Difference between DIRAC and PYSCF GTO normalization factor is (by comparing them numerically):
      ! norm_PYSCF = (2.D0 * dsqrt(pi/double_factorial(2*l+1))) * norm_DIRAC

      ! DIRAC norm:
      norm_cons = (2.D0*expo/pi)**(0.75D0)*dsqrt(2.D0**l)*dsqrt(2.D0*expo)**l

      ! Addition for PSI4:
      norm_cons = norm_cons*dsqrt(1.0/double_factorial(2*l-1))
      
      ! Addition for PYSCF: get rid of the norm_cons from PSI4 and add the one from PYSCF.
      ! I don't know why, but one should not apply this additional factor to s and p shells.
      if (l .ge. 2) then
        norm_cons = norm_cons*2.D0*dsqrt(pi/(double_factorial(2*l+1)/double_factorial(2*l-1)))
      end if

     end function normalization_constant_pyscf

     function normalization_constant_dirac(l,expo) result(norm_cons)
      ! Compute the normalization constant of cartesian GTOs.
      implicit none
      integer,intent(in) :: l
      real(8),intent(in) :: expo
      real(8)            :: norm_cons
      real(8)            :: pi
      real(8)            :: double_factorial, factorial

      pi = 4.D0*datan(1.D0)

      ! DIRAC NORMALIZATION FACTOR:
      norm_cons = (2.D0*expo/pi)**(0.75D0)*dsqrt(2.D0**l)*dsqrt(2.D0*expo)**l

     end function normalization_constant_dirac

     integer function nfunctions(l_value,basis_angular)
      ! Computes number of Cartesian/spherical functions for a given l-value
          integer, intent(in) :: l_value
          integer, intent(in) :: basis_angular
          if (basis_angular == 1) then
!            Cartesian
             nfunctions = l_value * ( l_value + 1 ) / 2
          else if (basis_angular == 2) then
!            spherical
             nfunctions = 2 * (l_value - 1) + 1
          end if
     end function nfunctions

     subroutine make_RIBO_Fock(RIBO, CMO_energy, RIBO_Fock)

        use rose_ao
        use rose_mo
        use linear_algebra
      ! Construct the Fock Matrix in the RIBO basis for the TDDFT calculation later on

        real(8), allocatable       :: HF_Fock(:,:,:), transconjgRIBO(:,:,:)
        real(8), intent(inout),allocatable     :: RIBO_Fock(:,:,:)
        real(8), intent(in)        :: RIBO(:,:,:)
        integer                    :: i,nmo1,nact,rcq
        real(8), intent(in)        :: CMO_energy(:)

        nmo1 = size(RIBO,1)
        nact = size(RIBO,2)
        rcq = size(RIBO,3)

        allocate(HF_Fock(nmo1,nmo1,rcq))
        HF_Fock = 0.D0
        forall (i=1:nmo1) HF_Fock(i,i,1) = CMO_energy(i)

        allocate(transconjgRIBO(nact,nmo1,rcq))
        call transpose_conjg(RIBO,transconjgRIBO)

        allocate(RIBO_Fock(nact,nact,rcq))
        call matrix_mul(transconjgRIBO,HF_Fock,RIBO,RIBO_Fock)

        deallocate(transconjgRIBO,HF_Fock)

     end subroutine make_RIBO_Fock     

end module rose_utils

     ! I tried to put this function inside the module but it gives a symbol error not found, why ?
     recursive function double_factorial(n) result(m)
      implicit none
      integer, intent(in) :: n
      real(8)             :: m

      if (n .le. 0) then
         m = 1.D0
      else
         m = 1.D0*n*double_factorial(n-2)
      end if

     end function double_factorial

   function factorial (n) result (fac)

    integer, intent (in) :: n
    real(8) :: fac
    integer :: i

    if (n<0) then
       fac = 0.D0
       return
    else
       fac = 1.D0
    end if
    do i = 2, n
       fac = fac * real(i,8)
    end do

   end function factorial
