! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Module:       Kramers unrestricted SCF module  
!
! Description:  
!
! Contains:     
!
! Licensing:
!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Logs:         
!
! Comments:     
!                            
! ------------------------------------------------------------------------------------

MODULE SCF_SOLVER 

  use scf_data_types
  use scf_matrix_operations
  use scf_interest_interface

  implicit none

  !-- Public routines --!

  public kramers_unrestricted_scf


  !-- Public variables --!


  !-- Private variables --!

  private

  integer,      parameter :: lmx1 = 6 
  logical,      parameter :: rkb = .true.
  character(*), parameter :: module_name = 'kramers_unrestricted_scf_solver'

 
! type(complex_matrix) :: Fock
! type(complex_matrix) :: Pmat 
! type(complex_matrix) :: Lowdin

  type(real_matrix)    :: Fock
  type(real_matrix)    :: Cmat 
  type(real_matrix)    :: CS12 
  type(real_matrix)    :: Lowdin


  interface diagonalization 
    module procedure diagonalization_1comp 
    module procedure diagonalization_4comp 
  end interface

  interface mo_occupation 
    module procedure mo_occupation_1comp 
    module procedure mo_occupation_4comp 
  end interface

  interface diis 
    module procedure diis_1comp 
   !module procedure diis_4comp 
  end interface

CONTAINS
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE kramers_unrestricted_scf() 

    !> local
    integer :: i
    integer :: j
    integer :: k
    integer :: ier
    real(8) :: rab 
    integer :: itrscf           = 0
    real(8) :: TrGP             = 0.0d0
    real(8) :: Enuc             = 0.0d0
    real(8) :: Etotal           = 0.0d0
    real(8) :: diis_error       = 0.0d0
    real(8) :: Etotal_last      = 0.0d0
    real(8) :: sum_orb_energies = 0.0d0 
    logical :: if_scf_converged = .false.

    external :: interest_eri    


    !> interface a data 
       call module_interface()
    !> end
    

    !> initialize InteRest integral package 
       call interest_initialize()
    !> end


    !> initialize Lowdin matrix used for symmetric AO-orthonormalization 
       call matrix_allocate( Lowdin )
       call process_overlap_integrals( S12 = Lowdin )
       !> test association
          !call test_association( Lowdin, 'Lowdin' )
       !> end
    !> end


    !> initialize MO expansion coefficients
       call matrix_allocate( Cmat )
    !> end


    !> initialize transformed MO coefficients: C'=S^{1/2}*C 
       call matrix_allocate( CS12 )
    !> end


    !> allocate Fock matrix 
       call matrix_allocate( Fock )
    !> end


    !> calculate the nuclear repulsion energy 
       Enuc = 0.0d0
       do i=1,size(atom)
         do j=(i+1),size(atom),1
           rab  = (atom(i)%coordinate_x - atom(j)%coordinate_x)**2
           rab  = (atom(i)%coordinate_y - atom(j)%coordinate_y)**2 + rab 
           rab  = (atom(i)%coordinate_z - atom(j)%coordinate_z)**2 + rab 
           Enuc = Enuc + atom(i)%charge*atom(j)%charge/dsqrt(rab) 
         enddo
       enddo
       print *,' Nuclear repulsion energy = ',Enuc
    !> end

   
    !> SCF cycles
       do 
         itrscf = itrscf + 1


         !> process AO-density matrix
            !> set MOs occupation  
               call mo_occupation( Cmat )
            !> end
            !> calculate AO-density matrix, P = C*dag(C)
               call matrix_multiplication( A=Cmat, opb='c', B=Cmat, C=Fock )
            !> end
         !> end
         
         
         !> two-electron Fock matrix contribution G (a direct mode)
            call two_electron_integrals( interest_eri, Fock, TrGP )
         !> end
         
         
         !> add one-electron Fock matrix contributions
            !> add kinetic energy contribution
               call process_overlap_integrals( Fock = Fock ) 
            !> end
            !> add nuclear attraction contribution
               call process_nuclear_attraction_integrals( Fock = Fock ) 
            !> end
            !> test association
               !*call test_association( Fock, 'Fock' )
            !> end
         !> end
         
         
         !> transform the Fock matrix to the orthonormal AO-basis (Lowdin basis)
         !> F' = S^{-1/2}*F*S^{-1/2}
            call matrix_multiplication( A=Fock,   B=Lowdin, C=Cmat )
            call matrix_multiplication( A=Lowdin, B=Cmat  , C=Fock )
            !> test: association
               !*call test_association( Fock, 'Fock_before_diagonalization' )
            !> end
         !> end
         

         !> diis procedure
            call diis( Fock, CS12, Cmat, itrscf, diis_error )
         !> end


         !> diagonalization
            call diagonalization( Fock, Cmat, CS12, Lowdin, sum_orb_energies )
         !> end
         

         !> scf cycle and convergence criteria report 
            Etotal  = sum_orb_energies - 0.5d0*TrGP + Enuc
            write(6,'(2x,a,i3,f20.10,2(3x,d8.2))') 'Cycle: ',itrscf,Etotal,abs(Etotal-Etotal_last),diis_error
            !write(6,*)

            if( itrscf ==  scf%max_scf_iterations )exit
            if( diis_error <  scf%diis_error_treshold .and. abs(Etotal-Etotal_last)  <  scf%energy_change_treshold )exit

            Etotal_last = Etotal
         !> end
         

       enddo
    !> end


    call matrix_deallocate( Cmat   )
    call matrix_deallocate( CS12   )
    call matrix_deallocate( Fock   )
    call matrix_deallocate( Lowdin )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE diis_1comp( Fock, CS12, Work, itrscf, diis_error ) 

    !> input/output 
       type(real_matrix), intent(inout) :: Fock
       type(real_matrix), intent(inout) :: CS12
       type(real_matrix), intent(inout) :: Work
    !> end

    !> input
       integer, intent(in) :: itrscf
    !> end

    !> output
       real(8), intent(out) :: diis_error
    !> end

    !> local
       !> a maximal dimension of B-matrix
       integer, parameter :: max_diis_dimension = 20

       logical :: do_diis    = .true.
       logical :: do_damping = .false.

       real(8)              :: ddot 
       integer              :: nrow
       integer              :: nshl
       integer              :: info 
       integer              :: itrdiis
       integer              :: i, j, k, ij
       integer              :: record_number 
       integer              :: record_length
       integer              :: record_number_old 
       real(8)              :: errvec_product
       real(8)              :: Bdiis(max_diis_dimension,max_diis_dimension) 
       integer, allocatable :: ipi(:)
       real(8), allocatable :: xxx(:)
       real(8), allocatable :: wrk(:)
       real(8), allocatable :: bij(:,:)

       integer              :: iij
       real(8)              :: BMCOND
       real(8), allocatable :: grc(:)
    !> end




    !> first iteration initialization
       if( itrscf == 1 )then
         !> initiate damping convergence acceleration
            if( do_damping )then
              open ( 16, file='Fock_matrix', form='unformatted' )
              write( 16 ) Fock%elms(:,:)  
              close( 16 )
            endif
         !> initiate DIIS convergence acceleration
            if( do_diis )then 
              Bdiis(1,1) = 0.0d0
              Bdiis(1,2:max_diis_dimension) = -1.0d0
            endif
            return
       endif
    !> end

    !> damping convergence acceleration
       if( do_damping )then
         open ( 16, file='Fock_matrix', form='unformatted' )
         read ( 16 ) Work%elms(:,:)  
                     Fock%elms(:,:) = (1.0d0 - scf%damping_factor)*Fock%elms(:,:)
                     Fock%elms(:,:) = Fock%elms(:,:) + scf%damping_factor*Work%elms(:,:)
         write( 16 ) Fock%elms(:,:)  
         close( 16 )
       endif
    !> end


    !> calculate diis error vector = F'*P' - dag(F'*P') and
    !> report the maximum absolute value within the occupied-virtual block

       nrow =       Fock%nrow
       nshl = size( Fock%elms, 2 )

       !> set MOs occupation for C' = S^{1/2} * C 
          call mo_occupation( CS12 )
       !> end
       
       !> calculate the density matrix (in the Lowdin basis), P' = C' * dag(C')
          call matrix_multiplication( A=CS12, opb='c', B=CS12, C=Work )
       !> end
       
       !> calculate F' * P' (virtual-occupied block)
       !> where F' = S^{-1/2} * F^{ao} * S^{-1/2}
          call matrix_multiplication( A=Fock, B=Work, C=CS12 )
       !> end
       
       !> anti-Hermitian matrix
          diis_error = 0.0d0 
          do k=1,nshl
            ij = 0
            do j=1,nrow 
              do i=j,nrow
                ij = ij + 1
                Work%elms(ij,k) = CS12%elms(index(i,j,nrow),k) - CS12%elms(index(j,i,nrow),k)
              enddo
            enddo
            !> find the maximum element in the error vector 
               do i=1,ij
                 if( abs(Work%elms(i,k)) > diis_error ) diis_error = abs(Work%elms(i,k)) 
               enddo
            !> end
          enddo
       !> end
    !> end
 

    if( .not. do_diis )return

    itrdiis       = min ( itrscf,   max_diis_dimension     ) 
    record_number = mod ( itrscf-2, max_diis_dimension - 1 ) + 1
    record_length = nshl*8*nrow*(nrow+1)/2

    open( 17, file='diis_error_vectors_alpha', access='direct', recl=record_length, form='unformatted' )
    open( 18, file='diis_error_vectors_beta',  access='direct', recl=record_length, form='unformatted' )
    open( 61, file='diis_fock_matrices_alpha', access='direct', recl=record_length, form='unformatted' )
    open( 62, file='diis_fock_matrices_beta',  access='direct', recl=record_length, form='unformatted' )

    !> shift the upper triangle of B-matrix in the case of max_diis_dimension
       if( itrdiis == max_diis_dimension )then
         do j=3,itrdiis
           do i=3,j
             Bdiis(i-1,j-1) = Bdiis(i,j) 
           enddo
         enddo
       endif
    !> end

    !> read diis error vectors from disk and calculate upper triangle of B-matrix 
       Bdiis(2:itrdiis,itrdiis) = 0.0d0
       do k=1,nshl
         do i=1,(itrdiis-2),1
           record_number_old = mod( (itrscf-itrdiis-1+i), (max_diis_dimension-1) ) + 1
           if( k==1 ) read( 17, rec=record_number_old ) CS12%elms(1:ij,k)
           if( k==2 ) read( 18, rec=record_number_old ) CS12%elms(1:ij,k)
           Bdiis((i+1),itrdiis) = Bdiis((i+1),itrdiis) + 2.0d0*ddot( ij, Work%elms(:,k), 1, CS12%elms(:,k), 1 ) 
         enddo
         Bdiis(itrdiis,itrdiis) = Bdiis(itrdiis,itrdiis) + 2.0d0*ddot( ij, Work%elms(:,k), 1, Work%elms(:,k), 1 ) 
       enddo 
    !> end

    !> write the current diis error vector and Fock matrix on disk
       do k=1,nshl
         if( k==1 ) write( 17, rec=record_number )Work%elms(1:ij,1)
         if( k==2 ) write( 18, rec=record_number )Work%elms(1:ij,2)
         ij = 0
         do j=1,nrow 
           do i=j,nrow
             ij = ij + 1
             Work%elms(ij,k) = Fock%elms(index(i,j,nrow),k)
           enddo
         enddo
         if( k==1 ) write( 61, rec=record_number )Work%elms(1:ij,k)
         if( k==2 ) write( 62, rec=record_number )Work%elms(1:ij,k)
       enddo
    !> end

    !> exit if there are not enough matrices for DIIS interpolation
       if( itrscf < 4 )then
         close( 17, status='keep' )
         close( 18, status='keep' )
         close( 61, status='keep' )
         close( 62, status='keep' )
         return
       endif
    !> end

    !> turn off damping convergence acceleration
       do_damping = .false.
    !> end

    !> solve the DIIS equation (using Lapack dsysv routine) and construct an diis-averaged Fock matrix
    !> note that Bdiis is upper-triangle stored
       allocate( ipi(itrdiis)         )
       allocate( xxx(itrdiis)         )
       allocate( wrk(itrdiis)         )
       allocate( bij(itrdiis,itrdiis) )

         xxx(1)         = -1.0d0
         xxx(2:itrdiis) =  0.0d0 
         bij            =  Bdiis(1:itrdiis,1:itrdiis)
         !> xxx = on input B, on output X for A*X = B 
            call dsysv( 'U', itrdiis, 1, bij,itrdiis, & 
                                         ipi,         &
                                         xxx,itrdiis, &
                                         wrk,itrdiis, &
                                         info         ) 
            if( info > 0 )then
               write(*,*)'Bdiis matrix is singular => the solution could not be computed and DIIS would not be applied!'
               stop 
            end if
         !> end
    !> end

    !> calculate an averaged Fock matrix
       do k=1,nshl
         Work%elms(1:ij,k) = Work%elms(1:ij,k)*xxx(itrdiis) 
         do i=1,(itrdiis-2),1
           record_number_old = mod( (itrscf-itrdiis-1+i), (max_diis_dimension-1) ) + 1
           if( k==1 ) read( 61, rec=record_number_old ) Fock%elms(1:ij,k)
           if( k==2 ) read( 62, rec=record_number_old ) Fock%elms(1:ij,k)
           Work%elms(1:ij,k) = Work%elms(1:ij,k) + Fock%elms(1:ij,k)*xxx(i+1) 
         enddo
       enddo
    !> end

    !> make full Fock matrix 
       do k=1,nshl
         ij = 0
         do j=1,nrow 
           do i=j,nrow
             ij = ij + 1
             Fock%elms(index(i,j,nrow),k) = Work%elms(ij,k)
           enddo
         enddo
         do j=1,nrow 
           do i=j,nrow
             Fock%elms(index(j,i,nrow),k) = Fock%elms(index(i,j,nrow),k) 
           enddo
         enddo
       enddo
    !> end

    deallocate( ipi )
    deallocate( xxx  )
    deallocate( wrk )
    deallocate( bij  )

    close( 17, status='keep' )
    close( 18, status='keep' )
    close( 61, status='keep' )
    close( 62, status='keep' )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE diagonalization_1comp( Fock, Cmat, CS12, Lowdin, sum_orb_energies ) 

    !> input/output 
       type(real_matrix) :: Fock
       type(real_matrix) :: Cmat
       type(real_matrix) :: CS12
       type(real_matrix) :: Lowdin
       real(8)           :: sum_orb_energies
    !> end

    !> local
       integer              :: i
       integer              :: j
       integer              :: ier
       integer              :: idw
       integer              :: nrow
       real(8), allocatable :: Work(:)
       real(8), allocatable :: Eval(:,:)

       type(real_matrix)    :: Tmp
       type(real_matrix)    :: Smat 
    !> end


    !> initialization
       nrow = Fock%nrow
       idw  = max( 1, 3*nrow-1 )
       allocate( Work( idw ) )
       allocate( Eval( nrow, size(Fock%elms,2) ) )
    !> end

   
    !> diagonalization using lapack routines
       do i=1,size(Fock%elms,2)
          call dsyev( 'v', 'l', nrow, Fock%elms(1,i), nrow, Eval(1,i), Work, idw, ier )
          if( ier /= 0 )then
            write(6,'(2x,a)') 'Error in the diagonalization procedure DSYEV => stop'
            stop 
          endif
          !> test electronic eigenvalues
            !if( i==1 )then
            !  write(6,*)' alpha orbital energies: '
            !  write(6,*) Eval(:,1)
            !else
            !  write(6,*)' beta  orbital energies: '
            !  write(6,*) Eval(:,2)
            !endif
          !> end
       enddo
    !> end


    !> calculate the sum of orbital energies of the occupied MOs  
       sum_orb_energies = 0.0d0
       do i=1,size(Fock%elms,2)
         do j=1,nocc(i)
           sum_orb_energies = sum_orb_energies + Eval(j,i)
         enddo
       enddo
    !> end

    !> test MO-density matrix; dag(C')*C' = 1_occ
       !*call matrix_allocate( Tmp, nrow )
       !*call matrix_multiplication( opa='c', A=Fock, opb='n', B=Fock, C=Tmp  )
       !*call test_association( Tmp, 'MO density matrix' )
       !*call matrix_deallocate( Tmp )
       !*stop
    !> end

    !> backup MO coefficients 
       CS12%nrow = Fock%nrow
       CS12%ncol = Fock%ncol
       CS12%elms = Fock%elms
    !> end

    !> back-transformation of MO coefficients to non-orthonormal AOs; C = S^{-1/2}*C'
       call matrix_multiplication( A=Lowdin, B=Fock, C=Cmat )
    !> end

    !> test the number of electrons = Tr(P*S) 
      !call matrix_allocate( Tmp, nrow  )
      !call matrix_allocate( Smat, nrow )
      !!> construct AO-density matrix; P = C*dag(C) 
      !   call matrix_multiplication( A=Cmat, opb='c', B=Cmat, C=Tmp )
      !!> end
      !!> number_of_electrons = Tr(P*S)
      !   call process_overlap_integrals( S = Smat ) 
      !   write(6,*) ' Number of electrons = ',matrix_dotproduct( A=Tmp, B=Smat )
      !!> end
      !call matrix_deallocate( Tmp  )
      !call matrix_deallocate( Smat )
    !> end

    deallocate( Eval ) 
    deallocate( Work ) 
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE diagonalization_4comp( Evec, Pmat, Lowdin ) 

    !> input/output 
       type(complex_matrix) :: Evec
       type(complex_matrix) :: Pmat
       type(complex_matrix) :: Lowdin
    !> end

    !> local
       integer              :: i
       integer              :: ier
       integer              :: idw
       integer              :: nrow
       type(complex_matrix) :: Smat 
       real(8), allocatable :: Work(:,:)
       real(8), allocatable :: Eval(:,:)
    !> end


    allocate( Eval( naoc(1), 4  ) )
    allocate( Work( naoc(1), 12 ) )

    !> forward-reordering of Fock matrix elements prior to diagonalisation (on output Pmat)
       call matrix_elements_reordering( 'f', Evec%elms, Pmat%elms, naoc(1) ) 
    !> end

    !> diagonalization using lapack routines
       call zheev('v', 'u', naoc(4), Pmat%elms, naoc(4), Eval, Evec%elms, naoc(4)*naoc(4), Work, ier)
       if( ier /= 0 )then
         write(6,'(2x,a)') 'Error in the diagonalization procedure ZHEEV => stop'
         stop 
       endif
    !> end

    !> test electronic eigenvalues
       write(6,*) Eval( :, 3:4 )
    !> end

    !> set MO occupation  
    !> fixme: 10 je docasny pocet elektronov, zmen!!!
       !old call mo_occupation( Pmat%elms, naoc(4), 10 )  
       call mo_occupation( Pmat )  
    !> end

    !> backward-reordering of eigenvector matrix elements (on output Evec)
       call matrix_elements_reordering( 'b', Evec%elms, Pmat%elms, naoc(1) ) 
    !> end

    !> reconsider the blocks association 
       call deassociate( A=Evec, full=.true. )
       call   associate( A=Evec, abll=.true., &
                                 absl=.true., &
                                 abls=.true., &
                                 abss=.true., &
                                 bbll=.true., &
                                 bbsl=.true., &
                                 bbls=.true., &
                                 bbss=.true.  )
    !> end

    !> test: MO-density matrix; dag(C')*C' = 1_occ
    !  allocate( TMP%elms( naoc(1)*naoc(1), 4, 4 ) )
    !  call matrix_multiplication( opa='c', A=Fock, opb='n', B=Fock, C=TMP  )
    !  call test_association( TMP, 'MO density matrix' )
    !  write(6,*) Evec%elms(index(1,1),3,3)
    !  deallocate( TMP%elms )
    !  stop
    !> end

    !> back-transformation of MO coefficients to non-orthonormal AOs; C = S^{-1/2}*C'
       call matrix_multiplication( A=Lowdin, B=Evec, C=Pmat )
    !> end

    !> construct AO-density matrix; P = C*dag(C) 
       call matrix_multiplication( A=Pmat, opb='c', B=Pmat, C=Evec  )
    !> end

    !> test number_of_electrons = Tr(P*S)
       call matrix_allocate( Smat )
       call process_overlap_integrals( S = Smat ) 
       write(6,*) ' Number of electrons = ',matrix_dotproduct( A=Evec, B=Smat )
       call matrix_deallocate( Smat )
    !> end

    deallocate( Eval ) 
    deallocate( Work ) 

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE matrix_elements_reordering( direction, A, B, n ) 

    !> input 
    integer,          intent(in) :: n 
    character(len=1), intent(in) :: direction

    !> forward:  A-input  B-output
    !> backward: A-output B-input
    complex(8) :: A(n,n,4,4)
    complex(8) :: B(n,4,n,4)

    !> local
    integer :: i
    integer :: j
    integer :: k


    select case( direction )
    !> forward reordering 
    case('f')
      do k=1,4
        do j=1,n
          do i=1,4
            B(:,i,j,k) = A(:,j,i,k) 
          enddo
        enddo
      enddo
    !> backward reordering 
    case('b')
      do k=1,4
        do j=1,n
          do i=1,4
            A(:,j,i,k) = B(:,i,j,k) 
          enddo
        enddo
      enddo
    case default
      write(6,'(2x,a)') 'Unknown reordering required => stop' 
      stop
    end select
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE mo_occupation_1comp( C ) 

    !> input 
       type(real_matrix), intent(inout) :: C
    !> end

    !> local
       integer :: min, max
    !> end

    
    max = C%nrow*C%ncol

    !> alpha MOs occupation
       min = C%nrow*nocc(1) + 1
       C%elms( min:max, 1 ) = 0.0d0
    !> end

    if( size(C%elms,2) == 1 )return

    !> beta MOs occupation
       min = C%nrow*nocc(2) + 1
       C%elms( min:max, 2 ) = 0.0d0
    !> end
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>
  SUBROUTINE mo_occupation_4comp( C ) 

    !> input 
       type(complex_matrix), intent(inout) :: C
    !> end

    !> local
       integer :: i, imin, imax
    !> end



    !> fixme: fix this routine
        write(6,*) 'Warrning: mo_occupation_4comp not finished yet ... => stop'   
        stop
    !> end
    
    !> alpha MOs occupation
       !code here ....
    !> end

    !> beta MOs occupation
       !code here ....
    !> end

    !> unoccupying the positronic states
      !old do j=1,(n/2)
      !old   do i=1,n
      !old     C(i,j) = dcmplx(0.0d0,0.0d0)
      !old   enddo
      !old enddo
    !> end

    !> occupying the electronic states 
      !old do j=((n/2)+num_electrons+1),n
      !old   do i=1,n
      !old     C(i,j) = dcmplx(0.0d0,0.0d0)
      !old   enddo
      !old enddo
    !> end
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE associate( A, aall, aals, aasl, aass, &
                           abll, abls, absl, abss, &
                           ball, bals, basl, bass, & 
                           bbll, bbls, bbsl, bbss, &
                           full                    )
    !> input
    type(complex_matrix), target, intent(inout) :: A
    logical, optional,            intent(in)    :: full 
    logical, optional,            intent(in)    :: aall, aals, aasl, aass 
    logical, optional,            intent(in)    :: abll, abls, absl, abss
    logical, optional,            intent(in)    :: ball, bals, basl, bass
    logical, optional,            intent(in)    :: bbll, bbls, bbsl, bbss

    !> local
    integer :: i
    integer :: j


    if( present(full) )then
      do j=1,4
        do i=1,4
          A%info(i,j)%p => A%elms(1,i,j)
        enddo
      enddo
      return
    endif

    if( present(aall) ) A%info(1,1)%p => A%elms(1,1,1)
    if( present(aals) ) A%info(1,2)%p => A%elms(1,1,2)
    if( present(aasl) ) A%info(2,1)%p => A%elms(1,2,1)
    if( present(aass) ) A%info(2,2)%p => A%elms(1,2,2)

    if( present(abll) ) A%info(1,3)%p => A%elms(1,1,3)
    if( present(abls) ) A%info(1,4)%p => A%elms(1,1,4)
    if( present(absl) ) A%info(2,3)%p => A%elms(1,2,3)
    if( present(abss) ) A%info(2,4)%p => A%elms(1,2,4)

    if( present(ball) ) A%info(3,1)%p => A%elms(1,3,1)
    if( present(bals) ) A%info(3,2)%p => A%elms(1,3,2)
    if( present(basl) ) A%info(4,1)%p => A%elms(1,4,1)
    if( present(bass) ) A%info(4,2)%p => A%elms(1,4,2)

    if( present(bbll) ) A%info(3,3)%p => A%elms(1,3,3)
    if( present(bbls) ) A%info(3,4)%p => A%elms(1,3,4)
    if( present(bbsl) ) A%info(4,3)%p => A%elms(1,4,3)
    if( present(bbss) ) A%info(4,4)%p => A%elms(1,4,4)

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE deassociate( A, aall, aals, aasl, aass, &
                             abll, abls, absl, abss, &
                             ball, bals, basl, bass, & 
                             bbll, bbls, bbsl, bbss, &
                             full                    )
    !> input
    type(complex_matrix), target, intent(inout) :: A
    logical, optional,            intent(in)    :: full 
    logical, optional,            intent(in)    :: aall, aals, aasl, aass 
    logical, optional,            intent(in)    :: abll, abls, absl, abss
    logical, optional,            intent(in)    :: ball, bals, basl, bass
    logical, optional,            intent(in)    :: bbll, bbls, bbsl, bbss

    !> local
    integer :: i
    integer :: j


    if( present(full) )then
      do j=1,4
        do i=1,4
          A%info(i,j)%p => null() 
        enddo
      enddo
      return
    endif

    if( present(aall) ) A%info(1,1)%p => null() 
    if( present(aals) ) A%info(1,2)%p => null() 
    if( present(aasl) ) A%info(2,1)%p => null() 
    if( present(aass) ) A%info(2,2)%p => null() 

    if( present(abll) ) A%info(1,3)%p => null() 
    if( present(abls) ) A%info(1,4)%p => null() 
    if( present(absl) ) A%info(2,3)%p => null() 
    if( present(abss) ) A%info(2,4)%p => null() 

    if( present(ball) ) A%info(3,1)%p => null() 
    if( present(bals) ) A%info(3,2)%p => null() 
    if( present(basl) ) A%info(4,1)%p => null() 
    if( present(bass) ) A%info(4,2)%p => null() 

    if( present(bbll) ) A%info(3,3)%p => null() 
    if( present(bbls) ) A%info(3,4)%p => null() 
    if( present(bbsl) ) A%info(4,3)%p => null() 
    if( present(bbss) ) A%info(4,4)%p => null() 

  END SUBROUTINE
! ------------------------------------------------------------------------------------
END MODULE
