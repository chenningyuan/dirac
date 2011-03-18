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
! Licensing:    dirac_copyright_start
!                     Copyright (c) 2010 by the authors of DIRAC.
!                     All Rights Reserved.
!                     
!                     This source code is part of the DIRAC program package.
!                     It is provided under a written license and may be used,
!                     copied, transmitted, or stored only in accordance to the
!                     conditions of that written license.
!                     
!                     In particular, no part of the source code or compiled modules may
!                     be distributed outside the research group of the license holder.
!                     This means also that persons (e.g. post-docs) leaving the research
!                     group of the license holder may not take any part of Dirac,
!                     including modified files, with him/her, unless that person has
!                     obtained his/her own license.
!                     
!                     For information on how to get a license, as well as the
!                     author list and the complete list of contributors to the
!                     DIRAC program, see: http://dirac.chem.vu.nl
!               dirac_copyright_end
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


  !> naos(:) - number of atomic orbitals, spherical coordinates
  !> naoc(:) - number of atomic orbitals, cartesian coordinates 
  integer :: naos(1:4)
  integer :: naoc(1:4)


  !> fixme: temporary 5 alpha-electrons and 5 beta-electrons
  integer :: nocc(2) = (/5,5/)

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
    real(8) :: Etotal_last      = 0.0d0
    real(8) :: diis_error(2)    = 0.0d0
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
       call matrix_allocate( Lowdin, naoc(1) )
       call process_overlap_integrals( S12 = Lowdin )
       !> test association
          !call test_association( Lowdin, 'Lowdin' )
       !> end
    !> end


    !> initialize MO expansion coefficients
       call matrix_allocate( Cmat, naoc(1) )
    !> end


    !> initialize transformed MO coefficients: C'=S^{1/2}*C 
       call matrix_allocate( CS12, naoc(1) )
    !> end


    !> allocate Fock matrix 
       call matrix_allocate( Fock, naoc(1) )
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
            write(6,'(2x,a,i3,f20.10,2(3x,d8.2))') 'Cycle: ',itrscf,Etotal,abs(Etotal-Etotal_last),diis_error(1)

            if( itrscf ==  scf%max_scf_iterations )exit
            if( diis_error(1) <  scf%diis_error_treshold .and. abs(Etotal-Etotal_last)  <  scf%energy_change_treshold )exit

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
       real(8), intent(out) :: diis_error(2) 
    !> end

    !> local
       !> on which scf iteration diis error vector starts to be accumulated
       integer, parameter :: diis_start_on_itrscf = 2
       !> a minimal dimension of B-matrix to start diis 
       integer, parameter :: min_dimension = 2
       !> a maximal dimension of B-matrix
       integer, parameter :: max_dimension = 20

       real(8)              :: xmax
       real(8)              :: ddot 
       integer              :: nrow
       integer              :: nshl
       integer              :: info 
       integer              :: itrdiis = 0
       integer              :: i, j, k, ij
       real(8)              :: errvec_product
       real(8)              :: Bdiis(max_dimension,max_dimension,2) 
       character*9          :: flfock 
       integer, allocatable :: ipi(:)
       real(8), allocatable :: xxx(:)
       real(8), allocatable :: wrk(:)
       real(8), allocatable :: bij(:,:)
    !> end


    
    if( itrscf < diis_start_on_itrscf )return

    nrow = Fock%nrow
    nshl = size(Fock%elms,2)

    !> set MOs occupation for C' = S^{1/2} * C 
       call mo_occupation( CS12 )
    !> end

    !> calculate density matrix (in the Lowdin basis), P' = C' * dag(C')
       call matrix_multiplication( A=CS12, opb='c', B=CS12, C=Work )
    !> end

    !> calculate F' * P' (virtual-occupied block)
    !> where F' = S^{-1/2} * F^{ao} * S^{-1/2}
       call matrix_multiplication( A=Fock, B=Work, C=CS12 )
    !> end

       itrdiis = itrdiis + 1

    !> calculate diis error vector = F'*P' - dag(F'*P')
    !> the error vector has components of anti-Hermitian matrix 
       do k=1,nshl
         ij = 0
         do j=1,nrow 
           do i=j,nrow
             ij = ij + 1
             Work%elms(ij,k) = CS12%elms(index(i,j,nrow),k) - CS12%elms(index(j,i,nrow),k)
           enddo
         enddo

         !> find the maximum element in the error vector 
            diis_error(k) = 0.0d0 
            do i=1,ij
              xmax = Work%elms(i,k)
              xmax = dsqrt(xmax*xmax)
              if( xmax > diis_error(k) ) diis_error(k) = xmax
            enddo
         !> end

         !> write diis error vector on disk
            write(flfock,'(a7,i2.2)')'FL1DIIS',itrdiis 
            open(17,file=flfock,form='unformatted')
            write(17)Work%elms(1:ij,k)
            close(17,status='keep')
         !> end

         !> add a new column/row to the upper triangle of B-matrix
            if( itrdiis == 1 )then
              Bdiis(1,1,k) = 0.0d0
              do i=2,max_dimension
                Bdiis(1,i,k) = -1.0d0
              enddo
            else 
              do i=2,itrdiis
                write(flfock,'(a7,i2.2)')'FL1DIIS',(i-1)
                open(17,file=flfock,form='unformatted')
                read(17)CS12%elms(1:ij,k)
                close(17,status='keep')
                errvec_product = 2.0d0*ddot( ij, Work%elms(:,k), 1, CS12%elms(:,k), 1 )
                Bdiis(i,(itrdiis+1),k) = errvec_product 
              enddo
            endif
            Bdiis((itrdiis+1),(itrdiis+1),k) = 2.0d0*ddot( ij, Work%elms(:,k), 1, Work%elms(:,k), 1 )
         !> end

         !> write the current Fock matrix on disk 
            ij = 0
            do j=1,nrow 
              do i=j,nrow              !fixme: triangle????
                ij = ij + 1
                Work%elms(ij,k) = Fock%elms(index(i,j,nrow),k)
              enddo
            enddo
            write(flfock,'(a7,i2.2)')'FLFOCKA',itrdiis 
            open(61,file=flfock,form='unformatted')
            write(61)Work%elms(1:ij,k)
            close(61,status='keep')
         !> end
       enddo
    !> end


    !> exit if not having enough matrices for DIIS interpolation
       if( itrdiis < (diis_start_on_itrscf + min_dimension - 1) )return 
    !> end


    !> solve the DIIS equation (using Lapack dsysv routine) and construct an diis-averaged Fock matrix
    !> note that Bdiis is lower-triangle stored
       itrdiis = itrdiis + 1
       allocate( ipi(itrdiis)         )
       allocate( xxx(itrdiis)         )
       allocate( wrk(itrdiis)         )
       allocate( bij(itrdiis,itrdiis) )

       do k=1,nshl
         xxx(1)         = -1.0d0
         xxx(2:itrdiis) =  0.0d0 
         bij            =  Bdiis(1:itrdiis,1:itrdiis,k)
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

         !> calculate an averaged Fock matrix
            Work%elms(1:ij,k) = Work%elms(1:ij,k)*xxx(itrdiis) 
            do i=2,(itrdiis-1),1
              write(flfock,'(a7,i2.2)')'FLFOCKA',(i-1)
              open(61,file=flfock,form='unformatted')
              read(61)Fock%elms(1:ij,k)
              close(61,status='keep')
              Work%elms(1:ij,k) = Work%elms(1:ij,k) + Fock%elms(1:ij,k)*xxx(i)
            enddo
         !> end
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


    itrdiis = itrdiis - 1
    deallocate( ipi )
    deallocate( xxx  )
    deallocate( wrk )
    deallocate( bij  )

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
       call matrix_allocate( Smat, naoc(1) )
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
!> Interface to DIRAC
  SUBROUTINE module_interface()

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"
#include "dcbdhf.h"
#include "nucdata.h"

    !> initialization of local variables
    integer :: i
    integer :: j
    integer :: ier
    type( type_gto ) :: tmpgto
  
  
    !> test if the large component basis is uncontracted
    if( nplrg.ne.nlarge )then
      write(6,*) ' Error: contracted basis set found => stop!'
      write(6,*) ' Note:  uncontracted basis sets are only supported for kramers_unrestricted_scf_solver'
      stop
    endif
  

    !> interface scf settings 
    !scf = type_scf( maxitr, mxdiis, scfcnv(1), scfcnv(1) )
     scf = type_scf( maxitr, mxdiis, 1.d-9, 1.d-9 )

  
    !> interface nuclear data
    do i=1,109 
      element(i) = type_element( 'not_set', 0.0d0, nucmass(i), nuclabel(i) )
    enddo
  
  
    !> interface molecular data
    allocate( atom(nucind), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: atom(:)'
    do i=1,size(atom)
      atom(i) = type_atom( charge(i), cord(1,i), cord(2,i), cord(3,i), gnuexp(i) )
    enddo
  
  
    !> interface basis set data
    allocate( gto(nlrgsh), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: gto(:)'
    i = 1
    do while( i <= size(gto) )  
      do j=1,nrco(i)
        gto(i)  = type_gto( i, 0, nhkt(i), ncent(i), 0, 0, priccf(i,j), priexp(i), priccf(i,j) )
        i = i + 1
      enddo
    enddo


    !> reorder basis functions to the angular momentum 
    !> ascending order to increase the integral performance
    do i=1,size(gto)-1
      do j=(i+1),size(gto)
        if( gto(j)%lvalue >= gto(i)%lvalue )cycle
        tmpgto = gto(i)
        gto(i) = gto(j)
        gto(j) = tmpgto
      enddo
    enddo


    !> set some data + AO-vector offsets 
    naos(1) = 0
    naoc(1) = 0
    do i=1,size(gto)
      gto(i)%offset = naoc(1) 
      gto(i)%sdegen = 2*gto(i)%lvalue - 1 
      gto(i)%cdegen = ( gto(i)%lvalue*(gto(i)%lvalue+1) )/2 
      naos(1) = naos(1) + gto(i)%sdegen 
      naoc(1) = naoc(1) + gto(i)%cdegen
    enddo
    do i=2,4
      naos(i) = i*naos(1)
      naoc(i) = i*naoc(1)
    enddo



    do i=1,size(gto)
      write(6,'(1x,5(a,i3,3x),3(a,d15.5,3x))') ' i = ',i,                            &
                                               ' index  = ',gto(i)%index,            &
                                               ' offset = ',gto(i)%offset,           &
                                               ' lvalue = ',gto(i)%lvalue,           &
                                               ' origin = ',gto(i)%origin,           &
                                               ' norm = ',  gto(i)%norm(1),          &
                                               ' exponent = ',gto(i)%exponent(1),    &
                                               ' coefficient = ',gto(i)%coefficient(1)
    enddo

  
!!! write(6,*) 
!!! write(6,'(2x,a,i5)') 'Basis set informations:'
!!! write(6,'(2x,a,i5)') '-----------------------'
!!! write(6,*) 
!!! write(6,*) 
!!! write(6,'(2x,a)')    'label    atom   charge   prim   cont   L-basis   RKB     S-basis'
!!! write(6,'(2x,70a1)') ('-',i=1,70)
!!!
!!! write(6,'(2x,70a1)') ('-',i=1,70)
    write(6,*) 
    write(6,'(2x,a,i5  )') 'Total number of basis function shells:    ',size(gto)
    write(6,'(2x,a,i5,a)') 'Total number of spherical basis functions:',naos(1),' (used for calculation)' 
    write(6,'(2x,a,i5  )') 'Total number of cartesian basis functions:',naoc(1)
    write(6,*) 
  
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
