! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Module:       module_interest_interface.f90 
!
! Description:  - interfacing basic data between Dirac & InteRest programs 
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

MODULE MODULE_INTEREST_INTERFACE

    use module_interest_eri
    use module_interest_peri

#ifdef USE_MPI_MOD_F90
  use mpi
#endif
  implicit none

#if defined (VAR_MPI) && !defined(USE_MPI_MOD_F90)
#include "mpif.h"
#endif

  public transfer_gto_info
  public initialize_interest 
  public initialize_interest_herd 
  public process_overlap_integrals_2comp 
  public process_overlap_integrals_4comp 
  public process_kinetic_energy_integrals_2comp 
  public process_kinetic_energy_integrals_4comp 
  public process_nuclear_attraction_integrals_2comp
  public process_nuclear_attraction_integrals_4comp

  public process_two_electron_integrals_2comp
  public process_two_electron_integrals_4comp
  public process_two_electron_integrals_2comp_slaves
  public process_two_electron_integrals_4comp_slaves

  private

  !> matrix dimensions
     !> note: by definition nrow, ncol refer to the number of cartesian GTOs
     integer, save :: nshl = 1
     integer, save :: nrow = 0
     integer, save :: ncol = 0
     integer, save :: nr_shells = 0 
     real(8), save :: cvalue
     real(8), save :: screening_threshold = 1.d-20
  !> end

  !> definition of "constant" data type 
     type type_constant
          real(8) :: alpha       !=(1/c)
          real(8) :: alpha2      !=(1/c)^2 
     end  type
     type(type_constant) :: constant 
  !> end

  !> definition of "symmetry" data type 
     !todo: later implementation 
     type type_symmetry
          integer :: dummy 
     end  type
     type(type_symmetry) :: symmetry 
  !> end

  !> definition of "atom" data type 
     type type_atom
          real(8) :: charge 
          real(8) :: coordinate_x
          real(8) :: coordinate_y
          real(8) :: coordinate_z
          real(8) :: gnu_exponent 
     end  type
     type(type_atom), allocatable :: atom(:)
  !> end

  !> definition of "gto" data type
     type type_gto
          integer :: index
          integer :: offset
          integer :: lvalue
          integer :: origin
          integer :: sdegen 
          integer :: cdegen 
          real(8) :: exponent(1)
          real(8) :: coefficient(1)
     end  type
     type(type_gto), allocatable, target :: gto(:)
  !> end

  integer, save :: myrank  = 0
  integer, save :: nr_cpus = 1 

CONTAINS



! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE initialize_interest()

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"
#include "dcbdhf.h"
#include "nucdata.h"
#include "dcbgen.h"
#include "cbihr2.h"

    !> local variables
       integer :: i
       integer :: j
       integer :: ij
       integer :: ier
       type( type_gto ) :: tmpgto
    !> end
  

    !> test if the large component basis is uncontracted
       if( nplrg.ne.nlarge )then
         write(6,*) ' Error: contracted basis set found => stop!'
         write(6,*) ' Note: .MDIRAC does not support contracted basis sets'
         stop
       endif
    !> end
  
    !> interface constants (cval == speed of light used in Dirac [might be different from 137])
       constant = type_constant( (1.0d0/cval), (1.0d0/cval)**2 )
    !> end

    !> interface molecular data
       allocate( atom(nucind), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: atom(:)'
       do i=1,size(atom)
         atom(i) = type_atom( charge(i), &
                              cord(1,i), &
                              cord(2,i), &
                              cord(3,i), &
                              gnuexp(i)  )
       enddo
    !> end
  
    !> interface basis set data
    !> fixme: contracted basis sets
       allocate( gto(nlrgsh), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: gto(:)'
       i =0
       ij=0
       do while( i < size(gto) )  
         do j=1,nrco(i+1)
           i=i+1
           gto(i) = type_gto( ij,                      &
                              ij,                      &
                              nhkt(i),                 &
                              ncent(i),                & 
                              (2*nhkt(i)-1),           &
                              (nhkt(i)*(nhkt(i)+1))/2, &
                              priexp(i),               &
                              priccf(i,j)              )
           !> note: cartesian indexation
           ij=ij+(nhkt(i)*(nhkt(i)+1))/2
         enddo
       enddo
    !> end

    !> temporary print ... 
       do i=1,size(gto)
         write(6,'(1x,5(a,i3,3x),2(a,d15.5,3x))') ' i = ',i,                             &
                                                  ' index  = ',     gto(i)%index,        &
                                                  ' offset = ',     gto(i)%offset,       &
                                                  ' lvalue = ',     gto(i)%lvalue,       &
                                                  ' origin = ',     gto(i)%origin,       &
                                                  ' exponent = ',   gto(i)%exponent(1),  &
                                                  ' coefficient = ',gto(i)%coefficient(1)
       enddo

       !> initialize InteRest integral package 
       call interest_initialize(.true.)

       write(6,*) 
       write(6,'(2x,a,i5  )') 'Total number of basis function shells:    ',size(gto)
       write(6,'(2x,a,i5  )') 'Total number of spherical basis functions:',sum(gto(:)%sdegen)
       write(6,'(2x,a,i5,a)') 'Total number of cartesian basis functions:',sum(gto(:)%cdegen),' (used in InteRest calculations)'
       write(6,*) 
    !> end

    !> global module variables
    nshl      = 1
    nrow      = sum( gto(:)%cdegen)
    ncol      = sum( gto(:)%cdegen)
    nr_shells = size(gto)
    cvalue    = cval
    if( scrfck > 0.0d0 ) screening_threshold = scrfck


#ifdef VAR_MPI
    call dirac_parctl( 8 ) 
#endif

    call initialize_interest_herd

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE initialize_interest_herd()

#ifdef VAR_MPI
    !> local MPI related variables
    integer :: i, j, ij, ier, nr_atoms
    integer(kind=mpi_address_kind) :: disp(9)                      !mpi_address_kind is always = 8
    integer                        :: blocklen(9), typex(9)
    integer                        :: newtype
  

    call MPI_comm_rank( MPI_comm_world, myrank,  ier )
    call MPI_comm_size( MPI_comm_world, nr_cpus, ier )

    !> initialize InteRest integral package 
    call interest_initialize(.false.)

    !> broadcast constants 
       !> define the local absolute address of each element within a derived data type
       call mpi_bcast(constant%alpha, 1, mpi_real8,0,mpi_comm_world,ier)
       call mpi_bcast(constant%alpha2,1, mpi_real8,0,mpi_comm_world,ier)
    !> end

    !> broadcast atom 
       if( allocated(atom) ) nr_atoms = size(atom)
       call mpi_bcast(nr_atoms,1,mpi_integer,0,mpi_comm_world,ier)
       if( .not. allocated(atom) ) allocate( atom(nr_atoms), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: atom(:)'
       !> define the local absolute address of each element within a derived data type
       do i=1,nr_atoms
        call MPI_bcast( atom(i)%charge,       1, MPI_real8, 0, MPI_comm_world, ier)
        call MPI_bcast( atom(i)%coordinate_x, 1, MPI_real8, 0, MPI_comm_world, ier)
        call MPI_bcast( atom(i)%coordinate_y, 1, MPI_real8, 0, MPI_comm_world, ier)
        call MPI_bcast( atom(i)%coordinate_z, 1, MPI_real8, 0, MPI_comm_world, ier)
        call MPI_bcast( atom(i)%gnu_exponent, 1, MPI_real8, 0, MPI_comm_world, ier)
       end do
    !> end

    !> broadcast gto 
       if( allocated(gto) ) nr_shells = size(gto)
       call mpi_bcast(nr_shells,1,mpi_integer,0,mpi_comm_world,ier)
       if( .not. allocated(gto) ) allocate( gto(nr_shells), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: gto(:)'
       do i=1,nr_shells
        call MPI_bcast( gto(i)%index,          1, MPI_integer, 0, MPI_comm_world, ier)
        call MPI_bcast( gto(i)%offset,         1, MPI_integer, 0, MPI_comm_world, ier)
        call MPI_bcast( gto(i)%lvalue,         1, MPI_integer, 0, MPI_comm_world, ier)
        call MPI_bcast( gto(i)%origin,         1, MPI_integer, 0, MPI_comm_world, ier)
        call MPI_bcast( gto(i)%sdegen,         1, MPI_integer, 0, MPI_comm_world, ier)
        call MPI_bcast( gto(i)%cdegen,         1, MPI_integer, 0, MPI_comm_world, ier)
        call MPI_bcast( gto(i)%exponent(1),    1, MPI_real8,   0, MPI_comm_world, ier)
        call MPI_bcast( gto(i)%coefficient(1), 1, MPI_real8,   0, MPI_comm_world, ier)
       end do
    !> end

    !> final broadcasting (global module variables)
       call mpi_bcast(nshl,               1,mpi_integer,0,mpi_comm_world,ier)
       call mpi_bcast(nrow,               1,mpi_integer,0,mpi_comm_world,ier)
       call mpi_bcast(ncol,               1,mpi_integer,0,mpi_comm_world,ier)
       call mpi_bcast(cvalue,             1,mpi_real8,  0,mpi_comm_world,ier)
       call mpi_bcast(screening_threshold,1,mpi_real8,  0,mpi_comm_world,ier)
    !> end
#endif

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_overlap_integrals_2comp( S ) 

    !> input
    real(8) :: S(nrow,ncol) 

    !> local
    integer              :: i, j 
    real(8), allocatable :: Smat(:,:,:)

    !> external InteRest integral routines
    external :: interest_overlap 


    allocate( Smat( nrow,ncol,2 ) ); Smat = 0.0d0
    call two_center_integrals( interest_overlap, Smat, nrow )
    call symmetrize_matrix( A = Smat(1,1,1), ida = nrow )
    S = 0.0d0
    do j=1,ncol
      do i=1,nrow
        S(i,j) = Smat(i,j,1) 
      enddo
    enddo
    deallocate( Smat )

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_overlap_integrals_4comp( S, ss_scaling ) 

    !> input
    real(8)           :: S(2*nrow,2*ncol) 
    real(8), optional :: ss_scaling

    !> local
    integer              :: i, j 
    real(8), allocatable :: Smat(:,:,:)

    !> external InteRest integral routines
    external :: interest_overlap 


    allocate( Smat( nrow,ncol,2 ) ); Smat = 0.0d0
    call two_center_integrals( interest_overlap, Smat, nrow )
    call symmetrize_matrix(                                 A = Smat(1,1,1), ida = nrow )  !Re(A)
    call symmetrize_matrix( factor = constant%alpha2/4.0d0, A = Smat(1,1,2), ida = nrow )  !Re(A) = (xx + yy + zz)
    S = 0.0d0
    if(present(ss_scaling))then
      do j=1,ncol
        do i=1,nrow
          S(i     ,j     ) = Smat(i,j,1) 
          S(i+nrow,j+ncol) = Smat(i,j,2) * ss_scaling
        enddo
      enddo
    else
      do j=1,ncol
        do i=1,nrow
          S(i     ,j     ) = Smat(i,j,1) 
          S(i+nrow,j+ncol) = Smat(i,j,2) 
        enddo
      enddo
    end if
    deallocate( Smat )

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_kinetic_energy_integrals_2comp( T ) 

    !> input
    !> fixme: the case when nz!=4, since we always assume nz=4
    real(8) :: T(nrow,ncol,4) 

    !> local
    integer              :: i, j 
    real(8), allocatable :: Tmat(:,:,:)

    !> external InteRest integral routines
    external :: interest_overlap 


    allocate( Tmat( nrow,ncol,2 ) ); Tmat = 0.0d0
    call two_center_integrals( interest_overlap, Tmat, nrow )
    call symmetrize_matrix( factor = 0.5d0, A = Tmat(1,1,2), ida = nrow )  !Re(A) = (xx + yy + zz)
    do j=1,ncol
      do i=1,nrow
        T(i,j,1) = T(i,j,1) + Tmat(i,j,2)
      enddo
    enddo
    deallocate( Tmat )

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_kinetic_energy_integrals_4comp( T ) 

    !> input
    !> fixme: the case when nz!=4, since we always assume nz=4
    real(8) :: T(2*nrow,2*ncol,4) 

    !> local
    integer              :: i, j 
    real(8), allocatable :: Tmat(:,:,:)

    !> external InteRest integral routines
    external :: interest_overlap 


    allocate( Tmat( nrow,ncol,2 ) ); Tmat = 0.0d0
    call two_center_integrals( interest_overlap, Tmat, nrow )
    call symmetrize_matrix( factor = 0.5d0, A = Tmat(1,1,2), ida = nrow )  !Re(A) = (xx + yy + zz)
    do j=1,ncol
      do i=1,nrow
        T(i+nrow,j     ,1) = T(i+nrow,j     ,1) + Tmat(i,j,2)
        T(i     ,j+ncol,1) = T(i     ,j+ncol,1) + Tmat(i,j,2)
        T(i+nrow,j+ncol,1) = T(i+nrow,j+ncol,1) - Tmat(i,j,2)
      enddo
    enddo
    deallocate( Tmat )

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE process_nuclear_attraction_integrals_2comp( V ) 
 
    !> input
    !> fixme: the case when nz!=4, since we always assume nz=4
    real(8) :: V(nrow,ncol,4)

    !> local
    integer              :: i, j 
    real(8), allocatable :: Vmat(:,:,:)

    !> external InteRest integral routines
    external :: interest_nuclear_attraction_point_nucleus 


    allocate( Vmat( nrow,ncol,5 ) ); Vmat = 0.0d0
    call two_center_integrals( interest_nuclear_attraction_point_nucleus, Vmat, nrow )
    call symmetrize_matrix   ( A = Vmat(1,1,1), ida = nrow )
    do j=1,ncol
      do i=1,nrow
        V(i,j,1) = V(i,j,1) + Vmat(i,j,1)
      enddo
    enddo
    deallocate( Vmat )

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE process_nuclear_attraction_integrals_4comp( V ) 
 
    !> input
    !> fixme: the case when nz!=4, since we always assume nz=4
    real(8) :: V(2*nrow,2*ncol,4)

    !> local
    integer              :: i, j 
    real(8), allocatable :: Vmat(:,:,:)

    !> external InteRest integral routines
    external :: interest_nuclear_attraction_point_nucleus 


    allocate( Vmat( nrow,ncol,5 ) ); Vmat = 0.0d0
    call two_center_integrals( interest_nuclear_attraction_point_nucleus, Vmat, nrow )
    call symmetrize_matrix    (                                 A = Vmat(1,1,1), ida = nrow )  !Re(A)
    call symmetrize_matrix    ( factor = constant%alpha2/4.0d0, A = Vmat(1,1,2), ida = nrow )  !Re(A) = (xx + yy + zz)
    call antisymmetrize_matrix( factor = constant%alpha2/4.0d0, A = Vmat(1,1,3), ida = nrow )  !Im(A) = (xy - yx)      
    call antisymmetrize_matrix( factor = constant%alpha2/4.0d0, A = Vmat(1,1,4), ida = nrow )  !Re(B) = (zx - xz)      
    call antisymmetrize_matrix( factor = constant%alpha2/4.0d0, A = Vmat(1,1,5), ida = nrow )  !Im(B) = (yz - zy)     
    do j=1,ncol
      do i=1,nrow
        V(i     ,j     ,1) = V(i     ,j     ,1) + Vmat(i,j,1)
        V(i+nrow,j+ncol,1) = V(i+nrow,j+ncol,1) + Vmat(i,j,2)
        V(i+nrow,j+ncol,2) = V(i+nrow,j+ncol,2) + Vmat(i,j,3)
        V(i+nrow,j+ncol,3) = V(i+nrow,j+ncol,3) + Vmat(i,j,4)
        V(i+nrow,j+ncol,4) = V(i+nrow,j+ncol,4) + Vmat(i,j,5)
      enddo
    enddo
    deallocate( Vmat )

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE two_center_integrals( interest_routine, final_integrals, ida ) 

    !> input
       interface 
         subroutine interest_routine(rkb,                    &
                                     fint,na,nb,nrkb,        &
                                     la,alpha,ax,ay,az,anorm,&
                                     lb,beta, bx,by,bz,bnorm,&
                                     ncentr,atomic_data      )
       
           integer, intent(out) :: na 
           integer, intent(out) :: nb 
           logical, intent(in)  :: rkb
           integer, intent(out) :: nrkb
           integer, intent(in)  :: la,lb
           integer, intent(in)  :: ncentr
           real(8), intent(out) :: fint(*)
           real(8), intent(in)  :: alpha,ax,ay,az,anorm
           real(8), intent(in)  :: beta, bx,by,bz,bnorm
           real(8), intent(in)  :: atomic_data(ncentr,5)
         end subroutine
       end interface
    !> end

    !> input
    integer, intent(in)  :: ida           !todo: better definition
    logical              :: rkb = .true.  !todo: better definition

    !> output
    real(8) :: final_integrals(*)

    !> local
    integer :: nrkb
    integer :: i, ni
    integer :: j, nj
    real(8) :: gout(5*441)  !todo: better definition 


    do j=1,nr_shells
      do i=j,nr_shells
        call interest_routine(rkb,                                &
                              gout, ni, nj, nrkb,                 &
                              !> bra-function
                              gto(i)%lvalue,                      &
                              gto(i)%exponent(1),                 &
                              atom( gto(i)%origin )%coordinate_x, &
                              atom( gto(i)%origin )%coordinate_y, &
                              atom( gto(i)%origin )%coordinate_z, &
                              gto(i)%coefficient(1),              & 
                              !> ket-function
                              gto(j)%lvalue,                      &
                              gto(j)%exponent(1),                 &
                              atom( gto(j)%origin )%coordinate_x, &
                              atom( gto(j)%origin )%coordinate_y, &
                              atom( gto(j)%origin )%coordinate_z, &
                              gto(j)%coefficient(1),              &
                              !> atomic data 
                              size(atom),                         &
                              (/atom(1:size(atom))%charge,        &
                                atom(1:size(atom))%coordinate_x,  &
                                atom(1:size(atom))%coordinate_y,  &
                                atom(1:size(atom))%coordinate_z,  &
                                atom(1:size(atom))%gnu_exponent/) )

        call store_batch(gout, ni, nj,  nrkb,       &
                         final_integrals, ida, ida, &
                         gto(i)%offset,             &
                         gto(j)%offset              )

      enddo
    enddo    

    CONTAINS 

    subroutine store_batch( A, ida, jda, nrkb, B, idb, jdb, ioffset, joffset )
      !> input 
      integer :: nrkb
      integer :: ioffset 
      integer :: joffset 
      integer :: ida, jda 
      integer :: idb, jdb
      real(8) :: A(5,ida,jda)
      !> output
      real(8) :: B(idb,jdb,nrkb)
      !> local
      integer :: i
      integer :: j
      integer :: k
      !>
      do k=1,nrkb
        do j=1,jda
          do i=1,ida
            B(ioffset+i,joffset+j,k) = A(k,i,j) 
          enddo
        enddo
      enddo
    end subroutine

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> note: on input, G and P must be in the cartesian GTOs 
  SUBROUTINE process_two_electron_integrals_2comp( Gmat, Pmat, intflg, hfxfac ) 

    !> input
    !> fixme: the case when nz!=4, since we always assume nz=4
    integer, intent(in) :: intflg
    real(8), intent(in) :: hfxfac
    !> 1-electron density matrix, Dirac format
    real(8), intent(in) :: Pmat(nrow,nrow,4)  

    !> output
    !> 2-electron Fock matrix, Dirac format
    real(8), intent(out) :: Gmat(nrow,nrow,4) 

    !> local
    real(8)                 :: fll
    integer                 :: i, j, ier
    complex(8), allocatable :: dP(:,:,:,:)
    complex(8), allocatable :: dG(:,:,:,:)


    !> check the eri class
    select case( intflg )
      case( 0 );    return 
      case( 1 );
      case default; stop "The desired class of ERI integrals is not implemented in InteRest :-( => stop"
    end select

    !> data forward alignment 
    fll = 0.5d0
    allocate( dP(2,2,nrow,nrow) ); dP = dcmplx( 0.0d0, 0.0d0 ) 
    allocate( dG(2,2,nrow,nrow) ); dG = dcmplx( 0.0d0, 0.0d0 ) 
    do j=1,nrow
      do i=1,nrow
        !> alpha-alpha = "A"
        dP(1,1,i,j) =  dcmplx( Pmat(i,j,1)*fll, Pmat(i,j,2)*fll )
        !> alpha-beta  = "B"
        dP(1,2,i,j) =  dcmplx( Pmat(i,j,3)*fll, Pmat(i,j,4)*fll )
        !> beta-beta   = "A*"
        dP(2,2,i,j) =  dconjg( dP(1,1,i,j) ) 
        !> beta-alpha  = "-B*"
        dP(2,1,i,j) = -dconjg( dP(1,2,i,j) ) 
      enddo
    enddo

#ifdef VAR_MPI
    !> send the AO-density matrix
    call dirac_parctl( 10 ) 
    call mpi_bcast( hfxfac,    1,           mpi_real8,          0, mpi_comm_world, ier )
    call mpi_bcast( dP,        4*nrow*nrow, mpi_double_complex, 0, mpi_comm_world, ier )
#endif

    !> call core InteRest routine (sequential part)
    call two_electron_integrals_2comp_1x( dG, dP, hfxfac )

#ifdef VAR_MPI
    !> collect the final two-electron Fock matrix
    call mpi_reduce( mpi_in_place, dG, 4*nrow*nrow, mpi_double_complex, mpi_sum, 0, mpi_comm_world, ier )
#endif

    !> data backward alignment 
    do j=1,nrow
      do i=1,nrow
        !> quaternion(0) 
        Gmat(i,j,1) = Gmat(i,j,1) + dreal(dG(1,1,i,j))
        !> quaternion(i) 
        Gmat(i,j,2) = Gmat(i,j,2) + dimag(dG(1,1,i,j))
        !> quaternion(j) 
        Gmat(i,j,3) = Gmat(i,j,3) + dreal(dG(1,2,i,j))
        !> quaternion(k) 
        Gmat(i,j,4) = Gmat(i,j,4) + dimag(dG(1,2,i,j))
      enddo
    enddo

    deallocate( dP )
    deallocate( dG )

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_two_electron_integrals_2comp_slaves( ) 

    !> slave local
    integer :: ier, i
    real(8) :: hfxfac
    
    complex(8), allocatable :: dP(:,:,:,:)
    complex(8), allocatable :: dG(:,:,:,:)


#ifdef VAR_MPI
    allocate( dP(2,2,nrow,nrow) ); dP = dcmplx( 0.0d0, 0.0d0 ) 
    allocate( dG(2,2,nrow,nrow) ); dG = dcmplx( 0.0d0, 0.0d0 ) 

    !> receive the AO-density matrix
    call mpi_bcast( hfxfac,    1,           mpi_real8,          0, mpi_comm_world, ier )
    call mpi_bcast( dP,        4*nrow*nrow, mpi_double_complex, 0, mpi_comm_world, ier )

    !> call core InteRest routine (sequential part)
    call two_electron_integrals_2comp_1x( dG, dP, hfxfac )

    !> send the two-electron Fock matrix of a particular slave upwards to the master
    call mpi_reduce( dG, mpi_in_place, 4*nrow*nrow, mpi_double_complex, mpi_sum, 0, mpi_comm_world, ier )

    deallocate( dP )
    deallocate( dG )
#endif

  END SUBROUTINE



! ------------------------------------------------------------------------------------
  SUBROUTINE two_electron_integrals_2comp_1x( dG, dP, hfxfac ) 

    use module_interest_eri

    !> input
    real(8), intent(in) :: hfxfac
    complex(8), intent(in) :: dP(2,2,nrow,nrow)  

    !> output
    complex(8), intent(out) :: dG(2,2,nrow,nrow)

    !> local
    logical :: doK
    real(8) :: fijkl, fll
    integer :: li, ni, oi   
    integer :: lj, nj, oj   
    integer :: lk, nk, ok   
    integer :: ll, nl, ol   
    integer :: nint, lmax
    integer :: i, j, k, l, ibas
    real(8) :: gout(441*441)    !todo
    integer :: nij, nkl, ij
    real(8) :: ei, ci, xi, yi, zi 
    real(8) :: ej, cj, xj, yj, zj 
    real(8) :: ek, ck, xk, yk, zk 
    real(8) :: el, cl, xl, yl, zl 
    integer :: ierr, mytask
    integer :: myrank=0
#ifdef VAR_MPI
    integer, allocatable :: eri_task_list(:) 
    integer              :: cpu, nr_cpus 
#endif

    !> init for DFT
    doK = ( hfxfac /= 0.0d0 )

#ifdef VAR_MPI
    ij = 0
    call MPI_comm_rank( MPI_comm_world, myrank, ierr)
    call MPI_comm_size( MPI_comm_world, nr_cpus,ierr)
    !> creating a simple MPI eri task list
    allocate( eri_task_list(nr_shells*nr_shells), stat=ierr)
    if( ierr.ne.0 )stop ' Error in allocation: eri_task_list(:)'
    ij =1
    cpu=1
    mytask=0
    do i=1,nr_shells
      do j=1,nr_shells
        if( cpu > nr_cpus )cpu=1
        eri_task_list(ij)=cpu-1
        ij =ij +1
        cpu=cpu+1
      enddo
    enddo
#endif
    ij = 0

 iloop:do i=1,nr_shells
         li =       gto(i)%lvalue
         oi =       gto(i)%offset
         ni =       gto(i)%cdegen
         ei =       gto(i)%exponent(1)
         xi = atom( gto(i)%origin )%coordinate_x
         yi = atom( gto(i)%origin )%coordinate_y
         zi = atom( gto(i)%origin )%coordinate_z
         ci =       gto(i)%coefficient(1)
     
 jloop:  do j=1,nr_shells
           lj =       gto(j)%lvalue
           oj =       gto(j)%offset
           nj =       gto(j)%cdegen
           ej =       gto(j)%exponent(1)
           xj = atom( gto(j)%origin )%coordinate_x
           yj = atom( gto(j)%origin )%coordinate_y
           zj = atom( gto(j)%origin )%coordinate_z
           cj =       gto(j)%coefficient(1)

#ifdef VAR_MPI
          mytask=mytask+1
          if( myrank /= eri_task_list(mytask) )cycle
#endif
       
 kloop:    do k=1,nr_shells
             lk =       gto(k)%lvalue
             ok =       gto(k)%offset
             nk =       gto(k)%cdegen
             ek =       gto(k)%exponent(1)
             xk = atom( gto(k)%origin )%coordinate_x
             yk = atom( gto(k)%origin )%coordinate_y
             zk = atom( gto(k)%origin )%coordinate_z
             ck =       gto(k)%coefficient(1)
       
 lloop:      do l=1,nr_shells
               ll   =       gto(l)%lvalue
               nl   =       gto(l)%cdegen
               ol   =       gto(l)%offset
               el   =       gto(l)%exponent(1)
               xl   = atom( gto(l)%origin )%coordinate_x
               yl   = atom( gto(l)%origin )%coordinate_y
               zl   = atom( gto(l)%origin )%coordinate_z
               cl   =       gto(l)%coefficient(1)
               fijkl= 1.0d0 
       
               call interest_eri(fijkl,gout,       &
                                 li,ei,xi,yi,zi,ci,&
                                 lj,ej,xj,yj,zj,cj,&
                                 lk,ek,xk,yk,zk,ck,&
                                 ll,el,xl,yl,zl,cl )
               call process_dG(ni,nj,nk,nl,oi,oj,ok,ol,  &
                               nrow,doK,hfxfac,gout,dP,dG)
             enddo lloop
           enddo kloop
         enddo jloop
       enddo iloop
#ifdef VAR_MPI
       deallocate( eri_task_list )
#endif
    !> end

    CONTAINS

    subroutine process_dG(ic,jc,kc,lc,io,jo,ko,lo,nbas,doK,scaleK,gout,pmat,gmat) 
      !> input
      logical,    intent(in) :: doK
      integer,    intent(in) :: nbas
      real(8),    intent(in) :: scaleK
      integer,    intent(in) :: ic, jc, kc, lc
      integer,    intent(in) :: io, jo, ko, lo
      complex(8), intent(in) :: pmat(2,2,nbas,nbas)
      !> input/output
      real(8), intent(inout) :: gout(kc,lc,ic,jc)
      !> output
      complex(8), intent(out) :: gmat(2,2,nbas,nbas)
      !> local
      integer    :: i, j, k, l
      integer    :: ibas, jbas, kbas, lbas
      complex(8) :: gg, pll, paa, pab, pba, pbb
      !>
      !> coulomb part
      do j=1,jc
        jbas=jo+j
        do i=1,ic
          ibas=io+i
          do l=1,lc
            lbas=lo+l
            do k=1,kc
              kbas = ko+k
              gg   = dcmplx( gout(k,l,i,j), 0.0d0 )
              pll  = pmat(1,1,lbas,kbas) + pmat(2,2,lbas,kbas)
              gmat(1,1,ibas,jbas) = gmat(1,1,ibas,jbas) + gg*pll
              gmat(2,2,ibas,jbas) = gmat(2,2,ibas,jbas) + gg*pll
            enddo
          enddo
        enddo
      enddo
      !> exchange part
      if( .not. doK )return
      !> pre-multiplication for scalling K-contribution 
      if( scalek /= 1.0d0 )then
        do j=1,jc
          do i=1,ic
            do l=1,lc
              do k=1,kc
                gout(k,l,i,j) = gout(k,l,i,j)*scaleK
              enddo
            enddo
          enddo
        enddo
      endif
      !>
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          do j=1,jc
            jbas=jo+j
            paa = pmat(1,1,jbas,kbas)
            pba = pmat(2,1,jbas,kbas)
            pab = pmat(1,2,jbas,kbas)
            pbb = pmat(2,2,jbas,kbas)
            do i=1,ic
              ibas = io+i
              gg   = dcmplx( gout(k,l,i,j), 0.0d0 )
              gmat(1,1,ibas,lbas) = gmat(1,1,ibas,lbas) - paa*gg
              gmat(1,2,ibas,lbas) = gmat(1,2,ibas,lbas) - pab*gg
              gmat(2,1,ibas,lbas) = gmat(2,1,ibas,lbas) - pba*gg
              gmat(2,2,ibas,lbas) = gmat(2,2,ibas,lbas) - pbb*gg
            enddo
          enddo
        enddo
      enddo
    end subroutine

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> note:  on input, both G and P must be in the cartesian GTOs 
!> fixme: the case when nz!=4, since we always assume nz=4
  SUBROUTINE process_two_electron_integrals_4comp( Gmat, Pmat, intflg, hfxfac ) 

    !> input
    integer, intent(in) :: intflg
    real(8), intent(in) :: hfxfac
    !> 1-electron density matrix, Dirac format
    real(8), intent(in) :: Pmat(2*nrow,2*nrow,4)  

    !> output
    !> 2-electron Fock matrix, Dirac format
    real(8), intent(out) :: Gmat(2*nrow,2*nrow,4) 

    !> local
    integer                 :: i, j, ier
    real(8)                 :: fll, fsl, fss 
    character(len=4)        :: eri_class
    complex(8), allocatable :: dP(:,:,:,:)
    complex(8), allocatable :: dG(:,:,:,:)


    !> init for eri classes
    select case( intflg )
      case( 0 );    return 
      case( 1 );    eri_class = 'llll'
      case( 3 );    eri_class = 'ssll'
      case( 7 );    eri_class = 'ssss'
      case default; stop "The desired class of ERI integrals is not implemented in InteRest :-( => stop"
    end select

    !> data forward alignment 
    fll = 0.5d0
    fsl = constant%alpha /4.0d0 
    fss = constant%alpha2/8.0d0 
    allocate( dP(4,4,nrow,nrow) ); dP = dcmplx( 0.0d0, 0.0d0 ) 
    allocate( dG(4,4,nrow,nrow) ); dG = dcmplx( 0.0d0, 0.0d0 ) 
    do j=1,nrow
      do i=1,nrow
        !> alpha-alpha = "A"
        dP(1,1,i,j) = dcmplx( Pmat(i     ,j     ,1)*fll, Pmat(i     ,j     ,2)*fll ) ! LL
        dP(2,1,i,j) = dcmplx( Pmat(i+nrow,j     ,1)*fsl, Pmat(i+nrow,j     ,2)*fsl ) ! SL
        dP(1,2,i,j) = dcmplx( Pmat(i     ,j+nrow,1)*fsl, Pmat(i     ,j+nrow,2)*fsl ) ! LS
        dP(2,2,i,j) = dcmplx( Pmat(i+nrow,j+nrow,1)*fss, Pmat(i+nrow,j+nrow,2)*fss ) ! SS
        !> beta-beta  = "A*"
        dP(3,3,i,j) = dconjg( dP(1,1,i,j) )
        dP(4,3,i,j) = dconjg( dP(2,1,i,j) )
        dP(3,4,i,j) = dconjg( dP(1,2,i,j) )
        dP(4,4,i,j) = dconjg( dP(2,2,i,j) )
        !> alpha-beta = "B"
        dP(1,3,i,j) = dcmplx( Pmat(i     ,j     ,3)*fll, Pmat(i     ,j     ,4)*fll ) ! LL
        dP(2,3,i,j) = dcmplx( Pmat(i+nrow,j     ,3)*fsl, Pmat(i+nrow,j     ,4)*fsl ) ! SL
        dP(1,4,i,j) = dcmplx( Pmat(i     ,j+nrow,3)*fsl, Pmat(i     ,j+nrow,4)*fsl ) ! LS
        dP(2,4,i,j) = dcmplx( Pmat(i+nrow,j+nrow,3)*fss, Pmat(i+nrow,j+nrow,4)*fss ) ! SS
        !> beta-alpha =" -B*"
        dP(3,1,i,j) = -dconjg( dP(1,3,i,j) )
        dP(4,1,i,j) = -dconjg( dP(2,3,i,j) )
        dP(3,2,i,j) = -dconjg( dP(1,4,i,j) )
        dP(4,2,i,j) = -dconjg( dP(2,4,i,j) )
      enddo
    enddo

#ifdef VAR_MPI
    !> send the AO-density matrix
    call dirac_parctl( 9 ) 
    call mpi_bcast( eri_class, 4,            mpi_character,      0, mpi_comm_world, ier )
    call mpi_bcast( hfxfac,    1,            mpi_real8,          0, mpi_comm_world, ier )
    call mpi_bcast( dP,        16*nrow*nrow, mpi_double_complex, 0, mpi_comm_world, ier )
#endif

    !> call core InteRest routine (sequential part)
!nshl=2
!call two_electron_integrals_4comp_8x( dG, dP, eri_class, hfxfac )
!nshl=1
    call two_electron_integrals_4comp_1x( dG, dP, eri_class, hfxfac )

#ifdef VAR_MPI
    !> collect the final two-electron Fock matrix
    call mpi_reduce( mpi_in_place, dG, 16*nrow*nrow, mpi_double_complex, mpi_sum, 0, mpi_comm_world, ier )
#endif

    !> data backward alignment 
    fsl = constant%alpha /2.0d0 
    fss = constant%alpha2/4.0d0 
    do j=1,nrow
      do i=1,nrow
        !> quaternion(0) 
        Gmat(i     ,j     ,1) = Gmat(i     ,j     ,1) + dreal(dG(1,1,i,j))
        Gmat(i+nrow,j     ,1) = Gmat(i+nrow,j     ,1) + dreal(dG(2,1,i,j))*fsl
        Gmat(i     ,j+nrow,1) = Gmat(i     ,j+nrow,1) + dreal(dG(2,1,j,i))*fsl
        Gmat(i+nrow,j+nrow,1) = Gmat(i+nrow,j+nrow,1) + dreal(dG(2,2,i,j))*fss
        !> quaternion(i) 
        Gmat(i     ,j     ,2) = Gmat(i     ,j     ,2) + dimag(dG(1,1,i,j))
        Gmat(i+nrow,j     ,2) = Gmat(i+nrow,j     ,2) + dimag(dG(2,1,i,j))*fsl
        Gmat(i     ,j+nrow,2) = Gmat(i     ,j+nrow,2) - dimag(dG(2,1,j,i))*fsl
        Gmat(i+nrow,j+nrow,2) = Gmat(i+nrow,j+nrow,2) + dimag(dG(2,2,i,j))*fss
        !> quaternion(j) 
        Gmat(i     ,j     ,3) = Gmat(i     ,j     ,3) + dreal(dG(1,3,i,j))
        Gmat(i+nrow,j     ,3) = Gmat(i+nrow,j     ,3) + dreal(dG(2,3,i,j))*fsl
        Gmat(i     ,j+nrow,3) = Gmat(i     ,j+nrow,3) - dreal(dG(2,3,j,i))*fsl
        Gmat(i+nrow,j+nrow,3) = Gmat(i+nrow,j+nrow,3) + dreal(dG(2,4,i,j))*fss
        !> quaternion(k) 
        Gmat(i     ,j     ,4) = Gmat(i     ,j     ,4) + dimag(dG(1,3,i,j))
        Gmat(i+nrow,j     ,4) = Gmat(i+nrow,j     ,4) + dimag(dG(2,3,i,j))*fsl
        Gmat(i     ,j+nrow,4) = Gmat(i     ,j+nrow,4) - dimag(dG(2,3,j,i))*fsl
        Gmat(i+nrow,j+nrow,4) = Gmat(i+nrow,j+nrow,4) + dimag(dG(2,4,i,j))*fss
      enddo
    enddo

    deallocate( dP )
    deallocate( dG )

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE process_two_electron_integrals_4comp_slaves( ) 

    !> slave local
    integer :: ier, i
    real(8) :: hfxfac
    character(len=4) :: eri_class
    
    complex(8), allocatable :: dP(:,:,:,:)
    complex(8), allocatable :: dG(:,:,:,:)


#ifdef VAR_MPI
    allocate( dP(4,4,nrow,nrow) ); dP = dcmplx( 0.0d0, 0.0d0 ) 
    allocate( dG(4,4,nrow,nrow) ); dG = dcmplx( 0.0d0, 0.0d0 ) 

    !> receive the AO-density matrix
    call mpi_bcast( eri_class, 4,            mpi_character,      0, mpi_comm_world, ier )
    call mpi_bcast( hfxfac,    1,            mpi_real8,          0, mpi_comm_world, ier )
    call mpi_bcast( dP,        16*nrow*nrow, mpi_double_complex, 0, mpi_comm_world, ier )

    !> call core InteRest routine (sequential part)
!nshl=2
!call two_electron_integrals_4comp_8x( dG, dP, eri_class, hfxfac )
!nshl=1
    call two_electron_integrals_4comp_1x( dG, dP, eri_class, hfxfac )

    !> send the two-electron Fock matrix of a particular slave upwards to the master
    call mpi_reduce( dG, mpi_in_place, 16*nrow*nrow, mpi_double_complex, mpi_sum, 0, mpi_comm_world, ier )

    deallocate( dP )
    deallocate( dG )
#endif

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>
  SUBROUTINE two_electron_integrals_4comp_1x( dG, dP, eri_class, hfxfac ) 

    use module_interest_eri

    !> input
    real(8),          intent(in) :: hfxfac
    character(len=4), intent(in) :: eri_class
    complex(8),       intent(in) :: dP(4,4,nrow,nrow)  

    !> output
    complex(8), intent(out) :: dG(4,4,nrow,nrow) 

    !> local
    logical :: doK
    real(8) :: fijkl
    integer :: li, ni, oi   
    integer :: lj, nj, oj   
    integer :: lk, nk, ok   
    integer :: ll, nl, ol   
    integer :: nint, lmax
    integer :: i, j, k, l, ibas
    real(8) :: gout(441*441*25)    !todo
    integer :: nij, nkl, ij
    real(8) :: ei, ci, xi, yi, zi 
    real(8) :: ej, cj, xj, yj, zj 
    real(8) :: ek, ck, xk, yk, zk 
    real(8) :: el, cl, xl, yl, zl 
    integer :: ierr, mytask
    integer :: myrank=0
#ifdef VAR_MPI
    integer, allocatable :: eri_task_list(:) 
    integer              :: cpu, nr_cpus 
#endif

    !> init for DFT
    doK = ( hfxfac /= 0.0d0 )

#ifdef VAR_MPI
    ij = 0
    call MPI_comm_rank( MPI_comm_world, myrank, ierr)
    call MPI_comm_size( MPI_comm_world, nr_cpus,ierr)
    !> creating a simple MPI eri task list
    allocate( eri_task_list(nr_shells*nr_shells), stat=ierr)
    if( ierr.ne.0 )stop ' Error in allocation: eri_task_list(:)'
    ij =1
    cpu=1
    mytask=0
    do i=1,nr_shells
      do j=1,nr_shells
        if( cpu > nr_cpus )cpu=1
        eri_task_list(ij)=cpu-1
        ij =ij +1
        cpu=cpu+1
      enddo
    enddo
#endif
    ij = 0

 iloop:do i=1,nr_shells
         li =       gto(i)%lvalue
         oi =       gto(i)%offset
         ni =       gto(i)%cdegen
         ei =       gto(i)%exponent(1)
         xi = atom( gto(i)%origin )%coordinate_x
         yi = atom( gto(i)%origin )%coordinate_y
         zi = atom( gto(i)%origin )%coordinate_z
         ci =       gto(i)%coefficient(1)
     
 jloop:  do j=1,nr_shells
           lj =       gto(j)%lvalue
           oj =       gto(j)%offset
           nj =       gto(j)%cdegen
           ej =       gto(j)%exponent(1)
           xj = atom( gto(j)%origin )%coordinate_x
           yj = atom( gto(j)%origin )%coordinate_y
           zj = atom( gto(j)%origin )%coordinate_z
           cj =       gto(j)%coefficient(1)

#ifdef VAR_MPI
          mytask=mytask+1
          if( myrank /= eri_task_list(mytask) )cycle
#endif

       
 kloop:    do k=1,nr_shells
             lk =       gto(k)%lvalue
             ok =       gto(k)%offset
             nk =       gto(k)%cdegen
             ek =       gto(k)%exponent(1)
             xk = atom( gto(k)%origin )%coordinate_x
             yk = atom( gto(k)%origin )%coordinate_y
             zk = atom( gto(k)%origin )%coordinate_z
             ck =       gto(k)%coefficient(1)
       
 lloop:      do l=1,nr_shells
               ll   =       gto(l)%lvalue
               nl   =       gto(l)%cdegen
               ol   =       gto(l)%offset
               el   =       gto(l)%exponent(1)
               xl   = atom( gto(l)%origin )%coordinate_x
               yl   = atom( gto(l)%origin )%coordinate_y
               zl   = atom( gto(l)%origin )%coordinate_z
               cl   =       gto(l)%coefficient(1)
               fijkl= 1.0d0 
       
               SELECT CASE( eri_class )
                 CASE( 'llll' )
                     call interest_eri_rkb('llll',fijkl,gout,&
                                           li,ei,xi,yi,zi,ci,&
                                           lj,ej,xj,yj,zj,cj,&
                                           lk,ek,xk,yk,zk,ck,&
                                           ll,el,xl,yl,zl,cl )
                     call process_dG_llll(ni,nj,nk,nl,oi,oj,ok,ol,  &
                                          nrow,doK,hfxfac,gout,dP,dG)
                 CASE( 'ssll' )
                     call interest_eri_rkb('ssll',fijkl,gout,&
                                           li,ei,xi,yi,zi,ci,&
                                           lj,ej,xj,yj,zj,cj,&
                                           lk,ek,xk,yk,zk,ck,&
                                           ll,el,xl,yl,zl,cl )
                     call process_dG_ssll(ni,nj,nk,nl,oi,oj,ok,ol,  &
                                          nrow,doK,hfxfac,gout,dP,dG)
                 CASE( 'ssss' ) 
                     call interest_eri_rkb('ssss',fijkl,gout,&
                                           li,ei,xi,yi,zi,ci,&
                                           lj,ej,xj,yj,zj,cj,&
                                           lk,ek,xk,yk,zk,ck,&
                                           ll,el,xl,yl,zl,cl )
                     call process_dG_ssss(ni,nj,nk,nl,oi,oj,ok,ol,  &
                                          nrow,doK,hfxfac,gout,dP,dG)
                 CASE DEFAULT; stop 'Error: this ERI class is not implemented yet => stop'
               END SELECT
             enddo lloop
           enddo kloop
         enddo jloop
       enddo iloop
#ifdef VAR_MPI
       deallocate( eri_task_list )
#endif
    !> end

    CONTAINS

    subroutine process_dG_llll(ic,jc,kc,lc,io,jo,ko,lo,nbas,doK,scaleK,gout,pmat,gmat) 
      !> input
      logical,    intent(in) :: doK
      integer,    intent(in) :: nbas
      real(8),    intent(in) :: scaleK
      integer,    intent(in) :: ic, jc, kc, lc
      integer,    intent(in) :: io, jo, ko, lo
      complex(8), intent(in) :: pmat(4,4,nbas,nbas)
      !> input/output
      real(8), intent(inout) :: gout(5,kc,lc,5,ic,jc)
      !> output
      complex(8), intent(out) :: gmat(4,4,nbas,nbas)
      !> local
      integer    :: i, j, k, l
      integer    :: ibas, jbas, kbas, lbas
      complex(8) :: gg, pll, paa_ll, pab_ll, pba_ll, pbb_ll
      !>
      !> coulomb part
      do j=1,jc
        jbas=jo+j
        do i=1,ic
          ibas=io+i
          do l=1,lc
            lbas=lo+l
            do k=1,kc
              kbas = ko+k
              gg   = dcmplx( gout(1,k,l,1,i,j), 0.0d0 )
              pll  = pmat(1,1,lbas,kbas) + pmat(3,3,lbas,kbas)
              gg   = gg*pll
              gmat(1,1,ibas,jbas) = gmat(1,1,ibas,jbas) +        gg
              gmat(3,3,ibas,jbas) = gmat(3,3,ibas,jbas) + dconjg(gg)
            enddo
          enddo
        enddo
      enddo
      !> exchange part
      if( .not. doK )return
      !> pre-multiplication for scalling K-contribution 
      if( scalek /= 1.0d0 )then
        do j=1,jc
          do i=1,ic
            do l=1,lc
              do k=1,kc
                gout(1,k,l,1,i,j) = gout(1,k,l,1,i,j)*scaleK
              enddo
            enddo
          enddo
        enddo
      endif
      !>
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          do j=1,jc
            jbas=jo+j
            paa_ll = pmat(1,1,jbas,kbas)
            pba_ll = pmat(3,1,jbas,kbas)
            pab_ll = pmat(1,3,jbas,kbas)
            pbb_ll = pmat(3,3,jbas,kbas)
            do i=1,ic
              ibas = io+i
              gg   = dcmplx( gout(1,k,l,1,i,j), 0.0d0 )
              gmat(1,1,ibas,lbas) = gmat(1,1,ibas,lbas) - paa_ll*gg
              gmat(1,3,ibas,lbas) = gmat(1,3,ibas,lbas) - pab_ll*gg
              gmat(3,1,ibas,lbas) = gmat(3,1,ibas,lbas) - pba_ll*gg
              gmat(3,3,ibas,lbas) = gmat(3,3,ibas,lbas) - pbb_ll*gg
            enddo
          enddo
        enddo
      enddo
    end subroutine
    !-------------
    subroutine process_dG_ssll(ic,jc,kc,lc,io,jo,ko,lo,nbas,doK,scaleK,gout,pmat,gmat) 
      !> input
      logical,    intent(in) :: doK
      integer,    intent(in) :: nbas
      real(8),    intent(in) :: scaleK
      integer,    intent(in) :: ic, jc, kc, lc
      integer,    intent(in) :: io, jo, ko, lo
      complex(8), intent(in) :: pmat(4,4,nbas,nbas)
      !> input/output
      real(8), intent(inout) :: gout(5,kc,lc,5,ic,jc)
      !> output
      complex(8), intent(out) :: gmat(4,4,nbas,nbas)
      !> local
      integer    :: i, j, k, l, m
      integer    :: ibas, jbas, kbas, lbas
      complex(8) :: pll, paa, pab, pba, pbb
      complex(8) :: gg, gaag, gabg, gbag, gbbg
      complex(8) :: paa_ll, pab_ll, pba_ll, pbb_ll
      complex(8) :: paa_sl, pab_sl, pba_sl, pbb_sl
      !>
      !> coulomb part
      do j=1,jc
        jbas=jo+j
        do i=1,ic
          ibas=io+i
          do l=1,lc
            lbas=lo+l
            do k=1,kc
              kbas=ko+k
              !> integrals processing: [ij|kl]
                 !> large-large block 
                 gg = dcmplx( gout(1,k,l,1,i,j), 0.0d0 )
                 !> small-large block 
                 gaag = dcmplx( gout(1,k,l,2,i,j), gout(1,k,l,3,i,j) )
                 gabg = dcmplx( gout(1,k,l,4,i,j), gout(1,k,l,5,i,j) )
                 gbbg = dconjg( gaag )
                 gbag =-dconjg( gabg ) 
              !> end
              !> contraction [kl|ij]*ji
                 paa = pmat(2,2,jbas,ibas)
                 pba = pmat(4,2,jbas,ibas)
                 pab = pmat(2,4,jbas,ibas)
                 pbb = pmat(4,4,jbas,ibas)
                 pll = pmat(1,1,jbas,ibas) + pmat(3,3,jbas,ibas)
                 !>
                 gg = gg*pll + gaag*paa + gabg*pba + gbag*pab + gbbg*pbb
                 !>
                 gmat(1,1,kbas,lbas) = gmat(1,1,kbas,lbas) +        gg
                 gmat(3,3,kbas,lbas) = gmat(3,3,kbas,lbas) + dconjg(gg)
              !> contraction [ij|kl]*lk
                 pll = pmat(1,1,lbas,kbas) + pmat(3,3,lbas,kbas)
                 !>
                 gaag = gaag*pll
                 gabg = gabg*pll
                 !>
                 gmat(2,2,ibas,jbas) = gmat(2,2,ibas,jbas) +        gaag
                 gmat(4,2,ibas,jbas) = gmat(4,2,ibas,jbas) - dconjg(gabg) !+ gbag
                 gmat(2,4,ibas,jbas) = gmat(2,4,ibas,jbas) +        gabg
                 gmat(4,4,ibas,jbas) = gmat(4,4,ibas,jbas) + dconjg(gaag) !+ gbbg
              !> end
            enddo
          enddo
        enddo
      enddo
      !> exchange part 
      if( .not. doK )return
      !> pre-multiplication for scalling K-contribution 
      if( scalek /= 1.0d0 )then
        do j=1,jc
          do i=1,ic
            do m=1,5
              do l=1,lc
                do k=1,kc
                  gout(1,k,l,m,i,j) = gout(1,k,l,m,i,j)*scaleK
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
      !>
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          do j=1,jc
            jbas=jo+j
            !> large-large block
            paa_ll = pmat(1,1,jbas,kbas)
            pba_ll = pmat(3,1,jbas,kbas)
            pab_ll = pmat(1,3,jbas,kbas)
            pbb_ll = pmat(3,3,jbas,kbas)
            !> small-large block
            paa_sl = pmat(2,1,jbas,kbas)
            pba_sl = pmat(4,1,jbas,kbas)
            pab_sl = pmat(2,3,jbas,kbas)
            pbb_sl = pmat(4,3,jbas,kbas)
            do i=1,ic
              ibas=io+i
              !> large-large block 
              gg = dcmplx( gout(1,k,l,1,i,j), 0.0d0 )
              !> small-large block 
              gaag = dcmplx( gout(1,k,l,2,i,j), gout(1,k,l,3,i,j) )
              gabg = dcmplx( gout(1,k,l,4,i,j), gout(1,k,l,5,i,j) )
              gbbg = dconjg( gaag )
              gbag =-dconjg( gabg ) 
              !> large-large block
              gmat(1,1,ibas,lbas) = gmat(1,1,ibas,lbas) - paa_ll*gg
              gmat(1,3,ibas,lbas) = gmat(1,3,ibas,lbas) - pab_ll*gg
              gmat(3,1,ibas,lbas) = gmat(3,1,ibas,lbas) - pba_ll*gg
              gmat(3,3,ibas,lbas) = gmat(3,3,ibas,lbas) - pbb_ll*gg
              !> small-large block
              gmat(2,1,ibas,lbas) = gmat(2,1,ibas,lbas) - gaag*paa_sl - gabg*pba_sl 
              gmat(2,3,ibas,lbas) = gmat(2,3,ibas,lbas) - gaag*pab_sl - gabg*pbb_sl 
              gmat(4,1,ibas,lbas) = gmat(4,1,ibas,lbas) - gbag*paa_sl - gbbg*pba_sl 
              gmat(4,3,ibas,lbas) = gmat(4,3,ibas,lbas) - gbag*pab_sl - gbbg*pbb_sl 
            enddo
          enddo
        enddo
      enddo
    end subroutine
    !-------------
    subroutine process_dG_ssss(ic,jc,kc,lc,io,jo,ko,lo,nbas,doK,scaleK,gout,pmat,gmat) 
      !> input
      logical,    intent(in) :: doK
      integer,    intent(in) :: nbas
      real(8),    intent(in) :: scaleK
      integer,    intent(in) :: ic, jc, kc, lc
      integer,    intent(in) :: io, jo, ko, lo
      complex(8), intent(in) :: pmat(4,4,nbas,nbas)
      !> input/output
      real(8), intent(inout) :: gout(5,kc,lc,5,ic,jc)
      !> output
      complex(8), intent(out) :: gmat(4,4,nbas,nbas)
      !> local
      integer    :: i, j, k, l, m
      integer    :: ibas, jbas, kbas, lbas
      complex(8) :: gg
      complex(8) :: gaag,   gabg,   gbag,   gbbg
      complex(8) :: ggaa,   ggab,   ggba,   ggbb   
      complex(8) :: gaagaa, gabgaa, gbagaa, gbbgaa
      complex(8) :: gaagab, gabgab, gbagab, gbbgab
      complex(8) :: gaagba, gabgba, gbagba, gbbgba
      complex(8) :: gaagbb, gabgbb, gbagbb, gbbgbb
      complex(8) :: pll
      complex(8) :: paa,    pab,    pba,    pbb
      complex(8) :: paa_ll, pab_ll, pba_ll, pbb_ll
      complex(8) :: paa_ls, pab_ls, pba_ls, pbb_ls
      complex(8) :: paa_sl, pab_sl, pba_sl, pbb_sl
      complex(8) :: paa_ss, pab_ss, pba_ss, pbb_ss
      !> 
      !> coulomb part
      do j=1,jc
        jbas=jo+j
        do i=1,ic
          ibas=io+i
          do l=1,lc
            lbas=lo+l
            do k=1,kc
              kbas=ko+k
              !> integrals processing: [ij|kl]
                 !> large-large block 
                 gg = dcmplx( gout(1,k,l,1,i,j), 0.0d0 )
                 !> large-small block 
                 ggaa = dcmplx( gout(2,k,l,1,i,j), gout(3,k,l,1,i,j) )
                 ggab = dcmplx( gout(4,k,l,1,i,j), gout(5,k,l,1,i,j) )
                 ggbb = dconjg( ggaa )
                 ggba =-dconjg( ggab )
                 !> small-large block 
                 gaag = dcmplx( gout(1,k,l,2,i,j), gout(1,k,l,3,i,j) )
                 gabg = dcmplx( gout(1,k,l,4,i,j), gout(1,k,l,5,i,j) )
                 gbbg = dconjg( gaag )
                 gbag =-dconjg( gabg ) 
                 !> small-small block
                 gaagaa = dcmplx( ( gout(2,k,l,2,i,j)-gout(3,k,l,3,i,j)), ( gout(3,k,l,2,i,j)+gout(2,k,l,3,i,j)) )
                 gaagba = dcmplx( (-gout(4,k,l,2,i,j)-gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)-gout(4,k,l,3,i,j)) )
                 gaagab = dcmplx( ( gout(4,k,l,2,i,j)-gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)+gout(4,k,l,3,i,j)) )
                 gaagbb = dcmplx( ( gout(2,k,l,2,i,j)+gout(3,k,l,3,i,j)), (-gout(3,k,l,2,i,j)+gout(2,k,l,3,i,j)) )
                 gabgaa = dcmplx( ( gout(2,k,l,4,i,j)-gout(3,k,l,5,i,j)), ( gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
                 gabgba = dcmplx( (-gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), ( gout(5,k,l,4,i,j)-gout(4,k,l,5,i,j)) )
                 gabgab = dcmplx( ( gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), ( gout(5,k,l,4,i,j)+gout(4,k,l,5,i,j)) )
                 gabgbb = dcmplx( ( gout(2,k,l,4,i,j)+gout(3,k,l,5,i,j)), (-gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
                !gbagbb =-dconjg( gabgaa )
                !gbagab = dconjg( gabgba )
                !gbagba = dconjg( gabgab )
                !gbagaa =-dconjg( gabgbb )
                !gbbgbb = dconjg( gaagaa )
                !gbbgab =-dconjg( gaagba ) 
                !gbbgba =-dconjg( gaagab )
                !gbbgaa = dconjg( gaagbb )
              !> end
              !> contraction [ij|kl]*lk = [ij|kl]*dconjg(kl)
                 paa = dconjg(pmat(2,2,kbas,lbas)                      )  !pmat(2,2,lbas,kbas)
                 pab = dconjg(pmat(4,2,kbas,lbas)                      )  !pmat(2,4,lbas,kbas)
                 pba = dconjg(pmat(2,4,kbas,lbas)                      )  !pmat(4,2,lbas,kbas)
                 pbb = dconjg(pmat(4,4,kbas,lbas)                      )  !pmat(4,4,lbas,kbas)
                 pll = dconjg(pmat(1,1,kbas,lbas) + pmat(3,3,kbas,lbas))  !pmat(1,1,lbas,kbas) + pmat(3,3,lbas,kbas)
                 !>
                 gg   =   gg*pll +   ggaa*paa +   ggab*pba +   ggba*pab +   ggbb*pbb
                 gaag = gaag*pll + gaagaa*paa + gaagab*pba + gaagba*pab + gaagbb*pbb
                 gabg = gabg*pll + gabgaa*paa + gabgab*pba + gabgba*pab + gabgbb*pbb
                !gbag = gbag*pll + gbagaa*paa + gbagab*pba + gbagba*pab + gbagbb*pbb
                !gbbg = gbbg*pll + gbbgaa*paa + gbbgab*pba + gbbgba*pab + gbbgbb*pbb
                 !>
                 gmat(1,1,ibas,jbas) = gmat(1,1,ibas,jbas) +        gg
                 gmat(3,3,ibas,jbas) = gmat(3,3,ibas,jbas) + dconjg(gg)
                 gmat(2,2,ibas,jbas) = gmat(2,2,ibas,jbas) +        gaag
                 gmat(4,2,ibas,jbas) = gmat(4,2,ibas,jbas) - dconjg(gabg) !+ gbag
                 gmat(2,4,ibas,jbas) = gmat(2,4,ibas,jbas) +        gabg
                 gmat(4,4,ibas,jbas) = gmat(4,4,ibas,jbas) + dconjg(gaag) !+ gbbg
              !> end
            enddo
          enddo
        enddo
      enddo
      !> exchange part  
      if( .not. doK )return
      !> pre-multiplication for scalling K-contribution 
      if( scalek /= 1.0d0 )then
        do j=1,jc
          do i=1,ic
            do m=1,5
              do l=1,lc
                do k=1,kc
                  gout(1,k,l,m,i,j) = gout(1,k,l,m,i,j)*scaleK
                  gout(2,k,l,m,i,j) = gout(2,k,l,m,i,j)*scaleK
                  gout(3,k,l,m,i,j) = gout(3,k,l,m,i,j)*scaleK
                  gout(4,k,l,m,i,j) = gout(4,k,l,m,i,j)*scaleK
                  gout(5,k,l,m,i,j) = gout(5,k,l,m,i,j)*scaleK
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
      !>
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          do j=1,jc
            jbas=jo+j
            !> large-large block
            paa_ll = pmat(1,1,jbas,kbas)
            pba_ll = pmat(3,1,jbas,kbas)
            pab_ll = pmat(1,3,jbas,kbas)
            pbb_ll = pmat(3,3,jbas,kbas)
            !> large-small block
            paa_ls = pmat(1,2,jbas,kbas)
            pba_ls = pmat(3,2,jbas,kbas)
            pab_ls = pmat(1,4,jbas,kbas)
            pbb_ls = pmat(3,4,jbas,kbas)
            !> small-large block
            paa_sl = pmat(2,1,jbas,kbas)
            pba_sl = pmat(4,1,jbas,kbas)
            pab_sl = pmat(2,3,jbas,kbas)
            pbb_sl = pmat(4,3,jbas,kbas)
            !> small-small block
            paa_ss = pmat(2,2,jbas,kbas)
            pba_ss = pmat(4,2,jbas,kbas)
            pab_ss = pmat(2,4,jbas,kbas)
            pbb_ss = pmat(4,4,jbas,kbas)
            do i=1,ic
              ibas=io+i
              !> large-large block 
              gg = dcmplx( gout(1,k,l,1,i,j), 0.0d0 )
              !> large-small block 
              ggaa = dcmplx( gout(2,k,l,1,i,j), gout(3,k,l,1,i,j) )
              ggab = dcmplx( gout(4,k,l,1,i,j), gout(5,k,l,1,i,j) )
              ggbb = dconjg( ggaa )
              ggba =-dconjg( ggab )
              !> small-large block 
              gaag = dcmplx( gout(1,k,l,2,i,j), gout(1,k,l,3,i,j) )
              gabg = dcmplx( gout(1,k,l,4,i,j), gout(1,k,l,5,i,j) )
              gbbg = dconjg( gaag )
              gbag =-dconjg( gabg ) 
              !> small-small block: aa
              gaagaa = dcmplx( ( gout(2,k,l,2,i,j)-gout(3,k,l,3,i,j)), ( gout(3,k,l,2,i,j)+gout(2,k,l,3,i,j)) )
              gaagba = dcmplx( (-gout(4,k,l,2,i,j)-gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)-gout(4,k,l,3,i,j)) )
              gabgaa = dcmplx( ( gout(2,k,l,4,i,j)-gout(3,k,l,5,i,j)), ( gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
              gabgba = dcmplx( (-gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), ( gout(5,k,l,4,i,j)-gout(4,k,l,5,i,j)) )
              !> small-small block: bb
              gbbgbb = dconjg( gaagaa )
              gbbgab =-dconjg( gaagba ) 
              gbagbb =-dconjg( gabgaa )
              gbagab = dconjg( gabgba )
             !gbbgbb = dcmplx( ( gout(2,k,l,2,i,j)-gout(3,k,l,3,i,j)), (-gout(3,k,l,2,i,j)-gout(2,k,l,3,i,j)) )
             !gbbgab = dcmplx( ( gout(4,k,l,2,i,j)+gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)-gout(4,k,l,3,i,j)) )
             !gbagbb = dcmplx( (-gout(2,k,l,4,i,j)+gout(3,k,l,5,i,j)), ( gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
             !gbagab = dcmplx( (-gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), (-gout(5,k,l,4,i,j)+gout(4,k,l,5,i,j)) )
              !> small-small block: ab
              gaagab = dcmplx( ( gout(4,k,l,2,i,j)-gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)+gout(4,k,l,3,i,j)) )
              gaagbb = dcmplx( ( gout(2,k,l,2,i,j)+gout(3,k,l,3,i,j)), (-gout(3,k,l,2,i,j)+gout(2,k,l,3,i,j)) )
              gabgab = dcmplx( ( gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), ( gout(5,k,l,4,i,j)+gout(4,k,l,5,i,j)) )
              gabgbb = dcmplx( ( gout(2,k,l,4,i,j)+gout(3,k,l,5,i,j)), (-gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
              !> small-small block: ba
              gbbgba =-dconjg( gaagab )
              gbbgaa = dconjg( gaagbb )
              gbagba = dconjg( gabgab )
              gbagaa =-dconjg( gabgbb )
             !gbbgba = dcmplx( (-gout(4,k,l,2,i,j)+gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)+gout(4,k,l,3,i,j)) )
             !gbbgaa = dcmplx( ( gout(2,k,l,2,i,j)+gout(3,k,l,3,i,j)), ( gout(3,k,l,2,i,j)-gout(2,k,l,3,i,j)) )
             !gbagba = dcmplx( ( gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), (-gout(5,k,l,4,i,j)-gout(4,k,l,5,i,j)) )
             !gbagaa = dcmplx( (-gout(2,k,l,4,i,j)-gout(3,k,l,5,i,j)), (-gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
              !> large-large block
              gmat(1,1,ibas,lbas) = gmat(1,1,ibas,lbas) - paa_ll*gg
              gmat(1,3,ibas,lbas) = gmat(1,3,ibas,lbas) - pab_ll*gg
              gmat(3,1,ibas,lbas) = gmat(3,1,ibas,lbas) - pba_ll*gg
              gmat(3,3,ibas,lbas) = gmat(3,3,ibas,lbas) - pbb_ll*gg
              !> large-small block
              gmat(1,2,ibas,lbas) = gmat(1,2,ibas,lbas) - paa_ls*ggaa - pab_ls*ggba  
              gmat(1,4,ibas,lbas) = gmat(1,4,ibas,lbas) - paa_ls*ggab - pab_ls*ggbb
              gmat(3,2,ibas,lbas) = gmat(3,2,ibas,lbas) - pba_ls*ggaa - pbb_ls*ggba
              gmat(3,4,ibas,lbas) = gmat(3,4,ibas,lbas) - pba_ls*ggab - pbb_ls*ggbb
              !> small-large block
              gmat(2,1,ibas,lbas) = gmat(2,1,ibas,lbas) - gaag*paa_sl - gabg*pba_sl 
              gmat(2,3,ibas,lbas) = gmat(2,3,ibas,lbas) - gaag*pab_sl - gabg*pbb_sl 
              gmat(4,1,ibas,lbas) = gmat(4,1,ibas,lbas) - gbag*paa_sl - gbbg*pba_sl 
              gmat(4,3,ibas,lbas) = gmat(4,3,ibas,lbas) - gbag*pab_sl - gbbg*pbb_sl 
              !> small-small block
              gmat(2,2,ibas,lbas) = gmat(2,2,ibas,lbas) - gaagaa*paa_ss - gaagba*pab_ss - gabgaa*pba_ss - gabgba*pbb_ss 
              gmat(2,4,ibas,lbas) = gmat(2,4,ibas,lbas) - gaagab*paa_ss - gaagbb*pab_ss - gabgab*pba_ss - gabgbb*pbb_ss 
              gmat(4,2,ibas,lbas) = gmat(4,2,ibas,lbas) - gbagaa*paa_ss - gbagba*pab_ss - gbbgaa*pba_ss - gbbgba*pbb_ss 
              gmat(4,4,ibas,lbas) = gmat(4,4,ibas,lbas) - gbagab*paa_ss - gbagbb*pab_ss - gbbgab*pba_ss - gbbgbb*pbb_ss 
            enddo
          enddo
        enddo
      enddo
    end subroutine

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE two_electron_integrals_4comp_8x( dG, dP, eri_class, scaleK ) 

    !> input
    real(8),          intent(in) :: scaleK
    character(len=4), intent(in) :: eri_class
    complex(8),       intent(in) :: dP(8*nshl*nrow*nrow)
    !> output
    complex(8), intent(inout) :: dG(8*nshl*nrow*nrow)
    !> local
    logical :: doK
    real(8) :: Fc,Fx
    real(8) :: c2,cc4
    real(8) :: fij, fijkl
    logical :: AB, AC, AD
    integer :: i, j, k, l
    integer :: li, ni, oi, ai, imax     
    integer :: lj, nj, oj, aj   
    integer :: lk, nk, ok, ak   
    integer :: ll, nl, ol, al, lupper  
    integer :: ij, kl, ik, il, jk, jl
    real(8) :: ei, ci, di, xi, yi, zi 
    real(8) :: ej, cj, dj, xj, yj, zj 
    real(8) :: ek, ck, dk, xk, yk, zk 
    real(8) :: el, cl, dl, xl, yl, zl 

    integer :: lmax
    integer :: nint
    real(8) :: gout(441*441*25)    !todo: rozumnejsia definicia
    !> coulomb
    real(8), allocatable :: dGij(:) 
    real(8), allocatable :: dPij(:) 
    !>exchange
    complex(8), allocatable :: dPi(:) 
    complex(8), allocatable :: dPj(:) 
    complex(8), allocatable :: dGi(:) 
    complex(8), allocatable :: dGj(:) 
    !> screening
    real(8)              :: kmax(2,2)
    real(8), allocatable :: smax(:,:,:)
    real(8), allocatable :: jmax(:,:,:)
    real(8), allocatable :: gmax(:,:,:)
    real(8), allocatable :: pmax(:,:,:,:)
    character(len=4)     :: eri_class_local
    !> MPI
    integer              :: mytask, cpu, ier
    integer, allocatable :: eri_task_list(:)    

    
    !> initialize
    c2   = 2.0d0*cvalue
    cc4  = 4.0d0*cvalue*cvalue
    lmax = maxval(gto(:)%lvalue) 
    nint = (lmax*(lmax+1))/2

    !> coulomb 
    allocate( dGij(5*nint*nint) )
    allocate( dPij(5*nint*nint) )

    !> exchange 
    doK = ( scaleK /= 0.0d0 )
    if( doK )then
      allocate( dPi(8*nshl*nint*nrow) )
      allocate( dPj(8*nshl*nint*nrow) )
      allocate( dGi(8*nshl*nint*nrow) )
      allocate( dGj(8*nshl*nint*nrow) )
    endif

    !> creating a simple MPI eri task list
    allocate( eri_task_list(nr_shells*(nr_shells+1)/2), stat=ier )
    if( ier.ne.0 )stop ' Error in allocation: eri_task_list(:)'
    ij =1
    cpu=1
    mytask=0
    do i=1,nr_shells
      do j=1,i
        if( cpu > nr_cpus )cpu=1
        eri_task_list(ij)=cpu-1
        ij =ij +1
        cpu=cpu+1
      enddo
    enddo

    allocate(gmax(  5,nr_shells,nr_shells))
    allocate(pmax(4,4,nr_shells,nr_shells))
    allocate(smax(  2,nr_shells,nr_shells))
    allocate(jmax(  2,nr_shells,nr_shells))
    call shell_guess_4comp( nr_shells, Jgmax=smax, Jrmax=jmax, Kgmax=gmax, pmax=pmax, cPmat=dP )


 iloop:do i=1,nr_shells
         ij   = (i-1)*nr_shells
         li   = gto(i)%lvalue
         oi   = gto(i)%offset
         ni   = gto(i)%cdegen
         ai   = gto(i)%origin
         ei   = gto(i)%exponent(1)
         ci   = gto(i)%coefficient(1)
         di   = 1.0d0 
         xi   = atom(ai)%coordinate_x
         yi   = atom(ai)%coordinate_y
         zi   = atom(ai)%coordinate_z
         imax = oi+ni

         if( doK )call get_dP_exchange( dP, nrow, oi, ni, dPi, dGi, imax )
          
 jloop:  do j=1,i
           ij=ij+1

           mytask=mytask+1 
           if( myrank /= eri_task_list(mytask) )cycle

           lj = gto(j)%lvalue
           oj = gto(j)%offset
           nj = gto(j)%cdegen
           aj = gto(j)%origin
           ej = gto(j)%exponent(1)
           cj = gto(j)%coefficient(1)
           dj = 1.0d0 
           xj = atom(aj)%coordinate_x
           yj = atom(aj)%coordinate_y
           zj = atom(aj)%coordinate_z
           AB = (ai==aj)

                     fij = 1.0d0 
           if( i==j )fij = 0.5d0
                  
                    call get_dP_coulomb ( dP, nrow, oi, ni, oj, nj, dPij, dGij      )
           if( doK )call get_dP_exchange( dP, nrow,         oj, nj, dPj,  dGj, imax )
           
 kloop:    do k=1,i
             kl = (k-1)*nr_shells
             lk = gto(k)%lvalue
             ok = gto(k)%offset
             nk = gto(k)%cdegen
             ak = gto(k)%origin
             ek = gto(k)%exponent(1)
             ck = gto(k)%coefficient(1)
             dk = 1.0d0 
             xk = atom(ak)%coordinate_x
             yk = atom(ak)%coordinate_y
             zk = atom(ak)%coordinate_z
             AC = (ai==ak)
       
                       lupper=k
             if( k==i )lupper=j
 lloop:      do l=1,lupper
               kl = kl+1
               ll = gto(l)%lvalue
               nl = gto(l)%cdegen
               ol = gto(l)%offset
               al = gto(l)%origin
               el = gto(l)%exponent(1)
               cl = gto(l)%coefficient(1)
               dl = 1.0d0 
               xl = atom(al)%coordinate_x
               yl = atom(al)%coordinate_y
               zl = atom(al)%coordinate_z
               AD = (ai==al)

               if( AB.and.AC.and.AD.and.mod(li+lj+lk+ll,2)/=0 )cycle

               fijkl=fij;if(  k==l  )fijkl=0.5d0*fijkl
                         if( ij==kl )fijkl=0.5d0*fijkl

               call estimate_dG_8x(eri_class_local,i,j,k,l,di,dj,dk,dl,kmax,pmax,gmax,nr_shells) 
               
               SELECT CASE( eri_class )
                 CASE( 'ssss' )
                     Fx = kmax(2,2)
                     Fc = dmax1(smax(2,j,i)*(jmax(2,l,k)-jmax(1,l,k)),&
                                smax(2,l,k)*(jmax(2,j,i)-jmax(1,j,i)) ) 
                     if( Fc<screening_threshold .and. Fx<screening_threshold )eri_class_local='ssll'
                     !>
                     Fx = dmax1(Fx,kmax(2,1),kmax(1,2))
                     Fc = dmax1(smax(2,j,i)*jmax(2,l,k),smax(2,l,k)*jmax(2,j,i)) 
                     if( Fc<screening_threshold .and. Fx<screening_threshold )eri_class_local='llll'
                     !>
                     Fx = dmax1(Fx,kmax(1,1) )
                     Fc = dmax1(Fc,smax(1,j,i)*jmax(2,l,k),smax(1,l,k)*jmax(2,j,i)) 
                     if( Fc<screening_threshold .and. Fx<screening_threshold )cycle  
                 CASE( 'ssll' )
                     Fx = dmax1(kmax(2,1),kmax(1,2))
                     Fc = dmax1(smax(2,j,i)*jmax(2,l,k),smax(2,l,k)*jmax(2,j,i)) 
                     if( Fc<screening_threshold .and. Fx<screening_threshold )eri_class_local='llll'
                     !>
                     Fx = dmax1(Fx,kmax(1,1) )
                     Fc = dmax1(Fc,smax(1,j,i)*jmax(2,l,k),smax(1,l,k)*jmax(2,j,i)) 
                     if( Fc<screening_threshold .and. Fx<screening_threshold )cycle  
                 CASE( 'llll' )
                     Fx = kmax(1,1)
                     Fc = dmax1( smax(1,j,i)*jmax(2,l,k),smax(1,l,k)*jmax(2,j,i) ) 
                     if( Fc<screening_threshold .and. Fx<screening_threshold )cycle  
               END SELECT
       
               SELECT CASE( eri_class_local )
                 CASE( 'llll' )
                     call interest_eri('llll',              &
                                       fijkl,gout,nint,     &
                                       li,ei,ai,xi,yi,zi,ci,&
                                       lj,ej,aj,xj,yj,zj,cj,&
                                       lk,ek,ak,xk,yk,zk,ck,&
                                       ll,el,al,xl,yl,zl,cl )
                     if( nshl==1 )call process_dG_llll_8x_trs(ni,nj,nk,nl,oi,oj,ok,ol,   &
                                       gout,dP,dG,nrow,doK,dPij,dGij,dPi,dPj,dGi,dGj,imax)
                    !if( nshl==2 )call process_dG_llll_8x_gen(ni,nj,nk,nl,oi,oj,ok,ol,   &
                    !                  gout,dP,dG,nrow,doK,dPij,dGij,dPi,dPj,dGi,dGj,imax)
                 CASE( 'ssll' )
                     call interest_eri('ssll',              &
                                       fijkl,gout,nint,     &
                                       li,ei,ai,xi,yi,zi,ci,&
                                       lj,ej,aj,xj,yj,zj,cj,&
                                       lk,ek,ak,xk,yk,zk,ck,&
                                       ll,el,al,xl,yl,zl,cl )
                     if( nshl==1 )call process_dG_ssll_8x_trs(ni,nj,nk,nl,oi,oj,ok,ol,   &
                                       gout,dP,dG,nrow,doK,dPij,dGij,dPi,dPj,dGi,dGj,imax)
                    !if( nshl==2 )call process_dG_ssll_8x_gen(ni,nj,nk,nl,oi,oj,ok,ol,   &
                    !                  gout,dP,dG,nrow,doK,dPij,dGij,dPi,dPj,dGi,dGj,imax)
                 CASE( 'ssss' )
                     call interest_eri('ssss',              &
                                       fijkl,gout,nint,     &
                                       li,ei,ai,xi,yi,zi,ci,&
                                       lj,ej,aj,xj,yj,zj,cj,&
                                       lk,ek,ak,xk,yk,zk,ck,&
                                       ll,el,al,xl,yl,zl,cl )
                     if( nshl==1 )call process_dG_ssss_8x_trs(ni,nj,nk,nl,oi,oj,ok,ol,   &
                                       gout,dP,dG,nrow,doK,dPij,dGij,dPi,dPj,dGi,dGj,imax)
                    !if( nshl==2 )call process_dG_ssss_8x_gen(ni,nj,nk,nl,oi,oj,ok,ol,   &
                    !                  gout,dP,dG,nrow,doK,dPij,dGij,dPi,dPj,dGi,dGj,imax)
                 CASE DEFAULT; stop 'Error: the ERI class is not implemented yet (4comp_8x) => stop'
               END SELECT

             enddo lloop
           enddo kloop
                     call add_dG_coulomb ( dG, nrow, oi, ni, oj, nj, dGij      ) 
           if( doK ) call add_dG_exchange( dG, nrow,         oj, nj, dGj, imax ) 
         enddo jloop
         if( doK ) call add_dG_exchange( dG, nrow, oi, ni, dGi, imax )
       enddo iloop
    !> end

    !> finalize
    deallocate( dGij )
    deallocate( dPij )
    if( doK )then
      deallocate( dPi )
      deallocate( dPj )
      deallocate( dGi )
      deallocate( dGj )
    endif
    deallocate( smax )
    deallocate( jmax )
    deallocate( gmax )
    deallocate( pmax )
    deallocate( eri_task_list )

    call hermitize_matrix( dG )
    
    CONTAINS
  
       subroutine get_dP_exchange( dP, nbas, io, ic, dP_block, dG_block, nbasr )
         !> input
         integer,    intent(in) :: io   
         integer,    intent(in) :: ic   
         integer,    intent(in) :: nbas 
         integer,    intent(in) :: nbasr 
         complex(8), intent(in) :: dP(2*nshl,4,nbas,nbas)
         !> output
         complex(8), intent(out) :: dP_block(2*nshl,4,ic,nbasr)
         complex(8), intent(out) :: dG_block(2*nshl,4,ic,nbasr)
         !> local
         integer :: ibas, i, j, k 
         !> 
         select case( nshl )
           case(1)
             do j=1,nbasr
               do i=1,ic
                 ibas = io+i
                 do k=1,4
                   dP_block(1,k,i,j) = dP(1,k,ibas,j)
                   dP_block(2,k,i,j) = dP(2,k,ibas,j)
                   dG_block(1,k,i,j) = dcmplx(0.0d0,0.0d0)
                   dG_block(2,k,i,j) = dcmplx(0.0d0,0.0d0)
                 enddo
               enddo
             enddo
           case(2)
             do j=1,nbasr
               do i=1,ic
                 ibas = io+i
                 do k=1,4
                   dP_block(1,k,i,j) = dP(1,k,ibas,j)
                   dP_block(2,k,i,j) = dP(2,k,ibas,j)
                   dP_block(3,k,i,j) = dP(3,k,ibas,j)
                   dP_block(4,k,i,j) = dP(4,k,ibas,j)
                   dG_block(1,k,i,j) = dcmplx(0.0d0,0.0d0)
                   dG_block(2,k,i,j) = dcmplx(0.0d0,0.0d0)
                   dG_block(3,k,i,j) = dcmplx(0.0d0,0.0d0)
                   dG_block(4,k,i,j) = dcmplx(0.0d0,0.0d0)
                 enddo
               enddo
             enddo
         end select
       end subroutine
  
       subroutine get_dP_coulomb( dP, nbas, io, ic, jo, jc, dP_block, dG_block )
         !> input
         integer,    intent(in) :: nbas 
         integer,    intent(in) :: io, jo   
         integer,    intent(in) :: ic, jc   
         complex(8), intent(in) :: dP(2*nshl,4,nbas,nbas)
         !> output
         real(8), intent(out) :: dP_block(5,ic,jc)
         real(8), intent(out) :: dG_block(5,ic,jc)
         !> local
         integer :: ibas, jbas, i, j, k 
         !> 
         select case( nshl )
           case(1)
             do j=1,jc
               jbas=jo+j
               do i=1,ic
                 ibas=io+i
                 dG_block(1,i,j) = 0.0d0
                 dG_block(2,i,j) = 0.0d0
                 dG_block(3,i,j) = 0.0d0
                 dG_block(4,i,j) = 0.0d0
                 dG_block(5,i,j) = 0.0d0
                 dP_block(1,i,j) = 4.0d0*(dreal(dP(1,1,ibas,jbas))) !Re(aa+bb)=2Re(aa)
                 dP_block(2,i,j) = 4.0d0*(dreal(dP(2,2,ibas,jbas))) !Re(aa+bb)=2Re(aa)
                 dP_block(3,i,j) = 4.0d0*(dimag(dP(2,2,ibas,jbas))) !Im(aa-bb)=2Im(aa)
                 dP_block(4,i,j) = 4.0d0*(dreal(dP(2,4,ibas,jbas))) !Re(ab-ba)=2Re(ab)
                 dP_block(5,i,j) = 4.0d0*(dimag(dP(2,4,ibas,jbas))) !Im(ab+ba)=2Im(ab)
               enddo
             enddo
           case(2)
             do j=1,jc
               jbas=jo+j
               do i=1,ic
                 ibas=io+i
                 dG_block(1,i,j) = 0.0d0
                 dG_block(2,i,j) = 0.0d0
                 dG_block(3,i,j) = 0.0d0
                 dG_block(4,i,j) = 0.0d0
                 dG_block(5,i,j) = 0.0d0
                 dP_block(1,i,j) = 2.0d0*(dreal(dP(1,1,ibas,jbas))+dreal(dP(3,3,ibas,jbas))) !Re(aa+bb)
                 dP_block(2,i,j) = 2.0d0*(dreal(dP(2,2,ibas,jbas))+dreal(dP(4,4,ibas,jbas))) !Re(aa+bb)
                 dP_block(3,i,j) = 2.0d0*(dimag(dP(2,2,ibas,jbas))-dimag(dP(4,4,ibas,jbas))) !Im(aa-bb)
                 dP_block(4,i,j) = 2.0d0*(dreal(dP(2,4,ibas,jbas))-dreal(dP(4,2,ibas,jbas))) !Re(ab-ba)
                 dP_block(5,i,j) = 2.0d0*(dimag(dP(2,4,ibas,jbas))+dimag(dP(4,2,ibas,jbas))) !Im(ab+ba)
               enddo
             enddo
           end select    
       end subroutine

       subroutine add_dG_exchange( dG, nbas, io, ic, dG_block, nbasr )
         !> note: dG is Hermitian on output
         !> input
         integer, intent(in) :: io   
         integer, intent(in) :: ic   
         integer, intent(in) :: nbas 
         integer, intent(in) :: nbasr 
         !> output
         complex(8), intent(inout) :: dG(2*nshl,4,nbas,nbas)
         complex(8), intent(inout) :: dG_block(2*nshl,4,ic,nbasr)
         !> local
         integer :: i, ibas, jbas, k 
         !> 
         select case( nshl )
           case(1)
             do jbas=1,nbasr
               do i=1,ic
                 dG_block(1,1,i,jbas) = dG_block(1,1,i,jbas)*dcmplx(scaleK,0.0d0)
                 dG_block(2,1,i,jbas) = dG_block(2,1,i,jbas)*dcmplx(scaleK,0.0d0)
                 dG_block(1,2,i,jbas) = dG_block(1,2,i,jbas)*dcmplx(scaleK,0.0d0)
                 dG_block(2,2,i,jbas) = dG_block(2,2,i,jbas)*dcmplx(scaleK,0.0d0)
                 dG_block(1,3,i,jbas) = dG_block(1,3,i,jbas)*dcmplx(scaleK,0.0d0)
                 dG_block(2,3,i,jbas) = dG_block(2,3,i,jbas)*dcmplx(scaleK,0.0d0)
                 dG_block(1,4,i,jbas) = dG_block(1,4,i,jbas)*dcmplx(scaleK,0.0d0)
                 dG_block(2,4,i,jbas) = dG_block(2,4,i,jbas)*dcmplx(scaleK,0.0d0)
               enddo
             enddo
             !>
             do jbas=1,nbasr
               do i=1,ic
                 ibas = io+i
                 dG(1,1,ibas,jbas) = dG(1,1,ibas,jbas) + dG_block(1,1,i,jbas)
                 dG(2,1,ibas,jbas) = dG(2,1,ibas,jbas) + dG_block(2,1,i,jbas)
                 dG(1,2,ibas,jbas) = dG(1,2,ibas,jbas) + dG_block(1,2,i,jbas)
                 dG(2,2,ibas,jbas) = dG(2,2,ibas,jbas) + dG_block(2,2,i,jbas)
                 dG(1,3,ibas,jbas) = dG(1,3,ibas,jbas) + dG_block(1,3,i,jbas)
                 dG(2,3,ibas,jbas) = dG(2,3,ibas,jbas) + dG_block(2,3,i,jbas)
                 dG(1,4,ibas,jbas) = dG(1,4,ibas,jbas) + dG_block(1,4,i,jbas)
                 dG(2,4,ibas,jbas) = dG(2,4,ibas,jbas) + dG_block(2,4,i,jbas)
               enddo
             enddo
             !> 
             do i=1,ic
               ibas = io+i
               do jbas=1,nbasr
                 dG(1,1,jbas,ibas) = dG(1,1,jbas,ibas) + dconjg(dG_block(1,1,i,jbas))
                 dG(2,1,jbas,ibas) = dG(2,1,jbas,ibas) + dconjg(dG_block(1,2,i,jbas))
                 dG(1,2,jbas,ibas) = dG(1,2,jbas,ibas) + dconjg(dG_block(2,1,i,jbas))
                 dG(2,2,jbas,ibas) = dG(2,2,jbas,ibas) + dconjg(dG_block(2,2,i,jbas))
                 dG(1,3,jbas,ibas) = dG(1,3,jbas,ibas) -        dG_block(1,3,i,jbas)  !AB(j,i)=-AB(i,j)
                 dG(2,3,jbas,ibas) = dG(2,3,jbas,ibas) -        dG_block(1,4,i,jbas)  !aB(j,i)=-Ab(i,j)
                 dG(1,4,jbas,ibas) = dG(1,4,jbas,ibas) -        dG_block(2,3,i,jbas)  !Ab(j,i)=-aB(i,j)
                 dG(2,4,jbas,ibas) = dG(2,4,jbas,ibas) -        dG_block(2,4,i,jbas)  !ab(j,i)=-ab(i,j)
               enddo
             enddo
           case(2)
             do jbas=1,nbasr
               do i=1,ic
                 do k=1,4
                   dG_block(1,k,i,jbas) = dG_block(1,k,i,jbas)*dcmplx(scaleK,0.0d0)
                   dG_block(2,k,i,jbas) = dG_block(2,k,i,jbas)*dcmplx(scaleK,0.0d0)
                   dG_block(3,k,i,jbas) = dG_block(3,k,i,jbas)*dcmplx(scaleK,0.0d0)
                   dG_block(4,k,i,jbas) = dG_block(4,k,i,jbas)*dcmplx(scaleK,0.0d0)
                 enddo
               enddo
             enddo
             !>
             do jbas=1,nbasr
               do i=1,ic
                 ibas = io+i
                 do k=1,4
                   dG(1,k,ibas,jbas) = dG(1,k,ibas,jbas) + dG_block(1,k,i,jbas)
                   dG(2,k,ibas,jbas) = dG(2,k,ibas,jbas) + dG_block(2,k,i,jbas)
                   dG(3,k,ibas,jbas) = dG(3,k,ibas,jbas) + dG_block(3,k,i,jbas)
                   dG(4,k,ibas,jbas) = dG(4,k,ibas,jbas) + dG_block(4,k,i,jbas)
                 enddo
               enddo
             enddo
             !> 
             do i=1,ic
               ibas = io+i
               do jbas=1,nbasr
                 do k=1,4
                   dG(1,k,jbas,ibas) = dG(1,k,jbas,ibas) + dconjg(dG_block(k,1,i,jbas))
                   dG(2,k,jbas,ibas) = dG(2,k,jbas,ibas) + dconjg(dG_block(k,2,i,jbas))
                   dG(3,k,jbas,ibas) = dG(3,k,jbas,ibas) + dconjg(dG_block(k,3,i,jbas))
                   dG(4,k,jbas,ibas) = dG(4,k,jbas,ibas) + dconjg(dG_block(k,4,i,jbas))
                 enddo
               enddo
             enddo
           end select
       end subroutine

       subroutine add_dG_coulomb( dG, nbas, io, ic, jo, jc, dG_block )
         !> note: J-matrix is Hermitian and time-reversal symmetric on output
         !> input
         integer, intent(in) :: nbas 
         integer, intent(in) :: io, jo   
         integer, intent(in) :: ic, jc   
         real(8), intent(in) :: dG_block(5,ic,jc)
         !> input/output
         complex(8), intent(inout) :: dG(2*nshl,4,nbas,nbas)
         !> local
         integer :: ibas, jbas, i, j
         !>  
         select case( nshl )
           case(1)
             do j=1,jc   
               jbas = jo+j
               do i=1,ic
                 ibas = io+i
                 dG(1,1,ibas,jbas) = dG(1,1,ibas,jbas) + dcmplx( dG_block(1,i,j), 0.0d0           )
                 dG(2,2,ibas,jbas) = dG(2,2,ibas,jbas) + dcmplx( dG_block(2,i,j), dG_block(3,i,j) )
                 dG(2,4,ibas,jbas) = dG(2,4,ibas,jbas) + dcmplx( dG_block(4,i,j), dG_block(5,i,j) )
               enddo
             enddo
             !>
             if( io/=jo )return
             !>
             do i=1,ic   
               ibas = io+i
               do j=1,jc
                 jbas = jo+j
                 dG(1,1,jbas,ibas) = dG(1,1,jbas,ibas) + dcmplx( dG_block(1,i,j), 0.0d0           ) 
                 dG(2,2,jbas,ibas) = dG(2,2,jbas,ibas) + dcmplx( dG_block(2,i,j),-dG_block(3,i,j) ) !aa(j,i)=+aa(i,j)*
                 dG(2,4,jbas,ibas) = dG(2,4,jbas,ibas) + dcmplx(-dG_block(4,i,j),-dG_block(5,i,j) ) !ab(j,i)=-ab(i,j)
               enddo
             enddo
           case(2)
             do j=1,jc   
               jbas = jo+j
               do i=1,ic
                 ibas = io+i
                 dG(1,1,ibas,jbas) = dG(1,1,ibas,jbas) + dcmplx( dG_block(1,i,j), 0.0d0           )
                 dG(3,3,ibas,jbas) = dG(3,3,ibas,jbas) + dcmplx( dG_block(1,i,j), 0.0d0           )
                 dG(2,2,ibas,jbas) = dG(2,2,ibas,jbas) + dcmplx( dG_block(2,i,j), dG_block(3,i,j) )
                 dG(4,2,ibas,jbas) = dG(4,2,ibas,jbas) + dcmplx(-dG_block(4,i,j), dG_block(5,i,j) )
                 dG(2,4,ibas,jbas) = dG(2,4,ibas,jbas) + dcmplx( dG_block(4,i,j), dG_block(5,i,j) )
                 dG(4,4,ibas,jbas) = dG(4,4,ibas,jbas) + dcmplx( dG_block(2,i,j),-dG_block(3,i,j) )
               enddo
             enddo
             !>
             if( io/=jo )return
             !>
             do i=1,ic   
               ibas = io+i
               do j=1,jc
                 jbas = jo+j
                 dG(1,1,jbas,ibas) = dG(1,1,jbas,ibas) + dcmplx( dG_block(1,i,j), 0.0d0           )
                 dG(3,3,jbas,ibas) = dG(3,3,jbas,ibas) + dcmplx( dG_block(1,i,j), 0.0d0           )
                 dG(2,2,jbas,ibas) = dG(2,2,jbas,ibas) + dcmplx( dG_block(2,i,j),-dG_block(3,i,j) ) !aa(j,i)=aa(i,j)*
                 dG(4,2,jbas,ibas) = dG(4,2,jbas,ibas) + dcmplx( dG_block(4,i,j),-dG_block(5,i,j) ) !ba(j,i)=ab(i,j)*
                 dG(2,4,jbas,ibas) = dG(2,4,jbas,ibas) + dcmplx(-dG_block(4,i,j),-dG_block(5,i,j) ) !ab(j,i)=ba(i,j)*
                 dG(4,4,jbas,ibas) = dG(4,4,jbas,ibas) + dcmplx( dG_block(2,i,j), dG_block(3,i,j) ) !bb(j,i)=bb(i,j)*
               enddo
             enddo
         end select
       end subroutine

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE shell_guess_4comp( nr_shells, Jgmax, Jrmax, Kgmax, pmax, rPmat, cPmat )
    !> on output, SS of Jgmax contains renormalized integral estimates
    !> on output, Jrmax contains contracted (P_KL*(I_KL)^{1/2}) estimates
    !> input
    integer,              intent(in) :: nr_shells
    real(8),    optional, intent(in) :: rPmat(       5,nrow,nrow)
    complex(8), optional, intent(in) :: cPmat(2*nshl,4,nrow,nrow)
    !> output
    real(8),           intent(inout) :: Jgmax(2,nr_shells,nr_shells)
    real(8),           intent(inout) :: Jrmax(2,nr_shells,nr_shells)
    real(8), optional, intent(inout) :: Kgmax(5,nr_shells,nr_shells)
    real(8), optional, intent(inout) :: Pmax(4,4,nr_shells,nr_shells)
    !> local 
    real(8) :: p(4,4)
    real(8) :: gout(441*441*25)
    integer :: k, ij, ijkl, nint
    real(8) :: ei, ci, di, xi, yi, zi 
    real(8) :: ej, cj, dj, xj, yj, zj 
    integer :: i, ii, li, ni, oi, ai, ibas
    integer :: j, jj, lj, nj, oj, aj, jbas
    real(8) :: rmax1, rmax2, rmax3, rmax4, rmax5
    real(8) :: gmax1, gmax2, gmax3, gmax4, gmax5
    !>
    ij=0
    do i=1,nr_shells
      li = gto(i)%lvalue
      oi = gto(i)%offset
      ni = gto(i)%cdegen
      ai = gto(i)%origin
      ei = gto(i)%exponent(1)
      ci = gto(i)%coefficient(1)
      di = 1.0d0 
      xi = atom(ai)%coordinate_x
      yi = atom(ai)%coordinate_y
      zi = atom(ai)%coordinate_z
      do j=1,i
        ij = ij+1
        lj = gto(j)%lvalue
        oj = gto(j)%offset
        nj = gto(j)%cdegen
        aj = gto(j)%origin
        ej = gto(j)%exponent(1)
        cj = gto(j)%coefficient(1)
        dj = 1.0d0
        xj = atom(aj)%coordinate_x
        yj = atom(aj)%coordinate_y
        zj = atom(aj)%coordinate_z
        !>
        call interest_eri('ssss',              &
                          1.0d0,gout,nint,     &
                          li,ei,ai,xi,yi,zi,ci,&
                          lj,ej,aj,xj,yj,zj,cj,&
                          li,ei,ai,xi,yi,zi,ci,&
                          lj,ej,aj,xj,yj,zj,cj )
        !>
        gmax1=0.0d0
        gmax2=0.0d0
        gmax3=0.0d0
        gmax4=0.0d0
        gmax5=0.0d0
        ijkl=0
        do jj=1,nj
          do ii=1,ni
            gmax1 = dmax1( gmax1, dsqrt(dabs(gout(index(ijkl+1,ijkl+1,5*ni*nj)))))
            gmax2 = dmax1( gmax2, dsqrt(dabs(gout(index(ijkl+2,ijkl+2,5*ni*nj)))))
            gmax3 = dmax1( gmax3, dsqrt(dabs(gout(index(ijkl+3,ijkl+3,5*ni*nj)))))
            gmax4 = dmax1( gmax4, dsqrt(dabs(gout(index(ijkl+4,ijkl+4,5*ni*nj)))))
            gmax5 = dmax1( gmax5, dsqrt(dabs(gout(index(ijkl+5,ijkl+5,5*ni*nj)))))
            ijkl=ijkl+5
          enddo
        enddo
        !>
        !> J-term
        !>
        rmax1=0.0d0
        rmax2=0.0d0
        rmax3=0.0d0
        rmax4=0.0d0
        rmax5=0.0d0
        if( present(rPmat) )then
          do jj=1,nj
            jbas=oj+jj
            do ii=1,ni
              ibas=oi+ii
              rmax1 = dmax1(rmax1,dabs(rPmat(1,ibas,jbas))) !Re(aa+bb)
              rmax2 = dmax1(rmax2,dabs(rPmat(2,ibas,jbas))) !Re(aa+bb)
              rmax3 = dmax1(rmax3,dabs(rPmat(3,ibas,jbas))) !Im(aa-bb)
              rmax4 = dmax1(rmax4,dabs(rPmat(4,ibas,jbas))) !Re(ab-ba)
              rmax5 = dmax1(rmax5,dabs(rPmat(5,ibas,jbas))) !Im(ab+ba)
            enddo
          enddo
        endif
        if( present(cPmat) )then
          select case( nshl )
            case( 1 )
              do jj=1,nj
                jbas=oj+jj
                do ii=1,ni
                  ibas=oi+ii
                  rmax1 = dmax1(rmax1,2.0d0*dabs(dreal(cPmat(1,1,ibas,jbas)))) !Re(aa+bb)
                  rmax2 = dmax1(rmax2,2.0d0*dabs(dreal(cPmat(2,2,ibas,jbas)))) !Re(aa+bb)
                  rmax3 = dmax1(rmax3,2.0d0*dabs(dimag(cPmat(2,2,ibas,jbas)))) !Im(aa-bb)
                  rmax4 = dmax1(rmax4,2.0d0*dabs(dreal(cPmat(2,4,ibas,jbas)))) !Re(ab-ba)
                  rmax5 = dmax1(rmax5,2.0d0*dabs(dimag(cPmat(2,4,ibas,jbas)))) !Im(ab+ba)
                enddo
              enddo
            case( 2 )
              do jj=1,nj
                jbas=oj+jj
                do ii=1,ni
                  ibas=oi+ii
                  rmax1 = dmax1(rmax1,dabs(dreal(cPmat(1,1,ibas,jbas))+dreal(cPmat(3,3,ibas,jbas)))) !Re(aa+bb)
                  rmax2 = dmax1(rmax2,dabs(dreal(cPmat(2,2,ibas,jbas))+dreal(cPmat(4,4,ibas,jbas)))) !Re(aa+bb)
                  rmax3 = dmax1(rmax3,dabs(dimag(cPmat(2,2,ibas,jbas))-dimag(cPmat(4,4,ibas,jbas)))) !Im(aa-bb)
                  rmax4 = dmax1(rmax4,dabs(dreal(cPmat(2,4,ibas,jbas))-dreal(cPmat(4,2,ibas,jbas)))) !Re(ab-ba)
                  rmax5 = dmax1(rmax5,dabs(dimag(cPmat(2,4,ibas,jbas))+dimag(cPmat(4,2,ibas,jbas)))) !Im(ab+ba)
                enddo
              enddo
          end select
        endif
        !> contract density matrix with integral estimates
        Jrmax(1,i,j) = rmax1*gmax1    !LL
        Jrmax(2,i,j) = rmax1*gmax1 &  !LL+SS
                     + rmax2*gmax2 &
                     + rmax3*gmax3 &
                     + rmax4*gmax4 &
                     + rmax5*gmax5
        Jrmax(1,j,i) = Jrmax(1,i,j)
        Jrmax(2,j,i) = Jrmax(2,i,j)
        !> estimate and renormalize SS block of integrals 
        !> (as they would contribute to the Fock matrix)
        Jgmax(1,i,j) = gmax1 
        Jgmax(2,i,j) = (dmax1(gmax2,gmax3,gmax4,gmax5))*di*dj
        Jgmax(1,j,i) = Jgmax(1,i,j)
        Jgmax(2,j,i) = Jgmax(2,i,j)
        !>
        !>
        !> K-term
        !>
        !>
        p(1:4,1)=0.0d0
        p(1:4,2)=0.0d0
        p(1:4,3)=0.0d0
        p(1:4,4)=0.0d0
        if( present(cPmat) )then
          select case( nshl )
            case( 1 )
              do jj=1,nj
                jbas=oj+jj
                do ii=1,ni
                  ibas=oi+ii
                  do k=1,4
                    !|z| = [Re(z)^2+Im(z)]^0.5 <= |Re(z)|+|Im(z)|
                    p(1,k) = dmax1(p(1,k),cdabs(cPmat(1,k,ibas,jbas)))
                    p(2,k) = dmax1(p(2,k),cdabs(cPmat(2,k,ibas,jbas)))
                  enddo
                enddo
              enddo
              p(3,3) = p(1,1)
              p(4,3) = p(2,1)
              p(3,4) = p(1,2)
              p(4,4) = p(2,2)
              p(3,1) = p(1,3)
              p(4,1) = p(2,3)
              p(3,2) = p(1,4)
              p(4,2) = p(2,4)
            case( 2 )
              do jj=1,nj
                jbas=oj+jj
                do ii=1,ni
                  ibas=oi+ii
                  do k=1,4
                    !|z| = [Re(z)^2+Im(z)]^0.5 <= |Re(z)|+|Im(z)|
                    p(1,k) = dmax1(p(1,k),cdabs(cPmat(1,k,ibas,jbas)))
                    p(2,k) = dmax1(p(2,k),cdabs(cPmat(2,k,ibas,jbas)))
                    p(3,k) = dmax1(p(3,k),cdabs(cPmat(3,k,ibas,jbas)))
                    p(4,k) = dmax1(p(4,k),cdabs(cPmat(4,k,ibas,jbas)))
                  enddo
                enddo
              enddo
          end select
        endif
        if( present(Pmax) )then
          do k=1,4
            Pmax(1,k,i,j) = p(1,k)
            Pmax(2,k,i,j) = p(2,k)
            Pmax(3,k,i,j) = p(3,k)
            Pmax(4,k,i,j) = p(4,k)
          enddo
          do k=1,4
            Pmax(1,k,j,i) = p(1,k)
            Pmax(2,k,j,i) = p(2,k)
            Pmax(3,k,j,i) = p(3,k)
            Pmax(4,k,j,i) = p(4,k)
          enddo
        endif
        if( present(Kgmax) )then
          Kgmax(1,i,j) = gmax1
          Kgmax(2,i,j) = gmax2
          Kgmax(3,i,j) = gmax3
          Kgmax(4,i,j) = gmax4
          Kgmax(5,i,j) = gmax5
          Kgmax(1,j,i) = gmax1
          Kgmax(2,j,i) = gmax2
          Kgmax(3,j,i) = gmax3
          Kgmax(4,j,i) = gmax4
          Kgmax(5,j,i) = gmax5
        endif
      enddo
    enddo

  END SUBROUTINE


! ------------------------------------------------------------------------------------
!>on input: lower triangle, on output: square symmetric matrix 
  SUBROUTINE symmetrize_matrix( factor, A, ida )

    !> input
    integer,           intent(in) :: ida
    real(8), optional, intent(in) :: factor

    !> output
    real(8), intent(inout) :: A(ida,ida)

    !> local
    integer :: i
    integer :: j


    if( present(factor) )then
      do j=1,ida
        do i=j,ida
          A(i,j) = A(i,j)*factor
          A(j,i) = A(i,j) 
        enddo
      enddo
    else
      do j=1,ida
        do i=j,ida
          A(j,i) = A(i,j) 
        enddo
      enddo
    endif

   !call backward_ao_resorting( A )
     
  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>on input: lower triangle, on output: square antisymmetric matrix 
  SUBROUTINE antisymmetrize_matrix( factor, A, ida )

    !> input
    integer,           intent(in) :: ida
    real(8), optional, intent(in) :: factor

    !> output
    real(8), intent(inout) :: A(ida,ida)

    !> local
    integer :: i
    integer :: j


    if( present(factor) )then
      do j=1,ida
        do i=j,ida
          A(i,j) =  A(i,j)*factor
          A(j,i) = -A(i,j) 
        enddo
      enddo
    else
      do j=1,ida
        do i=j,ida
          A(j,i) = -A(i,j) 
        enddo
      enddo
    endif

   !call backward_ao_resorting( A )
     
  END SUBROUTINE


! ------------------------------------------------------------------------------------
!>on input: the lower triangle of four-component Fock matrix, on output: full Hermite form 
!>
  SUBROUTINE hermitize_matrix( A )
    !> input/output
    complex(8), intent(inout) :: A(2*nshl,4,nrow,nrow)
    !> local
    integer :: i, j, m, n

    select case( nshl )
      case(1)
        do j=1,nrow
          do i=j,nrow
            A(1,1,j,i) = dconjg( A(1,1,i,j) )            
            A(1,2,j,i) = dconjg( A(2,1,i,j) )            
            A(2,1,j,i) = dconjg( A(1,2,i,j) )            
            A(2,2,j,i) = dconjg( A(2,2,i,j) )            
            A(1,3,j,i) =        -A(1,3,i,j)
            A(1,4,j,i) =        -A(2,3,i,j)
            A(2,3,j,i) =        -A(1,4,i,j)
            A(2,4,j,i) =        -A(2,4,i,j)
          enddo
        enddo
      case(2)
        do j=1,nrow
          do i=j,nrow
            do n=1,4
              do m=1,4
                A(n,m,j,i) = dconjg( A(m,n,i,j) ) 
              enddo
            enddo
          enddo
        enddo
      end select
     
  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>on input:  the matrix in "input"  basis order 
!>on output: the matrix in integral basis order 
  SUBROUTINE forward_ao_resorting( A )

    !> input/ouput
    real(8), intent(inout) :: A(nrow,ncol) 

    !> local
    real(8), allocatable :: B(:,:)
    integer :: i, ii, ic, io_new, io_old
    integer :: j, jj, jc, jo_new, jo_old


    allocate( B(nrow,ncol) ); B = A
    do j=1,nr_shells
      jc     = gto(j)%cdegen
      jo_new = gto(j)%offset
      jo_old = gto(j)%index 
      do i=1,nr_shells
        ic     = gto(i)%cdegen
        io_new = gto(i)%offset
        io_old = gto(i)%index 
        do jj=1,jc
          do ii=1,ic
            A(io_new+ii,jo_new+jj) = B(io_old+ii,jo_old+jj)
          enddo
        enddo
      enddo
    enddo
    deallocate( B )
        
  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>on input:  the matrix in "integral" basis order 
!>on output: the matrix in "input"    basis order 
  SUBROUTINE backward_ao_resorting( A )

    !> input/ouput
    real(8), intent(inout) :: A(nrow,ncol) 

    !> local
    real(8), allocatable :: B(:,:)
    integer :: i, ii, ic, io_new, io_old
    integer :: j, jj, jc, jo_new, jo_old


    allocate( B(nrow,ncol) ); B = A
    do j=1,nr_shells
      jc     = gto(j)%cdegen
      jo_new = gto(j)%index 
      jo_old = gto(j)%offset
      do i=1,nr_shells
        ic     = gto(i)%cdegen
        io_new = gto(i)%index 
        io_old = gto(i)%offset
        do jj=1,jc
          do ii=1,ic
            A(io_new+ii,jo_new+jj) = B(io_old+ii,jo_old+jj)
          enddo
        enddo
      enddo
    enddo
    deallocate( B )
        
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  FUNCTION index( i, j, n )

    !> input
    integer :: i
    integer :: j
    integer :: n

    !> output 
    integer :: index

    index = (j-1)*n + i
     
  END FUNCTION
! ------------------------------------------------------------------------------------
  SUBROUTINE transfer_gto_info(cart2lval)

    !> input/ouput
    integer, intent(out) :: cart2lval(*) 

    !> local
    integer :: i, j

    do i=1,size(gto)
      j = gto(i)%cdegen
      cart2lval((1+gto(i)%offset):(j+gto(i)%offset)) = gto(i)%lvalue
    enddo

  END SUBROUTINE
! ------------------------------------------------------------------------------------
END MODULE
