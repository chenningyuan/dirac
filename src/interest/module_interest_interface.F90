! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Module:       module_interest_interface.f90 
!
! Description:  - interfacing basic data between Dirac & InteRest program 
!               - providing two- and four-center cartesian integrals over RKB/RMB balanced basis
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

#ifdef VAR_MPI
  use mpi
#endif

  implicit none

  public initialize_interest 
  public initialize_interest_herd 
  public process_overlap_integrals_2comp 
  public process_overlap_integrals_4comp 
  public process_two_electron_integrals_4comp
  public process_two_electron_integrals_2comp
  public process_kinetic_energy_integrals_2comp 
  public process_kinetic_energy_integrals_4comp 
  public process_nuclear_attraction_integrals_2comp
  public process_nuclear_attraction_integrals_4comp

  public process_two_electron_integrals_4comp_slaves

  private

  !> matrix dimensions
     !> note: by definition nrow, ncol refer to the number of cartesian GTOs
     integer, save :: nshl = 1
     integer, save :: nrow = 0
     integer, save :: ncol = 0
     integer, save :: nr_shells = 0 
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
          integer :: worker
          integer :: lvalue
          integer :: origin
          integer :: sdegen 
          integer :: cdegen 
          real(8) :: exponent(1)
          real(8) :: coefficient(1)
     end  type
     type(type_gto), allocatable, target :: gto(:)
  !> end

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
       integer :: nr_cpus 
       type( type_gto ) :: tmpgto
    !> end
  

    !> test if the large component basis is uncontracted
       if( nplrg.ne.nlarge )then
         write(6,*) ' Error: contracted basis set found => stop!'
         write(6,*) ' Note:  uncontracted basis sets are only supported'
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
                              0,                       &
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
#ifdef VAR_MPI
       call mpi_comm_size(mpi_comm_world,nr_cpus,ier)
       do i=1,size(gto),nr_cpus    
         do j=0,nr_cpus-1
           if( i+j>size(gto) )exit
           gto(i+j)%worker = j
         enddo
       enddo
#endif
    !> end

    !> temporary print ... 
       do i=1,size(gto)
         write(6,'(1x,6(a,i3,3x),2(a,d15.5,3x))') ' i = ',i,                             &
                                                  ' index  = ',     gto(i)%index,        &
                                                  ' offset = ',     gto(i)%offset,       &
                                                  ' worker = ',     gto(i)%worker,       &
                                                  ' lvalue = ',     gto(i)%lvalue,       &
                                                  ' origin = ',     gto(i)%origin,       &
                                                  ' exponent = ',   gto(i)%exponent(1),  &
                                                  ' coefficient = ',gto(i)%coefficient(1)
       enddo

       write(6,*) 
       write(6,'(2x,a,i5  )') 'Total number of basis function shells:    ',size(gto)
       write(6,'(2x,a,i5  )') 'Total number of spherical basis functions:',sum(gto(:)%sdegen)
       write(6,'(2x,a,i5,a)') 'Total number of cartesian basis functions:',sum(gto(:)%cdegen),' (used in InteRest calculations)'
       write(6,*) 
    !> end

    !> global module variables
    nshl                = 1
    nrow                = sum( gto(:)%cdegen)
    ncol                = sum( gto(:)%cdegen)
    nr_shells           = size(gto)
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
    integer :: i,ier, nr_atoms, myrank
    integer(kind=mpi_address_kind) :: disp(9)                      !mpi_address_kind is always = 8
    integer                        :: blocklen(9), typex(9)
    integer                        :: newtype
#endif
  

    !> initialize InteRest integral package 
    call interest_initialize()

#ifdef VAR_MPI
    !> broadcast constants 
       !> define the local absolute address of each element within a derived data type
       call mpi_get_address(constant%alpha, disp(1),ier)
       call mpi_get_address(constant%alpha2,disp(2),ier)
       !> define the basic data type of each element within a derived data type
       typex(1:2) = mpi_real8
       !> define the length of each element within a derived data type
       blocklen(1:2) = 1
       !> redefine the absolute addresses to a relative address displacement (wrt to the address of the first element)
       disp(2) = disp(2) - disp(1) 
       disp(1) = disp(1) - disp(1) 
       !> create a MPI derived data type structure
       call mpi_type_create_struct(2,blocklen,disp,typex,newtype,ier)
       call mpi_type_commit(newtype,ier)
       call mpi_bcast(constant,1,newtype,0,mpi_comm_world,ier)
       call mpi_type_free(newtype,ier)
    !> end

    !> broadcast atom 
       if( allocated(atom) ) nr_atoms = size(atom)
       call mpi_bcast(nr_atoms,1,mpi_integer,0,mpi_comm_world,ier)
       if( .not. allocated(atom) ) allocate( atom(nr_atoms), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: atom(:)'
       !> define the local absolute address of each element within a derived data type
       call mpi_get_address(atom(1)%charge,      disp(1),ier)
       call mpi_get_address(atom(1)%coordinate_x,disp(2),ier)
       call mpi_get_address(atom(1)%coordinate_y,disp(3),ier)
       call mpi_get_address(atom(1)%coordinate_z,disp(4),ier)
       call mpi_get_address(atom(1)%gnu_exponent,disp(5),ier)
       !> define the basic data type of each element within a derived data type
       typex(1:5) = mpi_real8
       !> define the length of each element within a derived data type
       blocklen(1:5) = 1
       !> redefine the absolute addresses to a relative address displacement (wrt to the address of the first element)
       disp(5) = disp(5) - disp(1) 
       disp(4) = disp(4) - disp(1) 
       disp(3) = disp(3) - disp(1) 
       disp(2) = disp(2) - disp(1) 
       disp(1) = disp(1) - disp(1) 
       !> create a MPI derived data type structure
       call mpi_type_create_struct(5,blocklen,disp,typex,newtype,ier)
       call mpi_type_commit(newtype,ier)
       call mpi_bcast(atom,nr_atoms,newtype,0,mpi_comm_world,ier)
       call mpi_type_free(newtype,ier)
    !> end

    !> broadcast gto 
       if( allocated(gto) ) nr_shells = size(gto)
       call mpi_bcast(nr_shells,1,mpi_integer,0,mpi_comm_world,ier)
       if( .not. allocated(gto) ) allocate( gto(nr_shells), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: gto(:)'
       !> define the local absolute address of each element within a derived data type
       call mpi_get_address(gto(1)%index,         disp(1),ier)
       call mpi_get_address(gto(1)%offset,        disp(2),ier)
       call mpi_get_address(gto(1)%worker,        disp(3),ier)
       call mpi_get_address(gto(1)%lvalue,        disp(4),ier)
       call mpi_get_address(gto(1)%origin,        disp(5),ier)
       call mpi_get_address(gto(1)%sdegen,        disp(6),ier)
       call mpi_get_address(gto(1)%cdegen,        disp(7),ier)
       call mpi_get_address(gto(1)%exponent(1),   disp(8),ier)
       call mpi_get_address(gto(1)%coefficient(1),disp(9),ier)
       !> define the basic data type of each element within a derived data type
       typex(8:9) = mpi_real8
       typex(1:7) = mpi_integer
       !> define the length of each element within a derived data type
       blocklen(1:9) = 1
       !> redefine the absolute addresses to a relative address displacement (wrt to the address of the first element)
       disp(9) = disp(9) - disp(1) 
       disp(8) = disp(8) - disp(1) 
       disp(7) = disp(7) - disp(1) 
       disp(6) = disp(6) - disp(1) 
       disp(5) = disp(5) - disp(1) 
       disp(4) = disp(4) - disp(1) 
       disp(3) = disp(3) - disp(1) 
       disp(2) = disp(2) - disp(1) 
       disp(1) = disp(1) - disp(1) 
       !> create a MPI derived data type structure
       call mpi_type_create_struct(9,blocklen,disp,typex,newtype,ier)
       call mpi_type_commit(newtype,ier)
       call mpi_bcast(gto,nr_shells,newtype,0,mpi_comm_world,ier)
       call mpi_type_free(newtype,ier)
    !> end

    !> final broadcasting (global module variables)
       call mpi_bcast(nshl,               1,mpi_integer,0,mpi_comm_world,ier)
       call mpi_bcast(nrow,               1,mpi_integer,0,mpi_comm_world,ier)
       call mpi_bcast(ncol,               1,mpi_integer,0,mpi_comm_world,ier)
       call mpi_bcast(screening_threshold,1,mpi_real8,  0,mpi_comm_world,ier)
    !> end

    !> test printing
       call mpi_comm_rank( mpi_comm_world, myrank, ier )
       do i=1,size(gto)
         if( myrank /= gto(i)%worker )cycle
         write(6,*) ' Rank = ',myrank,gto(i)%exponent(1),gto(i)%lvalue
       enddo
      !call mpi_abort( mpi_comm_world, 1, ier )
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

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> Copy batch A -> matrix B 
  SUBROUTINE store_batch( A, ida, jda, nrkb, B, idb, jdb, ioffset, joffset )

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


    do k=1,nrkb
      do j=1,jda
        do i=1,ida
          B(ioffset+i,joffset+j,k) = A(k,i,j) 
        enddo
      enddo
    enddo
     
  END SUBROUTINE



! ------------------------------------------------------------------------------------
!> note: on input,  G and P must be in the cartesian GTOs 
  SUBROUTINE process_two_electron_integrals_2comp( Gmat, Pmat, intflg, hfxfac ) 

    use module_interest_eri

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

    complex(8), allocatable :: dP(:,:,:,:)
    complex(8), allocatable :: dG(:,:,:,:)


    !> check the eri class
    select case( intflg )
      case( 0 );    return 
      case( 1 );
      case default; stop "The desired class of ERI integrals is not implemented in InteRest :-( => stop"
    end select

    !> init for DFT
    doK = ( hfxfac /= 0.0d0 )

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
    !> end

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
    integer :: i, j, ier
    real(8) :: fll, fsl, fss 

    complex(8), allocatable :: dP(:,:,:,:)
    complex(8), allocatable :: dG(:,:,:,:)

    character(len=4) :: eri_class = 'ssss'


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
    call process_two_electron_integrals_4comp_4x( dG, dP, eri_class, hfxfac )

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
    integer :: ier
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
    call process_two_electron_integrals_4comp_4x( dG, dP, eri_class, hfxfac )

    !> send the two-electron Fock matrix of a particular slave upwards to the master
    call mpi_reduce( dG, mpi_in_place, 16*nrow*nrow, mpi_double_complex, mpi_sum, 0, mpi_comm_world, ier )

    deallocate( dP )
    deallocate( dG )
#endif

  END SUBROUTINE



! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_two_electron_integrals_4comp_1x( dG, dP, eri_class, hfxfac ) 

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


    !> init for DFT
    doK = ( hfxfac /= 0.0d0 )

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
  SUBROUTINE process_two_electron_integrals_4comp_4x( dG, dP, eri_class, hfxfac ) 

    use module_interest_eri

    !> input
    real(8),          intent(in) :: hfxfac
    character(len=4), intent(in) :: eri_class
    complex(8),       intent(in) :: dP(4,4,nrow,nrow)  

    !> output
    complex(8), intent(out) :: dG(4,4,nrow,nrow) 

    !> local
    logical    :: doK
    real(8)    :: a2, aa4
    real(8)    :: fij, fijkl
    integer    :: li, ni, oi   
    integer    :: lj, nj, oj   
    integer    :: lk, nk, ok   
    integer    :: ll, nl, ol   
    integer    :: i, j, k, l
    real(8)    :: ei, ci, xi, yi, zi 
    real(8)    :: ej, cj, xj, yj, zj 
    real(8)    :: ek, ck, xk, yk, zk 
    real(8)    :: el, cl, xl, yl, zl 

    integer    :: nint
    real(8)    :: gout(441*441*25)    !todo: rozumnejsia definicia
 
    complex(8), allocatable :: dGij(:,:,:) 
    complex(8), allocatable :: dPi (:,:,:) 
    complex(8), allocatable :: dPj (:,:,:) 
    complex(8), allocatable :: dGi (:,:,:) 
    complex(8), allocatable :: dGj (:,:,:) 

    integer              :: ierr 
    integer              :: myrank=0
    real(8)              :: Fc(3),Fx(3)
    real(8), allocatable :: gmax(:,:,:)
    real(8), allocatable :: pmax(:,:,:,:)
    character(len=4)     :: eri_class_local

    
    a2  = constant%alpha /2.0d0
    aa4 = constant%alpha2/4.0d0 
    
    !> init for DFT
    doK = ( hfxfac /= 0.0d0 )

    allocate(   dGij(4,4,21*21)  )   !todo: rozumnejsia definicia 
    if( doK )then
      allocate( dPi(4,4,21*nrow) )   !todo: rozumnejsia definicia 
      allocate( dPj(4,4,21*nrow) )   !todo: rozumnejsia definicia 
      allocate( dGi(4,4,21*nrow) )   !todo: rozumnejsia definicia 
      allocate( dGj(4,4,21*nrow) )   !todo: rozumnejsia definicia 
    endif

    allocate( gmax(  2,nr_shells,nr_shells) )
    allocate( pmax(2,2,nr_shells,nr_shells) )
    call max_shell_pair_contributions( dP, nrow, gmax, pmax, nr_shells )

#ifdef VAR_MPI
    call mpi_comm_rank(mpi_comm_world,myrank,ierr)
#endif

 iloop:do i=1,nr_shells
         if( myrank /= gto(i)%worker )cycle
         li =       gto(i)%lvalue
         oi =       gto(i)%offset
         ni =       gto(i)%cdegen
         ei =       gto(i)%exponent(1)
         xi = atom( gto(i)%origin )%coordinate_x
         yi = atom( gto(i)%origin )%coordinate_y
         zi = atom( gto(i)%origin )%coordinate_z
         ci =       gto(i)%coefficient(1)

         if( doK )call get_dP_exchange( dP, nrow, oi, ni, dPi, dGi )
     
 jloop:  do j=1,i
           lj =       gto(j)%lvalue
           oj =       gto(j)%offset
           nj =       gto(j)%cdegen
           ej =       gto(j)%exponent(1)
           xj = atom( gto(j)%origin )%coordinate_x
           yj = atom( gto(j)%origin )%coordinate_y
           zj = atom( gto(j)%origin )%coordinate_z
           cj =       gto(j)%coefficient(1)

                     fij = 1.0d0 
           if( i==j )fij = 0.5d0
                  
           dGij = dcmplx(0.0d0,0.0d0)

           if( doK )call get_dP_exchange( dP, nrow, oj, nj, dPj, dGj )
       
 kloop:    do k=1,nr_shells
             lk =       gto(k)%lvalue
             ok =       gto(k)%offset
             nk =       gto(k)%cdegen
             ek =       gto(k)%exponent(1)
             xk = atom( gto(k)%origin )%coordinate_x
             yk = atom( gto(k)%origin )%coordinate_y
             zk = atom( gto(k)%origin )%coordinate_z
             ck =       gto(k)%coefficient(1)
       
 lloop:      do l=1,k
               ll   =       gto(l)%lvalue
               nl   =       gto(l)%cdegen
               ol   =       gto(l)%offset
               el   =       gto(l)%exponent(1)
               xl   = atom( gto(l)%origin )%coordinate_x
               yl   = atom( gto(l)%origin )%coordinate_y
               zl   = atom( gto(l)%origin )%coordinate_z
               cl   =       gto(l)%coefficient(1)

                         fijkl = fij 
               if( k==l )fijkl = fij*0.5d0 
       
               SELECT CASE( eri_class )
                !CASE( 'llll' )
                !    eri_class_local = 'llll' 
                !    Fc = max( (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,l,k)), &
                !              (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,j,i))  )
                !    Fx = max( (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,j,k)), &
                !              (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,l,i))  )
                !    if( (Fc<screening_threshold).and.(Fx<screening_threshold) )cycle
                !CASE( 'ssll' )
                !    eri_class_local = 'ssll' 
                !    !> llll
                !    Fc = max( (gmax(1,i,j)*gmax(2,k,l)*pmax(2,2,l,k)),    & 
                !              (gmax(2,i,j)*gmax(1,k,l)*pmax(2,2,j,i)),    & 
                !              (gmax(2,i,j)*gmax(1,k,l)*pmax(1,1,l,k))*aa4,&  
                !              (gmax(1,i,j)*gmax(2,k,l)*pmax(1,1,j,i))*aa4 ) 
                !    Fx = max( (gmax(1,i,j)*gmax(2,k,l)*pmax(1,2,j,k)),    &
                !              (gmax(1,i,j)*gmax(2,k,l)*pmax(2,1,l,i)),    & 
                !              (gmax(2,i,j)*gmax(1,k,l)*pmax(2,1,j,k)),    & 
                !              (gmax(2,i,j)*gmax(1,k,l)*pmax(1,2,l,i))     )
                !    Fx=Fx*a2
                !    if( (Fc<screening_threshold).and.(Fx<screening_threshold) )eri_class_local='llll'
                !    !> exit
                !    Fc = max( (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,l,k)),&
                !              (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,j,i)) )
                !    Fx = max( (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,j,k)),&
                !              (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,l,i)) )
                !    if( (Fc<screening_threshold).and.(Fx<screening_threshold) )cycle
                 CASE( 'ssss' )
                     eri_class_local = 'ssss' 
                     !> ssll
                     Fc(1) = max( (gmax(2,i,j)*gmax(2,k,l)*pmax(2,2,l,k)),&
                                  (gmax(2,i,j)*gmax(2,k,l)*pmax(2,2,j,i)) ) 
                     Fx(1) = max( (gmax(2,i,j)*gmax(2,k,l)*pmax(2,2,j,k)),&
                                  (gmax(2,i,j)*gmax(2,k,l)*pmax(2,2,l,i)) ) 
                     Fc(1)=Fc(1)*aa4
                     Fx(1)=Fx(1)*aa4
                    !if( (Fc<screening_threshold).and.(Fx<screening_threshold) )eri_class_local='ssll'  
                     !> llll
                     Fc(2) = max( (gmax(1,i,j)*gmax(2,k,l)*pmax(2,2,l,k)),    & 
                                  (gmax(2,i,j)*gmax(1,k,l)*pmax(2,2,j,i)),    & 
                                  (gmax(2,i,j)*gmax(1,k,l)*pmax(1,1,l,k))*aa4,& 
                                  (gmax(1,i,j)*gmax(2,k,l)*pmax(1,1,j,i))*aa4 ) 
                     Fx(2) = max( (gmax(1,i,j)*gmax(2,k,l)*pmax(1,2,j,k)),    & 
                                  (gmax(1,i,j)*gmax(2,k,l)*pmax(2,1,l,i)),    & 
                                  (gmax(2,i,j)*gmax(1,k,l)*pmax(2,1,j,k)),    & 
                                  (gmax(2,i,j)*gmax(1,k,l)*pmax(1,2,l,i))     ) 
                     Fx(2)=Fx(2)*a2
                    !if( (Fc<screening_threshold).and.(Fx<screening_threshold) )eri_class_local='llll'
                     !> exit
                     Fc(3) = max( (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,l,k)),&
                                  (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,j,i)) )
                     Fx(3) = max( (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,j,k)),&
                                  (gmax(1,i,j)*gmax(1,k,l)*pmax(1,1,l,i)) )
                     if( (maxval(Fc)<screening_threshold).and.(maxval(Fx)<screening_threshold) )cycle
                 CASE DEFAULT; stop 'Error: the ERI class is not implemented yet => stop'
               END SELECT

               SELECT CASE( eri_class_local )
                 CASE( 'llll' )
                     write(6,*) 'stop: the LLLL-class not finished yet ...'
                     call interest_eri_rkb('llll',fijkl,gout,&
                                           li,ei,xi,yi,zi,ci,&
                                           lj,ej,xj,yj,zj,cj,&
                                           lk,ek,xk,yk,zk,ck,&
                                           ll,el,xl,yl,zl,cl )
                     call process_dG_llll_4x(ni,nj,nk,nl,oi,oj,ok,ol,                 &
                                             nrow,doK,gout,dGij,dPi,dPj,dGi,dGj,dP,dG )
                 CASE( 'ssll' )
                     write(6,*) 'stop: the SSLL-class not finished yet ...'
                     call interest_eri_rkb('ssll',fijkl,gout,&
                                           li,ei,xi,yi,zi,ci,&
                                           lj,ej,xj,yj,zj,cj,&
                                           lk,ek,xk,yk,zk,ck,&
                                           ll,el,xl,yl,zl,cl )
                     call process_dG_ssll_4x(ni,nj,nk,nl,oi,oj,ok,ol,                 &
                                             nrow,doK,gout,dGij,dPi,dPj,dGi,dGj,dP,dG )
                 CASE( 'ssss' ) 
                     call interest_eri_rkb('ssss',fijkl,gout,&
                                           li,ei,xi,yi,zi,ci,&
                                           lj,ej,xj,yj,zj,cj,&
                                           lk,ek,xk,yk,zk,ck,&
                                           ll,el,xl,yl,zl,cl )
                     call process_dG_ssss_4x(ni,nj,nk,nl,oi,oj,ok,ol,                 &
                                             nrow,doK,gout,dGij,dPi,dPj,dGi,dGj,dP,dG )
                 CASE DEFAULT; stop 'Error: the ERI class is not implemented yet => stop'
               END SELECT

             enddo lloop
           enddo kloop

                     call add_dG_coulomb ( dG, nrow, oi, ni, oj, nj, dGij        )
           if( doK ) call add_dG_exchange( dG, nrow,         oj, nj, dGj, hfxfac )

         enddo jloop

         if( doK ) call add_dG_exchange( dG, nrow, oi, ni, dGi, hfxfac )

       enddo iloop
    !> end

    deallocate( dGij )
    if( doK )then
      deallocate( dPi )
      deallocate( dPj )
      deallocate( dGi )
      deallocate( dGj )
    endif
  
    deallocate( gmax )
    deallocate( pmax )

    CONTAINS
  
    subroutine get_dP_exchange( dP, nbas, io, ic, dP_block, dG_block )
      !> input
      integer,    intent(in) :: io   
      integer,    intent(in) :: ic   
      integer,    intent(in) :: nbas 
      complex(8), intent(in) :: dP(4,4,nbas,nbas)
      !> output
      complex(8), intent(out) :: dP_block(4,4,ic,nbas)
      complex(8), intent(out) :: dG_block(4,4,ic,nbas)
      !> local
      integer :: ibas, i, j, k 
      !> 
      dG_block = dcmplx(0.0d0,0.0d0)
      do j=1,nbas
        do i=1,ic
          ibas = io+i
          do k=1,4
            dP_block(1,k,i,j) = dP(1,k,ibas,j)
            dP_block(2,k,i,j) = dP(2,k,ibas,j)
            dP_block(3,k,i,j) = dP(3,k,ibas,j)
            dP_block(4,k,i,j) = dP(4,k,ibas,j)
          enddo
        enddo
      enddo
    end subroutine

    subroutine add_dG_exchange( dG, nbas, io, ic, dG_block, scaleK )
      !> input
      integer,    intent(in) :: io   
      integer,    intent(in) :: ic   
      integer,    intent(in) :: nbas 
      real(8),    intent(in) :: scaleK
      complex(8), intent(in) :: dG_block(4,4,ic,nbas)
      !> output
      complex(8), intent(inout) :: dG(4,4,nbas,nbas)
      !> local
      integer :: ibas, i, j, k 
      !> 
      do j=1,nbas
        do i=1,ic
          ibas = io+i
          do k=1,4
            dG(1,k,ibas,j) = dG(1,k,ibas,j) + dG_block(1,k,i,j)*dcmplx(scaleK,0.0d0)
            dG(2,k,ibas,j) = dG(2,k,ibas,j) + dG_block(2,k,i,j)*dcmplx(scaleK,0.0d0)
            dG(3,k,ibas,j) = dG(3,k,ibas,j) + dG_block(3,k,i,j)*dcmplx(scaleK,0.0d0)
            dG(4,k,ibas,j) = dG(4,k,ibas,j) + dG_block(4,k,i,j)*dcmplx(scaleK,0.0d0)
          enddo
        enddo
      enddo
    end subroutine

    subroutine add_dG_coulomb( dG, nbas, io, ic, jo, jc, dG_block )
      !> input
      integer,    intent(in) :: nbas 
      integer,    intent(in) :: io, jo   
      integer,    intent(in) :: ic, jc   
      complex(8), intent(in) :: dG_block(4,4,ic,jc)
      !> input/output
      complex(8), intent(inout) :: dG(4,4,nbas,nbas)
      !> local
      integer :: ibas, jbas, i, j
      !> 
      do j=1,jc   
        jbas = jo+j
        do i=1,ic
          ibas = io+i
          dG(1,1,ibas,jbas) = dG(1,1,ibas,jbas) + dG_block(1,1,i,j)
          dG(3,3,ibas,jbas) = dG(3,3,ibas,jbas) + dG_block(3,3,i,j)
          dG(2,2,ibas,jbas) = dG(2,2,ibas,jbas) + dG_block(2,2,i,j)
          dG(4,2,ibas,jbas) = dG(4,2,ibas,jbas) + dG_block(4,2,i,j)
          dG(2,4,ibas,jbas) = dG(2,4,ibas,jbas) + dG_block(2,4,i,j)
          dG(4,4,ibas,jbas) = dG(4,4,ibas,jbas) + dG_block(4,4,i,j)
        enddo
      enddo
      !>
      do i=1,ic   
        ibas = io+i
        do j=1,jc
          jbas = jo+j
          dG(1,1,jbas,ibas) = dG(1,1,jbas,ibas) + dconjg(dG_block(1,1,i,j))
          dG(3,3,jbas,ibas) = dG(3,3,jbas,ibas) + dconjg(dG_block(3,3,i,j))
          dG(2,2,jbas,ibas) = dG(2,2,jbas,ibas) + dconjg(dG_block(2,2,i,j))
          dG(4,2,jbas,ibas) = dG(4,2,jbas,ibas) + dconjg(dG_block(2,4,i,j))
          dG(2,4,jbas,ibas) = dG(2,4,jbas,ibas) + dconjg(dG_block(4,2,i,j))
          dG(4,4,jbas,ibas) = dG(4,4,jbas,ibas) + dconjg(dG_block(4,4,i,j))
        enddo
      enddo
    end subroutine

    subroutine process_dG_ssss_4x(ic,jc,kc,lc,io,jo,ko,lo,nbas,doK,gout,dGij,dPi,dPj,dGi,dGj,dP,dG) 
      !> input
      logical,    intent(in) :: doK
      integer,    intent(in) :: nbas
      integer,    intent(in) :: ic, jc, kc, lc
      integer,    intent(in) :: io, jo, ko, lo
      complex(8), intent(in) :: dPi(4,4,ic,nbas)
      complex(8), intent(in) :: dPj(4,4,jc,nbas)
      complex(8), intent(in) :: dP (4,4,nbas,nbas)
      real(8),    intent(in) :: gout(5,kc,lc,5,ic,jc)
      !> input/output
      complex(8), intent(inout) :: dGij(4,4,ic,jc)
      complex(8), intent(inout) :: dGi (4,4,ic,nbas)
      complex(8), intent(inout) :: dGj (4,4,jc,nbas)
      complex(8), intent(inout) :: dG  (4,4,nbas,nbas)
      !> local
      integer    :: i, j, k, l, m
      integer    :: ibas, jbas, kbas, lbas
      complex(8) :: pll
      complex(8) :: paa, pab, pba, pbb
      complex(8) :: gg, fgg, fgaag, fgabg
      complex(8) :: gaag, gabg, gbag, gbbg
      complex(8) :: ggaa, ggab, ggba, ggbb   
      complex(8) :: gaagaa, gabgaa, gbagaa, gbbgaa
      complex(8) :: gaagab, gabgab, gbagab, gbbgab
      complex(8) :: gaagba, gabgba, gbagba, gbbgba
      complex(8) :: gaagbb, gabgbb, gbagbb, gbbgbb
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
                 gaagbb = dcmplx( ( gout(2,k,l,2,i,j)+gout(3,k,l,3,i,j)), (-gout(3,k,l,2,i,j)+gout(2,k,l,3,i,j)) )
                 gaagab = dcmplx( ( gout(4,k,l,2,i,j)-gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)+gout(4,k,l,3,i,j)) )
                 gaagba = dcmplx( (-gout(4,k,l,2,i,j)-gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)-gout(4,k,l,3,i,j)) )
                 gabgaa = dcmplx( ( gout(2,k,l,4,i,j)-gout(3,k,l,5,i,j)), ( gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
                 gabgbb = dcmplx( ( gout(2,k,l,4,i,j)+gout(3,k,l,5,i,j)), (-gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
                 gabgab = dcmplx( ( gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), ( gout(5,k,l,4,i,j)+gout(4,k,l,5,i,j)) )
                 gabgba = dcmplx( (-gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), ( gout(5,k,l,4,i,j)-gout(4,k,l,5,i,j)) )
              !> end
              !> contraction [ij] = [ij|kl]*lk = [ij|kl]*dconjg(kl)
                 paa = dconjg(dP(2,2,kbas,lbas)                    )  !=dP(2,2,lbas,kbas)
                 pab = dconjg(dP(4,2,kbas,lbas)                    )  !=dP(2,4,lbas,kbas)
                 pba = dconjg(dP(2,4,kbas,lbas)                    )  !=dP(4,2,lbas,kbas)
                 pbb = dconjg(dP(4,4,kbas,lbas)                    )  !=dP(4,4,lbas,kbas)
                 pll = dconjg(dP(1,1,kbas,lbas) + dP(3,3,kbas,lbas))  !=dP(1,1,lbas,kbas) + dP(3,3,lbas,kbas)
                 !>
                 fgg   =   gg*pll +   ggaa*paa +   ggab*pba +   ggba*pab +   ggbb*pbb
                 fgaag = gaag*pll + gaagaa*paa + gaagab*pba + gaagba*pab + gaagbb*pbb
                 fgabg = gabg*pll + gabgaa*paa + gabgab*pba + gabgba*pab + gabgbb*pbb
                 !>
                 dGij(1,1,i,j) = dGij(1,1,i,j) +        fgg
                 dGij(3,3,i,j) = dGij(3,3,i,j) + dconjg(fgg)
                 dGij(2,2,i,j) = dGij(2,2,i,j) +        fgaag
                 dGij(4,2,i,j) = dGij(4,2,i,j) - dconjg(fgabg)
                 dGij(2,4,i,j) = dGij(2,4,i,j) +        fgabg
                 dGij(4,4,i,j) = dGij(4,4,i,j) + dconjg(fgaag)
              !> end
              !> contraction [ji] = [ji|kl]*lk = dconjg([ij]) - transfered to add_dG_coulomb
              !  dG(1,1,jbas,ibas) = dG(1,1,jbas,ibas) + dconjg(fgg)
              !  dG(3,3,jbas,ibas) = dG(3,3,jbas,ibas) +        fgg    
              !  dG(2,2,jbas,ibas) = dG(2,2,jbas,ibas) + dconjg(fgaag) 
              !  dG(4,2,jbas,ibas) = dG(4,2,jbas,ibas) + dconjg(fgabg) 
              !  dG(2,4,jbas,ibas) = dG(2,4,jbas,ibas) -        fgabg
              !  dG(4,4,jbas,ibas) = dG(4,4,jbas,ibas) +        fgaag 
              !> end
              !> contraction [ij] = [ij|lk]*kl
                 fgg = pba
                 paa = dconjg(paa)        !=dP(2,2,kbas,lbas)                    
                 pba = dconjg(pab)        !=dP(4,2,kbas,lbas)                    
                 pab = dconjg(fgg)        !=dP(2,4,kbas,lbas)                    
                 pbb = dconjg(pbb)        !=dP(4,4,kbas,lbas)                    
                 pll = dconjg(pll)        !=dP(1,1,kbas,lbas) + dP(3,3,kbas,lbas)
                 !>
                 fgg   =   gg*pll +   ggbb*paa -   ggab*pba -   ggba*pab +   ggaa*pbb
                 fgaag = gaag*pll + gaagbb*paa - gaagab*pba - gaagba*pab + gaagaa*pbb
                 fgabg = gabg*pll + gabgbb*paa - gabgab*pba - gabgba*pab + gabgaa*pbb
                 !>
                 dGij(1,1,i,j) = dGij(1,1,i,j) +        fgg
                 dGij(3,3,i,j) = dGij(3,3,i,j) + dconjg(fgg)
                 dGij(2,2,i,j) = dGij(2,2,i,j) +        fgaag
                 dGij(4,2,i,j) = dGij(4,2,i,j) - dconjg(fgabg)
                 dGij(2,4,i,j) = dGij(2,4,i,j) +        fgabg
                 dGij(4,4,i,j) = dGij(4,4,i,j) + dconjg(fgaag)
              !> end
              !> contraction [ji] = [ji|lk]*kl = [ji|lk]*dconjg(lk) = dconjg([ij]) - transfered to add_dG_coulomb
              !  dG(1,1,jbas,ibas) = dG(1,1,jbas,ibas) + dconjg(fgg  )
              !  dG(3,3,jbas,ibas) = dG(3,3,jbas,ibas) +        fgg   
              !  dG(2,2,jbas,ibas) = dG(2,2,jbas,ibas) + dconjg(fgaag) 
              !  dG(4,2,jbas,ibas) = dG(4,2,jbas,ibas) + dconjg(fgabg) 
              !  dG(2,4,jbas,ibas) = dG(2,4,jbas,ibas) -        fgabg
              !  dG(4,4,jbas,ibas) = dG(4,4,jbas,ibas) +        fgaag 
              !> end
            enddo
          enddo
        enddo
      enddo
      !> exchange part
      if( .not. doK )return
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          do j=1,jc
            do i=1,ic
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
              !> end
              !> contraction [ij|kl]*jk
                 paa_ll = dPj(1,1,j,kbas)
                 pba_ll = dPj(3,1,j,kbas)
                 pab_ll = dPj(1,3,j,kbas)
                 pbb_ll = dPj(3,3,j,kbas)
                 paa_ls = dPj(1,2,j,kbas)
                 pba_ls = dPj(3,2,j,kbas)
                 pab_ls = dPj(1,4,j,kbas)
                 pbb_ls = dPj(3,4,j,kbas)
                 paa_sl = dPj(2,1,j,kbas)
                 pba_sl = dPj(4,1,j,kbas)
                 pab_sl = dPj(2,3,j,kbas)
                 pbb_sl = dPj(4,3,j,kbas)
                 paa_ss = dPj(2,2,j,kbas)
                 pba_ss = dPj(4,2,j,kbas)
                 pab_ss = dPj(2,4,j,kbas)
                 pbb_ss = dPj(4,4,j,kbas)
                 !> large-large block
                 dGi(1,1,i,lbas) = dGi(1,1,i,lbas) - paa_ll*gg
                 dGi(1,3,i,lbas) = dGi(1,3,i,lbas) - pab_ll*gg
                 dGi(3,1,i,lbas) = dGi(3,1,i,lbas) - pba_ll*gg
                 dGi(3,3,i,lbas) = dGi(3,3,i,lbas) - pbb_ll*gg
                 !> large-small block
                 dGi(1,2,i,lbas) = dGi(1,2,i,lbas) - paa_ls*ggaa - pab_ls*ggba  
                 dGi(1,4,i,lbas) = dGi(1,4,i,lbas) - paa_ls*ggab - pab_ls*ggbb
                 dGi(3,2,i,lbas) = dGi(3,2,i,lbas) - pba_ls*ggaa - pbb_ls*ggba
                 dGi(3,4,i,lbas) = dGi(3,4,i,lbas) - pba_ls*ggab - pbb_ls*ggbb
                 !> small-large block
                 dGi(2,1,i,lbas) = dGi(2,1,i,lbas) - gaag*paa_sl - gabg*pba_sl 
                 dGi(2,3,i,lbas) = dGi(2,3,i,lbas) - gaag*pab_sl - gabg*pbb_sl 
                 dGi(4,1,i,lbas) = dGi(4,1,i,lbas) - gbag*paa_sl - gbbg*pba_sl 
                 dGi(4,3,i,lbas) = dGi(4,3,i,lbas) - gbag*pab_sl - gbbg*pbb_sl 
                 !> small-small block
                 dGi(2,2,i,lbas) = dGi(2,2,i,lbas) - gaagaa*paa_ss - gaagba*pab_ss - gabgaa*pba_ss - gabgba*pbb_ss 
                 dGi(2,4,i,lbas) = dGi(2,4,i,lbas) - gaagab*paa_ss - gaagbb*pab_ss - gabgab*pba_ss - gabgbb*pbb_ss 
                 dGi(4,2,i,lbas) = dGi(4,2,i,lbas) - gbagaa*paa_ss - gbagba*pab_ss - gbbgaa*pba_ss - gbbgba*pbb_ss 
                 dGi(4,4,i,lbas) = dGi(4,4,i,lbas) - gbagab*paa_ss - gbagbb*pab_ss - gbbgab*pba_ss - gbbgbb*pbb_ss 
              !> end
              !> contraction [ij|lk]*jl
                 paa_ll = dPj(1,1,j,lbas)
                 pba_ll = dPj(3,1,j,lbas)
                 pab_ll = dPj(1,3,j,lbas)
                 pbb_ll = dPj(3,3,j,lbas)
                 paa_ls = dPj(1,2,j,lbas)
                 pba_ls = dPj(3,2,j,lbas)
                 pab_ls = dPj(1,4,j,lbas)
                 pbb_ls = dPj(3,4,j,lbas)
                 paa_sl = dPj(2,1,j,lbas)
                 pba_sl = dPj(4,1,j,lbas)
                 pab_sl = dPj(2,3,j,lbas)
                 pbb_sl = dPj(4,3,j,lbas)
                 paa_ss = dPj(2,2,j,lbas)
                 pba_ss = dPj(4,2,j,lbas)
                 pab_ss = dPj(2,4,j,lbas)
                 pbb_ss = dPj(4,4,j,lbas)
                 !> large-large block
                 dGi(1,1,i,kbas) = dGi(1,1,i,kbas) - paa_ll*gg
                 dGi(1,3,i,kbas) = dGi(1,3,i,kbas) - pab_ll*gg
                 dGi(3,1,i,kbas) = dGi(3,1,i,kbas) - pba_ll*gg
                 dGi(3,3,i,kbas) = dGi(3,3,i,kbas) - pbb_ll*gg
                 !> large-small block
                 dGi(1,2,i,kbas) = dGi(1,2,i,kbas) - paa_ls*ggbb + pab_ls*ggba  
                 dGi(1,4,i,kbas) = dGi(1,4,i,kbas) + paa_ls*ggab - pab_ls*ggaa
                 dGi(3,2,i,kbas) = dGi(3,2,i,kbas) - pba_ls*ggbb + pbb_ls*ggba
                 dGi(3,4,i,kbas) = dGi(3,4,i,kbas) + pba_ls*ggab - pbb_ls*ggaa
                 !> small-large block
                 dGi(2,1,i,kbas) = dGi(2,1,i,kbas) - gaag*paa_sl - gabg*pba_sl 
                 dGi(2,3,i,kbas) = dGi(2,3,i,kbas) - gaag*pab_sl - gabg*pbb_sl 
                 dGi(4,1,i,kbas) = dGi(4,1,i,kbas) - gbag*paa_sl - gbbg*pba_sl 
                 dGi(4,3,i,kbas) = dGi(4,3,i,kbas) - gbag*pab_sl - gbbg*pbb_sl 
                 !> small-small block
                 dGi(2,2,i,kbas) = dGi(2,2,i,kbas) - gaagbb*paa_ss + gaagba*pab_ss - gabgbb*pba_ss + gabgba*pbb_ss 
                 dGi(2,4,i,kbas) = dGi(2,4,i,kbas) + gaagab*paa_ss - gaagaa*pab_ss + gabgab*pba_ss - gabgaa*pbb_ss 
                 dGi(4,2,i,kbas) = dGi(4,2,i,kbas) - gbagbb*paa_ss + gbagba*pab_ss - gbbgbb*pba_ss + gbbgba*pbb_ss 
                 dGi(4,4,i,kbas) = dGi(4,4,i,kbas) + gbagab*paa_ss - gbagaa*pab_ss + gbbgab*pba_ss - gbbgaa*pbb_ss 
              !> end
              !> contraction [ji|kl]*ik
                 paa_ll = dPi(1,1,i,kbas)
                 pba_ll = dPi(3,1,i,kbas)
                 pab_ll = dPi(1,3,i,kbas)
                 pbb_ll = dPi(3,3,i,kbas)
                 paa_ls = dPi(1,2,i,kbas)
                 pba_ls = dPi(3,2,i,kbas)
                 pab_ls = dPi(1,4,i,kbas)
                 pbb_ls = dPi(3,4,i,kbas)
                 paa_sl = dPi(2,1,i,kbas)
                 pba_sl = dPi(4,1,i,kbas)
                 pab_sl = dPi(2,3,i,kbas)
                 pbb_sl = dPi(4,3,i,kbas)
                 paa_ss = dPi(2,2,i,kbas)
                 pba_ss = dPi(4,2,i,kbas)
                 pab_ss = dPi(2,4,i,kbas)
                 pbb_ss = dPi(4,4,i,kbas)
                 !> large-large block
                 dGj(1,1,j,lbas) = dGj(1,1,j,lbas) - paa_ll*gg
                 dGj(1,3,j,lbas) = dGj(1,3,j,lbas) - pab_ll*gg
                 dGj(3,1,j,lbas) = dGj(3,1,j,lbas) - pba_ll*gg
                 dGj(3,3,j,lbas) = dGj(3,3,j,lbas) - pbb_ll*gg
                 !> large-small block
                 dGj(1,2,j,lbas) = dGj(1,2,j,lbas) - paa_ls*ggaa - pab_ls*ggba  
                 dGj(1,4,j,lbas) = dGj(1,4,j,lbas) - paa_ls*ggab - pab_ls*ggbb
                 dGj(3,2,j,lbas) = dGj(3,2,j,lbas) - pba_ls*ggaa - pbb_ls*ggba
                 dGj(3,4,j,lbas) = dGj(3,4,j,lbas) - pba_ls*ggab - pbb_ls*ggbb
                 !> small-large block
                 dGj(2,1,j,lbas) = dGj(2,1,j,lbas) - gbbg*paa_sl + gabg*pba_sl 
                 dGj(2,3,j,lbas) = dGj(2,3,j,lbas) - gbbg*pab_sl + gabg*pbb_sl 
                 dGj(4,1,j,lbas) = dGj(4,1,j,lbas) + gbag*paa_sl - gaag*pba_sl 
                 dGj(4,3,j,lbas) = dGj(4,3,j,lbas) + gbag*pab_sl - gaag*pbb_sl 
                 !> small-small block
                 dGj(2,2,j,lbas) = dGj(2,2,j,lbas) - gbbgaa*paa_ss - gbbgba*pab_ss + gabgaa*pba_ss + gabgba*pbb_ss 
                 dGj(2,4,j,lbas) = dGj(2,4,j,lbas) - gbbgab*paa_ss - gbbgbb*pab_ss + gabgab*pba_ss + gabgbb*pbb_ss 
                 dGj(4,2,j,lbas) = dGj(4,2,j,lbas) + gbagaa*paa_ss + gbagba*pab_ss - gaagaa*pba_ss - gaagba*pbb_ss 
                 dGj(4,4,j,lbas) = dGj(4,4,j,lbas) + gbagab*paa_ss + gbagbb*pab_ss - gaagab*pba_ss - gaagbb*pbb_ss 
              !> end
              !> contraction [ji|lk]*il
                 paa_ll = dPi(1,1,i,lbas)
                 pba_ll = dPi(3,1,i,lbas)
                 pab_ll = dPi(1,3,i,lbas)
                 pbb_ll = dPi(3,3,i,lbas)
                 paa_ls = dPi(1,2,i,lbas)
                 pba_ls = dPi(3,2,i,lbas)
                 pab_ls = dPi(1,4,i,lbas)
                 pbb_ls = dPi(3,4,i,lbas)
                 paa_sl = dPi(2,1,i,lbas)
                 pba_sl = dPi(4,1,i,lbas)
                 pab_sl = dPi(2,3,i,lbas)
                 pbb_sl = dPi(4,3,i,lbas)
                 paa_ss = dPi(2,2,i,lbas)
                 pba_ss = dPi(4,2,i,lbas)
                 pab_ss = dPi(2,4,i,lbas)
                 pbb_ss = dPi(4,4,i,lbas)
                 !> large-large block
                 dGj(1,1,j,kbas) = dGj(1,1,j,kbas) - paa_ll*gg
                 dGj(1,3,j,kbas) = dGj(1,3,j,kbas) - pab_ll*gg
                 dGj(3,1,j,kbas) = dGj(3,1,j,kbas) - pba_ll*gg
                 dGj(3,3,j,kbas) = dGj(3,3,j,kbas) - pbb_ll*gg
                 !> large-small block
                 dGj(1,2,j,kbas) = dGj(1,2,j,kbas) - paa_ls*ggbb + pab_ls*ggba  
                 dGj(1,4,j,kbas) = dGj(1,4,j,kbas) + paa_ls*ggab - pab_ls*ggaa
                 dGj(3,2,j,kbas) = dGj(3,2,j,kbas) - pba_ls*ggbb + pbb_ls*ggba
                 dGj(3,4,j,kbas) = dGj(3,4,j,kbas) + pba_ls*ggab - pbb_ls*ggaa
                 !> small-large block
                 dGj(2,1,j,kbas) = dGj(2,1,j,kbas) - gbbg*paa_sl + gabg*pba_sl 
                 dGj(2,3,j,kbas) = dGj(2,3,j,kbas) - gbbg*pab_sl + gabg*pbb_sl 
                 dGj(4,1,j,kbas) = dGj(4,1,j,kbas) + gbag*paa_sl - gaag*pba_sl 
                 dGj(4,3,j,kbas) = dGj(4,3,j,kbas) + gbag*pab_sl - gaag*pbb_sl 
                 !> small-small block
                 dGj(2,2,j,kbas) = dGj(2,2,j,kbas) - gbbgbb*paa_ss + gbbgba*pab_ss + gabgbb*pba_ss - gabgba*pbb_ss 
                 dGj(2,4,j,kbas) = dGj(2,4,j,kbas) + gbbgab*paa_ss - gbbgaa*pab_ss - gabgab*pba_ss + gabgaa*pbb_ss 
                 dGj(4,2,j,kbas) = dGj(4,2,j,kbas) + gbagbb*paa_ss - gbagba*pab_ss - gaagbb*pba_ss + gaagba*pbb_ss 
                 dGj(4,4,j,kbas) = dGj(4,4,j,kbas) - gbagab*paa_ss + gbagaa*pab_ss + gaagab*pba_ss - gaagaa*pbb_ss 
              !> end
            enddo
          enddo
        enddo
      enddo
    end subroutine

    subroutine process_dG_ssll_4x(ic,jc,kc,lc,io,jo,ko,lo,nbas,doK,gout,dGij,dPi,dPj,dGi,dGj,dP,dG) 
      !> input
      logical,    intent(in) :: doK
      integer,    intent(in) :: nbas
      integer,    intent(in) :: ic, jc, kc, lc
      integer,    intent(in) :: io, jo, ko, lo
      complex(8), intent(in) :: dPi(4,4,ic,nbas)
      complex(8), intent(in) :: dPj(4,4,jc,nbas)
      complex(8), intent(in) :: dP (4,4,nbas,nbas)
      real(8),    intent(in) :: gout(5,kc,lc,5,ic,jc)
      !> input/output
      complex(8), intent(inout) :: dGij(4,4,ic,jc)
      complex(8), intent(inout) :: dGi (4,4,ic,nbas)
      complex(8), intent(inout) :: dGj (4,4,jc,nbas)
      complex(8), intent(inout) :: dG  (4,4,nbas,nbas)
      !> local
      integer    :: i, j, k, l, m
      integer    :: ibas, jbas, kbas, lbas
      complex(8) :: pll
      complex(8) :: paa, pab, pba, pbb
      complex(8) :: gg, fgg, fgaag, fgabg
      complex(8) :: gaag, gabg, gbag, gbbg
      complex(8) :: ggaa, ggab, ggba, ggbb   
      complex(8) :: gaagaa, gabgaa, gbagaa, gbbgaa
      complex(8) :: gaagab, gabgab, gbagab, gbbgab
      complex(8) :: gaagba, gabgba, gbagba, gbbgba
      complex(8) :: gaagbb, gabgbb, gbagbb, gbbgbb
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
              !> end
              !> contraction [ij] = [ij|kl]*lk = [ij|kl]*dconjg(kl)
                 paa = dconjg(dP(2,2,kbas,lbas)                    )  !=dP(2,2,lbas,kbas)
                 pab = dconjg(dP(4,2,kbas,lbas)                    )  !=dP(2,4,lbas,kbas)
                 pba = dconjg(dP(2,4,kbas,lbas)                    )  !=dP(4,2,lbas,kbas)
                 pbb = dconjg(dP(4,4,kbas,lbas)                    )  !=dP(4,4,lbas,kbas)
                 pll = dconjg(dP(1,1,kbas,lbas) + dP(3,3,kbas,lbas))  !=dP(1,1,lbas,kbas) + dP(3,3,lbas,kbas)
                 !>
                 fgg   =   gg*pll + ggaa*paa + ggab*pba + ggba*pab + ggbb*pbb
                 fgaag = gaag*pll
                 fgabg = gabg*pll
                 !>
                 dGij(1,1,i,j) = dGij(1,1,i,j) +        fgg
                 dGij(3,3,i,j) = dGij(3,3,i,j) + dconjg(fgg)
                 dGij(2,2,i,j) = dGij(2,2,i,j) +        fgaag
                 dGij(4,2,i,j) = dGij(4,2,i,j) - dconjg(fgabg)
                 dGij(2,4,i,j) = dGij(2,4,i,j) +        fgabg
                 dGij(4,4,i,j) = dGij(4,4,i,j) + dconjg(fgaag)
              !> end
              !> contraction [ij] = [ij|lk]*kl
                 fgg = pba
                 paa = dconjg(paa)        !=dP(2,2,kbas,lbas)                    
                 pba = dconjg(pab)        !=dP(4,2,kbas,lbas)                    
                 pab = dconjg(fgg)        !=dP(2,4,kbas,lbas)                    
                 pbb = dconjg(pbb)        !=dP(4,4,kbas,lbas)                    
                 pll = dconjg(pll)        !=dP(1,1,kbas,lbas) + dP(3,3,kbas,lbas)
                 !>
                 fgg   =   gg*pll + ggbb*paa - ggab*pba - ggba*pab + ggaa*pbb
                 fgaag = gaag*pll
                 fgabg = gabg*pll
                 !>
                 dGij(1,1,i,j) = dGij(1,1,i,j) +        fgg
                 dGij(3,3,i,j) = dGij(3,3,i,j) + dconjg(fgg)
                 dGij(2,2,i,j) = dGij(2,2,i,j) +        fgaag
                 dGij(4,2,i,j) = dGij(4,2,i,j) - dconjg(fgabg)
                 dGij(2,4,i,j) = dGij(2,4,i,j) +        fgabg
                 dGij(4,4,i,j) = dGij(4,4,i,j) + dconjg(fgaag)
              !> end
            enddo
          enddo
        enddo
      enddo
      !> exchange part
      if( .not. doK )return
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          do j=1,jc
            do i=1,ic
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
              !> end
              !> contraction [ij|kl]*jk
                 paa_ll = dPj(1,1,j,kbas)
                 pba_ll = dPj(3,1,j,kbas)
                 pab_ll = dPj(1,3,j,kbas)
                 pbb_ll = dPj(3,3,j,kbas)
                 paa_ls = dPj(1,2,j,kbas)
                 pba_ls = dPj(3,2,j,kbas)
                 pab_ls = dPj(1,4,j,kbas)
                 pbb_ls = dPj(3,4,j,kbas)
                 paa_sl = dPj(2,1,j,kbas)
                 pba_sl = dPj(4,1,j,kbas)
                 pab_sl = dPj(2,3,j,kbas)
                 pbb_sl = dPj(4,3,j,kbas)
                 !> large-large block
                 dGi(1,1,i,lbas) = dGi(1,1,i,lbas) - paa_ll*gg
                 dGi(1,3,i,lbas) = dGi(1,3,i,lbas) - pab_ll*gg
                 dGi(3,1,i,lbas) = dGi(3,1,i,lbas) - pba_ll*gg
                 dGi(3,3,i,lbas) = dGi(3,3,i,lbas) - pbb_ll*gg
                 !> large-small block
                 dGi(1,2,i,lbas) = dGi(1,2,i,lbas) - paa_ls*ggaa - pab_ls*ggba  
                 dGi(1,4,i,lbas) = dGi(1,4,i,lbas) - paa_ls*ggab - pab_ls*ggbb
                 dGi(3,2,i,lbas) = dGi(3,2,i,lbas) - pba_ls*ggaa - pbb_ls*ggba
                 dGi(3,4,i,lbas) = dGi(3,4,i,lbas) - pba_ls*ggab - pbb_ls*ggbb
                 !> small-large block
                 dGi(2,1,i,lbas) = dGi(2,1,i,lbas) - gaag*paa_sl - gabg*pba_sl 
                 dGi(2,3,i,lbas) = dGi(2,3,i,lbas) - gaag*pab_sl - gabg*pbb_sl 
                 dGi(4,1,i,lbas) = dGi(4,1,i,lbas) - gbag*paa_sl - gbbg*pba_sl 
                 dGi(4,3,i,lbas) = dGi(4,3,i,lbas) - gbag*pab_sl - gbbg*pbb_sl 
              !> end
              !> contraction [ij|lk]*jl
                 paa_ll = dPj(1,1,j,lbas)
                 pba_ll = dPj(3,1,j,lbas)
                 pab_ll = dPj(1,3,j,lbas)
                 pbb_ll = dPj(3,3,j,lbas)
                 paa_ls = dPj(1,2,j,lbas)
                 pba_ls = dPj(3,2,j,lbas)
                 pab_ls = dPj(1,4,j,lbas)
                 pbb_ls = dPj(3,4,j,lbas)
                 paa_sl = dPj(2,1,j,lbas)
                 pba_sl = dPj(4,1,j,lbas)
                 pab_sl = dPj(2,3,j,lbas)
                 pbb_sl = dPj(4,3,j,lbas)
                 !> large-large block
                 dGi(1,1,i,kbas) = dGi(1,1,i,kbas) - paa_ll*gg
                 dGi(1,3,i,kbas) = dGi(1,3,i,kbas) - pab_ll*gg
                 dGi(3,1,i,kbas) = dGi(3,1,i,kbas) - pba_ll*gg
                 dGi(3,3,i,kbas) = dGi(3,3,i,kbas) - pbb_ll*gg
                 !> large-small block
                 dGi(1,2,i,kbas) = dGi(1,2,i,kbas) - paa_ls*ggbb + pab_ls*ggba  
                 dGi(1,4,i,kbas) = dGi(1,4,i,kbas) + paa_ls*ggab - pab_ls*ggaa
                 dGi(3,2,i,kbas) = dGi(3,2,i,kbas) - pba_ls*ggbb + pbb_ls*ggba
                 dGi(3,4,i,kbas) = dGi(3,4,i,kbas) + pba_ls*ggab - pbb_ls*ggaa
                 !> small-large block
                 dGi(2,1,i,kbas) = dGi(2,1,i,kbas) - gaag*paa_sl - gabg*pba_sl 
                 dGi(2,3,i,kbas) = dGi(2,3,i,kbas) - gaag*pab_sl - gabg*pbb_sl 
                 dGi(4,1,i,kbas) = dGi(4,1,i,kbas) - gbag*paa_sl - gbbg*pba_sl 
                 dGi(4,3,i,kbas) = dGi(4,3,i,kbas) - gbag*pab_sl - gbbg*pbb_sl 
              !> end
              !> contraction [ji|kl]*ik
                 paa_ll = dPi(1,1,i,kbas)
                 pba_ll = dPi(3,1,i,kbas)
                 pab_ll = dPi(1,3,i,kbas)
                 pbb_ll = dPi(3,3,i,kbas)
                 paa_ls = dPi(1,2,i,kbas)
                 pba_ls = dPi(3,2,i,kbas)
                 pab_ls = dPi(1,4,i,kbas)
                 pbb_ls = dPi(3,4,i,kbas)
                 paa_sl = dPi(2,1,i,kbas)
                 pba_sl = dPi(4,1,i,kbas)
                 pab_sl = dPi(2,3,i,kbas)
                 pbb_sl = dPi(4,3,i,kbas)
                 !> large-large block
                 dGj(1,1,j,lbas) = dGj(1,1,j,lbas) - paa_ll*gg
                 dGj(1,3,j,lbas) = dGj(1,3,j,lbas) - pab_ll*gg
                 dGj(3,1,j,lbas) = dGj(3,1,j,lbas) - pba_ll*gg
                 dGj(3,3,j,lbas) = dGj(3,3,j,lbas) - pbb_ll*gg
                 !> large-small block
                 dGj(1,2,j,lbas) = dGj(1,2,j,lbas) - paa_ls*ggaa - pab_ls*ggba  
                 dGj(1,4,j,lbas) = dGj(1,4,j,lbas) - paa_ls*ggab - pab_ls*ggbb
                 dGj(3,2,j,lbas) = dGj(3,2,j,lbas) - pba_ls*ggaa - pbb_ls*ggba
                 dGj(3,4,j,lbas) = dGj(3,4,j,lbas) - pba_ls*ggab - pbb_ls*ggbb
                 !> small-large block
                 dGj(2,1,j,lbas) = dGj(2,1,j,lbas) - gbbg*paa_sl + gabg*pba_sl 
                 dGj(2,3,j,lbas) = dGj(2,3,j,lbas) - gbbg*pab_sl + gabg*pbb_sl 
                 dGj(4,1,j,lbas) = dGj(4,1,j,lbas) + gbag*paa_sl - gaag*pba_sl 
                 dGj(4,3,j,lbas) = dGj(4,3,j,lbas) + gbag*pab_sl - gaag*pbb_sl 
              !> end
              !> contraction [ji|lk]*il
                 paa_ll = dPi(1,1,i,lbas)
                 pba_ll = dPi(3,1,i,lbas)
                 pab_ll = dPi(1,3,i,lbas)
                 pbb_ll = dPi(3,3,i,lbas)
                 paa_ls = dPi(1,2,i,lbas)
                 pba_ls = dPi(3,2,i,lbas)
                 pab_ls = dPi(1,4,i,lbas)
                 pbb_ls = dPi(3,4,i,lbas)
                 paa_sl = dPi(2,1,i,lbas)
                 pba_sl = dPi(4,1,i,lbas)
                 pab_sl = dPi(2,3,i,lbas)
                 pbb_sl = dPi(4,3,i,lbas)
                 !> large-large block
                 dGj(1,1,j,kbas) = dGj(1,1,j,kbas) - paa_ll*gg
                 dGj(1,3,j,kbas) = dGj(1,3,j,kbas) - pab_ll*gg
                 dGj(3,1,j,kbas) = dGj(3,1,j,kbas) - pba_ll*gg
                 dGj(3,3,j,kbas) = dGj(3,3,j,kbas) - pbb_ll*gg
                 !> large-small block
                 dGj(1,2,j,kbas) = dGj(1,2,j,kbas) - paa_ls*ggbb + pab_ls*ggba  
                 dGj(1,4,j,kbas) = dGj(1,4,j,kbas) + paa_ls*ggab - pab_ls*ggaa
                 dGj(3,2,j,kbas) = dGj(3,2,j,kbas) - pba_ls*ggbb + pbb_ls*ggba
                 dGj(3,4,j,kbas) = dGj(3,4,j,kbas) + pba_ls*ggab - pbb_ls*ggaa
                 !> small-large block
                 dGj(2,1,j,kbas) = dGj(2,1,j,kbas) - gbbg*paa_sl + gabg*pba_sl 
                 dGj(2,3,j,kbas) = dGj(2,3,j,kbas) - gbbg*pab_sl + gabg*pbb_sl 
                 dGj(4,1,j,kbas) = dGj(4,1,j,kbas) + gbag*paa_sl - gaag*pba_sl 
                 dGj(4,3,j,kbas) = dGj(4,3,j,kbas) + gbag*pab_sl - gaag*pbb_sl 
              !> end
            enddo
          enddo
        enddo
      enddo
    end subroutine

    subroutine process_dG_llll_4x(ic,jc,kc,lc,io,jo,ko,lo,nbas,doK,gout,dGij,dPi,dPj,dGi,dGj,dP,dG) 
      !> input
      logical,    intent(in) :: doK
      integer,    intent(in) :: nbas
      integer,    intent(in) :: ic, jc, kc, lc
      integer,    intent(in) :: io, jo, ko, lo
      complex(8), intent(in) :: dPi(4,4,ic,nbas)
      complex(8), intent(in) :: dPj(4,4,jc,nbas)
      complex(8), intent(in) :: dP (4,4,nbas,nbas)
      real(8),    intent(in) :: gout(kc,lc,ic,jc)
      !> input/output
      complex(8), intent(inout) :: dGij(4,4,ic,jc)
      complex(8), intent(inout) :: dGi (4,4,ic,nbas)
      complex(8), intent(inout) :: dGj (4,4,jc,nbas)
      complex(8), intent(inout) :: dG  (4,4,nbas,nbas)
      !> local
      integer    :: i, j, k, l, m
      integer    :: ibas, jbas, kbas, lbas
      complex(8) :: pll
      complex(8) :: paa, pab, pba, pbb
      complex(8) :: gg, fgg, fgaag, fgabg
      complex(8) :: gaag, gabg, gbag, gbbg
      complex(8) :: ggaa, ggab, ggba, ggbb   
      complex(8) :: gaagaa, gabgaa, gbagaa, gbbgaa
      complex(8) :: gaagab, gabgab, gbagab, gbbgab
      complex(8) :: gaagba, gabgba, gbagba, gbbgba
      complex(8) :: gaagbb, gabgbb, gbagbb, gbbgbb
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
                 gg = dcmplx( gout(k,l,i,j), 0.0d0 )
              !> end
              !> contraction [ij] = [ij|kl]*lk = [ij|kl]*dconjg(kl)
                 pll = dconjg(dP(1,1,kbas,lbas) + dP(3,3,kbas,lbas))  !=dP(1,1,lbas,kbas) + dP(3,3,lbas,kbas)
                 !>
                 fgg = gg*pll 
                 !>
                 dGij(1,1,i,j) = dGij(1,1,i,j) +        fgg
                 dGij(3,3,i,j) = dGij(3,3,i,j) + dconjg(fgg)
              !> end
              !> contraction [ij] = [ij|lk]*kl
                 pll = dconjg(pll)        !=dP(1,1,kbas,lbas) + dP(3,3,kbas,lbas)
                 !>
                 fgg = gg*pll
                 !>
                 dGij(1,1,i,j) = dGij(1,1,i,j) +        fgg
                 dGij(3,3,i,j) = dGij(3,3,i,j) + dconjg(fgg)
              !> end
            enddo
          enddo
        enddo
      enddo
      !> exchange part
      if( .not. doK )return
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          do j=1,jc
            do i=1,ic
              !> integrals processing: [ij|kl]
                 !> large-large block 
                 gg = dcmplx( gout(k,l,i,j), 0.0d0 )
              !> end
              !> contraction [ij|kl]*jk
                 paa_ll = dPj(1,1,j,kbas)
                 pba_ll = dPj(3,1,j,kbas)
                 pab_ll = dPj(1,3,j,kbas)
                 pbb_ll = dPj(3,3,j,kbas)
                 dGi(1,1,i,lbas) = dGi(1,1,i,lbas) - paa_ll*gg
                 dGi(1,3,i,lbas) = dGi(1,3,i,lbas) - pab_ll*gg
                 dGi(3,1,i,lbas) = dGi(3,1,i,lbas) - pba_ll*gg
                 dGi(3,3,i,lbas) = dGi(3,3,i,lbas) - pbb_ll*gg
              !> end
              !> contraction [ij|lk]*jl
                 paa_ll = dPj(1,1,j,lbas)
                 pba_ll = dPj(3,1,j,lbas)
                 pab_ll = dPj(1,3,j,lbas)
                 pbb_ll = dPj(3,3,j,lbas)
                 dGi(1,1,i,kbas) = dGi(1,1,i,kbas) - paa_ll*gg
                 dGi(1,3,i,kbas) = dGi(1,3,i,kbas) - pab_ll*gg
                 dGi(3,1,i,kbas) = dGi(3,1,i,kbas) - pba_ll*gg
                 dGi(3,3,i,kbas) = dGi(3,3,i,kbas) - pbb_ll*gg
              !> end
              !> contraction [ji|kl]*ik
                 paa_ll = dPi(1,1,i,kbas)
                 pba_ll = dPi(3,1,i,kbas)
                 pab_ll = dPi(1,3,i,kbas)
                 pbb_ll = dPi(3,3,i,kbas)
                 dGj(1,1,j,lbas) = dGj(1,1,j,lbas) - paa_ll*gg
                 dGj(1,3,j,lbas) = dGj(1,3,j,lbas) - pab_ll*gg
                 dGj(3,1,j,lbas) = dGj(3,1,j,lbas) - pba_ll*gg
                 dGj(3,3,j,lbas) = dGj(3,3,j,lbas) - pbb_ll*gg
              !> end
              !> contraction [ji|lk]*il
                 paa_ll = dPi(1,1,i,lbas)
                 pba_ll = dPi(3,1,i,lbas)
                 pab_ll = dPi(1,3,i,lbas)
                 pbb_ll = dPi(3,3,i,lbas)
                 dGj(1,1,j,kbas) = dGj(1,1,j,kbas) - paa_ll*gg
                 dGj(1,3,j,kbas) = dGj(1,3,j,kbas) - pab_ll*gg
                 dGj(3,1,j,kbas) = dGj(3,1,j,kbas) - pba_ll*gg
                 dGj(3,3,j,kbas) = dGj(3,3,j,kbas) - pbb_ll*gg
              !> end
            enddo
          enddo
        enddo
      enddo
    end subroutine

    subroutine max_shell_pair_contributions( dP, nrow, gmax, pmax, nshl )
      !> input
      integer,    intent(in) :: nrow
      integer,    intent(in) :: nshl 
      complex(8), intent(in) :: dP(4,4,nrow,nrow)
      !> output
      real(8), intent(out) :: gmax(  2,nshl,nshl)
      real(8), intent(out) :: pmax(2,2,nshl,nshl)
      !> local 
      real(8) :: gout(441*441*25)    !todo: rozumnejsia definicia
      integer :: i, ii, li, ni, oi
      integer :: j, jj, lj, nj, oj
      real(8) :: ei, ci, xi, yi, zi 
      real(8) :: ej, cj, xj, yj, zj 
      real(8) :: praa, piaa, prab, piab
      real(8) :: prba, piba, prbb, pibb
      !>
      !> find max{dsqrt(abs[ij|ij])}
      do j=1,nshl
        lj =       gto(j)%lvalue
        oj =       gto(j)%offset
        nj =       gto(j)%cdegen
        ej =       gto(j)%exponent(1)
        xj = atom( gto(j)%origin )%coordinate_x
        yj = atom( gto(j)%origin )%coordinate_y
        zj = atom( gto(j)%origin )%coordinate_z
        cj =       gto(j)%coefficient(1)
        do i=1,nshl
          li =       gto(i)%lvalue
          oi =       gto(i)%offset
          ni =       gto(i)%cdegen
          ei =       gto(i)%exponent(1)
          xi = atom( gto(i)%origin )%coordinate_x
          yi = atom( gto(i)%origin )%coordinate_y
          zi = atom( gto(i)%origin )%coordinate_z
          ci =       gto(i)%coefficient(1)
          !>
          call interest_eri_rkb('ssss',1.0d0,gout,&
                                li,ei,xi,yi,zi,ci,&
                                lj,ej,xj,yj,zj,cj,&
                                li,ei,xi,yi,zi,ci,&
                                lj,ej,xj,yj,zj,cj )
          !>
          call get_max_elements( ni, nj, gout, gmax(1,i,j), gmax(2,i,j) )
        enddo
      enddo
      !> find max{abs(dP(i,j))}
      pmax = 0.0d0
      do j=1,nshl
        oj = gto(j)%offset
        nj = gto(j)%cdegen
        do i=1,nshl
          oi = gto(i)%offset
          ni = gto(i)%cdegen
          do jj=1,nj
            do ii=1,ni
              !> large-large
              praa = dreal(dP(1,1,oi+ii,oj+jj))
              piaa = dimag(dP(1,1,oi+ii,oj+jj))
              prab = dreal(dP(1,3,oi+ii,oj+jj))
              piab = dimag(dP(1,3,oi+ii,oj+jj))
              prba = dreal(dP(3,1,oi+ii,oj+jj))
              piba = dimag(dP(3,1,oi+ii,oj+jj))
              prbb = dreal(dP(3,3,oi+ii,oj+jj))
              pibb = dimag(dP(3,3,oi+ii,oj+jj))
              pmax(1,1,i,j) = max(pmax(1,1,i,j),             &
                                  dsqrt(praa*praa+piaa*piaa),&
                                  dsqrt(prab*prab+piab*piab),&
                                  dsqrt(prba*prba+piba*piba),&
                                  dsqrt(prbb*prbb+pibb*pibb) )
              !> small-large
              praa = dreal(dP(2,1,oi+ii,oj+jj))
              piaa = dimag(dP(2,1,oi+ii,oj+jj))
              prab = dreal(dP(2,3,oi+ii,oj+jj))
              piab = dimag(dP(2,3,oi+ii,oj+jj))
              prba = dreal(dP(4,1,oi+ii,oj+jj))
              piba = dimag(dP(4,1,oi+ii,oj+jj))
              prbb = dreal(dP(4,3,oi+ii,oj+jj))
              pibb = dimag(dP(4,3,oi+ii,oj+jj))
              pmax(2,1,i,j) = max(pmax(2,1,i,j),             &
                                  dsqrt(praa*praa+piaa*piaa),&
                                  dsqrt(prab*prab+piab*piab),&
                                  dsqrt(prba*prba+piba*piba),&
                                  dsqrt(prbb*prbb+pibb*pibb) )
              !> large-small
              praa = dreal(dP(1,2,oi+ii,oj+jj))
              piaa = dimag(dP(1,2,oi+ii,oj+jj))
              prab = dreal(dP(1,4,oi+ii,oj+jj))
              piab = dimag(dP(1,4,oi+ii,oj+jj))
              prba = dreal(dP(3,2,oi+ii,oj+jj))
              piba = dimag(dP(3,2,oi+ii,oj+jj))
              prbb = dreal(dP(3,4,oi+ii,oj+jj))
              pibb = dimag(dP(3,4,oi+ii,oj+jj))
              pmax(1,2,i,j) = max(pmax(1,2,i,j),             &
                                  dsqrt(praa*praa+piaa*piaa),&
                                  dsqrt(prab*prab+piab*piab),&
                                  dsqrt(prba*prba+piba*piba),&
                                  dsqrt(prbb*prbb+pibb*pibb) )
              !> small-small
              praa = dreal(dP(2,2,oi+ii,oj+jj))
              piaa = dimag(dP(2,2,oi+ii,oj+jj))
              prab = dreal(dP(2,4,oi+ii,oj+jj))
              piab = dimag(dP(2,4,oi+ii,oj+jj))
              prba = dreal(dP(4,2,oi+ii,oj+jj))
              piba = dimag(dP(4,2,oi+ii,oj+jj))
              prbb = dreal(dP(4,4,oi+ii,oj+jj))
              pibb = dimag(dP(4,4,oi+ii,oj+jj))
              pmax(2,2,i,j) = max(pmax(2,2,i,j),             &
                                  dsqrt(praa*praa+piaa*piaa),&
                                  dsqrt(prab*prab+piab*piab),&
                                  dsqrt(prba*prba+piba*piba),&
                                  dsqrt(prbb*prbb+pibb*pibb) )
            enddo
          enddo
        enddo
      enddo
    end subroutine
    
    subroutine get_max_elements( ni, nj, gout, gll, gss ) 
      !> input
      integer, intent(in) :: ni, nj
      real(8), intent(in) :: gout(5,ni,nj,5,ni,nj) 
      !> output
      real(8), intent(out) :: gll, gss
      !> local
      integer :: i, j
      !>
      gll = 0.0d0 
      gss = 0.0d0 
      do j=1,nj
        do i=1,ni
          gll = max(gll,dsqrt(abs(gout(1,i,j,1,i,j))))
          gss = max(gss,dsqrt(abs(gout(2,i,j,2,i,j)+gout(3,i,j,3,i,j))),&
                        dsqrt(abs(gout(4,i,j,4,i,j)+gout(5,i,j,5,i,j))) )
        enddo
      enddo
    end subroutine

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
END MODULE
