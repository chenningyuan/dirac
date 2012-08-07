! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Module:       module_interest_interface.f90 
!
! Description:  - interfacing basic data between Dirac & InteRest program 
!               - providing two- and four-center integrals in RKB/RMB balanced basis
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

  implicit none

  public initialize_interest 
  public process_overlap_integrals 
  public process_two_electron_integrals
  public process_kinetic_energy_integrals 
  public process_nuclear_attraction_integrals

  interface process_two_electron_integrals
    module procedure process_two_electron_integrals_4comp
  end interface

  private

  !> matrix dimensions
     !> note: by definition nrow, ncol refer to the number of cartesian GTOs
     integer :: nshl = 1
     integer :: nrow = 0
     integer :: ncol = 0
     integer :: nr_shells = 0 
     real(8) :: prescreening_treshold = 1.d-15
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

CONTAINS
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE initialize_interest()

  use codata

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"
#include "dcbdhf.h"
#include "nucdata.h"

    !> local variables
       integer :: i
       integer :: j
       integer :: ij
       integer :: ier
       type( type_gto ) :: tmpgto
    !> end
  

    !> initialize InteRest integral package 
       call interest_initialize()
    !> end

    !> test if the large component basis is uncontracted
       if( nplrg.ne.nlarge )then
         write(6,*) ' Error: contracted basis set found => stop!'
         write(6,*) ' Note:  uncontracted basis sets are only supported'
         stop
       endif
    !> end
  
    !> interface constants
       constant = type_constant( dsqrt(alpha2), alpha2 )
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
    !> end

    !> reorder basis functions to the angular momentum 
    !> ascending order to increase the integral performance
    !> set the AO-vector offsets and shell-degeneracy
       do i=1,size(gto)-1
         do j=(i+1),size(gto)
           if( gto(j)%lvalue >= gto(i)%lvalue )cycle
           tmpgto = gto(i)
           gto(i) = gto(j)
           gto(j) = tmpgto
         enddo
       enddo
       do i=1,size(gto)-1
         gto(i+1)%offset = gto(i)%offset + gto(i)%cdegen
       enddo
    !> end

    !> temporary print ... 
       do i=1,size(gto)
         write(6,'(1x,5(a,i3,3x),3(a,d15.5,3x))') ' i = ',i,                             &
                                                  ' index  = ',     gto(i)%index,        &
                                                  ' offset = ',     gto(i)%offset,       &
                                                  ' lvalue = ',     gto(i)%lvalue,       &
                                                  ' origin = ',     gto(i)%origin,       &
                                                  ' exponent = ',   gto(i)%exponent(1),  &
                                                  ' coefficient = ',gto(i)%coefficient(1)
       enddo
      !write(6,*) 
      !write(6,'(2x,a,i5)') 'Basis set informations:'
      !write(6,'(2x,a,i5)') '-----------------------'
      !write(6,*) 
      !write(6,*) 
      !write(6,'(2x,a)')    'label    atom   charge   prim   cont   L-basis   RKB     S-basis'
      !write(6,'(2x,70a1)') ('-',i=1,70)
      !write(6,'(2x,70a1)') ('-',i=1,70)

       nrow      = sum( gto(:)%cdegen)
       ncol      = sum( gto(:)%cdegen)
       nr_shells = size(gto)

       write(6,*) 
       write(6,'(2x,a,i5  )') 'Total number of basis function shells:    ',size(gto)
       write(6,'(2x,a,i5  )') 'Total number of spherical basis functions:',sum(gto(:)%sdegen)
       write(6,'(2x,a,i5,a)') 'Total number of cartesian basis functions:',sum(gto(:)%cdegen),' (used for calculation)'
       write(6,*) 
    !> end

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_overlap_integrals( S ) 

    !> input
       real(8) :: S(2*nrow,2*ncol) 
    !> end

    !> local
       integer              :: i, j 
       real(8), allocatable :: Smat(:,:,:)
    !> end

    !> external InteRest integral routines
       external :: interest_overlap 
    !> end


    allocate( Smat( nrow,ncol,2 ) ); Smat = 0.0d0
    call two_center_integrals( interest_overlap, Smat, nrow )
    call symmetrize_matrix(                                 A = Smat(1,1,1), ida = nrow )  !Re(A)
    call symmetrize_matrix( factor = constant%alpha2/4.0d0, A = Smat(1,1,2), ida = nrow )  !Re(A) = (xx + yy + zz)
    do j=1,ncol
      do i=1,nrow
        S(i     ,j     ) = Smat(i,j,1) 
        S(i+nrow,j+ncol) = Smat(i,j,2) 
      enddo
    enddo
    deallocate( Smat )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_kinetic_energy_integrals( T ) 

    !> input
       real(8) :: T(2*nrow,2*ncol) 
    !> end

    !> local
       integer              :: i, j 
       real(8), allocatable :: Tmat(:,:,:)
    !> end

    !> external InteRest integral routines
       external :: interest_overlap 
    !> end


    allocate( Tmat( nrow,ncol,2 ) ); Tmat = 0.0d0
    call two_center_integrals( interest_overlap, Tmat, nrow )
    call symmetrize_matrix( factor = 0.5d0, A = Tmat(1,1,2), ida = nrow )  !Re(A) = (xx + yy + zz)
    do j=1,ncol
      do i=1,nrow
        T(i+nrow,j     ) = T(i+nrow,j     ) + Tmat(i,j,2)
        T(i     ,j+ncol) = T(i     ,j+ncol) + Tmat(i,j,2)
        T(i+nrow,j+ncol) = T(i+nrow,j+ncol) - Tmat(i,j,2)
      enddo
    enddo
    deallocate( Tmat )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE process_nuclear_attraction_integrals( V ) 
 
    !> input
    !> note: we assume nz=4 !!!
       real(8) :: V(2*nrow,2*ncol,4)
    !> end

    !> local
       integer              :: i, j 
       real(8), allocatable :: Vmat(:,:,:)
    !> end

    !> external InteRest integral routines
       external :: interest_nuclear_attraction_point_nucleus 
    !> end


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
    !> end

    !> output
       real(8) :: final_integrals(*)
    !> end

    !> local
       integer :: nrkb
       integer :: i, ni
       integer :: j, nj
       real(8) :: gout(5*441)  !todo: better definition 
    !> end


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
    !> end

    !> output
       real(8) :: B(idb,jdb,nrkb)
    !> end

    !> local
       integer :: i
       integer :: j
       integer :: k
    !> end


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
  SUBROUTINE process_two_electron_integrals_4comp( Gmat, Pmat, intflg ) 

    use module_interest_eri

    !> input
    !> note:  we assume nz=4 !!!
    !> fixme: the case when nz!=4
    integer, intent(in) :: intflg
    real(8), intent(in) :: Pmat(2*nrow,2*nrow,4)  ! 1-electron density matrix, Dirac format

    !> output
    real(8), intent(out) :: Gmat(2*nrow,2*nrow,4) ! 2-electron Fock matrix, Dirac format, from Pmat

    !> local
    logical :: doK
    real(8) :: fijkl
    real(8) :: scaleK 
    integer :: li, ni, oi   
    integer :: lj, nj, oj   
    integer :: lk, nk, ok   
    integer :: ll, nl, ol   
    integer :: nint, lmax
    real(8) :: fll, fsl, fss 
    integer :: i, j, k, l, ibas
    real(8) :: gout(441*441*25)    !todo: rozumnejsia definicia
    integer :: nij, nkl, ij
    real(8) :: ei, ci, xi, yi, zi 
    real(8) :: ej, cj, xj, yj, zj 
    real(8) :: ek, ck, xk, yk, zk 
    real(8) :: el, cl, xl, yl, zl 
 
    !> fixme: later through the input ...
    character(len=4), parameter :: eri_class = 'ssss'

    complex(8), allocatable :: dP(:,:,:,:)
    complex(8), allocatable :: dG(:,:,:,:)


    !> initialize for DFT
    scaleK = 1.0d0
    doK    = ( scaleK /= 0.0d0 )

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
                 CASE( 'ssll' )
                     call interest_eri_rkb('ssll',fijkl,gout,&
                                           li,ei,xi,yi,zi,ci,&
                                           lj,ej,xj,yj,zj,cj,&
                                           lk,ek,xk,yk,zk,ck,&
                                           ll,el,xl,yl,zl,cl )
                 CASE( 'ssss' ) 
                     call interest_eri_rkb('ssss',fijkl,gout,&
                                           li,ei,xi,yi,zi,ci,&
                                           lj,ej,xj,yj,zj,cj,&
                                           lk,ek,xk,yk,zk,ck,&
                                           ll,el,xl,yl,zl,cl )
                 CASE default
                      stop 'Error: this ERI class is not implemented yet => stop'
               END SELECT
               call process_dG_4comp_test_ssss_1x(ni,nj,nk,nl,              &
                                                  oi,oj,ok,ol,              &
                                                  nrow,doK,scaleK,gout,dP,dG)

             enddo lloop
           enddo kloop
         enddo jloop
       enddo iloop
    !> end


    !> data backward alignment 
    fsl = constant%alpha /2.0d0 
    fss = constant%alpha2/4.0d0 
    do j=1,nrow
      do i=1,nrow
        !> quaternion(0) 
        Gmat(i     ,j     ,1) = Gmat(i     ,j     ,1) + dreal(dG(1,1,i,j))
        Gmat(i+nrow,j     ,1) = Gmat(i+nrow,j     ,1) + dreal(dG(2,1,i,j))*fsl
        Gmat(i     ,j+nrow,1) = Gmat(i     ,j+nrow,1) + dreal(dG(1,2,i,j))*fsl
        Gmat(i+nrow,j+nrow,1) = Gmat(i+nrow,j+nrow,1) + dreal(dG(2,2,i,j))*fss
        !> quaternion(i) 
        Gmat(i     ,j     ,2) = Gmat(i     ,j     ,2) + dimag(dG(1,1,i,j))
        Gmat(i+nrow,j     ,2) = Gmat(i+nrow,j     ,2) + dimag(dG(2,1,i,j))*fsl
        Gmat(i     ,j+nrow,2) = Gmat(i     ,j+nrow,2) + dimag(dG(1,2,i,j))*fsl
        Gmat(i+nrow,j+nrow,2) = Gmat(i+nrow,j+nrow,2) + dimag(dG(2,2,i,j))*fss
        !> quaternion(j) 
        Gmat(i     ,j     ,3) = Gmat(i     ,j     ,3) + dreal(dG(1,3,i,j))
        Gmat(i+nrow,j     ,3) = Gmat(i+nrow,j     ,3) + dreal(dG(2,3,i,j))*fsl
        Gmat(i     ,j+nrow,3) = Gmat(i     ,j+nrow,3) + dreal(dG(1,4,i,j))*fsl
        Gmat(i+nrow,j+nrow,3) = Gmat(i+nrow,j+nrow,3) + dreal(dG(2,4,i,j))*fss
        !> quaternion(k) 
        Gmat(i     ,j     ,4) = Gmat(i     ,j     ,4) + dimag(dG(1,3,i,j))
        Gmat(i+nrow,j     ,4) = Gmat(i+nrow,j     ,4) + dimag(dG(2,3,i,j))*fsl
        Gmat(i     ,j+nrow,4) = Gmat(i     ,j+nrow,4) + dimag(dG(1,4,i,j))*fsl
        Gmat(i+nrow,j+nrow,4) = Gmat(i+nrow,j+nrow,4) + dimag(dG(2,4,i,j))*fss
      enddo
    enddo

    deallocate( dP )
    deallocate( dG )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_dG_4comp_test_ssss_1x( ic, jc, kc, lc, io, jo, ko, lo, nbas, doK, scaleK, gout, pmat, gmat ) 

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
            !> contraction [ij|kl]*lk
               paa = pmat(2,2,lbas,kbas)
               pba = pmat(4,2,lbas,kbas)
               pab = pmat(2,4,lbas,kbas)
               pbb = pmat(4,4,lbas,kbas)
               pll = pmat(1,1,lbas,kbas) + pmat(3,3,lbas,kbas)
               !> large-large block
               gmat(1,1,ibas,jbas) = gmat(1,1,ibas,jbas) + gg*pll + ggaa*paa + ggab*pba + ggba*pab + ggbb*pbb
               gmat(3,3,ibas,jbas) = gmat(3,3,ibas,jbas) + gg*pll + ggaa*paa + ggab*pba + ggba*pab + ggbb*pbb
               !> small-small block
               gmat(2,2,ibas,jbas) = gmat(2,2,ibas,jbas) + gaag*pll + gaagaa*paa + gaagab*pba + gaagba*pab + gaagbb*pbb
               gmat(4,2,ibas,jbas) = gmat(4,2,ibas,jbas) + gbag*pll + gbagaa*paa + gbagab*pba + gbagba*pab + gbagbb*pbb
               gmat(2,4,ibas,jbas) = gmat(2,4,ibas,jbas) + gabg*pll + gabgaa*paa + gabgab*pba + gabgba*pab + gabgbb*pbb
               gmat(4,4,ibas,jbas) = gmat(4,4,ibas,jbas) + gbbg*pll + gbbgaa*paa + gbbgab*pba + gbbgba*pab + gbbgbb*pbb
            !> end
          enddo
        enddo
      enddo
    enddo

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

    !> exchange part
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

    !fixme
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

    !fixme
    !call backward_ao_resorting( A )
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>on input:  the matrix in "input"  AO order 
!>on output: the matrix in integral AO order 
!>fixme:     clsed-shell and LL-only case!!!
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
            A(io_new+ii,jo_new+jj) = B(io_old+ii,jo_old+jj)  !Re(A,LL) 
          enddo
        enddo
      enddo
    enddo
    deallocate( B )
        
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>on input:  the matrix in "integral" AO order 
!>on output: the matrix in "input"    AO order 
!>fixme:     clsed-shell and LL-only case!!!
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
            A(io_new+ii,jo_new+jj) = B(io_old+ii,jo_old+jj)  !Re(A,LL) 
          enddo
        enddo
      enddo
    enddo
    deallocate( B )
        
  END SUBROUTINE
! ------------------------------------------------------------------------------------
END MODULE
