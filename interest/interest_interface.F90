! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Module:       interest_interface 
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

MODULE INTEREST_INTERFACE

  implicit none

  public initialize_interest 
  public two_electron_integrals
  public process_overlap_integrals 
  public process_nuclear_attraction_integrals

  private

  !> definition of "constant" data type 
     type type_constant
          real(8) :: alpha2      !=(1/c)^2 
     end  type
     type(type_constant) :: constant 
  !> end

  !> definition of "symmetry" data type 
     !fixme: not initialized yet ...
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
          real(8) :: norm(1)
          real(8) :: exponent(1)
          real(8) :: coefficient(1)
     end  type
     type(type_gto), allocatable, target :: gto(:)
  !> end

  !> fixme: temporary
     !> definition of "real_vector" data type
        type real_vector
             real(8), pointer :: elms(:)
        end  type
     !> end

     !> definition of "real_matrix" data type 
        type real_matrix
             integer          :: nrow 
             integer          :: ncol  
             real(8), pointer :: elms(:,:) 
        end  type
     !> end
     
     !> definition of "complex_pointer" data type 
        type complex_pointer
             complex(8), pointer :: p => NULL()
        end  type
     !> end
     
     !> definition of "complex_matrix" data type 
        type complex_matrix
             integer                 :: nrow 
             integer                 :: ncol  
             type(complex_pointer)   :: info(4,4)  
             complex(8), allocatable :: elms(:,:,:) 
        end  type
     !> end
  !> end

 
  interface two_electron_integrals 
    module procedure two_electron_integrals_1comp
   !module procedure two_electron_integrals_4comp
  end interface

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
       constant = type_constant( alpha2 )
    !> end

    !> interface molecular data
       allocate( atom(nucind), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: atom(:)'
       do i=1,size(atom)
         atom(i) = type_atom( charge(i), cord(1,i), cord(2,i), cord(3,i), gnuexp(i) )
       enddo
    !> end
  
    !> interface basis set data
       allocate( gto(nlrgsh), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: gto(:)'
       i = 1
       do while( i <= size(gto) )  
         do j=1,nrco(i)
           gto(i)  = type_gto( i, 0, nhkt(i), ncent(i), 0, 0, priccf(i,j), priexp(i), priccf(i,j) )
           i = i + 1
         enddo
       enddo
    !> end

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
    !> end

    !> set some gto data + AO-vector offsets 
       j=0
       do i=1,size(gto)
         gto(i)%offset = j 
         gto(i)%sdegen = 2*gto(i)%lvalue - 1 
         gto(i)%cdegen = ( gto(i)%lvalue*(gto(i)%lvalue+1) )/2 
         j=j+gto(i)%cdegen !fixme:this is offset for cartesian AOs only!!!
       enddo
    !> end

    !> temporary print ... 
       do i=1,size(gto)
         write(6,'(1x,5(a,i3,3x),3(a,d15.5,3x))') ' i = ',i,                             &
                                                  ' index  = ',     gto(i)%index,        &
                                                  ' offset = ',     gto(i)%offset,       &
                                                  ' lvalue = ',     gto(i)%lvalue,       &
                                                  ' origin = ',     gto(i)%origin,       &
                                                  ' norm = ',       gto(i)%norm(1),      &
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
       real(8) :: S(*) 
    !> end

    !> local
       integer              :: i 
       integer              :: nrow
       real(8), allocatable :: Smat(:,:)
    !> end

    !> external InteRest integral routines
       external :: interest_overlap 
    !> end


    nrow = sum(gto(:)%cdegen) 
    allocate( Smat( nrow*nrow, 2 ) ); Smat = 0.0d0

    !> calculate overlap-Smat(:,1) and kinetic-Smat(:,2) integrals
    !> note: Smat(:,2) is not devided by the factor 2.0d0
       call two_center_integrals( interest_overlap, Smat, nrow )
       call symmetrize_matrix(                                 A = Smat(1,1), ida = nrow )
       call symmetrize_matrix( factor = 4.0d0*constant%alpha2, A = Smat(1,2), ida = nrow )
    !> end

    S(            1:  (nrow*nrow)) = Smat(1:nrow*nrow,1)
    S((nrow*nrow)+1:2*(nrow*nrow)) = Smat(1:nrow*nrow,2)

    deallocate( Smat )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE process_nuclear_attraction_integrals( V ) 
 
    !> input
       real(8) :: V(*) 
    !> end

    !> local
       integer              :: i 
       integer              :: nrow
       real(8), allocatable :: Vmat(:,:)
    !> end

    !> external InteRest integral routines
       external :: interest_nuclear_attraction_point_nucleus 
    !> end


    nrow = sum(gto(:)%cdegen) 
    allocate( Vmat( nrow*nrow, 5 ) ); Vmat = 0.0d0

    !> calculate nuclear attraction integrals in RKB basis
       allocate( Vmat( nrow*nrow, 5 ) ); Vmat = 0.0d0
       call two_center_integrals( interest_nuclear_attraction_point_nucleus, Vmat, nrow )
       call symmetrize_matrix    (                           A = Vmat(1,1), ida = nrow )   !Re(A )
       call symmetrize_matrix    ( factor = constant%alpha2, A = Vmat(1,2), ida = nrow )   !Re(A') = (xx + yy + zz)*alpha^{2}
       call antisymmetrize_matrix( factor = constant%alpha2, A = Vmat(1,3), ida = nrow )   !Im(A') = (xy - yx)     *alpha^{2} 
       call antisymmetrize_matrix( factor = constant%alpha2, A = Vmat(1,4), ida = nrow )   !Re(B') = (zx - xz)     *alpha^{2} 
       call antisymmetrize_matrix( factor = constant%alpha2, A = Vmat(1,5), ida = nrow )   !Im(B') = (yz - zy)     *alpha^{2}
    !> end

    V(              1:  (nrow*nrow)) = Vmat(1:nrow*nrow,1)
    V(  (nrow*nrow)+1:2*(nrow*nrow)) = Vmat(1:nrow*nrow,2)
    V(2*(nrow*nrow)+1:3*(nrow*nrow)) = Vmat(1:nrow*nrow,3)
    V(3*(nrow*nrow)+1:4*(nrow*nrow)) = Vmat(1:nrow*nrow,4)
    V(4*(nrow*nrow)+1:5*(nrow*nrow)) = Vmat(1:nrow*nrow,5)

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

    integer, intent(in)  :: ida           !fixme:umelo dodane
    logical              :: rkb = .true.  !fixme:umelo dodane

    !> output
    real(8) :: final_integrals(*)

    !> local
    integer :: nrkb
    integer :: i, ni
    integer :: j, nj
    real(8) :: gout(4410)  !todo: rozumnejsia definicia


    do j=1,size(gto)
      do i=j,size(gto)
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

        call store_batch(gout, ni, nj,  nrkb,              &
                         final_integrals, ida, ida,        &
                         gto(i)%offset,                    &
                         gto(j)%offset                     )

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
    real(8) :: A(ida,jda,nrkb)

    !> output
    real(8) :: B(idb,jdb,nrkb)

    !> local
    integer :: i
    integer :: j
    integer :: k


    do k=1,nrkb
      do j=1,jda
        do i=1,ida
          B(ioffset+i,joffset+j,k) = A(i,j,k) 
        enddo
      enddo
    enddo
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE two_electron_integrals_1comp( interest_routine, G, TrGP ) 

    !> input
       interface 
         subroutine interest_routine(rkb,                                            &
                                     factor,fint,nint,                               &
                                     la,alpha,ax,ay,az,anorm,lb,beta ,bx,by,bz,bnorm,&
                                     lc,gamma,cx,cy,cz,cnorm,ld,delta,dx,dy,dz,dnorm )
       
           integer, intent(out) :: nint
           real(8), intent(out) :: fint(*)
           logical, intent(in ) :: rkb 
           real(8), intent(in ) :: factor 
           integer, intent(in ) :: la,lb,lc,ld
           real(8), intent(in ) :: alpha,ax,ay,az,anorm
           real(8), intent(in ) :: beta, bx,by,bz,bnorm
           real(8), intent(in ) :: gamma,cx,cy,cz,cnorm
           real(8), intent(in ) :: delta,dx,dy,dz,dnorm
         end subroutine
       end interface
    !> end

    !> output
       type(real_matrix)    :: G
       real(8), intent(out) :: TrGP
    !> end

    !> local
       integer :: nint
       integer :: ni, nj, nk, nl
       integer :: i, j, k, l, ij, kl, lmax 
       integer :: nij,  nkl,  nik,  njk,  nil,  njl
       integer :: nij2, nkl2, nik2, njk2, nil2, njl2
       real(8) :: gout(441)                           !todo: rozumnejsia definicia
       real(8) :: fij, fkl, fijkl
       real(8), allocatable :: dpij(:) 
       real(8), allocatable :: dpkl(:) 
       real(8), allocatable :: dplj(:) 
       real(8), allocatable :: dpkj(:) 
       real(8), allocatable :: dpli(:) 
       real(8), allocatable :: dpki(:) 
       real(8), allocatable :: dgij(:) 
       real(8), allocatable :: dgkl(:) 
       real(8), allocatable :: dgjl(:) 
       real(8), allocatable :: dgjk(:) 
       real(8), allocatable :: dgil(:) 
       real(8), allocatable :: dgik(:) 
       type(real_vector), allocatable :: dP(:,:,:) 
    !> end
    
   

    !> allocate the working temporary fields
       allocate( dpij(2*441) )                   !todo: zmenit do buducna
       allocate( dpkl(2*441) )                   !todo: zmenit do buducna
       allocate( dplj(2*441) )                   !todo: zmenit do buducna
       allocate( dpkj(2*441) )                   !todo: zmenit do buducna
       allocate( dpli(2*441) )                   !todo: zmenit do buducna
       allocate( dpki(2*441) )                   !todo: zmenit do buducna
       allocate( dgij(2*441) )                   !todo: zmenit do buducna
       allocate( dgkl(2*441) )                   !todo: zmenit do buducna
       allocate( dgjl(2*441) )                   !todo: zmenit do buducna
       allocate( dgjk(2*441) )                   !todo: zmenit do buducna
       allocate( dgil(2*441) )                   !todo: zmenit do buducna
       allocate( dgik(2*441) )                   !todo: zmenit do buducna
       allocate( dP( 2, size(gto), size(gto) ) )
    !> end


    !> create shell-blocks of AO(-difference) density
       !> allocation 
          call doP_1comp( 0, G%elms(:,1), G%nrow, dP, size(gto) )
       !> alpha
          call doP_1comp( 1, G%elms(:,1), G%nrow, dP, size(gto) )
       !> beta
          call doP_1comp( 2, G%elms(:,2), G%nrow, dP, size(gto) )
    !> end


    !> loop over shells
iloop: do i=1,size(gto)
         ni = gto(i)%cdegen
         ij = (i-1)*size(gto)

jloop:   do j=1,i
           ij  = ij+1 
           fij = 1.0d0 
           nj  = gto(j)%cdegen
           if( j == i ) fij=0.5d0

           !> J-preparation
              nij  = ni*nj
              nij2 = 2*nij
              dgij(1:nij2)=0.0d0
             !dpij = dP(1,i,j)%elms  !komentovane pre trojuholnik
           !> end

kloop:     do k=1,size(gto)
             nk = gto(k)%cdegen
             kl = (k-1)*size(gto)

             !> K-preparation
                nik  = ni*nk
                njk  = nj*nk
                nik2 = 2*nik
                njk2 = 2*njk
                dgik(1:nik2) = 0.0d0
                dgjk(1:njk2) = 0.0d0
                dpkj = dP(1,k,j)%elms
                dpki = dP(1,k,i)%elms
             !> end

lloop:       do l=1,k
               kl  = kl+1
               fkl = 1.0d0 
               nl  = gto(l)%cdegen
               if(  k ==  l ) fkl=0.5d0
               !if( ij == kl )fijkl=0.5d0*fijkl  zakomentovany prefaktor lebo nerobim ij<->kl

               !> J-preparation
                  nkl  = nk*nl
                  nkl2 = 2*nkl
                 !dgkl(1:nkl2) = 0.0d0 !komentovane pre trojuholnik
                  dpkl = dP(1,k,l)%elms
               !> end

               !> K-preparation
                  nil  = ni*nl
                  njl  = nj*nl
                  nil2 = 2*nil
                  njl2 = 2*njl
                  dgil(1:nil2) = 0.0d0
                  dgjl(1:njl2) = 0.0d0
                  dplj = dP(1,l,j)%elms
                  dpli = dP(1,l,i)%elms
               !> end


               !> call InteRest library routine for a given batch
                  call interest_routine(.false.,                            &
                                        1.0d0,gout,nint,                    &
                                        !> bra-function(1)
                                        gto(i)%lvalue,                      &
                                        gto(i)%exponent(1),                 &
                                        atom( gto(i)%origin )%coordinate_x, &
                                        atom( gto(i)%origin )%coordinate_y, &
                                        atom( gto(i)%origin )%coordinate_z, &
                                        gto(i)%coefficient(1),              &
                                        !> ket-function(1)
                                        gto(j)%lvalue,                      &
                                        gto(j)%exponent(1),                 &
                                        atom( gto(j)%origin )%coordinate_x, &
                                        atom( gto(j)%origin )%coordinate_y, &
                                        atom( gto(j)%origin )%coordinate_z, &
                                        gto(j)%coefficient(1),              &
                                        !> bra-function(2)
                                        gto(k)%lvalue,                      &
                                        gto(k)%exponent(1),                 &
                                        atom( gto(k)%origin )%coordinate_x, &
                                        atom( gto(k)%origin )%coordinate_y, &
                                        atom( gto(k)%origin )%coordinate_z, &
                                        gto(k)%coefficient(1),              &
                                        !> ket-function(2)
                                        gto(l)%lvalue,                      &
                                        gto(l)%exponent(1),                 &
                                        atom( gto(l)%origin )%coordinate_x, &
                                        atom( gto(l)%origin )%coordinate_y, &
                                        atom( gto(l)%origin )%coordinate_z, &
                                        gto(l)%coefficient(1)               )
               !> end

               !> process integral batch with the density matrix to get JK terms 
                  call process_1comp_JK( gout, nij, nkl,       dpij, dgij, &
                                                               dpkl, dgkl, &
                                                               fkl,        &
                                               ni, nj, nk, nl, dplj, dgjl, &
                                                               dpkj, dgjk, &   
                                                               dpli, dgil, &
                                                               dpki, dgik, &
                                                               fij*fkl     )
               !> end 

               !> add K-contribution
                  dP(2,j,l)%elms(1:njl2) = dP(2,j,l)%elms(1:njl2) - dgjl(1:njl2)
                  dP(2,i,l)%elms(1:nil2) = dP(2,i,l)%elms(1:nil2) - dgil(1:nil2)
               !> end
             enddo lloop

             !> add K-contribution
                dP(2,j,k)%elms(1:njk2) = dP(2,j,k)%elms(1:njk2) - dgjk(1:njk2)
                dP(2,i,k)%elms(1:nik2) = dP(2,i,k)%elms(1:nik2) - dgik(1:nik2)
             !> end
           enddo kloop   

           !> add J-contribution
              dP(2,i,j)%elms(1:nij2) = dP(2,i,j)%elms(1:nij2) + 2.0d0*dgij(1:nij2)
           !> end 
         enddo jloop

       enddo iloop   
    !> end


    !> add G-contribution to the Fock matrix and calculate Tr(GP) 
       TrGP = 0.0d0
       !> alpha
          call doG_1comp( 1, dP, size(gto), G%elms(:,1), G%nrow, TrGP )
          call symmetrize_matrix( A = G%elms(:,1), ida = G%nrow       )
       !> beta
          call doG_1comp( 2, dP, size(gto), G%elms(:,2), G%nrow, TrGP )
          call symmetrize_matrix( A = G%elms(:,2), ida = G%nrow       )
       !> deallocation 
          call doG_1comp( 0, dP, size(gto), G%elms(:,2), G%nrow, TrGP )
    !> end
    

    deallocate( dpij )
    deallocate( dpkl )
    deallocate( dplj )
    deallocate( dpkj )
    deallocate( dpli )
    deallocate( dpki )
    deallocate( dgij )
    deallocate( dgkl )
    deallocate( dgjl )
    deallocate( dgjk )
    deallocate( dgil )
    deallocate( dgik )
    deallocate( dP   )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>note: doP_1comp a doP_4comp sa da pekne unifikovat, pokial budem mat optional GR a GI pre 1comp a 4comp data structure 
!>note: a pokial dlzku dat (v tomto pripade 2) nastavim ako vstupny parameter ktory mozem lahko menit z calling routine 
  SUBROUTINE doP_1comp( index, GR, idg, P, idp )

    !> input
       integer, intent(in) :: idg 
       integer, intent(in) :: idp 
       integer, intent(in) :: index 
       real(8), intent(in) :: GR(idg,idg)
    !> end

    !> output
       type(real_vector) :: P(2,idp,idp)
    !> end

    !> local
       integer :: ij
       integer :: i, ii, io, icc
       integer :: j, jj, jo, jcc
    !> end


    if( index == 0 )then
      do i=1,size(gto)
        icc = gto(i)%cdegen
        do j=1,size(gto)                           !fixme: reduce to lower/upper triangular form (i.le.j)
          jcc = gto(j)%cdegen
          allocate( P(1,i,j)%elms(2*icc*jcc) )     !fixme: always open-shell case 
          allocate( P(2,i,j)%elms(2*icc*jcc) )     !fixme: always open-shell case 
        enddo
      enddo 
      return
    endif

    do i=1,size(gto)
      io  = gto(i)%offset
      icc = gto(i)%cdegen
      do j=1,size(gto)                             !fixme: reduce to lower/upper triangular form (i.le.j)
        jo  = gto(j)%offset
        jcc = gto(j)%cdegen
        ij  =  index
        !> density-block
        do jj=1,jcc
          do ii=1,icc
            P(1,i,j)%elms(ij) = GR( jo+jj, io+ii ) !note:  P(ba) -> sa uklada ako (ab)
            ij = ij + 2                            !fixme: tu je 2 pre open-shell case 
          enddo
        enddo
        !> fock-block 
        ij =  index
        do jj=1,jcc
          do ii=1,icc
            P(2,i,j)%elms(ij) = 0.0d0
            ij = ij + 2                            !fixme: tu je 2 pre open-shell case 
          enddo
        enddo
      enddo
    enddo
 
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> Process Coulomb and Exchange contributions to Fock matrix 
  SUBROUTINE process_1comp_JK( gout, nij, nkl, pij, gij, pkl, gkl, fkl, &
                                               ni,  nj,  nk,  nl,       &
                           plj, gjl, pkj, gjk, pli, gil, pki, gik, fijkl)

    !> input 
       integer, intent(in) :: nij, nkl
       real(8), intent(in) :: fkl, fijkl 
       integer, intent(in) :: ni, nj, nk, nl 
       real(8), intent(in) :: gout(nkl,nij) 
    !> end

    !> input/output
       real(8), intent(inout) :: pij(2,nij), gij(2,nij)
       real(8), intent(inout) :: pkl(2,nkl), gkl(2,nkl)

       real(8), intent(inout) :: plj(2,nl,nj), gik(2,ni,nk)  !input P(l,j) -> output G(ik)
       real(8), intent(inout) :: pkj(2,nk,nj), gil(2,ni,nl)  !input P(k,j) -> output G(il)
       real(8), intent(inout) :: pli(2,nl,ni), gjk(2,nj,nk)  !input P(l,i) -> output G(jk)
       real(8), intent(inout) :: pki(2,nk,ni), gjl(2,nj,nl)  !input P(k,i) -> output G(jl)
    !> end

    !> local
       integer :: ij, kl
       integer :: i, j, k, l
       real(8) :: tk, pij1, pij2
    !> end


    !> Coulomb contribution
       !> contribution to G(ij)
          do ij=1,nij
            tk = 0.0d0
            do kl=1,nkl
              tk = tk + gout(kl,ij) * pkl(1,kl) !P_alpha
              tk = tk + gout(kl,ij) * pkl(2,kl) !P_beta
            enddo
            gij(1,ij) = gij(1,ij) + tk * fkl    !G_alpha
            gij(2,ij) = gij(2,ij) + tk * fkl    !G_beta
          enddo
       !> end
       
!      !> contribution to G(kl)
!         do ij=1,nij
!           pij1 = pij(1,ij,1) !P_alpha
!           pij2 = pij(2,ij,1) !P_beta
!           do kl=1,nkl
!             pkl(1,kl,2) = pkl(1,kl,2) + gout(kl,ij)*pij1 !G_alpha
!             pkl(2,kl,2) = pkl(2,kl,2) + gout(kl,ij)*pij2 !G_beta
!           enddo
!         enddo
!      !> end
    !> end


    !> Exchange contribution
       ij = 0
       do j=1,nj
         do i=1,ni
           kl = 0
           ij = ij + 1
           do l=1,nl
             do k=1,nk
               kl = kl + 1
               gil( 1,i,l ) = gil( 1,i,l ) + gout(kl,ij) * pkj( 1,k,j ) * fijkl
               gil( 2,i,l ) = gil( 2,i,l ) + gout(kl,ij) * pkj( 2,k,j ) * fijkl
               gjl( 1,j,l ) = gjl( 1,j,l ) + gout(kl,ij) * pki( 1,k,i ) * fijkl
               gjl( 2,j,l ) = gjl( 2,j,l ) + gout(kl,ij) * pki( 2,k,i ) * fijkl
               gik( 1,i,k ) = gik( 1,i,k ) + gout(kl,ij) * plj( 1,l,j ) * fijkl  
               gik( 2,i,k ) = gik( 2,i,k ) + gout(kl,ij) * plj( 2,l,j ) * fijkl
               gjk( 1,j,k ) = gjk( 1,j,k ) + gout(kl,ij) * pli( 1,l,i ) * fijkl
               gjk( 2,j,k ) = gjk( 2,j,k ) + gout(kl,ij) * pli( 2,l,i ) * fijkl
             enddo
           enddo
         enddo
       enddo
    !> end

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>tuto routinu je mozne generalizovat podla popisu postupu v doP routine ...
  SUBROUTINE doG_1comp( index, P, idp, GR, idg, TrGP )

    !> input
       integer, intent(in) :: idg 
       integer, intent(in) :: idp 
       integer, intent(in) :: index 
       type(real_vector) :: P(2,idp,idp)
    !> end

    !> output
       real(8), intent(  out) :: TrGP 
       real(8), intent(inout) :: GR(idg,idg)
    !> end

    !> local
       integer :: ij
       integer :: i, ii, io, icc
       integer :: j, jj, jo, jcc
    !> end


    if( index == 0 )then
      do j=1,size(gto)                    !fixme: reduce to lower/upper triangular form (i.le.j)
        do i=1,size(gto)
          deallocate( P(1,i,j)%elms )     !fixme: always open-shell case 
          deallocate( P(2,i,j)%elms )     !fixme: always open-shell case 
        enddo
      enddo 
      TrGP = 2.0d0*TrGP
      return
    endif


    do j=1,size(gto)
      jo  = gto(j)%offset
      jcc = gto(j)%cdegen
      !> energy / fock-block 
         ij =  index
         do jj=1,jcc
           do ii=1,jcc
             TrGP = TrGP + P(1,j,j)%elms(ij)*P(2,j,j)%elms(ij)*0.5d0
             GR( jo+ii, jo+jj ) = P(2,j,j)%elms(ij)
             ij = ij + 2                              !fixme: tu je 2 pre open-shell case 
           enddo
         enddo
      !> end
      do i=(j+1),size(gto),1
        io  = gto(i)%offset
        icc = gto(i)%cdegen
        !> energy / fock-block 
           ij =  index
           do jj=1,jcc
             do ii=1,icc
               TrGP = TrGP + P(1,i,j)%elms(ij)*P(2,i,j)%elms(ij)
               GR( io+ii, jo+jj ) = P(2,i,j)%elms(ij)
               ij = ij + 2                            !fixme: tu je 2 pre open-shell case 
             enddo
           enddo
        !> end
      enddo
    enddo
 
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE two_electron_integrals_4comp( interest_routine, G ) 

    !> input
    interface 
      subroutine interest_routine(rkb,                                            &
                                  fint,nint,                                      &
                                  la,alpha,ax,ay,az,anorm,lb,beta ,bx,by,bz,bnorm,&
                                  lc,gamma,cx,cy,cz,cnorm,ld,delta,dx,dy,dz,dnorm )

        integer, intent(out) :: nint
        real(8), intent(out) :: fint(*)
        logical, intent(in ) :: rkb 
        integer, intent(in ) :: la,lb,lc,ld
        real(8), intent(in ) :: alpha,ax,ay,az,anorm
        real(8), intent(in ) :: beta, bx,by,bz,bnorm
        real(8), intent(in ) :: gamma,cx,cy,cz,cnorm
        real(8), intent(in ) :: delta,dx,dy,dz,dnorm
      end subroutine
    end interface

    logical :: rkb = .true. !todo: zmenit do buducna, umelo dodane

    !> output
    type(complex_matrix) :: G

    !> local
    integer :: nint
    integer :: ni, nj, nk, nl
    integer :: i, j, k, l, lmax 
    real(8) :: gout(44100)                    !todo: rozumnejsia definicia
    real(8),           allocatable :: dpij(:) 
    type(real_vector), allocatable :: dp(:,:) 
    
   
    allocate( dpij(4410*2)               )    !todo: zmenit do buducna
    allocate( dp( size(gto), size(gto) ) )

    !> create shell-blocks of difference density (dp)
       !> Re(aa,ll)  
          call doP( index= 1, GR=G%elms(1,1,1), idg=G%nrow, P=dp, idp=size(gto) )
       !> Re(aa,ss)  
          call doP( index= 2, GR=G%elms(1,2,2), idg=G%nrow, P=dp, idp=size(gto) )
       !> Im(aa,ss)  
          call doP( index= 3, GI=G%elms(1,2,2), idg=G%nrow, P=dp, idp=size(gto) )
       !> Re(ab,ss)  
          call doP( index= 4, GR=G%elms(1,2,4), idg=G%nrow, P=dp, idp=size(gto) )
       !> Im(ab,ss)  
          call doP( index= 5, GI=G%elms(1,2,4), idg=G%nrow, P=dp, idp=size(gto) )
       !> Re(bb,ll)  
          call doP( index= 6, GR=G%elms(1,3,3), idg=G%nrow, P=dp, idp=size(gto) )
       !> Re(bb,ss)  
          call doP( index= 7, GR=G%elms(1,4,4), idg=G%nrow, P=dp, idp=size(gto) )
       !> Im(bb,ss)  
          call doP( index= 8, GI=G%elms(1,4,4), idg=G%nrow, P=dp, idp=size(gto) )
       !> Re(ba,ss)  
          call doP( index= 9, GR=G%elms(1,4,2), idg=G%nrow, P=dp, idp=size(gto) )
       !> Im(ba,ss)  
          call doP( index=10, GI=G%elms(1,4,2), idg=G%nrow, P=dp, idp=size(gto) )
    !> end


    !> loop over shells
iloop: do i=1,size(gto)
         ni = gto(i)%cdegen
jloop:   do j=1,i
           nj   = gto(j)%cdegen
           dpij = dp(i,j)%elms
kloop:     do k=1,size(gto)
             nk   = gto(k)%cdegen
lloop:       do l=1,k
               nl = gto(l)%cdegen


               !> call InteRest library routine for a given batch
                  call interest_routine(rkb,                                &
                                        gout, nint,                         &
                                        !> bra-function(1)
                                        gto(i)%lvalue,                      &
                                        gto(i)%exponent(1),                 &
                                        atom( gto(i)%origin )%coordinate_x, &
                                        atom( gto(i)%origin )%coordinate_y, &
                                        atom( gto(i)%origin )%coordinate_z, &
                                        gto(i)%coefficient(1),              &
                                        !> ket-function(1)
                                        gto(j)%lvalue,                      &
                                        gto(j)%exponent(1),                 &
                                        atom( gto(j)%origin )%coordinate_x, &
                                        atom( gto(j)%origin )%coordinate_y, &
                                        atom( gto(j)%origin )%coordinate_z, &
                                        gto(j)%coefficient(1),              &
                                        !> bra-function(2)
                                        gto(k)%lvalue,                      &
                                        gto(k)%exponent(1),                 &
                                        atom( gto(k)%origin )%coordinate_x, &
                                        atom( gto(k)%origin )%coordinate_y, &
                                        atom( gto(k)%origin )%coordinate_z, &
                                        gto(k)%coefficient(1),              &
                                        !> ket-function(2)
                                        gto(l)%lvalue,                      &
                                        gto(l)%exponent(1),                 &
                                        atom( gto(l)%origin )%coordinate_x, &
                                        atom( gto(l)%origin )%coordinate_y, &
                                        atom( gto(l)%origin )%coordinate_z, &
                                        gto(l)%coefficient(1)               )

               !> process JK-contribution
                  call process_JK( gout, ni*nj, nk*nl, dpij, dp(k,l)%elms, &
                                         ni,    nj,    nk,   nl,           &
                                         dp(i,l)%elms,                     &   
                                         dp(k,j)%elms,                     &
                                         dp(j,l)%elms,                     &
                                         dp(k,i)%elms                      )
               !> end 
             enddo lloop
           enddo kloop   
           !> add J-contribution
              call add_J( ni*nj, dpij, dp(i,j)%elms )
           !> end 
         enddo jloop
       enddo iloop   
       !> process G-contribution to the Fock matrix 
          call doG( 'aall', dp, size(gto), G%elms(1,1,1), G%nrow )
          call doG( 'aass', dp, size(gto), G%elms(1,2,2), G%nrow )
          call doG( 'abss', dp, size(gto), G%elms(1,2,4), G%nrow )
          call doG( 'bbll', dp, size(gto), G%elms(1,3,3), G%nrow )
          call doG( 'bbss', dp, size(gto), G%elms(1,4,4), G%nrow )
          call doG( 'bass', dp, size(gto), G%elms(1,4,2), G%nrow )
       !> end
    !> end

    deallocate( dp   )
    deallocate( dpij )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> Process Coulomb and Exchange contributions to Fock matrix 
  SUBROUTINE process_JK( gout, nij, nkl, pij, pkl, &
                               ni,  nj,  nk,  nl,  &
                               pil, pkj, pjl, pki  )

    !> input 
    integer :: nij, nkl
    integer :: ni, nj, nk, nl 
    real(8) :: gout(nkl,10,nij,10) 

    !> output
    real(8) :: pij(10,nij,2)  !batch data block
    real(8) :: pkl(10,nkl,2)  !full  data block

    real(8) :: pil(10,ni,nl,2)  !output G(il)
    real(8) :: pkj(10,nk,nj,2)  !input  P(kj)
    real(8) :: pjl(10,nj,nl,2)  !output G(jl)
    real(8) :: pki(10,nk,ni,2)  !input  P(ki)


    !> local
    integer :: m, i, j, k, l, ij, kl
    real(8) :: tk, gll, sum0, epsx, epsy, epsz

    !note: consider correct prefactors (alpha^2, alpha^4) for ss-blocks

!ok do k=1,10
    do k=1,1
      do ij=1,nij
        tk = 0.0d0      !fock(k,ij)
        do kl=1,nkl
      
          gll  = gout(kl, 1, ij, k)
!ok       epsz = gout(kl, 5, ij, k) - gout(kl, 3, ij, k)
!ok       epsy = gout(kl, 4, ij, k) - gout(kl, 8, ij, k)
!ok       epsx = gout(kl, 9, ij, k) - gout(kl, 7, ij, k)
!ok       sum0 = gout(kl, 2, ij, k) + gout(kl, 6, ij, k) + gout(kl,10, ij, k) 

          !> coulomb large-component only
          tk = tk + gll*( pkl( 1,kl,1) + pkl( 6,kl,1) )   !Re[(aa,ll) + (bb,ll)] 
!ok     
!ok       tk = tk + gll *( pkl( 1,kl,1) + pkl( 6,kl,1) ) &  !Re[(aa,ll) + (bb,ll)] 
!ok               + sum0*( pkl( 2,kl,1) + pkl( 7,kl,1) ) &  !Re[(aa,ss) + (bb,ss)]
!ok               - epsz*( pkl( 3,kl,1) - pkl( 8,kl,1) ) &  !Im[(aa,ss) - (bb,ss)] - Im*Im              
!ok               + epsy*( pkl( 4,kl,1) - pkl( 9,kl,1) ) &  !Re[(ab,ss) - (ba,ss)]
!ok               - epsx*( pkl( 5,kl,1) + pkl(10,kl,1) )    !Im[(ab,ss) + (ba,ss)] - Im*Im

        enddo
        pij(k,ij,2) = pij(k,ij,2) + tk
      enddo
    enddo

!   do k=1,10
!     do kl=1,nkl
!       tk = 0.0d0      !fock(k,kl)
!       do ij=1,nij
!     
!         gll  = gout(kl, k, ij,  1)
!         epsz = gout(kl, k, ij,  5) - gout(kl, k, ij, 3)
!         epsy = gout(kl, k, ij,  4) - gout(kl, k, ij, 8)
!         epsx = gout(kl, k, ij,  9) - gout(kl, k, ij, 7)
!         sum0 = gout(kl, k, ij,  2) + gout(kl, k, ij, 6) + gout(kl, k, ij, 10) 
!     
!         tk = tk + gll *( pij( 1,ij,1) + pij( 6,ij,1) ) &  !Re[(aa,ll) + (bb,ll)] 
!                 + sum0*( pij( 2,ij,1) + pij( 7,ij,1) ) &  !Re[(aa,ss) + (bb,ss)]
!                 - epsz*( pij( 3,ij,1) - pij( 8,ij,1) ) &  !Im[(aa,ss) - (bb,ss)] - Im*Im              
!                 + epsy*( pij( 4,ij,1) - pij( 9,ij,1) ) &  !Re[(ab,ss) - (ba,ss)]
!                 - epsx*( pij( 5,ij,1) + pij(10,ij,1) )    !Im[(ab,ss) + (ba,ss)] - Im*Im

!       enddo
!       pkl(k,kl,2) = pkl(k,kl,2) + tk
!     enddo
!   enddo



    !> HF exchange contribution

!ok do k=1,10
!   do m=1,1
!     ij = 0
!     do j=1,nj
!       do i=1,ni
!         kl = 0
!         ij = ij + 1
!         do l=1,nl
!           do k=1,nk
!             kl = kl + 1

!             gll = gout(kl, 1, ij, m)

!             pil( 1,i,l,2) = pil( 1,i,l,2) - gll*pkj( 1,k,j,1) !Re[(aa,ll)]
!             pil( 6,i,l,2) = pil( 6,i,l,2) - gll*pkj( 6,k,j,1) !Re[(bb,ll)]

!       enddo
!       pij(m,ij,2) = pij(m,ij,2) + tk
!     enddo
!   enddo

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>
  SUBROUTINE add_J( nij, pij, gij )  

    !> input 
    integer :: nij
    real(8) :: pij(10,nij,2) 

    !> output
    real(8) :: gij(10,nij,2)

    !> local
    integer :: ij 

    do ij=1,nij
      gij( 1,ij,2) = gij( 1,ij,2) + 2.0d0*pij( 1,ij,2)  
      gij( 2,ij,2) = gij( 2,ij,2) + 2.0d0*pij( 2,ij,2)  
      gij( 3,ij,2) = gij( 3,ij,2) + 2.0d0*pij( 3,ij,2)  
      gij( 4,ij,2) = gij( 4,ij,2) + 2.0d0*pij( 4,ij,2)  
      gij( 5,ij,2) = gij( 5,ij,2) + 2.0d0*pij( 5,ij,2)  
      gij( 6,ij,2) = gij( 6,ij,2) + 2.0d0*pij( 6,ij,2)  
      gij( 7,ij,2) = gij( 7,ij,2) + 2.0d0*pij( 7,ij,2)  
      gij( 8,ij,2) = gij( 8,ij,2) + 2.0d0*pij( 8,ij,2)  
      gij( 9,ij,2) = gij( 9,ij,2) + 2.0d0*pij( 9,ij,2)  
      gij(10,ij,2) = gij(10,ij,2) + 2.0d0*pij(10,ij,2)  
    enddo

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE doP( index, GR, GI, idg, P, idp )

    !> input
    integer,              intent(in) :: idg 
    integer,              intent(in) :: idp 
    integer,              intent(in) :: index 
    complex(8), optional, intent(in) :: GR(idg,idg)
    complex(8), optional, intent(in) :: GI(idg,idg)

    !> output
    type(real_vector) :: P(idp,idp)

    !> local
    integer :: ij
    integer :: i, ii, io, icc
    integer :: j, jj, jo, jcc


    if( index == 1 )then
      do i=1,size(gto)
        icc = gto(i)%cdegen
        do j=1,size(gto)
          jcc = gto(j)%cdegen
          allocate( P(i,j)%elms(10*icc*jcc*2) )     !todo: tu je 20 zatial pre coulomb cast iba
        enddo
      enddo 
    endif


    if( present(GR) )then
      do i=1,size(gto)
        io  = gto(i)%offset
        icc = gto(i)%cdegen
        do j=1,size(gto)
          jo  = gto(j)%offset
          jcc = gto(j)%cdegen
          ij  =  index
          !> density-block
          do jj=1,jcc
            do ii=1,icc
              P(i,j)%elms(ij) = dreal( GR( jo+jj, io+ii ) )      !note:  P(ba) -> sa uklada ako (ab)
              ij = ij + 10                                       !fixme: tu je 10 zatial pre coulomb cast iba
            enddo
          enddo
          !> fock-block 
          do jj=1,jcc
            do ii=1,icc
              P(i,j)%elms(ij) = 0.0d0
              ij = ij + 10                                       !fixme: tu je 10 zatial pre coulomb cast iba
            enddo
          enddo
        enddo
      enddo
    elseif( present(GI) )then
      do i=1,size(gto)
        io  = gto(i)%offset
        icc = gto(i)%cdegen
        do j=1,size(gto)
          jo  = gto(j)%offset
          jcc = gto(j)%cdegen
          ij  = index
          !> density-block
          do jj=1,jcc
            do ii=1,icc
              P(i,j)%elms(ij) = dimag( GI( jo+jj, io+ii ) )    !note:  P(ba) -> sa uklada ako (ab) 
              ij = ij + 10                                     !fixme: tu je 10 zatial pre coulomb cast iba
            enddo
          enddo
          !> fock-block
          do jj=1,jcc
            do ii=1,icc
              P(i,j)%elms(ij) = 0.0d0                        
              ij = ij + 10                                    !fixme: tu je 10 zatial pre coulomb cast iba
            enddo
          enddo
        enddo
      enddo
    endif
 
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE doG( block, p, idp, g, idg )

    !> input
    integer,          intent(in) :: idp 
    integer,          intent(in) :: idg 
    character(len=4), intent(in) :: block 
    type(real_vector)            :: p(idp,idp)

    !> output
    complex(8) :: g(idg,idg)

    !> local
    integer :: ij
    integer :: i, ii, io, icc
    integer :: j, jj, jo, jcc
    real(8) :: sum0, epsx, epsy, epsz


    select case( block )
      case('aall','bbll')
        do i=1,size(gto)
          io  = gto(i)%offset
          icc = gto(i)%cdegen
          do j=1,i
            jo  = gto(j)%offset
            jcc = gto(j)%cdegen
            ij  = 10*icc*jcc 
            do jj=1,jcc
              do ii=1,icc
                g( io+ii, jo+jj ) = dcmplx( p(i,j)%elms(ij+1), 0.0d0 ) 
                ij = ij + 10                                       
              enddo
            enddo
          enddo
        enddo

      case('aass')
        do i=1,size(gto)
          io  = gto(i)%offset
          icc = gto(i)%cdegen
          do j=1,i
            jo  = gto(j)%offset
            jcc = gto(j)%cdegen
            ij  = 10*icc*jcc 
            do jj=1,jcc
              do ii=1,icc
                epsz = p(i,j)%elms(ij+5) - p(i,j)%elms(ij+3) 
                sum0 = p(i,j)%elms(ij+2) + p(i,j)%elms(ij+6) + p(i,j)%elms(ij+10)
                g( io+ii, jo+jj ) = dcmplx( sum0, epsz ) 
                ij = ij + 10                                       
              enddo
            enddo
          enddo
        enddo

      case('bbss')
        do i=1,size(gto)
          io  = gto(i)%offset
          icc = gto(i)%cdegen
          do j=1,i
            jo  = gto(j)%offset
            jcc = gto(j)%cdegen
            ij  = 10*icc*jcc 
            do jj=1,jcc
              do ii=1,icc
                epsz = p(i,j)%elms(ij+5) - p(i,j)%elms(ij+3) 
                sum0 = p(i,j)%elms(ij+2) + p(i,j)%elms(ij+6) + p(i,j)%elms(ij+10)
                g( io+ii, jo+jj ) = dcmplx( sum0,-epsz ) 
                ij = ij + 10                                       
              enddo
            enddo
          enddo
        enddo

      case('abss')
        do i=1,size(gto)
          io  = gto(i)%offset
          icc = gto(i)%cdegen
          do j=1,i
            jo  = gto(j)%offset
            jcc = gto(j)%cdegen
            ij  = 10*icc*jcc 
            do jj=1,jcc
              do ii=1,icc
                epsy = p(i,j)%elms(ij+4) - p(i,j)%elms(ij+8) 
                epsx = p(i,j)%elms(ij+9) - p(i,j)%elms(ij+7) 
                g( io+ii, jo+jj ) = dcmplx( epsy, epsx ) 
                ij = ij + 10                                       
              enddo
            enddo
          enddo
        enddo

      case('bass')
        do i=1,size(gto)
          io  = gto(i)%offset
          icc = gto(i)%cdegen
          do j=1,i
            jo  = gto(j)%offset
            jcc = gto(j)%cdegen
            ij  = 10*icc*jcc 
            do jj=1,jcc
              do ii=1,icc
                epsy = p(i,j)%elms(ij+4) - p(i,j)%elms(ij+8) 
                epsx = p(i,j)%elms(ij+9) - p(i,j)%elms(ij+7) 
                g( io+ii, jo+jj ) = dcmplx(-epsy, epsx ) 
                ij = ij + 10                                       
              enddo
            enddo
          enddo
        enddo

      case default
        write(6,*) ' Storage machanism for this block is not implemented yet => stop'
        stop

      end select
       
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
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
END MODULE
