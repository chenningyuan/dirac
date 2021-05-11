!#define DEBUG_SOC
!dirac_copyright_start
!      Copyright (c) 2019 by the authors of DIRAC.
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
!
! !> this module contains utility routines related to the two-electron Fock matrix calculation
! !> within the X-AMFI / X2C workflow.
!
! written by stknecht july 2020
!
module x_fock

  use num_grid_gen
  use interface_mo_specific

  implicit none

  public get_scso_2e_fock

  private

contains

!----------------------------------------------------------------------
  subroutine get_scso_2e_fock(                       &
                             dmat,                   &
                             cmo,                    &
                             fmat,                   &
                             nrows,                  &
                             ncols,                  &
                             nr_2e_fock_matrices,    &
                             intflg,                 &
                             nz,                     &
                             ipqtoq,                 &
                             x_nopen,                &
                             x_dfopen,               &
                             x_n2tmt,                &
                             nzt_x,                  &
                             nr_fsym,                &
                             N2BBASXQ_dim_x,         &
                             nr_tmo,                 &
                             nr_ao_all,              &
                             ioff_tmot,              &
                             ioff_tmt,               &
                             ioff_aomat_x,           &
                             myAOMOMAT,              &
                             isDFT,                  &
                             lwork_ext,              &
                             print_lvl               &
                             )
!**********************************************************************
!    purpose: get 2e-Fock matrix with integral contributions according to
!             intflg (bit packed integer):
!             (LL|LL)                     == 1
!                       (SS|LL)           == 2
!             (LL|LL) + (SS|LL)           == 3
!                                 (SS|SS) == 4
!             (LL|LL)           + (SS|SS) == 5
!                       (SS|LL) + (SS|SS) == 6
!             (LL|LL) + (SS|LL) + (SS|SS) == 7
!
!             with GAUNT (LS|LS) entering the game (intflg=8) you know what to do by now...
!
! !> in : density matrix (inactive + active ones) in AO basis
  !> out: 2e-Fock matrix in AO-basis                            
!**********************************************************************
#ifdef VAR_MPI
#include "mpif.h"
#endif
     real(8), intent(inout)             :: dmat(nrows,ncols,nz,nr_2e_fock_matrices)
     real(8), intent(inout)             :: cmo(*)
     real(8), intent(out)               :: fmat(nrows,ncols,nz,nr_2e_fock_matrices)
     integer, intent(in)                :: nr_2e_fock_matrices
     integer, intent(in)                :: nrows
     integer, intent(in)                :: ncols
     integer, intent(in)                :: nz
     integer, intent(in)                :: intflg
     integer, intent(in)                :: ipqtoq(4,0:7)
     integer, INTENT(IN)                :: x_nopen
     real*8 , INTENT(IN)                :: x_dfopen(0:x_nopen)
     integer, INTENT(IN)                :: x_n2tmt
     integer, INTENT(IN)                :: nzt_x
     integer, INTENT(IN)                :: nr_fsym
     integer, INTENT(IN)                :: N2BBASXQ_dim_x
     integer, INTENT(IN)                :: nr_tmo(nr_fsym)
     integer, INTENT(IN)                :: nr_ao_all(nr_fsym)
     integer, INTENT(IN)                :: ioff_tmot(nr_fsym)
     integer, INTENT(IN)                :: ioff_tmt(nr_fsym)
     integer, INTENT(IN)                :: ioff_aomat_x(nr_fsym,nr_fsym)
     integer, intent(in)                :: lwork_ext
     integer, intent(in)                :: print_lvl
     character(len=10), intent(in)      :: myAOMOMAT
     logical, intent(in)                :: isDFT
!----------------------------------------------------------------------
     integer                            :: i, npos, local_lwork, j, kwork, kfree, lfree
     integer                            :: kdtemp, kdevec, kfmo, kfao, ktmat, ksaomo, islutmat
     real(8), allocatable               :: work(:)
     integer, allocatable               :: isymop(:), ifckop(:), ihrmop(:), pos(:)
     logical, allocatable               :: saveflags(:)
     logical                            :: parcal
!----------------------------------------------------------------------

       parcal = .false.
#ifdef VAR_MPI
       call mpi_initialized(parcal,i)
#endif
       call SetIntTaskArrayDimension(npos,parcal)

       local_lwork = lwork_ext - nrows*ncols*nz*nr_2e_fock_matrices - 4 - npos - 3*nr_2e_fock_matrices - 10000
       allocate(work(local_lwork))

       !> 2e-fock matrix construction
       allocate(isymop(nr_2e_fock_matrices))
       allocate(ifckop(nr_2e_fock_matrices))
       allocate(ihrmop(nr_2e_fock_matrices))

       ihrmop = 0
       ifckop = 0
       isymop = 0

       do i = 1,nr_2e_fock_matrices
!        totally symmetric operator
         isymop(i) = 1
!        fock matrix type
         ifckop(i) = 1
!        hermitian operator
         ihrmop(i) = 1
       end do

       allocate(saveflags(4))
       call SaveTaskDistribFlags(saveflags)
       call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /))
       allocate(pos(npos))


!      2e-fock matrices (Dirac driver routine)
       call twofck(isymop,ihrmop,ifckop,fmat,dmat,nr_2e_fock_matrices,   &
                   pos,intflg,print_lvl,work,local_lwork)

       if(parcal) call SetTaskDistribFlags(saveflags)
       deallocate(pos)
       deallocate(saveflags)
       deallocate(ihrmop)
       deallocate(ifckop)
       deallocate(isymop)

       if(nr_2e_fock_matrices == 1) goto 9999

       !> Add active two-electron Fock matrices to closed-shell Fock matrix
       do i = 1,x_nopen
         !> The factor is DF and not D2*DF because of the multiplication with D2 in TWOFCK.
         call daxpy(nrows*ncols*nz,x_dfopen(i),fmat(1,1,1,i+1),1,fmat(1,1,1,1),1)
       end do

       !> Transform AO Fock to MO Fock (taking care of coupling of open shells)
       !> ---------------------------------------------------------------------

       KWORK = 1; KFREE = KWORK; LFREE = local_lwork
       CALL MEMGET2('REAL','TMAT' ,KTMAT ,x_N2TMT,WORK,KFREE,LFREE)
       CALL MEMGET2('REAL','SAOMO',KSAOMO,x_N2TMT,WORK,KFREE,LFREE)

       CALL OPNFIL(isLUTMAT,trim(myAOMOMAT),'OLD','XAMFI ')

       !> read AO2MO transformation matrix TMAT and LUTMAT rec#2 with SAOMO
       CALL READT(isLUTMAT,x_N2TMT,work(ktmat)); CALL READT(isLUTMAT,x_N2TMT,work(KSAOMO))
       CLOSE(isLUTMAT,STATUS='KEEP')

       CALL MEMGET2('REAL','DTEMP',KDTEMP,N2BBASXQ_dim_x*nr_2e_fock_matrices,WORK,KFREE,LFREE)
       CALL MEMGET2('REAL','DEVEC',KDEVEC,N2BBASXQ_dim_x,WORK,KFREE,LFREE)
       CALL MEMGET2('REAL','FMO  ',KFMO  ,N2BBASXQ_dim_x,WORK,KFREE,LFREE)
       CALL MKMOFK(work(kFMO),fmat,dmat,work(kTMAT),cmo,work(KSAOMO),WORK(KDTEMP),WORK(KDEVEC),WORK(KFREE),LFREE)

       !> Backtransform MO Fock to AO Fock
       !> --------------------------------

       CALL OPNFIL(isLUTMAT,trim(myAOMOMAT),'OLD','XAMFI ')
       CALL READT(isLUTMAT,x_N2TMT,WORK(KTMAT)); CALL READT(isLUTMAT,x_N2TMT,work(KSAOMO))
       CLOSE(isLUTMAT,STATUS='KEEP')

       CALL MEMGET2('REAL','FAO  ',KFAO  ,nrows*ncols*nz,WORK,KFREE,LFREE)

       do i = 1,nr_fsym
         CALL QTRANS('MOAO','S',0.0d0,nr_ao_all(I),nr_ao_all(I),                              &
                     nr_tmo(I),nr_tmo(I),                                                     &
                     work(kfAO+ioff_aomat_x(I,I)),nrows,ncols,NZ,IPQTOQ(1,0),                 &
                     WORK(kfmo  +ioff_tmot(I)),nr_tmo(I)   ,nr_tmo(I),NZ,IPQTOQ(1,0),         &
                     WORK(ksaomo+ioff_tmt(I) ),nr_ao_all(I),nr_tmo(I),nzt_x,IPQTOQ(1,0),    &
                     WORK(ksaomo+ioff_tmt(I) ),nr_ao_all(I),nr_tmo(I),nzt_x,IPQTOQ(1,0),    &
                     WORK(KFREE),LFREE,print_lvl)
       end do

       fmat = 0; call dcopy(nrows*ncols*nz,work(kfao),1,fmat(1,1,1,1),1)

9999   deallocate(work)

       if(isDFT)then
          !print *, 'blubb DFT on'
          call interface_mo_write()
          call generate_num_grid(dmat)
#ifdef VAR_MPI
          if (parcal) call dirac_parctl( XCINT_PAR )
#endif
          call xcint_potential_rks(       &
                                   nrows, &
                                   dmat,  &
                                   fmat)
       end if 

  end subroutine get_scso_2e_fock
!----------------------------------------------------------------------

  

 end module
