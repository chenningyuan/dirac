!#define DEBUG_SOC
!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
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
  use interface_ao_specific
  use interface_mo_specific

  implicit none

  public get_scso_2e_fock
  public get_x_fock

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
                             nr_mo_all,              &
                             ioff_tmot,              &
                             ioff_tmt,               &
                             ioff_aomat_x,           &
                             myAOMOMAT,              &
                             isDFT,                  &
                             isLoA,                  &
                             xcmode,                 &
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
#include "dirac_partask.h"
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
     integer, INTENT(IN)                :: nr_mo_all(nr_fsym)
     integer, INTENT(IN)                :: ioff_tmot(nr_fsym)
     integer, INTENT(IN)                :: ioff_tmt(nr_fsym)
     integer, INTENT(IN)                :: ioff_aomat_x(nr_fsym,nr_fsym)
     integer, intent(in)                :: lwork_ext
     integer, intent(in)                :: xcmode
     integer, intent(in)                :: print_lvl
     character(len=10), intent(in)      :: myAOMOMAT
     logical, intent(in)                :: isDFT
     logical, intent(in)                :: isLoA
!----------------------------------------------------------------------
     integer                            :: i, local_lwork, j, k, kwork, kfree, lfree
     integer                            :: kdtemp, kdevec, kfmo, kfao, ktmat, ksaomo, islutmat
     real(8), allocatable               :: work(:)
     logical                            :: parcal
     character(len=12)                  :: flabel
     character(len=1)                   :: suffix
!----------------------------------------------------------------------

       parcal = .false.
#ifdef VAR_MPI
       call mpi_initialized(parcal,i)
#endif

      call get_x_fock(        &
      dmat,                   &
      fmat,                   &
      nrows,                  &
      ncols,                  &
      nr_2e_fock_matrices,    &
      intflg,                 &
      nz,                     &
      .false.,                &
      lwork_ext,              &
      print_lvl               &
      )


      if(nr_2e_fock_matrices > 1)then

        local_lwork = lwork_ext - 10000
        allocate(work(local_lwork))

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

        isLUTMAT = -1
        CALL OPNFIL(isLUTMAT,trim(myAOMOMAT),'OLD','XAMFI ')

        !> read AO2MO transformation matrix TMAT and LUTMAT rec#2 with SAOMO
        CALL READT(isLUTMAT,x_N2TMT,work(ktmat)); CALL READT(isLUTMAT,x_N2TMT,work(KSAOMO))
        CLOSE(isLUTMAT,STATUS='KEEP')

        CALL MEMGET2('REAL','DTEMP',KDTEMP,N2BBASXQ_dim_x*nr_2e_fock_matrices,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','DEVEC',KDEVEC,N2BBASXQ_dim_x,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','FMO  ',KFMO  ,N2BBASXQ_dim_x,WORK,KFREE,LFREE)
        call dzero(work(KFMO),N2BBASXQ_dim_x)
        call dzero(work(kDEVEC),N2BBASXQ_dim_x)
        call dzero(work(KDTEMP),N2BBASXQ_dim_x*nr_2e_fock_matrices)
        CALL MKMOFK(work(kFMO),fmat,dmat,work(kTMAT),cmo,work(KSAOMO),WORK(KDTEMP),WORK(KDEVEC),WORK(KFREE),LFREE)

        !> Backtransform MO Fock to AO Fock
        !> --------------------------------

        CALL MEMGET2('REAL','FAO  ',KFAO  ,nrows*ncols*nz,WORK,KFREE,LFREE)
        call dzero(work(kfAO),nrows*ncols*nz)

        if(xcmode == 2)then !> two-component mode

            CALL OPNFIL(isLUTMAT,trim(myAOMOMAT),'OLD','XAMFI ')
            CALL READT(isLUTMAT,x_N2TMT,WORK(KTMAT)); CALL READT(isLUTMAT,x_N2TMT,work(KSAOMO))
            CLOSE(isLUTMAT,STATUS='KEEP')

        else !> four-component mode

            ! open scratch file x2cscr
            open(104,file='x2cscr',status='old',form='unformatted',  &
            access='sequential',action="readwrite",position='rewind')

            if(isLoA)then !> linear or atomic symmetry

                call dzero(work(KDTEMP),N2BBASXQ_dim_x*nr_2e_fock_matrices)
                call dzero(work(KDEVEC),N2BBASXQ_dim_x)

                k = 1
                do i = 1, nr_fsym

                     j = nr_mo_all(i)
                     if(j == 0)cycle

                     write(flabel,'(a7,i4,i1)') '4cMOMOl',1,i
                     call x2c_read(flabel,work(KDTEMP+k-1),j**2 * nz,104)

!                    work(KDEVEC) (alias 4c-fock_onmo_nonlin) = work(KDTEMP) * work(kfmo) (alias 4c-fock_onmo) * work(KDTEMP)^+
                     CALL QTRANS('MOAO','S',0.0d0,j,j,                      &
                                  j,j,                                      &
                                  work(kdevec+k-1),j,j,NZ,IPQTOQ(1,0),      & !> out  "AO" matrix
                                  work(kfmo  +k-1),j,j,NZ,IPQTOQ(1,0),      & !> in   "MO" matrix
                                  work(KDTEMP+k-1),j,j,nz,IPQTOQ(1,0),      & !> TMAT
                                  work(KDTEMP+k-1),j,j,nz,IPQTOQ(1,0),      & !> TMAT
                                  WORK(KFREE),LFREE,print_lvl)

!                    update offsets for matrices
                     k = k + j**2 * nz
                end do

                call dzero(work(KFMO),N2BBASXQ_dim_x)
                call dcopy(k-1,work(kdevec),1,work(kfmo),1)

            end if

            !> read the half-transformed overlap matrix
            k = 1
            do i = 1, nr_fsym

               if(nr_mo_all(i) == 0)cycle

               write(flabel,'(a7,i4,i1)') '4cSaomo',1,i
               call x2c_read(flabel,work(ksaomo+k-1),nr_ao_all(i) * nr_mo_all(i) * nzt_x,104)

               k = k + nr_ao_all(i) * nr_mo_all(i) * nzt_x
            end do

            close(104,status="delete")

        end if

        do i = 1,nr_fsym
            CALL QTRANS('MOAO','S',0.0d0,nr_ao_all(I),nr_ao_all(I),                              &
                        nr_tmo(I),nr_tmo(I),                                                     &
                        work(kfAO+ioff_aomat_x(I,I)),nrows,ncols,NZ,IPQTOQ(1,0),                 &
                        WORK(kfmo  +ioff_tmot(I)),nr_tmo(I)   ,nr_tmo(I),NZ,IPQTOQ(1,0),         &
                        WORK(ksaomo+ioff_tmt(I) ),nr_ao_all(I),nr_tmo(I),nzt_x,IPQTOQ(1,0),      &
                        WORK(ksaomo+ioff_tmt(I) ),nr_ao_all(I),nr_tmo(I),nzt_x,IPQTOQ(1,0),      &
                        WORK(KFREE),LFREE,print_lvl)
        end do

        fmat = 0; call dcopy(nrows*ncols*nz,work(kfao),1,fmat(1,1,1,1),1)

        deallocate(work)
      end if

      if(isDFT)then
          !> delete and reset DFT grid
          call delete_num_grid(); call reset_num_grid()
          !> reset interface-mo
          call interface_mo_write()
          !> reset interface-ao
          call interface_ao_write()
          !> new grid
          call generate_num_grid(dmat)
#ifdef VAR_MPI
          if (parcal) call dirac_parctl( XCINT_PAR )
#endif
          allocate(work(nrows*ncols*nz)); work = 0
          call xcint_potential_rks(       &
                                   nrows, &
                                   dmat,  &
                                   work)
          call daxpy(nrows*ncols*nz,1.0d0,work,1,fmat,1)
          !> write F^{Xc,2e} aka V_xc
          write(suffix,'(i0)') xcmode
          open(103,file='amfVxc'//suffix,status='replace',form='unformatted',  &
          access='sequential',action="write",position='rewind')
          write(103) work(1:nrows*ncols*nz)
          close(103,status='keep')
          deallocate(work)
      end if

  end subroutine get_scso_2e_fock
!----------------------------------------------------------------------

  subroutine get_x_fock(    &
    dmat,                   &
    fmat,                   &
    nrows,                  &
    ncols,                  &
    nr_2e_fock_matrices,    &
    intflg,                 &
    nz,                     &
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
#include "dirac_partask.h"
#include "siripc.h"
#endif
real(8), intent(inout)             :: dmat(nrows,ncols,nz,nr_2e_fock_matrices)
real(8), intent(out)               :: fmat(nrows,ncols,nz,nr_2e_fock_matrices)
integer, intent(in)                :: nr_2e_fock_matrices
integer, intent(in)                :: nrows
integer, intent(in)                :: ncols
integer, intent(in)                :: nz
integer, intent(in)                :: intflg
integer, intent(in)                :: lwork_ext
integer, intent(in)                :: print_lvl
logical, intent(in)                :: isDFT
!----------------------------------------------------------------------
integer                            :: i, local_lwork
real(8), allocatable               :: work(:)
integer, allocatable               :: isymop(:), ifckop(:), ihrmop(:)
logical, allocatable               :: saveflags(:)
logical                            :: parcal
character(len=1)                   :: suffix
!----------------------------------------------------------------------
       parcal = .false.
#ifdef VAR_MPI
       call mpi_initialized(parcal,i)
#endif
       local_lwork = lwork_ext - nrows*ncols*nz*nr_2e_fock_matrices - 4 - 3*nr_2e_fock_matrices - 20000
       allocate(work(local_lwork)); work=0;

       !> set common block in Hermit
       CALL PAOVEC(work,local_lwork,0,0)

       !> Prepare for parallel calculation (set inforb)
       CALL RELINF
#ifdef VAR_MPI
       !> enforces a resynchronization of AO common blocks between master and co-workers
       NEWGEO = .true.
#endif

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

!      2e-fock matrices (Dirac driver routine)
       call twofck(isymop,ihrmop,ifckop,fmat,dmat,nr_2e_fock_matrices,   &
                   intflg,print_lvl,work,local_lwork)

       if(parcal) call SetTaskDistribFlags(saveflags)
       deallocate(saveflags)
       deallocate(ihrmop)
       deallocate(ifckop)
       deallocate(isymop)

       deallocate(work)

       if(isDFT)then
          !> delete and reset DFT grid
          call delete_num_grid(); call reset_num_grid()
          !> reset interface-mo
          call interface_mo_write()
          !> reset interface-ao
          call interface_ao_write()
          !> new grid
          call generate_num_grid(dmat)
#ifdef VAR_MPI
          if (parcal) call dirac_parctl( XCINT_PAR )
#endif
          allocate(work(nrows*ncols*nz)); work = 0
          call xcint_potential_rks(       &
                                   nrows, &
                                   dmat,  &
                                   work)
          call daxpy(nrows*ncols*nz,1.0d0,work,1,fmat,1)
          !> write F^{Xc,2e} aka V_xc
          if(intflg == 1)then
            write(suffix,'(i0)') 2
          else
            write(suffix,'(i0)') 4
          end if
          open(103,file='eamfVxc'//suffix,status='replace',form='unformatted',  &
          access='sequential',action="write",position='rewind')
          write(103) work(1:nrows*ncols*nz)
          close(103,status='keep')
          deallocate(work)

        end if

  end subroutine get_x_fock
!----------------------------------------------------------------------

 end module
