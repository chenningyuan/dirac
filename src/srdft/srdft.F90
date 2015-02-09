!dirac_copyright_start
!      Copyright (c) 2014 by the authors of DIRAC.
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

module srdft 

! Manu: - this module contains all necessary functionality 
!           to perform MCSCF-srDFT calculations  
!
!           written by Emmanuel Fromager, august 2010 

! space for some notes on srDFT: 
!
!
!
! ---------------------------------------------------------------------------------------
!                            "roadmap" of the srDFT code
! ---------------------------------------------------------------------------------------
!
!   
!
!    
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       1. 
!
!       2. 
!
!       3.
!
!
!   +++++++++++++++++++++++++++++++++++++++++++++
!
! ---------------------------------------------------------------------------------------
!                            end of "roadmap" of the srDFT code
! ---------------------------------------------------------------------------------------
!
  use memory_allocator
  use xcint_main 
  use num_grid_gen
! x2c modules
! use x2c_fio, only:      &
!     x2c_read (name of the subroutine)

  implicit none

  private

  save

  real(8), parameter, public :: srdft_d1     =  1.0D0
  real(8), parameter, public :: srdft_dm1    = -1.0D0

contains


#ifdef DEBUG_SRDFT

!**********************************************************************
  subroutine addvsrh(work,                   &
                     lwork,                  &
                     nao,                    &
                     ninactive,              &
                     fao,                    &
                     dcao,                   &
                     dvao)

!**********************************************************************
!
!    purpose: add the short-range Hartree potential V^{sr}_H[D]=J^{sr}*D to the Fock matrix "fao"
!             where D is the total density matrix (dcao+dvao)  
!
!**********************************************************************
     integer, intent(in)    :: lwork
     real(8), intent(inout) :: work(lwork)
     integer, intent(in)    :: nao               !  number of basis functions 
     integer, intent(in)    :: ninactive         !  number of inactive orbitals 
!----------------------------------------------------------------------
     integer                :: isymop
     integer                :: ihrmop
     integer                :: ifckop    
     integer                :: nao2              !  number of basis functions squared 
    

!**********************************************************************

       nao2= nao * nao


!      allocate memory for the srH potential J^{sr}*D

!      If there are inactive orbitals ...
!      compute the total density matrix and store it in dcao: (dcao, dvao=) --> (dcao=dcao+dvao,dvao) 
       if (ninactive.gt.0) then
           call daxpy(nao2,srdft_d1,dvao,1,dcao,1)
       end if

!      compute the srH potential J^{sr}*D 

       isymop=1
       ihrmop=1
       ifckop=1

       call twofck(isymop,ihrmop,ifckop,
     &             fao,WORK(KDAO),
     &            NFMAT,WORK(KPOS),IOPT_INTFLG,IPRTWO,WORK(KFREE),LFREE)



  end subroutine

#endif

!**********************************************************************
  subroutine addsrxc(nao,                    &
                     nz,                     &
                     ninactive,              &
                     fao,                    &
                     dcao,                   &
                     dvao,                   &  
                     energy) 

!**********************************************************************
!
!    purposes: 1. add the short-range XC (srXC) potential V^{sr}_xc[D] to the Fock matrix "fao"
!                 where D is the total density matrix (dcao+dvao)
!              2. add the srXC correction 
!                 E^{sr}_{xc}[D] - tr (V^{sr}_xc[D] * dvao) 
!                 to the energy   
!
!**********************************************************************
     real(8), intent(inout) :: fao(*)
     real(8), intent(inout) :: dcao(*)           ! will be used to store the
!                                                  total density matrix 
     real(8), intent(in)    :: dvao(*)
     real(8), intent(inout) :: energy 
     integer, intent(in)    :: nao               !  number of basis functions 
     integer, intent(in)    :: nz               
     integer, intent(in)    :: ninactive         !  number of inactive orbitals 
!----------------------------------------------------------------------

     integer                :: nao2z             !  number of basis functions squared * nz 
     real(8), allocatable   :: vsrxc(:)
     real(8)                :: vsrxcdv           !  srXC correction tr(V^{sr}_xc[D] * dvao) 
     real(8),      external :: ddot

!**********************************************************************


       nao2z = nao*nao*nz

!      allocate memory for the srXC potential vsrxc  
       call alloc(vsrxc,nao2z,"vsrxc-addsrxc")
      
!      initialize vsrxc to zero
       vsrxc=0

!      If there are inactive orbitals ...
!      compute the total density matrix and store it in dvao: (dcao, dvao) --> (dcao, dvao=dcao+dvao) 
       if (ninactive.gt.0) then
       call daxpy(nao2z,srdft_d1,dcao,1,dvao,1)
       end if
 
!!!!!! CHECK THE SCALING FACTOR: density matrices normalized to 1 ...


!      compute the srXC potential 
       call generate_num_grid(dvao)
       call integrate_xc(xc_mat_dim          = nao,  &
                              xc_nz           = nz,   &
                              xc_dmat_0       = dvao, &! contains the TOTAL density matrix
                              xc_nr_dmat      = 0,    &
                              xc_nr_fmat      = 1,    &
                              xc_fmat         = vsrxc,&
                              xc_do_potential = .true.)
      
!      add V^{sr}_xc[D] to the Fock operator fao
       call daxpy(nao2z,srdft_d1,vsrxc,1,fao,1)

!      add the srXC energy to "energy" 
       energy = energy + get_xc_energy()   
       
!      If there are inactive orbitals ...
!      restore the active density matrix dvao
       if (ninactive.gt.0) then
       call daxpy(nao2z,srdft_dm1,dcao,1,dvao,1)
       end if

!      CHECK THE SCALING FACTOR FOR DVAO !!!!!!!

!      compute tr (V^{sr}_xc[D] * dvao) ... 
       vsrxcdv = ddot(nao2z,dvao,1,vsrxc,1) 


!      ... and substract it from the energy   
       energy = energy - vsrxcdv

        
!      release memory
       call dealloc(vsrxc) 


  end subroutine

end module
