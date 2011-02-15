!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the
!      research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end
!
!
! module for import/export of matrices from/to C1* symmetry.
!
! written by sknecht feb 2011
!
module x2c_import_exportC1
  use memory_allocator
  use x2cmod_cfg
  use x2c_fio
  use x2c_utils, only: print_x2cmat

  implicit none

  public x2c_export_mat2C1

  private

contains

!**********************************************************************
  subroutine x2c_export_mat2C1(h1_saao)
!----------------------------------------------------------------------
!
!    purpose: export (defining h1) matrix in symmetry (quaternion packed)
!             AO basis to AO basis in C1* double group symmetry.
!
!----------------------------------------------------------------------
     real(8), intent(inout)          :: h1_saao(*)
!--------------------------------------------------------------
     real(8), allocatable            :: scr1_mat(:)
     integer                         :: i
     integer                         :: iq
     integer                         :: ipq
     integer                         :: iz
     integer                         :: irepd
     character (len=12)              :: flabel
!**********************************************************************

!      initialize dimensions and pointers which might be used outside this module in dirac calls
       call setdc2(0)

!      open file with 4c-h1
       open(x2c_fh1dirty,file='x2ch1dirty',status='unknown',form='unformatted', &
            access='sequential',action='readwrite',position='rewind')

       call dzero(h1_saao,nr_quat*nr_ao_total_x2c**2)
       i = 0
       write(flabel,'(a7,i4,i1)') 'h1_4cao',0,0
       call x2c_read(flabel,h1_saao,nr_ao_total_x2c*nr_ao_total_x2c*nr_quat,x2c_fh1dirty)

!      a. insert quaternion phase factors
!      ----------------------------------
       if(nr_quat < 4)then
         do iz = 1, nr_quat
           iq = x2c_cb_pq_to_uq(iz,0)
           call q2bphase('F',iq,1,h1_saao(1+(nr_ao_total_x2c*nr_ao_total_x2c*(iz-1))))
         end do
       end if

!      b. transform hamiltonian integral elements from DIRAC-sorted AO to Hermit-sorted AO basis
!      -----------------------------------------------------------------------------------------
       call bstobu_no_work(h1_saao,nr_quat)

!      allocate scratch matrix
       call alloc(scr1_mat,nr_ao_total_x2c*nr_ao_total_x2c*4)
       call dzero(scr1_mat,nr_ao_total_x2c*nr_ao_total_x2c*4)

!      c. transform from SO-AO basis to AO-basis (nosym)
!      -------------------------------------------------
       i = x2c_bs_to_fs(0,1)
       do iz = 1, 4
         irepd = x2c_bs_irrep_mat(iz,0)
         iq    = x2c_iqmult_trip_q(1,x2c_pointer_quat(irepd,i),iz)
         ipq   = x2c_cb_uq_to_pq(iq,0)
         call mtsoao(h1_saao( 1+(nr_ao_total_x2c*nr_ao_total_x2c*(ipq-1))),    &
                     scr1_mat(1+(nr_ao_total_x2c*nr_ao_total_x2c*(iz-1))),     &
                     nr_ao_total_x2c,                                          &
                     irepd,                                                    &
                     x2c_prt)
       end do
      
!      d. write h4c (nosym) in AO basis to file 
!      ----------------------------------------
       call x2c_write('h1_4caoNOSYM',scr1_mat,nr_ao_total_x2c*nr_ao_total_x2c*4,x2c_fh1dirty)

!      debug print
       if(x2c_prt > 2)then
         call print_x2cmat(scr1_mat,nr_ao_total_x2c,nr_ao_total_x2c,4,x2c_cb_pq_to_uq(1,0),'x2c - h1_4c-AO-nosym',6)
       end if

       call dealloc(scr1_mat)

       close(x2c_fh1dirty, status='keep')

  end subroutine x2c_export_mat2C1
end module
