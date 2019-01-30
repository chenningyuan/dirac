!This is a stand-alone version of exacorr that reads all information
!from an interface file.

   program exacorr_main

#if defined (VAR_MPI)
      use interface_to_mpi
#endif

      use exacorr_ao_to_mo,only     : exacorr_compute_mo_integrals
      use talsh_ao_to_mo,only       : talsh_compute_mo_integrals
      use talsh_cc,only             : talsh_cc_driver
      use exacorr_cc,only           : exacorr_cc_driver
      use exacc_cfg

      implicit none

      logical mo_transform, cc_calculation
      integer              :: numproc, mytid
      integer              :: imo, nmo(4)
      integer, allocatable :: mo_list(:)
!     Variables for reading of DFCOEF
      logical              :: input_found=.false. !For read_menu_input
      integer              :: file_unit,dfcoef_unit
      logical              :: tobe
!     Variables for storing DFCOEF info
      integer              :: kvec ! Dummy variable needed in first call to numlst
      integer              :: nesh(2) ! the second element is a dummy. Just to keep it consistent with the rest of the code
      character            :: text*74
      integer              :: nsym
      integer              :: idim(3,2)
      integer              :: i, j

!     Initialize MPI variables
      mytid   = 0
      numproc = 0

#if defined (VAR_MPI)
      call interface_mpi_init()
      call interface_mpi_comm_rank(global_communicator,mytid)
      call interface_mpi_comm_size(global_communicator,numproc)
#endif

      mo_transform   = .false.
      cc_calculation = .true.

      if (mytid == 0) then
        call get_free_fileunit(file_unit)
        call read_menu_input('exacorr.inp',file_unit,'**EXACC',input_found) 
        if (.not.input_found) call quit('Input was not found!')
      end if

#if defined (VAR_MPI)
      call interface_mpi_bcast_l0(exa_input%talsh,1,0,global_communicator)
      if (mytid /= 0 .and. exa_input%talsh) call interface_mpi_finalize()
#else
!     serial version of DIRAC
      if (.not.exa_input%talsh) then
!        EXATENSOR version does not work with serial version of DIRAC
         print*, "Attempting to use exatensor with a non MPI version"
         call quit ("error in setup: exatensor requires MPI")
      end if
#endif

       if (mo_transform) then
!      We have no input driver for this test code, so we hardwire this here.
       nmo = 8
       allocate(mo_list(nmo(1)+nmo(2)+nmo(3)+nmo(4)))
       do imo = 1, nmo(1)
           mo_list(imo)                      = imo
           mo_list(imo+nmo(1))               = imo
           mo_list(imo+nmo(1)+nmo(2))        = imo
           mo_list(imo+nmo(1)+nmo(2)+nmo(3)) = imo
       end do
 
!         call talsh_compute_mo_integrals (nmo, mo_list)
          call exacorr_compute_mo_integrals (nmo, mo_list)
       end if

      if (mytid == 0) then
!       open file and read dimension info
        inquire (file='DFCOEF',exist=tobe)
        if (.not.tobe) then
           call quit('error: MO coefficient file DFCOEF not found')
        else
           call get_free_fileunit(dfcoef_unit)
           open (dfcoef_unit,file='DFCOEF',status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')
           read (dfcoef_unit) text,nsym, ((idim(i,j),i=1,3),j=1,nsym)
        end if

!       check symmetry and get dimensions of the coefficient array
        if (nsym > 1) then
           call quit ('exacorr can not handle symmetry')
        else
           nesh(1) = idim(2,1) ! "electronic" (= positive energy) solutions (we only need these)
        end if 
        close (dfcoef_unit,status='KEEP')

        exa_input%nkr_occ = 0
        call numlst (string_occupied,kvec,1,1,nesh(1),1,exa_input%nkr_occ)
        print*, "nocc = ", exa_input%nkr_occ
        allocate(exa_input%mokr_occ(exa_input%nkr_occ))
        call numlst (string_occupied,exa_input%mokr_occ,exa_input%nkr_occ, &
                                            1,nesh(1),1,exa_input%nkr_occ)

        exa_input%nkr_vir = 0
        call numlst (string_virtual,kvec,1,1,nesh(1),1,exa_input%nkr_vir)
        allocate (exa_input%mokr_vir(exa_input%nkr_vir))
        call numlst (string_virtual,exa_input%mokr_vir,exa_input%nkr_vir, &
                                           1,nesh(1),1,exa_input%nkr_vir)
      end if

      if (cc_calculation) then
         if (mytid == 0) print*, "a CC calculation with",exa_input%nkr_occ," occupied"," and ",exa_input%nkr_vir," virtuals"
         if (exa_input%talsh) then
            if (mytid == 0) call talsh_cc_driver (exa_input)
         else
#if defined (VAR_MPI)
!           slaves need the same info to get started
            call interface_mpi_bcast(exa_input%nkr_occ,1,0,global_communicator)
            if (mytid /= 0) allocate(exa_input%mokr_occ(exa_input%nkr_occ))
            call interface_mpi_bcast(exa_input%nkr_vir,1,0,global_communicator)
            if (mytid /= 0) allocate(exa_input%mokr_vir(exa_input%nkr_vir))
            call exacc_sync_cw(0)
#endif
            call exacorr_cc_driver (exa_input)
         end if
      end if

      deallocate(exa_input%mokr_occ)
      deallocate(exa_input%mokr_vir)

#if defined (VAR_MPI) 
      call interface_mpi_finalize()
#endif

   end program exacorr_main
