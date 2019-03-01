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
      use exacorr_global

      implicit none

      logical mo_transform, cc_calculation
      integer              :: numproc, mytid
      integer              :: imo, nmo(4)
      logical              :: input_found=.false.
      integer              :: file_unit
      integer, allocatable :: mo_list(:)
      integer, allocatable :: mokr_occ(:), mokr_vir(:)
      integer, allocatable :: mobe_occ(:), mobe_vir(:)
      integer              :: nesh, kvec
      integer              :: i, ialpha, ibeta

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
!       print logo
        call printtitle
        call printlogo
        call printsubtitle

!       read input
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

!     Initialize global data (read from DFCOEF or similar by master)
      call set_talsh_only(exa_input%talsh)
      call initialize_global_data
      if (mytid == 0) print*," Initialized global data"

      if (mo_transform) then
!        We have no input driver for this test code, so we hardwire this here.
         nmo = 8
         allocate(mo_list(nmo(1)+nmo(2)+nmo(3)+nmo(4)))
         do imo = 1, nmo(1)
             mo_list(imo)                      = imo
             mo_list(imo+nmo(1))               = imo
             mo_list(imo+nmo(1)+nmo(2))        = imo
             mo_list(imo+nmo(1)+nmo(2)+nmo(3)) = imo
         end do

!        call talsh_compute_mo_integrals (nmo, mo_list)
         call exacorr_compute_mo_integrals (nmo, mo_list)
      end if

      if (mytid == 0) then
        nesh = get_nsolutions() / 2
        
        if (.not.exa_input%beta_occ) then
          if (.not.exa_input%beta_vir) then
            exa_input%nkr_occ = 0
            call numlst (string_occupied,kvec,1,1,nesh,1,exa_input%nkr_occ)
            allocate (mokr_occ(exa_input%nkr_occ))
            call numlst (string_occupied,mokr_occ,exa_input%nkr_occ,1,nesh,1,exa_input%nkr_occ)
            exa_input%nkr_occ = 2 * exa_input%nkr_occ
            allocate (exa_input%mokr_occ(exa_input%nkr_occ))
            do i = 1, exa_input%nkr_occ
              exa_input%mokr_occ(i) = (mokr_occ(1)-1)*2+i
            end do
            deallocate (mokr_occ)

            exa_input%nkr_vir = 0
            call numlst (string_virtual,kvec,1,1,nesh,1,exa_input%nkr_vir)
            allocate (mokr_vir(exa_input%nkr_vir))
            call numlst (string_virtual,mokr_vir,exa_input%nkr_vir,1,nesh,1,exa_input%nkr_vir)
            exa_input%nkr_vir = 2 * exa_input%nkr_vir
            allocate (exa_input%mokr_vir(exa_input%nkr_vir))
            do i = 1, exa_input%nkr_vir
              exa_input%mokr_vir(i) = exa_input%mokr_occ(exa_input%nkr_occ)+i
            end do
            deallocate (mokr_vir)

          else
            call quit ("also define .OCC_BETA in EXACORR")
          end if
        else
          if (.not.exa_input%beta_vir) then
            call quit ("also define .VIR_BETA in EXACORR")
          else
            call numlst (string_occupied,kvec,1,1,nesh,1,nmo(1))
            call numlst (string_occ_beta,kvec,1,1,nesh,1,nmo(2))
            call numlst (string_virtual,kvec,1,1,nesh,1,nmo(3))
            call numlst (string_vir_beta,kvec,1,1,nesh,1,nmo(4))

            allocate (mokr_occ(nmo(1)))
            allocate (mobe_occ(nmo(2)))
            allocate (mokr_vir(nmo(3)))
            allocate (mobe_vir(nmo(4)))

            call numlst (string_occupied,mokr_occ,nmo(1),1,nesh,1,nmo(1))
            call numlst (string_occ_beta,mobe_occ,nmo(2),1,nesh,1,nmo(2))
            call numlst (string_virtual,mokr_vir,nmo(3),1,nesh,1,nmo(3))
            call numlst (string_vir_beta,mobe_vir,nmo(4),1,nesh,1,nmo(4))

            if(mokr_occ(1).ne.mobe_occ(1)) then 
              call quit ("core orbital need to be closed shell, change .OCC_BETA in EXACORR")
            end if
            
            exa_input%nkr_occ=nmo(1)+nmo(2)
            exa_input%nkr_vir=nmo(3)+nmo(4)
            allocate (exa_input%mokr_occ(exa_input%nkr_occ))
            allocate (exa_input%mokr_vir(exa_input%nkr_vir))

            imo=0
            do i=1, max(nmo(1), nmo(2))
              if(i.le.nmo(1)) then
                imo=imo+1
                exa_input%mokr_occ(imo)=mokr_occ(i)*2-1
              end if
              if(i.le.nmo(2)) then
                imo=imo+1
                exa_input%mokr_occ(imo)=mobe_occ(i)*2
              end if
            end do

            ialpha=1
            ibeta=1
            imo=0
            do i=1, nmo(3)+nmo(4)
              if(ibeta.gt.nmo(4)) then
                if (ialpha.le.nmo(3)) then
                  imo=imo+1
                  exa_input%mokr_occ(imo)=mokr_occ(ialpha)*2-1
                  ialpha=ialpha+1
                end if
              else if (ialpha.gt.nmo(3)) then
                if (ibeta.le.nmo(4)) then
                  imo=imo+1
                  exa_input%mokr_occ(imo)=mobe_occ(ibeta)*2
                  ibeta=ibeta+1
                end if
              else
                if(mokr_occ(ialpha).eq.mobe_occ(ibeta)) then
                  imo=imo+1
                  exa_input%mokr_occ(imo)=mokr_occ(ialpha)*2-1
                  ialpha=ialpha+1
                  imo=imo+1
                  exa_input%mokr_occ(imo)=mobe_occ(ibeta)*2
                  ibeta=ibeta+1
                else if (mokr_occ(ialpha).lt.mobe_occ(ibeta)) then
                  imo=imo+1
                  exa_input%mokr_occ(imo)=mokr_occ(ialpha)*2-1
                  ialpha=ialpha+1
                else                   
                  imo=imo+1
                  exa_input%mokr_occ(imo)=mobe_occ(ibeta)*2
                  ibeta=ibeta+1
                end if
              end if
            end do

            deallocate (mokr_occ)
            deallocate (mobe_occ)
            deallocate (mokr_vir)
            deallocate (mobe_vir)
          end if
        end if
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
