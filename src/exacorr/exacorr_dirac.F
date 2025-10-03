!
!    Collection of routines needed for exacorr that use the DIRAC common
!    block infrastructure.
!
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXA_MOLTRA()
C***********************************************************************
C
C     Driver routine for exacorr MOLTRA implementation
C
C     Written by L.Visscher Jan 2017
C
C***********************************************************************

      use exacorr_cfg
      use exacorr_ao_to_mo,only     : exacorr_compute_mo_integrals
      use eti_ao_to_mo,only         : eti_compute_mo_integrals
      use, intrinsic:: ISO_C_BINDING
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
      use interface_to_mpi
#endif

      implicit none

#include "../moltra/dcbtra.h"
! to access tra4_indstr
#include "infpar.h"
! we need to broadcast tra4_indstr
#include "dcborb.h"
! to access nesh (number of positive energy "electronlike" orbitals)
#include "dgroup.h"
#include "dcbham.h"
!information about the group

      character (len=72)   :: orbital_string(4) ! String with range of MO's to be transformed for each index
      integer              :: kvec              ! Dummy variable needed in first call to numlst
      integer, allocatable :: mo_list(:), mo_list_tmp(:)
      integer              :: nmo(4)
      integer              :: i, j, imo, jmo, ntmo, rcw
      integer              :: numproc

      exa_input%print_level=IPRTRA

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
      ! parallel version of DIRAC
C     Transfer necessary information from common blocks
      call interface_mpi_comm_rank(global_communicator,mytid)
      call interface_mpi_comm_size(global_communicator,numproc)
      if (mytid == 0 .and. numproc < 4 .and. .not. exa_input%talsh) then
         print*," The EXATENSOR library needs ar least 4 mpi processes"
         print*," -> switched to serial mode using only TALSH"
         exa_input%talsh = .true.
      end if
      if (mytid == 0) then
         orbital_string = tra4_indstr(:,1) ! This is to be changed once we have an exacorr-specific input driver
         ! instruct slaves to go to the current routine (exa_moltra) as well
         call dirac_parctl( 19 )
         call interface_mpi_bcast_l0(exa_input%talsh,1,0,
     &                               global_communicator)
      else
         call interface_mpi_bcast_l0(exa_input%talsh,1,0,
     &                               global_communicator)
         if (exa_input%talsh) return ! TALSH has no use for MPI slaves
      end if
#else
      ! serial version of DIRAC
      if (exa_input%talsh) then
         ! TALSH version works fine with serial version of DIRAC
         orbital_string = tra4_indstr(:,1) ! This is to be changed once we have an exacorr-specific input driver
      else
         ! EXATENSOR version does not work with serial version of DIRAC
         print*, "Attempting to use exatensor with a non MPI version"
         call quit ("error in setup: exatensor requires MPI")
      end if
#endif

! next block is code is only for master
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
      if (mytid == 0) then
#endif
C
C     Loop over indices to convert string information to list of orbitals
C
C     Count number of orbitals to be transformed for each index 
      ntmo = 0
      do i = 1, 4
        nmo(i) = 0
        call numlst (orbital_string(i),kvec,1,1,nesh(1),1,nmo(i))
        ntmo = ntmo + nmo(i)
      end do

C     Fill the array with orbital indices
      allocate (mo_list_tmp(ntmo))
      imo = 1
      do i = 1, 4
       call numlst
     & (orbital_string(i),mo_list_tmp(imo),nmo(i),1,nesh(1),1,nmo(i))
       imo = imo + nmo(i)
      end do

C     Convert to Kramers unrestricted picture
      ntmo = ntmo * 2 ! spinors
      allocate (mo_list(ntmo))
      imo = 1
      jmo = 1
      do j = 1, 4
        do i = 1, 2*nmo(j)
            mo_list(imo) = (mo_list_tmp(jmo)-1)*2 + i
            imo = imo + 1
        end do
        jmo = jmo + nmo(j)
      end do
      deallocate(mo_list_tmp)
      nmo = nmo * 2 ! Kramers unrestricted

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
      end if
#endif
C
      IF (SPINFR.OR.NZ.EQ.1) THEN
         rcw = 1
      ELSE
         rcw = 2
      ENDIF

      if (.true.) then !Always use talsh for now
        call eti_compute_mo_integrals(nmo,mo_list,
     &      exa_input%talsh_buff,exa_input%executor,
     &      exa_input%print_level,rcw,TH_CHOLESKY)
      else
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
         !slaves need the same info to get started
         call interface_mpi_bcast(ntmo,1,mparid,global_communicator)
         if (mytid /= 0) allocate(mo_list(ntmo))
         call interface_mpi_bcast(nmo,4,mparid,global_communicator)
         call interface_mpi_bcast(mo_list,ntmo,mparid,
     &                            global_communicator)
#endif
         call exacorr_compute_mo_integrals (nmo,mo_list,
     &                            exa_input%print_level)
      end if
C
      !Only master should deallocate, since only master completes talsh
      if (mytid == 0) deallocate (mo_list)

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXA_CC()
C***********************************************************************
C
C     Driver routine for exacorr CC implementation
C
C     Written by L.Visscher Jul 2017
C
C***********************************************************************

      use exacorr_cfg
      use eti_cc,only             : eti_cc_driver
      use eti_cis,only            : eti_cis_driver
      use eti_mp2no,only          : eti_mp2no_driver
      use eti_mp2lap,only         : eti_mp2lap_driver
      use mp2no_driver,only         : get_no_driver
      use exacorr_cc,only           : exacorr_cc_driver
      use exacorr_mp2no,only        : exacorr_mp2no_driver
      use exacorr_printout,only     : print_exacorr_logo
      use exacorr_global
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
      use interface_to_mpi
#endif

      implicit none

      integer              :: mytid, numproc           ! MPI rank and # processors
      integer              :: n_mp2no,i                 

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
      ! parallel version of DIRAC
      ! the slaves do not have their input set yet, so call them
      ! and let them know whether we use talsh or exatensor
      call interface_mpi_comm_rank(global_communicator,mytid)
      call interface_mpi_comm_size(global_communicator,numproc)
      if (mytid == 0 .and. numproc < 4 .and. .not. exa_input%talsh) then
         print*," The EXATENSOR library needs ar least 4 mpi processes"
         print*," -> switched to serial mode using only TALSH"
         exa_input%talsh = .true.
      end if

      if (mytid == 0 ) then
         call print_exacorr_logo
         call dirac_parctl( 20 ) 
      end if
      call interface_mpi_bcast_l0(exa_input%talsh,1,0,
     &                               global_communicator)
      if (mytid /= 0 .and. exa_input%talsh) return
#else
      ! serial version of DIRAC
      if (.not. exa_input%talsh) then
         ! EXATENSOR version does not work with serial version of DIRAC
         print*, "Attempting to use exatensor with a non MPI version"
         call quit ("error in setup: exatensor requires MPI")
      end if
      call print_exacorr_logo
#endif
C
      call set_talsh_only(exa_input%talsh)
      call initialize_global_data
      if (mytid == 0) call print_date("Initialized global data")

! next block is code is only for master
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
      if (mytid == 0) then
#endif

       call make_spinor_list ( exa_input, string_occupied,
     & string_occ_beta, string_virtual, string_vir_beta)

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
      end if
#endif

      if (exa_input%talsh) then
       if (exa_input%do_cis) then
         call eti_cis_driver (exa_input)
       else if (exa_input%do_no) then
         call get_no_driver (exa_input)
       else if (exa_input%mp2no > 0 ) then
         call eti_mp2no_driver (exa_input,n_mp2no)
       elseif (exa_input%mp2no < 0) then
         call eti_mp2lap_driver (exa_input)
       else
         call eti_cc_driver (exa_input)
       end if
      else
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
         call exacc_sync_cw(0,mytid)
#endif
         !if lambda is activated set prop flags
       if (exa_input%lambda) then 
              if (mytid == 0) call set_prop_flags_ccsd()
       end if

       if (exa_input%mp2no > 0) then
              call exacorr_mp2no_driver (exa_input,n_mp2no)
       else
              call exacorr_cc_driver (exa_input) 
       end if 
      end if

      call exacc_deallocate_cw()

      end

       subroutine read_from_dirac (ao_basis,cspinor,section,ierr)

       use checkpoint
       use exacorr_mo
       use exacorr_ao
       use exacorr_global,only: nfunctions
       use exacorr_utils

       implicit none

!      read all positive energy mo coefficients present on CHECKPOINT file
!      check for incompatibilities (use of symmetry) as well and continue
!      by reading also the basis set information into gtos

       character(len=*),intent(in) :: section
       type(basis_set_info_t), intent(inout) :: ao_basis
       type(cmo), intent(inout):: cspinor
       integer, intent(out)   :: ierr 
       integer, parameter     :: NOERROR  = 0 ! no error
       integer, parameter     :: NOTFOUND = 1 ! file not found
       integer, parameter     :: CORRUPT  = 2 ! file unreadable
       integer, parameter     :: SYMMETRY = 3 ! calculation used symmetry
       integer, parameter     :: NOBASIS  = 4 ! no basis set information could be read (old CHECKPOINT?)
       integer, parameter     :: CONTRAC  = 5 ! basis set is contracted (Interest can not handle this)
       integer, parameter     :: INCONSIS = 6 ! CHECKPOINT does not contain consistent information wrt number of basis functions
       character              :: text*74      ! holds text read from CHECKPOINT, only used for printing this
       integer                :: maxopr,inversion
       real(8)                :: total_energy,e_corr
       integer                :: i, j, k, ishell, iprimstart, iprimend
       integer                :: nao, nao_bas, nshells, iao, ibas
       integer                :: nmo, npo, nto, npriexp, nzbuf
       integer, parameter     :: nz = 4 ! we use no symmetry, so always have this
       real(8), allocatable   :: coeff(:),eig(:),coeff_q(:,:,:)
       integer, allocatable   :: ibeig(:),index(:)
       integer, allocatable   :: orbmom(:)
       real(8), allocatable   :: priexp(:),priccf(:),coord(:)
       integer, allocatable   :: ao_map(:,:)
       character(8)           :: label(4)
       logical                :: labeled
       integer                :: n_mo(2),n_po(2),n_basis(2)
       integer,allocatable    :: n_prim(:),n_cont(:)
       character*(:),allocatable :: group
       logical                :: C1_info,sym_info

       ! check whether the MO coefficients were generated with symmetry, we then need to read the C1 form therof
       ! for backwards compatibility we first check whether the symmetry section itself is available
       call checkpoint_query ('/result/symmetry',sym_info)
       if (.not.sym_info) then
          print*, "WARNING: old CHECKPOINT, assuming NOSYM"
          maxopr = 0
          inversion = 1
       else
          call checkpoint_read ('/result/symmetry/maxopr',idata=maxopr)
          call checkpoint_read ('/result/symmetry/inversion',
     &         idata=inversion)
       end if

       group = '/result/wavefunctions/'//trim(section)//'/'
       ierr = NOERROR
       text = trim(section)//' data read from the checkpoint file '
       call checkpoint_read (group//'mobasis/nz',idata=nzbuf)
       call checkpoint_read (group//'mobasis/n_basis',idata=n_basis)
       call checkpoint_read (group//'mobasis/n_mo',idata=n_mo)
       call checkpoint_read (group//'mobasis/n_po',idata=n_po)
       call checkpoint_read (group//'energy',rdata=total_energy)
       if (section .ne. 'scf') then
          call checkpoint_read (group//'e_corr',rdata=e_corr)
          total_energy = total_energy - e_corr
       end if

       ! check symmetry and get dimensions of the coefficient array
       if (inversion == 2) then
          nao = sum(n_basis)
          npo = sum(n_po)
          nmo = sum(n_mo)
          nto = npo+nmo
       else
          nao = n_basis(1) ! number of ao's
          npo = n_po(1)    ! "positronic" (= negative energy) solutions (can be thrown away for our purposes)
          nmo = n_mo(1)    ! "electronic" (= positive energy) solutions (we only need these)
          nto = npo+nmo    ! total number of solutions on this file
       end if
       if (npo>0) print*, "warning: reading 4C CHECKPOINT file" ! this warning will go once we have the full index transform implemented

       allocate(coeff(nao*nto*nz))
       allocate(coeff_q(nao,nto,nz))
       allocate(eig(nto))
       allocate(ibeig(nto))
       allocate(index(nmo))
        
       ! read the mo coefficients and orbital eigenvalues
       if (maxopr == 0) then
        call checkpoint_read (group//'mobasis/orbitals',rdata=coeff)
        call checkpoint_read (group//'mobasis/eigenvalues',rdata=eig)
       else
        ! check whether the C1 formatted coefficients were written
        call checkpoint_query (group//'mobasis/orbitals_C1',C1_info)
        if (.not. C1_info) then
         print*," Activate the ACMOUT option in the SCF stage"
         call quit ("no C1 formatted MO coefficients available ")
        end if
        call checkpoint_read (group//'mobasis/orbitals_C1',rdata=coeff)
        call checkpoint_read (group//'mobasis/eigenvalues_C1',rdata=eig)
       endif
       call checkpoint_read (group//'mobasis/symmetry',idata=ibeig)
       ! fill in index values
       do i = 1, nmo
          index(i) = i ! we set the index compatible with the fact that we store all mo's
       end do

       ! put in cmo in quaternion format
       coeff_q = reshape(coeff,(/nao,nto,nz/))
       cspinor%algebra = 'q'
       call alloc_mo (cspinor,nao,nmo)
       cspinor%total_energy = total_energy
       cspinor%index        = index(npo+1:nto)
       cspinor%coeff_q      = coeff_q(:,npo+1:nto,:)
       cspinor%energy       = eig(npo+1:nto)
       cspinor%symmetry     = ibeig(npo+1:nto)

       deallocate(coeff)
       deallocate(eig)
       deallocate(ibeig)
       deallocate(index)
!      switch to unrestricted spinor picture
       call convert_mo_to_complex(cspinor)

       ! read general AO basis information
       if (maxopr == 0) then
          group = '/input/aobasis/1/' ! SCF was done in C1
       else
          group = '/input/aobasis/3/' ! we need the AO basis expanded from the symmetry-packed form
       end if
       call checkpoint_read 
     & (group//'angular',idata=ao_basis%basis_angular)
       call checkpoint_read (group//'n_ao',idata=nao_bas)
       call checkpoint_read (group//'n_shells',idata=nshells)
       allocate (orbmom(nshells))
       allocate (n_prim(nshells))
       allocate (n_cont(nshells))
       allocate (coord(3*nshells))
!      read information that is defined for each shell
       call checkpoint_read (group//'orbmom',idata=orbmom)
       call checkpoint_read (group//'n_prim',idata=n_prim)
       call checkpoint_read (group//'n_cont',idata=n_cont)
       call checkpoint_read (group//'center',rdata=coord)
!      values of the exponents and coefficients for all shells (packed into a 1-d arrays)
       npriexp = sum(n_prim)
       allocate (priexp(npriexp))
       allocate (priccf(npriexp))
       call checkpoint_read (group//'exponents',rdata=priexp)
       call checkpoint_read (group//'contractions',rdata=priccf)

       ! check for consistency with the MO first dimension
       if (nao < nao_bas) then
          ierr = INCONSIS
          return
       elseif (nao > nao_bas) then
          print*,"WARNING: assuming Levy-Leblond Hamiltonian"
          print*,"Reading fewer basis functions than MO-coefficients"
          allocate (ao_map(2,nao_bas))
          forall (i=1:nao_bas) ao_map(:,i)=i
          call compress_mo (cspinor,ao_map)
       end if

       ! put information in gtos
       ao_basis%nshells       = nshells
       ao_basis%nao           = nao_bas
       nullify(ao_basis%gtos)
       allocate(ao_basis%gtos(nshells))
       iPrimStart = 1
       do ishell = 1, nshells
          ao_basis%gtos(ishell)%orb_momentum = orbmom(ishell)
          ao_basis%gtos(ishell)%n_primitives = n_prim(ishell)
          nullify  (ao_basis%gtos(ishell)%exponent) 
          allocate (ao_basis%gtos(ishell)%exponent(n_prim(ishell)))
          nullify  (ao_basis%gtos(ishell)%coefficient)
          allocate (ao_basis%gtos(ishell)%coefficient(n_prim(ishell)))
          iPrimEnd   = iPrimStart+n_prim(ishell)-1
          ao_basis%gtos(ishell)%exponent   = priexp(iPrimStart:iPrimEnd)
          ao_basis%gtos(ishell)%coefficient= priccf(iPrimStart:iPrimEnd)
          ao_basis%gtos(ishell)%coord = coord(ishell*3-2:3*ishell)
!         In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
          call compress_unc_shell (ao_basis%gtos(ishell),ierr)
          if (ierr == 1) ierr = CONTRAC
          iPrimStart = iPrimStart + n_prim(ishell)
       end do

       deallocate(orbmom)
       deallocate(n_prim)
       deallocate(n_cont)
       deallocate(coord)
       deallocate(priexp)
       deallocate(priccf)

       ! make array with pointers from a basis function index to its shell
       nullify(ao_basis%shell_indices)
       allocate(ao_basis%shell_indices(nao))
       iao = 0
       do ishell = 1, nshells
          do ibas = 1,
     &       nfunctions(
     &       ao_basis%gtos(ishell)%orb_momentum,ao_basis%basis_angular)
             iao = iao + 1
             if (iao.gt.nao) then
                ierr = INCONSIS
                call quit ('Error in read_from_dirac')
             endif
             ao_basis%shell_indices(iao) = ishell
          end do
       end do
       if (iao.gt.nao) then
          ierr = INCONSIS
          call quit ('Error in read_from_dirac')
       endif

       call print_date(trim(text))
       return

  10   ierr = CORRUPT
       return

  20   ierr = NOBASIS
       return

       end subroutine read_from_dirac

       subroutine read_from_dfcoef (ao_basis,cspinor,file_name,ierr)

       use exacorr_mo
       use exacorr_ao
       use exacorr_global,only: nfunctions
       use exacorr_utils

       implicit none

!      read all positive energy mo coefficients present on CHECKPOINT file
!      check for incompatibilities (use of symmetry) as well and continue
!      by reading also the basis set information into gtos

       character(len=*),intent(in) :: file_name
       type(basis_set_info_t), intent(inout) :: ao_basis
       type(cmo), intent(inout):: cspinor
       integer, intent(out)   :: ierr 
       integer, parameter     :: NOERROR  = 0 ! no error
       integer, parameter     :: NOTFOUND = 1 ! file not found
       integer, parameter     :: CORRUPT  = 2 ! file unreadable
       integer, parameter     :: SYMMETRY = 3 ! calculation used symmetry
       integer, parameter     :: NOBASIS  = 4 ! no basis set information could be read (old CHECKPOINT?)
       integer, parameter     :: CONTRAC  = 5 ! basis set is contracted (Interest can not handle this)
       integer, parameter     :: INCONSIS = 6 ! CHECKPOINT does not contain consistent information wrt number of basis functions
       character              :: text*74      ! holds text read from CHECKPOINT, only used for printing this
       integer                :: nsym
       integer                :: idim(3,2)
       real(8)                :: total_energy
       logical                :: tobe
       integer                :: lucoef
       integer                :: i, j, k, ishell, iprimstart, iprimend
       integer                :: nao, nao_bas, nshells, iao, ibas
       integer                :: nmo, npo, nto, npriexp, ncnt, nzbuf
       integer, parameter     :: nz = 4 ! we use no symmetry, so always have this
       real(8), allocatable   :: coeff(:,:,:),eig(:)
       integer, allocatable   :: ibeig(:),index(:)
       integer, allocatable   :: ncent(:),nhkt(:),nuco(:),numcf(:)
       integer, allocatable   :: jstrt(:)
       real(8), allocatable   :: cent(:,:),priexp(:),priccf(:,:)
       integer, allocatable   :: ao_map(:,:)
       character(8)           :: label(4)
       logical                :: labeled

       ierr = NOERROR

       ! open file and read dimension info
       inquire (file=file_name,exist=tobe)
       if (.not.tobe) then
          print*, "error: MO coefficient file"//file_name//" not found"
          ierr = NOTFOUND
          return
       else
          call get_free_fileunit(lucoef)
          open (lucoef,file=file_name,status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
          read (lucoef) label
          if (label(1) == '********') then
             labeled = .true. ! new format
             read (lucoef,end=10,err=10) text,nsym,nzbuf,
     &         ((idim(i,j),i=1,3),j=1,nsym),total_energy
          else
             labeled = .false.
             backspace(lucoef) ! old format (DIRAC19 and earlier)
             read (lucoef,end=10,err=10) text,nsym,
     &         ((idim(i,j),i=1,3),j=1,nsym),total_energy
          endif
       end if

       ! check symmetry and get dimensions of the coefficient array
       if (nsym > 1) then
          ierr = SYMMETRY
          close (lucoef,status='KEEP')
          return
       else
          nao = idim(3,1) ! number of basis functions
          npo = idim(1,1) ! "positronic" (= negative energy) solutions (can be thrown away for our purposes)
          nmo = idim(2,1) ! "electronic" (= positive energy) solutions (we only need these)
          nto = npo+nmo   ! total number of solutions on this file
          if (npo > 0) print*, "warning: reading 4C CHECKPOINT file" ! this warning will go once we have the full index transform implemented
       end if

       allocate(coeff(nao,nto,nz))
       allocate(eig(nto))
       allocate(ibeig(nto))
       allocate(index(nmo))
        
       ! read the mo coefficients and orbital eigenvalues
       if (labeled) read (lucoef) label
       read (lucoef,end=10,err=10) coeff
       if (labeled) read (lucoef) label
       read (lucoef,end=10,err=10) eig
       if (labeled) read (lucoef) label
       read (lucoef,end=10,err=10) ibeig
       ! fill in index values
       do i = 1, nmo
          index(i) = i ! we set the index compatible with the fact that we store all mo's
       end do

       ! put in cmo in quaternion format
       cspinor%algebra = 'q'
       call alloc_mo (cspinor,nao,nmo)
       cspinor%total_energy = total_energy
       cspinor%index        = index(npo+1:nto)
       cspinor%coeff_q      = coeff(:,npo+1:nto,:)
       cspinor%energy       = eig(npo+1:nto)
       cspinor%symmetry     = ibeig(npo+1:nto)

       deallocate(coeff)
       deallocate(eig)
       deallocate(ibeig)
       deallocate(index)
!      switch to unrestricted spinor picture
       call convert_mo_to_complex(cspinor)

       ! read the basis set information and check for consistency
       if (labeled) read (lucoef) label
       read (LUCOEF,err=20,end=20) nao_bas,nshells,npriexp,ncnt 
       if (nao < nao_bas) then
          ierr = INCONSIS
          close (lucoef,status='KEEP')
          return
       elseif (nao > nao_bas) then
          print*,"WARNING: assuming Levy-Leblond Hamiltonian"
          print*,"Reading fewer basis functions than MO-coefficients"
          allocate (ao_map(2,nao_bas))
          forall (i=1:nao_bas) ao_map(:,i)=i
          call compress_mo (cspinor,ao_map)
       end if

       ! allocate the arrays to read the basis set information in DIRAC format
       allocate(ncent(nshells))
       allocate(nhkt(nshells))
       allocate(nuco(nshells))
       allocate(numcf(nshells))
       allocate(jstrt(nshells))
       allocate(cent(nshells,3))
       allocate(priexp(npriexp))
       allocate(priccf(npriexp,ncnt))

!      orb_momentum, atom_number, n_primives for all shells
       read (LUCOEF,err=20,end=20) nhkt,ncent,nuco
!      start of exponent values, number of contracted for all shells
       read (LUCOEF,err=20,end=20) jstrt,numcf
!      coordinates (x, y, z) of center for all shells
       read (LUCOEF,err=20,end=20) cent
!      values of the exponents for all shells (packed into a 1-d array)
       read (LUCOEF,err=20,end=20) priexp
!      values of the contraction coefficients for all shells (packed into a 1-d array)
       read (LUCOEF,err=20,end=20) priccf

       ! put in gtos
       ao_basis%nshells       = nshells
       ao_basis%nao           = nao_bas
       ao_basis%basis_angular = 1
       nullify(ao_basis%gtos)
       allocate(ao_basis%gtos(nshells))
       do ishell = 1, nshells
          ao_basis%gtos(ishell)%orb_momentum = nhkt(ishell)
          ao_basis%gtos(ishell)%atom_number  = ncent(ishell)
          ao_basis%gtos(ishell)%n_primitives = nuco(ishell)
          nullify  (ao_basis%gtos(ishell)%exponent) 
          allocate (ao_basis%gtos(ishell)%exponent(nuco(ishell)))
          nullify  (ao_basis%gtos(ishell)%coefficient)
          allocate (ao_basis%gtos(ishell)%coefficient(nuco(ishell)))
          iPrimStart = jstrt(ishell) + 1
          iPrimEnd   = jstrt(ishell) + nuco(ishell)
          k          = numcf(ishell)
          ao_basis%gtos(ishell)%exponent = priexp(iPrimStart:iPrimEnd)
          ao_basis%gtos(ishell)%coefficient = 
     &                                   priccf(iPrimStart:iPrimEnd,k)
          ao_basis%gtos(ishell)%coord = cent(ishell,1:3)
!         In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
          call compress_unc_shell (ao_basis%gtos(ishell),ierr)
          if (ierr == 1) ierr = CONTRAC
       end do

       deallocate(ncent)
       deallocate(nhkt)
       deallocate(nuco)
       deallocate(numcf)
       deallocate(jstrt)
       deallocate(cent)
       deallocate(priexp)
       deallocate(priccf)

       ! make array with pointers from a basis function index to its shell
       nullify(ao_basis%shell_indices)
       allocate(ao_basis%shell_indices(nao))
       iao = 0
       do ishell = 1, nshells
          do ibas = 1,
     &       nfunctions(
     &       ao_basis%gtos(ishell)%orb_momentum,ao_basis%basis_angular)
             iao = iao + 1
             if (iao.gt.nao) then
                ierr = INCONSIS
                call quit ('Error in read_from_dirac')
             endif
             ao_basis%shell_indices(iao) = ishell
          end do
       end do
       if (iao.gt.nao) then
          ierr = INCONSIS
          close (lucoef,status='KEEP')
          call quit ('Error in read_from_dirac')
       endif
       close (lucoef,status='KEEP')

       call print_date("read from CHECKPOINT:"//text)
       return

  10   ierr = CORRUPT
       close (lucoef,status='KEEP')
       return

  20   ierr = NOBASIS
       close (lucoef,status='KEEP')
       return

       end subroutine read_from_dfcoef

       subroutine write_mos_to_checkpoint(mobasis_id,cspinor,section)

       use exacorr_mo
       use checkpoint

!      write mo coefficients in new dirac format

       integer, intent(in)         :: mobasis_id
       character(len=*),intent(in) :: section
       type(cmo), intent(inout)    :: cspinor
       type(cmo)                   :: qorbital
       integer                     :: ncoeff
       character*(:), allocatable  :: group


       call copy_mo(cspinor,qorbital)
       call convert_mo_to_quaternion(qorbital)

       ncoeff = qorbital%nao * qorbital%nmo * qorbital%nz
       group = '/result/wavefunctions/'//trim(section)//'/mobasis'
       call checkpoint_write 
     & (trim(group)//'/mobasis_id',idata=mobasis_id)
       call checkpoint_write 
     & (trim(group)//'/eigenvalues',rdata=qorbital%energy)
       call checkpoint_write 
     & (trim(group)//'/nz',idata=qorbital%nz)
       call checkpoint_write 
     & (trim(group)//'/n_fsym',idata=1) ! we can not handle symmetry yet, so this is always 1
       call checkpoint_write 
     & (trim(group)//'/n_basis',idata=qorbital%nao)
       call checkpoint_write (trim(group)//'/n_mo',idata=qorbital%nmo)
       call checkpoint_write (trim(group)//'/n_po',idata=0)
       call checkpoint_write 
     & (trim(group)//'/symmetry',idata=qorbital%symmetry)
       call checkpoint_write 
     & (trim(group)//'/orbitals',rdata=
     & reshape(qorbital%coeff_q,(/ncoeff/)))

       call dealloc_mo(qorbital)

       end subroutine  write_mos_to_checkpoint

       subroutine write_mos_to_dfcoef(ao_basis,cspinor,file_name)

       use exacorr_mo
       use exacorr_ao
       use exacorr_utils
       use labeled_storage

       implicit none

!      write mo coefficients in dirac format

       character(len=*),intent(in) :: file_name
       type(basis_set_info_t), intent(inout) :: ao_basis
       type(cmo), intent(inout):: cspinor
       type(cmo)              :: qorbital
       character              :: text*74      ! holds text read from CHECKPOINT, only used for printing this
       integer                :: nsym
       real(8)                :: total_energy
       logical                :: tobe
       integer                :: lucoef
       integer                :: i, ishell, iprimstart, iprimend
       integer                :: nao, nshells, iao, ibas
       integer                :: nmo, npo, nto, npriexp, ncnt, ncoeff
       integer, allocatable   :: numcf(:),jstrt(:)
       real(8), allocatable   :: priexp(:),priccf(:)

!      switch to Kramers-restricted format
       text = "MO coefficients written by exacorr"
       nsym = 1 ! we are not using symmetry, so this is always one
       npo = 0  ! number of "positronic orbitals is hardwired to 0 for now
       call copy_mo(cspinor,qorbital)
       call convert_mo_to_quaternion(qorbital)
       ! print warning in case we overwrite a file and open the file
       call get_free_fileunit(lucoef)
       inquire (file=file_name,exist=tobe)
       if (tobe) then
          print*, "warning: overwriting file "//file_name
          open (lucoef,file=file_name,status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
       else
          open (lucoef,file=file_name,status='NEW',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
       end if

       ! convert basis set info to 1-d lists used by DIRAC (assuming uncontracted basis sets)
       ! TODO: change storage on DIRAC side to use the basis set type
       nshells = ao_basis%nshells
       ! we do not support general contraction so ncnt is always 1
       ncnt = 1
       ! first count the number of exponents (is also number of contraction coefficients) and set jstrt
       allocate(jstrt(nshells))
       allocate(numcf(nshells))
       npriexp = 0
       do ishell = 1, nshells
          jstrt(ishell) = npriexp
          npriexp = npriexp + ao_basis%gtos(ishell)%n_primitives
          numcf(ishell) = 1
       end do
       ! now copy the exponents and contraction coefficients into the 1-d lists
       allocate(priexp(npriexp))
       allocate(priccf(npriexp))
       do ishell = 1, nshells
        iPrimStart = jstrt(ishell) + 1
        iPrimEnd   = jstrt(ishell) + ao_basis%gtos(ishell)%n_primitives
        priexp(iPrimStart:iPrimEnd)=ao_basis%gtos(ishell)%exponent(:)
        priccf(iPrimStart:iPrimEnd)=ao_basis%gtos(ishell)%coefficient(:)
       end do
 
       ! write the header
       call lab_write(lucoef,'INFO    ')
       write (lucoef) text,nsym,qorbital%nz,npo,qorbital%nmo,
     &                qorbital%nao,qorbital%total_energy
       ! write the mo-coefficients
       ncoeff = qorbital%nao * qorbital%nmo * qorbital%nz
       call lab_write(lucoef,'COEFS   ',
     &      reshape(qorbital%coeff_q,(/ncoeff/)))
       ! write the orbital energies
       call lab_write(lucoef,'EVALS   ',qorbital%energy)
       ! write the symnmetry information
       call lab_write(lucoef,'SUPERSYM',idata=qorbital%symmetry)
       ! write the ao basis information header
       call lab_write(lucoef,'BASINFO ')
       write (lucoef) ao_basis%nao,ao_basis%nshells,npriexp,ncnt 
!      orb_momentum, atom_number, n_primives for all shells
       write (lucoef) ao_basis%gtos(:)%orb_momentum,
     &                ao_basis%gtos(:)%atom_number,
     &                ao_basis%gtos(:)%n_primitives
!      start of exponent values, number of contracted for all shells
       write (lucoef) jstrt,numcf
!      coordinates (x, y, z) of center for all shells
       write (lucoef) (ao_basis%gtos(:)%coord(i),i=1,3)
!      values of the exponents for all shells (packed into a 1-d array)
       write (lucoef) priexp
!      values of the contraction coefficients for all shells (packed into a 1-d array)
       write (lucoef) priccf
       call lab_write(lucoef,'EOFLABEL')

       deallocate(jstrt)
       deallocate(numcf)
       deallocate(priexp)
       deallocate(priccf)
       call dealloc_mo(qorbital)

       close (lucoef,status='KEEP')

       return

       end subroutine write_mos_to_dfcoef

       subroutine write_mos_to_dirac_formatted(cspinor,file_name)

       use exacorr_mo
       use exacorr_utils

       implicit none

!      write mo coefficients in dirac format

       character(len=*),intent(in) :: file_name
       type(cmo), intent(inout):: cspinor
       type(cmo)              :: qorbital
       character              :: text*74      ! holds text read from CHECKPOINT, only used for printing this
       integer                :: nsym
       real(8)                :: total_energy
       logical                :: tobe
       integer                :: lucoef
       integer                :: nao
       integer                :: nmo, npo

!      switch to Kramers-restricted format
       text = "MO coefficients produced by geniao "
       npo = 0  ! number of "positronic orbitals is hardwired to 0 for now
       nsym = 1 ! we do not use symmetry, so nsym is always 1
       call copy_mo(cspinor,qorbital)
       call convert_mo_to_quaternion(qorbital)
       ! print warning in case we overwrite a file and open the file
       call get_free_fileunit(lucoef)
       inquire (file=file_name,exist=tobe)
       if (tobe) then
          print*, "warning: overwriting file "//file_name
          open (lucoef,file=file_name,status='OLD',FORM='FORMATTED',
     &     access='SEQUENTIAL')
       else
          open (lucoef,file=file_name,status='NEW',FORM='FORMATTED',
     &     access='SEQUENTIAL')
       end if

       ! write the header
       write (lucoef,'(A74)') text
       write (lucoef,'(7(X,I0))') nsym,npo,qorbital%nmo,qorbital%nao
       write (lucoef,'(E24.16)') qorbital%total_energy
       ! write the mo-coefficients
       write (lucoef,'(6F22.16)') qorbital%coeff_q
       ! write the orbital energies
       write (lucoef,'(6E22.12)') qorbital%energy
       ! write the symnmetry information
       write (lucoef,'(66(X,I0))') qorbital%symmetry

       close (lucoef,status='KEEP')
       call dealloc_mo(qorbital)
       return
       
       end subroutine write_mos_to_dirac_formatted

       subroutine Write_To_Mdcint (nmo,mo_tens,rcw)

!      Utility routine to write to the DIRAC moltra format, makes it easier to test the index transform
 
!      Note that the integral matrix is transposed in moltra, this is why the 3 and 4th indices are written to the left below
!      This routine should therefore be rewritten if the matrix transpose I(ij,kl) .ne. I(kl,ij)

#include "implicit.h"

       integer rcw    ! this is used to toggle between complex (so, rcw=2) and real (scalar, rcw=1) tensors
       integer iclass, nclass
       integer nmo(4)
       real(8) mo_tens(2,nmo(3),nmo(4),nmo(1),nmo(2))
       real(8), allocatable :: fmat(:,:,:),rkl(:,:)
       integer, allocatable :: indk(:),indl(:)
       integer nkr(4),bar(8,4)

       iprint = 0
       irepij = 1
       nkr = nmo / 2 ! go to Kramers restricted
       nclass = 8    ! write the standard classes used in MOLTRA as well

!      Open the MDCINT file
       call mdinti (0) 

!      We only need to write a subset of integrals, the rest will be generated by permutation and/or Kramers' symmetry
 
!      Classes are defined to have a standard subset of integrals, make these consistent with the definition in mdintb (in moltra)
!      Index   :  3    4      1    2
!      Class 1 : (k    l    | i    j   )
!      Class 2 : (kbar lbar | i    j   )
!      Class 3 : (k    lbar | i    jbar)
!      Class 4 : (kbar l    | i    jbar)
!      Class 9 : (kbar l    | i    j   )
!      Class 10: (k    lbar | i    j   )
!      Class 11: (k    l    | ibar j   )
!      Class 12: (k    l    | i    jbar)

!      NB: classes 5 through 8 are not used, we therefore shift the final 4 classes in the array below to get a contiguous index

       bar = 1 ! set everything unbarred to initialize
!      define bars for class 2
       bar(2,3) = 2
       bar(2,4) = 2
!      define bars for class 3
       bar(3,2) = 2
       bar(3,4) = 2
!      define bars for class 4
       bar(4,2) = 2
       bar(4,3) = 2
!      define bars for class 5=9
       bar(5,3) = 2
!      define bars for class 6=10
       bar(6,4) = 2
!      define bars for class 7=11
       bar(7,1) = 2
!      define bars for class 8=12
       bar(8,2) = 2

!      Allocate memory needed by mdintb
       nkl = nkr(3) * nkr(4)
       allocate(rkl(rcw,nkl)) 
       allocate(indk(nkl))
       allocate(indl(nkl)) 
       allocate(fmat(nkr(3),nkr(4),rcw))

!      loop over classes to write out subsets of integrals
       do iclass = 1, nclass
          if (iclass.le.4) then
             jclass = iclass
          else
             jclass = iclass + 4  ! shift iclass to obey the convention used in mdintb (given above)
          end if

!         loop over blocks of integrals, taking i and j as outermost indices like is done in moltra
          do jj = 1, nkr(2)
             j = (jj-1)*2 + bar(iclass,2) ! Spinor index 2
             do ii = 1, nkr(1)
                i = (ii-1)*2 + bar(iclass,1) ! Spinor index 1

!              In moltra the real/imaginary index within a block is last rather than first, so we need to reorganize

!              If the particle/particle (index 12 <-> index 34) symmetry is not fulfilled we should reorganize further
!              this can be done here, but will be quite inefficient and not worth the trouble for this utility routine that
!              is to discarded later anyhow.
     
               do ll = 1, nkr(4)
                  l = (ll-1)*2 + bar(iclass,4) ! Spinor index 4
                  do kk = 1, nkr(3)
                     k = (kk-1)*2 + bar(iclass,3) ! Spinor index 3
                     fmat(kk,ll,1) = mo_tens(1,k,l,i,j)
                     if (rcw.gt.1) fmat(kk,ll,2) = mo_tens(2,k,l,i,j)
                  end do
               end do

!              write a block of integrals to file. mdintb takes care of thresholding and labelling
               call mdintb(iprint,rcw,jclass,irepij,nkr(3),nkr(4),
     &           ii,jj,fmat,nkl,indk,indl,rkl)

             end do
          end do
       end do

!      Close the MDCINT file
       call mdintf 
       deallocate(rkl)
       deallocate(indk)
       deallocate(indl)
       deallocate(fmat)

       return
       end

       subroutine read_from_mrconee (int_1el, ierr)

       use exacorr_1el
       use exacorr_utils, only: get_free_fileunit

       implicit none

!      routine to read the one electron integrals from MRCONEE

       type(one_el_t), intent(out) :: int_1el
       integer, intent(out)   :: ierr
       integer, parameter     :: NOERROR  = 0 ! no error
       integer, parameter     :: NOTFOUND = 1 ! file not found
       integer, parameter     :: CORRUPT  = 2 ! file unreadable
       integer, parameter     :: INCONSIS = 6 ! MRCONEE does not match information read from CHECKPOINT

       integer :: nsp, nkr
       logical :: breit
       real(8) :: e_core
       complex(8), allocatable :: h_core(:,:)

       logical                :: tobe
       integer                :: luone

       ierr=NOERROR

!      open file and read setup infromation
       inquire (file='MRCONEE',exist=tobe)
       if (.not.tobe) then
          print*, "error: file MRCONEE not found"
          ierr = NOTFOUND
          return
       else
          call get_free_fileunit(luone)
          open (luone,file='MRCONEE',status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
          read (luone,ERR=10,END=10) nsp,breit,e_core
       end if
!      save number of spinors
       int_1el%n_spinor=nsp
!      save core energy
       int_1el%e_core=e_core
!TODO  check consistency of e_core and nsp with the values read from CHECKPOINT !

!      skip symmetry information
       read (luone,ERR=10,END=10)
       read (luone,ERR=10,END=10)
       read (luone,ERR=10,END=10)
       read (luone,ERR=10,END=10)

!      allocate array to store the effective 1-body Hamiltonian
       nullify(int_1el%h_core)
       allocate(int_1el%h_core(NSP, NSP))

!      as the ordering on mrconee is different we need an auxilliary array to read in
       allocate(h_core(NSP, NSP))

!      read this "core" Hamiltonian
       read (luone,ERR=20,END=20) h_core

!      reorder and store in the data type
       nkr = nsp / 2 ! number of Kramers pairs is always half the number of spinors
       int_1el%h_core(1:nsp:2,1:nsp:2) = h_core(1:nkr,1:nkr)     ! unbar-unbar
       int_1el%h_core(2:nsp:2,1:nsp:2) = h_core(nkr+1:nsp,1:nkr)   ! bar-unbar
       int_1el%h_core(1:nsp:2,2:nsp:2) = h_core(1:nkr,nkr+1:nsp)   ! unbar-bar
       int_1el%h_core(2:nsp:2,2:nsp:2) = h_core(nkr+1:nsp,nkr+1:nsp) ! bar-bar

       deallocate (h_core)
       close (luone,status='KEEP')
       
       return

  10   ierr = CORRUPT
       return

  20   ierr = INCONSIS
       return

       end subroutine read_from_mrconee

       subroutine read_mdprop (int_1el, ierr)

       use exacorr_1el
       use exacorr_utils, only: get_free_fileunit

       implicit none

!      routine to read property integrals and add these to the int_1el data

       type(one_el_t), intent(inout) :: int_1el
       integer, intent(out)   :: ierr
       integer, parameter     :: NOERROR  = 0 ! no error
       integer, parameter     :: NOTFOUND = 1 ! file not found
       integer, parameter     :: CORRUPT  = 2 ! file unreadable

       integer, parameter     :: max_prop = 100 ! limit maximum of arrays to be read, should be high enough for the foreseeable future
       integer :: nsp, nkr, n_prop, i_prop
       complex(8), allocatable :: h_prop(:,:)

       character*32           :: label
       logical                :: tobe
       integer                :: luone

       ierr=NOERROR

!      open file and read setup infromation
       inquire (file='MDPROP',exist=tobe)
       if (.not.tobe) then
          print*, "error: file MDPROP not found"
          ierr = NOTFOUND
          return
       else
          call get_free_fileunit(luone)
          open (luone,file='MDPROP',status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
       end if
!      set number of spinors
       nsp=int_1el%n_spinor

!      allocate auxilliary array to read in
       allocate(h_prop(nsp, nsp))

!      count the number of matrices that we need to store
!      MDROP has a very simple format: record with label followed by record containing the data
       n_prop = 0
       do while (n_prop < max_prop)
          read (luone,err=10,end=10) label
          if (label(25:32) == 'EOFLABEL') exit
          read (luone,err=10,end=10) h_prop
          n_prop = n_prop + 1
       end do
       rewind (luone)

!      store property integrals
       nullify(int_1el%property_labels)
       nullify(int_1el%h_prop)
       allocate(int_1el%property_labels(n_prop))
       allocate(int_1el%h_prop(nsp,nsp,n_prop))
       int_1el%n_prop = n_prop
       nkr = nsp / 2 ! number of Kramers pairs is always half the number of spinors
       i_prop = 0
       do while (i_prop < n_prop)
          read (luone,err=10,end=10) label
          i_prop = i_prop + 1
          read (luone,err=10,end=10) h_prop
        ! reorder and store in the data type
          int_1el%property_labels(i_prop) = label(25:32)
          int_1el%h_prop(1:nsp:2,1:nsp:2,i_prop)
     &          = h_prop(1:nkr,1:nkr)     ! unbar-unbar
          int_1el%h_prop(2:nsp:2,1:nsp:2,i_prop)
     &         = h_prop(nkr+1:nsp,1:nkr)   ! bar-unbar
          int_1el%h_prop(1:nsp:2,2:nsp:2,i_prop)
     &         = h_prop(1:nkr,nkr+1:nsp)   ! unbar-bar
          int_1el%h_prop(2:nsp:2,2:nsp:2,i_prop)
     &         = h_prop(nkr+1:nsp,nkr+1:nsp) ! bar-bar
       end do

       close (luone,status='KEEP')
       deallocate (h_prop)
       return

  10   ierr = CORRUPT
       return

       end subroutine read_mdprop


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine set_laplace(num,print_level)
C
C---------------Description--------------------------------------------
C
C   set flags for the computation of laplace points
C
C   Johann Pototschnig, Fall 2019

C----------------------------------------------------------------------

      IMPLICIT NONE

C parameters
#include "priunit.h"
#include "lapdim.h"

C common blocks
#include "dcblap.h"


      integer  :: num
      integer  :: print_level

      integer  :: ILAP
      logical  :: LWGHTS, LXPNTS
      CHARACTER*7, CHRERR
      PARAMETER (CHRERR = 'LAPINP>')

C     Initialize /DCBLAP/
C     ===================
C
      IPRLAP = print_level
      NLAP = num
      MXITER = 50
      STEPMX = 0.3D0
      TOLRNG = 1.D-10
      TOLPAR = 1.D-15

      ! rather conservative accuracy threshold
      TOLLAP = 1.D-8

      LWGHTS = .FALSE.
      LXPNTS = .FALSE.

C     Processing
C     ====================

!     use pre-tabulated start values
      INILAP = (.NOT.(LWGHTS.AND.LXPNTS))
!     determine if points are given or not
      FNDLAP = (NLAP.LT.1)
!     print setting
      IF (IPRLAP.GT.1) THEN
       WRITE(LUPRI,*) CHRERR,"print level:",IPRLAP
       WRITE(LUPRI,*) CHRERR,"maximum iteratations:",MXITER
       WRITE(LUPRI,*) CHRERR,"number of quadrature points:",NLAP
       WRITE(LUPRI,*) CHRERR,
     &  "threshold for converging Newton-Maehly algorithms:",TOLRNG 
       WRITE(LUPRI,*) CHRERR,
     &  "threshold for converging Newton-Raphson algorithms:",TOLPAR
       WRITE(LUPRI,*) CHRERR,"accuracy threshold:",TOLLAP 
       WRITE(LUPRI,*) CHRERR,
     &  "step length in line search algorithm:",STEPMX 
       IF (LWGHTS) THEN
        WRITE(LUPRI,*) CHRERR, "start weights were given:"
        DO ILAP = 1,NLAP
         WRITE(LUPRI,*) CHRERR,ILAP,WGHTS(ILAP)
        END DO
       END IF
       IF (LXPNTS) THEN
        WRITE(LUPRI,*) CHRERR, "start exponents were given:"
        DO ILAP = 1,NLAP
         WRITE(LUPRI,*) CHRERR,ILAP,XPNTS(ILAP)
        END DO
       END IF
      END IF

      end subroutine set_laplace

      subroutine print_date(text) 
!      Prints time and date with a string
       character(len=10) :: datex
       character(len=8)  :: timex
       character(len=*)  :: text
       integer           :: lentext
       lentext = len(text)
       call daytime(datex,timex)
       write (*,'(1x,a10,1x,a8,1x,a)') datex,timex,text(1:lentext)
      end subroutine print_date
