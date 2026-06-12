submodule (hdf5_interface) hdf5_interface_standalone

contains

    module procedure write_to_checkpoint

       use labeled_storage
       use rose_mo, only: convert_mo_basis

       implicit none

!      write the rose section with orbitals to CHECKPOINT.h5 without DIRAC's checkpoint module

       type(file_info_t)      :: datafile
       character*(:), allocatable :: mo_label
       integer                :: nz

!      switch to Kramers-restricted format
       nz = 4
       call convert_mo_basis(spinors,nz)

       datafile%name = 'CHECKPOINT.h5'
       datafile%type = 2
       datafile%status = 0

       ! write generic data
       call lab_write (datafile,'/result/wavefunctions/rose/method',sdata='projection') ! to be generalized
       call lab_write (datafile,'/result/wavefunctions/rose/energy',rdata=spinors%total_energy)
       ! write the mo coefficients and associated data
       mo_label = '/result/wavefunctions/rose/mobasis/'
       call lab_write (datafile,mo_label//'orbitals',rdata=reshape(spinors%coeff_q,(/spinors%nao*spinors%nmo*nz/)))
       call lab_write (datafile,mo_label//'eigenvalues',rdata=spinors%energy)
       call lab_write (datafile,mo_label//'symmetry',idata=spinors%irrep)
       call lab_write (datafile,mo_label//'n_basis',idata=spinors%nao)
       call lab_write (datafile,mo_label//'n_mo',idata=spinors%nmo)
       call lab_write (datafile,mo_label//'nz',idata=nz)
       call lab_write (datafile,mo_label//'n_fsym',idata=1)

       call lab_write (datafile,mo_label//'shell_id',idata=spinors%shell_id)

    end procedure write_to_checkpoint

end submodule hdf5_interface_standalone
