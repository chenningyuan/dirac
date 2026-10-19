submodule (hdf5_interface) hdf5_interface_dirac

contains

    module procedure write_to_checkpoint

       use checkpoint
       use rose_mo, only: convert_mo_basis

       implicit none

!      write the rose section with orbitals to checkpoint

       character*(:), allocatable :: mo_label
       integer                :: nz, nmo_scf(2)
       integer, allocatable   :: shell_id(:)

!      switch to Kramers-restricted format
       nz = 4
       call convert_mo_basis(spinors,nz)

       ! write generic data
       call checkpoint_write ('/result/wavefunctions/rose/method',sdata='projection') ! to be generalized
       call checkpoint_write ('/result/wavefunctions/rose/energy',rdata=spinors%total_energy)
       ! write the mo coefficients and associated data
       mo_label = '/result/wavefunctions/rose/mobasis/'
       call checkpoint_write (mo_label//'orbitals',rdata=reshape(spinors%coeff_q,(/spinors%nao*spinors%nmo*nz/)))
       call checkpoint_write (mo_label//'eigenvalues',rdata=spinors%energy)
       call checkpoint_write (mo_label//'symmetry',idata=spinors%irrep)
       call checkpoint_write (mo_label//'n_basis',idata=spinors%nao)
       call checkpoint_write (mo_label//'n_mo',idata=spinors%nmo)
       call checkpoint_write (mo_label//'nz',idata=nz)
       call checkpoint_write (mo_label//'n_fsym',idata=1)
       ! shell_ids are copied from the scf section, a hack to keep the MO-type simple
       ! this hack should work as long as we do not reduce the occupied space and have no open shells
       call checkpoint_read  ('/result/wavefunctions/scf/mobasis/n_mo',idata=nmo_scf)
       allocate(shell_id(nmo_scf(1)))
       call checkpoint_read  ('/result/wavefunctions/scf/mobasis/shell_id',idata=shell_id)
       call checkpoint_write (mo_label//'shell_id',idata=shell_id(1:spinors%nmo))

    end procedure write_to_checkpoint

end submodule hdf5_interface_dirac
