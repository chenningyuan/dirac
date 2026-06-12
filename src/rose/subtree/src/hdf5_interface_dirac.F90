submodule (hdf5_interface) hdf5_interface_dirac

contains

    module procedure write_to_checkpoint

       use checkpoint
       use rose_mo, only: convert_mo_basis

       implicit none

!      write the rose section with orbitals to checkpoint

       character*(:), allocatable :: mo_label
       integer                :: nz

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
       call checkpoint_write (mo_label//'shell_id',idata=spinors%shell_id)

    end procedure write_to_checkpoint

end submodule hdf5_interface_dirac
