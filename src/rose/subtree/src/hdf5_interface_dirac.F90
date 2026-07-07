submodule (hdf5_interface) hdf5_interface_dirac

contains

    module procedure write_to_checkpoint

       use checkpoint
       use rose_mo, only: convert_mo_basis

       implicit none

!      write the rose section with orbitals to checkpoint

       character*(:), allocatable :: mo_label
       integer                :: nz,n_basis,n_po

!      switch to Kramers-restricted format
       nz = 4
       if (spinors%has_parent) then
          ! Parent coefficients are already stored in quaternion format; do not convert the active FW coefficients here.
          if (.not. associated(spinors%coeff_parent_q)) stop "parent metadata without parent coefficients in write_to_checkpoint"
          n_basis = spinors%nao_large + spinors%nao_small
       else
          call convert_mo_basis(spinors,nz)
          n_basis = spinors%nao
       endif
       n_po = 0

       ! write generic data
       call checkpoint_write ('/result/wavefunctions/rose/method',sdata='projection') ! to be generalized
       call checkpoint_write ('/result/wavefunctions/rose/energy',rdata=spinors%total_energy)
       ! write the mo coefficients and associated data
       mo_label = '/result/wavefunctions/rose/mobasis/'
       if (spinors%has_parent) then
          call checkpoint_write (mo_label//'orbitals',rdata=reshape(spinors%coeff_parent_q,(/n_basis*spinors%nmo*nz/)))
       else
          call checkpoint_write (mo_label//'orbitals',rdata=reshape(spinors%coeff_q,(/n_basis*spinors%nmo*nz/)))
       endif
       call checkpoint_write (mo_label//'eigenvalues',rdata=spinors%energy)
       call checkpoint_write (mo_label//'symmetry',idata=spinors%irrep)
       call checkpoint_write (mo_label//'n_basis',idata=n_basis)
       call checkpoint_write (mo_label//'n_mo',idata=spinors%nmo)
       call checkpoint_write (mo_label//'n_po',idata=n_po)
       call checkpoint_write (mo_label//'nz',idata=nz)
       call checkpoint_write (mo_label//'n_fsym',idata=1)
       call checkpoint_write (mo_label//'shell_id',idata=spinors%shell_id)

    end procedure write_to_checkpoint

end submodule hdf5_interface_dirac
