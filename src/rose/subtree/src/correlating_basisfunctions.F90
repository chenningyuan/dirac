module correlating_basisfunctions

    private

    public read_correlation_basis

    interface
      module subroutine read_correlation_basis(molecule,ao_basis)
        use read_xyzfile, only: molecule_info
        use rose_ao, only: basis_set_info_t
        implicit none
        type(molecule_info), intent(in) :: molecule
        type(basis_set_info_t), intent(inout) :: ao_basis(:)
      end subroutine read_correlation_basis
    end interface

end module correlating_basisfunctions
