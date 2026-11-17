subroutine change_aobasis(ao_basis_new,ao_basis_old,spinor_old,spinor_new)

        use rose_ao
        use rose_mo
        use linear_algebra

        implicit none
        type(basis_set_info_t), intent(in) :: ao_basis_new, ao_basis_old
        type(mo_basis), intent(in)         :: spinor_old
        type(mo_basis), intent(out)        :: spinor_new
        real(8), allocatable               :: Snew(:,:), Snew_inv(:,:) ! the AO overlap matrix of ao_basis_new
        real(8), allocatable               :: Snew_old(:,:) ! the AO overlap matrix between ao_basis_1 and ao_basis_2
        complex(8), allocatable            :: SC(:,:) ! Overlap times MO coefficients
        real(8),    allocatable            :: Snew_tmp(:,:,:), Snew_old_tmp(:,:,:)

        call alloc_mo_basis(spinor_new,ao_basis_new%nao,spinor_old%nmo,spinor_old%type)
        allocate (Snew(ao_basis_new%nao,ao_basis_new%nao))
        allocate (Snew_old(ao_basis_new%nao,ao_basis_old%nao))
        call get_ao_overlap (ao_basis_new,ao_basis_new,Snew)
        call get_ao_overlap (ao_basis_new,ao_basis_old,Snew_old)

        spinor_new%nao = ao_basis_new%nao

        ! C_new = Snew^-1 Snew_old C_old
        allocate(Snew_inv(spinor_new%nao,spinor_new%nao))
        allocate(Snew_tmp(spinor_new%nao,spinor_new%nao,1))
        Snew_tmp = 0.D0
        Snew_tmp(:,:,1) = Snew
        call matrix_inv(Snew_tmp)
        Snew_inv = Snew_tmp(:,:,1)
        deallocate(Snew_tmp)

        allocate(SC(spinor_new%nao,spinor_old%nmo))
        if (spinor_new%type .eq. 1) then
           SC = matmul(Snew_old,spinor_old%coeff_r(:,:,1))
           spinor_new%coeff_r(:,:,1) = matmul(Snew_inv,SC)
        else if (spinor_new%type .eq. 2) then
           SC = matmul(Snew_old,spinor_old%coeff_r(:,:,1))
           spinor_new%coeff_r(:,:,1) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_r(:,:,2))
           spinor_new%coeff_r(:,:,2) = matmul(Snew_inv,SC)
        else if (spinor_new%type .eq. 3) then
           SC = matmul(Snew_old,spinor_old%coeff_c(:,:,1))
           spinor_new%coeff_c(:,:,1) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_c(:,:,2))
           spinor_new%coeff_c(:,:,2) = matmul(Snew_inv,SC)
        else if (spinor_new%type .eq. 4) then
           SC = matmul(Snew_old,spinor_old%coeff_q(:,:,1))
           spinor_new%coeff_q(:,:,1) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_q(:,:,2))
           spinor_new%coeff_q(:,:,2) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_q(:,:,3))
           spinor_new%coeff_q(:,:,3) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_q(:,:,4))
           spinor_new%coeff_q(:,:,4) = matmul(Snew_inv,SC)
        else
           print*,"Error in the dimension defining the algebra."
           stop
        end if
        deallocate(SC)
        deallocate(Snew,Snew_inv,Snew_old)

end subroutine change_aobasis