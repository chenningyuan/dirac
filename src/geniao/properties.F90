module properties

contains

 subroutine partial_charges(Cocc_IAO,cspinor_fragment,n_fragments,nmo2,nocc)

        use read_xyz
        use exacorr_mo
        use exacorr_utils

        integer                 :: i,j,k,imo,nmoi,luxyz
        integer, intent(in)     :: n_fragments, nmo2, nocc
        complex(8), intent(in)  :: Cocc_IAO(nmo2,nocc)
        type(cmo), intent(in)   :: cspinor_fragment(n_fragments)
        type(atom), pointer     :: first_atom,current_atom ! linked list to store coordindate information
        complex(8), allocatable :: Cocc_IAO_frag(:,:)
        real(8)                 :: partial_charge(n_fragments)
        complex(8), allocatable :: vec(:)

        call get_free_fileunit(luxyz)
        open (luxyz,file='MOLECULE.XYZ',form='formatted')

        call get_atoms(luxyz,first_atom)
        close (luxyz,status='keep')

        write(*,*) " PARTIAL CHARGES:"
        imo=0
        current_atom => first_atom
        print*,
        write(*,'(2A12)') "FRAG.", "Par.Chrg."
        do i = 1, n_fragments
           nmoi = cspinor_fragment(i)%nmo
           partial_charge(i) = current_atom%bascharge
           do j=1,nocc
              allocate(vec(nmoi))
              vec = Cocc_IAO(imo+1:imo+nmoi,j)
              ! dot_product should be valid for complex vectors. Dot product of a complex vector with itself is real.
              partial_charge(i) = partial_charge(i) - real(dot_product(vec,vec))
              deallocate(vec)
           enddo
           write(*,'(A12,F12.4)'),current_atom%name,partial_charge(i)
           imo = imo + nmoi
           current_atom=>current_atom%next
        end do

 end subroutine partial_charges

 subroutine fragment_charges_IAO(IAO,cspinor_fragment,n_fragments,nmo1,nmo2,nocc)

        use read_xyz
        use exacorr_mo
        use exacorr_utils
        use utils

        complex(8), parameter   :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        integer                 :: i,j,k,imo,nmoi,luxyz
        integer, intent(in)     :: n_fragments, nmo1, nmo2, nocc
        complex(8), intent(in)  :: IAO(nmo1,nmo2)
        type(cmo), intent(in)   :: cspinor_fragment(n_fragments)
        type(atom), pointer     :: first_atom,current_atom ! linked list to store coordindate information
        complex(8), allocatable :: IAO_RDM_frag(:,:), IAO_RDM(:,:), HF_RDM(:,:)
        real(8)                 :: partial_charge(n_fragments), total_charge(n_fragments)
        real(8), allocatable    :: frag_occ(:)

        call get_free_fileunit(luxyz)
        open (luxyz,file='MOLECULE.XYZ',form='formatted')

        call get_atoms(luxyz,first_atom)
        close (luxyz,status='keep')

        write(*,*) " PARTIAL CHARGES (second version...):"
        ! 1RDM:
        allocate(HF_RDM(nmo1,nmo1)) ! In the MO basis, should be 0 everywhere except 1 on part of the diagonal (only for nocc)
        HF_RDM = ZERO
        forall (i=1:nocc) HF_RDM(i,i) = ONE
        allocate(IAO_RDM(nmo2,nmo2))
        IAO_RDM = matmul(conjg(transpose(IAO)),matmul(HF_RDM,IAO))
        imo=0
        current_atom => first_atom
        write(*,'(3A12)') "FRAG.", "Tot.Chrg.", "Par.Chrg."
        do i = 1, n_fragments
           nmoi = cspinor_fragment(i)%nmo
           allocate(IAO_RDM_frag(nmoi,nmoi))
           allocate(frag_occ(nmoi))
           IAO_RDM_frag = IAO_RDM(imo+1:imo+nmoi,imo+1:imo+nmoi)
           call diag_matrix(IAO_RDM_frag,nmoi,frag_occ)
           total_charge(i) = sum(frag_occ)
           partial_charge(i) = current_atom%bascharge - total_charge(i)
           write(*,'(A12,F12.4,F12.4)'),current_atom%name,total_charge(i),partial_charge(i)
           imo = imo + nmoi
           current_atom=>current_atom%next
           deallocate(IAO_RDM_frag,frag_occ)
        enddo
          
        deallocate(HF_RDM,IAO_RDM)

 end subroutine fragment_charges_IAO

end module
