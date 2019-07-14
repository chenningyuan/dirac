module properties

contains

 subroutine partial_charges(Cocc_IAO,cspinor_fragment,n_fragments,nmo2,nocc)

        use read_xyz
        use exacorr_mo
        use exacorr_utils

        integer, intent(in) :: n_fragments, nmo2, nocc
        complex(8), intent(in) :: Cocc_IAO(nmo2,nocc)
        type(cmo), intent(in)             :: cspinor_fragment(n_fragments)
        type(atom), pointer                :: first_atom,current_atom ! linked list to store coordindate information
        integer :: i,j,k,imo,nmoi,luxyz
        complex(8), allocatable :: Cocc_IAO_frag(:,:)
        real(8) :: partial_charge(n_fragments)
        real(8), allocatable :: vec(:)

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
!           allocate(Cocc_IAO_frag(nmoi,nocc))
!           Cocc_IAO_frag = Cocc_IAO(imo+1:imo+nmoi,:)
!           partial_charge(i) = current_atom%bascharge
!           do j=1,nocc
!            do k=1,nmoi
!              partial_charge(i) = partial_charge(i) - real(conjg(Cocc_IAO_frag(k,j))*Cocc_IAO_frag(k,j))
!            enddo
!           enddo
!           deallocate(Cocc_IAO_frag)
           partial_charge(i) = current_atom%bascharge
           do j=1,nocc
              allocate(vec(nmoi))
              vec = real(Cocc_IAO(imo+1:imo+nmoi,j))
              !  For real orbitals only. For complex, I should consider conjg(vec):
              partial_charge(i) = partial_charge(i) - dot_product(vec,vec)
              deallocate(vec)
           enddo
           write(*,'(A12,F12.2)'),current_atom%name,partial_charge(i)
           imo = imo + nmoi
           current_atom=>current_atom%next
        end do

 end subroutine partial_charges

end module
