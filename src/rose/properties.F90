module properties

 use rose_global

 interface print_orbital_delocalization
    module procedure print_orbital_delocalization_iao
    module procedure print_orbital_delocalization_cmo
 end interface

contains

 subroutine molecular_partial_charges(Cocc_IAO,spinor_fragment,fragments,restricted,atom_fragments)
        ! Population analysis. We use symmetrically orthogonal orbitals so this is like Lowdin population analysis, which is much less basis-dependent than Mulliken population analysis.

        use read_xyz
        use rose_mo
        use linear_algebra

        logical, intent(in)        :: restricted, atom_fragments
        integer                    :: i,j,imo,nmoi,nel_core_i
        integer                    :: n_fragments, nmo2, nocc, rcq
        real(8), intent(in)        :: Cocc_IAO(:,:,:) ! assumed shape array
        type(mo_basis), intent(in) :: spinor_fragment(:) ! assumed shape array
        type(molecule_info), intent(in) :: fragments(:) ! assumed shape array
        real(8), allocatable       :: partial_charge(:), population(:)
        real(8), allocatable       :: vec(:,:)

        n_fragments = size(fragments)
        rcq  = size(Cocc_IAO,3)
        nmo2 = size(Cocc_IAO,1)
        nocc = size(Cocc_IAO,2)
        
        write(luout,*)
        write(luout,*) "####################################"
        write(luout,*) "#  Partial charges of fragments    #"
        write(luout,*) "####################################"
        write(luout,*)

        allocate(partial_charge(n_fragments))
        imo=0
        write(luout,'(A15,A12)') "Fragment", "Charge"
        do i = 1, n_fragments
           nmoi = spinor_fragment(i)%nmo
           nel_core_i = fragments(i)%ncore
           partial_charge(i) = fragments(i)%Znumber - nel_core_i
           do j=1,nocc
              allocate(vec(nmoi,rcq))
              allocate(population(rcq)) ! number of electrons associated with orbital j
              vec = Cocc_IAO(imo+1:imo+nmoi,j,:)
              call dot_prod(vec,vec,population)
              if (restricted) then
                 partial_charge(i) = partial_charge(i) - 2*population(1) ! population(1) is the only non-zero contribution of the population array.
              else
                 partial_charge(i) = partial_charge(i) - population(1) ! population(1) is the only non-zero contribution of the population array.
              end if
              deallocate(vec,population)
           enddo
           if (atom_fragments) then
              write(luout,'(A15,F12.6)')trim(fragments(i)%firstatom%name),partial_charge(i)
           else
              write(luout,'(A15,F12.6)')trim(fragments(i)%name),partial_charge(i)
           endif
           imo = imo + nmoi
        end do
        write(luout,*)
        deallocate(partial_charge)

 end subroutine molecular_partial_charges

 subroutine molecular_fragment_charges_IAO(IAO,spinor_fragment,fragments,nocc,restricted)

        use read_xyz
        use rose_mo
        use linear_algebra
        use check

        logical, intent(in)        :: restricted
        integer                    :: i,imo,nmoi,nmo1,nmo2,rcq,n_fragments,nelec
        integer, intent(in)        :: nocc
        real(8), intent(in)        :: IAO(:,:,:) ! assumed shape array
        type(mo_basis), intent(in) :: spinor_fragment(:)
        type(molecule_info), intent(in) :: fragments(:)
        real(8), allocatable       :: IAO_RDM_frag(:,:,:), IAO_RDM(:,:,:), HF_RDM(:,:,:), transconjgIAO(:,:,:)
        real(8), allocatable       :: partial_charge(:), total_charge(:)
        real(8), allocatable       :: frag_occ(:)

        write(luout,*)
        write(luout,*) "------ PARTIAL CHARGES OF MOLECULAR FRAGMENTS (1-RDM method) ------"

        n_fragments = size(fragments)
        rcq  = size(IAO,3)
        nmo1 = size(IAO,1)
        nmo2 = size(IAO,2)
        nelec = nocc
        if (restricted) nelec = 2*nocc

        ! 1RDM:
        allocate(HF_RDM(nmo1,nmo1,rcq)) ! In the MO basis, should be 0 everywhere except 1 on part of the diagonal (only for nocc)
        HF_RDM = 0.D0
        forall (i=1:nocc) HF_RDM(i,i,1) = 1.D0 
 !       write(luout,*)"(HF)"
!        call purity_indicator(HF_RDM,nelec) ! should be 0 for idempotent matrices
        
        allocate(transconjgIAO(nmo2,nmo1,rcq))
        allocate(IAO_RDM(nmo2,nmo2,rcq))
        call transpose_conjg(IAO,transconjgIAO)
        call matrix_mul(transconjgIAO,HF_RDM,IAO,IAO_RDM)
        deallocate(transconjgIAO)
!        write(luout,*)"(IAO)"
!        call purity_indicator(IAO_RDM,nelec) ! should be 0 for idempotent matrices

        imo=0
        write(luout,'(3A12)') "FRAG.", "Tot.Chrg.", "Par.Chrg."
        allocate(partial_charge(n_fragments))
        allocate(total_charge(n_fragments))
        do i = 1, n_fragments
           nmoi = spinor_fragment(i)%nmo
           allocate(IAO_RDM_frag(nmoi,nmoi,rcq))
           allocate(frag_occ(nmoi))
           IAO_RDM_frag = IAO_RDM(imo+1:imo+nmoi,imo+1:imo+nmoi,:)
!           write(luout,*)"(IAO_frag)",i
!           call purity_indicator(IAO_RDM_frag,int(fragments(i)%Znumber)) ! Closer it is to zero, better we can treat this fragment as isolated wrt the whole molecule.
!           call diag_matrix(IAO_RDM_frag,frag_occ)
           call diag_matrix(IAO_RDM_frag,frag_occ)
           if (restricted) then
              total_charge(i) = 2*sum(frag_occ)
           else
              total_charge(i) = sum(frag_occ)
           end if
           partial_charge(i) = fragments(i)%Znumber - fragments(i)%ncore - total_charge(i)
           write(luout,'(A15,F12.6,F12.6)')fragments(i)%name,total_charge(i),partial_charge(i)
           imo = imo + nmoi
           deallocate(IAO_RDM_frag,frag_occ)
        enddo
        write(luout,*)  
        deallocate(HF_RDM,IAO_RDM)
        deallocate(partial_charge,total_charge)

 end subroutine molecular_fragment_charges_IAO

 subroutine oneRDM_fragment_purity_from_Cocc_IAO(Cocc_IAO,spinor_fragment,fragments,restricted)

        use read_xyz
        use rose_mo
        use linear_algebra
        use check

        logical, intent(in)        :: restricted
        real(8), intent(in)        :: Cocc_IAO(:,:,:) ! assumed shape array
        type(mo_basis), intent(in) :: spinor_fragment(:)
        type(molecule_info), intent(in) :: fragments(:)
        real(8), allocatable       :: IAO_RDM_frag(:,:,:),IAO_RDM(:,:,:),Cocc_IAO_bar(:,:,:)
        integer                    :: i,imo,nmoi,nmo2,rcq,n_fragments,nelec,nocc

        write(luout,*)
        write(luout,*) "------ DEVIATION FROM IDEMPOTENCY OF THE RDMs ------"

        n_fragments = size(fragments)
        nmo2  = size(Cocc_IAO,1)
        nocc  = size(Cocc_IAO,2)
        rcq   = size(Cocc_IAO,3)
        nelec = nocc
        if (restricted) nelec = 2*nocc

        allocate(Cocc_IAO_bar(nocc,nmo2,rcq))
        allocate(IAO_RDM(nmo2,nmo2,rcq))
        call transpose_conjg(Cocc_IAO,Cocc_IAO_bar)
        call matrix_mul(Cocc_IAO,Cocc_IAO_bar,IAO_RDM)
        if (restricted) IAO_RDM = 2*IAO_RDM
        deallocate(Cocc_IAO_bar)
        print*,"Full matrix:"
        call purity_indicator(IAO_RDM,nelec) ! should be 0 for idempotent matrices

        imo=0
        do i = 1, n_fragments
           nmoi = spinor_fragment(i)%nmo
           allocate(IAO_RDM_frag(nmoi,nmoi,rcq))
           IAO_RDM_frag = IAO_RDM(imo+1:imo+nmoi,imo+1:imo+nmoi,:)
           print*,"Fragment number",i
           call purity_indicator(IAO_RDM_frag,int(fragments(i)%Znumber-fragments(i)%ncore)) ! Closer it is to zero, better we can treat this fragment as isolated wrt the whole molecule.
           imo = imo + nmoi
           deallocate(IAO_RDM_frag)
        enddo
        write(luout,*)  
        deallocate(IAO_RDM)

 end subroutine oneRDM_fragment_purity_from_Cocc_IAO

 subroutine print_orbital_delocalization_cmo(LMO,IAO,energies,nmo_frag,text)
        use linear_algebra
        implicit none
        real(8), intent(in)        :: LMO(:,:,:),IAO(:,:,:),energies(:)
        integer, intent(in)        :: nmo_frag(:)
        character*(*), intent(in)  :: text

        real(8),allocatable        :: IAO_bar(:,:,:),LMO_IAO(:,:,:)
        integer                    :: nmo1, nmo2, rcq, nact

        nmo1 = size(IAO,1)
        nmo2 = size(IAO,2)
        rcq  = size(IAO,3)
        nact = size(LMO,2)

        ! Backtransform to IAO basis and then perform analysis
        allocate(IAO_bar(nmo2,nmo1,rcq))
        allocate(LMO_IAO(nmo2,nact,rcq))
        call transpose_conjg(IAO,IAO_bar)
        call matrix_mul(IAO_bar,LMO,LMO_IAO)

        call print_orbital_delocalization_iao(LMO_IAO,nmo_frag,text,energies)

 end subroutine print_orbital_delocalization_cmo

 subroutine print_orbital_delocalization_iao(Cocc,nmo_frag,text,energies)
        use linear_algebra
        implicit none
        real(8), intent(in)        :: Cocc(:,:,:)
        integer, intent(in)        :: nmo_frag(:)
        character*(*), intent(in)  :: text
        real(8), intent(in),optional  :: energies(:)
        integer                    :: n_fragments, nocc, nmoi, imo, i,j, rcq
        integer, allocatable       :: orb_deloc(:)
        real(8), allocatable       :: vec(:,:), population(:)
        character(len=52)          :: list_frag_deloc(size(Cocc,2))
        character(len=52)          :: list_population(size(Cocc,2))
        character(len=66)          :: header
        character(len=5)           :: frag,pop

        nocc = size(Cocc,2)
        rcq  = size(Cocc,3)
        n_fragments = size(nmo_frag)

        write(luout,*)
        write(luout,*)"----- Percentage of localization "//trim(text)//" ----- "
        write(luout,*)

        allocate(orb_deloc(nocc))
        orb_deloc = 0
        list_frag_deloc(:) = '['
        list_population(:) = '['
        do j=1,nocc
           imo = 0
           do i=1,n_fragments
              nmoi = nmo_frag(i)
              allocate(vec(nmoi,rcq))
              allocate(population(rcq))
              vec = Cocc(imo+1:imo+nmoi,j,:)
              call dot_prod(vec,vec,population)
              ! That's 1% of the orbital if restricted = False, otherwise that's 2% of the orbital. (is that true?)
              if (population(1) > 0.01D0) then ! population(1) is the only non-zero contribution of the population array.
                 orb_deloc(j) = orb_deloc(j) + 1
                 write(frag,'(I5)') i
                 write(pop,'(F5.1)') population(1)*100.D0
                 list_frag_deloc(j) = trim(list_frag_deloc(j))//frag
                 list_population(j) = trim(list_population(j))//pop
              endif
              deallocate(vec)
              deallocate(population)
              imo = imo + nmoi
           enddo
           list_frag_deloc(j) = trim(list_frag_deloc(j))//trim(']')
           list_population(j) = trim(list_population(j))//trim(']')
        enddo

        header = "Localized on fragments (with %) (thrs = 0.01, first ten are shown)"
        write(luout,'(A4,A)') " MO ",adjustl(header)
        do i=1,nocc
           if (present(energies)) then
              write(luout,'(I3,F12.4,A1,A52)') i,energies(i)," ",adjustl(list_frag_deloc(i))
              write(luout,'(A16,A52)') " ",adjustl(list_population(i))
           else
              write(luout,'(I3,A1,A52)') i," ",adjustl(list_frag_deloc(i))
              write(luout,'(A4,A52)') " ",adjustl(list_population(i))
           endif
        enddo
        deallocate(orb_deloc)
        
 end subroutine print_orbital_delocalization_iao

end module
