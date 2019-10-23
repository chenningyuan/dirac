module energy

contains

 subroutine HF_energy(E_HF,Hcore,Fock,rdm,nmo,nuc_rep)
      integer, intent(in   ) :: nmo
      real(8), intent(in   ) :: nuc_rep
      real(8), intent(  out) :: E_HF
      complex(8), intent(in) :: Hcore(nmo,nmo)
      complex(8), intent(in) :: Fock(nmo,nmo)
      complex(8), intent(in) :: rdm(nmo,nmo)
      complex(8),allocatable :: tmp(:,:)

      E_HF = 0.D0
      allocate(tmp(nmo,nmo))
      tmp=matmul(conjg(transpose(Hcore+Fock)),rdm)
      do i=1,nmo
        E_HF = E_HF + real(tmp(i,i))
      enddo
      E_HF = 0.5*E_HF + nuc_rep

      deallocate(tmp)
 end subroutine HF_energy

end module
