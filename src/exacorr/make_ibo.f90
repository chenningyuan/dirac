module make_ibo

contains

   subroutine MakeIBO(A,C,S11,cspinor_fragment,n_fragments,nmo1,nmo2,nocc)
   ! make intrinsic bond orbitals (IBOs)

        use exacorr_mo

        implicit none

        integer :: i,j
        integer, intent(in)     :: n_fragments, nmo1, nmo2, nocc
        complex(8), intent(in)  :: S11(nmo1,nmo1)
        complex(8), intent(in)  :: C(nmo1,nocc)
        complex(8), intent(in)  :: A(nmo1,nmo2)
        type(cmo),  intent(in)  :: cspinor_fragment(n_fragments) ! MO-coefficients of the fragments
        complex(8) :: Aocc(nmo2,nocc), Aocc2(nmo2,nocc)
         
        print*,
        print*,'-----Re-express occupied orbitals in terms of IAO basis-----(',nmo2,'x',nocc,')'
        Aocc = matmul(conjg(transpose(A)),matmul(S11,C))
        do i=1,nmo2
          write(*,'(100G10.2)') (REAL(Aocc(i,j)),j=1,nocc)
        enddo

        print*,
        print*,'-----Localize iteratively (check how to do that properly with complex entries?)-----'
        call LocalizeIbosInIaoBasis(Aocc,Aocc2,cspinor_fragment,n_fragments,nmo2,nocc)
        print*,
        print*,'-----Compute LIbOcc-----(',nmo2,'x',nocc,')'
        do i=1,nmo2
          write(*,'(100G10.2)') (REAL(Aocc2(i,j)),j=1,nocc)
        enddo
        stop 'MakeIBO subroutine not finished'

        ! Python ref:
        !CIbOcc = mdot(CIb.T, S1, COcc)
        ! localize iteratively (Pipek-Mezey-like) and order by energy.
        !LIbOcc = LocalizeIbosInIaoBasis(CIbOcc, MakeAtomIbOffsets(Atoms)[0])
        ! transform back localized orbitals in IAO basis to main basis,
        ! and order them by energy.
        !LOcc = OrderByEnergy(mdot(CIb, LIbOcc), Fock)

   end subroutine MakeIBO

   subroutine LocalizeIbosInIaoBasis(Aocc,Aocc2,cspinor_fragment,n_fragments,nmo2,nocc)
        ! Not working for now, I am a bit confused with the imaginary part.
        use exacorr_mo

        implicit none
        integer, parameter :: Expnt = 4
        logical :: conv
        integer :: i,j,k,l, imo, nmoi, it
        integer, intent(in) :: n_fragments, nmo2, nocc
        complex(8), intent(in) :: Aocc(nmo2,nocc)
        complex(8), intent(out) :: Aocc2(nmo2,nocc)
        type(cmo), intent(in) :: cspinor_fragment(n_fragments) ! MO-coefficients of the fragments
        real(8) :: phi,func, fGrad, Aij, Bij, Cii, Cij, Cjj
        real(8), allocatable :: C(:,:), Ai(:), Aj(:), vec(:)

        print*,"Entering LocalizeIbosInIaoBasis subroutine&
              &, not working yet"
        Aocc2 = Aocc
        conv = .false.
        do it=1,100
          func = 0.D0
          imo = 0
          do i=1,n_fragments
             nmoi = cspinor_fragment(i)%nmo
             do j=1,nocc
                allocate(vec(nmoi))
                vec = Aocc2(imo+1:imo+nmoi,j)
                func = func + real(dot_product(vec,vec))**Expnt
                deallocate(vec)
             enddo
             imo = imo + nmoi
          enddo
          print*,'func = ',func
  
          fGrad = 0.D0
          do i=1,nocc
             do j=1,i
                Aij = 0.D0
                Bij = 0.D0
                imo = 0
                do k=1,n_fragments
                   nmoi = cspinor_fragment(k)%nmo   
                   allocate(C(nmoi,nocc))
                   allocate(Ai(nmoi))
                   allocate(Aj(nmoi))
                   C = Aocc2(imo+1:imo+nmoi,:)
                   Ai = C(:,i)
                   Aj = C(:,j)
                   Cii = real(dot_product(Ai,Ai))
                   Cij = real(dot_product(Ai,Aj))
                   Cjj = real(dot_product(Aj,Aj))
                   if ((Expnt.eq.2) .or. (Expnt.eq.3)) then
                     Aij = Aij + 4.D0*Cij**2 - (Cii-Cjj)**2
                     Bij = Bij + 4*Cij*(Cii - Cjj)
                   else if (Expnt.eq.4) then
                     Bij = Bij + 4.D0*Cij*(Cii**3 - Cjj**3)
                     Aij = Aij - Cii**4 - Cjj**4 + &
                         & 6.D0*(Cii**2 + Cjj**2)*Cij**2 + Cjj*Cii**3 + &
                         & Cii*Cjj**3
                   else
                     write(*,*)"Expnt should be 2, 3 or 4."
                     stop
                   end if
                   deallocate(C,Ai,Aj)
                   imo = imo + nmoi
                enddo
                allocate(Ai(nmo2))
                allocate(Aj(nmo2))
                phi = 0.25*datan2(Bij,-Aij)
                fGrad = fGrad + Bij**2
                Ai = Aocc2(:,i)
                Aj = Aocc2(:,j)
                Aocc2(:,i) = dcos(phi)*Ai + dsin(phi)*Aj 
                Aocc2(:,j) = -dsin(phi)*Ai + dcos(phi)*Aj 
                deallocate(Ai,Aj)
             enddo
          enddo

          fGrad = fGrad**0.5
          if (abs(fGrad) .lt. 1.0D-8) then
             conv = .true.
             exit
          endif
        enddo
  
        if (.not. conv) then
           print*,"No convergence reached, value of fGrad = ", fGrad
        else
           print*,"Iteration ", it, " and fGrad = ", fGrad 
        endif

   end subroutine LocalizeIbosInIaoBasis

end module make_ibo
