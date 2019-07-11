module make_ibo

contains

   subroutine MakeIBO(CIbOcc,LIbocc,cspinor_fragment,n_fragments,nmo2,nocc,Expnt)
   ! make intrinsic bond orbitals (IBOs)

        use exacorr_mo

        implicit none

        integer :: i,j
        integer, intent(in)     :: n_fragments, nmo2, nocc, Expnt
        real(8), intent(in)     :: CIbOcc(nmo2,nocc)
        type(cmo),  intent(in)  :: cspinor_fragment(n_fragments) ! MO-coefficients of the fragments
        real(8), intent(  out)  :: LIbOcc(nmo2,nocc)
         
        call LocalizeIbosInIaoBasis(CIbOcc,LIbOcc,cspinor_fragment,n_fragments,nmo2,nocc,Expnt)
        ! In Knizia's program, another call is made to sort the new localized orbitals by energy.
        ! Not sure we need this here... If not, I can make a single subroutine MakeIBO instead of
        ! having an additional LocalizeIbosInIaoBasis subroutine...

   end subroutine MakeIBO

   subroutine LocalizeIbosInIaoBasis(CIbOcc,LIbOcc,cspinor_fragment,n_fragments,nmo2,nocc,Expnt)
        ! Not working for now, I am a bit confused with the imaginary part.
        use exacorr_mo

        implicit none
        integer, intent(in) :: Expnt
        logical :: conv
        integer :: i,j,k,l, imo, nmoi, it
        integer, intent(in) :: n_fragments, nmo2, nocc
        real(8), intent(in) :: CIbOcc(nmo2,nocc)
        real(8), intent(out) :: LIbOcc(nmo2,nocc)
        type(cmo), intent(in) :: cspinor_fragment(n_fragments) ! MO-coefficients of the fragments
        real(8) :: phi,func, fGrad, Aij, Bij, Cii, Cij, Cjj, Loc
        real(8), allocatable :: C(:,:), Ai(:), Aj(:), vec(:)

        LIbOcc = CIbOcc
        conv = .false.
        write(*,'(1A5,2A12)') "ITER","LOC(orbs)","Gradient"
        do it=1,2048
          func = 0.D0
          imo = 0
          do i=1,n_fragments
             nmoi = cspinor_fragment(i)%nmo
             do j=1,nocc
                allocate(vec(nmoi))
                vec = LIbOcc(imo+1:imo+nmoi,j)
                func = func + dot_product(vec,vec)**Expnt
                deallocate(vec)
             enddo
             imo = imo + nmoi
          enddo
  
          fGrad = 0.D0
!          write(*,'(3A5,5A8)') "i","j","Frag","Cii","Cij","Cjj","Aij","Bij"
          do i=2,nocc
             do j=1,i-1 ! j < i
                Aij = 0.D0
                Bij = 0.D0
                imo = 0
                do k=1,n_fragments
                   nmoi = cspinor_fragment(k)%nmo   
                   allocate(C(nmoi,nocc))
                   allocate(Ai(nmoi))
                   allocate(Aj(nmoi))
                   C = LIbOcc(imo+1:imo+nmoi,:)
                   Ai = C(:,i)
                   Aj = C(:,j)
                   Cii = dot_product(Ai,Ai)
                   Cij = dot_product(Ai,Aj)
                   Cjj = dot_product(Aj,Aj)
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
!                   write(*,'(3I5,5F8.4)') i,j,k,Cii,Cij,Cjj,Aij,Bij
                   deallocate(C,Ai,Aj)
                   imo = imo + nmoi
                enddo
                if ((Aij**2 + Bij**2) .lt. 1D-12) then
!                  print*,'no gradient for this rotation.'
                else
                  allocate(Ai(nmo2))
                  allocate(Aj(nmo2))
                  phi = 0.25*datan2(Bij,-Aij)
                  !Bij is the actual gradient. Aij is effectively
                  !the second derivative at phi=0.
                  fGrad = fGrad + Bij**2
                  Ai = LIbOcc(:,i)
                  Aj = LIbOcc(:,j)
                  LIbOcc(:,i) = dcos(phi)*Ai + dsin(phi)*Aj 
                  LIbOcc(:,j) = -dsin(phi)*Ai + dcos(phi)*Aj 
                  deallocate(Ai,Aj)
                endif
             enddo
          enddo

          fGrad = fGrad**0.5
          Loc = func**(1.D0/(1.D0*Expnt))
          write(*,"(I5,F12.6,ES12.2)") it,Loc,fGrad
          if (abs(fGrad) .lt. 1.0D-8) then
             conv = .true.
             exit
          endif
        enddo
  
        if (.not. conv) then
           print*,"No convergence reached, value of fGrad = ", fGrad
        else
           print*,"Final converged values (it, LOC, fGrad):"
           write(*,"(I5,F12.6,ES12.2)") it,Loc,fGrad
        endif

   end subroutine LocalizeIbosInIaoBasis

end module make_ibo
