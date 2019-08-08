module make_ibo

contains

   subroutine LocalizeIbosInIaoBasis(CIbOcc,LIbocc,cspinor_fragment,n_fragments,nmo2,nocc,Expnt)
   ! make intrinsic bond orbitals (IBOs) (localized occupied orbitals expressed in terms of IAOs)
   ! The localization procedure is in the spirit of Pipek-Mezey (J. Chem. Phys. 90, 4916 (1989))
   ! It consists in effectively minimizing the # of atoms an orbital is centered on.
   ! For this, unitary operations (rotations) are perform amongst the occupied MO
   ! (leaving the Slater Determinant invariant) until a function (that contains information
   ! about the # of electrons located on the IAO of each atoms) is maximised.

        use exacorr_mo

        implicit none

        logical                 :: conv
        integer                 :: i,j,k,l, imo, nmoi, it
        integer, intent(in)     :: n_fragments, nmo2, nocc, Expnt
        complex(8), intent(in)  :: CIbOcc(nmo2,nocc) ! The occupied MO in terms of IAOs
        type(cmo),  intent(in)  :: cspinor_fragment(n_fragments) ! MO-coefficients of the fragments
        complex(8), intent(out) :: LIbOcc(nmo2,nocc) ! The localized occupied MO in terms of IAOs
        real(8)                 :: phi, func, fGrad, Aij, Bij, Qii, Qij, Qjj, Loc
        complex(8), allocatable :: C(:,:), Ai(:), Aj(:), vec(:)

        LIbOcc = CIbOcc
        conv = .false.
        write(*,*)"Warning : real(Qij) is taken instead of complx(Qij) for now."
        write(*,'(1A5,2A12)') "ITER","LOC(orbs)","Gradient"

        do it=1,2048
          func = 0.D0
          imo = 0
          ! Compute the function in Eq.(4) of Knizia's paper [JCTC, 9, 4834 (2013)]
          ! The goal is to maximize this function. 
          ! Normally it requires between 5 and 10 iterations.
          do i=1,n_fragments
             nmoi = cspinor_fragment(i)%nmo
             do j=1,nocc
                allocate(vec(nmoi))
                vec = LIbOcc(imo+1:imo+nmoi,j)
                ! dot_product should work with complex matrices without any change:
                ! if vec is complex, dot_product = sum(conj(vec)*vec).
                ! dot_product of a complex vector with itself is real.
                func = func + real(dot_product(vec,vec))**dble(Expnt)
                deallocate(vec)
             enddo
             imo = imo + nmoi
          enddo
  
          ! Start the orbital localization procedure:
          fGrad = 0.D0
!          write(*,'(3A5,5A8)') "i","j","Frag","Qii","Qij","Qjj","Aij","Bij"
          ! Sum over all occupied orbital pairs i,j ; j < i
          do i=2,nocc
             do j=1,i-1
                Aij = 0.D0
                Bij = 0.D0
                imo = 0
                ! Increment Aij and Bij for each fragments:
                do k=1,n_fragments
                   nmoi = cspinor_fragment(k)%nmo   
                   allocate(C(nmoi,nocc))
                   allocate(Ai(nmoi))
                   allocate(Aj(nmoi))
                   C = LIbOcc(imo+1:imo+nmoi,:)
                   Ai = C(:,i)
                   Aj = C(:,j)
                   ! Compute the charge matrix elements of fragment k
                   Qii = real(dot_product(Ai,Ai)) ! dot product is real
                   Qij = real(dot_product(Ai,Aj)) ! dot product is complex....
                   Qjj = real(dot_product(Aj,Aj)) ! dot product is real
                   if ((Expnt.eq.2) .or. (Expnt.eq.3)) then
                   ! Pipek-Mezey localization
                     Aij = Aij + 4.D0*Qij**2.D0 - (Qii-Qjj)**2.D0
                     Bij = Bij + 4.D0*Qij*(Qii - Qjj)
                   else if (Expnt.eq.4) then
                   ! Appendix D of Knizia's paper
                     Bij = Bij + 4.D0*Qij*(Qii**3.D0 - Qjj**3.D0)
                     Aij = Aij - Qii**4.D0 - Qjj**4.D0 + &
                         & 6.D0*(Qii**2.D0 + Qjj**2.D0)*Qij**2.D0 + Qjj*Qii**3.D0 + &
                         & Qii*Qjj**3.D0
                   else
                     write(*,*)"STOP: exponent should be 2, 3 or 4."
                     stop
                   end if
!                   write(*,'(3I5,5F8.4)') i,j,k,Qii,Qij,Qjj,Aij,Bij
                   deallocate(C,Ai,Aj)
                   imo = imo + nmoi
                enddo

                ! Now perform the rotation:
                if ((Aij**2.D0 + Bij**2.D0) .gt. 1D-12) then
                  allocate(Ai(nmo2))
                  allocate(Aj(nmo2))
                  ! I don't know how to compute phi for complex Bij, Aij.
                  phi = 0.25D0*datan2(Bij,-Aij)
                  !Bij is the actual gradient. Aij is effectively
                  !the second derivative at phi=0.
                  fGrad = fGrad + Bij**2.D0
                  Ai = LIbOcc(:,i)
                  Aj = LIbOcc(:,j)
                  LIbOcc(:,i) = dcos(phi)*Ai + dsin(phi)*Aj 
                  LIbOcc(:,j) = -dsin(phi)*Ai + dcos(phi)*Aj 
                  deallocate(Ai,Aj)
                endif
             enddo
          enddo

          fGrad = fGrad**0.5D0
          Loc = func**(1.D0/(1.D0*Expnt))
          write(*,"(I5,F12.6,ES12.2)") it,Loc,fGrad
          if (abs(fGrad) .lt. 1.0D-12) then
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
