module make_ibo

   ! references : 
   ! [0] : Pipek and Mezey, J. Chem. Phys. 90, 4916 (1989)
   ! [1] : Knizia, J. Chem. Theory. Comput. 9, 4834 (2013)
   ! [2] : Dubillard et al. J. Chem. Phys. 124, 154307 (2006)

contains

   subroutine Localization(CIbOcc,LIbocc,cspinor_fragment,n_fragments,nmo2,nocc,Expnt)
   ! make REAL intrinsic bond orbitals (IBOs) (localized occupied orbitals expressed in terms of REAL IAOs)
   ! It consists in effectively minimizing the number of atoms an orbital is centered on.
   ! For this, unitary operations (rotations) are perform amongst the occupied real MO
   ! (leaving the Slater Determinant invariant) until a function (that contains information
   ! about the # of electrons located on the IAO of each atoms) is maximised.

        use exacorr_mo

        implicit none

        logical                 :: conv
        integer                 :: i,j,k,l, imo, nmoi, it, orb_deloc(nocc), n_ignored
        integer, intent(in)     :: n_fragments, nmo2, nocc, Expnt
        complex(8), intent(in)  :: CIbOcc(nmo2,nocc) ! The occupied MO in terms of IAOs
        type(cmo),  intent(in)  :: cspinor_fragment(n_fragments) ! MO-coefficients of the fragments
        complex(8), intent(out) :: LIbOcc(nmo2,nocc) ! The localized occupied MO in terms of IAOs
        real(8)                 :: Aij, Bij, phi, func, fGrad, Loc, func0, Delta_func, max_condition
        complex(8)              :: Qii, Qij, Qjj
        complex(8), allocatable :: Ai(:), Aj(:), vec(:)

        LIbOcc = CIbOcc
        conv = .false.
        write(*,*)"------- Localization Procedure (IBO) -------"
        write(*,*)
        write(*,'(1A11,I1)')"Exponent = ",Expnt
!        write(*,*)
!        write(*,*)"Function n_A(i): number of |i>'s electrons located on the IAOs of the fragments."
        imo = 0
        orb_deloc = 0
        func0 = 0.D0
        do i=1,n_fragments
!           write(*,"(A16,I2,A7)")"------ FRAGMENT ", i, " ------"
!           write(*,"(A12,A12)")"occ.MO","#.elec"
           nmoi = cspinor_fragment(i)%nmo
           do j=1,nocc
              allocate(vec(nmoi))
              vec = LIbOcc(imo+1:imo+nmoi,j)
              func0 = func0 + real(dot_product(vec,vec),8)**dble(Expnt)
!              write(*,"(I12,F12.6)") j,real(dot_product(vec,vec))
              if (real(dot_product(vec,vec)) > 0.01D0) then
                 orb_deloc(j) = orb_deloc(j) + 1
              endif
              deallocate(vec)
           enddo
           imo = imo + nmoi
        enddo

        write(*,*)
        write(*,"(A12,A12,A20)")"occ.MO","deloc.#frag","(thrs 0.01)"
        do i=1,nocc
           write(*,"(I12,I12)") i,orb_deloc(i)
        enddo

        write(*,*)
        write(*,'(1A5,2A12)') "ITER","Function","Gradient"
        func = 0.D0
        do it=1,2048
          ! Compute the function in Eq.(4) of [1]
          ! Note that it is the same as in Eq.(3) of [2] 
          ! because the IAOs are orthonormalized
          ! The goal is to maximize this function. 
          ! Normally it requires between 5 and 10 iterations.
          call get_func(LIbOcc,cspinor_fragment,n_fragments,nmo2,nocc,Expnt,func)
          Delta_func = func - func0
  
          ! Start the orbital localization procedure:
          fGrad = 0.D0
          ! Sum over all occupied orbital pairs i,j ; j < i
          do i=1,nocc
             do j=i+1,nocc
                ! Aij and Bij are given in Appendix D of [1]
                ! Aij is the second derivative at phi=0 (see Eq.(8) of [2]).
                ! Note that Aij is only the approximate diagonal elements of the Hessian, not the full one.
                ! Bij is the gradient (see g_ij in Eq.(7) of [2])
                Aij = 0.D0
                Bij = 0.D0
                imo = 0
                ! Increment Aij and Bij for each fragments:
                do k=1,n_fragments
                   nmoi = cspinor_fragment(k)%nmo   
                   allocate(Ai(nmoi))
                   allocate(Aj(nmoi))
                   Ai = LIbOcc(imo+1:imo+nmoi,i)
                   Aj = LIbOcc(imo+1:imo+nmoi,j)
                   ! Compute the charge matrix elements of fragment k
                   Qii = dot_product(Ai,Ai) ! correspond to Mii in [3] (or z_i, because I labels a general reciprocal lattice vector 
                   Qij = real(dot_product(Ai,Aj)) ! correspond to Mij in [3]
                   Qjj = dot_product(Aj,Aj) ! correspond to Mjj in [3]
                   if ((Expnt.eq.2) .or. (Expnt.eq.3)) then
                   ! Pipek-Mezey localization
                     Aij = Aij + real(conjg(Qij)*Qij - 0.25D0*conjg((Qii-Qjj))*(Qii-Qjj),8) ! Correspond to b in [3]
                     Bij = Bij + real(Qij*(Qii - Qjj),8) ! Correspond to a in [3]
                   else if (Expnt.eq.4) then
                   ! see Appendix D of [1].
                     Aij = Aij + real( - Qii*Qii*Qii*Qii - Qjj*Qjj*Qjj*Qjj + &
                         & 6.D0*(Qii*Qii + Qjj*Qjj)*conjg(Qij)*Qij + Qjj*Qii*Qii*Qii + &
                         & Qii*Qjj*Qjj*Qjj,8)
                     Bij = Bij + real(4.D0*Qij*(Qii*Qii*Qii - Qjj*Qjj*Qjj),8)
                   else
                     write(*,*)"STOP: exponent should be 2, 3 or 4."
                     stop
                   end if
                   deallocate(Ai,Aj)
                   imo = imo + nmoi
                enddo

                if ((Aij*Aij + Bij*Bij) .lt. 1.0D-12) then
                  write(*,'(A17,I2,A1,I2,A9)')"Rotation of pair ",i,"/",j," ignored."
                else
                  phi = 0.25D0*datan2(Bij,-Aij)
                  !Bij is the actual gradient. Aij is effectively
                  !the second derivative at phi=0.
                  fGrad = fGrad + Bij*Bij
                  allocate(Ai(nmo2))
                  allocate(Aj(nmo2))
                  Ai = LIbOcc(:,i)
                  Aj = LIbOcc(:,j)
                  LIbOcc(:,i) = dcos(phi)*Ai + dsin(phi)*Aj 
                  LIbOcc(:,j) = -dsin(phi)*Ai + dcos(phi)*Aj 
                  deallocate(Ai,Aj)
                  max_condition = 16.D0*Aij*dcos(4.D0*phi) - 16.D0*Bij*dsin(4.D0*phi)
                  if (max_condition > 0.D0) then
                     write(*,*)"Warning: The condition for the maximum is not fulfilled &
                              & (i.e. max_condition is not < 0) for orbital pair i,j = ",i,j
                     write(*,'(A20,ES12.4)') "  max_condition = ", max_condition
                  endif
                endif
             enddo
          enddo

          fGrad = fGrad**0.5D0
          write(*,"(I5,F12.6,ES12.2)") it,func,fGrad
          if (abs(fGrad) .lt. 1.0D-12) then
             conv = .true.
             exit
          endif
        enddo

        write(*,*)
        write(*,*)"*****************************************************"
        if (.not. conv) then
           print*,"No convergence reached, value of fGrad = ", fGrad
        else
           print*,"Final converged values (it, func, fGrad):"
           write(*,"(I5,F12.6,ES12.2)") it,func,fGrad
        endif
        write(*,*)"*****************************************************"

!        write(*,*)
!        write(*,*)"Function n_A(i): number of |i>'s electrons located on the IAOs of the fragments."
        imo = 0
        orb_deloc = 0
        do i=1,n_fragments
!           write(*,"(A16,I2,A7)")"------ FRAGMENT ", i, " ------"
!           write(*,"(A12,A12)")"occ.MO","#.elec"
           nmoi = cspinor_fragment(i)%nmo
           do j=1,nocc
              allocate(vec(nmoi))
              vec = LIbOcc(imo+1:imo+nmoi,j)
!              write(*,"(I12,F12.6)") j,real(dot_product(vec,vec))
              if (real(dot_product(vec,vec)) > 0.01D0) then
                 orb_deloc(j) = orb_deloc(j) + 1
              endif
              deallocate(vec)
           enddo
           imo = imo + nmoi
        enddo

        write(*,*)
        write(*,"(A12,A12,A20)")"occ.MO","deloc.#frag","(thrs 0.01)"
        do i=1,nocc
           write(*,"(I12,I12)") i,orb_deloc(i)
        enddo
        
        write(*,*)
        write(*,*)"------- End of Localization Procedure -------"
   end subroutine Localization

   subroutine get_func(LIbOcc,cspinor_fragment,n_fragments,nmo2,nocc,Expnt,func)
        ! Eq.(3) of [2]
        ! Actually, because the IAOs are orthonormalized, it reduces to
        ! Eq.(4) of [1].

        use exacorr_mo

        integer                 :: i,j,k,l, imo, nmoi
        integer, intent(in)     :: n_fragments, nmo2, nocc, Expnt
        type(cmo),  intent(in)  :: cspinor_fragment(n_fragments) ! MO-coefficients of the fragments
        complex(8), intent(in)  :: LIbOcc(nmo2,nocc) ! The localized occupied MO in terms of IAOs
        real(8), intent(out)    :: func
        complex(8), allocatable :: vec(:)

        func = 0.D0
        imo = 0
        do i=1,n_fragments
           nmoi = cspinor_fragment(i)%nmo
           do j=1,nocc
              allocate(vec(nmoi))
              vec = LIbOcc(imo+1:imo+nmoi,j)
              ! dot_product should work with complex matrices without any change:
              ! if vec is complex, dot_product = sum(conj(vec)*vec).
              ! dot_product of a complex vector with itself is real.
              func = func + real(dot_product(vec,vec),8)**dble(Expnt)
              deallocate(vec)
           enddo
           imo = imo + nmoi
        enddo
   end subroutine get_func

end module make_ibo
