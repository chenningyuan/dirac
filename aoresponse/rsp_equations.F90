!Copyright 2009 Andreas J. Thorvaldsen
!This file is made available under the terms of the
!GNU Lesser General Public License.

module rsp_equations
!calculates response densities, using integrals from module prop_integrals
!
!ajt: detailed description of contents to be added here ;-)

    use matrix_defop
    use matrix_genop
    use prop_integ

    implicit none

    public pert_dens    !perturbed (or response-) densities and Fock-matrices
    public pert_fock    !contract perturbed Fock-matrices
    public pert_scf_eq  !contract perturbed scf equation residuals
    !public pert_energy  !contract energy contribution to response functions (not yet impl.)

    private

    logical :: debug_me = .true.

contains

    subroutine pert_fock(p,D, F, c)
        character(*),intent(in)    :: p(:) !perturbation lables p(n)
        type(matrix),intent(in)    :: D(:) !un- and perturbed density matrices
        type(matrix),intent(inout) :: F(*) !output, perturbed Fock matrices
                                          !(*) to permit re-ranking
        integer,optional           :: c(:) !range within each p, c(2*n) (lower,upper indices)
        integer :: n,cc(2*size(p)),b(size(p)),i,pb,pb1,ne; n=size(p)
        call prop_check('PertFock',D(1),n,p, c=c, cc=cc,ne=ne)
        !call PropPrepare('PertFock',D(1),p, n,pp,cc, ci=c)
        if (size(D)/=product(cc(n+1:n+n)-cc(1:n)+2)) &
            call quit('PertFock: size(D) should be product 1+rank p')
        b = cc(n+1:n+n)-cc(1:n)+1
        pb = product(b); pb1 = product(b+1)
        !zero the Fock matrices
        do i=1,pb; F(i)=0.0d0*D(1); enddo
        if (n==0) then
            call quit('PertFock: unperturbed Fock matrix requested')
        elseif (all(p=='EL'.or.p=='MAGO')) then
            !lowest order contribution
            if (n==1.or.n==2.and.all(p=='MAGO')) &
                call prop_fock(p,D(1:1),F,c=c)
            !highest order contribution
            call prop_fock(p(1:0),D,F,b=b)
        elseif (n==1) then
            call prop_fock(p,(/D(1)/),F)
            call prop_fock(p(1:0),(/D(1),D(pb1-pb+1:pb1)/),F)
        elseif (n==2) then
            call prop_fock(p,(/D(1)/),F)
            call prop_fock(p(2:2),(/D(1),D(1+1:1+b(1))/),F)
            call prop_fock(p(1:1),(/D(1),D(1+b(1)+1:1+b(1)+b(2))/),F,o=(/2/))
            call prop_fock(p(1:0),(/D(1),D(1+1:1+b(1)),D(1+b(1)+1:1+b(1)+b(2)) &
                          ,D(pb1-pb+1:pb1)/),F,b=(/b(1),b(2)/)) !only works for HF (now also DFT)
        else !ajt: this needs to be done, at some point
            call quit('PertFock: General case not implemented')
        endif
    end subroutine


    subroutine pert_scf_eq(p,S0,D,F, DSD,FDS, c,w)
    !p'th derivative of idempotency DSD-D and SCF equation FDS-SDF
        character(*),intent(in):: p(:) !perturbation lables p(n)
        type(matrix),intent(in):: S0,D(:) !overlap, un- and perturbed density matrices
        type(matrix),intent(in):: F(:) !perturbed Fock matrices
        type(matrix),intent(inout):: DSD(*),FDS(*) !(*) to permit re-ranking
        integer,intent(in),optional:: c(:) !range within each p, c(2*n) (lower,upper indices)
        real*8,intent(in),optional:: w(:) !signed frequencies of each p
        integer:: n,cc(2*size(p)),b(size(p)),op(size(p)),ne
        integer:: i,j,k,i0,j0,k0,ij,ik,jk,ij0,ik0,jk0,ijk,ijk0,ji,jj
        integer:: imax
        real*8:: ww(size(p)+1); type(matrix):: Sp(3),Tp(3),X,SDF(1); n=size(p)
        call prop_check('pert_scf_eq',S0,n,p, c=c,w=w, op=op,cc=cc,ww=ww,ne=ne); ww(n+1)=0
        if (any(op/=(/(i,i=1,n)/))) call quit('PertSCFeq: General case not implemented')
        b = cc(n+1:n+n)-cc(1:n)+1
        if (size(D)/=product(b+1)) call quit( &
            'pert_scf_eq: size(D)/=product(b+1) wrong number of density matrices')
        if (size(F)/=product(b+1)) call quit( &
            'pert_scf_eq: size(F)/=product(b+1) wrong number of Fock matrices')
        !ajt fixme: This should be written optimized an for general orders



        X = -sum(ww/2.0d0)*S0

!#define AOR_ONLY_INDIVIDUAL_COMPONENTS
#ifdef AOR_ONLY_INDIVIDUAL_COMPONENTS
        if (all(p == 'EL' .or. p == 'MAGO')) then

          DSD(1) = tiny(1.0d0)*D(1)
          FDS(1) = tiny(1.0d0)*F(1)
          SDF(1) = tiny(1.0d0)*F(1)

          imax = size(D)

          do i = 1, imax
            j = imax - i + 1

            DSD(1) = DSD(1) + D(i)*S0*D(j)
            FDS(1) = FDS(1) + F(i)*D(j)
            SDF(1) = SDF(1) + D(i)*F(j)
          end do

          DSD(1) = DSD(1) - D(imax)

          FDS(1) = FDS(1) + X*D(imax)
          SDF(1) = SDF(1) + D(imax)*X

          FDS(1) = FDS(1)*S0 - S0*SDF(1)
#else



        if      (n == 1 .and. all(p == 'EL' .or. p == 'MAGO')) then
          do i = 1, b(1)
            DSD(i) = D(1+i)*S0*D(1  ) &
                   + D(1  )*S0*D(1+i) &
                   - D(1+i)
            FDS(i) = F(1+i)*D(1  )*S0 - S0*D(1  )*F(1+i) &
                   + F(1  )*D(1+i)*S0 - S0*D(1+i)*F(1  ) &
                   - ww(1)*S0*D(1+i)*S0
          end do
        else if (n == 2 .and. all(p == 'EL' .or. p == 'MAGO')) then
          do i0 = 1, b(1)
            do j0 = 1, b(2)
                ij0 = i0 + b(1)*(j0-1)
                i   = 1 + i0
                j   = 1 + b(1) + j0
                ij  = 1 + b(1) + b(2) + ij0

                DSD(ij0) = D(ij)*S0*D(1)+D(1)*S0*D(ij)-D(ij)+D(j)*S0*D(i)+D(i)*S0*D(j)
                FDS(ij0) = (F(ij)*D(1)+F(j)*D(i)+F(i)*D(j)+(F(1)-sum(ww/2)*S0)*D(ij))*S0 &
                      - S0*(D(1)*F(ij)+D(i)*F(j)+D(j)*F(i)+D(ij)*(F(1)+sum(ww/2)*S0))

             !  DSD(ij0) = D(ij)*S0*D(1) &
             !           + D(1 )*S0*D(ij) &
             !           + D(j )*S0*D(i) &
             !           + D(i )*S0*D(j) &
             !           - D(ij)
             !  FDS(ij0) = (F(ij)*D(1) &
             !           +  F(j )*D(i) &
             !           +  F(i )*D(j) &
             !           +  F(1) *D(ij))*S0 &
             !           -  sum(ww/2)*S0*D(ij)*S0 &
             !           - S0*(D(1 )*F(ij) &
             !           +     D(i )*F(j ) &
             !           +     D(j )*F(i ) &
             !           +     D(ij)*F(1)) &
             !           - S0*sum(ww/2)*S0
            end do
          end do
        else if (n == 3 .and. all(p == 'EL' .or. p == 'MAGO')) then
          do i0 = 0, b(1) - 1
          do j0 = 0, b(2) - 1
          do k0 = 0, b(3) - 1

            i    = i0 + 2
            j    = j0 + 2 + b(1)
            k    = k0 + 2 + b(1) + b(2)
                                     
            ij   = i0 + 2 + b(1) + b(2) + b(3) + b(1)*j0
            ik   = i0 + 2 + b(1) + b(2) + b(3) + b(1)*b(2) + b(1)*k0
            jk   = j0 + 2 + b(1) + b(2) + b(3) + b(1)*b(2) + b(1)*b(3) + b(2)*k0
                                            
            ijk  = i0 + 2 + b(1) + b(2) + b(3) + b(1)*b(2) + b(1)*b(3) + b(2)*b(3) + b(1)*j0 + b(1)*b(2)*k0
                                            
            ijk0 =  i0 + 1 +b(1)*j0 + b(1)*b(2)*k0
                                            
                DSD(ijk0) = D(ijk)*S0*D(1)*S0*D(ijk)-D(ijk)+D(jk)*S0*D(i)+D(i)*S0*D(jk) &
                          + D(ik)*S0*D(j)+D(j)*S0*D(ik)+D(ij)*S0*D(k)+D(k)*S0*D(ij)
                FDS(ijk0) = (F(ijk)*D(1)+F(jk)*D(i)+F(ik)*D(j)+F(ij)*D(k) &
                          +  F(k)*D(ij)+F(j)*D(ik)+F(i)*D(jk)+(F(1)-sum(ww/2)*S0)*D(ijk))*S0 &
                       - S0*(D(1)*F(ijk)+D(i)*F(jk)+D(j)*F(ik)+D(k)*F(ij) &
                          +  D(ij)*F(k)+D(ik)*F(j)+D(jk)*F(i)+D(ijk)*(F(1)+sum(ww/2)*S0))



          end do
          end do
          end do
#endif /* AOR_ONLY_INDIVIDUAL_COMPONENTS */
        elseif (n==1) then
            do i=1,b(1)
                call prop_overlap(p,S0,Sp,c=(/cc(1)+i-1/))
                call prop_overlap(p,S0,Tp,c=(/cc(1)+i-1/),w=(/-ww(1)/2,0d0/))
                DSD(i) = D(1)*Sp(1)*D(1) + D(1+i)*S0*D(1)+D(1)*S0*D(1+i)-D(1+i)
                FDS(i) = (F(1+i)-Tp(1))*D(1)*S0 - S0*D(1)*(F(1+i)+Tp(1)) &
                       + (F(1)-(ww(1)/2)*S0)*(D(1+i)*S0+D(1)*Sp(1)) &
                       - (Sp(1)*D(1)+S0*D(1+i))*(F(1)+(ww(1)/2)*S0); Sp(1)=0; Tp(1)=0
            enddo
        elseif (n==2) then
            do j=0,b(2)-1; do i=1,b(1); jj=1+b(1)+j+1; ij=i+b(1)*j; ji=1+b(1)+b(2)+ij
                call prop_overlap(p(1:1),S0,Sp(1),c=(/cc(1)+i-1/))
                call prop_overlap(p(1:1),S0,Tp(1),c=(/cc(1)+i-1/),w=(/-ww(1)/2,0d0/))
                call prop_overlap(p(2:2),S0,Sp(2),c=(/cc(2)+j/))
                call prop_overlap(p(2:2),S0,Tp(2),c=(/cc(2)+j/),w=(/-ww(2)/2,0d0/))
                call prop_overlap(p(1:2),S0,Sp(3),c=(/cc(1)+i-1,cc(2)+j/))
                call prop_overlap(p(1:2),S0,Tp(3),c=(/cc(1)+i-1,cc(2)+j/),w=(/-ww(1)/2,-ww(2)/2,0d0/))
                DSD(ij) = D(1)*Sp(3)*D(1)+D(1+i)*Sp(2)*D(1)+D(1)*Sp(2)*D(1+i) &
                        + (D(1+i)*S0+D(1)*Sp(1))*D(jj)+D(jj)*(S0*D(1+i)+Sp(1)*D(1))
                if (.not.iszero(D(ji))) DSD(ij) = DSD(ij)+D(ji)*S0*D(1)+D(1)*S0*D(ji)-D(ji)
                FDS(ij) = (F(ji)+Tp(3)-(ww(1)+ww(2))/2*Sp(3))*D(1)*S0+F(1)*(D(1)*Sp(3)) &
                        - S0*D(1)*(F(ji)+Tp(3)+(ww(1)+ww(2))/2*Sp(3))-(Sp(3)*D(1))*F(1) &
                        + (F(jj)+Tp(2)-ww(2)/2*Sp(2))*(D(1+i)*S0+D(1)*Sp(1)) &
                        - (S0*D(1+i)+Sp(1)*D(1))*(F(jj)+Tp(2)+ww(2)/2*Sp(2)) &
                        + (F(1+i)+Tp(1)-ww(1)/2*Sp(1))*D(jj)*S0-S0*D(jj)*(F(1+i)+Tp(1)+ww(1)/2*Sp(1))
                if (.not.iszero(D(ji))) FDS(ij) = FDS(ij) + (F(1)-(ww(1)+ww(2))/2*S0)*D(ji)*S0 &
                                                          - S0*D(ji)*(F(1)+(ww(1)+ww(2))/2*S0)
                if (.not.iszero(Sp(2))) FDS(ij) = FDS(ij) &
                        + ((F(1+i)+Tp(1)-ww(1)/2*Sp(1))*D(1)+(F(1)-ww(1)*S0)*D(1+i))*Sp(2) &
                        - Sp(2)*(D(1)*(F(1+i)+Tp(1)+ww(1)/2*Sp(1))+D(1+i)*(F(1)+ww(1)*S0))
                if (.not.iszero(Sp(1))) FDS(ij) = FDS(ij) &
                        + (F(1)-ww(2)*S0)*D(jj)*Sp(1)-Sp(1)*D(jj)*(F(1)+ww(2)*S0)
                Sp = 0; Tp = 0
            enddo; enddo
        else
            stop 'PertSCFeq: General case not implemented yet'
        endif
    end subroutine


    subroutine pert_dens(p,S0,D,F, DD,FF, c,w)
    !p'th derivative of idempotency DSD-D and SCF equation FDS-SDF
        character(*),intent(in):: p(:) !perturbation lables p(n)
        type(matrix),intent(in):: S0,D(:),F(:) !unperturbed overlap matrix,
                                  !lower-order density and Fock matrices
        type(matrix),intent(inout):: DD(*),FF(*) !highest order perturbed density and Fock
                                  !(*) is to permit variable ranks - size cannot be checked
        integer,intent(in),optional:: c(:) !range within each p, c(2*n) (lower,upper indices)
        real*8 ,intent(in),optional:: w(:) !signed frequencies of each p
        integer:: n,cc(2*size(p)),b(size(p)),op(size(p)),i,j,pb,pb1,ne
        real*8:: ww(size(p)+1); n=size(p)
        call prop_check('pert_dens',S0,n,p, c=c,w=w, op=op,cc=cc,ww=ww,ne=ne); ww(n+1)=0
        b = cc(n+1:n+n)-cc(1:n)+1
        pb = product(b); pb1 = product(b+1)
        if (size(D)/=pb1-pb) &
            call quit('pert_dens: D must have size product b+1 minus product b')
        if (size(F)/=pb1-pb) &
            call quit('pert_dens: F must have size product b+1 minus product b')
        !if (size(DD)/=pb) call QUIT('PertDensFock: DD must have size product b')
        !if (size(DD)/=pb) call QUIT('PertDensFock: FF must have size product b')
        call subr()
    contains
        subroutine subr() !called once pb, size of DSD and FDS is known
            type(matrix):: DSD(pb),FDS(pb)
            do i=1,pb; DD(i)=0d0*S0; FF(i)=0d0*S0; enddo
            if (all(p=='EL'.or.p=='MAGO').or.n==1.or.n==2) then
                call pert_scf_eq(p,S0,(/D,DD(1:pb)/),(/F,FF(1:pb)/),DSD,FDS,c=c,w=ww(1:n))
                call pert_fock(p,(/D,DD(1:pb)/),FF(1:pb),c=c)
                do i=1,pb
                    call solve_scf_eq(sum(ww(1:n)),S0,D(1),F(1),DSD(i),FDS(i),DD(i),FF(i))
                enddo
                print* !ajt newline after prints from SolveScfEq
            else
                stop 'pert_dens: General case not not implemented yet'
            endif
            !do i=1,product(b); DSD(i)=0; FDS(i)=0; enddo
        end subroutine
    end subroutine


    subroutine solve_scf_eq(w, S, D, F, DSD, FDS, Dp, Fp)
    !          in:          --------------------
    !          out:                               ------
    !ajt: should at some point also return the final residuals,
    !     at least FDS, for use in Sellers-type formulas

#ifdef PRG_DIRAC
        use dirac_interface
#else
        use newsolver
#endif

        real*8                      :: w
        type(matrix)                :: S, D, F
        type(matrix), intent(inout) :: DSD, FDS, Dp, Fp
        !----------------------------------------------------------------------------
        real*8                      :: x, ww(1)

        type(matrix)                :: DS, RHS(1), SOL(1)
#ifdef PRG_DIRAC               
        type(matrix)                :: R
#endif

        character*4                 :: nop(0)
        logical, save               :: first = .true.

#ifndef PRG_DIRAC
        if (first) call rsp2_init(1, 1, 1, 1, 1)
        first = .false.
#endif

        DS = D*S
        Dp = DSD - DSD*dag(DS) - DS*DSD
        RHS(1) = FDS*DS - dag(DS)*FDS + Fp
        if (debug_me) then
          call print_mat(Dp,     label = 'debug solve_scf_eq: Dp    ')
          call print_mat(Fp,     label = 'debug solve_scf_eq: Fp    ')
          call print_mat(DSD,    label = 'debug solve_scf_eq: DSD   ')
          call print_mat(FDS,    label = 'debug solve_scf_eq: FDS   ')
          call print_mat(RHS(1), label = 'debug solve_scf_eq: RHS(1)')
        end if
        FDS = 0

        call prop_fock(nop, (/D, Dp/), RHS)

        RHS(1) = -(RHS(1)*DS - dag(DS)*RHS(1)) !-(..) because RHS
        SOL(1) = tiny(1d0)*RHS(1)

#ifndef PRG_DIRAC
        x = norm(RHS(1))
        print *,'before response solver: norm(RHS) = ', x

        if (x > 1.0d-10) then
            ww(1) = -w !rsp_solver2 seems to solve for negative frequency
            call rsp_solver2(D, S, F, .true., 1, RHS, ww, SOL)
        else
            print *,'=> skipping this equation'
        end if

        DSD = SOL(1)*dag(DS) - DS*SOL(1)
        Dp  = Dp + 2.0d0*DSD
#else
!       h+ part
        R = -0.5d0*(RHS(1) - dag(RHS(1)))
        if (debug_me) then
          call print_mat(R, label = 'debug solve_scf_eq: R h+')
        end if
        if (real(dot(R, R)) > aor_cfg_threshold_rhs) then
            R%tr_sym = 1
            call response_solver(w, R, SOL(1))
            Dp = Dp + SOL(1)
        end if

!       h- part
        R = -0.5d0*(RHS(1) + dag(RHS(1)))
        if (debug_me) then
          call print_mat(R, label = 'debug solve_scf_eq: R h-')
        end if
        if (real(dot(R, R)) > aor_cfg_threshold_rhs) then
            R%tr_sym = -1
            call response_solver(w, R, SOL(1))
            Dp = Dp + SOL(1)
        end if

        R = 0
#endif /* ifndef PRG_DIRAC */

        DSD    = 0
        DS     = 0
        SOL(1) = 0

        RHS(1) = Fp
        call prop_fock(nop, (/D, Dp/), RHS)
        Fp = RHS(1)

        RHS(1) = 0

    end subroutine

end module
