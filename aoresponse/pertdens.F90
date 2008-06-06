module PertDens

    use matrix_mod
#ifdef PRG_DIRAC
    use aoresponse_input_m
    use dirac_interface_m
#endif /* PRG_DIRAC */
    use PropInt
    implicit none

    public PertFock,PertSCFeq,PertDensFock
    private

contains


    subroutine PertFock(p,D, F, c)
        character(*),intent(in):: p(:) !perturbation lables p(n)
        type(Matrix),intent(in):: D(:) !un- and perturbed density matrices
        type(Matrix),intent(inout):: F(*) !output, perturbed Fock matrices
                                          !(*) to permit re-ranking
        integer,optional:: c(:) !range within each p, c(2*n) (lower,upper indices)
        integer:: n,cc(2*size(p)),b(size(p)),i,pb,pb1,ne; n=size(p)
        call PropCheck('PertFock',D(1),n,p, c=c, cc=cc,ne=ne)
        !call PropPrepare('PertFock',D(1),p, n,pp,cc, ci=c)
        if (size(D)/=product(cc(n+1:n+n)-cc(1:n)+2)) &
            call QUIT('PertFock: size(D) should be product 1+rank p')
        b = cc(n+1:n+n)-cc(1:n)+1
        pb = product(b); pb1 = product(b+1)
        !zero the Fock matrices
        do i=1,pb; F(i)=0.0d0*D(1); enddo
        if (n==0) then
            call QUIT('PertFock: unperturbed Fock matrix requested')
        elseif (all(p=='EL'.or.p=='MAGO')) then
            !lowest order contribution
            if (n==1.or.n==2.and.all(p=='MAGO')) &
                call PropFock(p,D(1:1),F,c=c)
            !highest order contribution
            call PropFock(p(1:0),(/D(1),D(pb1-pb+1:pb1)/),F)
        else !ajt: this needs to be done, at some point
            call QUIT('PertFock: General case not implemented')
        endif
    end subroutine


    subroutine PertSCFeq(p,S0,D,F, DSD,FDS, c,w)
    !p'th derivative of idempotency DSD-D and SCF equation FDS-SDF
        character(*),intent(in):: p(:) !perturbation lables p(n)
        type(Matrix),intent(in):: S0,D(:) !overlap, un- and perturbed density matrices
        type(Matrix),intent(in):: F(:) !perturbed Fock matrices
        type(Matrix),intent(inout):: DSD(*),FDS(*) !(*) to permit re-ranking
        integer,intent(in),optional:: c(:) !range within each p, c(2*n) (lower,upper indices)
        real*8,intent(in),optional:: w(:) !signed frequencies of each p
        integer:: n,cc(2*size(p)),b(size(p)),op(size(p)),ne
        integer:: i,j,k,i0,j0,k0,ij,ik,jk,ij0,ik0,jk0,ijk,ijk0
        real*8:: ww(size(p)+1); n=size(p)
        call PropCheck('PertSCFeq',S0,n,p, c=c,w=w, op=op,cc=cc,ww=ww,ne=ne); ww(n+1)=0
        if (any(op/=(/(i,i=1,n)/))) call QUIT('PertSCFeq: General case not implemented')
        b = cc(n+1:n+n)-cc(1:n)+1
        if (size(D)/=product(b+1)) stop &
            'PertSCFeq: size(D)/=product(b+1) wrong number of density matrices'
        if (size(F)/=product(b+1)) stop &
            'PertSCFeq: size(F)/=product(b+1) wrong number of Fock matrices'
        if (n==1.and.all(p=='EL'.or.p=='MAGO')) then
            do i=1,b(1)
                DSD(i) = D(1+i)*S0*D(1)+D(1)*S0*D(1+i)-D(1+i)
                FDS(i) = F(1+i)*D(1)*S0-S0*D(1)*F(1+i) &
                       + F(1)*D(1+i)*S0-S0*D(1+i)*F(1)+ww(1)*S0*D(1+i)*S0
            enddo
        elseif (n==2.and.all(p=='EL'.or.p=='MAGO')) then
            do i0=1,b(1); do j0=1,b(2); ij0=i0+b(1)*(j0-1)
                i=1+i0; j=1+b(1)+j0; ij=1+b(1)+b(2)+ij0
                DSD(ij0) = D(ij)*S0*D(1)+D(1)*S0*D(ij)-D(ij)+D(j)*S0*D(i)+D(i)*S0*D(j)
                FDS(ij0) = (F(ij)*D(1)+F(j)*D(i)+F(i)*D(j)+(F(1)+sum(ww/2)*S0)*D(ij))*S0 &
                      - S0*(D(1)*F(ij)+D(i)*F(j)+D(j)*F(i)+D(ij)*(F(1)-sum(ww/2)*S0))
            enddo; enddo
        elseif (n==3.and.all(p=='EL'.or.p=='MAGO')) then
            do i0=1,b(1); do j0=1,b(2); do k0=1,b(3)
                ij0=i0+b(1)*(j0-1); ik0=i0+b(1)*(k0-1); jk0=j0+b(2)*(k0-1)
                ijk0=i0+b(1)*(j0-1)+b(1)*b(2)*(k0-1)
                i=1+i0; j=1+b(1)+j0; k=1+b(1)+b(2)+k0; ij=1+b(1)+b(2)+b(3)+ij0
                ik=1+b(1)+b(2)+b(3)+b(1)*b(2)+ik0; jk=1+b(1)+b(2)+b(3)+b(1)*b(2)+b(1)*b(3)+jk0
                ijk=1+b(1)+b(2)+b(3)+b(1)*b(2)+b(1)*b(3)+b(2)*b(3)+ijk0
                DSD(ijk0) = D(ijk)*S0*D(1)*S0*D(ijk)-D(ijk)+D(jk)*S0*D(i)+D(i)*S0*D(jk) &
                          + D(ik)*S0*D(j)+D(j)*S0*D(ik)+D(ij)*S0*D(k)+D(k)*S0*D(ij)
                FDS(ijk0) = (F(ijk)*D(1)+F(jk)*D(i)+F(ik)*D(j)+F(ij)*D(k) &
                          +  F(k)*D(ij)+F(j)*D(ik)+F(i)*D(jk)+(F(1)+sum(ww/2)*S0)*D(ijk))*S0 &
                       - S0*(D(1)*F(ijk)+D(i)*F(jk)+D(j)*F(ik)+D(k)*F(ij) &
                          +  D(ij)*F(k)+D(ik)*F(j)+D(jk)*F(i)+D(ijk)*(F(1)-sum(ww/2)*S0))
            enddo; enddo; enddo
        else
            stop 'PertSCFeq: General case not implemented yet'
        endif
    end subroutine


    subroutine PertDensFock(p,S0,D,F, DD,FF, c,w)
    !p'th derivative of idempotency DSD-D and SCF equation FDS-SDF
        character(*),intent(in):: p(:) !perturbation lables p(n)
        type(Matrix),intent(in):: S0,D(:),F(:) !unperturbed overlap matrix,
                                  !lower-order density and Fock matrices
        type(Matrix),intent(inout):: DD(*),FF(*) !highest order perturbed density and Fock
                                  !(*) is to permit variable ranks - size cannot be checked
        integer,intent(in),optional:: c(:) !range within each p, c(2*n) (lower,upper indices)
        real*8 ,intent(in),optional:: w(:) !signed frequencies of each p
        integer:: n,cc(2*size(p)),b(size(p)),op(size(p)),i,j,pb,pb1,ne
        real*8:: ww(size(p)+1); n=size(p)
        call PropCheck('PertDensFock',S0,n,p, c=c,w=w, op=op,cc=cc,ww=ww,ne=ne); ww(n+1)=0
        b = cc(n+1:n+n)-cc(1:n)+1
        pb = product(b); pb1 = product(b+1)
        if (size(D)/=pb1-pb) &
            call QUIT('PertDensFock: D must have size product b+1 minus product b')
        if (size(F)/=pb1-pb) &
            call QUIT('PertDensFock: F must have size product b+1 minus product b')
        !if (size(DD)/=pb) call QUIT('PertDensFock: DD must have size product b')
        !if (size(DD)/=pb) call QUIT('PertDensFock: FF must have size product b')
        call subr()
    contains
        subroutine subr() !called once pb, size of DSD and FDS is known
            type(Matrix):: DSD(pb),FDS(pb)
            if (all(p=='EL'.or.p=='MAGO')) then
                do i=1,pb; DD(i)=0.0d0*S0; FF(i)=0.0d0*S0; enddo
                call PertSCFeq(p,S0,(/D,DD(1:pb)/),(/F,FF(1:pb)/),DSD,FDS,c=c,w=w)
                call PertFock(p,(/D,DD(1:pb)/),FF(1:pb),c=c)
                do i=1,pb
                    call SolveSCFeq(sum(w),S0,D(1),F(1),DSD(i),FDS(i),DD(i),FF(i))
                enddo
            else
                stop 'PertDensFock: General case not not implemented yet'
            endif
            !do i=1,product(b); DSD(i)=0; FDS(i)=0; enddo
        end subroutine
    end subroutine

    subroutine SolveSCFeq(w,S0,D0,F0,DSD,FDS,D,F)
#ifndef PRG_DIRAC
        use NEWSOLVER
#endif /* PRG_DIRAC */
        real*8:: w; real*8:: ww(1)
        type(Matrix):: S0,D0,F0
        type(Matrix),intent(inout):: DSD,FDS,D,F
        type(Matrix):: DS,RHS(1),SOL(1)
#ifdef PRG_DIRAC
        type(matrix) :: T, R
#endif /* PRG_DIRAC */
        character*4:: NOP(0)
        integer:: NOI(0)
        logical,save:: first=.true.
#ifndef PRG_DIRAC
        if (first) call rsp2_init(1,1,1,1,1)
#endif /* PRG_DIRAC */
        first = .false.
        DS = D0*S0
        D = DSD - DSD*trps(DS) - DS*DSD
        RHS(1) = FDS*DS - trps(DS)*FDS + F; FDS=0
        call PropFock(NOP,(/D0,D/),RHS)
        RHS(1) = RHS(1)*DS - trps(DS)*RHS(1)
        ww(1) = w; SOL(1)=DSD; DSD=0
        if (iszero(SOL(1))) SOL(1) = 1d-32*S0 !trick to get it allocated
#ifndef PRG_DIRAC
        call rsp_solver2(D0,S0,F0,.true.,1,RHS,ww,SOL)
        DSD=SOL(1); RHS(1)=0; SOL(1)=0
        DSD = DSD*trps(DS) - DS*DSD; DS=0
        D = D + (-2)*DSD; DSD=0 !notice the -2
#else
        T = 1.0d-100*RHS(1)

        R = (RHS(1) - dag(RHS(1)))/2.0d0
        if (R**R > aoresponse_threshold_rhs) then
          R%tr_sym = 1
          call response_solver(ww, R, T)
          D = D + T
        end if

        R = (RHS(1) + dag(RHS(1)))/2.0d0
        if (R**R > aoresponse_threshold_rhs) then
          R%tr_sym = -1
          call response_solver(ww, R, T)
          D = D + T
        end if

        R = 0
        T = 0
#endif /* PRG_DIRAC */
        RHS(1)=F; F=0
        call PropFock(NOP,(/D0,D/),RHS); F=RHS(1); RHS(1)=0
    end subroutine


end module PertDens
