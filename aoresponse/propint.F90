module PropInt
!routines for property contractions and matrices

#ifndef PRG_DIRAC
    use TypeMatrix
    use matrix_operations,only: mat_init,mat_to_full,mat_set_from_full
    use dal_interface,only: di_GET_GbDs
#else
    use matrix_mod
    use dirac_interface_m
#endif /* PRG_DIRAC */

    implicit none
    public PropEnergy,PropReort,PropFock,PropOverlap
    public PropImag,PropRank,PropCheck
    private

    type prt !private: 'perturbation'
        character*4 :: code !four-letter abbreviation
        character*64:: name !long name
        integer     :: rank !number of components (when known, 0 otw.)
        logical:: imag !imaginary (1st ord.) perturbed integrals
        logical:: bas  !basis dependent (sa. GEO and MF)
        logical:: lin  !one-electron operator linear in field strength (EL)
        logical:: quad !one-electron operator quadratic in field strength (MAGO)
    end type

    !starting with two fields
    logical,parameter:: T=.true.,F=.false.
    type(prt):: tab(2) = (/ &                            !rk im ba ln qu
        prt('EL'  ,'Electric field'                      , 3, F, F, T, F), &
        prt('MAGO','Magnetic field w/o London orbitals'  , 3, T, F, F, T)/)

    integer:: nbas=0,natom=0 !unitialized
    integer,parameter:: lwork=6000000 !totally wrong, should be dynamic

    character*1,parameter:: xyz(3) = (/'X','Y','Z'/), &
            c09(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)

contains


    subroutine PropEnergy(p,D, e, c,b,o)
    !Front for 'erg' below, filling in optionals, checking sizes, and making permutations
    !NB: The result of the contraction is ADDED to 'e' - remember to zero before first call
        character(*),intent(in):: p(:) !p(n) perturbation lables
        type(Matrix),intent(in):: D(:) !density matrix expansion to contract integrals with
        real*8,intent(inout):: e(*) !result is added to 'e'. (*) to permit re-ranking
                                    !size(e) is not checked
        integer,intent(in),optional:: &
                 c(:), & !c(n) component indices, or c(n+n) components ranges of each p
                 b(:), & !b(m) ranks of D(:) expansion: size(D)=product(b+1)
                 o(:) !o(m) indices in e of each rank in D; default 1,2..m
        integer:: n,m,ne,op(size(p)),cc(2*size(p)),i
        logical:: zerg,imag; n=size(p)
        !process and check arguments
        call PropCheck('PropEnergy',D(1),n,p, c=c,nd=size(D),b=b,o=o &
                      ,op=op,cc=cc,m=m,ne=ne,zerg=zerg,imag=imag)
        if (zerg.or.m==0.and.imag) return !the result is zero
        call subr() !since bb(m) and oe(ne) and ee(ne)
    contains
        subroutine subr()
            integer:: bb(m),oe(ne); real*8:: ee(ne); ee=0
            call PropOrder(n,op,cc,m,size(D),ne, bb,oe, b=b,o=o)
            !call core function erg
            call erg(n,m,(/(p(op(i)),i=1,n)/) &
                    ,(/(cc(op(i)),i=1,n),(cc(n+op(i)),i=1,n)/),D,bb,ee)
            !add the reordered output to e
            e(1:ne) = e(1:ne) + (/(ee(oe(i)),i=1,ne)/)
        end subroutine
    end subroutine


    subroutine PropFock(p,D, F, c,b,o)
    !front for 'fck' below, checking dimensions and taking care permutations (not yet, actually)
    !Note that the result of the contraction is ADDED to F, ie. remeber to zero F before calling
        character(*),intent(in):: p(:) !p(n) perturbation lables
        type(Matrix),intent(in):: D(:) !density matrix expansion to contract integrals with
        type(Matrix):: F(*) !result is added to F. (*) to permit re-ranking, so size not checked
        integer,intent(in),optional:: &
                 c(:), & !c(n) component indices, or c(n+n) components ranges of each p
                 b(:), & !b(m) ranks of D(:) expansion: size(D)=product(b+1)
                 o(:)    !o(m) indices in e of each rank in D; default 1,2..m
        integer:: n,m,nf,op(size(p)),cc(2*size(p)),i
        logical:: zfck; n=size(p)
        call PropCheck('PropFock',D(1),n,p, c=c,nd=size(D),b=b,o=o &
                      ,op=op,cc=cc,m=m,ne=nf,zerg=zfck)
        if (zfck) return !the result is zero
        call subr() !since bb(m) and oe(ne) and FF(ne)
    contains
        subroutine subr()
            integer:: bb(m),of(nf); type(Matrix):: FF(nf)
            call PropOrder(n,op,cc,m,size(D),nf, bb,of, b=b,o=o)
            !order F into FF, call core function fck, then reorder back into F
            do i=1,nf; FF(of(i)) = F(i); F(i)=0; enddo
            call fck(n,m,(/(p(op(i)),i=1,n)/) &
                    ,(/(cc(op(i)),i=1,n),(cc(n+op(i)),i=1,n)/),D,bb,FF)
            do i=1,nf; F(i) = FF(of(i)); FF(of(i))=0; enddo
        end subroutine
    end subroutine


    subroutine PropOverlap(p,S0, S, c,w)
    !front for 'ovl' below, checking dimensions. Argument S0 is for context only
    !S will on exit contain the array of perturbed overlap matrices (which may all be zero (0*S0))
        character(*),intent(in):: p(:) !p(n) perturbation lables
        type(Matrix),intent(in):: S0 !overlap matrix, for context (shape, basis set etc.)
        type(Matrix):: S(*) !resulting perturbed overlap matrices.
                            !intent(inout) is dropped to permit re-ranking: (/S(:,:)/)
        integer,intent(in),optional:: &
                 c(:) !c(n) component indices, or c(n+n) components ranges of each p
        real*8,intent(in),optional:: w(:)    !w(n+1) frequences for each p
        integer:: n,op(size(p)),cc(2*size(p)),ns,i
        real*8:: ww(size(p)+1); logical:: zovl; n=size(p)
        !process and check arguments
        call PropCheck('PropOverlap',S0,n,p, c=c,w=w, op=op,cc=cc,ww=ww,ne=ns,zovl=zovl)
        do i=1,ns; S(i) = 0.0d0*S0; enddo
        if (zovl) return !the result is zero
        call subr() !since bb(m) and oe(ne) and FF(ne)
    contains
        subroutine subr()
            integer:: os(ns); type(Matrix):: SS(ns)
            call PropOrder(n,op,cc,0,0,ns, cc(1:0),os)
            !call core function ovl, then reorder back into S
            call ovl(n,(/(p(op(i)),i=1,n)/),(/(cc(op(i)),i=1,n),(cc(n+op(i)),i=1,n)/) &
                    ,(/(ww(op(i)),i=1,n),w(n+1)/),S0,SS)
            do i=1,ns; S(i) = SS(os(i)); SS(os(i))=0; enddo
        end subroutine
    end subroutine


    subroutine PropReort(p,c,S0,D,DFD,b,e, w)
    !front for 'reort' below, checking dimensions and taking care of permutations (not yet, actually)
        character(*):: p(:) !p(n) perturbation lables
        integer:: n,c(:) !c(n*2) component range for each p
        integer:: m,b(:) !b(m) ranks of D and DFD
        type(Matrix):: S0,D(product(b)) !S0 uneperturbed overlap and D(product(b)) perturbed density matrices
        type(Matrix):: DFD(product(b)) !energy-weighted density matrices
        real*8,intent(inout):: e(product(b)*product(c(size(p)+1:)-c(:size(p)))) !result is added to e
        real*8,optional:: w(:) !frequencies w(n+1)
        integer:: i,j,k,ip(size(p))
        real*8:: ww(size(p)+1)
        n = size(p); m = size(b)
        if (size(c)/=n*2) &
            stop 'PropReort: size(c)/=2*size(p)'
        if (present(w)) then
            if (size(w)/=n+1) stop 'PropReort: wrong number of frequencies, size(w)/=size(p)+1'
            ww = w
        endif
        call init(S0)
    end subroutine


    subroutine erg(n,m,p,c,D,b,e)
        integer:: n,m,c(2*n),b(m); character(*):: p(n) !perturbation lables
        type(Matrix):: D(product(b+1)) !density matrix expansion to contract over
        real*8,intent(out):: e(product(b)*product(c(n+1:n+n)-c(1:n)+1)) !result
        integer:: i,j,k,ii,jj,kk,pb,pb1
        type(Matrix):: A(3); real*8:: r(3*3*natom)
        pb = product(b); pb1 = product(b+1)
        if (n==0) then
            if (m==0) stop 'PropEnergy: unperturbed energy requested'
            !The highest order Ds multiply the unperturbed F, which is unknown here,
            !so insist that no highest-order Ds are present
            if (.not.all((/(iszero(D(i)),i=pb1-pb+1,pb1)/))) &
                stop 'PropEnergy: unperturbed Fock matrix needed'
            if (m==2) then !stopping at m=2, for now
                call mat_init(A(1),nbas,nbas)
                do j=0,b(2)-1
                    call di_GET_GbDs(D(1+b(1)+1+j),A(1))
                    do i=1,b(1)
                        e(i+b(1)*j) = trab(A(1),D(1+i)) !factor 2 hidden in trab
                    enddo
                enddo
                A(1) = 0
            elseif (m>2) then
                stop 'PropEnergy: m>2 not implemented at the moment'
            endif
        elseif (n==1.and.p(1)=='EL') then
            do j=0,c(2)-c(1); jj=c(1)+j
                A(1) = oneint(xyz(jj)//'DIPLEN ')
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*j) = e(i+pb*j) + A(1)**D(ii) !factor 2 hidden in **
                enddo
            enddo
#ifndef PRG_DIRAC
            if (m==0) call DIPNUC_ifc(natom,r(1:3),(/(0d0,i=1,9*natom)/))
#else
            !fixme
#endif /* PRG_DIRAC */
            if (m==0) e(1:c(2)-c(1)+1) = e(1:c(2)-c(1)+1) - r(c(1):c(2))
        elseif (n==1.and.p(1)=='MAGO') then
            do j=0,c(2)-c(1); jj=c(1)+j
                A(1) = oneint(xyz(jj)//'ANGMOM ')
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*j) = trab(A(1),D(ii))/2 !half, factor 2 hidden in trab
                enddo
            enddo
        elseif (n==2.and.all(p==(/'MAGO','MAGO'/))) then
            do k=0,c(4)-c(2); kk=c(2)+k
                do j=0,c(3)-c(1); jj=c(1)+j
                    A(1) = oneint(xyz(min(jj,kk))//xyz(max(jj,kk))//'SUSCGO')
                    do i=1,pb; ii=pb1-pb+i
                        e(i+pb*j+pb*3*k) = trab(A(1),D(ii)) !factor 2 hidden in trab
                    enddo
                enddo
            enddo
        else
            print *,'PropEnergy: no integrals for these perturbations: ', &
                    (p(i)//' ',i=1,n); stop
        endif
    end subroutine

    subroutine fck(n,m,p,c,D,b,F)
        integer:: n,m,c(n+n),b(m); character(*):: p(n) !perturbation lables
        type(Matrix):: D(product(b+1)) !density matrix expansion to contract over
        type(Matrix),intent(inout):: F(product(b)*product(c(n+1:n+n)-c(1:n)+1)) !result
        integer:: i,j,k,ii,jj,kk,pb,pb1; type(Matrix):: A(9)
        pb = product(b); pb1 = product(b+1)
        if (n==0) then
            if (m==0) stop 'PropFock: unperturbed Fock matrix requested'
            call mat_init(A(1),nbas,nbas)
            do i=1,pb; ii=pb1-pb+i
                if (iszero(D(ii))) cycle
                call di_GET_GbDs(D(ii),A(1))
                if (iszero(F(i))) F(i)=tiny(1d0)*D(1) !to get it allocated
                F(i) = F(i) + A(1)
            enddo
            A(1) = 0
        elseif (n==1.and.p(1)=='EL') then !only m==0
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                A(1) = oneint(xyz(ii)//'DIPLEN ')
                F(i) = F(i) + A(1)
            enddo
        elseif (n==1.and.p(1)=='MAGO') then !only m==0
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                A(1) = oneint(xyz(ii)//'ANGMOM ')
                F(i) = F(i) + A(1)/2.0d0 !
            enddo
        elseif (n==2.and.all(p==(/'MAGO','MAGO'/))) then !only m==0
            do j=0,c(4)-c(2); jj=c(2)+j
                do i=1,c(3)-c(1)+1; ii=c(1)+i-1
                    A(1) = oneint(xyz(min(ii,jj))//xyz(max(ii,jj))//'SUSCGO')
                    F(i+3*j) = F(i+3*j) + A(1)
                enddo
            enddo
        else
            print *,'PropFock: no integrals for these perturbations:', &
                    (' '//p(i),i=1,n); stop
        endif
!        do i=1,9; if (isdef(A(i))) call mat_free(A(i)); enddo
    end subroutine


    subroutine ovl(n,p,c,w,S0,S)
    !frequency generalized perturbed overlap matrices
        integer:: n,c(n*2)
        character(*):: p(n) !perturbation lables
        real*8:: w(n+1) !frequencies
        type(Matrix):: S0 !unperturbed overlap matrix, for context
        type(Matrix),intent(inout):: S(product(c(n+1:2*n)-c(1:n))) !result
        integer:: i,j,k,ii,jj,kk
        if (n==0) then
            stop 'PropOverlap: unperturbed overlap matrix requested'
        else
            print *,'PropOverlap: no integrals for these perturbations:', &
                    (' '//p(i),i=1,n); stop
        endif
    end subroutine


    subroutine init(D0)
        type(Matrix):: D0; integer:: NATOMS_ifc
        if (D0%nrow==nbas) return !stub
        nbas = D0%nrow
#ifndef PRG_DIRAC
        natom = NATOMS_ifc()
#else
        natom = 137 ! fixme
#endif /* PRG_DIRAC */
!no GEO yet        tab(idx('GEO'))%rank = 3*natom
        !AllPrt(idx('NS'))%rank = depends on which isotopes
    end subroutine

    function oneint(l) result(A)
#ifndef PRG_DIRAC
        use dal_interface
#endif /* PRG_DIRAC */
        character*8:: l; type(Matrix):: A
        if (nbas==0) stop 'PropInt, oneint: nbas not initialized'
        if (.not. A%is_defined) call mat_init(A,nbas,nbas)
#ifndef PRG_DIRAC
        call di_read_operator_int(l,A); A=A !to get it tracked by module TypeMatrix
#else
        call get_P(l, A)
        A = A
#endif /* PRG_DIRAC */
    end function

    subroutine twomatr(o,a,D,F)
    !o=1,a=0 -> three magnetic
    !o=2,a=0 -> six magnetic
    !o=1,a=1..natom -> three geometric
        integer:: o,a; type(Matrix):: D; type(Matrix),intent(inout):: F(3**o)
        integer:: lw,i,j,k; real*8:: Df(nbas**2),Ff(nbas**2,3*o),w(lwork)
        if (nbas==0) stop 'twomatr: nbas not initialized'
#ifndef PRG_DIRAC
        call mat_to_full(D,1d0,Df); lw = lwork
        call GRCONT(w,lw,Ff,nbas**2*3*o,a/=0,a==0,o,a,.false.,.true.,Df,1)
#else
        call quit('error in twomatr')
#endif /* PRG_DIRAC */
        do i=1,3**o; if (.not.F(i)%is_defined) call mat_init(F(i),nbas,nbas); enddo
        do i=1,3**o
            j = i
            if (o==2) then
                j=1+(i-1)/3; k=1+mod(i-1,3)
                j=max(j,k)*(max(j,k)-1)/2+min(j,k)
            endif
#ifndef PRG_DIRAC
            call mat_set_from_full(Ff(:,j),1d0,F(i))
#else
            call quit('error in twomatr')
#endif /* PRG_DIRAC */
        enddo
    end subroutine

    subroutine twocont(gm,D,DD,e)
    !gm=T -> 3*NATOMS geometriske
    !gm=F -> tre magnetiske
        logical:: gm; type(Matrix):: D,DD
        real*8,intent(out):: e(3*merge(natom,1,gm))
        integer:: lw
        real*8:: A(nbas,nbas,2),w(lwork)
        if (nbas==0) stop 'twocont: nbas not initialized'
#ifndef PRG_DIRAC
        call mat_to_full(D ,1d0,A(:,:,1))
        call mat_to_full(DD,1d0,A(:,:,2))
        lw=lwork; e=0
        call GRCONT(w,lw,e,size(e),gm,.not.gm,1,0,.true.,.false.,A,2)
#else
        call quit('error in twocont')
#endif /* PRG_DIRAC */
    end subroutine

    function idx(p)
        character(*):: p; integer:: idx
        do idx=1,size(tab)
            if (tab(idx)%code==p) return
        enddo
        print *,'Perturbation not found: '//p; stop
    end function

    function PropImag(p)
        character(*):: p; logical:: PropImag
        PropImag = tab(idx(p))%imag
    end function

    function PropRank(p)
        character(*):: p; integer:: PropRank
        PropRank = tab(idx(p))%rank
        if (PropRank==0) stop 'PropRank: rank not initialized'
    end function

    subroutine PropCheck(who,D0,n,p, c,w,nd,b,o, op,cc,ww,m,ne,zerg,zovl,imag)
    !Process, p,c,w,nd,b finding the order p(op(k)) which sorts p by idx,
    !fill in from c the component ranges cc(1:n+n), and the frequencies ww(1:n+1).
    !Determine m from nd and b
        integer,intent(in):: n
        character(*),intent(in):: who,p(n)
        type(Matrix),intent(in):: D0 !(or S0) for reference
        integer,intent(in),optional:: c(:),nd,b(:),o(:)
        real*8,intent(in),optional:: w(:)
        integer,intent(out),optional:: op(n),cc(n+n),m,ne
        real*8,intent(out),optional:: ww(n+1)
        logical,intent(out),optional:: zerg,zovl,imag
        integer:: ip(n),i,j,k
        call init(D0)
        ip = (/(idx(p(i)),i=1,n)/)
        !find zerg, zovl and imag
        if (present(zerg)) then
            j = count((/(tab(ip(i))%lin ,i=1,n)/))
            k = count((/(tab(ip(i))%quad,i=1,n)/))
            zerg = (j>1.or.k>2.or.j==1.and.k>=1)
            if (j==0.and.k==2) then
                j=0; do i=1,n; if (tab(i)%quad) then
                    if (j/=0) zerg = (ip(i)/=ip(j))
                j=i; endif; enddo
            endif
        endif
        !check whether the perturbed overlap integrals are zero
        if (present(zovl)) &
            zovl = (count((/(tab(ip(i))%bas,i=1,n)/))==n)
        !the integrals are imaginary when
        if (present(imag)) & 
            imag = (mod(count((/(tab(ip(i))%imag,i=1,n)/)),2)==1)
        !process optional c to produce cc
        if (present(c)) then
            if (size(c)/=n.and.size(c)/=n+n) &
                call QUIT(who//': c must be c(n) or c(n+n), with p(n)')
            if (any((/(1>c(i).or.c(i)>tab(ip(i))%rank,i=1,n)/))) &
                call QUIT(who//': c must be 1 <= c(1:n) <= p%comp')
            if (present(cc)) cc = (/c(1:n),c(1:n)/)
            if (size(c)==n+n) then
                if (any((/(c(i)>c(n+i).or.c(n+i)>tab(ip(i))%rank,i=1,n)/))) &
                    call QUIT(who//': c must be c(1:n) <= c(n+1:n+n) <= p%comp')
                if (present(cc)) cc = c
            endif
            if (present(ne).and.size(c)==n  ) ne = 1
            if (present(ne).and.size(c)==n+n) ne = product(c(n+1:n+n)-c(1:n)+1)
        else !ranges 1:p%comp
            if (present(cc)) cc = (/(1,i=1,n),(tab(ip(i))%rank,i=1,n)/)
            if (present(ne)) ne = product((/(tab(ip(i))%rank,i=1,n)/))
        endif
        !process optional w, produce ww
        if (present(w)) then
            if (size(w)/=n.and.size(w)/=n+1) &
                call QUIT(who//': w must be w(n) or w(n+1), with p(n)')
            if (present(ww)) ww(1:n+1) = (/w(1:n),0d0/)
            if (present(ww).and.size(w)==n+1) ww(n+1) = w(n+1)
        elseif (present(ww)) then
            ww(1:n) = 0; ww(n+1) = 1
        endif
        !find order of p(i) according to idx
        if (present(op)) then
            do i=1,n
                op(i:i) = minloc(ip)
                ip(op(i)) = huge(1) !not found again by minloc
            enddo
        endif
        !if nd and b, check product(b+1)=nd, set m
        if (present(nd).and.present(m)) then
            if (present(b)) then
                if (product(b+1)/=nd) call QUIT(who// &
                        ': wrong number of density matrices: product(b+1)/=nd')
                m = size(b)
                if (present(ne)) ne = ne*product(b)
            else
                m = merge(0,1,nd==1)
                if (present(ne)) ne = ne*merge(1,nd-1,nd==1)
            endif
            if (present(o)) then
                if (size(o)/=m) call QUIT(who//': o has wrong size')
                if (any(1>o.or.o>m+n)) call QUIT(who//': o out of range')
                if (any((/((o(i)==o(j),i=j+1,m),j=1,m-1)/))) &
                    call QUIT(who//': repeated element in o')
            endif
        elseif (present(b).or.present(o)) then
            call QUIT(who//': if b or o present, so must nd and m be')
        endif
    end subroutine


    subroutine PropOrder(n,op,cc,m,nd,ne, bb,oe, b,o)
        integer,intent(in):: n,op(n),cc(n+n),m,nd,ne
        integer,intent(out):: bb(m),oe(ne)
        integer,intent(in),optional:: b(:),o(:) !b(m),o(m) will break ifort
        integer:: kk(m+n,3),k(3),opi(n),i,j,l
        !cumulative product of dimensions (pitches)
        bb=nd-1; if (present(b)) bb = b
        kk(1:m+n,1) = (/bb,(cc(n+op(i))-cc(op(i))+1,i=1,n)/)
        kk(1,2)=1; do i=2,m+n
            kk(i,2) = kk(i-1,2)*kk(i-1,1); enddo
        !reorder dimensions and pitches into output order
        do i=1,n; opi(op(i)) = i; enddo !inverse permutation
        kk(m+1:m+n,1) = (/(kk(m+opi(i),1),i=1,n)/)
        kk(m+1:m+n,2) = (/(kk(m+opi(i),2),i=1,n)/)
        if (present(o)) then
            kk(1:m,3) = o(1:m)
            do i=m,1,-1
                k(1:1) = minloc(kk(1:i,3))
                j = k(1); k = kk(j,:)
                kk(j:k(3)-1,:) = kk(j+1:k(3),:)
                kk(k(3),:) = k
            enddo
        endif
        !construct oe from dimensions kk(:,1) and 'pitch'es kk(:,2)
        oe(1) = 1; j = 1
        do i=1,m+n
            oe(j+1:j*kk(i,1)) = (/(oe(1:j)+l*kk(i,2),l=1,kk(i,1)-1)/)
            j = j*kk(i,1)
        enddo
    end subroutine


end module PropInt
