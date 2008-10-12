module PropInt
!routines for property contractions and matrices

#ifndef PRG_DIRAC
    use matrix_autop
    use dal_interface,only: di_GET_GbDs,di_read_operator_int,di_get_sigma_xc_cont,di_get_T_xc_cont &
                     ,di_get_MagDeriv_FxD_DFT,di_get_MagDeriv_GxD_DFT &
                     ,di_get_geomDeriv_molgrad_DFT,di_get_geomDeriv_FxD_DFT,di_get_geomDeriv_GxD_DFT
    use matrix_operations,only: mat_init,mat_to_full,mat_set_from_full
#else
    use matrix_autop
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
    type(prt):: tab(6) = &                                 !rk im ba ln qu
        (/prt('EL'  ,'Electric field'                      , 3, F, F, T, F) &
         ,prt('MAGO','Magnetic field w/o. London orbitals' , 3, T, F, F, T) &
         ,prt('MAG' ,'Magnetic field with London orbitals' , 3, T, T, F, F) &
         ,prt('ELGR','Electric field gradient'             , 6, F, F, T, F) &
         ,prt('VEL' ,'Velocity'                            , 3, T, F, T, F) &
         ,prt('GEO' ,'Nuclear coordinates'                 , 0, F, T, F, F)/) !0=uninitialized

    integer:: nbas=0,natom=0 !unitialized

    character*1,parameter:: xyz(3) = (/'X','Y','Z'/), &
            c09(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)

#ifndef PRG_DIRAC
    logical,external:: do_dft
#endif /* PRG_DIRAC */

contains


    subroutine PropEnergy(p,D, e, c,b,o)
    !Front for 'erg' below, filling in optionals, checking sizes, and making permutations
    !NB: The result of the contraction is ADDED to 'e' - remember to zero before first call
        character(*),intent(in):: p(:) !p(n) perturbation lables
        type(matrix),intent(in):: D(:) !density matrix expansion to contract integrals with
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
        type(matrix),intent(in):: D(:) !density matrix expansion to contract integrals with
        type(matrix):: F(*) !result is added to F. (*) to permit re-ranking, so size not checked
        integer,intent(in),optional:: &
                 c(:), & !c(n) component indices, or c(n+n) components ranges of each p
                 b(:), & !b(m) ranks of D(:) expansion: size(D)=product(b+1)
                 o(:)    !o(m) indices in e of each rank in D; default 1,2..m
        integer:: n,m,nf,op(size(p)),cc(2*size(p)),i
        logical:: zfck; n=size(p)
        call PropCheck('PropFock',D(1),n,p, c=c,nd=size(D),b=b,o=o &
                      ,op=op,cc=cc,m=m,ne=nf,zfck=zfck)
        if (zfck) return !the result is zero
        call subr() !since bb(m) and oe(ne) and FF(ne)
    contains
        subroutine subr()
            integer:: bb(m),of(nf); type(matrix):: FF(nf)
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
        type(matrix),intent(in):: S0 !overlap matrix, for context (shape, basis set etc.)
        type(matrix):: S(*) !resulting perturbed overlap matrices.
                            !intent(inout) is dropped to permit re-ranking: (/S(:,:)/)
        integer,intent(in),optional:: &
                 c(:) !c(n) component indices, or c(n+n) components ranges of each p
        real*8,intent(in),optional:: w(:)    !w(n+1) frequences for each p
        integer:: n,op(size(p)),cc(2*size(p)),ns,i
        real*8:: ww(size(p)+1); logical:: zovl; n=size(p)
        !process and check arguments
        call PropCheck('PropOverlap',S0,n,p, c=c,w=w, op=op,cc=cc,ww=ww,ne=ns,zovl=zovl)
        do i=1,ns; S(i) = 0d0*S0; enddo
        if (zovl) return !the result is zero
        call subr() !since bb(m) and oe(ne) and FF(ne)
    contains
        subroutine subr()
            integer:: os(ns); type(matrix):: SS(ns)
            call PropOrder(n,op,cc,0,0,ns, cc(1:0),os)
            !call core function ovl, then reorder back into S
            call ovl(n,(/(p(op(i)),i=1,n)/),(/(cc(op(i)),i=1,n),(cc(n+op(i)),i=1,n)/) &
                    ,(/(ww(op(i)),i=1,n),ww(n+1)/),S0,SS)
            do i=1,ns; S(i) = SS(os(i)); SS(os(i))=0; enddo
        end subroutine
    end subroutine


    subroutine PropReort(p,c,S0,D,DFD,b,e, w)
    !front for 'reort' below, checking dimensions and taking care of permutations (not yet, actually)
        character(*):: p(:) !p(n) perturbation lables
        integer:: n,c(:) !c(n*2) component range for each p
        integer:: m,b(:) !b(m) ranks of D and DFD
        type(matrix):: S0,D(product(b)) !S0 uneperturbed overlap and D(product(b)) perturbed density matrices
        type(matrix):: DFD(product(b)) !energy-weighted density matrices
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
        integer:: n,m,c(n+n),b(m); character(*):: p(n) !perturbation lables
        type(matrix):: D(product(b+1)) !density matrix expansion to contract over
        real*8,intent(out):: e(product(b)*product(c(n+1:n+n)-c(1:n)+1)) !result
        integer:: i,j,k,l,ii,jj,kk,ll,i1,j1,k1,l1,pb,pb1
        type(matrix):: A(6); real*8:: r(3*natom*3) !maximum scratch needed
        pb = product(b); pb1 = product(b+1)
        if (n==0) then
            if (m==0) stop 'PropEnergy: unperturbed energy requested'
            !The highest order Ds multiply the unperturbed F, which is unknown here,
            !so insist that no highest-order Ds are present
            if (.not.all((/(iszero(D(i)),i=pb1-pb+1,pb1)/))) &
                stop 'PropEnergy: unperturbed Fock matrix needed'
            if (m==2) then !stopping at m=2, for now
                call mat_init(A(1),nbas,nbas)
                if (do_dft()) call mat_init(A(2),nbas,nbas)
                do j=0,b(2)-1
                    call di_GET_GbDs(D(1+b(1)+1+j),A(1))
                    if (do_dft()) then
!radovan please leave:
!  call xc_contribution(2, D(1), (/D(1+b(1)+1+j)/), (/A(1)/))
                        call di_get_sigma_xc_cont(D(1+b(1)+1+j),D(1),D(1),A(2))
                        A(1) = A(1) + A(2)
                    endif
                    do i=1,b(1)
                        e(i+b(1)*j) = trab(A(1),D(1+i)) !factor 2 hidden in trab
                    enddo
                enddo; A(1)=0; if (do_dft()) A(2)=0
            elseif (m>2) then
                call quit('PropEnergy: m>2 not implemented at the moment')
            endif
        elseif (n==1.and.p(1)=='EL') then
#ifndef PRG_DIRAC
            if (m==0) call DIPNUC_ifc(natom,r(1:3),(/(0d0,i=1,9*natom)/))
            if (m==0) e(1:c(2)-c(1)+1) = e(1:c(2)-c(1)+1) - r(c(1):c(2)) !sign change
#else
            !fixme
#endif /* PRG_DIRAC */
            do j=0,c(2)-c(1); jj=c(1)+j
                A(1) = oneint(xyz(jj)//'DIPLEN ')
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*j) = e(i+pb*j) + A(1)**D(ii) !factor 2 hidden in **
                enddo; A(1)=0
            enddo
        elseif (n==1.and.p(1)=='MAGO') then
            do j=0,c(2)-c(1); jj=c(1)+j
                A(1) = oneint(xyz(jj)//'ANGMOM ')
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*j) = trab(A(1),D(ii))/2 !half, factor 2 hidden in trab
                enddo; A(1)=0
            enddo
        elseif (n==2.and.all(p==(/'MAGO','MAGO'/))) then
            do k=0,c(4)-c(2); kk=c(2)+k
                do j=0,c(3)-c(1); jj=c(1)+j
                    A(1) = oneint(xyz(min(jj,kk))//xyz(max(jj,kk))//'SUSCGO')
                    do i=1,pb; ii=pb1-pb+i
                        e(i+pb*j+pb*3*k) = trab(A(1),D(ii)) !factor 2 hidden in trab
                    enddo; A(1)=0
                enddo
            enddo
        elseif (n==1.and.p(1)=='MAG'.and.m/=0) then !no m==0 because because MAG imag
            !one-electron integral contribution
            do k=0,c(2)-c(1); kk=c(1)+k
                A(1) = oneint('dh/dB'//xyz(kk)//'  ')
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*k) = e(i+pb*k) + trab(A(1),D(ii))
                enddo; A(1)=0
            enddo
            !two electron integral contributions
            if (.not.all((/(iszero(D(i)),i=pb1-pb+1,pb1)/))) then
                call twofck('M',D(1:1),A(1:3))
                if (do_dft()) then
                    do i=1,3; call mat_init(A(3+i),nbas,nbas); enddo
                    call di_get_MagDeriv_FxD_DFT(A(4:6),D(1))
                    do i=1,3; A(i) = A(i) + A(3+i); A(3+i)=0; enddo
                endif
                do i=1,pb; ii=pb1-pb+i
                    do k=0,c(2)-c(1); kk=c(1)+k
                        e(i+pb*k) = e(i+pb*k) + trab(A(kk),D(ii))
                    enddo
                enddo; A(1)=0; A(2)=0; A(3)=0
            endif
            if (m==1) then
                !nothing more
            elseif (m==2) then
                do i=1,b(1); if (iszero(D(1+i))) cycle
                    call twofck('M',D(1+i:1+i),A(1:3))
                    if (do_dft()) then
                        do k=1,3; call mat_init(A(3+k),nbas,nbas); enddo
                        call di_get_MagDeriv_GxD_DFT(D(1),D(1+i),A(4:6))
                        do k=1,3; A(k) = A(k) + A(3+k); A(3+k)=0; enddo
                    endif
                    do j=0,b(2)-1; jj=1+b(1)+1+j
                        do k=0,c(2)-c(1); kk=c(1)+k
                            e(i+b(1)*j+pb*k) = e(i+b(1)*j+pb*k) + trab(A(kk),D(jj))
                    enddo; enddo; A(1)=0; A(2)=0; A(3)=0
                enddo
            else
                call quit('PropEnergy: MAG, m>2 not implemented yet!')
            endif
        elseif (n==1.and.p(1)=='ELGR') then
#ifndef PRG_DIRAC
            if (m==0) call QDRNUC_ifc(r(1:6))
            if (m==0) e(1:c(2)-c(1)+1) = e(1:c(2)-c(1)+1) - r(c(1):c(2)) !sign change
#else
            !fixme
#endif /* PRG_DIRAC */
            do k=0,c(2)-c(1); kk=c(1)+k
                k1=kk-1; if (kk>1) k1=k1+2; if (kk>3) k1=k1+1
                A(1) = oneint(xyz(1+mod(k1,3))//xyz(1+k1/3)//'THETA ')
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*k) = e(i+pb*k) + trab(A(1),D(ii))
                enddo; A(1)=0
            enddo
        elseif (n==1.and.p(1)=='GEO') then
            if (m==0) then !nuclear repulsion contribution
                call GRADNN_ifc(natom,r(1:3*natom))
                do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                    e(i) = e(i) + r(ii); enddo
            endif
            !one-electron integral contribution
            do k=0,c(2)-c(1); kk=c(1)+k
                A(1) = oneint('1DHAM'//c09(kk/100)//c09(mod(kk,100)/10)//c09(mod(kk,10)))
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*k) = e(i+pb*k) + trab(A(1),D(ii))
                enddo; A(1)=0
            enddo
            !two-electron integral contribution
            do i=1,pb; ii=pb1-pb+i
                call twoctr('G',D(1),D(ii),r(1:3*natom))
                r(1:3*natom) = r(1:3*natom)*8 !factor 8 on these integrals
                if (m==0) r(1:3*natom) = r(1:3*natom)/2 !another half for molgra
                if (do_dft()) then
                    if (m==0) call di_get_geomDeriv_molgrad_DFT(r(3*natom+1:3*natom*2),natom,D(1))
                    if (m/=0) call di_get_geomDeriv_FxD_DFT(r(3*natom+1:3*natom*2),natom,D(1),D(ii))
                    r(1:3*natom) = r(1:3*natom) + r(3*natom+1:3*natom*2)
                endif
                do k=0,c(2)-c(1); kk=c(1)+k
                    e(i+pb*k) = e(i+pb*k) + r(kk)
            enddo; enddo
            if (m==0.or.m==1) then !nothing more
            elseif (m==2) then !integrals over products of first order densities
                do j=0,b(2)-1; do i=1,b(1); jj=1+b(1)+1+j
                    call twoctr('G',D(1+i),D(jj),r(1:3*natom))
                    r(1:3*natom) = r(1:3*natom)*8 !factor 8 on these integrals
                    if (do_dft()) then
                        call di_get_geomDeriv_GxD_DFT(r(3*natom+1:3*natom*2),natom,D(1),D(1+i),D(jj))
                        r(1:3*natom) = r(1:3*natom) + r(3*natom+1:3*natom*2)
                    endif
                    do k=0,c(2)-c(1); kk=c(1)+k
                        e(i+b(1)*j+pb*k) = e(i+b(1)*j+pb*k) + r(kk) !factor 8 here
                enddo; enddo; enddo
            else
                call quit('PropEnergy: GEO, m>2 not implemented')
            endif
            if (do_dft()) print* !after all the "...integrated to nn electrons..." prints
        elseif (n==1.and.p(1)=='VEL') then
            do j=0,c(2)-c(1); jj=c(1)+j
                A(1) = oneint(xyz(jj)//'DIPVEL ')
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*j) = e(i+pb*j) + trab(A(1),D(ii)) !factor 2 hidden in **
                enddo; A(1)=0
            enddo
        elseif (n==2.and.all(p==(/'EL','MAG'/)).and.m/=0) then !no m==0 because MAG imag
            do l=0,c(4)-c(2); ll=c(2)+l
                do k=0,c(3)-c(1); kk=c(1)+k
                    A(1) = oneint(xyz(kk)//'-CM1 '//xyz(ll)//' ')
                    kk=k+l*(c(3)-c(1)+1)
                    do i=1,pb; ii=pb1-pb+i
                        e(i+pb*kk) = e(i+pb*kk) + trab(A(1),D(ii))
            enddo; A(1)=0; enddo; enddo
        elseif (n==2.and.all(p==(/'MAG','MAG'/))) then
            !one-electron integral contribution
            do l=0,c(4)-c(2); ll=c(2)+l
                do k=0,c(3)-c(1); kk=c(1)+k
                    A(1) = oneint(xyz(min(kk,ll))//xyz(max(kk,ll))//'dh/dB2')
                    kk=k+l*(c(3)-c(1)+1)
                    do i=1,pb; ii=pb1-pb+i
                        e(i+pb*kk) = e(i+pb*kk) + trab(A(1),D(ii))
                    enddo; A(1)=0
            enddo; enddo
            if (m==0.or..not.all((/(iszero(D(i)),i=pb1-pb+1,pb1)/))) then
                call twofck('MM',D(1:1),A(1:6))
                do i=1,pb; ii=pb1-pb+i
                    do l=0,c(4)-c(2); ll=c(2)+l
                        do k=0,c(3)-c(1); kk=c(1)+k
                            kk = min(kk,ll)+max(kk,ll)*(max(kk,ll)-1)/2
                            k1 = k+(c(3)-c(1)+1)*l
                            e(i+pb*k1) = e(i+pb*k1) + trab(A(kk),D(ii))/merge(2,1,m==0)
                enddo; enddo; enddo
                do i=1,6; A(i)=0; enddo
            endif
            if (m==0.or.m==1) then
                !nothing more
            elseif (m==2) then
                do i=1,b(1); if (iszero(D(1+i))) cycle
                    call twofck('MM',D(1+i:1+i),A(1:6))
                    do j=0,b(2)-1; jj=1+b(1)+1+j
                        do l=0,c(4)-c(2); ll=c(2)+l
                            do k=0,c(3)-c(1); kk=c(1)+k
                                kk = min(kk,ll)+max(kk,ll)*(max(kk,ll)-1)/2
                                k1 = k+(c(4)-c(2)+1)*l
                                e(i+b(1)*j+pb*k1) = e(i+b(1)*j+pb*k1) + trab(A(kk),D(jj))
                    enddo; enddo; enddo
                    do k=1,6; A(k)=0; enddo
                enddo
            else
                call quit('PropEnergy: MAG-MAG, m>2 not implemented yet!')
            endif
        elseif (n==2.and.all(p==(/'EL','GEO'/))) then
            if (m==0) then !nuclear repulsion contribution
                call DIPNUC_ifc(natom,(/0d0,0d0,0d0/),r(1:3*3*natom))
                do j=0,c(4)-c(2); jj=c(2)+j
                    do i=1,c(3)-c(1)+1; ii=c(1)+i-1
                        i1 = i+(c(3)-c(1)+1)*j
                        e(i1) = e(i1) - r(ii+3*(jj-1))
                enddo; enddo
            endif
            !one-electron integral contribution
            do k=0,c(4)-c(2); kk=c(2)+k
                do j=0,c(3)-c(1); jj=c(1)+j
                    A(1) = oneint(c09(kk/10)//c09(mod(kk,10))//' DPG '//xyz(jj))
!ajt                    A(1) = oneint(c09(kk/100)//c09(mod(kk,100)/10) &
!ajt                                //c09(mod(kk,10))//'DPG '//xyz(jj))
                    do i=1,pb; ii=pb1-pb+i
                        i1 = i+pb*(j+(c(3)-c(1)+1)*k)
                        e(i1) = e(i1) - trab(A(1),D(ii)) !factor -1
                    enddo; A(1)=0
            enddo; enddo
        elseif (n==2.and.all(p==(/'MAGO','GEO'/)).and.m/=0) then
            do k=0,c(4)-c(2); kk=c(2)+k
                do j=0,c(3)-c(1); jj=c(1)+j
                    A(1) = oneint(c09(kk/100)//c09(mod(kk,100)/10) &
                                  //c09(mod(kk,10))//'AMDR'//xyz(jj))
                    do i=1,pb; ii=pb1-pb+i
                        i1 = i+pb*(j+(c(3)-c(1)+1)*k)
                        e(i1) = e(i1) + trab(A(1),D(ii))/2 !factor 1/2
                    enddo; A(1)=0
            enddo; enddo
        elseif (n==2.and.all(p==(/'ELGR','GEO'/))) then
            if (m==0) then !nuclear repulsion
                call quit('PropEnergy: ELGR,GEO - nucrep conrib not implemented')
            endif
            !one-electron integral contribution
            do k=0,c(4)-c(2); kk=c(2)+k
                do j=0,c(3)-c(1); jj=c(1)+j
                    j1=jj-1; if (jj>1) j1=j1+2; if (jj>3) j1=j1+1
                    A(1) = oneint(c09(kk/10)//c09(mod(kk,10))//'QDG ' &
                                //xyz(1+mod(j1,3))//xyz(1+j1/3))
!ajt                    call oneint(c09(jj/100)//c09(mod(jj,100)/10)//c09(mod(jj,10))// &
!ajt                                'QDG'//xyz(1+(ii-1)/3)//xyz(1+mod(ii-1,3)),A(1))
                    do i=1,pb; ii=pb1-pb+i
                        i1 = i+pb*(j+(c(3)-c(1)+1)*k)
                        e(i1) = e(i1) - 3/2d0*trab(A(1),D(ii)) !factor -3/2
                    enddo; A(1)=0
                enddo
                if (c(1)==1.and.c(3)==6) then !remove trace
                    do i=1,pb
                        i1 = i+pb*(c(3)-c(1)+1)*k
                        r(1) = e(i1+pb*0)+e(i1+pb*2)+e(i1+pb*5) !trace
                        e(i1+pb*0) = e(i1+pb*0) - r(1)/3
                        e(i1+pb*2) = e(i1+pb*2) - r(1)/3
                        e(i1+pb*5) = e(i1+pb*5) - r(1)/3
                    enddo
                endif
            enddo
        elseif (n==3.and.all(p==(/'EL','MAG','MAG'/)).and.(m==0.or.m==1)) then
            do l=0,c(6)-c(3); ll=c(3)+l
                do k=0,c(5)-c(2); kk=c(2)+k
                    do j=0,c(4)-c(1); jj=c(1)+j
                        A(1) = oneint(xyz(jj)//'-CM2'//xyz(min(kk,ll))//xyz(max(kk,ll))//' ')
                        jj = j+(c(4)-c(1)+1)*(k+(c(5)-c(2)+1)*l)
                        do i=1,pb; ii=pb1-pb+i
                            e(i+pb*jj) = e(i+pb*jj) + trab(A(1),D(ii))
            enddo; A(1)=0; enddo; enddo; enddo
        else
            print *,'PropEnergy: no integrals for these perturbations: ', &
                    (p(i)//' ',i=1,n)
            call quit('PropEnergy: no such integrals')
        endif
    end subroutine


    subroutine fck(n,m,p,c,D,b,F)
        integer:: n,m,c(n+n),b(m); character(*):: p(n) !perturbation lables
        type(matrix):: D(product(b+1)) !density matrix expansion to contract over
        type(matrix),intent(inout):: F(product(b)*product(c(n+1:n+n)-c(1:n)+1)) !result
        integer:: i,j,k,ii,jj,kk,i1,pb,pb1; type(Matrix):: A(9)
        pb = product(b); pb1 = product(b+1)
        if (n==0) then
            if (m==0) stop 'PropFock: unperturbed Fock matrix requested'
            call mat_init(A(1),nbas,nbas)
            do i=1,pb; ii=pb1-pb+i
                if (iszero(D(ii))) cycle
                call di_GET_GbDs(D(ii),A(1))
                if (iszero(F(i))) F(i)=tiny(1d0)*D(1) !to get it allocated
                F(i) = F(i) + A(1)
                if (do_dft()) then !DFT contribution
!radovan please leave:
!                    call xc_contribution(2, D(1), (/D(ii)/), (/F(i)/))
                    call di_get_sigma_xc_cont(D(ii),D(1),D(1),A(1))
                    F(i) = F(i) + A(1)
                endif
            enddo
            if (do_dft().and.m==1) then
                !nothing more
            elseif (do_dft().and.m==2) then
                do j=0,b(2)-1; jj=1+b(1)+1+j
                    if (iszero(D(jj))) cycle
                    do i=1,b(1); ii=1+i
                        if (iszero(D(ii))) cycle
!radovan please leave:
!                        A(1) = tiny(1d0)*A(1)
!                        A(1)%pg_sym = 1 + ieor(D(ii)%pg_sym - 1, D(jj)%pg_sym - 1)
!                        call xc_contribution(3, D(1), (/D(ii), D(jj)/), A(1:1))
                        call di_get_T_xc_cont(D(1),D(ii),D(jj),A(1))
                        F(i+b(1)*j) = F(i+b(1)*j) + A(1)
                    enddo
                enddo
            elseif (do_dft()) then
                call quit('PropFock: DFT m==3 not omplemented')
            endif
            A(1) = 0
        elseif (n==1.and.p(1)=='EL') then !only m==0
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                A(1) = oneint(xyz(ii)//'DIPLEN ')
                F(i) = F(i) + A(1)
            enddo; A(1)=0
        elseif (n==1.and.p(1)=='MAGO'.and.m==0) then !only m==0
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                A(1) = oneint(xyz(ii)//'ANGMOM ')
                F(i) = F(i) + (1/2d0)*A(1) !factor here
            enddo; A(1)=0
        elseif (n==1.and.p(1)=='MAG') then
            if (m==0) then
                do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                    A(1) = oneint('dh/dB'//xyz(ii)//'  ')
                    F(i) = F(i) + A(1)
                enddo; A(1)=0
            endif
            do i=1,pb; ii=pb1-pb+i
                if (iszero(D(ii))) cycle
                call twofck('M',D(ii:ii),A(1:3))
                do k=0,c(2)-c(1); kk=c(1)+k
                    F(i+pb*k) = F(i+pb*k) + A(kk)
                enddo
                if (do_dft()) then
                    if (m==0) call di_get_MagDeriv_FxD_DFT(A(1:3),D(1))
                    if (m/=0) call di_get_MagDeriv_GxD_DFT(D(1),D(ii),A(1:3))
                    do k=0,c(2)-c(1); kk=c(1)+k
                        F(i+pb*k) = F(i+pb*k) + A(kk)
                    enddo
                endif; do k=1,3; A(k)=0; enddo
            enddo
        elseif (n==1.and.p(1)=='ELGR'.and.m==0) then !only m==0
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                i1=ii-1; if (ii>1) i1=i1+2; if (ii>3) i1=i1+1
                A(1) = oneint(xyz(1+mod(i1,3))//xyz(1+i1/3)//'THETA ')
                F(i) = F(i) + A(1)
            enddo; A(1)=0
        elseif (n==1.and.p(1)=='GEO') then
            if (m==0) then
                do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                    A(1) = oneint('1DHAM'//c09(ii/100)//c09(mod(ii,100)/10)//c09(mod(ii,10)))
                    F(i) = F(i) + A(1)
                enddo; A(1)=0
            endif
            do i=1,pb; ii=pb1-pb+i
                if (iszero(D(ii))) cycle
                do k=0,c(2)-c(1); kk=c(1)+k
                    if (k==0.or.mod(kk-1,3)==0) &
                        call twofck('G',D(ii:ii),A(1:3), a=1+(kk-1)/3)
                    F(i+pb*k) = F(i+pb*k) + A(1+mod(kk-1,3))
                enddo
                if (do_dft()) call quit('PropFock: GEO with DFT not implemented')
                do k=1,3; A(k)=0; enddo
            enddo
        elseif (n==2.and.all(p==(/'EL','MAG'/)).and.m==0) then
            do j=0,c(4)-c(2); jj=c(2)+j; 
                do i=1,c(3)-c(1)+1; ii=c(1)+i-1; kk=i+(c(3)-c(1)+1)*j
                    A(1) = oneint(xyz(ii)//'-CM1 '//xyz(jj)//' ')
                    F(kk) = F(kk) + A(1); A(1)=0
            enddo; enddo
        elseif (n==2.and.all(p==(/'MAGO','MAGO'/))) then !only m==0
            do j=0,c(4)-c(2); jj=c(2)+j
                do i=1,c(3)-c(1)+1; ii=c(1)+i-1; kk=i+(c(3)-c(1)+1)*j
                    A(1) = oneint(xyz(min(ii,jj))//xyz(max(ii,jj))//'SUSCGO')
                    F(kk) = F(kk) + A(1); A(1)=0
                enddo
            enddo
        elseif (n==2.and.all(p==(/'MAG','MAG'/))) then
!            call twofck('MM',D(1),(/(D(1),i=1,9)/))
        else
            print *,'PropFock: no integrals for these perturbations:', &
                    (' '//p(i),i=1,n)
            call quit('PropFock: no such integrals')
        endif
!        do i=1,9; if (isdef(A(i))) call mat_free(A(i)); enddo
    end subroutine


    subroutine ovl(n,p,c,w,S0,S)
    !frequency generalized perturbed overlap matrices
        integer:: n,c(n+n)
        character(*):: p(n) !perturbation lables
        real*8:: w(n+1) !frequencies
        type(matrix):: S0 !unperturbed overlap matrix, for context
        type(matrix),intent(inout):: S(product(c(n+1:n+n)-c(1:n)+1)) !result
        integer:: i,j,k,ii,jj,kk
        if (n==0) then
            stop 'PropOverlap: unperturbed overlap matrix requested'
        elseif (n==1.and.p(1)=='MAG') then
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                if (w(1)==0.and.w(2)==1) then
                    S(i) = oneint('dS/dB'//xyz(ii)//'  '); cycle; endif
                S(i) = (w(2)-w(1))*oneint('d<S|/dB'//xyz(ii)//'  ') &
                     + (w(2)+w(1))*oneint('d|S>/dB'//xyz(ii)//'  ')
            enddo
        elseif (n==1.and.p(1)=='GEO') then
            if (.not.all(w==(/0d0,1d0/))) &
                call quit('PropInt.PropOverlap: unexpected w')
            do i=0,c(2)-c(1); ii=c(1)+i
                S(1+i) = -oneint('1DOVL'//c09(ii/100)//c09(mod(ii,100)/10) &
                                //c09(mod(ii,10))) !notice minus sign
            enddo
        elseif (n==2.and.p(1)=='MAG'.and.p(2)=='MAG') then
            if (.not.all(w==(/0d0,0d0,1d0/))) &
                call quit('PropInt.PropOverlap: unexpected w')
            do j=0,c(4)-c(2); jj=c(2)+j
                do i=1,c(3)-c(1)+1; ii=c(1)+i-1
                    S(i+(c(3)-c(1)+1)*j) = oneint('dS/dB2'//xyz(min(ii,jj))//xyz(max(ii,jj)))
            enddo; enddo
        else
            print *,'PropOverlap: no integrals for these perturbations:', &
                    (' '//p(i),i=1,n); call quit('')
        endif
    end subroutine


    subroutine init(SD)
        type(matrix):: SD; integer:: NATOMS_ifc
        if (SD%nrow==nbas) return !stub
        nbas = SD%nrow
#ifndef PRG_DIRAC
        natom = NATOMS_ifc()
#else
        natom = 137 ! fixme
#endif /* PRG_DIRAC */
        tab(idx('GEO'))%rank = 3*natom
        !AllPrt(idx('NS'))%rank = depends on which isotopes
    end subroutine


    function oneint(l) result(A)
        character*8:: l; type(matrix):: A
        if (nbas==0) stop 'PropInt, oneint: nbas not initialized'
        if (undef(A)) call mat_init(A,nbas,nbas)
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
        integer:: lw,i,j,k; real*8:: Df(nbas**2),Ff(nbas**2,3*o),w(600000)
        if (nbas==0) stop 'twomatr: nbas not initialized'
#ifndef PRG_DIRAC
        call mat_to_full(D,1d0,Df); lw = 600000
        call GRCONT(w,lw,Ff,nbas**2*3*o,a/=0,a==0,o,a,.false.,.true.,Df,1)
#else
        call quit('error in twomatr')
#endif /* PRG_DIRAC */
        do i=1,3**o; if (undef(F(i))) call mat_init(F(i),nbas,nbas); enddo
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


    subroutine twofck(what,D,F, a,b)
    !what='G' : 3 geometry-derivative Fock matrices; atom a
    !what='M' : 3 magnetic-derivative Fock matrices
    !what='GG': 9 2nd geom.-derivative Fock matrices; atoms a,b
    !what='MM': 9 2nd magn.-derivative Fock matrices
        character(*):: what; type(matrix):: D(:),F(:)
        integer,optional:: a,b; integer:: i,j,k,nd,nf,n2,lw
        real*8,allocatable:: wrk(:)
        nd=size(D); nf=size(F); n2=nbas**2
        if     (what=='M' .and.nf==3*nd) then
            lw = 50*n2+10000*nbas
        elseif (what=='G' .and.nf==3*nd) then
            if (.not.present(a)) &
                call quit("PropInt.twofck: Geometry Fock requires atom 'a' specified")
            lw = 50*n2+10000*nbas
        elseif (what=='MM'.and.nf==6*nd) then
            lw = 50*n2+10000*nbas
        elseif (what=='GG'.and.nf==3*3*nd) then
            if (.not.present(a).or..not.present(b)) &
                call quit("PropInt.twock: 2nd ord. geom Fock reqires atoms 'a','b' specified")
            lw = (nd+2*nf+50)*n2+10000*nbas
        else
            call quit('PropInt.twofck: wrong size(F) for what='//what)
        endif
        allocate(wrk(lw))
        do i=1,nd; wrk(1+n2*(i-1):n2*i) = D(i)%elms; enddo
        !ajt Redundant: wrk(nbas**2+1:nbas**2*(nf+1)) = 0
        call GRCONT(wrk(1+n2*(nd+nf):lw),lw-n2*(nd+nf),wrk(1+n2*nd:n2*(nd+nf)) &
                   ,n2*nf,what(1:1)=='G',what(1:1)=='M',len(what),0 &
                   ,.false.,.true.,wrk(1:n2*nd),nd)
!        call GRCONT(wrk(nbas**2*(nf+1)+1:),size(wrk)-nbas**2*(nf+1) &
!                   ,wrk(nbas**2+1:nbas**2*(size(F)+1)),nbas**2*nf &
!                   ,what(1:1)=='G',what(1:1)=='M',len(what),0 &
!                   ,.false.,.true.,wrk(1:nbas**2),1)
        do i=1,nf
            if (undef(F(i))) call mat_init(F(i),nbas,nbas)
            if (iszero(F(i))) call mat_init(F(i),nbas,nbas)
            F(i)%elms = wrk(1+n2*(nd+i-1):n2*(nd+i))
        enddo
        deallocate(wrk)
    end subroutine


    subroutine twoctr(what,D,DD,e)
    !gm=T -> 3*NATOMS geometriske
    !gm=F -> tre magnetiske
        character(*):: what; type(matrix):: D,DD
        real*8,intent(inout):: e(:)
        real*8,allocatable:: wrk(:)
        if (nbas==0) stop 'PropInt.twoctr: nbas not initialized'
#ifndef PRG_DIRAC
        !ajt MagSus contraction doesn't work, so I disabled this. Use twofck instead
        if (what=='MM'.or.what=='MM') &
            call quit("PropInt.twoctr: 'M' or 'MM' not yet implemented")
        if     (what=='G' .and.size(e)==3*natom) then
            allocate(wrk(10000*nbas+50*nbas**2))
        elseif (what=='M' .and.size(e)==3) then
            allocate(wrk(10000*nbas+50*nbas**2))
        elseif (what=='GG'.and.size(e)==3*3*natom**2) then
            allocate(wrk(10000*nbas+50*nbas**2))
        elseif (what=='MM'.and.size(e)==3*3) then
            allocate(wrk(10000*nbas+50*nbas**2))
        else
            call quit('PropInt.twoctr: wrong size(e) for what='//what)
        endif
        wrk(        1:  nbas**2) = D%elms !call mat_to_full(D ,1d0,A(:,:,1))
        wrk(nbas**2+1:2*nbas**2) = DD%elms !call mat_to_full(DD,1d0,A(:,:,2))
        e=0; call GRCONT(wrk(2*nbas**2+1:),(size(wrk)-2*nbas**2),e,size(e) &
                        ,what(1:1)=='G',what(1:1)=='M',len(what),0 &
                        ,.true.,.false.,wrk(1:2*nbas**2),2)
        deallocate(wrk)
#else
        call quit('error in twoctr')
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


    subroutine PropCheck(who,D0,n,p, c,w,nd,b,o, op,cc,ww,m,ne,zerg,zfck,zovl,imag)
    !Process, p,c,w,nd,b finding the order p(op(k)) which sorts p by idx,
    !fill in from c the component ranges cc(1:n+n), and the frequencies ww(1:n+1).
    !Determine m from nd and b
        integer,intent(in):: n
        character(*),intent(in):: who,p(n)
        type(matrix),intent(in):: D0 !(or S0) for reference
        integer,intent(in),optional:: c(:),nd,b(:),o(:)
        real*8,intent(in),optional:: w(:)
        integer,intent(out),optional:: op(n),cc(n+n),m,ne
        real*8,intent(out),optional:: ww(n+1)
        logical,intent(out),optional:: zerg,zfck,zovl,imag
        integer:: ip(n),i,j,k; logical:: one,zer
        call init(D0)
        ip = (/(idx(p(i)),i=1,n)/)
        !find out if the perturbed integrals are zero,
        !and if it's a one-electron operator
        j = count((/(tab(ip(i))%lin ,i=1,n)/))
        k = count((/(tab(ip(i))%quad,i=1,n)/))
        one = (j>0.or.k>0)
        zer = (j>1.or.k>2.or.j==1.and.k>=1)
        if (j==0.and.k==2) then
            j=0; do i=1,n; if (tab(i)%quad) then
                if (j/=0) zer = (ip(i)/=ip(j))
            j=i; endif; enddo
        endif
        if (present(zerg)) zerg = zer
        if (present(zfck)) zfck = zer
        !check whether the perturbed overlap integrals are zero
        if (present(zovl)) &
            zovl = any((/(.not.tab(ip(i))%bas,i=1,n)/))
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
            !knowing m, Fock contributions may be zero
            if (one.and.present(zfck).and.m>0) zfck = .true.
        elseif (present(b).or.present(o)) then
            call QUIT(who//': if b or o present, so must nd and m be')
        endif
        !is the energy contribution
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
!radovan: check whether this is not minloc (was in my version)
                k(1:1) = maxloc(kk(1:i,3))
                j = k(1); k = kk(j,:)
                kk(j:k(3)-1,:) = kk(j+1:k(3),:)
                kk(k(3),:) = k
            enddo
        endif
        !construct oe from dimensions kk(:,1) and pitches kk(:,2)
        oe(1) = 1; j = 1
        do i=1,m+n
            oe(j+1:j*kk(i,1)) = (/(oe(1:j)+l*kk(i,2),l=1,kk(i,1)-1)/)
            j = j*kk(i,1)
        enddo
    end subroutine


end module PropInt
