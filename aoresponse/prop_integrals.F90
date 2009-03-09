!Copyright 2009 Andreas J. Thorvaldsen
!This file is made available under the terms of the
!GNU Lesser General Public License.

module prop_integ
!routines for property contractions and matrices

!radovan: this file should not (but presently does) contain
!         any explicit (de)allocations (like in twoint and twoctr)
!         these are specific to programs and matrix/algebra dimensions
!         they should be moved to program specific interfaces
!ajt: mat_init(A,B,tb) is now supplied by matrix_defop and has
!     replaced all mat_init(A,nrow,ncol)

    use matrix_defop
#ifndef PRG_DIRAC
    use dal_interface,only: di_GET_GbDs,di_read_operator_int,di_get_sigma_xc_cont,di_get_T_xc_cont &
                     ,di_get_MagDeriv_FxD_DFT,di_get_MagDeriv_GxD_DFT &
                     ,di_get_geomDeriv_molgrad_DFT,di_get_geomDeriv_FxD_DFT,di_get_geomDeriv_GxD_DFT
    use matrix_operations,only: mat_to_full,mat_set_from_full
#else
    use dirac_interface_m
#endif /* PRG_DIRAC */

    implicit none
    public prop_energy,prop_reort,prop_fock,prop_overlap
    public prop_imag,prop_rank,prop_check
    private

    type prt !private: 'perturbation'
        character*4  :: code !four-letter abbreviation
        character*64 :: name !long name
        integer      :: rank !number of components (when known, 0 otw.)
        logical      :: imag !imaginary (1st ord.) perturbed integrals
        logical      :: bas  !basis dependent (sa. GEO and MF)
        logical      :: lin  !one-electron operator linear in field strength (EL)
        logical      :: quad !one-electron operator quadratic in field strength (MAGO)
    end type

    !starting with two fields
    logical,parameter:: T=.true.,F=.false.
    type(prt) :: tab(6) = &                                 !rk im ba ln qu
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


    subroutine prop_energy(p,D, e, c,b,o)
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
        call prop_check('prop_energy',D(1),n,p, c=c,nd=size(D),b=b,o=o &
                      ,op=op,cc=cc,m=m,ne=ne,zerg=zerg,imag=imag)
        if (zerg.or.m==0.and.imag) return !the result is zero
        call subr() !since bb(m) and oe(ne) and ee(ne)
    contains
        subroutine subr()
            integer:: bb(m),oe(ne); real*8:: ee(ne); ee=0
            call prop_order(n,op,cc,m,size(D),ne, bb,oe, b=b,o=o)
            !call core function erg
            call energy_core(n,m,(/(p(op(i)),i=1,n)/) &
                    ,(/(cc(op(i)),i=1,n),(cc(n+op(i)),i=1,n)/),D,bb,ee)
            !add the reordered output to e
            e(1:ne) = e(1:ne) + (/(ee(oe(i)),i=1,ne)/)
        end subroutine
    end subroutine


    subroutine prop_fock(p,D, F, c,b,o)
    !front for 'fock_core' below, checking dimensions and taking care permutations (not yet, actually)
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
        call prop_check('prop_fock',D(1),n,p, c=c,nd=size(D),b=b,o=o &
                      ,op=op,cc=cc,m=m,ne=nf,zfck=zfck)
        if (zfck) return !the result is zero
        call subr() !since bb(m) and oe(ne) and FF(ne)
    contains
        subroutine subr()
            integer:: bb(m),of(nf); type(matrix):: FF(nf)
            call prop_order(n,op,cc,m,size(D),nf, bb,of, b=b,o=o)
            !order F into FF, call core function fock_core, then reorder back into F
            do i=1,nf; FF(of(i)) = F(i); F(i)=0; end do
            call fock_core(n,m,(/(p(op(i)),i=1,n)/) &
                    ,(/(cc(op(i)),i=1,n),(cc(n+op(i)),i=1,n)/),D,bb,FF)
            do i=1,nf; F(i) = FF(of(i)); FF(of(i))=0; end do
        end subroutine
    end subroutine


    subroutine prop_overlap(p,S0, S, c,w)
    !front for 'overlap_core' below, checking dimensions. Argument S0 is for context only
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
        call prop_check('prop_overlap',S0,n,p, c=c,w=w, op=op,cc=cc,ww=ww,ne=ns,zovl=zovl)
        do i=1,ns; S(i) = 0d0*S0; end do !delete any previous S(:)
        if (zovl) return !the result is zero
        call subr() !since bb(m) and oe(ne) and FF(ne)
    contains
        subroutine subr()
            integer:: os(ns); type(matrix):: SS(ns)
            call prop_order(n,op,cc,0,0,ns, cc(1:0),os)
            !call core function ovl, then reorder back into S
            do i=1,ns; SS(i) = 0d0*S0; end do
            call overlap_core(n,(/(p(op(i)),i=1,n)/),(/(cc(op(i)),i=1,n),(cc(n+op(i)),i=1,n)/) &
                    ,(/(ww(op(i)),i=1,n),ww(n+1)/),S0,SS)
            S(1:ns) = (/(SS(os(i)),i=1,ns)/) !reorder (aliasing)
        end subroutine
    end subroutine


    subroutine prop_reort(p,c,S0,D,DFD,b,e, w)
    !front for 'reort' below, checking dimensions and taking care of permutations (not yet, actually)
        character(*):: p(:) !p(n) perturbation lables
        integer:: n,c(:) !c(n*2) component range for each p
        integer:: m,b(:) !b(m) ranks of D and DFD
        type(matrix):: S0,D(product(b)) !S0 uneperturbed overlap and
                                        !D(:) perturbed density matrices
        type(matrix):: DFD(product(b)) !energy-weighted density matrices
        real*8,intent(inout):: e(product(b)*product(c(size(p)+1:)-c(:size(p)))) !result is added to e
        real*8,optional:: w(:) !frequencies w(n+1)
        integer:: i,j,k,ip(size(p))
        real*8:: ww(size(p)+1)
        n = size(p); m = size(b)
        if (size(c)/=n*2) &
            stop 'prop_reort: size(c)/=2*size(p)'
        if (present(w)) then
            if (size(w)/=n+1) stop 'prop_reort: wrong number of frequencies, size(w)/=size(p)+1'
            ww = w
        end if
        call init(S0)
    end subroutine


    subroutine energy_core(n,m,p,c,D,b,e)
        integer:: n,m,c(n+n),b(m); character(*):: p(n) !perturbation lables
        type(matrix):: D(product(b+1)) !density matrix expansion to contract over
        real*8,intent(out):: e(product(b)*product(c(n+1:n+n)-c(1:n)+1)) !result
        integer:: i,j,k,l,ii,jj,kk,ll,i1,j1,k1,l1,pb,pb1
        type(matrix):: A(6); real*8:: r(3*natom*3) !maximum scratch needed
        do i=1,size(A); A(i)=0d0*D(1); end do
        pb = product(b); pb1 = product(b+1)
        if (n==0) then
            if (m==0) stop 'prop_energy: unperturbed energy requested'
            !The highest order Ds multiply the unperturbed F, which is unknown here,
            !so insist that no highest-order Ds are present
            if (.not.all((/(iszero(D(i)),i=pb1-pb+1,pb1)/))) &
                stop 'prop_energy: unperturbed Fock matrix needed'
            if (m==2) then !stopping at m=2, for now
                call mat_init(A(1),A(1)) !de_GET_GbDs expects initialized
                if (do_dft()) call mat_init(A(2),A(2)) !so does di_get_sigma..
                do j=0,b(2)-1
                    call di_GET_GbDs(D(1+b(1)+1+j),A(1))
                    if (do_dft()) then
!radovan please leave:
!  call xc_contribution(2, D(1), (/D(1+b(1)+1+j)/), (/A(1)/))
                        call di_get_sigma_xc_cont(D(1+b(1)+1+j),D(1),D(1),A(2))
                        A(1) = A(1) + A(2)
                    end if
                    do i=1,b(1)
                        e(i+b(1)*j) = tr(A(1),D(1+i)) !factor 2 hidden in tr
                    end do
                end do
            else if (m>2) then
                call quit('prop_energy: m>2 not implemented at the moment')
            end if
        else if (n==1.and.p(1)=='EL') then
            if (m==0) call DIPNUC_ifc(natom,r(1:3),(/(0d0,i=1,9*natom)/))
            if (m==0) e(1:c(2)-c(1)+1) = e(1:c(2)-c(1)+1) - r(c(1):c(2)) !sign change
            do j=0,c(2)-c(1); jj=c(1)+j
                call oneint(xyz(jj)//'DIPLEN ',A(1))
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*j) = e(i+pb*j) + dot(A(1),D(ii)) !factor 2 hidden in **
                end do
            end do
        else if (n==1.and.p(1)=='MAGO') then
            do j=0,c(2)-c(1); jj=c(1)+j
                call oneint(xyz(jj)//'ANGMOM ',A(1))
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*j) = tr(A(1),D(ii))/2 !half, factor 2 hidden in tr
                end do
            end do
        else if (n==2.and.all(p==(/'MAGO','MAGO'/))) then
            do k=0,c(4)-c(2); kk=c(2)+k
                do j=0,c(3)-c(1); jj=c(1)+j
                    call oneint(xyz(min(jj,kk))//xyz(max(jj,kk))//'SUSCGO',A(1))
                    do i=1,pb; ii=pb1-pb+i
                        e(i+pb*j+pb*3*k) = tr(A(1),D(ii)) !factor 2 hidden in tr
                    end do
                end do
            end do
        else if (n==1.and.p(1)=='MAG'.and.m/=0) then !no m==0 because because MAG imag
            !one-electron integral contribution
            do k=0,c(2)-c(1); kk=c(1)+k
                call oneint('dh/dB'//xyz(kk)//'  ',A(1))
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*k) = e(i+pb*k) + tr(A(1),D(ii))
                end do
            end do
            !two electron integral contributions
            if (.not.all((/(iszero(D(i)),i=pb1-pb+1,pb1)/))) then
                call twofck('M',D(1:1),A(1:3))
                if (do_dft()) then
                    do i=1,3; call mat_init(A(3+i),A(3+i)); end do
                    call di_get_MagDeriv_FxD_DFT(A(4:6),D(1))
                    do i=1,3; A(i) = A(i) + A(3+i); A(3+i)=0; end do
                end if
                do i=1,pb; ii=pb1-pb+i
                    do k=0,c(2)-c(1); kk=c(1)+k
                        e(i+pb*k) = e(i+pb*k) + tr(A(kk),D(ii))
                    end do
                end do
            end if
            if (m==1) then
                !nothing more
            else if (m==2) then
                do i=1,b(1); if (iszero(D(1+i))) cycle
                    call twofck('M',D(1+i:1+i),A(1:3))
                    if (do_dft()) then
                        do k=4,6 !since di_get_MagDeriv expects initialized
                            if (iszero(A(k))) call mat_init(A(k),A(k))
                        end do
                        call di_get_MagDeriv_GxD_DFT(D(1),D(1+i),A(4:6))
                        do k=1,3; A(k) = A(k) - A(3+k); end do !sign on KS contrib
                    end if
                    do j=0,b(2)-1; jj=1+b(1)+1+j
                        do k=0,c(2)-c(1); kk=c(1)+k
                            e(i+b(1)*j+pb*k) = e(i+b(1)*j+pb*k) + tr(A(kk),D(jj))
                        end do
                    end do
                end do
            else
                call quit('prop_energy: MAG, m>2 not implemented yet!')
            end if
        else if (n==1.and.p(1)=='ELGR') then
            if (m==0) call QDRNUC_ifc(r(1:6))
            if (m==0) e(1:c(2)-c(1)+1) = e(1:c(2)-c(1)+1) - r(c(1):c(2)) !sign change
            do k=0,c(2)-c(1); kk=c(1)+k
                k1=kk-1; if (kk>1) k1=k1+2; if (kk>3) k1=k1+1
                call oneint(xyz(1+mod(k1,3))//xyz(1+k1/3)//'THETA ',A(1))
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*k) = e(i+pb*k) + tr(A(1),D(ii))
                end do
            end do
        else if (n==1.and.p(1)=='GEO') then
            if (m==0) then !nuclear repulsion contribution
                call GRADNN_ifc(natom,r(1:3*natom))
                do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                    e(i) = e(i) + r(ii); end do
            end if
            !one-electron integral contribution
            do k=0,c(2)-c(1); kk=c(1)+k
#ifndef PRG_DIRAC
                call oneint('1DHAM'//c09(kk/100)//c09(mod(kk,100)/10)//c09(mod(kk,10)),A(1))
#else
                call oneint('GEO1KIN'              &
                            //c09(kk/100)          &
                            //c09(mod(kk, 100)/10) &
                            //c09(mod(kk,10)),     &
                            A(1))
!               add d1nuc
#endif /* PRG_DIRAC */
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*k) = e(i+pb*k) + tr(A(1),D(ii))
                end do
            end do
            !two-electron integral contribution
            do i=1,pb; ii=pb1-pb+i
                call twoctr('G',D(1),D(ii),r(1:3*natom))
                r(1:3*natom) = r(1:3*natom)*8 !factor 8 on these integrals
                if (m==0) r(1:3*natom) = r(1:3*natom)/2 !another half for molgra
                if (do_dft()) then
                    if (m==0) call di_get_geomDeriv_molgrad_DFT(r(3*natom+1:3*natom*2),natom,D(1))
                    if (m/=0) call di_get_geomDeriv_FxD_DFT(r(3*natom+1:3*natom*2),natom,D(1),D(ii))
                    r(1:3*natom) = r(1:3*natom) + r(3*natom+1:3*natom*2)
                end if
                do k=0,c(2)-c(1); kk=c(1)+k
                    e(i+pb*k) = e(i+pb*k) + r(kk)
                end do
            end do
            if (m==0.or.m==1) then !nothing more
            else if (m==2) then !integrals over products of first order densities
                do j=0,b(2)-1; jj=1+b(1)+1+j
                    do i=1,b(1)
                        call twoctr('G',D(1+i),D(jj),r(1:3*natom))
                        r(1:3*natom) = r(1:3*natom)*8 !factor 8 on these integrals
                        if (do_dft()) then
                            call di_get_geomDeriv_GxD_DFT(r(3*natom+1:3*natom*2) &
                                                        ,natom,D(1),D(1+i),D(jj))
                            r(1:3*natom) = r(1:3*natom) + r(3*natom+1:3*natom*2)
                        end if
                        do k=0,c(2)-c(1); kk=c(1)+k
                            e(i+b(1)*j+pb*k) = e(i+b(1)*j+pb*k) + r(kk) !factor 8 here
                        end do
                    end do
                end do
            else
                call quit('prop_energy: GEO, m>2 not implemented')
            end if
            if (do_dft()) print* !after all the "...integrated to nn electrons..." prints
        else if (n==1.and.p(1)=='VEL') then
            do j=0,c(2)-c(1); jj=c(1)+j
                call oneint(xyz(jj)//'DIPVEL ',A(1))
                do i=1,pb; ii=pb1-pb+i
                    e(i+pb*j) = e(i+pb*j) + tr(A(1),D(ii)) !factor 2 hidden in **
                end do
            end do
#ifdef FIXME_G95
        else if (n==2.and.all(p==(/'EL','MAG'/)).and.m/=0) then !no m==0 because MAG imag
            do l=0,c(4)-c(2); ll=c(2)+l
                do k=0,c(3)-c(1); kk=c(1)+k
                    call oneint(xyz(kk)//'-CM1 '//xyz(ll)//' ',A(1))
                    kk=k+l*(c(3)-c(1)+1)
                    do i=1,pb; ii=pb1-pb+i
                        e(i+pb*kk) = e(i+pb*kk) + tr(A(1),D(ii))
                    end do
                end do
            end do
#endif
        else if (n==2.and.all(p==(/'MAG','MAG'/))) then
            !one-electron integral contribution
            do l=0,c(4)-c(2); ll=c(2)+l
                do k=0,c(3)-c(1); kk=c(1)+k
                    call oneint(xyz(min(kk,ll))//xyz(max(kk,ll))//'dh/dB2',A(1))
                    kk=k+l*(c(3)-c(1)+1)
                    do i=1,pb; ii=pb1-pb+i
                        e(i+pb*kk) = e(i+pb*kk) + tr(A(1),D(ii))
                    end do
                end do
            end do
            if (m==0.or..not.all((/(iszero(D(i)),i=pb1-pb+1,pb1)/))) then
                call twofck('MM',D(1:1),A(1:6))
                do i=1,pb; ii=pb1-pb+i
                    do l=0,c(4)-c(2); ll=c(2)+l
                        do k=0,c(3)-c(1); kk=c(1)+k
                            kk = min(kk,ll)+max(kk,ll)*(max(kk,ll)-1)/2
                            k1 = k+(c(3)-c(1)+1)*l
                            e(i+pb*k1) = e(i+pb*k1) + tr(A(kk),D(ii))/merge(2,1,m==0)
                        end do
                    end do
                end do
            end if
            if (m==0.or.m==1) then
                !nothing more
            else if (m==2) then
                do i=1,b(1); if (iszero(D(1+i))) cycle
                    call twofck('MM',D(1+i:1+i),A(1:6))
                    do j=0,b(2)-1; jj=1+b(1)+1+j
                        do l=0,c(4)-c(2); ll=c(2)+l
                            do k=0,c(3)-c(1); kk=c(1)+k
                                kk = min(kk,ll)+max(kk,ll)*(max(kk,ll)-1)/2
                                k1 = k+(c(4)-c(2)+1)*l
                                e(i+b(1)*j+pb*k1) = e(i+b(1)*j+pb*k1) + tr(A(kk),D(jj))
                            end do
                        end do
                    end do
                end do
            else
                call quit('prop_energy: MAG-MAG, m>2 not implemented yet!')
            end if
#ifdef FIXME_G95
        else if (n==2.and.all(p==(/'EL','GEO'/))) then
            if (m==0) then !nuclear repulsion contribution
                call DIPNUC_ifc(natom,(/0d0,0d0,0d0/),r(1:3*3*natom))
                do j=0,c(4)-c(2); jj=c(2)+j
                    do i=1,c(3)-c(1)+1; ii=c(1)+i-1
                        i1 = i+(c(3)-c(1)+1)*j
                        e(i1) = e(i1) - r(ii+3*(jj-1))
                    end do
                end do
            end if
            !one-electron integral contribution
            do k=0,c(4)-c(2); kk=c(2)+k
                do j=0,c(3)-c(1); jj=c(1)+j
                    call oneint(c09(kk/10)//c09(mod(kk,10))//' DPG '//xyz(jj),A(1))
!ajt                    A(1) = oneint(c09(kk/100)//c09(mod(kk,100)/10) &
!ajt                                //c09(mod(kk,10))//'DPG '//xyz(jj))
                    do i=1,pb; ii=pb1-pb+i
                        i1 = i+pb*(j+(c(3)-c(1)+1)*k)
                        e(i1) = e(i1) - tr(A(1),D(ii)) !factor -1
                    end do
                end do
            end do
        else if (n==2.and.all(p==(/'MAGO','GEO'/)).and.m/=0) then
            do k=0,c(4)-c(2); kk=c(2)+k
                do j=0,c(3)-c(1); jj=c(1)+j
                    call oneint(c09(kk/100)//c09(mod(kk,100)/10) &
                              //c09(mod(kk,10))//'AMDR'//xyz(jj),A(1))
                    do i=1,pb; ii=pb1-pb+i
                        i1 = i+pb*(j+(c(3)-c(1)+1)*k)
                        e(i1) = e(i1) + tr(A(1),D(ii))/2 !factor 1/2
                    end do
                end do
            end do
        else if (n==2.and.all(p==(/'ELGR','GEO'/))) then
            if (m==0) then !nuclear repulsion
                call quit('prop_energy: ELGR,GEO - nucrep conrib not implemented')
            end if
            !one-electron integral contribution
            do k=0,c(4)-c(2); kk=c(2)+k
                do j=0,c(3)-c(1); jj=c(1)+j
                    j1=jj-1; if (jj>1) j1=j1+2; if (jj>3) j1=j1+1
                    call oneint(c09(kk/10)//c09(mod(kk,10))//'QDG ' &
                              //xyz(1+mod(j1,3))//xyz(1+j1/3),A(1))
!ajt                    call oneint(c09(jj/100)//c09(mod(jj,100)/10)//c09(mod(jj,10))// &
!ajt                                'QDG'//xyz(1+(ii-1)/3)//xyz(1+mod(ii-1,3)),A(1))
                    do i=1,pb; ii=pb1-pb+i
                        i1 = i+pb*(j+(c(3)-c(1)+1)*k)
                        e(i1) = e(i1) - 3/2d0*tr(A(1),D(ii)) !factor -3/2
                    end do
                end do
                if (c(1)==1.and.c(3)==6) then !remove trace
                    do i=1,pb
                        i1 = i+pb*(c(3)-c(1)+1)*k
                        r(1) = e(i1+pb*0)+e(i1+pb*2)+e(i1+pb*5) !trace
                        e(i1+pb*0) = e(i1+pb*0) - r(1)/3
                        e(i1+pb*2) = e(i1+pb*2) - r(1)/3
                        e(i1+pb*5) = e(i1+pb*5) - r(1)/3
                    end do
                end if
            end do
        else if (n==3.and.all(p==(/'EL','MAG','MAG'/)).and.(m==0.or.m==1)) then
            do l=0,c(6)-c(3); ll=c(3)+l
                do k=0,c(5)-c(2); kk=c(2)+k
                    do j=0,c(4)-c(1); jj=c(1)+j
                        call oneint(xyz(jj)//'-CM2'//xyz(min(kk,ll))//xyz(max(kk,ll))//' ',A(1))
                        jj = j+(c(4)-c(1)+1)*(k+(c(5)-c(2)+1)*l)
                        do i=1,pb; ii=pb1-pb+i
                            e(i+pb*jj) = e(i+pb*jj) + tr(A(1),D(ii))
                        end do
                    end do
                end do
            end do
#endif
        else
            print *,'prop_energy: no integrals for these perturbations: ', &
                    (p(i)//' ',i=1,n)
            call quit('prop_energy: no such integrals')
        end if
        A(:) = 0 !delete scratch matrices
    end subroutine


    subroutine fock_core(n,m,p,c,D,b,F)
        integer:: n,m,c(n+n),b(m); character(*):: p(n) !perturbation lables
        type(matrix):: D(product(b+1)) !density matrix expansion to contract over
        type(matrix),intent(inout):: F(product(b)*product(c(n+1:n+n)-c(1:n)+1)) !result
        integer:: i,j,k,ii,jj,kk,i1,pb,pb1; type(Matrix):: A(9)
        do i=1,size(A); A(i)=0d0*D(1); end do
        pb = product(b); pb1 = product(b+1)
        if (n==0) then
            if (m==0) stop 'prop_fock: unperturbed Fock matrix requested'
            call mat_init(A(1),A(1)) !di_GET_GbDs expects initalized
            do i=1,pb; ii=pb1-pb+i
                if (iszero(D(ii))) cycle
                call di_GET_GbDs(D(ii),A(1))
                F(i) = F(i) + A(1)
                if (do_dft()) then !DFT contribution
!radovan please leave:
!                    call xc_contribution(2, D(1), (/D(ii)/), (/F(i)/))
                    call di_get_sigma_xc_cont(D(ii),D(1),D(1),A(1))
                    F(i) = F(i) + A(1)
                end if
            end do
            if (do_dft().and.m==1) then
                !nothing more
            else if (do_dft().and.m==2) then
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
                    end do
                end do
            else if (do_dft()) then
                call quit('prop_fock: DFT m==3 not omplemented')
            end if
        else if (n==1.and.p(1)=='EL') then !only m==0
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                call oneint(xyz(ii)//'DIPLEN ',A(1))
                F(i) = F(i) + A(1)
            end do
        else if (n==1.and.p(1)=='MAGO'.and.m==0) then !only m==0
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                call oneint(xyz(ii)//'ANGMOM ',A(1))
                F(i) = F(i) + (1/2d0)*A(1) !factor 1/2 here
            end do
        else if (n==1.and.p(1)=='MAG') then
            if (m==0) then
                do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                    call oneint('dh/dB'//xyz(ii)//'  ',A(1))
                    F(i) = F(i) + A(1)
                end do
            end if
            do i=1,pb; ii=pb1-pb+i
                if (iszero(D(ii))) cycle
                call twofck('M',D(ii:ii),A(1:3))
                do k=0,c(2)-c(1); kk=c(1)+k
                    F(i+pb*k) = F(i+pb*k) + A(kk)
                end do
                if (do_dft()) then
                    if (m==0) call di_get_MagDeriv_FxD_DFT(A(1:3),D(1))
                    if (m/=0) call di_get_MagDeriv_GxD_DFT(D(1),D(ii),A(1:3))
                    do k=0,c(2)-c(1); kk=c(1)+k
                        F(i+pb*k) = F(i+pb*k) + A(kk)
                    end do
                end if
            end do
        else if (n==1.and.p(1)=='ELGR'.and.m==0) then !only m==0
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                i1=ii-1; if (ii>1) i1=i1+2; if (ii>3) i1=i1+1
                call oneint(xyz(1+mod(i1,3))//xyz(1+i1/3)//'THETA ',A(1))
                F(i) = F(i) + A(1)
            end do
        else if (n==1.and.p(1)=='GEO') then
            if (m==0) then
                do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                    call oneint('1DHAM'//c09(ii/100)//c09(mod(ii,100)/10)//c09(mod(ii,10)),A(1))
                    F(i) = F(i) + A(1)
                end do
            end if
            do i=1,3; call mat_init(A(i),D(1)); end do
            do i=1,pb; ii=pb1-pb+i
                if (iszero(D(ii))) cycle
                do k=0,c(2)-c(1); kk=c(1)+k
                    if (k==0.or.mod(kk-1,3)==0) &
                        call twofck('G',D(ii:ii),A(1:3), a=1+(kk-1)/3)
                    F(i+pb*k) = F(i+pb*k) + A(1+mod(kk-1,3))
                end do
                if (do_dft()) call quit('prop_fock: GEO with DFT not implemented')
            end do
#ifdef FIXME_G95
        else if (n==2.and.all(p==(/'EL','MAG'/)).and.m==0) then
            do j=0,c(4)-c(2); jj=c(2)+j;
                do i=1,c(3)-c(1)+1; ii=c(1)+i-1; kk=i+(c(3)-c(1)+1)*j
                    call oneint(xyz(ii)//'-CM1 '//xyz(jj)//' ',A(1))
                    F(kk) = F(kk) + A(1)
                end do
            end do
        else if (n==2.and.all(p==(/'MAGO','MAGO'/))) then !only m==0
            do j=0,c(4)-c(2); jj=c(2)+j
                do i=1,c(3)-c(1)+1; ii=c(1)+i-1; kk=i+(c(3)-c(1)+1)*j
                    call oneint(xyz(min(ii,jj))//xyz(max(ii,jj))//'SUSCGO',A(1))
                    F(kk) = F(kk) + A(1)
                end do
            end do
        else if (n==2.and.all(p==(/'MAG','MAG'/))) then
!            call twofck('MM',D(1),(/(D(1),i=1,9)/))
#endif
        else
            print *,'prop_fock: no integrals for these perturbations:', &
                    (' '//p(i),i=1,n)
            call quit('prop_fock: no such integrals')
        end if
        A(:) = 0
    end subroutine


    subroutine overlap_core(n,p,c,w,S0,S)
    !frequency generalized perturbed overlap matrices
        integer     ,intent(in)    :: n,c(n+n)
        character(*),intent(in)    :: p(n) !perturbation lables
        real*8      ,intent(in)    :: w(n+1) !frequencies
        type(matrix),intent(in)    :: S0 !unperturbed overlap matrix, for context
        type(matrix),intent(inout) :: S(product(c(n+1:n+n)-c(1:n)+1)) !result
        integer                    :: i,j,k,ii,jj,kk
        type(matrix)               :: A
        A = 0d0*S0 !scratch matrix
        if (n==0) then
            stop 'prop_overlap: unperturbed overlap matrix requested'
        else if (n==1.and.p(1)=='MAG') then
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
                if (w(1)==0.and.w(2)==1) then
                    call oneint('dS/dB'//xyz(ii)//'  ',S(i))
                else
                    call oneint('d<S|/dB'//xyz(ii),S(i))
                    call oneint('d|S>/dB'//xyz(ii),A)
                    S(i) = (w(2)-w(1))*S(i) + (w(2)+w(1))*A
                end if
            end do
        else if (n==1.and.p(1)=='GEO') then
            if (.not.all(w==(/0d0,1d0/))) &
                call quit('PropInt.prop_overlap: unexpected w')
            do i=1,c(2)-c(1)+1; ii=c(1)+i-1
#ifndef PRG_DIRAC
!fixme
                call oneint('1DOVL'//c09(ii/100)//c09(mod(ii,100)/10) &
                                   //c09(mod(ii,10)),S(i))
#endif /* PRG_DIRAC */
                S(i) = -S(i) !notice minus sign
            end do
        else if (n==2.and.p(1)=='MAG'.and.p(2)=='MAG') then
            do j=0,c(4)-c(2); jj=c(2)+j !loop over requested second components
                do i=1,c(3)-c(1)+1; ii=c(1)+i-1 !loop over first components
                    if (w(1)==0.and.w(2)==0.and.w(3)==1) then
                        call oneint('dS/dB2'//xyz(min(ii,jj))//xyz(max(ii,jj)) &
                                   ,S(i+(c(3)-c(1)+1)*j))
                    else
!ajt to dj: Does this look alright? Can you modify it so that the lables are correct?
!           I have assumed that '<<S/B2xy' is a label, and not '<<S/B2yx' (hence min before max)
!                     S(i+(c(3)-c(2))*j) =                                                    &
!                         (w(3)-w(2)-w(1))*oneint('<<S/B2'//xyz(min(ii,jj))//xyz(max(ii,jj))) &
!                       + (w(3)-w(2)+w(1))*oneint('<>S/B2'//xyz(ii)//xyz(jj))                 &
!                       + (w(3)+w(2)-w(1))*oneint('<>S/B2'//xyz(jj)//xyz(ii))                 &
!                       + (w(3)+w(2)+w(1))*oneint('>>S/B2'//xyz(min(ii,jj))//xyz(max(ii,jj)))
!ajt zero for now
                        S(i+(c(3)-c(2))*j) = 0d0*S0
                    end if
                end do
            end do
        else
            print *,'prop_overlap: no integrals for these perturbations:', &
                    (' '//p(i),i=1,n); call quit('')
        end if
        A = 0 !delete scratch matrix
    end subroutine


    subroutine init(SD)
    !ajt: this should be removed, somehow
        type(matrix):: SD
        if (SD%nrow==nbas) return !stub
        nbas = SD%nrow
        call get_natoms(natom)
        tab(idx('GEO'))%rank = 3*natom
        !AllPrt(idx('NS'))%rank = depends on which isotopes
    end subroutine

  subroutine oneint(l, A)

!   radovan: frontend to routine which reads AO integrals
!            (typically from AOPROPER)
!            with the label l into matrix A
!            please don't use fixed length of l, 
!            rather (*) to be more flexible

!   ----------------------------------------------------------------------------
    character(*), intent(in)    :: l
    type(matrix), intent(inout) :: A
!   ----------------------------------------------------------------------------

    if (iszero(A)) call mat_init(A, A)
    call di_read_operator_int(l, A)

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
            lw = 50*n2+10000*nbas+5000000
        else if (what=='G' .and.nf==3*nd) then
            if (.not.present(a)) &
                call quit("prop_integ-twofck: Geometry Fock requires atom 'a' specified")
            lw = 50*n2+10000*nbas+5000000
        else if (what=='MM'.and.nf==6*nd) then
            lw = 50*n2+10000*nbas+5000000
        else if (what=='GG'.and.nf==3*3*nd) then
            if (.not.present(a).or..not.present(b)) &
                call quit("prop_integ-twock: 2nd ord. geom Fock reqires atoms 'a','b' specified")
            lw = (nd+2*nf+50)*n2+10000*nbas+5000000
        else
            call quit('prop_integ-twofck: wrong size(F) for what='//what)
        end if
        allocate(wrk(lw))
        do i=1,nd; wrk(1+n2*(i-1):n2*i) = D(i)%elms; end do
        !ajt Redundant: wrk(nbas**2+1:nbas**2*(nf+1)) = 0
        call GRCONT(wrk(1+n2*(nd+nf):lw),lw-n2*(nd+nf),wrk(1+n2*nd:n2*(nd+nf)) &
                   ,n2*nf,what(1:1)=='G',what(1:1)=='M',len(what),0 &
                   ,.false.,.true.,wrk(1:n2*nd),nd)
        do i=1,nf
            if (iszero(F(i))) call mat_init(F(i),F(i))
            F(i)%elms = wrk(1+n2*(nd+i-1):n2*(nd+i))
        end do
        deallocate(wrk)
    end subroutine


    subroutine twoctr(what,D,DD,e)
    !gm=T -> 3*NATOMS geometriske
    !gm=F -> tre magnetiske
        character(*),intent(in)    :: what
        type(matrix),intent(in)    :: D,DD
        real*8,      intent(inout) :: e(:)
        real*8,allocatable         :: wrk(:)
        if (nbas==0) stop 'prop_integ-twoctr: nbas not initialized'
        !ajt MagSus contraction doesn't work, so I disabled this. Use twofck instead
        if (what=='MM'.or.what=='MM') &
            call quit("prop_integ-twoctr: 'M' or 'MM' not yet implemented")
        if     (what=='G' .and.size(e)==3*natom) then
            allocate(wrk(50*nbas**2+10000*nbas+5000000))
        else if (what=='M' .and.size(e)==3) then
            allocate(wrk(50*nbas**2+10000*nbas+5000000))
        else if (what=='GG'.and.size(e)==3*3*natom**2) then
            allocate(wrk(50*nbas**2+10000*nbas+5000000))
        else if (what=='MM'.and.size(e)==3*3) then
            allocate(wrk(50*nbas**2+10000*nbas+5000000))
        else
            call quit('prop_integ-twoctr: wrong size(e) for what='//what)
        end if
        wrk(        1:  nbas**2) = D%elms !call mat_to_full(D ,1d0,A(:,:,1))
        wrk(nbas**2+1:2*nbas**2) = DD%elms !call mat_to_full(DD,1d0,A(:,:,2))
        e=0; call GRCONT(wrk(2*nbas**2+1:),(size(wrk)-2*nbas**2),e,size(e) &
                        ,what(1:1)=='G',what(1:1)=='M',len(what),0 &
                        ,.true.,.false.,wrk(1:2*nbas**2),2)
        deallocate(wrk)
    end subroutine


    function idx(p)
        character(*):: p; integer:: idx
        do idx=1,size(tab)
            if (tab(idx)%code==p) return
        end do
        print *,'Perturbation not found: '//p
        call quit('Perturbation not found')
    end function


    function prop_imag(p)
        character(*):: p; logical:: prop_imag
        prop_imag = tab(idx(p))%imag
    end function


    function prop_rank(p)
        character(*):: p; integer:: prop_rank
        prop_rank = tab(idx(p))%rank
        if (prop_rank==0) stop 'prop_rank: rank not initialized'
    end function


    subroutine prop_check(who,D0,n,p, c,w,nd,b,o, op,cc,ww,m,ne,zerg,zfck,zovl,imag)
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
            j=i; end if; end do
        end if
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
            end if
            if (present(ne).and.size(c)==n  ) ne = 1
            if (present(ne).and.size(c)==n+n) ne = product(c(n+1:n+n)-c(1:n)+1)
        else !ranges 1:p%comp
            if (present(cc)) cc = (/(1,i=1,n),(tab(ip(i))%rank,i=1,n)/)
            if (present(ne)) ne = product((/(tab(ip(i))%rank,i=1,n)/))
        end if
        !process optional w, produce ww
        if (present(w)) then
            if (size(w)/=n.and.size(w)/=n+1) &
                call QUIT(who//': w must be w(n) or w(n+1), with p(n)')
            if (present(ww)) ww(1:n+1) = (/w(1:n),0d0/)
            if (present(ww).and.size(w)==n+1) ww(n+1) = w(n+1)
        else if (present(ww)) then
            ww(1:n) = 0; ww(n+1) = 1
        end if
        !find order of p(i) according to idx
        if (present(op)) then
            do i=1,n
                op(i:i) = minloc(ip)
                ip(op(i)) = huge(1) !not found again by minloc
            end do
        end if
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
            end if
            if (present(o)) then
                if (size(o)/=m) call QUIT(who//': o has wrong size')
                if (any(1>o.or.o>m+n)) call QUIT(who//': o out of range')
                if (any((/((o(i)==o(j),i=j+1,m),j=1,m-1)/))) &
                    call QUIT(who//': repeated element in o')
            end if
            !knowing m, Fock contributions may be zero
            if (one.and.present(zfck).and.m>0) zfck = .true.
        else if (present(b).or.present(o)) then
            call QUIT(who//': if b or o present, so must nd and m be')
        end if
        !is the energy contribution
    end subroutine


    subroutine prop_order(n,op,cc,m,nd,ne, bb,oe, b,o)
        integer,intent(in):: n,op(n),cc(n+n),m,nd,ne
        integer,intent(out):: bb(m),oe(ne)
        integer,intent(in),optional:: b(:),o(:) !b(m),o(m) will break ifort
        integer:: kk(m+n,3),k(3),opi(n),i,j,l
        !cumulative product of dimensions (pitches)
        bb=nd-1; if (present(b)) bb = b
        kk(1:m+n,1) = (/bb,(cc(n+op(i))-cc(op(i))+1,i=1,n)/)
        kk(1,2)=1; do i=2,m+n
            kk(i,2) = kk(i-1,2)*kk(i-1,1); end do
        !reorder dimensions and pitches into output order
        do i=1,n; opi(op(i)) = i; end do !inverse permutation
        kk(m+1:m+n,1) = (/(kk(m+opi(i),1),i=1,n)/)
        kk(m+1:m+n,2) = (/(kk(m+opi(i),2),i=1,n)/)
        if (present(o)) then
            kk(1:m,3) = o(1:m)
            do i=m,1,-1
                k(1:1) = maxloc(kk(1:i,3))
                j = k(1); k = kk(j,:)
                kk(j:k(3)-1,:) = kk(j+1:k(3),:)
                kk(k(3),:) = k
            end do
        end if
        !construct oe from dimensions kk(:,1) and pitches kk(:,2)
        oe(1) = 1; j = 1
        do i=1,m+n
            oe(j+1:j*kk(i,1)) = (/(oe(1:j)+l*kk(i,2),l=1,kk(i,1)-1)/)
            j = j*kk(i,1)
        end do
    end subroutine


end module
