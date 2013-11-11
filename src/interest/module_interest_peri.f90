! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Library:      InteRest
!
! Module:       module_interest_peri.f90 
!
! Description:  InteRest library routine 
!
! Contains:     
!
! Licensing:    dirac_copyright_start
!                     Copyright (c) 2010 by the authors of DIRAC.
!                     All Rights Reserved.
!                     
!                     This source code is part of the DIRAC program package.
!                     It is provided under a written license and may be used,
!                     copied, transmitted, or stored only in accordance to the
!                     conditions of that written license.
!                     
!                     In particular, no part of the source code or compiled modules may
!                     be distributed outside the research group of the license holder.
!                     This means also that persons (e.g. post-docs) leaving the research
!                     group of the license holder may not take any part of Dirac,
!                     including modified files, with him/her, unless that person has
!                     obtained his/her own license.
!                     
!                     For information on how to get a license, as well as the
!                     author list and the complete list of contributors to the
!                     DIRAC program, see: http://dirac.chem.vu.nl
!               dirac_copyright_end
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Logs:         
!
! Comments:     
!                            
! ------------------------------------------------------------------------------------
MODULE module_interest_peri

  !> InteRest 
  use module_interest_eri

  implicit none

  !> public
  public estimate_dG_8x
  public process_dG_llll_8x_trs
  public process_dG_ssll_8x_trs
  public process_dG_ssss_8x_trs

  !> private
  private

CONTAINS


! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_dG_llll_8x_trs(ic,jc,kc,lc,io,jo,ko,lo,gout,dP,dG,nbas,doK,dPij,dGij,dPi,dPj,dGi,dGj,nbasr) 
    !> input
    logical,    intent(in   ) :: doK
    integer,    intent(in   ) :: nbas
    integer,    intent(in   ) :: nbasr
    integer,    intent(in   ) :: ic, jc, kc, lc
    integer,    intent(in   ) :: io, jo, ko, lo
    real(8),    intent(in   ) :: gout(kc,lc,ic,jc)
    complex(8), intent(in   ) :: dP(2,4,nbas,nbas)
    complex(8), intent(inout) :: dG(2,4,nbas,nbas)
    !> input/output coulomb
    real(8),    intent(in   ) :: dPij(5,ic,jc)
    real(8),    intent(inout) :: dGij(5,ic,jc)
    !> input/output exchange
    complex(8), intent(in   ) :: dPi(2,4,ic,nbasr)
    complex(8), intent(in   ) :: dPj(2,4,jc,nbasr)
    complex(8), intent(inout) :: dGi(2,4,ic,nbasr)
    complex(8), intent(inout) :: dGj(2,4,jc,nbasr)
    !> local
    real(8)    :: pij, pkl(441)
    real(8)    :: gij, gkl(441)
    integer    :: ij, i, j, kl, k, l, m
    integer    :: ibas, jbas, kbas, lbas
    complex(8) :: gg, paa_ll, pab_ll
    !>
    !> coulomb part
    !>
    kl=0
    do l=1,lc
      lbas=lo+l
      do k=1,kc
        kbas=ko+k
        kl=kl+1
        gkl(kl) = 0.0d0
        !> density matrix (Hermitian and time-reversal symmetric) 
        pkl(kl) = 4.0d0*(dreal(dP(1,1,kbas,lbas))) !Re(aa+bb)=2Re(aa)
      enddo
    enddo
    !>
    do j=1,jc
      do i=1,ic
        gij = 0.0d0 
        pij = dPij(1,i,j) 
        kl=0
        do l=1,lc
          do k=1,kc
            kl=kl+1
            !> Gij = [ij|kl]*Plk + [ij|lk]*Pkl
            !> Gji = [ji|kl]*Plk + [ji|lk]*Pkl -> done outside 
            gij = gij + gout(k,l,i,j)*pkl(kl)
            !> Gkl = [kl|ij]*Pji + [kl|ji]*Pij
            !> Glk = [lk|ij]*Pji + [lk|ji]*Pij
            gkl(kl) = gkl(kl) + gout(k,l,i,j)*pij
          enddo
        enddo
        dGij(1,i,j) = dGij(1,i,j) + gij
      enddo
    enddo
    !>
    if( ko>=lo )then
      kl=0
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          kl=kl+1
          dG(1,1,kbas,lbas) = dG(1,1,kbas,lbas) + dcmplx( gkl(kl), 0.0d0 )
        enddo
      enddo
    endif
    if( lo>=ko )then
      kl=0
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          kl=kl+1
          dG(1,1,lbas,kbas) = dG(1,1,lbas,kbas) + dcmplx( gkl(kl), 0.0d0 )
        enddo
      enddo
    endif
    !>
    !> exchange part
    !>
    if( .not. doK )return
    do l=1,lc
      lbas=lo+l
      do k=1,kc
        kbas=ko+k
        do j=1,jc
          do i=1,ic
            !> integrals processing: [ij|kl]
               !> large-large block 
               gg = dcmplx( gout(k,l,i,j), 0.0d0 )
            !> end
            !> contraction [ij|kl]*jk
               paa_ll = dPj(1,1,j,kbas)
               pab_ll = dPj(1,3,j,kbas)
               !> large-large block
               dGi(1,1,i,lbas) = dGi(1,1,i,lbas) - paa_ll*gg
               dGi(1,3,i,lbas) = dGi(1,3,i,lbas) - pab_ll*gg
            !> end
            !> contraction [ij|lk]*jl
               paa_ll = dPj(1,1,j,lbas)
               pab_ll = dPj(1,3,j,lbas)
               !> large-large block
               dGi(1,1,i,kbas) = dGi(1,1,i,kbas) - paa_ll*gg
               dGi(1,3,i,kbas) = dGi(1,3,i,kbas) - pab_ll*gg
            !> end
            !> contraction [ji|kl]*ik
               paa_ll = dPi(1,1,i,kbas)
               pab_ll = dPi(1,3,i,kbas)
               !> large-large block
               dGj(1,1,j,lbas) = dGj(1,1,j,lbas) - paa_ll*gg
               dGj(1,3,j,lbas) = dGj(1,3,j,lbas) - pab_ll*gg
            !> end
            !> contraction [ji|lk]*il
               paa_ll = dPi(1,1,i,lbas)
               pab_ll = dPi(1,3,i,lbas)
               !> large-large block
               dGj(1,1,j,kbas) = dGj(1,1,j,kbas) - paa_ll*gg
               dGj(1,3,j,kbas) = dGj(1,3,j,kbas) - pab_ll*gg
            !> end
          enddo
        enddo
      enddo
    enddo
  END SUBROUTINE

! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE process_dG_ssll_8x_trs(ic,jc,kc,lc,io,jo,ko,lo,gout,dP,dG,nbas,doK,dPij,dGij,dPi,dPj,dGi,dGj,nbasr) 
    !> input
    logical,    intent(in   ) :: doK
    integer,    intent(in   ) :: nbas
    integer,    intent(in   ) :: nbasr
    integer,    intent(in   ) :: ic, jc, kc, lc
    integer,    intent(in   ) :: io, jo, ko, lo
    real(8),    intent(in   ) :: gout(5,kc,lc,5,ic,jc)
    complex(8), intent(in   ) :: dP(2,4,nbas,nbas)
    complex(8), intent(inout) :: dG(2,4,nbas,nbas)
    !> input/output coulomb
    real(8),    intent(in   ) :: dPij(5,ic,jc)
    real(8),    intent(inout) :: dGij(5,ic,jc)
    !> input/output exchange
    complex(8), intent(in   ) :: dPi(2,4,ic,nbasr)
    complex(8), intent(in   ) :: dPj(2,4,jc,nbasr)
    complex(8), intent(inout) :: dGi(2,4,ic,nbasr)
    complex(8), intent(inout) :: dGj(2,4,jc,nbasr)
    !> local
    integer    :: ij, i, j, kl, k, l, m
    integer    :: ibas, jbas, kbas, lbas
    real(8)    :: pij(5), pkl(5,441)
    real(8)    :: gij(5), gkl(5,441)
    complex(8) :: gg
    complex(8) :: gaag, gabg, gbag, gbbg
    complex(8) :: ggaa, ggab, ggba, ggbb   
    complex(8) :: paa_ll, pab_ll, pba_ll, pbb_ll
    complex(8) :: paa_ls, pab_ls, pba_ls, pbb_ls
    complex(8) :: paa_sl, pab_sl, pba_sl, pbb_sl
    !>
    !> coulomb part
    !>
    kl=0
    do l=1,lc
      lbas=lo+l
      do k=1,kc
        kbas=ko+k
        kl=kl+1
        gkl(1,kl) = 0.0d0
        gkl(2,kl) = 0.0d0
        gkl(3,kl) = 0.0d0
        gkl(4,kl) = 0.0d0
        gkl(5,kl) = 0.0d0
        !> density matrix (Hermitian and time-reversal symmetric) 
        pkl(1,kl) = 4.0d0*(dreal(dP(1,1,kbas,lbas))) !Re(aa+bb)=2Re(aa)
        pkl(2,kl) = 4.0d0*(dreal(dP(2,2,kbas,lbas))) !Re(aa+bb)=2Re(aa)
        pkl(3,kl) = 4.0d0*(dimag(dP(2,2,kbas,lbas))) !Im(aa-bb)=2Im(aa)
        pkl(4,kl) = 4.0d0*(dreal(dP(2,4,kbas,lbas))) !Re(ab-ba)=2Re(ab)
        pkl(5,kl) = 4.0d0*(dimag(dP(2,4,kbas,lbas))) !Im(ab+ba)=2Im(ab)
      enddo
    enddo
    !>
    do j=1,jc
      do i=1,ic
        gij(1) = 0.0d0 
        gij(2) = 0.0d0 
        gij(3) = 0.0d0 
        gij(4) = 0.0d0 
        gij(5) = 0.0d0 
        pij(1) = dPij(1,i,j) 
        pij(2) = dPij(2,i,j) 
        pij(3) = dPij(3,i,j) 
        pij(4) = dPij(4,i,j) 
        pij(5) = dPij(5,i,j) 
        kl=0
        do l=1,lc
          do k=1,kc
            kl=kl+1
            !> Gij = [ij|kl]*Plk + [ij|lk]*Pkl
            !> Gji = [ji|kl]*Plk + [ji|lk]*Pkl -> done outside 
            gij(1) = gij(1) + gout(1,k,l,1,i,j)*pkl(1,kl) &
                            + gout(2,k,l,1,i,j)*pkl(2,kl) &
                            + gout(3,k,l,1,i,j)*pkl(3,kl) &
                            + gout(4,k,l,1,i,j)*pkl(4,kl) &
                            + gout(5,k,l,1,i,j)*pkl(5,kl)
            gij(2) = gij(2) + gout(1,k,l,2,i,j)*pkl(1,kl)
            gij(3) = gij(3) + gout(1,k,l,3,i,j)*pkl(1,kl)
            gij(4) = gij(4) + gout(1,k,l,4,i,j)*pkl(1,kl)
            gij(5) = gij(5) + gout(1,k,l,5,i,j)*pkl(1,kl)
            !> Gkl = [kl|ij]*Pji + [kl|ji]*Pij
            !> Glk = [lk|ij]*Pji + [lk|ji]*Pij
            gkl(1,kl) = gkl(1,kl) + gout(1,k,l,1,i,j)*pij(1) &
                                  + gout(1,k,l,2,i,j)*pij(2) &
                                  + gout(1,k,l,3,i,j)*pij(3) &
                                  + gout(1,k,l,4,i,j)*pij(4) &
                                  + gout(1,k,l,5,i,j)*pij(5)
            gkl(2,kl) = gkl(2,kl) + gout(2,k,l,1,i,j)*pij(1)
            gkl(3,kl) = gkl(3,kl) + gout(3,k,l,1,i,j)*pij(1)
            gkl(4,kl) = gkl(4,kl) + gout(4,k,l,1,i,j)*pij(1)
            gkl(5,kl) = gkl(5,kl) + gout(5,k,l,1,i,j)*pij(1)
          enddo
        enddo
        dGij(1,i,j) = dGij(1,i,j) + gij(1)
        dGij(2,i,j) = dGij(2,i,j) + gij(2)
        dGij(3,i,j) = dGij(3,i,j) + gij(3)
        dGij(4,i,j) = dGij(4,i,j) + gij(4)
        dGij(5,i,j) = dGij(5,i,j) + gij(5)
      enddo
    enddo
    !>
    if( ko>=lo )then
      kl=0
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          kl=kl+1
          dG(1,1,kbas,lbas) = dG(1,1,kbas,lbas) + dcmplx( gkl(1,kl), 0.0d0     )
          dG(2,2,kbas,lbas) = dG(2,2,kbas,lbas) + dcmplx( gkl(2,kl), gkl(3,kl) )
          dG(2,4,kbas,lbas) = dG(2,4,kbas,lbas) + dcmplx( gkl(4,kl), gkl(5,kl) )
        enddo
      enddo
    endif
    if( lo>=ko )then
      kl=0
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          kl=kl+1
          dG(1,1,lbas,kbas) = dG(1,1,lbas,kbas) + dcmplx( gkl(1,kl), 0.0d0     )
          dG(2,2,lbas,kbas) = dG(2,2,lbas,kbas) + dcmplx( gkl(2,kl),-gkl(3,kl) )
          dG(2,4,lbas,kbas) = dG(2,4,lbas,kbas) + dcmplx(-gkl(4,kl),-gkl(5,kl) )
        enddo
      enddo
    endif
    !>
    !> exchange part
    !>
    if( .not. doK )return
    do l=1,lc
      lbas=lo+l
      do k=1,kc
        kbas=ko+k
        do j=1,jc
          do i=1,ic
            !> integrals processing: [ij|kl]
               !> large-large block 
               gg = dcmplx( gout(1,k,l,1,i,j), 0.0d0 )
               !> large-small block 
               ggaa = dcmplx( gout(2,k,l,1,i,j), gout(3,k,l,1,i,j) )
               ggab = dcmplx( gout(4,k,l,1,i,j), gout(5,k,l,1,i,j) )
               ggbb = dconjg( ggaa )
               ggba =-dconjg( ggab )
               !> small-large block 
               gaag = dcmplx( gout(1,k,l,2,i,j), gout(1,k,l,3,i,j) )
               gabg = dcmplx( gout(1,k,l,4,i,j), gout(1,k,l,5,i,j) )
               gbbg = dconjg( gaag )
               gbag =-dconjg( gabg ) 
            !> end
            !> contraction [ij|kl]*jk
               paa_ll = dPj(1,1,j,kbas)
               paa_sl = dPj(2,1,j,kbas)
               paa_ls = dPj(1,2,j,kbas)
               pab_ll = dPj(1,3,j,kbas)
               pab_sl = dPj(2,3,j,kbas)
               pab_ls = dPj(1,4,j,kbas)
               pba_ll =-dconjg(pab_ll )
               pba_sl =-dconjg(pab_sl )
               pba_ls =-dconjg(pab_ls )
               pbb_ll = dconjg(paa_ll )
               pbb_sl = dconjg(paa_sl )
               pbb_ls = dconjg(paa_ls )
               !> large-large block
               dGi(1,1,i,lbas) = dGi(1,1,i,lbas) - paa_ll*gg
               dGi(1,3,i,lbas) = dGi(1,3,i,lbas) - pab_ll*gg
               !> large-small block
               dGi(1,2,i,lbas) = dGi(1,2,i,lbas) - paa_ls*ggaa - pab_ls*ggba  
               dGi(1,4,i,lbas) = dGi(1,4,i,lbas) - paa_ls*ggab - pab_ls*ggbb
               !> small-large block
               dGi(2,1,i,lbas) = dGi(2,1,i,lbas) - gaag*paa_sl - gabg*pba_sl 
               dGi(2,3,i,lbas) = dGi(2,3,i,lbas) - gaag*pab_sl - gabg*pbb_sl 
            !> end
            !> contraction [ij|lk]*jl
               paa_ll = dPj(1,1,j,lbas)
               paa_sl = dPj(2,1,j,lbas)
               paa_ls = dPj(1,2,j,lbas)
               pab_ll = dPj(1,3,j,lbas)
               pab_sl = dPj(2,3,j,lbas)
               pab_ls = dPj(1,4,j,lbas)
               pba_ll =-dconjg(pab_ll )
               pba_sl =-dconjg(pab_sl )
               pba_ls =-dconjg(pab_ls )
               pbb_ll = dconjg(paa_ll )
               pbb_sl = dconjg(paa_sl )
               pbb_ls = dconjg(paa_ls )
               !> large-large block
               dGi(1,1,i,kbas) = dGi(1,1,i,kbas) - paa_ll*gg
               dGi(1,3,i,kbas) = dGi(1,3,i,kbas) - pab_ll*gg
               !> large-small block
               dGi(1,2,i,kbas) = dGi(1,2,i,kbas) - paa_ls*ggbb + pab_ls*ggba  
               dGi(1,4,i,kbas) = dGi(1,4,i,kbas) + paa_ls*ggab - pab_ls*ggaa
               !> small-large block
               dGi(2,1,i,kbas) = dGi(2,1,i,kbas) - gaag*paa_sl - gabg*pba_sl 
               dGi(2,3,i,kbas) = dGi(2,3,i,kbas) - gaag*pab_sl - gabg*pbb_sl 
            !> end
            !> contraction [ji|kl]*ik
               paa_ll = dPi(1,1,i,kbas)
               paa_sl = dPi(2,1,i,kbas)
               paa_ls = dPi(1,2,i,kbas)
               pab_ll = dPi(1,3,i,kbas)
               pab_sl = dPi(2,3,i,kbas)
               pab_ls = dPi(1,4,i,kbas)
               pba_ll =-dconjg(pab_ll )
               pba_sl =-dconjg(pab_sl )
               pba_ls =-dconjg(pab_ls )
               pbb_ll = dconjg(paa_ll )
               pbb_sl = dconjg(paa_sl )
               pbb_ls = dconjg(paa_ls )
               !> large-large block
               dGj(1,1,j,lbas) = dGj(1,1,j,lbas) - paa_ll*gg
               dGj(1,3,j,lbas) = dGj(1,3,j,lbas) - pab_ll*gg
               !> large-small block
               dGj(1,2,j,lbas) = dGj(1,2,j,lbas) - paa_ls*ggaa - pab_ls*ggba  
               dGj(1,4,j,lbas) = dGj(1,4,j,lbas) - paa_ls*ggab - pab_ls*ggbb
               !> small-large block
               dGj(2,1,j,lbas) = dGj(2,1,j,lbas) - gbbg*paa_sl + gabg*pba_sl 
               dGj(2,3,j,lbas) = dGj(2,3,j,lbas) - gbbg*pab_sl + gabg*pbb_sl 
            !> end
            !> contraction [ji|lk]*il
               paa_ll = dPi(1,1,i,lbas)
               paa_sl = dPi(2,1,i,lbas)
               paa_ls = dPi(1,2,i,lbas)
               pab_ll = dPi(1,3,i,lbas)
               pab_sl = dPi(2,3,i,lbas)
               pab_ls = dPi(1,4,i,lbas)
               pba_ll =-dconjg(pab_ll )
               pba_sl =-dconjg(pab_sl )
               pba_ls =-dconjg(pab_ls )
               pbb_ll = dconjg(paa_ll )
               pbb_sl = dconjg(paa_sl )
               pbb_ls = dconjg(paa_ls )
               !> large-large block
               dGj(1,1,j,kbas) = dGj(1,1,j,kbas) - paa_ll*gg
               dGj(1,3,j,kbas) = dGj(1,3,j,kbas) - pab_ll*gg
               !> large-small block
               dGj(1,2,j,kbas) = dGj(1,2,j,kbas) - paa_ls*ggbb + pab_ls*ggba  
               dGj(1,4,j,kbas) = dGj(1,4,j,kbas) + paa_ls*ggab - pab_ls*ggaa
               !> small-large block
               dGj(2,1,j,kbas) = dGj(2,1,j,kbas) - gbbg*paa_sl + gabg*pba_sl 
               dGj(2,3,j,kbas) = dGj(2,3,j,kbas) - gbbg*pab_sl + gabg*pbb_sl 
            !> end
          enddo
        enddo
      enddo
    enddo
  END SUBROUTINE


! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE process_dG_ssss_8x_trs(ic,jc,kc,lc,io,jo,ko,lo,gout,dP,dG,nbas,doK,dPij,dGij,dPi,dPj,dGi,dGj,nbasr) 
    !> input
    logical,    intent(in   ) :: doK
    integer,    intent(in   ) :: nbas
    integer,    intent(in   ) :: nbasr
    integer,    intent(in   ) :: ic, jc, kc, lc
    integer,    intent(in   ) :: io, jo, ko, lo
    real(8),    intent(in   ) :: gout(5,kc,lc,5,ic,jc)
    complex(8), intent(in   ) :: dP(2,4,nbas,nbas)
    complex(8), intent(inout) :: dG(2,4,nbas,nbas)
    !> input/output coulomb
    real(8),    intent(in   ) :: dPij(5,ic,jc)
    real(8),    intent(inout) :: dGij(5,ic,jc)
    !> input/output exchange
    complex(8), intent(in   ) :: dPi(2,4,ic,nbasr)
    complex(8), intent(in   ) :: dPj(2,4,jc,nbasr)
    complex(8), intent(inout) :: dGi(2,4,ic,nbasr)
    complex(8), intent(inout) :: dGj(2,4,jc,nbasr)
    !> local
    integer    :: ij, i, j, kl, k, l, m
    integer    :: ibas, jbas, kbas, lbas
    real(8)    :: pij(5), pkl(5,441)
    real(8)    :: gij(5), gkl(5,441)
    complex(8) :: gg
    complex(8) :: gaag, gabg, gbag, gbbg
    complex(8) :: ggaa, ggab, ggba, ggbb   
    complex(8) :: gaagaa, gabgaa, gbagaa, gbbgaa
    complex(8) :: gaagab, gabgab, gbagab, gbbgab
    complex(8) :: gaagba, gabgba, gbagba, gbbgba
    complex(8) :: gaagbb, gabgbb, gbagbb, gbbgbb
    complex(8) :: paa_ll, pab_ll, pba_ll, pbb_ll
    complex(8) :: paa_ls, pab_ls, pba_ls, pbb_ls
    complex(8) :: paa_sl, pab_sl, pba_sl, pbb_sl
    complex(8) :: paa_ss, pab_ss, pba_ss, pbb_ss
    !>
    !> coulomb part
    !>
    kl=0
    do l=1,lc
      lbas=lo+l
      do k=1,kc
        kbas=ko+k
        kl=kl+1
        gkl(1,kl) = 0.0d0
        gkl(2,kl) = 0.0d0
        gkl(3,kl) = 0.0d0
        gkl(4,kl) = 0.0d0
        gkl(5,kl) = 0.0d0
        !> density matrix (Hermitian and time-reversal symmetric) 
        pkl(1,kl) = 4.0d0*(dreal(dP(1,1,kbas,lbas))) !Re(aa+bb)=2Re(aa)
        pkl(2,kl) = 4.0d0*(dreal(dP(2,2,kbas,lbas))) !Re(aa+bb)=2Re(aa)
        pkl(3,kl) = 4.0d0*(dimag(dP(2,2,kbas,lbas))) !Im(aa-bb)=2Im(aa)
        pkl(4,kl) = 4.0d0*(dreal(dP(2,4,kbas,lbas))) !Re(ab-ba)=2Re(ab)
        pkl(5,kl) = 4.0d0*(dimag(dP(2,4,kbas,lbas))) !Im(ab+ba)=2Im(ab)
      enddo
    enddo
    !>
    do j=1,jc
      do i=1,ic
        gij(1) = 0.0d0 
        gij(2) = 0.0d0 
        gij(3) = 0.0d0 
        gij(4) = 0.0d0 
        gij(5) = 0.0d0 
        pij(1) = dPij(1,i,j) 
        pij(2) = dPij(2,i,j) 
        pij(3) = dPij(3,i,j) 
        pij(4) = dPij(4,i,j) 
        pij(5) = dPij(5,i,j) 
        do m=1,5  
          kl=0
          do l=1,lc
            do k=1,kc
              kl=kl+1
              !> Gij = [ij|kl]*Plk + [ij|lk]*Pkl
              !> Gji = [ji|kl]*Plk + [ji|lk]*Pkl
              gij(m) = gij(m) + gout(1,k,l,m,i,j)*pkl(1,kl) &
                              + gout(2,k,l,m,i,j)*pkl(2,kl) &
                              + gout(3,k,l,m,i,j)*pkl(3,kl) &
                              + gout(4,k,l,m,i,j)*pkl(4,kl) &
                              + gout(5,k,l,m,i,j)*pkl(5,kl)
              !> Gkl = [kl|ij]*Pji + [kl|ji]*Pij
              !> Glk = [lk|ij]*Pji + [lk|ji]*Pij
              gkl(1,kl) = gkl(1,kl) + gout(1,k,l,m,i,j)*pij(m)
              gkl(2,kl) = gkl(2,kl) + gout(2,k,l,m,i,j)*pij(m)
              gkl(3,kl) = gkl(3,kl) + gout(3,k,l,m,i,j)*pij(m)
              gkl(4,kl) = gkl(4,kl) + gout(4,k,l,m,i,j)*pij(m)
              gkl(5,kl) = gkl(5,kl) + gout(5,k,l,m,i,j)*pij(m)
            enddo
          enddo
        enddo
        dGij(1,i,j) = dGij(1,i,j) + gij(1)
        dGij(2,i,j) = dGij(2,i,j) + gij(2)
        dGij(3,i,j) = dGij(3,i,j) + gij(3)
        dGij(4,i,j) = dGij(4,i,j) + gij(4)
        dGij(5,i,j) = dGij(5,i,j) + gij(5)
      enddo
    enddo
    !>
    if( ko>=lo )then
      kl=0
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          kl=kl+1
          dG(1,1,kbas,lbas) = dG(1,1,kbas,lbas) + dcmplx( gkl(1,kl), 0.0d0     )
          dG(2,2,kbas,lbas) = dG(2,2,kbas,lbas) + dcmplx( gkl(2,kl), gkl(3,kl) )
          dG(2,4,kbas,lbas) = dG(2,4,kbas,lbas) + dcmplx( gkl(4,kl), gkl(5,kl) )
        enddo
      enddo
    endif
    if( lo>=ko )then
      kl=0
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          kl=kl+1
          dG(1,1,lbas,kbas) = dG(1,1,lbas,kbas) + dcmplx( gkl(1,kl), 0.0d0     )
          dG(2,2,lbas,kbas) = dG(2,2,lbas,kbas) + dcmplx( gkl(2,kl),-gkl(3,kl) )
          dG(2,4,lbas,kbas) = dG(2,4,lbas,kbas) + dcmplx(-gkl(4,kl),-gkl(5,kl) )
        enddo
      enddo
    endif
    !>
    !> exchange part
    !>
    if( .not. doK )return
    do l=1,lc
      lbas=lo+l
      do k=1,kc
        kbas=ko+k
        do j=1,jc
          do i=1,ic
            !> integrals processing: [ij|kl]
               !> large-large block 
               gg = dcmplx( gout(1,k,l,1,i,j), 0.0d0 )
               !> large-small block 
               ggaa = dcmplx( gout(2,k,l,1,i,j), gout(3,k,l,1,i,j) )
               ggab = dcmplx( gout(4,k,l,1,i,j), gout(5,k,l,1,i,j) )
               ggbb = dconjg( ggaa )
               ggba =-dconjg( ggab )
               !> small-large block 
               gaag = dcmplx( gout(1,k,l,2,i,j), gout(1,k,l,3,i,j) )
               gabg = dcmplx( gout(1,k,l,4,i,j), gout(1,k,l,5,i,j) )
               gbbg = dconjg( gaag )
               gbag =-dconjg( gabg ) 
               !> small-small block: aa
               gaagaa = dcmplx( ( gout(2,k,l,2,i,j)-gout(3,k,l,3,i,j)), ( gout(3,k,l,2,i,j)+gout(2,k,l,3,i,j)) )
               gaagba = dcmplx( (-gout(4,k,l,2,i,j)-gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)-gout(4,k,l,3,i,j)) )
               gabgaa = dcmplx( ( gout(2,k,l,4,i,j)-gout(3,k,l,5,i,j)), ( gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
               gabgba = dcmplx( (-gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), ( gout(5,k,l,4,i,j)-gout(4,k,l,5,i,j)) )
               !> small-small block: bb
               gbbgbb = dconjg( gaagaa )
               gbbgab =-dconjg( gaagba ) 
               gbagbb =-dconjg( gabgaa )
               gbagab = dconjg( gabgba )
               !> small-small block: ab
               gaagab = dcmplx( ( gout(4,k,l,2,i,j)-gout(5,k,l,3,i,j)), ( gout(5,k,l,2,i,j)+gout(4,k,l,3,i,j)) )
               gaagbb = dcmplx( ( gout(2,k,l,2,i,j)+gout(3,k,l,3,i,j)), (-gout(3,k,l,2,i,j)+gout(2,k,l,3,i,j)) )
               gabgab = dcmplx( ( gout(4,k,l,4,i,j)-gout(5,k,l,5,i,j)), ( gout(5,k,l,4,i,j)+gout(4,k,l,5,i,j)) )
               gabgbb = dcmplx( ( gout(2,k,l,4,i,j)+gout(3,k,l,5,i,j)), (-gout(3,k,l,4,i,j)+gout(2,k,l,5,i,j)) )
               !> small-small block: ba
               gbbgba =-dconjg( gaagab )
               gbbgaa = dconjg( gaagbb )
               gbagba = dconjg( gabgab )
               gbagaa =-dconjg( gabgbb )
            !> end
            !> contraction [ij|kl]*jk
               paa_ll = dPj(1,1,j,kbas)
               paa_sl = dPj(2,1,j,kbas)
               paa_ls = dPj(1,2,j,kbas)
               paa_ss = dPj(2,2,j,kbas)
               pab_ll = dPj(1,3,j,kbas)
               pab_sl = dPj(2,3,j,kbas)
               pab_ls = dPj(1,4,j,kbas)
               pab_ss = dPj(2,4,j,kbas)
               pba_ll =-dconjg(pab_ll )
               pba_sl =-dconjg(pab_sl )
               pba_ls =-dconjg(pab_ls )
               pba_ss =-dconjg(pab_ss )
               pbb_ll = dconjg(paa_ll )
               pbb_sl = dconjg(paa_sl )
               pbb_ls = dconjg(paa_ls )
               pbb_ss = dconjg(paa_ss )
               !> large-large block
               dGi(1,1,i,lbas) = dGi(1,1,i,lbas) - paa_ll*gg
               dGi(1,3,i,lbas) = dGi(1,3,i,lbas) - pab_ll*gg
               !> large-small block
               dGi(1,2,i,lbas) = dGi(1,2,i,lbas) - paa_ls*ggaa - pab_ls*ggba  
               dGi(1,4,i,lbas) = dGi(1,4,i,lbas) - paa_ls*ggab - pab_ls*ggbb
               !> small-large block
               dGi(2,1,i,lbas) = dGi(2,1,i,lbas) - gaag*paa_sl - gabg*pba_sl 
               dGi(2,3,i,lbas) = dGi(2,3,i,lbas) - gaag*pab_sl - gabg*pbb_sl 
               !> small-small block
               dGi(2,2,i,lbas) = dGi(2,2,i,lbas) - gaagaa*paa_ss - gaagba*pab_ss - gabgaa*pba_ss - gabgba*pbb_ss 
               dGi(2,4,i,lbas) = dGi(2,4,i,lbas) - gaagab*paa_ss - gaagbb*pab_ss - gabgab*pba_ss - gabgbb*pbb_ss 
            !> end
            !> contraction [ij|lk]*jl
               paa_ll = dPj(1,1,j,lbas)
               paa_sl = dPj(2,1,j,lbas)
               paa_ls = dPj(1,2,j,lbas)
               paa_ss = dPj(2,2,j,lbas)
               pab_ll = dPj(1,3,j,lbas)
               pab_sl = dPj(2,3,j,lbas)
               pab_ls = dPj(1,4,j,lbas)
               pab_ss = dPj(2,4,j,lbas)
               pba_ll =-dconjg(pab_ll )
               pba_sl =-dconjg(pab_sl )
               pba_ls =-dconjg(pab_ls )
               pba_ss =-dconjg(pab_ss )
               pbb_ll = dconjg(paa_ll )
               pbb_sl = dconjg(paa_sl )
               pbb_ls = dconjg(paa_ls )
               pbb_ss = dconjg(paa_ss )
               !> large-large block
               dGi(1,1,i,kbas) = dGi(1,1,i,kbas) - paa_ll*gg
               dGi(1,3,i,kbas) = dGi(1,3,i,kbas) - pab_ll*gg
               !> large-small block
               dGi(1,2,i,kbas) = dGi(1,2,i,kbas) - paa_ls*ggbb + pab_ls*ggba  
               dGi(1,4,i,kbas) = dGi(1,4,i,kbas) + paa_ls*ggab - pab_ls*ggaa
               !> small-large block
               dGi(2,1,i,kbas) = dGi(2,1,i,kbas) - gaag*paa_sl - gabg*pba_sl 
               dGi(2,3,i,kbas) = dGi(2,3,i,kbas) - gaag*pab_sl - gabg*pbb_sl 
               !> small-small block
               dGi(2,2,i,kbas) = dGi(2,2,i,kbas) - gaagbb*paa_ss + gaagba*pab_ss - gabgbb*pba_ss + gabgba*pbb_ss 
               dGi(2,4,i,kbas) = dGi(2,4,i,kbas) + gaagab*paa_ss - gaagaa*pab_ss + gabgab*pba_ss - gabgaa*pbb_ss 
            !> end
            !> contraction [ji|kl]*ik
               paa_ll = dPi(1,1,i,kbas)
               paa_sl = dPi(2,1,i,kbas)
               paa_ls = dPi(1,2,i,kbas)
               paa_ss = dPi(2,2,i,kbas)
               pab_ll = dPi(1,3,i,kbas)
               pab_sl = dPi(2,3,i,kbas)
               pab_ls = dPi(1,4,i,kbas)
               pab_ss = dPi(2,4,i,kbas)
               pba_ll =-dconjg(pab_ll )
               pba_sl =-dconjg(pab_sl )
               pba_ls =-dconjg(pab_ls )
               pba_ss =-dconjg(pab_ss )
               pbb_ll = dconjg(paa_ll )
               pbb_sl = dconjg(paa_sl )
               pbb_ls = dconjg(paa_ls )
               pbb_ss = dconjg(paa_ss )
               !> large-large block
               dGj(1,1,j,lbas) = dGj(1,1,j,lbas) - paa_ll*gg
               dGj(1,3,j,lbas) = dGj(1,3,j,lbas) - pab_ll*gg
               !> large-small block
               dGj(1,2,j,lbas) = dGj(1,2,j,lbas) - paa_ls*ggaa - pab_ls*ggba  
               dGj(1,4,j,lbas) = dGj(1,4,j,lbas) - paa_ls*ggab - pab_ls*ggbb
               !> small-large block
               dGj(2,1,j,lbas) = dGj(2,1,j,lbas) - gbbg*paa_sl + gabg*pba_sl 
               dGj(2,3,j,lbas) = dGj(2,3,j,lbas) - gbbg*pab_sl + gabg*pbb_sl 
               !> small-small block
               dGj(2,2,j,lbas) = dGj(2,2,j,lbas) - gbbgaa*paa_ss - gbbgba*pab_ss + gabgaa*pba_ss + gabgba*pbb_ss 
               dGj(2,4,j,lbas) = dGj(2,4,j,lbas) - gbbgab*paa_ss - gbbgbb*pab_ss + gabgab*pba_ss + gabgbb*pbb_ss 
            !> end
            !> contraction [ji|lk]*il
               paa_ll = dPi(1,1,i,lbas)
               paa_sl = dPi(2,1,i,lbas)
               paa_ls = dPi(1,2,i,lbas)
               paa_ss = dPi(2,2,i,lbas)
               pab_ll = dPi(1,3,i,lbas)
               pab_sl = dPi(2,3,i,lbas)
               pab_ls = dPi(1,4,i,lbas)
               pab_ss = dPi(2,4,i,lbas)
               pba_ll =-dconjg(pab_ll )
               pba_sl =-dconjg(pab_sl )
               pba_ls =-dconjg(pab_ls )
               pba_ss =-dconjg(pab_ss )
               pbb_ll = dconjg(paa_ll )
               pbb_sl = dconjg(paa_sl )
               pbb_ls = dconjg(paa_ls )
               pbb_ss = dconjg(paa_ss )
               !> large-large block
               dGj(1,1,j,kbas) = dGj(1,1,j,kbas) - paa_ll*gg
               dGj(1,3,j,kbas) = dGj(1,3,j,kbas) - pab_ll*gg
               !> large-small block
               dGj(1,2,j,kbas) = dGj(1,2,j,kbas) - paa_ls*ggbb + pab_ls*ggba  
               dGj(1,4,j,kbas) = dGj(1,4,j,kbas) + paa_ls*ggab - pab_ls*ggaa
               !> small-large block
               dGj(2,1,j,kbas) = dGj(2,1,j,kbas) - gbbg*paa_sl + gabg*pba_sl 
               dGj(2,3,j,kbas) = dGj(2,3,j,kbas) - gbbg*pab_sl + gabg*pbb_sl 
               !> small-small block
               dGj(2,2,j,kbas) = dGj(2,2,j,kbas) - gbbgbb*paa_ss + gbbgba*pab_ss + gabgbb*pba_ss - gabgba*pbb_ss 
               dGj(2,4,j,kbas) = dGj(2,4,j,kbas) + gbbgab*paa_ss - gbbgaa*pab_ss - gabgab*pba_ss + gabgaa*pbb_ss 
            !> end
          enddo
        enddo
      enddo
    enddo
  END SUBROUTINE


! ------------------------------------------------------------------------------------
!> note: Kmax sc-elements are already renormalized on output 
  SUBROUTINE estimate_dG_8x(class,i,j,k,l,ci,cj,ck,cl,Kmax,pmax,gmax,nshl) 
    !> input
    character(len=4), intent(in) :: class
    real(8),          intent(in) :: ci,cj,ck,cl
    integer,          intent(in) :: i,j,k,l,nshl
    real(8),          intent(in) :: gmax(  5,nshl,nshl)
    real(8),          intent(in) :: pmax(4,4,nshl,nshl)
    !> output
    real(8),          intent(out) :: Kmax(2,2)
    !> local 
    real(8) :: gg
    real(8) :: gaag, gabg, gbag, gbbg
    real(8) :: ggaa, ggab, ggba, ggbb   
    real(8) :: gaagaa, gabgaa, gbagaa, gbbgaa
    real(8) :: gaagab, gabgab, gbagab, gbbgab
    real(8) :: gaagba, gabgba, gbagba, gbbgba
    real(8) :: gaagbb, gabgbb, gbagbb, gbbgbb
    real(8) :: paa_ll, pab_ll, pba_ll, pbb_ll
    real(8) :: paa_ls, pab_ls, pba_ls, pbb_ls
    real(8) :: paa_sl, pab_sl, pba_sl, pbb_sl
    real(8) :: paa_ss, pab_ss, pba_ss, pbb_ss
    real(8) :: dGil(4,4), dGik(4,4), dGjl(4,4), dGjk(4,4)
    !>
    select case(class)
      case('llll')
        !> integrals processing: [ij|kl]
           !> large-large block 
           gg = gmax(1,i,j)*gmax(1,k,l) ![A|C]
        !> contraction [ij|kl]*jk
           paa_ll = pmax(1,1,j,k)
           pba_ll = pmax(3,1,j,k)
           pab_ll = pmax(1,3,j,k)
           pbb_ll = pmax(3,3,j,k)
           !> large-large block
           dGil(1,1) = paa_ll*gg
           dGil(1,3) = pab_ll*gg
           dGil(3,1) = pba_ll*gg
           dGil(3,3) = pbb_ll*gg
        !> end
        !> contraction [ij|lk]*jl
           paa_ll = pmax(1,1,j,l)
           pba_ll = pmax(3,1,j,l)
           pab_ll = pmax(1,3,j,l)
           pbb_ll = pmax(3,3,j,l)
           !> large-large block
           dGik(1,1) = paa_ll*gg
           dGik(1,3) = pab_ll*gg
           dGik(3,1) = pba_ll*gg
           dGik(3,3) = pbb_ll*gg
        !> end
        !> contraction [ji|kl]*ik
           paa_ll = pmax(1,1,i,k)
           pba_ll = pmax(3,1,i,k)
           pab_ll = pmax(1,3,i,k)
           pbb_ll = pmax(3,3,i,k)
           !> large-large block
           dGjl(1,1) = paa_ll*gg
           dGjl(1,3) = pab_ll*gg
           dGjl(3,1) = pba_ll*gg
           dGjl(3,3) = pbb_ll*gg
        !> end
        !> contraction [ji|lk]*il
           paa_ll = pmax(1,1,i,l)
           pba_ll = pmax(3,1,i,l)
           pab_ll = pmax(1,3,i,l)
           pbb_ll = pmax(3,3,i,l)
           !> large-large block
           dGjk(1,1) = paa_ll*gg
           dGjk(1,3) = pab_ll*gg
           dGjk(3,1) = pba_ll*gg
           dGjk(3,3) = pbb_ll*gg
        !> end
        !> K-estimate
        Kmax(1,1) = dmax1(dGil(1,1),dGik(1,1),dGjl(1,1),dGjk(1,1),&
                          dGil(3,1),dGik(3,1),dGjl(3,1),dGjk(3,1),&
                          dGil(1,3),dGik(1,3),dGjl(1,3),dGjk(1,3),&
                          dGil(3,3),dGik(3,3),dGjl(3,3),dGjk(3,3) )
        Kmax(2,1) = 0.0d0
        Kmax(1,2) = 0.0d0
        Kmax(2,2) = 0.0d0
        !>
        !>
        !>
      case('ssll','slaa','slcc')
        !> integrals processing: [ij|kl]
           !> large-large block 
           gg = gmax(1,i,j)*gmax(1,k,l) ![A|C]
           !> large-small block 
           ggaa = gmax(1,i,j)*(gmax(2,k,l)+gmax(3,k,l)) ![A|c]
           ggab = gmax(1,i,j)*(gmax(4,k,l)+gmax(5,k,l)) ![A|d]
           ggbb = ggaa
           ggba = ggab
           !> small-large block 
           gaag = gmax(1,k,l)*(gmax(2,i,j)+gmax(3,i,j)) ![a|C]
           gabg = gmax(1,k,l)*(gmax(4,i,j)+gmax(5,i,j)) ![b|C]
           gbbg = gaag
           gbag = gabg
        !> end
        !> contraction [ij|kl]*jk
           paa_ll = pmax(1,1,j,k)
           pba_ll = pmax(3,1,j,k)
           pab_ll = pmax(1,3,j,k)
           pbb_ll = pmax(3,3,j,k)
           paa_ls = pmax(1,2,j,k)
           pba_ls = pmax(3,2,j,k)
           pab_ls = pmax(1,4,j,k)
           pbb_ls = pmax(3,4,j,k)
           paa_sl = pmax(2,1,j,k)
           pba_sl = pmax(4,1,j,k)
           pab_sl = pmax(2,3,j,k)
           pbb_sl = pmax(4,3,j,k)
           !> large-large block
           dGil(1,1) = paa_ll*gg
           dGil(1,3) = pab_ll*gg
           dGil(3,1) = pba_ll*gg
           dGil(3,3) = pbb_ll*gg
           !> large-small block
           dGil(1,2) = (paa_ls*ggaa + pab_ls*ggba)*cl  
           dGil(1,4) = (paa_ls*ggab + pab_ls*ggbb)*cl  
           dGil(3,2) = (pba_ls*ggaa + pbb_ls*ggba)*cl  
           dGil(3,4) = (pba_ls*ggab + pbb_ls*ggbb)*cl  
           !> small-large block
           dGil(2,1) = (gaag*paa_sl + gabg*pba_sl)*ci 
           dGil(2,3) = (gaag*pab_sl + gabg*pbb_sl)*ci 
           dGil(4,1) = (gbag*paa_sl + gbbg*pba_sl)*ci 
           dGil(4,3) = (gbag*pab_sl + gbbg*pbb_sl)*ci 
        !> end
        !> contraction [ij|lk]*jl
           paa_ll = pmax(1,1,j,l)
           pba_ll = pmax(3,1,j,l)
           pab_ll = pmax(1,3,j,l)
           pbb_ll = pmax(3,3,j,l)
           paa_ls = pmax(1,2,j,l)
           pba_ls = pmax(3,2,j,l)
           pab_ls = pmax(1,4,j,l)
           pbb_ls = pmax(3,4,j,l)
           paa_sl = pmax(2,1,j,l)
           pba_sl = pmax(4,1,j,l)
           pab_sl = pmax(2,3,j,l)
           pbb_sl = pmax(4,3,j,l)
           !> large-large block
           dGik(1,1) = paa_ll*gg
           dGik(1,3) = pab_ll*gg
           dGik(3,1) = pba_ll*gg
           dGik(3,3) = pbb_ll*gg
           !> large-small block
           dGik(1,2) = (paa_ls*ggbb + pab_ls*ggba)*ck  
           dGik(1,4) = (paa_ls*ggab + pab_ls*ggaa)*ck
           dGik(3,2) = (pba_ls*ggbb + pbb_ls*ggba)*ck
           dGik(3,4) = (pba_ls*ggab + pbb_ls*ggaa)*ck
           !> small-large block
           dGik(2,1) = (gaag*paa_sl + gabg*pba_sl)*ci 
           dGik(2,3) = (gaag*pab_sl + gabg*pbb_sl)*ci 
           dGik(4,1) = (gbag*paa_sl + gbbg*pba_sl)*ci 
           dGik(4,3) = (gbag*pab_sl + gbbg*pbb_sl)*ci 
        !> end
        !> contraction [ji|kl]*ik
           paa_ll = pmax(1,1,i,k)
           pba_ll = pmax(3,1,i,k)
           pab_ll = pmax(1,3,i,k)
           pbb_ll = pmax(3,3,i,k)
           paa_ls = pmax(1,2,i,k)
           pba_ls = pmax(3,2,i,k)
           pab_ls = pmax(1,4,i,k)
           pbb_ls = pmax(3,4,i,k)
           paa_sl = pmax(2,1,i,k)
           pba_sl = pmax(4,1,i,k)
           pab_sl = pmax(2,3,i,k)
           pbb_sl = pmax(4,3,i,k)
           !> large-large block
           dGjl(1,1) = paa_ll*gg
           dGjl(1,3) = pab_ll*gg
           dGjl(3,1) = pba_ll*gg
           dGjl(3,3) = pbb_ll*gg
           !> large-small block
           dGjl(1,2) = (paa_ls*ggaa + pab_ls*ggba)*cl  
           dGjl(1,4) = (paa_ls*ggab + pab_ls*ggbb)*cl
           dGjl(3,2) = (pba_ls*ggaa + pbb_ls*ggba)*cl
           dGjl(3,4) = (pba_ls*ggab + pbb_ls*ggbb)*cl
           !> small-large block
           dGjl(2,1) = (gbbg*paa_sl + gabg*pba_sl)*cj 
           dGjl(2,3) = (gbbg*pab_sl + gabg*pbb_sl)*cj 
           dGjl(4,1) = (gbag*paa_sl + gaag*pba_sl)*cj 
           dGjl(4,3) = (gbag*pab_sl + gaag*pbb_sl)*cj 
        !> end
        !> contraction [ji|lk]*il
           paa_ll = pmax(1,1,i,l)
           pba_ll = pmax(3,1,i,l)
           pab_ll = pmax(1,3,i,l)
           pbb_ll = pmax(3,3,i,l)
           paa_ls = pmax(1,2,i,l)
           pba_ls = pmax(3,2,i,l)
           pab_ls = pmax(1,4,i,l)
           pbb_ls = pmax(3,4,i,l)
           paa_sl = pmax(2,1,i,l)
           pba_sl = pmax(4,1,i,l)
           pab_sl = pmax(2,3,i,l)
           pbb_sl = pmax(4,3,i,l)
           !> large-large block
           dGjk(1,1) = paa_ll*gg
           dGjk(1,3) = pab_ll*gg
           dGjk(3,1) = pba_ll*gg
           dGjk(3,3) = pbb_ll*gg
           !> large-small block
           dGjk(1,2) = (paa_ls*ggbb + pab_ls*ggba)*ck  
           dGjk(1,4) = (paa_ls*ggab + pab_ls*ggaa)*ck
           dGjk(3,2) = (pba_ls*ggbb + pbb_ls*ggba)*ck
           dGjk(3,4) = (pba_ls*ggab + pbb_ls*ggaa)*ck
           !> small-large block
           dGjk(2,1) = (gbbg*paa_sl + gabg*pba_sl)*cj 
           dGjk(2,3) = (gbbg*pab_sl + gabg*pbb_sl)*cj 
           dGjk(4,1) = (gbag*paa_sl + gaag*pba_sl)*cj 
           dGjk(4,3) = (gbag*pab_sl + gaag*pbb_sl)*cj 
        !> end
        !> K-estimate
        Kmax(1,1) = dmax1(dGil(1,1),dGik(1,1),dGjl(1,1),dGjk(1,1),&
                          dGil(3,1),dGik(3,1),dGjl(3,1),dGjk(3,1),&
                          dGil(1,3),dGik(1,3),dGjl(1,3),dGjk(1,3),&
                          dGil(3,3),dGik(3,3),dGjl(3,3),dGjk(3,3) )
        Kmax(2,1) = dmax1(dGil(2,1),dGik(2,1),dGjl(2,1),dGjk(2,1),&
                          dGil(4,1),dGik(4,1),dGjl(4,1),dGjk(4,1),&
                          dGil(2,3),dGik(2,3),dGjl(2,3),dGjk(2,3),&
                          dGil(4,3),dGik(4,3),dGjl(4,3),dGjk(4,3) )
        Kmax(1,2) = dmax1(dGil(1,2),dGik(1,2),dGjl(1,2),dGjk(1,2),&
                          dGil(3,2),dGik(3,2),dGjl(3,2),dGjk(3,2),&
                          dGil(1,4),dGik(1,4),dGjl(1,4),dGjk(1,4),&
                          dGil(3,4),dGik(3,4),dGjl(3,4),dGjk(3,4) )
        Kmax(2,2) = 0.0d0 
        !>
      case('ssss')
        !> integrals processing: [ij|kl]
           !> large-large block 
           gg = gmax(1,i,j)*gmax(1,k,l) ![A|C]
           !> large-small block 
           ggaa = gmax(1,i,j)*(gmax(2,k,l)+gmax(3,k,l)) ![A|c]
           ggab = gmax(1,i,j)*(gmax(4,k,l)+gmax(5,k,l)) ![A|d]
           ggbb = ggaa
           ggba = ggab
           !> small-large block 
           gaag = gmax(1,k,l)*(gmax(2,i,j)+gmax(3,i,j)) ![a|C]
           gabg = gmax(1,k,l)*(gmax(4,i,j)+gmax(5,i,j)) ![b|C]
           gbbg = gaag
           gbag = gabg
           !> small-small block: aa
           gaagaa = (gmax(2,i,j)+gmax(3,i,j))*(gmax(2,k,l)+gmax(3,k,l)) ![a|c]
           gaagba = (gmax(2,i,j)+gmax(3,i,j))*(gmax(4,k,l)+gmax(5,k,l)) ![a|d]
           gabgaa = (gmax(4,i,j)+gmax(5,i,j))*(gmax(2,k,l)+gmax(3,k,l)) ![b|c]
           gabgba = (gmax(4,i,j)+gmax(5,i,j))*(gmax(4,k,l)+gmax(5,k,l)) ![b|d]                                       
           !> 
           gbbgbb = gaagaa
           gaagbb = gaagaa
           gbbgaa = gaagaa
           !>
           gbbgab = gaagba
           gaagab = gaagba
           gbbgba = gaagba
           !>
           gbagbb = gabgaa
           gabgbb = gabgaa
           gbagaa = gabgaa
           !>
           gbagab = gabgba
           gabgab = gabgba
           gbagba = gabgba
        !> end
        !> contraction [ij|kl]*jk
           paa_ll = pmax(1,1,j,k)
           pba_ll = pmax(3,1,j,k)
           pab_ll = pmax(1,3,j,k)
           pbb_ll = pmax(3,3,j,k)
           paa_ls = pmax(1,2,j,k)
           pba_ls = pmax(3,2,j,k)
           pab_ls = pmax(1,4,j,k)
           pbb_ls = pmax(3,4,j,k)
           paa_sl = pmax(2,1,j,k)
           pba_sl = pmax(4,1,j,k)
           pab_sl = pmax(2,3,j,k)
           pbb_sl = pmax(4,3,j,k)
           paa_ss = pmax(2,2,j,k)
           pba_ss = pmax(4,2,j,k)
           pab_ss = pmax(2,4,j,k)
           pbb_ss = pmax(4,4,j,k)
           !> large-large block
           dGil(1,1) = paa_ll*gg
           dGil(1,3) = pab_ll*gg
           dGil(3,1) = pba_ll*gg
           dGil(3,3) = pbb_ll*gg
           !> large-small block
           dGil(1,2) = (paa_ls*ggaa + pab_ls*ggba)*cl  
           dGil(1,4) = (paa_ls*ggab + pab_ls*ggbb)*cl
           dGil(3,2) = (pba_ls*ggaa + pbb_ls*ggba)*cl
           dGil(3,4) = (pba_ls*ggab + pbb_ls*ggbb)*cl
           !> small-large block
           dGil(2,1) = (gaag*paa_sl + gabg*pba_sl)*ci 
           dGil(2,3) = (gaag*pab_sl + gabg*pbb_sl)*ci 
           dGil(4,1) = (gbag*paa_sl + gbbg*pba_sl)*ci 
           dGil(4,3) = (gbag*pab_sl + gbbg*pbb_sl)*ci 
           !> small-small block
           dGil(2,2) = (gaagaa*paa_ss + gaagba*pab_ss + gabgaa*pba_ss + gabgba*pbb_ss)*ci*cl 
           dGil(2,4) = (gaagab*paa_ss + gaagbb*pab_ss + gabgab*pba_ss + gabgbb*pbb_ss)*ci*cl 
           dGil(4,2) = (gbagaa*paa_ss + gbagba*pab_ss + gbbgaa*pba_ss + gbbgba*pbb_ss)*ci*cl 
           dGil(4,4) = (gbagab*paa_ss + gbagbb*pab_ss + gbbgab*pba_ss + gbbgbb*pbb_ss)*ci*cl 
        !> end
        !> contraction [ij|lk]*jl
           paa_ll = pmax(1,1,j,l)
           pba_ll = pmax(3,1,j,l)
           pab_ll = pmax(1,3,j,l)
           pbb_ll = pmax(3,3,j,l)
           paa_ls = pmax(1,2,j,l)
           pba_ls = pmax(3,2,j,l)
           pab_ls = pmax(1,4,j,l)
           pbb_ls = pmax(3,4,j,l)
           paa_sl = pmax(2,1,j,l)
           pba_sl = pmax(4,1,j,l)
           pab_sl = pmax(2,3,j,l)
           pbb_sl = pmax(4,3,j,l)
           paa_ss = pmax(2,2,j,l)
           pba_ss = pmax(4,2,j,l)
           pab_ss = pmax(2,4,j,l)
           pbb_ss = pmax(4,4,j,l)
           !> large-large block
           dGik(1,1) = paa_ll*gg
           dGik(1,3) = pab_ll*gg
           dGik(3,1) = pba_ll*gg
           dGik(3,3) = pbb_ll*gg
           !> large-small block
           dGik(1,2) = (paa_ls*ggbb + pab_ls*ggba)*ck  
           dGik(1,4) = (paa_ls*ggab + pab_ls*ggaa)*ck
           dGik(3,2) = (pba_ls*ggbb + pbb_ls*ggba)*ck
           dGik(3,4) = (pba_ls*ggab + pbb_ls*ggaa)*ck
           !> small-large block
           dGik(2,1) = (gaag*paa_sl + gabg*pba_sl)*ci 
           dGik(2,3) = (gaag*pab_sl + gabg*pbb_sl)*ci 
           dGik(4,1) = (gbag*paa_sl + gbbg*pba_sl)*ci 
           dGik(4,3) = (gbag*pab_sl + gbbg*pbb_sl)*ci 
           !> small-small block
           dGik(2,2) = (gaagbb*paa_ss + gaagba*pab_ss + gabgbb*pba_ss + gabgba*pbb_ss)*ci*ck 
           dGik(2,4) = (gaagab*paa_ss + gaagaa*pab_ss + gabgab*pba_ss + gabgaa*pbb_ss)*ci*ck 
           dGik(4,2) = (gbagbb*paa_ss + gbagba*pab_ss + gbbgbb*pba_ss + gbbgba*pbb_ss)*ci*ck 
           dGik(4,4) = (gbagab*paa_ss + gbagaa*pab_ss + gbbgab*pba_ss + gbbgaa*pbb_ss)*ci*ck 
        !> end
        !> contraction [ji|kl]*ik
           paa_ll = pmax(1,1,i,k)
           pba_ll = pmax(3,1,i,k)
           pab_ll = pmax(1,3,i,k)
           pbb_ll = pmax(3,3,i,k)
           paa_ls = pmax(1,2,i,k)
           pba_ls = pmax(3,2,i,k)
           pab_ls = pmax(1,4,i,k)
           pbb_ls = pmax(3,4,i,k)
           paa_sl = pmax(2,1,i,k)
           pba_sl = pmax(4,1,i,k)
           pab_sl = pmax(2,3,i,k)
           pbb_sl = pmax(4,3,i,k)
           paa_ss = pmax(2,2,i,k)
           pba_ss = pmax(4,2,i,k)
           pab_ss = pmax(2,4,i,k)
           pbb_ss = pmax(4,4,i,k)
           !> large-large block
           dGjl(1,1) = paa_ll*gg
           dGjl(1,3) = pab_ll*gg
           dGjl(3,1) = pba_ll*gg
           dGjl(3,3) = pbb_ll*gg
           !> large-small block
           dGjl(1,2) = (paa_ls*ggaa + pab_ls*ggba)*cl  
           dGjl(1,4) = (paa_ls*ggab + pab_ls*ggbb)*cl
           dGjl(3,2) = (pba_ls*ggaa + pbb_ls*ggba)*cl
           dGjl(3,4) = (pba_ls*ggab + pbb_ls*ggbb)*cl
           !> small-large block
           dGjl(2,1) = (gbbg*paa_sl + gabg*pba_sl)*cj 
           dGjl(2,3) = (gbbg*pab_sl + gabg*pbb_sl)*cj 
           dGjl(4,1) = (gbag*paa_sl + gaag*pba_sl)*cj 
           dGjl(4,3) = (gbag*pab_sl + gaag*pbb_sl)*cj 
           !> small-small block
           dGjl(2,2) = (gbbgaa*paa_ss + gbbgba*pab_ss + gabgaa*pba_ss + gabgba*pbb_ss)*cj*cl 
           dGjl(2,4) = (gbbgab*paa_ss + gbbgbb*pab_ss + gabgab*pba_ss + gabgbb*pbb_ss)*cj*cl 
           dGjl(4,2) = (gbagaa*paa_ss + gbagba*pab_ss + gaagaa*pba_ss + gaagba*pbb_ss)*cj*cl 
           dGjl(4,4) = (gbagab*paa_ss + gbagbb*pab_ss + gaagab*pba_ss + gaagbb*pbb_ss)*cj*cl 
        !> end
        !> contraction [ji|lk]*il
           paa_ll = pmax(1,1,i,l)
           pba_ll = pmax(3,1,i,l)
           pab_ll = pmax(1,3,i,l)
           pbb_ll = pmax(3,3,i,l)
           paa_ls = pmax(1,2,i,l)
           pba_ls = pmax(3,2,i,l)
           pab_ls = pmax(1,4,i,l)
           pbb_ls = pmax(3,4,i,l)
           paa_sl = pmax(2,1,i,l)
           pba_sl = pmax(4,1,i,l)
           pab_sl = pmax(2,3,i,l)
           pbb_sl = pmax(4,3,i,l)
           paa_ss = pmax(2,2,i,l)
           pba_ss = pmax(4,2,i,l)
           pab_ss = pmax(2,4,i,l)
           pbb_ss = pmax(4,4,i,l)
           !> large-large block
           dGjk(1,1) = paa_ll*gg
           dGjk(1,3) = pab_ll*gg
           dGjk(3,1) = pba_ll*gg
           dGjk(3,3) = pbb_ll*gg
           !> large-small block
           dGjk(1,2) = (paa_ls*ggbb + pab_ls*ggba)*ck  
           dGjk(1,4) = (paa_ls*ggab + pab_ls*ggaa)*ck
           dGjk(3,2) = (pba_ls*ggbb + pbb_ls*ggba)*ck
           dGjk(3,4) = (pba_ls*ggab + pbb_ls*ggaa)*ck
           !> small-large block
           dGjk(2,1) = (gbbg*paa_sl + gabg*pba_sl)*cj 
           dGjk(2,3) = (gbbg*pab_sl + gabg*pbb_sl)*cj 
           dGjk(4,1) = (gbag*paa_sl + gaag*pba_sl)*cj 
           dGjk(4,3) = (gbag*pab_sl + gaag*pbb_sl)*cj 
           !> small-small block
           dGjk(2,2) = (gbbgbb*paa_ss + gbbgba*pab_ss + gabgbb*pba_ss + gabgba*pbb_ss)*cj*ck 
           dGjk(2,4) = (gbbgab*paa_ss + gbbgaa*pab_ss + gabgab*pba_ss + gabgaa*pbb_ss)*cj*ck 
           dGjk(4,2) = (gbagbb*paa_ss + gbagba*pab_ss + gaagbb*pba_ss + gaagba*pbb_ss)*cj*ck 
           dGjk(4,4) = (gbagab*paa_ss + gbagaa*pab_ss + gaagab*pba_ss + gaagaa*pbb_ss)*cj*ck
        !> end
        !> K-estimate
        Kmax(1,1) = dmax1(dGil(1,1),dGik(1,1),dGjl(1,1),dGjk(1,1),&
                          dGil(3,1),dGik(3,1),dGjl(3,1),dGjk(3,1),&
                          dGil(1,3),dGik(1,3),dGjl(1,3),dGjk(1,3),&
                          dGil(3,3),dGik(3,3),dGjl(3,3),dGjk(3,3) )
        Kmax(2,1) = dmax1(dGil(2,1),dGik(2,1),dGjl(2,1),dGjk(2,1),&
                          dGil(4,1),dGik(4,1),dGjl(4,1),dGjk(4,1),&
                          dGil(2,3),dGik(2,3),dGjl(2,3),dGjk(2,3),&
                          dGil(4,3),dGik(4,3),dGjl(4,3),dGjk(4,3) )
        Kmax(1,2) = dmax1(dGil(1,2),dGik(1,2),dGjl(1,2),dGjk(1,2),&
                          dGil(3,2),dGik(3,2),dGjl(3,2),dGjk(3,2),&
                          dGil(1,4),dGik(1,4),dGjl(1,4),dGjk(1,4),&
                          dGil(3,4),dGik(3,4),dGjl(3,4),dGjk(3,4) )
        Kmax(2,2) = dmax1(dGil(2,2),dGik(2,2),dGjl(2,2),dGjk(2,2),&
                          dGil(4,2),dGik(4,2),dGjl(4,2),dGjk(4,2),&
                          dGil(2,4),dGik(2,4),dGjl(2,4),dGjk(2,4),&
                          dGil(4,4),dGik(4,4),dGjl(4,4),dGjk(4,4) )
        !>
      case default
        write(6,*) 'Error: unknown eri_class => stop'; stop
    end select
  END SUBROUTINE

! ------------------------------------------------------------------------------------
END MODULE
