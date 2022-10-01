module exacorr_eri

! Wrapper module for getting contracted AO electron repulsion integrals
! Lucas Visscher 2022

  use module_interest_eri

  implicit none

  public compute_eri

  private

  interface compute_eri
      module procedure compute_eri_uncontracted
      module procedure compute_eri_contracted
  end interface


contains

  subroutine compute_eri_uncontracted (eri_type,fijkl,gout,&
                                       lk,ek,xk,yk,zk,ck,&
                                       ll,el,xl,yl,zl,cl,&
                                       li,ei,xi,yi,zi,ci,&
                                       lj,ej,xj,yj,zj,cj )

     character*(*), intent(in) :: eri_type
     integer, intent(in)  :: li,lj,lk,ll
     real(8), intent(in)  :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
     real(8), intent(out) :: gout(:)
     real(8), intent(in)  :: fijkl

     call interest_eri(eri_type,fijkl,gout,&
                        lk,ek,xk,yk,zk,ck,&
                        ll,el,xl,yl,zl,cl,&
                        li,ei,xi,yi,zi,ci,&
                        lj,ej,xj,yj,zj,cj )


  end subroutine

  subroutine compute_eri_contracted (eri_type,fijkl,gout,&
                                     lk,ek,xk,yk,zk,ck,&
                                     ll,el,xl,yl,zl,cl,&
                                     li,ei,xi,yi,zi,ci,&
                                     lj,ej,xj,yj,zj,cj )

     character*(*), intent(in) :: eri_type
     integer, intent(in)  :: li,lj,lk,ll
     real(8), intent(in)  :: xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
     real(8), intent(in)  :: ei(:),ej(:),ek(:),el(:),ci(:),cj(:),ck(:),cl(:)
     real(8), intent(out) :: gout(:)
     real(8), intent(in)  :: fijkl
     integer :: iprim,jprim,kprim,lprim,nprim_i,nprim_j,nprim_k,nprim_l,ngout
     real(8), allocatable :: gout_part(:)

     ! determine the sizes, check them and go to the uncontracted code if possible
     nprim_i = size(ei)
     nprim_j = size(ej)
     nprim_k = size(ek)
     nprim_l = size(el)
     if (nprim_i /= size(ci)) stop "exacorr_eri can not handle general contractions"
     if (nprim_j /= size(cj)) stop "exacorr_eri can not handle general contractions"
     if (nprim_k /= size(ck)) stop "exacorr_eri can not handle general contractions"
     if (nprim_l /= size(cl)) stop "exacorr_eri can not handle general contractions"
     if (nprim_i*nprim_j*nprim_k*nprim_l == 1) then
        call compute_eri_uncontracted (eri_type,fijkl,gout,&
                                       lk,ek(1),xk,yk,zk,ck(1),&
                                       ll,el(1),xl,yl,zl,cl(1),&
                                       li,ei(1),xi,yi,zi,ci(1),&
                                       lj,ej(1),xj,yj,zj,cj(1) )
     else
       ! loop over primitives and gather partial contributions in gout
       ngout = size(gout)
       gout = 0.D0
       allocate(gout_part(ngout))
       do lprim = 1, nprim_l
          do kprim = 1, nprim_k
             do jprim = 1, nprim_j
                do iprim = 1, nprim_i
                   call interest_eri(eri_type,fijkl,gout_part,&
                                     lk,ek(kprim),xk,yk,zk,ck(kprim),&
                                     ll,el(lprim),xl,yl,zl,cl(lprim),&
                                     li,ei(iprim),xi,yi,zi,ci(iprim),&
                                     lj,ej(jprim),xj,yj,zj,cj(jprim) )
                    gout = gout + gout_part
                end do
             end do
          end do
       end do
       deallocate(gout_part)

     end if

  end subroutine

end module
