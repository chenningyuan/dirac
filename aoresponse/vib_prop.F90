module vib_prop
!Calculation and outputting of vibrational (optical) properties.
!This file is currently common to Dalton and Dirac

   use matrix_defop
   use prop_integrals
   use rsp_equations
#ifdef PRG_DIRAC
   use dirac_interface
#endif
   use aor_restart

   implicit none

   public roa_pol_gpri_aten_gra,cars_shyp_polgra
   public vib_ana_polari,print_shypol
   public dg_energy

   private

   !wavenumber, velocity of light, and nanometer (in au), pi
   real(8),parameter:: cm1=1/219474.631371d0,cvl=137.03599907d0, &
                       nm=10/0.52917706d0,pi=3.1415926535897931d0

   !field component lables for printing
   character*2:: fc(3)=(/'Fx','Fy','Fz'/)

contains


   subroutine print_tensor(dims, tensor, title, freqs, colwidth, unit)
   !ajt This is used for printing response function tensors
   !TODO: Add comp. lables argument. Make space between blocks of rows
      integer,                intent(in) :: dims(:)
      complex(8),             intent(in) :: tensor(*) !(product(dims))
      character(*), optional, intent(in) :: title
      integer,      optional, intent(in) :: unit, colwidth
      complex(8),   optional, intent(in) :: freqs(:)
      integer      :: uni, siz, dec, cwid, i
      character*9 :: fmt
      uni = 6
      if (present(unit)) uni = unit
      if (present(title)) then
         if (present(freqs)) then
            fmt = '(a,nf6.3)'
            write (fmt(4:4),'(i1)') size(freqs)
            write (uni,fmt) title//' w:',real(freqs)
         else
            write (uni,'(a)') title
         end if
      end if
      siz = product(dims)
      dec = max(2, 1+ceiling(log10(maxval(abs(tensor(1:siz))))))
      !dec = max(dec,ceiling(log10( maxval(tensor(1:siz)))))
      cwid = 15
      if (present(colwidth)) cwid = max(colwidth,dec+2)
      dec = cwid - dec - 1
      fmt = '(fww.dd)'
      write (fmt(3:7),'(i2,a1,i2)') cwid, '.', dec
      if (any( imag(tensor(1:siz))/=0 )) write (uni,'(a)') '    (real part)'
      call subr(dims, real(tensor(1:siz)), uni, cwid, fmt)
      if (any( imag(tensor(1:siz))/=0 )) then
         write (uni,'(a)') '    (imaginary part)'
         call subr(dims, imag(tensor(1:siz)), uni, cwid, fmt)
      end if
      write (uni,'()') !final blank line
   contains
      subroutine subr(dims, tnz, unit, cwid, fmt)
         integer,     intent(in) :: dims(:), unit, cwid
         real(8),     intent(in) :: tnz(*)
         character*8, intent(in) :: fmt
         integer :: i, j, l
         character*(dims(1)*(cwid+1)) :: line
         do i = 1, product(dims), dims(1)
            l = 1
            do j = 1, dims(1)
               line(l:l) = ' '
               write (line(l+1:l+cwid), fmt) tnz(i+j-1)
               l = l + 1 + cwid
            end do
            write (unit,'(a)') line
         end do
      end subroutine
   end subroutine


   subroutine load_vib_modes(nc,nq,w,Q)
   !load Hessian form DALTON.HES, mass-weight, diagonalize
   !and sqrt eigenvalues for frequencies (w<0 means imaginary)
      integer,intent(in)  :: nc !number of coordinates
      integer,intent(out) :: nq !number of vib. modes
      real(8),intent(out) :: w(nc),Q(nc,nc) !frequencies, normal modes Q(:,k)
      real(8)      :: M(nc/3),Mtot,Z(nc/3),G(nc),Rc(nc),H(nc,nc),orig(3)=0
! Bin Gao: change unit to -1 for Dalton main branch
      integer      :: IS(nc/3),i,internal=0,negative=0,unit=-1
      logical      :: ex
      character*22 :: fmt
      call NUCLEI_ifc(nc/3,Z,IS,G)
      call GPINQ('DALTON.HES','EXIST',ex)
      if (.not.ex) call quit('load_vib_modes: Hessian file DALTON.HES not found')
      call GPOPEN(unit,'DALTON.HES','OLD',' ','FORMATTED',0,.false.)
      read (unit,'(i12)') i
      if (i/=nc) call quit('load_vib_modes: DALTON.HES inconsistent with MOLECULE.INP')
      write (fmt,'(a,i4,a,i4,a)') '(',nc,'(/',nc,'(g24.17/)))'
      read (unit,fmt) H
      write (fmt(1:15),'(a,i3,a)') '(',nc,'(/g24.17))'
      read (unit,fmt,iostat=i) Rc !if no geometry at the end, reads EOF
      if (i==0.and..not.all(G==Rc)) & !if Rc read, compare
         call quit('load_vib_modes: Different geometry in DALTON.HES and MOLECULE.INP')
      call GPCLOSE(unit,'KEEP')
      call VIBHES(1,nc,G,nq,(/(0d0,i=1,nc)/),H,(/(0d0,i=1,6*nc)/) &
                 ,nint(Z),(/(0d0,i=1,2*nc*nc)/)); nq=nc-nq
      call VIBMAS(M,Mtot,IS,nint(Z),nc/3,G,orig,1)
      call VIBNOR(H,M,(/(0d0,i=1,nc)/),(/(0d0,i=1,nc*(nc+1)/2)/) &
                 ,w,Q,(/(0d0,i=1,nc*nc)/),(/(0d0,i=1,30*nc*nc+1)/) &
                 ,30*nc*nc+1,nc,internal,negative,.false.,1)
      call NUCLEI_ifc(nc/3,Z,IS,G)
      w = sign(sqrt(abs(w)),w)
   end subroutine


   subroutine print_vib_table(nq,nc,ne,lab,wc,ec,tab,o)
      integer           :: nq,nc,ne,wc(nc),ec(ne+1),o
      character(*)      :: lab(nc) !labels
      real(8)           :: tab(nc,nq) !values to be printed
      integer           :: i,j,k,l,w(nc),d(nc) !col width, decimals
      character*(32*nc) :: p !for pre-printing into
      character*8       :: fmt
      !count digits before decimal point
      do j=1,nc
         w(j) = 1+floor(log10(maxval(abs(tab(j,:)))))
         if (w(j)<0.and.w(j)<-wc(j)+3) then !zero column
            w(j) = 0
         elseif (all(tab(j,:)==nint(tab(j,:)))) then !integer column
            d(j) = -1
            if (any(tab(j,:)<0)) w(j) = w(j)+1
            w(j) = max(w(j)+1,wc(j))
         else !real column
            w(j) = 3+max(1,w(j))
            d(j) = max(0,wc(j)-w(j))
            w(j) = max(w(j),wc(j))
         end if
      end do
      !header
      l = 0
      do j=1,nc
         if (w(j)==0) cycle !skip column
         write (p(l+1:),'(a)') lab(j)(1:w(j))
         l = l+w(j)
      end do
      write (o,'(a)') trim(p(:l))
      !dashes ------
      write (o,'(1x,a)') repeat('-',l-1)
      !the table, line by line
      do i=1,nq
         l = 0
         do j=1,nc
            if (w(j)==0) cycle !skip column
            if (d(j)==-1) then
               write (fmt(1:5),'(a,i2,a)') '(i',w(j),')'
               write (p(l+1:),fmt(1:5)) nint(tab(j,i))
            else
               write (fmt,'(a,i2,a,i2,a)') '(f',w(j),'.',d(j),')'
               write (p(l+1:),fmt) tab(j,i)
            end if
            l = l+w(j)
         end do
         write (o,'(a)') p(:l)
      end do
      !dashes ------
      write (o,'(1x,a)') repeat('-',l-1)
      !column sums for 'energy' columns
      l = max(4,sum(w(:ec(1)-1)))
      write (p(1:),'(a)') repeat(' ',l-3)//'sum'
      do j=ec(1),ec(ne+1)-1
         if (w(j)==0.or.d(j)==-1) cycle
         write (fmt,'(a,i2,a,i2,a)') '(f',w(j),'.',d(j),')'
         write (p(l+1:),fmt) sum(tab(j,:))
         l = l+w(j)
      end do
      write (o,'(a)') p(:l)
      !total of 'energy' columns
      l = max(6,sum(w(:ec(1)-1)))
      write (p(1:),'(a)') repeat(' ',l-5)//'total'
      k = 1
      do j=ec(1),ec(ne)
         if (w(j)==0.or.d(j)==-1) cycle
         if (j==ec(k)) then
            write (fmt,'(a,i2,a,i2,a)') '(f',w(j),'.',d(j),')'
            write (p(l+1:),fmt) sum(tab(ec(k):ec(k+1)-1,:))
            k = k+1
         else
            write (p(l+1:),'(a)') repeat(' ',w(j))
         end if
         l = l+w(j)
      end do
      write (o,'(a/)') trim(p(:l))
   end subroutine


   subroutine print_polari(dax,w,pol,o)
   !Print the polarizability tensor pol, along with it's isotropic
   !and dipole-dipole components. dax is the NORMALIZED dipol axis.
      real(8) :: dax(3),w,pol(3,3)
      integer :: o,i,j
      write (o,'(//26x,a/26x,a/)') 'Polarizability (au)' &
                                  ,'-------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(2(8x,a,i1,a,f6.3,a,i7,a,i4,a/))') &
               'w',0,' = ', w,' au =', nint(w/cm1),' cm-1 ~ ', &
                                  nint((2*pi*cvl/w)/nm),' nm', &
               'w',1,' = ',-w,' au =',-nint(w/cm1),' cm-1 ~ ', &
                                 -nint((2*pi*cvl/w)/nm),' nm'
      write (o,'(a,13x,a2,18x,a2,18x,a2/)') ' w0 w1',fc
      write (o,'(3(1x,a,3f20.13/))') (fc(i),pol(i,:),i=1,3)
      write (o,'(a,g16.9)')  '@alpha    isotropic: ', &
               (pol(1,1)+pol(2,2)+pol(3,3))/3
      write (o,'(a,g16.9/)') '@alpha     dipole^2: ', &
               sum((/((pol(i,j)*dax(i)*dax(j),i=1,3),j=1,3)/))
   end subroutine


   subroutine print_shypol(w,dip,pol,hyp,shyp,o)
      real(8) :: w(4),dip(3),pol(3,3,4),hyp(3,3,3,6),shyp(3,3,3,3)
      integer :: o,i,j,k,l,m,ij
      real(8) :: a,b,dax(3)
      write (o,'(//28x,a/28x,a/)') 'Dipole moment (au)', &
                                   '------------------'
      write (o,'(1x,3(18x,a)//(3x,3f20.14/))') fc,dip
      if (sqrt(sum(dip*dip))<0.00001) then
         write (o,'(4x,a/)') 'The molecule has zero dipole moment' &
                           //' - using z as dipole axis'
         dax = (/0d0,0d0,1d0/)
      else
         dax = dip/sqrt(sum(dip*dip))
      end if
      do i=1,4; a=w(i)
         if (any(w(1:i-1)==a.or.w(1:i-1)==-a)) cycle
         if (all(pol(:,:,i)==0)) cycle
         call print_polari(dax,w(i),pol(:,:,i),o)
      end do
      do j=2,4; do i=1,j-1; a=w(i); b=w(j); ij=i+(j-1)*(j-2)/2
         !only print one beta for each unique freqency triple
         if (i==1.and.j==3.and.any(w(2)==(/-a-b,b/))) cycle
         if (i==2.and.j==3.and.any(w(1)==(/-a-b,a,b/))) cycle
         if (i==1.and.j==4.and.any(w(2)==(/-a-2*b,-a-b,-a,-b,0d0,b/))) cycle
         if (i==2.and.j==4.and.any(w(1)==(/-a-2*b,-a-b,-a,-b,0d0,b,a/))) cycle
         if (i==3.and.j==4.and.any(w(2)==(/-2*a-b,-a-2*b,-a-b,-a,-b,0d0,b,a/))) cycle
         if (all(hyp(:,:,:,ij)==0)) cycle
         write (o,'(//24x,a/24x,a/)') 'Hyperpolarizability (au)' &
                                     ,'------------------------'
         write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
         write (o,'(3(8x,a,i1,a,f6.3,a,i7,a,i4,a/))') &
                  'w',0,' = ',   a,' au =', nint(     a/cm1),' cm-1 ~ ', &
                                      nint((2*pi*cvl/a     )/nm),' nm',  &
                  'w',1,' = ',   b,' au =', nint(     b/cm1),' cm-1 ~ ', &
                                      nint((2*pi*cvl/b     )/nm),' nm',  &
                  'w',2,' = ',-a-b,' au =', nint((-a-b)/cm1),' cm-1 ~ ', &
                                      nint((2*pi*cvl/(-a-b))/nm),' nm'
         write (o,'(a,13x,a2,18x,a2,18x,a2)') ' w0w1 w2',fc
         write (o,'(3(/3(1x,2a2,3f20.12/)))') &
                 ((fc(k),fc(l),hyp(k,l,:,ij),l=1,3),k=1,3)
         write (o,'(a,g16.9)') '@beta     dip-iso: ',1/5d0*sum((/(( &
                 (hyp(k,k,l,ij)+hyp(k,l,k,ij)+hyp(l,k,k,ij))*dax(l), &
                  k=1,3),l=1,3)/))
         write (o,'(a,g16.9/)') '@beta    dipole^3: ', &
                  sum((/(((hyp(k,l,m,ij)*dax(k)*dax(l)*dax(m),k=1,3),l=1,3),m=1,3)/))
      end do; end do
      if (all(shyp==0)) return
      write (o,'(//22x,a/22x,a/)') 'Second hyperpolarizability (au)', &
                                   '-------------------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(4(8x,a,i1,a,f6.3,a,i7,a,i4,a/))')                   &
              ('w',i-1,' = ',w(i),' au =',nint(w(i)/cm1),' cm-1 ~ ', &
                             nint((2*pi*cvl/w(i))/nm),' nm',i=1,4)
      write (o,'(a,13x,a2,18x,a2,18x,a2)') ' w0w1w2 w3',fc
      write (o,'(3(3(/3(1x,3a2,3f20.10/))))') &
              (((fc(i),fc(j),fc(k),shyp(i,j,k,:),k=1,3),j=1,3),i=1,3)
      write (o,'(a,g16.9)') '@gamma    isotropic: ',1/15d0*sum((/(( &
              shyp(i,i,j,j)+shyp(i,j,i,j)+shyp(i,j,j,i),i=1,3),j=1,3)/))
      write (o,'(a,g16.9)') '@gamma    dip^2-iso: ',1/10d0*sum((/(((     &
              (shyp(i,i,k,l)+shyp(i,k,i,l)+shyp(k,i,i,l)                 &
              +shyp(i,k,l,i)+shyp(k,i,l,i)+shyp(i,i,k,l))*dax(k)*dax(l), &
               k=1,3),l=1,3),i=1,3)/))
      write (o,'(a,g16.9)') '@gamma        dip^4: ',                  &
               1*sum((/((((shyp(i,j,k,l)*dax(i)*dax(j)*dax(k)*dax(l), &
                                      i=1,3),j=1,3),k=1,3),l=1,3)/))
      write (o,'(4(a,1x,g13.7/)/)')                                      &
               '@gamma    Isotropic (AABB+ABAB+ABBA)/15:',               &
               sum((/((shyp(i,i,j,j)+shyp(i,j,i,j)+shyp(i,j,j,i),        &
                                            i=1,3),j=1,3)/))/15,         &
               '@gamma                   XXXX-component:',shyp(1,1,1,1), &
               '@gamma                   YYYY-component:',shyp(2,2,2,2), &
               '@gamma                   ZZZZ-component:',shyp(3,3,3,3)
   end subroutine


   subroutine vib_ana_polari(w,dip,nc,dipg,polg,hypg,o)
      real(8) :: dipg(3,nc,4),polg(3,3,nc,6),hypg(3,3,3,nc,4),w(4),dip(3)
      integer :: nc,o,nq,i,j,k,l,m,n,ij,ik,jk
      real(8) :: dax(3),dipq(3,nc*4),polq(3,3,nc*6),hypq(3,3,3,nc*4),wq(nc),Q(nc,nc)
      if (sqrt(sum(dip*dip))< 0.00001) dax = (/0d0,0d0,1d0/)
      if (sqrt(sum(dip*dip))>=0.00001) dax = dip/sqrt(sum(dip*dip))
      call load_vib_modes(nc,nq,wq,Q)
      write (o,'(//20x,a/20x,a/)') 'Mass-weighted normal coordinates (au)', &
                                   '-------------------------------------'
      call OUTPUT(Q(:,1:nq),1,nc,1,nq,nc,nq,1,o)
      do i=1,4; call DGEMM('N','N',3,nq,nc,1d0,dipg(:,:,i), &
              3,Q(:,1:nq),nc,0d0,dipq(:,1+nq*(i-1):nq*i),3); end do
      do i=1,6; call DGEMM('N','N',3*3,nq,nc,1d0,polg(:,:,:,i), &
              3*3,Q(:,1:nq),nc,0d0,polq(:,:,1+nq*(i-1):nq*i),3*3); end do
      do i=1,4; call DGEMM('N','N',3*3*3,nq,nc,1d0,hypg(:,:,:,:,i), &
              3*3*3,Q(:,1:nq),nc,0d0,hypq(:,:,:,1+nq*(i-1):nq*i),3*3*3); end do
      do i=1,4
         if (any(w(1:i-1)==w(i).or.w(1:i-1)==-w(i))) cycle
         if (all(dipq(:,1+nq*(i-1):nq*i)==0)) cycle
         call print_vib_alpha(dax,w(i),nq,wq(1:nq),dipq(:,1+nq*(i-1):nq*i),o)
      end do
      do n=4,1,-1; if (w(n)/=0) cycle
         i=merge(1,2,n>1); j=merge(2,3,n>2); k=merge(3,4,n>3)
         if (all(dipq(:,1+nq*(i-1):nq*i)==0).and. &
             all(dipq(:,1+nq*(j-1):nq*j)==0).and. &
             all(dipq(:,1+nq*(k-1):nq*k)==0)) cycle
         ij=i+(j-1)*(j-2)/2; ik=i+(k-1)*(k-2)/2; jk=j+(k-1)*(k-2)/2
         if (all(polq(:,:,1+nq*(ij-1):nq*ij)==0).and. &
             all(polq(:,:,1+nq*(ik-1):nq*ik)==0).and. &
             all(polq(:,:,1+nq*(jk-1):nq*jk)==0)) cycle
         call print_vib_beta(dax,(/w(i),w(j),w(k)/),nq,wq(1:nq), &
                 -(/dipq(:,1+nq*(i-1):nq*i),dipq(:,1+nq*(j-1):nq*j),dipq(:,1+nq*(k-1):nq*k)/), &
                 -(/polq(:,:,1+nq*(ij-1):nq*ij),polq(:,:,1+nq*(ik-1):nq*ik), &
                   polq(:,:,1+nq*(jk-1):nq*jk)/),o)
         exit !only one
      end do
      call print_vib_gamma(dax,w,nq,wq(1:nq),dipq(:,1:nq*4),    &
                           polq(:,:,1:nq*6),hypq(:,:,:,1:nq*4),o)
   end subroutine


   subroutine print_vib_alpha(dax,w,nq,wq,dipg,o)
      real(8) :: dax(3),w,wq(nq),dipg(3,nq)
      integer :: o,nq,i,j,k,m,q
      real(8) :: pol(3,3),R(3,3),t(5,nq)
      write (o,'(//21x,a/21x,a/)') 'Vibrational polarizability (au)', &
                                   '-------------------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(2(8x,a,i1,a,f6.3,a,i7,a,i4,a/))')             &
               'w',0,' = ', w,' au =', nint(w/cm1),' cm-1 ~ ', &
                                  nint((2*pi*cvl/w)/nm),' nm', &
               'w',1,' = ',-w,' au =',-nint(w/cm1),' cm-1 ~ ', &
                                 -nint((2*pi*cvl/w)/nm),' nm'
      write (o,'(4x,a,3f8.4/)') 'Normalized dipole axis: ',dax
      pol = 0
      do q=1,nq
         R(:,:) = reshape((/((dipg(j,q)*dipg(i,q),i=1,3),j=1,3)/),(/3,3/))
         pol = pol + R(:,:)/(wq(q)**2-w**2)
         t(1:3,q) = (/1d0*q,wq(q)/cm1,wq(q)/)
         t(4:5,q) = (/1d0/(wq(q)**2-w**2),1d0/)*1/3d0*sum((/(R(i,i),i=1,3)/))
      end do
      write (o,'(a,13x,a2,18x,a2,18x,a2/)') ' w0 w1',fc
      write (o,'(3(1x,a2,3f20.13/))') (fc(i),pol(i,:),i=1,3)
      write (o,'(a,g16.9)') '@alpha_vib    isotropic:',sum(t(4,:))
      write (o,'(a,g16.9/)') '@alpha_vib     dipole^2:', &
               sum((/((pol(i,j)*dax(i)*dax(j),i=1,3),j=1,3)/))
      if (abs(w)>0.03) t(5,:)=0
      write (o,'(4x,a/)') 'Table of mode contributions and dominant residues (au):'
      write (o,'(8x,a)') 'mm   :  dmu/dQ x dmu/dQ contribution (isotropic)'
      if (abs(w)<=0.03) write (o,'(8x,a)') 'm1   :  <<Q;F>>_w1 induced gradient with freq -w1 =w0'
      call print_vib_table(nq,5,1,(/'         ','freq/cm-1', &
               '  freq/au','    mm   ','    m0m1 '/),        &
               (/3,9,10,11,11/),(/4,5/),t,o)
   end subroutine


   subroutine print_vib_beta(dax,w,nq,wq,dipg,polg,o)
      real(8) :: dax(3),w(3),wq(nq),dipg(3,nq,3),polg(3,3,nq,3)
      real(8) :: hyp(3,3,3),R(3,3,3,4),t(15,nq)
      integer :: o,nq,i,j,k,n,q
      write (o,'(//20x,a/20x,a/)') 'Vibrational hyperpolarizability (au)', &
                                   '------------------------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(3(8x,a,i1,a,f6.3,a,i7,a,i4,a/))')                   &
              ('w',i-1,' = ',w(i),' au =',nint(w(i)/cm1),' cm-1 ~ ', &
                             nint((2*pi*cvl/w(i))/nm),' nm',i=1,3)
      write (o,'(4x,a,3f8.4/)') 'Normalized dipole axis: ',dax
      hyp = 0
      do q=1,nq
         R(:,:,:,2:4) = reshape((/ &
                 (((dipg(k,q,3)*polg(i,j,q,1),i=1,3),j=1,3),k=1,3), &
                 (((dipg(j,q,2)*polg(i,k,q,2),i=1,3),j=1,3),k=1,3), &
                 (((dipg(i,q,1)*polg(j,k,q,3),i=1,3),j=1,3),k=1,3)/),(/3,3,3,3/))
         R(:,:,:,1) = R(:,:,:,2)/(wq(q)**2-w(3)**2) &
                   + R(:,:,:,3)/(wq(q)**2-w(2)**2) &
                   + R(:,:,:,4)/(wq(q)**2-w(1)**2)
         hyp = hyp + R(:,:,:,1)
         t(1: 4,q) = (/(sum((/(((R(i,i,j,k)+R(i,j,i,k)+R(j,i,i,k))*dax(j),i=1,3),j=1,3)/))/5,k=1,4)/)
         t(5: 8,q) = (/(sum((/(((R(i,i,j,k)-R(i,j,i,k))*dax(j),i=1,3),j=1,3)/))/5,k=1,4)/)
         t(9:12,q) = (/(sum((/(((R(i,i,j,k)-R(j,i,i,k))*dax(j),i=1,3),j=1,3)/))/5,k=1,4)/)
         t(:,q) = (/1d0*q,wq(q)/cm1,wq(q),t(1,q),t(5,q),t(9,q),t(2:4,q),t(6:8,q),t(10:12,q)/)
      end do
      write (o,'(a,13x,a2,18x,a2,18x,a2)') ' w0w1 w2',Fc
      write (o,'(3(/3(1x,2a2,3f20.12/)))') ((Fc(i),Fc(j),hyp(i,j,:),j=1,3),i=1,3)
      write (o,'(a,g16.9)') '@beta_vib    dip-isotropic: ',sum(t(4,:))
      write (o,'(a,g16.9/)') '@beta_vib         dipole^3: ', &
               sum((/(((hyp(i,j,k)*dax(i)*dax(j)*dax(k),i=1,3),j=1,3),k=1,3)/))
      if (w(3)==0.or.maxval(abs(t(5,:)))<1d-8*maxval(abs(t(4,:)))) then
          t(5,:)=0; t(10:12,:)=0; end if
      if (w(2)==0.or.maxval(abs(t(6,:)))<1d-8*maxval(abs(t(4,:)))) then
          t(6,:)=0; t(13:15,:)=0; end if
      if (w(3)==w(2).or.w(3)==w(1).or.abs(w(3))>0.03) t(7,:)=0
      if (w(2)==w(1).or.abs(w(2))>0.03) t(8,:)=0
      if (abs(w(1))>0.03) t(9,:)=0
      do i=7,9; if (w(3)==0.and.all(t(i,:)==0)) t(i+3,:)=0; end do
      do i=7,9; if (w(2)==0.and.all(t(i,:)==0)) t(i+6,:)=0; end do
      write (o,'(4x,a/)') 'Table of mode contributions and dominant residues (au):'
      write (o,'(8x,a)') 'ma   :  dmu/dQ x dalpha/dQ contribution (dipole-isotropic)'
      if (any(t(5,:)/=0)) write (o,'(8x,a)') "ma'  :  as ma, but anisotropic: sum(iiz-izi)/5"
      if (any(t(6,:)/=0)) write (o,'(8x,a)') 'ma"  :  as ma, but anisotropic: sum(iiz-zii)/5'
      write (o,'(8x,a)') 'm2   :  <<Q;F>>_w2 induced gradient with freq -w2 =w0+w1'
      write (o,'(8x,a/)') 'a01  :  <<Q;F,F>>_w0,w1 ind. grad. with freq -w0-w1 =w2'
      call print_vib_table(nq,15,3,                                             &
               (/'         ','freq/cm-1','  freq/au','    ma   ',"    ma'  ",   &
                 '    ma"  ','   m2a01 ','   m1a02 ','   m0a12 ',"   m2a01'",   &
                 "   m1a02'","   m0a12'",'   m2a01"','   m1a02"','   m0a12"'/), &
               (/3,9,10,(11,i=1,12)/),(/4,5,6,7/),t,o)
   end subroutine


   subroutine print_vib_gamma(dax,w,nq,wq,dipq,polq,hypq,o)
      real(8) :: dax(3),w(4),wq(nq),dipq(3,nq,4),polq(3,3,nq,6),hypq(3,3,3,nq,4)
      real(8) :: shyp(3,3,3,3),R(3,3,3,3,9),t(30,nq)
      integer :: o,nq,i,j,k,l,n,q
      write (o,'(//19x,a/19x,a/)') 'Vibrational second hyperpolarizability (au)', &
                                   '-------------------------------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(4(8x,a,i1,a,f6.3,a,i7,a,i4,a/))')                   &
              ('w',i-1,' = ',w(i),' au =',nint(w(i)/cm1),' cm-1 ~ ', &
                             nint((2*pi*cvl/w(i))/nm),' nm',i=1,4)
      write (o,'(4x,a,3f8.4/)') 'Normalized dipole axis: ',dax
      shyp = 0
      do q=1,nq
         R(:,:,:,:,3:9) = reshape((/ &
                 ((((polq(j,k,q,3)*polq(i,l,q,4),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((polq(i,k,q,2)*polq(j,l,q,5),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((polq(i,j,q,1)*polq(k,l,q,6),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((dipq(l,q,4)*hypq(i,j,k,q,1),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((dipq(k,q,3)*hypq(i,j,l,q,2),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((dipq(j,q,2)*hypq(i,k,l,q,3),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((dipq(i,q,1)*hypq(j,k,l,q,4),i=1,3),j=1,3),k=1,3),l=1,3)/), &
                 (/3,3,3,3,7/))
         R(:,:,:,:,1) = R(:,:,:,:,3)/(wq(q)**2-(w(2)+w(3))**2) &
                      + R(:,:,:,:,4)/(wq(q)**2-(w(1)+w(3))**2) &
                      + R(:,:,:,:,5)/(wq(q)**2-(w(1)+w(2))**2)
         R(:,:,:,:,2) = R(:,:,:,:,6)/(wq(q)**2-w(4)**2) &
                      + R(:,:,:,:,7)/(wq(q)**2-w(3)**2) &
                      + R(:,:,:,:,8)/(wq(q)**2-w(2)**2) &
                      + R(:,:,:,:,9)/(wq(q)**2-w(1)**2)
         shyp = shyp + R(:,:,:,:,1) + R(:,:,:,:,2)
         t( 1: 9,q) = (/(sum((/((R(i,i,j,j,k)+R(i,j,i,j,k)+R(i,j,j,i,k),i=1,3),j=1,3)/)),k=1,9)/)/15
         t(10:18,q) = (/(sum((/((R(i,i,j,j,k)-R(i,j,j,i,k),i=1,3),j=1,3)/)),k=1,9)/)/10
         t(19:27,q) = (/(sum((/((R(i,i,j,j,k)+R(i,j,j,i,k)-2*R(i,j,i,j,k),i=1,3),j=1,3)/)),k=1,9)/)/10
         t(:,q) = (/1d0*q,wq(q)/cm1,wq(q),t(1:2,q),t(10:11,q), &
                    t(19:20,q),t(3:9,q),t(12:18,q),t(21:27,q)/)
      end do
      write (o,'(a,13x,a2,18x,a2,18x,a2)') ' w0w1w2 w3',Fc
      write (o,'(3(3(/3(1x,3a2,3f20.10/))))') (((Fc(k),Fc(j),Fc(i),shyp(k,j,i,:),i=1,3),j=1,3),k=1,3)
      write (o,'(a,g16.9)') '@gamma_vib    isotropic: ',sum(t(4:5,:))
      write (o,'(a,g16.9)') '@gamma_vib    dip^2-iso: ',                           &
               sum((/((((shyp(i,i,k,l)+shyp(i,k,i,l)+shyp(k,i,i,l)                 &
                        +shyp(i,k,l,i)+shyp(k,i,l,i)+shyp(i,i,k,l))*dax(k)*dax(l), &
                                                     k=1,3),l=1,3),i=1,3)/)) / 10
      write (o,'(a,g16.9)') '@gamma_vib        dip^4: ',            &
               sum((/((((shyp(i,j,k,l)*dax(i)*dax(j)*dax(k)*dax(l), &
                                   i=1,3),j=1,3),k=1,3),l=1,3)/))
      write (o,'(4(a,1x,g13.7/))')                                           &
               '@gamma_vib    Isotropic (AABB+ABAB+ABBA)/15:',sum(t(4:5,:)), &
               '@gamma_vib                   XXXX-component:',shyp(1,1,1,1), &
               '@gamma_vib                   YYYY-component:',shyp(2,2,2,2), &
               '@gamma_vib                   ZZZZ-component:',shyp(3,3,3,3)
      !zero (skip) repeated and high-frequency columns
      if (maxval(abs(t(6:7,:)))<1d-8*maxval(abs(t(4:5,:)))) then
          t(6:7,:)=0; t(17:23,:)=0; end if
      if (maxval(abs(t(8:9,:)))<1d-8*maxval(abs(t(4:5,:)))) then
          t(8:9,:)=0; t(24:30,:)=0; end if
      if (w(4)==w(3).or.w(4)==w(2).or.abs(w(2)+w(3))>0.03) t(10,:)=0
      if (w(3)==w(2).or.abs(w(1)+w(3))>0.03) t(11,:)=0
      if (abs(w(2)+w(3))>0.03) t(12,:)=0
      if (w(4)==w(3).or.w(4)==w(2).or.w(4)==w(1).or.abs(w(4))>0.03) t(13,:)=0
      if (w(3)==w(2).or.w(3)==w(1).or.abs(w(3))>0.03) t(14,:)=0
      if (w(2)==w(1).or.abs(w(2))>0.03) t(15,:)=0
      if (abs(w(1))>0.03) t(16,:)=0
      if (w(2)==w(3).or.w(2)==w(4)) then; do i=17,23
              if (all(t(i- 7,:)==0)) t(i,:)=0; end do; end if
      if (w(2)==w(4).or.w(1)==w(3)) then; do i=24,30
              if (all(t(i-14,:)==0)) t(i,:)=0; end do; end if
      !print the table
      write (o,'(4x,a/)') 'Table of mode contributions and dominant residues (au):'
      write (o,'(8x,a)') 'aa   :  dalpha/dQ x dalpha/dQ contribution (isotropic)'
      if (any(t(6:7,:)/=0)) write (o,'(8x,a)') "aa'  :  as aa, but anisotropic: sum(iijj-ijji)/10"
      if (any(t(8:9,:)/=0)) write (o,'(8x,a)') 'aa"  :  as aa, but anisotropic: sum(iijj-2ijij+ijji)/10'
      write (o,'(8x,a)') 'mb   :  dmu/dQ x dbeta/dQ  contribution (isotropic)'
      write (o,'(8x,a)') 'm3   :  <<Q;F>>_w3 induced gradient with freq -w3 =w0+w1+w2'
      write (o,'(8x,a)') 'a01  :  <<Q;F,F>>_w0,w1 ind. grad. with freq -w0-w1 =w2+w3'
      write (o,'(8x,a/)') 'b123 :  <<Q;F,F,F>>_w1,w2,w3 ind. grad. with freq -w1-w2-w3 =w0'
      call print_vib_table(nq,30,3,                                                         &
               (/'         ','freq/cm-1','  freq/au','    aa   ','    mb   ',"    aa'  ",   &
                 "    mb'  ",'    aa"  ','    mb"  ','  a12a03 ','  a02a13 ','  a01a23 ',   &
                 '  m3b012 ','  m2b013 ','  m1b023 ','  m0b123 ',"  a12a03'","  a02a13'",   &
                 "  a01a23'","  m3b012'","  m2b013'","  m1b023'","  m0b123'",'  a12a03"',   &
                 '  a02a13"','  a01a23"','  m3b012"','  m2b013"','  m1b023"','  m0b123"'/), &
               (/3,9,10,(11,i=1,27)/),(/4,6,8,10/),t,o)
   end subroutine


   subroutine roa_pol_gpri_aten_gra(S,D,F,nc,w,Eff,Efb,Efq,Effg,Efbg,Efqg)
   !For frequency w, Eff returns polarizability, Efb G-prime, Efq A-tensor,
   !Effg polarizability gradient, Efbg G-prime gradient, Efqg A-tensor gradient

     type(matrix), intent(in)  :: S, D, F
     integer,      intent(in)  :: nc
     complex(8),   intent(in)  :: w
     complex(8),   intent(out) :: Eff (3, 3)
     complex(8),   intent(out) :: Efb (3, 3)
     complex(8),   intent(out) :: Efq (3, 6)
     complex(8),   intent(out) :: Effg(3, 3, nc)
     complex(8),   intent(out) :: Efbg(3, 3, nc)
     complex(8),   intent(out) :: Efqg(3, 6, nc)

     character(4)              :: perturbation
     character(8)              :: restart_label
     character(10)             :: perturbation_print

!    tensor printout needs exact ranks
!    which are not known here
     complex(8),   allocatable :: Ef_j (:, :)
     complex(8),   allocatable :: Ef_jg(:, :, :)
                   
     integer                   :: r_i
     integer                   :: r_j
     integer,      parameter   :: r_f = 3
     integer,      parameter   :: r_e = 3
     integer,      parameter   :: r_b = 3
     integer,      parameter   :: r_q = 6

     integer                   :: i, j, k
                               
     type(matrix)              :: Sg(1), DFD3(3, 3), DFD6(3, 6), Df(3), Ff(3)
                               
     type(matrix)              :: D_j3(3), D_j6(6), Df_j3(3, 3), Df_j6(3, 6)
     type(matrix)              :: F_j3(3), F_j6(6), Ff_j3(3, 3), Ff_j6(3, 6)

!    open restart file
     call open_restart_file()

     r_i = r_f

     call pert_dens((/'EL'/), S, (/D/), (/F/), Df, Ff, w=(/-w/))


!    'EL'--'EL' equations
!    ======================

     r_j = r_e
     perturbation       = 'EL  '
     perturbation_print = 'polari    '

     allocate(Ef_j (r_i, r_j    ))
     allocate(Ef_jg(r_i, r_j, nc))

     Ef_j  = 0.0d0
     Ef_jg = 0.0d0

     restart_label = 'ff      '
     if (.not. find_restart_tensor(restart_label, (/-w, w/), size(Ef_j))) then

       call prop_energy((/perturbation/), (/D, Df/), Ef_j)
       call write_restart_tensor(restart_label, (/-w, w/), Ef_j)

     else

       call read_restart_tensor(Ef_j)
       write(*, *) 'successful restart from record '//restart_label

     end if

     call print_tensor((/r_i, r_j/), -Ef_j, perturbation_print, (/-w, w/))

     restart_label = 'ffg     '
     if (.not. find_restart_tensor(restart_label, (/-w, w, 0*w/), size(Ef_jg))) then
!      calculate, don't read from restart file

!      could be done like that:
!      call pert_dens((/perturbation/), S, (/D/), (/F/), D_j3, F_j3, w=(/w/))
!      call pert_dens((/'EL  ', perturbation/), &
!                       S,                      &
!                     (/D, Df, D_j3/),          &
!                     (/F, Ff, F_j3/),          &
!                       Df_j3,                  &
!                       Ff_j3,                  &
!                       w=(/-w, w/))

!      but here we save some time
       do j = 1, r_j
         D_j3(j) = dag(Df(j))
         F_j3(j) = dag(Ff(j))
      
         call pert_dens((/'EL  ', perturbation/), &
                          S,                      &
                        (/D, Df(1:j), D_j3(j)/),  &
                        (/F, Ff(1:j), F_j3(j)/),  &
                          Df_j3(1, j),            &
                          Ff_j3(1, j),            &
                          w=(/-w, w/),            &
                          c=(/1, j, j, j/))
      
         do i = 1, j - 1
           Df_j3(j, i) = dag(Df_j3(i, j))
           Ff_j3(j, i) = dag(Ff_j3(i, j))
         end do
       end do
       
       call prop_energy((/perturbation, 'GEO '/), (/D, Df/), Ef_jg)
       call prop_energy((/'EL ', 'GEO'/), (/D, D_j3/), Ef_jg, o=(/2/))
       call prop_energy((/'GEO'/), (/D, Df, D_j3, Df_j3/), Ef_jg, b=(/r_i, r_j/))
       
       do j = 1, r_j
         do i = 1, r_i
           DFD3(i, j) = Df_j3(i, j)*F*D         &
                      + D*F*Df_j3(i,j)          &
                      + D_j3(j)*Ff(i)*D         &
                      + D*Ff(i)*D_j3(j)         &
                      + Df(i)*F_j3(j)*D         &
                      + D*F_j3(j)*Df(i)         &
                      + D*Ff_j3(i, j)*D         &
                      + D_j3(j)*(F + w*S)*Df(i) &
                      + Df(i)*(F - w*S)*D_j3(j)
         end do
       end do
       
!      deallocate matrices
       D_j3  = 0
       Df_j3 = 0
       F_j3  = 0
       Ff_j3 = 0
       
       do k = 1, nc
         call prop_overlap((/'GEO'/), S, Sg, c=(/k/))
         do j = 1, r_j
           do i = 1, r_i
             Ef_jg(i, j, k) = Ef_jg(i, j, k) - tr(Sg(1), DFD3(i, j))
           end do
         end do
       end do
       Sg   = 0
       DFD3 = 0
       
       call write_restart_tensor(restart_label, (/-w, w, 0*w/), Ef_jg)

     else
!      don't calculate, read from restart file

       call read_restart_tensor(Ef_jg)
       write(*, *) 'successful restart from record '//restart_label

     end if

     call print_tensor((/r_i, r_j, nc/), -Ef_jg, 'd/dg '//perturbation_print, (/-w, w, 0*w/))

!    reverse signs of results
!    from quasi-energy derivative to property derivative
     Eff (1:r_i, 1:r_j   ) = -Ef_j (1:r_i, 1:r_j   )
     Effg(1:r_i, 1:r_j, :) = -Ef_jg(1:r_i, 1:r_j, :)

     deallocate(Ef_j)
     deallocate(Ef_jg)


!    'EL'--'MAGO' equations
!    ======================

     r_j = r_b
     perturbation       = 'MAGO'
     perturbation_print = 'G-prime   '

     allocate(Ef_j (r_i, r_j    ))
     allocate(Ef_jg(r_i, r_j, nc))

     Ef_j  = 0.0d0
     Ef_jg = 0.0d0

     restart_label = 'fb      '
     if (.not. find_restart_tensor(restart_label, (/-w, w/), size(Ef_j))) then

       call prop_energy((/perturbation/), (/D, Df/), Ef_j)
       call write_restart_tensor(restart_label, (/-w, w/), Ef_j)

     else

       call read_restart_tensor(Ef_j)
       write(*, *) 'successful restart from record '//restart_label

     end if

     call print_tensor((/r_i, r_j/), -Ef_j, perturbation_print, (/-w, w/))

     restart_label = 'fbg     '
     if (.not. find_restart_tensor(restart_label, (/-w, w, 0*w/), size(Ef_jg))) then
!      calculate, don't read from restart file

       call pert_dens((/perturbation/), S, (/D/), (/F/), D_j3, F_j3, w=(/w/))
       call pert_dens((/'EL  ', perturbation/), &
                        S,                      &
                      (/D, Df, D_j3/),          &
                      (/F, Ff, F_j3/),          &
                        Df_j3,                  &
                        Ff_j3,                  &
                        w=(/-w, w/))
       
       call prop_energy((/perturbation, 'GEO '/), (/D, Df/), Ef_jg)
       call prop_energy((/'EL ', 'GEO'/), (/D, D_j3/), Ef_jg, o=(/2/))
       call prop_energy((/'GEO'/), (/D, Df, D_j3, Df_j3/), Ef_jg, b=(/r_i, r_j/))
       
       do j = 1, r_j
         do i = 1, r_i
           DFD3(i, j) = Df_j3(i, j)*F*D         &
                      + D*F*Df_j3(i,j)          &
                      + D_j3(j)*Ff(i)*D         &
                      + D*Ff(i)*D_j3(j)         &
                      + Df(i)*F_j3(j)*D         &
                      + D*F_j3(j)*Df(i)         &
                      + D*Ff_j3(i, j)*D         &
                      + D_j3(j)*(F + w*S)*Df(i) &
                      + Df(i)*(F - w*S)*D_j3(j)
         end do
       end do
       
!      deallocate matrices
       D_j3  = 0
       Df_j3 = 0
       F_j3  = 0
       Ff_j3 = 0
       
       do k = 1, nc
         call prop_overlap((/'GEO'/), S, Sg, c=(/k/))
         do j = 1, r_j
           do i = 1, r_i
             Ef_jg(i, j, k) = Ef_jg(i, j, k) - tr(Sg(1), DFD3(i, j))
           end do
         end do
       end do
       Sg   = 0
       DFD3 = 0
       
       call write_restart_tensor(restart_label, (/-w, w, 0*w/), Ef_jg)

     else
!      don't calculate, read from restart file

       call read_restart_tensor(Ef_jg)
       write(*, *) 'successful restart from record '//restart_label

     end if

     call print_tensor((/r_i, r_j, nc/), -Ef_jg, 'd/dg '//perturbation_print, (/-w, w, 0*w/))

!    reverse signs of results
!    from quasi-energy derivative to property derivative
     Efb (1:r_i, 1:r_j   ) = -Ef_j (1:r_i, 1:r_j   )
     Efbg(1:r_i, 1:r_j, :) = -Ef_jg(1:r_i, 1:r_j, :)

     deallocate(Ef_j)
     deallocate(Ef_jg)


!    'EL'--'ELGR' equations
!    ======================

     r_j = r_q
     perturbation       = 'ELGR'
     perturbation_print = 'A-tensor  '

     allocate(Ef_j (r_i, r_j    ))
     allocate(Ef_jg(r_i, r_j, nc))

     Ef_j  = 0.0d0
     Ef_jg = 0.0d0

     restart_label = 'fq      '
     if (.not. find_restart_tensor(restart_label, (/-w, w/), size(Ef_j))) then

       call prop_energy((/perturbation/), (/D, Df/), Ef_j)
       call write_restart_tensor(restart_label, (/-w, w/), Ef_j)

     else

       call read_restart_tensor(Ef_j)
       write(*, *) 'successful restart from record '//restart_label

     end if

     call print_tensor((/r_i, r_j/), -Ef_j, perturbation_print, (/-w, w/))

     restart_label = 'fqg     '
     if (.not. find_restart_tensor(restart_label, (/-w, w, 0*w/), size(Ef_jg))) then
!      calculate, don't read from restart file

       call pert_dens((/perturbation/), S, (/D/), (/F/), D_j6, F_j6, w=(/w/))
       call pert_dens((/'EL  ', perturbation/), &
                        S,                      &
                      (/D, Df, D_j6/),          &
                      (/F, Ff, F_j6/),          &
                        Df_j6,                  &
                        Ff_j6,                  &
                        w=(/-w, w/))
       
       call prop_energy((/perturbation, 'GEO '/), (/D, Df/), Ef_jg)
       call prop_energy((/'EL ', 'GEO'/), (/D, D_j6/), Ef_jg, o=(/2/))
       call prop_energy((/'GEO'/), (/D, Df, D_j6, Df_j6/), Ef_jg, b=(/r_i, r_j/))
       
       do j = 1, r_j
         do i = 1, r_i
           DFD6(i, j) = Df_j6(i, j)*F*D         &
                      + D*F*Df_j6(i,j)          &
                      + D_j6(j)*Ff(i)*D         &
                      + D*Ff(i)*D_j6(j)         &
                      + Df(i)*F_j6(j)*D         &
                      + D*F_j6(j)*Df(i)         &
                      + D*Ff_j6(i, j)*D         &
                      + D_j6(j)*(F + w*S)*Df(i) &
                      + Df(i)*(F - w*S)*D_j6(j)
         end do
       end do
       
!      deallocate matrices
       D_j6  = 0
       Df_j6 = 0
       F_j6  = 0
       Ff_j6 = 0
       
       do k = 1, nc
         call prop_overlap((/'GEO'/), S, Sg, c=(/k/))
         do j = 1, r_j
           do i = 1, r_i
             Ef_jg(i, j, k) = Ef_jg(i, j, k) - tr(Sg(1), DFD6(i, j))
           end do
         end do
       end do
       Sg   = 0
       DFD6 = 0
       
       call write_restart_tensor(restart_label, (/-w, w, 0*w/), Ef_jg)

     else
!      don't calculate, read from restart file

       call read_restart_tensor(Ef_jg)
       write(*, *) 'successful restart from record '//restart_label

     end if

     call print_tensor((/r_i, r_j, nc/), -Ef_jg, 'd/dg '//perturbation_print, (/-w, w, 0*w/))

!    reverse signs of results
!    from quasi-energy derivative to property derivative
     Efq (1:r_i, 1:r_j   ) = -Ef_j (1:r_i, 1:r_j   )
     Efqg(1:r_i, 1:r_j, :) = -Ef_jg(1:r_i, 1:r_j, :)

     deallocate(Ef_j)
     deallocate(Ef_jg)


!    finally deallocate rest
     Df  = 0
     Ff  = 0

!    close restart file
     call close_restart_file()

   end subroutine

   subroutine cars_shyp_polgra(S,D,F,nc,w,Ef,Eff,Effff,Effg)
   !For frequency w, Eff returns polarizability alpha(-w,w), Effg polarizability
   !gradient d/dR alpha(-w,w), and Effff hyperpolarizability gamma(-w,w,-w,w)
      type(matrix),intent(in)  :: S,D,F
      integer,     intent(in)  :: nc
      complex(8),  intent(in)  :: w
      complex(8),  intent(out) ::Ef(3),Eff(3,3),Effff(3,3,3,3),Effg(3,3,nc)
      type(matrix) :: De(3),Df(3),Def(3,3),Dee(3,3),Lf(3),Sg(1),A
      type(matrix) :: Fe(3),Ff(3),Fef(3,3),Fee(3,3),Xf(3),DFD(3,3)
      integer      :: i,j,k,l
      !dipole moment Ef=-mu
      Ef=0; call prop_energy((/'EL'/),(/D/),Ef)
      call print_tensor((/3/),-Ef,'Mu = -Ef')
      !electric equations, Eff=-alpha
      call pert_dens((/'EL'/),S,(/D/),(/F/),Df,Ff,w=(/w/))
      do i=1,3; De(i)=dag(Df(i)); Fe(i)=dag(Ff(i)); end do
      Eff=0; call prop_energy((/'EL'/),(/D,De/),Eff)
      call print_tensor((/3,3/),-Eff,'Alpha = -Eff',(/-w,w/))
      !---------------------------------
      !electric-electric equations Def, e_-w, f_+w
      do j=1,3
         call pert_dens((/'EL','EL'/),S,(/D,De(1:j),Df(j)/),(/F,Fe(1:j),Ff(j)/), &
                          Def(1,j),Fef(1,j),w=(/-w,w/),c=(/1,j,j,j/))
         do i=1,j-1; Def(j,i)=dag(Def(i,j)); Fef(j,i)=dag(Fef(i,j)); end do
      end do
      !polarizability gradient
      Effg=0; call prop_energy((/'EL ','GEO'/),(/D,De/),Effg)
      call prop_energy((/'EL ','GEO'/),(/D,Df/),Effg,o=(/2/))
      call prop_energy((/'GEO'/),(/D,De,Df,Def/),Effg,b=(/3,3/))
      do j=1,3; do i=1,3
         DFD(i,j) = Def(i,j)*F*D +Df(j)*Fe(i)*D+Df(j)*(F+w*S)*De(i) &
                  + De(i)*Ff(j)*D+D*Fef(i,j)*D +D*Ff(j)*De(i)       &
                  + De(i)*(F-w*S)*Df(j)+D*Fe(i)*Df(j)+D*F*Def(i,j)
      end do; end do
      do k=1,nc
         call prop_overlap((/'GEO'/),S,Sg,c=(/k/))
         do j=1,3; do i=1,3
            Effg(i,j,k) = Effg(i,j,k) - tr(Sg(1),DFD(i,j))
         end do; end do; Sg(1)=0
      end do
      DFD = 0
      call print_tensor((/3,3,nc/),-Effg,'d/dg Alpha = -Effg',(/-w,w,0*w/))
      !---------------------------------
      !electric-electric equations Dee, e_-w
      do k=1,3
         call pert_dens((/'EL','EL'/),S,(/D,De(1:k),De(k)/),(/F,Fe(1:k),Fe(k)/), &
                          Dee(1,k),Fee(1,k),w=(/-w,-w/),c=(/1,k,k,k/))
         do i=1,k-1; Dee(k,i)=Dee(i,k); Fee(k,i)=Fee(i,k); end do
      end do
      !second hyperpolarizability
      Effff=0 !no energy or reorth contributions
      !stationarity multiplier contribution
      do l=1,3; Xf(l) = Df(l)*S*D-D*S*Df(l); end do
      do k=1,3; do j=1,3; do i=1,3
         A = (Fef(i,j)*De(k)+Fee(i,k)*Df(j)+Fef(k,j)*De(i)    &
             +Fe(k)*Def(i,j)+Ff(j)*Dee(i,k)+Fe(i)*Def(k,j))*S &
           - S*(Def(i,j)*Fe(k)+Dee(i,k)*Ff(j)+Def(k,j)*Fe(i)  &
               +De(k)*Fef(i,j)+Df(j)*Fee(i,k)+De(i)*Fef(k,j))
         do l=1,3; Effff(i,j,k,l) = Effff(i,j,k,l)-tr(Xf(l),A); end do; A=0
      end do; end do; end do
      Xf = 0
      !idempotency multiplier contribution
      do l=1,3; Lf(l) = Ff(l)*D*S+S*D*Ff(l)-Ff(l); end do
      do k=1,3; do j=1,3; do i=1,3
         A = Def(k,j)*S*De(i)+Dee(i,k)*S*Df(j)+Def(i,j)*S*De(k) &
           + De(i)*S*Def(k,j)+Df(j)*S*Dee(i,k)+De(k)*S*Def(i,j)
         do l=1,3; Effff(i,j,k,l) = Effff(i,j,k,l)-tr(Lf(l),A); end do; A=0
      end do; end do; end do; Lf=0
      call print_tensor((/3,3,3,3/),-Effff,'Gamma = -Effff',(/-w,w,-w,w/))
      !------------------------------
      !delete matrices, reverse signs (from quasi-energy- to property-derivatives)
      De=0; Fe=0; Df=0; Ff=0; Def=0; Fef=0; Dee=0; Fee=0
      Ef=-Ef; Eff=-Eff; Effff=-Effff; Effg=-Effg
   end subroutine


  subroutine dg_energy(S, D, F, n)

!   ----------------------------------------------------------------------------
    type(matrix)        :: S, D, F
    integer, intent(in) :: n
!   ----------------------------------------------------------------------------
    type(matrix)        :: W, Sg(1)
    integer             :: i, iatom, ix, iy, iz
    complex(8)          :: g(n)
!   ----------------------------------------------------------------------------

    g = 0.0d0

    call prop_energy((/'GEO'/), (/D/), g)

    W = D*F*D

    do i = 1, n

      call prop_overlap((/'GEO'/), S, Sg, c=(/i/))
      g(i) = g(i) - dot(W, Sg(1))
      Sg(1) = 0

    end do

    W = 0

    call header('molecular gradient', -1)

    do iatom = 1, nr_atoms
      ix = 3*(iatom - 1) + 1
      iy = 3*(iatom - 1) + 2
      iz = 3*(iatom - 1) + 3
      write(*, '(a6, 3e26.16)') sym_dep_nuc_name(iatom), real(g(ix)), real(g(iy)), real(g(iz))
    end do

    call prsymb(unit_output, '-', 46, 0)

  end subroutine



! subroutine DipGradi(S,D,F,nc)
!     type(matrix):: S,D,F,Df(3),Ff(3),DFD(3),Sg(1)
!     real(8):: Efg(3,nc); integer:: nc,i,j
!     !dipole moment
!     Efg(:,1)=0; call prop_energy((/'EL'/),(/D/),Efg)
!     i = max(1,ceiling(log10(maxval(abs(Efg(:,1))))))
!     print '("Dipole = -Ef"/3f19.<16-i>/)',-Efg(:,1)
!     !first order equations
!     call pert_dens((/'EL'/),S,(/D/),(/F/),Df,Ff,w=(/0d0/))
!     Efg=0; call prop_energy((/'EL','GEO'/),(/D/),Efg)
!     call prop_energy((/'GEO'/),(/D,Df/),Efg)
!     do i=1,3; DFD(i) = Df(i)*F*D+D*Ff(i)*D+D*F*Df(i); end do
!     do j=1,nc
!         call prop_overlap((/'GEO'/),S,Sg,c=(/j/))
!         do i=1,3
!             Efg(i,j) = Efg(i,j) - dot(DFD(i),Sg(1))
!         end do; Sg(1)=0
!     end do; do i=1,3; DFD(i)=0; end do
!     i = max(1,ceiling(log10(maxval(abs(Efg)))))
!     print '("d/dg Dipole = -Efg w:"2f6.3/<nc>(3f19.<16-i>/))',0d0,0d0,-Efg
!     do i=1,3; Df(i)=0; Ff(i)=0; end do
! end subroutine


end module
