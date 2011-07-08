  program hessian
! This program reads information extracted from a GOSCI run and
! sets up the electronic Hessian
  integer :: i,j,ii,ia,ji,ja,stat
  integer :: indi,indj,indk,indl,nxop
  integer :: norb,nocc,nvrt,nvar
  integer, allocatable :: jxop(:,:)
  real(8), allocatable :: orbeig(:)
  real(8), allocatable :: twoint(:,:,:,:)
  real(8), allocatable :: ablock(:,:),bblock(:,:)

! We start by getting the information extracted from GOSCI
! ========================================================
  open (5, Form='FORMATTED')
! First line contains total number of orbitals and occupied ones
! --------------------------------------------------------------
  read(5,*) norb,nocc
  nvrt=norb-nocc
  write(6,*) norb,nocc
! Next section is a list of orbital energies
! ------------------------------------------
  read(5,*)
  allocate (orbeig(norb))
  do i = 1,norb
    read(5,*) orbeig(i)
  enddo
  write(6,*) (orbeig(i),i=1,norb)

! Final section is the list of non-zero two electron integrals
! ------------------------------------------------------------
  allocate (twoint(norb,norb,norb,norb))
  twoint = 0
  read(5,*)
  stat = 0
  do while (stat.eq.0)
    read(5,*,iostat=stat) indi,indj,indk,indl,twoint(indi,indj,indk,indl)
    write(6,*) indi,indj,indk,indl,twoint(indi,indj,indk,indl)

  enddo
  close(5,status='KEEP')

! We set up orbital rotation indices 
! ========================================================================
  nvar=nvrt*nocc
  allocate (jxop(2,nvar))
  i = 0
  do ia = nocc+1,norb
    do ii = 1,nocc
      i = i+1
      jxop(1,i)=ia
      jxop(2,i)=ii
      write(6,*) i,ia,ii
    enddo
  enddo
! Set up A-block of Hessian (coupling terms only)
! ===============================================
  allocate (ablock(nvar,nvar))
  do j=1,nvar
    ja=jxop(1,j)
    ji=jxop(2,j)
    do i=1,nvar
      ia=jxop(1,i)
      ii=jxop(2,i)
      ablock(i,j)=twoint(ia,ii,ji,ja)-twoint(ia,ja,ji,ii)
    enddo
  enddo
  write(6,'(/A/)') '*** A block ***'
  do i = 1,nvar
    write(6,*) ablock(i,:)
  enddo

! Set up B-block of Hessian (coupling terms only)
! ===============================================
  allocate (bblock(nvar,nvar))
  do j=1,nvar
    ja=jxop(1,j)
    ji=jxop(2,j)
    do i=1,nvar
      ia=jxop(1,i)
      ii=jxop(2,i)
      bblock(i,j)=twoint(ia,ii,ja,ji)-twoint(ia,ji,ja,ii)
    enddo
  enddo
  write(6,'(/A/)') '*** B block ***'
  do i = 1,nvar
    write(6,*) bblock(i,:)
  enddo

  end
