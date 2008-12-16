module HEFFModule

   use Parameters
   use Options
   use OrbMod
   use IOmodule
   use VacMod
   !use SymBlock
   use SymBlocks
   use s2matele
   use TAmps
   use TIndexMod
   use GMSMod
   implicit none
   private

   ! Public Variables

   double precision, allocatable, dimension(:,:) :: HEFF
   double precision, allocatable, dimension(:,:) :: HEF0
   double precision, allocatable, dimension(:,:) :: HEFP
   public :: HEFF, HEF0, HEFP

   ! Private Vairables

   logical :: Init = .FALSE.

   ! Arrays for Calc_HEFF_ZerBody routine
   integer, allocatable, dimension(:), target :: ResP, ResQ
   double precision, allocatable, dimension(:) :: TMP_NPMOSQ
   double precision, allocatable, dimension(:,:) :: FINTHP

   ! Effective Hamiltonian block
   type (op_symblock), allocatable, target :: heff_op(:)

   ! Other arrays
   integer , allocatable, dimension(:,:) :: SHEF
   double precision, allocatable, dimension(:,:) :: HEFF_MATRIX

   ! Public Routines

   public :: InitHEFF, EndHEFF
   public :: DrvHEFF, DiagHEFF, PrintHEFF

contains

subroutine InitHEFF

   implicit none
   integer :: i, p, q, pq
   integer :: vac
   integer :: IntJunk(2)
   double precision :: DbleJunk
   integer :: isym_block, nvac

   if (Init) then
      write(UOUT,*) 'Attempt to Reinitialize HEFF module'
      call flush(UOUT)
      call abort
   endif

   if (NSTATE <= 0) then
      write(UOUT,*) 'Illegal Value of NSTATE'
      call flush(UOUT)
      call abort
   endif

   Init  = .TRUE.

   ! ... Allocations and Initializations for Calc_HEFF_ZerBody routine ...

   if (.not.allocated(TMP_NPMOSQ)) allocate(TMP_NPMOSQ(NPMOSQ))
   if (.not.allocated(ResP))       allocate(ResP(MAX(NPMOSQ,NHMOSQ)))
   if (.not.allocated(ResQ))       allocate(ResQ(MAX(NPMOSQ,NHMOSQ)))
   if (.not.allocated(FINTHP))     allocate(FINTHP(NPMO,NHMO))

   pq = 0
   do p = 2, MAX(NPMO,NHMO)
   do q = 1, (p-1)
      pq = pq + 1
      ResP(pq) = p
      ResQ(pq) = q
   enddo
   enddo

   ! ... Allocate public variables ...

   if (.not.allocated(HEFF))        allocate(HEFF(NSTATE,NSTATE))
   if (.not.allocated(HEF0))        allocate(HEF0(NSTATE,NSTATE))
   if (.not.allocated(HEFP))        allocate(HEFP(NSTATE,NSTATE))

   ! ... Allocate private arrays ...

   if (.not.allocated(SHEF))        allocate(SHEF(NSTATE,NSTATE))
   if (.not.allocated(HEFF_MATRIX)) allocate(HEFF_MATRIX(NSTATE,NSTATE))

   ! ... Allocate effective hamiltonian

   allocate(heff_op(NSYM_BLOCK))

   ! ... run over symmetry blocks ...
   do isym_block = 1, NSYM_BLOCK

      nvac = symblock_array(isym_block)%nvac
      allocate(heff_op(isym_block)%op(nvac,nvac))
      allocate(heff_op(isym_block)%eval(nvac))
      allocate(heff_op(isym_block)%rvec(nvac,nvac))
      allocate(heff_op(isym_block)%lvec(nvac,nvac))
      allocate(heff_op(isym_block)%s2p1(nvac))

   enddo

   ! ... Build Sign Matrix for Effective Hamiltonian  ...
   ! ... Effective Hamiltonian needs to be multiplied ...
   ! ... with this matrix (not matrix multiplication) ...
   ! ... before carrying out diagonalization          ...

   call SignHEFF

   ! ... Build HEF0

   call BuildHEF0

   ! ... Store HEF0 if it is fresh run ...

   if (.not.OPTION(RESTART)) then
      rewind(UHEFF)
      write(UHEFF) HEF0
   endif

   return

end subroutine InitHEFF

subroutine EndHEFF

   implicit none
   integer isym_block, nvac

   if (.not.Init) then
      write(UOUT,*) 'EndHEFF error'
      call flush(UOUT)
      call abort
   endif

   ! ... Deallocate ...

   if (allocated(TMP_NPMOSQ)) deallocate(TMP_NPMOSQ)
   if (allocated(ResP))       deallocate(ResP)
   if (allocated(ResQ))       deallocate(ResQ)
   if (allocated(FINTHP))     deallocate(FINTHP)

   if (allocated(HEFF))       deallocate(HEFF)
   if (allocated(HEF0))       deallocate(HEF0)
   if (allocated(HEFP))       deallocate(HEFP)

   if (allocated(SHEF))        deallocate(SHEF)
   if (allocated(HEFF_MATRIX)) deallocate(HEFF_MATRIX)

   do isym_block = 1, NSYM_BLOCK
      nvac = symblock_array(isym_block)%nvac
      deallocate(heff_op(isym_block)%op)
      deallocate(heff_op(isym_block)%eval)
      deallocate(heff_op(isym_block)%rvec)
      deallocate(heff_op(isym_block)%lvec)
      deallocate(heff_op(isym_block)%s2p1)
   enddo
   deallocate(heff_op)

   return

end subroutine EndHEFF

subroutine DrvHEFF(vac)

   use TimingMod
   implicit none
   integer, intent(in) :: vac

   if (.not.Init) then
      write(UOUT,*) 'HEFF module not Initialized'
      call flush(UOUT)
      call abort
   endif

   call Calc_HEFF_ZerBody(TimeTaken(:,HEFFEC),vac)
   call Calc_HEFF_OneBody(TimeTaken(:,HEFFEC),vac)
   call Calc_HEFF_TwoBody(TimeTaken(:,HEFFEC),vac)

   return

end subroutine DrvHEFF

! This routine calculates the Zero-body part of
! effective Hamiltonian for the current vacuum.
! It is just the diagonal element.

subroutine Calc_HEFF_ZerBody(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd
   integer :: i, j, k, l
   integer :: kl, ij, ix, jx
   integer :: ti
   integer, dimension(:), pointer :: aar, bar
   integer,          dimension(:), pointer :: VIND
   double precision, dimension(:), pointer :: VINT
   real :: dtim, dtime, tarray(2)
!   external :: dtime

   dtim = dtime(tarray)

   ! Add the vaccum expectation value HEF0 part

   HEFF(1,vac) = HEFF(1,vac) + HEF0(1,vac)

   ! + sum(i>j,a>b) [ <i,j||a,b>*(t(i,j->a,b)+t(i->a)*t(j->b)-t(i->b)*t(i->a)) ]

   call ReadFile(VHHPP)

   TMP_NPMOSQ = dble(0)
   aar => ResP
   bar => ResQ
   ij = 0
   do i = 2, NHMO
   do j = 1, (i-1)
      ij = ij + 1
      ti = TIND(ij)
      ix = (i-1)*NPMO
      jx = (j-1)*NPMO
      TMP_NPMOSQ = TOLD(ti+1:ti+NPMOSQ) + TOLD(ix+aar)*TOLD(jx+bar) &
                                        - TOLD(ix+bar)*TOLD(jx+aar)
      EndBufInd = 0
      kl = 0
      do k = 2, NHMO
      do l = 1, (k-1)
         kl = kl + 1
         StartBufInd = EndBufInd + 1
         EndBufInd   = EndBufInd + NZIntRec(kl)
         if ( EndBufInd < StartBufInd ) cycle
         if (kl/=ij) cycle
         VINT => IntBuf(StartBufInd:EndBufInd)
         VIND => IndBuf(StartBufInd:EndBufInd)
         HEFF(1,vac) = HEFF(1,vac) + DOT_PRODUCT(VINT,TMP_NPMOSQ(VIND))
      enddo
      enddo
   enddo
   enddo

   ! + sum(i,a)[<i|f|a>*t1(i->a)]

   read(UFHP,rec=vac) FINTHP
   do i = 1, NHMO
      ix = (i-1)*NPMO
      HEFF(1,vac) = HEFF(1,vac) + DOT_PRODUCT(FINTHP(:,i),TOLD(ix+1:ix+NPMO))
   enddo

   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Calc_HEFF_ZerBody

! This routines calculates one-body part of
! effective hamiltonian for the current vaccum.
! The contributions are in calculated in
! certain elements of TNEW when direct part
! calculation is carried out. Hence, they are
! simply collected into HEFF matrix.

subroutine Calc_HEFF_OneBody(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: iexc, tind
   integer :: lamda, gaama
   real :: dtim, dtime, tarray(2)
!   external :: dtime

   dtim = dtime(tarray)

   iexc = 1
   do lamda = 1, NVMO
      if (.not.IsHole(lamda+NCMO,vac)) cycle
      do gaama = 1, NVMO
         if (.not.IsPtcl(gaama,vac)) cycle
         iexc = iexc + 1
         tind = TINDEX(lamda+NCMO,gaama)
         HEFF(iexc,vac) = TNEW(tind)
         TNEW(tind) = dble(0)
       enddo
   enddo

   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Calc_HEFF_OneBody

! This routines calculates two-body part of
! effective hamiltonian for the current vaccum.
! The contributions are in calculated in
! certain elements of TNEW when direct part
! calculation is carried out. Hence, they are
! simply collected into HEFF matrix.

subroutine Calc_HEFF_TwoBody(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: iexc, tind
   integer :: lamda, gaama, sigma, delta
   real :: dtim, dtime, tarray(2)
!   external :: dtime

   dtim = dtime(tarray)

   iexc = 1 + NPTCL*NHOLE
   do lamda = 1, NVMO
      if (.not.IsHole(lamda+NCMO,vac)) cycle
      do sigma = 1, (lamda-1)
         if (.not.IsHole(sigma+NCMO,vac)) cycle
         ! lamda > sigma
         do gaama = 1, NVMO
            if (.not.IsPtcl(gaama,vac)) cycle
            do delta = 1, (gaama-1)
               if (.not.IsPtcl(delta,vac)) cycle
               ! gaama > delta
               iexc = iexc + 1
               tind = TINDEX(lamda+NCMO,sigma+NCMO,gaama,delta)
               HEFF(iexc,vac) = TNEW(tind)
               TNEW(tind) = dble(0)
            enddo
         enddo
      enddo
   enddo

   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Calc_HEFF_TwoBody

subroutine SignHEFF

   implicit none
   integer :: vac
   integer :: i, j, k, l, a, b, c, d
   integer :: nh, nperm, exc, ival
   integer, dimension(NHOLE) :: CurVac, ExctVac

   do vac = 1, NSTATE

      ! ... Construct Lexically Ordered               ...
      ! ... Determinant for Current Vaccum            ...
      ! ... A Lexically Ordered Determinant has all   ...
      ! ... Orbital indices in increasing order       ...

      nh = 0
      do ival = 1, NVMO
         if (IsActHole(ival,vac)) then
            nh = nh + 1
            CurVac(nh) = ival
         endif
      enddo
      if (nh/=NHOLE) stop 'Error in SignHEFF'

      ! ... CurVac is Lexically Ordered Determinant   ...
      ! ... Excitations out of this Determinant are   ...
      ! ... Constructed. The resulting determinant is ...
      ! ... not, in general, lexically ordered. The   ...
      ! ... number of permutations required to bring  ...
      ! ... this determinant to lexical order         ...
      ! ... determines the sign of the HEFF matrix    ...
      ! ... connecting these two vaccums              ...
      ! ...                                           ...
      ! ...     <ExctVac|HEFF|CurVac>                 ...
      ! ...                                           ...
      ! ... If the no of permutations is,             ...
      ! ...   Even => positive                        ...
      ! ...   Odd  => negative                        ...

      exc = 1

      ! ... No Excitation ...

      SHEF(exc,vac) = +1

      ! ... Single Excitations ...

      do i = 1, NHOLE
      do a = 1, NVMO
         if (.not.IsActPtcl(a,vac)) cycle
         exc = exc + 1
         ExctVac = CurVac
         ExctVac(i) = a
         ! ... Now Find the number of Permutations ...
         ! ... Required to Lexically Order this    ...
         nperm = NPermLexOrd(ExctVac)
         if (MOD(nperm,2)==0) then
            SHEF(exc,vac) = +1
         else
            SHEF(exc,vac) = -1
         endif
      enddo
      enddo

      ! ... Double Excitations ...

      do i = 1, NHOLE
      do j = 1, (i-1)
         do a = 1, NVMO
            if (.not.IsActPtcl(a,vac)) cycle
            do b = 1, (a-1)
               if (.not.IsActPtcl(b,vac)) cycle
               exc = exc + 1
               ExctVac = CurVac
               ExctVac(i) = a
               ExctVac(j) = b
               ! ... Now Find the number of Permutations ...
               ! ... Required to Lexically Order this    ...
               nperm = NPermLexOrd(ExctVac)
               if (MOD(nperm,2)==0) then
                  SHEF(exc,vac) = +1
               else
                  SHEF(exc,vac) = -1
               endif
            enddo
         enddo
      enddo
      enddo

      ! ... Triple Excitations ...

      do i = 1, NHOLE
      do j = 1, (i-1)
      do k = 1, (j-1)
         do a = 1, NVMO
            if (.not.IsActPtcl(a,vac)) cycle
            do b = 1, (a-1)
               if (.not.IsActPtcl(b,vac)) cycle
               do c = 1, (b-1)
                  if (.not.IsActPtcl(c,vac)) cycle
                  exc = exc + 1
                  ExctVac = CurVac
                  ExctVac(i) = a
                  ExctVac(j) = b
                  ExctVac(k) = c
                  ! ... Now Find the number of Permutations ...
                  ! ... Required to Lexically Order this    ...
                  nperm = NPermLexOrd(ExctVac)
                  if (MOD(nperm,2)==0) then
                     SHEF(exc,vac) = +1
                  else
                     SHEF(exc,vac) = -1
                  endif
               enddo
            enddo
         enddo
      enddo
      enddo
      enddo

      ! ... Quadruple Excitations ...

      do i = 1, NHOLE
      do j = 1, (i-1)
      do k = 1, (j-1)
      do l = 1, (k-1)
         do a = 1, NVMO
            if (.not.IsActPtcl(a,vac)) cycle
            do b = 1, (a-1)
               if (.not.IsActPtcl(b,vac)) cycle
               do c = 1, (b-1)
                  if (.not.IsActPtcl(c,vac)) cycle
                  do d = 1, (c-1)
                     if (.not.IsActPtcl(d,vac)) cycle
                     exc = exc + 1
                     ExctVac = CurVac
                     ExctVac(i) = a
                     ExctVac(j) = b
                     ExctVac(k) = c
                     ExctVac(l) = d
                     ! ... Now Find the number of Permutations ...
                     ! ... Required to Lexically Order this    ...
                     nperm = NPermLexOrd(ExctVac)
                     if (MOD(nperm,2)==0) then
                        SHEF(exc,vac) = +1
                     else
                        SHEF(exc,vac) = -1
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      enddo
      enddo
      enddo

      if (exc/=NSTATE) stop 'Error in SignHEFF'

   enddo

   return

end subroutine SignHEFF

! This functions returns 
! the number of Permutations
! required to bring a given
! Excitation to lexical order

function NPermLexOrd(Exct)

   implicit none
   integer, intent(in) :: Exct(NHOLE)
   integer :: NPermLexOrd
   integer :: LexOrd(NHOLE)
   integer :: i, j, k

   NPermLexOrd = 0
   LexOrd = Exct
   do i = 1, NHOLE
      do j = (i+1), NHOLE
         if (LexOrd(i) > LexOrd(j)) then
            k = LexOrd(i)
            LexOrd(i) = LexOrd(j)
            LexOrd(j) = k
            NPermLexOrd = NPermLexOrd + 1
         endif
      enddo
   enddo

   return

end function NPermLexOrd

!
!   This subroutine builds the vaccum expectation
!   value of the effective hamiltonian.

!   vaccum-vaccum expectation value of
!   effective hamiltonian.

!   HEF0(vac) = sigma(OrbEn(i),i is hole belongs to vac)
!               + 1/2 * sigma(<i,j||i,j>, i and j are holes
!                                         belonging to vac)
!               - sigma(<i,j||i,j>, i is hole in vac
!                                    j is hole in Ref vac)

subroutine BuildHEF0

   implicit none
   integer :: vac
   integer :: i, j, ij, kl
   integer :: StartBufInd, EndBufInd
   integer         , dimension(:), pointer :: VIND
   double precision, dimension(:), pointer :: VINT

   HEF0 = dble(0)

   ! orbital energy part add nuclear energy part also.

   do vac = 1, NSTATE
      HEF0(1,vac) = HEF0(1,vac) + NucEn
      do i = 1, NHMO
         if (IsHole(i,vac)) HEF0(1,vac) = HEF0(1,vac) + OrbEn(i)
      enddo
   enddo

   call ReadFile(VHHHH) ! ... read here two-el.integrals 
   do vac = 1, NSTATE
      EndBufInd = 0
      ij = 0
      do i = 2, NHMO
      do j = 1, (i-1)
         ij = ij + 1
         StartBufInd = EndBufInd + 1
         EndBufInd   = EndBufInd + NZIntRec(ij)
         if ( EndBufInd < StartBufInd ) cycle
         VINT => IntBuf(StartBufInd:EndBufInd)
         VIND => IndBuf(StartBufInd:EndBufInd)
         do kl = 1, SIZE(VIND)
            if (VIND(kl)/=ij) cycle
            if (IsHole(i,vac).and.IsHole(j,vac)) HEF0(1,vac) = HEF0(1,vac) + VINT(kl)
            if (IsHole(i,vac).and.IsHole(j,  1)) HEF0(1,vac) = HEF0(1,vac) - VINT(kl)
            if (IsHole(i,  1).and.IsHole(j,vac)) HEF0(1,vac) = HEF0(1,vac) - VINT(kl)
         enddo
      enddo
      enddo
   enddo

   return

end subroutine BuildHEF0


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine PrintHEFF
!----------------------------------------------------------
!
!   This subroutine prints HEFF.
!
!----------------------------------------------------------

   implicit none

   integer ir
   integer isym_block
   integer ms2, irp
   integer nvac

   write(UOUT,9011)

9011  format(/, &
      "**************************************************************",/, &
      "Effective Hamiltonian Matrix, its eigenvalues and eigenvectors",/, &
      "**************************************************************") 
            
   write(UOUT,9012) NSYM_BLOCK

9012  format("number of symmetry blocks = ",i3)

   do isym_block = 1, NSYM_BLOCK

      ms2  =  symblock_array(isym_block)%ms2
      irp  =  symblock_array(isym_block)%irp
      nvac =  symblock_array(isym_block)%nvac

      write(UOUT,9021) isym_block, ms2, irp, irrep_label(irp), nvac

      call prntmat("effective hamiltonian block", &
                   UOUT,nvac,nvac,heff_op(isym_block)%op,nvac, &
                   symblock_array(isym_block)%vacs,symblock_array(isym_block)%vacs)

      write(UOUT,9022)
      do ir = 1, nvac
         write(UOUT,9023) ir, ms2, irp, irrep_label(irp), heff_op(isym_block)%s2p1(ir), &
                          heff_op(isym_block)%eval(ir)
      enddo
      write(UOUT,9024)

      call prntmat("right eigenvectors (bi-orthogonalized)", &
                   UOUT,nvac,nvac,heff_op(isym_block)%rvec,nvac, &
                   symblock_array(isym_block)%vacs,symblock_array(isym_block)%vacs)
      call prntmat("left  eigenvectors (bi-orthogonalized)", &
                   UOUT,nvac,nvac,heff_op(isym_block)%lvec,nvac, &
                   symblock_array(isym_block)%vacs,symblock_array(isym_block)%vacs)

   enddo

9021  format(/,&
             "***********************************************************************",/,&
             "Symmetry Block = ",i3,1x,"Ms2 = ",i3,1x,"Irrep = ",i1,1x,"(",a4,")",1x,"Number of vacua = ",i5,/, &
             "***********************************************************************"/)
            
9022  format(/,&
             "**************************************",/, &
             "root ms2 irp label 2s+1     energy    ",/, &
             "**************************************"/)
            
9023  format(i4,1x,i3,1x,1x,i1,1x,1x,a4,1x,i3,1x,1x,F14.9,/)
9024  format("**************************************",/)

   return

end subroutine PrintHEFF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine DiagHEFF
!-----------------------------------------------------------------
!
!             Diagonalize the effective Hamiltonian
!
!  Called (3 times) from "hsmrcc.F90"/HSMRCC
!
!-----------------------------------------------------------------

   implicit none
   integer :: RightVac, LeftVac, Exct
   integer :: isym_block
   integer :: ms2, irp, nvac, nvac_ret
   integer :: ir, il, vacr, vacl
   double precision, allocatable :: tmpv(:)
   double precision :: dtmp
   double precision :: ddot


   ! ... Construct HEFF with proper sign  ...

   HEFF_MATRIX = dble(0)
   do RightVac = 1, NSTATE
      do Exct = 1, NSTATE
         ! LeftVac is the vaccum number
         ! of the Exctn with respect to RightVac.
         ! Sign for non-lexically ordered
         ! Excitations are taken here.
         LeftVac = VacTable(Exct,RightVac)
         HEFF_MATRIX(LeftVac,RightVac) = DBLE(SHEF(Exct,RightVac))*HEFF(Exct,RightVac)
      enddo
   enddo

   do isym_block = 1, NSYM_BLOCK

      ms2  =  symblock_array(isym_block)%ms2
      irp  =  symblock_array(isym_block)%irp
      nvac =  symblock_array(isym_block)%nvac

      do ir = 1, nvac
         vacr = symblock_array(isym_block)%vacs(ir)
         do il = 1, nvac
            vacl = symblock_array(isym_block)%vacs(il)
            heff_op(isym_block)%op(il,ir) = HEFF_MATRIX(vacl,vacr)
         enddo
      enddo

      ! ... for diagonalization use better routine ...
      ! ... obtain a set of biorthogonal vectors   ...
      ! ... abort if complex eigenvalues are found ...

      call diagonalize_hscc             &
           (                            &
              .false.,                  &
              nvac, nvac,               &
              heff_op(isym_block)%op,   &
              nvac_ret,                 &
              heff_op(isym_block)%eval, &
              heff_op(isym_block)%rvec, &
              heff_op(isym_block)%lvec, &
              .true., .false., .true.   &
           )

      allocate(tmpv(nvac))
      do ir = 1, nvac
         call dgemv('t',nvac,nvac,1.0d00,symblock_array(isym_block)%s2mat,nvac, &
                    heff_op(isym_block)%rvec(1,ir),1,0.0d00,tmpv,1)
         dtmp = ddot(nvac,heff_op(isym_block)%lvec(1,ir),1,tmpv,1)
         heff_op(isym_block)%s2p1(ir) = two_s_plus1(dtmp)
      enddo
      deallocate(tmpv)

   enddo

   return

end subroutine DiagHEFF

!include "tind.h"

end module HEFFModule
