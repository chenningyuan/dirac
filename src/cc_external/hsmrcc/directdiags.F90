! This module includes routines for
! computing HSMRCC direct diagrams.

module DirectDiags

   use Parameters
   use IOModule
   use TAmps
   implicit none
   private
   public :: InitDirect
   public :: DrvDirect
   public :: EndDirect

   logical :: Init = .FALSE.
   integer,          dimension(:), pointer :: VIND
   double precision, dimension(:), pointer :: VINT
   double precision, dimension(:), pointer :: TOLDP, TNEWP

   double precision, dimension(:)  , pointer :: TMP
   double precision, dimension(:,:), pointer :: TMP1, TMP2, TMP3
   double precision, dimension(:,:), pointer :: TMP4, TMP5
   double precision, dimension(:,:), pointer :: FINT

   integer         , allocatable, dimension(:)  , target :: ResP, ResQ
   double precision, allocatable, dimension(:)  , target :: TMP_NPSQ
   double precision, allocatable, dimension(:,:), target :: TMP_HH
   double precision, allocatable, dimension(:,:), target :: TMP_PH1
   double precision, allocatable, dimension(:,:), target :: TMP_PH2
   double precision, allocatable, dimension(:,:), target :: TMP_PP1
   double precision, allocatable, dimension(:,:), target :: TMP_PP2
   double precision, allocatable, dimension(:,:), target :: TMP_PP3
   double precision, allocatable, dimension(:,:), target :: TMP_PHP

contains

subroutine InitDirect

   implicit none
   integer :: p, q, pq

   if (Init) then
      write(UOUT,*) 'InitDirect Error'
      call flush(UOUT)
      call abort
   endif
   Init = .TRUE.

   if (.not.allocated(TMP_NPSQ)) allocate(TMP_NPSQ(NPMOSQ))
   if (.not.allocated(TMP_HH))   allocate(TMP_HH(NHMO,NHMO))
   if (.not.allocated(TMP_PH1))  allocate(TMP_PH1(NPMO,NHMO))
   if (.not.allocated(TMP_PP1))  allocate(TMP_PP1(NPMO,NPMO))
   if (.not.allocated(TMP_PH2))  allocate(TMP_PH2(NPMO,NHMO))
   if (.not.allocated(TMP_PP2))  allocate(TMP_PP2(NPMO,NPMO))
   if (.not.allocated(TMP_PHP))  allocate(TMP_PHP(NHMO*NPMO,NPMO))
   if (.not.allocated(TMP_PP3))  allocate(TMP_PP3(NPMO,NPMO))
   if (.not.allocated(ResP))     allocate(ResP(MAX(NPMOSQ,NHMOSQ)))
   if (.not.allocated(ResQ))     allocate(ResQ(MAX(NPMOSQ,NHMOSQ)))

   if (get_iprint() > 5 ) then
     write(UOUT,*) 'InitDirect: allocated stuff '
   endif


   pq = 0
   do p = 2, MAX(NPMO,NHMO)
   do q = 1, (p-1)
      pq = pq + 1
      ResP(pq) = p
      ResQ(pq) = q
   enddo
   enddo

   return

end subroutine InitDirect

subroutine EndDirect

   implicit none

   if (.not.Init) then
      write(UOUT,*) 'EndDirect Error'
      call flush(UOUT)
      call abort
   endif
   Init = .FALSE.

   if (allocated(TMP_NPSQ)) deallocate(TMP_NPSQ)
   if (allocated(TMP_HH))   deallocate(TMP_HH)
   if (allocated(TMP_PH1))  deallocate(TMP_PH1)
   if (allocated(TMP_PP1))  deallocate(TMP_PP1)
   if (allocated(TMP_PH2))  deallocate(TMP_PH2)
   if (allocated(TMP_PP2))  deallocate(TMP_PP2)
   if (allocated(TMP_PHP))  deallocate(TMP_PHP)
   if (allocated(TMP_PP3))  deallocate(TMP_PP3)
   if (allocated(ResP))     deallocate(ResP)
   if (allocated(ResQ))     deallocate(ResQ)

   return

end subroutine EndDirect

subroutine DrvDirect(vac)

   use TimingMod
   implicit none
   integer, intent(in) :: vac

   if (.not.Init) then
      write(UOUT,*) 'Direct Diagram module not Initialized'
      call flush(UOUT)
      call Quit('DrvDirect: Direct Diagram module not Initialized!')
      !call abort
   endif

   if (get_iprint() >= 9) then
     write(UOUT,'(2x,a,i3)') 'DrvDirect: beginning for vac=',vac
   endif

   call Dir_V2H2P(TimeTaken(:,V2H2P),vac)
   call Dir_VPPPP(TimeTaken(:,VPPPP),vac)
   call Dir_VHHHH(TimeTaken(:,VHHHH),vac)
   call Dir_VHHPP(TimeTaken(:,VHHPP),vac)
   call Dir_V2HPP(TimeTaken(:,V2HPP),vac)
   call Dir_V2PHH(TimeTaken(:,V2PHH),vac)
   call Dir_VPHHP(TimeTaken(:,VPHHP),vac)
   call Dir_VPH2P(TimeTaken(:,VPH2P),vac)
   call Dir_VPHPP(TimeTaken(:,VPHPP),vac)
   call Dir_V2HHP(TimeTaken(:,V2HHP),vac)
   call Dir_VHHHP(TimeTaken(:,VHHHP),vac)
   call Dir_FHH(TimeTaken(:,FHH),vac)
   call Dir_FHP(TimeTaken(:,FHP),vac)
   call Dir_FPP(TimeTaken(:,FPP),vac)

   if (get_iprint() >= 9) then
     write(UOUT,'(2x,a,i3)') '... DrvDirect: end for vac=',vac
   endif

   return

end subroutine DrvDirect

#include "maxmin.h" 
!include "tind.h"

subroutine Dir_V2H2P(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd, SaveEndBufInd
   integer :: i, j, a, b, k, l, c, d
   integer :: ix, jx, lx, kx
   integer :: ij, ab, ik, ki, lj
   integer :: kj, ac, ca, cb, db
   integer :: eq, ti
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers
   call ReadFile(V2H2P)

   ! ... TMP1(c,k) = sum(l,d)[<kl||cd>*t(l->d)] ...

   TMP1 => TMP_PH1
   EndBufInd = 0
   TMP1 = dble(0)
   do k = 1, NHMO
      kx = (k-1)*NPMO
      do c = 1, NPMO
         StartBufInd = EndBufInd + 1
         EndBufInd   = EndBufInd + NZIntRec(kx+c)
         if (EndBufInd < StartBufInd) cycle
         VINT => IntBuf(StartBufInd:EndBufInd)
         VIND => IndBuf(StartBufInd:EndBufInd)
         TMP1(c,k) = DOT_PRODUCT(VINT,TOLD(VIND))
      enddo
   enddo

   ! ... sum diagrams using TMP1 ...

   TMP2 => TMP_HH
   TOLDP => TOLD(1:NSING)
   TNEWP => TNEW(1:NSING)

   call DGEMM('T','N',NHMO,NHMO,NPMO,1.0d00,TOLDP,  &
              NPMO,TMP1,NPMO,0.0d00,TMP2,NHMO)

   ! 01. - <kl||cd>*t(l->d)*t(i->c)*t(k->a)

   call DGEMM('N','T',NPMO,NHMO,NHMO,-1.0d00,TOLDP, &
              NPMO,TMP2,NHMO,1.0d00,TNEWP,NPMO)

   ! 02. - <kl||cd>*t(l->d)*t(i->c)*t(k,j->a,b) + P(i|j)

   do j = 1, NHMO
      ! k < j case
      do k = 1, (j-1)
         kj = MAXMIN(k,j)
         ti = TIND(kj)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         do i = 1, (j-1)
            ij = MAXMIN(i,j)
            eq = TIND(ij)
            TNEWP => TNEW(eq+1:eq+NPMOSQ)
            TNEWP = TNEWP - TMP2(i,k)*TOLDP
         enddo
         do i = (j+1), NHMO
            ij = MAXMIN(i,j)
            eq = TIND(ij)
            TNEWP => TNEW(eq+1:eq+NPMOSQ)
            TNEWP = TNEWP + TMP2(i,k)*TOLDP
         enddo
      enddo
      ! k > j case
      do k = (j+1), NHMO
         kj = MAXMIN(k,j)
         ti = TIND(kj)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         do i = 1, (j-1)
            ij = MAXMIN(i,j)
            eq = TIND(ij)
            TNEWP => TNEW(eq+1:eq+NPMOSQ)
            TNEWP = TNEWP + TMP2(i,k)*TOLDP
         enddo
         do i = (j+1), NHMO
            ij = MAXMIN(i,j)
            eq = TIND(ij)
            TNEWP => TNEW(eq+1:eq+NPMOSQ)
            TNEWP = TNEWP - TMP2(i,k)*TOLDP
         enddo
      enddo
   enddo

   ! 03. + <kl||cd>*t(l->d)*t(k,i->c,a)

   TMP2 => TMP_PHP
   TMP3 => TMP_PP3

   do i = 1, NHMO
      ix = (i-1)*NPMO
      ! TMP2(kc,a) = t(k,i->c,a)
      ! k < i case
      do k = 1, (i-1)
         ki = MAXMIN(k,i)
         ti = TIND(ki)
         kx = (k-1)*NPMO
         ca = 0
         do c = 1, NPMO
            TMP3(c,c) = dble(0)
            do a = 1, (c-1)
               ca = ca + 1
               TMP3(c,a) = -TOLD(ti+ca)
               TMP3(a,c) = +TOLD(ti+ca)
            enddo
         enddo
         TMP2(kx+1:kx+NPMO,:) = TMP3
      enddo
      ! k = i case
      TMP2(ix+1:ix+NPMO,:) = dble(0)
      ! k > i case
      do k = (i+1), NHMO
         ki = MAXMIN(k,i)
         ti = TIND(ki)
         kx = (k-1)*NPMO
         ca = 0
         do c = 1, NPMO
            TMP3(c,c) = dble(0)
            do a = 1, (c-1)
               ca = ca + 1
               TMP3(c,a) = +TOLD(ti+ca)
               TMP3(a,c) = -TOLD(ti+ca)
            enddo
         enddo
         TMP2(kx+1:kx+NPMO,:) = TMP3
      enddo
      TNEWP => TNEW(ix+1:ix+NPMO)
      call DGEMV('T',SIZE(TMP2,1),SIZE(TMP2,2),1.0d00, &
                 TMP2,SIZE(TMP2,1),TMP1,1,1.0d00,TNEWP,1)
   enddo

   ! 04. - <kl||cd>*t(l->d)*t(k->a)*t(i,j->c,b) + P(a|b)

   TMP2 => TMP_PP1
   TMP3 => TMP_PH2
   TOLDP => TOLD(1:NSING)

   do ij = 1, NHMOSQ
      ti = TIND(ij)
      cb = 0
      do c = 1, NPMO
         TMP2(c,c) = dble(0)
         do b = 1, (c-1)
            cb = cb + 1
            TMP2(b,c) = +TOLD(ti+cb)
            TMP2(c,b) = -TOLD(ti+cb)
         enddo
      enddo
      call DGEMM('N','N',NPMO,NHMO,NPMO,1.0d00,TMP2,NPMO,  &
                 TMP1,NPMO,0.0d00,TMP3,NPMO)
      call DGEMM('N','T',NPMO,NPMO,NHMO,1.0d00,TOLDP,NPMO, &
                 TMP3,NPMO,0.0d00,TMP2,NPMO)
      ab = 0
      do a = 2, NPMO
      do b = 1, (a-1)
         ab = ab + 1
         eq = ti + ab
         TNEW(eq) = TNEW(eq) - TMP2(a,b) + TMP2(b,a)
      enddo
      enddo
   enddo

   ! 05. - <kl||cd>*t(l,j->d,b)*t(k->a)*t(i->c) + P(i|j)P(a|b)
   ! 06. + <kl||cd>*t(l,j->d,b)*t(i,k->a,c) + P(a|b) or P(i|j)

   TMP1 => TMP_PHP
   TMP2 => TMP_PP3
   TMP3 => TMP_PH2
   TMP4 => TMP_PP1
   TMP5 => TMP_PP2

   do j = 1, NHMO
      jx = (j-1)*NPMO
      ! Rearrangement
      ! TMP1(ld,b) = t(l,j->d,b)
      ! l < j case
      do l = 1, (j-1)
         lj = MAXMIN(l,j)
         ti = TIND(lj)
         lx = (l-1)*NPMO
         db = 0
         do d = 1, NPMO
            TMP2(d,d) = dble(0)
            do b = 1, (d-1)
               db = db + 1
               TMP2(d,b) = -TOLD(ti+db)
               TMP2(b,d) = +TOLD(ti+db)
            enddo
         enddo
         TMP1(lx+1:lx+NPMO,:) = TMP2
      enddo
      ! l = j case
      TMP1(jx+1:jx+NPMO,:) = dble(0)
      ! l > j case
      do l = (j+1), NHMO
         lj = MAXMIN(l,j)
         ti = TIND(lj)
         lx = (l-1)*NPMO
         db = 0
         do d = 1, NPMO
            TMP2(d,d) = dble(0)
            do b = 1, (d-1)
               db = db + 1
               TMP2(d,b) = +TOLD(ti+db)
               TMP2(b,d) = -TOLD(ti+db)
            enddo
         enddo
         TMP1(lx+1:lx+NPMO,:) = TMP2
      enddo

      ! Diagrams 05 & 06

      EndBufInd = 0
      do k = 1, NHMO
         kx = (k-1)*NPMO
         SaveEndBufInd = EndBufInd
         ! TMP2(c,b) = sum(l,d)[<kl||cd>*t(l,j->d,b)]
         TMP2 = dble(0)
         do b = 1, NPMO
            EndBufInd = SaveEndBufInd
            do c = 1, NPMO
               StartBufInd = EndBufInd + 1
               EndBufInd   = EndBufInd + NZIntRec(kx+c)
               if (EndBufInd < StartBufInd) cycle
               VINT => IntBuf(StartBufInd:EndBufInd)
               VIND => IndBuf(StartBufInd:EndBufInd)
               TMP2(c,b) = DOT_PRODUCT(VINT,TMP1(VIND,b))
            enddo
         enddo

         ! 05. - <kl||cd>*t(l,j->d,b)*t(k->a)*t(i->c) + P(i|j)P(a|b)
         ! TMP3(b,i) = sum(l,c,d)[<kl||cd>*t(l,j->d,b)*t(i->c)]

         TOLDP => TOLD(1:NSING)
         call DGEMM('T','N',NPMO,NHMO,NPMO,1.0d00,TMP2,NPMO, &
                         TOLDP,NPMO,0.0d00,TMP3,NPMO)
         do i = 1, (j-1)
            ij = MAXMIN(i,j)
            ti = TIND(ij)
            ab = 0
            do a = 2, NPMO
            do b = 1, (a-1)
               ab = ab + 1
               eq = ti + ab
               TNEW(eq) = TNEW(eq) + TMP3(b,i)*TOLD(kx+a) - TMP3(a,i)*TOLD(kx+b)
            enddo
            enddo
         enddo
         do i = (j+1), NHMO
            ij = MAXMIN(i,j)
            ti = TIND(ij)
            ab = 0
            do a = 2, NPMO
            do b = 1, (a-1)
               ab = ab + 1
               eq = ti + ab
               TNEW(eq) = TNEW(eq) - TMP3(b,i)*TOLD(kx+a) + TMP3(a,i)*TOLD(kx+b)
            enddo
            enddo
         enddo

         ! 06. + <kl||cd>*t(l,j->d,b)*t(i,k->a,c) + P(a|b) or P(i|j)
         !       [ here P(a|b) is being used => i>j ]

         do i = (j+1), NHMO
            if (i==k) cycle
            ij = MAXMIN(i,j)
            eq = TIND(ij)
            ik = MAXMIN(i,k)
            ti = TIND(ik)
            ! TMP4(c,a) = t(i,k->a,c)
            ac = 0
            do a = 1, NPMO
               TMP4(a,a) = dble(0)
               do c = 1, (a-1)
                  ac = ac + 1
                  TMP4(c,a) = +TOLD(ti+ac)
                  TMP4(a,c) = -TOLD(ti+ac)
               enddo
            enddo
            ! TMP5(b,a) = sum(l,c,d) [ <kl||cd>*t(l,j->d,b)*t(i,k->a,c) ]
            call DGEMM('T','N',NPMO,NPMO,NPMO,1.0d00,TMP2,NPMO, &
                        TMP4,NPMO,0.0d00,TMP5,NPMO)
            if (i<k) then
               ab = 0
               do a = 2, NPMO
               do b = 1, (a-1)
                  ab = ab + 1
                  TNEW(eq+ab) = TNEW(eq+ab) - TMP5(b,a) + TMP5(a,b)
               enddo
               enddo
            else
               ab = 0
               do a = 2, NPMO
               do b = 1, (a-1)
                  ab = ab + 1
                  TNEW(eq+ab) = TNEW(eq+ab) + TMP5(b,a) - TMP5(a,b)
               enddo
               enddo
            endif
         enddo
      enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Dir_V2H2P

subroutine Dir_VPPPP(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd
   integer :: Block, St, En
   integer :: i, j
   integer :: ix, jx, ij, ab
   integer :: eq, ti
   integer, dimension(:), pointer :: car, dar
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers

   ! + sum(c>d) <ab||cd>*[t(i,j->c,d)+t(i->c)*t(j->d)-t(i->d)*t(j->c)]

   TMP => TMP_NPSQ
   car => ResP
   dar => ResQ
   do Block = 1, NumBlock(VPPPP)
      call ReadFile(VPPPP,block,St,En)
      ij = 0
      do i = 2, NHMO
      do j = 1, (i-1)
         ij = ij + 1
         ti = TIND(ij)
         ix = (i-1)*NPMO
         jx = (j-1)*NPMO
         TMP = TOLD(ti+1:ti+NPMOSQ) + TOLD(ix+car)*TOLD(jx+dar) &
                                    - TOLD(ix+dar)*TOLD(jx+car)
         EndBufInd = 0
         do ab = St, En
            StartBufInd = EndBufInd + 1
            EndBufInd   = EndBufInd + NZIntRec(ab)
            if ( EndBufInd < StartBufInd ) cycle
            VINT => IntBuf(StartBufInd:EndBufInd)
            VIND => IndBuf(StartBufInd:EndBufInd)
            eq = ti + ab
            TNEW(eq) = TNEW(eq) + DOT_PRODUCT(VINT,TMP(VIND))
         enddo
      enddo
      enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return
       
end subroutine Dir_VPPPP

subroutine Dir_VHHHH(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd
   integer :: k, l, ind
   integer :: kx, lx, ij, kl
   integer :: eq, ti
   integer, dimension(:), pointer :: aar, bar
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers

   ! + sum(k>l) <kl||ij>*[t(k,l->a,b)+t(k->a)*t(l->b)-t(l->a)*t(k->b)]

   TMP => TMP_NPSQ
   aar => ResP
   bar => ResQ
   call ReadFile(VHHHH)
   EndBufInd = 0
   kl = 0
   do k = 2, NHMO
   do l = 1, (k-1)
      kl = kl + 1
      StartBufInd = EndBufInd + 1
      EndBufInd   = EndBufInd + NZIntRec(kl)
      if ( EndBufInd < StartBufInd ) cycle
      VINT => IntBuf(StartBufInd:EndBufInd)
      VIND => IndBuf(StartBufInd:EndBufInd)
      ti = TIND(kl)
      kx = (k-1)*NPMO
      lx = (l-1)*NPMO
      TMP = TOLD(ti+1:ti+NPMOSQ) + TOLD(kx+aar)*TOLD(lx+bar) &
                                 - TOLD(kx+bar)*TOLD(lx+aar)
      do ind = 1, SIZE(VIND)
         ij = VIND(ind)
         ti = TIND(ij)
         TNEWP => TNEW(ti+1:ti+NPMOSQ)
         TNEWP = TNEWP + VINT(ind)*TMP
      enddo
   enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Dir_VHHHH

subroutine Dir_VHHPP(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd
   integer :: i, j, k, l
   integer :: kl, ij, ix, jx, kx, lx
   integer :: ti, eq
   integer, dimension(:), pointer :: car, dar
   double precision, dimension(NHMOSQ,NHMOSQ) :: THSQ
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)

   call NullifyPointers
   call ReadFile(VHHPP)

   ! 01. + <i,j||a,b>

   EndBufInd = 0
   do ij = 1, NHMOSQ
      StartBufInd = EndBufInd + 1
      EndBufInd   = EndBufInd + NZIntRec(ij)
      if ( EndBufInd < StartBufInd ) cycle
      VINT => IntBuf(StartBufInd:EndBufInd)
      VIND => IndBuf(StartBufInd:EndBufInd)
      ti = TIND(ij)
      TNEW(ti+VIND) = TNEW(ti+VIND) + VINT
   enddo

   ! THSQ(ij,kl) = sum(c>d) [<kl||cd>*(t(i,j->c,d)+t(i->c)*t(j->d)-t(i->d)*t(j->c))]

   TMP => TMP_NPSQ
   car => ResP
   dar => ResQ
   ij = 0
   THSQ = dble(0)
   do i = 2, NHMO
   do j = 1, (i-1)
      ij = ij + 1
      ti = TIND(ij)
      ix = (i-1)*NPMO
      jx = (j-1)*NPMO
      TMP = TOLD(ti+1:ti+NPMOSQ) + TOLD(ix+car)*TOLD(jx+dar) &
                                 - TOLD(ix+dar)*TOLD(jx+car)
      EndBufInd = 0
      kl = 0
      do k = 2, NHMO
      do l = 1, (k-1)
         kl = kl + 1
         StartBufInd = EndBufInd + 1
         EndBufInd   = EndBufInd + NZIntRec(kl)
         if ( EndBufInd < StartBufInd ) cycle
         VINT => IntBuf(StartBufInd:EndBufInd)
         VIND => IndBuf(StartBufInd:EndBufInd)
         THSQ(kl,ij) = DOT_PRODUCT(VINT,TMP(VIND))
      enddo
      enddo
   enddo
   enddo
   THSQ = TRANSPOSE(THSQ)

   ! 02. + sum(k>l,c>d) [ <k,l||c,d>*T(i,j->c,d)*T(k,l->a,b) ]
   !       T(i,j->c,d) = t(i,j->c,d) + t(i->c)*t(j->d) - t(i->d)*t(j->c)
   !       T(k,l->a,b) = t(k,l->a,b) + t(k->a)*t(l->b) - t(l->a)*t(k->b)
   kl = 0
   do k = 2, NHMO
   do l = 1, (k-1)
      kl = kl + 1
      ti = TIND(kl)
      kx = (k-1)*NPMO
      lx = (l-1)*NPMO
      TMP = TOLD(ti+1:ti+NPMOSQ) + TOLD(kx+car)*TOLD(lx+dar) &
                                 - TOLD(kx+dar)*TOLD(lx+car)
      do ij = 1, NHMOSQ
         eq = TIND(ij)
         TNEWP => TNEW(eq+1:eq+NPMOSQ)
         TNEWP = TNEWP + THSQ(ij,kl)*TMP
      enddo
   enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Dir_VHHPP

subroutine Dir_V2HPP(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd
   integer :: k, l, i, j
   integer :: kx, il, ij, kj, ki
   integer :: ti, eq
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)

   call NullifyPointers
   call ReadFile(V2HPP)

   ! TMP1(k,i) = sum(l,c>d)[<kl||cd>*t(i,l->c,d)]

   TMP1 => TMP_HH
   TMP1 = dble(0)
   EndBufInd = 0
   do k = 1, NHMO
      kx = (k-1)*NHMO
      do l = 1, NHMO
         StartBufInd = EndBufInd + 1
         EndBufInd   = EndBufInd + NZIntRec(kx+l)
         if ( EndBufInd < StartBufInd ) cycle
         VINT => IntBuf(StartBufInd:EndBufInd)
         VIND => IndBuf(StartBufInd:EndBufInd)
         do i = 1, (l-1)
            il = MAXMIN(i,l)
            ti = TIND(il)
            TMP1(i,k) = TMP1(i,k) - DOT_PRODUCT(VINT,TOLD(ti+VIND))
         enddo
         do i = (l+1), NHMO
            il = MAXMIN(i,l)
            ti = TIND(il)
            TMP1(i,k) = TMP1(i,k) + DOT_PRODUCT(VINT,TOLD(ti+VIND))
         enddo
      enddo
   enddo
   TMP1 = TRANSPOSE(TMP1)

   ! 01. - sum(k,l,cd) [ <kl||cd>*t(i,l->c,d)*t(k->a) ]

   TOLDP => TOLD(1:NSING)
   TNEWP => TNEW(1:NSING)
   call DGEMM('N','N',NPMO,NHMO,NHMO,-1.0d00,TOLDP, &
              NPMO,TMP1,NHMO,1.0d00,TNEWP,NPMO)

   ! 02. - sum(k,l,cd) [ <kl||cd>*t(i,l->c,d)*t(k,j->a,b) ] + P(i,j)

   ij = 0
   do i = 2, NHMO
   do j = 1, (i-1)
      ij = ij + 1
      eq = TIND(ij)
      TNEWP => TNEW(eq+1:eq+NPMOSQ)
      ! i>j, k<j
      do k = 1, (j-1)
         kj = MAXMIN(k,j)
         ti = TIND(kj)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         TNEWP = TNEWP + TMP1(k,i)*TOLDP
      enddo
      ! i<j, k<j
      do k = 1, (i-1)
         ki = MAXMIN(k,i)
         ti = TIND(ki)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         TNEWP = TNEWP - TMP1(k,j)*TOLDP
      enddo
      ! i>j, k>j
      do k = (j+1), NHMO
         kj = MAXMIN(k,j)
         ti = TIND(kj)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         TNEWP = TNEWP - TMP1(k,i)*TOLDP
      enddo
      ! i<j, k>j
      do k = (i+1), NHMO
         ki = MAXMIN(k,i)
         ti = TIND(ki)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         TNEWP = TNEWP + TMP1(k,j)*TOLDP
      enddo
   enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Dir_V2HPP

subroutine Dir_V2PHH(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd
   integer :: a, b, c, d
   integer :: kl, ij, cx, ad, ab, cb
   integer :: ti, eq
   double precision, dimension(NHMOSQ,NPMO) :: THSQP
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)

   call NullifyPointers
   call ReadFile(V2PHH)

   ! TMP1(c,a) = + sum(d,k>l)[<cd||kl>*t(k,l->a,d)]

   TMP1 => TMP_PP1
   TMP1 = dble(0)
   do a = 1, NPMO
      ! THSQP(kl,d) => t(k,l->a,d)
      do kl = 1, NHMOSQ
         ti = TIND(kl)
         do d = 1, (a-1)
            ad = MAXMIN(a,d)
            THSQP(kl,d) = +TOLD(ti+ad)
         enddo
         THSQP(kl,a) = dble(0)
         do d = (a+1), NPMO
            ad = MAXMIN(a,d)
            THSQP(kl,d) = -TOLD(ti+ad)
         enddo
      enddo
      EndBufInd = 0
      do c = 1, NPMO
         cx = (c-1)*NPMO
         do d = 1, NPMO
            StartBufInd = EndBufInd + 1
            EndBufInd   = EndBufInd + NZIntRec(cx+d)
            if ( EndBufInd < StartBufInd ) cycle
            VINT => IntBuf(StartBufInd:EndBufInd)
            VIND => IndBuf(StartBufInd:EndBufInd)
            TMP1(c,a) = TMP1(c,a) + DOT_PRODUCT(VINT,THSQP(VIND,d))
         enddo
      enddo
   enddo
            
   ! 01. - sum (k>l,c,d) [ <cd||kl>*t(k,l->a,d)*t(i->c) ]

   TOLDP => TOLD(1:NSING)
   TNEWP => TNEW(1:NSING)
   call DGEMM('T','N',NPMO,NHMO,NPMO,-1.0d00,TMP1, &
              NPMO,TOLDP,NPMO,1.0d00,TNEWP,NPMO)

   ! 02. - sum (k>l,c,d) [ <cd||kl>*t(k,l->a,d)*t(i,j->c,b) ] + P(a,b)

   TMP2 => TMP_PP2
   TMP3 => TMP_PP3
   do ij = 1, NHMOSQ
      ti = TIND(ij)
      cb = 0
      do c = 1, NPMO
         TMP2(c,c) = dble(0)
         do b = 1, (c-1)
            cb = cb + 1
            TMP2(c,b) = +TOLD(ti+cb)
            TMP2(b,c) = -TOLD(ti+cb)
         enddo
      enddo
      call DGEMM('T','N',NPMO,NPMO,NPMO,1.0d00,TMP1, &
                 NPMO,TMP2,NPMO,0.0d00,TMP3,NPMO)
      ab = 0
      do a = 2, NPMO
      do b = 1, (a-1)
         ab = ab + 1
         eq = ti + ab
         TNEW(eq) = TNEW(eq) - TMP3(a,b) + TMP3(b,a)
      enddo
      enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

end subroutine Dir_V2PHH

subroutine Dir_VPHHP(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd, SaveEndBufInd
   integer :: i, j, a, b, l, d
   integer :: ix, jx, lx, lj, db
   integer :: ij, ab
   integer :: ti, eq
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers
   call ReadFile(VPHHP)

   ! 01. + sum(k,c) [ <ak||ic>*t(k->c) ]

   EndBufInd = 0
   do i = 1, NHMO
      ix = (i-1)*NPMO
      do a = 1, NPMO
         StartBufInd = EndBufInd + 1
         EndBufInd   = EndBufInd + NZIntRec(ix+a)
         if ( EndBufInd < StartBufInd ) cycle
         VINT => IntBuf(StartBufInd:EndBufInd)
         VIND => IndBuf(StartBufInd:EndBufInd)
         TNEW(ix+a) = TNEW(ix+a) + DOT_PRODUCT(VINT,TOLD(VIND))
      enddo
   enddo

   ! 02. + sum(l,d) [ <al||id>*{t(l,j->d,b)-t(l->b)*t(j->d)} ]

   TMP1 => TMP_PHP
   TMP2 => TMP_PP1
   do j = 1, NHMO
      jx = (j-1)*NPMO
      ! TMP1(ld,b) = [ t(l,j->d,b) - t(l->b)*t(j->d) ]
      ! l < j case
      do l = 1, (j-1)
         lj = MAXMIN(l,j)
         ti = TIND(lj)
         lx = (l-1)*NPMO
         db = 0
         do d = 1, NPMO
            TMP2(d,d) = -TOLD(lx+d)*TOLD(jx+d)
            do b = 1, (d-1)
               db = db + 1
               TMP2(d,b) = -TOLD(ti+db) - TOLD(lx+b)*TOLD(jx+d)
               TMP2(b,d) = +TOLD(ti+db) - TOLD(lx+d)*TOLD(jx+b)
            enddo
         enddo
         TMP1(lx+1:lx+NPMO,:) = TMP2
      enddo
      ! l = j case
      do b = 1, NPMO
         TMP1(jx+1:jx+NPMO,b) = -TOLD(jx+b)*TOLD(jx+1:jx+NPMO)
      enddo
      ! l > j case
      do l = (j+1), NHMO
         lj = MAXMIN(l,j)
         ti = TIND(lj)
         lx = (l-1)*NPMO
         db = 0
         do d = 1, NPMO
            TMP2(d,d) = -TOLD(lx+d)*TOLD(jx+d)
            do b = 1, (d-1)
               db = db + 1
               TMP2(d,b) = +TOLD(ti+db) - TOLD(lx+b)*TOLD(jx+d)
               TMP2(b,d) = -TOLD(ti+db) - TOLD(lx+d)*TOLD(jx+b)
            enddo
         enddo
         TMP1(lx+1:lx+NPMO,:) = TMP2
      enddo

      EndBufInd = 0
      do i = 1, NHMO
         ix = (i-1)*NPMO
         if (i==j) EndBufInd = EndBufInd + SUM(NZIntRec(ix+1:ix+NPMO))
         SaveEndBufInd = EndBufInd
         if (i==j) cycle
         ij = MAXMIN(i,j)
         ti = TIND(ij)
         if (i>j) then
            do b = 1, NPMO
               EndBufInd = SaveEndBufInd
               do a = 1, NPMO
                  StartBufInd = EndBufInd + 1
                  EndBufInd   = EndBufInd + NZIntRec(ix+a)
                  if ( EndBufInd < StartBufInd ) cycle
                  if (a==b) cycle
                  VINT => IntBuf(StartBufInd:EndBufInd)
                  VIND => IndBuf(StartBufInd:EndBufInd)
                  ab = MAXMIN(a,b)
                  eq = ti + ab
                  if (a>b) then
                     TNEW(eq) = TNEW(eq) + DOT_PRODUCT(VINT,TMP1(VIND,b))
                  else
                     TNEW(eq) = TNEW(eq) - DOT_PRODUCT(VINT,TMP1(VIND,b))
                  endif
               enddo
            enddo
         else
            do b = 1, NPMO
               EndBufInd = SaveEndBufInd
               do a = 1, NPMO
                  StartBufInd = EndBufInd + 1
                  EndBufInd   = EndBufInd + NZIntRec(ix+a)
                  if ( EndBufInd < StartBufInd ) cycle
                  if (a==b) cycle
                  VINT => IntBuf(StartBufInd:EndBufInd)
                  VIND => IndBuf(StartBufInd:EndBufInd)
                  ab = MAXMIN(a,b)
                  eq = ti + ab
                  if (a>b) then
                     TNEW(eq) = TNEW(eq) - DOT_PRODUCT(VINT,TMP1(VIND,b))
                  else
                     TNEW(eq) = TNEW(eq) + DOT_PRODUCT(VINT,TMP1(VIND,b))
                  endif
               enddo
            enddo
         endif
      enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Dir_VPHHP

subroutine Dir_V2HHP(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd, SaveEndBufInd
   integer :: i, j, a, b, k, l, d
   integer :: kj, ki, kx, ij, ab
   integer :: lj, lx, jx, db
   integer :: ti, eq
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers
   call ReadFile(V2HHP)

   ! TMP1(k,i) = sum(l,d) [ <kl||id>*t(l->d) ]

   TMP1 => TMP_HH
   TMP1 = dble(0)
   EndBufInd = 0
   do k = 1, NHMO
      kx = (k-1)*NHMO
      do i = 1, NHMO
         StartBufInd = EndBufInd + 1
         EndBufInd = EndBufInd + NZIntRec(kx+i)
         if ( EndBufInd < StartBufInd ) cycle
         VINT => IntBuf(StartBufInd:EndBufInd)
         VIND => IndBuf(StartBufInd:EndBufInd)
         TMP1(i,k) = DOT_PRODUCT(VINT,TOLD(VIND))
      enddo
   enddo
   TMP1 = TRANSPOSE(TMP1)

   ! 01. - sum(k,l,d) [ <k,l||i,d>*t(l->d)*t(k->a) ]

   TOLDP => TOLD(1:NSING)
   TNEWP => TNEW(1:NSING)
   call DGEMM('N','N',NPMO,NHMO,NHMO,-1.0d00,TOLDP, &
              NPMO,TMP1,NHMO,1.0d00,TNEWP,NPMO)

   ! 02. - sum(k,l,d) [ <k,l||i,d>*t(l->d)*t(k,j->a,b) ] + P(i|j)

   ij = 0
   do i = 2, NHMO
   do j = 1, (i-1)
      ij = ij + 1
      eq = TIND(ij)
      TNEWP => TNEW(eq+1:eq+NPMOSQ)
      ! i>j, k<j
      do k = 1, (j-1)
         kj = MAXMIN(k,j)
         ti = TIND(kj)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         TNEWP = TNEWP + TMP1(k,i)*TOLDP
      enddo
      ! i<j, k<j
      do k = 1, (i-1)
         ki = MAXMIN(k,i)
         ti = TIND(ki)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         TNEWP = TNEWP - TMP1(k,j)*TOLDP
      enddo
      ! i>j, k>j
      do k = (j+1), NHMO
         kj = MAXMIN(k,j)
         ti = TIND(kj)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         TNEWP = TNEWP - TMP1(k,i)*TOLDP
      enddo
      ! i<j, k>j
      do k = (i+1), NHMO
         ki = MAXMIN(k,i)
         ti = TIND(ki)
         TOLDP => TOLD(ti+1:ti+NPMOSQ)
         TNEWP = TNEWP + TMP1(k,j)*TOLDP
      enddo
   enddo
   enddo

   ! - sum(k,l,d) [ <k,l||i,d>*{t(l,j->d,b)-0.5*t(l->b)*t(j->d)}*t(k->a) ] + P(i|j)P(a,b)

   TMP1 => TMP_PHP
   TMP2 => TMP_PP1
   TMP3 => TMP_PH2
   do j = 1, NHMO
      jx = (j-1)*NPMO
      ! TMP1(ld,b) = [ t(l,j->d,b) - 0.5*t(l->b)*t(j->d) ]
      ! l < j case
      do l = 1, (j-1)
         lj = MAXMIN(l,j)
         ti = TIND(lj)
         lx = (l-1)*NPMO
         db = 0
         do d = 1, NPMO
            TMP2(d,d) = -0.5D00*TOLD(lx+d)*TOLD(jx+d)
            do b = 1, (d-1)
               db = db + 1
               TMP2(d,b) = -TOLD(ti+db) - 0.5D00*TOLD(lx+b)*TOLD(jx+d)
               TMP2(b,d) = +TOLD(ti+db) - 0.5D00*TOLD(lx+d)*TOLD(jx+b)
            enddo
         enddo
         TMP1(lx+1:lx+NPMO,:) = TMP2
      enddo
      ! l = j case
      do b = 1, NPMO
         TMP1(jx+1:jx+NPMO,b) = -0.5D00*TOLD(jx+b)*TOLD(jx+1:jx+NPMO)
      enddo
      ! l > j case
      do l = (j+1), NHMO
         lj = MAXMIN(l,j)
         ti = TIND(lj)
         lx = (l-1)*NPMO
         db = 0
         do d = 1, NPMO
            TMP2(d,d) = -0.5D00*TOLD(lx+d)*TOLD(jx+d)
            do b = 1, (d-1)
               db = db + 1
               TMP2(d,b) = +TOLD(ti+db) - 0.5D00*TOLD(lx+b)*TOLD(jx+d)
               TMP2(b,d) = -TOLD(ti+db) - 0.5D00*TOLD(lx+d)*TOLD(jx+b)
            enddo
         enddo
         TMP1(lx+1:lx+NPMO,:) = TMP2
      enddo

      EndBufInd = 0
      do k = 1, NHMO
         kx = (k-1)*NHMO
         SaveEndBufInd = EndBufInd
         ! TMP3(b,i) = sum(l,d) [ <k,l||i,d>*{t(l,j->d,b)-0.5*t(l->b)*t(j->d)}*t(k->a) ]
         TMP3 = dble(0)
         do b = 1, NPMO
            EndBufInd = SaveEndBufInd
            do i = 1, NHMO
               StartBufInd = EndBufInd + 1
               EndBufInd   = EndBufInd + NZIntRec(kx+i)
               if ( EndBufInd < StartBufInd ) cycle
               VINT => IntBuf(StartBufInd:EndBufInd)
               VIND => IndBuf(StartBufInd:EndBufInd)
               TMP3(b,i) = DOT_PRODUCT(VINT,TMP1(VIND,b))
            enddo
         enddo
         kx = (k-1)*NPMO
         do i = 1, (j-1)
            ij = MAXMIN(i,j)
            ti = TIND(ij)
            ab = 0
            do a = 2, NPMO
            do b = 1, (a-1)
               ab = ab + 1
               eq = ti + ab
               TNEW(eq) = TNEW(eq) + TMP3(b,i)*TOLD(kx+a) &
                                   - TMP3(a,i)*TOLD(kx+b)
            enddo
            enddo
         enddo
         do i = (j+1), NHMO
            ij = MAXMIN(i,j)
            ti = TIND(ij)
            ab = 0
            do a = 2, NPMO
            do b = 1, (a-1)
               ab = ab + 1
               eq = ti + ab
               TNEW(eq) = TNEW(eq) - TMP3(b,i)*TOLD(kx+a) &
                                   + TMP3(a,i)*TOLD(kx+b)
            enddo
            enddo
         enddo
      enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return
       
end subroutine Dir_V2HHP

subroutine Dir_VHHHP(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd
   integer :: i, j, a, d
   integer :: kl, ad, id, ix, ij
   integer :: ti, eq
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers
   call ReadFile(VHHHP)

   ! 01. - sum(k>l,d) [ <k,l||i,d>*T(k,l->a,d) ]
   ! 02. - sum(k) [ <i,j||k,b>*t(k->a) ] + P(a,b) => t(i,j->a,b)
   !     - sum(i) [ <k,l||i,d>*t(i->a) ] + P(a,d) => t(k,l->a,d)
   ! 03. + sum(k>l,d) [<k,l||i,d>*t(j->d)*t(k,l->a,b)] + P(i,j)

   TMP1 => TMP_PP1
   TMP2 => TMP_HH
   EndBufInd = 0
   do kl = 1, NHMOSQ
      StartBufInd = EndBufInd + 1
      EndBufInd   = EndBufInd + NZIntRec(kl)
      if ( EndBufInd < StartBufInd ) cycle
      VINT => IntBuf(StartBufInd:EndBufInd)
      VIND => IndBuf(StartBufInd:EndBufInd)
      ti = TIND(kl)
      TOLDP => TOLD(ti+1:ti+NPMOSQ)
      TMP2 = dble(0)
      ! TMP1(a,d) = t(k,l->a,d)
      ad = 0
      do a = 1, NPMO
         TMP1(a,a) = dble(0)
         do d = 1, (a-1)
            ad = ad + 1
            TMP1(a,d) = +TOLD(ti+ad)
            TMP1(d,a) = -TOLD(ti+ad)
         enddo
      enddo
      do id = 1, SIZE(VIND)
         i = (VIND(id)-1)/NPMO + 1
         d = VIND(id) - (i-1)*NPMO
         ix = (i-1)*NPMO
         ! 01. - sum(k>l,d) [ <k,l||i,d>*T(k,l->a,d) ]
         TNEWP => TNEW(ix+1:ix+NPMO)
         TNEWP = TNEWP - VINT(id)*TMP1(:,d)
         ! 02. - sum(k) [ <i,j||k,b>*t(k->a) ] + P(a,b) => t(i,j->a,b)
         !     - sum(i) [ <k,l||i,d>*t(i->a) ] + P(a,d) => t(k,l->a,d)
         do a = 1, (d-1)
            ad = MAXMIN(a,d)
            eq = ti + ad
            TNEW(eq) = TNEW(eq) + VINT(id)*TOLD(ix+a)
         enddo
         do a = (d+1), NPMO
            ad = MAXMIN(a,d)
            eq = ti + ad
            TNEW(eq) = TNEW(eq) - VINT(id)*TOLD(ix+a)
         enddo
         ! TMP2(j,i) = sum(d) [<k,l||i,d>*t(j->d)]
         TMP2(:,i) = TMP2(:,i) + VINT(id)*TOLD(d:NSING:NPMO)
      enddo
      ! 03. + sum(k>l,d) [<k,l||i,d>*t(j->d)*t(k,l->a,b)] + P(i,j)
      ij = 0
      do i = 2, NHMO
      do j = 1, (i-1)
         ij = ij + 1
         eq = TIND(ij)
         TNEWP => TNEW(eq+1:eq+NPMOSQ)
         TNEWP = TNEWP + (TMP2(j,i)-TMP2(i,j))*TOLDP
      enddo
      enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Dir_VHHHP

subroutine Dir_VPH2P(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd, SaveEndBufInd
   integer :: Block, St, En
   integer :: i, j, l, a, b, c, d
   integer :: ax, jx, lj, lx
   integer :: db, ij, ab, cb
   integer :: ti, eq
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers

   TMP1 => TMP_PP1
   TMP2 => TMP_PP2
   TMP3 => TMP_PHP
   TMP4 => TMP_PH2

   TMP1 = dble(0)
   do Block = 1, NumBlock(VPH2P)
      call ReadFile(VPH2P,Block,St,En)
      EndBufInd = 0
      do a = St, En
         ax = (a-1)*NPMO
         do c = 1, NPMO
            StartBufInd = EndBufInd + 1
            EndBufInd   = EndBufInd + NZIntRec(ax+c)
            if ( EndBufInd < StartBufInd ) cycle
            VINT => IntBuf(StartBufInd:EndBufInd)
            VIND => IndBuf(StartBufInd:EndBufInd)
            TMP1(c,a) = DOT_PRODUCT(VINT,TOLD(VIND))
         enddo
      enddo
      do j = 1, NHMO
         jx = (j-1)*NPMO
         ! TMP3(ld,b) = t(l,j->d,b) - 0.5*t(l->b)*t(j->d)
         ! l < j case
         do l = 1, (j-1)
            lj = MAXMIN(l,j)
            ti = TIND(lj)
            lx = (l-1)*NPMO
            db = 0
            do d = 1, NPMO
               TMP2(d,d) = - 0.5D00*TOLD(lx+d)*TOLD(jx+d)
               do b = 1, (d-1)
                  db = db + 1
                  TMP2(d,b) = -TOLD(ti+db) - 0.5D00*TOLD(lx+b)*TOLD(jx+d)
                  TMP2(b,d) = +TOLD(ti+db) - 0.5D00*TOLD(lx+d)*TOLD(jx+b)
               enddo
            enddo
            TMP3(lx+1:lx+NPMO,:) = TMP2
         enddo
         ! l = j case
         do b = 1, NPMO
            TMP3(jx+1:jx+NPMO,b) = - 0.5D00*TOLD(jx+b)*TOLD(jx+1:jx+NPMO)
         enddo
         ! l > j case
         do l = (j+1), NHMO
            lj = MAXMIN(l,j)
            ti = TIND(lj)
            lx = (l-1)*NPMO
            db = 0
            do d = 1, NPMO
               TMP2(d,d) = - 0.5D00*TOLD(lx+d)*TOLD(jx+d)
               do b = 1, (d-1)
                  db = db + 1
                  TMP2(d,b) = +TOLD(ti+db) - 0.5D00*TOLD(lx+b)*TOLD(jx+d)
                  TMP2(b,d) = -TOLD(ti+db) - 0.5D00*TOLD(lx+d)*TOLD(jx+b)
               enddo
            enddo
            TMP3(lx+1:lx+NPMO,:) = TMP2
         enddo
         EndBufInd = 0
         do a = St, En
            ax = (a-1)*NPMO
            SaveEndBufInd = EndBufInd
            ! TMP2(c,b) = sum(l,d) [ <al||cd>*(t(l,j->d,b) - 0.5*t(l->b)*t(j->d)) ]
            TMP2 = dble(0)
            do b = 1, NPMO
               EndBufInd = SaveEndBufInd
               do c = 1, NPMO
                  StartBufInd = EndBufInd + 1
                  EndBufInd = EndBufInd + NZIntRec(ax+c)
                  if ( EndBufInd < StartBufInd ) cycle
                  VINT => IntBuf(StartBufInd:EndBufInd)
                  VIND => IndBuf(StartBufInd:EndBufInd)
                  TMP2(c,b) = DOT_PRODUCT(VINT,TMP3(VIND,b))
               enddo
            enddo

            ! 03. + sum(l,c,d) [ <al||cd>*(t(l,j->d,b)-0.5*t(l->b)*t(j->d))*t(i->c) ] + P(i|j)P(a|b)

            TOLDP => TOLD(1:NSING)
            call DGEMM('T','N',NPMO,NHMO,NPMO,1.0d00,TMP2, &
                       NPMO,TOLDP,NPMO,0.0d00,TMP4,NPMO)

            do i = 1, (j-1)
               ij = MAXMIN(i,j)
               ti = TIND(ij)
               do b = 1, (a-1)
                  ab = MAXMIN(a,b)
                  eq = ti + ab
                  TNEW(eq) = TNEW(eq) - TMP4(b,i)
               enddo
               do b = (a+1), NPMO
                  ab = MAXMIN(a,b)
                  eq = ti + ab
                  TNEW(eq) = TNEW(eq) + TMP4(b,i)
               enddo
            enddo
            do i = (j+1), NHMO
               ij = MAXMIN(i,j)
               ti = TIND(ij)
               do b = 1, (a-1)
                  ab = MAXMIN(a,b)
                  eq = ti + ab
                  TNEW(eq) = TNEW(eq) + TMP4(b,i)
               enddo
               do b = (a+1), NPMO
                  ab = MAXMIN(a,b)
                  eq = ti + ab
                  TNEW(eq) = TNEW(eq) - TMP4(b,i)
               enddo
            enddo
         enddo
      enddo
   enddo

   ! 01. + sum(l,c,d) [ <a,l||c,d>*t(l->d)*t(i->c) ]

   TOLDP => TOLD(1:NSING)
   TNEWP => TNEW(1:NSING)
   call DGEMM('T','N',NPMO,NHMO,NPMO,1.0d00,TMP1, &
              NPMO,TOLDP,NPMO,1.0d00,TNEWP,NPMO)

   ! 02. + sum(l,c,d) [ <a,l||c,d>*t(l->d)*t(i,j->c,b) ] + P(a,b)

   TMP2 => TMP_PP2
   TMP3 => TMP_PP3
   do ij = 1, NHMOSQ
      ti = TIND(ij)
      cb = 0
      do c = 1, NPMO
         TMP2(c,c) = dble(0)
         do b = 1, (c-1)
            cb = cb + 1
            TMP2(c,b) = +TOLD(ti+cb)
            TMP2(b,c) = -TOLD(ti+cb)
         enddo
      enddo
      call DGEMM('T','N',NPMO,NPMO,NPMO,1.0d00,TMP1, &
                 NPMO,TMP2,NPMO,0.0d00,TMP3,NPMO)
      ab = 0
      do a = 2, NPMO
      do b = 1, (a-1)
         ab = ab + 1
         eq = ti + ab
         TNEW(eq) = TNEW(eq) + TMP3(a,b) - TMP3(b,a)
      enddo
      enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Dir_VPH2P

subroutine Dir_VPHPP(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: StartBufInd, EndBufInd
   integer :: Block, St, En
   integer :: i, j, a, b, l
   integer :: ix, ax, il, ij, ab
   integer :: ti, eq
   double precision, dimension(NHMO*NPMO,NHMOSQ) :: THPHSQ
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers

   THPHSQ = dble(0)
   do Block = 1, NumBlock(VPHPP)
      call ReadFile(VPHPP,Block,St,En)
      do i = 1, NHMO
         EndBufInd = 0
         ix = (i-1)*NPMO
         do a = St, En
            ax = (a-1)*NHMO
            do l = 1, NHMO
               StartBufInd = EndBufInd + 1
               EndBufInd   = EndBufInd + NZIntRec(ax+l)
               if ( EndBufInd < StartBufInd ) cycle
               VINT => IntBuf(StartBufInd:EndBufInd)
               VIND => IndBuf(StartBufInd:EndBufInd)
               if (i/=l) then
                  il = MAXMIN(i,l)
                  ti = TIND(il)
                  ! 01. + sum(c>d,l) [ <a,l||c,d>*t(i,l->c,d) ]
                  if (i<l) then
                     TNEW(ix+a) = TNEW(ix+a) - DOT_PRODUCT(VINT,TOLD(ti+VIND))
                  else
                     TNEW(ix+a) = TNEW(ix+a) + DOT_PRODUCT(VINT,TOLD(ti+VIND))
                  endif
                  ! 02. + sum(c) [ <c,j||a,b>*t(i->c) ] + P(i,j) => t(i,j->a,b)
                  !     + sum(a) [ <a,l||c,d>*t(i->a) ] + P(i,l) => t(i,l->c,d)
                  if (i<l) then
                     TNEW(ti+VIND) = TNEW(ti+VIND) - VINT*TOLD(ix+a)
                  else
                     TNEW(ti+VIND) = TNEW(ti+VIND) + VINT*TOLD(ix+a)
                  endif
               endif
               do j = 1, (i-1)
                  ij = MAXMIN(i,j)
                  ti = TIND(ij)
                  THPHSQ(ax+l,ij) = DOT_PRODUCT(VINT,TOLD(ti+VIND))
               enddo
            enddo
         enddo
      enddo
   enddo

   ! 03. - sum(l,c>d) [ <a,l||c,d>*t(i,j->c,d)*t(l->b) ] + P(a|b)

   TMP1 => TMP_PP1
   TOLDP => TOLD(1:NSING)
   do ij = 1, NHMOSQ
      ti = TIND(ij)
      call DGEMM('N','N',NPMO,NPMO,NHMO,1.0d00,TOLDP,NPMO, &
                 THPHSQ(:,ij),NHMO,0.0d00,TMP1,NPMO)
      ab = 0
      do a = 2, NPMO
      do b = 1, (a-1)
         ab = ab + 1
         eq = ti + ab
         TNEW(eq) = TNEW(eq) - TMP1(b,a) + TMP1(a,b)
      enddo
      enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Dir_VPHPP

subroutine Dir_FHH(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: i, j, k
   integer :: ij, kx, ki, kj
   integer :: ti, eq
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers

   FINT => TMP_HH
   read(UFHH,rec=vac) FINT

   ! 01. - sum(k)*[ <i|f|k>*t(k->a) ] 

   TOLDP => TOLD(1:NSING)
   TNEWP => TNEW(1:NSING)
   call DGEMM('N','N',NPMO,NHMO,NHMO,-1.0d00,TOLDP, &
              NPMO,FINT,NHMO,1.0d00,TNEWP,NPMO)

   ! 02. - sum(k)*[ <i|f|k>*t(k,j->a,b) ] + P(i|j)

   ij = 0
   do i = 2, NHMO
   do j = 1, (i-1)
      ij = ij + 1
      eq = TIND(ij)
      TNEWP => TNEW(eq+1:eq+NPMOSQ)
      do k = 1, NHMO
         if (k/=i) then
            ki = MAXMIN(k,i)
            ti = TIND(ki)
            TOLDP => TOLD(ti+1:ti+NPMOSQ)
            if (k>i) then
               TNEWP = TNEWP + FINT(k,j)*TOLDP
            else
               TNEWP = TNEWP - FINT(k,j)*TOLDP
            endif
         endif
         if (k/=j) then
            kj = MAXMIN(k,j)
            ti = TIND(kj)
            TOLDP => TOLD(ti+1:ti+NPMOSQ)
            if (k>j) then
               TNEWP = TNEWP - FINT(k,i)*TOLDP
            else
               TNEWP = TNEWP + FINT(k,i)*TOLDP
            endif
         endif
      enddo
   enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Dir_FHH

subroutine Dir_FPP(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: a, b, c
   integer :: cb, ab, ij
   integer :: ti, eq
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers

   FINT => TMP_PP1
   TMP1 => TMP_PP2
   TMP2 => TMP_PP3

   read(UFPP,rec=vac) FINT

   ! 01. + sum(c)*[ <a|f|c>*t(i->c) ] 

   TOLDP => TOLD(1:NSING)
   TNEWP => TNEW(1:NSING)
   call DGEMM('N','N',NPMO,NHMO,NPMO,1.0d00,FINT, &
              NPMO,TOLDP,NPMO,1.0d00,TNEWP,NPMO)

   ! 02. + sum(c)*[ <a|f|c>*t(i,j->c,b) ] + P(a|b) 

   do ij = 1, NHMOSQ
      ti = TIND(ij)
      ! TMP1(c,b) = t(i,j->c,b)
      cb = 0
      do c = 1, NPMO
         TMP1(c,c) = dble(0)
         do b = 1, (c-1)
            cb = cb + 1
            TMP1(c,b) = +TOLD(ti+cb)
            TMP1(b,c) = -TOLD(ti+cb)
         enddo
      enddo
      call DGEMM('N','N',NPMO,NPMO,NPMO,1.0d00,FINT, &
                 NPMO,TMP1,NPMO,0.0d00,TMP2,NPMO)
      ab = 0
      do a = 2, NPMO
      do b = 1, (a-1)
         ab = ab + 1
         eq = ti + ab
         TNEW(eq) = TNEW(eq) + TMP2(a,b) - TMP2(b,a)
      enddo
      enddo
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return
      
end subroutine Dir_FPP

subroutine Dir_FHP(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: i, j, a, b, k, c
   integer :: ij, ab, ki, ix
   integer :: kj, cb, ca, kx
   integer :: ti, eq
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)
   call NullifyPointers

   FINT => TMP_PH1
   read(UFHP,rec=vac) FINT

   ! 01. + <a|f|i>

   do i = 1, NHMO
      ix = (i-1)*NPMO
      TNEW(ix+1:ix+NPMO) = TNEW(ix+1:ix+NPMO) + FINT(:,i)
   enddo

   ! 02. - sum(k,c)[<k|f|c>*t(i->c)*t(k->a)]
   ! 03. - sum(k,c)[<k|f|c>*t(i->c)*t(k,j->a,b)] + P(i|j)

   TMP1 => TMP_HH
   TOLDP => TOLD(1:NSING)
   TNEWP => TNEW(1:NSING)
   call DGEMM('T','N',NHMO,NHMO,NPMO,1.0d00,FINT,  &
              NPMO,TOLDP,NPMO,0.0d00,TMP1,NHMO)

   ! 02. - sum(k,c) [ <k|f|c>*t(i->c)*t(k->a)

   call DGEMM('N','N',NPMO,NHMO,NHMO,-1.0d00,TOLDP, &
              NPMO,TMP1,NHMO,1.0d00,TNEWP,NPMO)

   ! 03. - sum(k,c) [ <k|f|c>*t(i->c)*t(k,j->a,b) + P(i|j)

   ij = 0
   do i = 2, NHMO
   do j = 1, (i-1)
      ij = ij + 1
      eq = TIND(ij)
      TNEWP => TNEW(eq+1:eq+NPMOSQ)
      do k = 1, NHMO
         if (k/=i) then
            ki = MAXMIN(k,i)
            ti = TIND(ki)
            TOLDP => TOLD(ti+1:ti+NPMOSQ)
            if (k>i) then
               TNEWP = TNEWP + TMP1(k,j)*TOLDP
            else
               TNEWP = TNEWP - TMP1(k,j)*TOLDP
            endif
         endif
         if (k/=j) then
            kj = MAXMIN(k,j)
            ti = TIND(kj)
            TOLDP => TOLD(ti+1:ti+NPMOSQ)
            if (k>j) then
               TNEWP = TNEWP - TMP1(k,i)*TOLDP
            else
               TNEWP = TNEWP + TMP1(k,i)*TOLDP
            endif
         endif
      enddo
   enddo
   enddo

   ! - 04. sum(k,c) [<k|f|c>*t(i,j->c,b)*t(k->a)] + P(a|b)

   TMP1 => TMP_PP1
   TMP2 => TMP_PP2
   TMP3 => TMP_PP3
   TOLDP => TOLD(1:NSING)
   call DGEMM('N','T',NPMO,NPMO,NHMO,1.0d00,FINT,NPMO,  &
              TOLDP,NPMO,0.0d00,TMP1,NPMO)
   do ij = 1, NHMOSQ
      ti = TIND(ij)
      ! TMP2(c,b) = t(i,j->c,b)
      cb = 0
      do c = 1, NPMO
         TMP2(c,c) = dble(0)
         do b = 1, (c-1)
            cb = cb + 1
            TMP2(c,b) = +TOLD(ti+cb)
            TMP2(b,c) = -TOLD(ti+cb)
         enddo
      enddo
      call DGEMM('T','N',NPMO,NPMO,NPMO,1.0d00,TMP1,NPMO,  &
                 TMP2,NPMO,0.0d00,TMP3,NPMO)
      ab = 0
      do a = 2, NPMO
      do b = 1, (a-1)
         ab = ab + 1
         eq = ti + ab
         TNEW(eq) = TNEW(eq) - TMP3(a,b) + TMP3(b,a)
      enddo
      enddo
   enddo

   ! 05. + sum(k,c) [ <k|f|c>*t(k,i->c,a) ]

   TMP1 => TMP_PHP
   TMP2 => TMP_PP1
   do i = 1, NHMO
      ix = (i-1)*NPMO
      ! TMP1(kc,a) = t(k,i->c,a)
      ! k < i case
      do k = 1, (i-1)
         ki = MAXMIN(k,i)
         ti = TIND(ki)
         kx = (k-1)*NPMO
         ca = 0
         do c = 1, NPMO
            TMP2(c,c) = dble(0)
            do a = 1, (c-1)
               ca = ca + 1
               TMP2(c,a) = -TOLD(ti+ca)
               TMP2(a,c) = +TOLD(ti+ca)
            enddo
         enddo
         TMP1(kx+1:kx+NPMO,:) = TMP2
      enddo
      ! k = i case
      TMP1(ix+1:ix+NPMO,:) = dble(0)
      ! k > i case
      do k = (i+1), NHMO
         ki = MAXMIN(k,i)
         ti = TIND(ki)
         kx = (k-1)*NPMO
         ca = 0
         do c = 1, NPMO
            TMP2(c,c) = dble(0)
            do a = 1, (c-1)
               ca = ca + 1
               TMP2(c,a) = +TOLD(ti+ca)
               TMP2(a,c) = -TOLD(ti+ca)
            enddo
         enddo
         TMP1(kx+1:kx+NPMO,:) = TMP2
      enddo
      TNEWP => TNEW(ix+1:ix+NPMO)
      call DGEMV('T',SIZE(TMP1,1),SIZE(TMP1,2),1.0d00,TMP1,  &
                 SIZE(TMP1,1),FINT,1,1.0d00,TNEWP,1)
   enddo

   call NullifyPointers
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return
      
end subroutine Dir_FHP

subroutine NullifyPointers

   implicit none

   nullify(TMP)
   nullify(TMP1)
   nullify(TMP2)
   nullify(TMP3)
   nullify(TMP4)
   nullify(TMP5)
   nullify(VINT)
   nullify(VIND)
   nullify(FINT)
   nullify(TNEWP)
   nullify(TOLDP)

   return

end subroutine NullifyPointers

end module directdiags
