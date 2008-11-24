c
c  RCS $Revision$  $Date$
c
*deck driver_argos
c
c*******************************************************************
c
c   this computer program contains work performed partially by the
c   argonne national laboratory theoretical chemistry group under
c   the auspices of the office of basic energy sciences,
c   division of chemical sciences, u.s. department of energy,
c   under contract w-31-109-eng-38.
c
c   these programs may not be (re)distributed without the
c   written consent of the argonne theoretical chemistry group.
c
c   since these programs are under development, correct results
c   are not guaranteed.
c
c*******************************************************************
c
      subroutine driver_argos( a, lcore, mem1, ifirst )
c
c            argonne - ohio state symmetry-adapted,
c       general-contraction integral program.  argos(i60439)
c  this program calculates one- and two-electron integrals over
c  symmetry orbitals of generally contracted gaussian orbitals and
c  forms the symmetry orbital integrals ordered by pk shell-blocks.
c  the point groups are limited to d2h and its subgroups.
c           by russell m. pitzer
c
c  version log:
c  16-aug-97 recur1 replaced by rad1. -rmp
c  21-may-97 many arrays allocated. -rmp
c            limits on quantum numbers removed. -rmp
c            tables replaced by code to generate them. -rmp/zyz
c  17-sep-96 ctss removed.
c  24-apr-92 rs6000 blocks added. -rls
c  11-sep-91 l2opxv() added. final parallel code added. -rls/rjh
c  01-sep-91 array declarations changed to workaround the ibm+convex
c            f77 equivalence bug. -rls/m. schueller/j. bentley
c  07-jul-91 preliminary parallel version using tcgmsg. -rjh
c  04-feb-91 npru changes in socfpd(). other cleanup. -rmp
c  07-dec-90 il(*) allocated from a(*) workspace. array(*,*) added
c            to oneint(). -rls
c  04-dec-90 title(1:5), only1e, and ndp<0 input options added. -rls
c  01-dec-90 SIFS output file written. -rls
c  28-nov-90 rad2() replaced. -rmp/rls
c  12-feb-90 stellar adaptation (rmp).
c  07-jan-90 change aoso2e, mpru test, eta 2-dim. (mn, rmp).
c  09-jul-89 mdc changes, miscellaneous changes (rmp).
c  07-mar-89 mdc changes, port to anl sun, alliant, titan (rls).
c  21-oct-88 dimension error of nprir and ica in seg1mn fixed (eas).
c  19-sep-88 sun version of code added (eas).
c  02-aug-88 rmp changes, mdc references, bummer calls added (rls).
c  05-jul-88 anl unicos version (rls).
c  09-may-88 timer() updated, bummer() included (rmp/rls).
c  20-apr-88 pseud2() and pseud3() modified (rmp/rls).
c  01-feb-88 cray-ctss and fps memory allocation added (rls).
c  17-jan-88 rest of rmp's corrections--common blocks (dcc)
c  17-dec-87 mblu defined in main, rmp's corrections (dcc,rmp)
c  12-oct-87 report max mblu, (*) arrays, cray pack (dcc)
c  21-sep-87 working under unicos (dcc)
c  15-sep-87 cmdc for osuibm and unicos (dcc)
c  11-aug-87 pkflag() and w*() routines added for integral io (rls).
c  07-aug-87 cmdc version (vax, fps, crayanl) from osu version (rls).
c  16-jul-86 dimension parameters added to main and seg1mn (rls/rmp).
c  02-apr-86 iperm(*,*) added to twoint (rls).
c  25-sep-84 a0 initialized in twoint and qpasy (rls/rab).
c  07-may-84 space allocation changes in socpfd (rls/rmp).
c  31-oct-83 fps-164 version changes (rls).
c  28-oct-83 tape read error corrected in routine recur (rls).
c
c  cmdc info:
c  keyword    description
c  -------   -------------
c  vax       vax code.
c  fps       fps code.
c  cray      cray code. includes "d" exponents in constants.
c  ibm       ibm (mainframe) code.
c  sun       sun code.
c  stellar   stellar specific code.
c  craycos   cos specific code.
c  unicos    cray unicos specific code.
c  mvsosu    uses osu local mvs routines.
c  crayanl   cray anl index packing code.
c  crayosu   cray osu index packing code.
c  ibmosu    ibm osu index packing code.
c  sunosu    sun osu index packing code.
c  vaxanl    vax anl index packing code.
c  vaxosu    vax osu index packing code.
c  crayxmp   cray xmp specific code.
c  parallel  tcgmsg parallel code.
c
c    *** things to do in this program ***
c  * eliminate hollerith strings in format statements.
c  * replace real constants with parameters.
c  * simplify the memory allocation method.
c  * improve vectorization (including so transformation).
c  * eliminate unnecessary statement labels (e.g. use if-then-else
c    blocks) to improve scalar and vector optimization.
c
      implicit logical(a-z)
c
c     # the following parameters should agree exactly with those in the
c     # seg1mn() routine.  (rls)
c
c     # msu    = max number of symmetry inequivalent types of atoms.
      integer    msup
      parameter (msup=203)
c
c     # mstu   = max number of irreps.
      integer    mstup
      parameter (mstup=8)
c
c     # mrcru  = max number of members in a contraction set.
      integer    mrcrup
      parameter (mrcrup=20)
c
c     # mconu  = max number of primitives in a contraction.
      integer    mconup
      parameter (mconup=24)
c
c     # mcu    = max number of symmetry equivalent atoms.
      integer    mcup
      parameter (mcup=12)
c
c     # mconsu = max number of contraction sets.
      integer    mconsp
      parameter (mconsp=65)
c
c     # msfu   = max number of function sets.
      integer    msfup
C      parameter (msfup=52)
      parameter (msfup=70)
c
      integer    mnsfup
      parameter (mnsfup=(msfup*(msfup+1))/2)
c
c     # mgu    = max number of operators in the
c     #          nuclear interchange group.
      integer    mgup
      parameter (mgup=49)
c
c     # mccu   = max number of symmetry unique center
c     #          combinations from 4 nuclei.
      integer    mccup
      parameter (mccup=182)
c
c     # mblu   = workspace array length.
c     #          (allocated in the main program)
c
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
      integer                  aoints
      equivalence ( iunits(2), aoints )
      integer                  aoint2
      equivalence ( iunits(3), aoint2 )
      integer                  ninput
      equivalence ( iunits(4), ninput )
c
c     # filenames.
      character*60    fnames
      common /cfname/ fnames(nunits)
c
      integer         only1e
      common /conly1/ only1e
c
c     # /bufout/ holds some output integral file parameters.
      integer
     & info,     ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & info(10), ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
c
      integer      l2rec
      equivalence( l2rec, info(4) )
      integer      n2max
      equivalence( n2max, info(5) )
c
      real*8         cutoff, tol
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      integer        ipq
      common /dim21/ ipq(256)
c
c     # dummy:
      integer lcore, mem1, ifirst
c
c     # local:
      integer iclock, iargos, ierr, h2o1e, h2o2e
c
      integer
     & nf(msup),      nc(msup),      ncon(mconsp),    nrcr(mconsp),
     & lmnp1(mconsp), ntl(msfup),    ntu(msfup),      nt(msfup),
     & mcons(msfup),  iprst(mnsfup), npair(2,mnsfup), icxst(2,mnsfup),
     & nfct(mccup),   ngw(mccup),    icb(4,24,mccup),
     & ica(mcup,msup,mgup),          nprir(2,mstup,mnsfup)
c
c     # ipt(*) is used for workspace stack pointers.
      integer ipt(55)
c
      real*8
     & a(lcore),       zet(mconup,mconsp),  eta(mrcrup,mconup,mconsp),
     & x(mcup,msup),   y(mcup,msup),        z(mcup,msup)
c
c     # timer operations.
      integer   tmin0
      parameter(tmin0=0)
      integer   tminit,  tmprt,  tmrein,  tmclr,  tmsusp,  tmresm
      parameter(tminit=1,tmprt=2,tmrein=3,tmclr=4,tmsusp=5,tmresm=6)
c
c     # bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
      integer  atebyt, forbyt
      external atebyt, forbyt
c
c      call timer_argos(' ', tmin0, iargos, nlist )
c
      msu   = msup
      mstu  = mstup
      mrcru = mrcrup
      mconu = mconup
      mcu   = mcup
      mccu  = mccup
      mblu  = lcore
c
c     # allocate space for lmnv
      ipt(1) = 1
c
      h2o1e = 0
      h2o2e = 0
c
c      call timer_argos(' ', tminit, iclock, nlist )
c
c     # read input, compute il(*), cx(*), and the 1-e integrals.
c
c     # initialize the high-water stack.
c     # it is not necessary to push, since no workspace has been
c     # allocated, but do it anyway in case allocation is performed
c     # at this point at a later time.
      call h2oini( 0 )
      call h2opsh( 0 )
c
            call seg1mn( a,      eta,       ica,    icb,    icxst,
     &   iprst,  ipt,    lmnp1,  a(ipt(1)), mcons,  nc,     ncon,
     &   nf,     nfct,   ngw,    npair,     nprir,  nrcr,   nt,
     &   ntl,    ntu,    x,      y,         z,      zet,    lcore,
     &   ifirst, mem1 )
c
c     # read the 1-e high-water mark and check for errors.
      call h2opop
      call h2omtr( h2o1e)
      call h2ofin( ierr )
      if ( ierr .ne. 0 ) then
        call bummer('driver: from 1-e h2ofin, ierr=',ierr,wrnerr)
      endif
c
c      call timer_argos('seg1mn required', tmrein, iclock, nlist )
c
*      if ( only1e .ne. 1 ) then
*c
*c       # compute the two-electron integrals.
*c
*c       # assign the output unit ntape, and open the 2-e file.
*c
*        call sifo2f( aoints, iunits(3), fnames(3), info, ntape, ierr )
*        if ( ierr .ne. 0 ) then
*          call bummer('driver: from sifo2f(), ierr=',ierr,faterr)
*        endif
*c
*c       # allocate output record and value buffers.
*c       # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(1:mcxu),
*c       #          4:buffer(1:l2rec), 5:values(1:n2max),
*c       #          6:labels(1:4,1:n2max), 7:ibitv(1:n2max), 8:h4(1:nblu
*c
*c       # note that ibitv(*) is rounded up to a multiple of 64.
*c
*        ipt(5) = ipt(4) + atebyt( l2rec )
*        ipt(6) = ipt(5) + atebyt( n2max )
*        ipt(7) = ipt(6) + forbyt( 4*n2max )
*        ipt(8) = ipt(7) + forbyt( ((n2max+63)/64)*64 )
*c
*c       # mblu is the remaining workspace.
*c       mblu    = lcore - ipt(8) + 1
*c
*        if ( mblu .le. 0 ) then
*          call bummer('driver: mblu=',mblu,faterr)
*        endif
*c
*c       # initialize the high-water stack, and push the current
*c       # workspace allocation onto the stack.
*        call h2oini( 0 )
*        call h2opsh( (ipt(8) - 1) )
*c
*            call twoint( a,         a(ipt(3)), eta,       a(ipt(8)),
*     &        ica,       icb,       icxst,     a(ipt(2)), iprst,
*     &        ipt,       lmnp1,     a(ipt(1)), mcons,     nc,
*     &        ncon,      nf,        nfct,      ngw,       npair,
*     &        nprir,     nrcr,      nt,        ntl,       ntu,
*     &        x,         y,         z,         zet,       a(ipt(4)),
*     &        a(ipt(5)), a(ipt(6)), a(ipt(7)) )
*c
*        call timer_argos('twoint required', tmclr, iclock, nlist )
*c
*c       # close the 2-e file.
*c
*        call sifc2f( ntape, info, ierr )
*        if ( ierr .ne. 0 ) then
*          call bummer('driver: from sifc2f(), ierr=',ierr,faterr)
*        endif
*c
*c       # pop back to the initial stack level, read the high-water
*c       # mark for the twoint() section, and check for errors.
*        call h2opop
*        call h2omtr( h2o2e )
*        call h2ofin( ierr )
*        if ( ierr .ne. 0 ) then
*          call bummer(' driver: from h2ofin, ierr=',ierr,wrnerr)
*        endif
*      endif
c
c     # close the 1-e file.
c
      close(unit=aoints)
c
c     # write out the core allocation high-water marks.
      write(nlist, * )
      write(nlist, 6010 ) '1-e integral ', h2o1e
      write(nlist, 6010 ) '2-e integral ', h2o2e
      write(nlist, 6010 ) 'overall argos', max( h2o1e, h2o2e )
c
c      call timer_argos('argos required', tmclr, iargos, nlist )
c
      return
6010  format(' driver: ',a,' workspace high-water mark =',i10)
      end
*deck seg1mn
      subroutine seg1mn( a,      eta,       ica,    icb,    icxst,
     &   iprst,  ipt,    lmnp1,  lmnv,      mcons,  nc,     ncon,
     &   nf,     nfct,   ngw,    npair,     nprir,  nrcr,   nt,
     &   ntl,    ntu,    x,      y,         z,      zet,    lcore,
     &   ifirst, mem1 )
c
c  read the user input and compute the 1-e integrals.
c
      implicit logical(a-z)
c
c     # the following parameters should agree exactly with those in the
c     # calling program.  (rls)
c
c     # msu    = max number of symmetry inequivalent types of atoms.
      integer    msup
      parameter (msup=203)
c
c     # mstu   = max number of irreps.
      integer    mstup
      parameter (mstup=8)
c
c     # kaords = max number of ao reduction sets.
      integer    kaordp
      parameter (kaordp=240)
c
c     # mconsu = max number of contraction sets.
      integer    mconsp
      parameter (mconsp=65)
c
c     # mgcsu  = max number of symmetry orbital transformation sets.
      integer    mgcsup
      parameter (mgcsup=240)
c
c     # mru    = max number of irreps in an ao reduction set.
      integer    mrup
      parameter (mrup=24)
c
c     # mcsu   = max number of so's in a transformation matrix.
      integer    mcsup
      parameter (mcsup=24)
c
c     # mctu   = max number of ao's in a transformation matrix.
      integer    mctup
      parameter (mctup=30)
c
c     # mcru   = max number of core potential expansion functions.
      integer    mcrup
      parameter (mcrup=140)
c
c     # msfu   = max number of function sets.
      integer    msfup
C      parameter (msfup=52)
      parameter (msfup=70)
c
      integer    mnsfup
      parameter (mnsfup=(msfup*(msfup+1))/2)
c
c     # mgu    = max number of operators in the nuclear
c     #          interchange group.
      integer    mgup
      parameter (mgup=49)
c
c     # msfru  = max number of so sets.
      integer    msfrup
      parameter (msfrup=312)
cbert parameter (msfrup=208)
c
c     # mnru   = max number of charge dists. from a function set pair.
      integer    mnrup
      parameter (mnrup=576)
c
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
      integer                  aoints
      equivalence ( iunits(2), aoints )
c
c     # filenames.
      character*60    fnames
      common /cfname/ fnames(nunits)
c
c     # /bufout/ holds some output integral file parameters.
      integer
     & info,     ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & info(10), ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
c
      integer      l1rec
      equivalence( l1rec, info(2) )
      integer      n1max
      equivalence( n1max, info(3) )
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)
c
      integer        ipq
      common /dim21/ ipq(256)
c
c     # dummy:
      integer lcore, ifirst, mem1
      integer
     & ica(mcu,msu,*), icb(4,24,*), icxst(2,*), iprst(*),
     & ipt(*),         lmnp1(*),    lmnv(3,*),
     & mcons(*),       nc(*),       ncon(*),    nf(*),
     & nfct(*),        ngw(*),      npair(2,*), nprir(2,mstu,*),
     & nrcr(*),        nt(*),       ntl(*),     ntu(*)
c
      real*8
     & a(*),     eta(mrcru,mconu,*), x(mcu,*), y(mcu,*),
     & z(mcu,*), zet(mconu,*)
c
c     # local:
      integer iclock, i, n, narray
      integer
     & lcr(msup),       lls(msup),      nkcrl(6,msup), nkcru(6,msup),
     & nklsl(4,msup),   nklsu(4,msup),  mcrs(msup),    nir(kaordp),
     & maords(mgcsup),  icxsv1(mgcsup), nd(mstup),     nso(mstup),
     & nsopr(mstup),    nblpr(mstup),   nopir(mstup),  ilxyz(3,mstup),
     & la(mrup,kaordp), ncr(mcrup),     nct(msfup),    mgcs(msfup),
     & lb(msfrup),      ms(msfrup),     mnl(msfrup),   mau(mnrup),
     & mics(mnrup),     mjcs(mnrup),    icxast(mnsfup),
     & icxsv2(mgcsup,mgcsup),           idp(mstup,mstup,mstup),
     & ihwsp,           igh,            lmn1,          nn,
     & lmn1u,           ncru,           lproju,        ndfac,
     & lmax,            l1max2,         lmnpwr,        l2m1
c
      real*8 chg(msup), zcr(mcrup), ccr(mcrup), eps
c
      character*3 ityp(mstup), mtype(msup)
c
c     # bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
c     # timer operations.
      integer   tminit,  tmprt,  tmrein,  tmclr,  tmsusp,  tmresm
      parameter(tminit=1,tmprt=2,tmrein=3,tmclr=4,tmsusp=5,tmresm=6)
c
      integer  atebyt, forbyt
      external atebyt, forbyt
c
      kaords = kaordp
      mconsu = mconsp
      mgcsu  = mgcsup
      mru    = mrup
      mcsu   = mcsup
      mctu   = mctup
      mcru   = mcrup
      msfu   = msfup
      mgu    = mgup
      msfru  = msfrup
      mnru   = mnrup
c
c     # lcore = total available workspace on entry to this routine.
      lcore  = mblu
c
      do 20 i = 1, 3
        inam(i) = i
   20 continue
c
c     # set up ipq(*).
c
      do 30 i = 1, 255
        ipq(i) = (i * (i - 1) ) / 2
   30 continue
c
c     # reserve space for c(mctu,mcsu,mgcsu) at the end of a(1:mblu).
c     # afterwards, a(1:mcxu) is the remaining available workspace.
c     # a(*)--> 1:...55:c(1:csize)
c
      ipt(55) = mblu + 1 - atebyt( mctu * mcsu * mgcsu )
c
      mcxu = ipt(55) - ipt(1)
c
      if ( mcxu .lt. 0 ) then
        call bummer('seg1mn: before syminp() mcxu=',mcxu,faterr)
      endif
c
c     # push the current allocation onto the stack.
      call h2opsh( (lcore - mcxu) )
c
c      call timer_argos(' ', tminit, iclock, nlist )
c
c     # input processing.
c
            call syminp_argos( a,      a(ipt(55)), ccr,    chg,    eta,
     &   ica,    idp,    ipt,    ityp,       la,     lb,     lcr,
     &   lls,    lmnp1,  lmnv,   lmn1u,      lproju, maords, mcons,
     &   mcrs,   mgcs,   mnl,    ms,         mtype,  nblpr,  nc,
     &   ncon,   ncr,    ncru,   nct,        nd,     nf,     nir,
     &   nkcrl,  nkcru,  nklsl,  nklsu,      nrcr,   nso,    nsopr,
     &   nt,     ntl,    ntu,    x,          y,      z,      zcr,
     &   zet,    lcore,  ifirst, mem1 )
c
c      call timer_argos('syminp required', tmrein, iclock, nlist )
c
c     # pop back to the current level.
      call h2opop
c
c     # allocate space for il(*), and compute the
c     # symmetry coefficient products.
c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(*)...,
c     #          55:c(*)
c
      ipt(3) = ipt(2) + forbyt( nnbft )
c
c     # mcxu is the remaining space.
      mcxu = ipt(55) - ipt(3)
c
c     # push the current allocation onto the high-water mark stack.
      call h2opsh( (lcore - mcxu) )
c
            call socfpd( a(ipt(55)), a(ipt(3)), icxast, icxst,
     &   icxsv1, icxsv2, idp,        a(ipt(2)), iprst,  la,
     &   lb,     maords, mau,        mcons,     mgcs,   mics,
     &   mjcs,   nblpr,  nc,         nd,        nf,     nir,
     &   npair,  nprir,  nrcr,       nsopr,     nt,     ntl,
     &   ntu  )
c
c     # pop back to the initial level.
      call h2opop
c
c      call timer_argos('socfpd required', tmrein, iclock, nlist )
c
c     # compute the one-electron integrals.
c
c     # on return from socfpd(*), cx(1:mcxu) was computed.
c     # allocate space.
c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(1:mcxu),
c     #          4:array(1:mpru,1:narray), 5:iarray(1:iamax,1:narray),
c     #          6:buffer(1:l1rec), 7:values(1:n1max),
c     #          8:labels(1:2,1:n1max), 9:hpt(1:ihwsp), 10:hwt(1:ihwsp),
c     #          11:dfac(1:ndfac), 12:binom(1:xyzdim), 13:lmf(1:l1max2),
c     #          14:lml(1:l1max2), 15:lmx(1:lmnpwr), 16:lmy(1:lmnpwr),
c     #          17:lmz(1:lmnpwr), 18:zlm(1:lmnpwr),
c     #          19:flmtx(1:3,1:lmnu2), 20:mc(1:3,1:l2m1),
c     #          21:mr(1:3,1:l2m1), 22:h2(1:nblu)
c
      if ( inam(nu) .eq. 5 ) then
c       # x,y,z components of spin-orbit integral are all computed
c       # together.
        narray = 3
      else
        narray = 1
      endif
c
      ihwsp = (lmn1u*(lmn1u+1))/2
      ipt(4)  = ipt(3)  + atebyt( mcxu )
      ipt(5)  = ipt(4)  + atebyt( narray * mpru )
      ipt(6)  = ipt(5)  + forbyt( narray * mpru )
      ipt(7)  = ipt(6)  + atebyt( l1rec )
      ipt(8)  = ipt(7)  + atebyt( n1max )
      ipt(9)  = ipt(8)  + forbyt( 2 * n1max )
      ipt(10) = ipt(9)  + atebyt( ihwsp )
      ipt(22) = ipt(10) + atebyt( ihwsp )
c
c     mblu = lcore - ipt(22) + 1
c
c     # compute gauss-hermite points and weights for s, t, v integrals.
c     # 10**(-12) needed to get 6-point values to sufficient accuracy.
      eps = 1.0d-12
      igh = 0
      do 40 lmn1 = 1, lmn1u
        call hermit(lmn1,a(ipt(9)+igh),a(ipt(10)+igh),eps)
        igh = igh + lmn1
   40 continue
c
      ntape = aoints
c
      write(nlist,*)
c
      do 60 n = 1, nu
c
        if( n.eq.4 ) then
c         # allocate space.
          ipt(10) = ipt(9)  + atebyt( 35 )
          ipt(11) = ipt(10) + atebyt( 35 )
          ndfac = max(4*lmn1u+2*lproju-3,6*lproju+3,4*lmn1u-1,
     &      2*lmn1u+2*lproju+1,4,ncru+4*lmn1u+2*lproju-1)
          ipt(12) = ipt(11) + atebyt( ndfac )
          ipt(13) = ipt(12) + atebyt( (lmn1u*(lmn1u+1))/2 )
          lmax =  max(1,lmn1u-1 + max(lmn1u-1,lproju))
          l1max2 = (lmax+1)**2
          ipt(14) = ipt(13) + forbyt( l1max2 )
          ipt(15) = ipt(14) + forbyt( l1max2 )
          lmnpwr = (((lmax*(lmax+2)*(lmax+4))/3)*(lmax+3)+
     &               (lmax+2)**2*(lmax+4))/16
          ipt(16) = ipt(15) + forbyt( lmnpwr )
          ipt(17) = ipt(16) + forbyt( lmnpwr )
          ipt(18) = ipt(17) + forbyt( lmnpwr )
          ipt(19) = ipt(18) + atebyt( lmnpwr )
          ipt(20) = ipt(19) + atebyt( 3*lproju**2 )
          l2m1 = 2*lproju-1
          ipt(21) = ipt(20) + forbyt( 3*l2m1 )
          ipt(22) = ipt(21) + forbyt( 3*l2m1 )
          call cortab(a(ipt(12)),a(ipt(11)),eps,a(ipt(19)),a(ipt(9)),
     &      a(ipt(10)),a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),
     &      a(ipt(17)),lmax,lmn1u,lproju,a(ipt(20)),a(ipt(21)),ndfac,
     &      a(ipt(18)))
        endif
c
c       # push the current allocation onto the high-water mark stack.
        call h2opsh( (ipt(22) - 1) )
c
        nnam  = inam(n)
c
            call oneint_argos( a,       ccr,     chg,        a(ipt(3)),
     &        eta,       a(ipt(22)), ica,
     &        icb,       icxast,     icxst,     a(ipt(2)),  ilxyz,
     &        iprst,     ipt,        lcr,       lls,        lmnp1,
     &        lproju,    mcons,      mcrs,      nc,         ncon,
     &        ncr,       nf,         nfct,      ngw,        nkcrl,
     &        nkcru,     nklsl,      nklsu,     nopir,      npair,
     &        nprir,     nrcr,       nt,        ntl,        ntu,
     &        x,         y,          z,         zcr,        zet,
     &        a(ipt(4)), a(ipt(5)),  a(ipt(6)), a(ipt(7)),  a(ipt(8)) )
c
c       # pop the high-water mark stack.
        call h2opop
   60 continue
      write(nlist,*)
c      call timer_argos('oneint required', tmclr, iclock, nlist )
c
      return
      end
*deck syminp_argos
      subroutine syminp_argos( a,      c,          ccr,    chg,    eta,
     &   ica,    idp,    ipt,    ityp,       la,     lb,     lcr,
     &   lls,    lmnp1,  lmnv,   lmn1u,      lproju, maords, mcons,
     &   mcrs,   mgcs,   mnl,    ms,         mtype,  nblpr,  nc,
     &   ncon,   ncr,    ncru,   nct,        nd,     nf,     nir,
     &   nkcrl,  nkcru,  nklsl,  nklsu,      nrcr,   nso,    nsopr,
     &   nt,     ntl,    ntu,    x,          y,      z,      zcr,
     &   zet,    lcore,  ifirst, mem1 )
c
c  read and verify the user input: symmetry info, basis sets,
c  geometry, etc.
c
c  11-sep-91 parallel code added. -rls/rjh
c  12-jun-91 l2opxv() call added. -rls
c  01-mar-91 aoint2,fsplit,l1rec,l2rec input order changed. -rls/rmp
c
      implicit logical(a-z)
c
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
      integer                  aoints
      equivalence ( iunits(2), aoints )
      integer                  aoint2
      equivalence ( iunits(3), aoint2 )
      integer                  ninput
      equivalence ( iunits(4), ninput )
c
c     # filenames.
      character*60    fnames
      common /cfname/ fnames(nunits)
c
c     # /bufout/ holds some output integral file parameters.
      integer
     & info,     ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & info(10), ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
c
      integer         only1e
      common /conly1/ only1e
c
      real*8         cutoff, tol
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)
c
      integer        ipq
      common /dim21/ ipq(256)
c
c     # dummy:
      integer lcore, ifirst, mem1
      integer
     & ica(mcu,msu,*), idp(mstu,mstu,*), ipt(*),     la(mru,*),
     & lb(*),          lcr(*),           lls(*),     lmnp1(*),
     & lmnv(3,*),      lmn1u,            lproju,     maords(*),
     & mcons(*),       mcrs(*),          mgcs(*),    mnl(*),
     & ms(*),          nblpr(*),         nc(*),      ncon(*),
     & ncr(*),         ncru,             nct(*),     nd(*),
     & nf(*),          nir(*),           nkcrl(6,*), nkcru(6,*),
     & nklsl(4,*),     nklsu(4,*),       nrcr(*),    nso(*),
     & nsopr(*),       nt(*),            ntl(*),     ntu(*)
      real*8
     & a(*),               c(mctu,mcsu,*), ccr(*),   chg(*),
     & eta(mrcru,mconu,*), x(mcu,*),       y(mcu,*), z(mcu,*),
     & zcr(*),             zet(mconu,*)
      character*3 ityp(*), mtype(*)
c
c     # local:
      integer i, ntitle, itol, icut, inrm, ncrs, fsplit, l1rec, l2rec,
     & ngen, ncons, naords, ist, ndpt, idpt, jst, kst, iprd, ierr,
     & iaords, iru, ir, igcs, icsu, ictu, ics, ict, ic, icons, iconu,
     & ircru, icon, ircr, kcru, icrs, lcru, llsu, lp1, nbfcr, kcrl,
     & kcr, l, ngenp1, isf, isfr, is, icu, ig, if, ilmnp1, lai,
     & nsf, itl, itu, jcts, it, jct, jt, ixyz, i2, i3, jconu, jsfr,
     & jcon, js, jcu, jc, ngp1, icsr, ncol, iso, icol, ibvtyp, l1recx,
     & l2recx, n1max, n2max, iau, ia, leig, lmn, lmnvmx, ircrl, ircrh
      integer ixyzir(3), lbla(16), lblso(16), iprdsm(3,7), ndptsm(0:9)
      real*8 twoc, prtint, sum, t, repnuc, chgprd, c1, ccc, eval
c
      character*80 title(5)
      character*3 ltyp(16)
      character*1 ibl, ipc, isk, lblsh(21)
c
      real*8            a0,       a1s2,       a1,       a2,
     & aln10,           a4,       a10
      parameter(        a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a2=2.0d0,
     & aln10=2.30258d0, a4=4.0d0, a10=10.0d0 )
c
c     # lrcmxp is the program default 1-e and 2-e record length.
c     # there is also an SIFS default, which may be different.
      integer    lrcmxp
      parameter( lrcmxp=4096 )
c
c     # bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
c     # for workspace allocation.
      integer  atebyt, forbyt
      external atebyt, forbyt
c
      integer  iodfac, nodeid
      external iodfac, nodeid
c
      data ibl, isk /' ','0'/, lblsh /'s','p','d','f','g','h','i','k',
     &  'l','m','n','o','q','r','t','u','v','w','x','y','z'/
c
c     # nst=4 and nst=8 product triples.
      data iprdsm / 2,3,4, 2,5,6, 2,7,8, 3,5,7, 3,6,8, 4,5,8, 4,6,7 /
c
c     # the number of product triples for each nst.
      data (ndptsm(i),i=0,9) /
     & -1, -1, 0, -1, 1, -1, -1, -1, 7, -1 /
c
c     # assign the default unit numbers and filenames.
      ninput = 5
      fnames(4)='argosin'
      nlist  = 6 
      fnames(1)='Dirac output file'
      aoints = 4
      fnames(2)='aoints'
      aoint2 = 8
      fnames(3)='aoints2'
c
c     # translate filenames if necessary.
c
*mdc*if parallel
*c     # Everyone has their own aoints_xxx, aoints2_xxx and argosls_xxx.
*c     # if not node 0 also rename argosin to avoid clashes.
*      call pfname( 3, fnames )
*      if ( nodeid() .ne. 0 ) call pfname(1,fnames(4))
*c     # copy the input file from node=0 to the other nodes.
*      call pfcopy( 999, 0, fnames(4) )
*c-debug      call flush(6)
*mdc*endif
c
c     # open the input file.
c
      open(unit=ninput,status='old',file=fnames(4))
c
c     # open the output listing file.
c
*mdc*if ibm
*c  use preconnected unit 6
*mdc*else
c      open(unit=nlist,status='unknown',file=fnames(1))
*mdc*endif
      call ibummr(nlist)
c
      write(nlist,'(/28x,a)')'================================='
      write(nlist,'(28x,a)') '!!       Starting  argos       !!'
      write(nlist,'(28x,a)') '!!             for             !!'
      write(nlist,'(28x,a)') '!!       ECP Calculation       !!'
      write(nlist,'(28x,a/)')'================================='

c     # echo the input cards.
      write(nlist,'(a)')'echo of the argos input file:'
      call echoin( ninput, nlist, ierr )
      if ( ierr .ne. 0 ) then
        call bummer('syminp: from echoin, ierr=', ierr, faterr )
      endif
c
      rewind ninput
c
c     # read the first title here.  remaining titles are read later.
c     # *** the next time the input to this program is changed, all
c     #     of the titles should be read together. ***
c
      read (ninput,'(a)') title(1)
c
      ntitle = 1
c
c     # intitalize some input variables.
c     # [01-dec-90 aoint2, fsplit, l1rec, l2rec added. -rls]
      itol   =  0
      icut   =  0
      only1e =  0
      inrm   =  0
      ncrs   =  0
      l1rec  = -1
      l2rec  = -1
      fsplit =  2
c
      read (ninput,*)  ngen,   ns,     naords, ncons,  ngcs,   itol,
     & icut,   aoints, only1e, inrm,   ncrs,   l1rec,  l2rec,  aoint2,
     & fsplit
c
      if ( ns .gt. msu ) then
        call bummer('change msup (two places) to ',ns,faterr)
      elseif ( naords .gt. kaords ) then
        call bummer('change kaordp (one place) to ',naords,faterr)
      elseif ( ncons .gt. mconsu ) then
        call bummer('change mconsp (two places) to ',ncons,faterr)
      elseif ( ngcs .gt. mgcsu ) then
        call bummer('change mgcsup (one place) to ',ngcs,faterr)
      endif
c
c     # reset aoints, aoint2, only1e, and fsplit if necessary.
      if ( aoints .le. 0 ) aoints = 4
      if ( aoint2 .le. 0 ) aoint2 = 8
      if ( fsplit .ne. 1 ) fsplit = 2
*mdc*if parallel
*      if ( fsplit .ne. 2 ) then
*        write(nlist, *) 'Need fsplit = 2 in parallel...resetting'
*        fsplit = 2
*      endif
*mdc*endif
      if ( only1e .ne. 0 ) then
c       # set only1e=0 for off, and only1e=1 for on.
        only1e = 1
        if ( fsplit .eq. 1 ) then
c         # this is inconsistent.  however, instead of treating
c         # this as an error, reset fsplit and continue.  the only
c         # harm is that the user may have to recompute the 1-e
c         # integrals later if this was not the right choice.
          write(nlist,*)'only1e is inconsistent with fsplit.'
     &     //' resetting fsplit=2 and continuing...'
          fsplit = 2
        endif
      endif
c
c     # read number of irreps, degeneracy, and character labels.
c
c     # *** the next time the input to this program is changed,
c     #     all of the items in this list should be read with
c     #     list-directed i/o.  ***
c
      read (ninput,'(i3,12(i3,a3))') nst,(nd(ist),ityp(ist),ist=1,nst)
c
      if ( nst .gt. mstu ) then
        call bummer('change mstup (two places) to ',nst,faterr)
      endif
c
      call izero_argos( nst,           nso, 1 )
      call izero_argos( mstu*mstu*nst, idp, 1 )
c
c     # assign the totally symmetric irrep binary products.
      do 104 ist = 1,nst
        idp(ist,ist,1) = 1
        idp(ist,1,ist) = 1
        idp(1,ist,ist) = 1
  104 continue
c
c     # read the number of product triples.
c
      ndpt = -1
      read(ninput,*) ndpt
c
      if ( ndpt .ge. 0 ) then
c
c       # read in the totally symmetric irrep product
c       # triples from the input file.
c
        do 108 idpt = 1, ndpt
c
          read (ninput,*) ist, jst, kst
c
          idp(ist,jst,kst) = 1
          idp(ist,kst,jst) = 1
          idp(jst,ist,kst) = 1
          idp(jst,kst,ist) = 1
          idp(kst,ist,jst) = 1
          idp(kst,jst,ist) = 1
  108   continue
      else
c
c       # assign default values based on nst.
c       # 04-dec-90 ndp<0 option added. -rls
c
        ndpt = ndptsm( min( max( 0, nst), 9 ) )
        if ( ndpt .lt. 0 ) then
          call bummer('syminp: unknown ndpt defaults for nst=',
     &     nst, faterr )
        else
          do 112 iprd = 1, ndpt
            ist = iprdsm(1,iprd)
            jst = iprdsm(2,iprd)
            kst = iprdsm(3,iprd)
c
            idp(ist,jst,kst) = 1
            idp(ist,kst,jst) = 1
            idp(jst,ist,kst) = 1
            idp(jst,kst,ist) = 1
            idp(kst,ist,jst) = 1
            idp(kst,jst,ist) = 1
  112     continue
        endif
      endif
c
      do 116 iaords = 1, naords
c
        read (ninput,*) iru, (la(ir,iaords), ir=1,iru)
c
        if ( iru .gt. mru ) then
          call bummer('change mrup (one place) to ',iru,faterr)
        endif
        nir(iaords) = iru
  116 continue
c
      do 128 igcs = 1,ngcs
c
        read (ninput,*) icsu, ictu, maords(igcs)
c
        if ( icsu .gt. mcsu ) then
          call bummer('change mcsup (one place) to ',icsu,faterr)
        elseif ( ictu .gt. mctu ) then
          call bummer('change mctup (one place) to ',ictu,faterr)
        endif
c
        do 124 ics = 1,icsu
c
          read (ninput,*) (c(ict,ics,igcs), ict=1,ictu)
c
          do 120 ict = 1,ictu
            c1 = c(ict,ics,igcs)
            ccc = abs(c1)
            if ( ccc .ne. a0 .and. ccc .ne. a1 ) then
               c(ict,ics,igcs) = sign(sqrt(ccc),c1)
            endif
  120     continue
  124   continue
c
  128 continue
c
c     # read in exponents and contraction coefficients.
c
      do 136 icons = 1,ncons
c
        read (ninput,*) iconu, lmnp1(icons), ircru
c
        if ( ircru .eq. 0 ) ircru = 1
        if ( ircru .gt. mrcru ) then
          call bummer('change mrcrup (one place) to ',ircru,faterr)
        elseif ( iconu .gt. mconu ) then
          call bummer('change mconup (one place) to ',iconu,faterr)
        endif
        nrcr(icons) = ircru
        ncon(icons) = iconu
c
        do 132 icon = 1,iconu
c
          read (ninput,*) zet(icon,icons), (eta(ircr,icon,icons),
     &     ircr = 1,ircru)
c
  132   continue
  136 continue
c
      if ( ncrs .eq. 0 ) go to 160
c
c     # read in expansions for core and spin-orbit potentials.
c
      ncru = -2
      lproju = 0
      kcru = 0
      do 156 icrs = 1,ncrs
c
        read (ninput,*) lcru, llsu
        lproju = max(lproju, lcru-1, llsu)
c
c       if ( lcru .gt. 5 ) then
c         call bummer('lcru too large ',lcru,faterr)
c       elseif ( llsu .gt. 4 ) then
c         call bummer('llsu too large ',llsu,faterr)
c       endif
        lcr(icrs) = lcru
        lls(icrs) = llsu
c
        do 144 lp1 = 1,(lcru+1)
c
          read (ninput,*) nbfcr
c
          kcrl = kcru+1
          kcru = kcru+nbfcr
          nkcrl(lp1,icrs) = kcrl
          nkcru(lp1,icrs) = kcru
          do 140 kcr = kcrl,kcru
c
            read (ninput,*) ncr(kcr), zcr(kcr), ccr(kcr)
            if ( ncr(kcr) .lt. 0 ) then
              call bummer('ncr(kcr) too small ',ncr(kcr),faterr)
            endif
            ncru = max(ncru,ncr(kcr))
c
  140     continue
  144   continue
c
        do 152 l = 1,llsu
c
          read (ninput,*) nbfcr
c
          kcrl = kcru+1
          kcru = kcru+nbfcr
          nklsl(l,icrs) = kcrl
          nklsu(l,icrs) = kcru
          do 148 kcr = kcrl,kcru
c
            read (ninput,*) ncr(kcr), zcr(kcr), ccr(kcr)
cisl-b for the original code
c      uncomment the following line
c
c            ccr(kcr)=ccr(kcr)/dfloat(l)
c
cisl-e
            if ( ncr(kcr) .lt. -2 ) then
              call bummer('ncr(kcr) too small ',ncr(kcr),faterr)
            endif
            ncru = max(ncru,ncr(kcr))
c
  148     continue
  152   continue
  156 continue
      if ( kcru .gt. mcru ) then
        call bummer('change mcrup (one place) to ',kcru,faterr)
      endif
c
  160 continue
c
      ngenp1 = ngen + 1
      ng     = ngenp1
c
c     # read in atomic labels, charges, coordinates, basis sets,
c     # and core and spin-orbit potentials
c
      call izero_argos(ns,mcrs,1)
      lmn1u = 0
      isf = 0
      isfr = 0
      do 184 is = 1,ns
c
c       # *** the next time the input to this program is changed,
c       #     all of the items in this list should be read with
c       #     list-directed i/o.  ***
c
c       read (ninput,'(a3,2i3,f3.0)') mtype(is),nf(is),nc(is),chg(is)
        read (ninput,*) mtype(is),nf(is),nc(is),chg(is)
c
        icu = nc(is)
        if ( icu .gt. mcu ) then
          call bummer('change mcup (one place) to ',icu,faterr)
        endif
        do 164 ic = 1,icu
c
          read (ninput,*) x(ic,is), y(ic,is), z(ic,is)
c
  164   continue
        if ( icu .gt. 1 ) then
          do 168 ig = 2, ng
c
            read (ninput,*) (ica(ic,is,ig), ic=1,icu)
c
  168     continue
        endif
c
        do 180 if = 1,nf(is)
          isf = isf+1
c
          read (ninput,*) icons, igcs
c
          mcons(isf) = icons
          iaords = maords(igcs)
          iru =  nir(iaords)
          mgcs(isf) = igcs
          ilmnp1 = lmnp1(icons)
          lmn1u = max(lmn1u,ilmnp1)
          nt(isf) = ipq(ilmnp1+1)
          ntu(isf) = (nt(isf)*(ilmnp1+2))/3
          ntl(isf) = ntu(isf)-nt(isf)+1
          nct(isf) = icu*nt(isf)
          ircru = nrcr(icons)
          do 176 ircr = 1,ircru
            do 172 ir = 1,iru
              isfr = isfr+1
              lai = la(ir,iaords)
              nso(lai) = nso(lai)+1
              lb(isfr) = nso(lai)
  172       continue
  176     continue
  180   continue
        if ( ncrs .gt. 0 ) then
          read(ninput,*) mcrs(is)
        endif
  184 continue
      lmnvmx = (lmn1u*(lmn1u + 1)*(lmn1u + 2))/6
      ipt(2) = ipt(1) + forbyt( 3*lmnvmx )
c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx)
      call lmnvgn(lmn1u,lmnv)
c
c     # read up to three more titles.
c     # stop reading at eof or a blank line.
c     #      04-dec-90 extra titles added. -rls
c
      do 188 i = 2, 4
c
        read(ninput,'(a)',iostat=ierr) title(i)
c
        if ( ierr .eq. 0 ) then
          if ( title(i) .ne. ' ') then
            ntitle = i
          else
c           # logical eof exit is ok.
            go to 192
          endif
        elseif ( ierr .gt. 0 ) then
          call bummer('syminp: title(2:4) input error, ierr=',
     &     ierr,faterr)
        elseif ( ierr .lt. 0 ) then
c         # eof exit is ok.
          go to 192
        endif
  188 continue
c
c     # all done with the input file.
  192 continue
      close(unit=ninput)
c
c     # generate an additional title with a job and tdy stamp.
c
      ntitle = ntitle + 1
#ifdef ARGOS_STANDALONE
      title(ntitle)(1:40) = 'aoints SIFS file created by argos.'
      call siftdy( title(ntitle)(41:) )
#endif
c
      nsf = isf
      if ( nsf .gt. msfu ) then
        call bummer('change msfup (two places) to ',nsf,faterr)
      elseif ( isfr .gt. msfru ) then
        call bummer('change msfrup (one place) to ',isfr,faterr)
      endif
c
c     # nbft = total number of symmetry orbital basis functions.
c
      nbft = isfr
c
c     # determine the number of types of one-electron integrals.
c
c     # initialize nu for s(*),t(*),v(*).
      nu = 3
c
      if ( ncrs .gt. 0 ) then
        do 196 is = 1,ns
          if ( mcrs(is) .ne. 0 ) then
            if ( lcr(mcrs(is)) .gt. 0 ) then
              nu = 4
              inam(4) = 4
              go to 200
            endif
          endif
  196   continue
c
  200   continue
        do 204 is = 1, ns
          if ( mcrs(is) .ne. 0 ) then
            if ( lls(mcrs(is)) .gt. 0 ) then
              nu = nu + 1
              inam(nu) = 5
              go to 208
            endif
          endif
  204   continue
c
  208   continue
      endif
c
      if ( inrm .ne. 0 ) then
c
c       # normalization of symmetry orbitals.
c
        do 244 igcs = 1,ngcs
          isf = 0
          do 240 is = 1,ns
            icu = nc(is)
            do 236 if = 1,nf(is)
              isf = isf+1
              if ( igcs .ne. mgcs(isf) ) go to 236
              itl = ntl(isf)
              itu = ntu(isf)
              iaords = maords(igcs)
              iru = nir(iaords)
              ilmnp1 = lmnp1(mcons(isf))
              ics = 0
              do 232 ir = 1,iru
                iau = nd(la(ir,iaords))
                do 228 ia = 1,iau
                  ics = ics+1
                  prtint = a0
                  ict = 0
                  do 220 ic = 1,icu
                    jcts = ict
                    do 216 it = itl,itu
                      ict = ict+1
                      if ( it .gt. itl ) then
                        twoc = a2*c(ict,ics,igcs)
                        jct = jcts
                        do 212 jt = itl,it-1
                          jct = jct+1
                          prtint = prtint + twoc*c(jct,ics,igcs)
     &                     * (iodfac(lmnv(1,it),lmnv(1,jt),lmnv(2,it),
     &                               lmnv(2,jt),lmnv(3,it),lmnv(3,jt)))
  212                   continue
                      endif
                      prtint = prtint + c(ict,ics,igcs)**2
     &                 * (iodfac(lmnv(1,it),lmnv(1,it),lmnv(2,it),
     &                           lmnv(2,it),lmnv(3,it),lmnv(3,it)))
  216                 continue
  220             continue
                  if ( prtint .ne. a1 ) then
                    prtint = sqrt(prtint)
                    do 224 jct = 1,ict
                      c(jct,ics,igcs) = c(jct,ics,igcs)/prtint
  224               continue
                  endif
  228           continue
  232         continue
              go to 244
  236       continue
  240     continue
  244   continue
      endif
c
      call izero_argos(3,ixyzir,1)
c
      if ( icut .le. 0 ) icut = 9
      cutoff = a1 / a10**icut
c
      if ( itol.le.0) itol = 20
      tol = aln10 * itol
c
c     # compute contracted orbital and lower-triangle-packed
c     # matrix symmetry offsets.
      nsopr(1) = 0
      nblpr(1) = 0
      do 288 ist = 2,nst
        nsopr(ist) = nsopr(ist-1) + nso(ist-1)
        nblpr(ist) = nblpr(ist-1) + ipq( nso(ist-1) + 1 )
  288 continue
c
c     # nnbft = total number of unique so function pairs.
c     # mpru = the number of elements in a diagonal,
c     #        symmetry-blocked, lower-triangle-packed matrix.  this
c     #        is also an upper bound to the size of any
c     #        symmetry-blocked array in which only unique
c     #        off-diagonal elements are stored and which is
c     #        associated with an operator that transforms as an
c     #        irreducible representation.
c
      nnbft = (nbft * (nbft + 1)) / 2
      mpru  = nblpr(nst) + ipq( nso(nst) + 1 )
c
c     # normalize the contraction coefficients.
c
      do 316 icons = 1, ncons
        iconu = ncon(icons)
        ircru = nrcr(icons)
        do 312 ircr =  1, ircru
          if ( iconu .eq. 1 ) then
            eta(1,1,icons) = a1
            go to 312
          endif
          prtint = eta(ircr,1,icons)**2
          do 304 icon = 2, iconu
            twoc  = a2 * eta(ircr,icon,icons)
            jconu = icon - 1
            do 300 jcon = 1, jconu
              t = (zet(icon,icons) + zet(jcon,icons)) * a1s2
              t = zet(icon,icons)/t * zet(jcon,icons)/t
              prtint = prtint + twoc * eta(ircr,jcon,icons)
     &         * sqrt(t**lmnp1(icons) * sqrt(t))
  300       continue
            prtint = prtint + eta(ircr,icon,icons)**2
  304     continue
          prtint = sqrt(prtint)
          do 308 icon = 1, iconu
            eta(ircr,icon,icons) = eta(ircr,icon,icons) / prtint
  308     continue
  312   continue
  316 continue
c
c     # calculate the nuclear repulsion energy.
c
      repnuc = a0
      do 332 is = 1, ns
        ic = nc(is)
        t = (ic) * chg(is)
c       # rearrangement theorem eliminates the ic loop in favor
c       # of the overall (ic) factor.
        do 328 js = 1, is
          chgprd = t * chg(js)
          if ( js .eq. is ) then
c           # repulsions within a set are double-counted, so
c           # compensate with a factor of 1/2.
            chgprd = a1s2 * chgprd
c           # don't include nc(is)=nc(js) self-repulsions.
            jcu = nc(js) - 1
          else
            jcu = nc(js)
          endif
          do 324 jc = 1, jcu
            repnuc = repnuc + chgprd / sqrt( (x(ic,is)-x(jc,js))**2
     1       + (y(ic,is)-y(jc,js))**2 + (z(ic,is)-z(jc,js))**2 )
  324     continue
  328   continue
  332 continue
c
      call gcentr( ica, nc )
c
      ngp1 = ng+1
      if ( ngp1 .gt. mgu ) then
        call bummer('change mgup (two places) to ',ngp1,faterr)
      endif
c
      write (nlist,600)
  600 format(30x,'program "argos" 4.1b2'/
     &  29x,'columbus program system'//
     &  13x,'this program computes integrals over symmetry orbitals'/
     &  17x,'of generally contracted gaussian atomic orbitals.'/
     &  25x,'programmed by russell m. pitzer'//
     &  29x,'version date: 01-sep-97'//' references:'//
     &  t5,'symmetry analysis (equal contributions):'/
     &  t5,'r. m. pitzer, j. chem. phys. 58, 3111 (1973).'//
     &  t5,'ao integral evaluation (hondo):'/
     &  t5,'m. dupuis, j. rys, and h. f. king, j. chem. phys. 65,',
     &     ' 111 (1976).'//
     &  t5,'general contraction of gaussian orbitals:'/
     &  t5,'r. c. raffenetti, j. chem. phys. 58, 4452 (1973).'//
     &  t5,'core potential ao integrals (meldps):'/
     &  t5,'l. e. mcmurchie and e. r. davidson, j. comput. phys.',
     &     ' 44, 289 (1981)'//
     &  t5,'spin-orbit and core potential integrals:'/
     &  t5,'r. m. pitzer and n. w. winter, int. j. quantum chem.',
     &     ' 40, 773 (1991)')
c
c     # print out the name and address of the local programmer
c     # responsible for maintaining this program.
c
      call who2c( 'ARGOS', nlist )
c
      write (nlist,604) lcore, mem1, ifirst
  604 format(/' workspace allocation parameters: lcore=',i10,
     & ' mem1=',i10,' ifirst=',i10)
c
      write (nlist,608) (iunits(i),fnames(i),i=1,nunits)
  608 format(/' filenames and unit numbers:'/
     & ' unit description',                      t37,'filename'/
     & 1x,4('-'),1x,11('-'),                     t37,10('-')/
     & 1x,i3,2x,'listing file:',                 t37,a/
     & 1x,i3,2x,'1-e integral file:',            t37,a/
     & 1x,i3,2x,'2-e integral file [fsplit=2]:', t37,a/
     & 1x,i3,2x,'input file:',                   t37,a )
c
      write (nlist,612) ngen,   ns,     naords, ncons,  ngcs,   itol,
     &  icut,   aoints, only1e, inrm,   ncrs,   l1rec,  l2rec,  aoint2,
     &  fsplit
  612 format(/' argos input parameters and titles:'/
     & ' ngen   =',i3,' ns     =',i3,' naords =',i3,' ncons  =',i3,
     & ' ngcs   =',i3,' itol   =',i3/
     & ' icut   =',i3,' aoints =',i3,' only1e =',i3,' inrm   =',i3,
     & ' ncrs   =',i3/
     & ' l1rec  =',i10,           5x,' l2rec  =',i10,           5x,
     & ' aoint2 =',i3,' fsplit =',i3/)
c
      write (nlist,'(a)') (title(i),i=1,ntitle)
      write (nlist,616) (ist, ist=1,nst)
  616 format(//'irrep',9x,12(i4,4x))
      write (nlist,620) (nd(ist), ist=1,nst)
  620 format('degeneracy',12i8)
      write (nlist,624) (ityp(ist), ist=1,nst)
  624 format('label',11x,12(a3,5x))
      write (nlist,628)
  628 format(//'direct product table')
      do 644 ist = 1,nst
        do 640 jst = 1,nst
          do 636 kst = 1,nst
            if ( idp(ist,jst,kst) .eq. 1 )
     1       write (nlist,632) ityp(ist), ityp(jst), ityp(kst)
  632        format(3x,'(',a3,') (',a3,') = ',a3)
  636     continue
  640   continue
  644 continue
c
      write (nlist,648) repnuc, itol, icut, aoints, aoint2
  648 format(//21x,'nuclear repulsion energy',f14.8///'primitive ao ',
     & 'integrals neglected if exponential factor below 10**(-',i2,')',/
     & 'contracted ao and so integrals neglected if value below 10**(-',
     & i2,')'/'symmetry orbital integrals written on units',2i3)
c
      isf = 0
      isfr = 0
      do 840 is = 1,ns
        icu = nc(is)
        write (nlist,656) mtype(is), chg(is), (ic, x(ic,is), y(ic,is),
     &   z(ic,is), ic=1,icu)
  656   format(//36x,a3,' atoms'//30x,'nuclear charge',f7.2//
     &    11x,'center',12x,'x',15x,'y',15x,'z'/(i14,7x,3f16.8))
        if ( icu .gt. 1 ) then
          write (nlist,660) (ica(ic,is,1), ic=1,icu)
  660     format(' operator',6x,'center interchanges'/
     &           '  1(id.)',6x,24i3)
          do 668 ig = 2,ngenp1
            write (nlist,664) ig, (ica(ic,is,ig), ic=1,icu)
  664       format(i3,'(gen.)',5x,24i3)
  668     continue
          do 672 ig = ngenp1+1,ng
            write (nlist,'(i3,11x,24i3)') ig, (ica(ic,is,ig), ic=1,icu)
  672     continue
        endif
        do 796 if = 1,nf(is)
          isf = isf+1
          itl = ntl(isf)
          itu = ntu(isf)
          igcs = mgcs(isf)
          iaords = maords(igcs)
          iru = nir(iaords)
          ictu = nct(isf)
          icons = mcons(isf)
          ircru = nrcr(icons)
          iconu = ncon(icons)
          icsu = 0
          do 704 ir = 1,iru
            icsu = icsu + nd(la(ir,iaords))
  704     continue
          write (nlist,708) lmnp1(icons), lblsh(lmnp1(icons))
  708     format(/i21,a1,' orbitals'//
     &      ' orbital exponents  contraction coefficients')
c         # use larger value than 4 for output lines longer than 80
          ircrh = min(4,ircru)
          do 712 icon = 1,iconu
            write (nlist,'(8(1pg16.7))') zet(icon,icons),
     &        (eta(ircr,icon,icons), ircr=1,ircrh)
  712     continue
          write (nlist,"(/21x,'symmetry orbital labels')")
          jsfr = isfr
          do 720 ir = 1,iru
            jsfr = jsfr + 1
            lai = la(ir,iaords)
            do 716 ia = 1,nd(lai)
              write (nlist,'(10x,7(i12,a3,i1))') (lb(jsfr+(ircr-1)*
     &          icsu),ityp(lai),ia, ircr=1,ircrh)
  716       continue
  720     continue
  724     if(ircrh.lt.ircru) then
            ircrl = ircrh + 1
            ircrh = min(ircrh+4,ircru)
            write (6,"(/20x,'contraction coefficients')")
            do 728 icon = 1,iconu
              write (nlist,'(16x,7(1pg16.7))')
     &          (eta(ircr,icon,icons), ircr=ircrl,ircrh)
  728       continue
            write (nlist,"(/21x,'symmetry orbital labels')")
            jsfr = isfr
            do 736 ir = 1,iru
              jsfr = jsfr + 1
              lai = la(ir,iaords)
              do 732 ia = 1,nd(lai)
                write (nlist,'(10x,7(i12,a3,i1))') (lb(jsfr+(ircr-1)*
     &            icsu),ityp(lai),ia, ircr=ircrl,ircrh)
  732         continue
  736       continue
            go to 724
          endif
          write (nlist,"(/11x,'symmetry orbitals')")
          icsr = 0
          ncol = 0
          do 760 ir = 1,iru
            lai = la(ir,iaords)
            iau = nd(lai)
            do 756 ia = 1,iau
              ncol = ncol+1
              ltyp(ncol) = ityp(lai)
              lbla(ncol) = ia
              if ( ncol .ne. 10 .and.
     &          (ir .ne. iru .or. ia .ne. iau) ) go to 756
              write (nlist,740) (ltyp(icol),lbla(icol), icol=1,ncol)
  740         format(' ctr, ao',16(3x,a3,i1))
              ict = 0
              do 752 ic = 1,icu
                do 748 it = itl,itu
                  ict = ict+1
                  write (nlist,744) ic, (lmnv(i,it), i=1,3),
     &              (c(ict,icsr+icol,igcs), icol=1,ncol)
  744             format(i3,', ',3i1,16f7.3)
  748           continue
  752         continue
              icsr = icsr+ncol
              ncol = 0
  756       continue
  760     continue
          do 792 ircr = 1,ircru
            ics = 0
            do 788 ir = 1,iru
              isfr = isfr+1
              lai = la(ir,iaords)
              iso = nsopr(lai)+lb(isfr)
              ms(iso) = is
c
c             # assign the mnl(*) values to the so functions.
c             # 1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g,...
c             #
c             # compute the L^2-operator expectation value and
c             # determine if it is an eigenvalue.
c
              lmn = lmnp1(icons) - 1
c             numxyz = ((lmn + 1) * (lmn + 2)) / 2
c             # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:scr(1:xyzdim)
              call l2opxv( lmn, c(1,ics+1,igcs), itl, itu, lmnv,
     &         ipt(55)-ipt(2), a(ipt(2)), eval, leig )
              if ( leig .lt. 0 ) then
c
c               # not an eigenvalue. print a warning.
c
                write(nlist,768) leig, iso, ityp(lai), eval
  768           format(' syminp() warning: from l2opxv(), ierr=',i3,
     &            ' s.o. basis function', i3, a, ' is not an ',
     &            'L^2 eigenvector. <v|L^2|v>/<v|v> =', f12.8 )
c
c               # just assign the maximum value.
                leig = lmn
              endif
              mnl(iso) = ( lmnp1(icons)**2 + 2 * leig + 5 ) / 4
              iau = nd(lai)
              do 784 ia = 1,iau
                ics = ics+1
c
                if ( inam(nu) .ge. 5 .and. leig .eq. 1 ) then
c
c                 # determine angular momentum symmetries.
c
                  do 780 ixyz = 1, 3
                    i2 = mod(ixyz,  3) + 1
                    i3 = mod(ixyz+1,3) + 1
                    sum = a0
                    ict = 0
                    do 776 ic = 1, icu
                      do 772 it = itl,itu
                        ict = ict + 1
                        if( mod(lmnv(ixyz,it),2) .eq. 1 .and.
     &                      mod(lmnv(i2,  it),2) .eq. 0 .and.
     &                      mod(lmnv(i3,  it),2) .eq. 0 ) then
                          sum = sum + c(ict,ics,igcs)
                        endif
  772                 continue
  776               continue
                    if (  (sum .eq. a0)
     &               .or. (lai .eq. ixyzir(ixyz)) ) go to 780
                    if ( ixyzir(ixyz) .ne. 0 ) then
                      write (nlist,*) 'axis choice is unsuitable for'
     &                 //' spin-orbit integrals'
                      call bummer('syminp: ixyzir(ixyz)=',
     &                 ixyzir(ixyz),faterr)
                    endif
                    ixyzir(ixyz) = lai
  780             continue
                endif
  784         continue
  788       continue
  792     continue
  796   continue
c
        icrs = mcrs(is)
        if ( icrs .eq. 0 ) go to 840
        lcru = lcr(icrs)
        if ( lcru .lt. 0 ) go to 816
        kcrl = nkcrl(1,icrs)
        kcru = nkcru(1,icrs)
        write (nlist,804) mtype(is),lblsh(lcru+1),
c    &   (ncr(kcr),kcr=kcrl,kcru)
     &   (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
  804   format(//31x,a3,' core potential'//35x,a1,' potential'/
c    &    7h0powers,6x,7(i8,7x))
     &    19x,'powers',6x,'exponentials',4x,'coefficients'/
     &    (i22,7x,2(1pg16.7)))
c       write (nlist,828) (zcr(kcr), kcr=kcrl,kcru)
c       write (nlist,832) (ccr(kcr), kcr=kcrl,kcru)
        if ( lcru .eq. 0 ) go to 816
        do 812 lp1 = 1,lcru
          kcrl = nkcrl(lp1+1,icrs)
          kcru = nkcru(lp1+1,icrs)
          write (nlist,808) lblsh(lp1), lblsh(lcru+1),
c    &     (ncr(kcr),kcr=kcrl,kcru)
     &     (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
c 808     format(/33x,a1,' - ',a1,' potential'/7h0powers,6x,7(i8,7x))
  808     format(/33x,a1,' - ',a1,' potential'/
     &      19x,'powers',6x,'exponentials',4x,'coefficients'/
     &      (i22,7x,2(1pg16.7)))
c         write (nlist,828) (zcr(kcr), kcr=kcrl,kcru)
c         write (nlist,832) (ccr(kcr), kcr=kcrl,kcru)
  812   continue
c
  816   continue
        llsu = lls(icrs)
        if ( llsu .le. 0 ) go to 840
        write (nlist,820) mtype(is)
  820   format(//28x,a3,' spin-orbit potential')
        do 836 l = 1,llsu
          kcrl = nklsl(l,icrs)
          kcru = nklsu(l,icrs)
c         write (nlist,824) lblsh(l+1), (ncr(kcr), kcr=kcrl,kcru)
          write (nlist,824) lblsh(l+1),
     &     (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
c 824     format(/10x,a1,10h potential/7h0powers,6x,7(i8,7x))
  824     format(/35x,a1,' potential'/
     &      19x,'powers',6x,'exponentials',4x,'coefficients'/
     &      (i22,7x,2(1pg16.7)))
c         write (nlist,828) (zcr(kcr), kcr=kcrl,kcru)
c 828     format(13h exponentials,7(1pg15.7))
c         write (nlist,832) (ccr(kcr), kcr=kcrl,kcru)
c 832     format(13h coefficients,7(1pg15.7))
  836   continue
  840 continue
c
      if ( inam(nu) .eq. 5 ) then
        do 860 ixyz = 1, 3
          i2 = mod(ixyz,  3) + 1
          i3 = mod(ixyz+1,3) + 1
          do 856 ist = 1, nst
            if ( idp(ist,ixyzir(i2),ixyzir(i3)) .ne. 0 ) then
              lxyzir(ixyz) = ist
              go to 860
            endif
  856     continue
  860   continue
        write (nlist,864) (ityp(lxyzir(ixyz)), ixyz=1,3)
  864 format(/'lx: ',a3,10x,'ly: ',a3,10x,'lz: ',a3)
      endif
c
c     # set up the 1-e record parameters.
      ibvtyp = 0
      if ( l1rec .eq. 0 ) then
c       # use the default value for this program.
        l1recx = lrcmxp
      else
c       # use the input value. [(-1) invokes the SIFS default.]
        l1recx = l1rec
      endif
c
      call sifcfg( 1, l1recx, nbft, ibvtyp, ifmt1, l1rec, n1max, ierr )
c
      if ( ierr .ne. 0 ) then
        call bummer('syminp: from 1-e sifcfg(), ierr=',ierr,faterr)
      endif
c
c     # set up the 2-e record parameters.
      ibvtyp = 1
      if ( l2rec .eq. 0 ) then
c       # use the default value for this program.
        l2recx = lrcmxp
      else
c       # use the input value. [(-1) envokes the SIFS default.]
        l2recx = l2rec
      endif
c
      call sifcfg( 2, l2recx, isfr, ibvtyp, ifmt2, l2rec, n2max, ierr )
c
      if ( ierr .ne. 0 ) then
        call bummer('syminp: from 2-e sifcfg(), ierr=',ierr,faterr)
      endif
c
c     # set up the info(*) array for the output file.
      info(1) = fsplit
      info(2) = l1rec
      info(3) = n1max
      info(4) = l2rec
      info(5) = n2max
c
c     # open the aoints file.
c
      open(unit=aoints,status='unknown',form='unformatted',
     & file=fnames(2))
c
      rewind aoints
c
c     # write out the header information
c
c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:map(1:2,1:nbft)
c     # push map(*) allocation onto the stack.
      call h2opsh( forbyt(2*nbft) )
c
      call wthead( nst,    ns,     isfr,   ntitle, title,  repnuc, ityp,
     &   nso,      mtype,  ms,     mnl,    info,   aoints, nlist,
     &   a(ipt(2)) )
c
c     # pop back to the initial level.
      call h2opop
c
c     # multiply normalization constants into contraction coefficients.
c
      do 876 icons = 1,ncons
        iconu = ncon(icons)
        ircru = nrcr(icons)
        do 872 icon = 1,iconu
          do 868 ircr = 1,ircru
            t = a4 * zet(icon,icons)
            eta(ircr,icon,icons) =
     &       sqrt( t**lmnp1(icons) / a2 * sqrt(t / a2) )
     &       * eta(ircr,icon,icons)
  868     continue
  872   continue
  876 continue
      write (nlist,880) isfr, (ityp(ist), nso(ist), ist=1,nst)
  880 format(//i9,' symmetry orbitals,',3x,4(3x,a3,':',i4),
     &                               (/31x,4(3x,a3,':',i4)))
c
      return
      end
*deck lmnvgn
      subroutine lmnvgn(lmn1u,lmnv)
c     # generate cartesian gaussian exponent array.
c     # lmnv(*,*) = exponents of the cartesian gaussian basis functions.
c     #             s  p  d   f   g   h   i
c     #       lmn = 0  1  2   3   4   5   6
c     #    numxyz = 1, 3, 6, 10, 15  21  28      =((lmn+1)*(lmn+2))/2
      implicit logical (a-z)
      integer i, j, k, lmn, lmn1u, lmnv(3,*), ndx
c
      ndx = 0
      do 40 lmn = 0,lmn1u-1
        do 30 i = lmn,(lmn+2)/3,-1
          do 20 j = min(i,lmn-i),(lmn-i+1)/2,-1
            k = lmn - i -j
            if( i.eq.j .and. j.eq.k ) then
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = i
              lmnv(3,ndx) = i
            elseif( i.eq.j .and.j.gt.k ) then
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = i
              lmnv(3,ndx) = k
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = k
              lmnv(3,ndx) = i
              ndx = ndx + 1
              lmnv(1,ndx) = k
              lmnv(2,ndx) = i
              lmnv(3,ndx) = i
            elseif( i.gt.j .and. j.eq.k ) then
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = j
              lmnv(3,ndx) = j
              ndx = ndx + 1
              lmnv(1,ndx) = j
              lmnv(2,ndx) = i
              lmnv(3,ndx) = j
              ndx = ndx + 1
              lmnv(1,ndx) = j
              lmnv(2,ndx) = j
              lmnv(3,ndx) = i
            elseif( i.gt.j .and. j.gt.k ) then
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = j
              lmnv(3,ndx) = k
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = k
              lmnv(3,ndx) = j
              ndx = ndx + 1
              lmnv(1,ndx) = j
              lmnv(2,ndx) = i
              lmnv(3,ndx) = k
              ndx = ndx + 1
              lmnv(1,ndx) = k
              lmnv(2,ndx) = i
              lmnv(3,ndx) = j
              ndx = ndx + 1
              lmnv(1,ndx) = j
              lmnv(2,ndx) = k
              lmnv(3,ndx) = i
              ndx = ndx + 1
              lmnv(1,ndx) = k
              lmnv(2,ndx) = j
              lmnv(3,ndx) = i
            endif
   20     continue
   30   continue
   40 continue
c
c  *********************
c     new ordering of 
c      d exponents
c  *********************
      lmnv(1,5) = 2     
      lmnv(2,5) = 0
      lmnv(3,5) = 0
      lmnv(1,6) = 1
      lmnv(2,6) = 1
      lmnv(3,6) = 0
      lmnv(1,7) = 1
      lmnv(2,7) = 0
      lmnv(3,7) = 1
      lmnv(1,8) = 0
      lmnv(2,8) = 2
      lmnv(3,8) = 0
      lmnv(1,9) = 0
      lmnv(2,9) = 1
      lmnv(3,9) = 1
      lmnv(1,10) = 0
      lmnv(2,10) = 0
      lmnv(3,10) = 2
c  *********************
c     new ordering of
c      f exponents  
c  *********************
      lmnv(1,11) = 3
      lmnv(2,11) = 0
      lmnv(3,11) = 0
      lmnv(1,12) = 2
      lmnv(2,12) = 1
      lmnv(3,12) = 0
      lmnv(1,13) = 2
      lmnv(2,13) = 0
      lmnv(3,13) = 1
      lmnv(1,14) = 1
      lmnv(2,14) = 2
      lmnv(3,14) = 0
      lmnv(1,15) = 1
      lmnv(2,15) = 1
      lmnv(3,15) = 1
      lmnv(1,16) = 1
      lmnv(2,16) = 0
      lmnv(3,16) = 2
      lmnv(1,17) = 0
      lmnv(2,17) = 3
      lmnv(3,17) = 0
      lmnv(1,18) = 0
      lmnv(2,18) = 2
      lmnv(3,18) = 1
      lmnv(1,19) = 0
      lmnv(2,19) = 1
      lmnv(3,19) = 2
      lmnv(1,20) = 0
      lmnv(2,20) = 0
      lmnv(3,20) = 3
c  *********************
c     new ordering of
c      g exponents  
c  *********************
      lmnv(1,21) = 4
      lmnv(2,21) = 0
      lmnv(3,21) = 0
      lmnv(1,22) = 3
      lmnv(2,22) = 1
      lmnv(3,22) = 0
      lmnv(1,23) = 3
      lmnv(2,23) = 0
      lmnv(3,23) = 1
      lmnv(1,24) = 2
      lmnv(2,24) = 2
      lmnv(3,24) = 0
      lmnv(1,25) = 2
      lmnv(2,25) = 1
      lmnv(3,25) = 1
      lmnv(1,26) = 2
      lmnv(2,26) = 0
      lmnv(3,26) = 2
      lmnv(1,27) = 1
      lmnv(2,27) = 3
      lmnv(3,27) = 0
      lmnv(1,28) = 1
      lmnv(2,28) = 2
      lmnv(3,28) = 1
      lmnv(1,29) = 1
      lmnv(2,29) = 1
      lmnv(3,29) = 2
      lmnv(1,30) = 1
      lmnv(2,30) = 0
      lmnv(3,30) = 3
      lmnv(1,31) = 0
      lmnv(2,31) = 4
      lmnv(3,31) = 0
      lmnv(1,32) = 0
      lmnv(2,32) = 3
      lmnv(3,32) = 1
      lmnv(1,33) = 0
      lmnv(2,33) = 2
      lmnv(3,33) = 2
      lmnv(1,34) = 0
      lmnv(2,34) = 1
      lmnv(3,34) = 3
      lmnv(1,35) = 0
      lmnv(2,35) = 0
      lmnv(3,35) = 4
c
      return
      end
*deck l2opxv
      subroutine l2opxv( lmn, v, itl, itu, lmnv, lenscr, scr, eval,
     & leig )
c
c  compute the matrix vector product, (L^2) * v, in an unnormalized
c  cartesian basis, and determine if v(*) is an eigenvector of the
c  total angular momentum operator.
c
c  input:
c  lmn = (l + m + n) where l, m, and n are the exponents of x, y, and z
c        respectively.  lmn=0 for cartesian "s" functions, 1 for
c        cartesian "p" functions, 2 for cartesian "d" functions, etc.
c  v(1:xyzdim) = input vector.  the elements are coefficients of the
c              cartesian basis functions (x**l * y**m * z**n).
c              where xyzdim = ((lmn + 1) * (lmn + 2))/2
c  lmnv(1:3,1:lmnvmx) = cartesian exponents of each basis function.
c  lenscr = scratch vector length.
c  scr(1:lenscr) = scratch work vector.  This must be at least as large
c                  as the cartesian subspace dimension.
c                  lenscr >= xyzdim
c
c  output:
c  eval = vector expectation value = <v| L^2 |v> / <v|v>
c  leig = eigenvector return code.
c       = -4 for vector norm error.
c       = -3 for lmnv(*,*) inconsistency.
c       = -2 for lenscr error.
c       = -1 if v(*) is not an eigenvector of L^2
c       = principal quantum number if v(*) is an eigenvector.
c         0 for s-type vectors, 1 for p-type vectors, etc.
c         the eigenvalue is given by (leig*(leig+1)) to within numerical
c         precision.
c  28-oct-96 account for nonorthogonality of cartesian gaussians (rmp)
c  11-jun-91 written by ron shepard with suggestions from r. m. pitzer.
c
      implicit logical(a-z)
c
c     # dummy:
c     integer lmn, lenv, lenscr, leig
      integer itl, itu, lmn, lenscr, leig
      integer lmnv(3,*)
      real*8 eval
      real*8 v(*), scr(*)
c
c     # local:
      integer l, m, n, i, j, it, jt, xyzdim, ll, mm, nn, p, q
      real*8 vi, rnorm, vnorm2
c
c     # small is used to determine if a vector is close enough to be
c     #       called an eigenvector.  this should be about 1000x the
c     #       machine epsilon to avoid sqrt() precision problems.
c
      real*8    zero,     one,     two,     fourth,        small
      parameter(zero=0d0, one=1d0, two=2d0, fourth=0.25d0, small=1d-10)
c
c     # bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
      integer  iodfac
      external iodfac
c
c     # the [l,m,n] basis function is assigned the internal index:
c     #  (lx*(lx+1))/2 + n + 1   with lx=(lmn-l)=(m+n)
c     # this convention is inconsistent with the rest of this program,
c     # but it is used anyway since it is convenient and since the
c     # computed results are local.
c
      integer llp1
      llp1(l) = (l * (l + 1)) / 2 + 1
c
      leig = -1
      eval = -one
c
c     # check to make sure scr(*) is large enough.
      xyzdim = llp1( lmn) + lmn
      if ( xyzdim .gt. lenscr ) then
         call bummer('l2opxv: need larger scr(*), (xyzdim-lenscr)=',
     &    (xyzdim-lenscr), wrnerr )
         leig = -2
         return
      endif
c
c     # initialize the scratch vector.
      call wzero(xyzdim,scr,1)
c
c     # compute the matrix-vector product and the vector norm.
c     # matrix elements of lop(*,*) are computed from the operator
c     # definition in cartesian space.
c
      vnorm2 = zero
      do 30 it = itl, itu
         i = it - (itl - 1)
         l = lmnv(1,it)
         m = lmnv(2,it)
         n = lmnv(3,it)
         vi = v(i)
         do 20 jt = itl, it-1
            j = jt - (itl - 1)
            vnorm2 = vnorm2 + (two * vi) * v(j) *
     &       iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
   20    continue
         vnorm2 = vnorm2 + vi * vi * iodfac(l,l,m,m,n,n)
c
         if ( (l+m+n) .ne. lmn ) then
c           # inconsistent exponents in the basis function.
            call bummer('l2opxv: (l+m+n-lmn)=', (l+m+n-lmn), wrnerr )
            leig = -3
            return
         endif
c
c        # L^2 is sparse in this representation.  each vi contributes
c        # to, at most, only 7 elements in the matrix-vector product.
c
         if ( l .ge. 2 ) then
c           # include -l*(l-1) * ( [l-2,m,n+2] + [l-2,m+2,n] ) terms.
            ll = l * (l - 1)
            p = llp1( lmn - l + 2 ) + n
            scr(p) = scr(p) - vi * (ll)
            p = p + 2
            scr(p) = scr(p) - vi * (ll)
         endif
c
         if ( m .ge. 2 ) then
c           # include -m*(m-1) * ( [l,m-2,n+2] + [l+2,m-2,n] ) terms.
            mm = m * (m - 1)
            p = llp1( lmn - l ) + n + 2
            scr(p) = scr(p) - vi * (mm)
            p = llp1( lmn - l - 2 ) + n
            scr(p) = scr(p) - vi * (mm)
         endif
c
         if ( n .ge. 2 ) then
c           # include -n*(n-1) * ( [l,m+2,n-2] + [l+2,m,n-2] ) terms.
            nn = n * (n - 1)
            p = llp1( lmn - l ) + n - 2
            scr(p) = scr(p) - vi * (nn)
            p = llp1( lmn - l - 2 ) + n - 2
            scr(p) = scr(p) - vi * (nn)
         endif
c
c        # include the 2*(l*m+l*n+m*n+l+m+n)*[l,m,n] diagonal term.
c
         p = llp1( lmn - l ) + n
         scr(p) = scr(p) + vi * (2 * (l * (m + n) + m * n + l + m + n))
   30 continue
c
      if ( vnorm2 .le. small ) then
         call bummer('l2mxv: small vector norm, lnm=', lmn, wrnerr)
         leig = -4
         return
      endif
c
c     # compute the expectation value.
c
      eval  = zero
      do 50 it = itl, itu
         i = it - (itl - 1)
         l = lmnv(1,it)
         m = lmnv(2,it)
         n = lmnv(3,it)
         p = llp1( lmn - l ) + n
         do 40 jt = itl, itu
            j = jt - (itl - 1)
            eval = eval + v(j) * scr(p) *
     &                  iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
   40    continue
   50 continue
      eval = eval / vnorm2
c
c     # compute the residual norm.
c     #
      rnorm = zero
c
      do 70 it = itl, itu
         i = it - (itl - 1)
         l = lmnv(1,it)
         m = lmnv(2,it)
         n = lmnv(3,it)
         p = llp1( lmn - l ) + n
         do 60 jt = itl, it-1
            j = jt - (itl - 1)
            q = llp1( lmn - lmnv(1,jt) ) + lmnv(3,jt)
            rnorm = rnorm + two * ( scr(p) - eval * v(i) ) *
     &                            ( scr(q) - eval * v(j) ) *
     &                iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
   60    continue
        rnorm = rnorm + (scr(p) - eval * v(i))**2 * iodfac(l,l,m,m,n,n)
   70 continue
c
c     # normalize w.r.t. |v|=1.
      rnorm = rnorm / vnorm2
c
      rnorm = sqrt( rnorm )
      if ( rnorm .gt. small ) then
c
c        # v(*) is not an eigenvector.
c
         leig = -1
      else
c
c        # v(*) is an eigenvector.
c
c        # determine leig such that eval = leig*(leig+1)
c
c        # the following assignment should truncate to the
c        # next smaller integer value...
         leig = ( sqrt( eval + fourth ) )
      endif
c
      return
      end
*deck socfpd
      subroutine socfpd( c,           cx,         icxast, icxst,
     &   icxsv1, icxsv2, idp,         il,         iprst,  la,
     &   lb,     maords, mau,         mcons,      mgcs,   mics,
     &   mjcs,   nblpr,  nc,          nd,         nf,     nir,
     &   npair,  nprir,  nrcr,        nsopr,      nt,     ntl,
     &   ntu  )
c
c  01-dec-90 unnecessary statement labels removed to allow
c            optimization. il(*) encoded with arithmetic. -rls
c
      implicit logical(a-z)
c
c     # il(*) packing factor. see also oneint() and wtint2().
      integer    ilfact
      parameter( ilfact=1024 )
c
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)
c
      integer        ipq
      common /dim21/ ipq(256)
c
c     # dummy:
      integer
     & icxast(*),        icxst(2,*),   icxsv1(*),    icxsv2(mgcsu,*),
     & idp(mstu,mstu,*), il(nnbft),    iprst(*),     la(mru,*),
     & lb(*),            maords(*),    mau(*),       mcons(*),
     & mgcs(*),          mics(*),      mjcs(*),      nblpr(*),
     & nc(*),            nd(*),        nf(*),        nir(*),
     & npair(2,*),       nrcr(*),      nsopr(*),     nprir(2,mstu,*),
     & nt(*),            ntl(*),       ntu(*)
      real*8 c(mctu,mcsu,*), cx(*)
c
c     # local:
      integer mcxu1, icx, ipair, ijsf, isf, isfr, igcs, jgcs, is,
     & ifu, icu, isfis, isfris, jsf, jsfr, js, jsfjs, jsfrjs, jfu,
     & jcu, if, isfrif, iaords, iru, ircru, itl, itu, jf, jaords,
     & jtl, jsfrjf, jtu, jru, ist, jrcru, ircr, isfrib, jrcr,
     & jsfrjb, iesfb, npr, ipr, ics, ir, lai, iau, jcs, jr, laj,
     & isv1, icxu, icxa, isv2, ict, ic, ictic, jct, jc, jctjc, it,
     & jt, ia, icxl, mnsf, npru
      logical ea, easfb, es, esf, esfb, esfbc, msfbct
c
      real*8     a0
      parameter (a0=0.0d0)
c
c     # bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
c     # check the il(*) packing factor.
      if ( nbft .ge. ilfact ) then
        call bummer('socfpd: larger ilfact required, nbft=',
     &   nbft, faterr )
      endif
c
      ea = inam(nu).eq.5
c
c     # mcxu1 is the original mcxu offset.
      mcxu1 = mcxu
      icx   = 0
      ipair = 0
c
c     # ijsf = one dimensionalization of is, js, if, jf
c     # isfr = one dimensionalization of is, if, ir
      ijsf = 0
      isf  = 0
      isfr = 0
      do 12 igcs = 1, ngcs
        icxsv1(igcs) = 0
        do 10 jgcs = 1, ngcs
          icxsv2(jgcs,igcs) = 0
10      continue
12    continue
      do 340 is = 1, ns
        ifu = nf(is)
        if(ifu.eq.0) go to 340
        icu = nc(is)
        isfis = isf
        isfris = isfr
        jsf = 0
        jsfr = 0
        do 330 js = 1,is
          jsfjs = jsf
          jsfrjs = jsfr
          es = is.eq.js
          if(.not.es) then
            jfu = nf(js)
            if(jfu.eq.0) go to 330
            jcu = nc(js)
          endif
          isf = isfis
          isfr = isfris
          do 324 if = 1,ifu
            isf = isf+1
            isfrif = isfr
            igcs = mgcs(isf)
            iaords = maords(igcs)
            iru = nir(iaords)
            ircru = nrcr(mcons(isf))
            itl = ntl(isf)
            itu = ntu(isf)
            jsf = jsfjs
            jsfr = jsfrjs
            if(es) jfu = if
            do 322 jf = 1,jfu
              jsf = jsf+1
              jgcs = mgcs(jsf)
              jaords = maords(jgcs)
              jtl = ntl(jsf)
              jsfrjf = jsfr
              ijsf = ijsf+1
              esf = isf.eq.jsf
              if(esf) then
                if(ircru.ge.2) then
                  jtu = ntu(jsf)
                  jru = nir(jaords)
                else
                  npair(2,ijsf) = 0
                  icxst(2,ijsf) = 0
                  do 20 ist = 1,nst
                    nprir(2,ist,ijsf) = 0
20                continue
                endif
              else
                jrcru = nrcr(mcons(jsf))
                jtu   = ntu(jsf)
                jru   = nir(jaords)
              endif
              iprst(ijsf) = ipair
              isfr = isfrif
              do 320 ircr = 1,ircru
                isfrib = isfr
                jsfr = jsfrjf
                if(esf) jrcru = ircr
                do 310 jrcr = 1,jrcru
                  jsfrjb = jsfr
                  esfb = esf.and.ircr.eq.jrcr
                  easfb = ea.and.esfb
                  if(esfb) then
                    iesfb = 1
                    npru = (iru * (iru + 1)) / 2
                  else
                    iesfb = 2
                    npru = iru * jru
                  endif
c
                  if ( npru .gt. mnru ) then
                    call bummer('change mnrup (one place) to ',
     &               npru,faterr)
                  endif
c
                  npr = 0
                  do 75 ist = 1, nst
                    ipr = npr
                    ics = 0
                    isfr = isfrib
                    do 70 ir = 1, iru
                      isfr = isfr + 1
                      lai = la(ir,iaords)
                      iau = nd(lai)
                      if(esfb) jru = ir
                      jsfr = jsfrjb
                      jcs = 0
                      do 65 jr = 1, jru
                        jsfr = jsfr + 1
                        laj = la(jr,jaords)
                        if ( idp(lai,laj,ist) .ne. 0 ) then
                          ipair = ipair + 1
                          il(ipair) = (nsopr(lai) + lb(isfr))
     &                     * ilfact + (nsopr(laj) + lb(jsfr))
c
                          npr = npr + 1
                          mau(npr)  = iau
                          mics(npr) = ics
                          mjcs(npr) = jcs
                        endif
                        jcs = jcs + nd(laj)
65                    continue
                      ics = ics + iau
70                  continue
                    nprir(iesfb,ist,ijsf) = npr - ipr
75                continue
                  npair(iesfb,ijsf) = npr
c
                  if ( npr .eq. 0 ) go to 310
c
                  if ( esfb ) then
                    isv1 = icxsv1(igcs)
                    if ( isv1 .ne. 0 ) then
                      icxst(1,ijsf) = icxst(1,isv1)
                      if(ea) icxast(ijsf) = icxast(isv1)
                      go to 310
                    endif
                    icxsv1(igcs) = ijsf
                    icxst(1,ijsf) = icx
                    icxu = icx+(npr*ipq(icu*nt(isf)+1))
                    if ( ea ) then
                      mcxu = mcxu-(npr*ipq(icu*nt(isf)+1))
                      icxa = mcxu
                      icxast(ijsf) = icxa
                    endif
                  else
                    isv2 = icxsv2(jgcs,igcs)
                    if ( isv2 .ne. 0 ) then
                      icxst(2,ijsf) = icxst(2,isv2)
                      go to 310
                    endif
                    icxsv2(jgcs,igcs) = ijsf
                    icxst(2,ijsf) = icx
                    icxu = icx+npr*icu*nt(isf)*jcu*nt(jsf)
                  endif
                  if ( icxu .gt. mcxu ) then
                    call bummer('socfpd: (icxu-mcxu)=',
     &               (icxu-mcxu), faterr )
                  endif
c
c                 # set the high-water mark.
c                 # account for both low-address allocations,
c                 # from (1 : icxu), and high-address allocations,
c                 # from (mcxu+1 : mcxu1).
                  call h2oset( (icxu + mcxu1 - mcxu) )
c
c                 # compute new cx block
c
                  ict = 0
                  do 298 ic = 1, icu
                    ictic = ict
                    if(esfb) jcu = ic
                    jct = 0
                    do 296 jc = 1, jcu
                      jctjc = jct
                      esfbc = esfb .and. ic .eq. jc
                      ict = ictic
                      do 294 it = itl, itu
                        ict = ict + 1
                        if(esfbc) jtu = it
                        jct = jctjc
                        do 292 jt = jtl, jtu
                          jct = jct + 1
                          msfbct =  .not. esfb .or.
     &                     (ic.eq.jc .and. it.eq.jt)
                          do 290 ipr = 1, npr
                            icx = icx + 1
                            cx(icx) = a0
                            if(easfb) then
                              icxa = icxa + 1
                              cx(icxa) = a0
                            endif
                            iau = mau(ipr)
                            ics = mics(ipr)
                            jcs = mjcs(ipr)
                            do 280 ia = 1, iau
                              ics = ics + 1
                              jcs = jcs + 1
                              cx(icx) = cx(icx)
     &                         + (c(ict,ics,igcs) * c(jct,jcs,jgcs))
                              if(easfb) cx(icxa) = cx(icxa)
     &                         + (c(ict,ics,igcs) * c(jct,jcs,jgcs))
                              if(msfbct) go to 280
                              cx(icx) = cx(icx)
     &                         + (c(jct,ics,igcs) * c(ict,jcs,jgcs))
                              if(easfb) cx(icxa) = cx(icxa)
     &                         - (c(jct,ics,igcs) * c(ict,jcs,jgcs))
280                         continue
290                       continue
292                     continue
294                   continue
296                 continue
298               continue
310             continue
320           continue
322         continue
324       continue
330     continue
340   continue
c
      if ( ea ) then
c
c       # move antisymmetric coefficient products down.
c
c       #       offsets: 0           mcxu1  mblu
c       #                |            |     |
c       # before: a(*)-->.............c(*)
c       #
c       #                0      mcxu mcxu1  mblu
c       #                |       |    |     |
c       # now:    a(*)-->cx(*)...x(*),c(*)
c       #
c       #                0   mcxu2 mcxu  (mblu+mcxu)
c       #                |     |   |        |
c       # after:  a(*)-->cx(*),x(*),........
c
c       # on return, mcxu2 is the offset for the s.o. coefficients.
        mcxu2 = icx
c
        icxl  = icx + 1
        icxu  = mcxu - icx
        mcxu  = mcxu1 - icxu
        do 370 icx = icxl, mcxu
          cx(icx) = cx(icx+icxu)
370     continue
c       # make pointers relative to mcxu2 instead of 0.
        do 380 mnsf = 1, ijsf
          icxast(mnsf) = icxast(mnsf) - icxu
380     continue
c       # on return, mcxu and mblu divide the workspace as:
c       # a(*)-->[cx//x](1:mcxu),w(1:mblu)
        mblu = mblu - mcxu
        write (nlist,6010) mcxu, mcxu2, mblu
      else
c       # on return, mcxu and mblu divide the workspace as:
c       # a(*)-->cx(1:mcxu),w(1:mblu)
        mcxu = icx
        mblu = mblu - mcxu
        write (nlist,6020) mcxu, mblu
      endif
c
      return
6010  format(/' socfpd: mcxu=',i9,' mcxu2=',i9,' left=',i9)
6020  format(/' socfpd: mcxu=',i9,' left=',i9)
      end
*deck hermit
      subroutine hermit(nn,x,a,eps)
c                  calculates the zeros  x(i)  of the nn-th order
c                hermite polynomial.  the largest zero will be
c                stored in x(1).  also calculates the corresponding
c                coefficients  a(i)  of the nn-th order gauss-hermite
c                quadrature formula of degree 2*nn-1.  the factor of
c                sqrt(pi) has been removed from the a(i).
c  a. h. stroud & d. secrest, gaussian quadrature formulas,
c  prentice-hall, 1966
c
      implicit real*8 (a-h,o-z)
      parameter (sixth = 1.0d0/6.0d0)
      dimension x(*), a(*)
c
      fn = (nn)
      n1 = nn - 1
      n2 = (nn+1)/2
      cc = 1.0d0
      s  = 0.0d0
      do 15 i=1,n1
        s  = s + 0.5d0
        cc = s*cc
   15 continue
      s  = (2.0d0*fn+1.0d0)**sixth
      do 30 i=1,n2
        if( i.eq.1 ) then
c         # largest zero
          xt = s**3 - 1.85575d0/s
        elseif( i.eq.2 ) then
c         # second zero
          xt = xt - 1.14d0*fn**0.426d0/xt
        elseif( i.eq.3 ) then
c         # third zero
          xt = 1.86d0*xt - 0.86d0*x(1)
        elseif( i.eq.4 ) then
c         # fourth zero
          xt = 1.91d0*xt - 0.91d0*x(2)
        else
c         # all other zeros
          xt = 2.0d0*xt - x(i-2)
        endif
c
        call hroot(xt,nn,dpn,pn1,eps)
        x(i) = xt
        a(i) = cc/dpn/pn1
c       write (6,'(2i4,2d25.17)') nn, i, xt, a(i)
        ni = nn-i+1
        x(ni) = -xt
        a(ni) = a(i)
   30 continue
      return
      end
*deck hroot
      subroutine hroot(x,nn,dpn,pn1,eps)
c                  improves the approximate root  x
c                in addition we also obtain
c                    dpn = derivative of h(n) at x
c                    pn1 = value of h(n-1) at x
c
      implicit real*8 (a-h,o-z)
c
c     # iter = 5 sufficient for 8-byte accuracy up to nn = 7
      do 14 iter=1,10
        call hrecur_argos(p,dp,pn1,x,nn)
        d  = p/dp
        x  = x - d
        if( abs(d).le.eps ) go to 16
   14 continue
   16 dpn = dp
      return
      end
*deck hrecur_argos
      subroutine hrecur_argos(pn,dpn,pn1,x,nn)
c
      implicit real*8 (a-h,o-z)
c
      p1 = 1.0d0
      p  = x
      dp1 = 0.0d0
      dp  = 1.0d0
      do 20 j=2,nn
        fj = (j)
        fj2 = 0.5d0*(fj-1.0d0)
        q  = x*p - fj2*p1
        dq = x*dp + p - fj2*dp1
        p1 = p
        p  = q
        dp1 = dp
        dp  = dq
   20 continue
      pn  = p
      dpn = dp
      pn1 = p1
      return
      end
*deck oneint_argos
      subroutine oneint_argos( a,          ccr,        chg,        cx,
     &       eta,        h2,         ica,
     &       icb,        icxast,     icxst,      il,         ilxyz,
     &       iprst,      ipt,        lcr,        lls,        lmnp1,
     &       lproju,     mcons,      mcrs,       nc,         ncon,
     &       ncr,        nf,         nfct,       ngw,        nkcrl,
     &       nkcru,      nklsl,      nklsu,      nopir,      npair,
     &       nprir,      nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zcr,        zet,
     &       array,      iarray,     buffer,     values,     labels )
c
c  04-dec-90 array(*,*) added for spin-orbit integrals. -rls
c  01-dec-90 SIFS version. il(*) decoded with arithmetic. -rls
c
      implicit logical(a-z)
c
      integer    ilfact
      parameter( ilfact=1024 )
c
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
c
      real*8         cutoff, tol
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)
c
      real*8
     & dxij,   dyij,   dzij,   fnfct,  rr,     xij,    xijm,
     & yij,    yijm,   zij,    zijm
      integer  ibl1,   ibl2,   icxi1,  icxi2,  ij,     ijsf,
     & ic,     icons,  igu,
     & ircru,  is,     isf,    itl,    itu,    jc,     jcons,
     & jgu,    jrcru,  js,     jsf,    jtl,    jtu,    lit,
     & ljt,    nblt1,  nc2,    nc1,    nop,    ntij1,  ntij2
      logical
     & esf,    esfc,   igueq1, jgueq1
      common /one/     dxij,   dyij,   dzij,   fnfct,  rr,
     & xij,    xijm,   yij,    yijm,   zij,    zijm,   ibl1,
     & ibl2,   icxi1,  icxi2,  ij,     ijsf,   ic,     icons,
     & igu,    ircru,  is,     isf,    itl,    itu,    jc,
     & jcons,  jgu,    jrcru,  js,     jsf,    jtl,    jtu,
     & lit,    ljt,    nblt1,  nc2,    nc1,    nop,    ntij1,
     & ntij2,  esf,    esfc,   igueq1, jgueq1
c
      integer        ipq
      common /dim21/ ipq(256)
c
c     # /bufout/ holds some output integral file parameters.
      integer
     & info,     ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & info(10), ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
c
      integer       l1rec
      equivalence ( l1rec, info(2) )
      integer       n1max
      equivalence ( n1max, info(3) )
c
      integer    nipv
      parameter( nipv=2 )
c
c     # dummy:
      integer           ica(mcu,msu,*), icb(4,24,*), icxast(*),
     & icxst(2,*),      il(nnbft),      ilxyz(3,*),  iprst(*),
     & ipt(*),          lcr(*),         lls(*),      lmnp1(*),
     & lproju,          mcons(*),       mcrs(*),     nc(*),
     & ncon(*),         ncr(*),         nf(*),       nfct(*),
     & ngw(*),          nkcrl(6,*),     nkcru(6,*),  nklsl(4,*),
     & nklsu(4,*),      nopir(*),       npair(2,*),  nprir(2,mstu,*),
     & nrcr(*),         nt(*),          ntl(*),      ntu(*),
     & iarray(mpru,*),  labels(nipv,*)
      real*8
     & a(*),         ccr(*),       chg(*),        cx(*),
     & h2(*),
     & x(mcu,*),     y(mcu,*),     z(mcu,*),      eta(mrcru,mconu,*),
     & zcr(*),       zet(mconu,*), array(mpru,*), buffer(*),
     & values(*)
c     # ibm+convex f77 bug; actual:
c     # integer labels(nipv,n1max)
c     # real*8 buffer(l1rec), values(n1max)
c
c     # local:
      integer ixyz, ist, iblk, ifu, icu, isfis, jsfjs, ncc, igwu,
     & ig, ic1, jc1, ic2, jc2, icc, igw, jfu, jcu, if, jf, iop,
     & ipr, itst, iprm, inx, jnx, nwt, ibl, ircr, jrcr, iopu, npairi,
     & lblop, ipair, ipr1, ipr2, num, iju, iiju
      logical  es
c
c     # bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
      integer  atebyt, forbyt
      external atebyt, forbyt
c
      integer    msame,   nmsame,   nomore
      parameter( msame=0, nmsame=1, nomore=2 )
c
      real*8     zero
      parameter( zero=0d0 )
c
      integer last
      integer itypan(7), itypbn(7), iapt(3)
      real*8 fcore
c
      character*8 chrtyp
c
c     #             s    t    v  veff so:x so:y so:z
      data itypan / 0,   0,   0,   0,   2,   2,   2 /
      data itypbn / 0,   1,   2,   3,   0,   1,   2 /
c
c     # operator symmetry
c
      call izero_argos( nst,   nopir, 1 )
      call izero_argos( 3*nst, ilxyz, 1 )
c
      if ( nnam .eq. 5 ) then
c       # x, y, and z components of spin-orbit integrals are all
c       # computed together.
        do 102 ixyz = 1, 3
          ist                   = lxyzir(ixyz)
          nopir(ist)            = nopir(ist) + 1
          ilxyz(nopir(ist),ist) = ixyz
102     continue
        nop     = 3
        iapt(1) = 0
        iapt(2) = 0
        iapt(3) = 0
      else
        nopir(1)   = 1
        ilxyz(1,1) = 1
        nop        = 1
        iapt(1)    = 0
      endif
c
      iblk = 0
      ijsf = 0
      isf  = 0
      do 712 is = 1,ns
      ifu = nf(is)
      if(ifu.eq.0) go to 712
      icu=nc(is)
      isfis=isf
      jsf = 0
      do 708 js = 1,is
      jsfjs=jsf
      es = is.eq.js
      ic=icu
      ncc = 0
      if(.not.es) go to 200
c
c  center combinations for is=js
c
      do 140 jc = 1, ic
      igwu=0
      do 132 ig = 1, ng
        ic1 = ica(ic,is,ig)
        jc1 = ica(jc,is,ig)
        if(ic1.ne.ic.and.jc1.ne.ic) go to 132
        ic2=max(ic1,jc1)
        jc2=min(ic1,jc1)
c
        do 112 icc=ncc,1,-1
          if(jc2-icb(2,1,icc)) 112,136,116
  112   continue
c
  116   if(jc2.ne.jc) go to 132
c
        do 124 igw=igwu,1,-1
          if(ic1.eq.icb(1,igw,ncc+1).and.
     &       jc1.eq.icb(2,igw,ncc+1)) go to 132
  124   continue
c
        igwu=igwu+1
        icb(1,igwu,ncc+1)=ic1
        icb(2,igwu,ncc+1)=jc1
  132 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc)=igwu
      go to 140
  136 nfct(icc) = nfct(icc) + ic
  140 continue
      do 144 icc = 1, ncc-1
        nfct(icc) = nfct(icc)/2
  144 continue
      go to 300
c
c  center combinations for is>js
c
  200 jfu = nf(js)
      if(jfu.eq.0) go to 708
      jcu = nc(js)
      do 240 jc = 1, jcu
      do 232 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 232
      jc1 = ica (jc,js,ig)
c
      do 212 icc=ncc,1,-1
        if(jc1-icb(2,1,icc)) 212,236,232
  212 continue
c
  232 continue
      ncc = ncc + 1
      icb(1,1,ncc)=ic
      icb(2,1,ncc)=jc
      nfct(ncc) = ic
      go to 240
  236 nfct(icc) = nfct(icc) + ic
  240 continue
c
c  evaluate integrals
c
  300 if(ncc.gt.mccu) then
        call bummer('change mccup (one place) to ',ncc,faterr)
      endif
      isf=isfis
      do 704 if = 1, ifu
      isf = isf + 1
      icons=mcons(isf)
      ircru=nrcr(icons)
      lit=lmnp1(icons)
      igu=ncon(icons)
      itl=ntl(isf)
      itu=ntu(isf)
      jsf=jsfjs
      if(es) jfu=if
      do 700 jf = 1, jfu
      jsf = jsf + 1
      esf=isf.eq.jsf
      jcons=mcons(jsf)
      jrcru=nrcr(jcons)
      ljt=lmnp1(jcons)
      ijsf = ijsf + 1
      ibl1=0
      ibl2=0
      do 308 ist=1,nst
        iop=nopir(ist)
        if(iop.ne.0) then
          if(esf) ibl1=ibl1+iop*nprir(1,ist,ijsf)
          ibl2=ibl2+iop*nprir(2,ist,ijsf)
        endif
  308 continue
      iju=nt(isf)*nt(jsf)
      ntij2=iju*npair(2,ijsf)
      icxi2=icxst(2,ijsf)
      if(esf) then
        nblt1=npair(1,ijsf)
        ntij1=iju*nblt1
        nblu=ircru*ibl1+ipq(ircru)*ibl2
        if(nnam.eq.5) then
          icxi1=icxast(ijsf)
        else
          icxi1=icxst(1,ijsf)
        endif
      else
        nblu=ircru*jrcru*ibl2
      endif
      if(nblu.eq.0) go to 696
c     # allocate space.
c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(1:mcxu),
c     #          4:array(1:mpru,1:narray), 5:iarray(1:iamax,1:narray),
c     #          6:buffer(1:l1rec), 7:values(1:n1max),
c     #          8:labels(1:2,1:n1max), 9:hpt(1:ihwsp), 10:hwt(1:ihwsp),
c     #          11:dfac(1:ndfac), 12:binom(1:xyzdim), 13:lmf(1,l1max2),
c     #          14:lml(1:l1max2), 15:lmx(1:lmnpwr), 16:lmy(1:lmnpwr),
c     #          17:lmz(1:lmnpwr), 18:zlm(1:lmnpwr),
c     #          19:flmtx(1:3,1:lmnu2), 20:mc(1:3,1:l2m1),
c     #          21:mr(1:3,1:l2m1), 22:h2(1:nblu), 23:ijx(1:iiju),
c     #          24:ijy(1,iiju), 25:ijz(1:iiju), 26:g2(1:ic1)
      ipt(23) = ipt(22) + atebyt( nblu )
      if( nnam.le.3 ) then
        iiju = iju
      else
        iiju = 0
      endif
      ipt(24) = ipt(23) + forbyt( iiju )
      ipt(25) = ipt(24) + forbyt( iiju )
      ipt(26) = ipt(25) + forbyt( iiju )
      ipr=iprst(ijsf)
c
      if ( (ipt(26)-1) .gt. mblu ) then
        call bummer('oneint: mblu too small ',(ipt(26)-1-mblu),faterr)
      endif
c
c     # push the current alloction onto the high-water mark stack.
      call h2opsh( ipt(26)-ipt(22) )
c
      jgu=ncon(jcons)
      jtl=ntl(jsf)
      jtu=ntu(jsf)
      call wzero(nblu,h2,1)
      if ( (.not. es) .or. esf ) then
c
        do 424 icc = 1, ncc
          ic = icb(1,1,icc)
          jc = icb(2,1,icc)
          fnfct = nfct(icc)
            call stvcz( a,          ccr,        chg,        cx,
     &      eta,        a(ipt(26)), h2,         a(ipt(23)), a(ipt(24)),
     &      a(ipt(25)), ilxyz,      ipt,        lcr,        lls,
     &      a(ipt(1)),  lproju,     mcrs,       nc,         ncr,
     &      nkcrl,      nkcru,      nklsl,      nklsu,      nopir,
     &      nprir,      nt,         x,          y,          z,
     &      zcr,        zet )
424     continue
      else
        do 524 icc = 1, ncc
          igwu = ngw(icc)
          itst = 0
          do 520 iprm = 1, 2
            if ( iprm .eq. 1 ) then
              inx=1
              jnx=2
            else
              inx=2
              jnx=1
            endif
            jc = icb(jnx,1,icc)
            if ( jc .eq. itst ) go to 520
            ic = icb(inx,1,icc)
            nwt = 1
            do 512 igw = 2, igwu
              if( icb(jnx,igw,icc) .lt. jc ) go to 520
              nwt = nwt + 1
512         continue
            fnfct = (nwt * nfct(icc))
            call stvcz( a,          ccr,        chg,        cx,
     &      eta,        a(ipt(26)), h2,         a(ipt(23)), a(ipt(24)),
     &      a(ipt(25)), ilxyz,      ipt,        lcr,        lls,
     &      a(ipt(1)),  lproju,     mcrs,       nc,         ncr,
     &      nkcrl,      nkcru,      nklsl,      nklsu,      nopir,
     &      nprir,      nt,         x,          y,          z,
     &      zcr,    zet )
            itst = jc
520       continue
524     continue
      endif
c
c     # pop back to the initial level.
      call h2opop
c
c     # put integrals and labels into the intermediate storage arrays.
c
      ibl = 0
      do 642 ircr = 1, ircru
        if(esf) jrcru = ircr
        do 640 jrcr = 1, jrcru
          do 636 ist = 1, nst
            iopu = nopir(ist)
            if ( esf .and. jrcr .eq. ircr ) then
              npairi = nprir(1,ist,ijsf)
            else
              npairi = nprir(2,ist,ijsf)
            endif
            if ( (iopu .ne. 0) .and. (npairi .ne. 0) ) then
              do 628 iop = 1, iopu
                lblop = ilxyz(iop,ist)
                do 624 ipair = 1, npairi
                  ibl = ibl + 1
c
c                  if ( abs(h2(ibl)) .gt. cutoff) then
                    iapt(lblop)                 = iapt(lblop) + 1
                    array(  iapt(lblop), lblop) = h2(ibl)
                    iarray( iapt(lblop), lblop) = (ipr+ipair)
c                  endif
624             continue
628           continue
            endif
            ipr = ipr + npairi
636       continue
640     continue
642   continue
696   continue
700   continue
704   continue
708   continue
712   continue
c
c     # write out the components of the operator arrays.
c
      do 860 lblop = 1, nop
c
c       # check that array(*) was correctly allocated.
        if ( iapt(lblop) .gt. mpru ) then
           call bummer('oneint: array(*) overflow, iapt=',
     &      iapt(lblop), faterr )
        endif
c
c       # initialize the integral record counters for this array.
c
        call w1stup
c
c       # SIFS frozen core for this array.
        fcore = zero
c
c       # SIFS integral type.
        itypea = itypan(nnam + lblop - 1)
        itypeb = itypbn(nnam + lblop - 1)
c
c       # loop over output records.
c       # do while ( ipr1 .le. iapt(lblop) )...
        ipr1 = 1
820     if ( ipr1 .le. iapt(lblop) ) then
          num  = min( (n1max - ibuf), (iapt(lblop) - ipr1 + 1) )
          ipr2 = ipr1 + num - 1
c
c         # fill up the record.
c
          do 840 ipr = ipr1, ipr2
            ibuf           = ibuf + 1
            labels(1,ibuf) =      il(iarray(ipr,lblop)) / ilfact
            labels(2,ibuf) = mod( il(iarray(ipr,lblop)),  ilfact )
            values(ibuf)   = array(ipr,lblop)
840       continue
          ipr1 = ipr2 + 1
c
c         # write out the record.
c
          if ( ipr2 .ne. iapt(lblop ) ) then
c           # more elements of this component still to go.
            last = msame
          elseif ( (nnam .eq. inam(nu)) .and. (lblop .eq. nop) ) then
c           # last component of the last 1-e operator.
            last = nomore
          else
c           # last record of this array, but more 1-e arrays to go.
            last = nmsame
          endif
          call wtlab1( last, fcore, buffer, values, labels )
c
c         # note: the "do while" construct is used since it is not
c         #       guaranteed that ibuf=0 on return from wtlab1().
c
          go to 820
        endif
c
        call siftyp( itypea, itypeb, chrtyp )
        write (nlist,6010) numout, chrtyp, nrec
c
860   continue
c
      if ( nnam .eq. 5 ) then
c       # all done with the antisymmetric coefficients.
c       # release the workspace, and adjust mcxu and mblu.
c       # this should probably be done in the higher-level
c       # calling program, but for now it's done down here.
        mblu = mblu + mcxu - mcxu2
        mcxu = mcxu2
      endif
c
      return
6010  format(' oneint:',i6,1x,a,' integrals were written in',
     & i3,' records.')
      end
      subroutine cortab(binom,dfac,eps,flmtx,hpt,hwt,lmf,lml,lmx,lmy,
     &  lmz,lmax,lmn1u,lproju,mc,mr,ndfac,zlm)
c     # tables for core potential and spin-orbit integrals.
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a2=2.0d0, a3=3.0d0)
      dimension binom(*),dfac(*),flmtx(3,*),hpt(*),hwt(*),lmf(*),
     &  lml(*),lmx(*),lmy(*),lmz(*),mc(3,*),mr(3,*),zlm(*)
c
c     # compute gauss-hermite points and weights for c, z integrals.
      igh =1
      nn = 5
      do 120 i = 1, 3
        call hermit(nn,hpt(igh),hwt(igh),eps)
        igh = igh + nn
        nn = 2*nn
  120 continue
c     # compute double factorials.
      dfac(1)=a1
      dfac(2)=a1
      fi=a0
      do 130 i=1,ndfac-2
        fi=fi+a1
        dfac(i+2)=fi * dfac(i)
  130 continue
c     # compute binomial coefficients.
      inew=1
      binom(1)=a1
      do 150 j=1,lmn1u-1
        inew=inew+1
        binom(inew)=a1
        do 140 i=1,j-1
          inew=inew+1
          binom(inew)=((j-i+1)*binom(inew-1))/(i)
  140   continue
        inew=inew+1
        binom(inew)=a1
  150 continue
c     # compute tables by recursion for real spherical harmonics.  they
c     # are indexed by l, m and sigma.  the sequence number of the
c     # harmonic with quantum numbers l, m and sigma is given by
c     #            l**2+2*m+1-sigma
c     # lmf(index) and lml(index) hold the positions of the first and
c     # last terms of the harmonic in the lmx, lmy, lmz, and zlm arrays.
c     # the harmonics with angular momentum l are generated from those
c     # with angular momenta l-1 and l-2.
c     # for m = 0,1,2,...,l-1, the recursion relation
c     z*Z(l-1,m,s) = sqrt(((l-m)*(l+m))/((2*l-1)*(2*l+1)))*Z(l,m,s)+
c                  sqrt(((l+m-1)*(l-m-1))/((2*l-3)*(2*l-1)))*Z(l-2,m,s)
c     # is used.
c     # for m = l, the recursion relation
c     x*Z(l-1,l-1,s)+(-1)**(1-s)*y*Z(l-1,l-1,1-s) =
c                  sqrt((2*l))/((2*l+1)))*Z(l,l,s)
c     # is used.
c     # l=0
      lmf(1) = 1
      lml(1) = 1
      lmx(1) = 0
      lmy(1) = 0
      lmz(1) = 0
      zlm(1) = a1
c     # l=1
      lmf(2) = 2
      lml(2) = 2
      lmx(2) = 0
      lmy(2) = 0
      lmz(2) = 1
      zlm(2) = sqrt(a3)
      lmf(3) = 3
      lml(3) = 3
      lmx(3) = 0
      lmy(3) = 1
      lmz(3) = 0
      zlm(3) = zlm(2)
      lmf(4) = 4
      lml(4) = 4
      lmx(4) = 1
      lmy(4) = 0
      lmz(4) = 0
      zlm(4) = zlm(2)
      nterm=4
      do 270 lang=2,lmax
        do 240 mang=0,lang-1
          anum = ((2*lang-1)*(2*lang+1))
          aden = ((lang-mang)*(lang+mang))
          coef1 = sqrt(anum/aden)
          anum = ((lang+mang-1)*(lang-mang-1)*(2*lang+1))
          aden = (2*lang-3)*aden
          coef2 = sqrt(anum/aden)
          nsigma=min(1,mang)
          do 230 isigma=nsigma,0,-1
            indexh=lang**2+2*mang+1-isigma
            lone=lang-1
            ltwo=lang-2
            ione=lone**2+2*mang+1-isigma
            itwo=ltwo**2+2*mang+1-isigma
            lmf(indexh)=lml(indexh-1)+1
            lml(indexh)=lml(indexh-1)
            nxy=(mang-isigma+2)/2
            iu=lmf(ione)+nxy-1
            do 200 i=lmf(ione),iu
              lml(indexh)=lml(indexh)+1
              j=lml(indexh)
              lmx(j)=lmx(i)
              lmy(j)=lmy(i)
              lmz(j)=lmz(i)+1
              zlm(j)=zlm(i)*coef1
              nterm=nterm+1
  200       end do
            if(ltwo.ge.mang) then
              il=iu+1
              do 210 i=il,lml(ione)
                lml(indexh)=lml(indexh)+1
                j=lml(indexh)
                k=lmf(itwo)+i-il
                lmx(j)=lmx(k)
                lmy(j)=lmy(k)
                lmz(j)=lmz(k)
                zlm(j)=zlm(i)*coef1-zlm(k)*coef2
                nterm=nterm+1
  210         enddo
              il=lml(itwo)-nxy+1
              if(mod(lang-mang,2).eq.0) then
                do 220 i=il,lml(itwo)
                  lml(indexh)=lml(indexh)+1
                  j=lml(indexh)
                  lmx(j)=lmx(i)
                  lmy(j)=lmy(i)
                  lmz(j)=lmz(i)
                  zlm(j)=-zlm(i)*coef2
                  nterm=nterm+1
  220           end do
              endif
            endif
  230     end do
  240   end do
        anum = (2*lang+1)
        aden = (2*lang)
        coef = sqrt(anum/aden)
        mang=lang
        isigma=1
        indexh=lang**2+2*mang+1-isigma
        lmf(indexh)=lml(indexh-1)+1
        lml(indexh)=lml(indexh-1)
c       # isig:  index of the harmonic (l-1),(m-1),sigma
c       # isigm: index of the harmonic (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        k=lmf(isigm)
        do 250 i=lmf(isig),lml(isig)
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(i)+1
          lmy(j)=lmy(i)
          lmz(j)=lmz(i)
          zlm(j)=(zlm(i)+zlm(k))*coef
          k=k+1
          nterm=nterm+1
  250   enddo
        if(mod(mang,2).eq.1) then
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(k)
          lmy(j)=lmy(k)+1
          lmz(j)=lmz(k)
          zlm(j)=zlm(k)*coef
          nterm=nterm+1
        endif
        isigma=0
        indexh=lang**2+2*mang+1-isigma
c       # isig:  index of the harmonic (l-1),(m-1),sigma
c       # isigm: index of the harmonuc (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        lmf(indexh)=lml(indexh-1)+1
        lml(indexh)=lmf(indexh)
        j=lml(indexh)
        i=lmf(isig)
        lmx(j)=lmx(i)+1
        lmy(j)=lmy(i)
        lmz(j)=lmz(i)
        zlm(j)=zlm(i)*coef
        nterm=nterm+1
        k=lmf(isigm)
        do 260 i=lmf(isig)+1,lml(isig)
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(i)+1
          lmy(j)=lmy(i)
          lmz(j)=lmz(i)
          zlm(j)=(zlm(i)-zlm(k))*coef
          k=k+1
          nterm=nterm+1
  260   enddo
        if(mod(mang,2).eq.0) then
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          k=lml(isigm)
          lmx(j)=lmx(k)
          lmy(j)=lmy(k)+1
          lmz(j)=lmz(k)
          zlm(j)=-zlm(k)*coef
          nterm=nterm+1
        endif
  270 end do
      ixy = 0
      iz = 0
      do 300 lang=1,lproju
        do 290 mang=0,lang-1
          nsigma=min(1,mang)
          ndelta=max(0,1-mang)
          anum = ((lang-mang)*(lang+mang+1))
          aden = (2*(2-ndelta))
          coef=sqrt(anum/aden)
          do 280 isigma=nsigma,0,-1
            isign=2*isigma-1
            ixy = ixy+1
            flmtx(1,ixy) = (isign)*coef
            flmtx(2,ixy) = coef
            if(mang.ne.0) then
              iz=iz+1
              flmtx(3,iz) = -(mang*isigma)
            endif
  280     enddo
  290   enddo
        iz=iz+1
        flmtx(3,iz) = -(lang)
  300 enddo
c     # column and row indices for angular momentum matrix elements.
      iadd = 1
      do 310 i=1,2*lproju-1
        mc(1,i) = i
        mc(2,i) = i
        mc(3,i) = i+1
        mr(1,i) = i+iadd
        mr(2,i) = i+2
        mr(3,i) = i+2
        iadd = 4-iadd
  310 continue
      return
      end
*deck stvcz
      subroutine stvcz( a,          ccr,        chg,        cx,
     &       eta,       g2,         h2,         ijx,        ijy,
     &       ijz,       ilxyz,      ipt,        lcr,        lls,
     &       lmnv,      lproju,     mcrs,       nc,         ncr,
     &       nkcrl,     nkcru,      nklsl,      nklsu,      nopir,
     &       nprir,     nt,         x,          y,          z,
     &       zcr,       zet )
c
      implicit real*8 (a-h,o-z)
      logical chsign,ec,esf,esfb,esfc,esfbc,igueq1,jgueq1
      parameter (a1s2=0.5d0)
c
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
c
      common /parmr/ cutoff, tol
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common /dim21/ ipq(256)
      dimension a(*),ccr(*),chg(*),cx(*),eta(mrcru,mconu,*),g2(*),h2(*),
     1  ijx(*),ijy(*),ijz(*),ilxyz(3,*),ipt(*),lcr(*),
     2  lls(*),lmnv(3,*),mcrs(*),nc(*),ncr(*),nkcrl(6,*),nkcru(6,*),
     3  nklsl(4,*),nklsu(4,*),nopir(*),nprir(2,mstu,*),nt(*),x(mcu,*),
     4  y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
c
c  bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
      integer  atebyt, forbyt
      external atebyt, forbyt
c
      lijt = lit*ljt
      ec   = ic.eq.jc
      esfc = esf.and.ec
      icx2 = icxi2+((ic-1)*nc(js)+jc-1)*ntij2
      if(.not.esfc.and.ircru.ne.1.and.esf) then
        ibld=ibl2-ibl1
        icx3=icx2-(ic-jc)*(nc(is)-1)*ntij2
      endif
      if(esf) icx1=icxi1+ipq((ic-1)*nt(isf)+1)*nblt1+(jc-1)*ntij1
c
c  ishell
c
      xi=x(ic,is)
      yi=y(ic,is)
      zi=z(ic,is)
c
c  jshell
c
      xj=x(jc,js)
      yj=y(jc,js)
      zj=z(jc,js)
      xij=a1s2*(xi+xj)
      dxij=xi-xj
      xijm=a1s2*dxij
      yij=a1s2*(yi+yj)
      dyij=yi-yj
      yijm=a1s2*dyij
      zij=a1s2*(zi+zj)
      dzij=zi-zj
      zijm=a1s2*dzij
      rr=dxij*dxij+dyij*dyij+dzij*dzij
c
c  prepare items for pairs of (i,j) functions
c
      ij=0
      do 20 it=itl,itu
        nx=ljt*lmnv(1,it)+1
        ny=ljt*lmnv(2,it)+1
        nz=ljt*lmnv(3,it)+1
        do 16 jt=jtl,jtu
          ij=ij+1
          ijx(ij)=nx+lmnv(1,jt)
          ijy(ij)=ny+lmnv(2,jt)
          ijz(ij)=nz+lmnv(3,jt)
   16   continue
   20 continue
      ij=nop*ij
      nc1=jrcru*ij
      if(esfc) then
        nc2=ipq(ircru+1)*ij
      else
        nc2=ircru*nc1
      endif
      igueq1=igu.eq.1
      jgueq1=jgu.eq.1
      if(igueq1) then
        ic1=0
      else
        ic1=nc2
      endif
      if(jgueq1) then
        ic0=0
      else
        ic0=nc1
      endif
c     # allocate space.
c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(1:mcxu),
c     #          4:array(1:mpru,1:narray), 5:iarray(1:iamax,1:narray),
c     #          6:buffer(1:l1rec), 7:values(1:n1max),
c     #          8:labels(1:2,1:n1max), 9:hpt(1:ihwsp), 10:hwt(1:ihwsp),
c     #          11:dfac(1:ndfac), 12:binom(1:xyzdim), 13:lmf(1,l1max2),
c     #          14:lml(1:l1max2), 15:lmx(1:lmnpwr), 16:lmy(1:lmnpwr),
c     #          17:lmz(1:lmnpwr), 18:zlm(1:lmnpwr),
c     #          19:flmtx(1:3,1:lmnu2), 20:mc(1:3,1:l2m1),
c     #          21:mr(1:3,1:l2m1), 22:h2(1:nblu), 23:ijx(1:iiju),
c     #          24:ijy(1,iiju), 25:ijz(1:iiju), 26:g2(1:ic1),
c     #          27:g1(1:ic0), 28:gout(1:ij),
c     #        if s, t, or v integrals, then
c     #          29:xin(1:nxyzin), 30:yin(1:nxyzin), 31:zin(1:nxyzin),
c     #          32:uf(1:nroots), 33:wf(1,nroots), 34:ff(1:nroots+n1),
c     #          (cf and sf must be kept together for droot)
c     #          35:cf(1:n1,1:n1), 36:sf(1:n1,1:n1), 37:af(1:n1),
c     #          38:rt(1:n1), 39:r(1:nroots,1:nroots),
c     #          40:w(1:nroots,1:nroots)
c     #        elseif c or z integrals, then
c     #          29:crda(1:lit,1:3), 30:crdb(1:ljt,1:3),
c     #          (for pseud1)
c     #          31:ang(1:ltot1,1:ltot1), 32:q(1:ltot1,1:ltot1),
c     #          33:qsum(1:ltot1,1:ltot1), 34:xab(1:ltot1),
c     #          35:yab:(1:ltot1), 36:zab(1:ltot1),
c     #          (for pseud2 and pseud3)
c     #          37:anga(1:lit,1:mproju,1:lamau),
c     #          38:angb(1:ljt,1:mproju,1:lambu),
c     #          39:qsum(1:ltot1,1:lambu,1:lamau), 40:apwr(1:ljt),
c     #          41:aterm1(1:ljt), 42:aterm2(1:ljt), 43:bpref(1:ljt),
c     #          44:bpwr(1:lit), 45:bterm1(1:lit), 46:ssi(ltot1+lproju),
c     #          47:abess(1:lamau), 48:bbess(1:lambu),
c     #          49:ptpow(1:ltot1), 50:q2(1:lambu,1:lamau)
c     #        endif
      ipt(27) = ipt(26) + atebyt( ic1 )
      ipt(28) = ipt(27) + atebyt( ic0 )
      ipt(29) = ipt(28) + atebyt( ij )
      go to (22,24,26,30,30), nnam
   22 nroots = 0
      nxyzin = lijt
      go to 28
   24 nroots = 0
      nxyzin = 2*lijt
      go to 28
   26 nroots = (lit+ljt)/2
      nxyzin = nroots*lijt
   28 n1 = nroots + 1
      ipt(30) = ipt(29) + atebyt( nxyzin )
      ipt(31) = ipt(30) + atebyt( nxyzin )
      ipt(32) = ipt(31) + atebyt( nxyzin )
      ipt(33) = ipt(32) + atebyt( nroots )
      ipt(34) = ipt(33) + atebyt( nroots )
      ipt(35) = ipt(34) + 2*atebyt( nroots + n1 )
      ipt(36) = ipt(35) + 2*atebyt( n1**2 )
      ipt(37) = ipt(36) + 2*atebyt( n1**2 )
      ipt(38) = ipt(37) + 2*atebyt( n1 )
      ipt(39) = ipt(38) + 2*atebyt( n1 )
      ipt(40) = ipt(39) + 2*atebyt( nroots**2 )
      ipt(51) = ipt(40) + 2*atebyt( nroots**2 )
      go to 32
   30 ipt(30) = ipt(29) + atebyt( 3*lit )
      ipt(31) = ipt(30) + atebyt( 3*ljt )
      ltot1 = lit+ljt-1
      ltot12 = ltot1**2
      ipt(32) = ipt(31) + atebyt( ltot12 )
      ipt(33) = ipt(32) + atebyt( ltot12 )
      ipt(34) = ipt(33) + atebyt( ltot12 )
      ipt(35) = ipt(34) + atebyt( ltot1 )
      ipt(36) = ipt(35) + atebyt( ltot1 )
      ipt(37) = ipt(36) + atebyt( ltot1 )
      mproju = 2*lproju+1
      lamau = lit+lproju
      ipt(38) = ipt(37) + atebyt( lit*mproju*lamau )
      lambu = ljt+lproju
      ipt(39) = ipt(38) + atebyt( ljt*mproju*lambu )
      ipt(40) = ipt(39) + atebyt( ltot1*lambu*lamau )
      ipt(41) = ipt(40) + atebyt( ljt )
      ipt(42) = ipt(41) + atebyt( ljt )
      ipt(43) = ipt(42) + atebyt( ljt )
      ipt(44) = ipt(43) + atebyt( ljt )
      ipt(45) = ipt(44) + atebyt( lit )
      ipt(46) = ipt(45) + atebyt( lit )
      ipt(47) = ipt(46) + atebyt( ltot1+lproju )
      ipt(48) = ipt(47) + atebyt( lamau )
      ipt(49) = ipt(48) + atebyt( lambu )
      ipt(50) = ipt(49) + atebyt( ltot1 )
      ipt(51) = ipt(50) + atebyt( lambu*lamau )
   32 continue
c
      if ( (ipt(51)-1) .gt. mblu ) then
        call bummer('oneint: mblu too small ',(ipt(51)-1-mblu),faterr)
      endif
c     if ( mgout .gt. mblu ) then
c       call bummer('stvcz: (mgout-mblu)=',(mgout-mblu),faterr)
c     endif
c
c     # set the high-water mark.
c     call h2oset( mgout )
      call h2oset( ipt(51)-ipt(26) )
c
      call wzero( nc2, g2, 1 )
c
      go to (36,38,40,42,54), nnam
   36       call sints( eta(1,1,icons), eta(1,1,jcons), g2,
     &      a(ipt(27)), a(ipt(28)),    	a(ipt(9)),      a(ipt(10)),
     &      ijx,        ijy,            ijz,            a(ipt(29)),
     &      a(ipt(30)), a(ipt(31)),     zet )
      go to 58
   38       call tints( eta(1,1,icons), eta(1,1,jcons), g2,
     &      a(ipt(27)), a(ipt(28)),     a(ipt(9)),      a(ipt(10)),
     &      ijx,        ijy,            ijz,            a(ipt(29)),
     &      a(ipt(30)), a(ipt(31)),     zet )
      go to 58
   40       call vints( a,              chg,            eta(1,1,icons),
     &  eta(1,1,jcons), g2,             a(ipt(27)),     a(ipt(28)),
     &  a(ipt(9)),      a(ipt(10)),     ijx,            ijy,
     &  ijz,            ipt,            nc,             nroots,
     &  n1,             a(ipt(32)),     a(ipt(33)),     x,
     &  a(ipt(29)),     y,              a(ipt(30)),     z,
     &  a(ipt(31)),     zet )
      go to 58
   42       call cints( a,              ccr,            a(ipt(29)),
     &                  a(ipt(30)),     eta(1,1,icons), eta(1,1,jcons),
     &                  g2,             a(ipt(27)),     a(ipt(28)),
     &                  ipt,            lcr,            lambu,
     &                  ltot1,          mcrs,           mproju,
     &                  nc,             ncr,            nkcrl,
     &                  nkcru,          x,              y,
     &                  z,              zcr,            zet )
      go to  58
   54       call lsints( a,              ccr,            a(ipt(29)),
     &                   a(ipt(30)),     eta(1,1,icons), eta(1,1,jcons),
     &                   g2,             a(ipt(27)),     a(ipt(28)),
     &                   ipt,            lls,            lambu,
     &                   ltot1,          mcrs,           mproju,
     &                   nc,             ncr,            nklsl,
     &                   nklsu,          x,              y,
     &                   z,              zcr,            zet )
   58 ndx=0
      ndxa=0
      do 92 ircr=1,ircru
      if(esfc) jrcru=ircr
      do 90 jrcr=1,jrcru
      esfb=esf.and.ircr.eq.jrcr
      esfbc=esfb.and.ec
      chsign=nnam.ge.5.and.esf.and.ircr.lt.jrcr
      if(.not.esfb) go to 60
      icx=icx1
      go to 66
   60 if((.not.esf).or.esfc) go to 64
      if(ircr.gt.jrcr) go to 62
      ndxa=(ipq(jrcr)+(ircr-1))*ibl2-(jrcr-1)*ibld
      icx=icx3
      ngti=nop
      ngtj=(nop*nt(isf))
      go to 70
   62 if(jrcr.eq.1) ndxa=(ipq(ircr)+(jrcr-1))*ibl2-(ircr-1)*ibld
   64 icx=icx2
   66 ngti=(nop*nt(jsf))
      ngtj=nop
   70 ndxi=ndx
      do 88 it=itl,itu
      ndxj=ndxi
      if(esfbc) jtu=it
      ndxi=ndxi+ngti
      do 86 jt=jtl,jtu
        index=ndxj
        ndxj=ndxj+ngtj
        ndxb=ndxa
        do 84 ist=1,nst
          iopu=nopir(ist)
          if(esfb) then
            npr=nprir(1,ist,ijsf)
          else
            npr=nprir(2,ist,ijsf)
          endif
          do 80 iop=1,iopu
            val=g2(index+ilxyz(iop,ist))
            if(chsign) val=-val
            val=fnfct*val
            do 76 ibl=1,npr
              h2(ndxb+ibl)=h2(ndxb+ibl)+cx(icx+ibl)*val
   76       continue
   78       ndxb=ndxb+npr
   80     continue
   82     icx=icx+npr
   84   continue
   86 continue
   88 continue
      ndxa=ndxb
      ndx=ndx+ij
   90 continue
   92 continue
      return
      end
*deck sints
      subroutine sints( etai, etaj, g2, g1, gout, hpt, hwt, ijx, ijy,
     &  ijz, xin, yin, zin, zet )
      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a1s2=0.5d0, a1=1.0d0)
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension etai(mrcru,*), etaj(mrcru,*), g2(*), g1(*), gout(*),
     1  hpt(*), hwt(*), ijx(*), ijy(*), ijz(*), xin(*), yin(*), zin(*),
     2  zet(mconu,*)
c
      if(igueq1) then
        fctr1=etai(1,1)
      else
        fctr1=a1
      endif
      if(jgueq1) fctr1=fctr1*etaj(1,1)
c
c  i primitive
c
      do 490 ig=1,igu
      if(.not.igueq1) call wzero(nc1,g1,1)
      ai=zet(ig,icons)
      arri=ai*rr
c
c  j primitive
c
      if(esfc) jgu=ig
      do 390 jg=1,jgu
      if(.not.jgueq1) call wzero(ij,gout,1)
      aj=zet(jg,jcons)
      aa=ai+aj
      aaa=(ai-aj)/aa
      dum=aj*arri/aa
      if(dum.gt.tol) go to 390
c
c  form overlap integrals
c
      x0=xij+aaa*xijm
      y0=yij+aaa*yijm
      z0=zij+aaa*zijm
      t=sqrt(aa)
      fctr2=fctr1*exp(-dum)/(t*aa)
      if((esfc.and..not.igueq1).and.ig.eq.jg) fctr2=a1s2*fctr2
      nij=0
      do 330 in=1,lit
        ni=in
        do 320 jn=1,ljt
          nj=jn
          nij=nij+1
          call stvint(hpt,hwt)
          xin(nij)=xint
          yin(nij)=yint
          zin(nij)=zint
  320   continue
  330 continue
      do 350 i=1,ij
        gout(i)=gout(i)+fctr2*xin(ijx(i))*yin(ijy(i))*zin(ijz(i))
  350 continue
c
c  j transformation
c
      if(jgueq1) go to 400
c
      j1=0
      do 380 jrcr=1,jrcru
        do 370 i=1,ij
          g1(j1+i)=g1(j1+i)+gout(i)*etaj(jrcr,jg)
  370   continue
        j1=j1+ij
  380 continue
  390 continue
c
c  i transformation
c
  400 if(igueq1) return
c
      if(esfc) i1=0
      ij2=0
      do 470 ircr=1,ircru
        j1=0
        if(esfc) then
          do 430 jrcr=1,ircr
            do 410 i=1,ij
              g2(ij2+i)=g2(ij2+i)+g1(j1+i)*etai(ircr,ig)+
     1                            g1(i1+i)*etai(jrcr,ig)
  410       continue
            ij2=ij2+ij
            j1=j1+ij
  430     continue
          i1=i1+ij
        else
          do 460 jrcr=1,jrcru
            do 450 i=1,ij
              g2(ij2+i)=g2(ij2+i)+g1(j1+i)*etai(ircr,ig)
  450       continue
            ij2=ij2+ij
            j1=j1+ij
  460     continue
        endif
  470 continue
  490 continue
      return
      end
*deck tints
      subroutine tints( etai, etaj, g2, g1, gout, hpt, hwt, ijx, ijy,
     &  ijz, xin, yin, zin, zet )
      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a2=2.0d0)
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension etai(mrcru,*), etaj(mrcru,*), g2(*), g1(*), gout(*),
     1  hpt(*), hwt(*), ijx(*), ijy(*), ijz(*), xin(*), yin(*), zin(*),
     2  zet(mconu,*)
c
      lijt = lit*ljt
      if(igueq1) then
        fctr1=etai(1,1)
      else
        fctr1=a1
      endif
      if(jgueq1) fctr1=fctr1*etaj(1,1)
c
c  i primitive
c
      do 490 ig=1,igu
      if(.not.igueq1) call wzero(nc1,g1,1)
      ai=zet(ig,icons)
      arri=ai*rr
c
c  j primitive
c
      if(esfc) jgu=ig
      do 390 jg=1,jgu
      if(.not.jgueq1) call wzero(ij,gout,1)
      aj=zet(jg,jcons)
      aa=ai+aj
      aaa=(ai-aj)/aa
      dum=aj*arri/aa
      if(dum.gt.tol) go to 390
c
c  form kinetic energy integrals
c
      x0=xij+aaa*xijm
      y0=yij+aaa*yijm
      z0=zij+aaa*zijm
      t=sqrt(aa)
      fctr2=fctr1*exp(-dum)/(t*aa**3)
      if((esfc.and..not.igueq1).and.ig.eq.jg) fctr2=a1s2*fctr2
      dum=ai*aj
      tx=dum*dxij
      ty=dum*dyij
      tz=dum*dzij
      dum=dum*aa
      txx=dum-a2*tx*tx
      tyy=dum-a2*ty*ty
      tzz=dum-a2*tz*tz
      nij=0
      tii=a0
      ti=a0
      do 330 in=1,lit
      ni=in
      tjj=a0
      tj=a0
      do 320 jn=1,ljt
      nj=jn
      nij=nij+1
      call stvint(hpt,hwt)
      xin(nij)=xint
      xin(nij+lijt)=txx*xint
      yin(nij)=yint
      yin(nij+lijt)=tyy*yint
      zin(nij)=zint
      zin(nij+lijt)=tzz*zint
      if(jn.eq.1) go to 210
      dum=tj+tj
      xin(nij+lijt)=xin(nij+lijt)+dum*tx*xin(nij-1)
      yin(nij+lijt)=yin(nij+lijt)+dum*ty*yin(nij-1)
      zin(nij+lijt)=zin(nij+lijt)+dum*tz*zin(nij-1)
      if(jn.eq.2) go to 210
      xin(nij+lijt)=xin(nij+lijt)+tjj*xin(nij-2)
      yin(nij+lijt)=yin(nij+lijt)+tjj*yin(nij-2)
      zin(nij+lijt)=zin(nij+lijt)+tjj*zin(nij-2)
  210 if(in.eq.1) go to 310
      dum=ti+ti
      xin(nij+lijt)=xin(nij+lijt)+dum*tx*xin(nij-ljt)
      yin(nij+lijt)=yin(nij+lijt)+dum*ty*yin(nij-ljt)
      zin(nij+lijt)=zin(nij+lijt)+dum*tz*zin(nij-ljt)
      if(jn.eq.1) go to 220
      dum=-tj*ti
      xin(nij+lijt)=xin(nij+lijt)+dum*xin(nij-ljt-1)
      yin(nij+lijt)=yin(nij+lijt)+dum*yin(nij-ljt-1)
      zin(nij+lijt)=zin(nij+lijt)+dum*zin(nij-ljt-1)
  220 if(in.eq.2) go to 310
      xin(nij+lijt)=xin(nij+lijt)+tii*xin(nij-2*ljt)
      yin(nij+lijt)=yin(nij+lijt)+tii*yin(nij-2*ljt)
      zin(nij+lijt)=zin(nij+lijt)+tii*zin(nij-2*ljt)
  310 tjj=tjj-tj*ai
      tj=tj+ai
  320 continue
      tii=tii+ti*aj
      ti=ti-aj
  330 continue
      do 350 i=1,ij
        nx=ijx(i)
        ny=ijy(i)
        nz=ijz(i)
        gout(i)=gout(i)+fctr2*(xin(nx+lijt)*yin(ny)*zin(nz)+xin(nx)*
     1    yin(ny+lijt)*zin(nz)+xin(nx)*yin(ny)*zin(nz+lijt))
  350 continue
c
c  j transformation
c
      if(jgueq1) go to 400
c
      j1=0
      do 380 jrcr=1,jrcru
        do 370 i=1,ij
          g1(j1+i)=g1(j1+i)+gout(i)*etaj(jrcr,jg)
  370   continue
        j1=j1+ij
  380 continue
  390 continue
c
c  i transformation
c
  400 if(igueq1) return
c
      if(esfc) i1=0
      ij2=0
      do 470 ircr=1,ircru
        j1=0
        if(esfc) then
          do 430 jrcr=1,ircr
            do 410 i=1,ij
              g2(ij2+i)=g2(ij2+i)+g1(j1+i)*etai(ircr,ig)+
     1                            g1(i1+i)*etai(jrcr,ig)
  410       continue
            ij2=ij2+ij
            j1=j1+ij
  430     continue
          i1=i1+ij
        else
          do 460 jrcr=1,jrcru
            do 450 i=1,ij
              g2(ij2+i)=g2(ij2+i)+g1(j1+i)*etai(ircr,ig)
  450       continue
            ij2=ij2+ij
            j1=j1+ij
  460     continue
        endif
  470 continue
  490 continue
      return
      end
*deck vints
      subroutine vints( a,              chg,            etai,
     &  etaj,           g2,             g1,             gout,
     &  hpt,            hwt,            ijx,            ijy,
     &  ijz,            ipt,            nc,             nroots,
     &  n1,             uf,             wf,             x,
     &  xin,            y,              yin,            z,
     &  zin,            zet )
      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a1s2=0.5d0, asrtpi=1.1283791670955126d0)
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension a(*), chg(*), etai(mrcru,*), etaj(mrcru,*), g2(*),
     1  g1(*), gout(*), hpt(*), hwt(*), ijx(*), ijy(*), ijz(*), ipt(*),
     2  nc(*), uf(*), wf(*), x(mcu,*), xin(*), y(mcu,*), yin(*),
     3  z(mcu,*), zin(*), zet(mconu,*)
c
      lijt = lit*ljt
      mmu=lijt*nroots
      if(igueq1) then
        fctr1=asrtpi*etai(1,1)
      else
        fctr1=asrtpi
      endif
      if(jgueq1) fctr1=fctr1*etaj(1,1)
c
c  i primitive
c
      do 490 ig=1,igu
      if(.not.igueq1) call wzero(nc1,g1,1)
      ai=zet(ig,icons)
      arri=ai*rr
c
c  j primitive
c
      if(esfc) jgu=ig
      do 390 jg=1,jgu
      if(.not.jgueq1) call wzero(ij,gout,1)
      aj=zet(jg,jcons)
      aa=ai+aj
      aaa=(ai-aj)/aa
      dum=aj*arri/aa
      if(dum.gt.tol) go to 390
c
c  form nuclear attraction integrals
c
      ax=xij+aaa*xijm
      ay=yij+aaa*yijm
      az=zij+aaa*zijm
      fctr2=fctr1*exp(-dum)/aa
      if((esfc.and..not.igueq1).and.ig.eq.jg) fctr2=a1s2*fctr2
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 356 ks=1,ns
        fctr3=-chg(ks)*fctr2
        do 354 kc=1,nc(ks)
          cx=x(kc,ks)
          cy=y(kc,ks)
          cz=z(kc,ks)
          xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
          if(nroots.le.3) then
            call rt123(nroots,xx,uf,wf)
          elseif(nroots.eq.4) then
            call root4(xx,uf,wf)
          elseif(nroots.eq.5) then
            call root5(xx,uf,wf)
          else
            call droot(nroots,n1,xx,a(ipt(34)),a(ipt(35)),a(ipt(36)),
     &        a(ipt(37)),a(ipt(38)),a(ipt(39)),a(ipt(40)),uf,wf)
          endif
          mm=0
          do 340 k=1,nroots
            uu=aa*uf(k)
            ww=wf(k)*fctr3
            tt=aa+uu
            t=sqrt(tt)
            x0=(aax+uu*cx)/tt
            y0=(aay+uu*cy)/tt
            z0=(aaz+uu*cz)/tt
            do 330 in=1,lit
              ni=in
              do 320 jn=1,ljt
                nj=jn
                mm=mm+1
                call stvint(hpt,hwt)
                xin(mm)=xint
                yin(mm)=yint
                zin(mm)=zint*ww
  320         continue
  330       continue
  340     continue
          do 352 i=1,ij
            nx=ijx(i)-lijt
            ny=ijy(i)-lijt
            nz=ijz(i)-lijt
            do 350 mm=lijt,mmu,lijt
              gout(i)=gout(i)+xin(nx+mm)*yin(ny+mm)*zin(nz+mm)
  350       continue
  352     continue
  354   continue
  356 continue
c
c  j transformation
c
      if(jgueq1) go to 400
c
      j1=0
      do 380 jrcr=1,jrcru
        do 370 i=1,ij
          g1(j1+i)=g1(j1+i)+gout(i)*etaj(jrcr,jg)
  370   continue
        j1=j1+ij
  380 continue
  390 continue
c
c  i transformation
c
  400 if(igueq1) return
c
      if(esfc) i1=0
      ij2=0
      do 470 ircr=1,ircru
        j1=0
        if(esfc) then
          do 430 jrcr=1,ircr
            do 410 i=1,ij
              g2(ij2+i)=g2(ij2+i)+g1(j1+i)*etai(ircr,ig)+
     1                            g1(i1+i)*etai(jrcr,ig)
  410       continue
            ij2=ij2+ij
            j1=j1+ij
  430     continue
          i1=i1+ij
        else
          do 460 jrcr=1,jrcru
            do 450 i=1,ij
              g2(ij2+i)=g2(ij2+i)+g1(j1+i)*etai(ircr,ig)
  450       continue
            ij2=ij2+ij
            j1=j1+ij
  460     continue
        endif
  470 continue
  490 continue
      return
      end
*deck iodfac
      function iodfac(l1,l2,m1,m2,n1,n2)
c
c  compute the product of odd number factorials which gives the
c  one-center overlaps of cartesian gaussian aos.
c
      if(mod(l1+l2,2).eq.1 .or.
     &   mod(m1+m2,2).eq.1 .or.
     &   mod(n1+n2,2).eq.1) then
        iodfac = 0
      else
        iprd = 1
        do 20 i=3,l1+l2-1,2
          iprd = i*iprd
   20   continue
        do 30 i=3,m1+m2-1,2
          iprd = i*iprd
   30   continue
        do 40 i=3,n1+n2-1,2
          iprd = i*iprd
   40   continue
        iodfac = iprd
      endif
      return
      end
