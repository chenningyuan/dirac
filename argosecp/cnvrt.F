*deck cnvrt
      subroutine cnvrt
c
c  this program transforms the integral file into one that is
c  compatible with the scf program by changing the indices or labels
c  of the one-electron integrals and by calculating the p and k/2
c  elements for the two-electron integrals.
c
c  version log:
c  11-sep-91 parallel version included in columbus distribution. -rjh
c  07-jul-91 preliminary parallel version using tcgmsg. -rjh
c  25-may-91 sifs version. minor cleanup. -rmp/rls
c  08-90 allocate arrays (rmp)
c  08-89 mdc update, zero-out improvement (rmp,mpn)
c  08-88 add sun, revise for distribution (eas,rmp)
c  1983 (cnvrt)   major revisions, r. m. pitzer
c  1981 (convert) major revisions, f. b. brown (i. shavitt group)
c  1979 (soconv)  y. yamaguchi & b. brooks (h. f. schaefer group)
c
c  cmdc info:
c  keyword   description
c  -------   -----------
c  ibm       ibm code.
c  vax       vax code.
c  sun       sun code.
c  stellar   stellar code.
c  harris    harris code.
c  cray      cray code. use compiler option for no double precision.
c  crayctss  ctss specific code.
c  parallel  tcgmsg parallel code.
c
      implicit logical(a-z)
c
      integer    ntitmx,    nbfmx,     ninfmx,    nengmx,    nmapmx
      parameter( ntitmx=20, nbfmx=255, ninfmx=20, nengmx=40, nmapmx=20 )
c
      integer     nlist, ifile, ifile2, nfile, nnam, idebug
      common /c1/ nlist, ifile, ifile2, nfile, nnam, idebug
c
      real*8 repnuc
      integer  nst,    ns,     isfr,   nu,     iau,    ntitle,
     & nso,    nsopr,          nblpr,          nnbft
      common /c2/
     & repnuc, nst,    ns,     isfr,   nu,     iau,    ntitle,
     & nso(8), nsopr(8),       nblpr(8),       nnbft
c
      character*80    title
      common /ctitle/ title(ntitmx)
c
      integer       ipq,        mapin,        mapout,        symb
      common /cipq/ ipq(nbfmx), mapin(nbfmx), mapout(nbfmx), symb(nbfmx)
c
c     # core(1:lencor) is used for miscellaneous  workspace.
      integer    lencor
      parameter( lencor=1 000 000 )
c
      real*8 core(lencor)
c
      character*8 bfnlab(nbfmx)
      real*8 energy(nengmx)
      integer ietype(20), imtype(nmapmx), info(20)
      character*4  slabel(8)
c
      integer i, buffer, values, labels, l1rec,
     & n1max, l2rec, n2max, mnl, ms, map, ninput, isto,
     & ninfo, ierr, nmap, nenrgy, isft, itotal, stvc, ibitv, itemp,
     & spo, lblij, plabs, psup, ksup, ilbl, klbl
      character*60 fname
c
c     # bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
c     # nbufo1 = number of 1-e values.
c     # nbufo2 = number of 2-e values.
      integer    nbufo1,      nbufo2
*mdc*if cray
*      parameter( nbufo1=4095, nbufo2=2730 )
*mdc*else
      parameter( nbufo1=2304, nbufo2=1152 )
*mdc*endif
c
      integer  atebyt, forbyt
      external atebyt, forbyt
c
      real*8   sifsce
      external sifsce
c
      integer  nodeid
      external nodeid
c
c     # parallel message passing initialization, or stub as appropriate.
      call pbginf
c
*mdc*if crayctss
*      call link('unit99=tty//')
*mdc*endif
c
c     # open the output listing file.  write to file 'cnvrtls'.
c
      nlist=6 
*mdc*if harris ibm
*c     # use preconnected unit 6
*mdc*elseif parallel
*      fname = 'cnvrtls'
*      call trnfln( 1, fname )
*      call pfname( 1, fname )
*      open ( nlist, file=fname, status='unknown' )
*mdc*else
c      fname = 'cnvrtls'
      call trnfln( 1, fname )
c      open ( nlist, file=fname, status='unknown' )
*mdc*endif
c
      call ibummr( nlist )
c
c     # open the input file.  read from file 'cnvrtin'.
c
      ninput=5
*mdc*if harris
*c     # use preconnected unit 5
*mdc*elseif parallel
*      fname = 'cnvrtin'
*      call trnfln( 1, fname )
*      if ( nodeid() .ne. 0 ) call pfname( 1, fname )
*c     # this distributes the input file to all nodes.
*      call pfcopy( 99, 0, fname )
*      open ( ninput, file=fname, status='old' )
*mdc*else
c     # translate the filename and open.
      fname = 'cnvrtin'
      call trnfln( 1, fname )
      open ( ninput, file=fname, status='old' )
*mdc*endif
c
      nmap   = 0
      ierr   = 0
      ninfo  = 0
      isft   = 0
      nenrgy = 0
c
      isto   = 0
      ifile  = 4
      nfile  = 3
      idebug = 0
      ifile2 = 8
c-      read (ninput,'(26i3)') isto, ifile, nfile, idebug
      read (ninput,*) isto, ifile, nfile, idebug, ifile2
c
cmckim start
c     close (ninput)
cmckim end
c
      if ( ifile  .eq. 0 ) ifile  = 4
      if ( ifile2 .eq. 0 ) ifile2 = 8
      if ( nfile  .eq. 0 ) nfile  = 3
c
      write (nlist,10)
   10 format(//28x,'program "cnvrt" 4.1b1'/
     1 28x,'columbus program system'/
     2 28x,'supermatrix formation program'//
     & 28x,'Programmed by Russell M. Pitzer'/
     3 28x,'version date: 11-sep-91')
c
      call who2c( 'CNVRT', nlist )
c
c     # open ifile.  read from file 'aoints'.
c     # open nfile.  write to file 'pkints'.
c
      fname = 'aoints'
      call trnfln( 1, fname )
*mdc*if parallel
*      call pfname( 1, fname )
*mdc*endif
      open (ifile, file=fname, status='old', form='unformatted')
c
      fname = 'pkints'
      call trnfln( 1, fname )
*mdc*if parallel
*      call pfname( 1, fname )
*mdc*endif
      open (nfile,file=fname, status='unknown', form='unformatted')
c
      write (nlist,40) ifile, nfile
   40 format(/' integrals are read in from file',i3,' (ifile).'/
     1  ' integrals are written out on file',i3,' (nfile).'//)
c
      call sifrh1(
     & ifile,  ntitle, nst,    isfr,
     & ninfo,  nenrgy, nmap,   ierr )
c
      if ( ierr .ne. 0 ) then
         call bummer('from sifrh1, ierr=', ierr, faterr )
      elseif ( ntitle .gt. ntitmx ) then
         call bummer('from sifrh1, ntitle=', ntitle, faterr )
      elseif ( isft .gt. nbfmx ) then
         call bummer('from sifrh1, isfr=', isfr, faterr )
      elseif ( ninfo .gt. ninfmx ) then
         call bummer('from sifrh1, ninfo=', ninfo, faterr )
      elseif ( nenrgy .gt. nengmx ) then
         call bummer('from sifrh1, nenrgy=', nenrgy, faterr )
      elseif ( nmap .gt. nmapmx ) then
         call bummer('from sifrh1, nmap=', nmap, faterr )
      endif
c
      call izero( 8, nso, 1 )
c
c     # read the second header record.
c     # core(*)==> mnl(1:isfr), ms(1:isfr), map(1:isfr,1:nmap)
c
      mnl = 1
      ms  = mnl + forbyt( isfr )
      map = ms  + forbyt( isfr )
c
      call sifrh2(
     & ifile,  ntitle, nst,    isfr,
     & ninfo,  nenrgy, nmap,   title,
     & nso,    slabel, info,   bfnlab,
     & ietype, energy, imtype, core(map),
     & ierr )
c
      if ( ierr .ne. 0 ) then
         call bummer('from sifrh2(), ierr=', ierr, faterr )
      endif
c
      write(nlist,6100) 'input SIFS file header information:'
      write(nlist,'(1x,a)') (title(i),i=1,ntitle)
c
      write(nlist,6100) 'input energy(*) values:'
      call sifpre( 6, nenrgy, energy, ietype )
c
      repnuc = sifsce( nenrgy, energy, ietype )
      write(nlist,6200) 'total core energy =', repnuc
c
      write(nlist,6040) nst, isfr
c
      write(nlist,6050) 'symmetry  =',(i,i=1,nst)
      write(nlist,6060) 'slabel(*) =',(slabel(i),i=1,nst)
      write(nlist,6050) 'nbpsy(*)  =',(nso(i),i=1,nst)
c
      write(nlist,6100) 'info(*) =',(info(i),i=1,ninfo)
c
      write(nlist,6100) 'input orbital labels, i:bfnlab(i)='
      write(nlist,6080) (i,bfnlab(i),i=1,isfr)
c
c     # initialize some addressing arrays.
c     # extract the bfn-to-center map vector, and
c     # the bfn-to-bfntyp map vector from map(*,*).
c
c     # core(*)==> mnl(1:isfr), ms(1:isfr), map(1:isfr,1:nmap)
c
      call getmap( nmap, imtype, core(map), core(mnl), core(ms) )
c
c     # extract some file and record length parameters.
      l1rec  = info(2)
      n1max  = info(3)
      l2rec  = info(4)
      n2max  = info(5)
c
c     # process the 1-e integral arrays.
c
c     # core(*)==> mnl(1:isfr), ms(1:isfr),
c     #            buffer(1:l1rec), values(1:n1max), labels(1:2,1:n1max)
c     #            stvc(1:nnbft,0:3), spo(1:nbufo1), lblij(1:nbufo1),
c     #            plabs(1:nbufo1)
      buffer = ms     + forbyt( isfr )
      values = buffer + atebyt( l1rec )
      labels = values + atebyt( n1max )
      stvc   = labels + forbyt( 2*n1max )
c      stvcso = labels + forbyt( 2*n1max )
cmckim start
c      spo    = stvc   + atebyt( 4*nnbft )
      spo    = stvc   + atebyt( 7*nnbft )
c      sposo  = stvcso + atebyt( 3*isfr*isfr )
cmckim end
      lblij  = spo    + atebyt( nbufo1 )
      plabs  = lblij  + forbyt( nbufo1 )
      itotal = plabs  + forbyt( nbufo1 ) - 1
c
      if ( itotal .gt. lencor ) then
         call bummer('before out1e(), itotal=', itotal, faterr )
      endif
c
      call out1e(
     & l1rec,        n1max,       core(buffer), core(values),
     & core(labels), core(stvc),  info, core(ms),
     & core(mnl),    slabel,      bfnlab,       nbufo1,
     & core(spo),    core(lblij), core(plabs),  isfr, nst )
c
      call pend
c
      return
6040  format(/' nsym =',i2,' nbft=',i4/)
6050  format(1x,a,8i5)
6060  format(1x,a,8(1x,a4))
6080  format(10(i4,':',a8))
6100  format(/1x,a,(1x,10i10))
6200  format(/1x,a,(1x,1p4e20.12))
c
      end
c
*deck getmap
      subroutine getmap( nmap, imtype, map, mnl, ms)
c
c  extract mnl(*) and ms(*) from map(*,*).
c  compute ipq(*), nsopr(*), and nblpr(*).
c
c  25-may-91 written by russ pitzer.
c
      implicit logical(a-z)
c
      real*8 repnuc
      integer  nst,    ns,     isfr,   nu,     iau,    ntitle,
     & nso,    nsopr,          nblpr,          nnbft
      common /c2/
     & repnuc, nst,    ns,     isfr,   nu,     iau,    ntitle,
     & nso(8), nsopr(8),       nblpr(8),       nnbft
c
      integer    nbfmx
      parameter( nbfmx=255 )
c
      integer       ipq,        mapin,        mapout,        symb
      common /cipq/ ipq(nbfmx), mapin(nbfmx), mapout(nbfmx), symb(nbfmx)
c
c     # dummy:
      integer nmap
      integer imtype(*), mnl(isfr), ms(isfr), map(isfr,*)
c
c     # local:
      integer    mnltyp,   mstyp
      parameter( mnltyp=4, mstyp=3 )
      integer i, nbft, imnl, ims, nbi, j
c
      integer nndxf
      nndxf(i) = (i * (i - 1)) / 2
c
c     # initialize the ipq(*) array, and set the
c     # default mapin(*) and mapout(*) values.
      do 100 i = 1, nbfmx
         ipq(i)    = nndxf(i)
         mapin(i)  = i
         mapout(i) = i
         symb(i)   = 0
100   continue
c
c     # initialize the symmetry offset arrays.
      nbft  = 0
      nnbft = 0
      do 300 i = 1, 8
         nbi      = nso(i)
         nsopr(i) = nbft
         nblpr(i) = nnbft
         do 200 j = (nbft+1), (nbft+nbi)
            symb(j) = i
200      continue
         nbft  = nbft + nbi
         nnbft = nnbft + nndxf( nbi + 1 )
300   continue
c
      imnl = 0
      ims  = 0
      do 400 i = 1, nmap
         if ( imtype(i) .eq. mnltyp ) then
            imnl = i
         elseif ( imtype(i) .eq. mstyp ) then
            ims = i
         endif
400   continue
c
c     # extract the bfntyp(*) array from map(*,*).
      if ( imnl .eq. 0 ) then
c        # default: set everything to an s-type function.
         call iset( isfr, 1, mnl, 1 )
      else
         call icopy( isfr, map(1,imnl), 1, mnl, 1 )
      endif
c
c     # extract the bfn-to-center array from map(*,*).
      if ( ims .eq. 0 ) then
c        # default: set everything to one center.
         call iset( isfr, 1, ms, 1 )
      else
         call icopy( isfr, map(1,ims), 1, ms, 1 )
      endif
c
      return
      end
c
*deck out1e
      subroutine out1e(
     & l1rec,  n1max,  buffer,  values,
     & labels, stvc,   info,  ms,
     & mnl,    slabel, bfnlab,  nbufo,
     & spo,    lblij,  plabs )
c
c  read in the stvc(*,*) array and write the output file header and
c  1-e integral arrays.
c
c  25-may-91 SIFS version. -rmp/rls
c
      implicit logical(a-z)
c
      integer     nlist, ifile, ifile2, nfile, nnam, idebug
      common /c1/ nlist, ifile, ifile2, nfile, nnam, idebug
c
      integer    ntitmx,    nbfmx
      parameter( ntitmx=20, nbfmx=255 )
c
      character*80    title
      common /ctitle/ title(ntitmx)
c
      integer       ipq,        mapin,        mapout,        symb
      common /cipq/ ipq(nbfmx), mapin(nbfmx), mapout(nbfmx), symb(nbfmx)
c
      real*8 repnuc
      integer  nst,    ns,     isfr,   nu,     iau,    ntitle,
     & nso,    nsopr,          nblpr,          nnbft
      common /c2/
     & repnuc, nst,    ns,     isfr,   nu,     iau,    ntitle,
     & nso(8), nsopr(8),       nblpr(8),       nnbft
c
c     # dummy:
      integer l1rec, n1max, nbufo
      integer info(*), ms(isfr), mnl(isfr), lblij(nbufo)
      real*8 buffer(l1rec), values(n1max), labels(2,n1max),
     & stvc(nnbft,0:3), spo(nbufo)
      character*4 slabel(nst)
      character*8 bfnlab(isfr)
c
c     # see dump1() for actual plabs(*) declaration.
      integer plabs(*)
c
c     # local:
      integer i, j, ij, ierr, ist, iso, is, itypeb, ibfld, isym, ni,
     & iblko, ibufo
      real*8 score, tcore, vcore, vecore
      character*8 chrtyp
c
      integer nd(8), lxyzir(3), inam(4)
      character*3 ityp(8), mtype(nbfmx)
c
      real*8     zero,     small
      parameter( zero=0d0, small=1.d-15 )
c
c     # bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
      integer    more,   nomore
      parameter( more=0, nomore=1 )
c
c     # dummy lxyzir(*) values, default inam(*) and nd(*) values:
      data lxyzir / 3*1 /
      data inam   / 1, 2, 3, 4 /
      data nd     / 8*1 /
c
      ierr   = 0
      score  = zero
      tcore  = zero
      vcore  = zero
      vecore = zero
c
c     # write out the ms(*) and mnl(*) map vectors:
c
      write(nlist,6100) 'bfn_to_center map(*), i:map(i)'
      write(nlist,6110) (i,ms(i),i=1,isfr)
c
      write(nlist,6100) 'bfn_to_orbital_type map(*), i:map(i)'
      write(nlist,6110) (i,mnl(i),i=1,isfr)
c
c     # read in the stvc(*,*) array.
c
      call rdstvc(
     & ifile,  info,   buffer, values,
     & labels, nst,    nso,    mapin,
     & nnbft,  stvc,   score,  tcore,
     & vcore,  vecore, symb,   ierr, 
     & isfr , nst)
c
      return
6010  format(1p/' out1e:',
     & t10,'score  =',e25.15/
     & t10,'tcore  =',e25.15/
     & t10,'vcore  =',e25.15/
     & t10,'vecore =',e25.15)
6100  format(/1x,a,(1x,10i10))
6110  format(10(i4,':',i3))
6200  format(' out1e:',i12,1x,a,' integrals written in',i3,' records.')
      end
c
*deck rdstvc
      subroutine rdstvc(
     & aoints,  info,    buffer,  values,
     & labels,  nsym,    nbpsy,   mapin,
     & nnbft,   stvc,    score,   tcore,
     & vcore,   vecore,  symb,    ierr,
     & isfr, nst )
c
c  read the overlap, kinetic energy, potential energy, and
c  effective  arrays.
c
c  on exit, the integral file is positioned after the
c  last 1-e integral record.
c
c  input:
c  aoints = input file unit number.
c  info(*) = info array for this file.
c  buffer(1:l1rec) = record buffer.
c  values(1:n1max) = value buffer.
c  labels(1:2,1:n1max) = orbital label buffer.
c  nsym = number of symmetry blocks.
c  nbpsy(*) = number of basis functions per symmetry block.
c  mapin(*) = input_ao-to-ao mapping vector.
c  nnbft = leading dimension of stvc(*,1:3).
c
c  output:
c  stvc(*) = the overlap s1(*) matrix is returned in stvc(*,1), the
c           total kinetic matrix T1_total(*) is returned in stvc(*,2),
c           the 1-e potential matrix V1_total(*) is returned in
c           stvc(*,3), and any remaining hamiltonian contributions are
c           returned in stvc(*,4).
c           all arrays are returned symmetry-blocked lower-triangle
c           packed by rows.
c           all 1-e contributions on the integral file are summed
c           into this array.  consequently, the entries on the file
c           must be only the distinct array elements.
c  score  = frozen core contribution. tr( s1 * dfc ) = nfrzct
c  tcore  = frozen core contribution. tr( t1 * dfc )
c  vcore  = frozen core contribution. tr( v1 * dfc )
c  vecore = frozen core contribution. tr( h1 * dfc )
c  symb(1:nbft) = symmetry index of each basis function
c  ierr = error return code.
c       =  0 for normal return.
c       = -1 if no arrays were found on the integral file.
c       = -n if n symmetry blocking errors were detected.
c       >  0 for iostat error.
c
c  25-may-91 rdstvc() created from modified sifstv(). -rmp/rls
c  31-oct-90 sifstv() created from sifrsh(). -rls
c  08-oct-90 (columbus day) 1-e fcore change.  sifr1n() interface
c            used. ierr added. -rls
c  04-oct-90 sifskh() call added. -rls
c  26-jul-89 written by ron shepard.
c
      implicit logical(a-z)
c
      integer ia,ib,ja
      integer  aoints, nsym,   nnbft,  ierr,  isfr, nst
      integer  info(*),        nbpsy(nsym),    labels(2,*),
     & mapin(*),       symb(*)
      real*8   score,  tcore,  vcore,  vecore
      real*8   buffer(*),      values(*), stvc(nnbft,7)
     &       ,  islso(isfr,isfr,5:7)
     &        , rtotbasis, sover, xnorm,   arppint(nnbft)
      integer nbft, ntotbasis, maxdistinct
      parameter( maxdistinct=10)
      integer nuarep, nusoc, ndistinct,
     &        nequivalent(maxdistinct), nbasis(maxdistinct)
      integer lensoc
      integer ntribasis
c     # local...
      integer    itypea,   btypmx
      parameter( itypea=0, btypmx=6 )
      integer i, j,   nntot,  isym,   nrec,   last,   lastb,  lasta
      integer symtot(36), idummy(1), btypes(0:btypmx)
      real*8  fcore(7), socx, socy, socz
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer nndxf
c
c     # set the btypes(*) array:
c     #            0:s1, 1:t1, 2:v1, 3:vec, 4:vfc, 5:vref, 6:generic_h1
      data btypes/ 1,    2,    3,    4,     4,     4,      4          /
c
      nndxf(i) = (i * (i - 1)) / 2
c
c     # nntot = the actual number of elements in the arrays.
      nntot=0
      do 20 isym=1,nsym
         nntot = nntot + nndxf( nbpsy(isym) + 1 )
20    continue
c
      if ( nntot .gt. nnbft ) then
c        # inconsistent nnbft value.
         call bummer('rdstvc: (nntot-nnbft)=',(nntot-nnbft),wrnerr)
         ierr = -2
         return
      endif
c
c     # initialize the output arrays...
c
      call wzero( nntot, stvc(1,1), 1 )
      call wzero( nntot, stvc(1,2), 1 )
      call wzero( nntot, stvc(1,3), 1 )
      call wzero( nntot, stvc(1,4), 1 )
      call wzero( nntot, stvc(1,5), 1 )
      call wzero( nntot, stvc(1,6), 1 )
      call wzero( nntot, stvc(1,7), 1 )
      fcore(1) = (0)
      fcore(2) = (0)
      fcore(3) = (0)
      fcore(4) = (0)
      fcore(5) = (0)
      fcore(6) = (0)
      fcore(7) = (0)
      call sifr1nisl(
     & aoints, info,   itypea, btypmx,
     & btypes, buffer, values, labels,
     & nsym,   nbpsy,  idummy, mapin,
     & nnbft,  stvc,   fcore,  symb,
     & symtot, lasta,  lastb,  last,
     & nrec,   ierr, nst, isfr, islso )
c
c      open(56,form='unformatted',file='argos_ovr')
c      open(57,form='unformatted',file='argos_kin')
c      open(58,form='unformatted',file='argos_v')
      open(59,form='unformatted',file='argos_vcore')
      open(60,form='unformatted',file='sox')
      open(61,form='unformatted',file='soy')
      open(62,form='unformatted',file='soz')
c      
c Write integrals to file
c      
c      write(56),(stvc(i,1),i=1,nnbft)
c      write(57),(stvc(i,2),i=1,nnbft)
c      write(58),(stvc(i,3),i=1,nnbft)
      write(59),(stvc(i,4),i=1,nnbft)
      write(60),((islso(ja,ia,5),ja=1,isfr),ia=1,isfr)
      write(61),((islso(ja,ia,6),ja=1,isfr),ia=1,isfr)
      write(62),((islso(ja,ia,7),ja=1,isfr),ia=1,isfr)
c
c     Debug
c
c      do j=1,7
c        i=0
c        if (j .eq. 1) then
c          write(6,'(/a)') 'S (including zero)'
c          do ia = 1, nnbft
c            i=i+1
c            write(6,1001) ia,stvc(i,j)
c          enddo
c        else if (j .eq. 2) then
c          write(6,'(/a)') 'T (including zero)'
c          do ia = 1, nnbft
c            i=i+1
c            write(6,1001) ia,stvc(i,j)
c          enddo
c        else if (j .eq. 3) then
c          write(6,'(/a)') 'V (including zero)'
c          do ia = 1, nnbft
c            i=i+1
c            write(6,1001) ia,stvc(i,j)
c          enddo
c        else if (j .eq. 4) then
c          write(6,'(/a)') 'Vcore (including zero)'
c          do ia = 1, nnbft
c            i=i+1
c            write(6,1001) ia,stvc(i,j)
c          enddo
c        endif
c      enddo
c
c      close(56)
c      close(57)
c      close(58)
      close(59)
      close(60)
      close(61)
      close(62)
  68  format(i5,2x,i5,10x,f20.16)
 1001 format(i5,10x,f20.16)
c
c     # save the appropriate core values.
      score  = fcore(1)
      tcore  = fcore(2)
      vcore  = fcore(3)
      vecore = fcore(4)
      socx   = fcore(5)
      socy   = fcore(6)
      socz   = fcore(7)
c
 9110 format(' Number of Distinct Atoms ', i3)
 9120 format(' # of Eq. Atoms of Ith Distinct Atom ', 10i3)
 9130 format(' # of Eq. Bases of Ith Distinct Atom ', 10i3)
 9140 format(' Total # of Bases ', i4)
 9150 format(' Error :     nnbft ', i4,
     &       ' not equal to ntribasis ', i4)
      return
      end
