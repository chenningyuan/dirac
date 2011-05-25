      integer function forbyt( i )
!  compute the number of working precision words required to hold
!  an (unpacked) integer array.
      integer i
      forbyt = i
      return
      end

      integer function atebyt( i )
!  compute the number of working precision words required to hold
!  a real (working precision) array.
      integer i
      atebyt = i
      return
      end


      subroutine izero_argos(n,ia,inca)
! zero the elements of ia
      integer n, inca
      integer ia(inca,*)
      integer i
      do i = 1, n
         ia(1,i) = 0
      enddo
      return
      end


      subroutine  icopy_argos(n,ivx,incx,ivy,incy)
!  blas-style integer vector copy.
      implicit integer(a-z)
      integer ivx(*),ivy(*)
      integer i,incx,incy,ix,iy,n
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)then
!        code for both increments equal to 1
         do i = 1,n
            ivy(i) = ivx(i)
         enddo
      else
!        code for unequal increments or equal increments
!        not equal to 1
         ix = 1
         iy = 1
         if(incx.lt.0)ix = (-n+1)*incx + 1
         if(incy.lt.0)iy = (-n+1)*incy + 1
         do i = 1, n
            ivy(iy) = ivx(ix)
            ix = ix + incx
            iy = iy + incy
         enddo
      endif
      return
      end

      subroutine bummer( text, ierr, errtyp )
c
c  process a program error.
c  input:
c  text  = character string to be printed.
c  ierr  = internal program error to be printed.
c  errtyp = 0 for warning.  traceback may be generated. execution
c             continues.
c         = 1 for nonfatal error.  traceback may be generated.
c             execution is stopped. jcl condition code is set to allow
c             subsequent program steps to continue if possible.
c         = 2 for fatal error.  traceback may be generated.
c             execution is stopped. jcl condition code is set to abort
c             subsequent program steps if possible.
c
c  entry ibummr must be called prior to bummer() to set the output
c  unit and to perform any additional initialization.
c
c  version log:
c  24-apr-92 %val() blocks added for ibm rs6000. -rls
c  11-sep-91 parerr() calls added. -rjh
c  13-mar-91 posix code added. -rls
c  01-mar-89 write(stderr,*) and call exit() for unix machines (rls).
c  05-jul-88 unicos version (rls).
c  03-nov-87 ibm version (dcc).
c  10-sep-87 written by ron shepard.
c
      implicit integer(a-z)
c
      character*(*) text
      integer iunit, errtyp
 
#include "argos_error.h"
 
      integer f77err
      integer nlist,        stderr
      save    nlist,        stderr
      data    nlist / 6 /,  stderr / 0 /
c
      if ( errtyp .eq. wrnerr ) then
c        # print a warning message and continue execution.
         write(nlist,6010) 'bummer (warning):', text, ierr
         write(stderr,6010) 'bummer (warning):', text, ierr
         return
c
      elseif ( errtyp .eq. nfterr ) then
c        # print a warning message, stop execution.
         write(nlist,6010) 'bummer (nonfatal):', text, ierr
         write(stderr,6010) 'bummer (nonfatal):', text, ierr
         call exit( 0 )
c
      elseif ( errtyp .eq. faterr ) then
c        # print an error message, stop execution, and abort job
c        # sequence.
         write(nlist,6010) 'bummer (fatal):', text, ierr
         write(stderr,6010) 'bummer (fatal):', text, ierr
         call exit( 1 )
c
      else
c        # unknown error level.  treat as a fatal error.
         write(nlist,6020) 'bummer (unknown): errtyp=',
     &    errtyp, text, ierr
         write(stderr,6020) 'bummer (unknown): errtyp=',
     &    errtyp, text, ierr
         call exit( 1 )
      endif
c
c     # this statement is not executed, it is included
c     # just to avoid compiler warnings. -rls
      stop 'bummer error'
c
c     # initialization...
c
      entry ibummr( iunit )
c
c     # save the listing unit for use later.
c
      nlist = iunit
c
      return
6010  format(1x,a,a,i10)
6020  format(1x,a,i10,a,i10)
      end

*deck sifcfg
      subroutine sifcfg(
     & nbft,  ifmt,    lrec,    nmax,    ierr )
c
c  return a set of consistent configuration parameters for a standard
c  integral file structure.  ( for 1-e integrals )
c
c  input:
c  nbft = total number of basis functions.
c
c  notes: (1) ifmt, lrec, and nmax may eventually have meaning on input
c         in future versions of this routine.  for now the input values
c         are ignored.
c         (2) for extensibility, the input variables may be passed
c         into this routine as array elements, instead of scalars, in
c         future versions of this routine.
c
c  output:
c  ifmt = ifmt parameter.
c  lrec = actual buffer length to be written.
c  nmax = maximum number of values in each record.
c  ierr = error return code.
c       =  0 for normal return.
c       = -1 for itype error.
c       = -2 for nbft error.
c       = -3 for lrec error.
c
      implicit none
      integer nbft, ifmt, lrec, nmax
      integer nbig, chunk, nword, space

#include "argos_error.h"
c
c     # nbftmx = largest nbft consistent with label packing methods.
c     # lrecmn = minimum allowed record length.
c     # lrecmx = maximum allowed record length.  this should be
c     #          consistent with the record "dword" encoding.
c     # lrcinc = record length block size increments.
c     # ndeflt = default number of block increments. such that
c     #          the default record length = (lrecmn + ndeflt * lrcinc)
c  
      integer    nbftmx
      parameter( nbftmx=2**16-1 )
      integer    lrecmn,     lrecmx,         lrcinc,      ndeflt
      parameter( lrecmn=512, lrecmx=2**16,   lrcinc=512,  ndeflt=3 )
c
c     # check the input parameters...
c
      ierr = 0
c
      if ( nbft .le. 0    .or.    nbft .gt. nbftmx ) then
         ierr = -2
         return
      endif
c
c     # adjust lrec for efficient I/O.
c     # machine-dependent quirks, magic record lengths, etc. should
c     # be localized here and in the above parameter definitions..
c
      lrec = lrecmn + (ndeflt * lrcinc)

      lrec = min( lrecmx, lrec )
      if ( lrec .lt. lrecmn ) then
         ierr = -3
         return
      endif
c
c     # for this record length, compute the number of values.
c
c     # 1-e records: dword // values // fcore // packed_labels // ibitv
c     # 2-e records: dword // values //          packed_labels // ibitv
c
      if ( nbft .lt. 2**8 ) then
c        # use standard 8-bit orbital label packing.
         ifmt=0
      else
c        # use standard 16-bit orbital label packing.
         ifmt=1
      endif
c
c     # compute a reasonable upper bound to nmax.
c     # chunk is used to constrain nmax so that ulab8() or ulab16() do
c     # not overwrite the labels(*) array.
c     #
c     # another choice would be to require the programmer to account
c     # for these overruns when allocating labels(*).  this would
c     # result in more efficient buffer(*) use, but the chunk choice
c     # is simpler.
c
c     # note that the programmer must explicitly account for ibitv(*)
c     # overruns; constraining chunk to 64 is too wasteful. -rls
c
      if ( ifmt .eq. 0 ) then
c        # n + n/4 <= lrec
         nbig  = ( 4 * lrec ) / 5
         chunk = 4
      elseif ( ifmt .eq. 1 ) then
c        # n + n/2 <= lrec
         nbig  = ( 2 * lrec ) / 3
         chunk = 2
      endif
c
c     # round up to a higher multiple of chunk.
      nbig = ( ( nbig + 2*chunk ) / chunk ) * chunk
c
c     # account for dword and fcore.
      nword = 2
c
c     # loop backwards from the upper bound until a valid nmax is found.
      do 10 nmax = nbig, 0, -chunk
         space = nword + nmax + ( ( nmax + chunk - 1) / chunk ) 
         if ( space .le. lrec ) goto 11
 10   continue
 11   continue
      if ( nmax .le. 0 ) then
         ierr = -4
         return
      endif
c
c     # ifmt, lrec, and nmax are all ok.
c
      return
      end

*deck sifpre
      subroutine sifpre( nlist, nenrgy, energy, ietype )
c
c  print out the energy(*) array with character labels.
c
c  24-apr-92 cnvginf values added. -rls
c  08-oct-90 (columbus day) written by ron shepard.
c
      implicit none
c
      integer  nlist,  nenrgy, ietype(*)
      real*8   energy(*)
c
      integer  i,      itypea, itypeb
      character*8      chrtyp
c
      do 10 i = 1, nenrgy
c
         itypea = ietype(i) / 1024
         itypeb = mod( ietype(i), 1024 )
         call siftyp( itypea, itypeb, chrtyp )
c
         if ( ietype(i) .ge. 0 ) then
c
c           # fcore from a 1-e hamiltonian array.
c
            write(nlist,6010) i, energy(i), ietype(i), 'fcore', chrtyp
c
         elseif ( itypea .eq. 0 ) then
c
c           # core energy value.
c
            write(nlist,6010) i, energy(i), ietype(i), 'core', chrtyp
c
         elseif ( itypea .eq. -1 ) then
c
c           # total energy value.
c
            write(nlist,6010) i, energy(i), ietype(i), 'total', chrtyp
c
         elseif ( itypea .eq. -2 ) then
c
c           # energy or wave function convergence value.
c
            write(nlist,6010) i, energy(i), ietype(i), 'cnvginf', chrtyp
c
         else
c
c           # undefined energy type.
c
            write(nlist,6010) i, energy(i), ietype(i), 'unknown', chrtyp
c
         endif
 10   continue
      return
c
 6010 format(' energy(', i2, ')=', 1pe20.12, ', ietype=', i5,
     & ', ', a7, ' energy of type: ', a )
      end

 

*deck sifsce
      function sifsce( nenrgy, energy, ietype )
c
c  select and sum the core energies in the energy(*) array.
c
c  usage: real*8   sifsce
c         external sifsce
c         total_core = sifsce(...)
c
c  input: nenrgy = number of energy(*) values.
c         ietype(1:nenrgy) = energy types.
c         energy(1:nenrgy) = energy array.
c
c  output: sifsce = total core energy such that
c                       total_potential = total_electronic + sifsce(...)
c                   is the total clamped-nucleus, born-oppenheimer
c                   potential.
c
      implicit none
c
      integer  nenrgy, ietype(nenrgy)
      real*8   sifsce, energy(nenrgy)
c
      integer  i,      itypea
      real*8   ecore
c
      ecore = (0)
      do i = 1, nenrgy
c
c        # 0 <= ietype elements are 1-e hamiltonian frozen core terms.
c        # ietype < 0  elements with itypea=0 are other core terms to
c                      be added.  e.g. ietype=-1 = nuclear repulsion.
c
         itypea = ietype(i) / 1024
         if ( (ietype(i) .ge. 0 ) .or.
     &    ( (ietype(i) .lt. 0) .and. (itypea .eq. 0) ) ) then
            ecore = ecore + energy(i)
         endif
      enddo   
c
      sifsce = ecore
      return
      end


*deck siftyp
      subroutine siftyp( itypea, itypeb, chrtyp )
c
c  return a character description of the integral type.
c
c  input:
c  itypea, itypeb = generic and specific integral or energy(*) types.
c
c  output:
c  chrtyp = character description. (this should be at least character*8
c           in the calling program.)
c
c  08-oct-90 total energy and convergence types added. -rls
c  01-jun-90 table version. -rls
c  11-aug-89 written by ron shepard.
c
      implicit none
c
c     # dummy...
      integer          itypea,  itypeb
      character*(*)    chrtyp
c
c     # local...
      integer  i
c
      integer  typ1e,  st1e,   end1e
      integer  typ2e,  st2e,   end2e
      integer  typc,   stc,    endc
      integer  typte,  stte,   endte
      integer  typcv,  stcv,   endcv
c
c     # typ1e = number of defined 1-e array types.
c     # typ2e = number of defined 2-e array types.
c     # typc  = number of core energy types.
c     # typte = number of total energy types.
c     # typcv = number of convergence types.
c
      parameter( typ1e = 16, st1e = 1,        end1e = typ1e        )
      parameter( typ2e =  2, st2e = end1e +1, end2e = end1e +typ2e )
      parameter( typc  =  1, stc  = end2e +1, endc  = end2e +typc  )
      parameter( typte = 14, stte = endc  +1, endte = endc  +typte )
      parameter( typcv = 10, stcv = endte +1, endcv = endte +typcv )
c
      integer    ntype
      parameter( ntype = endcv )
c
      integer ltypea(ntype), ltypeb(ntype)
      character*8 lctype(ntype)
c
c     # these array types can be in any order, so new ones can be added
c     # to the end or inserted into the middle as appropriate. -rls
c
c     # warning: case dependent code:  do not change the case of the
c     # following character strings.
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = st1e, end1e ) /
     & 0,    0, 'S1(*)',     0,    1, 'T1(*)',
     & 0,    2, 'V1(*)',     0,    3, 'Veff(*)',
     & 0,    4, 'VFC(*)',    0,    5, 'Vref(*)',
     & 0,    6, 'H1(*)',     0,    7, 'D1(*)',
     & 0,    8, 'F(*)',      0,    9, 'Q(*)',
     & 1,    0, 'X(*)',      1,    1, 'Y(*)',
     & 1,    2, 'Z(*)',
     & 2,    0, 'Im(SO:x)',  2,    1, 'Im(SO:y)',
     & 2,    2, 'Im(SO:z)'  /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = st2e, end2e ) /
     & 3,    0, '1/r12',     3,    1, 'd2(*)'     /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stc, endc ) /
     & 0,   -1, 'Nuc.Rep.' /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stte, endte ) /
     & -1,   0, 'SCF',      -1,   -1, 'MCSCF',
     & -1,  -2, 'MRSDCI',   -1,   -3, 'CPF',
     & -1,  -4, 'ACPF',     -1,   -5, 'LCC-SD',
     & -1,  -6, 'MRPT',     -1,   -7, 'Bk',
     & -1,  -8, 'DV1',      -1,   -9, 'DV2',
     & -1, -10, 'EPOPLE',   -1,  -11, 'S.O. CI',
     & -1, -12, 'SR-SDCI',  -1,  -13, 'UCEPA'  /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stcv, endcv ) /
     & -2,   0, 'SCF-D.E.', -2,   -1, 'SCF-D.D1',
     & -2,  -2, 'MC-D.E.',  -2,   -3, 'MC-Wnorm',
     & -2,  -4, 'MC-Knorm', -2,   -5, 'MC-ApxDE',
     & -2,  -6, 'Bk-Resid', -2,   -7, 'CI-Resid',
     & -2,  -8, 'CI-D.E.',  -2,   -9, 'CI-ApxDE' /
c
      do 10 i = 1, ntype
         if ( itypea .eq. ltypea(i) ) then
            if ( itypeb .eq. ltypeb(i) ) then
               chrtyp = lctype(i)
               return
            endif
         endif
 10   continue
c     # loop exit means unrecognized type.
      chrtyp = 'Unknown'
c
      return
      end



      subroutine setiarray(
     &   itypea, itypax, itypbx, btypmx, btypes, iarray)
!
!     # set iarray to point into array(*,iarray).
!
      integer itypea, itypax, itypbx, btypmx
      integer btypes(0:btypmx), iarray

         if ( itypea .eq. -1 ) then
!           # single array mode. accumulate into array(*,1) always.
            iarray = 1
         elseif ( itypax.eq.0 ) then
!           # correct generic array type for s, t, v, c.
!           # mask btypmx with btypes(*) to see if this array should
!           # be accumulated somewhere.
            if ( itypbx .le. btypmx) then
               iarray = btypes(itypbx)
            else
               iarray = -1
            endif
         elseif ( itypax.eq.2 ) then
!           # correct generic array type for soc.
!           # mask btypmx with btypes(*) to see if this array should
!           # be accumulated somewhere.
            if ( itypbx .le. btypmx) then
               iarray = 4 + btypes(itypbx)
            else
               iarray = -1
            endif
         else
!           # ignore this array.
            iarray = -1
         endif
      return
      end

*deck dcopy
*mdc*if ( cray .and. inline) or fujitsu
      subroutine  dcopy_argos(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses rolled loops for increments equal to one.
c     written by ron shepard, based on dcopy written by
c     jack dongarra, linpack, 3/11/78.
c
      real*8 dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if ( n.le.0 ) return
      if ( incx.eq.1 .and. incy.eq.1 ) go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if ( incx.lt.0 ) ix = (-n+1)*incx + 1
      if ( incy.lt.0 ) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         dy(iy) = dx(ix)
         ix = ix + incx
         iy = iy + incy
10    continue
      return
c
c        code for both increments equal to 1
c
20    continue
      do 50 i = 1, n
         dy(i) = dx(i)
50    continue
      return
      end

*deck wzero
      subroutine wzero(n,a,inca)
c
      implicit none
      integer n, inca
c
      integer ia, i
      real*8    zero
      parameter(zero=0.d0)

      real*8 a(*)
      ia = 1
      do 10 i = 1, n
         a(ia) = zero
         ia    = ia + inca
10    continue
c
      return
      end
