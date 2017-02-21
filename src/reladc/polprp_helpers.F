!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKEINVERSE(inverse)
      IMPLICIT none

      integer, dimension(:)            :: inverse

#include  "../relccsd/symm.inc"  !provides nrep,multb

      integer                          :: i,j

      write(*,*) 'lbound inverse:',lbound(inverse)
      write(*,*) 'ubound inverse:',ubound(inverse)
      do i=1,nrep
        inverse(i) = multb(i,nrep+1,2)
      enddo
      do i=1,nrep
        inverse(i+nrep) = multb(i+nrep,nrep+1,2)
      enddo

      Return
      End
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE property_handler(ioch,ioch2,lavail,nsp,irepspi)

      use polprp_cfg   !make some input variables available

      IMPLICIT none
!
! this routine reads in the MO-transformed
! property integrals from MDPROP via RDPROP_ADC maintaining
! complex arrays throughout.
! RDPROP_ADC is a modified version of the old relccsd/rdprop routine.
! The property integrals are analyzed for the CARITH=F cases where we still
! can have complex wave functions.
!
!---------------Calling variables--------------------------------------
!
      integer                          ::  ioch,ioch2  !file handles
      logical                          ::  lavail      !flag for property integrals
      integer                          ::  nsp         !number of active spinors
      integer, dimension(:)            ::  irepspi     !spinor translation array
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"     !provides MXREP,NXREP etc.
#include  "../relccsd/complex.inc"  !provides CARITH,RCW
!
!---------------Local variables--------------------------------------
!
      logical      :: existing,success
      integer      :: npaoo, npaov, npavo, npavv  !specific prop array dimensions

      integer      :: i,iprop
      integer      :: MDPROP

      character*8  :: nameofprop
      real*8, allocatable   :: CBUF(:)
      real*8, allocatable   :: CPVO(:), CPOV(:), CPOO(:), CPVV(:)
      complex*16            :: cphase(4)
      character*8,dimension(3)   :: naoplist
!
!---------------Interface region--------------------------------------
!
      interface

!  in order to be able to access the old F77 routine correctly
!  we do not provide an interface to this routine. Then the 
!  assumed size array transfer should work.

!       SUBROUTINE rdprop_adc(i1,c1,i2,l1,ra1,ia1
!         integer                  :: i1,i2,i3
!         logical                  :: l1
!         integer, dimension(:)    :: ia1
!       END SUBROUTINE

        SUBROUTINE propanaly(ch,l1,l2,ca1,
     &                       i1,i2,i3,i4,
     &                       ra1,ra2,ra3,ra4)
          character*8, intent(in)                :: ch
          logical, intent(in)                    :: l1,l2
          complex*16, dimension(:)               :: ca1
          integer                                :: i1,i2,i3,i4
          real*8, dimension(:)                   :: ra1,ra2,ra3,ra4
        END SUBROUTINE

        SUBROUTINE propflush(a1,i1,ca1,
     &                       i2,i3,i4,i5,ra1,ra2,ra3,ra4)
          character*1                      :: a1
          integer                          :: i1
          complex*16,dimension(:)          :: ca1
          integer                          :: i2,i3,i4,i5
          real*8, dimension(:)             :: ra1,ra2,ra3,ra4
        END SUBROUTINE

      end interface
!
!---------------Executable code --------------------------------------
!
      call PST('Entering property handler+')
!
!**  check for availability of MDPROP
!
      MDPROP = 15
      lavail = .false.
      INQUIRE(FILE='MDPROP',EXIST=existing)
      IF (.NOT.existing) THEN
         write(*,*) '  ** Warning ***'
         write(*,*) '  ** Property integrals (MDPROP) do not exist.'
         write(*,*) '  ** No transition moments are calculated.'
         write(*,*) '  ** Continuing with energy-only calculations.'
         RETURN
      ELSE
        lavail=.true.
        write(*,*) '  ** MDPROP file found.'
      END IF
!     write(*,*)
!     write(*,*) 'Print level is set to ',polprp_printlev
!     write(*,*)
!
! allocate proper space for the fou kinds of property integrals.
!
      npaoo = 0
      npaov = 0
      npavo = 0
      npavv = 0
      do i=1,nrep
        npaoo = npaoo + MOO(i)
        npaov = npaov + MOV(i)
        npavo = npavo + MVO(i)
        npavv = npavv + MVV(i)
      enddo
!
! check dimensions for consistency
!
      if( (npaoo+npaov+npavo+npavv) .ne. (nsp*nsp)) then
        write(*,*) '  ** Consistency check failed with:'
        write(*,*) '  ** property array has length',nsp*nsp
        write(*,*) '  ** Calculated dims are:',
     &                npaoo+npaov+npavo+npavv
        STOP
!     else
!       write(*,*) '  ** Consistency check passed with:',nsp*nsp
!       write(*,*) '  ** Size of OO array:',npaoo
!       write(*,*) '  ** Size of OV array:',npaov
!       write(*,*) '  ** Size of VO array:',npavo
!       write(*,*) '  ** Size of VV array:',npavv
      endif
!
! allocate proper space for property integrals
!
      allocate(CBUF(2*NSP*NSP))   !Always complex, holding raw integrals
      allocate(CPOO(RCW*npaoo))
      allocate(CPOV(RCW*npaov))
      allocate(CPVO(RCW*npavo))
      allocate(CPVV(RCW*npavv))
!
! set list of property names required for the TMs
!
      naoplist(1) = 'XDIPLEN'
      naoplist(2) = 'YDIPLEN'
      naoplist(3) = 'ZDIPLEN'
!
! loop over the three properties and store the integrals in the
! corresponding files.
!

      Do iprop = 1,3

        nameofprop=naoplist(iprop)

        success = .false.

        CALL read_flush_prop(success,nameofprop,NSP,IREPSPI,
     &                       npaoo,npaov,npavo,npavv,
     &                       CBUF,CPOO,CPOV,CPVO,CPVV,
     &                       polprp_printlev)

        if(.not.success) then
          call insuff_prop_return(nameofprop)
          lavail = .false.
          return
        endif

      Enddo
!
! releasing heap memory
!
      deallocate(CPVV)
      deallocate(CPVO)
      deallocate(CPOV)
      deallocate(CPOO)
      deallocate(CBUF)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE insuff_prop_return(nameofprop)
      character*8  :: nameofprop
      write(*,*)
      write(*,*) '  *** WARNING ***'
      write(*,*)
      write(*,*) 'Property file MDPROP available, but required dipole'
      write(*,*) 'moment properties are not available there.'
      write(*,*)
      write(*,*) 'Please add the following lines to your'
      write(*,*) '**MOLTRA section, in order to create the required'
      write(*,*) 'operator names for property processing:'
      write(*,*)
      write(*,*) '**MOLTRA'
      write(*,*) '.ACTIVE'
      write(*,*) 'energy ...  ...   ...'
      write(*,*) '.PRPTRA'
      write(*,*) '*PRPTRA'
      write(*,*) '.OPERATOR'
      write(*,*) ' ''X_DIPLEN'' '
      write(*,*) ' DIAGONAL'
      write(*,*) ' XDIPLEN'
      write(*,*) '.OPERATOR'
      write(*,*) ' ''Y_DIPLEN'' '
      write(*,*) ' DIAGONAL'
      write(*,*) ' YDIPLEN'
      write(*,*) '.OPERATOR'
      write(*,*) ' ''Z_DIPLEN'' '
      write(*,*) ' DIAGONAL'
      write(*,*) ' ZDIPLEN'
      write(*,*) ' COMFACTOR'
      write(*,*) ' 1.0000000000'
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) 'Calculation continues without transition moments.'
      write(*,*)
      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PROPFETCH(ioch,pph,
     &                     idoo,idov,idvo,idvv,
     &                     doox,dooy,dooz,
     &                     dovx,dovy,dovz,
     &                     dvox,dvoy,dvoz,
     &                     dvvx,dvvy,dvvz)
     &               
 
      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!     Reads all property matrices from disk
!
!---------------Calling variables--------------------------------------
!
      Integer                    :: ioch
      Complex*16,dimension(:)    :: pph
      Integer                    :: idoo,idov,idvo,idvv
      Real*8,dimension(:)        :: doox,dooy,dooz
      Real*8,dimension(:)        :: dovx,dovy,dovz
      Real*8,dimension(:)        :: dvox,dvoy,dvoz
      Real*8,dimension(:)        :: dvvx,dvvy,dvvz
!
!---------------Common blocks----------------------------------------
!
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      Character*1              :: co
      Character*8              :: filename
      Character*7              :: nameroot
      Integer                  :: icaoo,icaov,icavo,icavv,ixx,ipix,i
      Integer                  :: ist
      Complex*16               :: A0,A1,AI
      Real*8                   :: ur,ui
!
!---------------Executable code--------------------------------------
!
      A0=(0.0d0,0.0d0)
      A1=(1.0d0,0.0d0)
      AI=(0.0d0,1.0d0)
      nameroot = 'PROPINT'

!  ********** X

      co='X'
      ipix=1
      write(filename,'(A1,A7)') co,nameroot
      open(unit=ioch,file=filename,form='unformatted',
     &     access='sequential',status='unknown')
      read(ioch) icaoo,icaov,icavo,icavv
      if( (  (icaoo - idoo) +
     &       (icaov - idov) +
     &       (icavo - idvo) +
     &       (icavv - idvv) ).ne.0) then
        write(*,*) 'Dimension inconsistency in ',filename
        call quit('***')
      endif
      read(ioch) ur,ui
      read(ioch) (doox(ixx),ixx=1,rcw*idoo)
      read(ioch) (dovx(ixx),ixx=1,rcw*idov)
      read(ioch) (dvox(ixx),ixx=1,rcw*idvo)
      read(ioch) (dvvx(ixx),ixx=1,rcw*idvv)
      close(ioch)
      pph(ipix) = dcmplx(ur,ui)
      write(*,*) ' ** Phase for ',filename,'  ',pph(ipix)
      If(carith) then
          if(pph(ipix).ne.A0) then
            write(*,*) 'Phase inconsistency in ',filename
          endif
      endif
      write(*,*) ' ** Integrals from ',filename,' successfully read.'


!  ********** Y

      co='Y'
      ipix=2
      write(filename,'(A1,A7)') co,nameroot
      open(unit=ioch,file=filename,form='unformatted',
     &     access='sequential',status='unknown')
      read(ioch) icaoo,icaov,icavo,icavv
      if( (  (icaoo - idoo) +
     &       (icaov - idov) +
     &       (icavo - idvo) +
     &       (icavv - idvv) ).ne.0) then
        write(*,*) 'Dimension inconsistency in ',filename
        call quit('***')
      endif
      read(ioch) ur,ui
      read(ioch) (dooy(ixx),ixx=1,rcw*idoo)
      read(ioch) (dovy(ixx),ixx=1,rcw*idov)
      read(ioch) (dvoy(ixx),ixx=1,rcw*idvo)
      read(ioch) (dvvy(ixx),ixx=1,rcw*idvv)
      close(ioch)
      pph(ipix) = dcmplx(ur,ui)
      write(*,*) ' ** Phase for ',filename,'  ',pph(ipix)
      If(carith) then
          if(pph(ipix).ne.A0) then
            write(*,*) 'Phase inconsistency in ',filename
          endif
      endif
      write(*,*) ' ** Integrals from ',filename,' successfully read.'

!  ********** Z

      co='Z'
      ipix=3
      write(filename,'(A1,A7)') co,nameroot
      open(unit=ioch,file=filename,form='unformatted',
     &     access='sequential',status='unknown')
      read(ioch) icaoo,icaov,icavo,icavv
      if( (  (icaoo - idoo) +
     &       (icaov - idov) +
     &       (icavo - idvo) +
     &       (icavv - idvv) ).ne.0) then
        write(*,*) 'Dimension inconsistency in ',filename
        call quit('***')
      endif
      read(ioch) ur,ui
      read(ioch) (dooz(ixx),ixx=1,rcw*idoo)
      read(ioch) (dovz(ixx),ixx=1,rcw*idov)
      read(ioch) (dvoz(ixx),ixx=1,rcw*idvo)
      read(ioch) (dvvz(ixx),ixx=1,rcw*idvv)
      close(ioch)
      pph(ipix) = dcmplx(ur,ui)
      write(*,*) ' ** Phase for ',filename,'  ',pph(ipix)
      If(carith) then
          if(pph(ipix).ne.A0) then
            write(*,*) 'Phase inconsistency in ',filename
          endif
      endif
      write(*,*) ' ** Integrals from ',filename,' successfully read.'

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMVOTS(eps,aak)
 
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Divides the totally symmetric A_ak array by eps_a - eps_k
!
!---------------Calling variables--------------------------------------
!
      REAL*8,dimension(:)              :: EPS
      REAL*8,dimension(:)              :: AAK
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8                   :: epsa,epsk,ediff
C
C---------------Executable code--------------------------------------
C
      AK = 0
      DO KREP = 1, NREP
        AREP = KREP    ! array is totally symmetric
        KOFF = IO(KREP)
        AOFF = IO(NREP+1) + IV(AREP)
        DO K = 1,NO(KREP)
          DO A = 1,NV(AREP)
            epsa = EPS(AOFF + A)
            epsk = EPS(KOFF + K)
            ediff = epsa - epsk
            AK = AK+1
            AAK(AK) = AAK(AK)/ediff
!           write(*,*) 'DENOMVOTS: ak,aak(ak):',ak,aak(ak)
          ENDDO
        ENDDO
      ENDDO
      IF(AK.ne.MVO(1)) STOP 'Inconsistency in DENOMVOTS!'

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE rcaxpy(n,carith,phase,as,ad)
      Implicit none
!
!---------------- description ----------------
!
! this routine xaxpys from as to ad depending on phase
! if phase is 1.0,0.0 then ad = ad + as
! if phase is 0.0,1.0 then ad = ad + I*as
! ad is always complex, as is either real or complex (following carith)
!
!---------------- formal parameters ----------------
!
      integer                      :: n
      logical                      :: carith
      complex*16                   :: phase
      real*8,dimension(:)          :: as
      complex*16,dimension(:)      :: ad
!
!---------------- local variables ----------------
!
      integer                      :: rcw
      integer                      :: i
!
!---------------- executable code ----------------
!
      rcw = 1
      if(carith) rcw = 2

!     write(*,*) ' *** in rcaxpy ***',carith,rcw

      if(.not.carith) then       ! ** as = real case + phase **
        do i=1,n
          ad(i) = ad(i) + phase*dcmplx(as(i),0.0d0)
        enddo
      else                       ! ** as = complex case, no phase **
        do i=1,n
          ad(i) = ad(i) + dcmplx(as(2*i-1),as(2*i))
        enddo
      endif

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RCA_NEGATE(n,arr)
!
!  this routine changes sign of a real or a complex array
!  depending on the CARITH flag. In case of complex:
!  z -->  -z  not  z --> z*  !
!
      implicit none
      integer                      :: n
      real*8, dimension(:)         :: arr
!
! it is assumed that n counts the elements and not the array length!
!

#include  "../relccsd/complex.inc"

      integer                      :: i

      DO I=1,N
        IF(CARITH) THEN
          arr(2*i-1) = -1.0D0*arr(2*i-1)
          arr(2*i)   = -1.0D0*arr(2*i)
        ELSE
          arr(i) = -1.0D0*arr(i)
        ENDIF
      ENDDO
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TM_GEN(ioch,desrep,isatctl,iorder,
     &                  do2h2p,
     &                  faitx,faity,faitz,
     &                  fabijtx,fabijty,fabijtz)
      implicit none
!
!  this routine contracts the Davidson eigenvector with the 
!  F_ai/F_aibj matrix elements to form the three transition moments.
!  Remind that the F_ai and F_aibj elements together with the TMs
!  are kept complex!
!
!  due to high numerical effort in the calling routine a flag (do2h2p) 
!  determines if the 2h2p contributions should be included. In this
!  case the arrays fabijt[xyz] must be allocated! Otherwise they are not
!  referenced.
!
!  The routine can be called for each order of perturbation theory in the
!  F_ai matrix elements. It just needs access to the Davidson eigenvectors.
!
!  calling parameters
!
      integer                      :: ioch,desrep,isatctl
      character*30                 :: iorder
      logical                      :: do2h2p
      complex*16, dimension(:)     :: faitx,faity,faitz
      complex*16, dimension(:)     :: fabijtx,fabijty,fabijtz

!
!  includes
!
#include  "../relccsd/complex.inc"
!
!  locals
!
      character*8, parameter                  :: basename='DAVILONG'
      character*11                            :: davevcname
      logical                                 :: isthere

      real*8, allocatable, dimension(:)       :: rcoef
      complex*16, allocatable, dimension(:)   :: ccoef
      complex*16,allocatable,dimension(:)     :: tmx,tmy,tmz

      integer                                 :: ioresult
      integer                                 :: nroots,nmain,ladc
      real*8,allocatable,dimension(:)         :: evl,xnr
      real*8                                  :: ps
      real*8, parameter                       :: D0=0.0d0
      complex*16, parameter                   :: A0=(D0,D0)
      integer                                 :: i,k,ist
      complex*16                              :: ax,ay,az
      real*8, parameter                       :: autoev = 27.2113957D0
      integer                                 :: n2h2p
      character*7,dimension(2)                :: numtype
!
!  executable
!
      if(carith.and.RCW.ne.2) STOP 'Carith error'
      if(.not.carith.and.RCW.ne.1) STOP 'Carith error'
      numtype(1) = 'real'
      numtype(2) = 'complex'
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') basename,'.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') basename,'.0',desrep
      ENDIF
      INQUIRE(file=davevcname,exist=isthere)
      if(.not.isthere) then
        write(*,*) '***************************************'
        write(*,*) '*****  Warning:                   *****'
        write(*,*) '***************************************'
        write(*,*) 'No eigenvectors found in symmetry',desrep,'!'
        write(*,*) 'No TMs accessible!'
        write(*,*)
        return
      else
        write(*,*)
        write(*,*) 'Found ',numtype(rcw),' eigenvecs in symm.',desrep
      endif
      OPEN(ioch,file=davevcname,access='sequential',
     &     status='unknown',form='unformatted')
      REWIND(ioch)

!_________________________  perform the contractions
!|

      READ(ioch,iostat=ioresult) nroots
      READ(ioch,iostat=ioresult) nmain
      READ(ioch,iostat=ioresult) ladc

      write(*,*) 'Found',nroots,' states in the eigenvector file.'
      write(*,*) 'Length of ADC matrix (from file):',ladc
      write(*,*) 'Length of main space (from file):',nmain
      n2h2p = ladc-nmain
      write(*,*) 'Length of sat space:',n2h2p

      if(do2h2p) then
        write(*,*) 'Including 2h2p contributions to TM.'
        if (isatctl.ne.n2h2p) 
     &    call quit('Internal error for coupling block length')
      else
        write(*,*) 'Omitting 2h2p contributions to TM'
      endif
!
! allocate eigenvalue, transition moment and coefficient arrays
!
      allocate(evl(nroots))
      evl = D0

      allocate(tmx(nroots))
      allocate(tmy(nroots))
      allocate(tmz(nroots))
      tmx = A0
      tmy = A0
      tmz = A0
      if(carith) then
        allocate(ccoef(ladc))
      else
        allocate(rcoef(ladc))
      endif

      allocate(xnr(nroots))
      xnr = D0

      Do k=1,nroots
        read(ioch) ist
        read(ioch) evl(k)
        read(ioch) ps
        write(*,'(A,X,I3,X,A,F10.4,A,X,F10.4,A,F10.4)') 
     &    'State',ist,':exc. energy',evl(k),'(au)',evl(k)*autoev,
     &    ' (eV) and PS',ps
        if(carith) then
          read(ioch)(ccoef(i),i=1,ladc)
          do i=1,ladc
            xnr(k)=xnr(k) + 
     &         (real(ccoef(i))**2.0d0 + aimag(ccoef(i))**2.0d0)
          enddo
        else
          read(ioch)(rcoef(i),i=1,ladc)
          do i=1,ladc
            xnr(k)=xnr(k) + rcoef(i)**2.0d0
          enddo
        endif
        xnr(k) = dsqrt(xnr(k))
        
!
! the do2h2p flag controls the level of contraction.
!
        ax=A0
        ay=A0
        az=A0
        if(carith) then
          do i=1,nmain
            ax = ax + ccoef(i) * faitx(i)
            ay = ay + ccoef(i) * faity(i)
            az = az + ccoef(i) * faitz(i)
          enddo
          if(do2h2p) then
            do i=1,n2h2p
              ax = ax + ccoef(i+nmain) * fabijtx(i)
              ay = ay + ccoef(i+nmain) * fabijty(i)
              az = az + ccoef(i+nmain) * fabijtz(i)
            enddo
          endif
        else
          do i=1,nmain
            ax = ax + dcmplx(rcoef(i),D0)*faitx(i)
            ay = ay + dcmplx(rcoef(i),D0)*faity(i)
            az = az + dcmplx(rcoef(i),D0)*faitz(i)
!           write(*,'(F10.6,5X,2F10.6,5X,2F10.6)') rcoef(i),faitx(i),ax
          enddo
          if(do2h2p) then
            do i=1,n2h2p
              ax = ax + dcmplx(rcoef(i+nmain),D0)*fabijtx(i)
              ay = ay + dcmplx(rcoef(i+nmain),D0)*fabijty(i)
              az = az + dcmplx(rcoef(i+nmain),D0)*fabijtz(i)
            enddo
          endif
        endif
        tmx(k) = ax
        tmy(k) = ay
        tmz(k) = az

      Enddo   ! nroots

      CLOSE(ioch)
!|________________________ end reading section

!_________________________  print output
!|
!|
      write(*,*) '   Transition moments for states in symmetry',desrep
      write(*,*)
      write(*,*)
      write(*,*) '   *******************************'
      write(*,'(A,A)') '   Order of computation: ',iorder
      write(*,*) '   *******************************'
      write(*,*)
      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,'(A,4X,A,9X,A,20X,A,20X,A)')
     &       'State','Energy(eV)','TM(x)','TM(y)','TM(z)'
      write(*,*) '----------------------------------------------------'
      write(*,*)
      do k=1,nroots
        write(*,'(I4,F14.7,F8.4,2F15.9,A,2F15.9,A,2F15.9)') k,
     &     evl(k)*autoev,xnr(k),tmx(k),' | ',tmy(k),' | ',tmz(k)
      enddo
      write(*,*) '----------------------------------------------------'
!|
!|________________________ end print output

      if(carith) then
        deallocate(ccoef)
      else
        deallocate(rcoef)
      endif
      deallocate(tmx)
      deallocate(tmy)
      deallocate(tmz)

      deallocate(evl)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE read_flush_prop(SUCC,NAMEA,NSP,IREPSPI,
     &                           dimpoo,dimpov,dimpvo,dimpvv,
     &                           PROP,AOO,AOV,AVO,AVV,
     &                           printlevel)
!
      use interface_to_mpi
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!
!     Reads property integral files for operator A.
!     The routine was originally devised by LV and modified/extended
!     by MP for the RELADC purposes.
!     This version reads from the file created by PRTRAN
!
!     Important background for the use of MO-transformed property integrals:
!
!     Due to the possibility of spinfree calculations the wave function can
!     still be complex leading to complex property integrals also in the case
!     of CARITH=false. CARITH does not indicate if a calculation is quaternionic
!     or not but tells arithmetic only for calculations after the MOLTRA step.
!     A Spinfree calculation therefore leads to complex MOPIs even if CARITH=
!     false! In this case we can get purely imaginary property integrals for the
!     x,y,a dipole moment. Since carith=false, the integrals end up in the AXX
!     arrays which causes AOV = -AVO **seemingly violating Hermiticity**
!
!     We therefore need a canary indicating if we have purely imaginary integrals
!     The real array is therefore to be interpreted as the imaginary part of the
!     property integral in the case CARITH=F. The APHASE array transports the corresponding
!     phase information (purely real or purely imaginary) to the calling routine.
!
!     This allows for sticking at real arrays in the case of CARITH=F and postmultiplying
!     the F_ai matrix elements by the imaginary unit. The F_ai's therefore have to be kept
!     complex also in the SPINFREE/NONREL cases.
!
!     ATT: the absolute spinor numbers obtained through the irepspi array are *not*
!     the same ones as listed in the list of active occ/virt orbitals at the beginning
!     of CCSETI output! The irepspi absolute indices serve as pointers into the 
!     PROP array where they correspond to the proper indexing.
!
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/param.inc"
#include "../relccsd/symm.inc"
#include "../relccsd/complex.inc"
#include "../relccsd/ccpar.inc"
!
!---------------Calling variables--------------------------------------
!
      Logical                    SUCC
      Character*8                NAMEA
      INTEGER                    NSP
      INTEGER                    IREPSPI(NSP,MXREP,2)
      INTEGER                    dimpoo,dimpov,dimpvo,dimpvv
      REAL*8                     PROP(2,NSP,NSP)
      REAL*8                     AOO(*),AOV(*),AVO(*),AVV(*)
      integer                    printlevel
!
!---------------Local variables--------------------------------------
!
      Integer                 :: MDPROP, PROCPROP
      REAL*8,allocatable      :: DUMPROP(:,:,:)
      CHARACTER*32            :: ACHAR
      REAL*8                  :: APHASE(2)
      REAL*8                  :: AVOI,AOVI,AMAXR,AMAXI,ANORM
      LOGICAL                 :: PHASE
      REAL*8                  :: hdev
      COMPLEX*16              :: h1,h2
! --- flushing
      Character*1              :: co
      Character*8              :: filename
      Integer                  :: ixx,rawnum
      real*8, parameter        :: ZEROTHRS=1.0D-11
      integer,allocatable      :: acv(:,:)
      integer                  :: irow,icol,retcode
!
!---------------Executable code--------------------------------------
!
!
!  debugging:  check irepspi organization
!
!       WRITE(*,*)
!       WRITE(*,*) '|| IREPSPI organization: ||'
!       WRITE(*,*)
!       WRITE(*,*) 'occupied:'
!       Do qrep=1,nrep
!         Do n1=1,NO(qrep)
!           write(*,*) 'G,i,absi:',qrep,n1,irepspi(n1,qrep,1)
!         Enddo
!       Enddo
!       WRITE(*,*) 'virtual:'
!       Do qrep=1,nrep
!         Do n1=1,NV(qrep)
!           write(*,*) 'G,a,absa:',qrep,n1,irepspi(n1,qrep,2)
!         Enddo
!       Enddo

!
!     Read symmetry information (master node only) from MDPROP
!     The file is organized in a way that after a label identifying the
!     property a binary number stream without CRLF follows. If the label
!     is no match, we read til the next CRLF which corresponds to a skipping
!     of all the unwanted integrals.
!
      MDPROP = 15
      PROCPROP = MDPROP + 1

      IF(dimpov.ne.dimpvo) CALL QUIT('Inconsistent dimensions!')

      IF (MYPROC.EQ.MASTER) THEN

      allocate(DUMPROP(NSP,NSP,2))
      DUMPROP = 0.0d0
      OPEN (MDPROP,FILE='MDPROP',FORM='UNFORMATTED')
      INOP = 0
    1 READ (MDPROP,ERR=10,END=11) ACHAR
      IF (ACHAR(1:8).NE.'********'.OR.ACHAR(25:32).NE.NAMEA) GOTO 1
      WRITE (6,1000) NAMEA
      READ (MDPROP,ERR=10,END=11) 
     &    ((DUMPROP(I,J,1),DUMPROP(I,J,2),I=1,NSP),J=1,NSP)
      PROP=0.0d0
      Do I=1,NSP
        Do J=1,NSP
          PROP(1,I,J) = DUMPROP(I,J,1)
          PROP(2,I,J) = DUMPROP(I,J,2)
        Enddo
      Enddo
      deallocate(DUMPROP)
      CLOSE (MDPROP,STATUS='KEEP')
!
! debug C call
!
!     retcode = polprp_receiver(PROP,2*NSP*NSP)
!
! end debug
!

      write (6,*)
      hdev = 0.0d0
      rawnum = 0
      do i=1,nsp
        do j=1,nsp
          if ( (dabs(PROP(1,i,j)).GT.ZEROTHRS).OR.
     &         (dabs(PROP(2,i,j)).GT.ZEROTHRS) ) THEN
            if(printlevel.ge.3) then
              write(6,1013) I,J,PROP(1,i,j),PROP(2,i,j)
            endif
            rawnum = rawnum + 1
          endif
          h1=dcmplx(prop(1,j,i),prop(2,j,i))
          h2=dcmplx(prop(1,i,j),prop(2,i,j))
          hdev = hdev + abs(h1 - dconjg(h2))
        enddo
      enddo
      write (6,*) ' *** Raw property integrals for ',NAMEA,':',
     & ' being above ',ZEROTHRS,' : ',rawnum
      write (6,*) ' *** Number of spinors (squared)',nsp,nsp*nsp
      write (6,*) ' *** Deviation from Hermiticity:',hdev
      write (6,*)

      GOTO 12
   10 INOP = 1
      GOTO 12
   11 INOP = 2
   12 CONTINUE
 1013 FORMAT (2I5,F18.10,2X,F18.10)

      ENDIF !Myproc.eq.master
!
!     Distribute the information read from MDPROP
!
!if defined (VAR_MPI)
!     IF (NMPROC .GT. 1) THEN
!        call interface_mpi_BCAST(INOP,1,MASTER,
!    &                  global_communicator)
!        call interface_mpi_BCAST(PROP,2*NSP*NSP,MASTER,
!    &                  global_communicator)
!     END IF
!endif
!
!     Error exit if the integrals could not be read
!
      IF (INOP.EQ.1) THEN
        WRITE (6,*) ' Error reading property ',NAMEA,' in MDPROP'
        WRITE (6,*) ' Error reading property integrals'
        SUCC = .false.
        return
      ENDIF
      IF (INOP.EQ.2) THEN
        WRITE (6,*) ' Property ',NAMEA,' not found in MDPROP'
        WRITE (6,*) ' premature end of file or integrals missing!'
        SUCC = .false.
        return
      ELSE
        WRITE (6,*) ' Property ',NAMEA,' available in MDPROP'
      ENDIF

      if(carith) then
        write(*,*) 'We have complex arithmetic. RCW=',RCW
      else
        write(*,*) 'We have real arithmetic. RCW=',RCW
      endif
!
!__________________________________________________________________
!|   perform array coverage test for OO/OV/VO/VV property arrays
!|   only if printlevel >= 2
!
      if(printlevel.ge.2) then
        i=0
        write(*,*) 'Index coverage test....'
        call idxcoverage(irepspi,nsp,i)
        if(i.ne.0) stop 'Incomplete array coverage!'
        write(*,*) '... passed.'
      endif
!|
!|_________________________________________________________________
!

      APHASE(1) = D1
      APHASE(2) = D0
      PHASE = .FALSE.

!     if(printlevel.ge.2) then
!       write(*,*)
!       write(*,*) 'Symmetry ordering of PROP ints:'
!       write(*,*)
!     endif

    2 AMAXI = D0
      AAMAX = 0
      IIMAX = 0
      AI = 1
      DO AIREP = 1, NREP
         DO IREP = 1, NREP
            AREP = MULTB (IREP,AIREP+NREP,2)
            AREP = MULTB (AREP,NREP+1,2)
!           if(printlevel.ge.2)
!    &        write(*,*) AIREP,' (',repna(airep),') ',arep,
!    &            ' (',repna(arep),') ',irep,' (',repna(irep),')'
            DO I = 1, NO(IREP)
               II = IREPSPI(I,IREP,1)
               DO A = 1, NV(AREP)
                  AA = IREPSPI(A,AREP,2)
                  IF (CARITH) THEN
                      AVO(RCW*AI-1) = PROP(1,AA,II)
                      AVO(RCW*AI)   = PROP(2,AA,II)
                      AOV(RCW*AI-1) = PROP(1,II,AA)   !ov acquire same order as vo !
                      AOV(RCW*AI)   = PROP(2,II,AA)
                  ELSE
                      AVO(AI) = PROP(1,AA,II)*APHASE(1) -
     &                          PROP(2,AA,II)*APHASE(2)
                      AOV(AI) = PROP(1,II,AA)*APHASE(1) -
     &                          PROP(2,II,AA)*APHASE(2)
!       if(printlevel.ge.2)
!    &    write(*,*) '-avo- irepspi#:',aa,ii,'  ccseti#:',
!    &    io(nrep+1)+iv(arep)+a,  io(irep)+i,'  value:',
!    &    avo(ai)
                      AVOI = PROP(1,AA,II)*APHASE(2) +
     &                       PROP(2,AA,II)*APHASE(1)
                      AOVI = PROP(1,II,AA)*APHASE(2) +
     &                       PROP(2,II,AA)*APHASE(1)
                      IF (ABS(AVOI).GT.ABS(AMAXI)) THEN
                         AMAXR = AVO(AI) 
                         AMAXI = AVOI
                         AAMAX = AA
                         IIMAX = II
                      ENDIF
                      IF (ABS(AOVI).GT.ABS(AMAXI)) THEN
                         AMAXR = AOV(AI) 
                         AMAXI = AOVI 
                         AAMAX = II
                         IIMAX = AA
                      ENDIF
                  ENDIF
                   AI = AI + 1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      IF (.NOT.CARITH.AND.ABS(AMAXI).GT.ACCUR) THEN
         IF (PHASE) THEN
            WRITE (IW,1010) AA,II,AMAXI
            CALL QUIT('USE COMPLEX ARITHMETICS')
         ELSE
            APHASE(1) = AMAXR
            APHASE(2) = AMAXI 
            ANORM = SQRT(AMAXR*AMAXR+AMAXI*AMAXI)
            APHASE(1) = APHASE(1) / ANORM
            APHASE(2) = APHASE(2) / ANORM
            PHASE = .TRUE.
            WRITE (6,1020) NAMEA,APHASE
            WRITE (6,*) 'Phase correction applied, rescaling integrals.'
            WRITE (6,*) 'For this specific properties.'
!
! the following jump instruction leads to a second integral parsing
! and scaling with the proper phase.
!
            GOTO 2
         ENDIF
      ENDIF
C

      IF( (AI-1).ne.dimpov) Call Quit('AOV counting error!')

      write(6,*)
      write(6,'(A)') 'From the AOV integrals we determined a phase'
      write(6,'(A,2F10.5)') 'of ',APHASE(1),APHASE(2)
      write(6,'(A)') 'which is also applied to the AOO and AVV ints.'
      write(6,*)
!
!  generate AOO stream
!
      IJ = 1
      DO IJREP = 1,NREP
        DO JREP = 1, NREP
           IREP = MULTB (JREP,IJREP+NREP,2)
           IREP = MULTB (IREP,NREP+1,2) ! form gamma*_i
           DO J = 1, NO(JREP)
              JJ = IREPSPI(J,JREP,1)
              DO I = 1, NO(IREP)
                 II = IREPSPI(I,IREP,1)
                 IF (CARITH) THEN
                    AOO(RCW*IJ-1) = PROP(1,II,JJ)
                    AOO(RCW*IJ)   = PROP(2,II,JJ)
                 ELSE
                    AOO(IJ) = PROP(1,II,JJ)*APHASE(1) -
     &                        PROP(2,II,JJ)*APHASE(2)
                 ENDIF
                 IJ = IJ + 1
                 iaoocount = iaoocount + 1
              ENDDO
           ENDDO
        ENDDO
      ENDDO
      IF( (IJ-1).ne.dimpoo) Call Quit('AOO counting error!')
!
!  generate AVV stream
!
      AB = 1
      DO ABREP = 1,NREP
        DO BREP = 1, NREP
           AREP = MULTB (BREP,ABREP+NREP,2)
           AREP = MULTB (AREP,NREP+1,2) ! form gamma*_a
           DO B = 1, NV(BREP)
              BB = IREPSPI(B,BREP,2)
              DO A = 1, NV(AREP)
                 AA = IREPSPI(A,AREP,2)
                 IF (CARITH) THEN
                    AVV(RCW*AB-1) = PROP(1,AA,BB)
                    AVV(RCW*AB)   = PROP(2,AA,BB)
                 ELSE
                    AVV(AB) = PROP(1,AA,BB)*APHASE(1) -
     &                        PROP(2,AA,BB)*APHASE(2)
                 ENDIF
                 AB = AB + 1
                 iavvcount = iavvcount + 1
              ENDDO
           ENDDO
        ENDDO
      ENDDO
      IF( (AB-1).ne.dimpvv) Call Quit('AVV counting error!')

      write(*,*)
      write(*,*) 'Total number of AOO ints:     ',dimpoo
      write(*,*) 'Total number of AOV/AVO ints: ',dimpov
      write(*,*) 'Total number of AVV ints:     ',dimpvv

      SUCC = .true.
!
! ********  flush property integrals to disk
!
      co = NAMEA(1:1)
      if(carith) then         ! indicator for aphase that we have carith
        aphase(1) = 0.0d0
        aphase(2) = 0.0d0
      endif
      write(filename,'(A1,A7)') co,'PROPINT'
      write(*,*) 'Flushing properties to ',filename
      open(unit=PROCPROP,file=filename,form='unformatted',
     &     access='sequential',status='unknown')
      write(PROCPROP) dimpoo,dimpov,dimpvo,dimpvv
      write(PROCPROP) aphase(1),aphase(2)
      write(PROCPROP) (aoo(ixx),ixx=1,RCW*dimpoo)
      write(PROCPROP) (aov(ixx),ixx=1,RCW*dimpov)
      write(PROCPROP) (avo(ixx),ixx=1,RCW*dimpvo)
      write(PROCPROP) (avv(ixx),ixx=1,RCW*dimpvv)
      close(PROCPROP)
!
!********************************************
!***   testing the property streams     *****
!********************************************
!
!     write(*,*)
!     write(*,*) '----- Irep order in the property integrals'
!     write(*,*) '----- <A|I> == <V|O> == beta (bosonic)'
!     write(*,*)
!     write(*,*)
!     DO AIREP = 1, NREP
!        write(*,*) 'AIREP ==  *************',airep,' ***************'
!        write(*,*)
!        ipcount = 0
!        DO IREP = 1, NREP
!           AREP = MULTB (IREP,AIREP+NREP,2)
!           AREP = MULTB (AREP,NREP+1,2)
!           if(airep.ne.multb(arep,irep,2)) 
!    &         stop 'Property irep inconsistency!'
!           write(*,*) 'arep,irep',arep,irep
!           DO I = 1, NO(IREP)
!           DO A = 1, NV(AREP)
!             write(*,*) '       A,I:',A,I
!             ipcount = ipcount + 1
!           ENDDO
!           ENDDO
!        ENDDO
!        write (*,*) '<A|I> integrals in this irep:',ipcount
!        write (*,*)
!     ENDDO
!     write(*,*) '--------------------------------------'
!     write (*,*)
!
!********************************************
!***   printing the property integrals  *****
!********************************************
!
      if(printlevel.lt.6) return

      write(*,*) '---------------- AOO -----------------'
      IJ = 1
      DO IJREP = 1,NREP
        DO JREP = 1, NREP
           IREP = MULTB (JREP,IJREP+NREP,2)
           IREP = MULTB (IREP,NREP+1,2) ! form gamma*_i
           DO J = 1, NO(JREP)
              DO I = 1, NO(IREP)
                 irow = io(irep) + i
                 icol = io(jrep) + j
                 IF (CARITH) THEN
                    write(6,'(A,I6,X,2I4,X,2F16.10)') 'int #',ij,
     &                   irow,icol,AOO(2*IJ-1),AOO(2*IJ)
                 ELSE
                    write(6,'(A,I6,X,2I4,X,F16.10)') 'int #',ij,
     &                   irow,icol,AOO(IJ)
                 ENDIF
                 IJ = IJ + 1
              ENDDO
           ENDDO
        ENDDO
      ENDDO
!
!  AOV type
!
!     write(*,*) '---------------- AOV -----------------'
!     IA = 1
!     DO IAREP = 1, NREP
!        DO AREP = 1, NREP
!           IREP = MULTB (AREP,IAREP+NREP,2)
!           IREP = MULTB (IREP,NREP+1,2)
!           write(*,*) '   ----> IA,I,Arep:',iarep,irep,arep
!           DO A = 1, NV(AREP)
!             DO I = 1, NO(IREP)
!                irow = io(irep) + i
!                icol = io(nrep) + iv(arep) + a
!                 IF (CARITH) THEN
!                   write(6,'(A,I6,X,2I4,X,2F16.10)') 'j/a: ',IA,
!    &                    irow,icol,AOV(2*IA-1),AOV(2*IA)
!                 ELSE
!                   write(6,'(A,I6,X,2I4,X,F16.10)') 'j/a: ',IA,
!    &                    irow,icol,AOV(IA)
!                 ENDIF
!                 IA = IA + 1
!              ENDDO
!           ENDDO
!        ENDDO
!     ENDDO
!
!  AVO type
!
      write(*,*)
      write(*,*) '--------- AVO/AOV (same symmetry order!) --------'
      write(*,*) '--------- this makes eigenvector contraction'
      write(*,*) '--------- much easier!'
      write(*,*)
      AI = 1
      DO AIREP = 1, NREP
         DO IREP = 1, NREP
            AREP = MULTB (IREP,AIREP+NREP,2)
            AREP = MULTB (AREP,NREP+1,2)
!           write(*,*) '   ----> AI,A,Irep:',airep,arep,irep
            DO I = 1, NO(IREP)
               DO A = 1, NV(AREP)
                 irow = io(nrep) + iv(arep) + a
                 icol = io(irep) + i
                  IF (CARITH) THEN
                    write(6,'(A,I6,X2I4,X,2F16.10,A,2F16.10)') 'int # ',
     &              ai,irow,icol,AVO(2*AI-1),AVO(2*AI),'  ||  ',
     &              AOV(2*AI-1),AOV(2*AI)
                  ELSE
                    write(6,'(A,I6,X,2I4,X,F16.10,A,F16.10)') 'int # ',
     &              ai,irow,icol,AVO(AI),'  ||  ',AOV(AI)
                  ENDIF
                  AI = AI + 1
               ENDDO
            ENDDO
         ENDDO
      ENDDO

!
!  AVV type
!
      write(*,*)
      write(*,*) '------------- AVV ----------------'
      write(*,*)
      AB = 1
      DO ABREP = 1,NREP
        DO BREP = 1, NREP
           AREP = MULTB (BREP,ABREP+NREP,2)
           AREP = MULTB (AREP,NREP+1,2) ! form gamma*_a
           DO B = 1, NV(BREP)
              BB = IREPSPI(B,BREP,2)
              DO A = 1, NV(AREP)
                 AA = IREPSPI(A,AREP,2)
                 irow = io(nrep) + iv(arep) + a
                 icol = io(nrep) + iv(brep) + b
                 IF (CARITH) THEN
                    write(6,'(A,I6,X2I4,X,2F16.10)') 'int # ',
     &              ab,irow,icol,AVV(2*AB-1),AVV(2*AB)
                 ELSE
                    write(6,'(A,I6,X,2I4,X,F16.10)') 'int # ',
     &              ab,irow,icol,AVV(AB)
                 ENDIF
                 AB = AB + 1
              ENDDO
           ENDDO
        ENDDO
      ENDDO

      RETURN

 1000 FORMAT (/' Read integrals of type ',A8,
     & A8,' storage info : ',A8)
 1010 FORMAT (/' Largest imaginary part of matrix element',2I5,F10.2)
 1020 FORMAT (/' Property ',A8,' scaled with phase factor ',2F10.4)
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE IDXCOVERAGE(IREPSPI,NSP,ok)

      IMPLICIT none
!     IMPLICIT INTEGER (A-Z)
!
! perform array coverage test for the property integral arrays
! where absolute Kramers-unique spinor numbers are related to relative
! symmetry-grouped arrays.
!
!----------------Common Blocks-----------------------------------------
!
#include "../relccsd/symm.inc"    !provides MXREP and MULTB
!
!----------------Calling variables---------------------------------------------
!
      INTEGER                       :: NSP,ok
      INTEGER                       :: IREPSPI(NSP,MXREP,2)
!
!----------------Local variables---------------------------------------------
!
      integer,allocatable           :: acv(:,:)
      integer                       :: irep,jrep,arep,brep
      integer                       :: airep,ijrep,abrep
      integer                       :: i,ii,j,jj,a,aa,b,bb
      integer                       :: rem
!
!----------------Executable code ---------------------------------------------
!
      write(*,*) '*********'
      write(*,*) '*** Array coverage test triggered.'
      write(*,*) '*********'

      ok = 0

      allocate(acv(nsp,nsp))
      acv = 1

      write(*,*) '*** testing VO/OV...'
      DO AIREP = 1, NREP
        DO IREP = 1, NREP
          AREP = MULTB (IREP,AIREP+NREP,2)
          AREP = MULTB (AREP,NREP+1,2)
!         write(*,*) 'AIREP,AREP,IREP',airep,arep,irep
            DO I = 1, NO(IREP)
               II = IREPSPI(I,IREP,1)
               DO A = 1, NV(AREP)
                  AA = IREPSPI(A,AREP,2)
                  if(AA.gt.nsp) then
                    stop('Internal error IDXCoverage')
                  endif
                  if(II.gt.nsp) then
                    stop('Internal error IDXCoverage')
                  endif
                  acv(AA,II)=0
                  acv(II,AA)=0
               ENDDO
            ENDDO
        ENDDO
      ENDDO

      write(*,*) '*** testing OO...'
      DO IJREP = 1,NREP
        DO JREP = 1, NREP
           IREP = MULTB (JREP,IJREP+NREP,2)
           IREP = MULTB (IREP,NREP+1,2) ! form gamma*_i
           DO J = 1, NO(JREP)
              JJ = IREPSPI(J,JREP,1)
              DO I = 1, NO(IREP)
                 II = IREPSPI(I,IREP,1)
                 acv(II,JJ)=0
              ENDDO
           ENDDO
        ENDDO
      ENDDO

      write(*,*) '*** testing VV...'
      DO ABREP = 1,NREP
        DO BREP = 1, NREP
           AREP = MULTB (BREP,ABREP+NREP,2)
           AREP = MULTB (AREP,NREP+1,2) ! form gamma*_a
           DO B = 1, NV(BREP)
              BB = IREPSPI(B,BREP,2)
              DO A = 1, NV(AREP)
                 AA = IREPSPI(A,AREP,2)
                 acv(AA,BB)=0
              ENDDO
           ENDDO
        ENDDO
      ENDDO

      rem=sum(acv)
      if(rem.ne.0) then
         write(*,*) ' --->>> ERROR <<<---'
         write(*,*) ' --->>> index coverage: ',rem,
     &       ' should be zero !'
         ok = -1
      endif
      deallocate(acv)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE print_input_info(iw,nrep)
      use polprp_cfg     !these data are visible on master node only!

      implicit none

      integer               :: iw         ! I/O file handle for proper screen output
      integer               :: nrep       ! number of total symmetries in the calculation
!
      integer               :: i,n,icount

      write(iw,*) 'POLPRP control variables:'
      write(iw,*) '-------------------------'
      write(iw,*)
      write(iw,*) 'we do extended calculation: ',polprp_doextended
      write(iw,'(A,F10.5)') 
     &   'threshold for elements in ADC matrix: ',polprp_writethr
      write(iw,*) 'we perform a Davidson diagonalization:',polprp_dodiag
      write(iw,*) 'The way final state symmetries are chosen is:'
!
! determine symmetries to be calculated
!
      n = 0
      do i = 1,nrep
        n = n + polprp_statesym(i)
      enddo
      if(n.eq.0) then
        write(iw,*) 'Program checks symmetries automatically'
        do i=1,nrep
          polprp_statesym(i) = i
        enddo
      else
        write(iw,*) 'User has selected the following final'
        write(iw,*) 'state symmetries:'
        icount=1
        do i=1,nrep
          if(polprp_statesym(i).ne.0) then
            write(iw,*) '   No: ',icount,
     &                  '  symmetry: ',polprp_statesym(i)
            icount = icount + 1
          endif
        enddo
      endif
      write(iw,*)
!
! report diagonalization parameters
!
      write(iw,*) 'Davidson control variables:'
      write(iw,*) '---------------------------'
      write(iw,*)

      if(polprp_dodiag) then
        WRITE(IW,*) 'Number of eigenstates to be calculated (roots)',
     &     polprp_davroots
        WRITE(IW,*) 'Size of Krylov space: ',polprp_davmaxsp
        WRITE(IW,*) 'Number of Davidson macro iterations: ',
     &     polprp_davmaxit
        WRITE(IW,'(A,E20.4)') ' Requested energy conv.: ',polprp_davconv
        IF(polprp_davooc.eqv..false.) THEN
           WRITE(IW,*) 'Trying full incore storage of Davidson vectors.'
           WRITE(IW,*) 'Will default to out-of-core if resources are'
           WRITE(IW,*) 'not available.'
        ENDIF
      endif

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRECONTRACT(iounit,iw)
!
!
!  this routine is called by master only when parallel execution takes place.
!  Then the nodes provide integrals.
!
!

#if defined (VAR_MPI)
      use interface_to_mpi
#endif

      implicit none

      integer               :: iounit         ! I/O file handle for writing to file
      integer               :: iw             ! file handle for writing to screen
!
! local variables
!
      real*8, allocatable, dimension(:)        :: vvvvz
      real*8, allocatable, dimension(:)        :: vvoo1
      real*8, allocatable, dimension(:)        :: vvoop
      integer                                  :: n1,ixx
      integer                                  :: ialloc
      integer                                  :: signal
      integer*8                                :: n8,i8
      COMPLEX*16                               :: A0,A1
      character*9                              :: preconfn='VVVVXVVOO'
      logical                                  :: isthere


#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/ccpar.inc"
!
! execution
!
      if(carith) then
        CALL PST('Performing (complex) precontraction+')
        if(rcw.ne.2) call quit ('intenal error in PRECONTRACT')
      else
        CALL PST('Performing (real) precontraction+')
        if(rcw.ne.1) call quit ('intenal error in PRECONTRACT')
      endif
      A0=(0.0D00,0.0D00)
      A1=(1.0D00,0.0D00)


#if defined (VAR_MPI)
!
!  send the slaves a start signal for int supply only
!
      signal = 2  !start code.
      write(iw,*) 'Master indicates PRECONTRACTION (signal=2)'
      call interface_mpi_bcast(signal,1,MASTER,
     &                             global_communicator)
#endif

! allocate buffers and read in integrals

      n1 = ivvoott(nrep+1)
      allocate(vvoo1(n1*rcw))
      allocate(vvoop(n1*rcw))

      n8 = ivvvvtt(nrep+1)
      allocate(vvvvz(n8*rcw),stat=ialloc)
      if(ialloc.ne.0) then
        STOP 'VVVV buffer could not be allocated.'
      else
        write(*,*) 'Master allocated',n8,' words for VVVV buffer.'
      endif

#if defined (VAR_MPI)
      do i8=1,n8*rcw
        vvvvz(i8) = -4096.0d0
      enddo
      call rdvvvv(vvvvz)  ! master reads its own chunks
      write(*,*) 'masters chunk reading ... complete.'
      call mpi_vvvv_complete(vvvvz,n8*rcw)  !.... and completes the stream.
      write(*,*) 'VVVV check in PRECONTRACT...'
      do i8=1,n8
        if(vvvvz(i8).eq. -4096.0d0)
     &    call quit('Internal error: gap in the VVVV stream!')
      enddo
      write(*,*) '.... success. No leaks.'
#else
      call rdvvvv(vvvvz)  ! master reads complete stream in serial mode
#endif
      call getvvoo(vvoo1)
!
!  do the contraction over all symmetries
!
      CALL CNTRCT('N','N',NVVT,NOOT,NVVT,A1,vvvvz,vvoo1,A0,vvoop,NREP)
!
!  check if file already present. If yes, delete it.
!
      Inquire(file=preconfn,exist=isthere)
      If(isthere) then
         write(*,*) 'Found a precontraction file. Deleting...'
         open(unit=iounit, file=preconfn, status='old')
         close(iounit, status='delete')
      Endif
!
!  open corresponding file and write out contracted integrals.
!
      open(unit=iounit,file=preconfn, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
      write(iounit) (vvoop(ixx), ixx=1,n1*rcw)
      close(iounit)
!
!  free buffers immediately 
!
      deallocate(vvoop)
      deallocate(vvoo1)
      deallocate(vvvvz)

      write(*,*) 'PRECONTRACTION finished. Buffers freed'

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine compare_evc_r(io1,desrep,nroots,ladc,
     &                         davibase,fullbase)
      implicit none
!
! this routine evaluates eigenvectors for the available Davidson roots
! with respect to the full diagonalization results. It determines phases and
! absolute deviations
!
      integer        :: io1,desrep,nroots,ladc
      character*8    :: davibase
      character*7    :: fullbase 
!
! local variables
! 
      character*11                        :: daviname
      character*10                        :: fullname
      logical                             :: fullisthere
      integer                             :: io2
      integer                             :: i,j,k
      integer                             :: k1,k2
      integer                             :: nmain
      real*8                              :: lambda1,lambda2,dummy
      real*8,parameter                    :: autoev = 27.2113957D0
      real*8,allocatable,dimension (:)    :: e1,e2
      real*8                              :: devmod,devabs,devphs
!
! executable code
! 
      Call PST('Comparing D/F eigenvectors+')

      IF(desrep.GT.9) THEN
        WRITE(daviname,'(A8,A1,I2)') davibase,'.',desrep
        WRITE(fullname,'(A7,A1,I2)') fullbase,'.',desrep
      ELSE
        WRITE(daviname,'(A8,A2,I1)') davibase,'.0',desrep
        WRITE(fullname,'(A7,A2,I1)') fullbase,'.0',desrep
      ENDIF
!
! do we have an eigenvector file from fulldia ?
!
      fullisthere=.false.
      INQUIRE(FILE=fullname,EXIST=fullisthere)
      if(.not.fullisthere) then
        write(*,*) 'The eigenvector file from a full diagonalization'
        write(*,*) 'is not available. No comparison is made. Returning.'
        return
      endif
      write(*,*) 'Comparing eigenvectors for',nroots,' roots'
      write(*,*)
!
! comparing eigenvectors
!
      io2 = io1 + 1
      open(io1,file=daviname,access='sequential',
     &     status='unknown',form='unformatted')
      open(io2,file=fullname,access='sequential',
     &     status='unknown',form='unformatted')

      allocate(e1(ladc))
      allocate(e2(ladc))
      e1 = 0.0d0
      e2 = 0.0d0
!
! three integer reads from the Davidson file
!
      read(io1) k
      read(io1) nmain
      read(io1) k
      write(*,*) 'Main space has dimension',nmain
      write(*,*)

      write(*,'(18X,2A7,4X,3A8)') 'root1','root2',' devmod',
     &                          'devabs','devphs'
      do k=1,nroots
        read(io1) k1
        read(io2) k2
        if(.not.(k1.eq.k2 .and. k1.eq.k)) then
          write(*,*) 'Inconsistency in root numbering!'
          write(*,*) 'This is serious, continuation questionable!'
          write(*,*) 'k,k1,k2:',k,k1,k2
          exit  !leave this do loop
        endif
        read(io1) lambda1
        read(io2) lambda2
        read(io1) dummy
        read(io2) dummy
        read(io1) (e1(i),i=1,ladc)
        read(io2) (e2(i),i=1,ladc)
!
! compare vectors
!
        devmod = 0.0d0
        devabs = 0.0d0
        devphs = 0.0d0
        do i=1,ladc
          devmod = devmod + (abs(e1(i)) - abs(e2(i)))
          devabs = devabs + (e1(i) - e2(i))
          devphs = devphs + (e1(i) + e2(i))
        enddo
        write(*,'(A,I3,2X,2F12.6,2X,3F14.8)') 'Comparing root',k,
     &         lambda1,lambda2,devmod,devabs,devphs
        write(*,*)
        do i=1,min(10,nmain)
          write(*,*) '       ',e1(i),e2(i)
        enddo
        write(*,*) '..................'
        write(*,*)
      enddo
!
! release ressources
!
      deallocate(e1)
      deallocate(e2)

      close(io1)
      close(io2)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
