!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE property_handler(ioch,ioch2,lavail,nsp,irepspi)
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
      logical      :: existing,success,printints
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
      else
        write(*,*) '  ** Consistency check passed with:',nsp*nsp
        write(*,*) '  ** Size of OO array:',npaoo
        write(*,*) '  ** Size of OV array:',npaov
        write(*,*) '  ** Size of VO array:',npavo
        write(*,*) '  ** Size of VV array:',npavv
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
      naoplist(1) = 'X_DIPLEN'
      naoplist(2) = 'Y_DIPLEN'
      naoplist(3) = 'Z_DIPLEN'
!
! loop over the three properties and store the integrals in the
! corresponding files.
!

      printints = .false.
      Do iprop = 1,3

        nameofprop=naoplist(iprop)

        success = .false.

        CALL read_flush_prop(success,nameofprop,NSP,IREPSPI,
     &                       npaoo,npaov,npavo,npavv,
     &                       CBUF,CPOO,CPOV,CPVO,CPVV,
     &                       printints)

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
      SUBROUTINE BUILD_TRAMO(desrep,ioch,eps)
 
      IMPLICIT INTEGER (A-Z)
!     IMPLICIT none
!
!---------------Description--------------------------------------------
!
!     Calculates transition moments from the transformed
!     x,y,z dipole moment integrals and the Davidson eigenvectors
!
!     Important: Also for real two-electron integrals (high symmetry and/or 
!     NONREL,SPINFREE Hamiltonians) the property integrals are in complex form.
!     These cases will be handled by the corresponding phase factors stored at the
!     beginning of the property integral file.
!
!---------------Calling variables--------------------------------------
!
      INTEGER                          :: desrep,ioch
      REAL*8,dimension(:)              :: eps
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/complex.inc"
#include "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!

! filename of property integrals

      Character*8,dimension(3)               :: trafile

! constants

      Complex*16                             :: A0,A1,AI

! index arrays for the nontotally symmetric (NTS) property arrays

      INTEGER,allocatable,dimension(:)       :: JPOOXX,JPOVXX
      INTEGER,allocatable,dimension(:)       :: JPVOXX,JPVVXX
      INTEGER,allocatable,dimension(:)       :: JPVOOV

! offset arrays for NTS two-electron contributions
                                            
      INTEGER,allocatable,dimension(:)       :: ABIJTT
      INTEGER,allocatable,dimension(:,:)     :: JAVVOO
      INTEGER,allocatable,dimension(:,:)     :: JBVVOO
      INTEGER,allocatable,dimension(:)       :: MVVOOT
      INTEGER,allocatable,dimension(:)       :: BVVOOT

! phase information for the property x,y,z arrays in case of real arithmetic

      Complex*16,allocatable,dimension(:)    :: pphase
      Complex*16                             :: pphasex
      Complex*16                             :: pphasey
      Complex*16                             :: pphasez

! intermediate and final F_ai, F_aibj matrices (as vectors)

      REAL*8,allocatable,dimension(:)        :: faipx,faipy,faipz
      Complex*16,allocatable,dimension(:)    :: faitx,faity,faitz
      Complex*16,allocatable,dimension(:)    :: fabijtx,fabijty,fabijtz

! dipole moment integrals in MO basis for all three components

      REAL*8,allocatable,dimension(:)        :: doox,dooy,dooz
      REAL*8,allocatable,dimension(:)        :: dovx,dovy,dovz
      REAL*8,allocatable,dimension(:)        :: dvox,dvoy,dvoz
      REAL*8,allocatable,dimension(:)        :: dvvx,dvvy,dvvz
                                            
! auxiliary arrays for contractions

      REAL*8,allocatable,dimension(:)        :: aac,aki,aak,abj
                                            
! integral storage arrays

      REAL*8,allocatable,dimension(:)        :: vvoo1
      REAL*8,allocatable,dimension(:)        :: vvoo2
      REAL*8,allocatable,dimension(:)        :: vvoo3
      REAL*8,allocatable,dimension(:)        :: vvoo4
      REAL*8,allocatable,dimension(:)        :: vvoop
                                            
      REAL*8,allocatable,dimension(:)        :: vooo1,vooo2
      REAL*8,allocatable,dimension(:)        :: vovo1,vovo2
      REAL*8,allocatable,dimension(:)        :: vovv1,vovv2
                                            
      REAL*8,allocatable,dimension(:)        :: oooo

      REAL*8,allocatable,dimension(:)        :: vabikx
      REAL*8,allocatable,dimension(:)        :: vabiky
      REAL*8,allocatable,dimension(:)        :: vabikz

      Integer                                :: n1,n2,n3,n4
      Integer                                :: ialloc,iat

! some parameters

      Logical                                :: LFA=.false.
      Logical                                :: LTR=.true.

      integer*8                              :: n8
      character*9                            :: preconfn='VVVVXVVOO'
      character*30                           :: iorder
      integer                                :: isatctl
      logical                                :: okflag
!
!---------------interface region --------------------------------------
!
      interface

        SUBROUTINE denomvots(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE propfetch(i1,ca1,
     &                          i2,i3,i4,i5,
     &                          ra1,ra2,ra3,
     &                          ra4,ra5,ra6,
     &                          ra7,ra8,ra9,
     &                          ra10,ra11,ra12)
          integer                      :: i1
          complex*16, dimension(:)     :: ca1
          integer                      :: i2,i3,i4,i5
          real*8, dimension(:)         :: ra1,ra2,ra3
          real*8, dimension(:)         :: ra4,ra5,ra6
          real*8, dimension(:)         :: ra7,ra8,ra9
          real*8, dimension(:)         :: ra10,ra11,ra12
        END SUBROUTINE

        INTEGER FUNCTION XCOLDET(i1)
          integer                      :: i1
        END FUNCTION

        SUBROUTINE rcaxpy(i1,l1,c1,ra1,ca1)
          integer                      :: i1
          logical                      :: l1
          complex*16                   :: c1
          real*8, dimension(:)         :: ra1
          complex*16, dimension(:)     :: ca1
        END SUBROUTINE

        SUBROUTINE rca_negate(i1,ra1)
          integer                      :: i1
          real*8, dimension(:)         :: ra1
        END SUBROUTINE

        SUBROUTINE TM_GEN(i1,i2,i3,ch1,
     &                    l1,
     &                    ca1,ca2,ca3,
     &                    ca4,ca5,ca6)
          integer                       :: i1,i2,i3
          character*30                  :: ch1
          logical                       :: l1
          complex*16, dimension(:)      :: ca1,ca2,ca3
          complex*16, dimension(:)      :: ca4,ca5,ca6
        END SUBROUTINE

      end interface
!
!---------------Executable code--------------------------------------
!
      call pst('Calculating transition moments+')
      write(*,*)
      write(*,*) '   Remember: Even in SPINFREE/NONREL cases the TMs'
      write(*,*) '   can be complex due to the quaternionic structure'
      write(*,*) '   of the wave function.'
      write(*,*)
      write(*,*) '   **********************'
      write(*,*) '   Symmetry: ',desrep
      write(*,*) '   **********************'
      write(*,*)
!
! initialization of constant variables
!
      trafile = (/'XPROPINT','YPROPINT','ZPROPINT'/)
      do i=1,30
        iorder(i:i)=' '
      enddo
      A0=(0.0d0,0.0d0)
      A1=(1.0d0,0.0d0)
      AI=(0.0d0,1.0d0)
!
! Generate index arrays for OO,OV,VO,VV property arrays
!
      n1 = nrep + 1
      allocate(JPOOXX(n1))
      allocate(JPVOXX(n1))
      allocate(JPOVXX(n1))
      allocate(JPVVXX(n1))
      allocate(JPVOOV(n1))

      JPOOXX=0
      JPVOXX=0
      JPOVXX=0
      JPVVXX=0
      JPVOOV=0

      DO IREP = 1, NREP
        JPOOXX(IREP+1) = JPOOXX(IREP) + MOO(IREP)
        JPOVXX(IREP+1) = JPOVXX(IREP) + MOV(IREP)
        JPVOXX(IREP+1) = JPVOXX(IREP) + MVO(IREP)
        JPVVXX(IREP+1) = JPVVXX(IREP) + MVV(IREP)

        JPVOOV(IREP+1) = JPVOOV(IREP) + MVO(IREP)*MOV(IREP)
      ENDDO
!
! Generate index arrays for NTS two-particle contribs
! and check for correct coupling block lengths
!
      allocate(ABIJTT(NREP))

      allocate(JAVVOO(NREP,NREP))
      allocate(JBVVOO(NREP,NREP))
      allocate(MVVOOT(NREP))
      allocate(BVVOOT(NREP))

      ABIJTT  = 0
      JAVVOO  = 0
      JBVVOO  = 0
      MVVOOT  = 0
      BVVOOT  = 0

      DO JREP = 1, NREP
        DO IREP = 1, NREP
           IJREP = MULTB(IREP,JREP,1)
           IJBREP = MULTB(IREP+NREP,JREP+NREP,2)  !ABIJTT refers to bosonic ireps!

           ABIJTT(IJBREP) = ABIJTT(IJBREP) + NVVT(IREP)*NOOT(JREP)
           JAVVOO(IREP,JREP) = MVVOOT(IJREP)
           JBVVOO(IREP,JREP) = BVVOOT(IJREP)
           MVVOOT(IJREP) = MVVOOT(IJREP) + NVVOT(IREP)*NO(JREP)
           BVVOOT(IJREP) = BVVOOT(IJREP) + NV(IREP)*NVOOT(JREP)

        ENDDO
      ENDDO

      DO irep = 1, NREP
        if(ABIJTT(irep).ne.xcoldet(irep)) then
          write(*,*) 'Internal error in coupling block lengths!'
          write(*,*) 'IREP:',irep,' ABIJTT(irep):',ABIJTT(irep),
     &               ' XCOLDET(irep):',xcoldet(irep)
          Call quit('**********')
        endif
      ENDDO

!*****************************************
! complex array for the phases.
!*****************************************

      allocate(pphase(3))
      pphase = A0
!
!***************************************************
! ****  Transition element vectors. Attention
! ****  Even in a real calculation they will be kept
! ****  complex!
! ****  the final transition moment arrays are always complex.
! ****  the intermediate transition moment arrays follow CARITH.
!***************************************************
!
      n1 = mvo(desrep)

      allocate(faitx(n1))
      allocate(faity(n1))
      allocate(faitz(n1))
      faitx = A0
      faity = A0
      faitz = A0
      
      allocate(faipx(rcw*n1))
      allocate(faipy(rcw*n1))
      allocate(faipz(rcw*n1))
      faipx = 0.0d0
      faipy = 0.0d0
      faipz = 0.0d0
!
!***************************************************
! ****  Arrays holding the property integrals. 
! ****  They are either
! ****  complex if CARITH=T or real, if CARITH=F.
! ****  However, in the latter case a phase information
! ****  indicating purely imaginary property integrals
! ****  has to be additionally stored and processed.
!***************************************************
!
      idaoo = JPOOXX(NREP+1)   !  dimensions
      idaov = JPOVXX(NREP+1)
      idavo = JPVOXX(NREP+1)
      idavv = JPVVXX(NREP+1)
      n1 = idaoo*rcw   ! lengths
      n2 = idaov*rcw
      n3 = idavo*rcw
      n4 = idavv*rcw
      allocate(doox(n1),dooy(n1),dooz(n1))
      allocate(dovx(n2),dovy(n2),dovz(n2))
      allocate(dvox(n3),dvoy(n3),dvoz(n3))
      allocate(dvvx(n4),dvvy(n4),dvvz(n4))
      doox = 0.0d0
      dooy = 0.0d0
      dooz = 0.0d0
      dovx = 0.0d0
      dovy = 0.0d0
      dovz = 0.0d0
      dvox = 0.0d0
      dvoy = 0.0d0
      dvoz = 0.0d0
      dvvx = 0.0d0
      dvvy = 0.0d0
      dvvz = 0.0d0

!*********************************************************************
!  VVOO-type integral arrays, they always remain allocated since they
!  hold integrals repeatedly needed in nearly every step.
!*********************************************************************

      n1 = ivvoott(nrep+1)
      allocate(vvoo1(n1*rcw))
      allocate(vvoo2(n1*rcw))
      vvoo1 = 0.0d0
      vvoo2 = 0.0d0

      n1 = max0(lvvoot(nrep+1),jvovo(nrep+1),
     &          ivvoott(nrep+1),ivvoot(nrep+1),
     &          JPVOOV(nrep+1))

      allocate(vvoo3(n1*rcw))
      allocate(vvoo4(n1*rcw))
      vvoo3 = 0.0d0
      vvoo4 = 0.0d0
!
!  auxiliary arrays
!
      n1 = MVV(1)
      allocate(aac(n1*rcw))
      n1 = MOO(1)
      allocate(aki(n1*rcw))
      n1 = MVO(1)
      allocate(aak(n1*rcw))
      n1 = JPVOXX(NREP+1) 
      allocate(abj(n1*rcw))
      aac = 0.0d0
      aki = 0.0d0
      aak = 0.0d0
      abj = 0.0d0
!______________________________________________________
!|
!|
!|
!|
!|
!|       START CALCULATIONS
!|
!|
!|
!|
!|_____________________________________________________
!
!  Initialize VVOO integral arrays
!  v_abij ---> vvoo1  /  v*_cbjk ---> vvoo2
!
      n1 = ivvoott(nrep+1)
      call getvvoo(vvoo1)
      call xcopy(n1,vvoo1,1,vvoo2,1)
      if(carith) call conjuga(n1,vvoo2,1)
      call denomvvoo(eps,vvoo1,vvoo1)
      call denomvvoo(eps,vvoo2,vvoo2)
!
!  get the property integrals and phases from corresponding files.
!
      CALL propfetch(ioch,pphase,    !fetch properties, we transfer dimensions!
     &               idaoo,idaov,idavo,idavv,
     &               doox,dooy,dooz,
     &               dovx,dovy,dovz,
     &               dvox,dvoy,dvoz,
     &               dvvx,dvvy,dvvz)

      pphasex = pphase(1)
      pphasey = pphase(2)
      pphasez = pphase(3)
!
!------------------------
!------  F_ai(0)
!------------------------
!
      n1 = MVO(desrep)
      off = JPVOXX(desrep) + 1
      call xcopy(n1,dvox(off),1,faipx,1)
      call xcopy(n1,dvoy(off),1,faipy,1)
      call xcopy(n1,dvoz(off),1,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th order                     '
      isatctl = 0
      CALL TM_GEN(ioch,desrep,isatctl,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
!
!------------------------
!------  F_ai(1)
!------------------------
!
!  contract v_ab,ji and D_jb to F_ai
!  v_a>b,j>i  ---> -v_a>b,i>j
!  -v_a>b,i>j  ---> v_ab,i>j (1LS1)
!  v_ab,i>j  ---> C_ai,jb (SRT20M),(vvoo3)
!
      vvoo3 = 0.0d0
      vvoo4 = 0.0d0
      n1 = ivvoott(nrep+1)
      CALL RCA_NEGATE(n1,vvoo1)
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoo1,vvoo4)
      CALL RCA_NEGATE(n1,vvoo1)      !reverse negation for the vvoo1 array
      CALL SRT20M(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &            JPVOOV,JJVO,JJOV,vvoo4,vvoo3)


!  contract C_ai,jb (vvoo3)  and Djb over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JPVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)
!
! add result on top of the final arrays
!
      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)
!
!  calculate corresponding 0th + 1st order TMs
!
      iorder='0th + 1st order               '
      CALL TM_GEN(ioch,desrep,isatctl,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
!
!
!
!  enter calculation for second-order contributions
!
!------------------------
!------  F_ai(2)A
!------------------------
!
!  sort  v_ab,jk  -->  v_a,bjk    (vvoo3)
!  sort v*_cb,jk  -->  v*_c,bjk   (vvoo4)
!         v_a,bjk  *  v*_c,bjk   -->  aux1
!
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo2,vvoo4)
      call cntrct ('N','T',nv,nv,nvoot,a1,vvoo3,vvoo4,A0,aac,nrep)
!
!  contract A_ac and D_ci(x,y,z) over c in special loop because
!  D_ci is not totally symmetric and the order of Y_c in D_ci does not
!  follow the one in A_ac.
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         crep = multb(desrep+nrep,irep,2)
         arep = crep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NV(crep)
         ioff1 = JJVV(arep,crep) + 1
         ioff2 = JPVOXX(desrep) + JJVO(crep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,aac(ioff1),m,dvox(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aac(ioff1),m,dvoy(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aac(ioff1),m,dvoz(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo
      n1 = MVO(desrep)


      CALL XSCAL(n1,-A1/2.0d0,faipx,1)
      CALL XSCAL(n1,-A1/2.0d0,faipy,1)
      CALL XSCAL(n1,-A1/2.0d0,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2)A finished'
!
!------------------------
!------  F_ai(2)B
!------------------------
!
!  sort  v*_bc,jk (vvoo2)  -->  v*_bcj,k  (vvoo3)  (1T3)
!  sort  v_bc,ji  (vvoo1)  -->  v_bcj,i   (vvoo4)  (1T3)
!
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             vvoo2,vvoo3)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             vvoo1,vvoo4)
      CALL CNTRCT('T','N',NO,NO,NVVOT,A1,vvoo3,vvoo4,A0,aki,NREP)
!
!  contract D_ak(x,y,z) and A_ki over k in special loop because
!  D_ak is not totally symmetric!
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         arep = multb(desrep+nrep,irep,2)
         krep = irep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NO(krep)
         ioff1 = JPVOXX(desrep) + JJVO(arep,krep) + 1
         ioff2 = JJOO(krep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,dvox(ioff1),m,aki(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvoy(ioff1),m,aki(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvoz(ioff1),m,aki(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      CALL XSCAL(n1,-A1/2.0d0,faipx,1)
      CALL XSCAL(n1,-A1/2.0d0,faipy,1)
      CALL XSCAL(n1,-A1/2.0d0,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2)B finished'
!
!------------------------
!------  F_ai(2)C
!------------------------
!
!
!  sort  v_ab,ij  (vvoo1)  -->  v_ai,bj    (vvoo3)  (1TT4)
!  v*_cb,kj = v*_bc,jk !
!  sort  v*_bc,jk (vvoo2)  -->  v*_bj,ck   (vvoo4)  (1TT4)

!  additionally, allocate proper array
!
      vvoo3 = 0.0d0
      vvoo4 = 0.0d0
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,vvoo2,vvoo4)
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      vovo1 = 0.0d0
      if(ialloc.ne.0)  call quit('vovo1 allocation problem (tramo)')
!
!  do contraction of v_ai,bj  *  v*_bj,ck  =  A_ai,ck
!
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vvoo3(OFF1),M,vvoo4(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
!  do contraction of A_ai,ck  *  D_ck (determine proper symmetry!)
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      M = MVO(desrep)
      K = MVO(desrep)
      OFF1 = jvovo(desrep) + 1
      OFF2 = JPVOXX(desrep) + 1
      CALL XGEMV ('N',M,K,A1,vovo1(OFF1),M,dvox(OFF2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vovo1(OFF1),M,dvoy(OFF2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vovo1(OFF1),M,dvoz(OFF2),1,A0,faipz,1)
      deallocate(vovo1)
!
!  scale result by +0.5d0 and add to final values
!
      n1 = MVO(desrep)
      CALL XSCAL(n1,A1/2.0d0,faipx,1)
      CALL XSCAL(n1,A1/2.0d0,faipy,1)
      CALL XSCAL(n1,A1/2.0d0,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)
      
      write(*,*) '------  F_ai(2)C finished'
!
!------------------------
!------  F_ai(2,1)
!------------------------
!
!  only temporarily allocate space for the VOVV integrals.
!
      iat = 0
      n1 = ivovvt(nrep+1)
      allocate(vovv1(n1*rcw),stat=ialloc); iat = iat + ialloc
      n1 = lvovvt(nrep+1)
      allocate(vovv2(n1*rcw),stat=ialloc); iat = iat + ialloc
      if(ialloc.ne.0) STOP 'VOVV allocation problem in TRAMO!'
!
!  sort V_aj,bc  --->  Va,jbc   (1S2)    vovv1 --> vovv2
!  sort v_bc,kj  --->  v_jbc,k  (SRT22)  vvoo1 --> vvoo3
!  then contract Va,jbc  *  v_jbc,k  ---> A_ak (tot symm)
!  divide A_ak by (eps_a  -  eps_k), and multiply by -1.0
!  and contract A_ak with D_ki(x,y,z) (special loop!)
!
      CALL RDVOVV(vovv1); vovv2 = 0.0d0
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NV,NO,NVVT,LVOVVT,LLOVVT,
     &             vovv1,vovv2)
      CALL SRT22 (NREP,MULTB,LFA,NVVT,NO,NO,NOVVT,LOOVVT,LLOVVT,
     &            vvoo1,vvoo3)
      CALL CNTRCT ('N','N',NV,NO,NOVVT,A1,vovv2,vvoo3,A0,aak,NREP)
      CALL DENOMVOTS (EPS,aak)
      CALL XSCAL(MVO(1),-A1,aak,1)
      deallocate(vovv1,vovv2)
!
!  contract A_ak/(e_a - e_k) and D_ki(x,y,z) over k in special loop because
!  D_ki is not totally symmetric and the order of Y_k in D_ki does not
!  follow the one in A_ak.
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         krep = multb(desrep+nrep,irep,2)
         arep = krep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NO(krep)
         ioff1 = JJVO(arep,krep) + 1
         ioff2 = JPVOXX(desrep) + JJOO(krep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,doox(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooy(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooz(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,1) finished'
      
!
!------------------------
!------  F_ai(2,2)
!------------------------
!
      iat = 0
      n1 = ivooot(nrep+1)
      allocate(vooo1(n1*rcw),stat=ialloc); iat = iat + ialloc
      n1 = lovoot(nrep+1)
      allocate(vooo2(n1*rcw),stat=ialloc); iat = iat + ialloc
      if(ialloc.ne.0) STOP 'VOOO allocation problem in TRAMO!'
!
!  sort v_ab,jk  --->  v_a,bjk   (still in vvoo1 ---> vvoo3)
!  sort V*_bl,jk --->  V*_bjk,l  vooo1 ---> vooo2
!  contract  v_a,bjk  *  V*_bjk,l  ---> A_al (aak)
!  divide A_ak by (eps_a  -  eps_k), and multiply by -1.0
!  and contract A_ak with D_ki(x,y,z)
!
      call SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      CALL GETVOOO(vooo1);vooo2 = 0.0d0
      IF(CARITH) CALL CONJUGA(ivooot(nrep+1),vooo1,1)
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           vooo1,vooo2)
      CALL CNTRCT ('N','N',NV,NO,NVOOT,A1,vvoo3,vooo2,A0,aak,NREP)
      CALL DENOMVOTS (EPS,aak)
      CALL XSCAL(MVO(1),-A1,aak,1)
      deallocate(vooo1,vooo2)
!
!  special loop as above for contracting A_al and D_li (identical to
!  A_ak and D_ki)
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      do irep = 1, nrep
         krep = multb(desrep+nrep,irep,2)
         arep = krep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NO(krep)
         ioff1 = JJVO(arep,krep) + 1
         ioff2 = JPVOXX(desrep) + JJOO(krep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,doox(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooy(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooz(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)
      
      write(*,*) '------  F_ai(2,2) finished'
!
!-------------------------------------
!------  F_ai(2,3), a DVV contribution
!-------------------------------------
      n1 = ivooot(nrep+1)
      allocate(vooo1(n1*rcw),stat=ialloc)
      n1 = lovoot(nrep+1)
      allocate(vooo2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOOO allocation problem in TRAMO!'
!
!  sort v_bc,jk  --->  v_b,cjk   (still in vvoo1 ---> vvoo3)
!  sort V*_ci,jk --->  V*_cjk,i  vooo1 ---> vooo2
!  contract  v_b,cjk  *  V*_cjk,i  ---> A_bi (aak)
!  divide A_bi by (eps_b  -  eps_i)
!  and contract D_ab(x,y,z) with A_bi
!
      call SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      CALL GETVOOO(vooo1)
      IF(CARITH) CALL CONJUGA(ivooot(nrep+1),vooo1,1)
      vooo2 = 0.0d0
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           vooo1,vooo2)
      CALL CNTRCT ('N','N',NV,NO,NVOOT,A1,vvoo3,vooo2,A0,aak,NREP)
      deallocate(vooo1,vooo2)
      CALL DENOMVOTS (EPS,aak) !A_bi now in aak.
!
! do special contraction of D_ab and A_bi over b
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         arep = multb(desrep+nrep,irep,2)
         brep = irep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NV(brep)
         ioff1 = JPVVXX(desrep) + JJVV(arep,brep) + 1
         ioff2 = JJVO(brep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,dvvx(ioff1),m,aak(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvy(ioff1),m,aak(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvz(ioff1),m,aak(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,3) finished'
!
!-------------------------------------
!------  F_ai(2,4), a DVV contribution
!-------------------------------------
!
!  only temporarily allocate space for the VOVV integrals.
!
      n1 = ivovvt(nrep+1)
      allocate(vovv1(n1*rcw),stat=ialloc)
      n1 = lvovvt(nrep+1)
      allocate(vovv2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOVV allocation problem in TRAMO!'
!
!  sort v_bc,ij  --->  v_jbc,i  (SRT22)  vvoo1 --> vvoo3
!  sort V_dj,bc  --->  Vd,jbc   (1S2)    vovv1 --> vovv2
!  then contract Vd,jbc  *  v_jbc,i  ---> A_di (tot symm)
!  divide A_di by (eps_d  -  eps_i)
!  and contract D_ad(x,y,z) with Adi
!
      CALL SRT22 (NREP,MULTB,LFA,NVVT,NO,NO,NOVVT,LOOVVT,LLOVVT,
     &            vvoo1,vvoo3)
      CALL RDVOVV(vovv1); vovv2 = 0.0d0
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NV,NO,NVVT,LVOVVT,LLOVVT,
     &             vovv1,vovv2)
      CALL CNTRCT ('N','N',NV,NO,NOVVT,A1,vovv2,vvoo3,A0,aak,NREP)
      deallocate(vovv1, vovv2)
      CALL DENOMVOTS (EPS,aak)
!
! do special contraction of D_ad and A_di over d! (loop from above)
! is the same as of D_ab and A_bi over b
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         arep = multb(desrep+nrep,irep,2)
         brep = irep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NV(brep)
         ioff1 = JPVVXX(desrep) + JJVV(arep,brep) + 1
         ioff2 = JJVO(brep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,dvvx(ioff1),m,aak(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvy(ioff1),m,aak(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvz(ioff1),m,aak(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,4) finished'

!
!-------------------------------------
!------  F_ai(2,5), a DOV contribution
!-------------------------------------
!
!  v_b>c,k>j  --->  -v_b>c,j>k  (RCA_NEGATE) vvoo1 --> -vvoo1
!  -v_b>c,j>k  --->  -v_bj,ck   (SRT1TT4) -vvoo1 --> vvoo3
!  negate vvoo1 again in order to restore original values.

      n1 = ivvoott(nrep+1) 
      CALL RCA_NEGATE(n1,vvoo1)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)
      CALL RCA_NEGATE(n1,vvoo1)

!  read Vaj,bi ---> Vai,bj (SRT16) (vovo1 ---> vovo2)
!  in special symmetry contraction loop. Buffer not yet allocated.

      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
        if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      allocate(vovo2(n1*rcw),stat=ialloc)
        if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      CALL getvovo(vovo1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)

!  contract Vai,bj (vovo2) and -v_bj,ck (vvoo3) to Aai,ck (vovo1)
!  Aai,ck  attains  <||>,<||>  BKC  structure

      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vovo2(OFF1),M,vvoo3(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO

!  backsort Aai,ck  <||>,<||>, vovo1 
!      to Bac,ik (vvoo3)  (1TT4 inverse)
!  Bac,ik again attains <|<|, |>|> structure.
!  divide Bac,ik by the energy denominator (denomvvoo).
!  results in b_a>c,i>k (vvoo3).

      vvoo3 = 0.0d0
      CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vovo1,vvoo3)
      CALL DENOMVVOO (eps,vvoo3,vvoo3)
      deallocate(vovo1,vovo2)

!  sort ba>c,i>k (vvoo3) to bac,i>k (vvoo4)
!  and bac,i>k (vvoo4) ---> bai,kc (vvoo3) (SRT20M)

      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoo3,vvoo4)
      CALL SRT20M(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &            JPVOOV,JJVO,JJOV,vvoo4,vvoo3)

!  contract bai,kc (vvoo3)  and Dkc over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JPVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,5) finished'
!
!-------------------------------------
!------  F_ai(2,6): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      allocate(vovo2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
!
!  sort the -v_ab,ji (vvoo1) ---> -v_ai,bj (vvoo3)
!  and the V*_bk,cj ---> V*_bj,ck (SRT 16)
!  then contract -v_ai,bj and V*_bj,ck to  A_ai,ck
!       VO,VO  x  VO,VO  --->  VO,VO <||>,<||>
!
      n1 = ivvoott(nrep+1) 
      CALL RCA_NEGATE(n1,vvoo1)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)
      CALL RCA_NEGATE(n1,vvoo1)

      CALL getvovo(vovo1)
      n1 = ivovo(nrep+1) 
      if(carith) call conjuga(n1,vovo1,1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vvoo3(OFF1),M,vovo2(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
!  A_ai,ck now in vovo1
!  Do inverse sort to B_ac,ik into vvoo3
!  Divide B_ac,ik by denominators and proceed from here
!  as in F_ai(2,5)
!
      CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vovo1,vvoo3)
      CALL DENOMVVOO (eps,vvoo3,vvoo3)
      deallocate(vovo1,vovo2)

!  sort ba>c,i>k (vvoo3) to bac,i>k (vvoo4)
!  and bac,i>k (vvoo4) ---> bai,kc (vvoo3) (SRT20M)

      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoo3,vvoo4)
      CALL SRT20M(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &            JPVOOV,JJVO,JJOV,vvoo4,vvoo3)

!  contract bai,kc (vvoo3)  and Dkc over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JPVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,6) finished'
!
!-------------------------------------
!------  F_ai(2,7): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      allocate(vovo2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
!
!  sort the v_cb,ji (vvoo1) ---> v_cj,bi (1TT4)
!  and the V_aj,ck ---> V_ak,cj (SRT 16)
!  contract V_ak,cj  x  v_cj,bi  to A_ak,bi
!
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)
      CALL getvovo(vovo1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vovo2(OFF1),M,vvoo3(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
! A_ak,bi is now in vovo1
! Backsort A_ak,bi ---> B_ab,ki (1TT4 inv)
! Divide by energy denominator and negate in order to
! obtain b_ab,ik (note lowercase letter!)
!
      vvoo3 = 0.0d0
      CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vovo1,vvoo3)
      CALL DENOMVVOO (eps,vvoo3,vvoo3)
      CALL RCA_NEGATE(ivvoott(nrep+1),vvoo3)
      deallocate(vovo1,vovo2)
!
! resort b_a>b,i>k to b_ab,i>k to C_ai,kb
!
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoo3,vvoo4)
      CALL SRT20M(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &            JPVOOV,JJVO,JJOV,vvoo4,vvoo3)

!  contract C_ai,kb (vvoo3)  and Dkb over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JPVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,7) finished'
!
!-------------------------------------
!------  F_ai(2,8): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      allocate(vovo2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
!
!  sort the v_ac,kj (vvoo1) ---> v_ak,cj (1TT4)
!  and the V*_ci,bj ---> V_cj,bk (SRT 16)
!  contract v_ak,cj  and   V*_cj,bi to A_ak,bi
!
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)
      CALL getvovo(vovo1)
      if(carith) CALL conjuga(ivovo(nrep+1),vovo1,1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vvoo3(OFF1),M,vovo2(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
! A_ak,bi is now in vovo1, and we can proceed exactly in the same way
! as for F_ai(2,7)!
!
! Backsort A_ak,bi ---> B_ab,ki (1TT4 inv)
! Divide by energy denominator and negate in order to
! obtain b_ab,ik (note lowercase letter!)
!
      vvoo3 = 0.0d0
      CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vovo1,vvoo3)
      CALL DENOMVVOO (eps,vvoo3,vvoo3)
      CALL RCA_NEGATE(ivvoott(nrep+1),vvoo3)
      deallocate(vovo1,vovo2)
!
! resort b_a>b,i>k to b_ab,i>k to C_ai,kb
!
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoo3,vvoo4)
      CALL SRT20M(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &            JPVOOV,JJVO,JJOV,vvoo4,vvoo3)

!  contract C_ai,kb (vvoo3)  and Dkb over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JPVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,8) finished'
!
!-------------------------------------
!------  F_ai(2,9): a DOV contribution
!-------------------------------------
!
! contraction of V_ab,cd and v_cd,ij is done in precontraction section
! of main part. We read in the VVVV x VVOO contracted integrals.
!
!
! Read the A_a>b,i>j to vvoop
!
      n1 = ivvoott(nrep+1)
      allocate(vvoop(n1*rcw))
      open(ioch,file=preconfn,form='unformatted',
     &     access='sequential',status='unknown')
      n1 = ivvoott(nrep+1)
      read(ioch) (vvoop(K),K=1,n1*RCW)
      close(ioch)
      CALL DENOMVVOO (eps,vvoop,vvoop)
!
! resort b_a>b,i>j to b_ab,i>j to C_ai,jb,
! contract with D_jb and store.
!
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoop,vvoo4)
      CALL SRT20M(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &            JPVOOV,JJVO,JJOV,vvoo4,vvoo3)
      deallocate(vvoop)

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JPVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,9) finished'
!
!-------------------------------------
!-----  F_ai(2,10): a DOV contribution
!-------------------------------------
!
!  contract v_ab,kl and V_kl,ij to A_ab,ij
!  from here we go on as in F_ai(2,9)
!
      n1=ioooott(nrep+1)
      allocate(oooo(n1*rcw))
      call getoooo(oooo)
      CALL CNTRCT('N','N',NVVT,NOOT,NOOT,A1,vvoo1,oooo,A0,vvoo3,NREP)
      deallocate(oooo)
!
      CALL DENOMVVOO (eps,vvoo3,vvoo3)
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoo3,vvoo4)
      CALL SRT20M(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &            JPVOOV,JJVO,JJOV,vvoo4,vvoo3)

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JPVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,10) finished'
!
! ----------------------------------------------------------
! --- from here all second order h/p contributions are -----
! ---      done and we print the resulting TMs       -------
! ----------------------------------------------------------
!
      iorder='0th-2nd order (ph only)'
 7744 CALL TM_GEN(ioch,desrep,isatctl,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
!
!  ----------------------------------------------------------
!  -------     2h2p second order contributions to F  --------
!  ----------------------------------------------------------
!
!  allocate required space, faibj arrays have to be complex!
!
      write(*,*) '------  calculating the 2h2p 2nd order contributions'
      n1=ABIJTT(desrep)  !ABIJTT is non-cumulative!
      write(*,*) 'Allocating F_aibj buffers of length',n1
      allocate(fabijtx(n1))
      allocate(fabijty(n1))
      allocate(fabijtz(n1))
      fabijtx = A0
      fabijty = A0
      fabijtz = A0


      n1=mvvoot(desrep)  !mvvoot is non-cumulative!
      write(*,*) 'Allocating three A_abi,k buffers of length',n1
      allocate(vabikx(n1*rcw))
      allocate(vabiky(n1*rcw))
      allocate(vabikz(n1*rcw))
      vabikx = 0.0d0
      vabiky = 0.0d0
      vabikz = 0.0d0

      n1=ABIJTT(desrep)  !ABIJTT is non-cumulative!
      write(*,*) 'Allocating one more F_aibj buffer of length',n1
      allocate(vvoop(n1*rcw))
      vvoop = 0.0d0
!
!---------------------------------------------------
!-----  F_ab,ij(2,A+B and antisymmetrization) ------
!---------------------------------------------------
!
!  v_ab,ki = -v_ab,ik (vvoo1)
!  sort -v_ab,ik to -v_abi,k  (vvoo3)
!  contract v_abi,k and D(x,y,z)kj  ---> Babi,j **NTS**
!  backsort and antisymmetrize Babi,j ---> Bab,ij.
!
      n1=ivvoott(nrep+1)
      CALL RCA_NEGATE(n1,vvoo1)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             vvoo1,vvoo3)
      CALL RCA_NEGATE(n1,vvoo1)
!
! apply proper symmetry loop due to nontotally symmetric D_oo contributions
! 
      do jrep = 1, nrep
        krep = multb(desrep+nrep,jrep,2)
        if(MULTB(krep,jrep,2).ne.desrep) stop 'munch:-('
        m = NVVOT(krep)
        n = NO(jrep)
        k = NO(krep)
        ioff1 = KVVOOT(krep) + 1
        ioff2 = JPOOXX(desrep) + JJOO(krep,irep) + 1
        ioff3 = JAVVOO(krep,jrep) + 1
!x
        call xgemm ('N','N',m,n,k,a1,vvoo3(ioff1),m,doox(ioff2),k,
     &              a0,vabikx(ioff3),m)
!y
        call xgemm ('N','N',m,n,k,a1,vvoo3(ioff1),m,dooy(ioff2),k,
     &              a0,vabiky(ioff3),m)
!z
        call xgemm ('N','N',m,n,k,a1,vvoo3(ioff1),m,dooz(ioff2),k,
     &              a0,vabikz(ioff3),m)
      enddo


      n1 = ABIJTT(desrep)
      write(*,*) 'ICEIL 1:',n1
      CALL SRT1T3_NTSI (NREP,desrep,MULTB,NVVT,NO,NO,NVVOT,
     &                  JAVVOO,KKVVOT,vabikx,vvoop,n1)
      call rcaxpy(n1,carith,pphasex,vvoop,fabijtx)
      CALL SRT1T3_NTSI (NREP,desrep,MULTB,NVVT,NO,NO,NVVOT,
     &                  JAVVOO,KKVVOT,vabiky,vvoop,n1)
      call rcaxpy(n1,carith,pphasey,vvoop,fabijty)
      CALL SRT1T3_NTSI (NREP,desrep,MULTB,NVVT,NO,NO,NVVOT,
     &                  JAVVOO,KKVVOT,vabikz,vvoop,n1)
      call rcaxpy(n1,carith,pphasez,vvoop,fabijtz)

      deallocate(vvoop)
      deallocate(vabikx)
      deallocate(vabiky)
      deallocate(vabikz)

      write(*,*) 'F_ab,ij(2,A+B) contributions finished.'
!
!---------------------------------------------------
!-----  F_ab,ij(2,C+D and antisymmetrization)  -----
!---------------------------------------------------
!
      n1=BVVOOT(desrep)
      write(*,*) 'Allocating three B_a,bij buffers of length',n1
      allocate(vabikx(n1*rcw)) !still called vabik? for simplicity
      allocate(vabiky(n1*rcw))
      allocate(vabikz(n1*rcw))
      vabikx = 0.0d0
      vabiky = 0.0d0
      vabikz = 0.0d0

      n1=ABIJTT(desrep)
      write(*,*) 'Allocating one more F_aibj buffer of length',n1
      allocate(vvoop(n1*rcw))
      vvoop = 0.0d0
!
!  sort v_cb,ij to v_c,bij (1T2)  -->  vvoo3 (dim = lvvoot(nrep+1)! )
!  contract Dac with v_c,bij to B_a,bij, keeping proper symmetry
!  backsort B_a,bij to B_ab,ij  including antisymmetrization.
!
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      do crep = 1, nrep
         arep = multb(desrep+nrep,crep,2)
         if(MULTB(arep,crep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NVOOT(crep)
         k = NV(crep)
         ioff1 = JPVVXX(desrep) + JJVV(arep,crep) + 1
         ioff2 = LVVOOT(crep) + 1
         ioff3 = JBVVOO(arep,crep) + 1
!x
         call xgemm ('N','N',m,n,k,a1,dvvx(ioff1),m,vvoo3(ioff2),k,
     &               a0,vabikx(ioff3),m)
!y
         call xgemm ('N','N',m,n,k,a1,dvvy(ioff1),m,vvoo3(ioff2),k,
     &               a0,vabiky(ioff3),m)
!z
         call xgemm ('N','N',m,n,k,a1,dvvz(ioff1),m,vvoo3(ioff2),k,
     &               a0,vabikz(ioff3),m)
      enddo

      n1 = ABIJTT(desrep)
      write(*,*) 'ICEIL 2:',n1

! debug: original call leading to boundary violations
!     CALL SRT1T2_NTSI (NREP,desrep,MULTB,NVVT,NV,NV,NOOT,
!    &                  JBVVOO,vabikx,vvoop,n1)
! debug end

      CALL SRT1T2_NTSI (NREP,desrep,MULTB,NVVT,NV,NV,NOOT,
     &                  JBVVOO,LLVOOT,vabikx,vvoop,n1)
      call rcaxpy(n1,carith,pphasex,vvoop,fabijtx)
      CALL SRT1T2_NTSI (NREP,desrep,MULTB,NVVT,NV,NV,NOOT,
     &                  JBVVOO,LLVOOT,vabiky,vvoop,n1)
      call rcaxpy(n1,carith,pphasey,vvoop,fabijty)
      CALL SRT1T2_NTSI (NREP,desrep,MULTB,NVVT,NV,NV,NOOT,
     &                  JBVVOO,LLVOOT,vabikz,vvoop,n1)
      call rcaxpy(n1,carith,pphasez,vvoop,fabijtz)


      deallocate(vvoop)
      deallocate(vabikx)
      deallocate(vabiky)
      deallocate(vabikz)
!
!  calculate total transition moment including all contribs.
!
      iorder='0th-2nd order (ph + 2p2h)'
      isatctl = ABIJTT(desrep)

      CALL TM_GEN(ioch,desrep,isatctl,iorder,
     &            LTR,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      deallocate(fabijtx,fabijty,fabijtz)
!______________________________________________________
!|
!|
!|
!|
!|
!|  END CALCULATIONS, Deallocate all remaining arrays
!|
!|
!|
!|
!|_____________________________________________________
!
!
      deallocate(abj,aak,aki,aac)

      deallocate(vvoo3,vvoo4)
      deallocate(vvoo1,vvoo2)

      deallocate(doox,dooy,dooz)
      deallocate(dovx,dovy,dovz)
      deallocate(dvox,dvoy,dvoz)
      deallocate(dvvx,dvvy,dvvz)

      deallocate(faipx,faipy,faipz)
      deallocate(faitx,faity,faitz)

      deallocate(pphase)

      deallocate(ABIJTT)
      deallocate(JAVVOO)
      deallocate(JBVVOO)
      deallocate(MVVOOT)
      deallocate(BVVOOT)

      deallocate(JPVOOV)
      deallocate(JPOOXX,JPVOXX,JPOVXX,JPVVXX)

      RETURN
      END
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
      real*8,allocatable,dimension(:)         :: evl
      real*8                                  :: ps
      real*8, parameter                       :: D0=0.0d0
      complex*16, parameter                   :: A0=(D0,D0)
      integer                                 :: i,k,ist
      complex*16                              :: ax,ay,az
      real*8, parameter                       :: autoev = 27.2113957D0
      integer                                 :: n2h2p
!
!  executable
!
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
        write(*,*) 'Found eigenvectors in symmetry',desrep
      endif
      write(*,*) 'CARITH/RCW:',carith,RCW
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

      Do k=1,nroots
        read(ioch) ist
        read(ioch) evl(k)
        read(ioch) ps
        write(*,'(A,X,I3,X,A,F10.4,A,X,F10.4,A,F10.4)') 
     &    'State',ist,':exc. energy',evl(k),'(au)',evl(k)*autoev,
     &    ' (eV) and PS',ps
        if(carith) then
          read(ioch)(ccoef(i),i=1,ladc)
        else
          read(ioch)(rcoef(i),i=1,ladc)
        endif
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
        write(*,'(I4,F14.7,2F15.9,A,2F15.9,A,2F15.9)') k,evl(k)*autoev,
     &             tmx(k),' | ',tmy(k),' | ',tmz(k)
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
     &                           pri)
!
      use interface_to_mpi
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!
!     Reads property integral files for operator A. Original version by LV.
!     The diagonal blocks (OO and VV) are always read.
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
!     The routine was originally devised by LV and modified/extended
!     by MP for the RELADC purposes.
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
      Logical        SUCC
      Character*8    NAMEA
      INTEGER        NSP
      INTEGER        IREPSPI(NSP,MXREP,2)
      INTEGER        dimpoo,dimpov,dimpvo,dimpvv
      REAL*8         PROP(2,NSP,NSP)
      REAL*8         AOO(*),AOV(*),AVO(*),AVV(*)
      Logical        pri
!
!---------------Local variables--------------------------------------
!
      Integer MDPROP
      CHARACTER*32 ACHAR
      REAL*8 APHASE(2)
      REAL*8 AVOI,AOVI,AMAXR,AMAXI,ANORM
      LOGICAL PHASE
      REAL*8 hdev
      COMPLEX*16 h1,h2
! --- flushing
      Character*1              :: co
      Character*8              :: filename
      Integer                  :: ixx
!
!---------------Executable code--------------------------------------
!
!
!     Read symmetry information (master node only) from MDPROP
!     The file is organized in a way that after a label identifying the
!     property a binary number stream without CRLF follows. If the label
!     is no match, we read til the next CRLF which corresponds to a skipping
!     of all the unwanted integrals.
!
      MDPROP = 15

      IF(dimpov.ne.dimpvo) CALL QUIT('Inconsistent dimensions!')

      IF (MYPROC.EQ.MASTER) THEN

      OPEN (MDPROP,FILE='MDPROP',FORM='UNFORMATTED')
      INOP = 0
    1 READ (MDPROP,ERR=10,END=11) ACHAR
      IF (ACHAR(1:8).NE.'********'.OR.ACHAR(25:32).NE.NAMEA) GOTO 1
      WRITE (6,1000) NAMEA
      READ (MDPROP) PROP
      CLOSE (MDPROP,STATUS='KEEP')

      write (6,*)
      write (6,*) ' *** Raw property integrals:'
      hdev = 0.0d0
      do i=1,nsp
        do j=i,nsp
!         write(6,'(2I3,2F16.10)') j,i,prop(1,j,i),prop(2,j,i)
          h1=dcmplx(prop(1,j,i),prop(2,j,i))
          h2=dcmplx(prop(1,i,j),prop(2,i,j))
          hdev = hdev + abs(h1 - dconjg(h2))
        enddo
      enddo
      write (6,*) ' *** Deviation from Hermiticity:',hdev
      write (6,*)

      GOTO 12
   10 INOP = 1
      GOTO 12
   11 INOP = 2
   12 CONTINUE

      ENDIF
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
        WRITE (6,*) ' Error reading property ',NAMEA,' on file MDPROP'
        WRITE (6,*) ' Error reading property integrals'
        SUCC = .false.
        return
      ENDIF
      IF (INOP.EQ.2) THEN
        WRITE (6,*) ' Property ',NAMEA,' not found on file MDPROP'
        WRITE (6,*) ' Property integrals missing'
        SUCC = .false.
        return
      ENDIF

      if(carith) then
        write(*,*) 'We have complex arithmetic. RCW=',RCW
      else
        write(*,*) 'We have real arithmetic. RCW=',RCW
      endif
!
      APHASE(1) = D1
      APHASE(2) = D0
      PHASE = .FALSE.

    2 AMAXI = D0
      AAMAX = 0
      IIMAX = 0
      AI = 1
      DO AIREP = 1, NREP
         DO IREP = 1, NREP
            AREP = MULTB (IREP,AIREP+NREP,2)
            AREP = MULTB (AREP,NREP+1,2)
            DO I = 1, NO(IREP)
               II = IREPSPI(I,IREP,1)
               DO A = 1, NV(AREP)
                  AA = IREPSPI(A,AREP,2)
                  IF (CARITH) THEN
                      AVO(RCW*AI-1) = PROP(1,AA,II)
                      AVO(RCW*AI)   = PROP(2,AA,II)
                      AOV(RCW*AI-1) = PROP(1,II,AA)
                      AOV(RCW*AI)   = PROP(2,II,AA)
                  ELSE
                      AVO(AI) = PROP(1,AA,II)*APHASE(1) -
     &                          PROP(2,AA,II)*APHASE(2)
                      AOV(AI) = PROP(1,II,AA)*APHASE(1) -
     &                          PROP(2,II,AA)*APHASE(2)
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
! ********* from here all integrals are available and processed
! ********  we flush them to disk
!
      co = NAMEA(1:1)
      if(carith) then         ! indicator for aphase that we have carith
        aphase(1) = 0.0d0
        aphase(2) = 0.0d0
      endif
      write(filename,'(A1,A7)') co,'PROPINT'
      write(*,*) 'Flushing properties to ',filename
      open(unit=MDPROP,file=filename,form='unformatted',
     &     access='sequential',status='unknown')
      write(MDPROP) dimpoo,dimpov,dimpvo,dimpvv
      write(MDPROP) aphase(1),aphase(2)
      write(MDPROP) (aoo(ixx),ixx=1,RCW*dimpoo)
      write(MDPROP) (aov(ixx),ixx=1,RCW*dimpov)
      write(MDPROP) (avo(ixx),ixx=1,RCW*dimpvo)
      write(MDPROP) (avv(ixx),ixx=1,RCW*dimpvv)
      close(MDPROP)
!
! ** if pri=T we print the integrals
!
      if(pri) then
        write(*,*) '---- integral type AOO ----'
        do i=1,dimpoo
          if(carith) then
            write(6,'(A,I5,X,2F16.10)') 'Index: ',i,AOO(2*i-1),AOO(2*i)
          else
            write(6,'(A,I5,X,F16.10)') 'Index: ',i,AOO(i)
          endif
        enddo
        write(*,*) '---- integral type AOV ----'
        do i=1,dimpov
          if(carith) then
            write(6,'(A,I5,X,2F16.10)') 'Index: ',i,AOV(2*i-1),AOV(2*i)
          else
            write(6,'(A,I5,X,F16.10)') 'Index: ',i,AOV(i)
          endif
        enddo
        write(*,*) '---- integral type AVO ----'
        do i=1,dimpvo
          if(carith) then
            write(6,'(A,I5,X,2F16.10)') 'Index: ',i,AVO(2*i-1),AVO(2*i)
          else
            write(6,'(A,I5,X,F16.10)') 'Index: ',i,AVO(i)
          endif
        enddo
        write(*,*) '---- integral type AVV ----'
        do i=1,dimpvv
          if(carith) then
            write(6,'(A,I5,X,2F16.10)') 'Index: ',i,AVV(2*i-1),AVV(2*i)
          else
            write(6,'(A,I5,X,F16.10)') 'Index: ',i,AVV(i)
          endif
        enddo
      endif

      RETURN

 1000 FORMAT (/' Read integral type ',A8,' created ',
     & A8,' storage info : ',A8)
 1010 FORMAT (/' Largest imaginary part of matrix element',2I5,F10.2)
 1020 FORMAT (/' Property ',A8,' scaled with phase factor ',2F10.4)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C      symmetry testing loops for NTS two-electron parts.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     DO KREP=1,NREP
!       ABIREP = MULTB(KREP,DESREP+NREP,2)
!       write(*,*) 'ABIREP/KREP/A*K/K*A/DESREP',
!    &    abirep,krep,multb(abirep,krep,1),
!    &    multb(krep,abirep,1),desrep
!         DO KJREP=1,NREP
!           JREP = MULTB(KREP,KJREP+NREP,1)
!           write(*,*) '  -  ABIREP/JREP/A*J/DESREP',
!    &      ABIREP,JREP,MULTB(ABIREP,JREP,1),DESREP
!         ENDDO
!     ENDDO
