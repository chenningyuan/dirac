***********************************************************************

      SUBROUTINE Q_correction(C,NBLOCK,IBLOCK,IBLTP,LUC,
     &                        LUSCR,ICISTR,IRC)
*
*
* Analyze CI vector
*
* Jeppe Olsen, August 1995
* Driven By IBLOCK, May 1997
*           String occupations added, Feb. 98
*
*
*''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
* Modified for relativistic CI calculations
*
*        Timo Fleig, Feb. 2000
*''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!     in:
!         C: scratch space for ci vector
!         NBLOCK: # of ttsblocks
!         IBLOCK: array containing information about each ttsblock
!         IBLTP: ?
!         LUC: int filepointer
!         LUSCR: filepointer for scrathc file
!         ICISTR: 2
!         IRC: 1 (real) 2 (complex)
*
      use memory_allocator
      use symmetry_setup_krci
      use mospinor_info
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "strbas.inc"
#include "cicisp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "cands.inc"
#include "parluci.h"
#include "dgroup.h"
*
      DIMENSION IBLOCK(8,NBLOCK),IBLTP(*)
      dimension C(*)
      integer               :: mj_statistics(-32:32)
      integer, allocatable  :: mat1(:,:)
      integer, allocatable  :: astring(:,:)
      integer, allocatable  :: bstring(:,:)
      integer, allocatable  :: occupation_class(:,:)
      real(8), allocatable  :: mat2(:,:)
      character (len=2)     :: ferm_lab
      character (len=1)     :: inv_symbol
      character (len=200), allocatable   :: ref_wavefunction(:,:)
      real(8), allocatable  :: ref_wavefunction_coeff(:)
      real(8), allocatable  :: ref_wavefunction_coeff_i(:)       ! coeff von der ersten ci_iteration!!!
      integer :: nref = 1 ! for HF nref = 1, TODO: input parameter!!!
      integer :: nstates_ci
      real(8), allocatable :: e_corr(:), e_tot(:)
      real(8) :: q_corr, edci
      real(8) :: N = 8.0D0
*
      CALL QENTER('ANACI')
      allocate(ref_wavefunction(10,2))
      allocate(ref_wavefunction_coeff(10))
      allocate(ref_wavefunction_coeff_i(10))
      ref_wavefunction_coeff      = 0.0d0 ! initialize all to 0.0d0
      ref_wavefunction_coeff_i    = 0.0d0 ! initialize all to 0.0d0
      ref_wavefunction_coeff_i(1) = 1.0d0
      write(ref_wavefunction(1,1),"(4i4)"),1,2,3,4
      write(ref_wavefunction(1,2),"(4i4)"),1,2,3,4

*
** Specifications of internal space
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRWFC)
*
      LBLK = -1
      IBOFF = 0
      call REWINE(LUSCR,-1)
      call REWINE(LUC,-1)
      call itods(-1,1,-1,LUSCR)
      inv_symbol = ' '
*
* Loop over real and imaginary parts
      do IRIR = 1,IRC,1                                      !loop über call

        if (NTEST.ge.50) write(6,*) 'right r/i is  ',IRIR
        if (IRIR.EQ.2) then
*
*. Skip EOV mark between real and imaginary part
          call IFRMDS(IONEM,1,-1,LUC)
        end if
        do IMK2_R = 1,NMS2VAL,1
          MK2_R = MS2VAL(IMK2_R)
          IACTP = IST_FOR_DT(1,IMK2_R)
          IBCTP = IST_FOR_DT(2,IMK2_R)
*
*  Offset for IBLOCK(sym,type) information
          if (IMK2_R.eq.1) then
            IBOFF = 1
          else
            IBOFF = IBOFF + NBLK_MS2(IMK2_R-1)
          end if
*
          call REWINE(LUSCR,-1)
          NBL_C = NBLK_MS2(IMK2_R)
          CALL COPNBLKD(LUC,LUSCR,C,NBL_C,0,LBLK)
          CALL ITODS(-1,1,-1,LUSCR)
!     &         ' Contrib. Kramers proj., real/im  ',MK2_R,IRIR,
!     &         ' ----------------------------------------------'
          NAEL = NELEC(IACTP)
          NBEL = NELEC(IBCTP)
          call alloc(astring,NAEL,MXNSTR, id='astring')
          call alloc(bstring,NBEL,MXNSTR, id='bstring')

*
          NOCTPA   = NOCTYP(IACTP)
          NOCTPB   = NOCTYP(IBCTP)
*
          IOCTPA   = IBSPGPFTP(IACTP)
          IOCTPB   = IBSPGPFTP(IBCTP)
*
          IASPGPTP = IACTP
          IBSPGPTP = IBCTP
*
*. Number of occupation classes
          IWAY = 1
          NEL = NAEL + NBEL
          CALL OCCLS_REL(IWAY,NOCCLS,IOCCLS,NEL,NGAS,
     &                  IGSOCC(1,1),IGSOCC(1,2))
*
*. and then the occupation classes
          call alloc(occupation_class,ngas,noccls,id='occ_class')
          call alloc(mat1,10,NOCCLS, id='mat1')
          call alloc(mat2,10,NOCCLS, id='mat2')
          IWAY = 2
          CALL OCCLS_REL(IWAY,NOCCLS,occupation_class,NEL,NGAS,
     &                   IGSOCC(1,1),IGSOCC(1,2))

*
*. Occupation of strings of given sym and supergroup
          CALL get_Q_correction(C,LUSCR,NBL_C,WORK(KNSTSO(IACTP)),
     &                     WORK(KNSTSO2(IBCTP)),
     &                     NOCTPA,NOCTPB,MXPNGAS,IOCTPA,IOCTPB,
     &                     NBLOCK,IBLOCK(1,IBOFF),
     &                     THRES,MAXTRM,NAEL,NBEL,
     &                     astring, bstring,
     &                     IBLTP,NSMST,
     &                     mat1, mat2,
     &                     NELFSPGP,
     &                     NOCCLS,NGAS,occupation_class,
     &                     IASPGPTP,IBSPGPTP,
     &                     ICISTR,NTOOB,mj_statistics,inv_symbol,
     &                     IPRWFC, ref_wavefunction, 
     &                     ref_wavefunction_coeff,nref,
     &                     ref_wavefunction_coeff_i,N)
*
          call dealloc(mat2)
          call dealloc(mat1)
!         NAEL/NBEL may be zero and therefore no allocation of string arrays...
          if(allocated(bstring)) call dealloc(bstring)
          if(allocated(astring)) call dealloc(astring)
          call dealloc(occupation_class)
        end do
      end do

      !nref = 1       ! for HF nref = 1 TODO: input parameter!!!
      do i = 1, nref
      
      write(6,'(A)') '                                       '
      write(6,'(A)') ' Reference wavefunction & coefficient: ' 
      write(6,'(A)') ' ===================================== '

      print *, "alpha ref",ref_wavefunction(1,1)
      print *, "beta ref",ref_wavefunction(1,2)

      print *,"coefficient of reference wf #",i," == ",
     &         ref_wavefunction_coeff(i)
      print *,"coefficient of HF reference wf == ",
     &         ref_wavefunction_coeff_i(1)
      end do

      ! for +Q correction
      open(file="energies.CI",unit=10,status="old",
     &     form="unformatted",access="sequential")

      read(10) nstates_ci,reference_energy

      allocate(e_corr(nstates_ci))
      allocate(e_tot(nstates_ci))

      read(10) (e_tot(i),i=1,nstates_ci)
      close(10,status="keep")

      do i = 1, nstates_ci
        e_corr(i) = e_tot(i) - reference_energy
!        q_corr    = edc(e_corr(i),ref_wavefunction_coeff,
!     &                  nref)
        
        CALL E_Q(e_corr,ref_wavefunction_coeff,nref,
     &           ref_wavefunction_coeff_i,N,Edc_1,
     &           Edrc_1, Epc_1)

        write(6,"(/A,1f10.6,/A,1f10.6/)") " original ecorr = ",
     &        e_corr(i), " +Q correction Edc = ",Edc_1
        print *, "final energies: original ",e_tot(i), 
     &           "+Q corrected =",e_tot(i) + Edc_1
 
        print *, "original ecorr =",e_corr(i), "+Q correction Edrc_1=",
     &           Edrc_1
        print *, "final energies: original ",e_tot(i), 
     &           "+Q corrected =",e_tot(i) + Edrc_1
 
        print *, "original ecorr =",e_corr(i), "+Q correction Epc_1=",
     &           Epc_1
        print *, "final energies: original ",e_tot(i), 
     &           "+Q corrected =",e_tot(i) + Epc_1 
      end do


      deallocate(e_corr)
      deallocate(e_tot)

      deallocate(ref_wavefunction)
      deallocate(ref_wavefunction_coeff)
      deallocate(ref_wavefunction_coeff_i)

!
 
      CALL QEXIT('ANACI')
*
      END
***********************************************************************

      SUBROUTINE get_Q_correction(C,LUC,NBL_C,NSSOA,NSSOB,
     &                       NOCTPA,NOCTPB,
     &                       MXPNGAS,IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &                       THRES,MAXTRM,NAEL,NBEL,
     &                       IASTR,IBSTR,IBLTP,NSMST,
     &                       NCPMT,WCPMT,NELFSPGP,NOCCLS,
     &                       NGAS,IOCCLS,
     &                       IASPGPTP,IBSPGPTP,
     &                       ICISTR,NTOOB,mj_statistics,inv_symbol,
     &                       IPRNT, ref_wavefunction,
     &                       ref_wavefunction_coeff,nref,
     &                       ref_wavefunction_coeff_i,N)
*
***********************************************************************
*
      use mospinor_info
*  Orbitals
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION C(*)
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IBLTP(*)
      DIMENSION NELFSPGP(MXPNGAS,*)
      DIMENSION IOCCLS(NGAS,NOCCLS)
*. Specific input
      DIMENSION IBLOCK(8,NBL_C)
      character(len=200), intent(in)  :: ref_wavefunction(10,2)

*. Output
      DIMENSION NCPMT(10,NOCCLS)
      DIMENSION WCPMT(10,NOCCLS)
      real(8) , intent(inout) :: ref_wavefunction_coeff(10)
      real(8) , intent(inout) :: ref_wavefunction_coeff_i(10)       ! coeff von der ersten ci_iteration!!!
#include "dcborb.h"
#include "dgroup.h"
#include "maxash.h"
C
#include "maxorb.h"
#include "dcbidx.h"
!     scratch
      integer, allocatable              :: mjorb_vec(:)
      integer, intent(inout)            :: mj_statistics(-32:32)
      character (len=1) , intent(inout) :: inv_symbol
      character (len=1)                 :: inv_symbol_orb
      character (len=4),  allocatable   :: conf_string_a(:)
      character (len=4),  allocatable   :: conf_string_b(:)
      character (len=10), allocatable   :: stringa(:)
      character (len=10), allocatable   :: stringb(:)
      character (len=200), allocatable  :: ref_wavefunction_tmp(:,:)
      logical fndlab, ex, ex_open

      allocate(ref_wavefunction_tmp(10,2))

*
      NTESTL = 0000
      NTEST = max(IPRNT,NTESTL)
*
      MINPRT    = 0
      ITRM      = 0
      IDET      = 0
      IIDET     = 0
      ILOOP     = 0
      NCIVAR    = 0
      CNORM     = 0.0D0

      IF( ICISTR .GE. 2 ) CALL REWINE(LUC,-1)
      IIDET = 0
      ILOOP = ILOOP + 1
      
      IDET     = 0
*
*
      DO JBLOCK = 1, NBL_C
*
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)

*. Obtain alpha strings of sym IASM and type IATP
         CALL GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,
     &                              NASTR1,IASTR)
*. Obtain beta strings of sym IBSM and type IBTP
         CALL GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,
     &                              NBSTR1,IBSTR)
*
         IF(IBLTP(IASM).EQ.2) THEN
            IRESTR_ANA = 1
         ELSE
            IRESTR_ANA = 0
         END IF
*
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
         IMZERO = 0
         if (ICISTR.ge.2) then
*. Read in a Type-Type-symmetry block
*? Perhaps use GSTTBL here?
            CALL IFRMDS(IDET,1,-1,LUC)
            CALL FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
            IDET = 0
         END IF
!       write(6,*) ' IMZERO in gasanas:  ',IMZERO
!        if(linear) IMZERO = 0
         IF(IMZERO.NE.1) THEN
*
            IBBAS = 1
            IABAS = 1
*
            DO IB = IBBAS,IBBAS+NIB-1                                  ! loop über beta strings
               IF(IRESTR_ANA.EQ.1.AND.IATP.EQ.IBTP) THEN
                  MINIA = IB - IBBAS + IABAS
               ELSE
                  MINIA = IABAS
               END IF
               DO  IA = MINIA,IABAS+NIA-1                              ! loop über alpha strings 
*
                  IF (ILOOP .EQ. 1 ) NCIVAR = NCIVAR + 1
                    IDET = IDET + 1
                    ITRM = ITRM + 1
                    IIDET = IIDET + 1
                    IF( ITRM .LE. MAXTRM ) THEN
                       CNORM = CNORM + C(IDET) ** 2
                   

                    write(ref_wavefunction_tmp(1,1),"(60i4)")          ! save strings into array
     &              (iastr(iel,ia),iel=1,nael)                         
                    write(ref_wavefunction_tmp(1,2),"(60i4)") 
     &              (ibstr(iel,ib),iel=1,nael)

                    ilena   = len(trim(ref_wavefunction_tmp(1,1)))      
                    ilenb   = len(trim(ref_wavefunction_tmp(1,2)))      
                                                                        
                    iref_tmp = 0
                    !nref     = 1 ! TODO: schould be defined in first subroutine
                    do 

                      iref_tmp = iref_tmp + 1                           ! exit if iref_tmp > nref
                      if(iref_tmp > nref) exit

                      ilena2  = len(trim(ref_wavefunction(iref_tmp,1)))
                      ilenb2  = len(trim(ref_wavefunction(iref_tmp,2)))

                    if(ref_wavefunction_tmp(1,1)(1:ilena) ==            ! compare alpha strings 
     &                 ref_wavefunction(iref_tmp,1)(1:ilena2) .and.
     &                 ref_wavefunction_tmp(1,2)(1:ilenb) ==            ! compare beta strings 
     &                 ref_wavefunction(iref_tmp,2)(1:ilenb2)
     &                )then 
                       ref_wavefunction_coeff(iref_tmp) = C(IDET)       ! if true, save coeff
   
                    exit
                    end if
                    end do


                    END IF
               END DO
*     ^ End of loop over alpha strings
            END DO
*     ^ End of loop over beta strings
         END IF
*     ^ End of if statement for nonvanishing blocks
      END DO
*     ^ End of loop over blocks
!      IF(IIDET .EQ. 0 ) WRITE(6,*) '   ( no coefficients )'
*
!      WRITE(6,'(/A,E15.8/)')
!     &     '  Norm of printed CI vector .. ', CNORM                !ausschreiben des vektors



*
*.Size of CI coefficients
*
*
      IDET = 0
      IF(ICISTR .GE. 2 ) CALL REWINE(LUC,-1)
      CALL IZERO(NCPMT,10*NOCCLS)
      CALL DZERO(WCPMT,10*NOCCLS)
      DO JBLOCK = 1, NBL_C
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)
*
         IF(IBLTP(IASM).EQ.2) THEN
            IRESTR_ANA = 1
         ELSE
            IRESTR_ANA = 0
         END IF
*. Occupation class corresponding to given occupation
         JOCCLS = 0
         DO JJOCCLS = 1, NOCCLS
            IM_THE_ONE = 1
            DO IGAS = 1, NGAS
               IF(NELFSPGP(IGAS,IATP-1+IOCTPA)+
     &          NELFSPGP(IGAS,IBTP-1+IOCTPB).NE.IOCCLS(IGAS,JJOCCLS))
     &          IM_THE_ONE = 0
            END DO
            IF(IM_THE_ONE .EQ. 1 ) JOCCLS = JJOCCLS
         END DO
*
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
*
         IMZERO = 0
         IF( ICISTR.GE.2 ) THEN
*. Read in a Type-Type-symmetry block                               
            CALL IFRMDS(IDET,1,-1,LUC)
            CALL FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
            IDET = 0
         END IF
         IF(IMZERO.EQ.0) THEN

            DO IB = IBBAS,IBBAS+NIB-1                                  ! 2. do-loop ber a und b strings
               IF(IRESTR_ANA.EQ.1.AND.IATP.EQ.IBTP) THEN
                  MINIA = IB - IBBAS + IABAS
               ELSE
                  MINIA = IABAS
               END IF
               DO IA = MINIA,IABAS+NIA-1
                  IDET = IDET + 1
                  DO IPOT = 1, ipot_max
                     IF(10.0D0 ** (-IPOT+1).GE.ABS(C(IDET)).AND.
     &                 ABS(C(IDET)).GT. 10.0D0 ** ( - IPOT )) THEN
                       NCPMT(IPOT,JOCCLS)=
     &                       NCPMT(IPOT,JOCCLS)+ 1
                       WCPMT(IPOT,JOCCLS)=
     &                       WCPMT(IPOT,JOCCLS)+ C(IDET) ** 2
                     END IF
                  END DO
*     ^ End of loop over powers of ten
               END DO
*     ^ End of loop over alpha strings
            END DO
*     ^ End of loop over beta strings
         END IF
*     ^ End of test for novanishing blocks
      END DO
*     ^ End of loop over blocks
*

      END

***********************************************************************

      SUBROUTINE E_Q(e_corr,ref_wavefunction_coeff,nref,
     &               ref_wavefunction_coeff_i,N,Edc_1,
     &               Edrc_1, Epc_1)
*
*
***********************************************************************
*
      IMPLICIT NONE

      REAL*8,  intent(in)  :: e_corr
      real*8,  intent(in)  :: ref_wavefunction_coeff(nref)
      real*8,  intent(in)  :: ref_wavefunction_coeff_i(nref)
      real*8,  intent(out) :: Edc_1, Edrc_1, Epc_1
      real(8), intent(in)  :: N
      integer, intent(in)  :: nref

      REAL*8  :: c01 = 0         ! c01 = sum ( cp ) **2
      REAL*8  :: c02 = 0         ! co2 = (sum ( cref * cp )) **2
      REAL*8  :: c02_2 = 0
      real*8  :: Edc_2, Edrc_2, Epc_2
      integer :: i
      integer :: j
      

      do i = 1, nref
      c01 = c01 + ref_wavefunction_coeff(i) ** 2
      end do

      do j = 1, nref
      c02 = c02 + ref_wavefunction_coeff(j)*ref_wavefunction_coeff_i(j)
      end do

      c02_2 = c02 * c02
      print *, "c02_2 is ",c02_2
     
      Edc_1  = (1 - c01) * e_corr           ! TODO: Epc nicht pc', Emc, E_2, doc, 
      Edrc_1 = ((1-c01) / c01) * e_corr
      Epc_1  = (1 - (2/N)) * ((1 - c01) / c01) * e_corr    

!     if(N > 3.0d0) Emc = ...    

      Edc_2  = (1 - c02_2) * e_corr
      Edrc_2 = ((1-c02_2) / c02_2) * e_corr
      Epc_2  = (1 - (2/N)) * ((1 - c02_2) / c02_2) * e_corr

      end



