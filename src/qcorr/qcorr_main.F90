***********************************************************************

      subroutine Q_correction(C,NBLOCK,IBLOCK,IBLTP,LUC,
     &                        LUSCR,ICISTR,IRC,WORK,relQC)
!     incoming variables:
!                         C      - scratch space for ci vector
!                         NBLOCK - # of ttsblocks
!                         IBLOCK - array containing information about each ttsblock
!                         IBLTP  - block types
!                         LUC    - internal file pointer for CI vector file
!                         LUSCR  - file pointer for scratch vector file
!                         ICISTR - currently set to 2 by default
!                         IRC:   - algebra control value: 1 (real) or 2 (complex)
      use qcorr_cfg
      use qcorr_interface

      implicit none

!     input
      integer, intent(in)   :: IBLOCK(8,NBLOCK),IBLTP(*)
      real(8), intent(inout):: C(*), WORK(*)
      integer, intent(in)   :: nblock, luc, luscr, irc, icistr
      logical, intent(in)   :: relQC
!     scratch
      integer, allocatable  :: astring(:,:)
      integer, allocatable  :: bstring(:,:)
      integer, allocatable  :: occupation_class(:,:)
      real(8), allocatable  :: ref_wavefunction_coeff(:)
      real(8), allocatable  :: e_corr(:), e_tot(:)
      integer               :: nstates_ci
      real(8)               :: q_corr, edci
      real(8)               :: Edc_1, Erdc_1, Epc_1, Emc_1, Epcp_1
      real(8)               :: Edc_2, Erdc_2, Epc_2, Emc_2, Epcp_2
      real(8)               :: reference_energy
      integer               :: iactp, ibctp, nbl_c, nael, nbel
      integer               :: ntest, lblk, iboff, ionem
      integer               :: irir, imk2_r, i

      call QENTER('QCORR')

!     initialize
      call init_qcorr(relQC)

!     print flag
      ntest = 000
      ntest = max(ntest,print_qcorr)

      allocate(ref_wavefunction_coeff(nref_qcorr))
      ref_wavefunction_coeff = 0.0d0 ! initialize all to 0.0d0

      LBLK = -1
      IBOFF = 0
      call REWINE(LUSCR,-1)
      call REWINE(LUC,-1)
      call itods(-1,1,-1,LUSCR)

!     loop over real and imaginary parts
      do IRIR = 1,IRC

        if (IRIR.EQ.2) then
!         skip EOV mark between real and imaginary part
          call IFRMDS(IONEM,1,-1,LUC)
        end if

        do IMK2_R = 1,NMS2VAL_q,1

          if(relQC)then
            IACTP = IST_FOR_DT_q(1,IMK2_R)
            IBCTP = IST_FOR_DT_q(2,IMK2_R)
          else
            IACTP = 1
            IBCTP = 2
          end if
*  Offset for IBLOCK(sym,type) information
          if (IMK2_R.eq.1) then
            IBOFF = 1
          else
            IBOFF = IBOFF + NBLK_MS2_q(IMK2_R-1)
          end if
*
          if(relQC)then 
            NBL_C = NBLK_MS2_q(IMK2_R)
          else
            NBL_C = NBLOCK
          end if

          call REWINE(LUSCR,-1)
          call COPNBLKD(LUC,LUSCR,C,NBL_C,0,LBLK)
          call ITODS(-1,1,-1,LUSCR)

          NAEL = NELEC_q(IACTP)
          NBEL = NELEC_q(IBCTP)
          allocate(astring(NAEL,MXNSTR_q))
          allocate(bstring(NBEL,MXNSTR_q))

*
          call get_Q_correction(C,LUSCR,NBL_C,
     &                          WORK(KNSTSO_q(IACTP)),
     &                          WORK(KNSTSO2_q(IBCTP)),
     &                          IBLOCK(1,IBOFF),
     &                          NAEL,NBEL,astring, bstring,
     &                          IBLTP,NSMST_q,
     &                          IACTP,IBCTP,ICISTR,ntest,
     &                          ref_wavefunction_qcorr, 
     &                          ref_wavefunction_coeff,nref_qcorr,
     &                          NORB_q,relQC, wunit_q)
*
!         NAEL/NBEL may be zero and therefore no allocation of string arrays...
          if(allocated(bstring)) deallocate(bstring)
          if(allocated(astring)) deallocate(astring)
        end do
      end do

      if(ntest > 2)then
        write(wunit_q,'(A)') '                                       '
        write(wunit_q,'(A)') ' Reference wavefunction & coefficient: ' 
        write(wunit_q,'(A)') ' ===================================== '

        do i = 1, nref_qcorr

          write(wunit_q,*) "alpha string: ",ref_wavefunction_qcorr(i,1)
          write(wunit_q,*) "beta string:  ",ref_wavefunction_qcorr(i,2)

          write(wunit_q,*)" coefficient of the reference wf ",i,
     &                    " in the CI: ",ref_wavefunction_coeff(i)
          write(wunit_q,*)" original coefficient of the reference wf: ",
     &                                   ref_wavefunction_coeff_qcorr(i)
        end do
      end if

      ! read data for +Q correction
      open(file="energies.CI",unit=10,status="old",
     &     form="unformatted",access="sequential")

      read(10) nstates_ci,reference_energy

      allocate(e_corr(nstates_ci))
      allocate(e_tot(nstates_ci))

      read(10) (e_tot(i),i=1,nstates_ci)
      close(10,status="keep")

      write(wunit_q,"(/A,/A,/A,/A,/A,/A)")
     &                   " The acronyms for each Q correction below"//
     &                   " refer to: ",
     &  " Edc:  Davidson correction              [Eq. (34)]",
     &  " Erdc: renormalized Davidson correction [Eq. (36)]",
     &  " Epc:  Pople correction                 [Eq. (38)]",
     &  " Epc`: modified Pople correction        [Eq. (39)]",
     &  " Emc : Meissner correction              [Eq. (41)]"

      do i = 1, nstates_ci
        e_corr(i) = e_tot(i) - reference_energy
        
        call E_Q(e_corr,ref_wavefunction_coeff,nref_qcorr,
     &           ref_wavefunction_coeff_qcorr,NACTEL_q,
     &           Edc_1, Erdc_1, Epc_1, Emc_1, Epcp_1,
     &           Edc_2, Erdc_2, Epc_2, Emc_2, Epcp_2,
     &           ntest, wunit_q)

      ! output of the q-correction       
        write(wunit_q,"(/A,i4)")     ' Q correction for state #', i
        write(wunit_q,"( A,1f16.8)") ' correlation energy: ', e_corr(i)
        write(wunit_q,"( A,1f16.8)") ' original CI energy: ', e_tot(i)

        write(wunit_q,"(//A)")  ' Q corrections with c0^2 according'//
     &                          ' to Eq. (43)'
 
        write(wunit_q,"(/A)") " Edc correction :"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Edc_1
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Edc_1
        write(wunit_q,"(A )") " ----------------"
 
        write(wunit_q,"(/A)") " Erdc correction:"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Erdc_1
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Erdc_1
        write(wunit_q,"(A )") " ----------------"
 
        write(wunit_q,"(/A)") " Epc correction :"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Epc_1
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Epc_1
        write(wunit_q,"(A )") " ----------------"
 
        write(wunit_q,"(/A)") " Epc` correction:"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Epcp_1
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Epcp_1
        write(wunit_q,"(A )") " ----------------"
 
        write(wunit_q,"(/A)") " Emc correction :"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Emc_1
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Emc_1
        write(wunit_q,"(A )") " ----------------"
 
!       now using eq. (44) for c0^2: c0^2 == <ref|ci>^2 == ( sum\limits_{p=1}^{nref} c_p^{(0)} * c_p )^2 
!       where c_p^{(0)} are the coefficients in the normalized 0th-order reference (MCSCF) wave function

        write(wunit_q,"(//A)")  ' Q corrections with c0^2 according'//
     &                          ' to Eq. (44)'
!       write(wunit_q,"(/A)")  " Q corrections according to Eq. (44)"//
!    &                   " where c0^2 == <ref|ci>^2 == "//
!    &                   "( sum\limits_{p=1}^{nref} c_p^{(0)} * c_p )^2"
!       write(wunit_q,"(A/)") " Note: c_p^{(0)} are the reference"//
!    &                  " coefficients in the"//
!    &                  " normalized 0th-order reference (MCSCF) wave"//
!    &                  " function!"
 
        write(wunit_q,"(/A)") " Edc correction :"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Edc_2
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Edc_2
        write(wunit_q,"(A )") " ----------------"
 
        write(wunit_q,"(/A)") " Erdc correction:"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Erdc_2
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Erdc_2
        write(wunit_q,"(A )") " ----------------"
 
        write(wunit_q,"(/A)") " Epc correction :"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Epc_2
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Epc_2
        write(wunit_q,"(A )") " ----------------"
 
        write(wunit_q,"(/A)") " Epc` correction:"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Epcp_2
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Epcp_2
        write(wunit_q,"(A )") " ----------------"
 
        write(wunit_q,"(/A)") " Emc correction :"
        write(wunit_q,"(A )") " ----------------"
        write(wunit_q,"(A,1f16.8)") " Q correction:",Emc_2
        write(wunit_q,"( A,1f16.8)")" CI+Q energy: ", e_tot(i) + Emc_2
        write(wunit_q,"(A )") " ----------------"

      end do


      deallocate(e_corr)
      deallocate(e_tot)

      deallocate(ref_wavefunction_qcorr)
      deallocate(ref_wavefunction_coeff)
      deallocate(ref_wavefunction_coeff_qcorr)

!     finalize
      call exit_qcorr(relQC)

      call QEXIT('QCORR')

      end subroutine Q_correction
***********************************************************************

      subroutine get_Q_correction(C,LUC,NBL_C,
     &                            NSSOA,
     &                            NSSOB,
     &                            IBLOCK,
     &                            NAEL,NBEL,IASTR,IBSTR,
     &                            IBLTP,NSMST,
     &                            IASPGPTP,IBSPGPTP,ICISTR,IPRNT,
     &                            ref_wavefunction,
     &                            ref_wavefunction_coeff,nref,
     &                            NORB,relQC, wunit_q)
***********************************************************************

      implicit none

      real(8), intent(inout)          :: C(*)
      integer, intent(inout)          :: IASTR(NAEL,*), IBSTR(NBEL,*)
      integer, intent(in)             :: luc, nbl_c, nsmst
      integer, intent(in)             :: nael, nbel, icistr, iprnt
      integer, intent(in)             :: NSSOA(NSMST,*), NSSOB(NSMST,*)
      integer, intent(in)             :: IBLTP(*)
      integer, intent(in)             :: IBLOCK(8,NBL_C)
      integer, intent(in)             :: IASPGPTP, IBSPGPTP
      integer, intent(in)             :: norb
      integer, intent(in)             :: nref, wunit_q
      character(len=200), intent(in)  :: ref_wavefunction(nref,2)
      logical, intent(in)             :: relQC

!     output
      real(8) , intent(inout)         :: ref_wavefunction_coeff(nref)
     
!     scratch
      character (len=4),  allocatable :: conf_string_a(:)
      character (len=4),  allocatable :: conf_string_b(:)
      character (len=10), allocatable :: stringa(:)
      character (len=10), allocatable :: stringb(:)
      character (len=200), allocatable:: ref_wavefunction_tmp(:,:)
      integer                         :: iel
      integer                         :: ia, iabas, ib, ibbas
      integer                         :: iampack, imzero, idet, idum
      integer                         :: ilena, ilena2, ilenb, ilenb2
      integer                         :: iasm, iatp, ibsm, ibtp
      integer                         :: iref_tmp, irestr_ana, jblock
      integer                         :: minia, nastr1, nbstr1, nia, nib

      allocate(ref_wavefunction_tmp(nref,2))

      if( ICISTR .ge. 2 ) call REWINE(LUC,-1)
      
      IDET     = 0
      do JBLOCK = 1, NBL_C
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)

!        Obtain alpha/beta strings of sym IASM/IBSM and type IATP/IBTP
         if(relQC)then 
#ifdef PRG_DIRAC
           call GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,
     &                                NASTR1,IASTR)
           call GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,
     &                                NBSTR1,IBSTR)
#endif
         else

           IDUM = 0
           call GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                            NORB,0,IDUM,IDUM)
           call GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                            NORB,0,IDUM,IDUM)
         end if

         if(IBLTP(IASM).eq.2) then
            IRESTR_ANA = 1
         else
            IRESTR_ANA = 0
         end if
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
         IMZERO = 0
         if(ICISTR.ge.2) then
!          read in a Type-Type-symmetry block
           call IFRMDS(IDET,1,-1,LUC)
           call FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
           IDET = 0
         end if
         if(IMZERO.ne.1) then
            IBBAS = 1
            IABAS = 1
            do IB = IBBAS,IBBAS+NIB-1                                  ! loop over beta strings
               if(IRESTR_ANA.eq.1.and.IATP.eq.IBTP) then
                  MINIA = IB - IBBAS + IABAS
               else
                  MINIA = IABAS
               end if
               do  IA = MINIA,IABAS+NIA-1                              ! loop over alpha strings 
*
                    IDET = IDET + 1
                   
                    if(iprnt > 10)then
                      write(wunit_q,*) 'idet, c',idet,c(idet),nael,nbel
                      write(wunit_q,"(60i4)")          ! save strings into array
     &                (iastr(iel,ia),iel=1,nael)
                      write(wunit_q,"(60i4)") 
     &                (ibstr(iel,ib),iel=1,nbel)
                    end if

                    write(ref_wavefunction_tmp(1,1),"(60i4)")          ! save strings into array
     &              (iastr(iel,ia),iel=1,nael)                         
                    write(ref_wavefunction_tmp(1,2),"(60i4)") 
     &              (ibstr(iel,ib),iel=1,nbel)

                    ilena   = len(trim(ref_wavefunction_tmp(1,1)))      
                    ilenb   = len(trim(ref_wavefunction_tmp(1,2)))      
                                                                        
                    iref_tmp = 0

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


               end do
            end do
         end if
      end do

      end subroutine get_Q_correction

***********************************************************************

      subroutine E_Q(e_corr,ref_wavefunction_coeff,nref,
     &               ref_wavefunction_coeff_i,N,
     &               Edc_1, Edrc_1, Epc_1, Emc_1, Epcp_1,
     &               Edc_2, Edrc_2, Epc_2, Emc_2, Epcp_2,
     &               iprnt, wunit_q)
*
*
* Davidson Correction
*
* Equations from: Chemical Reviews 2012, 112 page 216-218
*
* Edc  = eq. 34 p. 116
* Edrc = eq. 36 p. 116
* Epc  = eq. 38 p. 116
* Emc  = eq. 41 p. 116
*
* c01  = eq. 43 p. 117
* c02  = eq. 44 p. 117
*
***********************************************************************
*
      implicit none

      real*8,  intent(in)  :: e_corr
      real*8,  intent(in)  :: ref_wavefunction_coeff(nref)
      real*8,  intent(in)  :: ref_wavefunction_coeff_i(nref)
      real*8,  intent(out) :: Edc_1, Edrc_1, Epc_1, Emc_1, Epcp_1
      real*8,  intent(out) :: Edc_2, Edrc_2, Epc_2, Emc_2, Epcp_2
      real*8,  intent(in)  :: N
      integer, intent(in)  :: nref
      integer, intent(in)  :: iprnt, wunit_q

      real*8  :: c01    = 0      ! c01    = sum ( cp ) **2
      real*8  :: c02    = 0      ! c02    = (sum ( cref * cp )) **2
      real*8  :: c02_2  = 0      ! c02_2  = c02 * c02, entspricht c01
      real*8  :: teta_1 = 0      ! teta_1 = arccos( sqrt(c01) )
      real*8  :: teta_2 = 0      ! teta_2 = arccos( sqrt(c02_2))
      integer :: i      = 0
      

      do i = 1, nref
      c01 = c01 + ref_wavefunction_coeff(i) ** 2
      end do

      teta_1 = acos(sqrt(c01))

      do i = 1, nref
      c02 = c02 + ref_wavefunction_coeff(i)*ref_wavefunction_coeff_i(i)
      end do

      c02_2  = c02 * c02
      teta_2 = acos(sqrt(c02_2))

      if(iprnt > 10)then
        write(wunit_q,*) "c01",         c01
        write(wunit_q,*) "sqrt (c01)",  sqrt(c01)
        write(wunit_q,*) "c02",         c02
        write(wunit_q,*) "c02_2",       c02_2
      end if
     
      Edc_1  = (1 - c01) * e_corr            
      Edrc_1 = ((1-c01) / c01) * e_corr
      Epc_1  = (sqrt(N **2 + 2*N * (tan(2*teta_1))**2) - N) / 
     &         (2*(1/(2*teta_1)) - 1) * e_corr
      Epcp_1  = (1 - (2/N)) * ((1-c01)/c01) * e_corr

      if(N > 3.0d0) Emc_1 = ((N - 2) * (N - 3) / (N * (N - 1))) *
     &                      ((1 - c01 **2) / c01 **2) * e_corr  

      Edc_2  = (1 - c02_2) * e_corr
      Edrc_2 = ((1-c02_2) / c02_2) * e_corr
      Epc_2  = (sqrt(N **2 + 2*N * (tan(2*teta_2))**2) - N) / 
     &         (2*(1/(2*teta_2)) - 1) * e_corr
      Epcp_2  = (1 - (2/N)) * ((1-c02_2)/c02_2) * e_corr

      if(N > 3.0d0) Emc_2 = ((N - 2) * (N - 3) / (N * (N - 1))) *
     &                      ((1 - c02_2 **2) / c02_2 **2) * e_corr  
      
      end subroutine E_Q



