* Driver routine for distribution of DIRAC integrals to
* double group GAS lists
*   Routine using excitation operator sequence for integral
*   storage definition (SIGDEN route)
*
*----------------------------------------------
*  Timo Fleig, Mar. 2001 (Ellen Tabea is here!)
*  Jeppe Olsen modifies it. (Apr. 2001)
*  In complete form, June 2001, Timo Fleig
*----------------------------------------------
*
      subroutine dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,ISOBEX,T,
     &                            ISTR_CA,ISTR_CB,ISTR_AA,ISTR_AB,
     &                            INDEX,NASHT,F1,F2,IBOSYM,
     &                            IBSPOBEX_TP,IHINDX,LABEL,CIRUN,IPRNT)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "ctcc.inc"
#include "symm.inc"
#include "cgas.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "clunit.inc"
*
      character*6 CIRUN
      dimension NICA_STR(MXPNGAS),NICB_STR(MXPNGAS),
     &          NIAA_STR(MXPNGAS),NIAB_STR(MXPNGAS)
      dimension ISOBEX(NGAS,4,NSPOBEX_TP),INDEX(4,NINT_CC),
     &          IBSPOBEX_TP(NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          LABEL(4,NSPOBEX_TP)
      dimension I_ADD(4)
*
* Group notation of the strings:
      dimension IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS),
     &          IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
      dimension T(NINT_CC)
*
* Space for the construction of SUBSTRINGS
* dimensioned outside as integers MX_ST_TSOSO_BLK_MX
      dimension ISTR_CA(*),ISTR_CB(*),ISTR_AA(*),ISTR_AB(*)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
* Double group type
      if (DOUGRP.le.3) IDGTYP = 1
      if (DOUGRP.ge.4.and.DOUGRP.le.6) IDGTYP = 2
      if (DOUGRP.ge.7.and.DOUGRP.le.8) IDGTYP = 3
*
* Check for inversion symmetry of group
      INVSYM = 0
      if (DOUGRP.eq.7.or.DOUGRP.eq.4.or.DOUGRP.eq.1) INVSYM = 1
*
* Get reordering info depending on environment
      if (INTIMP.eq.5) then
C FIXME: Write me
        call abend2('Not implemented yet.')
C       call reo_info_molcas()
        IRILP = 1
      else if (INTIMP.eq.6) then
        call reo_info_dirac(NDIM_TOT,NIRR_PN,SYMDIM,IBOSYM,
     &                      SYDI,DOUGRP,WORK(KISREO),WORK(KINVREO),
     &                      INVSYM,IDGTYP,ISPINFREE,IPRNT)
*
        IRILP = 2    ! Always in LUCIAREL. General r/i vectors.
      else if (INTIMP.eq.7) then
        call reo_info_dirac(NDIM_TOT,NIRR_PN,SYMDIM,IBOSYM,
     &                      SYDI,DOUGRP,WORK(KISREO),WORK(KINVREO),
     &                      INVSYM,IDGTYP,ISPINFREE,IPRNT)
        if (ISPINFREE.eq.1) then
          IRILP = 1
        else
          IRILP = 2
        end if
        IRILP = 2    ! enforced !
      end if
*
* Open ASCII unit for checking integral sorting information
      if (NTEST.ge.1) then
        if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
          LUINGR = LUSC94
          open(unit=LUINGR,file='INTEGRALS',
     &         status='UNKNOWN',form='FORMATTED')
          rewind(LUINGR)
        end if
      end if
*
* Loop over excitation types
      write(6,*) 'Number of types:',NSPOBEX_TP
*
* Loop over real/imaginary part
      ICOUNT = 1
      do IRI=1,IRILP,1
        if (NTEST.ge.1) then
          if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
            write(6,*)
            if (IRI.eq.1) write(6,*) 'REAL LOOP'
            if (IRI.eq.1) write(LUINGR,'(A)') 'Real  integrals.'
            if (IRI.eq.2) write(6,*) 'IMAG LOOP'
            if (IRI.eq.2) write(LUINGR,'(A)') 'Imag. integrals.'
          end if
        end if
        ITPSV = 0
        if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
          if (NTEST.ge.1) write(LUINGR,'(A)') 'INTEGRAL FILE'
        end if
        do IXTP=1,NSPOBEX_TP,1
*
* Check for correct pointers to type blocks
          if (IRI.eq.1.and.ICOUNT.ne.IBSPOBEX_TP(IXTP)) then
            write(6,*) 'Pointers to types do not match numbers'
            write(6,*) 'of integrals generated.'
            write(6,*) 'IHINDX array: ',(IHINDX(I,IXTP),I=1,4,1)
            write(6,*) 'LABEL array: ',(LABEL(I,IXTP),I=1,4,1)
            write(6,*) 'int.no.,offset: ',ICOUNT,IBSPOBEX_TP(IXTP)
            write(6,*) 'Current type: ',IXTP
            call abend2('Quitting.')
          end if
          if (NTEST.ge.2) then
            write(6,*) '======================'
            write(6,*) 'Processing type ',IXTP
            write(6,*) '======================'
            write(6,*) 'ISOBEX partition:'
            call iwrtma(ISOBEX(1,1,IXTP),NGAS,4,NGAS,4)
          end if
*
* write statement about next integral type (to file)
          if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
            if (NTEST.ge.1) 
     &          call integral_type(ISOBEX(1,1,IXTP),NGAS,ITPSV,LUINGR)
            if (NTEST.ge.1) write(LUINGR,'(A,I5)') 'Type ',IXTP
          end if
*
*. Change from Occupation to group notation of strings
          call occ_to_grp(ISOBEX(1,1,IXTP),IGRP_CA,1)
          call occ_to_grp(ISOBEX(1,2,IXTP),IGRP_CB,1)
          call occ_to_grp(ISOBEX(1,3,IXTP),IGRP_AA,1)
          call occ_to_grp(ISOBEX(1,4,IXTP),IGRP_AB,1)
          if (NTEST.ge.5.and.IRI.eq.1) then
            write(6,*) 'Group information:'
            write(6,'(A,4I6)') 'IGRP_CA : ',(IGRP_CA(I),I=1,NGAS,1)
            write(6,'(A,4I6)') 'IGRP_CB : ',(IGRP_CB(I),I=1,NGAS,1)
            write(6,'(A,4I6)') 'IGRP_AA : ',(IGRP_AA(I),I=1,NGAS,1)
            write(6,'(A,4I6)') 'IGRP_AB : ',(IGRP_AB(I),I=1,NGAS,1)
          end if
*
* Corresponding number of electrons
          NEL_CA = ielsum(ISOBEX(1,1,IXTP),NGAS)
          NEL_CB = ielsum(ISOBEX(1,2,IXTP),NGAS)
          NEL_AA = ielsum(ISOBEX(1,3,IXTP),NGAS)
          NEL_AB = ielsum(ISOBEX(1,4,IXTP),NGAS)
          if (NTEST.ge.2.and.IRI.eq.1) then
            write(6,*) 'NEL_CA = ',NEL_CA
            write(6,*) 'NEL_CB = ',NEL_CB
            write(6,*) 'NEL_AA = ',NEL_AA
            write(6,*) 'NEL_AB = ',NEL_AB
          end if
*
* Symmetry loops
* (Annihilator symm. is assumed in terms of corresponding creators!)
*  Symmetry restrictions according to number of electrons
          NEL_C = NEL_CA + NEL_CB
          NEL_A = NEL_AA + NEL_AB
          call evenodd(IEONELC,NEL_C)
          if (IEONELC.eq.1) then
            ISMCSTA = NIRR_PN + 1
            ISMCEND = NIRR_DG
          else if (IEONELC.eq.2) then
            ISMCSTA = 1
            ISMCEND = NIRR_PN
          end if
*
          call evenodd(IEONELCA,NEL_CA)
          if (IEONELCA.eq.1) then
            ISMCASTA = NIRR_PN + 1
            ISMCAEND = NIRR_DG
          else if (IEONELCA.eq.2) then
            ISMCASTA = 1
            ISMCAEND = NIRR_PN
          end if
*
          if (NTEST.ge.5.and.IRI.eq.1) then
            write(6,*) 'C  symmetry loop: ',ISMCSTA,ISMCEND
            write(6,*) 'CA symmetry loop: ',ISMCASTA,ISMCAEND
          end if
*
* Symmetry loops
          do ISM_C=ISMCSTA,ISMCEND,1
            ISM_A = IDBGMULT(ISYM_T,ISM_C)
            do ISM_CA=ISMCASTA,ISMCAEND,1
              ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
              do ISM_AA=1,NIRR_DG,1
                ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))
                if (NTEST.ge.10.and.IRI.eq.1) then
                  write(6,*)
                  write(6,*)
                  write(6,*) 'Printing symmetries:'
                  write(6,*) 'Creator symm.      :',ISM_C
                  write(6,*) 'Annihilator symm.  :',ISM_A
                  write(6,*) 'Creator alpha symm.:',ISM_CA
                  write(6,*) 'Creator beta symm. :',ISM_CB
                  write(6,*) 'Annihi. alpha symm.:',ISM_AA
                  write(6,*) 'Annihi. beta symm. :',ISM_AB
                end if
*
*. Generate strings of these symmetries
*. obtain strings
                IUB = 1
                call getstr2_totsm_spgp(IUB,IGRP_CA,NGAS,ISM_CA,NEL_CA,
     &                            NSTR_CA,ISTR_CA,NORBT,0,IDUM,IDUM)
                IUB = 2
                call getstr2_totsm_spgp(IUB,IGRP_CB,NGAS,ISM_CB,NEL_CB,
     &                            NSTR_CB,ISTR_CB,NORBT,0,IDUM,IDUM)
                IUB = 1
                call getstr2_totsm_spgp(IUB,IGRP_AA,NGAS,ISM_AA,NEL_AA,
     &                            NSTR_AA,ISTR_AA,NORBT,0,IDUM,IDUM)
                IUB = 2
                call getstr2_totsm_spgp(IUB,IGRP_AB,NGAS,ISM_AB,NEL_AB,
     &                            NSTR_AB,ISTR_AB,NORBT,0,IDUM,IDUM)
*
* Loops over strings for the given symmetries and types
* Loop over T elements as matrix T(I_CA,I_CB,I_AA,I_AB)
*  Determine possible connections
                if ((NSTR_CA.ge.1.or.NSTR_CB.ge.1).and.
     &              (NSTR_AA.ge.1.or.NSTR_AB.ge.1)) then
*
                  do I_AB=1,NSTR_AB,1
                    do I_AA=1,NSTR_AA,1
                      do I_CB=1,NSTR_CB,1
                        do I_CA=1,NSTR_CA,1
* Notes: 
*  1) Standard notation in second quantization for 2-el. ints is:
*     a+I a+K_ aL aJ_ (IJ|KL) = - a+I a+K_ aJ_ aL (IJ|KL)
*     The code proceeds according to the right hand side. The
*     signs are located in WORK(KSIGNNHX) for each excitation
*     class. So first indices 2 and 3 are interchanged, then
*     3 and 4 (implicit in the code).
*  2) Integral classes with 2 alpha or 2 beta creators (or both)
*     are treated differently: (due to exchange character)
*     a+I a+K aL aJ [(IJ|KL) - (IL|KJ)]
*  3) In the full spin-dependent formalism, extra classes of
*     integrals are required. 
* 
*. Indices of STRING I_CA are stored in IOCC_CA and first element is 
*  IOCC_CA(1+NEL_CA*(I_CA-1)) 
*
                          if (NEL_C.eq.1) then
                            IXCH = 0
*
                            if (NEL_CA.eq.1.and.NEL_CB.eq.0) then
                              INDEX(1,ICOUNT) 
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(2,ICOUNT) = 0
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.1) then
                              INDEX(1,ICOUNT) = 0
                              INDEX(2,ICOUNT) 
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else
                              write(6,*) 
     &                              'Impossible number of creators.'
                              write(6,*) 
     &                              'NEL_CA,NEL_CB = ',NEL_CA,NEL_CB
                              call abend2('Quitting.')
                            end if
*
                            if (NEL_AA.eq.1.and.NEL_AB.eq.0) then
                              INDEX(3,ICOUNT) 
     &                              = ISTR_AA(1+NEL_AA*(I_AA-1))
                              INDEX(4,ICOUNT) = 0
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.1) then
                              INDEX(3,ICOUNT) = 0
                              INDEX(4,ICOUNT) 
     &                              = ISTR_AB(1+NEL_AB*(I_AB-1))
                            else
                              write(6,*) 
     &                              'Impossible number of annihil.'
                              write(6,*) 
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
*
                          else if (NEL_C.eq.2) then
*
                            if (NEL_CA.eq.2.and.NEL_CB.eq.0) then
                              IXCH = 0
                              if (NEL_AA.eq.2) then
*                               Exchange integrals
                                IXCH = 1
                                do IEL=1,NEL_CA,1
                                  INDEX(IEL,ICOUNT) = 
     &                                ISTR_CA(IEL+NEL_CA*(I_CA-1))
                                end do
                              else if (NEL_AB.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(1,ICOUNT)
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                                INDEX(2,ICOUNT)
     &                              = ISTR_CA(2+NEL_CA*(I_CA-1))
                              else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
*                               Spin-orbit integral
                                IXCH = 0
                                do IEL=1,NEL_CA,1
                                  INDEX(IEL,ICOUNT) =
     &                                ISTR_CA(IEL+NEL_CA*(I_CA-1))
                                end do
                              end if
                            else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
                              IXCH = 0
                              INDEX(1,ICOUNT) 
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(3,ICOUNT) 
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.2) then
                              IXCH = 0
                              if (NEL_AB.eq.2) then
*                               Exchange integrals
                                IXCH = 1
                                do IEL=1,NEL_CB,1
                                  INDEX(IEL,ICOUNT) = 
     &                                ISTR_CB(IEL+NEL_CB*(I_CB-1))
                                end do
                              else if (NEL_AA.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(1,ICOUNT)
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                                INDEX(2,ICOUNT)
     &                              = ISTR_CB(2+NEL_CB*(I_CB-1))
                              else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
*                               Spin-orbit integral
                                INDEX(1,ICOUNT)
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                                INDEX(3,ICOUNT)
     &                              = ISTR_CB(2+NEL_CB*(I_CB-1))
                              end if
                            else
                              write(6,*) 
     &                              'Impossible number of creators.'
                              write(6,*) 
     &                              'NEL_CA,NEL_CB = ',NEL_CA,NEL_CB
                              call abend2('Quitting.')
                            end if
*
                            if (NEL_AA.eq.2.and.NEL_AB.eq.0) then
                              if (NEL_CA.eq.2) then
                                IXCH = 1
                                do IEL=1,NEL_AA,1
                                 INDEX(IEL+2,ICOUNT) = 
     &                                ISTR_AA(IEL+NEL_AA*(I_AA-1))
                                end do
                              else if (NEL_CB.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(3,ICOUNT)
     &                              = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(4,ICOUNT)
     &                              = ISTR_AA(2+NEL_AA*(I_AA-1))
                              else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spin-orbit integral
                                do IEL=1,NEL_AA,1
                                 INDEX(IEL*2,ICOUNT) =
     &                                ISTR_AA(IEL+NEL_AA*(I_AA-1))
                                end do
                              end if
                            else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
                              IXCH = 0
                              if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spinfree integral
                                INDEX(2,ICOUNT) 
     &                                = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(4,ICOUNT) 
     &                                = ISTR_AB(1+NEL_AB*(I_AB-1))
                              else
*                               Spin-orbit integral
                                INDEX(4,ICOUNT)
     &                                = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(2,ICOUNT)
     &                                = ISTR_AB(1+NEL_AB*(I_AB-1))
                              end if
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.2) then
                              if (NEL_CB.eq.2) then
                                IXCH = 1
                                do IEL=1,NEL_AB,1
                                 INDEX(IEL+2,ICOUNT) = 
     &                                ISTR_AB(IEL+NEL_AB*(I_AB-1))
                                end do
                              else if (NEL_CA.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                INDEX(3,ICOUNT)
     &                              = ISTR_AB(1+NEL_AB*(I_AB-1))
                                INDEX(4,ICOUNT)
     &                              = ISTR_AB(2+NEL_AB*(I_AB-1))
                              else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spin-orbit integral
                                do IEL=1,NEL_AB,1
                                 INDEX(IEL*2,ICOUNT) =
     &                                ISTR_AB(IEL+NEL_AB*(I_AB-1))
                                end do
                              end if
                            else
                              write(6,*) 
     &                              'Impossible number of annihil.'
                              write(6,*) 
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
*
                          else
                            write(6,*) 'Impossible number of creators.'
*  Apply only to 2-electron integrals:
                            write(6,*) 'NEL_C = ',NEL_C
                            call abend2('Quitting.')
                          end if
*
*  Modify if integrals are over relativistic spinors:
                          if ((INTIMP.eq.6.or.INTIMP.eq.7).and.
     &                         ISPINFREE.eq.0) then
                            if (INDEX(1,ICOUNT).ne.0.and.
     &                          INDEX(2,ICOUNT).ne.0) then
                              do IRUN = 1,4,1
                                  INDEX(IRUN,ICOUNT) = 
     &                                INDEX(IRUN,ICOUNT) * 
     &                                LABEL(IRUN,IXTP)
                              end do
*
* Correct internal ordering of two-bar exchange type integrals:
* (This might be implemented in some more elegant way later on.)
                              if (INDEX(1,ICOUNT).gt.0.and.
     &                            INDEX(2,ICOUNT).lt.0.and.
     &                            INDEX(3,ICOUNT).lt.0.and.
     &                            INDEX(4,ICOUNT).gt.0) then
                                J_SV = INDEX(2,ICOUNT)
                                L_SV = INDEX(4,ICOUNT)
                                INDEX(2,ICOUNT) = -L_SV
                                INDEX(4,ICOUNT) = -J_SV
                              end if
                            end if
                          end if
*
* Index reordering and second integral in exchange case:
                          YXINT = 0.D0
                          if (IXCH.eq.1) then
* Index sequence reads:
*  (IJ|KL) (which is not the operator sequence for the desired int!)
* Reorder for Coulomb integral:
*  (IJ|KL) --> (IK|LJ)
* Now create exchange integral:
*  (IJ|KL) --> (IJ|LK)
* (Overall, indices L and J are interchanged and a "-" introduced,
*  due to I+ K+ L J = - I+ K+ J L
* Coulomb:
                            I_SV = INDEX(1,ICOUNT)
                            L_SV = INDEX(2,ICOUNT)
                            K_SV = INDEX(3,ICOUNT)
                            J_SV = INDEX(4,ICOUNT)
                            INDEX(1,ICOUNT) = I_SV
                            INDEX(2,ICOUNT) = K_SV
                            INDEX(3,ICOUNT) = L_SV
                            INDEX(4,ICOUNT) = J_SV
* Exchange:
                            I_ADD(1) = I_SV
                            I_ADD(2) = J_SV
                            I_ADD(3) = L_SV
                            I_ADD(4) = K_SV
                            if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ')
     &                      call getint_sigden(YXINT,I_ADD,IRI,
     &                                    NDIM_TOT,NASHT,F1,F2,IPRNT)
                          end if
*
                          if (CIRUN.ne.'DENS1 '.and.
     &                        CIRUN.ne.'DENS2 ') then
                            call getint_sigden(YINT,INDEX(1,ICOUNT),
     &                                IRI,NDIM_TOT,NASHT,F1,F2,IPRNT)
                            T(ICOUNT) = YINT - YXINT
* FIXME !! Uncommented sign here!!
                            if (IXCH.eq.1) 
     &                          T(ICOUNT) = - T(ICOUNT)
* FIXME !! Uncommented sign here!!
*
*                           write integral to file (for debugging):
                            if (NTEST.ge.10 ) 
     &                          write(6,*) 'INTEGRAL no. ',ICOUNT
                            if (NTEST.ge.1) then
                              if (IXCH.eq.0) then
                                write(LUINGR,'(I4,4X,4I3,1F18.10)') 
     &                                ICOUNT,INDEX(1,ICOUNT),
     &                                       INDEX(2,ICOUNT),
     &                                       INDEX(3,ICOUNT),
     &                                       INDEX(4,ICOUNT),T(ICOUNT)
                              else
                          write(LUINGR,'(I4,4X,4I3,X,A,X,4I3,1F18.10)')
     &                                ICOUNT,INDEX(1,ICOUNT),
     &                                       INDEX(2,ICOUNT),
     &                                       INDEX(3,ICOUNT),
     &                                       INDEX(4,ICOUNT),'-',
     &                             I_ADD(1),I_ADD(2),I_ADD(3),I_ADD(4),
     &                                       T(ICOUNT)
                              end if
                            end if
                          end if
                          ICOUNT = ICOUNT + 1
*
                        end do
*                       ^ End loop over alpha creation strings
                      end do
*                     ^ End loop over beta creation strings
                    end do
*                   ^ End loop over alpha annihilation strings
                  end do
*                 ^ End loop over beta annihilation strings
                end if
*               ^ End if connection condition
              end do
*             ^ End loop over alpha annihilator symmetry
            end do
*           ^ End loop over alpha creator symmetry
          end do
*         ^ End loop over creation part symmetry
        end do
*       ^ End loop over excitation types
        print*,'ICOUNT=',ICOUNT
      end do
*     ^ End loop over real/imaginary
      if (NTEST.ge.1) close(LUINGR)
*
      return
      end
