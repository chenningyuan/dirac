* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*  1-particle: molfdir full spinor format
*  2-particle: molfdir full spinor format
*
*   November 12, 2001, Timo Fleig
*
      subroutine trnsfdens(F1,F2,DENS,INDEX,LABEL,IHINDX,SIGN,IOP_REO,
     &                     SIGN_OPREO,ISF,NASHT,IPRNT)
*
*  Input:  DENS
*  Output: F1,F2
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
      dimension F1(*),F2(*)
      dimension DENS(NINT_CC),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,NINT_CC),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP)
*
      NTESTL = 40
      NTEST = max(IPRNT,NTESTL)
      NTEST = 10  ! forced
*
      NINT_TOT = N1ELINT + N2ELINT
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHT2 = NASHT * 2
*
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) '=========================================='
        write(6,*) ' Transferring densities to MOLFDIR format.'
        write(6,*) '=========================================='
      end if
*
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
*
        if (NTEST.ge.4) then
          write(6,*)
          if (IRI.eq.1) write(6,*) ' Real part     : '
          if (IRI.eq.2) write(6,*) ' Imaginary part: '
          write(6,*) ' -------------------'
        end if
*
        do IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(6,*)
            write(6,'(A,I4)') ' Excitation type ... ',IXTP
            write(6,*)        ' -------------------------'
            write(6,*)
            write(6,*) '  LABEL array (for integral indices!):  '
            write(6,*) '   Internally, LABEL is reord. for operators.'
            call iwrtma(LABEL(1,IXTP),1,4,1,4)
            write(6,*) '   Input LUCIAREL operator order, output dens.'
            call iwrtma(IOP_REO(1,IXTP),1,4,1,4)
          end if
*
          NINT_TP = NELM_CC(WORK(KLIBSOBEX),IXTP,NSPOBEX_TP,NINT_TOT)
          do INT = 1,NINT_TP,1
*
* Determine type of density element
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
*
            if (NCT.eq.2) then
              I = INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHT
              J = INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHT
*
              IOUT =  (IRI-1)*NASHT2*NASHT2 +
     &                 (J-1) *NASHT2 +
     &                  I
              F1(IOUT) = DENS(IRUN+IRIOFF)
              if (NTEST.ge.10) then
                write(6,'(2I3,1F18.10,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
              end if
*
            else if (NCT.eq.4) then
              if (SIGN_OPREO(IXTP).eq.ONE.or.
     &            SIGN_OPREO(IXTP).eq.ONEM) then
*
*
*  INDEX array refers to integral indexing.
*  The corresponding second quantization operators then result
*  from the following resorting to CACA ordering:
*  CA AA CB AB
*  IHINDX relates integral and operator ordering of a particular
*  excitation class. LABEL has been defined for integrals (see
*  get_hx_rela.F), so it is reordered for operators here.
*
                IND_R(IOP_REO(1,IXTP)) = INDEX_L(1)  
     &                   - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHT
                if (NTEST.ge.50) then
                  write(6,*)
                  write(6,*) ' I '
                  write(6,*) ' INDEX_L(1) ',INDEX_L(1)
                  write(6,*) ' IHINDX(1,IXTP) ' ,IHINDX(1,IXTP)
                  write(6,*) ' LABEL(IHINDX(1,IXTP),IXTP) ',
     &                         LABEL(IHINDX(1,IXTP),IXTP)
                end if
*
                IND_R(IOP_REO(2,IXTP)) = INDEX_L(2)  
     &                   - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHT
                if (NTEST.ge.50) then
                  write(6,*)
                  write(6,*) ' L '
                  write(6,*) ' INDEX_L(2) ',INDEX_L(2)
                  write(6,*) ' IHINDX(2,IXTP) ' ,IHINDX(2,IXTP)
                  write(6,*) ' LABEL(IHINDX(2,IXTP),IXTP) ',
     &                         LABEL(IHINDX(2,IXTP),IXTP)
                end if
*
                IND_R(IOP_REO(3,IXTP)) = INDEX_L(3)  
     &                   - (LABEL(IHINDX(3,IXTP),IXTP) - 1)/2 * NASHT
                if (NTEST.ge.50) then
                  write(6,*)
                  write(6,*) ' K '
                  write(6,*) ' INDEX_L(3) ',INDEX_L(3)
                  write(6,*) ' IHINDX(3,IXTP) ' ,IHINDX(3,IXTP)
                  write(6,*) ' LABEL(IHINDX(3,IXTP),IXTP) ',
     &                         LABEL(IHINDX(3,IXTP),IXTP)
                end if
*
                IND_R(IOP_REO(4,IXTP)) = INDEX_L(4)  
     &                   - (LABEL(IHINDX(4,IXTP),IXTP) - 1)/2 * NASHT
                if (NTEST.ge.50) then
                  write(6,*)
                  write(6,*) ' J '
                  write(6,*) ' INDEX_L(4) ',INDEX_L(4)
                  write(6,*) ' IHINDX(4,IXTP) ' ,IHINDX(4,IXTP)
                  write(6,*) ' LABEL(IHINDX(4,IXTP),IXTP) ',
     &                         LABEL(IHINDX(4,IXTP),IXTP)
                end if
*
                IOUT = (IRI-1)*NASHT2*NASHT2*NASHT2*NASHT2 +
     &                  (IND_R(4)-1) *NASHT2*NASHT2*NASHT2 +
     &                  (IND_R(3)-1) *NASHT2*NASHT2 +
     &                  (IND_R(2)-1) *NASHT2 +
     &                   IND_R(1)
                F2(IOUT) = DENS(IRUN+IRIOFF) * SIGN(IXTP) 
     &                                       * SIGN_OPREO(IXTP)
                if (NTEST.ge.10) then
                  write(6,*)
                  write(6,'(4I3,1F18.10,A,I4)') 
     &                      IND_R(1),IND_R(2),IND_R(3),IND_R(4),
     &                      F2(IOUT),' to pos.',IOUT
                end if
*
* ... and copy to the entry for allowed index permutation
* (ab|cd) = -(ad|cb)
*
                I_R = IND_R(1)
                L_R = IND_R(2)
                K_R = IND_R(3)
                J_R = IND_R(4)
                SIGN_CP = -1 * ONE
                IOUT_CP = (IRI-1) *NASHT2*NASHT2*NASHT2*NASHT2 +
     &                    (L_R-1) *NASHT2*NASHT2*NASHT2 +
     &                    (K_R-1) *NASHT2*NASHT2 +
     &                    (J_R-1) *NASHT2 +
     &                     I_R
                F2(IOUT_CP) = F2(IOUT) * SIGN_CP
                if (NTEST.ge.10) then
                  write(6,'(4I3,1F18.10,A,I4)') 
     &                      I_R,J_R,K_R,L_R,F2(IOUT_CP),
     &                      ' to pos.',IOUT_CP
                end if
*
              end if
*             ^ Non-redundant density type (see get_hx_rela)
            else 
              write(6,*) 'No QFT implemented YET in LUCIAREL.'
              write(6,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting.')
            end if
            IRUN = IRUN + 1
          end do
*         ^ Loop over integrals in excitation type
        end do
*       ^ Loop over excitation types
      end do
*     ^ Loop over Real/Imaginary
*
      return
      end
