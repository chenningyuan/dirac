      SUBROUTINE WRTHD(LUOUT)
*
* Introduce LUCIAREL
*
      WRITE(6,'(1H ,8X,A)')
     &'****************************************************************'
      WRITE(6,'(1H ,8X,A)')
     &'*                                                              *'
      WRITE(6,'(1H ,8X,A)')
     &'*           Welcome to LUCIAREL, a direct relativistic         *'
      WRITE(6,'(1H ,8X,A)')
     &'*                                CI program                    *'
      WRITE(6,'(1H ,8X,A)')
     &'*                                                              *'
      WRITE(6,'(1H ,8X,A)')
     &'* LUCIA written by Jeppe Olsen, University of Lund, Sweden     *'
      WRITE(6,'(1H ,8X,A)')
     &'* Relativistic CI version by Timo Fleig and Jeppe Olsen        *'
      WRITE(6,'(1H ,8X,A)')
     &'*    - Double group symmetry implemented: JO and TF            *'
      WRITE(6,'(1H ,8X,A)')
     &'*    - Time-reversal symmetry implemented: JO and TF           *'
      WRITE(6,'(1H ,8X,A)')
     &'*    - Integral transformation section (one-el, spin-orbit): TF*'
      WRITE(6,'(1H ,8X,A)')
     &'*    - Integral transformation section (two-el): JO and TF     *'
      WRITE(6,'(1H ,8X,A)')
     &'*                                                              *'
      WRITE(6,'(1H ,8X,A)')
     &'* V 2.0                                                        *'
      WRITE(6,'(1H ,8X,A)')
     &'* Version of November 2000                                     *'
      WRITE(6,'(1H ,8X,A)')
     &'*   1.0       R-GAS-CI                                         *'
      WRITE(6,'(1H ,8X,A)')
     &'*   1.1       MOLCAS version                                   *'
      WRITE(6,'(1H ,8X,A)')
     &'*   1.2       External Fock matrices from MOTRA in use         *'
      WRITE(6,'(1H ,8X,A)')
     &'*   1.3       Density matrices can be obtained                 *'
      WRITE(6,'(1H ,8X,A)')
     &'*             (Implement C_1, C_i, and C_2h)                   *'
      WRITE(6,'(1H ,8X,A)')
     &'*    1.3.1     Integral transformation section generalized     *'
      WRITE(6,'(1H ,8X,A)')
     &'*   1.4       Resolution matrix batching active                *'
      WRITE(6,'(1H ,8X,A)')
     &'*   1.5       Analysis of optimized CI vectors possible        *'
      WRITE(6,'(1H ,8X,A)')
     &'*   1.6       Quaternionic and complex double groups active    *'
      WRITE(6,'(1H ,8X,A)')
     &'*   1.7       MOLCAS 5.0 version (input section changed)       *'
      WRITE(6,'(1H ,8X,A)')
     &'*    1.7.1     Abend calls for parallel runs included          *'
      WRITE(6,'(1H ,8X,A)')
     &'*    1.7.2     SEWARD-AMFI calls implemented for MOLCAS4.3 - 5.*'
      WRITE(6,'(1H ,8X,A)')
     &'*   1.8       DIRAC 4-component version                        *'
      WRITE(6,'(1H ,8X,A)')
     &'*    1.8.1     Spinfree formalism (Dyall, 1993) implemented    *'
      WRITE(6,'(1H ,8X,A)')
     &'*    1.8.2     Full spin-dependent formalism part. implemented *'
      WRITE(6,'(1H ,8X,A)')
     &'*                    ------------                              *'
      WRITE(6,'(1H ,8X,A)')
     &'*   2.0       Version for relativistic MCSCF                   *'
      WRITE(6,'(1H ,8X,A)')
     &'*              General sigma and density (from CC) implemented *'
      WRITE(6,'(1H ,8X,A)')
     &'****************************************************************'
*
      WRITE(6,*)
      WRITE(6,'(1H ,8X,A)')
     &' In case of trouble please contact : '
      WRITE(6,*)
      WRITE(6,'(1H ,12X,2A)')
     &' Timo Fleig, Theoretical Chemistry, Faculty of Sciences '
      WRITE(6,'(1H ,12X,2A)')
     &' Free University of Amsterdam                   '
      WRITE(6,'(1H ,12X,2A)')
     &' De Boelelaan 1083, NL-1081 HV Amsterdam, NL    '
      WRITE(6,'(1H ,12X,2A)')
     &' Telephone : +31 20 44 47627 '
      WRITE(6,'(1H ,12X,2A)')
     &' E-mail : fleig@chem.vu.nl '
      WRITE(6,*)
     &' or : '
      WRITE(6,*)
      WRITE(6,'(1H ,12X,2A)')
     &' Jeppe Olsen, Teoretisk Kemi, Aarhus Universitet '
      WRITE(6,'(1H ,12X,2A)')
     &' Langelandsgade 140, DK-8000 Aarhus C '
      WRITE(6,'(1H ,12X,2A)')
     &' E-mail : jeppe@signe.teokem.lu.se  or'
      WRITE(6,'(1H ,12X,2A)')
     &'          jeppe@aausgi.uni-c.dk'
      WRITE(6,*)
      WRITE(6,'(1H ,8X,A)')
      write(6,*)
      write(6,*)
      write(6,*) '-----------------------------------'
      write(6,*) '|  Stability on systems:          |'
      write(6,*) '-----------------------------------'
      write(6,*)
      write(6,*) ' SGI (Octane, Origin 2000) '
      write(6,*) ' IBM (RS6000)              '
      write(6,*) ' REDHAT LINUX              '
      write(6,*) ' SuSE Linux                '
      write(6,*)
      write(6,*)
*
      RETURN
      END
