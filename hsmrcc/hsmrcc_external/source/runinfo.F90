subroutine PrintRunInfo
!
!  print basic run info ...
!
   use Parameters
   use Options
   use IOModule
   use IntGenMod
   use VacMod
   use SymInf
   use GMSMod
   use OrbMod
   use Shift
   use ConvMod
   implicit none
   integer :: i, vac, ovac
   character*4 :: SpcStr, TypStr, SpnStr

   write(UOUT,9001)

9001 format(/,                                                                    &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|      PROGRAM FOR EFFECTIVE AND INTERMEDIATE HAMILTONIAN       |',/, &
            '|     HILBERT SPACE MULTI-REFERENCE COUPLED CLUSTER METHODS     |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '|           Shamasundar. K. R. and Sourav Pal                   |',/, &
            '|           Theoretical Chemistry Group                         |',/, &
            '|           Physical Chemistry Division                         |',/, &
            '|           National Chemical Laboratory (NCL)                  |',/, &
            '|           Pashon Road, Pune - 411 008                         |',/, &
            '|           INDIA                                               |',/, &
            '|                                                               |',/, &
            '|                       AND                                     |',/, &
            '|                                                               |',/, &
            '|           Ephraim Eliav and Uzi Kaldor                        |',/, &
            '|           School of Chemistry                                 |',/, &
            '|           Tel Aviv University                                 |',/, &
            '|           Tel Aviv, 69978                                     |',/, &
            '|           ISRAEL                                              |',/, &
            '|                                                               |',/, &
            '|    Note : This version is ported for opteron machines at      |',/, &
            '|           University of Waterloo and Tel Aviv University.     |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                INFORMATION ON CURRENT RUN                     |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'/)

!  print out basic parameters

   write(UOUT,9011) NTMO, NCMO, NVMO, NEMO, NHOLE, NPTCL, NELE, NSTATE, NTSIZE

9011 format('****************',/,&
            'Basic Parameters',/,&
            '****************',/,&
            'Number of spin-orbitals   (NTMO)   = ',i3,/, &
            '                   core   (NCMO)   = ',i3,/, &
            '                 active   (NVMO)   = ',i3,/, &
            '                virtual   (NEMO)   = ',i3,/, &
            '      valence     holes   (NHOLE)  = ',i3,/, &
            '      valence particles   (NPTCL)  = ',i3,/, &
            'Number of electrons       (NELE)   = ',i3,/, &
            'Dimension of model space  (NSTATE) = ',i8,/, &
            'No. of Cluster Amplitudes (NTSIZE) = ',i8,/)

! ...  print out options used for current run ...
   write(UOUT,'(a,i3,/)')  '**************** Program print level:',get_iprint()

   write(UOUT,9021) OnOff(OPTION(GENINT)),   OnOff(OPTION(PRNTINF)),  &
                    OnOff(OPTION(VERBOSE)),  OnOff(OPTION(RESTART)),  &
                    OnOff(OPTION(AVGVAL)),   OnOff(OPTION(REPTIME)),  &
                    OnOff(OPTION(QUAD)),     OnOff(OPTION(CONVTAMP)), &
                    OnOff(OPTION(NOHSMR)),   OnOff(OPTION(NOAPND)),   &
                    OnOff(OPTION(RDGUESS))


9021 format('************',/,&
            'Options used',/,&
            '************',/,&
            'Generate Integrals              [GENINT] = ',a3,/, &
            'Print Information              [PRNTINF] = ',a3,/, &
            'Verbose Output                 [VERBOSE] = ',a3,/, &
            'Restart an unfinished run      [RESTART] = ',a3,/, &
            'Average Valence MO energies     [AVGVAL] = ',a3,/, &
            'Report Subroutine Timings      [REPTIME] = ',a3,/, &
            'Quadratic Approximation           [QUAD] = ',a3,/, &
            'Use of T/HEFF for Convergence [CONVTAMP] = ',a3,/, &
            'Only Intergal Generation        [NOHSMR] = ',a3,/, &
            'Do not append results           [NOAPND] = ',a3,/, &
            'Read T-Amplitude guess         [RDGUESS] = ',a3,/ )  
           

!  print out other information

   write(UOUT,9031) CoreShiftVal, VirtShiftVal, ActvShiftVal, OnOff(OPTION(AVGVAL)), &
                    beeta, CONVTOL, CUTOFF, MAXITER

9031 format('****************',/,  &
            'Misc Information',/,  &
            '****************',/,  &
            'Orbital Energy Shifts for ... ',/, &
            '                   ...     core = ',F5.3,/, &
            '                   ...  virtual = ',F5.3,/, &
            '                   ...   active = ',F5.3,/, &
            'Average active orbital energies = ',a3,/, &
            'Extrapolation Parameter (Beeta) = ',F5.3,/, &
            'Convergence Tolerance           = ',F13.10,/, &
            'Integral Cutoff                 = ',E13.8 ,/, &
            'Maximum No of Iterations        = ',i3,/)


!   ...  information from gamess
! TODO: splitt for GAMESS / DIRAC

   write(UOUT,9041) TotEn, NucEn, TrfTol

9041 format('******************',/,  &
            'Gamess Information',/,  &
            '******************',/,  &
            'Total Electronic Energy         (TotEn) = ',F16.10,/, &
            'Nuclear Energy                  (NucEn) = ',F16.10,/, &
            'Integral Cutoff used in Gamess (TRFTOL) = ',E16.10,/)
           

!  Print out orbital information

   write(UOUT,9051)

   do i = 1, NTMO

      TypStr = "XXXX"
      if (OrbType(i).eq.CORE_ORB) TypStr = "CORE"
      if (OrbType(i).eq.ACTV_ORB) TypStr = "ACTV"
      if (OrbType(i).eq.VIRT_ORB) TypStr = "VIRT"

      SpcStr = "XXXX"
      if (OrbSpace(i).eq.Q_SPACE ) SpcStr = "Q   "
      if (OrbSpace(i).eq.PM_SPACE) SpcStr = "P-m "
      if (OrbSpace(i).eq.PI_SPACE) SpcStr = "P-i "

      SpnStr = "XXXX"
      if (SpinOfOrb(i).eq.ALPHA_SPIN) SpnStr = "ALPH"
      if (SpinOfOrb(i).eq.BEETA_SPIN) SpnStr = "BETA"

      write(UOUT,9052) i, TypStr, SpnStr, IrrepOrb(i),irrep_label(IrrepOrb(i)), &
                          OccGMSRef(i), OccRef(i), SpcStr, &
                          OrbEnergy(i), ShiftEn(i), OrbEnergy(i)-ShiftEn(i)

      if (i.eq.NCMO) write(UOUT,9053)
      if (i.eq.(NCMO+NVMO)) write(UOUT,9053)
      if (i.eq.NTMO) write(UOUT,9053)

   enddo

   write(UOUT,9054)

9051 format('**************************************************************************************',/,&
            '                                Spin-Orbital Information                              ',/,&
            '**************************************************************************************',/,&
            'Orbital Type Spin Irrep IrpLabel OccGMS OccRef Space    Energy     Shift     Shifed-En',/,&
            '**************************************************************************************',/)
           
9052 format(i4,3x,1x,a4,1x,a4,1x,i3,2x,2x,a4,2x,1x,i4,2x,1x,i4,2x,1x,1x,a4,1x,F10.6,1x,F10.6,1x,F10.6)
9053 format(1x)
9054 format('**************************************************************************************')

   ! Print out information on vacua

   write(UOUT,9061)

   do ovac = 1, NSTATE
      vac = RevOrderVac(ovac)
      write(UOUT,9062) ovac,vac,Ms2Vac(vac),IrpVac(vac),&
                       irrep_label(IrpVac(vac)), OccVac(1:NVMO,vac)
   enddo

   write(UOUT,9063)
   write(UOUT,9064)

9061 format('---------------------------------------------------------',/, &
            '           Information on Model Space Vacua              ',/, &
            '---------------------------------------------------------',/, &
            'Vacuum Orig-Order Ms2Vac Irrep IrpLabel Occupation String',/, &
            '---------------------------------------------------------',/)
           
9062 format(i4,2x,1x,2x,i4,4x,1x,i4,2x,1x,i3,2x,1x,2x,a4,2x,1x,20(i1))
9063 format('---------------------------------------------------------')
9064 format(1x)

   write(UOUT,9071)

9071 format('Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out-Out',/)

   return

contains

   character*3 function OnOff(flag)
! ... elegant print-out ...
   implicit none
   logical flag

   if (flag) then
      onoff = "ON "
   else
      onoff = "OFF"
   endif

   return

   end function OnOff

end subroutine PrintRunInfo
