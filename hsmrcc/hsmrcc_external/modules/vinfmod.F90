module VINFMod

   use Parameters
   use Options
   use IOModule
   use ConvMod
   use Shift
   use GMSMod
   use IntGenMod
   implicit none

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine OpenVINF
!------------------------------------------------------------------------------
!
!    Open GAMESS file  
!
!    called from InitHSMRBase
!
!------------------------------------------------------------------------------
   implicit none

   if (get_iprint() >=7) then
     write(UOUT,'(/,2x,a,a80)')   & 
     'OpenVINF:  going to open the Gamess input file FileVINF=',FileVINF
     write(UOUT,*) ' UVINF=',UVINF,' FileVINF=',FileVINF
     write(UOUT,'(2x,a,l1)') '...and OPTION(GENINT)=',OPTION(GENINT)
   endif

   if (OPTION(GENINT)) then
!     open(UVINF,file=FileVINF,access='SEQUENTIAL',form='UNFORMATTED',status='NEW')
      open(UVINF,file=FileVINF,access='SEQUENTIAL',form='UNFORMATTED',status='UNKNOWN')
   else
      open(UVINF,file=FileVINF,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')
   endif

   return

end subroutine OpenVINF

subroutine CloseVINF

   implicit none

   close(UVINF)
   return

end subroutine CloseVINF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine UpdateVINF
!---------------------------------------------------------------------------------
!
! ... called from HSMRCC (source/hsmrcc.F90) before iterations cycle starts
!
!---------------------------------------------------------------------------------

   implicit none

   ! ...............................................................................................
   ! ...                        Write Integral Information into FileVINF                         ...
   ! ...                        FileVINF contains following info on                              ...
   ! ...                        current run, this can be used for restart                        ...
   ! ...............................................................................................
   ! ...  Rec     Module    : Variables present in the record                                    ...
   ! ...............................................................................................
   ! ...  01 [ Parameters ] : NCMO, NHOLE, NPTCL, NEMO                                           ...
   ! ...  02 [ Options    ] : OPTION                                                             ...
   ! ...  03 [ ConvMod    ] : MAXITER, CONVTOL                                                   ...
   ! ...  04 [ IntGenMod  ] : CUTOFF                                                             ...
   ! ...  05 [ SymInf     ] : point_group, naxis, norder, numirrep, irrep_product,               ...
   ! ...                      totally_symmetric_irrep, irrep_label                               ...
   ! ...  06 [ Shift      ] : CoreShiftVal, VirtShiftVal                                         ...
   ! ...  07 [ GMSMod     ] : NAO, NSO, NFRZV, TotEn, NucEn,                                     ...
   !                          TrfTol, IUNIQ, LABELSIZE, NINTMX                                   ...
   ! ...  08 [ GMSMod     ] : SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb,                         ...
   !                          OrbType, OrbSpace, OrbEnergy, OrbCoeff, HCoreMatrix                ...
   ! ...  09 [ IOModule   ] : NumRec, RecLen                                                     ...
   ! ...  10 [ IOModule   ] : StartRec(V2H2P), EndngRec(V2H2P), TotNZInt(V2H2P), NZIntRec_V2H2P  ...
   ! ...  11 [ IOModule   ] : StartRec(VPPPP), EndngRec(VPPPP), TotNZInt(VPPPP), NZIntRec_VPPPP  ...
   ! ...  12 [ IOModule   ] : StartRec(VHHHH), EndngRec(VHHHH), TotNZInt(VHHHH), NZIntRec_VHHHH  ...
   ! ...  13 [ IOModule   ] : StartRec(VHHPP), EndngRec(VHHPP), TotNZInt(VHHPP), NZIntRec_VHHPP  ...
   ! ...  14 [ IOModule   ] : StartRec(V2HPP), EndngRec(V2HPP), TotNZInt(V2HPP), NZIntRec_V2HPP  ...
   ! ...  15 [ IOModule   ] : StartRec(V2PHH), EndngRec(V2PHH), TotNZInt(V2PHH), NZIntRec_V2PHH  ...
   ! ...  16 [ IOModule   ] : StartRec(VPHHP), EndngRec(VPHHP), TotNZInt(VPHHP), NZIntRec_VPHHP  ...
   ! ...  17 [ IOModule   ] : StartRec(VPH2P), EndngRec(VPH2P), TotNZInt(VPH2P), NZIntRec_VPH2P  ...
   ! ...  18 [ IOModule   ] : StartRec(VPHPP), EndngRec(VPHPP), TotNZInt(VPHPP), NZIntRec_VPHPP  ...
   ! ...  19 [ IOModule   ] : StartRec(V2HHP), EndngRec(V2HHP), TotNZInt(V2HHP), NZIntRec_V2HHP  ...
   ! ...  20 [ IOModule   ] : StartRec(VHHHP), EndngRec(VHHHP), TotNZInt(VHHHP), NZIntRec_VHHHP  ...
   ! ...  21 [ IOModule   ] : StartRec(FHH)  , EndngRec(FHH)  , TotNZInt(FHH)                    ...
   ! ...  22 [ IOModule   ] : StartRec(FHP)  , EndngRec(FHP)  , TotNZInt(FHP)                    ...
   ! ...  23 [ IOModule   ] : StartRec(FPP)  , EndngRec(FPP)  , TotNZInt(FPP)                    ...
   ! ...............................................................................................

   rewind(UVINF)

   write(UVINF) NCMO, NHOLE, NPTCL, NEMO 

   write(UVINF) OPTION

   write(UVINF) MAXITER, CONVTOL

   write(UVINF) CUTOFF

   write(UVINF) point_group, naxis, norder, numirrep, irrep_product, &
                totally_symmetric_irrep, irrep_label
   write(UVINF) CoreShiftVal, VirtShiftVal
   write(UVINF) NAO, NSO, NFRZV, TotEn, NucEn, TrfTol, IUNIQ, LABELSIZE, NINTMX

   write(UVINF) SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb, &
                OrbType, OrbSpace, OrbEnergy, OrbCoeff, HCoreMatrix

   write(UVINF) NumRec, RecLen

   write(UVINF) StartRec(V2H2P), EndngRec(V2H2P), TotNZInt(V2H2P), NZIntRec_V2H2P
   write(UVINF) StartRec(VPPPP), EndngRec(VPPPP), TotNZInt(VPPPP), NZIntRec_VPPPP
   write(UVINF) StartRec(VHHHH), EndngRec(VHHHH), TotNZInt(VHHHH), NZIntRec_VHHHH
   write(UVINF) StartRec(VHHPP), EndngRec(VHHPP), TotNZInt(VHHPP), NZIntRec_VHHPP
   write(UVINF) StartRec(V2HPP), EndngRec(V2HPP), TotNZInt(V2HPP), NZIntRec_V2HPP
   write(UVINF) StartRec(V2PHH), EndngRec(V2PHH), TotNZInt(V2PHH), NZIntRec_V2PHH
   write(UVINF) StartRec(VPHHP), EndngRec(VPHHP), TotNZInt(VPHHP), NZIntRec_VPHHP
   write(UVINF) StartRec(VPH2P), EndngRec(VPH2P), TotNZInt(VPH2P), NZIntRec_VPH2P
   write(UVINF) StartRec(VPHPP), EndngRec(VPHPP), TotNZInt(VPHPP), NZIntRec_VPHPP
   write(UVINF) StartRec(V2HHP), EndngRec(V2HHP), TotNZInt(V2HHP), NZIntRec_V2HHP
   write(UVINF) StartRec(VHHHP), EndngRec(VHHHP), TotNZInt(VHHHP), NZIntRec_VHHHP

   write(UVINF) StartRec(FHH)  , EndngRec(FHH)  , TotNZInt(FHH) 
   write(UVINF) StartRec(FHP)  , EndngRec(FHP)  , TotNZInt(FHP)
   write(UVINF) StartRec(FPP)  , EndngRec(FPP)  , TotNZInt(FPP)

   if (get_iprint() >= 10) then
     call Header('Output from UpdateVINF',-1)
     write(UOUT,'(/,2x,a)') 'UpdateVINF output - data set saved into <FileVINF>'
     write(UOUT,*) 'OPTION:',OPTION
     write(UOUT,*) 'MAXITER, CONVTOL:',MAXITER, CONVTOL
     write(UOUT,*) 'CUTOFF=',CUTOFF
     write(UOUT,*) 'point_group=',point_group,' naxis=',naxis
     write(UOUT,*) 
     write(UOUT,*) 'NAO=',NAO
     write(UOUT,*) 'NSO=',NSO
     write(UOUT,*) 'NFRZV=',NFRZV
     write(UOUT,*) 
     write(UOUT,*) 'TotEn=',TotEn
     write(UOUT,*) 'NucEn=',NucEn
     write(UOUT,*) 'IUNIQ=',IUNIQ
     write(UOUT,*) 'LABELSIZE=',LABELSIZE
     write(UOUT,*) 'NINTMX=',NINTMX
     write(UOUT,*) 
     write(UOUT,*) 'NumRec=',NumRec
     write(UOUT,*) 'RecLen=',RecLen
     write(UOUT,*) 
     write(UOUT,*) 'StartRec(FHH)=',StartRec(FHH)
     write(UOUT,*) 'EndngRec(FHH)=',EndngRec(FHH)
     write(UOUT,*) 'TotNZInt(FHH)=',TotNZInt(FHH)
     write(UOUT,*) 
     write(UOUT,*) 'StartRec(FHP)=',StartRec(FHP)
     write(UOUT,*) 'EndngRec(FHP)=',EndngRec(FHP)
     write(UOUT,*) 'TotNZInt(FHP)=',TotNZInt(FHP)
     write(UOUT,*) 
     write(UOUT,*) 'StartRec(FPP)=',StartRec(FPP)
     write(UOUT,*) 'EndngRec(FPP)=',EndngRec(FPP)
     write(UOUT,*) 'TotNZInt(FPP)=',TotNZInt(FPP)
   endif

   return

end subroutine UpdateVINF

end module VINFMod
