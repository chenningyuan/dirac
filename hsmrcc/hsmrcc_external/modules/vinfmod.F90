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
!    Open HSMRCC file VINF carrying important data of the run 
!
!   called from:
! ----------------  
!                    InitHSMRBase (source/inihsbase.F90)
!
!------------------------------------------------------------------------------
   implicit none
   integer :: ios
   logical :: existVINF

   if (get_iprint() >=7) then
     write(UOUT,'(/,2x,a,a80)')   & 
     'OpenVINF:  going to open the UNFORMATTED file FileVINF=',FileVINF
     write(UOUT,'(2x,a,i3,a,a)') 'where lu UVINF=',UVINF,'and file name, FileVINF=',FileVINF
     write(UOUT,'(2x,a,l1)') '...and OPTION(GENINT)=',OPTION(GENINT)
     write(UOUT,'(2x,a,l1)') '...and OPTION(DIRAC_MODE)=',OPTION(DIRAC_MODE)
     if (OPTION(GENINT).or.OPTION(DIRAC_MODE)) then
     !if (OPTION(GENINT).and..not.OPTION(DIRAC_MODE) )  then
       write(UOUT,'(6x,a)') 'going to open as UNKNOWN (new)'
     else
       write(UOUT,'(6x,a)') 'going to open as OLD - previously prepared file'
     endif
   endif


  ! if (OPTION(GENINT)) then
   if (OPTION(GENINT).or.OPTION(DIRAC_MODE) ) then
      open(UVINF,file=FileVINF,access='SEQUENTIAL',form='UNFORMATTED',status='NEW',iostat=ios)
      if (ios.ne.0) then
        call quit('OpenVINF: could not open NEW UVINF...')
      endif
   else
      inquire(file=FileVINF,exist=existVINF)

      if (existVINF) then
      ! ... restarting HSCC 
        open(UVINF,file=FileVINF,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')
      else
        write(UOUT,'(a,a)') FileVINF," file NOT found ! Error exit !"
        call quit('OpenVINF: FileVINF NOT found !')
      endif

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
! ... plus from others (source/endhsmr.F90)
!
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
   ! ...  05 [ SymInfHSCC ] : point_group, naxis, norder, numirrep, irrep_product,               ...
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

   ! ... symmetry info ...
   write(UVINF) point_group, naxis, norder, numirrep, irrep_product, &
                totally_symmetric_irrep, irrep_label

   write(UVINF) CoreShiftVal, VirtShiftVal

   write(UVINF) NAO, NSO, NFRZV, TotEn, NucEn, TrfTol, IUNIQ, LABELSIZE, NINTMX

   write(UVINF) SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb, &
                OrbType, OrbSpace, OrbEnergy, OrbCoeff, HCoreMatrix

   write(UVINF) NumRec, RecLen

 ! ... positions, lenghts,  number of integrals ...
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

   ! ... print out data ...
   call printVINFstuff

   return

end subroutine UpdateVINF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine printVINFstuff
!------------------------------------------------------------------------------------
!
!               ... print out data from update of the VINF file
!
!------------------------------------------------------------------------------------

   if (get_iprint() >= 10) then

     call Header('Output from UpdateVINF',-1)

     write(UOUT,'(/,2x,a)') 'UpdateVINF - data set saved into <FileVINF>'

     call print_Options

     write(UOUT,*) 'MAXITER, CONVTOL:',MAXITER, CONVTOL
     write(UOUT,*) 'CUTOFF=',CUTOFF
     write(UOUT,*) 
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

 !    write(UOUT,'(4(a,i5))') 'RecLen(V2H2P)=',RecLen(V2H2P),' NumRec(V2H2P)=', NumRec(V2H2P), & 
 !                ' TotNZInt(V2H2P)=',TotNZInt(V2H2P),' NZIntRec_V2H2P=',NZIntRec_V2H2P

  !   write(UOUT,'(4(a,i5))') 'RecLen(VPPPP)=',RecLen(VPPPP),' NumRec(VPPPP)=',NumRec(VPPPP), & 
  !                ' TotNZInt(VPPPP)=',TotNZInt(VPPPP),' NZIntRec_VPPPP=',NZIntRec_VPPPP

     call print_int_file_attributes
    ! ...
     write(UOUT,*) 
     write(UOUT,'(3(a,i5))') 'StartRec(FHH)=',StartRec(FHH),' RecLen(FHH)=',RecLen(FHH), & 
           ' TotNZInt(FHH)=',TotNZInt(FHH)

     write(UOUT,'(2x,a,i7)') 'RecLen(FHP)=',RecLen(FHP)
     write(UOUT,'(2x,a,i7)') 'RecLen(FPP)=',RecLen(FPP)
     write(UOUT,*)
 

 !    write(UOUT,'(4(a,i7))')  &  
 !            ' StartRec(VPHPP)=',StartRec(VPHPP),' EndngRec(VPHPP)=',EndngRec(VPHPP), &
 !            ' TotNZInt(VPHPP)=',TotNZInt(VPHPP),' NZIntRec_VPHPP =',NZIntRec_VPHPP

 !    write(UOUT,'(4(a,i7))')  &  
 !            ' StartRec(V2HHP)=',StartRec(V2HHP),' EndngRec(V2HHP)=',EndngRec(V2HHP), &
 !            ' TotNZInt(V2HHP)=',TotNZInt(V2HHP),' NZIntRec_V2HHP =',NZIntRec_V2HHP

!     write(UOUT,'(4(a,i7))')  &  
 !            ' StartRec(VHHHP)=',StartRec(VHHHP),' EndngRec(VHHHP)=',EndngRec(VHHHP), &
 !            ' TotNZInt(VHHHP)=',TotNZInt(VHHHP),' NZIntRec_VHHHP =',NZIntRec_VHHHP


     write(UOUT,'(/,3(a,i7))') 'StartRec(FHH)=',StartRec(FHH), & 
     ' EndngRec(FHH)=',EndngRec(FHH),' TotNZInt(FHH)=',TotNZInt(FHH)

     write(UOUT,'(3(a,i7))') 'StartRec(FHP)=',StartRec(FHP),  & 
     ' EndngRec(FHP)=',EndngRec(FHP),' TotNZInt(FHP)=',TotNZInt(FHP)

     write(UOUT,'(3(a,i7))') 'StartRec(FPP)=',StartRec(FPP), & 
     ' EndngRec(FPP)=',EndngRec(FPP),' TotNZInt(FPP)=',TotNZInt(FPP)

     write(UOUT,*) 

     call flushout

 endif

end subroutine printVINFstuff


end module VINFMod
