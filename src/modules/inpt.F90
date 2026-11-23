!-----------------------------------------
!  Module with variables previously defined
!  at src/relccsd/inpt.inc include file
!-----------------------------------------
   module inpt_mod
   implicit none
   integer, parameter :: MAXOP = 15
   integer, parameter :: MAXFR = 25
   real*8 :: FF_PROP_STRENGTHS(2,MAXOP), EFREQ(MAXFR)
   integer :: IPRNT, NELEC(16), NELEC_F1(16), NELEC_F2(16), NFROZ(16), MWORD, NFFOPER, &
              NTOLCC, MXDIMCC, MXITCC, NFOPER, NTOLGR, MXDIMGR, MXITGR, NAMEF(MAXOP), &
              NAOPER, NBOPER, NFREQ, NTOLRP, MXDIMRP, MXITRP, NEXC(16), NEOPER, NTOLEX, &
              MXDIMEX, MXITEX, FSSECT(6), MAXIT_SECT(6), NINACT(16), NACT(16), NTOLFS, &
              MXDIMFS, MXITFS, STATE_E_FSCC
   logical :: DEBUG, TIMING, BREIT, DOSORT, DOENER, DOFOPR, DOSOPR, DOFSPC, &
              ADD_FINITE_FIELD, DORESTART, DOMEMCOUNT, USEOE, NORECMP, LOWMEM, DOMP2, &
              DOCCSD, DOCCSDT, NOCCS, NOCCD, DOEOMCCSD, DOMP2G, DOMP2GOLD, DOCCSDG, &
              DOCCSDTG, NOZG, DONATORB, FS_SKIP_SECT(6), FS_REDO_SECT(6), &
              RESTART_UNCONVERGED, RESTART_REDO_SORT, RESTART_REDO_CCSD, &
              RESTART_IGNORE_CHECK, DOEOMH0DD, DOEOMPT2, DOEOMPT2CC
   character(len=10) :: INTERFACE
   character(len=8) :: FF_PROP_NAMES(MAXOP), NAMEA(MAXOP), NAMEB(MAXOP), NAMEE(MAXOP)
   end module inpt_mod
