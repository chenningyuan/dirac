!-----------------------------------------
! Module with variables previously defined
! at dcbpop.h include file (COMMON /CBCPOP/, /CBLPOP/, /CBIPOP/, /CBRPOP/)
!-----------------------------------------
!     Local common block for Mulliken Population analysis in dirana.F
!     (set in routine POPINP in dirrdn.F)
!     POPLAB - labels for population analysis
!     VECPOP(2) - orbital strings for population analysis
!     VECPOP_SAVE(2) - for backup copy of user input VECPOP
!     ADD_ALL   : collect large and small components for a given center
!                 into one label (AO-def): for DO_SCFPOP
!     ADDSML    : collect small components for a given center into one
!                 label (AO-def)
!                 collect small components for given center combination
!                 into one label (SO-def)
!     DONETP    : also print net populations
!     LABDEF    : use user defined labels in POPLAB(*)
!     DO_SCFPOP : do atomic gross mull. populations for SCF iterations
!     IPRPOP - print level
!     ILABDF = 1: base label definitions on primitive AO-labels
!            = 2: base label definitions on primitive SO-labels
!     NPOPLAB   : number of labels for population analysis
!     IPOPLAB(1:PNLAB(0)) - pointer from AO/SO-label to POPLAB
!     ICHRG(1:NPOPLAB) - nuclear charge of each entry in POPLAB (for DO_SCFPOP)
!     POP_DTOL : only mulpop contributions > POP_DTOL are printed
   module dcbpop_mod
   use dcblab_mod, only: MAXLAB
   implicit none
   real*8 :: CHRGNUC(MAXLAB), POP_DTOL
   integer :: IPRPOP, ILABDF, NPOPLAB, IPOPLAB(MAXLAB)
   logical :: ADD_ALL, ADDSML, DONETP, LABDEF, DO_SCFPOP
   character(len=12) :: POPLAB(MAXLAB)
   character(len=72) :: VECPOP(2), VECPOP_SAVE(2)
   end module dcbpop_mod
