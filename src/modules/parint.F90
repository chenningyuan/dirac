!-----------------------------------------
! Module with variables previously defined
! at parint.h include file
!-----------------------------------------
   module parint_mod
!
   use mxcent_mod, only: MXCENT
   implicit none
   save
   public
!
   integer I1_PARINT, JATOM
   integer JDOCOR(3, MXCENT), JDOREP(8), JDOSYM(8)
   integer JLUDAS, JLUINT, JLUONE, JLUSOL, JLUSUP
   integer JMXDIF, JMXREP, JNODES, JNOCNT, JNODV, JNOPV
   integer JPRINT, JRETUR, JTASK, JTIMIN, JTKTIM
   integer JTYPE, JSOLVN, J2TYP, JCEDIF, JRELCL
   integer JLEVYL, JTWOCOMP, JECPCALC, JNONREL, I2_PARINT
!
   ! Defined as common, as they are used for mpi parallelization
   common /PARINT/ I1_PARINT, JATOM, JDOCOR, JDOREP, JDOSYM,             &
                   JLUDAS, JLUINT, JLUONE, JLUSOL, JLUSUP, JMXDIF,      &
                   JMXREP, JNODES, JNOCNT, JNODV, JNOPV, JPRINT,        &
                   JRETUR, JTASK, JTIMIN, JTKTIM, JTYPE, JSOLVN,        &
                   J2TYP, JCEDIF, JRELCL, JLEVYL, JTWOCOMP, JECPCALC,   &
                   JNONREL, I2_PARINT
!
   end module parint_mod
