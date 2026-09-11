!-----------------------------------------
!
! Module with variables previously defined
!
! at twocom.h include file
!
!-----------------------------------------

module twocom_mod

   implicit none

   ! Logical variables

   logical :: TCONAB, TCONCD, TPRIAB, TPRICD, DIAGAB, DIAGCD, DIACAB, &
              DIACCD, SHAEQB, SHCEQD, SHABAB, GENAB, GENCD, RPRIAB,   &
              RPRICD, RCNTAB, RCNTCD, SPHRA, SPHRB, SPHRC, SPHRD,     &
              SPHRAB, SPHRCD, DIRAC, LSCALE_GAUNT

   ! Cartesian coordinates

   real*8 :: CORAX0, CORAY0, CORAZ0, CORBX0, CORBY0, CORBZ0, &
             CORCX0, CORCY0, CORCZ0, CORDX0, CORDY0, CORDZ0

   ! Integer variables

   integer :: ISHELA, ISHELB, ISHELC, ISHELD, NHKTA, NHKTB, NHKTC, &
              NHKTD, KHKTA, KHKTB, KHKTC, KHKTD, KHKTAB, KHKTCD,  &
              KCKTA, KCKTB, KCKTC, KCKTD, KCKTAB, KCKTCD, KCKMAX,  &
              NCENTA, NCENTB, NCENTC, NCENTD, MULA, MULB, MULC,    &
              MULD, MULTA, MULTB, MULTC, MULTD, NUCA, NUCB, NUCC,  &
              NUCD, NUCAB, NUCCD, MXUCAB, MXUCCD, NORBA, NORBB,    &
              NORBC, NORBD, NORBAB, NORBCD, NOABCD, NSETA, NSETB,  &
              NSETC, NSETD, NSTRA, NSTRB, NSTRC, NSTRD, JMAXA,     &
              JMAXB, JMAXC, JMAXD, IGENAB, IGENCD, MAXAB, MAXCD,   &
              IAB0X, IAB0Y, IAB0Z, ICD0X, ICD0Y, ICD0Z, NODCAB,     &
              NODCCD, NOPREP, NINTS, NINTMX, IASMAX, IBSMAX,       &
              ICSMAX, IDSMAX, IASTRT, IBSTRT, ICSTRT, IDSTRT, ICEN

   integer :: IPTREP(0:8,2), NINTSR(8)

end module twocom_mod