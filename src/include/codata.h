! ** COMMON block for CODATA
!     CODSET - Codata set of values.
!              Available as: CODATAXX, where XX possible values are
!                            86, 98, 02, 06, 10, 14, 18, 22
!     Fundamental and derived constants declared as real*8 values
!
!#include <codata.h>
      CHARACTER*8 CODSET
      REAL*8      xtang,echarge,hbar,xfmol,umass,pmass,emass,ccm,       &
     &            fermicc,s2thetaw,autk,planck,xtangm10,cvel,alphac,    &
     &            alpha2,xtj,xtkays,xthz,xtev,xkjmol,xkcmol,xtnm,       &
     &            xajoul,autime,tesla,debye,xtkmml,xfamu,xfmp,nmagn,    &
     &            nmagnau,efaumksa,cminv,nm,gfermi,gfac,g1fac
      COMMON/CODVAL/CODSET,xtang,echarge,hbar,xfmol,umass,pmass,emass,  &
     &              ccm,fermicc,s2thetaw,autk,planck,xtangm10,cvel,     &
     &              alphac,alpha2,xtj,xtkays,xthz,xtev,xkjmol,xkcmol,   &
     &              xtnm,xajoul,autime,tesla,debye,xtkmml,xfamu,xfmp,   &
     &              nmagn,nmagnau,efaumksa,cminv,nm,gfermi,gfac,g1fac          
