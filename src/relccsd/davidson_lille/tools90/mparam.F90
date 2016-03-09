   module mparam
!
!***********************************************************************
!     ----- module mparam : Parameters definitions -----
!=======================================================================
!     ----- If needed these parameters can be increased       -----
!     ----- An extensive description of the various variables -----
!     ----- is provided in the user manual epciso.ps          -----
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   use kindvar
   integer MAXAT,MAXAO,MAXMO,MAXACT,MXCLAS,MXGOP,MAXOP,MXET,MAXCF, &
           MXSYM,MXFOC,MXDP1,MXTR,MXSO,MAXCON,MAXCOR,MAXDAV, &
           MAXEXC,MXCRES,MAXLZ,MAXGRP,MAXGO,MAXSH, &
           MAXAZI,LXXX,LCTROU
!
   parameter(MAXAT = 20, &
             MAXMO = 600, &
             MAXGO = 2048, &
             MAXSH = 70,  &
             MAXOP = 9,     &
             MXGOP = 40,  &
             MXET  = 150,  &
             MAXDAV = 1500, &
             MAXCF = 1400000, &
             MAXCON = 700000,  &
             MAXCOR = 10,  &
             MXSYM = 20,  &
             MAXEXC = 8,  &
             MAXAZI = 4,  &
             LCTROU = 16,  &
!     ----- MXCRES = Combinaison(MAXEXC,MAXEXC/2) -----
             MXCRES = 1500,  &
             MAXGRP = 20,  &
             MAXLZ = 4096)
!
!     ----- name of symmetry operations and irreducible representations -----
!
   integer (kind=K_I4), parameter :: LENOP=16,LENIR=4
!
!     ----- molecule type (general, atom, linear) -----
!
   integer (kind=K_I4), parameter :: IMOGEN=0,IATOM=1,IMOLIN=2
!
!     ----- properties ------
!
!    PDIP = electric dipole moment  (r)
!    PQUA = electric quadrupole moment rr
!    PMAG = magentic dipole moment (l)
!    PLAN = Lande's factor g
!    POVL = overlap between bras and kets
!    PPOP = "Natural orbital" and Mulliken population analysis
!
   integer, parameter :: MXPROP=4,PDIP=1,PQUA=2,PMAG=3,PLAN=4
!
   parameter(MXFOC = MAXMO*(MAXMO+1)/2,   &
             MXDP1 = MAXCF+1,  &
             MXTR  = MAXEXC*MAXCF,  &
            MAXACT = MAXMO,  &
            MXCLAS = MAXMO,  &
             MXSO = MAXMO*2)
end module mparam
