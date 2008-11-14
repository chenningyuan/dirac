*   
*   All important symmetry variables not on pointers
*
      character*6 DBGSYM
      logical KRAMERS
      integer   NIRR_PN,DDIM,NIRR_DG,DOUGRP,
     &          IDGTYP,SYDI,SYMDIM,SYMDIM_F,ITOTDM,
     &          NDIM_TOT,IDIMGRP,NSDGIRR,NSPGIRR,
     &          N_SPINR_EXP,NGASPI,NGSOBPI,IMAXDIM,
     &          ITRIA,ITRIADB,ITRIASUM,ITRIADBSUM,
     &          ITRIA_F,ITRIADB_F,ITRIASUM_F,ITRIADBSUM_F,
     &          MORB_DIRAC,IOFPGSM,INVSYM
      dimension IDBGMULT(MXNDGIRR,MXNDGIRR),
     &          IPGMULT(MXNSPIR,MXNSPIR),
     &          IMOSF_SP(MXNDGIRR,2),
     &          ISPSF_MO(MXNDGIRR,2),IADJSYM(MXNDGIRR),
     &          INVELM(MXNDGIRR),IOFPGSM(MXNSPIR)
      dimension SYDI(MXNSPIR),SYMDIM(MXNSPIR),
     &          SYMDIM_F(MXNSPIR),
     &          NGSOBPI(MXPIRR,MXPNGAS),
     &          ITRIA(MXNSPIR),ITRIADB(MXNSPIR),
     &          ITRIA_F(MXNSPIR),ITRIADB_F(MXNSPIR),
     &          ITOTDM(MXNSPIR),NFRO_L(8)
      common /DGSYMM/KRAMERS,NIRR_PN,DDIM,NIRR_DG,DOUGRP,
     &          IDGTYP,SYDI,SYMDIM,SYMDIM_F,ITOTDM,
     &          NDIM_TOT,IDIMGRP,NSDGIRR,NSPGIRR,
     &          N_SPINR_EXP,NGASPI,NGSOBPI,IMAXDIM,
     &          ITRIA,ITRIADB,ITRIASUM,ITRIADBSUM,
     &          IDBGMULT,IPGMULT,IMOSF_SP,ISPSF_MO,
     &          IADJSYM,INVELM,NFRO_L,ITRIA_F,
     &          ITRIADB_F,ITRIASUM_F,ITRIADBSUM_F,
     &          MORB_DIRAC,IOFPGSM,INVSYM
      common /DGSYMMC/ DBGSYM

*
* KRAMERS	: true: time-reversal symmetry invariant formalism
* DDIM	 	: double dimension (of SYMDIM(x))
* NIRR_DG	: number of double group irreps
* DOUGRP	: double group
* IDGTYP	: type of double group
* NIRR_PN	: number of spatial irreps
* SYDI		: symmetry dimensions of spatial irreps
* SYMDIM	: two SYDI's added (depending on double group)
* NDIM_TOT	: Total number of spatial functions (orbitals)
* IDIMGRP	: dimension group for int trafo (1 SYMDIM)
* NSDGIRR	: number of sub double group irreps
* NSPGIRR	: number of sub point group irreps
* N_SPINR_EXP   : number of spinor expansions
* NGASPI	: number of GAS spaces (in picasso)
* NGSOBPI	: matrix with number orbitals in GAS and symmetries
* IMAXDIM	: largest dimension of spatial irreps (number from 
*                 ordering
*                    in Koster, Dimmock, Wheeler, Statz)
* ITRIA         : triangle lengths of point group blocks
* ITRIASUM      : sum of these
* ITRIADB       : triangle lengths of two connected blocks in double 
*                 group
* ITOTDM        : total orbital space dimensions per spatial irrep
* IDBGMULT      : Double Group multiplication table
* IPGMULT       : Point group mult. table
* IMOSF_SP      : Translation table from MO symmetry to spinor symmetry
*                 IMOSF_SP(IOBSM,IAB) Gives the spinor symmetry
*                 of an orbital of symmetry IOBSM times an alpha (IAB=1)
*                 or and beta spin-function (IAB = 2 )
* ISPSF_MO      : Translation table from spinor symmetry to MO symmetry
*                 ISPSF(ISPSM,IAB) gives the symmetry of an MO that
*                 combined with an alpha(IAB=1) or beta(IAB=2)
*                 spin function give a spinor with symmetry ISPSM
*
* IADJSYM       : IADJSYM(ISPSYM) is the adjoint symmetry of ISPSYM
* INVELM        : Inverse elements
* NFRO_L 	: Number of frozen (external Fock) orbitals per symmetry
* MORB_DIRAC	: Number of spinors per non-redundant fermion irrep 
*                 (DIRAC)
* IOFPGSM	: Offset for number of orbitals in spatial symmetry
* INVSYM        : 0 for groups without, 1 for groups with inversion s.
* _F (general)	: Full info with respect to Fock type calculations
*
* The irreps are organized in the usual way with the boson irreps first
