      CHARACTER*6 DBGSYM
      dimension IDBGMULT(8,8),IPGMULT(4,4),IMOSF_SP(8,2),ISPSF_MO(8,2),
     &          IADJSYM(8),INVELM(8)
      COMMON/DBGPSYM/IDBGMULT,IPGMULT,NPGIRR,NDBGIRR,IMOSF_SP,
     &              ISPSF_MO,IADJSYM,INVELM
* Contents of DBGPSYM
*    IDBGMULT : Double Group multiplication table 
*    IPGMULT  : Point group mult. table
*    NPGIRR   : Number of point group irreps
*    NDBGIRR  : Number of irreps in double group
*    IMOSF_SP : Translation table from MO symmetry to spinor symmetry
*               IMOSF_SP(IOBSM,IAB) Gives the spinor symmetry
*               of an orbital of symmetry IOBSM times an alpha (IAB=1)
*               or and beta spin-function (IAB = 2 )
*    ISPSF_MO : Translation table from spinor symmetry to MO symmetry
*               ISPSF(ISPSM,IAB) gives the symmetry of an MO that
*               combined with an alpha(IAB=1) or beta(IAB=2)
*               spin function give a spinor with symmetry ISPSM
*
*   IADJSYM   : IADJSYM(ISPSYM) is the adjoint symmetry of ISPSYM
*   INVELM    : Inverse elements
*
* The irreps are organized in the usual way with the boson irreps first 
