      SUBROUTINE INTTRA(I1234,I12,I34,ICCSM,IOPSMPN,IOPSMDG,NAOFSM,
     &                  NSMSP,NSMAO,MULTSP,MULTAO,
     &                  NLISTO,I12O,I34O,I1234O,ICCSMO,
     &                  ISPITRA,ISPITRA_RC,IAOSFSP,ISPSFAO,LSPINOR )
*
* General integral transformation for CI with general
* spinors
*
* On input we must have nformation about
*     1  The list of input ao. integrals
*     2  sp-ao transformation matrix( matrices)
*     3  The required list of transformed integrals
*     4  symmetry information,
*
* +++
*  1
* +++
*
* The list a.o. integrals are specified by
* a) The number of a.o per sym  NAOFSM
* b) The symmetry of the operator in the pointgroup ISMPN
* c) The permutational symmetry of the 4 indeces as specified by
*    I1234,I12,I34
*
* I1234 : Symmetry between particle 1 and 2
*
*     1 =>   (ij!kl) = (kl!ij)
*    -1 =>   (ij!kl) =-(kl!ij)
*     0 =>  No symmetry between particle 1 and 2
*
* I12   : Symmetry between first two indeces (indeces of particle 1 )
*
*    1 => (ij!kl) = (ji!kl)
*   -1 => (ij!kl) =-(ji!kl)
*    0 => No symmetry between the first two indeces
*
* I34   : Symmetry between last two indeces (indeces of particle 2)
*
*    1 => (ij!kl) = (ij!lk)
*   -1 => (ij!kl) =-(ij!lk)
*    0 => No symmetry between the last two indeces
*
*    We will not assume anything about the specific order of
*    the a.o integral list or the block sizes ...
*    We will instead invoke the a.o integrals through subroutine
*    calls where these details for specific implementations are
*    hidden
*
* +++
*  2
* +++
*
* There can be more than one sp-ao transformation matrix
* ( barred and unbarred spinors for example), and each
* transformation matrix is in general complex
*
*
* So we have as input
*  NSPAO : Number of mo-ao transformations
*  NSPFSMG(ISYM,IEXP) : Number of spinors of symmetry ISYM in
*  expansion IEXP, general expansion. The symmetry refers to
*  general spinor symmetry, say double groups.
*  NSPFSMR(ISYM,IEXP) : Number of spinors of symmetry ISYM in
*  expansion IEXP, restricted expansion,The symmetry refers to
*  general spinor symmetry, say double groups.
*  CSPAO : List of all SP-AO transformation matrices
*  ICSPAO : ICSPAO(IRI,IEXP)  Offset for Real (IRI = 1)
*  or imaginary part (IRI = 2 ) of SPAO expansion IEXP
*
* The  real or imaginary part of a given sp-ao expansion is
* ordered as
*
* Loop over symmetries
*  Loop over spinors of given symmetry
*   Loop over alpha functions * spatial functions of correct symmetry
*   End alpha functions * spatial functions of correct symmetry
*   Loop over beta  functions * spatial functions of correct symmetry
*   End beta  functions * spatial functions of correct symmetry
*  End of loop over spinors of given symmetry
* End of loop over symmetries
*
* The real/imaginary part of a given spao expansion  contains thus
* in each column the expansion of a given spinor. In the given column
* the are first alpha spin orbitals and the beta spin orbitals.
* The symmetries of the ao orbitals that combine with a given
* spin function to give a certain double group symmetry is
* given by the symmetry tables ( see these )
*
*
*
* +++
*  3
* +++
*
* I guess the following is obsolete, since the transformation
* to restricted sets is realized in a different way.
* One defines restricted spinor transformations,
* and includes in the input setup, the relevant types
* of reduced inegrals lsists
*
* N_GEN_IND
*
*    In general one is not interested in a complete list of integrals
*    rather, one or several indeces are restricted to a subset
*    This is specified by the N_GEN_IND parameter. Using the
*    notation that g stands for a general index, whereas
*    r stands for a restricted index, we have the following
*
*    0 => (rr!rr) integrals ( for a CI in restricted orbital space)
*    1 => (rr!rg) integral list ( for fisrst order casscf)
*    2 => (rr!gg)+(rg!rg) integral list ( for second order mcscf
*          and linear response, CASPT2)
*    3 => for quadratic response (not extremely relevant)
*    4 => for general CI in all orbitals
*
* Well, to make life (atleast programming) really complicated
* we are in general interested in more than one output integral list
* so another input parameter is
*
* NLISTO : Number of output integral lists
*
* The  output list are defined by
*
* I1234O(ILIST), I12O(ILIST), I34O(ILIST),and ICCSMO
* giving permutational
* symmetries of output list ILIST. The newcomer is ICCSMO
* that indicated whether a list ishould be restricted
* by complex conjugation restrictions, (12!34) = (21!43)*.
* This will define the index restrictions on the output integrals
*
*
* ISPITRA(INDEX,ILIST) : Index number INDEX(1,2,3,4) in output
* list ILIST refers to spao expansion ISPITRA (INDEX,ILIST)
*
* ISMDGG : Symmetry of the operator in the double group
* representation
*
* +++
*  4
* +++
*
* The information  about symmetry is
*
* NSMSP : Number of symmetries of spinors ( number of iireps in
*         double group )
* NSMAO : Number of symmetries of atomic orbitals ( Number of
*         irreps of pointgroup)
* ISPSFAO(ISPSM,IAB) : Spinor symmetry X spinfunction => ao symmetry
*                      The point group symmetry to be multiplied
*                      with alpha (IAB=1) or beta (IAB = 2)  to
*                      give spinor symmetry ISPSM
* IAOSFSP(IAOSM,IAB) : ao symmetry X spin function => spinor symmetry
*                      The spinor  symmetry to be multiplied
*                      with alpha (IAB=1) or beta (IAB = 2)  to
*                      give point group symmetry ISPSM
*
* Timo Fleig + Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. input
*.. Symmetry information
      INTEGER ISPSFAO(NSMSP,2), IAOSFSP(NSMAO,2)
*.. Information about output integral lists
      INTEGER ISPITRA(4,NLISTO)
      INTEGER I12O(NLISTO),I34O(NLISTO),I1234O(NLISTO)
*  NSPAO : Number of mo-ao transformations
*  NSPFSMG(ISYM,IEXP) : Number of spinors of symmetry ISYM in
*  expansion IEXP, general expansion. The symmetry refers to
*  general spinor symmetry, say double groups.
*  NSPFSMR(ISYM,IEXP) : Number of spinors of symmetry ISYM in
*  expansion IEXP, restricted expansion,The symmetry refers to
*  general spinor symmetry, say double groups.
*  CSPAO : List of all SP-AO transformation matrices
*  ICSPAO : ICSPAO(IRI,IEXP)  Offset for Real (IRI = 1)
*  or imaginary part (IRI = 2 ) of SPAO expansion IEXP
*
C      INTTRA(I1234,I12,I34,ICCSM,IOPSMPN,NAOFSM,
C     &                  NSMSP,NSMAO,MULTSP,MULTAO,
C     &                  NLISTO,I12O,I34O,I1234O,ICCSMO,
C     &                  ISPITRA,ISPITRA_RC,IAOSFSP,ISPSFAO         )
*
      WRITE(6,*) ' ISPITRA'
      CALL IWRTMA(ISPITRA,4,NLISTO,4,NLISTO)
*
      CALL SPITRA2(I1234,I12,I34,NSMAO,MULTAO,NSMSP,MULTSP,
     &             IOPSMPN,IOPSMDG,NLISTO,I1234O,I12O,I34O,ICCSMO,
     &             ISPITRA,ISPITRA_RC,IAOSFSP,ISPSFAO,LSPINOR)
      RETURN
      END
