*
* Lucia.f : GAS implementing (no pair) relativistic Theory (YES !!!!)
*
      SUBROUTINE DBG_INT_LISTS(NL1,IUB1,NL2,IUB2,ISM2,IDENSI,
     &                         ISPINFREE)
*
* Number of double group integral lists to be constructed
*
*     Jeppe Olsen, Feb. 98
*
* Reordered lists for convenience: We only need the list
* (uu|uu) in an initial or plain SOCI calculation. So this
* is going to be list number 1. Allocation modified as well.
* Modifications in GETINT_DBG incorporated (and elsewhere).
* We are furthermore using real wave functions. This implies
* That complex conjugation symmetry can also be used in the
* quaternion matrix groups C1 and Ci, because their irrep
* basis functions are not generally complex. But if we
* generalize the code to e.g. atomic complex 4- (or 2-)spinors
* or the like, these things have to be distinguished.
* So we will have 3 cases:
*
*   1) Initial or plain SOCI: 1 list
*   2) SOCI + spinor optimization in some way (INS or MCSCF): 3 lists
*   3) Complex wave function: 6 lists
*
* Here, time-reversal symmetry has already been accounted for in
* all cases, so these are non-redundant integral lists.
* At the moment, the first two cases are assumed to be implemented,
* so we will distinguish them by the input keyword DENSI, i.e.
* whether density matrices are to be calculated or not.
*
*     Timo Fleig, June 1999
*
* revised for DIRAC environment
*     Timo Fleig, December 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
*.output
      INTEGER IUB1(2,10),IUB2(4,10),ISM2(3,10)
*
*. ======================================
*. 1e Lists : all four lists constructed
*. ======================================
*
      NL1 = 4
*. 1 : unbarred, unbarred
      IUB1(1,1) = 1
      IUB1(2,1) = 1
*. 2 : unbarred, barred
      IUB1(1,2) = 1
      IUB1(2,2) = 2
*. 3 : barred, unbarred
      IUB1(1,3) = 2
      IUB1(2,3) = 1
*. 4 : barred, barred
      IUB1(1,4) = 2
      IUB1(2,4) = 2
*
*. =========
*. 2e Lists
*. =========
*
*  For the moment, no complex wave function:
      IWFCOMP = 0
*
* The lists of integrals are
*
* 1 : (u u | u u)
*--------------------
* 2 : (u b | u u)
* 3 : (b u | u b)
* 4 : (b u | b u)
*--------------------
* 5 : (u b | u b)
* 6 : (b u | u u)
*--------------------
*
*=====================
* List 1 ( u u ! u u )
*=====================
*
      IUB2(1,1) = 1
      IUB2(2,1) = 1
      IUB2(3,1) = 1
      IUB2(4,1) = 1
*. No Symmetry between 12 and 34
      ISM2(1,1) = 0
      ISM2(2,1) = 0
      ISM2(3,1) = 0
*
      if (IDENSI.ne.0.or.ISPINFREE.eq.0) then
*
*=====================
* List 2 ( u b ! u u )
*=====================
*
        IUB2(1,2) = 1
        IUB2(2,2) = 2
        IUB2(3,2) = 1
        IUB2(4,2) = 1
*. No Symmetry between 12 and 34
        ISM2(1,2) = 0
        ISM2(2,2) = 0
        ISM2(3,2) = 0
*
*=====================
* List 3 ( b u ! u b )
*=====================
*
        IUB2(1,3) = 2
        IUB2(2,3) = 1
        IUB2(3,3) = 1
        IUB2(4,3) = 2
*. No Symmetry between 12 and 34
        ISM2(1,3) = 0
        ISM2(2,3) = 0
        ISM2(3,3) = 0
*
*=====================
* List 4 ( b u ! b u )
*=====================
*
        IUB2(1,4) = 2
        IUB2(2,4) = 1
        IUB2(3,4) = 2
        IUB2(4,4) = 1
*. No Symmetry between 12 and 34
        ISM2(1,4) = 0
        ISM2(2,4) = 0
        ISM2(3,4) = 0
*
        if (IWFCOMP.eq.1) then
*
*=====================
* List 5 ( u b ! u b )
*=====================
*
          IUB2(1,5) = 1
          IUB2(2,5) = 2
          IUB2(3,5) = 1
          IUB2(4,5) = 2
*. Symmetry between 12 and 34
          ISM2(1,5) = 0
          ISM2(2,5) = 0
          ISM2(3,5) = 0
CTEMP     ISM2(3,5) = 1
*
*=====================
* List 6 ( b u ! u u )
*=====================
*
          IUB2(1,6) = 2
          IUB2(2,6) = 1
          IUB2(3,6) = 1
          IUB2(4,6) = 1
*. No Symmetry between 12 and 34
          ISM2(1,6) = 0
          ISM2(2,6) = 0
          ISM2(3,6) = 0
*
        end if
      end if
*
      if (IDENSI.eq.0.and.ISPINFREE.eq.1) then
        NL2 = 1
      else
        NL2 = 4
        if (IWFCOMP.eq.1) NL2 = 6
      end if
*
      RETURN
      END
