*   
*   Further important input information NOT on symm.inc
*
      integer IOPER,IPRINT,MOL4_3
      character*3 ARCHI
      logical ORBFILE
      common /INPINFO/IOPER,IPRINT,ORBFILE,ARCHI
*
*  IOPER:  Type of operator for integral transformation
*          (1) transform only spin-free ints
*          (2) transform only spin-orbit ints
*          (3) transform both types
*  IPRINT: Print level for transformation info
*          currently 0-3, 5-10 for massive debug output
*  ARCHI:  architecture flag for SGI/IBM (at the moment)
*  ORBFILE: For MOLCAS CMOAO vector
*             .true.  input file is formatted ASCII vector (INPORB)
*                                                        or LUMOUT
*             .false. input file is JOBIPH file from RASSCF
*
