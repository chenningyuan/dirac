      SUBROUTINE SET_HOP_DBG(IHTYPE,IPRNT)
*
* Construct orbital excitation types for different choices 
* of Hamiltonian operator, double group version
*
* 
* IHTYPE = 1 : Include only Delta MK = 0 terms
* IHTYPE = 2 : Include Delta Mk = 0 for two-electron part 
*              and complete one-electron part
* IHTYPE = 3  : Include all terms for general relativistic Hamiltonian
*
#include "implicit.inc"
*  
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "cprnt.inc"
#include "orbinp.inc"
*. Output
#include "ctcc.inc"
*
*. Number of spin-orbital exctation operators in H
      CALL MEMMAN(KLSOBEX_TEMP,4*NGAS,'ADDS  ',1,'SPOBET')
      IFLAG = 1
      CALL GET_HX_RELA(IFLAG,IHTYPE,NOBPT,NHX,
     &                 IFTONE,IFTTWO,ILTTWO,WORK(KLSOBEX_TEMP),
     &                 DUMMY,DUMMY,DUMMY,I2DIATP,DUMMY,DUMMY,DUMMY,
     &                 IPRNT)
      NSPOBEX_TP = NHX
*
*. and the actual spin-orbital excitations
      CALL MEMMAN(KLSOBEX,4*NGAS*NSPOBEX_TP,'ADDS  ',1,'SPOBEX')
      call memman(KLIHIND,4*NSPOBEX_TP,'ADDS  ',1,'IHINDX')
      call memman(KSIGNNHX,NSPOBEX_TP,'ADDS  ',2,'SIGNDX')
      call memman(KLABEXTP,4*NSPOBEX_TP,'ADDS  ',1,'LABXTP')
      call memman(KLOP_REO,4*NSPOBEX_TP,'ADDS  ',1,'OP_REO')
      call memman(KSIOPREO,NSPOBEX_TP,'ADDS  ',2,'SIOPRE')
      IFLAG = 0
      CALL GET_HX_RELA(IFLAG,IHTYPE,NOBPT,NHX,
     &                 IFTONE,IFTTWO,ILTTWO,WORK(KLSOBEX_TEMP),
     &                 WORK(KLSOBEX),WORK(KLIHIND),WORK(KSIGNNHX),
     &                 I2DIATP,WORK(KLABEXTP),WORK(KLOP_REO),
     &                 WORK(KSIOPREO),IPRNT)
*
*. Space for info in the different parts of Hamiltonian 
* KLLSOBEX : Length of this block
* KLIBSOBEX : Offset for this block
* KLSPOBEX_AC : Is this spinorbital excitation active 
*
      CALL MEMMAN(KLLSOBEX,NSPOBEX_TP,'ADDS  ',1,'LSPOBX')
      CALL MEMMAN(KLIBSOBEX,NSPOBEX_TP,'ADDS  ',1,'LSPOBY')
      CALL MEMMAN(KLSPOBEX_AC,NSPOBEX_TP,'ADDS  ',1,'SPOBAC')
*. ALl spinorbital excitations are initially active
      IONE = 1
      CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TP)
*. Dimensions for total symmetric operators, sym = 1
      ITOTSYM = 1
      MSCOMB_T = 0 
      CALL IDIM_TCC_DBG(WORK(KLSOBEX),NSPOBEX_TP,ITOTSYM,
     &           MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &           WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &           MSCOMB_T,MX_TBLK_AS,IFTONE,IFTTWO,ILTTWO,
     &           N1ELINT,N2ELINT,IPRSIG)
*
                 MX_TBLK_MX = MX_TBLK            
                 MX_ST_TSOSO_MX =  MX_ST_TSOSO 
                 MX_ST_TSOSO_BLK_MX = MX_ST_TSOSO_BLK
*
C     IDIM_TCC_DBG(ITSOSO_TP,NTSOSO_TP,ISYM,
C    &           MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
C    &           LTSOSO_TP,IBTSOSO_TP,IDIM_T,MSCOMB_CC,
C    &           MX_TBLK_AS,IFTONE,IFTTWO,ILTTWO,
C    &           N1ELINT,N2ELINT,IPRCC)
*
      RETURN 
      END
