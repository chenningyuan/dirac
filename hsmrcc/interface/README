  
  File with "know-how" for the HSMRCC program, which is to be connected to the DIRAC suite.
***********************************************************************************************


  DIRAC reads input and makes integral files for the external HSMRCC program of Sham.

 All codes are prefered to be in the Fotran90 language.

-------

       Gamess (after SCF) unformatted files with data/integrals:

      open(UMR1EINT,file=FileMR1EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')
      open(UMR2EINT,file=FileMR2EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')

      $GMS_SRC/"gmsmod.f90"

      and one formatted file with symmetry information:

      open(USYMINF,file=FileSYMINF)
 $GMS_SRC/"syminf.f90" 


     Know-how on hsmrcc stuff
==================================== 

  source/hsmrmain.c - main C code accepting parameters; 

-----------------------------------------------------------------------
 calling SetHSMRfiles(IFile,OFile)

   FileMR1EINT  = TRIM(IFile)//'.F295'   ! binary
   FileMR2EINT  = TRIM(IFile)//'.F296'   ! binary
   FileSYMINF   = TRIM(IFile)//'.F297'   ! text file...

  FileVINF  = TRIM(IFile)//'.VINF'

-----------------------------------------------------------------------

 calling source/hsmrcc.f90  main   HSMRCC (from C)

 calls DrvIntGen (in modules/intgenmod.f90)

-----------------------------------------------------------------------

modules/vinfmod.f90:     open(UVINF,file=FileVINF,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')

modules/gmsmod.f90:      open(UMR1EINT,file=FileMR1EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')
modules/gmsmod.f90:      open(UMR2EINT,file=FileMR2EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')

modules/intgenmod.f90:   open(UMR2EINT,file=FileMR2EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')

--------------------------------------------------------------------------------------------------

File name = LiH.scf_hsmrcc.gms.F297    formatted, sequential access   record = 1
 In source file modules/syminf.f90, at line number 96

 - important symmetry information from Gamess
--------------------------------------------------------------------------------------------------

  Files in working directory:
================================
hsmrcc.exe*              LiH.scf_hsmrcc.gms.FHH   LiH.scf_hsmrcc.gms.TAMP  LiH.scf_hsmrcc.gms.VINF_old
LiH.scf_hsmrcc.gms.COEF  LiH.scf_hsmrcc.gms.FHP   LiH.scf_hsmrcc.gms.TNEW  LiH.scf_hsmrcc.gms.VINT
LiH.scf_hsmrcc.gms.F295  LiH.scf_hsmrcc.gms.FPP   LiH.scf_hsmrcc.gms.TOLD  pgprof.out
LiH.scf_hsmrcc.gms.F296  LiH.scf_hsmrcc.gms.HEFF  LiH.scf_hsmrcc.gms.VIND
LiH.scf_hsmrcc.gms.F297  LiH.scf_hsmrcc.gms.inp   LiH.scf_hsmrcc.gms.VINF

Integrals generator within HSMRCC:

modules/intgenmod.f90

**********************************************************************************************

minmax.h and tind.h include files have been moved into common include/ directory


=====================================================================================


C     ----- INITIALIZE VARIABLES FOR NAMELIST $HSMR -----
C
C     THE FOLLOWING OPTIONS HAVE BEEN IMPLEMENTED
C
c     EXORB
c     *****
c              EXCHANGE ORBITAL PAIRS BEFORE MAKING SUBDIVISION OF
c              FULL ORBITAL INTO CORE, ACTIVE AND VIRTUAL SPACE.
c              THIS INFORMATION IS TO BE PROVIDED IN $EXCHMO CARD
c              AS FOLLOWS.
c
c        FOR RHF/ROHF
c        ------------
c
c              $EXCHMO
c              I1,I2
c              I3,I4
c              ..
c              ..
c              ..
c              $END
c
c        I1,I2 ARE MOS TO BE EXCHANGED
c
C        FOR UHF
C        -------
C
C              $EXCHMO
C              ORBTYP,I1,I2
C              ORBTYP,I3,I4
C              ..
C              ..
C              ..
C              $END
C
C        ORBTYP = +1 FOR ALPHA ORBITALS
C        ORBTYP = -1 FOR BETA  ORBITALS
C        I1,I2 ARE MOS TO BE EXCHANGED
C      
C     NFRZV
C     *****
C              FREEZE LAST NFRZV ORBITALS
C              AFTER CARRYING OUT ANY
C              ORBITAL EXCHANGE ACCORDING
C              TO EXORB OPTION ABOVE
C              (SAME VALUE FOR BOTH ALPHA AND BETA ORBITALS)
c              DEFAULT IS 0, MEANS NO FREEZING
C
C     TRFTOL
C     ******
C              CUTOFF TOLERANCE FOR TRANSFORMED
C              INTEGRALS (DEFAULT=1.0D-10)
c
C     NCOR
C     ****
C              NUMBER OF CORE ORBITALS
C              
C     NACT
C     ****
C              NUMBER OF ACTIVE ORBITALS
C              
C     NVIR
C     ****
C              NUMBER OF VIRTUAL ORBITALS
c
C     NACTEL
C     ******
C              NUMBER OF ACTIVE ELECTRONS (DETERMINES NHOLE AND NPTCL)
C
C     PMAIN
C     *****
c              AFTER MAKING SUBDIVISION OF ORBITAL SPACE INTO
c              CORE, ACTIVE AND VIRTUAL, WE CONSIDER ACTIVE
c              ORBITALS, AND FURTHER DIVIDE THEM INTO MAIN AND
C              INTERMEDIATE SPACE. THE ACTIVE ORBITALS ARE
c              RELABELLED FROM 1 TO NACT.
c
c              PMAIN IS THE LIST OF ACTIVE ORBITALS IN MAIN SPACE. REST OF
c              ACTIVE ORBITALS ARE IN INTERMEDIATE SPACE.
c
c              EXAMPLE :
c
c                 SUPPOSE THERE ARE 10 ACTIVE ORBITALS, AND WE WISH
c                 TO INCLUDE ORBITALS 2, 3 AND 5 AS MAIN, WE SPECIFY
c                 PMAIN INPUT IN $HSMR CARD AS FOLLOWS.
c
c                 $HSMR NCOR=.. NACT=.. NVIR=.. PMAIN(1)=2,3,5 $END
c
c     NOTE :: THE ORBITAL INDICES ABOVE REFER TO SPATIAL ORBITALS
c             DURING THE TRANSFORMATION THEY ARE CONVERTED TO SPIN ORBITALS
c             SUCH THAT FOR A GIVEN SPATIAL ORBITAL, CORRESPONDING ALPHA
c             ORBITAL COMES FIRST FOLLOWED BY BEETA ORBITAL. THIS ORDERING
c             IS ASSUMED IN HSMRCC PROGRAM. SEE SETUP_MRCCORBS FOR DETAILS.
c
c

