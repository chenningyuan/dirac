!    -- FILE: argoscom.h --
!    This is for the ECP routine
      INTEGER IGRPAG, ECPNUC, ECPNONT(600),                             &
     &        ndptag, MAXANG(60), irowmax(60,60), ireparg(60,60,100)   
      COMMON /argoscom/ IGRPAG, ECPNUC, ECPNONT,                        & 
     &                  ndptag, MAXANG, irowmax, ireparg           
!
      INTEGER         ECPINT_PRINT, ECPINT_DEBUG
      COMMON /ECPINT/ ECPINT_PRINT, ECPINT_DEBUG
!
      INTEGER ndptags(0:7)
      data ndptags/ 0,0,0,0,1,1,1,7/
!    -- end of argoscom.h --
