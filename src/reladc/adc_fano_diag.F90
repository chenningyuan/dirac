MODULE adc_fano_diag

CONTAINS

  SUBROUTINE fano_wcondat_s(LUN,FILENAME,IRECL,DESREP,lenmat, &
                            finalpos,fin_max,iw)

    IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!   Writes explicit irep names and spinor numbers for the various
!   final state configurations belonging to symmetry DESREP.
!   Avoiding too complicated things we fix IRECL in this routine
!   no matter which value is transferred to this routine. Therefore
!   IRECL is an OUTPUT parameter.
!   Configurations of the final states are kicked out, they won't
!   won't be used in FanoADC.
!
!---------------Calling variables--------------------------------------
!
    CHARACTER(8), INTENT(IN)                 :: FILENAME
    INTEGER, INTENT(IN)                     :: LUN,DESREP,lenmat
    INTEGER, INTENT(OUT)                    :: irecl
    INTEGER, INTENT(IN)                     :: fin_max
    INTEGER, INTENT(IN), DIMENSION(fin_max) :: finalpos
    INTEGER, INTENT(IN)                     :: iw
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!                                                                 
    INTEGER       :: lcount                                              
    INTEGER       :: satcount, fincount
    CHARACTER*44  :: FIELD44                                        
    INTEGER       :: bl_conf                                             
    CHARACTER*4   :: bl_irna,bl_abs                                  

    INTEGER       :: a, k , l
    INTEGER       :: arep, krep, lrep, klrep
    INTEGER       :: aoff, koff, loff
    INTEGER       :: kmin
                                                                  
!                                                                 
!---------------Executable code--------------------------------------
!                                                                 
!  Code for ** single ionization **:  1h/2h1p                     
!                                                                 
! fix record length for the direct access file and set some       
! blank variables. Then open file and initialize configuration    
! counter.                                                        
! *ATT* Is has turned out advantageous that storing the *absolute*
! spinor numbers is superior to storing the relative ones for the 
! subsequent postprocessing!                                      
!                                                                 
! this will be followed by the single ionization case for the moment
! (1h and 2h1p) configurations                       
!                                                                 
    IRECL = 44                                                  
    BL_CONF = 0                                                 
    BL_IRNA = '    '                                            
    BL_ABS  = 'abs.'                                            
    OPEN(LUN,FILE=FILENAME,ACCESS='DIRECT',RECL=IRECL,    &
   &     STATUS='UNKNOWN')                                      
    LCOUNT = 1                                                  
    satcount = 0
    !fincount = 1
                                                                
!_________________________ SIPS ______________________________    
!|                                                                
!|                                                                
!|                                                                
      WRITE(iw,*)
      WRITE(IW,*) 'Writing FanoSIP configs to ',FILENAME            
      WRITE(IW,*) 'Final state symmetry: ',DESREP               
      WRITE(iw,*)
                                                                  
! ** start with 1h configurations                                 
                                                                  
      DO k=1,no(DESREP)                                         
        koff=io(DESREP)+k                                       
        WRITE(FIELD44,381)                                &      
                k,repna(DESREP),                          &      
                koff,repna(DESREP),                       &      
                bl_conf,bl_irna,                          &      
                bl_conf,bl_irna                                 
        WRITE(LUN,REC=lcount) FIELD44                           
        lcount = lcount + 1                                     
      ENDDO                                                     
                                                                
! ** continue with 2h1p configurations                            
! ** write absolute spinor numbers                                
                                                                  
      klloop:DO  KLREP=1,NREP                                        
        AREP=MULTB(DESREP,KLREP+NREP,2)                         
        lloop:DO LREP=1,NREP                                       
          KREP=MULTB(LREP,KLREP+NREP,2)                          
          IF(KREP.LT.LREP) CYCLE lloop                           
          DO L=1,NO(LREP)                                        
            LOFF=IO(LREP)+L                                      
            KMIN=1                                               
            IF(KREP.EQ.LREP) KMIN = L + 1                        
            DO K=KMIN,NO(KREP)                                   
              KOFF=IO(KREP)+K                                    
              aloop:DO A=1,NV(AREP)                                    
                AOFF=IO(NREP+1)+IV(AREP)+A                       
                satcount = satcount + 1
                IF (ANY(finalpos == satcount)) THEN
                  CYCLE aloop
                END IF
!                WRITE(iw,*) 'satcount = ', satcount
!                 WRITE(FIELD44,381)                               
!    &              K,REPNA(KREP),                                 
!    &              L,REPNA(LREP),                                 
!    &              A,REPNA(AREP),                                 
!    &              BL_CONF,BL_IRNA                                
                WRITE(FIELD44,381)                         &      
   &              KOFF,REPNA(KREP),                        &      
   &              LOFF,REPNA(LREP),                        &      
   &              AOFF,REPNA(AREP),                        &      
   &              BL_CONF,BL_IRNA                                
                WRITE(LUN,REC=lcount) FIELD44                    
                LCOUNT=LCOUNT+1                                  
              END DO aloop                                             
            ENDDO                                                
          ENDDO                                                  
        END DO lloop
      END DO klloop
                                                                  
! ** file is open only within this subroutine                     
                                                                  
      CLOSE(LUN)                                                
                                                                  
! ** correct configuration counter                                
                                                                  
      LCOUNT = LCOUNT - 1                                       
      IF(LCOUNT.NE.lenmat) THEN                                
        WRITE(IW,*) 'LCOUNT,lenmat:',LCOUNT,lenmat            
        CALL QUIT('Error in fano_wcondat_s for single ionizations')  
      ENDIF                                                     
!|                                                                
!|                                                                
!|____________________________________________________________    
                                                                  
  381  FORMAT(4(I3,'(',A4,')  '))                                  
                                                                  
  END SUBROUTINE fano_wcondat_s             

END MODULE adc_fano_diag
