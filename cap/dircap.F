C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sajeev */
      SUBROUTINE SAJEEV(
     &   basis_n,ang_pow,ord_expon,
     &   ord_coe,r_ms,PE_AOr,PE_AOi,
     &   PEXR,PEXI,PEYR,PEYI,PEZR,PEZI,
     &   WORK,LWORK,NBASF,NANGPOW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Main program from Sajeev. Adapted into DIRAC enviroment
C  as the subroutine.
C
C On input: 
C              basis_n
C              ang_pow
C              ord_expon
C              ord_coe ...
C
C              NBASF
C              NANGPOW
C
C On output: 
C            complex CAP integrals written into formatted file, or test
C            calculation of overlap integrals only for comparison with
C            DIRAC
C
C  called from DOSAJEEV
C
C  Written by Sajeev, adapted to DIRAC by Miro Ilias, Tel Aviv, 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbcap.h"
#include "capint.h"
      
      DIMENSION WORK(LWORK)

CMI    ... arrays are entering subroutine, check the size !!!!!
C     integer *4 basis_n(200)
      integer *4 basis_n(NBASF)
C     integer *4 ang_pow(250,3)
      integer *4 ang_pow(NANGPOW,3)
      integer *4  N_B   
C     real *8 ord_expon(250),ord_coe(250)
      real *8 ord_expon(NANGPOW), ord_coe(NANGPOW)
      real *8 R_m(4),R_n(4)
C     real *8 R_ms(200,3)
      real *8 R_ms(NBASF,3)

C     COMPLEX *16 PE_AO(25000) ! former variable split into real and
C     imaginary parts
CMI   REAL*8 PE_AOr(25000),PE_AOi(25000)
      REAL*8 PE_AOr(NNBBASX),PE_AOi(NNBBASX)
      REAL*8 PEXR(*),PEXI(*),PEYR(*),PEYI(*),PEZR(*),PEZI(*)
      real *8 BETA0,PE0
      
CMI   LOGICAL ISFILE, FNDLAB,ARELABS
CMI   LOGICAL*1 ISFILE, FNDLAB,ARELABS,IS_OPEN
      LOGICAL ISFILE, FNDLAB,ARELABS,IS_OPEN

      COMPLEX *16 THETA,THETA2,CI,BETA

      CHARACTER*10 PEMAT_NAME
      CHARACTER*8 LBL

      CALL QENTER('SAJEEV')

      IF (IPRCAP.GE.2) THEN
        CALL HEADER('*** Output from SAJEEV ***',-1)
        CALL FLSHFO(LUPRI)
      ENDIF

      CALL INPREAD(N_B,ord_expon,ord_coe,R_ms,basis_n,
     &         ang_pow,CAPOVLP,IPRCAP,NBASF,NANGPOW)

CMI ... check the presence of PEMAT in order to decide whether to run 
CMI   or not the CAP integral evaluation
      ARELABS = .FALSE.
      INQUIRE(FILE='PEMAT',EXIST=ISFILE,OPENED=IS_OPEN)
      IF (ISFILE) THEN
        IF (.NOT.IS_OPEN) THEN
          CALL OPNFIL(LUPEMAT,'PEMAT','OLD','SAJEEV')
        ENDIF
        REWIND(LUPEMAT)
        IF (CAPOVLP) THEN
           LBL = 'CAP_OVLP'
        ELSE
         IF  (DONRCAP) THEN
           LBL = 'CAP_RE  '
         ELSE
           LBL = 'ZCAPD1I'
         ENDIF
        ENDIF
        IF (FNDLAB(LBL,LUPEMAT)) ARELABS = .TRUE.
      ENDIF
      IF (.NOT.ARELABS) THEN
CMI   ... evaluate CAP complex atomic untegrals
       IF (IPRCAP.GE.0) THEN
         WRITE(LUPRI,'(/2X,A)')
     &   'CAP integrals are not present, going'//
     &   ' to calculate and write them into the PEMAT file.'
         CALL FLSHFO(LUPRI)
       ENDIF

CMI   ... evaluate integrals - fill array with them
       CALL AOINTEGRAL(N_B,ord_expon,ord_coe,R_ms,basis_n,ang_pow,
     &         PW,PE_AOr,PE_AOi,PEXR,PEXI,PEYR,PEYI,PEZR,PEZI,
     &         NBASF,NANGPOW)
        
C ============================================================
CMI     ... write complex atomic integrals into the file
C
C ... alter the file name when calculating CAP overlap integrals !
C
C ============================================================
CMI      open(unit=LUPEMAT,file='PEMAT',
C        open(unit=LUPEMAT,file=PEMAT_NAME,
C    &    STATUS='UNKNOWN',FORM='FORMATTED')
        
C        CALL OPNFIL(LUPEMAT,'PEMAT','UNKNOWN','SAJEEV') 

CMI      do ii=1,n_b
C           do jj=1,n_b
C              IPU=IPOINTIJ(ii,jj)
C              if (jj.gt.ii)go to 301
C              write(LUPEMAT,*) ii,jj,CMPLX(PE_AOr(IPU),PE_AOi(IPU))
C           end do
C301        continue
C        end do

         REWIND(LUPEMAT)
         IF (CAPOVLP) THEN
            LBL='CAP_OVLP'
            CALL NEWLAB(LBL,LUPEMAT,LUPRI)
            CALL WRITT(LUPEMAT,NNBBASX,PE_AOr)
            WRITE(LUPRI,'(2X,A,A)')
     &      'CAP integrals generated, PEMAT label:',LBL
         ELSE
           IF (DONRCAP) THEN
            IF (FNDLAB('EOFLABEL',LUPEMAT)) THEN 
              BACKSPACE(LUPEMAT)
              LBL='CAP_RE  '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PE_AOr)
              WRITE(LUPRI,'(2X,A,A)')
     &      'CAP nonrel. integrals generated, PEMAT label:',LBL
              LBL='CAP_IM  '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PE_AOi)
             WRITE(LUPRI,'(2X,A,A)')
     &      'CAP nonrel. integrals generated, PEMAT label:',LBL
            ENDIF
           ELSE 
C             ... four-component CAP rel. operator 
            IF (FNDLAB('EOFLABEL',LUPEMAT)) THEN 
             BACKSPACE(LUPEMAT)
             IF (.NOT.CAPDIPVEL) THEN
              LBL='XCAPD1R '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PEXR)
              WRITE(LUPRI,'(2X,A,A)')
     &      'CAP relativistic integrals generated, PEMAT label:',LBL
              LBL='XCAPD1I '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PEXI)
              WRITE(LUPRI,'(2X,A,A)')
     &      'CAP relativistic integrals generated, PEMAT label:',LBL
              LBL='YCAPD1R '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PEYR)
              WRITE(LUPRI,'(2X,A,A)')
     &        'CAP relativistic integrals generated, PEMAT label:',LBL
              LBL='YCAPD1I '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PEYI)
              WRITE(LUPRI,'(2X,A,A)')
     &        'CAP relativistic integrals generated, PEMAT label:',LBL
              LBL='ZCAPD1R '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PEZR)
              WRITE(LUPRI,'(2X,A,A)')
     &        'CAP relativistic integrals generated, PEMAT label:',LBL
              LBL='ZCAPD1I '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PEZI)
              WRITE(LUPRI,'(2X,A,A)')
     &       'CAP relativistic integrals generated, PEMAT label:',LBL
             ELSE
CMI          ... CAP dipole velocity integrals
              LBL='XCAPDVE '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PEXR)
              WRITE(LUPRI,'(2X,A,A)')
     &      'CAP relativistic integrals generated, PEMAT label:',LBL
              LBL='YCAPDVE '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PEYR)
              WRITE(LUPRI,'(2X,A,A)')
     &      'CAP relativistic integrals generated, PEMAT label:',LBL
              LBL='ZCAPDVE '
              CALL NEWLAB(LBL,LUPEMAT,LUPRI)
              CALL WRITT(LUPEMAT,NNBBASX,PEZR)
              WRITE(LUPRI,'(2X,A,A)')
     &      'CAP relativistic integrals generated, PEMAT label:',LBL
             ENDIF
            ENDIF    
         ENDIF
         ENDIF
         CALL NEWLAB('EOFLABEL',LUPEMAT,LUPRI)

C        CLOSE(LUPEMAT,STATUS='KEEP')

         IF (IPRCAP.GE.1) THEN
          WRITE(LUPRI,'(/2X,A)')
     &    'PEMAT file with CAP integrals ready.'
         ENDIF

      ELSE
C================================================================================================
C
C================================================================================================
       WRITE(LUPRI,'(2X,A)')
     & 'SAJEEV: PEMAT file with desired CAP integrals is'//
     & ' present, not calculating them.'
      ENDIF

      CALL QEXIT('SAJEEV')
      RETURN
      END
      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck inpread */
      SUBROUTINE INPREAD(
     &  N_B,ord_expon,ord_coe,R_ms,basis_n,ang_pow,
     &  CAPOVLP,IPRINT,NBASF,NANGPOW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Read the GAMESS-like input file with the basis set
C
C  Parameters:
C
C           CAPOVLP - calculate only overlap integrals (for comparison
C           with DIRAC)
C
C  Called from: SAJEEV
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER(MSIZE=400)
#include "implicit.h"
#include "priunit.h"
      
CMI   static allocations are here !!!
      CHARACTER *1 shell_type,BLNK(6)
C     CHARACTER *2 shell_p(3),shell_d(6),Shell_arr(200)   
      CHARACTER *2 shell_p(3),shell_d(6)
      CHARACTER *3 shell_f(10)
      CHARACTER *4 shell_g(15), Shell_arr(MSIZE)   

      CHARACTER *8 IA

CMI   character *3 Shell_basis(200)
CMI   character *3 Shell_basis(MSIZE)
      character *4 Shell_basis(MSIZE)

C     integer *4 Shell_atom(200)
      integer *4 Shell_atom(MSIZE)
      
      integer *4 n_con_guass,nj,g_n,sh,n_c_g

!   ... powers on x,y,z
      integer *4 angu_pow_S(1,3),angu_pow_P(3,3)
      integer *4 angu_pow_D(6,3),angu_pow_F(10,3),angu_pow_G(15,3)

      integer *4 n_b,ia_p

C     integer *4 ang_pow(250,3)
      integer *4 ang_pow(NANGPOW,3)
      integer *4 n_pg_atom(50)
      integer *4 k_m(3),k_n(3)
C     integer *4 basis_n(200)
      integer *4 basis_n(NBASF)
      
      real *8 re(50,3),expon,coe
C     real *8 ord_expon(250),ord_coe(250)
      real *8 ord_expon(NANGPOW),ord_coe(NANGPOW)
      real *8 R_m(4),R_n(4)
C     real *8 R_ms(200,3)
      real *8 R_ms(NBASF,3)

      LOGICAL CAPOVLP,DIRAC_NRM

      DATA BLNK /'E','S','P','D','F','G' /

      DATA shell_p /'PX','PY','PZ'/
      DATA shell_d /'XX','YY','ZZ','XY','XZ','YZ'/
      DATA shell_f /'XXX','YYY','ZZZ','XXY','XXZ','YYX','YYZ',
     &              'ZZX','ZZY','XYZ'/

      DATA shell_g /'XXXX','YYYY','ZZZZ','XXXY','XXXZ','XYYY','YYYZ',
     &   'XZZZ','YZZZ','XXYY','XXZZ','YYZZ','XXYZ','XYYZ','XYZZ'/

      DATA angu_pow_S /0,0,0/

      DATA angu_pow_P /1,0,0,
     &                 0,1,0,
     &                 0,0,1/

      DATA angu_pow_D /2,0,0,1,1,0,
     &                 0,2,0,1,0,1,
     &                 0,0,2,0,1,1/

      DATA angu_pow_F /3,0,0,2,2,1,0,1,0,1,
     &                 0,3,0,1,0,2,2,0,1,1,
     &                 0,0,3,0,1,0,1,2,2,1 /

C    ... newest...not defined in Gamess
      DATA angu_pow_G /
     & 4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1,               
     & 0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1,
     & 0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2
     &  /

#include "dummy.h"
      CALL QENTER('INPREAD')

C     ... zero array for number of contracted basis function
C     do i=1,200
      do i=1,NBASF
        Basis_N(i)=0
      End do

C     ... open file for reading integrals
      LUinp = 57
      OPEN(unit=LUinp,file='CAPINT.inp',STATUS='OLD',FORM='FORMATTED')
      READ(LUinp,*) n_atoms

      IF (IPRINT.GE.3) THEN
       write(LUPRI,'(/2X,A)') 'Output for INPREAD'
       write(lupri,*) ' n_atoms= ',n_atoms
        CALL FLSHFO(LUPRI)
      ENDIF

      g_n=1
      sh=0
      iiii=1
      DO  k=1,n_atoms
        READ(LUinp,*) IA

        IF (IPRINT.GE.3) THEN
          write(lupri,*) 'label  IA=',IA
        ENDIF

        READ(LUinp,*)(re(k,i),i=1,3)

        IF (IPRINT.GE.3) THEN
           write(lupri,*) '> coordinate: ',(re(k,i),i=1,3)
        ENDIF
          
        Do ISNEW=1,1000  

          READ(LUinp,*) shell_type, n_con_guass 
           
          IF (IPRINT.GE.3) THEN
           write(lupri,*) 'shell_type=',shell_type,
     &     ' n_con_guass=',n_con_guass 
          ENDIF

CMI       'E', means end
            if(shell_type.eq.BLNK(1))go to 201
CMI       'S' shell
            if(shell_type.eq.blnk(2)) go to 200
CMI       'P' shell
            if(shell_type.eq.blnk(3)) go to 300
CMI       'D' shell
            if(shell_type.eq.blnk(4)) go to 400
CMI       'F' shell
            if(shell_type.eq.blnk(5)) go to 500
CMI       'G' shell
            if(shell_type.eq.blnk(6)) go to 600
              
 200        continue
          
            DO n_c_g=1,n_con_guass

             sh=sh+1

             READ(LUinp,*)nj,expon,coe

          IF (IPRINT.GE.3) THEN
           write(lupri,'(2X,A,I2,A,F12.5,A,F12.6)')
     & 'S: nj=',nj,' expon=',expon,' contr.coef=',coe
          ENDIF

               ord_expon(sh)=expon
               ord_coe(sh)=coe

               shell_arr(sh)=shell_type
               Shell_basis(g_n)=shell_type
               Shell_atom(g_n)=k

               do ia_p=1,3
                  ang_pow(sh,ia_p)=angu_pow_S(1,ia_p)
               end do
               basis_n(g_n)=n_con_guass

            END DO
              g_n=g_n+1
              go to 6000
             
 300         continue
             
             DO n_c_g=1,n_con_guass
                sh=sh+1
                READ(LUinp,*) nj,expon,coe

          IF (IPRINT.GE.3) THEN
C          write(lupri,*)
           write(lupri,'(2X,A,I2,A,F12.5,A,F12.6)')
     &     'P: nj=',nj,' expon=',expon,' contr.coef=',coe
          ENDIF

             do i=1,3
                   shell_arr(sh+(i-1)*n_con_guass)=shell_p(i)
                   Shell_basis(g_n+i-1)=shell_p(i)
                   Shell_atom(g_n+i-1)=k

                   ord_expon(sh+(i-1)*n_con_guass)=expon
                   ord_coe(sh+(i-1)*n_con_guass)=coe
                   do ia_p=1,3
                   ang_pow(sh+(i-1)*n_con_guass,ia_p)=angu_pow_p(i,ia_p)
                   end do
             end do                              
           END DO
          
          do  iji=1,3
             basis_n(g_n)=basis_n(g_n)+n_con_guass
             g_n=g_n+1
          end do
          sh=sh+2*n_con_guass
          go to 6000
             
 400      continue
             
        DO n_c_g=1,n_con_guass
             sh=sh+1
          READ(LUinp,*)nj,expon,coe

          IF (IPRINT.GE.3) THEN
           write(lupri,'(2X,A,I2,A,F12.5,A,F12.6)')
     & 'D: nj=',nj,' expon=',expon,' contr.coef=',coe
          ENDIF

                do i=1,6
                   shell_arr(sh+(i-1)*n_con_guass)=shell_d(i)
                   Shell_basis(g_n+i-1)=shell_d(i)
                   Shell_atom(g_n+i-1)=k

                   ord_expon(sh+(i-1)*n_con_guass)=expon
                   ord_coe(sh+(i-1)*n_con_guass)=coe

                   do ia_p=1,3
                 ang_pow(sh+(i-1)*n_con_guass,ia_p)= angu_pow_d(i,ia_p)
                   end do

                end do
            END DO
       
             do  iji=1,6
                basis_n(g_n)=basis_n(g_n)+n_con_guass
                g_n=g_n+1
             end do
             sh=sh+5*n_con_guass
             go to 6000
             
CMI  *********** F *****************************
 500         continue 
             
             DO  n_c_g=1,n_con_guass
                sh=sh+1

                READ(LUinp,*)nj,expon,coe

           IF (IPRINT.GE.3) THEN
           write(lupri,'(2X,A,I2,A,F12.5,A,F12.6)')
     & 'F: nj=',nj,' expon=',expon,' contr.coef=',coe
           ENDIF

                do i=1,10
                   shell_arr(sh+(i-1)*n_con_guass)=shell_f(i)
                   Shell_basis(g_n+i-1)=shell_f(i)
                   Shell_atom(g_n+i-1)=k

                   ord_expon(sh+(i-1)*n_con_guass)=expon
                   ord_coe(sh+(i-1)*n_con_guass)=coe
                   do ia_p=1,3
                   ang_pow(sh+(i-1)*n_con_guass,ia_p)=angu_pow_f(i,ia_p)
                   end do
                end do
             END DO
             
             do  iji=1,10
                basis_n(g_n)=basis_n(g_n)+n_con_guass
                g_n=g_n+1
             end do
             sh=sh+9*n_con_guass
             go to 6000
CMI  *********** G *****************************
 600         continue 
      DO   n_c_g=1,n_con_guass
              sh=sh+1
              READ(LUinp,*)nj,expon,coe
              do i=1,15
                 shell_arr(sh+(i-1)*n_con_guass)=shell_G(i)
                 Shell_basis(g_n+i-1)=shell_G(i)
                 Shell_atom(g_n+i-1)=k
                 ord_expon(sh+(i-1)*n_con_guass)=expon
                 ord_coe(sh+(i-1)*n_con_guass)=coe
                  do ia_p=1,3
                   ang_pow(sh+(i-1)*n_con_guass,ia_p)=angu_pow_G(i,ia_p)
                  end do
               end do
      END DO
  
             do  iji=1,15
                basis_n(g_n)=basis_n(g_n)+n_con_guass
                 g_n=g_n+1
             end do
             sh=sh+(14)*n_con_guass
             go to 6000

 6000       continue
     
              End  Do

 201       continue
                  n_pg_atom(1)=1
                  n_pg_atom(k+1)=sh

               do ii= n_pg_atom(k),sh
                     do isw=1,3
                       R_ms(iiii,isw)=re(k,isw)
                     end do
                     iiii=iiii+1
                end do
       END DO

       close(LUinp,status="keep")

CMI   ... get the total number of basis functions
             n_b=g_n-1
       
       WRITE(LUPRI,'(/2X,A,I3,A,I3)')
     & 'Total number of primitive guassians / total number of'//
     & ' contracted basis functions=',sh,'/',n_b

       LUCAP=0
       CALL GPOPEN(LUCAP,
     &  'CAPLBL.inp',' ',' ','FORMATTED',IDUMMY,.FALSE.)

       WRITE(LUCAP,*) 'CAP unique order of atomic basis functions'
       WRITE(LUPRI,*) 'CAP unique order of atomic basis functions'

       do i=1,n_b
          write(LUPRI,'(2i4,A6,A,I3)')i,SHELL_ATOM(i),Shell_basis(i),
     &   ' contr->',basis_n(i)  
C        ... write into file
         WRITE(LUCAP,'(i4,2x,i3,2x,A6)') i,SHELL_ATOM(i),Shell_basis(i)
       end do

       CALL GPCLOSE(LUCAP,'KEEP')

       IF (IPRINT.GE.1) THEN
         WRITE(LUPRI,'(2X,A)')
     &  'File "CAPLBL.inp" with CAP unique order of atomic basis'//
     &  ' functions ready.'
       ENDIF
       
       WRITE(LUPRI,'(/,2X,A)')
     & '>>> sh_ty, Pgua_expo,nor_pg_coe,angu_pow, cartesian'
       
       DIRAC_NRM = .TRUE.
       IF (DIRAC_NRM) THEN
C!NORMALIZATION PART
C!Only the primitives are normalized
C!The Contracted guassian is kept unnormalized
C!This normalization is compatible with the
C!NORMF=1 NORMP=0 of gamess US
C!For a different normalization procedure,
C!one may have to modify this
          do i=1,sh
           do isw=1,3
             R_m(isw)=R_ms(i,isw)
       k_m(isw)=(ang_pow(i,isw))
           end do

           R_m(4)=ord_expon(i)

           CALL PEM(R_m,R_m,K_m,K_m,PEINT)
           ord_coe(i)=ord_coe(i)/sqrt(PEINT)

           write (LUPRI,*)i,'fu=',shell_arr(i),ord_expon(i),
     &             ord_coe(i),' nrm PEINT=',PEINT

C           write (16,*)i,shell_arr(i),ord_expon(i),ord_coe(i)
C           write(16,*) (ang_pow(i,j),j=1,3)
C           write(16,*) (R_ms(i,j),j=1,3)
         end do

       ELSE
CMI  ... OLD normalization of CAP - not used in DIRAC
C ==========================================================
       do i=1,sh
          do isw=1,3
             R_m(isw)=R_ms(i,isw)
             R_n(isw)=R_ms(i,isw)
             k_m(isw)=(ang_pow(i,isw))
             k_n(isw)=(ang_pow(i,isw))
          end do
          
         R_m(4)=ord_expon(i)
         R_n(4)=ord_expon(i)
          
         CALL PEM(R_m,R_n,K_m,K_n,PEINT)

         unnormcoe=ord_coe(i)
         ord_coe(i)=ord_coe(i)/sqrt(PEINT)

         write(LUPRI,*) i,'. bfu: ',shell_arr(i),
     &   ' exponent:',ord_expon(i)
         write(LUPRI,*) 
     &   ' contr.coeff:',unnormcoe,' normal.contr. coeff:',ord_coe(i)
         write(LUPRI,*) 'ang_pow:',(ang_pow(i,j),j=1,3)
         write(LUPRI,*) 'R_ms:',(R_ms(i,j),j=1,3)

       end do
       ENDIF

       CALL FLSHFO(LUPRI)

       CALL QEXIT('INPREAD')
       RETURN
       END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck aointegral */
      SUBROUTINE AOINTEGRAL(N_B,ord_expon,ord_coe,R_ms,basis_n,ang_pow,
     & PW,PE_AOr,PE_AOi,PEXR,PEXI,PEYR,PEYI,PEZR,PEZI,NBASF,NANGPOW) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Core routine for calculating CAP integrals.
C
C  Returns filled arrays PE_AOr,PE_AOi,PEXR,... (if DONRCAP=.false.)
C
C    Parameters:
C -------------------
C           N_B          number of atomic basis functions
C           ord_expon    exponents of primitives
C           ord_coe      contraction coefficients
C           R_ms
C           basis_n
C           ang_pow
C           PW - theta
C           PE_AOr,PE_AOi
C
C    Called from SAJEEV
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbcap.h"

      integer *4 n_b
C     integer *4 ang_pow(250,3)
      integer *4 ang_pow(NANGPOW,3)
      integer *4 k_m(3),k_n(3)
C     integer *4 basis_n(200)
      integer *4 basis_n(NBASF)
C     integer *4 basis_ord(300)
      integer *4 basis_ord(400)
     
C     real *8 ord_expon(250),ord_coe(250)
      real *8 ord_expon(NANGPOW),ord_coe(NANGPOW)
      real *8 R_m(4),R_n(4)
C     real *8 R_ms(200,3)
      real *8 R_ms(NBASF,3)

      COMPLEX *16 PEINT,PEINT_D1X,PEINT_D1Y,PEINT_D1Z
CMI   REAL*8 PE_AOr(25000),PE_AOi(25000)
CMI   REAL*8 PE_AOr(NNBBASX),PE_AOi(NNBBASX)
      REAL*8 PE_AOr(*),PE_AOi(*)
      REAL*8 PEXR(*),PEXI(*),PEYR(*),PEYI(*),PEZR(*),PEZI(*)
      COMPLEX *16 BETA

      CALL QENTER('AOINTEGRAL') 

C     ... zero array for integrals
      do i=1,n_b
         do j=1,n_b
              IPU=IPOINTIJ(i,j)
            IF (CAPOVLP) THEN
              PE_AOr(IPU)=0.d0
            ELSE 
             IF (DONRCAP) THEN
              PE_AOr(IPU)=0.d0
              PE_AOi(IPU)=0.d0
             ELSE 
              PEXR(IPU)=0.d0
              PEXI(IPU)=0.d0
              PEYR(IPU)=0.d0
              PEYI(IPU)=0.d0
              PEZR(IPU)=0.d0
              PEZI(IPU)=0.d0
             ENDIF
            ENDIF
         end do
      end do

C ... prepare offsets for contracted basis functions
      do i=1,n_b
         if(i.eq.1)then
           iasp=0
         else
           iasp=basis_ord(i-1)
         endif
         basis_ord(i)=iasp+basis_n(i)
         IF (IPRCAP.GE.3) THEN
          WRITE(LUPRI,*) 'b.fu: ',i,' contr.shells:',basis_n(i),
     &    ' offset:',basis_ord(i)
         ENDIF
      end do
C     ... run over contracted basis functions
      DO ii=1,n_b
      DO jj=1,ii
         IPU=IPOINTIJ(ii,jj)
         if (ii.eq.1) then
            iangii=0
         else
            iangii=basis_ord(ii-1)
         endif
         if (jj.eq.1)then
            iangjj=0
         else
            iangjj=basis_ord(jj-1)
         endif
C        ... run over primitives in contracted basis functions
         Do i=iangii+1,basis_ord(ii)
         do j=iangjj+1,basis_ord(jj)
            do isw=1,3
               R_m(isw)=R_ms(i,isw)
               R_n(isw)=R_ms(j,isw)
               k_m(isw)=(ang_pow(i,isw))
               k_n(isw)=(ang_pow(j,isw))
            end do
            R_m(4)=ord_expon(i)
            R_n(4)=ord_expon(j)

            IF (CAPOVLP) THEN
CMI         ... get overlap integrals
             CALL PEM(R_m,R_n,K_m,K_n,PEINT)
            ELSE
CMI         ... get CAP integrals    
             CALL CAPM(R_m,R_n,K_m,K_n,PW,PEINT,
     &          PEINT_D1X,PEINT_D1Y,PEINT_D1Z)
            ENDIF

C           ... contraction coefficients of i,j
            COMUL=ord_coe(i)*ord_coe(j)
C           ... save real and imaginary part of integrals
            IF (CAPOVLP) THEN
              PE_AOr(IPU)= PE_AOr(IPU)+COMUL*DREAL(PEINT)
            ELSE 
             IF (DONRCAP) THEN
              PE_AOr(IPU)= PE_AOr(IPU)+COMUL*DREAL(PEINT)
              PE_AOi(IPU)= PE_AOi(IPU)+COMUL*DIMAG(PEINT)
             ELSE
              PEXR(IPU)= PEXR(IPU)+COMUL*DREAL(PEINT_D1X)
              PEXI(IPU)= PEXI(IPU)+COMUL*DIMAG(PEINT_D1X)

              PEYR(IPU)= PEYR(IPU)+COMUL*DREAL(PEINT_D1Y)
              PEYI(IPU)= PEYI(IPU)+COMUL*DIMAG(PEINT_D1Y)

              PEZR(IPU)= PEZR(IPU)+COMUL*DREAL(PEINT_D1Z)
              PEZI(IPU)= PEZI(IPU)+COMUL*DIMAG(PEINT_D1Z)
             ENDIF
            ENDIF
         END DO 
         END DO

C        ... control print out of integrals
         IF (IPRCAP.GE.7) THEN
          IF (CAPOVLP) THEN
          WRITE(LUPRI,'(2X,A,I3,A,I3,A,F22.8)') 
     &    'CAP overlap(',ii,';',jj,')=',PE_AOr(IPU)
          ENDIF
         ENDIF

      END DO
      END DO

      CALL QEXIT('AOINTEGRAL') 
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION IFACT(l)
      IMPLICIT REAL*8 (A-H,O-Z)
      ifact=1
      if((l-1).lt.2)go to 200
      do j=1,(l)-1
         IF ((J/2)*2.EQ.j)GO TO 100
         IFact=ifact*j
 100     CONTINUE
      end do
      go to 300
 200  ifact=1
 300  continue
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE F_INTC(k,PA,PB,k_m,K_n,FINTM) 
      
      IMPLICIT REAL*8 (A-H,O-Z)
      integer *4 k_m,k_n
      real *8 Int
      EXTERNAL JFACT
      COMPLEX *16 PA,PB,FINTM
      COMPLEX *16 PAMUL,PBMUL,RADD

      FINTM=(0.d0,0.d0)
      do i=0,(K_m)
         do j=0,k_n
            if ((i+j).eq.k)then
               PAMUL= PA**(K_m-i)
               PBMUL=PB**(K_n-j)
               JFKM=JFACT(k_m)
               JFKN=JFACT(k_n)
               JFI=JFACT(i)
               JFJ=JFACT(j)
               JFIM=JFACT(K_m-i)
               JFJN=JFACT(K_n-j)
               RMFACT=JFKM/(JFI*JFIM)
               RNFACT=JFKN/(JFJ*JFJN)
               RADD=PAMUL*PBMUL*RMFACT*RNFACT
               FINTM=FINTM+RADD
            endif
         end do 
      end do
      RETURN
      END
      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION JFACT(l)
      IMPLICIT REAL*8 (A-H,O-Z)
      jfact=1
      if((l).lt.1)go to 200
      do j=1,l
         jFact=jfact*j
      end do
      go to 300
 200  jfact=1
 300  continue
      END  
      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION IPOINTIJ(ii,jj)
      IMPLICIT REAL*8 (A-H,O-Z)
      IPOINTIJ=(ii*(ii-1)/2)+jj
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PEM */
      SUBROUTINE PEM(R_m,R_n,K_m,K_n,PEINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  On input: 
C           R_m
C           R_n
C           K_m
C           K_n
C
C  On output: PEINT
C
C  Called from: INPREAD 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      real *8 R_CM(2),R_CN(2),R_KM(2),R_KN(2),R_LM(2),R_LN(2)
      real *8 R_m(4),R_n(4)
      integer *4 K_m(3),K_n(3)
      real *8 OVRLPI,OVRLPJ,OVRLPK,PEINT
      CALL QENTER('PEM')

         PEINT=0.d0
         i=1
         k=2
         L=3
         R_CM(1)=R_m(i)
         R_KM(1)=R_m(k)
         R_LM(1)=R_m(L)

         R_CN(1)=R_n(i)
         R_KN(1)=R_n(k)
         R_LN(1)=R_n(L)

         R_CM(2)=R_m(4)
         R_KM(2)=R_m(4)
         R_LM(2)=R_m(4)

         R_CN(2)=R_n(4) 
         R_KN(2)=R_n(4)
         R_LN(2)=R_n(4)
  
         K_CM=K_m(i)
         K_KM=K_m(k)
         K_LM=K_m(L)

         K_CN=K_n(i)
         k_KN=K_n(k)     
         K_LN=k_n(L)

         CALL SMAT(R_KM,R_KN,K_KM,K_KN,0.d0,OVRLPI)

         CALL SMAT(R_LM,R_LN,K_LM,K_LN,0.d0,OVRLPJ)

         CALL SMAT(R_CM,R_CN,K_CM,K_CN,0.d0,OVRLPK)

         PEINT=PEINT+(OVRLPI*OVRLPJ*OVRLPK)

      CALL QEXIT('PEM')
      RETURN
      END
 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck SMAT */
      subroutine SMAT(R_m,R_n,K_m,K_n,BETA,Int_GR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Parameters:
C            R_m,R_n
C            K_m,K_n
C            BETA (usually zero)
C
C            Int_GR - returns overlap integrals
C
C  Called from: PEM
C               CAPM
C
C
C  Output:  Int_GR
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      integer *4 k_m,k_n
      real *8  R_m(2),R_n(2)
      real *8  BETA, INT_GR
      real *8  GAM,P,PA,PB,GAMP,RMULT,RMULTI,FINTM
      data PI/3.14/
      EXTERNAL  IFACT
   
      GAM=R_m(2)+(R_n(2)+BETA)
      P=(R_m(1)*R_m(2)+R_n(1)*(R_n(2)+BETA))/GAM
      PA=P-R_m(1)
      PB=P-R_n(1)
      AB=R_m(1)-R_n(1) 

      DO i=0,(K_m+k_n)/2
C         if (i.eq.1)write(*,*)i
         INT_GR=0.d0
         if (i.eq.0)go to 150
C         if ((i/2)*2.ne.i)go to 100
     
 150     continue
         IFI=IFACT(2*i)
C         if (i.eq.4)write(*,*)IFACT(i),'IFA'
         GAMP=(2*gam)**(-i)
         PIG=(PI/GAM)**0.5
         RMULT=IFI*PIG*gamp
         CALL F_INT(2*i,PA,PB,k_m,K_n,FINTM)
         RMULTI=RMULT*FINTM 
         INT_GR=INT_GR+RMULTI 
 100     continue
      END DO
      
      PREMULT=exp((-1.d0*R_m(2)*(R_n(2)+BETA)*((AB)**2))/gam)
      if (dABS(premult).lt.0.000000000d0)PREMULT=0.000000000d0
      INT_GR=INT_GR*PREMULT
      RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE F_INT(k,PA,PB,k_m,K_n,FINTM) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Called from SMAT
C
C
C   Returns FINTM
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      
CMI   IMPLICIT REAL*8 (A-H,O-Z)
      integer *4 k_m,k_n
      real *8 Int
C     REAL *8 PA,PB
      EXTERNAL JFACT
      real *8 PA,PB,FINTM
      real *8 PAMUL,PBMUL,RADD

      FINTM=(0.d0,0.d0)
      do i=0,(K_m)
         do j=0,k_n
            if ((i+j).eq.k)then
               PAMUL= PA**(K_m-i)
               PBMUL=PB**(K_n-j)
               JFKM=JFACT(k_m)
               JFKN=JFACT(k_n)
               JFI=JFACT(i)
               JFJ=JFACT(j)
               JFIM=JFACT(K_m-i)
               JFJN=JFACT(K_n-j)
               RMFACT=JFKM/(JFI*JFIM)
               RNFACT=JFKN/(JFJ*JFJN)
               RADD=PAMUL*PBMUL*RMFACT*RNFACT
               FINTM=FINTM+RADD
            endif
         end do 
      end do
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck capm */
      SUBROUTINE CAPM(R_m,R_n,K_m,K_n,theta,
     &         VCAP,VCAP_D1X,VCAP_D1Y,VCAP_D1Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Called from AOINTEGRAL (from within the 4-fold cycle)
C
C   entering parameter "theta" is from the INPQD file
C
C   returns: 
C          VCAP - the total nonrel. scalar CAP term for DONRCAP=.true.
C          VCAP_D1X,VCAP_D1Y,VCAP_D1Z relativ. terms  (DONRCAP=.false.)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbcap.h"
#include "capint.h"

      integer *4 NSTEP
C     integer *4 NGR
      
      real *8 R_CM(2),R_CN(2),R_KM(2),R_KN(2),R_LM(2),R_LN(2)
      real *8 R_m(4),R_n(4)
      integer *4 K_m(3),K_n(3)
      real *8 CAP_PAR(3)
CMI   REAL *8 PGRT(100),GINTPT(100)
      REAL *8 PGRT(100)
CMI   real *8 X0,theta,LAMDA,GRSTP(2)
      real *8 theta,GRSTP(2)
CMI   integer *4 NGRS(100)
      complex  *16 PE,PE1
      complex *16 VCAP1,VCAP,VCAPS(3),VCAP1C
      complex *16 VCAP_D1X,VCAP_D1Y,VCAP_D1Z,TEMP
      COMMON PGRT

      CALL QENTER('CAPM')
CMI   .... no need for reading this file !!!
C     open(unit=8,file='INPCAP',STATUS='UNKNOWN',FORM='FORMATTED')
C     REWIND (8)
     
      PE=(0.d0,0.d0)
      VCAP_D1X = (0.d0,0.d0)
      VCAP_D1Y = (0.d0,0.d0)
      VCAP_D1Z = (0.d0,0.d0)
C===================================================================
C         Run over x/y/z... i = 1/2/3 
C===================================================================
      DO i=1,3
CMI      ... replaced by X0(1:3), LAMDA(1:3) in capint.h
CMI      READ(8,*)X0,LAMDA
CMI      .... replaced by NGR(1:3) in capint.h
CMI      READ(8,*)NGR
CMI      .... replaced by GINTPT(3,MAXSPL) in capint.h
CMI      READ(8,*)(GINTPT(j),j=1,NGR)
CMI       .... replaced by NGRS(3,MAXSPL) in capint.h
CMI      READ(8,*)(NGRS(j),j=1,NGR-1)
   
      CAP_PAR(3)=THETA
CMI   CAP_PAR(2)=LAMDA
      CAP_PAR(2)=LAMDA(I)
CMI   CAP_PAR(1)=X0
      CAP_PAR(1)=X0(I)
      
         if(i.eq.1)then
            k=2
            L=3
         endif
         if(i.eq.2)then
            k=1
            L=3
         endif
         if(i.eq.3)then
            k=1
            L=2
         endif

         R_CM(1)=R_m(i)
         R_KM(1)=R_m(k)
         R_LM(1)=R_m(L)

         R_CN(1)=R_n(i)
         R_KN(1)=R_n(k)
         R_LN(1)=R_n(L)

         R_CM(2)=R_m(4)
         R_KM(2)=R_m(4)
         R_LM(2)=R_m(4)

         R_CN(2)=R_n(4) 
         R_KN(2)=R_n(4)
         R_LN(2)=R_n(4)

         K_CM=K_m(i)
         K_KM=K_m(k)
         K_LM=K_m(L)

         K_CN=K_n(i)
         k_KN=K_n(k)     
         K_LN=k_n(L)
       
CMI     ... initialize to zero !
         VCAPS(i)=(0.d0,0.d0)

CMI      ... run over various integration intervals
CMI      do j=1,NGR-1
         do j=1,NGR(i)-1
CMI         GRSTP(1)=GINTPT(j)
            GRSTP(1)=GINTPT(i,j)
CMI         GRSTP(2)=GINTPT(j+1)
            GRSTP(2)=GINTPT(i,j+1)
CMI         NSTEP=NGRS(j)
            NSTEP=NGRS(i,j)
CMI         write(LUPRI,*) 'i,j,NGRS(i,j):',i,j,NGRS(i,j)
            CALL CAPMAT1(CAP_PAR,GRSTP,NSTEP,R_CM,R_CN,K_CM,K_CN,
     &                    VCAP1C,VCAP1)

C           ... accumulate for each interval
            IF (DONRCAP) THEN
              VCAPS(i)=VCAPS(i)+VCAP1
            ELSE
              IF (I.EQ.1) THEN
                VCAP_D1X = VCAP_D1X + VCAP1C 
              ELSE IF (I.EQ.2) THEN
                VCAP_D1Y = VCAP_D1Y + VCAP1C
              ELSE IF (I.EQ.3) THEN
                VCAP_D1Z = VCAP_D1Z + VCAP1C
              ELSE
                write(LUPRI,*) 'I=',I
                CALL QUIT('CAPM: wrong I!')
              ENDIF
            ENDIF

         end do

         CALL SMAT(R_KM,R_KN,K_KM,K_KN,0.d0,OVRLPK)
         CALL SMAT(R_LM,R_LN,K_LM,K_LN,0.d0,OVRLPL)

         IF (DONRCAP) THEN
           VCAPS(i)=VCAPS(i)*OVRLPK*OVRLPL
         ELSE
           IF (I.EQ.1) THEN
              VCAP_D1X = VCAP_D1X*OVRLPK*OVRLPL
           ELSE IF (I.EQ.2) THEN
              VCAP_D1Y = VCAP_D1Y*OVRLPK*OVRLPL
           ELSE IF (I.EQ.3) THEN
              VCAP_D1Z = VCAP_D1Z*OVRLPK*OVRLPL
           ELSE
              write(LUPRI,*) 'I=',I
              CALL QUIT('CAPM 2: wrong I!')
           ENDIF
         ENDIF

      ENDDO ! of I=1,3
             
C   ... finally add all theree terms - for x,y and z - for DONRCAP=.true.!!!
      VCAP=VCAPS(1)+VCAPS(2)+VCAPS(3)
C ------------------------------------------------------------

      IF (IPRCAP.GE.19) THEN
       IF (.NOT.DONRCAP) THEN
        write(LUPRI,*) '** CAPM >>  VCAP_D1X,VCAP_D1Y,VCAP_D1Z:',
     &      VCAP_D1X,VCAP_D1Y,VCAP_D1Z   
       ELSE
        write(LUPRI,*) '** CAPM >>  VCAP=',VCAP
       ENDIF
      ENDIF

      CALL QEXIT('CAPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck capm */
CMI   SUBROUTINE CAPMAT1(CAP_PAR,GRSTP,NSTEP,R_m,R_n,K_m,K_n,VCAPT)
      SUBROUTINE CAPMAT1(CAP_PAR,GRSTP,NSTEP,R_m,R_n,K_m,K_n,
     &                   VCAP1C,VCAPT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  returns: 
C
C      VCAPT (either for x,y or z) when non-relativistic form of CAP is on
C      or  both VCAP1C, VCAPT (for x/y/z) when Dirac CAP is on
C
C  Called from CAPM
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbcap.h"
C#include "capint.h"

      integer *4 k_m,k_n
      integer *4 NSTEP
      real *8 R_m(2),R_n(2)
      real *8 CAP_PAR(3)
      real *8 M
      real *8 X0,theta,LAMDA,GRSTP(2)
      complex *16 FFx,fx,Dfx,D2fx
      complex *16 V0x,V1x,V2x
      complex *16 V0xy,V1xy,V2xy

      complex *16  TempV1x,TempV1xy

      complex *16 CI,ROT
      complex *16 PEM,XS,PEN

      COMPLEX *16 VCAPT,VCAPXT,VCAPYT,VCAPZT
      COMPLEX *16 VCAPX(4),VCAPX_D1C(4),VCAP1C
      COMPLEX *16 VCAP,VADD,VF
      COMPLEX *16 VCAP_D1C,VCAPy_D1C
      COMPLEX *16 VCAPTy,VCAPXTy,VCAPYTy,VCAPZTy
      COMPLEX *16 VCAPXy(4),VCAPXy_D1C(4)
      COMPLEX *16 VCAPy,VADDy,VFy
      COMPLEX *16 VCAPXT_D1C,VCAPXTy_D1C
      EXTERNAL V0,V1,V2,FF
      EXTERNAL V0D,V1D
CMI   ... declaration must be for gfortran
      double complex v0,v1,v2,V1D,FF,V0D
      
      data CI/(0.0d0,1.0d0)/
      data M/1.0d00/

      CALL QENTER('CAPMAT1')

      STEP=(GRSTP(2)-GRSTP(1))/(NSTEP)
      X0=CAP_PAR(1)
      LAMDA=CAP_PAR(2)
      THETA=CAP_PAR(3)
    
CMI   ... initialization of summation variables
      VCAPXT=(0.d0,0.d0)
      VCAPXTy=(0.d0,0.d0)

      VCAPXT_D1C=(0.d0,0.d0)
      VCAPXTy_D1C=(0.d0,0.d0)

CMI   ... starting position for the integration
      XI=GRSTP(1)

CMI  ... run over the integration interval with NSTEP nodes
      itx=0
      DO ix=1,NSTEP
         X=(XI+(ix-1)*STEP)
         Xy=X*(-1.0d0)

         itx=itx+1

         XR=X-R_n(1)
         XRM=X-R_m(1)
         DGr=K_n*(XR**(-1))

C         ... 1st derivative
         DG=DGr-(2*R_n(2)*(XR))
       
         XRy=Xy-R_n(1)
         XRMy=Xy-R_m(1)
         DGry=K_n*(XRy**(-1))

C         ... 1st derivative
         DGy=DGry-(2*R_n(2)*(XRy))

         D2G1X=K_n*(K_n-1)*(XR**(-2))
         D2G2X=2*R_n(2)*(2*K_n+1)
         D2G3X=4*(R_n(2)**2)*(XR**2)
         D2GX=D2G1X-D2G2X+D2G3X

C        ... 2nd derivative
         D2G=D2GX

         D2G1Xy=K_n*(K_n-1)*(XRy**(-2))
         D2G2Xy=2*R_n(2)*(2*K_n+1)
         D2G3Xy=4*(R_n(2)**2)*(XRy**2)
         D2GXy=D2G1Xy-D2G2Xy+D2G3Xy

C        ... 2nd derivative
         D2Gy=D2GXy

         GG1M=(XRM**2)
         GG1N=(XR**2)
         GG1=R_m(2)*GG1M+R_n(2)*GG1N

         GG1My=(XRMy**2)
         GG1Ny=(XRy**2)
         GG1y=R_m(2)*GG1My+R_n(2)*GG1Ny

         GGEX=(XR**K_n)*(XRM**K_m)
         GG=GGEX*DEXP(-GG1)

         GGEXy=(XRy**K_n)*(XRMy**K_m)
         GGy=GGEXy*DEXP(-GG1y)

C        VADD=DCMPLX(0.0d0,0.0d0)
         VADD=0.0d0
C        VADDy=DCMPLX(0.0d0,0.0d0)
         VADDy=0.0d0
       
         IF (DONRCAP) THEN
CMI      ... nonrelativistic CAP operator
           v0x=V0(M,THETA,LAMDA,X0,X)
           V1x=v1(M,THETA,LAMDA,X0,X)*DG
           v2x=v2(M,THETA,LAMDA,X0,X)*D2G
           v0xy=V0(M,THETA,LAMDA,X0,Xy)
           V1xy=v1(M,THETA,LAMDA,X0,Xy)*DGy
           v2xy=v2(M,THETA,LAMDA,X0,Xy)*D2Gy
           VCAP=(V0x+V1x+V2x+VADD)*GG
           VCAPy=(V0xy+V1xy+V2xy+VADDy)*GGy
           VCAPX(itx)=VCAP
           VCAPXy(itx)=VCAPy
           if(ix.eq.1) go to 321 
           VCAPXT=VCAPXT+0.5*(VCAPX(2)+VCAPX(1))*STEP
           VCAPXTy=VCAPXTy+0.5*(VCAPXy(2)+VCAPXy(1))*STEP
           if (ix.ne.1) VCAPX(2)=VCAPX(1)
           if (ix.ne.1) VCAPXy(2)=VCAPXy(1)
         ELSE
CMI       ... relativistic Hamiltonian
          IF (.NOT.CAPDIPVEL) THEN
           V0x  = V0D(M,THETA,LAMDA,X0,X )
           V0xy = V0D(M,THETA,LAMDA,X0,Xy)
           V1x  = V1D(M,THETA,LAMDA,X0,X )*DG
C          V1x  = V1D(M,THETA,LAMDA,X0,X )*DG*DG
           V1xy = V1D(M,THETA,LAMDA,X0,Xy)*DGy
C          V1xy = V1D(M,THETA,LAMDA,X0,Xy)*DGy*DGy
           VCAP_D1C  = (V0x +V1x +VADD )*GG
           VCAPy_D1C = (V0xy+V1xy+VADDy)*GGy
          ELSE
CMI       ... CAPDIPVEL = .TRUE.
C          V1x  =  DCMPLX(1.0D0,0.D0)*DG
           V1x  =  DG
C          V1xy =  DCMPLX(1.0D0,0.D0)*DGy
           V1xy =  DGy
           VCAP_D1C  = (V1x +VADD )*GG
C          VCAP_D1C  = (V1x +VADD )
           VCAPy_D1C = (V1xy+VADDy)*GGy
C          VCAPy_D1C = (V1xy+VADDy)
          ENDIF
           VCAPX_D1C(itx) = VCAP_D1C
           VCAPXy_D1C(itx)= VCAPy_D1C
           if(ix.eq.1) go to 321 
C         ... numerical (trapezoidal) integration 
           VCAPXT_D1C=VCAPXT_D1C  +
     &         0.50D0*(VCAPX_D1C(2)+VCAPX_D1C(1))*STEP
           VCAPXTy_D1C=VCAPXTy_D1C+
     &         0.50D0*(VCAPXy_D1C(2)+VCAPXy_D1C(1))*STEP
           if (ix.ne.1) then
             VCAPX_D1C(2)=VCAPX_D1C(1)
             VCAPXy_D1C(2)=VCAPXy_D1C(1)
           endif 
         ENDIF

         itx=0
C         if((zabs(vcap)).lt.(0.00000000001))exit

 321     CONTINUE

      END DO

      IF (DONRCAP) THEN
CMI    ... nonrelativistic CAP
        VCAPT=VCAPXT+VCAPXTy
CMI     write(LUPRI,*) 'VCAPXT,VCAPXTy:',VCAPXT,VCAPXTy
      ELSE
C      ... relativistic CAP
CMI    ... for V1D_C (x,y,z) term
        VCAP1C=VCAPXT_D1C+VCAPXTy_D1C
CMI     write(LUPRI,*) '*** CAPMAT1 > VCAP1C:',VCAP1C
      ENDIF

      CALL QEXIT('CAPMAT1')
      RETURN
      END

C************************************************************************
C     Implementation of basic functions from the paper.
C************************************************************************
C       ...big F(x);  Nimrod:1998, eq. 18
      DOUBLE COMPLEX FUNCTION FF(THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 THETA,LAMBDA,X0,X
      DOUBLE PRECISION AUX
      AUX = DCOSH(LAMBDA*(X-X0))/DCOSH(LAMBDA*(X+X0))
      AUX = X + DLOG(AUX)/(2.0D0*LAMBDA)
      FF = CDEXP(+(0.0D0,1.0D0)*DCMPLX(THETA))-(1.0D0,0.0D0)
      FF = DCMPLX(X) + FF*DCMPLX(AUX)
      END

C     ... small f(x);  eq.4 of Nimrod:1998 with g(x) defined in eq.17
      DOUBLE COMPLEX FUNCTION F0(THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 THETA,LAMBDA,X0,X
      DOUBLE PRECISION AUX
      AUX = DTANH(LAMBDA*(X-X0)) - DTANH(LAMBDA*(X+X0))
      AUX = 1.0D0 + AUX/2.0D0
      F0 = CDEXP(+(0.0D0,1.0D0)*DCMPLX(THETA))-(1.0D0,0.0D0)
      F0 = (1.0D0,0.0D0) + F0*DCMPLX(AUX)
      END

C        ... first derivative of small f(x), eq. 20 of Nimrod:1998
      DOUBLE COMPLEX FUNCTION F1(THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8  THETA,LAMBDA,X0,X
      DOUBLE PRECISION AUX
      AUX = (DTANH(LAMBDA*(X+X0))**2) - (DTANH(LAMBDA*(X-X0))**2)
      F1 = CDEXP(+(0.0D0,1.0D0)*DCMPLX(THETA))-(1.0D0,0.0D0)
      F1 = DCMPLX(LAMBDA)*F1*DCMPLX(AUX)/(2.0D0,0.0D0)
      END

C     ...second derivative of small f(x), eq.21 of Nimrod:1998
      DOUBLE COMPLEX FUNCTION F2(THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 THETA,LAMBDA,X0,X
      DOUBLE PRECISION AUX1,AUX2
      AUX1 = DTANH(LAMBDA*(X+X0))*(1.0D0-(DTANH(LAMBDA*(X+X0))**2))
      AUX2 = DTANH(LAMBDA*(X-X0))*(1.0D0-(DTANH(LAMBDA*(X-X0))**2))
      F2 = CDEXP(+(0.0D0,1.0D0)*DCMPLX(THETA))-(1.0D0,0.0D0)
      F2 = DCMPLX(LAMBDA**2)*F2*(DCMPLX(AUX1)-DCMPLX(AUX2)) 
      END

CMI     ... Eq. 15 of Sajeev:2006 and Eq. 16 of Nimrod:1998
      DOUBLE COMPLEX FUNCTION V0(M,THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8  M,THETA,LAMBDA,X0,X
      DOUBLE COMPLEX CA1,CA2
      EXTERNAL F0,F1,F2
      DOUBLE COMPLEX F0,F1,F2
      CA1 = F2(THETA,LAMBDA,X0,X)/(F0(THETA,LAMBDA,X0,X)**3)
      CA1 = CA1/(DCMPLX(4.0D0*M))
      CA2 = (F1(THETA,LAMBDA,X0,X)**2)/(F0(THETA,LAMBDA,X0,X)**4)
CMI   CA2 = CA2/DCMPLX(5.0D0*M/8.0D0)
      CA2 = CA2/DCMPLX(8.0D0*M/5.0D0)
      V0 = CA1 - CA2
      END

CMI   ...  new derived V0 for the DIRAC Hamiltonian
      DOUBLE COMPLEX FUNCTION V0D(M,THETA,LAMBDA,X0,X)
      IMPLICIT NONE
#include "priunit.h"
      real *8  M,THETA,LAMBDA,X0,X
      DOUBLE COMPLEX CA1,CA2
      EXTERNAL F0,F1
      DOUBLE COMPLEX F0,F1
      CA1 = -0.50D0*(F1(THETA,LAMBDA,X0,X)/(F0(THETA,LAMBDA,X0,X)**2))
C     CA1 = -0.50D0*( (F0(THETA,LAMBDA,X0,X)**(-2)) )
CMI   write(LUPRI,*) 'V0D function:',CA1
      V0D = CA1 
      END

CMI   .... Eq. 16 of Sajeev:2006
      DOUBLE COMPLEX FUNCTION V1(M,THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 M,THETA,LAMBDA,X0,X
      EXTERNAL F0,F1
      DOUBLE COMPLEX F0,F1
      V1 = F1(THETA,LAMBDA,X0,X)/(F0(THETA,LAMBDA,X0,X)**3)
      V1 = V1/DCMPLX(M)
      END

CMI  ... V1 for the DIRAC Hamiltonian
      DOUBLE COMPLEX FUNCTION V1D(M,THETA,LAMBDA,X0,X)
      IMPLICIT NONE
#include "priunit.h"
      real *8 M,THETA,LAMBDA,X0,X
      EXTERNAL F0
      DOUBLE COMPLEX F0,CA
      CA  = ( DCMPLX(M)/F0(THETA,LAMBDA,X0,X) ) - DCMPLX(M)
C     CA  = ( 1.0D0/F0(THETA,LAMBDA,X0,X) ) - 1.0D0
C     CA  = ( DCMPLX(M) - CA1 ) / CA1
CMI    write(LUPRI,*) 'V1D function:',CA
C      write(LUPRI,*) 'V1D function:',CA, ' and F1:',CA1,
C    &   'theta,la,xo=',THETA,LAMBDA,X0
      V1D = CA
CMI   V1D = 2.0D0*CA
      END

CMI  ... Eq. 17 of Sajeev:2006
      DOUBLE COMPLEX FUNCTION V2(M,THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 M,THETA,LAMBDA,X0,X
      EXTERNAL F0
      DOUBLE COMPLEX F0
      V2 = (1.0D0,0.0D0)/(F0(THETA,LAMBDA,X0,X)**2)
      V2 = ((1.0D0,0.0D0)-V2)/DCMPLX(2.0D0*M)
      END
