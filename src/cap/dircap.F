!dirac_copyright_start
!      Copyright (c) 2013 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sajeev */
      SUBROUTINE SAJEEV(
     &   basis_n,ang_pow,basis_ord,ord_expon,
     &   ord_coe,r_ms,PE_AOr,PE_AOi,
     &   PEXR,PEXI,PEYR,PEYI,PEZR,PEZI,
C    &   NBASF,NANGPOW)
     &   WORK,LWORK,NCBF,NUBF)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Main program from Sajeev. Adapted into DIRAC enviroment
C  as the subroutine.
C
C On input: 
C              basis_n
C              ang_pow
C              basis_ord
C              ord_expon
C              ord_coe ...
C
C              PE_AOr, PE_AOi, PEXR... allocated array
C
C              NCBF - # of contracted basis functions
C              NUBF - # of uncontracted (primitives) basis functions
C
C On output: 
C            complex CAP integrals written into formatted file, or test
C            calculation of overlap integrals only for comparison with
C            DIRAC
C
C  called from DOSAJEEV/prepcap.F
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
      integer *4 basis_n(NCBF),basis_ord(NCBF)

C     integer *4 ang_pow(250,3)
      integer *4 ang_pow(NUBF,3)

      integer *4  N_B   

C     real *8 ord_expon(250),ord_coe(250)
      real *8 ord_expon(NUBF), ord_coe(NUBF)

      real *8 R_m(4),R_n(4)

C     real *8 R_ms(200,3)
      real *8 R_ms(NUBF,3)

C     COMPLEX *16 PE_AO(25000) ! former variable split into real and
C     imaginary parts
CMI   REAL*8 PE_AOr(25000),PE_AOi(25000)
CMI   REAL*8 PE_AOr(NNBBASX),PE_AOi(NNBBASX)
      REAL*8 PE_AOr(*),PE_AOi(*)
      REAL*8 PEXR(*),PEXI(*),PEYR(*),PEYI(*),PEZR(*),PEZI(*)
      real *8 BETA0,PE0
      
      LOGICAL ISFILE, FNDLAB,CHECK_PEMAT

      COMPLEX *16 THETA,THETA2,CI,BETA

      CHARACTER*10 PEMAT_NAME
      CHARACTER*8 LBL

      CALL QENTER('SAJEEV')

      IF (IPRCAP.GE.2) THEN
        CALL HEADER('*** Output from SAJEEV ***',-1)
        write(LUPRI,'(2X,A,I3,A,I3)')
     &   '>> contracted NCBF=',NCBF,' uncontracted NUBF=',NUBF
        write(LUPRI,*) 'LWORK=',LWORK
        CALL FLSHFO(LUPRI)
      ENDIF

#include "memint.h"
      CALL MEMGET('REAL',KTEST,1,WORK,KFREE,LFREE)
      CALL MEMCHK('sajeev 1',WORK,1)
  
      CALL INPREAD(N_B,ord_expon,ord_coe,R_ms,basis_n,basis_ord,
     &         ang_pow,CAPOVLP,IPRCAP,NCBF,NUBF,WORK(KFREE),LFREE)

      CALL MEMCHK('sajeev 2',WORK,1)

      IF (.NOT.CHECK_PEMAT()) THEN

CMI   ... evaluate CAP complex atomic untegrals
       IF (IPRCAP.GE.0) THEN
         WRITE(LUPRI,'(/2X,A)')
     &   'CAP integrals are not present, going'//
     &   ' to calculate and write them into the PEMAT file.'
         CALL FLSHFO(LUPRI)
       ENDIF

CMI   ... evaluate integrals - fill array with them
       CALL AOINTEGRAL(N_B,ord_expon,ord_coe,R_ms,
     &      basis_ord,ang_pow,
     &      PW,PE_AOr,PE_AOi,PEXR,PEXI,PEYR,PEYI,PEZR,PEZI,
     &      NCBF,NUBF)

CMI   ... write integrals into the PEMAT file ...
       CALL  PEMAT_WRITE(PE_AOr,PE_AOi,
     &  PEXR,PEXI,PEYR,PEYI,PEZR,PEZI)
        
      ELSE
       WRITE(LUPRI,'(/2X,A)')
     & 'SAJEEV: PEMAT file with desired CAP integrals is'//
     & ' present, not calculating them.'
      ENDIF

      CALL QEXIT('SAJEEV')
      RETURN
      END

      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck inpread */
      SUBROUTINE INPREAD(
     &  N_B,ord_expon,ord_coe,R_ms,basis_n,basis_ord,ang_pow,
     &  CAPOVLP,IPRINT,NCBF,NUBF,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Read the GAMESS-like basis set input file 
C
C  Parameters:
C
C     CAPOVLP - calculate only overlap integrals (for comparison
C               with DIRAC)
C     
C   fills these arrays:
C -----------------------
C     ord_expon
C     ord_coe
C     basis_n
C     basis_ord
C     R_ms
C     ang_pow
C
C
C  Called from: SAJEEV
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
      PARAMETER(MSIZE=600)

      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     &           DP25 = 0.25D0, DP5 = 0.5D0, DP75 = 0.75D0,
     &           THRMIN = 1.D-17)

      REAL*8 WORK(LWORK)
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

!   ... powers on x,y,z = up to "g" 
      integer *4 angu_pow_S(1,3),angu_pow_P(3,3)
      integer *4 angu_pow_D(6,3),angu_pow_F(10,3),angu_pow_G(15,3)

      integer *4 n_b,ia_p

C     integer *4 ang_pow(250,3)
      integer *4 ang_pow(NUBF,3)
      integer *4 n_pg_atom(50)
      integer *4 k_m(3),k_n(3)

C     integer *4 basis_n(200)
      integer *4 basis_n(NCBF)
      
      integer *4 basis_ord(NCBF)

      real *8 re(50,3),expon,coe
C     real *8 ord_expon(250),ord_coe(250)
      real *8 ord_expon(NUBF),ord_coe(NUBF)

      real *8 R_m(4),R_n(4)
C     real *8 R_ms(200,3)
      real *8 R_ms(NUBF,3)

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

      IF (NCBF.GT.MSIZE) THEN
        WRITE(LUPRI,'(A,I4,A,I4)') 'INPREAD: number of basis fu:',NCBF,
     &                 ' while allocated space is MSIZE ',MSIZE
        CALL QUIT('INPREAD: increase parameter MSIZE !!!')
      ENDIF

C     ... zero array for number of contracted basis function
C     do i=1,200
      do i=1,NCBF
        Basis_N(i)=0
      End do

C     ... open file for reading integrals
      LUinp = 57
      OPEN(unit=LUinp,file='CAPINT.inp',STATUS='OLD',FORM='FORMATTED')
      READ(LUinp,*) n_atoms

      IF (IPRINT.GE.2) THEN
       CALL HEADER('*****  Output for INPREAD  *****',-1 )
       WRITE(LUPRI,'(2X,A,I3)') 'Number of atoms :',n_atoms
       CALL FLSHFO(LUPRI)
      ENDIF

      g_n=1
      sh=0
CMI ================================================
CMI       ... cycle over number of atoms ...
CMI ================================================
      DO k=1,n_atoms
        READ(LUinp,*) IA

        IF (IPRINT.GE.3) THEN
          write(lupri,*) '....label  IA=',IA
        ENDIF

CMI  ... next line - coordinates
        READ(LUinp,*)(re(k,i),i=1,3)

        IF (IPRINT.GE.3) THEN
           write(lupri,'(2X,A,3F13.7)')
     &    'nucleus coordinate: ',(re(k,i),i=1,3)
           CALL FLSHFO(LUPRI)
        ENDIF
          
        Do ISNEW=1,1000  

          READ(LUinp,*) shell_type, n_con_guass 
           
          IF (IPRINT.GE.3) THEN
           write(lupri,'(2X,A,A2,A,I3)') 'shell_type:',shell_type,
     &     ' # of contr. gauss:',n_con_guass 
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
           CALL FLSHFO(LUPRI)
          ENDIF

               ord_expon(sh)=expon
               ord_coe(sh)=coe

               shell_arr(sh)=shell_type
               Shell_basis(g_n)=shell_type
               Shell_atom(g_n)=k

               do ia_p=1,3
                  ang_pow(sh,ia_p)=angu_pow_S(1,ia_p)
               end do

               if (g_n.gt.NCBF) then
                write(LUPRI,*) 'INPREAD: g_n,NCBF=',g_n,NCBF
                call QUIT('INPREAD: g_n > NBASF ! Increase NBASF')
               endif

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
           CALL FLSHFO(LUPRI)
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
             if (g_n.gt.NCBF) then
               write(LUPRI,*) 'INPREAD: g_n,NBASF=',g_n,NBASF
               call QUIT('INPREAD: g_n > NBASF ! Increase NBASF')
             endif
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
     &     'D: nj=',nj,' expon=',expon,' contr.coef=',coe
           CALL FLSHFO(LUPRI)
          ENDIF

          do i=1,6
              shell_arr(sh+(i-1)*n_con_guass)=shell_d(i)
              Shell_basis(g_n+i-1)=shell_d(i)
              Shell_atom(g_n+i-1)=k
              ord_expon(sh+(i-1)*n_con_guass)=expon
              ord_coe(sh+(i-1)*n_con_guass)=coe
              do ia_p=1,3
                 ang_pow(sh+(i-1)*n_con_guass,ia_p)= angu_pow_d(i,ia_p)
              enddo
          enddo
        ENDDO
       
             do  iji=1,6
               if (g_n.gt.NCBF) then
                 write(LUPRI,*) 'INPREAD: g_n,NBASF=',g_n,NBASF
                 call QUIT('INPREAD: g_n > NBASF ! Increase NBASF')
               endif
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
                 enddo
              enddo
             ENDDO
             
             do  iji=1,10
               if (g_n.gt.NCBF) then
                 write(LUPRI,*) 'INPREAD: g_n,NBASF=',g_n,NCBF
                 call QUIT('INPREAD: g_n > NBASF ! Increase NBASF')
               endif
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
  
             do iji=1,15
               if (g_n.gt.NCBF) then
                 write(LUPRI,*) 'INPREAD: g_n,NCBF=',g_n,NCBF
                 call QUIT('INPREAD: g_n > NCBF ! Increase NCBF')
               endif
               basis_n(g_n)=basis_n(g_n)+n_con_guass
               g_n=g_n+1
             enddo
             sh=sh+(14)*n_con_guass
             go to 6000

 6000       continue
     
           EndDo ! of ISNEW

CMI *************************************************************
CMI ***              we are out of the main cycle....
CMI *************************************************************

 201       continue

CMI ... important changes from SAJEEV !!!!
           n_pg_atom(k)=sh

           If(k.eq.1)then
             INPG=0
           else
             INPG = n_pg_atom(k-1)
           endif

           do ii=INPG+1,sh
              do isw=1,3
                R_ms(ii,isw)=re(k,isw)
              enddo
              IF (IPRINT.GE.5) THEN
               write(LUPRI,'(2X,I3,A,3F14.5)')
     &         ii,'. R_ms x/y/z:',R_ms(ii,1),R_ms(ii,2),R_ms(ii,3)
              ENDIF
           enddo
C ...................................................................................
       ENDDO ! of k

       close(LUinp,status="keep")

C ===================================================================
CMI   ... get the total number of contracted basis functions
C ===================================================================
       n_b=g_n-1
       
       WRITE(LUPRI,'(/2X,A,I3,A,I3)')
     & 'Total number of all primitive guassians / total number of'//
     & ' contracted basis functions=',sh,'/',n_b
       CALL FLSHFO(LUPRI)

       IF (n_b.NE.NCBF) THEN
        WRITE(LUPRI,*) 'n_b,NCBF=',n_b,NCBF
        CALL QUIT('INPREAD: n_b <> NCBF !')
       ENDIF

       IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A,I3)')
     & 'Number of R_ms elements (sh)=',sh
        DO I=1,sh
         WRITE(LUPRI,'(2X,A,I3,A,3F12.6)')
     &  'R_ms(',I,') x,y,z:',R_ms(I,1),R_ms(I,2),R_ms(I,3)
        ENDDO
        WRITE(LUPRI,*)
        CALL FLSHFO(LUPRI)
       ENDIF

       LUCAP=-1
       CALL GPOPEN(LUCAP,
     &  'CAPLBL.inp',' ',' ','FORMATTED',IDUMMY,.FALSE.)

       WRITE(LUCAP,*) 'CAP unique order of atomic basis functions'

       IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,*) 'CAP unique order of atomic basis functions'
       ENDIF
       do i=1,n_b
        IF (IPRINT.GE.3) THEN
          write(LUPRI,'(2i4,A6,A,I3)')i,SHELL_ATOM(i),Shell_basis(i),
     &   ' contr->',basis_n(i)  
        ENDIF
C        ... write into file LUCAP !
        WRITE(LUCAP,'(i4,2x,i3,2x,A6)') i,SHELL_ATOM(i),Shell_basis(i)
       enddo

       CALL GPCLOSE(LUCAP,'KEEP')

       IF (IPRINT.GE.1) THEN
         WRITE(LUPRI,'(2X,A)')
     &  'File "CAPLBL.inp" with CAP unique order of atomic basis'//
     &  ' functions ready.'
       ENDIF

CMI  ... set up important offset -  basis_ord(i)
      do i=1,n_b
         if(i.eq.1)then
           iasp=0
         else
           iasp=basis_ord(i-1)
         endif
         basis_ord(i)=iasp+basis_n(i)
         IF (IPRINT.GE.3) THEN
          WRITE(LUPRI,'(2X,A,I3,A,I3,A,I3)')
     &    'b.fu: ',i,' # of contr.shells:',basis_n(i),
     &    ' offset:',basis_ord(i)
          CALL FLSHFO(LUPRI)
         ENDIF
      enddo

CMI  ... control print out of exponents/contraction coefficients
C       WRITE(LUPRI,*) 'basis_n(i=1,n_b)',(basis_n(i), i=1,n_b )
C       WRITE(LUPRI,*) 'contraction coefficients:'

        ix  = 1
        do i=1,n_b

          PE_AO=0.0d0

          j = i

          if (i.eq.1) then
            iangi=0
          else
            iangi=basis_ord(i-1)
          endif

          iangj = iangi
 
          IF (IPRINT.GE.3) THEN
          write(lupri,'(/2X,A,I4)')
     &    '*** contracted function #',i

          write(lupri,'(2X,A,2I4,A,2I4)')
     &    'contraction coefficients (indx from, to)',
     &     ix, ix+basis_n(i)-1,' :',iangi+1,basis_ord(i)

CMI       write(lupri,*) (ord_coe(ix+k), k=0, basis_n(i)-1 )
          write(lupri,*) (ord_coe(ii),ii=iangi+1,basis_ord(i) )

          write(lupri,'(2X,A,2I4)')
     &    'corresponding exponents of primitives'//
     &    ' functions (indx from, to):',iangi+1,basis_ord(i)

          write(lupri,*) (ord_expon(ii),ii=iangi+1,basis_ord(i) )
          ENDIF

C         do ii=iangi+1,basis_ord(i)
C         enddo

CMI      .... normalization part for primitive gaussians first
C ======================================================================
          do ii=iangi+1,basis_ord(i)
           do isw=1,3
             R_m(isw)=R_ms(ii,isw)
             K_m(isw)=ang_pow(ii,isw)
           enddo
           R_m(4)=ord_expon(ii)
           CALL PEM(R_m,R_m,K_m,K_m,PEINT)
           ord_coe(ii)=ord_coe(ii)/sqrt(PEINT)
          enddo

          IF (IPRINT.GE.3) THEN
          write(lupri,*) '1.normalized contraction coefficients:'
          write(lupri,*) (ord_coe(ix+k), k=0, basis_n(i)-1 )
          ENDIF

CMI  ... second normalization - get PE_AO first
C ====================================================================
          do ii=iangi+1,basis_ord(i)
          do jj=iangj+1,basis_ord(j)
 
           do isw=1,3
             R_m(isw)=R_ms(ii,isw)
             R_n(isw)=R_ms(jj,isw)
             K_m(isw)=ang_pow(ii,isw)
             K_n(isw)=ang_pow(jj,isw)
           enddo

           R_m(4)=ord_expon(ii)
           R_n(4)=ord_expon(jj)

           COMUL=(ord_coe(ii)*ord_coe(jj))
           CALL PEM(R_m,R_n,K_m,K_n,PEINT)
           PE_AO= PE_AO+(COMUL*PEINT)

           IF (IPRINT.GE.10) THEN
            write(LUPRI,'(2I3,A,2F12.8,2(A,F12.8))')
     &     ii,jj,' coeff:',ord_coe(ii),ord_coe(jj),
     &     ' PEINT=',PEINT,' PE_AO=',PE_AO
           ENDIF

          enddo
          enddo

CMI  ... do second normalization - now apply PE_AO as the norm
          do ii=iangi+1,basis_ord(i)
           unnormcoe =  ord_coe(ii)
           ord_coe(ii) = unnormcoe/sqrt(PE_AO)
           IF (IPRINT.GE.3) THEN
            write (LUPRI,'(i4,A,A2,A,F12.5,A,F12.5,A,F12.5,A,F12.5)')
     &      ii,'.fu=',shell_arr(ii),
     &      ' expon:', ord_expon(ii),
     &      ' unnorm coeff=',unnormcoe,
     &      ' PE_AO=',PE_AO,' normalized coeff=',
     &      ord_coe(ii)
           ENDIF
          enddo

          ix = ix + basis_n(i)

       enddo

       CALL FLSHFO(LUPRI)

       CALL QEXIT('INPREAD')
       RETURN
       END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck aointegral */
      SUBROUTINE AOINTEGRAL(
     &  N_B,
     &  ord_expon,
     &  ord_coe,
     &  R_ms,
     &  basis_ord,
     &  ang_pow,
     &  PW,
     &  PE_AOr,PE_AOi,
     &  PEXR,PEXI,PEYR,PEYI,PEZR,PEZI,
     &  NCBF,NUBF) 
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
C           R_ms         nuclear centers
C           basis_ord
C           ang_pow
C           PW - theta
C           PE_AOr,PE_AOi, PEXR... arrays for integrals
C
C  Called from SAJEEV
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbcap.h"

      integer *4 n_b

C     integer *4 ang_pow(250,3)
      integer *4 ang_pow(NUBF,3)

      integer *4 k_m(3),k_n(3)

      integer *4 basis_ord(NCBF)
     
CMI ... entering parameter
C     real *8 ord_expon(250),ord_coe(250)
      real *8 ord_expon(NUBF),ord_coe(NUBF)
      real *8 R_m(4),R_n(4)

CMI ... entering parameter
C     real *8 R_ms(200,3)
      real *8 R_ms(NUBF,3)

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
              IF( DCADD)  PE_AOr(IPU)=0.d0
             ENDIF
            ENDIF
         end do
      end do

      IF (IPRCAP.GE.3) THEN

       CALL HEADER(
     & 'Output from AOINTEGRAL - basic info on basis functions',-1)

       do ii=1,n_b

         if (ii.eq.1) then
            iangii=0
         else
            iangii=basis_ord(ii-1)
         endif

         WRITE(LUPRI,'(2X,A,I3,A,I4,A,I3)')
     &   '***  contr.bas.fu: ',ii,' # offset:',iangii,
     &   ' # of prim.fu:',basis_ord(ii)-iangii

         Do i=iangii+1,basis_ord(ii)
          write(LUPRI,'(4X,A,I3,A,3I2)')
     &    'subset.fu:',i,' ang.pow. x,y,z:',
     &    ang_pow(i,1),ang_pow(i,2),ang_pow(i,3)

          write(LUPRI,'(4X,A,3F15.7)')
     &    'R_ms:',R_ms(i,1),R_ms(i,2),R_ms(i,3)

          write(LUPRI,'(4X,A,F15.7,A,F15.7)')
     &    'expon:',ord_expon(i),
     &    ' normalized coeff:',ord_coe(i)
         EndDo
       enddo
       CALL FLSHFO(LUPRI)
      ENDIF

C ******************************************************************
C     ... run over basis functions and calculate own integrals ...
C ******************************************************************
      IF (IPRCAP.GE.4) THEN
       CALL HEADER(
     & 'Output from AOINTEGRAL - generation of integrals',-1)
      ENDIF

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
         IF (IPRCAP.GE.4) THEN
          WRITE(LUPRI,'(2X,A,2I3,A,I3,$)')
     &    '... calculating CAP integral(s) of ',ii,jj,'/',n_b
         ENDIF
C        ... run over primitives in contracted basis functions
C        ... i,j points to position of exponents !!!
         Do i=iangii+1,basis_ord(ii)
         do j=iangjj+1,basis_ord(jj)
            do isw=1,3
               R_m(isw)=R_ms(i,isw)
               R_n(isw)=R_ms(j,isw)
               k_m(isw)=(ang_pow(i,isw))
               k_n(isw)=(ang_pow(j,isw))
            enddo
            R_m(4)=ord_expon(i)
            R_n(4)=ord_expon(j)

            IF (CAPOVLP) THEN
CMI         ... get overlap integrals
             CALL PEM(R_m,R_n,K_m,K_n,PEINTr)
             PEINT=DCMPLX(PEINTr,0.0d0)
            ELSE
CMI         ... get CAP integrals    (nonrel+Vdc, relativ, Vdc)
             CALL CAPM(R_m,R_n,K_m,K_n,PW,PEINT,
     &          PEINT_D1X,PEINT_D1Y,PEINT_D1Z)
            ENDIF

C           ... contraction coefficients of i,j
            COMUL=ord_coe(i)*ord_coe(j)
C           ... cumulate real and imaginary part of integrals
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

              IF (DCADD) PE_AOr(IPU)= PE_AOr(IPU)+COMUL*DREAL(PEINT)

             ENDIF
            ENDIF
         END DO 
         END DO

         IF (IPRCAP.GE.4) THEN
           WRITE(LUPRI,'(A)') ' ... done.'
           CALL FLSHFO(LUPRI)
         ENDIF
C        ... control print out of integrals
         IF (IPRCAP.GE.5) THEN
          IF (CAPOVLP) THEN
           WRITE(LUPRI,'(2X,A,I3,A,I3,A,F22.8)') 
     &    'CAP overlap(',ii,';',jj,')=',PE_AOr(IPU)
          ELSE IF (DONRCAP) THEN
           WRITE(LUPRI,'(2X,A,I3,A,I3,A,2F22.8)') 
     &    'CAP nonrel re,im(',ii,';',jj,')=',PE_AOr(IPU),PE_AOi(IPU)
          ELSE IF (.NOT.DONRCAP) THEN
C         ... relativistic stuff ...
           WRITE(LUPRI,'(2X,A,I3,A,I3,A,2F22.8)') 
     &    'CAPd_X re,im(',ii,';',jj,')=',PEXR(IPU),PEXI(IPU)
           WRITE(LUPRI,'(2X,A,I3,A,I3,A,2F22.8)') 
     &    'CAPd_Y re,im(',ii,';',jj,')=',PEYR(IPU),PEYI(IPU)
           WRITE(LUPRI,'(2X,A,I3,A,I3,A,2F22.8)') 
     &    'CAPd_Z re,im(',ii,';',jj,')=',PEZR(IPU),PEZI(IPU)

           IF (DCADD) THEN
             WRITE(LUPRI,'(2X,A,I3,A,I3,A,F22.8)') 
     &      '..scalar Vdc (',ii,';',jj,')=',PE_AOr(IPU)
           ENDIF
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

      SUBROUTINE PEM(A,B,KA,KB,PEINT)
      IMPLICIT REAL*8 (A-H,O-Z)

      real *8 A(4),B(4)
      integer *4 KA(3),KB(3)
      real *8 PEINT,PREMULT
      real *8 GAMMA
      real *8 AB(3),P(3),PA(3),PB(3)
      REAL *8 OVINT(3)

      GAMMA=(A(4))+(B(4))
      DO ig=1,3
         AB(ig)=(A(ig))-(B(ig))
         P(ig)=((A(4))*(A(ig)) +(B(4))*(B(IG)))/GAMMA
         PA(IG)=((P(IG))-(A(IG)))
         PB(IG)=((P(IG))-(B(IG)))
      END DO

      ABSQ=(AB(1))*(AB(1))+(AB(2))*(AB(2))+(AB(3))*(AB(3))
      PREMULT=DEXP(((-A(4))*(B(4))*ABSQ)/GAMMA)

      DO I=1,3
         KAI=KA(I)
         KBI=KB(I)
         PAI=PA(I)
         PBI=PB(I)
         CALL OVMAT(GAMMA,KAI,KBI,PAI,PBI,OVINTI)
         OVINT(I)=OVINTI
      END DO

      PEINT=PREMULT*(OVINT(1))*(OVINT(2))*(OVINT(3))

      RETURN
      END

      SUBROUTINE OVMAT(GAMMA,L1,L2,AP,BP,OVINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER *4 LI,L2
      REAL *8 GAMMA,KA,KB,AP
      REAL *8 PI,GAMMA1,GAMMA2,PREMULT
      REAL *8 BINC
      INTEGER *4 IFI
CMI   EXTERNAL ISUMM
      INTEGER*4 ISUMM
      data PI/3.141592653589793/

      OVINT=0.0d0
      DO i=0,(L1+L2)/2
         IFI=ISUMM(2*i)
         GAMMA1=(2.0D0*GAMMA)**(i)
         GAMMA2=(PI/GAMMA)**(0.50d0)
         PREMULT=IFI*GAMMA2/GAMMA1
         CALL BINCOEF(2*i,L1,L2,AP,BP,BINC)
         OVINT=(OVINT)+(PREMULT*BINC)
      END DO
      RETURN
      END

      SUBROUTINE  BINCOEF(K,L1,L2,AP,BP,BINC)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER *4 K,L1,L2
      REAL *8 AP,BP,BINC
      REAL *8 COEL1,COEL2,PRM1,PRM2,PRM3
      REAL *8 PREMULT
      EXTERNAL BINEXP
      BINC=0.0D0

      DO I=0,L1
         DO J=0,L2
            IF((I+J).EQ.K)THEN
               COEL1=BINEXP(L1,I)
               COEL2=BINEXP(L2,J)
               PRM1=COEL1*COEL2
               IL=L1-i
               JL=L2-J
               PRM2=(AP)**(IL)
               PRM3=(BP)**(JL)
               PREMULT=PRM1*PRM2*PRM3
               BINC=BINC+PREMULT
            ENDIF
         END DO
      END DO
      RETURN
      END

      FUNCTION BINEXP(L,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 BINEXP
      INTEGER *4 L,I,LP1,LP2,L3
      INTEGER *4 NFACT
      EXTERNAL NFACT
      LP1=NFACT(L)
      LP2=NFACT(I)
      LP3=NFACT(L-I)
      BINEXP=(LP1)/(LP2*LP3)
      END

      FUNCTION NFACT(I)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER *4 I,K,KL,NFACT
      KL=1
      DO K=1,I
         KL=KL*K
      END DO
      NFACT=KL
      END

CMI ... this is INTEGER*4 function !!      
CMI   FUNCTION ISUMM(I)
      INTEGER*4 FUNCTION ISUMM(I)
      IMPLICIT REAL*8 (A-H,O-Z)
CMI   INTEGER *4 I,K,KL,M,ISUMM
      INTEGER *4 I,K,KL,M
      KL=1
      Do k=1,i
         m=(k-((k/2)*2))
         if(m.eq.0)go to 110
         KL=kL*k
 110     continue
      End Do
      ISUMM=KL
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PEM */
      SUBROUTINE PEMOLD(R_m,R_n,K_m,K_n,PEINT)
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

         CALL SMAT(R_KM,R_KN,K_KM,K_KN,OVRLPI)

         CALL SMAT(R_LM,R_LN,K_LM,K_LN,OVRLPJ)

         CALL SMAT(R_CM,R_CN,K_CM,K_CN,OVRLPK)

         PEINT=PEINT+(OVRLPI*OVRLPJ*OVRLPK)

      CALL QEXIT('PEM')
      RETURN
      END

      subroutine SMAT(A,B,L1,L2,OVINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer *4 L1, L2
      real *8  A(2),B(2)
      real *8 GAMMA,P,PA,PB, AB, ABSQ
      real *8 OVINT,RMULT,PREMULT, PI
      INTEGER *4 IFI
      real *8 GAMMA1, GAMMA2
      real *8 BINC
      data PI/3.141592653589793/
CMI   EXTERNAL  ISUMM
      INTEGER*4 ISUMM

      GAMMA=(A(2))+(B(2))
      P=((A(1))*(A(2))+(B(1))*(B(2)))/GAMMA
      PA=P-A(1)
      PB=P-B(1)
      AB=(A(1))-(B(1))
      OVINT=0.0d0

      DO i=0,(L1+l2)/2
         IFI=ISUMM(2*i)
         GAMMA1=(2.0D0*GAMMA)**(i)
         GAMMA2=(PI/GAMMA)**(0.50d0)
         RMULT=IFI*GAMMA2/GAMMA1
         CALL BINCOEF(2*i,L1,L2,PA,PB,BINC)
         OVINT=OVINT+(RMULT*BINC)
      END DO
      ABSQ=AB*AB
      PREMULT=exp((-(A(2))*(B(2))*ABSQ)/gamma)
      OVINT=OVINT*PREMULT

      RETURN
      END


 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck SMAT */
      subroutine SMAT_OLD(R_m,R_n,K_m,K_n,BETA,Int_GR)
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

CMI ... must be outside according to Sajeev !...
      INT_GR=0.d0

      DO i=0,(K_m+k_n)/2
C         if (i.eq.1)write(*,*)i
CMI      INT_GR=0.d0
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
     &           VCAP,VCAP_D1X,VCAP_D1Y,VCAP_D1Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Called from AOINTEGRAL (from within the 4-fold cycle)
C
C   R_m,R_n,K_m,K_n
C 
C   entering parameter "theta" was from the INPQD file
C
C   returns: 
C --------------
C          VCAP - the total nonrel. scalar CAP term (DONRCAP=.true.)
C
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
      
CMI    ... fix i,j,k so that is always 1,2,3 (like x,y,z)
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
     &                   VCAP1C,VCAP1,X0DC(i),EPSDC(i))

C           ... accumulate VCAPS(i) for each interval
            IF (DONRCAP) THEN
              VCAPS(i) = VCAPS(i) + VCAP1
            ELSE
C           ... relativistic integrals ...
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
              IF (DCADD) THEN
C               ... cumulate Vdc x,y,z scalar contributions
                VCAPS(i) = VCAPS(i) + VCAP1
              ENDIF
            ENDIF

         enddo

         CALL SMAT(R_KM,R_KN,K_KM,K_KN,OVRLPK)
         CALL SMAT(R_LM,R_LN,K_LM,K_LN,OVRLPL)

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
           IF (DCADD) THEN
              VCAPS(i)=VCAPS(i)*OVRLPK*OVRLPL
           ENDIF
         ENDIF

      ENDDO ! of I=1,3
             
C   ... finally add all theree terms - for x,y and z - for DONRCAP=.true., else
C     for relativistic case where DCADD=.true.
C --------------------------------------------------------------------------------
      VCAP=VCAPS(1)+VCAPS(2)+VCAPS(3)
C --------------------------------------------------------------------------------
      IF (IPRCAP.GE.10) THEN
       IF (.NOT.DONRCAP) THEN
        WRITE(LUPRI,'(/2X,A,3(2F15.7))')
     &  '** CAPM >> relativistic VCAP_D1X,VCAP_D1Y,VCAP_D1Z:',
     &  VCAP_D1X,VCAP_D1Y,VCAP_D1Z 
        IF (DCADD) THEN
        WRITE(LUPRI,'(2X,A,2F15.7)')
     &  '...scalar VDC=',VCAP
        ENDIF
       ELSE
        WRITE(LUPRI,'(2X,A,2F15.7)')
     &  '** CAPM >> nonrelativ. VCAP=',VCAP
       ENDIF
      ENDIF

      CALL QEXIT('CAPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck capmat1 */
      SUBROUTINE CAPMAT1(CAP_PAR,GRSTP,NSTEP,R_m,R_n,K_m,K_n,
     &                   VCAP1C,VCAPT, X_DC, EPS_DC)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C On input:
C
C   CAP_PAR, GRST, NSTEP
C   R_m, R_n, K_m, K_n
C
C   X_DC, EPS_DC ... if DCADD=.true., these are DC field components
C
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
      real *8  X_DC, EPS_DC, VF
      real *8 X0,theta,LAMDA,GRSTP(2)
      complex *16 FFx,fx,Dfx,D2fx
      complex *16 V0x,V1x,V2x
      complex *16 V0xy,V1xy,V2xy

      complex *16  TempV1x,TempV1xy

      complex *16 CI,ROT
      complex *16 PEM,XS,PEN

      COMPLEX *16 VCAPT,VCAPXT,VCAPYT,VCAPZT
      COMPLEX *16 VCAPX(4),VCAPX_D1C(4),VCAP1C
      COMPLEX *16 VCAP,VADD
      COMPLEX *16 VCAP_D1C,VCAPy_D1C
      COMPLEX *16 VCAPTy,VCAPXTy,VCAPYTy,VCAPZTy
      COMPLEX *16 VCAPXy(4),VCAPXy_D1C(4)
      COMPLEX *16 VCAPy,VADDy
      COMPLEX *16 VCAPXT_D1C,VCAPXTy_D1C
      EXTERNAL V0,V1,V2,FF
      EXTERNAL V0D,V1D
CMI   ... declaration must be for gfortran
      double complex v0,v1,v2,V1D,FF,V0D
      
      data CI/(0.0d0,1.0d0)/
      data M/1.0d00/
      data TOL/0.00000001/

      CALL QENTER('CAPMAT1')

CMI   ... calculate the integration step (per interval)
      STEP=(GRSTP(2)-GRSTP(1))/(NSTEP)

CMI   ... extract X0, lambda, theta parameters
      X0=CAP_PAR(1)
      LAMDA=CAP_PAR(2)
      THETA=CAP_PAR(3)
    
CMI   ... initialization of summation variables
      VCAPXT =(0.d0,0.d0)
      VCAPXTy=(0.d0,0.d0)

      VCAPXT_D1C =(0.d0,0.d0)
      VCAPXTy_D1C=(0.d0,0.d0)

CMI   ... starting position for the integration
      XI=GRSTP(1)

      IF (IPRCAP.GE.46) THEN
CMI    ... info for integration parameters  ...
       WRITE(LUPRI,'(2X,A,F12.5,A,I7)')
     & 'Trapezoidal 1-D integration starting from point',XI,
     & ' number of subintervals :',NSTEP
      ENDIF

CMI  ... run over the integration interval with NSTEP nodes
      itx=0
      DO ix=1,NSTEP

C        ... actual position for the integration in +/- ways...
         X=(XI+(ix-1)*STEP)
         Xy=X*(-1.0d0)

         itx=itx+1

         XR=X-R_n(1)
         XRM=X-R_m(1)

         IF(dabs(XR).lt.TOL) XR=TOL ! important for denominator
         DGr=K_n*(XR**(-1))

C         ... 1st derivative
         DG=DGr-(2*R_n(2)*(XR))
       
         XRy=Xy-R_n(1)
         IF(dabs(XRy).lt.TOL)XRy=TOL

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

CMI      ... additional electric field
         IF (DCADD) THEN
           IF (X.GE.X_DC) THEN
!          ... VF must be always positive !!!!
            VF = X - X_DC 
            VADD=DCMPLX(-VF*EPS_DC,0.0d0)
!           VADD=-VF*EPS_DC
           ELSE
            VADD=DCMPLX(0.0d0,0.0d0)
           ENDIF
         ELSE
            VADD=DCMPLX(0.0d0,0.0d0)
         ENDIF

         IF (DONRCAP) THEN
CMI        ... nonrelativistic CAP operator
           v0x=V0(M,THETA,LAMDA,X0,X)
           V1x=v1(M,THETA,LAMDA,X0,X)*DG
           v2x=v2(M,THETA,LAMDA,X0,X)*D2G
           v0xy=V0(M,THETA,LAMDA,X0,Xy)
           V1xy=v1(M,THETA,LAMDA,X0,Xy)*DGy
           v2xy=v2(M,THETA,LAMDA,X0,Xy)*D2Gy
           VCAP =(V0x +V1x +V2x +VADD)*GG   ! ... include Vdc contribution through VADD
           VCAPy=(V0xy+V1xy+V2xy-VADD)*GGy  !  
           VCAPX(itx)=VCAP
           VCAPXy(itx)=VCAPy
           if(ix.eq.1) go to 321 
CMI        ... calculate integral - trapezoidal method
           VCAPXT=VCAPXT+0.5*(VCAPX(2)+VCAPX(1))*STEP
           VCAPXTy=VCAPXTy+0.5*(VCAPXy(2)+VCAPXy(1))*STEP
           if (ix.ne.1) then
              VCAPX(2)=VCAPX(1)
             VCAPXy(2)=VCAPXy(1)
           endif
         ELSE
CMI       ... relativistic Hamiltonian with x,y,z components
CMI       ... Vdc comes in individual terms ...
          IF (.NOT.CAPDIPVEL) THEN
           V0x  = V0D(M,THETA,LAMDA,X0,X )
           V0xy = V0D(M,THETA,LAMDA,X0,Xy)

           V1x  = V1D(M,THETA,LAMDA,X0,X )*DG
           V1xy = V1D(M,THETA,LAMDA,X0,Xy)*DGy

           VCAP_D1C  = (V0x +V1x       )*GG
           VCAPy_D1C = (V0xy+V1xy      )*GGy
          ELSE
CMI       ... CAPDIPVEL = .TRUE.; comparison with DIRAC d/dx,d/dy,d/dz analytical integrals
           VCAP_D1C  = GG * DG
           VCAPy_D1C = GGy* DGy
C          write(LUPRI,*) 'VCAP_D1C, VCAPy_D1C:',VCAP_D1C, VCAPy_D1C
          ENDIF
          VCAPX_D1C(itx) = VCAP_D1C
          VCAPXy_D1C(itx)= VCAPy_D1C

          if(ix.eq.1) go to 321 

C           ... numerical (trapezoidal) integration 
           VCAPXT_D1C=VCAPXT_D1C  +
     &         0.50D0*(VCAPX_D1C(2)+VCAPX_D1C(1))*STEP
           VCAPXTy_D1C=VCAPXTy_D1C+
     &         0.50D0*(VCAPXy_D1C(2)+VCAPXy_D1C(1))*STEP

           IF (DCADD.AND..NOT.CAPDIPVEL) THEN
C          ... scalar contribution from DC field
             VCAP =(VADD)*GG
             VCAPy=(-VADD)*GGy
             VCAPX(itx)=VCAP
             VCAPXy(itx)=VCAPy
             if(ix.eq.1) go to 321 
             VCAPXT=VCAPXT+0.5*(VCAPX(2)+VCAPX(1))*STEP
             VCAPXTy=VCAPXTy+0.5*(VCAPXy(2)+VCAPXy(1))*STEP
             if (ix.ne.1) then
               VCAPX(2) =VCAPX(1)
               VCAPXy(2)=VCAPXy(1)
             endif
           ENDIF

C          ...rearrange for numerical integration 
           if (ix.ne.1) then
             VCAPX_D1C(2) =VCAPX_D1C(1)
             VCAPXy_D1C(2)=VCAPXy_D1C(1)
           endif 

         ENDIF

CMI     ... restore the counter
         itx=0

 321     CONTINUE

      ENDDO

      IF (DONRCAP) THEN
CMI     ... nonrelativistic CAP
        VCAPT=VCAPXT+VCAPXTy
      ELSE
CMI    ... summ of 2 integrations (0,+inf) + (0,-inf)
        VCAP1C=VCAPXT_D1C+VCAPXTy_D1C
        IF (DCADD) THEN
          VCAPT=VCAPXT+VCAPXTy ! take into account scalar Vdc term
        ENDIF
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

CMI   ...  newly derived V0 for the DIRAC Hamiltonian
! ============================================================
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

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck check_pemat */
      LOGICAL FUNCTION CHECK_PEMAT()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C      Check the presence of the PEMAT file and  of integral labels
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbcap.h"
      LOGICAL ISFILE,IS_OPEN,FNDLAB
      CHARACTER*8 LBL

      CHECK_PEMAT = .FALSE.
      INQUIRE(FILE='PEMAT',EXIST=ISFILE,OPENED=IS_OPEN)
      IF (ISFILE) THEN
        IF (.NOT.IS_OPEN) THEN
          CALL OPNFIL(LUPEMAT,'PEMAT','OLD','SAJEEV')
        ENDIF
        REWIND(LUPEMAT)
        IF (CAPOVLP) THEN
          LBL = 'CAP_OVLP'
          IF (FNDLAB(LBL,LUPEMAT)) CHECK_PEMAT = .TRUE.
        ELSE
         IF  (DONRCAP) THEN
          LBL = 'CAP_RE  '
         ELSE
          LBL = 'ZCAPD1I' ! relativistic case 
         ENDIF
        ENDIF
        IF (FNDLAB(LBL,LUPEMAT)) CHECK_PEMAT = .TRUE.
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pemat_write */
      SUBROUTINE PEMAT_WRITE(PE_AOr,PE_AOi,
     &            PEXR,PEXI,PEYR,PEYI,PEZR,PEZI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Save CAP like integrals (in entering arrays) into PEMAT file
C
C  Written by Miro Ilias, Sept. 2007, Tel Aviv
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbcap.h"
      CHARACTER*8 LBL
      LOGICAL FNDLAB
      REAL*8 PE_AOr(*),PE_AOi(*),PEXR(*),PEXI(*),
     &       PEYR(*),PEYI(*),PEZR(*),PEZI(*)

      CALL QENTER('PEMAT_WRITE')

C        do ii=1,n_b
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

              IF (DCADD) THEN
CMI            ... this is separate scalar contribution from Vdc term
               LBL='CAPD_VDC'
               CALL NEWLAB(LBL,LUPEMAT,LUPRI)
               CALL WRITT(LUPEMAT,NNBBASX,PE_AOr)
               WRITE(LUPRI,'(2X,A,A)')
     &        'CAP relativistic Vdc scalar contribution generated,'//
     &        ' PEMAT label:',LBL
              ENDIF

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

      CALL QEXIT('PEMAT_WRITE')
      RETURN
      END

************************************************************************
      DOUBLE PRECISION FUNCTION F(N)
C     ------------------------------------------------------------------
C     Returns the factorial N!. 
C     ------------------------------------------------------------------
      IMPLICIT NONE
C     ... input parameters ...
      INTEGER N
C     ... internal variables ...
      INTEGER I
C     ------------------------------------------------------------------
CMI   IF (N.LT.0) PAUSE 'ERROR: FACTORIAL UNDEFINED'
      IF (N.LT.0) CALL QUIT('ERROR: FACTORIAL UNDEFINED')
      F = 1.0D0
      DO I = 2,N
        F = F*DBLE(I)
      ENDDO
C     ------------------------------------------------------------------
C     Last modified: 05. 10. 2003
      END
************************************************************************
************************************************************************
      DOUBLE PRECISION FUNCTION NK(N,K)
C     ------------------------------------------------------------------
C     Returns the binomical coefficient {N \choose K}.
C     ------------------------------------------------------------------
      IMPLICIT NONE
C     ... input parameters ...
      INTEGER N,K
C     ... externals ...
      EXTERNAL F
      DOUBLE PRECISION F
C     ------------------------------------------------------------------
CMI   IF (N.LT.0) PAUSE 'ERROR (NK): N < 0'
      IF (N.LT.0) CALL QUIT('ERROR (NK): N < 0')
CMI   IF (K.LT.0) PAUSE 'ERROR (NK): K < 0'
      IF (K.LT.0) CALL QUIT('ERROR (NK): K < 0')
CMI   IF (K.GT.N) PAUSE 'ERROR (NK): K > N'
      IF (K.GT.N) CALL QUIT('ERROR (NK): K > N')
      NK = F(N)/(F(N-K)*F(K))
C     ------------------------------------------------------------------
C     Last modified: 05. 10. 2003
      END
************************************************************************
************************************************************************
      DOUBLE PRECISION FUNCTION TRIG(M,N,U,V)
C     ------------------------------------------------------------------
C     Returns the value of a trigonometric integral 
C     TRIG = int_{U}^{V} [sin x]^M [cos x]^N dx . 
C     ------------------------------------------------------------------
      IMPLICIT NONE
C     ... input parameters ...
      INTEGER M,N
      DOUBLE PRECISION U,V
C     ... constants ...
      DOUBLE PRECISION TOL
      PARAMETER(TOL=1.0D-15)
      DOUBLE COMPLEX C0,IM,C2
      PARAMETER(C0=(0.0D0,0.0D0))
      PARAMETER(IM=(0.0D0,1.0D0))
      PARAMETER(C2=(2.0D0,0.0D0))
C     ... internal variables ...
      INTEGER JM,JN,K
      DOUBLE COMPLEX CTRIG,CAUX1,CAUX2,CAUX3
C     ... externals ...
      EXTERNAL NK
      DOUBLE PRECISION NK
C     ------------------------------------------------------------------
CMI   IF (M.LT.0) PAUSE 'ERROR (TRIG): M < 0'
      IF (M.LT.0) CALL QUIT('ERROR (TRIG): M < 0')
CMI   IF (N.LT.0) PAUSE 'ERROR (TRIG): N < 0'
      IF (N.LT.0) CALL QUIT('ERROR (TRIG): N < 0')
C     ------------------------------------------------------------------
      CTRIG = C0
      DO JM = 0,M
      DO JN = 0,N
C       ----------------------------------------------------------------
        CAUX1 = DCMPLX(NK(M,JM)*NK(N,JN)*((-1.0D0)**JM))
C       ----------------------------------------------------------------
        K = M + N - 2*(JM+JN)
C       ----------------------------------------------------------------
        IF (K.EQ.0) THEN
          CAUX1 = CAUX1*DCMPLX(V-U)
        ELSE
          CAUX2 = IM*DCMPLX(DBLE(K))
          CAUX3 = CDEXP(CAUX2*DCMPLX(V)) - CDEXP(CAUX2*DCMPLX(U))
          CAUX3 = CAUX3/CAUX2
          CAUX1 = CAUX1*CAUX3
        ENDIF
C       ----------------------------------------------------------------
        CTRIG = CTRIG + CAUX1
C       ----------------------------------------------------------------
      ENDDO
      ENDDO
      CTRIG = CTRIG/((C2**M)*(IM**M)*(C2**N))
C     ------------------------------------------------------------------
      IF (DABS(DIMAG(CTRIG)).GT.TOL) 
CMI  &   PAUSE 'ERROR (TRIG): OUTPUT IS NOT A REAL NUMBER'
     &   CALL QUIT('ERROR (TRIG): OUTPUT IS NOT A REAL NUMBER')
C     ------------------------------------------------------------------
      TRIG = DREAL(CTRIG)
C     ------------------------------------------------------------------
C     Last modified: 05. 10. 2003
      END
************************************************************************
