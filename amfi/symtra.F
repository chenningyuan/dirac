C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
      Subroutine SYMTRAFO(LUPROP,LUPROP2,PROPRI,IPRINT,WRK,LFREE)
cbs   
cbs   Purpose: combine SO-integrals from amfi to symmetry-addapted   
cbs   integrals on one file AOPROPER_MF_SYM   
cbs   
cbs   LUPROP has the atomic integrals, LUPROP2 will have the molecular 
cbs   integrals    
cbs   
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <para.h>
      parameter(maxorbs=500)
      parameter(maxcent=40)
      character*8 xa,ya,za,xa2
      character*3 END                     
CBS   character*20 filename               
      logical  EX, PROPRI
CBS   namelist /SYMTRA/ none
      dimension xa(4),ya(4),za(4)
      dimension xa2(4)
      dimension WRK(LFREE)
      dimension ncent(maxorbs), Lval(maxorbs),mval(maxorbs), 
     *nadpt(maxorbs),nphase(8,maxorbs),jdummy(8),    
     *Lhighcent(maxcent),Lcent(maxorbs,maxcent),Mcent(maxorbs,maxcent), 
     *ncontcent(0:Lmax,maxcent),ipC(3,Maxcent),
     *numballcart(maxcent),ifirstLM(0:Lmax,-Lmax:Lmax,maxcent)  ,
     *ncount(maxorbs)
      Character Label*8
c###############################################################################
      IPNT(I,J)=(max(i,j)*max(i,j)-max(i,j))/2 +min(i,j) 
      zero=0d0
      xa(1)='********'
      ya(1)='********'
      za(1)='********'
      xa(2)='        '
      ya(2)='        '
      Za(2)='        '
      xa(3)='ANTISYMM'
      ya(3)='ANTISYMM'
      Za(3)='ANTISYMM'
      xa(4)='X1MNF-SO'
      ya(4)='Y1MNF-SO'
      ZA(4)='Z1MNF-SO'  
c
c     read information from SYMINFO   
      ISYMUNIT = 0
      inquire(FILE='SYMINFO',EXIST=EX) 

      if (.not.EX) CALL QUIT('SYMTRAFO: SYMINFO not present. Sorry')
CMI   if (.not.EX) RETURN

      CALL GPOPEN(ISYMUNIT,'SYMINFO','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      rewind(isymunit)  

      write(LUPRI,*) 'SYMTRAFO: Symmetry adapation of the SO-integrals'

      read(isymunit,*)   
      read(isymunit,*)   
      read(isymunit,*)   
      numboffunct=0
      do while(END.ne.'END') 
      numboffunct=numboffunct+1   
      read(isymunit,'(A3)') END   
      enddo    

      write(LUPRI,*) 'there are totally ',numboffunct,' functions' 

      if (numboffunct.gt.maxorbs) 
     *CALL QUIT('increase maxorbs in symtrafo')
      rewind isymunit 
      read(isymunit,*)   
      read(isymunit,*)
      numbofcent=0
CBS   read(5,SYMTRA,END=4711)
4711  do irun=1,numboffunct
      read(isymunit,*) index,ncent(irun),lval(irun),
     *mval(irun),ncount(irun),nadpt(irun)
     *,(nphase(I,irun),I=1,nadpt(irun)) 
      numbofcent=max(numbofcent,ncent(irun)) 
      if (index.ne.irun) CALL QUIT('weird numbering  on SYMINFO')
      enddo 
      CALL GPCLOSE(ISYMUNIT,'KEEP')

      write(LUPRI,*) 'number of unique centres' , numbofcent 
c
c     clean up arrays for new integrals  
      numboffunct3=(numboffunct*numboffunct+numboffunct)/2 
*
      nSCR=numboffunct3*3
      IPX = 1
      ipY=ipX+numboffunct3
      ipZ=ipY+numboffunct3
      IPSCR = IPZ + NUMBOFFUNCT3
      KLAST = IPZ + NSCR
      IF (KLAST .GT. LFREE) CALL STOPIT('SYMTRA',' ',KLAST,LFREE)
      Call DCopy(numboffunct3*3,Zero,0,WRK(ipx),1)
*
*
c     
c     loop over unique centres to read integrals and information   
c
      iunit=LUPROP
      nSCR=numboffunct3*3
      ipSCRX=ipSCR
      ipSCRY=ipSCRX+numboffunct3
      ipSCRZ=ipSCRY+numboffunct3    
      KLAST = IPZ + NSCR
      IF (KLAST .GT. LFREE) CALL STOPIT('SYMTRA',' ',KLAST,LFREE)
      length3_tot=0
*
      do icent=1,numbofcent    
*
      write(LUPRI,*) 'read integrals and info for centre ',icent   

              read(iunit)  xa2
     *        ,numbofsym,(jdummy(I),
     *        i=1,numbofsym),
     *        numballcart(icent),(Lcent(i,icent),
     *        I=1,numballcart(icent)),
     *        (mcent(i,icent),I=1,numballcart(icent)),
     *        Lhighcent(icent),(ncontcent(I,icent),I=0,Lhighcent(icent))

              write(LUPRI,*) numballcart(icent) , 
     &        ' functions on centre ',icent

              length3=ipnt(numballcart(icent),numballcart(icent)) 
              ipC(1,iCent)=ipSCRX
              ipC(2,iCent)=ipSCRY
              ipC(3,iCent)=ipSCRZ
              read(iunit) (Wrk(i),i=ipSCRX,ipSCRX+length3-1)
              read(iunit)  Ya
              read(iunit) (Wrk(i),i=ipSCRY,ipSCRY+length3-1)
              read(iunit)  Za
              read(iunit) (Wrk(i),i=ipSCRZ,ipSCRZ+length3-1)
              ipSCRX=ipSCRX+length3
              ipSCRY=ipSCRY+length3
              ipSCRZ=ipSCRZ+length3
              length3_tot=length3_tot+length3
culf
c      check if any L-value is missing
      LLhigh = Lhighcent(icent)
      do i=1,Lhighcent(icent)
         if(ncontcent(I,icent).eq.0) LLhigh=LLhigh-1
      enddo
      Lhighcent(icent)=LLhigh                                  
cbs   determize where the first function of a special type is..
      do Lrun=0,Lhighcent(icent) 
      do Mrun=-Lrun,Lrun 
      ifirstLM(Lrun,Mrun,icent)=ipnt(maxorbs,maxorbs)+1 
      enddo 
      enddo 
      do iorb=1,numballcart(icent)  
      Lrun=Lcent(iorb,icent) 
      Mrun=Mcent(iorb,icent) 
c     write(LUPRI,*) 'iorb,Lrun,mrun',iorb,Lrun,mrun 
      ifirstLM(Lrun,Mrun,icent)=min(iorb,ifirstLM(Lrun,Mrun,icent))   
      enddo 
cbs   determined..     
cbs   check if all of them were found 
      do Lrun=0,Lhighcent(icent) 
      do Mrun=-Lrun,Lrun 
      if(ifirstLM(Lrun,Mrun,icent).eq.(ipnt(maxorbs,maxorbs)+1)) then 
      write(LUPRI,*)
     &   'SYMTRAFO: problems for centre,L,M ',icent,Lrun,Mrun  
      CALL QUIT('SYMTRAFO: problems with L- and M-values')
      endif 
      enddo 
      enddo 
      enddo    !end of loop over centres   
cbs
cbs   Finally the transformation!!!! 
cbs
cbs
      lauf=0
      do irun=1,numboffunct 
      icent=ncent(irun)
      ilval=lval(irun)
      imval=mval(irun)
      do jrun=1,irun        
      lauf=lauf+1
      jcent=ncent(jrun)
      jlval=lval(jrun)
      jmval=mval(jrun)
cbs   check for same centers 
      if (ncent(irun).eq.ncent(jrun)) then 
      if (lval(irun).eq.lval(jrun).and.lval(irun).gt.0) then 
      if (iabs(iabs(mval(irun))-iabs(mval(jrun))).le.1) then 
      if (iabs(iabs(mval(irun))+iabs(mval(jrun))).ne.0) then 
      if (irun.ne.jrun) then 
*
*
*
cbs   the only cases  where non-zero integrals occur    
      if (nadpt(irun).eq.1) then 
      coeff=1d0
      else 
      icoeff=0
      do icc=1,nadpt(irun)
      icoeff=icoeff+nphase(icc,irun)*nphase(icc,jrun)
      enddo 
      coeff=icoeff
CBS   coeff=coeff/(nadpt(irun)*nadpt(irun))  ! this is HERMIT  
      endif  
cbs   determine indices of atomic integrals 
      indexi=ifirstLM(ilval,imval,icent)+ncount(irun)-1
      indexj=ifirstLM(jlval,jmval,jcent)+ncount(jrun)-1
      laufalt=ipnt(indexi,indexj) 
*
      ipSCRX=ipC(1,ncent(irun))-1+laufalt
      ipSCRY=ipC(2,ncent(irun))-1+laufalt
      ipSCRZ=ipC(3,ncent(irun))-1+laufalt
      if (indexi.gt.indexj) then 
         Wrk(ipX+lauf-1)=coeff*Wrk(ipSCRX)
         Wrk(ipY+lauf-1)=coeff*Wrk(ipSCRY)
         Wrk(ipZ+lauf-1)=coeff*Wrk(ipSCRZ)
      else 
         Wrk(ipX+lauf-1)=-coeff*Wrk(ipSCRX)
         Wrk(ipY+lauf-1)=-coeff*Wrk(ipSCRY)
         Wrk(ipZ+lauf-1)=-coeff*Wrk(ipSCRZ)
      endif 
CBS
CBS  the integrals are now on on the WRK array 
CBS  WRK(ipx+i) (i=1,numboffunct3)
CBS  WRK(ipy+i) (i=1,numboffunct3)
CBS  WRK(ipz+i) (i=1,numboffunct3)
CBS
*
*   
      endif 
      endif 
      endif 
      endif 
      endif 
      enddo   
      enddo   
      if (lauf.ne.numboffunct3) CALL QUIT('error in numbering ')
*
      IPRINT=5
      IF (PROPRI .OR. IPRINT .GT. 4) THEN
         CALL AROUND('Integrals of operator: '//XA(4))
         CALL OUTPAK(WRK(IPX),NUMBOFFUNCT,1,LUPRI)
         CALL AROUND('Integrals of operator: '//YA(4))
         CALL OUTPAK(WRK(IPY),NUMBOFFUNCT,1,LUPRI)
         CALL AROUND('Integrals of operator: '//ZA(4))
         CALL OUTPAK(WRK(IPZ),NUMBOFFUNCT,1,LUPRI)
      END IF
      IPRINT=0
      CALL WRTPRO(WRK(IPX),NUMBOFFUNCT3,XA(4),XA(2))
      CALL WRTPRO(WRK(IPY),NUMBOFFUNCT3,YA(4),YA(2))
      CALL WRTPRO(WRK(IPZ),NUMBOFFUNCT3,ZA(4),ZA(2))
      Return
      End 

  
      Subroutine Jacobi(A,C,N,LdC)
************************************************************************
*                                                                      *
*     (c) Copyright. All rights reserved                               *
*                                                                      *
*     No part of this code may be copied, redestributed or included    *
*     into any commercial product without the written permission of    *
*     the author. The use is restriced to research purposes only.      *
*                                                                      *
************************************************************************
*                                                                      *
*     Jacobi diagonalization                                           *
*                                                                      *
*     calling parameters:                                              *
*     A       : array of double precision real, input/output           *
*               On input this is the matrix to be diagonalized.        *
*               On ouput the diagonal matrix is returned.              *
*     C       : array of double precision real                         *
*               On input this are the basis vectors.                   *
*               The rotated basis vectors are returned.                *
*     N       : Integer, input                                         *
*               Dimension of the matrix to be diagonalized             *
*     LdC     : Integer, input                                         *
*               Leading dimension of matrix C.                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher, University of Lund, Sweden, 1993                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
*
      Dimension A(*)
      Dimension C(LdC,*)
*
      Parameter ( Eps=1.D-12 )
      Parameter ( EpsSqr=Eps*Eps )
*----------------------------------------------------------------------*
*     initialisation step:                                             *
*     - compute the largest subdiagonal element                        *
*----------------------------------------------------------------------*
      ASqrMax=0.0  
      Do i=2,N
         ii=i*(i-1)/2
         Do j=1,i-1
            Temp=A(ii+j)*A(ii+j)
            ASqrMax=Max(ASqrMax,Temp)
         End Do
      End Do
      If ( ASqrMax.le.EpsSqr ) Return  
      Thrs=0.1*ASqrMax
*----------------------------------------------------------------------*
*     start loop over sweeps:                                          *
*     - no pivot elementsearch is implemented. However, the rotation   *
*       is applied only to those submatrix elements which are bigger   *
*       than a threshold. The latter is based on the value of the      *
*       Largest subdiagonal element.                                   *
*----------------------------------------------------------------------*
      iter=0
100   iter=iter+1
      Do i=2,N
         ii=i*(i-1)/2
         Do j=1,i-1
            jj=j*(j-1)/2
            Aij=A(ii+j)
            AijSqr=Aij*Aij
            If ( AijSqr.ge.Thrs ) then
*---  compute the rotation angle --------------------------------------*
               Aii=A(ii+i)
               Ajj=A(jj+j)
               Diff=Aii-Ajj
               SigRot=1.0
               If ( Diff.lt.0.0 ) then
                  SigRot=-SigRot
                  Diff=-Diff
               End If
               Temp=Diff+SQRT(Diff*Diff+4.0*AijSqr)
               TanA=2.0*SigRot*Aij/Temp
               CosA=1.0/SQRT(1.0+TanA*TanA)
               SinA=CosA*TanA
*---  update rows/columnes of the matrix A ----------------------------*
               kmin=1  
               kmax=j-1
               kleft=Mod(kmax-kmin+1,4)
               If ( kleft.eq.1 ) then
                  Aii1=A(ii+1)
                  Ajj1=A(jj+1)
                  A(ii+1)=SinA*Ajj1+CosA*Aii1
                  A(jj+1)=CosA*Ajj1-SinA*Aii1
               Else If ( kleft.eq.2 ) then
                  Ajj1=A(jj+1)
                  Ajj2=A(jj+2)
                  Aii1=A(ii+1)
                  Aii2=A(ii+2)
                  A(ii+1)=SinA*Ajj1+CosA*Aii1
                  A(ii+2)=SinA*Ajj2+CosA*Aii2
                  A(jj+1)=CosA*Ajj1-SinA*Aii1
                  A(jj+2)=CosA*Ajj2-SinA*Aii2
               Else If ( kleft.eq.3 ) then
                  Ajj1=A(jj+1)
                  Ajj2=A(jj+2)
                  Ajj3=A(jj+3)
                  Aii1=A(ii+1)
                  Aii2=A(ii+2)
                  Aii3=A(ii+3)
                  A(ii+1)=SinA*Ajj1+CosA*Aii1
                  A(ii+2)=SinA*Ajj2+CosA*Aii2
                  A(ii+3)=SinA*Ajj3+CosA*Aii3
                  A(jj+1)=CosA*Ajj1-SinA*Aii1
                  A(jj+2)=CosA*Ajj2-SinA*Aii2
                  A(jj+3)=CosA*Ajj3-SinA*Aii3
               End If
               kmin=kmin+kleft
               Do k=kmin,kmax,4
                  Ajj0=A(jj+k+0)
                  Ajj1=A(jj+k+1)
                  Ajj2=A(jj+k+2)
                  Ajj3=A(jj+k+3)
                  Aii0=A(ii+k+0)
                  Aii1=A(ii+k+1)
                  Aii2=A(ii+k+2)
                  Aii3=A(ii+k+3)
                  A(ii+k+0)=SinA*Ajj0+CosA*Aii0
                  A(ii+k+1)=SinA*Ajj1+CosA*Aii1
                  A(ii+k+2)=SinA*Ajj2+CosA*Aii2
                  A(ii+k+3)=SinA*Ajj3+CosA*Aii3
                  A(jj+k+0)=CosA*Ajj0-SinA*Aii0
                  A(jj+k+1)=CosA*Ajj1-SinA*Aii1
                  A(jj+k+2)=CosA*Ajj2-SinA*Aii2
                  A(jj+k+3)=CosA*Ajj3-SinA*Aii3
               End Do
               kmin=j+1
               kmax=i-1
               kleft=Mod(kmax-kmin+1,4)
               kk=jj+j
               If ( kleft.eq.1 ) then
                  k0=kk
                  Ak0j=A(k0+j)
                  Aii0=A(ii+kmin+0)
                  A(ii+kmin+0)=SinA*Ak0j+CosA*Aii0
                  A(kk+j)=CosA*Ak0j-SinA*Aii0
                  kk=k0+kmin
               Else If ( kleft.eq.2 ) then
                  k0=kk
                  k1=k0+kmin
                  Ak0j=A(k0+j)
                  Ak1j=A(k1+j)
                  Aii0=A(ii+kmin+0)
                  Aii1=A(ii+kmin+1)
                  A(k0+j)=CosA*Ak0j-SinA*Aii0
                  A(k1+j)=CosA*Ak1j-SinA*Aii1
                  A(ii+kmin+0)=SinA*Ak0j+CosA*Aii0
                  A(ii+kmin+1)=SinA*Ak1j+CosA*Aii1
                  kk=k1+kmin+1
               Else If ( kleft.eq.3 ) then
                  k0=kk
                  k1=k0+kmin
                  k2=k1+kmin+1
                  Ak0j=A(k0+j)
                  Ak1j=A(k1+j)
                  Ak2j=A(k2+j)
                  Aii0=A(ii+kmin+0)
                  Aii1=A(ii+kmin+1)
                  Aii2=A(ii+kmin+2)
                  A(k0+j)=CosA*Ak0j-SinA*Aii0
                  A(k1+j)=CosA*Ak1j-SinA*Aii1
                  A(k2+j)=CosA*Ak2j-SinA*Aii2
                  A(ii+kmin+0)=SinA*Ak0j+CosA*Aii0
                  A(ii+kmin+1)=SinA*Ak1j+CosA*Aii1
                  A(ii+kmin+2)=SinA*Ak2j+CosA*Aii2
                  kk=k2+kmin+2
               End If
               kmin=kmin+kleft
               Do k=kmin,kmax,4
                  k0=kk
                  k1=k0+k
                  k2=k1+k+1
                  k3=k2+k+2
                  Ak0j=A(k0+j)
                  Ak1j=A(k1+j)
                  Ak2j=A(k2+j)
                  Ak3j=A(k3+j)
                  Aii0=A(ii+k+0)
                  Aii1=A(ii+k+1)
                  Aii2=A(ii+k+2)
                  Aii3=A(ii+k+3)
                  A(k0+j)=CosA*Ak0j-SinA*Aii0
                  A(k1+j)=CosA*Ak1j-SinA*Aii1
                  A(k2+j)=CosA*Ak2j-SinA*Aii2
                  A(k3+j)=CosA*Ak3j-SinA*Aii3
                  A(ii+k+0)=SinA*Ak0j+CosA*Aii0
                  A(ii+k+1)=SinA*Ak1j+CosA*Aii1
                  A(ii+k+2)=SinA*Ak2j+CosA*Aii2
                  A(ii+k+3)=SinA*Ak3j+CosA*Aii3
                  kk=k3+k+3
               End Do
               kmin=i+1
               kmax=N
               kleft=Mod(kmax-kmin+1,4)
               kk=ii+i
               If ( kleft.eq.1 ) then
                  k0=kk
                  Ak0j=A(k0+j)
                  Ak0i=A(k0+i)
                  A(k0+j)=CosA*Ak0j-SinA*Ak0i
                  A(k0+i)=SinA*Ak0j+CosA*Ak0i
                  kk=k0+kmin
               Else If ( kleft.eq.2 ) then
                  k0=kk
                  k1=k0+kmin
                  Ak0j=A(k0+j)
                  Ak0i=A(k0+i)
                  Ak1j=A(k1+j)
                  Ak1i=A(k1+i)
                  A(k0+j)=CosA*Ak0j-SinA*Ak0i
                  A(k0+i)=SinA*Ak0j+CosA*Ak0i
                  A(k1+j)=CosA*Ak1j-SinA*Ak1i
                  A(k1+i)=SinA*Ak1j+CosA*Ak1i
                  kk=k1+kmin+1
               Else If ( kleft.eq.3 ) then
                  k0=kk
                  k1=k0+kmin
                  k2=k1+kmin+1
                  Ak0j=A(k0+j)
                  Ak0i=A(k0+i)
                  Ak1j=A(k1+j)
                  Ak1i=A(k1+i)
                  Ak2j=A(k2+j)
                  Ak2i=A(k2+i)
                  A(k0+j)=CosA*Ak0j-SinA*Ak0i
                  A(k0+i)=SinA*Ak0j+CosA*Ak0i
                  A(k1+j)=CosA*Ak1j-SinA*Ak1i
                  A(k1+i)=SinA*Ak1j+CosA*Ak1i
                  A(k2+j)=CosA*Ak2j-SinA*Ak2i
                  A(k2+i)=SinA*Ak2j+CosA*Ak2i
                  kk=k2+kmin+2
               End If
               kmin=kmin+kleft
               Do k=kmin,kmax,4
                  k0=kk
                  k1=k0+k
                  k2=k1+k+1
                  k3=k2+k+2
                  Ak0j=A(k0+j)
                  Ak0i=A(k0+i)
                  Ak1j=A(k1+j)
                  Ak1i=A(k1+i)
                  Ak2j=A(k2+j)
                  Ak2i=A(k2+i)
                  Ak3j=A(k3+j)
                  Ak3i=A(k3+i)
                  A(k0+j)=CosA*Ak0j-SinA*Ak0i
                  A(k0+i)=SinA*Ak0j+CosA*Ak0i
                  A(k1+j)=CosA*Ak1j-SinA*Ak1i
                  A(k1+i)=SinA*Ak1j+CosA*Ak1i
                  A(k2+j)=CosA*Ak2j-SinA*Ak2i
                  A(k2+i)=SinA*Ak2j+CosA*Ak2i
                  A(k3+j)=CosA*Ak3j-SinA*Ak3i
                  A(k3+i)=SinA*Ak3j+CosA*Ak3i
                  kk=k3+k+3
               End Do
*---  update the diagonal elements of A -------------------------------*
               Temp=2.0*CosA*SinA*Aij
               CosA2=CosA*CosA
               SinA2=SinA*SinA
               A(jj+j)=SinA2*Aii+CosA2*Ajj-Temp
               A(ii+j)=0.0
               A(ii+i)=CosA2*Aii+SinA2*Ajj+Temp
*---  update rows/columnes of the eigenvectors C ----------------------*
               kmin=1  
               kmax=LdC
               kleft=Mod(kmax-kmin+1,4)
               If ( kleft.eq.1 ) then
                  C1j=C(1,j)
                  C1i=C(1,i)
                  C(1,i)=SinA*C1j+CosA*C1i
                  C(1,j)=CosA*C1j-SinA*C1i
               Else If ( kleft.eq.2 ) then
                  C1j=C(1,j)
                  C2j=C(2,j)
                  C1i=C(1,i)
                  C2i=C(2,i)
                  C(1,i)=SinA*C1j+CosA*C1i
                  C(2,i)=SinA*C2j+CosA*C2i
                  C(1,j)=CosA*C1j-SinA*C1i
                  C(2,j)=CosA*C2j-SinA*C2i
               Else If ( kleft.eq.3 ) then
                  C1j=C(1,j)
                  C2j=C(2,j)
                  C3j=C(3,j)
                  C1i=C(1,i)
                  C2i=C(2,i)
                  C3i=C(3,i)
                  C(1,i)=SinA*C1j+CosA*C1i
                  C(2,i)=SinA*C2j+CosA*C2i
                  C(3,i)=SinA*C3j+CosA*C3i
                  C(1,j)=CosA*C1j-SinA*C1i
                  C(2,j)=CosA*C2j-SinA*C2i
                  C(3,j)=CosA*C3j-SinA*C3i
               End If
               kmin=kmin+kleft
               Do k=kmin,kmax,4
                  C0j=C(k+0,j)
                  C1j=C(k+1,j)
                  C2j=C(k+2,j)
                  C3j=C(k+3,j)
                  C0i=C(k+0,i)
                  C1i=C(k+1,i)
                  C2i=C(k+2,i)
                  C3i=C(k+3,i)
                  C(k+0,i)=SinA*C0j+CosA*C0i
                  C(k+1,i)=SinA*C1j+CosA*C1i
                  C(k+2,i)=SinA*C2j+CosA*C2i
                  C(k+3,i)=SinA*C3j+CosA*C3i
                  C(k+0,j)=CosA*C0j-SinA*C0i
                  C(k+1,j)=CosA*C1j-SinA*C1i
                  C(k+2,j)=CosA*C2j-SinA*C2i
                  C(k+3,j)=CosA*C3j-SinA*C3i
               End Do
            End If
         End Do
      End Do
*---  find the value of the largest subdiagonal element ---------------*
      ASqrMax=0.0
      Do i=2,N
         ii=i*(i-1)/2
         Do j=1,i-1
            Temp=A(ii+j)*A(ii+j)
            ASqrMax=Max(ASqrMax,Temp)
         End Do
      End Do
*---  update the threshold --------------------------------------------*
      Thrs=0.1*ASqrMax
*---  check the accuracy ----------------------------------------------*
      If ( ASqrMax.gt.EpsSqr ) Goto 100
*----------------------------------------------------------------------*
*     Thats the end                                                    *
*----------------------------------------------------------------------*
200   Return
      End 
