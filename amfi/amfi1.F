     

      subroutine readbas(Lhigh,makemean,bonn,breit,
     *symmetry,sameorb,AIMP,oneonly,ncont4,numballcart,IN,ifinite,
     *EXP_FIN) 
cbs   suposed to read the maximum of l-values, the number of primitive and contracted 
cbs   functions, the exponents and contraction coefficients
#include <implicit.h>
#include <para.h>
#include <param.h>
#include <ired.h>
#include <priunit.h>
      character*4 WORD                         
      character*4 symmetry   
      character*13 Llimit 
      character*19 chcharge
      character*30 Nofprim
      character*28 addtext
      character*32 Nofcont
      character*76 Stars                  
      logical makemean,bonn,breit,
     *sameorb,AIMP,oneonly  
      common /nucleus/ charge,Exp_finite 
      Integer ibeginIRED(8),idelpersym(8) 
      dimension INOFT(Mxcart),INOFF(MxCart)
      stars='********************************************************'//
     * '*********************'
      Llimit='MAX. L-VALUE:'
      chcharge=' CHARGE OF NUCLEUS:'
      Nofprim='NUMBER OF PRIMITIVE FUNCTIONS:' 
      Nofcont=' NUMBER OF CONTRACTED FUNCTIONS:' 
      addtext='ADDITIONAL FUNCTIONS in IRS:'
CBS   write(LUPRI,*)                                                      
CBS   write(LUPRI,*) 'ATOMIC NO-PAIR SO-MF CODE starts'
CBS   write(LUPRI,*)                                                      
      bonn=.false.
      sameorb=.false.
      aimp=.false.
      oneonly=.false. 
      makemean=.true.
CBS   write(LUPRI,*) stars                                                
CBS   write(LUPRI,*) '2e-integrals for the mean-field only'
CBS   write(LUPRI,*) '    mean-field will be generated         ' 
CBS   write(LUPRI,*) stars                                                
        do i=0,Lmax 
        icore(i)=0
        enddo   
       if (ifinite.eq.1) Exp_finite=EXP_FIN
      if (BONN) then
CBS   write(LUPRI,*) 'Bonn-approach for spin-other-orbit part'
      endif 
      if (BREIT) then
CBS   write(LUPRI,*) ' Breit-Pauli-Approximation'
      else
CBS   write(LUPRI,*) 'Douglas-Kroll type operators '
      endif 
      if (ifinite.eq.0) then 
CBS   write(LUPRI,*) 'Point-nucleus '
      else
CBS   write(LUPRI,*) 'Finite Nucleus' 
      endif  
CBS   write(LUPRI,*) stars                                                
CBS   write(LUPRI,*) 'write out one-electron integrals in MOLCAS-style'
CBS   write(LUPRI,*) '   and with MOLCAS normalization ' 
CBS   write(LUPRI,*) stars                                                
CBS   write(LUPRI,*) stars                                                
CBS   write(LUPRI,*) 
      symmetry='D2H' 
CBS   write(LUPRI,*) 'Symmetry is D2H'
CBS   write(LUPRI,*) 'check whether order of IRs is correct!!!'
      numbofsym=8    
      if (SAMEORB) then 
CBS   write(LUPRI,*) 'SAME-ORBIT only'
      else
CBS   write(LUPRI,*) 'OTHER-ORBIT included' 
      endif 
      if (AIMP) then 
CBS   write(LUPRI,*) 'CORE removed for use with AIMP' 
      endif 
      read(IN,*) charge,Lhigh
      if (Lhigh.gt.Lmax) then 
      write(LUPRI,*) 'Sorry, so far this code deals only ',
     *'with maximum l-values of ',Lmax
      CALL QUIT('Too high angular momentum values in AMFI')            
      endif  
CBS   write(LUPRI,*) ' Functions will go up to an L-value of : ',Lhigh
CBS   write(LUPRI,'(A19,F5.2)') chcharge,charge    
      call initired
      Do iredrun=1,numbofsym
      do Lrun=0,Lhigh  
      nmbMperIRL(iredrun,Lrun)=0 
      enddo 
      enddo 
      do Lrun=0,Lhigh 
CBS   write(LUPRI,*) 'ANGULAR MOMENTUM ',LRUN  
         read(IN,*) nprimit(Lrun),ncontrac(Lrun) 
CBS   write(LUPRI,'(I3,I3)') nprimit(Lrun),ncontrac(Lrun) 
cbs   check keywords
cbs   check maximum numbers
         if (nprimit(Lrun).gt.MxprimL) then 
            write(LUPRI,*) 'To many primitives for L=',L,
     *           ' increase MxprimL in para.h or reduce ',
     *           'the number of primitives to at least ',MxprimL
            CALL QUIT('Too many primitive functions in AMFI')            
         endif   
         if (ncontrac(Lrun).gt.MxcontL) then 
            write(LUPRI,*) 'To many contracted fncts for L=',L,
     *           ' increase MxcontL in para.h or ',
     *           'reduce the number of contracted functions',
     *           'to at most ',MxcontL
            CALL QUIT('Too many contracted functions in AMFI')            
         endif   
         if (ncontrac(Lrun).gt.nprimit(Lrun)) then 
            write(LUPRI,*) 'You have more contracted than ',
     *           'uncontracted functions, I don''t believe ',
     *           'that. Sorry!! '
            CALL QUIT('Inconsistent input detected in AMFI')            
         endif
C     write(LUPRI,'(A7,I3,A15,I3,A33,I3,A24)') 'For L= ',Lrun,
C    *' there will be ',
C    *ncontrac(Lrun),' contracted functions, built from ',
C    *nprimit(Lrun),
C    *' uncontracted functions.' 
         do ILINE=1,nprimit(Lrun)
            read(IN,*) exponents(ILINE,Lrun), 
     *           (cntscrtch(ILINE,JRUN,Lrun),
     *           Jrun=1,ncontrac(Lrun))
         enddo 
ckr         read(IN,'(A76)') header
c
cbs   
cbs   end of reading for the current L-value 
cbs   
c     do  Irun=1,ncontrac(Lrun)
c     writE(LUPRI,*) 'orbital : ',irun   
c     write(LUPRI,'(6(X,E12.6))') 
c    *(cntscrtch(I,Irun,Lrun),I=1,nprimit(Lrun))
c     enddo 
c     write(LUPRI,*) ' ' 
cbs   setting the numbers of cartesians per IR
         do iredrun=1,numbofsym 
            nfunctions(iredrun,Lrun)=0
         enddo
         do mrun=-Lrun,Lrun
            nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),
     *           ipowxyz(2,mrun,Lrun),Ipowxyz(3,mrun,Lrun)),Lrun)=
     *           nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),
     *           ipowxyz(2,mrun,Lrun),
     *           ipowxyz(3,mrun,Lrun)),Lrun)+ncontrac(Lrun)
         enddo
         do mrun=-Lrun,Lrun
            nmbMperIRL(ipow2ired(ipowxyz(1,mrun,Lrun),
     *           ipowxyz(2,mrun,Lrun),Ipowxyz(3,mrun,Lrun)),lruN)=
     *           nmbMperIRL(ipOw2ired(ipowxyz(1,mrun,Lrun),
     *           ipowxyz(2,mrun,Lrun),IpowxYz(3,mrun,Lrun)),lruN)+1
         enddo
CBS   write(LUPRI,*) stars                                                       
CBS   write(LUPRI,'(A,8I4)') 
CBS  *'Number of functions per IR: ',(nfunctions(iredrun,Lrun),
CBS  *iredrun=1,numbofsym)
CBS   write(LUPRI,*) stars                                                          
      enddo                     ! enddo for loop over L-values 
C     write(LUPRI,*) 'distribution of M-values'
c     do Lrun=0,Lhigh
c     write(LUPRI,*) (nmbMperIRL(nsym,Lrun),nsym=1,numbofsym)
c     endDo 
      numbofcart=0
      do lrun=0,Lhigh 
      numbofcart=numbofcart+(Lrun+Lrun+1)*
     *ncontrac(Lrun)
      enddo
      do iredrun=1,numbofsym 
      nfunctperIRED(iredrun)=0
      enddo
      do Lrun=0,Lhigh
      do iredrun=1,numbofsym 
      nfunctperIRED(iredrun)=nfunctperIRED(iredrun)+
     *nfunctions(iredrun,Lrun)
      enddo
      enddo
CBS   write(LUPRI,*) stars                                                        
CBS   write(LUPRI,'(A,8I3)') 'total number of atomic functions per IRED ',
CBS  *(nfunctperIRED(iredrun),iredrun=1,numbofsym) 
CBS   write(LUPRI,*) stars                                                        
      isum=0
      do iredrun=1,numbofsym
      itotalperIR(iredrun)=nfunctperIRED(iredrun)
      isum=isum+itotalperIR(iredrun)
      enddo 
      numballcart=isum
      iorbrun=0
      do iredrun=1,numbofsym 
      do inired=1,itotalperIR(iredrun)
      iorbrun=iorbrun+1 
      IREDoffunctnew(Iorbrun)=iredrun 
      enddo 
      enddo 
CBS   write(LUPRI,*) stars                                                        
CBS   write(LUPRI,'(A,8I3)') 'including additional functions per IRED ',
CBS  *(itotalperIR(iredrun),iredrun=1,numbofsym) 
CBS   write(LUPRI,*) stars                                                          
      do iredrun=1,numbofsym  
      ibeginIRED(iredrun)=0                              
      enddo
      do lrun=0,Lhigh
      do mrun=-lrun,lrun
      iredLM(mrun,lrun)=ipow2ired(ipowxyz(1,mrun,Lrun),
     *ipowxyz(2,mrun,Lrun),
     *ipowxyz(3,mrun,Lrun))
      incrLM(mrun,lrun)=ibeginIRED(iredLM(mrun,lrun))
      ibeginIRED(iredLM(mrun,lrUn))=
     *ibeginIRED(iredLM(mrun,lrun))+ncontrac(lrun) 
      enddo
      enddo   
c     do lrun=0,Lhigh
c     write(LUPRI,'(A,I4,A,21I3)') 'L= ',lrun,
c    *' shifts inside the IRED',
c    *(incrLM(mrun,lrun),mrun=-lrun,lrun)
c     enddo
      shiftIRED(1)=0
      do iredrun=2,numbofsym  
      shiftIRED(iredrun)=shiftIRED(iredrun-1)
     *                   +itotalperIR(iredrun-1)
      enddo
c     write(LUPRI,'(A,8I4)') 'shifts for the IREDs ',
c    *(shiftIRED(iredrun),iredrun=1,numbofsym) 
cbs   test all orbital numbers
c     do lrun=0,Lhigh
c     do mrun=-Lrun,Lrun
c     do irun=1,ncontrac(lrun)
c     write(LUPRI,*) 'L,M,contr funct, absolute number ',
c    *lrun,mrun,irun,shiftired(iredLM(mrun,lrun))+
c    *incrLM(mrun,Lrun)+irun
c     enddo
c     enddo
c     enddo
      shiftIRIR(1)=0
      irun=1
      do ired1=2,numbofsym 
      do ired2=1,ired1 
      irun=irun+1
      if (ired2.eq.1) then
      shiftIRIR(irun)=shiftIRIR(irun-1)+
     *(itotalperIR(ired1-1)*itotalperIR(ired1-1)+
     *itotalperIR(ired1-1))/2
      else
      shiftIRIR(irun)=shiftIRIR(irun-1)+
     *itotalperIR(ired1)*itotalperIR(ired2-1)
      endif 
c     write(LUPRI,*) 'ired1,ired2 ',ired1,ired2,
c    *irun,shiftIRIR(irun)
      enddo
      enddo
cbs  
      do lrun=0,Lhigh
      do Mrun=-Lrun,Lrun
      ired=iredLM(Mrun,Lrun)
      ishifter=shiftIRED(ired)+incrLM(mrun,lrun)      
      do icart=1,ncontrac(Lrun)
      moffunction(ishifter+icart)=Mrun
      Loffunction(ishifter+icart)=Lrun
      IREDoffunction(ishifter+Icart)=ired  
      INOFT(ishifter+Icart)=icart     
      enddo
      enddo
      enddo
CBS   write(LUPRI,*) stars
CBS   write(LUPRI,*) 'SYMMETRY-INFORMATION ON FUNCTIONS '
CBS   write(LUPRI,*) stars
      do irun = 1, numbofcart
CBS   write(LUPRI,'(4(A,I3))') 'Number of function: ',
CBS  *irun,
CBS  *' IR of function: ',IREDoffunction(irun),
CBS  *' L-value: ',Loffunction(irun),
CBS  *' M-value: ',Moffunction(irun)
CBS   numboffunct(irun)=irun
      INOFF(irun)=irun
CBS   if (IREDoffunction(irun).ne.IREDoffunction(irun+1)) 
CBS  *write(LUPRI,*)   
      enddo
      do nsymrun=1,numbofsym
      idelpersym(nsymrun)=0
      enddo 
      do nsymrun=1,numbofsym
      nrtofiperIR(nsymrun)=itotalperIR(nsymrun) 
      enddo
      if (AIMP) then 
cbs   generate list of orbitals to be removed
      ikeeporb=0
      numbprev=0
      do irun=1,numbofcart
4712  if (irun.eq.1.or.(irun.ge.2.and.INOFF(irun).eq.  
     *numbprev+1)) then 
      Lval=Loffunction(irun)
      number=INOFF(irun)
      itype=INOFT(irun)  
      if (itype.le.icore(lval)) then
      write(LUPRI,777) number,itype,lval
      idelpersym(IREDoffunction(irun))=
     *               idelpersym(IREDoffunction(irun))+1
      numbprev=number
      else 
      ikeeporb=ikeeporb+1
      ikeeplist(ikeeporb)=number
      numbprev=number
      endif 
      else
      ikeeporb=ikeeporb+1
      ikeeplist(ikeeporb)=numbprev+1
      numbprev=numbprev+1     
      goto 4712
      endif 
      enddo 
      ikeeporb=0  
      do nsymrun=1,numbofsym
      nrtofiperIR(nsymrun)=itotalperIR(nsymrun)-idelpersym(nsymrun)  
      enddo
      do nsymrun=1,numbofsym
      ikeeporb=ikeeporb+nrtofiperIR(nsymrun)                        
      enddo
CBS   write(LUPRI,*) stars                                                          
      write(LUPRI,'(A,8I3)')'# of funct. per IRED after removing core ',
     *(nrtofiperIR(iredrun),iredrun=1,numbofsym) 
      write(LUPRI,*) ikeeporb,' orbitals left after deleting core' 
      endif 
CBS   write(LUPRI,*) stars
      nmax=max(6,ncontrac(0))  
      do lrun=1,Lhigh
      nmax=max(nmax,ncontrac(lrun))
      enddo
      ncont4=nmax*nmax*nmax*nmax
      return 
777   format('ORBITAL NUMBER ',I4,' IS THE ',I3,'TH of L-value ',I3,
     *' IT WILL BE REMOVED !!!')  
      end 
      double precision function  regge3j(
     *j1,     ! integer  2*j1
     *j2,     ! integer  2*j2
     *j3,     ! integer  2*j3
     *m1,     ! integer  2*m1
     *m2,     ! integer  2*m2
     *m3)     ! integer  2*m3
cbs   uses magic square of regge (see Lindner pp. 38-39)
cbs 
cbs    ---                                            ---
cbs   |                                                  |
cbs   | -j1+j2+j3     j1-j2+j3         j1+j2-j3          |
cbs   |                                                  |
cbs   |                                                  |
cbs   |  j1-m1        j2-m2            j3-m3             |
cbs   |                                                  |
cbs   |                                                  |
cbs   |  j1+m1        j2+m2            j3+m3             |
cbs   |                                                  |
cbs    ---                                            ---
cbs 
#include <implicit.h>
      dimension MAT(3,3)
      logical testup,testdown
#include <Regge.h>
cbs  facul,   integer array (nprim,0:mxLinRE) prime-expansion of factorials 
cbs  mxLinRE,    integer max. number for facul is given
cbs  nprim,   number of primes for expansion of factorials 
cbs  prim,    integer array with the first nprim prime numbers
cbs  iwork)   integer array of size nprim
      regge3j=0d0
c     write(6,'(A24,6I3)') '3J to be calculated for ',
c    *j1,j2,j3,m1,m2,m3
cbs   quick check  if =/= 0 at all
      icheck=m1+m2+m3 
      if (icheck.ne.0) then 
c     write(6,*) 'sum over m =/= 0'
      return    
      endif 
cbs   check triangular relation (|j1-j2|<= j3 <= j1+j2 )
      imini=iabs(j1-j2)
      imaxi=j1+j2 
      if (j3.lt.imini.or.j3.gt.imaxi) then 
c     write(6,*) 'triangular relation not fulfilled'
      return   
      endif
cbs   quick check  if =/= 0 at all  end 
cbs  
cbs   3J-symbol is not zero by simple rules 
cbs  
cbs   initialize MAT 
      MAT(1,1) =-j1+j2+j3
      MAT(2,1) =j1-m1    
      MAT(3,1) =j1+m1     
      MAT(1,2) =j1-j2+j3  
      MAT(2,2) =j2-m2      
      MAT(3,2) =j2+m2     
      MAT(1,3) =j1+j2-j3  
      MAT(2,3) =j3-m3      
      MAT(3,3) =j3+m3      
      do I=1,3
      do J=1,3
cbs   check for even numbers (2*integer) and positive or zero
      if (mod(MAT(J,I),2).ne.0.or.MAT(J,I).lt.0)  then 
c     write(6,*) 'J,I,MAT(J,I): ',J,I,MAT(J,I)
      return
      endif 
      MAT(J,I)=MAT(J,I)/2
      if (Mat(j,i).gt.mxLinRE)
     *CALL QUIT('increase mxLinRE for regge3j')
      enddo
      enddo
      Isigma=(j1+j2+j3)/2
cbs   check the magic sums
      do I=1,3
      IROW=0
      ICOL=0
      do J=1,3
      IROW=IROW+MAT(I,J)
      ICOL=ICOL+MAT(J,I)
      enddo
      if (IROW.ne.Isigma.or.ICOL.ne.Isigma) then 
c     write(6,*) 'I,IROW,ICOL ',I,IROW,ICOL  
      return
      endif
      enddo
cbs   if j1+j2+j3 is odd: check for equal rows or columns 
      Isign=1
      if (iabs(mod(Isigma,2)).eq.1) then 
      isign=-1
         do I=1,3
         do J=I+1,3
            if (MAT(1,I).eq.MAT(1,J).and. 
     *         MAT(2,I).eq.MAT(2,J).and.
     *         MAT(3,I).eq.MAT(3,J)) return 
            if (MAT(I,1).eq.MAT(J,1).and. 
     *         MAT(I,2).eq.MAT(J,2).and.
     *         MAT(I,3).eq.MAT(J,3)) return 
         enddo
         enddo
      endif 
cbs   look for the lowest element indices: IFIRST,ISECOND
      imini=MAT(1,1) 
      IFIRST=1
      ISECOND=1
      do I=1,3
      do J=1,3 
      if (MAT(J,I).lt.imini) then 
      IFIRST=J
      ISECOND=I
      imini=MAT(J,I)
      endif 
      enddo
      enddo
c     write(6,*) 'Matrix before commuting vectors'
      do ibm=1,3
c     write(6,'(3I5)') (Mat(ibm,j),j=1,3) 
      enddo
      if (IFIRST.ne.1) then  !interchange rows
c     write(6,*) 'IFIRST = ',ifirst
      do I=1,3
      IDUMMY=MAT(1,I) 
      MAT(1,I)=MAT(IFIRST,I)
      MAT(IFIRST,I)=IDUMMY
      enddo
      endif 
      if (ISECOND.ne.1) then  !interchange columns
c     write(6,*) 'ISECOND = ',isecond
      do I=1,3
      IDUMMY=MAT(I,1) 
      MAT(I,1)=MAT(I,ISECOND)
      MAT(I,ISECOND)=IDUMMY
      enddo
      endif 
cbs   lowest element is now on (1,1)
c     write(6,*) 'Matrix after commuting vectors'
c     do ibm=1,3
c     write(6,'(3I5)') (Mat(ibm,j),j=1,3) 
c     enddo
cbs   begin to calculate Sum over s_n
cbs   first the simple cases
      if (Mat(1,1).eq.0) then 
      isum=1 
      elseif (Mat(1,1).eq.1) then 
      isum=Mat(2,3)*Mat(3,2)-Mat(2,2)*Mat(3,3)
      elseif (Mat(1,1).eq.2) then 
      isum=Mat(2,3)*(Mat(2,3)-1)*Mat(3,2)*(Mat(3,2)-1)-
     *2*Mat(2,3)*Mat(3,2)*Mat(2,2)*Mat(3,3)+
     *Mat(2,2)*(Mat(2,2)-1)*Mat(3,3)*(Mat(3,3)-1)
      else !  all the cases with Mat(1,1) >= 3 
	      Icoeff=1
	      do Ibm=Mat(3,2)-Mat(1,1)+1,Mat(3,2)
	        icoeff=icoeff*ibm
	      enddo
	      do Ibm=Mat(2,3)-Mat(1,1)+1,Mat(2,3)
	        icoeff=icoeff*ibm
	      enddo
	      isum=icoeff
	      do Icount=1,MAT(1,1) 
	         icoeff=-icoeff*(Mat(1,1)+1-icount)*(Mat(2,2)+1-icount)*
     *           (Mat(3,3)+1-icount)
	         Idenom=icount*(Mat(2,3)-Mat(1,1)+icount)*
     *           (Mat(3,2)-Mat(1,1)+icount)
	         icoeff=icoeff/Idenom
	         isum=isum+icoeff
	      enddo
      endif   
cbs  additional sign from interchanging rows or columns
      if (ifirst.ne.1) isum=isum*isign
      if (isecond.ne.1) isum=isum*isign
c     write(6,*) 'isum = ',isum 
cbs       Mat(2,3)+Mat(3,2) 
cbs    (-) 
      if (iabs(mod((Mat(2,3)+Mat(3,2)),2)).eq.1) isum=-isum   
cbs   final factor
      LIMIT=ihigh(max(Mat(1,1),Mat(1,2),Mat(1,3),
     *Mat(2,1),Mat(2,2),Mat(2,3),Mat(3,1),Mat(3,2),
     *Mat(3,3),(Isigma+1)))
      do I=1,LIMIT 
      iwork(I)=facul(I,Mat(1,2))+facul(I,Mat(2,1))+
     *facul(I,Mat(3,1))+facul(I,Mat(1,3))-
     *facul(I,Mat(1,1))-facul(I,Mat(2,2))-
     *facul(I,Mat(3,3))-facul(I,(Isigma+1))-
     *facul(I,Mat(2,3))-facul(I,Mat(3,2))
      enddo
c     write(6,*) 'Iwork: ',(iwork(i),i=1,LIMIT)
      factor=1d0
      iup=1
      idown=1
      testup=.true.
      testdown=.true.
      do I=1,LIMIT   
      do J=1,iwork(I)
      iup=iup*prim(i)
      if (iup.lt.0) testup=.false. !check for Integer overflow
      enddo 
      Enddo 
      up=dfloat(iup) 
      if(.not.testup) then ! if the integers did not run correctly  
	      up=1d0
	      do I=1,LIMIT   
              do J=1,iwork(I)
              up=up*dfloat(prim(i))
              enddo 
	      enddo 
      endif 
      do I=1,LIMIT   
      do J=1,-iwork(I)
      idown=idown*prim(i)
      if (idown.lt.0) testdown=.false. 
      enddo 
      enddo 
      down=dfloat(idown) 
      if(.not.testdown) then 
	      down=1d0
	      do I=1,LIMIT   
              do J=1,-iwork(I)
              down=down*dfloat(prim(i))
              enddo 
	      enddo 
      endif 
c     if (.not.(testup.and.testdown)) then 
c     write(6,*) 'j1,j2,j3,m1,m2,m3 ',j1,j2,j3,m1,m2,m3
c     write(6,*) 'iup,idown ',iup,idown,'up,down ',up,down 
c     endif 
      factor=factor*up/down
cbs   final result
      regge3j=dsqrt(factor)*dfloat(isum)
      return 
      end 


      double precision function Tkinet(l,alpha1,alpha2)
cbs   calculates the matrix element of kinetic energy 
cbs   for primitive normalized functions with the same angular momentum l 
cbs   and exponents alpha1 and alpha2 
cbs   works only, if r**l is assumed for an l-value
cbs   formular obtained from the symmetric expression (d/dr's to (')
cbs   the left and to the right. 
cbs   Overlaps of the different powers are partially crossed out 
cbs   with  the overlap of functions with angular momentum l
cbs   final formula:
cbs   Tkinet=0.5*alpha12 (2l+3) (alpha1*alpha2/alpha12*alpha12)**((2L+7)/4)
cbs   with alpha12=0.5*(alpha1+alpha2)
cbs   as alpha12 has the dimensions 1/length**2, this can not be that bad...
      Implicit double precision (a-h,o-z)
Cbs   alpha12 is the effective exponent 
      Alpha12=0.5d0*(alpha1+alpha2)
      alphpro=alpha1*alpha2   
      ll3=l+l+3        
      ll7=l+l+7
      Tkinet=0.5d0*alpha12*ll3*(alphpro/
     *(alpha12*alpha12))**(0.25*dfloat(ll7))
      return 
      end    

CMI .......
      subroutine   tosigX(m1,m2,m3,m4,angint,
     *mcombina,ncontl1,ncontl2,ncontl3,
     *ncontl4,carteX,preXZ,interxyz,isgnprod,
     *cleaner) 
cbs   this subroutine combines the angular integrals 
cbs   to the integrals for the real-valued linear 
cbs   combinations for the sigma_X part 
cbs   definition of the real-valued linear combinations:
cbs
cbs
cbs   M=0  is the same as   Y(L,0)
cbs
cbs
cbs   M > 0
cbs   
cbs   | L,M,+> = 2**(-0.5) ( (-1)**M Y(L,M) + Y(L,-M)) 
cbs   
cbs   | L,M,-> = -i 2**(-0.5) ( (-1)**M Y(L,M) - Y(L,-M))  ($$$$) 
cbs
cbs
cbs   due to symmetry, there can be only integrals 
cbs   with indices one or three  (sigma_+ and sigma_-)- combinations 
cbs
#include <implicit.h>
#include <para.h>
#include <priunit.h>
      logical cleaner
      dimension mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *angint(ncontl1,ncontl2,ncontl3,ncontl4,*),
cbs  !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
     *carteX(ncontl1,ncontl3,ncontl2,ncontl4),
     *preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *interxyz(*),
     *isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),    
     *isgnM(-1:1,-1:1,-1:1,-1:1)
c     write(6,*) ' begin tosigx' 
cbs   cleaning up the integral-array
      irun=ncontl1*ncontl2*ncontl3*ncontl4 
      call dzero(cartex,irun)
cbs   set some signs
cbs   isgnM will give an additonal minus-sign if both m-values   
cbs   (cartesian and angular) are negative  see $$$$
      do irun4=-1,1
      do irun3=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,irun3,irun4)=1
      enddo
      enddo
      enddo
      enddo
      if (m1.lt.0) then  
      do irun4=-1,1
      do irun3=-1,1
      do irun2=-1,1
      isgnM(-1,irun2,irun3,irun4)=
     *-isgnM(-1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      endif 
      if (m2.lt.0) then
      do irun4=-1,1
      do irun3=-1,1
      do irun1=-1,1
      isgnM(irun1,-1,irun3,irun4)=
     *-isgnM(irun1,-1,irun3,irun4)
      enddo
      enddo
      enddo
      endif
      if (m3.lt.0) then 
      do irun4=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,-1,irun4)=
     *-isgnM(irun1,irun2,-1,irun4)
      enddo
      enddo
      enddo
      endif
      if (m4.lt.0) then 
      do irun3=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,irun3,-1)=
     *-isgnM(irun1,irun2,irun3,-1)
      enddo
      enddo
      enddo
      endif
cbs   define absolute m-values
      Mabs1=iabs(m1)
      Mabs2=iabs(m2)
      Mabs3=iabs(m3)
      Mabs4=iabs(m4)
      irun=0
      if (interxyz(1).eq.0) then 
      write(LUPRI,*) 'tosigx: no interaction: ',m1,m2,m3,m4
      CALL QUIT('Error in TOSIGX in AMFI')            
      endif 
      prexz1234=preXZ(m1,m2,m3,m4)
      do while (interxyz(irun+1).gt.0) 
      irun=irun+1
c     write(6,*) 'tosigx: ',irun,interxyz(irun)
c
cbs
cbs
cbs   This could be done with gotos, but I am biased to hate those..
cbs
cbs
         if (interxyz(irun).eq.1) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,1,1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.2) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,-1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.3) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,1,1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.4) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,-1,-1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.5) then  
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,1,-1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.6) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(-1,-1,1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.7) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,-1,1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.8) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.9) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,1,1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.10) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,-1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.11) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.12) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,-1,1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.13) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,-1,1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.14) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,1,-1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.15) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,-1,-1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.16) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(-1,1,1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
         endif  
       enddo 
        if (cleaner) then   
        do irun4=1,ncontl4
        do irun2=1,ncontl2
        do irun1=1,ncontl1
        cartex(irun1,irun1,irun2,irun4)=0d0
        enddo
        enddo
        enddo
        endif 
      return 
      end 
      subroutine   tosigY(m1,m2,m3,m4,angint,
     *mcombina,ncontl1,ncontl2,ncontl3,
     *ncontl4,carteY,preY,interxyz,isgnprod,
     *cleaner) 
cbs   this subroutine combines the angular integrals 
cbs   to the integrals for the real-valued linear 
cbs   combinations for the sigma_X part 
cbs   definition of the real-valued linear combinations:
cbs
cbs
cbs   M=0  is the same as   Y(L,0)
cbs
cbs
cbs   M > 0
cbs   
cbs   | L,M,+> = 2**(-0.5) ( (-1)**M Y(L,M) + Y(L,-M)) 
cbs   
cbs   | L,M,-> = -i 2**(-0.5) ( (-1)**M Y(L,M) - Y(L,-M)) ($$$$) 
cbs
cbs
cbs   due to symmetry, there can be only integrals 
cbs   with one or three (sigma_+ and sigma_-)  - combinations 
cbs
#include <implicit.h>
#include <priunit.h>
#include <para.h>
      dimension mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *angint(ncontl1,ncontl2,ncontl3,ncontl4,*),
cbs  !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
     *carteY(ncontl1,ncontl3,ncontl2,ncontl4),
     *preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *interxyz(*),
     *isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *isgnM(-1:1,-1:1,-1:1,-1:1)
      logical cleaner
c     write(6,*) 'begin tosigy '
cbs   cleaning up the integral-array
      irun=ncontl4*ncontl2*ncontl3*ncontl1
      call dzero(carteY,irun)
cbs   set some signs
cbs   isgnM will give an additonal minus-sign if both m-values   
cbs   (cartesian and angular) are negative  see $$$$
      do irun4=-1,1
      do irun3=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,irun3,irun4)=1
      enddo
      enddo
      enddo
      enddo
      if (m1.lt.0) then  
      do irun4=-1,1
      do irun3=-1,1
      do irun2=-1,1
      isgnM(-1,irun2,irun3,irun4)=
     *-isgnM(-1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      endif 
      if (m2.lt.0) then
      do irun4=-1,1
      do irun3=-1,1
      do irun1=-1,1
      isgnM(irun1,-1,irun3,irun4)=
     *-isgnM(irun1,-1,irun3,irun4)
      enddo
      enddo
      enddo
      endif
      if (m3.lt.0) then 
      do irun4=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,-1,irun4)=
     *-isgnM(irun1,irun2,-1,irun4)
      enddo
      enddo
      enddo
      endif
      if (m4.lt.0) then 
      do irun3=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,irun3,-1)=
     *-isgnM(irun1,irun2,irun3,-1)
      enddo
      enddo
      enddo
      endif
cbs   define absolute m-values
      Mabs1=iabs(m1)
      Mabs2=iabs(m2)
      Mabs3=iabs(m3)
      Mabs4=iabs(m4)
      irun=0
      if (interxyz(1).eq.0) then 
      write(LUPRI,*) 'tosigy: no interaction: ',m1,m2,m3,m4
      CALL QUIT('Error in TOSIGY in AMFI')
      endif
      prey1234=preY(m1,m2,m3,m4)
c     write(6,*) 'prey ',prey1234
      do while (interxyz(irun+1).gt.0) 
      irun=irun+1
c     write(6,*) 'tosigy: ',irun,interxyz(irun)
c
cbs
cbs
cbs   This could be done with gotos, but I am biased to hate those..
cbs
cbs
         if (interxyz(irun).eq.1) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,1,1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.2) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,-1,-1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.3) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,1,1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.4) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,-1,-1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.5) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,1,-1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.6) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(-1,-1,1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.7) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,-1,1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.8) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,1,-1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.9) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,1,1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.10) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,-1,-1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.11) then  
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,1,-1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.12) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,-1,1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.13) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,-1,1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.14) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,1,-1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.15) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,-1,-1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.16) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(-1,1,1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         endif  
      Enddo 
        if (cleaner) then   
        do irun4=1,ncontl4
        do irun2=1,ncontl2
        do irun1=1,ncontl1
        cartey(irun1,irun1,irun2,irun4)=0d0
        enddo
        enddo
        enddo
        endif
      return 
      end 
      subroutine   tosigZ(m1,m2,m3,m4,angint,
     *mcombina,ncontl1,ncontl2,ncontl3,
     *ncontl4,carteZ,preXZ,interxyz,isgnprod,
     *cleaner) 
cbs   this subroutine combines the angular integrals 
cbs   to the integrals for the real-valued linear 
cbs   combinations for the sigma_Z part 
cbs   definition of the real-valued linear combinations:
cbs
cbs
cbs   M=0  is the same as   Y(L,0)
cbs
cbs
cbs   M > 0
cbs   
cbs   | L,M,+> = 2**(-0.5) ( (-1)**M Y(L,M) + Y(L,-M)) 
cbs   
cbs   | L,M,-> = -i 2**(-0.5) ( (-1)**M Y(L,M) - Y(L,-M)) ($$$$) 
cbs
cbs   only angular integrals of type 2 (sigma_0) contribute  
cbs
#include <implicit.h>
#include <priunit.h>
#include <para.h>
      dimension mcombina(2,-Lmax:Lmax,-Lmax:Lmax,
     *-Lmax:Lmax,-Lmax:Lmax),
     *angint(ncontl1,ncontl2,ncontl3,ncontl4,*),
cbs  !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
     *carteZ(ncontl1,ncontl3,ncontl2,ncontl4),
     *preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *interxyz(*),
     *isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *isgnM(-1:1,-1:1,-1:1,-1:1)                          
      logical cleaner
cbs   cleaning up the integral-array
      irun=ncontl4*ncontl2*ncontl3*ncontl1
      call dzero(carteZ,irun)
c     write(6,*) 'begin tosigz'
cbs   set some signs 
cbs   isgnM will give an additonal minus-sign if both m-values   
cbs   (cartesian and angular) are negative  see $$$$
      do irun4=-1,1
      do irun3=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,irun3,irun4)=1
      enddo
      enddo
      enddo
      enddo
      if (m1.lt.0) then 
      do irun4=-1,1
      do irun3=-1,1
      do irun2=-1,1
      isgnM(-1,irun2,irun3,irun4)=
     *-isgnM(-1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      endif 
      if (m2.lt.0) then 
      do irun4=-1,1
      do irun3=-1,1
      do irun1=-1,1
      isgnM(irun1,-1,irun3,irun4)=
     *-isgnM(irun1,-1,irun3,irun4)
      enddo
      enddo
      enddo
      endif 
      if (m3.lt.0) then 
      do irun4=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,-1,irun4)=
     *-isgnM(irun1,irun2,-1,irun4)
      enddo
      enddo
      enddo
      endif 
      if (m4.lt.0) then 
      do irun3=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,irun3,-1)=
     *-isgnM(irun1,irun2,irun3,-1)
      enddo
      enddo
      enddo
      endif 
cbs   define absolute m-values
      Mabs1=iabs(m1)
      Mabs2=iabs(m2)
      Mabs3=iabs(m3)
      Mabs4=iabs(m4)
      irun=0
      if (interxyz(1).eq.0) then 
      write(LUPRI,*) 'tosigz: no interaction: ',m1,m2,m3,m4
      CALL QUIT('Error in TOSIGZ in AMFI')
      endif
      prexz1234=preXZ(m1,m2,m3,m4)
      do while (interxyz(irun+1).gt.0) 
      irun=irun+1
c
cbs
cbs
cbs   This could be done with gotos, but I am biased to hate those..
cbs
cbs
         if (interxyz(irun).eq.1) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,1,1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.2) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
         factor=-isgnM(-1,-1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.3) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,1,1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.4) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=-isgnM(-1,-1,-1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.5) then  
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,1,-1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.6) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=-isgnM(-1,-1,1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.7) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,-1,1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.8) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=-isgnM(-1,1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.9) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=-isgnM(-1,1,1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.10) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,-1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.11) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.12) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=-isgnM(-1,-1,1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.13) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,-1,1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.14) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=-isgnM(-1,1,-1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.15) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,-1,-1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.16) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=-isgnM(-1,1,1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         endif  
      enddo 
        if (cleaner) then   
        do irun4=1,ncontl4
        do irun2=1,ncontl2
        do irun1=1,ncontl1
        cartez(irun1,irun1,irun2,irun4)=0d0
        enddo
        enddo
        enddo
        endif
      return 
      end 
      subroutine trans(
cbs   makes the transformation for the ich-th index
     *coeffs, !(nolds(ith),nnew(ith)) modified contraction coefficients
     *idim1,  !  first dimension   
     *idim2,  !  second dimension  
     *ich,    ! index to be changed                         
     *nolds1,nolds2,nolds3,nolds4,  ! old dimensions                    
     *nnew1,nnew2,nnew3,nnew4, ! new dimensions      
     *array1, ! array of size (nolds1,nolds2,nolds3,nolds4)
     *array2  ! array of size (nnew1,nnew2,nnew3,nnew4)
     *)
#include <implicit.h>
      dimension coeffs(idim1,idim2),           
     *array1(nolds1,nolds2,nolds3,nolds4),
     *array2(nnew1,nnew2,nnew3,nnew4)
c     write(6,*) 'begin trans ' ,ich 
c     write(6,'(8I5)') nolds1,nolds2,nolds3,nolds4,
c    *nnew1,nnew2,nnew3,nnew4
      do ind4=1,nnew4
      do ind3=1,nnew3
      do ind2=1,nnew2
      do ind1=1,nnew1
      array2(ind1,ind2,ind3,ind4)=0d0
      enddo
      enddo
      enddo
      enddo
      if (ich.eq.1) then 
      do ind4=1,nnew4
      do ind3=1,nnew3
      do ind2=1,nnew2
      do ind5=1,nnew1  
      do ind1=1,nolds1
      array2(ind5,ind2,ind3,ind4)=array2(ind5,ind2,ind3,ind4)+
     *coeffs(ind1,ind5)*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
      elseif (ich.eq.2) then 
c     write(6,*) 'transform second index '
      do ind4=1,nnew4
      do ind3=1,nnew3
      do ind5=1,nnew2  
      do ind2=1,nolds2
      coeff=coeffs(ind2,ind5)
      do ind1=1,nnew1
      array2(ind1,ind5,ind3,ind4)=array2(ind1,ind5,ind3,ind4)+
     *coeff*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
c     write(6,*) 'end  to transform second index '
      elseif (ich.eq.3) then 
      do ind4=1,nnew4
      do ind5=1,nnew3  
      do ind3=1,nolds3
      coeff=coeffs(ind3,ind5)
      do ind2=1,nnew2
      do ind1=1,nnew1
      array2(ind1,ind2,ind5,ind4)=array2(ind1,ind2,ind5,ind4)+
     *coeff*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
      elseif (ich.eq.4) then 
      do ind5=1,nnew4  
      do ind4=1,nolds4
      coeff=coeffs(ind4,ind5)
      do ind3=1,nnew3
      do ind2=1,nnew2
      do ind1=1,nnew1
      array2(ind1,ind2,ind3,ind5)=array2(ind1,ind2,ind3,ind5)+
     *coeff*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
      endif  
c     write(6,*) 'end  trans ' 
      return
      end
      subroutine transcon(contold,idim1,idim2,ovlp,contnew,nprim,ncont)
#include <implicit.h>
      dimension contold(idim1,idim2),contnew(nprim,ncont),
     *ovlp(idim1,idim1)
c     write(6,*) 'begin transcon nprim,ncont ',nprim,ncont 
cbs   copy old contraction coefficients in dense form to common block 
      do Jrun=1,ncont
      do Irun=1,nprim
      contnew(Irun,Jrun)=contold(Irun,Jrun)
      enddo 
      enddo 
cbs   ensure normalization 
      do ICONT=1,ncont
      	xnorm=0d0   
      	do Jrun=1,nprim
      	do Irun=1,nprim
      	xnorm=xnorm+contnew(Irun,ICONT)*contnew(Jrun,ICONT)
     *	*ovlp(Irun,Jrun)
c       write(6,*) 'Icont,jrun,irun,xnorm ',
c    *  icont,jrun,irun,xnorm 
      	enddo 
      	enddo 
c       write(6,*) 'ICONT ',ICONT,xnorm   
      	xnorm=1d0/dsqrt(xnorm)
cbs   scale with normalization factor 
      	do Irun=1,nprim
      	contnew(Irun,ICONT)=xnorm*contnew(Irun,ICONT)
      	enddo 
      enddo 
c     write(6,*) 'end transcon nprim,ncont ',nprim,ncont 
      return 
      end 
      subroutine two2mean12a(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
#include <implicit.h>
#include <para.h>
      logical sameorb 
      dimension 
     *carteSO(ncontmf,norbsum,ncontmf,norbsum),
     *carteOO(ncontmf,norbsum,ncontmf,norbsum),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) THEN 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *carteSO(irun,icartleft,jrun,icartright)
      enddo
      enddo
      enddo
      enddo
      else 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *(carteSO(irun,icartleft,jrun,icartright)+
     *2d0*carteOO(irun,icartleft,jrun,icartright))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
      
      subroutine two2mean12b(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
#include <implicit.h>
#include <para.h>
      logical sameorb 
      dimension 
     *carteSO(ncontmf,norbsum,ncontmf,norbsum),
     *carteOO(ncontmf,norbsum,ncontmf,norbsum),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) then 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *carteSO(jrun,icartleft,irun,icartright)
      enddo
      enddo
      enddo
      enddo
      else 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *(carteSO(jrun,icartleft,irun,icartright)+
     *2d0*carteOO(jrun,icartleft,irun,icartright))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
      
      subroutine two2mean13(carteSO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb)
cbs   gives the two first contributions
cbs   < i M | j M >  with Malpha  and Mbeta 
cbs   the other orbit parts cancel    
#include <implicit.h>
#include <para.h>
      dimension carteSO(ncontmf,ncontmf,norbsum,norbsum),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb        
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *carteSO(irun,jrun,icartleft,icartright)
      enddo
      enddo
      enddo
      enddo
c     write(6,*) 'effective integrals' 
c     do jrun=1,ncontmf
c     write(6,'(4E20.14)') (onecart(irun,jrun),irun=1,ncontmf)
c     enddo
      return 
      end 
      
      subroutine two2mean34a(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
#include <implicit.h>
#include <para.h>
      logical sameorb 
      dimension 
     *carteSO(norbsum,ncontmf,norbsum,ncontmf),
     *carteOO(norbsum,ncontmf,norbsum,ncontmf),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) then 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *carteSO(icartleft,irun,icartright,jrun)
      enddo
      enddo
      enddo
      enddo
      else 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *(carteSO(icartleft,irun,icartright,jrun)+
     *2d0*carteOO(icartleft,irun,icartright,jrun))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
      
      subroutine two2mean34b(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
#include <implicit.h>
#include <para.h>
      logical sameorb 
      dimension 
     *carteSO(norbsum,ncontmf,norbsum,ncontmf),
     *carteOO(norbsum,ncontmf,norbsum,ncontmf),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) then 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5D0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *carteSO(icartleft,jrun,icartright,irun)
      enddo
      enddo
      enddo
      enddo
      else 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5D0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *(carteSO(icartleft,jrun,icartright,irun)+
     *2d0*carteOO(icartleft,jrun,icartright,irun))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
