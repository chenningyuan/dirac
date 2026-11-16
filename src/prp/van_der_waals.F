!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

      subroutine van_der_waals(nfreq,alpha,work,lfree)
C***********************************************************************
C
C     Written by R. Pilemalm and P. Norman 2008
C
C***********************************************************************
      use priunit_mod, only: LUPRI
#include "implicit.h"
C   
      real*8 C6
      dimension alpha(nfreq),work(lfree)
C
      kfree=1
      lwrk=lfree
C
      c6=0.0d0
      call gauss_legendre(nfreq,alpha,alpha,c6)
C
       CALL AROUND('Dispersion coefficient')
       WRITE(LUPRI,'(A)')
     &  'Please cite:',
     &  '   Sebastien Villaume, Trond Saue and Patrick Norman,',
     &  '   J. Chem. Phys. 133 (2010) 064105'
       WRITE(LUPRI,'(A)')
       WRITE(LUPRI,'(A)')'@'
       WRITE(LUPRI,'(A3,A10,F15.8,A5)') '@  ',
     &                         'C6       ', C6, 'a.u.'
       WRITE(LUPRI,'(A)')'@'
       WRITE(LUPRI,'(/A)')
     &  '----------------------------------------------'
       
C
      return 
      end
      subroutine gauss_legendre(N,A,B,C6)
C***********************************************************************
C
C     Written by R. Pilemalm and P. Norman 2008
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use pi_mod, only: PI
      use gauss_legendre_mod, only: ABSCISSAS_10PT, ABSCISSAS_12PT,
     &  ABSCISSAS_14PT, ABSCISSAS_16PT, ABSCISSAS_6PT, ABSCISSAS_8PT,
     &  MAX_GL_PT, WEIGHTS_10PT, WEIGHTS_12PT, WEIGHTS_14PT,
     &  WEIGHTS_16PT, WEIGHTS_6PT, WEIGHTS_8PT
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
C
C     Transformation of integration interval is made with w0=0.3 
C     in accordance with
C       R. D. Amos, N. C. Handy, P. J. Knowles, J. E. Rice and A. J. Stone
C       J. Phys. Chem. 89 (1985) 2186-2192
C     There is a slight typo in Eq. 9. The correct equation is Eq. 4 of
C       A. Jiemchooroj, P. Norman and B. E. Sernelius,
C       J. Chem. Phys. 125 (2006) 124306
C
      REAL(dp), PARAMETER :: w0 =0.3D0
C
      dimension A(N),B(N),abscissas(max_gl_pt),weights(max_gl_pt)
C
      if (N.eq.6) then
         do i=1,n
            abscissas(i)=abscissas_6pt(i)
            weights(i)=weights_6pt(i)
         end do
      else if (N.eq.8) then
         do i=1,n
            abscissas(i)=abscissas_8pt(i)
            weights(i)=weights_8pt(i)
         end do
      else if (N.eq.10) then
         do i=1,n
            abscissas(i)=abscissas_10pt(i)
            weights(i)=weights_10pt(i)
         end do
      else if (N.eq.12) then
         do i=1,n
            abscissas(i)=abscissas_12pt(i)
            weights(i)=weights_12pt(i)
         end do
      else if (N.eq.14) then
         do i=1,n
            abscissas(i)=abscissas_14pt(i)
            weights(i)=weights_14pt(i)
         end do
      else if (N.eq.16) then
         do i=1,n
            abscissas(i)=abscissas_16pt(i)
            weights(i)=weights_16pt(i)
         end do
      else
         write(lupri,*) 'Error in gauss-legendre integration'
         call quit('Error in gauss-legendre integration')
      endif
C     
      DO I=1,N
         C6 = C6 + weights(I)*A(I)*B(I)/(1 - abscissas(I))**2
      END DO
C
      C6 = 6.0D0/pi * w0 * C6
C     
      return
      end
