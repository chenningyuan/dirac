
/*

     subroutine ESRC_LDAERFSPIN (rhot,rhos,mu,dospin,e)
C*****************************************************************************
C     Short-range spin-dependent LDA correlation functional from 
C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
C
C     subroutine adapted from Molpro
C
C     Input: rhot   : total density
C            rhos   : spin density
C            mu     : Interation parameter
C            dospin : use spin density
C
C     Ouput: e      : energy
C
C     Created: 18-08-09, J. Toulouse
C*****************************************************************************
      implicit none

      double precision rhot, rhos, mu
      logical dospin
      double precision e

      double precision f13, pi, rsfac
      double precision rs
      double precision ec,ecd,ecz,eclr

      double precision rhoa, rhob, z
      double precision ea, eb

      parameter(f13 = 0.333333333333333d0)
      parameter(pi = 3.141592653589793d0)
      parameter(rsfac = 0.620350490899400d0)

C     Spin-unpolarized case
      if (dospin) then

        rs = rsfac/(rhot**f13)

        call ecPW(rs,0.0d0,ec,ecd,ecz)
        call ecorrlr(rs,0.0d0,mu,eclr)
        e=(ec-eclr)*rhot

C     Spin-polarized case
      else

        rs=rsfac/(rhot**f13)
        rhoa=max((rhot+rhos)*.5d0,1.0d-15)
        rhob=max((rhot-rhos)*.5d0,1.0d-15)
        z=(rhoa-rhob)/(rhoa+rhob)

        call ecPW(rs,z,ec,ecd,ecz)
        call ecorrlr(rs,z,mu,eclr)
        e=(ec-eclr)*rhot

      endif

      return
      end

*/

template<class num>
T esrc_ldaerfspin(const densvars<num> &d, double mu)
{
/*
C     Short-range spin-dependent LDA correlation functional from 
C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
C
C     subroutine adapted from Molpro
*/
  return d.n*(pw92_eps(d) - ecorrlr(d,mu));
}


template<class num>
T ecorrlr(const densvars<num> &d, double mu)
{
  /*
  double precision rs,z,mu,eclr,ec,ecd,ecz
    double precision alpha,cf,phi
    double precision g0f,dpol,d2anti,d3anti,Qrpa
    double precision coe2,coe3,coe4,coe5
    double precision a1,a2,a3,a4,b0
    double precision q1a,q2a,q3a,t1a,t2a,t3a,adib
  */
  const double alpha=pow(4.0/9.0/M_PI,1.0/3.0);
  double cf=1.d0/alpha; //TODO: The normal CF?
  num phi=(pow(1.0+d.zeta,2.0/3.0)+pow(1.0-d.zeta,2.0/3.0))/2.0;
  // cc parameters from the fit
  const double adib = 0.784949;
  const double q1a = -0.388;
  const double q2a = 0.676;
  const double q3a = 0.547;
  const double t1a = -4.95;
  const double t2a = 1.0;
  const double t3a = 0.31;
  
  num b0 = adib*d.r_s;
  num rs2 = d.r_s*d.r_s;
  num rs3 = rs2*d.r_s;

  num d2anti=(q1a*d.r_s+q2a*rs2)*exp(-q3a*d.r_s)/rs2;
  num d3anti=(t1a*d.r_s+t2a*rs2)*exp(-t3a*d.r_s)/rs3;

  num coe2=-3.0/8.0/rs3*(1.0-pow(d.zeta,2))*(g0f(rs)-0.5);

  num coe3=-(1.0-pow(d.zeta,2))*g0f(d.r_s)/(sqrt(2.0*M_PI)*rs3);

      if(abs(z).eq.1.0) then

        coe4=-9.0/64.0/rs**3*(dpol(d.r_s)
     $        -cf**2*2**(5.0/3.0)/5.0/rs**2)
        coe5=-9.0/40.0/(sqrt(2.0*M_PI)*rs**3)*dpol(d.r_s)

      else

         coe4=-9.0/64.0/rs**3*(((1.0+z)/2.0)**2*
     $        dpol(d.r_s*(2/(1.0+z))**(1.0/3.0))+((1.0-z)/2.0)**2
     $        *dpol(d.r_s*(2.0/(1.0-z))**(1.0/3.0))+
     $        (1.-z**2)*d2anti-cf**2/10.0*((1.0+z)**(8.0/3.0)
     $        +(1.-z)**(8.0/3.0))/rs**2)

         coe5=-9.0/40.0/(sqrt(2.0*M_PI)*rs**3)*(((1.0+z)/2.0)**2
     $        *dpol(d.r_s*(2.0/(1.0+z))**(1.0/3.0))+((1.0-z)/2.0)**2
     $        *dpol(d.r_s*(2.0/(1.0-z))**(1.0/3.0))+(1.0-z**2)*
     $        d3anti)
      end if

      call ecPW(d.r_s,z,ec,ecd,ecz)

      a1=4.0*b0**6*coe3+b0**8*coe5
      a2=4.0*b0**6*coe2+b0**8*coe4+6.0*b0**4*ec
      a3=b0**8*coe3
      a4=b0**6*(b0**2*coe2+4.0*ec)

      eclr=(phi**3*Qrpa(mu*sqrt(d.r_s)/phi)+a1*mu**3+a2*mu**4+a3*mu**5+
     $     a4*mu**6+b0**8*mu**8*ec)/((1.0+b0**2*mu**2)**4)

      return
      end

