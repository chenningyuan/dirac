#include "pw92_fun.h"
#include "pw9xx_fun.h"


namespace PBEx_internal
{
   using pw91_like_x_internal::prefactor;

   // enhancement factor F(S), common to PBEx and REVPBEx
  template<class num>
  static num pbex_enhancement(const double &R, 
			      const num &rho,
			      const num &grad)
   {
      using pw91_like_x_internal::S;

      double mu = 0.066725*M_PI*M_PI/3.0;
      // ulfek: mu from Daresbury implementation
    //double mu = 0.2195149727645171;
    //radovan: commented out above line because
    //  should be 0.21951645...
      num st = S(rho,grad);
      num t1 = 1 + mu*st*st/R;
   
      // TODO: veryfy 1/t1 does not overflow/underflow

      return 1 + R - R/t1; 
   }

  // enhancement factor F(S) for RPBE, see Hammer, Hansen and Norskov, PRB 1999, 59, 7413
  template<class num>
  static num rpbex_enhancement(const double &R, 
			       const num &rho,
			       const num &grad)
   {
      using pw91_like_x_internal::S;

      // double mu = 0.066725*M_PI*M_PI/3.0;
      // ulfek: mu from Daresbury implementation
      double mu = 0.2195149727645171;
      num st = S(rho,grad);

      return R*(1 - exp(-mu*st*st/R)); 
   }


   // original PBE exchange functional

   template<class num>
   static num pbe_exchange(const densvars<num> &d)
   {
     const double R = 0.804;

     return (  prefactor(d.a,d.gaa)*pbex_enhancement(R,d.a,d.gaa) 
             + prefactor(d.b,d.gbb)*pbex_enhancement(R,d.b,d.gbb) );
   }

   // REVPBE exchange functional

   template<class num>
   static num revpbe_exchange(const densvars<num> &d)
   {
     const double R = 1.245;

     return (  prefactor(d.a,d.gaa)*pbex_enhancement(R,d.a,d.gaa) 
             + prefactor(d.b,d.gbb)*pbex_enhancement(R,d.b,d.gbb) );
   }

   // RPBE exchange functional

   template<class num>
   static num rpbe_exchange(const densvars<num> &d)
   {
     const double R = 0.804;

     return (  prefactor(d.a,d.gaa)*rpbex_enhancement(R,d.a,d.gaa) 
             + prefactor(d.b,d.gbb)*rpbex_enhancement(R,d.b,d.gbb) );
   }

}


// PBEc implementation from MOLPRO manual, by ulfek.
namespace PBEc_internal
{
  //ulfek: I cannot find these constants directly in the PBE
  //paper
  //  const double iota = 0.0716;
  //MOLPRO value: const double kappa = 0.004235;
  //const double kappa = 0.00423488752945734;
  //const double nu = 16*pow(3*M_PI*M_PI,1.0/3)/M_PI;
  //const double lambda = nu*kappa;
  const double param_gamma = (1-log(2.0))/(M_PI*M_PI);
  //ulfek: beta comes from the slowly varying limit, see PBE paper
  const double param_beta_pbe_paper = 0.066725;
  const double param_beta_accurate  = 0.06672455060314922;
  //  const double param_omega = 0.046644;
  const double param_beta_gamma = param_beta_accurate/param_gamma;
  

  template<class num>
  static num u(const num &zeta)
  {
    // TODO: check if |zeta| ~= 1 
    return 0.5*(pow(1+zeta,2.0/3)+(pow(1-zeta,2.0/3)));
  }
  
  template<class num>
  static num A(const num &eps, const num &zeta)
  {
    return param_beta_gamma/(exp(-eps/(param_gamma*pow(u(zeta),3))) - 1);
  }

  template<class num>
  static num H(const num &d, const num &eps, const num &zeta)
  {
  num d2 = d*d;
  num d2A = d2*A(eps,zeta);
  return param_gamma*pow(u(zeta),3)*
    log(1+param_beta_gamma*d2*(1 + d2A)/(1+d2A*(1+d2A)));
  }

  template<class num>
  static num pbe_correlation(const densvars<num> &d)
  {
    const num &sigma = d.gnn;
    num dd = 1.0/12*(sqrt(sigma)*pow(3,5.0/6))/
      (u(d.zeta)*pow(M_PI,-1.0/6)*pow(d.n,7.0/6));
    num eps = pw92_eps_pbe(d);
    return d.n*(eps + H(dd,eps,d.zeta));
  }
}

using PBEc_internal::pbe_correlation;
using PBEx_internal::pbe_exchange;
using PBEx_internal::revpbe_exchange;
using PBEx_internal::rpbe_exchange;


