#include "m0xy_fun.h"


namespace M05_c_internal
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::m05_c_anti;
   using m0xy_metagga_xc_internal::ueg_c_anti;
   using m0xy_metagga_xc_internal::m05_c_para;
   using m0xy_metagga_xc_internal::ueg_c_para;

   // M05 correlation functional. 

   template<class num>
   static num m05_correlation(const densvars<num> &d)
   {
     // parameters for anti-parallel spin contributions
     const static double param_c_anti[5] =
       {  1.000000e+00,  3.785690e+00, -1.415261e+01, -7.465890e+00,  1.794491e+00 };

     // parameters for parallel spin contributions
     const static double param_c_para[5] =
       {  1.000000e+00,  3.773440e+00, -2.604463e+01,  3.069913e+01, -9.226950e+00 };

     num chi_a = chi(d.a, d.gaa);
     num chi_b = chi(d.b, d.gbb);
     num zet_a = zet(d.a, d.taua);
     num zet_b = zet(d.b, d.taub);

     num Ec_ab = ueg_c_anti(d)  *m05_c_anti(param_c_anti[5],chi_a,chi_b); 
     num Ec_aa = ueg_c_para(d.a)*m05_c_para(param_c_para[5],chi_a,zet_a);
     num Ec_bb = ueg_c_para(d.b)*m05_c_para(param_c_para[5],chi_b,zet_b);

     return Ec_ab + Ec_aa + Ec_bb;
   }

   // M05-2X correlation functional. 

   template<class num>
   static num m05_2x_correlation(const densvars<num> &d)
   {
     // parameters for anti-parallel spin contributions
     const static double param_c_anti[5] =
       {  1.000000e+00,  1.092970e+00, -3.791710e+00,  2.828100e+00, -1.058909e+01 };

     // parameters for parallel spin contributions
     const static double param_c_para[5] =
       {  1.000000e+00, -3.054300e+00,  7.618540e+00,  1.476650e+00, -1.192365e+01 };

     num chi_a = chi(d.a, d.gaa);
     num chi_b = chi(d.b, d.gbb);
     num zet_a = zet(d.a, d.taua);
     num zet_b = zet(d.b, d.taub);

     num Ec_ab = ueg_c_anti(d)  *m05_c_anti(param_c_anti[5],chi_a,chi_b); 
     num Ec_aa = ueg_c_para(d.a)*m05_c_para(param_c_para[5],chi_a,zet_a);
     num Ec_bb = ueg_c_para(d.b)*m05_c_para(param_c_para[5],chi_b,zet_b);

     return Ec_ab + Ec_aa + Ec_bb;
   }
}

using M05_c_internal::m05_correlation;
using M05_c_internal::m05_2x_correlation;

#include "m0xy_fun.h"

namespace M06_x_internal
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::fw;
   using m0xy_metagga_xc_internal::h;
   using PBEx_internal::pbe_exchange;
   using m0xy_metagga_xc_internal::alpha_x;

   // M06-L exchange functional. Does not use HF exchange. 

   template<class num>
   static num m06_l_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  3.987756e-01,  2.548219e-01,  2.923994e-01, -2.103655e+00, -6.302147e+00, 
          1.097615e+01,  3.097273e+01, -2.318489e+01, -5.673480e+01,  2.160364e+01,
          3.421814e+01, -9.049762e+00 }; 
     const static double param_d[5] = 
       {  6.012244e-01,  4.748822e-03, -8.635108e-04, -9.308062e-06,  4.482811e-05 }; 

     num chia2 = pow(chi(d.a, d.gaa), 2);
     num chib2 = pow(chi(d.b, d.gbb), 2);

     return (  (  pbe_exchange(d.a, d.gaa)*fw(param_a, d.a, d.taua) 
                + lsda_x(d.a)*h(param_d, alpha_x, chia2, zet(d.a, d.taua)))
              +
               (  pbe_exchange(d.b, d.gbb)*fw(param_a, d.b, d.taub) 
                + lsda_x(d.b)*h(param_d, alpha_x, chib2, zet(d.b, d.taub)))
            );
   } 

   // M06-HF exchange functional. Uses full HF exchange. 

   template<class num>
   static num m06_hf_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  1.179732e-01, -1.066708e+00, -1.462405e-01,  7.481848e+00,  3.776679e+00, 
         -4.436118e+01, -1.830962e+01,  1.003903e+02,  3.864360e+01, -9.806018e+01,
          2.557716e+01,  3.590404e+01 }; 
     const static double param_d[5] = 
       { -1.179732e-01, -2.500000e-03, -1.180056e-02,  0.000000e+00,  0.000000e+00 }; 

     num chia2 = pow(chi(d.a, d.gaa), 2);
     num chib2 = pow(chi(d.b, d.gbb), 2);

     return (  (  pbe_exchange(d.a, d.gaa)*fw(param_a, d.a, d.taua) 
                + lsda_x(d.a)*h(param_d, alpha_x, chia2, zet(d.a, d.taua)))
              +
               (  pbe_exchange(d.b, d.gbb)*fw(param_a, d.b, d.taub) 
                + lsda_x(d.b)*h(param_d, alpha_x, chib2, zet(d.b, d.taub)))
            );
   } 

}

using M06_x_internal::m06_l_exchange;
using M06_x_internal::m06_hf_exchange;


namespace M06_c_internal
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::m06_c_anti;
   using m0xy_metagga_xc_internal::ueg_c_anti;
   using m0xy_metagga_xc_internal::m06_c_para;
   using m0xy_metagga_xc_internal::ueg_c_para;

   // M06 correlation functional. 

   template<class num>
   static num m06_correlation(const densvars<num> &d)
   {
     // parameters for anti-parallel spin contributions
     const static double param_c_anti[5] =
       {  3.741539e+00,  2.187098e+02, -4.531252e+02,  2.936479e+02, -6.287470e+01 };
     const static double param_d_anti[5] =
       { -2.741539e+00, -6.720113e-01, -7.932688e-02,  1.918681e-03, -2.032902e-03 };

     // parameters for parallel spin contributions
     const static double param_c_para[5] =
       {  5.094055e-01, -1.491085e+00,  1.723922e+01, -3.859018e+01,  2.845044e+01 };
     const static double param_d_para[5] =
       {  4.905945e-01, -1.437348e-01,  2.357824e-01,  1.871015e-03, -3.788963e-03 };

     num chi_a = chi(d.a, d.gaa);
     num chi_b = chi(d.b, d.gbb);
     num zet_a = zet(d.a, d.taua);
     num zet_b = zet(d.b, d.taub);

     num Ec_ab = ueg_c_anti(d)*m06_c_anti(param_c_anti[5],param_d_anti[5],chi_a,chi_b,zet_a,zet_b); 
     num Ec_aa = ueg_c_para(d.a)*m06_c_para(param_c_para[5],param_d_para[5],chi_a,zet_a);
     num Ec_bb = ueg_c_para(d.b)*m06_c_para(param_c_para[5],param_d_para[5],chi_b,zet_b);

     return Ec_ab + Ec_aa + Ec_bb;
   }
}

using M06_c_internal::m06_correlation;

#ifndef M0xy_H
#define M0xy_H

// common functions for MO5 and M06 family of (hybrid) meta-gga functionals

namespace m0xy_metagga_xc_internal
{

  const char *m05_2x_x_reference = 
    "M05-2X Meta-GGA Exchange Functional\n"
      "Y Zhao, N. E. Schultz and D. G. Truhlar, "
      "J. Chem. Theory Comput. 2, 364 (2006)\n"; 
 
  const char *m05_x_reference = 
    "M05 Meta-GGA Exchange Functional\n"
      "Y Zhao, N. E. Schultz and D. G. Truhlar, "
      "J. Chem. Theory Comput. 2, 364 (2006)\n"; 
 
  const char *m06_l_x_reference = 
    "M06-L Meta-GGA Exchange Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "J. Chem. Phys. 125, 194101 (2006)\n"; 

  const char *m06_2x_x_reference = 
    "M06-2X Meta-GGA Exchange Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "Theor. Chem. Account 120, 215 (2008)\n"; 

  const char *m06_x_reference = 
    "M06 Meta-GGA Exchange Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "Theor. Chem. Account 120, 215 (2008)\n"; 


  const char *m05_2x_c_reference = 
    "M05-2X Meta-GGA Correlation Functional\n"
      "Y Zhao, N. E. Schultz and D. G. Truhlar, "
      "J. Chem. Theory Comput. 2, 364 (2006)\n"; 
 
  const char *m05_c_reference = 
    "M05 Meta-GGA Correlation Functional\n"
      "Y Zhao, N. E. Schultz and D. G. Truhlar, "
      "J. Chem. Theory Comput. 2, 364 (2006)\n"; 
 
  const char *m06_l_c_reference = 
    "M06-L Meta-GGA Correlation Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "J. Chem. Phys. 125, 194101 (2006)\n"; 

  const char *m06_2x_c_reference = 
    "M06-2X Meta-GGA Correlation Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "Theor. Chem. Account 120, 215 (2008)\n"; 

  const char *m06_c_reference = 
    "M06 Meta-GGA Correlation Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "Theor. Chem. Account 120, 215 (2008)\n"; 


  // additional information from the authors, and reference implementations 
  // in F77, can be found at
  //
  //    http://comp.chem.umn.edu/mfm/

  // values for alpha are taken from
  //   T. Van Voorhis and G. E. Scuseria, J. Chem. Phys. 129, 219901 (2008)
  //   DOI: 10.1063/1.3005348

  const static double alpha_x = 0.00186726;
  const static double alpha_c_parallel = 0.00515088;
  const static double alpha_c_antiparallel = 0.00304966;


  // zet calculates one of the working variables
  // see see TCA (2008) 120:215 eq 3
  // rho is the density, tau the kinetic energy density

  template<class num>
  static num zet(const num &rho, const num &tau)
  {
    using xc_constants::CF;

    return ((2.0*tau/pow(rho,5.0/3.0)) - CF);
  }

  // gamma calculates the quantity of TCA (2008) 120:215 eq 4, but we supply it with chi^2 as parameter

  template<class num>
  static num gamma(const static double alpha, const num &chi2, const num &zet)
  {
    return 1 + alpha*(chi2 + zet); 
  }

  // h calculates the quantity of TCA (2008) 120:215 eq 5, but here we provide it with the square of the x_\sigma
  // variable, since that is used always in terms of x^2

  template<class num>
  static num h(const static double d[], const static double alpha, const num &chi2, const num &zet)
  {
    num gam1 = gamma(alpha,chi2,zet);

    num t1 =  d[0]/gam1;
    num t2 = (d[1]*chi2 + d[2]*zet)/pow(gam1, 2);
    num t3 = (d[3]*pow(chi2,2) + d[4]*chi2*zet + d[5]*pow(zet, 2))/pow(gam1, 3);

    return t1+t2+t3;
  }


  // functions specific to the exchange functinals

  // fw is the spin kinetic enrgy density enhancement factor, see e.g.
  // TCA 2008 120:215 eq 8
  // here the number of parameters is hardcoded since all M05/M06 functionals
  // described at the time of coding (feb 2010) use 12 of them.

  template<class num>
  static num fw(const static double a[12], const num &rho, const num &tau)
  {
    using pw91_like_x_internal::pw91k_prefactor;

    num tau_lsda = pw91k_prefactor(rho);

    num  t = tau_lsda/tau;
    num  w = (t - 1)/(t + 1); 
    num fw = 0.0;

    for (int i = 0; i < 12; i++)
       fw += a[i] * pow(w,i);

    return fw;
  }


  // functions specific to the correlation functionals

  // Dsigma is the self-interaction correction factor, for parallel spins, used in M06 and M06-2X 
  // see for instance TCA 2008 120:215 eq 16
  // obviously, chi and zet are to have same spin: Dsigma(x_\alpha, z_\alpha) or Dsigma(x_\beta, z_\beta)
  // chi2 is chi*chi

  template<class num>
  static num Dsigma(const num &chi2, const num &zet)
  {
    using xc_constants::CF;

    return (1.0 - 0.25*chi2/(zet + CF)); 
  }

  // g is an auxiliary function that can be used in connection to both
  // parallel and antiparallel spins (provided in the latter that the variable
  // passed on is not chi_\sigma^2 but (chi_a^2 + chi_b^2)). 
  //
  // for compactness, we also pass (chi_\sigma^2 or (chi_a^2 + chi_b^2)) multiplied
  // by the corresponding gamma factor, since the two are always multiplied together.
  //
  // see TCA 2008 120:215 eqs 13 and 15 for the actual formulas
  //
  // here the number of parameters is hardcoded since all M05/M06 functionals
  // described at the time of coding (feb 2010) use 5 of them.

  template<class num>
  static num g(const static double param_c[5], const num &gamma_chi_squared)
  {
     num g = 0.0;
     num b   = gamma_chi_squared / (1.0 + gamma_chi_squared); 
  
     for (int i = 0; i < 5; i++)
       g += param_c[i] * pow(b,i);

    return g;
  }


  // m06_c_anti calculates the  (g + h) part of the integrand, for the M06y correlationg
  // functionals, for the antiparallel spins
  //
  // m06_c_para calculates the  (g + h)D part of the integrand, for the M06y correlationg
  // functionals, for parallel spins

  template<class num>
  static num m06_c_anti(const static double param_c[5], 
                        const static double param_d[5], 
                        const num &chi_a, const num &zet_a, 
                        const num &chi_b, const num &zet_b) 
  {
    const num gamma_c_anti = 0.0031;  // this is an "universal" constant for all M05/M06 functionals 

    const num zet_ab  = zet_a + zet_b;
    const num chi_ab2 = chi_a*chi_a + chi_b*chi_b;

    return g(param_c,gamma_c_anti*chi_ab2) + h(param_d,alpha_c_antiparallel,chi_ab2,zet_ab);
  }

  template<class num>
  static num m06_c_para(const static double param_c[5], 
                        const static double param_d[5], 
                        const num &chi, const num &zet) 
  {
    const num gamma_c_parallel = 0.06;  // this is an "universal" constant for all M05/M06 functionals 

    const num chi2 = chi*chi;

    return (g(param_c,gamma_c_parallel*chi2) + h(param_d,alpha_c_parallel,chi,zet))*Dsigma(chi2,zet);
  }


  // m05_c_anti calculates the  g part of the integrand, for the M05y correlationg
  // functionals, for the antiparallel spins
  //
  // m05_c_para calculates the  gD part of the integrand, for the M05y correlationg
  // functionals, for parallel spins

  template<class num>
  static num m05_c_anti(const static double param_c[5], const num &chi_a, const num &chi_b) 
  {
    const num gamma_c_anti = 0.0031;  // this is an "universal" constant for all M05/M06 functionals 

    const num chi_ab2 = chi_a*chi_a + chi_b*chi_b;

    return g(param_c,gamma_c_anti*chi_ab2);
  }

  template<class num>
  static num m05_c_para(const static double param_c[5], const num &chi, const num &zet) 
  {
    const num gamma_c_parallel = 0.06;  // this is an "universal" constant for all M05/M06 functionals 

    const num chi2 = chi*chi;

    return g(param_c,gamma_c_parallel*chi2)*Dsigma(chi2,zet);
  }


  // functions to calculate the \epsilon^{UEC}^{[\alpha\beta,\sigma\sigma]} factor in the
  // M05/M06 (and VSXC) functionals. The definitions of those can be found in 
  // 
  // Zhao, Schultz, Truhlar, J. Chem. Theory Comput. 2006, 2, 264   equations 7,8
  //  

  using pw92_internal::pw92_eps;
  using pw92_internal::pw92_eps_polarized;

  template<class num>
  static num ueg_c_anti(const densvars<num> &d)
  {
    return pw92_eps(d) - pw92_eps_polarized(d.a) - pw92_eps_polarized(d.b);
  }

  template<class num>
  static num ueg_c_para(const num &rho)
  {
     return pw92_eps_polarized(rho);
  }


}

// local spin density approximation for the exchange
template<class num>
static num lsda_x(const num &rho)
{
  return -(3.0/2.0)*pow(3.0/(4.0*M_PI), 1.0/3.0)*pow(rho, 4.0/3.0);
}


#endif
