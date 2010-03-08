// this file is ment to be included from xc_fun.cpp 

// KT exchange GGA correction

// reference:
// @article{keal:3015,
// author = {Thomas W. Keal and David J. Tozer},
// collaboration = {},
// title = {The exchange-correlation potential in Kohn--Sham nuclear magnetic resonance shielding calculations},
// publisher = {AIP},
// year = {2003},
// journal = {The Journal of Chemical Physics},
// volume = {119},
// number = {6},
// pages = {3015-3024},
// keywords = {eigenvalues and eigenfunctions; ab initio calculations; bond lengths; nuclear screening; nuclear magnetic resonance; density functional theory; ionisation potential; dissociation energies},
// url = {http://link.aip.org/link/?JCP/119/3015/1},
// doi = {10.1063/1.1590634}
// }

// autodiff version: Radovan Bast (radovan.bast@uit.no)
// tested against implementation in Dalton by Dave Wilson (davidwi@kjemi.uio.no)
// compared first derivatives only

template<class num>
static num kt_exchange(const densvars<num> &d)
{
  const double KT_THRESHOLD = 1e-14;
  const double DELTA = 0.1;  

  num ea = 0.0;
  if (d.a > KT_THRESHOLD) 
  {
    ea = d.gaa/(DELTA + pow(d.a, 4.0/3.0));
  }

  num eb = 0.0;
  if (d.b > KT_THRESHOLD) 
  {
    eb = d.gbb/(DELTA + pow(d.b, 4.0/3.0));
  }

  return (ea + eb);
}
