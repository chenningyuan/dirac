// Function table to pick the right left and right L value
void (*overlap_tab[4][4])(matrix<double> &S, 
			  const matrix<double> &c1, 
			  const shell &s1,
			  const matrix<double> &c2,
			  const shell &s2) = 
{
  {overlap<0,0>, overlap<0,1>, overlap<0,2>, overlap<0,3>},
  {overlap<1,0>, overlap<1,1>, overlap<1,2>, overlap<1,3>},
  {overlap<2,0>, overlap<2,1>, overlap<2,2>, overlap<2,3>},
  {overlap<3,0>, overlap<3,1>, overlap<3,2>, overlap<3,3>}
};

void (*kinetic_tab[4][4])(matrix<double> &T,
			  const matrix<double> &c1, 
			  const shell &s1,
			  const matrix<double> &c2,
			  const shell &s2) =
{
  {kinetic<0,0>, kinetic<0,1>, kinetic<0,2>, kinetic<0,3>},
  {kinetic<1,0>, kinetic<1,1>, kinetic<1,2>, kinetic<1,3>},
  {kinetic<2,0>, kinetic<2,1>, kinetic<2,2>, kinetic<2,3>},
  {kinetic<3,0>, kinetic<3,1>, kinetic<3,2>, kinetic<3,3>}
};

void (*coulomb_tab[7][7])(matrix<double> &J,
			  const matrix<double> &c1, 
			  const shell &s1,
			  const matrix<double> &c2,
			  const shell &s2) =
{
  {coulomb<0,0>, coulomb<0,1>, coulomb<0,2>, coulomb<0,3>,
    coulomb<0,4>, coulomb<0,5>, coulomb<0,6>},
  {coulomb<1,0>, coulomb<1,1>, coulomb<1,2>, coulomb<1,3>,
    coulomb<1,4>, coulomb<1,5>, coulomb<1,6>},
  {coulomb<2,0>, coulomb<2,1>, coulomb<2,2>, coulomb<2,3>,
    coulomb<2,4>, coulomb<2,5>, coulomb<2,6>},
  {coulomb<3,0>, coulomb<3,1>, coulomb<3,2>, coulomb<3,3>,
    coulomb<3,4>, coulomb<3,5>, coulomb<3,6>},
  {coulomb<4,0>, coulomb<4,1>, coulomb<4,2>, coulomb<4,3>,
    coulomb<4,4>, coulomb<4,5>, coulomb<4,6>},
  {coulomb<5,0>, coulomb<5,1>, coulomb<5,2>, coulomb<5,3>,
    coulomb<5,4>, coulomb<5,5>, coulomb<5,6>},
  {coulomb<6,0>, coulomb<6,1>, coulomb<6,2>, coulomb<6,3>,
    coulomb<6,4>, coulomb<6,5>, coulomb<6,6>}
};

void (*shell_product_tab[4][4])(matrix<double> &prod, shell &s12,
				const matrix<double> &c1, const shell &s1,
				const matrix<double> &c2, const shell &s2) =
{
  {shell_product<0,0>, shell_product<0,1>, shell_product<0,2>, shell_product<0,3>},
  {shell_product<1,0>, shell_product<1,1>, shell_product<1,2>, shell_product<1,3>},
  {shell_product<2,0>, shell_product<2,1>, shell_product<2,2>, shell_product<2,3>},
  {shell_product<3,0>, shell_product<3,1>, shell_product<3,2>, shell_product<3,3>}
};



void (*pointcharge_tab[7])(matrix<double> &C,
			   const matrix<double> &c1, 
			   const shell &s,
			   const double r[3],
			   double Z) =
{
  pointcharge<0>, pointcharge<1>, pointcharge<2>, pointcharge<3>,
  pointcharge<4>, pointcharge<5>, pointcharge<6>
};
