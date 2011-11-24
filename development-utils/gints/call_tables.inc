#define ARRAY88(n) \
{ {n<0,0>, n<0,1>, n<0,2>, n<0,3>, n<0,4>, n<0,5>, n<0,6>, n<0,7>},\
  {n<1,0>, n<1,1>, n<1,2>, n<1,3>, n<1,4>, n<1,5>, n<1,6>, n<1,7>},\
  {n<2,0>, n<2,1>, n<2,2>, n<2,3>, n<2,4>, n<2,5>, n<2,6>, n<2,7>},\
  {n<3,0>, n<3,1>, n<3,2>, n<3,3>, n<3,4>, n<3,5>, n<3,6>, n<3,7>},\
  {n<4,0>, n<4,1>, n<4,2>, n<4,3>, n<4,4>, n<4,5>, n<4,6>, n<4,7>},\
  {n<5,0>, n<5,1>, n<5,2>, n<5,3>, n<5,4>, n<5,5>, n<5,6>, n<5,7>},\
  {n<6,0>, n<6,1>, n<6,2>, n<6,3>, n<6,4>, n<6,5>, n<6,6>, n<6,7>},\
  {n<7,0>, n<7,1>, n<7,2>, n<7,3>, n<7,4>, n<7,5>, n<7,6>, n<7,7>}}

#define ARRAY44(n) \
{ {n<0,0>, n<0,1>, n<0,2>, n<0,3>, },\
  {n<1,0>, n<1,1>, n<1,2>, n<1,3>, },\
  {n<2,0>, n<2,1>, n<2,2>, n<2,3>, },\
  {n<3,0>, n<3,1>, n<3,2>, n<3,3>, }}	

// Function table to pick the right left and right L value
void (*overlap_tab[8][8])(matrix<double> &S, 
			  const matrix<double> &c1, 
			  const shell &s1,
			  const matrix<double> &c2,
			  const shell &s2) = 
  ARRAY88(overlap);

void (*kinetic_tab[8][8])(matrix<double> &T,
			  const matrix<double> &c1, 
			  const shell &s1,
			  const matrix<double> &c2,
			  const shell &s2) =
  ARRAY88(kinetic);


void (*coulomb_tab[8][8])(matrix<double> &J,
			  const matrix<double> &c1, 
			  const shell &s1,
			  const matrix<double> &c2,
			  const shell &s2) =
  ARRAY44(coulomb);

void (*shell_product_tab[8][8])(matrix<double> &prod, shell &s12,
				const matrix<double> &c1, const shell &s1,
				const matrix<double> &c2, const shell &s2) =
  ARRAY44(shell_product);


void (*pointcharge_tab[8])(matrix<double> &C,
			   const matrix<double> &c1, 
			   const shell &s,
			   const double r[3],
			   double Z) =
{
  pointcharge<0>, pointcharge<1>, pointcharge<2>, pointcharge<3>,
  pointcharge<4>, pointcharge<5>, pointcharge<6>, pointcharge<7>,
};
