// Function table to pick the right left and right L value
void (*overlaps[4][4])(const double R1[], double alpha1,
		       const double R2[], double alpha2) = 
{
  {print_overlap<0,0>, print_overlap<0,1>, print_overlap<0,2>, print_overlap<0,3>},
  {print_overlap<1,1>, print_overlap<1,1>, print_overlap<1,2>, print_overlap<1,3>},
  {print_overlap<2,2>, print_overlap<2,1>, print_overlap<2,2>, print_overlap<2,3>},
  {print_overlap<3,3>, print_overlap<3,1>, print_overlap<3,2>, print_overlap<3,3>}
};

void (*kinetics[4][4])(const double R1[], double alpha1,
		       const double R2[], double alpha2) = 
{
  {print_kinetic<0,0>, print_kinetic<0,1>, print_kinetic<0,2>, print_kinetic<0,3>},
  {print_kinetic<1,1>, print_kinetic<1,1>, print_kinetic<1,2>, print_kinetic<1,3>},
  {print_kinetic<2,2>, print_kinetic<2,1>, print_kinetic<2,2>, print_kinetic<2,3>},
  {print_kinetic<3,3>, print_kinetic<3,1>, print_kinetic<3,2>, print_kinetic<3,3>}
};

void (*coulombs[4][4][4][4])(const double R[], const double alpha[]) = 
{
  {
    {
      {print_coulomb<0,0,0,0>, print_coulomb<0,0,0,1>, print_coulomb<0,0,0,2>, print_coulomb<0,0,0,3>},
      {print_coulomb<0,0,1,0>, print_coulomb<0,0,1,1>, print_coulomb<0,0,1,2>, print_coulomb<0,0,1,3>},
      {print_coulomb<0,0,2,0>, print_coulomb<0,0,2,1>, print_coulomb<0,0,2,2>, print_coulomb<0,0,2,3>},
      {print_coulomb<0,0,3,0>, print_coulomb<0,0,3,1>, print_coulomb<0,0,3,2>, print_coulomb<0,0,3,3>},
    },
    {
      {print_coulomb<0,1,0,0>, print_coulomb<0,1,0,1>, print_coulomb<0,1,0,2>, print_coulomb<0,1,0,3>},
      {print_coulomb<0,1,1,0>, print_coulomb<0,1,1,1>, print_coulomb<0,1,1,2>, print_coulomb<0,1,1,3>},
      {print_coulomb<0,1,2,0>, print_coulomb<0,1,2,1>, print_coulomb<0,1,2,2>, print_coulomb<0,1,2,3>},
      {print_coulomb<0,1,3,0>, print_coulomb<0,1,3,1>, print_coulomb<0,1,3,2>, print_coulomb<0,1,3,3>},
    },
    {
      {print_coulomb<0,2,0,0>, print_coulomb<0,2,0,1>, print_coulomb<0,2,0,2>, print_coulomb<0,2,0,3>},
      {print_coulomb<0,2,1,0>, print_coulomb<0,2,1,1>, print_coulomb<0,2,1,2>, print_coulomb<0,2,1,3>},
      {print_coulomb<0,2,2,0>, print_coulomb<0,2,2,1>, print_coulomb<0,2,2,2>, print_coulomb<0,2,2,3>},
      {print_coulomb<0,2,3,0>, print_coulomb<0,2,3,1>, print_coulomb<0,2,3,2>, print_coulomb<0,2,3,3>},
    },
    {
      {print_coulomb<0,3,0,0>, print_coulomb<0,3,0,1>, print_coulomb<0,3,0,2>, print_coulomb<0,3,0,3>},
      {print_coulomb<0,3,1,0>, print_coulomb<0,3,1,1>, print_coulomb<0,3,1,2>, print_coulomb<0,3,1,3>},
      {print_coulomb<0,3,2,0>, print_coulomb<0,3,2,1>, print_coulomb<0,3,2,2>, print_coulomb<0,3,2,3>},
      {print_coulomb<0,3,3,0>, print_coulomb<0,3,3,1>, print_coulomb<0,3,3,2>, print_coulomb<0,3,3,3>},
    }},
  {
    {
      {print_coulomb<1,0,0,0>, print_coulomb<1,0,0,1>, print_coulomb<1,0,0,2>, print_coulomb<1,0,0,3>},
      {print_coulomb<1,0,1,0>, print_coulomb<1,0,1,1>, print_coulomb<1,0,1,2>, print_coulomb<1,0,1,3>},
      {print_coulomb<1,0,2,0>, print_coulomb<1,0,2,1>, print_coulomb<1,0,2,2>, print_coulomb<1,0,2,3>},
      {print_coulomb<1,0,3,0>, print_coulomb<1,0,3,1>, print_coulomb<1,0,3,2>, print_coulomb<1,0,3,3>},
    },
    {
      {print_coulomb<1,1,0,0>, print_coulomb<1,1,0,1>, print_coulomb<1,1,0,2>, print_coulomb<1,1,0,3>},
      {print_coulomb<1,1,1,0>, print_coulomb<1,1,1,1>, print_coulomb<1,1,1,2>, print_coulomb<1,1,1,3>},
      {print_coulomb<1,1,2,0>, print_coulomb<1,1,2,1>, print_coulomb<1,1,2,2>, print_coulomb<1,1,2,3>},
      {print_coulomb<1,1,3,0>, print_coulomb<1,1,3,1>, print_coulomb<1,1,3,2>, print_coulomb<1,1,3,3>},
    },
    {
      {print_coulomb<1,2,0,0>, print_coulomb<1,2,0,1>, print_coulomb<1,2,0,2>, print_coulomb<1,2,0,3>},
      {print_coulomb<1,2,1,0>, print_coulomb<1,2,1,1>, print_coulomb<1,2,1,2>, print_coulomb<1,2,1,3>},
      {print_coulomb<1,2,2,0>, print_coulomb<1,2,2,1>, print_coulomb<1,2,2,2>, print_coulomb<1,2,2,3>},
      {print_coulomb<1,2,3,0>, print_coulomb<1,2,3,1>, print_coulomb<1,2,3,2>, print_coulomb<1,2,3,3>},
    },
    {
      {print_coulomb<1,3,0,0>, print_coulomb<1,3,0,1>, print_coulomb<1,3,0,2>, print_coulomb<1,3,0,3>},
      {print_coulomb<1,3,1,0>, print_coulomb<1,3,1,1>, print_coulomb<1,3,1,2>, print_coulomb<1,3,1,3>},
      {print_coulomb<1,3,2,0>, print_coulomb<1,3,2,1>, print_coulomb<1,3,2,2>, print_coulomb<1,3,2,3>},
      {print_coulomb<1,3,3,0>, print_coulomb<1,3,3,1>, print_coulomb<1,3,3,2>, print_coulomb<1,3,3,3>},
    }},
  {
    {
      {print_coulomb<2,0,0,0>, print_coulomb<2,0,0,1>, print_coulomb<2,0,0,2>, print_coulomb<2,0,0,3>},
      {print_coulomb<2,0,1,0>, print_coulomb<2,0,1,1>, print_coulomb<2,0,1,2>, print_coulomb<2,0,1,3>},
      {print_coulomb<2,0,2,0>, print_coulomb<2,0,2,1>, print_coulomb<2,0,2,2>, print_coulomb<2,0,2,3>},
      {print_coulomb<2,0,3,0>, print_coulomb<2,0,3,1>, print_coulomb<2,0,3,2>, print_coulomb<2,0,3,3>},
    },
    {
      {print_coulomb<2,1,0,0>, print_coulomb<2,1,0,1>, print_coulomb<2,1,0,2>, print_coulomb<2,1,0,3>},
      {print_coulomb<2,1,1,0>, print_coulomb<2,1,1,1>, print_coulomb<2,1,1,2>, print_coulomb<2,1,1,3>},
      {print_coulomb<2,1,2,0>, print_coulomb<2,1,2,1>, print_coulomb<2,1,2,2>, print_coulomb<2,1,2,3>},
      {print_coulomb<2,1,3,0>, print_coulomb<2,1,3,1>, print_coulomb<2,1,3,2>, print_coulomb<2,1,3,3>},
    },
    {
      {print_coulomb<2,2,0,0>, print_coulomb<2,2,0,1>, print_coulomb<2,2,0,2>, print_coulomb<2,2,0,3>},
      {print_coulomb<2,2,1,0>, print_coulomb<2,2,1,1>, print_coulomb<2,2,1,2>, print_coulomb<2,2,1,3>},
      {print_coulomb<2,2,2,0>, print_coulomb<2,2,2,1>, print_coulomb<2,2,2,2>, print_coulomb<2,2,2,3>},
      {print_coulomb<2,2,3,0>, print_coulomb<2,2,3,1>, print_coulomb<2,2,3,2>, print_coulomb<2,2,3,3>},
    },
    {
      {print_coulomb<2,3,0,0>, print_coulomb<2,3,0,1>, print_coulomb<2,3,0,2>, print_coulomb<2,3,0,3>},
      {print_coulomb<2,3,1,0>, print_coulomb<2,3,1,1>, print_coulomb<2,3,1,2>, print_coulomb<2,3,1,3>},
      {print_coulomb<2,3,2,0>, print_coulomb<2,3,2,1>, print_coulomb<2,3,2,2>, print_coulomb<2,3,2,3>},
      {print_coulomb<2,3,3,0>, print_coulomb<2,3,3,1>, print_coulomb<2,3,3,2>, print_coulomb<2,3,3,3>},
    }},
  {
    {
      {print_coulomb<3,0,0,0>, print_coulomb<3,0,0,1>, print_coulomb<3,0,0,2>, print_coulomb<3,0,0,3>},
      {print_coulomb<3,0,1,0>, print_coulomb<3,0,1,1>, print_coulomb<3,0,1,2>, print_coulomb<3,0,1,3>},
      {print_coulomb<3,0,2,0>, print_coulomb<3,0,2,1>, print_coulomb<3,0,2,2>, print_coulomb<3,0,2,3>},
      {print_coulomb<3,0,3,0>, print_coulomb<3,0,3,1>, print_coulomb<3,0,3,2>, print_coulomb<3,0,3,3>},
    },
    {
      {print_coulomb<3,1,0,0>, print_coulomb<3,1,0,1>, print_coulomb<3,1,0,2>, print_coulomb<3,1,0,3>},
      {print_coulomb<3,1,1,0>, print_coulomb<3,1,1,1>, print_coulomb<3,1,1,2>, print_coulomb<3,1,1,3>},
      {print_coulomb<3,1,2,0>, print_coulomb<3,1,2,1>, print_coulomb<3,1,2,2>, print_coulomb<3,1,2,3>},
      {print_coulomb<3,1,3,0>, print_coulomb<3,1,3,1>, print_coulomb<3,1,3,2>, print_coulomb<3,1,3,3>},
    },
    {
      {print_coulomb<3,2,0,0>, print_coulomb<3,2,0,1>, print_coulomb<3,2,0,2>, print_coulomb<3,2,0,3>},
      {print_coulomb<3,2,1,0>, print_coulomb<3,2,1,1>, print_coulomb<3,2,1,2>, print_coulomb<3,2,1,3>},
      {print_coulomb<3,2,2,0>, print_coulomb<3,2,2,1>, print_coulomb<3,2,2,2>, print_coulomb<3,2,2,3>},
      {print_coulomb<3,2,3,0>, print_coulomb<3,2,3,1>, print_coulomb<3,2,3,2>, print_coulomb<3,2,3,3>},
    },
    {
      {print_coulomb<3,3,0,0>, print_coulomb<3,3,0,1>, print_coulomb<3,3,0,2>, print_coulomb<3,3,0,3>},
      {print_coulomb<3,3,1,0>, print_coulomb<3,3,1,1>, print_coulomb<3,3,1,2>, print_coulomb<3,3,1,3>},
      {print_coulomb<3,3,2,0>, print_coulomb<3,3,2,1>, print_coulomb<3,3,2,2>, print_coulomb<3,3,2,3>},
      {print_coulomb<3,3,3,0>, print_coulomb<3,3,3,1>, print_coulomb<3,3,3,2>, print_coulomb<3,3,3,3>},
    }},
};
