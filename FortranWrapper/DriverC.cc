#include <iostream>

// Utility functions
extern "C" 
{
  void initialiseevolution_(int const& PerturbativeOrder, double const& mu0, double const& AlphaQCDRef, double const& MuAlphaQCDRef, double const& mc, double const& mb, double* p);
  void evolvetransversities_(double const& x, double const& mu, double* xf);
}

//_____________________________________________________________________________
int main()
{
  // Initialisation
  const int PerturbativeOrder = 0;
  const double mu0 = 1;
  const double AlphaQCDRef = 0.13939;
  const double MuAlphaQCDRef = 91.1876;
  const double mc = 1.4;
  const double mb = 4.5;
  double *p = new double[10];
  for (int i = 0; i < 10; i++)
    p[i] = 10;
  initialiseevolution_(PerturbativeOrder, mu0, AlphaQCDRef, MuAlphaQCDRef, mc, mb, p);

  // Evolution
  const double x  = 0.1;
  const double mu = 10;
  double *xf = new double[13];
  evolvetransversities_(x, mu, xf);
  for(int i = 0; i < 13; i++)
    std::cout << std::scientific << i << "  " << xf[i] << std::endl;
  delete[] xf;
  return 0;
}
