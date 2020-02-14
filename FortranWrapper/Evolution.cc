//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

// Evolution object
std::unique_ptr<const apfel::Grid> g;
std::unique_ptr<const apfel::TabulateObject<apfel::Set<apfel::Distribution>>> TabulatedPDFs;

// Parameters
double pars[10] = {0};

// Initial scale functions
extern"C"
{
  double xh1up_  (double *x, double *norm, double *xpow, double *ceb1, double *ceb2, double *ceb3);
  double xh1down_(double *x, double *norm, double *xpow, double *ceb1, double *ceb2, double *ceb3);
}

//_________________________________________________________________________________
std::map<int, double> Dists(double const& x, double const&) 
{
  double xx = x;

  // Call all functions only once.
  const double dnv = xh1down_(&xx, &pars[0], &pars[1], &pars[2], &pars[3], &pars[4]);
  const double upv = xh1up_  (&xx, &pars[5], &pars[6], &pars[7], &pars[8], &pars[9]);


  // Construct QCD evolution basis conbinations.
  double const Gluon   = 0;
  double const Singlet = dnv + upv;
  double const T3      = upv - dnv;
  double const T8      = upv + dnv;
  double const T15     = upv + dnv;
  double const Valence = upv + dnv;
  double const V3      = upv - dnv;

  // Fill in map in the QCD evolution basis.
  std::map<int,double> QCDEvMap;
  QCDEvMap[0]  = Gluon;
  QCDEvMap[1]  = Singlet;
  QCDEvMap[2]  = Valence;
  QCDEvMap[3]  = T3;
  QCDEvMap[4]  = V3;
  QCDEvMap[5]  = T8;
  QCDEvMap[6]  = T15;
  QCDEvMap[7]  = Singlet;
  QCDEvMap[8]  = Valence;
  QCDEvMap[9]  = Singlet;
  QCDEvMap[10] = Valence;
  QCDEvMap[11] = Singlet;
  QCDEvMap[12] = Valence;
  return QCDEvMap;
}

extern "C" 
{
  //_____________________________________________________________________________
  void initialiseevolution_(int const& PerturbativeOrder, double const& mu0, double const& AlphaQCDRef, double const& MuAlphaQCDRef, double const& mc, double const& mb, double* p)
  {
    // Allocate parameters
    for (int i = 0; i < 10; i++)
      pars[i] = p[i];

    // Heavy quark mass vector
    const std::vector<double> Masses{0, 0, 0, mc, mb};

    // Running coupling object
    apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Masses, PerturbativeOrder};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 101, 3};
    const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

    // x-space grid
    g = std::unique_ptr<const apfel::Grid>(new apfel::Grid{{{100,1e-5,3}, {60,1e-1,3}, {50,6e-1,3}, {50,8e-1,3}}});

    // Initialize evolution objects
    const auto DglapObj = apfel::InitializeDglapObjectsQCDtrans(*g, Masses);  // Space-like (PDFs)
    //const auto DglapObj = apfel::InitializeDglapObjectsQCDTtrans(*g, Masses);  // Time-like (FFs)

    // Construct the DGLAP objects
    const auto EvolvedPDFs = apfel::BuildDglap(DglapObj, Dists, mu0, PerturbativeOrder, as);

    // Tabulate distributions
    TabulatedPDFs = std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>>(new apfel::TabulateObject<apfel::Set<apfel::Distribution>>{*EvolvedPDFs, 50, 1, 100, 3});
  }

  //_____________________________________________________________________________
  void evolvetransversities_(double const& x, double const& mu, double* xf)
  {
    // Evolve and rotate distributions
    for (auto const& f : apfel::QCDEvToPhys(TabulatedPDFs->EvaluateMapxQ(x, mu)))
      if (f.first != 21)
	xf[f.first + 6] = f.second;
  }
}
