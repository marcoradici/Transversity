#include <apfel/apfelxx.h>

// Parameters
double pars[10] = {0};

// Evolution parameters
const double mc = 1.4;
const double mb = 4.5;
const int PerturbativeOrder = 0;   // 0 = LO , 1 = NLO
const double mu0 = 1;   // Q0
const double AlphaQCDRef = 0.13939;   // alpha_S(MZ) from MSTW08: LO 0.13939 , NLO 0.12018
const double MuAlphaQCDRef = 91.1876;

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
  // Running coupling object
  apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, {0, 0, 0, mc, mb}, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 101, 3};
  double as(double const& mu) { return Alphas.Evaluate(mu); };

  // Evolution object
  const apfel::Grid g{{{100, 1e-5, 3}, {60, 1e-1, 3}, {50, 6e-1, 3}, {50, 8e-1, 3}}};
  const auto DglapObj = apfel::InitializeDglapObjectsQCDtrans(g, {0, 0, 0, mc, mb});  // Space-like (PDFs)
//  const auto DglapObj = apfel::InitializeDglapObjectsQCDTtrans(*g, Masses);  // Time-like (FFs)
  std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> TabulatedPDFs;

  //_____________________________________________________________________________
  void initialiseevolution_(double* p)
  {
    TabulatedPDFs = NULL;

    // Set silent mode for APFEL++
    apfel::SetVerbosityLevel(0);

    // Allocate parameters
    for (int i = 0; i < 10; i++)
      pars[i] = p[i];

    // Set initial-scale distributions
    std::unique_ptr<apfel::Dglap<apfel::Distribution>> EvolvedPDFs = apfel::BuildDglap(DglapObj, Dists, mu0, PerturbativeOrder, as);

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
