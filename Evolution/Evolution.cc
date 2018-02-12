//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/grid.h>
#include <apfel/alphaqcd.h>
#include <apfel/tabulateobject.h>
#include <apfel/dglapbuilder.h>
#include <apfel/lhtoypdfs.h>
#include <apfel/rotations.h>

using namespace apfel;
using namespace std;

int main()
{
  // x-space grid
  const Grid g{{SubGrid{100,1e-5,3}, SubGrid{60,1e-1,3}, SubGrid{50,6e-1,3}, SubGrid{50,8e-1,3}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vectors of masses and thresholds
  const vector<double> Masses = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 1;

  // Running coupling
  const double AlphaQCDRef = 0.35;
  const double MuAlphaQCDRef = sqrt(2);
  AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Masses, PerturbativeOrder};
  const TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Initialize QCD evolution objects
  //const auto DglapObj = InitializeDglapObjectsQCDtrans(g, Masses);  // Space-like (PDFs)
  const auto DglapObj = InitializeDglapObjectsQCDTtrans(g, Masses);  // Time-like (FFs)

  // Construct the DGLAP objects
  auto EvolvedPDFs = BuildDglap(DglapObj, LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const TabulateObject<Set<Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Final scale
  double mu = 100;

  // Print results
  cout << scientific;

  const vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  cout << "\n   x    "
       << "   u-ubar   "
       << "   d-dbar   "
       << " 2(ubr+dbr) "
       << "   c+cbar   "
       << "   gluon    "
       << endl;

  for (auto const& x : xlha)
    {
      const auto DistMap = QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x,mu));
      cout.precision(1);
      cout << x;
      cout.precision(4);
      cout << "  " <<
	DistMap.at(2) - DistMap.at(-2) << "  " <<
	DistMap.at(1) - DistMap.at(-1) << "  " <<
	2 * ( DistMap.at(-2) + DistMap.at(-1) ) << "  " <<
	DistMap.at(4) + DistMap.at(-4) << "  " <<
	DistMap.at(0) << "  "
	   << endl;
    }
  cout << "      " << endl;

  return 0;
}
