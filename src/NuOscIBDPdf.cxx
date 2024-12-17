#include "Riostream.h" 

#include "NuOscIBDPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
//#include <math.h> 
#include "TMath.h" 

#include <algorithm>

ClassImp(NuOscIBDPdf); 

NuOscIBDPdf::NuOscIBDPdf(const char *name, const char *title, RooAbsReal& x,
    RooAbsReal& l, RooAbsReal& sin13, RooAbsReal& dm13, RooAbsReal& sin14, RooAbsReal& dm14,
    RooAbsReal& fU235, RooAbsReal& fPu239, RooAbsReal& fPu241,
    const TGraph* grpU235, const TGraph* grpPu239, const TGraph* grpPu241,
    const TGraph* grpXsec):
  RooAbsPdf(name, title),
  x_("x", "x", this, x),
  l_("l", "l", this, l),
  sin13_("sin13", "sin13", this, sin13), dm13_("dm13", "dm13", this, dm13),
  sin14_("sin14", "sin14", this, sin14), dm14_("dm14", "dm14", this, dm14),
  fU235_("fU235", "fU235", this, fU235),
  fU238_("fU238", "fU238", this, fU235), // Note copy dummy value from U235
  fPu239_("fPu239", "fPu239", this, fPu239), fPu241_("fPu241", "fPu241", this, fPu241),
  hasU238_(false)
{
  loadFromTGraph(grpU235, xx_U235_, yy_U235_);
  loadFromTGraph(nullptr, xx_U238_, yy_U238_);
  loadFromTGraph(grpPu239, xx_Pu239_, yy_Pu239_);
  loadFromTGraph(grpPu241, xx_Pu241_, yy_Pu241_);
  loadFromTGraph(grpXsec, xx_Xsec_, yy_Xsec_);
}

NuOscIBDPdf::NuOscIBDPdf(const char *name, const char *title, RooAbsReal& x,
    RooAbsReal& l, RooAbsReal& sin13, RooAbsReal& dm13, RooAbsReal& sin14, RooAbsReal& dm14,
    RooAbsReal& fU235, RooAbsReal& fU238, RooAbsReal& fPu239, RooAbsReal& fPu241,
    const TGraph* grpU235, const TGraph* grpU238, const TGraph* grpPu239, const TGraph* grpPu241,
    const TGraph* grpXsec):
  RooAbsPdf(name, title),
  x_("x", "x", this, x),
  l_("l", "l", this, l),
  sin13_("sin13", "sin13", this, sin13), dm13_("dm13", "dm13", this, dm13),
  sin14_("sin14", "sin14", this, sin14), dm14_("dm14", "dm14", this, dm14),
  fU235_("fU235", "fU235", this, fU235), fU238_("fU238", "fU238", this, fU238),
  fPu239_("fPu239", "fPu239", this, fPu239), fPu241_("fPu241", "fPu241", this, fPu241),
  hasU238_(true)
{
  loadFromTGraph(grpU235, xx_U235_, yy_U235_);
  loadFromTGraph(grpU238, xx_U238_, yy_U238_);
  loadFromTGraph(grpPu239, xx_Pu239_, yy_Pu239_);
  loadFromTGraph(grpPu241, xx_Pu241_, yy_Pu241_);
  loadFromTGraph(grpXsec, xx_Xsec_, yy_Xsec_);
}

NuOscIBDPdf::NuOscIBDPdf(const NuOscIBDPdf& other, const char* name):
  RooAbsPdf(other, name),
  x_("x", this, other.x_),
  l_("l", this, other.l_),
  sin13_("sin13", this, other.sin13_), dm13_("dm13", this, other.dm13_),
  sin14_("sin14", this, other.sin14_), dm14_("dm14", this, other.dm14_),
  fU235_("fU235", this, other.fU235_), fU238_("fU238", this, other.fU238_),
  fPu239_("fPu239", this, other.fPu239_), fPu241_("fPu241", this, other.fPu241_),
  xx_U235_(other.xx_U235_), yy_U235_(other.yy_U235_),
  xx_U238_(other.xx_U238_), yy_U238_(other.yy_U238_),
  xx_Pu239_(other.xx_Pu239_), yy_Pu239_(other.yy_Pu239_),
  xx_Pu241_(other.xx_Pu241_), yy_Pu241_(other.yy_Pu241_),
  xx_Xsec_(other.xx_Xsec_), yy_Xsec_(other.yy_Xsec_),
  hasU238_(other.hasU238_)
{
}

void NuOscIBDPdf::loadFromTGraph(const TGraph* grp, std::vector<double>& xx, std::vector<double>& yy)
{
  if ( !grp ) {
    xx = {{0, 1}};
    yy = {{0, 0}};
    return;
  }

  xx.clear();
  yy.clear();

  for ( int i=0, n=grp->GetN(); i<n; ++i ) {
    xx.push_back(grp->GetX()[i]);
    yy.push_back(grp->GetY()[i]);
  }
}

double NuOscIBDPdf::interpolate(const double x, const std::vector<double>& xx, const std::vector<double>& yy) const
{
  auto itr = std::lower_bound(xx.begin(), xx.end(), x);
  if ( itr == xx.begin() or itr == xx.end() ) {
    if ( x == xx[0] ) return yy[0]; // At boundary
    else return 0; // Ouside of the boundary
  }

  const size_t idx = itr - xx.begin();
  const double x2 = *itr, x1 = *(itr-1);
  const double y2 = yy[idx], y1 = yy[idx-1];

  return y1 + (y2-y1)/(x2-x1)*(x-x1);
}

double NuOscIBDPdf::evaluate() const 
{ 
  const double x = x_->getVal();
  const double l = l_->getVal();
  if ( x <= 0 ) return 0; // safeguard for unphysical range

  const double sin13 = sin13_->getVal();
  const double dm13 = dm13_->getVal();
  const double sinD13 = std::sin(dm13*1.27*l/x); // 1.27 factor comes from 1/4 hbar c
  const double prob13 = 1 - sin13 * sinD13 * sinD13;

  const double sin14 = sin14_->getVal();
  const double dm14 = dm14_->getVal();
  const double sinD14 = std::sin(dm14*1.27*l/x); // 1.27 factor comes from 1/4 hbar c
  const double prob14 = prob13 - sin14 * sinD14 * sinD14;

  const double specU235 = fU235_->getVal() * interpolate(x, xx_U235_, yy_U235_);
  const double specU238 = hasU238_ ? fU238_->getVal() * interpolate(x, xx_U238_, yy_U238_) : 0;
  const double specPu239 = fPu239_->getVal() * interpolate(x, xx_Pu239_, yy_Pu239_);
  const double specPu241 = fPu241_->getVal() * interpolate(x, xx_Pu241_, yy_Pu241_);

  const double spec = specU235 + specU238 + specPu239 + specPu241;
  const double xsec = interpolate(x, xx_Xsec_, yy_Xsec_);

  return spec * xsec * prob14;
}

