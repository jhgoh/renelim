#include "Riostream.h" 

#include "NuOscIBDPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
//#include <math.h> 
#include "TMath.h" 

#include <algorithm>

ClassImp(NuOscIBDPdf); 

NuOscIBDPdf::NuOscIBDPdf(const char *name, const char *title, RooAbsReal& x,
    RooAbsReal& l, RooAbsReal& sin13, RooAbsReal& dm31, RooAbsReal& sin14, RooAbsReal& dm41,
    const RooArgList& elemFracs, const std::vector<const TGraph*> elemSpects,
    const TGraph* grpXsec):
  RooAbsPdf(name, title),
  x_("x", "x", this, x),
  l_("l", "l", this, l),
  sin13_("sin13", "sin13", this, sin13), dm31_("dm31", "dm31", this, dm31),
  sin14_("sin14", "sin14", this, sin14), dm41_("dm41", "dm41", this, dm41)
{
  assert(elemFracs.getSize() == elemSpects.size());

  elemFracs_.add(elemFracs);
  for ( int i=0; i<elemFracs_.getSize(); ++i ) {
    elemSpectsX_.push_back({});
    elemSpectsY_.push_back({});
    loadFromTGraph(elemSpects[i], elemSpectsX_[i], elemSpectsY_[i]);
  }

  loadFromTGraph(grpXsec, ibdXsecX_, ibdXsecY_);
}

NuOscIBDPdf::NuOscIBDPdf(const NuOscIBDPdf& other, const char* name):
  RooAbsPdf(other, name),
  x_("x", this, other.x_),
  l_("l", this, other.l_),
  sin13_("sin13", this, other.sin13_), dm31_("dm31", this, other.dm31_),
  sin14_("sin14", this, other.sin14_), dm41_("dm41", this, other.dm41_),
  elemSpectsX_(other.elemSpectsX_), elemSpectsY_(other.elemSpectsY_),
  ibdXsecX_(other.ibdXsecX_), ibdXsecY_(other.ibdXsecY_)
{
  elemFracs_.add(other.elemFracs_);
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
    else return 0; // Outside of the boundary
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
  const double dm31 = dm31_->getVal();
  const double sinD13 = std::sin(dm31*1.27*l/x); // 1.27 factor comes from 1/4 hbar c
  const double prob13 = 1 - sin13 * sinD13 * sinD13;

  const double sin14 = sin14_->getVal();
  const double dm41 = dm41_->getVal();
  const double sinD14 = std::sin(dm41*1.27*l/x); // 1.27 factor comes from 1/4 hbar c
  const double prob14 = prob13 - sin14 * sinD14 * sinD14;

  double spec = 0;
  for ( int i=0; i<elemFracs_.getSize(); ++i ) {
    const RooAbsReal& elemFrac = static_cast<const RooAbsReal&>(elemFracs_[i]);
    spec += elemFrac.getVal() * interpolate(x, elemSpectsX_[i], elemSpectsY_[i]);
  }
  const double xsec = interpolate(x, ibdXsecX_, ibdXsecY_);

  return spec * xsec * prob14;
}

