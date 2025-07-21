#include "Riostream.h"

#include "SmearedNuOscIBDPdf.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
//#include <math.h>
#include "TMath.h"

#include "RooRealVar.h"
#include <algorithm>

ClassImp(SmearedNuOscIBDPdf);

SmearedNuOscIBDPdf::SmearedNuOscIBDPdf(const char *name, const char *title,
    RooAbsReal& xr, RooAbsReal& xInt, RooAbsReal& l,
    RooAbsReal& sin13, RooAbsReal& dm31, RooAbsReal& sin14, RooAbsReal& dm41,
    const RooArgList& elemFracs, const std::vector<const TGraph*> elemSpects, const TGraph* grpXsec,
    const TH2* hResp):
  NuOscIBDPdf(name, title, xInt, l, sin13, dm31, sin14, dm41, elemFracs, elemSpects, grpXsec),
  xr_("xr", "xr", this, xr),
  respMat_(hResp->GetNbinsY()+1, hResp->GetNbinsX()+1)
{
  for (int ix=0; ix<=hResp->GetNbinsX(); ++ix ) {
    binsT_.push_back(hResp->GetXaxis()->GetBinLowEdge(ix+1));
  }
  for (int iy=0; iy<=hResp->GetNbinsY(); ++iy ) {
    binsR_.push_back(hResp->GetYaxis()->GetBinLowEdge(iy+1));
  }
  for (int ix=1; ix<=hResp->GetNbinsX(); ++ix) {
    for (int iy=1; iy<=hResp->GetNbinsY(); ++iy) {
      const double val = hResp->GetBinContent(ix, iy);
      respMat_(iy-1, ix-1) = val; // NOTE: different indexing
    }
  }
}

SmearedNuOscIBDPdf::SmearedNuOscIBDPdf(const char *name, const char *title,
    RooAbsReal& xr, RooAbsReal& xInt, RooAbsReal& l,
    RooAbsReal& sin13, RooAbsReal& dm31, RooAbsReal& sin14, RooAbsReal& dm41,
    const RooArgList& elemFracs, const std::vector<const TGraph*> elemSpects, const TGraph* grpXsec,
    const RooArgList& respPars):
  NuOscIBDPdf(name, title, xInt, l, sin13, dm31, sin14, dm41, elemFracs, elemSpects, grpXsec),
  xr_("xr", "xr", this, xr)
{
  respPars_.add(respPars);
}

SmearedNuOscIBDPdf::SmearedNuOscIBDPdf(const SmearedNuOscIBDPdf& other, const char* name):
  NuOscIBDPdf(other, name),
  xr_("xr", this, other.xr_),
  respMat_(other.respMat_), binsT_(other.binsT_), binsR_(other.binsR_)
{
  respPars_.add(other.respPars_);
}

int SmearedNuOscIBDPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if ( matchArgs(allVars, analVars, xr_) ) return 1;
  return 0;
}

double SmearedNuOscIBDPdf::analyticalIntegral(int code, const char* rangeName) const
{
  R__ASSERT(code == 1);

  return 1.0;
}

double SmearedNuOscIBDPdf::evaluate() const
{
  const RooRealVar* xInt_c = dynamic_cast<const RooRealVar*>(&NuOscIBDPdf::x_.arg());
  if (xInt_c == nullptr) return 0;
  RooRealVar* xInt = const_cast<RooRealVar*>(xInt_c);

  // Find the bin number of the given reco-energy, xr
  const double xr = xr_->getVal();
  if ( xr < binsR_.front() or xr >= binsR_.back() ) return 0;

  auto itr = std::upper_bound(binsR_.begin(), binsR_.end(), xr);
  const size_t idx = std::distance(binsR_.begin(), itr)-1;

  // Then integrate of the true-energy, for the given reco-energy index

  // FIXME: First version picks up only for the bin center.
  //        One needs perform integration within the range, for high freq. oscillations
  //        Maybe we have to consider cashing integrals, checking the parameter updates
  double sumW = 0;
  for (int i=0, n=binsT_.size()-1; i<n; ++i) {
    const double x0 = binsT_[i];
    /*xInt->setVal(x0);
    const double y0 = NoOscIBDPdf::evaluate();
    xx.push_back(x0);
    yy.push_back(y0);*/

    const double x1 = binsT_[i+1];
    const double xC = (x1+x0)/2;
    xInt->setVal(xC);
    const double yC = NuOscIBDPdf::evaluate();

    sumW += yC * respMat_(idx, i);
  }

  return sumW;
}

