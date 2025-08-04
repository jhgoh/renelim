#include "SmearedNuOscIBDPdf.h"

// We start from the reduced reco-energy formula:
//   P_i = \sum_j R_ij \int_{E_j}^{E_{j+1}) F(E)*S(E) (1-sin14*sin^2(K41/E)-sin13*...)
//     K41 = 1.27*dm41*L
//     R_ij = response matrix
// F(E) and S(E) Huber-Mueller and IBD cross section curve,
//   F(E) = F(E_j) + (dF/dE) * (E-E_j) = F_j + (F_{j+1}-F_j)/(E_{j+1}-E_j) * (E-E_j)
//   S(E) = S(E_j) + (dS/dE) * (E-E_j) = S_j + (S_{j+1}-S_j)/(E_{j+1}-E_j) * (E-E_j)
// Note that more precise evaluation can be done by introducing finer binning
// according to the data points of F(E)'s and S(E)'s original papers.
//
// Note: Some terms in the integral can be done in advance,
//   P_i = Penv_i - sin14*Posc_i(K41) - sin13*Posc_i(K31)
// The Penv_i can be done at the constructor, which was considered in the first-working version
// but the current version do the computation in every calls for code-reusability of
// the mother class, NuOscIBDPdf::subIntegral()

#include "Riostream.h"
#include "RooAbsCategory.h"
#include "RooAbsReal.h"
// #include <math.h>
#include "TMath.h"

#include "RooRealVar.h"
#include <algorithm>

ClassImp(SmearedNuOscIBDPdf);

SmearedNuOscIBDPdf::SmearedNuOscIBDPdf(const char *name, const char *title, RooAbsReal &xr,
                                       RooAbsReal &xInt, RooAbsReal &l, RooAbsReal &sin13,
                                       RooAbsReal &dm31, RooAbsReal &sin14, RooAbsReal &dm41,
                                       const RooArgList &elemFracs,
                                       const std::vector<const TGraph *> elemSpects,
                                       const TGraph *grpXsec, const TH2 *hResp)
    : NuOscIBDPdf(name, title, xInt, l, sin13, dm31, sin14, dm41, elemFracs, elemSpects, grpXsec),
      xr_("xr", "xr", this, xr), respMat_(hResp->GetNbinsY() + 1, hResp->GetNbinsX() + 1) {
  // Load the response matrix and its binning
  for (int ix = 0; ix <= hResp->GetNbinsX(); ++ix) {
    binsT_.push_back(hResp->GetXaxis()->GetBinLowEdge(ix + 1));
  }
  for (int iy = 0; iy <= hResp->GetNbinsY(); ++iy) {
    binsR_.push_back(hResp->GetYaxis()->GetBinLowEdge(iy + 1));
  }
  for (int ix = 1; ix <= hResp->GetNbinsX(); ++ix) { // Etrue axis
    double sumE = 0;
    for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy) { // Ereco axis
      const double val = hResp->GetBinContent(ix, iy);
      sumE += val;
    }
    if (sumE == 0)
      continue;
    for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy) { // Ereco axis
      const double val = hResp->GetBinContent(ix, iy);
      respMat_(iy - 1, ix - 1) = val / sumE; // NOTE: different indexing
    }
  }
}

SmearedNuOscIBDPdf::SmearedNuOscIBDPdf(const SmearedNuOscIBDPdf &other, const char *name)
    : NuOscIBDPdf(other, name), xr_("xr", this, other.xr_), respMat_(other.respMat_),
      binsT_(other.binsT_), binsR_(other.binsR_) {
}

int SmearedNuOscIBDPdf::getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                                              const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, xr_))
    return 1;
  return 0;
}

double SmearedNuOscIBDPdf::analyticalIntegral(int code, const char *rangeName) const {
  R__ASSERT(code == 1);

  return NuOscIBDPdf::analyticalIntegral(code, rangeName);
}

double SmearedNuOscIBDPdf::evaluate() const {
  // Find the bin number of the given reco-energy, xr
  const double xr = xr_->getVal();
  if (xr < binsR_.front() or xr >= binsR_.back())
    return 0;

  auto itr = std::upper_bound(binsR_.begin(), binsR_.end(), xr);
  const size_t idx = std::distance(binsR_.begin(), itr) - 1;

  // Compute the detector resolution convoluted, binned, energy distribution.
  // Details can be found in the header file.
  const double l = l_->getVal();

  const double sin13 = sin13_->getVal();
  const double dm31 = dm31_->getVal();
  const double sin14 = sin14_->getVal();
  const double dm41 = dm41_->getVal();

  const double k31 = 1.27 * dm31 * l;
  const double k41 = 1.27 * dm41 * l;

  std::vector<double> elemFracs;
  for (int i = 0; i < elemFracs_.getSize(); ++i) {
    const RooAbsReal &elemFrac = static_cast<const RooAbsReal &>(elemFracs_[i]);
    elemFracs.push_back(elemFrac.getVal());
  }

  double sumW = 0;
  for (int j = 0, nn = binsT_.size()-1; j < nn; ++j) {
    const double e0 = binsT_[j], e1 = binsT_[j + 1];

    double f0 = 0, f1 = 0;
    for (int j = 0; j < elemFracs.size(); ++j) {
      f0 += elemFracs[j] * interpolate(e0, elemSpectsX_[j], elemSpectsY_[j]);
      f1 += elemFracs[j] * interpolate(e1, elemSpectsX_[j], elemSpectsY_[j]);
    }
    const double s0 = interpolate(e0, ibdXsecX_, ibdXsecY_);
    const double s1 = interpolate(e1, ibdXsecX_, ibdXsecY_);

    sumW += respMat_(idx, j) * subIntegral(e0, e1, f0, f1, s0, s1, sin13, k31, sin14, k41);
  }

  return sumW;
}
