// Implementation of the smeared IBD oscillation spectrum. The starting point
// is the binned formula
//   P_i = \sum_j R_{ij} \int_{E_j}^{E_{j+1}} F(E) S(E)
//         [1 - \sin^2 2\theta_{14}\sin^2(K_{41}/E)
//          - \sin^2 2\theta_{13}\sin^2(K_{31}/E)],
// where $R_{ij}$ is the response matrix. Both the flux F(E) and cross section
// S(E) are approximated as linear functions within each bin.
// A more precise evaluation can be obtained by refining the binning according
// to the original data points of F(E) and S(E).
//
// The reconstructed-energy spectrum is cached across evaluate() calls.
// fillCache() computes the full spectrum once per unique set of oscillation
// parameters using NuOscIBDPdf::subIntegral(), so plotOn() only triggers one
// cache fill per drawn curve instead of one full evaluation per plotting bin.

#include "Riostream.h"
#include "BinnedNuOscIBDPdf.h"
#include "RooAbsCategory.h"
#include "RooAbsReal.h"
#include "TMath.h"
#include "TGraph.h"

#include "RooRealVar.h"
#include <algorithm>

ClassImp(BinnedNuOscIBDPdf);

BinnedNuOscIBDPdf::BinnedNuOscIBDPdf(const char *name, const char *title, RooAbsReal &xr,
                                       RooAbsReal &xInt, RooAbsReal &l, RooAbsReal &sin13,
                                       RooAbsReal &dm31, RooAbsReal &sin14, RooAbsReal &dm41,
                                       const RooArgList &elemFracs,
                                       const std::vector<std::vector<double>> &elemSpectsX,
                                       const std::vector<std::vector<double>> &elemSpectsY,
                                       const std::vector<double> &ibdXsecX,
                                       const std::vector<double> &ibdXsecY,
                                       const TH1* hBaseline, const TH2 *hResp)
    : NuOscIBDPdf(name, title, xInt, l, sin13, dm31, sin14, dm41, elemFracs,
                  elemSpectsX, elemSpectsY, ibdXsecX, ibdXsecY),
      xr_("xr", "xr", this, xr), respMat_(hResp->GetNbinsY() + 1, hResp->GetNbinsX() + 1) {
  // Load the baseline smearing kernel
  //l_->setVal(hBaseline->GetMean());
  for (int ix = 1; ix <= hBaseline->GetNbinsX(); ++ix ) {
    const double lw = hBaseline->GetBinContent(ix);
    if ( lw == 0 ) {
      if ( lws_.empty() )
        continue; // left-strip, skip empty entries front
      else
        break; // right-strip, skip empty entries back - caveat: ignoring tails
    }
    ls_.push_back(hBaseline->GetXaxis()->GetBinCenter(ix));
    lws_.push_back(hBaseline->GetBinContent(ix));
  }

  // Load the response matrix and normalise each true-energy slice.
  for (int ix = 0; ix <= hResp->GetNbinsX(); ++ix) {
    binsT_.push_back(hResp->GetXaxis()->GetBinLowEdge(ix + 1));
  }
  for (int iy = 0; iy <= hResp->GetNbinsY(); ++iy) {
    binsR_.push_back(hResp->GetYaxis()->GetBinLowEdge(iy + 1));
  }
  for (int ix = 1; ix <= hResp->GetNbinsX(); ++ix) { // true-E axis
    double sumE = 0;
    for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy) { // reco-E axis
      const double val = hResp->GetBinContent(ix, iy);
      sumE += val;
    }
    if (sumE == 0)
      continue;
    for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy) { // reco-E axis
      const double val = hResp->GetBinContent(ix, iy);
      // NOTE: different indexing: respMat_(iy - 1, ix - 1) stores (reco, true)
      // order, opposite to the (x, y) convention used by TH2D
      respMat_(iy - 1, ix - 1) = val / sumE; // normalised
    }
  }
}

BinnedNuOscIBDPdf::BinnedNuOscIBDPdf(const char *name, const char *title, RooAbsReal &xr,
                                       RooAbsReal &xInt, RooAbsReal &l, RooAbsReal &sin13,
                                       RooAbsReal &dm31, RooAbsReal &sin14, RooAbsReal &dm41,
                                       const RooArgList &elemFracs,
                                       const std::vector<const TGraph *> elemSpects,
                                       const TGraph *grpXsec, const TH1* hBaseline, const TH2 *hResp)
    : BinnedNuOscIBDPdf(name, title, xr, xInt, l, sin13, dm31, sin14, dm41, elemFracs,
                         xFromGraphs(elemSpects), yFromGraphs(elemSpects),
                         xFromGraph(grpXsec), yFromGraph(grpXsec), hBaseline, hResp) {}

BinnedNuOscIBDPdf::BinnedNuOscIBDPdf(const BinnedNuOscIBDPdf &other, const char *name)
    : NuOscIBDPdf(other, name), xr_("xr", this, other.xr_), respMat_(other.respMat_),
      binsT_(other.binsT_), binsR_(other.binsR_),
      ls_(other.ls_), lws_(other.lws_) {}

int BinnedNuOscIBDPdf::getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                                              const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, xr_))
    return 1;
  return 0;
}

double BinnedNuOscIBDPdf::analyticalIntegral(int code, const char *rangeName) const {
  R__ASSERT(code == 1);

  return NuOscIBDPdf::analyticalIntegral(code, rangeName);
}

// ---------------------------------------------------------------------------
// fillCache() – precompute the full reconstructed-energy spectrum for the
// current oscillation parameters and store it in cache_[iR].
//
// The original evaluate() did O(nL * nE) work per call, and plotOn() calls
// evaluate() once per plotting bin (1024 by default), so the total cost was
// O(n_curves * 1024 * nL * nE).  Here we observe that:
//   1. The flux/cross-section endpoint values (f0,f1,s0,s1) do not depend on
//      the baseline L, so they need only be computed once per true-energy bin.
//   2. The full reconstructed-energy spectrum can be written as
//        P[iR] = sum_iE  respMat(iR,iE) * trueSpec[iE]
//      where trueSpec[iE] = sum_iL lw[iL] * subIntegral(iE, l[iL]).
//      Factoring the loops this way lets us compute P[iR] for *all* iR in a
//      single pass of cost O(nE * nL + nE * nR), and then each subsequent
//      evaluate() call is a O(1) cache lookup.
// ---------------------------------------------------------------------------
void BinnedNuOscIBDPdf::fillCache() const {
  const double sin13 = sin13_->getVal();
  const double dm31  = dm31_->getVal();
  const double sin14 = sin14_->getVal();
  const double dm41  = dm41_->getVal();

  // Collect element fractions once.
  std::vector<double> elemFracs;
  elemFracs.reserve(elemFracs_.getSize());
  for (int iF = 0; iF < elemFracs_.getSize(); ++iF) {
    const RooAbsReal &ef = static_cast<const RooAbsReal &>(elemFracs_[iF]);
    elemFracs.push_back(ef.getVal());
  }

  const int nE = static_cast<int>(binsT_.size()) - 1;

  // Step 1 – precompute flux * cross-section endpoint values for each
  // true-energy bin.  These are independent of the baseline L.
  std::vector<double> f0v(nE, 0.0), f1v(nE, 0.0), s0v(nE), s1v(nE);
  for (int iE = 0; iE < nE; ++iE) {
    const double e0 = binsT_[iE], e1 = binsT_[iE + 1];
    s0v[iE] = interpolate(e0, ibdXsecX_, ibdXsecY_);
    s1v[iE] = interpolate(e1, ibdXsecX_, ibdXsecY_);
    for (size_t iF = 0; iF < elemFracs.size(); ++iF) {
      f0v[iE] += elemFracs[iF] * interpolate(e0, elemSpectsX_[iF], elemSpectsY_[iF]);
      f1v[iE] += elemFracs[iF] * interpolate(e1, elemSpectsX_[iF], elemSpectsY_[iF]);
    }
  }

  // Step 2 – compute the baseline-smeared true-energy spectrum.
  std::vector<double> trueSpec(nE, 0.0);
  for (size_t iL = 0; iL < ls_.size(); ++iL) {
    const double l  = ls_[iL];
    const double lw = lws_[iL];
    const double k31 = 1.27 * dm31 * l;
    const double k41 = 1.27 * dm41 * l;
    for (int iE = 0; iE < nE; ++iE) {
      trueSpec[iE] += lw * subIntegral(binsT_[iE], binsT_[iE + 1],
                                       f0v[iE], f1v[iE], s0v[iE], s1v[iE],
                                       sin13, k31, sin14, k41);
    }
  }

  // Step 3 – apply the response matrix to obtain the reconstructed-energy
  // spectrum for every reco bin.
  const int nR = static_cast<int>(binsR_.size()) - 1;
  cache_.assign(nR, 0.0);
  for (int iR = 0; iR < nR; ++iR) {
    for (int iE = 0; iE < nE; ++iE) {
      cache_[iR] += respMat_(iR, iE) * trueSpec[iE];
    }
  }

  cachedSin13_ = sin13;
  cachedDm31_  = dm31;
  cachedSin14_ = sin14;
  cachedDm41_  = dm41;
}

double BinnedNuOscIBDPdf::evaluate() const {
  // Determine the reconstructed-energy bin.
  const double xr = xr_->getVal();
  if (xr < binsR_.front() || xr >= binsR_.back())
    return 0;

  // Rebuild the cache whenever the oscillation parameters have changed.
  const double sin13 = sin13_->getVal();
  const double dm31  = dm31_->getVal();
  const double sin14 = sin14_->getVal();
  const double dm41  = dm41_->getVal();

  if (cache_.empty() ||
      sin13 != cachedSin13_ || dm31 != cachedDm31_ ||
      sin14 != cachedSin14_ || dm41 != cachedDm41_) {
    fillCache();
  }

  auto itr = std::upper_bound(binsR_.begin(), binsR_.end(), xr);
  const size_t idx = std::distance(binsR_.begin(), itr) - 1;

  return cache_[idx];
}
