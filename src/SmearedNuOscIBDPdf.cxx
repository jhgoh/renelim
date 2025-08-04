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
// Some terms in the integral can be done in advance,
//   P_i = Penv_i - sin14*Posc_i(K41) - sin13*Posc_i(K31)
// The Penv_i can be done at the constructor.

#include "Riostream.h"
#include "RooAbsCategory.h"
#include "RooAbsReal.h"
// #include <math.h>
#include "Math/SpecFunc.h"
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
      xr_("xr", "xr", this, xr), respMat_(hResp->GetNbinsY() + 1, hResp->GetNbinsX() + 1),
      a3Mat_(hResp->GetNbinsY() + 1, hResp->GetNbinsX() + 1),
      b2Mat_(hResp->GetNbinsY() + 1, hResp->GetNbinsX() + 1),
      c1Mat_(hResp->GetNbinsY() + 1, hResp->GetNbinsX() + 1) {
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

  // Cache element fractions for a convenience.
  // This is not really necessary since this part is in the constructor :-)
  std::vector<double> elemFracVals;
  for (int j = 0; j < elemFracs_.getSize(); ++j) {
    const RooAbsReal &elemFrac = static_cast<const RooAbsReal &>(elemFracs_[j]);
    elemFracVals.push_back(elemFrac.getVal());
  }

  // Read the flux and IBD cross section at the bin edges, to be used in pre-computed table
  std::vector<double> specAtBinEdge; // flux at the bin edges
  std::vector<double> xsecAtBinEdge; // IBD cross sections at the bin edges
  for (int i = 0; i < binsT_.size(); ++i) {
    const double x = binsT_[i];

    double spec = 0;
    for (int j = 0; j < elemFracVals.size(); ++j) {
      spec += elemFracVals[j] * interpolate(x, elemSpectsX_[j], elemSpectsY_[j]);
    }

    const double xsec = interpolate(x, ibdXsecX_, ibdXsecY_);

    specAtBinEdge.push_back(spec);
    xsecAtBinEdge.push_back(xsec);
  }

  // Compute the envelopment term
  //   Penv_i = \sum_j R_ij \int_{E_j}^{E_{j+1}) F(E)*S(E)
  // F(E) and S(E) Huber-Mueller and IBD cross section curve,
  //   F(E) = F(E_j) + (dF/dE) * (E-E_j) = F_j + (F_{j+1}-F_j)/(E_{j+1}-E_j) * (E-E_j)
  //   S(E) = S(E_j) + (dS/dE) * (E-E_j) = S_j + (S_{j+1}-S_j)/(E_{j+1}-E_j) * (E-E_j)
  // Therefore, the term inside of the integral becomes 2nd order polynomial
  //   F(E)*S(E) = A*E^2 + B*E + C
  // where the integral is elementary.
  // The A,B,C terms are need in the integral of the oscillating terms,
  //   \int_{E_j}^{E_{j+1}} dE F(E)*S(E) = [ A/3*E^3 + B/2*E^2 + CE ] |_{E_j}^{E_{j+1}}
  for (int i = 0; i < binsR_.size() - 1; ++i) {
    double pEnv = 0;
    for (int j = 0; j < binsT_.size() - 1; ++j) {
      const double e1 = binsT_[j + 1], e0 = binsT_[j];
      const double dE = e1 - e0;
      const double f1 = specAtBinEdge[j + 1], f0 = specAtBinEdge[j];
      const double s1 = xsecAtBinEdge[j + 1], s0 = xsecAtBinEdge[j];
      const double dFdE = (f1 - f0) / dE;
      const double dSdE = (s1 - s0) / dE;

      const double a3 = (dFdE * dSdE) / 3;
      const double b2 = (f0 * dSdE + dFdE * s0 - 2 * dFdE * dSdE * e0) / 2;
      const double c1 = (f0 - dFdE * e0) * (s0 - dSdE * e0);

      const double int0 = e0 * (a3 * e0 * e0 + b2 * e0 + c1);
      const double int1 = e1 * (a3 * e1 * e1 + b2 * e1 + c1);

      pEnv += respMat_(i, j) * (int1 - int0) / dE;

      a3Mat_(i, j) = a3;
      b2Mat_(i, j) = b2;
      c1Mat_(i, j) = c1;
    }

    pEnvs_.push_back(pEnv);
  }
}

SmearedNuOscIBDPdf::SmearedNuOscIBDPdf(const char *name, const char *title, RooAbsReal &xr,
                                       RooAbsReal &xInt, RooAbsReal &l, RooAbsReal &sin13,
                                       RooAbsReal &dm31, RooAbsReal &sin14, RooAbsReal &dm41,
                                       const RooArgList &elemFracs,
                                       const std::vector<const TGraph *> elemSpects,
                                       const TGraph *grpXsec, const RooArgList &respPars)
    : NuOscIBDPdf(name, title, xInt, l, sin13, dm31, sin14, dm41, elemFracs, elemSpects, grpXsec),
      xr_("xr", "xr", this, xr) {
  respPars_.add(respPars);
}

SmearedNuOscIBDPdf::SmearedNuOscIBDPdf(const SmearedNuOscIBDPdf &other, const char *name)
    : NuOscIBDPdf(other, name), xr_("xr", this, other.xr_), respMat_(other.respMat_),
      binsT_(other.binsT_), binsR_(other.binsR_), a3Mat_(other.a3Mat_), b2Mat_(other.b2Mat_),
      c1Mat_(other.c1Mat_), pEnvs_(other.pEnvs_) {
  respPars_.add(other.respPars_);
}

int SmearedNuOscIBDPdf::getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                                              const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, xr_))
    return 1;
  return 0;
}

double SmearedNuOscIBDPdf::analyticalIntegral(int code, const char *rangeName) const {
  R__ASSERT(code == 1);

  return 1.0;
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
  std::vector<double> int31As, int31Bs, int31Cs;
  std::vector<double> int41As, int41Bs, int41Cs;
  for (int j = 0, n = binsT_.size(); j < n; ++j) {
    const double e0 = binsT_[j];

    // Disappearance term for nu_3
    const double kk31e0 = 2 * k31 / e0;
    const double si31e0 = ROOT::Math::sinint(kk31e0);
    const double ci31e0 = ROOT::Math::cosint(kk31e0);
    const double sin31e0 = std::sin(kk31e0), cos31e0 = std::cos(kk31e0);

    const double int31A0 = 2 * k31 * k31 * k31 * si31e0 +
                           (k31 * k31 * e0 - e0 * e0 * e0 / 2) * cos31e0 + e0 * e0 * e0 / 2 +
                           k31 / 2 * e0 * e0 * sin31e0;
    const double int31B0 =
        -2 * k31 * k31 * ci31e0 + k31 * e0 * sin31e0 - e0 * e0 / 2 * cos31e0 + e0 * e0 / 2;
    const double int31C0 = -k31 * si31e0 + e0 / 2 * (1 - cos31e0);

    // Disappearance term for nu_4
    const double kk41e0 = 2 * k41 / e0;
    const double si41e0 = ROOT::Math::sinint(kk41e0);
    const double ci41e0 = ROOT::Math::cosint(kk41e0);
    const double sin41e0 = std::sin(kk41e0), cos41e0 = std::cos(kk41e0);

    const double int41A0 = 2 * k41 * k41 * k41 * si41e0 +
                           (k41 * k41 * e0 - e0 * e0 * e0 / 2) * cos41e0 + e0 * e0 * e0 / 2 +
                           k41 / 2 * e0 * e0 * sin41e0;
    const double int41B0 =
        -2 * k41 * k41 * ci41e0 + k41 * e0 * sin41e0 - e0 * e0 / 2 * cos41e0 + e0 * e0 / 2;
    const double int41C0 = -k41 * si41e0 + e0 / 2 * (1 - cos41e0);

    // Keep each terms
    int31As.push_back(e0 <= 0 ? 0 : int31A0);
    int31Bs.push_back(e0 <= 0 ? 0 : int31B0);
    int31Cs.push_back(e0 <= 0 ? 0 : int31C0);

    int41As.push_back(e0 <= 0 ? 0 : int41A0);
    int41Bs.push_back(e0 <= 0 ? 0 : int41B0);
    int41Cs.push_back(e0 <= 0 ? 0 : int41C0);
  }

  double pOsc31 = 0, pOsc41 = 0;
  for (int j = 0, nn = binsT_.size() - 1; j < nn; ++j) {
    const double dE = binsT_[j + 1] - binsT_[j];
    const double a3 = a3Mat_(idx, j);
    const double b2 = b2Mat_(idx, j);
    const double c1 = c1Mat_(idx, j);

    const double int31 = a3 * (int31As[j + 1] - int31As[j]) + b2 * (int31Bs[j + 1] - int31Bs[j]) +
                         c1 * (int31Cs[j + 1] - int31Cs[j]);
    const double int41 = a3 * (int41As[j + 1] - int41As[j]) + b2 * (int41Bs[j + 1] - int41Bs[j]) +
                         c1 * (int41Cs[j + 1] - int41Cs[j]);
    pOsc31 += respMat_(idx, j) * int31 / dE;
    pOsc41 += respMat_(idx, j) * int41 / dE;
  }

  return std::max(0.0, pEnvs_[idx] - sin13 * pOsc31 - sin14 * pOsc41);
  // return std::abs(pEnvs_[idx] - sin13*pOsc31 - sin14*pOsc41);
}
