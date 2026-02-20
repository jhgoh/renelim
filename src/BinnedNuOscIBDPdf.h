#ifndef BinnedNuOscIBDPdf_H
#define BinnedNuOscIBDPdf_H

// BinnedNuOscIBDPdf: energy spectrum of sterile neutrinos including detector resolution.
//
// A simple convolution cannot be used here because the detector resolution
// depends on energy. Instead the smeared spectrum is obtained by multiplying a
// response matrix and projecting onto the reconstructed-energy axis.
//
// The reconstructed-energy distribution is
//   P(E') = \int dE\, R(E',E) F(E) S(E)
//           [1 - \sin^2 2\theta_{14}\sin^2(K_{41}/E)
//            - \sin^2 2\theta_{13}\sin^2(K_{31}/E)],
// where E' is the reconstructed energy, R(E',E) the response matrix,
// F(E) the neutrino flux and S(E) the IBD cross section.
//
// This integral is numerically expensive. The class exploits the fact that the
// response matrix is provided as a finite histogram and the flux and cross
// section are piecewise linear, allowing analytic integrals over each bin using
// sine and cosine integral functions.

#include "RooAbsCategory.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooCategoryProxy.h"
#include "RooRealProxy.h"

#include "NuOscIBDPdf.h"

#include "TH2.h"
#include "TMatrixD.h"
#include <vector>

class BinnedNuOscIBDPdf : public NuOscIBDPdf {
public:
  BinnedNuOscIBDPdf() = default;
  BinnedNuOscIBDPdf(const char *name, const char *title, RooAbsReal &x, RooAbsReal &xInt,
                     RooAbsReal &l, RooAbsReal &sin13, RooAbsReal &dm31, RooAbsReal &sin14,
                     RooAbsReal &dm41, const RooArgList &elemFracs,
                     const std::vector<std::vector<double>> &elemSpectsX,
                     const std::vector<std::vector<double>> &elemSpectsY,
                     const std::vector<double> &ibdXsecX, const std::vector<double> &ibdXsecY,
                     const TH1* hBaseline, const TH2 *hResp);
  BinnedNuOscIBDPdf(const char *name, const char *title, RooAbsReal &x, RooAbsReal &xInt,
                     RooAbsReal &l, RooAbsReal &sin13, RooAbsReal &dm31, RooAbsReal &sin14,
                     RooAbsReal &dm41, const RooArgList &elemFracs,
                     const std::vector<const TGraph *> elemSpects, const TGraph *grpXsec,
                     const TH1* hBaseline, const TH2 *hResp);
  BinnedNuOscIBDPdf(const BinnedNuOscIBDPdf &other, const char *name = 0);
  virtual TObject *clone(const char *newname) const override {
    return new BinnedNuOscIBDPdf(*this, newname);
  }
  inline virtual ~BinnedNuOscIBDPdf() override = default;

protected:
  RooRealProxy xr_;                   //!< Reconstructed energy variable
  TMatrixD respMat_;                  //!< Normalised response matrix
  std::vector<double> binsT_, binsR_; //!< Bin edges in true and reconstructed energy

  std::vector<double> ls_, lws_;      //!< Bin centers and weights of baseline

  double evaluate() const override;
  int getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                            const char *rangeName = 0) const override;
  double analyticalIntegral(int code, const char *rangeName = 0) const override;

private:
  // Cache of the full reconstructed-energy spectrum for the current parameter
  // values.  fillCache() recomputes it whenever the oscillation parameters
  // change; evaluate() then returns the appropriate cached bin value.
  void fillCache() const;

  mutable std::vector<double> cache_;
  mutable double cachedSin13_ = -1.0, cachedDm31_ = -1.0;
  mutable double cachedSin14_ = -1.0, cachedDm41_ = -1.0;

  ClassDef(BinnedNuOscIBDPdf, 1) // RooFit class definition
};

#endif
