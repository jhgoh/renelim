#ifndef SmearedNuOscIBDPdf_H
#define SmearedNuOscIBDPdf_H

// SmearedNuOscIBDPdf: energy spectrum of sterile neutrinos including detector resolution.
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

class SmearedNuOscIBDPdf : public NuOscIBDPdf {
public:
  SmearedNuOscIBDPdf() = default;
  SmearedNuOscIBDPdf(const char *name, const char *title, RooAbsReal &x, RooAbsReal &xInt,
                     RooAbsReal &l, RooAbsReal &sin13, RooAbsReal &dm31, RooAbsReal &sin14,
                     RooAbsReal &dm41, const RooArgList &elemFracs,
                     const std::vector<const TGraph *> elemSpects, const TGraph *grpXsec,
                     const TH2 *hResp);
  SmearedNuOscIBDPdf(const SmearedNuOscIBDPdf &other, const char *name = 0);
  virtual TObject *clone(const char *newname) const override {
    return new SmearedNuOscIBDPdf(*this, newname);
  }
  inline virtual ~SmearedNuOscIBDPdf() override = default;

protected:
  RooRealProxy xr_;                   //!< Reconstructed energy variable
  TMatrixD respMat_;                  //!< Normalised response matrix
  std::vector<double> binsT_, binsR_; //!< Bin edges in true and reconstructed energy

  double evaluate() const override;
  int getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                            const char *rangeName = 0) const override;
  double analyticalIntegral(int code, const char *rangeName = 0) const override;

private:
  ClassDef(SmearedNuOscIBDPdf, 1) // RooFit class definition
};

#endif
