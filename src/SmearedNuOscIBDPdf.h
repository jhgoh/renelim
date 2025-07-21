#ifndef SmearedNuOscIBDPdf_H
#define SmearedNuOscIBDPdf_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include "NuOscIBDPdf.h"

#include "TH2.h"
#include "TMatrixD.h"
#include <vector>

class SmearedNuOscIBDPdf : public NuOscIBDPdf {
public:
  SmearedNuOscIBDPdf() = default;
  SmearedNuOscIBDPdf(const char *name, const char *title, RooAbsReal& x, RooAbsReal& xInt,
      RooAbsReal& l, // baseline L (meter)
      RooAbsReal& sin13, RooAbsReal& dm31, // sin^2(2theta_13), Delta m_13^2
      RooAbsReal& sim14, RooAbsReal& dm41, // sin^2(2theta_14), Delta m_14^2
      const RooArgList& elemFracs, const std::vector<const TGraph*> elemSpects, // Fuel composition and spectrum
      const TGraph* grpXsec, // IBD Cross-section curve
      const TH2* hResp // Response matrix, (x,y) = (ETrue, EReco)
  );
  SmearedNuOscIBDPdf(const char *name, const char *title, RooAbsReal& x, RooAbsReal& xInt,
      RooAbsReal& l, // baseline L (meter)
      RooAbsReal& sin13, RooAbsReal& dm31, // sin^2(2theta_13), Delta m_13^2
      RooAbsReal& sim14, RooAbsReal& dm41, // sin^2(2theta_14), Delta m_14^2
      const RooArgList& elemFracs, const std::vector<const TGraph*> elemSpects, // Fuel composition and spectrum
      const TGraph* grpXsec, // IBD Cross-section curve
      const RooArgList& respPars // Resolution parameters
  );
  SmearedNuOscIBDPdf(const SmearedNuOscIBDPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const override { return new SmearedNuOscIBDPdf(*this,newname); }
  inline virtual ~SmearedNuOscIBDPdf() override = default;

protected:
  RooRealProxy xr_;
  RooArgList respPars_;
  TMatrixD respMat_;
  std::vector<double> binsT_, binsR_; // bin edges, trueE and recoE

  double evaluate() const override;
  int getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const override;
  double analyticalIntegral(int code, const char* rangeName = 0) const override;

private:
  ClassDef(SmearedNuOscIBDPdf,1) // Your description goes here...
};

#endif
