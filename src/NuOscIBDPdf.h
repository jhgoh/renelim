#ifndef NuOscIBDPdf_H
#define NuOscIBDPdf_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include <vector>
 
class NuOscIBDPdf : public RooAbsPdf {
public:
  NuOscIBDPdf() = default;
  NuOscIBDPdf(const char *name, const char *title, RooAbsReal& x,
      RooAbsReal& l, // baseline L (meter)
      RooAbsReal& sin13, RooAbsReal& dm31, // sin^2(2theta_13), Delta m_13^2
      RooAbsReal& sim14, RooAbsReal& dm41, // sin^2(2theta_14), Delta m_14^2
      RooAbsReal& fU235, RooAbsReal& fPu239, RooAbsReal& fPu241, // Fuel composition (no U238)
      const TGraph* grpU235, const TGraph* grpPu239, const TGraph* grpPu241, // Spectrum (Huber-Mueller)
      const TGraph* grpXsec // IBD Cross-section curve
  );
  NuOscIBDPdf(const char *name, const char *title, RooAbsReal& x,
      RooAbsReal& l, // baseline L (meter)
      RooAbsReal& sin13, RooAbsReal& dm31, // sin^2(2theta_13), Delta m_13^2
      RooAbsReal& sim14, RooAbsReal& dm41, // sin^2(2theta_14), Delta m_14^2
      RooAbsReal& fU235, RooAbsReal& fU238, RooAbsReal& fPu239, RooAbsReal& fPu241, // Fuel composition
      const TGraph* grpU235, const TGraph* grpU238, const TGraph* grpPu239, const TGraph* grpPu241, // Spectrum (Huber-Mueller)
      const TGraph* grpXsec // IBD Cross-section curve
  );
  NuOscIBDPdf(const NuOscIBDPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new NuOscIBDPdf(*this,newname); }
  inline virtual ~NuOscIBDPdf() = default; 

protected:
  RooRealProxy x_; // Variable
  RooRealProxy l_, sin13_, dm31_, sin14_, dm41_; // Parameters
  RooRealProxy fU235_, fU238_, fPu239_, fPu241_;

  std::vector<double> xx_U235_, yy_U235_;
  std::vector<double> xx_U238_, yy_U238_;
  std::vector<double> xx_Pu239_, yy_Pu239_;
  std::vector<double> xx_Pu241_, yy_Pu241_;

  std::vector<double> xx_Xsec_, yy_Xsec_;

  bool hasU238_;
  
  double evaluate() const;
  int getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const override;
  double analyticalIntegral(int code, const char* rangeName = 0) const override;

  void loadFromTGraph(const TGraph* grp, std::vector<double>& xx, std::vector<double>& yy);
  double interpolate(const double x, const std::vector<double>& xx, const std::vector<double>& yy) const;

private:
  ClassDef(NuOscIBDPdf,1) // Your description goes here...
};
 
#endif
