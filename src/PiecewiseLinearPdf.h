#ifndef PiecewiseLinearPdf_H
#define PiecewiseLinearPdf_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include <vector>
 
class PiecewiseLinearPdf : public RooAbsPdf {
public:
  PiecewiseLinearPdf() = default;
  PiecewiseLinearPdf(const char *name, const char *title, RooAbsReal& x,
                     const int n, const double* xx, const double* yy);
  PiecewiseLinearPdf(const PiecewiseLinearPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new PiecewiseLinearPdf(*this,newname); }
  inline virtual ~PiecewiseLinearPdf() = default; 

protected:
  RooRealProxy x_;
  std::vector<double> xx_, yy_;
  
  void set(const int n, const double* xx, const double* yy);
  double evaluate() const;

private:
  ClassDef(PiecewiseLinearPdf,1) // Your description goes here...
};
 
#endif
