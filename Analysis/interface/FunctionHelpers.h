#ifndef FUNCTIONHELPERS_H
#define FUNCTIONHELPERS_H

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "ExtendedBinning.h"
#include "HelperFunctions.h"


namespace FunctionHelpers{

  class SimpleGaussian{
  public:
    enum RangeSetting{
      kHasLowHighRange,
      kHasLowRange,
      kHasHighRange,
      kHasInfRange
    };

  protected:
    double mean;
    double sigma;

    RangeSetting rangeset;
    double min;
    double max;

  public:

    SimpleGaussian(double mean_, double sigma_, RangeSetting rangeset_, double min_, double max_);

    double eval(double x);
    double integral(double xmin, double xmax);
    double norm();
    double evalNorm(double x);
    double integralNorm(double xmin, double xmax);

    void setRange(RangeSetting rangeset_, double min_, double max_);
    void setMean(double mean_);
    void setSigma(double sigma_);

  };

}

#endif
