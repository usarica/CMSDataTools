#ifndef FUNCTIONHELPERS_H
#define FUNCTIONHELPERS_H

#include <vector>


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

  // Piecewise polynomial that looks like
  // -- . -- ... -- . --
  class PiecewisePolynomial{
  protected:
    const int nfcn; // Number of piecewise functions
    const int polyndof; // Ndof of the polynomial (e.g. 4: cubic)
    const int nnodes; // How many nodes are in between
    const int ndof_endfcn; // Number of degrees of freedom in fcns in the middle of the nodes
    const int ndof_middlefcn; // Number of degrees of freedom in fcns outside the nodes

    // First [0,...,nnodes-1] parameters are nodes
    // There should be 2*ndof_endfcn+(nfcn-2)*ndof_middlefcn more parameters for the free dofs in the functions
    // where ndof_endfcn=polyndof-1 and ndof_middlefcn=polyndof-2
    std::vector<double> par;

  public:
    PiecewisePolynomial(const int nfcn_, const int polyndof_);
    PiecewisePolynomial(PiecewisePolynomial const& other);

    double eval(double x);

    void setParameters(std::vector<double> pars_);

  };

}

#endif
