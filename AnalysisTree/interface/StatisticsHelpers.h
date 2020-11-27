#ifndef STATISTICSHELPERS_H
#define STATISTICSHELPERS_H


namespace StatisticsHelpers{
  extern const double VAL_CL_1SIGMA;
  extern const double VAL_CL_2SIGMA;
  constexpr double VAL_CL_95PERCENT = 0.95;

  // chisq = 1, 3.84, 4 etc. for ndof=1
  double chisq_quantile(double CL, double ndof);

  // Actual range of counts. Subtract count to find difference from the central value.
  void getPoissonCountingConfidenceInterval_Frequentist(double count, double CL, double& vlow, double& vhigh);

  void getPoissonEfficiencyConfidenceInterval_Frequentist(double sw_total, double sw_passed, double swsq_total, double CL, double& vlow, double& vhigh);
  void getPoissonEfficiencyConfidenceInterval_Bayesian(double sw_total, double sw_passed, double swsq_total, double CL, double alpha, double beta, double& vlow, double& vhigh);

}

#endif
