#ifndef STATISTICSHELPERS_H
#define STATISTICSHELPERS_H


namespace StatisticsHelpers{
  extern const double VAL_CL_1SIGMA;
  extern const double VAL_CL_2SIGMA;
  extern const double VAL_CL_3SIGMA;
  extern const double VAL_CL_4SIGMA;
  extern const double VAL_CL_5SIGMA;
  constexpr double VAL_CL_95PERCENT = 0.95;

  // CL = 0.68, 0.95 etc. for ndof=1
  double getConfidenceLevelValue(double stddev, double ndof);

  // chisq = 1, 3.84, 4 etc. for ndof=1
  double chisq_quantile(double CL, double ndof);

  // Actual range of counts. Subtract count to find difference from the central value.
  void getPoissonCountingConfidenceInterval_Frequentist(double sw_total, double swsq_total, double CL, double& vlow, double& vhigh);
  void getPoissonCountingConfidenceInterval_Frequentist(double count, double CL, double& vlow, double& vhigh);

  void getPoissonEfficiencyConfidenceInterval_Frequentist(double sw_total, double sw_passed, double swsq_total, double CL, double& vlow, double& vhigh);
  void getPoissonEfficiencyConfidenceInterval_Bayesian(double sw_total, double sw_passed, double swsq_total, double CL, double alpha, double beta, double& vlow, double& vhigh);

}

#endif
