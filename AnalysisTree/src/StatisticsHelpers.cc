#include "Math/DistFunc.h"
#include "TEfficiency.h"
#include "StatisticsHelpers.h"


namespace StatisticsHelpers{
  const double VAL_CL_1SIGMA = ROOT::Math::chisquared_cdf(1., 1., 0.);
  const double VAL_CL_2SIGMA = ROOT::Math::chisquared_cdf(4., 1., 0.);

  double chisq_quantile(double CL, double ndof){ return ROOT::Math::chisquared_quantile(CL, ndof); }

  void getPoissonCountingConfidenceInterval_Frequentist(double count, double CL, double& vlow, double& vhigh){
    const double quant = (1. - CL) / 2.;
    vlow = (count == 0. ? 0. : ROOT::Math::chisquared_quantile(quant, 2. * count) / 2.);
    vhigh = ROOT::Math::chisquared_quantile_c(quant, 2 * (count + 1.)) / 2.;
  }

  void getPoissonEfficiencyConfidenceInterval_Frequentist(double sw_total, double sw_passed, double swsq_total, double CL, double& vlow, double& vhigh){
    double normval = sw_total/swsq_total;
    double total = sw_total*normval;
    double passed = sw_passed*normval;
    vlow = TEfficiency::ClopperPearson(total, passed, CL, false);
    vhigh = TEfficiency::ClopperPearson(total, passed, CL, true);
  }
  void getPoissonEfficiencyConfidenceInterval_Bayesian(double sw_total, double sw_passed, double swsq_total, double CL, double alpha, double beta, double& vlow, double& vhigh){
    double normval = sw_total/swsq_total;
    double passed = sw_passed*normval + alpha;
    double failed = (sw_total - sw_passed)*normval + beta;
    vlow = TEfficiency::BetaCentralInterval(CL, passed, failed, false);
    vhigh = TEfficiency::BetaCentralInterval(CL, passed, failed, true);
  }

}
