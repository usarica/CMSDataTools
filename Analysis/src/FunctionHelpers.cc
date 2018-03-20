#include "FunctionHelpers.h"
#include "TMath.h"


using namespace std;
using namespace HelperFunctions;


FunctionHelpers::SimpleGaussian::SimpleGaussian(double mean_, double sigma_, SimpleGaussian::RangeSetting rangeset_, double min_, double max_) :
  mean(mean_), sigma(sigma_), rangeset(rangeset_), min(min_), max(max_)
{}

double FunctionHelpers::SimpleGaussian::eval(double x){ return TMath::Gaus(x, mean, sigma, false); }
double FunctionHelpers::SimpleGaussian::norm(){
  double sigmasq=pow(sigma, 2);
  double norm=sqrt(2.*TMath::Pi())*sigma;
  if (rangeset==kHasInfRange) return norm;
  else if (rangeset==kHasLowHighRange){
    return norm*0.5*(TMath::Erfc(-(max-mean)/sqrt(2.*sigmasq))-TMath::Erfc(-(min-mean)/sqrt(2.*sigmasq)));
  }
  else if (rangeset==kHasLowRange){
    return norm*(1.-0.5*TMath::Erfc(-(min-mean)/sqrt(2.*sigmasq)));
  }
  else if (rangeset==kHasHighRange){
    return norm*(0.5*TMath::Erfc(-(max-mean)/sqrt(2.*sigmasq)));
  }
  else return 1;
}
double FunctionHelpers::SimpleGaussian::integral(double xmin, double xmax){
  if (rangeset==kHasLowHighRange || rangeset==kHasLowRange) xmin=std::max(xmin, min);
  if (rangeset==kHasLowHighRange || rangeset==kHasHighRange) xmax=std::min(xmax, max);
  double sigmasq=pow(sigma, 2);
  double norm=sqrt(2.*TMath::Pi())*sigma;
  return 0.5*norm*(TMath::Erfc(-(xmax-mean)/sqrt(2.*sigmasq))-TMath::Erfc(-(xmin-mean)/sqrt(2.*sigmasq)));
}
double FunctionHelpers::SimpleGaussian::evalNorm(double x){
  double val = eval(x);
  double n=norm();
  double res=0;
  if (n>0.) res=val/n;
  return res;
}
double FunctionHelpers::SimpleGaussian::integralNorm(double xmin, double xmax){
  double val = integral(xmin, xmax);
  double n=norm();
  double res=0;
  if (n>0.) res=val/n;
  return res;
}

void FunctionHelpers::SimpleGaussian::setRange(RangeSetting rangeset_, double min_, double max_){
  rangeset=rangeset_;
  min=min_;
  max=max_;
}
void FunctionHelpers::SimpleGaussian::setMean(double mean_){ mean=mean_; }
void FunctionHelpers::SimpleGaussian::setSigma(double sigma_){ sigma=sigma_; }
