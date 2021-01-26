#include <cmath>
#include "FunctionHelpers.h"
#include "ExtendedBinning.h"
#include "HelperFunctions.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace HelperFunctions;
using namespace MELAStreamHelpers;


FunctionHelpers::SimpleGaussian::SimpleGaussian(double mean_, double sigma_, SimpleGaussian::RangeSetting rangeset_, double min_, double max_) :
  mean(mean_), sigma(sigma_), rangeset(rangeset_), min(min_), max(max_)
{}

double FunctionHelpers::SimpleGaussian::eval(double x){
  if (sigma==0.){
    if (x==mean) return 1;
    else return 0;
  }
  else return TMath::Gaus(x, mean, sigma, false);
}
double FunctionHelpers::SimpleGaussian::norm(){
  if (sigma==0.){
    if (
      rangeset==kHasInfRange
      ||
      (rangeset==kHasLowHighRange && min<=mean && mean<=max)
      ||
      (rangeset==kHasLowRange && min<mean)
      ||
      (rangeset==kHasHighRange && mean<=max)
      ) return 1;
    else return 0;
  }
  else{
    double sigmasq = std::pow(sigma, 2);
    double norm = std::sqrt(2.*TMath::Pi())*sigma;
    if (rangeset==kHasInfRange) return norm;
    else if (rangeset==kHasLowHighRange){
      return norm*0.5*(TMath::Erfc(-(max-mean)/std::sqrt(2.*sigmasq))-TMath::Erfc(-(min-mean)/std::sqrt(2.*sigmasq)));
    }
    else if (rangeset==kHasLowRange){
      return norm*(1.-0.5*TMath::Erfc(-(min-mean)/std::sqrt(2.*sigmasq)));
    }
    else if (rangeset==kHasHighRange){
      return norm*(0.5*TMath::Erfc(-(max-mean)/std::sqrt(2.*sigmasq)));
    }
    else return 1;
  }
}
double FunctionHelpers::SimpleGaussian::integral(double xmin, double xmax){
  if (rangeset==kHasLowHighRange || rangeset==kHasLowRange) xmin=std::max(xmin, min);
  if (rangeset==kHasLowHighRange || rangeset==kHasHighRange) xmax=std::min(xmax, max);
  if (sigma==0.){
    if (
      rangeset==kHasInfRange
      ||
      (rangeset==kHasLowHighRange && xmin<=mean && mean<=xmax)
      ||
      (rangeset==kHasLowRange && xmin<mean)
      ||
      (rangeset==kHasHighRange && mean<=xmax)
      ) return 1;
    else return 0;
  }
  else{
    double sigmasq = std::pow(sigma, 2);
    double norm = std::sqrt(2.*TMath::Pi())*sigma;
    return 0.5*norm*(TMath::Erfc(-(xmax-mean)/std::sqrt(2.*sigmasq))-TMath::Erfc(-(xmin-mean)/std::sqrt(2.*sigmasq)));
  }
}
double FunctionHelpers::SimpleGaussian::evalNorm(double x){
  double val = eval(x);
  double n = norm();
  double res = 0;
  if (n>0.) res = val/n;
  return res;
}
double FunctionHelpers::SimpleGaussian::integralNorm(double xmin, double xmax){
  double val = integral(xmin, xmax);
  double n = norm();
  double res = 0;
  if (n>0.) res = val/n;
  return res;
}

void FunctionHelpers::SimpleGaussian::setRange(RangeSetting rangeset_, double min_, double max_){
  rangeset=rangeset_;
  min=min_;
  max=max_;
}
void FunctionHelpers::SimpleGaussian::setMean(double mean_){ mean=mean_; }
void FunctionHelpers::SimpleGaussian::setSigma(double sigma_){ sigma=sigma_; }


FunctionHelpers::PiecewisePolynomial::PiecewisePolynomial(const int nfcn_, const int polyndof_) :
  nfcn(nfcn_), polyndof(polyndof_),
  nnodes(nfcn-1), // How many nodes are in between
  ndof_endfcn(polyndof-1), // 1 degree for the function value
  ndof_middlefcn(polyndof-2) // +1 for slope of the other node
{
  assert((nfcn>2 && polyndof>=2) || (nfcn==2 && polyndof>=1));
}
FunctionHelpers::PiecewisePolynomial::PiecewisePolynomial(FunctionHelpers::PiecewisePolynomial const& other) :
  nfcn(other.nfcn), polyndof(other.polyndof),
  nnodes(other.nnodes),
  ndof_endfcn(other.ndof_endfcn),
  ndof_middlefcn(other.ndof_middlefcn),
  par(other.par)
{}

void FunctionHelpers::PiecewisePolynomial::setParameters(std::vector<double> par_){ par=par_; }

double FunctionHelpers::PiecewisePolynomial::eval(double x){
  // If we say the form of the polynomial is [0] + [1]*x + [2]*x2 + [3]*x3 + [4]*x4...,
  // use the highest two orders for matching at the nodes and free the rest.
  const double d_epsilon = 1e-14;
  if ((int) par.size()!=2*ndof_endfcn+(nfcn-2)*ndof_middlefcn+nnodes) assert(0);

  int npars_reduced[nfcn];
  for (int index=0; index<nfcn; index++){
    if (index==0 || index==nfcn-1) npars_reduced[index] = ndof_endfcn;
    else npars_reduced[index] = ndof_middlefcn;
  }

  std::vector<double> node(nnodes, 0); // First [0,...,nnodes-1] parameters are nodes
  std::vector<std::vector<double>> pars_full(nfcn, std::vector<double>(polyndof, 0)); // The full coefficients array

  for (int ip=0; ip<nnodes; ip++) node[ip] = par[ip];
  // Check for the nodes to be consecutive
  for (int ip=0; ip<nnodes; ip++){
    for (int ip2=ip+1; ip2<nnodes; ip2++){
      if (node[ip]>node[ip2]) return d_epsilon;
    }
  }
  int pos_ctr = nnodes;
  for (int index=0; index<nfcn; index++){
    for (int ipar=0; ipar<npars_reduced[index]; ipar++){
      if (!(index==(nfcn-1) && ipar==(npars_reduced[index]-1))) pars_full[index][ipar] = par[pos_ctr];
      else pars_full[index][ipar+1] = par[pos_ctr]; // Special case to avoid singular matrix. This corresponds to having the x^n contribution free instead of x^(n-1)
      pos_ctr++;
    }
  }

  std::vector<std::vector<double>> xton(nnodes, std::vector<double>(polyndof, 0)); // Array of node^power
  std::vector<std::vector<double>> nxtom(nnodes, std::vector<double>(polyndof, 0)); // Array of power*node^(power-1)
  for (int inode=0; inode<nnodes; inode++){
    for (int ipow=0; ipow<polyndof; ipow++){
      if (ipow==0) xton[inode][ipow]=1; // nxtom==0
      else if (ipow==1){
        xton[inode][ipow]=node[inode];
        nxtom[inode][ipow]=1;
      }
      else{
        xton[inode][ipow]=std::pow(node[inode], ipow);
        nxtom[inode][ipow]=((double) ipow)*std::pow(node[inode], ipow-1);
      }
    }
  }

  std::vector<double> ysbar_nodes(2*nnodes, 0);
  std::vector<std::vector<double>> coeff_ysbar(2*nnodes, std::vector<double>(2*nnodes, 0));
  int cstart=-1;
  for (int inode=0; inode<nnodes; inode++){
    int i=inode;
    int j=i+1;
    double sign_i = 1, sign_j=-1;
    for (int ip=0; ip<npars_reduced[i]; ip++){
      ysbar_nodes[inode] += sign_i*pars_full[i][ip]*xton[inode][ip];
      ysbar_nodes[nnodes+inode] += sign_i*pars_full[i][ip]*nxtom[inode][ip];
    }
    for (int ip=0; ip<npars_reduced[j]; ip++){
      if (!(j==(nfcn-1) && ip==(npars_reduced[j]-1))){
        ysbar_nodes[inode] += sign_j*pars_full[j][ip]*xton[inode][ip];
        ysbar_nodes[nnodes+inode] += sign_j*pars_full[j][ip]*nxtom[inode][ip];
      }
      else{
        ysbar_nodes[inode] += sign_j*pars_full[j][ip+1]*xton[inode][ip+1];
        ysbar_nodes[nnodes+inode] += sign_j*pars_full[j][ip+1]*nxtom[inode][ip+1];
      }
    }

    if (cstart>=0){
      coeff_ysbar[inode][cstart] = -sign_i*xton[inode][polyndof-2];
      coeff_ysbar[nnodes + inode][cstart] = -sign_i*nxtom[inode][polyndof-2];
    }
    coeff_ysbar[inode][cstart+1] = -sign_i*xton[inode][polyndof-1];
    coeff_ysbar[nnodes + inode][cstart+1] = -sign_i*nxtom[inode][polyndof-1];
    coeff_ysbar[inode][cstart+2] = -sign_j*xton[inode][polyndof-2];
    coeff_ysbar[nnodes + inode][cstart+2] = -sign_j*nxtom[inode][polyndof-2];
    if ((cstart+3)<2*nnodes){
      coeff_ysbar[inode][cstart+3] = -sign_j*xton[inode][polyndof-1];
      coeff_ysbar[nnodes + inode][cstart+3] = -sign_j*nxtom[inode][polyndof-1];
    }
    cstart+=2;
  }

  TVectorD polyvec(2*nnodes, ysbar_nodes.data());
  TMatrixD polycoeff(2*nnodes, 2*nnodes);
  for (int i=0; i<2*nnodes; i++){ for (int j=0; j<2*nnodes; j++) polycoeff(i, j)=coeff_ysbar[i][j]; }
  double testdet=0;
  TMatrixD polycoeff_inv = polycoeff.Invert(&testdet);
  if (testdet!=0){
    TVectorD unknowncoeffs = polycoeff_inv*polyvec;
    pos_ctr=0;
    for (int index=0; index<nfcn; index++){
      for (int ip=npars_reduced[index]; ip<polyndof; ip++){
        if (!(index==(nfcn-1) && ip==npars_reduced[index])) pars_full[index][ip] = unknowncoeffs[pos_ctr];
        else pars_full[index][ip-1] = unknowncoeffs[pos_ctr];
        pos_ctr++;
      }
    }

    int index_chosen=0;
    for (int index=0; index<nnodes-1; index++){
      if (x>=node[index] && x<node[index+1]){
        index_chosen = index+1;
        break;
      }
    }
    if (x>=node[nnodes-1]) index_chosen = nfcn-1;

    double res = 0;
    for (int ip=0; ip<polyndof; ip++) res += pars_full[index_chosen][ip]*std::pow(x, ip);
    return res;
  }
  else{
    MELAerr << "PiecewisePolynomial::eval: Something went wrong, and the determinant is 0!" << endl;
    return d_epsilon;
  }
}
