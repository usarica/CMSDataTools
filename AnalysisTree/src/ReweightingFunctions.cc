#include <cassert>
#include "ReweightingFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


std::vector<float*> ReweightingFunctions::getWeightRefs(BaseTree* tree, std::vector<TString> const& strWeights){
  std::vector<float*> res;
  if (!tree || strWeights.empty()) return res;
  res.reserve(strWeights.size());
  for (auto const& s:strWeights){
    float* v = nullptr;
    tree->getValRef<float>(s, v);
    if (!v) MELAerr << "ReweightingFunctions::getWeightRefs: Could not get the reference for weight " << s << endl;
    res.push_back(v);
  }
  return res;
}
float* ReweightingFunctions::getWeightRef(BaseTree* tree, TString const& strWeight){
  float* res=nullptr;
  if (!tree || strWeight=="") return res;
  tree->getValRef<float>(strWeight, res);
  if (!res) MELAerr << "ReweightingFunctions::getWeightRef: Could not get the reference for weight " << strWeight << endl;
  return res;
}

float ReweightingFunctions::getSimpleWeight(BaseTree* tree, std::vector<float*> const& vals){
  float res=0;
  if (!tree) return res;
  res=1;
  for (auto const& v:vals){ if (v) res *= *v; }
  return res;
}
float ReweightingFunctions::getA1PlusB1Weight(BaseTree* tree, std::vector<float*> const& vals){
  float res=0;
  if (!tree) return res;
  assert(vals.size()==2 && vals.front() && vals.back());
  res = *(vals.front()) + *(vals.back());
  return res;
}
float ReweightingFunctions::getA1MinusB1Weight(BaseTree* tree, std::vector<float*> const& vals){
  float res=0;
  if (!tree) return res;
  assert(vals.size()==2 && vals.front() && vals.back());
  res = *(vals.front()) - *(vals.back());
  return res;
}
float ReweightingFunctions::getA1MinusB1MinusC1Weight(BaseTree* tree, std::vector<float*> const& vals){
  float res=0;
  if (!tree) return res;
  assert(vals.size()==3 && vals.front() && vals.at(1) && vals.back());
  res = *(vals.front()) - *(vals.at(1)) - *(vals.back());
  return res;
}
float ReweightingFunctions::getOnePlusB1OverA1Weight(BaseTree* tree, std::vector<float*> const& vals){
  float res=0;
  if (!tree) return res;
  assert(vals.size()==2 && vals.front() && vals.back());
  res = 1. + *(vals.back()) / *(vals.front());
  return res;
}
float ReweightingFunctions::getOneMinusB1OverA1Weight(BaseTree* tree, std::vector<float*> const& vals){
  float res=0;
  if (!tree) return res;
  assert(vals.size()==2 && vals.front() && vals.back());
  res = 1. - *(vals.back()) / *(vals.front());
  return res;
}
float ReweightingFunctions::getA1OverB1Weight(BaseTree* tree, std::vector<float*> const& vals){
  float res=0;
  if (!tree) return res;
  res=1;
  assert(vals.size()==2 && vals.front() && vals.back());
  res = (*(vals.front()))/(*(vals.back()));
  return res;
}

float ReweightingFunctions::getAbsWeightThresholdByNeff(BaseTree* tree, std::vector<float*> const& vals, ReweightingFunction_t rule, double thr_Neff, TVar::VerbosityLevel verbosity){
  float res = -1;

  int nEntries = tree->getNEvents();
  thr_Neff = (thr_Neff>0. ? std::min(thr_Neff, double(nEntries)/3.*2.) : static_cast<double>(nEntries));
  unsigned int npos = 0;
  double Neff = 0;
  double sum_wgts[2]={ 0 }; // [0]: w, [1]: w^2
  std::vector<float> smallest_weights;
  if (verbosity>=TVar::ERROR) MELAout << "ReweightingFunctions::getAbsWeightThresholdByNeff: Determining the weight thresholds (Neff threshold = " << thr_Neff << ", number of events = " << nEntries << ")..." << endl;
  for (int ev=0; ev<nEntries; ev++){
    tree->getEvent(ev);
    HelperFunctions::progressbar(ev, nEntries);

    float wgt_combined = std::abs(rule(tree, vals));

    HelperFunctions::addByLowest(smallest_weights, wgt_combined, false);
    if (ev%100000 == 0 || ev == nEntries-1){
      sum_wgts[0] = sum_wgts[1] = 0;
      npos = 0;
      for (auto const& wgt:smallest_weights){
        sum_wgts[0] += wgt;
        sum_wgts[1] += wgt*wgt;
        Neff = std::pow(sum_wgts[0], 2) / sum_wgts[1];
        npos++;
        if (Neff>=thr_Neff) break;
      }
      if (verbosity>=TVar::ERROR) MELAout << "\t- Current Neff = " << Neff << " over " << ev+1 << " events..." << endl;
    }
    if (Neff>=thr_Neff) break;
  }

  if (!smallest_weights.empty()){
    res = (sum_wgts[0] + std::sqrt(sum_wgts[1]*Neff)) / (Neff-1.);
    if (verbosity>=TVar::ERROR) MELAout
      << "\t- " << res
      << " is the default weight threshold calculated from sN=" << sum_wgts[0] << ", vN=" << sum_wgts[1] << ", nN=" << Neff
      << " (N=" << npos << " / " << smallest_weights.size() << ", wN=" << smallest_weights.at(npos-1) << ", wLast=" << smallest_weights.back() << ")."
      << endl;
  }
  else{
    if (verbosity>=TVar::INFO) MELAout << "\t- No weight threshold is found." << endl;
  }

  return res;
}

float ReweightingFunctions::getSimpleVariable(BaseTree* tree, std::vector<float*> const& vals){
  assert(tree && vals.size()==1 && vals.front());
  return *(vals.front());
}
int ReweightingFunctions::getSimpleVariableBin(BaseTree* tree, ExtendedBinning const& binning, std::vector<float*> const& vals){
  return binning.getBin(getSimpleVariable(tree, vals));
}


std::vector<double> ReweightingFunctions::getSimpleNeffThrsPerBin(
  BaseTree* tree,
  ExtendedBinning const& binning, std::vector<float*> const& var_vals, ReweightingVariableBinFunction_t varbin_rule,
  std::vector<float*> const& wgt_vals, ReweightingFunction_t wgt_rule,
  double thr_Neff,
  TVar::VerbosityLevel verbosity
){
  unsigned int const nbins = (!binning.isValid() ? static_cast<unsigned int>(1) : binning.getNbins());

  std::vector<double> res(nbins, 0);

  int nEntries = tree->getNEvents();
  if (nEntries==0) return res;
  thr_Neff = (thr_Neff>0. ? std::min(thr_Neff, double(nEntries)/3.*2.) : static_cast<double>(nEntries));
  if (verbosity>=TVar::ERROR) MELAout << "ReweightingFunctions::getSimpleNeffThrsPerBin: Determining the distribution of Neff over the bins..." << endl;

  std::vector<unsigned int> nEntries_per_bin(nbins, 0);
  for (int ev=0; ev<nEntries; ev++){
    tree->getEvent(ev);
    HelperFunctions::progressbar(ev, nEntries);

    float wgt_combined = std::abs(wgt_rule(tree, wgt_vals));
    if (wgt_combined==0.f) continue;

    int ibin = varbin_rule(tree, binning, var_vals);
    if (ibin<0) ibin=0;
    else if (ibin>=(int) nbins) ibin = nbins-1;
    nEntries_per_bin.at(ibin)++;
  }

  // Determine the actual Neff thresholds per bin based on the raw event distribution
  for (unsigned int ibin=0; ibin<nbins; ibin++) res.at(ibin) = thr_Neff * ((double) nEntries_per_bin.at(ibin)) / ((double) nEntries);
  if (verbosity>=TVar::ERROR) MELAout << "\t- Neff thresholds = " << res << endl;

  return res;
}
std::vector<float> ReweightingFunctions::getAbsWeightThresholdsPerBinByNeff(
  BaseTree* tree,
  std::vector<float*> const& wgt_vals, ReweightingFunction_t wgt_rule,
  ExtendedBinning const& binning, std::vector<float*> const& var_vals, ReweightingVariableBinFunction_t varbin_rule,
  double tolerance,
  TVar::VerbosityLevel verbosity
){
  if (tolerance>1.) tolerance = 1.;
  unsigned int const nbins = (!binning.isValid() ? static_cast<unsigned int>(1) : binning.getNbins());

  std::vector<float> res(nbins, -1);

  int nEntries = tree->getNEvents();
  if (verbosity>=TVar::ERROR) MELAout << "ReweightingFunctions::getAbsWeightThresholdsPerBinByNeff: Determining the weight thresholds (number of events = " << nEntries << ")..." << endl;
  if (nEntries==0) return res;

  // A bit brute-force, but best to do this way...
  {
    std::vector< std::vector<float> > weights_all(nbins, std::vector<float>());
    for (int ev=0; ev<nEntries; ev++){
      tree->getEvent(ev);
      HelperFunctions::progressbar(ev, nEntries);

      float wgt_combined = std::abs(wgt_rule(tree, wgt_vals));
      if (wgt_combined==0.f) continue;
      int ibin = varbin_rule(tree, binning, var_vals);
      if (ibin<0) ibin=0;
      else if (ibin>=(int) nbins) ibin = nbins-1;
      HelperFunctions::addByLowest(weights_all.at(ibin), wgt_combined, false);
    }
    for (unsigned int ibin=0; ibin<nbins; ibin++){
      auto const& weights_bin = weights_all.at(ibin);
      unsigned int const nevts_bin = weights_bin.size();
      if (nevts_bin<=1) continue;
      unsigned int nevts_Neff = (tolerance>0. ? (static_cast<double>(nevts_bin)*tolerance + 0.5) : (static_cast<double>(nevts_bin)-std::sqrt(static_cast<double>(nevts_bin))+0.5));
      double last_smallest_weight = -1;
      double Neff = 0;
      std::pair<double, double> sumW(0, 0);
      unsigned int ipos=0;
      for (auto const& wgt:weights_bin){
        auto sumW_new = sumW;
        sumW_new.first += wgt;
        sumW_new.second += std::pow(wgt, 2);
        double Neff_new = std::pow(sumW_new.first, 2)/sumW_new.second;
        if (ipos>=nevts_Neff-1 && (Neff_new<Neff || (last_smallest_weight>0. && wgt>10.*last_smallest_weight))) break;
        Neff = Neff_new;
        sumW = sumW_new;
        last_smallest_weight = wgt;
        ipos++;
      }
      if (ipos==nevts_bin){
        if (verbosity>=TVar::ERROR) MELAout << "\t\t- Threshold for bin " << ibin << " is kept at " << res.at(ibin) << endl;
        continue;
      }
      auto const& last_good_wgt = weights_bin.at(ipos-1);
      auto const& first_bad_wgt = weights_bin.at(ipos);
      res.at(ibin) = (last_good_wgt + first_bad_wgt) / 2.;
      MELAout
        << "\t\t- " << res.at(ibin)
        << " is the weight threshold found from iterating over the sample. The calculated threshold would be " << (Neff>1. ? (sumW.first + std::sqrt(sumW.second*Neff)) / (Neff-1.) : -1.)
        << " from sN=" << sumW.first << ", vN=" << sumW.second << ", nN=" << Neff
        << " (N=" << ipos << " / " << weights_bin.size() << ", wN=" << weights_bin.at(ipos-1) << ", wLast=" << weights_bin.back()
        << "). Expected fraction to accept events: " << ((double) ipos) / ((double) weights_bin.size())
        << endl;
    }
  }

  return res;
}
std::vector<float> ReweightingFunctions::getAbsWeightThresholdsPerBinByFixedFractionalThreshold(
  BaseTree* tree,
  std::vector<float*> const& wgt_vals, ReweightingFunction_t wgt_rule,
  ExtendedBinning const& binning, std::vector<float*> const& var_vals, ReweightingVariableBinFunction_t varbin_rule,
  double frac, double tolerance,
  TVar::VerbosityLevel verbosity
){
  if (frac<=0.) frac = 0.9999;
  if (tolerance<=0.) tolerance = 5.;
  unsigned int const nbins = (!binning.isValid() ? static_cast<unsigned int>(1) : binning.getNbins());

  std::vector<float> res(nbins, -1);

  int nEntries = tree->getNEvents();
  if (verbosity>=TVar::ERROR) MELAout << "ReweightingFunctions::getAbsWeightThresholdsPerBinByNeff: Determining the weight thresholds (number of events = " << nEntries << ")..." << endl;
  if (nEntries==0) return res;

  // A bit brute-force, but best to do this way...
  {
    std::vector< std::vector<float> > weights_all(nbins, std::vector<float>());
    for (int ev=0; ev<nEntries; ev++){
      tree->getEvent(ev);
      HelperFunctions::progressbar(ev, nEntries);

      float wgt_combined = std::abs(wgt_rule(tree, wgt_vals));
      if (wgt_combined==0.f) continue;
      int ibin = varbin_rule(tree, binning, var_vals);
      if (ibin<0) ibin=0;
      else if (ibin>=(int) nbins) ibin = nbins-1;

      HelperFunctions::addByLowest(weights_all.at(ibin), wgt_combined, false);
    }
    for (unsigned int ibin=0; ibin<nbins; ibin++){
      if (verbosity>=TVar::ERROR) MELAout << "\t- Bin " << ibin << ":" << endl;
      auto const& weights_bin = weights_all.at(ibin);
      unsigned int const nevts_bin = weights_bin.size();

      if (verbosity>=TVar::ERROR) MELAout << "\t\t- Nevts = " << nevts_bin << "<3, skipping..." << endl;
      if (nevts_bin<3) continue;

      unsigned int index_entry = std::max((unsigned int) std::ceil(float(nevts_bin)*frac), (unsigned int) 1);
      unsigned int index_entry_prev=index_entry-1;
      float threshold = (weights_bin.at(index_entry_prev) + weights_bin.at(index_entry))*0.5;
      if (verbosity>=TVar::ERROR) MELAout << "\t\t- Raw threshold before checking tolerance: " << threshold << endl;
      if (weights_bin.back()<threshold*tolerance) threshold = -1; // Prevent false-positives
      res.at(ibin) = threshold;
      MELAout << "\t\t- Final threshold: " << res.at(ibin) << endl;
    }
  }

  return res;
}
