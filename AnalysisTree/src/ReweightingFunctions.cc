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
std::vector<std::vector<float>> ReweightingFunctions::getAbsWeightThresholdsPerBinByFixedFractionalThreshold(
  BaseTree* tree,
  std::vector<std::vector<float*>> const& wgt_vals_list, std::vector<ReweightingFunction_t> const& wgt_rule_list, std::vector< std::pair<double, double> > const& frac_tolerance_pair_list,
  ExtendedBinning const& binning, std::vector<float*> const& var_vals, ReweightingVariableBinFunction_t varbin_rule,
  TVar::VerbosityLevel verbosity
){
  constexpr unsigned int nevts_bin_skipThr = 3;

  unsigned int const nbins = (!binning.isValid() ? static_cast<unsigned int>(1) : binning.getNbins());
  unsigned int const nhypos = wgt_vals_list.size();
  assert(wgt_rule_list.size() == nhypos);
  assert(frac_tolerance_pair_list.size() == nhypos);

  std::vector<std::vector<float>> res(nhypos, std::vector<float>(nbins, -1));
  if (nhypos==0){
    if (verbosity>=TVar::ERROR) MELAerr << "ReweightingFunctions::getAbsWeightThresholdsPerBinByFixedFractionalThreshold: No hypotheses are passed." << endl;
    return res;
  }

  int nEntries = tree->getNEvents();
  if (verbosity>=TVar::ERROR) MELAout << "ReweightingFunctions::getAbsWeightThresholdsPerBinByFixedFractionalThreshold: Determining the weight thresholds (number of events = " << nEntries << ")..." << endl;
  if (nEntries==0) return res;

  std::vector<double> frac_list; frac_list.reserve(nhypos);
  std::vector<double> frac_rem_list; frac_rem_list.reserve(nhypos);
  std::vector<double> tolerance_list; tolerance_list.reserve(nhypos);
  std::vector<unsigned int> maxPrunedSize_list; maxPrunedSize_list.reserve(nhypos);
  std::vector<bool> skip_list(nhypos, false);
  // Vectors to accumulate the highest 'maxPrunedSize' number of weights
  std::vector< std::vector<unsigned int> > counts_all(nhypos, std::vector<unsigned int>(nbins, 0));
  std::vector< std::vector<double> > sumWgts_all(nhypos, std::vector<double>(nbins, 0));
  std::vector< std::vector< std::vector<float> > > weights_all(nhypos, std::vector< std::vector<float> >(nbins, std::vector<float>()));
  {
    unsigned int ihypo=0;
    for (auto const& pp:frac_tolerance_pair_list){
      double frac = pp.first;
      double tolerance = pp.second;
      assert(frac<=1. && (tolerance>=1. || tolerance<=0.));
      if (frac<0.) frac = 0.9999;
      if (tolerance<=0.) tolerance = 5.;
      // A bit brute-force, but best to do this way...
      double const frac_rem = 1. - frac;
      // Let the requested fraction threshold per bin to be f, number of total events to be N, and the number of events in bin i to be N_i such that sum_i(N_i)=N.
      // We don't know N_i before looping over the events, so we don't know f*N_i.
      // However, N_i<=N, so f*N_i<f*N, so we can use f*N in place of f*N_i to limit the number of floats we hold.
      unsigned int const maxPrunedSize = std::max((unsigned int) std::ceil(double(nEntries)*frac_rem), nevts_bin_skipThr+1);
      if (verbosity>=TVar::ERROR) MELAout << "\t- Maximum size of weight collection per bin for hypothesis " << ihypo << " is determined to be " << maxPrunedSize << "." << endl;
      auto& weights_hypo = weights_all.at(ihypo);
      for (auto& weights_hypo_bin:weights_hypo) weights_hypo_bin.reserve(maxPrunedSize);

      frac_list.push_back(frac);
      frac_rem_list.push_back(frac_rem);
      tolerance_list.push_back(tolerance);
      maxPrunedSize_list.push_back(maxPrunedSize);
      skip_list.at(ihypo) = (frac_rem==0.);

      ihypo++;
    }
  }

  for (int ev=0; ev<nEntries; ev++){
    tree->getEvent(ev);
    HelperFunctions::progressbar(ev, nEntries);

    int ibin = varbin_rule(tree, binning, var_vals);
    if (ibin<0) ibin=0;
    else if (ibin>=(int) nbins) ibin = nbins-1;

    for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
      if (skip_list.at(ihypo)) continue;

      float wgt_combined = std::abs(wgt_rule_list.at(ihypo)(tree, wgt_vals_list.at(ihypo)));
      if (wgt_combined==0.f) continue;

      counts_all.at(ihypo).at(ibin)++;
      sumWgts_all.at(ihypo).at(ibin) += wgt_combined;

      auto& weights = weights_all.at(ihypo).at(ibin);
      if (weights.size()<maxPrunedSize_list.at(ihypo)) HelperFunctions::addByHighest(weights, wgt_combined, false);
      else if (wgt_combined>weights.back()){
        weights.pop_back();
        HelperFunctions::addByHighest(weights, wgt_combined, false);
      }
    }
  }
  for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
    if (verbosity>=TVar::ERROR) MELAout << "Hypothesis " << ihypo << ":" << endl;
    double const& frac_rem = frac_rem_list.at(ihypo);
    double const& tolerance = tolerance_list.at(ihypo);
    if (skip_list.at(ihypo)) continue;

    bool isValidHypothesis = true;
    {
      double sumWgts_hypo = 0;
      for (unsigned int ibin=0; ibin<nbins; ibin++) sumWgts_hypo += sumWgts_all.at(ihypo).at(ibin);
      if (sumWgts_hypo==0.){
        if (verbosity>=TVar::ERROR) MELAout << "\t- Hypothesis " << ihypo << " probably does not apply to this sample because sums of weights after reweighting are all 0." << endl;
        isValidHypothesis = false;
      }
    }

    for (unsigned int ibin=0; ibin<nbins; ibin++){
      if (verbosity>=TVar::ERROR) MELAout << "\t- Bin " << ibin << ":" << endl;
      unsigned int const& nevts_bin = counts_all.at(ihypo).at(ibin);
      auto const& sumWgts_hypo_bin = sumWgts_all.at(ihypo).at(ibin);
      auto const& weights_hypo_bin = weights_all.at(ihypo).at(ibin);

      if (nevts_bin<nevts_bin_skipThr){
        if (verbosity>=TVar::ERROR) MELAout << "\t\t- Nevts = " << nevts_bin << "<3, skipping..." << endl;
        if (!isValidHypothesis) res.at(ihypo).at(ibin) = -99;
        continue;
      }

      unsigned int const nWeights_stored = weights_hypo_bin.size();
      unsigned int const nevts_bin_fracThr = std::min(nWeights_stored, std::max((unsigned int) std::ceil(double(nevts_bin)*frac_rem), nevts_bin_skipThr));

      unsigned int const index_entry = nevts_bin_fracThr-2;
      unsigned int const index_entry_prev = index_entry+1;
      float threshold = (weights_hypo_bin.at(index_entry_prev) + weights_hypo_bin.at(index_entry))*0.5;

      if (verbosity>=TVar::ERROR) MELAout << "\t\t- Raw threshold before checking tolerance: " << threshold << " / largest weight: " << weights_hypo_bin.front() << endl;
      if (weights_hypo_bin.front()<threshold*tolerance) threshold = -1; // Prevent false-positives
      res.at(ihypo).at(ibin) = threshold;
      MELAout << "\t\t- Final threshold: " << threshold << endl;
      MELAout << "\t\t- Number of events rejected: " << (threshold>0.f ? index_entry : (unsigned int) 0) << " / " << nevts_bin << endl;
      double sumWgts_hypo_bin_afterThr = sumWgts_hypo_bin;
      if (threshold>0.f){
        for (auto const& wgt:weights_hypo_bin){
          if (wgt<threshold) break;
          sumWgts_hypo_bin_afterThr -= double(wgt);
        }
      }
      MELAout
        << "\t\t- Sum of weights after / before thresholds: " << sumWgts_hypo_bin_afterThr << " / " << sumWgts_hypo_bin
        << " (fraction = " << sumWgts_hypo_bin_afterThr / sumWgts_hypo_bin << ")"
        << endl;
    }
  }

  return res;
}
