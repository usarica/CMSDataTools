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
  std::vector<double> const& thr_Neff_per_bin,
  TVar::VerbosityLevel verbosity
){
  unsigned int const nbins = (!binning.isValid() ? static_cast<unsigned int>(1) : binning.getNbins());

  std::vector<float> res(nbins, -1);

  int nEntries = tree->getNEvents();
  if (verbosity>=TVar::ERROR) MELAout << "ReweightingFunctions::getAbsWeightThresholdsPerBinByNeff: Determining the weight thresholds (number of events = " << nEntries << ")..." << endl;
  if (nEntries==0) return res;

  // Calculate the weight thresholds based on the Neff thresholds per bin
  std::vector<unsigned int> npos(nbins, 0);
  std::vector<double> Neff(nbins, 0);
  std::vector< std::pair<double, double> > sum_wgts(nbins, std::pair<double, double>(0, 0)); // first: w, second: w^2
  std::vector< std::vector<float> > smallest_weights(nbins, std::vector<float>());
  for (int ev=0; ev<nEntries; ev++){
    tree->getEvent(ev);
    HelperFunctions::progressbar(ev, nEntries);

    {
      float wgt_combined = std::abs(wgt_rule(tree, wgt_vals));
      int ibin = varbin_rule(tree, binning, var_vals);
      if (ibin<0) ibin=0;
      else if (ibin>=(int) nbins) ibin = nbins-1;
      HelperFunctions::addByLowest(smallest_weights.at(ibin), wgt_combined, false);
    }

    if (ev%100000 == 0 || ev == nEntries-1){
      // Reset sum_wgts and npos, and recalculate
      for (auto& sw:sum_wgts){ sw.first = sw.second = 0; }
      for (auto& nn:npos) nn = 0;
      for (unsigned int ibin=0; ibin<nbins; ibin++){
        auto const& smallest_weights_bin = smallest_weights.at(ibin);
        auto& Neff_bin = Neff.at(ibin);
        double const Neff_prev = Neff_bin;
        auto& npos_bin = npos.at(ibin);
        auto& sum_wgts_bin = sum_wgts.at(ibin);
        auto sum_wgts_bin_prev = sum_wgts_bin;

        for (auto const& wgt:smallest_weights_bin){
          sum_wgts_bin.first += wgt;
          sum_wgts_bin.second += wgt*wgt;
          if (sum_wgts_bin.second>0.) Neff_bin = std::pow(sum_wgts_bin.first, 2) / sum_wgts_bin.second;
          else Neff_bin = 0;
          npos_bin++;
          if (npos_bin>=thr_Neff_per_bin.at(ibin) && Neff_prev>Neff_bin){
            sum_wgts_bin = sum_wgts_bin_prev;
            Neff_bin = Neff_prev;
            npos_bin--;
            break;
          }
        }
      }
      if (verbosity>=TVar::ERROR) MELAout << "\t- Current Neff = " << Neff << " over " << ev+1 << " events..." << endl;
    }
    bool endEventLoop = true;
    for (unsigned int ibin=0; ibin<nbins; ibin++) endEventLoop &= (npos.at(ibin)>=thr_Neff_per_bin.at(ibin));
    if (endEventLoop) break;
  }
  for (unsigned int ibin=0; ibin<nbins; ibin++){
    auto const& smallest_weights_bin = smallest_weights.at(ibin);
    auto const& Neff_bin = Neff.at(ibin);
    auto const& npos_bin = npos.at(ibin);
    auto const& sum_wgts_bin = sum_wgts.at(ibin);

    if (verbosity>=TVar::ERROR) MELAout << "\t- Bin " << ibin << ":" << endl;
    if (!smallest_weights_bin.empty()){
      if (Neff_bin>1.) res.at(ibin) = (sum_wgts_bin.first + std::sqrt(sum_wgts_bin.second*Neff_bin)) / (Neff_bin-1.);
      if (verbosity>=TVar::ERROR){
        unsigned int nVeto = 0;
        for (auto const& wgt:smallest_weights_bin){
          if (wgt<res.at(ibin)) continue;
          nVeto++;
        }
        MELAout
          << "\t\t- " << res.at(ibin)
          << " is the default weight threshold calculated from sN=" << sum_wgts_bin.first << ", vN=" << sum_wgts_bin.second << ", nN=" << Neff_bin
          << " (N=" << npos_bin << " / " << smallest_weights_bin.size() << ", wN=" << smallest_weights_bin.at(npos_bin-1) << ", wLast=" << smallest_weights_bin.back()
          << "). Expected fraction of vetos: " << ((double) nVeto) / ((double) smallest_weights_bin.size())
          << endl;
      }
    }
    else{
      if (verbosity>=TVar::INFO) MELAout << "\t\t- No weight threshold is found." << endl;
    }
  }

  // For comparisons to brute-force procedure
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
      auto const& weights_bin=weights_all.at(ibin);
      double Neff=0;
      std::pair<double, double> sumW(0, 0);
      unsigned int ipos=0;
      for (auto const& wgt:weights_bin){
        auto sumW_new = sumW;
        sumW_new.first += wgt;
        sumW_new.second += std::pow(wgt, 2);
        double Neff_new = std::pow(sumW_new.first, 2)/sumW_new.second;
        if (ipos>=weights_bin.size()*2./3. && Neff_new<Neff) break;
        Neff = Neff_new;
        sumW = sumW_new;
        ipos++;
      }
      MELAout
        << "\t\t- " << (ipos>1 ? (weights_bin.at(ipos-1)+weights_bin.at(ipos-2))/2. : -1.)
        << " is the brute-force weight threshold calculated from sN=" << sumW.first << ", vN=" << sumW.second << ", nN=" << Neff
        << " (N=" << ipos << " / " << weights_bin.size() << ", wN=" << weights_bin.at(ipos-1) << ", wLast=" << weights_bin.back()
        << "). Expected fraction of accept: " << ((double) ipos) / ((double) weights_bin.size())
        << endl;
    }
  }

  return res;
}
std::vector<float> ReweightingFunctions::getAbsWeightThresholdsPerBinByNeff(
  BaseTree* tree,
  std::vector<float*> const& wgt_vals, ReweightingFunction_t wgt_rule,
  ExtendedBinning const& binning, std::vector<float*> const& var_vals, ReweightingVariableBinFunction_t varbin_rule,
  double thr_Neff,
  TVar::VerbosityLevel verbosity
){
  unsigned int const nbins = (!binning.isValid() ? static_cast<unsigned int>(1) : binning.getNbins());

  std::vector<float> res(nbins, -1);

  // Determine the actual Neff thresholds per bin based on the raw event distribution
  std::vector<double> thr_Neff_per_bin = getSimpleNeffThrsPerBin(tree, binning, var_vals, varbin_rule, wgt_vals, wgt_rule, thr_Neff, verbosity);

  int nEntries = tree->getNEvents();
  if (verbosity>=TVar::ERROR) MELAout << "ReweightingFunctions::getAbsWeightThresholdsPerBinByNeff: Determining the weight thresholds (number of events = " << nEntries << ")..." << endl;
  if (nEntries==0) return std::vector<float>(nbins, -1);

  // Calculate the weight thresholds based on the Neff thresholds per bin
  return getAbsWeightThresholdsPerBinByNeff(
    tree,
    wgt_vals, wgt_rule,
    binning, var_vals, varbin_rule,
    thr_Neff_per_bin, verbosity
  );
}
