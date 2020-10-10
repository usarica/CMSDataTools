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
  thr_Neff = std::min(thr_Neff, double(nEntries)/3.*2.);
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
