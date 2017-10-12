#include "ReweightingFunctions.h"
#include "SampleHelpers.h"


std::vector<float> ReweightingFunctions::getWeightValues(CJLSTTree* tree, const std::vector<TString>& strWeights){
  std::vector<float> res;
  if (!tree || strWeights.empty()) return res;
  res.reserve(strWeights.size());
  for (auto const& s:strWeights){
    float w = 1;
    tree->getVal(s, w);
    res.push_back(w);
  }
  return res;
}


float ReweightingFunctions::getSimpleWeight(CJLSTTree* tree, const std::vector<TString>& strWeights){
  float res=0;
  if (!tree) return res;
  std::vector<float> weights; getWeightValues(tree, strWeights);
  res=1;
  for(auto const& w:weights) res *= w;
  return res;
}

