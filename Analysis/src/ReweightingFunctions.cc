#include "ReweightingFunctions.h"
#include "SampleHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


std::vector<float*> ReweightingFunctions::getWeightRefs(CJLSTTree* tree, const std::vector<TString>& strWeights){
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


float ReweightingFunctions::getSimpleWeight(CJLSTTree* tree, const std::vector<float*>& vals){
  float res=0;
  if (!tree) return res;
  res=1;
  for (auto const& v:vals){ if (v) res *= *v; }
  return res;
}

