#include "SystematicsHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


SystematicsHelpers::YieldSystematic::YieldSystematic(const std::vector<ReweightingBuilder*>& inEvaluators, float(*infcn)(CJLSTTree*, const std::vector<ReweightingBuilder*>&)) : rule(infcn), evaluators(inEvaluators) {}
float SystematicsHelpers::YieldSystematic::eval(CJLSTTree* theTree) const{
  float weight=0;
  if (rule) weight=rule(theTree, evaluators);
  return weight;
}

float SystematicsHelpers::getRawSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders){
  float weight=1;
  for (auto const& builder:builders) weight *= builder->getPostThresholdWeight(theTree);
  return weight;
}
float SystematicsHelpers::getNormalizedSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders){
  float weight=0;
  assert(builders.size()==2 && builders.at(0) && builders.at(1));
  float sumWeights = builders.at(0)->getSumPostThresholdWeights(theTree);
  if (sumWeights!=0.) weight = builders.at(0)->getPostThresholdWeight(theTree)*builders.at(1)->getSumPostThresholdWeights(theTree)/sumWeights;
  return weight;
}


SystematicsHelpers::PerLeptonSystematic::PerLeptonSystematic(
  const TString inStrLepId, const std::vector<TString>& inStrVars, float(*infcn)(std::vector<short>* const&, std::vector<std::vector<float>*> const&)
) : rule(infcn), strLepId(inStrLepId), strVars(inStrVars)
{}
float SystematicsHelpers::PerLeptonSystematic::eval(CJLSTTree* theTree) const{
  auto it = componentRefs.find(theTree);
  if (it==componentRefs.cend()){
    MELAerr << "PerLeptonSystematic::eval: Could not find the variables for tree " << theTree->sampleIdentifier << ". Call PerLeptonSystematic::setup first!" << endl;
    return 0;
  }
  return rule(it->second.first, it->second.second);
}
void SystematicsHelpers::PerLeptonSystematic::setup(CJLSTTree* theTree){
  if (!theTree || strVars.empty() || strLepId=="") return;
  std::vector<short>* refId=nullptr; theTree->getValRef(strLepId, refId);
  std::vector<std::vector<float>*> refVar;
  for(TString const& s:strVars){
    vector<float>* v = nullptr;
    theTree->getValRef(s, v);
    if (!v) MELAerr << "PerLeptonSystematic::setup: Could not get the reference for " << s << endl;
    refVar.push_back(v);
  }
  if (refId && !refVar.empty()) componentRefs[theTree]=std::pair<std::vector<short>*, std::vector<std::vector<float>*>>(refId, refVar);
}

