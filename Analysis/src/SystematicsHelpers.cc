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
  const std::vector<TString>& inStrVars, std::pair<float, float> (*infcn)(short const& Z1Flav, short const& Z2Flav, std::vector<std::vector<float>*> const&), bool doUp_
) : rule(infcn), strVars(inStrVars), doUp(doUp_)
{}
float SystematicsHelpers::PerLeptonSystematic::eval(CJLSTTree* theTree) const{
  auto it = componentRefs.find(theTree);
  if (it==componentRefs.cend()){
    MELAerr << "PerLeptonSystematic::eval: Could not find the variables for tree " << theTree->sampleIdentifier << ". Call PerLeptonSystematic::setup first!" << endl;
    return 0;
  }
  std::pair<float, float> res = rule(*(it->second.first.at(0)), *(it->second.first.at(1)), it->second.second);
  return (doUp ? res.second : res.first);
}
void SystematicsHelpers::PerLeptonSystematic::setup(CJLSTTree* theTree){
  if (!theTree || strVars.empty()) return;
  std::vector<short*> refId; refId.assign(2, nullptr);
  theTree->getValRef("Z1Flav", refId.at(0));
  theTree->getValRef("Z2Flav", refId.at(1));
  std::vector<std::vector<float>*> refVar;
  for (TString const& s:strVars){
    vector<float>* v = nullptr;
    theTree->getValRef(s, v);
    if (!v) MELAerr << "PerLeptonSystematic::setup: Could not get the reference for " << s << endl;
    refVar.push_back(v);
  }
  if (!refVar.empty()) componentRefs[theTree]=std::pair<std::vector<short*>, std::vector<std::vector<float>*>>(refId, refVar);
}

std::pair<float, float> SystematicsHelpers::getLeptonSFSystematic(short const& Z1Flav, short const& Z2Flav, std::vector<std::vector<float>*> const& LepVars){
  std::pair<float, float> res(1, 1);
  float const& lep4pt = LepVars.at(0)->back();
  for (unsigned int ilep=0; ilep<LepVars.at(0)->size(); ilep++){
    res.first *= (1.-LepVars.at(2)->at(ilep)/LepVars.at(1)->at(ilep))*(1.-LepVars.at(4)->at(ilep)/LepVars.at(3)->at(ilep));
    res.second *= (1.+LepVars.at(2)->at(ilep)/LepVars.at(1)->at(ilep))*(1.+LepVars.at(4)->at(ilep)/LepVars.at(3)->at(ilep));
  }
  if (std::abs(Z1Flav)==169 && std::abs(Z2Flav)==169){
    if (lep4pt<7.){
      res.first *= 1.-0.032;
      res.second *= 1.+0.001;
    }
    else{
      res.first *= 1.-0.015;
      res.second *= 1.+0.001;
    }
  }
  else if (std::abs(Z1Flav)==121 && std::abs(Z2Flav)==121){
    if (lep4pt<12.){
      res.first *= 1.-0.011;
      res.second *= 1.+0.001;
    }
    else{
      res.first *= 1.-0.010;
      res.second *= 1.+0.005;
    }
  }
  else{
    if (lep4pt<7.){
      res.first *= 1.-0.080;
      res.second *= 1.+0.005;
    }
    else if (lep4pt<12.){
      res.first *= 1.-0.032;
      res.second *= 1.+0.005;
    }
    else{
      res.first *= 1.-0.010;
      res.second *= 1.+0.001;
    }
  }
  return res;
}
